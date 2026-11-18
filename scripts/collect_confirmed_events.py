#!/usr/bin/env python3
"""Collect CONFIRMED inversion / duplication / translocation events into a
publish-ready folder structure.

Writes VCF v4.3 (the de-facto standard, accepted by dbVar / DGVa) plus a flat
TSV with IS-pipeline metadata plus per-event evidence packs containing the
reference locus FASTA, target region FASTA, dnadiff .report, and progressiveMauve
XMFA — everything needed to reproduce the verdict.

CONFIRMED set (per category):
  inversion_only      : verdict == CONFIRMED_BOTH   (dnadiff + Mauve agree)
  duplication         : verdict == CONFIRMED_DNADIFF (Mauve n/a for dups)
  translocation_only  : verdict == CONFIRMED_BOTH AND deep-anchor (not artifact)

Inputs (per species):
  <results>/<sp>/rearr_validate/validation.tsv      — validator output
  <results>/<sp>/records/records.json               — reference IS loci
  <results>/<sp>/rearrangements_examples.tsv        — anchor positions / details
  <by_species>/<sp>/<assembly>.fna.gz               — assembly source

Usage:
  collect_confirmed_events.py --runs <results_dir> --src <by_species_dir> \\
      --out <output_dir>  [--species ...]
"""
import argparse, csv, gzip, json, os, re, shutil, subprocess, sys, tempfile
from collections import defaultdict

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_rearr as lr  # noqa: E402


CONFIRMED = {
    "inversion_only":     {"CONFIRMED_BOTH"},
    # CONFIRMED_DNADIFF (any dnadiff insertion in the window) is no longer a
    # confirmation: see validate_rearrangements_gold.py
    "duplication":        {"CONFIRMED_BOTH"},
    "translocation_only": {"CONFIRMED_BOTH"},   # already excludes contig-break artifacts
}
SVTYPE = {"inversion_only": "INV", "duplication": "DUP",
          "translocation_only": "BND"}
EVENT_PREFIX = {"inversion_only": "INV", "duplication": "DUP",
                "translocation_only": "TRANSLOC"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--runs", required=True, help="pipeline results dir (parent of <sp>/)")
    p.add_argument("--src", required=True, help="by_species dir with <asm>.fna.gz")
    p.add_argument("--out", required=True)
    p.add_argument("--species", nargs="*",
                   help="restrict to these species (default: all that have validation.tsv)")
    p.add_argument("--locus-window", type=int, default=30000,
                   help="bp of context around the IS locus written into evidence")
    p.add_argument("--mauve",
                   default="/global/home/users/kh36969/.conda/envs/mauve-env/bin/progressiveMauve",
                   help="progressiveMauve binary for regenerating evidence")
    p.add_argument("--regen-evidence", action="store_true",
                   help="re-run dnadiff + Mauve on each event and save outputs")
    p.add_argument("--threads", type=int, default=4)
    return p.parse_args()


def load_ref_loci(records_json):
    loci = {}
    with open(records_json) as f:
        recs = json.load(f)
    for r in recs:
        rid = r.get("ref_id") or r.get("is110_id")
        src = r.get("source", {})
        ie = r.get("is_element", {})
        contig = src.get("contig")
        s = ie.get("source_start") or src.get("transposase_start")
        e = ie.get("source_end") or src.get("transposase_end")
        strand = src.get("transposase_strand", "+")
        if rid and contig and s and e:
            local = contig.split("|", 1)[1] if "|" in contig else contig
            loci[rid] = {"contig": local, "start": int(s), "end": int(e), "strand": strand}
    return loci


def load_examples_details(tsv):
    """ref_id+tgt_asm+category -> details_string."""
    out = {}
    with open(tsv) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            key = (row["category"], row["ref_id"], row["assembly"])
            out[key] = row["details"]
    return out


def anchor_targets(details):
    d = lr.parse_details(details)
    return [(c, s0, e0) for c, s0, e0, st in d["up"] + d["down"]]


def decompress_assembly(src_dir, asm, work):
    for ext in (".fna.gz", ".fa.gz", ".fna", ".fa"):
        cand = os.path.join(src_dir, asm + ext)
        if os.path.exists(cand):
            out = os.path.join(work, asm + ".fa")
            if cand.endswith(".gz"):
                with gzip.open(cand, "rt") as fi, open(out, "w") as fo:
                    shutil.copyfileobj(fi, fo)
            else:
                shutil.copy(cand, out)
            return out
    return None


def extract_region(fa, contig, start, end, window, label, out_fh):
    try:
        clen = fa.get_reference_length(contig)
    except (KeyError, ValueError):
        return None
    s = max(0, start - window)
    e = min(clen, end + window)
    seq = fa.fetch(contig, s, e)
    if not seq:
        return None
    out_fh.write(f">{label}\n")
    for i in range(0, len(seq), 80):
        out_fh.write(seq[i:i+80] + "\n")
    return (s, e, clen)


def write_vcf_header(fh, contigs, species):
    fh.write("##fileformat=VCFv4.3\n")
    fh.write("##source=Cross_reference_IS pipeline + dnadiff(MUMmer4) + progressiveMauve\n")
    fh.write(f"##species={species}\n")
    fh.write("##ALT=<ID=INV,Description=\"Inversion\">\n")
    fh.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
    fh.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    fh.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"POS is the first base of the IS element at which the target genome is rearranged; the breakpoints of the rearrangement are not resolved\">\n")
    fh.write("##INFO=<ID=IS_END,Number=1,Type=Integer,Description=\"Last base of the IS element (1-based)\">\n")
    fh.write("##INFO=<ID=IS_STRAND,Number=1,Type=String,Description=\"Transposase strand\">\n")
    fh.write("##INFO=<ID=TGT_ANCHORS,Number=.,Type=String,Description=\"Anchor hit positions in the target assembly (contig:start-end, 1-based)\">\n")
    fh.write("##INFO=<ID=VAL_METRIC,Number=1,Type=String,Description=\"Second-evidence metric from validate_rearrangements_gold.py\">\n")
    fh.write("##INFO=<ID=IS_MEDIATED,Number=0,Type=Flag,Description=\"Event is IS-element-mediated\">\n")
    fh.write("##INFO=<ID=IS_FAMILY,Number=1,Type=String,Description=\"IS element family\">\n")
    fh.write("##INFO=<ID=REF_ID,Number=1,Type=String,Description=\"Cross_reference_IS pipeline ref_id\">\n")
    fh.write("##INFO=<ID=TGT_ASM,Number=1,Type=String,Description=\"Target assembly accession\">\n")
    fh.write("##INFO=<ID=VAL_DNADIFF,Number=1,Type=String,Description=\"MUMmer dnadiff verdict\">\n")
    fh.write("##INFO=<ID=VAL_MAUVE,Number=1,Type=String,Description=\"progressiveMauve verdict\">\n")
    fh.write("##INFO=<ID=MIN_ANCHOR_END_DIST,Number=1,Type=Integer,Description=\"Min bp from anchor to target contig end\">\n")
    fh.write("##INFO=<ID=CONFIRMED_BY,Number=.,Type=String,Description=\"Gold-standard tool(s) that independently confirmed this event\">\n")
    fh.write("##FILTER=<ID=PASS,Description=\"Confirmed by both gold-standard tools (or dnadiff-only for DUP)\">\n")
    for name, length in contigs:
        fh.write(f"##contig=<ID={name},length={length}>\n")
    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def group_target_events(flat, sp, slack=1000):
    """Mark records that describe the SAME event in the target genome.

    One row is one (category, reference element, target assembly). Copies of
    the same IS family at different reference loci can all point at one and
    the same rearranged place in the target, which then appears as several
    confirmed events. Rows of one category and one target assembly whose
    target anchors lie on the same contigs at overlapping positions (within
    `slack` bp) get one `target_event` id; count events with that column.
    Nothing is removed."""
    parent = list(range(len(flat)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    by = defaultdict(list)
    for i, r in enumerate(flat):
        by[(r["category"], r["tgt_asm"], tuple(sorted(r["_anchors"])))].append(i)
    for idx in by.values():
        for x in range(len(idx)):
            for y in range(x + 1, len(idx)):
                a, b = flat[idx[x]]["_anchors"], flat[idx[y]]["_anchors"]
                if all(a[c][0] - slack <= b[c][1] and b[c][0] - slack <= a[c][1] for c in a):
                    parent[find(idx[x])] = find(idx[y])
    members = defaultdict(list)
    for i in range(len(flat)):
        members[find(i)].append(i)
    for n, root in enumerate(sorted(members, key=lambda r: min(members[r])), 1):
        for i in members[root]:
            flat[i]["target_event"] = f"{sp}_T{n:04d}"
            flat[i]["n_records_same_target_event"] = len(members[root])


def collect_for_species(sp, args, out_root):
    runs = args.runs
    sp_dir = os.path.join(runs, sp)
    val_tsv = os.path.join(sp_dir, "rearr_validate", "validation.tsv")
    rec_json = os.path.join(sp_dir, "records", "records.json")
    ex_tsv = os.path.join(sp_dir, "rearrangements_examples.tsv")
    if not (os.path.exists(val_tsv) and os.path.exists(rec_json) and os.path.exists(ex_tsv)):
        print(f"[skip] {sp}: missing inputs", file=sys.stderr)
        return None

    ref_loci = load_ref_loci(rec_json)
    details_by = load_examples_details(ex_tsv)
    # the examples file holds 50 rows per category; validation.tsv may have
    # been made from the uncapped *_all.tsv files, whose events then had no
    # details -> empty evidence and no VCF record
    for extra in ("rearrangements_inversion_all.tsv", "rearrangements_duplication_all.tsv",
                  "rearrangements_translocation_all.tsv"):
        pth = os.path.join(sp_dir, extra)
        if os.path.exists(pth):
            for k, v in load_examples_details(pth).items():
                details_by.setdefault(k, v)

    sp_out = os.path.join(out_root, sp)
    os.makedirs(sp_out, exist_ok=True)
    ev_dir = os.path.join(sp_out, "evidence")
    os.makedirs(ev_dir, exist_ok=True)

    rows = list(csv.DictReader(open(val_tsv), delimiter="\t"))
    confirmed, seen_ev = [], set()
    for r in rows:
        if r["verdict"] not in CONFIRMED.get(r["category"], set()):
            continue
        k = (r["category"], r["ref_id"], r["tgt_asm"])
        if k in seen_ev:          # the same validation row twice is one event
            continue
        seen_ev.add(k)
        confirmed.append(r)
    print(f"  {sp}: {len(confirmed)} confirmed of {len(rows)} validated", file=sys.stderr)
    # no early return: the per-species files are always rewritten, so the
    # combined tables can never pick up a stale file from an earlier run
    skipped = defaultdict(int)

    src_dir = os.path.join(args.src, sp)
    work = tempfile.mkdtemp(prefix=f"collect_{sp}_", dir=out_root)

    flat = []
    contig_len = {}
    vcf_records = []
    cache_fa = {}

    def open_asm(asm):
        if asm in cache_fa:
            return cache_fa[asm]
        fa_path = decompress_assembly(src_dir, asm, work)
        if not fa_path:
            cache_fa[asm] = None; return None
        cache_fa[asm] = pysam.FastaFile(fa_path)
        cache_fa[f"{asm}__path"] = fa_path
        return cache_fa[asm]

    for idx, r in enumerate(confirmed, 1):
        cat = r["category"]
        ref_id = r["ref_id"]; ref_asm = r["ref_asm"]; tgt_asm = r["tgt_asm"]
        locus = ref_loci.get(ref_id)
        if not locus:
            skipped["ref_id not in records"] += 1
            continue
        details = details_by.get((cat, ref_id, tgt_asm), "")
        anchors = anchor_targets(details)
        if not anchors:
            skipped["no anchor details for this event"] += 1
            continue
        # species in the ID: INV_0001 used to exist once per species in the
        # combined table
        event_id = f"{sp}_{EVENT_PREFIX[cat]}_{idx:04d}"

        # --- per-event evidence pack ---
        ev = os.path.join(ev_dir, event_id)
        # start from an empty directory: IDs are positional, so a rerun with
        # other inputs left another event's dnadiff report under this name
        shutil.rmtree(ev, ignore_errors=True)
        os.makedirs(ev, exist_ok=True)
        ref_fa_obj = open_asm(ref_asm)
        qry_fa_obj = open_asm(tgt_asm)
        if not ref_fa_obj or not qry_fa_obj:
            skipped["assembly FASTA not found"] += 1
            shutil.rmtree(ev, ignore_errors=True)
            continue
        ref_region_fa = os.path.join(ev, "ref_locus.fa")
        with open(ref_region_fa, "w") as fh:
            ref_span = extract_region(ref_fa_obj, locus["contig"], locus["start"],
                                      locus["end"], args.locus_window, "REF_LOCUS", fh)
        if not ref_span:
            skipped["reference region not extractable"] += 1
            shutil.rmtree(ev, ignore_errors=True)
            continue
        qry_region_fa = os.path.join(ev, "target_region.fa")
        anchor_records = []
        with open(qry_region_fa, "w") as fh:
            seen_c = {}
            for (c, s, e) in anchors:
                seen_c.setdefault(c, [s, e])
                seen_c[c][0] = min(seen_c[c][0], s); seen_c[c][1] = max(seen_c[c][1], e)
            for c, (s, e) in seen_c.items():
                sp_span = extract_region(qry_fa_obj, c, s, e,
                                         args.locus_window, f"TGT_{c}", fh)
                if sp_span:
                    anchor_records.append({"contig": c, "anchor_start": s,
                                           "anchor_end": e, "extracted_span": sp_span})

        if not anchor_records:
            # none of the anchor contigs is in the target FASTA: nothing to
            # show for this event -- it was still written as confirmed
            skipped["target region not extractable"] += 1
            shutil.rmtree(ev, ignore_errors=True)
            continue

        # --- regenerate dnadiff + Mauve evidence (optional) ---
        dnadiff_report = ""
        mauve_xmfa = ""
        if args.regen_evidence:
            dd_pre = os.path.join(ev, "dnadiff")
            try:
                subprocess.run(["dnadiff", "-p", dd_pre, ref_region_fa, qry_region_fa],
                               check=True, stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL, timeout=300)
                if os.path.exists(dd_pre + ".report"):
                    dnadiff_report = dd_pre + ".report"
            except Exception:
                pass
            try:
                xmfa = os.path.join(ev, "mauve.xmfa")
                subprocess.run([args.mauve, f"--output={xmfa}",
                                ref_region_fa, qry_region_fa],
                               check=True, stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL, cwd=ev, timeout=900)
                if os.path.exists(xmfa):
                    mauve_xmfa = xmfa
            except Exception:
                pass

        # --- event.json (structured metadata) ---
        ev_meta = {
            "event_id": event_id,
            "species": sp,
            "category": cat,
            "svtype": SVTYPE[cat],
            "ref_assembly": ref_asm,
            "target_assembly": tgt_asm,
            "ref_id": ref_id,
            "is_family": "IS110",
            "is_locus": {"contig": locus["contig"],
                         "start": locus["start"], "end": locus["end"],
                         "strand": locus["strand"]},
            "target_anchors": anchor_records,
            "validation": {
                "dnadiff": r.get("dnadiff"),
                "mauve": r.get("mauve"),
                "mauve_metric": r.get("mauve_metric"),
                "min_anchor_end_dist_bp": int(r.get("min_anchor_end_dist", -1) or -1),
                "verdict": r["verdict"],
                "confirmed_by_software": [
                    n for n, ok in [
                        ("MUMmer4_dnadiff_v4.0.1",  r.get("dnadiff") == "YES"),
                        ("progressiveMauve_v2.4.0", r.get("mauve")   == "YES"),
                    ] if ok
                ],
                "tools": [
                    {"name": "MUMmer4 dnadiff", "version": "4.0.1",
                     "ref": "Marçais et al. 2018, PLOS Comput Biol"},
                    {"name": "progressiveMauve", "version": "2.4.0",
                     "ref": "Darling et al. 2010, PLOS ONE"},
                ],
                "params": {"locus_window_bp": args.locus_window,
                           "contig_break_window_bp": 5000,
                           "min_identity_pct": 95},
            },
        }
        with open(os.path.join(ev, "event.json"), "w") as f:
            json.dump(ev_meta, f, indent=2)

        # --- VCF record ---
        # What is known is WHERE THE ELEMENT IS in the reference assembly and
        # that the target genome is rearranged there. The breakpoints of the
        # rearrangement itself are not resolved here, so the record is
        # IMPRECISE and carries no END/SVLEN. It used to be written as a
        # precise <INV>/<DUP> spanning exactly the IS element (a ~1.3 kb
        # "inversion" with SVLEN one short), and translocations as two BND
        # records that were not valid mates.
        sv_chrom = f"{ref_asm}|{locus['contig']}"
        sv_pos = locus["start"]                       # 1-based first base of the element
        try:
            ref_base = ref_fa_obj.fetch(locus["contig"], sv_pos - 1, sv_pos).upper() or "N"
        except (KeyError, ValueError):
            ref_base = "N"
        _cby = []
        if r.get("dnadiff") == "YES": _cby.append("MUMmer4_dnadiff")
        if r.get("mauve")   == "YES": _cby.append("progressiveMauve")
        cby_vcf = ",".join(_cby) if _cby else "."
        # INFO values: no '=', ';', ',', or spaces (',' separates list items)
        clean = lambda v: (str(v) if v not in (None, "") else ".").replace(  # noqa: E731
            "=", ":").replace(";", "|").replace(",", "|").replace(" ", "_")
        if ref_base not in "ACGTN":
            ref_base = "N"                      # IUPAC codes are not allowed in REF
        tgt_pos = ",".join(f"{a['contig']}:{a['anchor_start'] + 1}-{a['anchor_end']}".replace(
                               ";", "_").replace("=", "_").replace(" ", "_")
                           for a in anchor_records) or "."
        info = (f"IMPRECISE;SVTYPE={SVTYPE[cat]};IS_END={locus['end']};IS_STRAND={locus['strand']};"
                f"IS_MEDIATED;IS_FAMILY=IS110;REF_ID={ref_id};TGT_ASM={tgt_asm};"
                f"TGT_ANCHORS={tgt_pos};VAL_DNADIFF={clean(r.get('dnadiff'))};"
                f"VAL_MAUVE={clean(r.get('mauve'))};"
                f"VAL_METRIC={clean(r.get('mauve_metric'))};"
                f"MIN_ANCHOR_END_DIST={clean(r.get('min_anchor_end_dist'))};"
                f"CONFIRMED_BY={cby_vcf}")
        # translocation: a single breakend with unknown mate ("N."), the VCF
        # form for "joined to something else here"; <TRA> is not a permitted
        # symbolic allele
        alt = f"{ref_base}." if cat == "translocation_only" else f"<{SVTYPE[cat]}>"
        vcf_records.append((sv_chrom, sv_pos, event_id, ref_base, alt, ".", "PASS", info))
        contig_len[sv_chrom] = ref_fa_obj.get_reference_length(locus["contig"])

        # which gold-standard tools backed this event
        confirmed_by = []
        if r.get("dnadiff") == "YES":
            confirmed_by.append("MUMmer4_dnadiff_v4.0.1")
        if r.get("mauve") == "YES":
            confirmed_by.append("progressiveMauve_v2.4.0")
        confirmed_by_str = ";".join(confirmed_by) if confirmed_by else ""

        # --- flat TSV row ---
        flat.append({
            "event_id": event_id, "species": sp, "category": cat,
            "svtype": SVTYPE[cat], "ref_id": ref_id, "ref_asm": ref_asm,
            "tgt_asm": tgt_asm, "is_contig": locus["contig"],
            "is_start": locus["start"], "is_end": locus["end"],
            "is_strand": locus["strand"],
            "target_anchor_contigs": ";".join(a["contig"] for a in anchor_records),
            "n_target_contigs": len(anchor_records),
            "val_dnadiff": r.get("dnadiff", ""),
            "val_mauve": r.get("mauve", ""),
            "mauve_metric": r.get("mauve_metric", ""),
            "min_anchor_end_dist_bp": r.get("min_anchor_end_dist", ""),
            "verdict": r["verdict"],
            "confirmed_by_software": confirmed_by_str,
            "evidence_dir": os.path.relpath(ev, sp_out),
            "_anchors": {a["contig"]: (a["anchor_start"], a["anchor_end"])
                         for a in anchor_records},
        })

    group_target_events(flat, sp)

    # --- write outputs for this species ---
    # contigs touched (for VCF header)
    contigs = sorted(contig_len.items())
    vcf_path = os.path.join(sp_out, "confirmed.vcf")
    with open(vcf_path, "w") as fh:
        write_vcf_header(fh, contigs, sp)
        for rec in sorted(vcf_records, key=lambda x: (x[0], x[1])):
            fh.write("\t".join(str(x) for x in rec) + "\n")

    tsv_path = os.path.join(sp_out, "confirmed.tsv")
    cols = ["event_id", "species", "category", "svtype", "ref_id", "ref_asm", "tgt_asm",
            "is_contig", "is_start", "is_end", "is_strand", "target_anchor_contigs",
            "n_target_contigs", "val_dnadiff", "val_mauve", "mauve_metric",
            "min_anchor_end_dist_bp", "verdict", "confirmed_by_software", "evidence_dir",
            "target_event", "n_records_same_target_event"]
    with open(tsv_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in flat:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")
    if skipped:
        print(f"  {sp}: {sum(skipped.values())} confirmed events NOT written: "
              f"{dict(skipped)}", file=sys.stderr)

    # close pysam handles + cleanup tmp
    for k, v in cache_fa.items():
        if hasattr(v, "close"):
            v.close()
    shutil.rmtree(work, ignore_errors=True)

    return {"species": sp, "n_confirmed": len(flat), "n_validated": len(rows),
            "n_target_events": len({r["target_event"] for r in flat}),
            "n_confirmed_in_validation": len(confirmed),
            "n_skipped": sum(skipped.values()),
            "by_category": {cat: sum(1 for r in flat if r["category"] == cat)
                            for cat in CONFIRMED}}


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    all_species = sorted(d for d in os.listdir(args.runs)
                         if os.path.isdir(os.path.join(args.runs, d, "rearr_validate")))
    target = args.species or all_species
    print(f"Collecting from {len(target)} species ...", file=sys.stderr)

    summary = []
    for sp in target:
        s = collect_for_species(sp, args, args.out)
        if s:
            summary.append(s)

    # all-species summary TSV
    summ_path = os.path.join(args.out, "all_species_summary.tsv")
    with open(summ_path, "w") as f:
        # n_target_events: records that describe the same place in the same
        # target genome counted once (see group_target_events)
        f.write("species\tn_validated\tn_confirmed\tn_inversion\tn_duplication\tn_translocation"
                "\tn_target_events\n")
        for s in summary:
            bc = s.get("by_category", {})
            f.write(f"{s['species']}\t{s['n_validated']}\t{s['n_confirmed']}\t"
                    f"{bc.get('inversion_only',0)}\t{bc.get('duplication',0)}\t"
                    f"{bc.get('translocation_only',0)}\t{s.get('n_target_events', 0)}\n")

    # combined VCF + TSV
    all_vcf = os.path.join(args.out, "all_confirmed.vcf")
    all_tsv = os.path.join(args.out, "all_confirmed.tsv")
    with open(all_vcf, "w") as fout:
        # full header (INFO/ALT/contig lines): with the bare three-line
        # header parsers read every INFO field as an untyped string
        all_contigs = {}
        for s_ in summary:
            spv = os.path.join(args.out, s_["species"], "confirmed.vcf")
            if os.path.exists(spv):
                for line in open(spv):
                    if line.startswith("##contig=<ID="):
                        nm, ln = line[len("##contig=<ID="):].rstrip(">\n").rsplit(",length=", 1)
                        all_contigs[nm] = int(ln)
        write_vcf_header(fout, sorted(all_contigs.items()), "all")
        for s in summary:
            sp = s["species"]
            spv = os.path.join(args.out, sp, "confirmed.vcf")
            if not os.path.exists(spv): continue
            for line in open(spv):
                if line.startswith("#"): continue
                fout.write(line)
    with open(all_tsv, "w") as fout:
        header_written = False
        for s in summary:
            sp = s["species"]
            spt = os.path.join(args.out, sp, "confirmed.tsv")
            if not os.path.exists(spt): continue
            with open(spt) as fin:
                hdr = fin.readline()
                if not header_written:
                    fout.write(hdr); header_written = True
                for line in fin:
                    fout.write(line)

    print(f"\nWrote: {summ_path}\n       {all_vcf}\n       {all_tsv}", file=sys.stderr)
    for s in summary:
        bc = s.get("by_category", {})
        print(f"  {s['species']:<28} confirmed={s['n_confirmed']:<4} "
              f"INV={bc.get('inversion_only',0)} DUP={bc.get('duplication',0)} "
              f"TRANSLOC={bc.get('translocation_only',0)}", file=sys.stderr)


if __name__ == "__main__":
    main()

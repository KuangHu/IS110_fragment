"""Shared anchor-placement, allele-decomposition and dedup helpers.

Ported from fna_based_mgefinder_project (70_allele_reconstructor.py,
lib_insert.py, lib/tsd.py), where each of these fixed a measured bug. The
reasons are kept next to the code so they are not "simplified" away again.

  place_pair()            one scored up/down placement per assembly, with a
                          runner-up margin, instead of every up x down pair
  decompose_alleles()     empty (short) vs filled (long) interval -> exact
                          insert coordinates, junction microhomology, lost bp
  canonical_insert_key()  orientation- AND junction-invariant insert key
  assembly_core() /       GCA_x.v and GCF_x.v are the same assembly; count it
  representative_asm()    once
"""
import hashlib
import re
from collections import namedtuple

# ---------------------------------------------------------------- sequences

# Full IUPAC table: with only ACGTN, ambiguity codes (R/Y/K/M/...) were left
# uncomplemented, so a minus-strand sequence carrying one was not the true
# reverse complement.
_RC = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn",
                    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def revcomp(seq):
    return seq.translate(_RC)[::-1]


def md5(s):
    return hashlib.md5(s.encode()).hexdigest()[:12] if s else "EMPTY"


def canon_key(s):
    """Orientation-invariant key: the same element on opposite strands must
    hash identically."""
    if not s:
        return ""
    r = revcomp(s)
    return hashlib.md5((s if s <= r else r).encode()).hexdigest()[:16]


def canonical_insert_key(long_seq, start, ilen):
    """Orientation- AND junction-invariant key for an insert of `ilen` bp that
    starts at `start` in `long_seq`.

    When the junction carries a direct repeat, the same insertion can be
    written several equally valid ways (the boundary may sit anywhere inside
    the repeat) and the aligner's choice is not symmetric under reverse
    complementation. Measured in the fna project: at junction overlap >= 1,
    0.0% of loci kept their key when re-decided in the revcomp frame, and
    1.67% of E. coli events were one event counted twice. So slide the insert
    across the whole repeat and take the smallest canonical key.
    """
    if ilen <= 0 or not long_seq:
        return ""
    lo = hi = start
    while lo > 0 and long_seq[lo - 1] == long_seq[lo - 1 + ilen]:
        lo -= 1
    while hi + ilen < len(long_seq) and long_seq[hi] == long_seq[hi + ilen]:
        hi += 1
    return min(canon_key(long_seq[s:s + ilen]) for s in range(lo, hi + 1))


# --------------------------------------------------------------- assemblies

_ASM_RE = re.compile(r"^GC[AF]_(\d+)\.(\d+)")


def split_target(tname):
    """'GCA_000005845.2|U00096.3' -> ('GCA_000005845.2', 'U00096.3')."""
    if "|" in tname:
        return tuple(tname.split("|", 1))
    return tname, tname


def assembly_core(asm):
    """GenBank and RefSeq carry every assembly twice (GCA_x / GCF_x, often at
    several versions). On the 13,027-assembly E. coli DB this collapses to
    7,723 distinct assemblies -- counting both inflates every empty/filled
    tally by ~1.7x. Non-NCBI names are returned unchanged."""
    m = _ASM_RE.match(asm)
    return m.group(1) if m else asm


def _asm_rank(asm):
    """Higher is preferred: RefSeq over GenBank, then the newest version."""
    m = _ASM_RE.match(asm)
    if not m:
        return (0, 0, asm)
    return (1 if asm.startswith("GCF_") else 0, int(m.group(2)), asm)


def representative_asm(assemblies):
    """{core: chosen assembly} -- deterministic, independent of input order."""
    best = {}
    for a in assemblies:
        c = assembly_core(a)
        if c not in best or _asm_rank(a) > _asm_rank(best[c]):
            best[c] = a
    return best


# ---------------------------------------------------------------- PAF / hits

Hit = namedtuple("Hit", "qname qlen qs qe strand tname tlen ts te nmatch block ident cov")


def parse_paf_line(line):
    c = line.rstrip("\n").split("\t")
    if len(c) < 12:
        return None
    qlen, qs, qe = int(c[1]), int(c[2]), int(c[3])
    nmatch, block = int(c[9]), int(c[10])
    return Hit(c[0], qlen, qs, qe, c[4], c[5], int(c[6]), int(c[7]), int(c[8]),
               nmatch, block,
               nmatch / block * 100 if block else 0.0,
               (qe - qs) / qlen * 100 if qlen else 0.0)


# <ref_id>__up<D> / <ref_id>__down<D>; legacy single-D names carry no D
ANCHOR_RE = re.compile(r"^(.+)__(up|down)(\d*)$")


def load_anchor_hits(paf_path, min_identity, min_coverage, row_counts=None):
    """PAF -> {(ref_id, side, D): {tname: [Hit]}}, D None for legacy names.

    D is part of the key: pooling up1000 with down80000 hits (what the old
    stage-4 parser did) pairs anchors from different distances.

    If `row_counts` (a dict) is given it is filled with raw PAF rows per
    anchor, so callers can tell which anchors hit minimap2's -N cap."""
    from collections import defaultdict
    hits = defaultdict(lambda: defaultdict(list))
    n_total = n_kept = 0
    with open(paf_path) as f:
        for line in f:
            h = parse_paf_line(line)
            if h is None:
                continue
            n_total += 1
            if row_counts is not None:
                row_counts[h.qname] = row_counts.get(h.qname, 0) + 1
            if h.ident < min_identity or h.cov < min_coverage:
                continue
            m = ANCHOR_RE.match(h.qname)
            if not m:
                continue
            n_kept += 1
            D = int(m.group(3)) if m.group(3) else None
            hits[(m.group(1), m.group(2), D)][h.tname].append(h)
    return hits, n_total, n_kept


def gap_tolerance(expected_v0, tol_bp=100, tol_frac=0.002):
    """Absolute tolerance plus a small divergence allowance.

    The old tolerance was max(200, 10% of the expected distance). At D=40 kb
    that is +/-8 kb and at 80 kb +/-16 kb -- wider than the element, so a
    filled 1.3 kb IS110 site classified as 'empty'. The tolerance must not
    scale with D faster than real indel divergence does."""
    return tol_bp + tol_frac * max(0, expected_v0)


def classify_gap(gap, expected_v0, tnp_len, tol, max_removed=None,
                 min_removed=None):
    """Classify an anchor gap against the source (filled) state.

    removed = expected_v0 - gap is how much DNA the target lacks relative to
    the source. An empty site lacks at least the transposase; the removed
    length at an empty site IS the element length."""
    delta = gap - expected_v0
    if abs(delta) <= tol:
        return "v0_filled"
    if delta > tol:
        return "v1plus_filled"
    removed = -delta
    if max_removed is not None and removed > max_removed + tol:
        return "very_short"
    if removed >= tnp_len - tol and (min_removed is None or removed >= min_removed):
        return "empty"
    return "intermediate"


def place_pair(ups, downs, min_gap=-50, max_gap=200000, margin_ratio=0.05):
    """Choose ONE up/down placement among all candidate pairs in one assembly.

    Every up x down combination on a contig used to be emitted and counted. On
    repeat-rich loci the two anchors then land on DIFFERENT copies of a repeat
    and fabricate an allele (fna project: a 42,367 bp 'allele'), and a genome
    with two anchor copies is counted several times. Instead:

      * each anchor's score is normalised by its own length, so a long weak
        hit cannot outvote two clean ones;
      * the gap must lie in [min_gap, max_gap]. min_gap is NEGATIVE on
        purpose: anchors that abut the element overlap by the target-site
        duplication at an empty site, and a clean empty site has gap 0 --
        both used to be discarded;
      * the runner-up pair is scored too; if the best is not clearly better
        (relative margin < margin_ratio) the placement is 'ambiguous' and is
        reported, not forced.

    `ups`/`downs` are Hit lists. Returns a dict or None.
    """
    def inner_edges(u, d):
        # Inward-facing anchor edges on the target, PROJECTED through any
        # unaligned query bases. Hits pass at >= min_coverage (80%), so up to
        # 100 bp of a 500-bp anchor may be clipped; using the raw alignment
        # ends added the clipped inner bases to the gap (a true 3,200-bp
        # filled gap measured 3,400 -> 'v1plus_filled') and to the fetched
        # interval.
        if u.strand == "+":
            return u.te + (u.qlen - u.qe), d.ts - d.qs
        # anchors swap sides on the minus strand
        return d.te + d.qs, u.ts - (u.qlen - u.qe)

    pairs = []
    for ui, u in enumerate(ups):
        for di, d in enumerate(downs):
            if u.tname != d.tname or u.strand != d.strand:
                continue
            s, e = inner_edges(u, d)
            gap = e - s
            if gap < min_gap or gap > max_gap:
                continue
            score = u.nmatch / max(1, u.qlen) + d.nmatch / max(1, d.qlen)
            pairs.append((score, gap, ui, di, u, d, s, e))
    if not pairs:
        return None
    # One candidate per LOCUS: keep a pair only if each anchor is the other's
    # nearest valid partner. Ranking every up x down combination by score
    # alone let a full-length up hit at copy 1 pair with a full-length down
    # hit at copy 2 and outscore both true (partly clipped) same-locus pairs:
    # a 'resolved' 12 kb allele that exists in no genome.
    near_u, near_d = {}, {}
    for c in pairs:
        _, gap, ui, di = c[:4]
        if ui not in near_u or gap < near_u[ui][0]:
            near_u[ui] = (gap, di)
        if di not in near_d or gap < near_d[di][0]:
            near_d[di] = (gap, ui)
    loci = [c for c in pairs
            if near_u[c[2]][1] == c[3] and near_d[c[3]][1] == c[2]]
    if not loci:
        return None
    in_locus_u = {c[2] for c in loci}
    in_locus_d = {c[3] for c in loci}
    # A stray extra copy of ONE anchor (not part of any locus of its own) is
    # still a rival for the anchor it could pair with: one up hit with two
    # equally good down hits must stay ambiguous, not be settled by distance.
    strays = [c for c in pairs
              if (c[2] in in_locus_u) != (c[3] in in_locus_d)]
    cands = [(c[0], c[4], c[5], c[6], c[7]) for c in loci + strays]
    # deterministic order: score, then target position
    cands.sort(key=lambda x: (-x[0], x[3], x[4]))
    # the placement itself is always a locus; strays only compete as runner-up
    best_locus = min(loci, key=lambda c: (-c[0], c[6], c[7]))
    score, u, d, s, e = best_locus[0], best_locus[4], best_locus[5], best_locus[6], best_locus[7]
    cands = [(score, u, d, s, e)] + [c for c in cands
                                      if not (c[1] is u and c[2] is d)]
    if len(cands) == 1:
        margin, status, second = 1.0, "unique", None
    else:
        second = cands[1][0]
        margin = (score - second) / abs(score) if score else 0.0
        status = "ambiguous" if margin < margin_ratio else "resolved"
    return {"tname": u.tname, "strand": u.strand, "up": u, "down": d,
            "start": s, "end": e, "gap": e - s,
            "n_pairs": len(cands), "score": score, "second_score": second,
            "margin": margin, "status": status}


def place_per_assembly(up_by_tname, down_by_tname, **kw):
    """{tname: [Hit]} x2 -> {core: placement}, one per distinct assembly.

    Groups contigs by assembly, places once per assembly, then keeps a single
    representative per GCA/GCF core. The placement dict gains 'assembly' and
    'core'."""
    by_asm_up, by_asm_down = {}, {}
    for t, hs in up_by_tname.items():
        by_asm_up.setdefault(split_target(t)[0], []).extend(hs)
    for t, hs in down_by_tname.items():
        by_asm_down.setdefault(split_target(t)[0], []).extend(hs)
    by_core = {}
    for asm in set(by_asm_up) & set(by_asm_down):
        by_core.setdefault(assembly_core(asm), []).append(asm)
    out = {}
    # sorted(): iteration order must not depend on the hash seed -- it decides
    # which observations survive a top-N cut downstream.
    for core in sorted(by_core):
        # Preferred twin first, but fall through to the other one: choosing
        # the representative BEFORE placing lost the site whenever the
        # preferred copy had its two anchors on different contigs.
        for asm in sorted(by_core[core], key=_asm_rank, reverse=True):
            pl = place_pair(by_asm_up[asm], by_asm_down[asm], **kw)
            if pl is not None:
                pl["assembly"], pl["core"] = asm, core
                out[core] = pl
                break
    return out


# ------------------------------------------------------------ decomposition

def decompose(short, long_):
    """Exact common prefix / suffix. Returns (lcp, lcs, inserted, amb) with
    amb = lcp + lcs - len(short); negative means target bases were lost."""
    n = min(len(short), len(long_))
    p = 0
    while p < n and short[p] == long_[p]:
        p += 1
    s = 0
    while s < n - p and short[len(short) - 1 - s] == long_[len(long_) - 1 - s]:
        s += 1
    return p, s, long_[p:len(long_) - s], p + s - len(short)


def _slide(long_, l0, l1):
    """How far the insert long_[l0:l1] can slide and still give the same two
    alleles = the junction microhomology (a TSD when there is one). Reported,
    never resolved: picking one placement would fake base-level precision."""
    left = 0
    while l0 - 1 - left >= 0 and long_[l0 - 1 - left] == long_[l1 - 1 - left]:
        left += 1
    right = 0
    while l1 + right < len(long_) and long_[l0 + right] == long_[l1 + right]:
        right += 1
    return left + right


def tolerant_decompose(short, long_, min_cov, min_ident, min_insert, dominance):
    """Alignment-based decomposition that tolerates strain-level SNPs/indels.

    Exact prefix matching stops at the first SNP (~1/100 bp between strains),
    which left most real insertions undecomposed (fna: median explained
    fraction 0.26 for failures). A strict gate keeps this from manufacturing
    insertions: coverage, identity, and ONE dominant extra block.
    """
    from Bio.Align import PairwiseAligner
    if not short or not long_ or len(long_) <= len(short):
        return None
    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                              open_gap_score=-12, extend_gap_score=-0.4)
    try:
        aln = aligner.align(short, long_)[0]
    except (ValueError, OverflowError, MemoryError):
        return None
    sb, lb = aln.aligned
    if len(sb) == 0:
        return None
    # native ints: Bio's coordinates are numpy.int64, which json.dump rejects
    # (Stage 6 crashed whenever the consensus call was alignment-based)
    blocks = [((int(s0), int(s1)), (int(l0), int(l1)))
              for (s0, s1), (l0, l1) in zip(sb, lb)]
    inserts, dels, cols, matches = [], [], 0, 0
    # Terminal gaps count too. Only gaps BETWEEN aligned blocks were tallied,
    # so a second large extra block at the start or end of the long allele
    # escaped the dominance test.
    prev_s, prev_l = 0, 0
    blocks.append(((len(short), len(short)), (len(long_), len(long_))))
    for (s0, s1), (l0, l1) in blocks:
        sg, lg = s0 - prev_s, l0 - prev_l
        if sg == 0 and lg > 0:
            inserts.append((lg, prev_s, prev_l, l0))
        elif lg == 0 and sg > 0:
            dels.append(sg)
        elif sg > 0 and lg > 0:
            inserts.append((lg, prev_s, prev_l, l0))
            dels.append(sg)
        matches += sum(1 for a, b in zip(short[s0:s1], long_[l0:l1]) if a == b)
        cols += s1 - s0
        prev_s, prev_l = s1, l1
    if cols == 0 or not inserts:
        return None
    cov, ident = cols / len(short), matches / cols
    inserts.sort(key=lambda x: -x[0])
    big, sp, bl0, bl1 = inserts[0]
    others = [x[0] for x in inserts[1:]] + dels
    second = max(others) if others else 0
    dom = big / max(1, second)
    ok = (cov >= min_cov and ident >= min_ident and big >= min_insert
          and (second == 0 or dom >= dominance))
    return {"ok": ok, "coverage": cov, "identity": ident,
            "largest_insert": big, "second_indel": second,
            "lcp": sp, "lcs": len(short) - sp,
            # sp is a SHORT-allele offset; bl0 is where the insert starts in
            # the LONG allele. They differ by upstream indels -- re-slicing the
            # long allele with sp was a 31% error rate in the fna project.
            "insert_start_long": bl0, "insert_len": bl1 - bl0,
            "target_lost": sum(dels),
            "ambiguity": _slide(long_, bl0, bl1)}


def hamming_split(short, long_, min_ident=0.95, min_insert=50,
                  window=20, max_window_mm=8, clean_junction=30):
    """SNP-tolerant, indel-free decomposition in O(len): find the offset where
    long = short[:off] + insert + short[off:] with the fewest mismatches.

    The alignment in tolerant_decompose() is quadratic, so it is skipped for
    long intervals (D >= 5 kb: ~10 kb x 11 kb cells). Exact matching then
    needs a 10 kb empty interval with zero SNPs against the source, which
    strain-level divergence (~1/100-1/1000 bp) almost never gives -- so
    elements whose ends lie > 1 kb from the transposase fell back to the
    symmetric guess. A single insertion needs no gaps elsewhere, only a split
    point, so mismatches can simply be counted on both sides of every
    candidate split.

    Refuses (returns None) when:
      * the best split is below `min_ident`;
      * any mismatch lies within `clean_junction` bp of the junction on
        either side. Both flanks are compared anchored at their OUTER ends,
        so an indel anywhere in a flank shifts every base between it and the
        junction: the shifted run always abuts the junction. A clean junction
        therefore rules out a net indel in either flank -- which would make
        the insert length wrong by the indel size (a 3-bp deletion 7 bp from
        the junction once gave a confident 1,503-bp call for a 1,500-bp
        element);
      * any `window`-bp stretch holds more than `max_window_mm` mismatches
        (compensating indels leave a shifted run that does not reach the
        junction).
    A real SNP next to the junction also refuses; another empty will vote.
    Ties take the rightmost split, the same convention as decompose().
    """
    n, ilen = len(short), len(long_) - len(short)
    if n == 0 or ilen < min_insert:
        return None
    # suf[i] = mismatches of short[i:] against the end of long_
    suf = [0] * (n + 1)
    for i in range(n - 1, -1, -1):
        suf[i] = suf[i + 1] + (short[i] != long_[i + ilen])
    best, off, pre = None, 0, 0
    for i in range(n + 1):
        tot = pre + suf[i]
        if best is None or tot <= best:
            best, off = tot, i
        if i < n:
            pre += short[i] != long_[i]
    ident = 1 - best / n
    if ident < min_ident:
        return None
    mm = [i for i in range(off) if short[i] != long_[i]]
    mm += [i for i in range(off, n) if short[i] != long_[i + ilen]]
    if any(off - clean_junction <= i < off + clean_junction for i in mm):
        return None
    for k in range(len(mm) - max_window_mm):
        if mm[k + max_window_mm] - mm[k] < window:
            return None
    return {"offset": off, "insert_start_long": off, "insert_len": ilen,
            "identity": ident, "n_mismatch": best,
            "ambiguity": _slide(long_, off, off + ilen)}


def decompose_alleles(short, long_, min_insert=50, min_cov=0.85,
                      min_ident=0.95, dominance=3.0, tolerant=True,
                      hamming_fallback=True):
    """Explain the long (filled) interval as the short (empty) one plus an
    insert. Exact first, alignment-tolerant fallback.

    When the alignment is skipped (`tolerant=False`, interval too long) and
    `hamming_fallback` is set, hamming_split() is tried instead.

    Returns a dict:
      method            exact | tolerant_alignment | hamming_split | none
      event_class       pure_insertion | insertion_target_retained |
                        replacement | below_min_insert | substitution_only |
                        undecomposed
      insert_start      insert start in the LONG allele (use this to slice)
      insert_len, insert_seq
      offset            insertion point in the SHORT allele, so that
                        short[:offset] + insert + short[offset:] ~= long
      junction_microhomology_bp   bp the junction can slide (TSD candidate)
      target_lost_bp    short-allele bases present on neither side
    """
    res = {"method": "none", "event_class": "undecomposed", "insert_start": -1,
           "insert_len": 0, "insert_seq": "", "offset": -1,
           "junction_microhomology_bp": 0, "target_lost_bp": 0,
           "identity": None, "coverage": None}
    if len(long_) == len(short):
        res["event_class"] = "substitution_only"
        return res
    if len(long_) < len(short):
        return res
    lcp, lcs, ins, amb = decompose(short, long_)
    if lcp + lcs >= len(short):
        # every base of the short allele is on one side of the junction
        start = lcp
        res.update(method="exact", insert_start=start, insert_len=len(ins),
                   insert_seq=ins, offset=lcp,
                   junction_microhomology_bp=_slide(long_, start, start + len(ins)),
                   target_lost_bp=0, identity=1.0, coverage=1.0,
                   event_class="pure_insertion" if not short
                   else "insertion_target_retained")
    elif tolerant:
        td = tolerant_decompose(short, long_, min_cov, min_ident, min_insert,
                                dominance)
        if td and td["ok"]:
            st, ln = td["insert_start_long"], td["insert_len"]
            res.update(method="tolerant_alignment", insert_start=st,
                       insert_len=ln, insert_seq=long_[st:st + ln],
                       offset=td["lcp"],
                       junction_microhomology_bp=td["ambiguity"],
                       target_lost_bp=td["target_lost"],
                       identity=td["identity"], coverage=td["coverage"],
                       event_class="replacement" if td["target_lost"]
                       else "insertion_target_retained")
    elif hamming_fallback:
        hs = hamming_split(short, long_, min_ident, min_insert)
        if hs:
            st, ln = hs["insert_start_long"], hs["insert_len"]
            res.update(method="hamming_split", insert_start=st, insert_len=ln,
                       insert_seq=long_[st:st + ln], offset=hs["offset"],
                       junction_microhomology_bp=hs["ambiguity"],
                       target_lost_bp=0, identity=hs["identity"], coverage=1.0,
                       event_class="insertion_target_retained")
    if res["method"] != "none" and res["insert_len"] < min_insert:
        res["event_class"] = "below_min_insert"
    return res


# ---------------------------------------------------------------------- TSD

def find_tsd(ins_seq, left_flank, right_flank, min_tsd=3, max_tsd=30):
    """Longest complexity-filtered direct repeat. Returns (len, seq, side, conf).

    RECORDED, NEVER A FILTER: IS110/IS1111 transpose without a TSD (fna
    project: 0/13 IS110 events carried one vs 43/170 for other families), so
    gating on TSD would discard exactly the family this pipeline was built for.
    """
    if not ins_seq:
        return 0, "", ".", "none"

    def complex_enough(sub):
        if len(set(sub)) < 3:
            return False
        return max(sub.count(b) for b in set(sub)) / len(sub) <= 0.75

    for n in range(min(max_tsd, len(ins_seq)), min_tsd - 1, -1):
        cands = []
        if len(right_flank) >= n and ins_seq[:n].upper() == right_flank[:n].upper():
            cands.append((ins_seq[:n].upper(), "right"))
        if len(left_flank) >= n and ins_seq[-n:].upper() == left_flank[-n:].upper():
            cands.append((ins_seq[-n:].upper(), "left"))
        for sub, side in cands:
            if complex_enough(sub):
                return n, sub, side, ("high" if n >= 6 else "medium" if n >= 4 else "low")
    return 0, "", ".", "none"

#!/bin/bash
# Completion markers for the resumable steps of the pipeline scripts.
#
# "Output file exists and is not empty" is not "this step finished with these
# inputs and options": a header-only .domtbl of a killed hmmsearch, a partial
# PAF, or the result of an earlier run with other inputs were all skipped as
# done. A step is done only if its marker file holds the signature of the
# current inputs and options; outputs are written to <name>.tmp and renamed.
#
#   SIG=$(stage_sig "step name" "$OPTIONS" "$(file_fp "$INPUT")" "$PREVIOUS_SIG")
#   if stage_done "$OUT.done" "$SIG"; then echo skip
#   else  cmd > "$OUT.tmp" && mv "$OUT.tmp" "$OUT" && stage_mark "$OUT.done" "$SIG"; fi
#
# Include the previous step's signature so that a change upstream reruns
# everything downstream.

stage_sig()  { printf '%s\n' "$@" | md5sum | cut -d' ' -f1; }
file_fp()    { stat -c '%n:%s:%Y' "$@" 2>/dev/null || echo "missing:$*"; }
# fingerprint (name, size, mtime) of the files listed on stdin, one per line:
# a list of PATHS is not a fingerprint of what the files contain
files_fp()   { xargs -d '\n' -r stat -c '%n:%s:%Y' 2>&1 | md5sum | cut -d' ' -f1; }
# stage_done MARKER SIG [PRODUCT ...]: the marker holds SIG and every product
# exists (an empty file is a product). Without PRODUCT arguments the product
# of X.done is X: a marker whose output had been removed was "done".
stage_done() {
    local marker="$1" sig="$2" p
    shift 2
    if [ ! -f "$marker" ] || [ "$(cat "$marker")" != "$sig" ]; then return 1; fi
    if [ $# -eq 0 ]; then
        case "$marker" in *.done) set -- "${marker%.done}" ;; esac
    fi
    for p in "$@"; do
        if [ ! -e "$p" ]; then return 1; fi
    done
    return 0
}
stage_mark() { printf '%s\n' "$2" > "$1"; }
# count FASTA records / data lines; 0 is a result, not an error (grep -c exits 1)
count_fasta() { grep -c '^>' "$1" || true; }
count_rows()  { grep -cv '^#' "$1" || true; }
# require_unique_seq_names DIR: sequence names must be unique across the
# *.fna files that are concatenated into one minimap2 target. Two genomes
# that both have a contig "ctg" become ONE target, and halves of an IS in two
# different genomes were joined into one "insertion".
require_unique_seq_names() {
    local dup
    dup=$(find "$1" -name "*.fna" -print0 | xargs -0 -r grep -h '^>' | awk '{print $1}' \
          | sort | uniq -d | head -3 || true)
    if [ -n "$dup" ]; then
        echo "ERROR: sequence names occur in more than one genome file under $1, e.g.:" >&2
        echo "$dup" >&2
        echo "Give the records unique names (prefix the assembly) before searching." >&2
        return 1
    fi
}

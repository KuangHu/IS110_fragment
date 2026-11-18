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
stage_done() { [ -f "$1" ] && [ "$(cat "$1")" = "$2" ]; }
stage_mark() { printf '%s\n' "$2" > "$1"; }
# count FASTA records / data lines; 0 is a result, not an error (grep -c exits 1)
count_fasta() { grep -c '^>' "$1" || true; }
count_rows()  { grep -cv '^#' "$1" || true; }

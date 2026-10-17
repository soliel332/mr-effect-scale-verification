#!/usr/bin/env bash
set -euo pipefail

root="${SCHEME_C_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
python_bin="${PYTHON_BIN:-python3}"
ldsc_software="${LDSC_SOFTWARE_DIR:-$root/ldsc_software}"
ldsc_pydeps="${LDSC_PYDEPS_DIR:-$root/ldsc_pydeps}"
export PYTHONPATH="$ldsc_pydeps:$ldsc_software"
ld_scores="${LD_SCORE_DIR:-$root/ldsc_reference_repo/GRCh37/results/ldscores/EUR/}"
sumstats="$root/ldsc_sumstats"
results="$root/ldsc_results"

run_anchor() {
  local anchor="$1"
  shift
  local rg_files="$sumstats/$anchor.sumstats.gz"
  local label
  for label in "$@"; do
    rg_files="$rg_files,$sumstats/$label.sumstats.gz"
  done
  "$python_bin" "$ldsc_software/ldsc.py" \
    --rg "$rg_files" \
    --ref-ld-chr "$ld_scores" \
    --w-ld-chr "$ld_scores" \
    --out "$results/rg_${anchor}_full"
}

run_anchor T1 T2 T3 ADHD ASD NEURO ASTHMA
run_anchor T2 T3 ADHD ASD NEURO ASTHMA
run_anchor T3 ADHD ASD NEURO ASTHMA
run_anchor ADHD ASD NEURO ASTHMA
run_anchor ASD NEURO ASTHMA
run_anchor NEURO ASTHMA

"$python_bin" "$root/parse_scheme_c_ldsc.py"

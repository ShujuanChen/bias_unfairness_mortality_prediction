#!/usr/bin/env bash
#
# Write the synthetic stand-in for every restricted input, overwriting whatever
# is there. No number produced from these files means anything.
#
#   bash data/regenerate_synthetic_inputs.sh
#
set -euo pipefail

D="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for s in UKB census PMR imd lookup HSE; do
  Rscript "$D/$s/make_synthetic_$s.R" >/dev/null || { echo "$s generator failed" >&2; exit 1; }
  echo "  $s"
done

echo "synthetic inputs written under $(basename "$D")/"

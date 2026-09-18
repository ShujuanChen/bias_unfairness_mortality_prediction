#!/usr/bin/env bash
#
# Multiplier bootstrap replicates of the whole pipeline, weights and models
# included.

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
CODE="$ROOT/code"

PYTHON="${REWEIGHTING_PYTHON:-python3}"

RUN="${REWEIGHTING_RUN_ID:-}"
if [ -z "$RUN" ]; then
  echo "No REWEIGHTING_RUN_ID. run_pipeline.R sets it, or set it to the run to write into." >&2
  exit 1
fi
TEMP="$ROOT/temp/$RUN"

cfg() { "$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1]))[sys.argv[2]][sys.argv[3]])' "$ROOT/framework_config.json" "$@"; }

FIRST="${1:-1}"
LAST="${2:-$(cfg phase5 replicates)}"
JOBS="${3:-2}"
CAUSES=$("$PYTHON" -c 'import json,sys; print(" ".join(json.load(open(sys.argv[1]))["phase2"]["outcomes"]))' "$ROOT/framework_config.json")

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

KEY=$(cfg phase3 weight_method)
case "$KEY" in
  census_*) REFERENCE="census" ;;
  *)        REFERENCE="hse" ;;
esac
case "$KEY" in
  hse_superlearner)       WEIGHT_CMD=(Rscript "$CODE/phase3_weighting/s1_estimate_hse.R" superlearner) ;;
  hse_lassologit)         WEIGHT_CMD=(Rscript "$CODE/phase3_weighting/s1_estimate_hse.R" lassologit) ;;
  hse_raking)             WEIGHT_CMD=(Rscript "$CODE/phase3_weighting/s3_moment_weights.R" raking) ;;
  hse_entropy_balancing)  WEIGHT_CMD=(Rscript "$CODE/phase3_weighting/s3_moment_weights.R" entropy_balancing) ;;
  census_superlearner)    WEIGHT_CMD=(Rscript "$CODE/phase3_weighting/s2_estimate_census.R" superlearner) ;;
  *) echo "phase3.weight_method is $KEY, which this runner does not know"; exit 1 ;;
esac

SOURCES="ukb,pmr,ukbw_$KEY"
BOOT="$TEMP/bootstrap_multiplier"
LOGS="$TEMP/bootstrap_logs/$KEY"
mkdir -p "$LOGS"

log() { echo "[$(date -u +'%F %H:%M:%SZ')] $*" | tee -a "$LOGS/status"; }

has_all_sources() {
  local f="$1"
  shift
  [ -s "$f" ] || return 1
  "$PYTHON" -c '
import csv, json, os, sys
try:
    with open(sys.argv[1], newline="") as fh:
        have = set(next(csv.reader(fh)))
except Exception:
    raise SystemExit(1)
with open(os.path.join(sys.argv[2], "framework_config.json")) as fh:
    years = json.load(fh)["phase1"]["horizons_years"]
want = {f"{s}_pred_{y}y" for s in sys.argv[3:] for y in years}
raise SystemExit(0 if want <= have else 1)
' "$f" "$ROOT" "$@" 2>/dev/null
}

cd "$ROOT"
[ -d "$TEMP/population" ] || { log "no population under temp/. Run phase 1 first."; exit 1; }

NCAUSE=$(set -- $CAUSES; echo $#)
NFOLD=$("$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1]))["phase2"]["splits"]["folds"])' "$ROOT/framework_config.json")
NFIT=$((NFOLD + 2))
log "=== multiplier bootstrap, replicates $FIRST to $LAST, $KEY"
log "causes: $CAUSES"
log "one replicate is 1 participation model and $((NCAUSE * NFIT)) outcome fits, $JOBS at a time"

FAILED=0
LINKED=""
for d in harmonised population folds; do
  [ -d "$TEMP/$d" ] && LINKED="$LINKED $TEMP/$d"
done
if [ -z "$LINKED" ]; then
  echo "temp/ has no harmonisation, population or folds to link." >&2; exit 1
fi
MARKER="$LOGS/.point_estimate_untouched"
: > "$MARKER"

for b in $(seq "$FIRST" "$LAST"); do
  D=$(printf "%s/replicate_%04d" "$BOOT" "$b")
  if [ -f "$D/rates.csv" ]; then
    log "replicate $b already done"
    continue
  fi

  BR=$(printf "%s/replicates/%s_%04d" "$TEMP" "$KEY" "$b")
  log "--- replicate $b, tree ${BR#$ROOT/}"

  if ! "$PYTHON" "$HERE/multiplier_bootstrap.py" \
      --replicate "$b" > "$LOGS/draw_$b.log" 2>&1; then
    log "  FAILED draw, see $LOGS/draw_$b.log"; FAILED=1; continue
  fi

  mkdir -p "$BR"
  for d in harmonised population folds; do
    mkdir -p "$BR/$d"
    for f in "$TEMP/$d"/*; do
      ln -sfn "$f" "$BR/$d/$(basename "$f")"
    done
  done
  rm -f "$BR/harmonised/$REFERENCE.csv"
  cp "$D/$REFERENCE.csv" "$BR/harmonised/$REFERENCE.csv"

  export REWEIGHTING_TEMP_ROOT="$BR"
  export REWEIGHTING_BOOTSTRAP_DIR="$D"

  if [ -f "$BR/weights/$KEY/ukb_weights.csv" ]; then
    log "  participation model already there, kept"
  elif ! "${WEIGHT_CMD[@]}" > "$LOGS/p3_$b.log" 2>&1; then
    log "  FAILED participation model, see $LOGS/p3_$b.log"; FAILED=1
    unset REWEIGHTING_TEMP_ROOT REWEIGHTING_BOOTSTRAP_DIR
    continue
  fi

  if [ "${REWEIGHTING_BOOTSTRAP_STAGE:-all}" = "participation" ]; then
    log "  participation stage only, weights in $BR/weights/$KEY"
    unset REWEIGHTING_TEMP_ROOT REWEIGHTING_BOOTSTRAP_DIR
    continue
  fi

  queue=()
  for c in $CAUSES; do
    for s in ukb "ukbw_$KEY"; do
      [ -f "$BR/models/$c/$s/single/model_state.pt" ] || queue+=("$c|$s|")
    done
    for f in $(seq 0 $((NFOLD - 1))); do
      [ -f "$BR/models/$c/pmr/fold_$f/model_state.pt" ] || queue+=("$c|pmr|$f")
    done
  done
  log "  ${#queue[@]} of $((NCAUSE * NFIT)) outcome fits to do"

  ok=1
  for (( k = 0; k < JOBS; k++ )); do
    spec=""
    ci=-1
    for c_deal in $CAUSES; do
      ci=$((ci + 1))
      [ "$((ci % JOBS))" -eq "$k" ] || continue
      for i in "${!queue[@]}"; do
        IFS='|' read -r c s f <<< "${queue[$i]}"
        [ "$c" = "$c_deal" ] || continue
        spec="${spec:+$spec,}${c}:${s}${f:+:$f}"
      done
    done
    [ -n "$spec" ] || continue
    "$PYTHON" "$CODE/phase2_prediction_bias/s2_train.py" --tasks "$spec" \
      > "$LOGS/train_${b}_slot${k}.log" 2>&1 &
    pids[$k]=$!
    labels[$k]="slot $k"
  done
  for j in "${!pids[@]}"; do
    wait "${pids[$j]}" || { log "  FAILED ${labels[$j]}, see $LOGS/train_${b}_slot${j}.log"; ok=0; }
  done
  unset pids labels
  if [ "$ok" -eq 0 ]; then
    FAILED=1
    unset REWEIGHTING_TEMP_ROOT REWEIGHTING_BOOTSTRAP_DIR
    continue
  fi

  for c in $CAUSES; do
    if has_all_sources "$BR/predictions/${c}_individual_predictions.csv" ${SOURCES//,/ }; then
      continue
    fi
    if ! "$PYTHON" "$CODE/phase2_prediction_bias/s3_assemble.py" \
        --cause "$c" --sources "$SOURCES" > "$LOGS/assemble_${b}_${c}.log" 2>&1; then
      log "  FAILED assemble $c, see $LOGS/assemble_${b}_${c}.log"; ok=0
    fi
  done
  if [ "$ok" -eq 0 ]; then
    FAILED=1
    unset REWEIGHTING_TEMP_ROOT REWEIGHTING_BOOTSTRAP_DIR
    continue
  fi

  if ! "$PYTHON" "$CODE/phase2_prediction_bias/s4_cumulative_incidence.py" \
      --sources "$SOURCES" > "$LOGS/cif_$b.log" 2>&1; then
    log "  FAILED cumulative incidence, see $LOGS/cif_$b.log"; FAILED=1
    unset REWEIGHTING_TEMP_ROOT REWEIGHTING_BOOTSTRAP_DIR
    continue
  fi

  if ! "$PYTHON" "$HERE/replicate_rates.py" \
      --sources "$SOURCES" --out "$D/rates.csv" \
      > "$LOGS/rates_$b.log" 2>&1; then
    log "  FAILED rates, see $LOGS/rates_$b.log"; FAILED=1
    unset REWEIGHTING_TEMP_ROOT REWEIGHTING_BOOTSTRAP_DIR
    continue
  fi

  if [ "${REWEIGHTING_BOOTSTRAP_KEEP:-1}" = "0" ]; then
    rm -rf "$BR"
  else
    log "  kept $BR"
  fi
  unset REWEIGHTING_TEMP_ROOT REWEIGHTING_BOOTSTRAP_DIR
  log "  replicate $b done"
done

TOUCHED=$(find $LINKED -type f -newer "$MARKER" 2>/dev/null || true)
if [ -n "$TOUCHED" ]; then
  log "STOPPING: a replicate wrote into the point estimate run, which it must only read:"
  echo "$TOUCHED" | while read -r f; do log "  $f"; done
  exit 1
fi
rm -f "$MARKER"

log "=== finished, failures=$FAILED"
log "replicates in $BOOT, one rates table each. Phase 5 stage 2 forms the intervals."
exit "$FAILED"

#!/bin/bash
# -----------------------------------------------------------------------------
# run_tests.sh
#
# Run selected cases of the ibScourFoam test suite and collect their
# metrics. Run directories go to $IBS_RUN_ROOT (default <repo>/run/tests,
# ignored by git).
#
# Usage:
#     tests/run_tests.sh [options] <case>[:variant,variant...] ... | all
#
# Options:
#     --submit        submit the run as a Slurm batch job ($IBS_SUBMIT,
#                     default sbatch; account and partition from
#                     IBS_ACCOUNT and IBS_PARTITION if set) instead of
#                     running here
#     --ntasks N      Slurm tasks / concurrent serial runs (default 8)
#     --time HH:MM:SS Slurm time limit (default 01:00:00)
#     --collect       only collect and print existing metrics
#
# Cases: pipelineScour ibChannelFlat sandSlideRelax humpMigration
#
# Examples:
#     tests/run_tests.sh --submit ibChannelFlat sandSlideRelax
#     tests/run_tests.sh ibChannelFlat:lam_f025,ras_tilt2
#     tests/run_tests.sh --collect all
# -----------------------------------------------------------------------------

TESTS_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=common/runFunctions.sh
. "$TESTS_DIR/common/runFunctions.sh"

ALL_CASES="pipelineScour ibChannelFlat sandSlideRelax humpMigration"
: "${IBS_SUBMIT:=sbatch}"
: "${IBS_ACCOUNT:=}"
: "${IBS_PARTITION:=}"

submit=false
collectOnly=false
ntasks=8
timeLimit=01:00:00
selection=()
while [ $# -gt 0 ]
do
    case "$1" in
        --submit) submit=true ;;
        --collect) collectOnly=true ;;
        --ntasks) ntasks=$2; shift ;;
        --time) timeLimit=$2; shift ;;
        -h|--help) sed -n '2,35p' "$0"; exit 0 ;;
        all) read -r -a selection <<< "$ALL_CASES" ;;
        *) selection+=("$1") ;;
    esac
    shift
done

if [ ${#selection[@]} -eq 0 ]
then
    echo "No case given (try --help)" >&2
    exit 1
fi

# Submit this script as a batch job from a clean environment
if $submit
then
    mkdir -p "$IBS_RUN_ROOT/jobs"
    cd "$IBS_RUN_ROOT/jobs" || exit 1
    site=()
    [ -n "$IBS_ACCOUNT" ] && site+=(--account="$IBS_ACCOUNT")
    [ -n "$IBS_PARTITION" ] && site+=(--partition="$IBS_PARTITION")
    jobId=$($IBS_SUBMIT \
        --job-name=ibs-tests "${site[@]}" --nodes=1 --ntasks="$ntasks" \
        --cpus-per-task=1 --mem=64G --time="$timeLimit" \
        --output="$IBS_RUN_ROOT/jobs/%x-%j.out" \
        --export=ALL,IBS_RUN_ROOT="$IBS_RUN_ROOT",IBS_MAX_JOBS="$ntasks",IBS_SOLVER="$IBS_SOLVER",IBS_ENV="$IBS_ENV",IBS_PYTHON="$IBS_PYTHON",IBS_TIMEOUT="${IBS_TIMEOUT:-}" \
        --wrap="bash $TESTS_DIR/run_tests.sh --ntasks $ntasks \
${selection[*]}" | tail -1)
    echo "Submitted job $jobId (log: $IBS_RUN_ROOT/jobs/ibs-tests-$jobId.out)"
    exit 0
fi

ibs_load_env || exit 1
set -e
# Job steps must inherit this environment
export SLURM_EXPORT_ENV=ALL
export IBS_MAX_JOBS=$ntasks

if ! $collectOnly
then
    echo "Running on $(hostname) at $(date); run root $IBS_RUN_ROOT"
    for item in "${selection[@]}"
    do
        case=${item%%:*}
        variants=()
        if [ "$item" != "$case" ]
        then
            IFS=, read -r -a variants <<< "${item#*:}"
        fi
        if [ ! -x "$TESTS_DIR/$case/Allrun" ]
        then
            echo "Unknown case $case" >&2
            continue
        fi
        echo "=== $case ${variants[*]} ($(date +%T))"
        tStart=$(date +%s)
        "$TESTS_DIR/$case/Allrun" "${variants[@]}" || \
            echo "$case: Allrun returned $?" >&2
        echo "=== $case done in $(( $(date +%s) - tStart )) s"
    done
fi

cases=()
for item in "${selection[@]}"
do
    cases+=("${item%%:*}")
done
$IBS_PYTHON "$TESTS_DIR/common/collect.py" "$IBS_RUN_ROOT" \
    --cases "${cases[@]}" --json "$IBS_RUN_ROOT/summary.json"

#------------------------------------------------------------------------------

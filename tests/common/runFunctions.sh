#!/bin/bash
# -----------------------------------------------------------------------------
# runFunctions.sh
#
# Shell helpers sourced by the Allrun/Allclean scripts of the ibScourFoam
# test suite. Source this file; do not execute it.
#
# Environment variables (all optional):
#     IBS_ENV       script that loads OpenFOAM v2606 and the ibScourFoam
#                   build, sourced if OpenFOAM is not loaded yet (default
#                   none: load OpenFOAM before running the tests)
#     IBS_RUN_ROOT  root of all run directories (default <repo>/run/tests)
#     IBS_MAX_JOBS  concurrent serial runs in a case's Allrun (default 4)
#     IBS_SOLVER    solver executable (default ibScourPimpleFoam)
#     IBS_PYTHON    Python >= 3.8 with numpy (default python3)
#     IBS_MPI_BIND  mpirun binding option (default "--bind-to core"; use
#                   "--bind-to none" if the batch system hands out hardware
#                   threads of fewer cores than tasks)
# -----------------------------------------------------------------------------

# Directory of the test suite and of the repository
IBS_TESTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IBS_REPO_DIR="$(cd "$IBS_TESTS_DIR/.." && pwd)"
export IBS_TESTS_DIR IBS_REPO_DIR

# Default locations outside the repository
: "${IBS_ENV:=}"
: "${IBS_RUN_ROOT:=$IBS_REPO_DIR/run/tests}"
: "${IBS_MAX_JOBS:=4}"
: "${IBS_SOLVER:=ibScourPimpleFoam}"
: "${IBS_PYTHON:=python3}"
: "${IBS_MPI_BIND:=--bind-to core}"
export IBS_ENV IBS_RUN_ROOT IBS_MAX_JOBS IBS_SOLVER IBS_PYTHON IBS_MPI_BIND


# Load OpenFOAM unless it is already loaded. Must be called before 'set -e'
# because OpenFOAM's etc/bashrc is not 'set -e' safe.
ibs_load_env()
{
    if [ -z "$WM_PROJECT_DIR" ] && [ -n "$IBS_ENV" ]
    then
        # shellcheck disable=SC1090
        source "$IBS_ENV" > /dev/null 2>&1 || true
    fi
    if ! command -v blockMesh > /dev/null 2>&1
    then
        echo "OpenFOAM environment not loaded (IBS_ENV=$IBS_ENV)" >&2
        return 1
    fi
}


# Run a command with its output in <case>/log.<name>.
# Usage: ibs_log <caseDir> <name> <command> [args...]
ibs_log()
{
    local caseDir=$1 name=$2
    shift 2
    ( cd "$caseDir" && "$@" > "log.$name" 2>&1 )
}


# Run a solver in <caseDir> serially or with mpirun. The solver exit
# status is recorded in <caseDir>/exitStatus.<solver> instead of aborting,
# because a crash is a test result. Parallel runs are reconstructed.
# If IBS_TIMEOUT is set (e.g. 45m) the solver is stopped after that time
# (exit status 124).
# Usage: ibs_solve <caseDir> <solver> <nProcs>
ibs_solve()
{
    local caseDir=$1 solver=$2 np=$3 status=0 famesh=no
    local limit=()
    if [ -n "$IBS_TIMEOUT" ]
    then
        limit=(timeout --signal=TERM "$IBS_TIMEOUT")
    fi
    # Height-field bed with system/finite-area/faMeshDefinition: the area
    # mesh is written by makeFaMesh, so the stock tools decompose and
    # reconstruct the bed fields
    [ -f "$caseDir/system/finite-area/faMeshDefinition" ] && famesh=yes
    if [ "$np" -gt 1 ]
    then
        foamDictionary -entry numberOfSubdomains -set "$np" \
            "$caseDir/system/decomposeParDict" > /dev/null
        if [ $famesh = yes ]
        then
            ibs_log "$caseDir" decomposePar decomposePar -force -no-finite-area
            ( cd "$caseDir" && mpirun -np "$np" $IBS_MPI_BIND makeFaMesh \
                -parallel > log.makeFaMesh 2>&1 )
        else
            ibs_log "$caseDir" decomposePar decomposePar -force
            # no area mesh on disk: decompose the bed fields by hand
            local f
            for f in "$caseDir"/0/finite-area/zb "$caseDir"/0/finite-area/zRigid
            do
                [ -f "$f" ] || continue
                ( cd "$caseDir" && $IBS_PYTHON \
                    "$IBS_TESTS_DIR/common/decompose_zb.py" . --time 0 \
                    --field "$(basename "$f")" >> log.decompose_zb 2>&1 )
            done
        fi
        # mpirun: srun --mpi=pmi2 does not work with this Open MPI build
        ( cd "$caseDir" && "${limit[@]}" mpirun -np "$np" $IBS_MPI_BIND \
            "$(command -v "$solver")" -parallel > "log.$solver" 2>&1 ) \
            || status=$?
        if [ $famesh = yes ]
        then
            # reconstructPar (v2606) corrupts the heap on finite-area
            # fields; redistributePar reconstructs volume and area fields
            ( cd "$caseDir" && mpirun -np "$np" $IBS_MPI_BIND \
                redistributePar -reconstruct -parallel \
                > log.redistributePar 2>&1 ) || true
        else
            # reconstructPar stops at the finite-area fields (no faMesh to
            # reconstruct onto): each time on its own, then the bed level
            local t
            : > "$caseDir/log.reconstructPar"
            for t in $(foamListTimes -case "$caseDir" -processor 2>/dev/null)
            do
                ( cd "$caseDir" && reconstructPar -time "$t" \
                    >> log.reconstructPar 2>&1 ) || true
            done
            ( cd "$caseDir" && $IBS_PYTHON \
                "$IBS_TESTS_DIR/common/reconstruct_zb.py" . \
                > log.reconstruct_zb 2>&1 ) || true
        fi
    else
        if [ $famesh = yes ]
        then
            ibs_log "$caseDir" makeFaMesh makeFaMesh
        fi
        ( cd "$caseDir" && "${limit[@]}" "$solver" > "log.$solver" 2>&1 ) \
            || status=$?
    fi
    echo "$status" > "$caseDir/exitStatus.$solver"
    return 0
}


# Run a ready-made case folder: copy <srcDir> to <dstDir> (below
# IBS_RUN_ROOT) and run its Allrun there; the solver exit status goes to
# <dstDir>/exitStatus.<application>.
# Usage: ibs_run_case <srcDir> <dstDir>
ibs_run_case()
{
    local src=$1 dst=$2 app status=0
    ibs_clean_dir "$dst"
    mkdir -p "$(dirname "$dst")"
    cp -r "$src" "$dst"
    ( cd "$dst" && ./Allrun > log.Allrun 2>&1 ) || status=$?
    app=$(foamDictionary -entry application -value "$dst/system/controlDict" 2>/dev/null)
    if [ -n "$app" ] && grep -q "^End" "$dst/log.$app" 2>/dev/null
    then
        echo 0 > "$dst/exitStatus.$app"
    else
        echo 1 > "$dst/exitStatus.${app:-solver}"
    fi
    return 0
}


# Prepare a run directory for ParaView: an empty <dirName>.foam file in the
# run directory and <dir>/view/ (make_view.py: links to the time
# directories that hold flow fields, plus view.foam).
# Usage: ibs_view <caseDir>
ibs_view()
{
    local caseDir=$1
    touch "$caseDir/$(basename "$caseDir").foam"
    python3 "$IBS_TESTS_DIR/common/make_view.py" "$caseDir" \
        > "$caseDir/log.make_view" 2>&1 || true
}


# Wait until fewer than IBS_MAX_JOBS background jobs are running.
ibs_throttle()
{
    while [ "$(jobs -rp | wc -l)" -ge "$IBS_MAX_JOBS" ]
    do
        wait -n || true
    done
}


# Remove a run directory below IBS_RUN_ROOT (refuses anything else).
# Usage: ibs_clean_dir <dir>
ibs_clean_dir()
{
    local dir=${1:?ibs_clean_dir: no directory}
    case "$dir" in
        "$IBS_RUN_ROOT"/?*) rm -rf -- "$dir" ;;
        *) echo "refusing to remove $dir (outside $IBS_RUN_ROOT)" >&2 ;;
    esac
}

# -----------------------------------------------------------------------------

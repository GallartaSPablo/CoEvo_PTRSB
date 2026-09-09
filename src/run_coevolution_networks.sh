#!/usr/bin/env bash
#
# Compiles src/coevolution_antagonistic_networks.c and runs it over a range
# of the 122 empirical networks published in /networks, in parallel.
#
# Usage (run from the repository root, so that ./networks/ resolves):
#   src/run_coevolution_networks.sh [-i NET_ID_INI] [-e NET_ID_END] [-t FLAG_TEMPORAL] [-j MAX_JOBS]
#
#   -i  first network id to run (default: 1)
#   -e  last network id to run, inclusive (default: 122)
#   -t  0 = write only the final state of each simulation (default)
#       1 = write the full time series of every simulation
#   -j  max number of networks simulated in parallel (default: number of CPU cores)
#
# Example: run only networks 1-10, keeping full time series, 4 at a time
#   src/run_coevolution_networks.sh -i 1 -e 10 -t 1 -j 4

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

if [ ! -d "networks" ]; then
    echo "Error: ./networks not found. Run this script from the repository root (found: $REPO_ROOT)." >&2
    exit 1
fi

NET_ID_INI=1
NET_ID_END=122
FLAG_TEMPORAL=0
MAX_JOBS="$( (command -v nproc >/dev/null 2>&1 && nproc) || sysctl -n hw.ncpu 2>/dev/null || echo 4 )"

while getopts "i:e:t:j:h" opt; do
    case "$opt" in
        i) NET_ID_INI=$OPTARG ;;
        e) NET_ID_END=$OPTARG ;;
        t) FLAG_TEMPORAL=$OPTARG ;;
        j) MAX_JOBS=$OPTARG ;;
        h)
            grep '^#' "$0" | sed 's/^#//'
            exit 0
            ;;
        *)
            echo "Usage: $0 [-i NET_ID_INI] [-e NET_ID_END] [-t FLAG_TEMPORAL] [-j MAX_JOBS]" >&2
            exit 1
            ;;
    esac
done

SOURCE="$SCRIPT_DIR/coevolution_antagonistic_networks.c"
BINARY="$SCRIPT_DIR/coevolution_antagonistic_networks"

echo "Compiling $SOURCE ..."
gcc -o "$BINARY" "$SOURCE" -O3 -lm -Wall -Wextra

LOG_DIR="files/networks/logs"
mkdir -p "files/networks/temporal" "$LOG_DIR"

echo "Running networks $NET_ID_INI..$NET_ID_END (flag_temporal=$FLAG_TEMPORAL, up to $MAX_JOBS in parallel)"

FAILED_LOG="$(mktemp)"
trap 'echo "Interrupted, stopping running jobs..."; kill 0; rm -f "$FAILED_LOG"' INT TERM

for (( NET_ID = NET_ID_INI; NET_ID <= NET_ID_END; NET_ID++ )); do
    (
        if ! "$BINARY" "$NET_ID" "$FLAG_TEMPORAL" > "$LOG_DIR/network_$(printf '%03d' "$NET_ID").log" 2>&1; then
            echo "$NET_ID" >> "$FAILED_LOG"
        fi
    ) &

    echo "Launched network $NET_ID (jobs running: $(jobs -rp | wc -l | tr -d ' '))"

    # Throttle: block until a slot frees up once MAX_JOBS are running.
    # (avoids "wait -n", not available in the bash 3.2 shipped on macOS)
    while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do
        sleep 1
    done
done

wait

if [ -s "$FAILED_LOG" ]; then
    echo "Finished with failures in networks: $(tr '\n' ' ' < "$FAILED_LOG")" >&2
    echo "See per-network logs under $LOG_DIR/" >&2
    rm -f "$FAILED_LOG"
    exit 1
fi

rm -f "$FAILED_LOG"
echo "Done: all networks $NET_ID_INI..$NET_ID_END finished successfully."

#!/bin/bash

set -e

# --- Defaults ---
INPUT=""
ENABLE_BIB=false
ENABLE_IP=false
NEVENTS=""

print_usage() {
    echo "Usage: $0 -i <input_file> [--bib] [--ip] [--nevents <number>]"
    echo "  -i <input_file>      Path to input file (must contain '_sim')"
    echo "  --bib                Enable BIB flag"
    echo "  --ip                 Enable IP flag"
    echo "  --nevents <number>   Number of events to run (positive integer)"
}

# --- Parse arguments ---
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT="$2"
            shift 2
            ;;
        --bib)
            ENABLE_BIB=true
            shift
            ;;
        --ip)
            ENABLE_IP=true
            shift
            ;;
        --nevents)
            NEVENTS="$2"
            if ! [[ "$NEVENTS" =~ ^[0-9]+$ ]]; then
                echo "Error: --nevents must be a positive integer."
                exit 1
            fi
            shift 2
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1"
            print_usage
            exit 1
            ;;
    esac
done

# --- Check input file ---
if [[ -z "$INPUT" ]]; then
    echo "Error: Input file not provided."
    print_usage
    exit 1
fi

if [[ "$INPUT" != *"_sim"* ]]; then
    echo "Error: Input file must contain '_sim' in its name: $INPUT"
    exit 1
fi

# --- Build argument list ---
ARGS=("steering/reco_steer.py" "--InFileName" "$INPUT")
$ENABLE_BIB && ARGS+=("--enableBIB")
$ENABLE_IP && ARGS+=("--enableIP")
[[ -n "$NEVENTS" ]] && ARGS+=("--NEvents" "$NEVENTS")

# --- Run ---
echo k4run "${ARGS[@]}"
k4run "${ARGS[@]}"

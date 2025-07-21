#!/bin/bash

# Usage message
usage() {
  echo "Usage: $0 --input path/to/input.slcio --output path/to/output_sim.slcio [--nevent N] [--skip N]"
  echo "Output filename must end with _sim.slcio or _sim_<number>.slcio"
  exit 1
}

# Defaults
NEVENT=1000
SKIP=0

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)
      INPUT_FILE="$2"
      shift 2
      ;;
    --output)
      OUTPUT_FILE="$2"
      shift 2
      ;;
    --nevent)
      NEVENT="$2"
      shift 2
      ;;
    --skip)
      SKIP="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
  esac
done

# Check required args
if [[ -z "$INPUT_FILE" || -z "$OUTPUT_FILE" ]]; then
  echo "Error: --input and --output are required"
  usage
fi

# Validate output filename format
if [[ ! "$OUTPUT_FILE" =~ _sim(\.slcio|_[0-9]+\.slcio)$ ]]; then
  echo "Error: Output filename must end with _sim.slcio or _sim_<number>.slcio"
  exit 1
fi

# Run ddsim
ddsim \
  --steeringFile steering/sim_steer.py \
  --inputFiles "$INPUT_FILE" \
  --numberOfEvents "$NEVENT" \
  --skipNEvents "$SKIP" \
  --outputFile "$OUTPUT_FILE" \
  --compactFile ../detector-simulation/geometries/MAIA_v0/MAIA_v0.xml

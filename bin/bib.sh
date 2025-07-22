#!/bin/bash

# Usage:
#   ./bib.sh --input path/to/summary1_DET_IP.dat --number 42.64 --output /path/to/outputdir


set -e

# Parse long options
TEMP=$(getopt -o '' \
--long input:,number:,output: \
-n 'bib.sh' -- "$@")

if [ $? != 0 ]; then
  echo "Error parsing arguments" >&2
  exit 1
fi

eval set -- "$TEMP"

# Initialize variables
input_dat=""
number=""
output_dir=""

# Extract options and values
while true; do
  case "$1" in
    --input)
      input_dat="$2"; shift 2 ;;
    --number)
      number="$2"; shift 2 ;;
    --output)
      output_dir="$2"; shift 2 ;;
    --)
      shift; break ;;
    *)
      echo "Unknown option: $1"; exit 1 ;;
  esac
done

# Validate required arguments
if [[ -z "$input_dat" || -z "$number" || -z "$output_dir" ]]; then
  echo "Usage: $0 --input INPUT_DAT --number NUMBER --output OUTPUT_DIR"
  exit 1
fi

if [[ ! -f "$input_dat" ]]; then
  echo "Error: Input file '$input_dat' not found"
  exit 1
fi

# Extract NUMBER from filename
basename=$(basename "$input_dat")
if [[ "$basename" =~ summary([0-9]+)_DET_IP\.dat ]]; then
    fileindex="${BASH_REMATCH[1]}"
else
  echo "Error: Filename must be in format summaryNUMBER_DET_IP.dat"
  exit 1
fi

mkdir -p "$output_dir/sim_mm"

output_file="$output_dir/sim_mm/BIBinput_gen_${fileindex}.slcio"

# Execute Python script
python3 $MY_MUCOLL_BASEDIR/../detector-simulation/utils/fluka_to_slcio_new.py -n "$number" "$input_dat" "$output_file"

mkdir -p "$output_dir/sim_mp"

output_file="$output_dir/sim_mp/BIBinput_gen_${fileindex}.slcio"

python3 $MY_MUCOLL_BASEDIR/../detector-simulation/utils/fluka_to_slcio_new.py -n "$number" "$input_dat" "$output_file" -i 1

sh bin/sim.sh --input ${output_dir}/sim_mm/BIBinput_gen_${fileindex}.slcio --output ${output_dir}/sim_mm/BIBinput_sim_${fileindex}.slcio --nevent 1
sh bin/sim.sh --input ${output_dir}/sim_mp/BIBinput_gen_${fileindex}.slcio --output ${output_dir}/sim_mp/BIBinput_sim_${fileindex}.slcio --nevent 1
 

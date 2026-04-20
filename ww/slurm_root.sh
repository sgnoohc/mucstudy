#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=22gb
#SBATCH --time=99:40:00

#SBATCH --job-name=ww_root
#SBATCH --output=logs/%x/%j.out
#SBATCH --error=logs/%x/%j.err

set -euo pipefail

pwd; hostname; date
echo "started running WW ROOT histogram generation"
apptainer exec \
  /cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-alma9:latest \
  bash -s << 'EOF'
    set -eo pipefail
    . /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/mucoll-stack-master-h2ssl2yh2yduqnhsv2i2zcjws74v7mcq/setup.sh # aliased by setup_mucoll
    export SAMPLEID=25620665
    export ODIR="/blue/avery/share/mucol" # NO TRAILING SLASH!!!!!!!!!!!!!!!!!!
    export MY_MUCOLL_BASEDIR="$PWD"
    export MARLIN_DLL="$MY_MUCOLL_BASEDIR/mucoll_software/MyBIBUtils/lib/libMyBIBUtils.so:${MARLIN_DLL:-}"

    mkdir -p "$ODIR/root"
    mkdir -p "$ODIR/bigroot"
    for i in $ODIR/slcio/MuMuToWW_reco_${SAMPLEID}_7*.slcio; do
      out=$ODIR/root/$(basename "${i%.slcio}")_batch.root
      python "$MY_MUCOLL_BASEDIR/ww/slcio_to_root.py" -i "$i" -o "$out"
    done
    hadd "$ODIR/bigroot/$WW_${SAMPLEID}.root" $ODIR/root/MuMuToWW_reco_${SAMPLEID}_*_batch.root
EOF

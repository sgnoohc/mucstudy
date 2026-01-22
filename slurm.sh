#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --time=02:10:00

#SBATCH --job-name=muc
#SBATCH --output=logs/%x/%A_%a.out
#SBATCH --error=logs/%x/%A_%a.err

set -euo pipefail

pwd; hostname; date
echo "started running muon collider simulation"
apptainer exec \
  /cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-alma9:latest \
  bash -s << 'EOF'
    set -eo pipefail

    N_EVENTS=3

    cd pythia
    (
      . setup.sh
      ./MuMuToZH "$N_EVENTS" "$SLURM_ARRAY_JOB_ID" "$SLURM_ARRAY_TASK_ID"
    )
    cd ..

    . /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/mucoll-stack-master-h2ssl2yh2yduqnhsv2i2zcjws74v7mcq/setup.sh # aliased by setup_mucoll
    export MY_MUCOLL_BASEDIR="$PWD"
    export MARLIN_DLL="$MY_MUCOLL_BASEDIR/mucoll_software/MyBIBUtils/lib/libMyBIBUtils.so:${MARLIN_DLL:-}"
    mkdir -p slcio
    ddsim \
      --inputFile  "pythia/hepmc/MuMuToZH_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.hepmc" \
      --steeringFile "steering/sim_steer.py" \
      --outputFile "slcio/MuMuToZH_sim_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.slcio" \
      --compactFile "mucoll_software/detector-simulation/geometries/MAIA_v0/MAIA_v0.xml" \
      --numberOfEvents "$N_EVENTS"
    k4run steering/reco_steer.py \
      --InFileName "slcio/MuMuToZH_sim_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.slcio" \
      --NEvents "$N_EVENTS" || true # continues the script after the crash

    mkdir -p root
    echo "CONVERTING SLCIO TO ROOT"
    python slcio_to_root.py \
      -i "slcio/MuMuToZH_reco_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.slcio" \
      -o "root/MuMuToZH_reco_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root"
EOF

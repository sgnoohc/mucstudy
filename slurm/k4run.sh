#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --time=04:10:00

#SBATCH --job-name=muc_k4run
#SBATCH --output=logs/%x/%j.out
#SBATCH --error=logs/%x/%j.err

set -euo pipefail
pwd; hostname; date
echo "started running muon collider simulation"
apptainer exec /cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-alma9:latest bash << 'EOF'
  set -eo pipefail
  . /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/mucoll-stack-master-h2ssl2yh2yduqnhsv2i2zcjws74v7mcq/setup.sh # aliased by setup_mucoll
  . ./setup.sh
  k4run steering/reco_steer.py --InFileName MuMuToZH_sim.slcio --NEvents 10
EOF
echo "finished running muon collider simulation"
date

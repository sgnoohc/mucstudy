#!/usr/bin/env bash
set -euo pipefail

MY_MUCOLL_BASEDIR="${MY_MUCOLL_BASEDIR:-$PWD}"

mkdir -p "${MY_MUCOLL_BASEDIR}/mucoll_software"
cd "${MY_MUCOLL_BASEDIR}/mucoll_software"
clone_if_missing() {
  local url="$1"
  local name="$2"
  if [ ! -d "$name/.git" ]; then
    git clone "$url" "$name"
  else
    echo "Repo $name already exists, skipping clone"
  fi
}
clone_if_missing https://github.com/MuonColliderSoft/ACTSTracking ACTSTracking
clone_if_missing https://github.com/madbaron/MyBIBUtils.git MyBIBUtils
clone_if_missing https://github.com/madbaron/LCIOmacros.git LCIOmacros
clone_if_missing https://github.com/madbaron/detector-simulation.git detector-simulation
clone_if_missing https://github.com/madbaron/SteeringMacros.git SteeringMacros
cd "${MY_MUCOLL_BASEDIR}"

# fixing hardcoded path in the SteeringMacros repo
SETTINGS="${MY_MUCOLL_BASEDIR}/mucoll_software/SteeringMacros/PandoraSettings/PandoraSettingsDefault.xml"
HISTO="${MY_MUCOLL_BASEDIR}/mucoll_software/SteeringMacros/PandoraSettings/PandoraLikelihoodData12EBin.xml"

if [[ ! -f "$SETTINGS" ]]; then
  echo "Error: SETTINGS file not found at $SETTINGS" >&2
  exit 1
fi

sed -i "s|<HistogramFile>.*</HistogramFile>|<HistogramFile>${HISTO}</HistogramFile>|" "${SETTINGS}"

apptainer exec /cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-alma9:latest bash << 'EOF'
  set -euo pipefail
  . /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/mucoll-stack-master-h2ssl2yh2yduqnhsv2i2zcjws74v7mcq/setup.sh # aliased by setup_mucoll

  cd pythia
  . setup.sh
  make
  ./MuMuToZH
  cd ..

  cd mucoll_software/MyBIBUtils
  mkdir -p build
  cd build/
  cmake ..
  make install
  #cd ../../..
EOF

#!/usr/bin/env bash
set -euo pipefail

MY_MUCOLL_BASEDIR="${MY_MUCOLL_BASEDIR:-$PWD}"

mkdir -p "${MY_MUCOLL_BASEDIR}/mucoll_software"
cd "${MY_MUCOLL_BASEDIR}/mucoll_software"
clone_if_missing() {
  local name="$1"
  shift
  if [ ! -d "$name/.git" ]; then
    git clone "$@" "$name"
  else
    echo "Repo $name already exists, skipping clone"
  fi
}
clone_if_missing ACTSTracking https://github.com/MuonColliderSoft/ACTSTracking
clone_if_missing MyBIBUtils https://github.com/madbaron/MyBIBUtils.git
clone_if_missing LCIOmacros https://github.com/madbaron/LCIOmacros.git
clone_if_missing detector-simulation -b KITP_10TeV https://github.com/madbaron/detector-simulation.git
clone_if_missing SteeringMacros https://github.com/madbaron/SteeringMacros.git
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
  set -eo pipefail
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

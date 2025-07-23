
## Installation

    singularity run -B /data:/data,/ceph:/ceph /cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-alma9:latest

    setup_mucoll

    mkdir mucoll_software
    cd mucoll_software

    export GIT_SSH_COMMAND="ssh -F /dev/null"

    git clone git@github.com:madbaron/MyBIBUtils.git
    git clone git@github.com:madbaron/LCIOmacros.git
    git clone -b KITP_10TeV git@github.com:madbaron/detector-simulation.git
    git clone git@github.com:madbaron/SteeringMacros.git

    cd MyBIBUtils/
    mkdir build
    cd build/
    cmake ..
    make install
    export MARLIN_DLL=/home/users/phchang/work/muc/workdir/MyBIBUtils/lib/libMyBIBUtils.so:$MARLIN_DLL
    cd ../../

    cd ../

## If git push has a problem

    GIT_SSH_COMMAND="ssh -F /dev/null"  git push

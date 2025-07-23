
## Installation

    singularity run -B /data:/data,/ceph:/ceph /cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-alma9:latest
    setup_mucoll
    source setup.sh

    mkdir mucoll_software
    cd mucoll_software
    export GIT_SSH_COMMAND="ssh -F /dev/null"
    git clone git@github.com:madbaron/MyBIBUtils.git
    git clone git@github.com:madbaron/LCIOmacros.git
    git clone -b KITP_10TeV git@github.com:madbaron/detector-simulation.git
    git clone git@github.com:madbaron/SteeringMacros.git
    git clone git@github.com:MuonColliderSoft/ACTSTracking.git
    cd MyBIBUtils/
    mkdir build
    cd build/
    cmake ..
    make install
    cd ../../../

## Next time

    singularity run -B /data:/data,/ceph:/ceph /cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-alma9:latest
    setup_mucoll
    source setup.sh


## If git push has a problem

    GIT_SSH_COMMAND="ssh -F /dev/null"  git push

## Error....

        PhotonReconstructionAlgorithm::InitialiseHistogramReading - Invalid xml file specified for pdf histograms.
    pLocalAlgorithm->ReadSettings(TiXmlHandle(pXmlElement)) throw STATUS_CODE_INVALID_PARAMETER
        in function: CreateAlgorithm
        in file:     /tmp/root/spack-stage/spack-stage-pandorasdk-3.4.2-33sx55dznpsqceud6wlroqaru7v4cpqo/spack-src/src/Managers/AlgorithmManager.cc line#: 135
    Failure in reading pandora settings, STATUS_CODE_INVALID_PARAMETER
    PandoraApi::ReadSettings(*m_pPandora, m_settings.m_pandoraSettingsXmlFile) throw STATUS_CODE_FAILURE
        in function: init
        in file:     /tmp/root/spack-stage/spack-stage-ddmarlinpandora-0.14-tr6xmfngr5hcva22kgokyrkyanrvuvah/spack-src/src/DDPandoraPFANewProcessor.cc line#: 180
    [ ERROR "DDMarlinPandora"] Failed to initialize marlin pandora: STATUS_CODE_FAILURE
    DDMarlinPandora     FATAL UNKNOWN Exception is caught
    EventLoopMgr        ERROR Unable to initialize Algorithm: DDMarlinPandora
    ServiceManager      ERROR Unable to initialize Service: EventLoopMgr
    ApplicationMgr      ERROR Application Manager Terminated with error code 1

universe=Vanilla
+DESIRED_Sites="T2_US_UCSD,UAF"
RequestMemory = 2048
RequestCpus = 1
executable=executable.sh
transfer_executable=True
transfer_input_files=package.tar.gz
transfer_output_files = ""
+Owner = undefined
+project_Name = "cmssurfandturf"
log=test.log
output=test.out
error=test.err
notification=Never
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
x509userproxy=/tmp/x509up_u31617
use_x509userproxy = True
Requirements = (HAS_SINGULARITY=?=True)
+SingularityImage = "/cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-alma9:latest"

arguments=

queue

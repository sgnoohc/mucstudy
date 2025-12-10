## Setup

### Optional

To get notifications on when SLURM jobs finish running, add this to your `~/.bashrc`:
```
export SBATCH_MAIL_USER="$USER@ufl.edu" # or whatever email you want
export SBATCH_MAIL_TYPE="END,FAIL"
```
Don't forget to source it `. ~/.bashrc`

## In the repo directory

```
./install_all.sh
sbatch slurm/ddsim.sh # to get MuMuToZH_sim.slcio

# after until the ddsim job is done and you got MuMuToZH_sim.slcio
sbatch slurm/k4run.sh # to get MuMuToZH_reco.slcio
```

## To view the output

```
apptainer shell /cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-alma9:latest
```
In the `Singularity>` shell
```
mucoll_setup
anajob MuMuToZH_sim.slcio
anajob MuMuToZH_reco.slcio
```

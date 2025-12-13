## Setup

In the repo directory
```
./install_all.sh
sbatch --mail-user="$USER@ufl.edu" --mail-type=END,FAIL --array=0-3 slurm.sh
```

## To view the output

```
apptainer run /cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-alma9:latest
```
In the `Singularity>` shell
```
mucoll_setup
cd slcio
anajob filename.slcio
```

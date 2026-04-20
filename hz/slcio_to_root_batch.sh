#!/usr/bin/env bash
JOBID="$1"
mkdir -p root

for i in slcio/MuMuToZH_reco_"$JOBID"_*.slcio; do
  out="root/$(basename "${i%.slcio}").root"
  python slcio_to_root.py -i "$in" -o "$out"
done

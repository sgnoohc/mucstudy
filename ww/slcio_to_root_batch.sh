#!/usr/bin/env bash
mkdir -p root

for i in slcio/MuMuToZH_reco_"$JOBID"_*.slcio; do
  out="$ODIR/root/$(basename "${i%.slcio}").root"
  python slcio_to_root.py -i "$in" -o "$out"
done

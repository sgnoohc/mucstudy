import os
import sys
from math import pi
from optparse import OptionParser

import pyLCIO as lcio
import ROOT

PDG_TO_NAME = {
    11: "e-",
    -11: "e+",
    13: "mu-",
    -13: "mu+",
    22: "gamma",
    12: "nu_e",
    -12: "nu_ebar",
    14: "nu_mu",
    -14: "nu_mubar",
    23: "Z0",
    24: "W+",
    -24: "W-",
    25: "H0",
    5: "b",
    -5: "bbar",
    6: "t",
    -6: "tbar",
}
NAME_TO_PDG = {name: pdg for pdg, name in PDG_TO_NAME.items()}


# Enable ROOT's automatic C++ STL vector handling
ROOT.gInterpreter.Declare("#include <vector>")


def get_PxPyPzE(mcp):
    p = mcp.getMomentum()
    px, py, pz = p[0], p[1], p[2]
    e = mcp.getEnergy()
    return ROOT.Math.PxPyPzEVector(px, py, pz, e)


def make_TH1D(name, nb, lo, hi, title=None):
    h = ROOT.TH1D(name, title or name, nb, lo, hi)
    h.SetDirectory(0)
    return h


# parse parameters
parser = OptionParser()
parser.add_option(
    "-i",
    "--inFile",
    help="--inFile output_reco.slcio",
    type=str,
    #default="output_reco.slcio",
)
parser.add_option(
    "-o",
    "--outFile",
    help="--outFile jet_study_output.root",
    type=str,
    default="jet_study_output.root",
)
(options, args) = parser.parse_args()

# create an LCIO reader and open an LCIO file
reader = lcio.IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# make higgs histograms
# NOTE: the histograms are empty and will be filled later
TH1D_PARAMS = {
    "pt": (180, 0, 6000),
    "eta": (180, -5, 5),
    "phi": (180, -pi, pi),
    "mass": (180, 0, 150),
    "energy": (180, 0, 6000),
}
# NOTE: *val unpacks a tuple to pass every item in it separately into make_TH1D()
hist_dict = {
    f"h_higgs_{key}": make_TH1D(f"h_higgs_{key}", *val)
    for key, val in TH1D_PARAMS.items()
}
hist_dict = hist_dict | {
    f"h_z_boson_{key}": make_TH1D(f"h_z_boson_{key}", *val)
    for key, val in TH1D_PARAMS.items()
}

all_PDGs = {}
for event in reader:
    mcps = event.getCollection("MCParticle")

    higgs_mcp = None
    z_boson_mcp = None
    for mcp in mcps:
        all_PDGs[mcp.getPDG()] = all_PDGs.get(mcp.getPDG(), 0) + 1
        if mcp.getPDG() == NAME_TO_PDG["H0"]:
            higgs_mcp = mcp
        if mcp.getPDG() == NAME_TO_PDG["Z0"]:
            z_boson_mcp = mcp

    if higgs_mcp is None or z_boson_mcp is None:
        continue
    higgs_PxPyPzE = get_PxPyPzE(higgs_mcp)
    z_boson_PxPyPzE = get_PxPyPzE(z_boson_mcp)

    # NOTE: will fail if you haven't created a histogram inside `hist_dict` beforehand
    hist_dict["h_higgs_pt"].Fill(higgs_PxPyPzE.Pt())
    hist_dict["h_higgs_eta"].Fill(higgs_PxPyPzE.Eta())
    hist_dict["h_higgs_phi"].Fill(higgs_PxPyPzE.Phi())
    hist_dict["h_higgs_mass"].Fill(higgs_PxPyPzE.M())
    hist_dict["h_higgs_energy"].Fill(higgs_PxPyPzE.E())
    hist_dict["h_z_boson_pt"].Fill(z_boson_PxPyPzE.Pt())
    hist_dict["h_z_boson_eta"].Fill(z_boson_PxPyPzE.Eta())
    hist_dict["h_z_boson_phi"].Fill(z_boson_PxPyPzE.Phi())
    hist_dict["h_z_boson_mass"].Fill(z_boson_PxPyPzE.M())
    hist_dict["h_z_boson_energy"].Fill(z_boson_PxPyPzE.E())
reader.close()


pdg_items = sorted(all_PDGs.items(), key=lambda kv: kv[1], reverse=True)
hist_dict["h_all_PDGs"] = ROOT.TH1I("h_all_PDGs", "PDG ID counts;PDG ID;Count", len(pdg_items), 0.5, len(pdg_items) + 0.5)
for i, (pdg, n) in enumerate(pdg_items, start=1):
    hist_dict["h_all_PDGs"].SetBinContent(i, n)
    hist_dict["h_all_PDGs"].GetXaxis().SetBinLabel(i, str(pdg))
hist_dict["h_all_PDGs"].LabelsOption("v", "X") # vertical labels
hist_dict["h_all_PDGs"].GetXaxis().SetLabelSize(0.03)


# write histograms to .root file
with ROOT.TFile(options.outFile, "recreate") as outfile:
    for name, hist in hist_dict.items():
        outfile.WriteObject(hist, name)

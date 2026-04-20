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
)
parser.add_option(
    "-o",
    "--outFile",
    help="--outFile jet_study_output.root",
    type=str,
)
(options, args) = parser.parse_args()

# create an LCIO reader and open an LCIO file
reader = lcio.IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# make higgs histograms
# NOTE: the histograms are empty and will be filled later
TH1D_PARAMS = {
    #"pt": (180, 0, 6000),
    #"eta": (180, -5, 5),
    #"phi": (180, -pi, pi),
    #"mass": (180, 0, 150),
    "energy": (180, 0, 6000),
}
hist_dict = {
    f"wp_boson_{key}": make_TH1D(f"wp_boson_{key}", *val)
    for key, val in TH1D_PARAMS.items()
}
hist_dict = {
    f"wn_boson_{key}": make_TH1D(f"wn_boson_{key}", *val)
    for key, val in TH1D_PARAMS.items()
}
# NOTE: *val unpacks a tuple to pass every item in it separately into make_TH1D()
hist_dict = {
    f"wp_boson_{key}": make_TH1D(f"wp_boson_{key}", *val)
    for key, val in TH1D_PARAMS.items()
}
hist_dict = hist_dict | {
    f"wn_boson_{key}": make_TH1D(f"wn_boson_{key}", *val)
    for key, val in TH1D_PARAMS.items()
}

all_PDGs = {}
for event in reader:
    mcps = event.getCollection("MCParticle")

    wp_boson_mcp = None # W positive boson
    wn_boson_mcp = None # W negative boson
    for mcp in mcps:
        if mcp.getPDG() == NAME_TO_PDG["W+"]:
            wp_boson_mcp = mcp
        if mcp.getPDG() == NAME_TO_PDG["W-"]:
            wn_boson_mcp = mcp

    if wp_boson_mcp is None or wn_boson_mcp is None:
        continue
    wp_boson_PxPyPzE = get_PxPyPzE(wp_boson_mcp)
    wn_boson_PxPyPzE = get_PxPyPzE(wn_boson_mcp)

    # NOTE: will fail if you haven't created a histogram inside `hist_dict` beforehand
    hist_dict["wp_boson_energy"].Fill(wp_boson_PxPyPzE.E())
    hist_dict["wn_boson_energy"].Fill(wn_boson_PxPyPzE.E())
reader.close()

# write histograms to .root file
print("hist outFile =", options.outFile)
with ROOT.TFile(options.outFile, "recreate") as outfile:
    for name, hist in hist_dict.items():
        outfile.WriteObject(hist, name)

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
    default="output_reco.slcio",
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
# NOTE: they are empty and will be filled later
HIGGS_1D_PARAMS = {
    "pt": (180, 0, 6000),
    "eta": (180, -5, 5),
    "phi": (180, -pi, pi),
    "mass": (180, 0, 150),
    "energy": (180, 0, 6000),
}
hist_dict = {
    # NOTE: *val unpacks a tuple to pass every item in it separately into book1d()
    f"h_higgs_{key}": make_TH1D(f"h_higgs_{key}", *val)
    for key, val in HIGGS_1D_PARAMS.items()
}

for ievent, event in enumerate(reader):
    mcps = event.getCollection("MCParticle")

    # gets first higgs mcp
    higgs_mcp = next((mcp for mcp in mcps if mcp.getPDG() == NAME_TO_PDG["H0"]), None)

    if higgs_mcp is None:
        continue
    higgs_PxPyPzE = get_PxPyPzE(higgs_mcp)

    # NOTE: will fail if you haven't created a histogram inside `hist_dict` beforehand
    hist_dict["h_higgs_pt"].Fill(higgs_PxPyPzE.Pt())
    hist_dict["h_higgs_eta"].Fill(higgs_PxPyPzE.Eta())
    hist_dict["h_higgs_phi"].Fill(higgs_PxPyPzE.Phi())
    hist_dict["h_higgs_mass"].Fill(higgs_PxPyPzE.M())
    hist_dict["h_higgs_energy"].Fill(higgs_PxPyPzE.E())


reader.close()
with ROOT.TFile(options.outFile, "recreate") as outfile:
    for hist in hist_dict.values():
        outfile.WriteObject(hist, "myhisto")

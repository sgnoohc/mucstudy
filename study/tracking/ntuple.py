from array import array
import os, sys
from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1D, TH2D, TFile, TLorentzVector, TTree, TMath, TVector3, std, gInterpreter
from math import *
from optparse import OptionParser

# Enable ROOT's automatic C++ STL vector handling
gInterpreter.Declare("#include <vector>")

# VERBOSE=True
VERBOSE=False
STOPEVENT=10

def vprint(*args, **kwargs):
    if VERBOSE:
        print(*args, **kwargs)

#########################
# parameters

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile ntuple.root',
                  type=str, default='ntuple.root')
(options, args) = parser.parse_args()

tree = TTree("tree", "tree")

## Create branches
col_i = {}
col_i["nhit"] = array('i', [0])
col_i["nmcp"] = array('i', [0])
for i_n in col_i:
    tree.Branch(i_n, col_i[i_n], f"{i_n}/I")
col_vf = {}
col_vf["hit_x"] = std.vector('float')()
col_vf["hit_y"] = std.vector('float')()
col_vf["hit_z"] = std.vector('float')()
col_vf["hit_system"] = std.vector('float')()
col_vf["hit_layer"] = std.vector('float')()
col_vf["hit_side"] = std.vector('float')()
col_vf["hit_module"] = std.vector('float')()
col_vf["hit_sensor"] = std.vector('float')()
col_vf["hit_mcp_id"] = std.vector('float')()
col_vf["mcp_id"] = std.vector('float')()
col_vf["mcp_pt"] = std.vector('float')()
col_vf["mcp_eta"] = std.vector('float')()
col_vf["mcp_phi"] = std.vector('float')()
col_vf["mcp_pdgId"] = std.vector('float')()
col_vf["mcp_q"] = std.vector('float')()
col_vf["mcp_isDenom"] = std.vector('float')()
for vf_n in col_vf:
    tree.Branch(vf_n, col_vf[vf_n])

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

tracker_systems = [
    "VB",
    "VE",
    "IB",
    "IE",
    "OB",
    "OE",
]

tracker_systems_map = {
    1: "VB",
    2: "VE",
    3: "IB",
    4: "IE",
    5: "OB",
    6: "OE",
}

# Note: VX first layer
# VXBarrel
# Layer 0-1 = 0-15 modules, 0-4 sensors
# Layer 2-3 = 0-14 modules, 0-4 sensors
# Layer 4-5 = 0-20 modules, 0-4 sensors
# Layer 6-7 = 0-28 modules, 0-4 sensors
# VXEndcap
# Layer ALL = 0 modules, 0-15 sensors, -1 or 1 side

def dPhiThreshold(lhit_v3):
    Bfield = 5
    kRinv1GeVf = (2.99792458e-3 * Bfield)
    k2Rinv1GeVf = kRinv1GeVf / 2.
    sinAlphaMax = 0.95
    ptCut = 1
    rt = lhit_v3.Perp()
    miniSlope = asin(min(rt * k2Rinv1GeVf / ptCut, sinAlphaMax))
    error = 0
    return miniSlope + sqrt(error**2)

# loop over all events in the file
for ievent, event in enumerate(reader):

    if ievent % 100 == 0:
        print("Processing event " + str(ievent))

    # Clear all the branches
    for i_n in col_i:
        col_i[i_n][0] = -999
    for vf_n in col_vf:
        col_vf[vf_n].clear()

    hitCollections = {}
    hitRelationCollections = {}
    hitRelations = {}

    for tracker_system in tracker_systems:
        print(f"tracker_system: {tracker_system} , ")
        hitCollections[f"{tracker_system}"] = event.getCollection(f"{tracker_system}TrackerHits")
        hitRelationCollections[f"{tracker_system}"] = event.getCollection(f"{tracker_system}TrackerHitsRelations")
        hitRelations[f"{tracker_system}"] = UTIL.LCRelationNavigator(hitRelationCollections[f"{tracker_system}"])
        encoding = hitCollections[f"{tracker_system}"].getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)
        for ihit, hit in enumerate(hitCollections[f"{tracker_system}"]):
            if ihit % 10000 == 0:
                print(f"ihit: {ihit} , ")
            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            system = decoder["system"].value()
            layer = decoder["layer"].value()
            side = decoder["side"].value()
            module = decoder["module"].value()
            sensor = decoder["sensor"].value()
            mcp_id = -999
            simhits = hitRelations[f"{tracker_systems_map[system]}"].getRelatedToObjects(hit)
            if len(simhits) > 0:
                mcp = simhits[0].getMCParticle()
                if mcp:
                    mcp_id = mcp.id()
            col_vf["hit_x"].push_back(hit.getPosition()[0])
            col_vf["hit_y"].push_back(hit.getPosition()[1])
            col_vf["hit_z"].push_back(hit.getPosition()[2])
            col_vf["hit_system"].push_back(system)
            col_vf["hit_layer"].push_back(layer)
            col_vf["hit_side"].push_back(side)
            col_vf["hit_module"].push_back(module)
            col_vf["hit_sensor"].push_back(sensor)
            col_vf["hit_mcp_id"].push_back(mcp_id)

    col_i["nhit"][0] = col_vf["hit_x"].size()

    # get all MC particle that I want to track
    mcps = event.getCollection("MCParticle")
    for mcp in mcps:
        mcp_id = mcp.id()
        mcp_px = mcp.getMomentum()[0]
        mcp_py = mcp.getMomentum()[1]
        mcp_pz = mcp.getMomentum()[2]
        mcp_energy = mcp.getEnergy()
        tlv = TLorentzVector()
        tlv.SetPxPyPzE(mcp_px, mcp_py, mcp_pz, mcp_energy)
        pt = tlv.Pt()
        eta = tlv.Eta()
        phi = tlv.Phi()
        q = mcp.getCharge()
        pdgid = mcp.getPDG()
        isdenom = pt > 1 and abs(eta) < 1.8
        vprint(f"mcp_id={mcp_id} pt={pt:.3f} eta={eta:.3f} phi={phi:.3f} pdgid={pdgid}")

        col_vf["mcp_id"].push_back(mcp_id)
        col_vf["mcp_pt"].push_back(pt)
        col_vf["mcp_eta"].push_back(eta)
        col_vf["mcp_phi"].push_back(phi)
        col_vf["mcp_pdgId"].push_back(pdgid)
        col_vf["mcp_q"].push_back(q)
        col_vf["mcp_isDenom"].push_back(isdenom)

    col_i["nmcp"][0] = col_vf["mcp_id"].size()

    tree.Fill()

    if VERBOSE and ievent == STOPEVENT:
        break

reader.close()
output_file = TFile(options.outFile, 'RECREATE')
tree.Write()
output_file.Close()

from array import array
import os, sys
from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1D, TH2D, TFile, TLorentzVector, TTree, TMath, TVector3, std, gInterpreter
from math import *
from optparse import OptionParser
from glob import glob

# Enable ROOT's automatic C++ STL vector handling
gInterpreter.Declare("#include <vector>")

# VERBOSE=True
VERBOSE=False
STOPEVENT=10

#########################
# parameters

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile mcp_scan_output.root',
                  type=str, default='jet_study_output.root')
parser.add_option('-e', '--event', help='--event 0',
                  type=int, default=0)
(options, args) = parser.parse_args()

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

#########################
h_higgs_pt = TH1D("h_higgs_pt", "h_higgs_pt", 180, 0, 6000)
h_higgs_eta = TH1D("h_higgs_eta", "h_higgs_eta", 180, -5, 5)
h_higgs_phi = TH1D("h_higgs_phi", "h_higgs_phi", 180, -3.1416, 3.1416)
h_higgs_mass = TH1D("h_higgs_mass", "h_higgs_mass", 180, 0, 150)
h_higgs_energy = TH1D("h_higgs_energy", "h_higgs_energy", 180, 0, 6000)
h_zboson_pt = TH1D("h_zboson_pt", "h_zboson_pt", 180, 0, 6000)
h_zboson_eta = TH1D("h_zboson_eta", "h_zboson_eta", 180, -5, 5)
h_zboson_phi = TH1D("h_zboson_phi", "h_zboson_phi", 180, -3.1416, 3.1416)
h_zboson_mass = TH1D("h_zboson_mass", "h_zboson_mass", 180, 0, 150)
h_zboson_energy = TH1D("h_zboson_energy", "h_zboson_energy", 180, 0, 6000)

h_reco_higgs_pt = TH1D("h_reco_higgs_pt", "h_reco_higgs_pt", 180, 0, 6000)
h_reco_higgs_eta = TH1D("h_reco_higgs_eta", "h_reco_higgs_eta", 180, -5, 5)
h_reco_higgs_phi = TH1D("h_reco_higgs_phi", "h_reco_higgs_phi", 180, -3.1416, 3.1416)
h_reco_higgs_mass = TH1D("h_reco_higgs_mass", "h_reco_higgs_mass", 180, 0, 150)
h_reco_higgs_energy = TH1D("h_reco_higgs_energy", "h_reco_higgs_energy", 180, 0, 6000)
h_reco_zboson_pt = TH1D("h_reco_zboson_pt", "h_reco_zboson_pt", 180, 0, 6000)
h_reco_zboson_eta = TH1D("h_reco_zboson_eta", "h_reco_zboson_eta", 180, -5, 5)
h_reco_zboson_phi = TH1D("h_reco_zboson_phi", "h_reco_zboson_phi", 180, -3.1416, 3.1416)
h_reco_zboson_mass = TH1D("h_reco_zboson_mass", "h_reco_zboson_mass", 180, 0, 150)
h_reco_zboson_energy = TH1D("h_reco_zboson_energy", "h_reco_zboson_energy", 180, 0, 6000)

h_higgs_pt_reco_v_truth = TH2D("h_higgs_pt_reco_v_truth", "h_higgs_pt_reco_v_truth", 50, 0, 6000, 50, 0, 6000)
h_zboson_pt_reco_v_truth = TH2D("h_zboson_pt_reco_v_truth", "h_zboson_pt_reco_v_truth", 50, 0, 6000, 50, 0, 6000)
h_higgs_energy_reco_v_truth = TH2D("h_higgs_energy_reco_v_truth", "h_higgs_energy_reco_v_truth", 50, 0, 6000, 50, 0, 6000)
h_zboson_energy_reco_v_truth = TH2D("h_zboson_energy_reco_v_truth", "h_zboson_energy_reco_v_truth", 50, 0, 6000, 50, 0, 6000)

nbins_2d = 520

cal_hit_rz = TH2D('CAL_rz', 'CAL_rz',  nbins_2d, -5000, 5000,  nbins_2d, 0, 5000)
ecal_hit_rz = TH2D('ECAL_rz', 'ECAL_rz',  nbins_2d, -5000, 5000,  nbins_2d, 0, 5000)
hcal_hit_rz = TH2D('HCAL_rz', 'HCAL_rz',  nbins_2d, -5000, 5000,  nbins_2d, 0, 5000)

cal_hit_xy = TH2D('CAL_xy', 'CAL_xy',  nbins_2d, -5000, 5000,  nbins_2d, -5000, 5000)
ecal_hit_xy = TH2D('ECAL_xy', 'ECAL_xy',  nbins_2d, -5000, 5000,  nbins_2d, -5000, 5000)
hcal_hit_xy = TH2D('HCAL_xy', 'HCAL_xy',  nbins_2d, -5000, 5000,  nbins_2d, -5000, 5000)

h_total_E = TH1D("h_total_E", "h_total_E", 180, 0, 12000)
h_total_E_rechit_v_mcp = TH2D("h_total_E_rechit_v_mcp", "h_total_E_rechit_v_mcp", 180, 0, 12000, 180, 0, 12000)

h_total_mcp_E = TH1D("h_total_mcp_E", "h_total_mcp_E", 180, 0, 12000)
h_total_mcp_E_notleft = TH1D("h_total_mcp_E_notleft", "h_total_mcp_E_notleft", 180, 0, 12000)

h_total_mcp_E_v_notleft = TH2D("h_total_mcp_E_v_notleft", "h_total_mcp_E_v_notleft", 180, 0, 12000, 180, 0, 12000)


hists = [h_higgs_pt, h_higgs_eta, h_higgs_phi, h_higgs_mass, h_higgs_energy,
         h_zboson_pt, h_zboson_eta, h_zboson_phi, h_zboson_mass, h_zboson_energy,
         h_higgs_pt_reco_v_truth, h_higgs_energy_reco_v_truth,
         h_zboson_pt_reco_v_truth, h_zboson_energy_reco_v_truth,
         h_reco_higgs_pt, h_reco_higgs_eta, h_reco_higgs_phi, h_reco_higgs_mass, h_reco_higgs_energy,
         h_reco_zboson_pt, h_reco_zboson_eta, h_reco_zboson_phi, h_reco_zboson_mass, h_reco_zboson_energy,
         h_total_E, h_total_E_rechit_v_mcp,
         cal_hit_rz, ecal_hit_rz, hcal_hit_rz,
         cal_hit_xy, ecal_hit_xy, hcal_hit_xy,
         h_total_mcp_E, h_total_mcp_E_notleft, h_total_mcp_E_v_notleft,
        ]

#_________________________________________________________________________________________________________________________________________
def get_tlv(mcp):
    mcp_px = mcp.getMomentum()[0]
    mcp_py = mcp.getMomentum()[1]
    mcp_pz = mcp.getMomentum()[2]
    mcp_energy = mcp.getEnergy()
    tlv = TLorentzVector()
    tlv.SetPxPyPzE(mcp_px, mcp_py, mcp_pz, mcp_energy)
    return tlv

#_________________________________________________________________________________________________________________________________________
def get_v3_hit(hit):
    v3 = TVector3()
    hit_x = hit.getPosition()[0]
    hit_y = hit.getPosition()[1]
    hit_z = hit.getPosition()[2]
    v3.SetXYZ(hit_x, hit_y, hit_z)
    return v3

#_________________________________________________________________________________________________________________________________________
def delta_r(hit, mcp):
    x = mcp.X()
    y = mcp.Y()
    z = mcp.Z()
    v3_mcp = TVector3()
    v3_mcp.SetXYZ(x, y, z)
    return v3_mcp.DeltaR(hit)

#_________________________________________________________________________________________________________________________________________
def print_daughters(particle, depth=0):
    indent = "  " * depth
    if particle.getGeneratorStatus() == 1:
        print(f"particle.getPDG(): {particle.getPDG()} , particle.getGeneratorStatus(): {particle.getGeneratorStatus()} , particle.getEnergy(): {particle.getEnergy()} , ")
    # print(f"{indent}PDG: {particle.getPDG()}, Energy: {particle.getEnergy():.2f} GeV")
    for daughter in particle.getDaughters():
        print_daughters(daughter, depth + 1)

#_________________________________________________________________________________________________________________________________________
def sum_daughter_tlv(particle):
    total_tlv = TLorentzVector()
    daughters = particle.getDaughters()
    for da in daughters:
        px = da.getMomentum()[0]
        py = da.getMomentum()[1]
        pz = da.getMomentum()[2]
        energy = da.getEnergy()
        status = da.getGeneratorStatus()
        da_tlv = TLorentzVector()
        da_tlv.SetPxPyPzE(px, py, pz, energy)
        if status == 0:
            total_tlv += da_tlv
        # Recursively add daughters of this daughter
        total_tlv += sum_daughter_tlv(da)
    return total_tlv

#_________________________________________________________________________________________________________________________________________
def print_content_relations(event):
    for collection_name in event.getCollectionNames():
        col = event.getCollection(collection_name)
        coltype = str(col.getTypeName())
        print(f"{str(collection_name):50} {coltype:50}")
        if coltype == "LCRelation":
            relation = UTIL.LCRelationNavigator(col)
            print(f"{relation.getFromType()} -> {relation.getToType()}")

#_________________________________________________________________________________________________________________________________________
def find_daughters_recursive(particle, pdg_id=None, status_code=None):
    """
    Recursively searches daughters of an MCParticle for matches by PDG ID and status code.
    Args:
        particle: MCParticle object (from pyLCIO).
        pdg_id: int or None. If set, filters by matching PDG ID.
        status_code: int or None. If set, filters by matching status code.
    Returns:
        List of matching daughter MCParticles.
    """
    matches = []
    def recurse(p):
        for daughter in p.getDaughters():
            match = True
            if pdg_id is not None and daughter.getPDG() != pdg_id:
                match = False
            if status_code is not None and daughter.getGeneratorStatus() != status_code:
                match = False
            if match:
                matches.append(daughter)
            recurse(daughter)
    recurse(particle)
    return matches

#_________________________________________________________________________________________________________________________________________
def find_first_daughter_recursive(particle, pdg_id=None, status_code=None):
    """
    Recursively searches daughters of an MCParticle and returns the first match
    by PDG ID and status code.
    Args:
        particle: MCParticle object (from pyLCIO).
        pdg_id: int or None. If set, filters by matching PDG ID.
        status_code: int or None. If set, filters by matching status code.
    Returns:
        The first matching daughter MCParticle, or None if no match is found.
    """
    def recurse(p):
        for daughter in p.getDaughters():
            match = True
            if pdg_id is not None and daughter.getPDG() != pdg_id:
                match = False
            if status_code is not None and daughter.getGeneratorStatus() != status_code:
                match = False
            if match:
                return daughter
            result = recurse(daughter)
            if result:
                return result
        return None
    return recurse(particle)


#_________________________________________________________________________________________________________________________________________
def print_mcparticle(mcp):
    def get_pdg_name(pdgid):
        pdg_names = {
            11: "e-", -11: "e+",
            13: "mu-", -13: "mu+",
            22: "gamma", 12: "nu_e", -12: "nu_ebar",
            14: "nu_mu", -14: "nu_mubar",
            23: "Z0", 24: "W+", -24: "W-",
            25: "H0", 5: "b", -5: "bbar", 6: "t", -6: "tbar"
        }
        return pdg_names.get(pdgid, f"({pdgid})")
    # header = (
    #     " i  id    name         stat mothers(d1,d2) daughters(d1,d2)"
    #     "       px       py       pz        E      vx      vy      vz"
    #     "     endx    endy    endz   left?"
    # )
    # print(header)
    # print("-" * len(header))
    p = mcp
    pdg = p.getPDG()
    status = p.getGeneratorStatus()
    name = get_pdg_name(pdg).ljust(13)
    momentum = p.getMomentum()
    tlv = get_tlv(mcp)
    pt = tlv.Pt()
    eta = tlv.Eta()
    phi = tlv.Phi()
    mass = tlv.M()
    energy = p.getEnergy()
    vertex = p.getVertex()
    endpoint = p.getEndpoint()
    left_detector = p.hasLeftDetector() if hasattr(p, "hasLeftDetector") else False
    print(f" {pdg:5} {name} {status:5}   "
          f"  {pt:12.2f} {eta:12.2f} {phi:12.2f} {mass:8.2f} {energy:12.2f}"
          f" {vertex[0]:6.2f} {vertex[1]:6.2f} {vertex[2]:6.2f}"
          f" {endpoint[0]:12.2f} {endpoint[1]:12.2f} {endpoint[2]:12.2f}"
          f"   {'yes' if left_detector else 'no ':>3}")

for ievent, event in enumerate(reader):

    # get collections
    mcps = event.getCollection("MCParticle")

    # get the Z and H and bb bb
    h = None
    z = None
    b_h = None
    bx_h = None
    b_z = None
    bx_z = None
    for mcp in mcps:
        status = mcp.getGeneratorStatus()
        pdgid = mcp.getPDG()
        if pdgid == 25 and status == 22:
            h = mcp
            b_h = find_first_daughter_recursive(h, pdg_id=5, status_code=23)
            bx_h = find_first_daughter_recursive(h, pdg_id=-5, status_code=23)
        if pdgid == 23 and status == 22:
            z = mcp
            b_z = find_first_daughter_recursive(z, pdg_id=5, status_code=23)
            bx_z = find_first_daughter_recursive(z, pdg_id=-5, status_code=23)

    print_mcparticle(h)
    print_mcparticle(z)
    print_mcparticle(b_h)
    print_mcparticle(bx_h)
    print_mcparticle(b_z)
    print_mcparticle(bx_z)

    tlv_h = get_tlv(h)
    tlv_z = get_tlv(z)
    h_higgs_pt.Fill(tlv_h.Pt())
    h_higgs_eta.Fill(tlv_h.Eta())
    h_higgs_phi.Fill(tlv_h.Phi())
    h_higgs_mass.Fill(tlv_h.M())
    h_higgs_energy.Fill(tlv_h.E())
    h_zboson_pt.Fill(tlv_z.Pt())
    h_zboson_eta.Fill(tlv_z.Eta())
    h_zboson_phi.Fill(tlv_z.Phi())
    h_zboson_mass.Fill(tlv_z.M())
    h_zboson_energy.Fill(tlv_z.E())

    # get all MC particle that I want to track
    mcp_status_1s = []
    mcp_status_1s_notleft = []
    for mcp in mcps:
        if mcp.getGeneratorStatus() == 1:
            mcp_status_1s.append(mcp)
            if not mcp.hasLeftDetector():
                mcp_status_1s_notleft.append(mcp)

    mcp_total_energy = sum(p.getEnergy() for p in mcp_status_1s)
    mcp_total_energy_notleft = sum(p.getEnergy() for p in mcp_status_1s_notleft)

    h_total_mcp_E.Fill(mcp_total_energy)
    h_total_mcp_E_notleft.Fill(mcp_total_energy_notleft)
    h_total_mcp_E_v_notleft.Fill(mcp_total_energy, mcp_total_energy_notleft)

    print(f"mcp_total_energy: {mcp_total_energy} , ")
    print(f"mcp_total_energy_notleft: {mcp_total_energy_notleft} , ")

    mcp_matched_to_h = []
    mcp_matched_to_z = []
    for mcp in mcp_status_1s_notleft:
        tlv_mcp = get_tlv(mcp)
        dr = tlv_h.DeltaR(tlv_mcp)
        if dr < 0.8:
            mcp_matched_to_h.append(mcp)
        dr = tlv_z.DeltaR(tlv_mcp)
        if dr < 0.8:
            mcp_matched_to_z.append(mcp)

    mcp_total_energy_h = sum(p.getEnergy() for p in mcp_matched_to_h)
    mcp_total_energy_z = sum(p.getEnergy() for p in mcp_matched_to_z)
    print(f"mcp_total_energy_h: {mcp_total_energy_h} , ")
    print(f"mcp_total_energy_z: {mcp_total_energy_z} , ")

    hit_matched_to_h = []
    hit_matched_to_z = []
    #-------------------------------------------------------------------
    total_E = 0
    # magic_scaling = 46.6
    magic_scaling = 1
    ECALrechitCollection = event.getCollection('EcalBarrelCollectionRec')
    for rechit in ECALrechitCollection:
        v3_hit = get_v3_hit(rechit)
        dr = delta_r(v3_hit, tlv_h)
        if dr < 0.8:
            hit_matched_to_h.append(rechit)
        dr = delta_r(v3_hit, tlv_z)
        if dr < 0.8:
            hit_matched_to_z.append(rechit)
        # print(f"pt: {pt} , eta: {eta} , phi: {phi} , mass: {mass} , ")
        total_E += rechit.getEnergy() * magic_scaling
    ECALrechitCollection = event.getCollection('EcalEndcapCollectionRec')
    for rechit in ECALrechitCollection:
        v3_hit = get_v3_hit(rechit)
        dr = delta_r(v3_hit, tlv_h)
        if dr < 0.8:
            hit_matched_to_h.append(rechit)
        dr = delta_r(v3_hit, tlv_z)
        if dr < 0.8:
            hit_matched_to_z.append(rechit)
        total_E += rechit.getEnergy() * magic_scaling
    HCALrechitCollection = event.getCollection('HcalBarrelCollectionRec')
    for rechit in HCALrechitCollection:
        v3_hit = get_v3_hit(rechit)
        dr = delta_r(v3_hit, tlv_h)
        if dr < 0.8:
            hit_matched_to_h.append(rechit)
        dr = delta_r(v3_hit, tlv_z)
        if dr < 0.8:
            hit_matched_to_z.append(rechit)
        total_E += rechit.getEnergy() * magic_scaling
    HCALrechitCollection = event.getCollection('HcalEndcapCollectionRec')
    for rechit in HCALrechitCollection:
        v3_hit = get_v3_hit(rechit)
        dr = delta_r(v3_hit, tlv_h)
        if dr < 0.8:
            hit_matched_to_h.append(rechit)
        dr = delta_r(v3_hit, tlv_z)
        if dr < 0.8:
            hit_matched_to_z.append(rechit)
        total_E += rechit.getEnergy() * magic_scaling
    print(f"total_E: {total_E} , ")
    hit_matched_energy_h = sum(p.getEnergy() for p in hit_matched_to_h)
    hit_matched_energy_z = sum(p.getEnergy() for p in hit_matched_to_z)
    print(f"hit_matched_energy_h: {hit_matched_energy_h} , ")
    print(f"hit_matched_energy_z: {hit_matched_energy_z} , ")
    h_total_E.Fill(total_E)
    h_total_E_rechit_v_mcp.Fill(total_E, mcp_total_energy_notleft)
    #-------------------------------------------------------------------

    simhit_matched_to_h = []
    simhit_matched_to_z = []
    #-------------------------------------------------------------------
    total_E = 0
    # magic_scaling = 46.6
    magic_scaling = 1
    ECALsimhitCollection = event.getCollection('ECalBarrelCollection')
    for simhit in ECALsimhitCollection:
        v3_hit = get_v3_hit(simhit)
        dr = delta_r(v3_hit, tlv_h)
        if dr < 0.8:
            simhit_matched_to_h.append(simhit)
        dr = delta_r(v3_hit, tlv_z)
        if dr < 0.8:
            simhit_matched_to_z.append(simhit)
        # print(f"pt: {pt} , eta: {eta} , phi: {phi} , mass: {mass} , ")
        total_E += simhit.getEnergy() * magic_scaling
    ECALsimhitCollection = event.getCollection('ECalEndcapCollection')
    for simhit in ECALsimhitCollection:
        v3_hit = get_v3_hit(simhit)
        dr = delta_r(v3_hit, tlv_h)
        if dr < 0.8:
            simhit_matched_to_h.append(simhit)
        dr = delta_r(v3_hit, tlv_z)
        if dr < 0.8:
            simhit_matched_to_z.append(simhit)
        total_E += simhit.getEnergy() * magic_scaling
    HCALsimhitCollection = event.getCollection('HCalBarrelCollection')
    for simhit in HCALsimhitCollection:
        v3_hit = get_v3_hit(simhit)
        dr = delta_r(v3_hit, tlv_h)
        if dr < 0.8:
            simhit_matched_to_h.append(simhit)
        dr = delta_r(v3_hit, tlv_z)
        if dr < 0.8:
            simhit_matched_to_z.append(simhit)
        total_E += simhit.getEnergy() * magic_scaling
    HCALsimhitCollection = event.getCollection('HCalEndcapCollection')
    for simhit in HCALsimhitCollection:
        v3_hit = get_v3_hit(simhit)
        dr = delta_r(v3_hit, tlv_h)
        if dr < 0.8:
            simhit_matched_to_h.append(simhit)
        dr = delta_r(v3_hit, tlv_z)
        if dr < 0.8:
            simhit_matched_to_z.append(simhit)
        total_E += simhit.getEnergy() * magic_scaling
    print(f"total_E: {total_E} , ")
    simhit_matched_energy_h = sum(p.getEnergy() for p in simhit_matched_to_h)
    simhit_matched_energy_z = sum(p.getEnergy() for p in simhit_matched_to_z)
    print(f"simhit_matched_energy_h: {simhit_matched_energy_h} , ")
    print(f"simhit_matched_energy_z: {simhit_matched_energy_z} , ")

    targets = [ ("h", get_tlv(h)), ("z", get_tlv(z)), ]
    dr_threshold = 0.4

    print("here")
    pfo_matched_to_h = []
    pfo_matched_to_z = []
    pfosCollection = event.getCollection("PandoraPFOs")
    for ipfo, pfo in enumerate(pfosCollection):
        px = pfo.getMomentum()[0]
        py = pfo.getMomentum()[1]
        pz = pfo.getMomentum()[2]
        E = pfo.getEnergy()
        tlv_pfo = get_tlv(pfo)
        dr = tlv_h.DeltaR(tlv_pfo)
        if dr < 0.8:
            pfo_matched_to_h.append(pfo)
        dr = tlv_z.DeltaR(tlv_pfo)
        if dr < 0.8:
            pfo_matched_to_z.append(pfo)
        # print(f"pt: {pt} , eta: {eta} , phi: {phi} , mass: {mass} , ")
    pfo_total_energy_h = sum(p.getEnergy() for p in pfo_matched_to_h)
    pfo_total_energy_z = sum(p.getEnergy() for p in pfo_matched_to_z)
    print(f"pfo_total_energy_h: {pfo_total_energy_h} , ")
    print(f"pfo_total_energy_z: {pfo_total_energy_z} , ")

    # For each target, find the closest jet
    for name, target in targets:
        best_jet = None
        best_dr = float("inf")

        jets = event.getCollection("JetOut")
        for i, jet in enumerate(jets):
            tlv = get_tlv(jet)
            dr = tlv.DeltaR(target)
            if dr < best_dr:
                best_dr = dr
                best_jet = jet

        if best_dr < dr_threshold:
            print(f"{name} matched to Jet {best_jet} with dR = {best_dr:.3f} E={best_jet.getEnergy()}")
            if name == "h":
                jet = best_jet
                tlv = get_tlv(jet)
                h_higgs_pt_reco_v_truth.Fill(target.Pt(), tlv.Pt())
                h_higgs_energy_reco_v_truth.Fill(target.E(), tlv.E())
                h_reco_higgs_pt.Fill(tlv.Pt())
                h_reco_higgs_eta.Fill(tlv.Eta())
                h_reco_higgs_phi.Fill(tlv.Phi())
                h_reco_higgs_mass.Fill(tlv.M())
                h_reco_higgs_energy.Fill(tlv.E())
            elif name == "z":
                jet = best_jet
                tlv = get_tlv(jet)
                h_zboson_pt_reco_v_truth.Fill(target.Pt(), tlv.Pt())
                h_zboson_energy_reco_v_truth.Fill(target.E(), tlv.E())
                h_reco_zboson_pt.Fill(tlv.Pt())
                h_reco_zboson_eta.Fill(tlv.Eta())
                h_reco_zboson_phi.Fill(tlv.Phi())
                h_reco_zboson_mass.Fill(tlv.M())
                h_reco_zboson_energy.Fill(tlv.E())
        else:
            print(f"{name} has no match within dR < {dr_threshold}")

    # ECAL barrel
    hitsCollection = event.getCollection("EcalBarrelCollectionRec")
    encoding = hitsCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
    decoder = UTIL.BitField64(encoding)
    for ihit, hit in enumerate(hitsCollection):
        r = sqrt(hit.getPosition()[0]*hit.getPosition()
                 [0] + hit.getPosition()[1]*hit.getPosition()[1])
        ecal_hit_rz.Fill(hit.getPosition()[2], r, hit.getEnergy())
        cal_hit_rz.Fill(hit.getPosition()[2], r, hit.getEnergy())
        ecal_hit_xy.Fill(hit.getPosition()[0], hit.getPosition()[1], hit.getEnergy())
        cal_hit_xy.Fill(hit.getPosition()[0], hit.getPosition()[1], hit.getEnergy())

    # HCAL barrel
    hitsCollection = event.getCollection("HcalBarrelCollectionRec")
    encoding = hitsCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
    decoder = UTIL.BitField64(encoding)
    for ihit, hit in enumerate(hitsCollection):
        r = sqrt(hit.getPosition()[0]*hit.getPosition()
                 [0] + hit.getPosition()[1]*hit.getPosition()[1])
        hcal_hit_rz.Fill(hit.getPosition()[2], r, hit.getEnergy())
        cal_hit_rz.Fill(hit.getPosition()[2], r, hit.getEnergy())
        hcal_hit_xy.Fill(hit.getPosition()[0], hit.getPosition()[1], hit.getEnergy())
        cal_hit_xy.Fill(hit.getPosition()[0], hit.getPosition()[1], hit.getEnergy())

    # ECAL barrel
    hitsCollection = event.getCollection("EcalEndcapCollectionRec")
    encoding = hitsCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
    decoder = UTIL.BitField64(encoding)
    for ihit, hit in enumerate(hitsCollection):
        r = sqrt(hit.getPosition()[0]*hit.getPosition()
                 [0] + hit.getPosition()[1]*hit.getPosition()[1])
        ecal_hit_rz.Fill(hit.getPosition()[2], r, hit.getEnergy())
        cal_hit_rz.Fill(hit.getPosition()[2], r, hit.getEnergy())
        ecal_hit_xy.Fill(hit.getPosition()[0], hit.getPosition()[1], hit.getEnergy())
        cal_hit_xy.Fill(hit.getPosition()[0], hit.getPosition()[1], hit.getEnergy())

    # HCAL barrel
    hitsCollection = event.getCollection("HcalEndcapCollectionRec")
    encoding = hitsCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
    decoder = UTIL.BitField64(encoding)
    for ihit, hit in enumerate(hitsCollection):
        r = sqrt(hit.getPosition()[0]*hit.getPosition()
                 [0] + hit.getPosition()[1]*hit.getPosition()[1])
        hcal_hit_rz.Fill(hit.getPosition()[2], r, hit.getEnergy())
        cal_hit_rz.Fill(hit.getPosition()[2], r, hit.getEnergy())
        hcal_hit_xy.Fill(hit.getPosition()[0], hit.getPosition()[1], hit.getEnergy())
        cal_hit_xy.Fill(hit.getPosition()[0], hit.getPosition()[1], hit.getEnergy())


reader.close()
output_file = TFile(options.outFile, 'RECREATE')
for hist in hists:
    hist.Write()
output_file.Close()































# ---------------------------------------------------------------------------
# COLLECTION NAME               COLLECTION TYPE          NUMBER OF ELEMENTS
# ===========================================================================
# AllTracks                     Track                          280
# ECalBarrelCollection          SimCalorimeterHit           125823
# ECalEndcapCollection          SimCalorimeterHit             2217
# EcalBarrelCollectionConed     CalorimeterHit               85597
# EcalBarrelCollectionDigi      CalorimeterHit               85673
# EcalBarrelCollectionRec       CalorimeterHit               85673
# EcalBarrelCollectionSel       CalorimeterHit               84247
# EcalBarrelRelationsSimConed   LCRelation                   85597
# EcalBarrelRelationsSimDigi    LCRelation                   85673
# EcalBarrelRelationsSimRec     LCRelation                   85673
# EcalBarrelRelationsSimSel     LCRelation                   84247
# EcalEndcapCollectionConed     CalorimeterHit                 218
# EcalEndcapCollectionDigi      CalorimeterHit                 283
# EcalEndcapCollectionRec       CalorimeterHit                 283
# EcalEndcapCollectionSel       CalorimeterHit                 163
# EcalEndcapRelationsSimConed   LCRelation                     218
# EcalEndcapRelationsSimDigi    LCRelation                     283
# EcalEndcapRelationsSimRec     LCRelation                     283
# EcalEndcapRelationsSimSel     LCRelation                     163
# HCalBarrelCollection          SimCalorimeterHit            89029
# HCalEndcapCollection          SimCalorimeterHit            13021
# HcalBarrelCollectionConed     CalorimeterHit               11064
# HcalBarrelCollectionDigi      CalorimeterHit               11070
# HcalBarrelCollectionRec       CalorimeterHit               11070
# HcalBarrelRelationsSimConed   LCRelation                   11064
# HcalBarrelRelationsSimDigi    LCRelation                   11070
# HcalBarrelRelationsSimRec     LCRelation                   11070
# HcalEndcapCollectionConed     CalorimeterHit                 588
# HcalEndcapCollectionDigi      CalorimeterHit                 589
# HcalEndcapCollectionRec       CalorimeterHit                 589
# HcalEndcapRelationsSimConed   LCRelation                     588
# HcalEndcapRelationsSimDigi    LCRelation                     589
# HcalEndcapRelationsSimRec     LCRelation                     589
# IBTrackerHits                 TrackerHitPlane                214
# IBTrackerHitsConed            TrackerHitPlane                214
# IBTrackerHitsRelations        LCRelation                     214
# IBTrackerHitsRelationsConed   LCRelation                     214
# IETrackerHits                 TrackerHitPlane                  0
# IETrackerHitsConed            TrackerHitPlane                  0
# IETrackerHitsRelations        LCRelation                       0
# IETrackerHitsRelationsConed   LCRelation                       0
# InnerTrackerBarrelCollection  SimTrackerHit                  252
# InnerTrackerBarrelCollectionConedSimTrackerHit                  214
# InnerTrackerEndcapCollection  SimTrackerHit                   49
# InnerTrackerEndcapCollectionConedSimTrackerHit                    0
# JetOut                        ReconstructedParticle            2
# MCParticle                    MCParticle                    3129
# MCParticle_SiTracks_Refitted  LCRelation                      39
# MUON                          CalorimeterHit                   2
# OBTrackerHits                 TrackerHitPlane                190
# OBTrackerHitsConed            TrackerHitPlane                190
# OBTrackerHitsRelations        LCRelation                     190
# OBTrackerHitsRelationsConed   LCRelation                     190
# OETrackerHits                 TrackerHitPlane                  7
# OETrackerHitsConed            TrackerHitPlane                  7
# OETrackerHitsRelations        LCRelation                       7
# OETrackerHitsRelationsConed   LCRelation                       7
# OuterTrackerBarrelCollection  SimTrackerHit                  375
# OuterTrackerBarrelCollectionConedSimTrackerHit                  190
# OuterTrackerEndcapCollection  SimTrackerHit                  311
# OuterTrackerEndcapCollectionConedSimTrackerHit                    7
# PandoraClusters               Cluster                        126
# PandoraPFOs                   ReconstructedParticle          126
# PandoraStartVertices          Vertex                         126
# RelationMuonHit               LCRelation                       2
# SeedTracks                    Track                          280
# SiTracks                      Track                           71
# SiTracks_Refitted             Track                           40
# VBTrackerHits                 TrackerHitPlane                183
# VBTrackerHitsConed            TrackerHitPlane                183
# VBTrackerHitsRelations        LCRelation                     183
# VBTrackerHitsRelationsConed   LCRelation                     183
# VETrackerHits                 TrackerHitPlane                 96
# VETrackerHitsConed            TrackerHitPlane                 96
# VETrackerHitsRelations        LCRelation                      96
# VETrackerHitsRelationsConed   LCRelation                      96
# VertexBarrelCollection        SimTrackerHit                  184
# VertexBarrelCollectionConed   SimTrackerHit                  183
# VertexEndcapCollection        SimTrackerHit                   97
# VertexEndcapCollectionConed   SimTrackerHit                   96
# YokeBarrelCollection          SimCalorimeterHit                2
# YokeEndcapCollection          SimCalorimeterHit                0
# ---------------------------------------------------------------------------

    # etap = []
    # etan = []

    # print("here")
    # pfosCollection = event.getCollection("PandoraPFOs")
    # for ipfo, pfo in enumerate(pfosCollection):
    #     px = pfo.getMomentum()[0]
    #     py = pfo.getMomentum()[1]
    #     pz = pfo.getMomentum()[2]
    #     E = pfo.getEnergy()
    #     tlv = TLorentzVector()
    #     tlv.SetPxPyPzE(px, py, pz, E)
    #     pt = tlv.Pt()
    #     eta = tlv.Eta()
    #     phi = tlv.Phi()
    #     mass = tlv.M()
    #     if eta > 0:
    #         etap.append(pfo)
    #     else:
    #         etan.append(pfo)
    #     # print(f"pt: {pt} , eta: {eta} , phi: {phi} , mass: {mass} , ")

    # tlv_p = TLorentzVector()
    # for pfo in etap:
    #     px = pfo.getMomentum()[0]
    #     py = pfo.getMomentum()[1]
    #     pz = pfo.getMomentum()[2]
    #     E = pfo.getEnergy()
    #     tlv = TLorentzVector()
    #     tlv.SetPxPyPzE(px, py, pz, E)
    #     tlv_p += tlv
    #     print(f"tlv.Pt(): {tlv.Pt()} , tlv.Eta(): {tlv.Eta()} , tlv.Phi(): {tlv.Phi()} , tlv.M(): {tlv.M()} , tlv.E(): {tlv.E()} , ")

    # print(f"tlv_p.Pt(): {tlv_p.Pt()} , tlv_p.Eta(): {tlv_p.Eta()} , tlv_p.Phi(): {tlv_p.Phi()} , tlv_p.M(): {tlv_p.M()} , tlv_p.E(): {tlv_p.E()} , ")

    # tlv_n = TLorentzVector()
    # for pfo in etan:
    #     px = pfo.getMomentum()[0]
    #     py = pfo.getMomentum()[1]
    #     pz = pfo.getMomentum()[2]
    #     E = pfo.getEnergy()
    #     tlv = TLorentzVector()
    #     tlv.SetPxPyPzE(px, py, pz, E)
    #     tlv_n += tlv
    #     print(f"tlv.Pt(): {tlv.Pt()} , tlv.Eta(): {tlv.Eta()} , tlv.Phi(): {tlv.Phi()} , tlv.M(): {tlv.M()} , tlv.E(): {tlv.E()} , ")

    # print(f"tlv_n.Pt(): {tlv_n.Pt()} , tlv_n.Eta(): {tlv_n.Eta()} , tlv_n.Phi(): {tlv_n.Phi()} , tlv_n.M(): {tlv_n.M()} , tlv_n.E(): {tlv_n.E()} , ")

    # tlv_all = tlv_p + tlv_n

    # print(f"tlv_all.Pt(): {tlv_all.Pt()} , tlv_all.Eta(): {tlv_all.Eta()} , tlv_all.Phi(): {tlv_all.Phi()} , tlv_all.M(): {tlv_all.M()} , tlv_all.E(): {tlv_all.E()} , ")


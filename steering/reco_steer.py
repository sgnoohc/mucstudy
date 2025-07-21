from Configurables import ApplicationMgr
from Gaudi.Configuration import *

from Configurables import LcioEvent, EventDataSvc, MarlinProcessorWrapper
from k4MarlinWrapper.parseConstants import *

import os

my_mucoll_basedir = os.getenv("MY_MUCOLL_BASEDIR")
runtag = os.getenv("RUNTAG")

from k4FWCore.parseArgs import parser

parser.add_argument("--enableBIB", action="store_true", default=False, help="Enable BIB overlay")
parser.add_argument("--enableIP", action="store_true", default=False, help="Enable IP overlay")
parser.add_argument("--InFileName", type=str, default="0", help="Input file name for the simulation")
parser.add_argument("--NEvents", type=int, default=-1, help="Input file name for the simulation")
the_args = parser.parse_args()

############################################################################################
# modifying file path
in_path = the_args.InFileName
dirname = os.path.dirname(in_path)
basename = os.path.basename(in_path)

# Check if '_sim' is present
if "_sim" not in basename:
    raise ValueError(f"Input filename must contain '_sim': {basename}")

# Replace _sim with _reco if present
if the_args.enableBIB:
    basename = basename.replace("_reco", "_recoBIB")
else:
    basename = basename.replace("_sim", "_reco")

# Prepend prefix
aida_basename = f"lctuple_actsseededckf_{basename}"

# Join back to full path
aida_path = os.path.join(dirname, aida_basename)
reco_path = os.path.join(dirname, basename)
############################################################################################

algList = []
evtsvc = EventDataSvc()

CONSTANTS = {
}

parseConstants(CONSTANTS)

read = LcioEvent()
read.OutputLevel = INFO
read.Files = [the_args.InFileName]
algList.append(read)

EventNumber = MarlinProcessorWrapper("EventNumber")
EventNumber.OutputLevel = INFO
EventNumber.ProcessorType = "Statusmonitor"
EventNumber.Parameters = {
    "HowOften": ["1"]
}


MyAIDAProcessor = MarlinProcessorWrapper("MyAIDAProcessor")
MyAIDAProcessor.OutputLevel = INFO
MyAIDAProcessor.ProcessorType = "AIDAProcessor"
MyAIDAProcessor.Parameters = {
    "FileName": [aida_path],
    "FileType": ["root"]
}

Output_REC = MarlinProcessorWrapper("Output_REC")
Output_REC.OutputLevel = INFO
Output_REC.ProcessorType = "LCIOOutputProcessor"
if not the_args.enableBIB:
    Output_REC.Parameters = {
        "DropCollectionTypes": [],
        "DropCollectionNames": [],
        "FullSubsetCollections": [],
        "KeepCollectionNames": ["MCParticle_SiTracks_Refitted"],
        "LCIOOutputFile": [reco_path],
        "LCIOWriteMode": ["WRITE_NEW"]
    }
else:
    Output_REC.Parameters = {
        "DropCollectionTypes": [
            # "SimTrackerHit", 
            # "CalorimeterHit", "TrackerHitPlane",
            # "LCRelation"
        ],
        "DropCollectionNames": [
            # "AllTracks", "SeedTracks", "SiTracks",
            # "MCPhysicsParticles", "MCPhysicsParticles_IP",
            # "EcalBarrelRelationsSimDigi", "EcalBarrelRelationsSimRec", "EcalBarrelRelationsSimSel", 
            # "EcalEndcapRelationsSimDigi", "EcalEndcapRelationsSimRec", "EcalEndcapRelationsSimSel",  
            # "HcalBarrelRelationsSimDigi", "HcalBarrelRelationsSimRec", "HcalBarrelRelationsSimSel",  
            # "HcalEndcapRelationsSimDigi", "HcalEndcapRelationsSimRec", "HcalEndcapRelationsSimSel"
        ],
        "FullSubsetCollections": [
            # "EcalBarrelCollectionSel", "EcalEndcapCollectionSel",
            # "HcalBarrelCollectionConed", "HcalEndcapCollectionConed",
            # "IBTrackerHitsConed", "IETrackerHitsConed",
            # "OBTrackerHitsConed", "OETrackerHitsConed", 
            # "VBTrackerHitsConed", "VETrackerHitsConed", 
            # "SiTracks_Refitted"
        ],
        "KeepCollectionNames": [
            # "EcalBarrelCollectionSel", "EcalEndcapCollectionSel", 
            # "HcalBarrelCollectionConed", "HcalEndcapCollectionConed",
            # "IBTrackerHitsConed", "IETrackerHitsConed",
            # "OBTrackerHitsConed", "OETrackerHitsConed", 
            # "VBTrackerHitsConed", "VETrackerHitsConed", 
            # "SiTracks_Refitted",
            "MCParticle_SiTracks_Refitted"
        ],
        "LCIOOutputFile": [reco_path],
        "LCIOWriteMode": ["WRITE_NEW"]
    }

InitDD4hep = MarlinProcessorWrapper("InitDD4hep")
InitDD4hep.OutputLevel = INFO
InitDD4hep.ProcessorType = "InitializeDD4hep"
InitDD4hep.Parameters = {
    "DD4hepXMLFile": [my_mucoll_basedir+"/../detector-simulation/geometries/MAIA_v0/MAIA_v0.xml"],
    #"DD4hepXMLFile": ["/code/detector-simulation/geometries/MuColl_10TeV_v0A/MuColl_10TeV_v0A.xml"],
    "EncodingStringParameterName": ["GlobalTrackerReadoutID"]
}

VXDBarrelDigitiser = MarlinProcessorWrapper("VXDBarrelDigitiser")
VXDBarrelDigitiser.OutputLevel = INFO
VXDBarrelDigitiser.ProcessorType = "DDPlanarDigiProcessor"
VXDBarrelDigitiser.Parameters = {
    "CorrectTimesForPropagation": ["true"],
    "IsStrip": ["false"],
    "ResolutionT": ["0.03"],
    "ResolutionU": ["0.005"],
    "ResolutionV": ["0.005"],
    "SimTrackHitCollectionName": ["VertexBarrelCollection"],
    "SimTrkHitRelCollection": ["VBTrackerHitsRelations"],
    "SubDetectorName": ["Vertex"],
    "TimeWindowMax": ["0.15"],
    "TimeWindowMin": ["-0.09"],
    "TrackerHitCollectionName": ["VBTrackerHits"],
    "UseTimeWindow": ["true"]
}

VXDEndcapDigitiser = MarlinProcessorWrapper("VXDEndcapDigitiser")
VXDEndcapDigitiser.OutputLevel = INFO
VXDEndcapDigitiser.ProcessorType = "DDPlanarDigiProcessor"
VXDEndcapDigitiser.Parameters = {
    "CorrectTimesForPropagation": ["true"],
    "IsStrip": ["false"],
    "ResolutionT": ["0.03"],
    "ResolutionU": ["0.005"],
    "ResolutionV": ["0.005"],
    "SimTrackHitCollectionName": ["VertexEndcapCollection"],
    "SimTrkHitRelCollection": ["VETrackerHitsRelations"],
    "SubDetectorName": ["Vertex"],
    "TimeWindowMax": ["0.15"],
    "TimeWindowMin": ["-0.09"],
    "TrackerHitCollectionName": ["VETrackerHits"],
    "UseTimeWindow": ["true"]
}

InnerPlanarDigiProcessor = MarlinProcessorWrapper("InnerPlanarDigiProcessor")
InnerPlanarDigiProcessor.OutputLevel = INFO
InnerPlanarDigiProcessor.ProcessorType = "DDPlanarDigiProcessor"
InnerPlanarDigiProcessor.Parameters = {
    "CorrectTimesForPropagation": ["true"],
    "IsStrip": ["false"],
    "ResolutionT": ["0.06"],
    "ResolutionU": ["0.007"],
    "ResolutionV": ["0.090"],
    "SimTrackHitCollectionName": ["InnerTrackerBarrelCollection"],
    "SimTrkHitRelCollection": ["IBTrackerHitsRelations"],
    "SubDetectorName": ["InnerTrackers"],
    "TimeWindowMax": ["0.3"],
    "TimeWindowMin": ["-0.18"],
    "TrackerHitCollectionName": ["IBTrackerHits"],
    "UseTimeWindow": ["true"]
}

InnerEndcapPlanarDigiProcessor = MarlinProcessorWrapper(
    "InnerEndcapPlanarDigiProcessor")
InnerEndcapPlanarDigiProcessor.OutputLevel = INFO
InnerEndcapPlanarDigiProcessor.ProcessorType = "DDPlanarDigiProcessor"
InnerEndcapPlanarDigiProcessor.Parameters = {
    "CorrectTimesForPropagation": ["true"],
    "IsStrip": ["false"],
    "ResolutionT": ["0.06"],
    "ResolutionU": ["0.007"],
    "ResolutionV": ["0.090"],
    "SimTrackHitCollectionName": ["InnerTrackerEndcapCollection"],
    "SimTrkHitRelCollection": ["IETrackerHitsRelations"],
    "SubDetectorName": ["InnerTrackers"],
    "TimeWindowMax": ["0.3"],
    "TimeWindowMin": ["-0.18"],
    "TrackerHitCollectionName": ["IETrackerHits"],
    "UseTimeWindow": ["true"]
}

OuterPlanarDigiProcessor = MarlinProcessorWrapper("OuterPlanarDigiProcessor")
OuterPlanarDigiProcessor.OutputLevel = INFO
OuterPlanarDigiProcessor.ProcessorType = "DDPlanarDigiProcessor"
OuterPlanarDigiProcessor.Parameters = {
    "CorrectTimesForPropagation": ["true"],
    "IsStrip": ["false"],
    "ResolutionT": ["0.06"],
    "ResolutionU": ["0.007"],
    "ResolutionV": ["0.090"],
    "SimTrackHitCollectionName": ["OuterTrackerBarrelCollection"],
    "SimTrkHitRelCollection": ["OBTrackerHitsRelations"],
    "SubDetectorName": ["OuterTrackers"],
    "TimeWindowMax": ["0.3"],
    "TimeWindowMin": ["-0.18"],
    "TrackerHitCollectionName": ["OBTrackerHits"],
    "UseTimeWindow": ["true"]
}

OuterEndcapPlanarDigiProcessor = MarlinProcessorWrapper(
    "OuterEndcapPlanarDigiProcessor")
OuterEndcapPlanarDigiProcessor.OutputLevel = INFO
OuterEndcapPlanarDigiProcessor.ProcessorType = "DDPlanarDigiProcessor"
OuterEndcapPlanarDigiProcessor.Parameters = {
    "CorrectTimesForPropagation": ["true"],
    "IsStrip": ["false"],
    "ResolutionT": ["0.06"],
    "ResolutionU": ["0.007"],
    "ResolutionV": ["0.090"],
    "SimTrackHitCollectionName": ["OuterTrackerEndcapCollection"],
    "SimTrkHitRelCollection": ["OETrackerHitsRelations"],
    "SubDetectorName": ["OuterTrackers"],
    "TimeWindowMax": ["0.3"],
    "TimeWindowMin": ["-0.18"],
    "TrackerHitCollectionName": ["OETrackerHits"],
    "UseTimeWindow": ["true"]
}

VXDBarrelConer = MarlinProcessorWrapper("VXDBarrelConer")
VXDBarrelConer.OutputLevel = INFO
VXDBarrelConer.ProcessorType = "FilterConeHits"
VXDBarrelConer.Parameters = {
    "MCParticleCollection": ["MCParticle"],
    "TrackerHitInputCollections": ["VBTrackerHits"],
    "TrackerSimHitInputCollections": ["VertexBarrelCollection"],
    "TrackerHitInputRelations": ["VBTrackerHitsRelations"],
    "TrackerHitOutputCollections": ["VBTrackerHitsConed"],
    "TrackerSimHitOutputCollections": ["VertexBarrelCollectionConed"],
    "TrackerHitOutputRelations": ["VBTrackerHitsRelationsConed"],
    "DeltaRCut": ["0.8"],
    "FillHistograms": ["false"]
}

VXDEndcapConer = MarlinProcessorWrapper("VXDEndcapConer")
VXDEndcapConer.OutputLevel = INFO
VXDEndcapConer.ProcessorType = "FilterConeHits"
VXDEndcapConer.Parameters = {
    "MCParticleCollection": ["MCParticle"],
    "TrackerHitInputCollections": ["VETrackerHits"],
    "TrackerSimHitInputCollections": ["VertexEndcapCollection"],
    "TrackerHitInputRelations": ["VETrackerHitsRelations"],
    "TrackerHitOutputCollections": ["VETrackerHitsConed"],
    "TrackerSimHitOutputCollections": ["VertexEndcapCollectionConed"],
    "TrackerHitOutputRelations": ["VETrackerHitsRelationsConed"],
    "DeltaRCut": ["0.8"],
    "FillHistograms": ["false"]
}

InnerPlanarConer = MarlinProcessorWrapper("InnerPlanarConer")
InnerPlanarConer.OutputLevel = INFO
InnerPlanarConer.ProcessorType = "FilterConeHits"
InnerPlanarConer.Parameters = {
    "MCParticleCollection": ["MCParticle"],
    "TrackerHitInputCollections": ["IBTrackerHits"],
    "TrackerSimHitInputCollections": ["InnerTrackerBarrelCollection"],
    "TrackerHitInputRelations": ["IBTrackerHitsRelations"],
    "TrackerHitOutputCollections": ["IBTrackerHitsConed"],
    "TrackerSimHitOutputCollections": ["InnerTrackerBarrelCollectionConed"],
    "TrackerHitOutputRelations": ["IBTrackerHitsRelationsConed"],
    "DeltaRCut": ["0.8"],
    "FillHistograms": ["false"]
}

InnerEndcapConer = MarlinProcessorWrapper("InnerEndcapConer")
InnerEndcapConer.OutputLevel = INFO
InnerEndcapConer.ProcessorType = "FilterConeHits"
InnerEndcapConer.Parameters = {
    "MCParticleCollection": ["MCParticle"],
    "TrackerHitInputCollections": ["IETrackerHits"],
    "TrackerSimHitInputCollections": ["InnerTrackerEndcapCollection"],
    "TrackerHitInputRelations": ["IETrackerHitsRelations"],
    "TrackerHitOutputCollections": ["IETrackerHitsConed"],
    "TrackerSimHitOutputCollections": ["InnerTrackerEndcapCollectionConed"],
    "TrackerHitOutputRelations": ["IETrackerHitsRelationsConed"],
    "DeltaRCut": ["0.8"],
    "FillHistograms": ["false"]
}

OuterPlanarConer = MarlinProcessorWrapper("OuterPlanarConer")
OuterPlanarConer.OutputLevel = INFO
OuterPlanarConer.ProcessorType = "FilterConeHits"
OuterPlanarConer.Parameters = {
    "MCParticleCollection": ["MCParticle"],
    "TrackerHitInputCollections": ["OBTrackerHits"],
    "TrackerSimHitInputCollections": ["OuterTrackerBarrelCollection"],
    "TrackerHitInputRelations": ["OBTrackerHitsRelations"],
    "TrackerHitOutputCollections": ["OBTrackerHitsConed"],
    "TrackerSimHitOutputCollections": ["OuterTrackerBarrelCollectionConed"],
    "TrackerHitOutputRelations": ["OBTrackerHitsRelationsConed"],
    "DeltaRCut": ["0.8"],
    "FillHistograms": ["false"]
}

OuterEndcapConer = MarlinProcessorWrapper("OuterEndcapConer")
OuterEndcapConer.OutputLevel = INFO
OuterEndcapConer.ProcessorType = "FilterConeHits"
OuterEndcapConer.Parameters = {
    "MCParticleCollection": ["MCParticle"],
    "TrackerHitInputCollections": ["OETrackerHits"],
    "TrackerSimHitInputCollections": ["OuterTrackerEndcapCollection"],
    "TrackerHitInputRelations": ["OETrackerHitsRelations"],
    "TrackerHitOutputCollections": ["OETrackerHitsConed"],
    "TrackerSimHitOutputCollections": ["OuterTrackerEndcapCollectionConed"],
    "TrackerHitOutputRelations": ["OETrackerHitsRelationsConed"],
    "DeltaRCut": ["0.8"],
    "FillHistograms": ["false"]
}

CKFTracking = MarlinProcessorWrapper("CKFTracking")
CKFTracking.OutputLevel = INFO
CKFTracking.ProcessorType = "ACTSSeededCKFTrackingProc"
CKFTracking.Parameters = {
    "CKF_Chi2CutOff": ["10"],
    "CKF_NumMeasurementsCutOff": ["1"],
    "CaloFace_Radius": ["1857"],
    "CaloFace_Z": ["2307"],
    "MatFile": [my_mucoll_basedir+"/../ACTSTracking/data/MAIA_v0_material.json"],
    "PropagateBackward": ["False"],
    "DetectorSchema": ["MAIA_v0"],
    "RunCKF": ["True"],
    "SeedFinding_CollisionRegion": ["6"],
    #"SeedFinding_DeltaRMax": ["60"],
    #"SeedFinding_DeltaRMin": ["2"],
    #"SeedFinding_DeltaRMaxBottom": ["50"],
    #"SeedFinding_DeltaRMaxTop": ["50"],
    #"SeedFinding_DeltaRMinBottom": ["5"],
    #"SeedFinding_DeltaRMinTop": ["2"],
    "SeedFinding_ImpactMax": ["3"],
    "SeedFinding_MinPt": ["500"],
    "SeedFinding_RMax": ["150"],
    "SeedFinding_ZMax": ["600"],
    "SeedFinding_RadLengthPerSeed": ["0.1"],
    #"SeedFinding_zBottomBinLen": ["1"],
    #"SeedFinding_zTopBinLen": ["2"],
    #"SeedFinding_phiBottomBinLen": ["1"],
    #"SeedFinding_phiTopBinLen": ["2"],
    "SeedFinding_SigmaScattering": ["50"],
    #"SeedFinding_zBottomBinLen": ["0"],
    #"SeedFinding_zTopBinLen": ["25"],
    #"SeedFinding_phiBottomBinLen": ["25"],
    #"SeedFinding_phiTopBinLen": ["50"],
    "SeedingLayers": ["13", "2", "13", "6", "13", "10", "13", "14",
                      "14", "2", "14", "6", "14", "8", "14", "10", 
                      "15", "2", "15", "6", "15", "10", "15", "14",
                      "8", "2",
                      "17", "2",
                      "18", "2"],
    "TGeoFile": [my_mucoll_basedir+"/../ACTSTracking/data/MAIA_v0.root"],
    "TGeoDescFile": [os.environ['ACTSTRACKING_DATA']+"/MAIA_v0.json"],
    "TrackCollectionName": ["AllTracks"],
    "TrackerHitCollectionNames": ["VBTrackerHitsConed", "IBTrackerHitsConed", "OBTrackerHitsConed", "VETrackerHitsConed", "IETrackerHitsConed", "OETrackerHitsConed"]
}

TrackDeduper = MarlinProcessorWrapper("TrackDeduper")
TrackDeduper.OutputLevel = INFO
TrackDeduper.ProcessorType = "ACTSDuplicateRemoval"
TrackDeduper.Parameters = {
    "InputTrackCollectionName": ["AllTracks"],
    "OutputTrackCollectionName": ["SiTracks"]
}

Refit = MarlinProcessorWrapper("Refit")
Refit.OutputLevel = WARNING
Refit.ProcessorType = "RefitFinal"
Refit.Parameters = {
    "DoCutsOnRedChi2Nhits": ["true"],
    "EnergyLossOn": ["true"],
    "InputRelationCollectionName": ["SiTrackRelations"],
    "InputTrackCollectionName": ["SiTracks"],
    "Max_Chi2_Incr": ["1.79769e+30"],
    "MinClustersOnTrackAfterFit": ["3"],
    "MultipleScatteringOn": ["true"],
    "NHitsCuts": ["1,2", "1", "3,4", "1", "5,6", "0"],
    "OutputRelationCollectionName": ["SiTracks_Refitted_Relation"],
    "OutputTrackCollectionName": ["SiTracks_Refitted"],
    "ReducedChi2Cut": ["10."],
    "ReferencePoint": ["-1"],
    "SmoothOn": ["false"],
    "extrapolateForward": ["true"]
}

MyTrackSelector = MarlinProcessorWrapper("MyTrackSelector")
MyTrackSelector.OutputLevel = INFO
MyTrackSelector.ProcessorType = "FilterTracks"
MyTrackSelector.Parameters = {
    "BarrelOnly": ["false"],
    "HasCaloState": ["true"],
    "NHitsTotal": ["5"],
    "NHitsVertex": ["2"],
    "NHitsInner": ["1"],
    "NHitsOuter": ["1"],
    "MinPt": ["0.5"],
    "Chi2Spatial": ["0"],
    "MaxHoles": ["5"],
    "InputTrackCollectionName": ["SiTracks"],
    "OutputTrackCollectionName": ["SelectedTracks"]
}

MyTrackTruth = MarlinProcessorWrapper("MyTrackTruth")
MyTrackTruth.OutputLevel = INFO
MyTrackTruth.ProcessorType = "TrackTruthProc"
MyTrackTruth.Parameters = {
    "MCParticleCollection": ["MCParticle"],
    "Particle2TrackRelationName": ["MCParticle_SiTracks_Refitted"],
    "TrackCollection": ["SiTracks_Refitted"],
    "TrackerHit2SimTrackerHitRelationName": ["VBTrackerHitsRelationsConed", "IBTrackerHitsRelationConed", "OBTrackerHitsRelationsConed", "VETrackerHitsRelationsConed", "IETrackerHitsRelationsConed", "OETrackerHitsRelationsConed"]
}

MyEcalBarrelDigi = MarlinProcessorWrapper("MyEcalBarrelDigi")
MyEcalBarrelDigi.OutputLevel = INFO
MyEcalBarrelDigi.ProcessorType = "RealisticCaloDigiSilicon"
MyEcalBarrelDigi.Parameters = {
    "CellIDLayerString": ["layer"],
    "calibration_mip": ["0.0001575"],
    "inputHitCollections": ["ECalBarrelCollection"],
    "outputHitCollections": ["EcalBarrelCollectionDigi"],
    "outputRelationCollections": ["EcalBarrelRelationsSimDigi"],
    #"threshold": ["0.002"],
    "threshold": ["5e-05"],
    "thresholdUnit": ["GeV"],
    "timingCorrectForPropagation": ["1"],
    "timingCut": ["1"],
    "timingResolution": ["0"],
    "timingWindowMax": ["10"],
    "timingWindowMin": ["-0.5"],
    "elec_range_mip": ["15000"]
}

MyEcalBarrelReco = MarlinProcessorWrapper("MyEcalBarrelReco")
MyEcalBarrelReco.OutputLevel = INFO
MyEcalBarrelReco.ProcessorType = "RealisticCaloRecoSilicon"
MyEcalBarrelReco.Parameters = {
    "CellIDLayerString": ["layer"],
    #    "calibration_factorsMipGev": ["0.00641222630095"],
    "calibration_factorsMipGev": ["0.0066150"], #used for v3
    #"calibration_factorsMipGev": ["0.00826875"],
    "calibration_layergroups": ["50"],
    "inputHitCollections": ["EcalBarrelCollectionDigi"],
    "inputRelationCollections": ["EcalBarrelRelationsSimDigi"],
    "outputHitCollections": ["EcalBarrelCollectionRec"],
    "outputRelationCollections": ["EcalBarrelRelationsSimRec"]
}

MyEcalEndcapDigi = MarlinProcessorWrapper("MyEcalEndcapDigi")
MyEcalEndcapDigi.OutputLevel = INFO
MyEcalEndcapDigi.ProcessorType = "RealisticCaloDigiSilicon"
MyEcalEndcapDigi.Parameters = {
    "CellIDLayerString": ["layer"],
    "calibration_mip": ["0.0001575"],
    "inputHitCollections": ["ECalEndcapCollection"],
    "outputHitCollections": ["EcalEndcapCollectionDigi"],
    "outputRelationCollections": ["EcalEndcapRelationsSimDigi"],
    #"threshold": ["0.002"],
    "threshold": ["5e-05"],
    "thresholdUnit": ["GeV"],
    "timingCorrectForPropagation": ["1"],
    "timingCut": ["1"],
    "timingResolution": ["0"],
    "timingWindowMax": ["10"],
    "timingWindowMin": ["-0.5"],
    "elec_range_mip": ["15000"]
}

MyEcalEndcapReco = MarlinProcessorWrapper("MyEcalEndcapReco")
MyEcalEndcapReco.OutputLevel = INFO
MyEcalEndcapReco.ProcessorType = "RealisticCaloRecoSilicon"
MyEcalEndcapReco.Parameters = {
    "CellIDLayerString": ["layer"],
    #    "calibration_factorsMipGev": ["0.00641222630095"],
    "calibration_factorsMipGev": ["0.0066150"], #used for v3
    #"calibration_factorsMipGev": ["0.00826875"],
    "calibration_layergroups": ["50"],
    "inputHitCollections": ["EcalEndcapCollectionDigi"],
    "inputRelationCollections": ["EcalEndcapRelationsSimDigi"],
    "outputHitCollections": ["EcalEndcapCollectionRec"],
    "outputRelationCollections": ["EcalEndcapRelationsSimRec"]
}

MyHcalBarrelDigi = MarlinProcessorWrapper("MyHcalBarrelDigi")
MyHcalBarrelDigi.OutputLevel = INFO
MyHcalBarrelDigi.ProcessorType = "RealisticCaloDigiScinPpd"
MyHcalBarrelDigi.Parameters = {
    "CellIDLayerString": ["layer"],
    "calibration_mip": ["0.0004725"],
    "inputHitCollections": ["HCalBarrelCollection"],
    "outputHitCollections": ["HcalBarrelCollectionDigi"],
    "outputRelationCollections": ["HcalBarrelRelationsSimDigi"],
    "ppd_mipPe": ["15"],
    "ppd_npix": ["2000"],
    "ppd_npix_uncert": ["0"],
    "ppd_pix_spread": ["0"],
    "threshold": ["0.5"],
    "thresholdUnit": ["MIP"],
    "timingCorrectForPropagation": ["1"],
    "timingCut": ["1"],
    "timingResolution": ["0"],
    "timingWindowMax": ["10"],
    "timingWindowMin": ["-0.5"]
}

MyHcalBarrelReco = MarlinProcessorWrapper("MyHcalBarrelReco")
MyHcalBarrelReco.OutputLevel = INFO
MyHcalBarrelReco.ProcessorType = "RealisticCaloRecoScinPpd"
MyHcalBarrelReco.Parameters = {
    "CellIDLayerString": ["layer"],
    #    "calibration_factorsMipGev": ["0.0287783798145"],
    "calibration_factorsMipGev": ["0.024625"],
    "calibration_layergroups": ["100"],
    "inputHitCollections": ["HcalBarrelCollectionDigi"],
    "inputRelationCollections": ["HcalBarrelRelationsSimDigi"],
    "outputHitCollections": ["HcalBarrelCollectionRec"],
    "outputRelationCollections": ["HcalBarrelRelationsSimRec"],
    "ppd_mipPe": ["15"],
    "ppd_npix": ["2000"]
}

MyHcalEndcapDigi = MarlinProcessorWrapper("MyHcalEndcapDigi")
MyHcalEndcapDigi.OutputLevel = INFO
MyHcalEndcapDigi.ProcessorType = "RealisticCaloDigiScinPpd"
MyHcalEndcapDigi.Parameters = {
    "CellIDLayerString": ["layer"],
    "calibration_mip": ["0.0004725"],
    "inputHitCollections": ["HCalEndcapCollection"],
    "outputHitCollections": ["HcalEndcapCollectionDigi"],
    "outputRelationCollections": ["HcalEndcapRelationsSimDigi"],
    "ppd_mipPe": ["15"],
    "ppd_npix": ["2000"],
    "ppd_npix_uncert": ["0"],
    "ppd_pix_spread": ["0"],
    "threshold": ["0.5"],
    "thresholdUnit": ["MIP"],
    "timingCorrectForPropagation": ["1"],
    "timingCut": ["1"],
    "timingResolution": ["0"],
    "timingWindowMax": ["10"],
    "timingWindowMin": ["-0.5"]
}

MyHcalEndcapReco = MarlinProcessorWrapper("MyHcalEndcapReco")
MyHcalEndcapReco.OutputLevel = INFO
MyHcalEndcapReco.ProcessorType = "RealisticCaloRecoScinPpd"
MyHcalEndcapReco.Parameters = {
    "CellIDLayerString": ["layer"],
    #   "calibration_factorsMipGev": ["0.0285819096797"],
    "calibration_factorsMipGev": ["0.024625"],
    "calibration_layergroups": ["100"],
    "inputHitCollections": ["HcalEndcapCollectionDigi"],
    "inputRelationCollections": ["HcalEndcapRelationsSimDigi"],
    "outputHitCollections": ["HcalEndcapCollectionRec"],
    "outputRelationCollections": ["HcalEndcapRelationsSimRec"],
    "ppd_mipPe": ["15"],
    "ppd_npix": ["2000"]
}

MyEcalBarrelConer = MarlinProcessorWrapper("MyEcalBarrelConer")
MyEcalBarrelConer.OutputLevel = INFO
MyEcalBarrelConer.ProcessorType = "CaloConer"
MyEcalBarrelConer.Parameters = {
    "MCParticleCollectionName": ["MCParticle"],
    "CaloHitCollectionName": ["EcalBarrelCollectionRec"],
    "CaloRelationCollectionName": ["EcalBarrelRelationsSimRec"],
    "GoodHitCollection": ["EcalBarrelCollectionConed"],
    "GoodRelationCollection": ["EcalBarrelRelationsSimConed"],
    "ConeWidth": ["0.4"]
}

MyEcalEndcapConer = MarlinProcessorWrapper("MyEcalEndcapConer")
MyEcalEndcapConer.OutputLevel = INFO
MyEcalEndcapConer.ProcessorType = "CaloConer"
MyEcalEndcapConer.Parameters = {
    "MCParticleCollectionName": ["MCParticle"],
    "CaloHitCollectionName": ["EcalEndcapCollectionRec"],
    "CaloRelationCollectionName": ["EcalEndcapRelationsSimRec"],
    "GoodHitCollection": ["EcalEndcapCollectionConed"],
    "GoodRelationCollection": ["EcalEndcapRelationsSimConed"],
    "ConeWidth": ["0.4"]
}

MyHcalBarrelConer = MarlinProcessorWrapper("MyHcalBarrelConer")
MyHcalBarrelConer.OutputLevel = INFO
MyHcalBarrelConer.ProcessorType = "CaloConer"
MyHcalBarrelConer.Parameters = {
    "MCParticleCollectionName": ["MCParticle"],
    "CaloHitCollectionName": ["HcalBarrelCollectionRec"],
    "CaloRelationCollectionName": ["HcalBarrelRelationsSimRec"],
    "GoodHitCollection": ["HcalBarrelCollectionConed"],
    "GoodRelationCollection": ["HcalBarrelRelationsSimConed"],
    "ConeWidth": ["0.4"]
}

MyHcalEndcapConer = MarlinProcessorWrapper("MyHcalEndcapConer")
MyHcalEndcapConer.OutputLevel = INFO
MyHcalEndcapConer.ProcessorType = "CaloConer"
MyHcalEndcapConer.Parameters = {
    "MCParticleCollectionName": ["MCParticle"],
    "CaloHitCollectionName": ["HcalEndcapCollectionRec"],
    "CaloRelationCollectionName": ["HcalEndcapRelationsSimRec"],
    "GoodHitCollection": ["HcalEndcapCollectionConed"],
    "GoodRelationCollection": ["HcalEndcapRelationsSimConed"],
    "ConeWidth": ["0.4"]
}


MyEcalBarrelSelector = MarlinProcessorWrapper("MyEcalBarrelSelector")
MyEcalBarrelSelector.OutputLevel = INFO
MyEcalBarrelSelector.ProcessorType = "CaloHitSelector"
MyEcalBarrelSelector.Parameters = {
    "CaloHitCollectionName": ["EcalBarrelCollectionConed"],
    "CaloRelationCollectionName": ["EcalBarrelRelationsSimConed"],
    "GoodHitCollection": ["EcalBarrelCollectionSel"],
    "GoodRelationCollection": ["EcalBarrelRelationsSimSel"],
    "ThresholdsFilePath": [my_mucoll_basedir+"/../MyBIBUtils/data/ECAL_Thresholds_10TeV.root"],
    "Nsigma": ["0"],
    "DoBIBsubtraction": ["false"]
}

MyEcalEndcapSelector = MarlinProcessorWrapper("MyEcalEndcapSelector")
MyEcalEndcapSelector.OutputLevel = INFO
MyEcalEndcapSelector.ProcessorType = "CaloHitSelector"
MyEcalEndcapSelector.Parameters = {
    "CaloHitCollectionName": ["EcalEndcapCollectionConed"],
    "CaloRelationCollectionName": ["EcalEndcapRelationsSimConed"],
    "GoodHitCollection": ["EcalEndcapCollectionSel"],
    "GoodRelationCollection": ["EcalEndcapRelationsSimSel"],
    "ThresholdsFilePath": [my_mucoll_basedir+"/../MyBIBUtils/data/ECAL_Thresholds_10TeV.root"],
    "Nsigma": ["0"],
    "DoBIBsubtraction": ["false"]
}

DDMarlinPandora = MarlinProcessorWrapper("DDMarlinPandora")
DDMarlinPandora.OutputLevel = INFO
DDMarlinPandora.ProcessorType = "DDPandoraPFANewProcessor"
DDMarlinPandora.Parameters = {
    "ClusterCollectionName": ["PandoraClusters"],
    "CreateGaps": ["false"],
    "CurvatureToMomentumFactor": ["0.00015"],
    "D0TrackCut": ["200"],
    "D0UnmatchedVertexTrackCut": ["5"],
    "DigitalMuonHits": ["0"],
    "ECalBarrelNormalVector": ["0", "0", "1"],
    "ECalCaloHitCollections": ["EcalBarrelCollectionSel", "EcalEndcapCollectionSel"],
    "ECalMipThreshold": ["0.5"],
    "ECalScMipThreshold": ["0"],
    "ECalScToEMGeVCalibration": ["1"],
    "ECalScToHadGeVCalibrationBarrel": ["1"],
    "ECalScToHadGeVCalibrationEndCap": ["1"],
    "ECalScToMipCalibration": ["1"],
    "ECalSiMipThreshold": ["0"],
    "ECalSiToEMGeVCalibration": ["1"],
    "ECalSiToHadGeVCalibrationBarrel": ["1"],
    "ECalSiToHadGeVCalibrationEndCap": ["1"],
    "ECalSiToMipCalibration": ["1"],
    "ECalToEMGeVCalibration": ["1.02373335516"],
    "ECalToHadGeVCalibrationBarrel": ["1.24223718397"],
    "ECalToHadGeVCalibrationEndCap": ["1.24223718397"],
    "ECalToMipCalibration": ["181.818"],
    "EMConstantTerm": ["0.01"], 
    "EMStochasticTerm": ["0.17"],
    "FinalEnergyDensityBin": ["110."],
    "HCalBarrelNormalVector": ["0", "0", "1"],
    "HCalCaloHitCollections": ["HcalBarrelCollectionConed", "HcalEndcapCollectionConed"],
    "HCalMipThreshold": ["0.3"],
    "HCalToEMGeVCalibration": ["1.02373335516"],
    "HCalToHadGeVCalibration": ["1.01799349172"],
    "HCalToMipCalibration": ["40.8163"],
    "HadConstantTerm": ["0.03"],
    "HadStochasticTerm": ["0.6"],
    "InputEnergyCorrectionPoints": [],
    "OutputEnergyCorrectionPoints": [],
    #"InputEnergyCorrectionPoints": ["1.166", "1.772", "1.468", "1.844", "2.384", "2.737", "3.085", "3.886", "3.97", "4.999", "5.64", "6.328", "6.909", "7.352", "7.893", "9.211", "8.783", "10.494", "9.644", "10.263", "10.754", "10.283", "12.955", "14.203", "15.32", "15.285", "17.759", "24.1", "34.602", "45.971", "66.333", "85.715", "100.868", "121.045", "141.793", "160.122", "190.229", "236.062", "272.852", "324.34", "367.599", "456.216", "508.0", "591.241", "677.222", "864.469", "1110.301", "1378.305", "1763.056", "2135.02", "2372.24"],
    #"OutputEnergyCorrectionPoints": ["5.0", "11.0", "13.0", "15.0", "17.0", "19.0", "22.5", "27.5", "32.5", "37.5", "42.5", "47.5", "52.5", "57.5", "62.5", "67.5", "72.5", "77.5", "82.5", "87.5", "92.5", "97.5", "105.0", "115.0", "125.0", "135.0", "145.0", "175.0", "225.0", "275.0", "325.0", "375.0", "425.0", "475.0", "525.0", "575.0", "650.0", "750.0", "850.0", "950.0", "1050.0", "1150.0", "1250.0", "1350.0", "1450.0", "1750.0", "2250.0", "2750.0", "3500.0", "4500.0", "5000."],
    "KinkVertexCollections": ["KinkVertices"],
    "LayersFromEdgeMaxRearDistance": ["250"],
    "MCParticleCollections": ["MCParticle"],
    "MaxBarrelTrackerInnerRDistance": ["200"],
    "MaxClusterEnergyToApplySoftComp": ["2000."],
    "MaxHCalHitHadronicEnergy": ["1000000"],
    "MaxTrackHits": ["5000"],
    "MaxTrackSigmaPOverP": ["0.15"],
    "MinBarrelTrackerHitFractionOfExpected": ["0"],
    "MinCleanCorrectedHitEnergy": ["0.1"],
    "MinCleanHitEnergy": ["0.5"],
    "MinCleanHitEnergyFraction": ["0.01"],
    "MinFtdHitsForBarrelTrackerHitFraction": ["0"],
    "MinFtdTrackHits": ["0"],
    "MinMomentumForTrackHitChecks": ["0"],
    "MinTpcHitFractionOfExpected": ["0"],
    "MinTrackECalDistanceFromIp": ["0"],
    "MinTrackHits": ["0"],
#    "MuonBarrelBField": ["5.0"],
    "MuonBarrelBField": ["0.0001"],
    "MuonCaloHitCollections": ["MUON"],
#    "MuonEndCapBField": ["5.0"],
    "MuonEndCapBField": ["0.0001"],
    "MuonHitEnergy": ["0.5"],
    "MuonToMipCalibration": ["19607.8"],
    "NEventsToSkip": ["0"],
    "NOuterSamplingLayers": ["3"],
    "PFOCollectionName": ["PandoraPFOs"],
    "PandoraSettingsXmlFile": [my_mucoll_basedir+"/../SteeringMacros/PandoraSettings/PandoraSettingsDefault.xml"],
    "ProngVertexCollections": ["ProngVertices"],
    "ReachesECalBarrelTrackerOuterDistance": ["-100"],
    "ReachesECalBarrelTrackerZMaxDistance": ["-50"],
    "ReachesECalFtdZMaxDistance": ["1"],
    "ReachesECalMinFtdLayer": ["0"],
    "ReachesECalNBarrelTrackerHits": ["0"],
    "ReachesECalNFtdHits": ["0"],
    "RelCaloHitCollections": ["EcalBarrelRelationsSimSel", "EcalEndcapRelationsSimSel", "HcalBarrelRelationsSimConed", "HcalEndcapRelationsSimConed", "RelationMuonHit"],
    "RelTrackCollections": ["SiTracks_Refitted_Relation"],
    "ShouldFormTrackRelationships": ["1"],
    "SoftwareCompensationEnergyDensityBins": ["0", "2.", "5.", "7.5", "9.5", "13.", "16.", "20.", "23.5", "28.", "33.", "40.", "50.", "75.", "100."],
    "SoftwareCompensationWeights": ["1.61741", "-0.00444385", "2.29683e-05", "-0.0731236", "-0.00157099", "-7.09546e-07", "0.868443", "1.0561", "-0.0238574"],
    # ECAL corrections w/ BIB w/ cell selection
    #"ECALInputEnergyCorrectionPoints": ["0.1", "11.894", "12.971", "13.166", "14.926", "15.26", "17.133", "20.106", "23.174", "25.777", 
    #                                    "29.132", "32.219", "34.876", "36.577", "39.751", "42.48", "46.22", "49.708", "53.274", "56.89", 
    #                                    "59.071", "62.913", "67.952", "75.322", "82.061", "89.277", "96.304", "116.911", "153.542", "189.536", 
    #                                    "227.668", "267.495", "308.292", "348.397", "386.6", "428.067", "488.265", "571.82", "655.051", "741.878", 
    #                                    "821.348", "913.81", "1000.185", "1089.007", "1169.009", "1415.481", "1859.895", "2315.64", "2941.297", "3851.56", 
    #                                    "4279.5"],
    #"ECALOutputEnergyCorrectionPoints": ["0.1", "11.0", "13.0", "15.0", "17.0", "19.0", "22.5", "27.5", "32.5", "37.5", 
    #                                     "42.5", "47.5", "52.5", "57.5", "62.5", "67.5", "72.5", "77.5", "82.5", "87.5", 
    #                                     "92.5", "97.5", "105.0", "115.0", "125.0", "135.0", "145.0", "175.0", "225.0", "275.0", 
    #                                     "325.0", "375.0", "425.0", "475.0", "525.0", "575.0", "650.0", "750.0", "850.0", "950.0", 
    #                                     "1050.0", "1150.0", "1250.0", "1350.0", "1450.0", "1750.0", "2250.0", "2750.0", "3500.0", "4500.0", 
    #                                     "5000."],
    # ECAL corrections w/o BIB w/cell selection
    #"ECALInputEnergyCorrectionPoints": ["0.1", 
    #                                    "35.871", "38.327", "41.386", "45.812", "50.01", "53.71", "58.548", "63.903", "66.884", "68.92", 
    #                                    "75.776", "79.581", "85.142", "91.092", "97.439", "103.42", "108.42", "112.653", "121.198", "133.57", 
    #                                    "142.56", "152.919", "163.913", "192.116", "242.636", "292.002", "342.592", "392.635", "446.347", 
    #                                    "495.021", "544.526", "597.864", "671.43", "774.283", "874.801", "983.224", "1082.058", "1186.351", 
    #                                    "1288.68", "1412.564", "1528.811", "1834.117", "2298.427", "2820.119", "3588.558", "4565.85", "5073.16"],
    #"ECALOutputEnergyCorrectionPoints": ["0.1",  
    #                                     "16.0", "18.5", "22.5", "27.5", "32.5", "37.5", "42.5", "47.5", "52.5", "57.5", "62.5", "67.5", 
    #                                     "72.5", "77.5", "82.5", "87.5", "92.5", "97.5", "105.0", "115.0", "125.0", "135.0", "145.0", 
    #                                     "175.0", "225.0", "275.0", "325.0", "375.0", "425.0", "475.0", "525.0", "575.0", "650.0", "750.0", 
    #                                     "850.0", "950.0", "1050.0", "1150.0", "1250.0", "1350.0", "1450.0", "1750.0", "2250.0", "2750.0", 
    #                                     "3500.0", "4500.0", "5000."],
    "SplitVertexCollections": ["SplitVertices"],
    "StartVertexAlgorithmName": ["PandoraPFANew"],
    "StartVertexCollectionName": ["PandoraStartVertices"],
    "StripSplittingOn": ["0"],
    "TrackCollections": ["SiTracks_Refitted"],
    "TrackCreatorName": ["DDTrackCreatorCLIC"],
    "TrackStateTolerance": ["0"],
    "TrackSystemName": ["DDKalTest"],
    "UnmatchedVertexTrackMaxEnergy": ["5"],
    "UseEcalScLayers": ["0"],
    "UseNonVertexTracks": ["1"],
    "UseOldTrackStateCalculation": ["0"],
    "UseUnmatchedNonVertexTracks": ["0"],
    "UseUnmatchedVertexTracks": ["1"],
    "V0VertexCollections": ["V0Vertices"],
    "YokeBarrelNormalVector": ["0", "0", "1"],
    "Z0TrackCut": ["200"],
    "Z0UnmatchedVertexTrackCut": ["5"],
    "ZCutForNonVertexTracks": ["250"]
}

FastJetProcessor = MarlinProcessorWrapper("FastJetProcessor")
FastJetProcessor.OutputLevel = INFO
FastJetProcessor.ProcessorType = "FastJetProcessor"
FastJetProcessor.Parameters = {
    "algorithm": ["antikt_algorithm", "0.4"],
    "clusteringMode": ["Inclusive", "5"],
    "jetOut": ["JetOut"],
    "recParticleIn": ["PandoraPFOs"],
    "recombinationScheme": ["E_scheme"]
}

ValenciaJetProcessor = MarlinProcessorWrapper("ValenciaJetProcessor")
ValenciaJetProcessor.OutputLevel = INFO
ValenciaJetProcessor.ProcessorType = "FastJetProcessor"
ValenciaJetProcessor.Parameters = {
    "algorithm": ["ValenciaPlugin", "1.2", "1.0", "0.7"],
    "clusteringMode": ["ExclusiveNJets", "2"],
    "jetOut": ["ValenciaJetOut"],
    "recParticleIn": ["PandoraPFOs"],
    "recombinationScheme": ["E_scheme"]
}

MyDDSimpleMuonDigi = MarlinProcessorWrapper("MyDDSimpleMuonDigi")
MyDDSimpleMuonDigi.OutputLevel = INFO
MyDDSimpleMuonDigi.ProcessorType = "DDSimpleMuonDigi"
MyDDSimpleMuonDigi.Parameters = {
    "CalibrMUON": ["70.1"],
    "MUONCollections": ["YokeBarrelCollection", "YokeEndcapCollection"],
    "MUONOutputCollection": ["MUON"],
    "MaxHitEnergyMUON": ["2.0"],
    "MuonThreshold": ["1e-06"],
    "RelationOutputCollection": ["RelationMuonHit"]
}

OverlayMIX = MarlinProcessorWrapper("OverlayMIX")
OverlayMIX.OutputLevel = INFO
OverlayMIX.ProcessorType = "OverlayTimingRandomMix"
OverlayMIX.Parameters = {
    "PathToMuPlus": ["/data/userdata/phchang/muc/BIB10TeV/BIB10TeV/sim_mm/"],
    "PathToMuMinus": ["/data/userdata/phchang/muc/BIB10TeV/BIB10TeV/sim_mp/"],
    "Collection_IntegrationTimes": [
        "VertexBarrelCollection", "-0.18", "0.18",
        "VertexEndcapCollection", "-0.18", "0.18",
        "InnerTrackerBarrelCollection", "-0.36", "0.36",
        "InnerTrackerEndcapCollection", "-0.36", "0.36",
        "OuterTrackerBarrelCollection", "-0.36", "0.36",
        "OuterTrackerEndcapCollection", "-0.36", "0.36",
        "ECalBarrelCollection", "-0.5", "15.",
        "ECalEndcapCollection", "-0.5", "15.",
        "HCalBarrelCollection", "-0.5", "15.",
        "HCalEndcapCollection", "-0.5", "15.",
        "YokeBarrelCollection", "-0.5", "15.",
        "YokeEndcapCollection", "-0.5", "15."
    ],
    "IntegrationTimeMin": ["-0.5"],
    "MCParticleCollectionName": ["MCParticle"],
    "MergeMCParticles": ["false"],
    "NumberBackground": ["1666"]
    # "NumberBackground": ["10"]
}


OverlayIP = MarlinProcessorWrapper("OverlayIP")
OverlayIP.OutputLevel = INFO
OverlayIP.ProcessorType = "OverlayTimingGeneric"
OverlayIP.Parameters = {
    "AllowReusingBackgroundFiles": ["true"],
    "BackgroundFileNames": [
        "/dataMuC/IPairs/sim/sim_pairs_cycle1.slcio",
        "/dataMuC/IPairs/sim/sim_pairs_cycle2.slcio",
        "/dataMuC/IPairs/sim/sim_pairs_cycle3.slcio",
        "/dataMuC/IPairs/sim/sim_pairs_cycle4.slcio"
    ],
    "Collection_IntegrationTimes": [
        "VertexBarrelCollection", "-0.18", "0.18",
        "VertexEndcapCollection", "-0.18", "0.18",
        "InnerTrackerBarrelCollection", "-0.36", "0.36",
        "InnerTrackerEndcapCollection", "-0.36", "0.36",
        "OuterTrackerBarrelCollection", "-0.36", "0.36",
        "OuterTrackerEndcapCollection", "-0.36", "0.36",
        "ECalBarrelCollection", "-0.5", "15.",
        "ECalEndcapCollection", "-0.5", "15.",
        "HCalBarrelCollection", "-0.5", "15.",
        "HCalEndcapCollection", "-0.5", "15.",
        "YokeBarrelCollection", "-0.5", "15.",
        "YokeEndcapCollection", "-0.5", "15."
    ],
    "Delta_t": ["10000"],
    "IntegrationTimeMin": ["-0.5"],
    "MCParticleCollectionName": ["MCParticle"],
    "MCPhysicsParticleCollectionName": ["MCPhysicsParticles_IP"],
    "MergeMCParticles": ["false"],
    "NBunchtrain": ["1"],
    "NumberBackground": ["1"],
    "PhysicsBX": ["1"],
    "Poisson_random_NOverlay": ["false"],
    "RandomBx": ["false"],
    "StartBackgroundFileIndex": ["0"],
    "TPCDriftvelocity": ["0.05"]
}

algList.append(MyAIDAProcessor)
algList.append(EventNumber)
algList.append(InitDD4hep)
if the_args.enableBIB:
    algList.append(OverlayMIX)
if the_args.enableIP:
    algList.append(OverlayIP)
algList.append(VXDBarrelDigitiser)
algList.append(VXDEndcapDigitiser)
algList.append(InnerPlanarDigiProcessor)
algList.append(InnerEndcapPlanarDigiProcessor)
algList.append(OuterPlanarDigiProcessor)
algList.append(OuterEndcapPlanarDigiProcessor)
algList.append(VXDBarrelConer)
algList.append(VXDEndcapConer)
algList.append(InnerPlanarConer)
algList.append(InnerEndcapConer)
algList.append(OuterPlanarConer)
algList.append(OuterEndcapConer)
#algList.append(CKFTracking)
#algList.append(TrackDeduper)
#algList.append(Refit)
##algList.append(MyTrackSelector)
#algList.append(MyTrackTruth)
#algList.append(MyEcalBarrelDigi)
#algList.append(MyEcalBarrelReco)
#algList.append(MyEcalEndcapDigi)
#algList.append(MyEcalEndcapReco)
#algList.append(MyHcalBarrelDigi)
#algList.append(MyHcalBarrelReco)
#algList.append(MyHcalEndcapDigi)
#algList.append(MyHcalEndcapReco)
#algList.append(MyEcalBarrelConer)
#algList.append(MyEcalEndcapConer)
#algList.append(MyHcalBarrelConer)
#algList.append(MyHcalEndcapConer)
#algList.append(MyEcalBarrelSelector)
#algList.append(MyEcalEndcapSelector)
#algList.append(MyDDSimpleMuonDigi)
#algList.append(DDMarlinPandora)
#algList.append(FastJetProcessor)
##algList.append(ValenciaJetProcessor)
algList.append(Output_REC)

ApplicationMgr(TopAlg=algList,
               EvtSel='NONE',
               EvtMax=the_args.NEvents,
               ExtSvc=[evtsvc],
               OutputLevel=INFO
               )

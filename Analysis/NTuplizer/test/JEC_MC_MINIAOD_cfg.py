# -*- coding: utf-8 -*-

import FWCore.ParameterSet.Config as cms
import sys

from RecoJets.Configuration.RecoPFJets_cff import *
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets, ak8GenJets
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff import *
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import selectedPatJets
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.patSequences_cff import *
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.slimming.metFilterPaths_cff import *


## Modified version of jetToolBox from https://github.com/cms-jet/jetToolbox
## Options for PUMethod: Puppi, CS, SK, CHS

# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("Combined")
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Conditions
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')
process.load('RecoJets.JetProducers.TrackJetParameters_cfi')
process.load('RecoJets.JetProducers.PileupJetIDParams_cfi')

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")


from RecoJets.Configuration.GenJetParticles_cff import *

process.GlobalTag.globaltag = "106X_upgrade2018_realistic_v15_L1v1"
#process.GlobalTag.globaltag = "102X_upgrade2018_realistic_v20"
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

from PhysicsTools.PatAlgos.tools.coreTools import *
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Input
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

inFiles = cms.untracked.vstring(
#'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/QCD_bEnriched_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/20000/105691C7-8434-2940-88D2-7369139692CB.root'   
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer19UL18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/260000/0112B4DB-39FA-4645-A39D-912087A8C335.root'
#'/store/mc/RunIISummer19UL18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/260000/00C28834-56C0-2343-B436-AA8521756E9E.root'
#'/store/mc/RunIISummer19UL18MiniAOD/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/100000/00BC2C64-EBF8-E843-ACC9-A1A2744FE11B.root'
#'file:/tmp/deroy/00C28834-56C0-2343-B436-AA8521756E9E.root'
#'root://xrootd-cms.infn.it//store/mc/RunIISummer19UL18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/260000/00C28834-56C0-2343-B436-AA8521756E9E.root'
   )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(3000))

#process.firstEvent = cms.untracked.PSet(input = cms.untracked.int32(5000))
process.source = cms.Source("PoolSource", fileNames = inFiles )


process.p = cms.Path()

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Services
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('CommonTools.UtilAlgos.TFileService_cfi')

process.TFileService = cms.Service("TFileService",
fileName = cms.string('hist.root')             #largest data till April5,2016 
)

process.patJets.addTagInfos = True

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# Electron IDs for AOD/MINIAOD
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

# define which IDs to produce
el_id_modules = [
##    "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff"
    "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff",
    "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff"
##    "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff"
]
# Add them to the VID producer
for iModule in el_id_modules:

	setupAllVIDIdsInModule(process, iModule, setupVIDElectronSelection)

switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)

pho_id_modules = [
	"RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff",
	"RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff"
]

for iModule in pho_id_modules:

	setupAllVIDIdsInModule(process, iModule, setupVIDPhotonSelection)


from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
deep_discriminators = ["pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD",
                       "pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD",
                       "pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD",
		       "pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHbbvsQCD",
		       "pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:bbvsLight"
                       
]

from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetJetTagsAll as pfParticleNetJetTagsAll
deep_discriminators += pfParticleNetJetTagsAll

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJetsAK8'),
   pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
   svSource = cms.InputTag('slimmedSecondaryVertices'),
   rParam = 0.8,
   labelName = 'SlimmedJetsAK8',
   jetCorrections = ('AK8PFPuppi', cms.vstring([]), 'None' ),
   btagDiscriminators = deep_discriminators 
)

# For prefire correction #

from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
process.prefiringweight = l1PrefiringWeightProducer.clone(
	TheJets = cms.InputTag("slimmedJets"), #"updatedPatJetsUpdatedJEC"), #this should be the slimmedJets collection with up to date JECs 
	DataEraECAL = cms.string("2017BtoF"), #Use 2016BtoH for 2016
	DataEraMuon = cms.string("20172018"), #Use 2016 for 2016
	UseJetEMPt = cms.bool(False),
	PrefiringRateSystematicUnctyECAL = cms.double(0.2),
	PrefiringRateSystematicUnctyMuon = cms.double(0.2)
)

# Analyzer #

process.mcjets =  cms.EDAnalyzer('Leptop',

	 Data =  cms.untracked.bool(False),
	 MonteCarlo =  cms.untracked.bool(True),
         YEAR = cms.untracked.int32(2018),
         UltraLegacy =  cms.untracked.bool(True),                        
	 isReco = cms.untracked.bool(True),
 	 ReRECO = cms.untracked.bool(True),
	 SoftDrop_ON =  cms.untracked.bool(True),
	 add_prefireweights =  cms.untracked.bool(False),

 	 RootFileName = cms.untracked.string('rootuple.root'),  #largest data till April5,2016
	
#	 PFJetsAK8 = cms.InputTag("slimmedJetsAK8"),
#        PFJetsAK8 = cms.InputTag("selectedPatJetsAK8PFPuppiSoftDropPacked","SubJets","Combined"),
#        PFJetsAK8 = cms.InputTag("selectedPatJetsAK8PFPuppi","","Combined"),
	 PFJetsAK8 = cms.InputTag("updatedPatJetsTransientCorrectedSlimmedJetsAK8"),#"","PAT"),
	 AK8PtCut = cms.untracked.double(180.),
         AK8GenPtCut = cms.untracked.double(150.),
	 softdropmass  = cms.untracked.string("ak8PFJetsSoftDropMass"),#ak8PFJetsPuppiSoftDropMass"),#('ak8PFJetsPuppiSoftDropMass'),
	 tau1  = cms.untracked.string("NjettinessAK8Puppi:tau1"),#'NjettinessAK8Puppi:tau1'),
	 tau2  = cms.untracked.string("NjettinessAK8Puppi:tau2"),#'NjettinessAK8Puppi:tau2'),
	 tau3  = cms.untracked.string("NjettinessAK8Puppi:tau3"),#'NjettinessAK8Puppi:tau3'),
	 subjets  = cms.untracked.string('SoftDropPuppi'),#("SoftDrop"),#'SoftDropPuppi'),
#	 subjets = cms.untracked.string("slimmedJetsAK8PFPuppiSoftDropPacked"),
	 toptagger = cms.untracked.string("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD"),
	 Wtagger = cms.untracked.string("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD"),
	 Ztagger = cms.untracked.string("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD"),
	 Htagger = cms.untracked.string("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHbbvsQCD"),
	 bbtagger = cms.untracked.string("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:bbvsLight"),	

	 PFJetsAK4 = cms.InputTag("slimmedJets"), 
	 minjPt = cms.untracked.double(25.),
	 maxEta = cms.untracked.double(3.),

	 GENJetAK8 = cms.InputTag("slimmedGenJetsAK8"),
	 GENJetAK4 = cms.InputTag("slimmedGenJets"),
         minGenPt = cms.untracked.double(15.),      
	 maxGenEta = cms.untracked.double(5.),         
         
	 Muons = cms.InputTag("slimmedMuons"),
	 minmuPt  = cms.untracked.double(10.),
	 EAFile_MiniIso = cms.FileInPath("PhysicsTools/NanoAOD/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_94X.txt"),

	 Electrons = cms.InputTag("slimmedElectrons"),
	 minePt   = cms.untracked.double(10.),
	 electronID_isowp90        = cms.string('mvaEleID-Fall17-iso-V2-wp90'),
         electronID_noisowp90      = cms.string('mvaEleID-Fall17-noIso-V2-wp90'),
         electronID_isowp80        = cms.string('mvaEleID-Fall17-iso-V2-wp80'),
         electronID_noisowp80      = cms.string('mvaEleID-Fall17-noIso-V2-wp80'),

	 Photons = cms.InputTag("slimmedPhotons"),
	 mingmPt  = cms.untracked.double(20.),
	 PhoID_FallV2_WP90 = cms.string("mvaPhoID-RunIIFall17-v2-wp90"),
         PhoID_FallV2_WP80 = cms.string("mvaPhoID-RunIIFall17-v2-wp80"),
         PhoID_SpringV1_WP90 = cms.string("mvaPhoID-Spring16-nonTrig-V1-wp90"),
         PhoID_SpringV1_WP80 = cms.string("mvaPhoID-Spring16-nonTrig-V1-wp80"),
	 label_mvaPhoID_FallV2_Value = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRunIIFall17v2Values"),
	 #label_mvaPhoID_FallV2_WP90 = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIFall17-v2-wp90"),
         #label_mvaPhoID_FallV2_WP80 = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIFall17-v2-wp80"),

	 Taus = cms.InputTag("slimmedTausUpdated"),
	 mintauPt = cms.untracked.double(25.),

	 PFMet = cms.InputTag("slimmedMETs"),
         PuppiMet = cms.InputTag("slimmedMETsPuppi"),
         GENMet  = cms.InputTag("genMetTrue","","SIM"),

	 GenParticles = cms.InputTag("prunedGenParticles"),#("prunedGenParticles"),#("packedGenParticles"),
	 jetFlavourInfos = cms.InputTag("slimmedGenJetsFlavourInfos"),
	 #jetFlavourInfos = cms.InputTag("genJetAK8FlavourAssociation"),

	 pfCands = cms.InputTag("packedPFCandidates"),

         Beamspot = cms.InputTag("offlineBeamSpot"),
         PrimaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
         SecondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
	 slimmedAddPileupInfo = cms.InputTag("slimmedAddPileupInfo"),
         PFRho = cms.InputTag("fixedGridRhoFastjetAll"),

	 Generator = cms.InputTag("generator"),
         LHEEventProductInputTag = cms.InputTag('externalLHEProducer'),
         GenEventProductInputTag = cms.InputTag('generator'),
	 nPDFsets = cms.untracked.uint32(103),

	 jecL1FastFileAK4          = cms.string('Summer19UL18_V5_MC/Summer19UL18_V5_MC_L1FastJet_AK4PFchs.txt'),
         jecL1FastFileAK8          = cms.string('Summer19UL18_V5_MC/Summer19UL18_V5_MC_L1FastJet_AK8PFPuppi.txt'),
         jecL2RelativeFileAK4      = cms.string('Summer19UL18_V5_MC/Summer19UL18_V5_MC_L2Relative_AK4PFchs.txt'),
         jecL2RelativeFileAK8      = cms.string('Summer19UL18_V5_MC/Summer19UL18_V5_MC_L2Relative_AK8PFPuppi.txt'),
         jecL3AbsoluteFileAK4      = cms.string('Summer19UL18_V5_MC/Summer19UL18_V5_MC_L3Absolute_AK4PFchs.txt'),
         jecL3AbsoluteFileAK8      = cms.string('Summer19UL18_V5_MC/Summer19UL18_V5_MC_L3Absolute_AK8PFPuppi.txt'),
         jecL2L3ResidualFileAK4    = cms.string('Summer19UL18_V5_MC/Summer19UL18_V5_MC_L2L3Residual_AK4PFchs.txt'),
         jecL2L3ResidualFileAK8    = cms.string('Summer19UL18_V5_MC/Summer19UL18_V5_MC_L2L3Residual_AK8PFPuppi.txt'),

	 PtResoFileAK4  = cms.string('Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt'),
         PtResoFileAK8  = cms.string('Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_PtResolution_AK8PFPuppi.txt'),
         PtSFFileAK4 = cms.string('Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_SF_AK4PFchs.txt'),
         PtSFFileAK8 = cms.string('Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_SF_AK8PFPuppi.txt'),

         JECUncFileAK4 = cms.string("Summer19UL18_V5_MC/Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.txt"),
	 JECUncFileAK8 = cms.string("Summer19UL18_V5_MC/Summer19UL18_V5_MC_UncertaintySources_AK8PFPuppi.txt"),

	 bits = cms.InputTag("TriggerResults","","HLT"),
         prescales = cms.InputTag("patTrigger","","RECO"),
         objects = cms.InputTag("slimmedPatTrigger")
)

#===== MET Filters ==

process.load('RecoMET.METFilters.primaryVertexFilter_cfi')
process.primaryVertexFilter.vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
process.load('RecoMET.METFilters.globalSuperTightHalo2016Filter_cfi')
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')
process.HBHENoiseFilterResultProducerNoMinZ = process.HBHENoiseFilterResultProducer.clone(minZeros = cms.int32(99999))
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonFilter.vtx = cms.InputTag("offlineSlimmedPrimaryVertices") 
process.BadPFMuonFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.BadPFMuonDzFilter_cfi')
process.BadPFMuonDzFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonDzFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonDzFilter.vtx = cms.InputTag("offlineSlimmedPrimaryVertices")
process.BadPFMuonDzFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag('reducedEgamma','reducedEERecHits')

process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('eventoutput.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.allMetFilterPaths=cms.Sequence(process.primaryVertexFilter*process.globalSuperTightHalo2016Filter*process.HBHENoiseFilter*process.HBHENoiseIsoFilter*process.EcalDeadCellTriggerPrimitiveFilter*process.BadPFMuonFilter*process.BadPFMuonDzFilter*process.ecalBadCalibFilter)

process.jetSeq=cms.Sequence(process.patJetCorrFactorsSlimmedJetsAK8+process.updatedPatJetsSlimmedJetsAK8+process.patJetCorrFactorsTransientCorrectedSlimmedJetsAK8+process.pfDeepBoostedJetTagInfosSlimmedJetsAK8+process.pfMassDecorrelatedDeepBoostedJetTagsSlimmedJetsAK8+process.pfMassDecorrelatedDeepBoostedDiscriminatorsJetTagsSlimmedJetsAK8+process.updatedPatJetsTransientCorrectedSlimmedJetsAK8)

import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
updatedTauName = "slimmedTausUpdated"
tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms, debug = True, updatedTauName = updatedTauName,
	toKeep = [ "2017v2",
        #       "dR0p32017v2", 
        #       "newDM2017v2", 
        "deepTau2017v2p1",
        "againstEle2018"\
        ]
)
tauIdEmbedder.runTauID()


process.p = cms.Path(process.egmPhotonIDSequence* 
 		     process.HBHENoiseFilterResultProducer*process.HBHENoiseFilterResultProducerNoMinZ*
		     process.allMetFilterPaths*
#		     process.egmGsfElectronIDSequence*
		     process.rerunMvaIsolationSequence*getattr(process,updatedTauName)* #process.slimmedTausUpdated*  # this also works
		     process.jetSeq *
		     process.prefiringweight *
		     process.mcjets)

process.schedule = cms.Schedule(process.p)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)

#process.options.numberOfThreads=cms.untracked.uint32(2)
#process.options.numberOfStreams=cms.untracked.uint32(0)

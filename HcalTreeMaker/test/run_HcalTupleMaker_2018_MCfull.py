#------------------------------------------------------------------------------------
# Imports
#------------------------------------------------------------------------------------
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
import FWCore.ParameterSet.VarParsing as VarParsing

#------------------------------------------------------------------------------------
# Declare the process and input variables
#------------------------------------------------------------------------------------
#process = cms.Process('NOISE',eras.Run2_50ns)#for 50ns 13 TeV data
#process = cms.Process('NOISE',eras.Run2_25ns)#for 25ns 13 TeV data
options = VarParsing.VarParsing ('analysis')
process = cms.Process("Trees",eras.Run2_2018)

options.register ('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "no of skipped events")
#
# Dataset
# /RelValSinglePiPt25Eta1p7_2p7/CMSSW_9_3_2-93X_upgrade2023_realistic_v2_2023D17noPU-v1/GEN-SIM-RECO
# /SinglePiPt*Eta1p6_2p8/PhaseIITDRFall17*93X_upgrade2023_realistic_v2*/GEN-SIM-RECO
#
# pt=50 GeV sample
options.inputFiles = '/store/relval/CMSSW_10_2_0_pre1/RelValTTbar_13/GEN-SIM-RECO/101X_upgrade2018_realistic_v7-v1/10000/BCE80CB8-F13C-E811-A56A-0CC47A4C8F06.root','/store/relval/CMSSW_10_2_0_pre1/RelValTTbar_13/GEN-SIM-RECO/101X_upgrade2018_realistic_v7-v1/10000/24C51D80-F23C-E811-ABC9-0025905B8582.root','/store/relval/CMSSW_10_2_0_pre1/RelValTTbar_13/GEN-SIM-RECO/101X_upgrade2018_realistic_v7-v1/10000/F04155FA-F73C-E811-87D9-0CC47A7C347E.root','/store/relval/CMSSW_10_2_0_pre1/RelValTTbar_13/GEN-SIM-RECO/101X_upgrade2018_realistic_v7-v1/10000/3C845010-F93C-E811-ADEB-0025905A6076.root'
options.secondaryInputFiles = '/store/relval/CMSSW_10_2_0_pre1/RelValTTbar_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/10000/08F9E721-EC3C-E811-B87E-0CC47A4D766C.root','/store/relval/CMSSW_10_2_0_pre1/RelValTTbar_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/10000/407C5724-EC3C-E811-BF1A-0CC47A7C347E.root','/store/relval/CMSSW_10_2_0_pre1/RelValTTbar_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/10000/04DE4232-EB3C-E811-AAC5-0025905A606A.root','/store/relval/CMSSW_10_2_0_pre1/RelValTTbar_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/10000/BA873832-EB3C-E811-B4E7-0025905A48F2.root','/store/relval/CMSSW_10_2_0_pre1/RelValTTbar_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/10000/E63A1431-EB3C-E811-AB84-0025905A48D8.root','/store/relval/CMSSW_10_2_0_pre1/RelValTTbar_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/10000/0CE933F9-F23C-E811-97A3-0025905A612E.root','/store/relval/CMSSW_10_2_0_pre1/RelValTTbar_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/10000/9277AFC4-F13C-E811-8AB8-0025905B85B6.root'
options.outputFile = 'relval_ttbar_2018_MCfull.root'
#
options.maxEvents = 10 # -1 means all events
#options.skipEvents = 0 # default is 0.

#------------------------------------------------------------------------------------
# Get and parse the command line arguments
#------------------------------------------------------------------------------------
options.parseArguments()
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(options.secondaryInputFiles),
    skipEvents = cms.untracked.uint32(options.skipEvents) # default is 0.
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string(options.outputFile)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    Rethrow = cms.untracked.vstring("ProductNotFound"), # make this exception fatal
    fileMode  =  cms.untracked.string('NOMERGE') # no ordering needed, but calls endRun/beginRun etc. at file boundaries
)

#------------------------------------------------------------------------------------
# import of standard configurations
#------------------------------------------------------------------------------------
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#KH
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#------------------------------------------------------------------------------------
# Set up our analyzer
#------------------------------------------------------------------------------------
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_Tree_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_Event_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_GenParticles_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HBHERecHits_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HcalSimHits_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HBHEDigis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HFDigis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HODigis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HcalTriggerPrimitives_cfi")

#------------------------------------------------------------------------------------
# Specify Global Tag
#------------------------------------------------------------------------------------
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')

#----------------------------
# Paths/Sequences Definitions
#----------------------------
#process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")
#process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")

process.digiPath = cms.Path(
    process.hcalDigis
)

#------------------------------------------------------------------------------------
# HcalTupleMaker sequence definition
#------------------------------------------------------------------------------------
process.tuple_step = cms.Sequence(
    # Make HCAL tuples: Event, run, ls number
    process.hcalTupleEvent*
    # Make HCAL tuples: digi info
    process.hcalTupleHBHEDigis*
    process.hcalTupleHODigis*
    process.hcalTupleHFDigis*
    #process.hcalTupleTriggerPrimitives*
    # Make HCAL tuples: reco info
    process.hcalTupleHBHERecHits*
    # Make HCAL tuples: simhit info
    process.hcalTupleHcalSimHits*
    process.hcalTupleGenParticles*
    #
    process.hcalTupleTree
)


#-----------------------------------------------------------------------------------
# Path and EndPath definitions
#-----------------------------------------------------------------------------------
process.preparation = cms.Path(
    process.tuple_step
)

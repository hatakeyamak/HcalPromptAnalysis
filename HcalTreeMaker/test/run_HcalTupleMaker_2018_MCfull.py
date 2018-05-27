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
# Dataset e.g.
# dasgoclient --query 'dataset dataset=/RelValTTbar_13/CMSSW_10_2_0_pre1-*realistic*/GEN-SIM-RECO'                 
# dasgoclient --query 'file dataset=/RelValTTbar_13/CMSSW_10_2_0_pre3-101X_upgrade2018_realistic_v7-v1/GEN-SIM-RECO'
#
# MinBias sample
options.inputFiles = 'root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-RECO/101X_upgrade2018_realistic_v7-v1/20000/3687EEB3-F952-E811-B04C-0025905A60F4.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-RECO/101X_upgrade2018_realistic_v7-v1/20000/92ED646F-FB52-E811-8D13-0CC47A7C3572.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-RECO/101X_upgrade2018_realistic_v7-v1/20000/3A0F2FA8-FC52-E811-8AAB-0CC47A4C8E20.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-RECO/101X_upgrade2018_realistic_v7-v1/20000/7C7C9FBA-FC52-E811-A146-0CC47A4C8E2A.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-RECO/101X_upgrade2018_realistic_v7-v1/20000/02FCE7FE-FB52-E811-9137-0CC47A745294.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-RECO/101X_upgrade2018_realistic_v7-v1/20000/D2B0744B-FD52-E811-9B60-0025905B85DA.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-RECO/101X_upgrade2018_realistic_v7-v1/20000/C2273091-FE52-E811-8E2B-0025905A610A.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-RECO/101X_upgrade2018_realistic_v7-v1/20000/CAE04C7A-0053-E811-A27A-0CC47A74524E.root'
options.secondaryInputFiles = 'root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/7E7ED8C4-F152-E811-B68E-0CC47A7C3638.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/6E3E8A31-F152-E811-A9B7-0CC47A4D7698.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/ACACA7CB-F152-E811-862B-0CC47A7C3612.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/D2F45456-F152-E811-A20D-0025905B85C0.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/F23F62AA-F252-E811-802C-0CC47A4D768C.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/E8086532-F452-E811-BED9-0CC47A4D767A.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/A0ECBA50-F352-E811-B557-0CC47A4D76C8.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/F2EB0C58-F352-E811-82DF-002618FDA265.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/9AFBAA52-F352-E811-A782-0CC47A78A3EC.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/DC521D51-F352-E811-9241-0CC47A4D768C.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/B20B32E0-F252-E811-A75F-0025905A60A6.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/FE2FB95B-F352-E811-BE41-0CC47A4D7668.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/367807E4-F152-E811-89F5-0CC47A78A3F8.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/D2F000D5-F252-E811-89C6-0CC47A74525A.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/F0351D69-F352-E811-A411-0CC47A78A360.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/D83FC9EF-F152-E811-B224-0025905B859A.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/5231F3FA-F152-E811-A26E-0025905B8604.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/4EB527AA-F752-E811-B678-0CC47A7C353E.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/1658DC6B-F352-E811-B7C4-0025905A611C.root','root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_2_0_pre3/RelValMinBias_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/1461CE62-F452-E811-9F47-0025905A4964.root'
options.outputFile = 'trees_relval_minbias_2018_MCfull.root'
#
#options.maxEvents = -1 # -1 means all events
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
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_QIE11Digis_cfi")
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
    process.hcalTupleQIE11Digis*
    #process.hcalTupleTriggerPrimitives*
    # Make HCAL tuples: reco info
    process.hcalTupleHBHERecHits*
    # Make HCAL tuples: gen & simhit info
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

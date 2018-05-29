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
#process = cms.Process("Trees",eras.Run2_2018)
process = cms.Process('Trees',eras.run2_HCAL_2017, eras.run2_HF_2017,eras.run2_HEPlan1_2017)

options.register ('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "no of skipped events")
#
# Dataset e.g.
# dasgoclient --query 'dataset dataset=/RelValTTbar_13/CMSSW_10_2_0_pre3-*realistic*/GEN-SIM-RECO'                 
# dasgoclient --query 'file dataset=/RelValTTbar_13/CMSSW_10_2_0_pre3-101X_upgrade2018_realistic_v7-v1/GEN-SIM-RECO'
#
# MinBias sample
options.inputFiles = 'root://kodiak-se.baylor.edu//store/data/Run2017A/HLTPhysics/RAW/v1/000/296/641/00000/1261D3DA-554F-E711-B1F5-02163E014342.root','root://kodiak-se.baylor.edu//store/data/Run2017A/HLTPhysics/RAW/v1/000/296/641/00000/4E56A0DC-554F-E711-AAA4-02163E011DD1.root','root://kodiak-se.baylor.edu//store/data/Run2017A/HLTPhysics/RAW/v1/000/296/642/00000/C6B5DAAB-604F-E711-87C8-02163E0145FE.root','root://kodiak-se.baylor.edu//store/data/Run2017A/HLTPhysics/RAW/v1/000/296/643/00000/F4E827C1-5B4F-E711-BD2C-02163E0144B1.root','root://kodiak-se.baylor.edu//store/data/Run2017A/HLTPhysics/RAW/v1/000/296/644/00000/705678DB-614F-E711-9531-02163E0146BC.root','root://kodiak-se.baylor.edu//store/data/Run2017A/HLTPhysics/RAW/v1/000/296/646/00000/B2506154-644F-E711-B4AF-02163E012A81.root','root://kodiak-se.baylor.edu//store/data/Run2017A/HLTPhysics/RAW/v1/000/296/647/00000/063A2991-724F-E711-80C1-02163E011B1A.root','root://kodiak-se.baylor.edu//store/data/Run2017A/HLTPhysics/RAW/v1/000/296/647/00000/624A76E2-714F-E711-9131-02163E014145.root','root://kodiak-se.baylor.edu//store/data/Run2017A/HLTPhysics/RAW/v1/000/296/647/00000/94B47B4D-6A4F-E711-B2A7-02163E01256B.root','root://kodiak-se.baylor.edu//store/data/Run2017A/HLTPhysics/RAW/v1/000/296/647/00000/C80F6034-6A4F-E711-A688-02163E013768.root'
#options.secondaryInputFiles = 
options.outputFile = 'trees_2017_HLTPhysics_DataFull.root'
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
#KH process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
#KH process.load('Configuration.StandardSequences.PATMC_cff')
#KH process.load('Configuration.StandardSequences.Validation_cff')
#KH process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
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

process.hcalTupleHBHERecHits.source = cms.untracked.InputTag("hbheprereco")
process.hcalTupleHBHEDigis.recHits = cms.untracked.InputTag("hbheprereco")

#------------------------------------------------------------------------------------
# Specify Global Tag
#------------------------------------------------------------------------------------
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data_promptlike', '')

#----------------------------
# Paths/Sequences Definitions
#----------------------------
#process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")
#process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")

process.digiPath = cms.Path(
    process.hcalDigis
)

process.recoPath = cms.Path(
    process.horeco
    *process.hfprereco
    *process.hfreco
    *process.hbheprereco
    *process.hbheplan1
    # *process.MEtoEDMConverter #KH
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
    #process.hcalTupleHcalSimHits*
    #process.hcalTupleGenParticles*
    #
    process.hcalTupleTree
)


#-----------------------------------------------------------------------------------
# Path and EndPath definitions
#-----------------------------------------------------------------------------------
process.preparation = cms.Path(
    process.tuple_step
)

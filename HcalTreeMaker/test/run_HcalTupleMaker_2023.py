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
process = cms.Process("Trees",eras.Phase2)

options.register ('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "no of skipped events")
options.inputFiles = '/store/relval/CMSSW_9_3_2/RelValSinglePiPt25Eta1p7_2p7/GEN-SIM-RECO/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/2637F672-C7A6-E711-B4EF-0025905A612A.root'
options.outputFile = 'results.root'
options.maxEvents = 10 # -1 means all events
#options.skipEvents = 0 # default is 0.

#------------------------------------------------------------------------------------
# Get and parse the command line arguments
#------------------------------------------------------------------------------------
options.parseArguments()
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
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
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')  # <=== to be checked
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
#process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_cfi") # Dont want to use this, load modules individually
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_Tree_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_Event_cfi")
#process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HBHEDigis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HBHERecHits_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HGCRecHits_cfi")
#process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_Trigger_cfi")

#------------------------------------------------------------------------------------
# Specify Global Tag
#------------------------------------------------------------------------------------
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

#------------------------------------------------------------------------------------
# HcalTupleMaker sequence definition
#------------------------------------------------------------------------------------
process.tuple_step = cms.Sequence(
    # Make HCAL tuples: Event, run, ls number
    process.hcalTupleEvent*
    #    # Make HCAL tuples: digi info
    #    process.hcalTupleHBHEDigis*
    #    process.hcalTupleHODigis*
    #    process.hcalTupleHFDigis*
    #    process.hcalTupleTriggerPrimitives*
    #    # Make HCAL tuples: reco info
    process.hcalTupleHBHERecHits*
    process.hcalTupleHGCRecHits*
    #
    process.hcalTupleTree
)


#-----------------------------------------------------------------------------------
# Path and EndPath definitions
#-----------------------------------------------------------------------------------
process.preparation = cms.Path(
    #process.my_hlt *
    #process.RawToDigi * #needed for RAW files
    #process.L1Reco *
    #rprocess.reconstruction * #needed for RAW files
    #process.caloglobalreco *
    #process.reconstructionCosmics *
    #
    #process.horeco *
    #process.hfreco *
    #
    #process.hbheprerecoMethod0 *
    #process.hbheprerecoMethod2 *
    #process.hbherecoMethod0 *
    #process.hbherecoMethod2 *
    #
    #process.towerMakerMethod0 *
    #process.towerMakerMethod2 *
    #
    #process.hcalnoiseMethod0 *
    #process.hcalnoiseMethod2 *
    #
    #process.HBHENoiseFilterResultProducerMethod0 *
    #process.HBHENoiseFilterResultProducerMethod2 *
    #
    #
    #process.hcalCosmicDigis *
    #process.hcalL1JetDigis *
    #
    #process.hcalnoise *  #needed for RAW files
    #process.HBHENoiseFilterResultProducer *
    #process.ApplyBaselineHBHENoiseFilter *
    #
    process.tuple_step
)

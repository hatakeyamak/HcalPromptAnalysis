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



##
## Setup command line options
##
options.register ('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "no of skipped events")
options.register ('isMINIAOD', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "MINIAODSIM input file(s)?")

##
## Default
##
options.maxEvents = 10 # -1 means all events
#options.skipEvents = 0 # default is 0.

##
## get and parse the command line arguments
##
options.parseArguments()
print("isMINIAOD: ", options.isMINIAOD)
print("maxEvents: ", options.maxEvents)

#
# Dataset e.g.
# dasgoclient --query 'dataset dataset=/RelValTTbar_13/CMSSW_10_2_0_pre3-*realistic*/GEN-SIM-RECO'                 
# dasgoclient --query 'file dataset=/RelValTTbar_13/CMSSW_10_2_0_pre3-101X_upgrade2018_realistic_v7-v1/GEN-SIM-RECO'
#
# TTbar sample
#
# MINIAODSIM
if options.isMINIAOD: 
    options.inputFiles = open('RelValNuGun_Default_MiniAOD.txt').readlines()
    #options.inputFiles = open('RelValMinBias_Calib.txt').readlines()
    #'/store/relval/CMSSW_10_3_0_pre4/RelValTTbar_13/MINIAODSIM/PUpmx25ns_103X_upgrade2018_realistic_v4-v1/20000/17E57223-3406-7340-B867-2DDC36E7C371.root'
    options.outputFile = 'relval_nugun_2018_pmx25ns_miniaodsim.root'
# GEN-SIM-RECO
else:
    options.inputFiles = open('RelValNuGun_Default.txt').readlines()
    #options.inputFiles = '/store/relval/CMSSW_10_4_0_pre3/RelValMinBias_13/GEN-SIM-RECO/103X_upgrade2018_realistic_v8-v1/20000/F5A3B356-D2B7-9D44-819E-AD2E300601D5.root','/store/relval/CMSSW_10_4_0_pre3/RelValMinBias_13/GEN-SIM-RECO/103X_upgrade2018_realistic_v8-v1/20000/EDF4FA12-68A1-2949-A4DF-314BE46B57D5.root','/store/relval/CMSSW_10_4_0_pre3/RelValMinBias_13/GEN-SIM-RECO/103X_upgrade2018_realistic_v8-v1/20000/A620CAFF-ACF6-1042-B6CB-8B2A5EE7AFEC.root','/store/relval/CMSSW_10_4_0_pre3/RelValMinBias_13/GEN-SIM-RECO/103X_upgrade2018_realistic_v8-v1/20000/97346745-2D06-FE49-9335-6BA4DF5E35F5.root','/store/relval/CMSSW_10_4_0_pre3/RelValMinBias_13/GEN-SIM-RECO/103X_upgrade2018_realistic_v8-v1/20000/763194EA-CD2C-5E49-8CC6-5FCD144DBEF7.root','/store/relval/CMSSW_10_4_0_pre3/RelValMinBias_13/GEN-SIM-RECO/103X_upgrade2018_realistic_v8-v1/20000/57829579-EA9D-B14D-A6D7-8364B31721B8.root','/store/relval/CMSSW_10_4_0_pre3/RelValMinBias_13/GEN-SIM-RECO/103X_upgrade2018_realistic_v8-v1/20000/571B56F0-275C-BB4B-918D-E484970D065C.root','/store/relval/CMSSW_10_4_0_pre3/RelValMinBias_13/GEN-SIM-RECO/103X_upgrade2018_realistic_v8-v1/20000/45C8B284-6695-2A44-9471-045DF2A560CC.root'
    options.outputFile = 'relval_nugun_pmix_2018.root'
#

#
#
#
print("maxEvents: ", options.maxEvents)
print("inputFiles: ", options.inputFiles)
print("outputFile: ", options.outputFile)

#------------------------------------------------------------------------------------
# Get and parse the command line arguments
#------------------------------------------------------------------------------------
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
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('Configuration.StandardSequences.PATMC_cff')
#process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
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
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HBHEDigis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HFDigis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HODigis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HcalTriggerPrimitives_cfi")

process.load("HcalPromptAnalysis.HcalTreeMaker.TupleMaker_PFCandidates_cfi")

#------------------------------------------------------------------------------------
# Specify Global Tag
#------------------------------------------------------------------------------------
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')

#------------------------------------------------------------------------------------
# HcalTupleMaker sequence definition
#------------------------------------------------------------------------------------
process.tuple_step = cms.Sequence(
    # Make HCAL tuples: Event, run, ls number
    process.hcalTupleEvent*
    # Make HCAL tuples: digi info
    #process.hcalTupleHBHEDigis*
    #process.hcalTupleHODigis*
    #process.hcalTupleHFDigis*
    #process.hcalTupleTriggerPrimitives*
    # Make HCAL tuples: reco info
    process.hcalTupleHBHERecHits*
    # Make HCAL tuples: gen info
    process.hcalTupleGenParticles*
    #
    process.tuplePFCandidates*
    #
    process.hcalTupleTree
)

#
# in case we are using MINIAOD files
#
if options.isMINIAOD: 
    process.tuple_step = cms.Sequence(
        # Make HCAL tuples: Event, run, ls number
        process.hcalTupleEvent*
        # Make HCAL tuples: gen info
        #process.hcalTupleGenParticles*
        #
        process.tuplePackedPFCandidates*
        #
        process.hcalTupleTree
    )

#-----------------------------------------------------------------------------------
# Path and EndPath definitions
#-----------------------------------------------------------------------------------
process.preparation = cms.Path(
    process.tuple_step
)

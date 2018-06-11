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
#
# Dataset
# /RelValSinglePiPt25Eta1p7_2p7/CMSSW_9_3_2-93X_upgrade2023_realistic_v2_2023D17noPU-v1/GEN-SIM-RECO
# /SinglePiPt*Eta1p6_2p8/PhaseIITDRFall17*93X_upgrade2023_realistic_v2*/GEN-SIM-RECO
#
# pt=5 GeV sample
#options.inputFiles = '/store/mc/PhaseIITDRFall17DR/SinglePiPt5Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/00000/18D5CC99-22AA-E711-8020-90B11C08CDC7.root','/store/mc/PhaseIITDRFall17DR/SinglePiPt5Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/00000/B8AEECCD-22AA-E711-82FF-1866DA7F8E98.root'
#options.outputFile = 'results_pt5.root'
#
# pt=10 GeV sample
#options.inputFiles = '/store/mc/PhaseIITDRFall17DR/SinglePiPt10Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/150000/0844EA3B-14A9-E711-9F1D-44A842CFD633.root','/store/mc/PhaseIITDRFall17DR/SinglePiPt10Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/150000/8CF1BCF2-3AA9-E711-B0BE-A4BF0112BD74.root'
#options.outputFile = 'results_pt10.root'
#
# pt=15 GeV sample
#options.inputFiles = '/store/mc/PhaseIITDRFall17DR/SinglePiPt15Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/00000/1E9C9D26-4EA9-E711-B775-6CC2173DA2F0.root','/store/mc/PhaseIITDRFall17DR/SinglePiPt15Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/00000/5C395DE1-ACA9-E711-B391-0CC47A7EEE32.root','/store/mc/PhaseIITDRFall17DR/SinglePiPt15Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/00000/6AD5611D-81A9-E711-92A7-3417EBE7009F.root'
#options.outputFile = 'results_pt15.root'
#
# pt=25 GeV sample *relval*
options.inputFiles = '/store/relval/CMSSW_9_3_2/RelValSinglePiPt25Eta1p7_2p7/GEN-SIM-RECO/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/2637F672-C7A6-E711-B4EF-0025905A612A.root','/store/relval/CMSSW_9_3_2/RelValSinglePiPt25Eta1p7_2p7/GEN-SIM-RECO/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/3C00B396-CBA6-E711-95E1-0025905A612A.root','/store/relval/CMSSW_9_3_2/RelValSinglePiPt25Eta1p7_2p7/GEN-SIM-RECO/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/60FC86A0-CDA6-E711-960D-0025905B856E.root','/store/relval/CMSSW_9_3_2/RelValSinglePiPt25Eta1p7_2p7/GEN-SIM-RECO/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/E0E29FFF-C6A6-E711-93A0-003048FFCC16.root'
options.outputFile = 'results_pt25.root'
#
# pt=50 GeV sample
#options.inputFiles = '/store/mc/PhaseIITDRFall17DR/SinglePiPt50Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/00000/B2807A02-50AD-E711-AF78-F01FAFDB45B7.root','/store/mc/PhaseIITDRFall17DR/SinglePiPt50Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/150000/1AB4309F-BAAE-E711-87BC-A0369FC5D904.root','/store/mc/PhaseIITDRFall17DR/SinglePiPt50Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/150000/1AE8206E-E1AD-E711-98C0-FA163E191258.root','/store/mc/PhaseIITDRFall17DR/SinglePiPt50Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/150000/20F00603-FEAD-E711-83F8-0090FAA59EE4.root'
#options.outputFile = 'results_pt50.root'
#
# pt=100 GeV sample
#options.inputFiles = '/store/mc/PhaseIITDRFall17DR/SinglePiPt100Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/150000/BE82B842-31AE-E711-9EA4-0026B94DBDA2.root','/store/mc/PhaseIITDRFall17DR/SinglePiPt100Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/150000/E20ADE15-A5AE-E711-9434-0023AEEEB55F.root','/store/mc/PhaseIITDRFall17DR/SinglePiPt100Eta1p6_2p8/GEN-SIM-RECO/noPUFEVT_93X_upgrade2023_realistic_v2-v1/150000/FC41ACC6-2FAE-E711-B431-F04DA2747854.root'
#options.outputFile = 'results_pt100.root'
#
#'/store/relval/CMSSW_9_3_2/RelValSinglePiPt25Eta1p7_2p7/GEN-SIM-RECO/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/2637F672-C7A6-E711-B4EF-0025905A612A.root'
options.maxEvents = -1 # -1 means all events
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
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_Tree_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_Event_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_GenParticles_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HBHERecHits_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HFRecHits_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HORecHits_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HcalSimHits_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HBHEDigis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HFDigis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HODigis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_QIE11Digis_cfi")

process.load("Validation.HGCalValidation.hgcalHitValidation_cfi")

#------------------------------------------------------------------------------------
# Specify Global Tag
#------------------------------------------------------------------------------------
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

#----------------------------
# Paths/Sequences Definitions
#----------------------------
#process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")
#process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")

process.digiPath = cms.Path(
    process.hcalDigis
)

# Aging models -  4500/fb
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_4500_ultimate 
process = customise_aging_4500_ultimate(process)
# 
process.es_hardcode.hbUpgrade.radiationDamage.depVsNeutrons = cms.vdouble(5.543e-10,8.012e-10)
process.es_hardcode.heUpgrade.radiationDamage.depVsNeutrons = cms.vdouble(5.543e-10,8.012e-10)
process.HBDarkeningEP.drdA = cms.double(2.7383)
process.HBDarkeningEP.drdB = cms.double(0.37471)

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
    process.hcalTupleHFRecHits*
    process.hcalTupleHORecHits*
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
    #process.hgcalHitValidation*
    process.tuple_step
)

#file = open('allDump_cfg.py','w')
#file.write(str(process.dumpPython()))
#file.close()


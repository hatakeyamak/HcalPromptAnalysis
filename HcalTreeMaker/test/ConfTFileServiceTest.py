import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

from Configuration.StandardSequences.Eras import eras

#process = cms.Process("RecHitsSpectraGen",eras.Phase2)
process = cms.Process("Trees",eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )                                              

#mylist = FileUtils.loadListFromFile ('dataset.txt')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport.reportEvery = 500

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'/store/relval/CMSSW_9_0_0/RelValTTbar_13/GEN-SIM-RECO/90X_upgrade2017_realistic_v20_resub-v1/00000/1485CC3F-CF10-E711-B1C3-0CC47A7C3636.root',
'/store/relval/CMSSW_9_0_0/RelValTTbar_13/GEN-SIM-RECO/90X_upgrade2017_realistic_v20_resub-v1/00000/7C57B2B6-CE10-E711-B0B6-0025905A60D2.root',
'/store/relval/CMSSW_9_0_0/RelValTTbar_13/GEN-SIM-RECO/90X_upgrade2017_realistic_v20_resub-v1/00000/A8960510-D010-E711-B02F-0CC47A745298.root'
#'/store/relval/CMSSW_9_0_0_pre6/RelValTTbar_13/GEN-SIM-RECO/90X_upgrade2017_realistic_v15-v1/00000/9446ED7D-5201-E711-BD98-0025905A48EC.root',
#'/store/relval/CMSSW_9_0_0_pre6/RelValTTbar_13/GEN-SIM-RECO/90X_upgrade2017_realistic_v15-v1/00000/B848B53A-5201-E711-92AF-0CC47A7C351E.root',
#'/store/relval/CMSSW_9_0_0_pre6/RelValTTbar_13/GEN-SIM-RECO/90X_upgrade2017_realistic_v15-v1/00000/BACB7F3B-5201-E711-A92D-0CC47A7C347A.root',
#'/store/relval/CMSSW_9_0_0_pre6/RelValTTbar_13/GEN-SIM-RECO/90X_upgrade2017_realistic_v15-v1/00000/C6DBB840-5201-E711-B44D-0025905B8596.root'
#	'/store/relval/CMSSW_8_1_0_pre6/RelValTTbar_14TeV/GEN-SIM-RECO/80X_mcRun2_asymptotic_v14_2023LReco-v1/00000/3C1D52ED-9629-E611-8E77-0025905B8600.root',
#	'/store/relval/CMSSW_8_1_0_pre6/RelValTTbar_14TeV/GEN-SIM-RECO/80X_mcRun2_asymptotic_v14_2023LReco-v1/00000/4E5C0F20-3329-E611-B827-0CC47A4D7690.root',
#	'/store/relval/CMSSW_8_1_0_pre6/RelValTTbar_14TeV/GEN-SIM-RECO/80X_mcRun2_asymptotic_v14_2023LReco-v1/00000/60B03426-6B28-E611-A89A-0CC47A4D7686.root',
#	'/store/relval/CMSSW_8_1_0_pre6/RelValTTbar_14TeV/GEN-SIM-RECO/80X_mcRun2_asymptotic_v14_2023LReco-v1/00000/98D2ECF0-9629-E611-A7DF-0025905B857E.root',
#	'/store/relval/CMSSW_8_1_0_pre6/RelValTTbar_14TeV/GEN-SIM-RECO/80X_mcRun2_asymptotic_v14_2023LReco-v1/00000/B0727767-7028-E611-ABA9-0025905B8572.root',
#	'/store/relval/CMSSW_8_1_0_pre6/RelValTTbar_14TeV/GEN-SIM-RECO/80X_mcRun2_asymptotic_v14_2023LReco-v1/00000/B82574C8-7F28-E611-BA7B-0CC47A4D76B6.root',
#	'/store/relval/CMSSW_8_1_0_pre6/RelValTTbar_14TeV/GEN-SIM-RECO/80X_mcRun2_asymptotic_v14_2023LReco-v1/00000/D4498450-5C28-E611-8C01-0CC47A4C8E8A.root',
#	'/store/relval/CMSSW_8_1_0_pre6/RelValTTbar_14TeV/GEN-SIM-RECO/80X_mcRun2_asymptotic_v14_2023LReco-v1/00000/FC8EE7BB-7028-E611-B6DC-0CC47A4D7698.root',
#	'/store/relval/CMSSW_8_1_0_pre6/RelValTTbar_14TeV/GEN-SIM-RECO/80X_mcRun2_asymptotic_v14_2023LReco-v1/00000/FE86EE41-7628-E611-8D0E-0CC47A4D75EC.root'
    ),
    secondaryFileNames=cms.untracked.vstring(
        '/store/relval/CMSSW_9_0_0/RelValTTbar_13/GEN-SIM-DIGI-RAW/90X_upgrade2017_realistic_v20_resub-v1/00000/0AC8EF8E-C410-E711-ACDA-0025905B8568.root',
'/store/relval/CMSSW_9_0_0/RelValTTbar_13/GEN-SIM-DIGI-RAW/90X_upgrade2017_realistic_v20_resub-v1/00000/505B7990-C410-E711-8420-0025905A60DE.root',
'/store/relval/CMSSW_9_0_0/RelValTTbar_13/GEN-SIM-DIGI-RAW/90X_upgrade2017_realistic_v20_resub-v1/00000/66A15436-C410-E711-9BA2-0CC47A4C8E3C.root',
'/store/relval/CMSSW_9_0_0/RelValTTbar_13/GEN-SIM-DIGI-RAW/90X_upgrade2017_realistic_v20_resub-v1/00000/9856BA9D-C410-E711-9896-0025905A6066.root',
'/store/relval/CMSSW_9_0_0/RelValTTbar_13/GEN-SIM-DIGI-RAW/90X_upgrade2017_realistic_v20_resub-v1/00000/ECEE2DF2-C410-E711-A86A-0CC47A78A3EE.root',
'/store/relval/CMSSW_9_0_0/RelValTTbar_13/GEN-SIM-DIGI-RAW/90X_upgrade2017_realistic_v20_resub-v1/00000/ECF5E993-C410-E711-AC63-0025905B85D8.root',
'/store/relval/CMSSW_9_0_0/RelValTTbar_13/GEN-SIM-DIGI-RAW/90X_upgrade2017_realistic_v20_resub-v1/00000/F84F5F36-C410-E711-85A5-0CC47A4C8E3C.root' 
    )  
)

########## Good run list ##########
#import PhysicsTools.PythonAnalysis.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_136033-149442_7TeV_Apr21ReReco_Collisions10_JSON_v2.txt').getVLuminosityBlockRange()
##########


process.TFileService = cms.Service("TFileService", fileName = cms.string("TFileServiceTest1.root") )

#process.load('RecoLocalCalo/HcalRecAlgos/hcalRecAlgoESProd_cfi')

#process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
#process.HBHENoiseFilterResultProducer.minZeros = cms.int32(9999)

#process.ApplyBaselineHBHENoiseFilter = cms.EDFilter(
#	'BooleanFlagFilter',
#	inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
#	reverseDecision = cms.bool(False)
#)

#process.goodVertices = cms.EDFilter(
#	'VertexSelector',
#	filter = cms.bool(True),
#	src = cms.InputTag("offlinePrimaryVertices"),
#	cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
#)

process.load('HcalPromptAnalysis/HcalTreeMaker/SimHitTree_cfi')

process.load('HcalPromptAnalysis/HcalTreeMaker/RecHitTree_cfi')
#process.hcalRecHitTree.rootOutputFile            = cms.string('TFileServiceTest2.root')
#process.hcalRecHitTree.HBHERecHitCollectionLabel = cms.untracked.InputTag("hbhereco")
#process.hcalRecHitTree.HFRecHitCollectionLabel   = cms.untracked.InputTag("hfreco")
process.hcalRecHitTree.HBHERecHitCollectionLabel = cms.untracked.InputTag("hbheprereco")
process.hcalRecHitTree.HFRecHitCollectionLabel   = cms.untracked.InputTag("hfreco")

#process.hcalRecHitTreePre = process.hcalRecHitTree.clone()
#process.hcalRecHitTreePre.treeName = cms.string('HcalRecHitPre')
#process.hcalRecHitTreePre.HBHERecHitCollectionLabel = cms.untracked.InputTag("hbheprereco")
#process.hcalRecHitTreePre.HFRecHitCollectionLabel   = cms.untracked.InputTag("hfreco")

process.plots = cms.Path(process.hcalSimHitTree*process.hcalRecHitTree)


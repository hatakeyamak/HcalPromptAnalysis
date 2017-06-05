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
mylist2 = FileUtils.loadListFromFile ('step2.txt')
mylist1 = FileUtils.loadListFromFile ('step1.txt')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport.reportEvery = 500

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(*mylist2),
      secondaryFileNames=cms.untracked.vstring(*mylist1)  
)

########## Good run list ##########
#import PhysicsTools.PythonAnalysis.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_136033-149442_7TeV_Apr21ReReco_Collisions10_JSON_v2.txt').getVLuminosityBlockRange()
##########

process.TFileService = cms.Service("TFileService", fileName = cms.string("TFileServiceOutputTree.root") )

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

#---------------------------------------------------------------------------
#Chris' Changes
process.load('HcalPromptAnalysis/HcalTreeMaker/HcalDigiTree_cfi')
#---------------------------------------------------------------------------

#process.hcalRecHitTree.rootOutputFile            = cms.string('TFileServiceTest2.root')
#process.hcalRecHitTree.HBHERecHitCollectionLabel = cms.untracked.InputTag("hbhereco")
#process.hcalRecHitTree.HFRecHitCollectionLabel   = cms.untracked.InputTag("hfreco")
process.hcalRecHitTree.HBHERecHitCollectionLabel = cms.untracked.InputTag("hbheprereco")
process.hcalRecHitTree.HFRecHitCollectionLabel   = cms.untracked.InputTag("hfreco")

#process.hcalRecHitTreePre = process.hcalRecHitTree.clone()
#process.hcalRecHitTreePre.treeName = cms.string('HcalRecHitPre')
#process.hcalRecHitTreePre.HBHERecHitCollectionLabel = cms.untracked.InputTag("hbheprereco")
#process.hcalRecHitTreePre.HFRecHitCollectionLabel   = cms.untracked.InputTag("hfreco")

#---------------------------------------------------------------------------
#Chris' Changes
#process.plots = cms.Path(process.hcalSimHitTree*process.hcalDigiTree)
process.plots = cms.Path(process.hcalDigiTree)
#---------------------------------------------------------------------------

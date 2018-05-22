import FWCore.ParameterSet.Config as cms

hcalTupleHcalSimHits = cms.EDProducer("HcalTupleMaker_HcalSimHits",
  source = cms.untracked.VInputTag(
        #cms.untracked.InputTag("g4SimHits","HGCHitsEE"),
        #cms.untracked.InputTag("g4SimHits","HGCHitsHEfront"),
        cms.untracked.InputTag("g4SimHits","HcalHits")
        ),
  geometrySource = cms.untracked.vstring(
        #'HGCalEESensitive',
        #'HGCalHESiliconSensitive',
        'HCal'
  ),
  Prefix = cms.untracked.string  ("HcalSimHits"),
  Suffix = cms.untracked.string  ("")
)



import FWCore.ParameterSet.Config as cms

hcalTupleHGCRecHits = cms.EDProducer("HcalTupleMaker_HGCRecHits",
  source = cms.untracked.VInputTag(
        cms.untracked.InputTag("HGCalRecHit","HGCEERecHits"),
        cms.untracked.InputTag("HGCalRecHit","HGCHEFRecHits"),
        cms.untracked.InputTag("HGCalRecHit","HGCHEBRecHits")
        ),
  geometrySource = cms.untracked.vstring(
        'HGCalEESensitive',
        'HGCalHESiliconSensitive',
        'HCal'
  ),
  Prefix = cms.untracked.string  ("HGCRecHit"),
  Suffix = cms.untracked.string  ("")
)


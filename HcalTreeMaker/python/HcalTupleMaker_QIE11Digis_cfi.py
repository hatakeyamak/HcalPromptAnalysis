import FWCore.ParameterSet.Config as cms

hcalTupleQIE11Digis = cms.EDProducer("HcalTupleMaker_QIE11Digis",
  tagQIE11  = cms.untracked.InputTag("hcalDigis"),
  taguMNio  = cms.untracked.InputTag("hcalDigis"),
  Prefix  = cms.untracked.string("HBHEQIE11Digi"),
  Suffix  = cms.untracked.string(""),
  StoreLaser = cms.untracked.bool(False)
)

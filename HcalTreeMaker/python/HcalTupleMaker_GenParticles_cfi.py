import FWCore.ParameterSet.Config as cms

hcalTupleGenParticles = cms.EDProducer("HcalTupleMaker_GenParticles",
  source    = cms.untracked.InputTag("genParticles"),
  HepMCProduct = cms.untracked.bool(False),
  Prefix    = cms.untracked.string  ("GenPar"),
  Suffix    = cms.untracked.string  ("")
)


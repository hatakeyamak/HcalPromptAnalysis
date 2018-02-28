import FWCore.ParameterSet.Config as cms

hcalTupleGenParticles = cms.EDProducer("HcalTupleMaker_GenParticles",
  source    = cms.untracked.InputTag("genParticles"),
  Prefix    = cms.untracked.string  ("GenPar"),
  Suffix    = cms.untracked.string  ("")
)


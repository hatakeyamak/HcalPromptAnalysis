import FWCore.ParameterSet.Config as cms

hcalSimHitTree = cms.EDFilter("SimHitTree",
    SimHitCollectionLabel = cms.untracked.InputTag("g4SimHits","HcalHits"),
    SubDetector   = cms.untracked.int32(0),
    TestNumbering = cms.bool(True)
    #SubDetector   = cms.untracked.int32(4),
    #TestNumbering = cms.bool(False)
)

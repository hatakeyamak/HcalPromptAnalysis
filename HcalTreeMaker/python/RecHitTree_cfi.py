import FWCore.ParameterSet.Config as cms

hcalRecHitTree = cms.EDAnalyzer("RecHitTree",
    rootOutputFile            = cms.string('HcalRecHit.root'),
    treeName                  = cms.string('HcalRecHit'),
    TestNumbering             = cms.bool(False),
    HBHERecHitCollectionLabel = cms.untracked.InputTag("hbheUpgradeReco"),
    HFRecHitCollectionLabel   = cms.untracked.InputTag("hfUpgradeReco"),
    HORecHitCollectionLabel   = cms.untracked.InputTag("horeco")                             
)

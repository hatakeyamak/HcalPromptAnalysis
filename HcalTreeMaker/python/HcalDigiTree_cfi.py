import FWCore.ParameterSet.Config as cms

hcalDigiTree = cms.EDAnalyzer("RecHitTree",
    rootOutputFile            = cms.string('HcalDigi.root'),
    treeName                  = cms.string('HcalRecHit'),
    TestNumbering             = cms.bool(True),
    HBHERecHitCollectionLabel = cms.untracked.InputTag("hbheUpgradeReco"),
    HFRecHitCollectionLabel   = cms.untracked.InputTag("hfUpgradeReco")
)

#hcalRecHitTree = cms.EDAnalyzer("RecHitTree",
#    rootOutputFile            = cms.string('HcalRecHitSpectra.root'),
#    treeName                  = cms.string('HcalRecHit'),
#    TestNumbering             = cms.bool(True),
#    HBHERecHitCollectionLabel = cms.untracked.InputTag("hbheUpgradeReco"),
#    HFRecHitCollectionLabel   = cms.untracked.InputTag("hfUpgradeReco")
#)

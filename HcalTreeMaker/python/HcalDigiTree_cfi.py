import FWCore.ParameterSet.Config as cms

hcalDigiTree = cms.EDAnalyzer("HcalDigiTree",
    rootOutputFile            = cms.string('HcalDigi.root'),
    treeName                  = cms.string('HcalDigi'),
    TestNumbering             = cms.bool(True),
    digiTag                   = cms.InputTag("hcalDigis"),
    QIE10digiTag              = cms.InputTag("hcalDigis"),
    QIE11digiTag              = cms.InputTag("hcalDigis")
)

#hcalRecHitTree = cms.EDAnalyzer("RecHitTree",
#    rootOutputFile            = cms.string('HcalRecHitSpectra.root'),
#    treeName                  = cms.string('HcalRecHit'),
#    TestNumbering             = cms.bool(True),
#    HBHERecHitCollectionLabel = cms.untracked.InputTag("hbheUpgradeReco"),
#    HFRecHitCollectionLabel   = cms.untracked.InputTag("hfUpgradeReco")
#)

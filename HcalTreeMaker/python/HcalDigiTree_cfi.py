import FWCore.ParameterSet.Config as cms

hcalDigiTree = cms.EDAnalyzer("HcalDigiTree",
    rootOutputFile            = cms.string('HcalDigiTree.root'),
    treeName                  = cms.string('HcalDigi'),
    TestNumbering             = cms.bool(True),

    digiTag                   = cms.InputTag("simHcalDigis"),
    QIE10digiTag              = cms.InputTag("simHcalDigis"),
    QIE11digiTag              = cms.InputTag("simHcalDigis"),

    #digiTag                   = cms.InputTag("hcalDigis"),
    #QIE10digiTag              = cms.InputTag("hcalDigis"),
    #QIE11digiTag              = cms.InputTag("simHcalUnsuppressedDigis"),

    HFRecHitCollectionLabel   = cms.untracked.InputTag("hfUpgradeReco")
)

#hcalRecHitTree = cms.EDAnalyzer("RecHitTree",
#    rootOutputFile            = cms.string('HcalRecHitSpectra.root'),
#    treeName                  = cms.string('HcalRecHit'),
#    TestNumbering             = cms.bool(True),
#    HBHERecHitCollectionLabel = cms.untracked.InputTag("hbheUpgradeReco"),
#    HFRecHitCollectionLabel   = cms.untracked.InputTag("hfUpgradeReco")
#)

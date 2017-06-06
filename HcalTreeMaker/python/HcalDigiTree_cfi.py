import FWCore.ParameterSet.Config as cms

hcalDigiTree = cms.EDAnalyzer("HcalDigiTree",
    rootOutputFile            = cms.string('HcalDigiTree.root'),
    treeName                  = cms.string('HcalDigiTree'),
    TestNumbering             = cms.bool(True),
    #digiTag                   = cms.InputTag("simHcalDigis"),
    #QIE10digiTag              = cms.InputTag("simHcalDigis","HFQIE10DigiCollection"),
    #QIE11digiTag              = cms.InputTag("simHcalDigis","HBHEQIE11DigiCollection"),
    digiTag                   = cms.InputTag("hcalDigis"),
    QIE10digiTag              = cms.InputTag("hcalDigis","HFQIE10DigiCollection"),
    QIE11digiTag              = cms.InputTag("hcalDigis","HBHEQIE11DigiCollection"),
    HFRecHitCollectionLabel   = cms.untracked.InputTag("hfUpgradeReco")
)

#hcalRecHitTree = cms.EDAnalyzer("RecHitTree",
#    rootOutputFile            = cms.string('HcalRecHitSpectra.root'),
#    treeName                  = cms.string('HcalRecHit'),
#    TestNumbering             = cms.bool(True),
#    HBHERecHitCollectionLabel = cms.untracked.InputTag("hbheUpgradeReco"),
#    HFRecHitCollectionLabel   = cms.untracked.InputTag("hfUpgradeReco")
#)

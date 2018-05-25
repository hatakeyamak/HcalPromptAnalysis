import sys
jobid=sys.argv[2]
nevt=sys.argv[3]
firstEvt=1+int(nevt)*int(jobid)
Nevt=int(nevt)
print jobid,nevt,firstEvt

import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('TEST',eras.Run2_2018)

### RANDOM setting (change last digit(s) to make runs different !)
process.load("Configuration.StandardSequences.SimulationRandomNumberGeneratorSeeds_cff")

process.load("Configuration.StandardSequences.Simulation_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load('Configuration/StandardSequences/DigiToRaw_cff')
process.load('Configuration/StandardSequences/RawToDigi_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# from Configuration.AlCa.autoCond import autoCond
# process.GlobalTag.globaltag = autoCond['run2_mc']
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')

process.load("IOMC.EventVertexGenerators.VtxSmearedGauss_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.g4SimHits.UseMagneticField = False

process.load("DQMServices.Core.DQMStore_cfi")
process.load("DQMServices.Components.MEtoEDMConverter_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(Nevt) 
)
# Input source
process.source = cms.Source("PoolSource",
    firstEvent = cms.untracked.uint32(firstEvt), 
    fileNames = cms.untracked.vstring('root://kodiak-se.baylor.edu//store/group/dpg_hcal/common/mc_pi50_eta05.root') 
) 

process.FEVT = cms.OutputModule("PoolOutputModule",
     outputCommands = cms.untracked.vstring('drop *', 'keep *_MEtoEDMConverter_*_*'),
     #outputCommands = cms.untracked.vstring('drop *_MEtoEDMXXX_*_*', 'keep *'),
     splitLevel = cms.untracked.int32(0),
     fileName = cms.untracked.string("output.root")
)

process.VtxSmeared.SigmaX = 0.00001
process.VtxSmeared.SigmaY = 0.00001
process.VtxSmeared.SigmaZ = 0.00001

process.load("Validation.HcalHits.HcalSimHitsValidation_cfi")
process.HcalSimHitsAnalyser.outputFile = cms.untracked.string('HcalSimHitsValidation.root')

process.load("Validation.HcalDigis.HcalDigisParam_cfi")
process.hcaldigisAnalyzer.outputFile = cms.untracked.string('HcalDigisValidationRelVal.root')

process.load("Validation.HcalRecHits.HcalRecHitParam_cfi")

process.load("Validation.CaloTowers.CaloTowersParam_cfi")
process.calotowersAnalyzer.outputFile = cms.untracked.string('CaloTowersValidationRelVal.root')

#--- replace hbhereco with hbheprereco
delattr(process,"hbhereco")
process.hbhereco = process.hbheprereco.clone()
process.hcalLocalRecoSequence = cms.Sequence(process.hbhereco+process.hfprereco+process.hfreco+process.horeco)

#--- post-LS1 customization 
process.mix.digitizers.hcal.minFCToDelay=cms.double(5.) # new TS model
process.mix.digitizers.hcal.ho.photoelectronsToAnalog = cms.vdouble([4.0]*16)
process.mix.digitizers.hcal.ho.siPMCode = cms.int32(1)
process.mix.digitizers.hcal.ho.pixels = cms.int32(2500)
process.mix.digitizers.hcal.ho.doSiPMSmearing = cms.bool(False)
process.mix.digitizers.hcal.hf1.samplingFactor = cms.double(0.67)
process.mix.digitizers.hcal.hf2.samplingFactor = cms.double(0.67)
process.g4SimHits.HFShowerLibrary.FileName = 'SimG4CMS/Calo/data/HFShowerLibrary_npmt_noatt_eta4_16en_v4.root'

#------------------------------------------------------------------------------------
# Set up our analyzer
#------------------------------------------------------------------------------------
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_Tree_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_Event_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_GenParticles_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HBHERecHits_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HcalSimHits_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HBHEDigis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HFDigis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HODigis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_QIE11Digis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HcalTriggerPrimitives_cfi")

#------------------------------------------------------------------------------------
# HcalTupleMaker sequence definition
#------------------------------------------------------------------------------------
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')
options.outputFile = 'pi50_trees_MCfull.root'
process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string(options.outputFile)
)
process.hcalTupleGenParticles.source = cms.untracked.InputTag("generatorSmeared")
process.hcalTupleGenParticles.HepMCProduct = cms.untracked.bool(True)

process.tuple_step = cms.Sequence(
    # Make HCAL tuples: Event, run, ls number
    process.hcalTupleEvent*
    # Make HCAL tuples: digi info
    process.hcalTupleHBHEDigis*
    process.hcalTupleHODigis*
    process.hcalTupleHFDigis*
    process.hcalTupleQIE11Digis*
    #process.hcalTupleTriggerPrimitives*
    # Make HCAL tuples: reco info
    process.hcalTupleHBHERecHits*
    # Make HCAL tuples: gen & simhit info
    process.hcalTupleHcalSimHits*
    process.hcalTupleGenParticles*
    #
    process.hcalTupleTree
)

#---------- PATH
# -- NB: for vertex smearing the Label should be: "unsmeared" 
# for GEN produced since 760pre6, for older GEN - just "": 
process.VtxSmeared.src = cms.InputTag("generator", "") 
process.generatorSmeared = cms.EDProducer("GeneratorSmearedProducer")
process.g4SimHits.Generator.HepMCProductLabel = 'VtxSmeared' 

process.p = cms.Path(
 process.VtxSmeared *
 process.generatorSmeared *
 process.g4SimHits *
 process.mix *
 process.ecalDigiSequence * 
 process.hcalDigiSequence *
 process.addPileupInfo *
 process.bunchSpacingProducer *
 process.ecalPacker *
 process.esDigiToRaw *
 process.hcalRawData *
 process.rawDataCollector *
 process.ecalDigis *
 process.ecalPreshowerDigis *
 process.hcalDigis *
 process.calolocalreco *
 process.caloTowersRec *
 process.hcalnoise *
 process.HcalSimHitsAnalyser *
 process.hcaldigisAnalyzer *
 process.calotowersAnalyzer *
 process.hcalRecoAnalyzer *
 process.tuple_step *
 process.MEtoEDMConverter
)

process.outpath = cms.EndPath(process.FEVT)


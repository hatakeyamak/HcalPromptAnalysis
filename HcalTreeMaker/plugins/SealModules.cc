#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_Tree.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_Event.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_HcalDigis.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_QIE10Digis.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_QIE11Digis.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_HcalRecHits.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_HGCRecHits.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_HcalTriggerPrimitives.h"

DEFINE_FWK_MODULE(HcalTupleMaker_Tree);
DEFINE_FWK_MODULE(HcalTupleMaker_Event);
DEFINE_FWK_MODULE(HcalTupleMaker_HBHERecHits);
DEFINE_FWK_MODULE(HcalTupleMaker_HORecHits);
DEFINE_FWK_MODULE(HcalTupleMaker_HFRecHits);
DEFINE_FWK_MODULE(HcalTupleMaker_HGCRecHits);
DEFINE_FWK_MODULE(HcalTupleMaker_HBHEDigis);
DEFINE_FWK_MODULE(HcalTupleMaker_HODigis);
DEFINE_FWK_MODULE(HcalTupleMaker_HFDigis);
DEFINE_FWK_MODULE(HcalTupleMaker_QIE10Digis);
DEFINE_FWK_MODULE(HcalTupleMaker_QIE11Digis);
DEFINE_FWK_MODULE(HcalTupleMaker_HcalTriggerPrimitives);

// -*- C++ -*-
//
// Class:      RecHitTree
// 
/**\class RecHitTree RecHitTree.cc  HcalPromptAnalysis/HcalTreeMaker/plugins/RecHitTree.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Kenneth Call
//         Created:  Wed, 15 Jul 2015 16:05:14 GMT
//
//

//Include files to use TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//All include files from HcalRecHitsValidation
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalSourcePositionData.h"
#include <DataFormats/EcalDetId/interface/EBDetId.h>
#include <DataFormats/EcalDetId/interface/EEDetId.h>
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "Geometry/HcalCommonData/interface/HcalHitRelabeller.h"
#include "DataFormats/HcalDetId/interface/HcalTestNumbering.h"
#include <vector>
#include <utility>
#include <ostream>
#include <string>
#include <algorithm>
#include <cmath>
#include "DataFormats/DetId/interface/DetId.h"
// channel status
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/HcalObjects/interface/HcalChannelQuality.h"
#include "CondFormats/DataRecord/interface/HcalChannelQualityRcd.h"
// severity level assignment for HCAL
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputer.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputerRcd.h"
// severity level assignment for ECAL
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
//All include files from HcalRecHitsValidation
// system include files
#include <memory>
// user include files
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalSourcePositionData.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "CondFormats/HcalObjects/interface/HcalChannelQuality.h"
#include "CondFormats/DataRecord/interface/HcalChannelQualityRcd.h"
// severity level assignment for HCAL
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputer.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputerRcd.h"
#include <math.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <vector>
#include <utility>
#include <ostream>
#include <string>
#include <algorithm>
//CaloTowers
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

//
// class declaration
//

class RecHitTree : public edm::EDAnalyzer {
   public:
      explicit RecHitTree(const edm::ParameterSet&);
      ~RecHitTree();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  // ----------member data ---------------------------

  std::string outputfile_;
  std::string treename_;

  edm::EDGetTokenT<HBHERecHitCollection> tok_hbhe_;
  edm::EDGetTokenT<HORecHitCollection> tok_ho_;
  edm::EDGetTokenT<HFRecHitCollection> tok_hf_;
  
  //TFile *tf1;
  TTree *tt1;

  // run:lumi:event
  int run;
  int lumi;
  int event;

  std::vector<float> recHitHB_En;
  std::vector<float> recHitHB_EnRAW;
  std::vector<float> recHitHB_EnM3;
  std::vector<int>   recHitHB_ieta;
  std::vector<int>   recHitHB_iphi;
  std::vector<float> recHitHB_time;
  std::vector<int>   recHitHB_depth;
  std::vector<float> chi2HB;

  std::vector<float> recHitHE_En;
  std::vector<float> recHitHE_EnRAW;
  std::vector<float> recHitHE_EnM3;
  std::vector<int>   recHitHE_ieta;
  std::vector<int>   recHitHE_iphi;
  std::vector<float> recHitHE_time;
  std::vector<int>   recHitHE_depth;
  std::vector<float> chi2HE;

  std::vector<float> recHitHF_En;
  std::vector<int>   recHitHF_ieta;
  std::vector<int>   recHitHF_iphi;
  std::vector<float> recHitHF_time;
  std::vector<int>   recHitHF_depth;

  std::vector<float> recHitHO_En;
  std::vector<int>   recHitHO_ieta;
  std::vector<int>   recHitHO_iphi;
  std::vector<float> recHitHO_time;

  bool testNumbering_;

};

//
// constructors and destructor
//
RecHitTree::RecHitTree(const edm::ParameterSet& iConfig)

{
  //outputfile_ = iConfig.getParameter<std::string>("rootOutputFile");
  treename_ = iConfig.getParameter<std::string>("treeName");
  testNumbering_ = iConfig.getParameter<bool>("TestNumbering");

  //Collections
  tok_hbhe_ = consumes<HBHERecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("HBHERecHitCollectionLabel"));
  tok_hf_  = consumes<HFRecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("HFRecHitCollectionLabel"));
  tok_ho_ = consumes<HORecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("HORecHitCollectionLabel"));

  //tf1 = new TFile(outputfile_.c_str(), "RECREATE");

  //Now we use the modification so that we can use the TFileService
  edm::Service<TFileService> fs;

  tt1 = fs->make<TTree>(treename_.c_str(),treename_.c_str());

  //branches
  tt1->Branch("run", &run, "run/I");
  tt1->Branch("lumi", &lumi, "lumi/I");
  tt1->Branch("event", &event, "event/I");

  tt1->Branch("recHitHB_En","std::vector<float>", &recHitHB_En, 32000, 0);
  tt1->Branch("recHitHB_EnRAW","std::vector<float>", &recHitHB_EnRAW, 32000, 0);
  tt1->Branch("recHitHB_EnM3","std::vector<float>", &recHitHB_EnM3, 32000, 0);
  tt1->Branch("recHitHB_ieta","std::vector<int>", &recHitHB_ieta, 32000, 0);
  tt1->Branch("recHitHB_iphi","std::vector<int>", &recHitHB_iphi, 32000, 0);
  tt1->Branch("recHitHB_depth","std::vector<int>", &recHitHB_depth, 32000, 0);
  tt1->Branch("recHitHB_time","std::vector<float>", &recHitHB_time, 32000, 0);
  tt1->Branch("chi2HB","std::vector<float>", &chi2HB, 32000, 0);

  tt1->Branch("recHitHE_En","std::vector<float>", &recHitHE_En, 32000, 0);
  tt1->Branch("recHitHE_EnRAW","std::vector<float>", &recHitHE_EnRAW, 32000, 0);
  tt1->Branch("recHitHE_EnM3","std::vector<float>", &recHitHE_EnM3, 32000, 0);
  tt1->Branch("recHitHE_ieta","std::vector<int>", &recHitHE_ieta, 32000, 0);
  tt1->Branch("recHitHE_iphi","std::vector<int>", &recHitHE_iphi, 32000, 0);
  tt1->Branch("recHitHE_depth","std::vector<int>", &recHitHE_depth, 32000, 0);
  tt1->Branch("recHitHE_time","std::vector<float>", &recHitHE_time, 32000, 0);
  tt1->Branch("chi2HE","std::vector<float>", &chi2HE, 32000, 0);

  tt1->Branch("recHitHF_En","std::vector<float>", &recHitHF_En, 32000, 0);
  tt1->Branch("recHitHF_ieta","std::vector<int>", &recHitHF_ieta, 32000, 0);
  tt1->Branch("recHitHF_iphi","std::vector<int>", &recHitHF_iphi, 32000, 0);
  tt1->Branch("recHitHF_depth","std::vector<int>", &recHitHF_depth, 32000, 0);
  tt1->Branch("recHitHF_time","std::vector<float>", &recHitHF_time, 32000, 0);

  tt1->Branch("recHitHO_En","std::vector<float>", &recHitHO_En, 32000, 0);
  tt1->Branch("recHitHO_ieta","std::vector<int>", &recHitHO_ieta, 32000, 0);
  tt1->Branch("recHitHO_iphi","std::vector<int>", &recHitHO_iphi, 32000, 0);
  tt1->Branch("recHitHO_time","std::vector<float>", &recHitHO_time, 32000, 0);

}

RecHitTree::~RecHitTree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  //tf1->cd();
  tt1->Write();
  //tf1->Write();
  //tf1->Close();

}

//
// member functions
//

// ------------ method called for each event  ------------
void
RecHitTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;

   edm::ESHandle<HcalDDDRecConstants> pHRNDC;
   iSetup.get<HcalRecNumberingRecord>().get( pHRNDC );
   const HcalDDDRecConstants* hcons = &(*pHRNDC);

   //Clear out the RecHit vectors
   recHitHB_En.clear();
   recHitHB_EnRAW.clear();
   recHitHB_EnM3.clear();
   recHitHB_ieta.clear();
   recHitHB_iphi.clear();
   recHitHB_depth.clear();
   recHitHB_time.clear();
   chi2HB.clear();

   recHitHE_En.clear();
   recHitHE_EnRAW.clear();
   recHitHE_EnM3.clear();
   recHitHE_ieta.clear();
   recHitHE_iphi.clear();
   recHitHE_depth.clear();
   recHitHE_time.clear();
   chi2HE.clear();

   recHitHF_En.clear();
   recHitHF_ieta.clear();
   recHitHF_iphi.clear();
   recHitHF_depth.clear();
   recHitHF_time.clear();

   recHitHO_En.clear();
   recHitHO_ieta.clear();
   recHitHO_iphi.clear();
   recHitHO_time.clear();

   //run:lumi:event
   run = iEvent.id().run();
   lumi = iEvent.id().luminosityBlock();
   event = iEvent.id().event();

   //HBHE RecHits
   edm::Handle<HBHERecHitCollection> hbhecoll;
   iEvent.getByToken(tok_hbhe_, hbhecoll);

   int depth = 0;
   for(HBHERecHitCollection::const_iterator j=hbhecoll->begin(); j != hbhecoll->end(); j++){
     HcalDetId cell;
     if (testNumbering_) cell = HcalHitRelabeller::relabel(j->id(),hcons);
     else cell = HcalDetId(j->id());
     depth = cell.depth();

     if(cell.subdet() == HcalBarrel){
       recHitHB_En.push_back(j->energy());
       recHitHB_EnRAW.push_back(j->eraw());
       recHitHB_EnM3.push_back(j->eaux());
       recHitHB_ieta.push_back(cell.ieta());
       recHitHB_iphi.push_back(cell.iphi());
       recHitHB_depth.push_back(depth);
       recHitHB_time.push_back(j->time());
       chi2HB.push_back(j->chi2());
     }//HB

     if(cell.subdet() == HcalEndcap){
       recHitHE_En.push_back(j->energy());
       recHitHE_EnRAW.push_back(j->eraw());
       recHitHE_EnM3.push_back(j->eaux());
       recHitHE_ieta.push_back(cell.ieta());
       recHitHE_iphi.push_back(cell.iphi());
       recHitHE_depth.push_back(depth);
       recHitHE_time.push_back(j->time());
       chi2HE.push_back(j->chi2());
     }//HE            
   } //HBHE

   //HF RecHits
   edm::Handle<HFRecHitCollection> hfcoll;
   iEvent.getByToken(tok_hf_, hfcoll);
      for( HFRecHitCollection::const_iterator j = hfcoll->begin(); j != hfcoll->end(); j++){
     HcalDetId cell(j->id());
     depth = cell.depth();
     recHitHF_En.push_back(j->energy());
     recHitHF_ieta.push_back(cell.ieta());
     recHitHF_iphi.push_back(cell.iphi());
     recHitHF_depth.push_back(depth);
     recHitHF_time.push_back(j->time());
   } //HF

   //HO RecHits
   edm::Handle<HORecHitCollection> hocoll;
   iEvent.getByToken(tok_ho_, hocoll);
   for( HORecHitCollection::const_iterator j = hocoll->begin(); j != hocoll->end(); j++){
     HcalDetId cell(j->id());
     //severityLevel = hcalSevLvl( (CaloRecHit*) &*j );
     //if(severityLevel > 8) continue;
     //nrechits++;
     recHitHO_En.push_back(j->energy());
     recHitHO_ieta.push_back(cell.ieta());
     recHitHO_iphi.push_back(cell.iphi());
     recHitHO_time.push_back(j->time());
   } //HO

   //Fill Tree
   tt1->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
RecHitTree::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RecHitTree::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecHitTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecHitTree);

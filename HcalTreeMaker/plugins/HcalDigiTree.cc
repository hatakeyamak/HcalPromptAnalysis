// -*- C++ -*-
//
// Package:    HcalPromptAnalysis/HcalTreeMaker/
// Class:      HcalDigiTree
// 
/**\class HcalDigiTree HcalDigiTree.cc HcalPromptAnalysis/HcalTreeMaker/plugins/HcalDigiTree.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Authors:  Kenneth Call and Christopher Madrid
//         Created:  Sun, 21 May 2017 16:05:14 GMT
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
#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"
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


//--KH
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "DQMServices/Core/interface/DQMStore.h"
//#include "DQMServices/Core/interface/MonitorElement.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "Geometry/Records/interface/HcalGeometryRecord.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"

#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CondFormats/HcalObjects/interface/HcalQIEShape.h"

#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"


#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "Geometry/CaloTopology/interface/HcalTopology.h"


/*TP Code*/
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
/*~TP Code*/
//--KH


//
// class declaration
//

class HcalDigiTree : public edm::EDAnalyzer {
   public:
      explicit HcalDigiTree(const edm::ParameterSet&);
      ~HcalDigiTree();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

  std::string outputfile_;
  std::string treename_;

  //edm::EDGetTokenT<HBHERecHitCollection> tok_hbhe_;
  //edm::EDGetTokenT<HFRecHitCollection> tok_hf_;

  //  edm::EDGetTokenT<reco::VertexCollection> tok_vertexc_;
  //  edm::EDGetTokenT<CaloTowerCollection> tok_calo_;

  TFile *tf1;
  TTree *tt1;

  //int nrechits;
  //int nvertx;

  // run:lumi:event
  int run;
  int lumi;
  int event;

  //std::vector<float> recHitEn_HB;
  //std::vector<float> recHitEn_HE;
  //std::vector<float> recHitEn_HF;
  //std::vector<float> recHitEn_HO;  

  std::vector<float> recHitHB_En;
  std::vector<float> recHitHB_EnRAW;
  std::vector<int>   recHitHB_ieta;
  std::vector<int>   recHitHB_iphi;
  std::vector<float> recHitHB_time;
  std::vector<int>   recHitHB_depth;

  std::vector<float> recHitHE_En;
  std::vector<float> recHitHE_EnRAW;
  std::vector<int>   recHitHE_ieta;
  std::vector<int>   recHitHE_iphi;
  std::vector<float> recHitHE_time;
  std::vector<int>   recHitHE_depth;

  std::vector<float> recHitHF_En;
  std::vector<float> recHitHF_EnRAW;
  std::vector<int>   recHitHF_ieta;
  std::vector<int>   recHitHF_iphi;
  std::vector<float> recHitHF_time;
  std::vector<int>   recHitHF_depth;

  bool testNumbering_;

  edm::InputTag inputTag_;
  edm::InputTag QIE10inputTag_;
  edm::InputTag QIE11inputTag_;

  edm::EDGetTokenT< HBHEDigiCollection > tok_hbhe_; 
  edm::EDGetTokenT< HODigiCollection > tok_ho_;
  edm::EDGetTokenT< HFDigiCollection > tok_hf_;

  edm::EDGetTokenT< QIE10DigiCollection > tok_qie10_hf_; 
  edm::EDGetTokenT< QIE11DigiCollection > tok_qie11_hbhe_; 

   // for checking the status of ECAL and HCAL channels stored in the DB 
  //const HcalChannelQuality* theHcalChStatus;
  // calculator of severety level for HCAL
  //const HcalSeverityLevelComputer* theHcalSevLvlComputer;
  //int hcalSevLvl(const CaloRecHit* hit);

  /*
  template<class Digi> void reco(const edm::Event& iEvent, const edm::EventSetup& iSetup, const edm::EDGetTokenT<edm::SortedCollection<Digi> > &tok);
  template<class dataFrameType> void reco(const edm::Event& iEvent, const edm::EventSetup& iSetup, const edm::EDGetTokenT<HcalDataFrameContainer<dataFrameType> > &tok);
  */

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HcalDigiTree::HcalDigiTree(const edm::ParameterSet& iConfig)

{
  outputfile_ = iConfig.getParameter<std::string>("rootOutputFile");
  treename_ = iConfig.getParameter<std::string>("treeName");
  testNumbering_ = iConfig.getParameter<bool>("TestNumbering");

  //subdet_ = iConfig.getUntrackedParameter<std::string > ("subdetector", "all");
  //outputFile_ = iConfig.getUntrackedParameter<std::string > ("outputFile", "");
  //    inputLabel_ = iConfig.getParameter<std::string > ("digiLabel");
  inputTag_   = iConfig.getParameter<edm::InputTag > ("digiTag");
  QIE10inputTag_   = iConfig.getParameter<edm::InputTag > ("QIE10digiTag");
  QIE11inputTag_   = iConfig.getParameter<edm::InputTag > ("QIE11digiTag");
  //emulTPsTag_ = iConfig.getParameter<edm::InputTag > ("emulTPs");
  //dataTPsTag_ = iConfig.getParameter<edm::InputTag > ("dataTPs");
  //mc_ = iConfig.getUntrackedParameter<std::string > ("mc", "no");
  //mode_ = iConfig.getUntrackedParameter<std::string > ("mode", "multi");
  //dirName_ = iConfig.getUntrackedParameter<std::string > ("dirName", "HcalDigisV/HcalDigiTask");
  //testNumber_= iConfig.getParameter<bool>("TestNumber");
  //hep17_     = iConfig.getParameter<bool>("hep17");

  //Collections
  //tok_hbhe_ = consumes<HBHERecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("HBHERecHitCollectionLabel"));
  //tok_hf_  = consumes<HFRecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("HFRecHitCollectionLabel"));
  //tok_ho_ = consumes<HORecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("HORecHitCollectionLabel"));
  //tok_vertexc_ = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCollectionLabel"));
  //tok_calo_ = consumes<CaloTowerCollection>(iConfig.getUntrackedParameter<edm::InputTag>("CaloTowerCollectionLabel"));

  tok_hbhe_ = consumes< HBHEDigiCollection >(inputTag_);
  tok_ho_ = consumes< HODigiCollection >(inputTag_);
  tok_hf_ = consumes< HFDigiCollection >(inputTag_);

  tok_qie10_hf_ = consumes< QIE10DigiCollection >(QIE10inputTag_);
  tok_qie11_hbhe_ = consumes< QIE11DigiCollection >(QIE11inputTag_);
  
  //now do what ever initialization is needed

  tf1 = new TFile(outputfile_.c_str(), "RECREATE");

  //Now we use the modification so that we can use the TFileService
  edm::Service<TFileService> fs;

  tt1 = fs->make<TTree>(treename_.c_str(),treename_.c_str());

  //branches

  tt1->Branch("run", &run, "run/I");
  tt1->Branch("lumi", &lumi, "lumi/I");
  tt1->Branch("event", &event, "event/I");

  tt1->Branch("recHitHB_En","std::vector<float>", &recHitHB_En, 32000, 0);
  tt1->Branch("recHitHB_ieta","std::vector<int>", &recHitHB_ieta, 32000, 0);
  tt1->Branch("recHitHB_iphi","std::vector<int>", &recHitHB_iphi, 32000, 0);
  tt1->Branch("recHitHB_depth","std::vector<int>", &recHitHB_depth, 32000, 0);
  tt1->Branch("recHitHB_time","std::vector<float>", &recHitHB_time, 32000, 0);

  tt1->Branch("recHitHE_En","std::vector<float>", &recHitHE_En, 32000, 0);
  //tt1->Branch("recHitHE_EnRAW","std::vector<float>", &recHitHE_EnRAW, 32000, 0);
  tt1->Branch("recHitHE_ieta","std::vector<int>", &recHitHE_ieta, 32000, 0);
  tt1->Branch("recHitHE_iphi","std::vector<int>", &recHitHE_iphi, 32000, 0);
  tt1->Branch("recHitHE_depth","std::vector<int>", &recHitHE_depth, 32000, 0);
  tt1->Branch("recHitHE_time","std::vector<float>", &recHitHE_time, 32000, 0);

  tt1->Branch("recHitHF_En","std::vector<float>", &recHitHF_En, 32000, 0);
  //tt1->Branch("recHitHF_EnRAW","std::vector<float>", &recHitHF_EnRAW, 32000, 0);
  tt1->Branch("recHitHF_ieta","std::vector<int>", &recHitHF_ieta, 32000, 0);
  tt1->Branch("recHitHF_iphi","std::vector<int>", &recHitHF_iphi, 32000, 0);
  tt1->Branch("recHitHF_depth","std::vector<int>", &recHitHF_depth, 32000, 0);
  tt1->Branch("recHitHF_time","std::vector<float>", &recHitHF_time, 32000, 0);

}


HcalDigiTree::~HcalDigiTree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  tf1->cd();
  tt1->Write();
  tf1->Write();
  tf1->Close();

}

//
// member functions
//

// ------------ method called for each event  ------------
void HcalDigiTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;

   edm::ESHandle<HcalDDDRecConstants> pHRNDC;
   iSetup.get<HcalRecNumberingRecord>().get( pHRNDC );
   const HcalDDDRecConstants* hcons = &(*pHRNDC);

   // HCAL channel status map ****************************************
   
   recHitHB_En.clear();
   //recHitHB_EnRAW.clear();
   recHitHB_ieta.clear();
   recHitHB_iphi.clear();
   recHitHB_depth.clear();
   recHitHB_time.clear();

   recHitHE_En.clear();
   //recHitHE_EnRAW.clear();
   recHitHE_ieta.clear();
   recHitHE_iphi.clear();
   recHitHE_depth.clear();
   recHitHE_time.clear();

   recHitHF_En.clear();
   //recHitHF_EnRAW.clear();
   recHitHF_ieta.clear();
   recHitHF_iphi.clear();
   recHitHF_depth.clear();
   recHitHF_time.clear();

   //run:lumi:event
   run = iEvent.id().run();
   lumi = iEvent.id().luminosityBlock();
   event = iEvent.id().event();

   //Clear the nvertx
   //nvertx = 0;
   //nrechits = 0;

   //Good Vertices
   //edm::Handle<VertexCollection> vertexcoll;
   //iEvent.getByToken(tok_vertexc_, vertexcoll);

   //for(VertexCollection::const_iterator vitr = vertexcoll->begin(); vitr != vertexcoll->end(); vitr++){
     //if(vitr->isFake()) continue;
     //if(vitr->ndof() <= 4) continue;
     //if(vitr->z() > 24) continue;
     //if(vitr->z() < -24) continue;
     //if(vitr->position().Rho()>=2) continue;
     //if(vitr->position().Rho()<=-2) continue;

     //nvertx++;
     
   //} //Good Vertices

   // CaloTowers
   //edm::Handle<CaloTowerCollection> towers;
   //iEvent.getByToken(tok_calo_,towers);
   //CaloTowerCollection::const_iterator cal;

   //for(cal = towers->begin(); cal != towers->end(); cal++){
   //  CaloTowerDetId idT = cal->id();
   //  caloTower_HadEt.push_back(cal->hadEt());
   //  caloTower_EmEt.push_back(cal->emEt());
   //  caloTower_ieta.push_back(idT.ieta());
   //  caloTower_iphi.push_back(idT.iphi());
   //}


   //-------------------------------------------------------------------------------------
   //HCAL DIGIS 
   //-------------------------------------------------------------------------------------

   //-------------------------------------------------------------------------------------
   //HBHE digis
   //-------------------------------------------------------------------------------------
   edm::Handle< HBHEDigiCollection > digiTag;
   iEvent.getByToken(tok_hbhe_, digiTag);
   typename edm::Handle<edm::SortedCollection<HBHEDigiCollection> > digiCollection;
   typename edm::SortedCollection<HBHEDigiCollection>::const_iterator digiItr;

   //reco<HBHEDataFrame > (iEvent, iSetup, tok_hbhe_);

   if (!digiTag.isValid()) return;
   ////   const HcalDataFrameContainer<HBHEDataFrame > *digiCollection = digiTag.product();   

   /*
   for (digiItr = digiTag->begin(); digiItr != digiTag->end(); digiItr++) { 
       HcalDetId cell(digiItr->id());
       int depth = cell.depth();
       int iphi = cell.iphi();
       int ieta = cell.ieta();
       int sub = cell.subdet();
       
       //std::cout << depth << std::endl;
       //std::cout << iphi << std::endl;
       //std::cout << ieta << std::endl;
       //std::cout << sub << std::endl;
   }
   */


   // //------------------------------------------------------------------------------------
   // //HO digis
   // //------------------------------------------------------------------------------------
   // edm::Handle< HODigiCollection > digiTag;
   // iEvent.getByToken(tok_ho_, digiTag);

   // //reco<HODataFrame > (iEvent, iSetup, tok_ho_);

   // //------------------------------------------------------------------------------------
   // //HF digis
   // //------------------------------------------------------------------------------------
   // edm::Handle< HFDigiCollection > digiTag;
   // iEvent.getByToken(tok_hf_, digiTag);

   // //reco<HFDataFrame > (iEvent, iSetup, tok_hf_);

   // //------------------------------------------------------------------------------------
   // //HF QIE10Digis
   // //------------------------------------------------------------------------------------
   // edm::Handle< QIE10DigiCollection > QIE10digiTag;
   // iEvent.getByToken(tok_qie10_hf_, QIE10digiTag);

   // //reco<QIE10DataFrame>(iEvent, iSetup, tok_qie10_hf_);

   // //------------------------------------------------------------------------------------
   // //HBHE QIE11Digis
   // //------------------------------------------------------------------------------------
   // edm::Handle< QIE11DigiCollection > QIE11digiTag;
   // iEvent.getByToken(tok_qie11_hbhe_, QIE11digiTag);

   //reco<QIE11DataFrame>(iEvent, iSetup, tok_qie11_hbhe_);

   //------------------------------------------------------------------------------------

   // //HBHE RecHits
   // edm::Handle<HBHERecHitCollection> hbhecoll;
   // iEvent.getByToken(tok_hbhe_, hbhecoll);

   // int depth = 0;
   // //int severityLevel = 0;

   // std::cout << "HBHE size: " << hbhecoll->size() << std::endl;

   // for(HBHERecHitCollection::const_iterator j=hbhecoll->begin(); j != hbhecoll->end(); j++){

   //    std::cout << "HBHE size: " << hbhecoll->size() << std::endl;

   //    HcalDetId cell;
   //    if (testNumbering_) cell = HcalHitRelabeller::relabel(j->id(),hcons);
   //    else cell = HcalDetId(j->id());
   //    depth = cell.depth();

   //    std::cout << depth << " " << cell.det()<< " " << cell.subdet() << " " << cell.ieta() << std::endl;
   //    if (testNumbering_) std::cout << "testNumbering_: " << depth << " " << cell.subdet() << std::endl;

   //    //HcalDetId cell(j->id());
   //    //depth = cell.depth();
   //    //severityLevel = hcalSevLvl( (CaloRecHit*) &*j );
   //    //if(severityLevel > 8) continue;

   //    if(cell.subdet() == HcalBarrel){
   // 	//nrechits++;
   // 	//recHitEn_HB.push_back(j->energy());
   // 	  recHitHB_En.push_back(j->energy());
   // 	 //recHitHB1_EnRAW.push_back(j->eraw());
   // 	  recHitHB_ieta.push_back(cell.ieta());
   // 	  recHitHB_iphi.push_back(cell.iphi());
   // 	  recHitHB_depth.push_back(depth);
   // 	  recHitHB_time.push_back(j->time());

   //    }//HB

   //    if(cell.subdet() == HcalEndcap){

   // 	std::cout << depth << " " << cell.subdet() << " " << j->energy() << std::endl;

   // 	//nrechits++;
   // 	//recHitEn_HE.push_back(j->energy());
   // 	  recHitHE_En.push_back(j->energy());
   // 	  //recHitHE_EnRAW.push_back(j->eraw());
   // 	  recHitHE_ieta.push_back(cell.ieta());
   // 	  recHitHE_iphi.push_back(cell.iphi());
   // 	  recHitHE_depth.push_back(depth);
   // 	  recHitHE_time.push_back(j->time());
   //    }//HE     

   //  } //HBHE

   //  //HF RecHits
   //  edm::Handle<HFRecHitCollection> hfcoll;
   //  iEvent.getByToken(tok_hf_, hfcoll);

   //  for( HFRecHitCollection::const_iterator j = hfcoll->begin(); j != hfcoll->end(); j++){
   //   HcalDetId cell(j->id());
   //   depth = cell.depth();
   //   //severityLevel = hcalSevLvl( (CaloRecHit*) &*j );

   //   //ZS emulation
   //   //int auxwd = j->aux();
   //   //bool reject = true;

   //   //for(int TSidx = 0; TSidx < 3; TSidx++){
   //     //int TS2 = (auxwd >> (TSidx*7)) & 0x7F;
   //     //int TS3 = (auxwd >> (TSidx*7+7)) & 0x7F;
   //     //if(TS2+TS3 >= 10) reject = false;
   //   //}
   //   //if(reject) continue;
   //   //ZS emulation

   //   //if(severityLevel > 8) continue;

   //   //nrechits++;
   //   //recHitEn_HF.push_back(j->energy());
   // 	 recHitHF_En.push_back(j->energy());
   // 	 //recHitHF1_EnRAW.push_back(j->eraw());
   // 	 recHitHF_ieta.push_back(cell.ieta());
   // 	 recHitHF_iphi.push_back(cell.iphi());
   // 	 recHitHF_depth.push_back(depth);
   // 	 recHitHF_time.push_back(j->time());
     
   // } //HF

   // //HO RecHits
   // //edm::Handle<HORecHitCollection> hocoll;
   // //iEvent.getByToken(tok_ho_, hocoll);
   // //for( HORecHitCollection::const_iterator j = hocoll->begin(); j != hocoll->end(); j++){
   //   //HcalDetId cell(j->id());
   //   //severityLevel = hcalSevLvl( (CaloRecHit*) &*j );
   //   //if(severityLevel > 8) continue;
   //   //nrechits++;
   //   //recHitEn_HO.push_back(j->energy());
   // 	 //recHitHO_En.push_back(j->energy());
   // 	 //recHitHO_EnRAW.push_back(j->eraw());
   // 	 //recHitHO_ieta.push_back(cell.ieta());
   // 	 //recHitHO_iphi.push_back(cell.iphi());
   // 	 //recHitHO_time.push_back(j->time());
       
   // //} //HO

   //Fill Tree
   tt1->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void HcalDigiTree::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void HcalDigiTree::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
HcalDigiTree::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
HcalDigiTree::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
HcalDigiTree::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
HcalDigiTree::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HcalDigiTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//int 
//HcalDigiTree::hcalSevLvl(const CaloRecHit* hit){
 
//   const DetId id = hit->detid();

//   const uint32_t recHitFlag = hit->flags();
//   const uint32_t dbStatusFlag = theHcalChStatus->getValues(id)->getValue();

//   int severityLevel = theHcalSevLvlComputer->getSeverityLevel(id, recHitFlag, dbStatusFlag);

//   return severityLevel;

//}

//define this as a plug-in
DEFINE_FWK_MODULE(HcalDigiTree);

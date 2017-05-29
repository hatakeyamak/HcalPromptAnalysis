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
// Original Author:  Christopher Madrid
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

  TFile *tf1;
  TTree *tt1;
  

  // run:lumi:event
  int run;
  int lumi;
  int event;

  std::vector<int> DigiHB_ieta;
  std::vector<int> DigiHB_iphi;
  std::vector<int> DigiHB_depth; 
  std::vector<int> DigiHB_sub;

  std::vector<int> DigiHE_ieta;
  std::vector<int> DigiHE_iphi;
  std::vector<int> DigiHE_depth; 
  std::vector<int> DigiHE_sub;

  bool testNumbering_;

  edm::InputTag inputTag_;
  edm::InputTag QIE10inputTag_;
  edm::InputTag QIE11inputTag_;

  edm::EDGetTokenT< HBHEDigiCollection > tok_hbhe_; 
  edm::EDGetTokenT< HODigiCollection > tok_ho_;
  edm::EDGetTokenT< HFDigiCollection > tok_hf_;

  edm::EDGetTokenT< QIE10DigiCollection > tok_qie10_hf_; 
  edm::EDGetTokenT< QIE11DigiCollection > tok_qie11_hbhe_; 

   
  
  template<class Digi> void reco(const edm::Event& iEvent, const edm::EventSetup& iSetup, const edm::EDGetTokenT<edm::SortedCollection<Digi> > &tok);
  template<class dataFrameType> void reco(const edm::Event& iEvent, const edm::EventSetup& iSetup, const edm::EDGetTokenT<HcalDataFrameContainer<dataFrameType> > &tok);
   
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
  tok_hbhe_ = consumes< HBHEDigiCollection >(inputTag_);
  tok_ho_ = consumes< HODigiCollection >(inputTag_);
  tok_hf_ = consumes< HFDigiCollection >(inputTag_);

  tok_qie10_hf_ = consumes< QIE10DigiCollection >(QIE10inputTag_);
  tok_qie11_hbhe_ = consumes< QIE11DigiCollection >(QIE11inputTag_);
 
  //Now do what ever initialization is needed

  tf1 = new TFile(outputfile_.c_str(), "RECREATE");

  //Now we use the modification so that we can use the TFileService
  edm::Service<TFileService> fs;
  ///tt1 = fs->make<TTree>("HcalDigiTree","HcalDigiTree");
  tt1 = fs->make<TTree>(treename_.c_str(),treename_.c_str());

  //Branches
  tt1->Branch("run", &run, "run/I");
  tt1->Branch("lumi", &lumi, "lumi/I");
  tt1->Branch("event", &event, "event/I");

  tt1->Branch("DigiHB_ieta","std::vector<int>", &DigiHB_ieta, 32000, 0);
  tt1->Branch("DigiHB_iphi","std::vector<int>", &DigiHB_iphi, 32000, 0);
  tt1->Branch("DigiHB_depth","std::vector<int>", &DigiHB_depth, 32000, 0);
  tt1->Branch("DigiHB_sub","std::vector<int>", &DigiHB_sub, 32000, 0);

  tt1->Branch("DigiHE_ieta","std::vector<int>", &DigiHE_ieta, 32000, 0);
  tt1->Branch("DigiHE_iphi","std::vector<int>", &DigiHE_iphi, 32000, 0);
  tt1->Branch("DigiHE_depth","std::vector<int>", &DigiHE_depth, 32000, 0);
  tt1->Branch("DigiHE_sub","std::vector<int>", &DigiHE_sub, 32000, 0);
  
  std::cout << "Made it Here 1" << std::endl;
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
void 
HcalDigiTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;

  ///edm::ESHandle<HcalDDDRecConstants> pHRNDC;
  ///iSetup.get<HcalRecNumberingRecord>().get( pHRNDC );
   ///const HcalDDDRecConstants* hcons = &(*pHRNDC);

   // HCAL channel status map ****************************************

   DigiHB_ieta.clear();
   DigiHB_iphi.clear();
   DigiHB_depth.clear();
   DigiHB_sub.clear();

   DigiHE_ieta.clear();
   DigiHE_iphi.clear();
   DigiHE_depth.clear();
   DigiHE_sub.clear();

   //run:lumi:event
   run = iEvent.id().run();
   lumi = iEvent.id().luminosityBlock();
   event = iEvent.id().event();

   std::cout << "Made it Here 2" << std::endl;

   //-------------------------------------------------------------------------------------
   //HCAL DIGIS 
   //-------------------------------------------------------------------------------------

   //-------------------------------------------------------------------------------------
   //HBHE digis
   //-------------------------------------------------------------------------------------
  
   
   edm::Handle< HBHEDigiCollection > digiTag;
   iEvent.getByToken(tok_hbhe_, digiTag);
  
   /*
   // ADC2fC
   HcalCalibrations calibrations;
   CaloSamples tool;
   iEvent.getByToken(tok, digiCollection);
   */

   std::cout << "Made it Here 3" << std::endl;

   for(HBHEDigiCollection::const_iterator j=digiTag->begin(); j != digiTag->end(); j++){
     HcalDetId cell;
     cell = HcalDetId(j->id());
     int ieta = cell.ieta();
     int iphi = cell.iphi();
     int depth = cell.depth();
     int sub = cell.subdet();
     
     ///if(cell.subdet() == HcalBarrel){
       //DigiHB_ieta.push_back(j->energy());
       DigiHB_ieta.push_back(ieta);
       DigiHB_iphi.push_back(iphi);
       DigiHB_depth.push_back(depth);
       DigiHB_sub.push_back(sub);
       
       std::cout << "Made it Here 4" << std::endl;

     ///}//HB
     /*
     if(cell.subdet() == HcalEndcap){
       //DigiHE_ieta.push_back(j->energy());
       DigiHE_ieta.push_back(ieta);
       DigiHE_iphi.push_back(iphi);
       DigiHE_depth.push_back(depth);
       DigiHE_sub.push_back(sub);

     }//HE 
     */

     /*
     for (int i = 0; i < 10; i++) {
       HcalCalibrations calibrations = conditions->getHcalCalibrations(cell);
       const HcalQIECoder* channelCoder = conditions->getHcalCoder(cell);
       const HcalQIEShape* shape = conditions->getHcalShape(channelCoder);
       HcalCoderDb coder(*channelCoder, *shape);
       coder.adc2fC(*digiItr, tool);

     }//Loop over Time Steps 
     */

   }//Loop over HBHE Digis 

   /*
   //------------------------------------------------------------------------------------
   //HO digis
   //------------------------------------------------------------------------------------
   edm::Handle< HODigiCollection > digiTag;
   iEvent.getByToken(tok_ho_, digiTag);

   //reco<HODataFrame > (iEvent, iSetup, tok_ho_);

   //------------------------------------------------------------------------------------
   //HF digis
   //------------------------------------------------------------------------------------
   edm::Handle< HFDigiCollection > digiTag;
   iEvent.getByToken(tok_hf_, digiTag);

   //reco<HFDataFrame > (iEvent, iSetup, tok_hf_);

   //------------------------------------------------------------------------------------
   //HF QIE10Digis
   //------------------------------------------------------------------------------------
   edm::Handle< QIE10DigiCollection > QIE10digiTag;
   iEvent.getByToken(tok_qie10_hf_, QIE10digiTag);

   //reco<QIE10DataFrame>(iEvent, iSetup, tok_qie10_hf_);

   //------------------------------------------------------------------------------------
   //HBHE QIE11Digis
   //------------------------------------------------------------------------------------
   edm::Handle< QIE11DigiCollection > QIE11digiTag;
   iEvent.getByToken(tok_qie11_hbhe_, QIE11digiTag);

   reco<QIE11DataFrame>(iEvent, iSetup, tok_qie11_hbhe_);

   //------------------------------------------------------------------------------------
   */

   //Fill the tree
   tt1->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void
HcalDigiTree::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HcalDigiTree::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void 
HcalDigiTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HcalDigiTree);

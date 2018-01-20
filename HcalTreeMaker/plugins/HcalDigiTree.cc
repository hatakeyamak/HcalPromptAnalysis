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

#include <vector>
#include <utility>
#include <ostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <memory> 

#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/HcalObjects/interface/HcalQIEShape.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

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
  int bx;

  std::vector<int> DigiHBHE_ieta;
  std::vector<int> DigiHBHE_iphi;
  std::vector<int> DigiHBHE_depth; 
  std::vector<int> DigiHBHE_sub;
  std::vector<float> DigiHBHE_charge0;
  std::vector<float> DigiHBHE_charge1;
  std::vector<float> DigiHBHE_charge2;
  std::vector<float> DigiHBHE_charge3;
  std::vector<float> DigiHBHE_charge4;
  std::vector<float> DigiHBHE_charge5;
  std::vector<float> DigiHBHE_charge6;
  std::vector<float> DigiHBHE_charge7;
  std::vector<float> DigiHBHE_charge8;
  std::vector<float> DigiHBHE_charge9;
  std::vector<float> DigiHBHE_adc0;
  std::vector<float> DigiHBHE_adc1;
  std::vector<float> DigiHBHE_adc2;
  std::vector<float> DigiHBHE_adc3;
  std::vector<float> DigiHBHE_adc4;
  std::vector<float> DigiHBHE_adc5;
  std::vector<float> DigiHBHE_adc6;
  std::vector<float> DigiHBHE_adc7;
  std::vector<float> DigiHBHE_adc8;
  std::vector<float> DigiHBHE_adc9;

  std::vector<int> DigiHBHE_QIE11_ieta;
  std::vector<int> DigiHBHE_QIE11_iphi;
  std::vector<int> DigiHBHE_QIE11_depth; 
  std::vector<int> DigiHBHE_QIE11_sub;
  std::vector<float> DigiHBHE_QIE11_charge0;
  std::vector<float> DigiHBHE_QIE11_charge1;
  std::vector<float> DigiHBHE_QIE11_charge2;
  std::vector<float> DigiHBHE_QIE11_charge3;
  std::vector<float> DigiHBHE_QIE11_charge4;
  std::vector<float> DigiHBHE_QIE11_charge5;
  std::vector<float> DigiHBHE_QIE11_charge6;
  std::vector<float> DigiHBHE_QIE11_charge7;
  std::vector<float> DigiHBHE_QIE11_charge8;
  std::vector<float> DigiHBHE_QIE11_charge9;
  std::vector<float> DigiHBHE_QIE11_adc0;
  std::vector<float> DigiHBHE_QIE11_adc1;
  std::vector<float> DigiHBHE_QIE11_adc2;
  std::vector<float> DigiHBHE_QIE11_adc3;
  std::vector<float> DigiHBHE_QIE11_adc4;
  std::vector<float> DigiHBHE_QIE11_adc5;
  std::vector<float> DigiHBHE_QIE11_adc6;
  std::vector<float> DigiHBHE_QIE11_adc7;
  std::vector<float> DigiHBHE_QIE11_adc8;
  std::vector<float> DigiHBHE_QIE11_adc9;

  std::vector<int> DigiHO_ieta;
  std::vector<int> DigiHO_iphi;
  std::vector<int> DigiHO_depth; 
  std::vector<int> DigiHO_sub;
  std::vector<float> DigiHO_charge0;
  std::vector<float> DigiHO_charge1;
  std::vector<float> DigiHO_charge2;
  std::vector<float> DigiHO_charge3;
  std::vector<float> DigiHO_charge4;
  std::vector<float> DigiHO_charge5;
  std::vector<float> DigiHO_charge6;
  std::vector<float> DigiHO_charge7;
  std::vector<float> DigiHO_charge8;
  std::vector<float> DigiHO_charge9;
  std::vector<float> DigiHO_adc0;
  std::vector<float> DigiHO_adc1;
  std::vector<float> DigiHO_adc2;
  std::vector<float> DigiHO_adc3;
  std::vector<float> DigiHO_adc4;
  std::vector<float> DigiHO_adc5;
  std::vector<float> DigiHO_adc6;
  std::vector<float> DigiHO_adc7;
  std::vector<float> DigiHO_adc8;
  std::vector<float> DigiHO_adc9;

  std::vector<int> DigiHF_ieta;
  std::vector<int> DigiHF_iphi;
  std::vector<int> DigiHF_depth; 
  std::vector<int> DigiHF_sub;
  std::vector<float> DigiHF_charge0;
  std::vector<float> DigiHF_charge1;
  std::vector<float> DigiHF_charge2;
  std::vector<float> DigiHF_charge3;
  std::vector<float> DigiHF_charge4;
  std::vector<float> DigiHF_charge5;
  std::vector<float> DigiHF_charge6;
  std::vector<float> DigiHF_charge7;
  std::vector<float> DigiHF_charge8;
  std::vector<float> DigiHF_charge9;
  std::vector<float> DigiHF_adc0;
  std::vector<float> DigiHF_adc1;
  std::vector<float> DigiHF_adc2;
  std::vector<float> DigiHF_adc3;
  std::vector<float> DigiHF_adc4;
  std::vector<float> DigiHF_adc5;
  std::vector<float> DigiHF_adc6;
  std::vector<float> DigiHF_adc7;
  std::vector<float> DigiHF_adc8;
  std::vector<float> DigiHF_adc9;
  
  std::vector<int> DigiHF_QIE10_ieta;
  std::vector<int> DigiHF_QIE10_iphi;
  std::vector<int> DigiHF_QIE10_depth; 
  std::vector<int> DigiHF_QIE10_sub;
  std::vector<float> DigiHF_QIE10_charge0;
  std::vector<float> DigiHF_QIE10_charge1;
  std::vector<float> DigiHF_QIE10_charge2;
  std::vector<float> DigiHF_QIE10_charge3;
  std::vector<float> DigiHF_QIE10_charge4;
  std::vector<float> DigiHF_QIE10_charge5;
  std::vector<float> DigiHF_QIE10_charge6;
  std::vector<float> DigiHF_QIE10_charge7;
  std::vector<float> DigiHF_QIE10_charge8;
  std::vector<float> DigiHF_QIE10_charge9;
  std::vector<float> DigiHF_QIE10_adc0;
  std::vector<float> DigiHF_QIE10_adc1;
  std::vector<float> DigiHF_QIE10_adc2;
  std::vector<float> DigiHF_QIE10_adc3;
  std::vector<float> DigiHF_QIE10_adc4;
  std::vector<float> DigiHF_QIE10_adc5;
  std::vector<float> DigiHF_QIE10_adc6;
  std::vector<float> DigiHF_QIE10_adc7;
  std::vector<float> DigiHF_QIE10_adc8;
  std::vector<float> DigiHF_QIE10_adc9;

  bool testNumbering_;
 
  edm::InputTag inputTag_;
  edm::InputTag QIE10inputTag_;
  edm::InputTag QIE11inputTag_;

  edm::EDGetTokenT< HBHEDigiCollection > tok_hbhe_; 
  edm::EDGetTokenT< HODigiCollection > tok_ho_;
  edm::EDGetTokenT< HFDigiCollection > tok_hf_;
  edm::EDGetTokenT< QIE10DigiCollection > tok_qie10_hf_; 
  edm::EDGetTokenT< QIE11DigiCollection > tok_qie11_hbhe_; 

  edm::ESHandle<HcalDbService> conditions;
  
  //template<class Digi> void reco(const edm::Event& iEvent, const edm::EventSetup& iSetup, const edm::EDGetTokenT<edm::SortedCollection<Digi> > &tok);
  //template<class dataFrameType> void reco(const edm::Event& iEvent, const edm::EventSetup& iSetup, const edm::EDGetTokenT<HcalDataFrameContainer<dataFrameType> > &tok);
   
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
  //inputLabel_ = iConfig.getParameter<std::string > ("digiLabel");
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

  //KH tf1 = new TFile(outputfile_.c_str(), "RECREATE");

  //Now we use the modification so that we can use the TFileService
  edm::Service<TFileService> fs;
  tt1 = fs->make<TTree>(treename_.c_str(),treename_.c_str());

  //Branches
  tt1->Branch("run", &run, "run/I");
  tt1->Branch("lumi", &lumi, "lumi/I");
  tt1->Branch("event", &event, "event/I");

  tt1->Branch("bx", &bx, "bx/I");

  tt1->Branch("DigiHBHE_ieta","std::vector<int>", &DigiHBHE_ieta, 32000, 0);
  tt1->Branch("DigiHBHE_iphi","std::vector<int>", &DigiHBHE_iphi, 32000, 0);
  tt1->Branch("DigiHBHE_depth","std::vector<int>", &DigiHBHE_depth, 32000, 0);
  tt1->Branch("DigiHBHE_sub","std::vector<int>", &DigiHBHE_sub, 32000, 0);
  tt1->Branch("DigiHBHE_charge0","std::vector<float>", &DigiHBHE_charge0, 32000, 0);
  tt1->Branch("DigiHBHE_charge1","std::vector<float>", &DigiHBHE_charge1, 32000, 0);
  tt1->Branch("DigiHBHE_charge2","std::vector<float>", &DigiHBHE_charge2, 32000, 0);
  tt1->Branch("DigiHBHE_charge3","std::vector<float>", &DigiHBHE_charge3, 32000, 0);
  tt1->Branch("DigiHBHE_charge4","std::vector<float>", &DigiHBHE_charge4, 32000, 0);
  tt1->Branch("DigiHBHE_charge5","std::vector<float>", &DigiHBHE_charge5, 32000, 0);
  tt1->Branch("DigiHBHE_charge6","std::vector<float>", &DigiHBHE_charge6, 32000, 0);
  tt1->Branch("DigiHBHE_charge7","std::vector<float>", &DigiHBHE_charge7, 32000, 0);
  tt1->Branch("DigiHBHE_charge8","std::vector<float>", &DigiHBHE_charge8, 32000, 0);
  tt1->Branch("DigiHBHE_charge9","std::vector<float>", &DigiHBHE_charge9, 32000, 0);
  tt1->Branch("DigiHBHE_adc0","std::vector<float>", &DigiHBHE_adc0, 32000, 0);
  tt1->Branch("DigiHBHE_adc1","std::vector<float>", &DigiHBHE_adc1, 32000, 0);
  tt1->Branch("DigiHBHE_adc2","std::vector<float>", &DigiHBHE_adc2, 32000, 0);
  tt1->Branch("DigiHBHE_adc3","std::vector<float>", &DigiHBHE_adc3, 32000, 0);
  tt1->Branch("DigiHBHE_adc4","std::vector<float>", &DigiHBHE_adc4, 32000, 0);
  tt1->Branch("DigiHBHE_adc5","std::vector<float>", &DigiHBHE_adc5, 32000, 0);
  tt1->Branch("DigiHBHE_adc6","std::vector<float>", &DigiHBHE_adc6, 32000, 0);
  tt1->Branch("DigiHBHE_adc7","std::vector<float>", &DigiHBHE_adc7, 32000, 0);
  tt1->Branch("DigiHBHE_adc8","std::vector<float>", &DigiHBHE_adc8, 32000, 0);
  tt1->Branch("DigiHBHE_adc9","std::vector<float>", &DigiHBHE_adc9, 32000, 0);
 
  tt1->Branch("DigiHBHE_QIE11_ieta","std::vector<int>", &DigiHBHE_QIE11_ieta, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_iphi","std::vector<int>", &DigiHBHE_QIE11_iphi, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_depth","std::vector<int>", &DigiHBHE_QIE11_depth, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_sub","std::vector<int>", &DigiHBHE_QIE11_sub, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_charge0","std::vector<float>", &DigiHBHE_QIE11_charge0, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_charge1","std::vector<float>", &DigiHBHE_QIE11_charge1, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_charge2","std::vector<float>", &DigiHBHE_QIE11_charge2, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_charge3","std::vector<float>", &DigiHBHE_QIE11_charge3, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_charge4","std::vector<float>", &DigiHBHE_QIE11_charge4, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_charge5","std::vector<float>", &DigiHBHE_QIE11_charge5, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_charge6","std::vector<float>", &DigiHBHE_QIE11_charge6, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_charge7","std::vector<float>", &DigiHBHE_QIE11_charge7, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_charge8","std::vector<float>", &DigiHBHE_QIE11_charge8, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_charge9","std::vector<float>", &DigiHBHE_QIE11_charge9, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_adc0","std::vector<float>", &DigiHBHE_QIE11_adc0, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_adc1","std::vector<float>", &DigiHBHE_QIE11_adc1, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_adc2","std::vector<float>", &DigiHBHE_QIE11_adc2, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_adc3","std::vector<float>", &DigiHBHE_QIE11_adc3, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_adc4","std::vector<float>", &DigiHBHE_QIE11_adc4, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_adc5","std::vector<float>", &DigiHBHE_QIE11_adc5, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_adc6","std::vector<float>", &DigiHBHE_QIE11_adc6, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_adc7","std::vector<float>", &DigiHBHE_QIE11_adc7, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_adc8","std::vector<float>", &DigiHBHE_QIE11_adc8, 32000, 0);
  tt1->Branch("DigiHBHE_QIE11_adc9","std::vector<float>", &DigiHBHE_QIE11_adc9, 32000, 0);

  tt1->Branch("DigiHO_ieta","std::vector<int>", &DigiHO_ieta, 32000, 0);
  tt1->Branch("DigiHO_iphi","std::vector<int>", &DigiHO_iphi, 32000, 0);
  tt1->Branch("DigiHO_depth","std::vector<int>", &DigiHO_depth, 32000, 0);
  tt1->Branch("DigiHO_sub","std::vector<int>", &DigiHO_sub, 32000, 0);
  tt1->Branch("DigiHO_charge0","std::vector<float>", &DigiHO_charge0, 32000, 0);
  tt1->Branch("DigiHO_charge1","std::vector<float>", &DigiHO_charge1, 32000, 0);
  tt1->Branch("DigiHO_charge2","std::vector<float>", &DigiHO_charge2, 32000, 0);
  tt1->Branch("DigiHO_charge3","std::vector<float>", &DigiHO_charge3, 32000, 0);
  tt1->Branch("DigiHO_charge4","std::vector<float>", &DigiHO_charge4, 32000, 0);
  tt1->Branch("DigiHO_charge5","std::vector<float>", &DigiHO_charge5, 32000, 0);
  tt1->Branch("DigiHO_charge6","std::vector<float>", &DigiHO_charge6, 32000, 0);
  tt1->Branch("DigiHO_charge7","std::vector<float>", &DigiHO_charge7, 32000, 0);
  tt1->Branch("DigiHO_charge8","std::vector<float>", &DigiHO_charge8, 32000, 0);
  tt1->Branch("DigiHO_charge9","std::vector<float>", &DigiHO_charge9, 32000, 0);
  tt1->Branch("DigiHO_adc0","std::vector<float>", &DigiHO_adc0, 32000, 0);
  tt1->Branch("DigiHO_adc1","std::vector<float>", &DigiHO_adc1, 32000, 0);
  tt1->Branch("DigiHO_adc2","std::vector<float>", &DigiHO_adc2, 32000, 0);
  tt1->Branch("DigiHO_adc3","std::vector<float>", &DigiHO_adc3, 32000, 0);
  tt1->Branch("DigiHO_adc4","std::vector<float>", &DigiHO_adc4, 32000, 0);
  tt1->Branch("DigiHO_adc5","std::vector<float>", &DigiHO_adc5, 32000, 0);
  tt1->Branch("DigiHO_adc6","std::vector<float>", &DigiHO_adc6, 32000, 0);
  tt1->Branch("DigiHO_adc7","std::vector<float>", &DigiHO_adc7, 32000, 0);
  tt1->Branch("DigiHO_adc8","std::vector<float>", &DigiHO_adc8, 32000, 0);
  tt1->Branch("DigiHO_adc9","std::vector<float>", &DigiHO_adc9, 32000, 0);

  /*
  tt1->Branch("DigiHF_ieta","std::vector<int>", &DigiHF_ieta, 32000, 0);
  tt1->Branch("DigiHF_iphi","std::vector<int>", &DigiHF_iphi, 32000, 0);
  tt1->Branch("DigiHF_depth","std::vector<int>", &DigiHF_depth, 32000, 0);
  tt1->Branch("DigiHF_sub","std::vector<int>", &DigiHF_sub, 32000, 0);
  tt1->Branch("DigiHF_charge0","std::vector<float>", &DigiHF_charge0, 32000, 0);
  tt1->Branch("DigiHF_charge1","std::vector<float>", &DigiHF_charge1, 32000, 0);
  tt1->Branch("DigiHF_charge2","std::vector<float>", &DigiHF_charge2, 32000, 0);
  tt1->Branch("DigiHF_charge3","std::vector<float>", &DigiHF_charge3, 32000, 0);
  tt1->Branch("DigiHF_charge4","std::vector<float>", &DigiHF_charge4, 32000, 0);
  tt1->Branch("DigiHF_charge5","std::vector<float>", &DigiHF_charge5, 32000, 0);
  tt1->Branch("DigiHF_charge6","std::vector<float>", &DigiHF_charge6, 32000, 0);
  tt1->Branch("DigiHF_charge7","std::vector<float>", &DigiHF_charge7, 32000, 0);
  tt1->Branch("DigiHF_charge8","std::vector<float>", &DigiHF_charge8, 32000, 0);
  tt1->Branch("DigiHF_charge9","std::vector<float>", &DigiHF_charge9, 32000, 0);
  tt1->Branch("DigiHF_adc0","std::vector<float>", &DigiHF_adc0, 32000, 0);
  tt1->Branch("DigiHF_adc1","std::vector<float>", &DigiHF_adc1, 32000, 0);
  tt1->Branch("DigiHF_adc2","std::vector<float>", &DigiHF_adc2, 32000, 0);
  tt1->Branch("DigiHF_adc3","std::vector<float>", &DigiHF_adc3, 32000, 0);
  tt1->Branch("DigiHF_adc4","std::vector<float>", &DigiHF_adc4, 32000, 0);
  tt1->Branch("DigiHF_adc5","std::vector<float>", &DigiHF_adc5, 32000, 0);
  tt1->Branch("DigiHF_adc6","std::vector<float>", &DigiHF_adc6, 32000, 0);
  tt1->Branch("DigiHF_adc7","std::vector<float>", &DigiHF_adc7, 32000, 0);
  tt1->Branch("DigiHF_adc8","std::vector<float>", &DigiHF_adc8, 32000, 0);
  tt1->Branch("DigiHF_adc9","std::vector<float>", &DigiHF_adc9, 32000, 0);
  */

  tt1->Branch("DigiHF_QIE10_ieta","std::vector<int>", &DigiHF_QIE10_ieta, 32000, 0);
  tt1->Branch("DigiHF_QIE10_iphi","std::vector<int>", &DigiHF_QIE10_iphi, 32000, 0);
  tt1->Branch("DigiHF_QIE10_depth","std::vector<int>", &DigiHF_QIE10_depth, 32000, 0);
  tt1->Branch("DigiHF_QIE10_sub","std::vector<int>", &DigiHF_QIE10_sub, 32000, 0);
  tt1->Branch("DigiHF_QIE10_charge0","std::vector<float>", &DigiHF_QIE10_charge0, 32000, 0);
  tt1->Branch("DigiHF_QIE10_charge1","std::vector<float>", &DigiHF_QIE10_charge1, 32000, 0);
  tt1->Branch("DigiHF_QIE10_charge2","std::vector<float>", &DigiHF_QIE10_charge2, 32000, 0);
  tt1->Branch("DigiHF_QIE10_charge3","std::vector<float>", &DigiHF_QIE10_charge3, 32000, 0);
  tt1->Branch("DigiHF_QIE10_charge4","std::vector<float>", &DigiHF_QIE10_charge4, 32000, 0);
  tt1->Branch("DigiHF_QIE10_charge5","std::vector<float>", &DigiHF_QIE10_charge5, 32000, 0);
  tt1->Branch("DigiHF_QIE10_charge6","std::vector<float>", &DigiHF_QIE10_charge6, 32000, 0);
  tt1->Branch("DigiHF_QIE10_charge7","std::vector<float>", &DigiHF_QIE10_charge7, 32000, 0);
  tt1->Branch("DigiHF_QIE10_charge8","std::vector<float>", &DigiHF_QIE10_charge8, 32000, 0);
  tt1->Branch("DigiHF_QIE10_charge9","std::vector<float>", &DigiHF_QIE10_charge9, 32000, 0);
  tt1->Branch("DigiHF_QIE10_adc0","std::vector<float>", &DigiHF_QIE10_adc0, 32000, 0);
  tt1->Branch("DigiHF_QIE10_adc1","std::vector<float>", &DigiHF_QIE10_adc1, 32000, 0);
  tt1->Branch("DigiHF_QIE10_adc2","std::vector<float>", &DigiHF_QIE10_adc2, 32000, 0);
  tt1->Branch("DigiHF_QIE10_adc3","std::vector<float>", &DigiHF_QIE10_adc3, 32000, 0);
  tt1->Branch("DigiHF_QIE10_adc4","std::vector<float>", &DigiHF_QIE10_adc4, 32000, 0);
  tt1->Branch("DigiHF_QIE10_adc5","std::vector<float>", &DigiHF_QIE10_adc5, 32000, 0);
  tt1->Branch("DigiHF_QIE10_adc6","std::vector<float>", &DigiHF_QIE10_adc6, 32000, 0);
  tt1->Branch("DigiHF_QIE10_adc7","std::vector<float>", &DigiHF_QIE10_adc7, 32000, 0);
  tt1->Branch("DigiHF_QIE10_adc8","std::vector<float>", &DigiHF_QIE10_adc8, 32000, 0);
  tt1->Branch("DigiHF_QIE10_adc9","std::vector<float>", &DigiHF_QIE10_adc9, 32000, 0);

  //std::cout << "Made it Here 1" << std::endl;
}

HcalDigiTree::~HcalDigiTree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  
  //tf1->cd();
  //tt1->Write();
  //tf1->Write();
  //tf1->Close();
 
}

//
// member functions
//

// ------------ method called for each event  ------------
void 
HcalDigiTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  //using namespace reco;

  //edm::ESHandle<HcalDDDRecConstants> pHRNDC;
  //iSetup.get<HcalRecNumberingRecord>().get( pHRNDC );
  //const HcalDDDRecConstants* hcons = &(*pHRNDC);

   // HCAL channel status map ****************************************

   DigiHBHE_ieta.clear();
   DigiHBHE_iphi.clear();
   DigiHBHE_depth.clear();
   DigiHBHE_sub.clear();
   DigiHBHE_charge0.clear();
   DigiHBHE_charge1.clear();
   DigiHBHE_charge2.clear();
   DigiHBHE_charge3.clear();
   DigiHBHE_charge4.clear();
   DigiHBHE_charge5.clear();
   DigiHBHE_charge6.clear();
   DigiHBHE_charge7.clear();
   DigiHBHE_charge8.clear();
   DigiHBHE_charge9.clear();
   DigiHBHE_adc0.clear();
   DigiHBHE_adc1.clear();
   DigiHBHE_adc2.clear();
   DigiHBHE_adc3.clear();
   DigiHBHE_adc4.clear();
   DigiHBHE_adc5.clear();
   DigiHBHE_adc6.clear();
   DigiHBHE_adc7.clear();
   DigiHBHE_adc8.clear();
   DigiHBHE_adc9.clear();

   DigiHBHE_QIE11_ieta.clear();
   DigiHBHE_QIE11_iphi.clear();
   DigiHBHE_QIE11_depth.clear();
   DigiHBHE_QIE11_sub.clear();
   DigiHBHE_QIE11_charge0.clear();
   DigiHBHE_QIE11_charge1.clear();
   DigiHBHE_QIE11_charge2.clear();
   DigiHBHE_QIE11_charge3.clear();
   DigiHBHE_QIE11_charge4.clear();
   DigiHBHE_QIE11_charge5.clear();
   DigiHBHE_QIE11_charge6.clear();
   DigiHBHE_QIE11_charge7.clear();
   DigiHBHE_QIE11_charge8.clear();
   DigiHBHE_QIE11_charge9.clear();
   DigiHBHE_QIE11_adc0.clear();
   DigiHBHE_QIE11_adc1.clear();
   DigiHBHE_QIE11_adc2.clear();
   DigiHBHE_QIE11_adc3.clear();
   DigiHBHE_QIE11_adc4.clear();
   DigiHBHE_QIE11_adc5.clear();
   DigiHBHE_QIE11_adc6.clear();
   DigiHBHE_QIE11_adc7.clear();
   DigiHBHE_QIE11_adc8.clear();
   DigiHBHE_QIE11_adc9.clear();

   DigiHO_ieta.clear();
   DigiHO_iphi.clear();
   DigiHO_depth.clear();
   DigiHO_sub.clear();
   DigiHO_charge0.clear();
   DigiHO_charge1.clear();
   DigiHO_charge2.clear();
   DigiHO_charge3.clear();
   DigiHO_charge4.clear();
   DigiHO_charge5.clear();
   DigiHO_charge6.clear();
   DigiHO_charge7.clear();
   DigiHO_charge8.clear();
   DigiHO_charge9.clear();
   DigiHO_adc0.clear();
   DigiHO_adc1.clear();
   DigiHO_adc2.clear();
   DigiHO_adc3.clear();
   DigiHO_adc4.clear();
   DigiHO_adc5.clear();
   DigiHO_adc6.clear();
   DigiHO_adc7.clear();
   DigiHO_adc8.clear();
   DigiHO_adc9.clear();

   DigiHF_ieta.clear();
   DigiHF_iphi.clear();
   DigiHF_depth.clear();
   DigiHF_sub.clear();
   DigiHF_charge0.clear();
   DigiHF_charge1.clear();
   DigiHF_charge2.clear();
   DigiHF_charge3.clear();
   DigiHF_charge4.clear();
   DigiHF_charge5.clear();
   DigiHF_charge6.clear();
   DigiHF_charge7.clear();
   DigiHF_charge8.clear();
   DigiHF_charge9.clear();
   DigiHF_adc0.clear();
   DigiHF_adc1.clear();
   DigiHF_adc2.clear();
   DigiHF_adc3.clear();
   DigiHF_adc4.clear();
   DigiHF_adc5.clear();
   DigiHF_adc6.clear();
   DigiHF_adc7.clear();
   DigiHF_adc8.clear();
   DigiHF_adc9.clear();
   
   DigiHF_QIE10_ieta.clear();
   DigiHF_QIE10_iphi.clear();
   DigiHF_QIE10_depth.clear();
   DigiHF_QIE10_sub.clear();
   DigiHF_QIE10_charge0.clear();
   DigiHF_QIE10_charge1.clear();
   DigiHF_QIE10_charge2.clear();
   DigiHF_QIE10_charge3.clear();
   DigiHF_QIE10_charge4.clear();
   DigiHF_QIE10_charge5.clear();
   DigiHF_QIE10_charge6.clear();
   DigiHF_QIE10_charge7.clear();
   DigiHF_QIE10_charge8.clear();
   DigiHF_QIE10_charge9.clear();
   DigiHF_QIE10_adc0.clear();
   DigiHF_QIE10_adc1.clear();
   DigiHF_QIE10_adc2.clear();
   DigiHF_QIE10_adc3.clear();
   DigiHF_QIE10_adc4.clear();
   DigiHF_QIE10_adc5.clear();
   DigiHF_QIE10_adc6.clear();
   DigiHF_QIE10_adc7.clear();
   DigiHF_QIE10_adc8.clear();
   DigiHF_QIE10_adc9.clear();

   //run:lumi:event
   run = iEvent.id().run();
   lumi = iEvent.id().luminosityBlock();
   event = iEvent.id().event();

   //Bunch Crossing??
   edm::EventBase const & eventbase = iEvent;
   //unsigned int bx = (unsigned int)eventbase.bunchCrossing();
   bx = eventbase.bunchCrossing();

   //std::cout << "Made it Here 2" << std::endl;

   //-------------------------------------------------------------------------------------
   //HCAL DIGIS 
   //-------------------------------------------------------------------------------------

   // ADC2fC                                                                                                 
   iSetup.get<HcalDbRecord > ().get(conditions);
   HcalCalibrations calibrations;
   CaloSamples tool;
   //iEvent.getByToken(tok, digiCollection);

   //-------------------------------------------------------------------------------------
   //HBHE digis
   //-------------------------------------------------------------------------------------
   edm::Handle< HBHEDigiCollection >HBHEdigiTag;
   iEvent.getByToken(tok_hbhe_, HBHEdigiTag);
    
   //std::cout << "Made it Here 3" << std::endl;

   for(HBHEDigiCollection::const_iterator j=HBHEdigiTag->begin(); j != HBHEdigiTag->end(); j++){
     HcalDetId cell;
     cell = HcalDetId(j->id());
     int ieta = cell.ieta();
     int iphi = cell.iphi();
     int depth = cell.depth();
     int sub = cell.subdet();

     DigiHBHE_ieta.push_back(ieta);
     DigiHBHE_iphi.push_back(iphi);
     DigiHBHE_depth.push_back(depth);
     DigiHBHE_sub.push_back(sub);

     //std::cout << ieta << " " << iphi << " " << depth << " " << sub << std::endl;
       
     HcalCalibrations calibrations = conditions->getHcalCalibrations(cell);
     const HcalQIECoder* channelCoder = conditions->getHcalCoder(cell);
     const HcalQIEShape* shape = conditions->getHcalShape(channelCoder);
     HcalCoderDb coder(*channelCoder, *shape);
     coder.adc2fC(*j, tool);

     for (int ii = 0; ii < tool.size(); ii++) {
       int capid = (*j)[ii].capid();
       // single ts amplitude                                                                               
       double val = (tool[ii] - calibrations.pedestal(capid));
       //std::cout << "val: " << ii << ": " << val << std::endl;

       if(ii==0){ 
	 DigiHBHE_charge0.push_back(val);      
	 DigiHBHE_adc0.push_back((*j)[ii].adc());
       }
       if(ii==1){
	 DigiHBHE_charge1.push_back(val);
	 DigiHBHE_adc1.push_back((*j)[ii].adc());
       }
       if(ii==2){ 
	 DigiHBHE_charge2.push_back(val);
	 DigiHBHE_adc2.push_back((*j)[ii].adc());
       }
       if(ii==3){ 
	 DigiHBHE_charge3.push_back(val);
	 DigiHBHE_adc3.push_back((*j)[ii].adc());
       }
       if(ii==4){ 
	 DigiHBHE_charge4.push_back(val);
	 DigiHBHE_adc4.push_back((*j)[ii].adc());
       }
       if(ii==5){ 
	 DigiHBHE_charge5.push_back(val);
	 DigiHBHE_adc5.push_back((*j)[ii].adc());
	 if(val > 50){ 
	   //std::cout << "ieta,iphi,depth,pedestal  : " << ieta<<"   " << iphi <<"   " << depth <<"   " <<calibrations.pedestal(capid) << std::endl;
	   //std::cout << "tool-pedestal=val_TS5 : "<< tool[ii] << " - " << calibrations.pedestal(capid) <<" = "<< val << std::endl;  
	 }
       }
       if(ii==6){ 
	 DigiHBHE_charge6.push_back(val);
	 DigiHBHE_adc6.push_back((*j)[ii].adc());
       }
       if(ii==7){ 
	 DigiHBHE_charge7.push_back(val);
	 DigiHBHE_adc7.push_back((*j)[ii].adc());
       }
       if(ii==8){ 
	 DigiHBHE_charge8.push_back(val);
	 DigiHBHE_adc8.push_back((*j)[ii].adc());
       }
       if(ii==9){ 
	 DigiHBHE_charge9.push_back(val);
	 DigiHBHE_adc9.push_back((*j)[ii].adc());
       }
       
     }//Loop to get Charge

   }//Loop over HBHE Digis 

   
   //------------------------------------------------------------------------------------
   //HBHE QIE11Digis
   //------------------------------------------------------------------------------------
   edm::Handle< QIE11DigiCollection > QIE11digiTag;
   iEvent.getByToken(tok_qie11_hbhe_, QIE11digiTag);

   for(QIE11DigiCollection::const_iterator j=QIE11digiTag->begin(); j != QIE11digiTag->end(); j++){  

     QIE11DataFrame dataFrame = *j;

     HcalDetId cell;
     cell = HcalDetId(j->id());
     int ieta = cell.ieta();
     int iphi = cell.iphi();
     int depth = cell.depth();
     int sub = cell.subdet();

     DigiHBHE_QIE11_ieta.push_back(ieta);
     DigiHBHE_QIE11_iphi.push_back(iphi);
     DigiHBHE_QIE11_depth.push_back(depth);
     DigiHBHE_QIE11_sub.push_back(sub);

     //std::cout << ieta << " " << iphi << " " << depth << " " << sub << std::endl;

     HcalCalibrations calibrations = conditions->getHcalCalibrations(cell);
     const HcalQIECoder* channelCoder = conditions->getHcalCoder(cell);
     const HcalQIEShape* shape = conditions->getHcalShape(channelCoder);
     HcalCoderDb coder(*channelCoder, *shape);
     coder.adc2fC(dataFrame, tool);

     for (int ii = 0; ii < tool.size(); ii++) {
       int capid = (dataFrame)[ii].capid();
       // single ts amplitude                                                                               
       double val = (tool[ii] - calibrations.pedestal(capid));

       //std::cout << tool[ii] <<"-" << calibrations.pedestal(capid) << "=" << val << std::endl;

       if(ii==0){ 
	 DigiHBHE_QIE11_charge0.push_back(val);      
	 DigiHBHE_QIE11_adc0.push_back((dataFrame)[ii].adc());
       }
       if(ii==1){
	 DigiHBHE_QIE11_charge1.push_back(val);
	 DigiHBHE_QIE11_adc1.push_back((dataFrame)[ii].adc());
       }
       if(ii==2){ 
	 DigiHBHE_QIE11_charge2.push_back(val);
	 DigiHBHE_QIE11_adc2.push_back((dataFrame)[ii].adc());
       }
       if(ii==3){ 
	 DigiHBHE_QIE11_charge3.push_back(val);
	 DigiHBHE_QIE11_adc3.push_back((dataFrame)[ii].adc());
       }
       if(ii==4){ 
	 DigiHBHE_QIE11_charge4.push_back(val);
	 DigiHBHE_QIE11_adc4.push_back((dataFrame)[ii].adc());
       }
       if(ii==5){ 
	 DigiHBHE_QIE11_charge5.push_back(val);
	 DigiHBHE_QIE11_adc5.push_back((dataFrame)[ii].adc());
       }
       if(ii==6){ 
	 DigiHBHE_QIE11_charge6.push_back(val);
	 DigiHBHE_QIE11_adc6.push_back((dataFrame)[ii].adc());
       }
       if(ii==7){ 
	 DigiHBHE_QIE11_charge7.push_back(val);
	 DigiHBHE_QIE11_adc7.push_back((dataFrame)[ii].adc());
       }
       if(ii==8){ 
	 DigiHBHE_QIE11_charge8.push_back(val);
	 DigiHBHE_QIE11_adc8.push_back((dataFrame)[ii].adc());
       }
       if(ii==9){ 
	 DigiHBHE_QIE11_charge9.push_back(val);
	 DigiHBHE_QIE11_adc9.push_back((dataFrame)[ii].adc());
       }
       
     }//Loop to get Charge

   }//Loop over HBHE_QIE11 Digis
   

   //------------------------------------------------------------------------------------
   //HO digis
   //------------------------------------------------------------------------------------
   edm::Handle< HODigiCollection > HOdigiTag;
   iEvent.getByToken(tok_ho_, HOdigiTag);
   //reco<HODataFrame > (iEvent, iSetup, tok_ho_);

   for(HODigiCollection::const_iterator j=HOdigiTag->begin(); j != HOdigiTag->end(); j++){  
     HcalDetId cell;
     cell = HcalDetId(j->id());
     int ieta = cell.ieta();
     int iphi = cell.iphi();
     int depth = cell.depth();
     int sub = cell.subdet();

     DigiHO_ieta.push_back(ieta);
     DigiHO_iphi.push_back(iphi);
     DigiHO_depth.push_back(depth);
     DigiHO_sub.push_back(sub);

     //std::cout << ieta << " " << iphi << " " << depth << " " << sub << std::endl;

     HcalCalibrations calibrations = conditions->getHcalCalibrations(cell);
     const HcalQIECoder* channelCoder = conditions->getHcalCoder(cell);
     const HcalQIEShape* shape = conditions->getHcalShape(channelCoder);
     HcalCoderDb coder(*channelCoder, *shape);
     coder.adc2fC(*j, tool);

     for (int ii = 0; ii < tool.size(); ii++) {
       int capid = (*j)[ii].capid();
       // single ts amplitude                                                                               
       double val = (tool[ii] - calibrations.pedestal(capid));
       //std::cout << "val: " << ii << ": " << val << std::endl;
      
       if(ii==0){ 
	 DigiHO_charge0.push_back(val);      
	 DigiHO_adc0.push_back((*j)[ii].adc());
       }
       if(ii==1){
	 DigiHO_charge1.push_back(val);
	 DigiHO_adc1.push_back((*j)[ii].adc());
       }
       if(ii==2){ 
	 DigiHO_charge2.push_back(val);
	 DigiHO_adc2.push_back((*j)[ii].adc());
       }
       if(ii==3){ 
	 DigiHO_charge3.push_back(val);
	 DigiHO_adc3.push_back((*j)[ii].adc());
       }
       if(ii==4){ 
	 DigiHO_charge4.push_back(val);
	 DigiHO_adc4.push_back((*j)[ii].adc());
       }
       if(ii==5){ 
	 DigiHO_charge5.push_back(val);
	 DigiHO_adc5.push_back((*j)[ii].adc());
       }
       if(ii==6){ 
	 DigiHO_charge6.push_back(val);
	 DigiHO_adc6.push_back((*j)[ii].adc());
       }
       if(ii==7){ 
	 DigiHO_charge7.push_back(val);
	 DigiHO_adc7.push_back((*j)[ii].adc());
       }
       if(ii==8){ 
	 DigiHO_charge8.push_back(val);
	 DigiHO_adc8.push_back((*j)[ii].adc());
       }
       if(ii==9){ 
	 DigiHO_charge9.push_back(val);
	 DigiHO_adc9.push_back((*j)[ii].adc());
       }

     }//Loop to get Charge

   }//Loop over HO Digis

   /*
   //------------------------------------------------------------------------------------
   //HF digis
   //------------------------------------------------------------------------------------
   edm::Handle< HFDigiCollection > HFdigiTag;
   iEvent.getByToken(tok_hf_, HFdigiTag);
   //reco<HFDataFrame > (iEvent, iSetup, tok_hf_);

   for(HFDigiCollection::const_iterator j=HFdigiTag->begin(); j != HFdigiTag->end(); j++){  
     HcalDetId cell;
     cell = HcalDetId(j->id());
     int ieta = cell.ieta();
     int iphi = cell.iphi();
     int depth = cell.depth();
     int sub = cell.subdet();

     DigiHF_ieta.push_back(ieta);
     DigiHF_iphi.push_back(iphi);
     DigiHF_depth.push_back(depth);
     DigiHF_sub.push_back(sub);

     //std::cout << ieta << " " << iphi << " " << depth << " " << sub << std::endl;

     HcalCalibrations calibrations = conditions->getHcalCalibrations(cell);
     const HcalQIECoder* channelCoder = conditions->getHcalCoder(cell);
     const HcalQIEShape* shape = conditions->getHcalShape(channelCoder);
     HcalCoderDb coder(*channelCoder, *shape);
     coder.adc2fC(*j, tool);

     for (int ii = 0; ii < tool.size(); ii++) {
       int capid = (*j)[ii].capid();
       // single ts amplitude                                                                               
       double val = (tool[ii] - calibrations.pedestal(capid));
       //std::cout << "val: " << ii << ": " << val << std::endl;
      
       if(ii==0){ 
	 DigiHF_charge0.push_back(val);      
	 DigiHF_adc0.push_back((*j)[ii].adc());
       }
       if(ii==1){
	 DigiHF_charge1.push_back(val);
	 DigiHF_adc1.push_back((*j)[ii].adc());
       }
       if(ii==2){ 
	 DigiHF_charge2.push_back(val);
	 DigiHF_adc2.push_back((*j)[ii].adc());
       }
       if(ii==3){ 
	 DigiHF_charge3.push_back(val);
	 DigiHF_adc3.push_back((*j)[ii].adc());
       }
       if(ii==4){ 
	 DigiHF_charge4.push_back(val);
	 DigiHF_adc4.push_back((*j)[ii].adc());
       }
       if(ii==5){ 
	 DigiHF_charge5.push_back(val);
	 DigiHF_adc5.push_back((*j)[ii].adc());
       }
       if(ii==6){ 
	 DigiHF_charge6.push_back(val);
	 DigiHF_adc6.push_back((*j)[ii].adc());
       }
       if(ii==7){ 
	 DigiHF_charge7.push_back(val);
	 DigiHF_adc7.push_back((*j)[ii].adc());
       }
       if(ii==8){ 
	 DigiHF_charge8.push_back(val);
	 DigiHF_adc8.push_back((*j)[ii].adc());
       }
       if(ii==9){ 
	 DigiHF_charge9.push_back(val);
	 DigiHF_adc9.push_back((*j)[ii].adc());
       }

     }//Loop to get Charge

   }//Loop over HF Digis
   */

   //------------------------------------------------------------------------------------
   //HF QIE10Digis
   //------------------------------------------------------------------------------------
   edm::Handle< QIE10DigiCollection > QIE10digiTag;
   iEvent.getByToken(tok_qie10_hf_, QIE10digiTag);
   //reco<QIE10DataFrame>(iEvent, iSetup, tok_qie10_hf_);   
   
   for(QIE10DigiCollection::const_iterator j=QIE10digiTag->begin(); j != QIE10digiTag->end(); j++){  

     QIE10DataFrame dataFrame = *j;

     HcalDetId cell;
     cell = HcalDetId(j->id());
     int ieta = cell.ieta();
     int iphi = cell.iphi();
     int depth = cell.depth();
     int sub = cell.subdet();

     DigiHF_QIE10_ieta.push_back(ieta);
     DigiHF_QIE10_iphi.push_back(iphi);
     DigiHF_QIE10_depth.push_back(depth);
     DigiHF_QIE10_sub.push_back(sub);

     //std::cout << ieta << " " << iphi << " " << depth << " " << sub << std::endl;

     HcalCalibrations calibrations = conditions->getHcalCalibrations(cell);
     const HcalQIECoder* channelCoder = conditions->getHcalCoder(cell);
     const HcalQIEShape* shape = conditions->getHcalShape(channelCoder);
     HcalCoderDb coder(*channelCoder, *shape);
     coder.adc2fC(dataFrame, tool);

     for (int ii = 0; ii < tool.size(); ii++) {
       int capid = (dataFrame)[ii].capid();
       // single ts amplitude                                                                               
       double val = (tool[ii] - calibrations.pedestal(capid));

       //std::cout << "val: " << ii << ": " << val << std::endl;

       if(ii==0){ 
	 DigiHF_QIE10_charge0.push_back(val);      
	 DigiHF_QIE10_adc0.push_back((dataFrame)[ii].adc());
       }
       if(ii==1){
	 DigiHF_QIE10_charge1.push_back(val);
	 DigiHF_QIE10_adc1.push_back((dataFrame)[ii].adc());
       }
       if(ii==2){ 
	 DigiHF_QIE10_charge2.push_back(val);
	 DigiHF_QIE10_adc2.push_back((dataFrame)[ii].adc());
       }
       if(ii==3){ 
	 DigiHF_QIE10_charge3.push_back(val);
	 DigiHF_QIE10_adc3.push_back((dataFrame)[ii].adc());
       }
       if(ii==4){ 
	 DigiHF_QIE10_charge4.push_back(val);
	 DigiHF_QIE10_adc4.push_back((dataFrame)[ii].adc());
       }
       if(ii==5){ 
	 DigiHF_QIE10_charge5.push_back(val);
	 DigiHF_QIE10_adc5.push_back((dataFrame)[ii].adc());
       }
       if(ii==6){ 
	 DigiHF_QIE10_charge6.push_back(val);
	 DigiHF_QIE10_adc6.push_back((dataFrame)[ii].adc());
       }
       if(ii==7){ 
	 DigiHF_QIE10_charge7.push_back(val);
	 DigiHF_QIE10_adc7.push_back((dataFrame)[ii].adc());
       }
       if(ii==8){ 
	 DigiHF_QIE10_charge8.push_back(val);
	 DigiHF_QIE10_adc8.push_back((dataFrame)[ii].adc());
       }
       if(ii==9){ 
	 DigiHF_QIE10_charge9.push_back(val);
	 DigiHF_QIE10_adc9.push_back((dataFrame)[ii].adc());
       }

     }//Loop to get Charge

   }//Loop over HF_QIE10 Digis

   //Fill the tree
   tt1->Fill();

//std::cout << "Made it Here 4" << std::endl;

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


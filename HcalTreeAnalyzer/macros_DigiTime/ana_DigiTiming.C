// ------------------------------------------------------------------------------------
//  ROOT macro that produces average RecHit energy from PFG ntuples
//
//  Author : Ken H
//  Written on May 24, 2018
// ------------------------------------------------------------------------------------
//  
// Pre-requisite :
//
//   You should have the PFG ntuple for the Run from which you want to do a measurement. 
//   Instruction on how to make PFG ntuples can be found here : FIXME link here 
//
//   You should have "Fig" directory for plots 
//
// Usage : 
//
//   $ root -b  
//   root> .L ana_DigiTiming.C++ 
//   root> ana_DigiTiming("root://kodiak-se.baylor.edu//store/user/hatake/HCAL/ntuples/10_2_x/pi50_trees_MCfull.root","hcal_histograms_pt50.root",-1)
//    
// -----------------------------------------------------------------------------------
// 

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip> // for setw()
#include <algorithm> 

#include "TROOT.h"
#include "TF1.h"
#include "TMath.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TBranch.h"
#include "TString.h"
#include "TStyle.h"
#include "TInterpreter.h"
#include "TStyle.h"
#include "TLorentzVector.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// In order to use vector of vectors : vector<vector<data type> >
// ACLiC makes dictionary for this
// [ref] http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=10236&p=44117#p44117
#ifdef __MAKECINT__
#pragma link C++ class std::vector < std::vector<int> >+;
#pragma link C++ class std::vector < std::vector<float> >+;
#endif

using namespace std;

bool DRAWPLOTS  = false;  // draw plots or not (make "Fig" directory first before turning this on)
bool VERBOSE    = false;  // print out mean +/- sigma for each channel or not

//
// h2 cosmetics
//
void h2cosmetic(TH2F* &h2, char* title, TString Xvar, TString Yvar, TString Zvar);

//
// Det Name 
//
const char* GetDetName(int Subdet); 

//
// Book 1D histograms
//
void book1D(TList *v_hist, std::string name, int n, double min, double max);
void bookHistograms(TList *v_hist);

//
// Fill 1D histograms
//
void fill1D(TList *v_hist, std::string name, double value);

//
// Main analyzer
//
void HCALCheckRun(TString rootfile, TString outfile, int maxevents=-1, int option=2) 
{ 

   cout << "[Hcal analyzer] Running option " << option << " for " << endl; 

   // fit pannel display option
   gStyle->SetOptFit(1011);

   //
   // Get the tree from the PFG ntuple 
   //
   TChain *ch = new TChain("hcalTupleTree/tree");
   ch->Add(rootfile);
   printf("%d;\n",ch->GetNtrees());
   printf("%lld;\n",ch->GetEntries());

   TTreeReader     fReader(ch);  //!the tree reader

   //
   // Set up TTreeReader's
   // -- use MakeSelector of root
   //
   // Readers to access the data (delete the ones you do not need).
   TTreeReaderArray<double> GenParEta = {fReader, "GenParEta"};
   TTreeReaderArray<double> GenParM = {fReader, "GenParM"};
   TTreeReaderArray<double> GenParPhi = {fReader, "GenParPhi"};
   TTreeReaderArray<double> GenParPt = {fReader, "GenParPt"};
   TTreeReaderArray<vector<double>> QIE11DigiFC = {fReader, "QIE11DigiFC"};
   TTreeReaderArray<float> HBHEDigiEta = {fReader, "HBHEDigiEta"};
   TTreeReaderArray<float> HBHEDigiPhi = {fReader, "HBHEDigiPhi"};
   TTreeReaderArray<float> HBHEDigiRecEnergy = {fReader, "HBHEDigiRecEnergy"};
   TTreeReaderArray<float> HBHEDigiRecTime = {fReader, "HBHEDigiRecTime"};
   TTreeReaderArray<float> HBHERecHitEnergy = {fReader, "HBHERecHitEnergy"};
   TTreeReaderArray<float> HBHERecHitEta = {fReader, "HBHERecHitEta"};
   TTreeReaderArray<float> HBHERecHitPhi = {fReader, "HBHERecHitPhi"};
   TTreeReaderArray<float> HBHERecHitTime = {fReader, "HBHERecHitTime"};
   TTreeReaderArray<float> HFDigiEta = {fReader, "HFDigiEta"};
   TTreeReaderArray<float> HFDigiPhi = {fReader, "HFDigiPhi"};
   TTreeReaderArray<float> HFDigiRecEnergy = {fReader, "HFDigiRecEnergy"};
   TTreeReaderArray<float> HFDigiRecTime = {fReader, "HFDigiRecTime"};
   TTreeReaderArray<float> HODigiEta = {fReader, "HODigiEta"};
   TTreeReaderArray<float> HODigiPhi = {fReader, "HODigiPhi"};
   TTreeReaderArray<float> HODigiRecEnergy = {fReader, "HODigiRecEnergy"};
   TTreeReaderArray<float> HODigiRecTime = {fReader, "HODigiRecTime"};
   TTreeReaderArray<float> HcalSimHitsEnergy = {fReader, "HcalSimHitsEnergy"};
   TTreeReaderArray<float> HcalSimHitsEta = {fReader, "HcalSimHitsEta"};
   TTreeReaderArray<float> HcalSimHitsPhi = {fReader, "HcalSimHitsPhi"};
   TTreeReaderArray<float> HcalSimHitsPosx = {fReader, "HcalSimHitsPosx"};
   TTreeReaderArray<float> HcalSimHitsPosy = {fReader, "HcalSimHitsPosy"};
   TTreeReaderArray<float> HcalSimHitsPosz = {fReader, "HcalSimHitsPosz"};
   TTreeReaderArray<float> HcalSimHitsTime = {fReader, "HcalSimHitsTime"};
   TTreeReaderArray<float> HcalSimHitsTimeTOF = {fReader, "HcalSimHitsTimeTOF"};
   TTreeReaderArray<vector<float>> HBHEDigiAllFC = {fReader, "HBHEDigiAllFC"};
   TTreeReaderArray<vector<float>> HBHEDigiEnergy = {fReader, "HBHEDigiEnergy"};
   TTreeReaderArray<vector<float>> HBHEDigiFC = {fReader, "HBHEDigiFC"};
   TTreeReaderArray<vector<float>> HBHEDigiGain = {fReader, "HBHEDigiGain"};
   TTreeReaderArray<vector<float>> HBHEDigiNomFC = {fReader, "HBHEDigiNomFC"};
   TTreeReaderArray<vector<float>> HBHEDigiPedFC = {fReader, "HBHEDigiPedFC"};
   TTreeReaderArray<vector<float>> HBHEDigiRCGain = {fReader, "HBHEDigiRCGain"};
   TTreeReaderArray<vector<float>> HFDigiAllFC = {fReader, "HFDigiAllFC"};
   TTreeReaderArray<vector<float>> HFDigiEnergy = {fReader, "HFDigiEnergy"};
   TTreeReaderArray<vector<float>> HFDigiFC = {fReader, "HFDigiFC"};
   TTreeReaderArray<vector<float>> HFDigiGain = {fReader, "HFDigiGain"};
   TTreeReaderArray<vector<float>> HFDigiNomFC = {fReader, "HFDigiNomFC"};
   TTreeReaderArray<vector<float>> HFDigiPedFC = {fReader, "HFDigiPedFC"};
   TTreeReaderArray<vector<float>> HFDigiRCGain = {fReader, "HFDigiRCGain"};
   TTreeReaderArray<vector<float>> HODigiAllFC = {fReader, "HODigiAllFC"};
   TTreeReaderArray<vector<float>> HODigiEnergy = {fReader, "HODigiEnergy"};
   TTreeReaderArray<vector<float>> HODigiFC = {fReader, "HODigiFC"};
   TTreeReaderArray<vector<float>> HODigiGain = {fReader, "HODigiGain"};
   TTreeReaderArray<vector<float>> HODigiNomFC = {fReader, "HODigiNomFC"};
   TTreeReaderArray<vector<float>> HODigiPedFC = {fReader, "HODigiPedFC"};
   TTreeReaderArray<vector<float>> HODigiRCGain = {fReader, "HODigiRCGain"};
   TTreeReaderArray<int> GenParPdgId = {fReader, "GenParPdgId"};
   TTreeReaderArray<int> GenParStatus = {fReader, "GenParStatus"};
   TTreeReaderArray<int> HBHEDigiDepth = {fReader, "HBHEDigiDepth"};
   TTreeReaderArray<int> HBHEDigiElectronicsID = {fReader, "HBHEDigiElectronicsID"};
   TTreeReaderArray<int> HBHEDigiFiberIdleOffset = {fReader, "HBHEDigiFiberIdleOffset"};
   TTreeReaderArray<int> HBHEDigiIEta = {fReader, "HBHEDigiIEta"};
   TTreeReaderArray<int> HBHEDigiIPhi = {fReader, "HBHEDigiIPhi"};
   TTreeReaderArray<int> HBHEDigiPresamples = {fReader, "HBHEDigiPresamples"};
   TTreeReaderArray<int> HBHEDigiRawID = {fReader, "HBHEDigiRawID"};
   TTreeReaderArray<int> HBHEDigiSize = {fReader, "HBHEDigiSize"};
   TTreeReaderArray<int> HBHEDigiSubdet = {fReader, "HBHEDigiSubdet"};
   TTreeReaderArray<int> HBHERecHitAux = {fReader, "HBHERecHitAux"};
   TTreeReaderArray<int> HBHERecHitDepth = {fReader, "HBHERecHitDepth"};
   TTreeReaderArray<int> HBHERecHitFlags = {fReader, "HBHERecHitFlags"};
   TTreeReaderArray<int> HBHERecHitHPDid = {fReader, "HBHERecHitHPDid"};
   TTreeReaderArray<int> HBHERecHitIEta = {fReader, "HBHERecHitIEta"};
   TTreeReaderArray<int> HBHERecHitIPhi = {fReader, "HBHERecHitIPhi"};
   TTreeReaderArray<int> HBHERecHitRBXid = {fReader, "HBHERecHitRBXid"};
   TTreeReaderArray<int> HBHEDigiSOI = {fReader, "HBHEDigiSOI"};
   TTreeReaderArray<int> HFDigiDepth = {fReader, "HFDigiDepth"};
   TTreeReaderArray<int> HFDigiElectronicsID = {fReader, "HFDigiElectronicsID"};
   TTreeReaderArray<int> HFDigiFiberIdleOffset = {fReader, "HFDigiFiberIdleOffset"};
   TTreeReaderArray<int> HFDigiIEta = {fReader, "HFDigiIEta"};
   TTreeReaderArray<int> HFDigiIPhi = {fReader, "HFDigiIPhi"};
   TTreeReaderArray<int> HFDigiPresamples = {fReader, "HFDigiPresamples"};
   TTreeReaderArray<int> HFDigiRawID = {fReader, "HFDigiRawID"};
   TTreeReaderArray<int> HFDigiSize = {fReader, "HFDigiSize"};
   TTreeReaderArray<int> HFDigiSubdet = {fReader, "HFDigiSubdet"};
   TTreeReaderArray<int> HODigiDepth = {fReader, "HODigiDepth"};
   TTreeReaderArray<int> HODigiElectronicsID = {fReader, "HODigiElectronicsID"};
   TTreeReaderArray<int> HODigiFiberIdleOffset = {fReader, "HODigiFiberIdleOffset"};
   TTreeReaderArray<int> HODigiIEta = {fReader, "HODigiIEta"};
   TTreeReaderArray<int> HODigiIPhi = {fReader, "HODigiIPhi"};
   TTreeReaderArray<int> HODigiPresamples = {fReader, "HODigiPresamples"};
   TTreeReaderArray<int> HODigiRawID = {fReader, "HODigiRawID"};
   TTreeReaderArray<int> HODigiSize = {fReader, "HODigiSize"};
   TTreeReaderArray<int> HODigiSubdet = {fReader, "HODigiSubdet"};
   TTreeReaderArray<int> HcalSimHitsDepth = {fReader, "HcalSimHitsDepth"};
   TTreeReaderArray<int> HcalSimHitsIeta = {fReader, "HcalSimHitsIeta"};
   TTreeReaderArray<int> HcalSimHitsIndex = {fReader, "HcalSimHitsIndex"};
   TTreeReaderArray<int> HcalSimHitsIphi = {fReader, "HcalSimHitsIphi"};
   TTreeReaderArray<int> HcalSimHitsSubdet = {fReader, "HcalSimHitsSubdet"};
   TTreeReaderArray<int> QIE11DigiCapIDError = {fReader, "QIE11DigiCapIDError"};
   TTreeReaderArray<int> QIE11DigiDepth = {fReader, "QIE11DigiDepth"};
   TTreeReaderArray<int> QIE11DigiFlags = {fReader, "QIE11DigiFlags"};
   TTreeReaderArray<int> QIE11DigiIEta = {fReader, "QIE11DigiIEta"};
   TTreeReaderArray<int> QIE11DigiIPhi = {fReader, "QIE11DigiIPhi"};
   TTreeReaderArray<int> QIE11DigiLinkError = {fReader, "QIE11DigiLinkError"};
   TTreeReaderArray<int> QIE11DigiRawID = {fReader, "QIE11DigiRawID"};
   TTreeReaderArray<int> QIE11DigiSubdet = {fReader, "QIE11DigiSubdet"};
   TTreeReaderArray<int> QIE11DigiSOI = {fReader, "QIE11DigiSOI"};   
   TTreeReaderArray<vector<int>> HBHEDigiADC = {fReader, "HBHEDigiADC"};
   TTreeReaderArray<vector<int>> HBHEDigiCapID = {fReader, "HBHEDigiCapID"};
   TTreeReaderArray<vector<int>> HBHEDigiDV = {fReader, "HBHEDigiDV"};
   TTreeReaderArray<vector<int>> HBHEDigiER = {fReader, "HBHEDigiER"};
   TTreeReaderArray<vector<int>> HBHEDigiFiber = {fReader, "HBHEDigiFiber"};
   TTreeReaderArray<vector<int>> HBHEDigiFiberChan = {fReader, "HBHEDigiFiberChan"};
   TTreeReaderArray<vector<int>> HBHEDigiLADC = {fReader, "HBHEDigiLADC"};
   TTreeReaderArray<vector<int>> HBHEDigiRaw = {fReader, "HBHEDigiRaw"};
   TTreeReaderArray<vector<int>> HFDigiADC = {fReader, "HFDigiADC"};
   TTreeReaderArray<vector<int>> HFDigiCapID = {fReader, "HFDigiCapID"};
   TTreeReaderArray<vector<int>> HFDigiDV = {fReader, "HFDigiDV"};
   TTreeReaderArray<vector<int>> HFDigiER = {fReader, "HFDigiER"};
   TTreeReaderArray<vector<int>> HFDigiFiber = {fReader, "HFDigiFiber"};
   TTreeReaderArray<vector<int>> HFDigiFiberChan = {fReader, "HFDigiFiberChan"};
   TTreeReaderArray<vector<int>> HFDigiLADC = {fReader, "HFDigiLADC"};
   TTreeReaderArray<vector<int>> HFDigiRaw = {fReader, "HFDigiRaw"};
   TTreeReaderArray<vector<int>> HODigiADC = {fReader, "HODigiADC"};
   TTreeReaderArray<vector<int>> HODigiCapID = {fReader, "HODigiCapID"};
   TTreeReaderArray<vector<int>> HODigiDV = {fReader, "HODigiDV"};
   TTreeReaderArray<vector<int>> HODigiER = {fReader, "HODigiER"};
   TTreeReaderArray<vector<int>> HODigiFiber = {fReader, "HODigiFiber"};
   TTreeReaderArray<vector<int>> HODigiFiberChan = {fReader, "HODigiFiberChan"};
   TTreeReaderArray<vector<int>> HODigiLADC = {fReader, "HODigiLADC"};
   TTreeReaderArray<vector<int>> HODigiRaw = {fReader, "HODigiRaw"};
   TTreeReaderArray<vector<int>> QIE11DigiADC = {fReader, "QIE11DigiADC"};
   TTreeReaderArray<vector<int>> QIE11DigiCapID = {fReader, "QIE11DigiCapID"};
   TTreeReaderArray<vector<int>> QIE11DigiTDC = {fReader, "QIE11DigiTDC"};
   TTreeReaderValue<UInt_t> bx = {fReader, "bx"};
   TTreeReaderValue<UInt_t> event = {fReader, "event"};
   TTreeReaderValue<UInt_t> ls = {fReader, "ls"};
   TTreeReaderValue<UInt_t> orbit = {fReader, "orbit"};
   TTreeReaderValue<UInt_t> run = {fReader, "run"};

   //
   // Define histograms to fill
   //
   TList *v_hist = new TList();
   
   TH1F *h_RecHitEtGenPt = new TH1F("h_RecHitEtGenPt","h_RecHitEtGenPt",100,0.,2.);
   v_hist->Add(h_RecHitEtGenPt);
   v_hist->FindObject("h_RecHitEtGenPt")->Print();

   bookHistograms(v_hist); // most of histograms booked here
   
   //
   // Loop over entries
   //
   unsigned int nentries = (Int_t)ch->GetEntries();
   cout << "[Hcal analyzer] The number of entries is: " << nentries << endl;

   //---------------------------------------------------------------------------------------------------------
   // main event loop
   //---------------------------------------------------------------------------------------------------------

   int ievent=0;
   while (fReader.Next()) {
  
     // Progress indicator 
     ievent++;
     if(ievent%100==0) cout << "[HCAL analyzer] Processed " << ievent << " out of " << nentries << " events" << endl; 
     if (maxevents>0 && ievent>maxevents) break;
     
     //std::cout << *event << std::endl;
	
     //--------------------
     // Loop over pions
     //------------------
     for (int iGenPar = 0, nGenPar =  GenParPt.GetSize(); iGenPar < nGenPar; ++iGenPar) {
       //std::cout << GenParPdgId[iGenPar] << std::endl;
       TLorentzVector TLVPion; TLVPion.SetPtEtaPhiM(GenParPt[iGenPar],GenParEta[iGenPar],GenParPhi[iGenPar],GenParM[iGenPar]);
       //TLVPion.Print();
       //if (fabs(TLVPion.Eta())<1.8 || fabs(TLVPion.Eta())>2.4) continue;
       
       // Loop over HBHERecHits
       double SumEt=0.;
       for (int irc = 0, nrc =  HBHERecHitEnergy.GetSize(); irc < nrc; ++irc) {
	 TLorentzVector TLVRecHit; 
	 double RecHitPt=HBHERecHitEnergy[irc]/cosh(HBHERecHitEta[irc]); //p=E for rechits
	 TLVRecHit.SetPtEtaPhiE(RecHitPt,HBHERecHitEta[irc],HBHERecHitPhi[irc],HBHERecHitEnergy[irc]);	  
	 double dR=TLVRecHit.DeltaR(TLVPion);
	 if (dR<0.3) SumEt+=TLVRecHit.Pt();	  
       }
       //std::cout << SumEt/TLVPion.Pt() <<std::endl;
       h_RecHitEtGenPt->Fill( SumEt/TLVPion.Pt());
     } // Loop over pions ends	

     //--------------------
     // Loop over digis 
     //--------------------
     for (int idigi = 0, ndigi =  HBHEDigiFC.GetSize(); idigi < ndigi; ++idigi) {
       int SOI=HBHEDigiSOI[idigi];

       //
       // Define SOI charge fraction
       //
       double v_ampl=0.;
       double fbinSOI = HBHEDigiFC[idigi][SOI];
       for (int iTS = SOI, nTS = HBHEDigiFC[idigi].size(); iTS < nTS; iTS++) v_ampl += HBHEDigiFC[idigi][iTS];
       double fbinPS = v_ampl - fbinSOI;
       if (v_ampl>0.){
	 fbinSOI /= v_ampl;
	 fbinPS /= v_ampl;
       }

       std::string strtmp;
       std::string subdet_="HB";
       
       // 
       if (v_ampl>30.){ // Amplitude selection

	 strtmp = "HcalDigiTask_SOI_frac_" + subdet_;
	 fill1D(v_hist, strtmp, fbinSOI);
	 strtmp = "HcalDigiTask_postSOI_frac_" + subdet_;
	 fill1D(v_hist, strtmp, fbinPS);
	 
       } // v_ampl threshold

     }

     //--------------------
     // Loop over digis 
     //--------------------
     for (int idigi = 0, ndigi =  QIE11DigiFC.GetSize(); idigi < ndigi; ++idigi) {
       int SOI=QIE11DigiSOI[idigi];

       int ieta=QIE11DigiIEta[idigi];
       int iphi=QIE11DigiIPhi[idigi];
       int depth=QIE11DigiDepth[idigi];
       
       //
       // Define SOI charge fraction
       //
       double v_ampl=0.;
       double fbinSOI = QIE11DigiFC[idigi][SOI];
       for (int iTS = SOI, nTS = QIE11DigiFC[idigi].size(); iTS < nTS; iTS++) v_ampl += QIE11DigiFC[idigi][iTS];
       double fbinPS = v_ampl - fbinSOI;
       if (v_ampl>0.){
	 fbinSOI /= v_ampl;
	 fbinPS /= v_ampl;
       }

       bool goodtest = true; // define this much ealier so that we can use it at different locations
       double chargeSOI = QIE11DigiFC[idigi][SOI];
       for (int iTS = SOI, nTS = QIE11DigiFC[idigi].size(); iTS < nTS; iTS++){
	 if (QIE11DigiFC[idigi][iTS]>chargeSOI) goodtest = false;
       }

       std::string strtmp;
       std::string subdet_="HE";
       
       // 
       if (v_ampl>300.){ // Amplitude selection

	 strtmp = "HcalDigiTask_SOI_frac_" + subdet_;
	 fill1D(v_hist, strtmp, fbinSOI);
	 strtmp = "HcalDigiTask_postSOI_frac_" + subdet_;
	 fill1D(v_hist, strtmp, fbinPS);
	 // 
	 //KH starts
	 if (v_ampl > 1000.) {
	   strtmp = "HcalDigiTask_SOI_frac_1000_" + subdet_;
	   fill1D(v_hist, strtmp, fbinSOI);
	   strtmp = "HcalDigiTask_postSOI_frac_1000_" + subdet_;
	   fill1D(v_hist, strtmp, fbinPS);
	 }
	 // 
	 if (v_ampl > 2000.) {
	   strtmp = "HcalDigiTask_SOI_frac_2000_" + subdet_;
	   fill1D(v_hist, strtmp, fbinSOI);
	   strtmp = "HcalDigiTask_postSOI_frac_2000_" + subdet_;
	   fill1D(v_hist, strtmp, fbinPS);
	 }
	 // 
	 if (v_ampl > 6000.) {
	   strtmp = "HcalDigiTask_SOI_frac_6000_" + subdet_;
	   fill1D(v_hist, strtmp, fbinSOI);
	   strtmp = "HcalDigiTask_postSOI_frac_6000_" + subdet_;
	   fill1D(v_hist, strtmp, fbinPS);
	 }

	 // Another test
	 // goodtest=true; if (fbinSOI<0.1) goodtest=false;
	 if (goodtest) fill1D(v_hist,"HcalDigiTask_Charge_Prompt_HE", v_ampl);  // total charge for SOI-lastTS
	 else          fill1D(v_hist,"HcalDigiTask_Charge_Delayed_HE", v_ampl); // total charge for SOI-lastTS
	 
	 double aveSimTime=0.;
	 double sumSimHitE=0.;
	 
	 for (int isim = 0, nsim =  HcalSimHitsEnergy.GetSize(); isim < nsim; ++isim) {
	   
	   int ieta_sim=HcalSimHitsIeta[isim];
	   int iphi_sim=HcalSimHitsIphi[isim];
	   int depth_sim=HcalSimHitsDepth[isim];
	   double ttime=HcalSimHitsTime[isim];
	   double ten=HcalSimHitsEnergy[isim];
	   
	   if (ieta==ieta_sim && iphi==iphi_sim && depth==depth_sim){
	   if (ttime<135.){  // simhits beyond 135 ns won't contribute to digis
	     aveSimTime += ten*ttime; // KH
	     sumSimHitE += ten; // KH	     
	   }	            
	   }

	 } // loop over simhits

	 aveSimTime /= sumSimHitE; //KH
	 if (goodtest) fill1D(v_hist,"Simhit_AveTime_PromptHits_HE", aveSimTime);
	 else          fill1D(v_hist,"Simhit_AveTime_DelayedHits_HE", aveSimTime);
	 
       } // v_ampl threshold
       
     } // Loop over HE digis    
     
   }   // Event loop ends
   //---------------------------------------------------------------------------------------------------------
   // main event loop ends
   //---------------------------------------------------------------------------------------------------------

   // output file for histograms
   TFile file_out(outfile,"RECREATE");
   
   h_RecHitEtGenPt->Fit("gaus");
   //h_RecHitEtGenPt->Write();
   v_hist->Write();
   
   file_out.ls();
   file_out.Close();

}

//
// Main function
//
void ana_DigiTiming(TString rootfile="/cms/data/store/user/hatake/HCAL/ntuples/10_2_x/pi50_trees_MCfull_CMSSW_10_2_0_pre3_*.root",TString outfile="hcal_timestudy_histograms.root",int maxevents=-1)
{
  HCALCheckRun(rootfile, outfile, maxevents, 0);
}

//
// --- Aux ---
//

//
// H2 cosmetics
//
void h2cosmetic(TH2F* &h2, char* title, TString Xvar="", TString Yvar="", TString Zvar="Events/bin")
{
    h2->SetTitle(title);
    h2->SetXTitle(Xvar);
    h2->SetYTitle(Yvar);
    h2->SetZTitle(Zvar);
    h2->SetStats(0);
}

//
// Det Name 
//
const char* GetDetName(int Subdet); 
const char* GetDetName(int Subdet) 
{ 
    const char* DetName="";
    if(Subdet==1) DetName = "HB"; 
    if(Subdet==2) DetName = "HE"; 
    if(Subdet==3) DetName = "HO"; 
    if(Subdet==4) DetName = "HF"; 
    return DetName;
}
//
// Book 1D histograms
//
void book1D(TList *v_hist, std::string name, int n, double min, double max)
{
  TH1D *htemp = new TH1D(name.c_str(), name.c_str(), n, min, max);
  v_hist->Add(htemp);
}
//
// Book histograms
//
void bookHistograms(TList *v_hist)
{

  Char_t histo[100];

  //
  // Booking histograms
  // 
  for (int idet=0; idet<=1; idet++){
    const char * sub = GetDetName(idet+1);
    
    //KH --
    sprintf(histo, "HcalDigiTask_SOI_frac_%s", sub);
    book1D(v_hist, histo, 80, -0.20, 1.40);
    sprintf(histo, "HcalDigiTask_postSOI_frac_%s", sub);
    book1D(v_hist, histo, 80, -0.20, 1.40);

    //-----
    if (idet==1){ 
    //KH --
    sprintf(histo, "HcalDigiTask_SOI_frac_pass_%s", sub);
    book1D(v_hist, histo, 80, -0.20, 1.40);
    sprintf(histo, "HcalDigiTask_postSOI_frac_pass_%s", sub);
    book1D(v_hist, histo, 80, -0.20, 1.40);
    sprintf(histo, "HcalDigiTask_SOI_frac_fail_%s", sub);
    book1D(v_hist, histo, 80, -0.20, 1.40);
    sprintf(histo, "HcalDigiTask_postSOI_frac_fail_%s", sub);
    book1D(v_hist, histo, 80, -0.20, 1.40);
    //KH --
    sprintf(histo, "HcalDigiTask_SOI_frac_1000_%s", sub);
    book1D(v_hist, histo, 80, -0.20, 1.40);
    sprintf(histo, "HcalDigiTask_SOI_frac_2000_%s", sub);
    book1D(v_hist, histo, 80, -0.20, 1.40);
    sprintf(histo, "HcalDigiTask_SOI_frac_6000_%s", sub);
    book1D(v_hist, histo, 80, -0.20, 1.40);
    //KH --
    sprintf(histo, "HcalDigiTask_postSOI_frac_1000_%s", sub);
    book1D(v_hist, histo, 80, -0.20, 1.40);
    sprintf(histo, "HcalDigiTask_postSOI_frac_2000_%s", sub);
    book1D(v_hist, histo, 80, -0.20, 1.40);
    sprintf(histo, "HcalDigiTask_postSOI_frac_6000_%s", sub);
    book1D(v_hist, histo, 80, -0.20, 1.40);
    } // HE only
  } // idet loop

  sprintf(histo, "HcalDigiTask_Charge_Prompt_HE");
  book1D(v_hist, histo, 200,0.,20000.);
  sprintf(histo, "HcalDigiTask_Charge_Delayed_HE");
  book1D(v_hist, histo, 200,0.,20000.);

  sprintf(histo, "Simhit_AveTime_PromptHits_HE");
  book1D(v_hist, histo, 200,1.,200.); // TEST plot - simhit time (simhit E weighted average) for prompt digis 
  sprintf(histo, "Simhit_AveTime_DelayedHits_HE");
  book1D(v_hist, histo, 200,1.,200.); // TEST plot - simhit time (simhit E weighted average) for delayed digis
    
}
//
// Fill 1D histograms
//
void fill1D(TList *v_hist, std::string name, double value)
{
  TH1F* htemp = (TH1F*) v_hist->FindObject(name.c_str());
  htemp->Fill(value);
}

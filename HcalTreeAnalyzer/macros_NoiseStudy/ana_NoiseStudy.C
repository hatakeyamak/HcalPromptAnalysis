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
//   root> .L ana_NoiseStudy.C+
//   root> ana_NoiseStudy("/cms/data/store/user/hatake/HCAL/ntuples/10_2_x/pi50_trees_MCfull_CMSSW_10_2_0_pre3_*.root","hcal_timestudy_pi50_histograms.root")
//   or
//   root> ana_NoiseStudy("list_trees_pi50_MCfull_CMSSW_10_2_0_pre3.txt","hcal_timestudy_pi50_histograms.root")
//   or
//   from command line:
/*
     root.exe -b -q 'ana_NoiseStudy.C++("trees_relval_ttbar_phase2_age_new2_4500ultimate.root","hcal_noisestudy_histograms_age_new2_4500ultimate.root")'
     root.exe -b -q 'ana_NoiseStudy.C++("trees_relval_ttbar_phase2_age_org.root","hcal_noisestudy_histograms_age_org.root")'
     root.exe -b -q 'ana_NoiseStudy.C++("trees_relval_ttbar_phase2_noage.root","hcal_noisestudy_histograms_noage.root")'
 */
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

   std::string filename(rootfile);
   std::string::size_type idx;
   idx = filename.rfind('.');
   std::string extension = filename.substr(idx+1);
   std::string line;
   
   if(idx != std::string::npos && extension=="txt")
     {
       std::cout << rootfile << " " << extension << std::endl;
       std::ifstream in(rootfile);
       while (std::getline(in, line)) {     // Process line
	 if (line.size()>0) ch->Add(line.c_str());
       }
     }
   else
     {
       // No extension found
       ch->Add(rootfile);
     }

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
   TTreeReaderArray<vector<double>> QIE11DigiRawFC = {fReader, "QIE11DigiRawFC"};
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
   TTreeReaderValue<Int_t> laserType = {fReader, "laserType"};
   TTreeReaderArray<int> GenParPdgId = {fReader, "GenParPdgId"};
   TTreeReaderArray<int> GenParStatus = {fReader, "GenParStatus"};
   TTreeReaderArray<int> HBHEDigiDepth = {fReader, "HBHEDigiDepth"};
   TTreeReaderArray<int> HBHEDigiElectronicsID = {fReader, "HBHEDigiElectronicsID"};
   TTreeReaderArray<int> HBHEDigiFiberIdleOffset = {fReader, "HBHEDigiFiberIdleOffset"};
   TTreeReaderArray<int> HBHEDigiIEta = {fReader, "HBHEDigiIEta"};
   TTreeReaderArray<int> HBHEDigiIPhi = {fReader, "HBHEDigiIPhi"};
   TTreeReaderArray<int> HBHEDigiPresamples = {fReader, "HBHEDigiPresamples"};
   TTreeReaderArray<int> HBHEDigiRawID = {fReader, "HBHEDigiRawID"};
   TTreeReaderArray<int> HBHEDigiSOI = {fReader, "HBHEDigiSOI"};
   TTreeReaderArray<int> HBHEDigiSize = {fReader, "HBHEDigiSize"};
   TTreeReaderArray<int> HBHEDigiSubdet = {fReader, "HBHEDigiSubdet"};
   TTreeReaderArray<int> HBHERecHitAux = {fReader, "HBHERecHitAux"};
   TTreeReaderArray<int> HBHERecHitDepth = {fReader, "HBHERecHitDepth"};
   TTreeReaderArray<int> HBHERecHitFlags = {fReader, "HBHERecHitFlags"};
   TTreeReaderArray<int> HBHERecHitHPDid = {fReader, "HBHERecHitHPDid"};
   TTreeReaderArray<int> HBHERecHitIEta = {fReader, "HBHERecHitIEta"};
   TTreeReaderArray<int> HBHERecHitIPhi = {fReader, "HBHERecHitIPhi"};
   TTreeReaderArray<int> HBHERecHitRBXid = {fReader, "HBHERecHitRBXid"};
   TTreeReaderArray<int> HFDigiDepth = {fReader, "HFDigiDepth"};
   TTreeReaderArray<int> HFDigiElectronicsID = {fReader, "HFDigiElectronicsID"};
   TTreeReaderArray<int> HFDigiFiberIdleOffset = {fReader, "HFDigiFiberIdleOffset"};
   TTreeReaderArray<int> HFDigiIEta = {fReader, "HFDigiIEta"};
   TTreeReaderArray<int> HFDigiIPhi = {fReader, "HFDigiIPhi"};
   TTreeReaderArray<int> HFDigiPresamples = {fReader, "HFDigiPresamples"};
   TTreeReaderArray<int> HFDigiRawID = {fReader, "HFDigiRawID"};
   TTreeReaderArray<int> HFDigiSOI = {fReader, "HFDigiSOI"};
   TTreeReaderArray<int> HFDigiSize = {fReader, "HFDigiSize"};
   TTreeReaderArray<int> HFDigiSubdet = {fReader, "HFDigiSubdet"};
   TTreeReaderArray<int> HODigiDepth = {fReader, "HODigiDepth"};
   TTreeReaderArray<int> HODigiElectronicsID = {fReader, "HODigiElectronicsID"};
   TTreeReaderArray<int> HODigiFiberIdleOffset = {fReader, "HODigiFiberIdleOffset"};
   TTreeReaderArray<int> HODigiIEta = {fReader, "HODigiIEta"};
   TTreeReaderArray<int> HODigiIPhi = {fReader, "HODigiIPhi"};
   TTreeReaderArray<int> HODigiPresamples = {fReader, "HODigiPresamples"};
   TTreeReaderArray<int> HODigiRawID = {fReader, "HODigiRawID"};
   TTreeReaderArray<int> HODigiSOI = {fReader, "HODigiSOI"};
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
   TTreeReaderArray<int> QIE11DigiSOI = {fReader, "QIE11DigiSOI"};
   TTreeReaderArray<int> QIE11DigiSubdet = {fReader, "QIE11DigiSubdet"};
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
     if(ievent%10==0) cout << "[HCAL analyzer] Processed " << ievent << " out of " << nentries << " events" << endl; 
     if (maxevents>0 && ievent>maxevents) break;
     
     //--------------------
     // Loop over rechits
     //--------------------
     for (int irechit = 0, nrechit =  HBHERecHitEnergy.GetSize(); irechit < nrechit; ++irechit) {
       
       std::string strtmp;
       std::string subdet_="HB";

       double rcenergy = HBHERecHitEnergy[irechit];
       int depth=HBHERecHitDepth[irechit];
       
       strtmp = "RecHitEnergy_depth" + std::to_string(depth) + "_" + subdet_;
       fill1D(v_hist, strtmp, rcenergy);

     }

     //--------------------
     // Loop over digis 
     //--------------------
     for (int idigi = 0, ndigi =  QIE11DigiFC.GetSize(); idigi < ndigi; ++idigi) {
       int SOI=QIE11DigiSOI[idigi];

       int ieta=QIE11DigiIEta[idigi];
       int iphi=QIE11DigiIPhi[idigi];
       int depth=QIE11DigiDepth[idigi];
       
       std::string strtmp;
       std::string subdet_="HB";

       //
       // Define SOI charge fraction
       //
       double v_ampl=0.;
       double v_ampl2=0.;
       for (int iTS = SOI, nTS = QIE11DigiFC[idigi].size(); iTS < min(nTS,SOI+2); iTS++){
	 v_ampl  += QIE11DigiFC[idigi][iTS];
	 v_ampl2 += QIE11DigiRawFC[idigi][iTS];
	 //std::cout << v_ampl << std::endl;
       }
	 
       strtmp = "HcalDigiTask_FC_depth" + std::to_string(depth) + "_" + subdet_;
       fill1D(v_hist, strtmp, v_ampl);
       strtmp = "HcalDigiTask_RawFC_depth" + std::to_string(depth) + "_" + subdet_;
       fill1D(v_hist, strtmp, v_ampl2);
       //std::cout << v_ampl << std::endl;

       if (depth<=2 || (depth==3 && abs(ieta)==16) ){
	 strtmp = "HcalDigiTask_FC_2p8mmSiPM_" + subdet_;
	 fill1D(v_hist,strtmp, v_ampl);
	 strtmp = "HcalDigiTask_RawFC_2p8mmSiPM_" + subdet_;
	 fill1D(v_hist,strtmp, v_ampl2);
       }
       else {
	 strtmp = "HcalDigiTask_FC_3p3mmSiPM_" + subdet_;
	 fill1D(v_hist,strtmp, v_ampl);
	 strtmp = "HcalDigiTask_RawFC_3p3mmSiPM_" + subdet_;
	 fill1D(v_hist,strtmp, v_ampl2);
       }
        
     } // Loop over HB digis    
     
   }   // Event loop ends
   //---------------------------------------------------------------------------------------------------------
   // main event loop ends
   //---------------------------------------------------------------------------------------------------------

   // output file for histograms
   TFile file_out(outfile,"RECREATE");
   
   v_hist->Write();
   
   file_out.ls();
   file_out.Close();

}

//
// Main function
//
void ana_NoiseStudy(TString rootfile="trees_relval_ttbar_phase2_age_new.root",TString outfile="hcal_aging_histograms.root",int maxevents=-1)
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
  std::string strtmp;
  
  //
  // Booking histograms
  // 
  for (int idet=0; idet<=1; idet++){
    const char * sub = GetDetName(idet+1);
    
    //KH --
    sprintf(histo, "HcalDigiTask_FC_2p8mmSiPM_%s", sub);
    book1D(v_hist, histo, 275, -10000., 100000.);
    sprintf(histo, "HcalDigiTask_FC_3p3mmSiPM_%s", sub);
    book1D(v_hist, histo, 275, -10000., 100000.);
    //
    sprintf(histo, "HcalDigiTask_FC_depth1_%s", sub);
    book1D(v_hist, histo, 250, 0., 100000.);
    sprintf(histo, "HcalDigiTask_FC_depth2_%s", sub);
    book1D(v_hist, histo, 250, 0., 100000.);
    sprintf(histo, "HcalDigiTask_FC_depth3_%s", sub);
    book1D(v_hist, histo, 250, 0., 100000.);
    sprintf(histo, "HcalDigiTask_FC_depth4_%s", sub);
    book1D(v_hist, histo, 250, 0., 100000.);
    //KH --
    sprintf(histo, "HcalDigiTask_RawFC_2p8mmSiPM_%s", sub);
    book1D(v_hist, histo, 275, -10000., 100000.);
    sprintf(histo, "HcalDigiTask_RawFC_3p3mmSiPM_%s", sub);
    book1D(v_hist, histo, 275, -10000., 100000.);
    //
    sprintf(histo, "HcalDigiTask_RawFC_depth1_%s", sub);
    book1D(v_hist, histo, 250, 0., 100000.);
    sprintf(histo, "HcalDigiTask_RawFC_depth2_%s", sub);
    book1D(v_hist, histo, 250, 0., 100000.);
    sprintf(histo, "HcalDigiTask_RawFC_depth3_%s", sub);
    book1D(v_hist, histo, 250, 0., 100000.);
    sprintf(histo, "HcalDigiTask_RawFC_depth4_%s", sub);
    book1D(v_hist, histo, 250, 0., 100000.);
    
    //KH
    sprintf(histo, "RecHitEnergy_depth1_%s", sub);
    book1D(v_hist, histo, 50, 0., 5.);
    sprintf(histo, "RecHitEnergy_depth2_%s", sub);
    book1D(v_hist, histo, 50, 0., 5.);
    sprintf(histo, "RecHitEnergy_depth3_%s", sub);
    book1D(v_hist, histo, 50, 0., 5.);
    sprintf(histo, "RecHitEnergy_depth4_%s", sub);
    book1D(v_hist, histo, 50, 0., 5.);
    //
    
  }
    
}
//
// Fill 1D histograms
//
void fill1D(TList *v_hist, std::string name, double value)
{
  TH1F* htemp = (TH1F*) v_hist->FindObject(name.c_str());
  htemp->Fill(value);
}

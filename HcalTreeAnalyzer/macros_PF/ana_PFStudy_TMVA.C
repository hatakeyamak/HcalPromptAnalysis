// ------------------------------------------------------------------------------------
//  ROOT macro example
//
//  Author : Ken H
//  Written on May 24, 2018
// ------------------------------------------------------------------------------------
//  
// Usage : 
//
//   $ root -b  
//   root> .L ana_PFStudy_TMVA.C+
//   root> ana_PFStudy_TMVA("relval_ttbar_2019.root","pfstudy_histograms_2019.root");
//   or
//   root> ana_PFStudy_TMVA();
//   or
//   from command line:
/*
     root.exe -b -q 'ana_PFStudy_TMVA.C+("relval_ttbar_2019.root","pfstudy_histograms_2019.root")'
 */
//    
// -----------------------------------------------------------------------------------
// 

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
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
#include "TProfile.h"
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

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVARegGui.h"

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
// Book histograms
//
void book1D(TList *v_hist, std::string name, int n, double min, double max);
void book1DProf(TList *v_hist, std::string name, int n, double min, double max, double ymin, double ymax, Option_t *option);
void bookHistograms(TList *v_hist);

//
// Fill histograms
//
void fill1D(TList *v_hist, std::string name, double value);
void fill1DProf(TList *v_hist, std::string name, double value, double valuey);

//
// Main analyzer
//
void PFCheckRun(TString rootfile, TString outfile, int maxevents=-1, int option=2) 
{ 

   cout << "[PF analyzer] Running option " << option << " for " << endl; 

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
   TTreeReaderArray<double> PFParEta = {fReader, "PFParEta"};
   TTreeReaderArray<double> PFParM = {fReader, "PFParM"};
   TTreeReaderArray<double> PFParPhi = {fReader, "PFParPhi"};
   TTreeReaderArray<double> PFParPt = {fReader, "PFParPt"};
   /*
   TTreeReaderArray<float> HBHERecHitEnergy = {fReader, "HBHERecHitEnergy"};
   TTreeReaderArray<float> HBHERecHitEta = {fReader, "HBHERecHitEta"};
   TTreeReaderArray<float> HBHERecHitPhi = {fReader, "HBHERecHitPhi"};
   TTreeReaderArray<float> HBHERecHitTime = {fReader, "HBHERecHitTime"};
   */
   TTreeReaderArray<float> PFParEcalEnergyFrac = {fReader, "PFParEcalEnergyFrac"};
   TTreeReaderArray<float> PFParHOEnergyFrac = {fReader, "PFParHOEnergyFrac"};
   TTreeReaderArray<float> PFParHcalEnergyFrac = {fReader, "PFParHcalEnergyFrac"};
   TTreeReaderArray<float> PFParHcalFrac1 = {fReader, "PFParHcalFrac1"};
   TTreeReaderArray<float> PFParHcalFrac2 = {fReader, "PFParHcalFrac2"};
   TTreeReaderArray<float> PFParHcalFrac3 = {fReader, "PFParHcalFrac3"};
   TTreeReaderArray<float> PFParHcalFrac4 = {fReader, "PFParHcalFrac4"};
   TTreeReaderArray<float> PFParHcalFrac5 = {fReader, "PFParHcalFrac5"};
   TTreeReaderArray<float> PFParHcalFrac6 = {fReader, "PFParHcalFrac6"};
   TTreeReaderArray<float> PFParHcalFrac7 = {fReader, "PFParHcalFrac7"};
   TTreeReaderArray<float> PFParTrackPt = {fReader, "PFParTrackPt"};
   TTreeReaderArray<int> GenParPdgId = {fReader, "GenParPdgId"};
   TTreeReaderArray<int> GenParStatus = {fReader, "GenParStatus"};
   /*
   TTreeReaderArray<int> HBHERecHitAux = {fReader, "HBHERecHitAux"};
   TTreeReaderArray<int> HBHERecHitDepth = {fReader, "HBHERecHitDepth"};
   TTreeReaderArray<int> HBHERecHitFlags = {fReader, "HBHERecHitFlags"};
   TTreeReaderArray<int> HBHERecHitHPDid = {fReader, "HBHERecHitHPDid"};
   TTreeReaderArray<int> HBHERecHitIEta = {fReader, "HBHERecHitIEta"};
   TTreeReaderArray<int> HBHERecHitIPhi = {fReader, "HBHERecHitIPhi"};
   TTreeReaderArray<int> HBHERecHitRBXid = {fReader, "HBHERecHitRBXid"};
   */
   TTreeReaderArray<int> PFParPdgId = {fReader, "PFParPdgId"};
   TTreeReaderArray<int> PFParStatus = {fReader, "PFParStatus"};
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

   //
   // TMVA related
   //
   // https://root.cern.ch/doc/v606/TMVARegression_8C_source.html
   // https://root.cern.ch/doc/v608/TMVARegression_8C.html
   // 
   
   // Create a new root output file
   TString outfileName( "TMVAReg.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory will
   // then run the performance analysis for you.
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );
   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //
   //     (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //     (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";  

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   dataloader->AddVariable( "var1", "Variable 1", "units", 'F' );
   dataloader->AddVariable( "var2", "Variable 2", "units", 'F' );

   // Add the variable carrying the regression target
   dataloader->AddTarget( "fvalue" );

   // Read training and test data (see TMVAClassification for reading ASCII files)
   // load the signal and background event samples from ROOT trees
   TFile *input(0);
   TString fname = "./tmva_reg_example.root";
   if (!gSystem->AccessPathName( fname )) 
     input = TFile::Open( fname ); // check if file in local directory exists
   else 
     input = TFile::Open( "http://root.cern.ch/files/tmva_reg_example.root" ); // if not: download from ROOT server
   
   if (!input) {
     std::cout << "ERROR: could not open data file" << std::endl;
     exit(1);
   }
   std::cout << "--- TMVARegression           : Using input file: " << input->GetName() << std::endl;
   
   // --- Register the regression tree
   
   TTree *regTree = (TTree*)input->Get("TreeR");
   
   // global event weights per tree (see below for setting event-wise weights)
   Double_t regWeight  = 1.0;   
   
   // You can add an arbitrary number of regression trees
   dataloader->AddRegressionTree( regTree, regWeight );
   
   // This would set individual event weights (the variables defined in the 
   // expression need to exist in the original TTree)
   dataloader->SetWeightExpression( "var1", "Regression" );
   
   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycut = ""; // for example: TCut mycut = "abs(var1)<0.5 && abs(var2-0.5)<1";
   
   // tell the factory to use all remaining events in the trees after training for testing:
   dataloader->PrepareTrainingAndTestTree( mycut, 
                                            "nTrain_Regression=1000:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V" );
   // factory->PrepareTrainingAndTestTree( mycut, 
   //                                      "nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V" );
   
   // If no numbers of events are given, half of the events in the tree are used 
   // for training, and the other half for testing:
   //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );  

   // PDE - RS method
   //if (Use["PDERS"])
   //factory->BookMethod( TMVA::Types::kPDERS, "PDERS", 
   //			"!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=40:NEventsMax=60:VarTransform=None" );
   // And the options strings for the MinMax and RMS methods, respectively:
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );   
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );   

   factory->BookMethod( dataloader,  TMVA::Types::kLD, "LD",
			"!H:!V:VarTransform=None" );
   
   // --------------------------------------------------------------------------------------------------
  
   // ---- Now you can tell the factory to train, test, and evaluate the MVAs
  
   // Train MVAs using the set of training events
   factory->TrainAllMethods();
   
   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();
   
   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();    
   
   // --------------------------------------------------------------
  
   // Save the output
   outputFile->Close();
  
   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVARegression is done!" << std::endl;      
  
   delete factory;
   
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
     // Loop over PF candidates
     //--------------------
     for (int ipfcand = 0, npfcand =  PFParPt.GetSize(); ipfcand < npfcand; ++ipfcand) {
       
       std::string strtmp;

       strtmp = "PFTask_hcalFrac1Zero_vs_pt";
       float zero1=0.;
       if (PFParHcalFrac1[ipfcand]==0.) zero1=1.; 
       fill1DProf(v_hist, strtmp, log10(PFParPt[ipfcand]),zero1);

       strtmp = "PFTask_hcalFrac2Zero_vs_pt";
       float zero2=0.;
       if (PFParHcalFrac2[ipfcand]==0.) zero2=1.; 
       fill1DProf(v_hist, strtmp, log10(PFParPt[ipfcand]),zero2);

       strtmp = "PFTask_hcalFrac3Zero_vs_pt";
       float zero3=0.;
       if (PFParHcalFrac3[ipfcand]==0.) zero3=1.; 
       fill1DProf(v_hist, strtmp, log10(PFParPt[ipfcand]),zero3);

       strtmp = "PFTask_hcalFrac4Zero_vs_pt";
       float zero4=0.;
       if (PFParHcalFrac4[ipfcand]==0.) zero4=1.; 
       fill1DProf(v_hist, strtmp, log10(PFParPt[ipfcand]),zero4);

       strtmp = "PFTask_hcalFrac5Zero_vs_pt";
       float zero5=0.;
       if (PFParHcalFrac5[ipfcand]==0.) zero5=1.; 
       fill1DProf(v_hist, strtmp, log10(PFParPt[ipfcand]),zero5);

       strtmp = "PFTask_hcalFrac6Zero_vs_pt";
       float zero6=0.;
       if (PFParHcalFrac6[ipfcand]==0.) zero6=1.;
       fill1DProf(v_hist, strtmp, log10(PFParPt[ipfcand]),zero6);

       strtmp = "PFTask_hcalFrac7Zero_vs_pt";
       float zero7=0.;
       if (PFParHcalFrac7[ipfcand]==0.) zero7=1.; 
       fill1DProf(v_hist, strtmp, log10(PFParPt[ipfcand]),zero7);

       strtmp = "PFTask_hcalFracAllZero_vs_pt";
       float zeroAll=0.;
       if (PFParHcalFrac1[ipfcand]==0.
	   && PFParHcalFrac2[ipfcand]==0. 
	   && PFParHcalFrac3[ipfcand]==0. 
	   && PFParHcalFrac4[ipfcand]==0. 
	   && PFParHcalFrac5[ipfcand]==0. 
	   && PFParHcalFrac6[ipfcand]==0. 
	   && PFParHcalFrac7[ipfcand]==0.) zeroAll=1.; 
       fill1DProf(v_hist, strtmp, log10(PFParPt[ipfcand]), zeroAll);

       //
       // Endcap
       // 
       if ( fabs(PFParEta[ipfcand])<2.9&&fabs(PFParEta[ipfcand])>1.5){

	 //
	 // Charged hadrons
	 //
	 if ( PFParPdgId[ipfcand]==1&&fabs(PFParEta[ipfcand])<2.5 ){ 
	   if ( PFParPt[ipfcand]>5.){
	     /*
	     std::cout << PFParPt[ipfcand] << " " << PFParEta[ipfcand] << " " << PFParPdgId[ipfcand] << std::endl;
	     std::cout << PFParTrackPt[ipfcand] << " " << PFParEcalEnergyFrac[ipfcand] << " "
		       << PFParHcalEnergyFrac[ipfcand] << " " << PFParHOEnergyFrac[ipfcand] << std::endl;
	     std::cout << PFParHcalFrac1[ipfcand] << " "
		       << PFParHcalFrac2[ipfcand] << " " 
		       << PFParHcalFrac3[ipfcand] << " " 
		       << PFParHcalFrac4[ipfcand] << " " 
		       << PFParHcalFrac5[ipfcand] << " " 
		       << PFParHcalFrac6[ipfcand] << " " 
		       << PFParHcalFrac7[ipfcand] << std::endl;
	     std::cout << std::endl;
	     */
	   }

	   //--- 
	   if ( PFParPt[ipfcand]>5. ){
	     strtmp = "PFTask_Profile_ChargedHadron_Endcap_PtAbove5";
	     fill1DProf(v_hist, strtmp, -3, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);
	     fill1DProf(v_hist, strtmp, -2, PFParEcalEnergyFrac[ipfcand]);
	     fill1DProf(v_hist, strtmp, -1, PFParHcalEnergyFrac[ipfcand]);
	     fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	     fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	     fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	     fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	     fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	     fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	     fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);

	     /*
	     std::cout << PFParPt[ipfcand] << " " << PFParEta[ipfcand] << " " << PFParPdgId[ipfcand] << std::endl;	   
	     std::cout << PFParTrackPt[ipfcand] << " " << PFParEcalEnergyFrac[ipfcand] << " "
		       << PFParHcalEnergyFrac[ipfcand] << " " << PFParHOEnergyFrac[ipfcand] << std::endl;
	     std::cout << PFParHcalFrac1[ipfcand] << " "
		       << PFParHcalFrac2[ipfcand] << " " 
		       << PFParHcalFrac3[ipfcand] << " " 
		       << PFParHcalFrac4[ipfcand] << " " 
		       << PFParHcalFrac5[ipfcand] << " " 
		       << PFParHcalFrac6[ipfcand] << " " 
		     << PFParHcalFrac7[ipfcand] << std::endl;
	     std::cout << std::endl;
	     */

	     strtmp = "PFTask_ChargedHadron_TrackPtRatio_PtAbove5";
	     fill1D(v_hist, strtmp, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);

	     if (PFParTrackPt[ipfcand]!=PFParPt[ipfcand]){
	       strtmp = "PFTask_Profile_ChargedHadron_Endcap_PtAbove5_Special";
	       fill1DProf(v_hist, strtmp, -3, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);
	       fill1DProf(v_hist, strtmp, -2, PFParEcalEnergyFrac[ipfcand]);
	       fill1DProf(v_hist, strtmp, -1, PFParHcalEnergyFrac[ipfcand]);
	       fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	       fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	       fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	       fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	       fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	       fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	       fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);
	     }
	     
	   } else if ( PFParPt[ipfcand]>1.){
	     
	     strtmp = "PFTask_Profile_ChargedHadron_Endcap_Pt1To5";
	     fill1DProf(v_hist, strtmp, -3, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);
	     fill1DProf(v_hist, strtmp, -2, PFParEcalEnergyFrac[ipfcand]);
	     fill1DProf(v_hist, strtmp, -1, PFParHcalEnergyFrac[ipfcand]);
	     fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	     fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	     fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	     fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	     fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	     fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	     fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);

	     strtmp = "PFTask_ChargedHadron_TrackPtRatio_Pt1To5";
	     fill1D(v_hist, strtmp, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);
	     
	     if (PFParTrackPt[ipfcand]!=PFParPt[ipfcand]){
	       strtmp = "PFTask_Profile_ChargedHadron_Endcap_Pt1To5_Special";
	       fill1DProf(v_hist, strtmp, -3, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);
	       fill1DProf(v_hist, strtmp, -2, PFParEcalEnergyFrac[ipfcand]);
	       fill1DProf(v_hist, strtmp, -1, PFParHcalEnergyFrac[ipfcand]);
	       fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	       fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	       fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	       fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	       fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	       fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	       fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);
	     }

	   } else {
	     
	     strtmp = "PFTask_Profile_ChargedHadron_Endcap_PtBelow1";
	     fill1DProf(v_hist, strtmp, -3, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);
	     fill1DProf(v_hist, strtmp, -2, PFParEcalEnergyFrac[ipfcand]);
	     fill1DProf(v_hist, strtmp, -1, PFParHcalEnergyFrac[ipfcand]);
	     fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	     fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	     fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	     fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	     fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	     fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	     fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);
	     
	     strtmp = "PFTask_ChargedHadron_TrackPtRatio_PtBelow1";
	     fill1D(v_hist, strtmp, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);

	     if (PFParTrackPt[ipfcand]!=PFParPt[ipfcand]){
	       strtmp = "PFTask_Profile_ChargedHadron_Endcap_PtBelow1_Special";
	       fill1DProf(v_hist, strtmp, -3, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);
	       fill1DProf(v_hist, strtmp, -2, PFParEcalEnergyFrac[ipfcand]);
	       fill1DProf(v_hist, strtmp, -1, PFParHcalEnergyFrac[ipfcand]);
	       fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	       fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	       fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	       fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	       fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	       fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	       fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);
	     }

	   }
	   
	 //
	 // Neutral hadrons
	 //
	 } else if ( PFParPdgId[ipfcand]==5 ){

	   if ( PFParPt[ipfcand]>5.){
	     strtmp = "PFTask_hcalProfile_NeutralHadron_Endcap_PtAbove5";
	     fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	     fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	     fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	     fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	     fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	     fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	     fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);
	     
	     /*
	     std::cout << PFParPt[ipfcand] << " " << PFParEta[ipfcand] << " " << PFParPdgId[ipfcand] << std::endl;	   
	     std::cout << PFParTrackPt[ipfcand] << " " << PFParEcalEnergyFrac[ipfcand] << " "
		       << PFParHcalEnergyFrac[ipfcand] << " " << PFParHOEnergyFrac[ipfcand] << std::endl;
	     std::cout << PFParHcalFrac1[ipfcand] << " "
		       << PFParHcalFrac2[ipfcand] << " " 
		       << PFParHcalFrac3[ipfcand] << " " 
		       << PFParHcalFrac4[ipfcand] << " " 
		       << PFParHcalFrac5[ipfcand] << " " 
		       << PFParHcalFrac6[ipfcand] << " " 
		     << PFParHcalFrac7[ipfcand] << std::endl;
	     std::cout << std::endl;
	     */
	     
	   } else if ( PFParPt[ipfcand]>1.){
	     
	     strtmp = "PFTask_hcalProfile_NeutralHadron_Endcap_Pt1To5";
	     fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	     fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	     fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	     fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	     fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	     fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	     fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);

	   } else {
	     
	     strtmp = "PFTask_hcalProfile_NeutralHadron_Endcap_PtBelow1";
	     fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	     fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	     fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	     fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	     fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	     fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	     fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);
	     
	   }
	   
	 } // Neutral Hadron
       } // Endcap
     } // PF candidiate loop
     
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
//void ana_PFStudy_TMVA(TString rootfile="relval_ttbar_2018_pmx25ns.root",TString outfile="pfstudy_histograms.root",int maxevents=-1)
void ana_PFStudy_TMVA(TString rootfile="relval_ttbar_2018_pmx25ns.root",TString outfile="pfstudy_histograms.root",int maxevents=100)
{
  PFCheckRun(rootfile, outfile, maxevents, 0);
}

//
// --- Aux ---
//

//
// Book 1D histograms
//
void book1D(TList *v_hist, std::string name, int n, double min, double max)
{
  TH1D *htemp = new TH1D(name.c_str(), name.c_str(), n, min, max);
  v_hist->Add(htemp);
}
//
// Book 1D profile histograms
//
void book1DProf(TList *v_hist, std::string name, int n, double min, double max, double ymin, double ymax, Option_t *option="")
{
  TProfile *htemp = new TProfile(name.c_str(), name.c_str(), n, min, max, ymin, ymax, option);
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
  sprintf(histo, "PFTask_hcalFrac1Zero_vs_pt");
  book1DProf(v_hist, histo, 50., -1, 4., -1., 2.);
  sprintf(histo, "PFTask_hcalFrac2Zero_vs_pt");
  book1DProf(v_hist, histo, 50., -1, 4., -1., 2.);
  sprintf(histo, "PFTask_hcalFrac3Zero_vs_pt");
  book1DProf(v_hist, histo, 50., -1, 4., -1., 2.);
  sprintf(histo, "PFTask_hcalFrac4Zero_vs_pt");
  book1DProf(v_hist, histo, 50., -1, 4., -1., 2.);
  sprintf(histo, "PFTask_hcalFrac5Zero_vs_pt");
  book1DProf(v_hist, histo, 50., -1, 4., -1., 2.);
  sprintf(histo, "PFTask_hcalFrac6Zero_vs_pt");
  book1DProf(v_hist, histo, 50., -1, 4., -1., 2.);
  sprintf(histo, "PFTask_hcalFrac7Zero_vs_pt");
  book1DProf(v_hist, histo, 50., -1, 4., -1., 2.);
  sprintf(histo, "PFTask_hcalFracAllZero_vs_pt");
  book1DProf(v_hist, histo, 50., -1, 4., -1., 2.);

  sprintf(histo, "PFTask_ChargedHadron_TrackPtRatio_PtAbove5");
  book1D(v_hist, histo, 70., -0.2, 1.2);
  sprintf(histo, "PFTask_ChargedHadron_TrackPtRatio_Pt1To5");
  book1D(v_hist, histo, 70., -0.2, 1.2);
  sprintf(histo, "PFTask_ChargedHadron_TrackPtRatio_PtBelow1");
  book1D(v_hist, histo, 70., -0.2, 1.2);
  
  sprintf(histo,"PFTask_Profile_ChargedHadron_Endcap_PtAbove5");
  book1DProf(v_hist, histo, 11, -3.5, 7.5, -1., 2., "S");
  TProfile* htemp = (TProfile*) v_hist->FindObject(histo);
  htemp->GetXaxis()->SetBinLabel(1,"Track");
  htemp->GetXaxis()->SetBinLabel(2,"ECAL");
  htemp->GetXaxis()->SetBinLabel(3,"HCAL");
  htemp->GetXaxis()->SetBinLabel(5,"HCAL d1");
  htemp->GetXaxis()->SetBinLabel(6,"HCAL d2");
  htemp->GetXaxis()->SetBinLabel(7,"HCAL d3");
  htemp->GetXaxis()->SetBinLabel(8,"HCAL d4");
  htemp->GetXaxis()->SetBinLabel(9,"HCAL d5");
  htemp->GetXaxis()->SetBinLabel(10,"HCAL d6");
  htemp->GetXaxis()->SetBinLabel(11,"HCAL d7");
  htemp->Print("all");
  sprintf(histo,"PFTask_Profile_ChargedHadron_Endcap_Pt1To5");
  book1DProf(v_hist, histo, 11, -3.5, 7.5, -1., 2., "S");
  htemp = (TProfile*) v_hist->FindObject(histo);
  htemp->GetXaxis()->SetBinLabel(1,"Track");
  htemp->GetXaxis()->SetBinLabel(2,"ECAL");
  htemp->GetXaxis()->SetBinLabel(3,"HCAL");
  htemp->GetXaxis()->SetBinLabel(5,"HCAL d1");
  htemp->GetXaxis()->SetBinLabel(6,"HCAL d2");
  htemp->GetXaxis()->SetBinLabel(7,"HCAL d3");
  htemp->GetXaxis()->SetBinLabel(8,"HCAL d4");
  htemp->GetXaxis()->SetBinLabel(9,"HCAL d5");
  htemp->GetXaxis()->SetBinLabel(10,"HCAL d6");
  htemp->GetXaxis()->SetBinLabel(11,"HCAL d7");
  htemp->Print("all");
  sprintf(histo,"PFTask_Profile_ChargedHadron_Endcap_PtBelow1");
  book1DProf(v_hist, histo, 11, -3.5, 7.5, -1., 2., "S");
  htemp = (TProfile*) v_hist->FindObject(histo);
  htemp->GetXaxis()->SetBinLabel(1,"Track");
  htemp->GetXaxis()->SetBinLabel(2,"ECAL");
  htemp->GetXaxis()->SetBinLabel(3,"HCAL");
  htemp->GetXaxis()->SetBinLabel(5,"HCAL d1");
  htemp->GetXaxis()->SetBinLabel(6,"HCAL d2");
  htemp->GetXaxis()->SetBinLabel(7,"HCAL d3");
  htemp->GetXaxis()->SetBinLabel(8,"HCAL d4");
  htemp->GetXaxis()->SetBinLabel(9,"HCAL d5");
  htemp->GetXaxis()->SetBinLabel(10,"HCAL d6");
  htemp->GetXaxis()->SetBinLabel(11,"HCAL d7");
  htemp->Print("all");

  // PF candidate Pt not 
  sprintf(histo,"PFTask_Profile_ChargedHadron_Endcap_PtAbove5_Special");
  book1DProf(v_hist, histo, 11, -3.5, 7.5, -1., 2., "S");
  htemp = (TProfile*) v_hist->FindObject(histo);
  htemp->GetXaxis()->SetBinLabel(1,"Track");
  htemp->GetXaxis()->SetBinLabel(2,"ECAL");
  htemp->GetXaxis()->SetBinLabel(3,"HCAL");
  htemp->GetXaxis()->SetBinLabel(5,"HCAL d1");
  htemp->GetXaxis()->SetBinLabel(6,"HCAL d2");
  htemp->GetXaxis()->SetBinLabel(7,"HCAL d3");
  htemp->GetXaxis()->SetBinLabel(8,"HCAL d4");
  htemp->GetXaxis()->SetBinLabel(9,"HCAL d5");
  htemp->GetXaxis()->SetBinLabel(10,"HCAL d6");
  htemp->GetXaxis()->SetBinLabel(11,"HCAL d7");
  htemp->Print("all");
  sprintf(histo,"PFTask_Profile_ChargedHadron_Endcap_Pt1To5_Special");
  book1DProf(v_hist, histo, 11, -3.5, 7.5, -1., 2., "S");
  htemp = (TProfile*) v_hist->FindObject(histo);
  htemp->GetXaxis()->SetBinLabel(1,"Track");
  htemp->GetXaxis()->SetBinLabel(2,"ECAL");
  htemp->GetXaxis()->SetBinLabel(3,"HCAL");
  htemp->GetXaxis()->SetBinLabel(5,"HCAL d1");
  htemp->GetXaxis()->SetBinLabel(6,"HCAL d2");
  htemp->GetXaxis()->SetBinLabel(7,"HCAL d3");
  htemp->GetXaxis()->SetBinLabel(8,"HCAL d4");
  htemp->GetXaxis()->SetBinLabel(9,"HCAL d5");
  htemp->GetXaxis()->SetBinLabel(10,"HCAL d6");
  htemp->GetXaxis()->SetBinLabel(11,"HCAL d7");
  htemp->Print("all");
  sprintf(histo,"PFTask_Profile_ChargedHadron_Endcap_PtBelow1_Special");
  book1DProf(v_hist, histo, 11, -3.5, 7.5, -1., 2., "S");
  htemp = (TProfile*) v_hist->FindObject(histo);
  htemp->GetXaxis()->SetBinLabel(1,"Track");
  htemp->GetXaxis()->SetBinLabel(2,"ECAL");
  htemp->GetXaxis()->SetBinLabel(3,"HCAL");
  htemp->GetXaxis()->SetBinLabel(5,"HCAL d1");
  htemp->GetXaxis()->SetBinLabel(6,"HCAL d2");
  htemp->GetXaxis()->SetBinLabel(7,"HCAL d3");
  htemp->GetXaxis()->SetBinLabel(8,"HCAL d4");
  htemp->GetXaxis()->SetBinLabel(9,"HCAL d5");
  htemp->GetXaxis()->SetBinLabel(10,"HCAL d6");
  htemp->GetXaxis()->SetBinLabel(11,"HCAL d7");
  htemp->Print("all");
  
  sprintf(histo,"PFTask_hcalProfile_NeutralHadron_Endcap_PtAbove5");
  book1DProf(v_hist, histo, 7, 0.5,7.5, -1., 2., "S");
  sprintf(histo,"PFTask_hcalProfile_NeutralHadron_Endcap_Pt1To5");
  book1DProf(v_hist, histo, 7, 0.5,7.5, -1., 2., "S");
  sprintf(histo,"PFTask_hcalProfile_NeutralHadron_Endcap_PtBelow1");
  book1DProf(v_hist, histo, 7, 0.5,7.5, -1., 2., "S");
  
}
//
// Fill 1D histograms
//
void fill1D(TList *v_hist, std::string name, double value)
{
  TH1F* htemp = (TH1F*) v_hist->FindObject(name.c_str());
  htemp->Fill(value);
}
//
// Fill 1D Profile histograms
//
void fill1DProf(TList *v_hist, std::string name, double value, double valuey)
{
  TProfile* htemp = (TProfile*) v_hist->FindObject(name.c_str());
  htemp->Fill(value,valuey);
}

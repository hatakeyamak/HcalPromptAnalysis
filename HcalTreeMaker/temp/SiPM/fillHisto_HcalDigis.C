#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TList.h"
#include "TLegendEntry.h"
#include "string.h"
#include <iostream>
#include "TTreeReader.h"
#include "TChain.h"
#include "TProfile.h"
#include <cmath>
#include "TProfile2D.h"
#include "TH2D.h"
#include "TTree.h"

void fillHisto_HcalDigis(const char *infile, const char *outfile, const char *flag = "single"){
   
  TTreeReader *tReader;
 
  if(!strcmp(flag,"single")){

    TFile *f1 = new TFile(infile);
    if(f1->IsZombie()){
      cout << "Root file: " << infile << " not found!" << endl;
      return;
    }
    tReader = new TTreeReader("hcalDigiTree/HcalDigiTree",f1);

  } else if(!strcmp(flag,"list")){

    ifstream in;
    in.open(infile);
    if(!in.is_open()){
      cout << "Cannot open list file: " << infile << endl;
      return;  
    }

    TChain *chain = new TChain("HcalDigiTree");
    
    string line;
    while(in.good()){
      if(!std::getline(in,line)) break; // We read a line from the file
      if(!chain->Add(line.c_str())){
	cout << "Problem loading tree from " << line << endl;
      }else{
	cout << "Adding file: " << line << "..." << endl;
      }
      
    }
    
    in.close();
    
    tReader = new TTreeReader((TTree *)chain);
    cout << "Finished loading files." << endl;
  } else {
    cout << "Unknown option: " << flag << endl;
    return;
  }

  TFile *f2 = new TFile(outfile, "NEW");
  if(f2->IsZombie()){
    cout << "Root file: " << outfile << " cannot be created!" << endl;
    return;
  }
  f2->cd();
      
  //Define the branchs to be read
  TTreeReaderValue<Int_t> eventID(*tReader, "event");
  TTreeReaderValue<Int_t> lumiID(*tReader, "lumi");
  TTreeReaderValue<Int_t> runID(*tReader, "run");  

  TTreeReaderValue<std::vector<float>> DigiHBHE_charge0(*tReader, "DigiHBHE_charge0");
  TTreeReaderValue<std::vector<float>> DigiHBHE_charge1(*tReader, "DigiHBHE_charge1");
  TTreeReaderValue<std::vector<float>> DigiHBHE_charge2(*tReader, "DigiHBHE_charge2");
  TTreeReaderValue<std::vector<float>> DigiHBHE_charge3(*tReader, "DigiHBHE_charge3");
  TTreeReaderValue<std::vector<float>> DigiHBHE_charge4(*tReader, "DigiHBHE_charge4");
  TTreeReaderValue<std::vector<float>> DigiHBHE_charge5(*tReader, "DigiHBHE_charge5");
  TTreeReaderValue<std::vector<float>> DigiHBHE_charge6(*tReader, "DigiHBHE_charge6");
  TTreeReaderValue<std::vector<float>> DigiHBHE_charge7(*tReader, "DigiHBHE_charge7");
  TTreeReaderValue<std::vector<float>> DigiHBHE_charge8(*tReader, "DigiHBHE_charge8");
  TTreeReaderValue<std::vector<float>> DigiHBHE_charge9(*tReader, "DigiHBHE_charge9");

  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_charge0(*tReader, "DigiHBHE_QIE11_charge0");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_charge1(*tReader, "DigiHBHE_QIE11_charge1");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_charge2(*tReader, "DigiHBHE_QIE11_charge2");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_charge3(*tReader, "DigiHBHE_QIE11_charge3");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_charge4(*tReader, "DigiHBHE_QIE11_charge4");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_charge5(*tReader, "DigiHBHE_QIE11_charge5");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_charge6(*tReader, "DigiHBHE_QIE11_charge6");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_charge7(*tReader, "DigiHBHE_QIE11_charge7");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_charge8(*tReader, "DigiHBHE_QIE11_charge8");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_charge9(*tReader, "DigiHBHE_QIE11_charge9");
  
  TTreeReaderValue<std::vector<float>> DigiHO_charge0(*tReader, "DigiHO_charge0");
  TTreeReaderValue<std::vector<float>> DigiHO_charge1(*tReader, "DigiHO_charge1");
  TTreeReaderValue<std::vector<float>> DigiHO_charge2(*tReader, "DigiHO_charge2");
  TTreeReaderValue<std::vector<float>> DigiHO_charge3(*tReader, "DigiHO_charge3");
  TTreeReaderValue<std::vector<float>> DigiHO_charge4(*tReader, "DigiHO_charge4");
  TTreeReaderValue<std::vector<float>> DigiHO_charge5(*tReader, "DigiHO_charge5");
  TTreeReaderValue<std::vector<float>> DigiHO_charge6(*tReader, "DigiHO_charge6");
  TTreeReaderValue<std::vector<float>> DigiHO_charge7(*tReader, "DigiHO_charge7");
  TTreeReaderValue<std::vector<float>> DigiHO_charge8(*tReader, "DigiHO_charge8");
  TTreeReaderValue<std::vector<float>> DigiHO_charge9(*tReader, "DigiHO_charge9");

  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_charge0(*tReader, "DigiHF_QIE10_charge0");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_charge1(*tReader, "DigiHF_QIE10_charge1");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_charge2(*tReader, "DigiHF_QIE10_charge2");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_charge3(*tReader, "DigiHF_QIE10_charge3");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_charge4(*tReader, "DigiHF_QIE10_charge4");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_charge5(*tReader, "DigiHF_QIE10_charge5");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_charge6(*tReader, "DigiHF_QIE10_charge6");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_charge7(*tReader, "DigiHF_QIE10_charge7");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_charge8(*tReader, "DigiHF_QIE10_charge8");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_charge9(*tReader, "DigiHF_QIE10_charge9");

  //Define the Histograms to be filled
  TProfile *pulseShape_DigiHBHE = new TProfile("pulseShape_DigiHBHE","pulseShape_DigiHBHE",10,-0.5,9.5,-1,1);
  TProfile *pulseShape_DigiHBHE_QIE11 = new TProfile("pulseShape_DigiHBHE_QIE11","pulseShape_DigiHBHE_QIE11",10,-0.5,9.5,-1,1);
  TProfile *pulseShape_DigiHO = new TProfile("pulseShape_DigiHO","pulseShape_DigiHO",10,-0.5,9.5,-1,1);
  TProfile *pulseShape_DigiHF_QIE10 = new TProfile("pulseShape_DigiHF_QIE10","pulseShape_DigiHF_QIE10",10,-0.5,9.5,-1,1);

  TH1D *DigiHBHE_charge0_hist = new TH1D("DigiHBHE_charge0_hist","DigiHBHE_charge0_hist",60,-20,20);

  cout << "Filling histograms..." << endl;
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                           While loop over all of the events and fills the histograms                                         //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  while(tReader->Next()){

    double n0 = 0;
    double n1 = 1;
    double n2 = 2;
    double n3 = 3;
    double n4 = 4;
    double n5 = 5;
    double n6 = 6;
    double n7 = 7;
    double n8 = 8;
    double n9 = 9;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                         HBHE_PulseShape                                                                    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double sumHBHE0 = 0;
    double sumHBHE1 = 0;
    double sumHBHE2 = 0;
    double sumHBHE3 = 0;
    double sumHBHE4 = 0;
    double sumHBHE5 = 0;
    double sumHBHE6 = 0;
    double sumHBHE7 = 0;
    double sumHBHE8 = 0;
    double sumHBHE9 = 0;

    for(int digi = 0; digi < (*DigiHBHE_charge0).size(); digi++){
      sumHBHE0 = sumHBHE0 + (*DigiHBHE_charge0).at(digi);
      DigiHBHE_charge0_hist->Fill((*DigiHBHE_charge0).at(digi));
    }//HBHE T0
    for(int digi = 0; digi < (*DigiHBHE_charge1).size(); digi++){
      sumHBHE1 = sumHBHE1 + (*DigiHBHE_charge1).at(digi);
    }//HBHE T1
    for(int digi = 0; digi < (*DigiHBHE_charge2).size(); digi++){
      sumHBHE2 = sumHBHE2 + (*DigiHBHE_charge2).at(digi);
    }//HBHE T2
    for(int digi = 0; digi < (*DigiHBHE_charge3).size(); digi++){
      sumHBHE3 = sumHBHE3 + (*DigiHBHE_charge3).at(digi);
    }//HBHE T3
    for(int digi = 0; digi < (*DigiHBHE_charge4).size(); digi++){
      sumHBHE4 = sumHBHE4 + (*DigiHBHE_charge4).at(digi);
    }//HBHE T4
    for(int digi = 0; digi < (*DigiHBHE_charge5).size(); digi++){
      sumHBHE5 = sumHBHE5 + (*DigiHBHE_charge5).at(digi);
    }//HBHE T5
    for(int digi = 0; digi < (*DigiHBHE_charge6).size(); digi++){
      sumHBHE6 = sumHBHE6 + (*DigiHBHE_charge6).at(digi);
    }//HBHE T6
    for(int digi = 0; digi < (*DigiHBHE_charge7).size(); digi++){
      sumHBHE7 = sumHBHE7 + (*DigiHBHE_charge7).at(digi);
    }//HBHE T7
    for(int digi = 0; digi < (*DigiHBHE_charge8).size(); digi++){
      sumHBHE8 = sumHBHE8 + (*DigiHBHE_charge8).at(digi);
    }//HBHE T8
    for(int digi = 0; digi < (*DigiHBHE_charge9).size(); digi++){
      sumHBHE9 = sumHBHE9 + (*DigiHBHE_charge9).at(digi);
    }//HBHE T9
    
    double totsumHBHE = sumHBHE0 + sumHBHE1 + sumHBHE2 + sumHBHE3 + sumHBHE4 + sumHBHE5 + sumHBHE6 + sumHBHE7 + sumHBHE8 + sumHBHE9;
    
    ///Total charge must be greater than 5000 fC
    //if(totsumHBHE > 5000){
      pulseShape_DigiHBHE->Fill(n0,sumHBHE0/totsumHBHE);
      pulseShape_DigiHBHE->Fill(n1,sumHBHE1/totsumHBHE);
      pulseShape_DigiHBHE->Fill(n2,sumHBHE2/totsumHBHE);
      pulseShape_DigiHBHE->Fill(n3,sumHBHE3/totsumHBHE);
      pulseShape_DigiHBHE->Fill(n4,sumHBHE4/totsumHBHE);
      pulseShape_DigiHBHE->Fill(n5,sumHBHE5/totsumHBHE);
      pulseShape_DigiHBHE->Fill(n6,sumHBHE6/totsumHBHE);
      pulseShape_DigiHBHE->Fill(n7,sumHBHE7/totsumHBHE);
      pulseShape_DigiHBHE->Fill(n8,sumHBHE8/totsumHBHE);
      pulseShape_DigiHBHE->Fill(n9,sumHBHE9/totsumHBHE);
    //}

    /*
    cout <<  "x=" << 0 << "y=" << sumHBHE0/totsumHBHE << endl;
    cout <<  "x=" << 1 << "y=" << sumHBHE1/totsumHBHE << endl;
    cout <<  "x=" << 2 << "y=" << sumHBHE2/totsumHBHE << endl;
    cout <<  "x=" << 3 << "y=" << sumHBHE3/totsumHBHE << endl;
    cout <<  "x=" << 4 << "y=" << sumHBHE4/totsumHBHE << endl;
    cout <<  "x=" << 5 << "y=" << sumHBHE5/totsumHBHE << endl;
    cout <<  "x=" << 6 << "y=" << sumHBHE6/totsumHBHE << endl;
    cout <<  "x=" << 7 << "y=" << sumHBHE7/totsumHBHE << endl;
    cout <<  "x=" << 8 << "y=" << sumHBHE8/totsumHBHE << endl;
    cout <<  "x=" << 9 << "y=" << sumHBHE9/totsumHBHE << endl;
    */
    
    //int size = (*DigiHBHE_charge0).size();
    //cout << size  << endl; 
     
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                         HBHE_QIE11_PulseShape                                                              //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double sumHBHE_QIE110 = 0;
    double sumHBHE_QIE111 = 0;
    double sumHBHE_QIE112 = 0;
    double sumHBHE_QIE113 = 0;
    double sumHBHE_QIE114 = 0;
    double sumHBHE_QIE115 = 0;
    double sumHBHE_QIE116 = 0;
    double sumHBHE_QIE117 = 0;
    double sumHBHE_QIE118 = 0;
    double sumHBHE_QIE119 = 0;

    for(int digi = 0; digi < (*DigiHBHE_QIE11_charge0).size(); digi++){
      sumHBHE_QIE110 = sumHBHE_QIE110 + (*DigiHBHE_QIE11_charge0).at(digi);
    }//HBHE_QIE11 T0
    for(int digi = 0; digi < (*DigiHBHE_QIE11_charge1).size(); digi++){
      sumHBHE_QIE111 = sumHBHE_QIE111 + (*DigiHBHE_QIE11_charge1).at(digi);
    }//HBHE_QIE11 T1
    for(int digi = 0; digi < (*DigiHBHE_QIE11_charge2).size(); digi++){
      sumHBHE_QIE112 = sumHBHE_QIE112 + (*DigiHBHE_QIE11_charge2).at(digi);
    }//HBHE_QIE11 T2
    for(int digi = 0; digi < (*DigiHBHE_QIE11_charge3).size(); digi++){
      sumHBHE_QIE113 = sumHBHE_QIE113 + (*DigiHBHE_QIE11_charge3).at(digi);
    }//HBHE_QIE11 T3
    for(int digi = 0; digi < (*DigiHBHE_QIE11_charge4).size(); digi++){
      sumHBHE_QIE114 = sumHBHE_QIE114 + (*DigiHBHE_QIE11_charge4).at(digi);
    }//HBHE_QIE11 T4
    for(int digi = 0; digi < (*DigiHBHE_QIE11_charge5).size(); digi++){
      sumHBHE_QIE115 = sumHBHE_QIE115 + (*DigiHBHE_QIE11_charge5).at(digi);
    }//HBHE_QIE11 T5
    for(int digi = 0; digi < (*DigiHBHE_QIE11_charge6).size(); digi++){
      sumHBHE_QIE116 = sumHBHE_QIE116 + (*DigiHBHE_QIE11_charge6).at(digi);
    }//HBHE_QIE11 T6
    for(int digi = 0; digi < (*DigiHBHE_QIE11_charge7).size(); digi++){
      sumHBHE_QIE117 = sumHBHE_QIE117 + (*DigiHBHE_QIE11_charge7).at(digi);
    }//HBHE_QIE11 T7
    for(int digi = 0; digi < (*DigiHBHE_QIE11_charge8).size(); digi++){
      sumHBHE_QIE118 = sumHBHE_QIE118 + (*DigiHBHE_QIE11_charge8).at(digi);
    }//HBHE_QIE11 T8
    for(int digi = 0; digi < (*DigiHBHE_QIE11_charge9).size(); digi++){
      sumHBHE_QIE119 = sumHBHE_QIE119 + (*DigiHBHE_QIE11_charge9).at(digi);
    }//HBHE_QIE11 T9
    
    double totsumHBHE_QIE11 = sumHBHE_QIE110 + sumHBHE_QIE111 + sumHBHE_QIE112 + sumHBHE_QIE113 + sumHBHE_QIE114 + sumHBHE_QIE115 + sumHBHE_QIE116 + sumHBHE_QIE117 + sumHBHE_QIE118 + sumHBHE_QIE119;
    
    ///Total charge must be greater than 5000 fC
    //if(totsumHBHE_QIE11 > 5000){
      pulseShape_DigiHBHE_QIE11->Fill(n0,sumHBHE_QIE110/totsumHBHE_QIE11);
      pulseShape_DigiHBHE_QIE11->Fill(n1,sumHBHE_QIE111/totsumHBHE_QIE11);
      pulseShape_DigiHBHE_QIE11->Fill(n2,sumHBHE_QIE112/totsumHBHE_QIE11);
      pulseShape_DigiHBHE_QIE11->Fill(n3,sumHBHE_QIE113/totsumHBHE_QIE11);
      pulseShape_DigiHBHE_QIE11->Fill(n4,sumHBHE_QIE114/totsumHBHE_QIE11);
      pulseShape_DigiHBHE_QIE11->Fill(n5,sumHBHE_QIE115/totsumHBHE_QIE11);
      pulseShape_DigiHBHE_QIE11->Fill(n6,sumHBHE_QIE116/totsumHBHE_QIE11);
      pulseShape_DigiHBHE_QIE11->Fill(n7,sumHBHE_QIE117/totsumHBHE_QIE11);
      pulseShape_DigiHBHE_QIE11->Fill(n8,sumHBHE_QIE118/totsumHBHE_QIE11);
      pulseShape_DigiHBHE_QIE11->Fill(n9,sumHBHE_QIE119/totsumHBHE_QIE11);
    //}

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                         HO_PulseShape                                                                      //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double sumHO0 = 0;
    double sumHO1 = 0;
    double sumHO2 = 0;
    double sumHO3 = 0;
    double sumHO4 = 0;
    double sumHO5 = 0;
    double sumHO6 = 0;
    double sumHO7 = 0;
    double sumHO8 = 0;
    double sumHO9 = 0;

    for(int digi = 0; digi < (*DigiHO_charge0).size(); digi++){
      sumHO0 = sumHO0 + (*DigiHO_charge0).at(digi);
    }//HO T0
    for(int digi = 0; digi < (*DigiHO_charge1).size(); digi++){
      sumHO1 = sumHO1 + (*DigiHO_charge1).at(digi);
    }//HO T1
    for(int digi = 0; digi < (*DigiHO_charge2).size(); digi++){
      sumHO2 = sumHO2 + (*DigiHO_charge2).at(digi);
    }//HO T2
    for(int digi = 0; digi < (*DigiHO_charge3).size(); digi++){
      sumHO3 = sumHO3 + (*DigiHO_charge3).at(digi);
    }//HO T3
    for(int digi = 0; digi < (*DigiHO_charge4).size(); digi++){
      sumHO4 = sumHO4 + (*DigiHO_charge4).at(digi);
    }//HO T4
    for(int digi = 0; digi < (*DigiHO_charge5).size(); digi++){
      sumHO5 = sumHO5 + (*DigiHO_charge5).at(digi);
    }//HO T5
    for(int digi = 0; digi < (*DigiHO_charge6).size(); digi++){
      sumHO6 = sumHO6 + (*DigiHO_charge6).at(digi);
    }//HO T6
    for(int digi = 0; digi < (*DigiHO_charge7).size(); digi++){
      sumHO7 = sumHO7 + (*DigiHO_charge7).at(digi);
    }//HO T7
    for(int digi = 0; digi < (*DigiHO_charge8).size(); digi++){
      sumHO8 = sumHO8 + (*DigiHO_charge8).at(digi);
    }//HO T8
    for(int digi = 0; digi < (*DigiHO_charge9).size(); digi++){
      sumHO9 = sumHO9 + (*DigiHO_charge9).at(digi);
    }//HO T9
    
    double totsumHO = sumHO0 + sumHO1 + sumHO2 + sumHO3 + sumHO4 + sumHO5 + sumHO6 + sumHO7 + sumHO8 + sumHO9;
    
    ///Total charge must be greater than 5000 fC
    //if(totsumHO > 5000){
      pulseShape_DigiHO->Fill(n0,sumHO0/totsumHO);
      pulseShape_DigiHO->Fill(n1,sumHO1/totsumHO);
      pulseShape_DigiHO->Fill(n2,sumHO2/totsumHO);
      pulseShape_DigiHO->Fill(n3,sumHO3/totsumHO);
      pulseShape_DigiHO->Fill(n4,sumHO4/totsumHO);
      pulseShape_DigiHO->Fill(n5,sumHO5/totsumHO);
      pulseShape_DigiHO->Fill(n6,sumHO6/totsumHO);
      pulseShape_DigiHO->Fill(n7,sumHO7/totsumHO);
      pulseShape_DigiHO->Fill(n8,sumHO8/totsumHO);
      pulseShape_DigiHO->Fill(n9,sumHO9/totsumHO);
    //}

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                         HF_QIE10_PulseShape                                                                //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double sumHF_QIE100 = 0;
    double sumHF_QIE101 = 0;
    double sumHF_QIE102 = 0;
    double sumHF_QIE103 = 0;
    double sumHF_QIE104 = 0;
    double sumHF_QIE105 = 0;
    double sumHF_QIE106 = 0;
    double sumHF_QIE107 = 0;
    double sumHF_QIE108 = 0;
    double sumHF_QIE109 = 0;

    for(int digi = 0; digi < (*DigiHF_QIE10_charge0).size(); digi++){
      sumHF_QIE100 = sumHF_QIE100 + (*DigiHF_QIE10_charge0).at(digi);
    }//HF_QIE10 T0
    for(int digi = 0; digi < (*DigiHF_QIE10_charge1).size(); digi++){
      sumHF_QIE101 = sumHF_QIE101 + (*DigiHF_QIE10_charge1).at(digi);
    }//HF_QIE10 T1
    for(int digi = 0; digi < (*DigiHF_QIE10_charge2).size(); digi++){
      sumHF_QIE102 = sumHF_QIE102 + (*DigiHF_QIE10_charge2).at(digi);
    }//HF_QIE10 T2
    for(int digi = 0; digi < (*DigiHF_QIE10_charge3).size(); digi++){
      sumHF_QIE103 = sumHF_QIE103 + (*DigiHF_QIE10_charge3).at(digi);
    }//HF_QIE10 T3
    for(int digi = 0; digi < (*DigiHF_QIE10_charge4).size(); digi++){
      sumHF_QIE104 = sumHF_QIE104 + (*DigiHF_QIE10_charge4).at(digi);
    }//HF_QIE10 T4
    for(int digi = 0; digi < (*DigiHF_QIE10_charge5).size(); digi++){
      sumHF_QIE105 = sumHF_QIE105 + (*DigiHF_QIE10_charge5).at(digi);
    }//HF_QIE10 T5
    for(int digi = 0; digi < (*DigiHF_QIE10_charge6).size(); digi++){
      sumHF_QIE106 = sumHF_QIE106 + (*DigiHF_QIE10_charge6).at(digi);
    }//HF_QIE10 T6
    for(int digi = 0; digi < (*DigiHF_QIE10_charge7).size(); digi++){
      sumHF_QIE107 = sumHF_QIE107 + (*DigiHF_QIE10_charge7).at(digi);
    }//HF_QIE10 T7
    for(int digi = 0; digi < (*DigiHF_QIE10_charge8).size(); digi++){
      sumHF_QIE108 = sumHF_QIE108 + (*DigiHF_QIE10_charge8).at(digi);
    }//HF_QIE10 T8
    for(int digi = 0; digi < (*DigiHF_QIE10_charge9).size(); digi++){
      sumHF_QIE109 = sumHF_QIE109 + (*DigiHF_QIE10_charge9).at(digi);
    }//HF_QIE10 T9
    
    double totsumHF_QIE10 = sumHF_QIE100 + sumHF_QIE101 + sumHF_QIE102 + sumHF_QIE103 + sumHF_QIE104 + sumHF_QIE105 + sumHF_QIE106 + sumHF_QIE107 + sumHF_QIE108 + sumHF_QIE109;
    
    ///Total charge must be greater than 5000 fC
    //if(totsumHF_QIE10 > 5000){
      pulseShape_DigiHF_QIE10->Fill(n0,sumHF_QIE100/totsumHF_QIE10);
      pulseShape_DigiHF_QIE10->Fill(n1,sumHF_QIE101/totsumHF_QIE10);
      pulseShape_DigiHF_QIE10->Fill(n2,sumHF_QIE102/totsumHF_QIE10);
      pulseShape_DigiHF_QIE10->Fill(n3,sumHF_QIE103/totsumHF_QIE10);
      pulseShape_DigiHF_QIE10->Fill(n4,sumHF_QIE104/totsumHF_QIE10);
      pulseShape_DigiHF_QIE10->Fill(n5,sumHF_QIE105/totsumHF_QIE10);
      pulseShape_DigiHF_QIE10->Fill(n6,sumHF_QIE106/totsumHF_QIE10);
      pulseShape_DigiHF_QIE10->Fill(n7,sumHF_QIE107/totsumHF_QIE10);
      pulseShape_DigiHF_QIE10->Fill(n8,sumHF_QIE108/totsumHF_QIE10);
      pulseShape_DigiHF_QIE10->Fill(n9,sumHF_QIE109/totsumHF_QIE10);
    //}
 
  }//End of While Loop

  //Fill Histograms
  f2->Write();

  cout << "End of Macro" << endl;
}

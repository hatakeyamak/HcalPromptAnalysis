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

  TTreeReaderValue<std::vector<Int_t>> DigiHBHE_sub(*tReader, "DigiHBHE_sub");
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
  TTreeReaderValue<std::vector<float>> DigiHBHE_adc0(*tReader, "DigiHBHE_adc0");
  TTreeReaderValue<std::vector<float>> DigiHBHE_adc1(*tReader, "DigiHBHE_adc1");
  TTreeReaderValue<std::vector<float>> DigiHBHE_adc2(*tReader, "DigiHBHE_adc2");
  TTreeReaderValue<std::vector<float>> DigiHBHE_adc3(*tReader, "DigiHBHE_adc3");
  TTreeReaderValue<std::vector<float>> DigiHBHE_adc4(*tReader, "DigiHBHE_adc4");
  TTreeReaderValue<std::vector<float>> DigiHBHE_adc5(*tReader, "DigiHBHE_adc5");
  TTreeReaderValue<std::vector<float>> DigiHBHE_adc6(*tReader, "DigiHBHE_adc6");
  TTreeReaderValue<std::vector<float>> DigiHBHE_adc7(*tReader, "DigiHBHE_adc7");
  TTreeReaderValue<std::vector<float>> DigiHBHE_adc8(*tReader, "DigiHBHE_adc8");
  TTreeReaderValue<std::vector<float>> DigiHBHE_adc9(*tReader, "DigiHBHE_adc9");

  TTreeReaderValue<std::vector<Int_t>> DigiHBHE_QIE11_sub(*tReader, "DigiHBHE_QIE11_sub");
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
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_adc0(*tReader, "DigiHBHE_QIE11_adc0");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_adc1(*tReader, "DigiHBHE_QIE11_adc1");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_adc2(*tReader, "DigiHBHE_QIE11_adc2");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_adc3(*tReader, "DigiHBHE_QIE11_adc3");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_adc4(*tReader, "DigiHBHE_QIE11_adc4");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_adc5(*tReader, "DigiHBHE_QIE11_adc5");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_adc6(*tReader, "DigiHBHE_QIE11_adc6");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_adc7(*tReader, "DigiHBHE_QIE11_adc7");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_adc8(*tReader, "DigiHBHE_QIE11_adc8");
  TTreeReaderValue<std::vector<float>> DigiHBHE_QIE11_adc9(*tReader, "DigiHBHE_QIE11_adc9");

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
  TTreeReaderValue<std::vector<float>> DigiHO_adc0(*tReader, "DigiHO_adc0");
  TTreeReaderValue<std::vector<float>> DigiHO_adc1(*tReader, "DigiHO_adc1");
  TTreeReaderValue<std::vector<float>> DigiHO_adc2(*tReader, "DigiHO_adc2");
  TTreeReaderValue<std::vector<float>> DigiHO_adc3(*tReader, "DigiHO_adc3");
  TTreeReaderValue<std::vector<float>> DigiHO_adc4(*tReader, "DigiHO_adc4");
  TTreeReaderValue<std::vector<float>> DigiHO_adc5(*tReader, "DigiHO_adc5");
  TTreeReaderValue<std::vector<float>> DigiHO_adc6(*tReader, "DigiHO_adc6");
  TTreeReaderValue<std::vector<float>> DigiHO_adc7(*tReader, "DigiHO_adc7");
  TTreeReaderValue<std::vector<float>> DigiHO_adc8(*tReader, "DigiHO_adc8");
  TTreeReaderValue<std::vector<float>> DigiHO_adc9(*tReader, "DigiHO_adc9");

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
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_adc0(*tReader, "DigiHF_QIE10_adc0");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_adc1(*tReader, "DigiHF_QIE10_adc1");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_adc2(*tReader, "DigiHF_QIE10_adc2");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_adc3(*tReader, "DigiHF_QIE10_adc3");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_adc4(*tReader, "DigiHF_QIE10_adc4");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_adc5(*tReader, "DigiHF_QIE10_adc5");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_adc6(*tReader, "DigiHF_QIE10_adc6");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_adc7(*tReader, "DigiHF_QIE10_adc7");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_adc8(*tReader, "DigiHF_QIE10_adc8");
  TTreeReaderValue<std::vector<float>> DigiHF_QIE10_adc9(*tReader, "DigiHF_QIE10_adc9");

  //Define the Histograms to be filled
  TProfile *pulseShape_DigiHB_HPD = new TProfile("pulseShape_DigiHB_HPD","pulseShape_DigiHB_HPD",10,-0.5,9.5,-1,1);
  TProfile *pulseShape_DigiHE_HPD = new TProfile("pulseShape_DigiHE_HPD","pulseShape_DigiHE_HPD",10,-0.5,9.5,-1,1);
  TProfile *pulseShape_adcHB_HPD = new TProfile("pulseShape_adcHB_HPD","pulseShape_adcHB_HPD",10,-0.5,9.5,-1,1);
  TProfile *pulseShape_adcHE_HPD = new TProfile("pulseShape_adcHE_HPD","pulseShape_adcHE_HPD",10,-0.5,9.5,-1,1);

  TProfile *pulseShape_DigiHB_SiPM = new TProfile("pulseShape_DigiHB_SiPM","pulseShape_DigiHB_SiPM",10,-0.5,9.5,-1,1);
  TProfile *pulseShape_DigiHE_SiPM = new TProfile("pulseShape_DigiHE_SiPM","pulseShape_DigiHE_SiPM",10,-0.5,9.5,-1,1);
  TProfile *pulseShape_adcHB_SiPM = new TProfile("pulseShape_adcHB_SiPM","pulseShape_adcHB_SiPM",10,-0.5,9.5,-1,1);
  TProfile *pulseShape_adcHE_SiPM = new TProfile("pulseShape_adcHE_SiPM","pulseShape_adcHE_SiPM",10,-0.5,9.5,-1,1);

  TProfile *pulseShape_DigiHO = new TProfile("pulseShape_DigiHO","pulseShape_DigiHO",10,-0.5,9.5,-1,1);
  TProfile *pulseShape_adcHO = new TProfile("pulseShape_adcHO","pulseShape_adcHO",10,-0.5,9.5,-1,1);

  TProfile *pulseShape_DigiHF_QIE10 = new TProfile("pulseShape_DigiHF_QIE10","pulseShape_DigiHF_QIE10",10,-0.5,9.5,-1,1);
  TProfile *pulseShape_adcHF_QIE10 = new TProfile("pulseShape_adcHF_QIE10","pulseShape_adcHF_QIE10",10,-0.5,9.5,-1,1);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                           While loop over all of the events and fills the histograms                                         //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  while(tReader->Next()){

    double min_totcharge_SiPM = 5000; //32000; //5000fC (~2GeV)
    double min_totcharge_HPD  = 50; //600;  //Not Sure where to put this 
    double min_totcharge_HF   = 50; //3200;
    
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
    double totsumHBHE = 0;

    double sumHBHE_adc0 = 0;
    double sumHBHE_adc1 = 0;
    double sumHBHE_adc2 = 0;
    double sumHBHE_adc3 = 0;
    double sumHBHE_adc4 = 0;
    double sumHBHE_adc5 = 0;
    double sumHBHE_adc6 = 0;
    double sumHBHE_adc7 = 0;
    double sumHBHE_adc8 = 0;
    double sumHBHE_adc9 = 0;
    double totsumHBHE_adc = 0;

    //cout  << (*DigiHBHE_charge0).size() << endl;

    for(int digi = 0; digi < (*DigiHBHE_charge0).size(); digi++){
      if((*DigiHBHE_charge0).size() > 0) sumHBHE0 = (*DigiHBHE_charge0).at(digi);
      if((*DigiHBHE_charge1).size() > 0) sumHBHE1 = (*DigiHBHE_charge1).at(digi);
      if((*DigiHBHE_charge2).size() > 0) sumHBHE2 = (*DigiHBHE_charge2).at(digi);
      if((*DigiHBHE_charge3).size() > 0) sumHBHE3 = (*DigiHBHE_charge3).at(digi);
      if((*DigiHBHE_charge4).size() > 0) sumHBHE4 = (*DigiHBHE_charge4).at(digi);
      if((*DigiHBHE_charge5).size() > 0) sumHBHE5 = (*DigiHBHE_charge5).at(digi);
      if((*DigiHBHE_charge6).size() > 0) sumHBHE6 = (*DigiHBHE_charge6).at(digi);
      if((*DigiHBHE_charge7).size() > 0) sumHBHE7 = (*DigiHBHE_charge7).at(digi);
      if((*DigiHBHE_charge8).size() > 0) sumHBHE8 = (*DigiHBHE_charge8).at(digi);
      if((*DigiHBHE_charge9).size() > 0) sumHBHE9 = (*DigiHBHE_charge9).at(digi);
      totsumHBHE = sumHBHE4 + sumHBHE5;// + sumHBHE0 + sumHBHE1 + sumHBHE2 + sumHBHE3 + sumHBHE6 + sumHBHE7 + sumHBHE8 + sumHBHE9;

      if((*DigiHBHE_adc0).size() > 0) sumHBHE_adc0 = (*DigiHBHE_adc0).at(digi);
      if((*DigiHBHE_adc1).size() > 0) sumHBHE_adc1 = (*DigiHBHE_adc1).at(digi);
      if((*DigiHBHE_adc2).size() > 0) sumHBHE_adc2 = (*DigiHBHE_adc2).at(digi);
      if((*DigiHBHE_adc3).size() > 0) sumHBHE_adc3 = (*DigiHBHE_adc3).at(digi);
      if((*DigiHBHE_adc4).size() > 0) sumHBHE_adc4 = (*DigiHBHE_adc4).at(digi);
      if((*DigiHBHE_adc5).size() > 0) sumHBHE_adc5 = (*DigiHBHE_adc5).at(digi);
      if((*DigiHBHE_adc6).size() > 0) sumHBHE_adc6 = (*DigiHBHE_adc6).at(digi);
      if((*DigiHBHE_adc7).size() > 0) sumHBHE_adc7 = (*DigiHBHE_adc7).at(digi);
      if((*DigiHBHE_adc8).size() > 0) sumHBHE_adc8 = (*DigiHBHE_adc8).at(digi);
      if((*DigiHBHE_adc9).size() > 0) sumHBHE_adc9 = (*DigiHBHE_adc9).at(digi);
      totsumHBHE_adc = sumHBHE_adc4 + sumHBHE_adc5;// + sumHBHE_adc0 + sumHBHE_adc1 + sumHBHE_adc2 + sumHBHE_adc3 + sumHBHE_adc6 + sumHBHE_adc7 + sumHBHE_adc8 + sumHBHE_adc9;

      ///Total charge must be greater than 5000 fC
      if(totsumHBHE > min_totcharge_HPD){
	
	if((*DigiHBHE_sub).at(digi) == 1){
	  pulseShape_DigiHB_HPD->Fill(0.,sumHBHE0/totsumHBHE);
	  pulseShape_DigiHB_HPD->Fill(1.,sumHBHE1/totsumHBHE);
	  pulseShape_DigiHB_HPD->Fill(2.,sumHBHE2/totsumHBHE);
	  pulseShape_DigiHB_HPD->Fill(3.,sumHBHE3/totsumHBHE);
	  pulseShape_DigiHB_HPD->Fill(4.,sumHBHE4/totsumHBHE);
	  pulseShape_DigiHB_HPD->Fill(5.,sumHBHE5/totsumHBHE);
	  pulseShape_DigiHB_HPD->Fill(6.,sumHBHE6/totsumHBHE);
	  pulseShape_DigiHB_HPD->Fill(7.,sumHBHE7/totsumHBHE);
	  pulseShape_DigiHB_HPD->Fill(8.,sumHBHE8/totsumHBHE);
	  pulseShape_DigiHB_HPD->Fill(9.,sumHBHE9/totsumHBHE);

	  pulseShape_adcHB_HPD->Fill(0.,sumHBHE_adc0/totsumHBHE_adc);
	  pulseShape_adcHB_HPD->Fill(1.,sumHBHE_adc1/totsumHBHE_adc);
	  pulseShape_adcHB_HPD->Fill(2.,sumHBHE_adc2/totsumHBHE_adc);
	  pulseShape_adcHB_HPD->Fill(3.,sumHBHE_adc3/totsumHBHE_adc);
	  pulseShape_adcHB_HPD->Fill(4.,sumHBHE_adc4/totsumHBHE_adc);
	  pulseShape_adcHB_HPD->Fill(5.,sumHBHE_adc5/totsumHBHE_adc);
	  pulseShape_adcHB_HPD->Fill(6.,sumHBHE_adc6/totsumHBHE_adc);
	  pulseShape_adcHB_HPD->Fill(7.,sumHBHE_adc7/totsumHBHE_adc);
	  pulseShape_adcHB_HPD->Fill(8.,sumHBHE_adc8/totsumHBHE_adc);
	  pulseShape_adcHB_HPD->Fill(9.,sumHBHE_adc9/totsumHBHE_adc);
	}//HB HPD

	if((*DigiHBHE_sub).at(digi) == 2){
	  pulseShape_DigiHE_HPD->Fill(0.,sumHBHE0/totsumHBHE);
	  pulseShape_DigiHE_HPD->Fill(1.,sumHBHE1/totsumHBHE);
	  pulseShape_DigiHE_HPD->Fill(2.,sumHBHE2/totsumHBHE);
	  pulseShape_DigiHE_HPD->Fill(3.,sumHBHE3/totsumHBHE);
	  pulseShape_DigiHE_HPD->Fill(4.,sumHBHE4/totsumHBHE);
	  pulseShape_DigiHE_HPD->Fill(5.,sumHBHE5/totsumHBHE);
	  pulseShape_DigiHE_HPD->Fill(6.,sumHBHE6/totsumHBHE);
	  pulseShape_DigiHE_HPD->Fill(7.,sumHBHE7/totsumHBHE);
	  pulseShape_DigiHE_HPD->Fill(8.,sumHBHE8/totsumHBHE);
	  pulseShape_DigiHE_HPD->Fill(9.,sumHBHE9/totsumHBHE);

	  pulseShape_adcHE_HPD->Fill(0.,sumHBHE_adc0/totsumHBHE_adc);
	  pulseShape_adcHE_HPD->Fill(1.,sumHBHE_adc1/totsumHBHE_adc);
	  pulseShape_adcHE_HPD->Fill(2.,sumHBHE_adc2/totsumHBHE_adc);
	  pulseShape_adcHE_HPD->Fill(3.,sumHBHE_adc3/totsumHBHE_adc);
	  pulseShape_adcHE_HPD->Fill(4.,sumHBHE_adc4/totsumHBHE_adc);
	  pulseShape_adcHE_HPD->Fill(5.,sumHBHE_adc5/totsumHBHE_adc);
	  pulseShape_adcHE_HPD->Fill(6.,sumHBHE_adc6/totsumHBHE_adc);
	  pulseShape_adcHE_HPD->Fill(7.,sumHBHE_adc7/totsumHBHE_adc);
	  pulseShape_adcHE_HPD->Fill(8.,sumHBHE_adc8/totsumHBHE_adc);
	  pulseShape_adcHE_HPD->Fill(9.,sumHBHE_adc9/totsumHBHE_adc);
	}//HE HPD

      }//HBHE Charge Restriction 
 
    }//HBHE Digi Loop

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
    double totsumHBHE_QIE11 = 0;

    double sumHBHE_QIE11_adc0 = 0;
    double sumHBHE_QIE11_adc1 = 0;
    double sumHBHE_QIE11_adc2 = 0;
    double sumHBHE_QIE11_adc3 = 0;
    double sumHBHE_QIE11_adc4 = 0;
    double sumHBHE_QIE11_adc5 = 0;
    double sumHBHE_QIE11_adc6 = 0;
    double sumHBHE_QIE11_adc7 = 0;
    double sumHBHE_QIE11_adc8 = 0;
    double sumHBHE_QIE11_adc9 = 0;
    double totsumHBHE_QIE11_adc = 0;

    for(int digi = 0; digi < (*DigiHBHE_QIE11_charge0).size(); digi++){
      if((*DigiHBHE_QIE11_charge0).size() > 0) sumHBHE_QIE110 = (*DigiHBHE_QIE11_charge0).at(digi);
      if((*DigiHBHE_QIE11_charge1).size() > 0) sumHBHE_QIE111 = (*DigiHBHE_QIE11_charge1).at(digi);
      if((*DigiHBHE_QIE11_charge2).size() > 0) sumHBHE_QIE112 = (*DigiHBHE_QIE11_charge2).at(digi);
      if((*DigiHBHE_QIE11_charge3).size() > 0) sumHBHE_QIE113 = (*DigiHBHE_QIE11_charge3).at(digi);
      if((*DigiHBHE_QIE11_charge4).size() > 0) sumHBHE_QIE114 = (*DigiHBHE_QIE11_charge4).at(digi);
      if((*DigiHBHE_QIE11_charge5).size() > 0) sumHBHE_QIE115 = (*DigiHBHE_QIE11_charge5).at(digi);
      if((*DigiHBHE_QIE11_charge6).size() > 0) sumHBHE_QIE116 = (*DigiHBHE_QIE11_charge6).at(digi);
      if((*DigiHBHE_QIE11_charge7).size() > 0) sumHBHE_QIE117 = (*DigiHBHE_QIE11_charge7).at(digi);
      if((*DigiHBHE_QIE11_charge8).size() > 0) sumHBHE_QIE118 = (*DigiHBHE_QIE11_charge8).at(digi);
      if((*DigiHBHE_QIE11_charge9).size() > 0) sumHBHE_QIE119 = (*DigiHBHE_QIE11_charge9).at(digi);
      totsumHBHE_QIE11 = sumHBHE_QIE114 + sumHBHE_QIE115;// + sumHBHE_QIE110 + sumHBHE_QIE111 + sumHBHE_QIE112 + sumHBHE_QIE113 + sumHBHE_QIE116 + sumHBHE_QIE117 + sumHBHE_QIE118 + sumHBHE_QIE119;

      if((*DigiHBHE_QIE11_adc0).size() > 0) sumHBHE_QIE11_adc0 = (*DigiHBHE_QIE11_adc0).at(digi);
      if((*DigiHBHE_QIE11_adc1).size() > 0) sumHBHE_QIE11_adc1 = (*DigiHBHE_QIE11_adc1).at(digi);
      if((*DigiHBHE_QIE11_adc2).size() > 0) sumHBHE_QIE11_adc2 = (*DigiHBHE_QIE11_adc2).at(digi);
      if((*DigiHBHE_QIE11_adc3).size() > 0) sumHBHE_QIE11_adc3 = (*DigiHBHE_QIE11_adc3).at(digi);
      if((*DigiHBHE_QIE11_adc4).size() > 0) sumHBHE_QIE11_adc4 = (*DigiHBHE_QIE11_adc4).at(digi);
      if((*DigiHBHE_QIE11_adc5).size() > 0) sumHBHE_QIE11_adc5 = (*DigiHBHE_QIE11_adc5).at(digi);
      if((*DigiHBHE_QIE11_adc6).size() > 0) sumHBHE_QIE11_adc6 = (*DigiHBHE_QIE11_adc6).at(digi);
      if((*DigiHBHE_QIE11_adc7).size() > 0) sumHBHE_QIE11_adc7 = (*DigiHBHE_QIE11_adc7).at(digi);
      if((*DigiHBHE_QIE11_adc8).size() > 0) sumHBHE_QIE11_adc8 = (*DigiHBHE_QIE11_adc8).at(digi);
      if((*DigiHBHE_QIE11_adc9).size() > 0) sumHBHE_QIE11_adc9 = (*DigiHBHE_QIE11_adc9).at(digi);
      totsumHBHE_QIE11_adc = sumHBHE_QIE11_adc4 + sumHBHE_QIE11_adc5;// + sumHBHE_QIE11_adc0 + sumHBHE_QIE11_adc1 + sumHBHE_QIE11_adc2 + sumHBHE_QIE11_adc3 + sumHBHE_QIE11_adc6 + sumHBHE_QIE11_adc7 + sumHBHE_QIE11_adc8 + sumHBHE_QIE11_adc9;

      ///Total charge must be greater than 5000 fC
      if(totsumHBHE_QIE11 > min_totcharge_SiPM){

	if((*DigiHBHE_QIE11_sub).at(digi) == 1){
	  pulseShape_DigiHB_SiPM->Fill(0.,sumHBHE_QIE110/totsumHBHE_QIE11);
	  pulseShape_DigiHB_SiPM->Fill(1.,sumHBHE_QIE111/totsumHBHE_QIE11);
	  pulseShape_DigiHB_SiPM->Fill(2.,sumHBHE_QIE112/totsumHBHE_QIE11);
	  pulseShape_DigiHB_SiPM->Fill(3.,sumHBHE_QIE113/totsumHBHE_QIE11);
	  pulseShape_DigiHB_SiPM->Fill(4.,sumHBHE_QIE114/totsumHBHE_QIE11);
	  pulseShape_DigiHB_SiPM->Fill(5.,sumHBHE_QIE115/totsumHBHE_QIE11);
	  pulseShape_DigiHB_SiPM->Fill(6.,sumHBHE_QIE116/totsumHBHE_QIE11);
	  pulseShape_DigiHB_SiPM->Fill(7.,sumHBHE_QIE117/totsumHBHE_QIE11);
	  pulseShape_DigiHB_SiPM->Fill(8.,sumHBHE_QIE118/totsumHBHE_QIE11);
	  pulseShape_DigiHB_SiPM->Fill(9.,sumHBHE_QIE119/totsumHBHE_QIE11);

	  pulseShape_adcHB_SiPM->Fill(0.,sumHBHE_QIE11_adc0/totsumHBHE_QIE11_adc);
	  pulseShape_adcHB_SiPM->Fill(1.,sumHBHE_QIE11_adc1/totsumHBHE_QIE11_adc);
	  pulseShape_adcHB_SiPM->Fill(2.,sumHBHE_QIE11_adc2/totsumHBHE_QIE11_adc);
	  pulseShape_adcHB_SiPM->Fill(3.,sumHBHE_QIE11_adc3/totsumHBHE_QIE11_adc);
	  pulseShape_adcHB_SiPM->Fill(4.,sumHBHE_QIE11_adc4/totsumHBHE_QIE11_adc);
	  pulseShape_adcHB_SiPM->Fill(5.,sumHBHE_QIE11_adc5/totsumHBHE_QIE11_adc);
	  pulseShape_adcHB_SiPM->Fill(6.,sumHBHE_QIE11_adc6/totsumHBHE_QIE11_adc);
	  pulseShape_adcHB_SiPM->Fill(7.,sumHBHE_QIE11_adc7/totsumHBHE_QIE11_adc);
	  pulseShape_adcHB_SiPM->Fill(8.,sumHBHE_QIE11_adc8/totsumHBHE_QIE11_adc);
	  pulseShape_adcHB_SiPM->Fill(9.,sumHBHE_QIE11_adc9/totsumHBHE_QIE11_adc);
	}//HB SiPM

	if((*DigiHBHE_QIE11_sub).at(digi) == 2){
	  pulseShape_DigiHE_SiPM->Fill(0.,sumHBHE_QIE110/totsumHBHE_QIE11);
	  pulseShape_DigiHE_SiPM->Fill(1.,sumHBHE_QIE111/totsumHBHE_QIE11);
	  pulseShape_DigiHE_SiPM->Fill(2.,sumHBHE_QIE112/totsumHBHE_QIE11);
	  pulseShape_DigiHE_SiPM->Fill(3.,sumHBHE_QIE113/totsumHBHE_QIE11);
	  pulseShape_DigiHE_SiPM->Fill(4.,sumHBHE_QIE114/totsumHBHE_QIE11);
	  pulseShape_DigiHE_SiPM->Fill(5.,sumHBHE_QIE115/totsumHBHE_QIE11);
	  pulseShape_DigiHE_SiPM->Fill(6.,sumHBHE_QIE116/totsumHBHE_QIE11);
	  pulseShape_DigiHE_SiPM->Fill(7.,sumHBHE_QIE117/totsumHBHE_QIE11);
	  pulseShape_DigiHE_SiPM->Fill(8.,sumHBHE_QIE118/totsumHBHE_QIE11);
	  pulseShape_DigiHE_SiPM->Fill(9.,sumHBHE_QIE119/totsumHBHE_QIE11);

	  pulseShape_adcHE_SiPM->Fill(0.,sumHBHE_QIE11_adc0/totsumHBHE_QIE11_adc);
	  pulseShape_adcHE_SiPM->Fill(1.,sumHBHE_QIE11_adc1/totsumHBHE_QIE11_adc);
	  pulseShape_adcHE_SiPM->Fill(2.,sumHBHE_QIE11_adc2/totsumHBHE_QIE11_adc);
	  pulseShape_adcHE_SiPM->Fill(3.,sumHBHE_QIE11_adc3/totsumHBHE_QIE11_adc);
	  pulseShape_adcHE_SiPM->Fill(4.,sumHBHE_QIE11_adc4/totsumHBHE_QIE11_adc);
	  pulseShape_adcHE_SiPM->Fill(5.,sumHBHE_QIE11_adc5/totsumHBHE_QIE11_adc);
	  pulseShape_adcHE_SiPM->Fill(6.,sumHBHE_QIE11_adc6/totsumHBHE_QIE11_adc);
	  pulseShape_adcHE_SiPM->Fill(7.,sumHBHE_QIE11_adc7/totsumHBHE_QIE11_adc);
	  pulseShape_adcHE_SiPM->Fill(8.,sumHBHE_QIE11_adc8/totsumHBHE_QIE11_adc);
	  pulseShape_adcHE_SiPM->Fill(9.,sumHBHE_QIE11_adc9/totsumHBHE_QIE11_adc);
	}//HE SiPM
	
      }//HBHE_QIE11 Charge Restriction 
 
    }//HBHE_QIE11 Digi Loop

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
    double totsumHO = 0;

    double sumHO_adc0 = 0;
    double sumHO_adc1 = 0;
    double sumHO_adc2 = 0;
    double sumHO_adc3 = 0;
    double sumHO_adc4 = 0;
    double sumHO_adc5 = 0;
    double sumHO_adc6 = 0;
    double sumHO_adc7 = 0;
    double sumHO_adc8 = 0;
    double sumHO_adc9 = 0;
    double totsumHO_adc = 0;

    for(int digi = 0; digi < (*DigiHO_charge0).size(); digi++){
      if((*DigiHO_charge0).size() > 0) sumHO0 = (*DigiHO_charge0).at(digi);
      if((*DigiHO_charge1).size() > 0) sumHO1 = (*DigiHO_charge1).at(digi);
      if((*DigiHO_charge2).size() > 0) sumHO2 = (*DigiHO_charge2).at(digi);
      if((*DigiHO_charge3).size() > 0) sumHO3 = (*DigiHO_charge3).at(digi);
      if((*DigiHO_charge4).size() > 0) sumHO4 = (*DigiHO_charge4).at(digi);
      if((*DigiHO_charge5).size() > 0) sumHO5 = (*DigiHO_charge5).at(digi);
      if((*DigiHO_charge6).size() > 0) sumHO6 = (*DigiHO_charge6).at(digi);
      if((*DigiHO_charge7).size() > 0) sumHO7 = (*DigiHO_charge7).at(digi);
      if((*DigiHO_charge8).size() > 0) sumHO8 = (*DigiHO_charge8).at(digi);
      if((*DigiHO_charge9).size() > 0) sumHO9 = (*DigiHO_charge9).at(digi);
      totsumHO = sumHO4 + sumHO5;//+ sumHO0 + sumHO1 + sumHO2 + sumHO3 + sumHO6 + sumHO7 + sumHO8 + sumHO9;

      if((*DigiHO_adc0).size() > 0) sumHO_adc0 = (*DigiHO_adc0).at(digi);
      if((*DigiHO_adc1).size() > 0) sumHO_adc1 = (*DigiHO_adc1).at(digi);
      if((*DigiHO_adc2).size() > 0) sumHO_adc2 = (*DigiHO_adc2).at(digi);
      if((*DigiHO_adc3).size() > 0) sumHO_adc3 = (*DigiHO_adc3).at(digi);
      if((*DigiHO_adc4).size() > 0) sumHO_adc4 = (*DigiHO_adc4).at(digi);
      if((*DigiHO_adc5).size() > 0) sumHO_adc5 = (*DigiHO_adc5).at(digi);
      if((*DigiHO_adc6).size() > 0) sumHO_adc6 = (*DigiHO_adc6).at(digi);
      if((*DigiHO_adc7).size() > 0) sumHO_adc7 = (*DigiHO_adc7).at(digi);
      if((*DigiHO_adc8).size() > 0) sumHO_adc8 = (*DigiHO_adc8).at(digi);
      if((*DigiHO_adc9).size() > 0) sumHO_adc9 = (*DigiHO_adc9).at(digi);
      totsumHO_adc = sumHO_adc4 + sumHO_adc5;//+ sumHO_adc0 + sumHO_adc1 + sumHO_adc2 + sumHO_adc3 + sumHO_adc6 + sumHO_adc7 + sumHO_adc8 + sumHO_adc9;

      ///Total charge must be greater than 5000 fC
      if(totsumHO > min_totcharge_SiPM){
	pulseShape_DigiHO->Fill(0.,sumHO0/totsumHO);
	pulseShape_DigiHO->Fill(1.,sumHO1/totsumHO);
	pulseShape_DigiHO->Fill(2.,sumHO2/totsumHO);
	pulseShape_DigiHO->Fill(3.,sumHO3/totsumHO);
	pulseShape_DigiHO->Fill(4.,sumHO4/totsumHO);
	pulseShape_DigiHO->Fill(5.,sumHO5/totsumHO);
	pulseShape_DigiHO->Fill(6.,sumHO6/totsumHO);
	pulseShape_DigiHO->Fill(7.,sumHO7/totsumHO);
	pulseShape_DigiHO->Fill(8.,sumHO8/totsumHO);
	pulseShape_DigiHO->Fill(9.,sumHO9/totsumHO);

	pulseShape_adcHO->Fill(0.,sumHO_adc0/totsumHO_adc);
	pulseShape_adcHO->Fill(1.,sumHO_adc1/totsumHO_adc);
	pulseShape_adcHO->Fill(2.,sumHO_adc2/totsumHO_adc);
	pulseShape_adcHO->Fill(3.,sumHO_adc3/totsumHO_adc);
	pulseShape_adcHO->Fill(4.,sumHO_adc4/totsumHO_adc);
	pulseShape_adcHO->Fill(5.,sumHO_adc5/totsumHO_adc);
	pulseShape_adcHO->Fill(6.,sumHO_adc6/totsumHO_adc);
	pulseShape_adcHO->Fill(7.,sumHO_adc7/totsumHO_adc);
	pulseShape_adcHO->Fill(8.,sumHO_adc8/totsumHO_adc);
	pulseShape_adcHO->Fill(9.,sumHO_adc9/totsumHO_adc);
      }//HO Charge Restriction 
 
    }//HO Digi Loop

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
    double totsumHF_QIE10 = 0;

    double sumHF_QIE10_adc0 = 0;
    double sumHF_QIE10_adc1 = 0;
    double sumHF_QIE10_adc2 = 0;
    double sumHF_QIE10_adc3 = 0;
    double sumHF_QIE10_adc4 = 0;
    double sumHF_QIE10_adc5 = 0;
    double sumHF_QIE10_adc6 = 0;
    double sumHF_QIE10_adc7 = 0;
    double sumHF_QIE10_adc8 = 0;
    double sumHF_QIE10_adc9 = 0;
    double totsumHF_QIE10_adc = 0;

    //cout << "(*DigiHF_QIE10_charge0).size() = " <<  (*DigiHF_QIE10_charge0).size() << endl;
    //cout << "(*DigiHF_QIE10_charge4).size() = "  << (*DigiHF_QIE10_charge4).size() << endl;

    for(int digi = 0; digi < (*DigiHF_QIE10_charge0).size(); digi++){
      if((*DigiHF_QIE10_charge0).size() > 0) sumHF_QIE100 = (*DigiHF_QIE10_charge0).at(digi);
      if((*DigiHF_QIE10_charge1).size() > 0) sumHF_QIE101 = (*DigiHF_QIE10_charge1).at(digi);
      if((*DigiHF_QIE10_charge2).size() > 0) sumHF_QIE102 = (*DigiHF_QIE10_charge2).at(digi);
      if((*DigiHF_QIE10_charge3).size() > 0) sumHF_QIE103 = (*DigiHF_QIE10_charge3).at(digi);
      if((*DigiHF_QIE10_charge4).size() > 0) sumHF_QIE104 = (*DigiHF_QIE10_charge4).at(digi);
      if((*DigiHF_QIE10_charge5).size() > 0) sumHF_QIE105 = (*DigiHF_QIE10_charge5).at(digi);
      if((*DigiHF_QIE10_charge6).size() > 0) sumHF_QIE106 = (*DigiHF_QIE10_charge6).at(digi);
      if((*DigiHF_QIE10_charge7).size() > 0) sumHF_QIE107 = (*DigiHF_QIE10_charge7).at(digi);
      if((*DigiHF_QIE10_charge8).size() > 0) sumHF_QIE108 = (*DigiHF_QIE10_charge8).at(digi);
      if((*DigiHF_QIE10_charge9).size() > 0) sumHF_QIE109 = (*DigiHF_QIE10_charge9).at(digi);
      totsumHF_QIE10 = sumHF_QIE104 + sumHF_QIE105 + sumHF_QIE100 + sumHF_QIE101 + sumHF_QIE102 + sumHF_QIE103 + sumHF_QIE106 + sumHF_QIE107 + sumHF_QIE108 + sumHF_QIE109;

      if((*DigiHF_QIE10_adc0).size() > 0) sumHF_QIE10_adc0 = (*DigiHF_QIE10_adc0).at(digi);
      if((*DigiHF_QIE10_adc1).size() > 0) sumHF_QIE10_adc1 = (*DigiHF_QIE10_adc1).at(digi);
      if((*DigiHF_QIE10_adc2).size() > 0) sumHF_QIE10_adc2 = (*DigiHF_QIE10_adc2).at(digi);
      if((*DigiHF_QIE10_adc3).size() > 0) sumHF_QIE10_adc3 = (*DigiHF_QIE10_adc3).at(digi);
      if((*DigiHF_QIE10_adc4).size() > 0) sumHF_QIE10_adc4 = (*DigiHF_QIE10_adc4).at(digi);
      if((*DigiHF_QIE10_adc5).size() > 0) sumHF_QIE10_adc5 = (*DigiHF_QIE10_adc5).at(digi);
      if((*DigiHF_QIE10_adc6).size() > 0) sumHF_QIE10_adc6 = (*DigiHF_QIE10_adc6).at(digi);
      if((*DigiHF_QIE10_adc7).size() > 0) sumHF_QIE10_adc7 = (*DigiHF_QIE10_adc7).at(digi);
      if((*DigiHF_QIE10_adc8).size() > 0) sumHF_QIE10_adc8 = (*DigiHF_QIE10_adc8).at(digi);
      if((*DigiHF_QIE10_adc9).size() > 0) sumHF_QIE10_adc9 = (*DigiHF_QIE10_adc9).at(digi);
      totsumHF_QIE10_adc = sumHF_QIE10_adc4 + sumHF_QIE10_adc5 + sumHF_QIE10_adc0 + sumHF_QIE10_adc1 + sumHF_QIE10_adc2 + sumHF_QIE10_adc3 + sumHF_QIE10_adc6 + sumHF_QIE10_adc7 + sumHF_QIE10_adc8 + sumHF_QIE10_adc9;

      ///Total charge must be greater than 5000 fC
      if(totsumHF_QIE10 > min_totcharge_HF){
	pulseShape_DigiHF_QIE10->Fill(0.,sumHF_QIE100/totsumHF_QIE10);
	pulseShape_DigiHF_QIE10->Fill(1.,sumHF_QIE101/totsumHF_QIE10);
	pulseShape_DigiHF_QIE10->Fill(2.,sumHF_QIE102/totsumHF_QIE10);
	pulseShape_DigiHF_QIE10->Fill(3.,sumHF_QIE103/totsumHF_QIE10);
	pulseShape_DigiHF_QIE10->Fill(4.,sumHF_QIE104/totsumHF_QIE10);
	pulseShape_DigiHF_QIE10->Fill(5.,sumHF_QIE105/totsumHF_QIE10);
	pulseShape_DigiHF_QIE10->Fill(6.,sumHF_QIE106/totsumHF_QIE10);
	pulseShape_DigiHF_QIE10->Fill(7.,sumHF_QIE107/totsumHF_QIE10);
	pulseShape_DigiHF_QIE10->Fill(8.,sumHF_QIE108/totsumHF_QIE10);
	pulseShape_DigiHF_QIE10->Fill(9.,sumHF_QIE109/totsumHF_QIE10);

	pulseShape_adcHF_QIE10->Fill(0.,sumHF_QIE10_adc0/totsumHF_QIE10_adc);
	pulseShape_adcHF_QIE10->Fill(1.,sumHF_QIE10_adc1/totsumHF_QIE10_adc);
	pulseShape_adcHF_QIE10->Fill(2.,sumHF_QIE10_adc2/totsumHF_QIE10_adc);
	pulseShape_adcHF_QIE10->Fill(3.,sumHF_QIE10_adc3/totsumHF_QIE10_adc);
	pulseShape_adcHF_QIE10->Fill(4.,sumHF_QIE10_adc4/totsumHF_QIE10_adc);
	pulseShape_adcHF_QIE10->Fill(5.,sumHF_QIE10_adc5/totsumHF_QIE10_adc);
	pulseShape_adcHF_QIE10->Fill(6.,sumHF_QIE10_adc6/totsumHF_QIE10_adc);
	pulseShape_adcHF_QIE10->Fill(7.,sumHF_QIE10_adc7/totsumHF_QIE10_adc);
	pulseShape_adcHF_QIE10->Fill(8.,sumHF_QIE10_adc8/totsumHF_QIE10_adc);
	pulseShape_adcHF_QIE10->Fill(9.,sumHF_QIE10_adc9/totsumHF_QIE10_adc);
      }//HF_QIE10 Charge Restriction 
 
    }//HF_QIE10 Digi Loop

  }//End of While Loop

  //Fill Histograms
  f2->Write();

}//End of Macro

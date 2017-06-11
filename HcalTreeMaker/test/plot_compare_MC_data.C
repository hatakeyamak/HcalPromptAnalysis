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
#include "TArrayD.h"

void plot_compare_MC_data(const char *Data_File, const char *MC_File, const char *run_number = "<Run Number>", const char *MC_Type = "<MC Type>" ){

  ////////////////////////////////////////////////////////////////////
  //                 Define everything that is needed               //
  ////////////////////////////////////////////////////////////////////

  //Set the plots max and min y-value
  Double_t ymin =0.;
  Double_t ymax =0.75;

  Double_t ymin_SiPM = 0.;
  Double_t ymax_SiPM = 0.6;

  Double_t ymin_HF =0.;
  Double_t ymax_HF =1;

  //Define the canvas 
  TCanvas *c0 = new TCanvas("PulseShape_HB","",1000,800);

  //Define/Get the Histos to be plotted
  TFile *f1 = new TFile(Data_File);
  if(f1->IsZombie()){
    cout << "Root file: " << Data_File << " not found!" << endl;
    return;
  }
  TFile *f2 = new TFile(MC_File);
  if(f2->IsZombie()){
    cout << "Root file: " << MC_File << " not found!" << endl;
    return;
  }

  TH1F *h1 = (TH1F*)f1->Get("pulseShape_DigiHB_HPD");
  TH1F *h2 = (TH1F*)f2->Get("pulseShape_DigiHB_HPD");
  TH1F *h3 = (TH1F*)f1->Get("pulseShape_DigiHE_HPD");
  TH1F *h4 = (TH1F*)f2->Get("pulseShape_DigiHE_HPD");
  TH1F *h5 = (TH1F*)f1->Get("pulseShape_DigiHE_SiPM");
  TH1F *h6 = (TH1F*)f2->Get("pulseShape_DigiHE_SiPM");
  TH1F *h7 = (TH1F*)f1->Get("pulseShape_DigiHF_QIE10");
  TH1F *h8 = (TH1F*)f2->Get("pulseShape_DigiHF_QIE10");

  //Define the Legends
  TLegend* catLeg0 = new TLegend(0.55,0.7,0.9,0.9);
  TLegend* catLeg1 = new TLegend(0.55,0.7,0.9,0.9);
  TLegend* catLeg2 = new TLegend(0.55,0.7,0.9,0.9);
  TLegend* catLeg3 = new TLegend(0.55,0.7,0.9,0.9);
  //TLegend* catLeg1 = new TLegend(0.7,0.8,0.99,1.0);
  //TLegend* catLeg2 = new TLegend(0.7,0.8,0.99,1.0);
  //TLegend* catLeg3 = new TLegend(0.7,0.8,0.99,1.0);

  //Divide the Canvas as needed
  c0->Divide(2,2);

  ////////////////////////////////////////////////////////////////////
  //                   Plotting Data vs MC for  HB HPD              //
  ////////////////////////////////////////////////////////////////////

  c0->cd(1);

  //Data
  h1->SetLineColor(1);
  h1->Draw("hist"); 
  h1->SetMinimum(ymin);
  h1->SetMaximum(ymax);
  h1->GetXaxis()->SetTitle("Time Step");
  h1->GetYaxis()->SetTitle("A.U.");
  h1->SetTitleOffset(1.5,"Y");
  h1->SetTitle("Data vs MC HB_HPD");
  h1->SetName("");
  h1->SetStats(false);

  catLeg0->SetTextSize(0.03);
  catLeg0->AddEntry(h1, run_number ,"l");  
  catLeg0->Draw();

  //MC
  h2->SetLineColor(2);
  h2->Draw("hist same");
  catLeg0->AddEntry(h2,MC_Type,"l");

  ////////////////////////////////////////////////////////////////////
  //                   Plotting Data vs MC for  HE HPD              //
  ////////////////////////////////////////////////////////////////////
  
  c0->cd(2);
  
  //Data
  h3->SetLineColor(1);
  h3->Draw("hist"); 
  h3->SetMinimum(ymin);
  h3->SetMaximum(ymax);
  h3->GetXaxis()->SetTitle("Time Step");
  h3->GetYaxis()->SetTitle("A.U.");
  h3->SetTitleOffset(1.5,"Y");
  h3->SetTitle("Data vs MC HE_HPD");
  h3->SetName("");
  h3->SetStats(false);

  catLeg1->SetTextSize(0.03);
  catLeg1->AddEntry(h3,run_number,"l");  
  catLeg1->Draw();

  //MC
  h4->SetLineColor(2);
  h4->Draw("hist same");
  catLeg1->AddEntry(h4,MC_Type,"l");

  ////////////////////////////////////////////////////////////////////
  //                   Plotting Data vs MC for  HE SiPM              //
  ////////////////////////////////////////////////////////////////////

  c0->cd(3);
  
  //Data
  h5->SetLineColor(1);
  h5->Draw("hist"); 
  h5->SetMinimum(ymin_SiPM);
  h5->SetMaximum(ymax_SiPM);
  h5->GetXaxis()->SetTitle("Time Step");
  h5->GetYaxis()->SetTitle("A.U.");
  h5->SetTitleOffset(1.5,"Y");
  h5->SetTitle("Data vs MC HE_SiPM");
  h5->SetName("");
  h5->SetStats(false);

  catLeg2->SetTextSize(0.03);
  catLeg2->AddEntry(h5,run_number,"l");  
  catLeg2->Draw();

  //MC
  h6->SetLineColor(2);
  h6->Draw("hist same");
  catLeg2->AddEntry(h6,MC_Type,"l");

  ////////////////////////////////////////////////////////////////////
  //                   Plotting Data vs MC for  Hf_QIE10            //
  ////////////////////////////////////////////////////////////////////

  c0->cd(4);
  
  //Data
  h7->SetLineColor(1);
  h7->Draw("hist"); 
  h7->SetMinimum(ymin_HF);
  h7->SetMaximum(ymax_HF);
  h7->GetXaxis()->SetTitle("Time Step");
  h7->GetYaxis()->SetTitle("A.U.");
  h7->SetTitleOffset(1.5,"Y");
  h7->SetTitle("Data vs MC HF_QIE10");
  h7->SetName("");
  h7->SetStats(false);

  catLeg3->SetTextSize(0.03);
  catLeg3->AddEntry(h7,run_number,"l");  
  catLeg3->Draw();

  //MC
  h8->SetLineColor(2);
  //h8->Draw("hist same");
  //catLeg3->AddEntry(h8,MC_Type,"l");

  //Shift HF MC 
  TH1F *hnew = new TH1F("hnew2","hnew2",10,-0.5,9.5);
  int ishift=-1; //to be tuned.
  for(int j=0; j<=9; j++){
    hnew->SetBinContent(j+ishift,h8->GetBinContent(j));
  }

  ////////////////////////////////////////////////////////////////////
  //                Shifting HF, Normalizing, and Saving            //
  ////////////////////////////////////////////////////////////////////

  hnew->SetLineColor(4);
  hnew->Draw("hist same");
  catLeg3->AddEntry(hnew,"Shifted MC","l");
  h7->Draw("hist same");

  Double_t norm = 1;
  Double_t min  = -0.5;
  Double_t max  = 9.5;
  h1->Scale(norm/h1->Integral(min,max));
  h2->Scale(norm/h2->Integral(min,max));
  h3->Scale(norm/h3->Integral(min,max));
  h4->Scale(norm/h4->Integral(min,max));
  h5->Scale(norm/h5->Integral(min,max));
  h6->Scale(norm/h6->Integral(min,max));
  h7->Scale(norm/h7->Integral(min,max));
  h8->Scale(norm/h8->Integral(min,max));
  hnew->Scale(norm/hnew->Integral(min,max));

  //Save Plot
  //c0->SaveAs("Data_vs_MC.gif");
  //c0->SaveAs("Data_vs_MC.pdf");

}//End of Macro

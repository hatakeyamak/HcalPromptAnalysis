{

  gStyle->SetOptStat(0);
  
  TFile *_file0 = TFile::Open("hcal_noisestudy_histograms_age_new2_4500ultimate.root");
  _file0->ls();

  // - MC
  TH1F *FC_2p8mm_SiPM    = (TH1F*) _file0->Get("HcalDigiTask_FC_2p8mmSiPM_HB");
  TH1F *RawFC_2p8mm_SiPM = (TH1F*) _file0->Get("HcalDigiTask_RawFC_2p8mmSiPM_HB");
  TH1F *FC_3p3mm_SiPM    = (TH1F*) _file0->Get("HcalDigiTask_FC_3p3mmSiPM_HB");
  TH1F *RawFC_3p3mm_SiPM = (TH1F*) _file0->Get("HcalDigiTask_RawFC_3p3mmSiPM_HB");

  //----------
  
  TCanvas *c5 = new TCanvas("c5","c5",800,600);
  gPad->SetLogy();
  gPad->SetGrid(1);

  RawFC_2p8mm_SiPM->SetTitleSize(0.05,"XY");
  RawFC_2p8mm_SiPM->GetYaxis()->SetTitle("# of digis");
  RawFC_2p8mm_SiPM->GetXaxis()->SetTitle("Q(SOI)+Q(SOI+1) [fC]");
  RawFC_2p8mm_SiPM->GetYaxis()->SetTitleOffset(0.8);
  RawFC_2p8mm_SiPM->SetLabelSize(.03,"XY");
  RawFC_2p8mm_SiPM->SetMinimum(1.0);
  RawFC_2p8mm_SiPM->SetTitle("2.8mm SiPM [depth 1, 2, 3 (|ieta|=16)]");
  RawFC_2p8mm_SiPM->Draw();
  FC_2p8mm_SiPM->SetLineColor(2);
  FC_2p8mm_SiPM->SetLineStyle(2);
  FC_2p8mm_SiPM->Draw("same");
  TLegend *tl5 = new TLegend(0.7,0.7,0.9,0.9);
  tl5->AddEntry(RawFC_2p8mm_SiPM,"Raw Charge");
  tl5->AddEntry(FC_2p8mm_SiPM,"Pedestal subtracted");
  tl5->Draw();

  c5->Print("plots/plot_FC_2p8mmSiPM.png");
  c5->Print("plots/plot_FC_2p8mmSiPM.pdf");
  
  //----------
  
  TCanvas *c6 = new TCanvas("c6","c6",800,600);
  gPad->SetLogy();
  gPad->SetGrid(1);

  RawFC_3p3mm_SiPM->SetTitleSize(0.05,"XY");
  RawFC_3p3mm_SiPM->GetYaxis()->SetTitle("# of digis");
  RawFC_3p3mm_SiPM->GetXaxis()->SetTitle("Q(SOI)+Q(SOI+1) [fC]");
  RawFC_3p3mm_SiPM->GetYaxis()->SetTitleOffset(0.8);
  RawFC_3p3mm_SiPM->SetLabelSize(.03,"XY");
  RawFC_3p3mm_SiPM->SetMinimum(1.0);
  RawFC_3p3mm_SiPM->SetTitle("3.3mm SiPM [depth 3 (|ieta|!=16), depth 4]");
  RawFC_3p3mm_SiPM->Draw();
  FC_3p3mm_SiPM->SetLineColor(2);
  FC_3p3mm_SiPM->SetLineStyle(2);
  FC_3p3mm_SiPM->Draw("same");
  TLegend *tl6 = new TLegend(0.7,0.7,0.9,0.9);
  tl6->AddEntry(RawFC_3p3mm_SiPM,"Raw Charge");
  tl6->AddEntry(FC_3p3mm_SiPM,"Pedestal subtracted");
  tl6->Draw();

  c6->Print("plots/plot_FC_3p3mmSiPM.png");
  c6->Print("plots/plot_FC_3p3mmSiPM.pdf");

  

  
}

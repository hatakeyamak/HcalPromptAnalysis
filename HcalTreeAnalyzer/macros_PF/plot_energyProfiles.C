{

  gStyle->SetOptStat(0);
  
  TFile *_file0 = TFile::Open("pfstudy_histograms.root");
  _file0->ls();

  // - MC
  TH1F *HcalProfile_NeutralHadron_Endcap_PtAbove5 = (TH1F*) _file0->Get("PFTask_hcalProfile_NeutralHadron_Endcap_PtAbove5");
  TH1F *HcalProfile_NeutralHadron_Endcap_Pt1To5   = (TH1F*) _file0->Get("PFTask_hcalProfile_NeutralHadron_Endcap_Pt1To5");
  TH1F *HcalProfile_NeutralHadron_Endcap_PtBelow1 = (TH1F*) _file0->Get("PFTask_hcalProfile_NeutralHadron_Endcap_PtBelow1");

  std::cout << "aaa" << std::endl;
  TH1F *HcalProfile_ChargedHadron_Endcap_PtAbove5 = (TH1F*) _file0->Get("PFTask_Profile_ChargedHadron_Endcap_PtAbove5");
  TH1F *HcalProfile_ChargedHadron_Endcap_Pt1To5   = (TH1F*) _file0->Get("PFTask_Profile_ChargedHadron_Endcap_Pt1To5");
  TH1F *HcalProfile_ChargedHadron_Endcap_PtBelow1 = (TH1F*) _file0->Get("PFTask_Profile_ChargedHadron_Endcap_PtBelow1");
  std::cout << "bbb" << std::endl;
  
  //----------
  
  TCanvas *c5 = new TCanvas("c5","c5",1200,400);
  c5->Divide(3,1);
  gPad->SetLogy();
  gPad->SetGrid(1);

  TH1F *tframe = new TH1F("tframe","tframe",7,0.5,7.5);
  tframe->SetMinimum(0.);
  tframe->SetMaximum(1.3);
  tframe->SetTitleSize(0.06,"XY");
  tframe->GetYaxis()->SetTitle("Fraction");
  tframe->GetXaxis()->SetTitle("HE Depth");
  tframe->GetYaxis()->SetTitleOffset(0.8);
  tframe->GetXaxis()->SetTitleOffset(0.8);
  tframe->SetLabelSize(.04,"XY");

  c5->cd(1);
  gPad->SetGrid();
  tframe->SetTitle("Neutral hadron, pt>5 GeV, 1.5<|#eta|<2.9");
  tframe->DrawCopy();
  HcalProfile_NeutralHadron_Endcap_PtAbove5->Draw("same");
  /*
  TLegend *tl5 = new TLegend(0.5,0.7,0.9,0.9);
  tl5->AddEntry(HcalFrac1Zero,"t#bar{t} PU 2018");
  tl5->AddEntry(HcalFrac1Zero_2019,"t#bar{t} noPU 2019");
  tl5->Draw();
  */

  c5->cd(2);
  gPad->SetGrid();
  tframe->SetTitle("Neutral hadron, 1<pt<5 GeV, 1.5<|#eta|<2.9");
  tframe->DrawCopy();
  HcalProfile_NeutralHadron_Endcap_Pt1To5->Draw("same");
  /*
  TLegend *tl5 = new TLegend(0.5,0.7,0.9,0.9);
  tl5->AddEntry(HcalFrac1Zero,"t#bar{t} PU 2018");
  tl5->AddEntry(HcalFrac1Zero_2019,"t#bar{t} noPU 2019");
  tl5->Draw();
  */

  c5->cd(3);
  gPad->SetGrid();
  tframe->SetTitle("Neutral hadron, pt<1 GeV, 1.5<|#eta|<2.9");
  tframe->DrawCopy();
  HcalProfile_NeutralHadron_Endcap_PtBelow1->Draw("same");
  /*
  TLegend *tl5 = new TLegend(0.5,0.7,0.9,0.9);
  tl5->AddEntry(HcalFrac1Zero,"t#bar{t} PU 2018");
  tl5->AddEntry(HcalFrac1Zero_2019,"t#bar{t} noPU 2019");
  tl5->Draw();
  */
    
  /*
  FC_2p8mm_SiPM->SetLineColor(2);
  FC_2p8mm_SiPM->SetLineStyle(2);
  FC_2p8mm_SiPM->Draw("same");
  TLegend *tl5 = new TLegend(0.7,0.7,0.9,0.9);
  tl5->AddEntry(RawFC_2p8mm_SiPM,"Raw Charge");
  tl5->AddEntry(FC_2p8mm_SiPM,"Pedestal subtracted");
  tl5->Draw();
  */

  c5->Print("plots/plot_PFCand_HcalEnergyProfile.png");
  c5->Print("plots/plot_PFCand_HcalEnergyProfile.pdf");
  
  //----------
  
  TCanvas *c6 = new TCanvas("c6","c6",1200,400);
  c6->Divide(3,1);
  gPad->SetLogy();
  gPad->SetGrid(1);

  TH1F *tframe2 = new TH1F("tframe2","tframe2",11,-3.5,7.5);
  tframe2->SetMinimum(0.);
  tframe2->SetMaximum(1.3);
  tframe2->SetTitleSize(0.06,"XY");
  tframe2->GetYaxis()->SetTitle("Fraction");
  tframe2->GetXaxis()->SetTitle("");
  tframe2->GetYaxis()->SetTitleOffset(0.8);
  tframe2->GetXaxis()->SetTitleOffset(0.8);
  tframe2->SetLabelSize(.04,"XY");
  tframe2->GetXaxis()->SetBinLabel(1,"Track");
  tframe2->GetXaxis()->SetBinLabel(2,"ECAL");
  tframe2->GetXaxis()->SetBinLabel(3,"HCAL");
  tframe2->GetXaxis()->SetBinLabel(5,"HCAL d1");
  tframe2->GetXaxis()->SetBinLabel(6,"d2");
  tframe2->GetXaxis()->SetBinLabel(7,"d3");
  tframe2->GetXaxis()->SetBinLabel(8,"d4");
  tframe2->GetXaxis()->SetBinLabel(9,"d5");
  tframe2->GetXaxis()->SetBinLabel(10,"d6");
  tframe2->GetXaxis()->SetBinLabel(11,"d7");
  
  c6->cd(1);
  gPad->SetGrid();
  tframe2->SetTitle("Charged hadron, pt>5 GeV, 1.5<|#eta|<2.5");
  tframe2->DrawCopy();
  HcalProfile_ChargedHadron_Endcap_PtAbove5->Draw("same");
  /*
  TLegend *tl5 = new TLegend(0.5,0.7,0.9,0.9);
  tl5->AddEntry(HcalFrac1Zero,"t#bar{t} PU 2018");
  tl5->AddEntry(HcalFrac1Zero_2019,"t#bar{t} noPU 2019");
  tl5->Draw();
  */

  c6->cd(2);
  gPad->SetGrid();
  tframe2->SetTitle("Charged hadron, 1<pt<5 GeV, 1.5<|#eta|<2.5");
  tframe2->DrawCopy();
  HcalProfile_ChargedHadron_Endcap_Pt1To5->Draw("same");
  /*
  TLegend *tl5 = new TLegend(0.5,0.7,0.9,0.9);
  tl5->AddEntry(HcalFrac1Zero,"t#bar{t} PU 2018");
  tl5->AddEntry(HcalFrac1Zero_2019,"t#bar{t} noPU 2019");
  tl5->Draw();
  */

  c6->cd(3);
  gPad->SetGrid();
  tframe2->SetTitle("Charged hadron, pt<1 GeV, 1.5<|#eta|<2.5");
  tframe2->DrawCopy();
  HcalProfile_ChargedHadron_Endcap_PtBelow1->Draw("same");
  /*
  TLegend *tl5 = new TLegend(0.5,0.7,0.9,0.9);
  tl5->AddEntry(HcalFrac1Zero,"t#bar{t} PU 2018");
  tl5->AddEntry(HcalFrac1Zero_2019,"t#bar{t} noPU 2019");
  tl5->Draw();
  */
    
  /*
  FC_2p8mm_SiPM->SetLineColor(2);
  FC_2p8mm_SiPM->SetLineStyle(2);
  FC_2p8mm_SiPM->Draw("same");
  TLegend *tl5 = new TLegend(0.7,0.7,0.9,0.9);
  tl5->AddEntry(RawFC_2p8mm_SiPM,"Raw Charge");
  tl5->AddEntry(FC_2p8mm_SiPM,"Pedestal subtracted");
  tl5->Draw();
  */

  c6->Print("plots/plot_PFCandChargedHadron_EnergyProfile.png");
  c6->Print("plots/plot_PFCandChargedHadron_EnergyProfile.pdf");
  
  //----------
  
}

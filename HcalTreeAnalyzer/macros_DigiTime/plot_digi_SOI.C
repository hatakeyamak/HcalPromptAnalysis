{

  gStyle->SetOptStat(0);
  
  TFile *_file0 = TFile::Open("hcal_histograms_pt50.root");
  TFile *_file1 = TFile::Open("hcal_histograms_pt50.root");
  _file0->ls();

  // - MC
  TH1F *SOIfrac_HB = (TH1F*) _file0->Get("HcalDigiTask_SOI_frac_HB");
  TH1F *SOIfrac_HE = (TH1F*) _file0->Get("HcalDigiTask_SOI_frac_HE");
  //TH1F *SOIfrac_HF = (TH1F*) _file0->Get("HcalDigiTask_SOI_frac_HF");
  //TH1F *SOIfrac_HO = (TH1F*) _file0->Get("HcalDigiTask_SOI_frac_HO");

  TH1F *PostSOIfrac_HB = (TH1F*) _file0->Get("HcalDigiTask_postSOI_frac_HB");
  TH1F *PostSOIfrac_HE = (TH1F*) _file0->Get("HcalDigiTask_postSOI_frac_HE");
  //TH1F *PostSOIfrac_HF = (TH1F*) _file0->Get("HcalDigiTask_postSOI_frac_HF");
  //TH1F *PostSOIfrac_HO = (TH1F*) _file0->Get("HcalDigiTask_postSOI_frac_HO");

  TH1F *charge_fail = (TH1F*) _file0->Get("HcalDigiTask_Charge_Delayed_HE");
  TH1F *charge_pass = (TH1F*) _file0->Get("HcalDigiTask_Charge_Prompt_HE");

  TH1F *SOIfrac_pass_HE = (TH1F*) _file0->Get("HcalDigiTask_SOI_frac_pass_HE");
  TH1F *SOIfrac_fail_HE = (TH1F*) _file0->Get("HcalDigiTask_SOI_frac_fail_HE");

  TH1F *SOIfrac_1000_HE = (TH1F*) _file0->Get("HcalDigiTask_SOI_frac_1000_HE");
  TH1F *SOIfrac_2000_HE = (TH1F*) _file0->Get("HcalDigiTask_SOI_frac_2000_HE");
  TH1F *SOIfrac_6000_HE = (TH1F*) _file0->Get("HcalDigiTask_SOI_frac_6000_HE");

  TH1F *SOIfrac_HE_normalized      = (TH1F*) SOIfrac_HE->Clone("SOIfrac_HE_normalized");
  TH1F *SOIfrac_1000_HE_normalized = (TH1F*) SOIfrac_1000_HE->Clone("SOIfrac_1000_HE_normalized");
  TH1F *SOIfrac_2000_HE_normalized = (TH1F*) SOIfrac_2000_HE->Clone("SOIfrac_2000_HE_normalized");
  TH1F *SOIfrac_6000_HE_normalized = (TH1F*) SOIfrac_6000_HE->Clone("SOIfrac_6000_HE_normalized");
  SOIfrac_HE_normalized->Scale(1./SOIfrac_HE->GetSum());
  SOIfrac_1000_HE_normalized->Scale(1./SOIfrac_1000_HE->GetSum());
  SOIfrac_2000_HE_normalized->Scale(1./SOIfrac_2000_HE->GetSum());
  SOIfrac_6000_HE_normalized->Scale(1./SOIfrac_6000_HE->GetSum());
  
  TH1F *simhit_avetime_prompt = (TH1F*) _file0->Get("Simhit_AveTime_PromptHits_HE");
  TH1F *simhit_avetime_delayed = (TH1F*) _file0->Get("Simhit_AveTime_DelayedHits_HE");

  // - Data
  TH1F *SOIfrac_HB_data = (TH1F*) _file1->Get("HcalDigiTask_SOI_frac_HB");
  TH1F *SOIfrac_HE_data = (TH1F*) _file1->Get("HcalDigiTask_SOI_frac_HE");
  //TH1F *SOIfrac_HF_data = (TH1F*) _file1->Get("HcalDigiTask_SOI_frac_HF");
  //TH1F *SOIfrac_HO_data = (TH1F*) _file1->Get("HcalDigiTask_SOI_frac_HO");

  TH1F *PostSOIfrac_HB_data = (TH1F*) _file1->Get("HcalDigiTask_postSOI_frac_HB");
  TH1F *PostSOIfrac_HE_data = (TH1F*) _file1->Get("HcalDigiTask_postSOI_frac_HE");
  //TH1F *PostSOIfrac_HF_data = (TH1F*) _file1->Get("HcalDigiTask_postSOI_frac_HF");
  //TH1F *PostSOIfrac_HO_data = (TH1F*) _file1->Get("HcalDigiTask_postSOI_frac_HO");

  TH1F *charge_fail_data = (TH1F*) _file1->Get("HcalDigiTask_Charge_Delayed_HE");
  TH1F *charge_pass_data = (TH1F*) _file1->Get("HcalDigiTask_Charge_Prompt_HE");

  //---------- 

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  TPad *c1_1 = new TPad("c1_1", "c1_1",0.0,0.32,1.0,1.0);
  c1_1->Draw();
  c1_1->cd();

  gPad->SetLogy();
  gPad->SetGrid(1);

  charge_pass->SetTitleSize(0.07,"XY");
  charge_pass->GetYaxis()->SetTitle("Number of hits");
  charge_pass->GetXaxis()->SetTitle("Charge (fC)");
  charge_pass->GetYaxis()->SetTitleOffset(0.6);
  charge_pass->SetLabelSize(.04,"XY");
  charge_pass->SetMinimum(1.0);
  charge_pass->SetTitle("All Charge");
  charge_pass->Draw();
  charge_fail->SetLineColor(2);
  charge_fail->SetLineStyle(2);  
  charge_fail->Draw("same"); 
  TLegend *tl1 = new TLegend(0.6,0.65,0.9,0.9);
  tl1->SetHeader("50 GeV pion MC");
  tl1->AddEntry(charge_fail,"Delayed Charge");
  tl1->AddEntry(charge_pass,"Prompt Charge");
  tl1->Draw();
  c1_1->Modified();
  c1->cd();
  
  TPad *c1_2 = new TPad("c1_2", "c1_2",0.0,0.0,1.0,0.35);
  c1_2->Draw();
  c1_2->cd();
  c1_2->SetTopMargin(0.01926445);
  c1_2->SetBottomMargin(0.2802102);
 
  TH1F* charge_ratio = (TH1F*) charge_fail->Clone("charge_ratio");
  charge_ratio->SetTitleSize(0.12,"XY");
  charge_ratio->Divide(charge_pass);
  charge_ratio->GetYaxis()->SetTitle("Delayed / Prompt");
  charge_ratio->GetXaxis()->SetTitle("Charge (fC)");
  charge_ratio->SetLabelSize(.08,"XY");
  //charge_ratio->SetLabelSize(32,"X");
  charge_ratio->GetYaxis()->SetTitleOffset(0.3);
  charge_ratio->SetTitle("");
  charge_ratio->SetLineColor(4);
  charge_ratio->SetLineStyle(1);
  charge_ratio->Draw();

  c1->Print("plot_charge_prompt_vs_delayed_pi50.png");
  c1->Print("plot_charge_prompt_vs_delayed_pi50.pdf");

  //----------

  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  TPad *c2_1 = new TPad("c2_1", "c2_1",0.0,0.32,1.0,1.0);
  c2_1->Draw();
  c2_1->cd();

  gPad->SetLogy();
  gPad->SetGrid(1);

  charge_pass_data->SetTitleSize(0.07,"XY");
  charge_pass_data->GetYaxis()->SetTitle("Number of hits");
  charge_pass_data->GetXaxis()->SetTitle("Charge (fC)");
  charge_pass_data->GetYaxis()->SetTitleOffset(0.6);
  charge_pass_data->SetLabelSize(.04,"XY");
  charge_pass_data->SetMinimum(1.0);
  charge_pass_data->SetTitle("All Charge");
  charge_pass_data->Draw();
  charge_fail_data->SetLineColor(2);
  charge_fail_data->SetLineStyle(2);  
  charge_fail_data->Draw("same"); 
  TLegend *tl2 = new TLegend(0.55,0.65,0.9,0.9);
  tl2->SetTextSize(0.05);
  tl2->SetHeader("2017A HEP17 Isobunch Data");
  tl2->AddEntry(charge_fail_data,"Delayed Charge");
  tl2->AddEntry(charge_pass_data,"Prompt Charge");
  tl2->Draw();
  c2_1->Modified();
  c2->cd();
  
  TPad *c2_2 = new TPad("c2_2", "c2_2",0.0,0.0,1.0,0.35);
  c2_2->Draw();
  c2_2->cd();
  c2_2->SetTopMargin(0.01926445);
  c2_2->SetBottomMargin(0.2802102);
 
  TH1F* charge_ratio_data = (TH1F*) charge_fail_data->Clone("charge_ratio_data");
  charge_ratio_data->SetTitleSize(0.12,"XY");
  charge_ratio_data->Divide(charge_pass_data);
  charge_ratio_data->GetYaxis()->SetTitle("Delayed / Prompt");
  charge_ratio_data->GetXaxis()->SetTitle("Charge (fC)");
  charge_ratio_data->SetLabelSize(.08,"XY");
  //charge_ratio_data->SetLabelSize(32,"X");
  charge_ratio_data->GetYaxis()->SetTitleOffset(0.3);
  charge_ratio_data->SetTitle("");
  charge_ratio_data->SetLineColor(4);
  charge_ratio_data->SetLineStyle(1);
  charge_ratio_data->Draw();

  c2->Print("plot_charge_prompt_vs_delayed_HLTPhysics2017A.png");
  c2->Print("plot_charge_prompt_vs_delayed_HLTPhysics2017A.pdf");

  //----------
  
  TCanvas *c3 = new TCanvas("c3","c3",800,800);
  gPad->SetLogy();
  gPad->SetGrid(1);

  TPad *c3_1 = new TPad("c3_1", "c3_1",0.0,0.5,0.5,1.0);
  c3_1->Draw();
  c3_1->cd();
  c3_1->SetLeftMargin(0.2);
  c3_1->SetBottomMargin(0.2);
  c3_1->SetRightMargin(0.05);

  SOIfrac_HB->SetTitleSize(0.06,"XY");
  SOIfrac_HB->GetYaxis()->SetTitle("Number of hits");
  SOIfrac_HB->GetXaxis()->SetTitle("Q(SOI)/Q(SOI:LastTS)");
  SOIfrac_HB->GetYaxis()->SetTitleOffset(1.2);
  SOIfrac_HB->SetLabelSize(.04,"XY");
  SOIfrac_HB->SetMinimum(1.0);
  SOIfrac_HB->SetTitle("All Charge");
  SOIfrac_HB->SetTitle("HB");
  SOIfrac_HB->Draw();

  c3_1->Modified();
  c3->cd();
  
  TPad *c3_2 = new TPad("c3_2", "c3_2",0.5,0.5,1.0,1.0);
  c3_2->Draw();
  c3_2->cd();
  c3_2->SetLeftMargin(0.2);
  c3_2->SetBottomMargin(0.2);
  c3_2->SetRightMargin(0.05);
  SOIfrac_HE->SetTitleSize(0.06,"XY");
  SOIfrac_HE->GetYaxis()->SetTitle("Number of hits");
  SOIfrac_HE->GetXaxis()->SetTitle("Q(SOI)/Q(SOI:LastTS)");
  SOIfrac_HE->GetYaxis()->SetTitleOffset(1.2);
  SOIfrac_HE->SetLabelSize(.04,"XY");
  SOIfrac_HE->SetMinimum(1.0);
  SOIfrac_HE->SetTitle("All Charge");
  SOIfrac_HE->SetTitle("HE");
  SOIfrac_HE->Draw();

  c3_2->Modified();
  c3->cd();

  TPad *c3_3 = new TPad("c3_3", "c3_3",0.0,0.0,0.5,0.5);
  c3_3->Draw();
  c3_3->cd();
  c3_3->SetLeftMargin(0.2);
  c3_3->SetBottomMargin(0.2);
  c3_3->SetRightMargin(0.05);

  PostSOIfrac_HB->SetTitleSize(0.06,"XY");
  PostSOIfrac_HB->GetYaxis()->SetTitle("Number of hits");
  PostSOIfrac_HB->GetXaxis()->SetTitle("Q(SOI+1)/Q(SOI:LastTS)");
  PostSOIfrac_HB->GetYaxis()->SetTitleOffset(1.2);
  PostSOIfrac_HB->SetLabelSize(.04,"XY");
  PostSOIfrac_HB->SetMinimum(1.0);
  PostSOIfrac_HB->SetTitle("All Charge");
  PostSOIfrac_HB->SetTitle("HB");
  PostSOIfrac_HB->Draw();

  c3_3->Modified();
  c3->cd();
  
  TPad *c3_4 = new TPad("c3_4", "c3_4",0.5,0.0,1.0,0.5);
  c3_4->Draw();
  c3_4->cd();
  c3_4->SetLeftMargin(0.2);
  c3_4->SetBottomMargin(0.2);
  c3_4->SetRightMargin(0.05);
  PostSOIfrac_HE->SetTitleSize(0.06,"XY");
  PostSOIfrac_HE->GetYaxis()->SetTitle("Number of hits");
  PostSOIfrac_HE->GetXaxis()->SetTitle("Q(SOI+1)/Q(SOI:LastTS)");
  PostSOIfrac_HE->GetYaxis()->SetTitleOffset(1.2);
  PostSOIfrac_HE->SetLabelSize(.04,"XY");
  PostSOIfrac_HE->SetMinimum(1.0);
  PostSOIfrac_HE->SetTitle("All Charge");
  PostSOIfrac_HE->SetTitle("HE");
  PostSOIfrac_HE->Draw();

  c3_4->Modified();
  c3->cd();

  c3->Print("plot_SOIfrac_pi50.png");
  c3->Print("plot_SOIfrac_pi50.pdf");
 
  //----------
  
  TCanvas *c4 = new TCanvas("c4","c4",800,800);
  gPad->SetLogy();
  gPad->SetGrid(1);

  TPad *c4_1 = new TPad("c4_1", "c4_1",0.0,0.5,0.5,1.0);
  c4_1->Draw();
  c4_1->cd();
  c4_1->SetLeftMargin(0.2);
  c4_1->SetBottomMargin(0.2);
  c4_1->SetRightMargin(0.05);

  SOIfrac_HB_data->SetTitleSize(0.06,"XY");
  SOIfrac_HB_data->GetYaxis()->SetTitle("Number of hits");
  SOIfrac_HB_data->GetXaxis()->SetTitle("Q(SOI)/Q(SOI:LastTS)");
  SOIfrac_HB_data->GetYaxis()->SetTitleOffset(1.2);
  SOIfrac_HB_data->SetLabelSize(.04,"XY");
  SOIfrac_HB_data->SetMinimum(1.0);
  SOIfrac_HB_data->SetTitle("All Charge");
  SOIfrac_HB_data->SetTitle("HB");
  SOIfrac_HB_data->Draw();

  c4_1->Modified();
  c4->cd();
  
  TPad *c4_2 = new TPad("c4_2", "c4_2",0.5,0.5,1.0,1.0);
  c4_2->Draw();
  c4_2->cd();
  c4_2->SetLeftMargin(0.2);
  c4_2->SetBottomMargin(0.2);
  c4_2->SetRightMargin(0.05);
  SOIfrac_HE_data->SetTitleSize(0.06,"XY");
  SOIfrac_HE_data->GetYaxis()->SetTitle("Number of hits");
  SOIfrac_HE_data->GetXaxis()->SetTitle("Q(SOI)/Q(SOI:LastTS)");
  SOIfrac_HE_data->GetYaxis()->SetTitleOffset(1.2);
  SOIfrac_HE_data->SetLabelSize(.04,"XY");
  SOIfrac_HE_data->SetMinimum(1.0);
  SOIfrac_HE_data->SetTitle("All Charge");
  SOIfrac_HE_data->SetTitle("HEP17");
  SOIfrac_HE_data->Draw();

  c4_2->Modified();
  c4->cd();

  TPad *c4_3 = new TPad("c4_3", "c4_3",0.0,0.0,0.5,0.5);
  c4_3->Draw();
  c4_3->cd();
  c4_3->SetLeftMargin(0.2);
  c4_3->SetBottomMargin(0.2);
  c4_3->SetRightMargin(0.05);

  PostSOIfrac_HB_data->SetTitleSize(0.06,"XY");
  PostSOIfrac_HB_data->GetYaxis()->SetTitle("Number of hits");
  PostSOIfrac_HB_data->GetXaxis()->SetTitle("Q(SOI+1)/Q(SOI:LastTS)");
  PostSOIfrac_HB_data->GetYaxis()->SetTitleOffset(1.2);
  PostSOIfrac_HB_data->SetLabelSize(.04,"XY");
  PostSOIfrac_HB_data->SetMinimum(1.0);
  PostSOIfrac_HB_data->SetTitle("All Charge");
  PostSOIfrac_HB_data->SetTitle("HB");
  PostSOIfrac_HB_data->Draw();

  c4_3->Modified();
  c4->cd();
  
  TPad *c4_4 = new TPad("c4_4", "c4_4",0.5,0.0,1.0,0.5);
  c4_4->Draw();
  c4_4->cd();
  c4_4->SetLeftMargin(0.2);
  c4_4->SetBottomMargin(0.2);
  c4_4->SetRightMargin(0.05);
  PostSOIfrac_HE_data->SetTitleSize(0.06,"XY");
  PostSOIfrac_HE_data->GetYaxis()->SetTitle("Number of hits");
  PostSOIfrac_HE_data->GetXaxis()->SetTitle("Q(SOI+1)/Q(SOI:LastTS)");
  PostSOIfrac_HE_data->GetYaxis()->SetTitleOffset(1.2);
  PostSOIfrac_HE_data->SetLabelSize(.04,"XY");
  PostSOIfrac_HE_data->SetMinimum(1.0);
  PostSOIfrac_HE_data->SetTitle("All Charge");
  PostSOIfrac_HE_data->SetTitle("HEP17");
  PostSOIfrac_HE_data->Draw();

  c4_4->Modified();
  c4->cd();

  c4->Print("plot_SOIfrac_HLTPhysics2017A.png");
  c4->Print("plot_SOIfrac_HLTPhysics2017A.pdf");
 
  //----------
  
  TCanvas *c5 = new TCanvas("c5","c5",800,600);
  gPad->SetLogy();
  gPad->SetGrid(1);

  SOIfrac_HE->SetTitleSize(0.05,"XY");
  SOIfrac_HE->GetYaxis()->SetTitle("# of hits");
  SOIfrac_HE->GetXaxis()->SetTitle("Q(SOI)/Q(SOI:LastTS)");
  SOIfrac_HE->GetYaxis()->SetTitleOffset(0.8);
  SOIfrac_HE->SetLabelSize(.03,"XY");
  SOIfrac_HE->SetMinimum(1.0);
  SOIfrac_HE->SetTitle("HE");
  SOIfrac_HE->Draw();
  SOIfrac_1000_HE->SetLineColor(2);
  SOIfrac_1000_HE->SetLineStyle(2);
  SOIfrac_1000_HE->Draw("same");
  //SOIfrac_2000_HE->SetLineColor(4);
  //SOIfrac_2000_HE->SetLineStyle(2);
  //SOIfrac_2000_HE->Draw("same");
  SOIfrac_6000_HE->SetLineColor(6);
  SOIfrac_6000_HE->SetLineStyle(4);
  SOIfrac_6000_HE->Draw("same");
  TLegend *tl5 = new TLegend(0.7,0.7,0.9,0.9);
  tl5->AddEntry(SOIfrac_HE,"Q>300 (fC)");
  tl5->AddEntry(SOIfrac_1000_HE,"Q>1000 (fC)");
  tl5->AddEntry(SOIfrac_6000_HE,"Q>6000 (fC)");
  tl5->Draw();

  c5->Print("plot_SOIfrac_pi50_Qdep.png");
  c5->Print("plot_SOIfrac_pi50_Qdep.pdf");
  
  //----------
  
  TCanvas *c5b = new TCanvas("c5b","c5b",800,600);
  //gPad->SetLogy();
  gPad->SetGrid(1);

  SOIfrac_HE_normalized->SetTitleSize(0.05,"XY");
  SOIfrac_HE_normalized->GetYaxis()->SetTitle("a.u.");
  SOIfrac_HE_normalized->GetXaxis()->SetTitle("Q(SOI)/Q(SOI:LastTS)");
  SOIfrac_HE_normalized->GetYaxis()->SetTitleOffset(0.8);
  SOIfrac_HE_normalized->SetLabelSize(.03,"XY");
  SOIfrac_HE_normalized->SetMinimum(1.0);
  SOIfrac_HE_normalized->SetTitle("HE");
  SOIfrac_HE_normalized->Draw("hist");
  SOIfrac_1000_HE_normalized->SetLineColor(2);
  SOIfrac_1000_HE_normalized->SetLineStyle(2);
  SOIfrac_1000_HE_normalized->Draw("same,hist");
  //SOIfrac_2000_HE_normalized->SetLineColor(4);
  //SOIfrac_2000_HE_normalized->SetLineStyle(2);
  //SOIfrac_2000_HE_normalized->Draw("same");
  SOIfrac_6000_HE_normalized->SetLineColor(6);
  SOIfrac_6000_HE_normalized->SetLineStyle(4);
  SOIfrac_6000_HE_normalized->Draw("same,hist");
  //TLegend *tl5 = new TLegend(0.7,0.7,0.9,0.9);
  //tl5->AddEntry(SOIfrac_HE,"Q>300 (fC)");
  //tl5->AddEntry(SOIfrac_1000_HE,"Q>1000 (fC)");
  //tl5->AddEntry(SOIfrac_6000_HE,"Q>6000 (fC)");
  tl5->Draw();

  c5b->Print("plot_SOIfrac_pi50_Qdep_normalized.png");
  c5b->Print("plot_SOIfrac_pi50_Qdep_normalized.pdf");
  
  //----------
  
  TCanvas *c6 = new TCanvas("c6","c6",800,600);
  //gPad->SetLogy();
  gPad->SetGrid(1);

  SOIfrac_pass_HE->SetTitleSize(0.05,"XY");
  SOIfrac_pass_HE->GetYaxis()->SetTitle("# of hits");
  SOIfrac_pass_HE->GetXaxis()->SetTitle("Q(SOI)/Q(SOI:LastTS)");
  SOIfrac_pass_HE->GetYaxis()->SetTitleOffset(0.8);
  SOIfrac_pass_HE->SetLabelSize(.03,"XY");
  SOIfrac_pass_HE->SetMinimum(1.0);
  SOIfrac_pass_HE->SetTitle("HE");
  SOIfrac_pass_HE->Draw();
  SOIfrac_fail_HE->SetLineColor(2);
  SOIfrac_fail_HE->SetLineStyle(2);
  SOIfrac_fail_HE->Draw("same");
  TLegend *tl6 = new TLegend(0.6,0.7,0.9,0.9);
  tl6->AddEntry(SOIfrac_pass_HE,"Q(SOI)>Q(any non-SOI)");
  tl6->AddEntry(SOIfrac_fail_HE,"Q(any non-SOI)>Q(SOI)");
  tl6->Draw();

  c6->Print("plot_SOIfrac_pi50_passfail.png");
  c6->Print("plot_SOIfrac_pi50_passfail.pdf");
  
  //----------
  
  TCanvas *c7 = new TCanvas("c7","c7",800,600);
  gPad->SetLogy();
  gPad->SetGrid(1);

  simhit_avetime_prompt->SetTitleSize(0.05,"XY");
  simhit_avetime_prompt->GetYaxis()->SetTitle("# of digis");
  simhit_avetime_prompt->GetXaxis()->SetTitle("Energy-weighted simhit time (ns)");
  simhit_avetime_prompt->GetYaxis()->SetTitleOffset(0.8);
  simhit_avetime_prompt->SetLabelSize(.03,"XY");
  simhit_avetime_prompt->SetMinimum(1.0);
  simhit_avetime_prompt->SetTitle("HE");
  simhit_avetime_prompt->GetXaxis()->SetRange(0.,135.); 
  simhit_avetime_prompt->Draw();
  simhit_avetime_delayed->SetLineColor(2);
  simhit_avetime_delayed->SetLineStyle(2);
  simhit_avetime_delayed->Draw("same");
  TLegend *tl7 = new TLegend(0.6,0.7,0.9,0.9);
  tl7->AddEntry(simhit_avetime_prompt,"Q(SOI)>Q(any non-SOI)");
  tl7->AddEntry(simhit_avetime_delayed,"Q(any non-SOI)>Q(SOI)");
  tl7->Draw();

  c7->Print("plot_simhit_avetime_pi50.png");
  c7->Print("plot_simhit_avetime_pi50.pdf");
  
  // //----------

  // TCanvas *c5 = new TCanvas("c5","c5",800,600);
  // gPad->SetLogy();
  // gPad->SetGrid(1);

  // hits_promptE->SetTitleSize(0.04,"XY");
  // hits_promptE->GetYaxis()->SetTitle("Energy of hits");
  // hits_promptE->GetXaxis()->SetTitle("Timing (ns)");
  // hits_promptE->GetYaxis()->SetTitleOffset(0.9);
  // chargef->SetLabelSize(.03,"XY");
  // //hits_promptE->SetMinimum(100.0);
  // hits_promptE->SetTitle("All Hits");
  // hits_promptE->Draw();
  // hits_delayedE->SetLineColor(2);
  // hits_delayedE->SetLineStyle(2);
  // hits_delayedE->Draw("same");
  // TLegend *tl3 = new TLegend(0.7,0.7,0.9,0.9);
  // tl3->AddEntry(charge_fail,"Delayed Hits");
  // tl3->AddEntry(charge_pass,"Prompt Hits");
  // tl3->Draw();

  // //c2->Print("plot2.png");
  // //c3->Print("plot3.png");
  // c4->Print("plot4.png");
  // c5->Print("plot5.png"); 

}

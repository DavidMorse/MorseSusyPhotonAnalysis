void temp(){

  gStyle->SetOptStat(0);

  TFile f1("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Photon_handMade_DiPho_Born_Pt20_doubleEMEnriched_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July10.root","READ");
  TFile f2("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Photon_handMade_DiPho_Box_Pt20_doubleEMEnriched_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July10_2.root","READ");
  TFile f3("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Photon_handMade_QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July12.root","READ");
  TFile f4("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Photon_handMade_QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July12.root","READ");
  TFile f5("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Photon_handMade_GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July10.root","READ");
  TFile f6("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Photon_handMade_GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July10.root","READ");
  
  TH1F* f1_0Jet=(TH1F*)f1.Get("ggffDiJetPtRatio_0Jet");TH1F* f1_1Jet=(TH1F*)f1.Get("ggffDiJetPtRatio_1Jet");TH1F* f1_2Jet=(TH1F*)f1.Get("ggffDiJetPtRatio_2Jet");
  TH1F* f2_0Jet=(TH1F*)f2.Get("ggffDiJetPtRatio_0Jet");TH1F* f2_1Jet=(TH1F*)f2.Get("ggffDiJetPtRatio_1Jet");TH1F* f2_2Jet=(TH1F*)f2.Get("ggffDiJetPtRatio_2Jet");
  TH1F* f3_0Jet=(TH1F*)f3.Get("ggffDiJetPtRatio_0Jet");TH1F* f3_1Jet=(TH1F*)f3.Get("ggffDiJetPtRatio_1Jet");TH1F* f3_2Jet=(TH1F*)f3.Get("ggffDiJetPtRatio_2Jet");
  TH1F* f4_0Jet=(TH1F*)f4.Get("ggffDiJetPtRatio_0Jet");TH1F* f4_1Jet=(TH1F*)f4.Get("ggffDiJetPtRatio_1Jet");TH1F* f4_2Jet=(TH1F*)f4.Get("ggffDiJetPtRatio_2Jet");
  TH1F* f5_0Jet=(TH1F*)f5.Get("ggffDiJetPtRatio_0Jet");TH1F* f5_1Jet=(TH1F*)f5.Get("ggffDiJetPtRatio_1Jet");TH1F* f5_2Jet=(TH1F*)f5.Get("ggffDiJetPtRatio_2Jet");
  TH1F* f6_0Jet=(TH1F*)f6.Get("ggffDiJetPtRatio_0Jet");TH1F* f6_1Jet=(TH1F*)f6.Get("ggffDiJetPtRatio_1Jet");TH1F* f6_2Jet=(TH1F*)f6.Get("ggffDiJetPtRatio_2Jet");

  f1_0Jet->SetLineColor(kBlue);f1_0Jet->SetMarkerColor(kBlue);f1_1Jet->SetLineColor(kBlue);f1_1Jet->SetMarkerColor(kBlue);f1_2Jet->SetLineColor(kBlue);f1_2Jet->SetMarkerColor(kBlue);
  f2_0Jet->SetLineColor(kRed);f2_0Jet->SetMarkerColor(kRed);f2_1Jet->SetLineColor(kRed);f2_1Jet->SetMarkerColor(kRed);f2_2Jet->SetLineColor(kRed);f2_2Jet->SetMarkerColor(kRed);
  f3_0Jet->SetLineColor(kCyan);f3_0Jet->SetMarkerColor(kCyan);f3_1Jet->SetLineColor(kCyan);f3_1Jet->SetMarkerColor(kCyan);f3_2Jet->SetLineColor(kCyan);f3_2Jet->SetMarkerColor(kCyan);
  f4_0Jet->SetLineColor(kViolet);f4_0Jet->SetMarkerColor(kViolet);f4_1Jet->SetLineColor(kViolet);f4_1Jet->SetMarkerColor(kViolet);f4_2Jet->SetLineColor(kViolet);f4_2Jet->SetMarkerColor(kViolet);
  f5_0Jet->SetLineColor(kOrange);f5_0Jet->SetMarkerColor(kOrange);f5_1Jet->SetLineColor(kOrange);f5_1Jet->SetMarkerColor(kOrange);f5_2Jet->SetLineColor(kOrange);f5_2Jet->SetMarkerColor(kOrange);
  f6_0Jet->SetLineColor(kBlack);f6_0Jet->SetMarkerColor(kBlack);f6_1Jet->SetLineColor(kBlack);f6_1Jet->SetMarkerColor(kBlack);f6_2Jet->SetLineColor(kBlack);f6_2Jet->SetMarkerColor(kBlack);

  TCanvas *c1 = new TCanvas("c1","c1",800,600);c1->cd();

  f1_0Jet->GetXaxis()->SetRangeUser(0,149);f1_1Jet->GetXaxis()->SetRangeUser(0,199);f1_2Jet->GetXaxis()->SetRangeUser(0,199);

  f1_0Jet->Draw("PE");f2_0Jet->Draw("PESAMES");f3_0Jet->Draw("PESAMES");f4_0Jet->Draw("PESAMES");f5_0Jet->Draw("PESAMES");f6_0Jet->Draw("PESAMES");
  c1->Print("Plots/Closure/handmade/moreEvents/AllMCdijetpt_0jet.png");

  f1_1Jet->Draw("PE");f2_1Jet->Draw("PESAMES");f3_1Jet->Draw("PESAMES");f4_1Jet->Draw("PESAMES");f5_1Jet->Draw("PESAMES");f6_1Jet->Draw("PESAMES");
  c1->Print("Plots/Closure/handmade/moreEvents/AllMCdijetpt_1jet.png");
  
  f1_2Jet->Draw("PE");f2_2Jet->Draw("PESAMES");f3_2Jet->Draw("PESAMES");f4_2Jet->Draw("PESAMES");f5_2Jet->Draw("PESAMES");f6_2Jet->Draw("PESAMES");
  c1->Print("Plots/Closure/handmade/moreEvents/AllMCdijetpt_2jet.png");

}

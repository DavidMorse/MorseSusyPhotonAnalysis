#include "TPaveStats.h"

void Jan22Compare(){

  gStyle->SetOptStat(0);

  TFile f_Jan22("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Jan22ReReco_2012A_Filter.root");
  TFile f_Old("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Data2012Aonly_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_15GevCHfakeCut_PixelVetoOnFakes_NewDiEMPtBins_19499pb_NewJetMatching_2JetReq_PUjetID_NewDoGG_May28.root");

  TCanvas c1("c1","c1",800,600);
  c1.cd();

  //c1.SetLogy(1);
  TH1F* hist = (TH1F*)f_Old.Get("ggMet");hist->Sumw2();
  float x = 1./hist->Integral();hist->Scale(x);
  TH1F* histJan22 = (TH1F*)f_Jan22.Get("ggMet");histJan22->Sumw2();
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,89);hist->SetTitle("ggMet");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ggMet.png");

  hist = (TH1F*)f_Old.Get("eeMet");hist->Sumw2();
  float x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("eeMet");histJan22->Sumw2();
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,99);hist->SetTitle("eeMet");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/eeMet.png");

  hist = (TH1F*)f_Old.Get("ffMet");hist->Sumw2();
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ffMet");histJan22->Sumw2();
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,74);hist->SetTitle("ffMet");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ffMet.png");

  hist = (TH1F*)f_Old.Get("met");
  float x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("met");
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,89);hist->SetTitle("met");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/met.png");

  hist = (TH1F*)f_Old.Get("Type1CorrMet");hist->Sumw2();
  float x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("Type1CorrMet");histJan22->Sumw2();
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,89);hist->SetTitle("Type1CorrMet");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/Type1CorrMet.png");

  hist = (TH1F*)f_Old.Get("ggPtLead");hist->Sumw2();
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ggPtLead");histJan22->Sumw2();
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.5);hist->GetXaxis()->SetRangeUser(35,299);hist->SetTitle("ggPtLead");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ggPtLead.png");

  hist = (TH1F*)f_Old.Get("ggPtTrail");hist->Sumw2();
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ggPtTrail");histJan22->Sumw2();
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(20,199);hist->SetTitle("ggPtTrail");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ggPtTrail.png");
  c1.SetLogy(0);
  hist = (TH1F*)f_Old.Get("ggDiEMPt");
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ggDiEMPt");
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,249);hist->SetTitle("ggDiEMPt");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ggDiEMPt.png");

  hist = (TH1F*)f_Old.Get("ggDiJetPt");
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ggDiJetPt");
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,249);hist->SetTitle("ggDiJetPt");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ggDiJetPt.png");
  
  c1.SetLogy(0);
  hist = (TH1F*)f_Old.Get("ggDiJetPt_0Jet");
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ggDiJetPt_0Jet");
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,119);hist->SetTitle("ggDiJetPt_0Jet");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ggDiJetPt_0Jet.png");

  hist = (TH1F*)f_Old.Get("ggDiJetPt_1Jet");hist->Rebin(2);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ggDiJetPt_1Jet");histJan22->Rebin(2);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,199);hist->SetTitle("ggDiJetPt_1Jet");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ggDiJetPt_1Jet.png");
  
  hist = (TH1F*)f_Old.Get("ggDiJetPt_2Jet");hist->Rebin(4);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ggDiJetPt_2Jet");histJan22->Rebin(4);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,299);hist->SetTitle("ggDiJetPt_2Jet");//hist->GetYaxis()->SetRangeUser(0,.08);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ggDiJetPt_2Jet.png");

  //c1.SetLogy(1);
  hist = (TH1F*)f_Old.Get("ffDiEMPt");
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ffDiEMPt");
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,149);hist->SetTitle("ffDiEMPt");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ffDiEMPt.png");

  hist = (TH1F*)f_Old.Get("ffDiJetPt");hist->Rebin(2);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ffDiJetPt");histJan22->Rebin(2);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,199);hist->SetTitle("ffDiJetPt");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ffDiJetPt.png");
  
  c1.SetLogy(0);
  hist = (TH1F*)f_Old.Get("ffDiJetPt_0Jet");
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ffDiJetPt_0Jet");
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,139);hist->SetTitle("ffDiJetPt_0Jet");hist->GetYaxis()->SetRangeUser(0,.095);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ffDiJetPt_0Jet.png");

  hist = (TH1F*)f_Old.Get("ffDiJetPt_1Jet");hist->Rebin(2);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ffDiJetPt_1Jet");histJan22->Rebin(2);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,199);hist->SetTitle("ffDiJetPt_1Jet");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ffDiJetPt_1Jet.png");
  
  hist = (TH1F*)f_Old.Get("ffDiJetPt_2Jet");hist->Rebin(5);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ffDiJetPt_2Jet");histJan22->Rebin(5);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,229);hist->SetTitle("ffDiJetPt_2Jet");//hist->GetYaxis()->SetRangeUser(0,.08);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ffDiJetPt_2Jet.png");

  //c1.SetLogy(1);
  hist = (TH1F*)f_Old.Get("eeDiEMPt");
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("eeDiEMPt");
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,249);hist->SetTitle("eeDiEMPt");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/eeDiEMPt.png");
  
  hist = (TH1F*)f_Old.Get("eeDiJetPt");
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("eeDiJetPt");
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,249);hist->SetTitle("eeDiJetPt");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/eeDiJetPt.png");
  
  c1.SetLogy(0);
  hist = (TH1F*)f_Old.Get("eeDiJetPt_0Jet");
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("eeDiJetPt_0Jet");
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,120);hist->SetTitle("eeDiJetPt_0Jet");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/eeDiJetPt_0Jet.png");

  hist = (TH1F*)f_Old.Get("eeDiJetPt_1Jet");hist->Rebin(2);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("eeDiJetPt_1Jet");histJan22->Rebin(2);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,349);hist->SetTitle("eeDiJetPt_1Jet");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/eeDiJetPt_1Jet.png");
  
  hist = (TH1F*)f_Old.Get("eeDiJetPt_2Jet");hist->Rebin(4);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("eeDiJetPt_2Jet");histJan22->Rebin(4);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,399);hist->SetTitle("eeDiJetPt_2Jet");//hist->GetYaxis()->SetRangeUser(0,.08);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/eeDiJetPt_2Jet.png");
  
  hist = (TH1F*)f_Old.Get("PfJet_phoCleaned_pT");hist->Sumw2();//hist->Rebin(4);
  cout<<"PhoCleaned Integral Old Files: "<<hist->Integral()<<endl;
  TH1F* histRat=(TH1F*)hist->Clone();histRat->Rebin(2);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("PfJet_phoCleaned_pT");histJan22->Sumw2();//histJan22->Rebin(4);
  cout<<"PhoCleaned Integral Jan22    : "<<histJan22->Integral()<<endl;
  TH1F* histRat22=(TH1F*)histJan22->Clone();histRat22->Rebin(2);histRat22->Divide(histRat);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(25,199);hist->SetTitle("PfJet_phoCleaned_pT");hist->GetYaxis()->SetRangeUser(0,.2);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/PfJet_phoCleaned_pT.png");
  histRat22->GetXaxis()->SetRangeUser(30,299);histRat22->SetMarkerSize(0.75);histRat22->Draw("PE");TLine ratLine(30,1,300,1);ratLine.Draw();
  c1.Print("Plots/Jan22Compare/PfJet_phoCleaned_pT_ratio.png");


  hist = (TH1F*)f_Old.Get("PfJet_phoCleaned_eta");hist->Sumw2();hist->Rebin(5);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("PfJet_phoCleaned_eta");histJan22->Sumw2();histJan22->Rebin(5);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(-2.7,2.7);hist->SetTitle("PfJet_phoCleaned_eta");//hist->GetYaxis()->SetRangeUser(0,.2);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/PfJet_phoCleaned_eta.png");
  
  hist = (TH1F*)f_Old.Get("PfJet_phoCleaned_phi");hist->Sumw2();hist->Rebin(2);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("PfJet_phoCleaned_phi");histJan22->Sumw2();histJan22->Rebin(2);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);//hist->GetXaxis()->SetRangeUser(-1.5,1.5);hist->SetTitle("PfJet_phoCleaned_phi");//hist->GetYaxis()->SetRangeUser(0,.2);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/PfJet_phoCleaned_phi.png");
  
  hist = (TH1F*)f_Old.Get("ggEta");hist->Sumw2();hist->Rebin(5);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ggEta");histJan22->Sumw2();histJan22->Rebin(5);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(-1.5,1.5);hist->SetTitle("ggEta");//hist->GetYaxis()->SetRangeUser(0,.2);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ggEta.png");
  
  hist = (TH1F*)f_Old.Get("ggPhi");hist->Sumw2();hist->Rebin(2);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ggPhi");histJan22->Sumw2();histJan22->Rebin(2);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(-1.5,1.5);hist->SetTitle("ggPhi");//hist->GetYaxis()->SetRangeUser(0,.2);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ggPhi.png");
  
  hist = (TH1F*)f_Old.Get("ffEta");hist->Sumw2();hist->Rebin(5);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ffEta");histJan22->Sumw2();histJan22->Rebin(5);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(-1.5,1.5);hist->SetTitle("ffEta");//hist->GetYaxis()->SetRangeUser(0,.2);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ffEta.png");
  
  hist = (TH1F*)f_Old.Get("ffPhi");hist->Sumw2();hist->Rebin(2);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ffPhi");histJan22->Sumw2();histJan22->Rebin(2);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(-1.5,1.5);hist->SetTitle("ffPhi");//hist->GetYaxis()->SetRangeUser(0,.2);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/ffPhi.png");
  
  hist = (TH1F*)f_Old.Get("eeEta");hist->Sumw2();hist->Rebin(5);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("eeEta");histJan22->Sumw2();histJan22->Rebin(5);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(-1.5,1.5);hist->SetTitle("eeEta");//hist->GetYaxis()->SetRangeUser(0,.2);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/eeEta.png");
  
  hist = (TH1F*)f_Old.Get("eePhi");hist->Sumw2();hist->Rebin(2);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("eePhi");histJan22->Sumw2();histJan22->Rebin(2);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(-1.5,1.5);hist->SetTitle("eePhi");//hist->GetYaxis()->SetRangeUser(0,.2);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/eePhi.png");

  
  hist = (TH1F*)f_Old.Get("PfJet_phoMatched_gg_pT");hist->Sumw2();//hist->Rebin(4);
  TH1F* histRat=(TH1F*)hist->Clone();histRat->Rebin(2);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("PfJet_phoMatched_gg_pT");histJan22->Sumw2();//histJan22->Rebin(4);
  TH1F* histRat22=(TH1F*)histJan22->Clone();histRat22->Rebin(2);histRat22->Divide(histRat);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,199);hist->SetTitle("PfJet_phoMatched_gg_pT");hist->GetYaxis()->SetRangeUser(0,.158);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/PfJet_phoMatched_gg_pT.png");
  histRat22->GetXaxis()->SetRangeUser(20,249);histRat22->SetMarkerSize(0.75);histRat22->Draw("PE");ratLine.DrawLine(20,1,250,1);
  c1.Print("Plots/Jan22Compare/PfJet_phoMatched_pT_ratio.png");
  
  hist = (TH1F*)f_Old.Get("PfJet_phoMatched_gg_eta");hist->Sumw2();hist->Rebin(5);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("PfJet_phoMatched_gg_eta");histJan22->Sumw2();histJan22->Rebin(5);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(-1.5,1.5);hist->SetTitle("PfJet_phoMatched_gg_eta");//hist->GetYaxis()->SetRangeUser(0,.2);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/PfJet_phoMatched_gg_eta.png");
  
  hist = (TH1F*)f_Old.Get("PfJet_phoMatched_gg_phi");hist->Sumw2();hist->Rebin(2);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("PfJet_phoMatched_gg_phi");histJan22->Sumw2();histJan22->Rebin(2);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);//hist->GetXaxis()->SetRangeUser(-1.5,1.5);hist->SetTitle("PfJet_phoMatched_gg_phi");//hist->GetYaxis()->SetRangeUser(0,.2);
  histJan22->SetLineColor(kRed);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/PfJet_phoMatched_gg_phi.png");

  
  Double_t DiEMPtBins[]={0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,80.,100.,150.,200.};
  int numBins=sizeof(DiEMPtBins)/sizeof(Double_t)-1;

  hist = (TH1F*)f_Old.Get("ggffDiJetPtRatio_0Jet");
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ggffDiJetPtRatio_0Jet");
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.75);hist->GetXaxis()->SetRangeUser(0,199);hist->SetTitle("ggffDiJetPtRatio_0Jet");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);histJan22->SetMarkerColor(kRed);histJan22->SetMarkerSize(0.75);
  histrebin = (TH1F*)hist->Rebin(numBins,"histrebin",DiEMPtBins);histrebin->GetXaxis()->SetRangeUser(0,199);histrebin->Draw("PE");//hist->Draw("histoSAMES");
  histJan22rebin = (TH1F*)histJan22->Rebin(numBins,"histJan22rebin",DiEMPtBins);histJan22rebin->Draw("peSAMES");
  c1.Print("Plots/Jan22Compare/ggffDiJetPtRatio_0Jet.png");

  hist = (TH1F*)f_Old.Get("ggffDiJetPtRatio_1Jet");
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ggffDiJetPtRatio_1Jet");
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.75);hist->GetXaxis()->SetRangeUser(0,199);hist->SetTitle("ggffDiJetPtRatio_1Jet");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);histJan22->SetMarkerColor(kRed);histJan22->SetMarkerSize(0.75);
  histrebin = (TH1F*)hist->Rebin(numBins,"histrebin",DiEMPtBins);histrebin->GetXaxis()->SetRangeUser(0,199);histrebin->GetYaxis()->SetRangeUser(0.02,.17);histrebin->Draw("PE");//hist->Draw("histoSAMES");
  histJan22rebin = (TH1F*)histJan22->Rebin(numBins,"histJan22rebin",DiEMPtBins);histJan22rebin->Draw("peSAMES");
  c1.Print("Plots/Jan22Compare/ggffDiJetPtRatio_1Jet.png");

  hist = (TH1F*)f_Old.Get("ggffDiJetPtRatio_2Jet");
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ggffDiJetPtRatio_2Jet");
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.75);hist->GetXaxis()->SetRangeUser(0,199);hist->SetTitle("ggffDiJetPtRatio_2Jet");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);histJan22->SetMarkerColor(kRed);histJan22->SetMarkerSize(0.75);
  histrebin = (TH1F*)hist->Rebin(numBins,"histrebin",DiEMPtBins);histrebin->GetXaxis()->SetRangeUser(0,199);histrebin->Draw("PE");//hist->Draw("histoSAMES");
  histJan22rebin = (TH1F*)histJan22->Rebin(numBins,"histJan22rebin",DiEMPtBins);histJan22rebin->Draw("peSAMES");
  c1.Print("Plots/Jan22Compare/ggffDiJetPtRatio_2Jet.png");

  hist = (TH1F*)f_Old.Get("ggeeDiJetPtRatio_0Jet");
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ggeeDiJetPtRatio_0Jet");
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.75);hist->GetXaxis()->SetRangeUser(0,199);hist->SetTitle("ggeeDiJetPtRatio_0Jet");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);histJan22->SetMarkerColor(kRed);histJan22->SetMarkerSize(0.75);
  histrebin = (TH1F*)hist->Rebin(numBins,"histrebin",DiEMPtBins);histrebin->GetXaxis()->SetRangeUser(0,199);histrebin->Draw("PE");//hist->Draw("histoSAMES");
  histJan22rebin = (TH1F*)histJan22->Rebin(numBins,"histJan22rebin",DiEMPtBins);histJan22rebin->Draw("peSAMES");
  c1.Print("Plots/Jan22Compare/ggeeDiJetPtRatio_0Jet.png");

  hist = (TH1F*)f_Old.Get("ggeeDiJetPtRatio_1Jet");
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ggeeDiJetPtRatio_1Jet");
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.75);hist->GetXaxis()->SetRangeUser(0,199);hist->SetTitle("ggeeDiJetPtRatio_1Jet");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);histJan22->SetMarkerColor(kRed);histJan22->SetMarkerSize(0.75);
  histrebin = (TH1F*)hist->Rebin(numBins,"histrebin",DiEMPtBins);histrebin->GetXaxis()->SetRangeUser(0,199);histrebin->GetYaxis()->SetRangeUser(0.02,.17);histrebin->Draw("PE");//hist->Draw("histoSAMES");
  histJan22rebin = (TH1F*)histJan22->Rebin(numBins,"histJan22rebin",DiEMPtBins);histJan22rebin->Draw("peSAMES");
  c1.Print("Plots/Jan22Compare/ggeeDiJetPtRatio_1Jet.png");

  hist = (TH1F*)f_Old.Get("ggeeDiJetPtRatio_2Jet");
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ggeeDiJetPtRatio_2Jet");
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.75);hist->GetXaxis()->SetRangeUser(0,199);hist->SetTitle("ggeeDiJetPtRatio_2Jet");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);histJan22->SetMarkerColor(kRed);histJan22->SetMarkerSize(0.75);
  histrebin = (TH1F*)hist->Rebin(numBins,"histrebin",DiEMPtBins);histrebin->GetXaxis()->SetRangeUser(0,199);histrebin->Draw("PE");//hist->Draw("histoSAMES");
  histJan22rebin = (TH1F*)histJan22->Rebin(numBins,"histJan22rebin",DiEMPtBins);histJan22rebin->Draw("peSAMES");
  c1.Print("Plots/Jan22Compare/ggeeDiJetPtRatio_2Jet.png");

  hist = (TH1F*)f_Old.Get("ggInvarMass");hist->Sumw2();hist->Rebin(3);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("ggInvarMass");histJan22->Sumw2();histJan22->Rebin(3);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.75);hist->GetXaxis()->SetRangeUser(15,199);hist->SetTitle("ggInvarMass");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);histJan22->SetMarkerColor(kRed);histJan22->SetMarkerSize(0.75);
  hist->Draw("PE");
  histJan22->Draw("PESAMES");
  c1.Print("Plots/Jan22Compare/ggInvarMass.png");

  gStyle->SetOptStat(2200);
  hist = (TH1F*)f_Old.Get("eeInvarMass");hist->Sumw2();//hist->Rebin(3);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("eeInvarMass");histJan22->Sumw2();//histJan22->Rebin(3);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.75);hist->GetXaxis()->SetRangeUser(81,100.9);hist->SetTitle("eeInvarMass");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);histJan22->SetMarkerColor(kRed);histJan22->SetMarkerSize(0.75);
  hist->Draw("PE");
  histJan22->Draw("PESAMES");
  c1.Update();
  TPaveStats* sbff=(TPaveStats*)(hist->GetListOfFunctions()->FindObject("stats"));
  sbff->SetX1NDC(.62);sbff->SetX2NDC(.99);sbff->SetY1NDC(.75);sbff->SetY2NDC(.99);sbff->SetTextColor(kBlue);
  TPaveStats* sbgf=(TPaveStats*)(histJan22->GetListOfFunctions()->FindObject("stats"));
  sbgf->SetX1NDC(.62);sbgf->SetX2NDC(.99);sbgf->SetY1NDC(.5);sbgf->SetY2NDC(.749);sbgf->SetTextColor(kRed);
  c1.Print("Plots/Jan22Compare/eeInvarMass.png");

  hist = (TH1F*)f_Old.Get("eeInvarMassFullRange");hist->Sumw2();hist->Rebin(3);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("eeInvarMassFullRange");histJan22->Sumw2();histJan22->Rebin(3);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.75);hist->GetXaxis()->SetRangeUser(15,199);hist->SetTitle("eeInvarMassFullRange");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);histJan22->SetMarkerColor(kRed);histJan22->SetMarkerSize(0.75);
  hist->Draw("PE");
  histJan22->Draw("PESAMES");
  c1.Print("Plots/Jan22Compare/eeInvarMassFullRange.png");


  gStyle->SetOptStat(0);
  hist = (TH1F*)f_Jan22.Get("met");
  float x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("Type1CorrMet");
  TH1F* histJan222 = (TH1F*)f_Jan22.Get("Type01CorrMet");histJan222->Sumw2();
  x=1./histJan22->Integral();histJan22->Scale(x);
  x=1./histJan222->Integral();histJan222->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,89);hist->SetTitle("");//hist->GetYaxis()->SetRangeUser(3e-5,3e-1);
  histJan22->SetLineColor(kRed);
  histJan222->SetLineColor(kGreen);
  hist->Draw("PE");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  histJan222->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/AllMets.png");

  
  hist = (TH1F*)f_Jan22.Get("PfJet_phoMatched_gg_pTlead");hist->Sumw2();//hist->Rebin(4);
  x = 1./hist->Integral();hist->Scale(x);
  histJan22 = (TH1F*)f_Jan22.Get("PfJet_phoMatched_gg_pTtrail");histJan22->Sumw2();//histJan22->Rebin(4);
  x=1./histJan22->Integral();histJan22->Scale(x);
  hist->SetLineColor(kBlue);hist->SetMarkerColor(kBlue);hist->SetMarkerSize(0.8);hist->GetXaxis()->SetRangeUser(0,199);//hist->SetTitle("PfJet_phoMatched_gg_pT");hist->GetYaxis()->SetRangeUser(0,.158);
  histJan22->SetLineColor(kRed);
  hist->Draw("histo");//hist->Draw("histoSAMES");
  histJan22->Draw("histoSAMES");
  c1.Print("Plots/Jan22Compare/PfJet_phoMatched_gg_pTleadANDtrail.png");

}

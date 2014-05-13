#include "TH1F.h"
#include "TH1.h"

void compare(){

  gStyle->SetOptStat(0);


  TFile f_old("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_15GevCHfakeCut_PixelVetoOnFakes_NewDiEMPtBins_19499pb_NewJetMatching_2JetReq_May7.root","READ");
  TFile f_new("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_15GevCHfakeCut_PixelCutOnFakes_NewDiEMPtBins_19499pb_NewJetMatching_2JetReq_GFsignalHLTonlyApr23.root");

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->cd();

  TH1F *hist;

  TLegend leg(.25,.6,.6,.8);
  leg.SetFillColor(kWhite);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);

  TLegend leg2(.5,.6,.85,.8);
  leg2.SetFillColor(kWhite);
  leg2.SetFillStyle(0);
  leg2.SetBorderSize(0);

  hist = (TH1F*)f_new.Get("gggammafakeDiJetPtRatio_0Jet");
  hist->GetXaxis()->SetRangeUser(0,399);hist->GetYaxis()->SetRangeUser(0.5,2);
  hist->SetLineColor(kRed);hist->SetMarkerColor(kRed);
  hist->SetMarkerSize(.5);
  hist->Draw("PE");
  leg.AddEntry(hist,"signal HLT only");
  hist = (TH1F*)f_old.Get("gggammafakeDiJetPtRatio_0Jet"); 
  hist->SetMarkerSize(.5);
  hist->Draw("PEsames");  
  leg.AddEntry(hist,"any HLT path");
  leg.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_gf_0Jet.png");

  hist = (TH1F*)f_new.Get("gggammafakeDiJetPtRatio_1Jet");
  hist->GetXaxis()->SetRangeUser(0,599);hist->GetYaxis()->SetRangeUser(0.5,4.2);
  hist->SetLineColor(kRed);hist->SetMarkerColor(kRed);
  hist->SetMarkerSize(.5);
  hist->Draw("PE");
  hist = (TH1F*)f_old.Get("gggammafakeDiJetPtRatio_1Jet"); 
  hist->SetMarkerSize(.5);
  hist->Draw("PEsames");
  leg.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_gf_1Jet.png");

  hist = (TH1F*)f_new.Get("gggammafakeDiJetPtRatio_2Jet");
  hist->GetXaxis()->SetRangeUser(0,599);hist->GetYaxis()->SetRangeUser(0.6,3);
  hist->SetLineColor(kRed);hist->SetMarkerColor(kRed);
  hist->SetMarkerSize(.5);
  hist->Draw("PE");
  hist = (TH1F*)f_old.Get("gggammafakeDiJetPtRatio_2Jet"); 
  hist->SetMarkerSize(.5);
  hist->Draw("PEsames");
  leg.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_gf_2Jet.png");


  hist = (TH1F*)f_old.Get("gammafakeMet");TH1F* gammafakeMetNew = (TH1F*)hist->Clone();
  hist->SetTitle("#gammaf Raw E_{T}^{miss}");
  hist->GetXaxis()->SetRangeUser(0,50);hist->GetYaxis()->SetRangeUser(0,60000);
  float x = hist->Integral();
  leg2.AddEntry(hist,"any HLT path","l");
  hist->Draw("histo");
  hist = (TH1F*)f_new.Get("gammafakeMet_reweightJet_binned");
  hist->SetLineColor(kRed);hist->SetMarkerColor(kRed);hist->SetMarkerSize(0.5);
  x/=hist->Integral();hist->Scale(x);
  leg2.AddEntry(hist,"signal HLT only");
  hist->Draw("pesames");leg2.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_gf_met.png");


  hist = (TH1F*)f_old.Get("gammafakeMet_reweightJet_binned");
  hist->SetTitle("#gammaf diJetP_{T} reweighted E_{T}^{miss}");
  hist->GetXaxis()->SetRangeUser(0,50);hist->GetYaxis()->SetRangeUser(0,60000);
  float x = hist->Integral();
  hist->Draw("histo");
  hist = (TH1F*)f_new.Get("gammafakeMet_reweightJet_binned");
  hist->SetLineColor(kRed);hist->SetMarkerColor(kRed);hist->SetMarkerSize(0.5);
  x/=hist->Integral();hist->Scale(x);
  hist->Draw("pesames");leg2.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_gf_met_reweight.png");

  f_old.cd();
  TH1F* ff = (TH1F*)f_old.Get("ffMet");ff->Sumw2();TH1F* ffMetNew = (TH1F*)ff->Clone();
  TH1F* ff_ci_ci = (TH1F*)f_old.Get("ffMet_CaloIso_CaloIso");ff_ci_ci->Sumw2();
  TH1F* ff_ci_r9 = (TH1F*)f_old.Get("ffMet_CaloIso_R9Id");ff_ci_r9->Sumw2();
  TH1F* ff_r9_ci = (TH1F*)f_old.Get("ffMet_R9Id_CaloIso");ff_r9_ci->Sumw2();
  TH1F* ff_r9_r9 = (TH1F*)f_old.Get("ffMet_R9Id_R9Id");ff_r9_r9->Sumw2();
  float marksize=0.5;
  ff->SetTitle("");ff->GetXaxis()->SetRangeUser(0,49);ff->SetMarkerSize(marksize);
  ff->Draw();
  x=ff->Integral()/ff_ci_ci->Integral();ff_ci_ci->Scale(x);
  x=ff->Integral()/ff_ci_r9->Integral();ff_ci_r9->Scale(x);
  x=ff->Integral()/ff_r9_ci->Integral();ff_r9_ci->Scale(x);
  x=ff->Integral()/ff_r9_r9->Integral();ff_r9_r9->Scale(x);
  ff_ci_ci->SetLineColor(kBlue);ff_ci_ci->SetMarkerColor(kBlue);ff_ci_ci->SetMarkerSize(marksize);
  ff_ci_r9->SetLineColor(kGreen);ff_ci_r9->SetMarkerColor(kGreen);ff_ci_r9->SetMarkerSize(marksize);
  ff_r9_ci->SetLineColor(kOrange);ff_r9_ci->SetMarkerColor(kOrange);ff_r9_ci->SetMarkerSize(marksize);
  ff_r9_r9->SetLineColor(kRed);ff_r9_r9->SetMarkerColor(kRed);ff_r9_r9->SetMarkerSize(marksize);
  ff_ci_ci->Draw("sames");
  ff_ci_r9->Draw("sames");
  ff_r9_ci->Draw("sames");
  ff_r9_r9->Draw("sames");
  TLegend legMet(.29,.21,.59,.51);
  legMet.SetFillColor(kWhite);
  legMet.SetFillStyle(0);
  legMet.SetBorderSize(0);
  legMet.AddEntry(ff,"ff_all");legMet.AddEntry(ff_ci_ci,"ff_Caloiso_CaloIso");legMet.AddEntry(ff_ci_r9,"ff_Caloiso_R9Id");legMet.AddEntry(ff_r9_ci,"ff_R9Id_CaloIso");legMet.AddEntry(ff_r9_r9,"ff_R9Id_R9Id");
  legMet.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_met_raw.png");

  float top=79.;
  ff_ci_r9->GetYaxis()->SetRangeUser(0,3);
  ff_ci_r9->GetXaxis()->SetRangeUser(0,top);
  TH1F* ffcici = (TH1F*)ff_ci_ci->Clone();TH1F* ffcir9 = (TH1F*)ff_ci_r9->Clone();TH1F* ffr9ci = (TH1F*)ff_r9_ci->Clone();TH1F* ffr9r9 = (TH1F*)ff_r9_r9->Clone();
  //ff_ci_ci->Divide(ff);ff_ci_ci->Draw("PE");
  ff_ci_r9->Divide(ff_ci_ci);ff_ci_r9->Draw("PE");
  ff_r9_ci->Divide(ff_ci_ci);ff_r9_ci->Draw("PEsames");
  ff_r9_r9->Divide(ff_ci_ci);ff_r9_r9->Draw("PEsames");
  TLine line(0,1,top,1);line.Draw("sames");
  TLegend legRat(.23,.53,.58,.83);
  legRat.SetFillColor(kWhite);
  legRat.SetFillStyle(0);
  legRat.SetBorderSize(0);
  legRat.AddEntry(ff_ci_r9,"CaloIso_R9Id / CaloIso_CaloIso");
  legRat.AddEntry(ff_r9_ci,"R9Id_CaloIso / CaloIso_CaloIso");
  legRat.AddEntry(ff_r9_r9,"R9Id_R9Id / CaloIso_CaloIso");
  legRat.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_metRatio_raw.png");
  ffcir9->Add(ffr9ci);ffcir9->Add(ffr9r9);
  x = ff_ci_ci->Integral()/ffcir9->Integral();ffcir9->Scale(x);
  ffcir9->Divide(ff_ci_ci);ffcir9->SetLineColor(kViolet);ffcir9->SetMarkerColor(kViolet);ffcir9->GetYaxis()->SetRangeUser(0.8,2);
  ffcir9->Draw();line.Draw("SAMES");
  legRat.Clear();
  legRat.AddEntry(ffcir9,"R9 triggers / CaloIso_CaloIso trigger");
  legRat.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_metRatioComb_raw.png");

  //now exclusive
  //ff = (TH1F*)f_old.Get("ffMet");ff->Sumw2();
  ff_ci_ci = (TH1F*)f_old.Get("ffMet_CaloIso_CaloIso_Only");ff_ci_ci->Sumw2();
  ff_ci_r9 = (TH1F*)f_old.Get("ffMet_CaloIso_R9Id_Only");ff_ci_r9->Sumw2();
  ff_r9_ci = (TH1F*)f_old.Get("ffMet_R9Id_CaloIso_Only");ff_r9_ci->Sumw2();
  ff_r9_r9 = (TH1F*)f_old.Get("ffMet_R9Id_R9Id_Only");ff_r9_r9->Sumw2();
  float marksize=0.5;
  ff->SetTitle("");ff->GetXaxis()->SetRangeUser(0,49);ff->SetMarkerSize(marksize);
  ff->Draw();
  x=ff->Integral()/ff_ci_ci->Integral();ff_ci_ci->Scale(x);
  x=ff->Integral()/ff_ci_r9->Integral();ff_ci_r9->Scale(x);
  x=ff->Integral()/ff_r9_ci->Integral();ff_r9_ci->Scale(x);
  x=ff->Integral()/ff_r9_r9->Integral();ff_r9_r9->Scale(x);
  ff_ci_ci->SetLineColor(kBlue);ff_ci_ci->SetMarkerColor(kBlue);ff_ci_ci->SetMarkerSize(marksize);
  ff_ci_r9->SetLineColor(kGreen);ff_ci_r9->SetMarkerColor(kGreen);ff_ci_r9->SetMarkerSize(marksize);
  ff_r9_ci->SetLineColor(kOrange);ff_r9_ci->SetMarkerColor(kOrange);ff_r9_ci->SetMarkerSize(marksize);
  ff_r9_r9->SetLineColor(kRed);ff_r9_r9->SetMarkerColor(kRed);ff_r9_r9->SetMarkerSize(marksize);
  ff_ci_ci->Draw("sames");
  ff_ci_r9->Draw("sames");
  ff_r9_ci->Draw("sames");
  ff_r9_r9->Draw("sames");
  TLegend legMet(.29,.21,.59,.51);
  legMet.SetFillColor(kWhite);
  legMet.SetFillStyle(0);
  legMet.SetBorderSize(0);
  legMet.AddEntry(ff,"ff_all");legMet.AddEntry(ff_ci_ci,"ff_Caloiso_CaloIso");legMet.AddEntry(ff_ci_r9,"ff_Caloiso_R9Id");legMet.AddEntry(ff_r9_ci,"ff_R9Id_CaloIso");legMet.AddEntry(ff_r9_r9,"ff_R9Id_R9Id");
  legMet.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_exclusive_met_raw.png");

  //float top=179.;
  ff_ci_r9->GetYaxis()->SetRangeUser(0,3);
  ff_ci_r9->GetXaxis()->SetRangeUser(0,top);
  TH1F* ffcici = (TH1F*)ff_ci_ci->Clone();TH1F* ffcir9 = (TH1F*)ff_ci_r9->Clone();TH1F* ffr9ci = (TH1F*)ff_r9_ci->Clone();TH1F* ffr9r9 = (TH1F*)ff_r9_r9->Clone();
  //ff_ci_ci->Divide(ff);ff_ci_ci->Draw("PE");
  ff_ci_r9->Divide(ff_ci_ci);ff_ci_r9->Draw("PE");
  ff_r9_ci->Divide(ff_ci_ci);ff_r9_ci->Draw("PEsames");
  ff_r9_r9->Divide(ff_ci_ci);ff_r9_r9->Draw("PEsames");
  TLine line(0,1,top,1);line.Draw("sames");
  TLegend legRat(.23,.53,.58,.83);
  legRat.SetFillColor(kWhite);
  legRat.SetFillStyle(0);
  legRat.SetBorderSize(0);
  legRat.AddEntry(ff_ci_r9,"CaloIso_R9Id / CaloIso_CaloIso");
  legRat.AddEntry(ff_r9_ci,"R9Id_CaloIso / CaloIso_CaloIso");
  legRat.AddEntry(ff_r9_r9,"R9Id_R9Id / CaloIso_CaloIso");
  legRat.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_exclusive_metRatio_raw.png");
  ffcir9->Add(ffr9ci);ffcir9->Add(ffr9r9);
  x = ff_ci_ci->Integral()/ffcir9->Integral();ffcir9->Scale(x);
  ffcir9->Divide(ff_ci_ci);ffcir9->SetLineColor(kViolet);ffcir9->SetMarkerColor(kViolet);ffcir9->GetYaxis()->SetRangeUser(0.8,2);
  ffcir9->Draw();line.Draw("SAMES");
  legRat.Clear();
  legRat.AddEntry(ffcir9,"R9 triggers / CaloIso_CaloIso trigger");
  legRat.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_exclusive_metRatioComb_raw.png");

  //end exclusive raw, now re-weighted

  ff = (TH1F*)f_old.Get("ffMet_reweightJet_binned");ff->Sumw2();
  ff_ci_ci = (TH1F*)f_old.Get("ffMet_CaloIso_CaloIso_reweightJet_binned");ff_ci_ci->Sumw2();
  ff_ci_r9 = (TH1F*)f_old.Get("ffMet_CaloIso_R9Id_reweightJet_binned");ff_ci_r9->Sumw2();
  ff_r9_ci = (TH1F*)f_old.Get("ffMet_R9Id_CaloIso_reweightJet_binned");ff_r9_ci->Sumw2();
  ff_r9_r9 = (TH1F*)f_old.Get("ffMet_R9Id_R9Id_reweightJet_binned");ff_r9_r9->Sumw2();
  float marksize=0.5;
  ff->SetTitle("");ff->GetXaxis()->SetRangeUser(0,49);ff->SetMarkerSize(marksize);
  ff->Draw();
  x=ff->Integral()/ff_ci_ci->Integral();ff_ci_ci->Scale(x);
  x=ff->Integral()/ff_ci_r9->Integral();ff_ci_r9->Scale(x);
  x=ff->Integral()/ff_r9_ci->Integral();ff_r9_ci->Scale(x);
  x=ff->Integral()/ff_r9_r9->Integral();ff_r9_r9->Scale(x);
  ff_ci_ci->SetLineColor(kBlue);ff_ci_ci->SetMarkerColor(kBlue);ff_ci_ci->SetMarkerSize(marksize);
  ff_ci_r9->SetLineColor(kGreen);ff_ci_r9->SetMarkerColor(kGreen);ff_ci_r9->SetMarkerSize(marksize);
  ff_r9_ci->SetLineColor(kOrange);ff_r9_ci->SetMarkerColor(kOrange);ff_r9_ci->SetMarkerSize(marksize);
  ff_r9_r9->SetLineColor(kRed);ff_r9_r9->SetMarkerColor(kRed);ff_r9_r9->SetMarkerSize(marksize);
  ff_ci_ci->Draw("sames");
  ff_ci_r9->Draw("sames");
  ff_r9_ci->Draw("sames");
  ff_r9_r9->Draw("sames");
  TLegend legMet(.29,.21,.59,.51);
  legMet.SetFillColor(kWhite);
  legMet.SetFillStyle(0);
  legMet.SetBorderSize(0);
  legMet.AddEntry(ff,"ff_all");legMet.AddEntry(ff_ci_ci,"ff_Caloiso_CaloIso");legMet.AddEntry(ff_ci_r9,"ff_Caloiso_R9Id");legMet.AddEntry(ff_r9_ci,"ff_R9Id_CaloIso");legMet.AddEntry(ff_r9_r9,"ff_R9Id_R9Id");
  legMet.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_met_reweight.png");

  ff_ci_r9->GetXaxis()->SetRangeUser(0,top);
  //ff_ci_ci->Divide(ff);ff_ci_ci->Draw("PE");
  ff_ci_r9->Divide(ff_ci_ci);ff_ci_r9->Draw("PE");
  ff_r9_ci->Divide(ff_ci_ci);ff_r9_ci->Draw("PEsames");
  ff_r9_r9->Divide(ff_ci_ci);ff_r9_r9->Draw("PEsames");
  TLine line(0,1,top,1);line.Draw("sames");
  TLegend legRat(.23,.53,.58,.83);
  legRat.SetFillColor(kWhite);
  legRat.SetFillStyle(0);
  legRat.SetBorderSize(0);
  legRat.AddEntry(ff_ci_r9,"CaloIso_R9Id / CaloIso_CaloIso");
  legRat.AddEntry(ff_r9_ci,"R9Id_CaloIso / CaloIso_CaloIso");
  legRat.AddEntry(ff_r9_r9,"R9Id_R9Id / CaloIso_CaloIso");
  legRat.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_metRatio_reweight.png");
  //

  gammafake = (TH1F*)f_old.Get("gammafakeMet");gammafake->Sumw2();
  gammafake_ci_ci = (TH1F*)f_old.Get("gammafakeMet_CaloIso_CaloIso");gammafake_ci_ci->Sumw2();
  gammafake_ci_r9 = (TH1F*)f_old.Get("gammafakeMet_CaloIso_R9Id");gammafake_ci_r9->Sumw2();
  gammafake_r9_ci = (TH1F*)f_old.Get("gammafakeMet_R9Id_CaloIso");gammafake_r9_ci->Sumw2();
  gammafake_r9_r9 = (TH1F*)f_old.Get("gammafakeMet_R9Id_R9Id");gammafake_r9_r9->Sumw2();
  float marksize=0.5;
  gammafake->SetTitle("");gammafake->GetXaxis()->SetRangeUser(0,49);gammafake->SetMarkerSize(marksize);
  gammafake->Draw();
  x=gammafake->Integral()/gammafake_ci_ci->Integral();gammafake_ci_ci->Scale(x);
  x=gammafake->Integral()/gammafake_ci_r9->Integral();gammafake_ci_r9->Scale(x);
  x=gammafake->Integral()/gammafake_r9_ci->Integral();gammafake_r9_ci->Scale(x);
  x=gammafake->Integral()/gammafake_r9_r9->Integral();gammafake_r9_r9->Scale(x);
  gammafake_ci_ci->SetLineColor(kBlue);gammafake_ci_ci->SetMarkerColor(kBlue);gammafake_ci_ci->SetMarkerSize(marksize);
  gammafake_ci_r9->SetLineColor(kGreen);gammafake_ci_r9->SetMarkerColor(kGreen);gammafake_ci_r9->SetMarkerSize(marksize);
  gammafake_r9_ci->SetLineColor(kOrange);gammafake_r9_ci->SetMarkerColor(kOrange);gammafake_r9_ci->SetMarkerSize(marksize);
  gammafake_r9_r9->SetLineColor(kRed);gammafake_r9_r9->SetMarkerColor(kRed);gammafake_r9_r9->SetMarkerSize(marksize);
  gammafake_ci_ci->Draw("sames");
  gammafake_ci_r9->Draw("sames");
  gammafake_r9_ci->Draw("sames");
  gammafake_r9_r9->Draw("sames");
  TLegend legMet(.29,.21,.59,.51);
  legMet.SetFillColor(kWhite);
  legMet.SetFillStyle(0);
  legMet.SetBorderSize(0);
  legMet.AddEntry(gammafake,"gammafake_all");legMet.AddEntry(gammafake_ci_ci,"gammafake_Caloiso_CaloIso");legMet.AddEntry(gammafake_ci_r9,"gammafake_Caloiso_R9Id");legMet.AddEntry(gammafake_r9_ci,"gammafake_R9Id_CaloIso");legMet.AddEntry(gammafake_r9_r9,"gammafake_R9Id_R9Id");
  legMet.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_gammafake_met_raw.png");

  gammafake_ci_r9->GetXaxis()->SetRangeUser(0,top);
  //gammafake_ci_ci->Divide(gammafake);gammafake_ci_ci->Draw("PE");
  gammafake_ci_r9->Divide(gammafake_ci_ci);gammafake_ci_r9->Draw("PE");
  gammafake_r9_ci->Divide(gammafake_ci_ci);gammafake_r9_ci->Draw("PEsames");
  gammafake_r9_r9->Divide(gammafake_ci_ci);gammafake_r9_r9->Draw("PEsames");
  TLine line(0,1,top,1);line.Draw("sames");
  TLegend legRat(.23,.53,.58,.83);
  legRat.SetFillColor(kWhite);
  legRat.SetFillStyle(0);
  legRat.SetBorderSize(0);
  legRat.AddEntry(gammafake_ci_r9,"CaloIso_R9Id / CaloIso_CaloIso");
  legRat.AddEntry(gammafake_r9_ci,"R9Id_CaloIso / CaloIso_CaloIso");
  legRat.AddEntry(gammafake_r9_r9,"R9Id_R9Id / CaloIso_CaloIso");
  legRat.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_gammafake_metRatio_raw.png");

  gammafake = (TH1F*)f_old.Get("gammafakeMet_reweightJet_binned");gammafake->Sumw2();
  gammafake_ci_ci = (TH1F*)f_old.Get("gammafakeMet_CaloIso_CaloIso_reweightJet_binned");gammafake_ci_ci->Sumw2();
  gammafake_ci_r9 = (TH1F*)f_old.Get("gammafakeMet_CaloIso_R9Id_reweightJet_binned");gammafake_ci_r9->Sumw2();
  gammafake_r9_ci = (TH1F*)f_old.Get("gammafakeMet_R9Id_CaloIso_reweightJet_binned");gammafake_r9_ci->Sumw2();
  gammafake_r9_r9 = (TH1F*)f_old.Get("gammafakeMet_R9Id_R9Id_reweightJet_binned");gammafake_r9_r9->Sumw2();
  float marksize=0.5;
  gammafake->SetTitle("");gammafake->GetXaxis()->SetRangeUser(0,49);gammafake->SetMarkerSize(marksize);
  gammafake->Draw();
  x=gammafake->Integral()/gammafake_ci_ci->Integral();gammafake_ci_ci->Scale(x);
  x=gammafake->Integral()/gammafake_ci_r9->Integral();gammafake_ci_r9->Scale(x);
  x=gammafake->Integral()/gammafake_r9_ci->Integral();gammafake_r9_ci->Scale(x);
  x=gammafake->Integral()/gammafake_r9_r9->Integral();gammafake_r9_r9->Scale(x);
  gammafake_ci_ci->SetLineColor(kBlue);gammafake_ci_ci->SetMarkerColor(kBlue);gammafake_ci_ci->SetMarkerSize(marksize);
  gammafake_ci_r9->SetLineColor(kGreen);gammafake_ci_r9->SetMarkerColor(kGreen);gammafake_ci_r9->SetMarkerSize(marksize);
  gammafake_r9_ci->SetLineColor(kOrange);gammafake_r9_ci->SetMarkerColor(kOrange);gammafake_r9_ci->SetMarkerSize(marksize);
  gammafake_r9_r9->SetLineColor(kRed);gammafake_r9_r9->SetMarkerColor(kRed);gammafake_r9_r9->SetMarkerSize(marksize);
  gammafake_ci_ci->Draw("sames");
  gammafake_ci_r9->Draw("sames");
  gammafake_r9_ci->Draw("sames");
  gammafake_r9_r9->Draw("sames");
  TLegend legMet(.29,.21,.59,.51);
  legMet.SetFillColor(kWhite);
  legMet.SetFillStyle(0);
  legMet.SetBorderSize(0);
  legMet.AddEntry(gammafake,"gammafake_all");legMet.AddEntry(gammafake_ci_ci,"gammafake_Caloiso_CaloIso");legMet.AddEntry(gammafake_ci_r9,"gammafake_Caloiso_R9Id");legMet.AddEntry(gammafake_r9_ci,"gammafake_R9Id_CaloIso");legMet.AddEntry(gammafake_r9_r9,"gammafake_R9Id_R9Id");
  legMet.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_gammafake_met_reweight.png");

  gammafake_ci_r9->GetXaxis()->SetRangeUser(0,top);
  //gammafake_ci_ci->Divide(gammafake);gammafake_ci_ci->Draw("PE");
  gammafake_ci_r9->Divide(gammafake_ci_ci);gammafake_ci_r9->Draw("PE");
  gammafake_r9_ci->Divide(gammafake_ci_ci);gammafake_r9_ci->Draw("PEsames");
  gammafake_r9_r9->Divide(gammafake_ci_ci);gammafake_r9_r9->Draw("PEsames");
  TLine line(0,1,top,1);line.Draw("sames");
  TLegend legRat(.23,.53,.58,.83);
  legRat.SetFillColor(kWhite);
  legRat.SetFillStyle(0);
  legRat.SetBorderSize(0);
  legRat.AddEntry(gammafake_ci_r9,"CaloIso_R9Id / CaloIso_CaloIso");
  legRat.AddEntry(gammafake_r9_ci,"R9Id_CaloIso / CaloIso_CaloIso");
  legRat.AddEntry(gammafake_r9_r9,"R9Id_R9Id / CaloIso_CaloIso");
  legRat.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_gammafake_metRatio_reweight.png");


  TH1F* ffs_lead = (TH1F*)f_old.Get("ffSigIetaIeta_vs_Met_lead");ffs_lead->Sumw2();
  TH1F* ffs_trail = (TH1F*)f_old.Get("ffSigIetaIeta_vs_Met_trail");ffs_lead->Sumw2();
  
  c1->SetLogz(1);
  ffs_lead->GetXaxis()->SetRangeUser(0,250);
  ffs_lead->GetYaxis()->SetRangeUser(0,.025);
  ffs_lead->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_sihih_lead.png");
  ffs_trail->GetXaxis()->SetRangeUser(0,250);
  ffs_trail->GetYaxis()->SetRangeUser(0,.025);
  ffs_trail->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_sihih_trail.png");

  ffs_lead = (TH1F*)f_old.Get("ffSigIetaIeta_vs_Met_lead_CaloIso_CaloIso_Only");ffs_lead->Sumw2();
  ffs_trail = (TH1F*)f_old.Get("ffSigIetaIeta_vs_Met_trail_CaloIso_CaloIso_Only");ffs_lead->Sumw2();
  ffs_lead->GetXaxis()->SetRangeUser(0,250);
  ffs_lead->GetYaxis()->SetRangeUser(0,.025);
  ffs_lead->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_sihih_CaloIso_CaloIso_Only_lead.png");
  ffs_trail->GetXaxis()->SetRangeUser(0,250);
  ffs_trail->GetYaxis()->SetRangeUser(0,.025);
  ffs_trail->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_sihih_CaloIso_CaloIso_Only_trail.png");

  ffs_lead = (TH1F*)f_old.Get("ffSigIetaIeta_vs_Met_lead_CaloIso_R9Id_Only");ffs_lead->Sumw2();
  ffs_trail = (TH1F*)f_old.Get("ffSigIetaIeta_vs_Met_trail_CaloIso_R9Id_Only");ffs_lead->Sumw2();
  ffs_lead->GetXaxis()->SetRangeUser(0,250);
  ffs_lead->GetYaxis()->SetRangeUser(0,.025);
  ffs_lead->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_sihih_CaloIso_R9Id_Only_lead.png");
  ffs_trail->GetXaxis()->SetRangeUser(0,250);
  ffs_trail->GetYaxis()->SetRangeUser(0,.025);
  ffs_trail->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_sihih_CaloIso_R9Id_Only_trail.png");

  ffs_lead = (TH1F*)f_old.Get("ffSigIetaIeta_vs_Met_lead_R9Id_CaloIso_Only");ffs_lead->Sumw2();
  ffs_trail = (TH1F*)f_old.Get("ffSigIetaIeta_vs_Met_trail_R9Id_CaloIso_Only");ffs_lead->Sumw2();
  ffs_lead->GetXaxis()->SetRangeUser(0,250);
  ffs_lead->GetYaxis()->SetRangeUser(0,.025);
  ffs_lead->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_sihih_R9Id_CaloIso_Only_lead.png");
  ffs_trail->GetXaxis()->SetRangeUser(0,250);
  ffs_trail->GetYaxis()->SetRangeUser(0,.025);
  ffs_trail->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_sihih_R9Id_CaloIso_Only_trail.png");

  ffs_lead = (TH1F*)f_old.Get("ffSigIetaIeta_vs_Met_lead_R9Id_R9Id_Only");ffs_lead->Sumw2();
  ffs_trail = (TH1F*)f_old.Get("ffSigIetaIeta_vs_Met_trail_R9Id_R9Id_Only");ffs_lead->Sumw2();
  ffs_lead->GetXaxis()->SetRangeUser(0,250);
  ffs_lead->GetYaxis()->SetRangeUser(0,.025);
  ffs_lead->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_sihih_R9Id_R9Id_Only_lead.png");
  ffs_trail->GetXaxis()->SetRangeUser(0,250);
  ffs_trail->GetYaxis()->SetRangeUser(0,.025);
  ffs_trail->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_sihih_R9Id_R9Id_Only_trail.png");

  TH2F* ffmetdphi = (TH2F*)f_old.Get("ffMETdPhi_vs_Met");ffmetdphi->Sumw2();
  ffmetdphi->GetXaxis()->SetRangeUser(0,250);
  ffmetdphi->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_metDphi.png");

  ffmetdphi = (TH2F*)f_old.Get("ffMETdPhi_vs_Met_CaloIso_CaloIso_Only");ffmetdphi->Sumw2();
  ffmetdphi->GetXaxis()->SetRangeUser(0,250);
  ffmetdphi->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_metDphi_CaloIso_CaloIso_Only.png");

  ffmetdphi = (TH2F*)f_old.Get("ffMETdPhi_vs_Met_CaloIso_R9Id_Only");ffmetdphi->Sumw2();
  ffmetdphi->GetXaxis()->SetRangeUser(0,250);
  ffmetdphi->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_metDphi_CaloIso_R9Id_Only.png");

  ffmetdphi = (TH2F*)f_old.Get("ffMETdPhi_vs_Met_R9Id_CaloIso_Only");ffmetdphi->Sumw2();
  ffmetdphi->GetXaxis()->SetRangeUser(0,250);
  ffmetdphi->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_metDphi_R9Id_CaloIso_Only.png");

  ffmetdphi = (TH2F*)f_old.Get("ffMETdPhi_vs_Met_R9Id_R9Id_Only");ffmetdphi->Sumw2();
  ffmetdphi->GetXaxis()->SetRangeUser(0,250);
  ffmetdphi->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_metDphi_R9Id_R9Id_Only.png");

  TH2F* ffjetmetdphi = (TH2F*)f_old.Get("ff_JetMETdPhi_vs_Met");ffjetmetdphi->Sumw2();
  ffjetmetdphi->GetXaxis()->SetRangeUser(0,250);
  ffjetmetdphi->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_JetMetDphi.png");

  ffjetmetdphi = (TH2F*)f_old.Get("ff_JetMETdPhi_vs_Met_CaloIso_CaloIso_Only");ffjetmetdphi->Sumw2();
  ffjetmetdphi->GetXaxis()->SetRangeUser(0,250);
  ffjetmetdphi->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_JetMetDphi_CaloIso_CaloIso_Only.png");
  
  ffjetmetdphi = (TH2F*)f_old.Get("ff_JetMETdPhi_vs_Met_CaloIso_R9Id_Only");ffjetmetdphi->Sumw2();
  ffjetmetdphi->GetXaxis()->SetRangeUser(0,250);
  ffjetmetdphi->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_JetMetDphi_CaloIso_R9Id_Only.png");
  
  ffjetmetdphi = (TH2F*)f_old.Get("ff_JetMETdPhi_vs_Met_R9Id_CaloIso_Only");ffjetmetdphi->Sumw2();
  ffjetmetdphi->GetXaxis()->SetRangeUser(0,250);
  ffjetmetdphi->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_JetMetDphi_R9Id_CaloIso_Only.png");
  
  ffjetmetdphi = (TH2F*)f_old.Get("ff_JetMETdPhi_vs_Met_R9Id_R9Id_Only");ffjetmetdphi->Sumw2();
  ffjetmetdphi->GetXaxis()->SetRangeUser(0,250);
  ffjetmetdphi->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtCompare_ff_JetMetDphi_R9Id_R9Id_Only.png");


  TH3F* th3ffall = (TH3F*)f_old.Get("ffSigIetaIeta_1v2_vs_Met");
  TH3F* th3ffcici = (TH3F*)f_old.Get("ffSigIetaIeta_1v2_vs_Met_CaloIso_CaloIso_Only");
  TH3F* th3ffcir9 = (TH3F*)f_old.Get("ffSigIetaIeta_1v2_vs_Met_CaloIso_R9Id_Only");
  TH3F* th3ffr9ci = (TH3F*)f_old.Get("ffSigIetaIeta_1v2_vs_Met_R9Id_CaloIso_Only");
  TH3F* th3ffr9r9 = (TH3F*)f_old.Get("ffSigIetaIeta_1v2_vs_Met_R9Id_R9Id_Only");

  TLine sihline(0,.014,250,.014);sihline.SetLineWidth(2);sihline.SetLineColor(kRed);

  th3ffall->Project3D("xze")->GetXaxis()->SetRangeUser(0,250);
  th3ffall->Project3D("xze")->GetYaxis()->SetRangeUser(.002,.025);
  th3ffall->Project3D("xz")->Draw("COLZ");sihline.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ff_sihih_1v2_vs_Met_lead.png");
  th3ffall->Project3D("yze")->GetXaxis()->SetRangeUser(0,250);
  th3ffall->Project3D("yz")->GetYaxis()->SetRangeUser(.002,.025);
  th3ffall->Project3D("yz")->Draw("COLZ");sihline.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ff_sihih_1v2_vs_Met_trail.png");

  th3ffcir9->Project3D("xze")->GetXaxis()->SetRangeUser(0,250);
  th3ffcir9->Project3D("xze")->GetYaxis()->SetRangeUser(.002,.025);
  th3ffcir9->Project3D("xz")->Draw("COLZ");sihline.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ff_sihih_1v2_vs_Met_lead_CaloIso_R9Id_Only.png");
  th3ffcir9->Project3D("yze")->GetXaxis()->SetRangeUser(0,250);
  th3ffcir9->Project3D("yz")->GetYaxis()->SetRangeUser(.002,.025);
  th3ffcir9->Project3D("yz")->Draw("COLZ");sihline.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ff_sihih_1v2_vs_Met_trail_CaloIso_R9Id_Only.png");
  
  th3ffr9ci->Project3D("xze")->GetXaxis()->SetRangeUser(0,250);
  th3ffr9ci->Project3D("xze")->GetYaxis()->SetRangeUser(.002,.025);
  th3ffr9ci->Project3D("xz")->Draw("COLZ");sihline.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ff_sihih_1v2_vs_Met_lead_R9Id_CaloIso_Only.png");
  th3ffr9ci->Project3D("yze")->GetXaxis()->SetRangeUser(0,250);
  th3ffr9ci->Project3D("yz")->GetYaxis()->SetRangeUser(.002,.025);
  th3ffr9ci->Project3D("yz")->Draw("COLZ");sihline.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ff_sihih_1v2_vs_Met_trail_R9Id_CaloIso_Only.png");
  
  th3ffr9r9->Project3D("xze")->GetXaxis()->SetRangeUser(0,250);
  th3ffr9r9->Project3D("xze")->GetYaxis()->SetRangeUser(.002,.025);
  th3ffr9r9->Project3D("xz")->Draw("COLZ");sihline.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ff_sihih_1v2_vs_Met_lead_R9Id_R9Id_Only.png");
  th3ffr9r9->Project3D("yze")->GetXaxis()->SetRangeUser(0,250);
  th3ffr9r9->Project3D("yz")->GetYaxis()->SetRangeUser(.002,.025);
  th3ffr9r9->Project3D("yz")->Draw("COLZ");sihline.Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ff_sihih_1v2_vs_Met_trail_R9Id_R9Id_Only.png");

  TH1F* ggmetdphil = (TH1F*)f_old.Get("ggMETdPhiLead");
  TH1F* ggmetdphit = (TH1F*)f_old.Get("ggMETdPhiTrail");ggmetdphit->SetLineColor(kRed);
  ggmetdphil->Draw();
  ggmetdphit->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/metdphi_gg.png");

  TH1F* ffmetdphil = (TH1F*)f_old.Get("ffMETdPhiLead");
  TH1F* ffmetdphit = (TH1F*)f_old.Get("ffMETdPhiTrail");ffmetdphit->SetLineColor(kRed);
  ffmetdphil->Draw();
  ffmetdphit->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/metdphi_ff.png");

  TH2F* ffmetdphi2 = (TH2F*)f_old.Get("ffMETdPhi_vs_Met");
  ffmetdphi2->Draw("COLZ");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/metdphi_ff_vsMet.png");

  ffmetdphi2->ProfileY()->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/metdphi_ff_vsMet_pfy.png");

  hist=(TH1F*)f_old.Get("ggMETdPhiLead");float itg=hist->Integral();
  hist->Draw();hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("eeMETdPhiLead");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("ffMETdPhiLead");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("gammafakeMETdPhiLead");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_METdPhiLead.png");

  hist=(TH1F*)f_old.Get("ggMETdPhiTrail");float itg=hist->Integral();
  hist->Draw();
  hist=(TH1F*)f_old.Get("ffMETdPhiTrail");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw();
  hist=(TH1F*)f_old.Get("ggMETdPhiTrail");hist->SetLineColor(kBlack);
  hist->Draw("SAMES");hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("eeMETdPhiTrail");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("gammafakeMETdPhiTrail");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_METdPhiTrail.png");

  hist=(TH1F*)f_old.Get("ggInvarMass");float itg=hist->Integral();hist->GetXaxis()->SetRangeUser(10,250);
  hist->Draw();hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("eeInvarMassFullRange");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);
  hist->Draw("SAMES");
  //hist=(TH1F*)f_old.Get("ggInvarMass");
  //hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("ffInvarMass");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("gammafakeInvarMass");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_InvarMass.png");

  c1->SetLogy(1);/*
  TH1F* hist1 = new TH1F("hist1","",11,-0.5,10.5);hist1->Sumw2();
  TH1F* hist2 = new TH1F("hist2","",11,-0.5,10.5);hist2->Sumw2();
  TH1F* hist3 = new TH1F("hist3","",11,-0.5,10.5);hist3->Sumw2();
  TH1F* hist4 = new TH1F("hist4","",11,-0.5,10.5);hist4->Sumw2();*/
  hist=(TH1F*)f_old.Get("nJets_gg");float itg=hist->Integral();hist->GetXaxis()->SetRangeUser(0,7.9);
  //for(int i=0;i<10;i++){hist1->SetBinContent(i+1,hist->GetBinContent(i*20+1));}
  hist->Draw();hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("nJets_ee");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);
  //for(int i=0;i<10;i++){hist2->SetBinContent(i+1,hist->GetBinContent(i*20+1));}hist2->SetLineColor(kRed);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("nJets_ff");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  //for(int i=0;i<10;i++){hist3->SetBinContent(i+1,hist->GetBinContent(i*20+1));}hist3->SetLineColor(kBlue);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("nJets_gammafake");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  //for(int i=0;i<10;i++){hist4->SetBinContent(i+1,hist->GetBinContent(i*20+1));}hist4->SetLineColor(kGreen);
  hist->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_nJets.png");
  c1->SetLogy(0);

  hist=(TH1F*)f_old.Get("NVertex_gg");float itg=hist->Integral();
  hist->Draw("histo");hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("NVertex_ee");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);
  hist->Draw("histoSAMES");
  hist=(TH1F*)f_old.Get("NVertex_ff");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("histoSAMES");
  hist=(TH1F*)f_old.Get("NVertex_gammafake");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("histoSAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_NVertex.png");

  c1->SetLogy(0);
  hist=(TH1F*)f_old.Get("ggPhotonLessHt");float itg=hist->Integral();hist->GetXaxis()->SetRangeUser(0,2000);
  hist->Draw();hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("eePhotonLessHt");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("ffPhotonLessHt");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("gammafakePhotonLessHt");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_PhotonLessHt.png");

  c1->SetLogz(1);
  TH2F* histth2(0);
  histth2=(TH2F*)f_old.Get("ggPhotonLessHtVsMET");float itg=1.,itgtemp=histth2->Integral();histth2->Scale(itg/itgtemp);histth2->GetXaxis()->SetRangeUser(0,150);histth2->GetYaxis()->SetRangeUser(0,2000);histth2->GetZaxis()->SetRangeUser(2e-7,5e-3);
  histth2->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_PhotonLessHtVsMET_gg.png");
  histth2=(TH2F*)f_old.Get("eePhotonLessHtVsMET");float itg2=histth2->Integral();histth2->Scale(itg/itg2);
  histth2->GetXaxis()->SetRangeUser(0,150);histth2->GetYaxis()->SetRangeUser(0,2000);histth2->GetZaxis()->SetRangeUser(2e-7,5e-3);
  histth2->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_PhotonLessHtVsMET_ee.png");
  histth2=(TH2F*)f_old.Get("ffPhotonLessHtVsMET");float itg2=histth2->Integral();histth2->Scale(itg/itg2);
  histth2->GetXaxis()->SetRangeUser(0,150);histth2->GetYaxis()->SetRangeUser(0,2000);histth2->GetZaxis()->SetRangeUser(2e-7,5e-3);
  histth2->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_PhotonLessHtVsMET_ff.png");
  histth2=(TH2F*)f_old.Get("gammafakePhotonLessHtVsMET");float itg2=histth2->Integral();histth2->Scale(itg/itg2);
  histth2->GetXaxis()->SetRangeUser(0,150);histth2->GetYaxis()->SetRangeUser(0,2000);histth2->GetZaxis()->SetRangeUser(2e-7,5e-3);
  histth2->Draw("colz");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_PhotonLessHtVsMET_gammafake.png");
  c1->SetLogz(0);

  c1->SetLogy(0);
  hist=(TH1F*)f_old.Get("sumEt_gg");float itg=hist->Integral();hist->GetXaxis()->SetRangeUser(0,2000);
  hist->Draw();hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("sumEt_ee");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("sumEt_ff");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("sumEt_gammafake");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_sumEt.png");
  c1->SetLogy(0);

  hist=(TH1F*)f_old.Get("ggdPhi");float itg=hist->Integral();
  hist->Draw();
  hist=(TH1F*)f_old.Get("eedPhi");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);
  hist->Draw("");
  hist=(TH1F*)f_old.Get("ggdPhi");
  hist->Draw("SAMES");hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("ffdPhi");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("gammafakedPhi");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_dPhi.png");
  c1->SetLogy(0);

  hist=(TH1F*)f_old.Get("ggdR");float itg=hist->Integral();
  hist->Draw();
  hist=(TH1F*)f_old.Get("eedR");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);
  hist->Draw("");
  hist=(TH1F*)f_old.Get("ggdR");
  hist->Draw("SAMES");hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("ffdR");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("gammafakedR");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_dR.png");

  hist=(TH1F*)f_old.Get("ggPhi");float itg=hist->Integral();hist->Rebin(3);
  hist->Draw();
  hist=(TH1F*)f_old.Get("ffPhi");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);hist->Rebin(3);
  hist->Draw();
  hist=(TH1F*)f_old.Get("eePhi");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);hist->Rebin(3);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("ggPhi");
  hist->Draw("SAMES");hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("gammafakePhi");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);hist->Rebin(3);
  hist->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_Phi.png");

  hist=(TH1F*)f_old.Get("ggEta");float itg=hist->Integral();hist->Rebin(4);
  hist->Draw();
  hist=(TH1F*)f_old.Get("eeEta");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);hist->Rebin(4);hist->GetXaxis()->SetRangeUser(-1.5,1.5);
  hist->Draw();
  hist=(TH1F*)f_old.Get("ggEta");
  hist->Draw("SAMES");hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("ffEta");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);hist->Rebin(4);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("gammafakeEta");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);hist->Rebin(4);
  hist->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_Eta.png");

  c1->SetLogy(1);
  hist=(TH1F*)f_old.Get("ggPt");float itg=hist->Integral();hist->GetXaxis()->SetRangeUser(0,500);
  hist->Draw();hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("eePt");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("ffPt");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("gammafakePt");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_Pt.png");

  hist=(TH1F*)f_old.Get("ggPtLead");float itg=hist->Integral();hist->GetXaxis()->SetRangeUser(0,500);
  hist->Draw();hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("eePtLead");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("ffPtLead");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("gammafakePtLead");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_PtLead.png");

  hist=(TH1F*)f_old.Get("ggPtTrail");float itg=hist->Integral();hist->GetXaxis()->SetRangeUser(0,500);
  hist->Draw();hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("eePtTrail");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("ffPtTrail");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("gammafakePtTrail");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_PtTrail.png");

  c1->SetLogy(0);
  hist=(TH1F*)f_old.Get("ggMet");float itg=hist->Integral();
  hist->Draw("histo");
  hist=(TH1F*)f_old.Get("eeMet");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);hist->GetXaxis()->SetRangeUser(0,79);
  hist->Draw("histo");
  hist=(TH1F*)f_old.Get("ggMet");
  hist->Draw("histoSAMES");hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)ffMetNew->Clone();float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("histoSAMES");
  hist=(TH1F*)gammafakeMetNew->Clone();float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("histoSAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_Met.png");

  c1->SetLogy(1);
  hist=(TH1F*)f_old.Get("ggDiEMPt");float itg=hist->Integral();
  hist->Draw("histo");
  hist=(TH1F*)f_old.Get("eeDiEMPt");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);hist->GetXaxis()->SetRangeUser(0,220);
  hist->Draw("histo");
  hist=(TH1F*)f_old.Get("ggDiEMPt");
  hist->Draw("histoSAMES");hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("ffDiEMPt");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("histoSAMES");
  hist=(TH1F*)f_old.Get("gammafakeDiEMPt");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("histoSAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_DiEMPt.png");

  hist=(TH1F*)f_old.Get("ggDiJetPt");float itg=hist->Integral();
  hist->Draw("histo");
  hist=(TH1F*)f_old.Get("eeDiJetPt");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);hist->GetXaxis()->SetRangeUser(0,220);
  hist->Draw("histo");
  hist=(TH1F*)f_old.Get("ggDiJetPt");
  hist->Draw("histoSAMES");hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("ffDiJetPt");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("histoSAMES");
  hist=(TH1F*)f_old.Get("gammafakeDiJetPt");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("histoSAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_DiJetPt.png");

  hist=(TH1F*)f_old.Get("ggDiJetPt_0Jet");float itg=hist->Integral();
  hist->Draw("histo");
  hist=(TH1F*)f_old.Get("eeDiJetPt_0Jet");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);hist->GetXaxis()->SetRangeUser(0,220);
  hist->Draw("histo");
  hist=(TH1F*)f_old.Get("ggDiJetPt_0Jet");
  hist->Draw("histoSAMES");hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("ffDiJetPt_0Jet");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("histoSAMES");
  hist=(TH1F*)f_old.Get("gammafakeDiJetPt_0Jet");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("histoSAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_DiJetPt_0Jet.png");

  hist=(TH1F*)f_old.Get("ggDiJetPt_1Jet");float itg=hist->Integral();
  hist->Draw("histo");
  hist=(TH1F*)f_old.Get("eeDiJetPt_1Jet");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);hist->GetXaxis()->SetRangeUser(0,220);
  hist->Draw("histo");
  hist=(TH1F*)f_old.Get("ggDiJetPt_1Jet");
  hist->Draw("histoSAMES");hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("ffDiJetPt_1Jet");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("histoSAMES");
  hist=(TH1F*)f_old.Get("gammafakeDiJetPt_1Jet");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("histoSAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_DiJetPt_1Jet.png");

  hist=(TH1F*)f_old.Get("ggDiJetPt_2Jet");float itg=hist->Integral();
  hist->Draw("histo");
  hist=(TH1F*)f_old.Get("eeDiJetPt_2Jet");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);hist->GetXaxis()->SetRangeUser(0,220);
  hist->Draw("histo");
  hist=(TH1F*)f_old.Get("ggDiJetPt_2Jet");
  hist->Draw("histoSAMES");hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  hist=(TH1F*)f_old.Get("ffDiJetPt_2Jet");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("histoSAMES");
  hist=(TH1F*)f_old.Get("gammafakeDiJetPt_2Jet");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  hist->Draw("histoSAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_DiJetPt_2Jet.png");

  c1->SetLogy(0);
  hist=(TH1F*)f_old.Get("gg_JetMETdPhi");float itg=hist->Integral();hist->GetXaxis()->SetRangeUser(0,500);
  hist->Draw();hist->SetMarkerSize(0.5);hist->Draw("SAMESpe");
  //hist=(TH1F*)f_old.Get("ee_JetMETdPhi");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kRed);
  //hist->Draw("SAMES");
  hist=(TH1F*)f_old.Get("ff_JetMETdPhi");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kBlue);
  hist->Draw("SAMES");
  //hist=(TH1F*)f_old.Get("gammafake_JetMETdPhi");float itg2=hist->Integral();hist->Scale(itg/itg2);hist->SetLineColor(kGreen);
  //hist->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/compare_gg_ee_ff_gf_JetMETdPhi.png");

}

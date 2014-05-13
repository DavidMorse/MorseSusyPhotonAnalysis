#include "TH2F.h"
#include "TH2.h"

float L_int = 19499.;//full dataset
Double_t xbins[]={0,5,10,15,20,25,30,35,40,50,65,80,100,150,250};
int NmetBins=(sizeof(xbins)/sizeof(Double_t))-1;


void WZ(){

 TFile *f_ggHgg = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Higgs_cms533v1_GluGluToHToGG_M-126_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");
  TFile *f_WZHgg = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Higgs_cms533v1_WH_ZH_HToGG_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");
  TFile *f_TTHgg = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Higgs_cms533v1_TTH_HToGG_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");
  TFile *f_VBFHgg = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Higgs_cms533v1_VBF_HToGG_M-126_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");

  TFile fout("SMhiggsBackground_lt2b_0lep_65-Mjj-120.root","RECREATE");

  TCanvas *c1 = new TCanvas("c1","",900,600);
  c1->cd();

  TH2F* h_ggHgg = (TH2F*)f_ggHgg->Get("ggMetVsInvarMass_Loose_lt2b_noMu_noEle_WZjets");if(h_ggHgg->GetSumw2N()==0)h_ggHgg->Sumw2();h_ggHgg->SetMarkerSize(.5);
  TH2F* h_WZHgg = (TH2F*)f_WZHgg->Get("ggMetVsInvarMass_Loose_lt2b_noMu_noEle_WZjets");if(h_WZHgg->GetSumw2N()==0)h_WZHgg->Sumw2();h_WZHgg->SetMarkerSize(.5);
  TH2F* h_TTHgg = (TH2F*)f_TTHgg->Get("ggMetVsInvarMass_Loose_lt2b_noMu_noEle_WZjets");if(h_TTHgg->GetSumw2N()==0)h_TTHgg->Sumw2();h_TTHgg->SetMarkerSize(.5);
  TH2F* h_VBFHgg = (TH2F*)f_VBFHgg->Get("ggMetVsInvarMass_Loose_lt2b_noMu_noEle_WZjets");if(h_VBFHgg->GetSumw2N()==0)h_VBFHgg->Sumw2();h_VBFHgg->SetMarkerSize(.5);



  TH1D* ggHggMet = (TH1D*)h_ggHgg->ProjectionY("ggHggMet",0,-1,"eo");
  TH1D* WZHggMet = (TH1D*)h_WZHgg->ProjectionY("WZHggMet",0,-1,"eo");
  TH1D* TTHggMet = (TH1D*)h_TTHgg->ProjectionY("TTHggMet",0,-1,"eo");
  TH1D* VBFHggMet = (TH1D*)h_VBFHgg->ProjectionY("VBFHggMet",0,-1,"eo");
 
 ggHggMet->Scale((L_int*2.29e-03*19.22)/99989.);//125GeV=/96290);  
  VBFHggMet->Scale((L_int*2.29e-03*1.544)/95677.);//125GeV=/99885); 
  TTHggMet->Scale((L_int*2.29e-03*.1271)/100048.);//125GeV=/100224);
  WZHggMet->Scale((L_int*2.29e-03*(.6782/**(.3257+.014)*/+.3843/**.2*/))/100320);//125 and 126 GeV have same # events
  ggHggMet->SetLineColor(kGreen);ggHggMet->SetMarkerColor(kGreen);ggHggMet->SetFillColor(kGreen);
  WZHggMet->SetLineColor(kCyan);WZHggMet->SetMarkerColor(kCyan);WZHggMet->SetFillColor(kCyan);
  TTHggMet->SetLineColor(kViolet);TTHggMet->SetMarkerColor(kViolet);TTHggMet->SetFillColor(kViolet);
  VBFHggMet->SetLineColor(kRed+3);VBFHggMet->SetMarkerColor(kRed+3);VBFHggMet->SetFillColor(kRed+3);
  //ggHggMet->SetFillStyle(0);WZHggMet->SetFillStyle(0);TTHggMet->SetFillStyle(0);VBFHggMet->SetFillStyle(0);
  ggHggMet->SetLineWidth(2);WZHggMet->SetLineWidth(2);VBFHggMet->SetLineWidth(2);TTHggMet->SetLineWidth(2);

  TH1F* ggHggMetRebin=(TH1F*)ggHggMet->Rebin(NmetBins,"ggHggMetRebin",xbins);
  TH1F* TTHggMetRebin=(TH1F*)TTHggMet->Rebin(NmetBins,"TTHggMetRebin",xbins);
  TH1F* WZHggMetRebin=(TH1F*)WZHggMet->Rebin(NmetBins,"WZHggMetRebin",xbins);
  TH1F* VBFHggMetRebin=(TH1F*)VBFHggMet->Rebin(NmetBins,"VBFHggMetRebin",xbins);
 
  c1->SetLogy(1);

  THStack *HiggsMetStack = new THStack("HiggsMetStack","");
  HiggsMetStack->Add(TTHggMet);HiggsMetStack->Add(WZHggMet);HiggsMetStack->Add(VBFHggMet);HiggsMetStack->Add(ggHggMet);
  HiggsMetStack->Draw("histo");
  HiggsMetStack->GetXaxis()->SetRangeUser(0,249);HiggsMetStack->SetMaximum(1);HiggsMetStack->SetMinimum(3e-4);
  HiggsMetStack->GetYaxis()->SetTitle("Events");HiggsMetStack->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  c1->Print("Plots/Higgs/ggMET_WZ_SMhiggsStack.png");

  THStack *HiggsMetStackRebin = new THStack("HiggsMetStackRebin","");
  HiggsMetStackRebin->Add(TTHggMetRebin);HiggsMetStackRebin->Add(WZHggMetRebin);HiggsMetStackRebin->Add(VBFHggMetRebin);HiggsMetStackRebin->Add(ggHggMetRebin);
  HiggsMetStackRebin->Draw("histo");
  HiggsMetStackRebin->GetXaxis()->SetRangeUser(0,249);HiggsMetStackRebin->SetMaximum(3);HiggsMetStackRebin->SetMinimum(8e-3);
  HiggsMetStackRebin->GetYaxis()->SetTitle("Events");HiggsMetStackRebin->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  c1->Print("Plots/Higgs/ggMET_WZ_SMhiggsRebinStack.png");


  TH1F* SMHiggsMet=(TH1F*)TTHggMet->Clone("SMHiggsMet");SMHiggsMet->Add(VBFHggMet);SMHiggsMet->Add(WZHggMet);SMHiggsMet->Add(ggHggMet);
  SMHiggsMet->SetFillColor(0);SMHiggsMet->SetLineColor(kBlack);SMHiggsMet->SetMarkerColor(kBlack);SMHiggsMet->SetMarkerSize(0.5);
SMHiggsMet->Draw("pe");
  c1->Print("Plots/Higgs/ggMET_WZ_SMhiggs.png");

  fout.cd();
  HiggsMetStack->Write("SMhiggsBackground_lt2b_0lep_65-Mjj-120_Stack");
  HiggsMetStackRebin->Write("SMhiggsBackground_lt2b_0lep_65-Mjj-120_StackRebin");
  SMHiggsMet->Write("SMhiggsBackground_lt2b_0lep_65-Mjj-120_raw");
  //fout.Write();
  fout.Close();




}

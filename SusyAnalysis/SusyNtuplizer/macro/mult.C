#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLine.h"

using namespace std;


void process(TH2F* &h){
  TLine lolo(0,0,0,0);lolo.SetLineColor(kRed);
  TLine lohi(0,0,0,0);lohi.SetLineColor(kRed);
  TLine hilo(0,0,0,0);hilo.SetLineColor(kRed);
  TLine hihi(0,0,0,0);hihi.SetLineColor(kRed);
  TLine siglo(0,0,0,0);siglo.SetLineColor(kBlue);
  TLine sighi(0,0,0,0);sighi.SetLineColor(kBlue);
  int xmax=0,ymax=0,xmin=0,ymin=0;
  h->SetTitle("");
  h->RebinX(2);
  h->RebinY(10);
  for(int i=h->GetXaxis()->FindBin(100);i<h->GetXaxis()->FindBin(200)-1;i++){
    if(h->Integral(i,i,0,-1)>0)xmax=h->GetBinLowEdge(i+1);
  }
  for(int i=0;i<h->GetNbinsY();i++){
    if(h->Integral(h->GetXaxis()->FindBin(100),h->GetXaxis()->FindBin(200)-1,i,i)>0)ymax=h->GetYaxis()->GetBinLowEdge(i+1);
  }
  for(int i=h->GetXaxis()->FindBin(200)-1;i>=h->GetXaxis()->FindBin(100);i--){
    if(h->Integral(i,i,0,-1)>0)xmin=h->GetBinLowEdge(i-1);
  }
  for(int i=h->GetNbinsY();i>0;i--){
    if(h->Integral(0,-1,i,i)>0)ymin=h->GetYaxis()->GetBinLowEdge(i-1);
  }
  //h->GetXaxis()->SetRangeUser(xmin,xmax);h->GetYaxis()->SetRangeUser(ymin,ymax);
  h->GetXaxis()->SetRangeUser(100,200);h->GetYaxis()->SetRangeUser(ymin,ymax);
  h->Draw("colz");
  lolo.DrawLine(103,h->GetBinLowEdge(h->GetYaxis()->FindBin(ymin)+1),103,ymax+h->GetYaxis()->GetBinWidth(0));
  lohi.DrawLine(118,h->GetBinLowEdge(h->GetYaxis()->FindBin(ymin)+1),118,ymax+h->GetYaxis()->GetBinWidth(0));
  hilo.DrawLine(133,h->GetBinLowEdge(h->GetYaxis()->FindBin(ymin)+1),133,ymax+h->GetYaxis()->GetBinWidth(0));
  hihi.DrawLine(163,h->GetBinLowEdge(h->GetYaxis()->FindBin(ymin)+1),163,ymax+h->GetYaxis()->GetBinWidth(0));
  siglo.DrawLine(120,h->GetBinLowEdge(h->GetYaxis()->FindBin(ymin)+1),120,ymax+h->GetYaxis()->GetBinWidth(0));
  sighi.DrawLine(131,h->GetBinLowEdge(h->GetYaxis()->FindBin(ymin)+1),131,ymax+h->GetYaxis()->GetBinWidth(0));
  return;
}

void mult(){

  gStyle->SetOptStat(10011);

  TFile f("hist_HiggsAna_Data2012A_13July_06Aug_recover_Data2012B_13July_Data2012C_Prompt_Runs198022-198903_Runs198941-203742_Data2012D_Prompt_Runs207920-209151_Runs203777-207905_Filter_Ana_HLT_JSON_Two40-25GeVbarrelPhotons_ALL-Higgsino-selectedEvents.root","READ");

  TCanvas *c1 = new TCanvas("c1","",800,600);c1->cd();

  TH2F* Met_1Mu0Ele = (TH2F*)f.Get("ggMetVsInvarMass_Loose_1Mu_0Ele_0_1Jets");
  TH2F* Met_2Mu0Ele = (TH2F*)f.Get("ggMetVsInvarMass_Loose_2Mu_0Ele_0_1Jets");
  TH2F* Met_1Mu1Ele = (TH2F*)f.Get("ggMetVsInvarMass_Loose_1Mu_1Ele_0_1Jets");
  //TH2F* Met_2Mu1Ele = (TH2F*)f.Get("ggMetVsInvarMass_Loose_2Mu_1Ele_0_1Jets");
  TH2F* Met_3Mu0Ele = (TH2F*)f.Get("ggMetVsInvarMass_Loose_3Mu_0Ele_0_1Jets");
  TH2F* Met_0Mu1Ele = (TH2F*)f.Get("ggMetVsInvarMass_Loose_exactly1Ele_0_1Jets");
  TH2F* Met_0Mu2Ele = (TH2F*)f.Get("ggMetVsInvarMass_Loose_exactly2Ele_0_1Jets");
  TH2F* Met_0Mu3Ele = (TH2F*)f.Get("ggMetVsInvarMass_Loose_exactly3Ele_0_1Jets");

  TH2F* MT_1Mu0Ele = (TH2F*)f.Get("ggMTvsInvarMass_Loose_1Mu_0Ele_0_1Jets");
  TH2F* MT_2Mu0Ele = (TH2F*)f.Get("ggMTvsInvarMass_Loose_2Mu_0Ele_0_1Jets");
  TH2F* MT_1Mu1Ele = (TH2F*)f.Get("ggMTvsInvarMass_Loose_1Mu_1Ele_0_1Jets");
  //TH2F* MT_2Mu1Ele = (TH2F*)f.Get("ggMTvsInvarMass_Loose_2Mu_1Ele_0_1Jets");
  TH2F* MT_3Mu0Ele = (TH2F*)f.Get("ggMTvsInvarMass_Loose_3Mu_0Ele_0_1Jets");
  TH2F* MT_0Mu1Ele = (TH2F*)f.Get("ggMTvsInvarMass_Loose_exactly1Ele_0_1Jets");
  TH2F* MT_0Mu2Ele = (TH2F*)f.Get("ggMTvsInvarMass_Loose_exactly2Ele_0_1Jets");
  TH2F* MT_0Mu3Ele = (TH2F*)f.Get("ggMTvsInvarMass_Loose_exactly3Ele_0_1Jets");


  process(Met_1Mu0Ele);
  c1->Print("Plots/Higgs/multiplicities_1Mu0Ele_Met.png");
  c1->Print("Plots/Higgs/multiplicities_1Mu0Ele_Met.pdf");
  process(MT_1Mu0Ele);
  c1->Print("Plots/Higgs/multiplicities_1Mu0Ele_MT.png");
  c1->Print("Plots/Higgs/multiplicities_1Mu0Ele_MT.pdf");

  process(Met_2Mu0Ele);
  c1->Print("Plots/Higgs/multiplicities_2Mu0Ele_Met.png");
  c1->Print("Plots/Higgs/multiplicities_2Mu0Ele_Met.pdf");
  process(MT_2Mu0Ele);
  c1->Print("Plots/Higgs/multiplicities_2Mu0Ele_MT.png");
  c1->Print("Plots/Higgs/multiplicities_2Mu0Ele_MT.pdf");

  process(Met_3Mu0Ele);
  c1->Print("Plots/Higgs/multiplicities_3Mu0Ele_Met.png");
  c1->Print("Plots/Higgs/multiplicities_3Mu0Ele_Met.pdf");
  process(MT_3Mu0Ele);
  c1->Print("Plots/Higgs/multiplicities_3Mu0Ele_MT.png");
  c1->Print("Plots/Higgs/multiplicities_3Mu0Ele_MT.pdf");

  process(Met_1Mu1Ele);
  c1->Print("Plots/Higgs/multiplicities_1Mu1Ele_Met.png");
  c1->Print("Plots/Higgs/multiplicities_1Mu1Ele_Met.pdf");
  process(MT_1Mu1Ele);
  c1->Print("Plots/Higgs/multiplicities_1Mu1Ele_MT.png");
  c1->Print("Plots/Higgs/multiplicities_1Mu1Ele_MT.pdf");
  /*
  process(Met_2Mu1Ele);
  c1->Print("Plots/Higgs/multiplicities_2Mu1Ele_Met.png");
  c1->Print("Plots/Higgs/multiplicities_2Mu1Ele_Met.pdf");
  process(MT_2Mu1Ele);
  c1->Print("Plots/Higgs/multiplicities_2Mu1Ele_MT.png");
  c1->Print("Plots/Higgs/multiplicities_2Mu1Ele_MT.pdf");
  */
  process(Met_0Mu1Ele);
  c1->Print("Plots/Higgs/multiplicities_0Mu1Ele_Met.png");
  c1->Print("Plots/Higgs/multiplicities_0Mu1Ele_Met.pdf");
  process(MT_0Mu1Ele);
  c1->Print("Plots/Higgs/multiplicities_0Mu1Ele_MT.png");
  c1->Print("Plots/Higgs/multiplicities_0Mu1Ele_MT.pdf");

  process(Met_0Mu2Ele);
  c1->Print("Plots/Higgs/multiplicities_0Mu2Ele_Met.png");
  c1->Print("Plots/Higgs/multiplicities_0Mu2Ele_Met.pdf");
  process(MT_0Mu2Ele);
  c1->Print("Plots/Higgs/multiplicities_0Mu2Ele_MT.png");
  c1->Print("Plots/Higgs/multiplicities_0Mu2Ele_MT.pdf");

  process(Met_0Mu3Ele);
  c1->Print("Plots/Higgs/multiplicities_0Mu3Ele_Met.png");
  c1->Print("Plots/Higgs/multiplicities_0Mu3Ele_Met.pdf");
  process(MT_0Mu3Ele);
  c1->Print("Plots/Higgs/multiplicities_0Mu3Ele_MT.png");
  c1->Print("Plots/Higgs/multiplicities_0Mu3Ele_MT.pdf");





}

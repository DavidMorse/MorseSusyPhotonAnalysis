#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "TString.h"
#include "TH2I.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include <map>

using namespace std;

Double_t MTbins[]={0,30,60,90,175};
int nMTbins=(sizeof(MTbins)/sizeof(Double_t))-1;

void AddOverflowToLastBin(TH1F* &hist){
  int overflowbin = hist->GetNbinsX()+1, lastBin=hist->GetNbinsX();
  double overflowErr(0.),lastErr(0.);
  double overflow = hist->IntegralAndError(overflowbin,-1,overflowErr);
  double last = hist->IntegralAndError(lastBin,lastBin,lastErr);
  double lastNew = last+overflow;
  double lastNewErr = sqrt(lastErr*lastErr + overflowErr*overflowErr);
  hist->SetBinContent(lastBin,lastNew);hist->SetBinError(lastBin,lastNewErr);
  hist->SetBinContent(overflowbin,0);hist->SetBinError(overflowbin,0);
  return;
}

void SMS_mt(){

  TFile *finWH=TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Photon_SMS_TChiWH_WincHgg_2J.root","READ");

  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","",1150,700);
 
  c1->cd();
  c1->Draw();
  c1->SetLogz(0);


 //Double_t xbins[]={130,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525};
  Double_t xbins[]={125,140,162.5,187.5,212.5,237.5,262.5,287.5,312.5,337.5,362.5,387.5,412.5,437.5,462.5,487.5,512.5};
  Double_t ybins[]={1,20,25,45,50,70,75,95,100,120,125,145,150,170,175,195,200,220,225,245,250,270,275,295,300,320,325,345,350,370,400};
  int nXbin=(sizeof(xbins)/sizeof(Double_t))-1;
  int nYbin=(sizeof(ybins)/sizeof(Double_t))-1;

  bool ScanBin[17][31];
  for(int i=0;i<17;i++){
    for(int j=0;j<31;j++){
      ScanBin[i][j]=false;
    }
  }


  float centrX[]={130,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500};
  int NcentrX = sizeof(centrX)/sizeof(float);
  int nBinsX=16;
  float low=112.5,high=512.5;

  float L_int=19499.;
  float brFracHgg=0.00229;
  float bf_HH_bbgg=2*0.561*0.00229, bf_HH_WWgg=2*0.231*0.00229, bf_HH_ZZgg=2*0.0289*0.00229, bf_HH_ttgg=2*0.0615*0.00229;
  float bf_WH_Wgg=0.00229, bf_ZH_Zgg=0.00229;
  std::map<float,float> x_secsHiggsino,x_secsWino;

  x_secsHiggsino[130]=3.248;
  x_secsHiggsino[150]=1.876;
  x_secsHiggsino[175]=1.027;
  x_secsHiggsino[200]=0.608;
  x_secsHiggsino[225]=0.377;
  x_secsHiggsino[250]=0.244;
  x_secsHiggsino[275]=0.162;
  x_secsHiggsino[300]=0.111;
  x_secsHiggsino[325]=0.0779;
  x_secsHiggsino[350]=0.0553;
  x_secsHiggsino[375]=0.0401;
  x_secsHiggsino[400]=0.0294;
  x_secsHiggsino[425]=0.0218;
  x_secsHiggsino[450]=0.0163;
  x_secsHiggsino[475]=0.0123;
  x_secsHiggsino[500]=0.0094;

  x_secsWino[130]=5.841;
  x_secsWino[150]=3.385;
  x_secsWino[175]=1.851;
  x_secsWino[200]=1.095;
  x_secsWino[225]=0.68;
  x_secsWino[250]=0.440;
  x_secsWino[275]=0.294;
  x_secsWino[300]=0.199;
  x_secsWino[325]=0.137;
  x_secsWino[350]=0.0967;
  x_secsWino[375]=0.0691;
  x_secsWino[400]=0.0500;
  x_secsWino[425]=0.0366;
  x_secsWino[450]=0.027;
  x_secsWino[475]=0.0201;
  x_secsWino[500]=0.0151;

  int nCats=7;

  TH3F* ggWe  = (TH3F*)finWH->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT");
  TH3F* ggWmu = (TH3F*)finWH->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT");
  TH2F* h_nEvents = (TH2F*)finWH->Get("SMS_mChi_mBino"); 	
  cout<<"before"<<endl;
  static const int xb=ggWe->GetNbinsX(), yb=ggWe->GetNbinsY();
  TH1F* histClone[xb][yb];
  /* for(int i=0;i<ggWe->GetNbinsX()*ggWe->GetNbinsY();i++){
TString 
    histClone[i] = new TH1F();
    }*/
  cout<<"after"<<endl;
  c1->cd();
  for(int i=1;i<=ggWe->GetNbinsX();i++){
    for(int j=1;j<=ggWe->GetNbinsY();j++){
      
      int mChi = ggWe->GetXaxis()->GetBinLowEdge(i);
      int mBino= ggWe->GetYaxis()->GetBinLowEdge(j);
      TString MCHI;ostringstream conv,conv2;conv<<mChi;MCHI=conv.str();
      TString MBINO;conv2<<mBino;MBINO=conv2.str();
      TH1F* RHO = (TH1F*)finWH->Get("rho");
      float nEvents=RHO->GetEntries();
      if(mChi<150)nEvents=120000.;
      else if(mChi<=175)nEvents=60000.;
      else if(mChi<425)nEvents=30000.;
      else if(mChi>=425)nEvents=60000.;
      int chiBin=h_nEvents->GetXaxis()->FindBin(mChi);
      int binoBin=h_nEvents->GetYaxis()->FindBin(mBino);
      float nEvents = h_nEvents->Integral(chiBin,chiBin,binoBin,binoBin);
      if(nEvents==0)continue;
      else ScanBin[i][j]=true;
      TH1F* ggWeMt=(TH1F*)ggWe->ProjectionZ("ggWeMt",i,i,j,j,"o");
      float val=ggWeMt->GetEntries();
      val/=nEvents;
      if(val>0 && val<lowest_gge)lowest_gge=val;
      if(val>highest_gge)highest_gge=val;
      cout<<"  highest_gge: "<<highest_gge<<endl;
      AcceptanceWe->Fill(mChi,mBino,val); 
      float scale = ggWeMt->GetEntries()/ggWeMt->Integral();
      float Norm = x_secsWino[mChi]*bf_WH_Wgg*L_int/nEvents;//wino model, WH but W inclusive so no W->enu branching fraction
      ggWeMt->Scale(scale);ggWeMt->Scale(Norm);
      ggWeMtJECup->Scale(scale);ggWeMtJECup->Scale(Norm);
      ggWeMtJECdown->Scale(scale);ggWeMtJECdown->Scale(Norm);
      TH1F* hist = (TH1F*)ggWeMt->Rebin(nMTbins,"hist",MTbins);
      AddOverflowToLastBin(hist);
      histClone[(i-1)*(j-1)] = (TH1F*)hist->Clone();
      if(i==1 && j==1)histCLone[(i-1)*(j-1)]->Draw("histo");
      else histCLone[(i-1)*(j-1)]->Draw("histo sames");
    }
  }
  c1->Print("Plots/Higgs/mT_SMS_TChiWH_WincHgg_2J_all.png");
  c1->Print("Plots/Higgs/mT_SMS_TChiWH_WincHgg_2J_all.pdf");
}

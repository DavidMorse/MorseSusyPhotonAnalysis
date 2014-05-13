#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "TString.h"
#include "TH2I.h"
#include "TFile.h"
#include <map>

using namespace std;

Double_t MTbins[]={0,30,60,90,180};
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

void HistScale(TH1F* &hist, float scale){
  int nbins = hist->GetNbinsX();
  for(int i=0;i<=nbins+1;i++){
    float val = hist->GetBinContent(i); val*=scale;
    float err = hist->GetBinError(i); err*=scale;
    hist->SetBinContent(i,val);hist->SetBinError(i,err);
  }
  return;
}

void AcceptanceSMS_Tchi_WH_ZH_HH_2J(){

  //TFile finWH("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_SMS_TChiWH_WincHgg_2J_Redo.root","READ");
  TFile *finWH=TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Photon_SMS_TChiWH_WincHgg_2J.root","READ");
  //TFile finWH("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_SMS_TChiWH_WincHgg_2J_3otherTriggers.root","READ");
  //TFile fLimitsSigSMS("results_SMS_TChiWH_WincHgg_2J.root","RECREATE");
  //TFile fLimitsSigSMS("results_aaz.root","RECREATE");
  //TFile fLimitsSigSMS("results_aaWW.root","RECREATE");

  //TFile fout_HH_2W2g_ggEle("SMS_HH_2W2g_ggEle.root");
  //TFile fout_HH_2W2g_ggMu("SMS_HH_2W2g_ggMu.root");



  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","",1150,700);
  TCanvas *c2 = new TCanvas("c2","",1400,700);
  //c1->SetFillColor(0);
  c1->cd();
  c1->Draw();
  c1->SetLogz(0);
  /*  TPad *p1 = new TPad("p1","",.04,.01,.9,1);
      p1->Draw();
      p1->SetLogz(0);
      p1->cd();
  */
  /*
    TText *titleRight = new TText(.98,.8,"Acceptance");
    titleRight->SetTextSize(1);
    //titleRight->SetFillColor(0);
    //titleRight->SetBorderSize(0);
    titleRight->SetTextAngle(90);
    //titleRight->AddText("Acceptance");
    */
  float lowest_gge=999999.,highest_gge=0.;
  float lowest_ggmu=999999.,highest_ggmu=0.;
  float lowest_ggLeftovers=999999.,highest_ggLeftovers=0.;

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
  float PhoEffScale2 = 0.994*0.994;
  std::map<float,float> x_secsHiggsino,x_secsWino;
  /* x_secsAAW[130]=4.12;
     x_secsAAW[150]=2.39;
     x_secsAAW[175]=1.32;
     x_secsAAW[200]=.785;
     x_secsAAW[225]=.491;
     x_secsAAW[250]=.318;
     x_secsAAW[275]=.213;
     x_secsAAW[300]=.146;
     x_secsAAW[325]=.102;
     x_secsAAW[350]=.0737;
     x_secsAAW[375]=.0541;
     x_secsAAW[400]=.0393;
     x_secsAAW[425]=.0296;
     x_secsAAW[450]=.0218;
     x_secsAAW[475]=.0167;
     x_secsAAW[500]=.0125;
  */
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

  vector<float> outValsWe,outValsWmu,outValsWWlnulnu,outValsZnunu,outValsInclusive,outValsLeftovers;
  int nCats=7;

  TH2F* AcceptanceWe = new TH2F("AcceptanceWe","",nXbin,xbins,nYbin,ybins);  			     
  //TH1F* AcceptanceWe = new TH1F("AcceptanceWe","",nBinsX,low,high);  			     
  AcceptanceWe->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceWe->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceWe->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceWe->GetYaxis()->SetTitleOffset(1.);
  AcceptanceWe->GetXaxis()->SetTitleOffset(0.9);
  AcceptanceWe->GetZaxis()->SetTitleOffset(.75);
  AcceptanceWe->GetXaxis()->SetLabelSize(0.05);
  AcceptanceWe->GetYaxis()->SetLabelSize(0.05);
  AcceptanceWe->GetZaxis()->SetLabelSize(0.04);
  //AcceptanceWe->SetLineColor(kBlue); AcceptanceWe->SetFillColor(kBlue);  			     
  TH2F* AcceptanceWmu = new TH2F("AcceptanceWmu","",nXbin,xbins,nYbin,ybins);  			     
  //TH1F* AcceptanceWmu = new TH1F("AcceptanceWmu","",nBinsX,low,high);  			     
  AcceptanceWmu->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceWmu->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceWmu->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceWmu->GetYaxis()->SetTitleOffset(1.);
  AcceptanceWmu->GetXaxis()->SetTitleOffset(.9);
  AcceptanceWmu->GetZaxis()->SetTitleOffset(.75);
  AcceptanceWmu->GetXaxis()->SetLabelSize(0.05);
  AcceptanceWmu->GetYaxis()->SetLabelSize(0.05);
  AcceptanceWmu->GetZaxis()->SetLabelSize(0.04);
  //AcceptanceWmu->SetLineColor(kBlue); AcceptanceWmu->SetFillColor(kBlue); 			     
  TH2F* AcceptanceLeftovers = new TH2F("AcceptanceLeftovers","",nXbin,xbins,nYbin,ybins);  			     
  AcceptanceLeftovers->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceLeftovers->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceLeftovers->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceLeftovers->GetYaxis()->SetTitleOffset(1.);
  AcceptanceLeftovers->GetXaxis()->SetTitleOffset(.9);
  AcceptanceLeftovers->GetZaxis()->SetTitleOffset(.75);
  AcceptanceLeftovers->GetXaxis()->SetLabelSize(0.05);
  AcceptanceLeftovers->GetYaxis()->SetLabelSize(0.05);
  AcceptanceLeftovers->GetZaxis()->SetLabelSize(0.04);

  //ifstream inputfilesAAW;
  //ifstream inputfilesAAZ;
  //ifstream inputfilesAAWW;
  //inputfilesAAW.open("AcceptanceFilesAAW.txt");
  //inputfilesAAZ.open("AcceptanceFilesAAZ.txt");
  //inputfilesAAWW.open("AcceptanceFilesAAWW.txt");

  TH3F* ggWe  = (TH3F*)finWH->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT");
  TH3F* ggWe_noSF  = (TH3F*)finWH->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor");
  TH3F* ggWe_noSF_noPU  = (TH3F*)finWH->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor_noPUweight");
  TH3F* ggWeJECup  = (TH3F*)finWH->Get("SMS_ggMT_JECup_Loose_1Ele_0_1Jets");
  TH3F* ggWeJECdown  = (TH3F*)finWH->Get("SMS_ggMT_JECdown_Loose_1Ele_0_1Jets");
  TH3F* ggWmu = (TH3F*)finWH->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT");
  TH3F* ggWmu_noSF = (TH3F*)finWH->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT_noScaleFactor");
  TH3F* ggWmu_noSF_noPU  = (TH3F*)finWH->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT_noScaleFactor_noPUweight");
  TH3F* ggWmuJECup  = (TH3F*)finWH->Get("SMS_ggMT_JECup_Loose_1Mu_0_1Jets");
  TH3F* ggWmuJECdown  = (TH3F*)finWH->Get("SMS_ggMT_JECdown_Loose_1Mu_0_1Jets");
  TH3F* ggLeftovers = (TH3F*)finWH->Get("gg_SMS_Loose_lt2b_noMu_noEle_lt2jets_mChi_mBino_met");
  TH2F* h_nEvents = (TH2F*)finWH->Get("SMS_mChi_mBino"); 
  ggWe->Sumw2();ggWe_noSF->Sumw2();ggWe_noSF_noPU->Sumw2();ggWeJECup->Sumw2();ggWeJECdown->Sumw2();
  ggWmu->Sumw2();ggWmu_noSF->Sumw2();ggWmu_noSF_noPU->Sumw2();ggWmuJECup->Sumw2();ggWmuJECdown->Sumw2();
  ggLeftovers->Sumw2();h_nEvents->Sumw2();
  h_nEvents->GetXaxis()->SetTitleOffset(0.9);		     
  h_nEvents->GetYaxis()->SetTitleOffset(0.75);		     
  h_nEvents->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  h_nEvents->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  c2->cd();
  h_nEvents->Draw("COLZtext");
  TPaveText *PrelimText = new TPaveText(.39,.86,.89,.9,"NDC");
  PrelimText->AddText("CMS Preliminary                  SMS W(#rightarrowInclusive)H(#gamma#gamma)+2Jet");
  PrelimText->SetFillStyle(0);
  PrelimText->SetFillColor(0);
  PrelimText->SetBorderSize(0);
  PrelimText->Draw();
  TPaveText *EventsText = new TPaveText(.24,.67,.44,.72,"NDC");
  EventsText->AddText("Event Counts");
  ////EventsText->AddText("");
  EventsText->SetFillStyle(0);
  EventsText->SetFillColor(0);
  EventsText->SetBorderSize(0);
  EventsText->Draw();
  c2->Print("Plots/Higgs/AcceptanceEvents_SMS_TChiWH_WincHgg_2J.png");
  c2->Print("Plots/Higgs/AcceptanceEvents_SMS_TChiWH_WincHgg_2J.pdf");
  c1->cd();
  //std::string filename;
  //if(inputfilesAAW.is_open()){
  for(int i=1;i<=ggWe->GetNbinsX();i++){
    for(int j=1;j<=ggWe->GetNbinsY();j++){
      //while(!inputfilesAAW.eof()){
      //std::getline(inputfilesAAW,filename);
      //cout<<"file: "<<filename<<endl;
      //TFile f(filename.c_str(),"READ");
      //f.cd();
      //TString str = f.GetName();
      //int one = str.Index("_chargino");
      //one+=9;
      //int two = str.Index("_",one+1);
      //TString MCHI (str(one,two-one));
    
      //int mChi = MCHI.Atof();
      
      int mChi = ggWe->GetXaxis()->GetBinLowEdge(i);
      int mBino= ggWe->GetYaxis()->GetBinLowEdge(j);
      TString MCHI;ostringstream conv,conv2;conv<<mChi;MCHI=conv.str();//= static_cast<ostringstream*>( &(ostringstream() << mChi) )->str();//sprintf(MCHI,"%i",mChi);
      TString MBINO;conv2<<mBino;MBINO=conv2.str();//sprintf(MBINO,"%i",mBino);
      TH1F* RHO = (TH1F*)finWH->Get("rho");RHO->Sumw2();
      float nEvents=RHO->GetEntries();
      if(mChi<150)nEvents=120000.;
      else if(mChi<=175)nEvents=60000.;
      else if(mChi<425)nEvents=30000.;
      else if(mChi>=425)nEvents=60000.;
      int chiBin=h_nEvents->GetXaxis()->FindBin(mChi);
      int binoBin=h_nEvents->GetYaxis()->FindBin(mBino);
      float nEventsHist = h_nEvents->Integral(chiBin,chiBin,binoBin,binoBin);
      if(nEventsHist==0)continue;
      else ScanBin[i][j]=true;
      //TH2F* ggWe = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Ele_0_1Jets");
      //if(ggWe){
      TH1F* ggWeMt=(TH1F*)ggWe->ProjectionZ("ggWeMt",i,i,j,j,"o");
      TH1F* ggWeMt_noSF=(TH1F*)ggWe_noSF->ProjectionZ("ggWeMt_noSF",i,i,j,j,"o");
      TH1F* ggWeMt_noSF_noPU=(TH1F*)ggWe_noSF_noPU->ProjectionZ("ggWeMt_noSF_noPU",i,i,j,j,"o");
      //TH1F* ggWeMt_noSF=(TH1F*)finWH->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor");ggWeMt_noSF->Sumw2();
      //TH1F* ggWeMt_noSF_noPU=(TH1F*)finWH->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor_noPUweight");ggWeMt_noSF_noPU->Sumw2();
      TH1F* ggWeMtJECup=(TH1F*)ggWeJECup->ProjectionZ("ggWeMtJECup",i,i,j,j,"o");
      TH1F* ggWeMtJECdown=(TH1F*)ggWeJECup->ProjectionZ("ggWeMtJECdown",i,i,j,j,"o");
      float val=ggWeMt->GetEntries();
      cout<<"mChi: "<<mChi<<"  mBino: "<<mBino<<"  nEvents: "<<nEvents<<"  gg We events: "<<val<<endl;
      val/=nEvents;
      //cout<<"highest_gge: "<<highest_gge<<" val: "<<val;
      if(val>0 && val<lowest_gge)lowest_gge=val;
      if(val>highest_gge)highest_gge=val;
      //cout<<"  highest_gge: "<<highest_gge<<endl;
      AcceptanceWe->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      cout<<"raw counts:         "<<ggWeMt->Integral(0,-1)<<endl;
      //float scale = ggWeMt_noSF->GetEntries()/ggWeMt_noSF->Integral(0,-1);
      float scale = ggWeMt_noSF_noPU->Integral(0,-1)/ggWeMt_noSF->Integral(0,-1);
      cout<<"scale:              "<<scale               <<"  counts after scale: "<<ggWeMt->Integral(0,-1)*scale<<endl;
      //cout<<"scale: "<<scale<<endl;
      //float Norm = PhoEffScale2*x_secsHiggsino[mChi]*brFracHgg*L_int/nEventsHist;
      float Norm = PhoEffScale2*x_secsWino[mChi]*bf_WH_Wgg*L_int/nEventsHist;//wino model, WH but W inclusive so no W->enu branching fraction
      cout<<"PhoEffScale2:       "<<PhoEffScale2        <<"  counts after PhoEffScale2: "<<ggWeMt->Integral(0,-1)*scale*PhoEffScale2<<endl;
      cout<<"x_secs["<<mChi<<"]: "<<x_secsWino[mChi]<<"  counts after x_sec: "<<ggWeMt->Integral(0,-1)*scale*PhoEffScale2*x_secsWino[mChi]<<endl;
      cout<<"bf_WH_Wgg:         "<<bf_WH_Wgg          <<"  counts after bf_WH_Wgg: "<<ggWeMt->Integral(0,-1)*scale*PhoEffScale2*x_secsWino[mChi]*bf_WH_Wgg<<endl;
      cout<<"L_int:              "<<L_int               <<"  counts after L_int: "<<ggWeMt->Integral(0,-1)*scale*PhoEffScale2*x_secsWino[mChi]*bf_WH_Wgg*L_int<<endl;
      cout<<"nEventsHist:        "<<nEventsHist         <<"  counts after nEventsHist: "<<ggWeMt->Integral(0,-1)*scale*PhoEffScale2*x_secsWino[mChi]*bf_WH_Wgg*L_int/nEventsHist<<endl;
      cout<<"Norm:               "<<Norm                <<"  counts after Norm: "<<ggWeMt->Integral(0,-1)*scale*Norm<<endl;
      //cout<<"x_secsWino["<<mChi<<"]: "<<x_secsWino[mChi]<<" bf_WH_Wgg: "<<bf_WH_Wgg<<"  nEvents: "<<nEvents<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      //ggWeMt->Scale(scale);ggWeMt->Scale(Norm);
      //ggWeMtJECup->Scale(scale);ggWeMtJECup->Scale(Norm);
      //ggWeMtJECdown->Scale(scale);ggWeMtJECdown->Scale(Norm);
      HistScale(ggWeMt,scale*Norm);HistScale(ggWeMtJECup,scale*Norm);HistScale(ggWeMtJECdown,scale*Norm);
      if(mChi==130 && mBino==1){cout<<"electron cutflow eventsHist: "<<nEventsHist<<endl;cout<<"electron cutflow histScale: "<<scale*Norm<<endl;}
      TH1F* hist = (TH1F*)ggWeMt->Rebin(nMTbins,"hist",MTbins);
      AddOverflowToLastBin(hist);
      TFile fLimitsSigSMS_WH_ggEle("HiggsLimitFiles/LimitPackage_mu_"+MCHI+"_mB_"+MBINO+"_WH_ggEle.root","RECREATE");
      fLimitsSigSMS_WH_ggEle.cd();
      //ggWeMt->Write("h_ggWe_WH_MT_mChi"+MCHI+"_mBino"+MBINO);
      //ggWeMt->Write("hMT_ggEle_MC");
      hist->Write("hMT_ggEle_MC");
      //ggWeMtJECup->Write("hMT_ggEle_JECUp");
      //ggWeMtJECdown->Write("hMT_ggEle_JECDown");
      //cout<<"Total expected events: "<<ggWeMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggWeMt->Integral(0,ggWeMt->FindBin(50))<<endl<<"mt>50 : "<<ggWeMt->Integral(ggWeMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggWeMt->Integral(ggWeMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggWeMt->Integral(ggWeMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggWeMt->Integral(ggWeMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggWeMt->Integral(ggWeMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggWeMt->Integral(ggWeMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggWeMt->Integral(ggWeMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggWeMt->FindBin(50),bin75=ggWeMt->FindBin(75),bin100=ggWeMt->FindBin(100),bin125=ggWeMt->FindBin(125),bin150=ggWeMt->FindBin(150),bin175=ggWeMt->FindBin(175),bin200=ggWeMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsWe.push_back(ggWeMt->Integral(0,-1));
      outValsWe.push_back(ggWeMt->Integral(0,bin50-1));
      outValsWe.push_back(ggWeMt->Integral(bin50,-1));
      outValsWe.push_back(ggWeMt->Integral(bin75,-1));
      outValsWe.push_back(ggWeMt->Integral(bin100,-1));
      //outValsWe.push_back(ggWeMt->Integral(bin125,-1));
      outValsWe.push_back(ggWeMt->Integral(bin150,-1));
      //outValsWe.push_back(ggWeMt->Integral(bin175,-1));
      outValsWe.push_back(ggWeMt->Integral(bin200,-1));
      
      //}
      //delete ggWe;
      fLimitsSigSMS_WH_ggEle.Close();

      //TH2F* ggWmu = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
      //if(ggWmu){
      TH1F* ggWmuMt=(TH1F*)ggWmu->ProjectionZ("ggWmuMt",i,i,j,j,"o");
      TH1F* ggWmuMt_noSF=(TH1F*)ggWmu_noSF->ProjectionZ("ggWmuMt_noSF",i,i,j,j,"o");
      TH1F* ggWmuMt_noSF_noPU=(TH1F*)ggWmu_noSF_noPU->ProjectionZ("ggWmuMt_noSF_noPU",i,i,j,j,"o");
      TH1F* ggWmuMtJECup=(TH1F*)ggWmuJECup->ProjectionZ("ggWmuMtJECup",i,i,j,j,"o");
      TH1F* ggWmuMtJECdown=(TH1F*)ggWmuJECdown->ProjectionZ("ggWmuMtJECdown",i,i,j,j,"o");
      float val=ggWmuMt->GetEntries();
      cout<<"mChi: "<<mChi<<"  mBino: "<<mBino<<"  nEvents: "<<nEventsHist<<"  gg+mu events: "<<val<<endl;
      val/=nEvents;
      //cout<<"highest_ggmu: "<<highest_ggmu<<" val: "<<val;
      if(val>0 && val<lowest_ggmu)lowest_ggmu=val;
      if(val>highest_ggmu)highest_ggmu=val;
      //cout<<"  highest_ggmu: "<<highest_ggmu<<endl<<endl;
      AcceptanceWmu->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      cout<<"raw counts:         "<<ggWmuMt->Integral(0,-1)<<endl;
      //float scale = ggWmuMt_noSF->GetEntries()/ggWmuMt_noSF->Integral(0,-1);
      float scale = ggWmuMt_noSF_noPU->Integral(0,-1)/ggWmuMt_noSF->Integral(0,-1);
      cout<<"scale:              "<<scale               <<"  counts after scale: "<<ggWmuMt->Integral(0,-1)*scale<<endl;
      //cout<<"scale: "<<scale<<endl;
      //float Norm = PhoEffScale2*x_secsWino[mChi]*brFracHgg*L_int/nEventsHist;
      float Norm = PhoEffScale2*x_secsWino[mChi]*bf_WH_Wgg*L_int/nEventsHist;//wino model, WH but W inclusive so no W->munu branching fraction
      cout<<"PhoEffScale2:       "<<PhoEffScale2        <<"  counts after PhoEffScale2: "<<ggWmuMt->Integral(0,-1)*scale*PhoEffScale2<<endl;
      cout<<"x_secs["<<mChi<<"]: "<<x_secsWino[mChi]<<"  counts after x_sec: "<<ggWmuMt->Integral(0,-1)*scale*PhoEffScale2*x_secsWino[mChi]<<endl;
      cout<<"bf_WH_Wgg:         "<<bf_WH_Wgg          <<"  counts after bf_WH_Wgg: "<<ggWmuMt->Integral(0,-1)*scale*PhoEffScale2*x_secsWino[mChi]*bf_WH_Wgg<<endl;
      cout<<"L_int:              "<<L_int               <<"  counts after L_int: "<<ggWmuMt->Integral(0,-1)*scale*PhoEffScale2*x_secsWino[mChi]*bf_WH_Wgg*L_int<<endl;
      cout<<"nEventsHist:        "<<nEventsHist         <<"  counts after nEventsHist: "<<ggWmuMt->Integral(0,-1)*scale*PhoEffScale2*x_secsWino[mChi]*bf_WH_Wgg*L_int/nEventsHist<<endl;
      cout<<"Norm:               "<<Norm                <<"  counts after Norm: "<<ggWmuMt->Integral(0,-1)*scale*Norm<<endl;
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      /*ggWmuMt->Scale(scale);ggWmuMt->Scale(Norm);
      ggWmuMtJECup->Scale(scale);ggWmuMtJECup->Scale(Norm);
      ggWmuMtJECdown->Scale(scale);ggWmuMtJECdown->Scale(Norm);*/
      HistScale(ggWmuMt,scale*Norm);HistScale(ggWmuMtJECup,scale*Norm);HistScale(ggWmuMtJECdown,scale*Norm);
      if(mChi==130 && mBino==1){cout<<"muon cutflow eventsHist: "<<nEventsHist<<endl;cout<<"muon cutflow histScale: "<<scale*Norm<<endl;}
      TH1F* hist = (TH1F*)ggWmuMt->Rebin(nMTbins,"hist",MTbins);
      AddOverflowToLastBin(hist);
      TFile fLimitsSigSMS_WH_ggMu("HiggsLimitFiles/LimitPackage_mu_"+MCHI+"_mB_"+MBINO+"_WH_ggMu.root","RECREATE");
      fLimitsSigSMS_WH_ggMu.cd();
      //ggWmuMt->Write("h_ggWmu_WH_MT_mChi"+MCHI+"_mBino"+MBINO);
      //ggWmuMt->Write("hMT_ggMu_MC");
      hist->Write("hMT_ggMu_MC");
      TH1F* histClone = (TH1F*)hist->Clone("histClone");
      float x = 1./histClone->Integral(0,-1);histClone->Scale(x);
      if(i==1 && j==1)histClone->Draw("histo");
      else histClone->Draw("histosames");
      if(i==ggWe->GetNbinsX() && j==ggWe->GetNbinsY()){c1->Print("Plots/Higgs/mT_SMS_TChiWH_WincHgg_2J_all.pdf");
	c1->Print("Plots/Higgs/mT_SMS_TChiWH_WincHgg_2J_all.png");}
      //ggWmuMtJECup->Write("hMT_ggMu_JECUp");
      //ggWmuMtJECdown->Write("hMT_ggMu_JECDown");
      //cout<<"Total expected events: "<<ggWmuMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggWmuMt->Integral(0,ggWmuMt->FindBin(50))<<endl<<"mt>50 : "<<ggWmuMt->Integral(ggWmuMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggWmuMt->Integral(ggWmuMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggWmuMt->Integral(ggWmuMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggWmuMt->Integral(ggWmuMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggWmuMt->Integral(ggWmuMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggWmuMt->Integral(ggWmuMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggWmuMt->Integral(ggWmuMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggWmuMt->FindBin(50),bin75=ggWmuMt->FindBin(75),bin100=ggWmuMt->FindBin(100),bin125=ggWmuMt->FindBin(125),bin150=ggWmuMt->FindBin(150),bin175=ggWmuMt->FindBin(175),bin200=ggWmuMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsWmu.push_back(ggWmuMt->Integral(0,-1));
      outValsWmu.push_back(ggWmuMt->Integral(0,bin50-1));
      outValsWmu.push_back(ggWmuMt->Integral(bin50,-1));
      outValsWmu.push_back(ggWmuMt->Integral(bin75,-1));
      outValsWmu.push_back(ggWmuMt->Integral(bin100,-1));
      //outValsWmu.push_back(ggWmuMt->Integral(bin125,-1));
      outValsWmu.push_back(ggWmuMt->Integral(bin150,-1));
      //outValsWmu.push_back(ggWmuMt->Integral(bin175,-1));
      outValsWmu.push_back(ggWmuMt->Integral(bin200,-1));
      // }
      // delete ggWmu;
      
      
      //f.Close();
      fLimitsSigSMS_WH_ggMu.Close();

      TH1F* ggLeftoversMt=(TH1F*)ggLeftovers->ProjectionZ("ggLeftoversMt",i,i,j,j,"o");
      float val=ggLeftoversMt->GetEntries();
      val/=nEvents;
      //cout<<"highest_ggLeftovers: "<<highest_ggLeftovers<<" val: "<<val;
      if(val>0 && val<lowest_ggLeftovers)lowest_ggLeftovers=val;
      if(val>highest_ggLeftovers)highest_ggLeftovers=val;
      //cout<<"  highest_ggLeftovers: "<<highest_ggLeftovers<<endl<<endl;
      AcceptanceLeftovers->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      //fLimitsSigSMS.cd();
      //ggLeftovers->Write("h_ggLeftovers_WH_MT_mChi"+MCHI+"_mBino"+MBINO);
      float scale = ggLeftoversMt->GetEntries()/ggLeftoversMt->Integral(0,-1);
      //cout<<"scale: "<<scale<<endl;
      float Norm = PhoEffScale2*x_secsWino[mChi]*brFracHgg*L_int/nEventsHist;
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      /*ggLeftoversMt->Scale(scale);
	ggLeftoversMt->Scale(Norm);*/
      HistScale(ggLeftoversMt,scale*Norm);
      //cout<<"Total expected events: "<<ggLeftoversMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggLeftoversMt->Integral(0,ggLeftoversMt->FindBin(50))<<endl<<"mt>50 : "<<ggLeftoversMt->Integral(ggLeftoversMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggLeftoversMt->Integral(ggLeftoversMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggLeftoversMt->Integral(ggLeftoversMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggLeftoversMt->Integral(ggLeftoversMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggLeftoversMt->Integral(ggLeftoversMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggLeftoversMt->Integral(ggLeftoversMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggLeftoversMt->Integral(ggLeftoversMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggLeftoversMt->FindBin(50),bin75=ggLeftoversMt->FindBin(75),bin100=ggLeftoversMt->FindBin(100),bin125=ggLeftoversMt->FindBin(125),bin150=ggLeftoversMt->FindBin(150),bin175=ggLeftoversMt->FindBin(175),bin200=ggLeftoversMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsLeftovers.push_back(ggLeftoversMt->Integral(0,-1));
      outValsLeftovers.push_back(ggLeftoversMt->Integral(0,bin50-1));
      outValsLeftovers.push_back(ggLeftoversMt->Integral(bin50,-1));
      outValsLeftovers.push_back(ggLeftoversMt->Integral(bin75,-1));
      outValsLeftovers.push_back(ggLeftoversMt->Integral(bin100,-1));
      //outValsLeftovers.push_back(ggLeftoversMt->Integral(bin125,-1));
      outValsLeftovers.push_back(ggLeftoversMt->Integral(bin150,-1));
      //outValsLeftovers.push_back(ggLeftoversMt->Integral(bin175,-1));
      outValsLeftovers.push_back(ggLeftoversMt->Integral(bin200,-1));
      // }
      // delete ggLeftovers;
      
      
      //f.Close();
    }
  }
  /*
    if(inputfilesAAZ.is_open()){
    
    while(!inputfilesAAZ.eof()){
    std::getline(inputfilesAAZ,filename);
    //cout<<"file: "<<filename<<endl;
    TFile f(filename.c_str(),"READ");
    f.cd();
    TString str = f.GetName();
    int one = str.Index("_chargino");
    one+=9;
    int two = str.Index("_",one+1);
    TString MCHI (str(one,two-one));
    
    int mChi = MCHI.Atof();
      
    TH1F* RHO = (TH1F*)f->Get("rho");
    float nEvents=RHO->GetEntries();
   
    TH2F* ggZnunu = (TH2F*)f->Get("ggMtVsInvarMass_Loose_0Lep_0_1Jets");
    if(ggZnunu){
    float val=ggZnunu->GetEntries();
    TH1F* ggZnunuMt=(TH1F*)ggZnunu->ProjectionY("ggZnunuMt");
    val/=nEvents;
    if(val<lowest)lowest=val;
    if(val>highest)highest=val;
    AcceptanceZnunu->Fill(mChi,val); 
    //cout<<mS<<endl<<mG<<endl<<val<<endl;
    fLimitsSigAAZ.cd();
    ggZnunu->Write("h_ggZnunu_WH_MT_mChi"+MCHI);
    float scale = ggZnunuMt->GetEntries()/ggZnunuMt->Integral(0,-1);
    //cout<<"scale: "<<scale<<endl;
    float Norm = PhoEffScale2*x_secsAAZ[mChi]*brFracHgg*L_int/nEventsHist;
    //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
    //cout<<"Norm: "<<Norm<<endl;
    ggZnunuMt->Scale(scale);
    ggZnunuMt->Scale(Norm);
    //cout<<"Total expected events: "<<ggZnunuMt->Integral(0,-1)<<endl;
    //cout<<"mt<50 : "<<ggZnunuMt->Integral(0,ggZnunuMt->FindBin(50))<<endl<<"mt>50 : "<<ggZnunuMt->Integral(ggZnunuMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggZnunuMt->Integral(ggZnunuMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggZnunuMt->Integral(ggZnunuMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggZnunuMt->Integral(ggZnunuMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggZnunuMt->Integral(ggZnunuMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggZnunuMt->Integral(ggZnunuMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggZnunuMt->Integral(ggZnunuMt->FindBin(200),-1)<<endl<<endl;

    float bin50=ggZnunuMt->FindBin(50),bin75=ggZnunuMt->FindBin(75),bin100=ggZnunuMt->FindBin(100),bin125=ggZnunuMt->FindBin(125),bin150=ggZnunuMt->FindBin(150),bin175=ggZnunuMt->FindBin(175),bin200=ggZnunuMt->FindBin(200);

    //if add or remove any, must change nCats at top
    outValsZnunu.push_back(ggZnunuMt->Integral(0,-1));
    outValsZnunu.push_back(ggZnunuMt->Integral(0,bin50-1));
    outValsZnunu.push_back(ggZnunuMt->Integral(bin50,-1));
    outValsZnunu.push_back(ggZnunuMt->Integral(bin75,-1));
    outValsZnunu.push_back(ggZnunuMt->Integral(bin100,-1));
    //outValsZnunu.push_back(ggZnunuMt->Integral(bin125,-1));
    outValsZnunu.push_back(ggZnunuMt->Integral(bin150,-1));
    //outValsZnunu.push_back(ggZnunuMt->Integral(bin175,-1));
    outValsZnunu.push_back(ggZnunuMt->Integral(bin200,-1));
    }
    delete ggZnunu;
  
    f.Close();
    }
    }

    if(inputfilesAAWW.is_open()){
    
    while(!inputfilesAAWW.eof()){
    std::getline(inputfilesAAWW,filename);
    //cout<<"file: "<<filename<<endl;
    TFile f(filename.c_str(),"READ");
    f.cd();
    TString str = f.GetName();
    int one = str.Index("_chargino");
    one+=9;
    int two = str.Index("_",one+1);
    TString MCHI (str(one,two-one));
    
    int mChi = MCHI.Atof();
      
    TH1F* RHO = (TH1F*)f->Get("rho");
    float nEvents=RHO->GetEntries();
    TH2F* ggWWlnulnu = (TH2F*)f->Get("ggMtVsInvarMass_Loose_2LepOSoffZ_0_1Jets");
    if(ggWWlnulnu){
    float val=ggWWlnulnu->GetEntries();
    TH1F* ggWWlnulnuMt=(TH1F*)ggWWlnulnu->ProjectionY("ggWWlnulnuMt");
    val/=nEvents;
    if(val<lowest)lowest=val;
    if(val>highest)highest=val;
    AcceptanceWWlnulnu->Fill(mChi,val); 
    //cout<<mS<<endl<<mG<<endl<<val<<endl;
    fLimitsSigSMS.cd();
    ggWWlnulnu->Write("h_ggWWlnulnu_WH_MT_mChi"+MCHI);
    float scale = ggWWlnulnuMt->GetEntries()/ggWWlnulnuMt->Integral(0,-1);
    //cout<<"scale: "<<scale<<endl;
    float Norm = PhoEffScale2*x_secsAAWW[mChi]*brFracHgg*L_int/nEventsHist;
    //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
    //cout<<"Norm: "<<Norm<<endl;
    ggWWlnulnuMt->Scale(scale);
    ggWWlnulnuMt->Scale(Norm);
    //cout<<"Total expected events: "<<ggWWlnulnuMt->Integral(0,-1)<<endl;
    //cout<<"mt<50 : "<<ggWWlnulnuMt->Integral(0,ggWWlnulnuMt->FindBin(50))<<endl<<"mt>50 : "<<ggWWlnulnuMt->Integral(ggWWlnulnuMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggWWlnulnuMt->Integral(ggWWlnulnuMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggWWlnulnuMt->Integral(ggWWlnulnuMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggWWlnulnuMt->Integral(ggWWlnulnuMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggWWlnulnuMt->Integral(ggWWlnulnuMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggWWlnulnuMt->Integral(ggWWlnulnuMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggWWlnulnuMt->Integral(ggWWlnulnuMt->FindBin(200),-1)<<endl<<endl;

    float bin50=ggWWlnulnuMt->FindBin(50),bin75=ggWWlnulnuMt->FindBin(75),bin100=ggWWlnulnuMt->FindBin(100),bin125=ggWWlnulnuMt->FindBin(125),bin150=ggWWlnulnuMt->FindBin(150),bin175=ggWWlnulnuMt->FindBin(175),bin200=ggWWlnulnuMt->FindBin(200);

    //if add or remove any, must change nCats at top
    outValsWWlnulnu.push_back(ggWWlnulnuMt->Integral(0,-1));
    outValsWWlnulnu.push_back(ggWWlnulnuMt->Integral(0,bin50-1));
    outValsWWlnulnu.push_back(ggWWlnulnuMt->Integral(bin50,-1));
    outValsWWlnulnu.push_back(ggWWlnulnuMt->Integral(bin75,-1));
    outValsWWlnulnu.push_back(ggWWlnulnuMt->Integral(bin100,-1));
    //outValsWWlnulnu.push_back(ggWWlnulnuMt->Integral(bin125,-1));
    outValsWWlnulnu.push_back(ggWWlnulnuMt->Integral(bin150,-1));
    //outValsWWlnulnu.push_back(ggWWlnulnuMt->Integral(bin175,-1));
    outValsWWlnulnu.push_back(ggWWlnulnuMt->Integral(bin200,-1));
    }
    delete ggWWlnulnu;
    
    f.Close();
    }
    }

  */


  /*
    for(int i=1;i<AcceptanceWe->GetNbinsX();i++){
    for(int j=1;j<AcceptanceWe->GetNbinsY();j++){
    int bin = AcceptanceWe->GetBin(i,j);
    if(AcceptanceWe->GetBinContent(bin)==0){
    int binRight = AcceptanceWe->GetBin(i+1,j);
    int binLeft = AcceptanceWe->GetBin(i-1,j);
    int binUp = AcceptanceWe->GetBin(i,j+1);
    int binDown = AcceptanceWe->GetBin(i,j-1);
    float SchmearedVal = (AcceptanceWe->GetBinContent(binRight)+AcceptanceWe->GetBinContent(binLeft)+AcceptanceWe->GetBinContent(binUp)+AcceptanceWe->GetBinContent(binDown))/4.;
    AcceptanceWe->SetBinContent(bin,SchmearedVal);
    //cout<<"No Jet Req Bin "<<bin<<" has 0 value , mS="<<AcceptanceWe->GetXaxis()->GetBinLowEdge(i)<<" mG="<<AcceptanceWe->GetYaxis()->GetBinLowEdge(j)<<"  Value set to:"<<SchmearedVal<<endl;
    }
    }
    }
  */
  /*
  //highest_gge=.012;
  //AcceptanceWe->SetMinimum(0);
  //AcceptanceWe->SetMaximum(highest_gge*1.1);
  AcceptanceWe->GetZaxis()->SetRangeUser(lowest_gge,highest_gge*1.1);//SetMaximum(highest_gge*1.1);
  //AcceptanceWe->SetMinimum(.22);
  AcceptanceWe->Draw("COLZ");
  TPaveText *PrelimText2 = new TPaveText(.25,.86,.87,.9,"NDC");
  PrelimText2->AddText("CMS Preliminary                  SMS W(#rightarrowInclusive)H(#gamma#gamma)+2Jet");
  PrelimText2->SetFillStyle(0);
  PrelimText2->SetFillColor(0);
  PrelimText2->SetBorderSize(0);
  PrelimText2->Draw();
  PrelimText2->Draw();
  TPaveText *WeText = new TPaveText(.24,.67,.64,.72,"NDC");
  WeText->AddText("#gamma#gamma+electron, #leq1 b-jet");
  ////WeText->AddText("");
  WeText->SetFillStyle(0);
  WeText->SetFillColor(0);
  WeText->SetBorderSize(0);
  WeText->Draw();
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_gge.png");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_gge.pdf");
  
  //c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_gge_3otherTrigs.png");
  //c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_gge_3otherTrigs.pdf");
  
  //AcceptanceWmu->SetMinimum(0);
  //AcceptanceWmu->SetMaximum(highest_ggmu*1.1);
  AcceptanceWmu->GetZaxis()->SetRangeUser(lowest_ggmu,highest_ggmu*1.1);//SetMaximum(highest_ggmu*1.1);
  //AcceptanceWmu->SetMinimum(.22);
  AcceptanceWmu->Draw("COLZ");
  PrelimText2->Draw();
  TPaveText *WmuText = new TPaveText(.24,.67,.64,.72,"NDC");
  WmuText->AddText("#gamma#gamma+#mu, #leq1 b-jet");
  ////WmuText->AddText("");
  WmuText->SetFillStyle(0);
  WmuText->SetFillColor(0);
  WmuText->SetBorderSize(0);
  WmuText->Draw();
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_ggmu.png");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_ggmu.pdf");
  */
  /*
  //c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_ggmu_3otherTrigs.png");
  //c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_ggmu_3otherTrigs.pdf");
  
  AcceptanceZnunu->SetMinimum(0);
  AcceptanceZnunu->SetMaximum(highest*1.1);
  //AcceptanceZnunu->SetMinimum(.22);
  AcceptanceZnunu->Draw("histo");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_Znunu.png");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_Znunu.pdf");
  
  //c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_Znunu_3otherTrigs.png");
  //c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_Znunu_3otherTrigs.pdf");
  
  AcceptanceWWlnulnu->SetMinimum(0);
  AcceptanceWWlnulnu->SetMaximum(highest*1.1);
  //AcceptanceWWlnulnu->SetMinimum(.22);
  AcceptanceWWlnulnu->Draw("histo");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_WWlnulnu.png");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_WWlnulnu.pdf");
  
  //c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_WWlnulnu_3otherTrigs.png");
  //c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_WWlnulnu_3otherTrigs.pdf");
  */
  int j=0,k=0;


  cout<<"\\documentclass[11pt]{article}"<<endl;
  cout<<"\\usepackage{calc}"<<endl;
  cout<<"\\usepackage{multirow}"<<endl;
  cout<<"\\usepackage{verbatim}"<<endl;
  cout<<"\\usepackage{changepage}"<<endl;
  cout<<"\\usepackage{tabularx}"<<endl;
  cout<<"\\begin{document}"<<endl;

  cout<<" \\begin{center}"<<endl;
  cout<<"  \\begin{small}"<<endl;
  cout<<"  \\hspace*{-3cm}"<<endl;
  cout<<"   \\begin{tabularx}{1.45\\textwidth}{ | c | c | c | c | c | c | c | c | c | c |}\n      \\hline"<<endl;
  cout<<"%\n%\n%\n";
  //cout<<"      \\bf{SMS WH e+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$125} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$175} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
  cout<<"      \\bf{SMS WH e+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
  
  for(int i=0;i<outValsWe.size();i++){
    if(i>0 && i%nCats==0)cout<<" \\\\ \\hline"<<endl;
    if(i%nCats==0 && k<=NcentrX){cout<<"      chargino"<<centrX[k];k++;}
    cout<<" & "<<outValsWe[i];
  }
  cout<<" \\\\ \\hline"<<endl;
  cout<<"%\n%\n%\n";
  cout<<"   \\end{tabularx}"<<endl;
  cout<<"  \\end{small}"<<endl;
  cout<<"%\n%\n";

  cout<<"  \\begin{small}"<<endl;
  cout<<"  \\hspace*{-3cm}"<<endl;
  cout<<"   \\begin{tabularx}{1.45\\textwidth}{ | c | c | c | c | c | c | c | c | c | c |}\n      \\hline"<<endl;
  cout<<"%\n%\n%\n";
  //cout<<"      \\bf{SMS WH $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$125} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$175} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
  cout<<"      \\bf{SMS WH $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
  k=0;
  for(int i=0;i<outValsWmu.size();i++){
    if(i>0 && i%nCats==0)cout<<" \\\\ \\hline"<<endl;
    if(i%nCats==0 && k<=NcentrX){cout<<"      chargino"<<centrX[k];k++;}
    cout<<" & "<<outValsWmu[i];
  }
  cout<<" \\\\ \\hline"<<endl;
  cout<<"%\n%\n%\n";
  cout<<"   \\end{tabularx}"<<endl;
  cout<<"  \\end{small}"<<endl;
  cout<<"%\n%\n";

  cout<<"  \\begin{small}"<<endl;
  cout<<"  \\hspace*{-3cm}"<<endl;
  cout<<"   \\begin{tabularx}{1.45\\textwidth}{ | c | c | c | c | c | c | c | c | c | c |}\n      \\hline"<<endl;
  cout<<"%\n%\n%\n";
  //cout<<"      \\bf{SMS WH $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$125} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$175} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
  cout<<"      \\bf{SMS WH $leftovers$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
  k=0;
  for(int i=0;i<outValsLeftovers.size();i++){
    if(i>0 && i%nCats==0)cout<<" \\\\ \\hline"<<endl;
    if(i%nCats==0 && k<=NcentrX){cout<<"      chargino"<<centrX[k];k++;}
    cout<<" & "<<outValsLeftovers[i];
  }
  cout<<" \\\\ \\hline"<<endl;
  cout<<"%\n%\n%\n";
  cout<<"   \\end{tabularx}"<<endl;
  cout<<"  \\end{small}"<<endl;
  cout<<"%\n%\n";
  /*
    cout<<"  \\begin{small}"<<endl;
    cout<<"  \\hspace*{-3cm}"<<endl;
    cout<<"   \\begin{tabularx}{1.45\\textwidth}{ | c | c | c | c | c | c | c | c | c | c |}\n      \\hline"<<endl;
    cout<<"%\n%\n%\n";
    //cout<<"      \\bf{SMS WZ $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$125} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$175} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
    cout<<"      \\bf{SMS WZ $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
    k=0;
    for(int i=0;i<outValsZnunu.size();i++){
    if(i>0 && i%nCats==0)cout<<" \\\\ \\hline"<<endl;
    if(i%nCats==0){cout<<"      chargino"<<centrX[k];k++;}
    cout<<" & "<<outValsZnunu[i];
    }
    cout<<" \\\\ \\hline"<<endl;
    cout<<"%\n%\n%\n";
    cout<<"   \\end{tabularx}"<<endl;
    cout<<"  \\end{small}"<<endl;
    cout<<"%\n%\n";

    cout<<"  \\begin{small}"<<endl;
    cout<<"  \\hspace*{-3cm}"<<endl;
    cout<<"   \\begin{tabularx}{1.45\\textwidth}{ | c | c | c | c | c | c | c | c | c | c |}\n      \\hline"<<endl;
    cout<<"%\n%\n%\n";
    //cout<<"      \\bf{SMS WZ $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$125} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$175} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
    cout<<"      \\bf{SMS WZ $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
    k=0;
    for(int i=0;i<outValsWWlnulnu.size();i++){
    if(i>0 && i%nCats==0)cout<<" \\\\ \\hline"<<endl;
    if(i%nCats==0){cout<<"      chargino"<<centrX[k];k++;}
    cout<<" & "<<outValsWWlnulnu[i];
    }
    cout<<" \\\\ \\hline"<<endl;
    cout<<"%\n%\n%\n";
    cout<<"   \\end{tabularx}"<<endl;
    cout<<"  \\end{small}"<<endl;
    cout<<"%\n%\n";
  */

  cout<<" \\end{center}"<<endl;
  cout<<"\\end{document}"<<endl;  

  //now do ZH

  //TFile finZH("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_SMS_TChiZH_ZincHgg_2J.root","READ");
  TFile *finZH=TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Photon_SMS_TChiZH_ZincHgg_2J.root","READ");

  //lowest=999999.;highest=0.;

 
  vector<float> outValsWe_ZH,outValsWmu_ZH,outValsLeftovers_ZH;
  int nCats=7;

  TH2F* AcceptanceWe_ZH = new TH2F("AcceptanceWe_ZH","",nXbin,xbins,nYbin,ybins);  			     
  AcceptanceWe_ZH->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceWe_ZH->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceWe_ZH->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceWe_ZH->GetYaxis()->SetTitleOffset(1.);
  AcceptanceWe_ZH->GetXaxis()->SetTitleOffset(0.9);
  AcceptanceWe_ZH->GetZaxis()->SetTitleOffset(.75);
  AcceptanceWe_ZH->GetXaxis()->SetLabelSize(0.05);
  AcceptanceWe_ZH->GetYaxis()->SetLabelSize(0.05);
  AcceptanceWe_ZH->GetZaxis()->SetLabelSize(0.04);
  TH2F* AcceptanceWmu_ZH = new TH2F("AcceptanceWmu_ZH","",nXbin,xbins,nYbin,ybins);  			     
  AcceptanceWmu_ZH->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceWmu_ZH->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceWmu_ZH->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceWmu_ZH->GetYaxis()->SetTitleOffset(1.);
  AcceptanceWmu_ZH->GetXaxis()->SetTitleOffset(.9);
  AcceptanceWmu_ZH->GetZaxis()->SetTitleOffset(.75);
  AcceptanceWmu_ZH->GetXaxis()->SetLabelSize(0.05);
  AcceptanceWmu_ZH->GetYaxis()->SetLabelSize(0.05);
  AcceptanceWmu_ZH->GetZaxis()->SetLabelSize(0.04);			     
  TH2F* AcceptanceLeftovers_ZH = new TH2F("AcceptanceLeftovers_ZH","",nXbin,xbins,nYbin,ybins);  			     
  AcceptanceLeftovers_ZH->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceLeftovers_ZH->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceLeftovers_ZH->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceLeftovers_ZH->GetYaxis()->SetTitleOffset(1.);
  AcceptanceLeftovers_ZH->GetXaxis()->SetTitleOffset(.9);
  AcceptanceLeftovers_ZH->GetZaxis()->SetTitleOffset(.75);
  AcceptanceLeftovers_ZH->GetXaxis()->SetLabelSize(0.05);
  AcceptanceLeftovers_ZH->GetYaxis()->SetLabelSize(0.05);
  AcceptanceLeftovers_ZH->GetZaxis()->SetLabelSize(0.04);

  TH3F* ggWe_ZH  = (TH3F*)finZH->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT");
  TH3F* ggWe_ZH_noSF  = (TH3F*)finZH->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor");
  TH3F* ggWe_ZHJECup  = (TH3F*)finZH->Get("SMS_ggMT_JECup_Loose_1Ele_0_1Jets");
  TH3F* ggWe_ZHJECdown  = (TH3F*)finZH->Get("SMS_ggMT_JECdown_Loose_1Ele_0_1Jets");
  TH3F* ggWmu_ZH = (TH3F*)finZH->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT");
  TH3F* ggWmu_ZH_noSF = (TH3F*)finZH->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT_noScaleFactor");
  TH3F* ggWmu_ZHJECup  = (TH3F*)finZH->Get("SMS_ggMT_JECup_Loose_1Mu_0_1Jets");
  TH3F* ggWmu_ZHJECdown  = (TH3F*)finZH->Get("SMS_ggMT_JECdown_Loose_1Mu_0_1Jets");
  TH3F* ggLeftovers_ZH = (TH3F*)finZH->Get("gg_SMS_Loose_lt2b_noMu_noEle_lt2jets_mChi_mBino_met");
  TH2F* h_nEvents_ZH = (TH2F*)finZH->Get("SMS_mChi_mBino"); 	
  ggWe_ZH->Sumw2();ggWe_ZH_noSF->Sumw2();ggWe_ZHJECup->Sumw2();ggWe_ZHJECdown->Sumw2();
  ggWmu_ZH->Sumw2();ggWmu_ZH_noSF->Sumw2();ggWmu_ZHJECup->Sumw2();ggWmu_ZHJECdown->Sumw2();
  ggLeftovers_ZH->Sumw2();h_nEvents_ZH->Sumw2();
  h_nEvents_ZH->GetXaxis()->SetTitleOffset(0.9);		     
  h_nEvents_ZH->GetYaxis()->SetTitleOffset(0.75);		     
  h_nEvents_ZH->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  h_nEvents_ZH->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  c2->cd();
  h_nEvents_ZH->Draw("COLZtext");
  TPaveText *PrelimText_ZH = new TPaveText(.39,.86,.89,.9,"NDC");
  PrelimText_ZH->AddText("CMS Preliminary                  SMS Z(#rightarrowInclusive)H(#gamma#gamma)+2Jet");
  PrelimText_ZH->SetFillStyle(0);
  PrelimText_ZH->SetFillColor(0);
  PrelimText_ZH->SetBorderSize(0);
  PrelimText_ZH->Draw();
  EventsText->Draw();
  c2->Print("Plots/Higgs/AcceptanceEvents_SMS_TChiZH_ZincHgg_2J.png");
  c2->Print("Plots/Higgs/AcceptanceEvents_SMS_TChiZH_ZincHgg_2J.pdf");
  c1->cd();
  //std::string filename;
  //if(inputfilesAAW.is_open()){
  for(int i=1;i<=ggWe_ZH->GetNbinsX();i++){
    for(int j=1;j<=ggWe_ZH->GetNbinsY();j++){

      //while(!inputfilesAAW.eof()){
      //std::getline(inputfilesAAW,filename);
      //cout<<"file: "<<filename<<endl;
      //TFile f(filename.c_str(),"READ");
      //f.cd();
      //TString str = f.GetName();
      //int one = str.Index("_chargino");
      //one+=9;
      //int two = str.Index("_",one+1);
      //TString MCHI (str(one,two-one));
    
      //int mChi = MCHI.Atof();
      
      int mChi = ggWe_ZH->GetXaxis()->GetBinLowEdge(i);
      int mBino= ggWe_ZH->GetYaxis()->GetBinLowEdge(j);
      TString MCHI;ostringstream conv,conv2;conv<<mChi;MCHI=conv.str();//= static_cast<ostringstream*>( &(ostringstream() << mChi) )->str();//sprintf(MCHI,"%i",mChi);
      TString MBINO;conv2<<mBino;MBINO=conv2.str();//sprintf(MBINO,"%i",mBino);
      TH1F* RHO = (TH1F*)finZH->Get("rho");RHO->Sumw2();
      float nEvents=RHO->GetEntries();
      if(mChi<150)nEvents=120000.;
      else if(mChi<=175)nEvents=60000.;
      else if(mChi<250)nEvents=16000.;
      else if(mChi>=250)nEvents=30000.;
      int chiBin=h_nEvents_ZH->GetXaxis()->FindBin(mChi);
      int binoBin=h_nEvents_ZH->GetYaxis()->FindBin(mBino);
      float nEventsHist = h_nEvents_ZH->Integral(chiBin,chiBin,binoBin,binoBin);
      if(nEventsHist==0 && ScanBin[i][j]==false)continue;
      //TH2F* ggWe_ZH = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Ele_0_1Jets");
      //if(ggWe_ZH){
      TH1F*  ggWe_ZHMt=(TH1F*)ggWe_ZH->ProjectionZ("ggWe_ZHMt",i,i,j,j,"o");
      TH1F*  ggWe_ZHMt_noSF=(TH1F*)ggWe_ZH_noSF->ProjectionZ("ggWe_ZHMt_noSF",i,i,j,j,"o");
      TH1F*  ggWe_ZHMtJECup=(TH1F*)ggWe_ZHJECup->ProjectionZ("ggWe_ZHMtJECup",i,i,j,j,"o");
      TH1F*  ggWe_ZHMtJECdown=(TH1F*)ggWe_ZHJECdown->ProjectionZ("ggWe_ZHMtJECdown",i,i,j,j,"o");
      float val=ggWe_ZHMt->GetEntries();
      //cout<<"mChi: "<<mChi<<"  mBino: "<<mBino<<"  nEvents: "<<nEvents<<"  gg We events: "<<val<<endl;
      val/=nEvents;
      //cout<<"highest_gge: "<<highest_gge<<" val: "<<val;
      if(val>0 && val<lowest_gge)lowest_gge=val;
      if(val>highest_gge)highest_gge=val;
      //cout<<"  highest_gge: "<<highest_gge<<endl;
      AcceptanceWe_ZH->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      float scale = ggWe_ZHMt_noSF->GetEntries()/ggWe_ZHMt_noSF->Integral(0,-1);
      //cout<<"scale: "<<scale<<endl;
      //float Norm = PhoEffScale2*x_secsHiggsino[mChi]*brFracHgg*L_int/nEventsHist;
      float Norm = PhoEffScale2*x_secsHiggsino[mChi]*bf_ZH_Zgg*L_int/nEventsHist;//higgsino model, ZH but Z inclusive so no Z->ee branching fraction
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      /*ggWe_ZHMt->Scale(scale);ggWe_ZHMt->Scale(Norm);
      ggWe_ZHMtJECup->Scale(scale);ggWe_ZHMtJECup->Scale(Norm);
      ggWe_ZHMtJECdown->Scale(scale);ggWe_ZHMtJECdown->Scale(Norm);*/
      HistScale(ggWe_ZHMt,scale*Norm);HistScale(ggWe_ZHMtJECup,scale*Norm);HistScale(ggWe_ZHMtJECdown,scale*Norm);
      TH1F* hist = (TH1F*)ggWe_ZHMt->Rebin(nMTbins,"hist",MTbins);
      AddOverflowToLastBin(hist);
      TFile fLimitsSigSMS_ZH_ggEle("HiggsLimitFiles/LimitPackage_mu_"+MCHI+"_mB_"+MBINO+"_ZH_ggEle.root","RECREATE");
      fLimitsSigSMS_ZH_ggEle.cd();
      //ggWe_ZHMt->Write("h_ggWe_ZH_MT_mChi"+MCHI+"_mBino"+MBINO);
      //ggWe_ZHMt->Write("hMT_ggEle_MC");
      hist->Write("hMT_ggEle_MC");
      //ggWe_ZHMtJECup->Write("hMT_ggEle_JECup");
      //ggWe_ZHMtJECdown->Write("hMT_ggEle_JECdown");
      //cout<<"Total expected events: "<<ggWe_ZHMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggWe_ZHMt->Integral(0,ggWe_ZHMt->FindBin(50))<<endl<<"mt>50 : "<<ggWe_ZHMt->Integral(ggWe_ZHMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggWe_ZHMt->Integral(ggWe_ZHMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggWe_ZHMt->Integral(ggWe_ZHMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggWe_ZHMt->Integral(ggWe_ZHMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggWe_ZHMt->Integral(ggWe_ZHMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggWe_ZHMt->Integral(ggWe_ZHMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggWe_ZHMt->Integral(ggWe_ZHMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggWe_ZHMt->FindBin(50),bin75=ggWe_ZHMt->FindBin(75),bin100=ggWe_ZHMt->FindBin(100),bin125=ggWe_ZHMt->FindBin(125),bin150=ggWe_ZHMt->FindBin(150),bin175=ggWe_ZHMt->FindBin(175),bin200=ggWe_ZHMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsWe_ZH.push_back(ggWe_ZHMt->Integral(0,-1));
      outValsWe_ZH.push_back(ggWe_ZHMt->Integral(0,bin50-1));
      outValsWe_ZH.push_back(ggWe_ZHMt->Integral(bin50,-1));
      outValsWe_ZH.push_back(ggWe_ZHMt->Integral(bin75,-1));
      outValsWe_ZH.push_back(ggWe_ZHMt->Integral(bin100,-1));
      //outValsWe_ZH.push_back(ggWe_ZHMt->Integral(bin125,-1));
      outValsWe_ZH.push_back(ggWe_ZHMt->Integral(bin150,-1));
      //outValsWe_ZH.push_back(ggWe_ZHMt->Integral(bin175,-1));
      outValsWe_ZH.push_back(ggWe_ZHMt->Integral(bin200,-1));
      
      //}
      //delete ggWe_ZH;
      fLimitsSigSMS_ZH_ggEle.Close();

      TFile fLimitsSigSMS_ZH_ggMu("HiggsLimitFiles/LimitPackage_mu_"+MCHI+"_mB_"+MBINO+"_ZH_ggMu.root","RECREATE");
      //TH2F* ggWmu_ZH = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
      //if(ggWmu_ZH){
      TH1F* ggWmu_ZHMt=(TH1F*)ggWmu_ZH->ProjectionZ("ggWmu_ZHMt",i,i,j,j,"o");
      TH1F* ggWmu_ZHMt_noSF=(TH1F*)ggWmu_ZH_noSF->ProjectionZ("ggWmu_ZHMt_noSF",i,i,j,j,"o");
      TH1F* ggWmu_ZHMtJECup=(TH1F*)ggWmu_ZHJECup->ProjectionZ("ggWmu_ZHMtJECup",i,i,j,j,"o");
      TH1F* ggWmu_ZHMtJECdown=(TH1F*)ggWmu_ZHJECdown->ProjectionZ("ggWmu_ZHMtJECdown",i,i,j,j,"o");
      float val=ggWmu_ZHMt->GetEntries();
      val/=nEvents;
      //cout<<"highest_ggmu: "<<highest_ggmu<<" val: "<<val;
      if(val>0 && val<lowest_ggmu)lowest_ggmu=val;
      if(val>highest_ggmu)highest_ggmu=val;
      //cout<<"  highest_ggmu: "<<highest_ggmu<<endl<<endl;
      AcceptanceWmu_ZH->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      float scale = ggWmu_ZHMt_noSF->GetEntries()/ggWmu_ZHMt_noSF->Integral(0,-1);
      //cout<<"scale: "<<scale<<endl;
      //float Norm = PhoEffScale2*x_secsHiggsino[mChi]*brFracHgg*L_int/nEventsHist;
      float Norm = PhoEffScale2*x_secsHiggsino[mChi]*bf_ZH_Zgg*L_int/nEventsHist;//higgsino model, ZH but Z inclusive so no Z->mumu branching fraction
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      /*ggWmu_ZHMt->Scale(scale);ggWmu_ZHMt->Scale(Norm);
      ggWmu_ZHMtJECup->Scale(scale);ggWmu_ZHMtJECup->Scale(Norm);
      ggWmu_ZHMtJECdown->Scale(scale);ggWmu_ZHMtJECdown->Scale(Norm);*/
      HistScale(ggWmu_ZHMt,scale*Norm);HistScale(ggWmu_ZHMtJECup,scale*Norm);HistScale(ggWmu_ZHMtJECdown,scale*Norm);
      TH1F* hist = (TH1F*)ggWmu_ZHMt->Rebin(nMTbins,"hist",MTbins);
      AddOverflowToLastBin(hist);
      fLimitsSigSMS_ZH_ggMu.cd();
      //ggWmu_ZHMt->Write("h_ggWmu_ZH_MT_mChi"+MCHI+"_mBino"+MBINO);
      //ggWmu_ZHMt->Write("hMT_ggMu_MC");
      hist->Write("hMT_ggMu_MC");
      //ggWmu_ZHMtJECup->Write("hMT_ggMu_JECup");
      //ggWmu_ZHMtJECdown->Write("hMT_ggMu_JECdown");
      //cout<<"Total expected events: "<<ggWmu_ZHMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggWmu_ZHMt->Integral(0,ggWmu_ZHMt->FindBin(50))<<endl<<"mt>50 : "<<ggWmu_ZHMt->Integral(ggWmu_ZHMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggWmu_ZHMt->Integral(ggWmu_ZHMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggWmu_ZHMt->Integral(ggWmu_ZHMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggWmu_ZHMt->Integral(ggWmu_ZHMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggWmu_ZHMt->Integral(ggWmu_ZHMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggWmu_ZHMt->Integral(ggWmu_ZHMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggWmu_ZHMt->Integral(ggWmu_ZHMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggWmu_ZHMt->FindBin(50),bin75=ggWmu_ZHMt->FindBin(75),bin100=ggWmu_ZHMt->FindBin(100),bin125=ggWmu_ZHMt->FindBin(125),bin150=ggWmu_ZHMt->FindBin(150),bin175=ggWmu_ZHMt->FindBin(175),bin200=ggWmu_ZHMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsWmu_ZH.push_back(ggWmu_ZHMt->Integral(0,-1));
      outValsWmu_ZH.push_back(ggWmu_ZHMt->Integral(0,bin50-1));
      outValsWmu_ZH.push_back(ggWmu_ZHMt->Integral(bin50,-1));
      outValsWmu_ZH.push_back(ggWmu_ZHMt->Integral(bin75,-1));
      outValsWmu_ZH.push_back(ggWmu_ZHMt->Integral(bin100,-1));
      //outValsWmu_ZH.push_back(ggWmu_ZHMt->Integral(bin125,-1));
      outValsWmu_ZH.push_back(ggWmu_ZHMt->Integral(bin150,-1));
      //outValsWmu_ZH.push_back(ggWmu_ZHMt->Integral(bin175,-1));
      outValsWmu_ZH.push_back(ggWmu_ZHMt->Integral(bin200,-1));
      // }
      

      fLimitsSigSMS_ZH_ggMu.Close(); 
      /*
      TFile fLimitsSigSMS_ZH_ggLeftovers("HiggsLimitFiles/LimitPackage_mu_"+MCHI+"_mB_"+MBINO+"_ZH_ggLeftovers.root","RECREATE");
      fLimitsSigSMS_ZH_ggLeftovers.cd();
      */
      //TH2F* ggLeftovers_ZH = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
      //if(ggLeftovers_ZH){
      TH1F* ggLeftovers_ZHMt=(TH1F*)ggLeftovers_ZH->ProjectionZ("ggLeftovers_ZHMt",i,i,j,j,"o");
      float val=ggLeftovers_ZHMt->GetEntries();
      val/=nEvents;
      //cout<<"highest_ggLeftovers: "<<highest_ggLeftovers<<" val: "<<val;
      if(val>0 && val<lowest_ggLeftovers)lowest_ggLeftovers=val;
      if(val>highest_ggLeftovers)highest_ggLeftovers=val;
      //cout<<"  highest_ggLeftovers: "<<highest_ggLeftovers<<endl<<endl;
      AcceptanceLeftovers_ZH->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      //fLimitsSigSMS_ZH.cd();
      //ggLeftovers_ZH->Write("h_ggLeftovers_ZH_MT_mChi"+MCHI+"_mBino"+MBINO);
      float scale = ggLeftovers_ZHMt->GetEntries()/ggLeftovers_ZHMt->Integral(0,-1);
      //cout<<"scale: "<<scale<<endl;
      float Norm = PhoEffScale2*x_secsHiggsino[mChi]*brFracHgg*L_int/nEventsHist;
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      /*ggLeftovers_ZHMt->Scale(scale);
	ggLeftovers_ZHMt->Scale(Norm);*/    
      HistScale(ggLeftovers_ZHMt,scale*Norm);
      //cout<<"Total expected events: "<<ggLeftovers_ZHMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggLeftovers_ZHMt->Integral(0,ggLeftovers_ZHMt->FindBin(50))<<endl<<"mt>50 : "<<ggLeftovers_ZHMt->Integral(ggLeftovers_ZHMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggLeftovers_ZHMt->Integral(ggLeftovers_ZHMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggLeftovers_ZHMt->Integral(ggLeftovers_ZHMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggLeftovers_ZHMt->Integral(ggLeftovers_ZHMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggLeftovers_ZHMt->Integral(ggLeftovers_ZHMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggLeftovers_ZHMt->Integral(ggLeftovers_ZHMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggLeftovers_ZHMt->Integral(ggLeftovers_ZHMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggLeftovers_ZHMt->FindBin(50),bin75=ggLeftovers_ZHMt->FindBin(75),bin100=ggLeftovers_ZHMt->FindBin(100),bin125=ggLeftovers_ZHMt->FindBin(125),bin150=ggLeftovers_ZHMt->FindBin(150),bin175=ggLeftovers_ZHMt->FindBin(175),bin200=ggLeftovers_ZHMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsLeftovers_ZH.push_back(ggLeftovers_ZHMt->Integral(0,-1));
      outValsLeftovers_ZH.push_back(ggLeftovers_ZHMt->Integral(0,bin50-1));
      outValsLeftovers_ZH.push_back(ggLeftovers_ZHMt->Integral(bin50,-1));
      outValsLeftovers_ZH.push_back(ggLeftovers_ZHMt->Integral(bin75,-1));
      outValsLeftovers_ZH.push_back(ggLeftovers_ZHMt->Integral(bin100,-1));
      //outValsLeftovers_ZH.push_back(ggLeftovers_ZHMt->Integral(bin125,-1));
      outValsLeftovers_ZH.push_back(ggLeftovers_ZHMt->Integral(bin150,-1));
      //outValsLeftovers_ZH.push_back(ggLeftovers_ZHMt->Integral(bin175,-1));
      outValsLeftovers_ZH.push_back(ggLeftovers_ZHMt->Integral(bin200,-1));
      // }
      //fLimitsSigSMS_ZH_ggLeftovers.Close();
    }
  }

  //now do HH_2b2g

  //TFile finHH_2b2g("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_SMS_TChiHH_2b2g_2J.root","READ");
  TFile *finHH_2b2g=TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Photon_SMS_TChiHH_2b2g_2J.root","READ");
  //TFile fLimitsSigSMS_HH_2b2g("results_SMS_TChiHH_2b2g_2J.root","RECREATE");

  //lowest=999999.;highest=0.;

 
  vector<float> outValsWe_HH_2b2g,outValsWmu_HH_2b2g,outValsLeftovers_HH_2b2g;
  int nCats=7;

  TH2F* AcceptanceWe_HH_2b2g = new TH2F("AcceptanceWe_HH_2b2g","",nXbin,xbins,nYbin,ybins);  			     
  AcceptanceWe_HH_2b2g->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceWe_HH_2b2g->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceWe_HH_2b2g->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceWe_HH_2b2g->GetYaxis()->SetTitleOffset(1.);
  AcceptanceWe_HH_2b2g->GetXaxis()->SetTitleOffset(0.9);
  AcceptanceWe_HH_2b2g->GetZaxis()->SetTitleOffset(.75);
  AcceptanceWe_HH_2b2g->GetXaxis()->SetLabelSize(0.05);
  AcceptanceWe_HH_2b2g->GetYaxis()->SetLabelSize(0.05);
  AcceptanceWe_HH_2b2g->GetZaxis()->SetLabelSize(0.04);
  TH2F* AcceptanceWmu_HH_2b2g = new TH2F("AcceptanceWmu_HH_2b2g","",nXbin,xbins,nYbin,ybins);  			     
  AcceptanceWmu_HH_2b2g->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceWmu_HH_2b2g->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceWmu_HH_2b2g->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceWmu_HH_2b2g->GetYaxis()->SetTitleOffset(1.);
  AcceptanceWmu_HH_2b2g->GetXaxis()->SetTitleOffset(.9);
  AcceptanceWmu_HH_2b2g->GetZaxis()->SetTitleOffset(.75);
  AcceptanceWmu_HH_2b2g->GetXaxis()->SetLabelSize(0.05);
  AcceptanceWmu_HH_2b2g->GetYaxis()->SetLabelSize(0.05);
  AcceptanceWmu_HH_2b2g->GetZaxis()->SetLabelSize(0.04);			     
  TH2F* AcceptanceLeftovers_HH_2b2g = new TH2F("AcceptanceLeftovers_HH_2b2g","",nXbin,xbins,nYbin,ybins);  			     
  AcceptanceLeftovers_HH_2b2g->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceLeftovers_HH_2b2g->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceLeftovers_HH_2b2g->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceLeftovers_HH_2b2g->GetYaxis()->SetTitleOffset(1.);
  AcceptanceLeftovers_HH_2b2g->GetXaxis()->SetTitleOffset(.9);
  AcceptanceLeftovers_HH_2b2g->GetZaxis()->SetTitleOffset(.75);
  AcceptanceLeftovers_HH_2b2g->GetXaxis()->SetLabelSize(0.05);
  AcceptanceLeftovers_HH_2b2g->GetYaxis()->SetLabelSize(0.05);
  AcceptanceLeftovers_HH_2b2g->GetZaxis()->SetLabelSize(0.04);

  TH3F* ggWe_HH_2b2g  = (TH3F*)finHH_2b2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT");
  TH3F* ggWe_HH_2b2g_noSF  = (TH3F*)finHH_2b2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor");
  TH3F* ggWe_HH_2b2gJECup  = (TH3F*)finHH_2b2g->Get("SMS_ggMT_JECup_Loose_1Ele_0_1Jets");
  TH3F* ggWe_HH_2b2gJECdown  = (TH3F*)finHH_2b2g->Get("SMS_ggMT_JECdown_Loose_1Ele_0_1Jets");
  TH3F* ggWmu_HH_2b2g = (TH3F*)finHH_2b2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT");
  TH3F* ggWmu_HH_2b2g_noSF = (TH3F*)finHH_2b2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT_noScaleFactor");
  TH3F* ggWmu_HH_2b2gJECup  = (TH3F*)finHH_2b2g->Get("SMS_ggMT_JECup_Loose_1Mu_0_1Jets");
  TH3F* ggWmu_HH_2b2gJECdown  = (TH3F*)finHH_2b2g->Get("SMS_ggMT_JECdown_Loose_1Mu_0_1Jets");
  TH3F* ggLeftovers_HH_2b2g = (TH3F*)finHH_2b2g->Get("gg_SMS_Loose_lt2b_noMu_noEle_lt2jets_mChi_mBino_met");
  TH2F* h_nEvents_HH_2b2g = (TH2F*)finHH_2b2g->Get("SMS_mChi_mBino"); 	 
  ggWe_HH_2b2g->Sumw2();ggWe_HH_2b2g_noSF->Sumw2();ggWe_HH_2b2gJECup->Sumw2();ggWe_HH_2b2gJECdown->Sumw2();
  ggWmu_HH_2b2g->Sumw2();ggWmu_HH_2b2g_noSF->Sumw2();ggWmu_HH_2b2gJECup->Sumw2();ggWmu_HH_2b2gJECdown->Sumw2();
  ggLeftovers_HH_2b2g->Sumw2();h_nEvents_HH_2b2g->Sumw2();
  h_nEvents_HH_2b2g->GetXaxis()->SetTitleOffset(0.9);		     
  h_nEvents_HH_2b2g->GetYaxis()->SetTitleOffset(0.75);		     
  h_nEvents_HH_2b2g->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  h_nEvents_HH_2b2g->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  c2->cd();
  h_nEvents_HH_2b2g->Draw("COLZtext");
  TPaveText *PrelimText_HH_2b2g = new TPaveText(.39,.86,.89,.9,"NDC");
  PrelimText_HH_2b2g->AddText("CMS Preliminary                  SMS H(bb)H(#gamma#gamma)+2Jet");
  PrelimText_HH_2b2g->SetFillStyle(0);
  PrelimText_HH_2b2g->SetFillColor(0);
  PrelimText_HH_2b2g->SetBorderSize(0);
  PrelimText_HH_2b2g->Draw();
  EventsText->Draw();
  c2->Print("Plots/Higgs/AcceptanceEvents_SMS_TChiHH_2b2g_2J.png");
  c2->Print("Plots/Higgs/AcceptanceEvents_SMS_TChiHH_2b2g_2J.pdf");
  c1->cd();
  //std::string filename;
  //if(inputfilesAAW.is_open()){
  for(int i=1;i<=ggWe_HH_2b2g->GetNbinsX();i++){
    for(int j=1;j<=ggWe_HH_2b2g->GetNbinsY();j++){
      //while(!inputfilesAAW.eof()){
      //std::getline(inputfilesAAW,filename);
      //cout<<"file: "<<filename<<endl;
      //TFile f(filename.c_str(),"READ");
      //f.cd();
      //TString str = f.GetName();
      //int one = str.Index("_chargino");
      //one+=9;
      //int two = str.Index("_",one+1);
      //TString MCHI (str(one,two-one));
    
      //int mChi = MCHI.Atof();
      
      int mChi = ggWe_HH_2b2g->GetXaxis()->GetBinLowEdge(i);
      int mBino= ggWe_HH_2b2g->GetYaxis()->GetBinLowEdge(j);
      TString MCHI;ostringstream conv,conv2;conv<<mChi;MCHI=conv.str();//= static_cast<ostringstream*>( &(ostringstream() << mChi) )->str();//sprintf(MCHI,"%i",mChi);
      TString MBINO;conv2<<mBino;MBINO=conv2.str();//sprintf(MBINO,"%i",mBino);
      TH1F* RHO = (TH1F*)finHH_2b2g->Get("rho");RHO->Sumw2();
      float nEvents=RHO->GetEntries();
      if(mChi<150)nEvents=60000.;
      else nEvents=30000.;
      int chiBin=h_nEvents_HH_2b2g->GetXaxis()->FindBin(mChi);
      int binoBin=h_nEvents_HH_2b2g->GetYaxis()->FindBin(mBino);
      float nEventsHist = h_nEvents_HH_2b2g->Integral(chiBin,chiBin,binoBin,binoBin);
      if(nEventsHist==0 && ScanBin[i][j]==false){/*cout<<"skipping mChi: "<<mChi<<"  mBino: "<<mBino<<endl;*/continue;}
      //TH2F* ggWe_HH_2b2g = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Ele_0_1Jets");
      //if(ggWe_HH_2b2g){
      TH1F*  ggWe_HH_2b2gMt=(TH1F*)ggWe_HH_2b2g->ProjectionZ("ggWe_HH_2b2gMt",i,i,j,j,"o");
      TH1F*  ggWe_HH_2b2gMt_noSF=(TH1F*)ggWe_HH_2b2g_noSF->ProjectionZ("ggWe_HH_2b2gMt_noSF",i,i,j,j,"o");
      TH1F*  ggWe_HH_2b2gMtJECup=(TH1F*)ggWe_HH_2b2gJECup->ProjectionZ("ggWe_HH_2b2gMtJECup",i,i,j,j,"o");
      TH1F*  ggWe_HH_2b2gMtJECdown=(TH1F*)ggWe_HH_2b2gJECdown->ProjectionZ("ggWe_HH_2b2gMtJECdown",i,i,j,j,"o");
      float val=ggWe_HH_2b2gMt->GetEntries();
      float nEventsWe = ggWe_HH_2b2gMt->Integral(0,-1);
      //cout<<"mChi: "<<mChi<<"  mBino: "<<mBino<<"  nEvents: "<<nEventsHist<<"  gg We events: "<<val<<endl;
      val/=nEvents;
      //cout<<"highest_gge: "<<highest_gge<<" val: "<<val;
      if(val>0 && val<lowest_gge)lowest_gge=val;
      if(val>highest_gge)highest_gge=val;
      //cout<<"  highest_gge: "<<highest_gge<<endl;
      AcceptanceWe_HH_2b2g->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      float scale = 0.;
      if(nEventsWe==0)scale=0.;
      else scale = ggWe_HH_2b2gMt_noSF->GetEntries()/ggWe_HH_2b2gMt_noSF->Integral(0,-1);
      //cout<<"scale: "<<scale<<endl;
      //float Norm = PhoEffScale2*x_secsHiggsino[mChi]*brFracHgg*L_int/nEventsHist;
      float Norm = 0;
      if(nEventsWe==0)Norm=0.;
      else Norm = PhoEffScale2*x_secsHiggsino[mChi]*bf_HH_bbgg*L_int/nEventsHist;//higgsino model, HH_2b2g
      //cout<<"x_secs["<<mChi<<"]: "<<x_secsHiggsino[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      /*ggWe_HH_2b2gMt->Scale(scale);ggWe_HH_2b2gMt->Scale(Norm);
      ggWe_HH_2b2gMtJECup->Scale(scale);ggWe_HH_2b2gMtJECup->Scale(Norm);
      ggWe_HH_2b2gMtJECdown->Scale(scale);ggWe_HH_2b2gMtJECdown->Scale(Norm);*/
      HistScale(ggWe_HH_2b2gMt,scale*Norm);HistScale(ggWe_HH_2b2gMtJECup,scale*Norm);HistScale(ggWe_HH_2b2gMtJECdown,scale*Norm);
      TH1F* hist = (TH1F*)ggWe_HH_2b2gMt->Rebin(nMTbins,"hist",MTbins);
      AddOverflowToLastBin(hist);
      if(nEventsWe==0){
	for(int bin=0;bin<=hist->GetNbinsX()+1;bin++){
	  hist->SetBinContent(bin,0);hist->SetBinError(bin,0);
	}
      }
      TFile fLimitsSigSMS_HH_2b2g_ggEle("HiggsLimitFiles/LimitPackage_mu_"+MCHI+"_mB_"+MBINO+"_HH_2b2g_ggEle.root","RECREATE");
      fLimitsSigSMS_HH_2b2g_ggEle.cd();
      //ggWe_HH_2b2gMt->Write("h_ggWe_HH_2b2g_MT_mChi"+MCHI+"_mBino"+MBINO);
      //ggWe_HH_2b2gMt->Write("hMT_ggEle_MC");
      hist->Write("hMT_ggEle_MC");
      //ggWe_HH_2b2gMtJECup->Write("hMT_ggEle_JECUp");
      //ggWe_HH_2b2gMtJECdown->Write("hMT_ggEle_JECDown");
      //cout<<"Total expected events: "<<ggWe_HH_2b2gMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggWe_HH_2b2gMt->Integral(0,ggWe_HH_2b2gMt->FindBin(50))<<endl<<"mt>50 : "<<ggWe_HH_2b2gMt->Integral(ggWe_HH_2b2gMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggWe_HH_2b2gMt->Integral(ggWe_HH_2b2gMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggWe_HH_2b2gMt->Integral(ggWe_HH_2b2gMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggWe_HH_2b2gMt->Integral(ggWe_HH_2b2gMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggWe_HH_2b2gMt->Integral(ggWe_HH_2b2gMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggWe_HH_2b2gMt->Integral(ggWe_HH_2b2gMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggWe_HH_2b2gMt->Integral(ggWe_HH_2b2gMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggWe_HH_2b2gMt->FindBin(50),bin75=ggWe_HH_2b2gMt->FindBin(75),bin100=ggWe_HH_2b2gMt->FindBin(100),bin125=ggWe_HH_2b2gMt->FindBin(125),bin150=ggWe_HH_2b2gMt->FindBin(150),bin175=ggWe_HH_2b2gMt->FindBin(175),bin200=ggWe_HH_2b2gMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsWe_HH_2b2g.push_back(ggWe_HH_2b2gMt->Integral(0,-1));
      outValsWe_HH_2b2g.push_back(ggWe_HH_2b2gMt->Integral(0,bin50-1));
      outValsWe_HH_2b2g.push_back(ggWe_HH_2b2gMt->Integral(bin50,-1));
      outValsWe_HH_2b2g.push_back(ggWe_HH_2b2gMt->Integral(bin75,-1));
      outValsWe_HH_2b2g.push_back(ggWe_HH_2b2gMt->Integral(bin100,-1));
      //outValsWe_HH_2b2g.push_back(ggWe_HH_2b2gMt->Integral(bin125,-1));
      outValsWe_HH_2b2g.push_back(ggWe_HH_2b2gMt->Integral(bin150,-1));
      //outValsWe_HH_2b2g.push_back(ggWe_HH_2b2gMt->Integral(bin175,-1));
      outValsWe_HH_2b2g.push_back(ggWe_HH_2b2gMt->Integral(bin200,-1));
      
      //}
      //delete ggWe_HH_2b2g;
      fLimitsSigSMS_HH_2b2g_ggEle.Close();

      //TH2F* ggWmu_HH_2b2g = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
      //if(ggWmu_HH_2b2g){
      TH1F* ggWmu_HH_2b2gMt=(TH1F*)ggWmu_HH_2b2g->ProjectionZ("ggWmu_HH_2b2gMt",i,i,j,j,"o");
      TH1F* ggWmu_HH_2b2gMt_noSF=(TH1F*)ggWmu_HH_2b2g_noSF->ProjectionZ("ggWmu_HH_2b2gMt_noSF",i,i,j,j,"o");
      TH1F* ggWmu_HH_2b2gMtJECup=(TH1F*)ggWmu_HH_2b2gJECup->ProjectionZ("ggWmu_HH_2b2gMtJECup",i,i,j,j,"o");
      TH1F* ggWmu_HH_2b2gMtJECdown=(TH1F*)ggWmu_HH_2b2gJECdown->ProjectionZ("ggWmu_HH_2b2gMtJECdown",i,i,j,j,"o");
      float val=ggWmu_HH_2b2gMt->GetEntries();
      float nEventsWmu = ggWmu_HH_2b2gMt->Integral(0,-1);
      val/=nEvents;
      //cout<<"highest_ggmu: "<<highest_ggmu<<" val: "<<val;
      if(val>0 && val<lowest_ggmu)lowest_ggmu=val;
      if(val>highest_ggmu)highest_ggmu=val;
      //cout<<"  highest_ggmu: "<<highest_ggmu<<endl<<endl;
      AcceptanceWmu_HH_2b2g->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      float scale = 0;
      if(nEventsWmu==0)scale=0.;
      else scale = ggWmu_HH_2b2gMt_noSF->GetEntries()/ggWmu_HH_2b2gMt_noSF->Integral(0,-1);
      //cout<<"scale: "<<scale<<endl;
      //float Norm = PhoEffScale2*x_secsHiggsino[mChi]*brFracHgg*L_int/nEventsHist;
      float Norm = 0.;
      if(nEventsWmu==0)Norm=0.;
      else Norm = PhoEffScale2*x_secsHiggsino[mChi]*bf_HH_bbgg*L_int/nEventsHist;//higgsino model, HH_2b2g
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      /*ggWmu_HH_2b2gMt->Scale(scale);ggWmu_HH_2b2gMt->Scale(Norm);
      ggWmu_HH_2b2gMtJECup->Scale(scale);ggWmu_HH_2b2gMtJECup->Scale(Norm);
      ggWmu_HH_2b2gMtJECdown->Scale(scale);ggWmu_HH_2b2gMtJECdown->Scale(Norm);*/
      HistScale(ggWmu_HH_2b2gMt,scale*Norm);HistScale(ggWmu_HH_2b2gMtJECup,scale*Norm);HistScale(ggWmu_HH_2b2gMtJECdown,scale*Norm);
      TH1F* hist = (TH1F*)ggWmu_HH_2b2gMt->Rebin(nMTbins,"hist",MTbins);
      AddOverflowToLastBin(hist);
      if(nEventsHist==0){
	for(int binM=0;binM<=hist->GetNbinsX()+1;binM++){
	  hist->SetBinContent(binM,0);hist->SetBinError(binM,0);
	}
      }
      TFile fLimitsSigSMS_HH_2b2g_ggMu("HiggsLimitFiles/LimitPackage_mu_"+MCHI+"_mB_"+MBINO+"_HH_2b2g_ggMu.root","RECREATE");
      fLimitsSigSMS_HH_2b2g_ggMu.cd();
      //ggWmu_HH_2b2gMt->Write("h_ggWmu_HH_2b2g_MT_mChi"+MCHI+"_mBino"+MBINO);
      //ggWmu_HH_2b2gMt->Write("hMT_ggMu_MC");
      hist->Write("hMT_ggMu_MC");
      //ggWmu_HH_2b2gMtJECup->Write("hMT_ggMu_JECUp");
      //ggWmu_HH_2b2gMtJECdown->Write("hMT_ggMu_JECDown");
      //cout<<"Total expected events: "<<ggWmu_HH_2b2gMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggWmu_HH_2b2gMt->Integral(0,ggWmu_HH_2b2gMt->FindBin(50))<<endl<<"mt>50 : "<<ggWmu_HH_2b2gMt->Integral(ggWmu_HH_2b2gMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggWmu_HH_2b2gMt->Integral(ggWmu_HH_2b2gMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggWmu_HH_2b2gMt->Integral(ggWmu_HH_2b2gMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggWmu_HH_2b2gMt->Integral(ggWmu_HH_2b2gMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggWmu_HH_2b2gMt->Integral(ggWmu_HH_2b2gMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggWmu_HH_2b2gMt->Integral(ggWmu_HH_2b2gMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggWmu_HH_2b2gMt->Integral(ggWmu_HH_2b2gMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggWmu_HH_2b2gMt->FindBin(50),bin75=ggWmu_HH_2b2gMt->FindBin(75),bin100=ggWmu_HH_2b2gMt->FindBin(100),bin125=ggWmu_HH_2b2gMt->FindBin(125),bin150=ggWmu_HH_2b2gMt->FindBin(150),bin175=ggWmu_HH_2b2gMt->FindBin(175),bin200=ggWmu_HH_2b2gMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsWmu_HH_2b2g.push_back(ggWmu_HH_2b2gMt->Integral(0,-1));
      outValsWmu_HH_2b2g.push_back(ggWmu_HH_2b2gMt->Integral(0,bin50-1));
      outValsWmu_HH_2b2g.push_back(ggWmu_HH_2b2gMt->Integral(bin50,-1));
      outValsWmu_HH_2b2g.push_back(ggWmu_HH_2b2gMt->Integral(bin75,-1));
      outValsWmu_HH_2b2g.push_back(ggWmu_HH_2b2gMt->Integral(bin100,-1));
      //outValsWmu_HH_2b2g.push_back(ggWmu_HH_2b2gMt->Integral(bin125,-1));
      outValsWmu_HH_2b2g.push_back(ggWmu_HH_2b2gMt->Integral(bin150,-1));
      //outValsWmu_HH_2b2g.push_back(ggWmu_HH_2b2gMt->Integral(bin175,-1));
      outValsWmu_HH_2b2g.push_back(ggWmu_HH_2b2gMt->Integral(bin200,-1));
      // }
      fLimitsSigSMS_HH_2b2g_ggMu.Close();

      //TH2F* ggLeftovers_HH_2b2g = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
      //if(ggLeftovers_HH_2b2g){
      TH1F* ggLeftovers_HH_2b2gMt=(TH1F*)ggLeftovers_HH_2b2g->ProjectionZ("ggLeftovers_HH_2b2gMt",i,i,j,j,"o");
      float val=ggLeftovers_HH_2b2gMt->GetEntries();
      val/=nEvents;
      //cout<<"highest_ggLeftovers: "<<highest_ggLeftovers<<" val: "<<val;
      if(val>0 && val<lowest_ggLeftovers)lowest_ggLeftovers=val;
      if(val>highest_ggLeftovers)highest_ggLeftovers=val;
      //cout<<"  highest_ggLeftovers: "<<highest_ggLeftovers<<endl<<endl;
      AcceptanceLeftovers_HH_2b2g->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      //fLimitsSigSMS_HH_2b2g.cd();
      //ggLeftovers_HH_2b2g->Write("h_ggLeftovers_HH_2b2g_MT_mChi"+MCHI+"_mBino"+MBINO);
      float scale = ggLeftovers_HH_2b2gMt->GetEntries()/ggLeftovers_HH_2b2gMt->Integral(0,-1);
      //cout<<"scale: "<<scale<<endl;
      float Norm = PhoEffScale2*x_secsHiggsino[mChi]*brFracHgg*L_int/nEventsHist;
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      /*ggLeftovers_HH_2b2gMt->Scale(scale);
	ggLeftovers_HH_2b2gMt->Scale(Norm);*/
      HistScale(ggLeftovers_HH_2b2gMt,scale*Norm);
      //cout<<"Total expected events: "<<ggLeftovers_HH_2b2gMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggLeftovers_HH_2b2gMt->Integral(0,ggLeftovers_HH_2b2gMt->FindBin(50))<<endl<<"mt>50 : "<<ggLeftovers_HH_2b2gMt->Integral(ggLeftovers_HH_2b2gMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggLeftovers_HH_2b2gMt->Integral(ggLeftovers_HH_2b2gMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggLeftovers_HH_2b2gMt->Integral(ggLeftovers_HH_2b2gMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggLeftovers_HH_2b2gMt->Integral(ggLeftovers_HH_2b2gMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggLeftovers_HH_2b2gMt->Integral(ggLeftovers_HH_2b2gMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggLeftovers_HH_2b2gMt->Integral(ggLeftovers_HH_2b2gMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggLeftovers_HH_2b2gMt->Integral(ggLeftovers_HH_2b2gMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggLeftovers_HH_2b2gMt->FindBin(50),bin75=ggLeftovers_HH_2b2gMt->FindBin(75),bin100=ggLeftovers_HH_2b2gMt->FindBin(100),bin125=ggLeftovers_HH_2b2gMt->FindBin(125),bin150=ggLeftovers_HH_2b2gMt->FindBin(150),bin175=ggLeftovers_HH_2b2gMt->FindBin(175),bin200=ggLeftovers_HH_2b2gMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsLeftovers_HH_2b2g.push_back(ggLeftovers_HH_2b2gMt->Integral(0,-1));
      outValsLeftovers_HH_2b2g.push_back(ggLeftovers_HH_2b2gMt->Integral(0,bin50-1));
      outValsLeftovers_HH_2b2g.push_back(ggLeftovers_HH_2b2gMt->Integral(bin50,-1));
      outValsLeftovers_HH_2b2g.push_back(ggLeftovers_HH_2b2gMt->Integral(bin75,-1));
      outValsLeftovers_HH_2b2g.push_back(ggLeftovers_HH_2b2gMt->Integral(bin100,-1));
      //outValsLeftovers_HH_2b2g.push_back(ggLeftovers_HH_2b2gMt->Integral(bin125,-1));
      outValsLeftovers_HH_2b2g.push_back(ggLeftovers_HH_2b2gMt->Integral(bin150,-1));
      //outValsLeftovers_HH_2b2g.push_back(ggLeftovers_HH_2b2gMt->Integral(bin175,-1));
      outValsLeftovers_HH_2b2g.push_back(ggLeftovers_HH_2b2gMt->Integral(bin200,-1));
      // }
    }
  }

  //now do HH_2Z2g

  //TFile finHH_2Z2g("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_SMS_TChiHH_2Z2g_2J.root","READ");
  TFile *finHH_2Z2g=TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Photon_SMS_TChiHH_2Z2g_2J.root","READ");
  //TFile fLimitsSigSMS_HH_2Z2g("results_SMS_TChiHH_2Z2g_2J.root","RECREATE");

  //lowest=999999.;highest=0.;

 
  vector<float> outValsWe_HH_2Z2g,outValsWmu_HH_2Z2g,outValsLeftovers_HH_2Z2g;
  int nCats=7;

  TH2F* AcceptanceWe_HH_2Z2g = new TH2F("AcceptanceWe_HH_2Z2g","",nXbin,xbins,nYbin,ybins);  			     
  AcceptanceWe_HH_2Z2g->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceWe_HH_2Z2g->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceWe_HH_2Z2g->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceWe_HH_2Z2g->GetYaxis()->SetTitleOffset(1.);
  AcceptanceWe_HH_2Z2g->GetXaxis()->SetTitleOffset(0.9);
  AcceptanceWe_HH_2Z2g->GetZaxis()->SetTitleOffset(.75);
  AcceptanceWe_HH_2Z2g->GetXaxis()->SetLabelSize(0.05);
  AcceptanceWe_HH_2Z2g->GetYaxis()->SetLabelSize(0.05);
  AcceptanceWe_HH_2Z2g->GetZaxis()->SetLabelSize(0.04);
  TH2F* AcceptanceWmu_HH_2Z2g = new TH2F("AcceptanceWmu_HH_2Z2g","",nXbin,xbins,nYbin,ybins);  			     
  AcceptanceWmu_HH_2Z2g->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceWmu_HH_2Z2g->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceWmu_HH_2Z2g->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceWmu_HH_2Z2g->GetYaxis()->SetTitleOffset(1.);
  AcceptanceWmu_HH_2Z2g->GetXaxis()->SetTitleOffset(.9);
  AcceptanceWmu_HH_2Z2g->GetZaxis()->SetTitleOffset(.75);
  AcceptanceWmu_HH_2Z2g->GetXaxis()->SetLabelSize(0.05);
  AcceptanceWmu_HH_2Z2g->GetYaxis()->SetLabelSize(0.05);
  AcceptanceWmu_HH_2Z2g->GetZaxis()->SetLabelSize(0.04);			     
  TH2F* AcceptanceLeftovers_HH_2Z2g = new TH2F("AcceptanceLeftovers_HH_2Z2g","",nXbin,xbins,nYbin,ybins);  			     
  AcceptanceLeftovers_HH_2Z2g->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceLeftovers_HH_2Z2g->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceLeftovers_HH_2Z2g->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceLeftovers_HH_2Z2g->GetYaxis()->SetTitleOffset(1.);
  AcceptanceLeftovers_HH_2Z2g->GetXaxis()->SetTitleOffset(.9);
  AcceptanceLeftovers_HH_2Z2g->GetZaxis()->SetTitleOffset(.75);
  AcceptanceLeftovers_HH_2Z2g->GetXaxis()->SetLabelSize(0.05);
  AcceptanceLeftovers_HH_2Z2g->GetYaxis()->SetLabelSize(0.05);
  AcceptanceLeftovers_HH_2Z2g->GetZaxis()->SetLabelSize(0.04);

  TH3F* ggWe_HH_2Z2g  = (TH3F*)finHH_2Z2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT");
  TH3F* ggWe_HH_2Z2g_noSF  = (TH3F*)finHH_2Z2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor");
  TH3F* ggWe_HH_2Z2gJECup  = (TH3F*)finHH_2Z2g->Get("SMS_ggMT_JECup_Loose_1Ele_0_1Jets");
  TH3F* ggWe_HH_2Z2gJECdown  = (TH3F*)finHH_2Z2g->Get("SMS_ggMT_JECdown_Loose_1Ele_0_1Jets");
  TH3F* ggWmu_HH_2Z2g = (TH3F*)finHH_2Z2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT");
  TH3F* ggWmu_HH_2Z2g_noSF = (TH3F*)finHH_2Z2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT_noScaleFactor");
  TH3F* ggWmu_HH_2Z2gJECup  = (TH3F*)finHH_2Z2g->Get("SMS_ggMT_JECup_Loose_1Mu_0_1Jets");
  TH3F* ggWmu_HH_2Z2gJECdown  = (TH3F*)finHH_2Z2g->Get("SMS_ggMT_JECdown_Loose_1Mu_0_1Jets");
  TH3F* ggLeftovers_HH_2Z2g = (TH3F*)finHH_2Z2g->Get("gg_SMS_Loose_lt2b_noMu_noEle_lt2jets_mChi_mBino_met");
  TH2F* h_nEvents_HH_2Z2g = (TH2F*)finHH_2Z2g->Get("SMS_mChi_mBino"); 	 
  ggWe_HH_2Z2g->Sumw2();ggWe_HH_2Z2g_noSF->Sumw2();ggWe_HH_2Z2gJECup->Sumw2();ggWe_HH_2Z2gJECdown->Sumw2();
  ggWmu_HH_2Z2g->Sumw2();ggWmu_HH_2Z2g_noSF->Sumw2();ggWmu_HH_2Z2gJECup->Sumw2();ggWmu_HH_2Z2gJECdown->Sumw2();
  ggLeftovers_HH_2Z2g->Sumw2();h_nEvents_HH_2Z2g->Sumw2();
  h_nEvents_HH_2Z2g->GetXaxis()->SetTitleOffset(0.9);		     
  h_nEvents_HH_2Z2g->GetYaxis()->SetTitleOffset(0.75);		     
  h_nEvents_HH_2Z2g->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  h_nEvents_HH_2Z2g->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  c2->cd();
  h_nEvents_HH_2Z2g->Draw("COLZtext");
  TPaveText *PrelimText_HH_2Z2g = new TPaveText(.39,.86,.89,.9,"NDC");
  PrelimText_HH_2Z2g->AddText("CMS Preliminary                  SMS H(ZZ)H(#gamma#gamma)+2Jet");
  PrelimText_HH_2Z2g->SetFillStyle(0);
  PrelimText_HH_2Z2g->SetFillColor(0);
  PrelimText_HH_2Z2g->SetBorderSize(0);
  PrelimText_HH_2Z2g->Draw();
  EventsText->Draw();
  c2->Print("Plots/Higgs/AcceptanceEvents_SMS_TChiHH_2Z2g_2J.png");
  c2->Print("Plots/Higgs/AcceptanceEvents_SMS_TChiHH_2Z2g_2J.pdf");
  c1->cd();
  //std::string filename;
  //if(inputfilesAAW.is_open()){
  for(int i=1;i<=ggWe_HH_2Z2g->GetNbinsX();i++){
    for(int j=1;j<=ggWe_HH_2Z2g->GetNbinsY();j++){
      //while(!inputfilesAAW.eof()){
      //std::getline(inputfilesAAW,filename);
      //cout<<"file: "<<filename<<endl;
      //TFile f(filename.c_str(),"READ");
      //f.cd();
      //TString str = f.GetName();
      //int one = str.Index("_chargino");
      //one+=9;
      //int two = str.Index("_",one+1);
      //TString MCHI (str(one,two-one));
    
      //int mChi = MCHI.Atof();
      
      int mChi = ggWe_HH_2Z2g->GetXaxis()->GetBinLowEdge(i);
      int mBino= ggWe_HH_2Z2g->GetYaxis()->GetBinLowEdge(j);
      TString MCHI;ostringstream conv,conv2;conv<<mChi;MCHI=conv.str();//= static_cast<ostringstream*>( &(ostringstream() << mChi) )->str();//sprintf(MCHI,"%i",mChi);
      TString MBINO;conv2<<mBino;MBINO=conv2.str();//sprintf(MBINO,"%i",mBino);
      TH1F* RHO = (TH1F*)finHH_2Z2g->Get("rho");RHO->Sumw2();
      float nEvents=RHO->GetEntries();
      if(mChi<150)nEvents=60000.;
      else nEvents=30000.;
      int chiBin=h_nEvents_HH_2Z2g->GetXaxis()->FindBin(mChi);
      int binoBin=h_nEvents_HH_2Z2g->GetYaxis()->FindBin(mBino);
      float nEventsHist = h_nEvents_HH_2Z2g->Integral(chiBin,chiBin,binoBin,binoBin);
      if(nEventsHist==0 && ScanBin[i][j]==false)continue;
      //TH2F* ggWe_HH_2Z2g = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Ele_0_1Jets");
      //if(ggWe_HH_2Z2g){
      TH1F*  ggWe_HH_2Z2gMt=(TH1F*)ggWe_HH_2Z2g->ProjectionZ("ggWe_HH_2Z2gMt",i,i,j,j,"o");
      TH1F*  ggWe_HH_2Z2gMt_noSF=(TH1F*)ggWe_HH_2Z2g_noSF->ProjectionZ("ggWe_HH_2Z2gMt_noSF",i,i,j,j,"o");
      TH1F*  ggWe_HH_2Z2gMtJECup=(TH1F*)ggWe_HH_2Z2gJECup->ProjectionZ("ggWe_HH_2Z2gMtJECup",i,i,j,j,"o");
      TH1F*  ggWe_HH_2Z2gMtJECdown=(TH1F*)ggWe_HH_2Z2gJECdown->ProjectionZ("ggWe_HH_2Z2gMtJECdown",i,i,j,j,"o");
      float val=ggWe_HH_2Z2gMt->GetEntries();
      //cout<<"mChi: "<<mChi<<"  mBino: "<<mBino<<"  nEvents: "<<nEvents<<"  gg We events: "<<val<<endl;
      val/=nEvents;
      //cout<<"highest_gge: "<<highest_gge<<" val: "<<val;
      if(val>0 && val<lowest_gge)lowest_gge=val;
      if(val>highest_gge)highest_gge=val;
      //cout<<"  highest_gge: "<<highest_gge<<endl;
      AcceptanceWe_HH_2Z2g->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      float scale = ggWe_HH_2Z2gMt_noSF->GetEntries()/ggWe_HH_2Z2gMt_noSF->Integral(0,-1);
      //cout<<"scale: "<<scale<<endl;
      //float Norm = PhoEffScale2*x_secsHiggsino[mChi]*brFracHgg*L_int/nEventsHist;
      float Norm = PhoEffScale2*x_secsHiggsino[mChi]*bf_HH_ZZgg*L_int/nEventsHist;//higgsino model, HH_2Z2g
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      /*ggWe_HH_2Z2gMt->Scale(scale);ggWe_HH_2Z2gMt->Scale(Norm);
      ggWe_HH_2Z2gMtJECup->Scale(scale);ggWe_HH_2Z2gMtJECdown->Scale(Norm);
      ggWe_HH_2Z2gMtJECup->Scale(scale);ggWe_HH_2Z2gMtJECdown->Scale(Norm);*/
      HistScale(ggWe_HH_2Z2gMt,scale*Norm);HistScale(ggWe_HH_2Z2gMtJECup,scale*Norm);HistScale(ggWe_HH_2Z2gMtJECdown,scale*Norm);
      TH1F* hist = (TH1F*)ggWe_HH_2Z2gMt->Rebin(nMTbins,"hist",MTbins);
      AddOverflowToLastBin(hist);
      TFile fLimitsSigSMS_HH_2Z2g_ggEle("HiggsLimitFiles/LimitPackage_mu_"+MCHI+"_mB_"+MBINO+"_HH_2Z2g_ggEle.root","RECREATE");
      fLimitsSigSMS_HH_2Z2g_ggEle.cd();
      //ggWe_HH_2Z2gMt->Write("h_ggWe_HH_2Z2g_MT_mChi"+MCHI+"_mBino"+MBINO);
      //ggWe_HH_2Z2gMt->Write("hMT_ggEle_MC");
      hist->Write("hMT_ggEle_MC");
      //ggWe_HH_2Z2gMtJECup->Write("hMT_ggEle_JECup");
      //ggWe_HH_2Z2gMtJECdown->Write("hMT_ggEle_JECdown");
      //cout<<"Total expected events: "<<ggWe_HH_2Z2gMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggWe_HH_2Z2gMt->Integral(0,ggWe_HH_2Z2gMt->FindBin(50))<<endl<<"mt>50 : "<<ggWe_HH_2Z2gMt->Integral(ggWe_HH_2Z2gMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggWe_HH_2Z2gMt->Integral(ggWe_HH_2Z2gMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggWe_HH_2Z2gMt->Integral(ggWe_HH_2Z2gMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggWe_HH_2Z2gMt->Integral(ggWe_HH_2Z2gMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggWe_HH_2Z2gMt->Integral(ggWe_HH_2Z2gMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggWe_HH_2Z2gMt->Integral(ggWe_HH_2Z2gMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggWe_HH_2Z2gMt->Integral(ggWe_HH_2Z2gMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggWe_HH_2Z2gMt->FindBin(50),bin75=ggWe_HH_2Z2gMt->FindBin(75),bin100=ggWe_HH_2Z2gMt->FindBin(100),bin125=ggWe_HH_2Z2gMt->FindBin(125),bin150=ggWe_HH_2Z2gMt->FindBin(150),bin175=ggWe_HH_2Z2gMt->FindBin(175),bin200=ggWe_HH_2Z2gMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsWe_HH_2Z2g.push_back(ggWe_HH_2Z2gMt->Integral(0,-1));
      outValsWe_HH_2Z2g.push_back(ggWe_HH_2Z2gMt->Integral(0,bin50-1));
      outValsWe_HH_2Z2g.push_back(ggWe_HH_2Z2gMt->Integral(bin50,-1));
      outValsWe_HH_2Z2g.push_back(ggWe_HH_2Z2gMt->Integral(bin75,-1));
      outValsWe_HH_2Z2g.push_back(ggWe_HH_2Z2gMt->Integral(bin100,-1));
      //outValsWe_HH_2Z2g.push_back(ggWe_HH_2Z2gMt->Integral(bin125,-1));
      outValsWe_HH_2Z2g.push_back(ggWe_HH_2Z2gMt->Integral(bin150,-1));
      //outValsWe_HH_2Z2g.push_back(ggWe_HH_2Z2gMt->Integral(bin175,-1));
      outValsWe_HH_2Z2g.push_back(ggWe_HH_2Z2gMt->Integral(bin200,-1));
      
      //}
      //delete ggWe_HH_2Z2g;
      fLimitsSigSMS_HH_2Z2g_ggEle.Close();


      //TH2F* ggWmu_HH_2Z2g = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
      //if(ggWmu_HH_2Z2g){
      TH1F* ggWmu_HH_2Z2gMt=(TH1F*)ggWmu_HH_2Z2g->ProjectionZ("ggWmu_HH_2Z2gMt",i,i,j,j,"o");
      TH1F* ggWmu_HH_2Z2gMt_noSF=(TH1F*)ggWmu_HH_2Z2g_noSF->ProjectionZ("ggWmu_HH_2Z2gMt_noSF",i,i,j,j,"o");
      TH1F* ggWmu_HH_2Z2gMtJECup=(TH1F*)ggWmu_HH_2Z2gJECup->ProjectionZ("ggWmu_HH_2Z2gMtJECup",i,i,j,j,"o");
      TH1F* ggWmu_HH_2Z2gMtJECdown=(TH1F*)ggWmu_HH_2Z2gJECdown->ProjectionZ("ggWmu_HH_2Z2gMtJECdown",i,i,j,j,"o");
      float val=ggWmu_HH_2Z2gMt->GetEntries();
      val/=nEvents;
      //cout<<"highest_ggmu: "<<highest_ggmu<<" val: "<<val;
      if(val>0 && val<lowest_ggmu)lowest_ggmu=val;
      if(val>highest_ggmu)highest_ggmu=val;
      //cout<<"  highest_ggmu: "<<highest_ggmu<<endl<<endl;
      AcceptanceWmu_HH_2Z2g->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      float scale = ggWmu_HH_2Z2gMt_noSF->GetEntries()/ggWmu_HH_2Z2gMt_noSF->Integral(0,-1);
      //cout<<"scale: "<<scale<<endl;
      //float Norm = PhoEffScale2*x_secsHiggsino[mChi]*brFracHgg*L_int/nEventsHist;
      float Norm = PhoEffScale2*x_secsHiggsino[mChi]*bf_HH_ZZgg*L_int/nEventsHist;//higgsino model, HH_2Z2g
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      /*ggWmu_HH_2Z2gMt->Scale(scale);ggWmu_HH_2Z2gMt->Scale(Norm);
      ggWmu_HH_2Z2gMtJECup->Scale(scale);ggWmu_HH_2Z2gMtJECup->Scale(Norm);
      ggWmu_HH_2Z2gMtJECdown->Scale(scale);ggWmu_HH_2Z2gMtJECdown->Scale(Norm);*/
      HistScale(ggWmu_HH_2Z2gMt,scale*Norm);HistScale(ggWmu_HH_2Z2gMtJECup,scale*Norm);HistScale(ggWmu_HH_2Z2gMtJECdown,scale*Norm);
      TH1F* hist = (TH1F*)ggWmu_HH_2Z2gMt->Rebin(nMTbins,"hist",MTbins);
      AddOverflowToLastBin(hist);
      TFile fLimitsSigSMS_HH_2Z2g_ggMu("HiggsLimitFiles/LimitPackage_mu_"+MCHI+"_mB_"+MBINO+"_HH_2Z2g_ggMu.root","RECREATE");
      fLimitsSigSMS_HH_2Z2g_ggMu.cd();
      //ggWmu_HH_2Z2gMt->Write("h_ggWmu_HH_2Z2g_MT_mChi"+MCHI+"_mBino"+MBINO);
      //ggWmu_HH_2Z2gMt->Write("hMT_ggMu_MC");
      hist->Write("hMT_ggMu_MC");
      //ggWmu_HH_2Z2gMtJECup->Write("hMT_ggMu_JECUp");
      //ggWmu_HH_2Z2gMtJECdown->Write("hMT_ggMu_JECDown");
      //cout<<"Total expected events: "<<ggWmu_HH_2Z2gMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggWmu_HH_2Z2gMt->Integral(0,ggWmu_HH_2Z2gMt->FindBin(50))<<endl<<"mt>50 : "<<ggWmu_HH_2Z2gMt->Integral(ggWmu_HH_2Z2gMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggWmu_HH_2Z2gMt->Integral(ggWmu_HH_2Z2gMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggWmu_HH_2Z2gMt->Integral(ggWmu_HH_2Z2gMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggWmu_HH_2Z2gMt->Integral(ggWmu_HH_2Z2gMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggWmu_HH_2Z2gMt->Integral(ggWmu_HH_2Z2gMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggWmu_HH_2Z2gMt->Integral(ggWmu_HH_2Z2gMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggWmu_HH_2Z2gMt->Integral(ggWmu_HH_2Z2gMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggWmu_HH_2Z2gMt->FindBin(50),bin75=ggWmu_HH_2Z2gMt->FindBin(75),bin100=ggWmu_HH_2Z2gMt->FindBin(100),bin125=ggWmu_HH_2Z2gMt->FindBin(125),bin150=ggWmu_HH_2Z2gMt->FindBin(150),bin175=ggWmu_HH_2Z2gMt->FindBin(175),bin200=ggWmu_HH_2Z2gMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsWmu_HH_2Z2g.push_back(ggWmu_HH_2Z2gMt->Integral(0,-1));
      outValsWmu_HH_2Z2g.push_back(ggWmu_HH_2Z2gMt->Integral(0,bin50-1));
      outValsWmu_HH_2Z2g.push_back(ggWmu_HH_2Z2gMt->Integral(bin50,-1));
      outValsWmu_HH_2Z2g.push_back(ggWmu_HH_2Z2gMt->Integral(bin75,-1));
      outValsWmu_HH_2Z2g.push_back(ggWmu_HH_2Z2gMt->Integral(bin100,-1));
      //outValsWmu_HH_2Z2g.push_back(ggWmu_HH_2Z2gMt->Integral(bin125,-1));
      outValsWmu_HH_2Z2g.push_back(ggWmu_HH_2Z2gMt->Integral(bin150,-1));
      //outValsWmu_HH_2Z2g.push_back(ggWmu_HH_2Z2gMt->Integral(bin175,-1));
      outValsWmu_HH_2Z2g.push_back(ggWmu_HH_2Z2gMt->Integral(bin200,-1));
      // }
      fLimitsSigSMS_HH_2Z2g_ggMu.Close();

      //TH2F* ggLeftovers_HH_2Z2g = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
      //if(ggLeftovers_HH_2Z2g){
      TH1F* ggLeftovers_HH_2Z2gMt=(TH1F*)ggLeftovers_HH_2Z2g->ProjectionZ("ggLeftovers_HH_2Z2gMt",i,i,j,j,"o");
      float val=ggLeftovers_HH_2Z2gMt->GetEntries();
      val/=nEvents;
      //cout<<"highest_ggLeftovers: "<<highest_ggLeftovers<<" val: "<<val;
      if(val>0 && val<lowest_ggLeftovers)lowest_ggLeftovers=val;
      if(val>highest_ggLeftovers)highest_ggLeftovers=val;
      //cout<<"  highest_ggLeftovers: "<<highest_ggLeftovers<<endl<<endl;
      AcceptanceLeftovers_HH_2Z2g->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      float scale = ggLeftovers_HH_2Z2gMt->GetEntries()/ggLeftovers_HH_2Z2gMt->Integral(0,-1);
      //cout<<"scale: "<<scale<<endl;
      float Norm = PhoEffScale2*x_secsHiggsino[mChi]*brFracHgg*L_int/nEventsHist;
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      /*ggLeftovers_HH_2Z2gMt->Scale(scale);
	ggLeftovers_HH_2Z2gMt->Scale(Norm);*/
      HistScale(ggLeftovers_HH_2Z2gMt,scale*Norm);
      //fLimitsSigSMS_HH_2Z2g.cd();
      //ggLeftovers_HH_2Z2g->Write("h_ggLeftovers_HH_2Z2g_MT_mChi"+MCHI+"_mBino"+MBINO);
      //cout<<"Total expected events: "<<ggLeftovers_HH_2Z2gMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggLeftovers_HH_2Z2gMt->Integral(0,ggLeftovers_HH_2Z2gMt->FindBin(50))<<endl<<"mt>50 : "<<ggLeftovers_HH_2Z2gMt->Integral(ggLeftovers_HH_2Z2gMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggLeftovers_HH_2Z2gMt->Integral(ggLeftovers_HH_2Z2gMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggLeftovers_HH_2Z2gMt->Integral(ggLeftovers_HH_2Z2gMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggLeftovers_HH_2Z2gMt->Integral(ggLeftovers_HH_2Z2gMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggLeftovers_HH_2Z2gMt->Integral(ggLeftovers_HH_2Z2gMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggLeftovers_HH_2Z2gMt->Integral(ggLeftovers_HH_2Z2gMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggLeftovers_HH_2Z2gMt->Integral(ggLeftovers_HH_2Z2gMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggLeftovers_HH_2Z2gMt->FindBin(50),bin75=ggLeftovers_HH_2Z2gMt->FindBin(75),bin100=ggLeftovers_HH_2Z2gMt->FindBin(100),bin125=ggLeftovers_HH_2Z2gMt->FindBin(125),bin150=ggLeftovers_HH_2Z2gMt->FindBin(150),bin175=ggLeftovers_HH_2Z2gMt->FindBin(175),bin200=ggLeftovers_HH_2Z2gMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsLeftovers_HH_2Z2g.push_back(ggLeftovers_HH_2Z2gMt->Integral(0,-1));
      outValsLeftovers_HH_2Z2g.push_back(ggLeftovers_HH_2Z2gMt->Integral(0,bin50-1));
      outValsLeftovers_HH_2Z2g.push_back(ggLeftovers_HH_2Z2gMt->Integral(bin50,-1));
      outValsLeftovers_HH_2Z2g.push_back(ggLeftovers_HH_2Z2gMt->Integral(bin75,-1));
      outValsLeftovers_HH_2Z2g.push_back(ggLeftovers_HH_2Z2gMt->Integral(bin100,-1));
      //outValsLeftovers_HH_2Z2g.push_back(ggLeftovers_HH_2Z2gMt->Integral(bin125,-1));
      outValsLeftovers_HH_2Z2g.push_back(ggLeftovers_HH_2Z2gMt->Integral(bin150,-1));
      //outValsLeftovers_HH_2Z2g.push_back(ggLeftovers_HH_2Z2gMt->Integral(bin175,-1));
      outValsLeftovers_HH_2Z2g.push_back(ggLeftovers_HH_2Z2gMt->Integral(bin200,-1));
      // }
    }
  }
  
  //now do HH_2W2g

  //TFile finHH_2W2g("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_SMS_TChiHH_2W2g_2J.root","READ");
  TFile *finHH_2W2g=TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Photon_SMS_TChiHH_2W2g_2J.root","READ");
  TFile fLimitsSigSMS_HH_2W2g("results_SMS_TChiHH_2W2g_2J.root","RECREATE");

  //lowest=999999.;highest=0.;

 
  vector<float> outValsWe_HH_2W2g,outValsWmu_HH_2W2g,outValsLeftovers_HH_2W2g;
  int nCats=7;

  TH2F* AcceptanceWe_HH_2W2g = new TH2F("AcceptanceWe_HH_2W2g","",nXbin,xbins,nYbin,ybins);  			     
  AcceptanceWe_HH_2W2g->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceWe_HH_2W2g->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceWe_HH_2W2g->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceWe_HH_2W2g->GetYaxis()->SetTitleOffset(1.);
  AcceptanceWe_HH_2W2g->GetXaxis()->SetTitleOffset(0.9);
  AcceptanceWe_HH_2W2g->GetZaxis()->SetTitleOffset(.75);
  AcceptanceWe_HH_2W2g->GetXaxis()->SetLabelSize(0.05);
  AcceptanceWe_HH_2W2g->GetYaxis()->SetLabelSize(0.05);
  AcceptanceWe_HH_2W2g->GetZaxis()->SetLabelSize(0.04);
  TH2F* AcceptanceWmu_HH_2W2g = new TH2F("AcceptanceWmu_HH_2W2g","",nXbin,xbins,nYbin,ybins);  			     
  AcceptanceWmu_HH_2W2g->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceWmu_HH_2W2g->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceWmu_HH_2W2g->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceWmu_HH_2W2g->GetYaxis()->SetTitleOffset(1.);
  AcceptanceWmu_HH_2W2g->GetXaxis()->SetTitleOffset(.9);
  AcceptanceWmu_HH_2W2g->GetZaxis()->SetTitleOffset(.75);
  AcceptanceWmu_HH_2W2g->GetXaxis()->SetLabelSize(0.05);
  AcceptanceWmu_HH_2W2g->GetYaxis()->SetLabelSize(0.05);
  AcceptanceWmu_HH_2W2g->GetZaxis()->SetLabelSize(0.04);			     
  TH2F* AcceptanceLeftovers_HH_2W2g = new TH2F("AcceptanceLeftovers_HH_2W2g","",nXbin,xbins,nYbin,ybins);  			     
  AcceptanceLeftovers_HH_2W2g->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceLeftovers_HH_2W2g->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceLeftovers_HH_2W2g->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceLeftovers_HH_2W2g->GetYaxis()->SetTitleOffset(1.);
  AcceptanceLeftovers_HH_2W2g->GetXaxis()->SetTitleOffset(.9);
  AcceptanceLeftovers_HH_2W2g->GetZaxis()->SetTitleOffset(.75);
  AcceptanceLeftovers_HH_2W2g->GetXaxis()->SetLabelSize(0.05);
  AcceptanceLeftovers_HH_2W2g->GetYaxis()->SetLabelSize(0.05);
  AcceptanceLeftovers_HH_2W2g->GetZaxis()->SetLabelSize(0.04);

  TH3F* ggWe_HH_2W2g  = (TH3F*)finHH_2W2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT");
  TH3F* ggWe_HH_2W2g_noSF  = (TH3F*)finHH_2W2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor");
  TH3F* ggWe_HH_2W2gJECup  = (TH3F*)finHH_2W2g->Get("SMS_ggMT_JECup_Loose_1Ele_0_1Jets");
  TH3F* ggWe_HH_2W2gJECdown  = (TH3F*)finHH_2W2g->Get("SMS_ggMT_JECdown_Loose_1Ele_0_1Jets");
  TH3F* ggWmu_HH_2W2g = (TH3F*)finHH_2W2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT");
  TH3F* ggWmu_HH_2W2g_noSF = (TH3F*)finHH_2W2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT_noScaleFactor");
  TH3F* ggWmu_HH_2W2gJECup  = (TH3F*)finHH_2W2g->Get("SMS_ggMT_JECup_Loose_1Mu_0_1Jets");
  TH3F* ggWmu_HH_2W2gJECdown  = (TH3F*)finHH_2W2g->Get("SMS_ggMT_JECdown_Loose_1Mu_0_1Jets");
  TH3F* ggLeftovers_HH_2W2g = (TH3F*)finHH_2W2g->Get("gg_SMS_Loose_lt2b_noMu_noEle_lt2jets_mChi_mBino_met");
  TH2F* h_nEvents_HH_2W2g = (TH2F*)finHH_2W2g->Get("SMS_mChi_mBino"); 	 
  ggWe_HH_2W2g->Sumw2();ggWe_HH_2W2g_noSF->Sumw2();ggWe_HH_2W2gJECup->Sumw2();ggWe_HH_2W2gJECdown->Sumw2();
  ggWmu_HH_2W2g->Sumw2();ggWmu_HH_2W2g_noSF->Sumw2();ggWmu_HH_2W2gJECup->Sumw2();ggWmu_HH_2W2gJECdown->Sumw2();
  ggLeftovers_HH_2W2g->Sumw2();h_nEvents_HH_2W2g->Sumw2();
  h_nEvents_HH_2W2g->GetXaxis()->SetTitleOffset(0.9);		     
  h_nEvents_HH_2W2g->GetYaxis()->SetTitleOffset(0.75);		     
  h_nEvents_HH_2W2g->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  h_nEvents_HH_2W2g->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  c2->cd();
  h_nEvents_HH_2W2g->Draw("COLZtext");
  TPaveText *PrelimText_HH_2W2g = new TPaveText(.39,.86,.89,.9,"NDC");
  PrelimText_HH_2W2g->AddText("CMS Preliminary                  SMS H(WW)H(#gamma#gamma)+2Jet");
  PrelimText_HH_2W2g->SetFillStyle(0);
  PrelimText_HH_2W2g->SetFillColor(0);
  PrelimText_HH_2W2g->SetBorderSize(0);
  PrelimText_HH_2W2g->Draw();
  EventsText->Draw();
  c2->Print("Plots/Higgs/AcceptanceEvents_SMS_TChiHH_2W2g_2J.png");
  c2->Print("Plots/Higgs/AcceptanceEvents_SMS_TChiHH_2W2g_2J.pdf");
  c1->cd();
  //std::string filename;
  //if(inputfilesAAW.is_open()){
  for(int i=1;i<=ggWe_HH_2W2g->GetNbinsX();i++){
    for(int j=1;j<=ggWe_HH_2W2g->GetNbinsY();j++){
      //while(!inputfilesAAW.eof()){
      //std::getline(inputfilesAAW,filename);
      //cout<<"file: "<<filename<<endl;
      //TFile f(filename.c_str(),"READ");
      //f.cd();
      //TString str = f.GetName();
      //int one = str.Index("_chargino");
      //one+=9;
      //int two = str.Index("_",one+1);
      //TString MCHI (str(one,two-one));
    
      //int mChi = MCHI.Atof();
      
      int mChi = ggWe_HH_2W2g->GetXaxis()->GetBinLowEdge(i);
      int mBino= ggWe_HH_2W2g->GetYaxis()->GetBinLowEdge(j);
      TString MCHI;ostringstream conv,conv2;conv<<mChi;MCHI=conv.str();//= static_cast<ostringstream*>( &(ostringstream() << mChi) )->str();//sprintf(MCHI,"%i",mChi);
      TString MBINO;conv2<<mBino;MBINO=conv2.str();//sprintf(MBINO,"%i",mBino);
      TH1F* RHO = (TH1F*)finHH_2W2g->Get("rho");RHO->Sumw2();
      float nEvents=RHO->GetEntries();
      if(mChi<150)nEvents=60000.;
      else nEvents=30000.;
      int chiBin=h_nEvents_HH_2W2g->GetXaxis()->FindBin(mChi);
      int binoBin=h_nEvents_HH_2W2g->GetYaxis()->FindBin(mBino);
      float nEventsHist = h_nEvents_HH_2W2g->Integral(chiBin,chiBin,binoBin,binoBin);
      if(nEventsHist==0 && ScanBin[i][j]==false)continue;
      //TH2F* ggWe_HH_2W2g = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Ele_0_1Jets");
      //if(ggWe_HH_2W2g){
      TH1F*  ggWe_HH_2W2gMt=(TH1F*)ggWe_HH_2W2g->ProjectionZ("ggWe_HH_2W2gMt",i,i,j,j,"o");
      TH1F*  ggWe_HH_2W2gMt_noSF=(TH1F*)ggWe_HH_2W2g_noSF->ProjectionZ("ggWe_HH_2W2gMt_noSF",i,i,j,j,"o");
      TH1F*  ggWe_HH_2W2gMtJECup=(TH1F*)ggWe_HH_2W2gJECup->ProjectionZ("ggWe_HH_2W2gMtJECup",i,i,j,j,"o");
      TH1F*  ggWe_HH_2W2gMtJECdown=(TH1F*)ggWe_HH_2W2gJECdown->ProjectionZ("ggWe_HH_2W2gMtJECdown",i,i,j,j,"o");
      float val=ggWe_HH_2W2gMt->GetEntries();
      //cout<<"mChi: "<<mChi<<"  mBino: "<<mBino<<"  nEvents: "<<nEventsHist<<"  gg+e events: "<<val<<endl;
      val/=nEvents;
      //cout<<"highest_gge: "<<highest_gge<<" val: "<<val;
      if(val>0 && val<lowest_gge)lowest_gge=val;
      if(val>highest_gge)highest_gge=val;
      //cout<<"  highest_gge: "<<highest_gge<<endl;
      AcceptanceWe_HH_2W2g->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      //cout<<"raw counts:         "<<ggWe_HH_2W2gMt->Integral(0,-1)<<endl;
      float scale = ggWe_HH_2W2gMt_noSF->GetEntries()/ggWe_HH_2W2gMt_noSF->Integral(0,-1);
      //cout<<"scale:              "<<scale               <<"  counts after scale: "<<ggWe_HH_2W2gMt->Integral(0,-1)*scale<<endl;
      //float Norm = PhoEffScale2*x_secsHiggsino[mChi]*brFracHgg*L_int/nEventsHist;
      float Norm = PhoEffScale2*x_secsHiggsino[mChi]*bf_HH_WWgg*L_int/nEventsHist;//higgsino model, HH_2W2g
      //cout<<"PhoEffScale2:       "<<PhoEffScale2        <<"  counts after PhoEffScale2: "<<ggWe_HH_2W2gMt->Integral(0,-1)*scale*PhoEffScale2<<endl;
      //cout<<"x_secs["<<mChi<<"]: "<<x_secsHiggsino[mChi]<<"  counts after x_sec: "<<ggWe_HH_2W2gMt->Integral(0,-1)*scale*PhoEffScale2*x_secsHiggsino[mChi]<<endl;
      //cout<<"bf_HH_WWgg:         "<<bf_HH_WWgg          <<"  counts after bf_HH_WWgg: "<<ggWe_HH_2W2gMt->Integral(0,-1)*scale*PhoEffScale2*x_secsHiggsino[mChi]*bf_HH_WWgg<<endl;
      //cout<<"L_int:              "<<L_int               <<"  counts after L_int: "<<ggWe_HH_2W2gMt->Integral(0,-1)*scale*PhoEffScale2*x_secsHiggsino[mChi]*bf_HH_WWgg*L_int<<endl;
      //cout<<"nEventsHist:        "<<nEventsHist         <<"  counts after nEventsHist: "<<ggWe_HH_2W2gMt->Integral(0,-1)*scale*PhoEffScale2*x_secsHiggsino[mChi]*bf_HH_WWgg*L_int/nEventsHist<<endl;
      //cout<<"Norm:               "<<Norm                <<"  counts after Norm: "<<ggWe_HH_2W2gMt->Integral(0,-1)*scale*Norm<<endl;
      /*ggWe_HH_2W2gMt->Scale(scale);ggWe_HH_2W2gMt->Scale(Norm);
      ggWe_HH_2W2gMtJECup->Scale(scale);ggWe_HH_2W2gMtJECup->Scale(Norm);
      ggWe_HH_2W2gMtJECdown->Scale(scale);ggWe_HH_2W2gMtJECdown->Scale(Norm);*/
      HistScale(ggWe_HH_2W2gMt,scale*Norm);HistScale(ggWe_HH_2W2gMtJECup,scale*Norm);HistScale(ggWe_HH_2W2gMtJECdown,scale*Norm);
      TH1F* hist = (TH1F*)ggWe_HH_2W2gMt->Rebin(nMTbins,"hist",MTbins);
      AddOverflowToLastBin(hist);
      TFile fLimitsSigSMS_HH_2W2g_ggEle("HiggsLimitFiles/LimitPackage_mu_"+MCHI+"_mB_"+MBINO+"_HH_2W2g_ggEle.root","RECREATE");
      fLimitsSigSMS_HH_2W2g_ggEle.cd();
      //hist->Write("h_ggEle_HH_2W2g_MT_mChi"+MCHI+"_mBino"+MBINO);
      //ggWe_HH_2W2gMt->Write("hMT_ggEle_MC");
      //ggWe_HH_2W2gMtJECup->Write("hMT_ggEle_JECup");
      //ggWe_HH_2W2gMtJECdown->Write("hMT_ggEle_JECdown");
      hist->Write("hMT_ggEle_MC");
      //cout<<"Total expected events: "<<ggWe_HH_2W2gMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggWe_HH_2W2gMt->Integral(0,ggWe_HH_2W2gMt->FindBin(50))<<endl<<"mt>50 : "<<ggWe_HH_2W2gMt->Integral(ggWe_HH_2W2gMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggWe_HH_2W2gMt->Integral(ggWe_HH_2W2gMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggWe_HH_2W2gMt->Integral(ggWe_HH_2W2gMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggWe_HH_2W2gMt->Integral(ggWe_HH_2W2gMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggWe_HH_2W2gMt->Integral(ggWe_HH_2W2gMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggWe_HH_2W2gMt->Integral(ggWe_HH_2W2gMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggWe_HH_2W2gMt->Integral(ggWe_HH_2W2gMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggWe_HH_2W2gMt->FindBin(50),bin75=ggWe_HH_2W2gMt->FindBin(75),bin100=ggWe_HH_2W2gMt->FindBin(100),bin125=ggWe_HH_2W2gMt->FindBin(125),bin150=ggWe_HH_2W2gMt->FindBin(150),bin175=ggWe_HH_2W2gMt->FindBin(175),bin200=ggWe_HH_2W2gMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsWe_HH_2W2g.push_back(ggWe_HH_2W2gMt->Integral(0,-1));
      outValsWe_HH_2W2g.push_back(ggWe_HH_2W2gMt->Integral(0,bin50-1));
      outValsWe_HH_2W2g.push_back(ggWe_HH_2W2gMt->Integral(bin50,-1));
      outValsWe_HH_2W2g.push_back(ggWe_HH_2W2gMt->Integral(bin75,-1));
      outValsWe_HH_2W2g.push_back(ggWe_HH_2W2gMt->Integral(bin100,-1));
      //outValsWe_HH_2W2g.push_back(ggWe_HH_2W2gMt->Integral(bin125,-1));
      outValsWe_HH_2W2g.push_back(ggWe_HH_2W2gMt->Integral(bin150,-1));
      //outValsWe_HH_2W2g.push_back(ggWe_HH_2W2gMt->Integral(bin175,-1));
      outValsWe_HH_2W2g.push_back(ggWe_HH_2W2gMt->Integral(bin200,-1));
      
      //}
      //delete ggWe_HH_2W2g;
      fLimitsSigSMS_HH_2W2g_ggEle.Close();

      //TH2F* ggWmu_HH_2W2g = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
      //if(ggWmu_HH_2W2g){
      TH1F* ggWmu_HH_2W2gMt=(TH1F*)ggWmu_HH_2W2g->ProjectionZ("ggWmu_HH_2W2gMt",i,i,j,j,"o");
      TH1F* ggWmu_HH_2W2gMt_noSF=(TH1F*)ggWmu_HH_2W2g_noSF->ProjectionZ("ggWmu_HH_2W2gMt_noSF",i,i,j,j,"o");
      TH1F* ggWmu_HH_2W2gMtJECup=(TH1F*)ggWmu_HH_2W2gJECup->ProjectionZ("ggWmu_HH_2W2gMtJECup",i,i,j,j,"o");
      TH1F* ggWmu_HH_2W2gMtJECdown=(TH1F*)ggWmu_HH_2W2gJECdown->ProjectionZ("ggWmu_HH_2W2gMtJECdown",i,i,j,j,"o");
      float val=ggWmu_HH_2W2gMt->GetEntries();
      cout<<"mChi: "<<mChi<<"  mBino: "<<mBino<<"  nEvents: "<<nEventsHist<<"  gg+mu events: "<<val<<endl;
      val/=nEvents;
      //cout<<"highest_ggmu: "<<highest_ggmu<<" val: "<<val;
      if(val>0 && val<lowest_ggmu)lowest_ggmu=val;
      if(val>highest_ggmu)highest_ggmu=val;
      //cout<<"  highest_ggmu: "<<highest_ggmu<<endl<<endl;
      AcceptanceWmu_HH_2W2g->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      cout<<"raw counts:         "<<ggWmu_HH_2W2gMt->Integral(0,-1)<<endl;
      float scale = ggWmu_HH_2W2gMt_noSF->GetEntries()/ggWmu_HH_2W2gMt_noSF->Integral(0,-1);
      cout<<"scale:              "<<scale               <<"  counts after scale: "<<ggWmu_HH_2W2gMt->Integral(0,-1)*scale<<endl;
      //cout<<"scale: "<<scale<<endl;
      //float Norm = PhoEffScale2*x_secsHiggsino[mChi]*brFracHgg*L_int/nEventsHist;
      float Norm = PhoEffScale2*x_secsHiggsino[mChi]*bf_HH_WWgg*L_int/nEventsHist;//higgsino model, HH_2W2g
      cout<<"PhoEffScale2:       "<<PhoEffScale2        <<"  counts after PhoEffScale2: "<<ggWmu_HH_2W2gMt->Integral(0,-1)*scale*PhoEffScale2<<endl;
      cout<<"x_secs["<<mChi<<"]: "<<x_secsHiggsino[mChi]<<"  counts after x_sec: "<<ggWmu_HH_2W2gMt->Integral(0,-1)*scale*PhoEffScale2*x_secsHiggsino[mChi]<<endl;
      cout<<"bf_HH_WWgg:         "<<bf_HH_WWgg          <<"  counts after bf_HH_WWgg: "<<ggWmu_HH_2W2gMt->Integral(0,-1)*scale*PhoEffScale2*x_secsHiggsino[mChi]*bf_HH_WWgg<<endl;
      cout<<"L_int:              "<<L_int               <<"  counts after L_int: "<<ggWmu_HH_2W2gMt->Integral(0,-1)*scale*PhoEffScale2*x_secsHiggsino[mChi]*bf_HH_WWgg*L_int<<endl;
      cout<<"nEventsHist:        "<<nEventsHist         <<"  counts after nEventsHist: "<<ggWmu_HH_2W2gMt->Integral(0,-1)*scale*PhoEffScale2*x_secsHiggsino[mChi]*bf_HH_WWgg*L_int/nEventsHist<<endl;
      cout<<"Norm:               "<<Norm                <<"  counts after Norm: "<<ggWmu_HH_2W2gMt->Integral(0,-1)*scale*Norm<<endl;
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      /*ggWmu_HH_2W2gMt->Scale(scale);ggWmu_HH_2W2gMt->Scale(Norm);
      ggWmu_HH_2W2gMtJECup->Scale(scale);ggWmu_HH_2W2gMtJECup->Scale(Norm);
      ggWmu_HH_2W2gMtJECdown->Scale(scale);ggWmu_HH_2W2gMtJECdown->Scale(Norm);*/
      HistScale(ggWmu_HH_2W2gMt,scale*Norm);HistScale(ggWmu_HH_2W2gMtJECup,scale*Norm);HistScale(ggWmu_HH_2W2gMtJECdown,scale*Norm);
      TH1F* hist = (TH1F*)ggWmu_HH_2W2gMt->Rebin(nMTbins,"hist",MTbins);
      AddOverflowToLastBin(hist);
      TFile fLimitsSigSMS_HH_2W2g_ggMu("HiggsLimitFiles/LimitPackage_mu_"+MCHI+"_mB_"+MBINO+"_HH_2W2g_ggMu.root","RECREATE");
      fLimitsSigSMS_HH_2W2g_ggMu.cd();
      //hist->Write("h_ggMu_HH_2W2g_MT_mChi"+MCHI+"_mBino"+MBINO);
      //ggWmu_HH_2W2gMt->Write("hMT_ggMu_MC");
      //ggWmu_HH_2W2gMtJECup->Write("hMT_ggMu_JECup");
      //ggWmu_HH_2W2gMtJECdown->Write("hMT_ggMu_JECdown");
      hist->Write("hMT_ggMu_MC");
      //cout<<"Total expected events: "<<ggWmu_HH_2W2gMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggWmu_HH_2W2gMt->Integral(0,ggWmu_HH_2W2gMt->FindBin(50))<<endl<<"mt>50 : "<<ggWmu_HH_2W2gMt->Integral(ggWmu_HH_2W2gMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggWmu_HH_2W2gMt->Integral(ggWmu_HH_2W2gMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggWmu_HH_2W2gMt->Integral(ggWmu_HH_2W2gMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggWmu_HH_2W2gMt->Integral(ggWmu_HH_2W2gMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggWmu_HH_2W2gMt->Integral(ggWmu_HH_2W2gMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggWmu_HH_2W2gMt->Integral(ggWmu_HH_2W2gMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggWmu_HH_2W2gMt->Integral(ggWmu_HH_2W2gMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggWmu_HH_2W2gMt->FindBin(50),bin75=ggWmu_HH_2W2gMt->FindBin(75),bin100=ggWmu_HH_2W2gMt->FindBin(100),bin125=ggWmu_HH_2W2gMt->FindBin(125),bin150=ggWmu_HH_2W2gMt->FindBin(150),bin175=ggWmu_HH_2W2gMt->FindBin(175),bin200=ggWmu_HH_2W2gMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsWmu_HH_2W2g.push_back(ggWmu_HH_2W2gMt->Integral(0,-1));
      outValsWmu_HH_2W2g.push_back(ggWmu_HH_2W2gMt->Integral(0,bin50-1));
      outValsWmu_HH_2W2g.push_back(ggWmu_HH_2W2gMt->Integral(bin50,-1));
      outValsWmu_HH_2W2g.push_back(ggWmu_HH_2W2gMt->Integral(bin75,-1));
      outValsWmu_HH_2W2g.push_back(ggWmu_HH_2W2gMt->Integral(bin100,-1));
      //outValsWmu_HH_2W2g.push_back(ggWmu_HH_2W2gMt->Integral(bin125,-1));
      outValsWmu_HH_2W2g.push_back(ggWmu_HH_2W2gMt->Integral(bin150,-1));
      //outValsWmu_HH_2W2g.push_back(ggWmu_HH_2W2gMt->Integral(bin175,-1));
      outValsWmu_HH_2W2g.push_back(ggWmu_HH_2W2gMt->Integral(bin200,-1));
      // }
      fLimitsSigSMS_HH_2W2g_ggMu.Close();

      //TH2F* ggLeftovers_HH_2W2g = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
      //if(ggLeftovers_HH_2W2g){
      TH1F* ggLeftovers_HH_2W2gMt=(TH1F*)ggLeftovers_HH_2W2g->ProjectionZ("ggLeftovers_HH_2W2gMt",i,i,j,j,"o");
      float val=ggLeftovers_HH_2W2gMt->GetEntries();
      val/=nEvents;
      //cout<<"highest_ggLeftovers: "<<highest_ggLeftovers<<" val: "<<val;
      if(val>0 && val<lowest_ggLeftovers)lowest_ggLeftovers=val;
      if(val>highest_ggLeftovers)highest_ggLeftovers=val;
      //cout<<"  highest_ggLeftovers: "<<highest_ggLeftovers<<endl<<endl;
      AcceptanceLeftovers_HH_2W2g->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      //fLimitsSigSMS_HH_2W2g.cd();
      //ggLeftovers_HH_2W2g->Write("h_ggLeftovers_HH_2W2g_MT_mChi"+MCHI+"_mBino"+MBINO);
      float scale = ggLeftovers_HH_2W2gMt->GetEntries()/ggLeftovers_HH_2W2gMt->Integral(0,-1);
      //cout<<"scale: "<<scale<<endl;
      float Norm = PhoEffScale2*x_secsHiggsino[mChi]*brFracHgg*L_int/nEventsHist;
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      /*ggLeftovers_HH_2W2gMt->Scale(scale);
	ggLeftovers_HH_2W2gMt->Scale(Norm);*/
      HistScale(ggLeftovers_HH_2W2gMt,scale*Norm);
      //cout<<"Total expected events: "<<ggLeftovers_HH_2W2gMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggLeftovers_HH_2W2gMt->Integral(0,ggLeftovers_HH_2W2gMt->FindBin(50))<<endl<<"mt>50 : "<<ggLeftovers_HH_2W2gMt->Integral(ggLeftovers_HH_2W2gMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggLeftovers_HH_2W2gMt->Integral(ggLeftovers_HH_2W2gMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggLeftovers_HH_2W2gMt->Integral(ggLeftovers_HH_2W2gMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggLeftovers_HH_2W2gMt->Integral(ggLeftovers_HH_2W2gMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggLeftovers_HH_2W2gMt->Integral(ggLeftovers_HH_2W2gMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggLeftovers_HH_2W2gMt->Integral(ggLeftovers_HH_2W2gMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggLeftovers_HH_2W2gMt->Integral(ggLeftovers_HH_2W2gMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggLeftovers_HH_2W2gMt->FindBin(50),bin75=ggLeftovers_HH_2W2gMt->FindBin(75),bin100=ggLeftovers_HH_2W2gMt->FindBin(100),bin125=ggLeftovers_HH_2W2gMt->FindBin(125),bin150=ggLeftovers_HH_2W2gMt->FindBin(150),bin175=ggLeftovers_HH_2W2gMt->FindBin(175),bin200=ggLeftovers_HH_2W2gMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsLeftovers_HH_2W2g.push_back(ggLeftovers_HH_2W2gMt->Integral(0,-1));
      outValsLeftovers_HH_2W2g.push_back(ggLeftovers_HH_2W2gMt->Integral(0,bin50-1));
      outValsLeftovers_HH_2W2g.push_back(ggLeftovers_HH_2W2gMt->Integral(bin50,-1));
      outValsLeftovers_HH_2W2g.push_back(ggLeftovers_HH_2W2gMt->Integral(bin75,-1));
      outValsLeftovers_HH_2W2g.push_back(ggLeftovers_HH_2W2gMt->Integral(bin100,-1));
      //outValsLeftovers_HH_2W2g.push_back(ggLeftovers_HH_2W2gMt->Integral(bin125,-1));
      outValsLeftovers_HH_2W2g.push_back(ggLeftovers_HH_2W2gMt->Integral(bin150,-1));
      //outValsLeftovers_HH_2W2g.push_back(ggLeftovers_HH_2W2gMt->Integral(bin175,-1));
      outValsLeftovers_HH_2W2g.push_back(ggLeftovers_HH_2W2gMt->Integral(bin200,-1));
      // }
    }
  }
  
  //now do HH_2tau2g

  //TFile finHH_2tau2g("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_SMS_TChiHH_2tau2g_2J.root","READ");
  TFile *finHH_2tau2g=TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Photon_SMS_TChiHH_2tau2g_2J.root","READ");
  TFile fLimitsSigSMS_HH_2tau2g("results_SMS_TChiHH_2tau2g_2J.root","RECREATE");

  //lowest=999999.;highest=0.;

 
  vector<float> outValsWe_HH_2tau2g,outValsWmu_HH_2tau2g,outValsLeftovers_HH_2tau2g;
  int nCats=7;

  TH2F* AcceptanceWe_HH_2tau2g = new TH2F("AcceptanceWe_HH_2tau2g","",nXbin,xbins,nYbin,ybins);  			     
  AcceptanceWe_HH_2tau2g->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceWe_HH_2tau2g->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceWe_HH_2tau2g->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceWe_HH_2tau2g->GetYaxis()->SetTitleOffset(1.);
  AcceptanceWe_HH_2tau2g->GetXaxis()->SetTitleOffset(0.9);
  AcceptanceWe_HH_2tau2g->GetZaxis()->SetTitleOffset(.75);
  AcceptanceWe_HH_2tau2g->GetXaxis()->SetLabelSize(0.05);
  AcceptanceWe_HH_2tau2g->GetYaxis()->SetLabelSize(0.05);
  AcceptanceWe_HH_2tau2g->GetZaxis()->SetLabelSize(0.04);
  TH2F* AcceptanceWmu_HH_2tau2g = new TH2F("AcceptanceWmu_HH_2tau2g","",nXbin,xbins,nYbin,ybins);  			     
  AcceptanceWmu_HH_2tau2g->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceWmu_HH_2tau2g->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceWmu_HH_2tau2g->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceWmu_HH_2tau2g->GetYaxis()->SetTitleOffset(1.);
  AcceptanceWmu_HH_2tau2g->GetXaxis()->SetTitleOffset(.9);
  AcceptanceWmu_HH_2tau2g->GetZaxis()->SetTitleOffset(.75);
  AcceptanceWmu_HH_2tau2g->GetXaxis()->SetLabelSize(0.05);
  AcceptanceWmu_HH_2tau2g->GetYaxis()->SetLabelSize(0.05);
  AcceptanceWmu_HH_2tau2g->GetZaxis()->SetLabelSize(0.04);			     
  TH2F* AcceptanceLeftovers_HH_2tau2g = new TH2F("AcceptanceLeftovers_HH_2tau2g","",nXbin,xbins,nYbin,ybins);  			     
  AcceptanceLeftovers_HH_2tau2g->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceLeftovers_HH_2tau2g->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  AcceptanceLeftovers_HH_2tau2g->GetZaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceLeftovers_HH_2tau2g->GetYaxis()->SetTitleOffset(1.);
  AcceptanceLeftovers_HH_2tau2g->GetXaxis()->SetTitleOffset(.9);
  AcceptanceLeftovers_HH_2tau2g->GetZaxis()->SetTitleOffset(.75);
  AcceptanceLeftovers_HH_2tau2g->GetXaxis()->SetLabelSize(0.05);
  AcceptanceLeftovers_HH_2tau2g->GetYaxis()->SetLabelSize(0.05);
  AcceptanceLeftovers_HH_2tau2g->GetZaxis()->SetLabelSize(0.04);

  TH3F* ggWe_HH_2tau2g  = (TH3F*)finHH_2tau2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT");
  TH3F* ggWe_HH_2tau2g_noSF  = (TH3F*)finHH_2tau2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor");
  TH3F* ggWe_HH_2tau2gJECup  = (TH3F*)finHH_2tau2g->Get("SMS_ggMT_JECup_Loose_1Ele_0_1Jets");
  TH3F* ggWe_HH_2tau2gJECdown  = (TH3F*)finHH_2tau2g->Get("SMS_ggMT_JECdown_Loose_1Ele_0_1Jets");
  TH3F* ggWmu_HH_2tau2g = (TH3F*)finHH_2tau2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT");
  TH3F* ggWmu_HH_2tau2g_noSF = (TH3F*)finHH_2tau2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT_noScaleFactor");
  TH3F* ggWmu_HH_2tau2gJECup  = (TH3F*)finHH_2tau2g->Get("SMS_ggMT_JECup_Loose_1Mu_0_1Jets");
  TH3F* ggWmu_HH_2tau2gJECdown  = (TH3F*)finHH_2tau2g->Get("SMS_ggMT_JECdown_Loose_1Mu_0_1Jets");
  TH3F* ggLeftovers_HH_2tau2g = (TH3F*)finHH_2tau2g->Get("gg_SMS_Loose_lt2b_noMu_noEle_lt2jets_mChi_mBino_met");
  TH2F* h_nEvents_HH_2tau2g = (TH2F*)finHH_2tau2g->Get("SMS_mChi_mBino"); 	 
  ggWe_HH_2tau2g->Sumw2();ggWe_HH_2tau2g_noSF->Sumw2();ggWe_HH_2tau2gJECup->Sumw2();ggWe_HH_2tau2gJECdown->Sumw2();
  ggWmu_HH_2tau2g->Sumw2();ggWmu_HH_2tau2g_noSF->Sumw2();ggWmu_HH_2tau2gJECup->Sumw2();ggWmu_HH_2tau2gJECdown->Sumw2();
  ggLeftovers_HH_2tau2g->Sumw2();h_nEvents_HH_2tau2g->Sumw2();
  h_nEvents_HH_2tau2g->GetXaxis()->SetTitleOffset(0.9);		     
  h_nEvents_HH_2tau2g->GetYaxis()->SetTitleOffset(0.75);		     
  h_nEvents_HH_2tau2g->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  h_nEvents_HH_2tau2g->GetYaxis()->SetTitle("m_{ #tilde{B}} (Gev/c^{2})"); 
  c2->cd();
  h_nEvents_HH_2tau2g->Draw("COLZtext");
  TPaveText *PrelimText_HH_2tau2g = new TPaveText(.39,.86,.89,.9,"NDC");
  PrelimText_HH_2tau2g->AddText("CMS Preliminary                  SMS H(#tau#tau)H(#gamma#gamma)+2Jet");
  PrelimText_HH_2tau2g->SetFillStyle(0);
  PrelimText_HH_2tau2g->SetFillColor(0);
  PrelimText_HH_2tau2g->SetBorderSize(0);
  PrelimText_HH_2tau2g->Draw();
  EventsText->Draw();
  c2->Print("Plots/Higgs/AcceptanceEvents_SMS_TChiHH_2tau2g_2J.png");
  c2->Print("Plots/Higgs/AcceptanceEvents_SMS_TChiHH_2tau2g_2J.pdf");
  c1->cd();
  //std::string filename;
  //if(inputfilesAAW.is_open()){
  for(int i=1;i<=ggWe_HH_2tau2g->GetNbinsX();i++){
    for(int j=1;j<=ggWe_HH_2tau2g->GetNbinsY();j++){
      //while(!inputfilesAAW.eof()){
      //std::getline(inputfilesAAW,filename);
      //cout<<"file: "<<filename<<endl;
      //TFile f(filename.c_str(),"READ");
      //f.cd();
      //TString str = f.GetName();
      //int one = str.Index("_chargino");
      //one+=9;
      //int two = str.Index("_",one+1);
      //TString MCHI (str(one,two-one));
    
      //int mChi = MCHI.Atof();
      
      int mChi = ggWe_HH_2tau2g->GetXaxis()->GetBinLowEdge(i);
      int mBino= ggWe_HH_2tau2g->GetYaxis()->GetBinLowEdge(j);
      TString MCHI;ostringstream conv,conv2;conv<<mChi;MCHI=conv.str();//= static_cast<ostringstream*>( &(ostringstream() << mChi) )->str();//sprintf(MCHI,"%i",mChi);
      TString MBINO;conv2<<mBino;MBINO=conv2.str();//sprintf(MBINO,"%i",mBino);
      TH1F* RHO = (TH1F*)finHH_2tau2g->Get("rho");RHO->Sumw2();
      float nEvents=RHO->GetEntries();
      if(mChi<150)nEvents=60000.;
      else nEvents=30000.;
      int chiBin=h_nEvents_HH_2tau2g->GetXaxis()->FindBin(mChi);
      int binoBin=h_nEvents_HH_2tau2g->GetYaxis()->FindBin(mBino);
      float nEventsHist = h_nEvents_HH_2tau2g->Integral(chiBin,chiBin,binoBin,binoBin);
      if(nEventsHist==0 && ScanBin[i][j]==false)continue;
      //TH2F* ggWe_HH_2tau2g = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Ele_0_1Jets");
      //if(ggWe_HH_2tau2g){
      TH1F*  ggWe_HH_2tau2gMt=(TH1F*)ggWe_HH_2tau2g->ProjectionZ("ggWe_HH_2tau2gMt",i,i,j,j,"o");
      TH1F*  ggWe_HH_2tau2gMt_noSF=(TH1F*)ggWe_HH_2tau2g_noSF->ProjectionZ("ggWe_HH_2tau2gMt_noSF",i,i,j,j,"o");
      TH1F*  ggWe_HH_2tau2gMtJECup=(TH1F*)ggWe_HH_2tau2gJECup->ProjectionZ("ggWe_HH_2tau2gMtJECup",i,i,j,j,"o");
      TH1F*  ggWe_HH_2tau2gMtJECdown=(TH1F*)ggWe_HH_2tau2gJECdown->ProjectionZ("ggWe_HH_2tau2gMtJECdown",i,i,j,j,"o");
      float val=ggWe_HH_2tau2gMt->GetEntries();
      //cout<<"mChi: "<<mChi<<"  mBino: "<<mBino<<"  nEvents: "<<nEvents<<"  gg We events: "<<val<<endl;
      val/=nEvents;
      //cout<<"highest_gge: "<<highest_gge<<" val: "<<val;
      if(val>0 && val<lowest_gge)lowest_gge=val;
      if(val>highest_gge)highest_gge=val;
      //cout<<"  highest_gge: "<<highest_gge<<endl;
      AcceptanceWe_HH_2tau2g->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      float scale = ggWe_HH_2tau2gMt_noSF->GetEntries()/ggWe_HH_2tau2gMt_noSF->Integral(0,-1);
      //cout<<"scale: "<<scale<<endl;
      //float Norm = PhoEffScale2*x_secsHiggsino[mChi]*brFracHgg*L_int/nEventsHist;
      float Norm = PhoEffScale2*x_secsHiggsino[mChi]*bf_HH_ttgg*L_int/nEventsHist;//higgsino model, HH_2tau2g
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      /*ggWe_HH_2tau2gMt->Scale(scale);ggWe_HH_2tau2gMt->Scale(Norm);
      ggWe_HH_2tau2gMtJECup->Scale(scale);ggWe_HH_2tau2gMtJECup->Scale(Norm);
      ggWe_HH_2tau2gMtJECdown->Scale(scale);ggWe_HH_2tau2gMtJECdown->Scale(Norm);*/
      HistScale(ggWe_HH_2tau2gMt,scale*Norm);HistScale(ggWe_HH_2tau2gMtJECup,scale*Norm);HistScale(ggWe_HH_2tau2gMtJECdown,scale*Norm);
      TH1F* hist = (TH1F*)ggWe_HH_2tau2gMt->Rebin(nMTbins,"hist",MTbins);
      AddOverflowToLastBin(hist);
      TFile fLimitsSigSMS_HH_2tau2g_ggEle("HiggsLimitFiles/LimitPackage_mu_"+MCHI+"_mB_"+MBINO+"_HH_2tau2g_ggEle.root","RECREATE");
      fLimitsSigSMS_HH_2tau2g_ggEle.cd();
      //ggWe_HH_2tau2gMt->Write("h_ggWe_HH_2tau2g_MT_mChi"+MCHI+"_mBino"+MBINO);
      //ggWe_HH_2tau2gMt->Write("hMT_ggEle_MC");
      hist->Write("hMT_ggEle_MC");
      //ggWe_HH_2tau2gMtJECup->Write("hMT_ggEle_JECUp");
      //ggWe_HH_2tau2gMtJECdown->Write("hMT_ggEle_JECDown");
      //cout<<"Total expected events: "<<ggWe_HH_2tau2gMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggWe_HH_2tau2gMt->Integral(0,ggWe_HH_2tau2gMt->FindBin(50))<<endl<<"mt>50 : "<<ggWe_HH_2tau2gMt->Integral(ggWe_HH_2tau2gMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggWe_HH_2tau2gMt->Integral(ggWe_HH_2tau2gMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggWe_HH_2tau2gMt->Integral(ggWe_HH_2tau2gMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggWe_HH_2tau2gMt->Integral(ggWe_HH_2tau2gMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggWe_HH_2tau2gMt->Integral(ggWe_HH_2tau2gMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggWe_HH_2tau2gMt->Integral(ggWe_HH_2tau2gMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggWe_HH_2tau2gMt->Integral(ggWe_HH_2tau2gMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggWe_HH_2tau2gMt->FindBin(50),bin75=ggWe_HH_2tau2gMt->FindBin(75),bin100=ggWe_HH_2tau2gMt->FindBin(100),bin125=ggWe_HH_2tau2gMt->FindBin(125),bin150=ggWe_HH_2tau2gMt->FindBin(150),bin175=ggWe_HH_2tau2gMt->FindBin(175),bin200=ggWe_HH_2tau2gMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsWe_HH_2tau2g.push_back(ggWe_HH_2tau2gMt->Integral(0,-1));
      outValsWe_HH_2tau2g.push_back(ggWe_HH_2tau2gMt->Integral(0,bin50-1));
      outValsWe_HH_2tau2g.push_back(ggWe_HH_2tau2gMt->Integral(bin50,-1));
      outValsWe_HH_2tau2g.push_back(ggWe_HH_2tau2gMt->Integral(bin75,-1));
      outValsWe_HH_2tau2g.push_back(ggWe_HH_2tau2gMt->Integral(bin100,-1));
      //outValsWe_HH_2tau2g.push_back(ggWe_HH_2tau2gMt->Integral(bin125,-1));
      outValsWe_HH_2tau2g.push_back(ggWe_HH_2tau2gMt->Integral(bin150,-1));
      //outValsWe_HH_2tau2g.push_back(ggWe_HH_2tau2gMt->Integral(bin175,-1));
      outValsWe_HH_2tau2g.push_back(ggWe_HH_2tau2gMt->Integral(bin200,-1));
      
      //}
      //delete ggWe_HH_2tau2g;
      fLimitsSigSMS_HH_2tau2g_ggEle.Close();

      //TH2F* ggWmu_HH_2tau2g = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
      //if(ggWmu_HH_2tau2g){
      TH1F* ggWmu_HH_2tau2gMt=(TH1F*)ggWmu_HH_2tau2g->ProjectionZ("ggWmu_HH_2tau2gMt",i,i,j,j,"o");
      TH1F* ggWmu_HH_2tau2gMt_noSF=(TH1F*)ggWmu_HH_2tau2g_noSF->ProjectionZ("ggWmu_HH_2tau2gMt_noSF",i,i,j,j,"o");
      TH1F* ggWmu_HH_2tau2gMtJECup=(TH1F*)ggWmu_HH_2tau2gJECup->ProjectionZ("ggWmu_HH_2tau2gMtJECup",i,i,j,j,"o");
      TH1F* ggWmu_HH_2tau2gMtJECdown=(TH1F*)ggWmu_HH_2tau2gJECdown->ProjectionZ("ggWmu_HH_2tau2gMtJECdown",i,i,j,j,"o");
      float val=ggWmu_HH_2tau2gMt->GetEntries();
      val/=nEvents;
      //cout<<"highest_ggmu: "<<highest_ggmu<<" val: "<<val;
      if(val>0 && val<lowest_ggmu)lowest_ggmu=val;
      if(val>highest_ggmu)highest_ggmu=val;
      //cout<<"  highest_ggmu: "<<highest_ggmu<<endl<<endl;
      AcceptanceWmu_HH_2tau2g->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      float scale = ggWmu_HH_2tau2gMt_noSF->GetEntries()/ggWmu_HH_2tau2gMt_noSF->Integral(0,-1);
      //cout<<"scale: "<<scale<<endl;
      //float Norm = PhoEffScale2*x_secsHiggsino[mChi]*brFracHgg*L_int/nEventsHist;
      float Norm = PhoEffScale2*x_secsHiggsino[mChi]*bf_HH_ttgg*L_int/nEventsHist;//higgsino model, HH_2tau2g
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      /*ggWmu_HH_2tau2gMt->Scale(scale);ggWmu_HH_2tau2gMt->Scale(Norm);
      ggWmu_HH_2tau2gMtJECup->Scale(scale);ggWmu_HH_2tau2gMtJECup->Scale(Norm);
      ggWmu_HH_2tau2gMtJECdown->Scale(scale);ggWmu_HH_2tau2gMtJECdown->Scale(Norm);*/
      HistScale(ggWmu_HH_2tau2gMt,scale*Norm);HistScale(ggWmu_HH_2tau2gMtJECup,scale*Norm);HistScale(ggWmu_HH_2tau2gMtJECdown,scale*Norm);
      TH1F* hist = (TH1F*)ggWmu_HH_2tau2gMt->Rebin(nMTbins,"hist",MTbins);
      AddOverflowToLastBin(hist);
      TFile fLimitsSigSMS_HH_2tau2g_ggMu("HiggsLimitFiles/LimitPackage_mu_"+MCHI+"_mB_"+MBINO+"_HH_2tau2g_ggMu.root","RECREATE");
      fLimitsSigSMS_HH_2tau2g_ggMu.cd();
      //ggWmu_HH_2tau2gMt->Write("h_ggWmu_HH_2tau2g_MT_mChi"+MCHI+"_mBino"+MBINO);
      //ggWmu_HH_2tau2gMt->Write("hMT_ggMu_MC");
      hist->Write("hMT_ggMu_MC");
      //ggWmu_HH_2tau2gMtJECup->Write("hMT_ggMu_JECUp");
      //ggWmu_HH_2tau2gMtJECdown->Write("hMT_ggMu_JECDown");
      fLimitsSigSMS_HH_2tau2g.cd();
      //cout<<"Total expected events: "<<ggWmu_HH_2tau2gMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggWmu_HH_2tau2gMt->Integral(0,ggWmu_HH_2tau2gMt->FindBin(50))<<endl<<"mt>50 : "<<ggWmu_HH_2tau2gMt->Integral(ggWmu_HH_2tau2gMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggWmu_HH_2tau2gMt->Integral(ggWmu_HH_2tau2gMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggWmu_HH_2tau2gMt->Integral(ggWmu_HH_2tau2gMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggWmu_HH_2tau2gMt->Integral(ggWmu_HH_2tau2gMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggWmu_HH_2tau2gMt->Integral(ggWmu_HH_2tau2gMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggWmu_HH_2tau2gMt->Integral(ggWmu_HH_2tau2gMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggWmu_HH_2tau2gMt->Integral(ggWmu_HH_2tau2gMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggWmu_HH_2tau2gMt->FindBin(50),bin75=ggWmu_HH_2tau2gMt->FindBin(75),bin100=ggWmu_HH_2tau2gMt->FindBin(100),bin125=ggWmu_HH_2tau2gMt->FindBin(125),bin150=ggWmu_HH_2tau2gMt->FindBin(150),bin175=ggWmu_HH_2tau2gMt->FindBin(175),bin200=ggWmu_HH_2tau2gMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsWmu_HH_2tau2g.push_back(ggWmu_HH_2tau2gMt->Integral(0,-1));
      outValsWmu_HH_2tau2g.push_back(ggWmu_HH_2tau2gMt->Integral(0,bin50-1));
      outValsWmu_HH_2tau2g.push_back(ggWmu_HH_2tau2gMt->Integral(bin50,-1));
      outValsWmu_HH_2tau2g.push_back(ggWmu_HH_2tau2gMt->Integral(bin75,-1));
      outValsWmu_HH_2tau2g.push_back(ggWmu_HH_2tau2gMt->Integral(bin100,-1));
      //outValsWmu_HH_2tau2g.push_back(ggWmu_HH_2tau2gMt->Integral(bin125,-1));
      outValsWmu_HH_2tau2g.push_back(ggWmu_HH_2tau2gMt->Integral(bin150,-1));
      //outValsWmu_HH_2tau2g.push_back(ggWmu_HH_2tau2gMt->Integral(bin175,-1));
      outValsWmu_HH_2tau2g.push_back(ggWmu_HH_2tau2gMt->Integral(bin200,-1));
      // }
      
      //TH2F* ggLeftovers_HH_2tau2g = (TH2F*)f->Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
      //if(ggLeftovers_HH_2tau2g){
      TH1F* ggLeftovers_HH_2tau2gMt=(TH1F*)ggLeftovers_HH_2tau2g->ProjectionZ("ggLeftovers_HH_2tau2gMt",i,i,j,j,"o");
      float val=ggLeftovers_HH_2tau2gMt->GetEntries();
      val/=nEvents;
      //cout<<"highest_ggLeftovers: "<<highest_ggLeftovers<<" val: "<<val;
      if(val>0 && val<lowest_ggLeftovers)lowest_ggLeftovers=val;
      if(val>highest_ggLeftovers)highest_ggLeftovers=val;
      //cout<<"  highest_ggLeftovers: "<<highest_ggLeftovers<<endl<<endl;
      AcceptanceLeftovers_HH_2tau2g->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      //fLimitsSigSMS_HH_2tau2g.cd();
      //ggLeftovers_HH_2tau2gMt->Write("h_ggLeftovers_HH_2tau2g_MT_mChi"+MCHI+"_mBino"+MBINO);
      float scale = ggLeftovers_HH_2tau2gMt->GetEntries()/ggLeftovers_HH_2tau2gMt->Integral(0,-1);
      //cout<<"scale: "<<scale<<endl;
      float Norm = PhoEffScale2*x_secsHiggsino[mChi]*brFracHgg*L_int/nEventsHist;
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      /*ggLeftovers_HH_2tau2gMt->Scale(scale);
	ggLeftovers_HH_2tau2gMt->Scale(Norm);*/
      HistScale(ggLeftovers_HH_2tau2gMt,scale*Norm);  
      //cout<<"Total expected events: "<<ggLeftovers_HH_2tau2gMt->Integral(0,-1)<<endl;
      //cout<<"mt<50 : "<<ggLeftovers_HH_2tau2gMt->Integral(0,ggLeftovers_HH_2tau2gMt->FindBin(50))<<endl<<"mt>50 : "<<ggLeftovers_HH_2tau2gMt->Integral(ggLeftovers_HH_2tau2gMt->FindBin(50),-1)<<endl<<"mt>75 : "<<ggLeftovers_HH_2tau2gMt->Integral(ggLeftovers_HH_2tau2gMt->FindBin(75),-1)<<endl<<"mt>100 : "<<ggLeftovers_HH_2tau2gMt->Integral(ggLeftovers_HH_2tau2gMt->FindBin(100),-1)<<endl<<"mt>125 : "<<ggLeftovers_HH_2tau2gMt->Integral(ggLeftovers_HH_2tau2gMt->FindBin(125),-1)<<endl<<"mt>150 : "<<ggLeftovers_HH_2tau2gMt->Integral(ggLeftovers_HH_2tau2gMt->FindBin(150),-1)<<endl<<"mt>175 : "<<ggLeftovers_HH_2tau2gMt->Integral(ggLeftovers_HH_2tau2gMt->FindBin(175),-1)<<endl<<"mt>200 : "<<ggLeftovers_HH_2tau2gMt->Integral(ggLeftovers_HH_2tau2gMt->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggLeftovers_HH_2tau2gMt->FindBin(50),bin75=ggLeftovers_HH_2tau2gMt->FindBin(75),bin100=ggLeftovers_HH_2tau2gMt->FindBin(100),bin125=ggLeftovers_HH_2tau2gMt->FindBin(125),bin150=ggLeftovers_HH_2tau2gMt->FindBin(150),bin175=ggLeftovers_HH_2tau2gMt->FindBin(175),bin200=ggLeftovers_HH_2tau2gMt->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsLeftovers_HH_2tau2g.push_back(ggLeftovers_HH_2tau2gMt->Integral(0,-1));
      outValsLeftovers_HH_2tau2g.push_back(ggLeftovers_HH_2tau2gMt->Integral(0,bin50-1));
      outValsLeftovers_HH_2tau2g.push_back(ggLeftovers_HH_2tau2gMt->Integral(bin50,-1));
      outValsLeftovers_HH_2tau2g.push_back(ggLeftovers_HH_2tau2gMt->Integral(bin75,-1));
      outValsLeftovers_HH_2tau2g.push_back(ggLeftovers_HH_2tau2gMt->Integral(bin100,-1));
      //outValsLeftovers_HH_2tau2g.push_back(ggLeftovers_HH_2tau2gMt->Integral(bin125,-1));
      outValsLeftovers_HH_2tau2g.push_back(ggLeftovers_HH_2tau2gMt->Integral(bin150,-1));
      //outValsLeftovers_HH_2tau2g.push_back(ggLeftovers_HH_2tau2gMt->Integral(bin175,-1));
      outValsLeftovers_HH_2tau2g.push_back(ggLeftovers_HH_2tau2gMt->Integral(bin200,-1));
      // }
    }
  }
  
// make Acceptance plots


  
//AcceptanceWe->GetZaxis()->SetRangeUser(lowest_gge,highest_gge*1.1);
//AcceptanceWe->GetZaxis()->SetRangeUser(.017,.039);
AcceptanceWe->Draw("COLZ");
TPaveText *PrelimText2 = new TPaveText(.25,.86,.87,.9,"NDC");
PrelimText2->AddText("CMS Preliminary                  SMS W(#rightarrowInclusive)H(#gamma#gamma)+2Jet");
PrelimText2->SetFillStyle(0);
PrelimText2->SetFillColor(0);
PrelimText2->SetBorderSize(0);
PrelimText2->Draw();
PrelimText2->Draw();
TPaveText *WeText = new TPaveText(.24,.67,.64,.72,"NDC");
WeText->AddText("#gamma#gamma+electron, #leq1 b-jet");
WeText->SetFillStyle(0);
WeText->SetFillColor(0);
WeText->SetBorderSize(0);
WeText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_gge.png");
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_gge.pdf");
  
//AcceptanceWmu->GetZaxis()->SetRangeUser(lowest_ggmu,highest_ggmu*1.1);
//AcceptanceWmu->GetZaxis()->SetRangeUser(.025,.047);
AcceptanceWmu->Draw("COLZ");
PrelimText2->Draw();
TPaveText *WmuText = new TPaveText(.24,.67,.64,.72,"NDC");
WmuText->AddText("#gamma#gamma+#mu, #leq1 b-jet");
WmuText->SetFillStyle(0);
WmuText->SetFillColor(0);
WmuText->SetBorderSize(0);
WmuText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_ggmu.png");
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_ggmu.pdf");
  
//AcceptanceLeftovers->GetZaxis()->SetRangeUser(lowest_ggLeftovers,highest_ggLeftovers*1.1);
//AcceptanceLeftovers->GetZaxis()->SetRangeUser(.18,.25);
AcceptanceLeftovers->Draw("COLZ");
PrelimText2->Draw();
TPaveText *LeftoversText = new TPaveText(.24,.67,.64,.72,"NDC");
LeftoversText->AddText("#gamma#gamma, #leq1 b-jet, 0leptons, no 70<m_{jj}<110");
LeftoversText->SetFillStyle(0);
LeftoversText->SetFillColor(0);
LeftoversText->SetBorderSize(0);
LeftoversText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_ggLeftovers.png");
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_ggLeftovers.pdf");

 
//AcceptanceWe_ZH->GetZaxis()->SetRangeUser(lowest_gge,highest_gge*1.1);
//AcceptanceWe_ZH->GetZaxis()->SetRangeUser(.02,.045);
AcceptanceWe_ZH->Draw("COLZ");
TPaveText *PrelimText_ZH2 = new TPaveText(.25,.86,.87,.9,"NDC");
PrelimText_ZH2->AddText("CMS Preliminary                  SMS Z(#rightarrowInclusive)H(#gamma#gamma)+2Jet");
PrelimText_ZH2->SetFillStyle(0);
PrelimText_ZH2->SetFillColor(0);
PrelimText_ZH2->SetBorderSize(0);
PrelimText_ZH2->Draw();
PrelimText_ZH2->Draw();
WeText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiZH_ZincHgg_2J_gge.png");
 c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiZH_ZincHgg_2J_gge.pdf");
 
 //AcceptanceWmu_ZH->GetZaxis()->SetRangeUser(lowest_ggmu,highest_ggmu*1.1);
 //AcceptanceWmu_ZH->GetZaxis()->SetRangeUser(.022,.045);
 AcceptanceWmu_ZH->Draw("COLZ");
 PrelimText_ZH2->Draw();
 WmuText->Draw();
 c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiZH_ZincHgg_2J_ggmu.png");
 c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiZH_ZincHgg_2J_ggmu.pdf");
 
 //AcceptanceLeftovers_ZH->GetZaxis()->SetRangeUser(lowest_ggLeftovers,highest_ggLeftovers*1.1);
 //AcceptanceLeftovers_ZH->GetZaxis()->SetRangeUser(.47,.7);
AcceptanceLeftovers_ZH->Draw("COLZ");
PrelimText_ZH2->Draw();
LeftoversText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiZH_ZincHgg_2J_ggLeftovers.png");
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiZH_ZincHgg_2J_ggLeftovers.pdf");

//
//AcceptanceWe_HH_2b2g->GetZaxis()->SetRangeUser(lowest_gge,highest_gge*1.1);
//AcceptanceWe_HH_2b2g->GetZaxis()->SetRangeUser(0.,.0012);
AcceptanceWe_HH_2b2g->Draw("COLZ");
TPaveText *PrelimText_HH_2b2g2 = new TPaveText(.25,.86,.87,.9,"NDC");
PrelimText_HH_2b2g2->AddText("CMS Preliminary                  SMS H(bb)H(#gamma#gamma)+2Jet");
PrelimText_HH_2b2g2->SetFillStyle(0);
PrelimText_HH_2b2g2->SetFillColor(0);
PrelimText_HH_2b2g2->SetBorderSize(0);
PrelimText_HH_2b2g2->Draw();
PrelimText_HH_2b2g2->Draw();
WeText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2b2g_2J_gge.png");
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2b2g_2J_gge.pdf");

//AcceptanceWmu_HH_2b2g->GetZaxis()->SetRangeUser(lowest_ggmu,highest_ggmu*1.1);
//AcceptanceWmu_HH_2b2g->GetZaxis()->SetRangeUser(0,.002);
AcceptanceWmu_HH_2b2g->Draw("COLZ");
PrelimText_HH_2b2g2->Draw();
WmuText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2b2g_2J_ggmu.png");
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2b2g_2J_ggmu.pdf");

//AcceptanceLeftovers_HH_2b2g->GetZaxis()->SetRangeUser(lowest_ggLeftovers,highest_ggLeftovers*1.1);
//AcceptanceLeftovers_HH_2b2g->GetZaxis()->SetRangeUser(.44,.61);
AcceptanceLeftovers_HH_2b2g->Draw("COLZ");
PrelimText_HH_2b2g2->Draw();
LeftoversText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2b2g_2J_ggLeftovers.png");
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2b2g_2J_ggLeftovers.pdf");

//
//AcceptanceWe_HH_2Z2g->GetZaxis()->SetRangeUser(lowest_gge,highest_gge*1.1);
//AcceptanceWe_HH_2Z2g->GetZaxis()->SetRangeUser(.03,.06);
AcceptanceWe_HH_2Z2g->Draw("COLZ");
TPaveText *PrelimText_HH_2Z2g2 = new TPaveText(.25,.86,.87,.9,"NDC");
PrelimText_HH_2Z2g2->AddText("CMS Preliminary                  SMS H(ZZ)H(#gamma#gamma)+2Jet");
PrelimText_HH_2Z2g2->SetFillStyle(0);
PrelimText_HH_2Z2g2->SetFillColor(0);
PrelimText_HH_2Z2g2->SetBorderSize(0);
PrelimText_HH_2Z2g2->Draw();
PrelimText_HH_2Z2g2->Draw();
WeText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2Z2g_2J_gge.png");
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2Z2g_2J_gge.pdf");

//AcceptanceWmu_HH_2Z2g->GetZaxis()->SetRangeUser(lowest_ggmu,highest_ggmu*1.1);
//AcceptanceWmu_HH_2Z2g->GetZaxis()->SetRangeUser(.039,.075);
AcceptanceWmu_HH_2Z2g->Draw("COLZ");
PrelimText_HH_2Z2g2->Draw();
WmuText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2Z2g_2J_ggmu.png");
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2Z2g_2J_ggmu.pdf");

//AcceptanceLeftovers_HH_2Z2g->GetZaxis()->SetRangeUser(lowest_ggLeftovers,highest_ggLeftovers*1.1);
//AcceptanceLeftovers_HH_2Z2g->GetZaxis()->SetRangeUser(.465,.67);
AcceptanceLeftovers_HH_2Z2g->Draw("COLZ");
PrelimText_HH_2Z2g2->Draw();
LeftoversText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2Z2g_2J_ggLeftovers.png");
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2Z2g_2J_ggLeftovers.pdf");

//

//AcceptanceWe_HH_2W2g->GetZaxis()->SetRangeUser(lowest_gge,highest_gge*1.1);
//AcceptanceWe_HH_2W2g->GetZaxis()->SetRangeUser(.065,.119);
AcceptanceWe_HH_2W2g->Draw("COLZ");
TPaveText *PrelimText_HH_2W2g2 = new TPaveText(.25,.86,.87,.9,"NDC");
PrelimText_HH_2W2g2->AddText("CMS Preliminary                  SMS H(WW)H(#gamma#gamma)+2Jet");
PrelimText_HH_2W2g2->SetFillStyle(0);
PrelimText_HH_2W2g2->SetFillColor(0);
PrelimText_HH_2W2g2->SetBorderSize(0);
PrelimText_HH_2W2g2->Draw();
PrelimText_HH_2W2g2->Draw();
WeText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2W2g_2J_gge.png");
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2W2g_2J_gge.pdf");

//AcceptanceWmu_HH_2W2g->GetZaxis()->SetRangeUser(lowest_ggmu,highest_ggmu*1.1);
//AcceptanceWmu_HH_2W2g->GetZaxis()->SetRangeUser(.09,.16);
AcceptanceWmu_HH_2W2g->Draw("COLZ");
PrelimText_HH_2W2g2->Draw();
WmuText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2W2g_2J_ggmu.png");
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2W2g_2J_ggmu.pdf");

//AcceptanceLeftovers_HH_2W2g->GetZaxis()->SetRangeUser(lowest_ggLeftovers,highest_ggLeftovers*1.1);
//AcceptanceLeftovers_HH_2W2g->GetZaxis()->SetRangeUser(.41,.57);
AcceptanceLeftovers_HH_2W2g->Draw("COLZ");
PrelimText_HH_2W2g2->Draw();
LeftoversText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2W2g_2J_ggLeftovers.png");
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2W2g_2J_ggLeftovers.pdf");


//

//AcceptanceWe_HH_2tau2g->GetZaxis()->SetRangeUser(lowest_gge,highest_gge*1.1);
//AcceptanceWe_HH_2tau2g->GetZaxis()->SetRangeUser(.09,.18);
AcceptanceWe_HH_2tau2g->Draw("COLZ");
TPaveText *PrelimText_HH_2tau2g2 = new TPaveText(.25,.86,.87,.9,"NDC");
PrelimText_HH_2tau2g2->AddText("CMS Preliminary                  SMS H(#tau#tau)H(#gamma#gamma)+2Jet");
PrelimText_HH_2tau2g2->SetFillStyle(0);
PrelimText_HH_2tau2g2->SetFillColor(0);
PrelimText_HH_2tau2g2->SetBorderSize(0);
PrelimText_HH_2tau2g2->Draw();
PrelimText_HH_2tau2g2->Draw();
WeText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2tau2g_2J_gge.png");
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2tau2g_2J_gge.pdf");

//AcceptanceWmu_HH_2tau2g->GetZaxis()->SetRangeUser(lowest_ggmu,highest_ggmu*1.1);
//AcceptanceWmu_HH_2tau2g->GetZaxis()->SetRangeUser(.12,.23);
AcceptanceWmu_HH_2tau2g->Draw("COLZ");
PrelimText_HH_2tau2g2->Draw();
WmuText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2tau2g_2J_ggmu.png");
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2tau2g_2J_ggmu.pdf");

//AcceptanceLeftovers_HH_2tau2g->GetZaxis()->SetRangeUser(lowest_ggLeftovers,highest_ggLeftovers*1.1);
//AcceptanceLeftovers_HH_2tau2g->GetZaxis()->SetRangeUser(.41,.57);
AcceptanceLeftovers_HH_2tau2g->Draw("COLZ");
PrelimText_HH_2tau2g2->Draw();
LeftoversText->Draw();
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2tau2g_2J_ggLeftovers.png");
c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiHH_2tau2g_2J_ggLeftovers.pdf");

int j=0,k=0;


cout<<"\\documentclass[11pt]{article}"<<endl;
cout<<"\\usepackage{calc}"<<endl;
cout<<"\\usepackage{multirow}"<<endl;
cout<<"\\usepackage{verbatim}"<<endl;
cout<<"\\usepackage{changepage}"<<endl;
cout<<"\\usepackage{tabularx}"<<endl;
cout<<"\\begin{document}"<<endl;

cout<<" \\begin{center}"<<endl;
cout<<"  \\begin{small}"<<endl;
cout<<"  \\hspace*{-3cm}"<<endl;
cout<<"   \\begin{tabularx}{1.45\\textwidth}{ | c | c | c | c | c | c | c | c | c | c |}\n      \\hline"<<endl;
cout<<"%\n%\n%\n";
//cout<<"      \\bf{SMS ZH e+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$125} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$175} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
cout<<"      \\bf{SMS ZH e+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
  
for(int i=0;i<outValsWe_ZH.size();i++){
  if(i>0 && i%nCats==0)cout<<" \\\\ \\hline"<<endl;
  if(i%nCats==0 && k<=NcentrX){cout<<"      chargino"<<centrX[k];k++;}
  cout<<" & "<<outValsWe_ZH[i];
 }
cout<<" \\\\ \\hline"<<endl;
cout<<"%\n%\n%\n";
cout<<"   \\end{tabularx}"<<endl;
cout<<"  \\end{small}"<<endl;
cout<<"%\n%\n";

cout<<"  \\begin{small}"<<endl;
cout<<"  \\hspace*{-3cm}"<<endl;
cout<<"   \\begin{tabularx}{1.45\\textwidth}{ | c | c | c | c | c | c | c | c | c | c |}\n      \\hline"<<endl;
cout<<"%\n%\n%\n";
//cout<<"      \\bf{SMS ZH $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$125} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$175} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
cout<<"      \\bf{SMS ZH $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
k=0;
for(int i=0;i<outValsWmu_ZH.size();i++){
  if(i>0 && i%nCats==0)cout<<" \\\\ \\hline"<<endl;
  if(i%nCats==0 && k<=NcentrX){cout<<"      chargino"<<centrX[k];k++;}
  cout<<" & "<<outValsWmu_ZH[i];
 }
cout<<" \\\\ \\hline"<<endl;
cout<<"%\n%\n%\n";
cout<<"   \\end{tabularx}"<<endl;
cout<<"  \\end{small}"<<endl;
cout<<"%\n%\n";

cout<<"  \\begin{small}"<<endl;
cout<<"  \\hspace*{-3cm}"<<endl;
cout<<"   \\begin{tabularx}{1.45\\textwidth}{ | c | c | c | c | c | c | c | c | c | c |}\n      \\hline"<<endl;
cout<<"%\n%\n%\n";
//cout<<"      \\bf{SMS ZH $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$125} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$175} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
cout<<"      \\bf{SMS ZH $Leftovers$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
k=0;
for(int i=0;i<outValsLeftovers_ZH.size();i++){
  if(i>0 && i%nCats==0)cout<<" \\\\ \\hline"<<endl;
  if(i%nCats==0 && k<=NcentrX){cout<<"      chargino"<<centrX[k];k++;}
  cout<<" & "<<outValsLeftovers_ZH[i];
 }
cout<<" \\\\ \\hline"<<endl;
cout<<"%\n%\n%\n";
cout<<"   \\end{tabularx}"<<endl;
cout<<"  \\end{small}"<<endl;
cout<<"%\n%\n";
/*
  cout<<"  \\begin{small}"<<endl;
  cout<<"  \\hspace*{-3cm}"<<endl;
  cout<<"   \\begin{tabularx}{1.45\\textwidth}{ | c | c | c | c | c | c | c | c | c | c |}\n      \\hline"<<endl;
  cout<<"%\n%\n%\n";
  //cout<<"      \\bf{SMS ZH $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$125} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$175} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
  cout<<"      \\bf{SMS ZH $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
  k=0;
  for(int i=0;i<outValsZnunu.size();i++){
  if(i>0 && i%nCats==0)cout<<" \\\\ \\hline"<<endl;
  if(i%nCats==0){cout<<"      chargino"<<centrX[k];k++;}
  cout<<" & "<<outValsZnunu[i];
  }
  cout<<" \\\\ \\hline"<<endl;
  cout<<"%\n%\n%\n";
  cout<<"   \\end{tabularx}"<<endl;
  cout<<"  \\end{small}"<<endl;
  cout<<"%\n%\n";

  cout<<"  \\begin{small}"<<endl;
  cout<<"  \\hspace*{-3cm}"<<endl;
  cout<<"   \\begin{tabularx}{1.45\\textwidth}{ | c | c | c | c | c | c | c | c | c | c |}\n      \\hline"<<endl;
  cout<<"%\n%\n%\n";
  //cout<<"      \\bf{SMS ZH $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$125} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$175} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
  cout<<"      \\bf{SMS ZH $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MT$\\le$50} & \\bf{MT$\\ge$50} & \\bf{MT$\\ge$75} & \\bf{MT$\\ge$100} & \\bf{MT$\\ge$150} & \\bf{MT$\\ge$200}    \\\\ \\hline"<<endl;
  k=0;
  for(int i=0;i<outValsWWlnulnu.size();i++){
  if(i>0 && i%nCats==0)cout<<" \\\\ \\hline"<<endl;
  if(i%nCats==0){cout<<"      chargino"<<centrX[k];k++;}
  cout<<" & "<<outValsWWlnulnu[i];
  }
  cout<<" \\\\ \\hline"<<endl;
  cout<<"%\n%\n%\n";
  cout<<"   \\end{tabularx}"<<endl;
  cout<<"  \\end{small}"<<endl;
  cout<<"%\n%\n";
*/

cout<<" \\end{center}"<<endl;
cout<<"\\end{document}"<<endl;  





}

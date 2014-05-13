#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "TString.h"
#include "TH2I.h"
#include "TFile.h"
#include <map>

using namespace std;

void AcceptanceSMS_TchiWH_WincHgg_2J(){

  TFile fin("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_SMS_TChiWH_WincHgg_2J_Redo.root");
  //TFile fin("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_SMS_TChiWH_WincHgg_2J_3otherTriggers.root");
  TFile fLimitsSigSMS("signal_contamination_SMS_TChiWH_WincHgg_2J.root","RECREATE");
  //TFile fLimitsSigSMS("signal_contamination_aaz.root","RECREATE");
  //TFile fLimitsSigSMS("signal_contamination_aaWW.root","RECREATE");

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
  float lowest=999999.,highest=0.;

  //Double_t xbins[]={130,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525};
  Double_t xbins[]={125,140,162.5,187.5,212.5,237.5,262.5,287.5,312.5,337.5,362.5,387.5,412.5,437.5,462.5,487.5,512.5};
  Double_t ybins[]={1,20,25,45,50,70,75,95,100,120,125,145,150,170,175,195,200,220,225,245,250,270,275,295,300,320,325,345,350,370,400};
  int nXbin=(sizeof(xbins)/sizeof(Double_t))-1;
  int nYbin=(sizeof(ybins)/sizeof(Double_t))-1;

  float centrX[]={130,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500};
  int NcentrX = sizeof(centrX)/sizeof(float);
  int nBinsX=16;
  float low=112.5,high=512.5;

  float L_int=19499.;
  float brRat=0.00229;
  std::map<float,float> x_secsHiggsino;
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
  x_secsHiggsino[150]=1.846;
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

  vector<float> outValsWe,outValsWmu,outValsWWlnulnu,outValsZnunu,outValsInclusive;
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

  //ifstream inputfilesAAW;
  //ifstream inputfilesAAZ;
  //ifstream inputfilesAAWW;
  //inputfilesAAW.open("AcceptanceFilesAAW.txt");
  //inputfilesAAZ.open("AcceptanceFilesAAZ.txt");
  //inputfilesAAWW.open("AcceptanceFilesAAWW.txt");

  TH3F* ggWe  = (TH3F*)fin.Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_met");
  TH3F* ggWmu = (TH3F*)fin.Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_met");
  TH2F* h_nEvents = (TH2F*)fin.Get("SMS_mChi_mBino"); 	
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
      TString MCHI;ostringstream conv;conv<<mChi;MCHI=conv.str();//= static_cast<ostringstream*>( &(ostringstream() << mChi) )->str();//sprintf(MCHI,"%i",mChi);
      TString MBINO;conv<<mChi;MBINO=conv.str();//sprintf(MBINO,"%i",mBino);
      TH1F* RHO = (TH1F*)fin.Get("rho");
      float nEvents=RHO->GetEntries();
      if(mChi<150)nEvents=120000.;
      else if(mChi<=175)nEvents=60000.;
      else if(mChi<425)nEvents=30000.;
      else if(mChi>=425)nEvents=60000.;
      int chiBin=h_nEvents->GetXaxis()->FindBin(mChi);
      int binoBin=h_nEvents->GetYaxis()->FindBin(mBino);
      float nEvents = h_nEvents->Integral(chiBin,chiBin,binoBin,binoBin);
      if(nEvents==0)continue;
      //TH2F* ggWe = (TH2F*)f.Get("ggMetVsInvarMass_Loose_1Ele_0_1Jets");
      //if(ggWe){
      TH1F*  ggWeMet=(TH1F*)ggWe->ProjectionZ(" ggWeMet",i,i,j,j,"o");
      float val=ggWeMet->GetEntries();
      cout<<"mChi: "<<mChi<<"  mBino: "<<mBino<<"  nEvents: "<<nEvents<<"  gg We events: "<<val<<endl;
      val/=nEvents;
      cout<<"highest: "<<highest<<" val: "<<val;
      if(val>0 && val<lowest)lowest=val;
      if(val>highest)highest=val;
      cout<<"  highest: "<<highest<<endl;
      AcceptanceWe->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      fLimitsSigSMS.cd();
      ggWe->Write("h_ggWe_WH_met_mChi"+MCHI+"_mBino"+MBINO);
      float scale = ggWeMet->GetEntries()/ggWeMet->Integral();
      //cout<<"scale: "<<scale<<endl;
      float Norm = x_secsHiggsino[mChi]*brRat*L_int/nEvents;
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      ggWeMet->Scale(scale);
      ggWeMet->Scale(Norm);
      //cout<<"Total expected events: "<<ggWeMet->Integral()<<endl;
      //cout<<"met<50 : "<<ggWeMet->Integral(0,ggWeMet->FindBin(50))<<endl<<"met>50 : "<<ggWeMet->Integral(ggWeMet->FindBin(50),-1)<<endl<<"met>75 : "<<ggWeMet->Integral(ggWeMet->FindBin(75),-1)<<endl<<"met>100 : "<<ggWeMet->Integral(ggWeMet->FindBin(100),-1)<<endl<<"met>125 : "<<ggWeMet->Integral(ggWeMet->FindBin(125),-1)<<endl<<"met>150 : "<<ggWeMet->Integral(ggWeMet->FindBin(150),-1)<<endl<<"met>175 : "<<ggWeMet->Integral(ggWeMet->FindBin(175),-1)<<endl<<"met>200 : "<<ggWeMet->Integral(ggWeMet->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggWeMet->FindBin(50),bin75=ggWeMet->FindBin(75),bin100=ggWeMet->FindBin(100),bin125=ggWeMet->FindBin(125),bin150=ggWeMet->FindBin(150),bin175=ggWeMet->FindBin(175),bin200=ggWeMet->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsWe.push_back(ggWeMet->Integral(0,-1));
      outValsWe.push_back(ggWeMet->Integral(0,bin50-1));
      outValsWe.push_back(ggWeMet->Integral(bin50,-1));
      outValsWe.push_back(ggWeMet->Integral(bin75,-1));
      outValsWe.push_back(ggWeMet->Integral(bin100,-1));
      //outValsWe.push_back(ggWeMet->Integral(bin125,-1));
      outValsWe.push_back(ggWeMet->Integral(bin150,-1));
      //outValsWe.push_back(ggWeMet->Integral(bin175,-1));
      outValsWe.push_back(ggWeMet->Integral(bin200,-1));
      
      //}
      //delete ggWe;
      
      //TH2F* ggWmu = (TH2F*)f.Get("ggMetVsInvarMass_Loose_1Mu_0_1Jets");
      //if(ggWmu){
      TH1F* ggWmuMet=(TH1F*)ggWmu->ProjectionZ("ggWmuMet",i,i,j,j,"o");
      float val=ggWmuMet->GetEntries();
      val/=nEvents;
      cout<<"highest: "<<highest<<" val: "<<val;
      if(val>0 && val<lowest)lowest=val;
      if(val>highest)highest=val;
      cout<<"  highest: "<<highest<<endl<<endl;
      AcceptanceWmu->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      fLimitsSigSMS.cd();
      ggWmu->Write("h_ggWmu_WH_met_mChi"+MCHI+"_mBino"+MBINO);
      float scale = ggWmuMet->GetEntries()/ggWmuMet->Integral();
      //cout<<"scale: "<<scale<<endl;
      float Norm = x_secsHiggsino[mChi]*brRat*L_int/nEvents;
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      ggWmuMet->Scale(scale);
      ggWmuMet->Scale(Norm);
      //cout<<"Total expected events: "<<ggWmuMet->Integral()<<endl;
      //cout<<"met<50 : "<<ggWmuMet->Integral(0,ggWmuMet->FindBin(50))<<endl<<"met>50 : "<<ggWmuMet->Integral(ggWmuMet->FindBin(50),-1)<<endl<<"met>75 : "<<ggWmuMet->Integral(ggWmuMet->FindBin(75),-1)<<endl<<"met>100 : "<<ggWmuMet->Integral(ggWmuMet->FindBin(100),-1)<<endl<<"met>125 : "<<ggWmuMet->Integral(ggWmuMet->FindBin(125),-1)<<endl<<"met>150 : "<<ggWmuMet->Integral(ggWmuMet->FindBin(150),-1)<<endl<<"met>175 : "<<ggWmuMet->Integral(ggWmuMet->FindBin(175),-1)<<endl<<"met>200 : "<<ggWmuMet->Integral(ggWmuMet->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggWmuMet->FindBin(50),bin75=ggWmuMet->FindBin(75),bin100=ggWmuMet->FindBin(100),bin125=ggWmuMet->FindBin(125),bin150=ggWmuMet->FindBin(150),bin175=ggWmuMet->FindBin(175),bin200=ggWmuMet->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsWmu.push_back(ggWmuMet->Integral(0,-1));
      outValsWmu.push_back(ggWmuMet->Integral(0,bin50-1));
      outValsWmu.push_back(ggWmuMet->Integral(bin50,-1));
      outValsWmu.push_back(ggWmuMet->Integral(bin75,-1));
      outValsWmu.push_back(ggWmuMet->Integral(bin100,-1));
      //outValsWmu.push_back(ggWmuMet->Integral(bin125,-1));
      outValsWmu.push_back(ggWmuMet->Integral(bin150,-1));
      //outValsWmu.push_back(ggWmuMet->Integral(bin175,-1));
      outValsWmu.push_back(ggWmuMet->Integral(bin200,-1));
      // }
      // delete ggWmu;
      
      
      //f.Close();
    }
  }
  
  /*
    if(inputfilesAAZ.is_open()){
    
    while(!inputfilesAAZ.eof()){
    std::getline(inputfilesAAZ,filename);
    cout<<"file: "<<filename<<endl;
    TFile f(filename.c_str(),"READ");
    f.cd();
    TString str = f.GetName();
    int one = str.Index("_chargino");
    one+=9;
    int two = str.Index("_",one+1);
    TString MCHI (str(one,two-one));
    
    int mChi = MCHI.Atof();
      
    TH1F* RHO = (TH1F*)f.Get("rho");
    float nEvents=RHO->GetEntries();
   
    TH2F* ggZnunu = (TH2F*)f.Get("ggMetVsInvarMass_Loose_0Lep_0_1Jets");
    if(ggZnunu){
    float val=ggZnunu->GetEntries();
    TH1F* ggZnunuMet=(TH1F*)ggZnunu->ProjectionY("ggZnunuMet");
    val/=nEvents;
    if(val<lowest)lowest=val;
    if(val>highest)highest=val;
    AcceptanceZnunu->Fill(mChi,val); 
    //cout<<mS<<endl<<mG<<endl<<val<<endl;
    fLimitsSigAAZ.cd();
    ggZnunu->Write("h_ggZnunu_WH_met_mChi"+MCHI);
    float scale = ggZnunuMet->GetEntries()/ggZnunuMet->Integral();
    //cout<<"scale: "<<scale<<endl;
    float Norm = x_secsAAZ[mChi]*brRat*L_int/nEvents;
    //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
    //cout<<"Norm: "<<Norm<<endl;
    ggZnunuMet->Scale(scale);
    ggZnunuMet->Scale(Norm);
    //cout<<"Total expected events: "<<ggZnunuMet->Integral()<<endl;
    //cout<<"met<50 : "<<ggZnunuMet->Integral(0,ggZnunuMet->FindBin(50))<<endl<<"met>50 : "<<ggZnunuMet->Integral(ggZnunuMet->FindBin(50),-1)<<endl<<"met>75 : "<<ggZnunuMet->Integral(ggZnunuMet->FindBin(75),-1)<<endl<<"met>100 : "<<ggZnunuMet->Integral(ggZnunuMet->FindBin(100),-1)<<endl<<"met>125 : "<<ggZnunuMet->Integral(ggZnunuMet->FindBin(125),-1)<<endl<<"met>150 : "<<ggZnunuMet->Integral(ggZnunuMet->FindBin(150),-1)<<endl<<"met>175 : "<<ggZnunuMet->Integral(ggZnunuMet->FindBin(175),-1)<<endl<<"met>200 : "<<ggZnunuMet->Integral(ggZnunuMet->FindBin(200),-1)<<endl<<endl;

    float bin50=ggZnunuMet->FindBin(50),bin75=ggZnunuMet->FindBin(75),bin100=ggZnunuMet->FindBin(100),bin125=ggZnunuMet->FindBin(125),bin150=ggZnunuMet->FindBin(150),bin175=ggZnunuMet->FindBin(175),bin200=ggZnunuMet->FindBin(200);

    //if add or remove any, must change nCats at top
    outValsZnunu.push_back(ggZnunuMet->Integral(0,-1));
    outValsZnunu.push_back(ggZnunuMet->Integral(0,bin50-1));
    outValsZnunu.push_back(ggZnunuMet->Integral(bin50,-1));
    outValsZnunu.push_back(ggZnunuMet->Integral(bin75,-1));
    outValsZnunu.push_back(ggZnunuMet->Integral(bin100,-1));
    //outValsZnunu.push_back(ggZnunuMet->Integral(bin125,-1));
    outValsZnunu.push_back(ggZnunuMet->Integral(bin150,-1));
    //outValsZnunu.push_back(ggZnunuMet->Integral(bin175,-1));
    outValsZnunu.push_back(ggZnunuMet->Integral(bin200,-1));
    }
    delete ggZnunu;
  
    f.Close();
    }
    }

    if(inputfilesAAWW.is_open()){
    
    while(!inputfilesAAWW.eof()){
    std::getline(inputfilesAAWW,filename);
    cout<<"file: "<<filename<<endl;
    TFile f(filename.c_str(),"READ");
    f.cd();
    TString str = f.GetName();
    int one = str.Index("_chargino");
    one+=9;
    int two = str.Index("_",one+1);
    TString MCHI (str(one,two-one));
    
    int mChi = MCHI.Atof();
      
    TH1F* RHO = (TH1F*)f.Get("rho");
    float nEvents=RHO->GetEntries();
    TH2F* ggWWlnulnu = (TH2F*)f.Get("ggMetVsInvarMass_Loose_2LepOSoffZ_0_1Jets");
    if(ggWWlnulnu){
    float val=ggWWlnulnu->GetEntries();
    TH1F* ggWWlnulnuMet=(TH1F*)ggWWlnulnu->ProjectionY("ggWWlnulnuMet");
    val/=nEvents;
    if(val<lowest)lowest=val;
    if(val>highest)highest=val;
    AcceptanceWWlnulnu->Fill(mChi,val); 
    //cout<<mS<<endl<<mG<<endl<<val<<endl;
    fLimitsSigSMS.cd();
    ggWWlnulnu->Write("h_ggWWlnulnu_WH_met_mChi"+MCHI);
    float scale = ggWWlnulnuMet->GetEntries()/ggWWlnulnuMet->Integral();
    //cout<<"scale: "<<scale<<endl;
    float Norm = x_secsAAWW[mChi]*brRat*L_int/nEvents;
    //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
    //cout<<"Norm: "<<Norm<<endl;
    ggWWlnulnuMet->Scale(scale);
    ggWWlnulnuMet->Scale(Norm);
    //cout<<"Total expected events: "<<ggWWlnulnuMet->Integral()<<endl;
    //cout<<"met<50 : "<<ggWWlnulnuMet->Integral(0,ggWWlnulnuMet->FindBin(50))<<endl<<"met>50 : "<<ggWWlnulnuMet->Integral(ggWWlnulnuMet->FindBin(50),-1)<<endl<<"met>75 : "<<ggWWlnulnuMet->Integral(ggWWlnulnuMet->FindBin(75),-1)<<endl<<"met>100 : "<<ggWWlnulnuMet->Integral(ggWWlnulnuMet->FindBin(100),-1)<<endl<<"met>125 : "<<ggWWlnulnuMet->Integral(ggWWlnulnuMet->FindBin(125),-1)<<endl<<"met>150 : "<<ggWWlnulnuMet->Integral(ggWWlnulnuMet->FindBin(150),-1)<<endl<<"met>175 : "<<ggWWlnulnuMet->Integral(ggWWlnulnuMet->FindBin(175),-1)<<endl<<"met>200 : "<<ggWWlnulnuMet->Integral(ggWWlnulnuMet->FindBin(200),-1)<<endl<<endl;

    float bin50=ggWWlnulnuMet->FindBin(50),bin75=ggWWlnulnuMet->FindBin(75),bin100=ggWWlnulnuMet->FindBin(100),bin125=ggWWlnulnuMet->FindBin(125),bin150=ggWWlnulnuMet->FindBin(150),bin175=ggWWlnulnuMet->FindBin(175),bin200=ggWWlnulnuMet->FindBin(200);

    //if add or remove any, must change nCats at top
    outValsWWlnulnu.push_back(ggWWlnulnuMet->Integral(0,-1));
    outValsWWlnulnu.push_back(ggWWlnulnuMet->Integral(0,bin50-1));
    outValsWWlnulnu.push_back(ggWWlnulnuMet->Integral(bin50,-1));
    outValsWWlnulnu.push_back(ggWWlnulnuMet->Integral(bin75,-1));
    outValsWWlnulnu.push_back(ggWWlnulnuMet->Integral(bin100,-1));
    //outValsWWlnulnu.push_back(ggWWlnulnuMet->Integral(bin125,-1));
    outValsWWlnulnu.push_back(ggWWlnulnuMet->Integral(bin150,-1));
    //outValsWWlnulnu.push_back(ggWWlnulnuMet->Integral(bin175,-1));
    outValsWWlnulnu.push_back(ggWWlnulnuMet->Integral(bin200,-1));
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
    cout<<"No Jet Req Bin "<<bin<<" has 0 value , mS="<<AcceptanceWe->GetXaxis()->GetBinLowEdge(i)<<" mG="<<AcceptanceWe->GetYaxis()->GetBinLowEdge(j)<<"  Value set to:"<<SchmearedVal<<endl;
    }
    }
    }
  */
  /*
  //highest=.012;
  //AcceptanceWe->SetMinimum(0);
  //AcceptanceWe->SetMaximum(highest*1.1);
  AcceptanceWe->GetZaxis()->SetRangeUser(lowest,highest*1.1);//SetMaximum(highest*1.1);
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
  WeText->AddText("#gamma#gamma+electron selection");
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
  //AcceptanceWmu->SetMaximum(highest*1.1);
  AcceptanceWmu->GetZaxis()->SetRangeUser(lowest,highest*1.1);//SetMaximum(highest*1.1);
  //AcceptanceWmu->SetMinimum(.22);
  AcceptanceWmu->Draw("COLZ");
  PrelimText2->Draw();
  TPaveText *WmuText = new TPaveText(.24,.67,.64,.72,"NDC");
  WmuText->AddText("#gamma#gamma+#mu selection");
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
  //cout<<"      \\bf{SMS WZ e+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$125} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$175} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
  cout<<"      \\bf{SMS WZ e+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
  
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
  //cout<<"      \\bf{SMS WZ $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$125} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$175} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
  cout<<"      \\bf{SMS WZ $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
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
  /*
    cout<<"  \\begin{small}"<<endl;
    cout<<"  \\hspace*{-3cm}"<<endl;
    cout<<"   \\begin{tabularx}{1.45\\textwidth}{ | c | c | c | c | c | c | c | c | c | c |}\n      \\hline"<<endl;
    cout<<"%\n%\n%\n";
    //cout<<"      \\bf{SMS WZ $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$125} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$175} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
    cout<<"      \\bf{SMS WZ $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
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
    //cout<<"      \\bf{SMS WZ $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$125} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$175} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
    cout<<"      \\bf{SMS WZ $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
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

  TFile finZH("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_SMS_TChiZH_ZincHgg_2J.root");
  TFile fLimitsSigSMS_ZH("signal_contamination_SMS_TChiZH_ZincHgg_2J.root","RECREATE");

  //lowest=999999.;highest=0.;

 
  vector<float> outValsWe_ZH,outValsWmu_ZH;
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

  TH3F* ggWe_ZH  = (TH3F*)finZH.Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_met");
  TH3F* ggWmu_ZH = (TH3F*)finZH.Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_met");
  TH2F* h_nEvents_ZH = (TH2F*)finZH.Get("SMS_mChi_mBino"); 	
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
      TString MCHI;ostringstream conv;conv<<mChi;MCHI=conv.str();//= static_cast<ostringstream*>( &(ostringstream() << mChi) )->str();//sprintf(MCHI,"%i",mChi);
      TString MBINO;conv<<mChi;MBINO=conv.str();//sprintf(MBINO,"%i",mBino);
      TH1F* RHO = (TH1F*)finZH.Get("rho");
      float nEvents=RHO->GetEntries();
      if(mChi<150)nEvents=120000.;
      else if(mChi<=175)nEvents=60000.;
      else if(mChi<250)nEvents=16000.;
      else if(mChi>=250)nEvents=30000.;
      int chiBin=h_nEvents_ZH->GetXaxis()->FindBin(mChi);
      int binoBin=h_nEvents_ZH->GetYaxis()->FindBin(mBino);
      float nEvents = h_nEvents_ZH->Integral(chiBin,chiBin,binoBin,binoBin);
      if(nEvents==0)continue;
      //TH2F* ggWe_ZH = (TH2F*)f.Get("ggMetVsInvarMass_Loose_1Ele_0_1Jets");
      //if(ggWe_ZH){
      TH1F*  ggWe_ZHMet=(TH1F*)ggWe_ZH->ProjectionZ(" ggWe_ZHMet",i,i,j,j,"o");
      float val=ggWe_ZHMet->GetEntries();
      cout<<"mChi: "<<mChi<<"  mBino: "<<mBino<<"  nEvents: "<<nEvents<<"  gg We events: "<<val<<endl;
      val/=nEvents;
      cout<<"highest: "<<highest<<" val: "<<val;
      if(val>0 && val<lowest)lowest=val;
      if(val>highest)highest=val;
      cout<<"  highest: "<<highest<<endl;
      AcceptanceWe_ZH->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      fLimitsSigSMS_ZH.cd();
      ggWe_ZH->Write("h_ggWe_ZH_met_mChi"+MCHI+"_mBino"+MBINO);
      float scale = ggWe_ZHMet->GetEntries()/ggWe_ZHMet->Integral();
      //cout<<"scale: "<<scale<<endl;
      float Norm = x_secsHiggsino[mChi]*brRat*L_int/nEvents;
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      ggWe_ZHMet->Scale(scale);
      ggWe_ZHMet->Scale(Norm);
      //cout<<"Total expected events: "<<ggWe_ZHMet->Integral()<<endl;
      //cout<<"met<50 : "<<ggWe_ZHMet->Integral(0,ggWe_ZHMet->FindBin(50))<<endl<<"met>50 : "<<ggWe_ZHMet->Integral(ggWe_ZHMet->FindBin(50),-1)<<endl<<"met>75 : "<<ggWe_ZHMet->Integral(ggWe_ZHMet->FindBin(75),-1)<<endl<<"met>100 : "<<ggWe_ZHMet->Integral(ggWe_ZHMet->FindBin(100),-1)<<endl<<"met>125 : "<<ggWe_ZHMet->Integral(ggWe_ZHMet->FindBin(125),-1)<<endl<<"met>150 : "<<ggWe_ZHMet->Integral(ggWe_ZHMet->FindBin(150),-1)<<endl<<"met>175 : "<<ggWe_ZHMet->Integral(ggWe_ZHMet->FindBin(175),-1)<<endl<<"met>200 : "<<ggWe_ZHMet->Integral(ggWe_ZHMet->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggWe_ZHMet->FindBin(50),bin75=ggWe_ZHMet->FindBin(75),bin100=ggWe_ZHMet->FindBin(100),bin125=ggWe_ZHMet->FindBin(125),bin150=ggWe_ZHMet->FindBin(150),bin175=ggWe_ZHMet->FindBin(175),bin200=ggWe_ZHMet->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsWe_ZH.push_back(ggWe_ZHMet->Integral(0,-1));
      outValsWe_ZH.push_back(ggWe_ZHMet->Integral(0,bin50-1));
      outValsWe_ZH.push_back(ggWe_ZHMet->Integral(bin50,-1));
      outValsWe_ZH.push_back(ggWe_ZHMet->Integral(bin75,-1));
      outValsWe_ZH.push_back(ggWe_ZHMet->Integral(bin100,-1));
      //outValsWe_ZH.push_back(ggWe_ZHMet->Integral(bin125,-1));
      outValsWe_ZH.push_back(ggWe_ZHMet->Integral(bin150,-1));
      //outValsWe_ZH.push_back(ggWe_ZHMet->Integral(bin175,-1));
      outValsWe_ZH.push_back(ggWe_ZHMet->Integral(bin200,-1));
      
      //}
      //delete ggWe_ZH;
      
      //TH2F* ggWmu_ZH = (TH2F*)f.Get("ggMetVsInvarMass_Loose_1Mu_0_1Jets");
      //if(ggWmu_ZH){
      TH1F* ggWmu_ZHMet=(TH1F*)ggWmu_ZH->ProjectionZ("ggWmu_ZHMet",i,i,j,j,"o");
      float val=ggWmu_ZHMet->GetEntries();
      val/=nEvents;
      cout<<"highest: "<<highest<<" val: "<<val;
      if(val>0 && val<lowest)lowest=val;
      if(val>highest)highest=val;
      cout<<"  highest: "<<highest<<endl<<endl;
      AcceptanceWmu_ZH->Fill(mChi,mBino,val); 
      //cout<<mS<<endl<<mG<<endl<<val<<endl;
      fLimitsSigSMS_ZH.cd();
      ggWmu_ZH->Write("h_ggWmu_ZH_met_mChi"+MCHI+"_mBino"+MBINO);
      float scale = ggWmu_ZHMet->GetEntries()/ggWmu_ZHMet->Integral();
      //cout<<"scale: "<<scale<<endl;
      float Norm = x_secsHiggsino[mChi]*brRat*L_int/nEvents;
      //cout<<"x_secs["<<mChi<<"]: "<<x_secs[mChi]<<endl;
      //cout<<"Norm: "<<Norm<<endl;
      ggWmu_ZHMet->Scale(scale);
      ggWmu_ZHMet->Scale(Norm);
      //cout<<"Total expected events: "<<ggWmu_ZHMet->Integral()<<endl;
      //cout<<"met<50 : "<<ggWmu_ZHMet->Integral(0,ggWmu_ZHMet->FindBin(50))<<endl<<"met>50 : "<<ggWmu_ZHMet->Integral(ggWmu_ZHMet->FindBin(50),-1)<<endl<<"met>75 : "<<ggWmu_ZHMet->Integral(ggWmu_ZHMet->FindBin(75),-1)<<endl<<"met>100 : "<<ggWmu_ZHMet->Integral(ggWmu_ZHMet->FindBin(100),-1)<<endl<<"met>125 : "<<ggWmu_ZHMet->Integral(ggWmu_ZHMet->FindBin(125),-1)<<endl<<"met>150 : "<<ggWmu_ZHMet->Integral(ggWmu_ZHMet->FindBin(150),-1)<<endl<<"met>175 : "<<ggWmu_ZHMet->Integral(ggWmu_ZHMet->FindBin(175),-1)<<endl<<"met>200 : "<<ggWmu_ZHMet->Integral(ggWmu_ZHMet->FindBin(200),-1)<<endl<<endl;
      
      float bin50=ggWmu_ZHMet->FindBin(50),bin75=ggWmu_ZHMet->FindBin(75),bin100=ggWmu_ZHMet->FindBin(100),bin125=ggWmu_ZHMet->FindBin(125),bin150=ggWmu_ZHMet->FindBin(150),bin175=ggWmu_ZHMet->FindBin(175),bin200=ggWmu_ZHMet->FindBin(200);
      
      //if add or remove any, must change nCats at top
      outValsWmu_ZH.push_back(ggWmu_ZHMet->Integral(0,-1));
      outValsWmu_ZH.push_back(ggWmu_ZHMet->Integral(0,bin50-1));
      outValsWmu_ZH.push_back(ggWmu_ZHMet->Integral(bin50,-1));
      outValsWmu_ZH.push_back(ggWmu_ZHMet->Integral(bin75,-1));
      outValsWmu_ZH.push_back(ggWmu_ZHMet->Integral(bin100,-1));
      //outValsWmu_ZH.push_back(ggWmu_ZHMet->Integral(bin125,-1));
      outValsWmu_ZH.push_back(ggWmu_ZHMet->Integral(bin150,-1));
      //outValsWmu_ZH.push_back(ggWmu_ZHMet->Integral(bin175,-1));
      outValsWmu_ZH.push_back(ggWmu_ZHMet->Integral(bin200,-1));
      // }
      
      
    }
  }
  

  //highest=.012;
  //AcceptanceWe->SetMinimum(0);
  //AcceptanceWe->SetMaximum(highest*1.1);
  AcceptanceWe->GetZaxis()->SetRangeUser(lowest,highest*1.1);//SetMaximum(highest*1.1);
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
  WeText->AddText("#gamma#gamma+electron selection");
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
  //AcceptanceWmu->SetMaximum(highest*1.1);
  AcceptanceWmu->GetZaxis()->SetRangeUser(lowest,highest*1.1);//SetMaximum(highest*1.1);
  //AcceptanceWmu->SetMinimum(.22);
  AcceptanceWmu->Draw("COLZ");
  PrelimText2->Draw();
  TPaveText *WmuText = new TPaveText(.24,.67,.64,.72,"NDC");
  WmuText->AddText("#gamma#gamma+#mu selection");
  ////WmuText->AddText("");
  WmuText->SetFillStyle(0);
  WmuText->SetFillColor(0);
  WmuText->SetBorderSize(0);
  WmuText->Draw();
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_ggmu.png");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiWH_WincHgg_2J_ggmu.pdf");

  //highest=.012;
  //AcceptanceWe_ZH->SetMinimum(0);
  //AcceptanceWe_ZH->SetMaximum(highest*1.1);
  AcceptanceWe_ZH->GetZaxis()->SetRangeUser(lowest,highest*1.1);//SetMaximum(highest*1.1);
  //AcceptanceWe_ZH->SetMinimum(.22);
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
  /*
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiZH_ZincHgg_2J_gge_3otherTrigs.png");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiZH_ZincHgg_2J_gge_3otherTrigs.pdf");
  */
  //AcceptanceWmu_ZH->SetMinimum(0);
  //AcceptanceWmu_ZH->SetMaximum(highest*1.1);
  AcceptanceWmu_ZH->GetZaxis()->SetRangeUser(lowest,highest*1.1);//SetMaximum(highest*1.1);
  //AcceptanceWmu_ZH->SetMinimum(.22);
  AcceptanceWmu_ZH->Draw("COLZ");
  PrelimText_ZH2->Draw();
  WmuText->Draw();
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiZH_ZincHgg_2J_ggmu.png");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_SMS_TChiZH_ZincHgg_2J_ggmu.pdf");

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
  //cout<<"      \\bf{SMS ZH e+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$125} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$175} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
  cout<<"      \\bf{SMS ZH e+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
  
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
  //cout<<"      \\bf{SMS ZH $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$125} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$175} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
  cout<<"      \\bf{SMS ZH $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
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
  /*
    cout<<"  \\begin{small}"<<endl;
    cout<<"  \\hspace*{-3cm}"<<endl;
    cout<<"   \\begin{tabularx}{1.45\\textwidth}{ | c | c | c | c | c | c | c | c | c | c |}\n      \\hline"<<endl;
    cout<<"%\n%\n%\n";
    //cout<<"      \\bf{SMS ZH $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$125} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$175} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
    cout<<"      \\bf{SMS ZH $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
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
    //cout<<"      \\bf{SMS ZH $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$125} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$175} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
    cout<<"      \\bf{SMS ZH $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
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

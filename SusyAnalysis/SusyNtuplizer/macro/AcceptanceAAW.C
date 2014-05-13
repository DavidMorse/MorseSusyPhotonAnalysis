#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "TString.h"
#include "TH2I.h"
#include "TFile.h"
#include <map>

using namespace std;

void AcceptanceAAW(){

  TFile fLimitsSigAAW("signal_contamination_aaw.root","RECREATE");
  //TFile fLimitsSigAAW("signal_contamination_aaz.root","RECREATE");
  //TFile fLimitsSigAAW("signal_contamination_aaWW.root","RECREATE");

  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","",800,700);
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

  float centr[]={130,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500};
  int nBins=16;
  float low=112.5,high=512.5;

  float L_int=19499.;
  float brRat=0.00229;
  std::map<float,float> x_secsAAW;
  x_secsAAW[130]=4.12;
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

  vector<float> outValsWe,outValsWmu,outValsWWlnulnu,outValsZnunu,outValsInclusive;
  int nCats=7;

  //TH2F* AcceptanceWe = new TH2F("AcceptanceWe","",25,112.5,512.5,1,0,1);  			     
  TH1F* AcceptanceWe = new TH1F("AcceptanceWe","",nBins,low,high);  			     
  AcceptanceWe->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceWe->GetYaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceWe->GetYaxis()->SetTitleOffset(1.6);
  AcceptanceWe->GetXaxis()->SetTitleOffset(.9);
  AcceptanceWe->GetXaxis()->SetLabelSize(0.04);
  AcceptanceWe->GetYaxis()->SetLabelSize(0.05);
  AcceptanceWe->SetLineColor(kBlue); AcceptanceWe->SetFillColor(kBlue);  			     
  TH1F* AcceptanceWmu = new TH1F("AcceptanceWmu","",nBins,low,high);  			     
  AcceptanceWmu->GetXaxis()->SetTitle("m_{ #tilde{#Chi}} (Gev/c^{2})");
  AcceptanceWmu->GetYaxis()->SetTitle("Acceptance x Efficiency"); 
  AcceptanceWmu->GetYaxis()->SetTitleOffset(1.6);
  AcceptanceWmu->GetXaxis()->SetTitleOffset(.9);
  AcceptanceWmu->GetXaxis()->SetLabelSize(0.04);
  AcceptanceWmu->GetYaxis()->SetLabelSize(0.05);
  AcceptanceWmu->SetLineColor(kBlue); AcceptanceWmu->SetFillColor(kBlue); 

  ifstream inputfilesAAW;
  //ifstream inputfilesAAZ;
  //ifstream inputfilesAAWW;
  inputfilesAAW.open("AcceptanceFilesAAW.txt");
  //inputfilesAAZ.open("AcceptanceFilesAAZ.txt");
  //inputfilesAAWW.open("AcceptanceFilesAAWW.txt");

  std::string filename;
  if(inputfilesAAW.is_open()){
    
    while(!inputfilesAAW.eof()){
      std::getline(inputfilesAAW,filename);
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
      //if(mChi<150)nEvents=120000.;
      //else if(mChi<=175)nEvents=60000.;
      //else if(mChi<=400)nEvents=30000.;
      //else if(mChi<=500)nEvents=30000.;
      

      TH2F* ggWe = (TH2F*)f.Get("ggMetVsInvarMass_Loose_1Ele_0_1Jets");
      if(ggWe){
	float val=ggWe->GetEntries();
	//cout<<"mChi: "<<mChi<<"  nEvents: "<<nEvents<<"  gg We events: "<<val<<endl;
	TH1F* ggWeMet=(TH1F*)ggWe->ProjectionY("ggWeMet");
	val/=nEvents;
	if(val<lowest)lowest=val;
	if(val>highest)highest=val;
	AcceptanceWe->Fill(mChi,val); 
	//cout<<mS<<endl<<mG<<endl<<val<<endl;
	fLimitsSigAAW.cd();
	ggWe->Write("h_ggWe_met_mChi"+MCHI);
	float scale = ggWeMet->GetEntries()/ggWeMet->Integral();
	//cout<<"scale: "<<scale<<endl;
	float Norm = x_secsAAW[mChi]*brRat*L_int/nEvents;
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

      }
      delete ggWe;
  
      TH2F* ggWmu = (TH2F*)f.Get("ggMetVsInvarMass_Loose_1Mu_0_1Jets");
      if(ggWmu){
	float val=ggWmu->GetEntries();
	TH1F* ggWmuMet=(TH1F*)ggWmu->ProjectionY("ggWmuMet");
	val/=nEvents;
	if(val<lowest)lowest=val;
	if(val>highest)highest=val;
	AcceptanceWmu->Fill(mChi,val); 
	//cout<<mS<<endl<<mG<<endl<<val<<endl;
	fLimitsSigAAW.cd();
	ggWmu->Write("h_ggWmu_met_mChi"+MCHI);
	float scale = ggWmuMet->GetEntries()/ggWmuMet->Integral();
	//cout<<"scale: "<<scale<<endl;
	float Norm = x_secsAAW[mChi]*brRat*L_int/nEvents;
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
      }
      delete ggWmu;
  

      f.Close();
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
    ggZnunu->Write("h_ggZnunu_met_mChi"+MCHI);
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
    fLimitsSigAAW.cd();
    ggWWlnulnu->Write("h_ggWWlnulnu_met_mChi"+MCHI);
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

  AcceptanceWe->SetMinimum(0);
  AcceptanceWe->SetMaximum(highest*1.1);
  //AcceptanceWe->SetMinimum(.22);
  AcceptanceWe->Draw("histo");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_aaW_We.png");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_aaW_We.pdf");

  AcceptanceWmu->SetMinimum(0);
  AcceptanceWmu->SetMaximum(highest*1.1);
  //AcceptanceWmu->SetMinimum(.22);
  AcceptanceWmu->Draw("histo");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_aaW_Wmu.png");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_aaW_Wmu.pdf");
  /*
  AcceptanceZnunu->SetMinimum(0);
  AcceptanceZnunu->SetMaximum(highest*1.1);
  //AcceptanceZnunu->SetMinimum(.22);
  AcceptanceZnunu->Draw("histo");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_aaW_Znunu.png");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_aaW_Znunu.pdf");

  AcceptanceWWlnulnu->SetMinimum(0);
  AcceptanceWWlnulnu->SetMaximum(highest*1.1);
  //AcceptanceWWlnulnu->SetMinimum(.22);
  AcceptanceWWlnulnu->Draw("histo");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_aaW_WWlnulnu.png");
  c1->Print("Plots/Higgs/AcceptanceTimesEfficiency_aaW_WWlnulnu.pdf");
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
  //cout<<"      \\bf{AAW e+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$125} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$175} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
  cout<<"      \\bf{AAW e+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
  
  for(int i=0;i<outValsWe.size();i++){
    if(i>0 && i%nCats==0)cout<<" \\\\ \\hline"<<endl;
    if(i%nCats==0){cout<<"      chargino"<<centr[k];k++;}
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
  //cout<<"      \\bf{AAW $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$125} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$175} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
  cout<<"      \\bf{AAW $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
  k=0;
  for(int i=0;i<outValsWmu.size();i++){
    if(i>0 && i%nCats==0)cout<<" \\\\ \\hline"<<endl;
    if(i%nCats==0){cout<<"      chargino"<<centr[k];k++;}
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
    //cout<<"      \\bf{AAW $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$125} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$175} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
    cout<<"      \\bf{AAW $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
    k=0;
    for(int i=0;i<outValsZnunu.size();i++){
    if(i>0 && i%nCats==0)cout<<" \\\\ \\hline"<<endl;
    if(i%nCats==0){cout<<"      chargino"<<centr[k];k++;}
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
    //cout<<"      \\bf{AAW $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$125} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$175} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
    cout<<"      \\bf{AAW $\\mu$+$\\gamma\\gamma$ Point} & \\bf{Total Yield} & \\bf{MET$\\le$50} & \\bf{MET$\\ge$50} & \\bf{MET$\\ge$75} & \\bf{MET$\\ge$100} & \\bf{MET$\\ge$150} & \\bf{MET$\\ge$200}    \\\\ \\hline"<<endl;
    k=0;
    for(int i=0;i<outValsWWlnulnu.size();i++){
    if(i>0 && i%nCats==0)cout<<" \\\\ \\hline"<<endl;
    if(i%nCats==0){cout<<"      chargino"<<centr[k];k++;}
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

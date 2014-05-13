#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "TString.h"
#include "TH2I.h"
#include "TFile.h"

using namespace std;

void Acceptance(){

  TFile fLimitsSigBino("signal_contamination_bino_chi0375.root","RECREATE");
  TFile fLimitsSigWino("signal_contamination_wino_chi0375.root","RECREATE");


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


 //First Bino:

  TH2F* Acceptance = new TH2F("Acceptance","",17,400,2100,17,420,2120);  			     
  Acceptance->GetXaxis()->SetTitle("m_{ #tilde{q}} (Gev/c^{2})");
  Acceptance->GetYaxis()->SetTitle("m_{ #tilde{g}} (Gev/c^{2})"); 
  Acceptance->GetYaxis()->SetTitleOffset(1.6);
  Acceptance->GetXaxis()->SetTitleOffset(.9);
  Acceptance->GetXaxis()->SetLabelSize(0.04);
  Acceptance->GetYaxis()->SetLabelSize(0.05);

  TH2F* Acceptance_JetReq = new TH2F("Acceptance_JetReq","",17,400,2100,17,420,2120);  							     
  Acceptance_JetReq->GetXaxis()->SetTitle("m_{ #tilde{q}} (Gev/c^{2})");
  Acceptance_JetReq->GetYaxis()->SetTitle("m_{ #tilde{g}} (Gev/c^{2})");
  Acceptance_JetReq->GetYaxis()->SetTitleOffset(1.6);
  Acceptance_JetReq->GetXaxis()->SetTitleOffset(.9);
  Acceptance_JetReq->GetXaxis()->SetLabelSize(0.04);
  Acceptance_JetReq->GetYaxis()->SetLabelSize(0.05);

  vector<TFile> files;

  ifstream inputfiles;
  inputfiles.open("AcceptanceFilesBino.txt");

  std::string filename;
  if(inputfiles.is_open()){
    
    while(!inputfiles.eof()){
      std::getline(inputfiles,filename);
      cout<<"file: "<<filename<<endl;
      TFile f(filename.c_str(),"READ");
      f.cd();
      TString str = f.GetName();
      int one = str.Index("_Bino_");
      one+=6;
      int two = str.Index("_",one+1);
      TString MS (str(one,two-one));
      int three = str.Index("_",two+1);
      TString MG (str(two+1,three-two-1));
      int four = str.Index(".",three+1);
      TString MCHI (str(three+1,four-three-1));
      
      int mG = MG.Atof();
      int mS = MS.Atof();
      int mChi = MCHI.Atof();
      
      //cout<<mG<<endl<<mS<<endl<<mChi<<endl<<endl;
      
      TH1F* gg = (TH1F*)f.Get("ggMet");
      TH1F* gg_JetReq = (TH1F*)f.Get("ggMet_JetReq");
      TH1F* ff = (TH1F*)f.Get("ffMet");
      TH1F* ff_JetReq = (TH1F*)f.Get("ffMet_JetReq");
      if(gg){
	float val=gg->GetEntries();
	val/=10000.;
	float val_JetReq=gg_JetReq->GetEntries();
	val_JetReq/=10000.;
	if(val<lowest)lowest=val;
	if(val>highest)highest=val;
	if(val_JetReq<lowest)lowest=val_JetReq;
	if(val_JetReq>highest)highest=val_JetReq;
	Acceptance->Fill(mS,mG,val); 
	Acceptance_JetReq->Fill(mS,mG,val_JetReq); 
	//cout<<mS<<endl<<mG<<endl<<val<<endl;
	fLimitsSigBino.cd();
	gg->Write("h_gg_met_nojet_mS"+MS+"_mG"+MG+"_mN"+MCHI);
	gg_JetReq->Write("h_gg_met_1jet_mS"+MS+"_mG"+MG+"_mN"+MCHI);
	ff->Write("h_ff_met_nojet_mS"+MS+"_mG"+MG+"_mN"+MCHI);
	ff_JetReq->Write("h_ff_met_1jet_mS"+MS+"_mG"+MG+"_mN"+MCHI);
      }
      delete gg;f.Close();
    }
  }

  for(int i=1;i<Acceptance->GetNbinsX();i++){
    for(int j=1;j<Acceptance->GetNbinsY();j++){
      int bin = Acceptance->GetBin(i,j);
      if(Acceptance->GetBinContent(bin)==0){
	int binRight = Acceptance->GetBin(i+1,j);
	int binLeft = Acceptance->GetBin(i-1,j);
	int binUp = Acceptance->GetBin(i,j+1);
	int binDown = Acceptance->GetBin(i,j-1);
	float SchmearedVal = (Acceptance->GetBinContent(binRight)+Acceptance->GetBinContent(binLeft)+Acceptance->GetBinContent(binUp)+Acceptance->GetBinContent(binDown))/4.;
	Acceptance->SetBinContent(bin,SchmearedVal);
	cout<<"No Jet Req Bin "<<bin<<" has 0 value , mS="<<Acceptance->GetXaxis()->GetBinLowEdge(i)<<" mG="<<Acceptance->GetYaxis()->GetBinLowEdge(j)<<"  Value set to:"<<SchmearedVal<<endl;
      }
      if(Acceptance_JetReq->GetBinContent(bin)==0){
	int binRight_JetReq = Acceptance_JetReq->GetBin(i+1,j);
	int binLeft_JetReq = Acceptance_JetReq->GetBin(i-1,j);
	int binUp_JetReq = Acceptance_JetReq->GetBin(i,j+1);
	int binDown_JetReq = Acceptance_JetReq->GetBin(i,j-1);
	float SchmearedVal_JetReq = (Acceptance_JetReq->GetBinContent(binRight_JetReq)+Acceptance_JetReq->GetBinContent(binLeft_JetReq)+Acceptance_JetReq->GetBinContent(binUp_JetReq)+Acceptance_JetReq->GetBinContent(binDown_JetReq))/4.;
	Acceptance_JetReq->SetBinContent(bin,SchmearedVal_JetReq);
	cout<<"1+ Jet Req Bin "<<bin<<" has 0 value , mS="<<Acceptance_JetReq->GetXaxis()->GetBinLowEdge(i)<<" mG="<<Acceptance_JetReq->GetYaxis()->GetBinLowEdge(j)<<"  Value set to:"<<SchmearedVal_JetReq<<endl;
      }
    }
  }


  Acceptance->SetMinimum(lowest-.01);
  Acceptance->SetMaximum(highest+.01);
  //Acceptance->SetMinimum(.22);
  Acceptance->Draw("colz");
  c1->Print("Plots/Figure22_AcceptanceTimesEfficiency_Bino.png");
  c1->Print("Plots/Figure22_AcceptanceTimesEfficiency_Bino.pdf");
  Acceptance_JetReq->SetMinimum(lowest-.01);
  Acceptance_JetReq->SetMaximum(highest+.01);
  Acceptance_JetReq->Draw("colz");
  c1->Print("Plots/Figure22_AcceptanceTimesEfficiency_Bino_JetReq.png");
  c1->Print("Plots/Figure22_AcceptanceTimesEfficiency_Bino_JetReq.pdf");


  //now Wino:
  lowest=999999.;highest=0.;

  TH2F* AcceptanceWino = new TH2F("AcceptanceWino","",16,500,2100,17,420,2120);  			     
  AcceptanceWino->GetXaxis()->SetTitle("m_{ #tilde{q}} (Gev/c^{2})");
  AcceptanceWino->GetYaxis()->SetTitle("m_{ #tilde{g}} (Gev/c^{2})"); 
  AcceptanceWino->GetYaxis()->SetTitleOffset(1.5);
  AcceptanceWino->GetXaxis()->SetTitleOffset(.85);
  AcceptanceWino->GetXaxis()->SetLabelSize(0.04);
  AcceptanceWino->GetYaxis()->SetLabelSize(0.05);

  TH2F* AcceptanceWino_JetReq = new TH2F("AcceptanceWino_JetReq","",16,500,2100,17,420,2120);  						     
  AcceptanceWino_JetReq->GetXaxis()->SetTitle("m_{ #tilde{q}} (Gev/c^{2})");
  AcceptanceWino_JetReq->GetYaxis()->SetTitle("m_{ #tilde{g}} (Gev/c^{2})");
  AcceptanceWino_JetReq->GetYaxis()->SetTitleOffset(1.5);
  AcceptanceWino_JetReq->GetXaxis()->SetTitleOffset(.85);
  AcceptanceWino_JetReq->GetXaxis()->SetLabelSize(0.04);
  AcceptanceWino_JetReq->GetYaxis()->SetLabelSize(0.05);

  vector<TFile> filesWino;

  ifstream inputfilesWino;
  inputfilesWino.open("AcceptanceFilesWino.txt");

  std::string filenameWino;
  if(inputfilesWino.is_open()){
    
    while(!inputfilesWino.eof()){
      std::getline(inputfilesWino,filenameWino);
      cout<<"file: "<<filenameWino<<endl;
      TFile f(filenameWino.c_str(),"READ");
      f.cd();
      TString str = f.GetName();
      int one = str.Index("_Wino_");
      one+=6;
      int two = str.Index("_",one+1);
      TString MS (str(one,two-one));
      int three = str.Index("_",two+1);
      TString MG (str(two+1,three-two-1));
      int four = str.Index(".",three+1);
      TString MCHI (str(three+1,four-three-1));
      
      int mG = MG.Atof();
      int mS = MS.Atof();
      int mChi = MCHI.Atof();
      
      //cout<<mG<<endl<<mS<<endl<<mChi<<endl<<endl;
      
      TH1F* ggWino = (TH1F*)f.Get("ggMet");
      TH1F* ggWino_JetReq = (TH1F*)f.Get("ggMet_JetReq");
      TH1F* ffWino = (TH1F*)f.Get("ffMet");
      TH1F* ffWino_JetReq = (TH1F*)f.Get("ffMet_JetReq");
      if(ggWino){
	float val=ggWino->GetEntries();
	//	if(mG==1720 && mS==500)cout<<"Acceptance: 
	val/=60000.;
	float val_JetReq=ggWino_JetReq->GetEntries();
	val_JetReq/=60000.;
	if(val<lowest)lowest=val;
	if(val>highest)highest=val;
	if(val_JetReq<lowest)lowest=val_JetReq;
	if(val_JetReq>highest)highest=val_JetReq;
	AcceptanceWino->Fill(mS,mG,val); 
	AcceptanceWino_JetReq->Fill(mS,mG,val_JetReq); 
	//cout<<mS<<endl<<mG<<endl<<val<<endl;
	fLimitsSigWino.cd();
	ggWino->Write("h_gg_met_nojet_mS"+MS+"_mG"+MG+"_mN"+MCHI);
	ggWino_JetReq->Write("h_gg_met_1jet_mS"+MS+"_mG"+MG+"_mN"+MCHI);
	ffWino->Write("h_ff_met_nojet_mS"+MS+"_mG"+MG+"_mN"+MCHI);
	ffWino_JetReq->Write("h_ff_met_1jet_mS"+MS+"_mG"+MG+"_mN"+MCHI);
      }
      delete ggWino;f.Close();
    }
  }

  for(int i=1;i<AcceptanceWino->GetNbinsX();i++){
    for(int j=1;j<AcceptanceWino->GetNbinsY();j++){
      int bin = AcceptanceWino->GetBin(i,j);
      if(AcceptanceWino->GetBinContent(bin)==0){
	int binRight = AcceptanceWino->GetBin(i+1,j);
	int binLeft = AcceptanceWino->GetBin(i-1,j);
	int binUp = AcceptanceWino->GetBin(i,j+1);
	int binDown = AcceptanceWino->GetBin(i,j-1);
	float SchmearedVal = (AcceptanceWino->GetBinContent(binRight)+AcceptanceWino->GetBinContent(binLeft)+AcceptanceWino->GetBinContent(binUp)+AcceptanceWino->GetBinContent(binDown))/4.;
	AcceptanceWino->SetBinContent(bin,SchmearedVal);
	cout<<"No Jet Req Bin "<<bin<<" has 0 value , mS="<<AcceptanceWino->GetXaxis()->GetBinLowEdge(i)<<" mG="<<AcceptanceWino->GetYaxis()->GetBinLowEdge(j)<<"  Value set to:"<<SchmearedVal<<endl;
      }
      if(AcceptanceWino_JetReq->GetBinContent(bin)==0){
	int binRight_JetReq = AcceptanceWino_JetReq->GetBin(i+1,j);
	int binLeft_JetReq = AcceptanceWino_JetReq->GetBin(i-1,j);
	int binUp_JetReq = AcceptanceWino_JetReq->GetBin(i,j+1);
	int binDown_JetReq = AcceptanceWino_JetReq->GetBin(i,j-1);
	float SchmearedVal_JetReq = (AcceptanceWino_JetReq->GetBinContent(binRight_JetReq)+AcceptanceWino_JetReq->GetBinContent(binLeft_JetReq)+AcceptanceWino_JetReq->GetBinContent(binUp_JetReq)+AcceptanceWino_JetReq->GetBinContent(binDown_JetReq))/4.;
	AcceptanceWino_JetReq->SetBinContent(bin,SchmearedVal_JetReq);
	cout<<"1+ Jet Req Bin "<<bin<<" has 0 value , mS="<<AcceptanceWino_JetReq->GetXaxis()->GetBinLowEdge(i)<<" mG="<<AcceptanceWino_JetReq->GetYaxis()->GetBinLowEdge(j)<<"  Value set to:"<<SchmearedVal_JetReq<<endl;
      }
    }
  }


  AcceptanceWino->SetMinimum(lowest-.0005);
  AcceptanceWino->SetMaximum(highest+.001);
  //AcceptanceWino->SetMinimum(.22);
  AcceptanceWino->Draw("colz");
  //titleRight->Draw("SAME");
  //p1->Draw();
  c1->Print("Plots/Figure22_AcceptanceTimesEfficiency_Wino.png");
  c1->Print("Plots/Figure22_AcceptanceTimesEfficiency_Wino.pdf");
  AcceptanceWino_JetReq->SetMinimum(lowest-.0005);
  AcceptanceWino_JetReq->SetMaximum(highest+.001);
  AcceptanceWino_JetReq->Draw("colz");
  //titleRight->Draw();
  c1->Print("Plots/Figure22_AcceptanceTimesEfficiency_Wino_JetReq.png");
  c1->Print("Plots/Figure22_AcceptanceTimesEfficiency_Wino_JetReq.pdf");
  /*
  for(int i=0;i<AcceptanceWino->GetNbinsX()+1;i++){
    for(int j=0;j<AcceptanceWino->GetNbinsY()+1;j++){
      int bin=AcceptanceWino->GetBin(i,j);
      cout<<"i:"<<i<<" j:"<<j<<"  "<<AcceptanceWino->GetXaxis()->GetBinLowEdge(i)<<"  "<<AcceptanceWino->GetYaxis()->GetBinLowEdge(j)<<"  "<<AcceptanceWino->GetBinContent(i,j)<<endl;
    }
  }
  */


}

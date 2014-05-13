#include <TFile>
#include <TH1>
#include <TH1F>
#include <math.h>
#include <TMath.h>
#include <utility>

void DR03(){

  bool useData=1,useQCD=0;

  TFile fSig("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_SignalMC_1400_1420_375_2012PU_Ntuplized_RhoPileupCorr_8TeV_DR03.root");
  //if(useData)TFile fBack("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Data2012_Filter_JsonHLTnVertexTwo40-25GeVBarrelPhosWithHoverE-RhoPileupCorr_Photon_cms525v3_Analysis_Runs190456-195397_DR03.root");
  if(useData)TFile fBack("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Data2012_Filter_JsonHLTnVertexTwo40-25GeVBarrelPhosWithHoverE_Photon_cms525v3_jec2012_Analysis_2012_Runs190456-195947_DR03.root");
  if(useQCD)TFile fBack("/data/ndpc3/c/dmorse/RA3/AnalysisOutput/hist_Background_Efficiencies_QCD_Pt15to1800_TuneZ2_7TeV_pythia6_DR03_NoReweight_RequireNotPho.root");
  TFile fout("IsolationOptimization2012.root","RECREATE");

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->Draw();
  c1->SetLogy(0);

  gStyle->SetOptStat(0);

  //----first Ecal------
  TH1F *sig_ecalDR03 = fSig.Get("EcalIsoDR03Nminus1Signal");
  if(useData)TH1F *back_ecalDR03 = fBack.Get("BgroundEcalIsoDR03Nminus1");//data
  if(useQCD)TH1F *back_ecalDR03 = fBack.Get("EcalIsoDR03Nminus1Bground");//qcd


  sig_ecalDR03->SetLineColor(kBlue);
  back_ecalDR03->SetLineColor(kRed);

  TH1F * sig_ecalDR03New = sig_ecalDR03->Clone();
  TH1F* back_ecalDR03New = back_ecalDR03->Clone();

  sig_ecalDR03New->Scale(1/sig_ecalDR03New->Integral());
  back_ecalDR03New->Scale(1/back_ecalDR03New->Integral());
  sig_ecalDR03New->SetTitle("Fraction of Events VS EcalIsoDR03");
  sig_ecalDR03New->GetYaxis()->SetTitle("Fraction of Events");
  sig_ecalDR03New->Draw();
  back_ecalDR03New->Draw("SAME");
  TLegend *DR03Iso = new TLegend(.45,.36,.85,.66,"","brNDC");
  DR03Iso->AddEntry(sig_ecalDR03,"EcalIsoDR03 SignalMC_2000_1015_305","lf");
  if(useQCD)DR03Iso->AddEntry(back_ecalDR03,"EcalIsoDR03 QCD","lf");
  if(useData)DR03Iso->AddEntry(back_ecalDR03,"EcalIsoDR03 Data MET<30","lf");
  DR03Iso->SetFillColor(kWhite);
  // DR03Iso->SetTextSize(0.03);
  DR03Iso->Draw("SAME");
  c1->Print("Plots/DR03/EcalDR03Iso.png");


  TH1F *sigEffEcalDR03 = sig_ecalDR03->Clone();
  TH1F *backEffEcalDR03 = back_ecalDR03->Clone();
  sigEffEcalDR03->SetTitle("Integrated Fraction of Events VS EcalIsoDR03");
  sigEffEcalDR03->GetYaxis()->SetTitle("Integrated Fraction of Events");
  sigEffEcalDR03->GetXaxis()->SetTitle("EcalIsoDR03");
  sigEffEcalDR03->GetYaxis()->SetTitleOffset(0.75);
  sigEffEcalDR03->GetXaxis()->SetTitleOffset(0.75);
  sigEffEcalDR03->SetLineColor(kBlue);
  backEffEcalDR03->SetLineColor(kRed);

  for(int i=1;i<151;i++){
    sigEffEcalDR03->SetBinContent(i,sig_ecalDR03->Integral(0,i)/sig_ecalDR03->Integral());
    backEffEcalDR03->SetBinContent(i,back_ecalDR03->Integral(0,i)/back_ecalDR03->Integral());
  }


  sigEffEcalDR03->Draw();
  backEffEcalDR03->Draw("SAME");
  TLegend *DR03EffEcal = new TLegend(.35,.3,.8,.5,"","brNDC");
  DR03EffEcal->AddEntry(sigEffEcalDR03,"EcalIsoDR03 SignalMC_2000_1015_305","lf");
  if(useQCD)DR03EffEcal->AddEntry(backEffEcalDR03,"EcalIsoDR03 QCD","lf");
  if(useData)DR03EffEcal->AddEntry(backEffEcalDR03,"EcalIsoDR03 Data MET<30","lf");
  DR03EffEcal->SetFillColor(kWhite);
  //DR03EffEcal->SetTextSize(0.03);
  DR03EffEcal->Draw("SAME");
  c1->Print("Plots/DR03/EcalIsoEffDR03.png");

  float DR03Sig[150]={0},DR03Bg[150]={0},DR03SoverRootB[150]={0},DR03LogSoverLogB[150]={0};
  for(int i=1;i<151;i++){
    DR03Sig[i-1]=sigEffEcalDR03->GetBinContent(i);
    DR03Bg[i-1]=backEffEcalDR03->GetBinContent(i);
    DR03SoverRootB[i-1]=DR03Sig[i-1]/sqrt(DR03Bg[i-1]);
    DR03LogSoverLogB[i-1]=log(DR03Sig[i-1])/log(DR03Bg[i-1]);
  }

  TGraph *bVsSEcalDR03 = new TGraph(150,DR03Bg,DR03Sig);
  if(useQCD)bVsSEcalDR03->SetTitle("Integrated Spectra - SignalMC VS QCD - ECAL");
  if(useData)bVsSEcalDR03->SetTitle("Integrated Spectra - SignalMC VS Data MET<30 - ECAL");
  bVsSEcalDR03->GetYaxis()->SetTitle("Signal MC");
  if(useQCD)bVsSEcalDR03->GetXaxis()->SetTitle("Background QCD");
  if(useData)bVsSEcalDR03->GetXaxis()->SetTitle("Background Data MET<30");
  bVsSEcalDR03->SetMarkerColor(kBlue);
  bVsSEcalDR03->SetMarkerStyle(20);
  bVsSEcalDR03->SetMarkerSize(.5);
  bVsSEcalDR03->Draw("AP");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(bVsSEcalDR03,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalVsBackgroundEffEcal.png");

  TH1F *SoverRootBvsEcalIsoDR03 =  sig_ecalDR03->Clone();
  pair<int,float> max;max.first=0;max.second=0.;
  for(int i=1;i<151;i++){
    if(TMath::IsNaN(DR03SoverRootB[i-1]) || !TMath::Finite(DR03SoverRootB[i-1]) ) SoverRootBvsEcalIsoDR03->SetBinContent(i,1);
    else {
      SoverRootBvsEcalIsoDR03->SetBinContent(i,DR03SoverRootB[i-1]);
      if(DR03SoverRootB[i-1]>max.second){
	max.first=i;max.second=DR03SoverRootB[i-1];
      }
    }
  }
  cout<<"Ecal SoverRootB Max:  (Bin,LowEdge,NextlowEdge) ("<<max.first<<","<<SoverRootBvsEcalIsoDR03->GetBinLowEdge(max.first)<<","<<SoverRootBvsEcalIsoDR03->GetBinLowEdge(max.first+1)<<")   Value: "<<max.second<<endl;

  SoverRootBvsEcalIsoDR03->SetTitle("S/#sqrt{B} Vs  EcalIsoDR03");
  SoverRootBvsEcalIsoDR03->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  SoverRootBvsEcalIsoDR03->GetXaxis()->SetTitle("EcalIsoDR03"); 
  SoverRootBvsEcalIsoDR03->GetYaxis()->SetTitleOffset(0.75);
  SoverRootBvsEcalIsoDR03->GetXaxis()->SetTitleOffset(0.75);
  SoverRootBvsEcalIsoDR03->SetMarkerColor(kBlue);
  SoverRootBvsEcalIsoDR03->SetMarkerSize(.5);
  SoverRootBvsEcalIsoDR03->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(SoverRootBvsEcalIsoDR03,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalEffOverRootBgroundEffVsEcalIso.png");


  TH1F *LogSoverLogBvsEcalIsoDR03 =  sig_ecalDR03->Clone();
  for(int i=1;i<151;i++){
    //cout<<"bin:"<<i<<"   "<<DR03LogSoverLogB[i-1]<<endl;
    if(TMath::IsNaN(DR03LogSoverLogB[i-1]) || !TMath::Finite(DR03LogSoverLogB[i-1]) ) LogSoverLogBvsEcalIsoDR03->SetBinContent(i,1);
    else LogSoverLogBvsEcalIsoDR03->SetBinContent(i,DR03LogSoverLogB[i-1]);
  }
  LogSoverLogBvsEcalIsoDR03->SetTitle("log(S)/log(B) Vs  EcalIsoDR03");
  LogSoverLogBvsEcalIsoDR03->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  LogSoverLogBvsEcalIsoDR03->GetXaxis()->SetTitle("EcalIsoDR03");
  LogSoverLogBvsEcalIsoDR03->GetYaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsEcalIsoDR03->GetXaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsEcalIsoDR03->SetMarkerColor(kBlue);
  LogSoverLogBvsEcalIsoDR03->SetMarkerSize(.5);
  LogSoverLogBvsEcalIsoDR03->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(LogSoverLogBvsEcalIsoDR03,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/LogSignalEffOverLogBgroundEffVsEcalIso.png");
  
  //----now Hcal------
  TH1F *sig_hcalDR03 = fSig.Get("HcalIsoDR03Nminus1Signal");
  if(useData)TH1F *back_hcalDR03 = fBack.Get("BgroundHcalIsoDR03Nminus1");//data
  if(useQCD) TH1F *back_hcalDR03 = fBack.Get("HcalIsoDR03Nminus1Bground");//QCD


  sig_hcalDR03->SetLineColor(kBlue);
  back_hcalDR03->SetLineColor(kRed);

  TH1F * sig_hcalDR03New = sig_hcalDR03->Clone();
  TH1F* back_hcalDR03New = back_hcalDR03->Clone();

  sig_hcalDR03New->Scale(1/sig_hcalDR03New->Integral());
  back_hcalDR03New->Scale(1/back_hcalDR03New->Integral());
  back_hcalDR03New->SetTitle("Fraction of Events VS HcalIsoDR03");
  back_hcalDR03New->GetYaxis()->SetTitle("Fraction of Events");
  back_hcalDR03New->Draw();
  sig_hcalDR03New->Draw("SAME");
  TLegend *DR03Iso = new TLegend(.45,.36,.85,.66,"","brNDC");
  DR03Iso->AddEntry(sig_hcalDR03,"HcalIsoDR03 SignalMC_2000_1015_305","lf");
  if(useQCD)DR03Iso->AddEntry(back_hcalDR03,"HcalIsoDR03 QCD","lf");
  if(useData)DR03Iso->AddEntry(back_hcalDR03,"HcalIsoDR03 Data MET<30","lf");
  DR03Iso->SetFillColor(kWhite);
  // DR03Iso->SetTextSize(0.03);
  DR03Iso->Draw("SAME");
  c1->Print("Plots/DR03/HcalDR03Iso.png");


  TH1F *sigEffHcalDR03 = sig_hcalDR03->Clone();
  TH1F *backEffHcalDR03 =  back_hcalDR03->Clone();
  sigEffHcalDR03->SetTitle("Integrated Fraction of Events VS HcalIsoDR03");
  sigEffHcalDR03->GetYaxis()->SetTitle("Integrated Fraction of Events");
  sigEffHcalDR03->GetXaxis()->SetTitle("HcalIsoDR03");
  sigEffHcalDR03->GetYaxis()->SetTitleOffset(0.75);
  sigEffHcalDR03->GetXaxis()->SetTitleOffset(0.75);
  sigEffHcalDR03->SetLineColor(kBlue);
  backEffHcalDR03->SetLineColor(kRed);

  for(int i=1;i<151;i++){
    sigEffHcalDR03->SetBinContent(i,sig_hcalDR03->Integral(0,i)/sig_hcalDR03->Integral());
    backEffHcalDR03->SetBinContent(i,back_hcalDR03->Integral(0,i)/back_hcalDR03->Integral());
    /*  if(i<26){
	sigEffHcalDR03->SetBinContent(i,sig_hcalDR03->Integral(0,i)/(sig_hcalDR03->Integral()-sig_hcalDR03->GetBinContent(26)));
	backEffHcalDR03->SetBinContent(i,back_hcalDR03->Integral(0,i)/(back_hcalDR03->Integral()-back_hcalDR03->GetBinContent(26)));
	}  
	if(i>26){
	sigEffHcalDR03->SetBinContent(i,(sig_hcalDR03->Integral(0,i)-sig_hcalDR03->GetBinContent(26))/(sig_hcalDR03->Integral()-sig_hcalDR03->GetBinContent(26)));
	backEffHcalDR03->SetBinContent(i,(back_hcalDR03->Integral(0,i)-back_hcalDR03->GetBinContent(26))/(back_hcalDR03->Integral()-back_hcalDR03->GetBinContent(26)));
	}*/
  }
  /* sigEffHcalDR03->SetBinContent(26,0);
     backEffHcalDR03->SetBinContent(26,0);*/
  
  sigEffHcalDR03->Draw();
  backEffHcalDR03->Draw("SAME");
  TLegend *DR03EffHcal = new TLegend(.35,.3,.8,.5,"","brNDC");
  DR03EffHcal->AddEntry(sigEffHcalDR03,"HcalIsoDR03 SignalMC_2000_1015_305","lf");
  if(useQCD)DR03EffHcal->AddEntry(backEffHcalDR03,"HcalIsoDR03 QCD","lf");
  if(useData)DR03EffHcal->AddEntry(backEffHcalDR03,"HcalIsoDR03 Data MET<30","lf");
  DR03EffHcal->SetFillColor(kWhite);
  //DR03EffHcal->SetTextSize(0.03);
  DR03EffHcal->Draw("SAME");
  c1->Print("Plots/DR03/HcalIsoEffDR03.png");

  float HDR03Sig[150]={0},HDR03Bg[150]={0},HDR03SoverRootB[150]={0},HDR03LogSoverLogB[150]={0};
  for(int i=1;i<151;i++){
    HDR03Sig[i-1]=sigEffHcalDR03->GetBinContent(i);
    HDR03Bg[i-1]=backEffHcalDR03->GetBinContent(i);
    HDR03SoverRootB[i-1]=HDR03Sig[i-1]/sqrt(HDR03Bg[i-1]);
    HDR03LogSoverLogB[i-1]=log(HDR03Sig[i-1])/log(HDR03Bg[i-1]);
  }

  TGraph *bVsSHcalDR03 = new TGraph(150,HDR03Bg,HDR03Sig);
  if(useQCD)bVsSHcalDR03->SetTitle("Integrated Spectra - SignalMC VS QCD - HCAL");
  if(useData)bVsSHcalDR03->SetTitle("Integrated Spectra - SignalMC VS Data MET<30 - HCAL");
  bVsSHcalDR03->GetYaxis()->SetTitle("Signal MC");
  if(useQCD)bVsSHcalDR03->GetXaxis()->SetTitle("Background QCD");
  if(useData)bVsSHcalDR03->GetXaxis()->SetTitle("Background Data MET<30");
  bVsSHcalDR03->SetMarkerColor(kBlue);
  bVsSHcalDR03->SetMarkerStyle(20);
  bVsSHcalDR03->SetMarkerSize(.5);
  bVsSHcalDR03->Draw("AP");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(bVsSHcalDR03,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalVsBackgroundEffHcal.png");

  TH1F *SoverRootBvsHcalIsoDR03 =  sig_hcalDR03->Clone();
  pair<int,float> max;max.first=0;max.second=0.;
  for(int i=1;i<151;i++){
    if(TMath::IsNaN(HDR03SoverRootB[i-1]) || !TMath::Finite(HDR03SoverRootB[i-1])) SoverRootBvsHcalIsoDR03->SetBinContent(i,1);
    else {
      SoverRootBvsHcalIsoDR03->SetBinContent(i,HDR03SoverRootB[i-1]);
      if(HDR03SoverRootB[i-1]>max.second && i< 120){
	max.first=i;max.second=HDR03SoverRootB[i-1];
      }
    }
  }
  cout<<"Hcal SoverRootB Max:  (Bin,LowEdge,NextlowEdge) ("<<max.first<<","<<SoverRootBvsHcalIsoDR03->GetBinLowEdge(max.first)<<","<<SoverRootBvsHcalIsoDR03->GetBinLowEdge(max.first+1)<<")   Value: "<<max.second<<endl;

  SoverRootBvsHcalIsoDR03->SetTitle("S/#sqrt{B} Vs  HcalIsoDR03");
  SoverRootBvsHcalIsoDR03->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  SoverRootBvsHcalIsoDR03->GetXaxis()->SetTitle("HcalIsoDR03");
  SoverRootBvsHcalIsoDR03->GetYaxis()->SetTitleOffset(0.75);
  SoverRootBvsHcalIsoDR03->GetXaxis()->SetTitleOffset(0.75);
  SoverRootBvsHcalIsoDR03->SetMarkerColor(kBlue);
  SoverRootBvsHcalIsoDR03->SetMarkerSize(.5);
  SoverRootBvsHcalIsoDR03->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(SoverRootBvsHcalIsoDR03,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalEffOverRootBgroundEffVsHcalIso.png");
  
  TH1F *LogSoverLogBvsHcalIsoDR03 =  sig_hcalDR03->Clone();
  for(int i=1;i<151;i++){
    //cout<<"bin:"<<i<<"   "<<DR03LogSoverLogB[i-1]<<endl;
    if(TMath::IsNaN(HDR03LogSoverLogB[i-1]) || !TMath::Finite(HDR03LogSoverLogB[i-1]) ) LogSoverLogBvsHcalIsoDR03->SetBinContent(i,1);
    else LogSoverLogBvsHcalIsoDR03->SetBinContent(i,HDR03LogSoverLogB[i-1]);
  }
  LogSoverLogBvsHcalIsoDR03->SetTitle("log(S)/log(B) Vs  HcalIsoDR03");
  LogSoverLogBvsHcalIsoDR03->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  LogSoverLogBvsHcalIsoDR03->GetXaxis()->SetTitle("HcalIsoDR03");
  LogSoverLogBvsHcalIsoDR03->GetYaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsHcalIsoDR03->GetXaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsHcalIsoDR03->SetMarkerColor(kBlue);
  LogSoverLogBvsHcalIsoDR03->SetMarkerSize(.5);
  LogSoverLogBvsHcalIsoDR03->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(LogSoverLogBvsHcalIsoDR03,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/LogSignalEffOverLogBgroundEffVsHcalIso.png");
  
  //----now Track------
  TH1F *sig_trackDR03 = fSig.Get("TrackIsoDR03Nminus1Signal");
  if(useData)TH1F *back_trackDR03 = fBack.Get("BgroundTrackIsoDR03Nminus1");//data
  if(useQCD) TH1F *back_trackDR03 = fBack.Get("TrackIsoDR03Nminus1Bground");//qcd


  sig_trackDR03->SetLineColor(kBlue);
  back_trackDR03->SetLineColor(kRed);

  TH1F * sig_trackDR03New = sig_trackDR03->Clone();
  TH1F* back_trackDR03New = back_trackDR03->Clone();

  sig_trackDR03New->Scale(1/sig_trackDR03New->Integral());
  back_trackDR03New->Scale(1/back_trackDR03New->Integral());
  back_trackDR03New->SetTitle("Fraction of Events VS TrackIsoDR03");
  back_trackDR03New->GetYaxis()->SetTitle("Fraction of Events");
  back_trackDR03New->Draw();
  sig_trackDR03New->Draw("SAME");
  TLegend *DR03Iso = new TLegend(.45,.36,.85,.66,"","brNDC");
  DR03Iso->AddEntry(sig_trackDR03,"TrackIsoDR03 SignalMC_2000_1015_305","lf");
  if(useQCD)DR03Iso->AddEntry(back_trackDR03,"TrackIsoDR03 QCD","lf");
  if(useData)DR03Iso->AddEntry(back_trackDR03,"TrackIsoDR03 Data MET<30","lf");
  DR03Iso->SetFillColor(kWhite);
  // DR03Iso->SetTextSize(0.03);
  DR03Iso->Draw("SAME");
  c1->Print("Plots/DR03/TrackDR03Iso.png");


  TH1F *sigEffTrackDR03 = sig_trackDR03->Clone();
  TH1F *backEffTrackDR03 =  back_trackDR03->Clone();
  sigEffTrackDR03->SetTitle("Integrated Fraction of Events VS TrackIsoDR03");
  sigEffTrackDR03->GetYaxis()->SetTitle("Integrated Fraction of Events");
  sigEffTrackDR03->GetXaxis()->SetTitle("TrackIsoDR03");
  sigEffTrackDR03->GetYaxis()->SetTitleOffset(0.75);
  sigEffTrackDR03->GetXaxis()->SetTitleOffset(0.75);
  sigEffTrackDR03->SetLineColor(kBlue);
  backEffTrackDR03->SetLineColor(kRed);

  for(int i=1;i<151;i++){
    sigEffTrackDR03->SetBinContent(i,sig_trackDR03->Integral(0,i)/sig_trackDR03->Integral());
    backEffTrackDR03->SetBinContent(i,back_trackDR03->Integral(0,i)/back_trackDR03->Integral());
  }


  sigEffTrackDR03->Draw();
  backEffTrackDR03->Draw("SAME");
  TLegend *DR03EffTrack = new TLegend(.35,.3,.8,.5,"","brNDC");
  DR03EffTrack->AddEntry(sigEffTrackDR03,"TrackIsoDR03 SignalMC_2000_1015_305","lf");
  if(useQCD)DR03EffTrack->AddEntry(backEffTrackDR03,"TrackIsoDR03 QCD","lf");
  if(useData)DR03EffTrack->AddEntry(backEffTrackDR03,"TrackIsoDR03 Data MET<30","lf");
  DR03EffTrack->SetFillColor(kWhite);
  //DR03EffTrack->SetTextSize(0.03);
  DR03EffTrack->Draw("SAME");
  c1->Print("Plots/DR03/TrackIsoEffDR03.png");

  float TDR03Sig[150]={0},TDR03Bg[150]={0},TDR03SoverRootB[150]={0},TDR03LogSoverLogB[150]={0};
  for(int i=1;i<151;i++){
    TDR03Sig[i-1]=sigEffTrackDR03->GetBinContent(i);
    TDR03Bg[i-1]=backEffTrackDR03->GetBinContent(i);
    TDR03SoverRootB[i-1]=TDR03Sig[i-1]/sqrt(TDR03Bg[i-1]);
    TDR03LogSoverLogB[i-1]=log(TDR03Sig[i-1])/log(TDR03Bg[i-1]);
  }

  TGraph *bVsSTrackDR03 = new TGraph(150,TDR03Bg,TDR03Sig);
  if(useQCD)bVsSTrackDR03->SetTitle("Integrated Spectra - SignalMC VS QCD - TRACK");
  if(useData)bVsSTrackDR03->SetTitle("Integrated Spectra - SignalMC VS Data MET<30 - TRACK");
  bVsSTrackDR03->GetYaxis()->SetTitle("Signal MC");
  if(useQCD)bVsSTrackDR03->GetXaxis()->SetTitle("Background QCD");
  if(useData)bVsSTrackDR03->GetXaxis()->SetTitle("Background Data MET<30");
  bVsSTrackDR03->SetMarkerColor(kBlue);
  bVsSTrackDR03->SetMarkerStyle(20);
  bVsSTrackDR03->SetMarkerSize(.5);
  bVsSTrackDR03->Draw("AP");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(bVsSTrackDR03,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalVsBackgroundEffTrack.png");

  TH1F *SoverRootBvsTrackIsoDR03 =  sig_trackDR03->Clone();
  pair<int,float> max;max.first=0;max.second=0.;
  for(int i=1;i<151;i++){
    if(TMath::IsNaN(TDR03SoverRootB[i-1]) || !TMath::Finite(TDR03SoverRootB[i-1])) SoverRootBvsTrackIsoDR03->SetBinContent(i,1);
    else{
      SoverRootBvsTrackIsoDR03->SetBinContent(i,TDR03SoverRootB[i-1]);
      if(TDR03SoverRootB[i-1]>max.second){
	max.first=i;max.second=TDR03SoverRootB[i-1];
      }
    }
  }
  cout<<"Track SoverRootB Max:  (Bin,LowEdge,NextlowEdge) ("<<max.first<<","<<SoverRootBvsTrackIsoDR03->GetBinLowEdge(max.first)<<","<<SoverRootBvsTrackIsoDR03->GetBinLowEdge(max.first+1)<<")   Value: "<<max.second<<endl;

  SoverRootBvsTrackIsoDR03->SetTitle("S/#sqrt{B} Vs  TrackIsoDR03");
  SoverRootBvsTrackIsoDR03->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  SoverRootBvsTrackIsoDR03->GetXaxis()->SetTitle("TrackIsoDR03");
  SoverRootBvsTrackIsoDR03->GetYaxis()->SetTitleOffset(0.75);
  SoverRootBvsTrackIsoDR03->GetXaxis()->SetTitleOffset(0.75);
  SoverRootBvsTrackIsoDR03->SetMarkerColor(kBlue);
  SoverRootBvsTrackIsoDR03->SetMarkerSize(.5);
  SoverRootBvsTrackIsoDR03->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(SoverRootBvsTrackIsoDR03,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalEffOverRootBgroundEffVsTrackIso.png");
  
  TH1F *LogSoverLogBvsTrackIsoDR03 =  sig_trackDR03->Clone();
  for(int i=1;i<151;i++){
    //cout<<"bin:"<<i<<"   "<<DR03LogSoverLogB[i-1]<<endl;
    if(TMath::IsNaN(TDR03LogSoverLogB[i-1]) || !TMath::Finite(TDR03LogSoverLogB[i-1]) ) LogSoverLogBvsTrackIsoDR03->SetBinContent(i,1);
    else LogSoverLogBvsTrackIsoDR03->SetBinContent(i,TDR03LogSoverLogB[i-1]);
  }
  LogSoverLogBvsTrackIsoDR03->SetTitle("log(S)/log(B) Vs  TrackIsoDR03");
  LogSoverLogBvsTrackIsoDR03->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  LogSoverLogBvsTrackIsoDR03->GetXaxis()->SetTitle("TrackIsoDR03");
  LogSoverLogBvsTrackIsoDR03->GetYaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsTrackIsoDR03->GetXaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsTrackIsoDR03->SetMarkerColor(kBlue);
  LogSoverLogBvsTrackIsoDR03->SetMarkerSize(.5);
  LogSoverLogBvsTrackIsoDR03->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(LogSoverLogBvsTrackIsoDR03,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/LogSignalEffOverLogBgroundEffVsTrackIso.png");
  
  //----now Comb------
  TH1F *sig_combDR03 = fSig.Get("SignalCombIsoDR03Nminus3");
  TH1F *back_combDR03 = fBack.Get("BgroundCombIsoDR03Nminus3");//data and qcd


  //sig_combDR03->Rebin(2);
  //back_combDR03->Rebin(2);
  sig_combDR03->SetLineColor(kBlue);
  back_combDR03->SetLineColor(kRed);

  TH1F * sig_combDR03New = sig_combDR03->Clone();
  TH1F* back_combDR03New = back_combDR03->Clone();

  sig_combDR03New->Scale(1/sig_combDR03New->Integral());
  back_combDR03New->Scale(1/back_combDR03New->Integral());
  sig_combDR03New->SetTitle("Fraction of Events VS CombIsoDR03");
  sig_combDR03New->GetYaxis()->SetTitle("Fraction of Events");
  sig_combDR03New->GetYaxis()->SetTitleOffset(0.75);
  sig_combDR03New->GetXaxis()->SetTitleOffset(0.75);
  sig_combDR03New->Draw();
  back_combDR03New->Draw("SAME");
  TLegend *DR03Iso = new TLegend(.45,.36,.85,.66,"","brNDC");
  DR03Iso->AddEntry(sig_combDR03,"CombIsoDR03 SignalMC_2000_1015_305","lf");
  if(useQCD)DR03Iso->AddEntry(back_combDR03,"CombIsoDR03 QCD","lf");
  if(useData)DR03Iso->AddEntry(back_combDR03,"CombIsoDR03 Data MET<30","lf");
  DR03Iso->SetFillColor(kWhite);
  // DR03Iso->SetTextSize(0.03);
  DR03Iso->Draw("SAME");
  c1->Print("Plots/DR03/CombDR03Iso.png");


  TH1F *sigEffCombDR03 = sig_combDR03->Clone();
  TH1F *backEffCombDR03 =  back_combDR03->Clone();
  sigEffCombDR03->SetTitle("Integrated Fraction of Events VS CombIsoDR03");
  sigEffCombDR03->GetYaxis()->SetTitle("Integrated Fraction of Events");
  sigEffCombDR03->GetXaxis()->SetTitle("CombIsoDR03");
  sigEffCombDR03->GetYaxis()->SetTitleOffset(0.75);
  sigEffCombDR03->GetXaxis()->SetTitleOffset(0.75);
  sigEffCombDR03->SetLineColor(kBlue);
  backEffCombDR03->SetLineColor(kRed);

  for(int i=1;i<211;i++){
    sigEffCombDR03->SetBinContent(i,sig_combDR03->Integral(0,i)/sig_combDR03->Integral());
    backEffCombDR03->SetBinContent(i,back_combDR03->Integral(0,i)/back_combDR03->Integral());
  }


  sigEffCombDR03->Draw();
  backEffCombDR03->Draw("SAME");
  TLegend *DR03EffComb = new TLegend(.35,.3,.8,.5,"","brNDC");
  DR03EffComb->AddEntry(sigEffCombDR03,"CombIsoDR03 SignalMC_2000_1015_305","lf");
  if(useQCD)DR03EffComb->AddEntry(backEffCombDR03,"CombIsoDR03 QCD","lf");
  if(useData)DR03EffComb->AddEntry(backEffCombDR03,"CombIsoDR03 Data MET<30","lf");
  DR03EffComb->SetFillColor(kWhite);
  //DR03EffComb->SetTextSize(0.03);
  DR03EffComb->Draw("SAME");
  c1->Print("Plots/DR03/CombIsoEffDR03.png");

  float CDR03Sig[210]={0},CDR03Bg[210]={0},CDR03SoverRootB[210]={0},CDR03LogSoverLogB[210]={0};
  for(int i=1;i<211;i++){
    CDR03Sig[i-1]=sigEffCombDR03->GetBinContent(i);
    CDR03Bg[i-1]=backEffCombDR03->GetBinContent(i);
    CDR03SoverRootB[i-1]=CDR03Sig[i-1]/sqrt(CDR03Bg[i-1]);
    CDR03LogSoverLogB[i-1]=log(CDR03Sig[i-1])/log(CDR03Bg[i-1]);
    // cout<<CDR03LogSoverLogB[i-1]<<endl;
  }

  TGraph *bVsSCombDR03 = new TGraph(210,CDR03Bg,CDR03Sig);
  if(useQCD)bVsSCombDR03->SetTitle("Integrated Spectra - SignalMC VS QCD - Combined");
  if(useData)bVsSCombDR03->SetTitle("Integrated Spectra - SignalMC VS Data MET<30 - Combined");
  bVsSCombDR03->GetYaxis()->SetTitle("Signal MC");
  if(useQCD)bVsSCombDR03->GetXaxis()->SetTitle("Background QCD");
  if(useData)bVsSCombDR03->GetXaxis()->SetTitle("Background Data MET<30");
  bVsSCombDR03->SetMarkerColor(kBlue);
  bVsSCombDR03->SetMarkerStyle(20);
  bVsSCombDR03->SetMarkerSize(.5);
  bVsSCombDR03->Draw("AP");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(bVsSCombDR03,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalVsBackgroundEffCombDR03.png");

  TH1F *SoverRootBvsCombIsoDR03 =  sig_combDR03->Clone();
  pair<int,float> max;max.first=0;max.second=0.;
  for(int i=1;i<211;i++){
    if(TMath::IsNaN(CDR03SoverRootB[i-1]) || !TMath::Finite(CDR03SoverRootB[i-1])) SoverRootBvsCombIsoDR03->SetBinContent(i,1);
    else{
      SoverRootBvsCombIsoDR03->SetBinContent(i,CDR03SoverRootB[i-1]);
      if(CDR03SoverRootB[i-1]>max.second){
	max.first=i;max.second=CDR03SoverRootB[i-1];
      }
    }
  }

  cout<<"CombDR03 SoverRootB Max:  (Bin,LowEdge,NextlowEdge) ("<<max.first<<","<<SoverRootBvsCombIsoDR03->GetBinLowEdge(max.first)<<","<<SoverRootBvsCombIsoDR03->GetBinLowEdge(max.first+1)<<")   Value: "<<max.second<<endl;

  SoverRootBvsCombIsoDR03->SetTitle("S/#sqrt{B} Vs  CombIsoDR03");
  SoverRootBvsCombIsoDR03->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  SoverRootBvsCombIsoDR03->GetXaxis()->SetTitle("CombIsoDR03");
  SoverRootBvsCombIsoDR03->GetYaxis()->SetTitleOffset(0.75);
  SoverRootBvsCombIsoDR03->GetXaxis()->SetTitleOffset(0.75);
  SoverRootBvsCombIsoDR03->SetMarkerColor(kBlue);
  SoverRootBvsCombIsoDR03->SetMarkerStyle(20);
  SoverRootBvsCombIsoDR03->SetMarkerSize(.4);
  SoverRootBvsCombIsoDR03->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(SoverRootBvsCombIsoDR03,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalEffOverRootBgroundEffVsCombIsoDR03.png");
  //SoverRootBvsCombIsoDR03->Fit("gaus","","",0,10);
  //c1->Print("Plots/DR03/SignalEffOverRootBgroundEffVsCombIsoDR03_Fit.png");

  fout.cd();SoverRootBvsCombIsoDR03->Write("SignalEffOverRootBgroundEffVsCombIsoDR03_2012");

  TH1F *LogSoverLogBvsCombIsoDR03 =  sig_combDR03->Clone();
  for(int i=1;i<211;i++){
    //cout<<"bin:"<<i<<"   "<<DR03LogSoverLogB[i-1]<<endl;
    if(TMath::IsNaN(CDR03LogSoverLogB[i-1]) || !TMath::Finite(CDR03LogSoverLogB[i-1]) ) LogSoverLogBvsCombIsoDR03->SetBinContent(i,1);
    else LogSoverLogBvsCombIsoDR03->SetBinContent(i,CDR03LogSoverLogB[i-1]);
  }
  LogSoverLogBvsCombIsoDR03->SetTitle("log(S)/log(B) Vs  CombIsoDR03");
  LogSoverLogBvsCombIsoDR03->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  LogSoverLogBvsCombIsoDR03->GetXaxis()->SetTitle("CombIsoDR03");
  LogSoverLogBvsCombIsoDR03->GetYaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsCombIsoDR03->GetXaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsCombIsoDR03->SetMarkerColor(kBlue);
  LogSoverLogBvsCombIsoDR03->SetMarkerSize(.5);
  LogSoverLogBvsCombIsoDR03->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(LogSoverLogBvsCombIsoDR03,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/LogSignalEffOverLogBgroundEffVsCombIsoDR03.png");

  TH1F *sig_combDR04 = fSig.Get("SignalCombIsoDR04Nminus3");
  TH1F *back_combDR04 = fBack.Get("BgroundCombIsoDR04Nminus3");//data and qcd
  
  
  //sig_combDR04->Rebin(2);
  //back_combDR04->Rebin(2);
  sig_combDR04->SetLineColor(kBlue);
  back_combDR04->SetLineColor(kRed);
  
  TH1F * sig_combDR04New = sig_combDR04->Clone();
  TH1F* back_combDR04New = back_combDR04->Clone();
  
  sig_combDR04New->Scale(1/sig_combDR04New->Integral());
  back_combDR04New->Scale(1/back_combDR04New->Integral());
  sig_combDR04New->SetTitle("Fraction of Events VS CombIsoDR04");
  sig_combDR04New->GetYaxis()->SetTitle("Fraction of Events");
  sig_combDR04New->GetYaxis()->SetTitleOffset(0.75);
  sig_combDR04New->GetXaxis()->SetTitleOffset(0.75);
  sig_combDR04New->Draw();
  back_combDR04New->Draw("SAME");
  TLegend *DR04Iso = new TLegend(.45,.36,.85,.66,"","brNDC");
  DR04Iso->AddEntry(sig_combDR04,"CombIsoDR04 SignalMC_2000_1015_305","lf");
  if(useQCD)DR04Iso->AddEntry(back_combDR04,"CombIsoDR04 QCD","lf");
  if(useData)DR04Iso->AddEntry(back_combDR04,"CombIsoDR04 Data MET<30","lf");
  DR04Iso->SetFillColor(kWhite);
  // DR04Iso->SetTextSize(0.03);
  DR04Iso->Draw("SAME");
  c1->Print("Plots/DR03/CombDR04Iso.png");


  TH1F *sigEffCombDR04 = sig_combDR04->Clone();
  TH1F *backEffCombDR04 =  back_combDR04->Clone();
  sigEffCombDR04->SetTitle("Integrated Fraction of Events VS CombIsoDR04");
  sigEffCombDR04->GetYaxis()->SetTitle("Integrated Fraction of Events");
  sigEffCombDR04->GetXaxis()->SetTitle("CombIsoDR04");
  sigEffCombDR04->GetYaxis()->SetTitleOffset(0.75);
  sigEffCombDR04->GetXaxis()->SetTitleOffset(0.75);
  sigEffCombDR04->SetLineColor(kBlue);
  backEffCombDR04->SetLineColor(kRed);

  for(int i=1;i<211;i++){
    sigEffCombDR04->SetBinContent(i,sig_combDR04->Integral(0,i)/sig_combDR04->Integral());
    backEffCombDR04->SetBinContent(i,back_combDR04->Integral(0,i)/back_combDR04->Integral());
  }


  sigEffCombDR04->Draw();
  backEffCombDR04->Draw("SAME");
  TLegend *DR04EffComb = new TLegend(.35,.3,.8,.5,"","brNDC");
  DR04EffComb->AddEntry(sigEffCombDR04,"CombIsoDR04 SignalMC_2000_1015_305","lf");
  if(useQCD)DR04EffComb->AddEntry(backEffCombDR04,"CombIsoDR04 QCD","lf");
  if(useData)DR04EffComb->AddEntry(backEffCombDR04,"CombIsoDR04 Data MET<30","lf");
  DR04EffComb->SetFillColor(kWhite);
  //DR04EffComb->SetTextSize(0.03);
  DR04EffComb->Draw("SAME");
  c1->Print("Plots/DR03/CombIsoEffDR04.png");

  float CDR04Sig[210]={0},CDR04Bg[210]={0},CDR04SoverRootB[210]={0},CDR04LogSoverLogB[210]={0};
  for(int i=1;i<211;i++){
    CDR04Sig[i-1]=sigEffCombDR04->GetBinContent(i);
    CDR04Bg[i-1]=backEffCombDR04->GetBinContent(i);
    CDR04SoverRootB[i-1]=CDR04Sig[i-1]/sqrt(CDR04Bg[i-1]);
    CDR04LogSoverLogB[i-1]=log(CDR04Sig[i-1])/log(CDR04Bg[i-1]);
    // cout<<CDR04LogSoverLogB[i-1]<<endl;
  }

  TGraph *bVsSCombDR04 = new TGraph(210,CDR04Bg,CDR04Sig);
  if(useQCD)bVsSCombDR04->SetTitle("Integrated Spectra - SignalMC VS QCD - Combined");
  if(useData)bVsSCombDR04->SetTitle("Integrated Spectra - SignalMC VS Data MET<30 - Combined");
  bVsSCombDR04->GetYaxis()->SetTitle("Signal MC");
  if(useQCD)bVsSCombDR04->GetXaxis()->SetTitle("Background QCD");
  if(useData)bVsSCombDR04->GetXaxis()->SetTitle("Background Data MET<30");
  bVsSCombDR04->SetMarkerColor(kRed);
  bVsSCombDR04->SetMarkerStyle(20);
  bVsSCombDR04->SetMarkerSize(.5);
  bVsSCombDR04->Draw("AP");
  TLegend *Eff = new TLegend(.55,.22,.8,.42,"","brNDC");
  Eff->AddEntry(bVsSCombDR04,"DR04","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalVsBackgroundEffCombDR04.png");
  bVsSCombDR03->Draw("sameP");
  TLegend *Eff2 = new TLegend(.55,.37,.8,.57,"","brNDC");
  Eff2->AddEntry(bVsSCombDR03,"DR03","p");
  Eff2->SetFillColor(0);
  Eff2->Draw("SAME");
  c1->Print("Plots/DR03/SignalVsBackgroundEffCombBoth.png");

  TH1F *SoverRootBvsCombIsoDR04 =  sig_combDR04->Clone();
  pair<int,float> max;max.first=0;max.second=0.;
  for(int i=1;i<211;i++){
    if(TMath::IsNaN(CDR04SoverRootB[i-1]) || !TMath::Finite(CDR04SoverRootB[i-1])) SoverRootBvsCombIsoDR04->SetBinContent(i,1);
    else{
      SoverRootBvsCombIsoDR04->SetBinContent(i,CDR04SoverRootB[i-1]);
      if(CDR04SoverRootB[i-1]>max.second && i<13){
	max.first=i;max.second=CDR04SoverRootB[i-1];
      }
    }
  }

  cout<<"CombDR04 SoverRootB Max:  (Bin,LowEdge,NextlowEdge) ("<<max.first<<","<<SoverRootBvsCombIsoDR04->GetBinLowEdge(max.first)<<","<<SoverRootBvsCombIsoDR04->GetBinLowEdge(max.first+1)<<")   Value: "<<max.second<<endl;

  SoverRootBvsCombIsoDR04->SetTitle("S/#sqrt{B} Vs  CombIsoDR04");
  SoverRootBvsCombIsoDR03->SetTitle("S/#sqrt{B}");
  SoverRootBvsCombIsoDR04->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  SoverRootBvsCombIsoDR04->GetXaxis()->SetTitle("CombIsoDR04");
  SoverRootBvsCombIsoDR04->GetYaxis()->SetTitleOffset(0.75);
  SoverRootBvsCombIsoDR04->GetXaxis()->SetTitleOffset(0.75);
  SoverRootBvsCombIsoDR04->SetMarkerColor(kBlue);
  SoverRootBvsCombIsoDR04->SetMarkerColor(kRed);
  SoverRootBvsCombIsoDR04->SetMarkerStyle(20);
  SoverRootBvsCombIsoDR04->SetMarkerSize(.4);
  SoverRootBvsCombIsoDR04->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(SoverRootBvsCombIsoDR04,"DR04","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalEffOverRootBgroundEffVsCombIsoDR04.png");
  SoverRootBvsCombIsoDR03->Draw("P");
  SoverRootBvsCombIsoDR04->Draw("PSAME");
  Eff->AddEntry(SoverRootBvsCombIsoDR03,"DR03","p");
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalEffOverRootBgroundEffVsCombIsoBoth.png");
  /*
    TH1F *LogSoverLogBvsCombIsoDR04 =  sig_combDR04->Clone();
    for(int i=1;i<211;i++){
    //cout<<"bin:"<<i<<"   "<<DR04LogSoverLogB[i-1]<<endl;
    if(TMath::IsNaN(CDR04LogSoverLogB[i-1]) || !TMath::Finite(CDR04LogSoverLogB[i-1]) ) LogSoverLogBvsCombIsoDR04->SetBinContent(i,1);
    else LogSoverLogBvsCombIsoDR04->SetBinContent(i,CDR04LogSoverLogB[i-1]);
    }
    LogSoverLogBvsCombIsoDR04->SetTitle("log(S)/log(B) Vs  CombIsoDR04");
    LogSoverLogBvsCombIsoDR04->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
    LogSoverLogBvsCombIsoDR04->GetXaxis()->SetTitle("CombIsoDR04");
    LogSoverLogBvsCombIsoDR04->GetYaxis()->SetTitleOffset(0.75);
    LogSoverLogBvsCombIsoDR04->GetXaxis()->SetTitleOffset(0.75);
    LogSoverLogBvsCombIsoDR04->SetMarkerColor(kBlue);
    LogSoverLogBvsCombIsoDR04->SetMarkerSize(.5);
    LogSoverLogBvsCombIsoDR04->Draw("P");
    TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
    Eff->AddEntry(LogSoverLogBvsCombIsoDR04,"DR04","p");
    Eff->SetFillColor(kWhite);
    Eff->Draw("SAME");
    c1->Print("Plots/DR03/LogSignalEffOverLogBgroundEffVsCombIso.png");
  */

  //Now chargedHadronIso
  TH1F *sig_chargedHadronIso = fSig.Get("chargedHadronIsoNminus3Signal");
  TH1F *back_chargedHadronIso = fBack.Get("chargedHadronIsoNminus3Bground");//data and qcd
  //if(useQCD)TH1F *back_chargedHadronIso = fBack.Get("chargedHadronIsoNminus1Bground");//qcd


  sig_chargedHadronIso->SetLineColor(kBlue);
  back_chargedHadronIso->SetLineColor(kRed);

  TH1F * sig_chargedHadronIsoNew = sig_chargedHadronIso->Clone();
  TH1F* back_chargedHadronIsoNew = back_chargedHadronIso->Clone();

  sig_chargedHadronIsoNew->Scale(1/sig_chargedHadronIsoNew->Integral());
  back_chargedHadronIsoNew->Scale(1/back_chargedHadronIsoNew->Integral());
  sig_chargedHadronIsoNew->SetTitle("Fraction of Events VS chargedHadronIso");
  sig_chargedHadronIsoNew->GetYaxis()->SetTitle("Fraction of Events");
  sig_chargedHadronIsoNew->Draw();
  back_chargedHadronIsoNew->Draw("SAME");
  TLegend *DR03Iso = new TLegend(.45,.36,.85,.66,"","brNDC");
  DR03Iso->AddEntry(sig_chargedHadronIso,"chargedHadronIso SignalMC_2000_1015_305","lf");
  if(useQCD)DR03Iso->AddEntry(back_chargedHadronIso,"chargedHadronIso QCD","lf");
  if(useData)DR03Iso->AddEntry(back_chargedHadronIso,"chargedHadronIso Data MET<30","lf");
  DR03Iso->SetFillColor(kWhite);
  // DR03Iso->SetTextSize(0.03);
  DR03Iso->Draw("SAME");
  c1->Print("Plots/DR03/chargedHadronIso.png");


  TH1F *sigEffchargedHadronIso = sig_chargedHadronIso->Clone();
  TH1F *backEffchargedHadronIso = back_chargedHadronIso->Clone();
  sigEffchargedHadronIso->SetTitle("Integrated Fraction of Events VS chargedHadronIso");
  sigEffchargedHadronIso->GetYaxis()->SetTitle("Integrated Fraction of Events");
  sigEffchargedHadronIso->GetXaxis()->SetTitle("chargedHadronIso");
  sigEffchargedHadronIso->GetYaxis()->SetTitleOffset(0.75);
  sigEffchargedHadronIso->GetXaxis()->SetTitleOffset(0.75);
  sigEffchargedHadronIso->SetLineColor(kBlue);
  backEffchargedHadronIso->SetLineColor(kRed);

  for(int i=1;i<151;i++){
    sigEffchargedHadronIso->SetBinContent(i,sig_chargedHadronIso->Integral(0,i)/sig_chargedHadronIso->Integral());
    backEffchargedHadronIso->SetBinContent(i,back_chargedHadronIso->Integral(0,i)/back_chargedHadronIso->Integral());
  }


  sigEffchargedHadronIso->Draw();
  backEffchargedHadronIso->Draw("SAME");
  TLegend *DR03EffchargedHadron = new TLegend(.35,.3,.8,.5,"","brNDC");
  DR03EffchargedHadron->AddEntry(sigEffchargedHadronIso,"chargedHadronIso SignalMC_2000_1015_305","lf");
  if(useQCD)DR03EffchargedHadron->AddEntry(backEffchargedHadronIso,"chargedHadronIso QCD","lf");
  if(useData)DR03EffchargedHadron->AddEntry(backEffchargedHadronIso,"chargedHadronIso Data MET<30","lf");
  DR03EffchargedHadron->SetFillColor(kWhite);
  //DR03EffchargedHadron->SetTextSize(0.03);
  DR03EffchargedHadron->Draw("SAME");
  c1->Print("Plots/DR03/chargedHadronIsoEffDR03.png");

  float DR03Sig[150]={0},DR03Bg[150]={0},DR03SoverRootB[150]={0},DR03LogSoverLogB[150]={0};
  for(int i=1;i<151;i++){
    DR03Sig[i-1]=sigEffchargedHadronIso->GetBinContent(i);
    DR03Bg[i-1]=backEffchargedHadronIso->GetBinContent(i);
    DR03SoverRootB[i-1]=DR03Sig[i-1]/sqrt(DR03Bg[i-1]);
    DR03LogSoverLogB[i-1]=log(DR03Sig[i-1])/log(DR03Bg[i-1]);
  }

  TGraph *bVsSchargedHadronIso = new TGraph(150,DR03Bg,DR03Sig);
  if(useQCD)bVsSchargedHadronIso->SetTitle("Integrated Spectra - SignalMC VS QCD - chargedHadronIso");
  if(useData)bVsSchargedHadronIso->SetTitle("Integrated Spectra - SignalMC VS Data MET<30 - chargedHadronIso");
  bVsSchargedHadronIso->GetYaxis()->SetTitle("Signal MC");
  if(useQCD)bVsSchargedHadronIso->GetXaxis()->SetTitle("Background QCD");
  if(useData)bVsSchargedHadronIso->GetXaxis()->SetTitle("Background Data MET<30");
  bVsSchargedHadronIso->SetMarkerColor(kBlue);
  bVsSchargedHadronIso->SetMarkerStyle(20);
  bVsSchargedHadronIso->SetMarkerSize(.5);
  bVsSchargedHadronIso->Draw("AP");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(bVsSchargedHadronIso,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalVsBackgroundEffchargedHadron.png");

  TH1F *SoverRootBvschargedHadronIso =  sig_chargedHadronIso->Clone();
  pair<int,float> max;max.first=0;max.second=0.;
  for(int i=1;i<151;i++){
    if(TMath::IsNaN(DR03SoverRootB[i-1]) || !TMath::Finite(DR03SoverRootB[i-1]) ) SoverRootBvschargedHadronIso->SetBinContent(i,1);
    else {
      SoverRootBvschargedHadronIso->SetBinContent(i,DR03SoverRootB[i-1]);
      if(DR03SoverRootB[i-1]>max.second){
	max.first=i;max.second=DR03SoverRootB[i-1];
      }
    }
  }
  cout<<"chargedHadron SoverRootB Max:  (Bin,LowEdge,NextlowEdge) ("<<max.first<<","<<SoverRootBvschargedHadronIso->GetBinLowEdge(max.first)<<","<<SoverRootBvschargedHadronIso->GetBinLowEdge(max.first+1)<<")   Value: "<<max.second<<endl;

  SoverRootBvschargedHadronIso->SetTitle("S/#sqrt{B} Vs  chargedHadronIso");
  SoverRootBvschargedHadronIso->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  SoverRootBvschargedHadronIso->GetXaxis()->SetTitle("chargedHadronIso"); 
  SoverRootBvschargedHadronIso->GetYaxis()->SetTitleOffset(0.75);
  SoverRootBvschargedHadronIso->GetXaxis()->SetTitleOffset(0.75);
  SoverRootBvschargedHadronIso->SetMarkerColor(kBlue);
  SoverRootBvschargedHadronIso->SetMarkerSize(.5);
  SoverRootBvschargedHadronIso->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(SoverRootBvschargedHadronIso,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalEffOverRootBgroundEffVschargedHadronIso.png");
  
  TH1F *LogSoverLogBvschargedHadronIso =  sig_chargedHadronIso->Clone();
  for(int i=1;i<151;i++){
    //cout<<"bin:"<<i<<"   "<<DR03LogSoverLogB[i-1]<<endl;
    if(TMath::IsNaN(DR03LogSoverLogB[i-1]) || !TMath::Finite(DR03LogSoverLogB[i-1]) ) LogSoverLogBvschargedHadronIso->SetBinContent(i,1);
    else LogSoverLogBvschargedHadronIso->SetBinContent(i,DR03LogSoverLogB[i-1]);
  }
  LogSoverLogBvschargedHadronIso->SetTitle("log(S)/log(B) Vs  chargedHadronIso");
  LogSoverLogBvschargedHadronIso->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  LogSoverLogBvschargedHadronIso->GetXaxis()->SetTitle("chargedHadronIso");
  LogSoverLogBvschargedHadronIso->GetYaxis()->SetTitleOffset(0.75);
  LogSoverLogBvschargedHadronIso->GetXaxis()->SetTitleOffset(0.75);
  LogSoverLogBvschargedHadronIso->SetMarkerColor(kBlue);
  LogSoverLogBvschargedHadronIso->SetMarkerSize(.5);
  LogSoverLogBvschargedHadronIso->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(LogSoverLogBvschargedHadronIso,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/LogSignalEffOverLogBgroundEffVschargedHadronIso.png");

  //Now chargedHadronIsoDeposit
  TH1F *sig_chargedHadronIsoDeposit = fSig.Get("chargedHadronIsoDepositNminus3Signal");
  TH1F *back_chargedHadronIsoDeposit = fBack.Get("chargedHadronIsoDepositNminus3Bground");//data and qcd
  //if(useQCD)TH1F *back_chargedHadronIsoDeposit = fBack.Get("chargedHadronIsoDepositNminus1Bground");//qcd


  sig_chargedHadronIsoDeposit->SetLineColor(kBlue);
  back_chargedHadronIsoDeposit->SetLineColor(kRed);

  TH1F * sig_chargedHadronIsoDepositNew = sig_chargedHadronIsoDeposit->Clone();
  TH1F* back_chargedHadronIsoDepositNew = back_chargedHadronIsoDeposit->Clone();

  sig_chargedHadronIsoDepositNew->Scale(1/sig_chargedHadronIsoDepositNew->Integral());
  back_chargedHadronIsoDepositNew->Scale(1/back_chargedHadronIsoDepositNew->Integral());
  sig_chargedHadronIsoDepositNew->SetTitle("Fraction of Events VS chargedHadronIsoDeposit");
  sig_chargedHadronIsoDepositNew->GetYaxis()->SetTitle("Fraction of Events");
  sig_chargedHadronIsoDepositNew->Draw();
  back_chargedHadronIsoDepositNew->Draw("SAME");
  TLegend *chargedHadronIsoDeposit = new TLegend(.45,.36,.85,.66,"","brNDC");
  chargedHadronIsoDeposit->AddEntry(sig_chargedHadronIsoDeposit,"chargedHadronIsoDeposit SignalMC_2000_1015_305","lf");
  if(useQCD)chargedHadronIsoDeposit->AddEntry(back_chargedHadronIsoDeposit,"chargedHadronIsoDeposit QCD","lf");
  if(useData)chargedHadronIsoDeposit->AddEntry(back_chargedHadronIsoDeposit,"chargedHadronIsoDeposit Data MET<30","lf");
  chargedHadronIsoDeposit->SetFillColor(kWhite);
  // chargedHadronIsoDeposit->SetTextSize(0.03);
  chargedHadronIsoDeposit->Draw("SAME");
  c1->Print("Plots/DR03/chargedHadronIsoDeposit.png");


  TH1F *sigEffchargedHadronIsoDeposit = sig_chargedHadronIsoDeposit->Clone();
  TH1F *backEffchargedHadronIsoDeposit = back_chargedHadronIsoDeposit->Clone();
  sigEffchargedHadronIsoDeposit->SetTitle("Integrated Fraction of Events VS chargedHadronIsoDeposit");
  sigEffchargedHadronIsoDeposit->GetYaxis()->SetTitle("Integrated Fraction of Events");
  sigEffchargedHadronIsoDeposit->GetXaxis()->SetTitle("chargedHadronIsoDeposit");
  sigEffchargedHadronIsoDeposit->GetYaxis()->SetTitleOffset(0.75);
  sigEffchargedHadronIsoDeposit->GetXaxis()->SetTitleOffset(0.75);
  sigEffchargedHadronIsoDeposit->SetLineColor(kBlue);
  backEffchargedHadronIsoDeposit->SetLineColor(kRed);

  for(int i=1;i<151;i++){
    sigEffchargedHadronIsoDeposit->SetBinContent(i,sig_chargedHadronIsoDeposit->Integral(0,i)/sig_chargedHadronIsoDeposit->Integral());
    backEffchargedHadronIsoDeposit->SetBinContent(i,back_chargedHadronIsoDeposit->Integral(0,i)/back_chargedHadronIsoDeposit->Integral());
  }


  sigEffchargedHadronIsoDeposit->Draw();
  backEffchargedHadronIsoDeposit->Draw("SAME");
  TLegend *DR03EffchargedHadronDeposit = new TLegend(.35,.3,.8,.5,"","brNDC");
  DR03EffchargedHadronDeposit->AddEntry(sigEffchargedHadronIsoDeposit,"chargedHadronIsoDeposit SignalMC_2000_1015_305","lf");
  if(useQCD)DR03EffchargedHadronDeposit->AddEntry(backEffchargedHadronIsoDeposit,"chargedHadronIsoDeposit QCD","lf");
  if(useData)DR03EffchargedHadronDeposit->AddEntry(backEffchargedHadronIsoDeposit,"chargedHadronIsoDeposit Data MET<30","lf");
  DR03EffchargedHadronDeposit->SetFillColor(kWhite);
  //DR03EffchargedHadronDeposit->SetTextSize(0.03);
  DR03EffchargedHadronDeposit->Draw("SAME");
  c1->Print("Plots/DR03/chargedHadronIsoDepositEffDR03.png");

  float DR03Sig[150]={0},DR03Bg[150]={0},DR03SoverRootB[150]={0},DR03LogSoverLogB[150]={0};
  for(int i=1;i<151;i++){
    DR03Sig[i-1]=sigEffchargedHadronIsoDeposit->GetBinContent(i);
    DR03Bg[i-1]=backEffchargedHadronIsoDeposit->GetBinContent(i);
    DR03SoverRootB[i-1]=DR03Sig[i-1]/sqrt(DR03Bg[i-1]);
    DR03LogSoverLogB[i-1]=log(DR03Sig[i-1])/log(DR03Bg[i-1]);
  }

  TGraph *bVsSchargedHadronIsoDeposit = new TGraph(150,DR03Bg,DR03Sig);
  bVsSchargedHadronIsoDeposit->SetTitle("Integrated Spectra - SignalMC VS QCD - chargedHadronIsoDeposit");
  bVsSchargedHadronIsoDeposit->GetYaxis()->SetTitle("Signal MC");
  bVsSchargedHadronIsoDeposit->GetXaxis()->SetTitle("Background Data MET<30");
  bVsSchargedHadronIsoDeposit->SetMarkerColor(kBlue);
  bVsSchargedHadronIsoDeposit->SetMarkerStyle(20);
  bVsSchargedHadronIsoDeposit->SetMarkerSize(.5);
  bVsSchargedHadronIsoDeposit->Draw("AP");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(bVsSchargedHadronIsoDeposit,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalVsBackgroundEffchargedHadronIsoDeposit.png");

  TH1F *SoverRootBvschargedHadronIsoDeposit =  sig_chargedHadronIsoDeposit->Clone();
  pair<int,float> max;max.first=0;max.second=0.;
  for(int i=1;i<151;i++){
    if(TMath::IsNaN(DR03SoverRootB[i-1]) || !TMath::Finite(DR03SoverRootB[i-1]) ) SoverRootBvschargedHadronIsoDeposit->SetBinContent(i,1);
    else {
      SoverRootBvschargedHadronIsoDeposit->SetBinContent(i,DR03SoverRootB[i-1]);
      if(DR03SoverRootB[i-1]>max.second){
	max.first=i;max.second=DR03SoverRootB[i-1];
      }
    }
  }
  cout<<"chargedHadronIsoDeposit SoverRootB Max:  (Bin,LowEdge,NextlowEdge) ("<<max.first<<","<<SoverRootBvschargedHadronIsoDeposit->GetBinLowEdge(max.first)<<","<<SoverRootBvschargedHadronIsoDeposit->GetBinLowEdge(max.first+1)<<")   Value: "<<max.second<<endl;

  SoverRootBvschargedHadronIsoDeposit->SetTitle("S/#sqrt{B} Vs  chargedHadronIsoDeposit");
  SoverRootBvschargedHadronIsoDeposit->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  SoverRootBvschargedHadronIsoDeposit->GetXaxis()->SetTitle("chargedHadronIsoDeposit"); 
  SoverRootBvschargedHadronIsoDeposit->GetYaxis()->SetTitleOffset(0.75);
  SoverRootBvschargedHadronIsoDeposit->GetXaxis()->SetTitleOffset(0.75);
  SoverRootBvschargedHadronIsoDeposit->SetMarkerColor(kBlue);
  SoverRootBvschargedHadronIsoDeposit->SetMarkerSize(.5);
  SoverRootBvschargedHadronIsoDeposit->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(SoverRootBvschargedHadronIsoDeposit,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalEffOverRootBgroundEffVschargedHadronIsoDeposit.png");
  
  TH1F *LogSoverLogBvschargedHadronIsoDeposit =  sig_chargedHadronIsoDeposit->Clone();
  for(int i=1;i<151;i++){
    //cout<<"bin:"<<i<<"   "<<DR03LogSoverLogB[i-1]<<endl;
    if(TMath::IsNaN(DR03LogSoverLogB[i-1]) || !TMath::Finite(DR03LogSoverLogB[i-1]) ) LogSoverLogBvschargedHadronIsoDeposit->SetBinContent(i,1);
    else LogSoverLogBvschargedHadronIsoDeposit->SetBinContent(i,DR03LogSoverLogB[i-1]);
  }
  LogSoverLogBvschargedHadronIsoDeposit->SetTitle("log(S)/log(B) Vs  chargedHadronIsoDeposit");
  LogSoverLogBvschargedHadronIsoDeposit->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  LogSoverLogBvschargedHadronIsoDeposit->GetXaxis()->SetTitle("chargedHadronIsoDeposit");
  LogSoverLogBvschargedHadronIsoDeposit->GetYaxis()->SetTitleOffset(0.75);
  LogSoverLogBvschargedHadronIsoDeposit->GetXaxis()->SetTitleOffset(0.75);
  LogSoverLogBvschargedHadronIsoDeposit->SetMarkerColor(kBlue);
  LogSoverLogBvschargedHadronIsoDeposit->SetMarkerSize(.5);
  LogSoverLogBvschargedHadronIsoDeposit->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(LogSoverLogBvschargedHadronIsoDeposit,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/LogSignalEffOverLogBgroundEffVschargedHadronIsoDeposit.png");


  //Now neutralHadronIso
  TH1F *sig_neutralHadronIso = fSig.Get("neutralHadronIsoNminus3Signal");
  TH1F *back_neutralHadronIso = fBack.Get("neutralHadronIsoNminus3Bground");//data and qcd
  //if(useQCD)TH1F *back_neutralHadronIso = fBack.Get("neutralHadronIsoNminus1Bground");//qcd


  sig_neutralHadronIso->SetLineColor(kBlue);
  back_neutralHadronIso->SetLineColor(kRed);

  TH1F * sig_neutralHadronIsoNew = sig_neutralHadronIso->Clone();
  TH1F* back_neutralHadronIsoNew = back_neutralHadronIso->Clone();

  sig_neutralHadronIsoNew->Scale(1/sig_neutralHadronIsoNew->Integral());
  back_neutralHadronIsoNew->Scale(1/back_neutralHadronIsoNew->Integral());
  sig_neutralHadronIsoNew->SetTitle("Fraction of Events VS neutralHadronIso");
  sig_neutralHadronIsoNew->GetYaxis()->SetTitle("Fraction of Events");
  sig_neutralHadronIsoNew->Draw();
  back_neutralHadronIsoNew->Draw("SAME");
  TLegend *DR03Iso = new TLegend(.45,.36,.85,.66,"","brNDC");
  DR03Iso->AddEntry(sig_neutralHadronIso,"neutralHadronIso SignalMC_2000_1015_305","lf");
  DR03Iso->AddEntry(back_neutralHadronIso,"neutralHadronIso Data MET<30","lf");
  DR03Iso->SetFillColor(kWhite);
  // DR03Iso->SetTextSize(0.03);
  DR03Iso->Draw("SAME");
  c1->Print("Plots/DR03/neutralHadronIso.png");


  TH1F *sigEffneutralHadronIso = sig_neutralHadronIso->Clone();
  TH1F *backEffneutralHadronIso = back_neutralHadronIso->Clone();
  sigEffneutralHadronIso->SetTitle("Integrated Fraction of Events VS neutralHadronIso");
  sigEffneutralHadronIso->GetYaxis()->SetTitle("Integrated Fraction of Events");
  sigEffneutralHadronIso->GetXaxis()->SetTitle("neutralHadronIso");
  sigEffneutralHadronIso->GetYaxis()->SetTitleOffset(0.75);
  sigEffneutralHadronIso->GetXaxis()->SetTitleOffset(0.75);
  sigEffneutralHadronIso->SetLineColor(kBlue);
  backEffneutralHadronIso->SetLineColor(kRed);

  for(int i=1;i<151;i++){
    sigEffneutralHadronIso->SetBinContent(i,sig_neutralHadronIso->Integral(0,i)/sig_neutralHadronIso->Integral());
    backEffneutralHadronIso->SetBinContent(i,back_neutralHadronIso->Integral(0,i)/back_neutralHadronIso->Integral());
  }


  sigEffneutralHadronIso->Draw();
  backEffneutralHadronIso->Draw("SAME");
  TLegend *DR03EffneutralHadron = new TLegend(.35,.3,.8,.5,"","brNDC");
  DR03EffneutralHadron->AddEntry(sigEffneutralHadronIso,"neutralHadronIso SignalMC_2000_1015_305","lf");
  DR03EffneutralHadron->AddEntry(backEffneutralHadronIso,"neutralHadronIso Data MET<30","lf");
  DR03EffneutralHadron->SetFillColor(kWhite);
  //DR03EffneutralHadron->SetTextSize(0.03);
  DR03EffneutralHadron->Draw("SAME");
  c1->Print("Plots/DR03/neutralHadronIsoEffDR03.png");

  float DR03Sig[150]={0},DR03Bg[150]={0},DR03SoverRootB[150]={0},DR03LogSoverLogB[150]={0};
  for(int i=1;i<151;i++){
    DR03Sig[i-1]=sigEffneutralHadronIso->GetBinContent(i);
    DR03Bg[i-1]=backEffneutralHadronIso->GetBinContent(i);
    DR03SoverRootB[i-1]=DR03Sig[i-1]/sqrt(DR03Bg[i-1]);
    DR03LogSoverLogB[i-1]=log(DR03Sig[i-1])/log(DR03Bg[i-1]);
  }

  TGraph *bVsSneutralHadronIso = new TGraph(150,DR03Bg,DR03Sig);
  bVsSneutralHadronIso->SetTitle("Integrated Spectra - SignalMC VS QCD - neutralHadronIso");
  bVsSneutralHadronIso->GetYaxis()->SetTitle("Signal MC");
  bVsSneutralHadronIso->GetXaxis()->SetTitle("Background Data MET<30");
  bVsSneutralHadronIso->SetMarkerColor(kBlue);
  bVsSneutralHadronIso->SetMarkerStyle(20);
  bVsSneutralHadronIso->SetMarkerSize(.5);
  bVsSneutralHadronIso->Draw("AP");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(bVsSneutralHadronIso,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalVsBackgroundEffneutralHadron.png");

  TH1F *SoverRootBvsneutralHadronIso =  sig_neutralHadronIso->Clone();
  pair<int,float> max;max.first=0;max.second=0.;
  for(int i=1;i<151;i++){
    if(TMath::IsNaN(DR03SoverRootB[i-1]) || !TMath::Finite(DR03SoverRootB[i-1]) ) SoverRootBvsneutralHadronIso->SetBinContent(i,1);
    else {
      SoverRootBvsneutralHadronIso->SetBinContent(i,DR03SoverRootB[i-1]);
      if(DR03SoverRootB[i-1]>max.second){
	max.first=i;max.second=DR03SoverRootB[i-1];
      }
    }
  }
  cout<<"neutralHadron SoverRootB Max:  (Bin,LowEdge,NextlowEdge) ("<<max.first<<","<<SoverRootBvsneutralHadronIso->GetBinLowEdge(max.first)<<","<<SoverRootBvsneutralHadronIso->GetBinLowEdge(max.first+1)<<")   Value: "<<max.second<<endl;

  SoverRootBvsneutralHadronIso->SetTitle("S/#sqrt{B} Vs  neutralHadronIso");
  SoverRootBvsneutralHadronIso->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  SoverRootBvsneutralHadronIso->GetXaxis()->SetTitle("neutralHadronIso"); 
  SoverRootBvsneutralHadronIso->GetYaxis()->SetTitleOffset(0.75);
  SoverRootBvsneutralHadronIso->GetXaxis()->SetTitleOffset(0.75);
  SoverRootBvsneutralHadronIso->SetMarkerColor(kBlue);
  SoverRootBvsneutralHadronIso->SetMarkerSize(.5);
  SoverRootBvsneutralHadronIso->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(SoverRootBvsneutralHadronIso,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalEffOverRootBgroundEffVsneutralHadronIso.png");
  
  TH1F *LogSoverLogBvsneutralHadronIso =  sig_neutralHadronIso->Clone();
  for(int i=1;i<151;i++){
    //cout<<"bin:"<<i<<"   "<<DR03LogSoverLogB[i-1]<<endl;
    if(TMath::IsNaN(DR03LogSoverLogB[i-1]) || !TMath::Finite(DR03LogSoverLogB[i-1]) ) LogSoverLogBvsneutralHadronIso->SetBinContent(i,1);
    else LogSoverLogBvsneutralHadronIso->SetBinContent(i,DR03LogSoverLogB[i-1]);
  }
  LogSoverLogBvsneutralHadronIso->SetTitle("log(S)/log(B) Vs  neutralHadronIso");
  LogSoverLogBvsneutralHadronIso->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  LogSoverLogBvsneutralHadronIso->GetXaxis()->SetTitle("neutralHadronIso");
  LogSoverLogBvsneutralHadronIso->GetYaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsneutralHadronIso->GetXaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsneutralHadronIso->SetMarkerColor(kBlue);
  LogSoverLogBvsneutralHadronIso->SetMarkerSize(.5);
  LogSoverLogBvsneutralHadronIso->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(LogSoverLogBvsneutralHadronIso,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/LogSignalEffOverLogBgroundEffVsneutralHadronIso.png");

  //Now neutralHadronIsoDeposit
  TH1F *sig_neutralHadronIsoDeposit = fSig.Get("neutralHadronIsoDepositNminus3Signal");
  TH1F *back_neutralHadronIsoDeposit = fBack.Get("neutralHadronIsoDepositNminus3Bground");//data and qcd
  //if(useQCD)TH1F *back_neutralHadronIsoDeposit = fBack.Get("neutralHadronIsoDepositNminus1Bground");//qcd


  sig_neutralHadronIsoDeposit->SetLineColor(kBlue);
  back_neutralHadronIsoDeposit->SetLineColor(kRed);

  TH1F * sig_neutralHadronIsoDepositNew = sig_neutralHadronIsoDeposit->Clone();
  TH1F* back_neutralHadronIsoDepositNew = back_neutralHadronIsoDeposit->Clone();

  sig_neutralHadronIsoDepositNew->Scale(1/sig_neutralHadronIsoDepositNew->Integral());
  back_neutralHadronIsoDepositNew->Scale(1/back_neutralHadronIsoDepositNew->Integral());
  sig_neutralHadronIsoDepositNew->SetTitle("Fraction of Events VS neutralHadronIsoDeposit");
  sig_neutralHadronIsoDepositNew->GetYaxis()->SetTitle("Fraction of Events");
  sig_neutralHadronIsoDepositNew->Draw();
  back_neutralHadronIsoDepositNew->Draw("SAME");
  TLegend *neutralHadronIsoDeposit = new TLegend(.45,.36,.85,.66,"","brNDC");
  neutralHadronIsoDeposit->AddEntry(sig_neutralHadronIsoDeposit,"neutralHadronIsoDeposit SignalMC_2000_1015_305","lf");
  neutralHadronIsoDeposit->AddEntry(back_neutralHadronIsoDeposit,"neutralHadronIsoDeposit Data MET<30","lf");
  neutralHadronIsoDeposit->SetFillColor(kWhite);
  // neutralHadronIsoDeposit->SetTextSize(0.03);
  neutralHadronIsoDeposit->Draw("SAME");
  c1->Print("Plots/DR03/neutralHadronIsoDeposit.png");


  TH1F *sigEffneutralHadronIsoDeposit = sig_neutralHadronIsoDeposit->Clone();
  TH1F *backEffneutralHadronIsoDeposit = back_neutralHadronIsoDeposit->Clone();
  sigEffneutralHadronIsoDeposit->SetTitle("Integrated Fraction of Events VS neutralHadronIsoDeposit");
  sigEffneutralHadronIsoDeposit->GetYaxis()->SetTitle("Integrated Fraction of Events");
  sigEffneutralHadronIsoDeposit->GetXaxis()->SetTitle("neutralHadronIsoDeposit");
  sigEffneutralHadronIsoDeposit->GetYaxis()->SetTitleOffset(0.75);
  sigEffneutralHadronIsoDeposit->GetXaxis()->SetTitleOffset(0.75);
  sigEffneutralHadronIsoDeposit->SetLineColor(kBlue);
  backEffneutralHadronIsoDeposit->SetLineColor(kRed);

  for(int i=1;i<151;i++){
    sigEffneutralHadronIsoDeposit->SetBinContent(i,sig_neutralHadronIsoDeposit->Integral(0,i)/sig_neutralHadronIsoDeposit->Integral());
    backEffneutralHadronIsoDeposit->SetBinContent(i,back_neutralHadronIsoDeposit->Integral(0,i)/back_neutralHadronIsoDeposit->Integral());
  }


  sigEffneutralHadronIsoDeposit->Draw();
  backEffneutralHadronIsoDeposit->Draw("SAME");
  TLegend *DR03EffneutralHadronDeposit = new TLegend(.35,.3,.8,.5,"","brNDC");
  DR03EffneutralHadronDeposit->AddEntry(sigEffneutralHadronIsoDeposit,"neutralHadronIsoDeposit SignalMC_2000_1015_305","lf");
  DR03EffneutralHadronDeposit->AddEntry(backEffneutralHadronIsoDeposit,"neutralHadronIsoDeposit Data MET<30","lf");
  DR03EffneutralHadronDeposit->SetFillColor(kWhite);
  //DR03EffneutralHadronDeposit->SetTextSize(0.03);
  DR03EffneutralHadronDeposit->Draw("SAME");
  c1->Print("Plots/DR03/neutralHadronIsoDepositEffDR03.png");

  float DR03Sig[150]={0},DR03Bg[150]={0},DR03SoverRootB[150]={0},DR03LogSoverLogB[150]={0};
  for(int i=1;i<151;i++){
    DR03Sig[i-1]=sigEffneutralHadronIsoDeposit->GetBinContent(i);
    DR03Bg[i-1]=backEffneutralHadronIsoDeposit->GetBinContent(i);
    DR03SoverRootB[i-1]=DR03Sig[i-1]/sqrt(DR03Bg[i-1]);
    DR03LogSoverLogB[i-1]=log(DR03Sig[i-1])/log(DR03Bg[i-1]);
  }

  TGraph *bVsSneutralHadronIsoDeposit = new TGraph(150,DR03Bg,DR03Sig);
  bVsSneutralHadronIsoDeposit->SetTitle("Integrated Spectra - SignalMC VS QCD - neutralHadronIsoDeposit");
  bVsSneutralHadronIsoDeposit->GetYaxis()->SetTitle("Signal MC");
  bVsSneutralHadronIsoDeposit->GetXaxis()->SetTitle("Background Data MET<30");
  bVsSneutralHadronIsoDeposit->SetMarkerColor(kBlue);
  bVsSneutralHadronIsoDeposit->SetMarkerStyle(20);
  bVsSneutralHadronIsoDeposit->SetMarkerSize(.5);
  bVsSneutralHadronIsoDeposit->Draw("AP");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(bVsSneutralHadronIsoDeposit,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalVsBackgroundEffneutralHadronIsoDeposit.png");

  TH1F *SoverRootBvsneutralHadronIsoDeposit =  sig_neutralHadronIsoDeposit->Clone();
  pair<int,float> max;max.first=0;max.second=0.;
  for(int i=1;i<151;i++){
    if(TMath::IsNaN(DR03SoverRootB[i-1]) || !TMath::Finite(DR03SoverRootB[i-1]) ) SoverRootBvsneutralHadronIsoDeposit->SetBinContent(i,1);
    else {
      SoverRootBvsneutralHadronIsoDeposit->SetBinContent(i,DR03SoverRootB[i-1]);
      if(DR03SoverRootB[i-1]>max.second){
	max.first=i;max.second=DR03SoverRootB[i-1];
      }
    }
  }
  cout<<"neutralHadronIsoDeposit SoverRootB Max:  (Bin,LowEdge,NextlowEdge) ("<<max.first<<","<<SoverRootBvsneutralHadronIsoDeposit->GetBinLowEdge(max.first)<<","<<SoverRootBvsneutralHadronIsoDeposit->GetBinLowEdge(max.first+1)<<")   Value: "<<max.second<<endl;

  SoverRootBvsneutralHadronIsoDeposit->SetTitle("S/#sqrt{B} Vs  neutralHadronIsoDeposit");
  SoverRootBvsneutralHadronIsoDeposit->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  SoverRootBvsneutralHadronIsoDeposit->GetXaxis()->SetTitle("neutralHadronIsoDeposit"); 
  SoverRootBvsneutralHadronIsoDeposit->GetYaxis()->SetTitleOffset(0.75);
  SoverRootBvsneutralHadronIsoDeposit->GetXaxis()->SetTitleOffset(0.75);
  SoverRootBvsneutralHadronIsoDeposit->SetMarkerColor(kBlue);
  SoverRootBvsneutralHadronIsoDeposit->SetMarkerSize(.5);
  SoverRootBvsneutralHadronIsoDeposit->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(SoverRootBvsneutralHadronIsoDeposit,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalEffOverRootBgroundEffVsneutralHadronIsoDeposit.png");
  
  TH1F *LogSoverLogBvsneutralHadronIsoDeposit =  sig_neutralHadronIsoDeposit->Clone();
  for(int i=1;i<151;i++){
    //cout<<"bin:"<<i<<"   "<<DR03LogSoverLogB[i-1]<<endl;
    if(TMath::IsNaN(DR03LogSoverLogB[i-1]) || !TMath::Finite(DR03LogSoverLogB[i-1]) ) LogSoverLogBvsneutralHadronIsoDeposit->SetBinContent(i,1);
    else LogSoverLogBvsneutralHadronIsoDeposit->SetBinContent(i,DR03LogSoverLogB[i-1]);
  }
  LogSoverLogBvsneutralHadronIsoDeposit->SetTitle("log(S)/log(B) Vs  neutralHadronIsoDeposit");
  LogSoverLogBvsneutralHadronIsoDeposit->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  LogSoverLogBvsneutralHadronIsoDeposit->GetXaxis()->SetTitle("neutralHadronIsoDeposit");
  LogSoverLogBvsneutralHadronIsoDeposit->GetYaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsneutralHadronIsoDeposit->GetXaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsneutralHadronIsoDeposit->SetMarkerColor(kBlue);
  LogSoverLogBvsneutralHadronIsoDeposit->SetMarkerSize(.5);
  LogSoverLogBvsneutralHadronIsoDeposit->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(LogSoverLogBvsneutralHadronIsoDeposit,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/LogSignalEffOverLogBgroundEffVsneutralHadronIsoDeposit.png");


  //Now photonIso
  TH1F *sig_photonIso = fSig.Get("photonIsoNminus3Signal");
  TH1F *back_photonIso = fBack.Get("photonIsoNminus3Bground");//data and qcd
  //if(useQCD)TH1F *back_photonIso = fBack.Get("photonIsoNminus1Bground");//qcd


  sig_photonIso->SetLineColor(kBlue);
  back_photonIso->SetLineColor(kRed);

  TH1F * sig_photonIsoNew = sig_photonIso->Clone();
  TH1F* back_photonIsoNew = back_photonIso->Clone();

  sig_photonIsoNew->Scale(1/sig_photonIsoNew->Integral());
  back_photonIsoNew->Scale(1/back_photonIsoNew->Integral());
  sig_photonIsoNew->SetTitle("Fraction of Events VS photonIso");
  sig_photonIsoNew->GetYaxis()->SetTitle("Fraction of Events");
  sig_photonIsoNew->Draw();
  back_photonIsoNew->Draw("SAME");
  TLegend *DR03Iso = new TLegend(.45,.36,.85,.66,"","brNDC");
  DR03Iso->AddEntry(sig_photonIso,"photonIso SignalMC_2000_1015_305","lf");
  DR03Iso->AddEntry(back_photonIso,"photonIso Data MET<30","lf");
  DR03Iso->SetFillColor(kWhite);
  // DR03Iso->SetTextSize(0.03);
  DR03Iso->Draw("SAME");
  c1->Print("Plots/DR03/photonIso.png");


  TH1F *sigEffphotonIso = sig_photonIso->Clone();
  TH1F *backEffphotonIso = back_photonIso->Clone();
  sigEffphotonIso->SetTitle("Integrated Fraction of Events VS photonIso");
  sigEffphotonIso->GetYaxis()->SetTitle("Integrated Fraction of Events");
  sigEffphotonIso->GetXaxis()->SetTitle("photonIso");
  sigEffphotonIso->GetYaxis()->SetTitleOffset(0.75);
  sigEffphotonIso->GetXaxis()->SetTitleOffset(0.75);
  sigEffphotonIso->SetLineColor(kBlue);
  backEffphotonIso->SetLineColor(kRed);

  for(int i=1;i<151;i++){
    sigEffphotonIso->SetBinContent(i,sig_photonIso->Integral(0,i)/sig_photonIso->Integral());
    backEffphotonIso->SetBinContent(i,back_photonIso->Integral(0,i)/back_photonIso->Integral());
  }


  sigEffphotonIso->Draw();
  backEffphotonIso->Draw("SAME");
  TLegend *DR03Effphoton = new TLegend(.35,.3,.8,.5,"","brNDC");
  DR03Effphoton->AddEntry(sigEffphotonIso,"photonIso SignalMC_2000_1015_305","lf");
  DR03Effphoton->AddEntry(backEffphotonIso,"photonIso Data MET<30","lf");
  DR03Effphoton->SetFillColor(kWhite);
  //DR03Effphoton->SetTextSize(0.03);
  DR03Effphoton->Draw("SAME");
  c1->Print("Plots/DR03/photonIsoEffDR03.png");

  float DR03Sig[150]={0},DR03Bg[150]={0},DR03SoverRootB[150]={0},DR03LogSoverLogB[150]={0};
  for(int i=1;i<151;i++){
    DR03Sig[i-1]=sigEffphotonIso->GetBinContent(i);
    DR03Bg[i-1]=backEffphotonIso->GetBinContent(i);
    DR03SoverRootB[i-1]=DR03Sig[i-1]/sqrt(DR03Bg[i-1]);
    DR03LogSoverLogB[i-1]=log(DR03Sig[i-1])/log(DR03Bg[i-1]);
  }

  TGraph *bVsSphotonIso = new TGraph(150,DR03Bg,DR03Sig);
  bVsSphotonIso->SetTitle("Integrated Spectra - SignalMC VS QCD - photonIso");
  bVsSphotonIso->GetYaxis()->SetTitle("Signal MC");
  bVsSphotonIso->GetXaxis()->SetTitle("Background Data MET<30");
  bVsSphotonIso->SetMarkerColor(kBlue);
  bVsSphotonIso->SetMarkerStyle(20);
  bVsSphotonIso->SetMarkerSize(.5);
  bVsSphotonIso->Draw("AP");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(bVsSphotonIso,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalVsBackgroundEffphoton.png");

  TH1F *SoverRootBvsphotonIso =  sig_photonIso->Clone();
  pair<int,float> max;max.first=0;max.second=0.;
  for(int i=1;i<151;i++){
    if(TMath::IsNaN(DR03SoverRootB[i-1]) || !TMath::Finite(DR03SoverRootB[i-1]) ) SoverRootBvsphotonIso->SetBinContent(i,1);
    else {
      SoverRootBvsphotonIso->SetBinContent(i,DR03SoverRootB[i-1]);
      if(DR03SoverRootB[i-1]>max.second){
	max.first=i;max.second=DR03SoverRootB[i-1];
      }
    }
  }
  cout<<"photon SoverRootB Max:  (Bin,LowEdge,NextlowEdge) ("<<max.first<<","<<SoverRootBvsphotonIso->GetBinLowEdge(max.first)<<","<<SoverRootBvsphotonIso->GetBinLowEdge(max.first+1)<<")   Value: "<<max.second<<endl;

  SoverRootBvsphotonIso->SetTitle("S/#sqrt{B} Vs  photonIso");
  SoverRootBvsphotonIso->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  SoverRootBvsphotonIso->GetXaxis()->SetTitle("photonIso"); 
  SoverRootBvsphotonIso->GetYaxis()->SetTitleOffset(0.75);
  SoverRootBvsphotonIso->GetXaxis()->SetTitleOffset(0.75);
  SoverRootBvsphotonIso->SetMarkerColor(kBlue);
  SoverRootBvsphotonIso->SetMarkerSize(.5);
  SoverRootBvsphotonIso->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(SoverRootBvsphotonIso,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalEffOverRootBgroundEffVsphotonIso.png");
  
  TH1F *LogSoverLogBvsphotonIso =  sig_photonIso->Clone();
  for(int i=1;i<151;i++){
    //cout<<"bin:"<<i<<"   "<<DR03LogSoverLogB[i-1]<<endl;
    if(TMath::IsNaN(DR03LogSoverLogB[i-1]) || !TMath::Finite(DR03LogSoverLogB[i-1]) ) LogSoverLogBvsphotonIso->SetBinContent(i,1);
    else LogSoverLogBvsphotonIso->SetBinContent(i,DR03LogSoverLogB[i-1]);
  }
  LogSoverLogBvsphotonIso->SetTitle("log(S)/log(B) Vs  photonIso");
  LogSoverLogBvsphotonIso->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  LogSoverLogBvsphotonIso->GetXaxis()->SetTitle("photonIso");
  LogSoverLogBvsphotonIso->GetYaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsphotonIso->GetXaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsphotonIso->SetMarkerColor(kBlue);
  LogSoverLogBvsphotonIso->SetMarkerSize(.5);
  LogSoverLogBvsphotonIso->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(LogSoverLogBvsphotonIso,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/LogSignalEffOverLogBgroundEffVsphotonIso.png");

  //Now photonIsoDeposit
  TH1F *sig_photonIsoDeposit = fSig.Get("photonIsoDepositNminus3Signal");
  TH1F *back_photonIsoDeposit = fBack.Get("photonIsoDepositNminus3Bground");//data and qcd
  //if(useQCD)TH1F *back_photonIsoDeposit = fBack.Get("photonIsoDepositNminus1Bground");//qcd


  sig_photonIsoDeposit->SetLineColor(kBlue);
  back_photonIsoDeposit->SetLineColor(kRed);

  TH1F * sig_photonIsoDepositNew = sig_photonIsoDeposit->Clone();
  TH1F* back_photonIsoDepositNew = back_photonIsoDeposit->Clone();

  sig_photonIsoDepositNew->Scale(1/sig_photonIsoDepositNew->Integral());
  back_photonIsoDepositNew->Scale(1/back_photonIsoDepositNew->Integral());
  sig_photonIsoDepositNew->SetTitle("Fraction of Events VS photonIsoDeposit");
  sig_photonIsoDepositNew->GetYaxis()->SetTitle("Fraction of Events");
  sig_photonIsoDepositNew->Draw();
  back_photonIsoDepositNew->Draw("SAME");
  TLegend *photonIsoDeposit = new TLegend(.45,.36,.85,.66,"","brNDC");
  photonIsoDeposit->AddEntry(sig_photonIsoDeposit,"photonIsoDeposit SignalMC_2000_1015_305","lf");
  photonIsoDeposit->AddEntry(back_photonIsoDeposit,"photonIsoDeposit Data MET<30","lf");
  photonIsoDeposit->SetFillColor(kWhite);
  // photonIsoDeposit->SetTextSize(0.03);
  photonIsoDeposit->Draw("SAME");
  c1->Print("Plots/DR03/photonIsoDeposit.png");


  TH1F *sigEffphotonIsoDeposit = sig_photonIsoDeposit->Clone();
  TH1F *backEffphotonIsoDeposit = back_photonIsoDeposit->Clone();
  sigEffphotonIsoDeposit->SetTitle("Integrated Fraction of Events VS photonIsoDeposit");
  sigEffphotonIsoDeposit->GetYaxis()->SetTitle("Integrated Fraction of Events");
  sigEffphotonIsoDeposit->GetXaxis()->SetTitle("photonIsoDeposit");
  sigEffphotonIsoDeposit->GetYaxis()->SetTitleOffset(0.75);
  sigEffphotonIsoDeposit->GetXaxis()->SetTitleOffset(0.75);
  sigEffphotonIsoDeposit->SetLineColor(kBlue);
  backEffphotonIsoDeposit->SetLineColor(kRed);

  for(int i=1;i<151;i++){
    sigEffphotonIsoDeposit->SetBinContent(i,sig_photonIsoDeposit->Integral(0,i)/sig_photonIsoDeposit->Integral());
    backEffphotonIsoDeposit->SetBinContent(i,back_photonIsoDeposit->Integral(0,i)/back_photonIsoDeposit->Integral());
  }


  sigEffphotonIsoDeposit->Draw();
  backEffphotonIsoDeposit->Draw("SAME");
  TLegend *DR03EffphotonDeposit = new TLegend(.35,.3,.8,.5,"","brNDC");
  DR03EffphotonDeposit->AddEntry(sigEffphotonIsoDeposit,"photonIsoDeposit SignalMC_2000_1015_305","lf");
  DR03EffphotonDeposit->AddEntry(backEffphotonIsoDeposit,"photonIsoDeposit Data MET<30","lf");
  DR03EffphotonDeposit->SetFillColor(kWhite);
  //DR03EffphotonDeposit->SetTextSize(0.03);
  DR03EffphotonDeposit->Draw("SAME");
  c1->Print("Plots/DR03/photonIsoDepositEffDR03.png");

  float DR03Sig[150]={0},DR03Bg[150]={0},DR03SoverRootB[150]={0},DR03LogSoverLogB[150]={0};
  for(int i=1;i<151;i++){
    DR03Sig[i-1]=sigEffphotonIsoDeposit->GetBinContent(i);
    DR03Bg[i-1]=backEffphotonIsoDeposit->GetBinContent(i);
    DR03SoverRootB[i-1]=DR03Sig[i-1]/sqrt(DR03Bg[i-1]);
    DR03LogSoverLogB[i-1]=log(DR03Sig[i-1])/log(DR03Bg[i-1]);
  }

  TGraph *bVsSphotonIsoDeposit = new TGraph(150,DR03Bg,DR03Sig);
  bVsSphotonIsoDeposit->SetTitle("Integrated Spectra - SignalMC VS QCD - photonIsoDeposit");
  bVsSphotonIsoDeposit->GetYaxis()->SetTitle("Signal MC");
  bVsSphotonIsoDeposit->GetXaxis()->SetTitle("Background Data MET<30");
  bVsSphotonIsoDeposit->SetMarkerColor(kBlue);
  bVsSphotonIsoDeposit->SetMarkerStyle(20);
  bVsSphotonIsoDeposit->SetMarkerSize(.5);
  bVsSphotonIsoDeposit->Draw("AP");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(bVsSphotonIsoDeposit,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalVsBackgroundEffphotonIsoDeposit.png");

  TH1F *SoverRootBvsphotonIsoDeposit =  sig_photonIsoDeposit->Clone();
  pair<int,float> max;max.first=0;max.second=0.;
  for(int i=1;i<151;i++){
    if(TMath::IsNaN(DR03SoverRootB[i-1]) || !TMath::Finite(DR03SoverRootB[i-1]) ) SoverRootBvsphotonIsoDeposit->SetBinContent(i,1);
    else {
      SoverRootBvsphotonIsoDeposit->SetBinContent(i,DR03SoverRootB[i-1]);
      if(DR03SoverRootB[i-1]>max.second){
	max.first=i;max.second=DR03SoverRootB[i-1];
      }
    }
  }
  cout<<"photonIsoDeposit SoverRootB Max:  (Bin,LowEdge,NextlowEdge) ("<<max.first<<","<<SoverRootBvsphotonIsoDeposit->GetBinLowEdge(max.first)<<","<<SoverRootBvsphotonIsoDeposit->GetBinLowEdge(max.first+1)<<")   Value: "<<max.second<<endl;

  SoverRootBvsphotonIsoDeposit->SetTitle("S/#sqrt{B} Vs  photonIsoDeposit");
  SoverRootBvsphotonIsoDeposit->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  SoverRootBvsphotonIsoDeposit->GetXaxis()->SetTitle("photonIsoDeposit"); 
  SoverRootBvsphotonIsoDeposit->GetYaxis()->SetTitleOffset(0.75);
  SoverRootBvsphotonIsoDeposit->GetXaxis()->SetTitleOffset(0.75);
  SoverRootBvsphotonIsoDeposit->SetMarkerColor(kBlue);
  SoverRootBvsphotonIsoDeposit->SetMarkerSize(.5);
  SoverRootBvsphotonIsoDeposit->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(SoverRootBvsphotonIsoDeposit,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalEffOverRootBgroundEffVsphotonIsoDeposit.png");
  
  TH1F *LogSoverLogBvsphotonIsoDeposit =  sig_photonIsoDeposit->Clone();
  for(int i=1;i<151;i++){
    //cout<<"bin:"<<i<<"   "<<DR03LogSoverLogB[i-1]<<endl;
    if(TMath::IsNaN(DR03LogSoverLogB[i-1]) || !TMath::Finite(DR03LogSoverLogB[i-1]) ) LogSoverLogBvsphotonIsoDeposit->SetBinContent(i,1);
    else LogSoverLogBvsphotonIsoDeposit->SetBinContent(i,DR03LogSoverLogB[i-1]);
  }
  LogSoverLogBvsphotonIsoDeposit->SetTitle("log(S)/log(B) Vs  photonIsoDeposit");
  LogSoverLogBvsphotonIsoDeposit->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  LogSoverLogBvsphotonIsoDeposit->GetXaxis()->SetTitle("photonIsoDeposit");
  LogSoverLogBvsphotonIsoDeposit->GetYaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsphotonIsoDeposit->GetXaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsphotonIsoDeposit->SetMarkerColor(kBlue);
  LogSoverLogBvsphotonIsoDeposit->SetMarkerSize(.5);
  LogSoverLogBvsphotonIsoDeposit->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(LogSoverLogBvsphotonIsoDeposit,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/LogSignalEffOverLogBgroundEffVsphotonIsoDeposit.png");

 //Now PfCombinedIso
  TH1F *sig_PfCombinedIso = fSig.Get("PfCombinedIsoNminus3Signal");
  TH1F *back_PfCombinedIso = fBack.Get("PfCombinedIsoNminus3Bground");//data and qcd
  //if(useQCD)TH1F *back_PfCombinedIso = fBack.Get("PfCombinedIsoNminus1Bground");//qcd


  sig_PfCombinedIso->SetLineColor(kBlue);
  back_PfCombinedIso->SetLineColor(kRed);

  TH1F * sig_PfCombinedIsoNew = sig_PfCombinedIso->Clone();
  TH1F* back_PfCombinedIsoNew = back_PfCombinedIso->Clone();

  sig_PfCombinedIsoNew->Scale(1/sig_PfCombinedIsoNew->Integral());
  back_PfCombinedIsoNew->Scale(1/back_PfCombinedIsoNew->Integral());
  sig_PfCombinedIsoNew->SetTitle("Fraction of Events VS PfCombinedIso");
  sig_PfCombinedIsoNew->GetYaxis()->SetTitle("Fraction of Events");
  sig_PfCombinedIsoNew->Draw();
  back_PfCombinedIsoNew->Draw("SAME");
  TLegend *DR03Iso = new TLegend(.45,.36,.85,.66,"","brNDC");
  DR03Iso->AddEntry(sig_PfCombinedIso,"PfCombinedIso SignalMC_2000_1015_305","lf");
  DR03Iso->AddEntry(back_PfCombinedIso,"PfCombinedIso Data MET<30","lf");
  DR03Iso->SetFillColor(kWhite);
  // DR03Iso->SetTextSize(0.03);
  DR03Iso->Draw("SAME");
  c1->Print("Plots/DR03/PfCombinedIso.png");


  TH1F *sigEffPfCombinedIso = sig_PfCombinedIso->Clone();
  TH1F *backEffPfCombinedIso = back_PfCombinedIso->Clone();
  sigEffPfCombinedIso->SetTitle("Integrated Fraction of Events VS PfCombinedIso");
  sigEffPfCombinedIso->GetYaxis()->SetTitle("Integrated Fraction of Events");
  sigEffPfCombinedIso->GetXaxis()->SetTitle("PfCombinedIso");
  sigEffPfCombinedIso->GetYaxis()->SetTitleOffset(0.75);
  sigEffPfCombinedIso->GetXaxis()->SetTitleOffset(0.75);
  sigEffPfCombinedIso->SetLineColor(kBlue);
  backEffPfCombinedIso->SetLineColor(kRed);

  for(int i=1;i<151;i++){
    sigEffPfCombinedIso->SetBinContent(i,sig_PfCombinedIso->Integral(0,i)/sig_PfCombinedIso->Integral());
    backEffPfCombinedIso->SetBinContent(i,back_PfCombinedIso->Integral(0,i)/back_PfCombinedIso->Integral());
  }


  sigEffPfCombinedIso->Draw();
  backEffPfCombinedIso->Draw("SAME");
  TLegend *DR03EffPfCombined = new TLegend(.35,.3,.8,.5,"","brNDC");
  DR03EffPfCombined->AddEntry(sigEffPfCombinedIso,"PfCombinedIso SignalMC_2000_1015_305","lf");
  DR03EffPfCombined->AddEntry(backEffPfCombinedIso,"PfCombinedIso Data MET<30","lf");
  DR03EffPfCombined->SetFillColor(kWhite);
  //DR03EffPfCombined->SetTextSize(0.03);
  DR03EffPfCombined->Draw("SAME");
  c1->Print("Plots/DR03/PfCombinedIsoEffDR03.png");

  float DR03Sig[150]={0},DR03Bg[150]={0},DR03SoverRootB[150]={0},DR03LogSoverLogB[150]={0};
  for(int i=1;i<151;i++){
    DR03Sig[i-1]=sigEffPfCombinedIso->GetBinContent(i);
    DR03Bg[i-1]=backEffPfCombinedIso->GetBinContent(i);
    DR03SoverRootB[i-1]=DR03Sig[i-1]/sqrt(DR03Bg[i-1]);
    DR03LogSoverLogB[i-1]=log(DR03Sig[i-1])/log(DR03Bg[i-1]);
  }

  TGraph *bVsSPfCombinedIso = new TGraph(150,DR03Bg,DR03Sig);
  bVsSPfCombinedIso->SetTitle("Integrated Spectra - SignalMC VS QCD - PfCombinedIso");
  bVsSPfCombinedIso->GetYaxis()->SetTitle("Signal MC");
  bVsSPfCombinedIso->GetXaxis()->SetTitle("Background Data MET<30");
  bVsSPfCombinedIso->SetMarkerColor(kBlue);
  bVsSPfCombinedIso->SetMarkerStyle(20);
  bVsSPfCombinedIso->SetMarkerSize(.5);
  bVsSPfCombinedIso->Draw("AP");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(bVsSPfCombinedIso,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalVsBackgroundEffPfCombined.png");

  TH1F *SoverRootBvsPfCombinedIso =  sig_PfCombinedIso->Clone();
  pair<int,float> max;max.first=0;max.second=0.;
  for(int i=1;i<151;i++){
    if(TMath::IsNaN(DR03SoverRootB[i-1]) || !TMath::Finite(DR03SoverRootB[i-1]) ) SoverRootBvsPfCombinedIso->SetBinContent(i,1);
    else {
      SoverRootBvsPfCombinedIso->SetBinContent(i,DR03SoverRootB[i-1]);
      if(DR03SoverRootB[i-1]>max.second){
	max.first=i;max.second=DR03SoverRootB[i-1];
      }
    }
  }
  cout<<"PfCombined SoverRootB Max:  (Bin,LowEdge,NextlowEdge) ("<<max.first<<","<<SoverRootBvsPfCombinedIso->GetBinLowEdge(max.first)<<","<<SoverRootBvsPfCombinedIso->GetBinLowEdge(max.first+1)<<")   Value: "<<max.second<<endl;

  SoverRootBvsPfCombinedIso->SetTitle("S/#sqrt{B} Vs  PfCombinedIso");
  SoverRootBvsPfCombinedIso->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  SoverRootBvsPfCombinedIso->GetXaxis()->SetTitle("PfCombinedIso"); 
  SoverRootBvsPfCombinedIso->GetYaxis()->SetTitleOffset(0.75);
  SoverRootBvsPfCombinedIso->GetXaxis()->SetTitleOffset(0.75);
  SoverRootBvsPfCombinedIso->SetMarkerColor(kBlue);
  SoverRootBvsPfCombinedIso->SetMarkerSize(.5);
  SoverRootBvsPfCombinedIso->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(SoverRootBvsPfCombinedIso,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalEffOverRootBgroundEffVsPfCombinedIso.png");
  
  TH1F *LogSoverLogBvsPfCombinedIso =  sig_PfCombinedIso->Clone();
  for(int i=1;i<151;i++){
    //cout<<"bin:"<<i<<"   "<<DR03LogSoverLogB[i-1]<<endl;
    if(TMath::IsNaN(DR03LogSoverLogB[i-1]) || !TMath::Finite(DR03LogSoverLogB[i-1]) ) LogSoverLogBvsPfCombinedIso->SetBinContent(i,1);
    else LogSoverLogBvsPfCombinedIso->SetBinContent(i,DR03LogSoverLogB[i-1]);
  }
  LogSoverLogBvsPfCombinedIso->SetTitle("log(S)/log(B) Vs  PfCombinedIso");
  LogSoverLogBvsPfCombinedIso->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  LogSoverLogBvsPfCombinedIso->GetXaxis()->SetTitle("PfCombinedIso");
  LogSoverLogBvsPfCombinedIso->GetYaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsPfCombinedIso->GetXaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsPfCombinedIso->SetMarkerColor(kBlue);
  LogSoverLogBvsPfCombinedIso->SetMarkerSize(.5);
  LogSoverLogBvsPfCombinedIso->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(LogSoverLogBvsPfCombinedIso,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/LogSignalEffOverLogBgroundEffVsPfCombinedIso.png");

  //Now PfCombinedIsoDeposit
  TH1F *sig_PfCombinedIsoDeposit = fSig.Get("PfCombinedIsoDepositNminus3Signal");
  TH1F *back_PfCombinedIsoDeposit = fBack.Get("PfCombinedIsoDepositNminus3Bground");//data and qcd
  //if(useQCD)TH1F *back_PfCombinedIsoDeposit = fBack.Get("PfCombinedIsoDepositNminus1Bground");//qcd


  sig_PfCombinedIsoDeposit->SetLineColor(kBlue);
  back_PfCombinedIsoDeposit->SetLineColor(kRed);

  TH1F * sig_PfCombinedIsoDepositNew = sig_PfCombinedIsoDeposit->Clone();
  TH1F* back_PfCombinedIsoDepositNew = back_PfCombinedIsoDeposit->Clone();

  sig_PfCombinedIsoDepositNew->Scale(1/sig_PfCombinedIsoDepositNew->Integral());
  back_PfCombinedIsoDepositNew->Scale(1/back_PfCombinedIsoDepositNew->Integral());
  sig_PfCombinedIsoDepositNew->SetTitle("Fraction of Events VS PfCombinedIsoDeposit");
  sig_PfCombinedIsoDepositNew->GetYaxis()->SetTitle("Fraction of Events");
  sig_PfCombinedIsoDepositNew->Draw();
  back_PfCombinedIsoDepositNew->Draw("SAME");
  TLegend *PfCombinedIsoDeposit = new TLegend(.45,.36,.85,.66,"","brNDC");
  PfCombinedIsoDeposit->AddEntry(sig_PfCombinedIsoDeposit,"PfCombinedIsoDeposit SignalMC_2000_1015_305","lf");
  PfCombinedIsoDeposit->AddEntry(back_PfCombinedIsoDeposit,"PfCombinedIsoDeposit Data MET<30","lf");
  PfCombinedIsoDeposit->SetFillColor(kWhite);
  // PfCombinedIsoDeposit->SetTextSize(0.03);
  PfCombinedIsoDeposit->Draw("SAME");
  c1->Print("Plots/DR03/PfCombinedIsoDeposit.png");


  TH1F *sigEffPfCombinedIsoDeposit = sig_PfCombinedIsoDeposit->Clone();
  TH1F *backEffPfCombinedIsoDeposit = back_PfCombinedIsoDeposit->Clone();
  sigEffPfCombinedIsoDeposit->SetTitle("Integrated Fraction of Events VS PfCombinedIsoDeposit");
  sigEffPfCombinedIsoDeposit->GetYaxis()->SetTitle("Integrated Fraction of Events");
  sigEffPfCombinedIsoDeposit->GetXaxis()->SetTitle("PfCombinedIsoDeposit");
  sigEffPfCombinedIsoDeposit->GetYaxis()->SetTitleOffset(0.75);
  sigEffPfCombinedIsoDeposit->GetXaxis()->SetTitleOffset(0.75);
  sigEffPfCombinedIsoDeposit->SetLineColor(kBlue);
  backEffPfCombinedIsoDeposit->SetLineColor(kRed);

  for(int i=1;i<151;i++){
    sigEffPfCombinedIsoDeposit->SetBinContent(i,sig_PfCombinedIsoDeposit->Integral(0,i)/sig_PfCombinedIsoDeposit->Integral());
    backEffPfCombinedIsoDeposit->SetBinContent(i,back_PfCombinedIsoDeposit->Integral(0,i)/back_PfCombinedIsoDeposit->Integral());
  }


  sigEffPfCombinedIsoDeposit->Draw();
  backEffPfCombinedIsoDeposit->Draw("SAME");
  TLegend *DR03EffPfCombinedDeposit = new TLegend(.35,.3,.8,.5,"","brNDC");
  DR03EffPfCombinedDeposit->AddEntry(sigEffPfCombinedIsoDeposit,"PfCombinedIsoDeposit SignalMC_2000_1015_305","lf");
  DR03EffPfCombinedDeposit->AddEntry(backEffPfCombinedIsoDeposit,"PfCombinedIsoDeposit Data MET<30","lf");
  DR03EffPfCombinedDeposit->SetFillColor(kWhite);
  //DR03EffPfCombinedDeposit->SetTextSize(0.03);
  DR03EffPfCombinedDeposit->Draw("SAME");
  c1->Print("Plots/DR03/PfCombinedIsoDepositEffDR03.png");

  float DR03Sig[150]={0},DR03Bg[150]={0},DR03SoverRootB[150]={0},DR03LogSoverLogB[150]={0};
  for(int i=1;i<151;i++){
    DR03Sig[i-1]=sigEffPfCombinedIsoDeposit->GetBinContent(i);
    DR03Bg[i-1]=backEffPfCombinedIsoDeposit->GetBinContent(i);
    DR03SoverRootB[i-1]=DR03Sig[i-1]/sqrt(DR03Bg[i-1]);
    DR03LogSoverLogB[i-1]=log(DR03Sig[i-1])/log(DR03Bg[i-1]);
  }

  TGraph *bVsSPfCombinedIsoDeposit = new TGraph(150,DR03Bg,DR03Sig);
  bVsSPfCombinedIsoDeposit->SetTitle("Integrated Spectra - SignalMC VS QCD - PfCombinedIsoDeposit");
  bVsSPfCombinedIsoDeposit->GetYaxis()->SetTitle("Signal MC");
  bVsSPfCombinedIsoDeposit->GetXaxis()->SetTitle("Background Data MET<30");
  bVsSPfCombinedIsoDeposit->SetMarkerColor(kBlue);
  bVsSPfCombinedIsoDeposit->SetMarkerStyle(20);
  bVsSPfCombinedIsoDeposit->SetMarkerSize(.5);
  bVsSPfCombinedIsoDeposit->Draw("AP");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(bVsSPfCombinedIsoDeposit,"PfCombinedIsoDeposit","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalVsBackgroundEffPfCombinedIsoDeposit.png");

  bVsSPfCombinedIsoDeposit->Draw("AP"); 
  bVsSPfCombinedIso->SetMarkerColor(kRed);
  bVsSPfCombinedIso->Draw("Psame");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(bVsSPfCombinedIso,"PfCombinedIso","p");
  Eff->AddEntry(bVsSPfCombinedIsoDeposit,"PfCombinedIsoDeposit","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalVsBackgroundEffPfCombinedIsoAndIsoDeposit.png");
  bVsSCombDR03->SetMarkerColor(kGreen);
  bVsSCombDR03->Draw("sameP");
  Eff->AddEntry(bVsSCombDR03,"DR03","p");
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalVsBackgroundEffPfCombinedIsoAndIsoDepositAndDR03.png");


  TH1F *SoverRootBvsPfCombinedIsoDeposit =  sig_PfCombinedIsoDeposit->Clone();
  pair<int,float> max;max.first=0;max.second=0.;
  for(int i=1;i<151;i++){
    if(TMath::IsNaN(DR03SoverRootB[i-1]) || !TMath::Finite(DR03SoverRootB[i-1]) ) SoverRootBvsPfCombinedIsoDeposit->SetBinContent(i,1);
    else {
      SoverRootBvsPfCombinedIsoDeposit->SetBinContent(i,DR03SoverRootB[i-1]);
      if(DR03SoverRootB[i-1]>max.second){
	max.first=i;max.second=DR03SoverRootB[i-1];
      }
    }
  }
  cout<<"PfCombinedIsoDeposit SoverRootB Max:  (Bin,LowEdge,NextlowEdge) ("<<max.first<<","<<SoverRootBvsPfCombinedIsoDeposit->GetBinLowEdge(max.first)<<","<<SoverRootBvsPfCombinedIsoDeposit->GetBinLowEdge(max.first+1)<<")   Value: "<<max.second<<endl;

  SoverRootBvsPfCombinedIsoDeposit->SetTitle("S/#sqrt{B} Vs  PfCombinedIsoDeposit");
  SoverRootBvsPfCombinedIsoDeposit->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  SoverRootBvsPfCombinedIsoDeposit->GetXaxis()->SetTitle("PfCombinedIsoDeposit"); 
  SoverRootBvsPfCombinedIsoDeposit->GetYaxis()->SetTitleOffset(0.75);
  SoverRootBvsPfCombinedIsoDeposit->GetXaxis()->SetTitleOffset(0.75);
  SoverRootBvsPfCombinedIsoDeposit->SetMarkerColor(kBlue);
  SoverRootBvsPfCombinedIsoDeposit->SetMarkerSize(.5);
  SoverRootBvsPfCombinedIsoDeposit->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(SoverRootBvsPfCombinedIsoDeposit,"DR03","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalEffOverRootBgroundEffVsPfCombinedIsoDeposit.png");
  
  TH1F *LogSoverLogBvsPfCombinedIsoDeposit =  sig_PfCombinedIsoDeposit->Clone();
  for(int i=1;i<151;i++){
    //cout<<"bin:"<<i<<"   "<<DR03LogSoverLogB[i-1]<<endl;
    if(TMath::IsNaN(DR03LogSoverLogB[i-1]) || !TMath::Finite(DR03LogSoverLogB[i-1]) ) LogSoverLogBvsPfCombinedIsoDeposit->SetBinContent(i,1);
    else LogSoverLogBvsPfCombinedIsoDeposit->SetBinContent(i,DR03LogSoverLogB[i-1]);
  }
  LogSoverLogBvsPfCombinedIsoDeposit->SetTitle("log(S)/log(B) Vs  PfCombinedIsoDeposit");
  LogSoverLogBvsPfCombinedIsoDeposit->GetYaxis()->SetTitle("#epsilon_{S}/#sqrt{#epsilon_{B}}");
  LogSoverLogBvsPfCombinedIsoDeposit->GetXaxis()->SetTitle("PfCombinedIsoDeposit");
  LogSoverLogBvsPfCombinedIsoDeposit->GetYaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsPfCombinedIsoDeposit->GetXaxis()->SetTitleOffset(0.75);
  LogSoverLogBvsPfCombinedIsoDeposit->SetMarkerColor(kBlue);
  LogSoverLogBvsPfCombinedIsoDeposit->SetMarkerSize(.5);
  LogSoverLogBvsPfCombinedIsoDeposit->Draw("P");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(LogSoverLogBvsPfCombinedIsoDeposit,"PfCombinedIsoDeposit","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/LogSignalEffOverLogBgroundEffVsPfCombinedIsoDeposit.png");

  SoverRootBvsPfCombinedIsoDeposit->GetYaxis()->SetRangeUser(0,1.6);
  SoverRootBvsPfCombinedIsoDeposit->Draw("P");
  SoverRootBvsPfCombinedIso->SetMarkerColor(kRed);
  SoverRootBvsPfCombinedIso->Draw("Psame");
  TLegend *Eff = new TLegend(.55,.3,.8,.5,"","brNDC");
  Eff->AddEntry(SoverRootBvsPfCombinedIso,"PfCombinedIso","p");
  Eff->AddEntry(SoverRootBvsPfCombinedIsoDeposit,"PfCombinedIsoDeposit","p");
  Eff->SetFillColor(kWhite);
  Eff->Draw("SAME");
  c1->Print("Plots/DR03/SignalEffOverRootBgroundEffVsPfCombinedIsoAndIsoDeposit.png");


  return;
}

#define HiggsTools_cxx
#include "HiggsTools.h"
#include <iostream>
#include <TLegend.h>
#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include <vector>
#include "RooGlobalFunc.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooVoigtian.h"
#include "RooPlot.h"
#include "RooMath.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooNumConvPdf.h"
#include "RooBreitWigner.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooBernstein.h"
#include "TProfile.h"
#include "TH2.h"
#include "TH2F.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TString.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
#include "TF1.h"
#include "PhysicsTools/TagAndProbe/interface/RooCMSShape.h"
#include "THStack.h"
#include "TLatex.h"

using namespace std;
using namespace RooFit;

float L_int = 19499.;//full dataset

void HiggsTools::Loop(){

  TFile fin("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_15GevCHfakeCut_EleVeto_NewDiEMPtBins_19499pb_NewJetMatching_2JetReq_May16.root","READ");
  TFile f_ggHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_May7.root","READ");
  TFile f_WZHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_May7.root","READ");
  TFile f_TTHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_TTH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_May7.root","READ");
  TFile f_VBFHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_May7.root","READ");
  //TFile f_VBFHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");

  TCanvas *c1 = new TCanvas("c1","",900,600);
  c1->cd();

  TH1F* ggInvarMassMET = (TH1F*)fin.Get("ggInvarMassMVAcorrVertexCorr");ggInvarMassMET->Sumw2();;ggInvarMassMET->SetMarkerSize(.5);
  TH1F* ggInvarMassMET_JetReq=(TH1F*)fin.Get("ggInvarMass_JetReq");ggInvarMassMET_JetReq->Sumw2();
  TH1F* ggInvarMassMET30 = (TH1F*)fin.Get("ggInvarMassMET30MVAcorr");
  TH1F* ggInvarMassMET30_JetReq = (TH1F*)fin.Get("ggInvarMassMET30_JetReq");
  //ggInvarMassMET->Rebin(3);ggInvarMassMET_JetReq->Rebin(3);
  ggInvarMassMET30->Rebin(3);ggInvarMassMET30_JetReq->Rebin(3);

  TH1F* ggHggNoMET = (TH1F*)f_ggHgg.Get("ggInvarMassMVAcorrVertexCorr");
  TH1F* WZHggNoMET = (TH1F*)f_WZHgg.Get("ggInvarMassMVAcorrVertexCorr");
  TH1F* TTHggNoMET = (TH1F*)f_TTHgg.Get("ggInvarMassMVAcorrVertexCorr");
  TH1F* VBFHggNoMET = (TH1F*)f_VBFHgg.Get("ggInvarMassMVAcorrVertexCorr");
  TH1F* ggHggMET30 = (TH1F*)f_ggHgg.Get("ggInvarMassMET30MVAcorr");
  TH1F* WZHggMET30 = (TH1F*)f_WZHgg.Get("ggInvarMassMET30MVAcorr");
  TH1F* TTHggMET30 = (TH1F*)f_TTHgg.Get("ggInvarMassMET30MVAcorr");
  TH1F* VBFHggMET30 = (TH1F*)f_VBFHgg.Get("ggInvarMassMET30MVAcorr");

  //fix these to have mvacorrvertexcorr in jetreq - fixed, now all invarmass are mvacorrvertexcorr
  TH1F* ggHggNoMET_JetReq = (TH1F*)f_ggHgg.Get("ggInvarMass_JetReq");
  TH1F* WZHggNoMET_JetReq = (TH1F*)f_WZHgg.Get("ggInvarMass_JetReq");
  TH1F* TTHggNoMET_JetReq = (TH1F*)f_TTHgg.Get("ggInvarMass_JetReq");
  TH1F* VBFHggNoMET_JetReq = (TH1F*)f_VBFHgg.Get("ggInvarMass_JetReq");

  ggHggNoMET->Sumw2();WZHggNoMET->Sumw2();TTHggNoMET->Sumw2();VBFHggNoMET->Sumw2();
  ggHggMET30->Sumw2();WZHggMET30->Sumw2();TTHggMET30->Sumw2();VBFHggMET30->Sumw2();
  ggHggNoMET_JetReq->Sumw2();WZHggNoMET_JetReq->Sumw2();TTHggNoMET_JetReq->Sumw2();VBFHggNoMET_JetReq->Sumw2();


  ggHggNoMET->Scale((L_int*2.29e-03*19.52)/100000);  
  VBFHggNoMET->Scale((L_int*2.29e-03*1.559)/100000); 
  TTHggNoMET->Scale((L_int*2.29e-03*.1302)/100000);//had some failed jobs, using number ntuplized instead of 100224(number produced)
  WZHggNoMET->Scale((L_int*2.29e-03*(.6966/**(.3257+.014)*/+.3943/**.2*/))/100000);
  ggHggNoMET->SetLineColor(kGreen);ggHggNoMET->SetMarkerColor(kGreen);
  WZHggNoMET->SetLineColor(kCyan);WZHggNoMET->SetMarkerColor(kCyan);
  TTHggNoMET->SetLineColor(kViolet);TTHggNoMET->SetMarkerColor(kViolet);
  VBFHggNoMET->SetLineColor(kRed+3);VBFHggNoMET->SetMarkerColor(kRed+3);
  ggHggNoMET->SetFillStyle(0);WZHggNoMET->SetFillStyle(0);TTHggNoMET->SetFillStyle(0);VBFHggNoMET->SetFillStyle(0);
  ggHggNoMET->SetLineWidth(2);WZHggNoMET->SetLineWidth(2);VBFHggNoMET->SetLineWidth(2);TTHggNoMET->SetLineWidth(2);
  ggHggMET30->Scale((L_int*2.29e-03*19.52)/100000);  
  VBFHggMET30->Scale((L_int*2.29e-03*1.559)/100000); 
  TTHggMET30->Scale((L_int*2.29e-03*.1302)/100000);//had some failed jobs, using number ntuplized instead of 100224(number produced)
  WZHggMET30->Scale((L_int*2.29e-03*(.6966/**(.3257+.014)*/+.3943/**.2*/))/100000);
  ggHggMET30->SetLineColor(kGreen);ggHggMET30->SetMarkerColor(kGreen);
  WZHggMET30->SetLineColor(kCyan);WZHggMET30->SetMarkerColor(kCyan);
  TTHggMET30->SetLineColor(kViolet);TTHggMET30->SetMarkerColor(kViolet);
  VBFHggMET30->SetLineColor(kRed+3);VBFHggMET30->SetMarkerColor(kRed+3);
  ggHggMET30->SetFillStyle(0);WZHggMET30->SetFillStyle(0);TTHggMET30->SetFillStyle(0);VBFHggMET30->SetFillStyle(0);
  ggHggMET30->SetLineWidth(2);WZHggMET30->SetLineWidth(2);VBFHggMET30->SetLineWidth(2);TTHggMET30->SetLineWidth(2);
  ggHggNoMET_JetReq->Scale((L_int*2.29e-03*19.52)/100000);  
  VBFHggNoMET_JetReq->Scale((L_int*2.29e-03*1.559)/100000); 
  TTHggNoMET_JetReq->Scale((L_int*2.29e-03*.1302)/100000);//had some failed jobs, using number ntuplized instead of 100224(number produced)
  WZHggNoMET_JetReq->Scale((L_int*2.29e-03*(.6966/**(.3257+.014)*/+.3943/**.2*/))/100000);
  ggHggNoMET_JetReq->SetLineColor(kGreen);ggHggNoMET_JetReq->SetMarkerColor(kGreen);
  WZHggNoMET_JetReq->SetLineColor(kCyan);WZHggNoMET_JetReq->SetMarkerColor(kCyan);
  TTHggNoMET_JetReq->SetLineColor(kViolet);TTHggNoMET_JetReq->SetMarkerColor(kViolet);
  VBFHggNoMET_JetReq->SetLineColor(kRed+3);VBFHggNoMET_JetReq->SetMarkerColor(kRed+3);
  ggHggNoMET_JetReq->SetFillStyle(0);WZHggNoMET_JetReq->SetFillStyle(0);TTHggNoMET_JetReq->SetFillStyle(0);VBFHggNoMET_JetReq->SetFillStyle(0);
  ggHggNoMET_JetReq->SetLineWidth(2);WZHggNoMET_JetReq->SetLineWidth(2);VBFHggNoMET_JetReq->SetLineWidth(2);TTHggNoMET_JetReq->SetLineWidth(2);

  //ggHggNoMET->Rebin(3);WZHggNoMET->Rebin(3);TTHggNoMET->Rebin(3);VBFHggNoMET->Rebin(3);
  ggHggMET30->Rebin(3);WZHggMET30->Rebin(3);TTHggMET30->Rebin(3);VBFHggMET30->Rebin(3);
  //ggHggNoMET_JetReq->Rebin(3);WZHggNoMET_JetReq->Rebin(3);TTHggNoMET_JetReq->Rebin(3);VBFHggNoMET_JetReq->Rebin(3);


  TH1F* SMHiggs = (TH1F*)ggHggNoMET->Clone();
  SMHiggs->Add(WZHggNoMET);
  SMHiggs->Add(VBFHggNoMET);
  SMHiggs->Add(TTHggNoMET);
  RooRealVar xSM("xSM","m_{#gamma#gamma}",110,140,"GeV");xSM.setRange("sigSMrange",115,135);
  RooDataHist dataSM("dataSM","data SM higgs",xSM,SMHiggs);
  //Gaussian for peak
  RooRealVar meanSM("meanSM","meanSM",125/*.3,122,128*/);
  RooRealVar sigmaSM("sigmaSM","sigmaSM",0.06,0,1.2);
  RooRealVar sigmaSM2("sigmaSM2","sigmaSM2",0.,6);
  //Voigtian (Breit-Vigner x Gaussian) for peak
  RooRealVar widthSM("widthSM","widthSM",1.9,0.1,3.7);
  RooRealVar SMYield("sig1 yield","sig1 yield",20000,0,400000);
  RooRealVar SMYield2("sig2 yield","sig2 yield",20000,0,400000);
  RooGaussian sigSM2("sigSM2","gaussSM2",xSM,meanSM,sigmaSM2);
  //RooGaussian sigSM ("sigSM", "gaussSM", xSM,meanSM,sigmaSM);
  RooVoigtian sigSM ("sigSM", "gaussSM", xSM,meanSM,widthSM,sigmaSM);
  //RooAddPdf SMPdf("SMPdf","SMPdf",RooArgList(sigSM,sigSM2),RooArgList(SMYield,SMYield2));
  RooAddPdf SMPdf("SMPdf","SMPdf",RooArgList(sigSM),RooArgList(SMYield));
  RooFitResult *rSM = SMPdf.fitTo(dataSM,Extended(kTRUE),Save(),Range("sigSMrange"));
  RooPlot *xframeSM = xSM.frame(Title("Breit-Wigner#otimesGaussian signal, gg SM higgs"));
  SMPdf.paramOn(xframeSM, Format("NE",AutoPrecision(2)),Layout(0.6,0.995,0.9) );
  dataSM.plotOn(xframeSM,LineColor(kBlack));
  SMPdf.plotOn(xframeSM,LineColor(kBlue));
  SMPdf.plotOn(xframeSM,LineColor(kBlue),VisualizeError(*rSM,2,kTRUE),FillColor(kGreen));
  SMPdf.plotOn(xframeSM,LineColor(kBlue),VisualizeError(*rSM,1,kTRUE),FillColor(kOrange));
  //SMPdf.plotOn(xframeSM,Components(sigSM),LineColor(kRed),LineStyle(kDashed));
  //SMPdf.plotOn(xframeSM,Components(sigSM2),LineColor(kGreen),LineStyle(kDashed));
  dataSM.plotOn(xframeSM,LineColor(kBlack));
  xframeSM->Draw();
  THStack* StackMet = new THStack("StackMet","");
  StackMet->Add(TTHggNoMET);//->Draw("SAME");
  StackMet->Add(WZHggNoMET);//->Draw("SAME"); 
  StackMet->Add(VBFHggNoMET);//->Draw("SAME");
  StackMet->Add(ggHggNoMET);//->Draw("SAME");
  StackMet->Draw("histoSAMES");
  xframeSM->Draw("SAMES");
  TLegend *legHiggs = new TLegend(.33,.65,.56,.8);
  legHiggs->AddEntry(ggHggNoMET,"GluGlu->H->#gamma#gamma","l");
  legHiggs->AddEntry(WZHggNoMET,"W/ZH->#gamma#gamma","l");
  legHiggs->AddEntry(VBFHggNoMET,"VBFH->#gamma#gamma","l");
  legHiggs->AddEntry(TTHggNoMET,"TTH->#gamma#gamma","l");
  legHiggs->SetFillColor(kWhite);
  legHiggs->SetFillStyle(0);
  legHiggs->SetBorderSize(0);
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassSMhiggs.png");
 

 //data, no met cut

  RooRealVar xMet("xMet","m_{#gamma#gamma}",95,200,"GeV");
  xMet.setRange("full",95,200);
  RooDataHist dataMet("dataMet","dataset",xMet,ggInvarMassMET);
  //Gaussian for peak
  RooRealVar gmeanMet("gmeanMet","gmean",125.3,122,128);
  RooRealVar gsigmaMet("gsigmaMet","gsigma",1.6/*,0.1,5*/);
  //RooGaussian sigMet("sigMet","gauss",xMet,gmeanMet,gsigmaMet);
  //Crystal Ball for signal
  RooRealVar meanMet("cb_mean", "mean" , 125.3, 122, 128.) ;
  RooRealVar sigmaMet("cb_sigma", "sigma",1.6/*2.6, 0., 3.*/);
  RooRealVar nMet("n","n", 10.,-5.,25.);
  RooRealVar alphaMet("alpha","alpha",25.,0.,50.);
  //Voigtian (Breit-Vigner x Gaussian) for peak
  RooRealVar VmeanMet("Vmean", "Vmean" , 125.3,122.,128.) ;
  RooRealVar VsigmaMet("VsigmaMet","VsigmaMet",0.06);
  RooRealVar VwidthMet("VwidthMet","VwidthMet",1.98);
  RooVoigtian sigMet ("sigMet", "gaussMet", xMet,VmeanMet,VwidthMet,VsigmaMet);
  //RooCBShape sigMet("sigMet", "crystal ball", xMet, meanMet, sigmaMet, alphaMet, nMet);
  RooRealVar sigMetYield("signal yield","signal yield",100,0,600);
  RooRealVar Bern1("Bern1","Berstein 1",14.5,8.,22.);
  RooRealVar Bern2("Bern2","Berstein 2",5.,1.,9.);
  RooRealVar Bern3("Bern3","Berstein 3",4.,0.,8.);
  RooRealVar Bern4("Bern4","Berstein 4",2.,0.,4.);
  RooRealVar Bern5("Bern5","Berstein 5",5.,0.,10.);
  RooRealVar Pol1("Pol1","Pol1",0,-.1,.1);
  RooRealVar Pol2("Pol2","Pol2",0,-.001,.001);
  RooRealVar Pol3("Pol3","Pol3",0,-.001,.001);
  RooRealVar Pol4("Pol4","Pol4",0,-.001,.001);
  RooRealVar Pol5("Pol5","Pol5",0,-.0001,.0001);
  //RooBernstein Bern("Bern","4th order Bernstein Polynomial",xMet,RooArgList(Bern1,Bern2,Bern3,Bern4,Bern5));
  RooPolynomial Bern("Bern","4th Order Polynomial",xMet,RooArgList(Pol1,Pol2,Pol3,Pol4));
  RooRealVar BernYield("bkgd yield","bkgd yield",20000,0,400000);
  RooAddPdf MetPdf("MetPdf","MetPdf",RooArgList(Bern,sigMet),RooArgList(BernYield,sigMetYield));
  RooFitResult *rNoMetCut = MetPdf.fitTo(dataMet,Extended(kTRUE),Save());
  RooPlot *xframeMet = xMet.frame(Title("Voigtian Signal, 4th Order Polynomial Background, gg with no MET cut"));
  MetPdf.paramOn(xframeMet, Format("NE",AutoPrecision(2)),Layout(0.6,0.995,0.9) );
  dataMet.plotOn(xframeMet,LineColor(kBlack),MarkerSize(0.3),MarkerStyle(20));
  //MetPdf.plotOn(xframeMet,Components(Bern),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut,3,kTRUE),FillColor(kViolet));
  MetPdf.plotOn(xframeMet,Components(Bern),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut,2,kTRUE),FillColor(kGreen));
  MetPdf.plotOn(xframeMet,Components(Bern),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut,1,kTRUE),FillColor(kOrange));
  //MetPdf.plotOn(xframeMet,LineColor(kBlue),LineStyle(kDashed),VisualizeError(*rNoMetCut,1,kTRUE),FillColor(kOrange));
  dataMet.plotOn(xframeMet,LineColor(kBlack),MarkerSize(0.3),MarkerStyle(20));
  MetPdf.plotOn(xframeMet,Components(Bern),LineColor(kRed),LineStyle(kDashed),Range("full"));
  MetPdf.plotOn(xframeMet,LineColor(kBlue));
  MetPdf.plotOn(xframeMet,Components(sigMet),LineColor(kBlue));
  dataMet.plotOn(xframeMet,LineColor(kBlack),MarkerSize(0.3),MarkerStyle(20));
  xframeMet->SetAxisRange(105.5,149,"X");
  xframeMet->SetAxisRange(0,1950,"Y");
  xframeMet->Draw();
  //line125->DrawLine(125.3,0,125.3,xframeMet->GetMaximum());
  StackMet->Draw("histoSAME");
  xframeMet->Draw("SAME");
  legHiggs->Draw();
  Double_t sigSMyield=0.;
  Double_t sigSMErr=0.;
  ggHggNoMET->Add(WZHggNoMET);ggHggNoMET->Add(VBFHggNoMET);ggHggNoMET->Add(TTHggNoMET);
  sigSMyield = ggHggNoMET->IntegralAndError(0,999,sigSMErr,"");
  //sigSMyield/=10;sigSMErr/=10;
  Double_t sigYield=sigMetYield.getVal();
  Double_t sigYieldErr=sigMetYield.getError();
  //cout<<"sigSMyield: "<<sigSMyield<<"  sigSMErr: "<<sigSMErr<<endl;
  //cout<<"sigYield: "<<sigYield<<"  sigYieldErr: "<<sigYieldErr<<endl;
  //sigStrengthVsMet->SetBinContent(1,sigYield/sigSMyield);
  //sigStrengthVsMet->SetBinError(1,sqrt(sigYieldErr*sigYieldErr/(sigSMyield*sigSMyield) + (sigYield*sigYield*sigSMErr*sigSMErr)/(sigSMyield*sigSMyield)));
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithNoMETcut.png");
  
  RooRealVar xMet_bg("xMet_bg","m_{#gamma#gamma}",95,200,"GeV");
  xMet_bg.setRange("sb_lo",95,120);xMet_bg.setRange("sb_hi",130,200);xMet.setRange("full",95,200);xMet.setRange("sig",120,130);
  RooDataHist dataMet_bg("dataMet_bg","dataset_bg",xMet_bg,ggInvarMassMET);
  RooRealVar Bern1_bg("Bern1_bg","Berstein 1",14.5,8.,22.);
  RooRealVar Bern2_bg("Bern2_bg","Berstein 2",5.,1.,9.);
  RooRealVar Bern3_bg("Bern3_bg","Berstein 3",4.,0.,8.);
  RooRealVar Bern4_bg("Bern4_bg","Berstein 4",2.,0.,4.);
  RooRealVar Bern5_bg("Bern5_bg","Berstein 5",5.,0.,10.);
  //RooBernstein Bern_bg("Bern","4th order Bernstein Polynomial",xMet_bg,RooArgList(Bern1_bg,Bern2_bg,Bern3_bg,Bern4_bg,Bern5_bg));
  RooPolynomial Bern_bg("Bern","4th order Bernstein Polynomial",xMet_bg,RooArgList(Pol1,Pol2,Pol3,Pol4));
  RooRealVar BernYield_bg("bkgd yield_bg","bkgd yield_bg",70000,0,170000);
  //RooAddPdf MetPdf_bg("MetPdf_bg","MetPdf_bg",RooArgList(Bern_bg),RooArgList(BernYield_bg));
  RooFitResult *rNoMetCut_bg = Bern_bg.chi2FitTo(dataMet_bg/*,Extended(kTRUE)*/,Save(),Range("sb_lo,sb_hi"));
  // RooFitResult *rNoMetCut_bg_full = MetPdf_bg.fitTo(dataMet,Extended(kTRUE),Save());
  // Print fit results 
  cout << "result of fit on all data " << endl ;
  //rNoMetCut_bg_full->Print() ;  
  cout << "result of fit in in background region" << endl ;
  //rNoMetCut_bg->Print() ;
  cout << "result of fit in background region" << endl ;
  RooPlot *xframeMet_bg = xMet_bg.frame(Title("4th Order Polynomial Background, gg with no MET cut"));
  Bern_bg.paramOn(xframeMet_bg, Format("NE",AutoPrecision(2)),Layout(0.6,0.995,0.9) );
  dataMet_bg.plotOn(xframeMet_bg,LineColor(kBlack),MarkerStyle(20),MarkerSize(0.3));
  //Bern_bg.plotOn(xframeMet_bg,Components(Bern_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg,3,kTRUE),FillColor(kViolet));
  Bern_bg.plotOn(xframeMet_bg,Components(Bern_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg,2,kTRUE),FillColor(kGreen));
  Bern_bg.plotOn(xframeMet_bg,Components(Bern_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg,1,kTRUE),FillColor(kOrange));
  //Bern_bg.plotOn(xframeMet_bg,Components(Bern),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg_full,3,kFALSE),FillColor(kViolet));
  //Bern_bg.plotOn(xframeMet_bg,Components(Bern),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg_full,2,kFALSE),FillColor(kGreen));
  //Bern_bg.plotOn(xframeMet_bg,Components(Bern),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg_full,1,kFALSE),FillColor(kOrange));
  Bern_bg.plotOn(xframeMet_bg,Components(Bern_bg),LineColor(kRed),LineStyle(kDashed),Range(95,180));
  Bern_bg.plotOn(xframeMet_bg,Components(Bern_bg),LineColor(kRed));
  dataMet_bg.plotOn(xframeMet_bg,LineColor(kBlack),MarkerStyle(20),MarkerSize(0.3));
  //xframeMet_bg->SetAxisRange(106,149,"X");
  //xframeMet_bg->SetAxisRange(0,1900,"Y");
  //xframeMet_bg->SetMaximum(2000);
  xframeMet_bg->Draw();
  //line125->DrawLine(125.3,0,125.3,xframeMet_bg->GetMaximum());
  StackMet->Draw("histoSAME");
  xframeMet_bg->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithNoMETcut_bgOnly.png");
  
  //SMhiggs, MET>30 SM higgs

  TH1F* SM30Higgs = (TH1F*)ggHggMET30->Clone();
  SM30Higgs->Add(WZHggMET30);
  SM30Higgs->Add(VBFHggMET30);
  SM30Higgs->Add(TTHggMET30);
  RooRealVar xSM30("xSM30","m_{#gamma#gamma}",110,140,"GeV");xSM30.setRange("sigSM30range",115,135);
  RooDataHist dataSM30("dataSM30","data SM30 higgs",xSM30,SM30Higgs);
  //Gaussian for peak
  RooRealVar meanSM30("meanSM30","meanSM30",125/*.3,122,128*/);
  RooRealVar sigmaSM30("sigmaSM30","sigmaSM30",0.06,0,2);
  RooRealVar sigmaSM302("sigmaSM302","sigmaSM302",0.,6);
  //Voigtian (Breit-Vigner x Gaussian) for peak
  RooRealVar widthSM30("widthSM30","widthSM30",1.9,0.1,3.7);
  RooRealVar SM30Yield("sig1 yield","sig1 yield",20000,0,400000);
  RooRealVar SM30Yield2("sig2 yield","sig2 yield",20000,0,400000);
  RooGaussian sigSM302("sigSM302","gaussSM302",xSM30,meanSM30,sigmaSM302);
  //RooGaussian sigSM30 ("sigSM30", "gaussSM30", xSM30,meanSM30,sigmaSM30);
  RooVoigtian sigSM30 ("sigSM30", "gaussSM30", xSM30,meanSM30,widthSM30,sigmaSM30);
  //RooAddPdf SM30Pdf("SM30Pdf","SM30Pdf",RooArgList(sigSM30,sigSM302),RooArgList(SM30Yield,SM30Yield2));
  RooAddPdf SM30Pdf("SM30Pdf","SM30Pdf",RooArgList(sigSM30),RooArgList(SM30Yield));
  RooFitResult *rSM30 = SM30Pdf.fitTo(dataSM30,Extended(kTRUE),Save(),Range("sigSM30range"));
  RooPlot *xframeSM30 = xSM30.frame(Title("Breit-Wigner#otimesGaussian signal, gg SM30 higgs"));
  SM30Pdf.paramOn(xframeSM30, Format("NE",AutoPrecision(2)),Layout(0.6,0.995,0.9) );
  dataSM30.plotOn(xframeSM30,LineColor(kBlack));
  SM30Pdf.plotOn(xframeSM30,LineColor(kBlue));
  SM30Pdf.plotOn(xframeSM30,LineColor(kBlue),VisualizeError(*rSM30,2,kTRUE),FillColor(kGreen));
  SM30Pdf.plotOn(xframeSM30,LineColor(kBlue),VisualizeError(*rSM30,1,kTRUE),FillColor(kOrange));
  //SM30Pdf.plotOn(xframeSM30,Components(sigSM30),LineColor(kRed),LineStyle(kDashed));
  //SM30Pdf.plotOn(xframeSM30,Components(sigSM302),LineColor(kGreen),LineStyle(kDashed));
  dataSM30.plotOn(xframeSM30,LineColor(kBlack));
  xframeSM30->Draw();
  THStack* StackMet30 = new THStack("StackMet30","");
  StackMet30->Add(TTHggMET30);//->Draw("SAME");
  StackMet30->Add(WZHggMET30);//->Draw("SAME"); 
  StackMet30->Add(VBFHggMET30);//->Draw("SAME");
  StackMet30->Add(ggHggMET30);//->Draw("SAME");
  StackMet30->Draw("histoSAMES");
  xframeSM30->Draw("SAMES");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassSMhiggsMET30cut.png");

  //data, MET>30

  RooRealVar xMet30("xMet30","m_{#gamma#gamma}",95,200,"GeV");
  xMet30.setRange("full",95,200);
  RooDataHist dataMet30("dataMet30","dataset",xMet30,ggInvarMassMET30);
  //Gaussian for peak
  RooRealVar gmeanMet30("gmeanMet30","gmean",125.3,122,128);
  RooRealVar gsigmaMet30("gsigmaMet30","gsigma",1.6/*,0.1,5*/);
  //RooGaussian sigMet30("sigMet30","gauss",xMet30,gmeanMet30,gsigmaMet30);
  //Crystal Ball for signal
  RooRealVar meanMet30("cb_mean", "mean" , 125.3, 122, 128.) ;
  RooRealVar sigmaMet30("cb_sigma", "sigma",1.6/*2.6, 0., 3.*/);
  RooRealVar nMet30("n","n", 10.,-5.,25.);
  RooRealVar alphaMet30("alpha","alpha",25.,0.,50.);
  //Voigtian (Breit-Vigner x Gaussian) for peak
  RooRealVar VmeanMet30("Vmean", "Vmean" , 125.3,122.,128.) ;
  RooRealVar VsigmaMet30("VsigmaMet30","VsigmaMet30",0.06);
  RooRealVar VwidthMet30("VwidthMet30","VwidthMet30",1.9);
  RooVoigtian sigMet30 ("sigMet30", "gaussMet30", xMet30,VmeanMet30,VwidthMet30,VsigmaMet30);
  //RooCBShape sigMet30("sigMet30", "crystal ball", xMet30, meanMet30, sigmaMet30, alphaMet30, nMet30);
  RooRealVar sigMet30Yield("signal yield","signal yield",100,0,600);
  RooRealVar Bern301("Bern301","Berstein 1",14.5,8.,22.);
  RooRealVar Bern302("Bern302","Berstein 2",5.,1.,9.);
  RooRealVar Bern303("Bern303","Berstein 3",4.,0.,8.);
  RooRealVar Bern304("Bern304","Berstein 4",2.,0.,4.);
  RooRealVar Bern305("Bern305","Berstein 5",5.,0.,10.);
  RooRealVar Pol301("Pol301","Pol301",0,-.1,.1);
  RooRealVar Pol302("Pol302","Pol302",0,-.001,.001);
  RooRealVar Pol303("Pol303","Pol303",0,-.001,.001);
  RooRealVar Pol304("Pol304","Pol304",0,-.001,.001);
  RooRealVar Pol305("Pol305","Pol305",0,-.0001,.0001);
  //RooBern30stein Bern30("Bern30","4th order Bern30stein Polynomial",xMet30,RooArgList(Bern301,Bern302,Bern303,Bern304,Bern305));
  RooPolynomial Bern30("Bern30","4th Order Polynomial",xMet30,RooArgList(Pol301,Pol302,Pol303,Pol304));
  RooRealVar Bern30Yield("bkgd yield","bkgd yield",20000,0,400000);
  RooAddPdf Met30Pdf("Met30Pdf","Met30Pdf",RooArgList(Bern30,sigMet30),RooArgList(Bern30Yield,sigMet30Yield));
  RooFitResult *rMet30Cut = Met30Pdf.fitTo(dataMet30,Extended(kTRUE),Save());
  RooPlot *xframeMet30 = xMet30.frame(Title("Voigtian Signal, 4th Order Polynomial Background, gg with MET>30 cut"));
  Met30Pdf.paramOn(xframeMet30, Format("NE",AutoPrecision(2)),Layout(0.6,0.995,0.9) );
  dataMet30.plotOn(xframeMet30,LineColor(kBlack),MarkerSize(0.3),MarkerStyle(20));
  //Met30Pdf.plotOn(xframeMet30,Components(Bern30),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rMet30Cut,3,kTRUE),FillColor(kViolet));
  Met30Pdf.plotOn(xframeMet30,Components(Bern30),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rMet30Cut,2,kTRUE),FillColor(kGreen));
  Met30Pdf.plotOn(xframeMet30,Components(Bern30),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rMet30Cut,1,kTRUE),FillColor(kOrange));
  //Met30Pdf.plotOn(xframeMet30,LineColor(kBlue),LineStyle(kDashed),VisualizeError(*rMet30Cut,1,kTRUE),FillColor(kOrange));
  dataMet30.plotOn(xframeMet30,LineColor(kBlack),MarkerSize(0.3),MarkerStyle(20));
  Met30Pdf.plotOn(xframeMet30,Components(Bern30),LineColor(kRed),LineStyle(kDashed),Range("full"));
  Met30Pdf.plotOn(xframeMet30,LineColor(kBlue));
  Met30Pdf.plotOn(xframeMet30,Components(sigMet30),LineColor(kBlue));
  dataMet30.plotOn(xframeMet30,LineColor(kBlack),MarkerSize(0.3),MarkerStyle(20));
  xframeMet30->SetAxisRange(106,149,"X");
  xframeMet30->SetAxisRange(0,1000,"Y");
  xframeMet30->Draw();
  //line125->DrawLine(125.3,0,125.3,xframeMet30->GetMaximum());
  StackMet30->Draw("histoSAME");
  xframeMet30->Draw("SAME");
  legHiggs->Draw();
  Double_t sigSM30yield=0.;
  Double_t sigSM30Err=0.;
  ggHggMET30->Add(WZHggMET30);ggHggMET30->Add(VBFHggMET30);ggHggMET30->Add(TTHggMET30);
  sigSM30yield = ggHggMET30->IntegralAndError(0,999,sigSM30Err,"");
  //sigSM30yield/=10;sigSM30Err/=10;
  Double_t sig30Yield=sigMet30Yield.getVal();
  Double_t sig30YieldErr=sigMet30Yield.getError();
  //cout<<"sigSM30yield: "<<sigSM30yield<<"  sigSM30Err: "<<sigSM30Err<<endl;
  //cout<<"sig30Yield: "<<sig30Yield<<"  sig30YieldErr: "<<sig30YieldErr<<endl;
  //sigStrengthVsMet30->SetBinContent(1,sig30Yield/sigSM30yield);
  //sigStrengthVsMet30->SetBinError(1,sqrt(sig30YieldErr*sig30YieldErr/(sigSM30yield*sigSM30yield) + (sig30Yield*sig30Yield*sigSM30Err*sigSM30Err)/(sigSM30yield*sigSM30yield)));
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET30cut.png");
  
  RooRealVar xMet30_bg("xMet30_bg","m_{#gamma#gamma}",95,200,"GeV");
  xMet30_bg.setRange("sb_lo",95,120);xMet30_bg.setRange("sb_hi",130,200);xMet30.setRange("full",95,200);xMet30.setRange("sig",120,130);
  RooDataHist dataMet30_bg("dataMet30_bg","dataset_bg",xMet30_bg,ggInvarMassMET30);
  RooRealVar Bern301_bg("Bern301_bg","Berstein 1",14.5,8.,22.);
  RooRealVar Bern302_bg("Bern302_bg","Berstein 2",5.,1.,9.);
  RooRealVar Bern303_bg("Bern303_bg","Berstein 3",4.,0.,8.);
  RooRealVar Bern304_bg("Bern304_bg","Berstein 4",2.,0.,4.);
  RooRealVar Bern305_bg("Bern305_bg","Berstein 5",5.,0.,10.);
  //RooBern30stein Bern30_bg("Bern30","4th order Bern30stein Polynomial",xMet30_bg,RooArgList(Bern301_bg,Bern302_bg,Bern303_bg,Bern304_bg,Bern305_bg));
  RooPolynomial Bern30_bg("Bern30","4th order Bern30stein Polynomial",xMet30_bg,RooArgList(Pol301,Pol302,Pol303,Pol304));
  RooRealVar Bern30Yield_bg("bkgd yield_bg","bkgd yield_bg",70000,0,170000);
  //RooAddPdf Met30Pdf_bg("Met30Pdf_bg","Met30Pdf_bg",RooArgList(Bern30_bg),RooArgList(Bern30Yield_bg));
  RooFitResult *rMet30Cut_bg = Bern30_bg.chi2FitTo(dataMet30_bg/*,Extended(kTRUE)*/,Save(),Range("sb_lo,sb_hi"));
  // RooFitResult *rMet30Cut_bg_full = Met30Pdf_bg.fitTo(dataMet30,Extended(kTRUE),Save());
  // Print fit results 
  cout << "result of fit on all data " << endl ;
  //rMet30Cut_bg_full->Print() ;  
  cout << "result of fit in in background region" << endl ;
  //rMet30Cut_bg->Print() ;
  cout << "result of fit in background region" << endl ;
  RooPlot *xframeMet30_bg = xMet30_bg.frame(Title("4th Order Polynomial Background, gg with MET>30 cut"));
  Bern30_bg.paramOn(xframeMet30_bg, Format("NE",AutoPrecision(2)),Layout(0.6,0.995,0.9) );
  dataMet30_bg.plotOn(xframeMet30_bg,LineColor(kBlack),MarkerStyle(20),MarkerSize(0.3));
  //Bern30_bg.plotOn(xframeMet30_bg,Components(Bern30_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rMet30Cut_bg,3,kTRUE),FillColor(kViolet));
  Bern30_bg.plotOn(xframeMet30_bg,Components(Bern30_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rMet30Cut_bg,2,kTRUE),FillColor(kGreen));
  Bern30_bg.plotOn(xframeMet30_bg,Components(Bern30_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rMet30Cut_bg,1,kTRUE),FillColor(kOrange));
  //Bern30_bg.plotOn(xframeMet30_bg,Components(Bern30),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rMet30Cut_bg_full,3,kFALSE),FillColor(kViolet));
  //Bern30_bg.plotOn(xframeMet30_bg,Components(Bern30),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rMet30Cut_bg_full,2,kFALSE),FillColor(kGreen));
  //Bern30_bg.plotOn(xframeMet30_bg,Components(Bern30),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rMet30Cut_bg_full,1,kFALSE),FillColor(kOrange));
  Bern30_bg.plotOn(xframeMet30_bg,Components(Bern30_bg),LineColor(kRed),LineStyle(kDashed),Range(95,180));
  Bern30_bg.plotOn(xframeMet30_bg,Components(Bern30_bg),LineColor(kRed));
  dataMet30_bg.plotOn(xframeMet30_bg,LineColor(kBlack),MarkerStyle(20),MarkerSize(0.3));
  xframeMet30_bg->SetAxisRange(106,149,"X");
  xframeMet30_bg->SetAxisRange(0,1000,"Y");
  //xframeMet30_bg->SetMaximum(2000);
  xframeMet30_bg->Draw();
  //line125->DrawLine(125.3,0,125.3,xframeMet30_bg->GetMaximum());
  StackMet30->Draw("histoSAME");
  xframeMet30_bg->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET30cut_bgOnly.png");


  fin.Close();
  // return;    
}

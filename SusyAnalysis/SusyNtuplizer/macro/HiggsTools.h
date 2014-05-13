#ifndef HiggsTools_h
#define HiggsTools_h

#include <iostream>
#include <TLegend.h>
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH3.h"
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
#include "RooChebychev.h"
#include "RooExtendPdf.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooNumConvPdf.h"
#include "RooBreitWigner.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooBernstein.h"
#include "RooProduct.h"
#include "TProfile.h"
#include "TH2.h"
#include "TH2F.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TString.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "PhysicsTools/TagAndProbe/interface/RooCMSShape.h"
#include "THStack.h"
#include "TLatex.h"
#include "TMatrixDSym.h"


class HiggsTools {

 public:

  HiggsTools();
  ~HiggsTools();
  void Loop();
  void MakeBlindInvMass(TH1D *&hInvMass);
  void MakeBlindMet(TH1F *&hMet);
  void MakeBlindMet30(TH1F *&hMet);
 
 private:

};

#endif

#ifdef HiggsTools_cxx
HiggsTools::HiggsTools()
{
 
}

HiggsTools::~HiggsTools()
{

}
/*
void HiggsTools::Loop()
{

}
*/
void HiggsTools::MakeBlindInvMass(TH1D *&hInvMass){
  int bin=hInvMass->FindBin(118);
  while(bin<=hInvMass->FindBin(133)-1){
    hInvMass->SetBinContent(bin,0);
    hInvMass->SetBinError(bin,0);
    bin++;
  }
}

void HiggsTools::MakeBlindMet(TH1F *&hMet){
  int bin=hMet->FindBin(50);
  while(bin<hMet->GetNbinsX()+1){
    hMet->SetBinContent(bin,0);
    hMet->SetBinError(bin,0);
    bin++;
  }}

void HiggsTools::MakeBlindMet30(TH1F *&hMet){
  int bin=hMet->FindBin(30);
  while(bin<hMet->GetNbinsX()+1){
    hMet->SetBinContent(bin,0);
    hMet->SetBinError(bin,0);
    bin++;
  }
}

#endif //#ifdef HiggsTools_cxx

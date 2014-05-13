#define HiggsTools_cxx
#include "HiggsTools.h"

#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TRandom.h"
#include "TRandom3.h"

using namespace std;
using namespace RooFit;
using namespace ROOT::Math;

const double alpha = 1 - 0.6827;
float L_int = 19499.;//full dataset
Double_t xbins[]={0,5,10,15,20,25,30,35,40,50,65,80,100,150,250};
int NmetBins=(sizeof(xbins)/sizeof(Double_t))-1;
float metPlotXmax = xbins[NmetBins];
Double_t xbinsWide[]={0,10,20,30,40,50,60,80,100,280};
int NmetBinsWide=(sizeof(xbinsWide)/sizeof(Double_t))-1;
float metPlotXmaxWide = xbinsWide[NmetBinsWide];
Double_t xbinsXtraWide[]={0,15,30,50,180};
//Double_t xbinsXtraWide[]={0,10,20,30,50,200};
//Double_t xbinsXtraWide[]={0,20,30,40,60,200};//SUS-13-014 binning
int NmetBinsXtraWide=(sizeof(xbinsXtraWide)/sizeof(Double_t))-1;
float metPlotXmaxXtraWide = xbinsXtraWide[NmetBinsXtraWide];
Double_t MTbins[]={0,30,60,90,180};
int nMTbins=(sizeof(MTbins)/sizeof(Double_t))-1;
float MTplotXmax = MTbins[nMTbins];
double sigLo=120,sigHi=131;
double sbFitLoLo=103,sbFitHiHi=163;
double sbLoLo=103,sbLoHi=118;
double sbHiLo=133,sbHiHi=163;
bool reject=true;

Double_t fpow(Double_t *x, Double_t *par){
  //Fit function definition
  if(reject && x[0]>=sbLoHi && x[0]<=sbHiLo){
    TF1::RejectPoint();
    return 0;
  }
  return par[0]*pow(x[0],par[1]);
}
Double_t fexp(Double_t *x, Double_t *par){
  //Fit function definition
  if(reject && x[0]>=sbLoHi && x[0]<=sbHiLo){
    TF1::RejectPoint();
    return 0;
  }
  return exp(par[0]+par[1]*x[0]);
}
Double_t flin(Double_t *x, Double_t *par){
  //Fit function definition
  if(reject && x[0]>=sbLoHi && x[0]<=sbHiLo){
    TF1::RejectPoint();
    return 0;
  }
  return par[0]+par[1]*x[0];
}
Double_t fpol3(Double_t *x, Double_t *par){
  //Fit function definition
  if(reject && x[0]>=sbLoHi && x[0]<=sbHiLo){
    TF1::RejectPoint();
    return 0;
  }
  return par[0]*x[0]+par[1]*x[0]*x[0]+par[2]*x[0]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]*x[0];
}
Double_t fpol5(Double_t *x, Double_t *par){
  //Fit function definition
  if(reject && x[0]>=sbLoHi && x[0]<=sbHiLo){
    TF1::RejectPoint();
    return 0;
  }
  return par[0]*x[0]+par[1]*x[0]*x[0]+par[2]*x[0]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]*x[0]*x[0];
}
Double_t TotErr(Double_t binVal, Double_t binL, Double_t binLerr, Double_t binU, Double_t binUerr, Double_t totL, Double_t totLerr, Double_t totU, Double_t totUerr, Double_t fitInt, Double_t fitIntErr, Double_t fitIntSyst, Double_t& StatErr, Double_t& FitSystErr, Double_t& FitShapeSystErr, Double_t lValScale, Double_t uValScale, Double_t& halfDiffErr){
  /*
  Double_t L = binLerr*binLerr/totL/totL - totLerr*totLerr*binL*binL/totL/totL/totL/totL;
  Double_t U = binUerr*binUerr/totU/totU - totUerr*totUerr*binU*binU/totU/totU/totU/totU;
  */
  Double_t L = binL/totL/totL - binL*binL/totL/totL/totL;
  Double_t U = binU/totU/totU - binU*binU/totU/totU/totU;
  Double_t LplusU = binL/totL + binU/totU; 
  Double_t LplusU2 = LplusU*LplusU;
  //Double_t HalfDiff2 = (binL-binU)*(binL-binU)/2./2.;
  Double_t HalfDiff2 = (lValScale-uValScale)*(lValScale-uValScale)/2./2.;
  halfDiffErr = (binVal==0. ? 0. : sqrt(HalfDiff2));
  Double_t StatErr2=fabs((1./4.)*((U+L)*fitInt*fitInt));
  StatErr=(binVal==0. ? 0. : sqrt(StatErr2));
  Double_t FitShapeSystErr2=fabs((fitInt-fitIntSyst)/2.);
  Double_t fitErrTot2 = fitIntErr*fitIntErr;// + FitShapeSystErr2;
  Double_t FitStatSystErr2=fabs((1./4.)*(LplusU2*fitErrTot2));
  FitSystErr=(binVal==0. ? 0. : sqrt(FitStatSystErr2));
  FitShapeSystErr=0.;//(binVal==0. ? 0. : sqrt(FitShapeSystErr2));
  //Double_t FitShapeSystErr2=fabs((fitInt-fitIntSyst)/2.);
  //FitShapeSystErr=sqrt(FitShapeSystErr2);
  //Double_t sigma2 = StatErr2/* + FitStatSystErr2*/ + HalfDiff2;//for limits take out fitstatsysterr2, keep in for plot
  Double_t sigma2 = StatErr2 + FitStatSystErr2 + HalfDiff2;//for limits take out fitstatsysterr2, keep in for plot
  Double_t sigma = (binVal==0. ? 0. : sqrt(sigma2));
  return sigma;

}

void DivideBy30gev(TH1F* &hist){
  /*for(int i=0;i<=hist->GetNbinsX();i++){
    float val = hist->GetBinContent(i), err=hist->GetBinError(i);
    hist->SetBinContent(i,val/hist->GetBinWidth(i));
    hist->SetBinError(i,err/hist->GetBinWidth(i));
    }*/
  float val = hist->GetBinContent(4), err=hist->GetBinError(4);
  hist->SetBinContent(4,val/3.);
  hist->SetBinError(4,err/3.);
  return;
}
void DivideBy30gev(TH1D* &hist){
  /* for(int i=0;i<=hist->GetNbinsX();i++){
     float val = hist->GetBinContent(i), err=hist->GetBinError(i);
     hist->SetBinContent(i,val/hist->GetBinWidth(i));
     hist->SetBinError(i,err/hist->GetBinWidth(i));
     }*/
  float val = hist->GetBinContent(4), err=hist->GetBinError(4);
  hist->SetBinContent(4,val/3.);
  hist->SetBinError(4,err/3.);
  return;
}
void MultBy30gev(TH1F* &hist){
  /* for(int i=0;i<=hist->GetNbinsX();i++){
     float val = hist->GetBinContent(i), err=hist->GetBinError(i);
     //hist->SetBinContent(i,val*hist->GetBinWidth(i));
     //hist->SetBinError(i,err*hist->GetBinWidth(i));
     }*/
  float val = hist->GetBinContent(4), err=hist->GetBinError(4);
  hist->SetBinContent(4,val*3.);
  hist->SetBinError(4,err*3.);
  return;
}

void DivideBy15gev(TH1F* &hist){
  for(int i=0;i<=hist->GetNbinsX();i++){
    float val = hist->GetBinContent(i), err=hist->GetBinError(i), w = hist->GetBinWidth(i);
    hist->SetBinContent(i,val*(15./w));
    hist->SetBinError(i,err*(15./w));
  }
  return;
}
void DivideBy15gev(TH1D* &hist){
  for(int i=0;i<=hist->GetNbinsX();i++){
    float val = hist->GetBinContent(i), err=hist->GetBinError(i), w = hist->GetBinWidth(i);
    hist->SetBinContent(i,val*(15./w));
    hist->SetBinError(i,err*(15./w));
  }
  return;
}
void MultBy15gev(TH1F* &hist){
  for(int i=0;i<=hist->GetNbinsX();i++){
    float val = hist->GetBinContent(i), err=hist->GetBinError(i), w = hist->GetBinWidth(i);
    hist->SetBinContent(i,val/(15./w));
    hist->SetBinError(i,err/(15./w));
  }
  return;
}
void DivideByBinWidth(TH1F* &hist){
  for(int i=0;i<=hist->GetNbinsX();i++){
    float val = hist->GetBinContent(i), err=hist->GetBinError(i);
    hist->SetBinContent(i,val/hist->GetBinWidth(i));
    hist->SetBinError(i,err/hist->GetBinWidth(i));
  }
  return;
}
void DivideByBinWidth(TH1D* &hist){
  for(int i=0;i<=hist->GetNbinsX();i++){
    float val = hist->GetBinContent(i), err=hist->GetBinError(i);
    hist->SetBinContent(i,val/hist->GetBinWidth(i));
    hist->SetBinError(i,err/hist->GetBinWidth(i));
  }
  return;
}

void MultByBinWidth(TH1F* &hist){
  for(int i=0;i<=hist->GetNbinsX();i++){
    float val = hist->GetBinContent(i), err=hist->GetBinError(i);
    //hist->SetBinContent(i,val*hist->GetBinWidth(i));
    //hist->SetBinError(i,err*hist->GetBinWidth(i));
  }
  return;
}


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
void AddOverflowToLastBin(TH1D* &hist){
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
void AddHiggsSigmaUncert(TH1F* &hist){
  for(int i=0;i<=hist->GetNbinsX();i++){
    float val = hist->GetBinError(i)*1.13;
    hist->SetBinError(i,val);
  }
}
void ScaleHiggs(TH1F* &hist, TH1F* histNoSF, float xsec, float nEvents){

  float a = L_int, aErr=(/*.024*/0.*L_int);
  float b = xsec, bErr = (/*.13*/0.*xsec);
  float c = 2.29e-03;float cErr = (0.*c);//h->gg bf
  float a2=a*a, aErr2=aErr*aErr,b2=b*b, bErr2=bErr*bErr,c2=c*c, cErr2=cErr*cErr;
  float d = 0., dErr=0.;
  float phoSF2 = 0.994*0.994;
  float norm = histNoSF->GetEntries()/histNoSF->Integral(0,-1);

  for(int i=0;i<=hist->GetNbinsX();i++){
    d = hist->GetBinContent(i);
    dErr = hist->GetBinError(i);
    float d2=d*d, dErr2=dErr*dErr;
    float val = norm*phoSF2*a*b*c*d/nEvents;
    float err = sqrt(b*c*d*aErr*b*c*d*aErr + a*c*d*bErr*a*c*d*bErr + a*b*d*cErr*a*b*d*cErr + a*b*c*dErr*a*b*c*dErr)/nEvents;
    hist->SetBinError(i,err);
    hist->SetBinContent(i,val);
  }




}
void HiggsTools::Loop(){
  gStyle->SetOptStat(0);
  gStyle->SetTitleFont(42,"t");
  //gStyle->SetErrorX(0);
  //gStyle->SetEndError(0);
  //TFile fin("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_15GevCHfakeCut_EleVeto_NewDiEMPtBins_19499pb_NewJetMatching_2JetReq_May16.root","READ");
  //TFile fin("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_15GevCHfakeCut_PixelVetoOnFakes_NewDiEMPtBins_19499pb_NewJetMatching_2JetReq_PUjetID_NewDoGG_Sep16_allTrigs_bEta2p0.root","READ");
  //TFile *fin = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/hist_HiggsAna_Data2012A_13July_06Aug_recover_Data2012B_13July_Data2012C_Prompt_Runs198022-198903_Runs198941-203742_Data2012D_Prompt_Runs207920-209151_Runs203777-207905_Filter_Ana_HLT_JSON_Two40-25GeVbarrelPhotons_ALL-Higgsino-selectedEvents.root","READ");
  TFile *fin = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Data2012A_13July_06Aug_recover_Data2012B_13July_Data2012C_Prompt_Runs198022-198903_Runs198941-203742_Data2012D_Prompt_Runs207920-209151_Runs203777-207905_Filter_Ana_HLT_JSON_Two40-25GeVbarrelPhotons_ALL-Higgsino-selectedEvents.root","READ");
  // TFile fin("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Higgs_cms538vop1_WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_May7.root","READ");
  /*
    TFile f_ggHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Higgs_cms533v1_GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Jul31.root","READ");
    TFile f_WZHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Higgs_cms533v1_WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Jul31.root","READ");
    TFile f_TTHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Higgs_cms533v1_TTH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Jul31.root","READ");
    TFile f_VBFHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Higgs_cms533v1_VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Jul31.root","READ");
  */
  /*
  TFile f_ggHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Higgs_cms533v1_GluGluToHToGG_M-126_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");
  TFile f_WZHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Higgs_cms533v1_WH_ZH_HToGG_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");
  TFile f_TTHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Higgs_cms533v1_TTH_HToGG_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");
  TFile f_VBFHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Higgs_cms533v1_VBF_HToGG_M-126_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");
  */
  TFile *f_ggHgg = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Higgs_cms533v1_GluGluToHToGG_M-126_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");
  TFile *f_WZHgg = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Higgs_cms533v1_WH_ZH_HToGG_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");
  TFile *f_TTHgg = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Higgs_cms533v1_TTH_HToGG_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");
  TFile *f_VBFHgg = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Higgs_cms533v1_VBF_HToGG_M-126_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");
  /*
  TFile f_aaW130("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_cms533_WinoNLSP_chargino130_bino1_5_10_15_hw_aaw_full_Sep11.root","READ");
  TFile f_aaW275("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_cms533_WinoNLSP_chargino275_bino1_5_hw_aaw_full_Sep11.root","READ");
  */
  /*
  TFile f_SMS_WH("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_SMS_TChiWH_WincHgg_2J_Redo.root","READ");
  TFile f_SMS_ZH("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_SMS_TChiZH_ZincHgg_2J.root","READ");
  TFile f_SMS_HH_2b2g("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_SMS_TChiHH_2b2g_2J.root","READ");
  TFile f_SMS_HH_2W2g("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_SMS_TChiHH_2W2g_2J.root","READ");
  TFile f_SMS_HH_2Z2g("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_HiggsAna_Photon_SMS_TChiHH_2Z2g_2J.root","READ");
  */

  TFile *f_SMS_WH = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Photon_SMS_TChiWH_WincHgg_2J.root","READ");
  TFile *f_SMS_WH_InvMass = TFile::Open("hist_HiggsAna_Photon_SMS_TChiWH_WincHgg_2J_InvMass.root","READ");
  TFile *f_SMS_HH2W_InvMass = TFile::Open("hist_HiggsAna_Photon_SMS_TChiHH_2W2g_2J_InvMass.root","READ");
  TFile *f_SMS_HH2Z_InvMass = TFile::Open("hist_HiggsAna_Photon_SMS_TChiHH_2Z2g_2J_InvMass.root","READ");
  TFile *f_SMS_HH2tau_InvMass = TFile::Open("hist_HiggsAna_Photon_SMS_TChiHH_2tau2g_2J_InvMass.root","READ");
  TFile *f_SMS_ZH = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Photon_SMS_TChiZH_ZincHgg_2J.root","READ");
  TFile *f_SMS_HH_2b2g = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Photon_SMS_TChiHH_2b2g_2J.root","READ");
  TFile *f_SMS_HH_2W2g = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Photon_SMS_TChiHH_2W2g_2J.root","READ");
  TFile *f_SMS_HH_2Z2g = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Photon_SMS_TChiHH_2Z2g_2J.root","READ");
  TFile *f_SMS_HH_2tau2g = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Photon_SMS_TChiHH_2tau2g_2J.root","READ");

  TFile fout("Data_and_DataDrivenBkg.root","RECREATE");
  TFile fout_We("Data_and_DataDrivenBkg_ggEle.root","RECREATE");
  TFile fout_Wmu("Data_and_DataDrivenBkg_ggMu.root","RECREATE");
  TFile fout_SMHiggs("SMHiggsBkg.root","RECREATE");

  TCanvas *c1 = new TCanvas("c1","",800,750);
  c1->SetMargin(0.15,0.05,0.15,0.1);
  c1->cd();
  TCanvas *c1_rat = new TCanvas("c1_rat","",800,250);
  c1_rat->SetMargin(0.15,0.05,0.15,0.1);
  c1_rat->cd();
  TCanvas *c2 = new TCanvas("c2","",800,900);
  //c2->Divide(1,2,0,.05);  
  /*
  TPad *p1 = new TPad("p1","",0.,.3,1.,1.);
  p1->SetBottomMargin(0);
  p1->Draw();
  TPad *p2 = new TPad("p2","",0.,.125,1.,.25);
  p2->SetTopMargin(0);
  p2->Draw();
  */
  /*
  TPad *p1 = new TPad("p1","",0.,.08,1.,0.78);
  p1->SetTopMargin(0.02);
  p1->Draw();
  TPad *p2 = new TPad("p2","",0.,.78,1.,.99);
  p2->SetBottomMargin(0.05);
  p2->Draw();
  */

  TPad *p1 = new TPad("p1","",0.,.29,1.,1.);
  p1->SetMargin(0.135,0.025,0.02,0.07);
  //p1->SetBottomMargin(0.02);
  p1->Draw();
  TPad *p2 = new TPad("p2","",0.,.08,1.,.29);
  p2->SetMargin(0.135,0.025,0.115,0.06);
  //p2->SetBottomMargin(0.02);
  p2->Draw();
  TPad *p3 = new TPad("p3","",0.,0.,1.,.06);
  p3->SetMargin(0.135,0.025,0.,0.);
  //p3->SetBottomMargin(0);
  p3->Draw();
  c1->cd();

  

  c1->cd();
  TH1F* PhoEleDR_data = (TH1F*)fin->Get("dRgEle_Loose");
  TH1F* PhoEleDR_WZHgg = (TH1F*)f_SMS_WH->Get("dR_HiggsG_WZele_GenMatch");
  c1->SetLogy(1);
  PhoEleDR_data->SetTitle("");PhoEleDR_data->GetXaxis()->SetTitle("#DeltaR_{e,#gamma}");PhoEleDR_data->GetXaxis()->SetRangeUser(0,3.59);
  PhoEleDR_data->Draw();
  c1->Print("Plots/Higgs/PhoEleDR_data.png");
  c1->Print("Plots/Higgs/PhoEleDR_data.pdf");
  c1->SetLogy(0);
  PhoEleDR_WZHgg->SetTitle("");PhoEleDR_WZHgg->GetXaxis()->SetTitle("#DeltaR_{e,#gamma}");PhoEleDR_WZHgg->GetXaxis()->SetRangeUser(0,4.49);
  PhoEleDR_WZHgg->Draw();
  c1->Print("Plots/Higgs/PhoEleDR_WZHgg.png");
  c1->Print("Plots/Higgs/PhoEleDR_WZHgg.pdf");


  TH1F* h_eleRelIso_nMinus1 = (TH1F*)fin->Get("eleRelIso_nMinus1");
  TH1F* h_eleRelIso_passID = (TH1F*)fin->Get("eleRelIso_passID");
  TH2F* h_eleEtaVsRelIso_passID = (TH2F*)fin->Get("eleEtaVsRelIso_passID");
  TH1F* h_eleRelIso_1Ele_0_1Jets = (TH1F*)fin->Get("eleRelIso_1Ele_0_1Jets");
  c1->SetLogy(0);
  h_eleRelIso_nMinus1->GetXaxis()->SetTitleSize(.04);h_eleRelIso_nMinus1->GetXaxis()->SetTitleOffset(1.5);
  h_eleRelIso_nMinus1->GetXaxis()->SetTitle("#sumIso/p_{T}");h_eleRelIso_nMinus1->GetYaxis()->SetTitle("Events");
  h_eleRelIso_nMinus1->GetXaxis()->SetRangeUser(0,.159);h_eleRelIso_nMinus1->Draw("PE");
  c1->Print("Plots/Higgs/EleRelIso_Nminus1.png");
  c1->Print("Plots/Higgs/EleRelIso_Nminus1.pdf");
  c1->SetLogy(0);
  h_eleRelIso_passID->GetXaxis()->SetTitleSize(.04);h_eleRelIso_passID->GetXaxis()->SetTitleOffset(1.5);
  h_eleRelIso_passID->GetXaxis()->SetTitle("#sumIso/p_{T}");h_eleRelIso_passID->GetYaxis()->SetTitle("Events");
  h_eleRelIso_passID->GetXaxis()->SetRangeUser(0,.159);h_eleRelIso_passID->Draw("PE");
  c1->Print("Plots/Higgs/EleRelIso_passID.png");
  c1->Print("Plots/Higgs/EleRelIso_passID.pdf");
  h_eleRelIso_1Ele_0_1Jets->GetXaxis()->SetTitleSize(.04);h_eleRelIso_1Ele_0_1Jets->GetXaxis()->SetTitleOffset(1.5);
  h_eleRelIso_1Ele_0_1Jets->GetXaxis()->SetTitle("#sumIso/p_{T}");h_eleRelIso_1Ele_0_1Jets->GetYaxis()->SetTitle("Events");
  h_eleRelIso_1Ele_0_1Jets->GetXaxis()->SetRangeUser(0,.159);h_eleRelIso_1Ele_0_1Jets->Draw("PE");
  c1->Print("Plots/Higgs/EleRelIso_1Ele_0_1Jets.png");
  c1->Print("Plots/Higgs/EleRelIso_1Ele_0_1Jets.pdf");
  c1->SetLogy(0);
  h_eleEtaVsRelIso_passID->GetXaxis()->SetTitleSize(.04);h_eleEtaVsRelIso_passID->GetXaxis()->SetTitleOffset(1.5);
  h_eleEtaVsRelIso_passID->GetXaxis()->SetTitle("#sumIso/p_{T}");h_eleEtaVsRelIso_passID->GetYaxis()->SetTitle("| #eta |");
  h_eleEtaVsRelIso_passID->GetXaxis()->SetRangeUser(0,.159);h_eleEtaVsRelIso_passID->GetYaxis()->SetRangeUser(0,2.49);h_eleEtaVsRelIso_passID->Draw("COLZ");
  c1->Print("Plots/Higgs/EleEtaVsRelIso_passID.png");
  c1->Print("Plots/Higgs/EleEtaVsRelIso_passID.pdf");
  
  TH1F* h_eleRelIso_1Ele_0_1Jets_Z = (TH1F*)fin->Get("eleRelIso_1Ele_0_1Jets_Z");
  TH1F* h_eleRelIso_1Ele_0_1Jets_noZ = (TH1F*)fin->Get("eleRelIso_1Ele_0_1Jets_noZ");
  TH1F* h_eleSihih_1Ele_0_1Jets_Z = (TH1F*)fin->Get("eleSihih_1Ele_0_1Jets_Z");
  TH1F* h_eleSihih_1Ele_0_1Jets_noZ = (TH1F*)fin->Get("eleSihih_1Ele_0_1Jets_noZ");
  TH1F* h_eleOneOverEminusOneOverP_1Ele_0_1Jets_Z = (TH1F*)fin->Get("eleOneOverEminusOneOverP_1Ele_0_1Jets_Z");
  TH1F* h_eleOneOverEminusOneOverP_1Ele_0_1Jets_noZ = (TH1F*)fin->Get("eleOneOverEminusOneOverP_1Ele_0_1Jets_noZ");

  h_eleRelIso_1Ele_0_1Jets_Z->SetLineColor(kRed);h_eleRelIso_1Ele_0_1Jets_Z->SetMarkerColor(kRed);
  h_eleSihih_1Ele_0_1Jets_Z->SetLineColor(kRed);h_eleSihih_1Ele_0_1Jets_Z->SetMarkerColor(kRed);
  h_eleOneOverEminusOneOverP_1Ele_0_1Jets_Z->SetLineColor(kRed);h_eleOneOverEminusOneOverP_1Ele_0_1Jets_Z->SetMarkerColor(kRed);

  float eleScale = h_eleRelIso_1Ele_0_1Jets_Z->Integral()/h_eleRelIso_1Ele_0_1Jets_noZ->Integral();
  h_eleRelIso_1Ele_0_1Jets_noZ->Scale(eleScale);
  eleScale = h_eleSihih_1Ele_0_1Jets_Z->Integral()/h_eleSihih_1Ele_0_1Jets_noZ->Integral();
  h_eleSihih_1Ele_0_1Jets_noZ->Scale(eleScale);
  eleScale = h_eleOneOverEminusOneOverP_1Ele_0_1Jets_Z->Integral()/h_eleOneOverEminusOneOverP_1Ele_0_1Jets_noZ->Integral();
  h_eleOneOverEminusOneOverP_1Ele_0_1Jets_noZ->Scale(eleScale);

  h_eleRelIso_1Ele_0_1Jets_Z->Rebin(2);h_eleRelIso_1Ele_0_1Jets_noZ->Rebin(2);
  h_eleRelIso_1Ele_0_1Jets_Z->GetXaxis()->SetTitle("#sumIso/p_{T}");;h_eleRelIso_1Ele_0_1Jets_Z->GetXaxis()->SetRangeUser(0,.149);
  h_eleRelIso_1Ele_0_1Jets_Z->Draw("PE");h_eleRelIso_1Ele_0_1Jets_noZ->Draw("PEsames");
  c1->Print("Plots/Higgs/EleZcomp_RelIso.png");
  c1_rat->cd();
  h_eleRelIso_1Ele_0_1Jets_Z->Divide(h_eleRelIso_1Ele_0_1Jets_noZ);
  h_eleRelIso_1Ele_0_1Jets_Z->Draw("PE");
  TLine eleL(0,1,.15,1);eleL.Draw();
  c1_rat->Print("Plots/Higgs/EleZcomp_RelIso_ratio.png");
  c1->cd();
  h_eleSihih_1Ele_0_1Jets_Z->Rebin(20);h_eleSihih_1Ele_0_1Jets_noZ->Rebin(20);
  h_eleSihih_1Ele_0_1Jets_Z->GetXaxis()->SetTitle("#sigma_{i#etai#eta}");h_eleSihih_1Ele_0_1Jets_Z->GetXaxis()->SetRangeUser(.0061,.0289);
  h_eleSihih_1Ele_0_1Jets_Z->Draw("PE");h_eleSihih_1Ele_0_1Jets_noZ->Draw("PEsames");
  c1->Print("Plots/Higgs/EleZcomp_Sihih.png");
  c1_rat->cd();
  h_eleSihih_1Ele_0_1Jets_Z->Divide(h_eleSihih_1Ele_0_1Jets_noZ);
  h_eleSihih_1Ele_0_1Jets_Z->Draw("PE");eleL.DrawLine(.006,1,.029,1);
  c1_rat->Print("Plots/Higgs/EleZcomp_Sihih_ratio.png");
  c1->cd();
  h_eleOneOverEminusOneOverP_1Ele_0_1Jets_Z->Rebin(40);h_eleOneOverEminusOneOverP_1Ele_0_1Jets_noZ->Rebin(40);
  h_eleOneOverEminusOneOverP_1Ele_0_1Jets_Z->GetXaxis()->SetTitle("1/E - 1/p");h_eleOneOverEminusOneOverP_1Ele_0_1Jets_Z->GetXaxis()->SetRangeUser(0,.049);
  h_eleOneOverEminusOneOverP_1Ele_0_1Jets_Z->Draw("PE");h_eleOneOverEminusOneOverP_1Ele_0_1Jets_noZ->Draw("PEsames");
  c1->Print("Plots/Higgs/EleZcomp_OneOverEminusOneOverP.png");
  c1_rat->cd();
  h_eleOneOverEminusOneOverP_1Ele_0_1Jets_Z->Divide(h_eleOneOverEminusOneOverP_1Ele_0_1Jets_noZ);
  h_eleOneOverEminusOneOverP_1Ele_0_1Jets_Z->Draw("PE");eleL.DrawLine(0,1,.05,1);
  c1_rat->Print("Plots/Higgs/EleZcomp_OneOverEminusOneOverP_ratio.png");
  c1->cd();
  
  float SMSExcChi=130.4,SMSExcBino=1.4,SMSnotExcChi=400.4,SMSnotExcBino=150.4;
  float BR_aaW=0.00229;
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/Electrohiggs
  float bf_HH_bbgg=2*0.561*0.00229, bf_HH_WWgg=2*0.231*0.00229, bf_HH_ZZgg=2*0.0289*0.00229, bf_HH_ttgg=2*0.0615*0.00229;
  //float xsec_SMS_Excluded=3.248,xsec_SMS_notExcluded=0.0294;//change if change SMS points used above - old https://twiki.cern.ch/twiki/bin/viewauth/CMS/Electrohiggs
  float xsec_SMS_Excluded=3.15958,xsec_SMS_notExcluded=0.071515;//change if change SMS points used above - updated https://twiki.cern.ch/twiki/bin/view/CMSPublic/Higgsinos
  float xsec_SMS_Wino_Excluded=5.841;//change if change SMS points used above
  //float nEvents_SMS_WHExcluded=120000.,nEvents_SMS_WHnotExcluded=30000.,nEvents_SMS_ZHExcluded=120000.,nEvents_SMS_ZHnotExcluded=30000.;//change if change SMS points used above
  float nEvents_SMS_WHExcluded=53439.,nEvents_SMS_WHnotExcluded=30000.,nEvents_SMS_ZHExcluded=53327.,nEvents_SMS_ZHnotExcluded=53327.;//change if change SMS points used above
  //float nEvents_SMS_HHExcluded=60000.;//,nEvents_SMS_HH_2W2gExcluded=30000.,nEvents_SMS_HH_2Z2gExcluded=30000.;
  float nEvents_SMS_HHExcluded=60000.,nEvents_SMS_HH_2W2gExcluded=26929.,nEvents_SMS_HH_2Z2gExcluded=26781.,nEvents_SMS_HH_2b2gExcluded=26833.,nEvents_SMS_HH_2tau2gExcluded=26402.;

  TH3F* h_SMS_WH_gg = (TH3F*)f_SMS_WH->Get("gg_SMS_Loose_mChi_mBino_met");h_SMS_WH_gg->Sumw2();
  int mChiBin = h_SMS_WH_gg->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_WH_gg->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_WH_gg_Excluded = (TH1F*)h_SMS_WH_gg->ProjectionZ("SMS_WH_gg_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_WH_gg->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_WH_gg->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_WH_gg_notExcluded = (TH1F*)h_SMS_WH_gg->ProjectionZ("SMS_WH_gg_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_ZH_gg = (TH3F*)f_SMS_ZH->Get("gg_SMS_Loose_mChi_mBino_met");h_SMS_ZH_gg->Sumw2();
  mChiBin = h_SMS_ZH_gg->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_ZH_gg->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_ZH_gg_Excluded = (TH1F*)h_SMS_ZH_gg->ProjectionZ("SMS_ZH_gg_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_ZH_gg->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_ZH_gg->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_ZH_gg_notExcluded = (TH1F*)h_SMS_ZH_gg->ProjectionZ("SMS_ZH_gg_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2b2g_gg = (TH3F*)f_SMS_HH_2b2g->Get("gg_SMS_Loose_mChi_mBino_met");h_SMS_HH_2b2g_gg->Sumw2();
  mChiBin = h_SMS_HH_2b2g_gg->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2b2g_gg->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2b2g_gg_Excluded = (TH1F*)h_SMS_HH_2b2g_gg->ProjectionZ("SMS_HH_2b2g_gg_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2b2g_gg->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2b2g_gg->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_HH_2b2g_gg_notExcluded = (TH1F*)h_SMS_HH_2b2g_gg->ProjectionZ("SMS_HH_2b2g_gg_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2W2g_gg = (TH3F*)f_SMS_HH_2W2g->Get("gg_SMS_Loose_mChi_mBino_met");h_SMS_HH_2W2g_gg->Sumw2();
  mChiBin = h_SMS_HH_2W2g_gg->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2W2g_gg->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2W2g_gg_Excluded = (TH1F*)h_SMS_HH_2W2g_gg->ProjectionZ("SMS_HH_2W2g_gg_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2W2g_gg->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2W2g_gg->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_HH_2W2g_gg_notExcluded = (TH1F*)h_SMS_HH_2W2g_gg->ProjectionZ("SMS_HH_2W2g_gg_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2Z2g_gg = (TH3F*)f_SMS_HH_2Z2g->Get("gg_SMS_Loose_mChi_mBino_met");h_SMS_HH_2Z2g_gg->Sumw2();
  mChiBin = h_SMS_HH_2Z2g_gg->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2Z2g_gg->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2Z2g_gg_Excluded = (TH1F*)h_SMS_HH_2Z2g_gg->ProjectionZ("SMS_HH_2Z2g_gg_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2Z2g_gg->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2Z2g_gg->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_HH_2Z2g_gg_notExcluded = (TH1F*)h_SMS_HH_2Z2g_gg->ProjectionZ("SMS_HH_2Z2g_gg_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");


  TH3F* h_SMS_WH_gg_1Ele = (TH3F*)f_SMS_WH->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_met");h_SMS_WH_gg_1Ele->Sumw2();
  TH3F* h_SMS_WH_gg_1Ele_InvMass = (TH3F*)f_SMS_WH_InvMass->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_InvMass");h_SMS_WH_gg_1Ele_InvMass->Sumw2();
  TH3F* h_SMS_HH_2W2g_1Ele_InvMass = (TH3F*)f_SMS_HH2W_InvMass->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_InvMass");h_SMS_HH_2W2g_1Ele_InvMass->Sumw2();
  TH3F* h_SMS_HH_2Z2g_1Ele_InvMass = (TH3F*)f_SMS_HH2Z_InvMass->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_InvMass");h_SMS_HH_2Z2g_1Ele_InvMass->Sumw2();
  TH3F* h_SMS_HH_2tau2g_1Ele_InvMass = (TH3F*)f_SMS_HH2tau_InvMass->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_InvMass");h_SMS_HH_2tau2g_1Ele_InvMass->Sumw2();
  mChiBin = h_SMS_WH_gg_1Ele->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_WH_gg_1Ele->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_WH_gg_1Ele_Excluded = (TH1F*)h_SMS_WH_gg_1Ele->ProjectionZ("SMS_WH_gg_1Ele_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  TH1F* h_SMS_WH_gg_1Ele_Excluded_InvMass = (TH1F*)h_SMS_WH_gg_1Ele_InvMass->ProjectionZ("SMS_WH_gg_1Ele_Excluded_InvMass",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  TH1F* h_SMS_HH_2W2g_1Ele_Excluded_InvMass = (TH1F*)h_SMS_HH_2W2g_1Ele_InvMass->ProjectionZ("SMS_HH_2W2g_1Ele_Excluded_InvMass",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  TH1F* h_SMS_HH_2Z2g_1Ele_Excluded_InvMass = (TH1F*)h_SMS_HH_2Z2g_1Ele_InvMass->ProjectionZ("SMS_HH_2Z2g_1Ele_Excluded_InvMass",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  TH1F* h_SMS_HH_2tau2g_1Ele_Excluded_InvMass = (TH1F*)h_SMS_HH_2tau2g_1Ele_InvMass->ProjectionZ("SMS_HH_2tau2g_1Ele_Excluded_InvMass",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_WH_gg_1Ele->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_WH_gg_1Ele->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_WH_gg_1Ele_notExcluded = (TH1F*)h_SMS_WH_gg_1Ele->ProjectionZ("SMS_WH_gg_1Ele_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_ZH_gg_1Ele = (TH3F*)f_SMS_ZH->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_met");h_SMS_ZH_gg_1Ele->Sumw2();
  mChiBin = h_SMS_ZH_gg_1Ele->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_ZH_gg_1Ele->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_ZH_gg_1Ele_Excluded = (TH1F*)h_SMS_ZH_gg_1Ele->ProjectionZ("SMS_ZH_gg_1Ele_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_ZH_gg_1Ele->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_ZH_gg_1Ele->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_ZH_gg_1Ele_notExcluded = (TH1F*)h_SMS_ZH_gg_1Ele->ProjectionZ("SMS_ZH_gg_1Ele_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2b2g_gg_1Ele = (TH3F*)f_SMS_HH_2b2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_met");h_SMS_HH_2b2g_gg_1Ele->Sumw2();
  mChiBin = h_SMS_HH_2b2g_gg_1Ele->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2b2g_gg_1Ele->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2b2g_gg_1Ele_Excluded = (TH1F*)h_SMS_HH_2b2g_gg_1Ele->ProjectionZ("SMS_HH_2b2g_gg_1Ele_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2b2g_gg_1Ele->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2b2g_gg_1Ele->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_HH_2b2g_gg_1Ele_notExcluded = (TH1F*)h_SMS_HH_2b2g_gg_1Ele->ProjectionZ("SMS_HH_2b2g_gg_1Ele_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2W2g_gg_1Ele = (TH3F*)f_SMS_HH_2W2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_met");h_SMS_HH_2W2g_gg_1Ele->Sumw2();
  mChiBin = h_SMS_HH_2W2g_gg_1Ele->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2W2g_gg_1Ele->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2W2g_gg_1Ele_Excluded = (TH1F*)h_SMS_HH_2W2g_gg_1Ele->ProjectionZ("SMS_HH_2W2g_gg_1Ele_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2W2g_gg_1Ele->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2W2g_gg_1Ele->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_HH_2W2g_gg_1Ele_notExcluded = (TH1F*)h_SMS_HH_2W2g_gg_1Ele->ProjectionZ("SMS_HH_2W2g_gg_1Ele_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2Z2g_gg_1Ele = (TH3F*)f_SMS_HH_2Z2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_met");h_SMS_HH_2Z2g_gg_1Ele->Sumw2();
  mChiBin = h_SMS_HH_2Z2g_gg_1Ele->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2Z2g_gg_1Ele->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2Z2g_gg_1Ele_Excluded = (TH1F*)h_SMS_HH_2Z2g_gg_1Ele->ProjectionZ("SMS_HH_2Z2g_gg_1Ele_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2Z2g_gg_1Ele->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2Z2g_gg_1Ele->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_HH_2Z2g_gg_1Ele_notExcluded = (TH1F*)h_SMS_HH_2Z2g_gg_1Ele->ProjectionZ("SMS_HH_2Z2g_gg_1Ele_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_WH_gg_1Ele_MT = (TH3F*)f_SMS_WH->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT");h_SMS_WH_gg_1Ele_MT->Sumw2();
  mChiBin = h_SMS_WH_gg_1Ele_MT->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_WH_gg_1Ele_MT->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_WH_gg_1Ele_MT_Excluded = (TH1F*)h_SMS_WH_gg_1Ele_MT->ProjectionZ("SMS_WH_gg_1Ele_MT_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_WH_gg_1Ele_MT->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_WH_gg_1Ele_MT->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_WH_gg_1Ele_MT_notExcluded = (TH1F*)h_SMS_WH_gg_1Ele_MT->ProjectionZ("SMS_WH_gg_1Ele_MT_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_ZH_gg_1Ele_MT = (TH3F*)f_SMS_ZH->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT");h_SMS_ZH_gg_1Ele_MT->Sumw2();
  mChiBin = h_SMS_ZH_gg_1Ele_MT->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_ZH_gg_1Ele_MT->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_ZH_gg_1Ele_MT_Excluded = (TH1F*)h_SMS_ZH_gg_1Ele_MT->ProjectionZ("SMS_ZH_gg_1Ele_MT_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_ZH_gg_1Ele_MT->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_ZH_gg_1Ele_MT->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_ZH_gg_1Ele_MT_notExcluded = (TH1F*)h_SMS_ZH_gg_1Ele_MT->ProjectionZ("SMS_ZH_gg_1Ele_MT_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2b2g_gg_1Ele_MT = (TH3F*)f_SMS_HH_2b2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT");h_SMS_HH_2b2g_gg_1Ele_MT->Sumw2();
  mChiBin = h_SMS_HH_2b2g_gg_1Ele_MT->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2b2g_gg_1Ele_MT->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2b2g_gg_1Ele_MT_Excluded = (TH1F*)h_SMS_HH_2b2g_gg_1Ele_MT->ProjectionZ("SMS_HH_2b2g_gg_1Ele_MT_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2b2g_gg_1Ele_MT->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2b2g_gg_1Ele_MT->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_HH_2b2g_gg_1Ele_MT_notExcluded = (TH1F*)h_SMS_HH_2b2g_gg_1Ele_MT->ProjectionZ("SMS_HH_2b2g_gg_1Ele_MT_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2W2g_gg_1Ele_MT = (TH3F*)f_SMS_HH_2W2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT");h_SMS_HH_2W2g_gg_1Ele_MT->Sumw2();
  mChiBin = h_SMS_HH_2W2g_gg_1Ele_MT->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2W2g_gg_1Ele_MT->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2W2g_gg_1Ele_MT_Excluded = (TH1F*)h_SMS_HH_2W2g_gg_1Ele_MT->ProjectionZ("SMS_HH_2W2g_gg_1Ele_MT_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2W2g_gg_1Ele_MT->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2W2g_gg_1Ele_MT->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_HH_2W2g_gg_1Ele_MT_notExcluded = (TH1F*)h_SMS_HH_2W2g_gg_1Ele_MT->ProjectionZ("SMS_HH_2W2g_gg_1Ele_MT_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2Z2g_gg_1Ele_MT = (TH3F*)f_SMS_HH_2Z2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT");h_SMS_HH_2Z2g_gg_1Ele_MT->Sumw2();
  mChiBin = h_SMS_HH_2Z2g_gg_1Ele_MT->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2Z2g_gg_1Ele_MT->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2Z2g_gg_1Ele_MT_Excluded = (TH1F*)h_SMS_HH_2Z2g_gg_1Ele_MT->ProjectionZ("SMS_HH_2Z2g_gg_1Ele_MT_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2Z2g_gg_1Ele_MT->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2Z2g_gg_1Ele_MT->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_HH_2Z2g_gg_1Ele_MT_notExcluded = (TH1F*)h_SMS_HH_2Z2g_gg_1Ele_MT->ProjectionZ("SMS_HH_2Z2g_gg_1Ele_MT_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2tau2g_gg_1Ele_MT = (TH3F*)f_SMS_HH_2tau2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT");h_SMS_HH_2tau2g_gg_1Ele_MT->Sumw2();
  mChiBin = h_SMS_HH_2tau2g_gg_1Ele_MT->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2tau2g_gg_1Ele_MT->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2tau2g_gg_1Ele_MT_Excluded = (TH1F*)h_SMS_HH_2tau2g_gg_1Ele_MT->ProjectionZ("SMS_HH_2tau2g_gg_1Ele_MT_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2tau2g_gg_1Mu_MT = (TH3F*)f_SMS_HH_2tau2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT");h_SMS_HH_2tau2g_gg_1Mu_MT->Sumw2();
  mChiBin = h_SMS_HH_2tau2g_gg_1Mu_MT->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2tau2g_gg_1Mu_MT->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2tau2g_gg_1Mu_MT_Excluded = (TH1F*)h_SMS_HH_2tau2g_gg_1Mu_MT->ProjectionZ("SMS_HH_2tau2g_gg_1Mu_MT_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2tau2g_gg_1Ele = (TH3F*)f_SMS_HH_2tau2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_met");h_SMS_HH_2tau2g_gg_1Ele->Sumw2();
  mChiBin = h_SMS_HH_2tau2g_gg_1Ele->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2tau2g_gg_1Ele->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2tau2g_gg_1Ele_Excluded = (TH1F*)h_SMS_HH_2tau2g_gg_1Ele->ProjectionZ("SMS_HH_2tau2g_gg_1Ele_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2tau2g_gg_1Mu = (TH3F*)f_SMS_HH_2tau2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_met");h_SMS_HH_2tau2g_gg_1Mu->Sumw2();
  mChiBin = h_SMS_HH_2tau2g_gg_1Mu->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2tau2g_gg_1Mu->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2tau2g_gg_1Mu_Excluded = (TH1F*)h_SMS_HH_2tau2g_gg_1Mu->ProjectionZ("SMS_HH_2tau2g_gg_1Mu_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  




  TH3F* h_SMS_WH_gg_1Ele_MT_noSF = (TH3F*)f_SMS_WH->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor");h_SMS_WH_gg_1Ele_MT_noSF->Sumw2();
  float norm_SMS_WH_ele = h_SMS_WH_gg_1Ele_MT_noSF->GetEntries()/h_SMS_WH_gg_1Ele_MT_noSF->Integral(0,-1);
  TH3F* h_SMS_ZH_gg_1Ele_MT_noSF = (TH3F*)f_SMS_ZH->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor");h_SMS_ZH_gg_1Ele_MT_noSF->Sumw2();
  float norm_SMS_ZH_ele = h_SMS_ZH_gg_1Ele_MT_noSF->GetEntries()/h_SMS_ZH_gg_1Ele_MT_noSF->Integral(0,-1);
  TH3F* h_SMS_HH_2b2g_gg_1Ele_MT_noSF = (TH3F*)f_SMS_HH_2b2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor");h_SMS_HH_2b2g_gg_1Ele_MT_noSF->Sumw2();
  float norm_SMS_HH_2b2g_ele = h_SMS_HH_2b2g_gg_1Ele_MT_noSF->GetEntries()/h_SMS_HH_2b2g_gg_1Ele_MT_noSF->Integral(0,-1);
  TH3F* h_SMS_HH_2W2g_gg_1Ele_MT_noSF = (TH3F*)f_SMS_HH_2W2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor");h_SMS_HH_2W2g_gg_1Ele_MT_noSF->Sumw2();
  float norm_SMS_HH_2W2g_ele = h_SMS_HH_2W2g_gg_1Ele_MT_noSF->GetEntries()/h_SMS_HH_2W2g_gg_1Ele_MT_noSF->Integral(0,-1);
  TH3F* h_SMS_HH_2Z2g_gg_1Ele_MT_noSF = (TH3F*)f_SMS_HH_2Z2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor");h_SMS_HH_2Z2g_gg_1Ele_MT_noSF->Sumw2();
  float norm_SMS_HH_2Z2g_ele = h_SMS_HH_2Z2g_gg_1Ele_MT_noSF->GetEntries()/h_SMS_HH_2Z2g_gg_1Ele_MT_noSF->Integral(0,-1);
  TH3F* h_SMS_HH_2tau2g_gg_1Ele_MT_noSF = (TH3F*)f_SMS_HH_2tau2g->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor");h_SMS_HH_2tau2g_gg_1Ele_MT_noSF->Sumw2();
  float norm_SMS_HH_2tau2g_ele = h_SMS_HH_2tau2g_gg_1Ele_MT_noSF->GetEntries()/h_SMS_HH_2tau2g_gg_1Ele_MT_noSF->Integral(0,-1);



  float PhoEffScale2 = 0.994*0.994;
  float scale_SMS_WHExcluded = PhoEffScale2*L_int*xsec_SMS_Wino_Excluded*BR_aaW/nEvents_SMS_WHExcluded;
  float scale_SMS_ZHExcluded = PhoEffScale2*L_int*xsec_SMS_Excluded*BR_aaW/nEvents_SMS_ZHExcluded;
  float scale_SMS_WHnotExcluded = PhoEffScale2*L_int*xsec_SMS_notExcluded*BR_aaW/nEvents_SMS_WHnotExcluded;
  float scale_SMS_ZHnotExcluded = PhoEffScale2*L_int*xsec_SMS_notExcluded*BR_aaW/nEvents_SMS_ZHnotExcluded;
  float scale_SMS_HHExcluded = PhoEffScale2*L_int*xsec_SMS_Excluded*BR_aaW/nEvents_SMS_HHExcluded;
  float scale_SMS_HH_2W2gExcluded = PhoEffScale2*L_int*xsec_SMS_Excluded*bf_HH_WWgg/nEvents_SMS_HH_2W2gExcluded;
  float scale_SMS_HH_2Z2gExcluded = PhoEffScale2*L_int*xsec_SMS_Excluded*bf_HH_ZZgg/nEvents_SMS_HH_2Z2gExcluded;
  float scale_SMS_HH_2tau2gExcluded = PhoEffScale2*L_int*xsec_SMS_Excluded*bf_HH_ttgg/nEvents_SMS_HH_2tau2gExcluded;
  float scale_SMS_HH_2b2gExcluded = PhoEffScale2*L_int*xsec_SMS_Excluded*bf_HH_bbgg/nEvents_SMS_HH_2b2gExcluded;
  
  //now inclusive
  h_SMS_WH_gg_Excluded->Scale(scale_SMS_WHExcluded);
  h_SMS_WH_gg_notExcluded->Scale(scale_SMS_WHnotExcluded);
  h_SMS_ZH_gg_Excluded->Scale(scale_SMS_ZHExcluded);
  h_SMS_ZH_gg_notExcluded->Scale(scale_SMS_ZHnotExcluded);
  
  h_SMS_HH_2b2g_gg_Excluded->Scale(scale_SMS_HH_2b2gExcluded);
  h_SMS_HH_2W2g_gg_Excluded->Scale(scale_SMS_HH_2W2gExcluded);
  h_SMS_HH_2Z2g_gg_Excluded->Scale(scale_SMS_HH_2Z2gExcluded);


  h_SMS_WH_gg_Excluded->SetLineColor(kBlue);h_SMS_WH_gg_Excluded->SetLineWidth(3);
  h_SMS_ZH_gg_Excluded->SetLineColor(kRed-4);h_SMS_ZH_gg_Excluded->SetLineWidth(3);
  h_SMS_WH_gg_notExcluded->SetLineColor(kRed);h_SMS_WH_gg_notExcluded->SetLineWidth(3);
  h_SMS_ZH_gg_notExcluded->SetLineColor(kRed);h_SMS_ZH_gg_notExcluded->SetLineWidth(3);

  h_SMS_HH_2b2g_gg_Excluded->SetLineColor(kGreen);h_SMS_HH_2b2g_gg_Excluded->SetLineWidth(3);
  h_SMS_HH_2W2g_gg_Excluded->SetLineColor(kViolet);h_SMS_HH_2W2g_gg_Excluded->SetLineWidth(3);
  h_SMS_HH_2Z2g_gg_Excluded->SetLineColor(kOrange);h_SMS_HH_2Z2g_gg_Excluded->SetLineWidth(3);

  TH1F* h_SMS_WH_gg_Excluded_Rebin = (TH1F*)h_SMS_WH_gg_Excluded->Rebin(NmetBins,"h_SMS_WH_gg_Excluded_Rebin ",xbins);
  TH1F* h_SMS_WH_gg_notExcluded_Rebin = (TH1F*)h_SMS_WH_gg_notExcluded->Rebin(NmetBins,"h_SMS_WH_gg_notExcluded_Rebin ",xbins);
  TH1F* h_SMS_ZH_gg_Excluded_Rebin = (TH1F*)h_SMS_ZH_gg_Excluded->Rebin(NmetBins,"h_SMS_ZH_gg_Excluded_Rebin ",xbins);
  TH1F* h_SMS_ZH_gg_notExcluded_Rebin = (TH1F*)h_SMS_ZH_gg_notExcluded->Rebin(NmetBins,"h_SMS_ZH_gg_notExcluded_Rebin ",xbins);

  TH1F* h_SMS_HH_2b2g_gg_Excluded_Rebin = (TH1F*)h_SMS_HH_2b2g_gg_Excluded->Rebin(NmetBins,"h_SMS_HH_2b2g_gg_Excluded_Rebin ",xbins);
  TH1F* h_SMS_HH_2W2g_gg_Excluded_Rebin = (TH1F*)h_SMS_HH_2W2g_gg_Excluded->Rebin(NmetBins,"h_SMS_HH_2W2g_gg_Excluded_Rebin ",xbins);
  TH1F* h_SMS_HH_2Z2g_gg_Excluded_Rebin = (TH1F*)h_SMS_HH_2Z2g_gg_Excluded->Rebin(NmetBins,"h_SMS_HH_2Z2g_gg_Excluded_Rebin ",xbins);


  AddOverflowToLastBin(h_SMS_WH_gg_Excluded_Rebin);AddOverflowToLastBin(h_SMS_WH_gg_notExcluded_Rebin);
  AddOverflowToLastBin(h_SMS_ZH_gg_Excluded_Rebin);AddOverflowToLastBin(h_SMS_ZH_gg_notExcluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2b2g_gg_Excluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2W2g_gg_Excluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2Z2g_gg_Excluded_Rebin);
  DivideByBinWidth(h_SMS_WH_gg_Excluded_Rebin);DivideByBinWidth(h_SMS_WH_gg_notExcluded_Rebin);
  DivideByBinWidth(h_SMS_ZH_gg_Excluded_Rebin);DivideByBinWidth(h_SMS_ZH_gg_notExcluded_Rebin);
  DivideByBinWidth(h_SMS_HH_2b2g_gg_Excluded_Rebin);
  DivideByBinWidth(h_SMS_HH_2W2g_gg_Excluded_Rebin);
  DivideByBinWidth(h_SMS_HH_2Z2g_gg_Excluded_Rebin);
  //now ele
  h_SMS_WH_gg_1Ele_Excluded->Scale(norm_SMS_WH_ele*scale_SMS_WHExcluded);
  h_SMS_WH_gg_1Ele_Excluded_InvMass->Scale(norm_SMS_WH_ele*scale_SMS_WHExcluded);
  h_SMS_HH_2W2g_1Ele_Excluded_InvMass->Scale(norm_SMS_HH_2W2g_ele*scale_SMS_HH_2W2gExcluded);
  h_SMS_HH_2Z2g_1Ele_Excluded_InvMass->Scale(norm_SMS_HH_2Z2g_ele*scale_SMS_HH_2Z2gExcluded);
  h_SMS_HH_2tau2g_1Ele_Excluded_InvMass->Scale(norm_SMS_HH_2tau2g_ele*scale_SMS_HH_2tau2gExcluded);
  h_SMS_WH_gg_1Ele_notExcluded->Scale(norm_SMS_WH_ele*scale_SMS_WHnotExcluded);
  h_SMS_ZH_gg_1Ele_Excluded->Scale(norm_SMS_ZH_ele*scale_SMS_ZHExcluded);
  h_SMS_ZH_gg_1Ele_notExcluded->Scale(norm_SMS_ZH_ele*scale_SMS_ZHnotExcluded);
  h_SMS_HH_2b2g_gg_1Ele_Excluded->Scale(norm_SMS_HH_2b2g_ele*scale_SMS_HH_2b2gExcluded);
  h_SMS_HH_2W2g_gg_1Ele_Excluded->Scale(norm_SMS_HH_2W2g_ele*scale_SMS_HH_2W2gExcluded);
  h_SMS_HH_2Z2g_gg_1Ele_Excluded->Scale(norm_SMS_HH_2Z2g_ele*scale_SMS_HH_2Z2gExcluded);
  h_SMS_HH_2tau2g_gg_1Ele_Excluded->Scale(norm_SMS_HH_2tau2g_ele*scale_SMS_HH_2tau2gExcluded);

  h_SMS_WH_gg_1Ele_Excluded->SetLineColor(kBlue);h_SMS_WH_gg_1Ele_Excluded->SetLineWidth(3);
  h_SMS_WH_gg_1Ele_Excluded_InvMass->SetLineColor(kRed);h_SMS_WH_gg_1Ele_Excluded_InvMass->SetLineWidth(0);
  h_SMS_WH_gg_1Ele_Excluded_InvMass->SetFillColor(kRed);h_SMS_WH_gg_1Ele_Excluded_InvMass->SetFillStyle(3004);
  h_SMS_HH_2W2g_1Ele_Excluded_InvMass->SetLineColor(kRed);h_SMS_HH_2W2g_1Ele_Excluded_InvMass->SetLineWidth(0);
  h_SMS_HH_2W2g_1Ele_Excluded_InvMass->SetFillColor(kRed);h_SMS_HH_2W2g_1Ele_Excluded_InvMass->SetFillStyle(3004);
  h_SMS_ZH_gg_1Ele_Excluded->SetLineColor(kRed-4);h_SMS_ZH_gg_1Ele_Excluded->SetLineWidth(3);
  h_SMS_WH_gg_1Ele_notExcluded->SetLineColor(kRed);h_SMS_WH_gg_1Ele_notExcluded->SetLineWidth(3);
  h_SMS_ZH_gg_1Ele_notExcluded->SetLineColor(kRed);h_SMS_ZH_gg_1Ele_notExcluded->SetLineWidth(3);
  h_SMS_HH_2b2g_gg_1Ele_Excluded->SetLineColor(kGreen);h_SMS_HH_2b2g_gg_1Ele_Excluded->SetLineWidth(3);
  h_SMS_HH_2W2g_gg_1Ele_Excluded->SetLineColor(kViolet);h_SMS_HH_2W2g_gg_1Ele_Excluded->SetLineWidth(3);
  h_SMS_HH_2Z2g_gg_1Ele_Excluded->SetLineColor(kOrange);h_SMS_HH_2Z2g_gg_1Ele_Excluded->SetLineWidth(3);


  h_SMS_WH_gg_1Ele_MT_Excluded->Scale(norm_SMS_WH_ele*scale_SMS_WHExcluded);
  h_SMS_WH_gg_1Ele_MT_notExcluded->Scale(norm_SMS_WH_ele*scale_SMS_WHnotExcluded);
  h_SMS_ZH_gg_1Ele_MT_Excluded->Scale(norm_SMS_ZH_ele*scale_SMS_ZHExcluded);
  h_SMS_ZH_gg_1Ele_MT_notExcluded->Scale(norm_SMS_ZH_ele*scale_SMS_ZHnotExcluded);
  h_SMS_HH_2b2g_gg_1Ele_MT_Excluded->Scale(norm_SMS_HH_2b2g_ele*scale_SMS_HH_2b2gExcluded);
  h_SMS_HH_2W2g_gg_1Ele_MT_Excluded->Scale(norm_SMS_HH_2W2g_ele*scale_SMS_HH_2W2gExcluded);
  h_SMS_HH_2Z2g_gg_1Ele_MT_Excluded->Scale(norm_SMS_HH_2Z2g_ele*scale_SMS_HH_2Z2gExcluded);
  h_SMS_HH_2tau2g_gg_1Ele_MT_Excluded->Scale(norm_SMS_HH_2tau2g_ele*scale_SMS_HH_2tau2gExcluded);

  h_SMS_WH_gg_1Ele_MT_Excluded->SetLineColor(kBlue);h_SMS_WH_gg_1Ele_MT_Excluded->SetLineWidth(3);
  h_SMS_ZH_gg_1Ele_MT_Excluded->SetLineColor(kRed-4);h_SMS_ZH_gg_1Ele_MT_Excluded->SetLineWidth(3);
  h_SMS_WH_gg_1Ele_MT_notExcluded->SetLineColor(kRed);h_SMS_WH_gg_1Ele_MT_notExcluded->SetLineWidth(3);
  h_SMS_ZH_gg_1Ele_MT_notExcluded->SetLineColor(kRed);h_SMS_ZH_gg_1Ele_MT_notExcluded->SetLineWidth(3);
  h_SMS_HH_2b2g_gg_1Ele_MT_Excluded->SetLineColor(kGreen);h_SMS_HH_2b2g_gg_1Ele_MT_Excluded->SetLineWidth(3);
  h_SMS_HH_2W2g_gg_1Ele_MT_Excluded->SetLineColor(kViolet);h_SMS_HH_2W2g_gg_1Ele_MT_Excluded->SetLineWidth(3);
  h_SMS_HH_2Z2g_gg_1Ele_MT_Excluded->SetLineColor(kOrange);h_SMS_HH_2Z2g_gg_1Ele_MT_Excluded->SetLineWidth(3);



  for(int i=0;i<h_SMS_WH_gg_1Ele_Excluded_InvMass->GetNbinsX();i++){
    if(h_SMS_WH_gg_1Ele_Excluded_InvMass->GetBinLowEdge(i)<120 || h_SMS_WH_gg_1Ele_Excluded_InvMass->GetBinLowEdge(i)>=131)h_SMS_WH_gg_1Ele_Excluded_InvMass->SetBinContent(i,0);
    if(h_SMS_HH_2W2g_1Ele_Excluded_InvMass->GetBinLowEdge(i)<120 || h_SMS_HH_2W2g_1Ele_Excluded_InvMass->GetBinLowEdge(i)>=131)h_SMS_HH_2W2g_1Ele_Excluded_InvMass->SetBinContent(i,0);
    if(h_SMS_HH_2Z2g_1Ele_Excluded_InvMass->GetBinLowEdge(i)<120 || h_SMS_HH_2Z2g_1Ele_Excluded_InvMass->GetBinLowEdge(i)>=131)h_SMS_HH_2Z2g_1Ele_Excluded_InvMass->SetBinContent(i,0);
    if(h_SMS_HH_2tau2g_1Ele_Excluded_InvMass->GetBinLowEdge(i)<120 || h_SMS_HH_2tau2g_1Ele_Excluded_InvMass->GetBinLowEdge(i)>=131)h_SMS_HH_2tau2g_1Ele_Excluded_InvMass->SetBinContent(i,0);
  }
  h_SMS_HH_2W2g_1Ele_Excluded_InvMass->Add(h_SMS_HH_2Z2g_1Ele_Excluded_InvMass);
  h_SMS_HH_2W2g_1Ele_Excluded_InvMass->Add(h_SMS_HH_2tau2g_1Ele_Excluded_InvMass);
  /*
  float smsInvMassScale = 4.78/h_SMS_WH_gg_1Ele_Excluded_InvMass->Integral();
  h_SMS_WH_gg_1Ele_Excluded_InvMass->Scale(smsInvMassScale);
  cout<<"smsInvMassScale: "<<smsInvMassScale<<endl;
  */
  //h_SMS_WH_gg_1Ele_Excluded->Add(h_SMS_ZH_gg_1Ele_Excluded);
  //h_SMS_WH_gg_1Ele_notExcluded->Add(h_SMS_ZH_gg_1Ele_notExcluded);

  TH1F* h_SMS_WH_gg_1Ele_Excluded_Rebin = (TH1F*)h_SMS_WH_gg_1Ele_Excluded->Rebin(NmetBinsXtraWide,"h_SMS_WH_gg_1Ele_Excluded_Rebin ",xbinsXtraWide);
  TH1F* h_SMS_WH_gg_1Ele_notExcluded_Rebin = (TH1F*)h_SMS_WH_gg_1Ele_notExcluded->Rebin(NmetBinsXtraWide,"h_SMS_WH_gg_1Ele_notExcluded_Rebin ",xbinsXtraWide);
  TH1F* h_SMS_ZH_gg_1Ele_Excluded_Rebin = (TH1F*)h_SMS_ZH_gg_1Ele_Excluded->Rebin(NmetBinsXtraWide,"h_SMS_ZH_gg_1Ele_Excluded_Rebin ",xbinsXtraWide);
  TH1F* h_SMS_ZH_gg_1Ele_notExcluded_Rebin = (TH1F*)h_SMS_ZH_gg_1Ele_notExcluded->Rebin(NmetBinsXtraWide,"h_SMS_ZH_gg_1Ele_notExcluded_Rebin ",xbinsXtraWide);
  TH1F* h_SMS_HH_2b2g_gg_1Ele_Excluded_Rebin = (TH1F*)h_SMS_HH_2b2g_gg_1Ele_Excluded->Rebin(NmetBinsXtraWide,"h_SMS_HH_2b2g_gg_1Ele_Excluded_Rebin ",xbinsXtraWide);
  TH1F* h_SMS_HH_2W2g_gg_1Ele_Excluded_Rebin = (TH1F*)h_SMS_HH_2W2g_gg_1Ele_Excluded->Rebin(NmetBinsXtraWide,"h_SMS_HH_2W2g_gg_1Ele_Excluded_Rebin ",xbinsXtraWide);
  TH1F* h_SMS_HH_2Z2g_gg_1Ele_Excluded_Rebin = (TH1F*)h_SMS_HH_2Z2g_gg_1Ele_Excluded->Rebin(NmetBinsXtraWide,"h_SMS_HH_2Z2g_gg_1Ele_Excluded_Rebin ",xbinsXtraWide);
  TH1F* h_SMS_HH_2tau2g_gg_1Ele_Excluded_Rebin = (TH1F*)h_SMS_HH_2tau2g_gg_1Ele_Excluded->Rebin(NmetBinsXtraWide,"h_SMS_HH_2tau2g_gg_1Ele_Excluded_Rebin ",xbinsXtraWide);

  TH1F* h_SMS_WH_gg_1Ele_MT_Excluded_Rebin = (TH1F*)h_SMS_WH_gg_1Ele_MT_Excluded->Rebin(nMTbins,"h_SMS_WH_gg_1Ele_MT_Excluded_Rebin ",MTbins);
  TH1F* h_SMS_WH_gg_1Ele_MT_notExcluded_Rebin = (TH1F*)h_SMS_WH_gg_1Ele_MT_notExcluded->Rebin(nMTbins,"h_SMS_WH_gg_1Ele_MT_notExcluded_Rebin ",MTbins);
  TH1F* h_SMS_ZH_gg_1Ele_MT_Excluded_Rebin = (TH1F*)h_SMS_ZH_gg_1Ele_MT_Excluded->Rebin(nMTbins,"h_SMS_ZH_gg_1Ele_MT_Excluded_Rebin ",MTbins);
  TH1F* h_SMS_ZH_gg_1Ele_MT_notExcluded_Rebin = (TH1F*)h_SMS_ZH_gg_1Ele_MT_notExcluded->Rebin(nMTbins,"h_SMS_ZH_gg_1Ele_MT_notExcluded_Rebin ",MTbins);
  TH1F* h_SMS_HH_2b2g_gg_1Ele_MT_Excluded_Rebin = (TH1F*)h_SMS_HH_2b2g_gg_1Ele_MT_Excluded->Rebin(nMTbins,"h_SMS_HH_2b2g_gg_1Ele_MT_Excluded_Rebin ",MTbins);
  TH1F* h_SMS_HH_2W2g_gg_1Ele_MT_Excluded_Rebin = (TH1F*)h_SMS_HH_2W2g_gg_1Ele_MT_Excluded->Rebin(nMTbins,"h_SMS_HH_2W2g_gg_1Ele_MT_Excluded_Rebin ",MTbins);
  TH1F* h_SMS_HH_2Z2g_gg_1Ele_MT_Excluded_Rebin = (TH1F*)h_SMS_HH_2Z2g_gg_1Ele_MT_Excluded->Rebin(nMTbins,"h_SMS_HH_2Z2g_gg_1Ele_MT_Excluded_Rebin ",MTbins);
  TH1F* h_SMS_HH_2tau2g_gg_1Ele_MT_Excluded_Rebin = (TH1F*)h_SMS_HH_2tau2g_gg_1Ele_MT_Excluded->Rebin(nMTbins,"h_SMS_HH_2tau2g_gg_1Ele_MT_Excluded_Rebin ",MTbins);

  /*
  int SMSoverflowbin=h_SMS_WH_gg_1Ele_Excluded_Rebin->FindBin(metPlotXmaxXtraWide+1),SMSlastbin=h_SMS_WH_gg_1Ele_Excluded_Rebin->FindBin(metPlotXmaxXtraWide-1);
  double SMSoverflowerr=0,SMSlastbinerr=0;
  double SMSoverflow = h_SMS_WH_gg_1Ele_Excluded_Rebin->IntegralAndError(SMSoverflowbin,-1,SMSoverflowerr);
  double SMSlast = h_SMS_WH_gg_1Ele_Excluded_Rebin->IntegralAndError(SMSlastbin,SMSlastbin,SMSlastbinerr);
  h_SMS_WH_gg_1Ele_Excluded_Rebin->SetBinContent(SMSlastbin,SMSlast+SMSoverflow);
  h_SMS_WH_gg_1Ele_Excluded_Rebin->SetBinError(SMSlastbin,sqrt(SMSlastbinerr*SMSlastbinerr+SMSoverflow*SMSoverflow));
  h_SMS_WH_gg_1Ele_Excluded_Rebin->SetBinContent(SMSoverflowbin,0);h_SMS_WH_gg_1Ele_Excluded_Rebin->SetBinError(SMSoverflowbin,0);
  SMSoverflow = h_SMS_WH_gg_1Ele_notExcluded_Rebin->IntegralAndError(SMSoverflowbin,-1,SMSoverflowerr);
  SMSlast = h_SMS_WH_gg_1Ele_notExcluded_Rebin->IntegralAndError(SMSlastbin,SMSlastbin,SMSlastbinerr);
  h_SMS_WH_gg_1Ele_notExcluded_Rebin->SetBinContent(SMSlastbin,SMSlast+SMSoverflow);
  h_SMS_WH_gg_1Ele_notExcluded_Rebin->SetBinError(SMSlastbin,sqrt(SMSlastbinerr*SMSlastbinerr+SMSoverflow*SMSoverflow));
  h_SMS_WH_gg_1Ele_notExcluded_Rebin->SetBinContent(SMSoverflowbin,0);h_SMS_WH_gg_1Ele_notExcluded_Rebin->SetBinError(SMSoverflowbin,0);
  
  SMSoverflow = h_SMS_ZH_gg_1Ele_Excluded_Rebin->IntegralAndError(SMSoverflowbin,-1,SMSoverflowerr);
  SMSlast = h_SMS_ZH_gg_1Ele_Excluded_Rebin->IntegralAndError(SMSlastbin,SMSlastbin,SMSlastbinerr);
  h_SMS_ZH_gg_1Ele_Excluded_Rebin->SetBinContent(SMSlastbin,SMSlast+SMSoverflow);
  h_SMS_ZH_gg_1Ele_Excluded_Rebin->SetBinError(SMSlastbin,sqrt(SMSlastbinerr*SMSlastbinerr+SMSoverflow*SMSoverflow));
  h_SMS_ZH_gg_1Ele_Excluded_Rebin->SetBinContent(SMSoverflowbin,0);h_SMS_ZH_gg_1Ele_Excluded_Rebin->SetBinError(SMSoverflowbin,0);
  SMSoverflow = h_SMS_ZH_gg_1Ele_notExcluded_Rebin->IntegralAndError(SMSoverflowbin,-1,SMSoverflowerr);
  SMSlast = h_SMS_ZH_gg_1Ele_notExcluded_Rebin->IntegralAndError(SMSlastbin,SMSlastbin,SMSlastbinerr);
  h_SMS_ZH_gg_1Ele_notExcluded_Rebin->SetBinContent(SMSlastbin,SMSlast+SMSoverflow);
  h_SMS_ZH_gg_1Ele_notExcluded_Rebin->SetBinError(SMSlastbin,sqrt(SMSlastbinerr*SMSlastbinerr+SMSoverflow*SMSoverflow));
  h_SMS_ZH_gg_1Ele_notExcluded_Rebin->SetBinContent(SMSoverflowbin,0);h_SMS_ZH_gg_1Ele_notExcluded_Rebin->SetBinError(SMSoverflowbin,0);
  */
  AddOverflowToLastBin(h_SMS_WH_gg_1Ele_Excluded_Rebin);AddOverflowToLastBin(h_SMS_WH_gg_1Ele_notExcluded_Rebin);
  AddOverflowToLastBin(h_SMS_ZH_gg_1Ele_Excluded_Rebin);AddOverflowToLastBin(h_SMS_ZH_gg_1Ele_notExcluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2b2g_gg_1Ele_Excluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2W2g_gg_1Ele_Excluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2Z2g_gg_1Ele_Excluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2tau2g_gg_1Ele_Excluded_Rebin);
  DivideBy15gev(h_SMS_WH_gg_1Ele_Excluded_Rebin);DivideBy15gev(h_SMS_WH_gg_1Ele_notExcluded_Rebin);
  DivideBy15gev(h_SMS_ZH_gg_1Ele_Excluded_Rebin);DivideBy15gev(h_SMS_ZH_gg_1Ele_notExcluded_Rebin);
  DivideBy15gev(h_SMS_HH_2b2g_gg_1Ele_Excluded_Rebin);
  DivideBy15gev(h_SMS_HH_2W2g_gg_1Ele_Excluded_Rebin);
  DivideBy15gev(h_SMS_HH_2Z2g_gg_1Ele_Excluded_Rebin);
  DivideBy15gev(h_SMS_HH_2tau2g_gg_1Ele_Excluded_Rebin);
  
  AddOverflowToLastBin(h_SMS_WH_gg_1Ele_MT_Excluded_Rebin);AddOverflowToLastBin(h_SMS_WH_gg_1Ele_MT_notExcluded_Rebin);
  AddOverflowToLastBin(h_SMS_ZH_gg_1Ele_MT_Excluded_Rebin);AddOverflowToLastBin(h_SMS_ZH_gg_1Ele_MT_notExcluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2b2g_gg_1Ele_MT_Excluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2W2g_gg_1Ele_MT_Excluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2Z2g_gg_1Ele_MT_Excluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2tau2g_gg_1Ele_MT_Excluded_Rebin);
  DivideBy30gev(h_SMS_WH_gg_1Ele_MT_Excluded_Rebin);DivideBy30gev(h_SMS_WH_gg_1Ele_MT_notExcluded_Rebin);
  DivideBy30gev(h_SMS_ZH_gg_1Ele_MT_Excluded_Rebin);DivideBy30gev(h_SMS_ZH_gg_1Ele_MT_notExcluded_Rebin);
  DivideBy30gev(h_SMS_HH_2b2g_gg_1Ele_MT_Excluded_Rebin);
  DivideBy30gev(h_SMS_HH_2W2g_gg_1Ele_MT_Excluded_Rebin);
  DivideBy30gev(h_SMS_HH_2Z2g_gg_1Ele_MT_Excluded_Rebin);
  DivideBy30gev(h_SMS_HH_2tau2g_gg_1Ele_MT_Excluded_Rebin);

  //now muon
  
  TH3F* h_SMS_WH_gg_1Mu = (TH3F*)f_SMS_WH->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_met");h_SMS_WH_gg_1Mu->Sumw2();
  mChiBin = h_SMS_WH_gg_1Mu->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_WH_gg_1Mu->GetYaxis()->FindBin(SMSExcBino);
  TH3F* h_SMS_WH_gg_1Mu_InvMass = (TH3F*)f_SMS_WH_InvMass->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_InvMass");h_SMS_WH_gg_1Mu_InvMass->Sumw2();
  TH3F* h_SMS_HH_2W2g_1Mu_InvMass = (TH3F*)f_SMS_HH2Z_InvMass->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_InvMass");h_SMS_HH_2W2g_1Mu_InvMass->Sumw2();
  TH3F* h_SMS_HH_2Z2g_1Mu_InvMass = (TH3F*)f_SMS_HH2Z_InvMass->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_InvMass");h_SMS_HH_2Z2g_1Mu_InvMass->Sumw2();
  TH3F* h_SMS_HH_2tau2g_1Mu_InvMass = (TH3F*)f_SMS_HH2tau_InvMass->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_InvMass");h_SMS_HH_2tau2g_1Mu_InvMass->Sumw2();
  TH1F* h_SMS_WH_gg_1Mu_Excluded = (TH1F*)h_SMS_WH_gg_1Mu->ProjectionZ("SMS_WH_gg_1Mu_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  TH1F* h_SMS_WH_gg_1Mu_Excluded_InvMass = (TH1F*)h_SMS_WH_gg_1Mu_InvMass->ProjectionZ("SMS_WH_gg_1Mu_Excluded_InvMass",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  TH1F* h_SMS_HH_2W2g_1Mu_Excluded_InvMass = (TH1F*)h_SMS_HH_2W2g_1Mu_InvMass->ProjectionZ("SMS_HH_2W2g_1Mu_Excluded_InvMass",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  TH1F* h_SMS_HH_2Z2g_1Mu_Excluded_InvMass = (TH1F*)h_SMS_HH_2Z2g_1Mu_InvMass->ProjectionZ("SMS_HH_2Z2g_1Mu_Excluded_InvMass",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  TH1F* h_SMS_HH_2tau2g_1Mu_Excluded_InvMass = (TH1F*)h_SMS_HH_2tau2g_1Mu_InvMass->ProjectionZ("SMS_HH_2taug_1Mu_Excluded_InvMass",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_WH_gg_1Mu->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_WH_gg_1Mu->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_WH_gg_1Mu_notExcluded = (TH1F*)h_SMS_WH_gg_1Mu->ProjectionZ("SMS_WH_gg_1Mu_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_ZH_gg_1Mu = (TH3F*)f_SMS_ZH->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_met");h_SMS_ZH_gg_1Mu->Sumw2();
  mChiBin = h_SMS_ZH_gg_1Mu->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_ZH_gg_1Mu->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_ZH_gg_1Mu_Excluded = (TH1F*)h_SMS_ZH_gg_1Mu->ProjectionZ("SMS_ZH_gg_1Mu_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_ZH_gg_1Mu->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_ZH_gg_1Mu->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_ZH_gg_1Mu_notExcluded = (TH1F*)h_SMS_ZH_gg_1Mu->ProjectionZ("SMS_ZH_gg_1Mu_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
 
  TH3F* h_SMS_HH_2b2g_gg_1Mu = (TH3F*)f_SMS_HH_2b2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_met");h_SMS_HH_2b2g_gg_1Mu->Sumw2();
  mChiBin = h_SMS_HH_2b2g_gg_1Mu->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2b2g_gg_1Mu->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2b2g_gg_1Mu_Excluded = (TH1F*)h_SMS_HH_2b2g_gg_1Mu->ProjectionZ("SMS_HH_2b2g_gg_1Mu_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2b2g_gg_1Mu->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2b2g_gg_1Mu->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_HH_2b2g_gg_1Mu_notExcluded = (TH1F*)h_SMS_HH_2b2g_gg_1Mu->ProjectionZ("SMS_HH_2b2g_gg_1Mu_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2W2g_gg_1Mu = (TH3F*)f_SMS_HH_2W2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_met");h_SMS_HH_2W2g_gg_1Mu->Sumw2();
  mChiBin = h_SMS_HH_2W2g_gg_1Mu->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2W2g_gg_1Mu->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2W2g_gg_1Mu_Excluded = (TH1F*)h_SMS_HH_2W2g_gg_1Mu->ProjectionZ("SMS_HH_2W2g_gg_1Mu_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2W2g_gg_1Mu->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2W2g_gg_1Mu->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_HH_2W2g_gg_1Mu_notExcluded = (TH1F*)h_SMS_HH_2W2g_gg_1Mu->ProjectionZ("SMS_HH_2W2g_gg_1Mu_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2Z2g_gg_1Mu = (TH3F*)f_SMS_HH_2Z2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_met");h_SMS_HH_2Z2g_gg_1Mu->Sumw2();
  mChiBin = h_SMS_HH_2Z2g_gg_1Mu->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2Z2g_gg_1Mu->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2Z2g_gg_1Mu_Excluded = (TH1F*)h_SMS_HH_2Z2g_gg_1Mu->ProjectionZ("SMS_HH_2Z2g_gg_1Mu_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2Z2g_gg_1Mu->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2Z2g_gg_1Mu->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_HH_2Z2g_gg_1Mu_notExcluded = (TH1F*)h_SMS_HH_2Z2g_gg_1Mu->ProjectionZ("SMS_HH_2Z2g_gg_1Mu_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  /*
  TH3F* h_SMS_HH_2tau2g_gg_1Mu = (TH3F*)f_SMS_HH_2tau2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_met");h_SMS_HH_2tau2g_gg_1Mu->Sumw2();
  mChiBin = h_SMS_HH_2tau2g_gg_1Mu->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2tau2g_gg_1Mu->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2tau2g_gg_1Mu_Excluded = (TH1F*)h_SMS_HH_2tau2g_gg_1Mu->ProjectionZ("SMS_HH_2tau2g_gg_1Mu_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  */
 

  //
  TH3F* h_SMS_WH_gg_1Mu_MT = (TH3F*)f_SMS_WH->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT");h_SMS_WH_gg_1Mu_MT->Sumw2();
  mChiBin = h_SMS_WH_gg_1Mu_MT->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_WH_gg_1Mu_MT->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_WH_gg_1Mu_MT_Excluded = (TH1F*)h_SMS_WH_gg_1Mu_MT->ProjectionZ("SMS_WH_gg_1Mu_MT_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_WH_gg_1Mu_MT->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_WH_gg_1Mu_MT->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_WH_gg_1Mu_MT_notExcluded = (TH1F*)h_SMS_WH_gg_1Mu_MT->ProjectionZ("SMS_WH_gg_1Mu_MT_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_ZH_gg_1Mu_MT = (TH3F*)f_SMS_ZH->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT");h_SMS_ZH_gg_1Mu_MT->Sumw2();
  mChiBin = h_SMS_ZH_gg_1Mu_MT->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_ZH_gg_1Mu_MT->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_ZH_gg_1Mu_MT_Excluded = (TH1F*)h_SMS_ZH_gg_1Mu_MT->ProjectionZ("SMS_ZH_gg_1Mu_MT_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_ZH_gg_1Mu_MT->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_ZH_gg_1Mu_MT->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_ZH_gg_1Mu_MT_notExcluded = (TH1F*)h_SMS_ZH_gg_1Mu_MT->ProjectionZ("SMS_ZH_gg_1Mu_MT_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
 
  TH3F* h_SMS_HH_2b2g_gg_1Mu_MT = (TH3F*)f_SMS_HH_2b2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT");h_SMS_HH_2b2g_gg_1Mu_MT->Sumw2();
  mChiBin = h_SMS_HH_2b2g_gg_1Mu_MT->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2b2g_gg_1Mu_MT->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2b2g_gg_1Mu_MT_Excluded = (TH1F*)h_SMS_HH_2b2g_gg_1Mu_MT->ProjectionZ("SMS_HH_2b2g_gg_1Mu_MT_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2b2g_gg_1Mu_MT->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2b2g_gg_1Mu_MT->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_HH_2b2g_gg_1Mu_MT_notExcluded = (TH1F*)h_SMS_HH_2b2g_gg_1Mu_MT->ProjectionZ("SMS_HH_2b2g_gg_1Mu_MT_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2W2g_gg_1Mu_MT = (TH3F*)f_SMS_HH_2W2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_met");h_SMS_HH_2W2g_gg_1Mu_MT->Sumw2();
  mChiBin = h_SMS_HH_2W2g_gg_1Mu_MT->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2W2g_gg_1Mu_MT->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2W2g_gg_1Mu_MT_Excluded = (TH1F*)h_SMS_HH_2W2g_gg_1Mu_MT->ProjectionZ("SMS_HH_2W2g_gg_1Mu_MT_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2W2g_gg_1Mu_MT->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2W2g_gg_1Mu_MT->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_HH_2W2g_gg_1Mu_MT_notExcluded = (TH1F*)h_SMS_HH_2W2g_gg_1Mu_MT->ProjectionZ("SMS_HH_2W2g_gg_1Mu_MT_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_HH_2Z2g_gg_1Mu_MT = (TH3F*)f_SMS_HH_2Z2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT");h_SMS_HH_2Z2g_gg_1Mu_MT->Sumw2();
  mChiBin = h_SMS_HH_2Z2g_gg_1Mu_MT->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2Z2g_gg_1Mu_MT->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2Z2g_gg_1Mu_MT_Excluded = (TH1F*)h_SMS_HH_2Z2g_gg_1Mu_MT->ProjectionZ("SMS_HH_2Z2g_gg_1Mu_MT_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2Z2g_gg_1Mu_MT->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2Z2g_gg_1Mu_MT->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_HH_2Z2g_gg_1Mu_MT_notExcluded = (TH1F*)h_SMS_HH_2Z2g_gg_1Mu_MT->ProjectionZ("SMS_HH_2Z2g_gg_1Mu_MT_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  /*
  TH3F* h_SMS_HH_2tau2g_gg_1Mu_MT = (TH3F*)f_SMS_HH_2tau2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT");h_SMS_HH_2tau2g_gg_1Mu_MT->Sumw2();
  mChiBin = h_SMS_HH_2tau2g_gg_1Mu_MT->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2tau2g_gg_1Mu_MT->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2tau2g_gg_1Mu_MT_Excluded = (TH1F*)h_SMS_HH_2tau2g_gg_1Mu_MT->ProjectionZ("SMS_HH_2tau2g_gg_1Mu_MT_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  */
  TH3F* h_SMS_WH_gg_1Mu_MT_noSF = (TH3F*)f_SMS_WH->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT_noScaleFactor");h_SMS_WH_gg_1Mu_MT_noSF->Sumw2();
  float norm_SMS_WH_mu = h_SMS_WH_gg_1Mu_MT_noSF->GetEntries()/h_SMS_WH_gg_1Mu_MT_noSF->Integral(0,-1);
  TH3F* h_SMS_ZH_gg_1Mu_MT_noSF = (TH3F*)f_SMS_ZH->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT_noScaleFactor");h_SMS_ZH_gg_1Mu_MT_noSF->Sumw2();
  float norm_SMS_ZH_mu = h_SMS_ZH_gg_1Mu_MT_noSF->GetEntries()/h_SMS_ZH_gg_1Mu_MT_noSF->Integral(0,-1);
  TH3F* h_SMS_HH_2b2g_gg_1Mu_MT_noSF = (TH3F*)f_SMS_HH_2b2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT_noScaleFactor");h_SMS_HH_2b2g_gg_1Mu_MT_noSF->Sumw2();
  float norm_SMS_HH_2b2g_mu = h_SMS_HH_2b2g_gg_1Mu_MT_noSF->GetEntries()/h_SMS_HH_2b2g_gg_1Mu_MT_noSF->Integral(0,-1);
  TH3F* h_SMS_HH_2W2g_gg_1Mu_MT_noSF = (TH3F*)f_SMS_HH_2W2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT_noScaleFactor");h_SMS_HH_2W2g_gg_1Mu_MT_noSF->Sumw2();
  float norm_SMS_HH_2W2g_mu = h_SMS_HH_2W2g_gg_1Mu_MT_noSF->GetEntries()/h_SMS_HH_2W2g_gg_1Mu_MT_noSF->Integral(0,-1);
  TH3F* h_SMS_HH_2Z2g_gg_1Mu_MT_noSF = (TH3F*)f_SMS_HH_2Z2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT_noScaleFactor");h_SMS_HH_2Z2g_gg_1Mu_MT_noSF->Sumw2();
  float norm_SMS_HH_2Z2g_mu = h_SMS_HH_2Z2g_gg_1Mu_MT_noSF->GetEntries()/h_SMS_HH_2Z2g_gg_1Mu_MT_noSF->Integral(0,-1);
  TH3F* h_SMS_HH_2tau2g_gg_1Mu_MT_noSF = (TH3F*)f_SMS_HH_2tau2g->Get("gg_SMS_Loose_1Mu_0_1Jets_mChi_mBino_MT_noScaleFactor");h_SMS_HH_2tau2g_gg_1Mu_MT_noSF->Sumw2();
  float norm_SMS_HH_2tau2g_mu = h_SMS_HH_2tau2g_gg_1Mu_MT_noSF->GetEntries()/h_SMS_HH_2tau2g_gg_1Mu_MT_noSF->Integral(0,-1);


  h_SMS_WH_gg_1Mu_Excluded->Scale(norm_SMS_WH_mu*scale_SMS_WHExcluded);
  h_SMS_WH_gg_1Mu_Excluded_InvMass->Scale(norm_SMS_WH_mu*scale_SMS_WHExcluded);
  h_SMS_HH_2W2g_1Mu_Excluded_InvMass->Scale(norm_SMS_HH_2W2g_mu*scale_SMS_HH_2W2gExcluded);
  h_SMS_HH_2Z2g_1Mu_Excluded_InvMass->Scale(norm_SMS_HH_2Z2g_mu*scale_SMS_HH_2Z2gExcluded);
  h_SMS_HH_2tau2g_1Mu_Excluded_InvMass->Scale(norm_SMS_HH_2tau2g_mu*scale_SMS_HH_2tau2gExcluded);
  h_SMS_WH_gg_1Mu_notExcluded->Scale(norm_SMS_WH_mu*scale_SMS_WHnotExcluded);
  h_SMS_ZH_gg_1Mu_Excluded->Scale(norm_SMS_ZH_mu*scale_SMS_ZHExcluded);
  h_SMS_ZH_gg_1Mu_notExcluded->Scale(norm_SMS_ZH_mu*scale_SMS_ZHnotExcluded);
  h_SMS_HH_2b2g_gg_1Mu_Excluded->Scale(norm_SMS_HH_2b2g_mu*scale_SMS_HH_2b2gExcluded);
  h_SMS_HH_2W2g_gg_1Mu_Excluded->Scale(norm_SMS_HH_2W2g_mu*scale_SMS_HH_2W2gExcluded);
  h_SMS_HH_2Z2g_gg_1Mu_Excluded->Scale(norm_SMS_HH_2Z2g_mu*scale_SMS_HH_2Z2gExcluded);
  h_SMS_HH_2tau2g_gg_1Mu_Excluded->Scale(norm_SMS_HH_2tau2g_mu*scale_SMS_HH_2tau2gExcluded);

  h_SMS_WH_gg_1Mu_Excluded->SetLineColor(kBlue);h_SMS_WH_gg_1Mu_Excluded->SetLineWidth(3);
  h_SMS_WH_gg_1Mu_Excluded_InvMass->SetLineColor(kRed);h_SMS_WH_gg_1Mu_Excluded_InvMass->SetLineWidth(0);
  h_SMS_WH_gg_1Mu_Excluded_InvMass->SetFillColor(kRed); h_SMS_WH_gg_1Mu_Excluded_InvMass->SetFillStyle(3004);
  h_SMS_HH_2W2g_1Mu_Excluded_InvMass->SetLineColor(kRed);h_SMS_HH_2W2g_1Mu_Excluded_InvMass->SetLineWidth(0);
  h_SMS_HH_2W2g_1Mu_Excluded_InvMass->SetFillColor(kRed); h_SMS_HH_2W2g_1Mu_Excluded_InvMass->SetFillStyle(3004);

  h_SMS_ZH_gg_1Mu_Excluded->SetLineColor(kRed-4);h_SMS_ZH_gg_1Mu_Excluded->SetLineWidth(3);
  h_SMS_WH_gg_1Mu_notExcluded->SetLineColor(kRed);h_SMS_WH_gg_1Mu_notExcluded->SetLineWidth(3);
  h_SMS_ZH_gg_1Mu_notExcluded->SetLineColor(kRed);h_SMS_ZH_gg_1Mu_notExcluded->SetLineWidth(3);
  h_SMS_HH_2b2g_gg_1Mu_Excluded->SetLineColor(kGreen);h_SMS_HH_2b2g_gg_1Mu_Excluded->SetLineWidth(3);
  h_SMS_HH_2W2g_gg_1Mu_Excluded->SetLineColor(kViolet);h_SMS_HH_2W2g_gg_1Mu_Excluded->SetLineWidth(3);
  h_SMS_HH_2Z2g_gg_1Mu_Excluded->SetLineColor(kOrange);h_SMS_HH_2Z2g_gg_1Mu_Excluded->SetLineWidth(3);

  h_SMS_WH_gg_1Mu_MT_Excluded->Scale(norm_SMS_WH_mu*scale_SMS_WHExcluded);
  h_SMS_WH_gg_1Mu_MT_notExcluded->Scale(norm_SMS_WH_mu*scale_SMS_WHnotExcluded);
  h_SMS_ZH_gg_1Mu_MT_Excluded->Scale(norm_SMS_ZH_mu*scale_SMS_ZHExcluded);
  h_SMS_ZH_gg_1Mu_MT_notExcluded->Scale(norm_SMS_ZH_mu*scale_SMS_ZHnotExcluded);
  h_SMS_HH_2b2g_gg_1Mu_MT_Excluded->Scale(norm_SMS_HH_2b2g_mu*scale_SMS_HH_2b2gExcluded);
  h_SMS_HH_2W2g_gg_1Mu_MT_Excluded->Scale(norm_SMS_HH_2W2g_mu*scale_SMS_HH_2W2gExcluded);
  h_SMS_HH_2Z2g_gg_1Mu_MT_Excluded->Scale(norm_SMS_HH_2Z2g_mu*scale_SMS_HH_2Z2gExcluded);
  h_SMS_HH_2tau2g_gg_1Mu_MT_Excluded->Scale(norm_SMS_HH_2tau2g_mu*scale_SMS_HH_2tau2gExcluded);

  h_SMS_WH_gg_1Mu_MT_Excluded->SetLineColor(kBlue);h_SMS_WH_gg_1Mu_MT_Excluded->SetLineWidth(3);
  h_SMS_ZH_gg_1Mu_MT_Excluded->SetLineColor(kRed-4);h_SMS_ZH_gg_1Mu_MT_Excluded->SetLineWidth(3);
  h_SMS_WH_gg_1Mu_MT_notExcluded->SetLineColor(kRed);h_SMS_WH_gg_1Mu_MT_notExcluded->SetLineWidth(3);
  h_SMS_ZH_gg_1Mu_MT_notExcluded->SetLineColor(kRed);h_SMS_ZH_gg_1Mu_MT_notExcluded->SetLineWidth(3);
  h_SMS_HH_2b2g_gg_1Mu_MT_Excluded->SetLineColor(kGreen);h_SMS_HH_2b2g_gg_1Mu_MT_Excluded->SetLineWidth(3);
  h_SMS_HH_2W2g_gg_1Mu_MT_Excluded->SetLineColor(kViolet);h_SMS_HH_2W2g_gg_1Mu_MT_Excluded->SetLineWidth(3);
  h_SMS_HH_2Z2g_gg_1Mu_MT_Excluded->SetLineColor(kOrange);h_SMS_HH_2Z2g_gg_1Mu_MT_Excluded->SetLineWidth(3);


  for(int i=0;i<h_SMS_WH_gg_1Mu_Excluded_InvMass->GetNbinsX();i++){
    if(h_SMS_WH_gg_1Mu_Excluded_InvMass->GetBinLowEdge(i)<120 || h_SMS_WH_gg_1Mu_Excluded_InvMass->GetBinLowEdge(i)>=131)h_SMS_WH_gg_1Mu_Excluded_InvMass->SetBinContent(i,0);
    if(h_SMS_HH_2W2g_1Mu_Excluded_InvMass->GetBinLowEdge(i)<120 || h_SMS_HH_2W2g_1Mu_Excluded_InvMass->GetBinLowEdge(i)>=131)h_SMS_HH_2W2g_1Mu_Excluded_InvMass->SetBinContent(i,0);
    if(h_SMS_HH_2Z2g_1Mu_Excluded_InvMass->GetBinLowEdge(i)<120 || h_SMS_HH_2Z2g_1Mu_Excluded_InvMass->GetBinLowEdge(i)>=131)h_SMS_HH_2Z2g_1Mu_Excluded_InvMass->SetBinContent(i,0);
    if(h_SMS_HH_2tau2g_1Mu_Excluded_InvMass->GetBinLowEdge(i)<120 || h_SMS_HH_2tau2g_1Mu_Excluded_InvMass->GetBinLowEdge(i)>=131)h_SMS_HH_2tau2g_1Mu_Excluded_InvMass->SetBinContent(i,0);
  }
  h_SMS_HH_2W2g_1Mu_Excluded_InvMass->Add(h_SMS_HH_2Z2g_1Mu_Excluded_InvMass);
  h_SMS_HH_2W2g_1Mu_Excluded_InvMass->Add(h_SMS_HH_2tau2g_1Mu_Excluded_InvMass);
  /*
  smsInvMassScale = 6.752/h_SMS_WH_gg_1Mu_Excluded_InvMass->Integral();
  h_SMS_WH_gg_1Mu_Excluded_InvMass->Scale(smsInvMassScale);
  cout<<"smsInvMassScale mu: "<<smsInvMassScale<<"  int: "<<h_SMS_WH_gg_1Mu_Excluded_InvMass->Integral()<< endl;
  */

  //h_SMS_WH_gg_1Mu_Excluded->Add(h_SMS_ZH_gg_1Mu_Excluded);
  //h_SMS_WH_gg_1Mu_notExcluded->Add(h_SMS_ZH_gg_1Mu_notExcluded);

  TH1F* h_SMS_WH_gg_1Mu_Excluded_Rebin = (TH1F*)h_SMS_WH_gg_1Mu_Excluded->Rebin(NmetBinsXtraWide,"h_SMS_WH_gg_1Mu_Excluded_Rebin ",xbinsXtraWide);
  TH1F* h_SMS_WH_gg_1Mu_notExcluded_Rebin = (TH1F*)h_SMS_WH_gg_1Mu_notExcluded->Rebin(NmetBinsXtraWide,"h_SMS_WH_gg_1Mu_notExcluded_Rebin ",xbinsXtraWide);
  TH1F* h_SMS_ZH_gg_1Mu_Excluded_Rebin = (TH1F*)h_SMS_ZH_gg_1Mu_Excluded->Rebin(NmetBinsXtraWide,"h_SMS_ZH_gg_1Mu_Excluded_Rebin ",xbinsXtraWide);
  TH1F* h_SMS_ZH_gg_1Mu_notExcluded_Rebin = (TH1F*)h_SMS_ZH_gg_1Mu_notExcluded->Rebin(NmetBinsXtraWide,"h_SMS_ZH_gg_1Mu_notExcluded_Rebin ",xbinsXtraWide);
    TH1F* h_SMS_HH_2b2g_gg_1Mu_Excluded_Rebin = (TH1F*)h_SMS_HH_2b2g_gg_1Mu_Excluded->Rebin(NmetBinsXtraWide,"h_SMS_HH_2b2g_gg_1Mu_Excluded_Rebin ",xbinsXtraWide);
    TH1F* h_SMS_HH_2W2g_gg_1Mu_Excluded_Rebin = (TH1F*)h_SMS_HH_2W2g_gg_1Mu_Excluded->Rebin(NmetBinsXtraWide,"h_SMS_HH_2W2g_gg_1Mu_Excluded_Rebin ",xbinsXtraWide);
    TH1F* h_SMS_HH_2Z2g_gg_1Mu_Excluded_Rebin = (TH1F*)h_SMS_HH_2Z2g_gg_1Mu_Excluded->Rebin(NmetBinsXtraWide,"h_SMS_HH_2Z2g_gg_1Mu_Excluded_Rebin ",xbinsXtraWide);
    TH1F* h_SMS_HH_2tau2g_gg_1Mu_Excluded_Rebin = (TH1F*)h_SMS_HH_2tau2g_gg_1Mu_Excluded->Rebin(NmetBinsXtraWide,"h_SMS_HH_2tau2g_gg_1Mu_Excluded_Rebin ",xbinsXtraWide);

  TH1F* h_SMS_WH_gg_1Mu_MT_Excluded_Rebin = (TH1F*)h_SMS_WH_gg_1Mu_MT_Excluded->Rebin(nMTbins,"h_SMS_WH_gg_1Mu_MT_Excluded_Rebin ",MTbins);
  TH1F* h_SMS_WH_gg_1Mu_MT_notExcluded_Rebin = (TH1F*)h_SMS_WH_gg_1Mu_MT_notExcluded->Rebin(nMTbins,"h_SMS_WH_gg_1Mu_MT_notExcluded_Rebin ",MTbins);
  TH1F* h_SMS_ZH_gg_1Mu_MT_Excluded_Rebin = (TH1F*)h_SMS_ZH_gg_1Mu_MT_Excluded->Rebin(nMTbins,"h_SMS_ZH_gg_1Mu_MT_Excluded_Rebin ",MTbins);
  TH1F* h_SMS_ZH_gg_1Mu_MT_notExcluded_Rebin = (TH1F*)h_SMS_ZH_gg_1Mu_MT_notExcluded->Rebin(nMTbins,"h_SMS_ZH_gg_1Mu_MT_notExcluded_Rebin ",MTbins);
    TH1F* h_SMS_HH_2b2g_gg_1Mu_MT_Excluded_Rebin = (TH1F*)h_SMS_HH_2b2g_gg_1Mu_MT_Excluded->Rebin(nMTbins,"h_SMS_HH_2b2g_gg_1Mu_MT_Excluded_Rebin ",MTbins);
    TH1F* h_SMS_HH_2W2g_gg_1Mu_MT_Excluded_Rebin = (TH1F*)h_SMS_HH_2W2g_gg_1Mu_MT_Excluded->Rebin(nMTbins,"h_SMS_HH_2W2g_gg_1Mu_MT_Excluded_Rebin ",MTbins);
    TH1F* h_SMS_HH_2Z2g_gg_1Mu_MT_Excluded_Rebin = (TH1F*)h_SMS_HH_2Z2g_gg_1Mu_MT_Excluded->Rebin(nMTbins,"h_SMS_HH_2Z2g_gg_1Mu_MT_Excluded_Rebin ",MTbins);
    TH1F* h_SMS_HH_2tau2g_gg_1Mu_MT_Excluded_Rebin = (TH1F*)h_SMS_HH_2tau2g_gg_1Mu_MT_Excluded->Rebin(nMTbins,"h_SMS_HH_2tau2g_gg_1Mu_MT_Excluded_Rebin ",MTbins);
/*
  SMSoverflowbin=h_SMS_WH_gg_1Mu_Excluded_Rebin->FindBin(metPlotXmaxXtraWide+1);SMSlastbin=h_SMS_WH_gg_1Mu_Excluded_Rebin->FindBin(metPlotXmaxXtraWide-1);
  SMSoverflowerr=0;SMSlastbinerr=0;
  SMSoverflow = h_SMS_WH_gg_1Mu_Excluded_Rebin->IntegralAndError(SMSoverflowbin,-1,SMSoverflowerr);
  SMSlast = h_SMS_WH_gg_1Mu_Excluded_Rebin->IntegralAndError(SMSlastbin,SMSlastbin,SMSlastbinerr);
  h_SMS_WH_gg_1Mu_Excluded_Rebin->SetBinContent(SMSlastbin,SMSlast+SMSoverflow);
  h_SMS_WH_gg_1Mu_Excluded_Rebin->SetBinError(SMSlastbin,sqrt(SMSlastbinerr*SMSlastbinerr+SMSoverflow*SMSoverflow));
  h_SMS_WH_gg_1Mu_Excluded_Rebin->SetBinContent(SMSoverflowbin,0);h_SMS_WH_gg_1Mu_Excluded_Rebin->SetBinError(SMSoverflowbin,0);
  SMSoverflow = h_SMS_WH_gg_1Mu_notExcluded_Rebin->IntegralAndError(SMSoverflowbin,-1,SMSoverflowerr);
  SMSlast = h_SMS_WH_gg_1Mu_notExcluded_Rebin->IntegralAndError(SMSlastbin,SMSlastbin,SMSlastbinerr);
  h_SMS_WH_gg_1Mu_notExcluded_Rebin->SetBinContent(SMSlastbin,SMSlast+SMSoverflow);
  h_SMS_WH_gg_1Mu_notExcluded_Rebin->SetBinError(SMSlastbin,sqrt(SMSlastbinerr*SMSlastbinerr+SMSoverflow*SMSoverflow));
  h_SMS_WH_gg_1Mu_notExcluded_Rebin->SetBinContent(SMSoverflowbin,0);h_SMS_WH_gg_1Mu_notExcluded_Rebin->SetBinError(SMSoverflowbin,0);
  
  SMSoverflow = h_SMS_ZH_gg_1Mu_Excluded_Rebin->IntegralAndError(SMSoverflowbin,-1,SMSoverflowerr);
  SMSlast = h_SMS_ZH_gg_1Mu_Excluded_Rebin->IntegralAndError(SMSlastbin,SMSlastbin,SMSlastbinerr);
  h_SMS_ZH_gg_1Mu_Excluded_Rebin->SetBinContent(SMSlastbin,SMSlast+SMSoverflow);
  h_SMS_ZH_gg_1Mu_Excluded_Rebin->SetBinError(SMSlastbin,sqrt(SMSlastbinerr*SMSlastbinerr+SMSoverflow*SMSoverflow));
  h_SMS_ZH_gg_1Mu_Excluded_Rebin->SetBinContent(SMSoverflowbin,0);h_SMS_ZH_gg_1Mu_Excluded_Rebin->SetBinError(SMSoverflowbin,0);
  SMSoverflow = h_SMS_ZH_gg_1Mu_notExcluded_Rebin->IntegralAndError(SMSoverflowbin,-1,SMSoverflowerr);
  SMSlast = h_SMS_ZH_gg_1Mu_notExcluded_Rebin->IntegralAndError(SMSlastbin,SMSlastbin,SMSlastbinerr);
  h_SMS_ZH_gg_1Mu_notExcluded_Rebin->SetBinContent(SMSlastbin,SMSlast+SMSoverflow);
  h_SMS_ZH_gg_1Mu_notExcluded_Rebin->SetBinError(SMSlastbin,sqrt(SMSlastbinerr*SMSlastbinerr+SMSoverflow*SMSoverflow));
  h_SMS_ZH_gg_1Mu_notExcluded_Rebin->SetBinContent(SMSoverflowbin,0);h_SMS_ZH_gg_1Mu_notExcluded_Rebin->SetBinError(SMSoverflowbin,0);
  */
  AddOverflowToLastBin(h_SMS_WH_gg_1Mu_Excluded_Rebin);AddOverflowToLastBin(h_SMS_WH_gg_1Mu_notExcluded_Rebin);
  AddOverflowToLastBin(h_SMS_ZH_gg_1Mu_Excluded_Rebin);AddOverflowToLastBin(h_SMS_ZH_gg_1Mu_notExcluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2b2g_gg_1Mu_Excluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2W2g_gg_1Mu_Excluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2Z2g_gg_1Mu_Excluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2tau2g_gg_1Mu_Excluded_Rebin);
  DivideBy15gev(h_SMS_WH_gg_1Mu_Excluded_Rebin);DivideBy15gev(h_SMS_WH_gg_1Mu_notExcluded_Rebin);
  DivideBy15gev(h_SMS_ZH_gg_1Mu_Excluded_Rebin);DivideBy15gev(h_SMS_ZH_gg_1Mu_notExcluded_Rebin);
  DivideBy15gev(h_SMS_HH_2b2g_gg_1Mu_Excluded_Rebin);
  DivideBy15gev(h_SMS_HH_2W2g_gg_1Mu_Excluded_Rebin);
  DivideBy15gev(h_SMS_HH_2Z2g_gg_1Mu_Excluded_Rebin);
  DivideBy15gev(h_SMS_HH_2tau2g_gg_1Mu_Excluded_Rebin);
  
  AddOverflowToLastBin(h_SMS_WH_gg_1Mu_MT_Excluded_Rebin);AddOverflowToLastBin(h_SMS_WH_gg_1Mu_MT_notExcluded_Rebin);
  AddOverflowToLastBin(h_SMS_ZH_gg_1Mu_MT_Excluded_Rebin);AddOverflowToLastBin(h_SMS_ZH_gg_1Mu_MT_notExcluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2b2g_gg_1Mu_MT_Excluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2W2g_gg_1Mu_MT_Excluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2Z2g_gg_1Mu_MT_Excluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2tau2g_gg_1Mu_MT_Excluded_Rebin);
  DivideBy30gev(h_SMS_WH_gg_1Mu_MT_Excluded_Rebin);DivideBy30gev(h_SMS_WH_gg_1Mu_MT_notExcluded_Rebin);
  DivideBy30gev(h_SMS_ZH_gg_1Mu_MT_Excluded_Rebin);DivideBy30gev(h_SMS_ZH_gg_1Mu_MT_notExcluded_Rebin);
  DivideBy30gev(h_SMS_HH_2b2g_gg_1Mu_MT_Excluded_Rebin);
  DivideBy30gev(h_SMS_HH_2W2g_gg_1Mu_MT_Excluded_Rebin);
  DivideBy30gev(h_SMS_HH_2Z2g_gg_1Mu_MT_Excluded_Rebin);
  DivideBy30gev(h_SMS_HH_2tau2g_gg_1Mu_MT_Excluded_Rebin);
 
  //now leftovers
  TH3F* h_SMS_WH_gg_leftovers = (TH3F*)f_SMS_WH->Get("gg_SMS_Loose_lt2b_noMu_noEle_lt2jets_mChi_mBino_met");h_SMS_WH_gg_leftovers->Sumw2();
  mChiBin = h_SMS_WH_gg_leftovers->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_WH_gg_leftovers->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_WH_gg_leftovers_Excluded = (TH1F*)h_SMS_WH_gg_leftovers->ProjectionZ("SMS_WH_gg_leftovers_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_WH_gg_leftovers->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_WH_gg_leftovers->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_WH_gg_leftovers_notExcluded = (TH1F*)h_SMS_WH_gg_leftovers->ProjectionZ("SMS_WH_gg_leftovers_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");

  TH3F* h_SMS_ZH_gg_leftovers = (TH3F*)f_SMS_ZH->Get("gg_SMS_Loose_lt2b_noMu_noEle_lt2jets_mChi_mBino_met");h_SMS_ZH_gg_leftovers->Sumw2();
  mChiBin = h_SMS_ZH_gg_leftovers->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_ZH_gg_leftovers->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_ZH_gg_leftovers_Excluded = (TH1F*)h_SMS_ZH_gg_leftovers->ProjectionZ("SMS_ZH_gg_leftovers_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_ZH_gg_leftovers->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_ZH_gg_leftovers->GetYaxis()->FindBin(SMSnotExcBino);
  TH1F* h_SMS_ZH_gg_leftovers_notExcluded = (TH1F*)h_SMS_ZH_gg_leftovers->ProjectionZ("SMS_ZH_gg_leftovers_notExcluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
 
  TH3F* h_SMS_HH_2b2g_gg_leftovers = (TH3F*)f_SMS_HH_2b2g->Get("gg_SMS_Loose_lt2b_noMu_noEle_lt2jets_mChi_mBino_met");h_SMS_HH_2b2g_gg_leftovers->Sumw2();
  mChiBin = h_SMS_HH_2b2g_gg_leftovers->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2b2g_gg_leftovers->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2b2g_gg_leftovers_Excluded = (TH1F*)h_SMS_HH_2b2g_gg_leftovers->ProjectionZ("SMS_HH_2b2g_gg_leftovers_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2b2g_gg_leftovers->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2b2g_gg_leftovers->GetYaxis()->FindBin(SMSnotExcBino);

  TH3F* h_SMS_HH_2W2g_gg_leftovers = (TH3F*)f_SMS_HH_2W2g->Get("gg_SMS_Loose_lt2b_noMu_noEle_lt2jets_mChi_mBino_met");h_SMS_HH_2W2g_gg_leftovers->Sumw2();
  mChiBin = h_SMS_HH_2W2g_gg_leftovers->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2W2g_gg_leftovers->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2W2g_gg_leftovers_Excluded = (TH1F*)h_SMS_HH_2W2g_gg_leftovers->ProjectionZ("SMS_HH_2W2g_gg_leftovers_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2W2g_gg_leftovers->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2W2g_gg_leftovers->GetYaxis()->FindBin(SMSnotExcBino);

  TH3F* h_SMS_HH_2Z2g_gg_leftovers = (TH3F*)f_SMS_HH_2Z2g->Get("gg_SMS_Loose_lt2b_noMu_noEle_lt2jets_mChi_mBino_met");h_SMS_HH_2Z2g_gg_leftovers->Sumw2();
  mChiBin = h_SMS_HH_2Z2g_gg_leftovers->GetXaxis()->FindBin(SMSExcChi),mBinoBin = h_SMS_HH_2Z2g_gg_leftovers->GetYaxis()->FindBin(SMSExcBino);
  TH1F* h_SMS_HH_2Z2g_gg_leftovers_Excluded = (TH1F*)h_SMS_HH_2Z2g_gg_leftovers->ProjectionZ("SMS_HH_2Z2g_gg_leftovers_Excluded",mChiBin,mChiBin,mBinoBin,mBinoBin,"eo");
  mChiBin = h_SMS_HH_2Z2g_gg_leftovers->GetXaxis()->FindBin(SMSnotExcChi),mBinoBin = h_SMS_HH_2Z2g_gg_leftovers->GetYaxis()->FindBin(SMSnotExcBino);

  h_SMS_WH_gg_leftovers_Excluded->Scale(scale_SMS_WHExcluded);
  h_SMS_WH_gg_leftovers_notExcluded->Scale(scale_SMS_WHnotExcluded);
  h_SMS_ZH_gg_leftovers_Excluded->Scale(scale_SMS_ZHExcluded);
  h_SMS_ZH_gg_leftovers_notExcluded->Scale(scale_SMS_ZHnotExcluded);
  h_SMS_HH_2b2g_gg_leftovers_Excluded->Scale(scale_SMS_HH_2b2gExcluded);
  h_SMS_HH_2W2g_gg_leftovers_Excluded->Scale(scale_SMS_HH_2W2gExcluded);
  h_SMS_HH_2Z2g_gg_leftovers_Excluded->Scale(scale_SMS_HH_2Z2gExcluded);

  h_SMS_WH_gg_leftovers_Excluded->SetLineColor(kBlue);h_SMS_WH_gg_leftovers_Excluded->SetLineWidth(3);
  h_SMS_ZH_gg_leftovers_Excluded->SetLineColor(kRed-4);h_SMS_ZH_gg_leftovers_Excluded->SetLineWidth(3);
  h_SMS_WH_gg_leftovers_notExcluded->SetLineColor(kRed);h_SMS_WH_gg_leftovers_notExcluded->SetLineWidth(3);
  h_SMS_ZH_gg_leftovers_notExcluded->SetLineColor(kRed);h_SMS_ZH_gg_leftovers_notExcluded->SetLineWidth(3);
  h_SMS_HH_2b2g_gg_leftovers_Excluded->SetLineColor(kGreen);h_SMS_HH_2b2g_gg_leftovers_Excluded->SetLineWidth(3);
  h_SMS_HH_2W2g_gg_leftovers_Excluded->SetLineColor(kViolet);h_SMS_HH_2W2g_gg_leftovers_Excluded->SetLineWidth(3);
  h_SMS_HH_2Z2g_gg_leftovers_Excluded->SetLineColor(kOrange);h_SMS_HH_2Z2g_gg_leftovers_Excluded->SetLineWidth(3);


  //h_SMS_WH_gg_leftovers_Excluded->Add(h_SMS_ZH_gg_leftovers_Excluded);
  //h_SMS_WH_gg_leftovers_notExcluded->Add(h_SMS_ZH_gg_leftovers_notExcluded);


  TH1F* h_SMS_WH_gg_leftovers_Excluded_Rebin = (TH1F*)h_SMS_WH_gg_leftovers_Excluded->Rebin(NmetBinsWide,"h_SMS_WH_gg_leftovers_Excluded_Rebin ",xbinsWide);
  TH1F* h_SMS_WH_gg_leftovers_notExcluded_Rebin = (TH1F*)h_SMS_WH_gg_leftovers_notExcluded->Rebin(NmetBinsWide,"h_SMS_WH_gg_leftovers_notExcluded_Rebin ",xbinsWide);
  TH1F* h_SMS_ZH_gg_leftovers_Excluded_Rebin = (TH1F*)h_SMS_ZH_gg_leftovers_Excluded->Rebin(NmetBinsWide,"h_SMS_ZH_gg_leftovers_Excluded_Rebin ",xbinsWide);
  TH1F* h_SMS_ZH_gg_leftovers_notExcluded_Rebin = (TH1F*)h_SMS_ZH_gg_leftovers_notExcluded->Rebin(NmetBinsWide,"h_SMS_ZH_gg_leftovers_notExcluded_Rebin ",xbinsWide);
  TH1F* h_SMS_HH_2b2g_gg_leftovers_Excluded_Rebin = (TH1F*)h_SMS_HH_2b2g_gg_leftovers_Excluded->Rebin(NmetBinsWide,"h_SMS_HH_2b2g_gg_leftovers_Excluded_Rebin ",xbinsWide);
  TH1F* h_SMS_HH_2W2g_gg_leftovers_Excluded_Rebin = (TH1F*)h_SMS_HH_2W2g_gg_leftovers_Excluded->Rebin(NmetBinsWide,"h_SMS_HH_2W2g_gg_leftovers_Excluded_Rebin ",xbinsWide);
  TH1F* h_SMS_HH_2Z2g_gg_leftovers_Excluded_Rebin = (TH1F*)h_SMS_HH_2Z2g_gg_leftovers_Excluded->Rebin(NmetBinsWide,"h_SMS_HH_2Z2g_gg_leftovers_Excluded_Rebin ",xbinsWide);

  /*
  SMSoverflowbin=h_SMS_WH_gg_leftovers_Excluded_Rebin->FindBin(metPlotXmaxWide+1);SMSlastbin=h_SMS_WH_gg_leftovers_Excluded_Rebin->FindBin(metPlotXmaxWide-1);
  SMSoverflowerr=0;SMSlastbinerr=0;
  SMSoverflow = h_SMS_WH_gg_leftovers_Excluded_Rebin->IntegralAndError(SMSoverflowbin,-1,SMSoverflowerr);
  SMSlast = h_SMS_WH_gg_leftovers_Excluded_Rebin->IntegralAndError(SMSlastbin,SMSlastbin,SMSlastbinerr);
  h_SMS_WH_gg_leftovers_Excluded_Rebin->SetBinContent(SMSlastbin,SMSlast+SMSoverflow);
  h_SMS_WH_gg_leftovers_Excluded_Rebin->SetBinError(SMSlastbin,sqrt(SMSlastbinerr*SMSlastbinerr+SMSoverflow*SMSoverflow));
  h_SMS_WH_gg_leftovers_Excluded_Rebin->SetBinContent(SMSoverflowbin,0);h_SMS_WH_gg_leftovers_Excluded_Rebin->SetBinError(SMSoverflowbin,0);
  SMSoverflow = h_SMS_WH_gg_leftovers_notExcluded_Rebin->IntegralAndError(SMSoverflowbin,-1,SMSoverflowerr);
  SMSlast = h_SMS_WH_gg_leftovers_notExcluded_Rebin->IntegralAndError(SMSlastbin,SMSlastbin,SMSlastbinerr);
  h_SMS_WH_gg_leftovers_notExcluded_Rebin->SetBinContent(SMSlastbin,SMSlast+SMSoverflow);
  h_SMS_WH_gg_leftovers_notExcluded_Rebin->SetBinError(SMSlastbin,sqrt(SMSlastbinerr*SMSlastbinerr+SMSoverflow*SMSoverflow));
  h_SMS_WH_gg_leftovers_notExcluded_Rebin->SetBinContent(SMSoverflowbin,0);h_SMS_WH_gg_leftovers_notExcluded_Rebin->SetBinError(SMSoverflowbin,0);
  
  SMSoverflow = h_SMS_ZH_gg_leftovers_Excluded_Rebin->IntegralAndError(SMSoverflowbin,-1,SMSoverflowerr);
  SMSlast = h_SMS_ZH_gg_leftovers_Excluded_Rebin->IntegralAndError(SMSlastbin,SMSlastbin,SMSlastbinerr);
  h_SMS_ZH_gg_leftovers_Excluded_Rebin->SetBinContent(SMSlastbin,SMSlast+SMSoverflow);
  h_SMS_ZH_gg_leftovers_Excluded_Rebin->SetBinError(SMSlastbin,sqrt(SMSlastbinerr*SMSlastbinerr+SMSoverflow*SMSoverflow));
  h_SMS_ZH_gg_leftovers_Excluded_Rebin->SetBinContent(SMSoverflowbin,0);h_SMS_ZH_gg_leftovers_Excluded_Rebin->SetBinError(SMSoverflowbin,0);
  SMSoverflow = h_SMS_ZH_gg_leftovers_notExcluded_Rebin->IntegralAndError(SMSoverflowbin,-1,SMSoverflowerr);
  SMSlast = h_SMS_ZH_gg_leftovers_notExcluded_Rebin->IntegralAndError(SMSlastbin,SMSlastbin,SMSlastbinerr);
  h_SMS_ZH_gg_leftovers_notExcluded_Rebin->SetBinContent(SMSlastbin,SMSlast+SMSoverflow);
  h_SMS_ZH_gg_leftovers_notExcluded_Rebin->SetBinError(SMSlastbin,sqrt(SMSlastbinerr*SMSlastbinerr+SMSoverflow*SMSoverflow));
  h_SMS_ZH_gg_leftovers_notExcluded_Rebin->SetBinContent(SMSoverflowbin,0);h_SMS_ZH_gg_leftovers_notExcluded_Rebin->SetBinError(SMSoverflowbin,0);
  */

  AddOverflowToLastBin(h_SMS_WH_gg_leftovers_Excluded_Rebin);AddOverflowToLastBin(h_SMS_WH_gg_leftovers_notExcluded_Rebin);
  AddOverflowToLastBin(h_SMS_ZH_gg_leftovers_Excluded_Rebin);AddOverflowToLastBin(h_SMS_ZH_gg_leftovers_notExcluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2b2g_gg_leftovers_Excluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2W2g_gg_leftovers_Excluded_Rebin);
  AddOverflowToLastBin(h_SMS_HH_2Z2g_gg_leftovers_Excluded_Rebin);
  DivideByBinWidth(h_SMS_WH_gg_leftovers_Excluded_Rebin);DivideByBinWidth(h_SMS_WH_gg_leftovers_notExcluded_Rebin);
  DivideByBinWidth(h_SMS_ZH_gg_leftovers_Excluded_Rebin);DivideByBinWidth(h_SMS_ZH_gg_leftovers_notExcluded_Rebin);
  DivideByBinWidth(h_SMS_HH_2b2g_gg_leftovers_Excluded_Rebin);
  DivideByBinWidth(h_SMS_HH_2W2g_gg_leftovers_Excluded_Rebin);
  DivideByBinWidth(h_SMS_HH_2Z2g_gg_leftovers_Excluded_Rebin);
  /*
  //TH1F* h_aaW130_ggMet = (TH1F*)f_aaW130->Get("ggMet");if(h_aaW130_ggMet->GetSumw2N()==0)h_aaW130_ggMet->Sumw2();
  float nEvents_aaW130=120000.,nEvents_aaW275=30000.;
  //float xSec_aaW130=5.841,xSec_aaW275=0.294;//https://twiki.cern.ch/twiki/bin/viewauth/CMS/Electrohiggs
  float xSec_aaW130=4.12,xSec_aaW275=0.213;//https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections8TeVcharginoneutralinoCMS#Chargino_Neutralino
  float scale_aaW130=L_int*xSec_aaW130*BR_aaW/nEvents_aaW130,scale_aaW275=L_int*xSec_aaW275*BR_aaW/nEvents_aaW275;
  //h_aaW130_ggMet->Scale(scale_aaW130);
  //h_aaW130_ggMet->Draw();
  //c1->Print("Plots/Higgs/aaW130_ggMet.png");
  */
  TH2F* h_ggMetData = (TH2F*)fin->Get("ggMetVsInvarMass_Loose");h_ggMetData->Sumw2();h_ggMetData->SetMarkerSize(.5);
  TH2F* h_ggMetData_JetReq = (TH2F*)fin->Get("ggMetVsInvarMass_Loose_JetReq");h_ggMetData_JetReq->Sumw2();h_ggMetData_JetReq->SetMarkerSize(.5);
  int bin30=h_ggMetData->GetYaxis()->FindBin(30);

  TH1D* h_ggInvarMassMET = (TH1D*)h_ggMetData->ProjectionX("h_ggInvarMassMET",0,-1,"eo");
  TH1D* h_ggInvarMassMET_JetReq = (TH1D*)h_ggMetData_JetReq->ProjectionX("h_ggInvarMassMET_JetReq",0,-1,"eo");
  TH1D* h_ggInvarMassMET30 = (TH1D*)h_ggMetData->ProjectionX("h_ggInvarMassMET30",bin30,-1,"eo");
  TH1D* h_ggInvarMassMET30_JetReq = (TH1D*)h_ggMetData_JetReq->ProjectionX("h_ggInvarMassMET30_JetReq",bin30,-1,"eo");
  /*
    TH1F* h_ggInvarMassMET = (TH1F*)fin->Get("ggInvarMass");h_ggInvarMassMET->Sumw2();h_ggInvarMassMET->SetMarkerSize(.5);
    TH1F* h_ggInvarMassMET_JetReq=(TH1F*)fin->Get("ggInvarMass_JetReq");h_ggInvarMassMET_JetReq->Sumw2();
    TH1F* h_ggInvarMassMET30 = (TH1F*)fin->Get("h_ggInvarMassMET30MVAcorr");
    TH1F* h_ggInvarMassMET30_JetReq = (TH1F*)fin->Get("h_ggInvarMassMET30_JetReq");
  */  
  h_ggInvarMassMET->Rebin(2);h_ggInvarMassMET_JetReq->Rebin(2);
  h_ggInvarMassMET30->Rebin(2);h_ggInvarMassMET30_JetReq->Rebin(2);
  //MakeBlindInvMass(h_ggInvarMassMET);
  //MakeBlindInvMass(h_ggInvarMassMET30); 
  TH1F* h_ggInvMassSubtract = (TH1F*)h_ggInvarMassMET->Clone();//h_ggInvarMassMET calling ggInvarMassMVAcorrVertexCorr makes problems later...??
  TH1F* h_ggInvMassSubtract30 = (TH1F*)h_ggInvarMassMET30->Clone();//h_ggInvarMassMET30 calling ggInvarMassMVAcorrVertexCorr makes problems later...??
  
  h_ggInvarMassMET->GetXaxis()->SetRangeUser(20,250);h_ggInvarMassMET->GetYaxis()->SetRangeUser(0,7999);
  h_ggInvarMassMET->SetTitle("");
  h_ggInvarMassMET->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV]");h_ggInvarMassMET->GetYaxis()->SetTitle("Events");
  h_ggInvarMassMET->Draw();
  TLine low(101,0,101,100);low.SetLineColor(kRed);low.SetLineWidth(3);
  low.DrawLine(101,0,101,7999);low.DrawLine(118,0,118,7999);low.DrawLine(133,0,133,7999);low.DrawLine(163,0,163,7999);
  TLine high(101,0,101,100);high.SetLineColor(kBlue);high.SetLineWidth(3);
  high.DrawLine(120,0,120,7999);high.DrawLine(131,0,131,7999);
  c1->Print("Plots/Higgs/ggInvarMassData.png");
  c1->Print("Plots/Higgs/ggInvarMassData.pdf");

  h_ggInvarMassMET->GetXaxis()->SetRangeUser(100,164.9);h_ggInvarMassMET->GetYaxis()->SetRangeUser(0,4400);//was 6550
  h_ggInvarMassMET->Draw();
  low.DrawLine(101,0,101,6550);low.DrawLine(118,0,118,6550);low.DrawLine(133,0,133,6550);low.DrawLine(163,0,163,6550);
  high.DrawLine(120,0,120,6550);high.DrawLine(131,0,131,6550);
  c1->Print("Plots/Higgs/ggInvarMassData_zoom.png");
  c1->Print("Plots/Higgs/ggInvarMassData_zoom.pdf");
  h_ggInvarMassMET->GetXaxis()->SetRangeUser(0,1000);


  cout<<"h_ggInvarMassMET integral: "<<h_ggInvarMassMET->Integral()<<endl<<"h_ggInvMassSubtract integral: "<<h_ggInvMassSubtract->Integral()<<endl;
  cout<<"h_ggInvarMassMET30 integral: "<<h_ggInvarMassMET30->Integral()<<endl<<"h_ggInvMassSubtract30 integral: "<<h_ggInvMassSubtract30->Integral()<<endl;

  TH2F* h_ggHgg = (TH2F*)f_ggHgg->Get("ggMetVsInvarMass_Loose");if(h_ggHgg->GetSumw2N()==0)h_ggHgg->Sumw2();h_ggHgg->SetMarkerSize(.5);
  TH2F* h_WZHgg = (TH2F*)f_WZHgg->Get("ggMetVsInvarMass_Loose");if(h_WZHgg->GetSumw2N()==0)h_WZHgg->Sumw2();h_WZHgg->SetMarkerSize(.5);
  TH2F* h_TTHgg = (TH2F*)f_TTHgg->Get("ggMetVsInvarMass_Loose");if(h_TTHgg->GetSumw2N()==0)h_TTHgg->Sumw2();h_TTHgg->SetMarkerSize(.5);
  TH2F* h_VBFHgg = (TH2F*)f_VBFHgg->Get("ggMetVsInvarMass_Loose");if(h_VBFHgg->GetSumw2N()==0)h_VBFHgg->Sumw2();h_VBFHgg->SetMarkerSize(.5);
  TH2F* h_ggHgg_uncorr = (TH2F*)f_ggHgg->Get("ggMetVsUncorrectedInvarMass_Loose");if(h_ggHgg_uncorr->GetSumw2N()==0)h_ggHgg_uncorr->Sumw2();h_ggHgg_uncorr->SetMarkerSize(.5);
  TH2F* h_WZHgg_uncorr = (TH2F*)f_WZHgg->Get("ggMetVsUncorrectedInvarMass_Loose");if(h_WZHgg_uncorr->GetSumw2N()==0)h_WZHgg_uncorr->Sumw2();h_WZHgg_uncorr->SetMarkerSize(.5);
  TH2F* h_TTHgg_uncorr = (TH2F*)f_TTHgg->Get("ggMetVsUncorrectedInvarMass_Loose");if(h_TTHgg_uncorr->GetSumw2N()==0)h_TTHgg_uncorr->Sumw2();h_TTHgg_uncorr->SetMarkerSize(.5);
  TH2F* h_VBFHgg_uncorr = (TH2F*)f_VBFHgg->Get("ggMetVsUncorrectedInvarMass_Loose");if(h_VBFHgg_uncorr->GetSumw2N()==0)h_VBFHgg_uncorr->Sumw2();h_VBFHgg_uncorr->SetMarkerSize(.5);
  /*
  TH2F* aaW130 = (TH2F*)f_aaW130->Get("ggMetVsInvarMass_Loose");
  TH2F* aaW275 = (TH2F*)f_aaW275->Get("ggMetVsInvarMass_Loose");
  aaW130->Scale(scale_aaW130);
  aaW130->SetLineColor(kBlue);aaW130->SetLineWidth(3);
  aaW275->Scale(scale_aaW275);
  aaW275->SetLineColor(kRed);aaW275->SetLineWidth(3);
  */
  TH2F* h_ggHgg_JetReq = (TH2F*)f_ggHgg->Get("ggMetVsInvarMass_Loose_JetReq");if(h_ggHgg_JetReq->GetSumw2N()==0)h_ggHgg_JetReq->Sumw2();h_ggHgg_JetReq->SetMarkerSize(.5);
  TH2F* h_WZHgg_JetReq = (TH2F*)f_WZHgg->Get("ggMetVsInvarMass_Loose_JetReq");if(h_WZHgg_JetReq->GetSumw2N()==0)h_WZHgg_JetReq->Sumw2();h_WZHgg_JetReq->SetMarkerSize(.5);
  TH2F* h_TTHgg_JetReq = (TH2F*)f_TTHgg->Get("ggMetVsInvarMass_Loose_JetReq");if(h_TTHgg_JetReq->GetSumw2N()==0)h_TTHgg_JetReq->Sumw2();h_TTHgg_JetReq->SetMarkerSize(.5);
  TH2F* h_VBFHgg_JetReq = (TH2F*)f_VBFHgg->Get("ggMetVsInvarMass_Loose_JetReq");if(h_VBFHgg_JetReq->GetSumw2N()==0)h_VBFHgg_JetReq->Sumw2();h_VBFHgg_JetReq->SetMarkerSize(.5);

  TH1D* ggHggNoMET = (TH1D*)h_ggHgg->ProjectionX("ggHggNoMET",0,-1,"eo");
  TH1D* WZHggNoMET = (TH1D*)h_WZHgg->ProjectionX("WZHggNoMET",0,-1,"eo");
  TH1D* TTHggNoMET = (TH1D*)h_TTHgg->ProjectionX("TTHggNoMET",0,-1,"eo");
  TH1D* VBFHggNoMET = (TH1D*)h_VBFHgg->ProjectionX("VBFHggNoMET",0,-1,"eo");
  TH1D* ggHggNoMET_uncorr = (TH1D*)h_ggHgg_uncorr->ProjectionX("ggHggNoMET_uncorr",0,-1,"eo");
  TH1D* WZHggNoMET_uncorr = (TH1D*)h_WZHgg_uncorr->ProjectionX("WZHggNoMET_uncorr",0,-1,"eo");
  TH1D* TTHggNoMET_uncorr = (TH1D*)h_TTHgg_uncorr->ProjectionX("TTHggNoMET_uncorr",0,-1,"eo");
  TH1D* VBFHggNoMET_uncorr = (TH1D*)h_VBFHgg_uncorr->ProjectionX("VBFHggNoMET_uncorr",0,-1,"eo");
  TH1D* ggHggMET30 = (TH1D*)h_ggHgg->ProjectionX("ggHggMET30",bin30,-1,"eo");
  TH1D* WZHggMET30 = (TH1D*)h_WZHgg->ProjectionX("WZHggMET30",bin30,-1,"eo");
  TH1D* TTHggMET30 = (TH1D*)h_TTHgg->ProjectionX("TTHggMET30",bin30,-1,"eo");
  TH1D* VBFHggMET30 = (TH1D*)h_VBFHgg->ProjectionX("VBFHggMET30",bin30,-1,"eo");
  //TH1D* AAW130ProjX = (TH1D*)aaW130->ProjectionX("AAW130ProjX",0,-1,"eo");
  //TH1D* AAW275ProjX = (TH1D*)aaW275->ProjectionX("AAW275ProjX",0,-1,"eo");

  TH1D* ggHggNoMET_JetReq = (TH1D*)h_ggHgg_JetReq->ProjectionX("ggHggNoMET_JetReq",0,-1,"eo");
  TH1D* WZHggNoMET_JetReq = (TH1D*)h_WZHgg_JetReq->ProjectionX("WZHggNoMET_JetReq",0,-1,"eo");
  TH1D* TTHggNoMET_JetReq = (TH1D*)h_TTHgg_JetReq->ProjectionX("TTHggNoMET_JetReq",0,-1,"eo");
  TH1D* VBFHggNoMET_JetReq = (TH1D*)h_VBFHgg_JetReq->ProjectionX("VBFHggNoMET_JetReq",0,-1,"eo");

  TH1D* ggHggMet = (TH1D*)h_ggHgg->ProjectionY("ggHggMet",0,-1,"eo");
  TH1D* WZHggMet = (TH1D*)h_WZHgg->ProjectionY("WZHggMet",0,-1,"eo");
  TH1D* TTHggMet = (TH1D*)h_TTHgg->ProjectionY("TTHggMet",0,-1,"eo");
  TH1D* VBFHggMet = (TH1D*)h_VBFHgg->ProjectionY("VBFHggMet",0,-1,"eo");
  //TH1D* AAW130ProjY = (TH1D*)aaW130->ProjectionY("AAW130ProjY",0,-1,"eo");
  //TH1D* AAW275ProjY = (TH1D*)aaW275->ProjectionY("AAW275ProjY",0,-1,"eo");
  TH1D* ggHggMet_JetReq = (TH1D*)h_ggHgg_JetReq->ProjectionY("ggHggMet_JetReq",0,-1,"eo");
  TH1D* WZHggMet_JetReq = (TH1D*)h_WZHgg_JetReq->ProjectionY("WZHggMet_JetReq",0,-1,"eo");
  TH1D* TTHggMet_JetReq = (TH1D*)h_TTHgg_JetReq->ProjectionY("TTHggMet_JetReq",0,-1,"eo");
  TH1D* VBFHggMet_JetReq = (TH1D*)h_VBFHgg_JetReq->ProjectionY("VBFHggMet_JetReq",0,-1,"eo");
  /*
    TH1F* ggHggNoMET = (TH1F*)f_ggHgg->Get("ggInvarMassMVAcorrVertexCorr");
    TH1F* WZHggNoMET = (TH1F*)f_WZHgg->Get("ggInvarMassMVAcorrVertexCorr");
    TH1F* TTHggNoMET = (TH1F*)f_TTHgg->Get("ggInvarMassMVAcorrVertexCorr");
    TH1F* VBFHggNoMET = (TH1F*)f_VBFHgg->Get("ggInvarMassMVAcorrVertexCorr");
    TH1F* ggHggMET30 = (TH1F*)f_ggHgg->Get("ggInvarMassMET30MVAcorr");
    TH1F* WZHggMET30 = (TH1F*)f_WZHgg->Get("ggInvarMassMET30MVAcorr");
    TH1F* TTHggMET30 = (TH1F*)f_TTHgg->Get("ggInvarMassMET30MVAcorr");
    TH1F* VBFHggMET30 = (TH1F*)f_VBFHgg->Get("ggInvarMassMET30MVAcorr");

    //fix these to have mvacorrvertexcorr in jetreq - fixed, now all invarmass are mvacorrvertexcorr
    TH1F* ggHggNoMET_JetReq = (TH1F*)f_ggHgg->Get("ggInvarMass_JetReq");
    TH1F* WZHggNoMET_JetReq = (TH1F*)f_WZHgg->Get("ggInvarMass_JetReq");
    TH1F* TTHggNoMET_JetReq = (TH1F*)f_TTHgg->Get("ggInvarMass_JetReq");
    TH1F* VBFHggNoMET_JetReq = (TH1F*)f_VBFHgg->Get("ggInvarMass_JetReq");

    TH1F* ggHggMet = (TH1F*)f_ggHgg->Get("ggMet");
    TH1F* WZHggMet = (TH1F*)f_WZHgg->Get("ggMet");
    TH1F* TTHggMet = (TH1F*)f_TTHgg->Get("ggMet");
    TH1F* VBFHggMet = (TH1F*)f_VBFHgg->Get("ggMet");
    TH1F* ggHggMet_JetReq = (TH1F*)f_ggHgg->Get("ggMet_JetReq");
    TH1F* WZHggMet_JetReq = (TH1F*)f_WZHgg->Get("ggMet_JetReq");
    TH1F* TTHggMet_JetReq = (TH1F*)f_TTHgg->Get("ggMet_JetReq");
    TH1F* VBFHggMet_JetReq = (TH1F*)f_VBFHgg->Get("ggMet_JetReq");
  */
  //don't need to call sumw2 since it was called in th2
  /*
  ggHggNoMET->Sumw2();WZHggNoMET->Sumw2();TTHggNoMET->Sumw2();VBFHggNoMET->Sumw2();
  ggHggMet->Sumw2();WZHggMet->Sumw2();TTHggMet->Sumw2();VBFHggMet->Sumw2();
  ggHggMet_JetReq->Sumw2();WZHggMet_JetReq->Sumw2();TTHggMet_JetReq->Sumw2();VBFHggMet_JetReq->Sumw2();
  ggHggMET30->Sumw2();WZHggMET30->Sumw2();TTHggMET30->Sumw2();VBFHggMET30->Sumw2();
  ggHggNoMET_JetReq->Sumw2();WZHggNoMET_JetReq->Sumw2();TTHggNoMET_JetReq->Sumw2();VBFHggNoMET_JetReq->Sumw2();*/
  //8TeV SM higgs production cross sections: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt8TeV2012ICHEP
  ggHggNoMET->Scale((PhoEffScale2*L_int*2.29e-03*19.22)/99989.);//125GeV=/96290);  
  VBFHggNoMET->Scale((PhoEffScale2*L_int*2.29e-03*1.544)/95677.);//125GeV=/99885); 
  TTHggNoMET->Scale((PhoEffScale2*L_int*2.29e-03*.1271)/100048.);//125GeV=/100224);
  WZHggNoMET->Scale((PhoEffScale2*L_int*2.29e-03*(.6782/**(.3257+.014)*/+.3843/**.2*/))/100320);//125 and 126 GeV have same # events
  ggHggNoMET->SetLineColor(kGreen);ggHggNoMET->SetMarkerColor(kGreen);
  WZHggNoMET->SetLineColor(kCyan);WZHggNoMET->SetMarkerColor(kCyan);
  TTHggNoMET->SetLineColor(31);TTHggNoMET->SetMarkerColor(31);
  VBFHggNoMET->SetLineColor(kRed+3);VBFHggNoMET->SetMarkerColor(kRed+3);
  //ggHggNoMET->SetFillStyle(0);WZHggNoMET->SetFillStyle(0);TTHggNoMET->SetFillStyle(0);VBFHggNoMET->SetFillStyle(0);
  ggHggNoMET->SetFillColor(kGreen);WZHggNoMET->SetFillColor(kCyan);TTHggNoMET->SetFillColor(31);VBFHggNoMET->SetFillColor(kRed+3);
  ggHggNoMET->SetLineWidth(2);WZHggNoMET->SetLineWidth(2);VBFHggNoMET->SetLineWidth(2);TTHggNoMET->SetLineWidth(2);


  ggHggNoMET_uncorr->Scale((PhoEffScale2*L_int*2.29e-03*19.22)/99989.);//125GeV=/96290);  
  VBFHggNoMET_uncorr->Scale((PhoEffScale2*L_int*2.29e-03*1.544)/95677.);//125GeV=/99885); 
  TTHggNoMET_uncorr->Scale((PhoEffScale2*L_int*2.29e-03*.1271)/100048.);//125GeV=/100224);
  WZHggNoMET_uncorr->Scale((PhoEffScale2*L_int*2.29e-03*(.6782/**(.3257+.014)*/+.3843/**.2*/))/100320);//125 and 126 GeV have same # events
  ggHggNoMET_uncorr->SetLineColor(kGreen);ggHggNoMET_uncorr->SetMarkerColor(kGreen);
  WZHggNoMET_uncorr->SetLineColor(kCyan);WZHggNoMET_uncorr->SetMarkerColor(kCyan);
  TTHggNoMET_uncorr->SetLineColor(31);TTHggNoMET_uncorr->SetMarkerColor(31);
  VBFHggNoMET_uncorr->SetLineColor(kRed+3);VBFHggNoMET_uncorr->SetMarkerColor(kRed+3);
  ggHggNoMET_uncorr->SetFillColor(kGreen);WZHggNoMET_uncorr->SetFillColor(kCyan);TTHggNoMET_uncorr->SetFillColor(31);VBFHggNoMET_uncorr->SetFillColor(kRed+3);
  ggHggNoMET_uncorr->SetLineWidth(2);WZHggNoMET_uncorr->SetLineWidth(2);VBFHggNoMET_uncorr->SetLineWidth(2);TTHggNoMET_uncorr->SetLineWidth(2);
 
  ggHggMet->Scale((PhoEffScale2*L_int*2.29e-03*19.22)/99989.);//125GeV=/96290);  
  VBFHggMet->Scale((PhoEffScale2*L_int*2.29e-03*1.544)/95677.);//125GeV=/99885); 
  TTHggMet->Scale((PhoEffScale2*L_int*2.29e-03*.1271)/100048.);//125GeV=/100224);
  WZHggMet->Scale((PhoEffScale2*L_int*2.29e-03*(.6782/**(.3257+.014)*/+.3843/**.2*/))/100320);//125 and 126 GeV have same # events
  ggHggMet->SetLineColor(kGreen);ggHggMet->SetMarkerColor(kGreen);ggHggMet->SetFillColor(kGreen);
  WZHggMet->SetLineColor(kCyan);WZHggMet->SetMarkerColor(kCyan);WZHggMet->SetFillColor(kCyan);
  TTHggMet->SetLineColor(31);TTHggMet->SetMarkerColor(31);TTHggMet->SetFillColor(31);
  VBFHggMet->SetLineColor(kRed+3);VBFHggMet->SetMarkerColor(kRed+3);VBFHggMet->SetFillColor(kRed+3);
  //ggHggMet->SetFillStyle(0);WZHggMet->SetFillStyle(0);TTHggMet->SetFillStyle(0);VBFHggMet->SetFillStyle(0);
  ggHggMet->SetLineWidth(2);WZHggMet->SetLineWidth(2);VBFHggMet->SetLineWidth(2);TTHggMet->SetLineWidth(2);

  ggHggMET30->Scale((PhoEffScale2*L_int*2.29e-03*19.22)/99989.);//125GeV=/96290);  
  VBFHggMET30->Scale((PhoEffScale2*L_int*2.29e-03*1.544)/95677.);//125GeV=/99885); 
  TTHggMET30->Scale((PhoEffScale2*L_int*2.29e-03*.1271)/100048.);//125GeV=/100224);
  WZHggMET30->Scale((PhoEffScale2*L_int*2.29e-03*(.6782/**(.3257+.014)*/+.3843/**.2*/))/100320);//125 and 126 GeV have same # events
  ggHggMET30->SetLineColor(kGreen);ggHggMET30->SetMarkerColor(kGreen);
  WZHggMET30->SetLineColor(kCyan);WZHggMET30->SetMarkerColor(kCyan);
  TTHggMET30->SetLineColor(31);TTHggMET30->SetMarkerColor(31);
  VBFHggMET30->SetLineColor(kRed+3);VBFHggMET30->SetMarkerColor(kRed+3);
  ggHggMET30->SetFillStyle(0);WZHggMET30->SetFillStyle(0);TTHggMET30->SetFillStyle(0);VBFHggMET30->SetFillStyle(0);
  ggHggMET30->SetLineWidth(2);WZHggMET30->SetLineWidth(2);VBFHggMET30->SetLineWidth(2);TTHggMET30->SetLineWidth(2);
  ggHggNoMET_JetReq->Scale((PhoEffScale2*L_int*2.29e-03*19.22)/99989.);//125GeV=/96290);  
  VBFHggNoMET_JetReq->Scale((PhoEffScale2*L_int*2.29e-03*1.544)/95677.);//125GeV=/99885); 
  TTHggNoMET_JetReq->Scale((PhoEffScale2*L_int*2.29e-03*.1271)/100048.);//125GeV=/100224);
  WZHggNoMET_JetReq->Scale((PhoEffScale2*L_int*2.29e-03*(.6782/**(.3257+.014)*/+.3843/**.2*/))/100320);//125 and 126 GeV have same # events
  ggHggNoMET_JetReq->SetLineColor(kGreen);ggHggNoMET_JetReq->SetMarkerColor(kGreen);
  WZHggNoMET_JetReq->SetLineColor(kCyan);WZHggNoMET_JetReq->SetMarkerColor(kCyan);
  TTHggNoMET_JetReq->SetLineColor(31);TTHggNoMET_JetReq->SetMarkerColor(31);
  VBFHggNoMET_JetReq->SetLineColor(kRed+3);VBFHggNoMET_JetReq->SetMarkerColor(kRed+3);
  ggHggNoMET_JetReq->SetFillStyle(0);WZHggNoMET_JetReq->SetFillStyle(0);TTHggNoMET_JetReq->SetFillStyle(0);VBFHggNoMET_JetReq->SetFillStyle(0);
  ggHggNoMET_JetReq->SetLineWidth(2);WZHggNoMET_JetReq->SetLineWidth(2);VBFHggNoMET_JetReq->SetLineWidth(2);TTHggNoMET_JetReq->SetLineWidth(2);

  ggHggNoMET->GetXaxis()->SetRangeUser(110,139.9);
  ggHggNoMET->Draw("histo");
  WZHggNoMET->Draw("histoSAMES");
  TTHggNoMET->Draw("histoSAMES");
  VBFHggNoMET->Draw("histoSAMES");
  c1->Print("Plots/Higgs/ggInvarMassSMhiggsIndividual.png");
  ggHggNoMET->GetXaxis()->SetRangeUser(0,1000);

  ggHggNoMET_uncorr->GetXaxis()->SetRangeUser(110,139.9);
  ggHggNoMET_uncorr->Draw("histo");
  WZHggNoMET_uncorr->Draw("histoSAMES");
  TTHggNoMET_uncorr->Draw("histoSAMES");
  VBFHggNoMET_uncorr->Draw("histoSAMES");
  c1->Print("Plots/Higgs/ggInvarMassSMhiggsIndividual_uncorrected.png");
  ggHggNoMET_uncorr->GetXaxis()->SetRangeUser(0,1000);


  TH1F* ggHggMetRebin=(TH1F*)ggHggMet->Rebin(NmetBins,"ggHggMetRebin",xbins);
  TH1F* TTHggMetRebin=(TH1F*)TTHggMet->Rebin(NmetBins,"TTHggMetRebin",xbins);
  TH1F* WZHggMetRebin=(TH1F*)WZHggMet->Rebin(NmetBins,"WZHggMetRebin",xbins);
  TH1F* VBFHggMetRebin=(TH1F*)VBFHggMet->Rebin(NmetBins,"VBFHggMetRebin",xbins);
  //TH1F* AAW130ProjYRebin=(TH1F*)AAW130ProjY->Rebin(NmetBins,"AAW130ProjYRebin",xbins);
  //TH1F* AAW275ProjYRebin=(TH1F*)AAW275ProjY->Rebin(NmetBins,"AAW275ProjYRebin",xbins);

  AddOverflowToLastBin(ggHggMetRebin);
  AddOverflowToLastBin(TTHggMetRebin);
  AddOverflowToLastBin(WZHggMetRebin);
  AddOverflowToLastBin(VBFHggMetRebin);
  //AddOverflowToLastBin(AAW130ProjYRebin);
  //AddOverflowToLastBin(AAW275ProjYRebin);
  DivideByBinWidth(ggHggMetRebin);
  DivideByBinWidth(TTHggMetRebin);
  DivideByBinWidth(WZHggMetRebin);
  DivideByBinWidth(VBFHggMetRebin);
  //DivideByBinWidth(AAW130ProjYRebin);
  //DivideByBinWidth(AAW275ProjYRebin);

  TH1F* ggHggMetRebin30=(TH1F*)ggHggMetRebin->Clone();
  TH1F* TTHggMetRebin30=(TH1F*)TTHggMetRebin->Clone();
  TH1F* VBFHggMetRebin30=(TH1F*)VBFHggMetRebin->Clone();
  TH1F* WZHggMetRebin30=(TH1F*)WZHggMetRebin->Clone();
  int bin30higgs=ggHggMetRebin30->GetXaxis()->FindBin(30);
  for(int i=0;i<bin30higgs;i++){
    ggHggMetRebin30->SetBinContent(i,0);ggHggMetRebin30->SetBinError(i,0);
    TTHggMetRebin30->SetBinContent(i,0);TTHggMetRebin30->SetBinError(i,0);
    VBFHggMetRebin30->SetBinContent(i,0);VBFHggMetRebin30->SetBinError(i,0);
    WZHggMetRebin30->SetBinContent(i,0);WZHggMetRebin30->SetBinError(i,0);
  }
  c1->SetLogy(1);
  THStack *HiggsMetStackRaw = new THStack("HiggsMetStackRaw","");
  HiggsMetStackRaw->Add(TTHggMet);HiggsMetStackRaw->Add(WZHggMet);HiggsMetStackRaw->Add(VBFHggMet);HiggsMetStackRaw->Add(ggHggMet);
  HiggsMetStackRaw->Draw("histo");
  HiggsMetStackRaw->GetXaxis()->SetRangeUser(0,249);HiggsMetStackRaw->SetMaximum(14);HiggsMetStackRaw->SetMinimum(3e-3);
  c1->Print("Plots/Higgs/ggMET_SMhiggsStack.png");

  THStack *HiggsMetStack = new THStack("HiggsMetStack","");
  TTHggMetRebin->SetLineColor(kBlack);VBFHggMetRebin->SetLineColor(kBlack);ggHggMetRebin->SetLineColor(kBlack);WZHggMetRebin->SetLineColor(kBlack);
  TTHggMetRebin->SetLineWidth(1);VBFHggMetRebin->SetLineWidth(1);ggHggMetRebin->SetLineWidth(1);WZHggMetRebin->SetLineWidth(1);
  HiggsMetStack->Add(TTHggMetRebin);HiggsMetStack->Add(WZHggMetRebin);HiggsMetStack->Add(VBFHggMetRebin);HiggsMetStack->Add(ggHggMetRebin);
  HiggsMetStack->Draw("histo");
  HiggsMetStack->GetXaxis()->SetRangeUser(0,249);HiggsMetStack->SetMaximum(14);HiggsMetStack->SetMinimum(3e-3);
  HiggsMetStack->GetYaxis()->SetTitle("Events / GeV");HiggsMetStack->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  c1->Print("Plots/Higgs/ggMET_SMhiggsStackRebin.png");
  TH1F* SMHiggsMet=(TH1F*)TTHggMetRebin->Clone();SMHiggsMet->Add(VBFHggMetRebin);SMHiggsMet->Add(WZHggMetRebin);SMHiggsMet->Add(ggHggMetRebin);

  c1->SetLogy(0);
  TH1F* SMHiggs = (TH1F*)ggHggNoMET->Clone();
  SMHiggs->Add(WZHggNoMET);
  SMHiggs->Add(VBFHggNoMET);
  SMHiggs->Add(TTHggNoMET);SMHiggs->SetLineColor(kBlack);SMHiggs->SetMarkerColor(kBlack);
  RooRealVar xSM("xSM","m_{#gamma#gamma}",110,140,"GeV");xSM.setRange("sigSMrange",120,131);
  RooDataHist dataSM("dataSM","data SM higgs",xSM,SMHiggs);
  //Gaussian for peak
  RooRealVar meanSM("meanSM","meanSM",125.883/*126,124,128*/);
  RooRealVar sigmaSM("sigmaSM","sigmaSM",1./*,0.4,1.5*/);
  RooRealVar sigmaSM2("sigmaSM2","sigmaSM2",0.,6);
  //Voigtian (Breit-Vigner x Gaussian) for peak
  RooRealVar widthSM("widthSM","widthSM",2.2/*,1.,3.*/);
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
  SMPdf.plotOn(xframeSM,LineColor(kBlue),VisualizeError(*rSM,2,kTRUE),FillColor(kRed));
  SMPdf.plotOn(xframeSM,LineColor(kBlue),VisualizeError(*rSM,1,kTRUE),FillColor(kRed+2));
  SMPdf.plotOn(xframeSM,LineColor(kBlue));
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
  TLegend *legHiggs = new TLegend(.2,.6,.52,.8);
  legHiggs->AddEntry(ggHggNoMET,"GluGlu->h->#gamma#gamma","l");
  legHiggs->AddEntry(VBFHggNoMET,"VBFh->#gamma#gamma","l");
  legHiggs->AddEntry(WZHggNoMET,"W/Zh->#gamma#gamma","l");
  legHiggs->AddEntry(TTHggNoMET,"TTh->#gamma#gamma","l");
  //legHiggs->AddEntry(SMHiggs,"All h->#gamma#gamma","lp");
  legHiggs->SetFillColor(kWhite);
  legHiggs->SetFillStyle(0);
  legHiggs->SetBorderSize(0);
  TLegend *legHiggs2 = new TLegend(.53,.62,.85,.82);
  legHiggs2->AddEntry(ggHggNoMET,"GluGlu->h->#gamma#gamma","l");
  legHiggs2->AddEntry(VBFHggNoMET,"VBFh->#gamma#gamma","l");
  legHiggs2->AddEntry(WZHggNoMET,"W/Zh->#gamma#gamma","l");
  legHiggs2->AddEntry(TTHggNoMET,"TTh->#gamma#gamma","l");
  //legHiggs2->AddEntry(SMHiggs2,"All h->#gamma#gamma","lp");
  legHiggs2->SetFillColor(kWhite);
  legHiggs2->SetFillStyle(0);
  legHiggs2->SetBorderSize(0);
  legHiggs->Draw();
  c1->Print("Plots/Higgs/ggInvarMassSMhiggsFit.png");

  c1->SetLogy(0);
  TH1F* SMHiggs_uncorr = (TH1F*)ggHggNoMET_uncorr->Clone();
  SMHiggs_uncorr->Add(WZHggNoMET_uncorr);
  SMHiggs_uncorr->Add(VBFHggNoMET_uncorr);
  SMHiggs_uncorr->Add(TTHggNoMET_uncorr);
  RooRealVar xSM_uncorr("xSM_uncorr","m_{#gamma#gamma}",110,140,"GeV");xSM_uncorr.setRange("sigSMrange",120,131);
  RooDataHist dataSM_uncorr("dataSM_uncorr","data SM higgs",xSM_uncorr,SMHiggs_uncorr);
  //Gaussian for peak
  RooRealVar meanSM_uncorr("meanSM_uncorr","meanSM_uncorr",125.677/*126.,124.,128.*/);
  RooRealVar sigmaSM_uncorr("sigmaSM_uncorr","sigmaSM_uncorr",.67/*.7,0.4,1.*/);
  //RooRealVar sigmaSM_uncorr2("sigmaSM_uncorr2","sigmaSM_uncorr2",0.,6);
  //Voigtian (Breit-Vigner x Gaussian) for peak
  RooRealVar widthSM_uncorr("widthSM_uncorr","widthSM_uncorr",1.98/*2.,1.2,2.8*/);
  RooRealVar SMYield_uncorr("sig1 yield","sig1 yield",20000,0,400000);
  //RooRealVar SMYield2("sig2 yield","sig2 yield",20000,0,400000);
  //RooGaussian sigSM2("sigSM2","gaussSM2",xSM,meanSM,sigmaSM2);
  //RooGaussian sigSM ("sigSM", "gaussSM", xSM,meanSM,sigmaSM);
  RooVoigtian sigSM_uncorr ("sigSM_uncorr", "gaussSM_uncorr", xSM_uncorr,meanSM_uncorr,widthSM_uncorr,sigmaSM_uncorr);
  //RooAddPdf SMPdf("SMPdf","SMPdf",RooArgList(sigSM,sigSM2),RooArgList(SMYield,SMYield2));
  RooAddPdf SMPdf_uncorr("SMPdf_uncorr","SMPdf_uncorr",RooArgList(sigSM_uncorr),RooArgList(SMYield_uncorr));
  RooFitResult *rSM_uncorr = SMPdf_uncorr.fitTo(dataSM_uncorr,Extended(kTRUE),Save(),Range("sigSMrange"));
  RooPlot *xframeSM_uncorr = xSM_uncorr.frame(Title("Breit-Wigner#otimesGaussian signal, gg SM higgs"));
  SMPdf_uncorr.paramOn(xframeSM_uncorr, Format("NE",AutoPrecision(2)),Layout(0.6,0.995,0.9) );
  dataSM_uncorr.plotOn(xframeSM_uncorr,LineColor(kBlack));
  SMPdf_uncorr.plotOn(xframeSM_uncorr,LineColor(kBlue));
  SMPdf_uncorr.plotOn(xframeSM_uncorr,LineColor(kBlue),VisualizeError(*rSM_uncorr,2,kTRUE),FillColor(kRed));
  SMPdf_uncorr.plotOn(xframeSM_uncorr,LineColor(kBlue),VisualizeError(*rSM_uncorr,1,kTRUE),FillColor(kRed+2));
  SMPdf_uncorr.plotOn(xframeSM_uncorr,LineColor(kBlue));
  //SMPdf.plotOn(xframeSM,Components(sigSM),LineColor(kRed),LineStyle(kDashed));
  //SMPdf.plotOn(xframeSM,Components(sigSM2),LineColor(kGreen),LineStyle(kDashed));
  dataSM_uncorr.plotOn(xframeSM_uncorr,LineColor(kBlack));
  xframeSM_uncorr->Draw();
  THStack* StackMet_uncorr = new THStack("StackMet_uncorr","");
  StackMet_uncorr->Add(TTHggNoMET_uncorr);//->Draw("SAME");
  StackMet_uncorr->Add(WZHggNoMET_uncorr);//->Draw("SAME"); 
  StackMet_uncorr->Add(VBFHggNoMET_uncorr);//->Draw("SAME");
  StackMet_uncorr->Add(ggHggNoMET_uncorr);//->Draw("SAME");
  StackMet_uncorr->Draw("histoSAMES");
  xframeSM_uncorr->Draw("SAMES");
  legHiggs->Draw();
  c1->Print("Plots/Higgs/ggInvarMassSMhiggsFit_uncorrected.png");
 
  StackMet->Draw("histo");
  StackMet->GetXaxis()->SetRangeUser(110,139.9);StackMet->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV]");StackMet->GetYaxis()->SetTitle("Events");
  legHiggs->Draw();
  c1->Print("Plots/Higgs/ggInvarMassSMhiggs.png");
  c1->Print("Plots/Higgs/ggInvarMassSMhiggs.pdf");

  StackMet_uncorr->Draw("histo");
  StackMet_uncorr->GetXaxis()->SetRangeUser(110,139.9);
  legHiggs->Draw();
  c1->Print("Plots/Higgs/ggInvarMassSMhiggs_uncorrected.png");


  ggHggNoMET->Rebin(2);WZHggNoMET->Rebin(2);TTHggNoMET->Rebin(2);VBFHggNoMET->Rebin(2);
  ggHggMET30->Rebin(2);WZHggMET30->Rebin(2);TTHggMET30->Rebin(2);VBFHggMET30->Rebin(2);
  //ggHggNoMET_JetReq->Rebin(2);WZHggNoMET_JetReq->Rebin(2);TTHggNoMET_JetReq->Rebin(2);VBFHggNoMET_JetReq->Rebin(2);


  //data, no met cut

  RooRealVar xMet("xMet","m_{#gamma#gamma}",95,200,"GeV");
  xMet.setRange("full",sbLoLo,sbHiHi);
  RooDataHist dataMet("dataMet","dataset",xMet,h_ggInvarMassMET);
  //Gaussian for peak
  RooRealVar gmeanMet("gmeanMet","gmean",125.3,122,128);
  RooRealVar gsigmaMet("gsigmaMet","gsigma",1.6/*,0.1,5*/);
  //RooGaussian sigMet("sigMet","gauss",xMet,gmeanMet,gsigmaMet);
  //Crystal Ball for signal
  RooRealVar meanMet("cb_mean", "mean" , 125.3, 122, 128.) ;
  RooRealVar sigmaMet("cb_sigma", "sigma",1.6/*.6, 0., 3.*/);
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
  RooRealVar Pol1("Pol1","Pol1",0,-.09,.09);
  RooRealVar Pol2("Pol2","Pol2",0,-.0001,.0001);
  RooRealVar Pol3("Pol3","Pol3",0,-.00001,.00001);
  RooRealVar Pol4("Pol4","Pol4",0,-.000001,.000001);
  RooRealVar Pol5("Pol5","Pol5",0,-.00000001,.00000001);
  //RooBernstein Bern("Bern","4th order Bernstein Polynomial",xMet,RooArgList(Bern1,Bern2,Bern3,Bern4,Bern5));
  RooPolynomial Bern("Bern","5th Order Polynomial",xMet,RooArgList(Pol1,Pol2,Pol3,Pol4,Pol5));
  RooRealVar BernYield("bkgd yield","bkgd yield",20000,0,400000);
  RooAddPdf MetPdf("MetPdf","MetPdf",RooArgList(Bern,sigMet),RooArgList(BernYield,sigMetYield));
  RooFitResult *rNoMetCut = MetPdf.fitTo(dataMet,Extended(kTRUE),Save());
  RooPlot *xframeMet = xMet.frame(Title("Voigtian Signal, 5th Order Polynomial Background, gg with no MET cut"));
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
  xframeMet->SetAxisRange(0,5200,"Y");
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
  c1->Print("Plots/Higgs/ggInvarMassWithNoMETcut.png");
  
  RooRealVar xMet_bg("xMet_bg","m_{#gamma#gamma}",95,200,"GeV");
  xMet_bg.setRange("sb_lo",sbLoLo,sbLoHi);xMet_bg.setRange("sbFit_lo",sbFitLoLo,sbLoHi);xMet_bg.setRange("sb_hi",sbHiLo,sbHiHi);xMet.setRange("full",sbLoLo,sbHiHi);xMet.setRange("sig",sigLo,sigHi);
  RooDataHist dataMet_bg("dataMet_bg","dataset_bg",xMet_bg,h_ggInvarMassMET);
  RooRealVar Bern1_bg("Bern1_bg","Berstein 1",14.5,8.,22.);
  RooRealVar Bern2_bg("Bern2_bg","Berstein 2",5.,1.,9.);
  RooRealVar Bern3_bg("Bern3_bg","Berstein 3",4.,0.,8.);
  RooRealVar Bern4_bg("Bern4_bg","Berstein 4",2.,0.,4.);
  RooRealVar Bern5_bg("Bern5_bg","Berstein 5",5.,0.,10.);
  RooRealVar Pol1_bg("Pol1_bg","Pol1_bg",0,-.09,.09);
  RooRealVar Pol2_bg("Pol2_bg","Pol2_bg",0,-.0001,.0001);
  RooRealVar Pol3_bg("Pol3_bg","Pol3_bg",0,-.00001,.00001);
  RooRealVar Pol4_bg("Pol4_bg","Pol4_bg",0,-.000001,.000001);
  RooRealVar Pol5_bg("Pol5_bg","Pol5_bg",0,-.00000001,.00000001);
  //RooBernstein Bern_bg("Bern","4th order Bernstein Polynomial",xMet_bg,RooArgList(Bern1_bg,Bern2_bg,Bern3_bg,Bern4_bg,Bern5_bg));
  RooPolynomial Bern_bg("Bern","4th order Bernstein Polynomial",xMet_bg,RooArgList(Pol1_bg,Pol2_bg,Pol3_bg,Pol4_bg,Pol5_bg));
  RooRealVar BernYield_bg("bkgd yield_bg","bkgd yield_bg",70000,0,170000);
  RooAddPdf MetPdf_bg("MetPdf_bg","MetPdf_bg",RooArgList(Bern_bg),RooArgList(BernYield_bg));
  //RooFitResult *rNoMetCut_bg = MetPdf_bg.chi2FitTo(dataMet_bg/*,Extended(kTRUE)*/,Save(),Range("sb_lo,sb_hi"));
  RooFitResult *rNoMetCut_bg = MetPdf_bg.fitTo(dataMet_bg,Extended(kTRUE),Save(),Range("sb_lo,sb_hi"));
  // Print fit results 
  cout << "result of fit on all data " << endl ;
  //rNoMetCut_bg_full->Print() ;  
  cout << "result of fit in in background region" << endl ;
  //rNoMetCut_bg->Print() ;
  cout << "result of fit in background region" << endl ;
  RooPlot *xframeMet_bg = xMet_bg.frame(Title("5th Order Polynomial Background, gg with no MET cut"));
  MetPdf_bg.paramOn(xframeMet_bg, Format("NE",AutoPrecision(2)),Layout(0.6,0.995,0.9) );
  dataMet_bg.plotOn(xframeMet_bg,LineColor(kBlack),MarkerStyle(20),MarkerSize(0.3));
  //MetPdf_bg.plotOn(xframeMet_bg,Components(Bern_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg,3,kTRUE),FillColor(kViolet));
  MetPdf_bg.plotOn(xframeMet_bg,Components(Bern_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg,2,kTRUE),FillColor(kGreen));
  MetPdf_bg.plotOn(xframeMet_bg,Components(Bern_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg,1,kTRUE),FillColor(kOrange));
  //MetPdf_bg.plotOn(xframeMet_bg,Components(Bern_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg_full,3,kFALSE),FillColor(kViolet));
  //MetPdf_bg.plotOn(xframeMet_bg,Components(Bern_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg_full,2,kFALSE),FillColor(kGreen));
  //MetPdf_bg.plotOn(xframeMet_bg,Components(Bern_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg_full,1,kFALSE),FillColor(kOrange));
  MetPdf_bg.plotOn(xframeMet_bg,Components(Bern_bg),LineColor(kRed),LineStyle(kDashed),Range(102,180));
  MetPdf_bg.plotOn(xframeMet_bg,Components(Bern_bg),LineColor(kRed));
  dataMet_bg.plotOn(xframeMet_bg,LineColor(kBlack),MarkerStyle(20),MarkerSize(0.3));
  //xframeMet_bg->SetAxisRange(106,149,"X");
  //xframeMet_bg->SetAxisRange(0,1900,"Y");
  //xframeMet_bg->SetMaximum(2000);
  xframeMet_bg->Draw();
  //line125->DrawLine(125.3,0,125.3,xframeMet_bg->GetMaximum());
  StackMet->Draw("histoSAME");
  xframeMet_bg->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/Higgs/ggInvarMassWithNoMETcut_bgOnly.png");
  

  xMet_bg.setRange("SigFull",sigLo,sigHi);/*
					    xMet_bg.setRange("Sig",120,121);
					    xMet_bg.setRange("Sig1",121,122);
					    xMet_bg.setRange("Sig2",122,123);
					    xMet_bg.setRange("Sig3",123,124);
					    xMet_bg.setRange("Sig4",124,125);
					    xMet_bg.setRange("Sig5",125,126);
					    xMet_bg.setRange("Sig6",126,127);
					    xMet_bg.setRange("Sig7",127,128);
					    xMet_bg.setRange("Sig8",128,129);
					    xMet_bg.setRange("Sig9",129,130);
					    //xMet_bg.setRange("xMetLowSB",sbLoLo,110);
					    //xMet_bg.setRange("xMetHighSB",135,160);
					    RooAbsReal* ig = MetPdf_bg.createIntegral(xMet_bg,NormSet(xMet_bg),Range("Sig"));
					    RooAbsReal* ig1 = MetPdf_bg.createIntegral(xMet_bg,NormSet(xMet_bg),Range("Sig1"));
					    RooAbsReal* ig2 = MetPdf_bg.createIntegral(xMet_bg,NormSet(xMet_bg),Range("Sig2"));
					    RooAbsReal* ig3 = MetPdf_bg.createIntegral(xMet_bg,NormSet(xMet_bg),Range("Sig3"));
					    RooAbsReal* ig4 = MetPdf_bg.createIntegral(xMet_bg,NormSet(xMet_bg),Range("Sig4"));
					    RooAbsReal* ig5 = MetPdf_bg.createIntegral(xMet_bg,NormSet(xMet_bg),Range("Sig5"));
					    RooAbsReal* ig6 = MetPdf_bg.createIntegral(xMet_bg,NormSet(xMet_bg),Range("Sig6"));
					    RooAbsReal* ig7 = MetPdf_bg.createIntegral(xMet_bg,NormSet(xMet_bg),Range("Sig7"));
					    RooAbsReal* ig8 = MetPdf_bg.createIntegral(xMet_bg,NormSet(xMet_bg),Range("Sig8"));
					    RooAbsReal* ig9 = MetPdf_bg.createIntegral(xMet_bg,NormSet(xMet_bg),Range("Sig9"));*/
  RooAbsReal* igSig = MetPdf_bg.createIntegral(xMet_bg,NormSet(xMet_bg),Range("SigFull"));
  RooAbsReal* igSigNoNorm = MetPdf_bg.createIntegral(xMet_bg,Range("SigFull"));
  RooAbsReal* igFull = MetPdf_bg.createIntegral(xMet_bg);
  RooAbsReal* igFullMinusSig = MetPdf_bg.createIntegral(xMet_bg,Range("sb_lo,sb_hi"));
  RooAbsReal* igSBlo = MetPdf_bg.createIntegral(xMet_bg,NormSet(xMet_bg),Range("sb_lo"));
  RooAbsReal* igSBhi = MetPdf_bg.createIntegral(xMet_bg,NormSet(xMet_bg),Range("sb_hi"));

  float HiggSigNoNorm=igSig->getVal()*BernYield_bg.getVal();
  float HiggLowSBNoNorm=igSBlo->getVal()*BernYield_bg.getVal(),HiggHighSBNoNorm=igSBhi->getVal()*BernYield_bg.getVal();
  cout<<"HiggSigNoNorm   : "<<HiggSigNoNorm<<endl;
  cout<<"HiggLowSBNoNorm : "<<HiggLowSBNoNorm<<endl;
  cout<<"HiggHighSBNoNorm: "<<HiggHighSBNoNorm<<endl;
  vector<RooAbsReal*> integralInts(0);
  h_ggInvMassSubtract->Rebin(2);
  for(int i=1;i<h_ggInvMassSubtract->GetBinLowEdge(h_ggInvMassSubtract->GetNbinsX());i++){//does this only work if the bins are integers?
    TString name = "count";name+=i;
    xMet_bg.setRange(name,h_ggInvMassSubtract->GetBinLowEdge(i),h_ggInvMassSubtract->GetBinLowEdge(i+1));
    integralInts.push_back(MetPdf_bg.createIntegral(xMet_bg,NormSet(xMet_bg),Range(name)));
    cout<<"i=:"<<i<<endl;
  }

  cout << "gx = " << MetPdf_bg.getVal() << endl ;
  cout << "gxYield = " << BernYield_bg.getVal() << endl ;
  
  // Return value of gx normalized over x in range [-10,10]
  RooArgSet nset(xMet_bg) ;
  cout << "gx_Norm[x] = " << MetPdf_bg.getVal(&nset) << endl ;

  // Create object representing integral over gx
  // which is used to calculate  gx_Norm[x] == gx / gx_Int[x]
  RooAbsReal* igx = MetPdf_bg.createIntegral(xMet_bg) ;
  cout << "gx_Int[x] = " << igx->getVal() << endl ;
  
  xMet_bg.setRange("signal",sigLo,129) ;
  RooAbsReal* igx_sig = MetPdf_bg.createIntegral(xMet_bg,NormSet(xMet_bg),Range("signal")) ;
  cout << "gx_Int[x|signal]_Norm[x] = " << igx_sig->getVal() << endl ;
  cout << "gx_Int[x|signal]_Norm[x] * Yield = " << igx_sig->getVal()*BernYield_bg.getVal() << endl ;

  cout<<"igFull: "<<igFull->getVal()<<endl;

  h_ggInvMassSubtract->GetXaxis()->SetRangeUser(100,164.9);
  h_ggInvMassSubtract->Draw();
  c1->Print("Plots/Higgs/ggInvarMassWithNoMETcut_bgOnly_presubtract.png");
  int bin120=h_ggInvMassSubtract->FindBin(sigLo),bin130=h_ggInvMassSubtract->FindBin(sigHi)-1;
  cout<<"full int:"<<igSig->getVal()<<"  full int no norm:"<<igSigNoNorm->getVal()<<"  full int * yield:"<<igSig->getVal()*BernYield_bg.getVal()<<"  data integral: "<<h_ggInvMassSubtract->Integral(bin120,bin130)<<endl;
  for(int i=1;i<h_ggInvMassSubtract->GetNbinsX()-1;i++){
    float x=0,xe=0;
    x=h_ggInvMassSubtract->GetBinContent(i);xe=h_ggInvMassSubtract->GetBinError(i);
    h_ggInvMassSubtract->SetBinContent(i,x-integralInts[i-1]->getVal()*BernYield_bg.getVal());
    //x=h_ggInvMassSubtract->GetBinContent(i);
    h_ggInvMassSubtract->SetBinError(i,sqrt(xe*xe+sqrt(integralInts[i-1]->getVal()*BernYield_bg.getVal()*integralInts[i-1]->getVal()*BernYield_bg.getVal())));
  }
  //h_ggInvMassSubtract->Rebin(2);h_ggInvMassSubtract->Scale(0.5);
  h_ggInvMassSubtract->GetXaxis()->SetRangeUser(100,164.9);//h_ggInvMassSubtract->GetYaxis()->SetRangeUser(-60,200);
  h_ggInvMassSubtract->Draw();StackMet->Draw("histoSAMES");h_ggInvMassSubtract->Draw("SAMES");
  c1->Print("Plots/Higgs/ggInvarMassWithNoMETcut_bgOnly_subtract.png");

  //SMhiggs, MET>30 SM higgs

  TH1F* SM30Higgs = (TH1F*)ggHggMET30->Clone();
  SM30Higgs->Add(WZHggMET30);
  SM30Higgs->Add(VBFHggMET30);
  SM30Higgs->Add(TTHggMET30);
  RooRealVar xSM30("xSM30","m_{#gamma#gamma}",110,140,"GeV");xSM30.setRange("sigSM30range",115,135);
  RooDataHist dataSM30("dataSM30","data SM30 higgs",xSM30,SM30Higgs);
  //Gaussian for peak
  RooRealVar meanSM30("meanSM30","meanSM30",126/*.3,122,128*/);
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
  c1->Print("Plots/Higgs/ggInvarMassSMhiggsMET30cut.png");

  //data, MET>30

  RooRealVar xMet30("xMet30","m_{#gamma#gamma}",95,200,"GeV");
  xMet30.setRange("full",sbLoLo,sbHiHi);
  RooDataHist dataMet30("dataMet30","dataset",xMet30,h_ggInvarMassMET30);
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
  RooPolynomial Bern30("Bern30","5th Order Polynomial",xMet30,RooArgList(Pol301,Pol302,Pol303,Pol304,Pol305));
  RooRealVar Bern30Yield("bkgd yield","bkgd yield",20000,0,400000);
  RooAddPdf Met30Pdf("Met30Pdf","Met30Pdf",RooArgList(Bern30,sigMet30),RooArgList(Bern30Yield,sigMet30Yield));
  RooFitResult *rMet30Cut = Met30Pdf.fitTo(dataMet30,Extended(kTRUE),Save());
  RooPlot *xframeMet30 = xMet30.frame(Title("Voigtian Signal, 5th Order Polynomial Background, gg with MET>30 cut"));
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
  // xframeMet30->SetAxisRange(106,149,"X");
  //xframeMet30->SetAxisRange(0,1300,"Y");
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
  c1->Print("Plots/Higgs/ggInvarMassWithMET30cut.png");
  
  RooRealVar xMet30_bg("xMet30_bg","m_{#gamma#gamma}",95,200,"GeV");
  xMet30_bg.setRange("sb_30lo",sbLoLo,sbLoHi);xMet30_bg.setRange("sbFit_30lo",sbFitLoLo,sbLoHi);xMet30_bg.setRange("sb_30hi",sbHiLo,sbHiHi);xMet30_bg.setRange("full",sbLoLo,sbHiHi);xMet30_bg.setRange("sig",sigLo,sigHi);
  RooDataHist dataMet30_bg("dataMet30_bg","dataset_bg",xMet30_bg,h_ggInvarMassMET30);
  RooRealVar Bern301_bg("Bern301_bg","Berstein 1",14.5,8.,22.);
  RooRealVar Bern302_bg("Bern302_bg","Berstein 2",5.,1.,9.);
  RooRealVar Bern303_bg("Bern303_bg","Berstein 3",4.,0.,8.);
  RooRealVar Bern304_bg("Bern304_bg","Berstein 4",2.,0.,4.);
  RooRealVar Bern305_bg("Bern305_bg","Berstein 5",5.,0.,10.);
  RooRealVar Pol301_bg("Pol301_bg","Pol301_bg",0,-.1,.1);
  RooRealVar Pol302_bg("Pol302_bg","Pol302_bg",0,-.00001,.00001);
  RooRealVar Pol303_bg("Pol303_bg","Pol303_bg",0,-.000001,.000001);
  RooRealVar Pol304_bg("Pol304_bg","Pol304_bg",0,-.0000001,.0000001);
  RooRealVar Pol305_bg("Pol305_bg","Pol305_bg",0,-.0000001,.0000001);
  //RooBern30stein Bern30_bg("Bern30","4th order Bern30stein Polynomial",xMet30_bg,RooArgList(Bern301_bg,Bern302_bg,Bern303_bg,Bern304_bg,Bern305_bg));
  RooPolynomial Bern30_bg("Bern30","4th order Bern30stein Polynomial",xMet30_bg,RooArgList(Pol301_bg,Pol302_bg,Pol303_bg,Pol304_bg,Pol305_bg));
  RooRealVar Bern30Yield_bg("bkgd yield_bg","bkgd yield_bg",70000,0,170000);
  RooAddPdf Met30Pdf_bg("Met30Pdf_bg","Met30Pdf_bg",RooArgList(Bern30_bg),RooArgList(Bern30Yield_bg));
  RooFitResult *rMet30Cut_bg = Met30Pdf_bg.fitTo(dataMet30_bg/*,Extended(kTRUE)*/,Save(),Range("sb_30lo,sb_30hi"));
  // RooFitResult *rMet30Cut_bg_full = Met30Pdf_bg.fitTo(dataMet30,Extended(kTRUE),Save());
  // Print fit results 
  cout << "result of fit on all data " << endl ;
  //rMet30Cut_bg_full->Print() ;  
  cout << "result of fit in in background region" << endl ;
  //rMet30Cut_bg->Print() ;
  cout << "result of fit in background region" << endl ;
  RooPlot *xframeMet30_bg = xMet30_bg.frame(Title("5th Order Polynomial Background, gg with MET>30 cut"));
  Met30Pdf_bg.paramOn(xframeMet30_bg, Format("NE",AutoPrecision(2)),Layout(0.6,0.995,0.9) );
  dataMet30_bg.plotOn(xframeMet30_bg,LineColor(kBlack),MarkerStyle(20),MarkerSize(0.3));
  //Met30Pdf_bg.plotOn(xframeMet30_bg,Components(Bern30_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rMet30Cut_bg,3,kTRUE),FillColor(kViolet));
  Met30Pdf_bg.plotOn(xframeMet30_bg,Components(Bern30_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rMet30Cut_bg,2,kTRUE),FillColor(kGreen));
  Met30Pdf_bg.plotOn(xframeMet30_bg,Components(Bern30_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rMet30Cut_bg,1,kTRUE),FillColor(kOrange));
  //Met30Pdf_bg.plotOn(xframeMet30_bg,Components(Bern30),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rMet30Cut_bg_full,3,kFALSE),FillColor(kViolet));
  //Met30Pdf_bg.plotOn(xframeMet30_bg,Components(Bern30),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rMet30Cut_bg_full,2,kFALSE),FillColor(kGreen));
  //Met30Pdf_bg.plotOn(xframeMet30_bg,Components(Bern30),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rMet30Cut_bg_full,1,kFALSE),FillColor(kOrange));
  Met30Pdf_bg.plotOn(xframeMet30_bg,Components(Bern30_bg),LineColor(kRed),LineStyle(kDashed),Range(95,180));
  Met30Pdf_bg.plotOn(xframeMet30_bg,Components(Bern30_bg),LineColor(kRed));
  dataMet30_bg.plotOn(xframeMet30_bg,LineColor(kBlack),MarkerStyle(20),MarkerSize(0.3));
  //xframeMet30_bg->SetAxisRange(106,149,"X");
  //xframeMet30_bg->SetAxisRange(0,1300,"Y");
  //xframeMet30_bg->SetMaximum(2000);
  xframeMet30_bg->Draw();
  //line125->DrawLine(125.3,0,125.3,xframeMet30_bg->GetMaximum());
  StackMet30->Draw("histoSAME");
  xframeMet30_bg->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/Higgs/ggInvarMassWithMET30cut_bgOnly.png");

  //

  xMet30_bg.setRange("SigFull",sigLo,sigHi);
  RooAbsReal* ig30Sig = Met30Pdf_bg.createIntegral(xMet30_bg,NormSet(xMet30_bg),Range("SigFull"));
  RooAbsReal* ig30SigNoNorm = Met30Pdf_bg.createIntegral(xMet30_bg,Range("SigFull"));
  RooAbsReal* ig30Full = Met30Pdf_bg.createIntegral(xMet30_bg);
  RooAbsReal* ig30FullMinusSig = Met30Pdf_bg.createIntegral(xMet30_bg,Range("sb_30lo,sb_30hi"));
  RooAbsReal* ig30SBlo = Met30Pdf_bg.createIntegral(xMet30_bg,NormSet(xMet30_bg),Range("sb_30lo"));
  RooAbsReal* ig30SBhi = Met30Pdf_bg.createIntegral(xMet30_bg,NormSet(xMet30_bg),Range("sb_30hi"));

  float Higg30SigNoNorm=ig30Sig->getVal()*Bern30Yield_bg.getVal();
  float Higg30LowSBNoNorm=ig30SBlo->getVal()*Bern30Yield_bg.getVal(),Higg30HighSBNoNorm=ig30SBhi->getVal()*Bern30Yield_bg.getVal();
  cout<<"Higg30SigNoNorm   : "<<Higg30SigNoNorm<<endl;
  cout<<"Higg30LowSBNoNorm : "<<Higg30LowSBNoNorm<<endl;
  cout<<"Higg30HighSBNoNorm: "<<Higg30HighSBNoNorm<<endl;
  vector<RooAbsReal*> integralInts30(0);
  h_ggInvMassSubtract30->Rebin(2);
  for(int i=1;i<h_ggInvMassSubtract30->GetBinLowEdge(h_ggInvMassSubtract30->GetNbinsX());i++){//does this only work if the bins are integers?
    TString name = "count30";name+=i;
    xMet30_bg.setRange(name,h_ggInvMassSubtract30->GetBinLowEdge(i),h_ggInvMassSubtract30->GetBinLowEdge(i+1));
    integralInts30.push_back(Met30Pdf_bg.createIntegral(xMet30_bg,NormSet(xMet30_bg),Range(name)));
  }

  h_ggInvMassSubtract30->GetXaxis()->SetRangeUser(100,164.9);
  h_ggInvMassSubtract30->Draw();
  c1->Print("Plots/Higgs/ggInvarMassWithMET30cut_bgOnly_presubtract.png");
  int bin30_120=h_ggInvMassSubtract30->FindBin(sigLo),bin30_130=h_ggInvMassSubtract30->FindBin(sigHi)-1;
  cout<<"full int:"<<ig30Sig->getVal()<<"  full int no norm:"<<ig30SigNoNorm->getVal()<<"  full int * yield:"<<ig30Sig->getVal()*Bern30Yield_bg.getVal()<<"  data integral: "<<h_ggInvMassSubtract30->Integral(bin30_120,bin30_130)<<endl;
  for(int i=1;i<h_ggInvMassSubtract30->GetNbinsX()-1;i++){
    float x=0,xe=0;
    x=h_ggInvMassSubtract30->GetBinContent(i);xe=h_ggInvMassSubtract30->GetBinError(i);
    h_ggInvMassSubtract30->SetBinContent(i,x-integralInts30[i-1]->getVal()*Bern30Yield_bg.getVal());
    //x=h_ggInvMassSubtract30->GetBinContent(i);
    h_ggInvMassSubtract30->SetBinError(i,sqrt(xe*xe+sqrt(integralInts30[i-1]->getVal()*Bern30Yield_bg.getVal()*integralInts30[i-1]->getVal()*Bern30Yield_bg.getVal())));
  }
  h_ggInvMassSubtract30->GetXaxis()->SetRangeUser(100,164.9);//h_ggInvMassSubtract30->GetYaxis()->SetRangeUser(-60,200);
  h_ggInvMassSubtract30->Draw();StackMet30->Draw("histoSAMES");h_ggInvMassSubtract30->Draw("SAMES");
  c1->Print("Plots/Higgs/ggInvarMassWithMET30cut_bgOnly_subtract.png");

  //

  gStyle->SetOptFit(0);

  sbFitHiHi=163;

  TH1F* DataProjXpow = (TH1F*)h_ggInvarMassMET->Clone();
  TH1F* DataProjXlin = (TH1F*)h_ggInvarMassMET->Clone();
  TH1F* DataProjXexpo = (TH1F*)h_ggInvarMassMET->Clone();
  TH1F* DataProjXpol = (TH1F*)h_ggInvarMassMET->Clone();

  reject=true;
  TF1* fitCurve = new TF1("fitCurve",fpow,sbFitLoLo,sbFitHiHi,2);
  Double_t avg_l = DataProjXpow->Integral(DataProjXpow->FindBin(sbFitLoLo),DataProjXpow->FindBin(sbLoHi))/float(sbLoHi-sbFitLoLo),avg_u = DataProjXpow->Integral(DataProjXpow->FindBin(sbHiLo),DataProjXpow->FindBin(sbFitHiHi))/float(sbFitHiHi-sbHiLo),avgX_l=(sbLoHi-sbFitLoLo)/2.,avgX_u=(sbFitHiHi-sbHiLo)/2.;
  cout<<avg_l<<"  "<<avg_u<<"  "<<avgX_l<<"  "<<avgX_u<<endl;
  Double_t param1= (log(avg_l) - log(avg_u))/(log(avgX_l) - log(avgX_u));
  Double_t param0= /*7e14;*/avg_l/pow(avgX_l, param1);
  cout<<"param0: "<<param0<<"  param1: "<<param1<<endl;
  fitCurve->SetParameter(0,param0);
  fitCurve->SetParameter(1,param1);
  //first just to check the status
  int status = DataProjXpow->Fit(fitCurve,"L","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResult = DataProjXpow->Fit(fitCurve,"SLLMEV0","",sbFitLoLo,sbFitHiHi);
  TMatrixDSym cov = fitResult->GetCovarianceMatrix();
  fitResult->Print("V");
  DataProjXpow->GetXaxis()->SetRangeUser(100,164.9);
  reject=false;
  float YieldBinWidth=DataProjXpow->GetBinWidth(1);
  float PowYield = fitCurve->Integral(sbLoLo,sbHiHi)/YieldBinWidth;
  float PowYieldErr = fitCurve->IntegralError(sbLoLo,sbHiHi,fitResult->GetParams(),cov.GetMatrixArray() )/YieldBinWidth; 
  float PowYieldSig = fitCurve->Integral(sigLo,sigHi)/YieldBinWidth;
  float PowYieldSigErr = fitCurve->IntegralError(sigLo,sigHi,fitResult->GetParams(),cov.GetMatrixArray() )/YieldBinWidth; 
  reject=true;
  char str[50],strSig[50];
  sprintf(str,"Fit Yield: %4.2f #pm %4.2f",PowYield,PowYieldErr);
  sprintf(strSig,"Fit Signal Yield: %4.2f #pm %4.2f",PowYieldSig,PowYieldSigErr);
  //char strRF[50];
  //sprintf(strRF,"RooFit signal Yield: %4.2f #pm %4.2f",integ_,dinteg_);
  TPaveText *text= new TPaveText(.5,.65,.84,.84,"NDC");text->SetFillStyle(0);text->SetBorderSize(0);
  //text->AddText(str);
  text->AddText(strSig);text->SetTextFont(42);
  //text->AddText(strRF);
  DataProjXpow->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  DataProjXpow->SetTitle("");
  DataProjXpow->GetYaxis()->SetRangeUser(0,5199);
  DataProjXpow->GetYaxis()->SetTitle("Events / GeV");
  DataProjXpow->SetLabelFont(42,"xy");
  DataProjXpow->SetTitleFont(42,"xy");
  DataProjXpow->Draw();
  TF1 *fitCurve2 = new TF1("fitCurve2","[0]*pow(x,[1])",sbLoHi,sbHiLo);
  fitCurve2->SetParameters(fitCurve->GetParameters());
  fitCurve2->SetLineColor(kBlue);fitCurve2->Draw("SAMES");
  TF1 *fitCurveLeft = new TF1("fitCurveLeft","[0]*pow(x,[1])",sbLoLo,sbLoHi);fitCurveLeft->SetParameters(fitCurve->GetParameters());
  DataProjXpow->GetListOfFunctions()->Add(fitCurveLeft);//gROOT->GetListOfFunctions()->Remove(fitCurveLeft);
  fitCurveLeft->SetLineColor(kRed);fitCurveLeft->Draw("SAMES");
  TF1 *fitCurveRight = new TF1("fitCurveRight","[0]*pow(x,[1])",sbHiLo,sbHiHi);fitCurveRight->SetParameters(fitCurve->GetParameters());
  DataProjXpow->GetListOfFunctions()->Add(fitCurveRight);//gROOT->GetListOfFunctions()->Remove(fitCurveRight);
  fitCurveRight->SetLineColor(kRed);fitCurveRight->Draw("SAMES");
  DataProjXpow->Draw("SAMES");
  //text->Draw("SAMES");
  TPaveText *ggTextInvMass = new TPaveText(.3,.7,.61,.81,"NDC");
  ggTextInvMass->AddText("#gamma#gamma");
  ////ggText->AddText("");
  ggTextInvMass->SetFillStyle(0);ggTextInvMass->SetTextFont(42);
  ggTextInvMass->SetFillColor(0);
  ggTextInvMass->SetBorderSize(0);
  ggTextInvMass->Draw();
  TLegend *InvMassLeg = new TLegend(.56,.6,.85,.85,"","brNDC");
  InvMassLeg->SetFillColor(kWhite);InvMassLeg->SetFillStyle(0);InvMassLeg->SetBorderSize(0);
  InvMassLeg->AddEntry(DataProjXpow,"Data","pe");
  InvMassLeg->AddEntry(fitCurveLeft,"Power Law Fit","l");
  //InvMassLeg->AddEntry(fitCurve2,"Fit Extrapolation","l");
  InvMassLeg->Draw();
  TPaveText *PrelimTextInvMass = new TPaveText(0.06,.91,.98,.965,"NDC");
  //PrelimTextInvMass->AddText("CMS Preliminary                  #sqrt{s}=8 TeV, L = 19.5 fb^{-1}");
  PrelimTextInvMass->AddText("CMS        L = 19.5 fb^{-1}        #sqrt{s} = 8 TeV");
  PrelimTextInvMass->SetFillStyle(0);PrelimTextInvMass->SetTextFont(42);
  PrelimTextInvMass->SetFillColor(0);
  PrelimTextInvMass->SetBorderSize(0);
  PrelimTextInvMass->Draw();
  c1->Print("Plots/Higgs/Exclusive_DataInvMassFit_powLaw.png");
  c1->Print("Plots/Higgs/Exclusive_DataInvMassFit_powLaw.pdf");

  reject=true;
  TF1* fitCurveExpo = new TF1("fitCurveExpo",fexp,sbFitLoLo,sbFitHiHi,2);
  fitCurveExpo->SetParameter(0,40.);
  fitCurveExpo->SetParameter(1,-2);
  //first just to check the status
  status = DataProjXexpo->Fit(fitCurveExpo,"L","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResultExpo = DataProjXexpo->Fit(fitCurveExpo,"SL","",sbFitLoLo,sbFitHiHi);
  TMatrixDSym covExpo = fitResultExpo->GetCovarianceMatrix();
  fitResultExpo->Print("V");
  DataProjXexpo->GetXaxis()->SetRangeUser(100,164.9);
  reject=false;
  YieldBinWidth=DataProjXexpo->GetBinWidth(1);
  float ExpoYield = fitCurveExpo->Integral(sbLoLo,sbHiHi)/YieldBinWidth;
  float ExpoYieldErr = fitCurveExpo->IntegralError(sbLoLo,sbHiHi,fitResultExpo->GetParams(),covExpo.GetMatrixArray() )/YieldBinWidth; 
  float ExpoYieldSig = fitCurveExpo->Integral(sigLo,sigHi)/YieldBinWidth;
  float ExpoYieldSigErr = fitCurveExpo->IntegralError(sigLo,sigHi,fitResultExpo->GetParams(),covExpo.GetMatrixArray() )/YieldBinWidth; 
  reject=true;
  sprintf(str,"Fit Yield: %4.2f #pm %4.2f",ExpoYield,ExpoYieldErr);
  sprintf(strSig,"Fit Signal Yield: %4.2f #pm %4.2f",ExpoYieldSig,ExpoYieldSigErr);
  //sprintf(strRF,"RooFit signal Yield: %4.2f #pm %4.2f",integ_,dinteg_);
  TPaveText *textExpo= new TPaveText(.5,.65,.84,.84,"NDC");textExpo->SetFillStyle(0);textExpo->SetBorderSize(0);
  //textExpo->AddText(str);
  textExpo->AddText(strSig);textExpo->SetTextFont(42);
  DataProjXexpo->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  DataProjXexpo->SetTitle("");
  DataProjXexpo->Draw();
  textExpo->Draw("SAMES");
  c1->Print("Plots/Higgs/Exclusive_DataInvMassFit_expo.png");

  reject=true;
  TF1* fitCurveLin = new TF1("fitCurveLin",flin,sbFitLoLo,sbFitHiHi,2);
  //fitCurveLin->FixParameter(0,13129.);
  //fitCurveLin->FixParameter(1,-64.5);
  //first just to check the status
  status = DataProjXlin->Fit(fitCurveLin,"L","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResultLin = DataProjXlin->Fit(fitCurveLin,"SL","",sbFitLoLo,sbFitHiHi);
  TMatrixDSym covLin = fitResultLin->GetCovarianceMatrix();
  fitResultLin->Print("V");
  DataProjXlin->GetXaxis()->SetRangeUser(100,164.9);
  reject=false;
  YieldBinWidth=DataProjXlin->GetBinWidth(1);
  float LinYield = fitCurveLin->Integral(sbLoLo,sbHiHi)/YieldBinWidth;
  float LinYieldErr = fitCurveLin->IntegralError(sbLoLo,sbHiHi,fitResultLin->GetParams(),covLin.GetMatrixArray() )/YieldBinWidth; 
  float LinYieldSig = fitCurveLin->Integral(sigLo,sigHi)/YieldBinWidth;
  float LinYieldSigErr = fitCurveLin->IntegralError(sigLo,sigHi,fitResultLin->GetParams(),covLin.GetMatrixArray() )/YieldBinWidth; 
  reject=true;
  sprintf(str,"Fit Yield: %4.2f #pm %4.2f",LinYield,LinYieldErr);
  sprintf(strSig,"Fit Signal Yield: %4.2f #pm %4.2f",LinYieldSig,LinYieldSigErr);
  //sprintf(strRF,"RooFit signal Yield: %4.2f #pm %4.2f",integ_,dinteg_);
  TPaveText *textLin= new TPaveText(.5,.65,.84,.84,"NDC");textLin->SetFillStyle(0);textLin->SetBorderSize(0);
  //textLin->AddText(str);
  textLin->AddText(strSig);textLin->SetTextFont(42);
  DataProjXlin->Draw();
  textLin->Draw("SAMES");
  c1->Print("Plots/Higgs/Exclusive_DataInvMassFit_linear.png");
  /*
  reject=true;
  TF1* fitCurvePol2 = new TF1("fitCurvePol2","pol2",sbFitLoLo,sbFitHiHi);
  fitCurvePol2->SetParameter(0,1e4);
  fitCurvePol2->SetParameter(1,-100);
  fitCurvePol2->SetParameter(2,10);
  //first just to check the status
  status = DataProjXpol->Fit(fitCurvePol2,"L","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResultPol2 = DataProjXpol->Fit(fitCurvePol2,"SL","",sbFitLoLo,sbFitHiHi);
  TMatrixDSym covPol2 = fitResultPol2->GetCovarianceMatrix();
  fitResultPol2->Print("V");
  DataProjXpol->GetXaxis()->SetRangeUser(100,164.9);
  reject=false;
  YieldBinWidth=DataProjXpol->GetBinWidth(1);
  float Pol2Yield = fitCurvePol2->Integral(sbLoLo,sbHiHi)/YieldBinWidth;
  float Pol2YieldErr = fitCurvePol2->IntegralError(sbLoLo,sbHiHi,fitResultPol2->GetParams(),covPol2.GetMatrixArray() )/YieldBinWidth; 
  float Pol2YieldSig = fitCurvePol2->Integral(sigLo,sigHi)/YieldBinWidth;
  float Pol2YieldSigErr = fitCurvePol2->IntegralError(sigLo,sigHi,fitResultPol2->GetParams(),covPol2.GetMatrixArray() )/YieldBinWidth; 
  reject=true;
  sprintf(str,"Fit Yield: %4.2f #pm %4.2f",Pol2Yield,Pol2YieldErr);
  sprintf(strSig,"Fit Signal Yield: %4.2f #pm %4.2f",Pol2YieldSig,Pol2YieldSigErr);
  //sprintf(strRF,"RooFit signal Yield: %4.2f #pm %4.2f",integ_,dinteg_);
  TPaveText *textPol2= new TPaveText(.5,.65,.84,.84,"NDC");textPol2->SetFillStyle(0);textPol2->SetBorderSize(0);
  textPol2->AddText(str);
  textPol2->AddText(strSig);
  DataProjXpol->Draw();
  textPol2->Draw("SAMES");
  c1->Print("Plots/Higgs/Exclusive_DataInvMassFit_pol2.png");
  */
  reject=true;
  TF1* fitCurvePol3 = new TF1("fitCurvePol3",fpol3,sbFitLoLo,sbFitHiHi,4);
  fitCurvePol3->SetParameter(0,100.);
  fitCurvePol3->SetParameter(1,-2.);
  //fitCurvePol3->SetParameter(2,1.);
  //first just to check the status
  status = DataProjXpol->Fit(fitCurvePol3,"L","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResultPol3 = DataProjXpol->Fit(fitCurvePol3,"SL","",sbFitLoLo,sbFitHiHi);
  TMatrixDSym covPol3 = fitResultPol3->GetCovarianceMatrix();
  fitResultPol3->Print("V");
  DataProjXpol->GetXaxis()->SetRangeUser(100,164.9);
  reject=false;
  YieldBinWidth=DataProjXpol->GetBinWidth(1);
  float Pol3Yield = fitCurvePol3->Integral(sbLoLo,sbHiHi)/YieldBinWidth;
  float Pol3YieldErr = fitCurvePol3->IntegralError(sbLoLo,sbHiHi,fitResultPol3->GetParams(),covPol3.GetMatrixArray() )/YieldBinWidth; 
  float Pol3YieldSig = fitCurvePol3->Integral(sigLo,sigHi)/YieldBinWidth;
  float Pol3YieldSigErr = fitCurvePol3->IntegralError(sigLo,sigHi,fitResultPol3->GetParams(),covPol3.GetMatrixArray() )/YieldBinWidth; 
  reject=true;
  sprintf(str,"Fit Yield: %4.2f #pm %4.2f",Pol3Yield,Pol3YieldErr);
  sprintf(strSig,"Fit Signal Yield: %4.2f #pm %4.2f",Pol3YieldSig,Pol3YieldSigErr);
  //sprintf(strRF,"RooFit signal Yield: %4.2f #pm %4.2f",integ_,dinteg_);
  TPaveText *textPol3= new TPaveText(.5,.65,.84,.84,"NDC");textPol3->SetFillStyle(0);textPol3->SetBorderSize(0);
  //textPol3->AddText(str);
  textPol3->AddText(strSig);textPol3->SetTextFont(42);
  DataProjXpol->Draw();
  textPol3->Draw("SAMES");
  c1->Print("Plots/Higgs/Exclusive_DataInvMassFit_pol3.png");

  reject=true;
  TF1* fitCurvePol4 = new TF1("fitCurvePol4","pol4",sbFitLoLo,sbFitHiHi);
  //first just to check the status
  status = DataProjXpol->Fit(fitCurvePol4,"L","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResultPol4 = DataProjXpol->Fit(fitCurvePol4,"SL","",sbFitLoLo,sbFitHiHi);
  TMatrixDSym covPol4 = fitResultPol4->GetCovarianceMatrix();
  fitResultPol4->Print("V");
  DataProjXpol->GetXaxis()->SetRangeUser(100,164.9);
  reject=false;
  YieldBinWidth=DataProjXpol->GetBinWidth(1);
  float Pol4Yield = fitCurvePol4->Integral(sbLoLo,sbHiHi)/YieldBinWidth;
  float Pol4YieldErr = fitCurvePol4->IntegralError(sbLoLo,sbHiHi,fitResultPol4->GetParams(),covPol4.GetMatrixArray() )/YieldBinWidth; 
  float Pol4YieldSig = fitCurvePol4->Integral(sigLo,sigHi)/YieldBinWidth;
  float Pol4YieldSigErr = fitCurvePol4->IntegralError(sigLo,sigHi,fitResultPol4->GetParams(),covPol4.GetMatrixArray() )/YieldBinWidth; 
  reject=true;
  sprintf(str,"Fit Yield: %4.2f #pm %4.2f",Pol4Yield,Pol4YieldErr);
  sprintf(strSig,"Fit Signal Yield: %4.2f #pm %4.2f",Pol4YieldSig,Pol4YieldSigErr);
  //sprintf(strRF,"RooFit signal Yield: %4.2f #pm %4.2f",integ_,dinteg_);
  TPaveText *textPol4= new TPaveText(.5,.65,.84,.84,"NDC");textPol4->SetFillStyle(0);textPol4->SetBorderSize(0);
  //textPol4->AddText(str);
  textPol4->AddText(strSig);textPol4->SetTextFont(42);
  DataProjXpol->Draw();
  textPol4->Draw("SAMES");
  c1->Print("Plots/Higgs/Exclusive_DataInvMassFit_pol4.png");

  reject=true;
  TF1* fitCurvePol5 = new TF1("fitCurvePol5",fpol5,sbFitLoLo,sbFitHiHi,5);
  fitCurvePol5->SetParameter(0,300.);
  fitCurvePol5->SetParameter(1,-5.);
  //first just to check the status
  status = DataProjXpol->Fit(fitCurvePol5,"L","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResultPol5 = DataProjXpol->Fit(fitCurvePol5,"SL","",sbFitLoLo,sbFitHiHi);
  TMatrixDSym covPol5 = fitResultPol5->GetCovarianceMatrix();
  fitResultPol5->Print("V");
  DataProjXpol->GetXaxis()->SetRangeUser(100,164.9);
  reject=false;
  YieldBinWidth=DataProjXpol->GetBinWidth(1);
  float Pol5Yield = fitCurvePol5->Integral(sbLoLo,sbHiHi)/YieldBinWidth;
  float Pol5YieldErr = fitCurvePol5->IntegralError(sbLoLo,sbHiHi,fitResultPol5->GetParams(),covPol5.GetMatrixArray() )/YieldBinWidth; 
  float Pol5YieldSig = fitCurvePol5->Integral(sigLo,sigHi)/YieldBinWidth;
  float Pol5YieldSigErr = fitCurvePol5->IntegralError(sigLo,sigHi,fitResultPol5->GetParams(),covPol5.GetMatrixArray() )/YieldBinWidth; 
  reject=true;
  sprintf(str,"Fit Yield: %4.2f #pm %4.2f",Pol5Yield,Pol5YieldErr);
  sprintf(strSig,"Fit Signal Yield: %4.2f #pm %4.2f",Pol5YieldSig,Pol5YieldSigErr);
  //sprintf(strRF,"RooFit signal Yield: %4.2f #pm %4.2f",integ_,dinteg_);
  TPaveText *textPol5= new TPaveText(.5,.65,.84,.84,"NDC");textPol5->SetFillStyle(0);textPol5->SetBorderSize(0);
  //textPol5->AddText(str);
  textPol5->AddText(strSig);textPol5->SetTextFont(42);
  DataProjXpol->Draw();
  TF1 *fitCurvePol5_2 = new TF1("fitCurvePol5_2","[0]*x+[1]*x*x+[2]*x*x*x+[3]*x*x*x*x+[4]*x*x*x*x*x",sbLoLo,sbHiHi);
  fitCurvePol5_2->SetParameters(fitCurvePol5->GetParameters());
  fitCurvePol5_2->SetLineColor(kBlue);fitCurvePol5_2->Draw("SAMES");
  DataProjXpol->Draw("SAMES");
  textPol5->Draw("SAMES");
  c1->Print("Plots/Higgs/Exclusive_DataInvMassFit_pol5.png");
  float param0WePol5Fix=fitCurvePol5->GetParameter(0);
  float param1WePol5Fix=fitCurvePol5->GetParameter(1);
  float param2WePol5Fix=fitCurvePol5->GetParameter(2);
  float param3WePol5Fix=fitCurvePol5->GetParameter(3);
  float param4WePol5Fix=fitCurvePol5->GetParameter(4);

  //higgs plot with ratio
  //p1->cd();p1->SetLogy(1);
  int sbLoBin1=h_ggMetData->GetXaxis()->FindBin(sbLoLo),sbLoBin2=h_ggMetData->GetXaxis()->FindBin(sbLoHi)-1;  
  int sbHiBin1=h_ggMetData->GetXaxis()->FindBin(sbHiLo),sbHiBin2=h_ggMetData->GetXaxis()->FindBin(sbHiHi)-1;
  int sigBin1=h_ggMetData->GetXaxis()->FindBin(sigLo),sigBin2=h_ggMetData->GetXaxis()->FindBin(sigHi)-1;
  TH1F* ggMet_SBloInvMass=(TH1F*)h_ggMetData->ProjectionY("ggMet_SBloInvMass",sbLoBin1,sbLoBin2,"eo");
  TH1F* ggMet_SBhiInvMass=(TH1F*)h_ggMetData->ProjectionY("ggMet_SBhiInvMass",sbHiBin1,sbHiBin2,"eo");
  TH1F* ggMet_sigInvMass=(TH1F*)h_ggMetData->ProjectionY("ggMet_sigInvMass",sigBin1,sigBin2,"eo");
  /*
  overFlowBin=ggMet_sigInvMass->FindBin(metPlotXmax+1), lastBin=ggMet_sigInvMass->FindBin(metPlotXmax-1);
  overFlow=0.;overFlowErr=0.;
  val=ggMet_sigInvMass->GetBinContent(lastBin);valErr=ggMet_sigInvMass->GetBinError(lastBin);
  overFlow=ggMet_sigInvMass->IntegralAndError(overFlowBin,-1,overFlowErr);
  ggMet_sigInvMass->SetBinContent(lastBin,val+overFlow);
  ggMet_sigInvMass->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  overFlow=0.;overFlowErr=0.;
  val=ggMet_SBloInvMass->GetBinContent(lastBin);valErr=ggMet_SBloInvMass->GetBinError(lastBin);
  overFlow=ggMet_SBloInvMass->IntegralAndError(overFlowBin,-1,overFlowErr);
  ggMet_SBloInvMass->SetBinContent(lastBin,val+overFlow);
  ggMet_SBloInvMass->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  overFlow=0.;overFlowErr=0.;
  val=ggMet_SBhiInvMass->GetBinContent(lastBin);valErr=ggMet_SBhiInvMass->GetBinError(lastBin);
  overFlow=ggMet_SBhiInvMass->IntegralAndError(overFlowBin,-1,overFlowErr);
  ggMet_SBhiInvMass->SetBinContent(lastBin,val+overFlow);
  ggMet_SBhiInvMass->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  */
  TH1F* ggMet_SBloInvMassRebin = (TH1F*)ggMet_SBloInvMass->Rebin(NmetBins,"ggMet_SBloInvMassRebin",xbins);
  TH1F* ggMet_SBhiInvMassRebin = (TH1F*)ggMet_SBhiInvMass->Rebin(NmetBins,"ggMet_SBhiInvMassRebin",xbins);
  TH1F* ggMet_sigInvMassRebin = (TH1F*)ggMet_sigInvMass->Rebin(NmetBins,"ggMet_sigInvMassRebin",xbins);

  AddOverflowToLastBin(ggMet_SBloInvMassRebin);AddOverflowToLastBin(ggMet_SBhiInvMassRebin);
  AddOverflowToLastBin(ggMet_sigInvMassRebin);

  float higgslowscale = ggMet_SBloInvMass->Integral() ? PowYieldSig/*(HiggSigNoNorm)*//ggMet_SBloInvMass->Integral() : 0;
  float higgshighscale = ggMet_SBhiInvMass->Integral() ? PowYieldSig/*(HiggSigNoNorm)*//ggMet_SBhiInvMass->Integral() : 0;
  cout<<"Higgs signal integral: "<<ggMet_sigInvMass->Integral()<<endl;
  cout<<"Higgs SBlow  integral: "<<ggMet_SBloInvMass->Integral()<<endl;
  cout<<"Higgs SBhigh integral: "<<ggMet_SBhiInvMass->Integral()<<endl;
  //float FitSyst=(BernYield_bg.getError()*BernYield_bg.getError())/(BernYield_bg.getVal()*BernYield_bg.getVal());
  float FitSyst=(PowYieldSigErr*PowYieldSigErr)/(PowYieldSig*PowYieldSig);
  double bgStatErrLO=0.;double bgStatLO=ggMet_SBloInvMassRebin->IntegralAndError(0,-1,bgStatErrLO);
  float bgStatSystLO=(bgStatErrLO*bgStatErrLO)/(bgStatLO*bgStatLO);
  double bgStatErrHI=0.;double bgStatHI=ggMet_SBhiInvMassRebin->IntegralAndError(0,-1,bgStatErrHI);
  float bgStatSystHI=(bgStatErrHI*bgStatErrHI)/(bgStatHI*bgStatHI);
  /*
  for(int i=0;i<ggMet_SBloInvMassRebin->GetNbinsX();i++){
    float Value=ggMet_SBloInvMassRebin->GetBinContent(i),err=ggMet_SBloInvMassRebin->GetBinError(i);
    float ValueNew=Value*higgslowscale;
    float errNew=ValueNew*sqrt((err*err/Value/Value)+FitSyst+bgStatSystLO);
    ggMet_SBloInvMassRebin->SetBinContent(i,ValueNew);ggMet_SBloInvMassRebin->SetBinError(i,errNew);
    Value=ggMet_SBhiInvMassRebin->GetBinContent(i),err=ggMet_SBhiInvMassRebin->GetBinError(i);
    ValueNew=Value*higgshighscale;
    errNew=ValueNew*sqrt((err*err/Value/Value)+FitSyst+bgStatSystHI);
    ggMet_SBhiInvMassRebin->SetBinContent(i,ValueNew);ggMet_SBhiInvMassRebin->SetBinError(i,errNew);
    }*/
  TH1F* ggMet_InvMassSBcombNoRebin=(TH1F*)ggMet_SBloInvMass->Clone();//ggMet_InvMassSBcombNoRebin->Add(ggMet_SBhiInvMass);ggMet_InvMassSBcombNoRebin->Scale(0.5);
  TH1F *StatErrs       =(TH1F*)ggMet_SBloInvMass->Clone();
  TH1F *FitStatSystErrs =(TH1F*)ggMet_SBloInvMass->Clone();
  TH1F *FitShapeSystErrs  =(TH1F*)ggMet_SBloInvMass->Clone();
  TH1F *HalfDiffErrs  =(TH1F*)ggMet_SBloInvMass->Clone();
  TH1F* ggMet_SBloInvMassClone = (TH1F*)ggMet_SBloInvMass->Clone("ggMet_SBloInvMassClone");
  TH1F* ggMet_SBhiInvMassClone = (TH1F*)ggMet_SBhiInvMass->Clone("ggMet_SBhiInvMassClone");
  ggMet_SBloInvMassClone->Scale(higgslowscale);ggMet_SBhiInvMassClone->Scale(higgshighscale);
  TH1F* ggMet_SBloInvMassCloneRebin = (TH1F*)ggMet_SBloInvMassClone->Rebin(NmetBins,"ggMet_SBloInvMassCloneRebin",xbins);
  TH1F* ggMet_SBhiInvMassCloneRebin = (TH1F*)ggMet_SBhiInvMassClone->Rebin(NmetBins,"ggMet_SBhiInvMassCloneRebin",xbins);
 
  for(int i=1;i<=ggMet_InvMassSBcombNoRebin->GetNbinsX();i++){
    Double_t val = ((ggMet_SBloInvMass->GetBinContent(i)*higgslowscale)+(ggMet_SBhiInvMass->GetBinContent(i)*higgshighscale))/2.;
    //double Lerr=0.,Uerr=0.;
    //double Lval = ggMet_SBloInvMass->IntegralAndError(0,-1,Lerr);
    //double Uval = ggMet_SBhiInvMass->IntegralAndError(0,-1,Uerr);
    Double_t LvalScale = ggMet_SBloInvMass->GetBinContent(i)*higgslowscale,UvalScale = ggMet_SBhiInvMass->GetBinContent(i)*higgshighscale;
    Double_t StatErr(0.),FitStatSystErr(0.),FitShapeSystErr(0.),HalfDiffErr(0.);
    double err = TotErr(val,ggMet_SBloInvMass->GetBinContent(i),ggMet_SBloInvMass->GetBinError(i),ggMet_SBhiInvMass->GetBinContent(i),ggMet_SBhiInvMass->GetBinError(i),bgStatLO,bgStatErrLO,bgStatHI,bgStatErrHI,PowYieldSig,PowYieldSigErr,ExpoYieldSig,StatErr,FitStatSystErr,FitShapeSystErr,LvalScale,UvalScale,HalfDiffErr);
    ggMet_InvMassSBcombNoRebin->SetBinContent(i,val);ggMet_InvMassSBcombNoRebin->SetBinError(i,err);
    StatErrs->SetBinError(i,StatErr);FitStatSystErrs->SetBinError(i,FitStatSystErr);FitShapeSystErrs->SetBinError(i,FitShapeSystErr);HalfDiffErrs->SetBinError(i,HalfDiffErr);
  }
  TH1F* ggMet_InvMassSBcomb = (TH1F*)ggMet_InvMassSBcombNoRebin->Rebin(NmetBins,"ggMet_InvMassSBcomb",xbins);
  TH1F* ggMet_sigInvMassRebin_Clone = (TH1F*)ggMet_sigInvMassRebin->Clone("ggMet_sigInvMassRebin_Clone");
  AddOverflowToLastBin(ggMet_InvMassSBcomb);
  DivideByBinWidth(ggMet_InvMassSBcomb);
  DivideByBinWidth(ggMet_sigInvMassRebin);

  c1->cd();c1->SetLogy(1);
  ggMet_SBloInvMassCloneRebin->SetFillColor(kRed);ggMet_SBloInvMassCloneRebin->SetFillStyle(3004);ggMet_SBloInvMassCloneRebin->SetMarkerSize(0);
  ggMet_SBloInvMassCloneRebin->GetXaxis()->SetRangeUser(0,250);
  ggMet_SBloInvMassCloneRebin->SetTitle("");ggMet_SBloInvMassCloneRebin->GetXaxis()->SetTitle("E_{T}^{miss}");
  ggMet_SBloInvMassCloneRebin->GetYaxis()->SetTitle("Events");
  ggMet_SBloInvMassCloneRebin->Draw("E2");ggMet_SBloInvMassCloneRebin->Draw("pesames");
  ggMet_SBhiInvMassCloneRebin->SetFillColor(kBlue);ggMet_SBhiInvMassCloneRebin->SetFillStyle(3004);ggMet_SBhiInvMassCloneRebin->SetMarkerSize(0);
  ggMet_SBloInvMassCloneRebin->SetLineColor(kRed);ggMet_SBloInvMassCloneRebin->SetMarkerColor(kRed);
  ggMet_SBhiInvMassCloneRebin->SetLineColor(kBlue);ggMet_SBhiInvMassCloneRebin->SetMarkerColor(kBlue);
  ggMet_SBloInvMassCloneRebin->SetLineWidth(2);ggMet_SBhiInvMassCloneRebin->SetLineWidth(2);
  ggMet_SBhiInvMassCloneRebin->Draw("E2sames");ggMet_SBhiInvMassCloneRebin->Draw("pesames");
  TLegend *lowhighleg = new TLegend(.55,.5,.85,.8);
  lowhighleg->AddEntry(ggMet_SBloInvMassCloneRebin,"low sideband","plf");
  lowhighleg->AddEntry(ggMet_SBhiInvMassCloneRebin,"high sideband","plf");
  lowhighleg->AddEntry(ggMet_sigInvMassRebin_Clone,"higgs mass window","ple");
  lowhighleg->SetFillStyle(0);lowhighleg->SetBorderSize(0);
  ggMet_sigInvMassRebin_Clone->Draw("PESAMES");
  lowhighleg->Draw();
  c1->Print("Plots/Higgs/Inclusive_Met_SbLowAndHigh.png");
  c1->Print("Plots/Higgs/Inclusive_Met_SbLowAndHigh.pdf");
  c1->SetLogy(0);
  p1->cd();p1->SetLogy(1);

  /*
  for(int i=0;i<ggMet_InvMassSBcomb->GetNbinsX()+2;i++){
    float x = ggMet_InvMassSBcomb->GetBinContent(i)/ggMet_InvMassSBcomb->GetBinWidth(i);ggMet_InvMassSBcomb->SetBinContent(i,x);
    x = ggMet_InvMassSBcomb->GetBinError(i)/ggMet_InvMassSBcomb->GetBinWidth(i);ggMet_InvMassSBcomb->SetBinError(i,x);
    x = ggMet_sigInvMassRebin->GetBinContent(i)/ggMet_sigInvMassRebin->GetBinWidth(i);ggMet_sigInvMassRebin->SetBinContent(i,x);
    x = ggMet_sigInvMassRebin->GetBinError(i)/ggMet_sigInvMassRebin->GetBinWidth(i);ggMet_sigInvMassRebin->SetBinError(i,x);
  }
  */
  //ggMet_SBloInvMassRebin->Scale(higgslowscale);
  //ggMet_SBhiInvMassRebin->Scale(higgshighscale);
  cout<<"after scaling --------"<<endl;
  cout<<"Higgs signal integral: "<<ggMet_sigInvMassRebin->Integral()<<endl;
  cout<<"Higgs SBlow  integral: "<<ggMet_SBloInvMassRebin->Integral()<<endl;
  cout<<"Higgs SBhigh integral: "<<ggMet_SBhiInvMassRebin->Integral()<<endl;
  /*
  for(int i=1;i<=ggMet_sigInvMassRebin->GetNbinsX();i++){
    float x = ggMet_SBloInvMassRebin->GetBinContent(i)/ggMet_SBloInvMassRebin->GetBinWidth(i);ggMet_SBloInvMassRebin->SetBinContent(i,x);
    x = ggMet_SBloInvMassRebin->GetBinError(i)/ggMet_SBloInvMassRebin->GetBinWidth(i);ggMet_SBloInvMassRebin->SetBinError(i,x);
    x = ggMet_SBhiInvMassRebin->GetBinContent(i)/ggMet_SBhiInvMassRebin->GetBinWidth(i);ggMet_SBhiInvMassRebin->SetBinContent(i,x);
    x = ggMet_SBhiInvMassRebin->GetBinError(i)/ggMet_SBhiInvMassRebin->GetBinWidth(i);ggMet_SBhiInvMassRebin->SetBinError(i,x);
    x = ggMet_sigInvMassRebin->GetBinContent(i)/ggMet_sigInvMassRebin->GetBinWidth(i);ggMet_sigInvMassRebin->SetBinContent(i,x);
    x = ggMet_sigInvMassRebin->GetBinError(i)/ggMet_sigInvMassRebin->GetBinWidth(i);ggMet_sigInvMassRebin->SetBinError(i,x);
    }*/
  //TH1F* ggMet_InvMassSBcomb=(TH1F*)ggMet_SBloInvMassRebin->Clone();ggMet_InvMassSBcomb->Add(ggMet_SBhiInvMassRebin);ggMet_InvMassSBcomb->Scale(0.5);

  ggMet_sigInvMassRebin->SetMarkerColor(kBlack);ggMet_sigInvMassRebin->SetLineColor(kBlack);
  ggMet_SBloInvMassRebin->SetMarkerColor(kBlue);ggMet_SBloInvMassRebin->SetLineColor(kBlue);ggMet_SBloInvMassRebin->SetFillColor(kBlue);ggMet_SBloInvMassRebin->SetFillStyle(3004);
  ggMet_SBhiInvMassRebin->SetMarkerColor(kRed);ggMet_SBhiInvMassRebin->SetLineColor(kRed);ggMet_SBhiInvMassRebin->SetFillColor(kRed);ggMet_SBhiInvMassRebin->SetFillStyle(3004);
  ggMet_InvMassSBcomb->SetMarkerColor(42);ggMet_InvMassSBcomb->SetLineColor(42);ggMet_InvMassSBcomb->SetFillColor(42);//ggMet_InvMassSBcomb->SetFillStyle(3004);
  ggMet_SBloInvMassRebin->SetMarkerSize(0.75);ggMet_SBhiInvMassRebin->SetMarkerSize(0.75);ggMet_InvMassSBcomb->SetMarkerSize(0.75);
  
  //THStack *HiggsBernStack = new THStack("HiggsBernStack","");
  //HiggsBernStack->Add(TTHggMetNew2);
  HiggsMetStack->Add(ggMet_InvMassSBcomb);
  TH1F* ggMet_InvMassSBcomb_plus_SMhiggs = (TH1F*)ggMet_InvMassSBcomb->Clone();ggMet_InvMassSBcomb_plus_SMhiggs->Add(SMHiggsMet);
  ggMet_InvMassSBcomb_plus_SMhiggs->SetMarkerSize(0);ggMet_InvMassSBcomb_plus_SMhiggs->SetFillColor(kBlack);ggMet_InvMassSBcomb_plus_SMhiggs->SetLineColor(kBlack);ggMet_InvMassSBcomb_plus_SMhiggs->SetFillStyle(3004);

  ggMet_sigInvMassRebin->GetXaxis()->SetRangeUser(0,249);
  ggMet_sigInvMassRebin->GetXaxis()->SetTitle("");
  ggMet_sigInvMassRebin->GetYaxis()->SetRangeUser(1e-4,5000);
  ggMet_sigInvMassRebin->GetYaxis()->SetTitle("Events / GeV");
  ggMet_sigInvMassRebin->GetYaxis()->SetTitleSize(0.05);
  ggMet_sigInvMassRebin->SetTitle("");
  ggMet_sigInvMassRebin->GetXaxis()->SetLabelSize(0);
  cout<<"Inclusive signal region full integral: "<<ggMet_sigInvMassRebin->Integral(0,-1)<<endl;
  //MakeBlindMet(ggMet_sigInvMassRebin);
  ggMet_sigInvMassRebin->Draw("PE");
  HiggsMetStack->Draw("histoSAMES");
  //MakeBlindMet(ggMet_InvMassSBcomb_plus_SMhiggs);
  ggMet_InvMassSBcomb_plus_SMhiggs->Draw("E2SAMES");
  double SMest=0.,SMestErr=0.;SMest=ggMet_InvMassSBcomb_plus_SMhiggs->IntegralAndError(0,-1,SMestErr);
  cout<<"Inclusive background estimate full integral: "<<SMest<<"  and error: "<<SMestErr<<endl;
  cout<<"Inclusive SM higgs full integral: "<<SMHiggsMet->Integral(0,-1)<<endl;
  ggMet_sigInvMassRebin->Draw("PESAMES");
  //AAW130ProjYRebin->Draw("histoSAMES");
  //AAW275ProjYRebin->Draw("histoSAMES");
  /*h_SMS_WH_gg_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_ZH_gg_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_HH_2b2g_gg_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_HH_2W2g_gg_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_HH_2Z2g_gg_Excluded_Rebin->Draw("histoSAMES");*/

  /*
  //TPaveText *PrelimText = new TPaveText(.3,.86,.89,.9,"NDC");
  TPaveText *PrelimText = new TPaveText(0.,.86,.89,1.,"NDC");
  //PrelimText->AddText("CMS Preliminary                  #sqrt{s}=8 TeV, L = 19.5 fb^{-1}");
  PrelimText->AddText("CMS  L = 19.5 fb^{-1}  #sqrt{s} = 8 TeV");
  PrelimText->SetFillStyle(0);
  PrelimText->SetFillColor(0);
  PrelimText->SetBorderSize(0);
  PrelimText->Draw();
  */
  TPaveText *ggText = new TPaveText(.28,.79,.58,.9,"NDC");
  ggText->AddText("Inclusive #gamma#gamma");
  ////ggText->AddText("");
  ggText->SetFillStyle(0);ggText->SetTextFont(42);
  ggText->SetFillColor(0);
  ggText->SetBorderSize(0);
  ggText->Draw();
  //TLegend *MetLeg = new TLegend(.47,.35,.82,.76,"","brNDC");
  TLegend *MetLeg = new TLegend(.59,.66,.93,.91,"","brNDC");
  MetLeg->SetFillColor(kWhite);
  MetLeg->AddEntry(ggMet_sigInvMassRebin,"Data","elpz");
  MetLeg->AddEntry(ggMet_InvMassSBcomb,"Non-higgs SM bg","f");
  MetLeg->AddEntry(ggMet_InvMassSBcomb_plus_SMhiggs,"Full SM background error","f");
  /*
  MetLeg->AddEntry(ggHggMetRebin,"SM GluGluH, m_{h}=126 GeV","f");
  MetLeg->AddEntry(VBFHggMetRebin,"SM VBFH, m_{h}=126 GeV","f");
  MetLeg->AddEntry(WZHggMetRebin,"SM W/ZH, m_{h}=126 GeV","f");
  MetLeg->AddEntry(TTHggMetRebin,"SM TTH, m_{h}=126 GeV","f");
  */
  TLegend *MetLeg2 = new TLegend(.595,.495,.94,.675,"","brNDC");
  MetLeg2->SetNColumns(2);
  MetLeg2->SetFillColor(kWhite);
  MetLeg2->AddEntry(ggHggMetRebin,"SM GluGluH","f");
  MetLeg2->AddEntry(VBFHggMetRebin,"SM VBFH","f");
  MetLeg2->AddEntry(WZHggMetRebin,"SM W/ZH","f");
  MetLeg2->AddEntry(TTHggMetRebin,"SM TTH","f");
  /*MetLeg->AddEntry(h_SMS_WH_gg_Excluded_Rebin,"SMS WH m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  MetLeg->AddEntry(h_SMS_ZH_gg_Excluded_Rebin,"SMS ZH m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  MetLeg->AddEntry(h_SMS_HH_2b2g_gg_Excluded_Rebin,"SMS H(bb)H m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  MetLeg->AddEntry(h_SMS_HH_2W2g_gg_Excluded_Rebin,"SMS H(WW)H m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  MetLeg->AddEntry(h_SMS_HH_2Z2g_gg_Excluded_Rebin,"SMS H(ZZ)H m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");*/
  //MetLeg->AddEntry(AAW130MetRebin,"AAW m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  //MetLeg->AddEntry(AAW275ProjYRebin,"AAW m_{#tilde{H}}=275, m_{#tilde{B}}=1","l");
  MetLeg->SetFillStyle(0);MetLeg->SetBorderSize(0);
  MetLeg2->SetFillStyle(0);MetLeg2->SetBorderSize(0);
  MetLeg->Draw();MetLeg2->Draw();
  //TPaveText *PrelimText = new TPaveText(.3,.86,.89,.9,"NDC");
  TPaveText *PrelimText = new TPaveText(0.08,.935,1.01,.99,"NDC");
  //PrelimText->AddText("CMS Preliminary                  #sqrt{s}=8 TeV, L = 19.5 fb^{-1}");
  PrelimText->AddText("CMS        L = 19.5 fb^{-1}        #sqrt{s} = 8 TeV");
  PrelimText->SetTextFont(42);
  PrelimText->SetFillStyle(0);
  PrelimText->SetFillColor(0);
  PrelimText->SetBorderSize(0);
  PrelimText->Draw();
  p1->RedrawAxis();
  p2->cd();
  TH1F* ggMet_sigInvMassRebin_Div_low = (TH1F*)ggMet_sigInvMassRebin->Clone();
  TH1F* ggMet_sigInvMassRebin_Div_high = (TH1F*)ggMet_sigInvMassRebin->Clone();
  TH1F* ggMet_sigInvMassRebin_Div_comb = (TH1F*)ggMet_sigInvMassRebin->Clone();
  TH1F* ggMet_sigInvMassRebin_Div_comb_plusSMhiggs = (TH1F*)ggMet_sigInvMassRebin->Clone();
  ggMet_sigInvMassRebin_Div_low->Divide(ggMet_SBloInvMassRebin);
  ggMet_sigInvMassRebin_Div_high->Divide(ggMet_SBhiInvMassRebin);
  ggMet_sigInvMassRebin_Div_comb->Divide(ggMet_InvMassSBcomb);
  //ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->Divide(ggMet_InvMassSBcomb_plus_SMhiggs);
  TH1F* h_SystErr = (TH1F*)ggMet_InvMassSBcomb_plus_SMhiggs->Clone();h_SystErr->SetFillColor(kBlack);h_SystErr->SetFillStyle(3004);h_SystErr->SetMarkerSize(0);
  for(int i=1;i<ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->GetNbinsX()+1;i++){
    float Value = ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->GetBinContent(i);float StatErr = ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->GetBinError(i);
    Value/=ggMet_InvMassSBcomb_plus_SMhiggs->GetBinContent(i);StatErr/=ggMet_InvMassSBcomb_plus_SMhiggs->GetBinContent(i);
    ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->SetBinContent(i,Value);ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->SetBinError(i,StatErr);
    float SystErr=h_SystErr->GetBinError(i)/h_SystErr->GetBinContent(i);
    cout<<"h_SystErr->GetBinError(i): "<<h_SystErr->GetBinError(i)<<"  /  h_SystErr->GetBinContent(i): "<<h_SystErr->GetBinContent(i)<<" = SystErr: "<<SystErr<<endl;
    h_SystErr->SetBinContent(i,1); h_SystErr->SetBinError(i,SystErr);
  }
  ggMet_sigInvMassRebin_Div_low->SetLineColor(kBlue);ggMet_sigInvMassRebin_Div_high->SetMarkerColor(kBlue);ggMet_sigInvMassRebin_Div_high->SetMarkerSize(0.75);
  ggMet_sigInvMassRebin_Div_comb->SetLineColor(kViolet);ggMet_sigInvMassRebin_Div_comb->SetMarkerColor(kViolet);ggMet_sigInvMassRebin_Div_comb->SetMarkerSize(0.75);
  ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->SetLineColor(kBlack);ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->SetMarkerColor(kBlack);//ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->SetMarkerSize(0.75);
  ggMet_sigInvMassRebin_Div_high->SetLineColor(kBlack);ggMet_sigInvMassRebin_Div_high->SetMarkerColor(kBlack);ggMet_sigInvMassRebin_Div_high->SetMarkerSize(0.75);
  ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->GetYaxis()->SetTitle("");
  ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->SetTitle("");
  ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->GetYaxis()->SetRangeUser(0.48,2.);
  ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->SetTitle("");
  ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->GetYaxis()->SetTitle("#frac{Data}{Prediction}");
  ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->GetYaxis()->SetTitleOffset(0.5);
  ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->GetYaxis()->SetTitleSize(0.12);
  ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->GetXaxis()->SetTitleSize(0.2);
  ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->GetYaxis()->SetLabelSize(0.12);
  ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->GetXaxis()->SetLabelSize(0.13);
  //MakeBlindMet(ggMet_sigInvMassRebin_Div_comb_plusSMhiggs);
  //MakeBlindMet(h_SystErr);
  //ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->GetXaxis()->SetLabelSize(0.);
  ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->GetXaxis()->SetTitle("");
  ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->Draw("PE");
  h_SystErr->Draw("E2SAMES");
  ggMet_sigInvMassRebin_Div_comb_plusSMhiggs->Draw("PESAMES");
  //ggMet_sigInvMassRebin_Div_high->Draw("PEsames");
  //ggMet_sigInvMassRebin_Div_low->Draw("PEsames");
  TLine l1(0,1,250,1);l1.DrawLine(0,1,250,1);
  //PrelimText->Draw();
  p3->cd();p3->Clear();
  TPaveText *met2 = new TPaveText(.77,0,.97,1,"NDC");
  met2->AddText("E_{T}^{miss} [GeV]");
  met2->SetFillColor(0);met2->SetTextFont(42);
  met2->SetBorderSize(0);
  met2->SetTextSize(.6);
  met2->Draw();
  
  c2->Print("Plots/Higgs/ggMet_InvarMassCuts_BERNfit.png");
  c2->Print("Plots/Higgs/ggMet_InvarMassCuts_BERNfit.pdf");

  //higgs plot with ratio, met>30 cut
  p1->cd();p1->SetLogy(1);
  TH1F* ggMET30_SBloInvMass=(TH1F*)h_ggMetData->ProjectionY("ggMET30_SBloInvMass",sbLoBin1,sbLoBin2,"eo");
  TH1F* ggMET30_SBhiInvMass=(TH1F*)h_ggMetData->ProjectionY("ggMET30_SBhiInvMass",sbHiBin1,sbHiBin2,"eo");
  TH1F* ggMET30_sigInvMass=(TH1F*)h_ggMetData->ProjectionY("ggMET30_sigInvMass",sigBin1,sigBin2,"eo");
  TH1F* ggMet30_SBloInvMassRebin = (TH1F*)ggMET30_SBloInvMass->Rebin(NmetBins,"ggMet30_SBloInvMassRebin",xbins);
  TH1F* ggMet30_SBhiInvMassRebin = (TH1F*)ggMET30_SBhiInvMass->Rebin(NmetBins,"ggMet30_SBhiInvMassRebin",xbins);
  TH1F* ggMet30_sigInvMassRebin = (TH1F*)ggMET30_sigInvMass->Rebin(NmetBins,"ggMet30_sigInvMassRebin",xbins);
  cout<<"Higgs signal integral, MET>30 cut: "<<ggMET30_sigInvMass->Integral()<<endl;
  cout<<"Higgs SBlow  integral, MET>30 cut: "<<ggMET30_SBloInvMass->Integral()<<endl;
  cout<<"Higgs SBhigh integral, MET>30 cut: "<<ggMET30_SBhiInvMass->Integral()<<endl;
  int bin30cut=ggMet30_sigInvMassRebin->GetXaxis()->FindBin(30);
  for(int i=0;i<bin30cut;i++){
    ggMet30_sigInvMassRebin->SetBinContent(i,0);ggMet30_sigInvMassRebin->SetBinError(i,0);
    ggMet30_SBloInvMassRebin->SetBinContent(i,0);ggMet30_SBloInvMassRebin->SetBinError(i,0);
    ggMet30_SBhiInvMassRebin->SetBinContent(i,0);ggMet30_SBhiInvMassRebin->SetBinError(i,0);
  }
  float higgslowscale30 = ggMet30_SBloInvMassRebin->Integral() ? (HiggSigNoNorm)/ggMet30_SBloInvMassRebin->Integral() : 0;
  float higgshighscale30 = ggMet30_SBhiInvMassRebin->Integral() ? (HiggSigNoNorm)/ggMet30_SBhiInvMassRebin->Integral() : 0;
  ggMet30_SBloInvMassRebin->Scale(higgslowscale30);
  ggMet30_SBhiInvMassRebin->Scale(higgshighscale30);
  cout<<"after scaling --------"<<endl;
  cout<<"Higgs signal integral, MET>30 cut: "<<ggMet30_sigInvMassRebin->Integral()<<endl;
  cout<<"Higgs SBlow  integral, MET>30 cut: "<<ggMet30_SBloInvMassRebin->Integral()<<endl;
  cout<<"Higgs SBhigh integral, MET>30 cut: "<<ggMet30_SBhiInvMassRebin->Integral()<<endl;
  /*
  for(int i=1;i<ggMet30_sigInvMassRebin->GetNbinsX()+2;i++){
    float x = ggMet30_SBloInvMassRebin->GetBinContent(i)/ggMet30_SBloInvMassRebin->GetBinWidth(i);ggMet30_SBloInvMassRebin->SetBinContent(i,x);
    x = ggMet30_SBloInvMassRebin->GetBinError(i)/ggMet30_SBloInvMassRebin->GetBinWidth(i);ggMet30_SBloInvMassRebin->SetBinError(i,x);
    x = ggMet30_SBhiInvMassRebin->GetBinContent(i)/ggMet30_SBhiInvMassRebin->GetBinWidth(i);ggMet30_SBhiInvMassRebin->SetBinContent(i,x);
    x = ggMet30_SBhiInvMassRebin->GetBinError(i)/ggMet30_SBhiInvMassRebin->GetBinWidth(i);ggMet30_SBhiInvMassRebin->SetBinError(i,x);
    x = ggMet30_sigInvMassRebin->GetBinContent(i)/ggMet30_sigInvMassRebin->GetBinWidth(i);ggMet30_sigInvMassRebin->SetBinContent(i,x);
    x = ggMet30_sigInvMassRebin->GetBinError(i)/ggMet30_sigInvMassRebin->GetBinWidth(i);ggMet30_sigInvMassRebin->SetBinError(i,x);
  }
  */
  DivideByBinWidth(ggMet30_SBloInvMassRebin);
  DivideByBinWidth(ggMet30_SBhiInvMassRebin);
  DivideByBinWidth(ggMet30_sigInvMassRebin);
  TH1F* ggMET30_InvMassSBcomb=(TH1F*)ggMet30_SBloInvMassRebin->Clone();ggMET30_InvMassSBcomb->Add(ggMet30_SBhiInvMassRebin);ggMET30_InvMassSBcomb->Scale(0.5);

  ggMet30_sigInvMassRebin->SetMarkerColor(kBlack);ggMet30_sigInvMassRebin->SetLineColor(kBlack);
  ggMet30_SBloInvMassRebin->SetMarkerColor(kBlue);ggMet30_SBloInvMassRebin->SetLineColor(kBlue);ggMet30_SBloInvMassRebin->SetFillColor(kBlue);ggMet30_SBloInvMassRebin->SetFillStyle(3004);
  ggMet30_SBhiInvMassRebin->SetMarkerColor(kRed);ggMet30_SBhiInvMassRebin->SetLineColor(kRed);ggMet30_SBhiInvMassRebin->SetFillColor(kRed);ggMet30_SBhiInvMassRebin->SetFillStyle(3004);
  ggMET30_InvMassSBcomb->SetMarkerColor(42);ggMET30_InvMassSBcomb->SetLineColor(42);ggMET30_InvMassSBcomb->SetFillColor(42);//ggMET30_InvMassSBcomb->SetFillStyle(3004);
  ggMet30_SBloInvMassRebin->SetMarkerSize(0.75);ggMet30_SBhiInvMassRebin->SetMarkerSize(0.75);ggMET30_InvMassSBcomb->SetMarkerSize(0.75);
  
  //THStack *HiggsBernStack = new THStack("HiggsBernStack","");
  //HiggsBernStack->Add(TTHggMet30New2);

  THStack *HiggsMetStack30 = new THStack("HiggsMetStack30","");
  HiggsMetStack30->Add(TTHggMetRebin30);HiggsMetStack30->Add(WZHggMetRebin30);HiggsMetStack30->Add(VBFHggMetRebin30);HiggsMetStack30->Add(ggHggMetRebin30);
  TH1F* SMHiggsMet30=(TH1F*)TTHggMetRebin30->Clone();SMHiggsMet30->Add(VBFHggMetRebin30);SMHiggsMet30->Add(WZHggMetRebin30);SMHiggsMet->Add(ggHggMetRebin30);
 
  HiggsMetStack30->Add(ggMET30_InvMassSBcomb);
  TH1F* ggMET30_InvMassSBcomb_plus_SMhiggs = (TH1F*)ggMET30_InvMassSBcomb->Clone();ggMET30_InvMassSBcomb_plus_SMhiggs->Add(SMHiggsMet30);
  ggMET30_InvMassSBcomb_plus_SMhiggs->SetMarkerSize(0);ggMET30_InvMassSBcomb_plus_SMhiggs->SetFillColor(kRed);ggMET30_InvMassSBcomb_plus_SMhiggs->SetFillStyle(3004);

  ggMet30_sigInvMassRebin->GetXaxis()->SetRangeUser(0,249);
  ggMet30_sigInvMassRebin->GetYaxis()->SetRangeUser(1e-5,1000);
  ggMet30_sigInvMassRebin->GetYaxis()->SetTitle("Events / GeV");
  ggMet30_sigInvMassRebin->Draw("PE");
  HiggsMetStack30->Draw("histoSAMES");
  ggMET30_InvMassSBcomb_plus_SMhiggs->Draw("E2SAMES");
  ggMet30_sigInvMassRebin->Draw("PESAMES");
  p1->RedrawAxis();
  p2->cd();
  TH1F* ggMet30_sigInvMassRebin_Div_low = (TH1F*)ggMet30_sigInvMassRebin->Clone();
  TH1F* ggMet30_sigInvMassRebin_Div_high = (TH1F*)ggMet30_sigInvMassRebin->Clone();
  TH1F* ggMet30_sigInvMassRebin_Div_comb = (TH1F*)ggMet30_sigInvMassRebin->Clone();
  TH1F* ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs = (TH1F*)ggMet30_sigInvMassRebin->Clone();
  ggMet30_sigInvMassRebin_Div_low->Divide(ggMet30_SBloInvMassRebin);
  ggMet30_sigInvMassRebin_Div_high->Divide(ggMet30_SBhiInvMassRebin);
  ggMet30_sigInvMassRebin_Div_comb->Divide(ggMET30_InvMassSBcomb);
  ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs->Divide(ggMET30_InvMassSBcomb_plus_SMhiggs);
  ggMet30_sigInvMassRebin_Div_low->SetLineColor(kBlue);ggMet30_sigInvMassRebin_Div_high->SetMarkerColor(kBlue);ggMet30_sigInvMassRebin_Div_high->SetMarkerSize(0.75);
  ggMet30_sigInvMassRebin_Div_comb->SetLineColor(kViolet);ggMet30_sigInvMassRebin_Div_comb->SetMarkerColor(kViolet);ggMet30_sigInvMassRebin_Div_comb->SetMarkerSize(0.75);
  ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs->SetLineColor(kRed);ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs->SetMarkerColor(kRed);//ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs->SetMarkerSize(0.75);
  ggMet30_sigInvMassRebin_Div_high->SetLineColor(kRed);ggMet30_sigInvMassRebin_Div_high->SetMarkerColor(kRed);ggMet30_sigInvMassRebin_Div_high->SetMarkerSize(0.75);
  ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs->GetYaxis()->SetTitle("");
  ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs->SetTitle("");
  ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs->GetYaxis()->SetRangeUser(0.,2.);
  ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs->SetTitle("");
  ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs->GetYaxis()->SetTitle("#frac{Data}{Prediction}");
  ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs->GetYaxis()->SetTitleOffset(0.48);
  ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs->GetYaxis()->SetTitleSize(0.15);
  ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs->GetXaxis()->SetTitleSize(0.2);
  ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs->GetYaxis()->SetLabelSize(0.12);
  ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs->GetXaxis()->SetLabelSize(0.15);
  ggMet30_sigInvMassRebin_Div_comb_plusSMhiggs->Draw("PE");
  //ggMet30_sigInvMassRebin_Div_high->Draw("PEsames");
  //ggMet30_sigInvMassRebin_Div_low->Draw("PEsames");
  l1.DrawLine(0,1,250,1);
  p3->cd();p3->Clear();
  met2->Draw();
  c2->Print("Plots/Higgs/ggMET30_InvarMassCuts_BERNfit.png");

  c1->cd();
  /*
  //W jets
  c1->SetLogz(1);
  TH3F* WjetInvMassVquarkProbs = (TH3F*)fin->Get("ggJetsInvarMass_vs_QuarkProb_LD_Lead_vs_Trail_Loose");WjetInvMassVquarkProbs->Sumw2();
  TH3F* WjetInvMassVquarkProbsNoB = (TH3F*)fin->Get("ggJetsInvarMassNoB_vs_QuarkProb_LD_Lead_vs_Trail_Loose");WjetInvMassVquarkProbsNoB->Sumw2();
  //TH3F* WjetInvMassVquarkProbsWmom = (TH3F*)fin->Get("ggJetsInvarMass_vs_QuarkProb_LD_Lead_vs_Trail_Wmother_Loose");WjetInvMassVquarkProbsWmom->Sumw2();
  //TH3F* WjetInvMassVquarkProbsZmom = (TH3F*)fin->Get("ggJetsInvarMass_vs_QuarkProb_LD_Lead_vs_Trail_Zmother_Loose");WjetInvMassVquarkProbsZmom->Sumw2();

  TH2D* WjetQuarkProb = (TH2D*)WjetInvMassVquarkProbs->Project3D("xy");
  WjetQuarkProb->GetXaxis()->SetTitle("Quark Prob. Lead");WjetQuarkProb->GetYaxis()->SetTitle("Quark Prob. Trail");
  WjetQuarkProb->RebinX(5);WjetQuarkProb->RebinY(5);
  WjetQuarkProb->Draw("COLZ");
  l1.SetLineColor(kRed);l1.SetLineWidth(2);l1.DrawLine(0.75,0,0.75,0.75);l1.DrawLine(0,0.75,0.75,0.75);
  c1->Print("Plots/Higgs/WjetsQuarkProdLeadVsTrail.png");

  TH2D* WjetQuarkProbNoB = (TH2D*)WjetInvMassVquarkProbsNoB->Project3D("xy");
  WjetQuarkProbNoB->GetXaxis()->SetTitle("Quark Prob. Lead");WjetQuarkProbNoB->GetYaxis()->SetTitle("Quark Prob. Trail");
  WjetQuarkProbNoB->RebinX(5);WjetQuarkProbNoB->RebinY(5);
  WjetQuarkProbNoB->Draw("COLZ");
  l1.DrawLine(0.75,0,0.75,0.75);l1.DrawLine(0,0.75,0.75,0.75);
  c1->Print("Plots/Higgs/WjetsNoBQuarkProdLeadVsTrail.png");

  TH1D* WjetInvMass = (TH1D*)WjetInvMassVquarkProbs->ProjectionZ("WjetInvMass",0,999,0,999,"eo");
  WjetInvMass->Rebin(6);
  WjetInvMass->GetXaxis()->SetRangeUser(15,250);
  WjetInvMass->Draw();
  c1->Print("Plots/Higgs/WjetsInvarMass.png");
  
//    TH1D* WjetInvMassWmom = (TH1D*)WjetInvMassVquarkProbsWmom->ProjectionZ("WjetInvMassWmom",0,999,0,999,"eo");
//    WjetInvMassWmom->Rebin(6);
//    WjetInvMassWmom->GetXaxis()->SetRangeUser(15,250);
//    WjetInvMassWmom->Draw();
//    c1->Print("Plots/Higgs/WjetsInvarMassWmom.png");
//
//    TH1D* WjetInvMassZmom = (TH1D*)WjetInvMassVquarkProbsZmom->ProjectionZ("WjetInvMassZmom",0,999,0,999,"eo");
//    WjetInvMassZmom->Rebin(6);
//    WjetInvMassZmom->GetXaxis()->SetRangeUser(15,250);
//    WjetInvMassZmom->Draw();
//    c1->Print("Plots/Higgs/WjetsInvarMassZmom.png");
//  
//    WjetInvMassWmom->SetLineColor(kRed);WjetInvMassWmom->SetMarkerColor(kRed);WjetInvMassWmom->SetFillColor(kRed);
//    WjetInvMassZmom->SetLineColor(kBlue);WjetInvMassZmom->SetMarkerColor(kBlue);WjetInvMassZmom->SetFillColor(kBlue);
//    THStack *WZmoms = new THStack("WZmoms","");WZmoms->Add(WjetInvMassWmom);WZmoms->Add(WjetInvMassZmom);
//    WZmoms->Draw("histo");
//    WZmoms->GetXaxis()->SetRangeUser(60,139);
//    c1->Print("Plots/Higgs/WjetsInvarMassWandZmoms.png");
  
  int bin75=WjetQuarkProbNoB->GetXaxis()->FindBin(0.75);
  TH1D* WjetInvMassHighQprob = (TH1D*)WjetInvMassVquarkProbs->ProjectionZ("WjetInvMassHighQprob",bin75,999,bin75,999,"eo");
  WjetInvMassHighQprob->Rebin(6);
  WjetInvMassHighQprob->GetXaxis()->SetRangeUser(15,250);
  WjetInvMassHighQprob->Draw();
  c1->Print("Plots/Higgs/WjetsInvarMassHighQprob.png");

  TH1D* WjetInvMassNoB = (TH1D*)WjetInvMassVquarkProbsNoB->ProjectionZ("WjetInvMassNoB",0,999,0,999,"eo");
  WjetInvMassNoB->Rebin(6);
  WjetInvMassNoB->GetXaxis()->SetRangeUser(15,250);
  WjetInvMassNoB->Draw();
  c1->Print("Plots/Higgs/WjetsInvarMassNoB.png");

  TH1D* WjetInvMassNoBHighQprob = (TH1D*)WjetInvMassVquarkProbsNoB->ProjectionZ("WjetInvMassNoBHighQprob",bin75,999,bin75,999,"eo");
  WjetInvMassNoBHighQprob->Rebin(6);
  WjetInvMassNoBHighQprob->GetXaxis()->SetRangeUser(15,250);
  WjetInvMassNoBHighQprob->Draw();
  c1->Print("Plots/Higgs/WjetsInvarMassNoBHighQprob.png");

  WjetInvMass->Draw();
  WjetInvMassHighQprob->SetLineColor(kRed);WjetInvMassHighQprob->SetMarkerColor(kRed);
  WjetInvMassNoB->SetLineColor(kBlue);WjetInvMassNoB->SetMarkerColor(kBlue);
  WjetInvMassNoBHighQprob->SetLineColor(kGreen);WjetInvMassNoBHighQprob->SetMarkerColor(kGreen);
  WjetInvMassHighQprob->Draw("SAMES");
  WjetInvMassNoB->Draw("SAMES");
  WjetInvMassNoBHighQprob->Draw("SAMES");
  c1->Print("Plots/Higgs/WjetsInvarMassAll.png");

  float scale = 1./WjetInvMass->Integral();WjetInvMass->Scale(scale);
  scale = 1./WjetInvMassHighQprob->Integral();WjetInvMassHighQprob->Scale(scale);
  scale = 1./WjetInvMassNoB->Integral();WjetInvMassNoB->Scale(scale);
  scale = 1./WjetInvMassNoBHighQprob->Integral();WjetInvMassNoBHighQprob->Scale(scale);
  WjetInvMassNoBHighQprob->Draw();
  WjetInvMass->Draw("SAMES");
  WjetInvMassHighQprob->Draw("SAMES");
  WjetInvMassNoB->Draw("SAMES");
  c1->Print("Plots/Higgs/WjetsInvarMassAllScale.png");
  */


  // all cases
  //W->e nu -- 1 electron, <2 jets
  c1->cd();
  TH2F* WeDataPhoE = (TH2F*)fin->Get("ggMetVsInvarMassEPho_Loose_1Ele_0_1Jets");
  TH1F* WeDataPhoEInvMass = (TH1F*)WeDataPhoE->ProjectionX("WeDataPhoEInvMass",0,-1,"eo");
  WeDataPhoEInvMass->Rebin(2);  WeDataPhoEInvMass->GetXaxis()->SetRangeUser(10,179);
  /*
  for(int i=0;i<WeDataPhoEInvMass->GetNbinsX();i++){
    cout<<"bin low edge: "<<WeDataPhoEInvMass->GetBinLowEdge(i)<<"  content: "<<WeDataPhoEInvMass->GetBinContent(i)<<endl;
    }*/
  c1->SetLogy(0);
  WeDataPhoEInvMass->SetTitle("");WeDataPhoEInvMass->GetXaxis()->SetTitle("M_{e,#gamma} [GeV]");WeDataPhoEInvMass->GetYaxis()->SetRangeUser(0,22);
  WeDataPhoEInvMass->Draw();
  TLine epholine(86,0,86,22);epholine.SetLineColor(kRed);epholine.Draw();epholine.DrawLine(96,0,96,22);
  c1->Print("Plots/Higgs/Exclusive_WeDataEPhoInvMass.png");
  c1->Print("Plots/Higgs/Exclusive_WeDataEPhoInvMass.pdf");
  //WeDataPhoEInvMass->Draw("histo");
  //c1->Print("Plots/Higgs/Exclusive_WeDataEPhoInvMassHisto.png");
  c1->SetLogy(1);
  TH2F* WeData = (TH2F*)fin->Get("ggMetVsInvarMass_Loose_1Ele_0_1Jets");
  TH2F* WeggHgg = (TH2F*)f_ggHgg->Get("ggMetVsInvarMass_Loose_1Ele_0_1Jets");
  TH2F* WeWZHgg = (TH2F*)f_WZHgg->Get("ggMetVsInvarMass_Loose_1Ele_0_1Jets");
  TH2F* WeTTHgg = (TH2F*)f_TTHgg->Get("ggMetVsInvarMass_Loose_1Ele_0_1Jets");
  TH2F* WeVBFHgg = (TH2F*)f_VBFHgg->Get("ggMetVsInvarMass_Loose_1Ele_0_1Jets");
  //TH2F* We_aaW130 = (TH2F*)f_aaW130->Get("ggMetVsInvarMass_Loose_1Ele_0_1Jets");
  //TH2F* We_aaW275 = (TH2F*)f_aaW275->Get("ggMetVsInvarMass_Loose_1Ele_0_1Jets");
  //WeggHgg->Scale((PhoEffScale2*L_int*2.29e-03*19.22)/99989.);//125GeV=/96290);  
  //WeVBFHgg->Scale((PhoEffScale2*L_int*2.29e-03*1.544)/95677.);//125GeV=/99885); 
  //WeTTHgg->Scale((PhoEffScale2*L_int*2.29e-03*.1271)/100048.);//125GeV=/100224);
  //WeWZHgg->Scale((PhoEffScale2*L_int*2.29e-03*(.6782/**(.3257+.014)*/+.3843/**.2*/))/100320);//125 and 126 GeV have same # events
  WeggHgg->SetLineColor(kGreen);WeggHgg->SetMarkerColor(kGreen);WeggHgg->SetFillColor(kGreen);
  WeWZHgg->SetLineColor(kCyan);WeWZHgg->SetMarkerColor(kCyan);WeWZHgg->SetFillColor(kCyan);
  WeTTHgg->SetLineColor(31);WeTTHgg->SetMarkerColor(31);WeTTHgg->SetFillColor(31);
  WeVBFHgg->SetLineColor(kRed+3);WeVBFHgg->SetMarkerColor(kRed+3);WeVBFHgg->SetFillColor(kRed+3);
  //WeggHgg->SetFillStyle(0);WeWZHgg->SetFillStyle(0);WeTTHgg->SetFillStyle(0);WeVBFHgg->SetFillStyle(0);
  WeggHgg->SetLineWidth(2);WeWZHgg->SetLineWidth(2);WeVBFHgg->SetLineWidth(2);WeTTHgg->SetLineWidth(2);
  /*We_aaW130->Scale(scale_aaW130);
  We_aaW130->SetLineColor(kBlue);We_aaW130->SetLineWidth(3);
  We_aaW275->Scale(scale_aaW275);
  We_aaW275->SetLineColor(kRed);We_aaW275->SetLineWidth(3);*/
  TH1D* WeDataProjX = (TH1D*)WeData->ProjectionX("WeDataProjX",0,-1,"eo");
  int bin30Down = WeData->GetXaxis()->FindBin(29.9);
  TH1D* WeDataProjX_METlt30 = (TH1D*)WeData->ProjectionX("WeDataProjX_METlt30",0,bin30Down,"eo");
  TH1D* WeggHggProjX = (TH1D*)WeggHgg->ProjectionX("WeggHggProjX",0,-1,"eo");
  TH1D* WeWZHggProjX = (TH1D*)WeWZHgg->ProjectionX("WeWZHggProjX",0,-1,"eo");
  TH1D* WeTTHggProjX = (TH1D*)WeTTHgg->ProjectionX("WeTTHggProjX",0,-1,"eo");
  TH1D* WeVBFHggProjX = (TH1D*)WeVBFHgg->ProjectionX("WeVBFHggProjX",0,-1,"eo");
  //TH1D* WeAAW130ProjX = (TH1D*)We_aaW130->ProjectionX("WeAAW130ProjX",0,-1,"eo");
  //TH1D* WeAAW275ProjX = (TH1D*)We_aaW275->ProjectionX("WeAAW275ProjX",0,-1,"eo");
  TH1D* WeDataProjY = (TH1D*)WeData->ProjectionY("WeDataProjY",0,-1,"eo");
  TH1D* WeggHggProjY = (TH1D*)WeggHgg->ProjectionY("WeggHggProjY",0,-1,"eo");
  TH1D* WeWZHggProjY = (TH1D*)WeWZHgg->ProjectionY("WeWZHggProjY",0,-1,"eo");
  TH1D* WeTTHggProjY = (TH1D*)WeTTHgg->ProjectionY("WeTTHggProjY",0,-1,"eo");
  TH1D* WeVBFHggProjY = (TH1D*)WeVBFHgg->ProjectionY("WeVBFHggProjY",0,-1,"eo");
  //TH1D* WeAAW130ProjY = (TH1D*)We_aaW130->ProjectionY("WeAAW130ProjY",0,-1,"eo");
  //TH1D* WeAAW275ProjY = (TH1D*)We_aaW275->ProjectionY("WeAAW275ProjY",0,-1,"eo");

  c1->SetLogy(0);/*
  WeAAW130ProjX->SetLineColor(kBlue);WeAAW275ProjX->SetLineColor(kRed);
  WeAAW130ProjX->GetXaxis()->SetRangeUser(110,140);WeAAW130ProjX->Draw("histo");WeAAW275ProjX->Draw("histoSAMES");
  c1->Print("Plots/Higgs/Exclusive_WeInvarMassAAW.png");
		 */
  WeDataProjX->Rebin(2);WeDataProjX_METlt30->Rebin(2);WeggHggProjX->Rebin(2);WeTTHggProjX->Rebin(2);WeWZHggProjX->Rebin(2);WeVBFHggProjX->Rebin(2); 
  int WeBin103 = WeDataProjX->FindBin(102.9),WeBin118 = WeDataProjX->FindBin(117.9),WeBin120 = WeDataProjX->FindBin(119.9),WeBin131 = WeDataProjX->FindBin(130.9),WeBin133 = WeDataProjX->FindBin(132.9),WeBin163 = WeDataProjX->FindBin(162.9);
  cout<<"We events invmass less than 103: "<<WeDataProjX->Integral(0,WeBin103)<<" 103<invmass<163: "<<WeDataProjX->Integral(WeBin103+1,WeBin163)<<" invmass<163: "<<WeDataProjX->Integral(WeBin163+1,-1)<<" 103<invmass<118: "<<WeDataProjX->Integral(WeBin103+1,WeBin118)<<" 120<invmass<131: "<<WeDataProjX->Integral(WeBin120+1,WeBin131)*0<<" 133<invmass<163: "<<WeDataProjX->Integral(WeBin133+1,WeBin163)<<endl;
  //MakeBlindInvMass(WeDataProjX);
  c1->SetLogy(1);
  WeDataProjY->GetXaxis()->SetRangeUser(0,200);WeDataProjY->Draw();
  c1->Print("Plots/Higgs/Exclusive_WeDataMet.png");

  THStack *WeMetStack = new THStack("WeMetStack","");
  WeMetStack->Add(WeTTHggProjY);WeMetStack->Add(WeVBFHggProjY);WeMetStack->Add(WeWZHggProjY);WeMetStack->Add(WeggHggProjY);
  WeMetStack->Draw("histo");WeMetStack->GetXaxis()->SetRangeUser(0,249.9);
  c1->Print("Plots/Higgs/Exclusive_WeSMHiggsMet.png");
  c1->Print("Plots/Higgs/Exclusive_WeSMHiggsMet.pdf");
  TH1F* WeggHggProjYRebin=(TH1F*)WeggHggProjY->Rebin(NmetBinsXtraWide,"WeggHggProjYRebin",xbinsXtraWide);
  TH1F* WeTTHggProjYRebin=(TH1F*)WeTTHggProjY->Rebin(NmetBinsXtraWide,"WeTTHggProjYRebin",xbinsXtraWide);
  TH1F* WeWZHggProjYRebin=(TH1F*)WeWZHggProjY->Rebin(NmetBinsXtraWide,"WeWZHggProjYRebin",xbinsXtraWide);
  TH1F* WeVBFHggProjYRebin=(TH1F*)WeVBFHggProjY->Rebin(NmetBinsXtraWide,"WeVBFHggProjYRebin",xbinsXtraWide);
  //TH1F* WeAAW130ProjYRebin=(TH1F*)WeAAW130ProjY->Rebin(NmetBinsXtraWide,"WeAAW130ProjYRebin",xbinsXtraWide);
  //TH1F* WeAAW275ProjYRebin=(TH1F*)WeAAW275ProjY->Rebin(NmetBinsXtraWide,"WeAAW275ProjYRebin",xbinsXtraWide);

  TH1F* met_Phi_ele_ggH = (TH1F*)f_ggHgg->Get("met_Phi_ggLoose_Ele");
  TH1F* met_Phi_ele_WZH = (TH1F*)f_WZHgg->Get("met_Phi_ggLoose_Ele");
  TH1F* met_Phi_ele_TTH = (TH1F*)f_TTHgg->Get("met_Phi_ggLoose_Ele");
  TH1F* met_Phi_ele_VBFH = (TH1F*)f_VBFHgg->Get("met_Phi_ggLoose_Ele");
  
  AddOverflowToLastBin(WeggHggProjYRebin);
  AddOverflowToLastBin(WeTTHggProjYRebin);
  AddOverflowToLastBin(WeWZHggProjYRebin);
  AddOverflowToLastBin(WeVBFHggProjYRebin);

  ScaleHiggs(WeggHggProjYRebin,met_Phi_ele_ggH,19.22,99989.);
  ScaleHiggs(WeTTHggProjYRebin,met_Phi_ele_TTH,.1271,100048.);
  ScaleHiggs(WeWZHggProjYRebin,met_Phi_ele_WZH,(.6782/**(.3257+.014)*/+.3843/**.2*/),100320.);
  ScaleHiggs(WeVBFHggProjYRebin,met_Phi_ele_VBFH,1.544,95677.);


  double SMHiggsMetIntWe_ggH=0.,SMHiggsMetIntErrWe_ggH=0.,SMHiggsMetIntWe_WZH=0.,SMHiggsMetIntErrWe_WZH=0.,SMHiggsMetIntWe_TTH=0.,SMHiggsMetIntErrWe_TTH=0.,SMHiggsMetIntWe_VBFH=0.,SMHiggsMetIntErrWe_VBFH=0.;
  SMHiggsMetIntWe_ggH=WeggHggProjYRebin->IntegralAndError(0,-1,SMHiggsMetIntErrWe_ggH);
  SMHiggsMetIntWe_WZH=WeWZHggProjYRebin->IntegralAndError(0,-1,SMHiggsMetIntErrWe_WZH);
  SMHiggsMetIntWe_TTH=WeTTHggProjYRebin->IntegralAndError(0,-1,SMHiggsMetIntErrWe_TTH);
  SMHiggsMetIntWe_VBFH=WeVBFHggProjYRebin->IntegralAndError(0,-1,SMHiggsMetIntErrWe_VBFH);
  //float SMHiggsMetIntWe=WeTTHggProjYRebin->Integral()+WeVBFHggProjYRebin->Integral()+WeWZHggProjYRebin->Integral()+WeggHggProjYRebin->Integral();
  float SMHiggsMetIntWe=SMHiggsMetIntWe_ggH+SMHiggsMetIntWe_WZH+SMHiggsMetIntWe_TTH+SMHiggsMetIntWe_VBFH;
  float SMHiggsMetIntErrWe=sqrt(SMHiggsMetIntErrWe_ggH*SMHiggsMetIntErrWe_ggH+SMHiggsMetIntErrWe_WZH*SMHiggsMetIntErrWe_WZH+SMHiggsMetIntErrWe_TTH*SMHiggsMetIntErrWe_TTH+SMHiggsMetIntErrWe_VBFH*SMHiggsMetIntErrWe_VBFH);

  TH1F* SMHiggsMetWeNoDivByWidth = (TH1F*)WeTTHggProjYRebin->Clone();SMHiggsMetWeNoDivByWidth->Add(WeVBFHggProjYRebin);SMHiggsMetWeNoDivByWidth->Add(WeggHggProjYRebin);SMHiggsMetWeNoDivByWidth->Add(WeWZHggProjYRebin);

  DivideBy15gev(WeggHggProjYRebin);
  DivideBy15gev(WeTTHggProjYRebin);
  DivideBy15gev(WeWZHggProjYRebin);
  DivideBy15gev(WeVBFHggProjYRebin);
  
  THStack *WeMetStackRebin = new THStack("WeMetStackRebin","");
  THStack *WeMetStackRebin2 = new THStack("WeMetStackRebin2","");
  WeTTHggProjYRebin->SetLineColor(kBlack);WeVBFHggProjYRebin->SetLineColor(kBlack);WeWZHggProjYRebin->SetLineColor(kBlack);WeggHggProjYRebin->SetLineColor(kBlack);
  WeTTHggProjYRebin->SetLineWidth(1);WeVBFHggProjYRebin->SetLineWidth(1);WeWZHggProjYRebin->SetLineWidth(1);WeggHggProjYRebin->SetLineWidth(1);
  TH1F* SMHiggsAdditionForWe = (TH1F*)WeTTHggProjYRebin->Clone();SMHiggsAdditionForWe->Add(WeVBFHggProjYRebin);SMHiggsAdditionForWe->Add(WeWZHggProjYRebin);SMHiggsAdditionForWe->Add(WeggHggProjYRebin);
  WeMetStackRebin2->Add(WeTTHggProjYRebin);WeMetStackRebin2->Add(WeVBFHggProjYRebin);WeMetStackRebin2->Add(WeWZHggProjYRebin);WeMetStackRebin2->Add(WeggHggProjYRebin);
  WeMetStackRebin->Add(SMHiggsAdditionForWe);
  WeMetStackRebin2->Draw("histo");WeMetStackRebin2->GetXaxis()->SetRangeUser(0,249.9);WeMetStackRebin2->SetMaximum(.4);WeMetStackRebin2->SetMinimum(3.1E-3);
  WeMetStackRebin2->GetYaxis()->SetTitle("Events / GeV");WeMetStackRebin2->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  legHiggs2->Draw();
  c1->Print("Plots/Higgs/Exclusive_WeSMHiggsMetRebin.png");
  c1->Print("Plots/Higgs/Exclusive_WeSMHiggsMetRebin.pdf");
  TH1F* SMHiggsMetWe = (TH1F*)WeTTHggProjYRebin->Clone();SMHiggsMetWe->Add(WeVBFHggProjYRebin);SMHiggsMetWe->Add(WeggHggProjYRebin);SMHiggsMetWe->Add(WeWZHggProjYRebin);
  c1->SetLogy(0);
  THStack *WeInvMassStack = new THStack("WeInvMassStack","");
  WeInvMassStack->Add(WeTTHggProjX);WeInvMassStack->Add(WeVBFHggProjX);WeInvMassStack->Add(WeWZHggProjX);WeInvMassStack->Add(WeggHggProjX);
  WeInvMassStack->Draw("histo");WeInvMassStack->GetXaxis()->SetRangeUser(110,139.9);
  c1->Print("Plots/Higgs/Exclusive_WeSMHiggsInvMass.png");
  cout<<"SM higgs Expected Yield, Wenu : "<<SMHiggsMetIntWe<<" +- "<<SMHiggsMetIntErrWe<<endl;
  WeDataProjX->GetXaxis()->SetRangeUser(20,289);WeDataProjX->SetTitle("");
  WeDataProjX->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");WeDataProjX->GetYaxis()->SetTitle("Events / GeV");
  WeDataProjX->GetXaxis()->SetTitleOffset(1.1);WeDataProjX->GetYaxis()->SetTitleOffset(1.4);
  WeDataProjX->GetXaxis()->SetLabelOffset(0.02);
  WeDataProjX->GetYaxis()->SetLabelOffset(0.02);
  WeDataProjX->Draw();
  PrelimTextInvMass->Draw();
  TPaveText *ggTextWeInvMass = new TPaveText(.29,.78,.6,.89,"NDC");
  ggTextWeInvMass->AddText("#gamma#gamma + e");
  ////ggTextWe->AddText("");
  ggTextWeInvMass->SetTextSize(.08);
  ggTextWeInvMass->SetTextFont(42);
  ggTextWeInvMass->SetFillStyle(0);
  ggTextWeInvMass->SetFillColor(0);
  ggTextWeInvMass->SetBorderSize(0);
  ggTextWeInvMass->Draw();
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMass.png");
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMass.pdf");
  WeDataProjX_METlt30->GetXaxis()->SetRangeUser(20,400);
  WeDataProjX_METlt30->Draw();
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMass_METlt30.png");
  //invariant mass roofit
  RooRealVar xMetWe_bg("xMetWe_bg","m_{#gamma#gamma}",95,200,"GeV");
  xMetWe_bg.setRange("sb_loWe_bg",sbLoLo,sbLoHi);xMetWe_bg.setRange("sb_hiWe_bg",sbHiLo,sbHiHi);xMetWe_bg.setRange("fullWe_bg",sbLoLo,sbHiHi);xMetWe_bg.setRange("sigWe_bg",sigLo,sigHi);
  RooDataHist dataMetWe_bg("dataMetWe_bg","dataset_bg",xMetWe_bg,WeDataProjX);
  RooRealVar lambdaWe("lambdaWe","lambdaWe",1.,-10.,5.);
  RooRealVar PolWe1_bg("PolWe1_bg","PolWe1_bg",0,-.09,.09);
  RooRealVar PolWe2_bg("PolWe2_bg","PolWe2_bg",0,-.00001,.00001);
  RooRealVar PolWe3_bg("PolWe3_bg","PolWe3_bg",0,-.000001,.000001);
  RooRealVar PolWe4_bg("PolWe4_bg","PolWe4_bg",0,-.0000001,.0000001);
  RooRealVar PolWe5_bg("PolWe5_bg","PolWe5_bg",0,-.000000001,.000000001);
  //RooPolynomial PolWe_bg("PolWe_bg","5th order Polynomial",xMetWe_bg,RooArgList(PolWe1_bg,PolWe2_bg,PolWe3_bg,PolWe4_bg,PolWe5_bg));
  RooExponential PolWe_bg("ExpoWe_bg","ExpoWe_bg",xMetWe_bg,lambdaWe);
  RooRealVar WeYield_bg("WeYield_bg","WeYield_bg",70000,0,170000);
  //RooAddPdf MetWePdf_bg("MetWePdf_bg","MetWePdf_bg",RooArgList(PolWe_bg),RooArgList(WeYield_bg));
  RooAddPdf MetWePdf_bg("MetWePdf_bg","MetWePdf_bg",RooArgList(PolWe_bg),RooArgList(WeYield_bg));
  //RooFitResult *rNoMetCutWe_bg = MetWePdf_bg.chi2FitTo(dataMetWe_bg/*,Extended(kTRUE)*/,Save(),Range("sb_lo,sb_hi"));
  RooFitResult *rNoMetCutWe_bg = MetWePdf_bg.fitTo(dataMetWe_bg,Extended(kTRUE),Save(),Range(/*"fullWe_bg"*/"sb_loWe_bg,sb_hiWe_bg"));
  RooArgList pars_We(*MetWePdf_bg.getParameters(RooArgSet(xMetWe_bg) ) );
  RooArgSet prodSet_We(MetWePdf_bg); //prodSet.add(nsig);
  RooProduct unNormPdf_We("fitted Function", "fitted Function", prodSet_We);
  TF1 * fit_We = unNormPdf_We.asTF(RooArgList(xMetWe_bg), pars_We, RooArgList(xMetWe_bg));
  float nsig_We = ((RooRealVar*) pars_We.find("WeYield_bg"))->getVal();
  Double_t integ_full_We = fit_We->Integral(sbLoLo, sbHiHi);
  Double_t integ_We = nsig_We*fit_We->Integral(sigLo, sigHi, 0)/integ_full_We;
  Double_t dinteg_We = nsig_We*fit_We->IntegralError(sigLo, sigHi, 0, rNoMetCutWe_bg->covarianceMatrix().GetMatrixArray())/integ_full_We;
  RooPlot *xframeMetWe_bg = xMetWe_bg.frame(Title("5th Order Polynomial Background, gg We with no MET cut"));
  MetWePdf_bg.paramOn(xframeMetWe_bg, Format("NE",AutoPrecision(2)),Layout(0.6,0.995,0.9) );
  dataMetWe_bg.plotOn(xframeMetWe_bg,LineColor(kBlack),MarkerStyle(20),MarkerSize(0.3));
  //MetWePdf_bg.plotOn(xframeMetWe_bg,Components(PolWe_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCutWe_bg,3,kTRUE),FillColor(kViolet));
  MetWePdf_bg.plotOn(xframeMetWe_bg,Components(PolWe_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCutWe_bg,2,kTRUE),FillColor(kGreen));
  MetWePdf_bg.plotOn(xframeMetWe_bg,Components(PolWe_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCutWe_bg,1,kTRUE),FillColor(kOrange));
  MetWePdf_bg.plotOn(xframeMetWe_bg,Components(PolWe_bg),LineColor(kRed),LineStyle(kDashed),Range(102,180));
  MetWePdf_bg.plotOn(xframeMetWe_bg,Components(PolWe_bg),LineColor(kRed));
  dataMetWe_bg.plotOn(xframeMetWe_bg,LineColor(kBlack),MarkerStyle(20),MarkerSize(0.3));
  //xframeMetWe_bg->SetAxisRange(106,149,"X");
  //xframeMetWe_bg->SetAxisRange(0,1900,"Y");
  //xframeMetWe_bg->SetMaximum(2000);
  xframeMetWe_bg->Draw();
  //line125->DrawLine(125.3,0,125.3,xframeMetWe_bg->GetMaximum());
  WeInvMassStack->Draw("histoSAME");
  xframeMetWe_bg->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit.png");

  TH1F* WeDataProjXpow = (TH1F*)WeDataProjX->Clone();
  TH1F* WeDataProjXlin = (TH1F*)WeDataProjX->Clone();
  TH1F* WeDataProjXexpo = (TH1F*)WeDataProjX->Clone();
  TH1F* WeDataProjXpol = (TH1F*)WeDataProjX->Clone();
  TH1F* WeDataProjX_METlt30pow = (TH1F*)WeDataProjX_METlt30->Clone();

  reject=true;
  TF1* fitCurveWe = new TF1("fitCurveWe",fpow,sbFitLoLo,sbFitHiHi,2);
  avg_l = WeDataProjXpow->Integral(WeDataProjXpow->FindBin(sbFitLoLo),WeDataProjXpow->FindBin(sigLo))/float(sbLoHi-sbFitLoLo);avg_u = WeDataProjXpow->Integral(WeDataProjXpow->FindBin(sbHiLo),WeDataProjXpow->FindBin(sbFitHiHi)-1)/float(sbFitHiHi-sbHiLo);avgX_l=(sbLoHi-sbFitLoLo)/2.;avgX_u=(sbFitHiHi-sbHiLo)/2.;
  cout<<avg_l<<"  "<<avg_u<<"  "<<avgX_l<<"  "<<avgX_u<<endl;
  Double_t param1We= (log(avg_l) - log(avg_u))/(log(avgX_l) - log(avgX_u));
  Double_t param0We= avg_l/pow(avgX_l, param1We);
  cout<<"param0We: "<<param0We<<"  param1We: "<<param1We<<endl;
  fitCurveWe->SetParameter(0,param0We);
  fitCurveWe->SetParameter(1,param1We);
  //first just to check the status
  status = WeDataProjXpow->Fit(fitCurveWe,"L0","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResultWe = WeDataProjXpow->Fit(fitCurveWe,"SLLMEV0","",sbFitLoLo,sbFitHiHi);
  TMatrixDSym covWe = fitResultWe->GetCovarianceMatrix();
  fitResultWe->Print("V");
  WeDataProjXpow->GetXaxis()->SetRangeUser(100,164.9);
  reject=false;
  YieldBinWidth=WeDataProjXpow->GetBinWidth(1);
  Double_t WePowYield = fitCurveWe->Integral(sbLoLo,sbHiHi)/YieldBinWidth;
  Double_t WePowYieldErr = fitCurveWe->IntegralError(sbLoLo,sbHiHi,fitResultWe->GetParams(),covWe.GetMatrixArray() )/YieldBinWidth; 
  Double_t WePowYieldSig = fitCurveWe->Integral(sigLo,sigHi)/YieldBinWidth;
  Double_t WePowYieldSigErr = fitCurveWe->IntegralError(sigLo,sigHi,fitResultWe->GetParams(),covWe.GetMatrixArray() )/YieldBinWidth; 
  Double_t WePowSBloYield = WeDataProjXpow->Integral(WeDataProjXpow->FindBin(sbLoLo),WeDataProjXpow->FindBin(sbLoHi-.1))/YieldBinWidth;
  Double_t WePowSBloYieldFit = fitCurveWe->Integral(sbLoLo,sbLoHi)/YieldBinWidth;
  Double_t WePowSBloYieldFitErr = fitCurveWe->IntegralError(sbLoLo,sbLoHi,fitResultWe->GetParams(),covWe.GetMatrixArray() )/YieldBinWidth; 
  cout<<"We low sideband yield from histo: "<<WePowSBloYield<<"  and from fit: "<<WePowSBloYieldFit<<" +- "<<WePowSBloYieldFitErr<<endl; 
  cout<<"We higgs window yield from fit: "<<WePowYieldSig<<" +- "<<WePowYieldSigErr<<endl; 
  Double_t WePowSBhiYield = WeDataProjXpow->Integral(WeDataProjXpow->FindBin(sbHiLo),WeDataProjXpow->FindBin(sbHiHi-.1))/YieldBinWidth;
  Double_t WePowSBhiYieldFit = fitCurveWe->Integral(sbHiLo,sbHiHi)/YieldBinWidth;
  Double_t WePowSBhiYieldFitErr = fitCurveWe->IntegralError(sbHiLo,sbHiHi,fitResultWe->GetParams(),covWe.GetMatrixArray() )/YieldBinWidth; 
  cout<<"We high sideband yield from histo: "<<WePowSBhiYield<<"  and from fit: "<<WePowSBhiYieldFit<<" +- "<<WePowSBhiYieldFitErr<<endl;
  reject=true;
  char strWe[50],strWeSig[50];
  sprintf(strWe,"Fit Yield: %4.2f #pm %4.2f",WePowYield,WePowYieldErr);
  sprintf(strWeSig,"Fit Signal Yield: %4.2f #pm %4.2f",WePowYieldSig,WePowYieldSigErr);
  char strWeRF[50];
  sprintf(strWeRF,"RooFit signal Yield: %4.2f #pm %4.2f",integ_We,dinteg_We);
  TPaveText *textWe= new TPaveText(.5,.65,.84,.84,"NDC");textWe->SetFillStyle(0);textWe->SetBorderSize(0);
  //textWe->AddText(strWe);
  textWe->AddText(strWeSig);textWe->SetTextFont(42);
  //textWe->AddText(strWeRF);
  WeDataProjXpow->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  WeDataProjXpow->GetYaxis()->SetTitle("Events / GeV");
  WeDataProjXpow->GetYaxis()->SetTitleOffset(0.8);
  WeDataProjXpow->GetXaxis()->SetTitleOffset(1.2);
  WeDataProjXpow->SetTitleFont(42,"xy");
  WeDataProjXpow->SetLabelFont(42,"xy");
  WeDataProjXpow->SetTitle("");
  WeDataProjXpow->GetListOfFunctions()->Remove(fitCurveWe);
  WeDataProjXpow->GetYaxis()->SetRangeUser(0,6);WeDataProjXpow->GetYaxis()->SetNdivisions(206);
  WeDataProjXpow->SetMarkerSize(1.5);
  WeDataProjXpow->Draw();
  h_SMS_HH_2W2g_1Ele_Excluded_InvMass->Draw("histosames");
  TF1 *fitCurveWeNew = new TF1("fitCurveWeNew","[0]*pow(x,[1])",sbLoLo,sbHiHi);fitCurveWeNew->SetParameters(fitCurveWe->GetParameters());
  fitCurveWeNew->SetLineColor(kBlue);fitCurveWeNew->SetLineWidth(3);fitCurveWeNew->Draw("SAMES");
  /*
  TF1 *fitCurveWe2 = new TF1("fitCurveWe2","[0]*pow(x,[1])",sbLoHi,sbHiLo);fitCurveWe2->SetParameters(fitCurveWe->GetParameters());
  //fitCurveWe2->SetLineColor(kBlue);fitCurveWe2->Draw("SAMES");
  TF1 *fitCurveWeLeft = new TF1("fitCurveWeLeft","[0]*pow(x,[1])",sbLoLo,sbLoHi);fitCurveWeLeft->SetParameters(fitCurveWe->GetParameters());
  //WeDataProjXpow->GetListOfFunctions()->Add(fitCurveWeLeft);//gROOT->GetListOfFunctions()->Remove(fitCurveWeLeft);
  //fitCurveWeLeft->SetLineColor(kRed);fitCurveWeLeft->Draw("SAMES");
  TF1 *fitCurveWeRight = new TF1("fitCurveWeRight","[0]*pow(x,[1])",sbHiLo,sbHiHi);fitCurveWeRight->SetParameters(fitCurveWe->GetParameters());
  //WeDataProjXpow->GetListOfFunctions()->Add(fitCurveWeRight);//gROOT->GetListOfFunctions()->Remove(fitCurveWeRight);
  //fitCurveWeRight->SetLineColor(kRed);fitCurveWeRight->Draw("SAMES");
  */
  WeDataProjXpow->Draw("SAMES");
  //textWe->Draw("SAMES");
  ggTextWeInvMass->Draw();
  TLegend *WeInvMassLeg = new TLegend(.56,.59,.85,.89,"","brNDC");
  WeInvMassLeg->SetTextSize(.06);WeInvMassLeg->SetTextFont(42);
  WeInvMassLeg->SetFillColor(kWhite);WeInvMassLeg->SetFillStyle(0);WeInvMassLeg->SetBorderSize(0);
  WeInvMassLeg->AddEntry(WeDataProjXpow,"Data","ep");
  WeInvMassLeg->AddEntry(fitCurveWeNew,"Sideband Fit","l");
  //WeInvMassLeg->AddEntry(fitCurveWe2,"Fit Extrapolation","l");
  WeInvMassLeg->AddEntry(h_SMS_HH_2W2g_1Ele_Excluded_InvMass,"Signal hh","f");
  WeInvMassLeg->AddEntry(h_SMS_HH_2W2g_1Ele_Excluded_InvMass,"m_{#tilde{#chi}_{1}^{0}}=130 GeV","");
  WeInvMassLeg->Draw();
  PrelimTextInvMass->Draw();
  TLine vertErr(136.27,5.4,136.27,5.85);vertErr.SetLineWidth(2);
  vertErr.Draw();
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_powLaw.png");
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_powLaw.pdf");

  //TRandom rand;
  TH1F* InvMassFits = new TH1F("InvMassFits","",80,100,180);
  for(int i=0;i<1000;i++){
    InvMassFits->Fill(fitCurveWe->GetRandom(100,180));
  }
  InvMassFits->Draw();
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_powLaw_randomFill.png");
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_powLaw_randomFill.pdf");

  WeDataProjX->GetXaxis()->SetRangeUser(0,1002.);
  TH1F* h_ggEle_fitIntegral = new TH1F("h_ggEle_fitIntegral","h_ggEle_fitIntegral",200,0,20);
  vector<TH1F*> ggEle_InvMass_toys;
  vector<TGraphAsymmErrors*> ggEle_InvMass_toys_tgraph;
  vector<TGraphAsymmErrors*> ggEle_InvMass_toys_tgraph2;
  int numToy=1;
  TRandom3 rr;
  for(int i=0;i<numToy;i++){
    TString name = "name";name+=i;
    ggEle_InvMass_toys.push_back(new TH1F(name,name,WeDataProjX->GetNbinsX(),WeDataProjX->GetXaxis()->GetXmin(),WeDataProjX->GetXaxis()->GetXmax()));
    ggEle_InvMass_toys[i]=(TH1F*)WeDataProjX->Clone();
    ggEle_InvMass_toys_tgraph.push_back(new TGraphAsymmErrors(ggEle_InvMass_toys[i]));
    
    for(int j=0;j<ggEle_InvMass_toys[i]->GetNbinsX();j++){
      double N = ggEle_InvMass_toys_tgraph[i]->GetY()[j];
      double L=0.,U=0.;
      L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
      U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1);
      double newVal = rr.Poisson(rr.Gaus(N,U));
      ggEle_InvMass_toys[i]->SetBinContent(j+1,newVal);
    }
    ggEle_InvMass_toys_tgraph2.push_back(new TGraphAsymmErrors(ggEle_InvMass_toys[i]));

    for(int j=0;j<ggEle_InvMass_toys[i]->GetNbinsX();j++){
      double N = ggEle_InvMass_toys_tgraph2[i]->GetY()[j];
      double L=0.,U=0.;
      L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
      U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1);		     
      ggEle_InvMass_toys_tgraph2[i]->SetPointEYlow(j, N-L);
      ggEle_InvMass_toys_tgraph2[i]->SetPointEYhigh(j, U-N);
    }

    reject=true;
    TF1* fitCurveWeToy = new TF1("fitCurveWeToy",fpow,sbFitLoLo,sbFitHiHi,2);
    avg_l = ggEle_InvMass_toys_tgraph2[i]->Integral(ggEle_InvMass_toys[i]->FindBin(sbFitLoLo),ggEle_InvMass_toys[i]->FindBin(sigLo))/float(sbLoHi-sbFitLoLo);avg_u = ggEle_InvMass_toys[i]->Integral(ggEle_InvMass_toys[i]->FindBin(sbHiLo),ggEle_InvMass_toys[i]->FindBin(sbFitHiHi)-1)/float(sbFitHiHi-sbHiLo);avgX_l=(sbLoHi-sbFitLoLo)/2.;avgX_u=(sbFitHiHi-sbHiLo)/2.;
    cout<<"averages: i="<<i<<"  "<<avg_l<<"  "<<avg_u<<"  "<<avgX_l<<"  "<<avgX_u<<endl;
    Double_t param1WeToy= (log(avg_l) - log(avg_u))/(log(avgX_l) - log(avgX_u));
    Double_t param0WeToy= avg_l/pow(avgX_l, param1WeToy);
    cout<<"param0WeToy: "<<param0WeToy<<"  param1WeToy: "<<param1WeToy<<endl;
    fitCurveWeToy->SetParameter(0,param0WeToy);
    fitCurveWeToy->SetParameter(1,param1WeToy);
    //first just to check the status
    status = ggEle_InvMass_toys_tgraph2[i]->Fit(fitCurveWeToy,"L0","",sbFitLoLo,sbFitHiHi);
    printf("fit status: %i \n",status);
    //Then to get the result
    TFitResultPtr fitResultWeToy = ggEle_InvMass_toys_tgraph2[i]->Fit(fitCurveWeToy,"SLLMEV0","",sbFitLoLo,sbFitHiHi);
    ggEle_InvMass_toys_tgraph2[i]->Draw();
    c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_powLaw_toys_fit.png");

    TMatrixDSym covWeToy = fitResultWeToy->GetCovarianceMatrix();
    fitResultWeToy->Print("V");
    ggEle_InvMass_toys_tgraph2[i]->GetXaxis()->SetRangeUser(100,164.9);
    reject=false;
    YieldBinWidth=ggEle_InvMass_toys[i]->GetBinWidth(1);
    Double_t WeToyPowYield = fitCurveWeToy->Integral(sbLoLo,sbHiHi)/YieldBinWidth;
    Double_t WeToyPowYieldErr = fitCurveWeToy->IntegralError(sbLoLo,sbHiHi,fitResultWeToy->GetParams(),covWeToy.GetMatrixArray() )/YieldBinWidth; 
    Double_t WeToyPowYieldSig = ggEle_InvMass_toys_tgraph2[i]->Integral(ggEle_InvMass_toys[i]->FindBin(sigLo),ggEle_InvMass_toys[i]->FindBin(sigHi-.1))/YieldBinWidth;
    Double_t WeToyPowYieldSigFit = fitCurveWeToy->Integral(sigLo,sigHi)/YieldBinWidth;
    Double_t WeToyPowYieldSigFitErr = fitCurveWeToy->IntegralError(sigLo,sigHi,fitResultWeToy->GetParams(),covWeToy.GetMatrixArray() )/YieldBinWidth; 
    Double_t WeToyPowSBloYield = ggEle_InvMass_toys_tgraph2[i]->Integral(ggEle_InvMass_toys[i]->FindBin(sbLoLo),ggEle_InvMass_toys[i]->FindBin(sbLoHi-.1))/YieldBinWidth;
    Double_t WeToyPowSBloYieldFit = fitCurveWeToy->Integral(sbLoLo,sbLoHi)/YieldBinWidth;
    Double_t WeToyPowSBloYieldFitErr = fitCurveWeToy->IntegralError(sbLoLo,sbLoHi,fitResultWeToy->GetParams(),covWeToy.GetMatrixArray() )/YieldBinWidth; 
    cout<<"WeToy low sideband yield from histo: "<<WeToyPowSBloYield<<"  and from fit: "<<WeToyPowSBloYieldFit<<" +- "<<WeToyPowSBloYieldFitErr<<endl; 
    cout<<"WeToy higgs window yield from histo: "<<WeToyPowYieldSig <<"  and from fit: "<<WeToyPowYieldSigFit<<" +- "<<WeToyPowYieldSigFitErr<<endl; 
    Double_t WeToyPowSBhiYield = ggEle_InvMass_toys_tgraph2[i]->Integral(ggEle_InvMass_toys[i]->FindBin(sbHiLo),ggEle_InvMass_toys[i]->FindBin(sbHiHi-.1))/YieldBinWidth;
    Double_t WeToyPowSBhiYieldFit = fitCurveWeToy->Integral(sbHiLo,sbHiHi)/YieldBinWidth;
    Double_t WeToyPowSBhiYieldFitErr = fitCurveWeToy->IntegralError(sbHiLo,sbHiHi,fitResultWeToy->GetParams(),covWeToy.GetMatrixArray() )/YieldBinWidth; 
    cout<<"WeToy high sideband yield from histo: "<<WeToyPowSBhiYield<<"  and from fit: "<<WeToyPowSBhiYieldFit<<" +- "<<WeToyPowSBhiYieldFitErr<<endl;
    reject=true;
    h_ggEle_fitIntegral->Fill(WeToyPowYieldSigFit);
  }
  h_ggEle_fitIntegral->Draw("histo");
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_powLaw_toys_integral.png");
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_powLaw_toys_integral.pdf");
  int toyColor=1;
  for(int i=0;i<numToy;i++){
    if(i!=0)toyColor+=3;
    ggEle_InvMass_toys_tgraph2[i]->SetLineColor(toyColor);ggEle_InvMass_toys_tgraph2[i]->SetMarkerColor(toyColor);
    if(i==0)ggEle_InvMass_toys_tgraph2[i]->Draw("PE");
    else ggEle_InvMass_toys_tgraph2[i]->Draw("PESAMES");
  }
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_powLaw_toys.png");
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_powLaw_toys.pdf");
  

  reject=true;
  TF1* fitCurveWE_METlt30 = new TF1("fitCurveWE_METlt30",fpow,sbFitLoLo,sbFitHiHi,2);
  avg_l = WeDataProjX_METlt30pow->Integral(WeDataProjX_METlt30pow->FindBin(sbFitLoLo),WeDataProjX_METlt30pow->FindBin(sigLo))/float(sbLoHi-sbFitLoLo);avg_u = WeDataProjX_METlt30pow->Integral(WeDataProjX_METlt30pow->FindBin(sbHiLo),WeDataProjX_METlt30pow->FindBin(sbFitHiHi)-1)/float(sbFitHiHi-sbHiLo);avgX_l=(sbLoHi-sbFitLoLo)/2.;avgX_u=(sbFitHiHi-sbHiLo)/2.;
  Double_t param1WE_METlt30= (log(avg_l) - log(avg_u))/(log(avgX_l) - log(avgX_u));
  Double_t param0WE_METlt30= avg_l/pow(avgX_l, param1WE_METlt30);
  fitCurveWE_METlt30->SetParameter(0,param0WE_METlt30);
  fitCurveWE_METlt30->SetParameter(1,param1WE_METlt30);
  //first just to check the status
  status = WeDataProjX_METlt30pow->Fit(fitCurveWE_METlt30,"L","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResultWE_METlt30 = WeDataProjX_METlt30pow->Fit(fitCurveWE_METlt30,"SLLMEV","",sbFitLoLo,sbFitHiHi);
  //TFitResultPtr fitResultWE_METlt30 = WeDataProjX_METlt30pow->Fit(fitCurveWE_METlt30,"SLLMEV","",110,165);
  TMatrixDSym covWE_METlt30 = fitResultWE_METlt30->GetCovarianceMatrix();
  fitResultWE_METlt30->Print("V");
  WeDataProjX_METlt30pow->GetXaxis()->SetRangeUser(45,199);
  reject=false;
  YieldBinWidth=WeDataProjX_METlt30pow->GetBinWidth(1);
  Double_t WePowYieldSig_METlt30 = fitCurveWE_METlt30->Integral(sigLo,sigHi)/YieldBinWidth;
  Double_t WePowYieldSBlo_METlt30 = fitCurveWE_METlt30->Integral(sbLoLo,sbLoHi)/YieldBinWidth;
  Double_t WePowYieldSBhi_METlt30 = fitCurveWE_METlt30->Integral(sbHiLo,sbHiHi)/YieldBinWidth;
  Double_t WePowYieldSigErr_METlt30 = fitCurveWE_METlt30->IntegralError(sigLo,sigHi,fitResultWE_METlt30->GetParams(),covWE_METlt30.GetMatrixArray() )/YieldBinWidth; 
  reject=true;
  Double_t WeYieldErr_METlt30=0.;
  Double_t WeYield_METlt30 = WeDataProjX_METlt30pow->IntegralAndError(WeDataProjX_METlt30pow->FindBin(sigLo),WeDataProjX_METlt30pow->FindBin(sigHi)-1,WeYieldErr_METlt30);
  Double_t WeYield_METlt30_sbLo = WeDataProjX_METlt30pow->Integral(WeDataProjX_METlt30pow->FindBin(sbLoLo),WeDataProjX_METlt30pow->FindBin(sbLoHi)-1);
  Double_t WeYield_METlt30_sbHi = WeDataProjX_METlt30pow->Integral(WeDataProjX_METlt30pow->FindBin(sbHiLo),WeDataProjX_METlt30pow->FindBin(sbHiHi)-1);
  cout<<"lolo: "<<WeDataProjX_METlt30pow->GetBinLowEdge(WeDataProjX_METlt30pow->FindBin(sbLoLo))<<"  lohi: "<<WeDataProjX_METlt30pow->GetBinLowEdge(WeDataProjX_METlt30pow->FindBin(sbLoHi))<<"  hilo: "<<WeDataProjX_METlt30pow->GetBinLowEdge(WeDataProjX_METlt30pow->FindBin(sbHiLo))<<"  hihi: "<<WeDataProjX_METlt30pow->GetBinLowEdge(WeDataProjX_METlt30pow->FindBin(sbHiHi))<<endl;
  sprintf(strWe,"Data Signal Yield: %4.2f",WeYield_METlt30);
  sprintf(strWeSig,"Fit Signal Yield: %4.2f #pm %4.2f",WePowYieldSig_METlt30,WePowYieldSigErr_METlt30);
  char strSbloWe[50],strSbloFitWe[50],strSbhiWe[50],strSbhiFitWe[50];
  sprintf(strSbloWe,"Data lowSB Yield: %4.2f",WeYield_METlt30_sbLo);
  sprintf(strSbloFitWe,"Fit lowSB Yield: %4.2f",WePowYieldSBlo_METlt30);
  sprintf(strSbhiWe,"Data highSB Yield: %4.2f",WeYield_METlt30_sbHi);
  sprintf(strSbhiFitWe,"Fit highSB Yield: %4.2f",WePowYieldSBhi_METlt30);
  textWe->Clear();
  textWe->AddText(strWe);
  textWe->AddText(strWeSig);
  textWe->AddText(strSbloWe);
  textWe->AddText(strSbloFitWe);
  textWe->AddText(strSbhiWe);
  textWe->AddText(strSbhiFitWe);
  WeDataProjX_METlt30pow->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  WeDataProjX_METlt30pow->SetTitle("");
  WeDataProjX_METlt30pow->Draw();
  textWe->Draw();
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_powLaw_METlt30.png");

  reject=true;
  TF1* fitCurveWeFix = new TF1("fitCurveWeFix",fpow,sbFitLoLo,sbFitHiHi,2);
  avg_l = WeDataProjXpow->Integral(WeDataProjXpow->FindBin(sbFitLoLo),WeDataProjXpow->FindBin(sbLoHi))/20.;avg_u = WeDataProjXpow->Integral(WeDataProjXpow->FindBin(sbHiLo),WeDataProjXpow->FindBin(sbFitHiHi))/69.5;avgX_l=110.;avgX_u=140.25;
  cout<<avg_l<<"  "<<avg_u<<"  "<<avgX_l<<"  "<<avgX_u<<endl;
  Double_t param1WeFix= (log(avg_l) - log(avg_u))/(log(avgX_l) - log(avgX_u));
  Double_t param0WeFix= avg_l/pow(avgX_l, param1WeFix);
  cout<<"param0WeFix: "<<param0WeFix<<"  param1WeFix: "<<param1WeFix<<endl;
  //fitCurveWeFix->SetParameter(0,param0WeFix);
  fitCurveWeFix->SetParameter(0,1.295e+08);
  fitCurveWeFix->FixParameter(1,-4.393);
  //first just to check the status
  status = WeDataProjXpow->Fit(fitCurveWeFix,"L","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResultWeFix = WeDataProjXpow->Fit(fitCurveWeFix,"SL","",sbFitLoLo,sbFitHiHi);
  TMatrixDSym covWeFix = fitResultWeFix->GetCovarianceMatrix();
  fitResultWeFix->Print("V");
  WeDataProjXpow->GetXaxis()->SetRangeUser(100,164.9);
  reject=false;
  YieldBinWidth=WeDataProjXpow->GetBinWidth(1);
  Double_t WeFixPowYield = fitCurveWeFix->Integral(sbLoLo,sbHiHi)/YieldBinWidth;
  Double_t WeFixPowYieldErr = fitCurveWeFix->IntegralError(sbLoLo,sbHiHi,fitResultWeFix->GetParams(),covWeFix.GetMatrixArray() )/YieldBinWidth; 
  Double_t WeFixPowYieldSig = fitCurveWeFix->Integral(sigLo,sigHi)/YieldBinWidth;
  Double_t WeFixPowYieldSigErr = fitCurveWeFix->IntegralError(sigLo,sigHi,fitResultWeFix->GetParams(),covWeFix.GetMatrixArray() )/YieldBinWidth; 
  reject=true;
  strWe[50];strWeSig[50];
  sprintf(strWe,"Fit Yield: %4.2f #pm %4.2f",WeFixPowYield,WeFixPowYieldErr);
  sprintf(strWeSig,"Fit Signal Yield: %4.2f #pm %4.2f",WeFixPowYieldSig,WeFixPowYieldSigErr);
  TPaveText *textWeFix= new TPaveText(.5,.65,.84,.84,"NDC");textWeFix->SetFillStyle(0);textWeFix->SetBorderSize(0);
  //textWeFix->AddText(strWe);
  textWeFix->AddText(strWeSig);textWeFix->SetTextFont(42);
  WeDataProjXpow->Draw();
  textWeFix->Draw("SAMES");
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_powLawFixed.png");

  reject=true;
  TF1* fitCurveWePol5Fix = new TF1("fitCurveWePol5Fix",fpol5,sbFitLoLo,sbFitHiHi,2);
  //avg_l = WeDataProjXpol->Integral(WeDataProjXpol->FindBin(sbFitLoLo),WeDataProjXpol->FindBin(sbLoHi))/20.;avg_u = WeDataProjXpol->Integral(WeDataProjXpol->FindBin(sbHiLo),WeDataProjXpol->FindBin(sbFitHiHi))/69.5;avgX_l=110.;avgX_u=140.25;
  //cout<<avg_l<<"  "<<avg_u<<"  "<<avgX_l<<"  "<<avgX_u<<endl;
  //Double_t param1WePol5Fix= (log(avg_l) - log(avg_u))/(log(avgX_l) - log(avgX_u));
  //Double_t param0WePol5Fix= avg_l/pow(avgX_l, param1WePol5Fix);
  //cout<<"param0WePol5Fix: "<<param0WePol5Fix<<"  param1WePol5Fix: "<<param1WePol5Fix<<endl;
  fitCurveWePol5Fix->SetParameter(0,param0WePol5Fix);
  fitCurveWePol5Fix->FixParameter(1,param1WePol5Fix);
  fitCurveWePol5Fix->FixParameter(2,param2WePol5Fix);
  fitCurveWePol5Fix->FixParameter(3,param3WePol5Fix);
  fitCurveWePol5Fix->FixParameter(4,param4WePol5Fix);
  //first just to check the status
  status = WeDataProjXpol->Fit(fitCurveWePol5Fix,"L","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResultWePol5Fix = WeDataProjXpol->Fit(fitCurveWePol5Fix,"SL","",sbFitLoLo,sbFitHiHi);
  TMatrixDSym covWePol5Fix = fitResultWePol5Fix->GetCovarianceMatrix();
  fitResultWePol5Fix->Print("V");
  WeDataProjXpol->GetXaxis()->SetRangeUser(100,164.9);
  reject=false;
  Double_t WePol5FixPolYield = fitCurveWePol5Fix->Integral(sbLoLo,sbHiHi);
  Double_t WePol5FixPolYieldErr = fitCurveWePol5Fix->IntegralError(sbLoLo,sbHiHi,fitResultWePol5Fix->GetParams(),covWePol5Fix.GetMatrixArray() ); 
  Double_t WePol5FixPolYieldSig = fitCurveWePol5Fix->Integral(sigLo,sigHi);
  Double_t WePol5FixPolYieldSigErr = fitCurveWePol5Fix->IntegralError(sigLo,sigHi,fitResultWePol5Fix->GetParams(),covWePol5Fix.GetMatrixArray() ); 
  reject=true;
  strWe[50],strWeSig[50];
  sprintf(strWe,"Fit Yield: %4.2f #pm %4.2f",WePol5FixPolYield,WePol5FixPolYieldErr);
  sprintf(strWeSig,"Fit Signal Yield: %4.2f #pm %4.2f",WePol5FixPolYieldSig,WePol5FixPolYieldSigErr);
  TPaveText *textWePol5Fix= new TPaveText(.5,.65,.84,.84,"NDC");textWePol5Fix->SetFillStyle(0);textWePol5Fix->SetBorderSize(0);
  //textWePol5Fix->AddText(strWe);
  textWePol5Fix->AddText(strWeSig);textWePol5Fix->SetTextFont(42);
  WeDataProjXpol->Draw();
  textWePol5Fix->Draw("SAMES");
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_Pol5Fixed.png");

  reject=true;
  TF1* fitCurveExpoWe = new TF1("fitCurveExpoWe",fexp,sbFitLoLo,sbFitHiHi,2);
  fitCurveExpoWe->SetParameter(0,-2.);
  //first just to check the status
  status = WeDataProjXexpo->Fit(fitCurveExpoWe,"L","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResultExpoWe = WeDataProjXexpo->Fit(fitCurveExpoWe,"SL","",sbFitLoLo,sbFitHiHi);
  TMatrixDSym covExpoWe = fitResultExpoWe->GetCovarianceMatrix();
  fitResultExpoWe->Print("V");
  WeDataProjXexpo->GetXaxis()->SetRangeUser(100,164.9);
  reject=false;
  YieldBinWidth=WeDataProjXexpo->GetBinWidth(1);
  Double_t WeExpoYield = fitCurveExpoWe->Integral(sbLoLo,sbHiHi)/YieldBinWidth;
  Double_t WeExpoYieldErr = fitCurveExpoWe->IntegralError(sbLoLo,sbHiHi,fitResultExpoWe->GetParams(),covExpoWe.GetMatrixArray() )/YieldBinWidth; 
  Double_t WeExpoYieldSig = fitCurveExpoWe->Integral(sigLo,sigHi)/YieldBinWidth;
  Double_t WeExpoYieldSigErr = fitCurveExpoWe->IntegralError(sigLo,sigHi,fitResultExpoWe->GetParams(),covExpoWe.GetMatrixArray() )/YieldBinWidth; 
  reject=true;
  sprintf(strWe,"Fit Yield: %4.2f #pm %4.2f",WeExpoYield,WeExpoYieldErr);
  sprintf(strWeSig,"Fit Signal Yield: %4.2f #pm %4.2f",WeExpoYieldSig,WeExpoYieldSigErr);
  sprintf(strWeRF,"RooFit signal Yield: %4.2f #pm %4.2f",integ_We,dinteg_We);
  TPaveText *textExpoWe= new TPaveText(.5,.65,.84,.84,"NDC");textExpoWe->SetFillStyle(0);textExpoWe->SetBorderSize(0);
  //textExpoWe->AddText(strWe);
  textExpoWe->AddText(strWeSig);
  WeDataProjXexpo->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  WeDataProjXexpo->SetTitle("");textExpoWe->SetTextFont(42);
  WeDataProjXexpo->Draw();
  textExpoWe->Draw("SAMES");
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_expo.png");

  reject=true;
  TF1* fitCurveLinWe = new TF1("fitCurveLinWe",flin,sbFitLoLo,sbFitHiHi,2);
  //first just to check the status
  status = WeDataProjXlin->Fit(fitCurveLinWe,"L","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResultLinWe = WeDataProjXlin->Fit(fitCurveLinWe,"SL","",sbFitLoLo,sbFitHiHi);
  TMatrixDSym covLinWe = fitResultLinWe->GetCovarianceMatrix();
  fitResultLinWe->Print("V");
  WeDataProjXlin->GetXaxis()->SetRangeUser(100,164.9);
  reject=false;
  Double_t WeLinYield = fitCurveLinWe->Integral(sbLoLo,sbHiHi);
  Double_t WeLinYieldErr = fitCurveLinWe->IntegralError(sbLoLo,sbHiHi,fitResultLinWe->GetParams(),covLinWe.GetMatrixArray() ); 
  Double_t WeLinYieldSig = fitCurveLinWe->Integral(sigLo,sigHi);
  Double_t WeLinYieldSigErr = fitCurveLinWe->IntegralError(sigLo,sigHi,fitResultLinWe->GetParams(),covLinWe.GetMatrixArray() ); 
  reject=true;
  sprintf(strWe,"Fit Yield: %4.2f #pm %4.2f",WeLinYield,WeLinYieldErr);
  sprintf(strWeSig,"Fit Signal Yield: %4.2f #pm %4.2f",WeLinYieldSig,WeLinYieldSigErr);
  sprintf(strWeRF,"RooFit signal Yield: %4.2f #pm %4.2f",integ_We,dinteg_We);
  TPaveText *textLinWe= new TPaveText(.5,.65,.84,.84,"NDC");textLinWe->SetFillStyle(0);textLinWe->SetBorderSize(0);
  //textLinWe->AddText(strWe);
  textLinWe->AddText(strWeSig);textLinWe->SetTextFont(42);
  WeDataProjXlin->Draw();
  textLinWe->Draw("SAMES");
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_linear.png");
  
  reject=true;
  TF1* fitCurvePol3We = new TF1("fitCurvePol3We",fpol3,sbFitLoLo,sbFitHiHi,4);
  //fitCurvePol3We->SetParameter(0,.03);
  //fitCurvePol3We->SetParameter(1,-.01);
  //fitCurvePol3We->SetParameter(2,-.01);
  //fitCurvePol3We->SetParameter(3,-.01);
  //first just to check the status
  status = WeDataProjXpol->Fit(fitCurvePol3We,"L","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResultPol3We = WeDataProjXpol->Fit(fitCurvePol3We,"SL","",sbFitLoLo,sbFitHiHi);
  TMatrixDSym covPol3We = fitResultPol3We->GetCovarianceMatrix();
  fitResultPol3We->Print("V");
  WeDataProjXpol->GetXaxis()->SetRangeUser(100,164.9);
  reject=false;
  Double_t WePol3Yield = fitCurvePol3We->Integral(sbLoLo,sbHiHi);
  Double_t WePol3YieldErr = fitCurvePol3We->IntegralError(sbLoLo,sbHiHi,fitResultPol3We->GetParams(),covPol3We.GetMatrixArray() ); 
  Double_t WePol3YieldSig = fitCurvePol3We->Integral(sigLo,sigHi);
  Double_t WePol3YieldSigErr = fitCurvePol3We->IntegralError(sigLo,sigHi,fitResultPol3We->GetParams(),covPol3We.GetMatrixArray() ); 
  reject=true;
  sprintf(strWe,"Fit Yield: %4.2f #pm %4.2f",WePol3Yield,WePol3YieldErr);
  sprintf(strWeSig,"Fit Signal Yield: %4.2f #pm %4.2f",WePol3YieldSig,WePol3YieldSigErr);
  sprintf(strWeRF,"RooFit signal Yield: %4.2f #pm %4.2f",integ_We,dinteg_We);
  TPaveText *textPol3We= new TPaveText(.5,.65,.84,.84,"NDC");textPol3We->SetFillStyle(0);textPol3We->SetBorderSize(0);
  //textPol3We->AddText(strWe);
  textPol3We->AddText(strWeSig);textPol3We->SetTextFont(42);
  WeDataProjXpol->Draw();
  textPol3We->Draw("SAMES");
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_pol3.png");

  reject=true;
  TF1* fitCurvePol5We = new TF1("fitCurvePol5We",fpol5,sbFitLoLo,sbFitHiHi,5);
  //fitCurvePol5We->SetParameter(0,300.);
  //fitCurvePol5We->SetParameter(1,-5.);
  //first just to check the status
  status = WeDataProjXpol->Fit(fitCurvePol5We,"L","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResultPol5We = WeDataProjXpol->Fit(fitCurvePol5We,"SL","",sbFitLoLo,sbFitHiHi);
  TMatrixDSym covPol5We = fitResultPol5We->GetCovarianceMatrix();
  fitResultPol5We->Print("V");
  WeDataProjXpol->GetXaxis()->SetRangeUser(100,164.9);
  reject=false;
  Double_t Pol5WeYield = fitCurvePol5We->Integral(sbLoLo,sbHiHi);
  Double_t Pol5WeYieldErr = fitCurvePol5We->IntegralError(sbLoLo,sbHiHi,fitResultPol5We->GetParams(),covPol5We.GetMatrixArray() ); 
  Double_t Pol5WeYieldSig = fitCurvePol5We->Integral(sigLo,sigHi);
  Double_t Pol5WeYieldSigErr = fitCurvePol5We->IntegralError(sigLo,sigHi,fitResultPol5We->GetParams(),covPol5We.GetMatrixArray() ); 
  reject=true;
  sprintf(str,"Fit Yield: %4.2f #pm %4.2f",Pol5WeYield,Pol5WeYieldErr);
  sprintf(strSig,"Fit Signal Yield: %4.2f #pm %4.2f",Pol5WeYieldSig,Pol5WeYieldSigErr);
  //sprintf(strRF,"RooFit signal Yield: %4.2f #pm %4.2f",integ_,dinteg_);
  TPaveText *textPol5We= new TPaveText(.5,.65,.84,.84,"NDC");textPol5We->SetFillStyle(0);textPol5We->SetBorderSize(0);
  //textPol5We->AddText(str);
  textPol5We->AddText(strSig);textPol5We->SetTextFont(42);
  WeDataProjXpol->Draw();
  textPol5We->Draw("SAMES");
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_pol5.png");

  

  //met plot with ratio
  RooAbsReal* igWeSig = MetWePdf_bg.createIntegral(xMetWe_bg,NormSet(xMetWe_bg),Range("sigWe_bg"));
  RooAbsReal* igWeSig2 = MetWePdf_bg.createIntegral(xMetWe_bg,Range("sigWe_bg"));
  RooAbsReal* igWeFULL = MetWePdf_bg.createIntegral(xMetWe_bg,NormSet(xMetWe_bg),Range("fullWe_bg"));
  RooAbsReal* igWeSBlo = MetWePdf_bg.createIntegral(xMetWe_bg,NormSet(xMetWe_bg),Range("sb_loWe_bg"));
  RooAbsReal* igWeSBhi = MetWePdf_bg.createIntegral(xMetWe_bg,NormSet(xMetWe_bg),Range("sb_hiWe_bg"));
  RooAbsReal* igWetest = MetWePdf_bg.createIntegral(xMetWe_bg);
  float HiggSigNoNormWe=igWeSig->getVal()*WeYield_bg.getVal();
  cout<<"poly 1 fit value: "<<PolWe1_bg.getVal()<<endl;//use these to make stupid signal fit to get yield.
  cout<<"signal integration old way : "<<HiggSigNoNormWe<<endl;
  cout<<"signal integration old way2: "<<igWeSig2->getVal()<<endl;
  cout<<"signal integration new way : "<<integ_We<<" +- "<<dinteg_We<<endl;
  cout<<" full integral: "<<igWeFULL->getVal()<<endl;
  cout<<" full integral test: "<<igWetest->getVal()<<endl;
  cout<<" full integral test2: "<<MetWePdf_bg.getVal()<<endl;
  float HiggLowSBNoNormWe=igWeSBlo->getVal()*WeYield_bg.getVal(),HiggHighSBNoNormWe=igWeSBhi->getVal()*WeYield_bg.getVal();
  p1->cd();
  sbLoBin1=WeData->GetXaxis()->FindBin(sbLoLo);sbLoBin2=WeData->GetXaxis()->FindBin(sbLoHi)-1;  
  sbHiBin1=WeData->GetXaxis()->FindBin(sbHiLo);sbHiBin2=WeData->GetXaxis()->FindBin(sbHiHi)-1;
  sigBin1=WeData->GetXaxis()->FindBin(sigLo);sigBin2=WeData->GetXaxis()->FindBin(sigHi)-1;
  TH1F* ggMetWe_SBloInvMass=(TH1F*)WeData->ProjectionY("ggMetWe_SBloInvMass",sbLoBin1,sbLoBin2,"eo");
  TH1F* ggMetWe_SBhiInvMass=(TH1F*)WeData->ProjectionY("ggMetWe_SBhiInvMass",sbHiBin1,sbHiBin2,"eo");
  TH1F* ggMetWe_sigInvMass=(TH1F*)WeData->ProjectionY("ggMetWe_sigInvMass",sigBin1,sigBin2,"eo");


  TH1F* ggMetWe_SBloInvMassRebin = (TH1F*)ggMetWe_SBloInvMass->Rebin(NmetBinsXtraWide,"ggMetWe_SBloInvMassRebin",xbinsXtraWide);
  TH1F* ggMetWe_SBhiInvMassRebin = (TH1F*)ggMetWe_SBhiInvMass->Rebin(NmetBinsXtraWide,"ggMetWe_SBhiInvMassRebin",xbinsXtraWide);
  TH1F* ggMetWe_sigInvMassRebin = (TH1F*)ggMetWe_sigInvMass->Rebin(NmetBinsXtraWide,"ggMetWe_sigInvMassRebin",xbinsXtraWide);

  AddOverflowToLastBin(ggMetWe_SBloInvMassRebin);
  AddOverflowToLastBin(ggMetWe_SBhiInvMassRebin);
  AddOverflowToLastBin(ggMetWe_sigInvMassRebin);

  for(int i=0;i<ggMetWe_sigInvMass->FindBin(31);i++){
    cout<<"We "<<ggMetWe_sigInvMass->GetBinLowEdge(i)<<"<"<<ggMetWe_sigInvMass->GetBinLowEdge(i+1)<<"  events: "<<ggMetWe_sigInvMass->GetBinContent(i)<<" +- "<<ggMetWe_sigInvMass->GetBinError(i)<<endl;
  }
  for(int i=0;i<ggMetWe_sigInvMassRebin->FindBin(31);i++){
    cout<<"We "<<ggMetWe_sigInvMassRebin->GetBinLowEdge(i)<<"<"<<ggMetWe_sigInvMassRebin->GetBinLowEdge(i+1)<<"  events: "<<ggMetWe_sigInvMassRebin->GetBinContent(i)<<" +- "<<ggMetWe_sigInvMassRebin->GetBinError(i)<<endl;
  }
  float higgslowscaleWe = ggMetWe_SBloInvMassRebin->Integral() ? WePowYieldSig/*(HiggSigNoNormWe)*//ggMetWe_SBloInvMassRebin->Integral() : 0;
  float higgshighscaleWe = ggMetWe_SBhiInvMassRebin->Integral() ? WePowYieldSig/*(HiggSigNoNormWe)*//ggMetWe_SBhiInvMassRebin->Integral() : 0;
  //cout<<"WePowYieldSig: "<<WePowYieldSig<<"  ggMetWe_SBloInvMassRebin->Integral(): "<<ggMetWe_SBloInvMassRebin->Integral()<<"  ggMetWe_SBhiInvMassRebin->Integral(): "<<ggMetWe_SBhiInvMassRebin->Integral()<<"  higgslowscaleWe: "<<higgslowscaleWe<<"  higgshighscaleWe: "<<higgshighscaleWe<<endl;
  //float WeFitSyst=(WeYield_bg.getError()*WeYield_bg.getError())/(WeYield_bg.getVal()*WeYield_bg.getVal());
  float WeFitSyst=(WePowYieldSigErr*WePowYieldSigErr)/(WePowYieldSig*WePowYieldSig);
  double WebgStatErrLO=0.;double WebgStatLO=ggMetWe_SBloInvMassRebin->IntegralAndError(0,-1,WebgStatErrLO);
  float WebgStatSystLO=(WebgStatErrLO*WebgStatErrLO)/(WebgStatLO*WebgStatLO);
  double WebgStatErrHI=0.;double WebgStatHI=ggMetWe_SBhiInvMassRebin->IntegralAndError(0,-1,WebgStatErrHI);
  float WebgStatSystHI=(WebgStatErrHI*WebgStatErrHI)/(WebgStatHI*WebgStatHI);
  /*  for(int i=0;i<ggMetWe_SBloInvMassRebin->GetNbinsX();i++){
      float Value=ggMetWe_SBloInvMassRebin->GetBinContent(i),err=ggMetWe_SBloInvMassRebin->GetBinError(i);
      float ValueNew=Value*higgslowscaleWe;
      float errNew=ValueNew*sqrt((err*err/Value/Value)+WeFitSyst+WebgStatSystLO);
      ggMetWe_SBloInvMassRebin->SetBinContent(i,ValueNew);ggMetWe_SBloInvMassRebin->SetBinError(i,errNew);
      Value=ggMetWe_SBhiInvMassRebin->GetBinContent(i),err=ggMetWe_SBhiInvMassRebin->GetBinError(i);
      ValueNew=Value*higgshighscaleWe;
      errNew=ValueNew*sqrt((err*err/Value/Value)+WeFitSyst+WebgStatSystHI);
      ggMetWe_SBhiInvMassRebin->SetBinContent(i,ValueNew);ggMetWe_SBhiInvMassRebin->SetBinError(i,errNew);
      }
      
      //ggMetWe_SBloInvMassRebin->Scale(higgslowscaleWe);
      //ggMetWe_SBhiInvMassRebin->Scale(higgshighscaleWe);
      for(int i=1;i<ggMetWe_sigInvMassRebin->GetNbinsX()+2;i++){
      float x = ggMetWe_SBloInvMassRebin->GetBinContent(i)/ggMetWe_SBloInvMassRebin->GetBinWidth(i);ggMetWe_SBloInvMassRebin->SetBinContent(i,x);
      x = ggMetWe_SBloInvMassRebin->GetBinError(i)/ggMetWe_SBloInvMassRebin->GetBinWidth(i);ggMetWe_SBloInvMassRebin->SetBinError(i,x);
      x = ggMetWe_SBhiInvMassRebin->GetBinContent(i)/ggMetWe_SBhiInvMassRebin->GetBinWidth(i);ggMetWe_SBhiInvMassRebin->SetBinContent(i,x);
      x = ggMetWe_SBhiInvMassRebin->GetBinError(i)/ggMetWe_SBhiInvMassRebin->GetBinWidth(i);ggMetWe_SBhiInvMassRebin->SetBinError(i,x);
      x = ggMetWe_sigInvMassRebin->GetBinContent(i)/ggMetWe_sigInvMassRebin->GetBinWidth(i);ggMetWe_sigInvMassRebin->SetBinContent(i,x);
      x = ggMetWe_sigInvMassRebin->GetBinError(i)/ggMetWe_sigInvMassRebin->GetBinWidth(i);ggMetWe_sigInvMassRebin->SetBinError(i,x);
      }
  */

  TH1F* ggMetWe_SBloInvMassClone = (TH1F*)ggMetWe_SBloInvMass->Clone("ggMetWe_SBloInvMassClone");
  TH1F* ggMetWe_SBhiInvMassClone = (TH1F*)ggMetWe_SBhiInvMass->Clone("ggMetWe_SBhiInvMassClone");
  ggMetWe_SBloInvMassClone->Scale(higgslowscale);ggMetWe_SBhiInvMassClone->Scale(higgshighscale);
  TH1F* ggMetWe_SBloInvMassCloneRebin = (TH1F*)ggMetWe_SBloInvMassClone->Rebin(NmetBinsXtraWide,"ggMetWe_SBloInvMassCloneRebin",xbinsXtraWide);
  TH1F* ggMetWe_SBhiInvMassCloneRebin = (TH1F*)ggMetWe_SBhiInvMassClone->Rebin(NmetBinsXtraWide,"ggMetWe_SBhiInvMassCloneRebin",xbinsXtraWide);
 
  TH1F* ggMetWe_InvMassSBcombNoRebin=(TH1F*)ggMetWe_SBloInvMass->Clone();//ggMetWe_InvMassSBcombNoRebin->Add(ggMetWe_SBhiInvMass);ggMetWe_InvMassSBcombNoRebin->Scale(0.5);
  TH1F *StatErrsWe       =(TH1F*)ggMetWe_SBloInvMass->Clone();
  TH1F *FitStatSystErrsWe =(TH1F*)ggMetWe_SBloInvMass->Clone();
  TH1F *FitShapeSystErrsWe  =(TH1F*)ggMetWe_SBloInvMass->Clone();
  TH1F *HalfDiffErrsWe  =(TH1F*)ggMetWe_SBloInvMass->Clone();
  double valTot=0.;
  for(int i=1;i<=ggMetWe_InvMassSBcombNoRebin->GetNbinsX();i++){
    Double_t val = ((ggMetWe_SBloInvMass->GetBinContent(i)*higgslowscaleWe)+(ggMetWe_SBhiInvMass->GetBinContent(i)*higgshighscaleWe))/2.;
    //Double_t val = (ggMetWe_SBloInvMass->GetBinContent(i)*higgslowscaleWe);
    //double Lerr=0.,Uerr=0.;
    //double Lval = ggMetWe_SBloInvMass->IntegralAndError(0,-1,Lerr);
    //double Uval = ggMetWe_SBhiInvMass->IntegralAndError(0,-1,Uerr);
    Double_t LvalScaleWe = ggMetWe_SBloInvMass->GetBinContent(i)*higgslowscaleWe,UvalScaleWe = ggMetWe_SBhiInvMass->GetBinContent(i)*higgshighscaleWe;
    Double_t WeStatErr(0.),WeFitStatSystErr(0.),WeFitShapeSystErr(0.),WeHalfDiffErr(0.);
    double err = TotErr(val,ggMetWe_SBloInvMass->GetBinContent(i),ggMetWe_SBloInvMass->GetBinError(i),ggMetWe_SBhiInvMass->GetBinContent(i),ggMetWe_SBhiInvMass->GetBinError(i),WebgStatLO,WebgStatErrLO,WebgStatHI,WebgStatErrHI,WePowYieldSig,WePowYieldSigErr,WeExpoYieldSig,WeStatErr,WeFitStatSystErr,WeFitShapeSystErr,LvalScaleWe,UvalScaleWe,WeHalfDiffErr);
    ggMetWe_InvMassSBcombNoRebin->SetBinContent(i,val);ggMetWe_InvMassSBcombNoRebin->SetBinError(i,err);
    StatErrsWe->SetBinError(i,WeStatErr);FitStatSystErrsWe->SetBinError(i,WeFitStatSystErr);FitShapeSystErrsWe->SetBinError(i,WeFitShapeSystErr);HalfDiffErrsWe->SetBinError(i,WeHalfDiffErr);
    valTot+=val;
  }
  cout<<"WePowYieldSig: "<<WePowYieldSig<<"  valTot: "<<valTot<<endl;

  TH1F* ggMetWe_InvMassSBcomb = (TH1F*)ggMetWe_InvMassSBcombNoRebin->Rebin(NmetBinsXtraWide,"ggMetWe_InvMassSBcomb",xbinsXtraWide);
  cout<<"ggMetWe_InvMassSBcomb after rebin: "<<ggMetWe_InvMassSBcomb->Integral()<<" (0,-1): "<<ggMetWe_InvMassSBcomb->Integral(0,-1)<<endl;
  for(int i=0;i<ggMetWe_sigInvMassRebin->FindBin(31);i++){
    cout<<"We "<<ggMetWe_sigInvMassRebin->GetBinLowEdge(i)<<"<"<<ggMetWe_sigInvMassRebin->GetBinLowEdge(i+1)<<"  events: "<<ggMetWe_sigInvMassRebin->GetBinContent(i)<<" +- "<<ggMetWe_sigInvMassRebin->GetBinError(i)<<endl;
  }
  TH1F* ggMetWe_sigInvMassRebin_Clone = (TH1F*)ggMetWe_sigInvMassRebin->Clone("ggMetWe_sigInvMassRebin_Clone");
  AddOverflowToLastBin(ggMetWe_InvMassSBcomb);
  TH1F* ggMetWe_InvMassSBcomb_Clone = (TH1F*)ggMetWe_InvMassSBcomb->Clone("ggMetWe_InvMassSBcomb_Clone");
  ggMetWe_InvMassSBcomb_Clone->Add(SMHiggsMetWeNoDivByWidth);
  cout<<"ggMetWe_InvMassSBcomb after overflow add: "<<ggMetWe_InvMassSBcomb->Integral()<<" (0,-1): "<<ggMetWe_InvMassSBcomb->Integral(0,-1)<<endl;
  double clone2err=0.;ggMetWe_InvMassSBcomb_Clone->IntegralAndError(0,-1,clone2err);
  cout<<"ggMetWe_InvMassSBcomb_plusSMhiggs after overflow add: "<<ggMetWe_InvMassSBcomb_Clone->Integral()<<" +- "<<clone2err<<endl;
  DivideBy15gev(ggMetWe_InvMassSBcomb);
  DivideBy15gev(ggMetWe_sigInvMassRebin);
  
  /*
  for(int i=0;i<ggMetWe_InvMassSBcomb->GetNbinsX()+2;i++){
    float x = ggMetWe_InvMassSBcomb->GetBinContent(i)/ggMetWe_InvMassSBcomb->GetBinWidth(i);ggMetWe_InvMassSBcomb->SetBinContent(i,x);
    x = ggMetWe_InvMassSBcomb->GetBinError(i)/ggMetWe_InvMassSBcomb->GetBinWidth(i);ggMetWe_InvMassSBcomb->SetBinError(i,x);
    x = ggMetWe_sigInvMassRebin->GetBinContent(i)/ggMetWe_sigInvMassRebin->GetBinWidth(i);ggMetWe_sigInvMassRebin->SetBinContent(i,x);
    x = ggMetWe_sigInvMassRebin->GetBinError(i)/ggMetWe_sigInvMassRebin->GetBinWidth(i);ggMetWe_sigInvMassRebin->SetBinError(i,x);
    }*/
  
  c1->cd();c1->SetLogy(0);
  ggMetWe_SBloInvMassCloneRebin->SetFillColor(kRed);ggMetWe_SBloInvMassCloneRebin->SetFillStyle(3004);ggMetWe_SBloInvMassCloneRebin->SetMarkerSize(0);
  ggMetWe_SBloInvMassCloneRebin->GetXaxis()->SetRangeUser(0,250);
  ggMetWe_SBhiInvMassCloneRebin->SetFillColor(kBlue);ggMetWe_SBhiInvMassCloneRebin->SetFillStyle(3004);ggMetWe_SBhiInvMassCloneRebin->SetMarkerSize(0);
  ggMetWe_SBloInvMassCloneRebin->SetLineColor(kRed);ggMetWe_SBloInvMassCloneRebin->SetMarkerColor(kRed);
  ggMetWe_SBhiInvMassCloneRebin->SetLineColor(kBlue);ggMetWe_SBhiInvMassCloneRebin->SetMarkerColor(kBlue);
  ggMetWe_SBhiInvMassCloneRebin->SetLineWidth(2);ggMetWe_SBloInvMassCloneRebin->SetLineWidth(2);
  ggMetWe_sigInvMassRebin_Clone->SetTitle("");ggMetWe_sigInvMassRebin_Clone->GetXaxis()->SetTitle("E_{T}^{miss}");ggMetWe_sigInvMassRebin_Clone->GetYaxis()->SetTitle("Events / 15 GeV");
  ggMetWe_sigInvMassRebin_Clone->GetXaxis()->SetLabelSize(0.06);
  ggMetWe_sigInvMassRebin_Clone->GetYaxis()->SetTitleSize(0.05);
  ggMetWe_sigInvMassRebin_Clone->Draw("PE");
  ggMetWe_SBhiInvMassCloneRebin->Draw("E2sames");ggMetWe_SBhiInvMassCloneRebin->Draw("pesames");
  ggMetWe_SBloInvMassCloneRebin->Draw("E2sames");ggMetWe_SBloInvMassCloneRebin->Draw("pesames");
  ggMetWe_sigInvMassRebin_Clone->Draw("PEsames");
  lowhighleg->Draw();
  c1->Print("Plots/Higgs/Exclusive_WeMet_SbLowAndHigh.png");
  c1->Print("Plots/Higgs/Exclusive_WeMet_SbLowAndHigh.pdf");
  c1->SetLogy(1);

  p1->cd();

  ggMetWe_sigInvMassRebin->SetMarkerColor(kBlack);ggMetWe_sigInvMassRebin->SetLineColor(kBlack);
  ggMetWe_SBloInvMassRebin->SetMarkerColor(kBlue);ggMetWe_SBloInvMassRebin->SetLineColor(kBlue);ggMetWe_SBloInvMassRebin->SetFillColor(kBlue);ggMetWe_SBloInvMassRebin->SetFillStyle(3004);
  ggMetWe_SBhiInvMassRebin->SetMarkerColor(kRed);ggMetWe_SBhiInvMassRebin->SetLineColor(kRed);ggMetWe_SBhiInvMassRebin->SetFillColor(kRed);ggMetWe_SBhiInvMassRebin->SetFillStyle(3004);
  ggMetWe_InvMassSBcomb->SetMarkerColor(42);ggMetWe_InvMassSBcomb->SetLineColor(42);ggMetWe_InvMassSBcomb->SetFillColor(42);//ggMetWe_InvMassSBcomb->SetFillStyle(3004);
  ggMetWe_SBloInvMassRebin->SetMarkerSize(0.75);ggMetWe_SBhiInvMassRebin->SetMarkerSize(0.75);ggMetWe_InvMassSBcomb->SetMarkerSize(0.75);
  WeMetStackRebin->Add(ggMetWe_InvMassSBcomb);
  TH1F* ggMetWe_InvMassSBcomb_plus_SMhiggs = (TH1F*)ggMetWe_InvMassSBcomb->Clone();ggMetWe_InvMassSBcomb_plus_SMhiggs->Add(SMHiggsMetWe);
  ggMetWe_InvMassSBcomb_plus_SMhiggs->SetMarkerSize(0);ggMetWe_InvMassSBcomb_plus_SMhiggs->SetFillColor(kBlack);ggMetWe_InvMassSBcomb_plus_SMhiggs->SetFillStyle(3004);


  //const double alpha = 1 - 0.6827;
  TGraphAsymmErrors *ggMetWe_sigInvMassRebin_tgraph = new TGraphAsymmErrors(ggMetWe_sigInvMassRebin);
  //ggMetWe_sigInvMassRebin_tgraph->SetMarkerSize(0.5);
  ggMetWe_sigInvMassRebin_tgraph->SetMarkerStyle(20);
  for (int i = 0; i < ggMetWe_sigInvMassRebin_tgraph->GetN(); ++i) {
    double N = ggMetWe_sigInvMassRebin_tgraph->GetY()[i];
    double L=0.,U=0.;
    double w = ggMetWe_sigInvMassRebin->GetBinWidth(i+1);
    L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N/(15./w),1.));
    U =  ROOT::Math::gamma_quantile_c(alpha/2,(N/(15./w))+1,1);
    L/=w/15.;U/=w/15.;
    ggMetWe_sigInvMassRebin_tgraph->SetPointEYlow(i, N-L);
    ggMetWe_sigInvMassRebin_tgraph->SetPointEYhigh(i, U-N);
    
    //      double U =  (N==0) ?  ( ROOT::Math::gamma_quantile_c(alpha,N+1,1) ) :
    //         ( ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) );
    //ggMetWeMT_sigInvMassRebin_tgraph->SetPointEYlow(i, N-L);
    //ggMetWeMT_sigInvMassRebin_tgraph->SetPointEYhigh(i, U-N);
    cout<<"ele i met:"<<i<<"  N: "<<N<<"  L: "<<L<<"  N-L: "<<N-L<<"  U: "<<U<<"  U-N: "<<U-N<<"  w: "<<w<<endl;
  }

  ggMetWe_sigInvMassRebin_tgraph->GetXaxis()->SetRangeUser(0,179);
  ggMetWe_sigInvMassRebin_tgraph->GetXaxis()->SetLabelSize(0);
  ggMetWe_sigInvMassRebin_tgraph->GetYaxis()->SetRangeUser(1e-2,70);
  ggMetWe_sigInvMassRebin_tgraph->GetYaxis()->SetTitle("Events / 15 GeV");
  ggMetWe_sigInvMassRebin_tgraph->GetYaxis()->SetTitleSize(0.05);
  ggMetWe_sigInvMassRebin_tgraph->SetTitle("");
  cout<<"We signal region full integral: "<<ggMetWe_sigInvMassRebin->Integral(0,-1)<<endl;
  for(int i=0;i<ggMetWe_sigInvMassRebin->FindBin(31);i++){
    cout<<"We "<<ggMetWe_sigInvMassRebin->GetBinLowEdge(i)<<"<"<<ggMetWe_sigInvMassRebin->GetBinLowEdge(i+1)<<"  events: "<<ggMetWe_sigInvMassRebin->GetBinContent(i)<<" +- "<<ggMetWe_sigInvMassRebin->GetBinError(i)<<endl;
  }
  //MakeBlindMet30(ggMetWe_sigInvMassRebin);
  for(int i=0;i<ggMetWe_sigInvMassRebin->FindBin(31);i++){
    cout<<"We "<<ggMetWe_sigInvMassRebin->GetBinLowEdge(i)<<"<"<<ggMetWe_sigInvMassRebin->GetBinLowEdge(i+1)<<"  events: "<<ggMetWe_sigInvMassRebin->GetBinContent(i)<<" +- "<<ggMetWe_sigInvMassRebin->GetBinError(i)<<endl;
  }

  TGraphAsymmErrors* ggMetWe_sigInvMassRebin_tgraph_clone = (TGraphAsymmErrors*)ggMetWe_sigInvMassRebin_tgraph->Clone("ggMetWe_sigInvMassRebin_tgraph_clone");

  ggMetWe_sigInvMassRebin_tgraph_clone->SetMarkerSize(2);
  ggMetWe_sigInvMassRebin_tgraph_clone->SetLineWidth(3);
  for(int i=0;i<ggMetWe_sigInvMassRebin_tgraph_clone->GetN();++i){
    ggMetWe_sigInvMassRebin_tgraph_clone->SetPointEXlow(i,0);
    ggMetWe_sigInvMassRebin_tgraph_clone->SetPointEXhigh(i,0);
  }
  ggMetWe_sigInvMassRebin_tgraph_clone->Draw("APE");
  WeMetStackRebin->Draw("histoSAMES");
  //MakeBlindMet30(ggMetWe_InvMassSBcomb_plus_SMhiggs);
  ggMetWe_InvMassSBcomb_plus_SMhiggs->Draw("E2SAMES");
  double WeSMest=0.,WeSMestErr=0.;WeSMest=ggMetWe_InvMassSBcomb_plus_SMhiggs->IntegralAndError(0,-1,WeSMestErr);
  cout<<"We background estimate full integral: "<<WeSMest<<"  and error: "<<WeSMestErr<<endl;
  cout<<"We SM higgs full integral: "<<SMHiggsMetWe->Integral(0,-1)<<endl;
  ggMetWe_sigInvMassRebin_tgraph_clone->Draw("PESAMES");
  //cout<<"We signal mX=130 full integral: "<<WeAAW130ProjYRebin->Integral(0,-1)<<endl;
  //cout<<"We signal mX=275 full integral: "<<WeAAW275ProjYRebin->Integral(0,-1)<<endl;
  //WeAAW130ProjYRebin->Draw("histoSAMES");
  //WeAAW275ProjYRebin->Draw("histoSAMES");
  h_SMS_WH_gg_1Ele_Excluded_Rebin->SetLineStyle(7);
  h_SMS_ZH_gg_1Ele_Excluded_Rebin->SetLineStyle(9);
  TH1F* h_SMS_HH_all_1Ele_Excluded_Rebin = (TH1F*)h_SMS_HH_2W2g_gg_1Ele_Excluded_Rebin->Clone();h_SMS_HH_all_1Ele_Excluded_Rebin->Add(h_SMS_HH_2Z2g_gg_1Ele_Excluded_Rebin);h_SMS_HH_all_1Ele_Excluded_Rebin->Add(h_SMS_HH_2tau2g_gg_1Ele_Excluded_Rebin);
  h_SMS_HH_all_1Ele_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_WH_gg_1Ele_notExcluded_Rebin->Draw("histoSAMES");
  h_SMS_ZH_gg_1Ele_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_WH_gg_1Ele_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_HH_2Z2g_gg_1Ele_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_HH_2W2g_gg_1Ele_Excluded_Rebin->Draw("histoSAMES");
  ggMetWe_sigInvMassRebin_tgraph_clone->Draw("PESAMES");
  p1->RedrawAxis();
  PrelimText->Draw();
  TPaveText *ggTextWe = new TPaveText(.37,.77,.57,.88,"NDC");
  //ggTextWe->AddText("#gamma#gamma + e, 0 #mu,  #leq1 b-Jet");
  //ggTextWe->AddText("#gamma#gamma + e, 0 #mu");
  ggTextWe->AddText("#gamma#gamma + e");
  ////ggTextWe->AddText("");
  ggTextWe->SetFillStyle(0);ggTextWe->SetTextFont(42);
  ggTextWe->SetFillColor(0);
  ggTextWe->SetBorderSize(0);
  //ggTextWe->Draw();
  //TLegend *WeMetLeg = new TLegend(.52,.57,.92,.78,"","brNDC");
  ggMetWe_sigInvMassRebin->SetMarkerSize(2);
  ggMetWe_sigInvMassRebin->SetLineWidth(2);
  TLegend *WeMetLeg = new TLegend(.57,.3,.97,.88,"","brNDC");
  WeMetLeg->SetFillColor(kWhite);WeMetLeg->SetTextSize(.05);WeMetLeg->SetTextFont(42);
  WeMetLeg->AddEntry(ggMetWe_sigInvMassRebin,"Data","ep");
  WeMetLeg->AddEntry(ggMetWe_InvMassSBcomb,"Non-higgs SM bg","f");
  WeMetLeg->AddEntry(SMHiggsAdditionForWe,"SM higgs","f");
  WeMetLeg->AddEntry(h_SMS_WH_gg_1Ele_Excluded_Rebin,   "Signal","");
  WeMetLeg->AddEntry(h_SMS_WH_gg_1Ele_Excluded_Rebin,   "hW, m_{#tilde{#chi}_{1}^{#pm}}=130 GeV","l");
  WeMetLeg->AddEntry(h_SMS_HH_all_1Ele_Excluded_Rebin,"hh,  m_{#tilde{#chi}_{1}^{0}}=130 GeV","l");
  WeMetLeg->AddEntry(h_SMS_ZH_gg_1Ele_Excluded_Rebin,"hZ,  m_{#tilde{#chi}_{1}^{0}}=130 GeV","l");

  /*
  WeMetLeg->AddEntry(WeWZHggProjYRebin,"SM W/ZH, m_{h}=126 GeV","f");
  WeMetLeg->AddEntry(WeVBFHggProjYRebin,"SM VBFH, m_{h}=126 GeV","f");
  WeMetLeg->AddEntry(WeTTHggProjYRebin,"SM TTH, m_{h}=126 GeV","f");
  WeMetLeg->AddEntry(WeggHggProjYRebin,"SM gg->H, m_{h}=126 GeV","f");
  */
  //WeMetLeg->AddEntry(WeWZHggProjYRebin,"SM W/ZH","f");
  //WeMetLeg->AddEntry(WeVBFHggProjYRebin,"SM VBFH","f");
  // WeMetLeg->AddEntry(WeTTHggProjYRebin,"SM TTH","f");
  //WeMetLeg->AddEntry(WeggHggProjYRebin,"SM gg->H","f");
  ///////////WeMetLeg->AddEntry(ggMetWe_InvMassSBcomb_plus_SMhiggs,"Full SM background error","f");
  //WeMetLeg->AddEntry(WeAAW130ProjYRebin,"AAW m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  //WeMetLeg->AddEntry(WeAAW275ProjYRebin,"AAW m_{#tilde{H}}=275, m_{#tilde{B}}=1","l");
  //WeMetLeg->AddEntry(h_SMS_WH_gg_1Ele_Excluded_Rebin,   "SMS WH+ZH m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  //WeMetLeg->AddEntry(h_SMS_WH_gg_1Ele_notExcluded_Rebin,"SMS WH+ZH m_{#tilde{H}}=400, m_{#tilde{B}}=150","l");
  /*
  WeMetLeg->AddEntry(h_SMS_WH_gg_1Ele_Excluded_Rebin,   "SMS WH m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  WeMetLeg->AddEntry(h_SMS_HH_2W2g_gg_1Ele_Excluded_Rebin,"SMS H(WW)H m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  WeMetLeg->AddEntry(h_SMS_ZH_gg_1Ele_Excluded_Rebin,"SMS ZH m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  WeMetLeg->AddEntry(h_SMS_HH_2Z2g_gg_1Ele_Excluded_Rebin,"SMS H(ZZ)H m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  */
  /*
  WeMetLeg->AddEntry(h_SMS_WH_gg_1Ele_Excluded_Rebin,   "Signal hW","l");
  WeMetLeg->AddEntry(h_SMS_HH_2W2g_gg_1Ele_Excluded_Rebin,"SMS H(WW)H","l");
  WeMetLeg->AddEntry(h_SMS_ZH_gg_1Ele_Excluded_Rebin,"Signal hZ","l");
  WeMetLeg->AddEntry(h_SMS_HH_2Z2g_gg_1Ele_Excluded_Rebin,"SMS H(ZZ)H","l");
  */
  WeMetLeg->SetFillStyle(0);WeMetLeg->SetBorderSize(0);
  //WeMetLeg->Draw();
  TLegend *WeMetLeg2 = new TLegend(.53,.4,.89,.57,"","brNDC");
  WeMetLeg2->SetFillColor(kWhite);
  WeMetLeg2->SetNColumns(2);
  /*WeMetLeg2->AddEntry(WeWZHggProjYRebin,"SM W/ZH","f");
  WeMetLeg2->AddEntry(WeVBFHggProjYRebin,"SM VBFH","f");
  WeMetLeg2->AddEntry(WeTTHggProjYRebin,"SM TTH","f");
  WeMetLeg2->AddEntry(WeggHggProjYRebin,"SM gg->H","f");*/
  //WeMetLeg2->AddEntry(h_SMS_WH_gg_1Ele_Excluded_Rebin,   "Signal hW, m_{#tilde{#chi}_{2}^{0}}=m_{#tilde{#chi}_{1}^{#pm}}=130 GeV","l");
  //WeMetLeg2->AddEntry(h_SMS_HH_2W2g_gg_1Ele_Excluded_Rebin,"SMS H(WW)H","l");
  //WeMetLeg2->AddEntry(h_SMS_HH_all_1Ele_Excluded_Rebin,"Signal hh","l");
  //WeMetLeg2->AddEntry(h_SMS_ZH_gg_1Ele_Excluded_Rebin,"Signal hZ","l");

  TPaveText *textSMSlog = new TPaveText(.775,.485,.94,.575,"NDC");
  //textSMSlog->AddText("m_{#chi_{2}}=130, m_{#chi_{1}}=1 GeV");
  textSMSlog->AddText("m_{#tilde{#chi}_{1}}=130 GeV");
  textSMSlog->SetFillStyle(0);textSMSlog->SetBorderSize(0);textSMSlog->SetTextFont(42);
  textSMSlog->Draw();
  //TPaveText *textSMS = new TPaveText(.72,.395,.87,.475,"NDC");
  TPaveText *textSMS = new TPaveText(.742,.29,.892,.37,"NDC");
  //textSMS->AddText("m_{#chi_{1}}=130, m_{#chi_{0}}=1 GeV");
  //textSMS->AddText("m_{#tilde{#chi}_{1}}=130 GeV");
  //textSMS->AddText("m_{#tilde{#chi}_{1}^{0}}=1 GeV");
  textSMS->SetFillStyle(0);textSMS->SetBorderSize(0);textSMS->SetTextSize(0.035);textSMS->SetTextFont(42);
  //textSMS->Draw();


  //WeMetLeg2->AddEntry(h_SMS_HH_2Z2g_gg_1Ele_Excluded_Rebin,"SMS H(ZZ)H","l");
  WeMetLeg2->SetFillStyle(0);WeMetLeg2->SetBorderSize(0);
  //WeMetLeg2->Draw();
  TPaveText *ggTextWeLog = new TPaveText(.33,.79,.63,.9,"NDC");
  ggTextWeLog->AddText("#gamma#gamma + e");
  ////ggTextWe->AddText("");
  ggTextWeLog->SetFillStyle(0);ggTextWeLog->SetTextFont(42);
  ggTextWeLog->SetFillColor(0);
  ggTextWeLog->SetBorderSize(0);
  ggTextWeLog->Draw();
  TLegend *WeMetLegLog = new TLegend(.59,.66,.93,.91,"","brNDC");
  WeMetLegLog->SetFillColor(kWhite);
  WeMetLegLog->AddEntry(ggMetWe_sigInvMassRebin,"Data","elpz");
  WeMetLegLog->AddEntry(ggMetWe_InvMassSBcomb,"Non-higgs SM bg","f");
  WeMetLegLog->AddEntry(SMHiggsAdditionForWe,"SM higgs","f");
  WeMetLegLog->SetFillStyle(0);WeMetLegLog->SetBorderSize(0);
  WeMetLegLog->Draw();
  TLegend *WeMetLeg2Log = new TLegend(.595,.495,.94,.675,"","brNDC");
  WeMetLeg2Log->SetFillColor(kWhite);
  WeMetLeg2Log->SetNColumns(2);
  WeMetLeg2Log->AddEntry(h_SMS_WH_gg_1Ele_Excluded_Rebin,"Signal hW","l");
  WeMetLeg2Log->AddEntry(h_SMS_HH_all_1Ele_Excluded_Rebin,"Signal hh","l");
  WeMetLeg2Log->AddEntry(h_SMS_ZH_gg_1Ele_Excluded_Rebin,"Signal hZ","l");
  WeMetLeg2Log->SetFillStyle(0);WeMetLeg2Log->SetBorderSize(0);
  WeMetLeg2Log->Draw();
  p2->cd();
  TH1F* ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs = (TH1F*)ggMetWe_sigInvMassRebin->Clone();
  //ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs->Divide(ggMetWe_InvMassSBcomb_plus_SMhiggs);
  TH1F* h_SystErrWe = (TH1F*)ggMetWe_InvMassSBcomb_plus_SMhiggs->Clone();h_SystErrWe->SetFillColor(kBlack);h_SystErrWe->SetFillStyle(3004);h_SystErrWe->SetMarkerSize(0);
  for(int i=1;i<ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs->GetNbinsX()+1;i++){
    float Value = ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs->GetBinContent(i);float StatErr = ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs->GetBinError(i);
    Value/=ggMetWe_InvMassSBcomb_plus_SMhiggs->GetBinContent(i);StatErr/=ggMetWe_InvMassSBcomb_plus_SMhiggs->GetBinContent(i);
    ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs->SetBinContent(i,Value);ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs->SetBinError(i,StatErr);
    float SystErr=h_SystErrWe->GetBinError(i)/h_SystErrWe->GetBinContent(i);
    h_SystErrWe->SetBinContent(i,1); h_SystErrWe->SetBinError(i,SystErr);
  }


  TGraphAsymmErrors* ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph = new TGraphAsymmErrors(ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs);
  for (int i = 0; i < ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetN(); ++i) {

    float N = ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetY()[i];
    double U = ggMetWe_sigInvMassRebin_tgraph->GetErrorYhigh(i);
    double L = ggMetWe_sigInvMassRebin_tgraph->GetErrorYlow(i);
    double Utemp=U,Ltemp=L;
    U/=ggMetWe_InvMassSBcomb_plus_SMhiggs->GetBinContent(i+1);L/=ggMetWe_InvMassSBcomb_plus_SMhiggs->GetBinContent(i+1);
    ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetPointEYlow(i, L);
    ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetPointEYhigh(i, U);
    if(i==0)cout<<"Ratio:"<<endl;
    cout<<"i:"<<i<<"  N: "<<N<<"  Ltemp: "<<Ltemp<<"Lbin: "<<ggMetWe_InvMassSBcomb_plus_SMhiggs->GetBinContent(i+1)<<"  L: "<<L<<"  Utemp: "<<Utemp<<"  Ubin: "<<ggMetWe_InvMassSBcomb_plus_SMhiggs->GetBinContent(i+1)<<" U: "<<U<<endl;
  }


  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetLineColor(kBlack);//ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs->SetMarkerColor(kBlack);//ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs->SetMarkerSize(0.75);
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetTitle("");
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetTitle("");
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetRangeUser(0.,5);
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetTitle("");
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetTitle("#frac{Data}{Prediction}");
  //ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetTitle("");
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetTitleOffset(0.5);
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetTitleSize(0.12);
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetTitleSize(0.2);
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetLabelSize(0.12);
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetNdivisions(205,0);
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetLabelSize(0.13);
  //MakeBlindMet30(ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs);
  //MakeBlindMet30(h_SystErrWe);
  //ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs->GetXaxis()->SetLabelSize(0.);
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetRangeUser(0,179);
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetMarkerSize(2);
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetLineWidth(3);
  for(int i=0;i<ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetN();++i){
    ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetPointEXlow(i,0);
    ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetPointEXhigh(i,0);
  }
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Draw("APE");
  h_SystErrWe->Draw("E2SAMES");
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Draw("PESAMES");
  l1.DrawLine(0,1,metPlotXmaxXtraWide,1);
  //PrelimText->Draw();
  p3->cd();p3->Clear();
  met2->Draw();
  c2->Print("Plots/Higgs/Exclusive_WeDataMetWithRatio.png");
  c2->Print("Plots/Higgs/Exclusive_WeDataMetWithRatio.pdf");
  p1->cd();p1->SetLogy(0);
  ggMetWe_sigInvMassRebin_tgraph_clone->GetYaxis()->SetRangeUser(0,10);
  ggMetWe_sigInvMassRebin_tgraph_clone->SetMarkerSize(2);
  ggMetWe_sigInvMassRebin_tgraph_clone->SetLineWidth(3);
  ggMetWe_sigInvMassRebin_tgraph_clone->Draw("APE");
  WeMetStackRebin->Draw("histoSAMES");
  ggMetWe_InvMassSBcomb_plus_SMhiggs->Draw("E2SAMES");
  ggMetWe_sigInvMassRebin_tgraph_clone->Draw("PESAMES");
  h_SMS_HH_all_1Ele_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_ZH_gg_1Ele_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_WH_gg_1Ele_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_HH_2Z2g_gg_1Ele_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_HH_2W2g_gg_1Ele_Excluded_Rebin->Draw("histoSAMES");
  ggMetWe_sigInvMassRebin_tgraph_clone->Draw("PESAMES");
  p1->RedrawAxis();
  PrelimText->Draw();
  ggTextWe->Draw();
  WeMetLeg->Draw();
  WeMetLeg2->Draw();
  textSMS->Draw();
  TLine metVertErr(104.13,8.6,104.13,9.4);metVertErr.SetLineWidth(3);
  metVertErr.Draw();
  p2->cd();
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetMarkerSize(2);
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetLineWidth(3);
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Draw("APE");
  h_SystErrWe->Draw("E2SAMES");
  ggMetWe_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Draw("PESAMES");
  l1.DrawLine(0,1,metPlotXmaxXtraWide,1);
  //PrelimText->Draw();
  p3->cd();p3->Clear();
  met2->Draw();
  c2->Print("Plots/Higgs/Exclusive_WeDataMetWithRatio_linear.png");
  c2->Print("Plots/Higgs/Exclusive_WeDataMetWithRatio_linear.pdf");
  p1->SetLogy(1);
 
  c1->cd();
  int metCut = WeData->GetYaxis()->FindBin(50);
  TH1D* WeDataProjX_metCut = (TH1D*)WeData->ProjectionX("WeDataProjX_metCut",metCut,-1,"eo");
  TH1D* WeggHggProjX_metCut = (TH1D*)WeggHgg->ProjectionX("WeggHggProjX_metCut",metCut,-1,"eo");
  TH1D* WeWZHggProjX_metCut = (TH1D*)WeWZHgg->ProjectionX("WeWZHggProjX_metCut",metCut,-1,"eo");
  TH1D* WeTTHggProjX_metCut = (TH1D*)WeTTHgg->ProjectionX("WeTTHggProjX_metCut",metCut,-1,"eo");
  TH1D* WeVBFHggProjX_metCut = (TH1D*)WeVBFHgg->ProjectionX("WeVBFHggProjX_metCut",metCut,-1,"eo");
  //TH1D* WeAAW130ProjX_metCut = (TH1D*)We_aaW130->ProjectionX("WeAAW130ProjX_metCut",metCut,-1,"eo");
  //TH1D* WeAAW275ProjX_metCut = (TH1D*)We_aaW275->ProjectionX("WeAAW275ProjX_metCut",metCut,-1,"eo");
  WeDataProjX_metCut->Rebin(2);WeggHggProjX_metCut->Rebin(2);WeWZHggProjX_metCut->Rebin(2);WeTTHggProjX_metCut->Rebin(2);WeVBFHggProjX_metCut->Rebin(2);
  //WeAAW130ProjX_metCut->Rebin(2);WeAAW275ProjX_metCut->Rebin(2);
  //MakeBlindInvMass(WeDataProjX_metCut);
  //TH1F* ggMetWe_SB_metCut = (TH1F*)ggMetWe_SBloInvMass_metCut->Clone();
  //ggMetWe_SB_metCut->Add(ggMetWe_SBloInvMass_metCut,ggMetWe_SBhiInvMass_metCut,0.5,0.5);
  THStack *WeBGSMHiggs = new THStack("WeBGSMHiggs","");
  WeBGSMHiggs->Add(WeTTHggProjX_metCut);WeBGSMHiggs->Add(WeVBFHggProjX_metCut);WeBGSMHiggs->Add(WeWZHggProjX_metCut);WeBGSMHiggs->Add(WeggHggProjX_metCut);
  //WeBGSMHiggs->Add(ggMetWe_SB_metCut);
  WeDataProjX_metCut->GetXaxis()->SetRangeUser(0,300);
  WeDataProjX_metCut->Draw("PE");
  WeBGSMHiggs->Draw("histoSAMES");
  //WeAAW130ProjX_metCut->Draw("histoSAMES");
  //WeAAW275ProjX_metCut->Draw("histoSAMES");
  c1->Print("Plots/Higgs/Exclusive_WeInvMass.png");
  //W->e MT plot
  c1->cd();
  TH2F* WeMTData = (TH2F*)fin->Get("ggMTvsInvarMass_Loose_1Ele_0_1Jets");
  TH2F* WeMTggHgg = (TH2F*)f_ggHgg->Get("ggMTvsInvarMass_Loose_1Ele_0_1Jets");
  TH2F* WeMTWZHgg = (TH2F*)f_WZHgg->Get("ggMTvsInvarMass_Loose_1Ele_0_1Jets");
  TH2F* WeMTTTHgg = (TH2F*)f_TTHgg->Get("ggMTvsInvarMass_Loose_1Ele_0_1Jets");
  TH2F* WeMTVBFHgg = (TH2F*)f_VBFHgg->Get("ggMTvsInvarMass_Loose_1Ele_0_1Jets");

  float ggHnormSF_ele = met_Phi_ele_ggH->GetEntries()/met_Phi_ele_ggH->Integral(0,-1);
  float WZHnormSF_ele = met_Phi_ele_WZH->GetEntries()/met_Phi_ele_WZH->Integral(0,-1);
  float TTHnormSF_ele = met_Phi_ele_TTH->GetEntries()/met_Phi_ele_TTH->Integral(0,-1);
  float VBFHnormSF_ele = met_Phi_ele_VBFH->GetEntries()/met_Phi_ele_VBFH->Integral(0,-1);

  WeMTggHgg->Scale((ggHnormSF_ele*PhoEffScale2*L_int*2.29e-03*19.22)/99989.);//125GeV=/96290);  
  WeMTVBFHgg->Scale((WZHnormSF_ele*PhoEffScale2*L_int*2.29e-03*1.544)/95677.);//125GeV=/99885); 
  WeMTTTHgg->Scale((TTHnormSF_ele*PhoEffScale2*L_int*2.29e-03*.1271)/100048.);//125GeV=/100224);
  WeMTWZHgg->Scale((VBFHnormSF_ele*PhoEffScale2*L_int*2.29e-03*(.6782/**(.3257+.014)*/+.3843/**.2*/))/100320);//125 and 126 GeV have same # events
  WeMTggHgg->SetLineColor(kGreen);WeMTggHgg->SetMarkerColor(kGreen);WeMTggHgg->SetFillColor(kGreen);
  WeMTWZHgg->SetLineColor(kCyan);WeMTWZHgg->SetMarkerColor(kCyan);WeMTWZHgg->SetFillColor(kCyan);
  WeMTTTHgg->SetLineColor(31);WeMTTTHgg->SetMarkerColor(31);WeMTTTHgg->SetFillColor(31);
  WeMTVBFHgg->SetLineColor(kRed+3);WeMTVBFHgg->SetMarkerColor(kRed+3);WeMTVBFHgg->SetFillColor(kRed+3);
  //WeMTggHgg->SetFillStyle(0);WeMTWZHgg->SetFillStyle(0);WeMTTTHgg->SetFillStyle(0);WeMTVBFHgg->SetFillStyle(0);
  WeMTggHgg->SetLineWidth(2);WeMTWZHgg->SetLineWidth(2);WeMTVBFHgg->SetLineWidth(2);WeMTTTHgg->SetLineWidth(2);
  TH1D* WeMTDataProjY = (TH1D*)WeMTData->ProjectionY("WeMTDataProjY",0,-1,"eo");
  TH1D* WeMTggHggProjY = (TH1D*)WeMTggHgg->ProjectionY("WeMTggHggProjY",0,-1,"eo");
  TH1D* WeMTWZHggProjY = (TH1D*)WeMTWZHgg->ProjectionY("WeMTWZHggProjY",0,-1,"eo");
  TH1D* WeMTTTHggProjY = (TH1D*)WeMTTTHgg->ProjectionY("WeMTTTHggProjY",0,-1,"eo");
  TH1D* WeMTVBFHggProjY = (TH1D*)WeMTVBFHgg->ProjectionY("WeMTVBFHggProjY",0,-1,"eo");
  //MakeBlindInvMass(WeMTDataProjX);
  c1->SetLogy(1);
  WeMTDataProjY->GetXaxis()->SetRangeUser(0,200);WeMTDataProjY->Draw();
  c1->Print("Plots/Higgs/Exclusive_WeDataMT.png");
  c1->Print("Plots/Higgs/Exclusive_WeDataMT.pdf");
  /*
  overFlowBin=WeMTggHggProjY->FindBin(metPlotXmax+1);lastBin=WeMTggHggProjY->FindBin(metPlotXmax-1);
  overFlow=0.;overFlowErr=0.;
  val=WeMTggHggProjY->GetBinContent(lastBin);valErr=WeMTggHggProjY->GetBinError(lastBin);
  overFlow=WeMTggHggProjY->IntegralAndError(overFlowBin,-1,overFlowErr);
  WeMTggHggProjY->SetBinContent(lastBin,val+overFlow);
  WeMTggHggProjY->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  val=WeMTTTHggProjY->GetBinContent(lastBin);valErr=WeMTTTHggProjY->GetBinError(lastBin);
  overFlow=WeMTTTHggProjY->IntegralAndError(overFlowBin,-1,overFlowErr);
  WeMTTTHggProjY->SetBinContent(lastBin,val+overFlow);
  WeMTTTHggProjY->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  val=WeMTWZHggProjY->GetBinContent(lastBin);valErr=WeMTWZHggProjY->GetBinError(lastBin);
  overFlow=WeMTWZHggProjY->IntegralAndError(overFlowBin,-1,overFlowErr);
  WeMTWZHggProjY->SetBinContent(lastBin,val+overFlow);
  WeMTWZHggProjY->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  val=WeMTVBFHggProjY->GetBinContent(lastBin);valErr=WeMTVBFHggProjY->GetBinError(lastBin);
  overFlow=WeMTVBFHggProjY->IntegralAndError(overFlowBin,-1,overFlowErr);
  WeMTVBFHggProjY->SetBinContent(lastBin,val+overFlow);
  WeMTVBFHggProjY->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  */
  THStack *WeMTStack = new THStack("WeMTStack","");
  WeMTStack->Add(WeMTTTHggProjY);WeMTStack->Add(WeMTVBFHggProjY);WeMTStack->Add(WeMTWZHggProjY);WeMTStack->Add(WeMTggHggProjY);
  WeMTStack->Draw("histo");WeMTStack->GetXaxis()->SetRangeUser(0,249.9);
  float SMHiggsMetIntWeMT=WeMTTTHggProjY->Integral()+WeMTVBFHggProjY->Integral()+WeMTWZHggProjY->Integral()+WeMTggHggProjY->Integral();
  c1->Print("Plots/Higgs/Exclusive_WeSMHiggsMT.png");
  TH1F* WeMTggHggProjYRebin=(TH1F*)WeMTggHggProjY->Rebin(nMTbins,"WeMTggHggProjYRebin",MTbins);
  TH1F* WeMTTTHggProjYRebin=(TH1F*)WeMTTTHggProjY->Rebin(nMTbins,"WeMTTTHggProjYRebin",MTbins);
  TH1F* WeMTWZHggProjYRebin=(TH1F*)WeMTWZHggProjY->Rebin(nMTbins,"WeMTWZHggProjYRebin",MTbins);
  TH1F* WeMTVBFHggProjYRebin=(TH1F*)WeMTVBFHggProjY->Rebin(nMTbins,"WeMTVBFHggProjYRebin",MTbins);
 
  AddOverflowToLastBin(WeMTggHggProjYRebin);
  AddOverflowToLastBin(WeMTTTHggProjYRebin);
  AddOverflowToLastBin(WeMTWZHggProjYRebin);
  AddOverflowToLastBin(WeMTVBFHggProjYRebin);
  DivideBy30gev(WeMTggHggProjYRebin);
  DivideBy30gev(WeMTTTHggProjYRebin);
  DivideBy30gev(WeMTWZHggProjYRebin);
  DivideBy30gev(WeMTVBFHggProjYRebin);

  THStack *WeMTStackRebin = new THStack("WeMTStackRebin","");
  THStack *WeMTStackRebin2 = new THStack("WeMTStackRebin2","");
  WeMTTTHggProjYRebin->SetLineColor(kBlack);WeMTVBFHggProjYRebin->SetLineColor(kBlack);WeMTggHggProjYRebin->SetLineColor(kBlack);WeMTWZHggProjYRebin->SetLineColor(kBlack);
  WeMTTTHggProjYRebin->SetLineWidth(1);WeMTVBFHggProjYRebin->SetLineWidth(1);WeMTggHggProjYRebin->SetLineWidth(1);WeMTWZHggProjYRebin->SetLineWidth(1);
  TH1F* SMHiggsAdditionForWeMT = (TH1F*)WeMTTTHggProjYRebin->Clone();SMHiggsAdditionForWeMT->Add(WeMTVBFHggProjYRebin);SMHiggsAdditionForWeMT->Add(WeMTWZHggProjYRebin);SMHiggsAdditionForWeMT->Add(WeMTggHggProjYRebin);
  WeMTStackRebin2->Add(WeMTTTHggProjYRebin);WeMTStackRebin2->Add(WeMTVBFHggProjYRebin);WeMTStackRebin2->Add(WeMTWZHggProjYRebin);WeMTStackRebin2->Add(WeMTggHggProjYRebin);
  WeMTStackRebin->Add(SMHiggsAdditionForWeMT);
  WeMTStackRebin2->Draw("histo");WeMTStackRebin2->GetXaxis()->SetRangeUser(0,249.9);WeMTStackRebin2->SetMinimum(3.1E-3);WeMTStackRebin2->SetMaximum(6);
  WeMTStackRebin2->GetXaxis()->SetTitle("M_{T} [GeV]");WeMTStackRebin2->GetYaxis()->SetTitle("Events / 30 GeV");
  legHiggs2->Draw();
  c1->Print("Plots/Higgs/Exclusive_WeSMHiggsMTRebin.png");
  c1->Print("Plots/Higgs/Exclusive_WeSMHiggsMTRebin.pdf");
  TH1F* SMHiggsMTWeMT = (TH1F*)WeMTTTHggProjYRebin->Clone();SMHiggsMTWeMT->Add(WeMTVBFHggProjYRebin);SMHiggsMTWeMT->Add(WeMTggHggProjYRebin);SMHiggsMTWeMT->Add(WeMTWZHggProjYRebin);
  //MT plot with ratio
  p1->cd();
  sbLoBin1=WeMTData->GetXaxis()->FindBin(sbLoLo);sbLoBin2=WeMTData->GetXaxis()->FindBin(sbLoHi)-1;  
  sbHiBin1=WeMTData->GetXaxis()->FindBin(sbHiLo);sbHiBin2=WeMTData->GetXaxis()->FindBin(sbHiHi)-1;
  sigBin1=WeMTData->GetXaxis()->FindBin(sigLo);sigBin2=WeMTData->GetXaxis()->FindBin(sigHi)-1;
  TH1F* ggMetWeMT_SBloInvMass=(TH1F*)WeMTData->ProjectionY("ggMetWeMT_SBloInvMass",sbLoBin1,sbLoBin2,"eo");
  TH1F* ggMetWeMT_SBhiInvMass=(TH1F*)WeMTData->ProjectionY("ggMetWeMT_SBhiInvMass",sbHiBin1,sbHiBin2,"eo");
  TH1F* ggMetWeMT_sigInvMass=(TH1F*)WeMTData->ProjectionY("ggMetWeMT_sigInvMass",sigBin1,sigBin2,"eo");
  /*
  overFlowBin=ggMetWeMT_sigInvMass->FindBin(metPlotXmax+1), lastBin=ggMetWeMT_sigInvMass->FindBin(metPlotXmax-1);
  overFlow=0.;overFlowErr=0.;
  val=ggMetWeMT_sigInvMass->GetBinContent(lastBin);valErr=ggMetWeMT_sigInvMass->GetBinError(lastBin);
  overFlow=ggMetWeMT_sigInvMass->IntegralAndError(overFlowBin,-1,overFlowErr);
  ggMetWeMT_sigInvMass->SetBinContent(lastBin,val+overFlow);
  ggMetWeMT_sigInvMass->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  overFlow=0.;overFlowErr=0.;
  val=ggMetWeMT_SBloInvMass->GetBinContent(lastBin);valErr=ggMetWeMT_SBloInvMass->GetBinError(lastBin);
  overFlow=ggMetWeMT_SBloInvMass->IntegralAndError(overFlowBin,-1,overFlowErr);
  ggMetWeMT_SBloInvMass->SetBinContent(lastBin,val+overFlow);
  ggMetWeMT_SBloInvMass->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  overFlow=0.;overFlowErr=0.;
  val=ggMetWeMT_SBhiInvMass->GetBinContent(lastBin);valErr=ggMetWeMT_SBhiInvMass->GetBinError(lastBin);
  overFlow=ggMetWeMT_SBhiInvMass->IntegralAndError(overFlowBin,-1,overFlowErr);
  ggMetWeMT_SBhiInvMass->SetBinContent(lastBin,val+overFlow);
  ggMetWeMT_SBhiInvMass->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  */
  TH1F* ggMetWeMT_SBloInvMassRebin = (TH1F*)ggMetWeMT_SBloInvMass->Rebin(nMTbins,"ggMetWeMT_SBloInvMassRebin",MTbins);
  TH1F* ggMetWeMT_SBhiInvMassRebin = (TH1F*)ggMetWeMT_SBhiInvMass->Rebin(nMTbins,"ggMetWeMT_SBhiInvMassRebin",MTbins);
  TH1F* ggMetWeMT_sigInvMassRebin = (TH1F*)ggMetWeMT_sigInvMass->Rebin(nMTbins,"ggMetWeMT_sigInvMassRebin",MTbins);
 
  AddOverflowToLastBin(ggMetWeMT_SBloInvMassRebin);
  AddOverflowToLastBin(ggMetWeMT_SBhiInvMassRebin);
  AddOverflowToLastBin(ggMetWeMT_sigInvMassRebin);

  float higgslowscaleWeMT = ggMetWeMT_SBloInvMassRebin->Integral() ? WePowYieldSig/*(HiggSigNoNormWe)*//ggMetWeMT_SBloInvMassRebin->Integral() : 0;
  float higgshighscaleWeMT = ggMetWeMT_SBhiInvMassRebin->Integral() ? WePowYieldSig/*(HiggSigNoNormWe)*//ggMetWeMT_SBhiInvMassRebin->Integral() : 0;
  float WeMTFitSyst=(WePowYieldSigErr*WePowYieldSigErr)/(WePowYieldSig*WePowYieldSig);
  double WeMTbgStatErrLO=0.;double WeMTbgStatLO=ggMetWeMT_SBloInvMassRebin->IntegralAndError(0,-1,WeMTbgStatErrLO);
  float WeMTbgStatSystLO=(WeMTbgStatErrLO*WeMTbgStatErrLO)/(WeMTbgStatLO*WeMTbgStatLO);
  double WeMTbgStatErrHI=0.;double WeMTbgStatHI=ggMetWeMT_SBhiInvMassRebin->IntegralAndError(0,-1,WeMTbgStatErrHI);
  float WeMTbgStatSystHI=(WeMTbgStatErrHI*WeMTbgStatErrHI)/(WeMTbgStatHI*WeMTbgStatHI); 

  TH1F* ggMetWeMT_SBloInvMassClone = (TH1F*)ggMetWeMT_SBloInvMass->Clone("ggMetWeMT_SBloInvMassClone");
  TH1F* ggMetWeMT_SBhiInvMassClone = (TH1F*)ggMetWeMT_SBhiInvMass->Clone("ggMetWeMT_SBhiInvMassClone");
  ggMetWeMT_SBloInvMassClone->Scale(higgslowscale);ggMetWeMT_SBhiInvMassClone->Scale(higgshighscale);
  TH1F* ggMetWeMT_SBloInvMassCloneRebin = (TH1F*)ggMetWeMT_SBloInvMassClone->Rebin(nMTbins,"ggMetWeMT_SBloInvMassCloneRebin",MTbins);
  TH1F* ggMetWeMT_SBhiInvMassCloneRebin = (TH1F*)ggMetWeMT_SBhiInvMassClone->Rebin(nMTbins,"ggMetWeMT_SBhiInvMassCloneRebin",MTbins);
 
  TH1F* ggMetWeMT_InvMassSBcombNoRebin=(TH1F*)ggMetWeMT_SBloInvMass->Clone();//ggMetWeMT_InvMassSBcombNoRebin->Add(ggMetWeMT_SBhiInvMass);ggMetWeMT_InvMassSBcombNoRebin->Scale(0.5);
  TH1F *StatErrsWeMT       =(TH1F*)ggMetWeMT_SBloInvMass->Clone();
  TH1F *FitStatSystErrsWeMT =(TH1F*)ggMetWeMT_SBloInvMass->Clone();
  TH1F *FitShapeSystErrsWeMT  =(TH1F*)ggMetWeMT_SBloInvMass->Clone();
  TH1F *HalfDiffErrsWeMT  =(TH1F*)ggMetWeMT_SBloInvMass->Clone();
  valTot=0.;
  for(int i=1;i<=ggMetWeMT_InvMassSBcombNoRebin->GetNbinsX();i++){
    Double_t val = ((ggMetWeMT_SBloInvMass->GetBinContent(i)*higgslowscaleWeMT)+(ggMetWeMT_SBhiInvMass->GetBinContent(i)*higgshighscaleWeMT))/2.;
    //Double_t val = (ggMetWeMT_SBloInvMass->GetBinContent(i)*higgslowscaleWeMT);
    //double Lerr=0.,Uerr=0.;
    //double Lval = ggMetWeMT_SBloInvMass->IntegralAndError(0,-1,Lerr);
    //double Uval = ggMetWeMT_SBhiInvMass->IntegralAndError(0,-1,Uerr);
    Double_t LvalScaleWeMT = ggMetWeMT_SBloInvMass->GetBinContent(i)*higgslowscaleWeMT,UvalScaleWeMT = ggMetWeMT_SBhiInvMass->GetBinContent(i)*higgshighscaleWeMT;
    Double_t WeMTStatErr(0.),WeMTFitStatSystErr(0.),WeMTFitShapeSystErr(0.),WeMTHalfDiffErr(0.);
    double err = TotErr(val,ggMetWeMT_SBloInvMass->GetBinContent(i),ggMetWeMT_SBloInvMass->GetBinError(i),ggMetWeMT_SBhiInvMass->GetBinContent(i),ggMetWeMT_SBhiInvMass->GetBinError(i),WeMTbgStatLO,WeMTbgStatErrLO,WeMTbgStatHI,WeMTbgStatErrHI,WePowYieldSig,WePowYieldSigErr,WeExpoYieldSig,WeMTStatErr,WeMTFitStatSystErr,WeMTFitShapeSystErr,LvalScaleWeMT,UvalScaleWeMT,WeMTHalfDiffErr);
    ggMetWeMT_InvMassSBcombNoRebin->SetBinContent(i,val);ggMetWeMT_InvMassSBcombNoRebin->SetBinError(i,err);
    StatErrsWeMT->SetBinError(i,WeMTStatErr);FitStatSystErrsWeMT->SetBinError(i,WeMTFitStatSystErr);FitShapeSystErrsWeMT->SetBinError(i,WeMTFitShapeSystErr);HalfDiffErrsWeMT->SetBinError(i,WeMTHalfDiffErr);
    valTot+=val; 
 }
  cout<<"WePowYieldSig: "<<WePowYieldSig<<"  valTot MT: "<<valTot<<endl;
  TH1F* ggMetWeMT_InvMassSBcomb = (TH1F*)ggMetWeMT_InvMassSBcombNoRebin->Rebin(nMTbins,"ggMetWeMT_InvMassSBcomb",MTbins);
  cout<<"ggMetWe_InvMassSBcomb MT after rebin: "<<ggMetWeMT_InvMassSBcomb->Integral()<<" (0,-1): "<<ggMetWeMT_InvMassSBcomb->Integral(0,-1)<<endl;
  for(int i=0;i<ggMetWeMT_sigInvMassRebin->FindBin(31);i++){
    cout<<"WeMT "<<ggMetWeMT_sigInvMassRebin->GetBinLowEdge(i)<<"<"<<ggMetWeMT_sigInvMassRebin->GetBinLowEdge(i+1)<<"  events: "<<ggMetWeMT_sigInvMassRebin->GetBinContent(i)<<" +- "<<ggMetWeMT_sigInvMassRebin->GetBinError(i)<<endl;
  }
  AddOverflowToLastBin(ggMetWeMT_InvMassSBcomb);
  cout<<"ggMetWeMT_InvMassSBcomb after overflow add: "<<ggMetWeMT_InvMassSBcomb->Integral()<<" (0,-1): "<<ggMetWeMT_InvMassSBcomb->Integral(0,-1)<<endl;
  fout_We.cd();
  TH1F* ggMetWeMT_sigInvMassRebin_clone = (TH1F*)ggMetWeMT_sigInvMassRebin->Clone();
  //MakeBlindMet30(ggMetWeMT_sigInvMassRebin_clone);
  ggMetWeMT_sigInvMassRebin_clone->SetTitle("");ggMetWeMT_sigInvMassRebin_clone->GetXaxis()->SetTitle("M_{T} [GeV]");
  ggMetWeMT_sigInvMassRebin_clone->Write("hMT_ggEle_tag");
  fout.cd();
  ggMetWeMT_sigInvMassRebin_clone->Write("hMT_ggEle_tag");
  DivideBy30gev(ggMetWeMT_InvMassSBcomb);
  DivideBy30gev(ggMetWeMT_sigInvMassRebin);
  /*
  ggMetWeMT_SBloInvMassRebin->Scale(higgslowscaleWeMT);
  ggMetWeMT_SBhiInvMassRebin->Scale(higgshighscaleWeMT);

  AddOverflowToLastBin(ggMetWeMT_SBloInvMassRebin);
  AddOverflowToLastBin(ggMetWeMT_SBhiInvMassRebin);
  AddOverflowToLastBin(ggMetWeMT_sigInvMassRebin);
  DivideByBinWidth(ggMetWeMT_SBloInvMassRebin);
  DivideByBinWidth(ggMetWeMT_SBhiInvMassRebin);
  DivideByBinWidth(ggMetWeMT_sigInvMassRebin);
 
  TH1F* ggMetWeMT_InvMassSBcomb=(TH1F*)ggMetWeMT_SBloInvMassRebin->Clone();ggMetWeMT_InvMassSBcomb->Add(ggMetWeMT_SBhiInvMassRebin);ggMetWeMT_InvMassSBcomb->Scale(0.5);
*/
  c1->cd();c1->SetLogy(0);
  ggMetWeMT_SBloInvMassCloneRebin->SetFillColor(kRed);ggMetWeMT_SBloInvMassCloneRebin->SetFillStyle(3004);ggMetWeMT_SBloInvMassCloneRebin->SetMarkerSize(0);
  ggMetWeMT_SBloInvMassCloneRebin->GetXaxis()->SetRangeUser(0,250);
  ggMetWeMT_SBhiInvMassCloneRebin->SetFillColor(kBlue);ggMetWeMT_SBhiInvMassCloneRebin->SetFillStyle(3004);ggMetWeMT_SBhiInvMassCloneRebin->SetMarkerSize(0);
  ggMetWeMT_SBloInvMassCloneRebin->SetLineColor(kRed);ggMetWeMT_SBloInvMassCloneRebin->SetMarkerColor(kRed);
  ggMetWeMT_SBhiInvMassCloneRebin->SetLineColor(kBlue);ggMetWeMT_SBhiInvMassCloneRebin->SetMarkerColor(kBlue);
  ggMetWeMT_SBhiInvMassCloneRebin->SetLineWidth(2);ggMetWeMT_SBloInvMassCloneRebin->SetLineWidth(2);
  ggMetWeMT_sigInvMassRebin_clone->SetTitle("");ggMetWeMT_sigInvMassRebin_clone->GetXaxis()->SetTitle("M_{T}");ggMetWeMT_sigInvMassRebin_clone->GetYaxis()->SetTitle("Events");
  ggMetWeMT_sigInvMassRebin_clone->Draw("PE");
  ggMetWeMT_SBhiInvMassCloneRebin->Draw("E2sames");ggMetWeMT_SBhiInvMassCloneRebin->Draw("pesames");
  ggMetWeMT_SBloInvMassCloneRebin->Draw("E2sames");ggMetWeMT_SBloInvMassCloneRebin->Draw("pesames");
  ggMetWeMT_sigInvMassRebin_clone->Draw("PEsames");
  lowhighleg->Draw();
  c1->Print("Plots/Higgs/Exclusive_WeMT_SbLowAndHigh.png");
  c1->Print("Plots/Higgs/Exclusive_WeMT_SbLowAndHigh.pdf");
  c1->SetLogy(1);

  p1->cd();
  ggMetWeMT_sigInvMassRebin->SetMarkerColor(kBlack);ggMetWeMT_sigInvMassRebin->SetLineColor(kBlack);
  ggMetWeMT_SBloInvMassRebin->SetMarkerColor(kBlue);ggMetWeMT_SBloInvMassRebin->SetLineColor(kBlue);ggMetWeMT_SBloInvMassRebin->SetFillColor(kBlue);ggMetWeMT_SBloInvMassRebin->SetFillStyle(3004);
  ggMetWeMT_SBhiInvMassRebin->SetMarkerColor(kRed);ggMetWeMT_SBhiInvMassRebin->SetLineColor(kRed);ggMetWeMT_SBhiInvMassRebin->SetFillColor(kRed);ggMetWeMT_SBhiInvMassRebin->SetFillStyle(3004);
  ggMetWeMT_InvMassSBcomb->SetMarkerColor(42);ggMetWeMT_InvMassSBcomb->SetLineColor(42);ggMetWeMT_InvMassSBcomb->SetFillColor(42);//ggMetWeMT_InvMassSBcomb->SetFillStyle(3004);
  ggMetWeMT_SBloInvMassRebin->SetMarkerSize(0.75);ggMetWeMT_SBhiInvMassRebin->SetMarkerSize(0.75);ggMetWeMT_InvMassSBcomb->SetMarkerSize(0.75);
  ggMetWeMT_InvMassSBcomb->SetMarkerSize(0);
  WeMTStackRebin->Add(ggMetWeMT_InvMassSBcomb);
  TH1F* ggMetWeMT_InvMassSBcomb_plus_SMhiggs = (TH1F*)ggMetWeMT_InvMassSBcomb->Clone();ggMetWeMT_InvMassSBcomb_plus_SMhiggs->Add(SMHiggsMTWeMT);
  TH1F* ggMetWeMT_InvMassSBcomb_plus_SMhiggs_clone = (TH1F*)ggMetWeMT_InvMassSBcomb_plus_SMhiggs->Clone();
  TH1F* ggMetWeMT_InvMassSBcomb_LimitClone = (TH1F*)ggMetWeMT_InvMassSBcomb->Clone();
  TH1F* SMHiggsMTWeMT_LimitCLone = (TH1F*)SMHiggsMTWeMT->Clone();
  MultBy30gev(ggMetWeMT_InvMassSBcomb_plus_SMhiggs_clone);
  MultBy30gev(ggMetWeMT_InvMassSBcomb_LimitClone);
  MultBy30gev(SMHiggsMTWeMT_LimitCLone);

  fout_We.cd();
  ggMetWeMT_InvMassSBcomb_plus_SMhiggs_clone->SetTitle("");ggMetWeMT_InvMassSBcomb_plus_SMhiggs_clone->GetXaxis()->SetTitle("M_{T} [GeV]");
  ggMetWeMT_InvMassSBcomb_plus_SMhiggs_clone->Write("smBackground");
  fout.cd();
  ggMetWeMT_InvMassSBcomb_LimitClone->Write("hMT_ggEle_bkg");
  fout_SMHiggs.cd();
  SMHiggsMTWeMT_LimitCLone->Write("SMHiggsMET_ggEle");
  ggMetWeMT_InvMassSBcomb_plus_SMhiggs->SetMarkerSize(0.75);
  ggMetWeMT_InvMassSBcomb_plus_SMhiggs->SetMarkerSize(0);ggMetWeMT_InvMassSBcomb_plus_SMhiggs->SetFillColor(kBlack);ggMetWeMT_InvMassSBcomb_plus_SMhiggs->SetFillStyle(3004);

  //const double alpha = 1 - 0.6827;
  TGraphAsymmErrors *ggMetWeMT_sigInvMassRebin_tgraph = new TGraphAsymmErrors(ggMetWeMT_sigInvMassRebin);
  //ggMetWeMT_sigInvMassRebin_tgraph->SetMarkerSize(0.5);
  ggMetWeMT_sigInvMassRebin_tgraph->SetMarkerStyle(20);
  for (int i = 0; i < ggMetWeMT_sigInvMassRebin_tgraph->GetN(); ++i) {
    double N = ggMetWeMT_sigInvMassRebin_tgraph->GetY()[i];
    double L=0.,U=0.;
    if(i!=3){
      L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
      U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
      ggMetWeMT_sigInvMassRebin_tgraph->SetPointEYlow(i, N-L);
      ggMetWeMT_sigInvMassRebin_tgraph->SetPointEYhigh(i, U-N);
    }
    else{
      L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N*3.,1.));
      U =  ROOT::Math::gamma_quantile_c(alpha/2,3.*N+1,1);
      L/=3.;U/=3.;
      ggMetWeMT_sigInvMassRebin_tgraph->SetPointEYlow(i, N-L);
      ggMetWeMT_sigInvMassRebin_tgraph->SetPointEYhigh(i, U-N);
    }
    //      double U =  (N==0) ?  ( ROOT::Math::gamma_quantile_c(alpha,N+1,1) ) :
    //         ( ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) );
    //ggMetWeMT_sigInvMassRebin_tgraph->SetPointEYlow(i, N-L);
    //ggMetWeMT_sigInvMassRebin_tgraph->SetPointEYhigh(i, U-N);
    cout<<"ele i:"<<i<<"  N: "<<N<<"  L: "<<L<<"  N-L: "<<N-L<<"  U: "<<U<<"  U-N: "<<U-N<<endl;
  }

  ggMetWeMT_sigInvMassRebin_tgraph->GetXaxis()->SetRangeUser(0,179);
  ggMetWeMT_sigInvMassRebin_tgraph->GetXaxis()->SetLabelSize(0);
  ggMetWeMT_sigInvMassRebin_tgraph->GetYaxis()->SetRangeUser(9e-3,50);
  ggMetWeMT_sigInvMassRebin_tgraph->GetYaxis()->SetTitle("Events / 30 GeV");
  ggMetWeMT_sigInvMassRebin_tgraph->GetYaxis()->SetTitleSize(0.05);
  ggMetWeMT_sigInvMassRebin_tgraph->SetTitle("");
  ggMetWeMT_sigInvMassRebin_tgraph->SetMarkerSize(1);
  //MakeBlindMet30(ggMetWeMT_sigInvMassRebin);

  TGraphAsymmErrors* ggMetWeMT_sigInvMassRebin_tgraph_clone = (TGraphAsymmErrors*)ggMetWeMT_sigInvMassRebin_tgraph->Clone("ggMetWeMT_sigInvMassRebin_tgraph_clone");

  ggMetWeMT_sigInvMassRebin_tgraph_clone->SetMarkerSize(2);
  ggMetWeMT_sigInvMassRebin_tgraph_clone->SetLineWidth(3);
  for(int i=0;i<ggMetWeMT_sigInvMassRebin_tgraph_clone->GetN();++i){
    ggMetWeMT_sigInvMassRebin_tgraph_clone->SetPointEXlow(i,0);
    ggMetWeMT_sigInvMassRebin_tgraph_clone->SetPointEXhigh(i,0);
  }
  ggMetWeMT_sigInvMassRebin_tgraph_clone->Draw("APE");
  WeMTStackRebin->Draw("histoSAMES");
  //MakeBlindMet30(ggMetWeMT_InvMassSBcomb_plus_SMhiggs);
  ggMetWeMT_InvMassSBcomb_plus_SMhiggs->Draw("E2SAMES");
  ggMetWeMT_sigInvMassRebin_tgraph_clone->Draw("PESAMES");
  h_SMS_WH_gg_1Ele_MT_Excluded_Rebin->SetLineStyle(7);
  h_SMS_ZH_gg_1Ele_MT_Excluded_Rebin->SetLineStyle(9);
  TH1F* h_SMS_HH_all_1Ele_MT_Excluded_Rebin = (TH1F*)h_SMS_HH_2W2g_gg_1Ele_MT_Excluded_Rebin->Clone();h_SMS_HH_all_1Ele_MT_Excluded_Rebin->Add(h_SMS_HH_2Z2g_gg_1Ele_MT_Excluded_Rebin);h_SMS_HH_all_1Ele_MT_Excluded_Rebin->Add(h_SMS_HH_2tau2g_gg_1Ele_MT_Excluded_Rebin);
  h_SMS_HH_all_1Ele_MT_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_ZH_gg_1Ele_MT_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_WH_gg_1Ele_MT_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_HH_2b2g_gg_1Ele_MT_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_HH_2W2g_gg_1Ele_MT_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_HH_2Z2g_gg_1Ele_MT_Excluded_Rebin->Draw("histoSAMES");
  ggMetWeMT_sigInvMassRebin_tgraph_clone->Draw("PESAMES");
  p1->RedrawAxis();
  PrelimText->Draw();
  /*
  ggTextWe->Draw();
  WeMetLeg->Draw();
  WeMetLeg2->Draw();
  */
  ggTextWeLog->Draw();
  WeMetLegLog->Draw();
  WeMetLeg2Log->Draw();
  textSMSlog->Draw();
  p2->cd();
  TH1F* ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs = (TH1F*)ggMetWeMT_sigInvMassRebin->Clone();
  //ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs->Divide(ggMetWeMT_InvMassSBcomb_plus_SMhiggs);
  TH1F* h_WeSystErrMT = (TH1F*)ggMetWeMT_InvMassSBcomb_plus_SMhiggs->Clone();h_WeSystErrMT->SetFillColor(kBlack);h_WeSystErrMT->SetFillStyle(3004);h_WeSystErrMT->SetMarkerSize(0);
  for(int i=1;i<ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs->GetNbinsX()+1;i++){
    float Value = ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs->GetBinContent(i);float StatErr = ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs->GetBinError(i);
    Value/=ggMetWeMT_InvMassSBcomb_plus_SMhiggs->GetBinContent(i);StatErr/=ggMetWeMT_InvMassSBcomb_plus_SMhiggs->GetBinContent(i);
    ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs->SetBinContent(i,Value);ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs->SetBinError(i,StatErr);
    float SystErr=h_WeSystErrMT->GetBinError(i)/h_WeSystErrMT->GetBinContent(i);
    h_WeSystErrMT->SetBinContent(i,1); h_WeSystErrMT->SetBinError(i,SystErr);
  }

  TGraphAsymmErrors* ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph = new TGraphAsymmErrors(ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs);
  for (int i = 0; i < ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetN(); ++i) {

    float N = ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetY()[i];
    double U = ggMetWeMT_sigInvMassRebin_tgraph->GetErrorYhigh(i);
    double L = ggMetWeMT_sigInvMassRebin_tgraph->GetErrorYlow(i);
    double Utemp=U,Ltemp=L;
    U/=ggMetWeMT_InvMassSBcomb_plus_SMhiggs->GetBinContent(i+1);L/=ggMetWeMT_InvMassSBcomb_plus_SMhiggs->GetBinContent(i+1);
    ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetPointEYlow(i, L);
    ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetPointEYhigh(i, U);
    if(i==0)cout<<"Ratio:"<<endl;
    cout<<"i:"<<i<<"  N: "<<N<<"  Ltemp: "<<Ltemp<<"Lbin: "<<ggMetWeMT_InvMassSBcomb_plus_SMhiggs->GetBinContent(i+1)<<"  L: "<<L<<"  Utemp: "<<Utemp<<"  Ubin: "<<ggMetWeMT_InvMassSBcomb_plus_SMhiggs->GetBinContent(i+1)<<" U: "<<U<<endl;
  }


  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetLineColor(kBlack);//ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetMarkerColor(kBlack);ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetMarkerSize(0.75);
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetTitle("");
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetTitle("");
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetRangeUser(0.,6);
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetTitle("");
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetTitle("#frac{Data}{Prediction}");
  //ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetTitle("M_{T} [GeV]");
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetTitle("");
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetTitleOffset(0.4);
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetTitleSize(0.14);
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetTitleSize(0.2);
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetLabelSize(0.12);
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetLabelSize(.13);
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetNdivisions(206,0);

  //MakeBlindMet30(ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs);
  //MakeBlindMet30(h_WeSystErrMT);
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetRangeUser(0,179);
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetMarkerSize(2);
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetLineWidth(3);
  for(int i=0;i<ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetN();++i){
    ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetPointEXlow(i,0);
    ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetPointEXhigh(i,0);
  }
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Draw("APE");
  h_WeSystErrMT->Draw("E2SAMES");
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Draw("PESAMES");
  l1.DrawLine(0,1,180,1);
  //PrelimText->Draw();
  p3->cd();p3->Clear();
  TPaveText *mtText = new TPaveText(.77,0,.97,1,"NDC");
  mtText->AddText("M_{T} [GeV]");
  mtText->SetFillColor(0);mtText->SetTextFont(42);
  mtText->SetBorderSize(0);
  mtText->SetTextSize(.6);
  mtText->Draw();
  c2->Print("Plots/Higgs/Exclusive_WeDataMTWithRatio.png");
  c2->Print("Plots/Higgs/Exclusive_WeDataMTWithRatio.pdf");
  p1->cd();p1->SetLogy(0);
  ggMetWeMT_sigInvMassRebin_tgraph_clone->GetYaxis()->SetRangeUser(0,14);
  ggMetWeMT_sigInvMassRebin_tgraph_clone->GetYaxis()->SetNdivisions(207,0);
  ggMetWeMT_sigInvMassRebin_tgraph_clone->SetMarkerSize(2);
  ggMetWeMT_sigInvMassRebin_tgraph_clone->Draw("APE");
  WeMTStackRebin->Draw("histoSAMES");
  ggMetWeMT_InvMassSBcomb_plus_SMhiggs->Draw("E2SAMES");
  ggMetWeMT_sigInvMassRebin_tgraph_clone->Draw("PESAMES");
  h_SMS_HH_all_1Ele_MT_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_ZH_gg_1Ele_MT_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_WH_gg_1Ele_MT_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_HH_2W2g_gg_1Ele_MT_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_HH_2Z2g_gg_1Ele_MT_Excluded_Rebin->Draw("histoSAMES");
  ggMetWeMT_sigInvMassRebin_tgraph_clone->Draw("PESAMES");
  p1->RedrawAxis();
  PrelimText->Draw();
  ggTextWe->Draw();
  WeMetLeg->Draw();
  WeMetLeg2->Draw();
  textSMS->Draw();
  TLine metVertErr2(104.13,12.,104.13,13.2);metVertErr2.SetLineWidth(3);
  metVertErr2.Draw();
  p2->cd();
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetMarkerSize(2);
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Draw("APE");
  h_WeSystErrMT->Draw("E2SAMES");
  ggMetWeMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Draw("PESAMES");
  l1.DrawLine(0,1,180,1);
  //PrelimText->Draw();
  p3->cd();p3->Clear();
  mtText->Draw();
  c2->Print("Plots/Higgs/Exclusive_WeDataMTWithRatio_linear.png");
  c2->Print("Plots/Higgs/Exclusive_WeDataMTWithRatio_linear.pdf");
  p1->SetLogy(1);

  //W->muon nu -- 1 muon, <2 jets
  c1->cd();
  TH2F* WmuData = (TH2F*)fin->Get("ggMetVsInvarMass_Loose_1Mu_0_1Jets");
  TH2F* WmuggHgg = (TH2F*)f_ggHgg->Get("ggMetVsInvarMass_Loose_1Mu_0_1Jets");
  TH2F* WmuWZHgg = (TH2F*)f_WZHgg->Get("ggMetVsInvarMass_Loose_1Mu_0_1Jets");
  TH2F* WmuTTHgg = (TH2F*)f_TTHgg->Get("ggMetVsInvarMass_Loose_1Mu_0_1Jets");
  TH2F* WmuVBFHgg = (TH2F*)f_VBFHgg->Get("ggMetVsInvarMass_Loose_1Mu_0_1Jets");

  TH1F* met_Phi_mu_ggH = (TH1F*)f_ggHgg->Get("met_Phi_ggLoose_Mu");
  TH1F* met_Phi_mu_WZH = (TH1F*)f_WZHgg->Get("met_Phi_ggLoose_Mu");
  TH1F* met_Phi_mu_TTH = (TH1F*)f_TTHgg->Get("met_Phi_ggLoose_Mu");
  TH1F* met_Phi_mu_VBFH = (TH1F*)f_VBFHgg->Get("met_Phi_ggLoose_Ele");

  float ggHnormSF_mu = met_Phi_mu_ggH->GetEntries()/met_Phi_mu_ggH->Integral(0,-1);
  float WZHnormSF_mu = met_Phi_mu_WZH->GetEntries()/met_Phi_mu_WZH->Integral(0,-1);
  float TTHnormSF_mu = met_Phi_mu_TTH->GetEntries()/met_Phi_mu_TTH->Integral(0,-1);
  float VBFHnormSF_mu = met_Phi_mu_VBFH->GetEntries()/met_Phi_mu_VBFH->Integral(0,-1);

  //TH2F* Wmu_aaW130 = (TH2F*)f_aaW130->Get("ggMetVsInvarMass_Loose_1Mu_0_1Jets");
  //TH2F* Wmu_aaW275 = (TH2F*)f_aaW275->Get("ggMetVsInvarMass_Loose_1Mu_0_1Jets");
  WmuggHgg->Scale((ggHnormSF_mu*PhoEffScale2*L_int*2.29e-03*19.22)/99989.);//125GeV=/96290);  
  WmuVBFHgg->Scale((VBFHnormSF_mu*PhoEffScale2*L_int*2.29e-03*1.544)/95677.);//125GeV=/99885); 
  WmuTTHgg->Scale((TTHnormSF_mu*PhoEffScale2*L_int*2.29e-03*.1271)/100048.);//125GeV=/100224);
  WmuWZHgg->Scale((WZHnormSF_mu*PhoEffScale2*L_int*2.29e-03*(.6782/**(.3257+.014)*/+.3843/**.2*/))/100320);//125 and 126 GeV have same # events
  WmuggHgg->SetLineColor(kGreen);WmuggHgg->SetMarkerColor(kGreen);WmuggHgg->SetFillColor(kGreen);
  WmuWZHgg->SetLineColor(kCyan);WmuWZHgg->SetMarkerColor(kCyan);WmuWZHgg->SetFillColor(kCyan);
  WmuTTHgg->SetLineColor(31);WmuTTHgg->SetMarkerColor(31);WmuTTHgg->SetFillColor(31);
  WmuVBFHgg->SetLineColor(kRed+3);WmuVBFHgg->SetMarkerColor(kRed+3);WmuVBFHgg->SetFillColor(kRed+3);
  //WmuggHgg->SetFillStyle(0);WmuWZHgg->SetFillStyle(0);WmuTTHgg->SetFillStyle(0);WmuVBFHgg->SetFillStyle(0);
  WmuggHgg->SetLineWidth(2);WmuWZHgg->SetLineWidth(2);WmuVBFHgg->SetLineWidth(2);WmuTTHgg->SetLineWidth(2);
  /*Wmu_aaW130->Scale(scale_aaW130);
  Wmu_aaW130->SetLineColor(kBlue);Wmu_aaW130->SetLineWidth(3);
  Wmu_aaW275->Scale(scale_aaW275);
  Wmu_aaW275->SetLineColor(kRed);Wmu_aaW275->SetLineWidth(3);*/
  TH1D* WmuDataProjX = (TH1D*)WmuData->ProjectionX("WmuDataProjX",0,-1,"eo");
  TH1D* WmuggHggProjX = (TH1D*)WmuggHgg->ProjectionX("WmuggHggProjX",0,-1,"eo");
  TH1D* WmuWZHggProjX = (TH1D*)WmuWZHgg->ProjectionX("WmuWZHggProjX",0,-1,"eo");
  TH1D* WmuTTHggProjX = (TH1D*)WmuTTHgg->ProjectionX("WmuTTHggProjX",0,-1,"eo");
  TH1D* WmuVBFHggProjX = (TH1D*)WmuVBFHgg->ProjectionX("WmuVBFHggProjX",0,-1,"eo");
  //TH1D* WmuAAW130ProjX = (TH1D*)Wmu_aaW130->ProjectionX("WmuAAW130ProjX",0,-1,"eo");
  //TH1D* WmuAAW275ProjX = (TH1D*)Wmu_aaW275->ProjectionX("WmuAAW275ProjX",0,-1,"eo");
  TH1D* WmuDataProjY = (TH1D*)WmuData->ProjectionY("WmuDataProjY",0,-1,"eo");
  TH1D* WmuggHggProjY = (TH1D*)WmuggHgg->ProjectionY("WmuggHggProjY",0,-1,"eo");
  TH1D* WmuWZHggProjY = (TH1D*)WmuWZHgg->ProjectionY("WmuWZHggProjY",0,-1,"eo");
  TH1D* WmuTTHggProjY = (TH1D*)WmuTTHgg->ProjectionY("WmuTTHggProjY",0,-1,"eo");
  TH1D* WmuVBFHggProjY = (TH1D*)WmuVBFHgg->ProjectionY("WmuVBFHggProjY",0,-1,"eo");
  //TH1D* WmuAAW130ProjY = (TH1D*)Wmu_aaW130->ProjectionY("WmuAAW130ProjY",0,-1,"eo");
  //TH1D* WmuAAW275ProjY = (TH1D*)Wmu_aaW275->ProjectionY("WmuAAW275ProjY",0,-1,"eo");
  WmuDataProjX->Rebin(2);WmuggHggProjX->Rebin(2);WmuTTHggProjX->Rebin(2);WmuWZHggProjX->Rebin(2);WmuVBFHggProjX->Rebin(2);
  
  int WmuBin103 = WmuDataProjX->FindBin(102.9),WmuBin118 = WmuDataProjX->FindBin(117.9),WmuBin120 = WmuDataProjX->FindBin(119.9),WmuBin131 = WmuDataProjX->FindBin(130.9),WmuBin133 = WmuDataProjX->FindBin(132.9),WmuBin163 = WmuDataProjX->FindBin(162.9);
    cout<<"Wmu events invmass less than 103: "<<WmuDataProjX->Integral(0,WmuBin103)<<" 103<invmass<163: "<<WmuDataProjX->Integral(WmuBin103+1,WmuBin163)<<" invmass<163: "<<WmuDataProjX->Integral(WmuBin163+1,-1)<<" 103<invmass<118: "<<WmuDataProjX->Integral(WmuBin103+1,WmuBin118)<<" 120<invmass<131: "<<WmuDataProjX->Integral(WmuBin120+1,WmuBin131)*0<<" 133<invmass<163: "<<WmuDataProjX->Integral(WmuBin133+1,WmuBin163)<<endl;
    //MakeBlindInvMass(WmuDataProjX);
  c1->SetLogy(1);
  WmuDataProjY->GetXaxis()->SetRangeUser(0,200);WmuDataProjY->Draw();
  c1->Print("Plots/Higgs/Exclusive_WmuDataMet.png");

 
  THStack *WmuMetStack = new THStack("WmuMetStack","");
  WmuMetStack->Add(WmuTTHggProjY);WmuMetStack->Add(WmuVBFHggProjY);WmuMetStack->Add(WmuWZHggProjY);WmuMetStack->Add(WmuggHggProjY);
  WmuMetStack->Draw("histo");WmuMetStack->GetXaxis()->SetRangeUser(0,249.9);
  double SMHiggsMetIntWmu_ggH=0.,SMHiggsMetIntErrWmu_ggH=0.,SMHiggsMetIntWmu_WZH=0.,SMHiggsMetIntErrWmu_WZH=0.,SMHiggsMetIntWmu_TTH=0.,SMHiggsMetIntErrWmu_TTH=0.,SMHiggsMetIntWmu_VBFH=0.,SMHiggsMetIntErrWmu_VBFH=0.;
  SMHiggsMetIntWmu_ggH=WmuggHggProjY->IntegralAndError(0,-1,SMHiggsMetIntErrWmu_ggH);
  SMHiggsMetIntWmu_WZH=WmuWZHggProjY->IntegralAndError(0,-1,SMHiggsMetIntErrWmu_WZH);
  SMHiggsMetIntWmu_TTH=WmuTTHggProjY->IntegralAndError(0,-1,SMHiggsMetIntErrWmu_TTH);
  SMHiggsMetIntWmu_VBFH=WmuVBFHggProjY->IntegralAndError(0,-1,SMHiggsMetIntErrWmu_VBFH);
  //float SMHiggsMetIntWmu=WmuTTHggProjY->Integral()+WmuVBFHggProjY->Integral()+WmuWZHggProjY->Integral()+WmuggHggProjY->Integral();
  float SMHiggsMetIntWmu=SMHiggsMetIntWmu_ggH+SMHiggsMetIntWmu_WZH+SMHiggsMetIntWmu_TTH+SMHiggsMetIntWmu_VBFH;
  float SMHiggsMetIntErrWmu=sqrt(SMHiggsMetIntErrWmu_ggH*SMHiggsMetIntErrWmu_ggH+SMHiggsMetIntErrWmu_WZH*SMHiggsMetIntErrWmu_WZH+SMHiggsMetIntErrWmu_TTH*SMHiggsMetIntErrWmu_TTH+SMHiggsMetIntErrWmu_VBFH*SMHiggsMetIntErrWmu_VBFH);
  //float SMHiggsMetIntWmu=WmuTTHggProjY->Integral()+WmuVBFHggProjY->Integral()+WmuWZHggProjY->Integral()+WmuggHggProjY->Integral();
  //float SMHiggsMetIntWmu=WmuTTHggProjY->Integral()+WmuVBFHggProjY->Integral()+WmuWZHggProjY->Integral()+WmuggHggProjY->Integral();
  c1->Print("Plots/Higgs/Exclusive_WmuSMHiggsMet.png");
  c1->Print("Plots/Higgs/Exclusive_WmuSMHiggsMet.pdf");
  TH1F* WmuggHggProjYRebin=(TH1F*)WmuggHggProjY->Rebin(NmetBinsXtraWide,"WmuggHggProjYRebin",xbinsXtraWide);
  TH1F* WmuTTHggProjYRebin=(TH1F*)WmuTTHggProjY->Rebin(NmetBinsXtraWide,"WmuTTHggProjYRebin",xbinsXtraWide);
  TH1F* WmuWZHggProjYRebin=(TH1F*)WmuWZHggProjY->Rebin(NmetBinsXtraWide,"WmuWZHggProjYRebin",xbinsXtraWide);
  TH1F* WmuVBFHggProjYRebin=(TH1F*)WmuVBFHggProjY->Rebin(NmetBinsXtraWide,"WmuVBFHggProjYRebin",xbinsXtraWide);
  //TH1F* WmuAAW130ProjYRebin=(TH1F*)WmuAAW130ProjY->Rebin(NmetBinsXtraWide,"WmuAAW130ProjYRebin",xbinsXtraWide);
  //TH1F* WmuAAW275ProjYRebin=(TH1F*)WmuAAW275ProjY->Rebin(NmetBinsXtraWide,"WmuAAW275ProjYRebin",xbinsXtraWide);
  AddOverflowToLastBin(WmuggHggProjYRebin);
  AddOverflowToLastBin(WmuTTHggProjYRebin);
  AddOverflowToLastBin(WmuWZHggProjYRebin);
  AddOverflowToLastBin(WmuVBFHggProjYRebin);
  AddHiggsSigmaUncert(WmuggHggProjYRebin);
  AddHiggsSigmaUncert(WmuTTHggProjYRebin);
  AddHiggsSigmaUncert(WmuWZHggProjYRebin);
  AddHiggsSigmaUncert(WmuVBFHggProjYRebin);
  DivideBy15gev(WmuggHggProjYRebin);
  DivideBy15gev(WmuTTHggProjYRebin);
  DivideBy15gev(WmuWZHggProjYRebin);
  DivideBy15gev(WmuVBFHggProjYRebin);

  THStack *WmuMetStackRebin = new THStack("WmuMetStackRebin","");
  THStack *WmuMetStackRebin2 = new THStack("WmuMetStackRebin2","");
  WmuTTHggProjYRebin->SetLineColor(kBlack);WmuVBFHggProjYRebin->SetLineColor(kBlack);WmuWZHggProjYRebin->SetLineColor(kBlack);WmuggHggProjYRebin->SetLineColor(kBlack);
  WmuTTHggProjYRebin->SetLineWidth(1);WmuVBFHggProjYRebin->SetLineWidth(1);WmuWZHggProjYRebin->SetLineWidth(1);WmuggHggProjYRebin->SetLineWidth(1);
  TH1F* SMHiggsAdditionForWmu = (TH1F*)WmuTTHggProjYRebin->Clone();SMHiggsAdditionForWmu->Add(WmuVBFHggProjYRebin);SMHiggsAdditionForWmu->Add(WmuWZHggProjYRebin);SMHiggsAdditionForWmu->Add(WmuggHggProjYRebin);
  WmuMetStackRebin2->Add(WmuTTHggProjYRebin);WmuMetStackRebin2->Add(WmuVBFHggProjYRebin);WmuMetStackRebin2->Add(WmuWZHggProjYRebin);WmuMetStackRebin2->Add(WmuggHggProjYRebin);
  WmuMetStackRebin->Add(SMHiggsAdditionForWmu);
  WmuMetStackRebin2->Draw("histo");WmuMetStackRebin2->GetXaxis()->SetRangeUser(0,249.9);
  WmuMetStackRebin2->SetMaximum(.4);WmuMetStackRebin2->SetMinimum(3.1E-3);
  WmuMetStackRebin2->GetYaxis()->SetTitle("Events / 15 GeV");WmuMetStackRebin2->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  legHiggs2->Draw();
  c1->Print("Plots/Higgs/Exclusive_WmuSMHiggsMetRebin.png");
  c1->Print("Plots/Higgs/Exclusive_WmuSMHiggsMetRebin.pdf");
  TH1F* SMHiggsMetWmu = (TH1F*)WmuTTHggProjYRebin->Clone();SMHiggsMetWmu->Add(WmuVBFHggProjYRebin);SMHiggsMetWmu->Add(WmuggHggProjYRebin);SMHiggsMetWmu->Add(WmuWZHggProjYRebin);
  c1->SetLogy(0);
  THStack *WmuInvMassStack = new THStack("WmuInvMassStack","");
  WmuInvMassStack->Add(WmuTTHggProjX);WmuInvMassStack->Add(WmuVBFHggProjX);WmuInvMassStack->Add(WmuWZHggProjX);WmuInvMassStack->Add(WmuggHggProjX);
  WmuInvMassStack->Draw("histo");WmuInvMassStack->GetXaxis()->SetRangeUser(110,139.9);
  c1->Print("Plots/Higgs/Exclusive_WmuSMHiggsInvMass.png");
  cout<<"SM higgs Expected Yield, Wmunu : "<<SMHiggsMetIntWmu<<" +- "<<SMHiggsMetIntErrWmu<<endl;
  WmuDataProjX->GetXaxis()->SetRangeUser(20,269);WmuDataProjX->SetTitle("");
  WmuDataProjX->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");WmuDataProjX->GetYaxis()->SetTitle("Events / GeV");WmuDataProjX->GetYaxis()->SetTitleOffset(0.7);
  WmuDataProjX->Draw();
  PrelimTextInvMass->Draw();
  TPaveText *ggTextWmuInvMass = new TPaveText(.29,.78,.6,.89,"NDC");
  ggTextWmuInvMass->AddText("#gamma#gamma + #mu");
  ////ggTextWmu->AddText("");
  ggTextWmuInvMass->SetFillStyle(0);ggTextWmuInvMass->SetTextFont(42);
  ggTextWmuInvMass->SetFillColor(0);
  ggTextWmuInvMass->SetBorderSize(0);
  ggTextWmuInvMass->Draw();
  c1->Print("Plots/Higgs/Exclusive_WmuDataInvMass.png");
  c1->Print("Plots/Higgs/Exclusive_WmuDataInvMass.pdf");
  //invariant mass roofit
  RooRealVar xMetWmu_bg("xMetWmu_bg","m_{#gamma#gamma}",95,200,"GeV");
  xMetWmu_bg.setRange("sb_loWmu_bg",100,sbLoHi);xMetWmu_bg.setRange("sb_hiWmu_bg",sbHiLo,sbHiHi);xMetWmu_bg.setRange("fullWmu_bg",sbLoLo,sbHiHi);xMetWmu_bg.setRange("sigWmu_bg",sigLo,sigHi);
  RooDataHist dataMetWmu_bg("dataMetWmu_bg","dataset_bg",xMetWmu_bg,WmuDataProjX);
  RooRealVar lambdaWmu("lambdaWmu","lambdaWmu",-1.,-10.,5.);
  RooRealVar PolWmu1_bg("PolWmu1_bg","PolWmu1_bg",0,-.04,.04);
  RooRealVar PolWmu2_bg("PolWmu2_bg","PolWmu2_bg",0,-.000001,.000001);
  RooRealVar PolWmu3_bg("PolWmu3_bg","PolWmu3_bg",0,-.0000001,.0000001);
  RooRealVar PolWmu4_bg("PolWmu4_bg","PolWmu4_bg",0,-.00000001,.00000001);
  RooRealVar PolWmu5_bg("PolWmu5_bg","PolWmu5_bg",0,-.0000000001,.0000000001);
  //RooPolynomial PolWmu_bg("PolWmu_bg","5th order Polynomial",xMetWmu_bg,RooArgList(PolWmu1_bg,PolWmu2_bg,PolWmu3_bg,PolWmu4_bg,PolWmu5_bg));
  RooExponential PolWmu_bg("ExpoWmu_bg","ExpoWmu_bg",xMetWmu_bg,lambdaWmu);
  RooRealVar WmuYield_bg("WmuYield_bg","WmuYield_bg",70000,0,170000);
  RooAddPdf MetWmuPdf_bg("MetWmuPdf_bg","MetWmuPdf_bg",RooArgList(PolWmu_bg),RooArgList(WmuYield_bg));
  //RooFitResult *rNoMetCutWmu_bg = MetWmuPdf_bg.chi2FitTo(dataMetWmu_bg/*,Extended(kTRUE)*/,Save(),Range("sb_lo,sb_hi"));
  RooFitResult *rNoMetCutWmu_bg = MetWmuPdf_bg.fitTo(dataMetWmu_bg,Extended(kTRUE),Save(),Range("sb_loWmu_bg,sb_hiWmu_bg"));
  RooArgList pars_Wmu(*MetWmuPdf_bg.getParameters(RooArgSet(xMetWmu_bg) ) );
  RooArgSet prodSet_Wmu(MetWmuPdf_bg); //prodSet.add(nsig);
  RooProduct unNormPdf_Wmu("fitted Function", "fitted Function", prodSet_Wmu);
  TF1 * fit_Wmu = unNormPdf_Wmu.asTF(RooArgList(xMetWmu_bg), pars_Wmu, RooArgList(xMetWmu_bg));
  float nsig_Wmu = ((RooRealVar*) pars_Wmu.find("WmuYield_bg"))->getVal();
  Double_t integ_full_Wmu = fit_Wmu->Integral(sbLoLo, sbHiHi);
  Double_t integ_Wmu = nsig_Wmu*fit_Wmu->Integral(sigLo, sigHi, 0)/integ_full_Wmu;
  Double_t dinteg_Wmu = nsig_Wmu*fit_Wmu->IntegralError(sigLo, sigHi, 0, rNoMetCutWmu_bg->covarianceMatrix().GetMatrixArray())/integ_full_Wmu;
  RooPlot *xframeMetWmu_bg = xMetWmu_bg.frame(Title("5th Order Polynomial Background, gg Wmu with no MET cut"));
  MetWmuPdf_bg.paramOn(xframeMetWmu_bg, Format("NE",AutoPrecision(2)),Layout(0.6,0.995,0.9) );
  dataMetWmu_bg.plotOn(xframeMetWmu_bg,LineColor(kBlack),MarkerStyle(20),MarkerSize(0.3));
  //MetWmuPdf_bg.plotOn(xframeMetWmu_bg,Components(PolWmu_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCutWmu_bg,3,kTRUE),FillColor(kViolet));
  MetWmuPdf_bg.plotOn(xframeMetWmu_bg,Components(PolWmu_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCutWmu_bg,2,kTRUE),FillColor(kGreen));
  MetWmuPdf_bg.plotOn(xframeMetWmu_bg,Components(PolWmu_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCutWmu_bg,1,kTRUE),FillColor(kOrange));
  MetWmuPdf_bg.plotOn(xframeMetWmu_bg,Components(PolWmu_bg),LineColor(kRed),LineStyle(kDashed),Range(102,180));
  MetWmuPdf_bg.plotOn(xframeMetWmu_bg,Components(PolWmu_bg),LineColor(kRed));
  dataMetWmu_bg.plotOn(xframeMetWmu_bg,LineColor(kBlack),MarkerStyle(20),MarkerSize(0.3));
  //xframeMetWmu_bg->SetAxisRange(106,149,"X");
  //xframeMetWmu_bg->SetAxisRange(0,1900,"Y");
  //xframeMetWmu_bg->SetMaximum(2000);
  xframeMetWmu_bg->Draw();
  //line125->DrawLine(125.3,0,125.3,xframeMetWmu_bg->GetMaximum());
  WmuInvMassStack->Draw("histoSAME");
  xframeMetWmu_bg->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/Higgs/Exclusive_WmuDataInvMassFit.png");

  reject=true;
  TF1* fitCurveWmu = new TF1("mgg_Wmu_fit",fpow,sbFitLoLo,sbFitHiHi,2);
  avg_l = WmuDataProjX->Integral(WmuDataProjX->FindBin(100.),WmuDataProjX->FindBin(sbLoHi))/20.,avg_u = WmuDataProjX->Integral(WmuDataProjX->FindBin(sbHiLo),WmuDataProjX->FindBin(sbFitHiHi))/69.5,avgX_l=110.,avgX_u=148.;
  cout<<avg_l<<"  "<<avg_u<<"  "<<avgX_l<<"  "<<avgX_u<<endl;
  cout<<WmuDataProjX->Integral()<<"  "<<WmuDataProjX->Integral(0,sbHiLo)<<"  "<<WmuDataProjX->Integral(sbHiLo,sbFitHiHi)<<"  "<<WmuDataProjX->Integral(sbHiLo,sbLoHi)<<"  "<<WmuDataProjX->Integral(sbHiLo,sbFitHiHi)<<endl;
  Double_t param1Wmu= (log(avg_l) - log(avg_u))/(log(avgX_l) - log(avgX_u));
  Double_t param0Wmu= avg_l/pow(avgX_l, param1Wmu);
  cout<<"param0Wmu: "<<param0Wmu<<"  param1Wmu: "<<param1Wmu<<endl;
  fitCurveWmu->SetParameter(0,param0Wmu);
  fitCurveWmu->SetParameter(1,param1Wmu);
  //first just to check the status
  status = WmuDataProjX->Fit(fitCurveWmu,"L","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResultWmu = WmuDataProjX->Fit(fitCurveWmu,"SLLMEV0","",sbFitLoLo,sbFitHiHi);
  TMatrixDSym covWmu = fitResultWmu->GetCovarianceMatrix();
  fitResultWmu->Print("V");
  WmuDataProjX->GetXaxis()->SetRangeUser(100,164.9);
  reject=false;
  YieldBinWidth=WmuDataProjX->GetBinWidth(1);
  float WmuPowYield = fitCurveWmu->Integral(sbLoLo,sbHiHi)/YieldBinWidth;
  float WmuPowYieldErr = fitCurveWmu->IntegralError(sbLoLo,sbHiHi,fitResultWmu->GetParams(),covWmu.GetMatrixArray() )/YieldBinWidth; 
  float WmuPowYieldSig = fitCurveWmu->Integral(sigLo,sigHi)/YieldBinWidth;
  float WmuPowYieldSigErr = fitCurveWmu->IntegralError(sigLo,sigHi,fitResultWmu->GetParams(),covWmu.GetMatrixArray() )/YieldBinWidth; 
  Double_t WmuPowSBloYield = WmuDataProjX->Integral(WmuDataProjX->FindBin(sbLoLo),WmuDataProjX->FindBin(sbLoHi-.1));
  Double_t WmuPowSBloYieldFit = fitCurveWmu->Integral(sbLoLo,sbLoHi)/YieldBinWidth;
  Double_t WmuPowSBloYieldFitErr = fitCurveWmu->IntegralError(sbLoLo,sbLoHi,fitResultWmu->GetParams(),covWmu.GetMatrixArray() )/YieldBinWidth; 
  cout<<"Wmu low sideband yield from histo: "<<WmuPowSBloYield<<"  and from fit: "<<WmuPowSBloYieldFit<<" +- "<<WmuPowSBloYieldFitErr<<endl; 
  cout<<"Wmu higgs window yield from fit: "<<WmuPowYieldSig<<" +- "<<WmuPowYieldSigErr<<endl; 
  Double_t WmuPowSBhiYield = WmuDataProjX->Integral(WmuDataProjX->FindBin(sbHiLo),WmuDataProjX->FindBin(sbHiHi-.1));
  Double_t WmuPowSBhiYieldFit = fitCurveWmu->Integral(sbHiLo,sbHiHi)/YieldBinWidth;
  Double_t WmuPowSBhiYieldFitErr = fitCurveWmu->IntegralError(sbHiLo,sbHiHi,fitResultWmu->GetParams(),covWmu.GetMatrixArray() )/YieldBinWidth; 
  cout<<"Wmu high sideband yield from histo: "<<WmuPowSBhiYield<<"  and from fit: "<<WmuPowSBhiYieldFit<<" +- "<<WmuPowSBhiYieldFitErr<<endl;
  reject=true;
  char strWmu[50],strWmuSig[50];
  sprintf(strWmu,"Fit Yield: %4.2f #pm %4.2f",WmuPowYield,WmuPowYieldErr);
  sprintf(strWmuSig,"Fit Signal Yield: %4.2f #pm %4.2f",WmuPowYieldSig,WmuPowYieldSigErr);
  char strWmuRF[50];
  sprintf(strWmuRF,"RooFit signal Yield: %4.2f #pm %4.2f",integ_Wmu,dinteg_Wmu);
  TPaveText *textWmu= new TPaveText(.5,.65,.84,.84,"NDC");textWmu->SetFillStyle(0);textWmu->SetBorderSize(0);
  //textWmu->AddText(strWmu);
  textWmu->AddText(strWmuSig);textWmu->SetTextFont(42);
  //textWmu->AddText(strWmuRF);
  WmuDataProjX->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  WmuDataProjX->GetYaxis()->SetTitle("Events / GeV");
  WmuDataProjX->SetTitle("");
  WmuDataProjX->GetYaxis()->SetRangeUser(0,6);WmuDataProjX->GetYaxis()->SetNdivisions(206);
  WmuDataProjX->GetYaxis()->SetTitleOffset(0.8);
  WmuDataProjX->GetXaxis()->SetTitleOffset(1.2);
  WmuDataProjX->SetTitleFont(42,"xy");
  WmuDataProjX->SetLabelFont(42,"xy");
  WmuDataProjX->GetXaxis()->SetLabelOffset(0.02);
  WmuDataProjX->GetYaxis()->SetLabelOffset(0.02);  
  WmuDataProjX->GetListOfFunctions()->Remove(fitCurveWmu);
  WmuDataProjX->SetMarkerSize(1.5);
  WmuDataProjX->Draw();
  h_SMS_HH_2W2g_1Mu_Excluded_InvMass->Draw("histosames");
  TF1 *fitCurveWmuNew = new TF1("fitCurveWmuNew","[0]*pow(x,[1])",sbLoLo,sbHiHi);fitCurveWmuNew->SetParameters(fitCurveWmu->GetParameters());
  fitCurveWmuNew->SetLineColor(kBlue);fitCurveWmuNew->SetLineWidth(3);fitCurveWmuNew->Draw("SAMES");
  /*
  TF1 *fitCurveWmu2 = new TF1("fitCurveWmu2","[0]*pow(x,[1])",sbLoHi,sbHiLo);
  fitCurveWmu2->SetParameters(fitCurveWmu->GetParameters());
  //fitCurveWmu2->SetLineColor(kBlue);fitCurveWmu2->Draw("SAMES");
  TF1 *fitCurveWmuLeft = new TF1("fitCurveWmuLeft","[0]*pow(x,[1])",sbLoLo,sbLoHi);fitCurveWmuLeft->SetParameters(fitCurveWmu->GetParameters());
  WmuDataProjX->GetListOfFunctions()->Add(fitCurveWmuLeft);//gROOT->GetListOfFunctions()->Remove(fitCurveWmuLeft);
  //fitCurveWmuLeft->SetLineColor(kRed);fitCurveWmuLeft->Draw("SAMES");
  TF1 *fitCurveWmuRight = new TF1("fitCurveWmuRight","[0]*pow(x,[1])",sbHiLo,sbHiHi);fitCurveWmuRight->SetParameters(fitCurveWmu->GetParameters());
  WmuDataProjX->GetListOfFunctions()->Add(fitCurveWmuRight);//gROOT->GetListOfFunctions()->Remove(fitCurveWmuRight);
  //fitCurveWmuRight->SetLineColor(kRed);fitCurveWmuRight->Draw("SAMES");
  */
  WmuDataProjX->Draw("SAMES");//textWmu->Draw("SAMES");
  TPaveText *ggTextWmuLog = new TPaveText(.29,.78,.6,.89,"NDC");
  ggTextWmuLog->AddText("#gamma#gamma + #mu");
  ////ggTextWmu->AddText("");
  ggTextWmuInvMass->SetTextSize(.08);
  ggTextWmuInvMass->SetTextFont(42);
  ggTextWmuLog->SetFillStyle(0);
  ggTextWmuLog->SetFillColor(0);
  ggTextWmuLog->SetBorderSize(0);
  //ggTextWmuLog->Draw();
  ggTextWmuInvMass->Draw();
  TLegend *WmuInvMassLeg = new TLegend(.56,.59,.85,.89,"","brNDC");
  WmuInvMassLeg->SetFillColor(kWhite);WmuInvMassLeg->SetFillStyle(0);WmuInvMassLeg->SetBorderSize(0);
  WmuInvMassLeg->SetTextSize(.06);WmuInvMassLeg->SetTextFont(42);
  WmuInvMassLeg->AddEntry(WmuDataProjX,"Data","EP");
  WmuInvMassLeg->AddEntry(fitCurveWmuNew,"Sideband Fit","l");
  WmuInvMassLeg->AddEntry(h_SMS_HH_2W2g_1Mu_Excluded_InvMass,"Signal hh","f");
  WmuInvMassLeg->AddEntry(h_SMS_HH_2W2g_1Mu_Excluded_InvMass,"m_{#tilde{#chi}_{1}^{0}}=130 GeV","");
  //WmuInvMassLeg->AddEntry(fitCurveWmu2,"Fit Extrapolation","l");
  WmuInvMassLeg->Draw();
  PrelimTextInvMass->Draw();
  vertErr.Draw();
  c1->Print("Plots/Higgs/Exclusive_WmuDataInvMassFit_powLaw.png");
  c1->Print("Plots/Higgs/Exclusive_WmuDataInvMassFit_powLaw.pdf");

  reject=true;
  TF1* fitCurveExpoWmu = new TF1("fitCurveExpoWmu",fexp,sbFitLoLo,sbFitHiHi,2);
  fitCurveExpoWmu->SetParameter(0,-2.);
  //first just to check the status
  status = WmuDataProjX->Fit(fitCurveExpoWmu,"L","",sbFitLoLo,sbFitHiHi);
  printf("fit status: %i \n",status);
  //Then to get the result
  TFitResultPtr fitResultExpoWmu = WmuDataProjX->Fit(fitCurveExpoWmu,"SL","",sbFitLoLo,sbFitHiHi);
  TMatrixDSym covExpoWmu = fitResultExpoWmu->GetCovarianceMatrix();
  fitResultExpoWmu->Print("V");
  WmuDataProjX->GetXaxis()->SetRangeUser(100,164.9);
  reject=false;
  YieldBinWidth=WeDataProjXexpo->GetBinWidth(1);
  Double_t WmuExpoYield = fitCurveExpoWmu->Integral(sbLoLo,sbHiHi)/YieldBinWidth;
  Double_t WmuExpoYieldErr = fitCurveExpoWmu->IntegralError(sbLoLo,sbHiHi,fitResultExpoWmu->GetParams(),covExpoWmu.GetMatrixArray() )/YieldBinWidth; 
  Double_t WmuExpoYieldSig = fitCurveExpoWmu->Integral(sigLo,sigHi)/YieldBinWidth;
  Double_t WmuExpoYieldSigErr = fitCurveExpoWmu->IntegralError(sigLo,sigHi,fitResultExpoWmu->GetParams(),covExpoWmu.GetMatrixArray() )/YieldBinWidth; 
  reject=true;
  sprintf(strWmu,"Fit Yield: %4.2f #pm %4.2f",WmuExpoYield,WmuExpoYieldErr);
  sprintf(strWmuSig,"Fit Signal Yield: %4.2f #pm %4.2f",WmuExpoYieldSig,WmuExpoYieldSigErr);
  sprintf(strWmuRF,"RooFit signal Yield: %4.2f #pm %4.2f",integ_Wmu,dinteg_Wmu);
  TPaveText *textExpoWmu= new TPaveText(.5,.65,.84,.84,"NDC");textExpoWmu->SetFillStyle(0);textExpoWmu->SetBorderSize(0);
  //textExpoWmu->AddText(strWmu);
  textExpoWmu->AddText(strWmuSig);
  WmuDataProjX->Draw();textExpoWmu->SetTextFont(42);
  textExpoWmu->Draw("SAMES");
  c1->Print("Plots/Higgs/Exclusive_WmuDataInvMassFit_expo.png");

  //met plot with ratio
  RooAbsReal* igWmuSig = MetWmuPdf_bg.createIntegral(xMetWmu_bg,NormSet(xMetWmu_bg),Range("sigWmu_bg"));
  RooAbsReal* igWmuSBlo = MetWmuPdf_bg.createIntegral(xMetWmu_bg,NormSet(xMetWmu_bg),Range("sb_loWmu_bg"));
  RooAbsReal* igWmuSBhi = MetWmuPdf_bg.createIntegral(xMetWmu_bg,NormSet(xMetWmu_bg),Range("sb_hiWmu_bg"));
  float HiggSigNoNormWmu=igWmuSig->getVal()*WmuYield_bg.getVal();
  float HiggLowSBNoNormWmu=igWmuSBlo->getVal()*WmuYield_bg.getVal(),HiggHighSBNoNormWmu=igWmuSBhi->getVal()*WmuYield_bg.getVal();
  //p1->cd();
  sbLoBin1=WmuData->GetXaxis()->FindBin(sbLoLo);sbLoBin2=WmuData->GetXaxis()->FindBin(sbLoHi)-1;  
  sbHiBin1=WmuData->GetXaxis()->FindBin(sbHiLo);sbHiBin2=WmuData->GetXaxis()->FindBin(sbHiHi)-1;
  sigBin1=WmuData->GetXaxis()->FindBin(sigLo);sigBin2=WmuData->GetXaxis()->FindBin(sigHi)-1;
  TH1F* ggMetWmu_SBloInvMass=(TH1F*)WmuData->ProjectionY("ggMetWmu_SBloInvMass",sbLoBin1,sbLoBin2,"eo");
  TH1F* ggMetWmu_SBhiInvMass=(TH1F*)WmuData->ProjectionY("ggMetWmu_SBhiInvMass",sbHiBin1,sbHiBin2,"eo");
  TH1F* ggMetWmu_sigInvMass=(TH1F*)WmuData->ProjectionY("ggMetWmu_sigInvMass",sigBin1,sigBin2,"eo");

  TH1F* ggMetWmu_SBloInvMassRebin = (TH1F*)ggMetWmu_SBloInvMass->Rebin(NmetBinsXtraWide,"ggMetWmu_SBloInvMassRebin",xbinsXtraWide);
  TH1F* ggMetWmu_SBhiInvMassRebin = (TH1F*)ggMetWmu_SBhiInvMass->Rebin(NmetBinsXtraWide,"ggMetWmu_SBhiInvMassRebin",xbinsXtraWide);
  TH1F* ggMetWmu_sigInvMassRebin = (TH1F*)ggMetWmu_sigInvMass->Rebin(NmetBinsXtraWide,"ggMetWmu_sigInvMassRebin",xbinsXtraWide);

  AddOverflowToLastBin(ggMetWmu_SBloInvMassRebin);
  AddOverflowToLastBin(ggMetWmu_SBhiInvMassRebin);
  AddOverflowToLastBin(ggMetWmu_sigInvMassRebin);

  float higgslowscaleWmu = ggMetWmu_SBloInvMassRebin->Integral() ? WmuPowYieldSig/*HiggSigNoNormWmu*//ggMetWmu_SBloInvMassRebin->Integral() : 0;
  float higgshighscaleWmu = ggMetWmu_SBhiInvMassRebin->Integral() ? WmuPowYieldSig/*HiggSigNoNormWmu*//ggMetWmu_SBhiInvMassRebin->Integral() : 0;
  //float WmuFitSyst=(WmuYield_bg.getError()*WmuYield_bg.getError())/(WmuYield_bg.getVal()*WmuYield_bg.getVal());
  float WmuFitSyst=(WmuPowYieldSigErr*WmuPowYieldSigErr)/(WmuPowYieldSig*WmuPowYieldSig);
  double WmubgStatErrLO=0.;double WmubgStatLO=ggMetWmu_SBloInvMassRebin->IntegralAndError(0,-1,WmubgStatErrLO);
  float WmubgStatSystLO=(WmubgStatErrLO*WmubgStatErrLO)/(WmubgStatLO*WmubgStatLO);
  double WmubgStatErrHI=0.;double WmubgStatHI=ggMetWmu_SBhiInvMassRebin->IntegralAndError(0,-1,WmubgStatErrHI);
  float WmubgStatSystHI=(WmubgStatErrHI*WmubgStatErrHI)/(WmubgStatHI*WmubgStatHI);

  TH1F* ggMetWmu_SBloInvMassClone = (TH1F*)ggMetWmu_SBloInvMass->Clone("ggMetWmu_SBloInvMassClone");
  TH1F* ggMetWmu_SBhiInvMassClone = (TH1F*)ggMetWmu_SBhiInvMass->Clone("ggMetWmu_SBhiInvMassClone");
  ggMetWmu_SBloInvMassClone->Scale(higgslowscale);ggMetWmu_SBhiInvMassClone->Scale(higgshighscale);
  TH1F* ggMetWmu_SBloInvMassCloneRebin = (TH1F*)ggMetWmu_SBloInvMassClone->Rebin(NmetBinsXtraWide,"ggMetWmu_SBloInvMassCloneRebin",xbinsXtraWide);
  TH1F* ggMetWmu_SBhiInvMassCloneRebin = (TH1F*)ggMetWmu_SBhiInvMassClone->Rebin(NmetBinsXtraWide,"ggMetWmu_SBhiInvMassCloneRebin",xbinsXtraWide);

  TH1F* ggMetWmu_InvMassSBcombNoRebin=(TH1F*)ggMetWmu_SBloInvMass->Clone();//ggMetWmu_InvMassSBcomb->Add(ggMetWmu_SBhiInvMassRebin);ggMetWmu_InvMassSBcomb->Scale(0.5);
  TH1F *StatErrsWmu       =(TH1F*)ggMetWmu_SBloInvMass->Clone();
  TH1F *FitStatSystErrsWmu =(TH1F*)ggMetWmu_SBloInvMass->Clone();
  TH1F *FitShapeSystErrsWmu  =(TH1F*)ggMetWmu_SBloInvMass->Clone();
  TH1F *HalfDiffErrsWmu  =(TH1F*)ggMetWmu_SBloInvMass->Clone();
  for(int i=1;i<=ggMetWmu_InvMassSBcombNoRebin->GetNbinsX();i++){
    Double_t val = ((ggMetWmu_SBloInvMass->GetBinContent(i)*higgslowscaleWmu)+(ggMetWmu_SBhiInvMass->GetBinContent(i)*higgshighscaleWmu))/2.;
    //double Lerr=0.,Uerr=0.;
    //double Lval = ggMetWmu_SBloInvMass->IntegralAndError(0,-1,Lerr);
    //double Uval = ggMetWmu_SBhiInvMass->IntegralAndError(0,-1,Uerr);
    Double_t LvalScaleWmu = ggMetWmu_SBloInvMass->GetBinContent(i)*higgslowscaleWmu,UvalScaleWmu = ggMetWmu_SBhiInvMass->GetBinContent(i)*higgshighscaleWmu;
    Double_t WmuStatErr(0.),WmuFitStatSystErr(0.),WmuFitShapeSystErr(0.),WmuHalfDiffErr(0.);
    double err = TotErr(val,ggMetWmu_SBloInvMass->GetBinContent(i),ggMetWmu_SBloInvMass->GetBinError(i),ggMetWmu_SBhiInvMass->GetBinContent(i),ggMetWmu_SBhiInvMass->GetBinError(i),WmubgStatLO,WmubgStatErrLO,WmubgStatHI,WmubgStatErrHI,WmuPowYieldSig,WmuPowYieldSigErr,WmuExpoYieldSig,WmuStatErr,WmuFitStatSystErr,WmuFitShapeSystErr,LvalScaleWmu,UvalScaleWmu,WmuHalfDiffErr);
    ggMetWmu_InvMassSBcombNoRebin->SetBinContent(i,val);ggMetWmu_InvMassSBcombNoRebin->SetBinError(i,err);
    StatErrsWmu->SetBinError(i,WmuStatErr);FitStatSystErrsWmu->SetBinError(i,WmuFitStatSystErr);FitShapeSystErrsWmu->SetBinError(i,WmuFitShapeSystErr);HalfDiffErrsWmu->SetBinError(i,WmuHalfDiffErr);
  }
  TH1F* ggMetWmu_InvMassSBcomb = (TH1F*)ggMetWmu_InvMassSBcombNoRebin->Rebin(NmetBinsXtraWide,"ggMetWmu_InvMassSBcomb",xbinsXtraWide);
  AddOverflowToLastBin(ggMetWmu_InvMassSBcomb);
  TH1F* ggMetWmu_sigInvMassRebin_Clone = (TH1F*)ggMetWmu_sigInvMassRebin->Clone("ggMetWmu_sigInvMassRebin_Clone");
  DivideBy15gev(ggMetWmu_InvMassSBcomb);
  DivideBy15gev(ggMetWmu_sigInvMassRebin);

  c1->cd();c1->SetLogy(0);
  ggMetWmu_SBloInvMassCloneRebin->SetFillColor(kRed);ggMetWmu_SBloInvMassCloneRebin->SetFillStyle(3004);ggMetWmu_SBloInvMassCloneRebin->SetMarkerSize(0);
  ggMetWmu_SBloInvMassCloneRebin->GetXaxis()->SetRangeUser(0,250);
  ggMetWmu_SBhiInvMassCloneRebin->SetFillColor(kBlue);ggMetWmu_SBhiInvMassCloneRebin->SetFillStyle(3004);ggMetWmu_SBhiInvMassCloneRebin->SetMarkerSize(0);
  ggMetWmu_SBloInvMassCloneRebin->SetLineColor(kRed);ggMetWmu_SBloInvMassCloneRebin->SetMarkerColor(kRed);
  ggMetWmu_SBhiInvMassCloneRebin->SetLineColor(kBlue);ggMetWmu_SBhiInvMassCloneRebin->SetMarkerColor(kBlue);
  ggMetWmu_SBhiInvMassCloneRebin->SetLineWidth(2);ggMetWmu_SBloInvMassCloneRebin->SetLineWidth(2);
  ggMetWmu_SBhiInvMassCloneRebin->SetTitle("");ggMetWmu_SBhiInvMassCloneRebin->GetXaxis()->SetTitle("E_{T}^{miss}");ggMetWmu_SBhiInvMassCloneRebin->GetYaxis()->SetTitle("Events / 15 GeV");
  ggMetWmu_SBhiInvMassCloneRebin->Draw("E2");ggMetWmu_SBhiInvMassCloneRebin->Draw("pesames");
  ggMetWmu_SBloInvMassCloneRebin->Draw("E2sames");ggMetWmu_SBloInvMassCloneRebin->Draw("pesames");
  ggMetWmu_sigInvMassRebin_Clone->Draw("PESAMES");
  lowhighleg->Draw();
  c1->Print("Plots/Higgs/Exclusive_WmuMet_SbLowAndHigh.png");
  c1->Print("Plots/Higgs/Exclusive_WmuMet_SbLowAndHigh.pdf");
  c1->SetLogy(1);
  p1->cd();

  ggMetWmu_sigInvMassRebin->SetMarkerColor(kBlack);ggMetWmu_sigInvMassRebin->SetLineColor(kBlack);
  ggMetWmu_SBloInvMassRebin->SetMarkerColor(kBlue);ggMetWmu_SBloInvMassRebin->SetLineColor(kBlue);ggMetWmu_SBloInvMassRebin->SetFillColor(kBlue);ggMetWmu_SBloInvMassRebin->SetFillStyle(3004);
  ggMetWmu_SBhiInvMassRebin->SetMarkerColor(kRed);ggMetWmu_SBhiInvMassRebin->SetLineColor(kRed);ggMetWmu_SBhiInvMassRebin->SetFillColor(kRed);ggMetWmu_SBhiInvMassRebin->SetFillStyle(3004);
  ggMetWmu_InvMassSBcomb->SetMarkerColor(42);ggMetWmu_InvMassSBcomb->SetLineColor(42);ggMetWmu_InvMassSBcomb->SetFillColor(42);//ggMetWmu_InvMassSBcomb->SetFillStyle(3004);
  ggMetWmu_SBloInvMassRebin->SetMarkerSize(0.75);ggMetWmu_SBhiInvMassRebin->SetMarkerSize(0.75);ggMetWmu_InvMassSBcomb->SetMarkerSize(0.75);
  WmuMetStackRebin->Add(ggMetWmu_InvMassSBcomb);
  TH1F* ggMetWmu_InvMassSBcomb_plus_SMhiggs = (TH1F*)ggMetWmu_InvMassSBcomb->Clone();ggMetWmu_InvMassSBcomb_plus_SMhiggs->Add(SMHiggsMetWmu);
  ggMetWmu_InvMassSBcomb_plus_SMhiggs->SetMarkerSize(0);ggMetWmu_InvMassSBcomb_plus_SMhiggs->SetFillColor(kBlack);ggMetWmu_InvMassSBcomb_plus_SMhiggs->SetFillStyle(3004);


  TGraphAsymmErrors *ggMetWmu_sigInvMassRebin_tgraph = new TGraphAsymmErrors(ggMetWmu_sigInvMassRebin);
  //ggMetWmu_sigInvMassRebin_tgraph->SetMarkerSize(0.5);
  ggMetWmu_sigInvMassRebin_tgraph->SetMarkerStyle(20);
  for (int i = 0; i < ggMetWmu_sigInvMassRebin_tgraph->GetN(); ++i) {
    double N = ggMetWmu_sigInvMassRebin_tgraph->GetY()[i];
    double L=0.,U=0.;
    double w = ggMetWmu_sigInvMassRebin->GetBinWidth(i+1);
    L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N/(15./w),1.));
    U =  ROOT::Math::gamma_quantile_c(alpha/2,(N/(15./w))+1,1);
    L/=w/15.;U/=w/15.;
    ggMetWmu_sigInvMassRebin_tgraph->SetPointEYlow(i, N-L);
    ggMetWmu_sigInvMassRebin_tgraph->SetPointEYhigh(i, U-N);
    
    //      double U =  (N==0) ?  ( ROOT::Math::gamma_quantile_c(alpha,N+1,1) ) :
    //         ( ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) );
    //ggMetWmuMT_sigInvMassRebin_tgraph->SetPointEYlow(i, N-L);
    //ggMetWmuMT_sigInvMassRebin_tgraph->SetPointEYhigh(i, U-N);
    cout<<"mu i met:"<<i<<"  N: "<<N<<"  L: "<<L<<"  N-L: "<<N-L<<"  U: "<<U<<"  U-N: "<<U-N<<"  w: "<<w<<endl;
  }

  ggMetWmu_sigInvMassRebin_tgraph->GetXaxis()->SetRangeUser(0,179);
  ggMetWmu_sigInvMassRebin_tgraph->GetXaxis()->SetLabelSize(0);
  ggMetWmu_sigInvMassRebin_tgraph->GetYaxis()->SetRangeUser(1e-2,70);
  ggMetWmu_sigInvMassRebin_tgraph->GetYaxis()->SetTitle("Events / 15 GeV");
  ggMetWmu_sigInvMassRebin_tgraph->GetYaxis()->SetTitleSize(0.05);
  ggMetWmu_sigInvMassRebin_tgraph->SetTitle("");
  cout<<"Wmu signal region full integral: "<<ggMetWmu_sigInvMassRebin->Integral(0,-1)<<endl;
  //MakeBlindMet30(ggMetWmu_sigInvMassRebin);
  double WmuSMest=0.,WmuSMestErr=0.;WmuSMest=ggMetWmu_InvMassSBcomb_plus_SMhiggs->IntegralAndError(0,-1,WmuSMestErr);
  cout<<"Wmu background estimate full integral: "<<WmuSMest<<"  and error: "<<WmuSMestErr<<endl;
  cout<<"Wmu SM higgs full integral: "<<SMHiggsMetWmu->Integral()<<endl;
  //cout<<"Wmu signal mX=130 full integral: "<<WmuAAW130ProjYRebin->Integral(0,-1)<<endl;
  //cout<<"Wmu signal mX=275 full integral: "<<WmuAAW275ProjYRebin->Integral(0,-1)<<endl;

  TGraphAsymmErrors* ggMetWmu_sigInvMassRebin_tgraph_clone = (TGraphAsymmErrors*)ggMetWmu_sigInvMassRebin_tgraph->Clone("ggMetWmu_sigInvMassRebin_tgraph_clone");

  ggMetWmu_sigInvMassRebin_tgraph_clone->SetMarkerSize(2);
  ggMetWmu_sigInvMassRebin_tgraph_clone->SetLineWidth(3);
  for(int i=0;i<ggMetWmu_sigInvMassRebin_tgraph_clone->GetN();++i){
    ggMetWmu_sigInvMassRebin_tgraph_clone->SetPointEXlow(i,0);
    ggMetWmu_sigInvMassRebin_tgraph_clone->SetPointEXhigh(i,0);
  }
  ggMetWmu_sigInvMassRebin_tgraph_clone->Draw("APE");
  WmuMetStackRebin->Draw("histoSAMES");
  //MakeBlindMet30(ggMetWmu_InvMassSBcomb_plus_SMhiggs);
  ggMetWmu_InvMassSBcomb_plus_SMhiggs->Draw("E2SAMES");
  ggMetWmu_sigInvMassRebin_tgraph_clone->Draw("PESAMES");
  //WmuAAW130ProjYRebin->Draw("histoSAMES");
  //WmuAAW275ProjYRebin->Draw("histoSAMES");
  h_SMS_WH_gg_1Mu_Excluded_Rebin->SetLineStyle(7);
  h_SMS_ZH_gg_1Mu_Excluded_Rebin->SetLineStyle(9);
  TH1F* h_SMS_HH_all_1Mu_Excluded_Rebin = (TH1F*)h_SMS_HH_2W2g_gg_1Mu_Excluded_Rebin->Clone();h_SMS_HH_all_1Mu_Excluded_Rebin->Add(h_SMS_HH_2Z2g_gg_1Mu_Excluded_Rebin);h_SMS_HH_all_1Mu_Excluded_Rebin->Add(h_SMS_HH_2tau2g_gg_1Mu_Excluded_Rebin);
  h_SMS_HH_all_1Mu_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_WH_gg_1Mu_notExcluded_Rebin->Draw("histoSAMES");
  h_SMS_ZH_gg_1Mu_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_WH_gg_1Mu_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_HH_2Z2g_gg_1Mu_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_HH_2W2g_gg_1Mu_Excluded_Rebin->Draw("histoSAMES");
  ggMetWmu_sigInvMassRebin_tgraph_clone->Draw("PESAMES");
  p1->RedrawAxis();
  PrelimText->Draw();
  //TPaveText *ggTextWmu = new TPaveText(.45,.78,.65,.83,"NDC");
  TPaveText *ggTextWmu = new TPaveText(.37,.77,.57,.88,"NDC");
  //ggTextWmu->AddText("#gamma#gamma + #mu, #leq1 b-Jet");
  ggTextWmu->AddText("#gamma#gamma + #mu");
  ////ggTextWmu->AddText("");
  ggTextWmu->SetFillStyle(0);ggTextWmu->SetTextFont(42);
  ggTextWmu->SetFillColor(0);
  ggTextWmu->SetBorderSize(0);
  //ggTextWmu->Draw();
  //TLegend *WmuMetLeg = new TLegend(.44,.33,.81,.79,"","brNDC");
  //TLegend *WmuMetLeg = new TLegend(.52,.57,.92,.78,"","brNDC");
  ggMetWmu_sigInvMassRebin->SetMarkerSize(2);
  TLegend *WmuMetLeg = new TLegend(.57,.3,.97,.88,"","brNDC");
  WmuMetLeg->SetFillColor(kWhite);WmuMetLeg->SetTextSize(.05);WmuMetLeg->SetTextFont(42);
  WmuMetLeg->AddEntry(ggMetWmu_sigInvMassRebin,"Data","epz");
  WmuMetLeg->AddEntry(ggMetWmu_InvMassSBcomb,"Non-higgs SM bg","f");
  WmuMetLeg->AddEntry(SMHiggsAdditionForWmu,"SM higgs","f");
  WmuMetLeg->AddEntry(h_SMS_WH_gg_1Mu_Excluded_Rebin,   "Signal","");
  WmuMetLeg->AddEntry(h_SMS_WH_gg_1Mu_Excluded_Rebin,   "hW, m_{#tilde{#chi}_{1}^{#pm}}=130 GeV","l");
  WmuMetLeg->AddEntry(h_SMS_HH_all_1Mu_Excluded_Rebin,"hh,  m_{#tilde{#chi}_{1}^{0}}=130 GeV","l");
  WmuMetLeg->AddEntry(h_SMS_ZH_gg_1Mu_Excluded_Rebin,"hZ,  m_{#tilde{#chi}_{1}^{0}}=130 GeV","l");
  /*
    WmuMetLeg->AddEntry(WmuWZHggProjYRebin,"SM W/ZH, m_{h}=126 GeV","f");
    WmuMetLeg->AddEntry(WmuVBFHggProjYRebin,"SM VBFH, m_{h}=126 GeV","f");
    WmuMetLeg->AddEntry(WmuTTHggProjYRebin,"SM TTH, m_{h}=126 GeV","f");
    WmuMetLeg->AddEntry(WmuggHggProjYRebin,"SM gg->H, m_{h}=126 GeV","f");
  */		
  /*						       
  WmuMetLeg->AddEntry(WmuWZHggProjYRebin,"SM W/ZH","f");
  WmuMetLeg->AddEntry(WmuVBFHggProjYRebin,"SM VBFH","f");
  WmuMetLeg->AddEntry(WmuTTHggProjYRebin,"SM TTH","f");
  WmuMetLeg->AddEntry(WmuggHggProjYRebin,"SM gg->H","f");
  */
  //WmuMetLeg->AddEntry(ggMetWmu_InvMassSBcomb_plus_SMhiggs,"Full SM background error","f");
  //WmuMetLeg->AddEntry(WmuAAW130ProjYRebin,"AAW m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  //WmuMetLeg->AddEntry(WmuAAW275ProjYRebin,"AAW m_{#tilde{H}}=275, m_{#tilde{B}}=1","l");
  //WmuMetLeg->AddEntry(h_SMS_WH_gg_1Mu_Excluded_Rebin,   "SMS WH+ZH m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  //WmuMetLeg->AddEntry(h_SMS_WH_gg_1Mu_notExcluded_Rebin,"SMS WH+ZH m_{#tilde{H}}=400, m_{#tilde{B}}=150","l");
  /*
  WmuMetLeg->AddEntry(h_SMS_WH_gg_1Mu_Excluded_Rebin,"SMS WH m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  WmuMetLeg->AddEntry(h_SMS_HH_2W2g_gg_1Mu_Excluded_Rebin,"SMS H(WW)H m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  WmuMetLeg->AddEntry(h_SMS_ZH_gg_1Mu_Excluded_Rebin,"SMS ZH m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  WmuMetLeg->AddEntry(h_SMS_HH_2Z2g_gg_1Mu_Excluded_Rebin,"SMS H(ZZ)H m_{#tilde{H}}=130, m_{#tilde{B}}=1","l");
  */
  /*
  WmuMetLeg->AddEntry(h_SMS_WH_gg_1Mu_Excluded_Rebin,"Signal hW","l");
  WmuMetLeg->AddEntry(h_SMS_HH_2W2g_gg_1Mu_Excluded_Rebin,"SMS H(WW)H","l");
  WmuMetLeg->AddEntry(h_SMS_ZH_gg_1Mu_Excluded_Rebin,"Signal hZ","l");
  WmuMetLeg->AddEntry(h_SMS_HH_2Z2g_gg_1Mu_Excluded_Rebin,"SMS H(ZZ)H","l");
  */
  WmuMetLeg->SetFillStyle(0);WmuMetLeg->SetBorderSize(0);
  //WmuMetLeg->Draw();
  TLegend *WmuMetLeg2 = new TLegend(.53,.4,.89,.57,"","brNDC");
  WmuMetLeg2->SetFillColor(kWhite);WmuMetLeg2->SetTextFont(42);
  WmuMetLeg2->SetNColumns(2);
  /*
  WmuMetLeg2->AddEntry(WmuWZHggProjYRebin,"SM W/ZH","f");
  WmuMetLeg2->AddEntry(WmuVBFHggProjYRebin,"SM VBFH","f");
  WmuMetLeg2->AddEntry(WmuTTHggProjYRebin,"SM TTH","f");
  WmuMetLeg2->AddEntry(WmuggHggProjYRebin,"SM gg->H","f");
  */
  ////WmuMetLeg2->AddEntry(h_SMS_WH_gg_1Mu_Excluded_Rebin,"Signal hW","l");
  //WmuMetLeg2->AddEntry(h_SMS_HH_2W2g_gg_1Mu_Excluded_Rebin,"SMS H(WW)H","l");
  ////WmuMetLeg2->AddEntry(h_SMS_HH_all_1Mu_Excluded_Rebin,"Signal hh","l");
  ////WmuMetLeg2->AddEntry(h_SMS_ZH_gg_1Mu_Excluded_Rebin,"Signal hZ","l");
  //WmuMetLeg2->AddEntry(h_SMS_HH_2Z2g_gg_1Mu_Excluded_Rebin,"SMS H(ZZ)H","l");
  WmuMetLeg2->SetFillStyle(0);WmuMetLeg2->SetBorderSize(0);
  //WmuMetLeg2->Draw();

  ggTextWmuLog->Draw();
  ggMetWmu_sigInvMassRebin->SetMarkerSize(2);
  TLegend *WmuMetLegLog = new TLegend(.59,.66,.93,.91,"","brNDC");
  WmuMetLegLog->SetFillColor(kWhite);WmuMetLegLog->SetTextFont(42);
  WmuMetLegLog->AddEntry(ggMetWmu_sigInvMassRebin,"Data","elpz");
  WmuMetLegLog->AddEntry(ggMetWmu_InvMassSBcomb,"Non-higgs SM bg","f");
  WmuMetLegLog->AddEntry(SMHiggsAdditionForWmu,"SM higgs","f");
  WmuMetLegLog->SetFillStyle(0);WmuMetLegLog->SetBorderSize(0);
  WmuMetLegLog->Draw();
  TLegend *WmuMetLeg2Log = new TLegend(.595,.495,.94,.675,"","brNDC");
  WmuMetLeg2Log->SetFillColor(kWhite);WmuMetLeg2Log->SetTextFont(42);
  WmuMetLeg2Log->SetNColumns(2);
  WmuMetLeg2Log->AddEntry(h_SMS_WH_gg_1Mu_Excluded_Rebin,"Signal hW","l");
  WmuMetLeg2Log->AddEntry(h_SMS_HH_all_1Mu_Excluded_Rebin,"Signal hh","l");
  WmuMetLeg2Log->AddEntry(h_SMS_ZH_gg_1Mu_Excluded_Rebin,"Signal hZ","l");
  WmuMetLeg2Log->SetFillStyle(0);WmuMetLeg2Log->SetBorderSize(0);
  WmuMetLeg2Log->Draw();
  textSMSlog->Draw();
  p2->cd();
  TH1F* ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs = (TH1F*)ggMetWmu_sigInvMassRebin->Clone();
  //ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs->Divide(ggMetWmu_InvMassSBcomb_plus_SMhiggs);
  TH1F* h_SystErrWmu = (TH1F*)ggMetWmu_InvMassSBcomb_plus_SMhiggs->Clone();h_SystErrWmu->SetFillColor(kBlack);h_SystErrWmu->SetFillStyle(3004);h_SystErrWmu->SetMarkerSize(0);
  for(int i=1;i<ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs->GetNbinsX()+1;i++){
    float Value = ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs->GetBinContent(i);float StatErr = ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs->GetBinError(i);
    Value/=ggMetWmu_InvMassSBcomb_plus_SMhiggs->GetBinContent(i);StatErr/=ggMetWmu_InvMassSBcomb_plus_SMhiggs->GetBinContent(i);
    ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs->SetBinContent(i,Value);ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs->SetBinError(i,StatErr);
    float SystErr=h_SystErrWmu->GetBinError(i)/h_SystErrWmu->GetBinContent(i);
    h_SystErrWmu->SetBinContent(i,1); h_SystErrWmu->SetBinError(i,SystErr);
  }

  TGraphAsymmErrors* ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph = new TGraphAsymmErrors(ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs);
  for (int i = 0; i < ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetN(); ++i) {

    float N = ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetY()[i];
    double U = ggMetWmu_sigInvMassRebin_tgraph->GetErrorYhigh(i);
    double L = ggMetWmu_sigInvMassRebin_tgraph->GetErrorYlow(i);
    double Utemp=U,Ltemp=L;
    U/=ggMetWmu_InvMassSBcomb_plus_SMhiggs->GetBinContent(i+1);L/=ggMetWmu_InvMassSBcomb_plus_SMhiggs->GetBinContent(i+1);
    ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetPointEYlow(i, L);
    ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetPointEYhigh(i, U);
    if(i==0)cout<<"Ratio:"<<endl;
    cout<<"i:"<<i<<"  N: "<<N<<"  Ltemp: "<<Ltemp<<"Lbin: "<<ggMetWmu_InvMassSBcomb_plus_SMhiggs->GetBinContent(i+1)<<"  L: "<<L<<"  Utemp: "<<Utemp<<"  Ubin: "<<ggMetWmu_InvMassSBcomb_plus_SMhiggs->GetBinContent(i+1)<<" U: "<<U<<endl;
  }



  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetLineColor(kBlack);//ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs->SetMarkerColor(kBlack);//ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs->SetMarkerSize(0.75);
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetTitle("");
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetTitle("");
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetRangeUser(0.,2.5);
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetTitle("");
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetTitle("#frac{Data}{Prediction}");
  //ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetTitle("");
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetTitleOffset(0.5);
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetTitleSize(0.12);
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetTitleSize(0.2);
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetLabelSize(0.13);
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetLabelSize(0.12);
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetNdivisions(502,0);
  //MakeBlindMet30(ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs);
  //MakeBlindMet30(h_SystErrWmu);
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetRangeUser(0,179);
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetMarkerSize(2);
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetLineWidth(3);
  //ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetErrorX(0);
  for(int i=0;i<ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetN();++i){
    ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetPointEXlow(i,0);
    ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetPointEXhigh(i,0);
  }
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Draw("APE");
  h_SystErrWmu->Draw("E2SAMES");
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Draw("PESAMES");
  l1.DrawLine(0,1,metPlotXmaxXtraWide,1);
  //PrelimText->Draw();
  p3->cd();p3->Clear();
  met2->Draw();
  c2->Print("Plots/Higgs/Exclusive_WmuDataMetWithRatio.png");
  c2->Print("Plots/Higgs/Exclusive_WmuDataMetWithRatio.pdf");
 
  p1->cd();p1->SetLogy(0);
  ggMetWmu_sigInvMassRebin_tgraph_clone->GetYaxis()->SetRangeUser(0,5);
  ggMetWmu_sigInvMassRebin_tgraph_clone->Draw("APE");
  WmuMetStackRebin->Draw("histoSAMES");
  ggMetWmu_InvMassSBcomb_plus_SMhiggs->Draw("E2SAMES");
  ggMetWmu_sigInvMassRebin_tgraph_clone->Draw("PESAMES");
  h_SMS_HH_all_1Mu_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_ZH_gg_1Mu_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_WH_gg_1Mu_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_HH_2Z2g_gg_1Mu_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_HH_2W2g_gg_1Mu_Excluded_Rebin->Draw("histoSAMES");
  ggMetWmu_sigInvMassRebin_tgraph_clone->Draw("PESAMES");
  p1->RedrawAxis();
  PrelimText->Draw();
  ggTextWmu->Draw();
  WmuMetLeg->Draw();
  WmuMetLeg2->Draw();
  textSMS->Draw(); 
  TLine metVertErr3(104.13,4.3,104.13,4.7);metVertErr3.SetLineWidth(3);
  metVertErr3.Draw();
  p2->cd();
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Draw("APE");
  h_SystErrWmu->Draw("E2SAMES");
  ggMetWmu_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Draw("PESAMES");
  l1.DrawLine(0,1,metPlotXmaxXtraWide,1);
  //PrelimText->Draw();
  p3->cd();p3->Clear();
  met2->Draw();
  c2->Print("Plots/Higgs/Exclusive_WmuDataMetWithRatio_linear.png");
  c2->Print("Plots/Higgs/Exclusive_WmuDataMetWithRatio_linear.pdf");
  p1->SetLogy(1);
  //W->mu MT plot
  c1->cd();
  TH2F* WmuMTData = (TH2F*)fin->Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
  TH2F* WmuMTggHgg = (TH2F*)f_ggHgg->Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
  TH2F* WmuMTWZHgg = (TH2F*)f_WZHgg->Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
  TH2F* WmuMTTTHgg = (TH2F*)f_TTHgg->Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
  TH2F* WmuMTVBFHgg = (TH2F*)f_VBFHgg->Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
  WmuMTggHgg->Scale((ggHnormSF_mu*PhoEffScale2*L_int*2.29e-03*19.22)/99989.);//125GeV=/96290);  
  WmuMTVBFHgg->Scale((VBFHnormSF_mu*PhoEffScale2*L_int*2.29e-03*1.544)/95677.);//125GeV=/99885); 
  WmuMTTTHgg->Scale((TTHnormSF_mu*PhoEffScale2*L_int*2.29e-03*.1271)/100048.);//125GeV=/100224);
  WmuMTWZHgg->Scale((WZHnormSF_mu*PhoEffScale2*L_int*2.29e-03*(.6782/**(.3257+.014)*/+.3843/**.2*/))/100320);//125 and 126 GeV have same # events
  WmuMTggHgg->SetLineColor(kGreen);WmuMTggHgg->SetMarkerColor(kGreen);WmuMTggHgg->SetFillColor(kGreen);
  WmuMTWZHgg->SetLineColor(kCyan);WmuMTWZHgg->SetMarkerColor(kCyan);WmuMTWZHgg->SetFillColor(kCyan);
  WmuMTTTHgg->SetLineColor(31);WmuMTTTHgg->SetMarkerColor(31);WmuMTTTHgg->SetFillColor(31);
  WmuMTVBFHgg->SetLineColor(kRed+3);WmuMTVBFHgg->SetMarkerColor(kRed+3);WmuMTVBFHgg->SetFillColor(kRed+3);
  //WmuMTggHgg->SetFillStyle(0);WmuMTWZHgg->SetFillStyle(0);WmuMTTTHgg->SetFillStyle(0);WmuMTVBFHgg->SetFillStyle(0);
  WmuMTggHgg->SetLineWidth(2);WmuMTWZHgg->SetLineWidth(2);WmuMTVBFHgg->SetLineWidth(2);WmuMTTTHgg->SetLineWidth(2);
  TH1D* WmuMTDataProjY = (TH1D*)WmuMTData->ProjectionY("WmuMTDataProjY",0,-1,"eo");
  TH1D* WmuMTggHggProjY = (TH1D*)WmuMTggHgg->ProjectionY("WmuMTggHggProjY",0,-1,"eo");
  TH1D* WmuMTWZHggProjY = (TH1D*)WmuMTWZHgg->ProjectionY("WmuMTWZHggProjY",0,-1,"eo");
  TH1D* WmuMTTTHggProjY = (TH1D*)WmuMTTTHgg->ProjectionY("WmuMTTTHggProjY",0,-1,"eo");
  TH1D* WmuMTVBFHggProjY = (TH1D*)WmuMTVBFHgg->ProjectionY("WmuMTVBFHggProjY",0,-1,"eo");
  //MakeBlindInvMass(WmuMTDataProjX);
  c1->SetLogy(1);
  WmuMTDataProjY->GetXaxis()->SetRangeUser(0,200);WmuMTDataProjY->Draw();
  c1->Print("Plots/Higgs/Exclusive_WmuDataMT.png");
  /*
  overFlowBin=WmuMTggHggProjY->FindBin(metPlotXmax+1);lastBin=WmuMTggHggProjY->FindBin(metPlotXmax-1);
  overFlow=0.;overFlowErr=0.;
  val=WmuMTggHggProjY->GetBinContent(lastBin);valErr=WmuMTggHggProjY->GetBinError(lastBin);
  overFlow=WmuMTggHggProjY->IntegralAndError(overFlowBin,-1,overFlowErr);
  WmuMTggHggProjY->SetBinContent(lastBin,val+overFlow);
  WmuMTggHggProjY->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  val=WmuMTTTHggProjY->GetBinContent(lastBin);valErr=WmuMTTTHggProjY->GetBinError(lastBin);
  overFlow=WmuMTTTHggProjY->IntegralAndError(overFlowBin,-1,overFlowErr);
  WmuMTTTHggProjY->SetBinContent(lastBin,val+overFlow);
  WmuMTTTHggProjY->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  val=WmuMTWZHggProjY->GetBinContent(lastBin);valErr=WmuMTWZHggProjY->GetBinError(lastBin);
  overFlow=WmuMTWZHggProjY->IntegralAndError(overFlowBin,-1,overFlowErr);
  WmuMTWZHggProjY->SetBinContent(lastBin,val+overFlow);
  WmuMTWZHggProjY->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  val=WmuMTVBFHggProjY->GetBinContent(lastBin);valErr=WmuMTVBFHggProjY->GetBinError(lastBin);
  overFlow=WmuMTVBFHggProjY->IntegralAndError(overFlowBin,-1,overFlowErr);
  WmuMTVBFHggProjY->SetBinContent(lastBin,val+overFlow);
  WmuMTVBFHggProjY->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  */
  THStack *WmuMTStack = new THStack("WmuMTStack","");
  WmuMTStack->Add(WmuMTTTHggProjY);WmuMTStack->Add(WmuMTVBFHggProjY);WmuMTStack->Add(WmuMTWZHggProjY);WmuMTStack->Add(WmuMTggHggProjY);
  WmuMTStack->Draw("histo");WmuMTStack->GetXaxis()->SetRangeUser(0,249.9);
  float SMHiggsMetIntWmuMT=WmuMTTTHggProjY->Integral()+WmuMTVBFHggProjY->Integral()+WmuMTWZHggProjY->Integral()+WmuMTggHggProjY->Integral();
  c1->Print("Plots/Higgs/Exclusive_WmuSMHiggsMT.png");
  TH1F* WmuMTggHggProjYRebin=(TH1F*)WmuMTggHggProjY->Rebin(nMTbins,"WmuMTggHggProjYRebin",MTbins);
  TH1F* WmuMTTTHggProjYRebin=(TH1F*)WmuMTTTHggProjY->Rebin(nMTbins,"WmuMTTTHggProjYRebin",MTbins);
  TH1F* WmuMTWZHggProjYRebin=(TH1F*)WmuMTWZHggProjY->Rebin(nMTbins,"WmuMTWZHggProjYRebin",MTbins);
  TH1F* WmuMTVBFHggProjYRebin=(TH1F*)WmuMTVBFHggProjY->Rebin(nMTbins,"WmuMTVBFHggProjYRebin",MTbins);

  AddOverflowToLastBin(WmuMTggHggProjYRebin);
  AddOverflowToLastBin(WmuMTTTHggProjYRebin);
  AddOverflowToLastBin(WmuMTWZHggProjYRebin);
  AddOverflowToLastBin(WmuMTVBFHggProjYRebin);
  DivideBy30gev(WmuMTggHggProjYRebin);
  DivideBy30gev(WmuMTTTHggProjYRebin);
  DivideBy30gev(WmuMTWZHggProjYRebin);
  DivideBy30gev(WmuMTVBFHggProjYRebin);

  THStack *WmuMTStackRebin = new THStack("WmuMTStackRebin","");
  THStack *WmuMTStackRebin2 = new THStack("WmuMTStackRebin2","");
  WmuMTTTHggProjYRebin->SetLineColor(kBlack);WmuMTVBFHggProjYRebin->SetLineColor(kBlack);WmuMTWZHggProjYRebin->SetLineColor(kBlack);WmuMTggHggProjYRebin->SetLineColor(kBlack);
  WmuMTTTHggProjYRebin->SetLineWidth(1);WmuMTVBFHggProjYRebin->SetLineWidth(1);WmuMTWZHggProjYRebin->SetLineWidth(1);WmuMTggHggProjYRebin->SetLineWidth(1);
  TH1F* SMHiggsAdditionForWmuMT = (TH1F*)WmuMTTTHggProjYRebin->Clone();SMHiggsAdditionForWmuMT->Add(WmuMTVBFHggProjYRebin);SMHiggsAdditionForWmuMT->Add(WmuMTWZHggProjYRebin);SMHiggsAdditionForWmuMT->Add(WmuMTggHggProjYRebin);
  WmuMTStackRebin2->Add(WmuMTTTHggProjYRebin);WmuMTStackRebin2->Add(WmuMTVBFHggProjYRebin);WmuMTStackRebin2->Add(WmuMTWZHggProjYRebin);WmuMTStackRebin2->Add(WmuMTggHggProjYRebin);
  WmuMTStackRebin->Add(SMHiggsAdditionForWmuMT);
  WmuMTStackRebin2->Draw("histo");WmuMTStackRebin2->GetXaxis()->SetRangeUser(0,249.9);
  WmuMTStackRebin2->SetMaximum(6);WmuMTStackRebin2->SetMinimum(3.1E-3);
  WmuMTStackRebin2->GetYaxis()->SetTitle("Events / 30 GeV");WmuMTStackRebin2->GetXaxis()->SetTitle("M_{T} [GeV]");
  legHiggs2->Draw();
  c1->Print("Plots/Higgs/Exclusive_WmuSMHiggsMTRebin.png");
  c1->Print("Plots/Higgs/Exclusive_WmuSMHiggsMTRebin.pdf");
  TH1F* SMHiggsMTWmuMT = (TH1F*)WmuMTTTHggProjYRebin->Clone();SMHiggsMTWmuMT->Add(WmuMTVBFHggProjYRebin);SMHiggsMTWmuMT->Add(WmuMTggHggProjYRebin);SMHiggsMTWmuMT->Add(WmuMTWZHggProjYRebin);
  //MT plot with ratio
  //p1->cd();
  sbLoBin1=WmuMTData->GetXaxis()->FindBin(sbLoLo);sbLoBin2=WmuMTData->GetXaxis()->FindBin(sbLoHi)-1;  
  sbHiBin1=WmuMTData->GetXaxis()->FindBin(sbHiLo);sbHiBin2=WmuMTData->GetXaxis()->FindBin(sbHiHi)-1;
  sigBin1=WmuMTData->GetXaxis()->FindBin(sigLo);sigBin2=WmuMTData->GetXaxis()->FindBin(sigHi)-1;
  TH1F* ggMetWmuMT_SBloInvMass=(TH1F*)WmuMTData->ProjectionY("ggMetWmuMT_SBloInvMass",sbLoBin1,sbLoBin2,"eo");
  TH1F* ggMetWmuMT_SBhiInvMass=(TH1F*)WmuMTData->ProjectionY("ggMetWmuMT_SBhiInvMass",sbHiBin1,sbHiBin2,"eo");
  TH1F* ggMetWmuMT_sigInvMass=(TH1F*)WmuMTData->ProjectionY("ggMetWmuMT_sigInvMass",sigBin1,sigBin2,"eo");
  /*
  overFlowBin=ggMetWmuMT_sigInvMass->FindBin(metPlotXmax+1), lastBin=ggMetWmuMT_sigInvMass->FindBin(metPlotXmax-1);
  overFlow=0.;overFlowErr=0.;
  val=ggMetWmuMT_sigInvMass->GetBinContent(lastBin);valErr=ggMetWmuMT_sigInvMass->GetBinError(lastBin);
  overFlow=ggMetWmuMT_sigInvMass->IntegralAndError(overFlowBin,-1,overFlowErr);
  ggMetWmuMT_sigInvMass->SetBinContent(lastBin,val+overFlow);
  ggMetWmuMT_sigInvMass->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  overFlow=0.;overFlowErr=0.;
  val=ggMetWmuMT_SBloInvMass->GetBinContent(lastBin);valErr=ggMetWmuMT_SBloInvMass->GetBinError(lastBin);
  overFlow=ggMetWmuMT_SBloInvMass->IntegralAndError(overFlowBin,-1,overFlowErr);
  ggMetWmuMT_SBloInvMass->SetBinContent(lastBin,val+overFlow);
  ggMetWmuMT_SBloInvMass->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  overFlow=0.;overFlowErr=0.;
  val=ggMetWmuMT_SBhiInvMass->GetBinContent(lastBin);valErr=ggMetWmuMT_SBhiInvMass->GetBinError(lastBin);
  overFlow=ggMetWmuMT_SBhiInvMass->IntegralAndError(overFlowBin,-1,overFlowErr);
  ggMetWmuMT_SBhiInvMass->SetBinContent(lastBin,val+overFlow);
  ggMetWmuMT_SBhiInvMass->SetBinError(lastBin,sqrt(valErr*valErr+overFlowErr*overFlowErr));
  */
  TH1F* ggMetWmuMT_SBloInvMassRebin = (TH1F*)ggMetWmuMT_SBloInvMass->Rebin(nMTbins,"ggMetWmuMT_SBloInvMassRebin",MTbins);
  TH1F* ggMetWmuMT_SBhiInvMassRebin = (TH1F*)ggMetWmuMT_SBhiInvMass->Rebin(nMTbins,"ggMetWmuMT_SBhiInvMassRebin",MTbins);
  TH1F* ggMetWmuMT_sigInvMassRebin = (TH1F*)ggMetWmuMT_sigInvMass->Rebin(nMTbins,"ggMetWmuMT_sigInvMassRebin",MTbins);
 
  AddOverflowToLastBin(ggMetWmuMT_SBloInvMassRebin);
  AddOverflowToLastBin(ggMetWmuMT_SBhiInvMassRebin);
  AddOverflowToLastBin(ggMetWmuMT_sigInvMassRebin);

  float higgslowscaleWmuMT = ggMetWmuMT_SBloInvMassRebin->Integral() ? WmuPowYieldSig/*(HiggSigNoNormWmu)*//ggMetWmuMT_SBloInvMassRebin->Integral() : 0;
  float higgshighscaleWmuMT = ggMetWmuMT_SBhiInvMassRebin->Integral() ? WmuPowYieldSig/*(HiggSigNoNormWmu)*//ggMetWmuMT_SBhiInvMassRebin->Integral() : 0;
  float WmuMTFitSyst=(WmuPowYieldSigErr*WmuPowYieldSigErr)/(WmuPowYieldSig*WmuPowYieldSig);
  double WmuMTbgStatErrLO=0.;double WmuMTbgStatLO=ggMetWmuMT_SBloInvMassRebin->IntegralAndError(0,-1,WmuMTbgStatErrLO);
  float WmuMTbgStatSystLO=(WmuMTbgStatErrLO*WmuMTbgStatErrLO)/(WmuMTbgStatLO*WmuMTbgStatLO);
  double WmuMTbgStatErrHI=0.;double WmuMTbgStatHI=ggMetWmuMT_SBhiInvMassRebin->IntegralAndError(0,-1,WmuMTbgStatErrHI);
  float WmuMTbgStatSystHI=(WmuMTbgStatErrHI*WmuMTbgStatErrHI)/(WmuMTbgStatHI*WmuMTbgStatHI); 

  TH1F* ggMetWmuMT_SBloInvMassClone = (TH1F*)ggMetWmuMT_SBloInvMass->Clone("ggMetWmuMT_SBloInvMassClone");
  TH1F* ggMetWmuMT_SBhiInvMassClone = (TH1F*)ggMetWmuMT_SBhiInvMass->Clone("ggMetWmuMT_SBhiInvMassClone");
  ggMetWmuMT_SBloInvMassClone->Scale(higgslowscale);ggMetWmuMT_SBhiInvMassClone->Scale(higgshighscale);
  TH1F* ggMetWmuMT_SBloInvMassCloneRebin = (TH1F*)ggMetWmuMT_SBloInvMassClone->Rebin(nMTbins,"ggMetWmuMT_SBloInvMassCloneRebin",MTbins);
  TH1F* ggMetWmuMT_SBhiInvMassCloneRebin = (TH1F*)ggMetWmuMT_SBhiInvMassClone->Rebin(nMTbins,"ggMetWmuMT_SBhiInvMassCloneRebin",MTbins);
  TH1F* ggMetWmuMT_InvMassSBcombNoRebin=(TH1F*)ggMetWmuMT_SBloInvMass->Clone();//ggMetWmuMT_InvMassSBcombNoRebin->Add(ggMetWmuMT_SBhiInvMass);ggMetWmuMT_InvMassSBcombNoRebin->Scale(0.5);
  TH1F *StatErrsWmuMT       =(TH1F*)ggMetWmuMT_SBloInvMass->Clone();
  TH1F *FitStatSystErrsWmuMT =(TH1F*)ggMetWmuMT_SBloInvMass->Clone();
  TH1F *FitShapeSystErrsWmuMT  =(TH1F*)ggMetWmuMT_SBloInvMass->Clone();
  TH1F *HalfDiffErrsWmuMT  =(TH1F*)ggMetWmuMT_SBloInvMass->Clone();
  for(int i=1;i<=ggMetWmuMT_InvMassSBcombNoRebin->GetNbinsX();i++){
    Double_t val = ((ggMetWmuMT_SBloInvMass->GetBinContent(i)*higgslowscaleWmuMT)+(ggMetWmuMT_SBhiInvMass->GetBinContent(i)*higgshighscaleWmuMT))/2.;
    //double Lerr=0.,Uerr=0.;
    //double Lval = ggMetWmuMT_SBloInvMass->IntegralAndError(0,-1,Lerr);
    //double Uval = ggMetWmuMT_SBhiInvMass->IntegralAndError(0,-1,Uerr);
    Double_t LvalScaleWmuMT = ggMetWmuMT_SBloInvMass->GetBinContent(i)*higgslowscaleWmuMT,UvalScaleWmuMT = ggMetWmuMT_SBhiInvMass->GetBinContent(i)*higgshighscaleWmuMT;
    Double_t WmuMTStatErr(0.),WmuMTFitStatSystErr(0.),WmuMTFitShapeSystErr(0.),WmuMTHalfDiffErr(0.);
    double err = TotErr(val,ggMetWmuMT_SBloInvMass->GetBinContent(i),ggMetWmuMT_SBloInvMass->GetBinError(i),ggMetWmuMT_SBhiInvMass->GetBinContent(i),ggMetWmuMT_SBhiInvMass->GetBinError(i),WmuMTbgStatLO,WmuMTbgStatErrLO,WmuMTbgStatHI,WmuMTbgStatErrHI,WmuPowYieldSig,WmuPowYieldSigErr,WmuExpoYieldSig,WmuMTStatErr,WmuMTFitStatSystErr,WmuMTFitShapeSystErr,LvalScaleWmuMT,UvalScaleWmuMT,WmuMTHalfDiffErr);
    ggMetWmuMT_InvMassSBcombNoRebin->SetBinContent(i,val);ggMetWmuMT_InvMassSBcombNoRebin->SetBinError(i,err);
    StatErrsWmuMT->SetBinError(i,WmuMTStatErr);FitStatSystErrsWmuMT->SetBinError(i,WmuMTFitStatSystErr);FitShapeSystErrsWmuMT->SetBinError(i,WmuMTFitShapeSystErr);HalfDiffErrsWmuMT->SetBinError(i,WmuMTHalfDiffErr);
  }
  TH1F* ggMetWmuMT_InvMassSBcomb = (TH1F*)ggMetWmuMT_InvMassSBcombNoRebin->Rebin(nMTbins,"ggMetWmuMT_InvMassSBcomb",MTbins);
  for(int i=0;i<ggMetWmuMT_sigInvMassRebin->FindBin(31);i++){
    cout<<"WmuMT "<<ggMetWmuMT_sigInvMassRebin->GetBinLowEdge(i)<<"<"<<ggMetWmuMT_sigInvMassRebin->GetBinLowEdge(i+1)<<"  events: "<<ggMetWmuMT_sigInvMassRebin->GetBinContent(i)<<" +- "<<ggMetWmuMT_sigInvMassRebin->GetBinError(i)<<endl;
  }
  AddOverflowToLastBin(ggMetWmuMT_InvMassSBcomb);
  fout_Wmu.cd();
  TH1F* ggMetWmuMT_sigInvMassRebin_clone = (TH1F*)ggMetWmuMT_sigInvMassRebin->Clone();
  //MakeBlindMet30(ggMetWmuMT_sigInvMassRebin_clone);
  ggMetWmuMT_sigInvMassRebin_clone->SetTitle("");ggMetWmuMT_sigInvMassRebin_clone->GetXaxis()->SetTitle("M_{T} [GeV]");
  ggMetWmuMT_sigInvMassRebin_clone->Write("hMT_ggMu_tag");
  fout.cd();
  ggMetWmuMT_sigInvMassRebin_clone->Write("hMT_ggMu_tag");
  DivideBy30gev(ggMetWmuMT_InvMassSBcomb);
  DivideBy30gev(ggMetWmuMT_sigInvMassRebin);
  /*
  ggMetWmuMT_SBloInvMassRebin->Scale(higgslowscaleWmuMT);
  ggMetWmuMT_SBhiInvMassRebin->Scale(higgshighscaleWmuMT);

  AddOverflowToLastBin(ggMetWmuMT_SBloInvMassRebin);
  AddOverflowToLastBin(ggMetWmuMT_SBhiInvMassRebin);
  AddOverflowToLastBin(ggMetWmuMT_sigInvMassRebin);
  DivideByBinWidth(ggMetWmuMT_SBloInvMassRebin);
  DivideByBinWidth(ggMetWmuMT_SBhiInvMassRebin);
  DivideByBinWidth(ggMetWmuMT_sigInvMassRebin);
  
  TH1F* ggMetWmuMT_InvMassSBcomb=(TH1F*)ggMetWmuMT_SBloInvMassRebin->Clone();ggMetWmuMT_InvMassSBcomb->Add(ggMetWmuMT_SBhiInvMassRebin);ggMetWmuMT_InvMassSBcomb->Scale(0.5);
  */


  c1->cd();c1->SetLogy(0);
  ggMetWmuMT_SBloInvMassCloneRebin->SetFillColor(kRed);ggMetWmuMT_SBloInvMassCloneRebin->SetFillStyle(3004);ggMetWmuMT_SBloInvMassCloneRebin->SetMarkerSize(0);
  ggMetWmuMT_SBloInvMassCloneRebin->GetXaxis()->SetRangeUser(0,250);
  ggMetWmuMT_SBhiInvMassCloneRebin->SetFillColor(kBlue);ggMetWmuMT_SBhiInvMassCloneRebin->SetFillStyle(3004);ggMetWmuMT_SBhiInvMassCloneRebin->SetMarkerSize(0);
  ggMetWmuMT_SBloInvMassCloneRebin->SetLineColor(kRed);ggMetWmuMT_SBloInvMassCloneRebin->SetMarkerColor(kRed);
  ggMetWmuMT_SBhiInvMassCloneRebin->SetLineColor(kBlue);ggMetWmuMT_SBhiInvMassCloneRebin->SetMarkerColor(kBlue);
  ggMetWmuMT_SBhiInvMassCloneRebin->SetLineWidth(2);ggMetWmuMT_SBloInvMassCloneRebin->SetLineWidth(2);
  ggMetWmuMT_SBhiInvMassCloneRebin->SetTitle("");ggMetWmuMT_SBhiInvMassCloneRebin->GetXaxis()->SetTitle("M_{T}");ggMetWmuMT_SBhiInvMassCloneRebin->GetYaxis()->SetTitle("Events / 30 GeV");
  ggMetWmuMT_SBhiInvMassCloneRebin->Draw("E2");ggMetWmuMT_SBhiInvMassCloneRebin->Draw("pesames");
  ggMetWmuMT_SBloInvMassCloneRebin->Draw("E2sames");ggMetWmuMT_SBloInvMassCloneRebin->Draw("pesames");
  ggMetWmuMT_sigInvMassRebin_clone->Draw("PESAMES");
  lowhighleg->Draw();
  c1->Print("Plots/Higgs/Exclusive_WmuMT_SbLowAndHigh.png");
  c1->Print("Plots/Higgs/Exclusive_WmuMT_SbLowAndHigh.pdf");
  c1->SetLogy(1);
  p1->cd();

  ggMetWmuMT_sigInvMassRebin->SetMarkerColor(kBlack);ggMetWmuMT_sigInvMassRebin->SetLineColor(kBlack);
  ggMetWmuMT_SBloInvMassRebin->SetMarkerColor(kBlue);ggMetWmuMT_SBloInvMassRebin->SetLineColor(kBlue);ggMetWmuMT_SBloInvMassRebin->SetFillColor(kBlue);ggMetWmuMT_SBloInvMassRebin->SetFillStyle(3004);
  ggMetWmuMT_SBhiInvMassRebin->SetMarkerColor(kRed);ggMetWmuMT_SBhiInvMassRebin->SetLineColor(kRed);ggMetWmuMT_SBhiInvMassRebin->SetFillColor(kRed);ggMetWmuMT_SBhiInvMassRebin->SetFillStyle(3004);
  ggMetWmuMT_InvMassSBcomb->SetMarkerColor(42);ggMetWmuMT_InvMassSBcomb->SetLineColor(42);ggMetWmuMT_InvMassSBcomb->SetFillColor(42);//ggMetWmuMT_InvMassSBcomb->SetFillStyle(3004);
  ggMetWmuMT_SBloInvMassRebin->SetMarkerSize(0.75);ggMetWmuMT_SBhiInvMassRebin->SetMarkerSize(0.75);ggMetWmuMT_InvMassSBcomb->SetMarkerSize(0.75);
  WmuMTStackRebin->Add(ggMetWmuMT_InvMassSBcomb);
  TH1F* ggMetWmuMT_InvMassSBcomb_plus_SMhiggs = (TH1F*)ggMetWmuMT_InvMassSBcomb->Clone();ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->Add(SMHiggsMTWmuMT);
  TH1F* ggMetWmuMT_InvMassSBcomb_plus_SMhiggs_clone = (TH1F*)ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->Clone();
  TH1F* ggMetWmuMT_InvMassSBcomb_LimitClone = (TH1F*)ggMetWmuMT_InvMassSBcomb->Clone();
  TH1F* SMHiggsMTWmuMT_LimitCLone = (TH1F*)SMHiggsMTWmuMT->Clone();
  MultBy30gev(ggMetWmuMT_InvMassSBcomb_plus_SMhiggs_clone);
  MultBy30gev(ggMetWmuMT_InvMassSBcomb_LimitClone);
  MultBy30gev(SMHiggsMTWmuMT_LimitCLone);

  fout_Wmu.cd();
  ggMetWmuMT_InvMassSBcomb_plus_SMhiggs_clone->Write("smBackground");
  fout.cd();
  ggMetWmuMT_InvMassSBcomb_LimitClone->Write("hMT_ggMu_bkg");
  fout_SMHiggs.cd();
  SMHiggsMTWmuMT_LimitCLone->Write("SMHiggsMET_ggMu");
  ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->SetMarkerSize(0);ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->SetFillColor(kBlack);ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->SetFillStyle(3004);
  //const double alpha = 1 - 0.6827;
  TGraphAsymmErrors *ggMetWmuMT_sigInvMassRebin_tgraph = new TGraphAsymmErrors(ggMetWmuMT_sigInvMassRebin);
  //ggMetWmuMT_sigInvMassRebin_tgraph->SetMarkerSize(0.5);
  ggMetWmuMT_sigInvMassRebin_tgraph->SetMarkerStyle(20);
  for (int i = 0; i < ggMetWmuMT_sigInvMassRebin_tgraph->GetN(); ++i) {
    int N = ggMetWmuMT_sigInvMassRebin_tgraph->GetY()[i];
    double L=0.,U=0.;
    if(i!=3){
      L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
      U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
      ggMetWmuMT_sigInvMassRebin_tgraph->SetPointEYlow(i, N-L);
      ggMetWmuMT_sigInvMassRebin_tgraph->SetPointEYhigh(i, U-N);
    }
    else{
      L =  (ROOT::Math::gamma_quantile(alpha/2,1,1.));
      U =  ROOT::Math::gamma_quantile_c(alpha/2,1+1,1);
      L/=3.;U/=3.;
      ggMetWmuMT_sigInvMassRebin_tgraph->SetPointEYlow(i, 1./3.-L);
      ggMetWmuMT_sigInvMassRebin_tgraph->SetPointEYhigh(i, U-1./3.);
    }
    //      double U =  (N==0) ?  ( ROOT::Math::gamma_quantile_c(alpha,N+1,1) ) :
    //         ( ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) );
    //ggMetWmuMT_sigInvMassRebin_tgraph->SetPointEYlow(i, N-L);
    //ggMetWmuMT_sigInvMassRebin_tgraph->SetPointEYhigh(i, U-N);
    cout<<"mu i:"<<i<<"  N: "<<N<<"  L: "<<L<<"  N-L: "<<N-L<<"  U: "<<U<<"  U-N: "<<U-N<<endl;
  }
  

  ggMetWmuMT_sigInvMassRebin_tgraph->GetXaxis()->SetRangeUser(0,179);
  ggMetWmuMT_sigInvMassRebin_tgraph->GetXaxis()->SetLabelSize(0);
  ggMetWmuMT_sigInvMassRebin_tgraph->GetYaxis()->SetRangeUser(9e-3,50);
  ggMetWmuMT_sigInvMassRebin_tgraph->GetYaxis()->SetTitle("Events / 30 GeV");
  ggMetWmuMT_sigInvMassRebin_tgraph->GetYaxis()->SetTitleSize(0.05);
  ggMetWmuMT_sigInvMassRebin_tgraph->SetTitle("");
  ggMetWmuMT_sigInvMassRebin_tgraph->SetMarkerSize(1);
  //MakeBlindMet30(ggMetWmuMT_sigInvMassRebin);
  

  TGraphAsymmErrors* ggMetWmuMT_sigInvMassRebin_tgraph_clone = (TGraphAsymmErrors*)ggMetWmuMT_sigInvMassRebin_tgraph->Clone("ggMetWmuMT_sigInvMassRebin_tgraph_clone");

  //ggMetWmuMT_sigInvMassRebin_tgraph_clone->SetPointEYhigh(2,0);
  //ggMetWmuMT_sigInvMassRebin_tgraph_clone->SetPointEYlow(3,ggMetWmuMT_sigInvMassRebin_tgraph->GetY()[3]);
  ggMetWmuMT_sigInvMassRebin_tgraph_clone->SetMarkerSize(2);
  ggMetWmuMT_sigInvMassRebin_tgraph_clone->SetLineWidth(3);
  for(int i=0;i<ggMetWmuMT_sigInvMassRebin_tgraph_clone->GetN();++i){
    ggMetWmuMT_sigInvMassRebin_tgraph_clone->SetPointEXlow(i,0);
    ggMetWmuMT_sigInvMassRebin_tgraph_clone->SetPointEXhigh(i,0);
  }
  ggMetWmuMT_sigInvMassRebin_tgraph_clone->Draw("APE");
  WmuMTStackRebin->Draw("histoSAMES");
  //MakeBlindMet30(ggMetWmuMT_InvMassSBcomb_plus_SMhiggs);
  ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->Draw("E2SAMES");
  ggMetWmuMT_sigInvMassRebin_tgraph_clone->Draw("PESAMES");
  h_SMS_WH_gg_1Mu_MT_Excluded_Rebin->SetLineStyle(7);
  h_SMS_ZH_gg_1Mu_MT_Excluded_Rebin->SetLineStyle(9);
  //h_SMS_HH_2b2g_gg_1Mu_MT_Excluded_Rebin->Draw("histoSAMES");
  TH1F* h_SMS_HH_all_1Mu_MT_Excluded_Rebin = (TH1F*)h_SMS_HH_2W2g_gg_1Mu_MT_Excluded_Rebin->Clone();h_SMS_HH_all_1Mu_MT_Excluded_Rebin->Add(h_SMS_HH_2Z2g_gg_1Mu_MT_Excluded_Rebin);h_SMS_HH_all_1Mu_MT_Excluded_Rebin->Add(h_SMS_HH_2tau2g_gg_1Mu_MT_Excluded_Rebin);
  //h_SMS_HH_2W2g_gg_1Mu_MT_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_HH_2Z2g_gg_1Mu_MT_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_HH_all_1Mu_MT_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_ZH_gg_1Mu_MT_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_WH_gg_1Mu_MT_Excluded_Rebin->Draw("histoSAMES");
  ggMetWmuMT_sigInvMassRebin_tgraph_clone->Draw("PESAMES");
  p1->RedrawAxis();
  PrelimText->Draw();
  /*ggTextWmu->Draw();
  WmuMetLeg->Draw();
  WmuMetLeg2->Draw();*/
  ggTextWmuLog->Draw();
  WmuMetLegLog->Draw();
  WmuMetLeg2Log->Draw();
  textSMSlog->Draw();
  p2->cd();
  TH1F* ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs = (TH1F*)ggMetWmuMT_sigInvMassRebin->Clone();
  //TGraphAsymmErrors* ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph = new TGraphAsymmErrors(ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs);
  //ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Divide(ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs,ggMetWmuMT_InvMassSBcomb_plus_SMhiggs,"cp");
  //ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs->Divide(ggMetWmuMT_InvMassSBcomb_plus_SMhiggs);
  TH1F* h_WmuSystErrMT = (TH1F*)ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->Clone();h_WmuSystErrMT->SetFillColor(kBlack);h_WmuSystErrMT->SetFillStyle(3004);h_WmuSystErrMT->SetMarkerSize(0);



  for(int i=1;i<ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs->GetNbinsX()+1;i++){

   
    float Value = ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs->GetBinContent(i);float StatErr = ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs->GetBinError(i);
    Value/=ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->GetBinContent(i);StatErr/=ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->GetBinContent(i);
    ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs->SetBinContent(i,Value);ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs->SetBinError(i,StatErr);

    float SystErr=h_WmuSystErrMT->GetBinError(i)/h_WmuSystErrMT->GetBinContent(i);
    h_WmuSystErrMT->SetBinContent(i,1); h_WmuSystErrMT->SetBinError(i,SystErr);
  }

  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs->SetTitleFont(42,"xy");
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs->SetLabelFont(42,"xy");

  TGraphAsymmErrors* ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph = new TGraphAsymmErrors(ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs);
  for (int i = 0; i < ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetN(); ++i) {

    float N = ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetY()[i];
    double U = ggMetWmuMT_sigInvMassRebin_tgraph->GetErrorYhigh(i);
    double L = ggMetWmuMT_sigInvMassRebin_tgraph->GetErrorYlow(i);
    double Utemp=U,Ltemp=L;
    U/=ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->GetBinContent(i+1);L/=ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->GetBinContent(i+1);
    ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetPointEYlow(i, L);
    ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetPointEYhigh(i, U);
    if(i==0)cout<<"Ratio:"<<endl;
    cout<<"i:"<<i<<"  N: "<<N<<"  Ltemp: "<<Ltemp<<"Lbin: "<<ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->GetBinContent(i+1)<<"  L: "<<L<<"  Utemp: "<<Utemp<<"  Ubin: "<<ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->GetBinContent(i+1)<<" U: "<<U<<endl;
  }



  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetLineColor(kBlack);ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetMarkerColor(kBlack);//ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetMarkerSize(0.75);
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetTitle("");
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetTitle("");
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetRangeUser(0.,6.);
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetTitle("");
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetTitle("#frac{Data}{Prediction}");
  //ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetTitle("M_{T} [GeV]");
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetTitle("");
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetTitleOffset(0.4);
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetTitleSize(0.14);
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetTitleSize(0.);
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetLabelSize(0.12);
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetLabelSize(0.13);
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetYaxis()->SetNdivisions(206,0);
  //MakeBlindMet30(ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs);
  //MakeBlindMet30(h_WmuSystErrMT);
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetXaxis()->SetRangeUser(0,179);
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetMarkerSize(2);
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetLineWidth(3);
  for(int i=0;i<ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->GetN();++i){
    ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetPointEXlow(i,0);
    ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->SetPointEXhigh(i,0);
  }
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Draw("APE");
  h_WmuSystErrMT->Draw("E2SAMES");
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Draw("PESAMES");
  l1.DrawLine(0,1,180,1);
  //PrelimText->Draw();
  p3->cd();p3->Clear();
  mtText->Draw();
  c2->Print("Plots/Higgs/Exclusive_WmuDataMTWithRatio.png");
  c2->Print("Plots/Higgs/Exclusive_WmuDataMTWithRatio.pdf");

  p1->cd();p1->SetLogy(0);
  //ggMetWmuMT_sigInvMassRebin->GetYaxis()->SetRangeUser(0,6.49);ggMetWmuMT_sigInvMassRebin->Draw("PE");
  //ggMetWmuMT_sigInvMassRebin_tgraph_clone->RemovePoint(2);
  ggMetWmuMT_sigInvMassRebin_tgraph_clone->GetYaxis()->SetRangeUser(0,6.49);
  ggMetWmuMT_sigInvMassRebin_tgraph_clone->Draw("APE");
  WmuMTStackRebin->Draw("histoSAMES");
  ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->Draw("E2SAMES");
  ggMetWmuMT_sigInvMassRebin_tgraph_clone->Draw("PESAMES");
  h_SMS_HH_all_1Mu_MT_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_ZH_gg_1Mu_MT_Excluded_Rebin->Draw("histoSAMES");
  h_SMS_WH_gg_1Mu_MT_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_HH_2W2g_gg_1Mu_MT_Excluded_Rebin->Draw("histoSAMES");
  //h_SMS_HH_2Z2g_gg_1Mu_MT_Excluded_Rebin->Draw("histoSAMES");
  ggMetWmuMT_sigInvMassRebin_tgraph_clone->Draw("PESAMES");
  p1->RedrawAxis();
  PrelimText->Draw();
  ggTextWmu->Draw();
  WmuMetLeg->Draw();
  WmuMetLeg2->Draw();
  textSMS->Draw();
  TLine metVertErr4(104.13,5.58,104.13,6.13);metVertErr4.SetLineWidth(3);
  metVertErr4.Draw();
  p2->cd();
  float wmusig3 = 1., wmubg = ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->GetBinContent(3), wmubgerr = ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->GetBinError(3);
  float wmuratioerr3 = (wmusig3/wmubg);//*sqrt(1+wmubgerr*wmubgerr/wmubg/wmubg);
  //ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs->SetBinError(3,wmuratioerr3);
  //ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs->Draw("PE0");
  //ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs->Draw("PE");
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Draw("APE");
  h_WmuSystErrMT->Draw("E2SAMES");
  ggMetWmuMT_sigInvMassRebin_Div_comb_plusSMhiggs_tgraph->Draw("PESAMES");
  l1.DrawLine(0,1,180,1);
  //PrelimText->Draw();
  p3->cd();p3->Clear();
  mtText->Draw();
  c2->Print("Plots/Higgs/Exclusive_WmuDataMTWithRatio_linear.png");
  c2->Print("Plots/Higgs/Exclusive_WmuDataMTWithRatio_linear.pdf");
  p1->SetLogy(1);




  /*
  for(int i=0;i<StatErrs->GetNbinsX()/5;i++){
    if(StatErrs->GetBinLowEdge(i*5+1)<30){
      if(i==0)cout<<"Inclusive errors:"<<endl;
      Double_t errNewStat(0.),errNewSyst(0.),errNewShape(0.),errSigStat(0.);
      Double_t sigStat = ggMet_sigInvMass->IntegralAndError((i)*5+1,(i+1)*5,errSigStat);
      Double_t stat=StatErrs->IntegralAndError((i)*5+1,(i+1)*5,errNewStat);
      stat=FitStatSystErrs->IntegralAndError((i)*5+1,(i+1)*5,errNewSyst);
      stat=FitShapeSystErrs->IntegralAndError((i)*5+1,(i+1)*5,errNewShape);
      stat = ggMet_InvMassSBcombNoRebin->Integral((i)*5+1,(i+1)*5);
      cout<<i*5<<"<MET<"<<(i+1)*5<<" Candidate   Events: "<<sigStat<<" +- "<<errSigStat<<endl;
      cout<<i*5<<"<MET<"<<(i+1)*5<<" BG estimate Events: "<<stat<<"  StatErr: "<<errNewStat<<"  FitStatSystErr: "<<errNewSyst<<"  FitShapeSyst: "<<errNewShape<<endl;
    }
  }
 
  for(int i=0;i<StatErrsWe->GetNbinsX()/5;i++){
    if(StatErrsWe->GetBinLowEdge(i*5+1)<30){
      if(i==0)cout<<"W->e+nu errors:"<<endl;
      Double_t errNewStat(0.),errNewSyst(0.),errNewShape(0.),errSigStat(0.);
      Double_t sigStat = ggMetWe_sigInvMass->IntegralAndError((i)*5+1,(i+1)*5,errSigStat);
      Double_t stat=StatErrsWe->IntegralAndError((i)*5+1,(i+1)*5,errNewStat);
      stat=FitStatSystErrsWe->IntegralAndError((i)*5+1,(i+1)*5,errNewSyst);
      stat=FitShapeSystErrsWe->IntegralAndError((i)*5+1,(i+1)*5,errNewShape);
      stat = ggMetWe_InvMassSBcombNoRebin->Integral((i)*5+1,(i+1)*5);
      cout<<i*5<<"<MET<"<<(i+1)*5<<" Candidate   Events: "<<sigStat<<" +- "<<errSigStat<<endl;
      cout<<i*5<<"<MET<"<<(i+1)*5<<" BG estimate Events: "<<stat<<"  StatErr: "<<errNewStat<<"  FitStatSystErr: "<<errNewSyst<<"  FitShapeSyst: "<<errNewShape<<endl;
    }
  }
 
  for(int i=0;i<StatErrsWmu->GetNbinsX()/5;i++){
    if(StatErrsWmu->GetBinLowEdge(i*5+1)<30){
      if(i==0)cout<<"W->mu+nu errors:"<<endl;
      Double_t errNewStat(0.),errNewSyst(0.),errNewShape(0.),errSigStat(0.);
      Double_t sigStat = ggMetWmu_sigInvMass->IntegralAndError((i)*5+1,(i+1)*5,errSigStat);
      Double_t stat=StatErrsWmu->IntegralAndError((i)*5+1,(i+1)*5,errNewStat);
      stat=FitStatSystErrsWmu->IntegralAndError((i)*5+1,(i+1)*5,errNewSyst);
      stat=FitShapeSystErrsWmu->IntegralAndError((i)*5+1,(i+1)*5,errNewShape);
      stat = ggMetWmu_InvMassSBcombNoRebin->Integral((i)*5+1,(i+1)*5);
      cout<<i*5<<"<MET<"<<(i+1)*5<<" Candidate   Events: "<<sigStat<<" +- "<<errSigStat<<endl;
      cout<<i*5<<"<MET<"<<(i+1)*5<<" BG estimate Events: "<<stat<<"  StatErr: "<<errNewStat<<"  FitStatSystErr: "<<errNewSyst<<"  FitShapeSyst: "<<errNewShape<<endl;
    }
  }
 
  for(int i=0;i<StatErrsleftovers->GetNbinsX()/5;i++){
    if(StatErrsleftovers->GetBinLowEdge(i*5+1)<30){
      if(i==0)cout<<"leftovers errors:"<<endl;
      Double_t errNewStat(0.),errNewSyst(0.),errNewShape(0.),errSigStat(0.);
      Double_t sigStat = ggMetleftovers_sigInvMass->IntegralAndError((i)*5+1,(i+1)*5,errSigStat);
      Double_t stat=StatErrsleftovers->IntegralAndError((i)*5+1,(i+1)*5,errNewStat);
      stat=FitStatSystErrsleftovers->IntegralAndError((i)*5+1,(i+1)*5,errNewSyst);
      stat=FitShapeSystErrsleftovers->IntegralAndError((i)*5+1,(i+1)*5,errNewShape);
      stat = ggMetleftovers_InvMassSBcombNoRebin->Integral((i)*5+1,(i+1)*5);
      cout<<i*5<<"<MET<"<<(i+1)*5<<" Candidate   Events: "<<sigStat<<" +- "<<errSigStat<<endl;
      cout<<i*5<<"<MET<"<<(i+1)*5<<" BG estimate Events: "<<stat<<"  StatErr: "<<errNewStat<<"  FitStatSystErr: "<<errNewSyst<<"  FitShapeSyst: "<<errNewShape<<endl;
    }
  }
  
  */
  

  for(int i=1;i<=ggMet_InvMassSBcomb_plus_SMhiggs->GetNbinsX();i++){
      if(i==1)cout<<"Inclusive errors:"<<endl;
      Double_t errNewStat(0.),errNewSyst(0.),errNewShape(0.),errSigStat(0.),errHalfDiff(0.);
      Double_t sigStat = ggMet_sigInvMassRebin->IntegralAndError(i,i,errSigStat);
      sigStat*=ggMet_sigInvMassRebin->GetBinWidth(i);errSigStat*=ggMet_sigInvMassRebin->GetBinWidth(i);
      int errBinLow = StatErrs->FindBin(ggMet_sigInvMassRebin->GetBinLowEdge(i)),errBinHigh = StatErrs->FindBin(ggMet_sigInvMassRebin->GetBinLowEdge(i+1)-.1);
      Double_t stat=StatErrs->IntegralAndError(errBinLow,errBinHigh,errNewStat);
      stat=FitStatSystErrs->IntegralAndError(errBinLow,errBinHigh,errNewSyst);
      stat=FitShapeSystErrs->IntegralAndError(errBinLow,errBinHigh,errNewShape);
      stat=HalfDiffErrs->IntegralAndError(errBinLow,errBinHigh,errHalfDiff);
      Double_t fullErr(0.);
      stat = ggMet_InvMassSBcomb_plus_SMhiggs->IntegralAndError(i,i,fullErr);
      stat*=ggMet_InvMassSBcomb_plus_SMhiggs->GetBinWidth(i);fullErr*=ggMet_InvMassSBcomb_plus_SMhiggs->GetBinWidth(i);
      /*if(ggMet_InvMassSBcomb_plus_SMhiggs->GetBinLowEdge(i)<50)*/cout<<ggMet_sigInvMassRebin->GetBinLowEdge(i)<<"<MET<"<<ggMet_sigInvMassRebin->GetBinLowEdge(i+1)<<" Candidate   Events: "<<sigStat<<" +- "<<errSigStat<<endl;
      cout<<ggMet_sigInvMassRebin->GetBinLowEdge(i)<<"<MET<"<<ggMet_sigInvMassRebin->GetBinLowEdge(i+1)<<" BG estimate Events: "<<stat<<" +- "<<fullErr<<"  StatErr: "<<errNewStat<<"  FitStatSystErr: "<<errNewSyst<</*"  FitShapeSyst: "<<errNewShape<<*/"  HalfDiffErr: "<<errHalfDiff<<endl;
      cout<<ggMet_sigInvMassRebin->GetBinLowEdge(i)<<"<MET<"<<ggMet_sigInvMassRebin->GetBinLowEdge(i+1)<<" HH_2W2g 130_1 signal Events: "<<h_SMS_HH_2W2g_gg_1Ele_Excluded_Rebin->GetBinContent(i)*h_SMS_HH_2W2g_gg_1Ele_Excluded_Rebin->GetBinWidth(i)<<endl;
      cout<<ggMet_sigInvMassRebin->GetBinLowEdge(i)<<"<MET<"<<ggMet_sigInvMassRebin->GetBinLowEdge(i+1)<<" ZH 130_1 signal Events: "<<h_SMS_ZH_gg_1Ele_Excluded_Rebin->GetBinContent(i)*h_SMS_ZH_gg_1Ele_Excluded_Rebin->GetBinWidth(i)<<endl;
  }


  for(int i=1;i<=ggMetWe_InvMassSBcomb_plus_SMhiggs->GetNbinsX();i++){
      if(i==1)cout<<"W->e+nu errors:"<<endl;
      Double_t errNewStat(0.),errNewSyst(0.),errNewShape(0.),errSigStat(0.),errHalfDiff(0.);
      Double_t sigStat = ggMetWe_sigInvMassRebin->IntegralAndError(i,i,errSigStat);
      sigStat*=ggMetWe_sigInvMassRebin->GetBinWidth(i);errSigStat*=ggMetWe_sigInvMassRebin->GetBinWidth(i);
      int errBinLow = StatErrsWe->FindBin(ggMetWe_sigInvMassRebin->GetBinLowEdge(i)),errBinHigh = StatErrsWe->FindBin(ggMetWe_sigInvMassRebin->GetBinLowEdge(i+1)-.1);
      Double_t stat=StatErrsWe->IntegralAndError(errBinLow,errBinHigh,errNewStat);
      stat=FitStatSystErrsWe->IntegralAndError(errBinLow,errBinHigh,errNewSyst);
      stat=FitShapeSystErrsWe->IntegralAndError(errBinLow,errBinHigh,errNewShape);
      stat=HalfDiffErrsWe->IntegralAndError(errBinLow,errBinHigh,errHalfDiff);
      Double_t fullErr(0.);
      stat = ggMetWe_InvMassSBcomb_plus_SMhiggs->IntegralAndError(i,i,fullErr);
      stat*=ggMetWe_InvMassSBcomb_plus_SMhiggs->GetBinWidth(i);fullErr*=ggMetWe_InvMassSBcomb_plus_SMhiggs->GetBinWidth(i);
      /*if(ggMetWe_InvMassSBcomb_plus_SMhiggs->GetBinLowEdge(i)<30)*/cout<<ggMetWe_sigInvMassRebin->GetBinLowEdge(i)<<"<MET<"<<ggMetWe_sigInvMassRebin->GetBinLowEdge(i+1)<<" Candidate   Events: "<<sigStat<<" +- "<<errSigStat<<endl;
      cout<<ggMetWe_sigInvMassRebin->GetBinLowEdge(i)<<"<MET<"<<ggMetWe_sigInvMassRebin->GetBinLowEdge(i+1)<<" BG estimate Events: "<<stat<<" +- "<<fullErr<<"  StatErr: "<<errNewStat<<"  FitStatSystErr: "<<errNewSyst<</*"  FitShapeSyst: "<<errNewShape<<*/"  HalfDiffErr: "<<errHalfDiff<<endl;
      cout<<ggMetWe_sigInvMassRebin->GetBinLowEdge(i)<<"<MET<"<<ggMetWe_sigInvMassRebin->GetBinLowEdge(i+1)<<" HH_2W2g 130_1 signal Events: "<<h_SMS_HH_2W2g_gg_1Ele_Excluded_Rebin->GetBinContent(i)*h_SMS_HH_2W2g_gg_1Ele_Excluded_Rebin->GetBinWidth(i)<<" +- "<<h_SMS_HH_2W2g_gg_1Ele_Excluded_Rebin->GetBinError(i)*h_SMS_HH_2W2g_gg_1Ele_Excluded_Rebin->GetBinWidth(i)<<endl;
      cout<<ggMetWe_sigInvMassRebin->GetBinLowEdge(i)<<"<MET<"<<ggMetWe_sigInvMassRebin->GetBinLowEdge(i+1)<<" ZH 130_1 signal Events: "<<h_SMS_ZH_gg_1Ele_Excluded_Rebin->GetBinContent(i)*h_SMS_ZH_gg_1Ele_Excluded_Rebin->GetBinWidth(i)<<endl;
  }

  for(int i=1;i<=ggMetWmu_InvMassSBcomb_plus_SMhiggs->GetNbinsX();i++){
      if(i==1)cout<<"W->mu+nu errors:"<<endl;
      Double_t errNewStat(0.),errNewSyst(0.),errNewShape(0.),errSigStat(0.),errHalfDiff(0.);
      Double_t sigStat = ggMetWmu_sigInvMassRebin->IntegralAndError(i,i,errSigStat);
      sigStat*=ggMetWmu_sigInvMassRebin->GetBinWidth(i);errSigStat*=ggMetWmu_sigInvMassRebin->GetBinWidth(i);
      int errBinLow = StatErrsWmu->FindBin(ggMetWmu_sigInvMassRebin->GetBinLowEdge(i)),errBinHigh = StatErrsWmu->FindBin(ggMetWmu_sigInvMassRebin->GetBinLowEdge(i+1)-.1);
      Double_t stat=StatErrsWmu->IntegralAndError(errBinLow,errBinHigh,errNewStat);
      stat=FitStatSystErrsWmu->IntegralAndError(errBinLow,errBinHigh,errNewSyst);
      stat=FitShapeSystErrsWmu->IntegralAndError(errBinLow,errBinHigh,errNewShape);
      stat=HalfDiffErrsWmu->IntegralAndError(errBinLow,errBinHigh,errHalfDiff);
      Double_t fullErr(0.);
      stat = ggMetWmu_InvMassSBcomb_plus_SMhiggs->IntegralAndError(i,i,fullErr);
      stat*=ggMetWmu_InvMassSBcomb_plus_SMhiggs->GetBinWidth(i);fullErr*=ggMetWmu_InvMassSBcomb_plus_SMhiggs->GetBinWidth(i);
      /*if(ggMetWmu_InvMassSBcomb_plus_SMhiggs->GetBinLowEdge(i)<30)*/cout<<ggMetWmu_sigInvMassRebin->GetBinLowEdge(i)<<"<MET<"<<ggMetWmu_sigInvMassRebin->GetBinLowEdge(i+1)<<" Candidate   Events: "<<sigStat<<" +- "<<errSigStat<<endl;
      cout<<ggMetWmu_sigInvMassRebin->GetBinLowEdge(i)<<"<MET<"<<ggMetWmu_sigInvMassRebin->GetBinLowEdge(i+1)<<" BG estimate Events: "<<stat<<" +- "<<fullErr<<"  StatErr: "<<errNewStat<<"  FitStatSystErr: "<<errNewSyst<</*"  FitShapeSyst: "<<errNewShape<<*/"  HalfDiffErr: "<<errHalfDiff<<endl;
      cout<<ggMetWmu_sigInvMassRebin->GetBinLowEdge(i)<<"<MET<"<<ggMetWmu_sigInvMassRebin->GetBinLowEdge(i+1)<<" HH_2W2g 130_1 signal Events: "<<h_SMS_HH_2W2g_gg_1Mu_Excluded_Rebin->GetBinContent(i)*h_SMS_HH_2W2g_gg_1Mu_Excluded_Rebin->GetBinWidth(i)<<" +- "<<h_SMS_HH_2W2g_gg_1Mu_Excluded_Rebin->GetBinError(i)*h_SMS_HH_2W2g_gg_1Mu_Excluded_Rebin->GetBinWidth(i)<<endl;
      cout<<ggMetWmu_sigInvMassRebin->GetBinLowEdge(i)<<"<MET<"<<ggMetWmu_sigInvMassRebin->GetBinLowEdge(i+1)<<" ZH 130_1 signal Events: "<<h_SMS_ZH_gg_1Mu_Excluded_Rebin->GetBinContent(i)*h_SMS_ZH_gg_1Mu_Excluded_Rebin->GetBinWidth(i)<<endl;
  }


  //MT
  for(int i=1;i<=ggMetWeMT_InvMassSBcomb_plus_SMhiggs->GetNbinsX();i++){
      if(i==1)cout<<"W->e+nu MT errors:"<<endl;
      Double_t errNewStat(0.),errNewSyst(0.),errNewShape(0.),errSigStat(0.),errHalfDiff(0.);
      Double_t sigStat = ggMetWeMT_sigInvMassRebin->IntegralAndError(i,i,errSigStat);
      sigStat*=ggMetWeMT_sigInvMassRebin->GetBinWidth(i);errSigStat*=ggMetWeMT_sigInvMassRebin->GetBinWidth(i);
      int errBinLow = StatErrsWeMT->FindBin(ggMetWeMT_sigInvMassRebin->GetBinLowEdge(i)),errBinHigh = StatErrsWeMT->FindBin(ggMetWeMT_sigInvMassRebin->GetBinLowEdge(i+1)-.1);
      Double_t stat=StatErrsWeMT->IntegralAndError(errBinLow,errBinHigh,errNewStat);
      stat=FitStatSystErrsWeMT->IntegralAndError(errBinLow,errBinHigh,errNewSyst);
      stat=FitShapeSystErrsWeMT->IntegralAndError(errBinLow,errBinHigh,errNewShape);
      stat=HalfDiffErrsWeMT->IntegralAndError(errBinLow,errBinHigh,errHalfDiff);
      Double_t fullErr(0.);
      stat = ggMetWeMT_InvMassSBcomb_plus_SMhiggs->IntegralAndError(i,i,fullErr);
      stat*=ggMetWeMT_InvMassSBcomb_plus_SMhiggs->GetBinWidth(i);fullErr*=ggMetWeMT_InvMassSBcomb_plus_SMhiggs->GetBinWidth(i);
      /*if(ggMetWeMT_InvMassSBcomb_plus_SMhiggs->GetBinLowEdge(i)<30)*/cout<<ggMetWeMT_sigInvMassRebin->GetBinLowEdge(i)<<"<MT<"<<ggMetWeMT_sigInvMassRebin->GetBinLowEdge(i+1)<<" Candidate   Events: "<<sigStat<<" +- "<<errSigStat<<endl;
      cout<<ggMetWeMT_sigInvMassRebin->GetBinLowEdge(i)<<"<MT<"<<ggMetWeMT_sigInvMassRebin->GetBinLowEdge(i+1)<<" BG estimate Events: "<<stat<<" +- "<<fullErr<<"  StatErr: "<<errNewStat<<"  FitStatSystErr: "<<errNewSyst<</*"  FitShapeSyst: "<<errNewShape<<*/"  HalfDiffErr: "<<errHalfDiff<<endl;
      cout<<ggMetWeMT_sigInvMassRebin->GetBinLowEdge(i)<<"<MT<"<<ggMetWeMT_sigInvMassRebin->GetBinLowEdge(i+1)<<" HH_2W2g 130_1 signal Events: "<<h_SMS_HH_2W2g_gg_1Ele_MT_Excluded_Rebin->GetBinContent(i)*h_SMS_HH_2W2g_gg_1Ele_MT_Excluded_Rebin->GetBinWidth(i)<<" +- "<<h_SMS_HH_2W2g_gg_1Ele_MT_Excluded_Rebin->GetBinError(i)*h_SMS_HH_2W2g_gg_1Ele_MT_Excluded_Rebin->GetBinWidth(i)<<endl;
      cout<<ggMetWeMT_sigInvMassRebin->GetBinLowEdge(i)<<"<MT<"<<ggMetWeMT_sigInvMassRebin->GetBinLowEdge(i+1)<<" ZH 130_1 signal Events: "<<h_SMS_ZH_gg_1Ele_MT_Excluded_Rebin->GetBinContent(i)*h_SMS_ZH_gg_1Ele_MT_Excluded_Rebin->GetBinWidth(i)<<endl;
  }

  for(int i=1;i<=ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->GetNbinsX();i++){
      if(i==1)cout<<"W->mu+nu MTerrors:"<<endl;
      Double_t errNewStat(0.),errNewSyst(0.),errNewShape(0.),errSigStat(0.),errHalfDiff(0.);
      Double_t sigStat = ggMetWmuMT_sigInvMassRebin->IntegralAndError(i,i,errSigStat);
      sigStat*=ggMetWmuMT_sigInvMassRebin->GetBinWidth(i);errSigStat*=ggMetWmuMT_sigInvMassRebin->GetBinWidth(i);
      int errBinLow = StatErrsWmuMT->FindBin(ggMetWmuMT_sigInvMassRebin->GetBinLowEdge(i)),errBinHigh = StatErrsWmuMT->FindBin(ggMetWmuMT_sigInvMassRebin->GetBinLowEdge(i+1)-.1);
      Double_t stat=StatErrsWmuMT->IntegralAndError(errBinLow,errBinHigh,errNewStat);
      stat=FitStatSystErrsWmuMT->IntegralAndError(errBinLow,errBinHigh,errNewSyst);
      stat=FitShapeSystErrsWmuMT->IntegralAndError(errBinLow,errBinHigh,errNewShape);
      stat=HalfDiffErrsWmuMT->IntegralAndError(errBinLow,errBinHigh,errHalfDiff);
      Double_t fullErr(0.);
      stat = ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->IntegralAndError(i,i,fullErr);
      stat*=ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->GetBinWidth(i);fullErr*=ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->GetBinWidth(i);
      /*if(ggMetWmuMT_InvMassSBcomb_plus_SMhiggs->GetBinLowEdge(i)<30)*/cout<<ggMetWmuMT_sigInvMassRebin->GetBinLowEdge(i)<<"<MT<"<<ggMetWmuMT_sigInvMassRebin->GetBinLowEdge(i+1)<<" Candidate   Events: "<<sigStat<<" +- "<<errSigStat<<endl;
      cout<<ggMetWmuMT_sigInvMassRebin->GetBinLowEdge(i)<<"<MT<"<<ggMetWmuMT_sigInvMassRebin->GetBinLowEdge(i+1)<<" BG estimate Events: "<<stat<<" +- "<<fullErr<<"  StatErr: "<<errNewStat<<"  FitStatSystErr: "<<errNewSyst<</*"  FitShapeSyst: "<<errNewShape<<*/"  HalfDiffErr: "<<errHalfDiff<<endl;
      cout<<ggMetWmuMT_sigInvMassRebin->GetBinLowEdge(i)<<"<MT<"<<ggMetWmuMT_sigInvMassRebin->GetBinLowEdge(i+1)<<" HH_2W2g 130_1 signal Events: "<<h_SMS_HH_2W2g_gg_1Mu_MT_Excluded_Rebin->GetBinContent(i)*h_SMS_HH_2W2g_gg_1Mu_MT_Excluded_Rebin->GetBinWidth(i)<<" +- "<<h_SMS_HH_2W2g_gg_1Mu_MT_Excluded_Rebin->GetBinError(i)*h_SMS_HH_2W2g_gg_1Mu_MT_Excluded_Rebin->GetBinWidth(i)<<endl;
      cout<<ggMetWmuMT_sigInvMassRebin->GetBinLowEdge(i)<<"<MT<"<<ggMetWmuMT_sigInvMassRebin->GetBinLowEdge(i+1)<<" ZH 130_1 signal Events: "<<h_SMS_ZH_gg_1Mu_MT_Excluded_Rebin->GetBinContent(i)*h_SMS_ZH_gg_1Mu_MT_Excluded_Rebin->GetBinWidth(i)<<endl;
  }



   fin->Close(); 
   fout_We.Close();
   fout_Wmu.Close();
   fout.Close();
   fout_SMHiggs.Close();
   // return;    
}  

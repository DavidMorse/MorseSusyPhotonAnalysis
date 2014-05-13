#define DavePlots_cxx
#include "DavePlots.h"
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

//float L_int = 12162.;
float L_int = 19499.;//full dataset
//float L_int = 19388.;//full dataset
//float L_int = 4041.;//ICHEP dataset
//float ffPercent=0.4372,gammafakePercent=0.5628,eePercent=0.;//these numbers are from jet req min chi2 fraction fit
//float ffPercent=0.228237,eePercent=0.,gfPercent=0.62523,fgPercent=0.146533,gammafakePercent=gfPercent+fgPercent;//these numbers are from purity study
float ffPercent=0.118914,eePercent=0.452198,gfPercent=0.328314,fgPercent=0.100574,gammafakePercent=gfPercent+fgPercent;//these numbers are from purity study
float FakeRateSyst=0.01;
float eeScale=0.,eeScale_JetReq=0.,eeScale_2JetReq=0.,ffScale=0.,ffScale_JetReq=0.,ffScale_2JetReq=0.;
float eeScaleErr=0.,eeScaleErr_JetReq=0.,eeScaleErr_2JetReq=0.,ffScaleErr=0.,ffScaleErr_JetReq=0.,ffScaleErr_2JetReq=0.;
//to get correct numbers for AN should always have the following low bin edges:0,5,10,15,30,50,100
/*Double_t xbins[15]={ 
  0,   //1
  5,   //2
  10,  //3
  15,  //4
  20,  //5
  25,  //6
  30,  //7
  35,  //8
  40,  //9
  //45,  //10
  50,  //11
  //60,  //12
  70,  //13
  //80,  //14
  100, //15
  //120,
  180, //16
  280,
  //250,
  400 //17
  };*/
Double_t xbins[14]={ 
  0,   //1
  5,   //2
  10,  //3
  15,  //4
  20,  //5
  25,  //6
  30,  //7
  35,  //8
  40,  //9
  //45,  //10
  50,  //11
  65,  //12
  80,  //13
  //80,  //14
  100, //15
  //170,  //300,
  250, //16
  //280,
  //250,
  //400 //17
};///ICHEP binning
int NmetBins=(sizeof(xbins)/sizeof(Double_t))-1;
float metPlotXmax = xbins[NmetBins];

Double_t xbinsGGonly[15]={ 
  0,   //1
  5,   //2
  10,  //3
  15,  //4
  20,  //5
  25,  //6
  30,  //7
  35,  //8
  40,  //9
  //45,  //10
  50,  //11
  65,  //12
  80,  //13
  //80,  //14
  100, //15
  //170,  //300,
  250, //16
  //280,
  //250,
  400 //17
};
int NmetBinsGGonly=(sizeof(xbinsGGonly)/sizeof(Double_t))-1;
float metPlotXmaxGGonly = xbinsGGonly[NmetBinsGGonly];

Double_t DiEMPtBins[]={0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100.,110.,120.,130.,140.,150.,200.,300.,400.,600.,1000.};
int numBins=sizeof(DiEMPtBins)/sizeof(Double_t)-1;

void GetFakeRateAndErrorFullRange(float eeVal, float eeErr, float egVal, float egErr, float &fakeRate, float &fakeRateErr, float &fakeRateErrForLimits){

  fakeRate = (egVal / (2*eeVal + egVal));
  float fakeRateErrSquared = ( egVal/( 2*eeVal + egVal ) )*( egVal/( 2*eeVal + egVal ) )*( ( egErr/egVal )*( egErr/egVal ) + ( (4*eeErr*eeErr+egErr*egErr)/((2*eeVal+egVal)*(2*eeVal+egVal)) )  );
  fakeRateErrForLimits=sqrt(fakeRateErrSquared);
  fakeRateErrSquared += FakeRateSyst*FakeRateSyst;
  fakeRateErr=sqrt(fakeRateErrSquared);
  //fakeRate = 0.018;
  //fakeRate = 0.080;
  fakeRate = 0.0209;
  //fakeRate = 0.0199;
  //fakeRateErr = 0.0009;
  //fakeRateErrForLimits = 0.0003;
  fakeRateErr = 0.0006;
  fakeRateErrForLimits = 0.0001;
  return;
}

void GetFakeRateAndError(float eeVal, float eeErr, float egVal, float egErr, float &fakeRate, float &fakeRateErr){

  fakeRate = (egVal / (eeVal + egVal));
  float fakeRateErrSquared = ( egVal/( eeVal + egVal ) )*( egVal/( eeVal + egVal ) )*( ( egErr/egVal )*( egErr/egVal ) + ( (eeErr*eeErr+egErr*egErr)/((eeVal+egVal)*(eeVal+egVal)) )  );
  fakeRateErr=sqrt(fakeRateErrSquared);
  return;
}

float AddInQuadrature(vector<float> entries,int firstEntry,int lastEntry){
  if(lastEntry==-1){
    float sum=0.;
    for(int i=firstEntry-1;i<=(int)entries.size()-1;i++){
      sum+=entries[i]*entries[i];
    }
    float err=sqrt(sum);
    return err;
  }
  else if((int)entries.size()>=lastEntry){
    float sum=0.;
    for(int i=firstEntry-1;i<=lastEntry-1;i++){
      sum+=entries[i]*entries[i];
    }
    float err=sqrt(sum);
    return err;
  }
  else return 0.;
}
float AddUp(vector<float> entries,int firstEntry,int lastEntry){
  float sum=0.; 
  if(lastEntry==-1){
    for(int i=firstEntry-1;i<=(int)entries.size()-1;i++){
      sum+=entries[i];
    }
    //float err=sqrt(sum);
    return sum;
  }
  else if((int)entries.size()>=lastEntry){
    for(int i=firstEntry-1;i<=lastEntry-1;i++){
      sum+=entries[i];
    }
    //float err=sqrt(sum);
    return sum;
  }
  else return 0.;
}
/*
  float AddInQuadrature(TH1F* ErrHist,int firstBinLowEdge,int lastBinHighEdge){
  if(lastBin==-1){
  float sum=0.;
  int firstBin=ErrHist->FindBin(firstBinLowEdge);
  int lastBin=ErrHist->FindBin(lastBinHighEdge-1);
  for(int i=firstBin;i<=ErrHist->GetNbinsX()+1;i++){
  sum+=ErrHist->GetBinContent[i]*ErrHist->GetBinContent[i];
  }
  float err=sqrt(sum);
  return err;
  }
  else if(ErrHist->GetNbinsX()>=lastBin){
  float sum=0.;
  for(int i=firstBin;i<=lastBin;i++){
  sum+=ErrHist->GetBinContent[i]*ErrHist->GetBinContent[i];
  }
  float err=sqrt(sum);
  return err;
  }
  else return 0.;
  }
*/
void GetAndSetFFerrors(TH1F* ggHist,TH1F* egHist,TH1F* eeHist,/*TH1F *&eeHistTemp,*/ TH1F *&ffHist,TH1F *&ffHistRebin,TH1F *&ffHistForLimits,TH1F *&ffHistForLimitsStatOnly,float fakerate,float fakerateErr,float fakerateErrForLimits,float STDDEV[],vector<float> &FFnormErr,vector<float> &FFreweightErr,vector<float> &FFstatErr,vector<float> &DeltaEEerr,vector<float> &DeltaFFGFerr,bool useffgferr,float &ffscale, float &ffscaleErr){

  cout<<"at beginning of function"<<endl;

  double eg04Err=0.,eg04=0.,gg04Err=0.,gg04=0.,ff04Err=0.,ff04=0.,ee04=0.,ee04Err=0.,scaleEE=0.,scaleEEErr=0.,egPrime04=0.,egPrime04Err=0.,egPrime04ErrForLimits=0.,scaleFF=0.,scaleFFErr=0.,scaleFFErrSquared=0.,scaleFFErrForLimits=0.,scaleFFErrSquaredForLimits=0.,Nbin=0.,NbinErr=0.,NbinErrSquared=0.,binErr=0.,binErrForLimits=0.,binErrSquared=0.,reweightErr=0.,reweightErrSquared=0.,deltaEEerr=0.,eePrime04=0.,eetemp=0.,eetempErr=0.,eescale;

  FFnormErr.clear();FFreweightErr.clear();FFstatErr.clear();DeltaEEerr.clear();

  gg04=ggHist->IntegralAndError(0,4,gg04Err);
  ff04=ffHist->IntegralAndError(0,4,ff04Err);
  ee04=eeHist->IntegralAndError(0,4,ee04Err);

  eePrime04=fakerate*fakerate*ee04;

  eg04=egHist->IntegralAndError(0,4,eg04Err);
  egPrime04=/*(fakerate/(1-fakerate))**/eg04;
  //egPrime04=fakerate*(eg04-eePrime04);
  //egPrime04Err=sqrt(eg04*eg04*fakerateErr*fakerateErr + fakerate*fakerate*eg04Err*eg04Err);
  //egPrime04ErrForLimits=sqrt(eg04*eg04*fakerateErrForLimits*fakerateErrForLimits + fakerate*fakerate*eg04Err*eg04Err);
  egPrime04Err=egPrime04*sqrt(fakerateErr*fakerateErr*(1./(fakerate*fakerate)+1./((1-fakerate)*(1-fakerate)))+(eg04Err*eg04Err)/(eg04*eg04));
  egPrime04ErrForLimits=egPrime04*sqrt(fakerateErrForLimits*fakerateErrForLimits*(1./(fakerate*fakerate)+1./((1-fakerate)*(1-fakerate)))+(eg04Err*eg04Err)/(eg04*eg04));


  scaleFF=(gg04-egPrime04)/ff04;
  //scaleFF=(gg04-egPrime04-eePrime04)/ff04;
  scaleEE=(gg04-egPrime04)/ee04;
  //scaleEEErr=scaleEE*sqrt( (gg04Err*gg04Err+egPrime04Err*egPrime04Err)/((gg04-egPrime04)*(gg04-egPrime04)) + (ee04Err*ee04Err)/(ee04*ee04) );
  scaleEEErr=scaleEE*sqrt( (gg04Err*gg04Err)/(eePrime04*eePrime04) + (gg04*gg04)*(ee04Err*ee04Err)/(eePrime04*eePrime04) );
  scaleFFErr=scaleFF*sqrt(  ( sqrt(gg04Err*gg04Err+egPrime04Err*egPrime04Err)/((gg04-egPrime04)*(gg04-egPrime04)) )  +  (ff04Err*ff04Err)/(ff04*ff04)  );
  scaleFFErrForLimits=scaleFF*sqrt(  ( sqrt(gg04Err*gg04Err+egPrime04ErrForLimits*egPrime04ErrForLimits)/((gg04-egPrime04)*(gg04-egPrime04)) )  +  (ff04Err*ff04Err)/(ff04*ff04)  );
  ffscale=scaleFF;ffscaleErr=scaleFFErr;

  for(int i=1;i<ffHist->GetNbinsX()+1;i++){
    double ee=0.,eeErr=0.,eePrime=0.,eeErrFull=0.,eeErrFullForLimits=0.;
    ee=eeHist->IntegralAndError(i,i,eeErr);
    eePrime=ee*fakerate*fakerate;
    //eeHist->SetBinContent(i,eePrime);
    eeErrFull=sqrt(fakerate*fakerate*fakerate*fakerate*eeErr*eeErr + 2*ee*ee*fakerate*fakerate*fakerateErr*fakerateErr);
    //cout<<"Bin: "<<i<<" eeVal: "<<eePrime<<" eeErr: "<<eeErrFull<<endl;
    //eeHist->SetBinError(i,eeErrFull);
    Nbin=ffHist->IntegralAndError(i,i,NbinErr);
    reweightErr=std::fabs(STDDEV[i-1]);
    float ScaledFFbin=Nbin*scaleFF;
    float ScaledEEbin=eeHist->GetBinContent(i)*scaleEE;
    if(useffgferr && DeltaFFGFerr.size()>i-1)deltaEEerr=DeltaFFGFerr[i-1];
    else deltaEEerr=std::fabs(ScaledFFbin - ScaledEEbin);
    //cout<<ScaledFFbin<<endl<<eeHist->GetBinContent(i)<<endl<<deltaEEerr<<endl<<endl;
    scaleFFErrSquared=ScaledFFbin*ScaledFFbin*scaleFFErr*scaleFFErr/(scaleFF*scaleFF);
    scaleFFErrSquaredForLimits=ScaledFFbin*ScaledFFbin*scaleFFErrForLimits*scaleFFErrForLimits/(scaleFF*scaleFF);
    if(Nbin==0){reweightErrSquared=0;NbinErrSquared=0;}
    else{
      reweightErrSquared=ScaledFFbin*ScaledFFbin*reweightErr*reweightErr/(Nbin*Nbin);
      NbinErrSquared=ScaledFFbin*ScaledFFbin*NbinErr*NbinErr/(Nbin*Nbin);
    }
    binErrForLimits=sqrt(reweightErrSquared + NbinErrSquared + scaleFFErrSquaredForLimits);
    binErrSquared=reweightErrSquared + NbinErrSquared + scaleFFErrSquared/* + deltaEEerr*deltaEEerr*/;//deltaEE error needs to be applied after all the others
    binErr=sqrt(binErrSquared);
    //cout<<reweightErr<<endl<<NbinErr<<endl<<scaleFFErr<<endl<<scaleFF<<endl<<deltaEEerr<<endl<<sqrt( (reweightErr*reweightErr + NbinErr*NbinErr)/(Nbin*Nbin) + (scaleFFErr*scaleFFErr)/(scaleFF*scaleFF) + deltaEEerr*deltaEEerr )<<endl<<binErr<<endl<<endl;
    //cout<<scaleFFErr<<endl;
    FFnormErr.push_back(sqrt(scaleFFErrSquared));
    FFreweightErr.push_back(sqrt(reweightErrSquared));
    FFstatErr.push_back(sqrt(NbinErrSquared));
    if(!useffgferr)DeltaEEerr.push_back(deltaEEerr);
    else DeltaFFGFerr.push_back(deltaEEerr);
    ffHist->SetBinContent(i,Nbin*scaleFF);
    ffHist->SetBinError(i,binErr);
    ffHistForLimits->SetBinContent(i,Nbin*scaleFF);
    ffHistForLimits->SetBinError(i,binErrForLimits);
    ffHistForLimitsStatOnly->SetBinContent(i,Nbin*scaleFF);
    ffHistForLimitsStatOnly->SetBinError(i,sqrt(NbinErrSquared + scaleFFErrSquared));
    eeHist->SetBinContent(i,eeHist->GetBinContent(i)*scaleEE);
    eeHist->SetBinError(i,eeHist->GetBinError(i)*scaleEE);
  }
  TH1F* eeHistTemp2 = (TH1F*)eeHist->Rebin(NmetBins,"eeHistTemp2",xbins);
  //cout<<"OldOverFlow: "<<eeHist->Integral(eeHist->FindBin(275.1),-1)<<"  NewOverFlow: "<<eeHistTemp->Integral(eeHistTemp->FindBin(275.1),-1)<<endl;
  ffHistRebin = (TH1F*)ffHist->Rebin(NmetBins,ffHistRebin->GetName(),xbins);
  //include overflow in last bin
  Double_t ffOverFlow=0.,ffOverFlowErr=0.;int overflowBinff=ffHist->FindBin(ffHistRebin->GetBinLowEdge(NmetBins+1)+.01);
  ffOverFlow=ffHist->IntegralAndError(overflowBinff,999,ffOverFlowErr);
  //cout<<"ffOverFlowBinff:"<<overflowBinff<<" lowEdge:"<<ffHist->GetBinLowEdge(overflowBinff)<<" ffOverFlow:"<<ffOverFlow<<" ffOverFlowErr:"<<ffOverFlow<<endl;
  //cout<<"ffHistRebin->GetMaximumBin():"<<ffHistRebin->GetMaximumBin()<<"ffHistRebin->GetNbinsX():"<<ffHistRebin->GetNbinsX()<<" ffHistRebin->GetBinLowEdge(ffHistRebin->GetMaximumBin()+1):"<<ffHistRebin->GetBinLowEdge(ffHistRebin->GetMaximumBin()+1)<<" ffHist->FindBin(ffHistRebin->GetBinLowEdge(ffHistRebin->GetMaximumBin()+1)+.01)"<<ffHist->FindBin(ffHistRebin->GetBinLowEdge(ffHistRebin->GetMaximumBin()+1)+.01)<<endl;
  float ffLastBinNew=ffHistRebin->GetBinContent(NmetBins)+ffOverFlow;
  float ffLastBinNewErr=sqrt(ffHistRebin->GetBinError(NmetBins)*ffHistRebin->GetBinError(NmetBins)+ffOverFlowErr*ffOverFlowErr);
  Double_t eeOverFlow=0.,eeOverFlowErr=0.;int overflowBinee=eeHist->FindBin(eeHistTemp2->GetBinLowEdge(NmetBins+1)+.01);
  eeOverFlow=eeHist->IntegralAndError(overflowBinee,999,eeOverFlowErr);
  float eeLastBinNew=eeHistTemp2->GetBinContent(NmetBins)+eeOverFlow;
  float eeLastBinNewErr=sqrt(eeHistTemp2->GetBinError(NmetBins)*eeHistTemp2->GetBinError(NmetBins)+eeOverFlowErr*eeOverFlowErr);

  ffHistRebin->SetBinContent(NmetBins,ffLastBinNew);
  ffHistRebin->SetBinError(NmetBins,ffLastBinNewErr);
  eeHistTemp2->SetBinContent(NmetBins,eeLastBinNew);
  eeHistTemp2->SetBinError(NmetBins,eeLastBinNewErr);
  for(int i=1;i<ffHistRebin->GetNbinsX()+2;i++){
    float errRaw = ffHistRebin->GetBinError(i);
    float errRawSquared=errRaw*errRaw;
    float diff = 0.;
    int which=0;
    if(!useffgferr){diff=std::fabs(eeHistTemp2->GetBinContent(i)-ffHistRebin->GetBinContent(i));which=1;}
    else if(DeltaFFGFerr.size()>i){diff=DeltaFFGFerr[i-1];which=2;}
    else{diff=0;which=3;}
    cout<<"bin:"<<i<<" of:"<<ffHistRebin->GetNbinsX()+2<<" method:"<<which<<" errRaw="<<errRaw<<" Diff:"<<diff<<endl;
    float errNewSquared=errRawSquared+diff*diff;
    float errNew=sqrt(errNewSquared);
    ffHistRebin->SetBinError(i,errNew);
  }
  cout<<"at end of function"<<endl;
  return;
}

void GetAndSetEEerrors(TH1F* ggHist,TH1F *&egHist,TH1F *&egHistForLimits,TH1F *&eeHist,TH1F *&eeHistRebin,TH1F *&eeHistForLimits,TH1F *&eeHistForLimitsStatOnly,TH1F *eeSBhighHist,TH1F *eeSBlowHist,TH1F *ffHist,float fakerate,float fakerateErr,float fakerateErrForLimits,float STDDEV[],float STDDEVsbh[],float STDDEVsbl[],vector<float> &EEnormErr,vector<float> &EEreweightErr,vector<float> &EEstatErr,vector<float> &DeltaFFerr,float &eescale, float &eescaleErr, vector<float> &egStatErr, vector<float> &egNormErr, TH1F* eeMetFromZZ, TH1F* eeMetFromWZ){

  double eg04Err=0.,eg04=0.,gg04Err=0.,gg04=0.,ee04Err=0.,ee04=0.,ff04=0.,eeSBhigh04Err=0.,eeSBhigh04=0.,eeSBlow04Err=0.,eeSBlow04=0.,egPrime04=0.,egPrime04Err=0.,egPrime04ErrForLimits=0.,scaleFF=0.,scaleEE=0.,scaleEEErr=0.,scaleEEErrSquaredForLimits=0.,scaleEEErrForLimits=0.,scaleEEErrSquared=0.,Nbin=0.,NbinErr=0.,NbinErrFullSquared=0.,NbinSBh=0.,NbinSBhErr=0.,NbinSBl=0.,NbinSBlErr=0.,binErr=0.,binErrForLimits=0.,binErrSquared=0.,reweightErr=0.,reweightErrFullSquared=0.,reweightErrSBh=0.,reweightErrSBl=0.,deltaFFerr=0.;
 
  EEnormErr.clear();EEreweightErr.clear();EEstatErr.clear();DeltaFFerr.clear();egStatErr.clear();egNormErr.clear();

  gg04       =ggHist->IntegralAndError(0,4,gg04Err);
  ee04       =eeHist->IntegralAndError(0,4,ee04Err);
  ff04       =ffHist->Integral(0,4);
  eeSBhigh04 =eeSBhighHist->IntegralAndError(0,4,eeSBhigh04Err);
  eeSBlow04  =eeSBlowHist->IntegralAndError(0,4,eeSBlow04Err);
  float eeZZ04 =eeMetFromZZ->Integral(0,4);
  float eeWZ04 =eeMetFromWZ->Integral(0,4);
  float eeMinusSB04=ee04-eeSBhigh04-eeSBlow04-eeZZ04-eeWZ04;
  eg04=egHist->IntegralAndError(0,4,eg04Err);
  egPrime04=/*(fakerate/(1-fakerate))**/eg04;
  //egPrime04=fakerate*eg04;
  //egPrime04Err=sqrt(eg04*eg04*fakerateErr*fakerateErr + fakerate*fakerate*eg04Err*eg04Err);
  //egPrime04ErrForLimits=sqrt(eg04*eg04*fakerateErrForLimits*fakerateErrForLimits + fakerate*fakerate*eg04Err*eg04Err);
  egPrime04Err=egPrime04*sqrt(fakerateErr*fakerateErr*(1./(fakerate*fakerate)+1./((1-fakerate)*(1-fakerate)))+(eg04Err*eg04Err)/(eg04*eg04));
  egPrime04ErrForLimits=egPrime04*sqrt(fakerateErrForLimits*fakerateErrForLimits*(1./(fakerate*fakerate)+1./((1-fakerate)*(1-fakerate)))+(eg04Err*eg04Err)/(eg04*eg04));

  scaleEE=(gg04-egPrime04)/(ee04-eeSBhigh04-eeSBlow04-eeZZ04-eeWZ04);  
  scaleFF=(gg04-egPrime04)/(ff04);
  scaleEEErr=scaleEE*sqrt(  ( sqrt(gg04Err*gg04Err+egPrime04Err*egPrime04Err)/((gg04-egPrime04)*(gg04-egPrime04)) )  +  (ee04Err*ee04Err+eeSBhigh04Err*eeSBhigh04Err+eeSBlow04Err*eeSBlow04Err)/(eeMinusSB04*eeMinusSB04)  );
  scaleEEErrForLimits=scaleEE*sqrt(  ( sqrt(gg04Err*gg04Err+egPrime04ErrForLimits*egPrime04ErrForLimits)/((gg04-egPrime04)*(gg04-egPrime04)) )  +  (ee04Err*ee04Err+eeSBhigh04Err*eeSBhigh04Err+eeSBlow04Err*eeSBlow04Err)/(eeMinusSB04*eeMinusSB04)  );
  eescale=scaleEE;eescaleErr=scaleEEErr;

  for(int i=1;i<eeHist->GetNbinsX()+1;i++){
    double eg=0.,egErr=0.,egPrime=0.,egErrFull=0.,egErrFullForLimits=0.;
    eg=egHist->IntegralAndError(i,i,egErr);
    egPrime=eg/**(fakerate/(1-fakerate))*/;
    egHist->SetBinContent(i,egPrime);
    egHistForLimits->SetBinContent(i,egPrime);
    //egErrFull=sqrt(fakerate*egErr*fakerate*egErr + fakerateErr*eg*fakerateErr*eg);//this is using only egPrime=fakerate*eg
    //egErrFullForLimits=sqrt(fakerate*egErr*fakerate*egErr + fakerateErrForLimits*eg*fakerateErr*eg);
    //cout<<"bin:"<<i<<"  egPrime:"<<egPrime<<"  fakerate:"<<fakerate<<"  fakerateErr:"<<fakerateErr<<"  eg:"<<eg<<"  egErr:"<<egErr<<endl;
    if(eg==0){egErrFull=0;egErrFullForLimits=0;}//need this or error blows up
    else{egErrFull=egPrime*sqrt(fakerateErr*fakerateErr*(1./(fakerate*fakerate)+1./((1-fakerate)*(1-fakerate)))+(egErr*egErr)/(eg*eg));
      egErrFullForLimits=egPrime*sqrt(fakerateErrForLimits*fakerateErrForLimits*(1./(fakerate*fakerate)+1./((1-fakerate)*(1-fakerate)))+(egErr*egErr)/(eg*eg));}
    egHist->SetBinError(i,egErrFull);
    egHistForLimits->SetBinError(i,egErrFullForLimits);
    egStatErr.push_back(fakerate*egErr);
    egNormErr.push_back(fakerateErr*eg);
    Nbin=eeHist->IntegralAndError(i,i,NbinErr);
    NbinSBh=eeSBhighHist->IntegralAndError(i,i,NbinSBhErr);
    NbinSBl=eeSBlowHist->IntegralAndError(i,i,NbinSBlErr);
    double NbinZZErr=0.,NbinWZErr=0.;
    float NbinZZ=eeMetFromZZ->IntegralAndError(i,i,NbinZZErr);
    float NbinWZ=eeMetFromWZ->IntegralAndError(i,i,NbinWZErr);
    float NbinErrFull=sqrt(NbinErr*NbinErr+NbinSBhErr*NbinSBhErr+NbinSBlErr*NbinSBlErr+NbinZZErr*NbinZZErr+NbinWZErr*NbinWZErr);
    //cout<<NbinErr<<endl<<NbinSBhErr<<endl<<NbinSBlErr<<endl<<NbinErrFull<<endl<<endl;
    float NbinEEminusSB=Nbin-NbinSBh-NbinSBl;
    float NbinEEminusSBandZZandWZ=NbinEEminusSB-NbinZZ-NbinWZ;
    reweightErr=std::fabs(STDDEV[i-1]);
    reweightErrSBh=std::fabs(STDDEVsbh[i-1]);
    reweightErrSBl=std::fabs(STDDEVsbl[i-1]);
    float reweightErrFull=sqrt(reweightErr*reweightErr+reweightErrSBh*reweightErrSBh+reweightErrSBl*reweightErrSBl);
    float ScaledEEbin=NbinEEminusSBandZZandWZ*scaleEE;
    float ScaledFFbin=ffHist->GetBinContent(i)*scaleFF;
    deltaFFerr=std::fabs(ScaledEEbin - ScaledFFbin);
    //cout<<ScaledEEbin<<endl<<eeHist->GetBinContent(i)<<endl<<deltaEEerr<<endl<<endl;
    scaleEEErrSquared=ScaledEEbin*ScaledEEbin*scaleEEErr*scaleEEErr/(scaleEE*scaleEE);
    scaleEEErrSquaredForLimits=ScaledEEbin*ScaledEEbin*scaleEEErrForLimits*scaleEEErrForLimits/(scaleEE*scaleEE);
    if(NbinEEminusSBandZZandWZ==0){reweightErrFullSquared=0.;NbinErrFull=0.;}
    else{
      reweightErrFullSquared=ScaledEEbin*ScaledEEbin*reweightErrFull*reweightErrFull/(NbinEEminusSBandZZandWZ*NbinEEminusSBandZZandWZ);
      NbinErrFullSquared=ScaledEEbin*ScaledEEbin*NbinErrFull*NbinErrFull/(NbinEEminusSBandZZandWZ*NbinEEminusSBandZZandWZ);
    }
    binErrSquared=reweightErrFullSquared + NbinErrFullSquared + scaleEEErrSquared/* + deltaFFerr*deltaFFerr*/;//deltaFF error needs to be applied after all the others
    binErr= sqrt(binErrSquared);
    binErrForLimits=sqrt(reweightErrFullSquared + NbinErrFullSquared + scaleEEErrSquaredForLimits);
    //cout<<reweightErrFull<<endl<<NbinErrFull<<endl<<scaleEEErr<<endl<<scaleEE<<endl<<deltaFFerr<<endl<<binErr<<endl<<endl;

    EEnormErr.push_back(sqrt(scaleEEErrSquared));
    EEreweightErr.push_back(sqrt(reweightErrFullSquared));
    EEstatErr.push_back(sqrt(NbinErrFullSquared));
    DeltaFFerr.push_back(deltaFFerr);
    eeHist->SetBinContent(i,NbinEEminusSBandZZandWZ*scaleEE);
    eeHist->SetBinError(i,binErr);
    eeHistForLimits->SetBinContent(i,NbinEEminusSBandZZandWZ*scaleEE);
    eeHistForLimits->SetBinError(i,binErrForLimits);
    eeHistForLimitsStatOnly->SetBinContent(i,NbinEEminusSBandZZandWZ*scaleEE);
    eeHistForLimitsStatOnly->SetBinError(i,sqrt(NbinErrFullSquared+scaleEEErrSquared));
  }
  ffHist->Scale(scaleFF);
  TH1F* ffHistTemp = (TH1F*)ffHist->Rebin(NmetBins,"ffHistTemp",xbins);
  eeHistRebin = (TH1F*)eeHist->Rebin(NmetBins,eeHistRebin->GetName(),xbins);
  //include overflow in last bin
  Double_t eeOverFlow=0.,eeOverFlowErr=0.;int overflowBin=eeHist->FindBin(eeHistRebin->GetBinLowEdge(NmetBins+1)+.01);
  eeOverFlow=eeHist->IntegralAndError(overflowBin,999,eeOverFlowErr);
  float eeLastBinNew=eeHistRebin->GetBinContent(NmetBins)+eeOverFlow;
  float eeLastBinNewErr=sqrt(eeHistRebin->GetBinError(NmetBins)*eeHistRebin->GetBinError(NmetBins)+eeOverFlowErr*eeOverFlowErr);
  Double_t ffOverFlow=0.,ffOverFlowErr=0.;int overflowBinff=ffHist->FindBin(ffHistTemp->GetBinLowEdge(NmetBins+1)+.01);
  ffOverFlow=ffHist->IntegralAndError(overflowBinff,999,ffOverFlowErr);
  float ffLastBinNew=ffHistTemp->GetBinContent(NmetBins)+ffOverFlow;
  float ffLastBinNewErr=sqrt(ffHistTemp->GetBinError(NmetBins)*ffHistTemp->GetBinError(NmetBins)+ffOverFlowErr*ffOverFlowErr);
  for(int i=1;i<eeHistRebin->GetNbinsX()+2;i++){
    float errRaw = eeHistRebin->GetBinError(i);
    float errRawSquared=errRaw*errRaw;
    float diff = std::fabs(ffHistTemp->GetBinContent(i)-eeHistRebin->GetBinContent(i));
    //cout<<"bin:"<<i<<" errRaw="<<errRaw<<" Diff:"<<diff<<endl;
    float errNewSquared=errRawSquared+diff*diff;
    float errNew=sqrt(errNewSquared);
    //eeHistRebin->SetBinError(i,errNew);
  }
  return;
}

double gg0_20=0.,gg0_20Error=0.,eg0_20=0.,eg0_20Error=0.,ee0_20=0.,ee0_20Error=0.,ff0_20=0.,ff0_20Error=0.,QCDee0_20=0.,QCDee0_20Error=0.,QCDff0_20=0.,QCDff0_20Error=0.;
double gg30_50=0.,gg30_50Error=0.,eg30_50=0.,eg30_50Error=0.,ee30_50=0.,ee30_50Error=0.,ff30_50=0.,ff30_50Error=0.,QCDee30_50=0.,QCDee30_50Error=0.,QCDff30_50=0.,QCDff30_50Error=0.;
double gg50up=0.,gg50upError=0.,eg50up=0.,eg50upError=0.,ee50up=0.,ee50upError=0.,ff50up=0.,ff50upError=0.,QCDee50up=0.,QCDee50upError=0.,QCDff50up=0.,QCDff50upError=0.;
double gg100up=0.,gg100upError=0.,eg100up=0.,eg100upError=0.,ee100up=0.,ee100upError=0.,ff100up=0.,ff100upError=0.,QCDee100up=0.,QCDee100upError=0.,QCDff100up=0.,QCDff100upError=0.,comb100up=0.;
double gg0_20_JetReq=0.,gg0_20Error_JetReq=0.,eg0_20_JetReq=0.,eg0_20Error_JetReq=0.,ee0_20_JetReq=0.,ee0_20Error_JetReq=0.,ff0_20_JetReq=0.,ff0_20Error_JetReq=0.,QCDee0_20_JetReq=0.,QCDee0_20Error_JetReq=0.,QCDff0_20_JetReq=0.,QCDff0_20Error_JetReq=0.;
double gg30_50_JetReq=0.,gg30_50Error_JetReq=0.,eg30_50_JetReq=0.,eg30_50Error_JetReq=0.,ee30_50_JetReq=0.,ee30_50Error_JetReq=0.,ff30_50_JetReq=0.,ff30_50Error_JetReq=0.,QCDee30_50_JetReq=0.,QCDee30_50Error_JetReq=0.,QCDff30_50_JetReq=0.,QCDff30_50Error_JetReq=0.;
double gg50up_JetReq=0.,gg50upError_JetReq=0.,eg50up_JetReq=0.,eg50upError_JetReq=0.,ee50up_JetReq=0.,ee50upError_JetReq=0.,ff50up_JetReq=0.,ff50upError_JetReq=0.,QCDee50up_JetReq=0.,QCDee50upError_JetReq=0.,QCDff50up_JetReq=0.,QCDff50upError_JetReq=0.;
double gg100up_JetReq=0.,gg100upError_JetReq=0.,eg100up_JetReq=0.,eg100upError_JetReq=0.,ee100up_JetReq=0.,ee100upError_JetReq=0.,ff100up_JetReq=0.,ff100upError_JetReq=0.,QCDee100up_JetReq=0.,QCDee100upError_JetReq=0.,QCDff100up_JetReq=0.,QCDff100upError_JetReq=0.,comb100up_JetReq=0.;
double gg50_60_JetReq=0.,gg50_60Error_JetReq=0.,eg50_60_JetReq=0.,eg50_60Error_JetReq=0.,ee50_60_JetReq=0.,ee50_60Error_JetReq=0.,ff50_60_JetReq=0.,ff50_60Error_JetReq=0.,QCDee50_60_JetReq=0.,QCDee50_60Error_JetReq=0.,QCDff50_60_JetReq=0.,QCDff50_60Error_JetReq=0.;
double gg60_70_JetReq=0.,gg60_70Error_JetReq=0.,eg60_70_JetReq=0.,eg60_70Error_JetReq=0.,ee60_70_JetReq=0.,ee60_70Error_JetReq=0.,ff60_70_JetReq=0.,ff60_70Error_JetReq=0.,QCDee60_70_JetReq=0.,QCDee60_70Error_JetReq=0.,QCDff60_70_JetReq=0.,QCDff60_70Error_JetReq=0.;
double gg70_80_JetReq=0.,gg70_80Error_JetReq=0.,eg70_80_JetReq=0.,eg70_80Error_JetReq=0.,ee70_80_JetReq=0.,ee70_80Error_JetReq=0.,ff70_80_JetReq=0.,ff70_80Error_JetReq=0.,QCDee70_80_JetReq=0.,QCDee70_80Error_JetReq=0.,QCDff70_80_JetReq=0.,QCDff70_80Error_JetReq=0.;
double gg80_100_JetReq=0.,gg80_100Error_JetReq=0.,eg80_100_JetReq=0.,eg80_100Error_JetReq=0.,ee80_100_JetReq=0.,ee80_100Error_JetReq=0.,ff80_100_JetReq=0.,ff80_100Error_JetReq=0.,QCDee80_100_JetReq=0.,QCDee80_100Error_JetReq=0.,QCDff80_100_JetReq=0.,QCDff80_100Error_JetReq=0.;
double gg0_20_2JetReq=0.,gg0_20Error_2JetReq=0.,eg0_20_2JetReq=0.,eg0_20Error_2JetReq=0.,ee0_20_2JetReq=0.,ee0_20Error_2JetReq=0.,ff0_20_2JetReq=0.,ff0_20Error_2JetReq=0.,QCDee0_20_2JetReq=0.,QCDee0_20Error_2JetReq=0.,QCDff0_20_2JetReq=0.,QCDff0_20Error_2JetReq=0.;
double gg30_50_2JetReq=0.,gg30_50Error_2JetReq=0.,eg30_50_2JetReq=0.,eg30_50Error_2JetReq=0.,ee30_50_2JetReq=0.,ee30_50Error_2JetReq=0.,ff30_50_2JetReq=0.,ff30_50Error_2JetReq=0.,QCDee30_50_2JetReq=0.,QCDee30_50Error_2JetReq=0.,QCDff30_50_2JetReq=0.,QCDff30_50Error_2JetReq=0.;
double gg50up_2JetReq=0.,gg50upError_2JetReq=0.,eg50up_2JetReq=0.,eg50upError_2JetReq=0.,ee50up_2JetReq=0.,ee50upError_2JetReq=0.,ff50up_2JetReq=0.,ff50upError_2JetReq=0.,QCDee50up_2JetReq=0.,QCDee50upError_2JetReq=0.,QCDff50up_2JetReq=0.,QCDff50upError_2JetReq=0.;
double gg100up_2JetReq=0.,gg100upError_2JetReq=0.,eg100up_2JetReq=0.,eg100upError_2JetReq=0.,ee100up_2JetReq=0.,ee100upError_2JetReq=0.,ff100up_2JetReq=0.,ff100upError_2JetReq=0.,QCDee100up_2JetReq=0.,QCDee100upError_2JetReq=0.,QCDff100up_2JetReq=0.,QCDff100upError_2JetReq=0.,comb100up_2JetReq=0.;
double gg50_60_2JetReq=0.,gg50_60Error_2JetReq=0.,eg50_60_2JetReq=0.,eg50_60Error_2JetReq=0.,ee50_60_2JetReq=0.,ee50_60Error_2JetReq=0.,ff50_60_2JetReq=0.,ff50_60Error_2JetReq=0.,QCDee50_60_2JetReq=0.,QCDee50_60Error_2JetReq=0.,QCDff50_60_2JetReq=0.,QCDff50_60Error_2JetReq=0.;
double gg60_70_2JetReq=0.,gg60_70Error_2JetReq=0.,eg60_70_2JetReq=0.,eg60_70Error_2JetReq=0.,ee60_70_2JetReq=0.,ee60_70Error_2JetReq=0.,ff60_70_2JetReq=0.,ff60_70Error_2JetReq=0.,QCDee60_70_2JetReq=0.,QCDee60_70Error_2JetReq=0.,QCDff60_70_2JetReq=0.,QCDff60_70Error_2JetReq=0.;
double gg70_80_2JetReq=0.,gg70_80Error_2JetReq=0.,eg70_80_2JetReq=0.,eg70_80Error_2JetReq=0.,ee70_80_2JetReq=0.,ee70_80Error_2JetReq=0.,ff70_80_2JetReq=0.,ff70_80Error_2JetReq=0.,QCDee70_80_2JetReq=0.,QCDee70_80Error_2JetReq=0.,QCDff70_80_2JetReq=0.,QCDff70_80Error_2JetReq=0.;
double gg80_100_2JetReq=0.,gg80_100Error_2JetReq=0.,eg80_100_2JetReq=0.,eg80_100Error_2JetReq=0.,ee80_100_2JetReq=0.,ee80_100Error_2JetReq=0.,ff80_100_2JetReq=0.,ff80_100Error_2JetReq=0.,QCDee80_100_2JetReq=0.,QCDee80_100Error_2JetReq=0.,QCDff80_100_2JetReq=0.,QCDff80_100Error_2JetReq=0.;
double gg50_60=0.,gg50_60Error=0.,eg50_60=0.,eg50_60Error=0.,ee50_60=0.,ee50_60Error=0.,ff50_60=0.,ff50_60Error=0.,QCDee50_60=0.,QCDee50_60Error=0.,QCDff50_60=0.,QCDff50_60Error=0.;
double gg60_70=0.,gg60_70Error=0.,eg60_70=0.,eg60_70Error=0.,ee60_70=0.,ee60_70Error=0.,ff60_70=0.,ff60_70Error=0.,QCDee60_70=0.,QCDee60_70Error=0.,QCDff60_70=0.,QCDff60_70Error=0.;
double gg70_80=0.,gg70_80Error=0.,eg70_80=0.,eg70_80Error=0.,ee70_80=0.,ee70_80Error=0.,ff70_80=0.,ff70_80Error=0.,QCDee70_80=0.,QCDee70_80Error=0.,QCDff70_80=0.,QCDff70_80Error=0.;
double gg80_100=0.,gg80_100Error=0.,eg80_100=0.,eg80_100Error=0.,ee80_100=0.,ee80_100Error=0.,ff80_100=0.,ff80_100Error=0.,QCDee80_100=0.,QCDee80_100Error=0.,QCDff80_100=0.,QCDff80_100Error=0.;

float CBsigEE=0.,CBsigEG=0.,CBsigGG=0.,CBsigEE25to40=0.,CBsigEE40to45=0.,CBsigEE45to50=0.,CBsigEE50to60=0.,CBsigEE60to80=0.,CBsigEE80=0.,CBsigEG25to40=0.,CBsigEG40to45=0.,CBsigEG45to50=0.,CBsigEG50to60=0.,CBsigEG60to80=0.,CBsigEG80=0./*,CBsigEE_JetReq=0.,CBsigEG_JetReq=0.*/;
float CBsigEEerror=0.,CBsigEGerror=0.,CBsigGGerror=0.,CBsigEEerror25to40=0.,CBsigEEerror40to45=0.,CBsigEEerror45to50=0.,CBsigEEerror50to60=0.,CBsigEEerror60to80=0.,CBsigEEerror80=0.,CBsigEGerror25to40=0.,CBsigEGerror40to45=0.,CBsigEGerror45to50=0.,CBsigEGerror50to60=0.,CBsigEGerror60to80=0.,CBsigEGerror80=0./*,CBsigEEerror_JetReq=0.,CBsigEGerror_JetReq=0.*/;

vector<float> reweightErree,reweightErrff,normErreg,statErree,statErrff,statErreg,normErree,normErrff;
vector<float> reweightErree_JetReq,reweightErrff_JetReq,normErreg_JetReq,statErree_JetReq,statErrff_JetReq,statErreg_JetReq,normErree_JetReq,normErrff_JetReq;
vector<float> reweightErree_2JetReq,reweightErrff_2JetReq,normErreg_2JetReq,statErree_2JetReq,statErrff_2JetReq,statErreg_2JetReq,normErree_2JetReq,normErrff_2JetReq;
vector<float> diffFromeeErrorff,diffFromeeErrorff_JetReq,diffFromeeErrorff_2JetReq,diffFromffErroree,diffFromffErroree_JetReq,diffFromffErroree_2JetReq,diffFromfgErrorff,diffFromfgErrorff_JetReq,diffFromfgErrorff_2JetReq;

TH1F *totalBGee, *totalBGee_JetReq, *totalBGee_2JetReq, *totalBGff, *totalBGff_JetReq,*totalBGff_2JetReq;
TH1F *eeMetMinusSideBandForErrors,*ffMet_reweightJet_binnedForErrors;
TH1F *eeMetMinusSideBandForErrors_JetReq,*ffMet_reweightJet_binnedForErrors_JetReq;
TH1F *eeMetMinusSideBandForErrors_2JetReq,*ffMet_reweightJet_binnedForErrors_2JetReq;
TH1F *eeMetNew,*ffMetNew;
void DavePlots::Loop(){
  
  bool doMC=false, doPileup=false;
  float FakeRate=0.0181,FakeRateErr=0.0003,FakeRateErrForLimits=0.;


  //TFile fin("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Data2012_Filter_JsonHLTnVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ICHEPdata_2012IDLoose_15GevFakeHighCut_sihih012_pixelVetoNoEleVeto.root","READ");
  //TFile fin("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_15GevCHfakeCut_PixelCutOnFakes_NewDiEMPtBins_19499pb_NewJetMatching_2JetReq_Apr17.root","READ");
  TFile fin("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_15GevCHfakeCut_PixelCutOnFakes_NewDiEMPtBins_19499pb_NewJetMatching_2JetReq_May7.root","READ");
  //TFile fin("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_15GevCHfakeCut_EleVetoOnFakes_19499pb.root","READ");
  //TFile fin("reweights.root","READ");
  TFile fSig("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_SignalMC_2000_1015_305_ALL_PU_Ntuplized_RhoPileupCorr_8TeV_Photon.root","READ");
  //TFile fSig2012("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_SignalMCAcceptance_1400_1420_375.root","READ");
  //TFile fSig2012_2("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_SignalMCAcceptance_1100_720_375.root","READ");
  //TFile fSig2012("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_SignalMCAcceptance_Bino_2012IDloose_15GevCHfakeCut_sihih012_pixelCutNoEleVeto_MS_1400_MG_1720_MN_375.root");
  //TFile fSig2012_2("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_SignalMCAcceptance_Bino_2012IDloose_15GevCHfakeCut_sihih012_pixelCutNoEleVeto_MS_1100_MG_720_MN_375.root");
  TFile fSig2012("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_SignalMCAcceptance_Bino_2012IDloose_15GevCHfakeCut_sihih012_pixelCutNoEleVeto_2Jetreq_Apr17_MS_1400_MG_1720_MN_375.root");
  TFile fSig2012_2("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_SignalMCAcceptance_Bino_2012IDloose_15GevCHfakeCut_sihih012_pixelCutNoEleVeto_2Jetreq_Apr17_MS_1100_MG_720_MN_375.root");
  if(doMC){  
    //TFile fRhoCorr("/data/ndpc3/c/dmorse/RA3/AnalysisOutput/hist_Data2011A_ToRun167913_Filter_RhoPileupCorr_Photon_NEW2.root","READ");
    TFile fSig720("/data/ndpc3/c/dmorse/RA3/AnalysisOutput/hist_SignalMC_720_640_375_ALL_PU_Ntuplized_NoPileupCorr_Photon.root","READ");
    TFile fSig800("/data/ndpc3/c/dmorse/RA3/AnalysisOutput/hist_SignalMC_800_960_375_ALL_PU_Ntuplized_NoPileupCorr_Photon.root","READ");
  }
  TFile fout("TRASH.root","RECREATE");
  TFile fLimits_nojet("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/met_reweighted_nojet.root","RECREATE");
  TFile fLimits_1jet("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/met_reweighted_1jet.root","RECREATE");
  TFile fLimits_2jet("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/met_reweighted_2jet.root","RECREATE");
  //TFile fLimitsSig("signal_contamination_bino_chi0375.root","RECREATE");

  //TFile f_ZZ("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/SignalMC/hist_ForDiBosonSubtraction_ZZ_cms533v1.root","READ");
  //TFile f_WZ("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/SignalMC/hist_ForDiBosonSubtraction_WZ_cms533v1.root","READ");
  TFile f_ZZ("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_cms533v1_ZZTo2L2Nu_TuneZ2star_8TeV_pythia6_tauola_Summer12-PU_S7_START52_V9-v1_May6.root");
  TFile f_WZ("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_cms533v1_WZTo3LNu_TuneZ2star_8TeV_pythia6_tauola_Summer12-PU_S7_START52_V9-v1_May6.root");
  /*
  TFile f_ggHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");
  TFile f_WZHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");
  TFile f_TTHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_TTH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");
  TFile f_VBFHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root","READ");
  */
  TFile f_ggHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_May7.root","READ");
  TFile f_WZHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_May7.root","READ");
  TFile f_TTHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_TTH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_May7.root","READ");
  TFile f_VBFHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_May7.root","READ");

  TFile f_MetAndMassShapes("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MetAndMassShapes.root","RECREATE");

  //This toggles whether or not to use ff/gf combination with ff syst background, or ff with gf systematic background (though second case sends ee to limits...need to change this)
  bool doffgfcomb=true;

  TCanvas *c4 = new TCanvas("c4","",1700,800);
  c4->Divide(2,1);

  TCanvas *c3 = new TCanvas("c3","",1350,900);
  c3->Divide(3,2);

  TCanvas *c2 = new TCanvas("c2","",900,850);
  //c2->Divide(1,2,0,.05);  
  TPad *p1 = new TPad("p1","",0.,.3,1.,1.);
  p1->SetBottomMargin(0);
  p1->Draw();
  TPad *p2 = new TPad("p2","",0.,.125,1.,.3);
  p2->SetTopMargin(0);
  p2->Draw();
  TPad *p3 = new TPad("p3","",0.,0,1.,.125);
  p3->SetTopMargin(0);
  p3->Draw();

  TCanvas *c1 = new TCanvas("c1","",900,600);
  c1->cd();
  c1->Draw();

  TPaveText *ProtoText;
  ProtoText = new TPaveText(.53,.58,.83,.79,"NDC");
  ProtoText->AddText("CMS Preliminary");
  ProtoText->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  //ProtoText->AddText("==0 Jet Requirement");
  ProtoText->SetFillStyle(4000);
  ProtoText->SetFillColor(0);
  ProtoText->SetBorderSize(0);
  TPaveText *InvarMassTextee;
  InvarMassTextee = new TPaveText(.185,.58,.485,.8,"NDC");
  InvarMassTextee->AddText("CMS Preliminary");
  InvarMassTextee->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  InvarMassTextee->AddText("e^{+}e^{-} sample");
  InvarMassTextee->SetFillStyle(4000);
  InvarMassTextee->SetFillColor(0);
  InvarMassTextee->SetBorderSize(0);
  TPaveText *InvarMassTexteg;
  InvarMassTexteg = new TPaveText(.185,.58,.485,.8,"NDC");
  InvarMassTexteg->AddText("CMS Preliminary");
  InvarMassTexteg->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  InvarMassTexteg->AddText("e#gamma sample");
  InvarMassTexteg->SetFillStyle(4000);
  InvarMassTexteg->SetFillColor(0);
  InvarMassTexteg->SetBorderSize(0);
  TPaveText *Text_noJet;
  Text_noJet = new TPaveText(.23,.58,.53,.79,"NDC");
  Text_noJet->AddText("CMS Preliminary");
  Text_noJet->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  Text_noJet->AddText("no Jet Requirement");
  Text_noJet->SetFillStyle(4000);
  Text_noJet->SetFillColor(0);
  Text_noJet->SetBorderSize(0);
  TPaveText *Text_Jet;
  Text_Jet = new TPaveText(.53,.58,.83,.79,/*.23,.58,.53,.79,*/"NDC");
  Text_Jet->AddText("CMS Preliminary");
  Text_Jet->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  Text_Jet->AddText(">=1 Jet Requirement");
  Text_Jet->SetFillStyle(4000);
  Text_Jet->SetFillColor(0);
  Text_Jet->SetBorderSize(0);
  TPaveText *Text_0Jet;
  Text_0Jet = new TPaveText(.53,.58,.83,.79,"NDC");
  Text_0Jet->AddText("CMS Preliminary");
  Text_0Jet->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  Text_0Jet->AddText("==0 Jet Requirement");
  Text_0Jet->SetFillStyle(4000);
  Text_0Jet->SetFillColor(0);
  Text_0Jet->SetBorderSize(0);
  TPaveText *Text_1Jet;
  Text_1Jet = new TPaveText(.53,.58,.83,.79,"NDC");
  Text_1Jet->AddText("CMS Preliminary");
  Text_1Jet->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  Text_1Jet->AddText("==1 Jet Requirement");
  Text_1Jet->SetFillStyle(4000);
  Text_1Jet->SetFillColor(0);
  Text_1Jet->SetBorderSize(0);
  TPaveText *Text_2Jet;
  Text_2Jet = new TPaveText(.53,.58,.83,.79,"NDC");
  Text_2Jet->AddText("CMS Preliminary");
  Text_2Jet->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  Text_2Jet->AddText(">=2 Jet Requirement");
  Text_2Jet->SetFillStyle(4000);
  Text_2Jet->SetFillColor(0);
  Text_2Jet->SetBorderSize(0);
  TPaveText *Textee = new TPaveText(.53,.5,.83,.58,"NDC");
  Textee->AddText("e^{+}e^{-} sample");
  Textee->SetFillStyle(4000);
  Textee->SetFillColor(0);
  Textee->SetBorderSize(0);
  TPaveText *Textff = new TPaveText(.53,.5,.83,.58,"NDC");
  Textff->AddText("fake fake sample");
  Textff->SetFillStyle(4000);
  Textff->SetFillColor(0);
  Textff->SetBorderSize(0);
  TPaveText *Textgammafake = new TPaveText(.53,.5,.83,.58,"NDC");
  Textgammafake->AddText("#gamma fake sample");
  Textgammafake->SetFillStyle(4000);
  Textgammafake->SetFillColor(0);
  Textgammafake->SetBorderSize(0);
  TPaveText *Text25to40 = new TPaveText(.185,.5,.485,.58,"NDC");
  Text25to40->AddText("25<p_{T}<40 GeV");
  Text25to40->SetFillStyle(4000);
  Text25to40->SetFillColor(0);
  Text25to40->SetBorderSize(0);
  TPaveText *Text40to45 = new TPaveText(.185,.5,.485,.58,"NDC");
  Text40to45->AddText("40<p_{T}<45 GeV");
  Text40to45->SetFillStyle(4000);
  Text40to45->SetFillColor(0);
  Text40to45->SetBorderSize(0);
  TPaveText *Text45to50 = new TPaveText(.185,.5,.485,.58,"NDC");
  Text45to50->AddText("45<p_{T}<50 GeV");
  Text45to50->SetFillStyle(4000);
  Text45to50->SetFillColor(0);
  Text45to50->SetBorderSize(0);
  TPaveText *Text50to60 = new TPaveText(.185,.5,.485,.58,"NDC");
  Text50to60->AddText("50<p_{T}<60 GeV");
  Text50to60->SetFillStyle(4000);
  Text50to60->SetFillColor(0);
  Text50to60->SetBorderSize(0);
  TPaveText *Text60to80 = new TPaveText(.185,.5,.485,.58,"NDC");
  Text60to80->AddText("60<p_{T}<80 GeV");
  Text60to80->SetFillStyle(4000);
  Text60to80->SetFillColor(0);
  Text60to80->SetBorderSize(0);
  TPaveText *Text80 = new TPaveText(.185,.5,.485,.58,"NDC");
  Text80->AddText("p_{T}>80 GeV");
  Text80->SetFillStyle(4000);
  Text80->SetFillColor(0);
  Text80->SetBorderSize(0);

  TLine *line125 = new TLine(125.3,0,125.3,740);
  //line125->SetLineColor(kCyan);
  gStyle->SetLineStyleString(11,"6 20");

  TH1F* ggInvarMassMET = (TH1F*)fin.Get("ggInvarMassMVAcorrVertexCorr");
  TH1F* ggInvarMassMET_JetReq = (TH1F*)fin.Get("ggInvarMass_JetReq");//fix to have mvacorrvertexcorr for jetreq
  TH1F* ggInvarMass = (TH1F*)ggInvarMassMET->Clone();
  c1->SetLogy(0); gStyle->SetOptStat(0);
  TH1F* ggInvarMassMET30 = (TH1F*)fin.Get("ggInvarMassMET30MVAcorr");
  TH1F* ggInvarMassMET30Clone = (TH1F*)ggInvarMassMET30->Clone();
  TH1F* ggInvarMassMET30_JetReq = (TH1F*)fin.Get("ggInvarMassMET30_JetReq");
  TH1F* egInvarMassMET30 = (TH1F*)fin.Get("egInvarMassMET30MVAcorr");
  TH1F* ffInvarMassMET30 = (TH1F*)fin.Get("ffInvarMassMET30MVAcorr");
  TH1F* gammafakeInvarMassMET30 = (TH1F*)fin.Get("gammafakeInvarMassMET30MVAcorr");
  TH1F* ggInvarMassMET40 = (TH1F*)fin.Get("ggInvarMassMET40MVAcorr");
  TH1F* ggInvarMassMET40_JetReq = (TH1F*)fin.Get("ggInvarMassMET40_JetReq");
  TH1F* ggInvarMassMET50 = (TH1F*)fin.Get("ggInvarMassMET50MVAcorr");
  TH1F* ggInvarMassMET50_JetReq = (TH1F*)fin.Get("ggInvarMassMET50_JetReq");
  TH1F* ggInvarMassMET60 = (TH1F*)fin.Get("ggInvarMassMET60MVAcorr");
  TH1F* ggInvarMassMET60_JetReq = (TH1F*)fin.Get("ggInvarMassMET60_JetReq");
  TH1F* ggInvarMassMET70 = (TH1F*)fin.Get("ggInvarMassMET70MVAcorr");
  TH1F* ggInvarMassMET70_JetReq = (TH1F*)fin.Get("ggInvarMassMET70_JetReq");
  TH1F* ggInvarMassMET80 = (TH1F*)fin.Get("ggInvarMassMET80MVAcorr");
  TH1F* ggInvarMassMET100 = (TH1F*)fin.Get("ggInvarMassMET100MVAcorr");

  //Higgs dataset:
  //ggHgg:96290 events
  //WZHgg:100320
  //VBFHgg:99885
  //TTHgg:100224 but jobs fail so use number ntuplized (94424)

  TH1F* ggHggNoMET = (TH1F*)f_ggHgg.Get("ggInvarMassMVAcorrVertexCorr");
  TH1F* WZHggNoMET = (TH1F*)f_WZHgg.Get("ggInvarMassMVAcorrVertexCorr");
  TH1F* TTHggNoMET = (TH1F*)f_TTHgg.Get("ggInvarMassMVAcorrVertexCorr");
  TH1F* VBFHggNoMET = (TH1F*)f_VBFHgg.Get("ggInvarMassMVAcorrVertexCorr");

  //fix these to have mvacorrvertexcorr in jetreq - fixed, now all invarmass are mvacorrvertexcorr
  TH1F* ggHggNoMET_JetReq = (TH1F*)f_ggHgg.Get("ggInvarMass_JetReq");
  TH1F* WZHggNoMET_JetReq = (TH1F*)f_WZHgg.Get("ggInvarMass_JetReq");
  TH1F* TTHggNoMET_JetReq = (TH1F*)f_TTHgg.Get("ggInvarMass_JetReq");
  TH1F* VBFHggNoMET_JetReq = (TH1F*)f_VBFHgg.Get("ggInvarMass_JetReq");

  TH1F* ggHggMET30 = (TH1F*)f_ggHgg.Get("ggInvarMassMET30MVAcorr");
  TH1F* WZHggMET30 = (TH1F*)f_WZHgg.Get("ggInvarMassMET30MVAcorr");
  TH1F* TTHggMET30 = (TH1F*)f_TTHgg.Get("ggInvarMassMET30MVAcorr");
  TH1F* VBFHggMET30 = (TH1F*)f_VBFHgg.Get("ggInvarMassMET30MVAcorr");
  
  TH1F* ggHggMET30_JetReq = (TH1F*)f_ggHgg.Get("ggInvarMassMET30_JetReq");
  TH1F* WZHggMET30_JetReq = (TH1F*)f_WZHgg.Get("ggInvarMassMET30_JetReq");
  TH1F* TTHggMET30_JetReq = (TH1F*)f_TTHgg.Get("ggInvarMassMET30_JetReq");
  TH1F* VBFHggMET30_JetReq = (TH1F*)f_VBFHgg.Get("ggInvarMassMET30_JetReq");
  
  TH1F* ggHggMET40 = (TH1F*)f_ggHgg.Get("ggInvarMassMET40MVAcorr");
  TH1F* WZHggMET40 = (TH1F*)f_WZHgg.Get("ggInvarMassMET40MVAcorr");
  TH1F* TTHggMET40 = (TH1F*)f_TTHgg.Get("ggInvarMassMET40MVAcorr");
  TH1F* VBFHggMET40 = (TH1F*)f_VBFHgg.Get("ggInvarMassMET40MVAcorr");

  TH1F* ggHggMET50 = (TH1F*)f_ggHgg.Get("ggInvarMassMET50MVAcorr");
  TH1F* WZHggMET50 = (TH1F*)f_WZHgg.Get("ggInvarMassMET50MVAcorr");
  TH1F* TTHggMET50 = (TH1F*)f_TTHgg.Get("ggInvarMassMET50MVAcorr");
  TH1F* VBFHggMET50 = (TH1F*)f_VBFHgg.Get("ggInvarMassMET50MVAcorr");

  TH1F* ggHggMET60 = (TH1F*)f_ggHgg.Get("ggInvarMassMET60MVAcorr");
  TH1F* WZHggMET60 = (TH1F*)f_WZHgg.Get("ggInvarMassMET60MVAcorr");
  TH1F* TTHggMET60 = (TH1F*)f_TTHgg.Get("ggInvarMassMET60MVAcorr");
  TH1F* VBFHggMET60 = (TH1F*)f_VBFHgg.Get("ggInvarMassMET60MVAcorr");

  TH1F* ggHggMET70 = (TH1F*)f_ggHgg.Get("ggInvarMassMET70MVAcorr");
  TH1F* WZHggMET70 = (TH1F*)f_WZHgg.Get("ggInvarMassMET70MVAcorr");
  TH1F* TTHggMET70 = (TH1F*)f_TTHgg.Get("ggInvarMassMET70MVAcorr");
  TH1F* VBFHggMET70 = (TH1F*)f_VBFHgg.Get("ggInvarMassMET70MVAcorr");

  TH1F* ggHggMET80 = (TH1F*)f_ggHgg.Get("ggInvarMassMET80MVAcorr");
  TH1F* WZHggMET80 = (TH1F*)f_WZHgg.Get("ggInvarMassMET80MVAcorr");
  TH1F* TTHggMET80 = (TH1F*)f_TTHgg.Get("ggInvarMassMET80MVAcorr");
  TH1F* VBFHggMET80 = (TH1F*)f_VBFHgg.Get("ggInvarMassMET80MVAcorr");

  TH1F* ggHggMET100 = (TH1F*)f_ggHgg.Get("ggInvarMassMET100MVAcorr");
  TH1F* WZHggMET100 = (TH1F*)f_WZHgg.Get("ggInvarMassMET100MVAcorr");
  TH1F* TTHggMET100 = (TH1F*)f_TTHgg.Get("ggInvarMassMET100MVAcorr");
  TH1F* VBFHggMET100 = (TH1F*)f_VBFHgg.Get("ggInvarMassMET100MVAcorr");


  TH1F* ggHggMet = (TH1F*)f_ggHgg.Get("ggMet");
  TH1F* WZHggMet = (TH1F*)f_WZHgg.Get("ggMet");
  TH1F* TTHggMet = (TH1F*)f_TTHgg.Get("ggMet");
  TH1F* VBFHggMet = (TH1F*)f_VBFHgg.Get("ggMet");

  TH1F* ggHggInvarMassMET30 = (TH1F*)ggHggMET30->Clone();//(TH1F*)f_ggHgg.Get("ggInvarMassMET30MVAcorr");
  TH1F* WZHggInvarMassMET30 = (TH1F*)WZHggMET30->Clone();//(TH1F*)f_WZHgg.Get("ggInvarMassMET30MVAcorr");
  TH1F* TTHggInvarMassMET30 = (TH1F*)TTHggMET30->Clone();//(TH1F*)f_TTHgg.Get("ggInvarMassMET30MVAcorr");
  TH1F* VBFHggInvarMassMET30 = (TH1F*)VBFHggMET30->Clone();//(TH1F*)f_VBFHgg.Get("ggInvarMassMET30MVAcorr");

  TH1F* ggHggInvarMassMET40 = (TH1F*)ggHggMET40->Clone();//(TH1F*)f_ggHgg.Get("ggInvarMassMET40MVAcorr");
  TH1F* WZHggInvarMassMET40 = (TH1F*)WZHggMET40->Clone();//(TH1F*)f_WZHgg.Get("ggInvarMassMET40MVAcorr");
  TH1F* TTHggInvarMassMET40 = (TH1F*)TTHggMET40->Clone();//(TH1F*)f_TTHgg.Get("ggInvarMassMET40MVAcorr");
  TH1F* VBFHggInvarMassMET40 = (TH1F*)VBFHggMET40->Clone();//(TH1F*)f_VBFHgg.Get("ggInvarMassMET40MVAcorr");

  ggHggMet->Sumw2();WZHggMet->Sumw2();TTHggMet->Sumw2();VBFHggMet->Sumw2();
  ggHggInvarMassMET30->Sumw2();WZHggInvarMassMET30->Sumw2();TTHggInvarMassMET30->Sumw2();VBFHggInvarMassMET30->Sumw2();
  ggHggInvarMassMET40->Sumw2();WZHggInvarMassMET40->Sumw2();TTHggInvarMassMET40->Sumw2();VBFHggInvarMassMET40->Sumw2();
  ggHggNoMET->Sumw2();WZHggNoMET->Sumw2();TTHggNoMET->Sumw2();VBFHggNoMET->Sumw2();
  ggHggNoMET_JetReq->Sumw2();WZHggNoMET_JetReq->Sumw2();TTHggNoMET_JetReq->Sumw2();VBFHggNoMET_JetReq->Sumw2();
  ggHggMET30->Sumw2();WZHggMET30->Sumw2();TTHggMET30->Sumw2();VBFHggMET30->Sumw2();
  ggHggMET40->Sumw2();WZHggMET40->Sumw2();TTHggMET40->Sumw2();VBFHggMET40->Sumw2();
  ggHggMET50->Sumw2();WZHggMET50->Sumw2();TTHggMET50->Sumw2();VBFHggMET50->Sumw2();
  ggHggMET60->Sumw2();WZHggMET60->Sumw2();TTHggMET60->Sumw2();VBFHggMET60->Sumw2();
  ggHggMET70->Sumw2();WZHggMET70->Sumw2();TTHggMET70->Sumw2();VBFHggMET70->Sumw2();
  ggHggMET80->Sumw2();WZHggMET80->Sumw2();TTHggMET80->Sumw2();VBFHggMET80->Sumw2();
  ggHggMET100->Sumw2();WZHggMET100->Sumw2();TTHggMET100->Sumw2();VBFHggMET100->Sumw2();

  TH1F* sigStrengthVsMet = new TH1F("sigStrengthVsMet","Signal Yield from fit / SM yield",11,0,110);
  TH1F* sigStrengthVsMet_JetReq = new TH1F("sigStrengthVsMet_JetReq","Signal Yield from fit / SM yield with >=1 jet requirement",11,0,110);


  //from DAS: ggHgg: 96290, VBFHgg:100118 , WZHgg:100320 , TTHgg:100224
  //from prep: ggHgg:100000 , VBFHgg:100000 , WZHgg:100000 , TTHgg:100000
  //used before May14: ggHgg: 96290, VBFHgg:99885 , WZHgg:100320 , TTHgg:94424
  //x-sections: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt8TeV2012ICHEP
  //L_int*(H->gg BF)*(individual xSection (pb))/nEvents
  ggHggMet->Scale((L_int*2.29e-03*19.52)/100000);//Scale(1./ggHggMet->Integral());
  WZHggMet->Scale((L_int*2.29e-03*(.6966/**(.3257+.014)*/+.3943/**.2*/))/100000);//Scale(1./WZHggMet->Integral());
  TTHggMet->Scale((L_int*2.29e-03*.1302)/100000);//Scale(1./TTHggMet->Integral());
  VBFHggMet->Scale((L_int*2.29e-03*1.559)/100000);//Scale(1./VBFHggMet->Integral());

  ggHggMet->SetLineColor(kGreen);ggHggMet->SetMarkerColor(kGreen);
  WZHggMet->SetLineColor(kCyan);WZHggMet->SetMarkerColor(kCyan);
  TTHggMet->SetLineColor(kViolet);TTHggMet->SetMarkerColor(kViolet);
  VBFHggMet->SetLineColor(kRed+3);VBFHggMet->SetMarkerColor(kRed+3);
  ggHggMet->SetLineWidth(2);WZHggMet->SetLineWidth(2);VBFHggMet->SetLineWidth(2);TTHggMet->SetLineWidth(2);

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

  ggHggMET30->Scale((L_int*2.29e-03*19.52)/100000);  
  VBFHggMET30->Scale((L_int*2.29e-03*1.559)/100000); 
  TTHggMET30->Scale((L_int*2.29e-03*.1302)/100000);
  WZHggMET30->Scale((L_int*2.29e-03*(.6966/**(.3257+.014)*/+.3943/**.2*/))/100000);
  ggHggMET30->SetLineColor(kGreen);ggHggMET30->SetMarkerColor(kGreen);
  WZHggMET30->SetLineColor(kCyan);WZHggMET30->SetMarkerColor(kCyan);
  TTHggMET30->SetLineColor(kViolet);TTHggMET30->SetMarkerColor(kViolet);
  VBFHggMET30->SetLineColor(kRed+3);VBFHggMET30->SetMarkerColor(kRed+3);
  ggHggMET30->SetLineWidth(2);WZHggMET30->SetLineWidth(2);VBFHggMET30->SetLineWidth(2);TTHggMET30->SetLineWidth(2);
  
  ggHggMET30_JetReq->Scale((L_int*2.29e-03*19.52)/100000);  
  VBFHggMET30_JetReq->Scale((L_int*2.29e-03*1.559)/100000); 
  TTHggMET30_JetReq->Scale((L_int*2.29e-03*.1302)/100000);
  WZHggMET30_JetReq->Scale((L_int*2.29e-03*(.6966+.3943))/100000);
  ggHggMET30_JetReq->SetLineColor(kGreen);ggHggMET30_JetReq->SetMarkerColor(kGreen);
  WZHggMET30_JetReq->SetLineColor(kCyan);WZHggMET30_JetReq->SetMarkerColor(kCyan);
  TTHggMET30_JetReq->SetLineColor(kViolet);TTHggMET30_JetReq->SetMarkerColor(kViolet);
  VBFHggMET30_JetReq->SetLineColor(kRed+3);VBFHggMET30_JetReq->SetMarkerColor(kRed+3);
  ggHggMET30_JetReq->SetLineWidth(2);WZHggMET30_JetReq->SetLineWidth(2);VBFHggMET30_JetReq->SetLineWidth(2);TTHggMET30_JetReq->SetLineWidth(2);

  ggHggMET40->Scale((L_int*2.29e-03*19.52)/100000);  
  VBFHggMET40->Scale((L_int*2.29e-03*1.559)/100000); 
  TTHggMET40->Scale((L_int*2.29e-03*.1302)/100000);
  WZHggMET40->Scale((L_int*2.29e-03*(.6966/**(.3257+.014)*/+.3943/**.2*/))/100000);
  ggHggMET40->SetLineColor(kGreen);ggHggMET40->SetMarkerColor(kGreen);
  WZHggMET40->SetLineColor(kCyan);WZHggMET40->SetMarkerColor(kCyan);
  TTHggMET40->SetLineColor(kViolet);TTHggMET40->SetMarkerColor(kViolet);
  VBFHggMET40->SetLineColor(kRed+3);VBFHggMET40->SetMarkerColor(kRed+3);
  ggHggMET40->SetFillStyle(0);WZHggMET40->SetFillStyle(0);TTHggMET40->SetFillStyle(0);VBFHggMET40->SetFillStyle(0);
  ggHggMET40->SetLineWidth(2);WZHggMET40->SetLineWidth(2);VBFHggMET40->SetLineWidth(2);TTHggMET40->SetLineWidth(2);

  ggHggMET50->Scale((L_int*2.29e-03*19.52)/100000);  
  VBFHggMET50->Scale((L_int*2.29e-03*1.559)/100000); 
  TTHggMET50->Scale((L_int*2.29e-03*.1302)/100000);
  WZHggMET50->Scale((L_int*2.29e-03*(.6966/**(.3257+.014)*/+.3943/**.2*/))/100000);
  ggHggMET50->SetLineColor(kGreen);ggHggMET50->SetMarkerColor(kGreen);
  WZHggMET50->SetLineColor(kCyan);WZHggMET50->SetMarkerColor(kCyan);
  TTHggMET50->SetLineColor(kViolet);TTHggMET50->SetMarkerColor(kViolet);
  VBFHggMET50->SetLineColor(kRed+3);VBFHggMET50->SetMarkerColor(kRed+3);
  ggHggMET50->SetLineWidth(2);WZHggMET50->SetLineWidth(2);VBFHggMET50->SetLineWidth(2);TTHggMET50->SetLineWidth(2);

  ggHggMET60->Scale((L_int*2.29e-03*19.52)/100000);  
  VBFHggMET60->Scale((L_int*2.29e-03*1.559)/100000); 
  TTHggMET60->Scale((L_int*2.29e-03*.1302)/100000);
  WZHggMET60->Scale((L_int*2.29e-03*(.6966/**(.3257+.014)*/+.3943/**.2*/))/100000);
  ggHggMET60->SetLineColor(kGreen);ggHggMET60->SetMarkerColor(kGreen);
  WZHggMET60->SetLineColor(kCyan);WZHggMET60->SetMarkerColor(kCyan);
  TTHggMET60->SetLineColor(kViolet);TTHggMET60->SetMarkerColor(kViolet);
  VBFHggMET60->SetLineColor(kRed+3);VBFHggMET60->SetMarkerColor(kRed+3);
  ggHggMET60->SetLineWidth(2);WZHggMET60->SetLineWidth(2);VBFHggMET60->SetLineWidth(2);TTHggMET60->SetLineWidth(2);

  ggHggMET70->Scale((L_int*2.29e-03*19.52)/100000);  
  VBFHggMET70->Scale((L_int*2.29e-03*1.559)/100000); 
  TTHggMET70->Scale((L_int*2.29e-03*.1302)/100000);
  WZHggMET70->Scale((L_int*2.29e-03*(.6966/**(.3257+.014)*/+.3943/**.2*/))/100000);
  ggHggMET70->SetLineColor(kGreen);ggHggMET70->SetMarkerColor(kGreen);
  WZHggMET70->SetLineColor(kCyan);WZHggMET70->SetMarkerColor(kCyan);
  TTHggMET70->SetLineColor(kViolet);TTHggMET70->SetMarkerColor(kViolet);
  VBFHggMET70->SetLineColor(kRed+3);VBFHggMET70->SetMarkerColor(kRed+3);
  ggHggMET70->SetLineWidth(2);WZHggMET70->SetLineWidth(2);VBFHggMET70->SetLineWidth(2);TTHggMET70->SetLineWidth(2);

  ggHggMET80->Scale((L_int*2.29e-03*19.52)/100000);  
  VBFHggMET80->Scale((L_int*2.29e-03*1.559)/100000); 
  TTHggMET80->Scale((L_int*2.29e-03*.1302)/100000);
  WZHggMET80->Scale((L_int*2.29e-03*(.6966/**(.3257+.014)*/+.3943/**.2*/))/100000);
  ggHggMET80->SetLineColor(kGreen);ggHggMET80->SetMarkerColor(kGreen);
  WZHggMET80->SetLineColor(kCyan);WZHggMET80->SetMarkerColor(kCyan);
  TTHggMET80->SetLineColor(kViolet);TTHggMET80->SetMarkerColor(kViolet);
  VBFHggMET80->SetLineColor(kRed+3);VBFHggMET80->SetMarkerColor(kRed+3);
  ggHggMET80->SetLineWidth(2);WZHggMET80->SetLineWidth(2);VBFHggMET80->SetLineWidth(2);TTHggMET80->SetLineWidth(2);

  ggHggMET100->Scale((L_int*2.29e-03*19.52)/100000);  
  VBFHggMET100->Scale((L_int*2.29e-03*1.559)/100000); 
  TTHggMET100->Scale((L_int*2.29e-03*.1302)/100000);
  WZHggMET100->Scale((L_int*2.29e-03*(.6966/**(.3257+.014)*/+.3943/**.2*/))/100000);
  ggHggMET100->SetLineColor(kGreen);ggHggMET100->SetMarkerColor(kGreen);
  WZHggMET100->SetLineColor(kCyan);WZHggMET100->SetMarkerColor(kCyan);
  TTHggMET100->SetLineColor(kViolet);TTHggMET100->SetMarkerColor(kViolet);
  VBFHggMET100->SetLineColor(kRed+3);VBFHggMET100->SetMarkerColor(kRed+3);
  ggHggMET100->SetLineWidth(2);WZHggMET100->SetLineWidth(2);VBFHggMET100->SetLineWidth(2);TTHggMET100->SetLineWidth(2);
  
  //ggHggMet->Rebin(2);WZHggMet->Rebin(2);VBFHggMet->Rebin(2);TTHggMet->Rebin(2);
  //ggHggNoMET->Rebin(2);WZHggNoMET->Rebin(2);VBFHggNoMET->Rebin(2);TTHggNoMET->Rebin(2);
  //ggHggNoMET_JetReq->Rebin(2);WZHggNoMET_JetReq->Rebin(2);VBFHggNoMET_JetReq->Rebin(2);TTHggNoMET_JetReq->Rebin(2);
  //ggHggMET30->Rebin(2);WZHggMET30->Rebin(2);VBFHggMET30->Rebin(2);TTHggMET30->Rebin(2);
  //ggHggMET30_JetReq->Rebin(2);WZHggMET30_JetReq->Rebin(2);VBFHggMET30_JetReq->Rebin(2);TTHggMET30_JetReq->Rebin(2);
  //ggHggMET40->Rebin(2);WZHggMET40->Rebin(2);VBFHggMET40->Rebin(2);TTHggMET40->Rebin(2);
  //ggHggMET50->Rebin(2);WZHggMET50->Rebin(2);VBFHggMET50->Rebin(2);TTHggMET50->Rebin(2);
  //ggHggMET60->Rebin(2);WZHggMET60->Rebin(2);VBFHggMET60->Rebin(2);TTHggMET60->Rebin(2);
  //ggHggMET70->Rebin(2);WZHggMET70->Rebin(2);VBFHggMET70->Rebin(2);TTHggMET70->Rebin(2);
  //ggHggMET80->Rebin(2);WZHggMET80->Rebin(2);VBFHggMET80->Rebin(2);TTHggMET80->Rebin(2);
  //ggHggMET100->Rebin(2);WZHggMET100->Rebin(2);VBFHggMET100->Rebin(2);TTHggMET100->Rebin(2);

  ggHggInvarMassMET30->Scale((L_int*2.29e-03*19.52)/100000);  
  VBFHggInvarMassMET30->Scale((L_int*2.29e-03*1.559)/100000); 
  TTHggInvarMassMET30->Scale((L_int*2.29e-03*.1302)/100000);
  WZHggInvarMassMET30->Scale((L_int*2.29e-03*(.6966/**(.3257+.014)*/+.3943/**.2*/))/100000);

  ggHggInvarMassMET40->Scale((L_int*2.29e-03*19.52)/100000);  
  VBFHggInvarMassMET40->Scale((L_int*2.29e-03*1.559)/100000); 
  TTHggInvarMassMET40->Scale((L_int*2.29e-03*.1302)/100000);
  WZHggInvarMassMET40->Scale((L_int*2.29e-03*(.6966/**(.3257+.014)*/+.3943/**.2*/))/100000);

  TH1F* HiggsMetForShape = (TH1F*)TTHggMet->Clone();HiggsMetForShape->Add(VBFHggMet);HiggsMetForShape->Add(ggHggMet);HiggsMetForShape->Add(WZHggMet);
  TH1F* HiggsMetForShape30 = (TH1F*)TTHggMET30->Clone();HiggsMetForShape30->Add(VBFHggMET30);HiggsMetForShape30->Add(ggHggMET30);HiggsMetForShape30->Add(WZHggMET30);
  TH1F* HiggsMetForShape30_2 = (TH1F*)TTHggInvarMassMET30->Clone();HiggsMetForShape30_2->Add(VBFHggInvarMassMET30);HiggsMetForShape30_2->Add(ggHggInvarMassMET30);HiggsMetForShape30_2->Add(WZHggInvarMassMET30);
  float SMhiggsYield=HiggsMetForShape->Integral(0,999);
  cout<<"Standard model Higgs expected yield: "<<HiggsMetForShape->Integral()<<endl;
  cout<<"Standard model Higgs expected yield, MET>30 from full hist: "<<HiggsMetForShape->Integral(7,999)<<endl;
  cout<<"Standard model Higgs expected yield, MET>30: "<<HiggsMetForShape30->Integral()<<endl;
  cout<<"Standard model Higgs expected yield, MET>30 from clone: "<<HiggsMetForShape30_2->Integral()<<endl;
  TH1F* HiggsInvarMassForShape = (TH1F*)TTHggInvarMassMET40->Clone();HiggsInvarMassForShape->Add(VBFHggInvarMassMET40);HiggsInvarMassForShape->Add(ggHggInvarMassMET40);HiggsInvarMassForShape->Add(WZHggInvarMassMET40);
  HiggsMetForShape->Rebin(4);
  float HiggsScale = 1./HiggsMetForShape->Integral(0,999);
  HiggsMetForShape->Scale(HiggsScale);
  HiggsScale = 1./HiggsInvarMassForShape->Integral(0,999);
  HiggsInvarMassForShape->Scale(HiggsScale);
  f_MetAndMassShapes.cd();
  HiggsMetForShape->Write("StandardModelHiggsMet");
  HiggsInvarMassForShape->Write("StandardModelHiggsInvarMass");
  fin.cd();
  TH1F* TTHggMetNew = (TH1F*)TTHggMet->Rebin(NmetBins,"TTHggMetNew",xbins);
  TH1F* VBFHggMetNew = (TH1F*)VBFHggMet->Rebin(NmetBins,"VBFHggMetNew",xbins);
  TH1F* ggHggMetNew = (TH1F*)ggHggMet->Rebin(NmetBins,"ggHggMetNew",xbins);
  TH1F* WZHggMetNew = (TH1F*)WZHggMet->Rebin(NmetBins,"WZHggMetNew",xbins);
  TTHggMetNew->SetFillColor(0);VBFHggMetNew->SetFillColor(0);ggHggMetNew->SetFillColor(0);WZHggMetNew->SetFillColor(0);
  for(int i=0;i<=TTHggMetNew->GetNbinsX();i++){
    float x =  TTHggMetNew->GetBinContent(i)/TTHggMetNew->GetBinWidth(i);
    TTHggMetNew->SetBinContent(i,x);
    x = TTHggMetNew->GetBinError(i)/TTHggMetNew->GetBinWidth(i);
    TTHggMetNew->SetBinError(i,x);
    x =  VBFHggMetNew->GetBinContent(i)/VBFHggMetNew->GetBinWidth(i);
    VBFHggMetNew->SetBinContent(i,x);
    x = VBFHggMetNew->GetBinError(i)/VBFHggMetNew->GetBinWidth(i);
    VBFHggMetNew->SetBinError(i,x);
    x =  ggHggMetNew->GetBinContent(i)/ggHggMetNew->GetBinWidth(i);
    ggHggMetNew->SetBinContent(i,x);
    x = ggHggMetNew->GetBinError(i)/ggHggMetNew->GetBinWidth(i);
    ggHggMetNew->SetBinError(i,x);
    x =  WZHggMetNew->GetBinContent(i)/WZHggMetNew->GetBinWidth(i);
    WZHggMetNew->SetBinContent(i,x);
    x = WZHggMetNew->GetBinError(i)/WZHggMetNew->GetBinWidth(i);
    WZHggMetNew->SetBinError(i,x);
  }
  cout<<"check here"<<endl;
  THStack *metstackRebin = new THStack("metstackRebin","");
  metstackRebin->Add(TTHggMetNew);
  metstackRebin->Add(VBFHggMetNew);
  metstackRebin->Add(WZHggMetNew);
  metstackRebin->Add(ggHggMetNew);
  TH1F* TTHggMetNew2 = (TH1F*)TTHggMetNew->Clone();
  TTHggMetNew2->Add(VBFHggMetNew);TTHggMetNew2->Add(ggHggMetNew);TTHggMetNew2->Add(WZHggMetNew);TTHggMetNew2->SetLineColor(kGreen);TTHggMetNew2->SetMarkerColor(kGreen);TTHggMetNew2->SetFillColor(kGreen);

  THStack *metstack = new THStack("metstack","");
  metstack->Add(TTHggMet);
  metstack->Add(VBFHggMet);
  metstack->Add(WZHggMet);
  metstack->Add(ggHggMet);
  metstack->Draw("histo");
  metstack->GetHistogram()->GetXaxis()->SetRangeUser(0,349);
  metstack->SetMinimum(1e-3);
  metstack->SetMaximum(150);
  metstack->Draw("histo");
  TLegend *legHiggs = new TLegend(.36,.65,.59,.8);
  legHiggs->AddEntry(ggHggMet,"GluGlu->H->#gamma#gamma","l");
  legHiggs->AddEntry(WZHggMet,"W/ZH->#gamma#gamma","l");
  legHiggs->AddEntry(VBFHggMet,"VBFH->#gamma#gamma","l");
  legHiggs->AddEntry(TTHggMet,"TTH->#gamma#gamma","l");
  legHiggs->SetFillColor(kWhite);
  legHiggs->SetFillStyle(0);
  legHiggs->SetBorderSize(0);
  legHiggs->Draw();
  c1->SetLogy(1);
  line125->SetLineColor(kRed);
  line125->DrawLine(30.,0.,30.,300);
  cout<<"ggHggMet integral:"<<ggHggMet->Integral()<<endl;
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggMetFromHiggsWithNoMETcutStack.png");
  c1->SetLogy(0);

  THStack *met30stack = new THStack("met30stack","");
  met30stack->Add(TTHggMET30);
  met30stack->Add(VBFHggMET30);
  met30stack->Add(WZHggMET30);
  met30stack->Add(ggHggMET30);
  met30stack->Draw("histo");
  met30stack->GetXaxis()->SetRangeUser(100,180);
  TLegend *legHiggs2 = new TLegend(.46,.65,.69,.8);
  legHiggs2->AddEntry(ggHggMET30,"GluGlu->H->#gamma#gamma","l");
  legHiggs2->AddEntry(WZHggMET30,"W/ZH->#gamma#gamma","l");
  legHiggs2->AddEntry(VBFHggMET30,"VBFH->#gamma#gamma","l");
  legHiggs2->AddEntry(TTHggMET30,"TTH->#gamma#gamma","l");
  legHiggs2->SetFillColor(kWhite);
  legHiggs2->SetFillStyle(0);
  legHiggs2->SetBorderSize(0);
  legHiggs2->Draw();
  cout<<"ggHggMET30 integral:"<<ggHggMET30->Integral()<<endl;
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassFromHiggsWithMET30cut.png");
  c1->SetLogy(0);

  //c1->SetLogy(1);
  THStack *metNoCutstack = new THStack("metNoCutstack","");
  metNoCutstack->Add(TTHggNoMET);
  metNoCutstack->Add(WZHggNoMET);
  metNoCutstack->Add(VBFHggNoMET);
  metNoCutstack->Add(ggHggNoMET);
  metNoCutstack->Draw("histo");
  metNoCutstack->GetXaxis()->SetRangeUser(100,180);
  legHiggs2->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassFromHiggsWithNoMETcutStack.png");

  c1->SetLogy(1);
  ggHggMet->SetTitle("ggMet");
  ggHggMet->GetYaxis()->SetTitle("Events");
  ggHggMet->GetXaxis()->SetRangeUser(0,520);
  ggHggMet->Draw();
  WZHggMet->Draw("SAME");
  VBFHggMet->Draw("SAME");
  TTHggMet->Draw("SAME");
  TLegend *legHiggsMet = new TLegend(.45,.6,.8,.8);
  legHiggsMet->AddEntry(TTHggMet,"TTH->#gamma#gamma","l");
  legHiggsMet->AddEntry(WZHggMet,"W/ZH->#gamma#gamma","l");
  legHiggsMet->AddEntry(VBFHggMet,"VBFH->#gamma#gamma","l");
  legHiggsMet->AddEntry(ggHggMet,"GluGlu->H->#gamma#gamma","l");
  legHiggsMet->SetFillColor(kWhite);
  legHiggsMet->SetBorderSize(0);
  legHiggsMet->Draw();
  line125->SetLineColor(kRed);
  line125->DrawLine(30.,0.,30.,1050);
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggMetFromHiggsWithNoMETcutIndividual.png");
  c1->SetLogy(0);
  line125->SetLineColor(kBlack);
  line125->SetLineStyle(11);
  line125->SetLineWidth(11);

  ggInvarMassMET30->SetTitle("");ggInvarMassMET30_JetReq->SetTitle("");ggInvarMassMET->SetTitle("");ggInvarMassMET_JetReq->SetTitle("");egInvarMassMET30->SetTitle("");ffInvarMassMET30->SetTitle("");gammafakeInvarMassMET30->SetTitle("");ggInvarMassMET40->SetTitle("");ggInvarMassMET50_JetReq->SetTitle("");ggInvarMassMET50->SetTitle("");ggInvarMassMET50_JetReq->SetTitle("");ggInvarMassMET60->SetTitle("");ggInvarMassMET60_JetReq->SetTitle("");ggInvarMassMET70->SetTitle("");ggInvarMassMET70_JetReq->SetTitle("");ggInvarMassMET80->SetTitle("");ggInvarMassMET100->SetTitle("");
  ggInvarMassMET->GetXaxis()->SetTitle("Invariant Mass (GeV)");ggInvarMassMET_JetReq->GetXaxis()->SetTitle("Invariant Mass (GeV)");ggInvarMassMET30->GetXaxis()->SetTitle("Invariant Mass (GeV)");ggInvarMassMET30_JetReq->GetXaxis()->SetTitle("Invariant Mass (GeV)");egInvarMassMET30->GetXaxis()->SetTitle("Invariant Mass (GeV)");ffInvarMassMET30->GetXaxis()->SetTitle("Invariant Mass (GeV)");ggInvarMassMET40->GetXaxis()->SetTitle("Invariant Mass (GeV)");ggInvarMassMET40_JetReq->GetXaxis()->SetTitle("Invariant Mass (GeV)");ggInvarMassMET50->GetXaxis()->SetTitle("Invariant Mass (GeV)");ggInvarMassMET50_JetReq->GetXaxis()->SetTitle("Invariant Mass (GeV)");ggInvarMassMET60->GetXaxis()->SetTitle("Invariant Mass (GeV)");ggInvarMassMET60_JetReq->GetXaxis()->SetTitle("Invariant Mass (GeV)");ggInvarMassMET70->GetXaxis()->SetTitle("Invariant Mass (GeV)");ggInvarMassMET70_JetReq->GetXaxis()->SetTitle("Invariant Mass (GeV)");ggInvarMassMET80->GetXaxis()->SetTitle("Invariant Mass (GeV)");ggInvarMassMET100->GetXaxis()->SetTitle("Invariant Mass (GeV)");
  ggInvarMassMET->Sumw2();ggInvarMassMET_JetReq->Sumw2();ggInvarMassMET30->Sumw2();ggInvarMassMET30_JetReq->Sumw2();egInvarMassMET30->Sumw2();ffInvarMassMET30->Sumw2();ggInvarMassMET40->Sumw2();ggInvarMassMET40_JetReq->Sumw2();ggInvarMassMET50->Sumw2();ggInvarMassMET50_JetReq->Sumw2();ggInvarMassMET60->Sumw2();ggInvarMassMET60_JetReq->Sumw2();ggInvarMassMET70->Sumw2();ggInvarMassMET70_JetReq->Sumw2();ggInvarMassMET80->Sumw2();ggInvarMassMET100->Sumw2();
  ggInvarMassMET->SetMarkerSize(0.5);ggInvarMassMET_JetReq->SetMarkerSize(0.5);ggInvarMassMET30->SetMarkerSize(0.5);ggInvarMassMET30_JetReq->SetMarkerSize(0.5);egInvarMassMET30->SetMarkerSize(0.5);ffInvarMassMET30->SetMarkerSize(0.5);ggInvarMassMET40->SetMarkerSize(0.5);ggInvarMassMET40_JetReq->SetMarkerSize(0.5);ggInvarMassMET50->SetMarkerSize(0.5);ggInvarMassMET50_JetReq->SetMarkerSize(0.5);ggInvarMassMET60->SetMarkerSize(0.5);ggInvarMassMET60_JetReq->SetMarkerSize(0.5);ggInvarMassMET70->SetMarkerSize(0.5);ggInvarMassMET70_JetReq->SetMarkerSize(0.5);ggInvarMassMET80->SetMarkerSize(0.5);ggInvarMassMET100->SetMarkerSize(0.5);

  ggInvarMassMET->Rebin(3);ggInvarMassMET_JetReq->Rebin(3);ggInvarMassMET30->Rebin(3);ggInvarMassMET30_JetReq->Rebin(3);egInvarMassMET30->Rebin(3);ffInvarMassMET30->Rebin(3);ggInvarMassMET40->Rebin(3);ggInvarMassMET40_JetReq->Rebin(3);ggInvarMassMET50->Rebin(3);ggInvarMassMET50_JetReq->Rebin(3);ggInvarMassMET60->Rebin(3);ggInvarMassMET60_JetReq->Rebin(3);ggInvarMassMET70->Rebin(3);ggInvarMassMET70_JetReq->Rebin(3);ggInvarMassMET80->Rebin(3);ggInvarMassMET100->Rebin(3);gammafakeInvarMassMET30->Rebin(3);

  TH1F* egInvarMassMET30_ana = (TH1F*)egInvarMassMET30->Clone();egInvarMassMET30_ana->Sumw2();
  TH1F* ggInvarMassMET30_ana = (TH1F*)ggInvarMassMET30->Clone();ggInvarMassMET30_ana->Sumw2();
  TH1F* ffInvarMassMET30_ana = (TH1F*)ffInvarMassMET30->Clone();ffInvarMassMET30_ana->Sumw2();
  TH1F* gammafakeInvarMassMET30_ana = (TH1F*)gammafakeInvarMassMET30->Clone();gammafakeInvarMassMET30_ana->Sumw2();
  //float ScaleegInvMass=(FakeRate/(1-FakeRate));
  float ScaleegInvMass=(.02/(1-.02));
  egInvarMassMET30_ana->Scale(ScaleegInvMass);

  float gfffscaleinvmass = ffInvarMassMET30_ana->Integral()/gammafakeInvarMassMET30_ana->Integral();
  gammafakeInvarMassMET30_ana->Scale(gfffscaleinvmass);
  ffInvarMassMET30_ana->Scale(ffPercent);gammafakeInvarMassMET30_ana->Scale(gammafakePercent);
  //ffInvarMassMET30_ana->Add(gammafakeInvarMassMET30_ana);

  int bin1=ggInvarMassMET30_ana->FindBin(100),bin2=ggInvarMassMET30_ana->FindBin(110),bin3=ggInvarMassMET30_ana->FindBin(140),bin4=ggInvarMassMET30_ana->FindBin(180);
  //float ScaleffInvMass=(ggInvarMassMET30_ana->Integral(bin1,bin2)-egInvarMassMET30_ana->Integral(bin1,bin2))/ffInvarMassMET30_ana->Integral(bin1,bin2);
  float ggadd=ggInvarMassMET30_ana->Integral(bin1,bin2)+ggInvarMassMET30_ana->Integral(bin3,bin4);
  float egadd=egInvarMassMET30_ana->Integral(bin1,bin2)+egInvarMassMET30_ana->Integral(bin3,bin4);
  float ffadd=ffInvarMassMET30_ana->Integral(bin1,bin2)+ffInvarMassMET30_ana->Integral(bin3,bin4);
  float ggInt_ana = ggInvarMassMET30_ana->Integral() - ggInvarMassMET30_ana->Integral(ggInvarMassMET30_ana->FindBin(120),ggInvarMassMET30_ana->FindBin(130))-ggInvarMassMET30_ana->Integral(ggInvarMassMET30_ana->FindBin(80),ggInvarMassMET30_ana->FindBin(100));
  float egInt_ana = egInvarMassMET30_ana->Integral() - egInvarMassMET30_ana->Integral(egInvarMassMET30_ana->FindBin(120),egInvarMassMET30_ana->FindBin(130))-egInvarMassMET30_ana->Integral(egInvarMassMET30_ana->FindBin(80),egInvarMassMET30_ana->FindBin(100));
  float ffInt_ana = ffInvarMassMET30_ana->Integral() - ffInvarMassMET30_ana->Integral(ffInvarMassMET30_ana->FindBin(120),ffInvarMassMET30_ana->FindBin(130))-ffInvarMassMET30_ana->Integral(ffInvarMassMET30_ana->FindBin(80),ffInvarMassMET30_ana->FindBin(100));
  float gammafakeInt_ana = gammafakeInvarMassMET30_ana->Integral() - gammafakeInvarMassMET30_ana->Integral(gammafakeInvarMassMET30_ana->FindBin(120),gammafakeInvarMassMET30_ana->FindBin(130))-gammafakeInvarMassMET30_ana->Integral(gammafakeInvarMassMET30_ana->FindBin(80),gammafakeInvarMassMET30_ana->FindBin(100));
  float ScaleffInvMass=(ggInt_ana-egInt_ana)/(ffInt_ana+gammafakeInt_ana);//(gg-eg)/ff excluding Z peak (80-100) and Higgs window (120-130)
  //float ScaleffInvMass=(ggadd-egadd)/ffadd;//(gg-eg)/ff using only 100-110 and 140-180
  ffInvarMassMET30_ana->Scale(ScaleffInvMass);gammafakeInvarMassMET30_ana->Scale(ScaleffInvMass);
  
  //egInvarMassMET30_ana->Rebin(2);ggInvarMassMET30_ana->Rebin(2);ffInvarMassMET30_ana->Rebin(2);

  ffInvarMassMET30_ana->SetMarkerSize(0);
  TH1F* bgInvMass = (TH1F*)ffInvarMassMET30_ana->Clone();
  bgInvMass->Add(egInvarMassMET30_ana);
  bgInvMass->Add(gammafakeInvarMassMET30_ana);
  egInvarMassMET30_ana->SetLineColor(kGreen);
  egInvarMassMET30_ana->SetFillColor(kGreen);
  ffInvarMassMET30_ana->SetLineColor(kAzure);
  ffInvarMassMET30_ana->SetFillColor(kAzure);
  gammafakeInvarMassMET30_ana->SetLineColor(kGray+1);
  gammafakeInvarMassMET30_ana->SetFillColor(kGray+1);
  bgInvMass->SetFillColor(kRed);
  bgInvMass->SetFillStyle(3154);
  THStack *InvMassStack = new THStack("InvMassStack",";;Number of Events");
  InvMassStack->Add(egInvarMassMET30_ana);
  InvMassStack->Add(ffInvarMassMET30_ana);
  InvMassStack->Add(gammafakeInvarMassMET30_ana);
  InvMassStack->Draw("histo");
  InvMassStack->SetMaximum(820);
  InvMassStack->GetHistogram()->GetXaxis()->SetTitle("Invariant Mass [GeV]");
  InvMassStack->GetHistogram()->GetXaxis()->SetRangeUser(105,150);
  InvMassStack->Draw("histo");
  bgInvMass->Draw("E2SAME");
  ggInvarMassMET30_ana->Draw("PEsame");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET30cut_RA3bgMethod.png");
  //ggInvarMassMET30->GetXaxis()->SetRangeUser(0,300);
  
  RooRealVar xMet("xMet","m_{#gamma#gamma}",95,200,"GeV");
  xMet.setRange("full",95,200);
  RooDataHist dataMet("dataMet","dataset",xMet,ggInvarMassMET);
  //Gaussian for peak
  RooRealVar gmeanMet("gmeanMet","gmean",125.3,122,128);
  RooRealVar gsigmaMet("gsigmaMet","gsigma",1.6/*,0.1,5*/);
  RooGaussian sigMet("sigMet","gauss",xMet,gmeanMet,gsigmaMet);
  //Crystal Ball for signal
  RooRealVar cbmeanMet("cb_mean", "cbmean" , 125.3, 122, 128.) ;
  RooRealVar cbsigmaMet("cb_sigma", "cbsigma",1.6/*2.6, 0., 3.*/);
  RooRealVar nMet("n","n", 15/*10.,-5.,25.*/);
  RooRealVar alphaMet("alpha","alpha",25./*,0.,5.*/);
  //RooCBShape sigMet("sigMet", "crystal ball", xMet, cbmeanMet, cbsigmaMet, alphaMet, nMet);
  RooRealVar sigMetYield("signal yield","signal yield",100,0,600);
   RooRealVar Bern1("Bern1","Berstein 1",12.,0.,19.);
  RooRealVar Bern2("Bern2","Berstein 2",5.,1.,9.);
  RooRealVar Bern3("Bern3","Berstein 3",4.,0.,8.);
  RooRealVar Bern4("Bern4","Berstein 4",2.,0.,4.);
  RooRealVar Bern5("Bern5","Berstein 5",5.,0.,10.);
  RooRealVar Pol1("Pol1","Pol1",0,-1,1);
  RooRealVar Pol2("Pol2","Pol2",0,-1,1);
  RooRealVar Pol3("Pol3","Pol3",0,-1,1);
  RooRealVar Pol4("Pol4","Pol4",0,-1,1);
  RooRealVar Pol5("Pol5","Pol5",0,-1,1);
  RooBernstein Bern("Bern","4th order Bernstein Polynomial",xMet,RooArgList(Bern1,Bern2,Bern3,Bern4,Bern5));
  //RooPolynomial Bern("Bern","4th Order Polynomial",xMet,RooArgList(Pol1,Pol2,Pol3,Pol4));
  RooRealVar BernYield("bkgd yield","bkgd yield",20000,0,400000);
  RooAddPdf MetPdf("MetPdf","MetPdf",RooArgList(Bern,sigMet),RooArgList(BernYield,sigMetYield));
  RooFitResult *rNoMetCut = MetPdf.fitTo(dataMet,Extended(kTRUE),Save());
  RooPlot *xframeMet = xMet.frame(Title("Crystal Ball Signal, 5th Order Bernstein Polynomial, gg with no MET cut"));
  MetPdf.paramOn(xframeMet, Format("NE",AutoPrecision(2)),Layout(0.5,0.995,0.9) );
  dataMet.plotOn(xframeMet,LineColor(kBlack));
  //MetPdf.plotOn(xframeMet,Components(Bern),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut,3,kTRUE),FillColor(kViolet));
  MetPdf.plotOn(xframeMet,Components(Bern),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut,2,kFALSE),FillColor(kGreen));
  MetPdf.plotOn(xframeMet,Components(Bern),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut,1,kFALSE),FillColor(kOrange));
  dataMet.plotOn(xframeMet,LineColor(kBlack));
  MetPdf.plotOn(xframeMet,Components(Bern),LineColor(kRed),LineStyle(kDashed),Range("full"));
  MetPdf.plotOn(xframeMet,LineColor(kBlue));
  MetPdf.plotOn(xframeMet,Components(sigMet),LineColor(kBlue));
  dataMet.plotOn(xframeMet,LineColor(kBlack));
  xframeMet->SetAxisRange(105.5,154,"X");
  //xframeMet->SetAxisRange(1550,3200,"Y");
  xframeMet->SetMaximum(5000);
  xframeMet->Draw();
  line125->DrawLine(125.3,0,125.3,xframeMet->GetMaximum());
  THStack* StackMet = new THStack("StackMet","");
  StackMet->Add(TTHggNoMET);//->Draw("SAME");
  StackMet->Add(WZHggNoMET);//->Draw("SAME"); 
  StackMet->Add(VBFHggNoMET);//->Draw("SAME");
  StackMet->Add(ggHggNoMET);//->Draw("SAME");
  StackMet->Draw("histoSAME");
  xframeMet->Draw("SAME");
  legHiggs->Draw();
  Double_t sigSM=0.;
  Double_t sigSMErr=0.;
  ggHggNoMET->Add(WZHggNoMET);ggHggNoMET->Add(VBFHggNoMET);ggHggNoMET->Add(TTHggNoMET);
  sigSM = ggHggNoMET->IntegralAndError(0,999,sigSMErr,"");
  //sigSM/=10;sigSMErr/=10;
  Double_t sigYield=sigMetYield.getVal();
  Double_t sigYieldErr=sigMetYield.getError();
  //cout<<"sigSM: "<<sigSM<<"  sigSMErr: "<<sigSMErr<<endl;
  //cout<<"sigYield: "<<sigYield<<"  sigYieldErr: "<<sigYieldErr<<endl;
  sigStrengthVsMet->SetBinContent(1,sigYield/sigSM);
  sigStrengthVsMet->SetBinError(1,sqrt(sigYieldErr*sigYieldErr/(sigSM*sigSM) + (sigYield*sigYield*sigSMErr*sigSMErr)/(sigSM*sigSM)));
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithNoMETcut.png");

  xMet.setRange("sb_lo",95,120);xMet.setRange("sb_hi",130,200);;xMet.setRange("full",95,200);
  RooRealVar Bern1_bg("Bern1","Berstein 1",12.,0.,19.);
  RooRealVar Bern2_bg("Bern2","Berstein 2",5.,1.,9.);
  RooRealVar Bern3_bg("Bern3","Berstein 3",4.,0.,8.);
  RooRealVar Bern4_bg("Bern4","Berstein 4",2.,0.,4.);
  RooRealVar Bern5_bg("Bern5","Berstein 5",5.,0.,10.);
  RooBernstein Bern_bg("Bern","4th order Bernstein Polynomial",xMet,RooArgList(Bern1,Bern2,Bern3,Bern4,Bern5));
  RooRealVar BernYield_bg("bkgd yield","bkgd yield",70000,0,170000);
  RooAddPdf MetPdf_bg("MetPdf_bg","MetPdf_bg",RooArgList(Bern),RooArgList(BernYield));
  //RooFitResult *rNoMetCut_bg = MetPdf_bg.fitTo(dataMet,Extended(kTRUE),Save(),Range("sb_lo,sb_hi"));
  RooFitResult *rNoMetCut_bg_full = MetPdf_bg.fitTo(dataMet,Extended(kTRUE),Save());
  // Print fit results 
   cout << "result of fit on all data " << endl ;
  rNoMetCut_bg_full->Print() ;  
  cout << "result of fit in in background region" << endl ;
  //rNoMetCut_bg->Print() ;
  cout << "result of fit in background region" << endl ;
  RooPlot *xframeMet_bg = xMet.frame(Title("5th Order Polynomial Background, gg with no MET cut"));
  MetPdf.paramOn(xframeMet_bg, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.9) );
  dataMet.plotOn(xframeMet_bg,LineColor(kBlack));
  //MetPdf_bg.plotOn(xframeMet_bg,Components(Bern),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg,3,kFALSE),FillColor(kViolet),);
  //MetPdf_bg.plotOn(xframeMet_bg,Components(Bern),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg,2,kFALSE),FillColor(kGreen));
  //MetPdf_bg.plotOn(xframeMet_bg,Components(Bern),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg,1,kFALSE),FillColor(kOrange));
  //MetPdf_bg.plotOn(xframeMet_bg,Components(Bern),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg_full,3,kFALSE),FillColor(kViolet));
  //MetPdf_bg.plotOn(xframeMet_bg,Components(Bern),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg_full,2,kFALSE),FillColor(kGreen));
  //MetPdf_bg.plotOn(xframeMet_bg,Components(Bern),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_bg_full,1,kFALSE),FillColor(kOrange));
  MetPdf.plotOn(xframeMet_bg,Components(Bern),LineColor(kRed),LineStyle(kDashed));
  MetPdf.plotOn(xframeMet_bg,Components(Bern),LineColor(kRed),Range("sb_lo,sb_hi"));
  dataMet.plotOn(xframeMet_bg,LineColor(kBlack));
  xframeMet_bg->SetAxisRange(105.5,154,"X");
  //xframeMet_bg->SetAxisRange(1550,3200,"Y");
  xframeMet_bg->SetMaximum(5000);
  xframeMet_bg->Draw();
  //line125->DrawLine(125.3,0,125.3,xframeMet_bg->GetMaximum());
  StackMet->Draw("histoSAME");
  xframeMet_bg->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithNoMETcut_bgOnly.png");

  //no met cut jetreq
  RooRealVar xMet_JetReq("xMet_JetReq","gg Invariant Mass with no MET cut",90,180,"GeV");
  RooDataHist dataMet_JetReq("dataMet_JetReq","dataset",xMet_JetReq,ggInvarMassMET_JetReq);
  //Gaussian for peak
  RooRealVar gmeanMet_JetReq("gmeanMet_JetReq","gmean",125.3,122,128);
  RooRealVar gsigmaMet_JetReq("gsigmaMet_JetReq","gsigma",1.6/*,0.1,5*/);
  RooGaussian sigMet_JetReq("sigMet_JetReq","gauss",xMet_JetReq,gmeanMet_JetReq,gsigmaMet_JetReq);
  //Crystal Ball for signal
  RooRealVar cbmeanMet_JetReq("cb_mean_JetReq", "cbmean" , 125.3, 122, 128.) ;
  RooRealVar cbsigmaMet_JetReq("cb_sigma_JetReq", "cbsigma",1.6/*2.6, 0., 3.*/);
  RooRealVar nMet_JetReq("n_JetReq","n", 15/*10.,-5.,25.*/);
  RooRealVar alphaMet_JetReq("alpha_JetReq","alpha",25./*,0.,5.*/);
  //RooCBShape sigMet_JetReq("sigMet_JetReq", "crystal ball", xMet_JetReq, cbmeanMet_JetReq, cbsigmaMet_JetReq, alphaMet_JetReq, nMet_JetReq);
  RooRealVar sigMetYield_JetReq("signal yield","signal yield",10,0,500);
  RooRealVar Bern1_JetReq("Bern1_JetReq","Berstein 1",12.,4.,20.);
  RooRealVar Bern2_JetReq("Bern2_JetReq","Berstein 2",5.,1.,9.);
  RooRealVar Bern3_JetReq("Bern3_JetReq","Berstein 3",4.,0.,8.);
  RooRealVar Bern4_JetReq("Bern4_JetReq","Berstein 4",2.,0.,4.);
  RooRealVar Bern5_JetReq("Bern5_JetReq","Berstein 5",5.,0.,10.);
  RooRealVar Pol1_JetReq("Pol1_JetReq","Pol1",0,-.1,.1);
  RooRealVar Pol2_JetReq("Pol2_JetReq","Pol2",0,-.1,.1);
  RooRealVar Pol3_JetReq("Pol3_JetReq","Pol3",0,-.1,.1);
  RooRealVar Pol4_JetReq("Pol4_JetReq","Pol4",0,-.1,.1);
  RooRealVar Pol5_JetReq("Pol5_JetReq","Pol5",0,-.1,.1);
  //RooBernstein Bern_JetReq("Bern_JetReq","4th order Bernstein Polynomial",xMet_JetReq,RooArgList(Bern1_JetReq,Bern2_JetReq,Bern3_JetReq,Bern4_JetReq/*,Bern5_JetReq*/));
  RooPolynomial Bern_JetReq("Bern_JetReq","4th Order Polynomial",xMet_JetReq,RooArgList(Pol1_JetReq,Pol2_JetReq,Pol3_JetReq,Pol4_JetReq));
  RooRealVar BernYield_JetReq("bkgd yield_JetReq","bkgd yield",20000,10,90000);
  RooAddPdf MetPdf_JetReq("MetPdf_JetReq","MetPdf_JetReq",RooArgList(Bern_JetReq,sigMet_JetReq),RooArgList(BernYield_JetReq,sigMetYield_JetReq));
  RooFitResult *rNoMetCut_JetReq = MetPdf.fitTo(dataMet_JetReq,Extended(kTRUE),Save());
  RooPlot *xframeMet_JetReq = xMet_JetReq.frame(Title("Crystal Ball Signal, 4th Order Bernstein Polynomial Background, gg_JetReq with no MET cut"));
  MetPdf_JetReq.paramOn(xframeMet_JetReq, Format("NE",AutoPrecision(2)),Layout(0.5,0.995,0.9) );
  dataMet_JetReq.plotOn(xframeMet_JetReq,LineColor(kBlack));
  //MetPdf_JetReq.plotOn(xframeMet_JetReq,Components(Bern_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_JetReq,3,kFALSE),FillColor(kViolet));
  //MetPdf_JetReq.plotOn(xframeMet_JetReq,Components(Bern_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_JetReq,2,kFALSE),FillColor(kGreen));
  //MetPdf_JetReq.plotOn(xframeMet_JetReq,Components(Bern_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_JetReq,1,kFALSE),FillColor(kOrange));
  dataMet_JetReq.plotOn(xframeMet_JetReq,LineColor(kBlack));
  MetPdf_JetReq.plotOn(xframeMet_JetReq,Components(Bern_JetReq),LineColor(kRed),LineStyle(kDashed));
  MetPdf_JetReq.plotOn(xframeMet_JetReq,LineColor(kBlue));
  MetPdf_JetReq.plotOn(xframeMet_JetReq,Components(sigMet_JetReq),LineColor(kBlue));
  //xframeMet_JetReq->SetAxisRange(105,154,"X");
  //xframeMet_JetReq->SetAxisRange(1550,3200,"Y");
  xframeMet_JetReq->Draw();
  line125->DrawLine(125.3,0,125.3,xframeMet_JetReq->GetMaximum());
  THStack* StackMet_JetReq = new THStack("StackMet_JetReq","");
  StackMet_JetReq->Add(TTHggNoMET_JetReq);//->Draw("SAME");
  StackMet_JetReq->Add(VBFHggNoMET_JetReq);//->Draw("SAME");
  StackMet_JetReq->Add(WZHggNoMET_JetReq);//->Draw("SAME"); 
  StackMet_JetReq->Add(ggHggNoMET_JetReq);//->Draw("SAME");
  StackMet_JetReq->Draw("histoSAME");
  xframeMet_JetReq->Draw("SAME");
  legHiggs->Draw();
  Double_t sigSM_JetReq=0.;
  Double_t sigSMErr_JetReq=0.;
  ggHggNoMET_JetReq->Add(WZHggNoMET_JetReq);ggHggNoMET_JetReq->Add(VBFHggNoMET_JetReq);ggHggNoMET_JetReq->Add(TTHggNoMET_JetReq);
  sigSM_JetReq = ggHggNoMET_JetReq->IntegralAndError(0,999,sigSMErr_JetReq,"");
  //sigSM_JetReq/=10;sigSMErr_JetReq/=10;
  Double_t sigYield_JetReq=sigMetYield_JetReq.getVal();
  Double_t sigYieldErr_JetReq=sigMetYield_JetReq.getError();
  //cout<<"sigSM_JetReq: "<<sigSM_JetReq<<"  sigSMErr_JetReq: "<<sigSMErr_JetReq<<endl;
  //cout<<"sigYield_JetReq: "<<sigYield_JetReq<<"  sigYieldErr_JetReq: "<<sigYieldErr_JetReq<<endl;
  sigStrengthVsMet_JetReq->SetBinContent(1,sigYield_JetReq/sigSM_JetReq);
  sigStrengthVsMet_JetReq->SetBinError(1,sqrt(sigYieldErr_JetReq*sigYieldErr_JetReq/(sigSM_JetReq*sigSM_JetReq) + (sigYield_JetReq*sigYield_JetReq*sigSMErr_JetReq*sigSMErr_JetReq)/(sigSM_JetReq*sigSM_JetReq)));
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithNoMETcut_JetReq.png");
  
  RooAddPdf MetPdf_bg_JetReq("MetPdf_bg_JetReq","MetPdf_bg",RooArgList(Bern_JetReq),RooArgList(BernYield_JetReq));
  RooFitResult *rNoMetCut_bg_JetReq = MetPdf_bg_JetReq.fitTo(dataMet_JetReq,Extended(kTRUE),Save());
  RooPlot *xframeMet_bg_JetReq = xMet_JetReq.frame(Title("4th Order Bernstein Polynomial Background, gg_JetReq with no MET cut"));
  MetPdf_bg_JetReq.paramOn(xframeMet_bg_JetReq, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.9) );
  dataMet_JetReq.plotOn(xframeMet_bg_JetReq,LineColor(kBlack));
  //MetPdf_JetReq.plotOn(xframeMet_JetReq,Components(Bern_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_JetReq,3,kFALSE),FillColor(kViolet));
  //MetPdf_JetReq.plotOn(xframeMet_JetReq,Components(Bern_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_JetReq,2,kFALSE),FillColor(kGreen));
  //MetPdf_JetReq.plotOn(xframeMet_JetReq,Components(Bern_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*rNoMetCut_JetReq,1,kFALSE),FillColor(kOrange));
  MetPdf_bg_JetReq.plotOn(xframeMet_bg_JetReq,Components(Bern_JetReq),LineColor(kRed),LineStyle(kDashed));
  dataMet_JetReq.plotOn(xframeMet_bg_JetReq,LineColor(kBlack));
  //xframeMet_bg_JetReq->SetAxisRange(105,154,"X");
  //xframeMet_bg_JetReq->SetAxisRange(1550,3200,"Y");
  xframeMet_bg_JetReq->Draw();
  line125->DrawLine(125.3,0,125.3,xframeMet_bg_JetReq->GetMaximum());
  StackMet_JetReq->Draw("histoSAME");
  xframeMet_bg_JetReq->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithNoMETcut_bgOnly_JetReq.png");

  xMet.setRange("xMetSigFull",120,130);
  xMet.setRange("xMetSig",120,121);
  xMet.setRange("xMetSig1",121,122);
  xMet.setRange("xMetSig2",122,123);
  xMet.setRange("xMetSig3",123,124);
  xMet.setRange("xMetSig4",124,125);
  xMet.setRange("xMetSig5",125,126);
  xMet.setRange("xMetSig6",126,127);
  xMet.setRange("xMetSig7",127,128);
  xMet.setRange("xMetSig8",128,129);
  xMet.setRange("xMetSig9",129,130);
  xMet.setRange("xMetLowSB",95,110);
  xMet.setRange("xMetHighSB",135,160);
  RooAbsReal* ig = MetPdf_bg.createIntegral(xMet,NormSet(xMet),Range("xMetSig"));
  RooAbsReal* ig1 = MetPdf_bg.createIntegral(xMet,NormSet(xMet),Range("xMetSig1"));
  RooAbsReal* ig2 = MetPdf_bg.createIntegral(xMet,NormSet(xMet),Range("xMetSig2"));
  RooAbsReal* ig3 = MetPdf_bg.createIntegral(xMet,NormSet(xMet),Range("xMetSig3"));
  RooAbsReal* ig4 = MetPdf_bg.createIntegral(xMet,NormSet(xMet),Range("xMetSig4"));
  RooAbsReal* ig5 = MetPdf_bg.createIntegral(xMet,NormSet(xMet),Range("xMetSig5"));
  RooAbsReal* ig6 = MetPdf_bg.createIntegral(xMet,NormSet(xMet),Range("xMetSig6"));
  RooAbsReal* ig7 = MetPdf_bg.createIntegral(xMet,NormSet(xMet),Range("xMetSig7"));
  RooAbsReal* ig8 = MetPdf_bg.createIntegral(xMet,NormSet(xMet),Range("xMetSig8"));
  RooAbsReal* ig9 = MetPdf_bg.createIntegral(xMet,NormSet(xMet),Range("xMetSig9"));
  RooAbsReal* igSig = MetPdf_bg.createIntegral(xMet,NormSet(xMet),Range("xMetSigFull"));
  RooAbsReal* igLowSB  = MetPdf_bg.createIntegral(xMet,NormSet(xMet),Range("xMetLowSB"));
  RooAbsReal* igHighSB = MetPdf_bg.createIntegral(xMet,NormSet(xMet),Range("xMetHighSB"));
  RooAbsReal* igLowSBNoNorm  = MetPdf_bg.createIntegral(xMet,Range("xMetLowSB"));
  RooAbsReal* igHighSBNoNorm = MetPdf_bg.createIntegral(xMet,Range("xMetHighSB"));
  RooAbsReal* igFull = MetPdf_bg.createIntegral(xMet,NormSet(xMet));
  RooAbsReal* igSigNoNorm = MetPdf_bg.createIntegral(xMet,Range("xMetSigFull"));
  float HiggSigNoNorm = igSig->getVal()*BernYield.getVal();
  float HiggSig = igSig->getVal();
  float HiggHighSBNoNorm = igHighSB->getVal()*BernYield.getVal();
  float HiggLowSBNoNorm = igLowSB->getVal()*BernYield.getVal();
  float HiggSigOverHighSB = igSig->getVal()/igHighSB->getVal();
  float HiggSigOverLowSB  = igSig->getVal()/igLowSB ->getVal();
  cout<<"HiggSig:"<<HiggSig<<"  HiggSig no norm.:"<<HiggSigNoNorm<<"  High SB:"<<igHighSB->getVal()<<"  High SB no norm:"<<HiggHighSBNoNorm<<"  HiggSigOverHighSB:"<<HiggSigOverHighSB<<"  Low SB:"<<igLowSB->getVal()<<"  Low SB no norm:"<<HiggLowSBNoNorm<<"  HiggSigOverLowSB:"<<HiggSigOverLowSB<<endl;
  cout<<"getVal(): "<<ig->getVal()<<endl;
  //cout<<"getError(): "<<ig->getError()<<endl;
  cout<<"getValFull(): "<<igFull->getVal()<<endl;
  TH1F* dataFitComp = new TH1F("dataFitComp","dataFitComp",10,120,130);
  dataFitComp->SetLineColor(kRed);dataFitComp->SetMarkerColor(kRed);
  //dataFitComp->GetYaxis()->SetRangeUser(900,1400);
  dataFitComp->SetBinContent(1,ig->getVal()*BernYield.getVal());
  dataFitComp->SetBinContent(2,ig1->getVal()*BernYield.getVal());
  dataFitComp->SetBinContent(3,ig2->getVal()*BernYield.getVal());
  dataFitComp->SetBinContent(4,ig3->getVal()*BernYield.getVal());
  dataFitComp->SetBinContent(5,ig4->getVal()*BernYield.getVal());
  dataFitComp->SetBinContent(6,ig5->getVal()*BernYield.getVal());
  dataFitComp->SetBinContent(7,ig6->getVal()*BernYield.getVal());
  dataFitComp->SetBinContent(8,ig7->getVal()*BernYield.getVal());
  dataFitComp->SetBinContent(9,ig8->getVal()*BernYield.getVal());
  dataFitComp->SetBinContent(10,ig9->getVal()*BernYield.getVal());
  dataFitComp->Draw("PE");
  TH1F* ggInvarMassClone = (TH1F*)ggInvarMass->Clone();
  ggInvarMassClone->Draw("PESAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMass_gg_FITbg.png");

  RooRealVar xMet30("xMet30","gg Invariant Mass with MET>30",90,180,"GeV");
  RooDataHist dataMet30("dataMet30","dataset",xMet30,ggInvarMassMET30);
  //Gaussian for peak
  RooRealVar gmeanMet30("gmeanMet30","gmean",125.3,122,128);
  RooRealVar gsigmaMet30("gsigmaMet30","gsigma",1,0.1,5);
  //RooGaussian sigMet30("sigMet30","gauss",xMet30,gmeanMet30,gsigmaMet30);
  //Crystal Ball for signal
  RooRealVar cbmeanMet30("cb_mean", "cbmean" , 125.3, 120, 130.) ;
  RooRealVar cbsigmaMet30("cb_sigma", "cbsigma" ,1.8/*,0.,3.*/);
  RooRealVar nMet30("n","n", 2.5,-1.,25.);
  RooRealVar alphaMet30("alpha","alpha",1.35,.5,5.);
  RooCBShape sigMet30("sigMet30", "crystal ball", xMet30, cbmeanMet30, cbsigmaMet30, alphaMet30, nMet30);
  RooRealVar sigMet30Yield("signal yield","signal yield",40,0,1000);
  RooRealVar Bern301("Bern301","Berstein 1",19.,15.,23.);
  RooRealVar Bern302("Bern302","Berstein 2",3.,1.5,4.5);
  RooRealVar Bern303("Bern303","Berstein 3",4.,3.1,4.9);
  RooRealVar Bern304("Bern304","Berstein 4",1.4,1.,1.8);
  RooRealVar Bern305("Bern305","Berstein 5",5.,0.,10.);
  RooRealVar Pol301("Pol301","Pol301",0,-1,1);
  RooRealVar Pol302("Pol302","Pol302",0,-1,1);
  RooRealVar Pol303("Pol303","Pol303",0,-1,1);
  RooRealVar Pol304("Pol304","Pol304",0,-1,1);
  RooRealVar Pol305("Pol305","Pol305",0,-1,1);
  RooBernstein Bern30("Bern30","4th order Bernstein Polynomial",xMet30,RooArgList(Bern301,Bern302,Bern303,Bern304/*,Bern305*/));
  //RooPolynomial Bern30("Bern30","4th Order Polynomial",xMet30,RooArgList(Pol301,Pol302,Pol303,Pol304));
  RooRealVar Bern30Yield("bkgd yield","bkgd yield",7000,0,100000);
  RooAddPdf Met30Pdf("Met30Pdf","Met30Pdf",RooArgList(Bern30,sigMet30),RooArgList(Bern30Yield,sigMet30Yield));
  RooFitResult *r30 = Met30Pdf.fitTo(dataMet30,Extended(kTRUE),Save());
  RooPlot *xframeMet30 = xMet30.frame(Title("Crystal Ball Signal, 4th Order Bernstein Polynomial Background, gg MET>30"));
  Met30Pdf.paramOn(xframeMet30, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.9) );
  dataMet30.plotOn(xframeMet30,LineColor(kBlack));
  //Met30Pdf.plotOn(xframeMet30,Components(Bern30),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r30,3,kFALSE),FillColor(kViolet));
  //Met30Pdf.plotOn(xframeMet30,Components(Bern30),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r30,2,kTRUE),FillColor(kGreen));
  //Met30Pdf.plotOn(xframeMet30,Components(Bern30),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r30,1,kTRUE),FillColor(kOrange));
  dataMet30.plotOn(xframeMet30,LineColor(kBlack));
  Met30Pdf.plotOn(xframeMet30,LineColor(kBlue));
  Met30Pdf.plotOn(xframeMet30,Components(Bern30),LineColor(kRed),LineStyle(kDashed));
  Met30Pdf.plotOn(xframeMet30,Components(sigMet30),LineColor(kBlue));
  xframeMet30->SetMinimum(0);
  xframeMet30->SetAxisRange(105,154,"X");
  xframeMet30->Draw();
  //line125->DrawLine(125.3,0,125.3,xframeMet30->GetMaximum());
  double chi2Met30=xframeMet30->chiSquare();
  //cout<<"MET30 Chi2: "<<chi2Met30<<endl;
  THStack* StackMet30 = new THStack("StackMet30","");
  StackMet30->Add(TTHggMET30);//->Draw("SAME");
  StackMet30->Add(VBFHggMET30);//->Draw("SAME");
  StackMet30->Add(WZHggMET30);//->Draw("SAME"); 
  StackMet30->Add(ggHggMET30);//->Draw("SAME");
  StackMet30->Draw("histoSAME");
  xframeMet30->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET30cut.png");

  RooRealVar xMet30_bg("xMet30_bg","gg Invariant Mass with MET>30",90,180,"GeV");
  RooDataHist dataMet30_bg("dataMet30_bg","dataset",xMet30_bg,ggInvarMassMET30);
  xMet30.setRange("sb_lo30",100,120);xMet30.setRange("sb_hi30",130,180);
  RooRealVar Bern301_bg("Bern301","Berstein 1",19.,15.,23.);
  RooRealVar Bern302_bg("Bern302","Berstein 2",3.,1.5,4.5);
  RooRealVar Bern303_bg("Bern303","Berstein 3",4.,3.1,4.9);
  RooRealVar Bern304_bg("Bern304","Berstein 4",1.4,1.,1.8);
  RooRealVar Bern305_bg("Bern305","Berstein 5",5.,0.,10.);
  RooBernstein Bern30_bg("Bern30","4th order Bernstein Polynomial",xMet30_bg,RooArgList(Bern301_bg,Bern302_bg,Bern303_bg,Bern304_bg/*,Bern305_bg*/));
  RooRealVar Bern30Yield_bg("bkgd yield","bkgd yield",7000,0,100000);
  RooAddPdf Met30Pdf_bg("Met30Pdf_bg","Met30Pdf_bg",RooArgList(Bern30_bg),RooArgList(Bern30Yield_bg));
  RooFitResult *r30Cut_bg = Met30Pdf_bg.fitTo(dataMet30_bg,Extended(kTRUE),Save(),Range("sb_lo30,sb_hi30"));
  RooFitResult *r30Cut_bg_full = Met30Pdf_bg.fitTo(dataMet30_bg,Extended(kTRUE),Save());
  RooPlot *xframeMet30_bg = xMet30_bg.frame(Title("4th Order Bernstein Polynomial Background, gg with MET>30 cut"));
  Met30Pdf_bg.paramOn(xframeMet30_bg, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.9) );
  dataMet30_bg.plotOn(xframeMet30_bg,LineColor(kBlack));
  //Met30Pdf_bg.plotOn(xframeMet30_bg,Components(Bern30_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r30Cut_bg,3,kFALSE),FillColor(kViolet),Range(100,180));
  Met30Pdf_bg.plotOn(xframeMet30_bg,Components(Bern30_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r30Cut_bg,2,kFALSE),FillColor(kGreen),Range(100,180));
  Met30Pdf_bg.plotOn(xframeMet30_bg,Components(Bern30_bg),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r30Cut_bg,1,kFALSE),FillColor(kOrange),Range(100,180));
  /*Met30Pdf_bg.plotOn(xframeMet30_bg,Components(Bern30),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r30Cut_bg_full,3,kFALSE),FillColor(kViolet));
  Met30Pdf_bg.plotOn(xframeMet30_bg,Components(Bern30),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r30Cut_bg_full,2,kFALSE),FillColor(kGreen));
  Met30Pdf_bg.plotOn(xframeMet30_bg,Components(Bern30),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r30Cut_bg_full,1,kFALSE),FillColor(kOrange));*/
  Met30Pdf_bg.plotOn(xframeMet30_bg,Components(Bern30_bg),LineColor(kBlue),LineStyle(kDashed),Range(100,180));
  Met30Pdf_bg.plotOn(xframeMet30_bg,Components(Bern30_bg),LineColor(kBlue));
  dataMet30_bg.plotOn(xframeMet30_bg,LineColor(kBlack));
  xframeMet30_bg->SetAxisRange(105.5,154,"X");
  xframeMet30_bg->SetMaximum(875);
  //xframeMet30_bg->SetAxisRange(1550,3200,"Y");
  xframeMet30_bg->Draw();
  //line125->DrawLine(125.3,0,125.3,xframeMet30_bg->GetMaximum());
  StackMet30->Draw("histoSAME");
  xframeMet30_bg->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET30cut_bgOnly.png");
  sigSM=0.;
  sigSMErr=0.;
  ggHggMET30->Add(WZHggMET30);ggHggMET30->Add(VBFHggMET30);ggHggMET30->Add(TTHggMET30);
  sigSM = ggHggMET30->IntegralAndError(0,999,sigSMErr,"");
  //sigSM/=10;sigSMErr/=10;
  sigYield=sigMet30Yield.getVal();
  sigYieldErr=sigMet30Yield.getError();
  //cout<<"sigSM: "<<sigSM<<"  sigSMErr: "<<sigSMErr<<endl;
  //cout<<"sigYield: "<<sigYield<<"  sigYieldErr: "<<sigYieldErr<<endl;
  sigStrengthVsMet->SetBinContent(4,sigYield/sigSM);
  sigStrengthVsMet->SetBinError(4,sqrt(sigYieldErr*sigYieldErr/(sigSM*sigSM) + (sigYield*sigYield*sigSMErr*sigSMErr)/(sigSM*sigSM)));
  xMet30.setRange("xMet30SigFull",120,130);
  xMet30.setRange("xMet30Sig",120,121);
  xMet30.setRange("xMet30Sig1",121,122);
  xMet30.setRange("xMet30Sig2",122,123);
  xMet30.setRange("xMet30Sig3",123,124);
  xMet30.setRange("xMet30Sig4",124,125);
  xMet30.setRange("xMet30Sig5",125,126);
  xMet30.setRange("xMet30Sig6",126,127);
  xMet30.setRange("xMet30Sig7",127,128);
  xMet30.setRange("xMet30Sig8",128,129);
  xMet30.setRange("xMet30Sig9",129,130);
  RooAbsReal* ig30 = Bern30.createIntegral(xMet30,NormSet(xMet30),Range("xMet30Sig"));
  RooAbsReal* ig301 = Bern30.createIntegral(xMet30,NormSet(xMet30),Range("xMet30Sig1"));
  RooAbsReal* ig302 = Bern30.createIntegral(xMet30,NormSet(xMet30),Range("xMet30Sig2"));
  RooAbsReal* ig303 = Bern30.createIntegral(xMet30,NormSet(xMet30),Range("xMet30Sig3"));
  RooAbsReal* ig304 = Bern30.createIntegral(xMet30,NormSet(xMet30),Range("xMet30Sig4"));
  RooAbsReal* ig305 = Bern30.createIntegral(xMet30,NormSet(xMet30),Range("xMet30Sig5"));
  RooAbsReal* ig306 = Bern30.createIntegral(xMet30,NormSet(xMet30),Range("xMet30Sig6"));
  RooAbsReal* ig307 = Bern30.createIntegral(xMet30,NormSet(xMet30),Range("xMet30Sig7"));
  RooAbsReal* ig308 = Bern30.createIntegral(xMet30,NormSet(xMet30),Range("xMet30Sig8"));
  RooAbsReal* ig309 = Bern30.createIntegral(xMet30,NormSet(xMet30),Range("xMet30Sig9"));
  RooAbsReal* ig30Full = Bern30.createIntegral(xMet30,NormSet(xMet30));
  cout<<"getVal(): "<<ig30->getVal()<<endl;
  //cout<<"getError(): "<<ig30->getError()<<endl;
  cout<<"getValFull(): "<<ig30Full->getVal()<<endl;
  TH1F* dataFitComp30 = new TH1F("dataFitComp30","dataFitComp30",10,120,130);
  dataFitComp30->SetLineColor(kRed);dataFitComp30->SetMarkerColor(kRed);
  //dataFitComp30->GetYaxis()->SetRangeUser(0,500);
  dataFitComp30->SetBinContent(1,ig30->getVal()*Bern30Yield.getVal());
  dataFitComp30->SetBinContent(2,ig301->getVal()*Bern30Yield.getVal());
  dataFitComp30->SetBinContent(3,ig302->getVal()*Bern30Yield.getVal());
  dataFitComp30->SetBinContent(4,ig303->getVal()*Bern30Yield.getVal());
  dataFitComp30->SetBinContent(5,ig304->getVal()*Bern30Yield.getVal());
  dataFitComp30->SetBinContent(6,ig305->getVal()*Bern30Yield.getVal());
  dataFitComp30->SetBinContent(7,ig306->getVal()*Bern30Yield.getVal());
  dataFitComp30->SetBinContent(8,ig307->getVal()*Bern30Yield.getVal());
  dataFitComp30->SetBinContent(9,ig308->getVal()*Bern30Yield.getVal());
  dataFitComp30->SetBinContent(10,ig309->getVal()*Bern30Yield.getVal());
  dataFitComp30->Draw("PE");
  //TH1F* ggInvarMassCloneMET30 = (TH1F*)ggInvarMassMET30->Clone();
  ggInvarMassMET30Clone->Draw("PESAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassMET30_gg_FITbg.png");
 

  RooRealVar xMet30_JetReq("xMet30","gg Invariant Mass with MET>30",90,180,"GeV");
  RooDataHist dataMet30_JetReq("dataMet30","dataset",xMet30_JetReq,ggInvarMassMET30_JetReq);
  //Gaussian for peak
  RooRealVar gmeanMet30_JetReq("gmeanMet30","gmean",125.3,122,128);
  RooRealVar gsigmaMet30_JetReq("gsigmaMet30","gsigma",1,0.1,5);
  //RooGaussian sigMet30("sigMet30","gauss",xMet30,gmeanMet30,gsigmaMet30);
  //Crystal Ball for signal
  RooRealVar cbmeanMet30_JetReq("cb_mean", "cbmean" , 125.3, 122, 128) ;
  RooRealVar cbsigmaMet30_JetReq("cb_sigma", "cbsigma" , 1., 0.0, 3.);
  RooRealVar nMet30_JetReq("n","n", 2.2/*5.,-5.,15.*/);
  RooRealVar alphaMet30_JetReq("alpha","alpha",1.35/*,0.,3.*/);
  RooCBShape sigMet30_JetReq("sigMet30", "crystal ball", xMet30_JetReq, cbmeanMet30_JetReq, cbsigmaMet30_JetReq, alphaMet30_JetReq, nMet30_JetReq);
  RooRealVar sigMet30Yield_JetReq("signal yield","signal yield",30,0,500);
  RooRealVar Bern301_JetReq("Bern301","Berstein 1",18.,14.,22.);
  RooRealVar Bern302_JetReq("Bern302","Berstein 2",5.3,2.,7.6);
  RooRealVar Bern303_JetReq("Bern303","Berstein 3",5.,2.,8.);
  RooRealVar Bern304_JetReq("Bern304","Berstein 4",2.,0.,5.);
  RooRealVar Bern305_JetReq("Bern305","Berstein 5",5.,0.,10.);
  RooBernstein Bern30_JetReq("Bern30_JetReq","4th order Bernstein Polynomial",xMet30_JetReq,RooArgList(Bern301_JetReq,Bern302_JetReq,Bern303_JetReq,Bern304_JetReq/*,Bern305*/));
  RooRealVar Bern30Yield_JetReq("bkgd yield","bkgd yield",2500,0,100000);
  RooAddPdf Met30Pdf_JetReq("Met30Pdf_JetReq","Met30Pdf_JetReq",RooArgList(Bern30_JetReq,sigMet30_JetReq),RooArgList(Bern30Yield_JetReq,sigMet30Yield_JetReq));
  RooFitResult *r30_JetReq = Met30Pdf_JetReq.fitTo(dataMet30_JetReq,Extended(kTRUE),Save());
  RooPlot *xframeMet30_JetReq = xMet30_JetReq.frame(Title("Crystal Ball Signal, 4th Order Bernstein Polynomial Background, gg_JetReq MET>30"));
  Met30Pdf_JetReq.paramOn(xframeMet30_JetReq, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.9) );
  dataMet30_JetReq.plotOn(xframeMet30_JetReq,LineColor(kBlack));
  //Met30Pdf_JetReq.plotOn(xframeMet30_JetReq,Components(Bern30_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r30_JetReq,3,kFALSE),FillColor(kViolet));
  //Met30Pdf_JetReq.plotOn(xframeMet30_JetReq,Components(Bern30_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r30_JetReq,2,kFALSE),FillColor(kGreen));
  //Met30Pdf_JetReq.plotOn(xframeMet30_JetReq,Components(Bern30_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r30_JetReq,1,kFALSE),FillColor(kOrange));
  dataMet30_JetReq.plotOn(xframeMet30_JetReq,LineColor(kBlack));
  Met30Pdf_JetReq.plotOn(xframeMet30_JetReq,Components(Bern30_JetReq),LineColor(kRed),LineStyle(kDashed));
  Met30Pdf_JetReq.plotOn(xframeMet30_JetReq,LineColor(kBlue));
  Met30Pdf_JetReq.plotOn(xframeMet30_JetReq,Components(sigMet30_JetReq),LineColor(kBlue));
  xframeMet30_JetReq->SetMinimum(0);
  xframeMet30_JetReq->SetAxisRange(105,154,"X");
  xframeMet30_JetReq->Draw();
  //line125->DrawLine(125.3,0,125.3,xframeMet30_JetReq->GetMaximum());
  double chi2Met30_JetReq=xframeMet30_JetReq->chiSquare();
  //cout<<"MET30_JetReq Chi2: "<<chi2Met30_JetReq<<endl;
   THStack* StackMet30_JetReq = new THStack("StackMet30_JetReq","");
   StackMet30_JetReq->Add(TTHggMET30_JetReq);//->Draw("SAME");
   StackMet30_JetReq->Add(VBFHggMET30_JetReq);//->Draw("SAME");
   StackMet30_JetReq->Add(WZHggMET30_JetReq);//->Draw("SAME"); 
   StackMet30_JetReq->Add(ggHggMET30_JetReq);//->Draw("SAME");
   StackMet30_JetReq->Draw("histoSAME");
   xframeMet30_JetReq->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET30cut_JetReq.png");

 
  RooRealVar xMet30_JetReq_bg("xMet30","gg Invariant Mass with MET>30",90,180,"GeV");
  RooDataHist dataMet30_JetReq_bg("dataMet30","dataset",xMet30_JetReq_bg,ggInvarMassMET30_JetReq);
  xMet30_JetReq_bg.setRange("sb_lo30jb",100,120);xMet30_JetReq_bg.setRange("sb_hi30jb",130,180);
  RooRealVar Bern301_JetReq_bg("Bern301","Berstein 1",18.,14.,22.);
  RooRealVar Bern302_JetReq_bg("Bern302","Berstein 2",5.3,2.,7.6);
  RooRealVar Bern303_JetReq_bg("Bern303","Berstein 3",5.,2.,8.);
  RooRealVar Bern304_JetReq_bg("Bern304","Berstein 4",2.,0.,5.);
  RooRealVar Bern305_JetReq_bg("Bern305","Berstein 5",5.,0.,10.);
  RooBernstein Bern30_JetReq_bg("Bern30_JetReq","4th order Bernstein Polynomial",xMet30_JetReq_bg,RooArgList(Bern301_JetReq_bg,Bern302_JetReq_bg,Bern303_JetReq_bg,Bern304_JetReq_bg/*,Bern305_JetReq_bg*/));
  RooRealVar Bern30Yield_JetReq_bg("bkgd yield","bkgd yield",2500,0,100000);
  RooAddPdf Met30Pdf_JetReq_bg("Met30Pdf_JetReq_bg","Met30Pdf_JetReq_bg",RooArgList(Bern30_JetReq_bg),RooArgList(Bern30Yield_JetReq_bg));
  RooFitResult *r30Cut_JetReq_bg = Met30Pdf_JetReq_bg.fitTo(dataMet30_JetReq_bg,Extended(kTRUE),Save(),Range("sb_lojb,sb_hijb"));
  RooPlot *xframeMet30_JetReq_bg = xMet30_JetReq_bg.frame(Title("4th Order Bernstein Polynomial Background, gg_JetReq with MET>30 cut"));
  Met30Pdf_JetReq_bg.paramOn(xframeMet30_JetReq_bg, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.9) );
  dataMet30_JetReq_bg.plotOn(xframeMet30_JetReq_bg,LineColor(kBlack));
  //Met30Pdf_JetReq_bg.plotOn(xframeMet30_JetReq_bg,Components(Bern30_JetReq_bg),LineColor(kBlue),LineStyle(kDashed),VisualizeError(*r30Cut_JetReq_bg,3,kFALSE),FillColor(kViolet));
  Met30Pdf_JetReq_bg.plotOn(xframeMet30_JetReq_bg,Components(Bern30_JetReq_bg),LineColor(kBlue),LineStyle(kDashed),VisualizeError(*r30Cut_JetReq_bg,2,kFALSE),FillColor(kGreen));
  Met30Pdf_JetReq_bg.plotOn(xframeMet30_JetReq_bg,Components(Bern30_JetReq_bg),LineColor(kBlue),LineStyle(kDashed),VisualizeError(*r30Cut_JetReq_bg,1,kFALSE),FillColor(kOrange));
  Met30Pdf_JetReq_bg.plotOn(xframeMet30_JetReq_bg,Components(Bern30_JetReq_bg),LineColor(kBlue),LineStyle(kDashed),Range(100,180));
  Met30Pdf_JetReq_bg.plotOn(xframeMet30_JetReq_bg,Components(Bern30_JetReq_bg),LineColor(kBlue));
  dataMet30_JetReq_bg.plotOn(xframeMet30_JetReq_bg,LineColor(kBlack));
  xframeMet30_JetReq_bg->SetAxisRange(105,154,"X");
  //xframeMet30_bg->SetAxisRange(1550,3200,"Y");
  xframeMet30_JetReq_bg->SetMaximum(250);
  xframeMet30_JetReq_bg->Draw();
  //line125->DrawLine(125.3,0,125.3,xframeMet30_JetReq_bg->GetMaximum());
  StackMet30_JetReq->Draw("histoSAME");
  xframeMet30_JetReq_bg->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET30cut_JetReq_bgOnly.png");

  RooRealVar xMet40("xMet40","gg Invariant Mass with MET>40",90,180,"GeV");
  RooDataHist dataMet40("dataMet40","dataset",xMet40,ggInvarMassMET40);
  //Gaussian for peak
  RooRealVar gmeanMet40("gmeanMet40","gmean",125.3,122,128);
  RooRealVar gsigmaMet40("gsigmaMet40","gsigma",1,0.1,5);
  //RooGaussian sigMet40("sigMet40","gauss",xMet40,gmeanMet40,gsigmaMet40);
  //Crystal Ball for signal
  RooRealVar cbmeanMet40("cb_mean","cbmean",125.3,122.,128.) ;
  RooRealVar cbsigmaMet40("cb_sigma","cbsigma",1.5,0.5,5.5) ;
  RooRealVar nMet40("n","n",2.5/*,-5.,20.*/);
  RooRealVar alphaMet40("alpha","alpha",1.35,.1,5);
  RooCBShape sigMet40("sigMet40", "crystal ball", xMet40, cbmeanMet40, cbsigmaMet40, alphaMet40, nMet40);
  RooRealVar sigMet40Yield("signal yield","signal yield",20,0,500);
  RooRealVar Bern401("Bern401","Berstein 1",12.,6.,18.);
  RooRealVar Bern402("Bern402","Berstein 2",4.,1.,7.);
  RooRealVar Bern403("Bern403","Berstein 3",3.,1.,5.);
  RooRealVar Bern404("Bern404","Berstein 4",2.,0.,4.);
  RooRealVar Bern405("Bern405","Berstein 5",5.,0.,10.);
  RooBernstein Bern40("Bern40","4th order Bernstein Polynomial",xMet40,RooArgList(Bern401,Bern402,Bern403,Bern404/*,Bern405*/));
  RooRealVar Bern40Yield("bkgd yield","bkgd yield",1000,0,100000);
  RooAddPdf Met40Pdf("Met40Pdf","Met40Pdf",RooArgList(Bern40,sigMet40),RooArgList(Bern40Yield,sigMet40Yield));
  RooFitResult *r40 = Met40Pdf.fitTo(dataMet40,Extended(kTRUE),Save());
  RooPlot *xframeMet40 = xMet40.frame(Title("Crystal Ball Signal, 4th Order Bernstein Polynomial Background, gg MET>40"));
  Met40Pdf.paramOn(xframeMet40, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.9) );
  dataMet40.plotOn(xframeMet40,LineColor(kBlack));
  //Met40Pdf.plotOn(xframeMet40,Components(Bern40),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r40,3,kFALSE),FillColor(kViolet));
  //Met40Pdf.plotOn(xframeMet40,Components(Bern40),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r40,2,kFALSE),FillColor(kGreen));
  //Met40Pdf.plotOn(xframeMet40,Components(Bern40),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r40,1,kFALSE),FillColor(kOrange));
  dataMet40.plotOn(xframeMet40,LineColor(kBlack));
  Met40Pdf.plotOn(xframeMet40,Components(Bern40),LineColor(kRed),LineStyle(kDashed));
  Met40Pdf.plotOn(xframeMet40,LineColor(kBlue));
  Met40Pdf.plotOn(xframeMet40,Components(sigMet40),LineColor(kBlue));
  xframeMet40->SetMinimum(0);
  xframeMet40->Draw();
  line125->DrawLine(125.3,0,125.3,xframeMet40->GetMaximum());
  double chi2Met40=xframeMet40->chiSquare();
  //cout<<"MET40 Chi2: "<<chi2Met40<<endl;
  THStack* StackMet40 = new THStack("StackMet40","");
  StackMet40->Add(TTHggMET40);//->Draw("SAME");
  StackMet40->Add(VBFHggMET40);//->Draw("SAME");
  StackMet40->Add(WZHggMET40);//->Draw("SAME"); 
  StackMet40->Add(ggHggMET40);//->Draw("SAME");
  StackMet40->Draw("histoSAME");
  xframeMet40->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET40cut.png");
  ggHggMET40->Add(WZHggMET40);ggHggMET40->Add(VBFHggMET40);ggHggMET40->Add(TTHggMET40);
  sigSM = ggHggMET40->IntegralAndError(0,999,sigSMErr,"");
  //sigSM/=10;sigSMErr/=10;
  sigYield=sigMet40Yield.getVal();
  sigYieldErr=sigMet40Yield.getError();
  //cout<<"sigSM: "<<sigSM<<"  sigSMErr: "<<sigSMErr<<endl;
  //cout<<"sigYield: "<<sigYield<<"  sigYieldErr: "<<sigYieldErr<<endl;
  sigStrengthVsMet->SetBinContent(5,sigYield/sigSM);
  sigStrengthVsMet->SetBinError(5,sqrt(sigYieldErr*sigYieldErr/(sigSM*sigSM) + (sigYield*sigYield*sigSMErr*sigSMErr)/(sigSM*sigSM)));

  RooAddPdf Met40Pdf_bg("Met40Pdf_bg","Met40Pdf_bg",RooArgList(Bern40),RooArgList(Bern40Yield));
  RooFitResult *r40_bg = Met40Pdf_bg.fitTo(dataMet40,Extended(kTRUE),Save());
  RooPlot *xframeMet40_bg = xMet40.frame(Title("4th Order Bernstein Polynomial Background, gg MET>40 bground only model"));
  Met40Pdf_bg.paramOn(xframeMet40_bg, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.9) );
  dataMet40.plotOn(xframeMet40_bg,LineColor(kBlack));
  //Met40Pdf_bg.plotOn(xframeMet40_bg,LineColor(kRed),LineStyle(kDashed),VisualizeError(*r40_bg,3,kFALSE),FillColor(kViolet));
  //Met40Pdf_bg.plotOn(xframeMet40_bg,LineColor(kRed),LineStyle(kDashed),VisualizeError(*r40_bg,2,kFALSE),FillColor(kGreen)); 
  //Met40Pdf_bg.plotOn(xframeMet40_bg,LineColor(kRed),LineStyle(kDashed),VisualizeError(*r40_bg,1,kFALSE),FillColor(kOrange));
  dataMet40.plotOn(xframeMet40_bg,LineColor(kBlack));
  xframeMet40_bg->SetMinimum(0);
  xframeMet40_bg->Draw();
  line125->DrawLine(125.3,0,125.3,xframeMet40_bg->GetMaximum());
  xframeMet40_bg->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET40cut_bgOnly.png");

  RooRealVar xMet40_JetReq("xMet40","gg Invariant Mass with MET>40",90,180,"GeV");
  RooDataHist dataMet40_JetReq("dataMet40","dataset",xMet40_JetReq,ggInvarMassMET40_JetReq);
  //Gaussian for peak
  RooRealVar gmeanMet40_JetReq("gmeanMet40","gmean",125.3,122,128);
  RooRealVar gsigmaMet40_JetReq("gsigmaMet40","gsigma",1,0.1,5);
  //RooGaussian sigMet40("sigMet40","gauss",xMet40,gmeanMet40,gsigmaMet40);
  //Crystal Ball for signal
  RooRealVar cbmeanMet40_JetReq("cb_mean", "cbmean" , 125.3, 122, 128) ;
  RooRealVar cbsigmaMet40_JetReq("cb_sigma", "cbsigma" , 1.3/*1.5, 0.0, 3.*/);
  RooRealVar nMet40_JetReq("n","n", 3./*5.,-5.,15.*/);
  RooRealVar alphaMet40_JetReq("alpha","alpha",1.35/*,0.,3.*/);
  RooCBShape sigMet40_JetReq("sigMet40", "crystal ball", xMet40_JetReq, cbmeanMet40_JetReq, cbsigmaMet40_JetReq, alphaMet40_JetReq, nMet40_JetReq);
  RooRealVar sigMet40Yield_JetReq("signal yield","signal yield",40,0,500);
  RooRealVar Bern401_JetReq("Bern401","Berstein 1",16.,9.,23.);
  RooRealVar Bern402_JetReq("Bern402","Berstein 2",5.3,2.,7.6);
  RooRealVar Bern403_JetReq("Bern403","Berstein 3",5.,2.,8.);
  RooRealVar Bern404_JetReq("Bern404","Berstein 4",2.,0.,5.);
  RooRealVar Bern405_JetReq("Bern405","Berstein 5",5.,0.,10.);
  RooBernstein Bern40_JetReq("Bern40_JetReq","4th order Bernstein Polynomial",xMet40_JetReq,RooArgList(Bern401_JetReq,Bern402_JetReq,Bern403_JetReq,Bern404_JetReq/*,Bern405*/));
  RooRealVar Bern40Yield_JetReq("bkgd yield","bkgd yield",2500,0,100000);
  RooAddPdf Met40Pdf_JetReq("Met40Pdf_JetReq","Met40Pdf_JetReq",RooArgList(Bern40_JetReq,sigMet40_JetReq),RooArgList(Bern40Yield_JetReq,sigMet40Yield_JetReq));
  RooFitResult *r40_JetReq = Met40Pdf_JetReq.fitTo(dataMet40_JetReq,Extended(kTRUE),Save());
  RooPlot *xframeMet40_JetReq = xMet40_JetReq.frame(Title("Crystal Ball Signal, 4th Order Bernstein Polynomial Background, gg_JetReq MET>40"));
  Met40Pdf_JetReq.paramOn(xframeMet40_JetReq, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.9) );
  dataMet40_JetReq.plotOn(xframeMet40_JetReq,LineColor(kBlack));
  //Met40Pdf_JetReq.plotOn(xframeMet40_JetReq,Components(Bern40_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r40_JetReq,3,kFALSE),FillColor(kViolet));
  //Met40Pdf_JetReq.plotOn(xframeMet40_JetReq,Components(Bern40_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r40_JetReq,2,kFALSE),FillColor(kGreen));
  //Met40Pdf_JetReq.plotOn(xframeMet40_JetReq,Components(Bern40_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r40_JetReq,1,kFALSE),FillColor(kOrange));
  dataMet40_JetReq.plotOn(xframeMet40_JetReq,LineColor(kBlack));
  Met40Pdf_JetReq.plotOn(xframeMet40_JetReq,Components(Bern40_JetReq),LineColor(kRed),LineStyle(kDashed));
  Met40Pdf_JetReq.plotOn(xframeMet40_JetReq,LineColor(kBlue));
  //Met40Pdf_JetReq.plotOn(xframeMet40_JetReq,Components(sigMet40_JetReq),LineColor(kBlue));
  xframeMet40_JetReq->SetMinimum(0);
  xframeMet40_JetReq->Draw();
  line125->DrawLine(125.3,0,125.3,xframeMet40_JetReq->GetMaximum());
  double chi2Met40_JetReq=xframeMet40_JetReq->chiSquare();
  //cout<<"MET40_JetReq Chi2: "<<chi2Met40_JetReq<<endl;
  /* THStack* StackMet40_JetReq = new THStack("StackMet40_JetReq","");
     StackMet40_JetReq->Add(TTHggMET40_JetReq);//->Draw("SAME");
     StackMet40_JetReq->Add(VBFHggMET40_JetReq);//->Draw("SAME");
     StackMet40_JetReq->Add(WZHggMET40_JetReq);//->Draw("SAME"); 
     StackMet40_JetReq->Add(ggHggMET40_JetReq);//->Draw("SAME");
     StackMet40_JetReq->Draw("histoSAME");
     xframeMet40_JetReq->Draw("SAME");*/
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET40cut_JetReq.png");
  RooAddPdf Met40Pdf_JetReq_bg("Met40Pdf_JetReq_bg","Met40Pdf_JetReq_bg",RooArgList(Bern40_JetReq),RooArgList(Bern40Yield_JetReq));
  Met40Pdf_JetReq_bg.fitTo(dataMet40_JetReq,Extended(kTRUE));
  RooPlot *xframeMet40_JetReq_bg = xMet40_JetReq.frame(Title("4th Order Bernstein Polynomial Background, gg_JetReq MET>40"));
  Met40Pdf_JetReq_bg.paramOn(xframeMet40_JetReq_bg, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.9) );
  dataMet40_JetReq.plotOn(xframeMet40_JetReq_bg,LineColor(kBlack));
  Met40Pdf_JetReq_bg.plotOn(xframeMet40_JetReq_bg,LineColor(kRed),LineStyle(kDashed));
  xframeMet40_JetReq_bg->SetMinimum(0);
  xframeMet40_JetReq_bg->Draw();
  line125->DrawLine(125.3,0,125.3,xframeMet40_JetReq_bg->GetMaximum());
  xframeMet40_JetReq_bg->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET40cut_JetReq_bgOnly.png");

  RooRealVar xMet50("xMet50","gg Invariant Mass with MET>50",90,180,"GeV");
  RooDataHist dataMet50("dataMet50","dataset",xMet50,ggInvarMassMET50);
  //Gaussian for peak
  RooRealVar gmeanMet50("gmeanMet50","gmean",125.3,122,128);
  RooRealVar gsigmaMet50("gsigmaMet50","gsigma",1,0.1,5);
  //RooGaussian sigMet50("sigMet50","gauss",xMet50,gmeanMet50,gsigmaMet50);
  //Crystal Ball for signal
  RooRealVar cbmeanMet50("cb_mean", "cbmean" , 125.3, 122., 128.) ;
  RooRealVar cbsigmaMet50("cb_sigma", "cbsigma" , 1.5, 0.0, 6.) ;
  RooRealVar nMet50("n","n",2.5/*,-5.,20.*/);
  RooRealVar alphaMet50("alpha","alpha",1.8,1.0,3.6);
  RooCBShape sigMet50("sigMet50", "crystal ball", xMet50, cbmeanMet50, cbsigmaMet50, alphaMet50, nMet50);
  RooRealVar sigMet50Yield("signal yield","signal yield",20,0,500);
  RooRealVar Bern501("Bern501","Berstein 1",12.,8.,16.);
  RooRealVar Bern502("Bern502","Berstein 2",2.,0.,5.);
  RooRealVar Bern503("Bern503","Berstein 3",3.,1.,5.);
  RooRealVar Bern504("Bern504","Berstein 4",1.5,0.,3.);
  RooRealVar Bern505("Bern505","Berstein 5",5.,0.,10.);
  RooBernstein Bern50("Bern50","4th order Bernstein Polynomial",xMet50,RooArgList(Bern501,Bern502,Bern503,Bern504/*,Bern505*/));
  RooRealVar Bern50Yield("bkgd yield","bkgd yield",500,0,100000);
  RooAddPdf Met50Pdf("Met50Pdf","Met50Pdf",RooArgList(Bern50,sigMet50),RooArgList(Bern50Yield,sigMet50Yield));
  RooFitResult *r50 = Met50Pdf.fitTo(dataMet50,Extended(kTRUE),Save());
  RooPlot *xframeMet50 = xMet50.frame(Title("Crystal Ball Signal, 4th Order Bernstein Polynomial Background, gg MET>50"));
  Met50Pdf.paramOn(xframeMet50, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.9) );
  dataMet50.plotOn(xframeMet50,LineColor(kBlack));
  /*Met50Pdf.plotOn(xframeMet50,Components(Bern50),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r50,5),FillColor(kRed));
    Met50Pdf.plotOn(xframeMet50,Components(Bern50),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r50,4),FillColor(kCyan));*/
  //Met50Pdf.plotOn(xframeMet50,Components(Bern50),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r50,3,kFALSE),FillColor(kViolet));
  //Met50Pdf.plotOn(xframeMet50,Components(Bern50),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r50,2,kFALSE),FillColor(kGreen));
  //Met50Pdf.plotOn(xframeMet50,Components(Bern50),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r50,1,kFALSE),FillColor(kOrange));
  dataMet50.plotOn(xframeMet50,LineColor(kBlack));
  Met50Pdf.plotOn(xframeMet50,Components(Bern50),LineColor(kRed),LineStyle(kDashed));
  Met50Pdf.plotOn(xframeMet50,LineColor(kBlue));
  //Met50Pdf.plotOn(xframeMet50,Components(sigMet50),LineColor(kBlue));
  xframeMet50->SetMinimum(0);
  xframeMet50->Draw();
  line125->DrawLine(125.3,0,125.3,xframeMet50->GetMaximum());
  double chi2Met50=xframeMet50->chiSquare();
  //cout<<"MET50 Chi2: "<<chi2Met50<<endl;
  THStack* StackMet50 = new THStack("StackMet50","");
  StackMet50->Add(TTHggMET50);//->Draw("SAME");
  StackMet50->Add(VBFHggMET50);//->Draw("SAME");
  StackMet50->Add(WZHggMET50);//->Draw("SAME"); 
  StackMet50->Add(ggHggMET50);//->Draw("SAME");
  StackMet50->Draw("histoSAME");
  xframeMet50->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET50cut.png");
  ggHggMET50->Add(WZHggMET50);ggHggMET50->Add(VBFHggMET50);ggHggMET50->Add(TTHggMET50);
  sigSM = ggHggMET50->IntegralAndError(0,999,sigSMErr,"");
  //sigSM/=10;sigSMErr/=10;
  sigYield=sigMet50Yield.getVal();
  sigYieldErr=sigMet50Yield.getError();
  //cout<<"sigSM: "<<sigSM<<"  sigSMErr: "<<sigSMErr<<endl;
  //cout<<"sigYield: "<<sigYield<<"  sigYieldErr: "<<sigYieldErr<<endl;
  sigStrengthVsMet->SetBinContent(6,sigYield/sigSM);
  sigStrengthVsMet->SetBinError(6,sqrt(sigYieldErr*sigYieldErr/(sigSM*sigSM) + (sigYield*sigYield*sigSMErr*sigSMErr)/(sigSM*sigSM)));


  RooRealVar xMet50_JetReq("xMet50_JetReq","gg_JetReq Invariant Mass with MET>50",90,180,"GeV");
  RooDataHist dataMet50_JetReq("dataMet50_JetReq","dataset",xMet50_JetReq,ggInvarMassMET50_JetReq);
  //Gaussian for peak
  RooRealVar gmeanMet50_JetReq("gmeanMet50_JetReq","gmean",125.3,122,128);
  RooRealVar gsigmaMet50_JetReq("gsigmaMet50_JetReq","gsigma",1,0.1,5);
  //RooGaussian sigMet50("sigMet50","gauss",xMet50,gmeanMet50,gsigmaMet50);
  //Crystal Ball for signal
  RooRealVar cbmeanMet50_JetReq("cb_mean", "cbmean" , 125.3, 122., 128.) ;
  RooRealVar cbsigmaMet50_JetReq("cb_sigma", "cbsigma",2., 0.0, 4.) ;
  RooRealVar nMet50_JetReq("n","n",5.,-5.,20.);
  RooRealVar alphaMet50_JetReq("alpha","alpha",2.5,1.0,4.0);
  RooCBShape sigMet50_JetReq("sigMet50_JetReq", "crystal ball", xMet50_JetReq, cbmeanMet50_JetReq, cbsigmaMet50_JetReq, alphaMet50_JetReq, nMet50_JetReq);
  RooRealVar sigMet50Yield_JetReq("signal yield","signal yield",20,0,500);
  RooRealVar Bern501_JetReq("Bern501","Berstein 1",12.,8.,16.);
  RooRealVar Bern502_JetReq("Bern502","Berstein 2",2.,0.,5.);
  RooRealVar Bern503_JetReq("Bern503","Berstein 3",3.,1.,5.);
  RooRealVar Bern504_JetReq("Bern504","Berstein 4",1.5,0.,3.);
  RooRealVar Bern505_JetReq("Bern505","Berstein 5",5.,0.,10.);
  RooBernstein Bern50_JetReq("Bern50_JetReq","4th order Bernstein Polynomial",xMet50_JetReq,RooArgList(Bern501_JetReq,Bern502_JetReq,Bern503_JetReq,Bern504_JetReq/*,Bern505_JetReq*/));
  RooRealVar Bern50Yield_JetReq("bkgd yield","bkgd yield",500,0,100000);
  RooAddPdf Met50Pdf_JetReq("Met50Pdf_JetReq","Met50Pdf_JetReq",RooArgList(Bern50_JetReq,sigMet50_JetReq),RooArgList(Bern50Yield_JetReq,sigMet50Yield_JetReq));
  RooFitResult *r50_JetReq = Met50Pdf_JetReq.fitTo(dataMet50_JetReq,Extended(kTRUE),Save());
  RooPlot *xframeMet50_JetReq = xMet50_JetReq.frame(Title("Crystal Ball Signal, 4th Order Bernstein Polynomial Background, gg_JetReq MET>50"));
  Met50Pdf_JetReq.paramOn(xframeMet50_JetReq, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.9) );
  dataMet50_JetReq.plotOn(xframeMet50_JetReq,LineColor(kBlack));
  //Met50Pdf_JetReq.plotOn(xframeMet50_JetReq,Components(Bern50_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r50_JetReq,3,kFALSE),FillColor(kViolet));
  //Met50Pdf_JetReq.plotOn(xframeMet50_JetReq,Components(Bern50_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r50_JetReq,2,kFALSE),FillColor(kGreen));
  //Met50Pdf_JetReq.plotOn(xframeMet50_JetReq,Components(Bern50_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r50_JetReq,1,kFALSE),FillColor(kOrange)); 
  dataMet50_JetReq.plotOn(xframeMet50_JetReq,LineColor(kBlack));
  Met50Pdf_JetReq.plotOn(xframeMet50_JetReq,Components(Bern50_JetReq),LineColor(kRed),LineStyle(kDashed));
  Met50Pdf_JetReq.plotOn(xframeMet50_JetReq,LineColor(kBlue));
  //Met50Pdf_JetReq.plotOn(xframeMet50_JetReq,Components(sigMet50_JetReq),LineColor(kBlue));
  xframeMet50_JetReq->SetMinimum(0);
  xframeMet50_JetReq->Draw();
  line125->DrawLine(125.3,0,125.3,xframeMet50_JetReq->GetMaximum());
  double chi2Met50_JetReq=xframeMet50_JetReq->chiSquare();
  //cout<<"MET50_JetReq Chi2: "<<chi2Met50_JetReq<<endl;
  THStack* StackMet50_JetReq = new THStack("StackMet50_JetReq","");
  /* StackMet50_JetReq->Add(TTHggMET50);//->Draw("SAME");
     StackMet50_JetReq->Add(VBFHggMET50);//->Draw("SAME");
     StackMet50_JetReq->Add(WZHggMET50);//->Draw("SAME"); 
     StackMet50_JetReq->Add(ggHggMET50);//->Draw("SAME");
     StackMet50_JetReq->Draw("histoSAME");*/
  xframeMet50_JetReq->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET50cut_JetReq.png");


  RooRealVar xMet60("xMet60","gg Invariant Mass with MET>60",100,164,"GeV");
  RooDataHist dataMet60("dataMet60","dataset",xMet60,ggInvarMassMET60);
  //Gaussian for peak
  RooRealVar gmeanMet60("gmeanMet60","gmean",125.3,122,128);
  RooRealVar gsigmaMet60("gsigmaMet60","gsigma",1,0.1,5);
  //RooGaussian sigMet60("sigMet60","gauss",xMet60,gmeanMet60,gsigmaMet60);
  //Crystal Ball for signal
  RooRealVar cbmeanMet60("cb_mean", "cbmean" , 125.3, 121., 130.) ;
  RooRealVar cbsigmaMet60("cb_sigma", "cbsigma" , 1.2/*, 0.0, 2.2*/) ;
  RooRealVar nMet60("n","n", 9.,0.,35.);
  RooRealVar alphaMet60("alpha","alpha",2.5,1.,4.);
  RooCBShape sigMet60("sigMet50", "crystal ball", xMet60, cbmeanMet60, cbsigmaMet60, alphaMet60, nMet60);
  RooRealVar sigMet60Yield("signal yield","signal yield",20,0,500);
  RooRealVar Bern601("Bern601","Berstein 1",9.,4.,20.);
  RooRealVar Bern602("Bern602","Berstein 2",8.,1.,12.);
  RooRealVar Bern603("Bern603","Berstein 3",6.,1.,8.);
  RooRealVar Bern604("Bern604","Berstein 4",2.,0.,4.);
  RooRealVar Bern605("Bern605","Berstein 5",5.,0.,10.);
  RooBernstein Bern60("Bern60","4th order Bernstein Polynomial",xMet60,RooArgList(Bern601,Bern602,Bern603,Bern604/*,Bern605*/));
  RooRealVar Bern60Yield("bkgd yield","bkgd yield",200,0,10000);
  RooAddPdf Met60Pdf("Met60Pdf","Met60Pdf",RooArgList(Bern60,sigMet60),RooArgList(Bern60Yield,sigMet60Yield));
  RooFitResult *r60 = Met60Pdf.fitTo(dataMet60,Extended(kTRUE),Save());
  RooPlot *xframeMet60 = xMet60.frame(Title("Crystal Ball Signal, 4th Order Bernstein Polynomial Background, gg MET>60"));
  Met60Pdf.paramOn(xframeMet60, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.9) );
  dataMet60.plotOn(xframeMet60,LineColor(kBlack));
  //Met60Pdf.plotOn(xframeMet60,Components(Bern60),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r60,3,kFALSE),FillColor(kViolet));
  //Met60Pdf.plotOn(xframeMet60,Components(Bern60),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r60,2,kFALSE),FillColor(kGreen));
  //Met60Pdf.plotOn(xframeMet60,Components(Bern60),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r60,1,kFALSE),FillColor(kOrange));
  dataMet60.plotOn(xframeMet60,LineColor(kBlack));
  Met60Pdf.plotOn(xframeMet60,Components(Bern60),LineColor(kRed),LineStyle(kDashed));
  Met60Pdf.plotOn(xframeMet60,LineColor(kBlue));
  Met60Pdf.plotOn(xframeMet60,Components(sigMet60),LineColor(kBlue));
  xframeMet60->SetMinimum(0);
  xframeMet60->Draw();
  line125->DrawLine(125.3,0,125.3,xframeMet60->GetMaximum());
  double chi2Met60=xframeMet60->chiSquare();
  //cout<<"MET60 Chi2: "<<chi2Met60<<endl;
  THStack* StackMet60 = new THStack("StackMet60","");
  StackMet60->Add(TTHggMET60);//->Draw("SAME");
  StackMet60->Add(VBFHggMET60);//->Draw("SAME");
  StackMet60->Add(WZHggMET60);//->Draw("SAME"); 
  StackMet60->Add(ggHggMET60);//->Draw("SAME");
  StackMet60->Draw("histoSAME");
  xframeMet60->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET60cut.png");
  ggHggMET60->Add(WZHggMET60);ggHggMET60->Add(VBFHggMET60);ggHggMET60->Add(TTHggMET60);
  sigSM = ggHggMET60->IntegralAndError(0,999,sigSMErr,"");
  //sigSM/=10;sigSMErr/=10;
  sigYield=sigMet60Yield.getVal();
  sigYieldErr=sigMet60Yield.getError();
  //cout<<"sigSM: "<<sigSM<<"  sigSMErr: "<<sigSMErr<<endl;
  //cout<<"sigYield: "<<sigYield<<"  sigYieldErr: "<<sigYieldErr<<endl;
  sigStrengthVsMet->SetBinContent(7,sigYield/sigSM);
  sigStrengthVsMet->SetBinError(7,sqrt(sigYieldErr*sigYieldErr/(sigSM*sigSM) + (sigYield*sigYield*sigSMErr*sigSMErr)/(sigSM*sigSM)));

  RooRealVar xMet60_JetReq("xMet60_JetReq","gg_JetReq Invariant Mass with MET>60",90,180,"GeV");
  RooDataHist dataMet60_JetReq("dataMet60_JetReq","dataset",xMet60_JetReq,ggInvarMassMET60_JetReq);
  //Gaussian for peak
  RooRealVar gmeanMet60_JetReq("gmeanMet60_JetReq","gmean",125.3,122,128);
  RooRealVar gsigmaMet60_JetReq("gsigmaMet60_JetReq","gsigma",1,0.1,5);
  //RooGaussian sigMet60("sigMet60","gauss",xMet60,gmeanMet60,gsigmaMet60);
  //Crystal Ball for signal
  RooRealVar cbmeanMet60_JetReq("cb_mean", "cbmean" , 125.3, 122., 128.) ;
  RooRealVar cbsigmaMet60_JetReq("cb_sigma", "cbsigma",1/*2., 0.0, 4.*/) ;
  RooRealVar nMet60_JetReq("n","n", 4.5/*5.,-5.,20.*/);
  RooRealVar alphaMet60_JetReq("alpha","alpha",2.5,1.0,4.0);
  RooCBShape sigMet60_JetReq("sigMet60_JetReq", "crystal ball", xMet60_JetReq, cbmeanMet60_JetReq, cbsigmaMet60_JetReq, alphaMet60_JetReq, nMet60_JetReq);
  RooRealVar sigMet60Yield_JetReq("signal yield","signal yield",20,0,600);
  RooRealVar Bern601_JetReq("Bern601","Berstein 1",12.,8.,16.);
  RooRealVar Bern602_JetReq("Bern602","Berstein 2",2.,0.,5.);
  RooRealVar Bern603_JetReq("Bern603","Berstein 3",3.,1.,5.);
  RooRealVar Bern604_JetReq("Bern604","Berstein 4",1.5,0.,3.);
  RooRealVar Bern605_JetReq("Bern605","Berstein 5",5.,0.,10.);
  RooBernstein Bern60_JetReq("Bern60_JetReq","4th order Bernstein Polynomial",xMet60_JetReq,RooArgList(Bern601_JetReq,Bern602_JetReq,Bern603_JetReq,Bern604_JetReq/*,Bern605_JetReq*/));
  RooRealVar Bern60Yield_JetReq("bkgd yield","bkgd yield",500,0,100000);
  RooAddPdf Met60Pdf_JetReq("Met60Pdf_JetReq","Met60Pdf_JetReq",RooArgList(Bern60_JetReq,sigMet60_JetReq),RooArgList(Bern60Yield_JetReq,sigMet60Yield_JetReq));
  RooFitResult *r60_JetReq = Met60Pdf_JetReq.fitTo(dataMet60_JetReq,Extended(kTRUE),Save());
  RooPlot *xframeMet60_JetReq = xMet60_JetReq.frame(Title("Crystal Ball Signal, 4th Order Bernstein Polynomial Background, gg_JetReq MET>60"));
  Met60Pdf_JetReq.paramOn(xframeMet60_JetReq, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.9) );
  dataMet60_JetReq.plotOn(xframeMet60_JetReq,LineColor(kBlack));
  //Met60Pdf_JetReq.plotOn(xframeMet60_JetReq,Components(Bern60_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r60_JetReq,3,kFALSE),FillColor(kViolet));
  //Met60Pdf_JetReq.plotOn(xframeMet60_JetReq,Components(Bern60_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r60_JetReq,2,kFALSE),FillColor(kGreen));
  //Met60Pdf_JetReq.plotOn(xframeMet60_JetReq,Components(Bern60_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r60_JetReq,1,kFALSE),FillColor(kOrange)); 
  dataMet60_JetReq.plotOn(xframeMet60_JetReq,LineColor(kBlack));
  Met60Pdf_JetReq.plotOn(xframeMet60_JetReq,Components(Bern60_JetReq),LineColor(kRed),LineStyle(kDashed));
  Met60Pdf_JetReq.plotOn(xframeMet60_JetReq,LineColor(kBlue));
  //Met60Pdf_JetReq.plotOn(xframeMet60_JetReq,Components(sigMet60_JetReq),LineColor(kBlue));
  xframeMet60_JetReq->SetMinimum(0);
  xframeMet60_JetReq->Draw();
  line125->DrawLine(125.3,0,125.3,xframeMet60_JetReq->GetMaximum());
  double chi2Met60_JetReq=xframeMet60_JetReq->chiSquare();
  //cout<<"MET60_JetReq Chi2: "<<chi2Met60_JetReq<<endl;
  THStack* StackMet60_JetReq = new THStack("StackMet60_JetReq","");
  /* StackMet60_JetReq->Add(TTHggMET60);//->Draw("SAME");
     StackMet60_JetReq->Add(VBFHggMET60);//->Draw("SAME");
     StackMet60_JetReq->Add(WZHggMET60);//->Draw("SAME"); 
     StackMet60_JetReq->Add(ggHggMET60);//->Draw("SAME");
     StackMet60_JetReq->Draw("histoSAME");*/
  xframeMet60_JetReq->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET60cut_JetReq.png");

  RooRealVar xMet70("xMet70","gg Invariant Mass with MET>70",90,180,"GeV");
  RooDataHist dataMet70("dataMet70","dataset",xMet70,ggInvarMassMET70);
  //Gaussian for peak
  RooRealVar gmeanMet70("gmeanMet70","gmean",125.3,122,128);
  RooRealVar gsigmaMet70("gsigmaMet70","gsigma",1,0.1,5);
  //RooGaussian sigMet70("sigMet70","gauss",xMet70,gmeanMet70,gsigmaMet70);
  //Crystal Ball for signal
  RooRealVar cbmeanMet70("cb_mean", "cbmean" , 124.3, 120., 125.) ;
  RooRealVar cbsigmaMet70("cb_sigma", "cbsigma" , .9/*1.0, 0.0, 2.2*/) ;
  RooRealVar nMet70("n","n", 5.,1.,25.);
  RooRealVar alphaMet70("alpha","alpha",2.0,0.0,4.);
  RooCBShape sigMet70("sigMet70", "crystal ball", xMet70, cbmeanMet70, cbsigmaMet70, alphaMet70, nMet70);
  RooRealVar sigMet70Yield("signal yield","signal yield",20,0,500);
  RooRealVar Bern701("Bern701","Berstein 1",12.,4.,20.);
  RooRealVar Bern702("Bern702","Berstein 2",6.,-2.,15.);
  RooRealVar Bern703("Bern703","Berstein 3",3.,-3.,10.);
  RooRealVar Bern704("Bern704","Berstein 4",2.,-2.,6.);
  RooRealVar Bern705("Bern705","Berstein 5",5.,0.,10.);
  RooBernstein Bern70("Bern70","4th order Bernstein Polynomial",xMet70,RooArgList(Bern701,Bern702,Bern703,Bern704/*,Bern705*/));
  RooRealVar Bern70Yield("bkgd yield","bkgd yield",100,0,10000);
  RooAddPdf Met70Pdf("Met70Pdf","Met70Pdf",RooArgList(Bern70,sigMet70),RooArgList(Bern70Yield,sigMet70Yield));
  RooFitResult *r70 = Met70Pdf.fitTo(dataMet70,Extended(kTRUE));
  //RooFitResult *r70bg = Bern70.fitTo(dataMet70);
  //Bern70.fitTo(dataMet70,Extended(kTRUE));
  RooPlot *xframeMet70 = xMet70.frame(Title("Crystal Ball Signal, 4th Order Bernstein Polynomial Background, gg MET>70"));
  Met70Pdf.paramOn(xframeMet70, Format("NE",AutoPrecision(1)),Layout(0.68,0.995,0.9) );
  //Bern70.paramOn(xframeMet70, Format("NE",AutoPrecision(1)),Layout(0.68,0.995,0.9) );
  dataMet70.plotOn(xframeMet70,LineColor(kBlack));
  //Met70Pdf.plotOn(xframeMet70_JetReq,Components(Bern70),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r70,3,kFALSE),FillColor(kViolet));
  //Met70Pdf.plotOn(xframeMet70,Components(Bern70),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r70,2,kFALSE),FillColor(kGreen));
  //Met70Pdf.plotOn(xframeMet70,Components(Bern70),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r70,1,kFALSE),FillColor(kOrange)); 
  dataMet70.plotOn(xframeMet70,LineColor(kBlack));
  Met70Pdf.plotOn(xframeMet70,Components(Bern70),LineColor(kRed),LineStyle(kDashed));
  //Met70Pdf.plotOn(xframeMet70,LineColor(kBlue));
  //Bern70.plotOn(xframeMet70,LineColor(kBlue));
  //Met70Pdf.plotOn(xframeMet70,Components(sigMet70),LineColor(kBlue));
  xframeMet70->Draw();
  line125->DrawLine(125.3,0,125.3,xframeMet70->GetMaximum());
  double chi2Met70=xframeMet70->chiSquare();
  //cout<<"MET70 Chi2: "<<chi2Met70<<endl;
  THStack* StackMet70 = new THStack("StackMet70","");
  StackMet70->Add(TTHggMET70);//->Draw("SAME");
  StackMet70->Add(VBFHggMET70);//->Draw("SAME");
  StackMet70->Add(WZHggMET70);//->Draw("SAME"); 
  StackMet70->Add(ggHggMET70);//->Draw("SAME");
  StackMet70->Draw("histoSAME");
  xframeMet70->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET70cut.png");
  ggHggMET70->Add(WZHggMET70);ggHggMET70->Add(VBFHggMET70);ggHggMET70->Add(TTHggMET70);
  sigSM = ggHggMET70->IntegralAndError(0,999,sigSMErr,"");
  //sigSM/=10;sigSMErr/=10;
  sigYield=sigMet70Yield.getVal();
  sigYieldErr=sigMet70Yield.getError();
  //cout<<"sigSM: "<<sigSM<<"  sigSMErr: "<<sigSMErr<<endl;
  //cout<<"sigYield: "<<sigYield<<"  sigYieldErr: "<<sigYieldErr<<endl;
  sigStrengthVsMet->SetBinContent(8,sigYield/sigSM);
  sigStrengthVsMet->SetBinError(8,sqrt(sigYieldErr*sigYieldErr/(sigSM*sigSM) + (sigYield*sigYield*sigSMErr*sigSMErr)/(sigSM*sigSM)));

  RooRealVar xMet70_JetReq("xMet70_JetReq","gg_JetReq Invariant Mass with MET>70",90,180,"GeV");
  RooDataHist dataMet70_JetReq("dataMet70_JetReq","dataset",xMet70_JetReq,ggInvarMassMET70_JetReq);
  //Gaussian for peak
  RooRealVar gmeanMet70_JetReq("gmeanMet70_JetReq","gmean",125.3,122,128);
  RooRealVar gsigmaMet70_JetReq("gsigmaMet70_JetReq","gsigma",1,0.1,5);
  //RooGaussian sigMet70("sigMet70","gauss",xMet70,gmeanMet70,gsigmaMet70);
  //Crystal Ball for signal
  RooRealVar cbmeanMet70_JetReq("cb_mean", "cbmean" , 125.3, 122., 128.) ;
  RooRealVar cbsigmaMet70_JetReq("cb_sigma", "cbsigma",1/*2., 0.0, 4.*/) ;
  RooRealVar nMet70_JetReq("n","n", 4.5/*5.,-5.,20.*/);
  RooRealVar alphaMet70_JetReq("alpha","alpha",2.5,1.0,4.0);
  RooCBShape sigMet70_JetReq("sigMet70_JetReq", "crystal ball", xMet70_JetReq, cbmeanMet70_JetReq, cbsigmaMet70_JetReq, alphaMet70_JetReq, nMet70_JetReq);
  RooRealVar sigMet70Yield_JetReq("signal yield","signal yield",20,0,600);
  RooRealVar Bern701_JetReq("Bern701","Berstein 1",12.,8.,16.);
  RooRealVar Bern702_JetReq("Bern702","Berstein 2",2.,0.,5.);
  RooRealVar Bern703_JetReq("Bern703","Berstein 3",3.,1.,5.);
  RooRealVar Bern704_JetReq("Bern704","Berstein 4",1.5,0.,3.);
  RooRealVar Bern705_JetReq("Bern705","Berstein 5",5.,0.,10.);
  RooBernstein Bern70_JetReq("Bern70_JetReq","4th order Bernstein Polynomial",xMet70_JetReq,RooArgList(Bern701_JetReq,Bern702_JetReq,Bern703_JetReq,Bern704_JetReq/*,Bern705_JetReq*/));
  RooRealVar Bern70Yield_JetReq("bkgd yield","bkgd yield",500,0,100000);
  RooAddPdf Met70Pdf_JetReq("Met70Pdf_JetReq","Met70Pdf_JetReq",RooArgList(Bern70_JetReq,sigMet70_JetReq),RooArgList(Bern70Yield_JetReq,sigMet70Yield_JetReq));
  RooFitResult *r70_JetReq = Met70Pdf_JetReq.fitTo(dataMet70_JetReq,Extended(kTRUE),Save());
  RooPlot *xframeMet70_JetReq = xMet70_JetReq.frame(Title("Crystal Ball Signal, 4th Order Bernstein Polynomial Background, gg_JetReq MET>70"));
  Met70Pdf_JetReq.paramOn(xframeMet70_JetReq, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.9) );
  dataMet70_JetReq.plotOn(xframeMet70_JetReq,LineColor(kBlack));
  //Met70Pdf_JetReq.plotOn(xframeMet70_JetReq,Components(Bern70_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r70_JetReq,3,kFALSE),FillColor(kViolet));
  //Met70Pdf_JetReq.plotOn(xframeMet70_JetReq,Components(Bern70_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r70_JetReq,2,kFALSE),FillColor(kGreen));
  //Met70Pdf_JetReq.plotOn(xframeMet70_JetReq,Components(Bern70_JetReq),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r70_JetReq,1,kFALSE),FillColor(kOrange)); 
  dataMet70_JetReq.plotOn(xframeMet70_JetReq,LineColor(kBlack));
  Met70Pdf_JetReq.plotOn(xframeMet70_JetReq,Components(Bern70_JetReq),LineColor(kRed),LineStyle(kDashed));
  Met70Pdf_JetReq.plotOn(xframeMet70_JetReq,LineColor(kBlue));
  //Met70Pdf_JetReq.plotOn(xframeMet70_JetReq,Components(sigMet70_JetReq),LineColor(kBlue));
  xframeMet70_JetReq->SetMinimum(0);
  xframeMet70_JetReq->Draw();
  line125->DrawLine(125.3,0,125.3,xframeMet70_JetReq->GetMaximum());
  double chi2Met70_JetReq=xframeMet70_JetReq->chiSquare();
  //cout<<"MET70_JetReq Chi2: "<<chi2Met70_JetReq<<endl;
  THStack* StackMet70_JetReq = new THStack("StackMet70_JetReq","");
  /* StackMet70_JetReq->Add(TTHggMET70);//->Draw("SAME");
     StackMet70_JetReq->Add(VBFHggMET70);//->Draw("SAME");
     StackMet70_JetReq->Add(WZHggMET70);//->Draw("SAME"); 
     StackMet70_JetReq->Add(ggHggMET70);//->Draw("SAME");
     StackMet70_JetReq->Draw("histoSAME");*/
  xframeMet70_JetReq->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET70cut_JetReq.png");

  RooRealVar xMet80("xMet80","gg Invariant Mass with MET>80",90,180,"GeV");
  RooDataHist dataMet80("dataMet80","dataset",xMet80,ggInvarMassMET80);
  //Gaussian for peak
  RooRealVar gmeanMet80("gmeanMet80","gmean",125.3,122,128);
  RooRealVar gsigmaMet80("gsigmaMet80","gsigma",1,0.1,5);
  //RooGaussian sigMet80("sigMet80","gauss",xMet80,gmeanMet80,gsigmaMet80);
  //Crystal Ball for signal
  RooRealVar cbmeanMet80("cb_mean", "cbmean" , 124.3, 120., 125.) ;
  RooRealVar cbsigmaMet80("cb_sigma", "cbsigma" , .9/*1.0, 0.0, 2.2*/) ;
  RooRealVar nMet80("n","n", 5.,1.,25.);
  RooRealVar alphaMet80("alpha","alpha",2.0,0.0,4.);
  RooCBShape sigMet80("sigMet80", "crystal ball", xMet80, cbmeanMet80, cbsigmaMet80, alphaMet80, nMet80);
  RooRealVar sigMet80Yield("signal yield","signal yield",20,0,500);
  RooRealVar Bern801("Bern801","Berstein 1",12.,4.,20.);
  RooRealVar Bern802("Bern802","Berstein 2",6.,-2.,15.);
  RooRealVar Bern803("Bern803","Berstein 3",3.,-3.,10.);
  RooRealVar Bern804("Bern804","Berstein 4",2.,-2.,6.);
  RooRealVar Bern805("Bern805","Berstein 5",5.,0.,10.);
  RooBernstein Bern80("Bern80","4th order Bernstein Polynomial",xMet80,RooArgList(Bern801,Bern802,Bern803,Bern804/*,Bern805*/));
  RooRealVar Bern80Yield("bkgd yield","bkgd yield",100,0,10000);
  RooAddPdf Met80Pdf("Met80Pdf","Met80Pdf",RooArgList(Bern80,sigMet80),RooArgList(Bern80Yield,sigMet80Yield));
  RooFitResult *r80 = Met80Pdf.fitTo(dataMet80,Extended(kTRUE));
  //RooFitResult *r80bg = Bern80.fitTo(dataMet80);
  //Bern80.fitTo(dataMet80,Extended(kTRUE));
  RooPlot *xframeMet80 = xMet80.frame(Title("Crystal Ball Signal, 4th Order Bernstein Polynomial Background, gg MET>80"));
  Met80Pdf.paramOn(xframeMet80, Format("NE",AutoPrecision(1)),Layout(0.68,0.995,0.9) );
  //Bern80.paramOn(xframeMet80, Format("NE",AutoPrecision(1)),Layout(0.68,0.995,0.9) );
  dataMet80.plotOn(xframeMet80,LineColor(kBlack));
  //Met80Pdf.plotOn(xframeMet80_JetReq,Components(Bern80),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r80,3,kFALSE),FillColor(kViolet));
  //Met80Pdf.plotOn(xframeMet80,Components(Bern80),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r80,2,kFALSE),FillColor(kGreen));
  //Met80Pdf.plotOn(xframeMet80,Components(Bern80),LineColor(kRed),LineStyle(kDashed),VisualizeError(*r80,1,kFALSE),FillColor(kOrange)); 
  dataMet80.plotOn(xframeMet80,LineColor(kBlack));
  Met80Pdf.plotOn(xframeMet80,Components(Bern80),LineColor(kRed),LineStyle(kDashed));
  //Met80Pdf.plotOn(xframeMet80,LineColor(kBlue));
  //Bern80.plotOn(xframeMet80,LineColor(kBlue));
  //Met80Pdf.plotOn(xframeMet80,Components(sigMet80),LineColor(kBlue));
  xframeMet80->Draw();
  line125->DrawLine(125.3,0,125.3,xframeMet80->GetMaximum());
  double chi2Met80=xframeMet80->chiSquare();
  //cout<<"MET80 Chi2: "<<chi2Met80<<endl;
  THStack* StackMet80 = new THStack("StackMet80","");
  StackMet80->Add(TTHggMET80);//->Draw("SAME");
  StackMet80->Add(VBFHggMET80);//->Draw("SAME");
  StackMet80->Add(WZHggMET80);//->Draw("SAME"); 
  StackMet80->Add(ggHggMET80);//->Draw("SAME");
  StackMet80->Draw("histoSAME");
  xframeMet80->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET80cut.png");
  ggHggMET80->Add(WZHggMET80);ggHggMET80->Add(VBFHggMET80);ggHggMET80->Add(TTHggMET80);
  sigSM = ggHggMET80->IntegralAndError(0,999,sigSMErr,"");
  //sigSM/=10;sigSMErr/=10;
  sigYield=sigMet80Yield.getVal();
  sigYieldErr=sigMet80Yield.getError();
  //cout<<"sigSM: "<<sigSM<<"  sigSMErr: "<<sigSMErr<<endl;
  //cout<<"sigYield: "<<sigYield<<"  sigYieldErr: "<<sigYieldErr<<endl;
  sigStrengthVsMet->SetBinContent(9,sigYield/sigSM);
  sigStrengthVsMet->SetBinError(9,sqrt(sigYieldErr*sigYieldErr/(sigSM*sigSM) + (sigYield*sigYield*sigSMErr*sigSMErr)/(sigSM*sigSM)));

  RooRealVar xMet100("xMet100","gg Invariant Mass with MET>100",90,180,"GeV");
  RooDataHist dataMet100("dataMet100","dataset",xMet100,ggInvarMassMET100);
  //Gaussian for peak
  RooRealVar gmeanMet100("gmeanMet100","gmean",125.3,122,128);
  RooRealVar gsigmaMet100("gsigmaMet100","gsigma",1,0.1,5);
  //RooGaussian sigMet100("sigMet100","gauss",xMet100,gmeanMet100,gsigmaMet100);
  //Crystal Ball for signal
  RooRealVar cbmeanMet100("cb_mean", "cbmean" , 124.3, 120., 125.) ;
  RooRealVar cbsigmaMet100("cb_sigma", "cbsigma" , .9/*1.0, 0.0, 2.2*/) ;
  RooRealVar nMet100("n","n", 5.,1.,25.);
  RooRealVar alphaMet100("alpha","alpha",2.0,0.0,4.);
  RooCBShape sigMet100("sigMet50", "crystal ball", xMet100, cbmeanMet100, cbsigmaMet100, alphaMet100, nMet100);
  RooRealVar sigMet100Yield("signal yield","signal yield",20,0,500);
  RooRealVar Bern1001("Bern1001","Berstein 1",12.,4.,20.);
  RooRealVar Bern1002("Bern1002","Berstein 2",6.,-2.,15.);
  RooRealVar Bern1003("Bern1003","Berstein 3",3.,-3.,10.);
  RooRealVar Bern1004("Bern1004","Berstein 4",2.,-2.,6.);
  RooRealVar Bern1005("Bern1005","Berstein 5",5.,0.,10.);
  RooBernstein Bern100("Bern100","4th order Bernstein Polynomial",xMet100,RooArgList(Bern1001,Bern1002,Bern1003,Bern1004/*,Bern1005*/));
  RooRealVar Bern100Yield("bkgd yield","bkgd yield",100,0,10000);
  RooAddPdf Met100Pdf("Met100Pdf","Met100Pdf",RooArgList(Bern100,sigMet100),RooArgList(Bern100Yield,sigMet100Yield));
  Met100Pdf.fitTo(dataMet100,Extended(kTRUE));
  //Bern100.fitTo(dataMet100,Extended(kTRUE));
  RooPlot *xframeMet100 = xMet100.frame(Title("Crystal Ball Signal, 4th Order Bernstein Polynomial Background, gg MET>100"));
  Met100Pdf.paramOn(xframeMet100, Format("NE",AutoPrecision(1)),Layout(0.68,0.995,0.9) );
  //Bern100.paramOn(xframeMet100, Format("NE",AutoPrecision(1)),Layout(0.68,0.995,0.9) );
  dataMet100.plotOn(xframeMet100,LineColor(kBlack));
  Met100Pdf.plotOn(xframeMet100,Components(Bern100),LineColor(kRed),LineStyle(kDashed));
  Met100Pdf.plotOn(xframeMet100,LineColor(kBlue));
  //Bern100.plotOn(xframeMet100,LineColor(kBlue));
  //Met100Pdf.plotOn(xframeMet100,Components(sigMet100),LineColor(kBlue));
  xframeMet100->Draw();
  line125->DrawLine(125.3,0,125.3,xframeMet100->GetMaximum());
  double chi2Met100=xframeMet100->chiSquare();
  //cout<<"MET100 Chi2: "<<chi2Met100<<endl;
  THStack* StackMet100 = new THStack("StackMet100","");
  StackMet100->Add(TTHggMET100);//->Draw("SAME");
  StackMet100->Add(VBFHggMET100);//->Draw("SAME");
  StackMet100->Add(WZHggMET100);//->Draw("SAME"); 
  StackMet100->Add(ggHggMET100);//->Draw("SAME");
  StackMet100->Draw("histoSAME");
  xframeMet100->Draw("SAME");
  legHiggs->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMET100cut.png");
  ggHggMET100->Add(WZHggMET100);ggHggMET100->Add(VBFHggMET100);ggHggMET100->Add(TTHggMET30);
  sigSM = ggHggMET100->IntegralAndError(0,999,sigSMErr,"");
  //sigSM/=10;sigSMErr/=10;
  sigYield=sigMet100Yield.getVal();
  sigYieldErr=sigMet100Yield.getError();
  //cout<<"sigSM: "<<sigSM<<"  sigSMErr: "<<sigSMErr<<endl;
  //cout<<"sigYield: "<<sigYield<<"  sigYieldErr: "<<sigYieldErr<<endl;
  sigStrengthVsMet->SetBinContent(11,sigYield/sigSM);
  sigStrengthVsMet->SetBinError(11,sqrt(sigYieldErr*sigYieldErr/(sigSM*sigSM) + (sigYield*sigYield*sigSMErr*sigSMErr)/(sigSM*sigSM)));
 
  sigStrengthVsMet->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  sigStrengthVsMet->GetYaxis()->SetTitle("Signal Yield / SM yield");  
  sigStrengthVsMet->Fit("pol0","","",30,110);
  sigStrengthVsMet->SetMinimum(0);
  sigStrengthVsMet->SetMaximum(15);
  sigStrengthVsMet->Draw();
  TLine lineSigStrength(0,1,110,1);lineSigStrength.SetLineWidth(3);lineSigStrength.Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/signalStrengthVsMet.png");

  ggInvarMassMET40->SetLineColor(kRed);ggInvarMassMET40->SetMarkerColor(kRed);
  ggInvarMassMET50->SetLineColor(kBlue);ggInvarMassMET50->SetMarkerColor(kBlue);
  ggInvarMassMET60->SetLineColor(kGreen);ggInvarMassMET60->SetMarkerColor(kGreen);
  ggInvarMassMET70->SetLineColor(kCyan);ggInvarMassMET70->SetMarkerColor(kCyan);
  ggInvarMassMET80->SetLineColor(kViolet);ggInvarMassMET80->SetMarkerColor(kViolet);
  ggInvarMassMET100->SetLineColor(kYellow);ggInvarMassMET100->SetMarkerColor(kYellow);

  ggInvarMassMET30->GetXaxis()->SetRangeUser(10,200);
  ggInvarMassMET40->GetXaxis()->SetRangeUser(10,200);
  ggInvarMassMET50->GetXaxis()->SetRangeUser(10,200);
  ggInvarMassMET60->GetXaxis()->SetRangeUser(10,200);
  ggInvarMassMET70->GetXaxis()->SetRangeUser(10,200);
  ggInvarMassMET80->GetXaxis()->SetRangeUser(10,200);
  ggInvarMassMET100->GetXaxis()->SetRangeUser(10,200);

  ggInvarMassMET30->Draw("PE");
  ggInvarMassMET40->Draw("PESAME");
  ggInvarMassMET50->Draw("PESAME");
  ggInvarMassMET60->Draw("PESAME");
  //ggInvarMassMET70->Draw("PESAME");
  //ggInvarMassMET80->Draw("PESAME");
  //ggInvarMassMET100->Draw("PESAME");
  line125->DrawLine(125.3,0,125.3,ggInvarMassMET30->GetMaximum());
  TLegend *legInvMass = new TLegend(.6,.6,.8,.8);
  legInvMass->AddEntry(ggInvarMassMET30,"gg, MET>30");
  legInvMass->AddEntry(ggInvarMassMET40,"gg, MET>40");
  legInvMass->AddEntry(ggInvarMassMET50,"gg, MET>50");
  legInvMass->AddEntry(ggInvarMassMET60,"gg, MET>60");
  legInvMass->SetFillColor(kWhite);
  legInvMass->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMETcut.png");

  ggInvarMassMET30->GetXaxis()->SetRangeUser(100,150);
  ggInvarMassMET40->GetXaxis()->SetRangeUser(100,150);
  ggInvarMassMET50->GetXaxis()->SetRangeUser(100,150);
  ggInvarMassMET60->GetXaxis()->SetRangeUser(100,150);
  ggInvarMassMET70->GetXaxis()->SetRangeUser(100,150);
  ggInvarMassMET80->GetXaxis()->SetRangeUser(100,150);
  ggInvarMassMET100->GetXaxis()->SetRangeUser(100,150);
  ggInvarMassMET30->GetYaxis()->SetRangeUser(0,375);
  ggInvarMassMET40->GetYaxis()->SetRangeUser(0,375);
  ggInvarMassMET50->GetYaxis()->SetRangeUser(0,375);
  ggInvarMassMET60->GetYaxis()->SetRangeUser(0,375);
  ggInvarMassMET70->GetYaxis()->SetRangeUser(0,375);
  ggInvarMassMET80->GetYaxis()->SetRangeUser(0,375);
  ggInvarMassMET100->GetYaxis()->SetRangeUser(0,375);

  ggInvarMassMET30->Draw("PE");
  ggInvarMassMET40->Draw("PESAME");
  ggInvarMassMET50->Draw("PESAME");
  ggInvarMassMET60->Draw("PESAME");
  //ggInvarMassMET70->Draw("PESAME");
  //ggInvarMassMET80->Draw("PESAME");
  //ggInvarMassMET100->Draw("PESAME");
  line125->DrawLine(125.3,0,125.3,378);
  legInvMass->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMETcutZOOM.png");

  float ggMetX = 1./ggInvarMassMET30->Integral();ggInvarMassMET30->Scale(ggMetX);
  ggMetX = 1./ggInvarMassMET40->Integral();ggInvarMassMET40->Scale(ggMetX);
  ggMetX = 1./ggInvarMassMET50->Integral();ggInvarMassMET50->Scale(ggMetX);
  ggMetX = 1./ggInvarMassMET60->Integral();ggInvarMassMET60->Scale(ggMetX);
  ggMetX = 1./ggInvarMassMET70->Integral();ggInvarMassMET70->Scale(ggMetX);
  ggMetX = 1./ggInvarMassMET80->Integral();ggInvarMassMET80->Scale(ggMetX);
  ggMetX = 1./ggInvarMassMET100->Integral();ggInvarMassMET100->Scale(ggMetX);

  ggInvarMassMET30->GetXaxis()->SetRangeUser(10,200);
  ggInvarMassMET40->GetXaxis()->SetRangeUser(10,200);
  ggInvarMassMET50->GetXaxis()->SetRangeUser(10,200);
  ggInvarMassMET60->GetXaxis()->SetRangeUser(10,200);
  ggInvarMassMET70->GetXaxis()->SetRangeUser(10,200);
  ggInvarMassMET80->GetXaxis()->SetRangeUser(10,200);
  ggInvarMassMET100->GetXaxis()->SetRangeUser(10,200);
  ggInvarMassMET30->Draw("PE");
  ggInvarMassMET40->Draw("PESAME");
  ggInvarMassMET50->Draw("PESAME");
  ggInvarMassMET60->Draw("PESAME");
  //ggInvarMassMET70->Draw("PESAME");
  //ggInvarMassMET80->Draw("PESAME");
  //ggInvarMassMET100->Draw("PESAME");
  line125->DrawLine(125.3,0,125.3,ggInvarMassMET30->GetMaximum());
  legInvMass->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMETcut_Scaled.png");

  ggInvarMassMET30->GetXaxis()->SetRangeUser(100,150);
  ggInvarMassMET40->GetXaxis()->SetRangeUser(100,150);
  ggInvarMassMET50->GetXaxis()->SetRangeUser(100,150);
  ggInvarMassMET60->GetXaxis()->SetRangeUser(100,150);
  ggInvarMassMET70->GetXaxis()->SetRangeUser(100,150);
  ggInvarMassMET80->GetXaxis()->SetRangeUser(100,150);
  ggInvarMassMET100->GetXaxis()->SetRangeUser(100,150);

  ggInvarMassMET30->Draw("PE");
  ggInvarMassMET40->Draw("PESAME");
  ggInvarMassMET50->Draw("PESAME");
  ggInvarMassMET60->Draw("PESAME");
  //ggInvarMassMET70->Draw("PESAME");
  //ggInvarMassMET80->Draw("PESAME");
  //ggInvarMassMET100->Draw("PESAME");
  line125->DrawLine(125.3,.0092,125.3,.08);
  legInvMass->Draw("SAME");
  ggInvarMassMET30->Draw("PESAME");
  ggInvarMassMET40->Draw("PESAME");
  ggInvarMassMET50->Draw("PESAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggInvarMassWithMETcut_ScaledZOOM.png");

  c2->cd();
  p1->cd();p1->SetLogy(1);
  TH1F* ggMETInvarMass123_128 = (TH1F*)fin.Get("ggMETInvarMass120_130");
  TH1F* ggMETInvarMass95_110 = (TH1F*)fin.Get("ggMETInvarMass95_110");
  TH1F* ggMETInvarMass135_160 = (TH1F*)fin.Get("ggMETInvarMass135_160");

  ggMETInvarMass123_128->Sumw2();ggMETInvarMass95_110->Sumw2();ggMETInvarMass135_160->Sumw2();

  //ggMETInvarMass123_128->Scale(1./ggMETInvarMass123_128->Integral());
  //ggMETInvarMass95_110->Scale(1./ggMETInvarMass95_110->Integral());
  TH1F*ggMETInvarMass135_160_forBern=(TH1F*)ggMETInvarMass135_160->Clone();
  ggMETInvarMass135_160->Scale(ggMETInvarMass95_110->Integral()/ggMETInvarMass135_160->Integral());

  TH1F* ggMETInvarMassSideBand = (TH1F*)ggMETInvarMass95_110->Clone();
  ggMETInvarMassSideBand->Add(ggMETInvarMass135_160);
  float scaleHiggsTag = (ggMETInvarMass123_128->Integral(1,4)-TTHggMetNew2->Integral(1,4))/ggMETInvarMassSideBand->Integral(1,4);
  ggMETInvarMassSideBand->Scale(scaleHiggsTag);

  //ggMETInvarMass123_128->Rebin(5);ggMETInvarMass95_110->Rebin(5);ggMETInvarMass135_160->Rebin(5);ggMETInvarMassSideBand->Rebin(5);
  //ggMETInvarMass123_128->SetFillStyle(3154);ggMETInvarMass123_128->SetFillColor(kRed);ggMETInvarMass123_128->SetMarkerSize(1);
  //ggMETInvarMass123_128->SetMarkerColor(kRed);ggMETInvarMass123_128->SetLineColor(kRed);
  ggMETInvarMass95_110->SetLineColor(kOrange);ggMETInvarMass135_160->SetLineColor(kGreen);
  ggMETInvarMassSideBand->SetLineColor(kGray+1);ggMETInvarMassSideBand->SetMarkerSize(0);
  ggMETInvarMassSideBand->SetFillColor(kGray+1);
  //ggMETInvarMass123_128->GetXaxis()->SetRangeUser(0,150);
  //ggMETInvarMass123_128->GetYaxis()->SetRangeUser(2e-5,.3);
  ggMETInvarMass123_128->GetYaxis()->SetTitle("Normalized Events");
  ggMETInvarMass123_128->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  ggMETInvarMass123_128->SetTitle("ggMet with invariant mass cuts");
  ggMETInvarMass123_128->SetLineWidth(2);
  ggMETInvarMassSideBand->SetLineWidth(2);
  //c1->SetLogy(1);
  //gStyle->SetErrorX(1);
  //ggMETInvarMassSideBand->GetYaxis()->SetRangeUser(1e-6,.3);
  //ggMETInvarMass123_128->GetXaxis()->SetRangeUser(0,250);
  //  TH1F* newHiggsMettemp = (TH1F*)ggMETInvarMass123_128->Clone();
  //  TH1F* newHiggsMetSBtemp = (TH1F*)ggMETInvarMassSideBand->Clone();
  TH1F* newHiggsMet=(TH1F*)ggMETInvarMass123_128->Rebin(NmetBins,"newHiggsMet",xbins);
  TH1F* newHiggsMetSB=(TH1F*)ggMETInvarMassSideBand->Rebin(NmetBins,"newHiggsMetSB",xbins);
  //newHiggsMet->Rebin(20);newHiggsMetSB->Rebin(20);
  newHiggsMet->Add(newHiggsMetSB,-1);
  TH1F* sig123_128 = (TH1F*)ggMETInvarMass123_128->Rebin(NmetBins,"sig123_128",xbins);
  TH1F* bg95_110_and_135_160 = (TH1F*)ggMETInvarMassSideBand->Rebin(NmetBins,"bg95_110_and_135_160",xbins);
  for(int i=1;i<sig123_128->GetNbinsX()+2;i++){
    float x = sig123_128->GetBinContent(i)/sig123_128->GetBinWidth(i);
    sig123_128->SetBinContent(i,x);
    x = sig123_128->GetBinError(i)/sig123_128->GetBinWidth(i);
    sig123_128->SetBinError(i,x);
    x = bg95_110_and_135_160->GetBinContent(i)/bg95_110_and_135_160->GetBinWidth(i);
    bg95_110_and_135_160->SetBinContent(i,x);
    x = bg95_110_and_135_160->GetBinError(i)/bg95_110_and_135_160->GetBinWidth(i);
    bg95_110_and_135_160->SetBinError(i,x);
  }
  sig123_128->GetYaxis()->SetRangeUser(2e-3,1000);
  gStyle->SetErrorX();
  //sig123_128->GetXaxis()->SetRangeUser(0,79);
  sig123_128->Draw("PE");
  TH1F* ratioLowMet = (TH1F*)bg95_110_and_135_160->Clone();ratioLowMet->Add(TTHggMetNew2);
  ratioLowMet->SetLineColor(kRed);ratioLowMet->SetMarkerColor(kRed);ratioLowMet->SetFillColor(kRed);ratioLowMet->SetFillStyle(3254);//ratioLowMet->SetMarkerSize(0);
  THStack *HiggsScaleLowMET = new THStack("HiggsScaleLowMET","");
  HiggsScaleLowMET->Add(TTHggMetNew2);HiggsScaleLowMET->Add(bg95_110_and_135_160);
  HiggsScaleLowMET->Draw("histoSAMES");
  //ggMETInvarMass95_110->Draw("SAME");
  //ggMETInvarMass135_160->Draw("SAME");
  ////bg95_110_and_135_160->Draw("PESAMES");
  //metstackRebin->Draw("histoSAME");
  ////TTHggMetNew2->Draw("histoSAMES");
  sig123_128->Draw("PEsames");
  TLegend *legMetInvMass = new TLegend(.4,.5,.8,.8);
  legMetInvMass->SetBorderSize(0);
  legMetInvMass->SetFillStyle(0);
  legMetInvMass->AddEntry(ggMETInvarMass123_128,"ggMet, 123<InvMass<128","pl");
  //legMetInvMass->AddEntry(ggMETInvarMass95_110,"ggMet, 95<InvMass<110","l");
  //legMetInvMass->AddEntry(ggMETInvarMass135_160,"ggMet, 135<InvMass<160","l");
  legMetInvMass->AddEntry(ggMETInvarMassSideBand,"#splitline{ggMet, InvMass SideBand}{95<InvMass110 and 135<InvMass<160}","f");
  legMetInvMass->AddEntry(TTHggMetNew2,"#splitline{SM Higgs}{ggHgg+VBFHgg+W/ZHgg+TTHgg}","f");
  legMetInvMass->SetFillColor(kWhite);
  ratioLowMet->Draw("E2SAMES");
  legMetInvMass->Draw("SAME");  
  p2->cd();
  TH1F* sig123_128_Div = (TH1F*)sig123_128->Clone();
  sig123_128_Div->SetLineColor(kRed);sig123_128_Div->SetMarkerColor(kRed);
  //sig123_128_Div->Divide(bg95_110_and_135_160);
  sig123_128_Div->Divide(ratioLowMet);
  sig123_128_Div->GetYaxis()->SetRangeUser(0,1.8);
  sig123_128_Div->SetTitle("");
  sig123_128_Div->GetYaxis()->SetTitle("Data/Prediction");
  sig123_128_Div->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  sig123_128_Div->GetYaxis()->SetTitleOffset(0.48);
  sig123_128_Div->GetYaxis()->SetTitleSize(0.15);
  sig123_128_Div->GetXaxis()->SetTitleSize(0.2);
  sig123_128_Div->GetYaxis()->SetLabelSize(0.12);
  sig123_128_Div->GetXaxis()->SetLabelSize(0.15);
  sig123_128_Div->Draw("PE");
  TLine l1(0,1,250,1);l1.Draw();
  p3->cd();
  TPaveText *met2 = new TPaveText(.65,.5,.85,.99,"NDC");
  met2->AddText("E_{T}^{miss} [GeV]");
  met2->SetFillColor(0);
  met2->SetBorderSize(0);
  met2->SetTextSize(.3);
  met2->Draw();
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggMET_InvarMassCuts.png");
  
  p1->cd();
  //TH1F* ggMETInvarMass95_110_bern2 = (TH1F*)ggMETInvarMass95_110->Clone();
  //TH1F* ggMETInvarMass135_160_bern2 = (TH1F*)ggMETInvarMass135_160->Clone();
  //HiggSigOverLowSB/=ggMETInvarMass95_110_bern->Integral();
  //ggMETInvarMass95_110_bern2->Scale(HiggSigOverLowSB);
  //HiggSigOverHighSB/=ggMETInvarMass135_160_bern->Integral();
  //ggMETInvarMass135_160_bern2->Scale(HiggSigOverHighSB);
  TH1F* ggMETInvarMass95_110_bern = (TH1F*)ggMETInvarMass95_110->Rebin(NmetBins,"ggMETInvarMass95_110_bern",xbins);
  TH1F* ggMETInvarMass135_160_bern = (TH1F*)ggMETInvarMass135_160_forBern->Rebin(NmetBins,"ggMETInvarMass135_160_bern",xbins);
  //float higgslowscale=HiggSigNoNorm/ggMETInvarMass95_110_bern->Integral();
  //float higgshighscale=HiggSigNoNorm/ggMETInvarMass135_160_bern->Integral();
  float higgslowscale = (HiggSigNoNorm/*-SMhiggsYield*/)/HiggLowSBNoNorm;//(HiggSigOverLowSB*ggMETInvarMass95_110_bern->Integral()-SMhiggsYield)/ggMETInvarMass95_110_bern->Integral();
  float higgshighscale = (HiggSigNoNorm/*-SMhiggsYield*/)/HiggHighSBNoNorm;//(HiggSigOverHighSB*ggMETInvarMass135_160_bern->Integral()-SMhiggsYield)/ggMETInvarMass135_160_bern->Integral();
  //ggMETInvarMass95_110_bern->Scale(HiggSigOverLowSB);
  //ggMETInvarMass135_160_bern->Scale(HiggSigOverHighSB);
  ggMETInvarMass95_110_bern->Scale(higgslowscale);
  ggMETInvarMass135_160_bern->Scale(higgshighscale);
  for(int i=1;i<sig123_128->GetNbinsX()+2;i++){
    float x = ggMETInvarMass95_110_bern->GetBinContent(i)/ggMETInvarMass95_110_bern->GetBinWidth(i);
    ggMETInvarMass95_110_bern->SetBinContent(i,x);
    x = ggMETInvarMass95_110_bern->GetBinError(i)/ggMETInvarMass95_110_bern->GetBinWidth(i);
    ggMETInvarMass95_110_bern->SetBinError(i,x);
    x = ggMETInvarMass135_160_bern->GetBinContent(i)/ggMETInvarMass135_160_bern->GetBinWidth(i);
    ggMETInvarMass135_160_bern->SetBinContent(i,x);
    x = ggMETInvarMass135_160_bern->GetBinError(i)/ggMETInvarMass135_160_bern->GetBinWidth(i);
    ggMETInvarMass135_160_bern->SetBinError(i,x);
  }
  TH1F* ggMETInvarMass_comb_bern=(TH1F*)ggMETInvarMass95_110_bern->Clone();ggMETInvarMass_comb_bern->Add(ggMETInvarMass135_160_bern);ggMETInvarMass_comb_bern->Scale(0.5);

  sig123_128->SetMarkerColor(kBlack);sig123_128->SetLineColor(kBlack);
  ggMETInvarMass95_110_bern->SetMarkerColor(kBlue);ggMETInvarMass95_110_bern->SetLineColor(kBlue);ggMETInvarMass95_110_bern->SetFillColor(kBlue);ggMETInvarMass95_110_bern->SetFillStyle(3254);
  ggMETInvarMass135_160_bern->SetMarkerColor(kRed);ggMETInvarMass135_160_bern->SetLineColor(kRed);ggMETInvarMass135_160_bern->SetFillColor(kRed);ggMETInvarMass135_160_bern->SetFillStyle(3254);
  ggMETInvarMass_comb_bern->SetMarkerColor(kGray+1);ggMETInvarMass_comb_bern->SetLineColor(kGray+1);ggMETInvarMass_comb_bern->SetFillColor(kGray+1);//ggMETInvarMass_comb_bern->SetFillStyle(3254);
  ggMETInvarMass95_110_bern->SetMarkerSize(0.75);ggMETInvarMass135_160_bern->SetMarkerSize(0.75);ggMETInvarMass_comb_bern->SetMarkerSize(0.75);
  
  THStack *HiggsBernStack = new THStack("HiggsBernStack","");
  HiggsBernStack->Add(TTHggMetNew2);HiggsBernStack->Add(ggMETInvarMass_comb_bern);
  TH1F* ggMETInvarMass_comb_bern_plus_SMhiggs = (TH1F*)ggMETInvarMass_comb_bern->Clone();ggMETInvarMass_comb_bern_plus_SMhiggs->Add(TTHggMetNew2);
  ggMETInvarMass_comb_bern_plus_SMhiggs->SetMarkerSize(0);ggMETInvarMass_comb_bern_plus_SMhiggs->SetFillColor(kRed);ggMETInvarMass_comb_bern_plus_SMhiggs->SetFillStyle(3254);

  sig123_128->GetXaxis()->SetRangeUser(0,79);
  sig123_128->Draw("PE");
  HiggsBernStack->Draw("histoSAMES");
  ggMETInvarMass_comb_bern_plus_SMhiggs->Draw("E2SAMES");
  sig123_128->Draw("PESAMES");
  //ggMETInvarMass_comb_bern->Draw("E2SAMES");
  //ggMETInvarMass95_110_bern->Draw("E2SAMES");
  //ggMETInvarMass135_160_bern->Draw("E2SAMES");
  //TTHggMetNew2->Draw("histoSAMES");
  //legMetInvMass->Draw("SAME");  
  p2->cd();
  TH1F* sig123_128_Div_low = (TH1F*)sig123_128->Clone();
  TH1F* sig123_128_Div_high = (TH1F*)sig123_128->Clone();
  TH1F* sig123_128_Div_comb = (TH1F*)sig123_128->Clone();
  TH1F* sig123_128_Div_comb_plusSMhiggs = (TH1F*)sig123_128->Clone();
  sig123_128_Div_low->Divide(ggMETInvarMass95_110_bern);
  sig123_128_Div_high->Divide(ggMETInvarMass135_160_bern);
  sig123_128_Div_comb->Divide(ggMETInvarMass_comb_bern);
  sig123_128_Div_comb_plusSMhiggs->Divide(ggMETInvarMass_comb_bern_plus_SMhiggs);
  sig123_128_Div_low->SetLineColor(kBlue);sig123_128_Div_high->SetMarkerColor(kBlue);sig123_128_Div_high->SetMarkerSize(0.75);
  sig123_128_Div_comb->SetLineColor(kViolet);sig123_128_Div_comb->SetMarkerColor(kViolet);sig123_128_Div_comb->SetMarkerSize(0.75);
  sig123_128_Div_comb_plusSMhiggs->SetLineColor(kRed);sig123_128_Div_comb_plusSMhiggs->SetMarkerColor(kRed);//sig123_128_Div_comb_plusSMhiggs->SetMarkerSize(0.75);
  sig123_128_Div_high->SetLineColor(kRed);sig123_128_Div_high->SetMarkerColor(kRed);sig123_128_Div_high->SetMarkerSize(0.75);
  sig123_128_Div_comb_plusSMhiggs->GetYaxis()->SetTitle("");
  sig123_128_Div_comb_plusSMhiggs->SetTitle("");
  sig123_128_Div_comb_plusSMhiggs->GetYaxis()->SetRangeUser(0.,1.8);
  sig123_128_Div_comb_plusSMhiggs->SetTitle("");
  sig123_128_Div_comb_plusSMhiggs->GetYaxis()->SetTitle("Data/Prediction");
  sig123_128_Div_comb_plusSMhiggs->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  sig123_128_Div_comb_plusSMhiggs->GetYaxis()->SetTitleOffset(0.48);
  sig123_128_Div_comb_plusSMhiggs->GetYaxis()->SetTitleSize(0.15);
  sig123_128_Div_comb_plusSMhiggs->GetXaxis()->SetTitleSize(0.2);
  sig123_128_Div_comb_plusSMhiggs->GetYaxis()->SetLabelSize(0.12);
  sig123_128_Div_comb_plusSMhiggs->GetXaxis()->SetLabelSize(0.15);
  sig123_128_Div_comb_plusSMhiggs->Draw("PE");
  //sig123_128_Div_high->Draw("PEsames");
  //sig123_128_Div_low->Draw("PEsames");
  l1.DrawLine(0,1,80,1);
  p3->cd();
  met2->Draw();
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggMET_InvarMassCuts_BERNfit.png");
  

  c1->cd();
  c1->SetLogy(0);
  ggMETInvarMass123_128->Rebin(10);
  TH1F* ggMETInvarMassSideBand_new = (TH1F*)ggMETInvarMassSideBand->Clone();
  ggMETInvarMassSideBand->Rebin(10);
  RooRealVar SigPDF("SigPDF","123-128",20,150,"GeV");
  RooRealVar BgPDF("BgPDF","95-110 and 135-160",20,150,"GeV");
  //RooRealVar SigPDFfit("SigPDFfit","123-128",25,200,"GeV");
  //RooRealVar BgPDFfit("BgPDFfit","95-110 and 135-160",25,200,"GeV");
  RooDataHist dataSigPDF("dataSigPDF","dataSigPDF",SigPDF,ggMETInvarMass123_128);
  RooDataHist dataBgPDF("dataBgPDF","dataBgPDF",BgPDF,ggMETInvarMassSideBand);
  RooRealVar Polsig1("Polsig1","Polsig1",10,0,20.);
  RooRealVar Polsig2("Polsig2","Polsig2",0,-10,10);
  RooRealVar Polsig3("Polsig3","Polsig3",0,-10,10);
  RooRealVar Polsig4("Polsig4","Polsig4",0,-10,10);
  RooRealVar Polsig5("Polsig5","Polsig5",0,-10,10);  
  RooRealVar Polbg1("Polbg1","Polbg1",10,0,20);
  RooRealVar Polbg2("Polbg2","Polbg2",0,-10,10);
  RooRealVar Polbg3("Polbg3","Polbg3",0,-10,10);
  RooRealVar Polbg4("Polbg4","Polbg4",0,-10,10);
  RooRealVar Polbg5("Polbg5","Polbg5",0,-10,10);  
  RooBernstein PolySig("PolySig","4th Order Polynomial",SigPDF,RooArgList(Polsig1,Polsig2,Polsig3,Polsig4/*,Polsig5*/));
  RooBernstein PolyBg("PolyBg","4th Order Polynomial",BgPDF,RooArgList(Polbg1,Polbg2,Polbg3,Polbg4/*,Polbg5*/));
  // RooPolynomial PolySig("PolySig","4th Order Polynomial",SigPDF,RooArgList(Polsig1,Polsig2,Polsig3,Polsig4/*,Polsig5*/));
  // RooPolynomial PolyBg("PolyBg","4th Order Polynomial",BgPDF,RooArgList(Polbg1,Polbg2,Polbg3,Polbg4/*,Polbg5*/));
  RooFitResult *polyFitSig = PolySig.fitTo(dataSigPDF,Extended(kFALSE),Save());
  RooFitResult *polyFitBg = PolyBg.fitTo(dataBgPDF,Extended(kFALSE),Save());
  RooPlot *xframeSig = SigPDF.frame(Title("MET PDFs"));
  RooPlot *xframeBg = BgPDF.frame(Title("MET PDFs"));
  PolySig.paramOn(xframeSig, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.9) );
  PolyBg.paramOn(xframeSig, Format("NE",AutoPrecision(2)),Layout(0.68,0.995,0.6) );
  dataSigPDF.plotOn(xframeSig,LineColor(kRed));
  PolySig.plotOn(xframeSig,LineColor(kRed));
  dataBgPDF.plotOn(xframeBg,LineColor(kBlue));
  PolyBg.plotOn(xframeBg,LineColor(kBlue));
  //xframeSig->SetAxisRange(.001,3000,"Y");
  xframeSig->Draw();
  xframeBg->Draw("SAMES");
  //sig123_128->Draw();
  //bg95_110_and_135_160->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggMET_InvarMassCutsFitModels.png");
  TH1F* sig123_128_Diff = (TH1F*)sig123_128->Clone();
  TH1F* bg95_110_and_135_160_Diff = (TH1F*)bg95_110_and_135_160->Clone();
  //float dif = sig123_128_Diff->Integral(1,4)/bg95_110_and_135_160_Diff->Integral(1,4);
  //bg95_110_and_135_160_Diff->Scale(dif);
  sig123_128_Diff->Add(bg95_110_and_135_160_Diff,-1);
  sig123_128_Diff->SetTitle("InvarMass signal minus background");
  c1->SetLogy(1);
  sig123_128_Diff->GetYaxis()->SetRangeUser(1e-4,11);
  sig123_128_Diff->Draw("PE");
  TTHggMetNew2->Draw("histoSAME");
  TLine l0(0,0,300,0);l0.Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggMET_InvarMassCutsDiff.png");

  c1->SetLogy(0);
  newHiggsMet->GetYaxis()->SetRangeUser(-10,10);
  newHiggsMet->Draw("PE");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/newHiggsMet_before.png");
  //HiggsScale=1./newHiggsMet->Integral();newHiggsMet->Scale(HiggsScale);
  //newHiggsMet->GetYaxis()->SetRangeUser(-10,10);
  //newHiggsMet->Draw("PE");
  //c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/newHiggsMet_after.png");
  c1->SetLogy(1);

  f_MetAndMassShapes.cd();
  float low=ggInvarMass->FindBin(95),high=ggInvarMass->FindBin(160);
  HiggsScale=1./ggInvarMassMET30->Integral();ggInvarMassMET30->Scale(HiggsScale);ggInvarMassMET30->Write("BackgroundInvarMass");
  ggMETInvarMassSideBand_new->Rebin(2);
  HiggsScale=1./ggMETInvarMassSideBand_new->Integral(0,999);ggMETInvarMassSideBand_new->Scale(HiggsScale);ggMETInvarMassSideBand_new->Write("BackgroundMet");
  //newHiggsMet->Rebin(7);
  for(int i=1;i<newHiggsMet->GetNbinsX()+1;i++){
    if(newHiggsMet->GetBinContent(i)<=0 || newHiggsMet->GetBinLowEdge(i)<80){newHiggsMet->SetBinContent(i,1e-8);newHiggsMet->SetBinError(i,1e-9);}
    float x = newHiggsMet->GetBinContent(i)/newHiggsMet->GetBinWidth(i);
    float y = newHiggsMet->GetBinError(i)/newHiggsMet->GetBinWidth(i);
    newHiggsMet->SetBinContent(i,x);newHiggsMet->SetBinError(i,y);
  }
  //newHiggsMet->Rebin(7);
  HiggsScale=1./newHiggsMet->Integral();newHiggsMet->Scale(HiggsScale);
  newHiggsMet->Write("NewHiggsMet");
  newHiggsMet->GetYaxis()->SetRangeUser(1e-8,1);
  c1->SetLogy(1);
  newHiggsMet->Draw("PE");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/newHiggsMet_after.png");
  f_MetAndMassShapes.Close();
  fin.cd();  
  c1->SetLogy(1);

  TH1F *eeInvarMassFullRange = (TH1F*)fin.Get("eeInvarMassFullRange");
  TH1F *eeInvarMassFullRange_Pt25to40 = (TH1F*)fin.Get("eeInvarMassFullRange_Pt25to40");
  TH1F *eeInvarMassFullRange_Pt40to45 = (TH1F*)fin.Get("eeInvarMassFullRange_Pt40to45");
  TH1F *eeInvarMassFullRange_Pt45to50 = (TH1F*)fin.Get("eeInvarMassFullRange_Pt45to50");
  TH1F *eeInvarMassFullRange_Pt50to60 = (TH1F*)fin.Get("eeInvarMassFullRange_Pt50to60");
  TH1F *eeInvarMassFullRange_Pt60to80 = (TH1F*)fin.Get("eeInvarMassFullRange_Pt60to80");
  TH1F *eeInvarMassFullRange_Pt80 = (TH1F*)fin.Get("eeInvarMassFullRange_Pt80");
  TH1F* egInvarMass = (TH1F*)fin.Get("egInvarMass");
  TH1F *egInvarMass_Pt25to40 = (TH1F*)fin.Get("egInvarMass_Pt25to40");
  TH1F *egInvarMass_Pt40to45 = (TH1F*)fin.Get("egInvarMass_Pt40to45");
  TH1F *egInvarMass_Pt45to50 = (TH1F*)fin.Get("egInvarMass_Pt45to50");
  TH1F *egInvarMass_Pt50to60 = (TH1F*)fin.Get("egInvarMass_Pt50to60");
  TH1F *egInvarMass_Pt60to80 = (TH1F*)fin.Get("egInvarMass_Pt60to80");
  TH1F *egInvarMass_Pt80 = (TH1F*)fin.Get("egInvarMass_Pt80");

  for(int i=1;i<7;i++){
    c3->cd(i)->SetLogy(0);
  }
  c1->cd();
  RooRealVar x("x","ee Invariant Mass",61,121,"GeV");
  RooRealVar x2("x2","ee Invariant Mass",81,101,"GeV");
  RooDataHist data("data","dataset",x,eeInvarMassFullRange);
  /*  //Gaussian for peak
      RooRealVar gmean("gmean","gmean",91.2);
      RooRealVar gsigma("gsigma","gsigma",3.4);
      RooGaussian gauss("gauss","gauss",x,gmean,gsigma);
      //RooVoigtian
      RooRealVar vmean("vmean","vmean",90,85,95);
      RooRealVar vwidth("vwidth","vwidth",3,0,5);
      RooRealVar vsigma("vsigma","vsigma",2,0,5);
      RooVoigtian Vgauss("Vgauss","Vgauss",x,vmean,vwidth,vsigma);
      RooRealVar Vgauss_yield("Vgauss_yield","signal events",10000,0,1000000);*/
  //linear for background
  RooRealVar bkgd_poly_c1("bkgd_poly_c1","slope of background",0,-10.,10.);
  RooRealVar bkgd_poly_c2("bkgd_poly_c2","coefficient of x^2 term",0.,-10.,10.);
  //RooRealVar bkgd_poly_c3("bkgd_poly_c3","coefficient of x^3 term",0.,-10.,10.);
  //RooRealVar bkgd_poly_c4("bkgd_poly_c4","coefficient of x^4 term",0.,-10.,10.);
  //RooPolynomial background("background", "linear function for background",x,RooArgList(bkgd_poly_c1/*,bkgd_poly_c2/*,bkgd_poly_c3,bkgd_poly_c4*/) );
  //RooRealVar backgroundYield("bkgd yield","bkgd yield",6000,0,50000);
  //RooCMSShape for background
  RooRealVar aBkg("aBkg", "aBkg",/* 93*/88.0, 71.0, 111.0, "GeV"); //RooCMSShape (fixed)
  RooRealVar bBkg("bBkg", "bBkg",.14/*.2, 0., 0.6*/);            //RooCMSShape (fixed)
  RooRealVar c("c", "c", .14, 0., 2.0);                   //RooCMSShape (fixed)
  RooRealVar ZMass("bkgd peak", "ZMass",100./*71., 55.0, 111.0, "GeV"*/);            //Z mass 91.2=pdg
  RooCMSShape background("background", "background", x, aBkg, bBkg, c, ZMass);
  RooRealVar backgroundYield("bkgd events","bkgd events",200000,50000,700000);
  //aBkg.setConstant();
  //bBkg.setConstant();
  //c.setConstant();
  //Try other background models
  RooRealVar  bgtau("a_{BG}", "Background Shape", -0.15, -10.0, 10.0, "1/GeV/c^{2}");
  RooRealVar  bgaf("a_{BG}", "Background Shape", 0.0);
  //RooExponential background("bg", "Background Distribution", x, bgtau);
  //Breit-Wigner for peak
  RooRealVar  mRes("M_{Z^{0}}", "Z0 Resonance  Mass", 90.9, 88.0, 92.0);//,"GeV/c^{2}"); 
  RooRealVar  Gamma("#Gamma", "#Gamma", 2.0, 0.0,3.2);//,"GeV/c^{2}"); 
  //RooRealVar  Gamma("#Gamma", "#Gamma", 2.4952, 2.0,3.0);//,"GeV/c^{2}"); //default
  //mRes.setConstant();
  //Gamma.setConstant();
  RooBreitWigner bw("bw","A Breit-Wigner Distribution",x,mRes,Gamma);
  //Crystal Ball for resolution
  RooRealVar cbmean("cb_mean", "cbmean" , 0.0,-.25,.2) ;
  RooRealVar cbsigma("cb_sigma", "cbsigma" ,1.6,.8,2.4) ;
  //cbsigma.setConstant();
  //RooRealVar cbsig("signal events", "cbsignal", 75000., 50000, 100000);
  RooRealVar n("n","n", 5,0.,25.);
  RooRealVar alpha("alpha","alpha",0.6,0.2,1.0);
  RooCBShape cball("cball", "crystal ball", x, cbmean, cbsigma, alpha, n);
  //Convolution of Breit-Wigner and CrystalBall
  RooFFTConvPdf bw_cball("bw_cball","Convolution",x,bw,cball);
  RooRealVar bw_cball_yield("signal events","signal events",1500000,500000,2500000);
  //Addition of crystall ball  and Linear
  //RooAddPdf totalPDF("totalPDF","totalPDF",RooArgList(bkgd_linear, cball),RooArgList(bkgd_yield, cbsig));
  //RooAddPdf totalPDF("totalPDF","totalPDF",RooArgList(background, cball),RooArgList(backgroundYield, cbsig));
  RooAddPdf totalPDF("totalPDF","totalPDF",RooArgList(background,bw_cball),RooArgList(backgroundYield,bw_cball_yield));
  //fit to data
  totalPDF.fitTo(data,Extended(kTRUE));
  //Plot the whole thing
  RooPlot *xframe = x.frame(Title("Breit Wigner x Crystal Ball Signal Plus a 2nd Order Polynomial Background -- ee"));
  xframe->SetTitle("");
  totalPDF.paramOn(xframe, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  //totalPDF.paramOn(xframe, Format("NELU",AutoPrecision(1)),Layout(0.6,0.9,0.9) );
  //totalPDF.fitTo(data);
  //Vgauss.fitTo(data);
  data.plotOn(xframe,LineColor(kBlack));
  totalPDF.plotOn(xframe,LineColor(kBlue));
  //Vgauss.plotOn(xframe,LineColor(kBlue));
  //bwgauss.plotOn(xframe,LineColor(kBlue));
  //bw.plotOn(xframe,LineColor(kGreen));
  //totalPDF.plotOn(xframe,Components(bw_cball),LineColor(kCyan));
  totalPDF.plotOn(xframe, Components(background),/*LineStyle(Dashed),*/LineColor(kRed),DrawOption("F"),FillColor(kRed));
  totalPDF.plotOn(xframe, Components(bw_cball),LineColor(kGreen)/*,FillColor(kRed),DrawOption("F")*/);
  //totalPDF.plotOn(xframe, Components(bkgd_linear),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  //totalPDF.plotOn(xframe, Components(Vgauss),LineColor(kGreen));
  //totalPDF.plotOn(xframe, Components(cball),LineColor(kGreen));
  //bkgd_linear.plotOn(xframe,LineColor(kYellow));
  data.plotOn(xframe,LineColor(kBlack));
  c1->SetLogy(0);
  xframe->Draw();
  InvarMassTextee->Draw();
  CBsigEE=bw_cball_yield.getVal();
  CBsigEEerror=bw_cball_yield.getError();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_ee.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_ee.pdf");
  c1->SetLogy(1);
  xframe->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eeLOG.png");
  fout.cd();xframe->Write("InvarMassCrystalBall_ee");fin.cd();
    
  c1->SetLogy(0);
  // now pt dependent
  //25-40
  RooRealVar x25to40("x","ee Invariant Mass",61,121,"GeV");
  RooDataHist data25to40("data","dataset",x25to40,eeInvarMassFullRange_Pt25to40);
  //RooCMSShape for background
  RooRealVar aBkg25to40("aBkg", "aBkg", 98/*89.0, 71.0, 111.0, "GeV"*/); //RooCMSShape (fixed)
  RooRealVar bBkg25to40("bBkg", "bBkg", 0.14/*0.12, 0.08, 0.5*/);            //RooCMSShape (fixed)
  RooRealVar c25to40("c", "c", 0.21, .10, 0.5);                   //RooCMSShape (fixed)
  RooRealVar ZMass25to40("bkgd peak", "ZMass", 94/*79.0, 71.0, 111.0, "GeV"*/);            //Z mass 91.2=pdg
  RooRealVar backgroundYield25to40("bkgd events", "backgroundYield25to40", 60000,8000, 120000);
  RooCMSShape background25to40("background", "background", x25to40, aBkg25to40, bBkg25to40, c25to40, ZMass25to40);
  //linear bground
  RooRealVar bkgd_poly_c1_25to40("bkgd_poly_c1","coefficient of x^1 term",0.,-10,10);
  RooRealVar bkgd_poly_c2_25to40("bkgd_poly_c2","coefficient of x^2 term",0,-10,10);
  //RooPolynomial background25to40("bkgd_linear", "linear function for background",x25to40,RooArgList(bkgd_poly_c1_25to40,bkgd_poly_c2_25to40/*,bkgd_poly_c3,bkgd_poly_c4*/) );
  //RooRealVar backgroundYield25to40("bkg events","bkg events",6000,0,50000);
  //try other background shapes
  RooRealVar  bgtau25to40("a_{BG}", "Background Shape", -0.15, -10.0, 10.0, "1/GeV/c^{2}");
  RooRealVar  bgaf25to40("a_{BG}", "Background Shape", 0.0);
  //RooExponential background25to40("bg", "Background Distribution", x, bgtau25to40);
  //Breit-Wigner for peak
  RooRealVar  mRes25to40("M_{Z^{0}}", "Z0 Resonance  Mass", 90.9, 88.0, 92.0);//,"GeV/c^{2}"); 
  RooRealVar  Gamma25to40("#Gamma", "#Gamma", 2.0, 1.0,3.0);//,"GeV/c^{2}"); 
  RooBreitWigner bw25to40("bw","A Breit-Wigner Distribution",x25to40,mRes25to40,Gamma25to40);
  //Crystal Ball for resolution
  RooRealVar cbmean25to40("cb_mean", "cbmean" , 0.0,-1.0,1.0) ;
  RooRealVar cbsigma25to40("cb_sigma", "cbsigma" ,1.,0.1,5.0) ;
  //RooRealVar cbsig25to40("signal events", "cbsignal", 75000., 50000, 100000);
  RooRealVar n25to40("n","n", 50,0.,150.);
  RooRealVar alpha25to40("alpha","alpha",0.4,0.2,1.0);
  RooCBShape cball25to40("cball", "crystal ball", x25to40, cbmean25to40, cbsigma25to40, alpha25to40, n25to40);
  //Convolution of Breit-Wigner and CrystalBall
  RooFFTConvPdf bw_cball25to40("bw_cball","Convolution",x25to40,bw25to40,cball25to40);
  RooRealVar bw_cball_yield_25to40("signal events","signal events",600000,200000,1100000);
  //Crystal Ball for signal
  /*RooRealVar cbmean25to40("cb_mean", "cbmean" , 90., 88., 92.) ;
    RooRealVar cbsigma25to40("cb_sigma", "cbsigma" , 2.3, 1.0, 3.6) ;
    RooRealVar cbsig25to40("signal events", "cbsignal", 25000., 2000, 300000);
    RooRealVar n25to40("n","n", 5,1.,25);
    RooRealVar alpha25to40("alpha","alpha",0.9,0.5,1.3);
    RooCBShape cball25to40("cball", "crystal ball", x25to40, cbmean25to40, cbsigma25to40, alpha25to40, n25to40);*/
  //RooAddPdf totalPDF25to40("totalPDF25to40","totalPDF25to40",RooArgList(background25to40, cball25to40),RooArgList(backgroundYield25to40, cbsig25to40));
  //RooAddPdf totalPDF25to40("totalPDF25to40","totalPDF25to40",RooArgList(bkgd_linear_25to40, cball25to40),RooArgList(bkgd_yield_25to40, cbsig25to40));   
  RooAddPdf totalPDF25to40("totalPDF25to40","totalPDF25to40",RooArgList(background25to40, bw_cball25to40),RooArgList(backgroundYield25to40, bw_cball_yield_25to40));
  //fit to data
  totalPDF25to40.fitTo(data25to40,Extended());
  //Plot the whole thing
  RooPlot *xframe25to40 = x25to40.frame(Title("Breit Wigner x Crystal Ball Signal Plus a 2nd Order Polynomial Background -- ee Pt 25to40"));
  xframe25to40->SetTitle("");
  totalPDF25to40.paramOn(xframe25to40, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  data25to40.plotOn(xframe25to40,LineColor(kBlack));
  //totalPDF25to40.plotOn(xframe25to40, Components(background25to40),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDF25to40.plotOn(xframe25to40, Components(background25to40),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDF25to40.plotOn(xframe25to40, Components(bw_cball25to40),LineColor(kGreen)/*,FillColor(kRed),DrawOption("F")*/);
  totalPDF25to40.plotOn(xframe25to40,LineColor(kBlue));
  xframe25to40->Draw();
  CBsigEE25to40=bw_cball_yield_25to40.getVal();
  CBsigEEerror25to40=bw_cball_yield_25to40.getError();
  InvarMassTextee->Draw();
  Text25to40->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt25to40.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt25to40.pdf");
  c1->SetLogy(1);xframe25to40->Draw();c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt25to40LOG.png");c1->SetLogy(0);
  fout.cd();xframe25to40->Write("InvarMassCrystalBall_eePt25to40");fin.cd();
  c3->cd(1);xframe25to40->Draw();
  c1->cd();

  //40-45
  RooRealVar x40to45("x","ee Invariant Mass",61,121,"GeV");
  RooDataHist data40to45("data","dataset",x40to45,eeInvarMassFullRange_Pt40to45);
  //RooCMSShape for background
  RooRealVar aBkg40to45("aBkg", "aBkg", 70/*82.0, 71.0, 91.0, "GeV"*/); //RooCMSShape (fixed)
  RooRealVar bBkg40to45("bBkg", "bBkg40to45",.17/*0.17, 0.05, 0.3*/);            //RooCMSShape (fixed)
  RooRealVar c40to45("c", "c40to45", 0.25, 0.0, 0.5);                   //RooCMSShape (fixed)
  RooRealVar ZMass40to45("bkgd peak", "ZMass40to45", 91/*79.0, 71.0, 111.0, "GeV"*/);            //Z mass 91.2=pdg
  RooRealVar backgroundYield40to45("bkgd events", "backgroundYield", 4000,800, 10000);
  RooCMSShape background40to45("background", "background", x40to45, aBkg40to45, bBkg40to45, c40to45, ZMass40to45);
  //linear bground
  RooRealVar bkgd_poly_c1_40to45("bkgd_poly_c1","coefficient of x^1 term",0.,-10,10);
  RooRealVar bkgd_poly_c2_40to45("bkgd_poly_c2","coefficient of x^2 term",0,-100,100);
  //RooPolynomial background40to45("bkgd_linear", "linear function for background",x40to45,RooArgList(bkgd_poly_c1_40to45,bkgd_poly_c2_40to45/*,bkgd_poly_c3,bkgd_poly_c4*/) );
  //RooRealVar backgroundYield40to45("bkg events","bkg events",6000,0,50000);
  //Breit-Wigner for peak
  RooRealVar  mRes40to45("M_{Z^{0}}", "Z0 Resonance  Mass", 90., 87.0, 92.0);//,"GeV/c^{2}"); 
  RooRealVar  Gamma40to45("#Gamma", "#Gamma", 2.0, 0.5,3.0);//,"GeV/c^{2}"); 
  RooBreitWigner bw40to45("bw","A Breit-Wigner Distribution",x40to45,mRes40to45,Gamma40to45);
  //Crystal Ball for resolution
  RooRealVar cbmean40to45("cb_mean", "cbmean" , 0.0,-1.0,1.0) ;
  RooRealVar cbsigma40to45("cb_sigma", "cbsigma" ,1.,0.1,5.0) ;
  RooRealVar cbsig40to45("signal events", "cbsignal", 25000., 50000, 200000);
  RooRealVar n40to45("n","n", 50,0.,150.);
  RooRealVar alpha40to45("alpha","alpha",0.6,0.2,1.0);
  RooCBShape cball40to45("cball", "crystal ball", x40to45, cbmean40to45, cbsigma40to45, alpha40to45, n40to45);
  //Convolution of Breit-Wigner and CrystalBall
  RooFFTConvPdf bw_cball40to45("bw_cball","Convolution",x40to45,bw40to45,cball40to45);
  RooRealVar bw_cball_yield_40to45("signal events","signal events",400000,100000,800000);
  //Crystal Ball for signal
  /*RooRealVar cbmean40to45("cb_mean", "cbmean40to45" , 91., 88., 93.) ;
    RooRealVar cbsigma40to45("cb_sigma", "cbsigma40to45" , 2.0, 0.0, 5.0) ;
    RooRealVar cbsig40to45("signal events", "cbsignal", 50000., 2000, 80000);
    RooRealVar n40to45("n","n", 10,0.,180.);
    RooRealVar alpha40to45("alpha","alpha",1.4,0.8,2.);
    RooCBShape cball40to45("cball", "crystal ball", x40to45, cbmean40to45, cbsigma40to45, alpha40to45, n40to45);*/
  //RooAddPdf totalPDF40to45("totalPDF40to45","totalPDF40to45",RooArgList(background40to45, cball40to45),RooArgList(backgroundYield40to45, cbsig40to45));
  RooAddPdf totalPDF40to45("totalPDF40to45","totalPDF40to45",RooArgList(background40to45, bw_cball40to45),RooArgList(backgroundYield40to45, bw_cball_yield_40to45));
  //fit to data
  totalPDF40to45.fitTo(data40to45,Extended());
  //Plot the whole thing
  RooPlot *xframe40to45 = x40to45.frame(Title("Breit Wigner x Crystal Ball Signal Plus a 2nd Order Polynomial Background -- ee Pt 40to45"));
  xframe40to45->SetTitle("");
  totalPDF40to45.paramOn(xframe40to45, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  data40to45.plotOn(xframe40to45,LineColor(kBlack));
  //totalPDF40to45.plotOn(xframe40to45, Components(background40to45),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDF40to45.plotOn(xframe40to45, Components(background40to45),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDF40to45.plotOn(xframe40to45, Components(bw_cball40to45),LineColor(kGreen)/*,FillColor(kRed),DrawOption("F")*/);
  totalPDF40to45.plotOn(xframe40to45,LineColor(kBlue));
  xframe40to45->Draw();
  CBsigEE40to45=bw_cball_yield_40to45.getVal();
  CBsigEEerror40to45=bw_cball_yield_40to45.getError();InvarMassTextee->Draw();
  Text40to45->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt40to45.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt40to45.pdf");
  c1->SetLogy(1);xframe40to45->Draw();c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt40to45LOG.png");c1->SetLogy(0);
  fout.cd();xframe40to45->Write("InvarMassCrystalBall_eePt40to45");fin.cd();
  c3->cd(2);xframe40to45->Draw();
  c1->cd();
    
  //45-50
  RooRealVar x45to50("x","ee Invariant Mass",61,121,"GeV");
  RooDataHist data45to50("data","dataset",x45to50,eeInvarMassFullRange_Pt45to50);
  //RooCMSShape for background
  RooRealVar aBkg45to50("aBkg", "aBkg", 73/*89.0, 71.0, 111.0, "GeV"*/); //RooCMSShape (fixed)
  RooRealVar bBkg45to50("bBkg", "bBkg",.12/*0.12, 0.07, 0.33*/);            //RooCMSShape (fixed)
  RooRealVar c45to50("c", "c",0.18, 0.0, 0.4);                   //RooCMSShape (fixed)
  RooRealVar ZMass45to50("bkgd peak", "ZMass", 91/*80.0, 61.0, 111.0, "GeV"*/);            //Z mass 91.2=pdg
  RooRealVar backgroundYield45to50("bkgd events", "backgroundYield", 2000,0, 75000);
  RooCMSShape background45to50("background", "background", x45to50, aBkg45to50, bBkg45to50, c45to50, ZMass45to50);
  //linear bground
  RooRealVar bkgd_poly_c1_45to50("bkgd_poly_c1","coefficient of x^1 term",0.,-1,1);
  RooRealVar bkgd_poly_c2_45to50("bkgd_poly_c2","coefficient of x^2 term",0,-10,10);
  //RooPolynomial background45to50("bkgd_linear", "linear function for background",x45to50,RooArgList(bkgd_poly_c1_45to50,bkgd_poly_c2_45to50/*,bkgd_poly_c3,bkgd_poly_c4*/) );
  //    RooRealVar backgroundYield45to50("bkg events","bkg events",1500,0,5000);
  //Breit-Wigner for peak
  RooRealVar  mRes45to50("M_{Z^{0}}", "Z0 Resonance  Mass", 91, 88.0, 93.0);//,"GeV/c^{2}"); 
  RooRealVar  Gamma45to50("#Gamma", "#Gamma", 1.0, 0.0,2.0);//,"GeV/c^{2}"); 
  RooBreitWigner bw45to50("bw","A Breit-Wigner Distribution",x45to50,mRes45to50,Gamma45to50);
  //Crystal Ball for resolution
  RooRealVar cbmean45to50("cb_mean", "cbmean" , 0.0,-1.0,1.0) ;
  RooRealVar cbsigma45to50("cb_sigma", "cbsigma" ,1.,0.1,5.0) ;
  RooRealVar cbsig45to50("signal events", "cbsignal", 75000., 50000, 200000);
  RooRealVar n45to50("n","n", 50,0.,150.);
  RooRealVar alpha45to50("alpha","alpha",0.6,0.2,1.2);
  RooCBShape cball45to50("cball", "crystal ball", x45to50, cbmean45to50, cbsigma45to50, alpha45to50, n45to50);
  //Convolution of Breit-Wigner and CrystalBall
  RooFFTConvPdf bw_cball45to50("bw_cball","Convolution",x45to50,bw45to50,cball45to50);
  RooRealVar bw_cball_yield_45to50("signal events","signal events",500000,100000,1000000);
  //Crystal Ball for signal
  /* RooRealVar cbmean45to50("cb_mean", "cbmean" , 91., 88., 93.) ;
     RooRealVar cbsigma45to50("cb_sigma", "cbsigma" , 2.0, 0.0, 5.0) ;
     RooRealVar cbsig45to50("signal events", "cbsignal", 25000., 2000, 80000);
     RooRealVar n45to50("n","n", 20,0.,80.);
     RooRealVar alpha45to50("alpha","alpha",1.4,0.8,2.);
     RooCBShape cball45to50("cball", "crystal ball", x45to50, cbmean45to50, cbsigma45to50, alpha45to50, n45to50);*/
  //RooAddPdf totalPDF45to50("totalPDF45to50","totalPDF45to50",RooArgList(background45to50, cball45to50),RooArgList(backgroundYield45to50, cbsig45to50));
  RooAddPdf totalPDF45to50("totalPDF45to50","totalPDF45to50",RooArgList(background45to50, bw_cball45to50),RooArgList(backgroundYield45to50, bw_cball_yield_45to50));
  //fit to data
  totalPDF45to50.fitTo(data45to50,Extended());
  //Plot the whole thing
  RooPlot *xframe45to50 = x45to50.frame(Title("Breit Wigner x Crystal Ball Signal Plus a 2nd Order Polynomial Background -- ee Pt 45to50"));
  xframe45to50->SetTitle("");
  totalPDF45to50.paramOn(xframe45to50, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  data45to50.plotOn(xframe45to50,LineColor(kBlack));
  //totalPDF45to50.plotOn(xframe45to50, Components(background45to50),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDF45to50.plotOn(xframe45to50, Components(background45to50),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDF45to50.plotOn(xframe45to50, Components(bw_cball45to50),LineColor(kGreen)/*,FillColor(kRed),DrawOption("F")*/);
  totalPDF45to50.plotOn(xframe45to50,LineColor(kBlue));
  xframe45to50->Draw();
  CBsigEE45to50=bw_cball_yield_45to50.getVal();
  CBsigEEerror45to50=bw_cball_yield_45to50.getError();InvarMassTextee->Draw();
  Text45to50->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt45to50.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt45to50.pdf");
  c1->SetLogy(1);xframe45to50->Draw();c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt45to50LOG.png");c1->SetLogy(0);
  fout.cd();xframe45to50->Write("InvarMassCrystalBall_eePt45to50");fin.cd();
  c3->cd(3);xframe45to50->Draw();
  c1->cd();
  
  //50to60
  RooRealVar x50to60("x","ee Invariant Mass",61,121,"GeV");
  RooDataHist data50to60("data","dataset",x50to60,eeInvarMassFullRange_Pt50to60);
  //RooCMSShape for background
  RooRealVar aBkg50to60("aBkg", "aBkg", 97/*89.0, 71.0, 111.0, "GeV"*/); //RooCMSShape (fixed)
  RooRealVar bBkg50to60("bBkg", "bBkg",.165/*0.2, 0.0, 0.5*/);            //RooCMSShape (fixed)
  RooRealVar c50to60("c", "c",0.1, 0.0, 0.22);                   //RooCMSShape (fixed)
  RooRealVar ZMass50to60("bkgd peak", "ZMass", 90.9/*75.0, 65.0, 85.0, "GeV"*/);            //Z mass 91.2=pdg
  RooRealVar backgroundYield50to60("bkgd events", "backgroundYield", 50000,10000, 100000);
  RooCMSShape background50to60("background", "background", x50to60, aBkg50to60, bBkg50to60, c50to60, ZMass50to60);
  //linear bground
  RooRealVar bkgd_poly_c1_50to60("bkgd_poly_c1","coefficient of x^1 term",0.,-1,1);
  RooRealVar bkgd_poly_c2_50to60("bkgd_poly_c2","coefficient of x^2 term",0,-2,2);
  //RooPolynomial background50to60("bkgd_linear", "linear function for background",x50to60,RooArgList(bkgd_poly_c1_50to60,bkgd_poly_c2_50to60/*,bkgd_poly_c3,bkgd_poly_c4*/) );
  //    RooRealVar backgroundYield50to60("bkg events","bkg events",6000,0,50000);
  //Breit-Wigner for peak
  RooRealVar  mRes50to60("M_{Z^{0}}", "Z0 Resonance  Mass", 90.9, 88.0, 92.0);//,"GeV/c^{2}"); 
  RooRealVar  Gamma50to60("#Gamma", "#Gamma", 2.0, 0.8,4);//,"GeV/c^{2}"); 
  RooBreitWigner bw50to60("bw","A Breit-Wigner Distribution",x50to60,mRes50to60,Gamma50to60);
  //Crystal Ball for resolution
  RooRealVar cbmean50to60("cb_mean", "cbmean" , 0.0,-1.0,1.0) ;
  RooRealVar cbsigma50to60("cb_sigma", "cbsigma" ,1.,0.1,5.0) ;
  //RooRealVar cbsig50to60("signal events", "cbsignal", 85000., 50000, 200000);
  RooRealVar n50to60("n","n", 50,0.,150.);
  RooRealVar alpha50to60("alpha","alpha",0.6,0.2,1.2);
  RooCBShape cball50to60("cball", "crystal ball", x50to60, cbmean50to60, cbsigma50to60, alpha50to60, n50to60);
  //Convolution of Breit-Wigner and CrystalBall
  RooFFTConvPdf bw_cball50to60("bw_cball","Convolution",x50to60,bw50to60,cball50to60);
  RooRealVar bw_cball_yield_50to60("signal events","signal events",350000,100000,800000);
  //Crystal Ball for signal
  /* RooRealVar cbmean50to60("cb_mean", "cbmean" , 91., 88., 93.) ;
     RooRealVar cbsigma50to60("cb_sigma", "cbsigma" , 2.0, 0.0, 5.0) ;
     RooRealVar cbsig50to60("signal events", "cbsignal", 5000., 2000, 80000);
     RooRealVar n50to60("n","n", 10,0.,200.);
     RooRealVar alpha50to60("alpha","alpha",1.4,0.8,5.);
     RooCBShape cball50to60("cball", "crystal ball", x50to60, cbmean50to60, cbsigma50to60, alpha50to60, n50to60);*/
  //RooAddPdf totalPDF50to60("totalPDF50to60","totalPDF50to60",RooArgList(background50to60, cball50to60),RooArgList(backgroundYield50to60, cbsig50to60));
  RooAddPdf totalPDF50to60("totalPDF50to60","totalPDF50to60",RooArgList(background50to60, bw_cball50to60),RooArgList(backgroundYield50to60, bw_cball_yield_50to60));
  //fit to data
  totalPDF50to60.fitTo(data50to60,Extended());
  //Plot the whole thing
  RooPlot *xframe50to60 = x50to60.frame(Title("Breit Wigner x Crystal Ball Signal Plus a 2nd Order Polynomial Background -- ee Pt 50to60"));
  xframe50to60->SetTitle("");
  totalPDF50to60.paramOn(xframe50to60, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  data50to60.plotOn(xframe50to60,LineColor(kBlack));
  //totalPDF50to60.plotOn(xframe50to60, Components(background50to60),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDF50to60.plotOn(xframe50to60, Components(background50to60),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDF50to60.plotOn(xframe50to60, Components(bw_cball50to60),LineColor(kGreen)/*,FillColor(kRed),DrawOption("F")*/);
  totalPDF50to60.plotOn(xframe50to60,LineColor(kBlue));
  xframe50to60->Draw();
  CBsigEE50to60=bw_cball_yield_50to60.getVal();
  CBsigEEerror50to60=bw_cball_yield_50to60.getError();InvarMassTextee->Draw();
  Text50to60->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt50to60.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt50to60.pdf");
  c1->SetLogy(1);xframe50to60->Draw();c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt50to60LOG.png");c1->SetLogy(0);
  fout.cd();xframe50to60->Write("InvarMassCrystalBall_eePt50to60");fin.cd();
  c3->cd(4);xframe50to60->Draw();
  c1->cd();

  //60to80
  RooRealVar x60to80("x","ee Invariant Mass",61,121,"GeV");
  RooDataHist data60to80("data","dataset",x60to80,eeInvarMassFullRange_Pt60to80);
  //RooCMSShape for background
  RooRealVar aBkg60to80("aBkg", "aBkg",95/*89.0, 71.0, 111.0, "GeV"*/); //RooCMSShape (fixed)
  RooRealVar bBkg60to80("bBkg", "bBkg",.14/*0.17, 0.1, 0.5*/);            //RooCMSShape (fixed)
  RooRealVar c60to80("c", "c",0.18, -0.1, 0.3);                   //RooCMSShape (fixed)
  RooRealVar ZMass60to80("bkgd peak", "ZMass", 91.2/*79.0, 71.0, 111.0, "GeV"*/);            //Z mass 91.2=pdg
  RooRealVar backgroundYield60to80("bkgd events", "backgroundYield", 3500,1000, 7000);
  RooCMSShape background60to80("background", "background", x60to80, aBkg60to80, bBkg60to80, c60to80, ZMass60to80);
  //linear bground
  RooRealVar bkgd_poly_c1_60to80("bkgd_poly_c1","coefficient of x^1 term",1.5,-0.5,3.5);
  RooRealVar bkgd_poly_c2_60to80("bkgd_poly_c2","coefficient of x^2 term",0,-10,10);
  //RooPolynomial background60to80("bkgd_linear", "linear function for background",x60to80,RooArgList(bkgd_poly_c1_60to80,bkgd_poly_c2_60to80/*,bkgd_poly_c3,bkgd_poly_c4*/) );
  //    RooRealVar backgroundYield60to80("bkg events","bkg events",1900,1000,5000);
  //Breit-Wigner for peak
  RooRealVar  mRes60to80("M_{Z^{0}}", "Z0 Resonance  Mass", 90.9, 88.0, 92.0);//,"GeV/c^{2}"); 
  RooRealVar  Gamma60to80("#Gamma", "#Gamma", 2.0, 1.0,3.5);//,"GeV/c^{2}"); 
  RooBreitWigner bw60to80("bw","A Breit-Wigner Distribution",x60to80,mRes60to80,Gamma60to80);
  //Crystal Ball for resolution
  RooRealVar cbmean60to80("cb_mean", "cbmean" , 0.0,-1.0,1.0) ;
  RooRealVar cbsigma60to80("cb_sigma", "cbsigma" ,2.,0.1,5.0) ;
  RooRealVar cbsig60to80("signal events", "cbsignal", 75000., 50000, 1000000);
  RooRealVar n60to80("n","n", 50,0.,160.);
  RooRealVar alpha60to80("alpha","alpha",0.7,0.2,1.3);
  RooCBShape cball60to80("cball", "crystal ball", x60to80, cbmean60to80, cbsigma60to80, alpha60to80, n60to80);
  //Convolution of Breit-Wigner and CrystalBall
  RooFFTConvPdf bw_cball60to80("bw_cball","Convolution",x60to80,bw60to80,cball60to80);
  RooRealVar bw_cball_yield_60to80("signal events","signal events",40000,10000,100000);
  //Crystal Ball for signal
  /*  RooRealVar cbmean60to80("cb_mean", "cbmean" , 91., 88., 93.) ;
      RooRealVar cbsigma60to80("cb_sigma", "cbsigma" , 2.0, 0.0, 5.0) ;
      RooRealVar cbsig60to80("signal events", "cbsignal", 8000., 2000, 20000);
      RooRealVar n60to80("n","n", 75,0.,250.);
      RooRealVar alpha60to80("alpha","alpha",1.4,0.8,2.);
      RooCBShape cball60to80("cball", "crystal ball", x60to80, cbmean60to80, cbsigma60to80, alpha60to80, n60to80);*/
  //RooAddPdf totalPDF60to80("totalPDF60to80","totalPDF60to80",RooArgList(background60to80, cball60to80),RooArgList(backgroundYield60to80, cbsig60to80));
  RooAddPdf totalPDF60to80("totalPDF60to80","totalPDF60to80",RooArgList(background60to80, bw_cball60to80),RooArgList(backgroundYield60to80, bw_cball_yield_60to80));
  //fit to data
  totalPDF60to80.fitTo(data60to80,Extended());
  //Plot the whole thing
  RooPlot *xframe60to80 = x60to80.frame(Title("Breit Wigner x Crystal Ball Signal Plus a 2nd Order Polynomial Background -- ee Pt 60to80"));
  xframe60to80->SetTitle("");
  totalPDF60to80.paramOn(xframe60to80, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  data60to80.plotOn(xframe60to80,LineColor(kBlack));
  //totalPDF60to80.plotOn(xframe60to80, Components(background60to80),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDF60to80.plotOn(xframe60to80, Components(background60to80),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDF60to80.plotOn(xframe60to80, Components(bw_cball60to80),LineColor(kGreen)/*,FillColor(kRed),DrawOption("F")*/);
  totalPDF60to80.plotOn(xframe60to80,LineColor(kBlue));
  data60to80.plotOn(xframe60to80,LineColor(kBlack));
  xframe60to80->Draw();
  CBsigEE60to80=bw_cball_yield_60to80.getVal();
  CBsigEEerror60to80=bw_cball_yield_60to80.getError();InvarMassTextee->Draw();
  Text60to80->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt60to80.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt60to80.pdf");
  c1->SetLogy(1);xframe60to80->Draw();c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt60to80LOG.png");c1->SetLogy(0);
  fout.cd();xframe60to80->Write("InvarMassCrystalBall_eePt60to80");fin.cd();
  c3->cd(5);xframe60to80->Draw();
  c1->cd();
    
  //over80    
  RooRealVar x80("x","ee Invariant Mass",61,121,"GeV");
  RooDataHist data80("data","dataset",x,eeInvarMassFullRange_Pt80);
  //RooCMSShape for background
  RooRealVar aBkg80("aBkg", "aBkg",97/*89.0, 71.0, 111.0, "GeV"*/); //RooCMSShape (fixed)
  RooRealVar bBkg80("bBkg", "bBkg",.17/*0.17, 0.1, 0.8*/);            //RooCMSShape (fixed)
  RooRealVar c80("c", "c",0.18,-0.1,0.3);                   //RooCMSShape (fixed)
  RooRealVar ZMass80("bkgd peak", "ZMass",90.5/*79.0, 71.0, 111.0, "GeV"*/);            //Z mass 91.2=pdg
  RooRealVar backgroundYield80("bkgd events", "backgroundYield", 1500,250, 40000);
  RooCMSShape background80("background", "background", x80, aBkg80, bBkg80, c80, ZMass80);
  //linear bground
  RooRealVar bkgd_poly_c1_80("bkgd_poly_c1","coefficient of x^1 term",0,-8,8);
  RooRealVar bkgd_poly_c2_80("bkgd_poly_c2","coefficient of x^2 term",0,-10,10);
  //RooPolynomial background80("bkgd_linear", "linear function for background",x80,RooArgList(bkgd_poly_c1_80/*,bkgd_poly_c2_80/*,bkgd_poly_c3,bkgd_poly_c4*/) );
  //    RooRealVar backgroundYield80("bkg events","bkg events",800,500,2000);
  //Breit-Wigner for peak
  RooRealVar  mRes80("M_{Z^{0}}", "Z0 Resonance  Mass", 90.9, 88.0, 92.0);//,"GeV/c^{2}"); 
  RooRealVar  Gamma80("#Gamma", "#Gamma", 2.0, 1.0,3.5);//,"GeV/c^{2}"); 
  RooBreitWigner bw80("bw","A Breit-Wigner Distribution",x80,mRes80,Gamma80);
  //Crystal Ball for resolution
  RooRealVar cbmean80("cb_mean", "cbmean" , 0.0,-1.0,1.0) ;
  RooRealVar cbsigma80("cb_sigma", "cbsigma" ,1.,0.1,5.0) ;
  //RooRealVar cbsig80("signal events", "cbsignal", 75000., 50000, 100000);
  RooRealVar n80("n","n", 75,0.,175.);
  RooRealVar alpha80("alpha","alpha",1.,0.2,1.8);
  RooCBShape cball80("cball", "crystal ball", x80, cbmean80, cbsigma80, alpha80, n80);
  //Convolution of Breit-Wigner and CrystalBall
  RooFFTConvPdf bw_cball80("bw_cball","Convolution",x80,bw80,cball80);
  RooRealVar bw_cball_yield_80("signal events","signal events",15000,8000,50000);
  //Crystal Ball for signal
  /* RooRealVar cbmean80("cb_mean", "cbmean" , 91., 88., 93.) ;
     RooRealVar cbsigma80("cb_sigma", "cbsigma" , 2.0, 0.0, 5.0) ;
     RooRealVar cbsig80("signal events", "cbsignal", 3500., 2000, 8000);
     RooRealVar n80("n","n", 150,0.,250.);
     RooRealVar alpha80("alpha","alpha",1.4,0.8,2.);
     RooCBShape cball80("cball", "crystal ball", x80, cbmean80, cbsigma80, alpha80, n80);*/
  //RooAddPdf totalPDF80("totalPDF80","totalPDF80",RooArgList(background80, cball80),RooArgList(backgroundYield80, cbsig80));
  RooAddPdf totalPDF80("totalPDF80","totalPDF80",RooArgList(background80, bw_cball80),RooArgList(backgroundYield80, bw_cball_yield_80));
  //fit to data
  totalPDF80.fitTo(data80,Extended());
  //Plot the whole thing
  RooPlot *xframe80 = x80.frame(Title("Breit Wigner x Crystal Ball Signal Plus a 2nd Order Polynomial Background -- ee Pt > 80"));
  xframe80->SetTitle("");
  totalPDF80.paramOn(xframe80, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  data80.plotOn(xframe80,LineColor(kBlack));
  //totalPDF80.plotOn(xframe80, Components(background80),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDF80.plotOn(xframe80, Components(background80),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDF80.plotOn(xframe80, Components(bw_cball80),LineColor(kGreen)/*,FillColor(kRed),DrawOption("F")*/);
  totalPDF80.plotOn(xframe80,LineColor(kBlue));
  data80.plotOn(xframe80,LineColor(kBlack));
  xframe80->Draw();
  CBsigEE80=bw_cball_yield_80.getVal();
  CBsigEEerror80=bw_cball_yield_80.getError();InvarMassTextee->Draw();
  Text80->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt80.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt80.pdf");
  c1->SetLogy(1);xframe80->Draw();c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt80LOG.png");c1->SetLogy(0);
  fout.cd();xframe80->Write("InvarMassCrystalBall_eePt80");fin.cd();
  c3->cd(6);xframe80->Draw();
  c3->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt.png");
  c3->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eePt.pdf");
  c1->cd();

  //now eg

  RooRealVar xeg("x","eg Invariant Mass",61,121,"GeV");
  RooDataHist dataeg("dataeg","dataset",xeg,egInvarMass);
  //linear for background
  /*RooRealVar bkgd_poly_c1("slope of bkgd","slope of background",7,0,1500);
    RooPolynomial bkgd_linear("bkgd_linear", "linear function for background",xeg,RooArgList(bkgd_poly_c1) );
    RooRealVar bkgd_yield("bkg events","bkg events",100,0,1000000);*/

  //RooCMSShape for background
  RooRealVar aBkgeg("aBkg", "aBkgeg", /*75*/70., 61.0, 80.0, "GeV"); //RooCMSShape (fixed)
  RooRealVar bBkgeg("bBkg", "bBkgeg", 0.1, 0.005, 0.4);            //RooCMSShape (fixed)
  RooRealVar ceg("c", "ceg", 0.05, 0.01, 0.1);                   //RooCMSShape (fixed)
  RooRealVar ZMasseg("bkgd peak", "ZMasseg", 91.2/*75.0, 61.0, 111.0, "GeV"*/);            //Z mass 91.2=pdg
  //RooRealVar aBkg("aBkg", "aBkg", 52.9, 40.0, 100.0, "GeV"); //RooCMSShape (fixed)
  //RooRealVar bBkg("bBkg", "bBkg", 0.47, 0.1, 0.6);            //RooCMSShape (fixed)
  //RooRealVar c("c", "c", 1.4, 0.0, 5.0);                   //RooCMSShape (fixed)
  RooRealVar backgroundYieldeg("bkgd events", "backgroundYieldeg", 30000,15000, 80000);
  RooCMSShape backgroundeg("backgroundeg", "backgroundeg", xeg, aBkgeg, bBkgeg, ceg, ZMasseg);

  //Breit-Wigner for peak
  RooRealVar  mReseg("M_{Z^{0}}", "Z0 Resonance  Mass", 90.9, 88.0, 92.0);//,"GeV/c^{2}"); 
  RooRealVar  Gammaeg("#Gamma", "#Gamma", 1.2, 0.1,3.5);//,"GeV/c^{2}"); 
  RooBreitWigner bweg("bw","A Breit-Wigner Distribution",xeg,mReseg,Gammaeg);
  //Crystal Ball for resolution
  RooRealVar cbmeaneg("cb_mean", "cbmean" , 0.0,-1.0,1.0) ;
  RooRealVar cbsigmaeg("cb_sigma", "cbsigma" ,2.,0.1,4.0) ;
  //RooRealVar cbsigeg("signal events", "cbsignal", 75000., 50000, 100000);
  RooRealVar neg("n","n", 80,5.,150.);
  RooRealVar alphaeg("alpha","alpha",1.,0.2,1.8);
  RooCBShape cballeg("cball", "crystal ball", xeg, cbmeaneg, cbsigmaeg, alphaeg, neg);
  //Convolution of Breit-Wigner and CrystalBall
  RooFFTConvPdf bw_cballeg("bw_cball","Convolution",xeg,bweg,cballeg);
  RooRealVar bw_cball_yield_eg("signal events","signal events",150000,65000,400000);
  //Crystal Ball
  /*  RooRealVar cbmeaneg("mass", "cbmeaneg" , 91., 87., 95.) ;
      RooRealVar cbsigmaeg("cb_sigma", "cbsigmaeg" , 2.5, 1., 4.) ;
      RooRealVar cbsigeg("signal events", "cbsignaleg", 11000, 6000, 40000);
      RooRealVar neg("n","neg", 35.,2.,100.);
      RooRealVar alphaeg("alpha","alphaeg", 1.1,0.,2.);
      RooCBShape cballeg("cball", "crystal ball eg", xeg, cbmeaneg, cbsigmaeg, alphaeg, neg);*/
  //Addition of crystall ball  and Linear
  //RooAddPdf totalPDFeg("totalPDFeg","totalPDFeg",RooArgList(backgroundeg, cballeg),RooArgList(backgroundYieldeg, cbsigeg));    
  RooAddPdf totalPDFeg("totalPDFeg","totalPDFeg",RooArgList(backgroundeg, bw_cballeg),RooArgList(backgroundYieldeg, bw_cball_yield_eg));

  //fit to data
  totalPDFeg.fitTo(dataeg,Extended());
  //Plot the whole thing
  RooPlot *xframeeg = xeg.frame(Title("Crystal Ball Signal Plus a RooCMSShape Background -- eg"));
  xframeeg->SetTitle("");
  totalPDFeg.paramOn(xframeeg, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  dataeg.plotOn(xframeeg,LineColor(kBlack));
  totalPDFeg.plotOn(xframeeg, Components(backgroundeg),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDFeg.plotOn(xframeeg, Components(bw_cballeg),LineColor(kGreen)/*,FillColor(kRed),DrawOption("F")*/);
  totalPDFeg.plotOn(xframeeg,LineColor(kBlue));
  dataeg.plotOn(xframeeg,LineColor(kBlack));
  //bkgd_linear.plotOn(xframeeg,LineColor(kYellow));
  xframeeg->Draw();
  InvarMassTexteg->Draw();
  //CBsigEG=cbsigeg.getVal();
  //CBsigEGerror=cbsigeg.getError();
  CBsigEG=bw_cball_yield_eg.getVal();
  CBsigEGerror=bw_cball_yield_eg.getError();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eg.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eg.pdf");
  fout.cd();xframeeg->Write("InvarMassCrystalBall_eg");fin.cd();
  c4->cd(1);
  xframe->Draw();
  InvarMassTextee->Draw();
  c4->cd(2);
  xframeeg->Draw();
  InvarMassTexteg->Draw();
  c4->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eeANDeg.png");
  c4->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eeANDeg.pdf");
  c1->cd();
  //Pt25to40
  RooRealVar xeg25to40("x","eg Invariant Mass",61,121,"GeV");
  RooDataHist dataeg25to40("dataeg25to40","dataset",xeg25to40,egInvarMass_Pt25to40);
  //RooCMSShape for background
  RooRealVar aBkgeg25to40("aBkg", "aBkgeg", 76./*, 61.0, 111.0, "GeV"*/); //RooCMSShape (fixed)
  RooRealVar bBkgeg25to40("bBkg", "bBkgeg", .08/*0.04, 0.01, 0.15*/);            //RooCMSShape (fixed)
  RooRealVar ceg25to40("c", "ceg",0.045, 0.02, 0.1);                   //RooCMSShape (fixed)
  RooRealVar ZMasseg25to40("bkgd peak", "ZMasseg", 91.0/*, 87.0, 93.0, "GeV"*/);            //Z mass 91.2=pdg
  RooRealVar backgroundYieldeg25to40("bkgd events", "backgroundYieldeg", 7000,600, 120000);
  RooCMSShape backgroundeg25to40("backgroundeg", "backgroundeg", xeg25to40, aBkgeg25to40, bBkgeg25to40, ceg25to40, ZMasseg25to40);
  //Breit-Wigner for peak
  RooRealVar  mReseg25to40("M_{Z^{0}}", "Z0 Resonance  Mass", 90.9, 88.0, 92.5);//,"GeV/c^{2}"); 
  RooRealVar  Gammaeg25to40("#Gamma", "#Gamma",1,-1,5);//,"GeV/c^{2}"); 
  RooBreitWigner bweg25to40("bw","A Breit-Wigner Distribution",xeg25to40,mReseg25to40,Gammaeg25to40);
  //Crystal Ball for resolution
  RooRealVar cbmeaneg25to40("cb_mean", "cbmean" , 0.0,-5.0,2.0) ;
  RooRealVar cbsigmaeg25to40("cb_sigma", "cbsigma" ,1.,0.1,5.0) ;
  RooRealVar cbsigeg25to40("signal events", "cbsignal", 75000., 50000, 100000);
  RooRealVar neg25to40("n","n", 50,0.,150.);
  RooRealVar alphaeg25to40("alpha","alpha",0.6,0.2,3);
  RooCBShape cballeg25to40("cball", "crystal ball", xeg25to40, cbmeaneg25to40, cbsigmaeg25to40, alphaeg25to40, neg25to40);
  //Convolution of Breit-Wigner and CrystalBall
  RooFFTConvPdf bw_cballeg25to40("bw_cball","Convolution",xeg25to40,bweg25to40,cballeg25to40);
  RooRealVar bw_cball_yield_eg25to40("signal events","signal events",4000,2000,100000);   
  /*  //Crystal Ball
      RooRealVar cbmeaneg25to40("mass", "cbmeaneg" , 91., 87., 95.) ;
      RooRealVar cbsigmaeg25to40("cb_sigma", "cbsigmaeg" , 2.5, 1., 4.) ;
      RooRealVar cbsigeg25to40("signal events", "cbsignaleg", 2000, 500, 10000);
      RooRealVar neg25to40("n","neg", 25.,20.,45.);
      RooRealVar alphaeg25to40("alpha","alphaeg", 1.1,0,2);
      RooCBShape cballeg25to40("cball", "crystal ball eg", xeg25to40, cbmeaneg25to40, cbsigmaeg25to40, alphaeg25to40, neg25to40);*/
  //Addition of crystall ball  and RooCMSShape
  //RooAddPdf totalPDFeg25to40("totalPDFeg25to40","totalPDFeg25to40",RooArgList(backgroundeg25to40, cballeg25to40),RooArgList(backgroundYieldeg25to40, cbsigeg25to40));   
  RooAddPdf totalPDFeg25to40("totalPDFeg25to40","totalPDFeg25to40",RooArgList(backgroundeg25to40, bw_cballeg25to40),RooArgList(backgroundYieldeg25to40, bw_cball_yield_eg25to40));
  //fit to data
  totalPDFeg25to40.fitTo(dataeg25to40,Extended());
  //Plot the whole thing
  RooPlot *xframeeg25to40 = xeg25to40.frame(Title("Crystal Ball Signal Plus a RooCMSShape Background -- eg25to40"));
  xframeeg25to40->SetTitle("");
  totalPDFeg25to40.paramOn(xframeeg25to40, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  dataeg25to40.plotOn(xframeeg25to40,LineColor(kBlack));
  totalPDFeg25to40.plotOn(xframeeg25to40, Components(backgroundeg25to40),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDFeg25to40.plotOn(xframeeg25to40, Components(bw_cballeg25to40),LineColor(kGreen)/*,FillColor(kRed),DrawOption("F")*/);
  totalPDFeg25to40.plotOn(xframeeg25to40,LineColor(kBlue));
  dataeg25to40.plotOn(xframeeg25to40,LineColor(kBlack));
  xframeeg25to40->Draw();
  //CBsigEG25to40=cbsigeg25to40.getVal();
  //CBsigEGerror25to40=cbsigeg25to40.getError();
  CBsigEG25to40=bw_cball_yield_eg25to40.getVal();
  CBsigEGerror25to40=bw_cball_yield_eg25to40.getError();
  InvarMassTexteg->Draw();
  Text25to40->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eg25to40.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eg25to40.pdf");
  fout.cd();xframeeg25to40->Write("InvarMassCrystalBall_eg25to40");fin.cd();
  c3->cd(1);xframeeg25to40->Draw();
  c1->cd();
   
  //Pt40to45
  RooRealVar xeg40to45("x","eg Invariant Mass",61,121,"GeV");
  RooDataHist dataeg40to45("dataeg40to45","dataset",xeg40to45,egInvarMass_Pt40to45);
  //RooCMSShape for background
  RooRealVar aBkgeg40to45("aBkg", "aBkgeg", 70/*80., 61.0, 111.0, "GeV"*/); //RooCMSShape (fixed)
  RooRealVar bBkgeg40to45("bBkg", "bBkgeg",.12/*0.2, 0.01, 0.4*/);            //RooCMSShape (fixed)
  RooRealVar ceg40to45("c", "ceg", 0.06, 0.02, 0.12);                   //RooCMSShape (fixed)
  RooRealVar ZMasseg40to45("bkgd peak", "ZMasseg", 91/*91.0, 87.0, 93.0, "GeV"*/);            //Z mass 91.2=pdg
  RooRealVar backgroundYieldeg40to45("bkgd events", "backgroundYieldeg", 3000,0, 60000);
  RooCMSShape backgroundeg40to45("backgroundeg", "backgroundeg", xeg40to45, aBkgeg40to45, bBkgeg40to45, ceg40to45, ZMasseg40to45);
  //Breit-Wigner for peak
  RooRealVar  mReseg40to45("M_{Z^{0}}", "Z0 Resonance  Mass", 90.9, 88.0, 92.0);//,"GeV/c^{2}"); 
  RooRealVar  Gammaeg40to45("#Gamma", "#Gamma", 1.0,-1,3);//,"GeV/c^{2}"); 
  RooBreitWigner bweg40to45("bw","A Breit-Wigner Distribution",xeg40to45,mReseg40to45,Gammaeg40to45);
  //Crystal Ball for resolution
  RooRealVar cbmeaneg40to45("cb_mean", "cbmean" , 0.0,-1.0,2.0) ;
  RooRealVar cbsigmaeg40to45("cb_sigma", "cbsigma" ,1.,0.1,5.0) ;
  RooRealVar cbsigeg40to45("signal events", "cbsignal", 75000., 50000, 100000);
  RooRealVar neg40to45("n","n", 50,0.,150.);
  RooRealVar alphaeg40to45("alpha","alpha",0.6,0.01,3);
  RooCBShape cballeg40to45("cball", "crystal ball", xeg40to45, cbmeaneg40to45, cbsigmaeg40to45, alphaeg40to45, neg40to45);
  //Convolution of Breit-Wigner and CrystalBall
  RooFFTConvPdf bw_cballeg40to45("bw_cball","Convolution",xeg40to45,bweg40to45,cballeg40to45);
  RooRealVar bw_cball_yield_eg40to45("signal events","signal events",3000,1000,100000); 
  /*   //Crystal Ball
       RooRealVar cbmeaneg40to45("mass", "cbmeaneg" , 91., 87., 95.) ;
       RooRealVar cbsigmaeg40to45("cb_sigma", "cbsigmaeg" , 2.5, 1., 4.) ;
       RooRealVar cbsigeg40to45("signal events", "cbsignaleg", 1400, 500, 4000);
       RooRealVar neg40to45("n","neg", 40.,0.,90.);
       RooRealVar alphaeg40to45("alpha","alphaeg", 1.1,0,2);
       RooCBShape cballeg40to45("cball", "crystal ball eg", xeg40to45, cbmeaneg40to45, cbsigmaeg40to45, alphaeg40to45, neg40to45);*/
  //Addition of crystall ball  and RooCMSShape
  //RooAddPdf totalPDFeg40to45("totalPDFeg40to45","totalPDFeg40to45",RooArgList(backgroundeg40to45, cballeg40to45),RooArgList(backgroundYieldeg40to45, cbsigeg40to45)); 
  RooAddPdf totalPDFeg40to45("totalPDFeg40to45","totalPDFeg40to45",RooArgList(backgroundeg40to45, bw_cballeg40to45),RooArgList(backgroundYieldeg40to45, bw_cball_yield_eg40to45));
  //fit to data
  totalPDFeg40to45.fitTo(dataeg40to45,Extended());
  //Plot the whole thing
  RooPlot *xframeeg40to45 = xeg40to45.frame(Title("Crystal Ball Signal Plus a RooCMSShape Background -- eg40to45"));
  xframeeg40to45->SetTitle("");
  totalPDFeg40to45.paramOn(xframeeg40to45, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  dataeg40to45.plotOn(xframeeg40to45,LineColor(kBlack));
  totalPDFeg40to45.plotOn(xframeeg40to45, Components(backgroundeg40to45),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDFeg40to45.plotOn(xframeeg40to45, Components(bw_cballeg40to45),LineColor(kGreen)/*,FillColor(kRed),DrawOption("F")*/);
  totalPDFeg40to45.plotOn(xframeeg40to45,LineColor(kBlue));
  dataeg40to45.plotOn(xframeeg40to45,LineColor(kBlack));
  xframeeg40to45->Draw();
  //CBsigEG40to45=cbsigeg40to45.getVal();
  //CBsigEGerror40to45=cbsigeg40to45.getError();
  CBsigEG40to45=bw_cball_yield_eg40to45.getVal();
  CBsigEGerror40to45=bw_cball_yield_eg40to45.getError();
  InvarMassTexteg->Draw();
  Text40to45->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eg40to45.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eg40to45.pdf");
  fout.cd();xframeeg40to45->Write("InvarMassCrystalBall_eg40to45");fin.cd();
  c3->cd(2);xframeeg40to45->Draw();
  c1->cd();

  //Pt45to50
  RooRealVar xeg45to50("x","eg Invariant Mass",61,121,"GeV");
  RooDataHist dataeg45to50("dataeg45to50","dataset",xeg45to50,egInvarMass_Pt45to50);
  //RooCMSShape for background
  RooRealVar aBkgeg45to50("aBkg", "aBkgeg",74/*78., 61.0, 111.0, "GeV"*/); //RooCMSShape (fixed)
  RooRealVar bBkgeg45to50("bBkg", "bBkgeg",.115/*0.2, 0.01, 0.4*/);            //RooCMSShape (fixed)
  RooRealVar ceg45to50("c", "ceg", 0.1, 0.02, 0.3);                   //RooCMSShape (fixed)
  RooRealVar ZMasseg45to50("bkgd peak", "ZMasseg", 91/*91.0, 87.0, 93.5, "GeV"*/);            //Z mass 91.2=pdg
  RooRealVar backgroundYieldeg45to50("bkgd events", "backgroundYieldeg", 1500,100, 100000);
  RooCMSShape backgroundeg45to50("backgroundeg", "backgroundeg", xeg45to50, aBkgeg45to50, bBkgeg45to50, ceg45to50, ZMasseg45to50);
  //Breit-Wigner for peak
  RooRealVar  mReseg45to50("M_{Z^{0}}", "Z0 Resonance  Mass", 90.9, 88.0, 92.0);//,"GeV/c^{2}"); 
  RooRealVar  Gammaeg45to50("#Gamma", "#Gamma", 3.0,.1,5);//,"GeV/c^{2}"); 
  RooBreitWigner bweg45to50("bw","A Breit-Wigner Distribution",xeg45to50,mReseg45to50,Gammaeg45to50);
  //Crystal Ball for resolution
  RooRealVar cbmeaneg45to50("cb_mean", "cbmean" , 0.0,-1.0,2.0) ;
  RooRealVar cbsigmaeg45to50("cb_sigma", "cbsigma" ,1.,-0.1,5.0) ;
  RooRealVar cbsigeg45to50("signal events", "cbsignal", 75000., 50000, 100000);
  RooRealVar neg45to50("n","n", 5,0.,25.);
  RooRealVar alphaeg45to50("alpha","alpha",0.6,0.01,3);
  RooCBShape cballeg45to50("cball", "crystal ball", xeg45to50, cbmeaneg45to50, cbsigmaeg45to50, alphaeg45to50, neg45to50);
  //Convolution of Breit-Wigner and CrystalBall
  RooFFTConvPdf bw_cballeg45to50("bw_cball","Convolution",xeg45to50,bweg45to50,cballeg45to50);
  RooRealVar bw_cball_yield_eg45to50("signal events","signal events",4000,2000,100000); 
  /*   //Crystal Ball
       RooRealVar cbmeaneg45to50("mass", "cbmeaneg" , 91., 87., 95.) ;
       RooRealVar cbsigmaeg45to50("cb_sigma", "cbsigmaeg" , 2.5, 1., 4.) ;
       RooRealVar cbsigeg45to50("signal events", "cbsignaleg", 1400, 90, 4000);
       RooRealVar neg45to50("n","neg", 50.,0.,150.);
       RooRealVar alphaeg45to50("alpha","alphaeg", 1.1,0,2);
       RooCBShape cballeg45to50("cball", "crystal ball eg", xeg45to50, cbmeaneg45to50, cbsigmaeg45to50, alphaeg45to50, neg45to50);*/
  //Addition of crystall ball  and RooCMSShape
  // RooAddPdf totalPDFeg45to50("totalPDFeg45to50","totalPDFeg45to50",RooArgList(backgroundeg45to50, cballeg45to50),RooArgList(backgroundYieldeg45to50, cbsigeg45to50)); 
  RooAddPdf totalPDFeg45to50("totalPDFeg45to50","totalPDFeg45to50",RooArgList(backgroundeg45to50, bw_cballeg45to50),RooArgList(backgroundYieldeg45to50, bw_cball_yield_eg45to50));
  //fit to data
  totalPDFeg45to50.fitTo(dataeg45to50,Extended());
  //Plot the whole thing
  RooPlot *xframeeg45to50 = xeg45to50.frame(Title("Crystal Ball Signal Plus a RooCMSShape Background -- eg45to50"));
  xframeeg45to50->SetTitle("");
  totalPDFeg45to50.paramOn(xframeeg45to50, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  dataeg45to50.plotOn(xframeeg45to50,LineColor(kBlack));
  totalPDFeg45to50.plotOn(xframeeg45to50, Components(backgroundeg45to50),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDFeg45to50.plotOn(xframeeg45to50, Components(bw_cballeg45to50),LineColor(kGreen)/*,FillColor(kRed),DrawOption("F")*/);
  totalPDFeg45to50.plotOn(xframeeg45to50,LineColor(kBlue));
  dataeg45to50.plotOn(xframeeg45to50,LineColor(kBlack));
  xframeeg45to50->Draw();
  //CBsigEG45to50=cbsigeg45to50.getVal();
  //CBsigEGerror45to50=cbsigeg45to50.getError();
  CBsigEG45to50=bw_cball_yield_eg45to50.getVal();
  CBsigEGerror45to50=bw_cball_yield_eg45to50.getError();
  InvarMassTexteg->Draw();
  Text45to50->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eg45to50.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eg45to50.pdf");
  fout.cd();xframeeg45to50->Write("InvarMassCrystalBall_eg45to50");fin.cd();
  c3->cd(3);xframeeg45to50->Draw();
  c1->cd();

  //Pt50to60
  RooRealVar xeg50to60("x","eg Invariant Mass",61,121,"GeV");
  RooDataHist dataeg50to60("dataeg50to60","dataset",xeg50to60,egInvarMass_Pt50to60);
  //RooCMSShape for background
  RooRealVar aBkgeg50to60("aBkg", "aBkgeg",77/*78., 71.0, 111.0, "GeV"*/); //RooCMSShape (fixed)
  RooRealVar bBkgeg50to60("bBkg", "bBkgeg",.09/*0.04, 0.01, 0.11*/);            //RooCMSShape (fixed)
  RooRealVar ceg50to60("c", "ceg", 0.05, 0.02, 0.14);                   //RooCMSShape (fixed)
  RooRealVar ZMasseg50to60("bkgd peak", "ZMasseg",91/*91.0, 85.0, 95.0, "GeV"*/);            //Z mass 91.2=pdg
  RooRealVar backgroundYieldeg50to60("bkgd events", "backgroundYieldeg", 3200,100, 100000);
  RooCMSShape backgroundeg50to60("backgroundeg", "backgroundeg", xeg50to60, aBkgeg50to60, bBkgeg50to60, ceg50to60, ZMasseg50to60);
  //Breit-Wigner for peak
  RooRealVar  mReseg50to60("M_{Z^{0}}", "Z0 Resonance  Mass", 90.9, 88.0, 92.0);//,"GeV/c^{2}"); 
  RooRealVar  Gammaeg50to60("#Gamma", "#Gamma", 4.0,1.,7);//,"GeV/c^{2}"); 
  RooBreitWigner bweg50to60("bw","A Breit-Wigner Distribution",xeg50to60,mReseg50to60,Gammaeg50to60);
  //Crystal Ball for resolution
  RooRealVar cbmeaneg50to60("cb_mean", "cbmean" , 0.0,-1.0,2.0) ;
  RooRealVar cbsigmaeg50to60("cb_sigma", "cbsigma" ,1.,-1,5.0) ;
  RooRealVar cbsigeg50to60("signal events", "cbsignal", 75000., 50000, 100000);
  RooRealVar neg50to60("n","n", 5,0.,25.);
  RooRealVar alphaeg50to60("alpha","alpha",.4,0.01,3);
  RooCBShape cballeg50to60("cball", "crystal ball", xeg50to60, cbmeaneg50to60, cbsigmaeg50to60, alphaeg50to60, neg50to60);
  //Convolution of Breit-Wigner and CrystalBall
  RooFFTConvPdf bw_cballeg50to60("bw_cball","Convolution",xeg50to60,bweg50to60,cballeg50to60);
  RooRealVar bw_cball_yield_eg50to60("signal events","signal events",1500,800,70000); 
  /*    //Crystal Ball
	RooRealVar cbmeaneg50to60("mass", "cbmeaneg" , 91., 87., 95.) ;
	RooRealVar cbsigmaeg50to60("cb_sigma", "cbsigmaeg" , 2.5, 1., 4.) ;
	RooRealVar cbsigeg50to60("signal events", "cbsignaleg", 1400, 500, 40000);
	RooRealVar neg50to60("n","neg", 50.,0.,100.);
	RooRealVar alphaeg50to60("alpha","alphaeg", 2.,0,5.5);
	RooCBShape cballeg50to60("cball", "crystal ball eg", xeg50to60, cbmeaneg50to60, cbsigmaeg50to60, alphaeg50to60, neg50to60);*/
  //Addition of crystall ball  and RooCMSShape
  //RooAddPdf totalPDFeg50to60("totalPDFeg50to60","totalPDFeg50to60",RooArgList(backgroundeg50to60, cballeg50to60),RooArgList(backgroundYieldeg50to60, cbsigeg50to60));
  RooAddPdf totalPDFeg50to60("totalPDFeg50to60","totalPDFeg50to60",RooArgList(backgroundeg50to60, bw_cballeg50to60),RooArgList(backgroundYieldeg50to60, bw_cball_yield_eg50to60));
  //fit to data
  totalPDFeg50to60.fitTo(dataeg50to60,Extended());
  //Plot the whole thing
  RooPlot *xframeeg50to60 = xeg50to60.frame(Title("Crystal Ball Signal Plus a RooCMSShape Background -- eg50to60"));
  xframeeg50to60->SetTitle("");
  totalPDFeg50to60.paramOn(xframeeg50to60, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  dataeg50to60.plotOn(xframeeg50to60,LineColor(kBlack));
  totalPDFeg50to60.plotOn(xframeeg50to60, Components(backgroundeg50to60),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDFeg50to60.plotOn(xframeeg50to60, Components(bw_cballeg50to60),LineColor(kGreen)/*,FillColor(kRed),DrawOption("F")*/);
  totalPDFeg50to60.plotOn(xframeeg50to60,LineColor(kBlue));
  dataeg50to60.plotOn(xframeeg50to60,LineColor(kBlack));
  xframeeg50to60->Draw();
  //CBsigEG50to60=cbsigeg50to60.getVal();
  //CBsigEGerror50to60=cbsigeg50to60.getError();
  CBsigEG50to60=bw_cball_yield_eg50to60.getVal();
  CBsigEGerror50to60=bw_cball_yield_eg50to60.getError();
  InvarMassTexteg->Draw();
  Text50to60->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eg50to60.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eg50to60.pdf");
  fout.cd();xframeeg50to60->Write("InvarMassCrystalBall_eg50to60");fin.cd();
  c3->cd(4);xframeeg50to60->Draw();
  c1->cd();

  //Pt60to80
  RooRealVar xeg60to80("x","eg Invariant Mass",61,121,"GeV");
  RooDataHist dataeg60to80("dataeg60to80","dataset",xeg60to80,egInvarMass_Pt60to80);
  //RooCMSShape for background
  RooRealVar aBkgeg60to80("aBkg", "aBkgeg",95., 71.0, 141.0, "GeV"); //RooCMSShape (fixed)
  RooRealVar bBkgeg60to80("bBkg", "bBkgeg",0.04, -.1, 0.10);            //RooCMSShape (fixed)
  RooRealVar ceg60to80("c", "ceg",0.04, 0.02, 0.13);                   //RooCMSShape (fixed)
  RooRealVar ZMasseg60to80("bkgd peak", "ZMasseg",91/*91.0, 87.0, 93.0, "GeV"*/);            //Z mass 91.2=pdg
  RooRealVar backgroundYieldeg60to80("bkgd events", "backgroundYieldeg", 2700,100, 400000);
  RooCMSShape backgroundeg60to80("backgroundeg", "backgroundeg", xeg60to80, aBkgeg60to80, bBkgeg60to80, ceg60to80, ZMasseg60to80);
  //Breit-Wigner for peak
  RooRealVar  mReseg60to80("M_{Z^{0}}", "Z0 Resonance  Mass", 90.9, 88.0, 92.0);//,"GeV/c^{2}"); 
  RooRealVar  Gammaeg60to80("#Gamma", "#Gamma", 7,5,9);//,"GeV/c^{2}"); 
  RooBreitWigner bweg60to80("bw","A Breit-Wigner Distribution",xeg60to80,mReseg60to80,Gammaeg60to80);
  //Crystal Ball for resolution
  RooRealVar cbmeaneg60to80("cb_mean", "cbmean" , 0.0,-1.0,2.0) ;
  RooRealVar cbsigmaeg60to80("cb_sigma", "cbsigma" ,1.,0,5.0) ;
  RooRealVar cbsigeg60to80("signal events", "cbsignal", 75000., 50000, 100000);
  RooRealVar neg60to80("n","n", 5,0.,25.);
  RooRealVar alphaeg60to80("alpha","alpha",1,0.01,4);
  RooCBShape cballeg60to80("cball", "crystal ball", xeg60to80, cbmeaneg60to80, cbsigmaeg60to80, alphaeg60to80, neg60to80);
  //Convolution of Breit-Wigner and CrystalBall
  RooFFTConvPdf bw_cballeg60to80("bw_cball","Convolution",xeg60to80,bweg60to80,cballeg60to80);
  RooRealVar bw_cball_yield_eg60to80("signal events","signal events",300,200,100000); 
  /*    //Crystal Ball
	RooRealVar cbmeaneg60to80("mass", "cbmeaneg" , 91., 87., 95.) ;
	RooRealVar cbsigmaeg60to80("cb_sigma", "cbsigmaeg" , 2.5, 1., 4.) ;
	RooRealVar cbsigeg60to80("signal events", "cbsignaleg", 1400, 90, 40000);
	RooRealVar neg60to80("n","neg", 50.,0.,100.);
	RooRealVar alphaeg60to80("alpha","alphaeg", 1.1,0,2.2);
	RooCBShape cballeg60to80("cball", "crystal ball eg", xeg60to80, cbmeaneg60to80, cbsigmaeg60to80, alphaeg60to80, neg60to80);*/
  //Addition of crystall ball  and RooCMSShape
  //RooAddPdf totalPDFeg60to80("totalPDFeg60to80","totalPDFeg60to80",RooArgList(backgroundeg60to80, cballeg60to80),RooArgList(backgroundYieldeg60to80, cbsigeg60to80));
  RooAddPdf totalPDFeg60to80("totalPDFeg60to80","totalPDFeg60to80",RooArgList(backgroundeg60to80, bw_cballeg60to80),RooArgList(backgroundYieldeg60to80, bw_cball_yield_eg60to80));
  //fit to data
  totalPDFeg60to80.fitTo(dataeg60to80,Extended());
  //Plot the whole thing
  RooPlot *xframeeg60to80 = xeg60to80.frame(Title("Crystal Ball Signal Plus a RooCMSShape Background -- eg60to80"));
  xframeeg60to80->SetTitle("");
  totalPDFeg60to80.paramOn(xframeeg60to80, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  dataeg60to80.plotOn(xframeeg60to80,LineColor(kBlack));
  totalPDFeg60to80.plotOn(xframeeg60to80, Components(backgroundeg60to80),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDFeg60to80.plotOn(xframeeg60to80, Components(bw_cballeg60to80),LineColor(kGreen)/*,FillColor(kRed),DrawOption("F")*/);
  totalPDFeg60to80.plotOn(xframeeg60to80,LineColor(kBlue));
  dataeg60to80.plotOn(xframeeg60to80,LineColor(kBlack));
  xframeeg60to80->Draw();
  //CBsigEG60to80=cbsigeg60to80.getVal();
  //CBsigEGerror60to80=cbsigeg60to80.getError();
  CBsigEG60to80=bw_cball_yield_eg60to80.getVal();
  CBsigEGerror60to80=bw_cball_yield_eg60to80.getError();
  InvarMassTexteg->Draw();
  Text60to80->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eg60to80.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eg60to80.pdf");
  fout.cd();xframeeg60to80->Write("InvarMassCrystalBall_eg60to80");fin.cd();
  c3->cd(5);xframeeg60to80->Draw();
  c1->cd();

  //Pt80
  RooRealVar xeg80("x","eg Invariant Mass",61,121,"GeV");
  RooDataHist dataeg80("dataeg80","dataset",xeg80,egInvarMass_Pt80);
  //RooCMSShape for background
  RooRealVar aBkgeg80("aBkg", "aBkgeg",121/*100., 71.0, 131.0, "GeV"*/); //RooCMSShape (fixed)
  RooRealVar bBkgeg80("bBkg", "bBkgeg",.017/*0.04, -0.1, 0.3*/);            //RooCMSShape (fixed)
  RooRealVar ceg80("c", "ceg", 0.08, -0.1, 0.15);                   //RooCMSShape (fixed)
  RooRealVar ZMasseg80("bkgd peak", "ZMasseg",91.2/*91.0, 89.0, 93.0, "GeV"*/);            //Z mass 91.2=pdg
  RooRealVar backgroundYieldeg80("bkgd events", "backgroundYieldeg", 700,20, 8000);
  RooCMSShape backgroundeg80("backgroundeg", "backgroundeg", xeg80, aBkgeg80, bBkgeg80, ceg80, ZMasseg80);
  //Breit-Wigner for peak
  RooRealVar  mReseg80("M_{Z^{0}}", "Z0 Resonance  Mass", 90.9, 88.0, 92.0);//,"GeV/c^{2}"); 
  RooRealVar  Gammaeg80("#Gamma", "#Gamma", 2.5,-1,4);//,"GeV/c^{2}"); 
  RooBreitWigner bweg80("bw","A Breit-Wigner Distribution",xeg80,mReseg80,Gammaeg80);
  //Crystal Ball for resolution
  RooRealVar cbmeaneg80("cb_mean", "cbmean" , 0.0,-2,2) ;
  RooRealVar cbsigmaeg80("cb_sigma", "cbsigma" ,1.,-1,5.0) ;
  //RooRealVar cbsigeg80("signal events", "cbsignal", 75000., 50000, 100000);
  RooRealVar neg80("n","n", 5,-10.,50.);
  RooRealVar alphaeg80("alpha","alpha",2,-2,5.5);
  RooCBShape cballeg80("cball", "crystal ball", xeg80, cbmeaneg80, cbsigmaeg80, alphaeg80, neg80);
  //Convolution of Breit-Wigner and CrystalBall
  RooFFTConvPdf bw_cballeg80("bw_cball","Convolution",xeg80,bweg80,cballeg80);
  RooRealVar bw_cball_yield_eg80("signal events","signal events",1000,50,5000); 
  /*    //Crystal Ball
	RooRealVar cbmeaneg80("mass", "cbmeaneg" , 91., 87., 95.) ;
	RooRealVar cbsigmaeg80("cb_sigma", "cbsigmaeg" , 2.5, 1., 4.) ;
	RooRealVar cbsigeg80("signal events", "cbsignaleg", 1400, 90, 40000);
	RooRealVar neg80("n","neg", 50.,-10.,100.);
	RooRealVar alphaeg80("alpha","alphaeg", 1.1,0,3);
	RooCBShape cballeg80("cball", "crystal ball eg", xeg80, cbmeaneg80, cbsigmaeg80, alphaeg80, neg80);*/
  //Addition of crystall ball  and RooCMSShape
  //RooAddPdf totalPDFeg80("totalPDFeg80","totalPDFeg80",RooArgList(backgroundeg80, cballeg80),RooArgList(backgroundYieldeg80, cbsigeg80));
  RooAddPdf totalPDFeg80("totalPDFeg80","totalPDFeg80",RooArgList(backgroundeg80, bw_cballeg80),RooArgList(backgroundYieldeg80, bw_cball_yield_eg80));
  //fit to data
  totalPDFeg80.fitTo(dataeg80,Extended());
  //Plot the whole thing
  RooPlot *xframeeg80 = xeg80.frame(Title("Crystal Ball Signal Plus a RooCMSShape Background -- eg80"));
  xframeeg80->SetTitle("");
  totalPDFeg80.paramOn(xframeeg80, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  dataeg80.plotOn(xframeeg80,LineColor(kBlack));
  totalPDFeg80.plotOn(xframeeg80, Components(backgroundeg80),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDFeg80.plotOn(xframeeg80, Components(bw_cballeg80),LineColor(kGreen)/*,FillColor(kRed),DrawOption("F")*/);
  totalPDFeg80.plotOn(xframeeg80,LineColor(kBlue));
  dataeg80.plotOn(xframeeg80,LineColor(kBlack));
  xframeeg80->Draw();
  //CBsigEG80=cbsigeg80.getVal();
  //CBsigEGerror80=cbsigeg80.getError();
  CBsigEG80=bw_cball_yield_eg80.getVal();
  CBsigEGerror80=bw_cball_yield_eg80.getError();
  InvarMassTexteg->Draw();
  Text80->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eg80.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eg80.pdf");
  fout.cd();xframeeg80->Write("InvarMassCrystalBall_eg80");fin.cd();
  c3->cd(6);xframeeg80->Draw();
  c3->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_egPt.png");
  c1->cd();

  //now gg
  TH1F* ggInvarMassFakeRate = (TH1F*)fin.Get("ggInvarMass");
  RooRealVar xgg("x","gg Invariant Mass",61,121,"GeV");
  RooDataHist datagg("datagg","dataset",xgg,ggInvarMassFakeRate);
  //linear for background
  /*  RooRealVar bkgd_poly_c1("slope of bkgd","slope of background",7,-50,50);
      RooPolynomial bkgd_linear("bkgd_linear", "linear function for background",xgg,RooArgList(bkgd_poly_c1) );
      RooRealVar bkgd_yield("bkg events","bkg events",100,0,1000000);*/

  //RooCMSShape for background
  RooRealVar aBkggg("aBkg", "aBkggg", 68.0, 61.0, 75.0, "GeV"); //RooCMSShape (fixed)
  RooRealVar bBkggg("bBkg", "bBkggg", 0.4, -.1, 1.0);            //RooCMSShape (fixed)
  RooRealVar cgg("c", "cgg", 0.02, -0.01, 0.05);                   //RooCMSShape (fixed)
  RooRealVar ZMassgg("bkgd peak", "ZMassgg", 91/*83.0, 75.0, 89.0, "GeV"*/);            //Z mass 91.2=pdg
  //RooRealVar aBkg("aBkg", "aBkg", 52.9, 40.0, 100.0, "GeV"); //RooCMSShape (fixed)
  //RooRealVar bBkg("bBkg", "bBkg", 0.47, 0.1, 0.6);            //RooCMSShape (fixed)
  //RooRealVar c("c", "c", 1.4, 0.0, 5.0);                   //RooCMSShape (fixed)
  RooRealVar backgroundYieldgg("bkgd events", "backgroundYieldgg", 200000,100000, 300000);
  RooCMSShape backgroundgg("backgroundgg", "backgroundgg", xgg, aBkggg, bBkggg, cgg, ZMassgg);

  //Crystal Ball
  RooRealVar cbmeangg("mass", "cbmeangg" , 88., 78., 95.) ;
  RooRealVar cbsigmagg("cb_sigma", "cbsigmagg" , 6., 2., 10.) ;
  RooRealVar cbsiggg("signal events", "cbsignalgg", 20000, 5000, 60000);
  RooRealVar ngg("n","ngg", 5.,0.,60.);
  RooRealVar alphagg("alpha","alphagg", 2.7/*,1.6,4.0*/);
  RooCBShape cballgg("cballgg", "crystal ball gg", xgg, cbmeangg, cbsigmagg, alphagg, ngg);
  //Addition of crystall ball  and Linear
  RooAddPdf totalPDFgg("totalPDFgg","totalPDFgg",RooArgList(backgroundgg, cballgg),RooArgList(backgroundYieldgg, cbsiggg));
  //fit to data
  totalPDFgg.fitTo(datagg,Extended());
  //Plot the whole thing
  RooPlot *xframegg = xgg.frame(Title("Crystal Ball Signal Plus a RooCMSShape Background -- gg"));
  datagg.plotOn(xframegg,LineColor(kBlack));
  totalPDFgg.plotOn(xframegg, Components(backgroundgg),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDFgg.plotOn(xframegg,LineColor(kBlue));
  //bkgd_linear.plotOn(xframegg,LineColor(kYellow));
  totalPDFgg.paramOn(xframegg, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9),FillStyle(0) );
  datagg.plotOn(xframegg,LineColor(kBlack));
  xframegg->Draw();
  CBsigGG=cbsiggg.getVal();
  CBsigGGerror=cbsiggg.getError();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_gg.png");
  fout.cd();xframegg->Write("InvarMassCrystalBall_gg");fin.cd();
  /*
  //------------------Now One Jet Requirement Invariant Mass fits-----
  RooRealVar x("x","ee Invariant Mass",61,121,"GeV");
  RooDataHist data("data","dataset",x,eeInvarMassFullRange_JetReq);
  //linear for background
  RooRealVar bkgd_poly_c1("slope of bkgd","slope of background",7.5,0,1500);
  RooPolynomial bkgd_linear("bkgd_linear", "linear function for background",x,RooArgList(bkgd_poly_c1) );
  RooRealVar bkgd_yield("bkg events","bkg events",1000,0,1000000);
  //Crystal Ball
  RooRealVar cbmean("mass", "cbmean" , 90, 85, 95) ;
  RooRealVar cbsigma("cb_sigma", "cbsigma" , 4, 1, 10) ;
  RooRealVar cbsig_JetReq("signal events", "cbsignal_JetReq", 10000, 0, 1000000);
  RooRealVar n("n","n", 50,0,1000);
  RooRealVar alpha("alpha","alpha", 2,0,10);
  RooCBShape cball("cball", "crystal ball", x, cbmean, cbsigma, alpha, n);
  //Addition of crystall ball  and Linear
  RooAddPdf totalPDF("totalPDF","totalPDF",RooArgList(bkgd_linear, cball),RooArgList(bkgd_yield, cbsig_JetReq));
  //fit to data
  totalPDF.fitTo(data,Extended());
  //Plot the whole thing
  RooPlot *xframe = x.frame(Title("Crystal Ball Plus a Linear Background - One Loose Jet Required"));
  totalPDF.paramOn(xframe, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  data.plotOn(xframe,LineColor(kBlack));
  totalPDF.plotOn(xframe, Components(bkgd_linear),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDF.plotOn(xframe,LineColor(kBlue));
  xframe->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_ee_JetReq.png");
  fout.cd();xframe->Write("InvarMassCrystalBall_ee_JetReq");fin.cd();

  RooRealVar x("x","eg Invariant Mass",61,121,"GeV");
  RooDataHist data("data","dataset",x,egInvarMass_JetReq);
  //linear for background
  RooRealVar bkgd_poly_c1("slope of bkgd","slope of background",7,0,1500);
  RooPolynomial bkgd_linear("bkgd_linear", "linear function for background",x,RooArgList(bkgd_poly_c1) );
  RooRealVar bkgd_yield("bkg events","bkg events",100,0,1000000);
  //Crystal Ball
  RooRealVar cbmean("mass", "cbmean" , 90, 85, 95) ;
  RooRealVar cbsigma("cb_sigma", "cbsigma" , 4, 1, 10) ;
  RooRealVar cbsigeg_JetReq("signal events", "cbsignal2_JetReq", 1000, 0, 1000000);
  RooRealVar n("n","n", 5,0,1000);
  RooRealVar alpha("alpha","alpha", 2,0,10);
  RooCBShape cball("cball", "crystal ball", x, cbmean, cbsigma, alpha, n);
  //Addition of crystall ball  and Linear
  RooAddPdf totalPDF("totalPDF","totalPDF",RooArgList(bkgd_linear, cball),RooArgList(bkgd_yield, cbsigeg_JetReq));
  //fit to data
  totalPDF.fitTo(data,Extended());
  //Plot the whole thing
  RooPlot *xframe = x.frame(Title("Crystal Ball Plus a Linear Background - One Loose Jet Required"));
  totalPDF.paramOn(xframe, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  data.plotOn(xframe,LineColor(kBlack));
  totalPDF.plotOn(xframe, Components(bkgd_linear),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDF.plotOn(xframe,LineColor(kBlue));
  //bkgd_linear.plotOn(xframe,LineColor(kYellow));
  xframe->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_eg_JetReq.png");
  fout.cd();xframe->Write("InvarMassCrystalBall_eg_JetReq");fin.cd();

  RooRealVar x("x","gg Invariant Mass",61,121,"GeV");
  RooDataHist data("data","dataset",x,ggInvarMass_JetReq);
  //linear for background
  RooRealVar bkgd_poly_c1("slope of bkgd","slope of background",7,-50,50);
  RooPolynomial bkgd_linear("bkgd_linear", "linear function for background",x,RooArgList(bkgd_poly_c1) );
  RooRealVar bkgd_yield("bkg events","bkg events",100,0,1000000);
  //Crystal Ball
  RooRealVar cbmean("mass", "cbmean" , 90, 85, 95) ;
  RooRealVar cbsigma("cb_sigma", "cbsigma" , 4, 1, 10) ;
  RooRealVar cbsiggg_JetReq("signal events", "cbsignal3_JetReq", 1000, 0, 1000000);
  RooRealVar n("n","n", 5,-50,50);
  RooRealVar alpha("alpha","alpha", 2,0,10);
  RooCBShape cball("cball", "crystal ball", x, cbmean, cbsigma, alpha, n);
  //Addition of crystall ball  and Linear
  RooAddPdf totalPDF("totalPDF","totalPDF",RooArgList(bkgd_linear, cball),RooArgList(bkgd_yield, cbsiggg_JetReq));
  //fit to data
  totalPDF.fitTo(data,Extended());
  //Plot the whole thing
  RooPlot *xframe = x.frame(Title("Crystal Ball Plus a Linear Background - One Loose Jet Required"));
  totalPDF.paramOn(xframe, Format("NE",AutoPrecision(1)),Layout(0.6,0.98,0.9) );
  data.plotOn(xframe,LineColor(kBlack));
  totalPDF.plotOn(xframe, Components(bkgd_linear),LineColor(kRed),FillColor(kRed),DrawOption("F"));
  totalPDF.plotOn(xframe,LineColor(kBlue));
  //bkgd_linear.plotOn(xframe,LineColor(kYellow));
  xframe->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMassCrystalBall_gg_JetReq.png");
  fout.cd();xframe->Write("InvarMassCrystalBall_gg_JetReq");fin.cd();
  */

  GetFakeRateAndErrorFullRange(CBsigEE,CBsigEEerror,CBsigEG,CBsigEGerror,FakeRate,FakeRateErr,FakeRateErrForLimits);
   
  //----------------------end of Invariant Mass Fit plots


  TH1F* eeMetFromZZ=(TH1F*)f_ZZ.Get("eeMet_reweightJet_binned");
  TH1F* eeMetSBhighFromZZ=(TH1F*)f_ZZ.Get("eeSidebandHighMet_reweightJet_binned");
  TH1F* eeMetSBlowFromZZ=(TH1F*)f_ZZ.Get("eeSidebandLowMet_reweightJet_binned");
  TH1F* eeMetFromWZ=(TH1F*)f_WZ.Get("eeMet_reweightJet_binned");
  TH1F* eeMetSBhighFromWZ=(TH1F*)f_WZ.Get("eeSidebandHighMet_reweightJet_binned");
  TH1F* eeMetSBlowFromWZ=(TH1F*)f_WZ.Get("eeSidebandLowMet_reweightJet_binned");
  
  TH1F* eeMetFromZZ_JetReq=(TH1F*)f_ZZ.Get("eeMet_reweightJet_binned_JetReq");
  TH1F* eeMetSBhighFromZZ_JetReq=(TH1F*)f_ZZ.Get("eeSidebandHighMet_reweightJet_binned_JetReq");
  TH1F* eeMetSBlowFromZZ_JetReq=(TH1F*)f_ZZ.Get("eeSidebandLowMet_reweightJet_binned_JetReq");
  TH1F* eeMetFromWZ_JetReq=(TH1F*)f_WZ.Get("eeMet_reweightJet_binned_JetReq");
  TH1F* eeMetSBhighFromWZ_JetReq=(TH1F*)f_WZ.Get("eeSidebandHighMet_reweightJet_binned_JetReq");
  TH1F* eeMetSBlowFromWZ_JetReq=(TH1F*)f_WZ.Get("eeSidebandLowMet_reweightJet_binned_JetReq");
  
  TH1F* eeMetFromZZ_2JetReq=(TH1F*)f_ZZ.Get("eeMet_reweightJet_binned_JetReq");//temp
  TH1F* eeMetSBhighFromZZ_2JetReq=(TH1F*)f_ZZ.Get("eeSidebandHighMet_reweightJet_binned_JetReq");//temp
  TH1F* eeMetSBlowFromZZ_2JetReq=(TH1F*)f_ZZ.Get("eeSidebandLowMet_reweightJet_binned_JetReq");//temp
  TH1F* eeMetFromWZ_2JetReq=(TH1F*)f_WZ.Get("eeMet_reweightJet_binned_JetReq");//temp
  TH1F* eeMetSBhighFromWZ_2JetReq=(TH1F*)f_WZ.Get("eeSidebandHighMet_reweightJet_binned_JetReq");//temp
  TH1F* eeMetSBlowFromWZ_2JetReq=(TH1F*)f_WZ.Get("eeSidebandLowMet_reweightJet_binned_JetReq");//temp

  //(lumi * LOxsec * (kfactor) / nEvents)
  //float ZZscale = L_int*0.231*1.359/500000.;
  //float WZscale = L_int*0.412*1.618/(2956080*(500000/497920));//don't know how many events, have to estimate
  float ZZscale = L_int*0.1913*1.359/500000.;
  float WZscale = L_int*0.412*1.618/3000000.;//from prep
 
  c1->SetLogy(1);

  eeMetFromZZ->Draw();
  eeMetSBlowFromZZ->SetLineColor(kRed);
  eeMetSBhighFromZZ->SetLineColor(kBlue);
  eeMetSBlowFromZZ->Draw("SAME");
  eeMetSBhighFromZZ->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/eeMetfromZZwithSB.png");
  eeMetFromZZ_JetReq->Draw();
  eeMetSBlowFromZZ_JetReq->SetLineColor(kRed);
  eeMetSBhighFromZZ_JetReq->SetLineColor(kBlue);
  eeMetSBlowFromZZ_JetReq->Draw("SAME");
  eeMetSBhighFromZZ_JetReq->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/eeMetfromZZwithSB_JetReq.png");
  eeMetFromWZ->Draw();
  eeMetSBlowFromWZ->SetLineColor(kRed);
  eeMetSBhighFromWZ->SetLineColor(kBlue);
  eeMetSBlowFromWZ->Draw("SAME");
  eeMetSBhighFromWZ->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/eeMetfromWZwithSB.png");
  eeMetFromWZ_JetReq->Draw();
  eeMetSBlowFromWZ_JetReq->SetLineColor(kRed);
  eeMetSBhighFromWZ_JetReq->SetLineColor(kBlue);
  eeMetSBlowFromWZ_JetReq->Draw("SAME");
  eeMetSBhighFromWZ_JetReq->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/eeMetfromWZwithSB_JetReq.png");


  eeMetFromZZ->Add(eeMetSBhighFromZZ,-1);
  eeMetFromZZ->Add(eeMetSBlowFromZZ,-1);
  eeMetFromWZ->Add(eeMetSBhighFromWZ,-1);
  eeMetFromWZ->Add(eeMetSBlowFromWZ,-1);
  eeMetFromZZ_JetReq->Add(eeMetSBhighFromZZ_JetReq,-1);
  eeMetFromZZ_JetReq->Add(eeMetSBlowFromZZ_JetReq,-1);
  eeMetFromWZ_JetReq->Add(eeMetSBhighFromWZ_JetReq,-1);
  eeMetFromWZ_JetReq->Add(eeMetSBlowFromWZ_JetReq,-1);


  eeMetFromZZ->Scale(ZZscale);
  eeMetFromWZ->Scale(WZscale);
  eeMetFromZZ_JetReq->Scale(ZZscale);
  eeMetFromWZ_JetReq->Scale(WZscale);


  float eeMetRawFromZZ100up=eeMetFromZZ->Integral(eeMetFromZZ->FindBin(100.1),-1);
  float eeMetRawFromWZ100up=eeMetFromWZ->Integral(eeMetFromWZ->FindBin(100.1),-1);
  float eeMetRawFromZZ100up_JetReq=eeMetFromZZ_JetReq->Integral(eeMetFromZZ_JetReq->FindBin(100.1),-1);
  float eeMetRawFromWZ100up_JetReq=eeMetFromWZ_JetReq->Integral(eeMetFromWZ_JetReq->FindBin(100.1),-1);


  eeMetFromZZ->SetLineColor(kGreen);
  eeMetFromZZ->SetMarkerColor(kGreen);
  eeMetFromWZ->SetLineColor(kRed);
  eeMetFromWZ->SetMarkerColor(kRed);
  TH1F* Zz=(TH1F*)eeMetFromZZ->Rebin(NmetBins,"Zz",xbins);
  TH1F* Wz=(TH1F*)eeMetFromWZ->Rebin(NmetBins,"Wz",xbins);
  Wz->Draw();
  Zz->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/eeMetZZandWZbgSubtracted.png");
  fin.cd();



    
  //Plots that need log scale here:
  //need to do log scale DiEMPt here so they don't get screwed up by ratio plots
  c1->SetLogy(1);
  TH1F* ggDiEMPt=(TH1F*)fin.Get("ggDiEMPt");ggDiEMPt->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiEMPt_gg.png");
  fout.cd();ggDiEMPt->Write("DiEMPt_gg");fin.cd();
  TH1F* egDiEMPt=(TH1F*)fin.Get("egDiEMPt");egDiEMPt->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiEMPt_eg.png");
  fout.cd();egDiEMPt->Write("DiEMPt_eg");fin.cd();
  TH1F* eeDiEMPt=(TH1F*)fin.Get("eeDiEMPt");eeDiEMPt->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiEMPt_ee.png");
  fout.cd();eeDiEMPt->Write("DiEMPt_ee");fin.cd();
  TH1F* ffDiEMPt=(TH1F*)fin.Get("ffDiEMPt");ffDiEMPt->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiEMPt_ff.png");
  fout.cd();ffDiEMPt->Write("DiEMPt_ff");fin.cd();
  TH1F* ggDiEMPt_JetReq=(TH1F*)fin.Get("ggDiEMPt_JetReq");ggDiEMPt_JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiEMPt_gg_JetReq.png");
  fout.cd();ggDiEMPt_JetReq->Write("DiEMPt_gg_JetReq");fin.cd();
  TH1F* egDiEMPt_JetReq=(TH1F*)fin.Get("egDiEMPt_JetReq");egDiEMPt_JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiEMPt_eg_JetReq.png");
  fout.cd();egDiEMPt_JetReq->Write("DiEMPt_eg_JetReq");fin.cd();
  TH1F* eeDiEMPt_JetReq=(TH1F*)fin.Get("eeDiEMPt_JetReq");eeDiEMPt_JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiEMPt_ee_JetReq.png");
  fout.cd();eeDiEMPt_JetReq->Write("DiEMPt_ee_JetReq");fin.cd();
  TH1F* ffDiEMPt_JetReq=(TH1F*)fin.Get("ffDiEMPt_JetReq");ffDiEMPt_JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiEMPt_ff_JetReq.png");
  fout.cd();ffDiEMPt_JetReq->Write("DiEMPt_ff_JetReq");fin.cd();


  gStyle->SetOptStat(0);

  TH1F* eeDiJetPt = (TH1F*)fin.Get("eeDiJetPt");eeDiJetPt ->Sumw2();
  TH1F* eeDiJetPt_reweight_binned = (TH1F*)fin.Get("eeDiJetPt_reweight_binned");eeDiJetPt_reweight_binned->Sumw2();
  TH1F* ggDiJetPt = (TH1F*)fin.Get("ggDiJetPt");ggDiJetPt->Sumw2();
  //TH1F* ggDiJetPt_reweight_binned = (TH1F*)fin.Get("ggDiJetPt_reweight_binned");
  TH1F* ffDiJetPt = (TH1F*)fin.Get("ffDiJetPt");ffDiJetPt->Sumw2();
  TH1F* ffDiJetPt_reweight_binned = (TH1F*)fin.Get("ffDiJetPt_reweight_binned");ffDiJetPt_reweight_binned->Sumw2();
  TH1F* eeDiJetPt_JetReq = (TH1F*)fin.Get("eeDiJetPt_JetReq");eeDiJetPt_JetReq->Sumw2();
  TH1F* eeDiJetPt_reweight_binned_JetReq = (TH1F*)fin.Get("eeDiJetPt_reweight_binned_JetReq");eeDiJetPt_reweight_binned_JetReq->Sumw2();
  TH1F* ggDiJetPt_JetReq = (TH1F*)fin.Get("ggDiJetPt_JetReq");ggDiJetPt_JetReq->Sumw2();
  //TH1F* ggDiJetPt_reweight_binned_JetReq = (TH1F*)fin.Get("ggDiJetPt_reweightJet_binned_JetReq");
  TH1F* ffDiJetPt_JetReq = (TH1F*)fin.Get("ffDiJetPt_JetReq");ffDiJetPt_JetReq->Sumw2();
  TH1F* ffDiJetPt_reweight_binned_JetReq = (TH1F*)fin.Get("ffDiJetPt_reweight_binned_JetReq");ffDiJetPt_reweight_binned_JetReq->Sumw2();

  TH1F* ggDiJetPt_0Jet = (TH1F*)fin.Get("ggDiJetPt_0Jet");ggDiJetPt_0Jet->Sumw2();
  TH1F* ggDiJetPt_1Jet = (TH1F*)fin.Get("ggDiJetPt_1Jet");ggDiJetPt_1Jet->Sumw2();
  TH1F* ggDiJetPt_2Jet = (TH1F*)fin.Get("ggDiJetPt_2Jet");ggDiJetPt_2Jet->Sumw2();
  TH1F* ffDiJetPt_reweight_binned_0Jet = (TH1F*)fin.Get("ffDiJetPt_reweight_binned_0Jet");ffDiJetPt_reweight_binned_0Jet->Sumw2();
  TH1F* ffDiJetPt_reweight_binned_1Jet = (TH1F*)fin.Get("ffDiJetPt_reweight_binned_1Jet");ffDiJetPt_reweight_binned_1Jet->Sumw2();
  TH1F* ffDiJetPt_reweight_binned_2Jet = (TH1F*)fin.Get("ffDiJetPt_reweight_binned_2Jet");ffDiJetPt_reweight_binned_2Jet->Sumw2();
  TH1F* eeDiJetPt_reweight_binned_0Jet = (TH1F*)fin.Get("eeDiJetPt_reweight_binned_0Jet");eeDiJetPt_reweight_binned_0Jet->Sumw2();
  TH1F* eeDiJetPt_reweight_binned_1Jet = (TH1F*)fin.Get("eeDiJetPt_reweight_binned_1Jet");eeDiJetPt_reweight_binned_1Jet->Sumw2();
  TH1F* eeDiJetPt_reweight_binned_2Jet = (TH1F*)fin.Get("eeDiJetPt_reweight_binned_2Jet");eeDiJetPt_reweight_binned_2Jet->Sumw2();

  /*
    int numBins=15;
    Double_t DiEMPtBins[16]={0.,5.,10.,15.,20.,30.,40.,50.,60.,70.,80.,100.,120.,150.,200.,600.,};
    TH1F* eeRebin = (TH1F*)eeDiJetPt->Rebin(numBins,"eeRebin",DiEMPtBins);
    TH1F* ggRebin = (TH1F*)ggDiJetPt->Rebin(numBins,"ggRebin",DiEMPtBins);
    TH1F* ffRebin = (TH1F*)ffDiJetPt->Rebin(numBins,"ffRebin",DiEMPtBins);
    TH1F* eeRebin_JetReq = (TH1F*)eeDiJetPt_JetReq->Rebin(numBins,"eeRebin_JetReq",DiEMPtBins);
    TH1F* ggRebin_JetReq = (TH1F*)ggDiJetPt_JetReq->Rebin(numBins,"ggRebin_JetReq",DiEMPtBins);
    TH1F* ffRebin_JetReq = (TH1F*)ffDiJetPt_JetReq->Rebin(numBins,"ffRebin_JetReq",DiEMPtBins);
    TH1F* eeReweightRebin = (TH1F*)eeDiJetPt_reweight_binned->Rebin(numBins,"eeReweightRebin",DiEMPtBins);
    //TH1F* ggReweightRebin = (TH1F*)ggDiJetPt_reweight_binned->Rebin(numBins,"ggReweightRebin",DiEMPtBins);
    TH1F* ffReweightRebin = (TH1F*)ffDiJetPt_reweight_binned->Rebin(numBins,"ffReweightRebin",DiEMPtBins);
    TH1F* eeReweightRebin_JetReq = (TH1F*)eeDiJetPt_reweight_binned_JetReq->Rebin(numBins,"eeReweightRebin_JetReq",DiEMPtBins);
    //TH1F* ggReweightRebin_JetReq = (TH1F*)ggDiJetPt_reweight_binned_JetReq->Rebin(numBins,"ggReweightRebin_JetReq",DiEMPtBins);
    TH1F* ffReweightRebin_JetReq = (TH1F*)ffDiJetPt_reweight_binned_JetReq->Rebin(numBins,"ffReweightRebin_JetReq",DiEMPtBins);
  */
  c1->SetLogy(1);

  float dijet = 1./eeDiJetPt->Integral();eeDiJetPt->Scale(dijet);
  eeDiJetPt->GetYaxis()->SetTitle("Normalized Number of Events");
  eeDiJetPt->GetXaxis()->SetTitle("diJetPt (GeV)");
  eeDiJetPt->SetTitle("");
  ggDiJetPt->GetXaxis()->SetRangeUser(0,400);
  dijet = 1./ggDiJetPt->Integral();ggDiJetPt->Scale(dijet);
  dijet = 1./ggDiJetPt_0Jet->Integral();ggDiJetPt_0Jet->Scale(dijet);
  dijet = 1./ggDiJetPt_1Jet->Integral();ggDiJetPt_1Jet->Scale(dijet);
  dijet = 1./ggDiJetPt_2Jet->Integral();ggDiJetPt_2Jet->Scale(dijet);
  ggDiJetPt->Draw("histo");
  eeDiJetPt->SetLineColor(kRed);
  eeDiJetPt->Draw("histosames");
  dijet = 1./ffDiJetPt->Integral();ffDiJetPt->Scale(dijet);
  ffDiJetPt->SetLineColor(kBlue);
  ffDiJetPt->Draw("histoSAMEs");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt.png");
  float dijetMax=ggDiJetPt->GetMaximum(),dijetMin=ggDiJetPt->GetMinimum();
  dijet = 1./eeDiJetPt_JetReq->Integral();eeDiJetPt_JetReq->Scale(dijet);
  eeDiJetPt_JetReq->GetYaxis()->SetTitle("Normalized Number of Events");
  eeDiJetPt_JetReq->GetXaxis()->SetTitle("diJetPt (GeV)");
  eeDiJetPt_JetReq->SetTitle("");
  ggDiJetPt_JetReq->GetXaxis()->SetRangeUser(0,400);
  dijet = 1./ggDiJetPt_JetReq->Integral();ggDiJetPt_JetReq->Scale(dijet);
  ggDiJetPt_JetReq->Draw("histo");
  eeDiJetPt_JetReq->SetLineColor(kRed);
  eeDiJetPt_JetReq->Draw("histosames");
  dijet = 1./ffDiJetPt_JetReq->Integral();ffDiJetPt_JetReq->Scale(dijet);
  ffDiJetPt_JetReq->SetLineColor(kBlue);
  ffDiJetPt_JetReq->Draw("histoSAMEs");
  Text_Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_JetReq.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_JetReq.pdf");
  float dijetMax_JetReq=ggDiJetPt_JetReq->GetMaximum(),dijetMin_JetReq=ggDiJetPt_JetReq->GetMinimum();
  dijet = 1./eeDiJetPt_reweight_binned->Integral();eeDiJetPt_reweight_binned->Scale(dijet);
  eeDiJetPt_reweight_binned->GetYaxis()->SetTitle("Normalized Number of Events");
  eeDiJetPt_reweight_binned->GetXaxis()->SetTitle("diJetPt (GeV)");
  eeDiJetPt_reweight_binned->SetTitle("");
  eeDiJetPt_reweight_binned->SetLineColor(kRed);eeDiJetPt_reweight_binned->SetMarkerColor(kRed);
  dijet = 1./ffDiJetPt_reweight_binned->Integral();ffDiJetPt_reweight_binned->Scale(dijet);
  ffDiJetPt_reweight_binned->SetLineColor(kBlue);ffDiJetPt_reweight_binned->SetMarkerColor(kBlue);ffDiJetPt_reweight_binned->SetMarkerSize(.6);
  dijet = 1./ffDiJetPt_reweight_binned_0Jet->Integral();ffDiJetPt_reweight_binned_0Jet->Scale(dijet);
  ffDiJetPt_reweight_binned_0Jet->SetLineColor(kBlue);ffDiJetPt_reweight_binned_0Jet->SetMarkerColor(kBlue);ffDiJetPt_reweight_binned_0Jet->SetMarkerSize(.6);
  dijet = 1./ffDiJetPt_reweight_binned_1Jet->Integral();ffDiJetPt_reweight_binned_1Jet->Scale(dijet);
  ffDiJetPt_reweight_binned_1Jet->SetLineColor(kBlue);ffDiJetPt_reweight_binned_1Jet->SetMarkerColor(kBlue);ffDiJetPt_reweight_binned_1Jet->SetMarkerSize(.6);
  dijet = 1./ffDiJetPt_reweight_binned_2Jet->Integral();ffDiJetPt_reweight_binned_2Jet->Scale(dijet);
  ffDiJetPt_reweight_binned_2Jet->SetLineColor(kBlue);ffDiJetPt_reweight_binned_2Jet->SetMarkerColor(kBlue);ffDiJetPt_reweight_binned_2Jet->SetMarkerSize(.6);
  dijet = 1./eeDiJetPt_reweight_binned_0Jet->Integral();eeDiJetPt_reweight_binned_0Jet->Scale(dijet);
  eeDiJetPt_reweight_binned_0Jet->SetLineColor(kRed);eeDiJetPt_reweight_binned_0Jet->SetMarkerColor(kRed);eeDiJetPt_reweight_binned_0Jet->SetMarkerSize(.6);
  dijet = 1./eeDiJetPt_reweight_binned_1Jet->Integral();eeDiJetPt_reweight_binned_1Jet->Scale(dijet);
  eeDiJetPt_reweight_binned_1Jet->SetLineColor(kRed);eeDiJetPt_reweight_binned_1Jet->SetMarkerColor(kRed);eeDiJetPt_reweight_binned_1Jet->SetMarkerSize(.6);
  dijet = 1./eeDiJetPt_reweight_binned_2Jet->Integral();eeDiJetPt_reweight_binned_2Jet->Scale(dijet);
  eeDiJetPt_reweight_binned_2Jet->SetLineColor(kRed);eeDiJetPt_reweight_binned_2Jet->SetMarkerColor(kRed);eeDiJetPt_reweight_binned_2Jet->SetMarkerSize(.6);
  ggDiJetPt->Draw("histo");
  ffDiJetPt_reweight_binned->Draw("histosames");
  eeDiJetPt_reweight_binned->Draw("histosames");
  //dijet = 1./ggDiJetPt->Integral();ggDiJetPt->Scale(dijet);
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_reweight.png");
  dijet = 1./eeDiJetPt_reweight_binned_JetReq->Integral();eeDiJetPt_reweight_binned_JetReq->Scale(dijet);
  eeDiJetPt_reweight_binned_JetReq->GetYaxis()->SetTitle("Normalized Number of Events");
  eeDiJetPt_reweight_binned_JetReq->GetXaxis()->SetTitle("diJetPt (GeV)");
  eeDiJetPt_reweight_binned_JetReq->SetTitle("");
  eeDiJetPt_reweight_binned_JetReq->SetLineColor(kRed);eeDiJetPt_reweight_binned_JetReq->SetMarkerColor(kRed);
  dijet = 1./ffDiJetPt_reweight_binned_JetReq->Integral();ffDiJetPt_reweight_binned_JetReq->Scale(dijet);
  ffDiJetPt_reweight_binned_JetReq->SetLineColor(kBlue);ffDiJetPt_reweight_binned_JetReq->SetMarkerColor(kBlue);ffDiJetPt_reweight_binned_JetReq->SetMarkerSize(.6);
  ggDiJetPt_JetReq->Draw("histo");
  ffDiJetPt_reweight_binned_JetReq->Draw("histosames");
  eeDiJetPt_reweight_binned_JetReq->Draw("histosames");
  //dijet = 1./ggDiJetPt_JetReq->Integral();ggDiJetPt_JetReq->Scale(dijet);
  Text_Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_reweight_JetReq.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_reweight_JetReq.pdf");

  ggDiJetPt->GetXaxis()->SetRangeUser(0,600);
  ggDiJetPt_JetReq->GetXaxis()->SetRangeUser(0,600);

  TH1F* ffDiJetPt_reweight_binned2 = (TH1F*)ffDiJetPt_reweight_binned->Rebin(numBins,"ffDiJetPt_reweight_binned2",DiEMPtBins);
  TH1F* ffDiJetPt_reweight_binned2_JetReq = (TH1F*)ffDiJetPt_reweight_binned_JetReq->Rebin(numBins,"ffDiJetPt_reweight_binned2_JetReq",DiEMPtBins);
  TH1F* eeDiJetPt_reweight_binned2 = (TH1F*)eeDiJetPt_reweight_binned->Rebin(numBins,"eeDiJetPt_reweight_binned2",DiEMPtBins);
  TH1F* eeDiJetPt_reweight_binned2_JetReq = (TH1F*)eeDiJetPt_reweight_binned_JetReq->Rebin(numBins,"eeDiJetPt_reweight_binned2_JetReq",DiEMPtBins);
  TH1F* ggDiJetPt2 = (TH1F*)ggDiJetPt->Rebin(numBins,"ggDiJetPt2",DiEMPtBins);
  TH1F* ggDiJetPt2_JetReq = (TH1F*)ggDiJetPt_JetReq->Rebin(numBins,"ggDiJetPt2_JetReq",DiEMPtBins);

  TH1F* ggDiJetPt2_0Jet = (TH1F*)ggDiJetPt_0Jet->Rebin(numBins,"ggDiJetPt2_0Jet",DiEMPtBins);
  TH1F* ggDiJetPt2_1Jet = (TH1F*)ggDiJetPt_1Jet->Rebin(numBins,"ggDiJetPt2_1Jet",DiEMPtBins);
  TH1F* ggDiJetPt2_2Jet = (TH1F*)ggDiJetPt_2Jet->Rebin(numBins,"ggDiJetPt2_2Jet",DiEMPtBins);
  TH1F* ffDiJetPt_reweight_binned2_0Jet = (TH1F*)ffDiJetPt_reweight_binned_0Jet->Rebin(numBins,"ffDiJetPt_reweight_binned2_0Jet",DiEMPtBins);
  TH1F* ffDiJetPt_reweight_binned2_1Jet = (TH1F*)ffDiJetPt_reweight_binned_1Jet->Rebin(numBins,"ffDiJetPt_reweight_binned2_1Jet",DiEMPtBins);
  TH1F* ffDiJetPt_reweight_binned2_2Jet = (TH1F*)ffDiJetPt_reweight_binned_2Jet->Rebin(numBins,"ffDiJetPt_reweight_binned2_2Jet",DiEMPtBins);
  TH1F* eeDiJetPt_reweight_binned2_0Jet = (TH1F*)eeDiJetPt_reweight_binned_0Jet->Rebin(numBins,"eeDiJetPt_reweight_binned2_0Jet",DiEMPtBins);
  TH1F* eeDiJetPt_reweight_binned2_1Jet = (TH1F*)eeDiJetPt_reweight_binned_1Jet->Rebin(numBins,"eeDiJetPt_reweight_binned2_1Jet",DiEMPtBins);
  TH1F* eeDiJetPt_reweight_binned2_2Jet = (TH1F*)eeDiJetPt_reweight_binned_2Jet->Rebin(numBins,"eeDiJetPt_reweight_binned2_2Jet",DiEMPtBins);


  for(int i=1;i<ffDiJetPt_reweight_binned2->GetNbinsX()+1;i++){
    ggDiJetPt2->SetBinContent(i,ggDiJetPt2->GetBinContent(i)/ggDiJetPt2->GetBinWidth(i));
    ggDiJetPt2_JetReq->SetBinContent(i,ggDiJetPt2_JetReq->GetBinContent(i)/ggDiJetPt2_JetReq->GetBinWidth(i));
    ffDiJetPt_reweight_binned2->SetBinContent(i,ffDiJetPt_reweight_binned2->GetBinContent(i)/ffDiJetPt_reweight_binned2->GetBinWidth(i));
    ffDiJetPt_reweight_binned2_JetReq->SetBinContent(i,ffDiJetPt_reweight_binned2_JetReq->GetBinContent(i)/ffDiJetPt_reweight_binned2_JetReq->GetBinWidth(i));
    eeDiJetPt_reweight_binned2->SetBinContent(i,eeDiJetPt_reweight_binned2->GetBinContent(i)/eeDiJetPt_reweight_binned2->GetBinWidth(i));
    eeDiJetPt_reweight_binned2_JetReq->SetBinContent(i,eeDiJetPt_reweight_binned2_JetReq->GetBinContent(i)/eeDiJetPt_reweight_binned2_JetReq->GetBinWidth(i));
    ggDiJetPt2->SetBinError(i,ggDiJetPt2->GetBinError(i)/ggDiJetPt2->GetBinWidth(i));
    ggDiJetPt2_JetReq->SetBinError(i,ggDiJetPt2_JetReq->GetBinError(i)/ggDiJetPt2_JetReq->GetBinWidth(i));
    ffDiJetPt_reweight_binned2->SetBinError(i,ffDiJetPt_reweight_binned2->GetBinError(i)/ffDiJetPt_reweight_binned2->GetBinWidth(i));
    ffDiJetPt_reweight_binned2_JetReq->SetBinError(i,ffDiJetPt_reweight_binned2_JetReq->GetBinError(i)/ffDiJetPt_reweight_binned2_JetReq->GetBinWidth(i));
    eeDiJetPt_reweight_binned2->SetBinError(i,eeDiJetPt_reweight_binned2->GetBinError(i)/eeDiJetPt_reweight_binned2->GetBinWidth(i));
    eeDiJetPt_reweight_binned2_JetReq->SetBinError(i,eeDiJetPt_reweight_binned2_JetReq->GetBinError(i)/eeDiJetPt_reweight_binned2_JetReq->GetBinWidth(i));
 
    ggDiJetPt2_0Jet->SetBinContent(i,ggDiJetPt2_0Jet->GetBinContent(i)/ggDiJetPt2_0Jet->GetBinWidth(i));
    ggDiJetPt2_1Jet->SetBinContent(i,ggDiJetPt2_1Jet->GetBinContent(i)/ggDiJetPt2_1Jet->GetBinWidth(i));
    ggDiJetPt2_2Jet->SetBinContent(i,ggDiJetPt2_2Jet->GetBinContent(i)/ggDiJetPt2_2Jet->GetBinWidth(i));
    ffDiJetPt_reweight_binned2_0Jet->SetBinContent(i,ffDiJetPt_reweight_binned2_0Jet->GetBinContent(i)/ffDiJetPt_reweight_binned2_0Jet->GetBinWidth(i));
    ffDiJetPt_reweight_binned2_1Jet->SetBinContent(i,ffDiJetPt_reweight_binned2_1Jet->GetBinContent(i)/ffDiJetPt_reweight_binned2_1Jet->GetBinWidth(i));
    ffDiJetPt_reweight_binned2_2Jet->SetBinContent(i,ffDiJetPt_reweight_binned2_2Jet->GetBinContent(i)/ffDiJetPt_reweight_binned2_2Jet->GetBinWidth(i));
    eeDiJetPt_reweight_binned2_0Jet->SetBinContent(i,eeDiJetPt_reweight_binned2_0Jet->GetBinContent(i)/eeDiJetPt_reweight_binned2_0Jet->GetBinWidth(i));
    eeDiJetPt_reweight_binned2_1Jet->SetBinContent(i,eeDiJetPt_reweight_binned2_1Jet->GetBinContent(i)/eeDiJetPt_reweight_binned2_1Jet->GetBinWidth(i));
    eeDiJetPt_reweight_binned2_2Jet->SetBinContent(i,eeDiJetPt_reweight_binned2_2Jet->GetBinContent(i)/eeDiJetPt_reweight_binned2_2Jet->GetBinWidth(i));

  }

  ggDiJetPt2->GetXaxis()->SetRangeUser(0,399);
  ggDiJetPt2->GetYaxis()->SetTitle("Events/GeV");
  ggDiJetPt2->Draw("histo");
  ffDiJetPt_reweight_binned2->Draw("pesames");
  eeDiJetPt_reweight_binned2->Draw("pesames");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_reweight_rebin.png");
  ggDiJetPt2_JetReq->GetXaxis()->SetRangeUser(0,399);
  ggDiJetPt2_JetReq->GetYaxis()->SetTitle("Events/GeV");
  ggDiJetPt2_JetReq->Draw("histo");
  ffDiJetPt_reweight_binned2_JetReq->Draw("pesames");
  eeDiJetPt_reweight_binned2_JetReq->Draw("pesames");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_reweight_rebin_JetReq.png");

  ggDiJetPt2_0Jet->GetXaxis()->SetRangeUser(0,399);
  ggDiJetPt2_0Jet->GetYaxis()->SetTitle("Events/GeV");
  ggDiJetPt2_0Jet->Draw("histo");
  ffDiJetPt_reweight_binned2_0Jet->Draw("pesames");
  eeDiJetPt_reweight_binned2_0Jet->Draw("pesames");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_reweight_rebin_0Jet.png");
  ggDiJetPt2_1Jet->GetXaxis()->SetRangeUser(0,399);
  ggDiJetPt2_1Jet->GetYaxis()->SetTitle("Events/GeV");
  ggDiJetPt2_1Jet->Draw("histo");
  ffDiJetPt_reweight_binned2_1Jet->Draw("pesames");
  eeDiJetPt_reweight_binned2_1Jet->Draw("pesames");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_reweight_rebin_1Jet.png");
  ggDiJetPt2_2Jet->GetXaxis()->SetRangeUser(0,399);
  ggDiJetPt2_2Jet->GetYaxis()->SetTitle("Events/GeV");
  ggDiJetPt2_2Jet->Draw("histo");
  ffDiJetPt_reweight_binned2_2Jet->Draw("pesames");
  eeDiJetPt_reweight_binned2_2Jet->Draw("pesames");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_reweight_rebin_2Jet.png");
 

  c1->SetLogy(0);
  TLine l_diJet1(0,1,400,1);
  ffDiJetPt_reweight_binned2->Divide(ggDiJetPt2);
  ffDiJetPt_reweight_binned2->GetXaxis()->SetRangeUser(0,399);
  ffDiJetPt_reweight_binned2->GetYaxis()->SetRangeUser(.4,1.6);
  ffDiJetPt_reweight_binned2->Draw("PE");l_diJet1.Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_ffOVERgg.png");
  ffDiJetPt_reweight_binned2_JetReq->Divide(ggDiJetPt2_JetReq);
  ffDiJetPt_reweight_binned2_JetReq->GetXaxis()->SetRangeUser(0,399);
  ffDiJetPt_reweight_binned2_JetReq->GetYaxis()->SetRangeUser(.4,1.6);
  ffDiJetPt_reweight_binned2_JetReq->Draw("PE");l_diJet1.Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_ffOVERgg_JetReq.png");


  c1->SetLogy(0);
  ffDiJetPt_reweight_binned2_0Jet->Divide(ggDiJetPt2_0Jet);
  ffDiJetPt_reweight_binned2_0Jet->GetYaxis()->SetRangeUser(.4,1.6);
  ffDiJetPt_reweight_binned2_0Jet->GetXaxis()->SetRangeUser(0,399);
  ffDiJetPt_reweight_binned2_0Jet->Draw("PE");l_diJet1.Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_ffOVERgg_0Jet.png");
  ffDiJetPt_reweight_binned2_1Jet->Divide(ggDiJetPt2_1Jet);
  ffDiJetPt_reweight_binned2_1Jet->GetYaxis()->SetRangeUser(.4,1.6);
  ffDiJetPt_reweight_binned2_1Jet->GetXaxis()->SetRangeUser(0,399);
  ffDiJetPt_reweight_binned2_1Jet->Draw("PE");l_diJet1.Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_ffOVERgg_1Jet.png");
  ffDiJetPt_reweight_binned2_2Jet->Divide(ggDiJetPt2_2Jet);
  ffDiJetPt_reweight_binned2_2Jet->GetYaxis()->SetRangeUser(.4,1.6);
  ffDiJetPt_reweight_binned2_2Jet->GetXaxis()->SetRangeUser(0,399);
  ffDiJetPt_reweight_binned2_2Jet->Draw("PE");l_diJet1.Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_ffOVERgg_2Jet.png");
  


  eeDiJetPt_reweight_binned2->Divide(ggDiJetPt2);
  eeDiJetPt_reweight_binned2->GetYaxis()->SetRangeUser(.4,1.6);
  eeDiJetPt_reweight_binned2->Draw("PE");l_diJet1.Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_eeOVERgg.png");
  eeDiJetPt_reweight_binned2_JetReq->Divide(ggDiJetPt2_JetReq);
  eeDiJetPt_reweight_binned2_JetReq->GetYaxis()->SetRangeUser(.4,1.6);
  eeDiJetPt_reweight_binned2_JetReq->Draw("PE");l_diJet1.Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPt_eeOVERgg_JetReq.png");
  c1->SetLogy(1);

  TH1F* eeDiJetPt_0Jet=(TH1F*)fin.Get("eeDiJetPt_0Jet");
  dijet = 1./eeDiJetPt_0Jet->Integral();eeDiJetPt_0Jet->Scale(dijet);
  /*eeDiJetPt_0Jet->GetYaxis()->SetTitle("Normalized Number of Events");
    eeDiJetPt_0Jet->GetXaxis()->SetTitle("diJetPt (GeV)");
    eeDiJetPt_0Jet->SetTitle("");
    eeDiJetPt_0Jet->GetXaxis()->SetRangeUser(0,220);
    eeDiJetPt_0Jet->SetLineColor(kRed);
    eeDiJetPt_0Jet->Draw("histo");*/
  TH1F* gammafakeDiJetPt_0Jet=(TH1F*)fin.Get("gammafakeDiJetPt_0Jet");
  dijet = 1./gammafakeDiJetPt_0Jet->Integral();gammafakeDiJetPt_0Jet->Scale(dijet);
  gammafakeDiJetPt_0Jet->GetYaxis()->SetTitle("Normalized Number of Events");
  gammafakeDiJetPt_0Jet->GetXaxis()->SetTitle("diJetPt (GeV)");
  gammafakeDiJetPt_0Jet->SetTitle("");
  gammafakeDiJetPt_0Jet->GetXaxis()->SetRangeUser(0,220);
  gammafakeDiJetPt_0Jet->SetLineColor(kRed);
  gammafakeDiJetPt_0Jet->Draw("histo");
  TH1F* ggDiJetPt_0Jetnew=(TH1F*)fin.Get("ggDiJetPt_0Jet");
  dijet = 1./ggDiJetPt_0Jetnew->Integral();ggDiJetPt_0Jetnew->Scale(dijet);/*ggDiJetPt_0Jetnew->SetLineColor(kRed);*/ggDiJetPt_0Jetnew->Draw("histoSAME");
  TH1F* ffDiJetPt_0Jet=(TH1F*)fin.Get("ffDiJetPt_0Jet");
  dijet = 1./ffDiJetPt_0Jet->Integral();ffDiJetPt_0Jet->Scale(dijet);ffDiJetPt_0Jet->SetLineColor(kBlue);ffDiJetPt_0Jet->Draw("histoSAME");
  Text_0Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure14_DiJetPt_0Jet.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure14_DiJetPt_0Jet.pdf");
  TH1F* ggDiJetPt_1Jetnew=(TH1F*)fin.Get("ggDiJetPt_1Jet");
  dijet = 1./ggDiJetPt_1Jetnew->Integral();ggDiJetPt_1Jetnew->Scale(dijet);
  ggDiJetPt_1Jetnew->GetYaxis()->SetTitle("Normalized Number of Events");
  ggDiJetPt_1Jetnew->GetXaxis()->SetTitle("diJetPt (GeV)");
  ggDiJetPt_1Jetnew->SetTitle("");
  ggDiJetPt_1Jetnew->GetXaxis()->SetRangeUser(0,220);
  ggDiJetPt_1Jetnew->Draw("histo");
  TH1F* eeDiJetPt_1Jet=(TH1F*)fin.Get("eeDiJetPt_1Jet");
  //dijet = 1./eeDiJetPt_1Jet->Integral();eeDiJetPt_1Jet->Scale(dijet);eeDiJetPt_1Jet->SetLineColor(kRed);eeDiJetPt_1Jet->Draw("histoSAME");
  TH1F* gammafakeDiJetPt_1Jet=(TH1F*)fin.Get("gammafakeDiJetPt_1Jet");
  dijet = 1./gammafakeDiJetPt_1Jet->Integral();gammafakeDiJetPt_1Jet->Scale(dijet);gammafakeDiJetPt_1Jet->SetLineColor(kRed);gammafakeDiJetPt_1Jet->Draw("histoSAME");
  TH1F* ffDiJetPt_1Jet=(TH1F*)fin.Get("ffDiJetPt_1Jet");
  dijet = 1./ffDiJetPt_1Jet->Integral();ffDiJetPt_1Jet->Scale(dijet);ffDiJetPt_1Jet->SetLineColor(kBlue);ffDiJetPt_1Jet->Draw("histoSAME");
  Text_1Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure14_DiJetPt_1Jet.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure14_DiJetPt_1Jet.pdf");
  TH1F* ggDiJetPt_2Jetnew=(TH1F*)fin.Get("ggDiJetPt_2Jet");
  dijet = 1./ggDiJetPt_2Jetnew->Integral();ggDiJetPt_2Jetnew->Scale(dijet);
  ggDiJetPt_2Jetnew->GetYaxis()->SetTitle("Normalized Number of Events");
  ggDiJetPt_2Jetnew->GetXaxis()->SetTitle("diJetPt (GeV)");
  ggDiJetPt_2Jetnew->SetTitle("");
  ggDiJetPt_2Jetnew->GetXaxis()->SetRangeUser(0,220);
  ggDiJetPt_2Jetnew->Draw("histo");
  TH1F* eeDiJetPt_2Jet=(TH1F*)fin.Get("eeDiJetPt_2Jet");
  //dijet = 1./eeDiJetPt_2Jet->Integral();eeDiJetPt_2Jet->Scale(dijet);eeDiJetPt_2Jet->SetLineColor(kRed);eeDiJetPt_2Jet->Draw("histoSAME");
  TH1F* gammafakeDiJetPt_2Jet=(TH1F*)fin.Get("gammafakeDiJetPt_2Jet");
  dijet = 1./gammafakeDiJetPt_2Jet->Integral();gammafakeDiJetPt_2Jet->Scale(dijet);gammafakeDiJetPt_2Jet->SetLineColor(kRed);gammafakeDiJetPt_2Jet->Draw("histoSAME");
  TH1F* ffDiJetPt_2Jet=(TH1F*)fin.Get("ffDiJetPt_2Jet");
  dijet = 1./ffDiJetPt_2Jet->Integral();ffDiJetPt_2Jet->Scale(dijet);ffDiJetPt_2Jet->SetLineColor(kBlue);ffDiJetPt_2Jet->Draw("histoSAME");
  Text_2Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure14_DiJetPt_2Jet.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure14_DiJetPt_2Jet.pdf");
  //fin.cd()

  TH1F* ggMetData=(TH1F*)fin.Get("ggMet");c1->cd();ggMetData->Draw();
  //cout<<"Look Here"<<hex<<ggMetData<<dec<<fin.GetName()<<endl;
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_gg.png");
  fout.cd();ggMetData->Write("Met_gg");fin.cd();

  //TH1F* egMet=(TH1F*)fin.Get("egMet");egMet->Draw();egMet->Scale(0.0209);
  TH1F* egMet=(TH1F*)fin.Get("egMet_reweight_Pt_Nvtx_Eta");egMet->Draw();//egMet->Scale(0.01);
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_eg.png");
  fout.cd();egMet->Write("Met_eg");fin.cd();
  TH1F* eeMet=(TH1F*)fin.Get("eeMet");eeMet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_ee.png");
  fout.cd();eeMet->Write("Met_ee");fin.cd();
  TH1F* eeSideBandMet=(TH1F*)fin.Get("eeSideBandMet");eeSideBandMet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_eeSideband.png");
  fout.cd();eeSideBandMet->Write("Met_eeSideband");fin.cd();
  TH1F* ffMet=(TH1F*)fin.Get("ffMet");ffMet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_ff.png");
  fout.cd();ffMet->Write("Met_ff");fin.cd();
  TH1F* ggMet_JetReq=(TH1F*)fin.Get("ggMet_JetReq");ggMet_JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_gg_JetReq.png");
  fout.cd();ggMet_JetReq->Write("Met_gg_JetReq");fin.cd();
  TH1F* ggMet_2JetReq=(TH1F*)fin.Get("ggMet_2JetReq");ggMet_2JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_gg_2JetReq.png");
  fout.cd();ggMet_2JetReq->Write("Met_gg_2JetReq");fin.cd();
  //TH1F* egMet_JetReq=(TH1F*)fin.Get("egMet_JetReq");egMet_JetReq->Draw();egMet_JetReq->Scale(0.0209);
  TH1F* egMet_JetReq=(TH1F*)fin.Get("egMet_reweight_Pt_Nvtx_Eta_JetReq");egMet_JetReq->Draw();//egMet_JetReq->Scale(0.01);
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_eg_JetReq.png");
  fout.cd();egMet_JetReq->Write("Met_eg_JetReq");fin.cd();
  //TH1F* egMet_2JetReq=(TH1F*)fin.Get("egMet_2JetReq");egMet_2JetReq->Draw();egMet_2JetReq->Scale(0.0209);
  TH1F* egMet_2JetReq=(TH1F*)fin.Get("egMet_reweight_Pt_Nvtx_Eta_2JetReq");egMet_2JetReq->Draw();//egMet_2JetReq->Scale(0.01);
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_eg_2JetReq.png");
  fout.cd();egMet_2JetReq->Write("Met_eg_2JetReq");fin.cd();
  TH1F* eeMet_JetReq=(TH1F*)fin.Get("eeMet_JetReq");eeMet_JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_ee_JetReq.png");
  fout.cd();eeMet_JetReq->Write("Met_ee_JetReq");fin.cd();
  TH1F* eeSideBandMet_JetReq=(TH1F*)fin.Get("eeSideBandMet_JetReq");eeSideBandMet_JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_eeSideband_JetReq.png");
  fout.cd();eeSideBandMet_JetReq->Write("Met_eeSideband_JetReq");fin.cd();
  TH1F* ffMet_JetReq=(TH1F*)fin.Get("ffMet_JetReq");ffMet_JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_ff_JetReq.png");
  fout.cd();ffMet_JetReq->Write("Met_ff_JetReq");fin.cd();
  TH1F* eeMet_2JetReq=(TH1F*)fin.Get("eeMet_2JetReq");eeMet_2JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_ee_2JetReq.png");
  fout.cd();eeMet_2JetReq->Write("Met_ee_2JetReq");fin.cd();
  TH1F* eeSideBandMet_2JetReq=(TH1F*)fin.Get("eeSideBandMet_2JetReq");eeSideBandMet_2JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_eeSideband_2JetReq.png");
  fout.cd();eeSideBandMet_2JetReq->Write("Met_eeSideband_2JetReq");fin.cd();
  TH1F* ffMet_2JetReq=(TH1F*)fin.Get("ffMet_2JetReq");ffMet_2JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_ff_2JetReq.png");
  fout.cd();ffMet_2JetReq->Write("Met_ff_2JetReq");fin.cd();
  
  TH1F* egMet_reweight_Pt=(TH1F*)fin.Get("egMet_reweight_Pt");
  TH1F* egMet_reweight_Nvtx=(TH1F*)fin.Get("egMet_reweight_Nvtx");
  TH1F* egMet_reweight_Pt_Nvtx_Eta=(TH1F*)egMet->Clone();//(TH1F*)fin.Get("egMet_reweight_Pt_Nvtx_Eta");//egMet_reweight_Pt_Nvtx_Eta->Scale(0.01);//changes whether you use egMet or egMet_reweight
  float egxmax=250.;
  egMet_reweight_Pt->GetXaxis()->SetRangeUser(0,egxmax);egMet_reweight_Nvtx->GetXaxis()->SetRangeUser(0,egxmax);egMet_reweight_Pt_Nvtx_Eta->GetXaxis()->SetRangeUser(0,egxmax);

  egMet_reweight_Pt->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_eg_reweight_Pt.png");
   egMet_reweight_Nvtx->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_eg_reweight_Nvtx.png");
   egMet_reweight_Pt_Nvtx_Eta->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_eg_reweight_Pt_Nvtx.png");
  c1->SetLogy(0);  
  TH1F* egnew = (TH1F*)egMet->Rebin(NmetBins,"egnew",xbins);//egnew->Scale(0.0209);
  TH1F* egptnew = (TH1F*)egMet_reweight_Pt->Rebin(NmetBins,"egptnew",xbins);egptnew->SetLineColor(kBlue);egptnew->SetMarkerColor(kBlue);
  TH1F* egnvtxnew = (TH1F*)egMet_reweight_Nvtx->Rebin(NmetBins,"egnvtxnew",xbins);egnvtxnew->SetLineColor(kGreen);egnvtxnew->SetMarkerColor(kGreen);
  TH1F* egptnvtxnew = (TH1F*)egMet_reweight_Pt_Nvtx_Eta->Rebin(NmetBins,"egptnvtxnew",xbins);egptnvtxnew->SetLineColor(kRed);egptnvtxnew->SetMarkerColor(kRed);
  c1->SetLogy(1);
  float ptnvtxscale = egnew->Integral()/egptnvtxnew->Integral();egptnvtxnew->Scale(ptnvtxscale);
  egnew->GetYaxis()->SetRangeUser(1.5e1,1e3);
  egnew->Draw();//egptnew->Draw("SAMES");egnvtxnew->Draw("SAMES");
  egptnvtxnew->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_eg_and_reweight.png");
  c1->SetLogy(0);

  TLine egline(0,1,250,1);
  //TH1F* egMetRatio=(TH1F*)egMet->Rebin(NmetBins,"egMetRatio",xbins);
  //egMetRatio->Scale(0.0199);
  egptnew->Divide(egnew);
  egptnew->Draw();egline.Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_eg_reweight_Pt_over_eg.png");
  //egMetRatio=(TH1F*)egMet->Rebin(NmetBins,"egMetRatio",xbins);
  //egMetRatio->Scale(0.0199);
  egnvtxnew->Divide(egnew);
  egnvtxnew->Draw();egline.Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_eg_reweight_Nvtx_over_eg.png");
  //egMetRatio=(TH1F*)egMet->Rebin(NmetBins,"egMetRatio",xbins);
  //egMetRatio->Scale(0.0199);
  egptnvtxnew->Divide(egnew);
  egptnvtxnew->Draw();egline.Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_eg_reweight_Pt_Nvtx_over_eg.png");
  
  c1->SetLogy(1);


  //Now do reweighted MET plots
 
  
  /* TH1F *ggSig720 = fSig720.Get("ggMet");
     TH1F *ggSig800 = fSig800.Get("ggMet");
     ggSig720->Rebin(NmetBins,"ggSig720New",xbins);
     ggSig800->Rebin(NmetBins,"ggSig800New",xbins);
     float Sig720Scale = ggSig720New->Integral()/ggSig720New->Integral();
     float Sig800Scale = ggSig800New->Integral()/ggSig800New->Integral();
     ggSig720New->Scale(Sig720Scale);
     ggSig800New->Scale(Sig800Scale);*/
  TH1F *ggSig = (TH1F*)fSig.Get("ggMet");
  TH1F *ffSig = (TH1F*)fSig.Get("ffMet");
  //TH1F *ffSig = (TH1F*)fSig.Get("gfMet");
  TH1F *ggSig_JetReq = (TH1F*)fSig.Get("ggMet_JetReq");
  TH1F *ffSig_JetReq = (TH1F*)fSig.Get("ffMet_JetReq");
  //TH1F *ffSig_JetReq = (TH1F*)fSig.Get("gfMet_JetReq");
  TH1F *ggSig2012 = (TH1F*)fSig2012.Get("ggMet");
  TH1F *ggSig2012_2 = (TH1F*)fSig2012_2.Get("ggMet");
  TH1F *ggSig2012_JetReq = (TH1F*)fSig2012.Get("ggMet_JetReq");
  TH1F *ggSig2012_JetReq_2 = (TH1F*)fSig2012_2.Get("ggMet_JetReq");
  TH1F *ggSig2012_2JetReq = (TH1F*)fSig2012.Get("ggMet_2JetReq");
  TH1F *ggSig2012_2JetReq_2 = (TH1F*)fSig2012_2.Get("ggMet_2JetReq");
  /* fLimitsSig.cd();
     ggSig->Write("h_gg_met_nojet_mS2000_mG1015_mN375");ffSig->Write("h_ff_met_nojet_mS2000_mG1015_mN375");
     ggSig_JetReq->Write("h_gg_met_1jet_mS2000_mG1015_mN375");ffSig_JetReq->Write("h_ff_met_1jet_mS2000_mG1015_mN375");
     fin.cd();*/
  //apply scale factor and trigger efficiency and number of events and luminosity and cross section to signal MET to over-plot
  float scaleFactorSquared=1.005*1.005;
  float trigEff=.883;
  float nSignalEvents=10000.;
  //TH1F *rhoSig = (TH1F*)fSig.Get("rho");
  TH1F* ggSigNew = (TH1F*)ggSig->Rebin(NmetBins,"ggSigNew",xbins);
  ggSigNew->Scale(0.034*L_int/10000);//get the right x-sec
  TH1F* ggSig2012Rebin = (TH1F*)ggSig2012->Rebin(NmetBins,"ggSig2012Rebin",xbins);
  //include overflow in last bin
  Double_t sigOverFlow=0.,sigOverFlowErr=0.;int overflowBinsig=ggSig2012->FindBin(ggSig2012Rebin->GetBinLowEdge(NmetBins+1)+.01);
  sigOverFlow=ggSig2012->IntegralAndError(overflowBinsig,999,sigOverFlowErr);
  float sigLastBinNew=ggSig2012Rebin->GetBinContent(NmetBins)+sigOverFlow;
  float sigLastBinNewErr=sqrt(ggSig2012Rebin->GetBinError(NmetBins)*ggSig2012Rebin->GetBinError(NmetBins)+sigOverFlowErr*sigOverFlowErr);
  ggSig2012Rebin->SetBinContent(NmetBins,sigLastBinNew);
  ggSig2012Rebin->SetBinError(NmetBins,sigLastBinNewErr);
  //have to use this tempScale for now - otherwise no MET plot - need to understand
  float tempScale = ggSig2012Rebin->GetEntries()/ggSig2012Rebin->Integral();ggSig2012Rebin->Scale(tempScale);
  ggSig2012Rebin->Scale(scaleFactorSquared*trigEff*.0006669*L_int/nSignalEvents);
  TH1F* ggSig2012Rebin_2 = (TH1F*)ggSig2012_2->Rebin(NmetBins,"ggSig2012Rebin_2",xbins);
  //include overflow in last bin
  sigOverFlow=0.,sigOverFlowErr=0.;overflowBinsig=ggSig2012_2->FindBin(ggSig2012Rebin_2->GetBinLowEdge(NmetBins+1)+.01);
  sigOverFlow=ggSig2012_2->IntegralAndError(overflowBinsig,999,sigOverFlowErr);
  sigLastBinNew=ggSig2012Rebin_2->GetBinContent(NmetBins)+sigOverFlow;
  sigLastBinNewErr=sqrt(ggSig2012Rebin_2->GetBinError(NmetBins)*ggSig2012Rebin_2->GetBinError(NmetBins)+sigOverFlowErr*sigOverFlowErr);
  ggSig2012Rebin_2->SetBinContent(NmetBins,sigLastBinNew);
  ggSig2012Rebin_2->SetBinError(NmetBins,sigLastBinNewErr);
  tempScale = ggSig2012Rebin_2->GetEntries()/ggSig2012Rebin_2->Integral();ggSig2012Rebin_2->Scale(tempScale);
  ggSig2012Rebin_2->Scale(scaleFactorSquared*trigEff*.1151*L_int/nSignalEvents);
  TH1F* ggSig2012Rebin_JetReq = (TH1F*)ggSig2012_JetReq->Rebin(NmetBins,"ggSig2012Rebin_JetReq",xbins);
  //include overflow in last bin
  sigOverFlow=0.,sigOverFlowErr=0.;overflowBinsig=ggSig2012_JetReq->FindBin(ggSig2012Rebin_JetReq->GetBinLowEdge(NmetBins+1)+.01);
  sigOverFlow=ggSig2012_JetReq->IntegralAndError(overflowBinsig,999,sigOverFlowErr);
  sigLastBinNew=ggSig2012Rebin_JetReq->GetBinContent(NmetBins)+sigOverFlow;
  sigLastBinNewErr=sqrt(ggSig2012Rebin_JetReq->GetBinError(NmetBins)*ggSig2012Rebin_JetReq->GetBinError(NmetBins)+sigOverFlowErr*sigOverFlowErr);
  ggSig2012Rebin_JetReq->SetBinContent(NmetBins,sigLastBinNew);
  ggSig2012Rebin_JetReq->SetBinError(NmetBins,sigLastBinNewErr);
  tempScale = ggSig2012Rebin_JetReq->GetEntries()/ggSig2012Rebin_JetReq->Integral();ggSig2012Rebin_JetReq->Scale(tempScale);
  ggSig2012Rebin_JetReq->Scale(scaleFactorSquared*trigEff*.0006669*L_int/nSignalEvents);
  TH1F* ggSig2012Rebin_JetReq_2 = (TH1F*)ggSig2012_JetReq_2->Rebin(NmetBins,"ggSig2012Rebin_JetReq_2",xbins);
  //include overflow in last bin
  sigOverFlow=0.,sigOverFlowErr=0.;overflowBinsig=ggSig2012_JetReq_2->FindBin(ggSig2012Rebin_JetReq_2->GetBinLowEdge(NmetBins+1)+.01);
  sigOverFlow=ggSig2012_JetReq_2->IntegralAndError(overflowBinsig,999,sigOverFlowErr);
  sigLastBinNew=ggSig2012Rebin_JetReq_2->GetBinContent(NmetBins)+sigOverFlow;
  sigLastBinNewErr=sqrt(ggSig2012Rebin_JetReq_2->GetBinError(NmetBins)*ggSig2012Rebin_JetReq_2->GetBinError(NmetBins)+sigOverFlowErr*sigOverFlowErr);
  ggSig2012Rebin_JetReq_2->SetBinContent(NmetBins,sigLastBinNew);
  ggSig2012Rebin_JetReq_2->SetBinError(NmetBins,sigLastBinNewErr);
  tempScale = ggSig2012Rebin_JetReq_2->GetEntries()/ggSig2012Rebin_JetReq_2->Integral();ggSig2012Rebin_JetReq_2->Scale(tempScale);
  ggSig2012Rebin_JetReq_2->Scale(scaleFactorSquared*trigEff*.1151*L_int/nSignalEvents);
  TH1F* ggSig2012Rebin_2JetReq = (TH1F*)ggSig2012_2JetReq->Rebin(NmetBins,"ggSig2012Rebin_2JetReq",xbins);
  //include overflow in last bin
  sigOverFlow=0.,sigOverFlowErr=0.;overflowBinsig=ggSig2012_2JetReq->FindBin(ggSig2012Rebin_2JetReq->GetBinLowEdge(NmetBins+1)+.01);
  sigOverFlow=ggSig2012_2JetReq->IntegralAndError(overflowBinsig,999,sigOverFlowErr);
  sigLastBinNew=ggSig2012Rebin_2JetReq->GetBinContent(NmetBins)+sigOverFlow;
  sigLastBinNewErr=sqrt(ggSig2012Rebin_2JetReq->GetBinError(NmetBins)*ggSig2012Rebin_2JetReq->GetBinError(NmetBins)+sigOverFlowErr*sigOverFlowErr);
  ggSig2012Rebin_2JetReq->SetBinContent(NmetBins,sigLastBinNew);
  ggSig2012Rebin_2JetReq->SetBinError(NmetBins,sigLastBinNewErr);
  tempScale = ggSig2012Rebin_2JetReq->GetEntries()/ggSig2012Rebin_2JetReq->Integral();ggSig2012Rebin_2JetReq->Scale(tempScale);
  ggSig2012Rebin_2JetReq->Scale(scaleFactorSquared*trigEff*.0006669*L_int/nSignalEvents);
  TH1F* ggSig2012Rebin_2JetReq_2 = (TH1F*)ggSig2012_2JetReq_2->Rebin(NmetBins,"ggSig2012Rebin_2JetReq_2",xbins);
  //include overflow in last bin
  sigOverFlow=0.,sigOverFlowErr=0.;overflowBinsig=ggSig2012_2JetReq_2->FindBin(ggSig2012Rebin_2JetReq_2->GetBinLowEdge(NmetBins+1)+.01);
  sigOverFlow=ggSig2012_2JetReq_2->IntegralAndError(overflowBinsig,999,sigOverFlowErr);
  sigLastBinNew=ggSig2012Rebin_2JetReq_2->GetBinContent(NmetBins)+sigOverFlow;
  sigLastBinNewErr=sqrt(ggSig2012Rebin_2JetReq_2->GetBinError(NmetBins)*ggSig2012Rebin_2JetReq_2->GetBinError(NmetBins)+sigOverFlowErr*sigOverFlowErr);
  ggSig2012Rebin_2JetReq_2->SetBinContent(NmetBins,sigLastBinNew);
  ggSig2012Rebin_2JetReq_2->SetBinError(NmetBins,sigLastBinNewErr);
  tempScale = ggSig2012Rebin_2JetReq_2->GetEntries()/ggSig2012Rebin_2JetReq_2->Integral();ggSig2012Rebin_2JetReq_2->Scale(tempScale);
  ggSig2012Rebin_JetReq_2->Scale(scaleFactorSquared*trigEff*.1151*L_int/nSignalEvents);

  ggSig2012Rebin->SetLineColor(kGreen);ggSig2012Rebin->SetLineWidth(2);
  ggSig2012Rebin_2->SetLineColor(kViolet);ggSig2012Rebin_2->SetLineWidth(2);
  ggSig2012Rebin_JetReq->SetLineColor(kGreen);ggSig2012Rebin_JetReq->SetLineWidth(2);
  ggSig2012Rebin_JetReq_2->SetLineColor(kViolet);ggSig2012Rebin_JetReq_2->SetLineWidth(2);
  ggSig2012Rebin_2JetReq->SetLineColor(kGreen);ggSig2012Rebin_2JetReq->SetLineWidth(2);
  ggSig2012Rebin_2JetReq_2->SetLineColor(kViolet);ggSig2012Rebin_2JetReq_2->SetLineWidth(2);

  //-------------First ee------------   
  TString title="Toys/eeMet_toy_1";
  TString titleSBlow="Toys/eeSidebandLowMet_toy_1";TString titleSBhigh="Toys/eeSidebandHighMet_toy_1";
  TH1F* hist=(TH1F*)fin.Get(title);hist->Sumw2();
  TH1F* histSBlow=(TH1F*)fin.Get(titleSBlow);histSBlow->Sumw2();
  TH1F* histSBhigh=(TH1F*)fin.Get(titleSBhigh);histSBhigh->Sumw2();
    
  const int HistBins=hist->GetNbinsX();

  float StdDev[HistBins]/*={0.}*/,StdDevSBLow[HistBins]/*={0.}*/,StdDevSBHigh[HistBins]/*={0.}*/;
  float mean[HistBins]/*={0.}*/,meanSBlow[HistBins]/*={0.}*/,meanSBhigh[HistBins]/*={0.}*/;

  for(int i=0;i<HistBins;i++){
    StdDev[i]=0.;StdDevSBLow[i]=0.;StdDevSBHigh[i]=0.;mean[i]=0.;meanSBlow[i]=0.;meanSBhigh[i]=0.;
  }

  for(int bin=1;bin<HistBins+1;bin++){
    for(int i=1;i<1001;i++){
      title="Toys/eeMet_toy_";title+=i;
      titleSBlow="Toys/eeSidebandLowMet_toy_";titleSBlow+=i;
      titleSBhigh="Toys/eeSidebandHighMet_toy_";titleSBhigh+=i;
      hist=(TH1F*)fin.Get(title);if(hist->GetSumw2()==0)hist->Sumw2();
      histSBlow=(TH1F*)fin.Get(titleSBlow);if(histSBlow->GetSumw2()==0)histSBlow->Sumw2();
      histSBhigh=(TH1F*)fin.Get(titleSBhigh);if(histSBhigh->GetSumw2()==0)histSBhigh->Sumw2();
      mean[bin-1]+=hist->GetBinContent(bin);
      //cout<<hist->GetBinContent(bin)<<endl;
      meanSBlow[bin-1]+=histSBlow->GetBinContent(bin);
      meanSBhigh[bin-1]+=histSBhigh->GetBinContent(bin);
    }
    //cout<<"Bin:"<<bin<<"  MeanSum:"<<mean[bin-1]<<"  Mean:"<<mean[bin-1]/1000<<endl;
    mean[bin-1]=mean[bin-1]/1000;
    meanSBlow[bin-1]=meanSBlow[bin-1]/1000;
    meanSBhigh[bin-1]=meanSBhigh[bin-1]/1000;
  }
    
  float StdDev_JetReq[HistBins]/*={0.}*/,StdDevSBLow_JetReq[HistBins]/*={0.}*/,StdDevSBHigh_JetReq[HistBins]/*={0.}*/;
  float mean_JetReq[HistBins]/*={0.}*/,meanSBlow_JetReq[HistBins]/*={0.}*/,meanSBhigh_JetReq[HistBins]/*={0.}*/;
    
  for(int i=0;i<HistBins;i++){
    StdDev_JetReq[i]=0.;StdDevSBLow_JetReq[i]=0.;StdDevSBHigh_JetReq[i]=0.;mean_JetReq[i]=0.;meanSBlow_JetReq[i]=0.;meanSBhigh_JetReq[i]=0.;
  }

  TString title_JetReq="Toys/ee_JetReqMet_toy_1";
  TString titleSBlow_JetReq="Toys/eeSideband_JetReqLowMet_toy_1";TString titleSBhigh_JetReq="Toys/eeSideband_JetReqHighMet_toy_1";
  TH1F* hist_JetReq=(TH1F*)fin.Get(title_JetReq);hist_JetReq->Sumw2();
  TH1F* histSBlow_JetReq=(TH1F*)fin.Get(titleSBlow_JetReq);histSBlow_JetReq->Sumw2();
  TH1F* histSBhigh_JetReq=(TH1F*)fin.Get(titleSBhigh_JetReq);histSBhigh_JetReq->Sumw2();
    
  for(int bin=1;bin<HistBins+1;bin++){
    for(int i=1;i<1001;i++){
      title_JetReq="Toys/ee_JetReqMet_toy_";title_JetReq+=i;
      titleSBlow_JetReq="Toys/eeSideband_JetReqLowMet_toy_";titleSBlow_JetReq+=i;
      titleSBhigh_JetReq="Toys/eeSideband_JetReqHighMet_toy_";titleSBhigh_JetReq+=i;
      hist_JetReq=(TH1F*)fin.Get(title_JetReq);if(hist_JetReq->GetSumw2()==0)hist_JetReq->Sumw2();
      histSBlow_JetReq=(TH1F*)fin.Get(titleSBlow_JetReq);if(histSBlow_JetReq->GetSumw2()==0)histSBlow_JetReq->Sumw2();
      histSBhigh_JetReq=(TH1F*)fin.Get(titleSBhigh_JetReq);if(histSBhigh_JetReq->GetSumw2()==0)histSBhigh_JetReq->Sumw2();
      mean_JetReq[bin-1]+=hist_JetReq->GetBinContent(bin);
      //cout<<hist->GetBinContent(bin)<<endl;
      meanSBlow_JetReq[bin-1]+=histSBlow_JetReq->GetBinContent(bin);
      meanSBhigh_JetReq[bin-1]+=histSBhigh_JetReq->GetBinContent(bin);
    }
    //cout<<"Bin:"<<bin<<"  MeanSum:"<<mean_JetReq[bin-1]<<"  Mean:"<<mean_JetReq[bin-1]/1000<<endl;
    mean_JetReq[bin-1]=mean_JetReq[bin-1]/1000;
    meanSBlow_JetReq[bin-1]=meanSBlow_JetReq[bin-1]/1000;
    meanSBhigh_JetReq[bin-1]=meanSBhigh_JetReq[bin-1]/1000;
  } 

  float StdDev_2JetReq[HistBins]/*={0.}*/,StdDevSBLow_2JetReq[HistBins]/*={0.}*/,StdDevSBHigh_2JetReq[HistBins]/*={0.}*/;
  float mean_2JetReq[HistBins]/*={0.}*/,meanSBlow_2JetReq[HistBins]/*={0.}*/,meanSBhigh_2JetReq[HistBins]/*={0.}*/;
    
  for(int i=0;i<HistBins;i++){
    StdDev_2JetReq[i]=0.;StdDevSBLow_2JetReq[i]=0.;StdDevSBHigh_2JetReq[i]=0.;mean_2JetReq[i]=0.;meanSBlow_2JetReq[i]=0.;meanSBhigh_2JetReq[i]=0.;
  }

  TString title_2JetReq="Toys/ee_2JetReqMet_toy_1";
  TString titleSBlow_2JetReq="Toys/eeSideband_2JetReqLowMet_toy_1";TString titleSBhigh_2JetReq="Toys/eeSideband_2JetReqHighMet_toy_1";
  TH1F* hist_2JetReq=(TH1F*)fin.Get(title_2JetReq);hist_2JetReq->Sumw2();
  TH1F* histSBlow_2JetReq=(TH1F*)fin.Get(titleSBlow_2JetReq);histSBlow_2JetReq->Sumw2();
  TH1F* histSBhigh_2JetReq=(TH1F*)fin.Get(titleSBhigh_2JetReq);histSBhigh_2JetReq->Sumw2();
    
  for(int bin=1;bin<HistBins+1;bin++){
    for(int i=1;i<1001;i++){
      title_2JetReq="Toys/ee_2JetReqMet_toy_";title_2JetReq+=i;
      titleSBlow_2JetReq="Toys/eeSideband_2JetReqLowMet_toy_";titleSBlow_2JetReq+=i;
      titleSBhigh_2JetReq="Toys/eeSideband_2JetReqHighMet_toy_";titleSBhigh_2JetReq+=i;
      hist_2JetReq=(TH1F*)fin.Get(title_2JetReq);if(hist_2JetReq->GetSumw2()==0)hist_2JetReq->Sumw2();
      histSBlow_2JetReq=(TH1F*)fin.Get(titleSBlow_2JetReq);if(histSBlow_2JetReq->GetSumw2()==0)histSBlow_2JetReq->Sumw2();
      histSBhigh_2JetReq=(TH1F*)fin.Get(titleSBhigh_2JetReq);if(histSBhigh_2JetReq->GetSumw2()==0)histSBhigh_2JetReq->Sumw2();
      mean_2JetReq[bin-1]+=hist_2JetReq->GetBinContent(bin);
      //cout<<hist->GetBinContent(bin)<<endl;
      meanSBlow_2JetReq[bin-1]+=histSBlow_2JetReq->GetBinContent(bin);
      meanSBhigh_2JetReq[bin-1]+=histSBhigh_2JetReq->GetBinContent(bin);
    }
    //cout<<"Bin:"<<bin<<"  MeanSum:"<<mean[bin-1]<<"  Mean:"<<mean[bin-1]/1000<<endl;
    mean_2JetReq[bin-1]=mean_2JetReq[bin-1]/1000;
    meanSBlow_2JetReq[bin-1]=meanSBlow_2JetReq[bin-1]/1000;
    meanSBhigh_2JetReq[bin-1]=meanSBhigh_2JetReq[bin-1]/1000;
  }
    

  for(int i=1;i<1001;i++){
    title="Toys/eeMet_toy_";title+=i;
    titleSBlow="Toys/eeSidebandLowMet_toy_";titleSBlow+=i;
    titleSBhigh="Toys/eeSidebandHighMet_toy_";titleSBhigh+=i;
    hist = (TH1F*)fin.Get(title);
    histSBlow = (TH1F*)fin.Get(titleSBlow);
    histSBhigh = (TH1F*)fin.Get(titleSBhigh);
    //
    title_JetReq="Toys/ee_JetReqMet_toy_";title_JetReq+=i;
    titleSBlow_JetReq="Toys/eeSideband_JetReqLowMet_toy_";titleSBlow_JetReq+=i;
    titleSBhigh_JetReq="Toys/eeSideband_JetReqHighMet_toy_";titleSBhigh_JetReq+=i;
    //cout<<"title: "<<title<<endl;
    //fin.Get(title)->Draw("SAME");
    hist_JetReq = (TH1F*)fin.Get(title_JetReq);
    histSBlow_JetReq = (TH1F*)fin.Get(titleSBlow_JetReq);
    histSBhigh_JetReq = (TH1F*)fin.Get(titleSBhigh_JetReq);
    for(int bin=1;bin<HistBins+1;bin++){
      StdDev[bin-1]+= (hist->GetBinContent(bin)-mean[bin-1])*(hist->GetBinContent(bin)-mean[bin-1]);
      StdDevSBLow[bin-1]+=(histSBlow->GetBinContent(bin)-meanSBlow[bin-1])*(histSBlow->GetBinContent(bin)-meanSBlow[bin-1]);
      StdDevSBHigh[bin-1]+=(histSBhigh->GetBinContent(bin)-meanSBhigh[bin-1])*(histSBhigh->GetBinContent(bin)-meanSBhigh[bin-1]);
      //if(hist->GetBinContent(bin)<binMin[bin-1])binMin[bin-1]=hist->GetBinContent(bin);
      //if(hist->GetBinContent(bin)>binMax[bin-1])binMax[bin-1]=hist->GetBinContent(bin);
      StdDev_JetReq[bin-1]+= (hist_JetReq->GetBinContent(bin)-mean_JetReq[bin-1])*(hist_JetReq->GetBinContent(bin)-mean_JetReq[bin-1]);
      StdDevSBLow_JetReq[bin-1]+=(histSBlow_JetReq->GetBinContent(bin)-meanSBlow_JetReq[bin-1])*(histSBlow_JetReq->GetBinContent(bin)-meanSBlow_JetReq[bin-1]);
      StdDevSBHigh_JetReq[bin-1]+=(histSBhigh_JetReq->GetBinContent(bin)-meanSBhigh_JetReq[bin-1])*(histSBhigh_JetReq->GetBinContent(bin)-meanSBhigh_JetReq[bin-1]);
      //if(bin==1)cout<<"StdDev: "<<StdDev[bin-1]<<endl;
    }
  }
  for(int bin=1;bin<HistBins+1;bin++){
    StdDev[bin-1]=sqrt(StdDev[bin-1]/1000);
    StdDevSBLow[bin-1]=sqrt(StdDevSBLow[bin-1]/1000);
    StdDevSBHigh[bin-1]=sqrt(StdDevSBHigh[bin-1]/1000);
    StdDev_JetReq[bin-1]=sqrt(StdDev_JetReq[bin-1]/1000);
    StdDevSBLow_JetReq[bin-1]=sqrt(StdDevSBLow_JetReq[bin-1]/1000);
    StdDevSBHigh_JetReq[bin-1]=sqrt(StdDevSBHigh_JetReq[bin-1]/1000);
    //cout<<"Bin:" <<bin<<endl<<" BinMin:"<<binMin[bin-1]<<endl<< " BinMax:"<<binMax[bin-1]<<endl<<" STDDEV:"<<StdDev[bin-1]<<endl;
  }

  title="Toys/eeMet_toy_1";
  hist = (TH1F*)fin.Get(title);
  hist->SetTitle("eeMet Toys");
  hist->Draw("histo");
  for(int i=2;i<1001;i++){
    title="Toys/eeMet_toy_";title+=i;
    //cout<<"title: "<<title<<endl;
    //fin.Get(title)->Draw("SAME");
    hist = (TH1F*)fin.Get(title);
    hist->Draw("histo SAME");
  }
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MetToys_ee.png");

  titleSBlow="Toys/eeSidebandLowMet_toy_1";
  histSBlow = (TH1F*)fin.Get(titleSBlow);
  histSBlow->SetTitle("eeSideBandLowMet Toys");
  histSBlow->Draw("histo");
  for(int i=2;i<1001;i++){
    titleSBlow="Toys/eeSidebandLowMet_toy_";titleSBlow+=i;
    histSBlow = (TH1F*)fin.Get(titleSBlow);
    histSBlow->Draw("histo SAME");
  }
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MetToys_eeSideBandLow.png");

  titleSBhigh="Toys/eeSidebandHighMet_toy_1";
  histSBhigh = (TH1F*)fin.Get(titleSBhigh);
  histSBhigh->SetTitle("eeSideBandHighMet Toys");
  histSBhigh->Draw("histo");
  for(int i=2;i<1001;i++){
    titleSBhigh="Toys/eeSidebandHighMet_toy_";titleSBhigh+=i;
    histSBhigh = (TH1F*)fin.Get(titleSBhigh);
    histSBhigh->Draw("histo SAME");
  }
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MetToys_eeSideBandHigh.png");

  title="Toys/ee_JetReqMet_toy_1";
  hist = (TH1F*)fin.Get(title);
  hist->SetTitle("ee_JetReqMet Toys");
  hist->Draw("histo");
  for(int i=2;i<1001;i++){
    title="Toys/ee_JetReqMet_toy_";title+=i;
    //cout<<"title: "<<title<<endl;
    //fin.Get(title)->Draw("SAME");
    hist = (TH1F*)fin.Get(title);
    hist->Draw("histo SAME");
  }
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MetToys_ee_JetReq.png");

  titleSBlow_JetReq="Toys/eeSideband_JetReqLowMet_toy_1";
  histSBlow_JetReq = (TH1F*)fin.Get(titleSBlow_JetReq);
  histSBlow_JetReq->SetTitle("eeSideBand_JetReqLowMet Toys");
  histSBlow_JetReq->Draw("histo");
  for(int i=2;i<1001;i++){
    titleSBlow_JetReq="Toys/eeSideband_JetReqLowMet_toy_";titleSBlow_JetReq+=i;
    histSBlow_JetReq = (TH1F*)fin.Get(titleSBlow_JetReq);
    histSBlow_JetReq->Draw("histo SAME");
  }
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MetToys_eeSideBandLow_JetReq.png");

  titleSBhigh="Toys/eeSideband_JetReqHighMet_toy_1";
  histSBhigh = (TH1F*)fin.Get(titleSBhigh);
  histSBhigh->SetTitle("eeSideBand_JetReqHighMet Toys");
  histSBhigh->Draw("histo");
  for(int i=2;i<1001;i++){
    titleSBhigh="Toys/eeSideband_JetReqHighMet_toy_";titleSBhigh+=i;
    histSBhigh = (TH1F*)fin.Get(titleSBhigh);
    histSBhigh->Draw("histo SAME");
  }
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MetToys_eeSideBandHigh_JetReq.png");

  fin.cd();
  //FakeRate=0.0145;
  //Now do 
  gStyle->SetErrorX();
  fLimits_nojet.cd();ggMetData->Write("met_gg_nojet");fin.cd();
  TH1F* ggMetNew = (TH1F*)ggMetData->Rebin(NmetBins,"ggMetNew",xbins);
  //include overflow in last bin
  Double_t ggOverFlow=0.,ggOverFlowErr=0.;int overflowBingg=ggMetData->FindBin(ggMetNew->GetBinLowEdge(NmetBins+1)+.01);
  ggOverFlow=ggMetData->IntegralAndError(overflowBingg,999,ggOverFlowErr);
  float ggLastBinNew=ggMetNew->GetBinContent(NmetBins)+ggOverFlow;
  float ggLastBinNewErr=sqrt(ggMetNew->GetBinError(NmetBins)*ggMetNew->GetBinError(NmetBins)+ggOverFlowErr*ggOverFlowErr);
  ggMetNew->SetBinContent(NmetBins,ggLastBinNew);
  ggMetNew->SetBinError(NmetBins,ggLastBinNewErr);
  //float egScale = ggMetNew->Integral(0,4)/egMetNew->Integral(0,4);
  //egMetNew->Scale(egScale,"");
  TH1F* egMetForLimit = (TH1F*)egMet->Clone();
  //egMetForLimit->Scale(FakeRate,"");
  TH1F* eeMet_reweightJet_binned = (TH1F*)fin.Get("eeMet_reweightJet_binned");
  TH1F* eeSidebandLowMet_reweightJet_binned = (TH1F*)fin.Get("eeSidebandLowMet_reweightJet_binned");
  TH1F* eeSidebandHighMet_reweightJet_binned = (TH1F*)fin.Get("eeSidebandHighMet_reweightJet_binned");
  /*
    TH1F* eeMetNew=(TH1F*)eeMet_reweightJet_binned->Rebin(NmetBins,"eeMetNew",xbins);
    TH1F* eeSidebandLowMetNew=(TH1F*)eeSidebandLowMet_reweightJet_binned->Rebin(NmetBins,"eeSidebandLowMetNew",xbins);
    TH1F* eeSidebandHighMetNew=(TH1F*)eeSidebandHighMet_reweightJet_binned->Rebin(NmetBins,"eeSidebandHighMetNew",xbins);
  */
  //Try unbinned
  eeMetNew=(TH1F*)eeMet_reweightJet_binned->Clone();//Rebin(NmetBins,"eeMetNew",xbins);
  TH1F* eeSidebandLowMetNew=(TH1F*)eeSidebandLowMet_reweightJet_binned->Clone();//Rebin(NmetBins,"eeSidebandLowMetNew",xbins);
  TH1F* eeSidebandHighMetNew=(TH1F*)eeSidebandHighMet_reweightJet_binned->Clone();//Rebin(NmetBins,"eeSidebandHighMetNew",xbins);


  //TH1F* eeMetMinusSideBand = new TH1F("eeMetMinusSideBand","",NmetBins,xbins);
  TH1F* eeMetMinusSideBand = new TH1F("eeMetMinusSideBand","",200,0.,1000.);//NmetBins,xbins);

  eeMetMinusSideBand->Sumw2();
  eeMetMinusSideBand->Add(eeMetNew,eeSidebandLowMetNew,1,-1);
  eeMetMinusSideBand->Add(eeSidebandHighMetNew,-1);
  eeMetMinusSideBand->GetXaxis()->SetRangeUser(0,400);
  eeMetMinusSideBand->GetYaxis()->SetRangeUser(1e-3,1.2e5);
  eeMetMinusSideBand->Draw();
  eeMetFromZZ->SetLineColor(kRed);
  eeMetFromZZ->Draw("SAMES");
  eeMetFromWZ->SetLineColor(kBlue);
  eeMetFromWZ->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/eeMetAndZZandWZ.png");
  eeMetMinusSideBand->Add(eeMetFromZZ,-1);
  eeMetMinusSideBand->Add(eeMetFromWZ,-1);

  eeMetMinusSideBandForErrors=(TH1F*)eeMetMinusSideBand->Rebin(NmetBins,"eeMetMinusSideBandForErrors",xbins);
  eeMetMinusSideBand->Draw();
  ggMetData->SetLineColor(kRed);
  ggMetData->Draw("SAME");
  egMetForLimit->SetLineColor(kBlue);
  egMetForLimit->Draw("SAME");
  c1->Print("test2.png");
  TH1F* ffMetForee = (TH1F*)fin.Get("ffMet_reweightJet_binned");
  //TH1F* ffMetForee = (TH1F*)fin.Get("gfMet_reweightJet_binned");
  TH1F* egMet_ForLimits = (TH1F*)egMet->Clone();
  TH1F* eeMetNew_ForLimits = (TH1F*)eeMetNew->Clone();
  TH1F* eeMetNew_ForLimitsStatOnly = (TH1F*)eeMetNew->Clone();
  TH1F* egMetForff=(TH1F*)egMet->Clone();
  TH1F* eeMetMinusSideBand2=(TH1F*)eeMetNew->Clone();

  GetAndSetEEerrors(ggMetData,egMet,egMet_ForLimits,eeMetNew,eeMetMinusSideBand2,eeMetNew_ForLimits,eeMetNew_ForLimitsStatOnly,eeSidebandHighMetNew,eeSidebandLowMetNew,ffMetForee,FakeRate,FakeRateErr,FakeRateErrForLimits,StdDev,StdDevSBHigh,StdDevSBLow,normErree,reweightErree,statErree,diffFromffErroree,eeScale,eeScaleErr,statErreg,normErreg,eeMetFromZZ,eeMetFromWZ);

  fLimits_nojet.cd();
  egMet_ForLimits->Write("met_eg_nojet");
  if(doffgfcomb==false)eeMetNew_ForLimits->Write("met_ee_nojet");
  fin.cd();


  TH1F* egMetNew = (TH1F*)egMet->Rebin(NmetBins,"egMetNew",xbins);
  //include overflow in last bin
  Double_t egOverFlow=0.,egOverFlowErr=0.;int overflowBineg=egMet->FindBin(egMetNew->GetBinLowEdge(NmetBins+1)+.01);
  egOverFlow=egMet->IntegralAndError(overflowBineg,999,egOverFlowErr);
  float egLastBinNew=egMetNew->GetBinContent(NmetBins)+egOverFlow;
  float egLastBinNewErr=sqrt(egMetNew->GetBinError(NmetBins)*egMetNew->GetBinError(NmetBins)+egOverFlowErr*egOverFlowErr);
  egMetNew->SetBinContent(NmetBins,egLastBinNew);
  egMetNew->SetBinError(NmetBins,egLastBinNewErr);
  //  TH1F* eeMetMinusSideBand2=(TH1F*)eeMetNew->Rebin(NmetBins,"eeMetMinusSideBand2",xbins);
    
  //eeMetNew->Scale( ggMetNew->Integral(0,4)/eeMetNew->Integral(0,4) , "");
  TH1F *totalBG=new TH1F("totalBG","",NmetBins,xbins);
  totalBG->Sumw2();
  totalBG->Add(egMetNew,eeMetMinusSideBand2,1,1);

  gg0_20=ggMetNew->IntegralAndError(0,4,gg0_20Error);
  eg0_20=egMetNew->IntegralAndError(0,4,eg0_20Error);
  ee0_20=eeMetMinusSideBand2->IntegralAndError(0,4,ee0_20Error);
  QCDee0_20=totalBG->IntegralAndError(0,4,QCDee0_20Error);
  gg30_50=ggMetNew->IntegralAndError(ggMetNew->FindBin(30),ggMetNew->FindBin(50-1),gg30_50Error);
  eg30_50 = egMetNew->IntegralAndError(egMetNew->FindBin(30),egMetNew->FindBin(50-1),eg30_50Error);
  ee30_50=eeMetMinusSideBand2->IntegralAndError(eeMetMinusSideBand2->FindBin(30),eeMetMinusSideBand2->FindBin(50-1),ee30_50Error);
  QCDee30_50=totalBG->IntegralAndError(totalBG->FindBin(30),totalBG->FindBin(50-1),QCDee30_50Error);
  //cout<<"Bin(30):"<<totalBG->FindBin(30)<<" Bin(50-1):"<<totalBG->FindBin(50-1)<<" Bin(50):"<<totalBG->FindBin(50)<<" Bin(100):"<<totalBG->FindBin(100)<<endl;
  gg50up=ggMetNew->IntegralAndError(ggMetNew->FindBin(50),-1,gg50upError);
  eg50up=egMetNew->IntegralAndError(egMetNew->FindBin(50),-1,eg50upError);
  ee50up=eeMetMinusSideBand2->IntegralAndError(eeMetMinusSideBand2->FindBin(50),-1,ee50upError);
  QCDee50up=totalBG->IntegralAndError(totalBG->FindBin(50),-1,QCDee50upError);
  gg100up=ggMetNew->IntegralAndError(ggMetNew->FindBin(100),-1,gg100upError);
  eg100up=egMetNew->IntegralAndError(egMetNew->FindBin(100),-1,eg100upError);
  ee100up=eeMetMinusSideBand2->IntegralAndError(eeMetMinusSideBand2->FindBin(100),-1,ee100upError);
  QCDee100up=totalBG->IntegralAndError(totalBG->FindBin(100),-1,QCDee100upError);
  ee50_60=eeMetMinusSideBand2->IntegralAndError(eeMetMinusSideBand2->FindBin(50),eeMetMinusSideBand2->FindBin(60-1),ee50_60Error);
  QCDee50_60=totalBG->IntegralAndError(totalBG->FindBin(50),totalBG->FindBin(60-1),QCDee50_60Error);
  ee60_70=eeMetMinusSideBand2->IntegralAndError(eeMetMinusSideBand2->FindBin(60),eeMetMinusSideBand2->FindBin(70-1),ee60_70Error);
  QCDee60_70=totalBG->IntegralAndError(totalBG->FindBin(60),totalBG->FindBin(70-1),QCDee60_70Error);
  ee70_80=eeMetMinusSideBand2->IntegralAndError(eeMetMinusSideBand2->FindBin(70),eeMetMinusSideBand2->FindBin(80-1),ee70_80Error);
  QCDee70_80=totalBG->IntegralAndError(totalBG->FindBin(70),totalBG->FindBin(80-1),QCDee70_80Error);
  ee80_100=eeMetMinusSideBand2->IntegralAndError(eeMetMinusSideBand2->FindBin(80),eeMetMinusSideBand2->FindBin(100-1),ee80_100Error);
  QCDee80_100=totalBG->IntegralAndError(totalBG->FindBin(80),totalBG->FindBin(100-1),QCDee80_100Error);
  gg50_60=ggMetNew->IntegralAndError(ggMetNew->FindBin(50),ggMetNew->FindBin(60-1),gg50_60Error);
  gg60_70=ggMetNew->IntegralAndError(ggMetNew->FindBin(60),ggMetNew->FindBin(70-1),gg60_70Error);
  gg70_80=ggMetNew->IntegralAndError(ggMetNew->FindBin(70),ggMetNew->FindBin(80-1),gg70_80Error);
  gg80_100=ggMetNew->IntegralAndError(ggMetNew->FindBin(80),ggMetNew->FindBin(100-1),gg80_100Error);
  eg50_60=egMetNew->IntegralAndError(egMetNew->FindBin(50),egMetNew->FindBin(60-1),eg50_60Error);
  eg60_70=egMetNew->IntegralAndError(egMetNew->FindBin(60),egMetNew->FindBin(70-1),eg60_70Error);
  eg70_80=egMetNew->IntegralAndError(egMetNew->FindBin(70),egMetNew->FindBin(80-1),eg70_80Error);
  eg80_100=egMetNew->IntegralAndError(egMetNew->FindBin(80),egMetNew->FindBin(100-1),eg80_100Error);
   
  TH1F *eeOverBinWidth = new TH1F("eeOverBinWidth","",NmetBins,xbins);eeOverBinWidth->Sumw2();
  TH1F *egOverBinWidth = new TH1F("egOverBinWidth","",NmetBins,xbins);egOverBinWidth->Sumw2();
  TH1F *ggOverBinWidth = new TH1F("ggOverBinWidth","",NmetBins,xbins);ggOverBinWidth->Sumw2();
  TH1F *ggSigOverBinWidth = new TH1F("ggSigOverBinWidth","",NmetBins,xbins);ggSigOverBinWidth->Sumw2();
  for(int i=1;i<NmetBins+1;i++){
    float eg = egMetNew->GetBinContent(i)/egMetNew->GetBinWidth(i);
    float ee = eeMetMinusSideBand2->GetBinContent(i)/eeMetMinusSideBand2->GetBinWidth(i);
    float gg = ggMetNew->GetBinContent(i)/ggMetNew->GetBinWidth(i);
    float ggE = ggMetNew->GetBinError(i)/ggMetNew->GetBinWidth(i);
    float bg = totalBG->GetBinContent(i)/totalBG->GetBinWidth(i);
    float bgE = totalBG->GetBinError(i)/totalBG->GetBinWidth(i);
    float ggS = ggSigNew->GetBinContent(i)/ggSigNew->GetBinWidth(i);
    float ggSE = ggSigNew->GetBinError(i)/ggSigNew->GetBinWidth(i);
    float sig=ggSig2012Rebin->GetBinContent(i)/ggSig2012Rebin->GetBinWidth(i);
    ggSig2012Rebin->SetBinContent(i,sig);
    float sig2=ggSig2012Rebin_2->GetBinContent(i)/ggSig2012Rebin_2->GetBinWidth(i);
    ggSig2012Rebin_2->SetBinContent(i,sig2);
    egOverBinWidth->SetBinContent(i,eg);
    eeOverBinWidth->SetBinContent(i,ee);
    ggOverBinWidth->SetBinContent(i,gg);
    ggOverBinWidth->SetBinError(i,ggE);
    totalBG->SetBinContent(i,bg);
    totalBG->SetBinError(i,bgE);
    ggSigOverBinWidth->SetBinContent(i,ggS);
    ggSigOverBinWidth->SetBinError(i,ggSE);
  }
  ggSigOverBinWidth->SetLineColor(kPink);
  ggOverBinWidth->SetMarkerColor(kBlack);
  ggOverBinWidth->SetLineColor(kBlack);
  ggOverBinWidth->SetLineWidth(2.2);
  ggOverBinWidth->SetMarkerSize(.6);
  ggOverBinWidth->SetStats(0);
  float Max = ggOverBinWidth->GetMaximum()+70000;
  float Min = .008;
  egOverBinWidth->SetMaximum(Max);
  egOverBinWidth->GetYaxis()->SetTitle("Number of Events / GeV");
  egOverBinWidth->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  egOverBinWidth->SetTitle("");
  egOverBinWidth->SetLineColor(kAzure);
  egOverBinWidth->SetFillColor(kAzure);
  //egOverBinWidth->SetFillStyle(3005);
  egOverBinWidth->SetStats(0);
  eeOverBinWidth->SetLineColor(kGray+1);
  eeOverBinWidth->SetFillColor(kGray+1);
  //eeOverBinWidth->SetFillStyle(3004);
  eeOverBinWidth->SetMarkerSize(0);
  eeOverBinWidth->SetStats(0);
  eeOverBinWidth->SetTitle("");
  totalBG->SetMarkerSize(0);
  totalBG->SetFillColor(kRed);
  //gStyle->SetHatchesSpacing(5);
  totalBG->SetFillStyle(3154);
  totalBG->SetStats(0);
  THStack *eeStack = new THStack("eeStack",";;Number of Events / GeV");
  eeStack->SetMaximum(Max);
  eeStack->SetMinimum(Min);
  eeStack->Add(egOverBinWidth);
  eeStack->Add(eeOverBinWidth); 
  eeStack->Draw();
  totalBG->Draw("E2SAME");
  ggOverBinWidth->Draw("PESAME");
  //ggSigOverBinWidth->Draw("histo SAME");
  // ggSig720New->Draw("SAME");
  //ggSig800New->Draw("SAME");
  TLegend *legEE = new TLegend(.5,.45,.8,.65,"","brNDC");
  legEE->AddEntry(ggOverBinWidth,"#gamma#gamma Candidate Sample","lpf");
  legEE->AddEntry(totalBG,"QCD + Electroweak Error","f");
  legEE->AddEntry(eeOverBinWidth,"QCD (e^{+}e^{-} sample)","f");
  legEE->AddEntry(egOverBinWidth,"Electroweak","f");
  legEE->SetFillColor(kWhite);
  legEE->Draw("SAME");
  TPaveText *Text;
  Text = new TPaveText(.5,.66,.8,.87,"NDC");
  Text->AddText("CMS Preliminary");
  Text->AddText("");
  Text->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  Text->AddText("No Jet Requirement");
  Text->SetFillStyle(4000);
  Text->SetFillColor(0);
  Text->SetBorderSize(0);
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_ggANDeeBinnedReweighted.png");

  TLine *line1 = new TLine(0,1,metPlotXmax,1);//line1->SetLineColor(kBlue);

  TH1F* ggOveree=(TH1F*)ggOverBinWidth->Clone();ggOveree->SetTitle("ggOveree");ggOveree->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  ggOveree->Divide(eeOverBinWidth);
  c1->SetLogy(0);
  ggOveree->Draw();
  line1->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggOveree.png");
  fout.cd();ggOveree->Write("ggOveree");fin.cd();
  totalBGee = (TH1F*)totalBG->Clone();
  TH1F* ggOvereeQCD=(TH1F*)ggOverBinWidth->Clone();ggOvereeQCD->SetTitle("ggOvereeQCD");ggOvereeQCD->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  TH1F* ggOvereeQCDerr=(TH1F*)totalBG->Clone();
  for(int i=1;i<ggOvereeQCD->GetNbinsX()+1;i++){
    float Value = ggOvereeQCD->GetBinContent(i);float StatErr = ggOvereeQCD->GetBinError(i);
    Value/=totalBG->GetBinContent(i);StatErr/=totalBG->GetBinContent(i);
    ggOvereeQCD->SetBinContent(i,Value);ggOvereeQCD->SetBinError(i,StatErr);
    float SystErr = totalBG->GetBinError(i)/totalBG->GetBinContent(i);
    ggOvereeQCDerr->SetBinContent(i,1.);ggOvereeQCDerr->SetBinError(i,SystErr);
  }
  c1->SetLogy(0);   
  ggOvereeQCDerr->SetFillColor(kRed);
  ggOvereeQCDerr->GetYaxis()->SetRangeUser(0.0,2.0);
  ggOvereeQCDerr->Draw("E2");
  ggOvereeQCD->Draw("SAME");
  line1->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggOvereeQCD.png");
  fout.cd();ggOvereeQCD->Write("ggOvereeQCD");ggOvereeQCDerr->Write("ggOvereeQCDerr");fin.cd();
  c1->SetLogy(1);

  c2->cd();
  p1->cd();
  p1->SetLogy(1);
  /*  egOverBinWidth->Draw("B");
      eeOverBinWidth->Draw("BSAME");
      //eeOverBinWidth->Draw("SAME");
      egOverBinWidth->Draw("BSAME");
      //egOverBinWidth->Draw("SAME");*/
  eeStack->Draw();
  totalBG->Draw("E2SAME");
  ggOverBinWidth->Draw("PESAME");
  ggSig2012Rebin->Draw("histo SAME");
  ggSig2012Rebin_2->Draw("histo SAME");
  TLegend *legEEwr = new TLegend(.5,.35,.83,.61,"","brNDC");
  legEEwr->AddEntry(ggOverBinWidth,"#gamma#gamma Candidate Sample","lpf");
  legEEwr->AddEntry(totalBG,"QCD + Electroweak Error","f");
  legEEwr->AddEntry(eeOverBinWidth,"QCD (e^{+}e^{-} sample)","f");
  legEEwr->AddEntry(egOverBinWidth,"Electroweak","f");
  legEEwr->AddEntry(ggSig2012Rebin_2,"GGM #gamma#gamma (1100_720_375)","l");
  legEEwr->AddEntry(ggSig2012Rebin,"GGM #gamma#gamma (1400_1720_375)","l");
  legEEwr->SetFillColor(kWhite);
  legEEwr->Draw("SAME");
  TPaveText *Textwr;
  Textwr = new TPaveText(.5,.62,.8,.83,"NDC");
  Textwr->AddText("CMS Preliminary");
  Textwr->AddText("");
  Textwr->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  Textwr->AddText("No Jet Requirement");
  Textwr->SetFillColor(0);
  Textwr->SetBorderSize(0);
  Textwr->Draw();
  p2->cd();
  ggOvereeQCDerr->SetTitle("");
  ggOvereeQCDerr->GetYaxis()->SetTitle("Data/Prediction");
  //ggOvereeQCDerr->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  //ggOvereeQCDerr->GetXaxis()->SetTitleOffset(0.);
  ggOvereeQCDerr->GetYaxis()->SetTitleOffset(0.48);
  ggOvereeQCDerr->GetYaxis()->SetTitleSize(0.15);
  ggOvereeQCDerr->GetXaxis()->SetTitleSize(0.2);
  ggOvereeQCDerr->GetYaxis()->SetLabelSize(0.12);
  ggOvereeQCDerr->GetXaxis()->SetLabelSize(0.15);
  float ratMin=0.6,ratMax=1.3;
  ggOvereeQCDerr->GetYaxis()->SetRangeUser(ratMin,ratMax);
  ggOvereeQCDerr->Draw("E2");
  line1->Draw("SAME");
  ggOvereeQCD->Draw("SAME");
  p3->cd();
  TPaveText *met = new TPaveText(.65,.5,.85,.99,"NDC");
  met->AddText("E_{T}^{miss} [GeV]");
  met->SetFillColor(0);
  met->SetBorderSize(0);
  met->SetTextSize(.3);
  met->Draw();
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDeeBinnedReweighted_WithRatio.png");
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDeeBinnedReweighted_WithRatio.pdf");
 
  fout.cd();eeOverBinWidth->Write("eeMetBinnedReweightedAndScaled");ggOverBinWidth->Write("ggMetReweightedAndScaled");fin.cd();
  c1->cd();
  //---------------Now ee Jet_Req----------------- //fix sideband
  TH1F* ggMet_JetReqForLimit=(TH1F*)fin.Get("ggMet_JetReq");
  //TH1F* egMet_JetReq=(TH1F*)fin.Get("egMet_JetReq");
  //float egScale_JetReq = ggMetNew_JetReq->Integral(0,4)/egMetNew_JetReq->Integral(0,4);
  //egMetNew_JetReq->Scale(egScale_JetReq,"");
  //egMet_JetReqForLimit->Scale(FakeRate,"");
  fLimits_1jet.cd();
  ggMet_JetReqForLimit->Write("met_gg_1jet");
  fin.cd();

  TH1F* ggMetNew_JetReq=(TH1F*)ggMet_JetReqForLimit->Rebin(NmetBins,"ggMetNew_JetReq",xbins);
  //include overflow in last bin
  ggOverFlow=0.,ggOverFlowErr=0.;overflowBingg=ggMet_JetReqForLimit->FindBin(ggMetNew_JetReq->GetBinLowEdge(NmetBins+1)+.01);
  ggOverFlow=ggMet_JetReqForLimit->IntegralAndError(overflowBingg,999,ggOverFlowErr);
  ggLastBinNew=ggMetNew_JetReq->GetBinContent(NmetBins)+ggOverFlow;
  ggLastBinNewErr=sqrt(ggMetNew_JetReq->GetBinError(NmetBins)*ggMetNew_JetReq->GetBinError(NmetBins)+ggOverFlowErr*ggOverFlowErr);
  ggMetNew_JetReq->SetBinContent(NmetBins,ggLastBinNew);
  ggMetNew_JetReq->SetBinError(NmetBins,ggLastBinNewErr);

  TH1F* eeMet_reweightJet_binned_JetReq=(TH1F*)fin.Get("eeMet_reweightJet_binned_JetReq");
  TH1F* eeSidebandLowMet_reweightJet_binned_JetReq = (TH1F*)fin.Get("eeSidebandLowMet_reweightJet_binned_JetReq");
  TH1F* eeSidebandHighMet_reweightJet_binned_JetReq = (TH1F*)fin.Get("eeSidebandHighMet_reweightJet_binned_JetReq");
  /* TH1F* eeMetNew_JetReq=(TH1F*)eeMet_reweightJet_binned_JetReq->Clone();Rebin(NmetBins,"eeMetNew_JetReq",xbins);
     TH1F* eeSidebandLowMetNew_JetReq=(TH1F*)eeSidebandLowMet_reweightJet_binned_JetReq->Rebin(NmetBins,"eeSidebandLowMetNew_JetReq",xbins);
     TH1F* eeSidebandHighMetNew_JetReq=(TH1F*)eeSidebandHighMet_reweightJet_binned_JetReq->Rebin(NmetBins,"eeSidebandHighMetNew_JetReq",xbins);*/
  //try unbinned

  TH1F* eeMetNew_JetReq=(TH1F*)eeMet_reweightJet_binned_JetReq->Clone();//;Rebin(NmetBins,"eeMetNew_JetReq",xbins);
  TH1F* eeSidebandLowMetNew_JetReq=(TH1F*)eeSidebandLowMet_reweightJet_binned_JetReq->Clone();//Rebin(NmetBins,"eeSidebandLowMetNew_JetReq",xbins);
  TH1F* eeSidebandHighMetNew_JetReq=(TH1F*)eeSidebandHighMet_reweightJet_binned_JetReq->Clone();//Rebin(NmetBins,"eeSidebandHighMetNew_JetReq",xbins);


  //TH1F* eeSidebandMetNew_JetReq=(TH1F*)eeSidebandMet_reweightJet_binned_JetReq->Rebin(NmetBins,"eeSidebandMetNew_JetReq",xbins);
  TH1F* eeMetMinusSideBand_JetReq = new TH1F("eeMetMinusSideBand_JetReq","",200,0.,1000.);//NmetBins,xbins);

  eeMetMinusSideBand_JetReq->Sumw2();
  eeMetMinusSideBand_JetReq->Add(eeMetNew_JetReq,eeSidebandLowMetNew_JetReq,1,-1);
  eeMetMinusSideBand_JetReq->Add(eeSidebandHighMetNew_JetReq,-1);
  eeMetMinusSideBand_JetReq->GetXaxis()->SetRangeUser(0,400);
  eeMetMinusSideBand_JetReq->GetYaxis()->SetRangeUser(1e-3,1.2e5);
  eeMetMinusSideBand_JetReq->Draw();
  eeMetFromZZ_JetReq->SetLineColor(kRed);
  eeMetFromZZ_JetReq->Draw("SAMES");
  eeMetFromWZ_JetReq->SetLineColor(kBlue);
  eeMetFromWZ_JetReq->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/eeMetAndZZandWZ_JetReq.png");
  eeMetMinusSideBand_JetReq->Add(eeMetFromZZ_JetReq,-1);
  eeMetMinusSideBand_JetReq->Add(eeMetFromWZ_JetReq,-1);
  TH1F* eeMetNew_ForLimits_JetReq=(TH1F*)eeMetNew_JetReq->Clone();
  TH1F* eeMetNew_ForLimitsStatOnly_JetReq=(TH1F*)eeMetNew_JetReq->Clone();
  TH1F* ffMetForee_JetReq = (TH1F*)fin.Get("ffMet_reweightJet_binned_JetReq");
  //TH1F* ffMetForee_JetReq = (TH1F*)fin.Get("gfMet_reweightJet_binned_JetReq");
  TH1F* egMet_ForLimits_JetReq=(TH1F*)egMet_JetReq->Clone();
  TH1F* egMetForff_JetReq=(TH1F*)egMet_JetReq->Clone();
  TH1F* eeMetMinusSideBand2_JetReq=(TH1F*)eeMetNew_JetReq->Clone();
  GetAndSetEEerrors(ggMet_JetReq,egMet_JetReq,egMet_ForLimits_JetReq,eeMetNew_JetReq,eeMetMinusSideBand2_JetReq,eeMetNew_ForLimits_JetReq,eeMetNew_ForLimitsStatOnly_JetReq,eeSidebandHighMetNew_JetReq,eeSidebandLowMetNew_JetReq,ffMetForee_JetReq,FakeRate,FakeRateErr,FakeRateErrForLimits,StdDev_JetReq,StdDevSBHigh_JetReq,StdDevSBLow_JetReq,normErree_JetReq,reweightErree_JetReq,statErree_JetReq,diffFromffErroree_JetReq,eeScale_JetReq,eeScaleErr_JetReq,statErreg_JetReq,normErreg_JetReq,eeMetFromZZ_JetReq,eeMetFromWZ_JetReq);
 
  eeMetMinusSideBandForErrors_JetReq=(TH1F*)eeMetMinusSideBand_JetReq->Rebin(NmetBins,"eeMetMinusSideBandForErrors_JetReq",xbins);
 
  //cout<<"Integral Over 100, ee JetReq.  from eeMetNew_JetReq:"<<eeMetNew_JetReq->Integral(eeMetNew_JetReq->FindBin(100.1),-1)<<"  from eeMetMinusSideBand2:"<<eeMetMinusSideBand2_JetReq->Integral(eeMetMinusSideBand2_JetReq->FindBin(100.1),-1)<<endl;   

  fLimits_1jet.cd();
  egMet_ForLimits_JetReq->Write("met_eg_1jet");
  if(doffgfcomb==false)eeMetNew_ForLimits_JetReq->Write("met_ee_1jet");
  fin.cd();
 
  TH1F* egMetNew_JetReq=(TH1F*)egMet_JetReq->Rebin(NmetBins,"egMetNew_JetReq",xbins);
  //include overflow in last bin
  egOverFlow=0.,egOverFlowErr=0.;overflowBineg=egMet_JetReq->FindBin(egMetNew_JetReq->GetBinLowEdge(NmetBins+1)+.01);
  egOverFlow=egMet_JetReq->IntegralAndError(overflowBineg,999,egOverFlowErr);
  egLastBinNew=egMetNew_JetReq->GetBinContent(NmetBins)+egOverFlow;
  egLastBinNewErr=sqrt(egMetNew_JetReq->GetBinError(NmetBins)*egMetNew_JetReq->GetBinError(NmetBins)+egOverFlowErr*egOverFlowErr);
  egMetNew_JetReq->SetBinContent(NmetBins,egLastBinNew);
  egMetNew_JetReq->SetBinError(NmetBins,egLastBinNewErr);
  //TH1F* eeMetMinusSideBand2_JetReq=(TH1F*)eeMetNew_JetReq->Rebin(NmetBins,"eeMetMinusSideBand2_JetReq",xbins);

  //eeMetNew_JetReq->Scale( ggMetNew_JetReq->Integral(0,4)/eeMetNew_JetReq->Integral(0,4) , "");
  TH1F *totalBG_JetReq=new TH1F("totalBG_JetReq","",NmetBins,xbins);
  totalBG_JetReq->Sumw2();
  totalBG_JetReq->Add(egMetNew_JetReq,eeMetMinusSideBand2_JetReq,1,1);
   
  gg0_20_JetReq=ggMetNew_JetReq->IntegralAndError(0,4,gg0_20Error_JetReq);
  eg0_20_JetReq=egMetNew_JetReq->IntegralAndError(0,4,eg0_20Error_JetReq);
  ee0_20_JetReq=eeMetMinusSideBand2_JetReq->IntegralAndError(0,4,ee0_20Error_JetReq);
  QCDee0_20_JetReq=totalBG_JetReq->IntegralAndError(0,4,QCDee0_20Error_JetReq);
  gg30_50_JetReq=ggMetNew_JetReq->IntegralAndError(ggMetNew_JetReq->FindBin(30),ggMetNew_JetReq->FindBin(50-1),gg30_50Error_JetReq);
  eg30_50_JetReq=egMetNew_JetReq->IntegralAndError(egMetNew_JetReq->FindBin(30),egMetNew_JetReq->FindBin(50-1),eg30_50Error_JetReq);
  ee30_50_JetReq=eeMetMinusSideBand2_JetReq->IntegralAndError(eeMetMinusSideBand2_JetReq->FindBin(30),eeMetMinusSideBand2_JetReq->FindBin(50-1),ee30_50Error_JetReq);
  QCDee30_50_JetReq=totalBG_JetReq->IntegralAndError(totalBG_JetReq->FindBin(30),totalBG_JetReq->FindBin(50-1),QCDee30_50Error_JetReq);
  gg50up_JetReq=ggMetNew_JetReq->IntegralAndError(ggMetNew_JetReq->FindBin(50),-1,gg50upError_JetReq);
  eg50up_JetReq=egMetNew_JetReq->IntegralAndError(egMetNew_JetReq->FindBin(50),-1,eg50upError_JetReq);
  ee50up_JetReq=eeMetMinusSideBand2_JetReq->IntegralAndError(eeMetMinusSideBand2_JetReq->FindBin(50),-1,ee50upError_JetReq);
  QCDee50up_JetReq=totalBG_JetReq->IntegralAndError(totalBG_JetReq->FindBin(50),-1,QCDee50upError_JetReq);
  gg100up_JetReq=ggMetNew_JetReq->IntegralAndError(ggMetNew_JetReq->FindBin(100),-1,gg100upError_JetReq);
  eg100up_JetReq=egMetNew_JetReq->IntegralAndError(egMetNew_JetReq->FindBin(100),-1,eg100upError_JetReq);
  ee100up_JetReq=eeMetMinusSideBand2_JetReq->IntegralAndError(eeMetMinusSideBand2_JetReq->FindBin(100),-1,ee100upError_JetReq);
  QCDee100up_JetReq=totalBG_JetReq->IntegralAndError(totalBG_JetReq->FindBin(100),-1,QCDee100upError_JetReq);
  ee50_60_JetReq=eeMetMinusSideBand2_JetReq->IntegralAndError(eeMetMinusSideBand2_JetReq->FindBin(50),eeMetMinusSideBand2_JetReq->FindBin(60-1),ee50_60Error_JetReq);
  QCDee50_60_JetReq=totalBG_JetReq->IntegralAndError(totalBG_JetReq->FindBin(50),totalBG_JetReq->FindBin(60-1),QCDee50_60Error_JetReq);
  ee60_70_JetReq=eeMetMinusSideBand2_JetReq->IntegralAndError(eeMetMinusSideBand2_JetReq->FindBin(60),eeMetMinusSideBand2_JetReq->FindBin(70-1),ee60_70Error_JetReq);
  QCDee60_70_JetReq=totalBG_JetReq->IntegralAndError(totalBG_JetReq->FindBin(60),totalBG_JetReq->FindBin(70-1),QCDee60_70Error_JetReq);
  ee70_80_JetReq=eeMetMinusSideBand2_JetReq->IntegralAndError(eeMetMinusSideBand2_JetReq->FindBin(70),eeMetMinusSideBand2_JetReq->FindBin(80-1),ee70_80Error_JetReq);
  QCDee70_80_JetReq=totalBG_JetReq->IntegralAndError(totalBG_JetReq->FindBin(70),totalBG_JetReq->FindBin(80-1),QCDee70_80Error_JetReq);
  ee80_100_JetReq=eeMetMinusSideBand2_JetReq->IntegralAndError(eeMetMinusSideBand2_JetReq->FindBin(80),eeMetMinusSideBand2_JetReq->FindBin(100-1),ee80_100Error_JetReq);
  QCDee80_100_JetReq=totalBG_JetReq->IntegralAndError(totalBG_JetReq->FindBin(80),totalBG_JetReq->FindBin(100-1),QCDee80_100Error_JetReq);
  gg50_60_JetReq=ggMetNew_JetReq->IntegralAndError(ggMetNew_JetReq->FindBin(50),ggMetNew_JetReq->FindBin(60-1),gg50_60Error_JetReq);
  gg60_70_JetReq=ggMetNew_JetReq->IntegralAndError(ggMetNew_JetReq->FindBin(60),ggMetNew_JetReq->FindBin(70-1),gg60_70Error_JetReq);
  gg70_80_JetReq=ggMetNew_JetReq->IntegralAndError(ggMetNew_JetReq->FindBin(70),ggMetNew_JetReq->FindBin(80-1),gg70_80Error_JetReq);
  gg80_100_JetReq=ggMetNew_JetReq->IntegralAndError(ggMetNew_JetReq->FindBin(80),ggMetNew_JetReq->FindBin(100-1),gg80_100Error_JetReq);
  eg50_60_JetReq=egMetNew_JetReq->IntegralAndError(egMetNew_JetReq->FindBin(50),egMetNew_JetReq->FindBin(60-1),eg50_60Error_JetReq);
  eg60_70_JetReq=egMetNew_JetReq->IntegralAndError(egMetNew_JetReq->FindBin(60),egMetNew_JetReq->FindBin(70-1),eg60_70Error_JetReq);
  eg70_80_JetReq=egMetNew_JetReq->IntegralAndError(egMetNew_JetReq->FindBin(70),egMetNew_JetReq->FindBin(80-1),eg70_80Error_JetReq);
  eg80_100_JetReq=egMetNew_JetReq->IntegralAndError(egMetNew_JetReq->FindBin(80),egMetNew_JetReq->FindBin(100-1),eg80_100Error_JetReq);
  
  TH1F *eeOverBinWidth_JetReq = new TH1F("eeOverBinWidth_JetReq","",NmetBins,xbins);eeOverBinWidth_JetReq->Sumw2();
  TH1F *egOverBinWidth_JetReq = new TH1F("egOverBinWidth_JetReq","",NmetBins,xbins);egOverBinWidth_JetReq->Sumw2();
  TH1F *ggOverBinWidth_JetReq = new TH1F("ggOverBinWidth_JetReq","",NmetBins,xbins);ggOverBinWidth_JetReq->Sumw2();
  for(int i=1;i<NmetBins+1;i++){
    float eg = egMetNew_JetReq->GetBinContent(i)/egMetNew_JetReq->GetBinWidth(i);
    float ee = eeMetMinusSideBand2_JetReq->GetBinContent(i)/eeMetMinusSideBand2_JetReq->GetBinWidth(i);
    float gg = ggMetNew_JetReq->GetBinContent(i)/ggMetNew_JetReq->GetBinWidth(i);
    float ggE = ggMetNew_JetReq->GetBinError(i)/ggMetNew_JetReq->GetBinWidth(i);
    float bg = totalBG_JetReq->GetBinContent(i)/totalBG_JetReq->GetBinWidth(i);
    float bgE = totalBG_JetReq->GetBinError(i)/totalBG_JetReq->GetBinWidth(i);
    egOverBinWidth_JetReq->SetBinContent(i,eg);
    eeOverBinWidth_JetReq->SetBinContent(i,ee);
    ggOverBinWidth_JetReq->SetBinContent(i,gg);
    ggOverBinWidth_JetReq->SetBinError(i,ggE);
    totalBG_JetReq->SetBinContent(i,bg);
    totalBG_JetReq->SetBinError(i,bgE);
    float sig_JetReq=ggSig2012Rebin_JetReq->GetBinContent(i)/ggSig2012Rebin_JetReq->GetBinWidth(i);
    ggSig2012Rebin_JetReq->SetBinContent(i,sig_JetReq);
    float sig_JetReq2=ggSig2012Rebin_JetReq_2->GetBinContent(i)/ggSig2012Rebin_JetReq_2->GetBinWidth(i);
    ggSig2012Rebin_JetReq_2->SetBinContent(i,sig_JetReq2);
  }
  ggOverBinWidth_JetReq->SetMarkerColor(kBlack);
  ggOverBinWidth_JetReq->SetLineColor(kBlack);
  ggOverBinWidth_JetReq->SetLineWidth(2.2);
  ggOverBinWidth_JetReq->SetMarkerSize(0.6);
  ggOverBinWidth_JetReq->SetStats(0);    
  egOverBinWidth_JetReq->SetMaximum(Max);
  egOverBinWidth_JetReq->GetYaxis()->SetTitle("Number of Events / GeV");
  egOverBinWidth_JetReq->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  egOverBinWidth_JetReq->SetTitle("");
  egOverBinWidth_JetReq->SetLineColor(kAzure);
  egOverBinWidth_JetReq->SetFillColor(kAzure);
  //egOverBinWidth_JetReq->SetFillStyle(3005);
  egOverBinWidth_JetReq->SetStats(0);
  eeOverBinWidth_JetReq->SetLineColor(kGray+1);
  eeOverBinWidth_JetReq->SetFillColor(kGray+1);
  //eeOverBinWidth_JetReq->SetFillStyle(3004);
  eeOverBinWidth_JetReq->SetMarkerSize(0);
  eeOverBinWidth_JetReq->SetStats(0);
  eeOverBinWidth_JetReq->SetTitle("");
  totalBG_JetReq->SetMarkerSize(0);
  totalBG_JetReq->SetFillColor(kRed);
  totalBG_JetReq->SetFillStyle(3154);
  totalBG_JetReq->SetStats(0);
  THStack *eeStack_JetReq = new THStack("eeStack_JetReq",";;Number of Events / GeV");
  eeStack_JetReq->SetMaximum(Max);
  eeStack_JetReq->SetMinimum(Min);
  eeStack_JetReq->Add(egOverBinWidth_JetReq);
  eeStack_JetReq->Add(eeOverBinWidth_JetReq);
  eeStack_JetReq->Draw();
  totalBG_JetReq->Draw("E2SAME");
  ggOverBinWidth_JetReq->Draw("PESAME");
  TLegend *legEE_JetReq = new TLegend(.5,.45,.8,.65,"","brNDC");
  legEE_JetReq->AddEntry(ggOverBinWidth_JetReq,"#gamma#gamma Candidate Sample","lpf");
  legEE_JetReq->AddEntry(totalBG_JetReq,"QCD + Electroweak Error","f");
  legEE_JetReq->AddEntry(eeOverBinWidth_JetReq,"QCD (e^{+}e^{-} sample)","f");
  legEE_JetReq->AddEntry(egOverBinWidth_JetReq,"Electroweak","f");
  legEE_JetReq->SetFillColor(kWhite);
  legEE_JetReq->Draw("SAME");
  TPaveText *Text_JetReq;
  Text_JetReq= new TPaveText(.5,.66,.8,.87,"NDC");
  Text_JetReq->AddText("CMS Preliminary");
  Text_JetReq->AddText("");
  Text_JetReq->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  Text_JetReq->AddText(">=1 Jet Requirement");
  Text_JetReq->SetFillColor(0);
  Text_JetReq->SetBorderSize(0);
  Text_JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_ggANDeeBinnedReweighted_JetReq.png");
  fout.cd();eeOverBinWidth_JetReq->Write("eeMetBinnedReweightedAndScaled_JetReq");ggOverBinWidth_JetReq->Write("ggMetReweightedAndScaled_JetReq");fin.cd();

  TH1F* ggOveree_JetReq=(TH1F*)ggOverBinWidth_JetReq->Clone();ggOveree_JetReq->SetTitle("ggOveree_JetReq");ggOveree_JetReq->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  ggOveree_JetReq->Divide(eeOverBinWidth_JetReq);
  c1->SetLogy(0);
  ggOveree_JetReq->Draw();
  line1->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggOveree_JetReq.png");
  fout.cd();ggOveree->Write("ggOveree_JetReq.png");fin.cd();
  totalBGee_JetReq = (TH1F*)totalBG_JetReq->Clone();
  TH1F* ggOvereeQCD_JetReq=(TH1F*)ggOverBinWidth_JetReq->Clone();ggOvereeQCD_JetReq->SetTitle("ggOvereeQCD_JetReq");ggOvereeQCD_JetReq->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  TH1F* ggOvereeQCDerr_JetReq=(TH1F*)totalBG_JetReq->Clone();
  for(int i=1;i<ggOvereeQCD_JetReq->GetNbinsX()+1;i++){
    float Value = ggOvereeQCD_JetReq->GetBinContent(i);float StatErr = ggOvereeQCD_JetReq->GetBinError(i);
    Value/=totalBG_JetReq->GetBinContent(i);StatErr/=totalBG_JetReq->GetBinContent(i);
    ggOvereeQCD_JetReq->SetBinContent(i,Value);ggOvereeQCD_JetReq->SetBinError(i,StatErr);
    float SystErr = totalBG_JetReq->GetBinError(i)/totalBG_JetReq->GetBinContent(i);
    ggOvereeQCDerr_JetReq->SetBinContent(i,1.);ggOvereeQCDerr_JetReq->SetBinError(i,SystErr);
  }
  //ggOvereeQCD->Divide(totalBG);
  c1->SetLogy(0);
  ggOvereeQCDerr_JetReq->SetFillColor(kRed);
  ggOvereeQCDerr_JetReq->GetYaxis()->SetRangeUser(0.0,2.0);
  ggOvereeQCDerr_JetReq->Draw("E2");
  ggOvereeQCD_JetReq->Draw("SAME");
  line1->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggOvereeQCD_JetReq.png");
  fout.cd();ggOvereeQCD_JetReq->Write("ggOvereeQCD_JetReq");ggOvereeQCDerr_JetReq->Write("ggOvereeQCDerr_JetReq");fin.cd();
  c1->SetLogy(1);

  c2->cd();
  p1->cd();
  p1->SetLogy(1);
  /* egOverBinWidth_JetReq->Draw("B");
     eeOverBinWidth_JetReq->Draw("BSAME");
     eeOverBinWidth_JetReq->Draw("SAME");
     egOverBinWidth_JetReq->Draw("BSAME");
     egOverBinWidth_JetReq->Draw("SAME");*/
  eeStack_JetReq->Draw();
  totalBG_JetReq->Draw("E2SAME");
  ggOverBinWidth_JetReq->Draw("PESAME");
  ggSig2012Rebin_JetReq->Draw("histo SAME");
  ggSig2012Rebin_JetReq_2->Draw("histo SAME");
  TLegend *legEEwr_JetReq = new TLegend(.5,.35,.83,.61,"","brNDC");
  legEEwr_JetReq->AddEntry(ggOverBinWidth_JetReq,"#gamma#gamma Candidate Sample","lpf");
  legEEwr_JetReq->AddEntry(totalBG_JetReq,"QCD + Electroweak Error","f");
  legEEwr_JetReq->AddEntry(eeOverBinWidth_JetReq,"QCD (e^{+}e^{-} sample)","f");
  legEEwr_JetReq->AddEntry(egOverBinWidth_JetReq,"Electroweak","f");
  legEEwr_JetReq->AddEntry(ggSig2012Rebin_JetReq_2,"GGM #gamma#gamma (1100_720_375)","l");
  legEEwr_JetReq->AddEntry(ggSig2012Rebin_JetReq,"GGM #gamma#gamma (1400_1720_375)","l");
  legEEwr_JetReq->SetFillColor(kWhite);
  legEEwr_JetReq->Draw("SAME");
  TPaveText *Textwr_JetReq;
  Textwr_JetReq = new TPaveText(.5,.62,.8,.83,"NDC");
  Textwr_JetReq->AddText("CMS Preliminary");
  Textwr_JetReq->AddText("");
  Textwr_JetReq->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  Textwr_JetReq->AddText(">=1 Jet Requirement");
  Textwr_JetReq->SetFillColor(0);
  Textwr_JetReq->SetBorderSize(0);
  Textwr_JetReq->Draw();
  p2->cd();
  ggOvereeQCDerr_JetReq->SetTitle("");
  ggOvereeQCDerr_JetReq->GetYaxis()->SetTitle("Data/Prediction");
  ggOvereeQCDerr_JetReq->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  ggOvereeQCDerr_JetReq->GetYaxis()->SetTitleOffset(0.48);
  ggOvereeQCDerr_JetReq->GetYaxis()->SetTitleSize(0.15);
  ggOvereeQCDerr_JetReq->GetXaxis()->SetTitleSize(0.2);
  ggOvereeQCDerr_JetReq->GetYaxis()->SetLabelSize(0.12);
  ggOvereeQCDerr_JetReq->GetXaxis()->SetLabelSize(0.15);
  ggOvereeQCDerr_JetReq->GetYaxis()->SetRangeUser(ratMin,ratMax);
  //ggOvereeQCD->GetYaxis()->SetTitleSize(3);
  ggOvereeQCDerr_JetReq->Draw("E2");
  ggOvereeQCD_JetReq->Draw("SAME");
  line1->Draw("SAME");
  p3->cd();
  met->Draw();
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDeeBinnedReweighted_WithRatio_JetReq.png");
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDeeBinnedReweighted_WithRatio_JetReq.pdf");
  c1->cd();
  //---------Now ee 2JetReq-------
  TH1F* ggMet_2JetReqForLimit=(TH1F*)fin.Get("ggMet_2JetReq");
  //TH1F* egMet_2JetReq=(TH1F*)fin.Get("egMet_2JetReq");
  //float egScale_2JetReq = ggMetNew_2JetReq->Integral(0,4)/egMetNew_2JetReq->Integral(0,4);
  //egMetNew_2JetReq->Scale(egScale_2JetReq,"");
  //egMet_2JetReqForLimit->Scale(FakeRate,"");
  fLimits_1jet.cd();
  ggMet_2JetReqForLimit->Write("met_gg_1jet");
  fin.cd();

  TH1F* ggMetNew_2JetReq=(TH1F*)ggMet_2JetReqForLimit->Rebin(NmetBins,"ggMetNew_2JetReq",xbins);
  //include overflow in last bin
  ggOverFlow=0.,ggOverFlowErr=0.;overflowBingg=ggMet_2JetReqForLimit->FindBin(ggMetNew_2JetReq->GetBinLowEdge(NmetBins+1)+.01);
  ggOverFlow=ggMet_2JetReqForLimit->IntegralAndError(overflowBingg,999,ggOverFlowErr);
  ggLastBinNew=ggMetNew_2JetReq->GetBinContent(NmetBins)+ggOverFlow;
  ggLastBinNewErr=sqrt(ggMetNew_2JetReq->GetBinError(NmetBins)*ggMetNew_2JetReq->GetBinError(NmetBins)+ggOverFlowErr*ggOverFlowErr);
  ggMetNew_2JetReq->SetBinContent(NmetBins,ggLastBinNew);
  ggMetNew_2JetReq->SetBinError(NmetBins,ggLastBinNewErr);

  TH1F* eeMet_reweightJet_binned_2JetReq=(TH1F*)fin.Get("eeMet_reweightJet_binned_2JetReq");
  TH1F* eeSidebandLowMet_reweightJet_binned_2JetReq = (TH1F*)fin.Get("eeSidebandLowMet_reweightJet_binned_2JetReq");
  TH1F* eeSidebandHighMet_reweightJet_binned_2JetReq = (TH1F*)fin.Get("eeSidebandHighMet_reweightJet_binned_2JetReq");
  /* TH1F* eeMetNew_2JetReq=(TH1F*)eeMet_reweightJet_binned_2JetReq->Clone();Rebin(NmetBins,"eeMetNew_2JetReq",xbins);
     TH1F* eeSidebandLowMetNew_2JetReq=(TH1F*)eeSidebandLowMet_reweightJet_binned_2JetReq->Rebin(NmetBins,"eeSidebandLowMetNew_2JetReq",xbins);
     TH1F* eeSidebandHighMetNew_2JetReq=(TH1F*)eeSidebandHighMet_reweightJet_binned_2JetReq->Rebin(NmetBins,"eeSidebandHighMetNew_2JetReq",xbins);*/
  //try unbinned

  TH1F* eeMetNew_2JetReq=(TH1F*)eeMet_reweightJet_binned_2JetReq->Clone();//;Rebin(NmetBins,"eeMetNew_2JetReq",xbins);
  TH1F* eeSidebandLowMetNew_2JetReq=(TH1F*)eeSidebandLowMet_reweightJet_binned_2JetReq->Clone();//Rebin(NmetBins,"eeSidebandLowMetNew_2JetReq",xbins);
  TH1F* eeSidebandHighMetNew_2JetReq=(TH1F*)eeSidebandHighMet_reweightJet_binned_2JetReq->Clone();//Rebin(NmetBins,"eeSidebandHighMetNew_2JetReq",xbins);


  //TH1F* eeSidebandMetNew_2JetReq=(TH1F*)eeSidebandMet_reweightJet_binned_2JetReq->Rebin(NmetBins,"eeSidebandMetNew_2JetReq",xbins);
  TH1F* eeMetMinusSideBand_2JetReq = new TH1F("eeMetMinusSideBand_2JetReq","",200,0.,1000.);//NmetBins,xbins);

  eeMetMinusSideBand_2JetReq->Sumw2();
  eeMetMinusSideBand_2JetReq->Add(eeMetNew_2JetReq,eeSidebandLowMetNew_2JetReq,1,-1);
  eeMetMinusSideBand_2JetReq->Add(eeSidebandHighMetNew_2JetReq,-1);
  eeMetMinusSideBand_2JetReq->GetXaxis()->SetRangeUser(0,400);
  eeMetMinusSideBand_2JetReq->GetYaxis()->SetRangeUser(1e-3,1.2e5);
  eeMetMinusSideBand_2JetReq->Draw();
  //eeMetFromZZ_2JetReq->SetLineColor(kRed);//temp
  eeMetFromZZ_JetReq->Draw("SAMES");//temp
  //eeMetFromWZ_2JetReq->SetLineColor(kBlue);//temp
  eeMetFromWZ_JetReq->Draw("SAMES");//temp
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/eeMetAndZZandWZ_2JetReq.png");
  eeMetMinusSideBand_2JetReq->Add(eeMetFromZZ_JetReq,-1);//temp
  eeMetMinusSideBand_2JetReq->Add(eeMetFromWZ_JetReq,-1);//temp
  TH1F* eeMetNew_ForLimits_2JetReq=(TH1F*)eeMetNew_2JetReq->Clone();
  TH1F* eeMetNew_ForLimitsStatOnly_2JetReq=(TH1F*)eeMetNew_2JetReq->Clone();
  TH1F* ffMetForee_2JetReq = (TH1F*)fin.Get("ffMet_reweightJet_binned_2JetReq");
  //TH1F* ffMetForee_2JetReq = (TH1F*)fin.Get("gfMet_reweightJet_binned_2JetReq");
  TH1F* egMet_ForLimits_2JetReq=(TH1F*)egMet_2JetReq->Clone();
  TH1F* egMetForff_2JetReq=(TH1F*)egMet_2JetReq->Clone();
  TH1F* eeMetMinusSideBand2_2JetReq=(TH1F*)eeMetNew_2JetReq->Clone();
  GetAndSetEEerrors(ggMet_2JetReq,egMet_2JetReq,egMet_ForLimits_2JetReq,eeMetNew_2JetReq,eeMetMinusSideBand2_2JetReq,eeMetNew_ForLimits_2JetReq,eeMetNew_ForLimitsStatOnly_2JetReq,eeSidebandHighMetNew_2JetReq,eeSidebandLowMetNew_2JetReq,ffMetForee_2JetReq,FakeRate,FakeRateErr,FakeRateErrForLimits,StdDev_2JetReq,StdDevSBHigh_2JetReq,StdDevSBLow_2JetReq,normErree_2JetReq,reweightErree_2JetReq,statErree_2JetReq,diffFromffErroree_2JetReq,eeScale_2JetReq,eeScaleErr_2JetReq,statErreg_2JetReq,normErreg_2JetReq,eeMetFromZZ_2JetReq,eeMetFromWZ_2JetReq);
 
  eeMetMinusSideBandForErrors_2JetReq=(TH1F*)eeMetMinusSideBand_2JetReq->Rebin(NmetBins,"eeMetMinusSideBandForErrors_2JetReq",xbins);
 
  //cout<<"Integral Over 100, ee 2JetReq.  from eeMetNew_2JetReq:"<<eeMetNew_2JetReq->Integral(eeMetNew_2JetReq->FindBin(100.1),-1)<<"  from eeMetMinusSideBand2:"<<eeMetMinusSideBand2_2JetReq->Integral(eeMetMinusSideBand2_2JetReq->FindBin(100.1),-1)<<endl;   

  fLimits_1jet.cd();
  egMet_ForLimits_2JetReq->Write("met_eg_1jet");
  if(doffgfcomb==false)eeMetNew_ForLimits_2JetReq->Write("met_ee_1jet");
  fin.cd();
 
  TH1F* egMetNew_2JetReq=(TH1F*)egMet_2JetReq->Rebin(NmetBins,"egMetNew_2JetReq",xbins);
  //include overflow in last bin
  egOverFlow=0.,egOverFlowErr=0.;overflowBineg=egMet_2JetReq->FindBin(egMetNew_2JetReq->GetBinLowEdge(NmetBins+1)+.01);
  egOverFlow=egMet_2JetReq->IntegralAndError(overflowBineg,999,egOverFlowErr);
  egLastBinNew=egMetNew_2JetReq->GetBinContent(NmetBins)+egOverFlow;
  egLastBinNewErr=sqrt(egMetNew_2JetReq->GetBinError(NmetBins)*egMetNew_2JetReq->GetBinError(NmetBins)+egOverFlowErr*egOverFlowErr);
  egMetNew_2JetReq->SetBinContent(NmetBins,egLastBinNew);
  egMetNew_2JetReq->SetBinError(NmetBins,egLastBinNewErr);
  //TH1F* eeMetMinusSideBand2_2JetReq=(TH1F*)eeMetNew_2JetReq->Rebin(NmetBins,"eeMetMinusSideBand2_2JetReq",xbins);

  //eeMetNew_2JetReq->Scale( ggMetNew_2JetReq->Integral(0,4)/eeMetNew_2JetReq->Integral(0,4) , "");
  TH1F *totalBG_2JetReq=new TH1F("totalBG_2JetReq","",NmetBins,xbins);
  totalBG_2JetReq->Sumw2();
  totalBG_2JetReq->Add(egMetNew_2JetReq,eeMetMinusSideBand2_2JetReq,1,1);
   
  gg0_20_2JetReq=ggMetNew_2JetReq->IntegralAndError(0,4,gg0_20Error_2JetReq);
  eg0_20_2JetReq=egMetNew_2JetReq->IntegralAndError(0,4,eg0_20Error_2JetReq);
  ee0_20_2JetReq=eeMetMinusSideBand2_2JetReq->IntegralAndError(0,4,ee0_20Error_2JetReq);
  QCDee0_20_2JetReq=totalBG_2JetReq->IntegralAndError(0,4,QCDee0_20Error_2JetReq);
  gg30_50_2JetReq=ggMetNew_2JetReq->IntegralAndError(ggMetNew_2JetReq->FindBin(30),ggMetNew_2JetReq->FindBin(50-1),gg30_50Error_2JetReq);
  eg30_50_2JetReq=egMetNew_2JetReq->IntegralAndError(egMetNew_2JetReq->FindBin(30),egMetNew_2JetReq->FindBin(50-1),eg30_50Error_2JetReq);
  ee30_50_2JetReq=eeMetMinusSideBand2_2JetReq->IntegralAndError(eeMetMinusSideBand2_2JetReq->FindBin(30),eeMetMinusSideBand2_2JetReq->FindBin(50-1),ee30_50Error_2JetReq);
  QCDee30_50_2JetReq=totalBG_2JetReq->IntegralAndError(totalBG_2JetReq->FindBin(30),totalBG_2JetReq->FindBin(50-1),QCDee30_50Error_2JetReq);
  gg50up_2JetReq=ggMetNew_2JetReq->IntegralAndError(ggMetNew_2JetReq->FindBin(50),-1,gg50upError_2JetReq);
  eg50up_2JetReq=egMetNew_2JetReq->IntegralAndError(egMetNew_2JetReq->FindBin(50),-1,eg50upError_2JetReq);
  ee50up_2JetReq=eeMetMinusSideBand2_2JetReq->IntegralAndError(eeMetMinusSideBand2_2JetReq->FindBin(50),-1,ee50upError_2JetReq);
  QCDee50up_2JetReq=totalBG_2JetReq->IntegralAndError(totalBG_2JetReq->FindBin(50),-1,QCDee50upError_2JetReq);
  gg100up_2JetReq=ggMetNew_2JetReq->IntegralAndError(ggMetNew_2JetReq->FindBin(100),-1,gg100upError_2JetReq);
  eg100up_2JetReq=egMetNew_2JetReq->IntegralAndError(egMetNew_2JetReq->FindBin(100),-1,eg100upError_2JetReq);
  ee100up_2JetReq=eeMetMinusSideBand2_2JetReq->IntegralAndError(eeMetMinusSideBand2_2JetReq->FindBin(100),-1,ee100upError_2JetReq);
  QCDee100up_2JetReq=totalBG_2JetReq->IntegralAndError(totalBG_2JetReq->FindBin(100),-1,QCDee100upError_2JetReq);
  ee50_60_2JetReq=eeMetMinusSideBand2_2JetReq->IntegralAndError(eeMetMinusSideBand2_2JetReq->FindBin(50),eeMetMinusSideBand2_2JetReq->FindBin(60-1),ee50_60Error_2JetReq);
  QCDee50_60_2JetReq=totalBG_2JetReq->IntegralAndError(totalBG_2JetReq->FindBin(50),totalBG_2JetReq->FindBin(60-1),QCDee50_60Error_2JetReq);
  ee60_70_2JetReq=eeMetMinusSideBand2_2JetReq->IntegralAndError(eeMetMinusSideBand2_2JetReq->FindBin(60),eeMetMinusSideBand2_2JetReq->FindBin(70-1),ee60_70Error_2JetReq);
  QCDee60_70_2JetReq=totalBG_2JetReq->IntegralAndError(totalBG_2JetReq->FindBin(60),totalBG_2JetReq->FindBin(70-1),QCDee60_70Error_2JetReq);
  ee70_80_2JetReq=eeMetMinusSideBand2_2JetReq->IntegralAndError(eeMetMinusSideBand2_2JetReq->FindBin(70),eeMetMinusSideBand2_2JetReq->FindBin(80-1),ee70_80Error_2JetReq);
  QCDee70_80_2JetReq=totalBG_2JetReq->IntegralAndError(totalBG_2JetReq->FindBin(70),totalBG_2JetReq->FindBin(80-1),QCDee70_80Error_2JetReq);
  ee80_100_2JetReq=eeMetMinusSideBand2_2JetReq->IntegralAndError(eeMetMinusSideBand2_2JetReq->FindBin(80),eeMetMinusSideBand2_2JetReq->FindBin(100-1),ee80_100Error_2JetReq);
  QCDee80_100_2JetReq=totalBG_2JetReq->IntegralAndError(totalBG_2JetReq->FindBin(80),totalBG_2JetReq->FindBin(100-1),QCDee80_100Error_2JetReq);
  gg50_60_2JetReq=ggMetNew_2JetReq->IntegralAndError(ggMetNew_2JetReq->FindBin(50),ggMetNew_2JetReq->FindBin(60-1),gg50_60Error_2JetReq);
  gg60_70_2JetReq=ggMetNew_2JetReq->IntegralAndError(ggMetNew_2JetReq->FindBin(60),ggMetNew_2JetReq->FindBin(70-1),gg60_70Error_2JetReq);
  gg70_80_2JetReq=ggMetNew_2JetReq->IntegralAndError(ggMetNew_2JetReq->FindBin(70),ggMetNew_2JetReq->FindBin(80-1),gg70_80Error_2JetReq);
  gg80_100_2JetReq=ggMetNew_2JetReq->IntegralAndError(ggMetNew_2JetReq->FindBin(80),ggMetNew_2JetReq->FindBin(100-1),gg80_100Error_2JetReq);
  eg50_60_2JetReq=egMetNew_2JetReq->IntegralAndError(egMetNew_2JetReq->FindBin(50),egMetNew_2JetReq->FindBin(60-1),eg50_60Error_2JetReq);
  eg60_70_2JetReq=egMetNew_2JetReq->IntegralAndError(egMetNew_2JetReq->FindBin(60),egMetNew_2JetReq->FindBin(70-1),eg60_70Error_2JetReq);
  eg70_80_2JetReq=egMetNew_2JetReq->IntegralAndError(egMetNew_2JetReq->FindBin(70),egMetNew_2JetReq->FindBin(80-1),eg70_80Error_2JetReq);
  eg80_100_2JetReq=egMetNew_2JetReq->IntegralAndError(egMetNew_2JetReq->FindBin(80),egMetNew_2JetReq->FindBin(100-1),eg80_100Error_2JetReq);
  
  TH1F *eeOverBinWidth_2JetReq = new TH1F("eeOverBinWidth_2JetReq","",NmetBins,xbins);eeOverBinWidth_2JetReq->Sumw2();
  TH1F *egOverBinWidth_2JetReq = new TH1F("egOverBinWidth_2JetReq","",NmetBins,xbins);egOverBinWidth_2JetReq->Sumw2();
  TH1F *ggOverBinWidth_2JetReq = new TH1F("ggOverBinWidth_2JetReq","",NmetBins,xbins);ggOverBinWidth_2JetReq->Sumw2();
  for(int i=1;i<NmetBins+1;i++){
    float eg = egMetNew_2JetReq->GetBinContent(i)/egMetNew_2JetReq->GetBinWidth(i);
    float ee = eeMetMinusSideBand2_2JetReq->GetBinContent(i)/eeMetMinusSideBand2_2JetReq->GetBinWidth(i);
    float gg = ggMetNew_2JetReq->GetBinContent(i)/ggMetNew_2JetReq->GetBinWidth(i);
    float ggE = ggMetNew_2JetReq->GetBinError(i)/ggMetNew_2JetReq->GetBinWidth(i);
    float bg = totalBG_2JetReq->GetBinContent(i)/totalBG_2JetReq->GetBinWidth(i);
    float bgE = totalBG_2JetReq->GetBinError(i)/totalBG_2JetReq->GetBinWidth(i);
    egOverBinWidth_2JetReq->SetBinContent(i,eg);
    eeOverBinWidth_2JetReq->SetBinContent(i,ee);
    ggOverBinWidth_2JetReq->SetBinContent(i,gg);
    ggOverBinWidth_2JetReq->SetBinError(i,ggE);
    totalBG_2JetReq->SetBinContent(i,bg);
    totalBG_2JetReq->SetBinError(i,bgE);
    float sig_2JetReq=ggSig2012Rebin_2JetReq->GetBinContent(i)/ggSig2012Rebin_2JetReq->GetBinWidth(i);
    ggSig2012Rebin_2JetReq->SetBinContent(i,sig_2JetReq);
    float sig_2JetReq2=ggSig2012Rebin_2JetReq_2->GetBinContent(i)/ggSig2012Rebin_2JetReq_2->GetBinWidth(i);
    ggSig2012Rebin_2JetReq_2->SetBinContent(i,sig_2JetReq2);
  }
  ggOverBinWidth_2JetReq->SetMarkerColor(kBlack);
  ggOverBinWidth_2JetReq->SetLineColor(kBlack);
  ggOverBinWidth_2JetReq->SetLineWidth(2.2);
  ggOverBinWidth_2JetReq->SetMarkerSize(0.6);
  ggOverBinWidth_2JetReq->SetStats(0);    
  egOverBinWidth_2JetReq->SetMaximum(Max);
  egOverBinWidth_2JetReq->GetYaxis()->SetTitle("Number of Events / GeV");
  egOverBinWidth_2JetReq->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  egOverBinWidth_2JetReq->SetTitle("");
  egOverBinWidth_2JetReq->SetLineColor(kAzure);
  egOverBinWidth_2JetReq->SetFillColor(kAzure);
  //egOverBinWidth_2JetReq->SetFillStyle(3005);
  egOverBinWidth_2JetReq->SetStats(0);
  eeOverBinWidth_2JetReq->SetLineColor(kGray+1);
  eeOverBinWidth_2JetReq->SetFillColor(kGray+1);
  //eeOverBinWidth_2JetReq->SetFillStyle(3004);
  eeOverBinWidth_2JetReq->SetMarkerSize(0);
  eeOverBinWidth_2JetReq->SetStats(0);
  eeOverBinWidth_2JetReq->SetTitle("");
  totalBG_2JetReq->SetMarkerSize(0);
  totalBG_2JetReq->SetFillColor(kRed);
  totalBG_2JetReq->SetFillStyle(3154);
  totalBG_2JetReq->SetStats(0);
  THStack *eeStack_2JetReq = new THStack("eeStack_2JetReq",";;Number of Events / GeV");
  eeStack_2JetReq->SetMaximum(Max);
  eeStack_2JetReq->SetMinimum(Min);
  eeStack_2JetReq->Add(egOverBinWidth_2JetReq);
  eeStack_2JetReq->Add(eeOverBinWidth_2JetReq);
  eeStack_2JetReq->Draw();
  totalBG_2JetReq->Draw("E2SAME");
  ggOverBinWidth_2JetReq->Draw("PESAME");
  TLegend *legEE_2JetReq = new TLegend(.5,.45,.8,.65,"","brNDC");
  legEE_2JetReq->AddEntry(ggOverBinWidth_2JetReq,"#gamma#gamma Candidate Sample","lpf");
  legEE_2JetReq->AddEntry(totalBG_2JetReq,"QCD + Electroweak Error","f");
  legEE_2JetReq->AddEntry(eeOverBinWidth_2JetReq,"QCD (e^{+}e^{-} sample)","f");
  legEE_2JetReq->AddEntry(egOverBinWidth_2JetReq,"Electroweak","f");
  legEE_2JetReq->SetFillColor(kWhite);
  legEE_2JetReq->Draw("SAME");
  TPaveText *Text_2JetReq;
  Text_2JetReq= new TPaveText(.5,.66,.8,.87,"NDC");
  Text_2JetReq->AddText("CMS Preliminary");
  Text_2JetReq->AddText("");
  Text_2JetReq->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  Text_2JetReq->AddText(">=1 Jet Requirement");
  Text_2JetReq->SetFillColor(0);
  Text_2JetReq->SetBorderSize(0);
  Text_2JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_ggANDeeBinnedReweighted_2JetReq.png");
  fout.cd();eeOverBinWidth_2JetReq->Write("eeMetBinnedReweightedAndScaled_2JetReq");ggOverBinWidth_2JetReq->Write("ggMetReweightedAndScaled_2JetReq");fin.cd();

  TH1F* ggOveree_2JetReq=(TH1F*)ggOverBinWidth_2JetReq->Clone();ggOveree_2JetReq->SetTitle("ggOveree_2JetReq");ggOveree_2JetReq->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  ggOveree_2JetReq->Divide(eeOverBinWidth_2JetReq);
  c1->SetLogy(0);
  ggOveree_2JetReq->Draw();
  line1->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggOveree_2JetReq.png");
  fout.cd();ggOveree->Write("ggOveree_2JetReq.png");fin.cd();
  totalBGee_2JetReq = (TH1F*)totalBG_2JetReq->Clone();
  TH1F* ggOvereeQCD_2JetReq=(TH1F*)ggOverBinWidth_2JetReq->Clone();ggOvereeQCD_2JetReq->SetTitle("ggOvereeQCD_2JetReq");ggOvereeQCD_2JetReq->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  TH1F* ggOvereeQCDerr_2JetReq=(TH1F*)totalBG_2JetReq->Clone();
  for(int i=1;i<ggOvereeQCD_2JetReq->GetNbinsX()+1;i++){
    float Value = ggOvereeQCD_2JetReq->GetBinContent(i);float StatErr = ggOvereeQCD_2JetReq->GetBinError(i);
    Value/=totalBG_2JetReq->GetBinContent(i);StatErr/=totalBG_2JetReq->GetBinContent(i);
    ggOvereeQCD_2JetReq->SetBinContent(i,Value);ggOvereeQCD_2JetReq->SetBinError(i,StatErr);
    float SystErr = totalBG_2JetReq->GetBinError(i)/totalBG_2JetReq->GetBinContent(i);
    ggOvereeQCDerr_2JetReq->SetBinContent(i,1.);ggOvereeQCDerr_2JetReq->SetBinError(i,SystErr);
  }
  //ggOvereeQCD->Divide(totalBG);
  c1->SetLogy(0);
  ggOvereeQCDerr_2JetReq->SetFillColor(kRed);
  ggOvereeQCDerr_2JetReq->GetYaxis()->SetRangeUser(0.0,2.0);
  ggOvereeQCDerr_2JetReq->Draw("E2");
  ggOvereeQCD_2JetReq->Draw("SAME");
  line1->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggOvereeQCD_2JetReq.png");
  fout.cd();ggOvereeQCD_2JetReq->Write("ggOvereeQCD_2JetReq");ggOvereeQCDerr_2JetReq->Write("ggOvereeQCDerr_2JetReq");fin.cd();
  c1->SetLogy(1);

  c2->cd();
  p1->cd();
  p1->SetLogy(1);
  /* egOverBinWidth_2JetReq->Draw("B");
     eeOverBinWidth_2JetReq->Draw("BSAME");
     eeOverBinWidth_2JetReq->Draw("SAME");
     egOverBinWidth_2JetReq->Draw("BSAME");
     egOverBinWidth_2JetReq->Draw("SAME");*/
  eeStack_2JetReq->Draw();
  totalBG_2JetReq->Draw("E2SAME");
  ggOverBinWidth_2JetReq->Draw("PESAME");
  ggSig2012Rebin_2JetReq->Draw("histo SAME");
  ggSig2012Rebin_2JetReq_2->Draw("histo SAME");
  TLegend *legEEwr_2JetReq = new TLegend(.5,.35,.83,.61,"","brNDC");
  legEEwr_2JetReq->AddEntry(ggOverBinWidth_2JetReq,"#gamma#gamma Candidate Sample","lpf");
  legEEwr_2JetReq->AddEntry(totalBG_2JetReq,"QCD + Electroweak Error","f");
  legEEwr_2JetReq->AddEntry(eeOverBinWidth_2JetReq,"QCD (e^{+}e^{-} sample)","f");
  legEEwr_2JetReq->AddEntry(egOverBinWidth_2JetReq,"Electroweak","f");
  legEEwr_2JetReq->AddEntry(ggSig2012Rebin_2JetReq_2,"GGM #gamma#gamma (1100_720_375)","l");
  legEEwr_2JetReq->AddEntry(ggSig2012Rebin_2JetReq,"GGM #gamma#gamma (1400_1720_375)","l");
  legEEwr_2JetReq->SetFillColor(kWhite);
  legEEwr_2JetReq->Draw("SAME");
  TPaveText *Textwr_2JetReq;
  Textwr_2JetReq = new TPaveText(.5,.62,.8,.83,"NDC");
  Textwr_2JetReq->AddText("CMS Preliminary");
  Textwr_2JetReq->AddText("");
  Textwr_2JetReq->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  Textwr_2JetReq->AddText(">=1 Jet Requirement");
  Textwr_2JetReq->SetFillColor(0);
  Textwr_2JetReq->SetBorderSize(0);
  Textwr_2JetReq->Draw();
  p2->cd();
  ggOvereeQCDerr_2JetReq->SetTitle("");
  ggOvereeQCDerr_2JetReq->GetYaxis()->SetTitle("Data/Prediction");
  ggOvereeQCDerr_2JetReq->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  ggOvereeQCDerr_2JetReq->GetYaxis()->SetTitleOffset(0.48);
  ggOvereeQCDerr_2JetReq->GetYaxis()->SetTitleSize(0.15);
  ggOvereeQCDerr_2JetReq->GetXaxis()->SetTitleSize(0.2);
  ggOvereeQCDerr_2JetReq->GetYaxis()->SetLabelSize(0.12);
  ggOvereeQCDerr_2JetReq->GetXaxis()->SetLabelSize(0.15);
  ggOvereeQCDerr_2JetReq->GetYaxis()->SetRangeUser(ratMin,ratMax);
  //ggOvereeQCD->GetYaxis()->SetTitleSize(3);
  ggOvereeQCDerr_2JetReq->Draw("E2");
  ggOvereeQCD_2JetReq->Draw("SAME");
  line1->Draw("SAME");
  p3->cd();
  met->Draw();
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDeeBinnedReweighted_WithRatio_2JetReq.png");
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDeeBinnedReweighted_WithRatio_2JetReq.pdf");
  c1->cd();
  //---------Now ff---------
  //Get Errors From Toys
  //float binMin[12]={0.};
  //float binMax[12]={0.};
  float StdDevff[HistBins];//={0.};
  float meanff[HistBins];//={0.};
    
  for(int i=0;i<HistBins;i++){
    StdDevff[i]=0.;meanff[i]=0.;
  }

  TString titleff="Toys/ffMet_toy_1";
  hist=(TH1F*)fin.Get(titleff);hist->Sumw2();
    
  for(int bin=1;bin<HistBins+1;bin++){
    for(int i=1;i<1001;i++){
      titleff="Toys/ffMet_toy_";titleff+=i;
      hist=(TH1F*)fin.Get(titleff);if(hist->GetSumw2()==0)hist->Sumw2();
      meanff[bin-1]+=hist->GetBinContent(bin);
      //cout<<hist->GetBinContent(bin)<<endl;
    }
    //cout<<"Bin:"<<bin<<"  MeanSum:"<<meanff[bin-1]<<"  Meanff:"<<meanff[bin-1]/1000<<endl;
    meanff[bin-1]=meanff[bin-1]/1000;
  }
    
  for(int i=1;i<1001;i++){
    titleff="Toys/ffMet_toy_";titleff+=i;
    //cout<<"titleff: "<<titleff<<endl;
    //fin.Get(titleff)->Draw("SAME");
    hist = (TH1F*)fin.Get(titleff);
    for(int bin=1;bin<HistBins+1;bin++){
      //if(hist->GetBinContent(bin)<binMin[bin-1])binMin[bin-1]=hist->GetBinContent(bin);
      //if(hist->GetBinContent(bin)>binMax[bin-1])binMax[bin-1]=hist->GetBinContent(bin);
      StdDevff[bin-1]+= (hist->GetBinContent(bin)-meanff[bin-1])*(hist->GetBinContent(bin)-meanff[bin-1]);
      //if(bin==1)cout<<"StdDev: "<<StdDevff[bin-1]<<endl;
    }
  }
  for(int bin=1;bin<HistBins+1;bin++){
    StdDevff[bin-1]=sqrt(StdDevff[bin-1]/1000);
    //cout<<"Bin:" <<bin<<endl<<" BinMin:"<<binMin[bin-1]<<endl<< " BinMax:"<<binMax[bin-1]<<endl<<" STDDEV:"<<StdDev[bin-1]<<endl;
  }
  //gammafake no jet
  float StdDevgammafake[HistBins];//={0.};
  float meangammafake[HistBins];//={0.};
    
  for(int i=0;i<HistBins;i++){
    StdDevgammafake[i]=0.;meangammafake[i]=0.;
  }

  TString titlegammafake="Toys/gammafakeMet_toy_1";
  hist=(TH1F*)fin.Get(titlegammafake);hist->Sumw2();
    
  for(int bin=1;bin<HistBins+1;bin++){
    for(int i=1;i<1001;i++){
      titlegammafake="Toys/gammafakeMet_toy_";titlegammafake+=i;
      hist=(TH1F*)fin.Get(titlegammafake);if(hist->GetSumw2()==0)hist->Sumw2();
      meangammafake[bin-1]+=hist->GetBinContent(bin);
      //cout<<hist->GetBinContent(bin)<<endl;
    }
    //cout<<"Bin:"<<bin<<"  MeanSum:"<<meangammafake[bin-1]<<"  Meangammafake:"<<meangammafake[bin-1]/1000<<endl;
    meangammafake[bin-1]=meangammafake[bin-1]/1000;
  }
    
  for(int i=1;i<1001;i++){
    titlegammafake="Toys/gammafakeMet_toy_";titlegammafake+=i;
    //cout<<"titlegammafake: "<<titlegammafake<<endl;
    //fin.Get(titlegammafake)->Draw("SAME");
    hist = (TH1F*)fin.Get(titlegammafake);
    for(int bin=1;bin<HistBins+1;bin++){
      //if(hist->GetBinContent(bin)<binMin[bin-1])binMin[bin-1]=hist->GetBinContent(bin);
      //if(hist->GetBinContent(bin)>binMax[bin-1])binMax[bin-1]=hist->GetBinContent(bin);
      StdDevgammafake[bin-1]+= (hist->GetBinContent(bin)-meangammafake[bin-1])*(hist->GetBinContent(bin)-meangammafake[bin-1]);
      //if(bin==1)cout<<"StdDev: "<<StdDevgammafake[bin-1]<<endl;
    }
  }
  for(int bin=1;bin<HistBins+1;bin++){
    StdDevgammafake[bin-1]=sqrt(StdDevgammafake[bin-1]/1000);
    //cout<<"Bin:" <<bin<<endl<<" BinMin:"<<binMin[bin-1]<<endl<< " BinMax:"<<binMax[bin-1]<<endl<<" STDDEV:"<<StdDev[bin-1]<<endl;
  }
  //done gammafake
  //gammafake ff combination
  float StdDevgammafakeffcomb[HistBins];//={0.};
  for(int bin=1;bin<HistBins+1;bin++){
    StdDevgammafakeffcomb[bin-1]=sqrt(gammafakePercent*gammafakePercent*StdDevgammafake[bin-1]*StdDevgammafake[bin-1]+ffPercent*ffPercent*StdDevff[bin-1]*StdDevff[bin-1]);
    //cout<<"Bin:" <<bin<<endl<<" BinMin:"<<binMin[bin-1]<<endl<< " BinMax:"<<binMax[bin-1]<<endl<<" STDDEV:"<<StdDev[bin-1]<<endl;
  }
  //end comb

  titleff="Toys/ffMet_toy_1";
  hist = (TH1F*)fin.Get(titleff);
  hist->SetTitle("ffMet Toys");
  hist->Draw("histo");
  for(int i=2;i<1001;i++){
    titleff="Toys/ffMet_toy_";titleff+=i;
    //cout<<"titleff: "<<titleff<<endl;
    //fin.Get(titleff)->Draw("SAME");
    hist = (TH1F*)fin.Get(titleff);
    hist->Draw("histo SAME");
  }
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MetToys_ff.png");

  titleff="Toys/ff_JetReqMet_toy_1";
  hist = (TH1F*)fin.Get(titleff);
  hist->SetTitle("ffMet Toys JetReq");
  hist->Draw("histo");
  for(int i=2;i<1001;i++){
    titleff="Toys/ff_JetReqMet_toy_";titleff+=i;
    //cout<<"titleff: "<<titleff<<endl;
    //fin.Get(titleff)->Draw("SAME");
    hist = (TH1F*)fin.Get(titleff);
    hist->Draw("histo SAME");
  }
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MetToys_ff_JetReq.png");

  titleff="Toys/ff_2JetReqMet_toy_1";
  hist = (TH1F*)fin.Get(titleff);
  hist->SetTitle("ffMet Toys 2JetReq");
  hist->Draw("histo");
  for(int i=2;i<1001;i++){
    titleff="Toys/ff_2JetReqMet_toy_";titleff+=i;
    //cout<<"titleff: "<<titleff<<endl;
    //fin.Get(titleff)->Draw("SAME");
    hist = (TH1F*)fin.Get(titleff);
    hist->Draw("histo SAME");
  }
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MetToys_ff_2JetReq.png");


  fin.cd();
  /*
    TH1F* ffMet_reweightJet_binned1 = (TH1F*)fin.Get("ffMet_reweightJet_binned");ffMet_reweightJet_binned1->Sumw2();
    TH1F* ffMet_reweightJet_binned2 = (TH1F*)fin.Get("gfMet_reweightJet_binned");ffMet_reweightJet_binned2->Sumw2();
    TH1F* ffMet_reweightJet_binned = (TH1F*)ffMet_reweightJet_binned1->Clone();
    ffMet_reweightJet_binned->Add(ffMet_reweightJet_binned2);ffMet_reweightJet_binned->Scale(.5);*/
  TH1F* fgMet_reweightJet_binnedComp = (TH1F*)fin.Get("fgMet_reweightJet_binned");fgMet_reweightJet_binnedComp->Sumw2();
  TH1F* gfMet_reweightJet_binnedComp = (TH1F*)fin.Get("gfMet_reweightJet_binned");gfMet_reweightJet_binnedComp->Sumw2();
  TH1F* gammafakeMet_reweightJet_binnedComp = (TH1F*)fin.Get("gammafakeMet_reweightJet_binned");gammafakeMet_reweightJet_binnedComp->Sumw2();
  TH1F* ffMet_reweightJet_binnedComp = (TH1F*)fin.Get("ffMet_reweightJet_binned");ffMet_reweightJet_binnedComp->Sumw2();
  float scale = ffMet_reweightJet_binnedComp->Integral()/gfMet_reweightJet_binnedComp->Integral();
  gfMet_reweightJet_binnedComp->Scale(scale);gfMet_reweightJet_binnedComp->SetLineColor(kRed);gfMet_reweightJet_binnedComp->SetMarkerColor(kRed);
  scale = ffMet_reweightJet_binnedComp->Integral()/fgMet_reweightJet_binnedComp->Integral();
  fgMet_reweightJet_binnedComp->Scale(scale);fgMet_reweightJet_binnedComp->SetLineColor(kBlue);fgMet_reweightJet_binnedComp->SetMarkerColor(kBlue);
  scale = ffMet_reweightJet_binnedComp->Integral()/gammafakeMet_reweightJet_binnedComp->Integral();
  gammafakeMet_reweightJet_binnedComp->Scale(scale);gammafakeMet_reweightJet_binnedComp->SetLineColor(kGreen);gammafakeMet_reweightJet_binnedComp->SetMarkerColor(kGreen);
  ffMet_reweightJet_binnedComp->GetXaxis()->SetRangeUser(0,99);
  ffMet_reweightJet_binnedComp->Draw("PE");gfMet_reweightJet_binnedComp->Draw("PEsames");fgMet_reweightJet_binnedComp->Draw("PEsames");gammafakeMet_reweightJet_binnedComp->Draw("PEsames");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_ffANDgfANDfgBinnedReweighted.png");
  TH1F* fgMet_reweightJet_binnedCompRebin = (TH1F*)fgMet_reweightJet_binnedComp->Rebin(NmetBins,"fgMet_reweightJet_binnedCompRebin",xbins);
  TH1F* ffMet_reweightJet_binnedCompRebin = (TH1F*)ffMet_reweightJet_binnedComp->Rebin(NmetBins,"ffMet_reweightJet_binnedCompRebin",xbins);
  TH1F* gfMet_reweightJet_binnedCompRebin = (TH1F*)gfMet_reweightJet_binnedComp->Rebin(NmetBins,"gfMet_reweightJet_binnedCompRebin",xbins);
  TH1F* gammafakeMet_reweightJet_binnedCompRebin = (TH1F*)gammafakeMet_reweightJet_binnedComp->Rebin(NmetBins,"gammafakeMet_reweightJet_binnedCompRebin",xbins);
  for(int i=1;i<ffMet_reweightJet_binnedCompRebin->GetNbinsX()+1;i++){
    float Value = ffMet_reweightJet_binnedCompRebin->GetBinContent(i);Value/=ffMet_reweightJet_binnedCompRebin->GetBinWidth(i);
    ffMet_reweightJet_binnedCompRebin->SetBinContent(i,Value);
    Value = ffMet_reweightJet_binnedCompRebin->GetBinError(i);Value/=ffMet_reweightJet_binnedCompRebin->GetBinWidth(i);
    ffMet_reweightJet_binnedCompRebin->SetBinError(i,Value);
    Value = gfMet_reweightJet_binnedCompRebin->GetBinContent(i);Value/=gfMet_reweightJet_binnedCompRebin->GetBinWidth(i);
    gfMet_reweightJet_binnedCompRebin->SetBinContent(i,Value);
    Value = gfMet_reweightJet_binnedCompRebin->GetBinError(i);Value/=gfMet_reweightJet_binnedCompRebin->GetBinWidth(i);
    gfMet_reweightJet_binnedCompRebin->SetBinError(i,Value);
    Value = fgMet_reweightJet_binnedCompRebin->GetBinContent(i);Value/=fgMet_reweightJet_binnedCompRebin->GetBinWidth(i);
    fgMet_reweightJet_binnedCompRebin->SetBinContent(i,Value);
    Value = fgMet_reweightJet_binnedCompRebin->GetBinError(i);Value/=fgMet_reweightJet_binnedCompRebin->GetBinWidth(i);
    fgMet_reweightJet_binnedCompRebin->SetBinError(i,Value);
    Value = gammafakeMet_reweightJet_binnedCompRebin->GetBinContent(i);Value/=gammafakeMet_reweightJet_binnedCompRebin->GetBinWidth(i);
    gammafakeMet_reweightJet_binnedCompRebin->SetBinContent(i,Value);
    Value = gammafakeMet_reweightJet_binnedCompRebin->GetBinError(i);Value/=gammafakeMet_reweightJet_binnedCompRebin->GetBinWidth(i);
    gammafakeMet_reweightJet_binnedCompRebin->SetBinError(i,Value);
  }
  ffMet_reweightJet_binnedCompRebin->Draw("PE");gfMet_reweightJet_binnedCompRebin->Draw("PEsames");fgMet_reweightJet_binnedCompRebin->Draw("PEsames");gammafakeMet_reweightJet_binnedCompRebin->Draw("PEsames");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_ffANDgfANDfgBinnedReweightedRebin.png");
  gfMet_reweightJet_binnedCompRebin->Divide(ffMet_reweightJet_binnedCompRebin);
  fgMet_reweightJet_binnedCompRebin->Divide(ffMet_reweightJet_binnedCompRebin);
  gammafakeMet_reweightJet_binnedCompRebin->Divide(ffMet_reweightJet_binnedCompRebin);
  c1->SetLogy(0);
  gfMet_reweightJet_binnedCompRebin->Draw("PE");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_gfOVERff.png");
  fgMet_reweightJet_binnedCompRebin->Draw("PE");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_fgOVERff.png");
  gammafakeMet_reweightJet_binnedCompRebin->Draw("PE");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_gammafakeOVERff.png");
  gfMet_reweightJet_binnedCompRebin->Draw("PE");
  fgMet_reweightJet_binnedCompRebin->Draw("PEsames");
  gammafakeMet_reweightJet_binnedCompRebin->Draw("PEsames");
  line1->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_gammafakesOVERff.png");

  ffMet_reweightJet_binnedComp->GetXaxis()->SetRangeUser(0,1000);

  c1->SetLogy(1);
  TH1F* ffMet_reweightJet_binned = (TH1F*)fin.Get("ffMet_reweightJet_binned");ffMet_reweightJet_binned->Sumw2();
  TH1F* ffMet_reweightJet_binnedSyst = (TH1F*)ffMet_reweightJet_binned->Clone();ffMet_reweightJet_binnedSyst->Sumw2();
  TH1F* gfMet_reweightJet_binned = (TH1F*)fin.Get("gfMet_reweightJet_binned");gfMet_reweightJet_binned->Sumw2();
  TH1F* gfMet_reweightJet_binnedSyst = (TH1F*)gfMet_reweightJet_binned->Clone();gfMet_reweightJet_binnedSyst->Sumw2();
  TH1F* fgMet_reweightJet_binned = (TH1F*)fin.Get("fgMet_reweightJet_binned");fgMet_reweightJet_binned->Sumw2();
  TH1F* gammafakeMet_reweightJet_binned = (TH1F*)fin.Get("gammafakeMet_reweightJet_binned");gammafakeMet_reweightJet_binned->Sumw2();
  TH1F* Met_reweightJet_binned_ffSyst = (TH1F*)gammafakeMet_reweightJet_binned->Clone();Met_reweightJet_binned_ffSyst->Sumw2();
  ffMet_reweightJet_binnedForErrors=(TH1F*)ffMet_reweightJet_binned->Rebin(NmetBins,"ffMet_reweightJet_binnedForErrors",xbins);
  TH1F* ffMet_ForLimits = (TH1F*)ffMet_reweightJet_binned->Clone();
  TH1F* ffMet_ForLimitsStatOnly = (TH1F*)ffMet_reweightJet_binned->Clone();
  ffMetNew = (TH1F*)ffMet_reweightJet_binned->Clone();
  TH1F *eeOverBinWidthForFF=(TH1F*)eeMetMinusSideBand->Clone();eeOverBinWidthForFF->Sumw2();

  //float ffCont = 0.29412,gammafakeCont = 0.70588;
  float ffCont = 0.26,gammafakeCont = 0.46;
  //float ffPercent = ffCont/(ffCont+gammafakeCont),gammafakePercent = gammafakeCont/(ffCont+gammafakeCont);
  //float ffPercent=0.3635,gammafakePercent=0.6365;//these numbers are from jet req min chi2 fraction fit
  if(doffgfcomb){
    //this for systematic error, grab ff
    GetAndSetFFerrors(ggMetData,egMetForff,Met_reweightJet_binned_ffSyst/*eeMetMinusSideBand/*,eeOverBinWidthForFF*/,ffMet_reweightJet_binnedSyst,ffMetNew,ffMet_ForLimits,ffMet_ForLimitsStatOnly,FakeRate,FakeRateErr,FakeRateErrForLimits,StdDevff,normErrff,reweightErrff,statErrff,diffFromeeErrorff,diffFromfgErrorff,0,ffScale,ffScaleErr);
    fLimits_nojet.cd();
    ffMet_ForLimits->Write("met_ee_nojet");
    fin.cd();
    ee0_20=ffMetNew->IntegralAndError(0,4,ee0_20Error);
    ee30_50=ffMetNew->IntegralAndError(ffMetNew->FindBin(30),ffMetNew->FindBin(50-1),ee30_50Error);
    ee50up=ffMetNew->IntegralAndError(ffMetNew->FindBin(50),-1,ee50upError);
    ee100up=ffMetNew->IntegralAndError(ffMetNew->FindBin(100),-1,ee100upError);
    ee50_60=ffMetNew->IntegralAndError(ffMetNew->FindBin(50),ffMetNew->FindBin(60-1),ee50_60Error);
    ee60_70=ffMetNew->IntegralAndError(ffMetNew->FindBin(60),ffMetNew->FindBin(70-1),ee60_70Error);
    ee70_80=ffMetNew->IntegralAndError(ffMetNew->FindBin(70),ffMetNew->FindBin(80-1),ee70_80Error);
    ee80_100=ffMetNew->IntegralAndError(ffMetNew->FindBin(80),ffMetNew->FindBin(100-1),ee80_100Error);
    TH1F *totalBGtemp=new TH1F("totalBGtemp","",NmetBins,xbins);totalBGtemp->Sumw2();
    totalBGtemp->Add(ffMetNew,egMetNew,1,1);
    QCDee0_20=totalBGtemp->IntegralAndError(0,4,QCDee0_20Error);
    QCDee30_50=totalBGtemp->IntegralAndError(totalBGtemp->FindBin(30),totalBGtemp->FindBin(50-1),QCDee30_50Error);
    QCDee50up=totalBGtemp->IntegralAndError(totalBGtemp->FindBin(50),-1,QCDee50upError);
    QCDee100up=totalBGtemp->IntegralAndError(totalBGtemp->FindBin(100),-1,QCDee100upError);
    QCDee50_60=totalBGtemp->IntegralAndError(totalBGtemp->FindBin(50),totalBGtemp->FindBin(60-1),QCDee50_60Error);
    QCDee60_70=totalBGtemp->IntegralAndError(totalBGtemp->FindBin(60),totalBGtemp->FindBin(70-1),QCDee60_70Error);
    QCDee70_80=totalBGtemp->IntegralAndError(totalBGtemp->FindBin(70),totalBGtemp->FindBin(80-1),QCDee70_80Error);
    QCDee80_100=totalBGtemp->IntegralAndError(totalBGtemp->FindBin(80),totalBGtemp->FindBin(100-1),QCDee80_100Error);
    //make ff, gf combination.  Need to clone instead of get from file because getting points both to same histo
    //TH1F* ffMet_reweightJet_binned_comb = (TH1F*)fin.Get("ffMet_reweightJet_binned");ffMet_reweightJet_binned_comb->Sumw2();
    //TH1F* ffMet_reweightJet_binned_comb = (TH1F*)ffMet_reweightJet_binned->Clone();ffMet_reweightJet_binned_comb->Sumw2();
    //TH1F* gfMet_reweightJet_binned_comb = (TH1F*)gfMet_reweightJet_binned->Clone();gfMet_reweightJet_binned_comb->Sumw2();
    Met_reweightJet_binned_ffSyst = (TH1F*)ffMet_reweightJet_binned->Clone();
    //gfMet_reweightJet_binned->Add(fgMet_reweightJet_binned);


    float gammafakeTOffScale = ffMet_reweightJet_binned->Integral()/gammafakeMet_reweightJet_binned->Integral();
    gammafakeMet_reweightJet_binned->Scale(gammafakeTOffScale);
    gammafakeMet_reweightJet_binned->Scale(gammafakePercent);

    float gfTOffScale = ffMet_reweightJet_binned->Integral()/gfMet_reweightJet_binned->Integral();
    gfMet_reweightJet_binned->Scale(gfTOffScale);
    gfMet_reweightJet_binned->Scale(gfPercent);

    float fgTOffScale = ffMet_reweightJet_binned->Integral()/fgMet_reweightJet_binned->Integral();
    fgMet_reweightJet_binned->Scale(fgTOffScale);
    fgMet_reweightJet_binned->Scale(fgPercent);

    float eeTOffScale = ffMet_reweightJet_binned->Integral()/eeMetNew_ForLimitsStatOnly->Integral();
    eeMetNew_ForLimitsStatOnly->Scale(eeTOffScale);
    eeMetNew_ForLimitsStatOnly->Scale(eePercent);

    ffMet_reweightJet_binned->Scale(ffPercent);

    ffMet_reweightJet_binned->Add(eeMetNew_ForLimitsStatOnly);
    //ffMet_reweightJet_binned->Add(gammafakeMet_reweightJet_binned);
    ffMet_reweightJet_binned->Add(gfMet_reweightJet_binned);
    ffMet_reweightJet_binned->Add(fgMet_reweightJet_binned);
    ffMet_ForLimits = (TH1F*)ffMet_reweightJet_binned->Clone();
    ffMet_ForLimitsStatOnly = (TH1F*)ffMet_reweightJet_binned->Clone();
    for(int i=0;i<sizeof(StdDevff)/sizeof(StdDevff[0]);i++){StdDevff[i]=StdDevgammafakeffcomb[i];}
  }
  
  GetAndSetFFerrors(ggMetData,egMetForff,Met_reweightJet_binned_ffSyst/*eeMetMinusSideBand/*,eeOverBinWidthForFF*/,ffMet_reweightJet_binned,ffMetNew,ffMet_ForLimits,ffMet_ForLimitsStatOnly,FakeRate,FakeRateErr,FakeRateErrForLimits,StdDevff,normErrff,reweightErrff,statErrff,diffFromeeErrorff,diffFromfgErrorff,0,ffScale,ffScaleErr);
  
  fLimits_nojet.cd();
  ffMet_ForLimits->Write("met_ff_nojet");
  fin.cd();

  //ffMetNew = (TH1F*)ffMet_reweightJet_binned->Rebin(NmetBins,"ffMetNew",xbins);
  totalBG->Add(egMetNew,ffMetNew,1,1);
  //totalBG->Add(eeOverBinWidthForFF,1);

  ff0_20=ffMetNew->IntegralAndError(0,4,ff0_20Error);
  QCDff0_20=totalBG->IntegralAndError(0,4,QCDff0_20Error);
  ff30_50=ffMetNew->IntegralAndError(ffMetNew->FindBin(30),ffMetNew->FindBin(50-1),ff30_50Error);
  QCDff30_50=totalBG->IntegralAndError(totalBG->FindBin(30),totalBG->FindBin(50-1),QCDff30_50Error);
  ff50up=ffMetNew->IntegralAndError(ffMetNew->FindBin(50),-1,ff50upError);
  QCDff50up=totalBG->IntegralAndError(totalBG->FindBin(50),-1,QCDff50upError);
  ff100up=ffMetNew->IntegralAndError(ffMetNew->FindBin(100),-1,ff100upError);
  QCDff100up=totalBG->IntegralAndError(totalBG->FindBin(100),-1,QCDff100upError);
  ff50_60=ffMetNew->IntegralAndError(ffMetNew->FindBin(50),ffMetNew->FindBin(60-1),ff50_60Error);
  QCDff50_60=totalBG->IntegralAndError(totalBG->FindBin(50),totalBG->FindBin(60-1),QCDff50_60Error);
  ff60_70=ffMetNew->IntegralAndError(ffMetNew->FindBin(60),ffMetNew->FindBin(70-1),ff60_70Error);
  QCDff60_70=totalBG->IntegralAndError(totalBG->FindBin(60),totalBG->FindBin(70-1),QCDff60_70Error);
  ff70_80=ffMetNew->IntegralAndError(ffMetNew->FindBin(70),ffMetNew->FindBin(80-1),ff70_80Error);
  QCDff70_80=totalBG->IntegralAndError(totalBG->FindBin(70),totalBG->FindBin(80-1),QCDff70_80Error);
  ff80_100=ffMetNew->IntegralAndError(ffMetNew->FindBin(80),ffMetNew->FindBin(100-1),ff80_100Error);
  QCDff80_100=totalBG->IntegralAndError(totalBG->FindBin(80),totalBG->FindBin(100-1),QCDff80_100Error);
   
  TH1F *ffOverBinWidth = new TH1F("ffOverBinWidth","",NmetBins,xbins);ffOverBinWidth->Sumw2();
  for(int i=1;i<NmetBins+1;i++){
    float ff = ffMetNew->GetBinContent(i)/ffMetNew->GetBinWidth(i);
    float ee = eeOverBinWidthForFF->GetBinContent(i)/eeOverBinWidthForFF->GetBinWidth(i);
    float eeE = eeOverBinWidthForFF->GetBinError(i)/eeOverBinWidthForFF->GetBinWidth(i);
    float bg = totalBG->GetBinContent(i)/totalBG->GetBinWidth(i);
    float bgE = totalBG->GetBinError(i)/totalBG->GetBinWidth(i);
    ffOverBinWidth->SetBinContent(i,ff);
    eeOverBinWidthForFF->SetBinContent(i,ee);
    eeOverBinWidthForFF->SetBinError(i,eeE);
    totalBG->SetBinContent(i,bg);
    totalBG->SetBinError(i,bgE);
  }


  eeOverBinWidthForFF->SetMaximum(Max);
  eeOverBinWidthForFF->GetYaxis()->SetTitle("Number of Events / GeV");
  eeOverBinWidthForFF->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  eeOverBinWidthForFF->SetTitle("");
  eeOverBinWidthForFF->SetLineColor(kRed);
  eeOverBinWidthForFF->SetFillColor(kRed);
  eeOverBinWidthForFF->SetStats(0);

  //eeOverBinWidthForFF->SetMarkerSize(0);
  ffOverBinWidth->SetLineColor(kGray+1);
  ffOverBinWidth->SetFillColor(kGray+1);
  //ffOverBinWidth->SetFillStyle(3004);
  ffOverBinWidth->SetMarkerSize(0);
  ffOverBinWidth->SetStats(0);
  ffOverBinWidth->SetTitle("");
  THStack *ffStack = new THStack("ffStack",";;Number of Events / GeV");
  ffStack->SetMaximum(Max);
  ffStack->SetMinimum(Min);
  //ffStack->Add(eeOverBinWidthForFF);
  ffStack->Add(egOverBinWidth);
  ffStack->Add(ffOverBinWidth);
  /*  egOverBinWidth->Draw("B");
      ffOverBinWidth->Draw("BSAME");
      ffOverBinWidth->Draw("SAME");
      egOverBinWidth->Draw("SAME");*/
  ffStack->Draw("histo");
  //totalBG->Draw("E2SAME");
  ggOverBinWidth->Draw("PESAME");
  TLegend *legFF = new TLegend(.5,.45,.8,.65,"","brNDC");
  legFF->AddEntry(ggOverBinWidth,"#gamma#gamma Candidate Sample","lpf");
  legFF->AddEntry(totalBG,"QCD + Electroweak Error","f");
  legFF->AddEntry(ffOverBinWidth,"QCD","f");
  legFF->AddEntry(egOverBinWidth,"Electroweak","f");
  legFF->SetFillColor(kWhite);
  legFF->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_ggANDffBinnedReweighted.png");

  TH1F* ggOverff=(TH1F*)ggOverBinWidth->Clone();ggOverff->SetTitle("ggOverff");ggOverff->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  ggOverff->Divide(ffOverBinWidth);
  c1->SetLogy(0);
  ggOverff->Draw();
  line1->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggOverff.png");
  fout.cd();ggOverff->Write("ggOverff");fin.cd();
  c1->SetLogy(1);
  fout.cd();ffOverBinWidth->Write("ffMetBinnedReweightedAndScaled");fin.cd();
   
  totalBGff = (TH1F*)totalBG->Clone();
  TH1F* ggOverffQCD=(TH1F*)ggOverBinWidth->Clone();ggOverffQCD->SetTitle("ggOverffQCD");ggOverffQCD->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  TH1F* ggOverffQCDerr=(TH1F*)totalBG->Clone();
  ggOverffQCDerr->SetFillColor(kRed);
  //ggOverffQCDerr->SetFillColor(kBlack-1);
  for(int i=1;i<ggOverffQCD->GetNbinsX()+1;i++){
    float Value = ggOverffQCD->GetBinContent(i);float StatErr = ggOverffQCD->GetBinError(i);
    Value/=totalBG->GetBinContent(i);StatErr/=totalBG->GetBinContent(i);
    ggOverffQCD->SetBinContent(i,Value);ggOverffQCD->SetBinError(i,StatErr);
    float SystErr = totalBG->GetBinError(i)/totalBG->GetBinContent(i);
    ggOverffQCDerr->SetBinContent(i,1.);ggOverffQCDerr->SetBinError(i,SystErr);
  }
  c1->SetLogy(0);
  ggOverffQCD->GetYaxis()->SetRangeUser(0.7,1.17);
  ggOverffQCD->Draw("PE");
  ggOverffQCDerr->Draw("E2same");
  line1->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggOverffQCD.png");
  fout.cd();ggOverffQCD->Write("ggOverffQCD");ggOverffQCDerr->Write("ggOverffQCDerr");fin.cd();
  c1->SetLogy(1);


  c2->cd();
  p1->cd();
  p1->SetLogy(1);
  //p1->SetLogx(1);
  /*  egOverBinWidth->Draw("B");
      ffOverBinWidth->Draw("BSAME");
      ffOverBinWidth->Draw("SAME");
      egOverBinWidth->Draw("BSAME");
      egOverBinWidth->Draw("SAME");*/
  ffStack->Draw("histo");
  totalBG->Draw("E2SAME");
  ggOverBinWidth->Draw("PESAME");
  ggSig2012Rebin->Draw("histo SAME");
  ggSig2012Rebin_2->Draw("histo SAME");
  TLegend *legFFwr = new TLegend(.5,.35,.83,.61,"","brNDC");
  legFFwr->AddEntry(ggOverBinWidth,"#gamma#gamma Candidate Sample","lpf");
  legFFwr->AddEntry(totalBG,"QCD + Electroweak Error","f");
  legFFwr->AddEntry(ffOverBinWidth,"QCD","f");
  legFFwr->AddEntry(egOverBinWidth,"Electroweak","f");
  legFFwr->AddEntry(ggSig2012Rebin_2,"GGM #gamma#gamma (1100_720_375)","l");
  legFFwr->AddEntry(ggSig2012Rebin,"GGM #gamma#gamma (1400_1720_375)","l");
  legFFwr->SetFillColor(kWhite);
  legFFwr->Draw("SAME");
  TPaveText *TextFFwr;
  TextFFwr = new TPaveText(.5,.62,.8,.83,"NDC");
  TextFFwr->AddText("CMS Preliminary");
  TextFFwr->AddText("");
  TextFFwr->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  TextFFwr->AddText("No Jet Requirement");
  TextFFwr->SetFillColor(0);
  TextFFwr->SetBorderSize(0);
  TextFFwr->Draw();
  p2->cd();
  //p2->SetLogx(1);
  ggOverffQCDerr->SetTitle("");
  ggOverffQCDerr->GetYaxis()->SetTitle("Data/Prediction");
  ggOverffQCDerr->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  ggOverffQCDerr->GetYaxis()->SetTitleOffset(0.48);
  ggOverffQCDerr->GetYaxis()->SetTitleSize(0.15);
  ggOverffQCDerr->GetXaxis()->SetTitleSize(0.2);
  ggOverffQCDerr->GetYaxis()->SetLabelSize(0.12);
  ggOverffQCDerr->GetXaxis()->SetLabelSize(0.15);
  //ggOverffQCD->GetYaxis()->SetTitleSize(3);
  ggOverffQCDerr->GetYaxis()->SetRangeUser(ratMin,ratMax);
  ggOverffQCDerr->Draw("E2");
  ggOverffQCD->Draw("SAME");
  line1->Draw("SAME");
  p3->cd();
  met->Draw();
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDffBinnedReweighted_WithRatio.png");
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDffBinnedReweighted_WithRatio.pdf");
  c1->cd();

  c2->cd();
  p1->cd();
  ffStack->GetXaxis()->SetRangeUser(0,49);
  ffStack->Draw("histo");
  totalBG->Draw("E2SAME");
  ggOverBinWidth->Draw("PESAME");
  ggSig2012Rebin->Draw("histo SAME");
  ggSig2012Rebin_2->Draw("histo SAME");
  //legFFwr->Draw("SAME");
  //TextFFwr->Draw();
  p2->cd();
  ggOverffQCDerr->GetXaxis()->SetRangeUser(0,49);
  ggOverffQCDerr->GetYaxis()->SetRangeUser(0.8,1.2);
  ggOverffQCDerr->Draw("E2");
  ggOverffQCD->Draw("SAME");
  line1->DrawLine(0,1,50,1);
  p3->cd();
  met->Draw();
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDffBinnedReweighted_WithRatio_Zoom0-50.png");
  //c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDffBinnedReweighted_WithRatio_Zoom0-50.pdf");
  c1->cd();

  //---------------Now ff Jet_Req-----------------
    

  float StdDevff_JetReq[HistBins];//={0.};
  float meanff_JetReq[HistBins];//={0.};
    
  for(int i=0;i<HistBins;i++){
    StdDevff_JetReq[i]=0.;meanff_JetReq[i]=0.;
  }

  TString titleff_JetReq="Toys/ff_JetReqMet_toy_1";
  hist_JetReq=(TH1F*)fin.Get(titleff_JetReq);hist_JetReq->Sumw2();
    
  for(int bin=1;bin<HistBins+1;bin++){
    for(int i=1;i<1001;i++){
      titleff_JetReq="Toys/ff_JetReqMet_toy_";titleff_JetReq+=i;
      hist_JetReq=(TH1F*)fin.Get(titleff_JetReq);if(hist_JetReq->GetSumw2()==0)hist_JetReq->Sumw2();
      meanff_JetReq[bin-1]+=hist_JetReq->GetBinContent(bin);
    }
    meanff_JetReq[bin-1]=meanff_JetReq[bin-1]/1000;
  }
    
  for(int i=1;i<1001;i++){
    titleff_JetReq="Toys/ff_JetReqMet_toy_";titleff_JetReq+=i;
    hist_JetReq = (TH1F*)fin.Get(titleff_JetReq);
    for(int bin=1;bin<HistBins+1;bin++){
      StdDevff_JetReq[bin-1]+= (hist_JetReq->GetBinContent(bin)-meanff_JetReq[bin-1])*(hist_JetReq->GetBinContent(bin)-meanff_JetReq[bin-1]);
    }
  }
  for(int bin=1;bin<HistBins+1;bin++){
    StdDevff_JetReq[bin-1]=sqrt(StdDevff_JetReq[bin-1]/1000);
  }
  //gammafake jetreq
  float StdDevgammafake_JetReq[HistBins];//={0.};
  float meangammafake_JetReq[HistBins];//={0.};
    
  for(int i=0;i<HistBins;i++){
    StdDevgammafake_JetReq[i]=0.;meangammafake_JetReq[i]=0.;
  }

  TString titlegammafake_JetReq="Toys/gammafake_JetReqMet_toy_1";
  hist_JetReq=(TH1F*)fin.Get(titlegammafake_JetReq);hist_JetReq->Sumw2();
    
  for(int bin=1;bin<HistBins+1;bin++){
    for(int i=1;i<1001;i++){
      titlegammafake_JetReq="Toys/gammafake_JetReqMet_toy_";titlegammafake_JetReq+=i;
      hist_JetReq=(TH1F*)fin.Get(titlegammafake_JetReq);if(hist_JetReq->GetSumw2()==0)hist_JetReq->Sumw2();
      meangammafake_JetReq[bin-1]+=hist_JetReq->GetBinContent(bin);
    }
    meangammafake_JetReq[bin-1]=meangammafake_JetReq[bin-1]/1000;
  }
    
  for(int i=1;i<1001;i++){
    titlegammafake_JetReq="Toys/gammafake_JetReqMet_toy_";titlegammafake_JetReq+=i;
    hist_JetReq = (TH1F*)fin.Get(titlegammafake_JetReq);
    for(int bin=1;bin<HistBins+1;bin++){
      StdDevgammafake_JetReq[bin-1]+= (hist_JetReq->GetBinContent(bin)-meangammafake_JetReq[bin-1])*(hist_JetReq->GetBinContent(bin)-meangammafake_JetReq[bin-1]);
    }
  }
  for(int bin=1;bin<HistBins+1;bin++){
    StdDevgammafake_JetReq[bin-1]=sqrt(StdDevgammafake_JetReq[bin-1]/1000);
  }
  //end gammafake
  //gammafake ff combination
  float StdDevgammafakeffcomb_JetReq[HistBins];//={0.};
  for(int bin=1;bin<HistBins+1;bin++){
    StdDevgammafakeffcomb_JetReq[bin-1]=sqrt(gammafakePercent*gammafakePercent*StdDevgammafake_JetReq[bin-1]*StdDevgammafake_JetReq[bin-1]+ffPercent*ffPercent*StdDevff_JetReq[bin-1]*StdDevff_JetReq[bin-1]);
    //cout<<"Bin:" <<bin<<endl<<" BinMin:"<<binMin[bin-1]<<endl<< " BinMax:"<<binMax[bin-1]<<endl<<" STDDEV:"<<StdDev[bin-1]<<endl;
  }
  //end comb
  fin.cd();

  TH1F* ffMet_reweightJet_binned_JetReq = (TH1F*)fin.Get("ffMet_reweightJet_binned_JetReq");ffMet_reweightJet_binned_JetReq->Sumw2();
  TH1F* ffMet_reweightJet_binned_JetReqSyst = (TH1F*)ffMet_reweightJet_binned_JetReq->Clone();ffMet_reweightJet_binned_JetReqSyst->Sumw2();
  TH1F* gfMet_reweightJet_binned_JetReq = (TH1F*)fin.Get("gfMet_reweightJet_binned_JetReq");gfMet_reweightJet_binned_JetReq->Sumw2();
  TH1F* gfMet_reweightJet_binned_JetReqSyst = (TH1F*)gfMet_reweightJet_binned_JetReq->Clone();gfMet_reweightJet_binned_JetReqSyst->Sumw2();
  TH1F* fgMet_reweightJet_binned_JetReq = (TH1F*)fin.Get("fgMet_reweightJet_binned_JetReq");fgMet_reweightJet_binned_JetReq->Sumw2();
  TH1F* gammafakeMet_reweightJet_binned_JetReq = (TH1F*)fin.Get("gammafakeMet_reweightJet_binned_JetReq");gammafakeMet_reweightJet_binned_JetReq->Sumw2();
  TH1F* Met_reweightJet_binned_JetReqffSyst = (TH1F*)gfMet_reweightJet_binned_JetReq->Clone();Met_reweightJet_binned_JetReqffSyst->Sumw2();
  ffMet_reweightJet_binnedForErrors_JetReq=(TH1F*)ffMet_reweightJet_binned_JetReq->Rebin(NmetBins,"ffMet_reweightJet_binnedForErrors_JetReq",xbins);
  TH1F* ffMet_ForLimits_JetReq = (TH1F*)ffMet_reweightJet_binned_JetReq->Clone();
  TH1F* ffMet_ForLimitsStatOnly_JetReq = (TH1F*)ffMet_reweightJet_binned_JetReq->Clone();
  TH1F* ffMetNew_JetReq=(TH1F*)ffMet_reweightJet_binned_JetReq->Clone();

  if(doffgfcomb){
    //this for systematic error, grab ff
    GetAndSetFFerrors(ggMet_JetReq,egMetForff_JetReq,Met_reweightJet_binned_JetReqffSyst/*eeMetMinusSideBand_JetReq/*,eeOverBinWidthForFF*/,ffMet_reweightJet_binned_JetReqSyst,ffMetNew_JetReq,ffMet_ForLimits_JetReq,ffMet_ForLimitsStatOnly_JetReq,FakeRate,FakeRateErr,FakeRateErrForLimits,StdDevff_JetReq,normErrff_JetReq,reweightErrff_JetReq,statErrff_JetReq,diffFromeeErrorff_JetReq,diffFromfgErrorff_JetReq,0,ffScale_JetReq,ffScaleErr_JetReq);
    fLimits_1jet.cd();
    ffMet_ForLimits_JetReq->Write("met_ee_1jet");
    fin.cd();
    ee0_20_JetReq=ffMetNew_JetReq->IntegralAndError(0,4,ee0_20Error_JetReq);
    ee30_50_JetReq=ffMetNew_JetReq->IntegralAndError(ffMetNew_JetReq->FindBin(30),ffMetNew_JetReq->FindBin(50-1),ee30_50Error_JetReq);
    ee50up_JetReq=ffMetNew_JetReq->IntegralAndError(ffMetNew_JetReq->FindBin(50),-1,ee50upError_JetReq);
    ee100up_JetReq=ffMetNew_JetReq->IntegralAndError(ffMetNew_JetReq->FindBin(100),-1,ee100upError_JetReq);
    ee50_60_JetReq=ffMetNew_JetReq->IntegralAndError(ffMetNew_JetReq->FindBin(50),ffMetNew_JetReq->FindBin(60-1),ee50_60Error_JetReq);
    ee60_70_JetReq=ffMetNew_JetReq->IntegralAndError(ffMetNew_JetReq->FindBin(60),ffMetNew_JetReq->FindBin(70-1),ee60_70Error_JetReq);
    ee70_80_JetReq=ffMetNew_JetReq->IntegralAndError(ffMetNew_JetReq->FindBin(70),ffMetNew_JetReq->FindBin(80-1),ee70_80Error_JetReq);
    ee80_100_JetReq=ffMetNew_JetReq->IntegralAndError(ffMetNew_JetReq->FindBin(80),ffMetNew_JetReq->FindBin(100-1),ee80_100Error_JetReq);
    TH1F *totalBGtemp_JetReq=new TH1F("totalBGtemp_JetReq","",NmetBins,xbins);totalBGtemp_JetReq->Sumw2();
    totalBGtemp_JetReq->Add(ffMetNew_JetReq,egMetNew_JetReq,1,1);
    QCDee0_20_JetReq=totalBGtemp_JetReq->IntegralAndError(0,4,QCDee0_20Error_JetReq);
    QCDee30_50_JetReq=totalBGtemp_JetReq->IntegralAndError(totalBGtemp_JetReq->FindBin(30),totalBGtemp_JetReq->FindBin(50-1),QCDee30_50Error_JetReq);
    QCDee50up_JetReq=totalBGtemp_JetReq->IntegralAndError(totalBGtemp_JetReq->FindBin(50),-1,QCDee50upError_JetReq);
    QCDee100up_JetReq=totalBGtemp_JetReq->IntegralAndError(totalBGtemp_JetReq->FindBin(100),-1,QCDee100upError_JetReq);
    QCDee50_60_JetReq=totalBGtemp_JetReq->IntegralAndError(totalBGtemp_JetReq->FindBin(50),totalBGtemp_JetReq->FindBin(60-1),QCDee50_60Error_JetReq);
    QCDee60_70_JetReq=totalBGtemp_JetReq->IntegralAndError(totalBGtemp_JetReq->FindBin(60),totalBGtemp_JetReq->FindBin(70-1),QCDee60_70Error_JetReq);
    QCDee70_80_JetReq=totalBGtemp_JetReq->IntegralAndError(totalBGtemp_JetReq->FindBin(70),totalBGtemp_JetReq->FindBin(80-1),QCDee70_80Error_JetReq);
    QCDee80_100_JetReq=totalBGtemp_JetReq->IntegralAndError(totalBGtemp_JetReq->FindBin(80),totalBGtemp_JetReq->FindBin(100-1),QCDee80_100Error_JetReq); 
    //make ff, gf combination.  Need to clone instead of get from file because getting points both to same histo
    //TH1F* ffMet_reweightJet_binned_JetReq_comb = (TH1F*)fin.Get("ffMet_reweightJet_binned_JetReq");ffMet_reweightJet_binned_JetReq_comb->Sumw2();
    //TH1F* ffMet_reweightJet_binned_JetReq_comb = (TH1F*)ffMet_reweightJet_binned_JetReq->Clone();ffMet_reweightJet_binned_JetReq_comb->Sumw2();
    //TH1F* gfMet_reweightJet_binned_JetReq_comb = (TH1F*)gfMet_reweightJet_binned_JetReq->Clone();gfMet_reweightJet_binned_JetReq_comb->Sumw2();
    Met_reweightJet_binned_JetReqffSyst = (TH1F*)ffMet_reweightJet_binned_JetReq->Clone();
    //gfMet_reweightJet_binned_JetReq->Add(fgMet_reweightJet_binned_JetReq);
    float gfTOffScale_JetReq = ffMet_reweightJet_binned_JetReq->Integral()/gammafakeMet_reweightJet_binned_JetReq->Integral();
    gammafakeMet_reweightJet_binned_JetReq->Scale(gfTOffScale_JetReq);
    ffMet_reweightJet_binned_JetReq->Scale(ffPercent);
    gammafakeMet_reweightJet_binned_JetReq->Scale(gammafakePercent);

    float eeTOffScale_JetReq = ffMet_reweightJet_binned_JetReq->Integral()/eeMetNew_ForLimitsStatOnly_JetReq->Integral();
    eeMetNew_ForLimitsStatOnly_JetReq->Scale(eeTOffScale_JetReq);
    eeMetNew_ForLimitsStatOnly_JetReq->Scale(eePercent);

    ffMet_reweightJet_binned_JetReq->Add(eeMetNew_ForLimitsStatOnly_JetReq);
    ffMet_reweightJet_binned_JetReq->Add(gammafakeMet_reweightJet_binned_JetReq);
    ffMet_ForLimits_JetReq = (TH1F*)ffMet_reweightJet_binned_JetReq->Clone();
    ffMet_ForLimitsStatOnly_JetReq = (TH1F*)ffMet_reweightJet_binned_JetReq->Clone();
    for(int i=0;i<sizeof(StdDevff_JetReq)/sizeof(StdDevff_JetReq[0]);i++){StdDevff_JetReq[i]=StdDevgammafakeffcomb_JetReq[i];}
  }
  
  GetAndSetFFerrors(ggMet_JetReq,egMetForff_JetReq,Met_reweightJet_binned_JetReqffSyst/*eeMetMinusSideBand_JetReq/*,eeOverBinWidthForFF*/,ffMet_reweightJet_binned_JetReq,ffMetNew_JetReq,ffMet_ForLimits_JetReq,ffMet_ForLimitsStatOnly_JetReq,FakeRate,FakeRateErr,FakeRateErrForLimits,StdDevff_JetReq,normErrff_JetReq,reweightErrff_JetReq,statErrff_JetReq,diffFromeeErrorff_JetReq,diffFromfgErrorff_JetReq,0,ffScale_JetReq,ffScaleErr_JetReq);
  

  fLimits_1jet.cd();
  ffMet_ForLimits_JetReq->Write("met_ff_1jet");
  fin.cd();

  //TH1F* ffMetNew_JetReq=(TH1F*)ffMet_reweightJet_binned_JetReq->Rebin(NmetBins,"ffMetNew_JetReq",xbins);



  totalBG_JetReq->Add(egMetNew_JetReq,ffMetNew_JetReq,1,1);

  ff0_20_JetReq=ffMetNew_JetReq->IntegralAndError(0,4,ff0_20Error_JetReq);
  QCDff0_20_JetReq=totalBG_JetReq->IntegralAndError(0,4,QCDff0_20Error_JetReq);
  ff30_50_JetReq=ffMetNew_JetReq->IntegralAndError(ffMetNew_JetReq->FindBin(30),ffMetNew_JetReq->FindBin(50-1),ff30_50Error_JetReq);
  QCDff30_50_JetReq=totalBG_JetReq->IntegralAndError(totalBG_JetReq->FindBin(30),totalBG_JetReq->FindBin(50-1),QCDff30_50Error_JetReq);
  ff50up_JetReq=ffMetNew_JetReq->IntegralAndError(ffMetNew_JetReq->FindBin(50),-1,ff50upError_JetReq);
  QCDff50up_JetReq=totalBG_JetReq->IntegralAndError(totalBG_JetReq->FindBin(50),-1,QCDff50upError_JetReq);
  ff100up_JetReq=ffMetNew_JetReq->IntegralAndError(ffMetNew_JetReq->FindBin(100),-1,ff100upError_JetReq);
  QCDff100up_JetReq=totalBG_JetReq->IntegralAndError(totalBG_JetReq->FindBin(100),-1,QCDff100upError_JetReq);
  ff50_60_JetReq=ffMetNew_JetReq->IntegralAndError(ffMetNew_JetReq->FindBin(50),ffMetNew_JetReq->FindBin(60-1),ff50_60Error_JetReq);
  QCDff50_60_JetReq=totalBG_JetReq->IntegralAndError(totalBG_JetReq->FindBin(50),totalBG_JetReq->FindBin(60-1),QCDff50_60Error_JetReq);
  ff60_70_JetReq=ffMetNew_JetReq->IntegralAndError(ffMetNew_JetReq->FindBin(60),ffMetNew_JetReq->FindBin(70-1),ff60_70Error_JetReq);
  QCDff60_70_JetReq=totalBG_JetReq->IntegralAndError(totalBG_JetReq->FindBin(60),totalBG_JetReq->FindBin(70-1),QCDff60_70Error_JetReq);
  ff70_80_JetReq=ffMetNew_JetReq->IntegralAndError(ffMetNew_JetReq->FindBin(70),ffMetNew_JetReq->FindBin(80-1),ff70_80Error_JetReq);
  QCDff70_80_JetReq=totalBG_JetReq->IntegralAndError(totalBG_JetReq->FindBin(70),totalBG_JetReq->FindBin(80-1),QCDff70_80Error_JetReq);
  ff80_100_JetReq=ffMetNew_JetReq->IntegralAndError(ffMetNew_JetReq->FindBin(80),ffMetNew_JetReq->FindBin(100-1),ff80_100Error_JetReq);
  QCDff80_100_JetReq=totalBG_JetReq->IntegralAndError(totalBG_JetReq->FindBin(80),totalBG_JetReq->FindBin(100-1),QCDff80_100Error_JetReq);

  
  TH1F *ffOverBinWidth_JetReq = new TH1F("ffOverBinWidth_JetReq","",NmetBins,xbins);ffOverBinWidth_JetReq->Sumw2();
  for(int i=1;i<NmetBins+1;i++){
    float ff_JetReq = ffMetNew_JetReq->GetBinContent(i)/ffMetNew_JetReq->GetBinWidth(i);
    float bg_JetReq = totalBG_JetReq->GetBinContent(i)/totalBG_JetReq->GetBinWidth(i);
    float bgE_JetReq = totalBG_JetReq->GetBinError(i)/totalBG_JetReq->GetBinWidth(i);
    ffOverBinWidth_JetReq->SetBinContent(i,ff_JetReq);
    totalBG_JetReq->SetBinContent(i,bg_JetReq);
    totalBG_JetReq->SetBinError(i,bgE_JetReq);
  }
  ffOverBinWidth_JetReq->SetLineColor(kGray+1);
  ffOverBinWidth_JetReq->SetFillColor(kGray+1);
  //ffOverBinWidth_JetReq->SetFillStyle(3004);
  ffOverBinWidth_JetReq->SetMarkerSize(0);
  ffOverBinWidth_JetReq->SetStats(0);
  ffOverBinWidth_JetReq->SetTitle("");
  THStack *ffStack_JetReq = new THStack("ffStack_JetReq",";;Number of Events / GeV");
  ffStack_JetReq->SetMaximum(Max);
  ffStack_JetReq->SetMinimum(Min);
  ffStack_JetReq->Add(egOverBinWidth_JetReq);
  ffStack_JetReq->Add(ffOverBinWidth_JetReq);
  /*  egOverBinWidth_JetReq->Draw("B");
      ffOverBinWidth_JetReq->Draw("BSAME");
      ffOverBinWidth_JetReq->Draw("SAME");
      egOverBinWidth_JetReq->Draw("SAME");*/
  ffStack_JetReq->Draw();
  totalBG_JetReq->Draw("E2SAME");
  ggOverBinWidth_JetReq->Draw("PESAME");
  legFF->Draw("SAME");
  Text_JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_ggANDffBinnedReweighted_JetReq.png");
  fout.cd();ffOverBinWidth_JetReq->Write("ffMetBinnedReweightedAndScaled_JetReq");fin.cd();

  TH1F* ggOverff_JetReq=(TH1F*)ggOverBinWidth_JetReq->Clone();ggOverff_JetReq->SetTitle("ggOverff_JetReq");ggOverff_JetReq->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  ggOverff_JetReq->Divide(ffOverBinWidth_JetReq);
  c1->SetLogy(0);
  ggOverff_JetReq->Draw();
  line1->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggOverff_JetReq.png");
  fout.cd();ggOverff_JetReq->Write("ggOverff_JetReq");fin.cd();
  c1->SetLogy(1);
  fout.cd();ffOverBinWidth_JetReq->Write("ffMetBinnedReweightedAndScaled_JetReq");fin.cd();
  /*
    TH1F* ggOverffQCD_JetReq=(TH1F*)ggOverBinWidth_JetReq->Clone();ggOverffQCD_JetReq->SetTitle("ggOverffQCD_JetReq");ggOverffQCD_JetReq->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
    ggOverffQCD_JetReq->Divide(totalBG_JetReq);
    c1->SetLogy(0);
    ggOverffQCD_JetReq->Draw();
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggOverffQCD_JetReq.png");
    fout.cd();ggOverffQCD->Write("ggOverffQCD_JetReq");fin.cd();
    c1->SetLogy(1);
  */
  totalBGff_JetReq = (TH1F*)totalBG_JetReq->Clone();
  TH1F* ggOverffQCD_JetReq=(TH1F*)ggOverBinWidth_JetReq->Clone();ggOverffQCD_JetReq->SetTitle("ggOverffQCD_JetReq");ggOverffQCD_JetReq->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  TH1F* ggOverffQCDerr_JetReq=(TH1F*)totalBG_JetReq->Clone();
  for(int i=1;i<ggOverffQCD_JetReq->GetNbinsX()+1;i++){
    float Value = ggOverffQCD_JetReq->GetBinContent(i);float StatErr = ggOverffQCD_JetReq->GetBinError(i);
    Value/=totalBG_JetReq->GetBinContent(i);StatErr/=totalBG_JetReq->GetBinContent(i);
    ggOverffQCD_JetReq->SetBinContent(i,Value);ggOverffQCD_JetReq->SetBinError(i,StatErr);
    float SystErr = totalBG_JetReq->GetBinError(i)/totalBG_JetReq->GetBinContent(i);
    ggOverffQCDerr_JetReq->SetBinContent(i,1.);ggOverffQCDerr_JetReq->SetBinError(i,SystErr);
  }
  //ggOverffQCD->Divide(totalBG);
  c1->SetLogy(0);
  ggOverffQCDerr_JetReq->SetFillColor(kRed);
  ggOverffQCD_JetReq->GetYaxis()->SetRangeUser(0.7,1.17);
  ggOverffQCD_JetReq->Draw("pe");
  ggOverffQCDerr_JetReq->Draw("E2same");
  line1->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggOverffQCD_JetReq.png");
  fout.cd();ggOverffQCD_JetReq->Write("ggOverffQCD_JetReq");ggOverffQCDerr_JetReq->Write("ggOverffQCDerr_JetReq");fin.cd();
  c1->SetLogy(1);



  c2->cd();
  p1->cd();
  p1->SetLogy(1);
  /*   egOverBinWidth_JetReq->Draw("B");
       ffOverBinWidth_JetReq->Draw("BSAME");
       ffOverBinWidth_JetReq->Draw("SAME");
       egOverBinWidth_JetReq->Draw("BSAME");
       egOverBinWidth_JetReq->Draw("SAME");*/
  ffStack_JetReq->Draw();
  totalBG_JetReq->Draw("E2SAME");
  ggOverBinWidth_JetReq->Draw("PESAME");
  ggSig2012Rebin_JetReq->Draw("histo SAME");
  ggSig2012Rebin_JetReq_2->Draw("histo SAME");
  TLegend *legFFwr_JetReq = new TLegend(.5,.35,.83,.61,"","brNDC");
  legFFwr_JetReq->AddEntry(ggOverBinWidth_JetReq,"#gamma#gamma Candidate Sample","lpf");
  legFFwr_JetReq->AddEntry(totalBG_JetReq,"QCD + Electroweak Error","f");
  legFFwr_JetReq->AddEntry(ffOverBinWidth_JetReq,"QCD","f");
  legFFwr_JetReq->AddEntry(egOverBinWidth_JetReq,"Electroweak","f");
  legFFwr_JetReq->AddEntry(ggSig2012Rebin_JetReq_2,"GGM #gamma#gamma (1100_720_375)","l");
  legFFwr_JetReq->AddEntry(ggSig2012Rebin_JetReq,"GGM #gamma#gamma (1400_1720_375)","l");
  legFFwr_JetReq->SetFillColor(kWhite);
  legFFwr_JetReq->Draw("SAME");
  TPaveText *TextFFwr_JetReq;
  TextFFwr_JetReq = new TPaveText(.5,.62,.8,.83,"NDC");
  TextFFwr_JetReq->AddText("CMS Preliminary");
  TextFFwr_JetReq->AddText("");
  TextFFwr_JetReq->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  TextFFwr_JetReq->AddText(">=1 Jet Requirement");
  TextFFwr_JetReq->SetFillColor(0);
  TextFFwr_JetReq->SetBorderSize(0);
  TextFFwr_JetReq->Draw();
  p2->cd();
  ggOverffQCDerr_JetReq->SetTitle("");
  ggOverffQCDerr_JetReq->GetYaxis()->SetTitle("Data/Prediction");
  ggOverffQCDerr_JetReq->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  ggOverffQCDerr_JetReq->GetYaxis()->SetTitleOffset(0.48);
  ggOverffQCDerr_JetReq->GetYaxis()->SetTitleSize(0.15);
  ggOverffQCDerr_JetReq->GetXaxis()->SetTitleSize(0.2);
  ggOverffQCDerr_JetReq->GetYaxis()->SetLabelSize(0.12);
  ggOverffQCDerr_JetReq->GetXaxis()->SetLabelSize(0.15);
  //ggOverffQCD_JetReq->GetYaxis()->SetTitleSize(3);
  ggOverffQCDerr_JetReq->GetYaxis()->SetRangeUser(ratMin,ratMax);
  ggOverffQCDerr_JetReq->Draw("E2");
  ggOverffQCD_JetReq->Draw("SAME");
  line1->Draw("SAME");
  p3->cd();
  met->Draw();
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDffBinnedReweighted_WithRatio_JetReq.png");
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDffBinnedReweighted_WithRatio_JetReq.pdf");
  c1->cd();

  c2->cd();
  p1->cd();
  ffStack_JetReq->GetXaxis()->SetRangeUser(0,49);
  ffStack_JetReq->Draw("histo");
  totalBG_JetReq->Draw("E2SAME");
  ggOverBinWidth_JetReq->Draw("PESAME");
  ggSig2012Rebin_JetReq->Draw("histo SAME");
  ggSig2012Rebin_JetReq_2->Draw("histo SAME");
  //legFFwr->Draw("SAME");
  //TextFFwr->Draw();
  p2->cd();
  ggOverffQCDerr_JetReq->GetXaxis()->SetRangeUser(0,49);
  ggOverffQCDerr_JetReq->GetYaxis()->SetRangeUser(0.8,1.2);
  ggOverffQCDerr_JetReq->Draw("E2");
  ggOverffQCD_JetReq->Draw("SAME");
  line1->DrawLine(0,1,50,1);
  p3->cd();
  met->Draw();
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDffBinnedReweighted_WithRatio_JetReq_Zoom0-50.png");
  //c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDffBinnedReweighted_WithRatio_Zoom0-50.pdf");
  c1->cd();


  //---------------Now ff 2JetReq-----------------
    

  float StdDevff_2JetReq[HistBins];//={0.};
  float meanff_2JetReq[HistBins];//={0.};
    
  for(int i=0;i<HistBins;i++){
    StdDevff_2JetReq[i]=0.;meanff_2JetReq[i]=0.;
  }

  TString titleff_2JetReq="Toys/ff_2JetReqMet_toy_1";
  hist_2JetReq=(TH1F*)fin.Get(titleff_2JetReq);hist_2JetReq->Sumw2();
    
  for(int bin=1;bin<HistBins+1;bin++){
    for(int i=1;i<1001;i++){
      titleff_2JetReq="Toys/ff_2JetReqMet_toy_";titleff_2JetReq+=i;
      hist_2JetReq=(TH1F*)fin.Get(titleff_2JetReq);if(hist_2JetReq->GetSumw2()==0)hist_2JetReq->Sumw2();
      meanff_2JetReq[bin-1]+=hist_2JetReq->GetBinContent(bin);
    }
    meanff_2JetReq[bin-1]=meanff_2JetReq[bin-1]/1000;
  }
    
  for(int i=1;i<1001;i++){
    titleff_2JetReq="Toys/ff_2JetReqMet_toy_";titleff_2JetReq+=i;
    hist_2JetReq = (TH1F*)fin.Get(titleff_2JetReq);
    for(int bin=1;bin<HistBins+1;bin++){
      StdDevff_2JetReq[bin-1]+= (hist_2JetReq->GetBinContent(bin)-meanff_2JetReq[bin-1])*(hist_2JetReq->GetBinContent(bin)-meanff_2JetReq[bin-1]);
    }
  }
  for(int bin=1;bin<HistBins+1;bin++){
    StdDevff_2JetReq[bin-1]=sqrt(StdDevff_2JetReq[bin-1]/1000);
  }
 //gammafake 2jetreq
  float StdDevgammafake_2JetReq[HistBins];//={0.};
  float meangammafake_2JetReq[HistBins];//={0.};
    
  for(int i=0;i<HistBins;i++){
    StdDevgammafake_2JetReq[i]=0.;meangammafake_2JetReq[i]=0.;
  }

  TString titlegammafake_2JetReq="Toys/gammafake_2JetReqMet_toy_1";
  hist_2JetReq=(TH1F*)fin.Get(titlegammafake_2JetReq);hist_2JetReq->Sumw2();
    
  for(int bin=1;bin<HistBins+1;bin++){
    for(int i=1;i<1001;i++){
      titlegammafake_2JetReq="Toys/gammafake_2JetReqMet_toy_";titlegammafake_2JetReq+=i;
      hist_2JetReq=(TH1F*)fin.Get(titlegammafake_2JetReq);if(hist_2JetReq->GetSumw2()==0)hist_2JetReq->Sumw2();
      meangammafake_2JetReq[bin-1]+=hist_2JetReq->GetBinContent(bin);
    }
    meangammafake_2JetReq[bin-1]=meangammafake_2JetReq[bin-1]/1000;
  }
    
  for(int i=1;i<1001;i++){
    titlegammafake_2JetReq="Toys/gammafake_2JetReqMet_toy_";titlegammafake_2JetReq+=i;
    hist_2JetReq = (TH1F*)fin.Get(titlegammafake_2JetReq);
    for(int bin=1;bin<HistBins+1;bin++){
      StdDevgammafake_2JetReq[bin-1]+= (hist_2JetReq->GetBinContent(bin)-meangammafake_2JetReq[bin-1])*(hist_2JetReq->GetBinContent(bin)-meangammafake_2JetReq[bin-1]);
    }
  }
  for(int bin=1;bin<HistBins+1;bin++){
    StdDevgammafake_2JetReq[bin-1]=sqrt(StdDevgammafake_2JetReq[bin-1]/1000);
  }
  //end gammafake
  //gammafake ff 2jetreq combination
  float StdDevgammafakeffcomb_2JetReq[HistBins];//={0.};
  for(int bin=1;bin<HistBins+1;bin++){
    StdDevgammafakeffcomb_2JetReq[bin-1]=sqrt(gammafakePercent*gammafakePercent*StdDevgammafake_2JetReq[bin-1]*StdDevgammafake_2JetReq[bin-1]+ffPercent*ffPercent*StdDevff_2JetReq[bin-1]*StdDevff_2JetReq[bin-1]);
    //cout<<"Bin:" <<bin<<endl<<" BinMin:"<<binMin[bin-1]<<endl<< " BinMax:"<<binMax[bin-1]<<endl<<" STDDEV:"<<StdDev[bin-1]<<endl;
  }
  //end comb
  fin.cd();

  TH1F* ffMet_reweightJet_binned_2JetReq = (TH1F*)fin.Get("ffMet_reweightJet_binned_2JetReq");ffMet_reweightJet_binned_2JetReq->Sumw2();
  TH1F* ffMet_reweightJet_binned_2JetReqSyst = (TH1F*)ffMet_reweightJet_binned_2JetReq->Clone();ffMet_reweightJet_binned_2JetReqSyst->Sumw2();
  TH1F* gfMet_reweightJet_binned_2JetReq = (TH1F*)fin.Get("gfMet_reweightJet_binned_2JetReq");gfMet_reweightJet_binned_2JetReq->Sumw2();
  TH1F* gfMet_reweightJet_binned_2JetReqSyst = (TH1F*)gfMet_reweightJet_binned_2JetReq->Clone();gfMet_reweightJet_binned_2JetReqSyst->Sumw2();
  TH1F* fgMet_reweightJet_binned_2JetReq = (TH1F*)fin.Get("fgMet_reweightJet_binned_2JetReq");fgMet_reweightJet_binned_2JetReq->Sumw2();
  TH1F* gammafakeMet_reweightJet_binned_2JetReq = (TH1F*)fin.Get("gammafakeMet_reweightJet_binned_2JetReq");gammafakeMet_reweightJet_binned_2JetReq->Sumw2();
  TH1F* Met_reweightJet_binned_2JetReqffSyst = (TH1F*)gfMet_reweightJet_binned_2JetReq->Clone();Met_reweightJet_binned_2JetReqffSyst->Sumw2();
  ffMet_reweightJet_binnedForErrors_2JetReq=(TH1F*)ffMet_reweightJet_binned_2JetReq->Rebin(NmetBins,"ffMet_reweightJet_binnedForErrors_2JetReq",xbins);
  TH1F* ffMet_ForLimits_2JetReq = (TH1F*)ffMet_reweightJet_binned_2JetReq->Clone();
  TH1F* ffMet_ForLimitsStatOnly_2JetReq = (TH1F*)ffMet_reweightJet_binned_2JetReq->Clone();
  TH1F* ffMetNew_2JetReq=(TH1F*)ffMet_reweightJet_binned_2JetReq->Clone();

  if(doffgfcomb){
    //this for systematic error, grab ff
    GetAndSetFFerrors(ggMet_2JetReq,egMetForff_2JetReq,Met_reweightJet_binned_2JetReqffSyst/*eeMetMinusSideBand_2JetReq/*,eeOverBinWidthForFF*/,ffMet_reweightJet_binned_2JetReqSyst,ffMetNew_2JetReq,ffMet_ForLimits_2JetReq,ffMet_ForLimitsStatOnly_2JetReq,FakeRate,FakeRateErr,FakeRateErrForLimits,StdDevff_2JetReq,normErrff_2JetReq,reweightErrff_2JetReq,statErrff_2JetReq,diffFromeeErrorff_2JetReq,diffFromfgErrorff_2JetReq,0,ffScale_2JetReq,ffScaleErr_2JetReq);
    fLimits_2jet.cd();
    ffMet_ForLimits_2JetReq->Write("met_ee_2jet");
    fin.cd();
    ee0_20_2JetReq=ffMetNew_2JetReq->IntegralAndError(0,4,ee0_20Error_2JetReq);
    ee30_50_2JetReq=ffMetNew_2JetReq->IntegralAndError(ffMetNew_2JetReq->FindBin(30),ffMetNew_2JetReq->FindBin(50-1),ee30_50Error_2JetReq);
    ee50up_2JetReq=ffMetNew_2JetReq->IntegralAndError(ffMetNew_2JetReq->FindBin(50),-1,ee50upError_2JetReq);
    ee100up_2JetReq=ffMetNew_2JetReq->IntegralAndError(ffMetNew_2JetReq->FindBin(100),-1,ee100upError_2JetReq);
    ee50_60_2JetReq=ffMetNew_2JetReq->IntegralAndError(ffMetNew_2JetReq->FindBin(50),ffMetNew_2JetReq->FindBin(60-1),ee50_60Error_2JetReq);
    ee60_70_2JetReq=ffMetNew_2JetReq->IntegralAndError(ffMetNew_2JetReq->FindBin(60),ffMetNew_2JetReq->FindBin(70-1),ee60_70Error_2JetReq);
    ee70_80_2JetReq=ffMetNew_2JetReq->IntegralAndError(ffMetNew_2JetReq->FindBin(70),ffMetNew_2JetReq->FindBin(80-1),ee70_80Error_2JetReq);
    ee80_100_2JetReq=ffMetNew_2JetReq->IntegralAndError(ffMetNew_2JetReq->FindBin(80),ffMetNew_2JetReq->FindBin(100-1),ee80_100Error_2JetReq);
    TH1F *totalBGtemp_2JetReq=new TH1F("totalBGtemp_2JetReq","",NmetBins,xbins);totalBGtemp_2JetReq->Sumw2();
    totalBGtemp_2JetReq->Add(ffMetNew_2JetReq,egMetNew_2JetReq,1,1);
    QCDee0_20_2JetReq=totalBGtemp_2JetReq->IntegralAndError(0,4,QCDee0_20Error_2JetReq);
    QCDee30_50_2JetReq=totalBGtemp_2JetReq->IntegralAndError(totalBGtemp_2JetReq->FindBin(30),totalBGtemp_2JetReq->FindBin(50-1),QCDee30_50Error_2JetReq);
    QCDee50up_2JetReq=totalBGtemp_2JetReq->IntegralAndError(totalBGtemp_2JetReq->FindBin(50),-1,QCDee50upError_2JetReq);
    QCDee100up_2JetReq=totalBGtemp_2JetReq->IntegralAndError(totalBGtemp_2JetReq->FindBin(100),-1,QCDee100upError_2JetReq);
    QCDee50_60_2JetReq=totalBGtemp_2JetReq->IntegralAndError(totalBGtemp_2JetReq->FindBin(50),totalBGtemp_2JetReq->FindBin(60-1),QCDee50_60Error_2JetReq);
    QCDee60_70_2JetReq=totalBGtemp_2JetReq->IntegralAndError(totalBGtemp_2JetReq->FindBin(60),totalBGtemp_2JetReq->FindBin(70-1),QCDee60_70Error_2JetReq);
    QCDee70_80_2JetReq=totalBGtemp_2JetReq->IntegralAndError(totalBGtemp_2JetReq->FindBin(70),totalBGtemp_2JetReq->FindBin(80-1),QCDee70_80Error_2JetReq);
    QCDee80_100_2JetReq=totalBGtemp_2JetReq->IntegralAndError(totalBGtemp_2JetReq->FindBin(80),totalBGtemp_2JetReq->FindBin(100-1),QCDee80_100Error_2JetReq); 
    //make ff, gf combination.  Need to clone instead of get from file because getting points both to same histo
    //TH1F* ffMet_reweightJet_binned_2JetReq_comb = (TH1F*)fin.Get("ffMet_reweightJet_binned_2JetReq");ffMet_reweightJet_binned_2JetReq_comb->Sumw2();
    //TH1F* ffMet_reweightJet_binned_2JetReq_comb = (TH1F*)ffMet_reweightJet_binned_2JetReq->Clone();ffMet_reweightJet_binned_2JetReq_comb->Sumw2();
    //TH1F* gfMet_reweightJet_binned_2JetReq_comb = (TH1F*)gfMet_reweightJet_binned_2JetReq->Clone();gfMet_reweightJet_binned_2JetReq_comb->Sumw2();
    Met_reweightJet_binned_2JetReqffSyst = (TH1F*)ffMet_reweightJet_binned_2JetReq->Clone();
    //gfMet_reweightJet_binned_2JetReq->Add(fgMet_reweightJet_binned_2JetReq);
    float gfTOffScale_2JetReq = ffMet_reweightJet_binned_2JetReq->Integral()/gammafakeMet_reweightJet_binned_2JetReq->Integral();
    gammafakeMet_reweightJet_binned_2JetReq->Scale(gfTOffScale_2JetReq);
    ffMet_reweightJet_binned_2JetReq->Scale(ffPercent);
    gammafakeMet_reweightJet_binned_2JetReq->Scale(gammafakePercent);
   
    float eeTOffScale_2JetReq = ffMet_reweightJet_binned_2JetReq->Integral()/eeMetNew_ForLimitsStatOnly_2JetReq->Integral();
    eeMetNew_ForLimitsStatOnly_2JetReq->Scale(eeTOffScale_2JetReq);
    eeMetNew_ForLimitsStatOnly_2JetReq->Scale(eePercent);

    ffMet_reweightJet_binned_2JetReq->Add(eeMetNew_ForLimitsStatOnly_2JetReq); 
    ffMet_reweightJet_binned_2JetReq->Add(gammafakeMet_reweightJet_binned_2JetReq);
    ffMet_ForLimits_2JetReq = (TH1F*)ffMet_reweightJet_binned_2JetReq->Clone();
    ffMet_ForLimitsStatOnly_2JetReq = (TH1F*)ffMet_reweightJet_binned_2JetReq->Clone();
    for(int i=0;i<sizeof(StdDevff_2JetReq)/sizeof(StdDevff_2JetReq[0]);i++){StdDevff_2JetReq[i]=StdDevgammafakeffcomb_2JetReq[i];}
  }
  
  GetAndSetFFerrors(ggMet_2JetReq,egMetForff_2JetReq,Met_reweightJet_binned_2JetReqffSyst/*eeMetMinusSideBand_2JetReq/*,eeOverBinWidthForFF*/,ffMet_reweightJet_binned_2JetReq,ffMetNew_2JetReq,ffMet_ForLimits_2JetReq,ffMet_ForLimitsStatOnly_2JetReq,FakeRate,FakeRateErr,FakeRateErrForLimits,StdDevff_2JetReq,normErrff_2JetReq,reweightErrff_2JetReq,statErrff_2JetReq,diffFromeeErrorff_2JetReq,diffFromfgErrorff_2JetReq,0,ffScale_2JetReq,ffScaleErr_2JetReq);
  

  fLimits_2jet.cd();
  ffMet_ForLimits_2JetReq->Write("met_ff_2jet");
  fin.cd();

  //TH1F* ffMetNew_2JetReq=(TH1F*)ffMet_reweightJet_binned_2JetReq->Rebin(NmetBins,"ffMetNew_2JetReq",xbins);



  totalBG_2JetReq->Add(egMetNew_2JetReq,ffMetNew_2JetReq,1,1);



  ff0_20_2JetReq=ffMetNew_2JetReq->IntegralAndError(0,4,ff0_20Error_2JetReq);
  QCDff0_20_2JetReq=totalBG_2JetReq->IntegralAndError(0,4,QCDff0_20Error_2JetReq);
  ff30_50_2JetReq=ffMetNew_2JetReq->IntegralAndError(ffMetNew_2JetReq->FindBin(30),ffMetNew_2JetReq->FindBin(50-1),ff30_50Error_2JetReq);
  QCDff30_50_2JetReq=totalBG_2JetReq->IntegralAndError(totalBG_2JetReq->FindBin(30),totalBG_2JetReq->FindBin(50-1),QCDff30_50Error_2JetReq);
  ff50up_2JetReq=ffMetNew_2JetReq->IntegralAndError(ffMetNew_2JetReq->FindBin(50),-1,ff50upError_2JetReq);
  QCDff50up_2JetReq=totalBG_2JetReq->IntegralAndError(totalBG_2JetReq->FindBin(50),-1,QCDff50upError_2JetReq);
  ff100up_2JetReq=ffMetNew_2JetReq->IntegralAndError(ffMetNew_2JetReq->FindBin(100),-1,ff100upError_2JetReq);
  QCDff100up_2JetReq=totalBG_2JetReq->IntegralAndError(totalBG_2JetReq->FindBin(100),-1,QCDff100upError_2JetReq);
  ff50_60_2JetReq=ffMetNew_2JetReq->IntegralAndError(ffMetNew_2JetReq->FindBin(50),ffMetNew_2JetReq->FindBin(60-1),ff50_60Error_2JetReq);
  QCDff50_60_2JetReq=totalBG_2JetReq->IntegralAndError(totalBG_2JetReq->FindBin(50),totalBG_2JetReq->FindBin(60-1),QCDff50_60Error_2JetReq);
  ff60_70_2JetReq=ffMetNew_2JetReq->IntegralAndError(ffMetNew_2JetReq->FindBin(60),ffMetNew_2JetReq->FindBin(70-1),ff60_70Error_2JetReq);
  QCDff60_70_2JetReq=totalBG_2JetReq->IntegralAndError(totalBG_2JetReq->FindBin(60),totalBG_2JetReq->FindBin(70-1),QCDff60_70Error_2JetReq);
  ff70_80_2JetReq=ffMetNew_2JetReq->IntegralAndError(ffMetNew_2JetReq->FindBin(70),ffMetNew_2JetReq->FindBin(80-1),ff70_80Error_2JetReq);
  QCDff70_80_2JetReq=totalBG_2JetReq->IntegralAndError(totalBG_2JetReq->FindBin(70),totalBG_2JetReq->FindBin(80-1),QCDff70_80Error_2JetReq);
  ff80_100_2JetReq=ffMetNew_2JetReq->IntegralAndError(ffMetNew_2JetReq->FindBin(80),ffMetNew_2JetReq->FindBin(100-1),ff80_100Error_2JetReq);
  QCDff80_100_2JetReq=totalBG_2JetReq->IntegralAndError(totalBG_2JetReq->FindBin(80),totalBG_2JetReq->FindBin(100-1),QCDff80_100Error_2JetReq);

  
  TH1F *ffOverBinWidth_2JetReq = new TH1F("ffOverBinWidth_2JetReq","",NmetBins,xbins);ffOverBinWidth_2JetReq->Sumw2();
  for(int i=1;i<NmetBins+1;i++){
    float ff_2JetReq = ffMetNew_2JetReq->GetBinContent(i)/ffMetNew_2JetReq->GetBinWidth(i);
    float bg_2JetReq = totalBG_2JetReq->GetBinContent(i)/totalBG_2JetReq->GetBinWidth(i);
    float bgE_2JetReq = totalBG_2JetReq->GetBinError(i)/totalBG_2JetReq->GetBinWidth(i);
    ffOverBinWidth_2JetReq->SetBinContent(i,ff_2JetReq);
    totalBG_2JetReq->SetBinContent(i,bg_2JetReq);
    totalBG_2JetReq->SetBinError(i,bgE_2JetReq);
  }
  ffOverBinWidth_2JetReq->SetLineColor(kGray+1);
  ffOverBinWidth_2JetReq->SetFillColor(kGray+1);
  //ffOverBinWidth_2JetReq->SetFillStyle(3004);
  ffOverBinWidth_2JetReq->SetMarkerSize(0);
  ffOverBinWidth_2JetReq->SetStats(0);
  ffOverBinWidth_2JetReq->SetTitle("");
  THStack *ffStack_2JetReq = new THStack("ffStack_2JetReq",";;Number of Events / GeV");
  ffStack_2JetReq->SetMaximum(Max);
  ffStack_2JetReq->SetMinimum(Min);
  ffStack_2JetReq->Add(egOverBinWidth_2JetReq);
  ffStack_2JetReq->Add(ffOverBinWidth_2JetReq);
  /*  egOverBinWidth_2JetReq->Draw("B");
      ffOverBinWidth_2JetReq->Draw("BSAME");
      ffOverBinWidth_2JetReq->Draw("SAME");
      egOverBinWidth_2JetReq->Draw("SAME");*/
  ffStack_2JetReq->Draw();
  totalBG_2JetReq->Draw("E2SAME");
  ggOverBinWidth_2JetReq->Draw("PESAME");
  legFF->Draw("SAME");
  Text_2JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_ggANDffBinnedReweighted_2JetReq.png");
  fout.cd();ffOverBinWidth_2JetReq->Write("ffMetBinnedReweightedAndScaled_2JetReq");fin.cd();

  TH1F* ggOverff_2JetReq=(TH1F*)ggOverBinWidth_2JetReq->Clone();ggOverff_2JetReq->SetTitle("ggOverff_2JetReq");ggOverff_2JetReq->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  ggOverff_2JetReq->Divide(ffOverBinWidth_2JetReq);
  c1->SetLogy(0);
  ggOverff_2JetReq->Draw();
  line1->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggOverff_2JetReq.png");
  fout.cd();ggOverff_2JetReq->Write("ggOverff_2JetReq");fin.cd();
  c1->SetLogy(1);
  fout.cd();ffOverBinWidth_2JetReq->Write("ffMetBinnedReweightedAndScaled_2JetReq");fin.cd();
  /*
    TH1F* ggOverffQCD_2JetReq=(TH1F*)ggOverBinWidth_2JetReq->Clone();ggOverffQCD_2JetReq->SetTitle("ggOverffQCD_2JetReq");ggOverffQCD_2JetReq->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
    ggOverffQCD_2JetReq->Divide(totalBG_2JetReq);
    c1->SetLogy(0);
    ggOverffQCD_2JetReq->Draw();
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggOverffQCD_2JetReq.png");
    fout.cd();ggOverffQCD->Write("ggOverffQCD_2JetReq");fin.cd();
    c1->SetLogy(1);
  */
  totalBGff_2JetReq = (TH1F*)totalBG_2JetReq->Clone();
  TH1F* ggOverffQCD_2JetReq=(TH1F*)ggOverBinWidth_2JetReq->Clone();ggOverffQCD_2JetReq->SetTitle("ggOverffQCD_2JetReq");ggOverffQCD_2JetReq->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  TH1F* ggOverffQCDerr_2JetReq=(TH1F*)totalBG_2JetReq->Clone();
  for(int i=1;i<ggOverffQCD_2JetReq->GetNbinsX()+1;i++){
    float Value = ggOverffQCD_2JetReq->GetBinContent(i);float StatErr = ggOverffQCD_2JetReq->GetBinError(i);
    Value/=totalBG_2JetReq->GetBinContent(i);StatErr/=totalBG_2JetReq->GetBinContent(i);
    ggOverffQCD_2JetReq->SetBinContent(i,Value);ggOverffQCD_2JetReq->SetBinError(i,StatErr);
    float SystErr = totalBG_2JetReq->GetBinError(i)/totalBG_2JetReq->GetBinContent(i);
    ggOverffQCDerr_2JetReq->SetBinContent(i,1.);ggOverffQCDerr_2JetReq->SetBinError(i,SystErr);
  }
  //ggOverffQCD->Divide(totalBG);
  c1->SetLogy(0);
  ggOverffQCDerr_2JetReq->SetFillColor(kRed);
  ggOverffQCD_2JetReq->GetYaxis()->SetRangeUser(0.7,1.17);
  ggOverffQCD_2JetReq->Draw("pe");
  ggOverffQCDerr_2JetReq->Draw("E2same");
  line1->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggOverffQCD_2JetReq.png");
  fout.cd();ggOverffQCD_2JetReq->Write("ggOverffQCD_2JetReq");ggOverffQCDerr_2JetReq->Write("ggOverffQCDerr_2JetReq");fin.cd();
  c1->SetLogy(1);



  c2->cd();
  p1->cd();
  p1->SetLogy(1);
  /*   egOverBinWidth_2JetReq->Draw("B");
       ffOverBinWidth_2JetReq->Draw("BSAME");
       ffOverBinWidth_2JetReq->Draw("SAME");
       egOverBinWidth_2JetReq->Draw("BSAME");
       egOverBinWidth_2JetReq->Draw("SAME");*/
  ffStack_2JetReq->Draw();
  totalBG_2JetReq->Draw("E2SAME");
  ggOverBinWidth_2JetReq->Draw("PESAME");
  ggSig2012Rebin_2JetReq->Draw("histo SAME");
  ggSig2012Rebin_2JetReq_2->Draw("histo SAME");
  TLegend *legFFwr_2JetReq = new TLegend(.5,.35,.83,.61,"","brNDC");
  legFFwr_2JetReq->AddEntry(ggOverBinWidth_2JetReq,"#gamma#gamma Candidate Sample","lpf");
  legFFwr_2JetReq->AddEntry(totalBG_2JetReq,"QCD + Electroweak Error","f");
  legFFwr_2JetReq->AddEntry(ffOverBinWidth_2JetReq,"QCD","f");
  legFFwr_2JetReq->AddEntry(egOverBinWidth_2JetReq,"Electroweak","f");
  legFFwr_2JetReq->AddEntry(ggSig2012Rebin_2JetReq_2,"GGM #gamma#gamma (1100_720_375)","l");
  legFFwr_2JetReq->AddEntry(ggSig2012Rebin_2JetReq,"GGM #gamma#gamma (1400_1720_375)","l");
  legFFwr_2JetReq->SetFillColor(kWhite);
  legFFwr_2JetReq->Draw("SAME");
  TPaveText *TextFFwr_2JetReq;
  TextFFwr_2JetReq = new TPaveText(.5,.62,.8,.83,"NDC");
  TextFFwr_2JetReq->AddText("CMS Preliminary");
  TextFFwr_2JetReq->AddText("");
  TextFFwr_2JetReq->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  TextFFwr_2JetReq->AddText(">=1 Jet Requirement");
  TextFFwr_2JetReq->SetFillColor(0);
  TextFFwr_2JetReq->SetBorderSize(0);
  TextFFwr_2JetReq->Draw();
  p2->cd();
  ggOverffQCDerr_2JetReq->SetTitle("");
  ggOverffQCDerr_2JetReq->GetYaxis()->SetTitle("Data/Prediction");
  ggOverffQCDerr_2JetReq->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  ggOverffQCDerr_2JetReq->GetYaxis()->SetTitleOffset(0.48);
  ggOverffQCDerr_2JetReq->GetYaxis()->SetTitleSize(0.15);
  ggOverffQCDerr_2JetReq->GetXaxis()->SetTitleSize(0.2);
  ggOverffQCDerr_2JetReq->GetYaxis()->SetLabelSize(0.12);
  ggOverffQCDerr_2JetReq->GetXaxis()->SetLabelSize(0.15);
  //ggOverffQCD_2JetReq->GetYaxis()->SetTitleSize(3);
  ggOverffQCDerr_2JetReq->GetYaxis()->SetRangeUser(ratMin,ratMax);
  ggOverffQCDerr_2JetReq->Draw("E2");
  ggOverffQCD_2JetReq->Draw("SAME");
  line1->Draw("SAME");
  p3->cd();
  met->Draw();
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDffBinnedReweighted_WithRatio_2JetReq.png");
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDffBinnedReweighted_WithRatio_2JetReq.pdf");
  c1->cd();

  c2->cd();
  p1->cd();
  ffStack_2JetReq->GetXaxis()->SetRangeUser(0,49);
  ffStack_2JetReq->Draw("histo");
  totalBG_2JetReq->Draw("E2SAME");
  ggOverBinWidth_2JetReq->Draw("PESAME");
  ggSig2012Rebin_2JetReq->Draw("histo SAME");
  ggSig2012Rebin_2JetReq_2->Draw("histo SAME");
  //legFFwr->Draw("SAME");
  //TextFFwr->Draw();
  p2->cd();
  ggOverffQCDerr_2JetReq->GetXaxis()->SetRangeUser(0,49);
  ggOverffQCDerr_2JetReq->GetYaxis()->SetRangeUser(0.8,1.2);
  ggOverffQCDerr_2JetReq->Draw("E2");
  ggOverffQCD_2JetReq->Draw("SAME");
  line1->DrawLine(0,1,50,1);
  p3->cd();
  met->Draw();
  c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDffBinnedReweighted_WithRatio_2JetReq_Zoom0-50.png");
  //c2->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_Met_ggANDffBinnedReweighted_WithRatio_Zoom0-50.pdf");
  c1->cd();

  //end 2JetReq

  c1->SetLogy(0);

  TH1F* eeOverff = (TH1F*)eeMetNew_ForLimitsStatOnly->Rebin(NmetBins,"eeOverff",xbins);
  TH1F* ff = (TH1F*)ffMet_ForLimitsStatOnly->Rebin(NmetBins,"ff",xbins);
  eeOverff->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  eeOverff->GetYaxis()->SetTitle("ee sample / ff sample");
  TH1F* eeOverffErr = (TH1F*)eeMetNew_ForLimits->Rebin(NmetBins,"eeOverffErr",xbins);//eeOverffErr->Sumw2();
  TH1F* ffErr = (TH1F*)ffMet_ForLimits->Rebin(NmetBins,"ffErr",xbins);//ffErr->Sumw2();
  eeOverff->Divide(ff);
  //TH1F* eeOverff = (TH1F*)eeOverffOld->Rebin(NmetBins,"eeOverff",xbins);
  eeOverffErr->Divide(ffErr);
  for(int i=1;i<eeOverff->GetNbinsX()+1;i++){
    //float temp=eeOverff->GetBinContent(i)/ff->GetBinContent(i);
    //float StatErr=temp*sqrt((1)/(eeOverff->GetBinContent(i))+(1)/(ff->GetBinContent(i)));
    //eeOverff->SetBinContent(i,temp);
    //eeOverff->SetBinError(i,StatErr);
    //float SystErr=eeOverff->GetBinContent(i)*sqrt((reweightErree[i-1]*reweightErree[i-1])/(eeMetMinusSideBandForErrors->GetBinContent(i)*eeMetMinusSideBandForErrors->GetBinContent(i))+(reweightErrff[i-1]*reweightErrff[i-1])/(ffMet_reweightJet_binnedForErrors->GetBinContent(i)*ffMet_reweightJet_binnedForErrors->GetBinContent(i)));
    eeOverffErr->SetBinContent(i,1.);//eeOverffErr->SetBinError(i,SystErr);
  }
  eeOverffErr->GetYaxis()->SetRangeUser(0,2);
  eeOverffErr->SetTitle("");
  eeOverffErr->SetFillStyle(3001);
  eeOverffErr->SetFillColor(kGray+2);
  eeOverffErr->SetMarkerSize(0);
  eeOverff->SetMarkerSize(0.75);
  eeOverffErr->Draw("E2");
  eeOverff->Draw("PEsame");
  Text_noJet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_eeOverff.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_eeOverff.pdf");

  TH1F* eeOverff_JetReq = (TH1F*)eeMetNew_ForLimitsStatOnly_JetReq->Rebin(NmetBins,"eeOverff_JetReq",xbins);
  TH1F* ff_JetReq = (TH1F*)ffMet_ForLimitsStatOnly_JetReq->Rebin(NmetBins,"ff_JetReq",xbins);
  eeOverff_JetReq->GetXaxis()->SetTitle("Missing E_{T} (GeV)");
  eeOverff_JetReq->GetYaxis()->SetTitle("ee sample / ff sample");
  TH1F* eeOverffErr_JetReq = (TH1F*)eeMetNew_ForLimits_JetReq->Rebin(NmetBins,"eeOverffErr_JetReq",xbins);
  TH1F* ffErr_JetReq = (TH1F*)ffMet_ForLimits_JetReq->Rebin(NmetBins,"ffErr_JetReq",xbins);
  eeOverff_JetReq->Divide(ff_JetReq);
  //TH1F* eeOverff_JetReq = (TH1F*)eeOverffOld_JetReq->Rebin(NmetBins,"eeOverff_JetReq",xbins);
  eeOverffErr_JetReq->Divide(ffErr_JetReq);
  for(int i=1;i<eeOverff_JetReq->GetNbinsX()+1;i++){
    //float temp_JetReq=eeOverff_JetReq->GetBinContent(i)/ff_JetReq->GetBinContent(i);
    //float StatErr_JetReq=temp_JetReq*sqrt((1)/(eeOverff->GetBinContent(i))+(1)/(ff->GetBinContent(i)));
    //eeOverff->SetBinContent(i,temp);
    //eeOverff->SetBinError(i,StatErr);
    //float SystErr=eeOverff->GetBinContent(i)*sqrt((reweightErree[i-1]*reweightErree[i-1])/(eeMetMinusSideBandForErrors->GetBinContent(i)*eeMetMinusSideBandForErrors->GetBinContent(i))+(reweightErrff[i-1]*reweightErrff[i-1])/(ffMet_reweightJet_binnedForErrors->GetBinContent(i)*ffMet_reweightJet_binnedForErrors->GetBinContent(i)));
    eeOverffErr_JetReq->SetBinContent(i,1.);//eeOverffErr->SetBinError(i,SystErr);
  }
  eeOverffErr_JetReq->GetYaxis()->SetRangeUser(0,2);
  eeOverffErr_JetReq->SetTitle("");
  eeOverffErr_JetReq->SetFillStyle(3001);
  eeOverffErr_JetReq->SetFillColor(kGray+2);
  eeOverffErr_JetReq->SetMarkerSize(0);
  eeOverff_JetReq->SetMarkerSize(0.75);
  eeOverffErr_JetReq->Draw("E2");
  eeOverff_JetReq->Draw("PEsame");
  Text_Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_eeOverff_JetReq.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure18_eeOverff_JetReq.pdf");

  totalBGee->SetLineWidth(2.2);totalBGee->SetMarkerSize(.6);
  totalBGff->SetLineWidth(2.2);totalBGff->SetMarkerSize(.6);
  totalBGee->SetLineColor(kRed);totalBGee->SetMarkerColor(kRed);
  totalBGff->SetLineColor(kBlue);totalBGff->SetMarkerColor(kBlue);
  totalBGff->GetXaxis()->SetRangeUser(0,39);
  totalBGff->Draw("PE");
  ggOverBinWidth->Draw("PESAMES");
  totalBGee->Draw("PESAMES");
  TLegend *lowMetLeg = new TLegend(.64,.6,.84,.84,"","brNDC");
  lowMetLeg->AddEntry(ggOverBinWidth,"gg","lp");
  lowMetLeg->AddEntry(totalBGee,"tot bg - ee","lp");
  lowMetLeg->AddEntry(totalBGff,"tot bg - ff","lp");
  lowMetLeg->SetFillStyle(0);
  lowMetLeg->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/gg_ee_ff_lowMET.png");

  totalBGee_JetReq->SetLineWidth(2.2);totalBGee_JetReq->SetMarkerSize(.6);
  totalBGff_JetReq->SetLineWidth(2.2);totalBGff_JetReq->SetMarkerSize(.6);
  totalBGee_JetReq->SetLineColor(kRed);totalBGee_JetReq->SetMarkerColor(kRed);
  totalBGff_JetReq->SetLineColor(kBlue);totalBGff_JetReq->SetMarkerColor(kBlue);
  totalBGff_JetReq->GetXaxis()->SetRangeUser(0,39);
  totalBGff_JetReq->Draw("PE");
  ggOverBinWidth_JetReq->Draw("PESAMES");
  totalBGee_JetReq->Draw("PESAMES");
  lowMetLeg->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/gg_ee_ff_lowMET_JetReq.png");


  TH1F* ggMetComp = (TH1F*)fin.Get("ggMet");    
  TH1F* ggType1MetComp = (TH1F*)fin.Get("ggType1CorrMet");    
  TH1F* ggMetCompNew = (TH1F*)ggMetComp->Rebin(NmetBins,"ggMetCompNew",xbins);
  TH1F* ggType1MetCompNew = (TH1F*)ggType1MetComp->Rebin(NmetBins,"ggMetCompNew",xbins);
  for(int i=0;i<ggMetCompNew->GetNbinsX()+1;i++){
    float temp = ggMetCompNew->GetBinContent(i)/ggMetCompNew->GetBinWidth(i);
    float tempErr = ggMetCompNew->GetBinError(i)/ggMetCompNew->GetBinWidth(i);
    ggMetCompNew->SetBinContent(i,temp);
    ggMetCompNew->SetBinError(i,tempErr);
    temp = ggType1MetCompNew->GetBinContent(i)/ggType1MetCompNew->GetBinWidth(i);
    tempErr = ggType1MetCompNew->GetBinError(i)/ggType1MetCompNew->GetBinWidth(i);
    ggType1MetCompNew->SetBinContent(i,temp);
    ggType1MetCompNew->SetBinError(i,tempErr);
  }
  ggType1MetCompNew->SetLineColor(kRed);ggType1MetCompNew->SetMarkerColor(kRed);
  ggMetCompNew->GetXaxis()->SetRangeUser(0,39);
  ggMetCompNew->Draw("PE");
  ggType1MetCompNew->Draw("PESAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggMET_ggType1CorrMet_Comp_METlt40.png");

  c1->SetLogy(1);

  ggMetCompNew->GetXaxis()->SetRangeUser(0,400);
  ggMetCompNew->Draw("PE");
  ggType1MetCompNew->Draw("PESAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggMET_ggType1CorrMet_Comp_Full_LOG.png");

  c1->SetLogy(0);

  TH1F* ffMetComp = (TH1F*)fin.Get("ffMet");    
  TH1F* ffType1MetComp = (TH1F*)fin.Get("ffType1CorrMet");    
  TH1F* ffMetCompNew = (TH1F*)ffMetComp->Rebin(NmetBins,"ffMetCompNew",xbins);
  TH1F* ffType1MetCompNew = (TH1F*)ffType1MetComp->Rebin(NmetBins,"ffMetCompNew",xbins);
  for(int i=0;i<ffMetCompNew->GetNbinsX()+1;i++){
    float temp = ffMetCompNew->GetBinContent(i)/ffMetCompNew->GetBinWidth(i);
    float tempErr = ffMetCompNew->GetBinError(i)/ffMetCompNew->GetBinWidth(i);
    ffMetCompNew->SetBinContent(i,temp);
    ffMetCompNew->SetBinError(i,tempErr);
    temp = ffType1MetCompNew->GetBinContent(i)/ffType1MetCompNew->GetBinWidth(i);
    tempErr = ffType1MetCompNew->GetBinError(i)/ffType1MetCompNew->GetBinWidth(i);
    ffType1MetCompNew->SetBinContent(i,temp);
    ffType1MetCompNew->SetBinError(i,tempErr);
  }
  ffType1MetCompNew->SetLineColor(kRed);ffType1MetCompNew->SetMarkerColor(kRed);
  ffMetCompNew->GetXaxis()->SetRangeUser(0,39);
  ffMetCompNew->Draw("PE");
  ffType1MetCompNew->Draw("PESAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ffMET_ffType1CorrMet_Comp_METlt40.png");

  c1->SetLogy(1);

  ffMetCompNew->GetXaxis()->SetRangeUser(0,400);
  ffMetCompNew->GetYaxis()->SetRangeUser(1e-4,6e3);
  ffMetCompNew->Draw("PE");
  ffType1MetCompNew->Draw("PESAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ffMET_ffType1CorrMet_Comp_Full_LOG.png");




  //-----------Now do combined----------------------
  TH1F *QCDcombined = new TH1F("QCDcombined","",NmetBins,xbins);
  QCDcombined->Sumw2();
  //ffMet_reweightJet_binned->Rebin(NmetBins,"ffMet_reweightJet_binned_new",xbins);
  //QCDcombined->Add(eeMetMinusSideBand,ffMet_reweightJet_binned_new,1,1);
  QCDcombined->Add(eeMetMinusSideBand2,ffMetNew,1,1);
  QCDcombined->Scale(0.5,"");
  QCDcombined->Scale((ggMetNew->Integral(0,4)-egMetNew->Integral(0,4))/QCDcombined->Integral(0,4) ,"");
  comb100up = QCDcombined->Integral(12,13);
  totalBG->Add(egMetNew,QCDcombined,1,1);
  TH1F *QCDOverBinWidth = new TH1F("QCDOverBinWidth","",NmetBins,xbins);QCDOverBinWidth->Sumw2();
  for(int i=1;i<NmetBins+1;i++){
    float QCD = QCDcombined->GetBinContent(i)/QCDcombined->GetBinWidth(i);
    float bg = totalBG->GetBinContent(i)/totalBG->GetBinWidth(i);
    float bgE = totalBG->GetBinError(i)/totalBG->GetBinWidth(i);
    QCDOverBinWidth->SetBinContent(i,QCD);
    totalBG->SetBinContent(i,bg);
    totalBG->SetBinError(i,bgE);
  }
  QCDOverBinWidth->SetFillColor(kBlue);
  QCDOverBinWidth->SetFillStyle(3004);
  QCDOverBinWidth->SetMarkerSize(0);
  QCDOverBinWidth->SetStats(0);
  QCDOverBinWidth->SetTitle("");
  egOverBinWidth->Draw("B");
  QCDOverBinWidth->Draw("BSAME");
  QCDOverBinWidth->Draw("SAME");
  egOverBinWidth->Draw("SAME");
  totalBG->Draw("E2SAME");
  ggOverBinWidth->Draw("PESAME");
  TLegend *legQCD = new TLegend(.5,.45,.8,.65,"","brNDC");
  legQCD->AddEntry(ggOverBinWidth,"#gamma#gamma Candidate Sample","lpf");
  legQCD->AddEntry(totalBG,"QCD + Electroweak","f");
  legQCD->AddEntry(QCDcombined,"QCD (ee and ff Combined)","f");
  legQCD->AddEntry(egOverBinWidth,"Electroweak","f");
  legQCD->SetFillColor(kWhite);
  legQCD->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_ggANDeeffCombinedBinnedReweighted.png");
  //------------Now combined _JetReq-----------------
  TH1F *QCDcombined_JetReq = new TH1F("QCDcombined_JetReq","",NmetBins,xbins);
  QCDcombined_JetReq->Sumw2();
  //ffMet_reweightJet_binned_JetReq->Rebin(NmetBins,"ffMet_reweightJet_binned_JetReq_new",xbins);
  //QCDcombined_JetReq->Add(eeMetMinusSideBand_JetReq,ffMet_reweightJet_binned_JetReq_new,1,1);
  QCDcombined_JetReq->Add(eeMetMinusSideBand2_JetReq,ffMetNew_JetReq,1,1);
  QCDcombined_JetReq->Scale(0.5,"");
  QCDcombined_JetReq->Scale((ggMetNew_JetReq->Integral(0,4)-egMetNew_JetReq->Integral(0,4))/QCDcombined_JetReq->Integral(0,4) ,"");
  comb100up_JetReq = QCDcombined_JetReq->Integral(12,13);
  totalBG_JetReq->Add(egMetNew_JetReq,QCDcombined_JetReq,1,1);
  TH1F *QCDOverBinWidth_JetReq = new TH1F("QCDOverBinWidth_JetReq","",NmetBins,xbins);QCDOverBinWidth_JetReq->Sumw2();
  for(int i=1;i<NmetBins+1;i++){
    float QCD_JetReq = QCDcombined_JetReq->GetBinContent(i)/QCDcombined_JetReq->GetBinWidth(i);
    float bg = totalBG_JetReq->GetBinContent(i)/totalBG_JetReq->GetBinWidth(i);
    float bgE = totalBG_JetReq->GetBinError(i)/totalBG_JetReq->GetBinWidth(i);
    QCDOverBinWidth_JetReq->SetBinContent(i,QCD_JetReq);
    totalBG_JetReq->SetBinContent(i,bg);
    totalBG_JetReq->SetBinError(i,bgE);
  }
  QCDOverBinWidth_JetReq->SetFillColor(kBlue);
  QCDOverBinWidth_JetReq->SetFillStyle(3004);
  QCDOverBinWidth_JetReq->SetMarkerSize(0);
  QCDOverBinWidth_JetReq->SetStats(0);
  QCDOverBinWidth_JetReq->SetTitle("");
  egOverBinWidth_JetReq->Draw("B");
  QCDOverBinWidth_JetReq->Draw("BSAME");
  QCDOverBinWidth_JetReq->Draw("SAME");
  egOverBinWidth_JetReq->Draw("SAME");
  totalBG_JetReq->Draw("E2SAME");
  ggOverBinWidth_JetReq->Draw("PESAME");
  TLegend *legQCD_JetReq = new TLegend(.5,.45,.8,.65,"","brNDC");
  legQCD_JetReq->AddEntry(ggOverBinWidth,"#gamma#gamma Candidate Sample","lpf");
  legQCD_JetReq->AddEntry(totalBG,"QCD + Electroweak","f");
  legQCD_JetReq->AddEntry(QCDcombined,"QCD (ee and ff Combined)","f");
  legQCD_JetReq->AddEntry(egOverBinWidth,"Electroweak","f");
  legQCD_JetReq->SetFillColor(kWhite);
  legQCD_JetReq->Draw("SAME");
  Text_JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_ggANDeeffCombinedBinnedReweighted_JetReq.png");
  //----------
    
  //now Razor MR and R2 variables
  //first ee
  /*
    float egScale = ggMR->Integral(1,60)/egMR->Integral(1,60);
    egMR->Scale(egScale,"");
    egMR->Scale(0.0145,"");
    TH1F* eeMRMinusSideBand = (TH1F*)eeMR->Clone();//new TH1F("eeMRMinusSideBand","",ggMR->GetNbinsX(),0.,1500.);
    eeMRMinusSideBand->Sumw2();
    eeMRMinusSideBand->Add(eeMR,eeSideBandMR,1,-2);
    eeMRMinusSideBand->Scale( (ggMR->Integral(1,60)-egMR->Integral(1,60))/eeMRMinusSideBand->Integral(1,60) , "");
    TH1F* totBG = (TH1F*)egMR->Clone();
    totBG->Add(egMR, eeMRMinusSideBand,1,1);
    ggMR->Rebin(2);totBG->Rebin(2);egMR->Rebin(2);eeMRMinusSideBand->Rebin(2);
    ggMR->SetMarkerColor(kBlack);
    ggMR->SetLineColor(kBlack);
    ggMR->SetLineWidth(2);
    ggMR->SetMarkerSize(.8);
    ggMR->SetStats(0);
    float Max = ggMR->GetMaximum()+20000;
    egMR->SetMaximum(Max);
    egMR->GetYaxis()->SetTitle("Number of Events");
    egMR->GetXaxis()->SetTitle("Razor M_{R}");
    egMR->SetTitle("");
    egMR->SetFillColor(kCyan);
    egMR->SetFillStyle(3005);
    egMR->SetStats(0);
    eeMRMinusSideBand->SetFillColor(kAzure);
    eeMRMinusSideBand->SetFillStyle(3004);
    eeMRMinusSideBand->SetMarkerSize(0);
    eeMRMinusSideBand->SetStats(0);
    eeMRMinusSideBand->SetTitle("");
    totBG->SetMarkerSize(0);
    totBG->SetFillColor(kRed);
    totBG->SetFillStyle(3003);
    totBG->SetStats(0);
    egMR->Draw("B");
    eeMRMinusSideBand->Draw("BSAME");
    eeMRMinusSideBand->Draw("SAME");
    egMR->Draw("BSAME");
    egMR->Draw("SAME");
    totBG->Draw("E2SAME");
    ggMR->Draw("PESAME");
    // ggSig720New->Draw("SAME");
    //ggSig800New->Draw("SAME");
    TLegend *legMR = new TLegend(.5,.45,.8,.65,"","brNDC");
    legMR->AddEntry(ggMR,"#gamma#gamma Candidate Sample","lpf");
    legMR->AddEntry(totBG,"QCD + Electroweak","f");
    legMR->AddEntry(eeMRMinusSideBand,"QCD (e^{+}e^{-} sample)","f");
    legMR->AddEntry(egMR,"Electroweak","f");
    legMR->SetFillColor(kWhite);
    legMR->Draw("SAME");
    Text = new TPaveText(.5,.66,.8,.87,"NDC");
    Text->AddText("CMS Preliminary");
    Text->AddText("");
    Text->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
    Text->AddText("No Jet Requirement");
    Text->SetFillColor(0);
    Text->SetBorderSize(0);
    Text->Draw();
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MR_ggANDeeNoReweight.png");
    fout.cd();eeMRMinusSideBand->Write("eeMRMinusSideBandBinnedReweightedAndScaled");ggMR->Write("ggMRReweightedAndScaled");fin.cd();

    egScale = ggR2->Integral(1,10)/egR2->Integral(1,10);
    egR2->Scale(egScale,"");
    egR2->Scale(0.0145,"");
    TH1F* eeR2MinusSideBand = (TH1F*)eeR2->Clone();
    eeR2MinusSideBand->Sumw2();
    eeR2MinusSideBand->Add(eeR2,eeSideBandR2,1,-2);
    eeR2MinusSideBand->Scale( (ggR2->Integral(1,10)-egR2->Integral(1,10))/eeR2MinusSideBand->Integral(1,10) , "");
    TH1F* totBGeR2 = (TH1F*)ggR2->Clone();
    totBGeR2->Add(egR2, eeR2MinusSideBand,1,1);
    ggR2->Rebin(2);totBGeR2->Rebin(2);egR2->Rebin(2);eeR2MinusSideBand->Rebin(2);
    ggR2->SetMarkerColor(kBlack);
    ggR2->SetLineColor(kBlack);
    ggR2->SetLineWidth(2);
    ggR2->SetMarkerSize(.8);
    ggR2->SetStats(0);
    float Max = ggR2->GetMaximum()+2000;
    egR2->SetMaximum(Max);
    egR2->GetYaxis()->SetTitle("Number of Events");
    egR2->GetXaxis()->SetTitle("Razor R^{2}");
    egR2->SetTitle("");
    egR2->SetFillColor(kCyan);
    egR2->SetFillStyle(3005);
    egR2->SetStats(0);
    eeR2->SetFillColor(kAzure);
    eeR2->SetFillStyle(3004);
    eeR2->SetMarkerSize(0);
    eeR2->SetStats(0);
    eeR2->SetTitle("");
    totBGeR2->SetMarkerSize(0);
    totBGeR2->SetFillColor(kRed);
    totBGeR2->SetFillStyle(3003);
    totBGeR2->SetStats(0);
    egR2->Draw("B");
    eeR2->Draw("BSAME");
    eeR2->Draw("SAME");
    egR2->Draw("BSAME");
    egR2->Draw("SAME");
    totBGeR2->Draw("E2SAME");
    ggR2->Draw("PESAME");
    // ggSig720New->Draw("SAME");
    //ggSig800New->Draw("SAME");
    TLegend *legR2 = new TLegend(.5,.45,.8,.65,"","brNDC");
    legR2->AddEntry(ggR2,"#gamma#gamma Candidate Sample","lpf");
    legR2->AddEntry(totBGeR2,"QCD + Electroweak","f");
    legR2->AddEntry(eeR2,"QCD (e^{+}e^{-} sample)","f");
    legR2->AddEntry(egR2,"Electroweak","f");
    legR2->SetFillColor(kWhite);
    legR2->Draw("SAME");
    Text = new TPaveText(.5,.66,.8,.87,"NDC");
    Text->AddText("CMS Preliminary");
    Text->AddText("");
    Text->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
    Text->AddText("No Jet Requirement");
    Text->SetFillColor(0);
    Text->SetBorderSize(0);
    Text->Draw();
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/R2_ggANDeeNoReweight.png");
    fout.cd();eeR2->Write("eeR2BinnedReweightedAndScaled");ggR2->Write("ggR2ReweightedAndScaled");fin.cd();
  */
  //Now ff  
  //c1->SetLogy(0);
  //TH1F* ggMR = fSig.Get("MR");ggMR->SetStats(0);

  TH2F* ggR2vsMR = (TH2F*)fin.Get("ggR2vsMR");
  TH2F* egR2vsMR = (TH2F*)fin.Get("egR2vsMR");
  TH2F* eeR2vsMR = (TH2F*)fin.Get("eeR2vsMR_reweight_binned");
  TH2F* eeSideBandR2vsMR = (TH2F*)fin.Get("eeSideBandR2vsMR_reweight_binned");
  TH2F* ffR2vsMR = (TH2F*)fin.Get("ffR2vsMR_reweight_binned");
  TH1F* ggMr_R2point8cut = new TH1F("ggMr_R2point8cut","",ggR2vsMR->GetNbinsX(),ggR2vsMR->GetXaxis()->GetBinLowEdge(1),ggR2vsMR->GetXaxis()->GetBinLowEdge(ggR2vsMR->GetNbinsX()+1));
  //cout <<ggR2vsMR->GetXaxis()->GetBinLowEdge(1)<<endl<<ggR2vsMR->GetXaxis()->GetBinLowEdge(ggR2vsMR->GetNbinsX()+1)<<endl;
  TH1F* egMr_R2point8cut = new TH1F("egMr_R2point8cut","",ggR2vsMR->GetNbinsX(),ggR2vsMR->GetXaxis()->GetBinLowEdge(1),ggR2vsMR->GetXaxis()->GetBinLowEdge(ggR2vsMR->GetNbinsX()+1));
  TH1F* eeMr_R2point8cut = new TH1F("eeMr_R2point8cut","",ggR2vsMR->GetNbinsX(),ggR2vsMR->GetXaxis()->GetBinLowEdge(1),ggR2vsMR->GetXaxis()->GetBinLowEdge(ggR2vsMR->GetNbinsX()+1));
  TH1F* eeSideBandMr_R2point8cut = new TH1F("eeSideBandMr_R2point8cut","",ggR2vsMR->GetNbinsX(),ggR2vsMR->GetXaxis()->GetBinLowEdge(1),ggR2vsMR->GetXaxis()->GetBinLowEdge(ggR2vsMR->GetNbinsX()+1));
  TH1F* ffMr_R2point8cut = new TH1F("ffMr_R2point8cut","",ggR2vsMR->GetNbinsX(),ggR2vsMR->GetXaxis()->GetBinLowEdge(1),ggR2vsMR->GetXaxis()->GetBinLowEdge(ggR2vsMR->GetNbinsX()+1));
  //ggMr_R2point8cut->Sumw2();    


  //int j=0;
  for(int i=0;i<ggR2vsMR->GetNbinsX();i++){
    //cout<<"bin:"<<i<<" LowEdge:"<<ggR2vsMR->GetYaxis()->GetBinLowEdge(i)<<endl;
    //cout<<ggR2vsMR->GetNbinsY()<<endl;
    ggMr_R2point8cut->SetBinContent(i,ggR2vsMR->Integral(i,i,0,ggR2vsMR->GetNbinsY()+1));
    egMr_R2point8cut->SetBinContent(i,egR2vsMR->Integral(i,i,0,ggR2vsMR->GetNbinsY()+1));
    eeMr_R2point8cut->SetBinContent(i,eeR2vsMR->Integral(i,i,0,ggR2vsMR->GetNbinsY()+1));
    eeSideBandMr_R2point8cut->SetBinContent(i,eeSideBandR2vsMR->Integral(i,i,4,ggR2vsMR->GetNbinsY()+1));
    ffMr_R2point8cut->SetBinContent(i,ffR2vsMR->Integral(i,i,0,ggR2vsMR->GetNbinsY()+1));
    //cout<<"bin:"<<i<<" binlowedge: "<<ggMr_R2point8cut->GetBinLowEdge(i)<<" Integral:"<<ggR2vsMR->Integral(i,i,0,ggR2vsMR->GetNbinsY()+1)<<endl;
  }
  float egScale = ggMr_R2point8cut->Integral(1,60)/egMr_R2point8cut->Integral(1,60);
  egMr_R2point8cut->Scale(egScale,"");
  egMr_R2point8cut->Scale(0.0145,"");
  ffMr_R2point8cut->Scale( (ggMr_R2point8cut->Integral(1,60)-egMr_R2point8cut->Integral(1,60))/ffMr_R2point8cut->Integral(1,60) , "");
  TH1F* totBGMR = (TH1F*)ffMr_R2point8cut->Clone();
  totBGMR->Add(egMr_R2point8cut,ffMr_R2point8cut,1,1);
  /*  for(int i=0;i<totBGMR->GetNbinsX();i++){
      cout<<"bin:"<<i<<" binlowedge: "<<totBGMR->GetBinLowEdge(i)<<endl;
      }*/
  Double_t MrBins[20]={
    0,  //1
    20,
    40,
    60,
    80,
    100,
    120,
    140,
    160,
    180,  //10
    200,
    225,
    250,
    300,
    350,
    400,
    500, 
    600,
    800,
    1500};//20
  ffMr_R2point8cut->Draw();
  int numBinsMr = sizeof(MrBins)/sizeof(Double_t)-1;
  //ffMr_R2point8cut->Rebin(4);totBGMR->Rebin(4);egMr_R2point8cut->Rebin(4);ggMr_R2point8cut->Rebin(4);
  TH1F* ffMr_R2point8cut_new=(TH1F*)ffMr_R2point8cut->Rebin(numBinsMr,"ffMr_R2point8cut_new",MrBins);
  TH1F* totBGMR_new=(TH1F*)totBGMR->Rebin(numBinsMr,"totBGMR_new",MrBins);
  TH1F* egMr_R2point8cut_new=(TH1F*)egMr_R2point8cut->Rebin(numBinsMr,"egMr_R2point8cut_new",MrBins);
  TH1F* ggMr_R2point8cut_new=(TH1F*)ggMr_R2point8cut->Rebin(numBinsMr,"ggMr_R2point8cut_new",MrBins);

  TH1F* ggMr_R2point8cut_overBinWidth = (TH1F*)ggMr_R2point8cut_new->Clone();
  TH1F* egMr_R2point8cut_overBinWidth = (TH1F*)egMr_R2point8cut_new->Clone();
  TH1F* ffMr_R2point8cut_overBinWidth = (TH1F*)ffMr_R2point8cut_new->Clone();
  TH1F* totBGMR_overBinWidth = (TH1F*)totBGMR_new->Clone();

  ggMr_R2point8cut_overBinWidth->Sumw2();
  egMr_R2point8cut_overBinWidth->Sumw2();
  ffMr_R2point8cut_overBinWidth->Sumw2();
  totBGMR_overBinWidth->Sumw2();
 
  for(int i=0;i<ggMr_R2point8cut_new->GetNbinsX()+1;i++){
    ggMr_R2point8cut_overBinWidth->SetBinContent(i,ggMr_R2point8cut_overBinWidth->GetBinContent(i)/ggMr_R2point8cut_overBinWidth->GetBinWidth(i));
    egMr_R2point8cut_overBinWidth->SetBinContent(i,egMr_R2point8cut_overBinWidth->GetBinContent(i)/egMr_R2point8cut_overBinWidth->GetBinWidth(i));
    ffMr_R2point8cut_overBinWidth->SetBinContent(i,ffMr_R2point8cut_overBinWidth->GetBinContent(i)/ffMr_R2point8cut_overBinWidth->GetBinWidth(i));
    totBGMR_overBinWidth->SetBinContent(i,totBGMR_overBinWidth->GetBinContent(i)/totBGMR_overBinWidth->GetBinWidth(i));
    ggMr_R2point8cut_overBinWidth->SetBinError(i,ggMr_R2point8cut_overBinWidth->GetBinError(i)/ggMr_R2point8cut_overBinWidth->GetBinWidth(i));
    egMr_R2point8cut_overBinWidth->SetBinError(i,egMr_R2point8cut_overBinWidth->GetBinError(i)/egMr_R2point8cut_overBinWidth->GetBinWidth(i));
    ffMr_R2point8cut_overBinWidth->SetBinError(i,ffMr_R2point8cut_overBinWidth->GetBinError(i)/ffMr_R2point8cut_overBinWidth->GetBinWidth(i));
    totBGMR_overBinWidth->SetBinError(i,totBGMR_overBinWidth->GetBinError(i)/totBGMR_overBinWidth->GetBinWidth(i));
  }
  egMr_R2point8cut_overBinWidth->GetYaxis()->SetTitle("Number of Events / GeV");
  egMr_R2point8cut_overBinWidth->GetXaxis()->SetTitle("Razor M_{R}");
  egMr_R2point8cut_overBinWidth->SetTitle("");
  egMr_R2point8cut_overBinWidth->SetFillColor(kCyan);
  egMr_R2point8cut_overBinWidth->SetFillStyle(3005);
  egMr_R2point8cut_overBinWidth->SetStats(0);
  egMr_R2point8cut_overBinWidth->SetMarkerSize(0);
  ffMr_R2point8cut_overBinWidth->SetFillColor(kAzure);
  ffMr_R2point8cut_overBinWidth->SetFillStyle(3004);
  ffMr_R2point8cut_overBinWidth->SetMarkerSize(0);
  ffMr_R2point8cut_overBinWidth->SetStats(0);
  ffMr_R2point8cut_overBinWidth->SetTitle("");
  totBGMR_overBinWidth->SetMarkerSize(0);
  totBGMR_overBinWidth->SetFillColor(kRed);
  totBGMR_overBinWidth->SetFillStyle(3003);
  totBGMR_overBinWidth->SetStats(0);
  float max2 = ggMr_R2point8cut_overBinWidth->GetMaximum()+2500;
  //float max2 =1000;
  egMr_R2point8cut_overBinWidth->SetMaximum(max2);
  egMr_R2point8cut_overBinWidth->Draw("histo");
  //ffMr_R2point8cut_overBinWidth->Draw("BSAME");
  ffMr_R2point8cut_overBinWidth->Draw("histoSAME");
  //egMr_R2point8cut_overBinWidth->Draw("BSAME");
  egMr_R2point8cut_overBinWidth->Draw("histoSAME");
  totBGMR_overBinWidth->Draw("E2SAME");
  ggMr_R2point8cut_overBinWidth->Draw("PESAME");
  // ggSig720New->Draw("SAME");
  //ggSig800New->Draw("SAME");
  TLegend *legMR = new TLegend(.6,.55,.9,.75,"","brNDC");
  legMR->AddEntry(ggMr_R2point8cut_overBinWidth,"#gamma#gamma Candidate Sample","lpf");
  legMR->AddEntry(totBGMR_overBinWidth,"QCD + Electroweak","f");
  legMR->AddEntry(ffMr_R2point8cut_overBinWidth,"QCD (ff sample)","f");
  legMR->AddEntry(egMr_R2point8cut_overBinWidth,"Electroweak","f");
  legMR->SetFillColor(kWhite);
  legMR->Draw("SAME");
  Text = new TPaveText(.6,.76,.9,.97,"NDC");
  Text->AddText("CMS Preliminary");
  Text->AddText("");
  Text->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  Text->AddText("No Jet Requirement");
  Text->SetFillColor(0);
  Text->SetBorderSize(0);
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MR_ggANDffNoReweight.png");
  fout.cd();ffMr_R2point8cut_overBinWidth->Write("ffMRBinnedReweightedAndScaled");ggMr_R2point8cut_overBinWidth->Write("ggMRReweightedAndScaled");fin.cd();

  TH1F* ggR2 = (TH1F*)fin.Get("ggR2");//ggR2->SetStats(0);
  TH1F* egR2 = (TH1F*)fin.Get("egR2");
  TH1F* ffR2 = (TH1F*)fin.Get("ffR2");
  //TH1F* eeR2 = (TH1F*)fin.Get("eeR2");
  egScale = ggR2->Integral(1,10)/egR2->Integral(1,10);
  egR2->Scale(egScale,"");
  egR2->Scale(0.0145,"");
  ffR2->Scale( (ggR2->Integral(1,60)-egR2->Integral(1,60))/ffR2->Integral(1,60) , "");
  TH1F* totBGfR2 = (TH1F*)ffR2->Clone();
  totBGfR2->Add(egR2, ffR2,1,1);
  ffR2->Rebin(2);totBGfR2->Rebin(2);egR2->Rebin(2);ggR2->Rebin(2);
  ffR2->SetFillColor(kAzure);
  ffR2->SetFillStyle(3004);
  ffR2->SetMarkerSize(0);
  ffR2->SetStats(0);
  ffR2->SetTitle("");
  egR2->GetYaxis()->SetTitle("Number of Events");
  egR2->GetXaxis()->SetTitle("Razor R^{2}");
  egR2->SetTitle("");
  egR2->SetFillColor(kCyan);
  egR2->SetFillStyle(3005);
  egR2->SetStats(0);
  totBGfR2->SetMarkerSize(0);
  totBGfR2->SetFillColor(kRed);
  totBGfR2->SetFillStyle(3003);
  totBGfR2->SetStats(0);
  float maxR2 = ggR2->GetMaximum()+4000;
  egR2->SetMaximum(maxR2);
  egR2->Draw("B");
  ffR2->Draw("BSAME");
  ffR2->Draw("SAME");
  egR2->Draw("BSAME");
  egR2->Draw("SAME");
  totBGfR2->Draw("E2SAME");
  ggR2->Draw("PESAME");
  // ggSig720New->Draw("SAME");
  //ggSig800New->Draw("SAME");
  TLegend *legR2 = new TLegend(.55,.55,.85,.75,"","brNDC");
  legR2->AddEntry(ggR2,"#gamma#gamma Candidate Sample","lpf");
  legR2->AddEntry(totBGfR2,"QCD + Electroweak","f");
  legR2->AddEntry(ffR2,"QCD (ff sample)","f");
  legR2->AddEntry(egR2,"Electroweak","f");
  legR2->SetFillColor(kWhite);
  legR2->Draw("SAME");
  Text = new TPaveText(.55,.76,.85,.97,"NDC");
  Text->AddText("CMS Preliminary");
  Text->AddText("");
  Text->AddText("#sqrt{s} = 8 TeV, #intLdt = 19.5 fb^{-1}");
  Text->AddText("No Jet Requirement");
  Text->SetFillColor(0);
  Text->SetBorderSize(0);
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/R2_ggANDffNoReweight.png");
  fout.cd();ffR2->Write("ffR2BinnedReweightedAndScaled");ggR2->Write("ggR2ReweightedAndScaled");fin.cd();


  gStyle->SetErrorX(1);

  //-------------------end of log scale plots
  

  //Plots that use linear scale here:
  c1->cd();
  c1->SetLogy(0);
  gStyle->SetErrorX();
  gStyle->SetOptStat(0);
  //DiEMPt Ratio
  TH1F* ggeeDiEMPtRatio=(TH1F*)fin.Get("ggeeDiEMPtRatio");
  ggeeDiEMPtRatio->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiEMPtRatio_ggee.png");
  fout.cd();
  ggeeDiEMPtRatio->Write();
  fin.cd();

  TH1F* ggeeSideBandDiEMPtRatio=(TH1F*)fin.Get("ggeeSideBandDiEMPtRatio");
  ggeeSideBandDiEMPtRatio->Draw();
  ggeeSideBandDiEMPtRatio->GetXaxis()->SetRangeUser(0,400);
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiEMPtRatio_ggeeSideband.png");
  fout.cd();
  ggeeSideBandDiEMPtRatio->Write();
  fin.cd();

  TH1F* ggffDiEMPtRatio=(TH1F*)fin.Get("ggffDiEMPtRatio");
  ggffDiEMPtRatio->Draw();
  ggffDiEMPtRatio->GetXaxis()->SetRangeUser(0,400);
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiEMPtRatio_ggff.png");
  fout.cd();
  ggffDiEMPtRatio->Write();
  fin.cd();

  TH1F* ggeeDiEMPtRatio_JetReq=(TH1F*)fin.Get("ggeeDiEMPtRatio_JetReq");
  ggeeDiEMPtRatio_JetReq->Draw();
  ggeeDiEMPtRatio_JetReq->GetXaxis()->SetRangeUser(0,400);
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiEMPtRatio_ggee_JetReq.png");
  fout.cd();
  ggeeDiEMPtRatio_JetReq->Write();
  fin.cd();

  TH1F* ggeeSideBandDiEMPtRatio_JetReq=(TH1F*)fin.Get("ggeeSideBandDiEMPtRatio_JetReq");
  ggeeSideBandDiEMPtRatio_JetReq->Draw();
  ggeeSideBandDiEMPtRatio_JetReq->GetXaxis()->SetRangeUser(0,400);
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiEMPtRatio_ggeeSideband_JetReq.png");
  fout.cd();
  ggeeSideBandDiEMPtRatio_JetReq->Write();
  fin.cd();

  TH1F* ggffDiEMPtRatio_JetReq=(TH1F*)fin.Get("ggffDiEMPtRatio_JetReq");
  ggffDiEMPtRatio_JetReq->Draw();
  ggffDiEMPtRatio_JetReq->GetXaxis()->SetRangeUser(0,400);
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiEMPtRatio_ggff_JetReq.png");
  fout.cd();
  ggffDiEMPtRatio_JetReq->Write();
  fin.cd();


  //DiJetPt Ratio
  TH1F* ggeeDiJetPtRatio=(TH1F*)fin.Get("ggeeDiJetPtRatio");
  ggeeDiJetPtRatio->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggee.png");
  fout.cd();
  ggeeDiJetPtRatio->Write();
  fin.cd();

  TH1F* ggeeDiJetPtRatio_0Jet=(TH1F*)fin.Get("ggeeDiJetPtRatio_0Jet");
  ggeeDiJetPtRatio_0Jet->GetXaxis()->SetRangeUser(0,599);
  ggeeDiJetPtRatio_0Jet->SetTitle("");
  ggeeDiJetPtRatio_0Jet->Draw();
  Text_0Jet->Draw();
  Textee->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_ggee_0Jet.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_ggee_0Jet.pdf");
  c3->cd(1);ggeeDiJetPtRatio_0Jet->Draw();c1->cd();
  fout.cd();
  ggeeDiJetPtRatio_0Jet->Write();
  fin.cd();

  TH1F* ggeeDiJetPtRatio_1Jet=(TH1F*)fin.Get("ggeeDiJetPtRatio_1Jet");
  ggeeDiJetPtRatio_1Jet->GetXaxis()->SetRangeUser(0,599);
  ggeeDiJetPtRatio_1Jet->SetTitle("");
  ggeeDiJetPtRatio_1Jet->Draw();
  Text_1Jet->Draw();
  Textee->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_ggee_1Jet.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_ggee_1Jet.pdf");
  c3->cd(2);ggeeDiJetPtRatio_1Jet->Draw();c1->cd();
  fout.cd();
  ggeeDiJetPtRatio_1Jet->Write();
  fin.cd();

  TH1F* ggeeDiJetPtRatio_2Jet=(TH1F*)fin.Get("ggeeDiJetPtRatio_2Jet");
  ggeeDiJetPtRatio_2Jet->GetXaxis()->SetRangeUser(0,599);
  ggeeDiJetPtRatio_2Jet->SetTitle("");
  ggeeDiJetPtRatio_2Jet->Draw();
  Text_2Jet->Draw();
  Textee->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_ggee_2Jet.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_ggee_2Jet.pdf");
  c3->cd(3);ggeeDiJetPtRatio_2Jet->Draw();c1->cd();
  fout.cd();
  ggeeDiJetPtRatio_2Jet->Write();
  fin.cd();

  TH1F* ggeeDiJetPtRatio_JetReq_0Jet=(TH1F*)fin.Get("ggeeDiJetPtRatio_JetReq_0Jet");
  ggeeDiJetPtRatio_JetReq_0Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_JetReq_ggee_0Jet.png");
  c3->cd(4);ggeeDiJetPtRatio_JetReq_0Jet->Draw();c1->cd();
  fout.cd();
  ggeeDiJetPtRatio_JetReq_0Jet->Write();
  fin.cd();

  TH1F* ggeeDiJetPtRatio_JetReq_1Jet=(TH1F*)fin.Get("ggeeDiJetPtRatio_JetReq_1Jet");
  ggeeDiJetPtRatio_JetReq_1Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_JetReq_ggee_1Jet.png");
  c3->cd(5);ggeeDiJetPtRatio_JetReq_1Jet->Draw();c1->cd();
  fout.cd();
  ggeeDiJetPtRatio_JetReq_1Jet->Write();
  fin.cd();

  TH1F* ggeeDiJetPtRatio_JetReq_2Jet=(TH1F*)fin.Get("ggeeDiJetPtRatio_JetReq_2Jet");
  ggeeDiJetPtRatio_JetReq_2Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_JetReq_ggee_2Jet.png");
  c3->cd(6);ggeeDiJetPtRatio_JetReq_2Jet->Draw();
  c3->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggee_all.png");
  c1->cd();
  fout.cd();
  ggeeDiJetPtRatio_JetReq_2Jet->Write();
  fin.cd();

  TH1F* ggeeSideBandLowDiJetPtRatio=(TH1F*)fin.Get("ggeeSideBandLowDiJetPtRatio");
  ggeeSideBandLowDiJetPtRatio->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggeeSideBandLow.png");
  fout.cd();
  ggeeSideBandLowDiJetPtRatio->Write();
  fin.cd();

  TH1F* ggeeSideBandLowDiJetPtRatio_0Jet=(TH1F*)fin.Get("ggeeSideBandLowDiJetPtRatio_0Jet");
  ggeeSideBandLowDiJetPtRatio_0Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggeeSideBandLow_0Jet.png");
  c3->cd(1);ggeeSideBandLowDiJetPtRatio_0Jet->Draw();c1->cd();
  fout.cd();
  ggeeSideBandLowDiJetPtRatio_0Jet->Write();
  fin.cd();

  TH1F* ggeeSideBandLowDiJetPtRatio_1Jet=(TH1F*)fin.Get("ggeeSideBandLowDiJetPtRatio_1Jet");
  ggeeSideBandLowDiJetPtRatio_1Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggeeSideBandLow_1Jet.png");
  c3->cd(2);ggeeSideBandLowDiJetPtRatio_1Jet->Draw();c1->cd();
  fout.cd();
  ggeeSideBandLowDiJetPtRatio_1Jet->Write();
  fin.cd();

  TH1F* ggeeSideBandLowDiJetPtRatio_2Jet=(TH1F*)fin.Get("ggeeSideBandLowDiJetPtRatio_2Jet");
  ggeeSideBandLowDiJetPtRatio_2Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggeeSideBandLow_2Jet.png");
  c3->cd(3);ggeeSideBandLowDiJetPtRatio_2Jet->Draw();c1->cd();
  fout.cd();
  ggeeSideBandLowDiJetPtRatio_2Jet->Write();
  fin.cd();

  TH1F* ggeeSideBandLowDiJetPtRatio_JetReq_0Jet=(TH1F*)fin.Get("ggeeSideBandLowDiJetPtRatio_JetReq_0Jet");
  ggeeSideBandLowDiJetPtRatio_JetReq_0Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_JetReq_ggeeSideBandLow_0Jet.png");
  c3->cd(4);ggeeSideBandLowDiJetPtRatio_JetReq_0Jet->Draw();c1->cd();
  fout.cd();
  ggeeSideBandLowDiJetPtRatio_JetReq_0Jet->Write();
  fin.cd();

  TH1F* ggeeSideBandLowDiJetPtRatio_JetReq_1Jet=(TH1F*)fin.Get("ggeeSideBandLowDiJetPtRatio_JetReq_1Jet");
  ggeeSideBandLowDiJetPtRatio_JetReq_1Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_JetReq_ggeeSideBandLow_1Jet.png");
  c3->cd(5);ggeeSideBandLowDiJetPtRatio_JetReq_1Jet->Draw();c1->cd();
  fout.cd();
  ggeeSideBandLowDiJetPtRatio_JetReq_1Jet->Write();
  fin.cd();

  TH1F* ggeeSideBandLowDiJetPtRatio_JetReq_2Jet=(TH1F*)fin.Get("ggeeSideBandLowDiJetPtRatio_JetReq_2Jet");
  ggeeSideBandLowDiJetPtRatio_JetReq_2Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_JetReq_ggeeSideBandLow_2Jet.png");
  c3->cd(6);ggeeSideBandLowDiJetPtRatio_JetReq_2Jet->Draw();
  c3->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggeeSideBandLow_all.png");
  c1->cd();
  fout.cd();
  ggeeSideBandLowDiJetPtRatio_JetReq_2Jet->Write();
  fin.cd();

  TH1F* ggeeSideBandHighDiJetPtRatio=(TH1F*)fin.Get("ggeeSideBandHighDiJetPtRatio");
  ggeeSideBandHighDiJetPtRatio->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggeeSideBandHigh.png");
  fout.cd();
  ggeeSideBandHighDiJetPtRatio->Write();
  fin.cd();

  TH1F* ggeeSideBandHighDiJetPtRatio_0Jet=(TH1F*)fin.Get("ggeeSideBandHighDiJetPtRatio_0Jet");
  ggeeSideBandHighDiJetPtRatio_0Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggeeSideBandHigh_0Jet.png");
  c3->cd(1);ggeeSideBandHighDiJetPtRatio_0Jet->Draw();c1->cd();
  fout.cd();
  ggeeSideBandHighDiJetPtRatio_0Jet->Write();
  fin.cd();

  TH1F* ggeeSideBandHighDiJetPtRatio_1Jet=(TH1F*)fin.Get("ggeeSideBandHighDiJetPtRatio_1Jet");
  ggeeSideBandHighDiJetPtRatio_1Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggeeSideBandHigh_1Jet.png");
  c3->cd(2);ggeeSideBandHighDiJetPtRatio_1Jet->Draw();c1->cd();
  fout.cd();
  ggeeSideBandHighDiJetPtRatio_1Jet->Write();
  fin.cd();

  TH1F* ggeeSideBandHighDiJetPtRatio_2Jet=(TH1F*)fin.Get("ggeeSideBandHighDiJetPtRatio_2Jet");
  ggeeSideBandHighDiJetPtRatio_2Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggeeSideBandHigh_2Jet.png");
  c3->cd(3);ggeeSideBandHighDiJetPtRatio_2Jet->Draw();c1->cd();
  fout.cd();
  ggeeSideBandHighDiJetPtRatio_2Jet->Write();
  fin.cd();

  TH1F* ggeeSideBandHighDiJetPtRatio_JetReq_0Jet=(TH1F*)fin.Get("ggeeSideBandHighDiJetPtRatio_JetReq_0Jet");
  ggeeSideBandHighDiJetPtRatio_JetReq_0Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_JetReq_ggeeSideBandHigh_0Jet.png");
  c3->cd(4);ggeeSideBandHighDiJetPtRatio_JetReq_0Jet->Draw();c1->cd();
  fout.cd();
  ggeeSideBandHighDiJetPtRatio_JetReq_0Jet->Write();
  fin.cd();

  TH1F* ggeeSideBandHighDiJetPtRatio_JetReq_1Jet=(TH1F*)fin.Get("ggeeSideBandHighDiJetPtRatio_JetReq_1Jet");
  ggeeSideBandHighDiJetPtRatio_JetReq_1Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_JetReq_ggeeSideBandHigh_1Jet.png");
  c3->cd(5);ggeeSideBandHighDiJetPtRatio_JetReq_1Jet->Draw();c1->cd();
  fout.cd();
  ggeeSideBandHighDiJetPtRatio_JetReq_1Jet->Write();
  fin.cd();

  TH1F* ggeeSideBandHighDiJetPtRatio_JetReq_2Jet=(TH1F*)fin.Get("ggeeSideBandHighDiJetPtRatio_JetReq_2Jet");
  ggeeSideBandHighDiJetPtRatio_JetReq_2Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_JetReq_ggeeSideBandHigh_2Jet.png");
  c3->cd(6);ggeeSideBandHighDiJetPtRatio_JetReq_2Jet->Draw();
  c3->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggeeSideBandHigh_all.png");
  c1->cd();
  fout.cd();
  ggeeSideBandHighDiJetPtRatio_JetReq_2Jet->Write();
  fin.cd();

  //now ff


  TH1F* ggffDiJetPtRatio_0Jet=(TH1F*)fin.Get("ggffDiJetPtRatio_0Jet");
  ggffDiJetPtRatio_0Jet->GetXaxis()->SetRangeUser(0,599);
  //ggffDiJetPtRatio_0Jet->GetYaxis()->SetRangeUser(.5,1.7);
  ggffDiJetPtRatio_0Jet->SetTitle("");
  ggffDiJetPtRatio_0Jet->Draw();
  Text_0Jet->Draw();
  Textff->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_ggff_0Jet.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_ggff_0Jet.pdf");
  c3->cd(1);ggffDiJetPtRatio_0Jet->GetXaxis()->SetRangeUser(0,599);ggffDiJetPtRatio_0Jet->Draw();c1->cd();
  fout.cd();
  ggffDiJetPtRatio_0Jet->Write();
  fin.cd();

  TH1F* ggffDiJetPtRatio_1Jet=(TH1F*)fin.Get("ggffDiJetPtRatio_1Jet");
  ggffDiJetPtRatio_1Jet->GetXaxis()->SetRangeUser(0,599);
  //ggffDiJetPtRatio_1Jet->GetYaxis()->SetRangeUser(.5,1.7);
  ggffDiJetPtRatio_1Jet->SetTitle("");
  ggffDiJetPtRatio_1Jet->Draw();
  Text_1Jet->Draw();
  Textff->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_ggff_1Jet.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_ggff_1Jet.pdf");
  c3->cd(2);ggffDiJetPtRatio_1Jet->GetXaxis()->SetRangeUser(0,599);ggffDiJetPtRatio_1Jet->Draw();c1->cd();
  fout.cd();
  ggffDiJetPtRatio_1Jet->Write();
  fin.cd();

  TH1F* ggffDiJetPtRatio_2Jet=(TH1F*)fin.Get("ggffDiJetPtRatio_2Jet");
  ggffDiJetPtRatio_2Jet->GetXaxis()->SetRangeUser(0,599);
  //ggffDiJetPtRatio_2Jet->GetYaxis()->SetRangeUser(.5,1.7);
  ggffDiJetPtRatio_2Jet->SetTitle("");
  ggffDiJetPtRatio_2Jet->Draw();
  Text_2Jet->Draw();
  Textff->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_ggff_2Jet.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_ggff_2Jet.pdf");
  c3->cd(3);ggffDiJetPtRatio_2Jet->GetXaxis()->SetRangeUser(0,599);ggffDiJetPtRatio_2Jet->Draw();c1->cd();
  fout.cd();
  ggffDiJetPtRatio_2Jet->Write();
  fin.cd();

  //now gammafake
  TH1F* gggammafakeDiJetPtRatio_0Jet=(TH1F*)fin.Get("gggammafakeDiJetPtRatio_0Jet");
  gggammafakeDiJetPtRatio_0Jet->GetXaxis()->SetRangeUser(0,599);
  //gggammafakeDiJetPtRatio_0Jet->GetYaxis()->SetRangeUser(.5,1.7);
  gggammafakeDiJetPtRatio_0Jet->SetTitle("");
  gggammafakeDiJetPtRatio_0Jet->Draw();
  Text_0Jet->Draw();
  Textgammafake->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_gggammafake_0Jet.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_gggammafake_0Jet.pdf");
  c3->cd(1);gggammafakeDiJetPtRatio_0Jet->GetXaxis()->SetRangeUser(0,599);gggammafakeDiJetPtRatio_0Jet->Draw();c1->cd();
  fout.cd();
  gggammafakeDiJetPtRatio_0Jet->Write();
  fin.cd();

  TH1F* gggammafakeDiJetPtRatio_1Jet=(TH1F*)fin.Get("gggammafakeDiJetPtRatio_1Jet");
  gggammafakeDiJetPtRatio_1Jet->GetXaxis()->SetRangeUser(0,599);
  //gggammafakeDiJetPtRatio_1Jet->GetYaxis()->SetRangeUser(.5,1.7);
  gggammafakeDiJetPtRatio_1Jet->SetTitle("");
  gggammafakeDiJetPtRatio_1Jet->Draw();
  Text_1Jet->Draw();
  Textgammafake->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_gggammafake_1Jet.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_gggammafake_1Jet.pdf");
  c3->cd(2);gggammafakeDiJetPtRatio_1Jet->GetXaxis()->SetRangeUser(0,599);gggammafakeDiJetPtRatio_1Jet->Draw();c1->cd();
  fout.cd();
  gggammafakeDiJetPtRatio_1Jet->Write();
  fin.cd();

  TH1F* gggammafakeDiJetPtRatio_2Jet=(TH1F*)fin.Get("gggammafakeDiJetPtRatio_2Jet");
  gggammafakeDiJetPtRatio_2Jet->GetXaxis()->SetRangeUser(0,599);
  //gggammafakeDiJetPtRatio_2Jet->GetYaxis()->SetRangeUser(.5,1.7);
  gggammafakeDiJetPtRatio_2Jet->SetTitle("");
  gggammafakeDiJetPtRatio_2Jet->Draw();
  Text_2Jet->Draw();
  Textgammafake->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_gggammafake_2Jet.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure15_DiJetPtRatio_gggammafake_2Jet.pdf");
  c3->cd(3);gggammafakeDiJetPtRatio_2Jet->GetXaxis()->SetRangeUser(0,599);gggammafakeDiJetPtRatio_2Jet->Draw();c1->cd();
  fout.cd();
  gggammafakeDiJetPtRatio_2Jet->Write();
  fin.cd();
  //

  TH1F* ggffDiJetPtRatio_JetReq_0Jet=(TH1F*)fin.Get("ggffDiJetPtRatio_JetReq_0Jet");
  ggffDiJetPtRatio_JetReq_0Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_JetReq_ggff_0Jet.png");
  c3->cd(4);ggffDiJetPtRatio_JetReq_0Jet->Draw();c1->cd();
  fout.cd();
  ggffDiJetPtRatio_JetReq_0Jet->Write();
  fin.cd();

  TH1F* ggffDiJetPtRatio_JetReq_1Jet=(TH1F*)fin.Get("ggffDiJetPtRatio_JetReq_1Jet");
  ggffDiJetPtRatio_JetReq_1Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_JetReq_ggff_1Jet.png");
  c3->cd(5);ggffDiJetPtRatio_JetReq_1Jet->Draw();c1->cd();
  fout.cd();
  ggffDiJetPtRatio_JetReq_1Jet->Write();
  fin.cd();

  TH1F* ggffDiJetPtRatio_JetReq_2Jet=(TH1F*)fin.Get("ggffDiJetPtRatio_JetReq_2Jet");
  ggffDiJetPtRatio_JetReq_2Jet->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_JetReq_ggff_2Jet.png");
  c3->cd(6);ggffDiJetPtRatio_JetReq_2Jet->Draw();
  c3->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggff_all.png");
  c1->cd();
  fout.cd();
  ggffDiJetPtRatio_JetReq_2Jet->Write();
  fin.cd();
  //
  TH1F* ggeeSideBandDiJetPtRatio=(TH1F*)fin.Get("ggeeSideBandDiJetPtRatio");
  ggeeSideBandDiJetPtRatio->Draw();
  ggeeSideBandDiJetPtRatio->GetXaxis()->SetRangeUser(0,400);
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggeeSideband.png");
  fout.cd();
  ggeeSideBandDiJetPtRatio->Write();
  fin.cd();

  TH1F* ggffDiJetPtRatio=(TH1F*)fin.Get("ggffDiJetPtRatio");
  ggffDiJetPtRatio->Draw();
  ggffDiJetPtRatio->GetXaxis()->SetRangeUser(0,400);
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggff.png");
  fout.cd();
  ggffDiJetPtRatio->Write();
  fin.cd();

  TH1F* ggeeDiJetPtRatio_JetReq=(TH1F*)fin.Get("ggeeDiJetPtRatio_JetReq");
  ggeeDiJetPtRatio_JetReq->Draw();
  ggeeDiJetPtRatio_JetReq->GetXaxis()->SetRangeUser(0,400);
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggee_JetReq.png");
  fout.cd();
  ggeeDiJetPtRatio_JetReq->Write();
  fin.cd();

  TH1F* ggeeSideBandDiJetPtRatio_JetReq=(TH1F*)fin.Get("ggeeSideBandDiJetPtRatio_JetReq");
  ggeeSideBandDiJetPtRatio_JetReq->Draw();
  ggeeSideBandDiJetPtRatio_JetReq->GetXaxis()->SetRangeUser(0,400);
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggeeSideband_JetReq.png");
  fout.cd();
  ggeeSideBandDiJetPtRatio_JetReq->Write();
  fin.cd();

  TH1F* ggffDiJetPtRatio_JetReq=(TH1F*)fin.Get("ggffDiJetPtRatio_JetReq");
  ggffDiJetPtRatio_JetReq->Draw();
  ggffDiJetPtRatio_JetReq->GetXaxis()->SetRangeUser(0,400);
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/DiJetPtRatio_ggff_JetReq.png");
  fout.cd();
  ggffDiJetPtRatio_JetReq->Write();
  fin.cd();

  c1->SetLogy(1);

  //BDT
  TH1F* ggBDT = (TH1F*)fin.Get("ggMva_BDT");ggBDT->Sumw2();
  TH1F* ffBDT = (TH1F*)fin.Get("ffMva_BDT");ffBDT->Sumw2();
  TH1F* gammafakeBDT = (TH1F*)fin.Get("gammafakeMva_BDT");gammafakeBDT->Sumw2();
  TH1F* egBDT = (TH1F*)fin.Get("egMva_BDT");egBDT->Sumw2();

  egBDT->SetFillColor(kGreen);
  gammafakeBDT->SetFillColor(kCyan);
  ffBDT->SetFillColor(kGray+1);
  ffPercent=1;gammafakePercent=0.;
  float bdtscale = ffBDT->Integral()/gammafakeBDT->Integral();
  gammafakeBDT->Scale(bdtscale);
  ffBDT->Scale(ffPercent);gammafakeBDT->Scale(gammafakePercent);
  cout<<"fakerate:"<<FakeRate<<endl;
  egBDT->Scale(FakeRate);

  int binLow=ggBDT->FindBin(-0.8), binHigh=ggBDT->FindBin(-0.2);
  bdtscale = (ggBDT->Integral(binLow,binHigh)-egBDT->Integral(binLow,binHigh))/(ffBDT->Integral(binLow,binHigh)+gammafakeBDT->Integral(binLow,binHigh));

  ffBDT->Scale(bdtscale);gammafakeBDT->Scale(bdtscale);
  TH1F *bdtBG = (TH1F*)ffBDT->Clone();bdtBG->Add(gammafakeBDT);bdtBG->Add(egBDT);bdtBG->SetFillStyle(3154);bdtBG->SetFillColor(kRed);bdtBG->SetMarkerStyle(0);

  THStack *bdtStack = new THStack("bdtStack","");
  bdtStack->Add(egBDT);
  bdtStack->Add(gammafakeBDT);
  bdtStack->Add(ffBDT);
  bdtStack->Draw("histo");
  bdtStack->SetMaximum(4e4);bdtStack->SetMinimum(.12);
  bdtStack->GetHistogram()->GetXaxis()->SetRangeUser(-0.8,0.8);
  bdtStack->Draw("histo");
  bdtBG->Draw("E2SAMES");
  ggBDT->Draw("PESAMES");
  TLegend *legBDT = new TLegend(.65,.52,.8,.82,"","brNDC");legBDT->SetFillStyle(0);legBDT->SetBorderSize(0);
  legBDT->AddEntry(ggBDT,"#gamma#gamma","lp");
  legBDT->AddEntry(bdtBG,"BG err","f");
  legBDT->AddEntry(ffBDT,"ff","f");
  legBDT->AddEntry(gammafakeBDT,"#gammaf","f");
  legBDT->AddEntry(egBDT,"e#gamma","f");
  legBDT->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MVA_BDT_ggANDcomb.png");

  //SVM
  TH1F* ggSVM = (TH1F*)fin.Get("ggMva_SVM");ggSVM->Sumw2();
  TH1F* ffSVM = (TH1F*)fin.Get("ffMva_SVM");ffSVM->Sumw2();
  TH1F* gammafakeSVM = (TH1F*)fin.Get("gammafakeMva_SVM");gammafakeSVM->Sumw2();
  TH1F* egSVM = (TH1F*)fin.Get("egMva_SVM");egSVM->Sumw2();

  egSVM->SetFillColor(kGreen);
  gammafakeSVM->SetFillColor(kCyan);
  ffSVM->SetFillColor(kGray+1);
  
  float svmscale = ffSVM->Integral()/gammafakeSVM->Integral();
  gammafakeSVM->Scale(svmscale);
  ffSVM->Scale(ffPercent);gammafakeSVM->Scale(gammafakePercent);
  egSVM->Scale(FakeRate);

  binLow=ggSVM->FindBin(0.0), binHigh=ggSVM->FindBin(0.3);
  svmscale = (ggSVM->Integral(binLow,binHigh)-egSVM->Integral(binLow,binHigh))/(ffSVM->Integral(binLow,binHigh)+gammafakeSVM->Integral(binLow,binHigh));

  ffSVM->Scale(svmscale);gammafakeSVM->Scale(svmscale);
  TH1F *svmBG = (TH1F*)ffSVM->Clone();svmBG->Add(gammafakeSVM);svmBG->Add(egSVM);svmBG->SetFillStyle(3154);svmBG->SetFillColor(kRed);svmBG->SetMarkerStyle(0);

  THStack *svmStack = new THStack("svmStack","");
  svmStack->Add(egSVM);
  svmStack->Add(gammafakeSVM);
  svmStack->Add(ffSVM);
  svmStack->Draw("histo");
  svmStack->SetMaximum(6e4);svmStack->SetMinimum(.12);
  svmStack->GetHistogram()->GetXaxis()->SetRangeUser(0.,1.0);
  svmStack->Draw("histo");
  svmBG->Draw("E2SAMES");
  ggSVM->Draw("PESAMES");
  TLegend *legSVM = new TLegend(.65,.52,.8,.82,"","brNDC");legSVM->SetFillStyle(0);legSVM->SetBorderSize(0);
  legSVM->AddEntry(ggSVM,"#gamma#gamma","lp");
  legSVM->AddEntry(svmBG,"BG err","f");
  legSVM->AddEntry(ffSVM,"ff","f");
  legSVM->AddEntry(gammafakeSVM,"#gammaf","f");
  legSVM->AddEntry(egSVM,"e#gamma","f");
  legSVM->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MVA_SVM_ggANDcomb.png");

  //LD
  TH1F* ggLD = (TH1F*)fin.Get("ggMva_LD");ggLD->Sumw2();
  TH1F* ffLD = (TH1F*)fin.Get("ffMva_LD");ffLD->Sumw2();
  TH1F* gammafakeLD = (TH1F*)fin.Get("gammafakeMva_LD");gammafakeLD->Sumw2();
  TH1F* egLD = (TH1F*)fin.Get("egMva_LD");egLD->Sumw2();

  egLD->SetFillColor(kGreen);
  gammafakeLD->SetFillColor(kCyan);
  ffLD->SetFillColor(kGray+1);
  
  float ldscale = ffLD->Integral()/gammafakeLD->Integral();
  gammafakeLD->Scale(ldscale);
  ffLD->Scale(ffPercent);gammafakeLD->Scale(gammafakePercent);
  egLD->Scale(FakeRate);

  binLow=ggLD->FindBin(-0.6), binHigh=ggLD->FindBin(-0.1);
  ldscale = (ggLD->Integral(binLow,binHigh)-egLD->Integral(binLow,binHigh))/(ffLD->Integral(binLow,binHigh)+gammafakeLD->Integral(binLow,binHigh));

  ffLD->Scale(ldscale);gammafakeLD->Scale(ldscale);
  TH1F *ldBG = (TH1F*)ffLD->Clone();ldBG->Add(gammafakeLD);ldBG->Add(egLD);ldBG->SetFillStyle(3154);ldBG->SetFillColor(kRed);ldBG->SetMarkerStyle(0);

  THStack *ldStack = new THStack("ldStack","");
  ldStack->Add(egLD);
  ldStack->Add(gammafakeLD);
  ldStack->Add(ffLD);
  ldStack->Draw("histo");
  ldStack->SetMaximum(4e4);ldStack->SetMinimum(.02);
  ldStack->GetHistogram()->GetXaxis()->SetRangeUser(-0.6,1.);
  ldStack->Draw("histo");
  ldBG->Draw("E2SAMES");
  ggLD->Draw("PESAMES");
  TLegend *legLD = new TLegend(.65,.52,.8,.82,"","brNDC");legLD->SetFillStyle(0);legLD->SetBorderSize(0);
  legLD->AddEntry(ggLD,"#gamma#gamma","lp");
  legLD->AddEntry(ldBG,"BG err","f");
  legLD->AddEntry(ffLD,"ff","f");
  legLD->AddEntry(gammafakeLD,"#gammaf","f");
  legLD->AddEntry(egLD,"e#gamma","f");
  legLD->Draw("SAMES");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MVA_LD_ggANDcomb.png");


  c1->SetLogy(0);

/*
    TH1F *ggeeDiEMPtRatio = new TH1F("ggeeDiEMPtRatio","ggDiEMPt/eeDiEMPt - No Jet Requirement;DiEMP_{T};DiEMP_{T} Ratio",200,0.,1000.);
    TH1F *ggffDiEMPtRatio = new TH1F("ggffDiEMPtRatio","ggDiEMPt/ffDiEMPt - No Jet Requirement;DiEMP_{T};DiEMP_{T} Ratio",200,0.,1000.);
    TH1F *ggeeDiEMPtRatio_JetReq = new TH1F("ggeeDiEMPtRatio_JetReq","ggDiEMPt/eeDiEMPt - >=1 Jet Requirement;DiEMP_{T};DiEMP_{T} Ratio",200,0.,1000.);
    TH1F *ggffDiEMPtRatio_JetReq = new TH1F("ggffDiEMPtRatio_JetReq","ggDiEMPt/ffDiEMPt - >=1 Jet Requirement;DiEMP_{T};DiEMP_{T} Ratio",200,0.,1000.);
    ggDiEMPt->Rebin(2);
    eeDiEMPt->Rebin(2);
    ffDiEMPt->Rebin(2);
    ggDiEMPt_JetReq->Rebin(2);
    eeDiEMPt_JetReq->Rebin(2);
    ffDiEMPt_JetReq->Rebin(2);
    if(ggDiEMPt->GetSumw2N()==0)ggDiEMPt->Sumw2();
    if(ffDiEMPt->GetSumw2N()==0)ffDiEMPt->Sumw2();
    if(eeDiEMPt->GetSumw2N()==0)eeDiEMPt->Sumw2();
    if(ggDiEMPt_JetReq->GetSumw2N()==0)ggDiEMPt_JetReq->Sumw2();
    if(ffDiEMPt_JetReq->GetSumw2N()==0)ffDiEMPt_JetReq->Sumw2();
    if(eeDiEMPt_JetReq->GetSumw2N()==0)eeDiEMPt_JetReq->Sumw2();
    float ggInt=ggDiEMPt->Integral(),eeInt=eeDiEMPt->Integral(),ffInt=ffDiEMPt->Integral();
    float scaleggee=ggInt/eeInt;
    float scaleggff=ggInt/ffInt;
    eeDiEMPt->Scale(scaleggee);
    ffDiEMPt->Scale(scaleggff);
    float ggInt_JetReq=ggDiEMPt_JetReq->Integral(),eeInt_JetReq=eeDiEMPt_JetReq->Integral(),ffInt_JetReq=ffDiEMPt_JetReq->Integral();
    float scaleggee_JetReq=ggInt_JetReq/eeInt_JetReq;
    float scaleggff_JetReq=ggInt_JetReq/ffInt_JetReq;
    eeDiEMPt_JetReq->Scale(scaleggee_JetReq);
    ffDiEMPt_JetReq->Scale(scaleggff_JetReq);

    ggeeDiEMPtRatio->Divide(ggDiEMPt,eeDiEMPt,1,1,"");
    ggffDiEMPtRatio->Divide(ggDiEMPt,ffDiEMPt,1,1,"");
    ggeeDiEMPtRatio_JetReq->Divide(ggDiEMPt_JetReq,eeDiEMPt_JetReq,1,1,"");
    ggffDiEMPtRatio_JetReq->Divide(ggDiEMPt_JetReq,ffDiEMPt_JetReq,1,1,"");
    ggeeDiEMPtRatio->Draw();
    //ggeeDiEMPtRatio->SetAxisRange(0,150,"");
    c1->Update();
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggeeDiEMPtRatio.png");
    ggffDiEMPtRatio->Draw();
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggffDiEMPtRatio.png");
    ggeeDiEMPtRatio_JetReq->Draw();
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggeeDiEMPtRatio_JetReq.png");
    ggffDiEMPtRatio_JetReq->Draw();
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ggffDiEMPtRatio_JetReq.png");
  */
  //end DiEMPt Ratio Plots


  //Invariant mass plots
  //TH1F* ggInvarMass=(TH1F*)fin.Get("ggInvarMass");
  ggInvarMass->GetXaxis()->SetRangeUser(75,175);
  ggInvarMass->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMass_gg.png");
  fout.cd();ggInvarMass->Write("InvarMass_gg");fin.cd();

  //TH1F* egInvarMass=(TH1F*)fin.Get("egInvarMass");
  egInvarMass->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMass_eg.png");
  fout.cd();egInvarMass->Write("InvarMass_eg");fin.cd();
  TH1F* eeInvarMass=(TH1F*)fin.Get("eeInvarMass");
  eeInvarMass->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMass_ee.png");
  fout.cd();eeInvarMass->Write("InvarMass_ee");fin.cd();
  //TH1F* eeInvarMassFullRange=(TH1F*)fin.Get("eeInvarMassFullRange");
  eeInvarMassFullRange->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMass_eeFullRange.png");
  fout.cd();eeInvarMassFullRange->Write("InvarMass_eeFullRange");fin.cd();
  TH1F* ffInvarMass=(TH1F*)fin.Get("ffInvarMass");
  ffInvarMass->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMass_ff.png");
  fout.cd();ffInvarMass->Write("InvarMass_ff");fin.cd();
  TH1F* ggInvarMass_JetReq=(TH1F*)fin.Get("ggInvarMass_JetReq");
  ggInvarMass_JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMass_gg_JetReq.png");
  fout.cd();ggInvarMass_JetReq->Write("InvarMass_gg_JetReq");fin.cd();
  TH1F* egInvarMass_JetReq=(TH1F*)fin.Get("egInvarMass_JetReq");
  egInvarMass_JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMass_eg_JetReq.png");
  fout.cd();egInvarMass_JetReq->Write("InvarMass_eg_JetReq");fin.cd();
  TH1F* eeInvarMass_JetReq=(TH1F*)fin.Get("eeInvarMass_JetReq");
  eeInvarMass_JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMass_ee_JetReq.png");
  fout.cd();eeInvarMass_JetReq->Write("InvarMass_ee_JetReq");fin.cd();
  TH1F* ffInvarMass_JetReq=(TH1F*)fin.Get("ffInvarMass_JetReq");
  ffInvarMass_JetReq->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InvarMass_ff_JetReq.png");
  fout.cd();ffInvarMass_JetReq->Write("InvarMass_ff_JetReq");fin.cd();
  
  /*
    if(doPileup){
    //Pileup Effective Areas
    //normal 
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    AvgEcalIsoVsNVertex->GetYaxis()->SetTitleOffset(.8);
    AvgEcalIsoVsNVertex->Draw();
    AvgEcalIsoVsNVertex->Fit("pol1","","",0,20);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ToRun180252/40_25AvgEcalIsoVsNVertex.png");
    fout.cd();AvgEcalIsoVsNVertex->Write("AvgEcalIsoVsNVertex");fin.cd();
    AvgEcalIsoVsRho->Draw();
    AvgEcalIsoVsRho->Fit("pol1","","",0,20);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ToRun180252/40_25AvgEcalIsoVsRho.png");
    fout.cd();AvgEcalIsoVsRho->Write("AvgEcalIsoVsRho");fin.cd();

    AvgHcalIsoVsNVertex->Draw();
    AvgHcalIsoVsNVertex->Fit("pol1","","",0,20);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ToRun180252/40_25AvgHcalIsoVsNVertex.png");
    fout.cd();AvgHcalIsoVsNVertex->Write("AvgHcalIsoVsNVertex");fin.cd();
    AvgHcalIsoVsRho->Draw();
    AvgHcalIsoVsRho->Fit("pol1","","",0,20);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ToRun180252/40_25AvgHcalIsoVsRho.png");
    fout.cd();AvgHcalIsoVsRho->Write("AvgHcalIsoVsRho");fin.cd();

    AvgHOverEVsNVertex->Draw();
    AvgHOverEVsNVertex->Fit("pol1","","",0,20);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ToRun180252/40_25AvgHOverEVsNVertex.png");
    fout.cd();AvgHOverEVsNVertex->Write("AvgHOverEVsNVertex");fin.cd();
    AvgHOverEVsRho->Draw();
    AvgHOverEVsRho->Fit("pol1","","",0,20);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ToRun180252/40_25AvgHOverEVsRho.png");
    fout.cd();AvgHOverEVsRho->Write("AvgHOverEVsRho");fin.cd();
  
    //    
    //    AvgTrackVsNVertex->Draw();
    //    AvgTrackVsNVertex->Fit("pol1","","",0,20);
    //    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ToRun180252/40_25AvgTrackVsNVertex.png");
    //   fout.cd();AvgTrackVsNVertex->Write("AvgTrackVsNVertex");fin.cd();
    //    AvgTrackVsRho->Draw();
    //    AvgTrackVsRho->Fit("pol1","","",0,20);
    //    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ToRun180252/40_25AvgTrackVsRho.png");
    //    fout.cd();AvgTrackVsRho->Write("AvgTrackVsRho");fin.cd();
    //    
    
    AvgEcalIsoVsNVertex_JetReq->Draw();
    AvgEcalIsoVsNVertex_JetReq->Fit("pol1","","",0,20);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ToRun180252/40_25AvgEcalIsoVsNVertex_JetReq.png");
    fout.cd();AvgEcalIsoVsNVertex_JetReq->Write("AvgEcalIsoVsNVertex_JetReq");fin.cd();
    AvgEcalIsoVsRho_JetReq->Draw();
    AvgEcalIsoVsRho_JetReq->Fit("pol1","","",0,20);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ToRun180252/40_25AvgEcalIsoVsRho_JetReq.png");
    fout.cd();AvgEcalIsoVsRho_JetReq->Write("AvgEcalIsoVsRho_JetReq");fin.cd();
    
    AvgHcalIsoVsNVertex_JetReq->Draw();
    AvgHcalIsoVsNVertex_JetReq->Fit("pol1","","",0,20);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ToRun180252/40_25AvgHcalIsoVsNVertex_JetReq.png");
    fout.cd();AvgHcalIsoVsNVertex_JetReq->Write("AvgHcalIsoVsNVertex_JetReq");fin.cd();
    AvgHcalIsoVsRho_JetReq->Draw();
    AvgHcalIsoVsRho_JetReq->Fit("pol1","","",0,20);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ToRun180252/40_25AvgHcalIsoVsRho_JetReq.png");
    fout.cd();AvgHcalIsoVsRho_JetReq->Write("AvgHcalIsoVsRho_JetReq");fin.cd();
    
    AvgHOverEVsNVertex_JetReq->Draw();
    AvgHOverEVsNVertex_JetReq->Fit("pol1","","",0,20);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ToRun180252/40_25AvgHOverEVsNVertex_JetReq.png");
    fout.cd();AvgHOverEVsNVertex_JetReq->Write("AvgHOverEVsNVertex_JetReq");fin.cd();
    AvgHOverEVsRho_JetReq->Draw();
    AvgHOverEVsRho_JetReq->Fit("pol1","","",0,20);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ToRun180252/40_25AvgHOverEVsRho_JetReq.png");
    fout.cd();AvgHOverEVsRho_JetReq->Write("AvgHOverEVsRho_JetReq");fin.cd();
    //  
    //    AvgTrackVsNVertex_JetReq->Draw();
    //    AvgTrackVsNVertex_JetReq->Fit("pol1","","",0,20);
    //    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ToRun180252/40_25AvgTrackVsNVertex_JetReq.png");
    //    fout.cd();AvgTrackVsNVertex_JetReq->Write("AvgTrackVsNVertex_JetReq");fin.cd();
    //    AvgTrackVsRho_JetReq->Draw();
    //    AvgTrackVsRho_JetReq->Fit("pol1","","",0,20);
    //    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ToRun180252/40_25AvgTrackVsRho_JetReq.png");
    //    fout.cd();AvgTrackVsRho_JetReq->Write("AvgTrackVsRho_JetReq");fin.cd();
    //  
    gStyle->SetOptStat(1);
    
    
    //photon from neutralino
    AvgEcalIsoVsNVertex_gFromN->Draw();
    AvgEcalIsoVsNVertex_gFromN->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsNVertex_gFromN.png");
    fout.cd();AvgEcalIsoVsNVertex_gFromN->Write("AvgEcalIsoVsNVertex_gFromN");fin.cd();
    AvgEcalIsoVsRho_gFromN->Draw();
    AvgEcalIsoVsRho_gFromN->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsRho_gFromN.png");
    fout.cd();AvgEcalIsoVsRho_gFromN->Write("AvgEcalIsoVsRho_gFromN");fin.cd();
    AvgHcalIsoVsNVertex_gFromN->Draw();
    AvgHcalIsoVsNVertex_gFromN->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsNVertex_gFromN.png");
    fout.cd();AvgHcalIsoVsNVertex_gFromN->Write("AvgHcalIsoVsNVertex_gFromN");fin.cd();
    AvgHcalIsoVsRho_gFromN->Draw();
    AvgHcalIsoVsRho_gFromN->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsRho_gFromN.png");
    fout.cd();AvgHcalIsoVsRho_gFromN->Write("AvgHcalIsoVsRho_gFromN");fin.cd();
    AvgEcalIsoVsNVertex_gFromN_JetReq->Draw();
    AvgEcalIsoVsNVertex_gFromN_JetReq->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsNVertex_gFromN_JetReq.png");
    fout.cd();AvgEcalIsoVsNVertex_gFromN_JetReq->Write("AvgEcalIsoVsNVertex_gFromN_JetReq");fin.cd();
    AvgEcalIsoVsRho_gFromN_JetReq->Draw();
    AvgEcalIsoVsRho_gFromN_JetReq->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsRho_gFromN_JetReq.png");
    fout.cd();AvgEcalIsoVsRho_gFromN->Write("AvgEcalIsoVsRho_gFromN_JetReq");fin.cd();
    AvgHcalIsoVsNVertex_gFromN_JetReq->Draw();
    AvgHcalIsoVsNVertex_gFromN_JetReq->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsNVertex_gFromN_JetReq.png");
    fout.cd();AvgHcalIsoVsNVertex_gFromN->Write("AvgHcalIsoVsNVertex_gFromN_JetReq");fin.cd();
    AvgHcalIsoVsRho_gFromN_JetReq->Draw();
    AvgHcalIsoVsRho_gFromN_JetReq->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsRho_gFromN_JetReq.png");
    fout.cd();AvgHcalIsoVsRho_gFromN_JetReq->Write("AvgHcalIsoVsRho_gFromN_JetReq");fin.cd();
    //photon from NOT neutralino
    AvgEcalIsoVsNVertex_gFromJ->Draw();
    AvgEcalIsoVsNVertex_gFromJ->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsNVertex_gFromJ.png");
    fout.cd();AvgEcalIsoVsNVertex_gFromJ->Write("AvgEcalIsoVsNVertex_gFromJ");fin.cd();
    AvgEcalIsoVsRho_gFromJ->Draw();
    AvgEcalIsoVsRho_gFromJ->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsRho_gFromJ.png");
    fout.cd();AvgEcalIsoVsRho_gFromJ->Write("AvgEcalIsoVsRho_gFromJ");fin.cd();
    AvgHcalIsoVsNVertex_gFromJ->Draw();
    AvgHcalIsoVsNVertex_gFromJ->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsNVertex_gFromJ.png");
    fout.cd();AvgHcalIsoVsNVertex_gFromJ->Write("AvgHcalIsoVsNVertex_gFromJ");fin.cd();
    AvgHcalIsoVsRho_gFromJ->Draw();
    AvgHcalIsoVsRho_gFromJ->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsRho_gFromJ.png");
    fout.cd();AvgHcalIsoVsRho_gFromJ->Write("AvgHcalIsoVsRho_gFromJ");fin.cd();
    AvgEcalIsoVsNVertex_gFromJ_JetReq->Draw();
    AvgEcalIsoVsNVertex_gFromJ_JetReq->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsNVertex_gFromJ_JetReq.png");
    fout.cd();AvgEcalIsoVsNVertex_gFromJ_JetReq->Write("AvgEcalIsoVsNVertex_gFromJ_JetReq");fin.cd();
    AvgEcalIsoVsRho_gFromJ_JetReq->Draw();
    AvgEcalIsoVsRho_gFromJ_JetReq->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsRho_gFromJ_JetReq.png");
    fout.cd();AvgEcalIsoVsRho_gFromJ_JetReq->Write("AvgEcalIsoVsRho_gFromJ_JetReq");fin.cd();
    AvgHcalIsoVsNVertex_gFromJ_JetReq->Draw();
    AvgHcalIsoVsNVertex_gFromJ_JetReq->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsNVertex_gFromJ_JetReq.png");
    fout.cd();AvgHcalIsoVsNVertex_gFromJ_JetReq->Write("AvgHcalIsoVsNVertex_gFromJ_JetReq");fin.cd();
    AvgHcalIsoVsRho_gFromJ_JetReq->Draw();
    AvgHcalIsoVsRho_gFromJ_JetReq->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsRho_gFromJ_JetReq.png");
    fout.cd();AvgHcalIsoVsRho_gFromJ_JetReq->Write("AvgHcalIsoVsRho_gFromJ_JetReq");fin.cd();
    //Fail SigmaIetaIeta or TrackIso
    AvgEcalIsoVsNVertex_FailSigIetaOrTrackIso->Draw();
    AvgEcalIsoVsNVertex_FailSigIetaOrTrackIso->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsNVertex_FailSigIetaOrTrackIso.png");
    fout.cd();AvgEcalIsoVsNVertex_FailSigIetaOrTrackIso->Write("AvgEcalIsoVsNVertex_FailSigIetaOrTrackIso");fin.cd();
    AvgEcalIsoVsRho_FailSigIetaOrTrackIso->Draw();
    AvgEcalIsoVsRho_FailSigIetaOrTrackIso->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsRho_FailSigIetaOrTrackIso.png");
    fout.cd();AvgEcalIsoVsRho_FailSigIetaOrTrackIso->Write("AvgEcalIsoVsRho_FailSigIetaOrTrackIso");fin.cd();
    AvgHcalIsoVsNVertex_FailSigIetaOrTrackIso->Draw();
    AvgHcalIsoVsNVertex_FailSigIetaOrTrackIso->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsNVertex_FailSigIetaOrTrackIso.png");
    fout.cd();AvgHcalIsoVsNVertex_FailSigIetaOrTrackIso->Write("AvgHcalIsoVsNVertex_FailSigIetaOrTrackIso");fin.cd();
    AvgHcalIsoVsRho_FailSigIetaOrTrackIso->Draw();
    AvgHcalIsoVsRho_FailSigIetaOrTrackIso->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsRho_FailSigIetaOrTrackIso.png");
    fout.cd();AvgHcalIsoVsRho_FailSigIetaOrTrackIso->Write("AvgHcalIsoVsRho_FailSigIetaOrTrackIso");fin.cd();
    AvgEcalIsoVsNVertex_JetReq_FailSigIetaOrTrackIso->Draw();
    AvgEcalIsoVsNVertex_JetReq_FailSigIetaOrTrackIso->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsNVertex_JetReq_FailSigIetaOrTrackIso.png");
    fout.cd();AvgEcalIsoVsNVertex_JetReq_FailSigIetaOrTrackIso->Write("AvgEcalIsoVsNVertex_JetReq_FailSigIetaOrTrackIso");fin.cd();
    AvgEcalIsoVsRho_JetReq_FailSigIetaOrTrackIso->Draw();
    AvgEcalIsoVsRho_JetReq_FailSigIetaOrTrackIso->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsRho_JetReq_FailSigIetaOrTrackIso.png");
    fout.cd();AvgEcalIsoVsRho_JetReq_FailSigIetaOrTrackIso->Write("AvgEcalIsoVsRho_JetReq_FailSigIetaOrTrackIso");fin.cd();
    AvgHcalIsoVsNVertex_JetReq_FailSigIetaOrTrackIso->Draw();
    AvgHcalIsoVsNVertex_JetReq_FailSigIetaOrTrackIso->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsNVertex_JetReq_FailSigIetaOrTrackIso.png");
    fout.cd();AvgHcalIsoVsNVertex_JetReq_FailSigIetaOrTrackIso->Write("AvgHcalIsoVsNVertex_JetReq_FailSigIetaOrTrackIso");fin.cd();
    AvgHcalIsoVsRho_JetReq_FailSigIetaOrTrackIso->Draw();
    AvgHcalIsoVsRho_JetReq_FailSigIetaOrTrackIso->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsRho_JetReq_FailSigIetaOrTrackIso.png");
    fout.cd();AvgHcalIsoVsRho_JetReq_FailSigIetaOrTrackIso->Write("AvgHcalIsoVsRho_JetReq_FailSigIetaOrTrackIso");fin.cd();
    //Signal Cuts, No pixel seed
    AvgEcalIsoVsNVertex_SignalCutsNoPixel->Draw();
    AvgEcalIsoVsNVertex_SignalCutsNoPixel->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsNVertex_SignalCutsNoPixel.png");
    fout.cd();AvgEcalIsoVsNVertex_SignalCutsNoPixel->Write("AvgEcalIsoVsNVertex_SignalCutsNoPixel");fin.cd();
    AvgEcalIsoVsRho_SignalCutsNoPixel->Draw();
    AvgEcalIsoVsRho_SignalCutsNoPixel->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsRho_SignalCutsNoPixel.png");
    fout.cd();AvgEcalIsoVsRho_SignalCutsNoPixel->Write("AvgEcalIsoVsRho_SignalCutsNoPixel");fin.cd();
    AvgHcalIsoVsNVertex_SignalCutsNoPixel->Draw();
    AvgHcalIsoVsNVertex_SignalCutsNoPixel->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsNVertex_SignalCutsNoPixel.png");
    fout.cd();AvgHcalIsoVsNVertex_SignalCutsNoPixel->Write("AvgHcalIsoVsNVertex_SignalCutsNoPixel");fin.cd();
    AvgHcalIsoVsRho_SignalCutsNoPixel->Draw();
    AvgHcalIsoVsRho_SignalCutsNoPixel->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsRho_SignalCutsNoPixel.png");
    fout.cd();AvgHcalIsoVsRho_SignalCutsNoPixel->Write("AvgHcalIsoVsRho_SignalCutsNoPixel");fin.cd();
    AvgEcalIsoVsNVertex_JetReq_SignalCutsNoPixel->Draw();
    AvgEcalIsoVsNVertex_JetReq_SignalCutsNoPixel->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsNVertex_JetReq_SignalCutsNoPixel.png");
    fout.cd();AvgEcalIsoVsNVertex_JetReq_SignalCutsNoPixel->Write("AvgEcalIsoVsNVertex_JetReq_SignalCutsNoPixel");fin.cd();
    AvgEcalIsoVsRho_JetReq_SignalCutsNoPixel->Draw();
    AvgEcalIsoVsRho_JetReq_SignalCutsNoPixel->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsRho_JetReq_SignalCutsNoPixel.png");
    fout.cd();AvgEcalIsoVsRho_JetReq_SignalCutsNoPixel->Write("AvgEcalIsoVsRho_JetReq_SignalCutsNoPixel");fin.cd();
    AvgHcalIsoVsNVertex_JetReq_SignalCutsNoPixel->Draw();
    AvgHcalIsoVsNVertex_JetReq_SignalCutsNoPixel->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsNVertex_JetReq_SignalCutsNoPixel.png");
    fout.cd();AvgHcalIsoVsNVertex_JetReq_SignalCutsNoPixel->Write("AvgHcalIsoVsNVertex_JetReq_SignalCutsNoPixel");fin.cd();
    AvgHcalIsoVsRho_JetReq_SignalCutsNoPixel->Draw();
    AvgHcalIsoVsRho_JetReq_SignalCutsNoPixel->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsRho_JetReq_SignalCutsNoPixel.png");
    fout.cd();AvgHcalIsoVsRho_JetReq_SignalCutsNoPixel->Write("AvgHcalIsoVsRho_JetReq_SignalCutsNoPixel");fin.cd();
    //Signal Cuts, With pixel seed
    AvgEcalIsoVsNVertex_SignalCutsPixel->Draw();
    AvgEcalIsoVsNVertex_SignalCutsPixel->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsNVertex_SignalCutsPixel.png");
    fout.cd();AvgEcalIsoVsNVertex_SignalCutsPixel->Write("AvgEcalIsoVsNVertex_SignalCutsPixel");fin.cd();
    AvgEcalIsoVsRho_SignalCutsPixel->Draw();
    AvgEcalIsoVsRho_SignalCutsPixel->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsRho_SignalCutsPixel.png");
    fout.cd();AvgEcalIsoVsRho_SignalCutsPixel->Write("AvgEcalIsoVsRho_SignalCutsPixel");fin.cd();
    AvgHcalIsoVsNVertex_SignalCutsPixel->Draw();
    AvgHcalIsoVsNVertex_SignalCutsPixel->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsNVertex_SignalCutsPixel.png");
    fout.cd();AvgHcalIsoVsNVertex_SignalCutsPixel->Write("AvgHcalIsoVsNVertex_SignalCutsPixel");fin.cd();
    AvgHcalIsoVsRho_SignalCutsPixel->Draw();
    AvgHcalIsoVsRho_SignalCutsPixel->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsRho_SignalCutsPixel.png");
    fout.cd();AvgHcalIsoVsRho_SignalCutsPixel->Write("AvgHcalIsoVsRho_SignalCutsPixel");fin.cd();
    AvgEcalIsoVsNVertex_JetReq_SignalCutsPixel->Draw();
    AvgEcalIsoVsNVertex_JetReq_SignalCutsPixel->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsNVertex_JetReq_SignalCutsPixel.png");
    fout.cd();AvgEcalIsoVsNVertex_JetReq_SignalCutsPixel->Write("AvgEcalIsoVsNVertex_JetReq_SignalCutsPixel");fin.cd();
    AvgEcalIsoVsRho_JetReq_SignalCutsPixel->Draw();
    AvgEcalIsoVsRho_JetReq_SignalCutsPixel->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsRho_JetReq_SignalCutsPixel.png");
    fout.cd();AvgEcalIsoVsRho_JetReq_SignalCutsPixel->Write("AvgEcalIsoVsRho_JetReq_SignalCutsPixel");fin.cd();
    AvgHcalIsoVsNVertex_JetReq_SignalCutsPixel->Draw();
    AvgHcalIsoVsNVertex_JetReq_SignalCutsPixel->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgHcalIsoVsNVertex_JetReq_SignalCutsPixel.png");
    fout.cd();AvgHcalIsoVsNVertex_JetReq_SignalCutsPixel->Write("AvgHcalIsoVsNVertex_JetReq_SignalCutsPixel");fin.cd();
    AvgHcalIsoVsRho_JetReq_SignalCutsPixel->Draw();
    AvgHcalIsoVsRho_JetReq_SignalCutsPixel->Fit("pol1","","",1,16);
    c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/AvgEcalIsoVsRho_JetReq_SignalCutsPixel.png");
    fout.cd();AvgHcalIsoVsRho_JetReq_SignalCutsPixel->Write("AvgHcalIsoVsRho_JetReq_SignalCutsPixel");fin.cd();*/
  /*
  //-------Background Efficiencies for Pileup study-------
  TGraphAsymmErrors *BgroundEffVsNVertex = new TGraphAsymmErrors(NumerVsNVertex,DenomVsNVertex_afterHE,"");
  BgroundEffVsNVertex->SetMarkerSize(0.6);
  BgroundEffVsNVertex->GetXaxis()->SetTitle("NVertex");
  //BgroundEffVsNVertex->GetYaxis()->SetTitleOffset(.8);
  BgroundEffVsNVertex->GetYaxis()->SetTitle("Efficiency");
  //BgroundEffVsNVertex->Write("BgroundEffVsNVertex");
  TGraphAsymmErrors *RhoCorrectedBgroundEffVsNVertex = new TGraphAsymmErrors(RhoCorrectedNumerVsNVertex,DenomVsNVertex_afterHE,"");
  RhoCorrectedBgroundEffVsNVertex->SetMarkerSize(0.6);
  RhoCorrectedBgroundEffVsNVertex->GetXaxis()->SetTitle("NVertex");
  //RhoCorrectedBgroundEffVsNVertex->GetYaxis()->SetTitleOffset(.8);
  RhoCorrectedBgroundEffVsNVertex->GetYaxis()->SetTitle("Efficiency");
  RhoCorrectedBgroundEffVsNVertex->SetLineColor(kRed);
  RhoCorrectedBgroundEffVsNVertex->SetMarkerColor(kRed);
  //RhoCorrectedBgroundEffVsNVertex->Write("RhoCorrectedBgroundEffVsNVertex");
  TGraphAsymmErrors *NVertexCorrectedBgroundEffVsNVertex = new TGraphAsymmErrors(NVertexCorrectedNumerVsNVertex,DenomVsNVertex_afterHE,"");
  NVertexCorrectedBgroundEffVsNVertex->SetMarkerSize(0.6);
  NVertexCorrectedBgroundEffVsNVertex->GetXaxis()->SetTitle("NVertex");
  //NVertexCorrectedBgroundEffVsNVertex->GetYaxis()->SetTitleOffset(.8);
  NVertexCorrectedBgroundEffVsNVertex->GetYaxis()->SetTitle("Efficiency");
  NVertexCorrectedBgroundEffVsNVertex->SetLineColor(kBlue);
  NVertexCorrectedBgroundEffVsNVertex->SetMarkerColor(kBlue);
  //NVertexCorrectedBgroundEffVsNVertex->Write("NVertexCorrectedBgroundEffVsNVertex");
  BgroundEffVsNVertex->SetMaximum(0.04);
  BgroundEffVsNVertex->Draw("AP");
  RhoCorrectedBgroundEffVsNVertex->Draw("PSAME");
  NVertexCorrectedBgroundEffVsNVertex->Draw("PSAME");
  TLegend *leg= new TLegend(.3,.5,.5,.7,"","brNDC");
  leg->SetTextSize(0.04);
  leg->AddEntry(BgroundEffVsNVertex,"No Correction");
  leg->AddEntry(RhoCorrectedBgroundEffVsNVertex,"Rho Correction");
  leg->AddEntry(NVertexCorrectedBgroundEffVsNVertex,"NVertex Correction");
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/BgroundEffVsNVertex.png");

  NumerVsEt->Rebin(2);DenomVsEt_afterHE->Rebin(2);
  TGraphAsymmErrors *BgroundEffVsEt = new TGraphAsymmErrors(NumerVsEt,DenomVsEt_afterHE,"");
  BgroundEffVsEt->SetMarkerSize(0.6);
  BgroundEffVsEt->GetXaxis()->SetTitle("Et");
  //BgroundEffVsEt->GetYaxis()->SetTitleOffset(.8);
  BgroundEffVsEt->GetYaxis()->SetTitle("Efficiency");
  //BgroundEffVsEt->Write("BgroundEffVsEt");
  RhoCorrectedNumerVsEt->Rebin(2);
  TGraphAsymmErrors *RhoCorrectedBgroundEffVsEt = new TGraphAsymmErrors(RhoCorrectedNumerVsEt,DenomVsEt_afterHE,"");
  RhoCorrectedBgroundEffVsEt->SetMarkerSize(0.6);
  RhoCorrectedBgroundEffVsEt->GetXaxis()->SetTitle("Et");
  //RhoCorrectedBgroundEffVsEt->GetYaxis()->SetTitleOffset(.8);
  RhoCorrectedBgroundEffVsEt->GetYaxis()->SetTitle("Efficiency");
  RhoCorrectedBgroundEffVsEt->SetLineColor(kRed);
  RhoCorrectedBgroundEffVsEt->SetMarkerColor(kRed);
  //RhoCorrectedBgroundEffVsEt->Write("RhoCorrectedBgroundEffVsEt");
  NVertexCorrectedNumerVsEt->Rebin(2);
  TGraphAsymmErrors *NVertexCorrectedBgroundEffVsEt = new TGraphAsymmErrors(NVertexCorrectedNumerVsEt,DenomVsEt_afterHE,"");
  NVertexCorrectedBgroundEffVsEt->SetMarkerSize(0.6);
  NVertexCorrectedBgroundEffVsEt->GetXaxis()->SetTitle("Et");
  //NVertexCorrectedBgroundEffVsEt->GetYaxis()->SetTitleOffset(.8);
  NVertexCorrectedBgroundEffVsEt->GetYaxis()->SetTitle("Efficiency");
  NVertexCorrectedBgroundEffVsEt->SetLineColor(kBlue);
  NVertexCorrectedBgroundEffVsEt->SetMarkerColor(kBlue);
  //NVertexCorrectedBgroundEffVsEt->Write("NVertexCorrectedBgroundEffVsEt");
  BgroundEffVsEt->SetMinimum(0);
  BgroundEffVsEt->SetMaximum(0.04);
  BgroundEffVsEt->GetXaxis()->SetLimits(0,300);
  BgroundEffVsEt->Draw("AP");
  RhoCorrectedBgroundEffVsEt->Draw("PSAME");
  NVertexCorrectedBgroundEffVsEt->Draw("PSAME");
  TLegend *leg= new TLegend(.25,.6,.45,.8,"","brNDC");
  leg->SetTextSize(0.04);
  leg->AddEntry(BgroundEffVsEt,"No Correction");
  leg->AddEntry(RhoCorrectedBgroundEffVsEt,"Rho Correction");
  leg->AddEntry(NVertexCorrectedBgroundEffVsEt,"NVertex Correction");
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/BgroundEffVsEt.png");

  TGraphAsymmErrors *BgroundEffVsEta = new TGraphAsymmErrors(NumerVsEta,DenomVsEta_afterHE,"");
  BgroundEffVsEta->SetMarkerSize(0.6);
  BgroundEffVsEta->GetXaxis()->SetTitle("Eta");
  //BgroundEffVsEta->GetYaxis()->SetTitleOffset(.8);
  BgroundEffVsEta->GetYaxis()->SetTitle("Efficiency");
  //BgroundEffVsEta->Write("BgroundEffVsEta");
  TGraphAsymmErrors *RhoCorrectedBgroundEffVsEta = new TGraphAsymmErrors(RhoCorrectedNumerVsEta,DenomVsEta_afterHE,"");
  RhoCorrectedBgroundEffVsEta->SetMarkerSize(0.6);
  RhoCorrectedBgroundEffVsEta->GetXaxis()->SetTitle("Eta");
  //RhoCorrectedBgroundEffVsEta->GetYaxis()->SetTitleOffset(.8);
  RhoCorrectedBgroundEffVsEta->GetYaxis()->SetTitle("Efficiency");
  RhoCorrectedBgroundEffVsEta->SetLineColor(kRed);
  RhoCorrectedBgroundEffVsEta->SetMarkerColor(kRed);
  //RhoCorrectedBgroundEffVsEta->Write("RhoCorrectedBgroundEffVsEta");
  TGraphAsymmErrors *NVertexCorrectedBgroundEffVsEta = new TGraphAsymmErrors(NVertexCorrectedNumerVsEta,DenomVsEta_afterHE,"");
  NVertexCorrectedBgroundEffVsEta->SetMarkerSize(0.6);
  NVertexCorrectedBgroundEffVsEta->GetXaxis()->SetTitle("Eta");
  //NVertexCorrectedBgroundEffVsEta->GetYaxis()->SetTitleOffset(.8);
  NVertexCorrectedBgroundEffVsEta->GetYaxis()->SetTitle("Efficiency");
  NVertexCorrectedBgroundEffVsEta->SetLineColor(kBlue);
  NVertexCorrectedBgroundEffVsEta->SetMarkerColor(kBlue);
  //NVertexCorrectedBgroundEffVsEta->Write("NVertexCorrectedBgroundEffVsEta");
  BgroundEffVsEta->SetMinimum(0.015);
  BgroundEffVsEta->SetMaximum(.02);
  BgroundEffVsEta->Draw("AP");
  RhoCorrectedBgroundEffVsEta->Draw("PSAME");
  NVertexCorrectedBgroundEffVsEta->Draw("PSAME");
  TLegend *leg= new TLegend(.45,.2,.6,.4,"","brNDC");
  leg->SetTextSize(0.04);
  leg->AddEntry(BgroundEffVsEta,"No Correction");
  leg->AddEntry(RhoCorrectedBgroundEffVsEta,"Rho Correction");
  leg->AddEntry(NVertexCorrectedBgroundEffVsEta,"NVertex Correction");
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/BgroundEffVsEta.png");

  TGraphAsymmErrors *BgroundEffVsMET = new TGraphAsymmErrors(NumerVsMET,DenomVsMET_afterHE,"");
  BgroundEffVsMET->SetMarkerSize(0.6);
  BgroundEffVsMET->GetXaxis()->SetTitle("MET");
  //BgroundEffVsMET->GetYaxis()->SetTitleOffset(.8);
  BgroundEffVsMET->GetYaxis()->SetTitle("Efficiency");
  //BgroundEffVsMET->Write("BgroundEffVsMET");
  TGraphAsymmErrors *RhoCorrectedBgroundEffVsMET = new TGraphAsymmErrors(RhoCorrectedNumerVsMET,DenomVsMET_afterHE,"");
  RhoCorrectedBgroundEffVsMET->SetMarkerSize(0.6);
  RhoCorrectedBgroundEffVsMET->GetXaxis()->SetTitle("MET");
  //RhoCorrectedBgroundEffVsMET->GetYaxis()->SetTitleOffset(.8);
  RhoCorrectedBgroundEffVsMET->GetYaxis()->SetTitle("Efficiency");
  RhoCorrectedBgroundEffVsMET->SetLineColor(kRed);
  RhoCorrectedBgroundEffVsMET->SetMarkerColor(kRed);
  //RhoCorrectedBgroundEffVsMET->Write("RhoCorrectedBgroundEffVsMET");
  TGraphAsymmErrors *NVertexCorrectedBgroundEffVsMET = new TGraphAsymmErrors(NVertexCorrectedNumerVsMET,DenomVsMET_afterHE,"");
  NVertexCorrectedBgroundEffVsMET->SetMarkerSize(0.6);
  NVertexCorrectedBgroundEffVsMET->GetXaxis()->SetTitle("MET");
  //NVertexCorrectedBgroundEffVsMET->GetYaxis()->SetTitleOffset(.8);
  NVertexCorrectedBgroundEffVsMET->GetYaxis()->SetTitle("Efficiency");
  NVertexCorrectedBgroundEffVsMET->SetLineColor(kBlue);
  NVertexCorrectedBgroundEffVsMET->SetMarkerColor(kBlue);
  BgroundEffVsMET->GetXaxis()->SetLimits(0,200);
  BgroundEffVsMET->SetMaximum(0.2);
  BgroundEffVsMET->Draw("AP");
  RhoCorrectedBgroundEffVsMET->Draw("PSAME");
  NVertexCorrectedBgroundEffVsMET->Draw("PSAME");
  TLegend *leg= new TLegend(.3,.5,.5,.7,"","brNDC");
  leg->SetTextSize(0.04);
  leg->AddEntry(BgroundEffVsMET,"No Correction");
  leg->AddEntry(RhoCorrectedBgroundEffVsMET,"Rho Correction");
  leg->AddEntry(NVertexCorrectedBgroundEffVsMET,"NVertex Correction");
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/BgroundEffVsMET.png");
  // end comment
  }
  */
  
  //-------------end of linear scale plots
  float fakeRate[6]={0.,0.,0.,0.,0.,0.},fakeRateErr[6]={0.,0.,0.,0.,0.,0.},fakeRateFull=0.,fakeRateErrFull=0.;

  cout<<"ee Invariant Mass Fit:"<<endl<<" signal events: "<<CBsigEE<<" +- "<<CBsigEEerror<<endl;
  cout<<"eg Invariant Mass Fit:"<<endl<<" signal events: "<<CBsigEG<<" +- "<<CBsigEGerror<<endl;
  cout<<"e->g fake rate f(e->g) = "<<FakeRate<< " +- "<< FakeRateErrForLimits<<" (stat.) +- "<<FakeRateSyst<<" (syst.)"<<endl;
  //Pt25to40
  cout<<"|---------Pt25to40---------|"<<endl;
  cout<<"ee Invariant Mass Fit Pt25to40:"<<endl<<" signal events: "<<CBsigEE25to40<<" +- "<<CBsigEEerror25to40<<endl;
  cout<<"eg Invariant Mass Fit Pt25to40:"<<endl<<" signal events: "<<CBsigEG25to40<<" +- "<<CBsigEGerror25to40<<endl;
  GetFakeRateAndError(CBsigEE25to40,CBsigEEerror25to40,CBsigEG25to40,CBsigEGerror25to40,fakeRate[0],fakeRateErr[0]);
  cout<<"e->g fake rate f(e->g) = "<<fakeRate[0]<< " +- "<<fakeRateErr[0]<<endl;
  //Pt40to45
  cout<<"|---------Pt40to45---------|"<<endl;
  cout<<"ee Invariant Mass Fit Pt40to45:"<<endl<<" signal events: "<<CBsigEE40to45<<" +- "<<CBsigEEerror40to45<<endl;
  cout<<"eg Invariant Mass Fit Pt40to45:"<<endl<<" signal events: "<<CBsigEG40to45<<" +- "<<CBsigEGerror40to45<<endl;
  GetFakeRateAndError(CBsigEE40to45,CBsigEEerror40to45,CBsigEG40to45,CBsigEGerror40to45,fakeRate[1],fakeRateErr[1]);
  cout<<"e->g fake rate f(e->g) = "<<fakeRate[1]<< " +- "<<fakeRateErr[1]<<endl;
  //Pt45to50
  cout<<"|---------Pt45to50---------|"<<endl;
  cout<<"ee Invariant Mass Fit Pt45to50:"<<endl<<" signal events: "<<CBsigEE45to50<<" +- "<<CBsigEEerror45to50<<endl;
  cout<<"eg Invariant Mass Fit Pt45to50:"<<endl<<" signal events: "<<CBsigEG45to50<<" +- "<<CBsigEGerror45to50<<endl;
  GetFakeRateAndError(CBsigEE45to50,CBsigEEerror45to50,CBsigEG45to50,CBsigEGerror45to50,fakeRate[2],fakeRateErr[2]);
  cout<<"e->g fake rate f(e->g) = "<<fakeRate[2]<< " +- "<<fakeRateErr[2]<<endl;
  //Pt50to60
  cout<<"|---------Pt50to60---------|"<<endl;
  cout<<"ee Invariant Mass Fit Pt50to60:"<<endl<<" signal events: "<<CBsigEE50to60<<" +- "<<CBsigEEerror50to60<<endl;
  cout<<"eg Invariant Mass Fit Pt50to60:"<<endl<<" signal events: "<<CBsigEG50to60<<" +- "<<CBsigEGerror50to60<<endl;
  GetFakeRateAndError(CBsigEE50to60,CBsigEEerror50to60,CBsigEG50to60,CBsigEGerror50to60,fakeRate[3],fakeRateErr[3]);
  cout<<"e->g fake rate f(e->g) = "<<fakeRate[3]<< " +- "<<fakeRateErr[3]<<endl;
  //Pt60to80
  cout<<"|---------Pt60to80---------|"<<endl;
  cout<<"ee Invariant Mass Fit Pt60to80:"<<endl<<" signal events: "<<CBsigEE60to80<<" +- "<<CBsigEEerror60to80<<endl;
  cout<<"eg Invariant Mass Fit Pt60to80:"<<endl<<" signal events: "<<CBsigEG60to80<<" +- "<<CBsigEGerror60to80<<endl;
  GetFakeRateAndError(CBsigEE60to80,CBsigEEerror60to80,CBsigEG60to80,CBsigEGerror60to80,fakeRate[4],fakeRateErr[4]);
  cout<<"e->g fake rate f(e->g) = "<<fakeRate[4]<< " +- "<<fakeRateErr[4]<<endl;
  //Pt80
  cout<<"|-----------Pt80-----------|"<<endl;
  cout<<"ee Invariant Mass Fit Pt80:"<<endl<<" signal events: "<<CBsigEE80<<" +- "<<CBsigEEerror80<<endl;
  cout<<"eg Invariant Mass Fit Pt80:"<<endl<<" signal events: "<<CBsigEG80<<" +- "<<CBsigEGerror80<<endl;
  GetFakeRateAndError(CBsigEE80,CBsigEEerror80,CBsigEG80,CBsigEGerror80,fakeRate[5],fakeRateErr[5]);
  cout<<"e->g fake rate f(e->g) = "<<fakeRate[5]<< " +- "<<fakeRateErr[5]<<endl;
  cout<<"|-----------Now using gg-----------|"<<endl;


  cout<<"gg Invariant Mass Fit:"<<endl<<" signal events: "<<CBsigGG<<endl;
  fakeRateFull=(-CBsigGG+sqrt(CBsigEE*CBsigGG))/(CBsigEE-CBsigGG);
  fakeRateErrFull=( (sqrt(CBsigEE*CBsigGG)-CBsigGG)/(CBsigEE-CBsigGG) ) *
    ( ( ( (sqrt(CBsigGG*CBsigGG*CBsigEEerror*CBsigEEerror + CBsigGGerror*CBsigGGerror*CBsigEE*CBsigEE)/(2*sqrt(CBsigGG*CBsigEE)) + CBsigGGerror) *
	  (sqrt(CBsigGG*CBsigGG*CBsigEEerror*CBsigEEerror + CBsigGGerror*CBsigGGerror*CBsigEE*CBsigEE)/(2*sqrt(CBsigGG*CBsigEE)) + CBsigGGerror) ) / 
	( (sqrt(CBsigEE-CBsigGG)-CBsigGG)*(sqrt(CBsigEE-CBsigGG)-CBsigGG) ) ) +
      ( ((CBsigEEerror+CBsigGGerror)*(CBsigEEerror+CBsigGGerror)) / ( (CBsigEE-CBsigGG)*(CBsigEE-CBsigGG) ) ) );
  cout<<"e->g fake rate f(e->g) = "<<fakeRateFull<<" +- "<<fakeRateErrFull<<endl;

  float misIdBins[8]={0,25,40,45,50,60,80,200};
  TH1F* PtDepFakeRate = new TH1F("PtDepFakeRate","",7,misIdBins);
  for(int i=1;i<7;i++){
    PtDepFakeRate->SetBinContent(i+1,fakeRate[i-1]);
    PtDepFakeRate->SetBinError(i+1,fakeRateErr[i-1]);
  }
  //gStyle->SetErrorX(0);
  c1->cd();
  PtDepFakeRate->GetYaxis()->SetTitle("Electron Fake Rate");
  PtDepFakeRate->GetXaxis()->SetTitle("p_{T} (GeV)");
  //PtDepFakeRate->GetYaxis()->SetRangeUser(0.0,0.035);
  //PtDepFakeRate->Fit("pol0","","",40,200);
  PtDepFakeRate->Draw();
  ProtoText->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure13_PtDepFakeRate.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure13_PtDepFakeRate.pdf");
  fout.cd();PtDepFakeRate->Write("PtDepFakeRate");fin.cd();

  fin.cd();
  c1->cd();
  c1->SetLogy(1);
  TH1F* gg010 = (TH1F*)fin.Get("ggMet_NV0_10");
  TH1F* eg010 = (TH1F*)fin.Get("egMet_NV0_10");
  TH1F* ee010 = (TH1F*)fin.Get("eeMet_NV0_10");
  TH1F* eeSB010 = (TH1F*)fin.Get("eeSBMet_NV0_10");
  TH1F* ff010 = (TH1F*)fin.Get("ffMet_NV0_10");
  TH1F* gg1015 = (TH1F*)fin.Get("ggMet_NV10_15");
  TH1F* eg1015 = (TH1F*)fin.Get("egMet_NV10_15");
  TH1F* ee1015 = (TH1F*)fin.Get("eeMet_NV10_15");
  TH1F* eeSB1015 = (TH1F*)fin.Get("eeSBMet_NV10_15");
  TH1F* ff1015 = (TH1F*)fin.Get("ffMet_NV10_15");
  TH1F* gg15up = (TH1F*)fin.Get("ggMet_NV15up");
  TH1F* eg15up = (TH1F*)fin.Get("egMet_NV15up");
  TH1F* ee15up = (TH1F*)fin.Get("eeMet_NV15up");
  TH1F* eeSB15up = (TH1F*)fin.Get("eeSBMet_NV15up");
  TH1F* ff15up = (TH1F*)fin.Get("ffMet_NV15up");

  gg010->Sumw2();ee010->Sumw2();eeSB010->Sumw2();ff010->Sumw2();
  gg1015->Sumw2();ee1015->Sumw2();eeSB1015->Sumw2();ff1015->Sumw2();
  gg15up->Sumw2();ee15up->Sumw2();eeSB15up->Sumw2();ff15up->Sumw2();

  TH1F* gg010n = (TH1F*)gg010->Rebin(NmetBins,"gg010n",xbins);
  TH1F* eg010n = (TH1F*)eg010->Rebin(NmetBins,"eg010n",xbins);
  TH1F* ee010n = (TH1F*)ee010->Rebin(NmetBins,"ee010n",xbins);
  TH1F* eeSB010n = (TH1F*)eeSB010->Rebin(NmetBins,"eeSB010n",xbins);
  TH1F* ff010n = (TH1F*)ff010->Rebin(NmetBins,"ff010n",xbins);
  TH1F* gg1015n = (TH1F*)gg1015->Rebin(NmetBins,"gg1015n",xbins);
  TH1F* eg1015n = (TH1F*)eg1015->Rebin(NmetBins,"eg1015n",xbins);
  TH1F* ee1015n = (TH1F*)ee1015->Rebin(NmetBins,"ee1015n",xbins);
  TH1F* eeSB1015n = (TH1F*)eeSB1015->Rebin(NmetBins,"eeSB1015n",xbins);
  TH1F* ff1015n = (TH1F*)ff1015->Rebin(NmetBins,"ff1015n",xbins);
  TH1F* gg15upn = (TH1F*)gg15up->Rebin(NmetBins,"gg15upn",xbins);
  TH1F* eg15upn = (TH1F*)eg15up->Rebin(NmetBins,"eg15upn",xbins);
  TH1F* ee15upn = (TH1F*)ee15up->Rebin(NmetBins,"ee15upn",xbins);
  TH1F* eeSB15upn = (TH1F*)eeSB15up->Rebin(NmetBins,"eeSB15upn",xbins);
  TH1F* ff15upn = (TH1F*)ff15up->Rebin(NmetBins,"ff15upn",xbins);
  
  ee010n->Add(eeSB010n,-1);
  ee1015n->Add(eeSB1015n,-1);
  ee15upn->Add(eeSB15upn,-1);
  /*
    float scale = 1./gg010n->Integral(0,4);
    gg010n->Scale(scale);
    scale = 1./gg1015n->Integral();
    gg1015n->Scale(scale);
    scale = 1./gg15upn->Integral();
    gg15upn->Scale(scale);*/
  scale = gg010n->Integral(0,4)/ee010n->Integral(0,4);
  ee010n->Scale(scale);
  scale = gg1015n->Integral(0,4)/ee1015n->Integral(0,4);
  ee1015n->Scale(scale);
  scale = gg15upn->Integral(0,4)/ee15upn->Integral(0,4);
  ee15upn->Scale(scale);
  scale = gg010n->Integral(0,4)/ff010n->Integral(0,4);
  ff010n->Scale(scale);
  scale = gg1015n->Integral(0,4)/ff1015n->Integral(0,4);
  ff1015n->Scale(scale);
  scale = gg15upn->Integral(0,4)/ff15upn->Integral(0,4);
  ff15upn->Scale(scale);

  for(int i=1;i<gg010n->GetNbinsX()+1;i++){
    gg010n->SetBinContent(i,gg010n->GetBinContent(i)/gg010n->GetBinWidth(i));
    eg010n->SetBinContent(i,eg010n->GetBinContent(i)/eg010n->GetBinWidth(i));
    ee010n->SetBinContent(i,ee010n->GetBinContent(i)/ee010n->GetBinWidth(i));
    ff010n->SetBinContent(i,ff010n->GetBinContent(i)/ff010n->GetBinWidth(i));
    gg1015n->SetBinContent(i,gg1015n->GetBinContent(i)/gg1015n->GetBinWidth(i));
    eg1015n->SetBinContent(i,eg1015n->GetBinContent(i)/eg1015n->GetBinWidth(i));
    ee1015n->SetBinContent(i,ee1015n->GetBinContent(i)/ee1015n->GetBinWidth(i));
    ff1015n->SetBinContent(i,ff1015n->GetBinContent(i)/ff1015n->GetBinWidth(i));
    gg15upn->SetBinContent(i,gg15upn->GetBinContent(i)/gg15upn->GetBinWidth(i));
    eg15upn->SetBinContent(i,eg15upn->GetBinContent(i)/eg15upn->GetBinWidth(i));
    ee15upn->SetBinContent(i,ee15upn->GetBinContent(i)/ee15upn->GetBinWidth(i));
    ff15upn->SetBinContent(i,ff15upn->GetBinContent(i)/ff15upn->GetBinWidth(i));
    gg010n->SetBinError(i,gg010n->GetBinError(i)/gg010n->GetBinWidth(i));
    eg010n->SetBinError(i,eg010n->GetBinError(i)/eg010n->GetBinWidth(i));
    ee010n->SetBinError(i,ee010n->GetBinError(i)/ee010n->GetBinWidth(i));
    ff010n->SetBinError(i,ff010n->GetBinError(i)/ff010n->GetBinWidth(i));
    gg1015n->SetBinError(i,gg1015n->GetBinError(i)/gg1015n->GetBinWidth(i));
    eg1015n->SetBinError(i,eg1015n->GetBinError(i)/eg1015n->GetBinWidth(i));
    ee1015n->SetBinError(i,ee1015n->GetBinError(i)/ee1015n->GetBinWidth(i));
    ff1015n->SetBinError(i,ff1015n->GetBinError(i)/ff1015n->GetBinWidth(i));
    gg15upn->SetBinError(i,gg15upn->GetBinError(i)/gg15upn->GetBinWidth(i));
    eg15upn->SetBinError(i,eg15upn->GetBinError(i)/eg15upn->GetBinWidth(i));
    ee15upn->SetBinError(i,ee15upn->GetBinError(i)/ee15upn->GetBinWidth(i));
    ff15upn->SetBinError(i,ff15upn->GetBinError(i)/ff15upn->GetBinWidth(i));
  }

  float min = 5e-3;//gg1015n->GetMinimum()/2.;
  float max = 2000;//ee010n->GetMaximum()+1000;

  gg010n->SetMarkerSize(.75);
  gg1015n->SetMarkerSize(.75);
  gg15upn->SetMarkerSize(.75);
  ee010n->SetMarkerSize(.75);
  ee1015n->SetMarkerSize(.75);
  ee15upn->SetMarkerSize(.75);
  ff010n->SetMarkerSize(.75);
  ff1015n->SetMarkerSize(.75);
  ff15upn->SetMarkerSize(.75);

  gg010n->SetTitle("");gg1015n->SetTitle("");gg15upn->SetTitle("");
  gg010n->GetYaxis()->SetTitle("Events / GeV");gg1015n->GetYaxis()->SetTitle("Events / GeV");gg15upn->GetYaxis()->SetTitle("Events / GeV");
  gg010n->SetMinimum(min);gg010n->SetMaximum(max);
  gg1015n->SetMinimum(min);gg1015n->SetMaximum(max);
  gg15upn->SetMinimum(min);gg15upn->SetMaximum(max);
  gg010n->Draw();
  gg1015n->SetMarkerColor(kRed);
  gg15upn->SetMarkerColor(kBlue);
  gg1015n->SetLineColor(kRed);
  gg15upn->SetLineColor(kBlue);
  gg1015n->Draw("SAME");
  gg15upn->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_NVertex_gg.png");
  ee010n->SetTitle("");ee010n->GetYaxis()->SetTitle("Events / GeV");
  ee010n->SetMinimum(min);ee010n->SetMaximum(max);
  ee010n->Draw();
  ee1015n->SetMarkerColor(kRed);
  ee15upn->SetMarkerColor(kBlue);
  ee1015n->SetLineColor(kRed);
  ee15upn->SetLineColor(kBlue);
  ee1015n->Draw("SAME");
  ee15upn->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_NVertex_ee.png");
  ff010n->SetTitle("");ff010n->GetYaxis()->SetTitle("Events / GeV");
  ff010n->SetMinimum(min);ff010n->SetMaximum(max);
  ff010n->Draw();
  ff1015n->SetMarkerColor(kRed);
  ff15upn->SetMarkerColor(kBlue);
  ff1015n->SetLineColor(kRed);
  ff15upn->SetLineColor(kBlue);
  ff1015n->Draw("SAME");
  ff15upn->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_NVertex_ff.png");

  gg010n->Draw();
  ee010n->SetMarkerColor(kRed);
  ff010n->SetMarkerColor(kBlue);
  ee010n->SetLineColor(kRed);
  ff010n->SetLineColor(kBlue);
  ee010n->Draw("SAME");
  ff010n->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_NVertex0_10.png");
  gg1015n->Draw();
  gg1015n->SetMarkerColor(kBlack);
  ee1015n->SetMarkerColor(kRed);
  ff1015n->SetMarkerColor(kBlue);
  gg1015n->SetLineColor(kBlack);
  ee1015n->SetLineColor(kRed);
  ff1015n->SetLineColor(kBlue);
  ee1015n->Draw("SAME");
  ff1015n->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_NVertex10_15.png");
  gg15upn->Draw();
  gg15upn->SetMarkerColor(kBlack);
  ee15upn->SetMarkerColor(kRed);
  ff15upn->SetMarkerColor(kBlue);
  gg15upn->SetLineColor(kBlack);
  ee15upn->SetLineColor(kRed);
  ff15upn->SetLineColor(kBlue);
  ee15upn->Draw("SAME");
  ff15upn->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Met_NVertex15up.png");

  c1->SetLogy(0);

  TH1* ggOVERee010=(TH1F*)gg010n->Clone();
  TH1* ggOVERee1015=(TH1F*)gg1015n->Clone();
  TH1* ggOVERee15up=(TH1F*)gg15upn->Clone();
  ggOVERee1015->Divide(ee1015n);
  ggOVERee010->Divide(ee010n);
  ggOVERee15up->Divide(ee15upn);
  ggOVERee010->GetYaxis()->SetRangeUser(0,3);
  ggOVERee1015->GetYaxis()->SetRangeUser(0,3);
  ggOVERee15up->GetYaxis()->SetRangeUser(0,3);
  ggOVERee1015->SetLineColor(kRed);
  ggOVERee1015->SetMarkerColor(kRed);
  ggOVERee15up->SetLineColor(kBlue);
  ggOVERee15up->SetMarkerColor(kBlue);
  ggOVERee010->SetMinimum(0);
  ggOVERee010->SetMaximum(5);
  ggOVERee010->Draw("PE");
  ggOVERee1015->Draw("PEsame");
  ggOVERee15up->Draw("PEsame");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MET_NVertexGGoverEE.png");


  TH1* ggOVERff010=(TH1F*)gg010n->Clone();
  TH1* ggOVERff1015=(TH1F*)gg1015n->Clone();
  TH1* ggOVERff15up=(TH1F*)gg15upn->Clone();
  ggOVERff1015->Divide(ff1015n);
  ggOVERff010->Divide(ff010n);
  ggOVERff15up->Divide(ff15upn);
  ggOVERff010->GetYaxis()->SetRangeUser(0,5);
  ggOVERff1015->GetYaxis()->SetRangeUser(0,5);
  ggOVERff15up->GetYaxis()->SetRangeUser(0,5);
  ggOVERff1015->SetLineColor(kRed);
  ggOVERff1015->SetMarkerColor(kRed);
  ggOVERff15up->SetLineColor(kBlue);
  ggOVERff15up->SetMarkerColor(kBlue);
  ggOVERff010->SetMinimum(0);
  ggOVERff010->SetMaximum(5);
  ggOVERff010->Draw("PE");
  ggOVERff1015->Draw("PEsame");
  ggOVERff15up->Draw("PEsame");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MET_NVertexGGoverFF.png");


  TH1* eeOVERff010=(TH1F*)ee010n->Clone();
  TH1* eeOVERff1015=(TH1F*)ee1015n->Clone();
  TH1* eeOVERff15up=(TH1F*)ee15upn->Clone();
  eeOVERff1015->Divide(ff1015n);
  eeOVERff010->Divide(ff010n);
  eeOVERff15up->Divide(ff15upn);
  eeOVERff010->GetYaxis()->SetRangeUser(0,5);
  eeOVERff1015->GetYaxis()->SetRangeUser(0,5);
  eeOVERff15up->GetYaxis()->SetRangeUser(0,5);
  eeOVERff010->SetLineColor(kBlack);
  eeOVERff010->SetMarkerColor(kBlack);
  eeOVERff1015->SetLineColor(kRed);
  eeOVERff1015->SetMarkerColor(kRed);
  eeOVERff15up->SetLineColor(kBlue);
  eeOVERff15up->SetMarkerColor(kBlue);
  eeOVERff010->SetMinimum(0);
  eeOVERff010->SetMaximum(5);
  eeOVERff010->Draw("PE");
  eeOVERff1015->Draw("PEsame");
  eeOVERff15up->Draw("PEsame");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/MET_NVertexEEoverFF.png");

  TH1F* ggLumi = (TH1F*)fin.Get("ggPerInstLumi");ggLumi->Sumw2();
  TH1F* egLumi = (TH1F*)fin.Get("egPerInstLumi");egLumi->Sumw2();
  TH1F* eeLumi = (TH1F*)fin.Get("eePerInstLumi");eeLumi->Sumw2();
  TH1F* ffLumi = (TH1F*)fin.Get("ffPerInstLumi");ffLumi->Sumw2();
  TH1F* Lumi = (TH1F*)fin.Get("InstLumi");Lumi->Sumw2();

  ggLumi->Rebin(2);
  egLumi->Rebin(2);
  eeLumi->Rebin(2);
  ffLumi->Rebin(2);
  Lumi->Rebin(2);

  ggLumi->Divide(Lumi);
  egLumi->Divide(Lumi);
  eeLumi->Divide(Lumi);
  ffLumi->Divide(Lumi);

  ggLumi->SetTitle(";Avg Inst Lumi (Hz/#mub);gg Events / Total Events");
  egLumi->SetTitle(";Avg Inst Lumi (Hz/#mub);eg Events / Total Events");
  eeLumi->SetTitle(";Avg Inst Lumi (Hz/#mub);ee Events / Total Events");
  ffLumi->SetTitle(";Avg Inst Lumi (Hz/#mub);ff Events / Total Events");

  ggLumi->GetYaxis()->SetRangeUser(.05,.07);
  egLumi->GetYaxis()->SetRangeUser(.035,.05);
  eeLumi->GetYaxis()->SetRangeUser(.3,.4);
  ffLumi->GetYaxis()->SetRangeUser(.03,.06);

  ggLumi->GetXaxis()->SetRangeUser(200,1250);
  egLumi->GetXaxis()->SetRangeUser(200,1250);
  eeLumi->GetXaxis()->SetRangeUser(200,1250);
  ffLumi->GetXaxis()->SetRangeUser(200,1250);

  ggLumi->Fit("pol1","","",300,1150);
  ggLumi->Draw("PE");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InstLumi_gg.png");
  egLumi->Draw("PE");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InstLumi_eg.png"); 
  eeLumi->Draw("PE");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InstLumi_ee.png"); 
  ffLumi->Fit("pol1","","",300,1150);
  ffLumi->Draw("PE");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/InstLumi_ff.png");
  /*
    cout<<"ee Invariant Mass Fit with One Loose Jet Requirement:"<<endl<<"signal events: "<<CBsigEE_JetReq<<endl;
    cout<<"eg Invariant Mass Fit with One Loose Jet Requirement:"<<endl<<"signal events: "<<CBsigEG_JetReq<<endl;
    cout<<"e->g fake rate with One Loose Jet Requirement f(e->g) = "<<CBsigEG_JetReq / (2*CBsigEE_JetReq + CBsigEG_JetReq)<< " +- "<< 
    ( CBsigEG_JetReq/( 2*CBsigEE_JetReq + CBsigEG_JetReq ) )*
    sqrt( ( CBsigEG_JetReq/CBsigEG_JetReq )*( CBsigEG_JetReq/CBsigEG_JetReq ) + ((2*CBsigEE_JetReq+CBsigEG_JetReq)/(2*CBsigEE_JetReq+CBsigEG_JetReq))*((2*CBsigEE_JetReq+CBsigEG_JetReq)/(2*CBsigEE_JetReq+CBsigEG_JetReq)) )
    <<endl;
  */
  cout<<"-----------------------------------"<<endl;
  cout<<"gg Events with MET>100GeV:"<<gg100up/*<< " +- "<<gg100upError*/<<endl;
  cout<<"EWK(eg) Events with MET>100GeV:"<<eg100up<< " +- "<<eg100upError<<endl;
  cout<<"ee Events with MET>100GeV:"<<ee100up<< " +- "<<ee100upError<<endl;
  cout<<"QCD(ee) Events with MET>100GeV:"<<ee100up+eg100up<< " +- "<<sqrt(ee100upError*ee100upError + eg100upError*eg100upError)<<endl;
  cout<<"ff Events with MET>100GeV:"<<ff100up<< " +- "<<ff100upError<<endl; 
  cout<<"QCD(ff) Events with MET>100GeV:"<<ff100up+eg100up<< " +- "<<sqrt(ff100upError*ff100upError + eg100upError*eg100upError)<<endl;
  //cout<<"QCDcombined Events with MET>100GeV:"<<comb100up<<endl;
  cout<<"-----------------------------------"<<endl; 
  cout<<"gg_JetReq Events with MET>100GeV:"<<gg100up_JetReq<</* " +- "<<gg100upError_JetReq<<*/endl;
  cout<<"EWK(eg_JetReq) Events with MET>100GeV:"<<eg100up_JetReq<< " +- "<<eg100upError_JetReq<<endl;
  cout<<"ee_JetReq Events with MET>100GeV:"<<ee100up_JetReq<< " +- "<<ee100upError_JetReq<<endl;
  cout<<"QCD(ee_JetReq) Events with MET>100GeV:"<<ee100up_JetReq+eg100up_JetReq<< " +- "<<sqrt(ee100upError_JetReq*ee100upError_JetReq + eg100upError_JetReq*eg100upError_JetReq)<<endl;
  cout<<"ff_JetReq Events with MET>100GeV:"<<ff100up_JetReq<< " +- "<<ff100upError_JetReq<<endl;
  cout<<"QCD(ff_JetReq) Events with MET>100GeV:"<<ff100up_JetReq+eg100up_JetReq<< " +- "<<sqrt(ff100upError_JetReq*ff100upError_JetReq + eg100upError_JetReq*eg100upError_JetReq)<<endl;
  //cout<<"QCDcombined_JetReq Events with MET>100GeV:"<<comb100up_JetReq<<endl;
  cout<<"-----------------------------------"<<endl; 
  cout<<"Normalization Factor ee        : "<<eeScale<<" +- "<<eeScaleErr<<endl;
  cout<<"Normalization Factor ee_JetReq : "<<eeScale_JetReq<<" +- "<<eeScaleErr_JetReq<<endl;
  cout<<"Normalization Factor ff        : "<<ffScale<<" +- "<<ffScaleErr<<endl;
  cout<<"Normalization Factor ff_JetReq : "<<ffScale_JetReq<<" +- "<<ffScaleErr_JetReq<<endl;
  cout<<"-----------------------------------"<<endl; 
 

  

  eeMetFromZZ->Scale(eeScale);
  eeMetFromWZ->Scale(eeScale);
  eeMetFromZZ_JetReq->Scale(eeScale_JetReq);
  eeMetFromWZ_JetReq->Scale(eeScale_JetReq);


  TH1F* ZZrebin = (TH1F*)eeMetFromZZ->Rebin(NmetBins,"ZZrebin",xbins);
  TH1F* WZrebin = (TH1F*)eeMetFromWZ->Rebin(NmetBins,"WZrebin",xbins);
  TH1F* ZZrebin_JetReq = (TH1F*)eeMetFromZZ_JetReq->Rebin(NmetBins,"ZZrebin_JetReq",xbins);
  TH1F* WZrebin_JetReq = (TH1F*)eeMetFromWZ_JetReq->Rebin(NmetBins,"WZrebin_JetReq",xbins);

  ZZrebin->SetTitle("");WZrebin->SetTitle("");ZZrebin_JetReq->SetTitle("");WZrebin_JetReq->SetTitle("");
  //ZZrebin->GetYaxis()->SetTitle("Events / GeV");WZrebin->GetYaxis()->SetTitle("Events / GeV");ZZrebin_JetReq->GetYaxis()->SetTitle("Events / GeV");WZrebin_JetReq->GetYaxis()->SetTitle("Events / GeV");
  /*
    for(int i=1;i<=eeOverBinWidth_JetReq->GetNbinsX();i++){
    WZrebin->SetBinContent(i,WZrebin->GetBinContent(i)/WZrebin->GetBinWidth(i));
    ZZrebin->SetBinContent(i,ZZrebin->GetBinContent(i)/ZZrebin->GetBinWidth(i));
    WZrebin->SetBinError(i,WZrebin->GetBinError(i)/WZrebin->GetBinWidth(i));
    ZZrebin->SetBinError(i,ZZrebin->GetBinError(i)/ZZrebin->GetBinWidth(i));
    WZrebin_JetReq->SetBinContent(i,WZrebin_JetReq->GetBinContent(i)/WZrebin_JetReq->GetBinWidth(i));
    ZZrebin_JetReq->SetBinContent(i,ZZrebin_JetReq->GetBinContent(i)/ZZrebin_JetReq->GetBinWidth(i));
    WZrebin_JetReq->SetBinError(i,WZrebin_JetReq->GetBinError(i)/WZrebin_JetReq->GetBinWidth(i));
    ZZrebin_JetReq->SetBinError(i,ZZrebin_JetReq->GetBinError(i)/ZZrebin_JetReq->GetBinWidth(i));
    }
  */
  eeMetFromWZ->Draw();
  eeMetFromZZ->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/eeMetFromZZandWZ_Scaled.png");
  eeMetFromWZ_JetReq->Draw();
  eeMetFromZZ_JetReq->Draw("SAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/eeMetFromZZandWZ_JetReq_Scaled.png");
  c1->SetLogy(1);  
  //eeOverBinWidth->Draw("hist");
  WZrebin->Draw("hist");
  ZZrebin->Draw("histSAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/eeMetFromZZandWZ_Scaledrebin.png");
  //eeOverBinWidth_JetReq->Draw("hist");
  WZrebin_JetReq->Draw("hist");
  ZZrebin_JetReq->Draw("histSAME");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/eeMetFromZZandWZ_JetReq_Scaledrebin.png");
  float eeMetScaledFromZZ100up=eeMetFromZZ->Integral(eeMetFromZZ->FindBin(100.1),-1);
  float eeMetScaledFromWZ100up=eeMetFromWZ->Integral(eeMetFromWZ->FindBin(100.1),-1);
  float eeMetScaledFromZZ100up_JetReq=eeMetFromZZ_JetReq->Integral(eeMetFromZZ_JetReq->FindBin(100.1),-1);
  float eeMetScaledFromWZ100up_JetReq=eeMetFromWZ_JetReq->Integral(eeMetFromWZ_JetReq->FindBin(100.1),-1);
  c1->SetLogy(0);  


  cout<<"Latex table:"<<endl<<endl;

  //cout<<"\\begin{center}"<<endl;
  cout<<"   \\begin{tabular}{ | l | c | c | c | c |}\n      \\hline"<<endl;
  cout<<"      Type & MET$<$20 & 30$<$MET$<$50 & MET$>$50 & MET$>$100 \\\\ \\hline"<<endl;
  cout<<"      $\\gamma\\gamma$ & "<<gg0_20<<" +- "<<gg0_20Error<<" & "<<gg30_50<<" +- "<<gg30_50Error<<" & "<<gg50up<<" +- "<<gg50upError<<" & "<<gg100up<<" +- "<<gg100upError<<" \\\\ \\hline"<<endl;
  cout<<"      ff QCD background & "<<ff0_20<<" +- "<<ff0_20Error<<" & "<<ff30_50<<" +- "<<ff30_50Error<<" & "<<ff50up<<" +- "<<ff50upError<<" & "<<ff100up<<" +- "<<ff100upError<<" \\\\ \\hline"<<endl;
  cout<<"      ee QCD background & "<<ee0_20<<" +- "<<ee0_20Error<<" & "<<ee30_50<<" +- "<<ee30_50Error<<" & "<<ee50up<<" +- "<<ee50upError<<" & "<<ee100up<<" +- "<<ee100upError<<" \\\\ \\hline"<<endl;
  cout<<"      EWK background & "<<eg0_20<<" +- "<<eg0_20Error<<" & "<<eg30_50<<" +- "<<eg30_50Error<<" & "<<eg50up<<" +- "<<eg50upError<<" & "<<eg100up<<" +- "<<eg100upError<<" \\\\ \\hline"<<endl;
  cout<<"      total background (ff) & "<<QCDff0_20<<" +- "<<QCDff0_20Error<<" & "<<QCDff30_50<<" +- "<<QCDff30_50Error<<" & "<<QCDff50up<<" +- "<<QCDff50upError<<" & "<<QCDff100up<<" +- "<<QCDff100upError<<" \\\\ \\hline"<<endl;
  cout<<"      total background (ee) & "<<QCDee0_20<<" +- "<<QCDee0_20Error<<" & "<<QCDee30_50<<" +- "<<QCDee30_50Error<<" & "<<QCDee50up<<" +- "<<QCDee50upError<<" & "<<QCDee100up<<" +- "<<QCDee100upError<<" \\\\ \\hline"<<endl;
  cout<<"      $\\gamma\\gamma$ $\\ge$1 Jet & "<<gg0_20_JetReq<<" +- "<<gg0_20Error_JetReq<<" & "<<gg30_50_JetReq<<" +- "<<gg30_50Error_JetReq<<" & "<<gg50up_JetReq<<" +- "<<gg50upError_JetReq<<" & "<<gg100up_JetReq<<" +- "<<gg100upError_JetReq<<" \\\\ \\hline"<<endl;
  cout<<"      ff QCD background $\\ge$1 Jet & "<<ff0_20_JetReq<<" +- "<<ff0_20Error_JetReq<<" & "<<ff30_50_JetReq<<" +- "<<ff30_50Error_JetReq<<" & "<<ff50up_JetReq<<" +- "<<ff50upError_JetReq<<" & "<<ff100up_JetReq<<" +- "<<ff100upError_JetReq<<" \\\\ \\hline"<<endl;
  cout<<"      ee QCD background $\\ge$1 Jet & "<<ee0_20_JetReq<<" +- "<<ee0_20Error_JetReq<<" & "<<ee30_50_JetReq<<" +- "<<ee30_50Error_JetReq<<" & "<<ee50up_JetReq<<" +- "<<ee50upError_JetReq<<" & "<<ee100up_JetReq<<" +- "<<ee100upError_JetReq<<" \\\\ \\hline"<<endl;
  cout<<"      EWK background $\\ge$1 Jet & "<<eg0_20_JetReq<<" +- "<<eg0_20Error_JetReq<<" & "<<eg30_50_JetReq<<" +- "<<eg30_50Error_JetReq<<" & "<<eg50up_JetReq<<" +- "<<eg50upError_JetReq<<" & "<<eg100up_JetReq<<" +- "<<eg100upError_JetReq<<" \\\\ \\hline"<<endl;
  cout<<"      total background (ff) $\\ge$1 Jet & "<<QCDff0_20_JetReq<<" +- "<<QCDff0_20Error_JetReq<<" & "<<QCDff30_50_JetReq<<" +- "<<QCDff30_50Error_JetReq<<" & "<<QCDff50up_JetReq<<" +- "<<QCDff50upError_JetReq<<" & "<<QCDff100up_JetReq<<" +- "<<QCDff100upError_JetReq<<" \\\\ \\hline"<<endl;
  cout<<"      total background (ee) $\\ge$1 Jet & "<<QCDee0_20_JetReq<<" +- "<<QCDee0_20Error_JetReq<<" & "<<QCDee30_50_JetReq<<" +- "<<QCDee30_50Error_JetReq<<" & "<<QCDee50up_JetReq<<" +- "<<QCDee50upError_JetReq<<" & "<<QCDee100up_JetReq<<" +- "<<QCDee100upError_JetReq<<" \\\\ \\hline"<<endl;
  cout<<"   \\end{tabular}"<<endl;
  //  cout<<"\\end{center}"<<endl;
  

  float ffReweightErr0_20 = AddInQuadrature(reweightErrff,1,4);
  float ffStatErr0_20 = AddInQuadrature(statErrff,1,4);
  float ffNormErr0_20 = AddInQuadrature(normErrff,1,4);
  float ffDiffFromeeErr0_20 = fabs(ff0_20-ee0_20);//AddUp(diffFromeeErrorff,1,4);
  vector<float> ffErr0_20;
  ffErr0_20.push_back(ffReweightErr0_20);
  ffErr0_20.push_back(ffStatErr0_20);
  ffErr0_20.push_back(ffNormErr0_20);
  ffErr0_20.push_back(ffDiffFromeeErr0_20);
  float fferr0_20 = AddInQuadrature(ffErr0_20,1,4);
  cout<<endl;
  cout <<"ffReweightErr0_20   : "<<ffReweightErr0_20<<endl;
  cout <<"ffStatErr0_20       : "<<ffStatErr0_20<<endl;
  cout <<"ffNormErr0_20       : "<<ffNormErr0_20<<endl;
  cout <<"ffDiffFromeeErr0_20 : "<<ffDiffFromeeErr0_20<<endl;
  cout <<"fferr0_20           : "<<fferr0_20<<endl;


  float ffReweightErr0_20_JetReq = AddInQuadrature(reweightErrff_JetReq,1,4);
  float ffStatErr0_20_JetReq = AddInQuadrature(statErrff_JetReq,1,4);
  float ffNormErr0_20_JetReq = AddInQuadrature(normErrff_JetReq,1,4);
  float ffDiffFromeeErr0_20_JetReq = ff0_20_JetReq-ee0_20_JetReq;//AddUp(diffFromeeErrorff_JetReq,1,4);
  cout<<endl;
  vector<float> ffErr0_20_JetReq;
  ffErr0_20_JetReq.push_back(ffReweightErr0_20_JetReq);
  ffErr0_20_JetReq.push_back(ffStatErr0_20_JetReq);
  ffErr0_20_JetReq.push_back(ffNormErr0_20_JetReq);
  ffErr0_20_JetReq.push_back(ffDiffFromeeErr0_20_JetReq);
  float fferr0_20_JetReq = AddInQuadrature(ffErr0_20_JetReq,1,4);
  cout <<"ffReweightErr0_20_JetReq   : "<<ffReweightErr0_20_JetReq<<endl;
  cout <<"ffStatErr0_20_JetReq       : "<<ffStatErr0_20_JetReq<<endl;
  cout <<"ffNormErr0_20_JetReq       : "<<ffNormErr0_20_JetReq<<endl;
  cout <<"ffDiffFromeeErr0_20_JetReq : "<<ffDiffFromeeErr0_20_JetReq<<endl;
  cout <<"fferr0_20_JetReq           : "<<fferr0_20_JetReq<<endl;


  float ffReweightErr0_20_2JetReq = AddInQuadrature(reweightErrff_2JetReq,1,4);
  float ffStatErr0_20_2JetReq = AddInQuadrature(statErrff_2JetReq,1,4);
  float ffNormErr0_20_2JetReq = AddInQuadrature(normErrff_2JetReq,1,4);
  float ffDiffFromeeErr0_20_2JetReq = ff0_20_2JetReq-ee0_20_2JetReq;//AddUp(diffFromeeErrorff_2JetReq,1,4);
  cout<<endl;
  vector<float> ffErr0_20_2JetReq;
  ffErr0_20_2JetReq.push_back(ffReweightErr0_20_2JetReq);
  ffErr0_20_2JetReq.push_back(ffStatErr0_20_2JetReq);
  ffErr0_20_2JetReq.push_back(ffNormErr0_20_2JetReq);
  ffErr0_20_2JetReq.push_back(ffDiffFromeeErr0_20_2JetReq);
  float fferr0_20_2JetReq = AddInQuadrature(ffErr0_20_2JetReq,1,4);
  cout <<"ffReweightErr0_20_2JetReq   : "<<ffReweightErr0_20_2JetReq<<endl;
  cout <<"ffStatErr0_20_2JetReq       : "<<ffStatErr0_20_2JetReq<<endl;
  cout <<"ffNormErr0_20_2JetReq       : "<<ffNormErr0_20_2JetReq<<endl;
  cout <<"ffDiffFromeeErr0_20_2JetReq : "<<ffDiffFromeeErr0_20_2JetReq<<endl;
  cout <<"fferr0_20_2JetReq           : "<<fferr0_20_2JetReq<<endl;



  float eeReweightErr0_20 = AddInQuadrature(reweightErree,1,4);
  float eeStatErr0_20 = AddInQuadrature(statErree,1,4);
  float eeNormErr0_20 = AddInQuadrature(normErree,1,4);
  float eeDiffFromffErr0_20 = fabs(ff0_20-ee0_20);//AddUp(diffFromffErroree,1,4);
  vector<float> eeErr0_20;
  eeErr0_20.push_back(eeReweightErr0_20);
  eeErr0_20.push_back(eeStatErr0_20);
  eeErr0_20.push_back(eeNormErr0_20);
  eeErr0_20.push_back(ffDiffFromeeErr0_20);
  float eeerr0_20 = AddInQuadrature(eeErr0_20,1,4);
  cout<<endl;
  cout <<"eeReweightErr0_20   : "<<eeReweightErr0_20<<endl;
  cout <<"eeStatErr0_20       : "<<eeStatErr0_20<<endl;
  cout <<"eeNormErr0_20       : "<<eeNormErr0_20<<endl;

  float eeReweightErr0_20_JetReq = AddInQuadrature(reweightErree_JetReq,1,4);
  float eeStatErr0_20_JetReq = AddInQuadrature(statErree_JetReq,1,4);
  float eeNormErr0_20_JetReq = AddInQuadrature(normErree_JetReq,1,4);
  float eeDiffFromffErr0_20_JetReq = fabs(ff0_20_JetReq-ee0_20_JetReq);//AddUp(diffFromffErroree_JetReq,1,4);
  vector<float> eeErr0_20_JetReq;
  eeErr0_20_JetReq.push_back(eeReweightErr0_20_JetReq);
  eeErr0_20_JetReq.push_back(eeStatErr0_20_JetReq);
  eeErr0_20_JetReq.push_back(eeNormErr0_20_JetReq);
  eeErr0_20_JetReq.push_back(ffDiffFromeeErr0_20_JetReq);
  float eeerr0_20_JetReq = AddInQuadrature(eeErr0_20_JetReq,1,4);
  cout<<endl;
  cout <<"eeReweightErr0_20_JetReq   : "<<eeReweightErr0_20_JetReq<<endl;
  cout <<"eeStatErr0_20_JetReq       : "<<eeStatErr0_20_JetReq<<endl;
  cout <<"eeNormErr0_20_JetReq       : "<<eeNormErr0_20_JetReq<<endl;

  float eeReweightErr0_20_2JetReq = AddInQuadrature(reweightErree_2JetReq,1,4);
  float eeStatErr0_20_2JetReq = AddInQuadrature(statErree_2JetReq,1,4);
  float eeNormErr0_20_2JetReq = AddInQuadrature(normErree_2JetReq,1,4);
  float eeDiffFromffErr0_20_2JetReq = fabs(ff0_20_2JetReq-ee0_20_2JetReq);//AddUp(diffFromffErroree_2JetReq,1,4);
  vector<float> eeErr0_20_2JetReq;
  eeErr0_20_2JetReq.push_back(eeReweightErr0_20_2JetReq);
  eeErr0_20_2JetReq.push_back(eeStatErr0_20_2JetReq);
  eeErr0_20_2JetReq.push_back(eeNormErr0_20_2JetReq);
  eeErr0_20_2JetReq.push_back(ffDiffFromeeErr0_20_2JetReq);
  float eeerr0_20_2JetReq = AddInQuadrature(eeErr0_20_2JetReq,1,4);
  cout<<endl;
  cout <<"eeReweightErr0_20_2JetReq   : "<<eeReweightErr0_20_2JetReq<<endl;
  cout <<"eeStatErr0_20_2JetReq       : "<<eeStatErr0_20_2JetReq<<endl;
  cout <<"eeNormErr0_20_2JetReq       : "<<eeNormErr0_20_2JetReq<<endl;

  float egStatErr0_20 = AddInQuadrature(statErreg,1,4);
  float egNormErr0_20 = AddInQuadrature(normErreg,1,4);
  cout<<endl;
  cout <<"egStatErr0_20       : "<<egStatErr0_20<<endl;
  cout <<"egNormErr0_20       : "<<egNormErr0_20<<endl;
  vector<float> egErr0_20;
  egErr0_20.push_back(egStatErr0_20);
  egErr0_20.push_back(egNormErr0_20);
  float egerr0_20 = AddInQuadrature(egErr0_20,1,2);
  float QcdFFerr0_20=sqrt(egerr0_20*egerr0_20+fferr0_20*fferr0_20);
  float QcdEEerr0_20=sqrt(egerr0_20*egerr0_20+eeerr0_20*eeerr0_20);

  float egStatErr0_20_JetReq = AddInQuadrature(statErreg_JetReq,1,4);
  float egNormErr0_20_JetReq = AddInQuadrature(normErreg_JetReq,1,4);
  cout<<endl;
  cout <<"egStatErr0_20_JetReq       : "<<egStatErr0_20_JetReq<<endl;
  cout <<"egNormErr0_20_JetReq       : "<<egNormErr0_20_JetReq<<endl;
  vector<float> egErr0_20_JetReq;
  egErr0_20_JetReq.push_back(egStatErr0_20_JetReq);
  egErr0_20_JetReq.push_back(egNormErr0_20_JetReq);
  float egerr0_20_JetReq = AddInQuadrature(egErr0_20_JetReq,1,2);
  float QcdFFerr0_20_JetReq=sqrt(egerr0_20_JetReq*egerr0_20_JetReq+fferr0_20_JetReq*fferr0_20_JetReq);
  float QcdEEerr0_20_JetReq=sqrt(egerr0_20_JetReq*egerr0_20_JetReq+eeerr0_20_JetReq*eeerr0_20_JetReq);

  float egStatErr0_20_2JetReq = AddInQuadrature(statErreg_2JetReq,1,4);
  float egNormErr0_20_2JetReq = AddInQuadrature(normErreg_2JetReq,1,4);
  cout<<endl;
  cout <<"egStatErr0_20_2JetReq       : "<<egStatErr0_20_2JetReq<<endl;
  cout <<"egNormErr0_20_2JetReq       : "<<egNormErr0_20_2JetReq<<endl;
  vector<float> egErr0_20_2JetReq;
  egErr0_20_2JetReq.push_back(egStatErr0_20_2JetReq);
  egErr0_20_2JetReq.push_back(egNormErr0_20_2JetReq);
  float egerr0_20_2JetReq = AddInQuadrature(egErr0_20_2JetReq,1,2);
  float QcdFFerr0_20_2JetReq=sqrt(egerr0_20_2JetReq*egerr0_20_2JetReq+fferr0_20_2JetReq*fferr0_20_2JetReq);
  float QcdEEerr0_20_2JetReq=sqrt(egerr0_20_2JetReq*egerr0_20_2JetReq+eeerr0_20_2JetReq*eeerr0_20_2JetReq);

  float ffReweightErr30_50 = AddInQuadrature(reweightErrff,7,10);
  float ffStatErr30_50 = AddInQuadrature(statErrff,7,10);
  float ffNormErr30_50 = AddInQuadrature(normErrff,7,10);
  float ffDiffFromeeErr30_50 = fabs(ff30_50 - ee30_50);//AddUp(diffFromeeErrorff,7,10);
  cout<<endl;
  cout <<"ffReweightErr30_50   : "<<ffReweightErr30_50<<endl;
  cout <<"ffStatErr30_50       : "<<ffStatErr30_50<<endl;
  cout <<"ffNormErr30_50       : "<<ffNormErr30_50<<endl;
  cout <<"ffDiffFromeeErr30_50 : "<<ffDiffFromeeErr30_50<<endl;
  vector<float> ffErr30_50;
  ffErr30_50.push_back(ffReweightErr30_50);
  ffErr30_50.push_back(ffStatErr30_50);
  ffErr30_50.push_back(ffNormErr30_50);
  ffErr30_50.push_back(ffDiffFromeeErr30_50);
  float fferr30_50 = AddInQuadrature(ffErr30_50,1,4);

  float ffReweightErr30_50_JetReq = AddInQuadrature(reweightErrff_JetReq,7,10);
  float ffStatErr30_50_JetReq = AddInQuadrature(statErrff_JetReq,7,10);
  float ffNormErr30_50_JetReq = AddInQuadrature(normErrff_JetReq,7,10);
  float ffDiffFromeeErr30_50_JetReq = fabs(ff30_50_JetReq - ee30_50_JetReq);//AddUp(diffFromeeErrorff_JetReq,7,10);
  cout <<"ffReweightErr30_50_JetReq   : "<<ffReweightErr30_50_JetReq<<endl;
  cout <<"ffStatErr30_50_JetReq       : "<<ffStatErr30_50_JetReq<<endl;
  cout <<"ffNormErr30_50_JetReq       : "<<ffNormErr30_50_JetReq<<endl;
  cout <<"ffDiffFromeeErr30_50_JetReq : "<<ffDiffFromeeErr30_50_JetReq<<endl;
  vector<float> ffErr30_50_JetReq;
  ffErr30_50_JetReq.push_back(ffReweightErr30_50_JetReq);
  ffErr30_50_JetReq.push_back(ffStatErr30_50_JetReq);
  ffErr30_50_JetReq.push_back(ffNormErr30_50_JetReq);
  ffErr30_50_JetReq.push_back(ffDiffFromeeErr30_50_JetReq);
  float fferr30_50_JetReq = AddInQuadrature(ffErr30_50_JetReq,1,4);

  float ffReweightErr30_50_2JetReq = AddInQuadrature(reweightErrff_2JetReq,7,10);
  float ffStatErr30_50_2JetReq = AddInQuadrature(statErrff_2JetReq,7,10);
  float ffNormErr30_50_2JetReq = AddInQuadrature(normErrff_2JetReq,7,10);
  float ffDiffFromeeErr30_50_2JetReq = fabs(ff30_50_2JetReq - ee30_50_2JetReq);//AddUp(diffFromeeErrorff_2JetReq,7,10);
  cout <<"ffReweightErr30_50_2JetReq   : "<<ffReweightErr30_50_2JetReq<<endl;
  cout <<"ffStatErr30_50_2JetReq       : "<<ffStatErr30_50_2JetReq<<endl;
  cout <<"ffNormErr30_50_2JetReq       : "<<ffNormErr30_50_2JetReq<<endl;
  cout <<"ffDiffFromeeErr30_50_2JetReq : "<<ffDiffFromeeErr30_50_2JetReq<<endl;
  vector<float> ffErr30_50_2JetReq;
  ffErr30_50_2JetReq.push_back(ffReweightErr30_50_2JetReq);
  ffErr30_50_2JetReq.push_back(ffStatErr30_50_2JetReq);
  ffErr30_50_2JetReq.push_back(ffNormErr30_50_2JetReq);
  ffErr30_50_2JetReq.push_back(ffDiffFromeeErr30_50_2JetReq);
  float fferr30_50_2JetReq = AddInQuadrature(ffErr30_50_2JetReq,1,4);

  float eeReweightErr30_50 = AddInQuadrature(reweightErree,7,10);
  float eeStatErr30_50 = AddInQuadrature(statErree,7,10);
  float eeNormErr30_50 = AddInQuadrature(normErree,7,10);
  float eeDiffFromffErr30_50 = fabs(ff30_50 - ee30_50);//AddUp(diffFromffErroree,7,10);
  cout<<endl;
  cout <<"eeReweightErr30_50   : "<<eeReweightErr30_50<<endl;
  cout <<"eeStatErr30_50       : "<<eeStatErr30_50<<endl;
  cout <<"eeNormErr30_50       : "<<eeNormErr30_50<<endl;
  vector<float> eeErr30_50;
  eeErr30_50.push_back(eeReweightErr30_50);
  eeErr30_50.push_back(eeStatErr30_50);
  eeErr30_50.push_back(eeNormErr30_50);
  eeErr30_50.push_back(ffDiffFromeeErr30_50);
  float eeerr30_50 = AddInQuadrature(eeErr30_50,1,4);

  float eeReweightErr30_50_JetReq = AddInQuadrature(reweightErree_JetReq,7,10);
  float eeStatErr30_50_JetReq = AddInQuadrature(statErree_JetReq,7,10);
  float eeNormErr30_50_JetReq = AddInQuadrature(normErree_JetReq,7,10);
  float eeDiffFromffErr30_50_JetReq = fabs(ff30_50_JetReq - ee30_50_JetReq);//AddUp(diffFromffErroree_JetReq,7,10);
  cout<<endl;
  cout <<"eeReweightErr30_50_JetReq   : "<<eeReweightErr30_50_JetReq<<endl;
  cout <<"eeStatErr30_50_JetReq       : "<<eeStatErr30_50_JetReq<<endl;
  cout <<"eeNormErr30_50_JetReq       : "<<eeNormErr30_50_JetReq<<endl;
  vector<float> eeErr30_50_JetReq;
  eeErr30_50_JetReq.push_back(eeReweightErr30_50_JetReq);
  eeErr30_50_JetReq.push_back(eeStatErr30_50_JetReq);
  eeErr30_50_JetReq.push_back(eeNormErr30_50_JetReq);
  eeErr30_50_JetReq.push_back(ffDiffFromeeErr30_50_JetReq);
  float eeerr30_50_JetReq = AddInQuadrature(eeErr30_50_JetReq,1,4);

  float eeReweightErr30_50_2JetReq = AddInQuadrature(reweightErree_2JetReq,7,10);
  float eeStatErr30_50_2JetReq = AddInQuadrature(statErree_2JetReq,7,10);
  float eeNormErr30_50_2JetReq = AddInQuadrature(normErree_2JetReq,7,10);
  float eeDiffFromffErr30_50_2JetReq = fabs(ff30_50_2JetReq - ee30_50_2JetReq);//AddUp(diffFromffErroree_2JetReq,7,10);
  cout<<endl;
  cout <<"eeReweightErr30_50_2JetReq   : "<<eeReweightErr30_50_2JetReq<<endl;
  cout <<"eeStatErr30_50_2JetReq       : "<<eeStatErr30_50_2JetReq<<endl;
  cout <<"eeNormErr30_50_2JetReq       : "<<eeNormErr30_50_2JetReq<<endl;
  vector<float> eeErr30_50_2JetReq;
  eeErr30_50_2JetReq.push_back(eeReweightErr30_50_2JetReq);
  eeErr30_50_2JetReq.push_back(eeStatErr30_50_2JetReq);
  eeErr30_50_2JetReq.push_back(eeNormErr30_50_2JetReq);
  eeErr30_50_2JetReq.push_back(ffDiffFromeeErr30_50_2JetReq);
  float eeerr30_50_2JetReq = AddInQuadrature(eeErr30_50_2JetReq,1,4);

  float egStatErr30_50 = AddInQuadrature(statErreg,7,10);
  float egNormErr30_50 = AddInQuadrature(normErreg,7,10);
  cout<<endl;
  cout <<"egStatErr30_50       : "<<egStatErr30_50<<endl;
  cout <<"egNormErr30_50       : "<<egNormErr30_50<<endl;
  vector<float> egErr30_50;
  egErr30_50.push_back(egStatErr30_50);
  egErr30_50.push_back(egNormErr30_50);
  float egerr30_50 = AddInQuadrature(egErr30_50,1,2);
  float QcdFFerr30_50=sqrt(egerr30_50*egerr30_50+fferr30_50*fferr30_50);
  float QcdEEerr30_50=sqrt(egerr30_50*egerr30_50+eeerr30_50*eeerr30_50);

  float egStatErr30_50_JetReq = AddInQuadrature(statErreg_JetReq,7,10);
  float egNormErr30_50_JetReq = AddInQuadrature(normErreg_JetReq,7,10);
  cout<<endl;
  cout <<"egStatErr30_50_JetReq       : "<<egStatErr30_50_JetReq<<endl;
  cout <<"egNormErr30_50_JetReq       : "<<egNormErr30_50_JetReq<<endl;
  vector<float> egErr30_50_JetReq;
  egErr30_50_JetReq.push_back(egStatErr30_50_JetReq);
  egErr30_50_JetReq.push_back(egNormErr30_50_JetReq);
  float egerr30_50_JetReq = AddInQuadrature(egErr30_50_JetReq,1,2);
  float QcdFFerr30_50_JetReq=sqrt(egerr30_50_JetReq*egerr30_50_JetReq+fferr30_50_JetReq*fferr30_50_JetReq);
  float QcdEEerr30_50_JetReq=sqrt(egerr30_50_JetReq*egerr30_50_JetReq+eeerr30_50_JetReq*eeerr30_50_JetReq);

  float egStatErr30_50_2JetReq = AddInQuadrature(statErreg_2JetReq,7,10);
  float egNormErr30_50_2JetReq = AddInQuadrature(normErreg_2JetReq,7,10);
  cout<<endl;
  cout <<"egStatErr30_50_2JetReq       : "<<egStatErr30_50_2JetReq<<endl;
  cout <<"egNormErr30_50_2JetReq       : "<<egNormErr30_50_2JetReq<<endl;
  vector<float> egErr30_50_2JetReq;
  egErr30_50_2JetReq.push_back(egStatErr30_50_2JetReq);
  egErr30_50_2JetReq.push_back(egNormErr30_50_2JetReq);
  float egerr30_50_2JetReq = AddInQuadrature(egErr30_50_2JetReq,1,2);
  float QcdFFerr30_50_2JetReq=sqrt(egerr30_50_2JetReq*egerr30_50_2JetReq+fferr30_50_2JetReq*fferr30_50_2JetReq);
  float QcdEEerr30_50_2JetReq=sqrt(egerr30_50_2JetReq*egerr30_50_2JetReq+eeerr30_50_2JetReq*eeerr30_50_2JetReq);

  float ffReweightErr50up = AddInQuadrature(reweightErrff,11,-1);
  float ffStatErr50up = AddInQuadrature(statErrff,11,-1);
  float ffNormErr50up = AddInQuadrature(normErrff,11,-1);
  float ffDiffFromeeErr50up = fabs(ff50up - ee50up);//AddUp(diffFromeeErrorff,11,-1);
  cout<<endl;
  cout <<"ffReweightErr50up   : "<<ffReweightErr50up<<endl;
  cout <<"ffStatErr50up       : "<<ffStatErr50up<<endl;
  cout <<"ffNormErr50up       : "<<ffNormErr50up<<endl;
  cout <<"ffDiffFromeeErr50up : "<<ffDiffFromeeErr50up<<endl;
  vector<float> ffErr50up;
  ffErr50up.push_back(ffReweightErr50up);
  ffErr50up.push_back(ffStatErr50up);
  ffErr50up.push_back(ffNormErr50up);
  ffErr50up.push_back(ffDiffFromeeErr50up);
  float fferr50up = AddInQuadrature(ffErr50up,1,4);

  float ffReweightErr50up_JetReq = AddInQuadrature(reweightErrff_JetReq,11,-1);
  float ffStatErr50up_JetReq = AddInQuadrature(statErrff_JetReq,11,-1);
  float ffNormErr50up_JetReq = AddInQuadrature(normErrff_JetReq,11,-1);
  float ffDiffFromeeErr50up_JetReq = fabs(ff50up_JetReq - ee50up_JetReq);//AddUp(diffFromeeErrorff_JetReq,11,-1);
  cout <<"ffReweightErr50up_JetReq   : "<<ffReweightErr50up_JetReq<<endl;
  cout <<"ffStatErr50up_JetReq       : "<<ffStatErr50up_JetReq<<endl;
  cout <<"ffNormErr50up_JetReq       : "<<ffNormErr50up_JetReq<<endl;
  cout <<"ffDiffFromeeErr50up_JetReq : "<<ffDiffFromeeErr50up_JetReq<<endl;
  vector<float> ffErr50up_JetReq;
  ffErr50up_JetReq.push_back(ffReweightErr50up_JetReq);
  ffErr50up_JetReq.push_back(ffStatErr50up_JetReq);
  ffErr50up_JetReq.push_back(ffNormErr50up_JetReq);
  ffErr50up_JetReq.push_back(ffDiffFromeeErr50up_JetReq);
  float fferr50up_JetReq = AddInQuadrature(ffErr50up_JetReq,1,4);

  float ffReweightErr50up_2JetReq = AddInQuadrature(reweightErrff_2JetReq,11,-1);
  float ffStatErr50up_2JetReq = AddInQuadrature(statErrff_2JetReq,11,-1);
  float ffNormErr50up_2JetReq = AddInQuadrature(normErrff_2JetReq,11,-1);
  float ffDiffFromeeErr50up_2JetReq = fabs(ff50up_2JetReq - ee50up_2JetReq);//AddUp(diffFromeeErrorff_2JetReq,11,-1);
  cout <<"ffReweightErr50up_2JetReq   : "<<ffReweightErr50up_2JetReq<<endl;
  cout <<"ffStatErr50up_2JetReq       : "<<ffStatErr50up_2JetReq<<endl;
  cout <<"ffNormErr50up_2JetReq       : "<<ffNormErr50up_2JetReq<<endl;
  cout <<"ffDiffFromeeErr50up_2JetReq : "<<ffDiffFromeeErr50up_2JetReq<<endl;
  vector<float> ffErr50up_2JetReq;
  ffErr50up_2JetReq.push_back(ffReweightErr50up_2JetReq);
  ffErr50up_2JetReq.push_back(ffStatErr50up_2JetReq);
  ffErr50up_2JetReq.push_back(ffNormErr50up_2JetReq);
  ffErr50up_2JetReq.push_back(ffDiffFromeeErr50up_2JetReq);
  float fferr50up_2JetReq = AddInQuadrature(ffErr50up_2JetReq,1,4);

  float eeReweightErr50up = AddInQuadrature(reweightErree,11,-1);
  float eeStatErr50up = AddInQuadrature(statErree,11,-1);
  float eeNormErr50up = AddInQuadrature(normErree,11,-1);
  float eeDiffFromffErr50up = fabs(ff50up - ee50up);//AddUp(diffFromffErroree,11,-1);
  cout<<endl;
  cout <<"eeReweightErr50up   : "<<eeReweightErr50up<<endl;
  cout <<"eeStatErr50up       : "<<eeStatErr50up<<endl;
  cout <<"eeNormErr50up       : "<<eeNormErr50up<<endl;
  vector<float> eeErr50up;
  eeErr50up.push_back(eeReweightErr50up);
  eeErr50up.push_back(eeStatErr50up);
  eeErr50up.push_back(eeNormErr50up);
  eeErr50up.push_back(ffDiffFromeeErr50up);
  float eeerr50up = AddInQuadrature(eeErr50up,1,4);

  float eeReweightErr50up_JetReq = AddInQuadrature(reweightErree_JetReq,11,-1);
  float eeStatErr50up_JetReq = AddInQuadrature(statErree_JetReq,11,-1);
  float eeNormErr50up_JetReq = AddInQuadrature(normErree_JetReq,11,-1);
  float eeDiffFromffErr50up_JetReq = fabs(ff50up_JetReq - ee50up_JetReq);//AddUp(diffFromffErroree_JetReq,11,-1);
  cout<<endl;
  cout <<"eeReweightErr50up_JetReq   : "<<eeReweightErr50up_JetReq<<endl;
  cout <<"eeStatErr50up_JetReq       : "<<eeStatErr50up_JetReq<<endl;
  cout <<"eeNormErr50up_JetReq       : "<<eeNormErr50up_JetReq<<endl;
  vector<float> eeErr50up_JetReq;
  eeErr50up_JetReq.push_back(eeReweightErr50up_JetReq);
  eeErr50up_JetReq.push_back(eeStatErr50up_JetReq);
  eeErr50up_JetReq.push_back(eeNormErr50up_JetReq);
  eeErr50up_JetReq.push_back(ffDiffFromeeErr50up_JetReq);
  float eeerr50up_JetReq = AddInQuadrature(eeErr50up_JetReq,1,4);

  float eeReweightErr50up_2JetReq = AddInQuadrature(reweightErree_2JetReq,11,-1);
  float eeStatErr50up_2JetReq = AddInQuadrature(statErree_2JetReq,11,-1);
  float eeNormErr50up_2JetReq = AddInQuadrature(normErree_2JetReq,11,-1);
  float eeDiffFromffErr50up_2JetReq = fabs(ff50up_2JetReq - ee50up_2JetReq);//AddUp(diffFromffErroree_2JetReq,11,-1);
  cout<<endl;
  cout <<"eeReweightErr50up_2JetReq   : "<<eeReweightErr50up_2JetReq<<endl;
  cout <<"eeStatErr50up_2JetReq       : "<<eeStatErr50up_2JetReq<<endl;
  cout <<"eeNormErr50up_2JetReq       : "<<eeNormErr50up_2JetReq<<endl;
  vector<float> eeErr50up_2JetReq;
  eeErr50up_2JetReq.push_back(eeReweightErr50up_2JetReq);
  eeErr50up_2JetReq.push_back(eeStatErr50up_2JetReq);
  eeErr50up_2JetReq.push_back(eeNormErr50up_2JetReq);
  eeErr50up_2JetReq.push_back(ffDiffFromeeErr50up_2JetReq);
  float eeerr50up_2JetReq = AddInQuadrature(eeErr50up_2JetReq,1,4);

  float egStatErr50up = AddInQuadrature(statErreg,11,-1);
  float egNormErr50up = AddInQuadrature(normErreg,11,-1);
  cout<<endl;
  cout <<"egStatErr50up       : "<<egStatErr50up<<endl;
  cout <<"egNormErr50up       : "<<egNormErr50up<<endl;
  vector<float> egErr50up;
  egErr50up.push_back(egStatErr50up);
  egErr50up.push_back(egNormErr50up);
  float egerr50up = AddInQuadrature(egErr50up,1,2);
  float QcdFFerr50up=sqrt(egerr50up*egerr50up+fferr50up*fferr50up);
  float QcdEEerr50up=sqrt(egerr50up*egerr50up+eeerr50up*eeerr50up);

  float egStatErr50up_JetReq = AddInQuadrature(statErreg_JetReq,11,-1);
  float egNormErr50up_JetReq = AddInQuadrature(normErreg_JetReq,11,-1);
  cout<<endl;
  cout <<"egStatErr50up_JetReq       : "<<egStatErr50up_JetReq<<endl;
  cout <<"egNormErr50up_JetReq       : "<<egNormErr50up_JetReq<<endl;
  vector<float> egErr50up_JetReq;
  egErr50up_JetReq.push_back(egStatErr50up_JetReq);
  egErr50up_JetReq.push_back(egNormErr50up_JetReq);
  float egerr50up_JetReq = AddInQuadrature(egErr50up_JetReq,1,2);
  float QcdFFerr50up_JetReq=sqrt(egerr50up_JetReq*egerr50up_JetReq+fferr50up_JetReq*fferr50up_JetReq);
  float QcdEEerr50up_JetReq=sqrt(egerr50up_JetReq*egerr50up_JetReq+eeerr50up_JetReq*eeerr50up_JetReq);

  float egStatErr50up_2JetReq = AddInQuadrature(statErreg_2JetReq,11,-1);
  float egNormErr50up_2JetReq = AddInQuadrature(normErreg_2JetReq,11,-1);
  cout<<endl;
  cout <<"egStatErr50up_2JetReq       : "<<egStatErr50up_2JetReq<<endl;
  cout <<"egNormErr50up_2JetReq       : "<<egNormErr50up_2JetReq<<endl;
  vector<float> egErr50up_2JetReq;
  egErr50up_2JetReq.push_back(egStatErr50up_2JetReq);
  egErr50up_2JetReq.push_back(egNormErr50up_2JetReq);
  float egerr50up_2JetReq = AddInQuadrature(egErr50up_2JetReq,1,2);
  float QcdFFerr50up_2JetReq=sqrt(egerr50up_2JetReq*egerr50up_2JetReq+fferr50up_2JetReq*fferr50up_2JetReq);
  float QcdEEerr50up_2JetReq=sqrt(egerr50up_2JetReq*egerr50up_2JetReq+eeerr50up_2JetReq*eeerr50up_2JetReq);

  float ffReweightErr50_60 = AddInQuadrature(reweightErrff,11,12);
  float ffStatErr50_60 = AddInQuadrature(statErrff,11,12);
  float ffNormErr50_60 = AddInQuadrature(normErrff,11,12);
  float ffDiffFromeeErr50_60 = fabs(ff50_60 - ee50_60);
  vector<float> ffErr50_60;
  ffErr50_60.push_back(ffReweightErr50_60);
  ffErr50_60.push_back(ffStatErr50_60);
  ffErr50_60.push_back(ffNormErr50_60);
  ffErr50_60.push_back(ffDiffFromeeErr50_60);
  float fferr50_60 = AddInQuadrature(ffErr50_60,1,4);

  float ffReweightErr50_60_JetReq = AddInQuadrature(reweightErrff_JetReq,11,12);
  float ffStatErr50_60_JetReq = AddInQuadrature(statErrff_JetReq,11,12);
  float ffNormErr50_60_JetReq = AddInQuadrature(normErrff_JetReq,11,12);
  float ffDiffFromeeErr50_60_JetReq = fabs(ff50_60_JetReq - ee50_60_JetReq);
  vector<float> ffErr50_60_JetReq;
  ffErr50_60_JetReq.push_back(ffReweightErr50_60_JetReq);
  ffErr50_60_JetReq.push_back(ffStatErr50_60_JetReq);
  ffErr50_60_JetReq.push_back(ffNormErr50_60_JetReq);
  ffErr50_60_JetReq.push_back(ffDiffFromeeErr50_60_JetReq);
  float fferr50_60_JetReq = AddInQuadrature(ffErr50_60_JetReq,1,4);

  float ffReweightErr50_60_2JetReq = AddInQuadrature(reweightErrff_2JetReq,11,12);
  float ffStatErr50_60_2JetReq = AddInQuadrature(statErrff_2JetReq,11,12);
  float ffNormErr50_60_2JetReq = AddInQuadrature(normErrff_2JetReq,11,12);
  float ffDiffFromeeErr50_60_2JetReq = fabs(ff50_60_2JetReq - ee50_60_2JetReq);
  vector<float> ffErr50_60_2JetReq;
  ffErr50_60_2JetReq.push_back(ffReweightErr50_60_2JetReq);
  ffErr50_60_2JetReq.push_back(ffStatErr50_60_2JetReq);
  ffErr50_60_2JetReq.push_back(ffNormErr50_60_2JetReq);
  ffErr50_60_2JetReq.push_back(ffDiffFromeeErr50_60_2JetReq);
  float fferr50_60_2JetReq = AddInQuadrature(ffErr50_60_2JetReq,1,4);

  float eeReweightErr50_60 = AddInQuadrature(reweightErree,11,12);
  float eeStatErr50_60 = AddInQuadrature(statErree,11,12);
  float eeNormErr50_60 = AddInQuadrature(normErree,11,12);
  float eeDiffFromffErr50_60 = fabs(ff50_60 - ee50_60);
  vector<float> eeErr50_60;
  eeErr50_60.push_back(eeReweightErr50_60);
  eeErr50_60.push_back(eeStatErr50_60);
  eeErr50_60.push_back(eeNormErr50_60);
  eeErr50_60.push_back(ffDiffFromeeErr50_60);
  float eeerr50_60 = AddInQuadrature(eeErr50_60,1,4);

  float eeReweightErr50_60_JetReq = AddInQuadrature(reweightErree_JetReq,11,12);
  float eeStatErr50_60_JetReq = AddInQuadrature(statErree_JetReq,11,12);
  float eeNormErr50_60_JetReq = AddInQuadrature(normErree_JetReq,11,12);
  float eeDiffFromffErr50_60_JetReq = fabs(ff50_60_JetReq - ee50_60_JetReq);
  vector<float> eeErr50_60_JetReq;
  eeErr50_60_JetReq.push_back(eeReweightErr50_60_JetReq);
  eeErr50_60_JetReq.push_back(eeStatErr50_60_JetReq);
  eeErr50_60_JetReq.push_back(eeNormErr50_60_JetReq);
  eeErr50_60_JetReq.push_back(ffDiffFromeeErr50_60_JetReq);
  float eeerr50_60_JetReq = AddInQuadrature(eeErr50_60_JetReq,1,4);

  float eeReweightErr50_60_2JetReq = AddInQuadrature(reweightErree_2JetReq,11,12);
  float eeStatErr50_60_2JetReq = AddInQuadrature(statErree_2JetReq,11,12);
  float eeNormErr50_60_2JetReq = AddInQuadrature(normErree_2JetReq,11,12);
  float eeDiffFromffErr50_60_2JetReq = fabs(ff50_60_2JetReq - ee50_60_2JetReq);
  vector<float> eeErr50_60_2JetReq;
  eeErr50_60_2JetReq.push_back(eeReweightErr50_60_2JetReq);
  eeErr50_60_2JetReq.push_back(eeStatErr50_60_2JetReq);
  eeErr50_60_2JetReq.push_back(eeNormErr50_60_2JetReq);
  eeErr50_60_2JetReq.push_back(ffDiffFromeeErr50_60_2JetReq);
  float eeerr50_60_2JetReq = AddInQuadrature(eeErr50_60_2JetReq,1,4);

  float egStatErr50_60 = AddInQuadrature(statErreg,11,12);
  float egNormErr50_60 = AddInQuadrature(normErreg,11,12);
  vector<float> egErr50_60;
  egErr50_60.push_back(egStatErr50_60);
  egErr50_60.push_back(egNormErr50_60);
  float egerr50_60 = AddInQuadrature(egErr50_60,1,2);
  float QcdFFerr50_60=sqrt(egerr50_60*egerr50_60+fferr50_60*fferr50_60);
  float QcdEEerr50_60=sqrt(egerr50_60*egerr50_60+eeerr50_60*eeerr50_60);

  float egStatErr50_60_JetReq = AddInQuadrature(statErreg_JetReq,11,12);
  float egNormErr50_60_JetReq = AddInQuadrature(normErreg_JetReq,11,12);
  vector<float> egErr50_60_JetReq;
  egErr50_60_JetReq.push_back(egStatErr50_60_JetReq);
  egErr50_60_JetReq.push_back(egNormErr50_60_JetReq);
  float egerr50_60_JetReq = AddInQuadrature(egErr50_60_JetReq,1,2);
  float QcdFFerr50_60_JetReq=sqrt(egerr50_60_JetReq*egerr50_60_JetReq+fferr50_60_JetReq*fferr50_60_JetReq);
  float QcdEEerr50_60_JetReq=sqrt(egerr50_60_JetReq*egerr50_60_JetReq+eeerr50_60_JetReq*eeerr50_60_JetReq);

  float egStatErr50_60_2JetReq = AddInQuadrature(statErreg_2JetReq,11,12);
  float egNormErr50_60_2JetReq = AddInQuadrature(normErreg_2JetReq,11,12);
  vector<float> egErr50_60_2JetReq;
  egErr50_60_2JetReq.push_back(egStatErr50_60_2JetReq);
  egErr50_60_2JetReq.push_back(egNormErr50_60_2JetReq);
  float egerr50_60_2JetReq = AddInQuadrature(egErr50_60_2JetReq,1,2);
  float QcdFFerr50_60_2JetReq=sqrt(egerr50_60_2JetReq*egerr50_60_2JetReq+fferr50_60_2JetReq*fferr50_60_2JetReq);
  float QcdEEerr50_60_2JetReq=sqrt(egerr50_60_2JetReq*egerr50_60_2JetReq+eeerr50_60_2JetReq*eeerr50_60_2JetReq);

  float ffReweightErr60_70 = AddInQuadrature(reweightErrff,13,14);
  float ffStatErr60_70 = AddInQuadrature(statErrff,13,14);
  float ffNormErr60_70 = AddInQuadrature(normErrff,13,14);
  float ffDiffFromeeErr60_70 = fabs(ff60_70 - ee60_70);
  vector<float> ffErr60_70;
  ffErr60_70.push_back(ffReweightErr60_70);
  ffErr60_70.push_back(ffStatErr60_70);
  ffErr60_70.push_back(ffNormErr60_70);
  ffErr60_70.push_back(ffDiffFromeeErr60_70);
  float fferr60_70 = AddInQuadrature(ffErr60_70,1,4);

  float ffReweightErr60_70_JetReq = AddInQuadrature(reweightErrff_JetReq,13,14);
  float ffStatErr60_70_JetReq = AddInQuadrature(statErrff_JetReq,13,14);
  float ffNormErr60_70_JetReq = AddInQuadrature(normErrff_JetReq,13,14);
  float ffDiffFromeeErr60_70_JetReq = fabs(ff60_70_JetReq - ee60_70_JetReq);
  vector<float> ffErr60_70_JetReq;
  ffErr60_70_JetReq.push_back(ffReweightErr60_70_JetReq);
  ffErr60_70_JetReq.push_back(ffStatErr60_70_JetReq);
  ffErr60_70_JetReq.push_back(ffNormErr60_70_JetReq);
  ffErr60_70_JetReq.push_back(ffDiffFromeeErr60_70_JetReq);
  float fferr60_70_JetReq = AddInQuadrature(ffErr60_70_JetReq,1,4);

  float ffReweightErr60_70_2JetReq = AddInQuadrature(reweightErrff_2JetReq,13,14);
  float ffStatErr60_70_2JetReq = AddInQuadrature(statErrff_2JetReq,13,14);
  float ffNormErr60_70_2JetReq = AddInQuadrature(normErrff_2JetReq,13,14);
  float ffDiffFromeeErr60_70_2JetReq = fabs(ff60_70_2JetReq - ee60_70_2JetReq);
  vector<float> ffErr60_70_2JetReq;
  ffErr60_70_2JetReq.push_back(ffReweightErr60_70_2JetReq);
  ffErr60_70_2JetReq.push_back(ffStatErr60_70_2JetReq);
  ffErr60_70_2JetReq.push_back(ffNormErr60_70_2JetReq);
  ffErr60_70_2JetReq.push_back(ffDiffFromeeErr60_70_2JetReq);
  float fferr60_70_2JetReq = AddInQuadrature(ffErr60_70_2JetReq,1,4);

  float eeReweightErr60_70 = AddInQuadrature(reweightErree,13,14);
  float eeStatErr60_70 = AddInQuadrature(statErree,13,14);
  float eeNormErr60_70 = AddInQuadrature(normErree,13,14);
  float eeDiffFromffErr60_70 = fabs(ff60_70 - ee60_70);
  vector<float> eeErr60_70;
  eeErr60_70.push_back(eeReweightErr60_70);
  eeErr60_70.push_back(eeStatErr60_70);
  eeErr60_70.push_back(eeNormErr60_70);
  eeErr60_70.push_back(ffDiffFromeeErr60_70);
  float eeerr60_70 = AddInQuadrature(eeErr60_70,1,4);

  float eeReweightErr60_70_JetReq = AddInQuadrature(reweightErree_JetReq,13,14);
  float eeStatErr60_70_JetReq = AddInQuadrature(statErree_JetReq,13,14);
  float eeNormErr60_70_JetReq = AddInQuadrature(normErree_JetReq,13,14);
  float eeDiffFromffErr60_70_JetReq = fabs(ff60_70_JetReq - ee60_70_JetReq);
  vector<float> eeErr60_70_JetReq;
  eeErr60_70_JetReq.push_back(eeReweightErr60_70_JetReq);
  eeErr60_70_JetReq.push_back(eeStatErr60_70_JetReq);
  eeErr60_70_JetReq.push_back(eeNormErr60_70_JetReq);
  eeErr60_70_JetReq.push_back(ffDiffFromeeErr60_70_JetReq);
  float eeerr60_70_JetReq = AddInQuadrature(eeErr60_70_JetReq,1,4);

  float eeReweightErr60_70_2JetReq = AddInQuadrature(reweightErree_2JetReq,13,14);
  float eeStatErr60_70_2JetReq = AddInQuadrature(statErree_2JetReq,13,14);
  float eeNormErr60_70_2JetReq = AddInQuadrature(normErree_2JetReq,13,14);
  float eeDiffFromffErr60_70_2JetReq = fabs(ff60_70_2JetReq - ee60_70_2JetReq);
  vector<float> eeErr60_70_2JetReq;
  eeErr60_70_2JetReq.push_back(eeReweightErr60_70_2JetReq);
  eeErr60_70_2JetReq.push_back(eeStatErr60_70_2JetReq);
  eeErr60_70_2JetReq.push_back(eeNormErr60_70_2JetReq);
  eeErr60_70_2JetReq.push_back(ffDiffFromeeErr60_70_2JetReq);
  float eeerr60_70_2JetReq = AddInQuadrature(eeErr60_70_2JetReq,1,4);

  float egStatErr60_70 = AddInQuadrature(statErreg,13,14);
  float egNormErr60_70 = AddInQuadrature(normErreg,13,14);
  vector<float> egErr60_70;
  egErr60_70.push_back(egStatErr60_70);
  egErr60_70.push_back(egNormErr60_70);
  float egerr60_70 = AddInQuadrature(egErr60_70,1,2);
  float QcdFFerr60_70=sqrt(egerr60_70*egerr60_70+fferr60_70*fferr60_70);
  float QcdEEerr60_70=sqrt(egerr60_70*egerr60_70+eeerr60_70*eeerr60_70);

  float egStatErr60_70_JetReq = AddInQuadrature(statErreg_JetReq,13,14);
  float egNormErr60_70_JetReq = AddInQuadrature(normErreg_JetReq,13,14);
  vector<float> egErr60_70_JetReq;
  egErr60_70_JetReq.push_back(egStatErr60_70_JetReq);
  egErr60_70_JetReq.push_back(egNormErr60_70_JetReq);
  float egerr60_70_JetReq = AddInQuadrature(egErr60_70_JetReq,1,2);
  float QcdFFerr60_70_JetReq=sqrt(egerr60_70_JetReq*egerr60_70_JetReq+fferr60_70_JetReq*fferr60_70_JetReq);
  float QcdEEerr60_70_JetReq=sqrt(egerr60_70_JetReq*egerr60_70_JetReq+eeerr60_70_JetReq*eeerr60_70_JetReq);

  float egStatErr60_70_2JetReq = AddInQuadrature(statErreg_2JetReq,13,14);
  float egNormErr60_70_2JetReq = AddInQuadrature(normErreg_2JetReq,13,14);
  vector<float> egErr60_70_2JetReq;
  egErr60_70_2JetReq.push_back(egStatErr60_70_2JetReq);
  egErr60_70_2JetReq.push_back(egNormErr60_70_2JetReq);
  float egerr60_70_2JetReq = AddInQuadrature(egErr60_70_2JetReq,1,2);
  float QcdFFerr60_70_2JetReq=sqrt(egerr60_70_2JetReq*egerr60_70_2JetReq+fferr60_70_2JetReq*fferr60_70_2JetReq);
  float QcdEEerr60_70_2JetReq=sqrt(egerr60_70_2JetReq*egerr60_70_2JetReq+eeerr60_70_2JetReq*eeerr60_70_2JetReq);


  float ffReweightErr70_80 = AddInQuadrature(reweightErrff,15,16);
  float ffStatErr70_80 = AddInQuadrature(statErrff,15,16);
  float ffNormErr70_80 = AddInQuadrature(normErrff,15,16);
  float ffDiffFromeeErr70_80 = fabs(ff70_80 - ee70_80);
  vector<float> ffErr70_80;
  ffErr70_80.push_back(ffReweightErr70_80);
  ffErr70_80.push_back(ffStatErr70_80);
  ffErr70_80.push_back(ffNormErr70_80);
  ffErr70_80.push_back(ffDiffFromeeErr70_80);
  float fferr70_80 = AddInQuadrature(ffErr70_80,1,4);

  float ffReweightErr70_80_JetReq = AddInQuadrature(reweightErrff_JetReq,15,16);
  float ffStatErr70_80_JetReq = AddInQuadrature(statErrff_JetReq,15,16);
  float ffNormErr70_80_JetReq = AddInQuadrature(normErrff_JetReq,15,16);
  float ffDiffFromeeErr70_80_JetReq = fabs(ff70_80_JetReq - ee70_80_JetReq);
  vector<float> ffErr70_80_JetReq;
  ffErr70_80_JetReq.push_back(ffReweightErr70_80_JetReq);
  ffErr70_80_JetReq.push_back(ffStatErr70_80_JetReq);
  ffErr70_80_JetReq.push_back(ffNormErr70_80_JetReq);
  ffErr70_80_JetReq.push_back(ffDiffFromeeErr70_80_JetReq);
  float fferr70_80_JetReq = AddInQuadrature(ffErr70_80_JetReq,1,4);

  float ffReweightErr70_80_2JetReq = AddInQuadrature(reweightErrff_2JetReq,15,16);
  float ffStatErr70_80_2JetReq = AddInQuadrature(statErrff_2JetReq,15,16);
  float ffNormErr70_80_2JetReq = AddInQuadrature(normErrff_2JetReq,15,16);
  float ffDiffFromeeErr70_80_2JetReq = fabs(ff70_80_2JetReq - ee70_80_2JetReq);
  vector<float> ffErr70_80_2JetReq;
  ffErr70_80_2JetReq.push_back(ffReweightErr70_80_2JetReq);
  ffErr70_80_2JetReq.push_back(ffStatErr70_80_2JetReq);
  ffErr70_80_2JetReq.push_back(ffNormErr70_80_2JetReq);
  ffErr70_80_2JetReq.push_back(ffDiffFromeeErr70_80_2JetReq);
  float fferr70_80_2JetReq = AddInQuadrature(ffErr70_80_2JetReq,1,4);

  float eeReweightErr70_80 = AddInQuadrature(reweightErree,15,16);
  float eeStatErr70_80 = AddInQuadrature(statErree,15,16);
  float eeNormErr70_80 = AddInQuadrature(normErree,15,16);
  float eeDiffFromffErr70_80 = fabs(ff70_80 - ee70_80);
  vector<float> eeErr70_80;
  eeErr70_80.push_back(eeReweightErr70_80);
  eeErr70_80.push_back(eeStatErr70_80);
  eeErr70_80.push_back(eeNormErr70_80);
  eeErr70_80.push_back(ffDiffFromeeErr70_80);
  float eeerr70_80 = AddInQuadrature(eeErr70_80,1,4);

  float eeReweightErr70_80_JetReq = AddInQuadrature(reweightErree_JetReq,15,16);
  float eeStatErr70_80_JetReq = AddInQuadrature(statErree_JetReq,15,16);
  float eeNormErr70_80_JetReq = AddInQuadrature(normErree_JetReq,15,16);
  float eeDiffFromffErr70_80_JetReq = fabs(ff70_80_JetReq - ee70_80_JetReq);
  vector<float> eeErr70_80_JetReq;
  eeErr70_80_JetReq.push_back(eeReweightErr70_80_JetReq);
  eeErr70_80_JetReq.push_back(eeStatErr70_80_JetReq);
  eeErr70_80_JetReq.push_back(eeNormErr70_80_JetReq);
  eeErr70_80_JetReq.push_back(ffDiffFromeeErr70_80_JetReq);
  float eeerr70_80_JetReq = AddInQuadrature(eeErr70_80_JetReq,1,4);

  float eeReweightErr70_80_2JetReq = AddInQuadrature(reweightErree_2JetReq,15,16);
  float eeStatErr70_80_2JetReq = AddInQuadrature(statErree_2JetReq,15,16);
  float eeNormErr70_80_2JetReq = AddInQuadrature(normErree_2JetReq,15,16);
  float eeDiffFromffErr70_80_2JetReq = fabs(ff70_80_2JetReq - ee70_80_2JetReq);
  vector<float> eeErr70_80_2JetReq;
  eeErr70_80_2JetReq.push_back(eeReweightErr70_80_2JetReq);
  eeErr70_80_2JetReq.push_back(eeStatErr70_80_2JetReq);
  eeErr70_80_2JetReq.push_back(eeNormErr70_80_2JetReq);
  eeErr70_80_2JetReq.push_back(ffDiffFromeeErr70_80_2JetReq);
  float eeerr70_80_2JetReq = AddInQuadrature(eeErr70_80_2JetReq,1,4);

  float egStatErr70_80 = AddInQuadrature(statErreg,15,16);
  float egNormErr70_80 = AddInQuadrature(normErreg,15,16);
  vector<float> egErr70_80;
  egErr70_80.push_back(egStatErr70_80);
  egErr70_80.push_back(egNormErr70_80);
  float egerr70_80 = AddInQuadrature(egErr70_80,1,2);
  float QcdFFerr70_80=sqrt(egerr70_80*egerr70_80+fferr70_80*fferr70_80);
  float QcdEEerr70_80=sqrt(egerr70_80*egerr70_80+eeerr70_80*eeerr70_80);

  float egStatErr70_80_JetReq = AddInQuadrature(statErreg_JetReq,15,16);
  float egNormErr70_80_JetReq = AddInQuadrature(normErreg_JetReq,15,16);
  vector<float> egErr70_80_JetReq;
  egErr70_80_JetReq.push_back(egStatErr70_80_JetReq);
  egErr70_80_JetReq.push_back(egNormErr70_80_JetReq);
  float egerr70_80_JetReq = AddInQuadrature(egErr70_80_JetReq,1,2);
  float QcdFFerr70_80_JetReq=sqrt(egerr70_80_JetReq*egerr70_80_JetReq+fferr70_80_JetReq*fferr70_80_JetReq);
  float QcdEEerr70_80_JetReq=sqrt(egerr70_80_JetReq*egerr70_80_JetReq+eeerr70_80_JetReq*eeerr70_80_JetReq);

  float egStatErr70_80_2JetReq = AddInQuadrature(statErreg_2JetReq,15,16);
  float egNormErr70_80_2JetReq = AddInQuadrature(normErreg_2JetReq,15,16);
  vector<float> egErr70_80_2JetReq;
  egErr70_80_2JetReq.push_back(egStatErr70_80_2JetReq);
  egErr70_80_2JetReq.push_back(egNormErr70_80_2JetReq);
  float egerr70_80_2JetReq = AddInQuadrature(egErr70_80_2JetReq,1,2);
  float QcdFFerr70_80_2JetReq=sqrt(egerr70_80_2JetReq*egerr70_80_2JetReq+fferr70_80_2JetReq*fferr70_80_2JetReq);
  float QcdEEerr70_80_2JetReq=sqrt(egerr70_80_2JetReq*egerr70_80_2JetReq+eeerr70_80_2JetReq*eeerr70_80_2JetReq);

  float ffReweightErr80_100 = AddInQuadrature(reweightErrff,17,20);
  float ffStatErr80_100 = AddInQuadrature(statErrff,17,20);
  float ffNormErr80_100 = AddInQuadrature(normErrff,17,20);
  float ffDiffFromeeErr80_100 = fabs(ff80_100 - ee80_100);
  vector<float> ffErr80_100;
  ffErr80_100.push_back(ffReweightErr80_100);
  ffErr80_100.push_back(ffStatErr80_100);
  ffErr80_100.push_back(ffNormErr80_100);
  ffErr80_100.push_back(ffDiffFromeeErr80_100);
  float fferr80_100 = AddInQuadrature(ffErr80_100,1,4);

  float ffReweightErr80_100_JetReq = AddInQuadrature(reweightErrff_JetReq,17,20);
  float ffStatErr80_100_JetReq = AddInQuadrature(statErrff_JetReq,17,20);
  float ffNormErr80_100_JetReq = AddInQuadrature(normErrff_JetReq,17,20);
  float ffDiffFromeeErr80_100_JetReq = fabs(ff80_100_JetReq - ee80_100_JetReq);
  vector<float> ffErr80_100_JetReq;
  ffErr80_100_JetReq.push_back(ffReweightErr80_100_JetReq);
  ffErr80_100_JetReq.push_back(ffStatErr80_100_JetReq);
  ffErr80_100_JetReq.push_back(ffNormErr80_100_JetReq);
  ffErr80_100_JetReq.push_back(ffDiffFromeeErr80_100_JetReq);
  float fferr80_100_JetReq = AddInQuadrature(ffErr80_100_JetReq,1,4);

  float ffReweightErr80_100_2JetReq = AddInQuadrature(reweightErrff_2JetReq,17,20);
  float ffStatErr80_100_2JetReq = AddInQuadrature(statErrff_2JetReq,17,20);
  float ffNormErr80_100_2JetReq = AddInQuadrature(normErrff_2JetReq,17,20);
  float ffDiffFromeeErr80_100_2JetReq = fabs(ff80_100_2JetReq - ee80_100_2JetReq);
  vector<float> ffErr80_100_2JetReq;
  ffErr80_100_2JetReq.push_back(ffReweightErr80_100_2JetReq);
  ffErr80_100_2JetReq.push_back(ffStatErr80_100_2JetReq);
  ffErr80_100_2JetReq.push_back(ffNormErr80_100_2JetReq);
  ffErr80_100_2JetReq.push_back(ffDiffFromeeErr80_100_2JetReq);
  float fferr80_100_2JetReq = AddInQuadrature(ffErr80_100_2JetReq,1,4);

  float eeReweightErr80_100 = AddInQuadrature(reweightErree,17,20);
  float eeStatErr80_100 = AddInQuadrature(statErree,17,20);
  float eeNormErr80_100 = AddInQuadrature(normErree,17,20);
  float eeDiffFromffErr80_100 = fabs(ff80_100 - ee80_100);
  vector<float> eeErr80_100;
  eeErr80_100.push_back(eeReweightErr80_100);
  eeErr80_100.push_back(eeStatErr80_100);
  eeErr80_100.push_back(eeNormErr80_100);
  eeErr80_100.push_back(ffDiffFromeeErr80_100);
  float eeerr80_100 = AddInQuadrature(eeErr80_100,1,4);

  float eeReweightErr80_100_JetReq = AddInQuadrature(reweightErree_JetReq,17,20);
  float eeStatErr80_100_JetReq = AddInQuadrature(statErree_JetReq,17,20);
  float eeNormErr80_100_JetReq = AddInQuadrature(normErree_JetReq,17,20);
  float eeDiffFromffErr80_100_JetReq = fabs(ff80_100_JetReq - ee80_100_JetReq);
  vector<float> eeErr80_100_JetReq;
  eeErr80_100_JetReq.push_back(eeReweightErr80_100_JetReq);
  eeErr80_100_JetReq.push_back(eeStatErr80_100_JetReq);
  eeErr80_100_JetReq.push_back(eeNormErr80_100_JetReq);
  eeErr80_100_JetReq.push_back(ffDiffFromeeErr80_100_JetReq);
  float eeerr80_100_JetReq = AddInQuadrature(eeErr80_100_JetReq,1,4);

  float eeReweightErr80_100_2JetReq = AddInQuadrature(reweightErree_2JetReq,17,20);
  float eeStatErr80_100_2JetReq = AddInQuadrature(statErree_2JetReq,17,20);
  float eeNormErr80_100_2JetReq = AddInQuadrature(normErree_2JetReq,17,20);
  float eeDiffFromffErr80_100_2JetReq = fabs(ff80_100_2JetReq - ee80_100_2JetReq);
  vector<float> eeErr80_100_2JetReq;
  eeErr80_100_2JetReq.push_back(eeReweightErr80_100_2JetReq);
  eeErr80_100_2JetReq.push_back(eeStatErr80_100_2JetReq);
  eeErr80_100_2JetReq.push_back(eeNormErr80_100_2JetReq);
  eeErr80_100_2JetReq.push_back(ffDiffFromeeErr80_100_2JetReq);
  float eeerr80_100_2JetReq = AddInQuadrature(eeErr80_100_2JetReq,1,4);

  float egStatErr80_100 = AddInQuadrature(statErreg,17,20);
  float egNormErr80_100 = AddInQuadrature(normErreg,17,20);
  vector<float> egErr80_100;
  egErr80_100.push_back(egStatErr80_100);
  egErr80_100.push_back(egNormErr80_100);
  float egerr80_100 = AddInQuadrature(egErr80_100,1,2);
  float QcdFFerr80_100=sqrt(egerr80_100*egerr80_100+fferr80_100*fferr80_100);
  float QcdEEerr80_100=sqrt(egerr80_100*egerr80_100+eeerr80_100*eeerr80_100);

  float egStatErr80_100_JetReq = AddInQuadrature(statErreg_JetReq,17,20);
  float egNormErr80_100_JetReq = AddInQuadrature(normErreg_JetReq,17,20);
  vector<float> egErr80_100_JetReq;
  egErr80_100_JetReq.push_back(egStatErr80_100_JetReq);
  egErr80_100_JetReq.push_back(egNormErr80_100_JetReq);
  float egerr80_100_JetReq = AddInQuadrature(egErr80_100_JetReq,1,2);
  float QcdFFerr80_100_JetReq=sqrt(egerr80_100_JetReq*egerr80_100_JetReq+fferr80_100_JetReq*fferr80_100_JetReq);
  float QcdEEerr80_100_JetReq=sqrt(egerr80_100_JetReq*egerr80_100_JetReq+eeerr80_100_JetReq*eeerr80_100_JetReq);

  float egStatErr80_100_2JetReq = AddInQuadrature(statErreg_2JetReq,17,20);
  float egNormErr80_100_2JetReq = AddInQuadrature(normErreg_2JetReq,17,20);
  vector<float> egErr80_100_2JetReq;
  egErr80_100_2JetReq.push_back(egStatErr80_100_2JetReq);
  egErr80_100_2JetReq.push_back(egNormErr80_100_2JetReq);
  float egerr80_100_2JetReq = AddInQuadrature(egErr80_100_2JetReq,1,2);
  float QcdFFerr80_100_2JetReq=sqrt(egerr80_100_2JetReq*egerr80_100_2JetReq+fferr80_100_2JetReq*fferr80_100_2JetReq);
  float QcdEEerr80_100_2JetReq=sqrt(egerr80_100_2JetReq*egerr80_100_2JetReq+eeerr80_100_2JetReq*eeerr80_100_2JetReq);

  float ffReweightErr100up = AddInQuadrature(reweightErrff,21,-1);
  float ffStatErr100up = AddInQuadrature(statErrff,21,-1);
  float ffNormErr100up = AddInQuadrature(normErrff,21,-1);
  float ffDiffFromeeErr100up = fabs(ff100up - ee100up);//AddUp(diffFromeeErrorff,21,-1);
  cout<<endl;
  cout <<"ffReweightErr100up   : "<<ffReweightErr100up<<endl;
  cout <<"ffStatErr100up       : "<<ffStatErr100up<<endl;
  cout <<"ffNormErr100up       : "<<ffNormErr100up<<endl;
  cout <<"ffDiffFromeeErr100up : "<<ffDiffFromeeErr100up<<endl;
  vector<float> ffErr100up;
  ffErr100up.push_back(ffReweightErr100up);
  ffErr100up.push_back(ffStatErr100up);
  ffErr100up.push_back(ffNormErr100up);
  ffErr100up.push_back(ffDiffFromeeErr100up);
  float fferr100up = AddInQuadrature(ffErr100up,1,4);


  float ffReweightErr100up_JetReq = AddInQuadrature(reweightErrff_JetReq,21,-1);
  float ffStatErr100up_JetReq = AddInQuadrature(statErrff_JetReq,21,-1);
  float ffNormErr100up_JetReq = AddInQuadrature(normErrff_JetReq,21,-1);
  float ffDiffFromeeErr100up_JetReq = fabs(ff100up_JetReq - ee100up_JetReq);//AddUp(diffFromeeErrorff_JetReq,21,-1);
  cout <<"ffReweightErr100up_JetReq   : "<<ffReweightErr100up_JetReq<<endl;
  cout <<"ffStatErr100up_JetReq       : "<<ffStatErr100up_JetReq<<endl;
  cout <<"ffNormErr100up_JetReq       : "<<ffNormErr100up_JetReq<<endl;
  cout <<"ffDiffFromeeErr100up_JetReq : "<<ffDiffFromeeErr100up_JetReq<<endl;
  vector<float> ffErr100up_JetReq;
  ffErr100up_JetReq.push_back(ffReweightErr100up_JetReq);
  ffErr100up_JetReq.push_back(ffStatErr100up_JetReq);
  ffErr100up_JetReq.push_back(ffNormErr100up_JetReq);
  ffErr100up_JetReq.push_back(ffDiffFromeeErr100up_JetReq);
  float fferr100up_JetReq = AddInQuadrature(ffErr100up_JetReq,1,4);

  float ffReweightErr100up_2JetReq = AddInQuadrature(reweightErrff_2JetReq,21,-1);
  float ffStatErr100up_2JetReq = AddInQuadrature(statErrff_2JetReq,21,-1);
  float ffNormErr100up_2JetReq = AddInQuadrature(normErrff_2JetReq,21,-1);
  float ffDiffFromeeErr100up_2JetReq = fabs(ff100up_2JetReq - ee100up_2JetReq);//AddUp(diffFromeeErrorff_2JetReq,21,-1);
  cout <<"ffReweightErr100up_2JetReq   : "<<ffReweightErr100up_2JetReq<<endl;
  cout <<"ffStatErr100up_2JetReq       : "<<ffStatErr100up_2JetReq<<endl;
  cout <<"ffNormErr100up_2JetReq       : "<<ffNormErr100up_2JetReq<<endl;
  cout <<"ffDiffFromeeErr100up_2JetReq : "<<ffDiffFromeeErr100up_2JetReq<<endl;
  vector<float> ffErr100up_2JetReq;
  ffErr100up_2JetReq.push_back(ffReweightErr100up_2JetReq);
  ffErr100up_2JetReq.push_back(ffStatErr100up_2JetReq);
  ffErr100up_2JetReq.push_back(ffNormErr100up_2JetReq);
  ffErr100up_2JetReq.push_back(ffDiffFromeeErr100up_2JetReq);
  float fferr100up_2JetReq = AddInQuadrature(ffErr100up_2JetReq,1,4);

  float eeReweightErr100up = AddInQuadrature(reweightErree,21,-1);
  float eeStatErr100up = AddInQuadrature(statErree,21,-1);
  float eeNormErr100up = AddInQuadrature(normErree,21,-1);
  float eeDiffFromffErr100up = fabs(ff100up - ee100up);//AddUp(diffFromffErroree,21,-1);
  cout<<endl;
  cout <<"eeReweightErr100up   : "<<eeReweightErr100up<<endl;
  cout <<"eeStatErr100up       : "<<eeStatErr100up<<endl;
  cout <<"eeNormErr100up       : "<<eeNormErr100up<<endl;
  vector<float> eeErr100up;
  eeErr100up.push_back(eeReweightErr100up);
  eeErr100up.push_back(eeStatErr100up);
  eeErr100up.push_back(eeNormErr100up);
  eeErr100up.push_back(ffDiffFromeeErr100up);
  float eeerr100up = AddInQuadrature(eeErr100up,1,4);

  float eeReweightErr100up_JetReq = AddInQuadrature(reweightErree_JetReq,21,-1);
  float eeStatErr100up_JetReq = AddInQuadrature(statErree_JetReq,21,-1);
  float eeNormErr100up_JetReq = AddInQuadrature(normErree_JetReq,21,-1);
  float eeDiffFromffErr100up_JetReq = fabs(ff100up_JetReq - ee100up_JetReq);//AddUp(diffFromffErroree_JetReq,21,-1);
  cout<<endl;
  cout <<"eeReweightErr100up_JetReq   : "<<eeReweightErr100up_JetReq<<endl;
  cout <<"eeStatErr100up_JetReq       : "<<eeStatErr100up_JetReq<<endl;
  cout <<"eeNormErr100up_JetReq       : "<<eeNormErr100up_JetReq<<endl;
  vector<float> eeErr100up_JetReq;
  eeErr100up_JetReq.push_back(eeReweightErr100up_JetReq);
  eeErr100up_JetReq.push_back(eeStatErr100up_JetReq);
  eeErr100up_JetReq.push_back(eeNormErr100up_JetReq);
  eeErr100up_JetReq.push_back(ffDiffFromeeErr100up_JetReq);
  float eeerr100up_JetReq = AddInQuadrature(eeErr100up_JetReq,1,4);

  float eeReweightErr100up_2JetReq = AddInQuadrature(reweightErree_2JetReq,21,-1);
  float eeStatErr100up_2JetReq = AddInQuadrature(statErree_2JetReq,21,-1);
  float eeNormErr100up_2JetReq = AddInQuadrature(normErree_2JetReq,21,-1);
  float eeDiffFromffErr100up_2JetReq = fabs(ff100up_2JetReq - ee100up_2JetReq);//AddUp(diffFromffErroree_2JetReq,21,-1);
  cout<<endl;
  cout <<"eeReweightErr100up_2JetReq   : "<<eeReweightErr100up_2JetReq<<endl;
  cout <<"eeStatErr100up_2JetReq       : "<<eeStatErr100up_2JetReq<<endl;
  cout <<"eeNormErr100up_2JetReq       : "<<eeNormErr100up_2JetReq<<endl;
  vector<float> eeErr100up_2JetReq;
  eeErr100up_2JetReq.push_back(eeReweightErr100up_2JetReq);
  eeErr100up_2JetReq.push_back(eeStatErr100up_2JetReq);
  eeErr100up_2JetReq.push_back(eeNormErr100up_2JetReq);
  eeErr100up_2JetReq.push_back(ffDiffFromeeErr100up_2JetReq);
  float eeerr100up_2JetReq = AddInQuadrature(eeErr100up_2JetReq,1,4);

  float egStatErr100up = AddInQuadrature(statErreg,21,-1);
  float egNormErr100up = AddInQuadrature(normErreg,21,-1);
  cout<<endl;
  cout <<"egStatErr100up       : "<<egStatErr100up<<endl;
  cout <<"egNormErr100up       : "<<egNormErr100up<<endl;
  vector<float> egErr100up;
  egErr100up.push_back(egStatErr100up);
  egErr100up.push_back(egNormErr100up);
  float egerr100up = AddInQuadrature(egErr100up,1,2);
  float QcdFFerr100up=sqrt(egerr100up*egerr100up+fferr100up*fferr100up);
  float QcdEEerr100up=sqrt(egerr100up*egerr100up+eeerr100up*eeerr100up);

  float egStatErr100up_JetReq = AddInQuadrature(statErreg_JetReq,21,-1);
  float egNormErr100up_JetReq = AddInQuadrature(normErreg_JetReq,21,-1);
  cout<<endl;
  cout <<"egStatErr100up_JetReq       : "<<egStatErr100up_JetReq<<endl;
  cout <<"egNormErr100up_JetReq       : "<<egNormErr100up_JetReq<<endl;
  vector<float> egErr100up_JetReq;
  egErr100up_JetReq.push_back(egStatErr100up_JetReq);
  egErr100up_JetReq.push_back(egNormErr100up_JetReq);
  float egerr100up_JetReq = AddInQuadrature(egErr100up_JetReq,1,2);
  float QcdFFerr100up_JetReq=sqrt(egerr100up_JetReq*egerr100up_JetReq+fferr100up_JetReq*fferr100up_JetReq);
  float QcdEEerr100up_JetReq=sqrt(egerr100up_JetReq*egerr100up_JetReq+eeerr100up_JetReq*eeerr100up_JetReq);

  float egStatErr100up_2JetReq = AddInQuadrature(statErreg_2JetReq,21,-1);
  float egNormErr100up_2JetReq = AddInQuadrature(normErreg_2JetReq,21,-1);
  cout<<endl;
  cout <<"egStatErr100up_2JetReq       : "<<egStatErr100up_2JetReq<<endl;
  cout <<"egNormErr100up_2JetReq       : "<<egNormErr100up_2JetReq<<endl;
  vector<float> egErr100up_2JetReq;
  egErr100up_2JetReq.push_back(egStatErr100up_2JetReq);
  egErr100up_2JetReq.push_back(egNormErr100up_2JetReq);
  float egerr100up_2JetReq = AddInQuadrature(egErr100up_2JetReq,1,2);
  float QcdFFerr100up_2JetReq=sqrt(egerr100up_2JetReq*egerr100up_2JetReq+fferr100up_2JetReq*fferr100up_2JetReq);
  float QcdEEerr100up_2JetReq=sqrt(egerr100up_2JetReq*egerr100up_2JetReq+eeerr100up_2JetReq*eeerr100up_2JetReq);

  cout<<"\\documentclass[11pt]{article}"<<endl;
  cout<<"\\usepackage{calc}"<<endl;
  cout<<"\\usepackage{multirow}"<<endl;
  cout<<"\\usepackage{verbatim}"<<endl;
  cout<<"\\usepackage{changepage}"<<endl;
  cout<<"\\usepackage{tabularx}"<<endl;
  cout<<"\\begin{document}"<<endl;

  cout<<"begin{table}[tpb]\\label{tbl:metcounts}"<<endl;
  cout<<" \\begin{center}"<<endl;
  cout<<"  \\begin{tiny}"<<endl;
  cout<<"   \\begin{tabularx}{1.0\\textwidth}{ | l | c | c | c | c | c | c |}\n      \\hline"<<endl;
  cout<<"%\n%\n%\n";
  cout<<"      \\bf{Type} & \\bf{no jet cut} & \\bf{$\\ge$1 jet}&\\multicolumn{4}{c|}{Individual Errors} \\\\ \\hline"<<endl;
  cout<<"      %met<20"<<endl;
  if(doffgfcomb){
    cout<<"       & \\bf{$<$20} & \\bf{$<$20} & \\bf{Stat} & \\bf{Norm} & \\bf{Reweight} & \\bf{$\\Delta$(comb,ff)} \\\\ \\hline"<<endl;
    printf("      $\\gamma\\gamma$ & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f &  & & \\\\ \\hline\n",gg0_20,gg0_20Error,gg0_20_JetReq,gg0_20Error_JetReq,gg0_20Error_JetReq);
    printf("      comb QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ff0_20,fferr0_20,ff0_20_JetReq,fferr0_20_JetReq,ffStatErr0_20_JetReq,ffNormErr0_20_JetReq,ffReweightErr0_20_JetReq,ffDiffFromeeErr0_20_JetReq);
    printf("      ff QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ee0_20,eeerr0_20,ee0_20_JetReq,eeerr0_20_JetReq,eeStatErr0_20_JetReq,eeNormErr0_20_JetReq,eeReweightErr0_20_JetReq,eeDiffFromffErr0_20_JetReq);
    printf("      EWK background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & & \\\\ \\hline\n",eg0_20,egerr0_20,eg0_20_JetReq,egerr0_20_JetReq,egStatErr0_20_JetReq,egNormErr0_20_JetReq);
    printf("      Total background (comb) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDff0_20,QcdFFerr0_20,QCDff0_20_JetReq,QcdFFerr0_20_JetReq);
    printf("      Total background (ff) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDee0_20,QcdEEerr0_20,QCDee0_20_JetReq,QcdEEerr0_20_JetReq);
  }
  else{
    cout<<"      \\bf{Type} & \\bf{no jet cut} & \\bf{$\\ge$1 jet} & \\bf{$\\ge$1 jet} & \\bf{$\\ge$1 jet} & \\bf{$\\ge$1 jet} & \\bf{$\\ge$1 jet} \\\\ \\hline"<<endl;
    cout<<"       & \\bf{$<$20} & \\bf{$<$20} & \\bf{Stat} & \\bf{Norm} & \\bf{Reweight} & \\bf{$\\Delta$(ee,ff)} \\\\ \\hline"<<endl;
    printf("      $\\gamma\\gamma$ & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f &  & & \\\\ \\hline\n",gg0_20,gg0_20Error,gg0_20_JetReq,gg0_20Error_JetReq,gg0_20Error_JetReq);
    printf("      ff QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ff0_20,fferr0_20,ff0_20_JetReq,fferr0_20_JetReq,ffStatErr0_20_JetReq,ffNormErr0_20_JetReq,ffReweightErr0_20_JetReq,ffDiffFromeeErr0_20_JetReq);
    printf("      ee QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ee0_20,eeerr0_20,ee0_20_JetReq,eeerr0_20_JetReq,eeStatErr0_20_JetReq,eeNormErr0_20_JetReq,eeReweightErr0_20_JetReq,eeDiffFromffErr0_20_JetReq);
    printf("      EWK background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & & \\\\ \\hline\n",eg0_20,egerr0_20,eg0_20_JetReq,egerr0_20_JetReq,egStatErr0_20_JetReq,egNormErr0_20_JetReq);
    printf("      Total background (ff) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDff0_20,QcdFFerr0_20,QCDff0_20_JetReq,QcdFFerr0_20_JetReq);
    printf("      Total background (ee) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDee0_20,QcdEEerr0_20,QCDee0_20_JetReq,QcdEEerr0_20_JetReq);
  }
  //
  if(doffgfcomb){
    cout<<"       %30<met<50"<<endl;
    cout<<"       & \\bf{30-50} & \\bf{30-50} & \\bf{Stat} & \\bf{Norm} & \\bf{Reweight} & \\bf{$\\Delta$(comb,ff)} \\\\ \\hline"<<endl;
    printf("      $\\gamma\\gamma$ & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f &  & & \\\\ \\hline\n",gg30_50,gg30_50Error,gg30_50_JetReq,gg30_50Error_JetReq,gg30_50Error_JetReq);
    printf("      comb QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ff30_50,fferr30_50,ff30_50_JetReq,fferr30_50_JetReq,ffStatErr30_50_JetReq,ffNormErr30_50_JetReq,ffReweightErr30_50_JetReq,ffDiffFromeeErr30_50_JetReq);
    printf("      ff QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ee30_50,eeerr30_50,ee30_50_JetReq,eeerr30_50_JetReq,eeStatErr30_50_JetReq,eeNormErr30_50_JetReq,eeReweightErr30_50_JetReq,eeDiffFromffErr30_50_JetReq);
    printf("      EWK background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & & \\\\ \\hline\n",eg30_50,egerr30_50,eg30_50_JetReq,egerr30_50_JetReq,egStatErr30_50_JetReq,egNormErr30_50_JetReq);
    printf("      Total background (comb) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDff30_50,QcdFFerr30_50,QCDff30_50_JetReq,QcdFFerr30_50_JetReq);
    printf("      Total background (ff) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDee30_50,QcdEEerr30_50,QCDee30_50_JetReq,QcdEEerr30_50_JetReq);
  }
  else{
    cout<<"       %30<met<50"<<endl;
    cout<<"       & \\bf{30-50} & \\bf{30-50} & \\bf{Stat} & \\bf{Norm} & \\bf{Reweight} & \\bf{$\\Delta$(ee,ff)} \\\\ \\hline"<<endl;
    printf("      $\\gamma\\gamma$ & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f &  & & \\\\ \\hline\n",gg30_50,gg30_50Error,gg30_50_JetReq,gg30_50Error_JetReq,gg30_50Error_JetReq);
    printf("      ff QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ff30_50,fferr30_50,ff30_50_JetReq,fferr30_50_JetReq,ffStatErr30_50_JetReq,ffNormErr30_50_JetReq,ffReweightErr30_50_JetReq,ffDiffFromeeErr30_50_JetReq);
    printf("      ee QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ee30_50,eeerr30_50,ee30_50_JetReq,eeerr30_50_JetReq,eeStatErr30_50_JetReq,eeNormErr30_50_JetReq,eeReweightErr30_50_JetReq,eeDiffFromffErr30_50_JetReq);
    printf("      EWK background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & & \\\\ \\hline\n",eg30_50,egerr30_50,eg30_50_JetReq,egerr30_50_JetReq,egStatErr30_50_JetReq,egNormErr30_50_JetReq);
    printf("      Total background (ff) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDff30_50,QcdFFerr30_50,QCDff30_50_JetReq,QcdFFerr30_50_JetReq);
    printf("      Total background (ee) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDee30_50,QcdEEerr30_50,QCDee30_50_JetReq,QcdEEerr30_50_JetReq);
  }
  //
  cout<<"       %50<met<60"<<endl;
  cout<<"        & \\bf{50-60} & \\bf{50-60} & \\bf{Stat} & \\bf{Norm} & \\bf{Reweight} & \\bf{$\\Delta$(ee,ff)} \\\\ \\hline"<<endl;
  printf("      $\\gamma\\gamma$ & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f &  & & \\\\ \\hline\n",gg50_60,gg50_60Error,gg50_60_JetReq,gg50_60Error_JetReq,gg50_60Error_JetReq);
  printf("      ff QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ff50_60,fferr50_60,ff50_60_JetReq,fferr50_60_JetReq,ffStatErr50_60_JetReq,ffNormErr50_60_JetReq,ffReweightErr50_60_JetReq,ffDiffFromeeErr50_60_JetReq);
  printf("      ee QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ee50_60,eeerr50_60,ee50_60_JetReq,eeerr50_60_JetReq,eeStatErr50_60_JetReq,eeNormErr50_60_JetReq,eeReweightErr50_60_JetReq,eeDiffFromffErr50_60_JetReq);
  printf("      EWK background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & & \\\\ \\hline\n",eg50_60,egerr50_60,eg50_60_JetReq,egerr50_60_JetReq,egStatErr50_60_JetReq,egNormErr50_60_JetReq);
  printf("      Total background (ff) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDff50_60,QcdFFerr50_60,QCDff50_60_JetReq,QcdFFerr50_60_JetReq);
  printf("      Total background (ee) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDee50_60,QcdEEerr50_60,QCDee50_60_JetReq,QcdEEerr50_60_JetReq);
  //
  cout<<"       %60<met<70"<<endl;
  cout<<"        & \\bf{60-70} & \\bf{60-70} & \\bf{Stat} & \\bf{Norm} & \\bf{Reweight} & \\bf{$\\Delta$(ee,ff)} \\\\ \\hline"<<endl;
  printf("      $\\gamma\\gamma$ & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f &  & & \\\\ \\hline\n",gg60_70,gg60_70Error,gg60_70_JetReq,gg60_70Error_JetReq,gg60_70Error_JetReq);
  printf("      ff QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ff60_70,fferr60_70,ff60_70_JetReq,fferr60_70_JetReq,ffStatErr60_70_JetReq,ffNormErr60_70_JetReq,ffReweightErr60_70_JetReq,ffDiffFromeeErr60_70_JetReq);
  printf("      ee QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ee60_70,eeerr60_70,ee60_70_JetReq,eeerr60_70_JetReq,eeStatErr60_70_JetReq,eeNormErr60_70_JetReq,eeReweightErr60_70_JetReq,eeDiffFromffErr60_70_JetReq);
  printf("      EWK background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & & \\\\ \\hline\n",eg60_70,egerr60_70,eg60_70_JetReq,egerr60_70_JetReq,egStatErr60_70_JetReq,egNormErr60_70_JetReq);
  printf("      Total background (ff) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDff60_70,QcdFFerr60_70,QCDff60_70_JetReq,QcdFFerr60_70_JetReq);
  printf("      Total background (ee) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDee60_70,QcdEEerr60_70,QCDee60_70_JetReq,QcdEEerr60_70_JetReq);
  //
  cout<<"       %70<met<80"<<endl;
  cout<<"        & \\bf{70-80} & \\bf{70-80} & \\bf{Stat} & \\bf{Norm} & \\bf{Reweight} & \\bf{$\\Delta$(ee,ff)} \\\\ \\hline"<<endl;
  printf("      $\\gamma\\gamma$ & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f &  & & \\\\ \\hline\n",gg70_80,gg70_80Error,gg70_80_JetReq,gg70_80Error_JetReq,gg70_80Error_JetReq);
  printf("      ff QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ff70_80,fferr70_80,ff70_80_JetReq,fferr70_80_JetReq,ffStatErr70_80_JetReq,ffNormErr70_80_JetReq,ffReweightErr70_80_JetReq,ffDiffFromeeErr70_80_JetReq);
  printf("      ee QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ee70_80,eeerr70_80,ee70_80_JetReq,eeerr70_80_JetReq,eeStatErr70_80_JetReq,eeNormErr70_80_JetReq,eeReweightErr70_80_JetReq,eeDiffFromffErr70_80_JetReq);
  printf("      EWK background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & & \\\\ \\hline\n",eg70_80,egerr70_80,eg70_80_JetReq,egerr70_80_JetReq,egStatErr70_80_JetReq,egNormErr70_80_JetReq);
  printf("      Total background (ff) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDff70_80,QcdFFerr70_80,QCDff70_80_JetReq,QcdFFerr70_80_JetReq);
  printf("      Total background (ee) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDee70_80,QcdEEerr70_80,QCDee70_80_JetReq,QcdEEerr70_80_JetReq);
  //
  cout<<"       %80<met<100"<<endl;
  if(doffgfcomb){
    cout<<"        & \\bf{80-100} & \\bf{80-100} & \\bf{Stat} & \\bf{Norm} & \\bf{Reweight} & \\bf{$\\Delta$(comb,ff)} \\\\ \\hline"<<endl;
    printf("      $\\gamma\\gamma$ & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f &  & & \\\\ \\hline\n",gg80_100,gg80_100Error,gg80_100_JetReq,gg80_100Error_JetReq,gg80_100Error_JetReq);
    printf("      comb QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ff80_100,fferr80_100,ff80_100_JetReq,fferr80_100_JetReq,ffStatErr80_100_JetReq,ffNormErr80_100_JetReq,ffReweightErr80_100_JetReq,ffDiffFromeeErr80_100_JetReq);
    printf("      ff QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ee80_100,eeerr80_100,ee80_100_JetReq,eeerr80_100_JetReq,eeStatErr80_100_JetReq,eeNormErr80_100_JetReq,eeReweightErr80_100_JetReq,eeDiffFromffErr80_100_JetReq);
    printf("      EWK background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & & \\\\ \\hline\n",eg80_100,egerr80_100,eg80_100_JetReq,egerr80_100_JetReq,egStatErr80_100_JetReq,egNormErr80_100_JetReq);
    printf("      Total background (comb) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDff80_100,QcdFFerr80_100,QCDff80_100_JetReq,QcdFFerr80_100_JetReq);
    printf("      Total background (ff) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDee80_100,QcdEEerr80_100,QCDee80_100_JetReq,QcdEEerr80_100_JetReq);
  }
  else{
    cout<<"        & \\bf{80-100} & \\bf{80-100} & \\bf{Stat} & \\bf{Norm} & \\bf{Reweight} & \\bf{$\\Delta$(ee,ff)} \\\\ \\hline"<<endl;
    printf("      $\\gamma\\gamma$ & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f &  & & \\\\ \\hline\n",gg80_100,gg80_100Error,gg80_100_JetReq,gg80_100Error_JetReq,gg80_100Error_JetReq);
    printf("      ff QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ff80_100,fferr80_100,ff80_100_JetReq,fferr80_100_JetReq,ffStatErr80_100_JetReq,ffNormErr80_100_JetReq,ffReweightErr80_100_JetReq,ffDiffFromeeErr80_100_JetReq);
    printf("      ee QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ee80_100,eeerr80_100,ee80_100_JetReq,eeerr80_100_JetReq,eeStatErr80_100_JetReq,eeNormErr80_100_JetReq,eeReweightErr80_100_JetReq,eeDiffFromffErr80_100_JetReq);
    printf("      EWK background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & & \\\\ \\hline\n",eg80_100,egerr80_100,eg80_100_JetReq,egerr80_100_JetReq,egStatErr80_100_JetReq,egNormErr80_100_JetReq);
    printf("      Total background (ff) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDff80_100,QcdFFerr80_100,QCDff80_100_JetReq,QcdFFerr80_100_JetReq);
    printf("      Total background (ee) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDee80_100,QcdEEerr80_100,QCDee80_100_JetReq,QcdEEerr80_100_JetReq);
  }
  //
  cout<<"       %met>50"<<endl;
  if(doffgfcomb){
    cout<<"        & \\bf{$>$50} & \\bf{$>$50} & \\bf{Stat} & \\bf{Norm} & \\bf{Reweight} & \\bf{$\\Delta$(comb,ff)} \\\\ \\hline"<<endl;
    printf("      $\\gamma\\gamma$ & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f &  & & \\\\ \\hline\n",gg50up,gg50upError,gg50up_JetReq,gg50upError_JetReq,gg50upError_JetReq);
    printf("      comb QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ff50up,fferr50up,ff50up_JetReq,fferr50up_JetReq,ffStatErr50up_JetReq,ffNormErr50up_JetReq,ffReweightErr50up_JetReq,ffDiffFromeeErr50up_JetReq);
    printf("      ff QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ee50up,eeerr50up,ee50up_JetReq,eeerr50up_JetReq,eeStatErr50up_JetReq,eeNormErr50up_JetReq,eeReweightErr50up_JetReq,eeDiffFromffErr50up_JetReq);
    printf("      EWK background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & & \\\\ \\hline\n",eg50up,egerr50up,eg50up_JetReq,egerr50up_JetReq,egStatErr50up_JetReq,egNormErr50up_JetReq);
    printf("      Total background (comb) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDff50up,QcdFFerr50up,QCDff50up_JetReq,QcdFFerr50up_JetReq);
    printf("      Total background (ff) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDee50up,QCDee50upError,QCDee50up_JetReq,QcdEEerr50up_JetReq);
  }
  else{
    cout<<"        & \\bf{$>$50} & \\bf{$>$50} & \\bf{Stat} & \\bf{Norm} & \\bf{Reweight} & \\bf{$\\Delta$(ee,ff)} \\\\ \\hline"<<endl;
    printf("      $\\gamma\\gamma$ & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f &  & & \\\\ \\hline\n",gg50up,gg50upError,gg50up_JetReq,gg50upError_JetReq,gg50upError_JetReq);
    printf("      ff QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ff50up,fferr50up,ff50up_JetReq,fferr50up_JetReq,ffStatErr50up_JetReq,ffNormErr50up_JetReq,ffReweightErr50up_JetReq,ffDiffFromeeErr50up_JetReq);
    printf("      ee QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ee50up,eeerr50up,ee50up_JetReq,eeerr50up_JetReq,eeStatErr50up_JetReq,eeNormErr50up_JetReq,eeReweightErr50up_JetReq,eeDiffFromffErr50up_JetReq);
    printf("      EWK background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & & \\\\ \\hline\n",eg50up,egerr50up,eg50up_JetReq,egerr50up_JetReq,egStatErr50up_JetReq,egNormErr50up_JetReq);
    printf("      Total background (ff) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDff50up,QcdFFerr50up,QCDff50up_JetReq,QcdFFerr50up_JetReq);
    printf("      Total background (ee) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDee50up,QCDee50upError,QCDee50up_JetReq,QcdEEerr50up_JetReq);
  }
  //
  cout<<"       %met>100"<<endl;
  if(doffgfcomb){
    cout<<"        & \\bf{$>$100} & \\bf{$>$100} & \\bf{Stat} & \\bf{Norm} & \\bf{Reweight} & \\bf{$\\Delta$(comb,ff)} \\\\ \\hline"<<endl;
    printf("      $\\gamma\\gamma$ & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f &  & & \\\\ \\hline\n",gg100up,gg100upError,gg100up_JetReq,gg100upError_JetReq,gg100upError_JetReq);
    printf("      comb QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ff100up,fferr100up,ff100up_JetReq,fferr100up_JetReq,ffStatErr100up_JetReq,ffNormErr100up_JetReq,ffReweightErr100up_JetReq,ffDiffFromeeErr100up_JetReq);
    printf("      ff QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ee100up,eeerr100up,ee100up_JetReq,eeerr100up_JetReq,eeStatErr100up_JetReq,eeNormErr100up_JetReq,eeReweightErr100up_JetReq,eeDiffFromffErr100up_JetReq);
    printf("      EWK background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & & \\\\ \\hline\n",eg100up,egerr100up,eg100up_JetReq,egerr100up_JetReq,egStatErr100up_JetReq,egNormErr100up_JetReq);
    printf("      Total background (comb) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDff100up,QcdFFerr100up,QCDff100up_JetReq,QcdFFerr100up_JetReq);
    printf("      Total background (ff) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDee100up,QcdEEerr100up,QCDee100up_JetReq,QcdEEerr100up_JetReq);
  }
  else{
    cout<<"        & \\bf{$>$100} & \\bf{$>$100} & \\bf{Stat} & \\bf{Norm} & \\bf{Reweight} & \\bf{$\\Delta$(ee,ff)} \\\\ \\hline"<<endl;
    printf("      $\\gamma\\gamma$ & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f &  & & \\\\ \\hline\n",gg100up,gg100upError,gg100up_JetReq,gg100upError_JetReq,gg100upError_JetReq);
    printf("      ff QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ff100up,fferr100up,ff100up_JetReq,fferr100up_JetReq,ffStatErr100up_JetReq,ffNormErr100up_JetReq,ffReweightErr100up_JetReq,ffDiffFromeeErr100up_JetReq);
    printf("      ee QCD background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & %4.1f & %4.1f \\\\ \\hline\n",ee100up,eeerr100up,ee100up_JetReq,eeerr100up_JetReq,eeStatErr100up_JetReq,eeNormErr100up_JetReq,eeReweightErr100up_JetReq,eeDiffFromffErr100up_JetReq);
    printf("      EWK background & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f & %4.1f & & \\\\ \\hline\n",eg100up,egerr100up,eg100up_JetReq,egerr100up_JetReq,egStatErr100up_JetReq,egNormErr100up_JetReq);
    printf("      Total background (ff) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDff100up,QcdFFerr100up,QCDff100up_JetReq,QcdFFerr100up_JetReq);
    printf("      Total background (ee) & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f &  &  & & \\\\ \\hline\n",QCDee100up,QcdEEerr100up,QCDee100up_JetReq,QcdEEerr100up_JetReq);
  }
  cout<<"%\n%\n%\n";
  cout<<"   \\end{tabularx}"<<endl;
  cout<<"  \\end{tiny}"<<endl;
  cout<<" \\end{center}"<<endl;
  cout<<"\\end{table}\n"<<endl;
  cout<<"\\end{document}"<<endl;

  fin.Close();
  // return;    
}

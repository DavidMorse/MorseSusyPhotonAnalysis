#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAttPad.h"
using namespace std;
Double_t fpow(Double_t *x, Double_t *par){
  //Fit function definition
  return par[0]*pow(x[0],par[1]);
}

void FitTest(){
  //gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","",1280,800);c1->cd();c1->Draw();
  TPad *p1 = new TPad("p1","",0.,0.,.33,1.);//p1->SetRightMargin(0);
  TPad *p2 = new TPad("p2","",.33,0.,.66,1.);//p2->SetLeftMargin(0);p2->SetRightMargin(0);
  TPad *p3 = new TPad("p3","",.66,0.,1.,1.);//p3->SetLeftMargin(0);
  p1->Draw();p2->Draw();p3->Draw();
  c1->cd();p1->cd();
  TH1F* hist = new TH1F("hist","No Rebin",200,0,100);hist->Sumw2();
  TH1F* hist2 = new TH1F("hist2","Rebin(2)",200,0,100);hist2->Sumw2();
  TH1F* hist3 = new TH1F("hist3","Rebin(4)",200,0,100);hist3->Sumw2();

  for(int i=0;i<=2*hist->GetNbinsX();i++){
    for(int j=0;j<int(500*pow(i+.1,-.08));j++){
      hist->Fill(i/2.+.4);
      hist2->Fill(i/2.+.4);
      hist3->Fill(i/2.+.4);
    }
  }
  TF1* fitCurve = new TF1("fitCurve",fpow,0,100,2);
  Double_t param1= log(48.)/log(350.);
  Double_t param0= 350./pow(48.,param1);
  fitCurve->SetParameter(0,param0);
  fitCurve->SetParameter(1,param1);

  hist->Fit(fitCurve,"LS","",0,100);
  float Int4_72_hist = hist->Integral(hist->FindBin(4),hist->FindBin(71.9));
  float Int4_72_fit = fitCurve->Integral(4,72);
  hist->Draw();
  char str[4],strFit[4];
  sprintf(str,"histo: %4.1f",Int4_72_hist);
  sprintf(strFit,"fit   : %4.1f",Int4_72_fit);
  TPaveText *text= new TPaveText(.3,.55,.84,.74,"NDC");text->SetFillStyle(0);text->SetBorderSize(0);
  text->AddText("Integral from 4-72:");
  text->AddText(str);
  text->AddText(strFit);
  text->Draw();
  p2->cd();
 
  hist2->Rebin(2);

  TF1* fitCurve2 = new TF1("fitCurve2",fpow,0,100,2);
  param1= log(50.)/log(350.);
  param0= 350./pow(48.,param1);
  fitCurve2->SetParameter(0,param0);
  fitCurve2->SetParameter(1,param1);
  hist2->Fit(fitCurve2,"LS","",0,100);
  float Int4_72_hist2 = hist2->Integral(hist2->FindBin(4),hist2->FindBin(71.9));
  float Int4_72_fit2 = fitCurve2->Integral(4,72);
  hist2->Draw();
  char str2[4],strFit2[4];
  sprintf(str2,"histo: %4.1f",Int4_72_hist2);
  sprintf(strFit2,"fit   : %4.1f",Int4_72_fit2);
  TPaveText *text2= new TPaveText(.3,.55,.84,.74,"NDC");text2->SetFillStyle(0);text2->SetBorderSize(0);
  text2->AddText("Integral from 4-72:");
  text2->AddText(str2);
  text2->AddText(strFit2);
  text2->Draw();
  p3->cd();

  
  hist3->Rebin(4);

  TF1* fitCurve3 = new TF1("fitCurve3",fpow,0,100,2);
  param1= log(48.)/log(350.);
  param0= 350./pow(48.,param1);
  fitCurve3->SetParameter(0,param0);
  fitCurve3->SetParameter(1,param1);

  hist3->Fit(fitCurve3,"LS","",0,100);  
  hist3->Draw();
  float Int4_72_hist3 = hist3->Integral(hist3->FindBin(4),hist3->FindBin(71.9));
  float Int4_72_fit3 = fitCurve3->Integral(4,72);
  char str3[4],strFit3[4];
  sprintf(str3,"histo: %4.1f",Int4_72_hist3);
  sprintf(strFit3,"fit   : %4.1f",Int4_72_fit3);
  TPaveText *text3= new TPaveText(.3,.55,.84,.74,"NDC");text3->SetFillStyle(0);text3->SetBorderSize(0);
  text3->AddText("Integral from 4-72:");
  text3->AddText(str3);
  text3->AddText(strFit3);
  text3->Draw();

  cout<<"Integral from 4-72 from:\nhisto: "<<Int4_72_hist<<endl<<"fit  : "<<Int4_72_fit<<endl;
  cout<<endl<<"Rebin(2)"<<endl;
  cout<<"Integral from 4-72 from:\nhisto: "<<Int4_72_hist2<<endl<<"fit  : "<<Int4_72_fit2<<endl;
  cout<<endl<<"Rebin(3)"<<endl;
  cout<<"Integral from 4-72 from:\nhisto: "<<Int4_72_hist3<<endl<<"fit  : "<<Int4_72_fit3<<endl;

}

#include "TRandom.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"

bool reject=true;
Double_t fpow(Double_t *x, Double_t *par){
  //Fit function definition
  if(reject && x[0]>=118 && x[0]<=133){
    TF1::RejectPoint();
    return 0;
  }
  return par[0]*pow(x[0],par[1]);
}

void fitToys2(){

  int count=0,countNotFail=0;

  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","c1",800,600);c1->cd();
  //TCanvas *c2 = new TCanvas("c2","c2",800,600);c2->cd();
  //c1->cd();

  TF1* fit = new TF1("fit","[0]*pow(x,[1])",103,163);
  fit->SetParameters(1.28542e+08,-3.90422);
  fit->SetParError(0,1.27534e+09);
  fit->SetParError(1,2.06166);

  TH1F* h_counts = new TH1F("h_counts","",2,0,2);
  TH1F* h_HiggsWindowFit = new TH1F("h_HiggsWindowFit","# events in Higgs window from fit",80,0,40);
  TH1F* h_pull = new TH1F("h_pull","(fit - generated)/(fit error)",160,-8,8);
  TH1F* h_pull_fail = new TH1F("h_pull_fail","(fit - generated)/(fit error)",160,-8,8);
  TH1F* h_pull_lowErr  = new TH1F("h_pull_lowErr","(fit - generated)/(fit error)",160,-8,8);
  TH1F* h_pull_highErr = new TH1F("h_pull_highErr","(fit - generated)/(fit error)",160,-8,8);
  TH1F* h_diff = new TH1F("h_diff","(fit - generated)",160,-40,40);
  TH1F* h_pullU = new TH1F("h_pullU","(fit - generated)/(fit error)",160,-8,8);
  TH1F* h_diffU = new TH1F("h_diffU","(fit - generated)",160,-40,40);
  TH1F* h_pullL = new TH1F("h_pullL","(fit - generated)/(fit error)",160,-8,8);
  TH1F* h_diffL = new TH1F("h_diffL","(fit - generated)",160,-40,40);
  h_pullL->SetLineColor(kRed);h_pullL->SetMarkerColor(kRed);
  h_diffL->SetLineColor(kRed);h_diffL->SetMarkerColor(kRed);
  h_pullU->SetLineColor(kBlue);h_pullU->SetMarkerColor(kBlue);
  h_diffU->SetLineColor(kBlue);h_diffU->SetMarkerColor(kBlue);
  h_pull_fail->SetFillColor(kRed);
  //  TRandom3 rr;

  for(int j=0;j<10;j++){
 
    TH1F* h = new TH1F("h","",60,103,163);
    for(int i=0;i<5200;i++){
      //for(int i=0;i<rr.Poisson(52);i++){
      h->Fill(fit->GetRandom(103,163));
    }
    
    //h->Draw("PE");
    
    reject=true;
    TF1* fitCurve = new TF1("fitCurve",fpow,103,163,2);
    
    Double_t avg_l = h->Integral(h->FindBin(103),h->FindBin(118))/float(118-103),avg_u = h->Integral(h->FindBin(133),h->FindBin(163))/float(163-133),avgX_l=(118-103)/2.,avgX_u=(163-133)/2.;
    cout<<avg_l<<"  "<<avg_u<<"  "<<avgX_l<<"  "<<avgX_u<<endl;
    Double_t param1= (log(avg_l) - log(avg_u))/(log(avgX_l) - log(avgX_u));
    Double_t param0= /*7e14;*/avg_l/pow(avgX_l, param1);
    cout<<"param0: "<<param0<<"  param1: "<<param1<<endl;
    fitCurve->SetParameter(0,param0);
    fitCurve->SetParameter(1,param1);
    int status = h->Fit(fitCurve,"L","",103,163);
    //Then to get the result
    TFitResultPtr fitResult = h->Fit(fitCurve,"SLLMEV0","",103,163);
    TMatrixDSym cov = fitResult->GetCovarianceMatrix();
    fitResult->Print("V");
    h->GetXaxis()->SetRangeUser(100,164.9);
    reject=false;
    float YieldBinWidth=h->GetBinWidth(1);
    Double_t PowYieldSig = h->Integral(h->FindBin(120),h->FindBin(131-.1))/YieldBinWidth;
    Double_t PowYieldSigFit = fitCurve->Integral(120,131)/YieldBinWidth;
    Double_t PowYieldSigFitErr = fitCurve->IntegralError(120,131,fitResult->GetParams(),cov.GetMatrixArray() )/YieldBinWidth; 
    Double_t PowSBloYield = h->Integral(h->FindBin(103),h->FindBin(118-.1))/YieldBinWidth;
    Double_t PowSBloYieldFit = fitCurve->Integral(103,118)/YieldBinWidth;
    Double_t PowSBloYieldFitErr = fitCurve->IntegralError(103,118,fitResult->GetParams(),cov.GetMatrixArray() )/YieldBinWidth; 
    Double_t PowSBhiYield = h->Integral(h->FindBin(133),h->FindBin(163-.1))/YieldBinWidth;
    Double_t PowSBhiYieldFit = fitCurve->Integral(133,163)/YieldBinWidth;
    Double_t PowSBhiYieldFitErr = fitCurve->IntegralError(133,163,fitResult->GetParams(),cov.GetMatrixArray() )/YieldBinWidth; 
    
    cout<<"gMinuit->fStatus : "<< gMinuit->fStatus <<"  gMinuit->fCstatu : "<< gMinuit->fCstatu<<endl;

    cout<<" low sideband yield from histo: "<<PowSBloYield<<"  and from fit: "<<PowSBloYieldFit<<" +- "<<PowSBloYieldFitErr<<endl; 
    cout<<" higgs window yield from histo: "<<PowYieldSig <<"  and from fit: "<<PowYieldSigFit<<" +- "<<PowYieldSigFitErr<<endl; 
    cout<<" high sideband yield from histo: "<<PowSBhiYield<<"  and from fit: "<<PowSBhiYieldFit<<" +- "<<PowSBhiYieldFitErr<<endl;
    reject=true;
    
    TString str = (TString)gMinuit->fCstatu;
    cout<<"str: "<<str;
    h_counts->Fill(0);
    if(str.Contains("FAILURE")){
      double pull_fail = (PowYieldSigFit-9.11)/PowYieldSigFitErr;
      h_pull_fail->Fill(pull_fail);
      //  continue;
    }
    h_counts->Fill(1);
    double pull = (PowYieldSigFit-9.11)/PowYieldSigFitErr;
    double pull_lowErr = (PowYieldSigFit-(9.11-1.6))/PowYieldSigFitErr;
    double pull_highErr = (PowYieldSigFit-(9.11+1.6))/PowYieldSigFitErr;
    double diff = (PowYieldSigFit-9.11);
    double pullL = (PowSBloYieldFit-PowSBloYield)/PowSBloYieldFitErr;
    double diffL = (PowSBloYieldFit-PowSBloYield);
    double pullU = (PowSBhiYieldFit-PowSBhiYield)/PowSBhiYieldFitErr;
    double diffU = (PowSBhiYieldFit-PowSBhiYield);
    h_pull->Fill(pull);
    h_pull_lowErr->Fill(pull_lowErr);
    h_pull_highErr->Fill(pull_highErr);
    h_diff->Fill(diff);
    h_pullU->Fill(pullU);
    h_diffU->Fill(diffU);
    h_pullL->Fill(pullL);
    h_diffL->Fill(diffL);
    h_HiggsWindowFit->Fill(PowYieldSigFit);
  }
  /*
  TFile fout("InvarMassFitToys.root","RECREATE");
  fout.cd();
  h_diff->Write();
  h_diffL->Write();
  h_diffU->Write();
  h_pull->Write();
  h_pull_lowErr->Write();
  h_pull_highErr->Write();
  h_pullL->Write();
  h_pullU->Write();
  h_HiggsWindowFit->Write();
  h_counts->Write();
  fout.Close();
  */
  //h->Draw();
  h_diffU->Draw();h_diff->Draw("SAMES");h_diffL->Draw("SAMES");
  TLegend *legd = new TLegend(.6,.55,.8,.85,"","brNDC");
  legd->SetFillStyle(0);legd->SetBorderSize(0);
  legd->AddEntry(h_diff,"Higgs window","l");
  legd->AddEntry(h_diffL,"lower sideband","l");
  legd->AddEntry(h_diffU,"upper sideband","l");
  legd->Draw();
  //c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_toys__diff.png");
  //c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_toys__diff.pdf");


  //c2->cd();
  h_pullU->Draw();h_pull->Draw("SAMES");h_pullL->Draw("SAMES");
  legd->Draw();
  //c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_toys__pull.png");
  //c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_toys__pull.pdf");


  h_diff->Fit("gaus");h_diff->SetTitle("Higgs window (fit - generated) with Gaussian fit");
  h_diff->Draw();
  //c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_toys__diff_gausFit.png");
  //c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_toys__diff_gausFit.pdf");

  h_pull->Fit("gaus");h_pull->SetTitle("Higgs window pull with Gaussian fit");//h_pull_fail->Fit("gaus");
  h_pull->Draw();//h_pull_fail->Draw("SAMES");
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_toys__pull_gausFit_withFail.png");
  c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_toys__pull_gausFit_withFail.pdf");

  h_HiggsWindowFit->Fit("gaus");
  h_HiggsWindowFit->Draw();
  //c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_toys__HiggsWindow_gausFit.png");
  //c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_toys__HiggsWindow_gausFit.pdf");
  
  h_counts->Draw();
  //c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_toys__counts.png");
  //c1->Print("Plots/Higgs/Exclusive_WeDataInvMassFit_toys__counts.pdf");

}

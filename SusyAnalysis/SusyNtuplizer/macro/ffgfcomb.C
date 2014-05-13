#include "TH1.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLatex.h"

using namespace std;

void ffgfcomb(){


  //TFile f("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_15GevCHfakeCut_PixelCutOnFakes_NewDiEMPtBins_19499pb_NewJetMatching_DissertationOutput.root","READ");
  //TFile f_nojet("Plots/Closure/handmade/moreEvents/met_reweighted_nojet_comb.root","READ");
  TFile f_nojet("Plots/Closure/handmade/moreEvents/met_closure_comb.root","READ");
  //TFile f_1jet("Plots/Closure/handmade/moreEvents/met_reweighted_1jet_comb.root","READ");
  //TFile f_2jet("Plots/Closure/handmade/moreEvents/met_reweighted_2jet_comb.root","READ");
  f_nojet.cd();

  gStyle->SetOptStat(0);gStyle->SetOptFit(0011);

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->cd();
  c1->SetLogz(0);

  TH2F* h_chi2 = new TH2F("h_chi2",";ff fraction  #omega_{ff};#gammaf fraction  #omega_{#gamma f}",100,0.,1.,100,0.,1.);

  int Nbins=5,NbinsScale=5;

  TH1F* gg = (TH1F*)f_nojet.Get("ggMet");gg->Sumw2();
  TH1F* ff = (TH1F*)f_nojet.Get("ffMet");ff->Sumw2();
  TH1F* gf = (TH1F*)f_nojet.Get("gammafakeMet");gf->Sumw2();
  //TH1F* fg = (TH1F*)f_nojet.Get("fgMet_reweightJet_binned");fg->Sumw2();

  TH1F* ggVec = (TH1F*)gg->Clone(); TH1F *ffVec[28],*gfVec[28];
  float scale=0.;
  for(int k=0;k<28;k++){
    ffVec[k]=(TH1F*)ff->Clone();
    gfVec[k]=(TH1F*)gf->Clone();
    //scale = ffVec[k]->Integral()/gfVec[k]->Integral();gfVec[k]->Scale(scale);//don't need since it gets scaled to gg anyway
    int i=k;
    if(k<2)i=2;
    else i=k;
    scale = ggVec->Integral(0,i+1)/ffVec[k]->Integral(0,i+1);ffVec[k]->Scale(scale);
    scale = ggVec->Integral(0,i+1)/gfVec[k]->Integral(0,i+1);gfVec[k]->Scale(scale);
    //scale = ggVec->Integral(0,NbinsScale)/ffVec[k]->Integral(0,NbinsScale);ffVec[k]->Scale(scale);
    //scale = ggVec->Integral(0,NbinsScale)/gfVec[k]->Integral(0,NbinsScale);gfVec[k]->Scale(scale);
  }

  TH1F* QCD = (TH1F*)ff->Clone();
  //gf->Add(fg);
  scale = ff->Integral()/gf->Integral();
  gf->Scale(scale);
  scale = gg->Integral(0,NbinsScale)/ff->Integral(0,NbinsScale);ff->Scale(scale);
  scale = gg->Integral(0,NbinsScale)/gf->Integral(0,NbinsScale);gf->Scale(scale);

  gg->SetLineWidth(2);gg->GetXaxis()->SetRangeUser(0,24);gg->Draw("histo");
  ff->SetLineColor(kBlue);ff->Draw("histoSAMES");
  gf->SetLineColor(kRed);gf->Draw("histoSAMES");
  TLegend leg2(.62,.52,.92,.72,"","brNDC");leg2.SetFillStyle(0);leg2.SetShadowColor(0);leg2.SetBorderSize(0);
  leg2.AddEntry(gg,"#gamma#gamma","l");leg2.AddEntry(gf,"#gammaf","l");leg2.AddEntry(ff,"ff","l");//leg2.Draw();
  c1->Print("Plots/Closure/handmade/moreEvents/ffgfcomb_met_ggANDffANDgf.png");
  gg->GetXaxis()->SetRangeUser(0,1000);

  float chi2=0.;float min=999.,max=0.,minx=0,miny=0;
  for(int i=1;i<=h_chi2->GetNbinsX();i++){
    for(int j=1;j<=h_chi2->GetNbinsY();j++){
      chi2=0.;
      if(h_chi2->GetBinCenter(i)+h_chi2->GetBinCenter(j)==1){
	//cout<<h_chi2->GetBinCenter(i)+h_chi2->GetBinCenter(j)<<endl;
	QCD->Add(ff,gf,h_chi2->GetBinCenter(i),h_chi2->GetBinCenter(j));
	for(int k=1;k<=Nbins;k++){
	  float chi2Temp =(QCD->GetBinContent(k)-gg->GetBinContent(k))*(QCD->GetBinContent(k)-gg->GetBinContent(k))/(QCD->GetBinError(k)*QCD->GetBinError(k)+gg->GetBinError(k)*gg->GetBinError(k));
	  chi2+=chi2Temp;
	  //cout<<"i:"<<i<<" j:"<<j<<" k:"<<k<<"  chi2Temp:"<<chi2Temp<<"  chi2:"<<chi2<<" min:"<<min<<endl;
	}
	chi2/=Nbins;
	//if(i==20)cout<<"i:"<<i<<" j:"<<j<<" chi2: "<<chi2<<endl;
	int bin=h_chi2->FindBin(h_chi2->GetBinCenter(i)+.0001,h_chi2->GetBinCenter(j)+.0001);
	h_chi2->SetBinContent(bin,chi2);
	//cout<<"i:"<<i<<" j:"<<j<<" chi2:"<<chi2<<" min:"<<min<<endl;
	if(chi2<min && h_chi2->GetBinCenter(i)>0.01 && h_chi2->GetBinCenter(j)>0.01){min=chi2;minx=h_chi2->GetBinCenter(i);miny=h_chi2->GetBinCenter(j);}
	if(chi2>max)max=chi2;
	//cout<<"  chi2:"<<chi2<<"  min:"<<min<<endl<<endl;
	//cout<<"minx: "<<minx<<"  "<<h_chi2->GetBinCenter(i)<<endl;
	if(i==75){
	  gg->GetXaxis()->SetRangeUser(0,24);
	  gg->Draw("histo");ff->Draw("histoSAMES");gf->Draw("histoSAMES");QCD->SetLineColor(kGreen);QCD->Draw("histoSAMES");
	  c1->Print("Plots/Closure/handmade/moreEvents/ffgfcomb_met_ggANDffANDgfANDqcd.png");
	  gg->GetXaxis()->SetRangeUser(0,1000);
	}
      } 
    }
  }
  cout<<"no jet requirement MET<25 min chi2: "<<min<<"  max chi2: "<<max<<"  minx:"<<minx<<"  miny:"<<miny<<endl;
  char str[50];sprintf(str,"Minimum #chi^{2}/ndof: %4.3f",min);
  char str2[50];sprintf(str2,"ff fraction: %4.3f",minx);
  char str3[50];sprintf(str3,"#gammaf fraction: %4.3f",miny);
  TPaveText *text = new TPaveText(.42,.55,.83,.82,"NDC");
  //text->AddText("Minimum chi2:"+str);
  text->AddText("No jet requirement");text->AddText("E_{T}^{miss}<25");text->AddText(str);text->AddText(str2);text->AddText(str3);
  text->SetFillStyle(4000);
  text->SetFillColor(0);
  text->SetBorderSize(0);

  h_chi2->GetZaxis()->SetRangeUser(min-.01,max+.01);
  h_chi2->Draw("colz");text->Draw();
  c1->Print("Plots/Closure/handmade/moreEvents/ffgfcomb_Chi2_025.png");
  c1->Print("Plots/Closure/handmade/moreEvents/ffgfcomb_Chi2_025.pdf");
  
  //f_nojet.Close();
  f_nojet.cd();
  c1->cd();
  //now jet requirement:
  
  TH2F* h_chi2_jet = new TH2F("h_chi2_jet",";ff fraction  #omega_{ff};#gammaf fraction  #omega_{#gamma f}",100,0.,1.,100,0.,1.);
  //TH1F* gg_jet = (TH1F*)f_nojet.Get("ggMet_JetReq");gg_jet->Sumw2();
  //TH1F* ff_jet = (TH1F*)f_nojet.Get("ffMet_reweightJet_binned_JetReq");ff_jet->Sumw2();
  //TH1F* gf_jet = (TH1F*)f_nojet.Get("gammafakeMet_reweightJet_binned_JetReq");gf_jet->Sumw2();
  //TH1F* fg_jet = (TH1F*)f_nojet.Get("fgMet_reweightJet_binned_JetReq");fg_jet->Sumw2();
  
  TH1F* gg_jet = (TH1F*)f_nojet.Get("ggMet_JetReq");gg_jet->Sumw2();
  TH1F* ff_jet = (TH1F*)f_nojet.Get("ffMet_JetReq");ff_jet->Sumw2();
  TH1F* gf_jet = (TH1F*)f_nojet.Get("gammafakeMet_JetReq");gf_jet->Sumw2();

  TH1F* ggVec_jet = (TH1F*)gg_jet->Clone(); TH1F *ffVec_jet[28],*gfVec_jet[28];
  for(int k=0;k<28;k++){
    ffVec_jet[k]=(TH1F*)ff_jet->Clone();
    gfVec_jet[k]=(TH1F*)gf_jet->Clone();
    //scale = gfVec_jet[k]->Integral()/ffVec_jet[k]->Integral();gfVec_jet[k]->Scale(scale);
    int i=k;
    if(k<4)i=3;
    else i=k;
    scale = ggVec_jet->Integral(0,i+1)/ffVec_jet[k]->Integral(0,i+1);ffVec_jet[k]->Scale(scale);
    scale = ggVec_jet->Integral(0,i+1)/gfVec_jet[k]->Integral(0,i+1);gfVec_jet[k]->Scale(scale);
    //scale = ggVec_jet->Integral(0,NbinsScale)/ffVec_jet[k]->Integral(0,NbinsScale);ffVec_jet[k]->Scale(scale);
    //scale = ggVec_jet->Integral(0,NbinsScale)/gfVec_jet[k]->Integral(0,NbinsScale);gfVec_jet[k]->Scale(scale);
  }

 TH1F* QCD_jet = (TH1F*)ff_jet->Clone();
  //gf_jet->Add(fg_jet);
  scale = ff_jet->Integral()/gf_jet->Integral();
  gf_jet->Scale(scale);
  scale = gg_jet->Integral(0,NbinsScale)/ff_jet->Integral(0,NbinsScale);ff_jet->Scale(scale);
  scale = gg_jet->Integral(0,NbinsScale)/gf_jet->Integral(0,NbinsScale);gf_jet->Scale(scale);

  gg_jet->SetLineWidth(2);gg_jet->GetXaxis()->SetRangeUser(0,49);gg_jet->Draw("histo");
  ff_jet->SetLineColor(kBlue);ff_jet->Draw("histoSAMES");
  gf_jet->SetLineColor(kRed);gf_jet->Draw("histoSAMES");
  leg2.Draw();
  c1->Print("Plots/Closure/handmade/moreEvents/ffgfcomb_met_ggANDffANDgf_jet.png");
  gg_jet->GetXaxis()->SetRangeUser(0,1000);

  chi2=0.;min=999.;max=0.;minx=0;miny=0;
  for(int i=1;i<=h_chi2_jet->GetNbinsX();i++){
    for(int j=1;j<=h_chi2_jet->GetNbinsY();j++){
      chi2=0.;
      if(h_chi2_jet->GetBinCenter(i)+h_chi2_jet->GetBinCenter(j)==1){
	//cout<<h_chi2_jet->GetBinCenter(i)+h_chi2_jet->GetBinCenter(j)<<endl;
	QCD_jet->Add(ff_jet,gf_jet,h_chi2_jet->GetBinCenter(i),h_chi2_jet->GetBinCenter(j));
	for(int k=1;k<=Nbins;k++){
	  chi2+=(QCD_jet->GetBinContent(k)-gg_jet->GetBinContent(k))*(QCD_jet->GetBinContent(k)-gg_jet->GetBinContent(k))/(QCD_jet->GetBinError(k)*QCD_jet->GetBinError(k)+gg_jet->GetBinError(k)*gg_jet->GetBinError(k));
	  //chi2Err2+=(2*(QCD_jet->GetBinContent(k)-gg_jet->GetBinContent(k))*sqrt(1/(QCD_jet->GetBinError(k)*QCD_jet->GetBinError(k)+gg_jet->GetBinError(k)*gg_jet->GetBinError(k))))*(2*(QCD_jet->GetBinContent(k)-gg_jet->GetBinContent(k))*sqrt(1/(QCD_jet->GetBinError(k)*QCD_jet->GetBinError(k)+gg_jet->GetBinError(k)*gg_jet->GetBinError(k))));
	}
	chi2/=Nbins;//cout<<"i:"<<i<<" j:"<<j<<" chi2: "<<chi2<<endl;
	//chi2Err=sqrt(chi2Err2);chi2Err/=Nbins;
	int bin=h_chi2_jet->FindBin(h_chi2_jet->GetBinCenter(i)+.0001,h_chi2_jet->GetBinCenter(j)+.0001);
	h_chi2_jet->SetBinContent(bin,chi2);
	if(chi2<min && h_chi2_jet->GetBinCenter(i)>0.01 && h_chi2_jet->GetBinCenter(j)>0.01){min=chi2;minx=h_chi2_jet->GetBinCenter(i);miny=h_chi2_jet->GetBinCenter(j);}
	if(chi2>max)max=chi2;
	//cout<<"minx: "<<minx<<h_chi2_jet->GetXaxis()_jet->GetBinCenter(i)<<endl;
      } 
    }
  }
  cout<<">=1 jet requirement MET<25 min chi2: "<<min<<"  max chi2: "<<max<<endl;
  sprintf(str,"Minimum #chi^{2}/ndof: %4.3f",min);
  sprintf(str2,"ff fraction: %4.3f",minx);
  sprintf(str3,"#gammaf fraction: %4.3f",miny);
  TPaveText *text_jet = new TPaveText(.42,.55,.83,.82,"NDC");
  //text_jet->AddText("Minimum chi2:"+str);
  text_jet->AddText("#geq 1 jet requirement");text_jet->AddText("E_{T}^{miss}<25");text_jet->AddText(str);text_jet->AddText(str2);text_jet->AddText(str3);
  text_jet->SetFillStyle(4000);
  text_jet->SetFillColor(0);
  text_jet->SetBorderSize(0);

  h_chi2_jet->GetZaxis()->SetRangeUser(min-.01,max+.01);
  h_chi2_jet->Draw("colz");text_jet->Draw();
  c1->Print("Plots/Closure/handmade/moreEvents/ffgfcomb_Chi2_025_jet.png");
  c1->Print("Plots/Closure/handmade/moreEvents/ffgfcomb_Chi2_025_jet.pdf");
  
  Nbins=28;
  //float fracYmin=.31,fracYmax=.8;
  float fracYmin=0,fracYmax=1;
  TH1F* h_ffFracVsNbins = new TH1F("h_ffFracVsNbins",";N_{bin};fraction #omega",Nbins,0.5,Nbins+.5);h_ffFracVsNbins->SetLineColor(kBlue);h_ffFracVsNbins->GetYaxis()->SetRangeUser(fracYmin,fracYmax); 
  TH1F* h_gfFracVsNbins = new TH1F("h_gfFracVsNbins",";N_{bin};fraction #omega",Nbins,0.5,Nbins+.5);h_gfFracVsNbins->SetLineColor(kRed);
  h_ffFracVsNbins->SetLineWidth(2);h_gfFracVsNbins->SetLineWidth(2);
  chi2=0.;min=999.;max=0.;minx=0;miny=0;
  float chi2Vec[28]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},minVec[28]={999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.},minxVec[28]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},minyVec[28]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  for(int i=1;i<=h_chi2->GetNbinsX();i++){//cout<<"i="<<i<<endl;
    for(int j=1;j<=h_chi2->GetNbinsY();j++){//cout<<"j="<<j<<endl;
      chi2=0.;
      for(int p=0;p<Nbins;p++){//cout<<"p="<<p<<endl;
	chi2Vec[p]=0.;
      }
      if(h_chi2->GetBinCenter(i)+h_chi2->GetBinCenter(j)==1){
	//cout<<h_chi2->GetBinCenter(i)+h_chi2->GetBinCenter(j)<<endl;
	for(int m=1;m<=Nbins;m++){
	  //if(m==1)cout<<" gg bin1:"<<ggVec->GetBinContent(1)<<" ff bin1:"<<ffVec[m-1]->GetBinContent(1)<<" gf bin1:"<<gfVec[m-1]->GetBinContent(1)<<endl;
	  QCD->Add(ffVec[m-1],gfVec[m-1],h_chi2->GetBinCenter(i),h_chi2->GetBinCenter(j));
	  chi2=0.;float chi2Temp=0;
	  for(int k=1;k<=Nbins;k++){
	    float ggVal=ggVec->GetBinContent(1),QCDval=(ggVec->GetBinContent(1)/ggVec->GetBinContent(k))*QCD->GetBinContent(k),ggValErr=(ggVec->GetBinContent(1)/ggVec->GetBinContent(k))*ggVec->GetBinError(k),QCDvalErr=(ggVec->GetBinContent(1)/ggVec->GetBinContent(k))*QCD->GetBinError(k);
	    chi2Temp=(QCD->GetBinContent(k)-ggVec->GetBinContent(k))*(QCD->GetBinContent(k)-ggVec->GetBinContent(k))/(QCD->GetBinError(k)*QCD->GetBinError(k)+ggVec->GetBinError(k)*ggVec->GetBinError(k));
	    //chi2Temp=(QCDval-ggVal)*(QCDval-ggVal)/(QCDvalErr*QCDvalErr+ggValErr*ggValErr);
	    //cout<<"ggVal: "<<ggVal<<"  QCDval: "<<QCDval<<"  ggCont: "<< ggVec->GetBinContent(k) << "  QCDcont: "<<QCD->GetBinContent(k)<<endl;
	    //cout<<"ggValErr: "<<ggValErr<<"  QCDvalErr: "<<QCDvalErr<<"  ggErr: "<< ggVec->GetBinError(k) << "  QCDerr: "<<QCD->GetBinError(k)<<endl;
	    chi2+=chi2Temp;
	    //if(k==m && i==1)cout<<"k:"<<k<<" i:"<<i<<" j:"<<j<<" chi2Temp: "<<chi2Temp<<"  chi2: "<<chi2<<endl;
	    if(k==m)chi2Vec[k-1]=chi2;
	  }
	}
	chi2/=Nbins;//cout<<"i:"<<i<<" j:"<<j<<" chi2: "<<chi2<<endl;
	for(int p=0;p<Nbins;p++){//cout<<"p="<<p<<endl;
	  chi2Vec[p]/=p+1;//element 0 = bin 1
	  //if(i==20){cout<<"i:"<<i<<" j:"<<j<<" p:"<<p<<" chi2Vec[p]:"<<chi2Vec[p]<<" minVec[p]:"<<minVec[p]<<endl;}
	  if(chi2Vec[p]<minVec[p] && h_chi2->GetBinCenter(i)>0.01 && h_chi2->GetBinCenter(j)>0.01){minVec[p]=chi2Vec[p];minxVec[p]=h_chi2->GetBinCenter(i);minyVec[p]=h_chi2->GetBinCenter(j);}
	}
	int bin=h_chi2->FindBin(h_chi2->GetBinCenter(i)+.0001,h_chi2->GetBinCenter(j)+.0001);
	h_chi2->SetBinContent(bin,chi2);
	//if(i==20)cout<<"i:"<<i<<" j:"<<j<<" chi2: "<<chi2<<endl;
	if(chi2<min && h_chi2->GetBinCenter(i)>0.01 && h_chi2->GetBinCenter(j)>0.01){min=chi2;minx=h_chi2->GetBinCenter(i);miny=h_chi2->GetBinCenter(j);}
	if(chi2>max)max=chi2;
	//cout<<"minx: "<<minx<<h_chi2->GetXaxis()->GetBinCenter(i)<<endl;
      } 
    }
  }
  for(int q=0;q<Nbins;q++){
    //cout<<"q="<<q<<"  minVec[q]="<<minVec[q]<<endl;
    h_ffFracVsNbins->SetBinContent(q+1,minxVec[q]);
    h_gfFracVsNbins->SetBinContent(q+1,minyVec[q]);
    h_ffFracVsNbins->SetBinError(q+1,0.1);
    h_gfFracVsNbins->SetBinError(q+1,0.1);
    //cout<<"q:"<<q<<"  ff:"<<h_ffFracVsNbins->GetBinContent(q+1)<<"  gf:"<<h_gfFracVsNbins->GetBinContent(q+1)<<endl;
  }
  TPaveText *text = new TPaveText(.65,.748,.83,.82,"NDC");text->AddText("no jet req");text->SetFillStyle(4000);text->SetFillColor(0);text->SetBorderSize(0);
  TLegend* leg = new TLegend(.24,.635,.54,.835,"","brNDC");leg->SetFillStyle(0);leg->SetShadowColor(0);leg->SetBorderSize(0);
  leg->AddEntry(h_ffFracVsNbins,"ff","l"); 
  leg->AddEntry(h_gfFracVsNbins,"#gammaf","l"); 
  h_ffFracVsNbins->Fit("pol0","","",.5,Nbins+.5);TF1* fffit = (TF1*)h_ffFracVsNbins->GetFunction("pol0");fffit->SetLineColor(kBlue);
  h_gfFracVsNbins->Fit("pol0","","",.5,Nbins+.5);TF1* gffit = (TF1*)h_gfFracVsNbins->GetFunction("pol0");
  h_ffFracVsNbins->Draw("histo");
  h_gfFracVsNbins->Draw("histoSAMES");
  fffit->Draw("SAMES");gffit->Draw("SAMES");
  c1->Update();
  TPaveStats* sbff=(TPaveStats*)(h_ffFracVsNbins->GetListOfFunctions()->FindObject("stats"));
  sbff->SetX1NDC(.39);sbff->SetX2NDC(.6);sbff->SetY1NDC(.71);sbff->SetY2NDC(.86);sbff->SetBorderSize(0);sbff->SetTextColor(kBlue);sbff->SetFillStyle(0);
  TPaveStats* sbgf=(TPaveStats*)(h_gfFracVsNbins->GetListOfFunctions()->FindObject("stats"));
  sbgf->SetX1NDC(.39);sbgf->SetX2NDC(.6);sbgf->SetY1NDC(.65);sbgf->SetY2NDC(.73);sbgf->SetBorderSize(0);sbgf->SetTextColor(kRed);sbgf->SetFillStyle(0);
  TPaveText *textEq = new TPaveText(.45,.635,.48,.82,"NDC");textEq->SetFillStyle(4000);textEq->SetFillColor(0);textEq->SetBorderSize(0);
  textEq->AddText("=");textEq->AddText("=");
  leg->Draw();text->Draw();textEq->Draw();
  c1->Print("Plots/Closure/handmade/moreEvents/ffgfcomb_ffFracVsNbins.png");
  c1->Print("Plots/Closure/handmade/moreEvents/ffgfcomb_ffFracVsNbins.pdf");
  
  cout<<"no jet requirement MET<50 min chi2: "<<min<<"  max chi2: "<<max<<"  minx: "<<minx<<"  miny: "<<miny<<endl;
  sprintf(str,"Minimum #chi^{2}/ndof: %4.3f",min);
  sprintf(str2,"ff fraction: %4.3f",minx);
  sprintf(str3,"#gammaf fraction: %4.3f",miny);
  TPaveText *text2 = new TPaveText(.42,.55,.83,.82,"NDC");
  //text->AddText("Minimum chi2:"+str);
  text2->AddText("No jet requirement");text2->AddText("E_{T}^{miss}<50");text2->AddText(str);text2->AddText(str2);text2->AddText(str3);
  text2->SetFillStyle(4000);
  text2->SetFillColor(0);
  text2->SetBorderSize(0);
  h_chi2->GetZaxis()->SetRangeUser(min-.01,max+.01);
  h_chi2->Draw("colz");text2->Draw();
  c1->Print("Plots/Closure/handmade/moreEvents/ffgfcomb_Chi2_050.png");
  c1->Print("Plots/Closure/handmade/moreEvents/ffgfcomb_Chi2_050.pdf");

  TH1F* h_ffFracVsNbins_jet = new TH1F("h_ffFracVsNbins_jet",";N_{bin};fraction #omega",Nbins,0.5,Nbins+.5);h_ffFracVsNbins_jet->SetLineColor(kBlue);  h_ffFracVsNbins_jet->GetYaxis()->SetRangeUser(fracYmin,fracYmax);
  TH1F* h_gfFracVsNbins_jet = new TH1F("h_gfFracVsNbins_jet",";N_{bin};fraction #omega",Nbins,0.5,Nbins+.5);h_gfFracVsNbins_jet->SetLineColor(kRed);
  h_ffFracVsNbins_jet->SetLineWidth(2);h_gfFracVsNbins_jet->SetLineWidth(2);
  chi2=0.;min=999.;max=0.;minx=0;miny=0;  
  float chi2Vec_jet[28]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},minVec_jet[28]={999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.},minxVec_jet[28]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},minyVec_jet[28]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  for(int i=1;i<=h_chi2_jet->GetNbinsX();i++){
    for(int j=1;j<=h_chi2_jet->GetNbinsY();j++){
      chi2=0.;
      for(int p=0;p<Nbins;p++){//cout<<"p="<<p<<endl;
	chi2Vec_jet[p]=0.;
      }
      if(h_chi2_jet->GetBinCenter(i)+h_chi2_jet->GetBinCenter(j)==1){
	//cout<<h_chi2_jet->GetBinCenter(i)+h_chi2_jet->GetBinCenter(j)<<endl;
	for(int m=1;m<=Nbins;m++){
	  //if(m==1)cout<<" gg_jet bin1:"<<ggVec_jet->GetBinContent(1)<<" ff_jet bin1:"<<ffVec_jet[m-1]->GetBinContent(1)<<" gf_jet bin1:"<<gfVec_jet[m-1]->GetBinContent(1)<<endl;
	  QCD_jet->Add(ffVec_jet[m-1],gfVec_jet[m-1],h_chi2_jet->GetBinCenter(i),h_chi2_jet->GetBinCenter(j));
	  chi2=0.;
	  for(int k=1;k<=Nbins;k++){
	    chi2+=(QCD_jet->GetBinContent(k)-ggVec_jet->GetBinContent(k))*(QCD_jet->GetBinContent(k)-ggVec_jet->GetBinContent(k))/(QCD_jet->GetBinError(k)*QCD_jet->GetBinError(k)+ggVec_jet->GetBinError(k)*ggVec_jet->GetBinError(k));
	    if(k==m)chi2Vec_jet[k-1]=chi2;
	  }
	}
	chi2/=Nbins;//cout<<"i:"<<i<<" j:"<<j<<" chi2: "<<chi2<<endl;
	for(int p=0;p<Nbins;p++){
	  //cout<<"p="<<p<<endl;
	  chi2Vec_jet[p]/=p+1;//element 0 = bin 1
	  if(chi2Vec_jet[p]<minVec_jet[p] && h_chi2_jet->GetBinCenter(i)>0.01 && h_chi2_jet->GetBinCenter(j)>0.01){minVec_jet[p]=chi2Vec_jet[p];minxVec_jet[p]=h_chi2_jet->GetBinCenter(i);minyVec_jet[p]=h_chi2_jet->GetBinCenter(j);}
	}
	int bin=h_chi2_jet->FindBin(h_chi2_jet->GetBinCenter(i)+.0001,h_chi2_jet->GetBinCenter(j)+.0001);
	h_chi2_jet->SetBinContent(bin,chi2);
	if(chi2<min && h_chi2_jet->GetBinCenter(i)>0.01 && h_chi2_jet->GetBinCenter(j)>0.01){min=chi2;minx=h_chi2_jet->GetXaxis()->GetBinCenter(i);miny=h_chi2_jet->GetYaxis()->GetBinCenter(j);}
	if(chi2>max)max=chi2;
	//cout<<"minx: "<<minx<<"  "<<h_chi2_jet->GetXaxis()_jet->GetBinCenter(i)<<endl;
      } 
    }
  }
  for(int q=0;q<Nbins;q++){
    //cout<<"q="<<q<<"  minVec_jet[q]="<<minVec_jet[q]<<endl;
    h_ffFracVsNbins_jet->SetBinContent(q+1,minxVec_jet[q]);
    h_gfFracVsNbins_jet->SetBinContent(q+1,minyVec_jet[q]);
    h_ffFracVsNbins_jet->SetBinError(q+1,0.1);
    h_gfFracVsNbins_jet->SetBinError(q+1,0.1);
    //cout<<"q:"<<q<<"  ff_jet:"<<h_ffFracVsNbins_jet->GetBinContent(q+1)<<"  gf_jet:"<<h_gfFracVsNbins_jet->GetBinContent(q+1)<<endl;
  }
  TPaveText *text_jet = new TPaveText(.65,.748,.83,.82,"NDC");text_jet->AddText("#geq1 jet req");text_jet->SetFillStyle(4000);text_jet->SetFillColor(0);text_jet->SetBorderSize(0);
  h_ffFracVsNbins_jet->Fit("pol0","","",.5,Nbins+.5);TF1* fffit_jet = (TF1*)h_ffFracVsNbins_jet->GetFunction("pol0");fffit_jet->SetLineColor(kBlue);
  h_gfFracVsNbins_jet->Fit("pol0","","",.5,Nbins+.5);TF1* gffit_jet = (TF1*)h_gfFracVsNbins_jet->GetFunction("pol0");
  h_ffFracVsNbins_jet->Draw("histo");h_gfFracVsNbins_jet->Draw("histoSAMES");
  fffit_jet->Draw("SAMES");gffit_jet->Draw("SAMES");
  c1->Update();
  TPaveStats* sbff=(TPaveStats*)(h_ffFracVsNbins_jet->GetListOfFunctions()->FindObject("stats"));
  sbff->SetX1NDC(.39);sbff->SetX2NDC(.6);sbff->SetY1NDC(.71);sbff->SetY2NDC(.86);sbff->SetBorderSize(0);sbff->SetTextColor(kBlue);sbff->SetFillStyle(0);
  TPaveStats* sbgf=(TPaveStats*)(h_gfFracVsNbins_jet->GetListOfFunctions()->FindObject("stats"));
  sbgf->SetX1NDC(.39);sbgf->SetX2NDC(.6);sbgf->SetY1NDC(.65);sbgf->SetY2NDC(.73);sbgf->SetBorderSize(0);sbgf->SetTextColor(kRed);sbgf->SetFillStyle(0);
  leg->Draw();text_jet->Draw();textEq->Draw();
  c1->Print("Plots/Closure/handmade/moreEvents/ffgfcomb_ffFracVsNbins_jet.png");
  c1->Print("Plots/Closure/handmade/moreEvents/ffgfcomb_ffFracVsNbins_jet.pdf");
  
  cout<<">=1 jet requirement MET<50 min chi2: "<<min<<"  max chi2: "<<max<<endl;
  sprintf(str,"Minimum #chi^{2}/ndof: %4.3f",min);
  sprintf(str2,"ff fraction: %4.3f",minx);
  sprintf(str3,"#gammaf fraction: %4.3f",miny);
  TPaveText *text2_jet = new TPaveText(.42,.55,.83,.82,"NDC");
  text2_jet->AddText("#geq 1 jet requirement");text2_jet->AddText("E_{T}^{miss}<50");text2_jet->AddText(str);text2_jet->AddText(str2);text2_jet->AddText(str3);
  text2_jet->SetFillStyle(4000);
  text2_jet->SetFillColor(0);
  text2_jet->SetBorderSize(0);
  h_chi2_jet->GetZaxis()->SetRangeUser(min-.01,max+.01);
  h_chi2_jet->Draw("colz");text2_jet->Draw();
  c1->Print("Plots/Closure/handmade/moreEvents/ffgfcomb_Chi2_050_jet.png");
  c1->Print("Plots/Closure/handmade/moreEvents/ffgfcomb_Chi2_050_jet.pdf");
  

  f_nojet.Close();
}

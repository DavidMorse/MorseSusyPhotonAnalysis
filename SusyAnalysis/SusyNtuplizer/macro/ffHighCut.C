#include "TH1.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveText.h"

void ffHighCut(){

  //TFile f("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_50GevCHfakeCut_PixelCutOnFakes_NewDiEMPtBins_19499pb.root");  
  TFile f("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_50GevCHfakeCut_PixelCutOnFakes_NewDiEMPtBins_19499pb_NewJetMatching_2JetReq_Apr23.root","READ");

  f.cd();

  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->cd();
  c1->SetLogy(0);

  TPaveText *Text;
  Text = new TPaveText(.26,.68,.56,.82,"NDC");
  Text->AddText("CMS Preliminary 2012");
  Text->AddText("#sqrt{s} = 8 TeV, #intL = 19.5 fb^{-1}");
  Text->SetFillStyle(4000);
  Text->SetFillColor(0);
  Text->SetBorderSize(0);

  TPaveText *TextRight = new TPaveText(.46,.68,.76,.82,"NDC");
  TextRight->AddText("CMS Preliminary 2012");
  TextRight->AddText("#sqrt{s} = 8 TeV, #intL = 19.5 fb^{-1}");
  TextRight->SetFillStyle(4000);
  TextRight->SetFillColor(0);
  TextRight->SetBorderSize(0);


  float Chi25[50]={0.},Chi50[50]={0.},Chi25_JetReq[50]={0.},Chi50_JetReq[50]={0.};
  float PfChi25[50]={0.},PfChi50[50]={0.},PfChi25_JetReq[50]={0.},PfChi50_JetReq[50]={0.};
  float PfChi25gf[50]={0.},PfChi50gf[50]={0.},PfChi25gf_JetReq[50]={0.},PfChi50gf_JetReq[50]={0.};
  float PfChi25gammafake[50]={0.},PfChi50gammafake[50]={0.},PfChi25gammafake_JetReq[50]={0.},PfChi50gammafake_JetReq[50]={0.};
  float PfChi25fg[50]={0.},PfChi50fg[50]={0.},PfChi25fg_JetReq[50]={0.},PfChi50fg_JetReq[50]={0.};
  TString fftitle="",gftitle="",fgtitle="",gammafaketitle="";
  TH1F *ff;
  TH1F *ff_JetReq;
  TH1F *gf;
  TH1F *gf_JetReq;
  TH1F *gammafake;
  TH1F *gammafake_JetReq;
  TH1F *fg;
  TH1F *fg_JetReq;
  
  for(int i=1;i<50;i++){
    fftitle="Toys/ffMet_Cut_";fftitle+=i;fftitle+="GeV";
    ff = (TH1F*)f.Get(fftitle);  
    ff->GetYaxis()->SetRangeUser(0,22000);
    ff->GetXaxis()->SetRangeUser(0,120);
    ff->SetLineColor(4*i);
    if(i==20){ff->SetLineWidth(5);cout<<"ffMet_Cut_20GeV entries: "<<ff->GetEntries()<<endl;}
    if(i==1)ff->Draw();
    else ff->Draw("SAME");
  }
  TextRight->Draw();
  //c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ffChi2.png");

  for(int i=1;i<50;i++){
    fftitle="Toys/ffMet_Cut_";fftitle+=i;fftitle+="GeV_JetReq";
    ff = (TH1F*)f.Get(fftitle);  
    ff->GetYaxis()->SetRangeUser(0,7200);
    ff->GetXaxis()->SetRangeUser(0,120);
    ff->SetLineColor(4*i);
    if(i==20){ff->SetLineWidth(5);cout<<"ffMet_Cut_20GeV_JetReq entries: "<<ff->GetEntries()<<endl;}
    if(i==1)ff->Draw();
    else ff->Draw("SAME");
  }
  TextRight->Draw();
  //c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ffChi2_JetReq.png");

for(int i=1;i<50;i++){
    fftitle="Toys/ffMet_Cut_";fftitle+=i;fftitle+="GeV_PF";
    ff = (TH1F*)f.Get(fftitle);  
    ff->GetYaxis()->SetRangeUser(0,21000);
    ff->GetXaxis()->SetRangeUser(0,120);
    ff->SetLineColor(4*i);
    if(i==27){cout<<"i:"<<i<<" Entries:"<<ff->GetEntries()<<endl;}
    if(i==20){ff->SetLineWidth(5);cout<<"ffMet_Cut_20GeV_PF entries: "<<ff->GetEntries()<<endl;}
    if(i==1)ff->Draw();
    else ff->Draw("SAME");
  }
  TextRight->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ffChi2_PF.png");

  for(int i=1;i<50;i++){
    fftitle="Toys/ffMet_Cut_";fftitle+=i;fftitle+="GeV_PF_JetReq";
    ff = (TH1F*)f.Get(fftitle);  
    ff->GetYaxis()->SetRangeUser(0,7000);
    ff->GetXaxis()->SetRangeUser(0,120);
    ff->SetLineColor(4*i);
    if(i==27){cout<<"i:"<<i<<" Entries_JetReq:"<<ff->GetEntries()<<endl;}
    if(i==20){ff->SetLineWidth(5);cout<<"ffMet_Cut_20GeV_PF_JetReq entries: "<<ff->GetEntries()<<endl;}
    if(i==1)ff->Draw();
    else ff->Draw("SAME");
  }
  TextRight->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/ffChi2_PF_JetReq.png");

  c1->SetLogy(0);

  float scale=0., scale_JetReq=0.;

  for(int i=5;i<50;i++){
    for(int j=1;j<=5;j++){
      fftitle="Toys/ffMet_Cut_";fftitle+=i;fftitle+="GeV";
      ff = (TH1F*)f.Get(fftitle);  
      fftitle="Toys/ffMet_Cut_";fftitle+=i;fftitle+="GeV_JetReq";
      ff_JetReq = (TH1F*)f.Get(fftitle);
      if(ff->Integral(0,5)>0)scale=ggMet->Integral(0,5)/ff->Integral(0,5);
      else{scale=0.;}
      if(ff_JetReq->Integral(0,5)>0)scale_JetReq=ggMet_JetReq->Integral(0,5)/ff_JetReq->Integral(0,5);
      else scale_JetReq=0.;
      ff->Scale(scale);ff_JetReq->Scale(scale_JetReq);
      Chi25[i]+=(((ff->GetBinContent(j)-ggMet->GetBinContent(j))*(ff->GetBinContent(j)-ggMet->GetBinContent(j)))/(ff->GetBinError(j)*ff->GetBinError(j)+ggMet->GetBinError(j)*ggMet->GetBinError(j)));
      Chi25_JetReq[i]+=(((ff_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j))*(ff_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j)))/(ff_JetReq->GetBinError(j)*ff_JetReq->GetBinError(j)+ggMet_JetReq->GetBinError(j)*ggMet_JetReq->GetBinError(j)));
      
    }
    
  }
 for(int i=5;i<50;i++){
    for(int j=1;j<=10;j++){
      fftitle="Toys/ffMet_Cut_";fftitle+=i;fftitle+="GeV";
      ff = (TH1F*)f.Get(fftitle);  
      fftitle="Toys/ffMet_Cut_";fftitle+=i;fftitle+="GeV_JetReq";
      ff_JetReq = (TH1F*)f.Get(fftitle);
      if(ff->Integral(0,10)>0)scale=ggMet->Integral(0,10)/ff->Integral(0,10);
      else scale=0.;
      if(ff_JetReq->Integral(0,10)>0)scale_JetReq=ggMet_JetReq->Integral(0,10)/ff_JetReq->Integral(0,10);
      else scale_JetReq=0.;
      ff->Scale(scale);ff_JetReq->Scale(scale_JetReq);
      Chi50[i]+=(((ff->GetBinContent(j)-ggMet->GetBinContent(j))*(ff->GetBinContent(j)-ggMet->GetBinContent(j)))/(ff->GetBinError(j)*ff->GetBinError(j)+ggMet->GetBinError(j)*ggMet->GetBinError(j)));
      Chi50_JetReq[i]+=(((ff_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j))*(ff_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j)))/(ff_JetReq->GetBinError(j)*ff_JetReq->GetBinError(j)+ggMet_JetReq->GetBinError(j)*ggMet_JetReq->GetBinError(j)));
    }
 
 }

 //Now PF

  for(int i=5;i<50;i++){
    for(int j=1;j<=5;j++){
      fftitle="Toys/ffMet_Cut_";fftitle+=i;fftitle+="GeV_PF";
      ff = (TH1F*)f.Get(fftitle);  
      fftitle="Toys/ffMet_Cut_";fftitle+=i;fftitle+="GeV_PF_JetReq";
      ff_JetReq = (TH1F*)f.Get(fftitle);
      if(ff->Integral(0,5)>0)scale=ggMet->Integral(0,5)/ff->Integral(0,5);
      else scale=0.;
      if(ff_JetReq->Integral(0,5)>0)scale_JetReq=ggMet_JetReq->Integral(0,5)/ff_JetReq->Integral(0,5);
      else scale_JetReq=0.;
      ff->Scale(scale);ff_JetReq->Scale(scale_JetReq);
      PfChi25[i]+=(((ff->GetBinContent(j)-ggMet->GetBinContent(j))*(ff->GetBinContent(j)-ggMet->GetBinContent(j)))/(ff->GetBinError(j)*ff->GetBinError(j)+ggMet->GetBinError(j)*ggMet->GetBinError(j)));
      PfChi25_JetReq[i]+=(((ff_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j))*(ff_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j)))/(ff_JetReq->GetBinError(j)*ff_JetReq->GetBinError(j)+ggMet_JetReq->GetBinError(j)*ggMet_JetReq->GetBinError(j)));
      
    }
    
  }
 for(int i=5;i<50;i++){
    for(int j=1;j<=10;j++){
      fftitle="Toys/ffMet_Cut_";fftitle+=i;fftitle+="GeV_PF";
      ff = (TH1F*)f.Get(fftitle);  
      fftitle="Toys/ffMet_Cut_";fftitle+=i;fftitle+="GeV_PF_JetReq";
      ff_JetReq = (TH1F*)f.Get(fftitle);
      if(ff->Integral(0,10)>0)scale=ggMet->Integral(0,10)/ff->Integral(0,10);
      else scale=0.;
      if(ff_JetReq->Integral(0,10)>0)scale_JetReq=ggMet_JetReq->Integral(0,10)/ff_JetReq->Integral(0,10);
      else scale_JetReq=0.;
      ff->Scale(scale);ff_JetReq->Scale(scale_JetReq);
      PfChi50[i]+=(((ff->GetBinContent(j)-ggMet->GetBinContent(j))*(ff->GetBinContent(j)-ggMet->GetBinContent(j)))/(ff->GetBinError(j)*ff->GetBinError(j)+ggMet->GetBinError(j)*ggMet->GetBinError(j)));
      PfChi50_JetReq[i]+=(((ff_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j))*(ff_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j)))/(ff_JetReq->GetBinError(j)*ff_JetReq->GetBinError(j)+ggMet_JetReq->GetBinError(j)*ggMet_JetReq->GetBinError(j)));
    }
 
 }
 c1->SetLogy(1);
for(int i=5;i<50;i++){
    for(int j=1;j<=5;j++){
      gftitle="Toys/gfMet_Cut_";gftitle+=i;gftitle+="GeV_PF";
      gf = (TH1F*)f.Get(gftitle);//if(i==15)cout<<endl<<"gf original integral: "<<gf->Integral()<<endl;  
      gftitle="Toys/fgMet_Cut_";gftitle+=i;gftitle+="GeV_PF";
      TH1F* gftemp=(TH1F*)gf->Clone();
      gftemp->Add((TH1F*)f.Get(gftitle));//if(i==15)cout<<"gf integral after adding fg: "<<gftemp->Integral()<<endl;
      gftitle="Toys/gfMet_Cut_";gftitle+=i;gftitle+="GeV_PF_JetReq";
      gf_JetReq = (TH1F*)f.Get(gftitle);
      if(gftemp->Integral(0,5)>0)scale=ggMet->Integral(0,5)/gftemp->Integral(0,5);
      else scale=0.;
      if(gf_JetReq->Integral(0,5)>0)scale_JetReq=ggMet_JetReq->Integral(0,5)/gf_JetReq->Integral(0,5);
      else scale_JetReq=0.;
      gftemp->Scale(scale);gf_JetReq->Scale(scale_JetReq);//if(i==15)cout<<"gf integral after scaling to gg: "<<gftemp->Integral()<<endl;
      //if(j==1)cout<<"gfInt: "<<gf->Integral(0,5)<<"  ggInt: "<<ggMet->Integral(0,5)<<endl;;
      if(j==1 && i==5){gftemp->GetXaxis()->SetRangeUser(0,300);gftemp->Draw("histo");}
      else{gftemp->SetLineColor(4*i);gftemp->Draw("histosames");}
      PfChi25gf[i]+=(((gftemp->GetBinContent(j)-ggMet->GetBinContent(j))*(gftemp->GetBinContent(j)-ggMet->GetBinContent(j)))/(gftemp->GetBinError(j)*gftemp->GetBinError(j)+ggMet->GetBinError(j)*ggMet->GetBinError(j)));
      PfChi25gf_JetReq[i]+=(((gf_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j))*(gf_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j)))/(gf_JetReq->GetBinError(j)*gf_JetReq->GetBinError(j)+ggMet_JetReq->GetBinError(j)*ggMet_JetReq->GetBinError(j)));
      cout<<"bin:"<<j<<"  UpperCut:"<<i<<"  PfChi25gf[i]: "<<PfChi25gf[i]<<"  gf: "<<gftemp->GetBinContent(j)<<"  gg: "<<ggMet->GetBinContent(j)<<"  diff: "<<gftemp->GetBinContent(j)-ggMet->GetBinContent(j)<<" bin error gf:"<<gftemp->GetBinError(j)<<"  bin error gg:"<<ggMet->GetBinError(j)<<endl;
      
    }
    
 }
 c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gfAll.png");
 c1->SetLogy(0);
 
 for(int i=5;i<50;i++){
    for(int j=1;j<=10;j++){
      gftitle="Toys/gfMet_Cut_";gftitle+=i;gftitle+="GeV_PF";
      gf = (TH1F*)f.Get(gftitle);//if(i==15)cout<<"gf integral after getting: "<<gf->Integral()<<endl<<endl;
      gftitle="Toys/fgMet_Cut_";gftitle+=i;gftitle+="GeV_PF";
      TH1F* gftemp=(TH1F*)gf->Clone();
      gftemp->Add((TH1F*)f.Get(gftitle));  
      gftitle="Toys/gfMet_Cut_";gftitle+=i;gftitle+="GeV_PF_JetReq";
      gf_JetReq = (TH1F*)f.Get(gftitle);
      if(gftemp->Integral(0,10)>0)scale=ggMet->Integral(0,10)/gftemp->Integral(0,10);
      else scale=0.;
      if(gf_JetReq->Integral(0,10)>0)scale_JetReq=ggMet_JetReq->Integral(0,10)/gf_JetReq->Integral(0,10);
      else scale_JetReq=0.;
      gftemp->Scale(scale);gf_JetReq->Scale(scale_JetReq);
      PfChi50gf[i]+=(((gftemp->GetBinContent(j)-ggMet->GetBinContent(j))*(gftemp->GetBinContent(j)-ggMet->GetBinContent(j)))/(gftemp->GetBinError(j)*gftemp->GetBinError(j)+ggMet->GetBinError(j)*ggMet->GetBinError(j)));
      PfChi50gf_JetReq[i]+=(((gf_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j))*(gf_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j)))/(gf_JetReq->GetBinError(j)*gf_JetReq->GetBinError(j)+ggMet_JetReq->GetBinError(j)*ggMet_JetReq->GetBinError(j)));
    }
 
 }

for(int i=5;i<50;i++){
    for(int j=1;j<=5;j++){
      gammafaketitle="Toys/gammafakeMet_Cut_";gammafaketitle+=i;gammafaketitle+="GeV_PF";
      gammafake = (TH1F*)f.Get(gammafaketitle);//if(i==15)cout<<endl<<"gammafake original integral: "<<gammafake->Integral()<<endl;  
      gammafaketitle="Toys/fgMet_Cut_";gammafaketitle+=i;gammafaketitle+="GeV_PF";
      TH1F* gammafaketemp=(TH1F*)gammafake->Clone();
      gammafaketemp->Add((TH1F*)f.Get(gammafaketitle));//if(i==15)cout<<"gammafake integral after adding fg: "<<gammafaketemp->Integral()<<endl;
      gammafaketitle="Toys/gammafakeMet_Cut_";gammafaketitle+=i;gammafaketitle+="GeV_PF_JetReq";
      gammafake_JetReq = (TH1F*)f.Get(gammafaketitle);
      if(gammafaketemp->Integral(0,5)>0)scale=ggMet->Integral(0,5)/gammafaketemp->Integral(0,5);
      else scale=0.;
      if(gammafake_JetReq->Integral(0,5)>0)scale_JetReq=ggMet_JetReq->Integral(0,5)/gammafake_JetReq->Integral(0,5);
      else scale_JetReq=0.;
      gammafaketemp->Scale(scale);gammafake_JetReq->Scale(scale_JetReq);//if(i==15)cout<<"gammafake integral after scaling to gg: "<<gammafaketemp->Integral()<<endl;
      //if(j==1)cout<<"gammafakeInt: "<<gammafake->Integral(0,5)<<"  ggInt: "<<ggMet->Integral(0,5)<<endl;;
      if(j==1 && i==5){gammafaketemp->GetXaxis()->SetRangeUser(0,300);gammafaketemp->Draw("histo");}
      else{gammafaketemp->SetLineColor(4*i);gammafaketemp->Draw("histosames");}
      PfChi25gammafake[i]+=(((gammafaketemp->GetBinContent(j)-ggMet->GetBinContent(j))*(gammafaketemp->GetBinContent(j)-ggMet->GetBinContent(j)))/(gammafaketemp->GetBinError(j)*gammafaketemp->GetBinError(j)+ggMet->GetBinError(j)*ggMet->GetBinError(j)));
      PfChi25gammafake_JetReq[i]+=(((gammafake_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j))*(gammafake_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j)))/(gammafake_JetReq->GetBinError(j)*gammafake_JetReq->GetBinError(j)+ggMet_JetReq->GetBinError(j)*ggMet_JetReq->GetBinError(j)));
      cout<<"bin:"<<j<<"  UpperCut:"<<i<<"  PfChi25gammafake[i]: "<<PfChi25gammafake[i]<<"  gammafake: "<<gammafaketemp->GetBinContent(j)<<"  gg: "<<ggMet->GetBinContent(j)<<"  diff: "<<gammafaketemp->GetBinContent(j)-ggMet->GetBinContent(j)<<" bin error gammafake:"<<gammafaketemp->GetBinError(j)<<"  bin error gg:"<<ggMet->GetBinError(j)<<endl;
      
    }
    
 }
 for(int i=5;i<50;i++){
    for(int j=1;j<=10;j++){
      gammafaketitle="Toys/gammafakeMet_Cut_";gammafaketitle+=i;gammafaketitle+="GeV_PF";
      gammafake = (TH1F*)f.Get(gammafaketitle);//if(i==15)cout<<"gammafake integral after getting: "<<gammafake->Integral()<<endl<<endl;
      gammafaketitle="Toys/fgMet_Cut_";gammafaketitle+=i;gammafaketitle+="GeV_PF";
      TH1F* gammafaketemp=(TH1F*)gammafake->Clone();
      gammafaketemp->Add((TH1F*)f.Get(gammafaketitle));  
      gammafaketitle="Toys/gammafakeMet_Cut_";gammafaketitle+=i;gammafaketitle+="GeV_PF_JetReq";
      gammafake_JetReq = (TH1F*)f.Get(gammafaketitle);
      if(gammafaketemp->Integral(0,10)>0)scale=ggMet->Integral(0,10)/gammafaketemp->Integral(0,10);
      else scale=0.;
      if(gammafake_JetReq->Integral(0,10)>0)scale_JetReq=ggMet_JetReq->Integral(0,10)/gammafake_JetReq->Integral(0,10);
      else scale_JetReq=0.;
      gammafaketemp->Scale(scale);gammafake_JetReq->Scale(scale_JetReq);
      PfChi50gammafake[i]+=(((gammafaketemp->GetBinContent(j)-ggMet->GetBinContent(j))*(gammafaketemp->GetBinContent(j)-ggMet->GetBinContent(j)))/(gammafaketemp->GetBinError(j)*gammafaketemp->GetBinError(j)+ggMet->GetBinError(j)*ggMet->GetBinError(j)));
      PfChi50gammafake_JetReq[i]+=(((gammafake_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j))*(gammafake_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j)))/(gammafake_JetReq->GetBinError(j)*gammafake_JetReq->GetBinError(j)+ggMet_JetReq->GetBinError(j)*ggMet_JetReq->GetBinError(j)));
    }
 
 }

for(int i=5;i<50;i++){
    for(int j=1;j<=5;j++){
      fgtitle="Toys/fgMet_Cut_";fgtitle+=i;fgtitle+="GeV_PF";
      fg = (TH1F*)f.Get(fgtitle);  
      fgtitle="Toys/fgMet_Cut_";fgtitle+=i;fgtitle+="GeV_PF_JetReq";
      fg_JetReq = (TH1F*)f.Get(fgtitle);
      if(fg->Integral(0,5)>0)scale=ggMet->Integral(0,5)/fg->Integral(0,5);
      else scale=0.;
      if(fg_JetReq->Integral(0,5)>0)scale_JetReq=ggMet_JetReq->Integral(0,5)/fg_JetReq->Integral(0,5);
      else scale_JetReq=0.;
      fg->Scale(scale);fg_JetReq->Scale(scale_JetReq);
      PfChi25fg[i]+=(((fg->GetBinContent(j)-ggMet->GetBinContent(j))*(fg->GetBinContent(j)-ggMet->GetBinContent(j)))/(fg->GetBinError(j)*fg->GetBinError(j)+ggMet->GetBinError(j)*ggMet->GetBinError(j)));
      PfChi25fg_JetReq[i]+=(((fg_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j))*(fg_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j)))/(fg_JetReq->GetBinError(j)*fg_JetReq->GetBinError(j)+ggMet_JetReq->GetBinError(j)*ggMet_JetReq->GetBinError(j)));
      //cout<<"bin:"<<j<<"  UpperCut:"<<i<<"  PfChi25fg[i]: "<<PfChi25fg[i]<<"  fg: "<<fg->GetBinContent(j)<<"  gg: "<<ggMet->GetBinContent(j)<<"  diff: "<<fg->GetBinContent(j)-ggMet->GetBinContent(j)<<endl;

    }
    
  }
 for(int i=5;i<50;i++){
    for(int j=1;j<=10;j++){
      fgtitle="Toys/fgMet_Cut_";fgtitle+=i;fgtitle+="GeV_PF";
      fg = (TH1F*)f.Get(fgtitle);  
      fgtitle="Toys/fgMet_Cut_";fgtitle+=i;fgtitle+="GeV_PF_JetReq";
      fg_JetReq = (TH1F*)f.Get(fgtitle);
      if(fg->Integral(0,10)>0)scale=ggMet->Integral(0,10)/fg->Integral(0,10);
      else scale=0.;
      if(fg_JetReq->Integral(0,10)>0)scale_JetReq=ggMet_JetReq->Integral(0,10)/fg_JetReq->Integral(0,10);
      else scale_JetReq=0.;
      fg->Scale(scale);fg_JetReq->Scale(scale_JetReq);
      PfChi50fg[i]+=(((fg->GetBinContent(j)-ggMet->GetBinContent(j))*(fg->GetBinContent(j)-ggMet->GetBinContent(j)))/(fg->GetBinError(j)*fg->GetBinError(j)+ggMet->GetBinError(j)*ggMet->GetBinError(j)));
      PfChi50fg_JetReq[i]+=(((fg_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j))*(fg_JetReq->GetBinContent(j)-ggMet_JetReq->GetBinContent(j)))/(fg_JetReq->GetBinError(j)*fg_JetReq->GetBinError(j)+ggMet_JetReq->GetBinError(j)*ggMet_JetReq->GetBinError(j)));
    }
 
 }
 //put it all together

  TH1F* chi2_25 = new TH1F("chi2_25","",50,0,50);  
  TH1F* chi2_50 = new TH1F("chi2_50","",50,0,50);
  TH1F* chi2_25_JetReq = new TH1F("chi2_25_JetReq","",50,0,50);  
  TH1F* chi2_50_JetReq = new TH1F("chi2_50_JetReq","",50,0,50);

  TH1F* PFchi2_25 = new TH1F("Pfchi2_25","",50,0,50);  
  TH1F* PFchi2_50 = new TH1F("Pfchi2_50","",50,0,50);
  TH1F* PFchi2_25_JetReq = new TH1F("Pfchi2_25_JetReq","",50,0,50);  
  TH1F* PFchi2_50_JetReq = new TH1F("Pfchi2_50_JetReq","",50,0,50);
  
  TH1F* PFchi2_25gf = new TH1F("Pfchi2_25gf","",50,0,50);  
  TH1F* PFchi2_50gf = new TH1F("Pfchi2_50gf","",50,0,50);
  TH1F* PFchi2_25gf_JetReq = new TH1F("Pfchi2_25gf_JetReq","",50,0,50);  
  TH1F* PFchi2_50gf_JetReq = new TH1F("Pfchi2_50gf_JetReq","",50,0,50);
  
  TH1F* PFchi2_25gammafake = new TH1F("Pfchi2_25gammafake","",50,0,50);  
  TH1F* PFchi2_50gammafake = new TH1F("Pfchi2_50gammafake","",50,0,50);
  TH1F* PFchi2_25gammafake_JetReq = new TH1F("Pfchi2_25gammafake_JetReq","",50,0,50);  
  TH1F* PFchi2_50gammafake_JetReq = new TH1F("Pfchi2_50gammafake_JetReq","",50,0,50);

  TH1F* PFchi2_25fg = new TH1F("Pfchi2_25fg","",50,0,50);  
  TH1F* PFchi2_50fg = new TH1F("Pfchi2_50fg","",50,0,50);
  TH1F* PFchi2_25fg_JetReq = new TH1F("Pfchi2_25fg_JetReq","",50,0,50);  
  TH1F* PFchi2_50fg_JetReq = new TH1F("Pfchi2_50fg_JetReq","",50,0,50);

  for(int i=1;i<50;i++){
    Chi25[i]/=5;Chi25_JetReq[i]/=5;
    Chi50[i]/=10;Chi50_JetReq[i]/=10;
    chi2_25->Fill(i,Chi25[i]);
    chi2_25_JetReq->Fill(i,Chi25_JetReq[i]);
    chi2_50->Fill(i,Chi50[i]);
    chi2_50_JetReq->Fill(i,Chi50_JetReq[i]);

    PfChi25[i]/=5;PfChi25_JetReq[i]/=5;
    PfChi50[i]/=10;PfChi50_JetReq[i]/=10;
    PFchi2_25->Fill(i,PfChi25[i]);
    PFchi2_25_JetReq->Fill(i,PfChi25_JetReq[i]);
    PFchi2_50->Fill(i,PfChi50[i]);
    PFchi2_50_JetReq->Fill(i,PfChi50_JetReq[i]);

    PfChi25gf[i]/=5;PfChi25gf_JetReq[i]/=5;
    PfChi50gf[i]/=10;PfChi50gf_JetReq[i]/=10;
    PFchi2_25gf->Fill(i,PfChi25gf[i]);
    PFchi2_25gf_JetReq->Fill(i,PfChi25gf_JetReq[i]);
    PFchi2_50gf->Fill(i,PfChi50gf[i]);
    PFchi2_50gf_JetReq->Fill(i,PfChi50gf_JetReq[i]);

    PfChi25gammafake[i]/=5;PfChi25gammafake_JetReq[i]/=5;
    PfChi50gammafake[i]/=10;PfChi50gammafake_JetReq[i]/=10;
    PFchi2_25gammafake->Fill(i,PfChi25gammafake[i]);
    PFchi2_25gammafake_JetReq->Fill(i,PfChi25gammafake_JetReq[i]);
    PFchi2_50gammafake->Fill(i,PfChi50gammafake[i]);
    PFchi2_50gammafake_JetReq->Fill(i,PfChi50gammafake_JetReq[i]);
  
    PfChi25fg[i]/=5;PfChi25fg_JetReq[i]/=5;
    PfChi50fg[i]/=10;PfChi50fg_JetReq[i]/=10;
    PFchi2_25fg->Fill(i,PfChi25fg[i]);
    PFchi2_25fg_JetReq->Fill(i,PfChi25fg_JetReq[i]);
    PFchi2_50fg->Fill(i,PfChi50fg[i]);
    PFchi2_50fg_JetReq->Fill(i,PfChi50fg_JetReq[i]);
  }

  chi2_25->GetYaxis()->SetTitle("#Chi^{2}");
  chi2_25->GetXaxis()->SetTitle("Upper Cut on ff DR03 Combined Isolation (GeV)");
  chi2_25->GetXaxis()->SetTitleSize(.04);
  chi2_25->GetXaxis()->SetTitleOffset(1.3);

  chi2_25_JetReq->GetYaxis()->SetTitle("#Chi^{2}");
  chi2_25_JetReq->GetXaxis()->SetTitle("Upper Cut on ff DR03 Combined Isolation (GeV)");
  chi2_25_JetReq->GetXaxis()->SetTitleSize(.04);
  chi2_25_JetReq->GetXaxis()->SetTitleOffset(1.3);

  chi2_25->SetLineColor(kBlue);
  chi2_25->Draw("L");
  chi2_50->SetLineColor(kRed);
  chi2_50->Draw("LSAME");
  TLegend* leg = new TLegend(.5,.5,.8,.7,"","brNDC");
  leg->AddEntry(chi2_50,"E_{T}^{miss} < 50GeV","l"); 
  leg->AddEntry(chi2_25,"E_{T}^{miss} < 25GeV","l"); 
  leg->SetFillColor(kWhite);
  leg->Draw("SAME");
  Text->Draw();
  //c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure8_Chi2.png");
  //c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure8_Chi2.pdf");

  chi2_25_JetReq->SetLineColor(kBlue);
  chi2_25_JetReq->Draw("L");
  chi2_50_JetReq->SetLineColor(kRed);
  chi2_50_JetReq->Draw("LSAME");
  TLegend* leg_JetReq = new TLegend(.4,.5,.8,.7,"","brNDC");
  leg_JetReq->AddEntry(chi2_25,"E_{T}^{miss} < 25GeV JetReq","l"); 
  leg_JetReq->AddEntry(chi2_50,"E_{T}^{miss} < 50GeV JetReq","l"); 
  leg_JetReq->SetFillColor(kWhite);
  leg_JetReq->Draw("SAME");
  Text->Draw();
  //c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure8_Chi2_JetReq.png");
  //c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure8_Chi2_JetReq.pdf");

  chi2_25->GetXaxis()->SetRangeUser(0,50);
  chi2_25->GetYaxis()->SetRangeUser(0,12);
  chi2_25->Draw("L");
  chi2_50->Draw("LSAME");
  TLegend* leg = new TLegend(.28,.5,.58,.7,"","brNDC");
  leg->AddEntry(chi2_50,"E_{T}^{miss} < 50GeV","l"); 
  leg->AddEntry(chi2_25,"E_{T}^{miss} < 25GeV","l"); 
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(4000);
  leg->Draw("SAME");
  Text->Draw();
  //c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure8_Chi2_ZOOM.png");
  //c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure8_Chi2_ZOOM.pdf");

  chi2_25_JetReq->GetXaxis()->SetRangeUser(0,50);
  chi2_25_JetReq->GetYaxis()->SetRangeUser(0.8,3.6);
  chi2_25_JetReq->Draw("L");
  chi2_50_JetReq->Draw("LSAME");
  TLegend* leg_JetReq = new TLegend(.45,.15,.85,.35,"","brNDC");
  leg_JetReq->AddEntry(chi2_25,"E_{T}^{miss} < 25GeV JetReq","l"); 
  leg_JetReq->AddEntry(chi2_50,"E_{T}^{miss} < 50GeV JetReq","l"); 
  leg_JetReq->SetFillColor(kWhite);
  leg_JetReq->Draw("SAME");
  Text->Draw();
  //c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure8_Chi2_ZOOM_JetReq.png");
  //c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/Figure8_Chi2_ZOOM_JetReq.pdf");

  //now pf
  PFchi2_50->GetYaxis()->SetTitle("#Chi^{2}/ndof");
  PFchi2_50->GetXaxis()->SetTitle("Upper Cut on fake PF Combined Isolation (GeV)");
  PFchi2_50->GetXaxis()->SetTitleSize(.04);
  PFchi2_50->GetXaxis()->SetTitleOffset(1.3);

  PFchi2_25_JetReq->GetYaxis()->SetTitle("#Chi^{2}");
  PFchi2_25_JetReq->GetXaxis()->SetTitle("Upper Cut on ff PF Combined Isolation (GeV)");
  PFchi2_25_JetReq->GetXaxis()->SetTitleSize(.04);
  PFchi2_25_JetReq->GetXaxis()->SetTitleOffset(1.3);
  
  PFchi2_50->SetLineColor(kRed);
  PFchi2_25->SetLineColor(kBlue);
  PFchi2_50->GetXaxis()->SetRangeUser(0,50);
  PFchi2_50->Draw("L");
  PFchi2_25->Draw("Lsame");
  TLegend* leg = new TLegend(.53,.36,.83,.56,"","brNDC");
  leg->AddEntry(PFchi2_50,"E_{T}^{miss} < 50GeV","l"); 
  leg->AddEntry(PFchi2_25,"E_{T}^{miss} < 25GeV","l"); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2.pdf");

  
  PFchi2_50_JetReq->SetLineColor(kRed);
  PFchi2_25_JetReq->SetLineColor(kBlue);
  PFchi2_50_JetReq->GetXaxis()->SetRangeUser(0,50);
  PFchi2_50_JetReq->Draw("L");
  PFchi2_25_JetReq->Draw("LSAME");
  TLegend* leg_JetReq = new TLegend(.45,.5,.85,.7,"","brNDC");
  leg_JetReq->AddEntry(PFchi2_25,"E_{T}^{miss} < 25GeV JetReq","l"); 
  leg_JetReq->AddEntry(PFchi2_50,"E_{T}^{miss} < 50GeV JetReq","l"); 
  leg_JetReq->SetFillColor(kWhite);
  leg_JetReq->SetBorderSize(0);
  leg_JetReq->SetFillStyle(0);
  leg_JetReq->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2_JetReq.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2_JetReq.pdf");

  PFchi2_50->GetXaxis()->SetRangeUser(0,50);
  //PFchi2_25->GetYaxis()->SetRangeUser(0,3.2);
  PFchi2_50->Draw("L");
  PFchi2_25->Draw("LSAME");
  TLegend* leg = new TLegend(.59,.15,.89,.35,"","brNDC");
  leg->AddEntry(PFchi2_50,"E_{T}^{miss} < 50GeV","l"); 
  leg->AddEntry(PFchi2_25,"E_{T}^{miss} < 25GeV","l"); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2_ZOOM.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2_ZOOM.pdf");

  PFchi2_50_JetReq->GetXaxis()->SetRangeUser(0,50);
  //PFchi2_25_JetReq->GetYaxis()->SetRangeUser(0,4);
  PFchi2_50_JetReq->Draw("L");
  PFchi2_25_JetReq->Draw("LSAME");
  TLegend* leg_JetReq = new TLegend(.5,.15,.85,.35,"","brNDC");
  leg_JetReq->AddEntry(PFchi2_25,"E_{T}^{miss} < 25GeV JetReq","l"); 
  leg_JetReq->AddEntry(PFchi2_50,"E_{T}^{miss} < 50GeV JetReq","l"); 
  leg_JetReq->SetFillColor(kWhite);
  leg_JetReq->SetBorderSize(0);
  leg_JetReq->SetFillStyle(0);
  leg_JetReq->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2_ZOOM_JetReq.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2_ZOOM_JetReq.pdf");

  //now gf

  PFchi2_25gf->GetYaxis()->SetTitle("#Chi^{2}");
  PFchi2_25gf->GetXaxis()->SetTitle("Upper Cut on gf PF Combined Isolation (GeV)");
  PFchi2_25gf->GetXaxis()->SetTitleSize(.04);
  PFchi2_25gf->GetXaxis()->SetTitleOffset(1.3);

  PFchi2_25gf_JetReq->GetYaxis()->SetTitle("#Chi^{2}");
  PFchi2_25gf_JetReq->GetXaxis()->SetTitle("Upper Cut on gf PF Combined Isolation (GeV)");
  PFchi2_25gf_JetReq->GetXaxis()->SetTitleSize(.04);
  PFchi2_25gf_JetReq->GetXaxis()->SetTitleOffset(1.3);
  
  PFchi2_50gf->SetLineColor(kRed);
  PFchi2_25gf->SetLineColor(kBlue);
  PFchi2_50gf->GetXaxis()->SetRangeUser(0,50);
  PFchi2_50gf->Draw("L");
  PFchi2_25gf->Draw("Lsame");
  TLegend* leg = new TLegend(.5,.2,.8,.4,"","brNDC");
  leg->AddEntry(PFchi2_50gf,"E_{T}^{miss} < 50GeV","l"); 
  leg->AddEntry(PFchi2_25gf,"E_{T}^{miss} < 25GeV","l"); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gf.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gf.pdf");

  
  PFchi2_50gf_JetReq->SetLineColor(kRed);
  PFchi2_25gf_JetReq->SetLineColor(kBlue);
  PFchi2_50gf_JetReq->GetXaxis()->SetRangeUser(0,50);
  PFchi2_50gf_JetReq->Draw("L");
  PFchi2_25gf_JetReq->Draw("LSAME");
  TLegend* leg_JetReq = new TLegend(.45,.5,.85,.7,"","brNDC");
  leg_JetReq->AddEntry(PFchi2_25gf,"E_{T}^{miss} < 25GeV JetReq","l"); 
  leg_JetReq->AddEntry(PFchi2_50gf,"E_{T}^{miss} < 50GeV JetReq","l"); 
  leg_JetReq->SetFillColor(kWhite);
  leg_JetReq->SetBorderSize(0);
  leg_JetReq->SetFillStyle(0);
  leg_JetReq->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gf_JetReq.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gf_JetReq.pdf");

  PFchi2_50gf->GetXaxis()->SetRangeUser(0,50);
  //PFchi2_25gf->GetYaxis()->SetRangeUser(0,3.2);
  PFchi2_50gf->Draw("L");
  PFchi2_25gf->Draw("LSAME");
  TLegend* leg = new TLegend(.59,.15,.89,.35,"","brNDC");
  leg->AddEntry(PFchi2_50gf,"E_{T}^{miss} < 50GeV","l"); 
  leg->AddEntry(PFchi2_25gf,"E_{T}^{miss} < 25GeV","l"); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gf_ZOOM.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gf_ZOOM.pdf");

  PFchi2_50gf_JetReq->GetXaxis()->SetRangeUser(0,50);
  //PFchi2_25gf_JetReq->GetYaxis()->SetRangeUser(0,4);
  PFchi2_50gf_JetReq->Draw("L");
  PFchi2_25gf_JetReq->Draw("LSAME");
  TLegend* leg_JetReq = new TLegend(.5,.15,.85,.35,"","brNDC");
  leg_JetReq->AddEntry(PFchi2_25gf,"E_{T}^{miss} < 25GeV JetReq","l"); 
  leg_JetReq->AddEntry(PFchi2_50gf,"E_{T}^{miss} < 50GeV JetReq","l"); 
  leg_JetReq->SetFillColor(kWhite);
  leg_JetReq->SetBorderSize(0);
  leg_JetReq->SetFillStyle(0);
  leg_JetReq->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gf_ZOOM_JetReq.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gf_ZOOM_JetReq.pdf");

  //now gammafake

  PFchi2_25gammafake->GetYaxis()->SetTitle("#Chi^{2}");
  PFchi2_25gammafake->GetXaxis()->SetTitle("Upper Cut on gammafake PF Combined Isolation (GeV)");
  PFchi2_25gammafake->GetXaxis()->SetTitleSize(.04);
  PFchi2_25gammafake->GetXaxis()->SetTitleOffset(1.3);

  PFchi2_25gammafake_JetReq->GetYaxis()->SetTitle("#Chi^{2}");
  PFchi2_25gammafake_JetReq->GetXaxis()->SetTitle("Upper Cut on gammafake PF Combined Isolation (GeV)");
  PFchi2_25gammafake_JetReq->GetXaxis()->SetTitleSize(.04);
  PFchi2_25gammafake_JetReq->GetXaxis()->SetTitleOffset(1.3);
  
  PFchi2_50gammafake->SetLineColor(kRed);
  PFchi2_25gammafake->SetLineColor(kBlue);
  PFchi2_50gammafake->GetXaxis()->SetRangeUser(0,50);
  PFchi2_50gammafake->Draw("L");
  PFchi2_25gammafake->Draw("Lsame");
  TLegend* leg = new TLegend(.5,.2,.8,.4,"","brNDC");
  leg->AddEntry(PFchi2_50gammafake,"E_{T}^{miss} < 50GeV","l"); 
  leg->AddEntry(PFchi2_25gammafake,"E_{T}^{miss} < 25GeV","l"); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gammafake.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gammafake.pdf");

  
  PFchi2_50gammafake_JetReq->SetLineColor(kRed);
  PFchi2_25gammafake_JetReq->SetLineColor(kBlue);
  PFchi2_50gammafake_JetReq->GetXaxis()->SetRangeUser(0,50);
  PFchi2_50gammafake_JetReq->Draw("L");
  PFchi2_25gammafake_JetReq->Draw("LSAME");
  TLegend* leg_JetReq = new TLegend(.45,.5,.85,.7,"","brNDC");
  leg_JetReq->AddEntry(PFchi2_25gammafake,"E_{T}^{miss} < 25GeV JetReq","l"); 
  leg_JetReq->AddEntry(PFchi2_50gammafake,"E_{T}^{miss} < 50GeV JetReq","l"); 
  leg_JetReq->SetFillColor(kWhite);
  leg_JetReq->SetBorderSize(0);
  leg_JetReq->SetFillStyle(0);
  leg_JetReq->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gammafake_JetReq.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gammafake_JetReq.pdf");

  PFchi2_50gammafake->GetXaxis()->SetRangeUser(0,50);
  //PFchi2_25gammafake->GetYaxis()->SetRangeUser(0,3.2);
  PFchi2_50gammafake->Draw("L");
  PFchi2_25gammafake->Draw("LSAME");
  TLegend* leg = new TLegend(.59,.15,.89,.35,"","brNDC");
  leg->AddEntry(PFchi2_50gammafake,"E_{T}^{miss} < 50GeV","l"); 
  leg->AddEntry(PFchi2_25gammafake,"E_{T}^{miss} < 25GeV","l"); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gammafake_ZOOM.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gammafake_ZOOM.pdf");

  PFchi2_50gammafake_JetReq->GetXaxis()->SetRangeUser(0,50);
  //PFchi2_25gammafake_JetReq->GetYaxis()->SetRangeUser(0,4);
  PFchi2_50gammafake_JetReq->Draw("L");
  PFchi2_25gammafake_JetReq->Draw("LSAME");
  TLegend* leg_JetReq = new TLegend(.5,.15,.85,.35,"","brNDC");
  leg_JetReq->AddEntry(PFchi2_25gammafake,"E_{T}^{miss} < 25GeV JetReq","l"); 
  leg_JetReq->AddEntry(PFchi2_50gammafake,"E_{T}^{miss} < 50GeV JetReq","l"); 
  leg_JetReq->SetFillColor(kWhite);
  leg_JetReq->SetBorderSize(0);
  leg_JetReq->SetFillStyle(0);
  leg_JetReq->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gammafake_ZOOM_JetReq.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2gammafake_ZOOM_JetReq.pdf");

//now fg

  PFchi2_25fg->GetYaxis()->SetTitle("#Chi^{2}");
  PFchi2_25fg->GetXaxis()->SetTitle("Upper Cut on fg PF Combined Isolation (GeV)");
  PFchi2_25fg->GetXaxis()->SetTitleSize(.04);
  PFchi2_25fg->GetXaxis()->SetTitleOffset(1.3);

  PFchi2_25fg_JetReq->GetYaxis()->SetTitle("#Chi^{2}");
  PFchi2_25fg_JetReq->GetXaxis()->SetTitle("Upper Cut on fg PF Combined Isolation (GeV)");
  PFchi2_25fg_JetReq->GetXaxis()->SetTitleSize(.04);
  PFchi2_25fg_JetReq->GetXaxis()->SetTitleOffset(1.3);
  
  PFchi2_50fg->SetLineColor(kRed);
  PFchi2_25fg->SetLineColor(kBlue);
  PFchi2_50fg->GetXaxis()->SetRangeUser(0,50);
  PFchi2_50fg->Draw("L");
  PFchi2_25fg->Draw("Lsame");
  TLegend* leg = new TLegend(.5,.2,.8,.4,"","brNDC");
  leg->AddEntry(PFchi2_50fg,"E_{T}^{miss} < 50GeV","l"); 
  leg->AddEntry(PFchi2_25fg,"E_{T}^{miss} < 25GeV","l"); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2fg.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2fg.pdf");

  
  PFchi2_50fg_JetReq->SetLineColor(kRed);
  PFchi2_25fg_JetReq->SetLineColor(kBlue);
  PFchi2_50fg_JetReq->GetXaxis()->SetRangeUser(0,50);
  PFchi2_50fg_JetReq->Draw("L");
  PFchi2_25fg_JetReq->Draw("LSAME");
  TLegend* leg_JetReq = new TLegend(.45,.5,.85,.7,"","brNDC");
  leg_JetReq->AddEntry(PFchi2_25fg,"E_{T}^{miss} < 25GeV JetReq","l"); 
  leg_JetReq->AddEntry(PFchi2_50fg,"E_{T}^{miss} < 50GeV JetReq","l"); 
  leg_JetReq->SetFillColor(kWhite);
  leg_JetReq->SetBorderSize(0);
  leg_JetReq->SetFillStyle(0);
  leg_JetReq->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2fg_JetReq.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2fg_JetReq.pdf");

  PFchi2_50fg->GetXaxis()->SetRangeUser(0,50);
  //PFchi2_25fg->GetYaxis()->SetRangeUser(0,3.2);
  PFchi2_50fg->Draw("L");
  PFchi2_25fg->Draw("LSAME");
  TLegend* leg = new TLegend(.59,.15,.89,.35,"","brNDC");
  leg->AddEntry(PFchi2_50fg,"E_{T}^{miss} < 50GeV","l"); 
  leg->AddEntry(PFchi2_25fg,"E_{T}^{miss} < 25GeV","l"); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2fg_ZOOM.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2fg_ZOOM.pdf");

  PFchi2_50fg_JetReq->GetXaxis()->SetRangeUser(0,50);
  //PFchi2_25fg_JetReq->GetYaxis()->SetRangeUser(0,4);
  PFchi2_50fg_JetReq->Draw("L");
  PFchi2_25fg_JetReq->Draw("LSAME");
  TLegend* leg_JetReq = new TLegend(.5,.15,.85,.35,"","brNDC");
  leg_JetReq->AddEntry(PFchi2_25fg,"E_{T}^{miss} < 25GeV JetReq","l"); 
  leg_JetReq->AddEntry(PFchi2_50fg,"E_{T}^{miss} < 50GeV JetReq","l"); 
  leg_JetReq->SetFillColor(kWhite);
  leg_JetReq->SetBorderSize(0);
  leg_JetReq->SetFillStyle(0);
  leg_JetReq->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2fg_ZOOM_JetReq.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2fg_ZOOM_JetReq.pdf");

  //PFchi2_50->GetXaxis()->SetRangeUser(0,50);PFchi2_50->GetYaxis()->SetRangeUser(0,25);
  //PFchi2_50->Draw("L");
  PFchi2_25->Draw("L");PFchi2_25->GetYaxis()->SetRangeUser(0,25);
  PFchi2_25gammafake->SetLineColor(kCyan);PFchi2_50gammafake->SetLineColor(kOrange);
  //PFchi2_50gammafake->Draw("Lsame");
  PFchi2_25gammafake->Draw("Lsame");
  //PFchi2_25fg->SetLineColor(kGreen);PFchi2_50fg->SetLineColor(kViolet);
  //PFchi2_50fg->Draw("Lsame");
  //PFchi2_25fg->Draw("Lsame");
  TLegend* leg2 = new TLegend(.21,.47,.52,.69,"","brNDC");
  //leg2->AddEntry(PFchi2_50,"ff  - E_{T}^{miss} < 50GeV","l"); 
  leg2->AddEntry(PFchi2_25,"ff  - E_{T}^{miss} < 25GeV","l"); 
  //leg2->AddEntry(PFchi2_50gammafake,"gammafake - E_{T}^{miss} < 50GeV","l"); 
  leg2->AddEntry(PFchi2_25gammafake,"#gammafake - E_{T}^{miss} < 25GeV","l"); 
  //leg2->AddEntry(PFchi2_50fg,"fg - E_{T}^{miss} < 50GeV","l"); 
  //leg2->AddEntry(PFchi2_25fg,"fg - E_{T}^{miss} < 25GeV","l"); 
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2_ff_and_gammafake.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2_ff_and_gammafake.pdf");

  // PFchi2_50_JetReq->GetXaxis()->SetRangeUser(0,50);PFchi2_50_JetReq->GetYaxis()->SetRangeUser(0,25);
  //PFchi2_50_JetReq->Draw("L");
  PFchi2_25_JetReq->Draw("L");PFchi2_25_JetReq->GetYaxis()->SetRangeUser(0,10);
  PFchi2_25gammafake_JetReq->SetLineColor(kCyan);PFchi2_50gammafake_JetReq->SetLineColor(kOrange);
  // PFchi2_50gammafake_JetReq->Draw("Lsame");
  PFchi2_25gammafake_JetReq->Draw("Lsame");
  //PFchi2_25fg_JetReq->SetLineColor(kGreen);PFchi2_50fg_JetReq->SetLineColor(kViolet);
  //PFchi2_50fg_JetReq->Draw("Lsame");
  //PFchi2_25fg_JetReq->Draw("Lsame");
  TLegend* leg2 = new TLegend(.21,.47,.52,.69,"","brNDC");
  //leg2->AddEntry(PFchi2_50,"ff  - E_{T}^{miss} < 50GeV JetReq","l"); 
  leg2->AddEntry(PFchi2_25,"ff  - E_{T}^{miss} < 25GeV JetReq","l"); 
  //leg2->AddEntry(PFchi2_50gammafake,"gammafake - E_{T}^{miss} < 50GeV JetReq","l"); 
  leg2->AddEntry(PFchi2_25gammafake,"#gammafake - E_{T}^{miss} < 25GeV JetReq","l"); 
  //leg2->AddEntry(PFchi2_50fg,"fg - E_{T}^{miss} < 50GeV JetReq","l"); 
  //leg2->AddEntry(PFchi2_25fg,"fg - E_{T}^{miss} < 25GeV JetReq","l"); 
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2_ff_and_gammafake_JetReq.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2_ff_and_gammafake_JetReq.pdf");

  PFchi2_50->GetXaxis()->SetRangeUser(0,25);PFchi2_50->GetYaxis()->SetRangeUser(0,25);
  PFchi2_50->Draw("L");
  PFchi2_25->Draw("Lsame");
  PFchi2_25gammafake->SetLineColor(kCyan);PFchi2_50gammafake->SetLineColor(kOrange);
  PFchi2_50gammafake->Draw("Lsame");
  PFchi2_25gammafake->Draw("Lsame");
  //PFchi2_25fg->SetLineColor(kGreen);PFchi2_50fg->SetLineColor(kViolet);
  //PFchi2_50fg->Draw("Lsame");
  //PFchi2_25fg->Draw("Lsame");
  leg2->Draw("SAME");
  Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2_ff_and_gammafake_ZOOM.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2_ff_and_gammafake_ZOOM.pdf");
 
  PFchi2_50->GetXaxis()->SetRangeUser(0,10);PFchi2_50->GetYaxis()->SetRangeUser(0,25);
  PFchi2_50->Draw("L");
  PFchi2_25->Draw("Lsame");
  PFchi2_25gammafake->SetLineColor(kCyan);PFchi2_50gammafake->SetLineColor(kOrange);
  PFchi2_50gammafake->Draw("Lsame");
  PFchi2_25gammafake->Draw("Lsame");
  //PFchi2_25fg->SetLineColor(kGreen);PFchi2_50fg->SetLineColor(kViolet);
  //PFchi2_50fg->Draw("Lsame");
  //PFchi2_25fg->Draw("Lsame");
  //leg2->Draw("SAME");
  //Text->Draw();
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2_ff_and_gammafake_SuperZOOM.png");
  c1->Print("Plots/LooseWP/15GeVFakeCut/sihih012/pixelVeto/19499pb/NewDiEMPtBins/NewJetMatch/TwoJetReqApr24/PFchi2_ff_and_gammafake_SuperZOOM.png");

}


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
void JEC(){

  gStyle->SetOptStat(0);

  Double_t MTbins[]={0,30,60,90,175};
  int nMTbins=(sizeof(MTbins)/sizeof(Double_t))-1;
  float MTplotXmax = MTbins[nMTbins];

  Double_t Metbins[]={0,20,40,65,100,150,250};
  int nMetbins=(sizeof(Metbins)/sizeof(Double_t))-1;
  float MetplotXmax = Metbins[nMetbins];

  TFile f("hist_HiggsAna_Photon_SMS_TChiWH_WincHgg_2J.root","READ");
  
  TCanvas *c1 = new TCanvas("c1","c1",800,600);c1->cd();
  c1->SetLogy(0);

  TH2F* mt_ele_Vinvmass = (TH2F*)f.Get("ggMTvsInvarMass_Loose_1Ele_0_1Jets");
  TH1F* mt_ele= (TH1F*)mt_ele_Vinvmass->ProjectionY("mt_ele");
  //mt_ele->Rebin(50);
  mt_ele->SetTitle("");
  mt_ele->SetMarkerSize(0.5);

  TH1F* mt_ele_NoJEC = f.Get("ggMT_NoJEC_Loose_1Ele_0_1Jets");
  //mt_ele_NoJEC->Rebin(50);
  mt_ele_NoJEC->SetLineColor(kGreen);

  TH1F* mt_ele_JECup = f.Get("ggMT_JECup_Loose_1Ele_0_1Jets");
  //mt_ele_JECup->Rebin(50);
  mt_ele_JECup->SetLineColor(kRed);
  mt_ele_JECup->SetMarkerColor(kRed);
  mt_ele_JECup->SetMarkerSize(0.5);

  TH1F* mt_ele_JECdown = f.Get("ggMT_JECdown_Loose_1Ele_0_1Jets");
  //mt_ele_JECdown->Rebin(50);
  mt_ele_JECdown->SetLineColor(kBlue);
  mt_ele_JECdown->SetMarkerColor(kBlue);
  mt_ele_JECdown->SetMarkerSize(0.5);
  cout<<"Events:  raw: "<<mt_ele->Integral(0,-1)<<"  NoJEC: "<<mt_ele_NoJEC->Integral(0,-1)<<"  JECup: "<<mt_ele_JECup->Integral(0,-1)<<"  JECdown: "<<mt_ele_JECdown->Integral(0,-1)<<endl;

  TH1F* mt_ele_Rebin = (TH1F*)mt_ele->Rebin(nMTbins,"mt_ele_Rebin",MTbins);
  TH1F* mt_ele_NoJECRebin = (TH1F*)mt_ele_NoJEC->Rebin(nMTbins,"mt_ele_NoJECRebin",MTbins);
  TH1F* mt_ele_JECupRebin = (TH1F*)mt_ele_JECup->Rebin(nMTbins,"mt_ele_JECupRebin",MTbins);
  TH1F* mt_ele_JECdownRebin = (TH1F*)mt_ele_JECdown->Rebin(nMTbins,"mt_ele_JECdownRebin",MTbins);

  cout<<"Events Rebin:  raw: "<<mt_ele_Rebin->Integral(0,-1)<<"  NoJEC: "<<mt_ele_NoJECRebin->Integral(0,-1)<<"  JECup: "<<mt_ele_JECupRebin->Integral(0,-1)<<"  JECdown: "<<mt_ele_JECdownRebin->Integral(0,-1)<<endl;


  AddOverflowToLastBin(mt_ele_Rebin);
  AddOverflowToLastBin(mt_ele_NoJECRebin);
  AddOverflowToLastBin(mt_ele_JECupRebin);
  AddOverflowToLastBin(mt_ele_JECdownRebin);

  mt_ele_Rebin->GetXaxis()->SetTitle("M_{T} [GeV]"); mt_ele_Rebin->GetXaxis()->SetTitleOffset(1.);
  //mt_ele_Rebin->GetYaxis()->SetRangeUser(300,2380);
  //mt_ele_Rebin->GetXaxis()->SetRangeUser(120,121);
  mt_ele_Rebin->Draw("pe");
  //mt_ele_NoJECRebin->Draw("pesames");
  mt_ele_JECupRebin->Draw("pesames");
  mt_ele_JECdownRebin->Draw("pesames");
  
  c1->Print("Plots/Higgs/JEC_up_down_mt_ele.png");
  c1->Print("Plots/Higgs/JEC_up_down_mt_ele.pdf");

  for(int i=1;i<=mt_ele_Rebin->GetNbinsX();i++){
    if(i==1)cout<<"ele mt"<<endl<<"No JEC  "<<endl;
    cout<<"bin:"<<i<<" events: "<<mt_ele_Rebin->GetBinContent(i)<<endl;
  }
  for(int i=1;i<=mt_ele_Rebin->GetNbinsX();i++){
    if(i==1)cout<<"JECup   "<<endl;
    cout<<"bin:"<<i<<" events: "<<mt_ele_JECupRebin->GetBinContent(i)<<endl;
  }
  for(int i=1;i<=mt_ele_Rebin->GetNbinsX();i++){
    if(i==1)cout<<"JEC down"<<endl;
    cout<<"bin:"<<i<<" events: "<<mt_ele_JECdownRebin->GetBinContent(i)<<endl;
  }
  
  cout<<"Events Rebin after Overflow add:  raw: "<<mt_ele_Rebin->Integral(0,-1)<<"  NoJEC: "<<mt_ele_NoJECRebin->Integral(0,-1)<<"  JECup: "<<mt_ele_JECupRebin->Integral(0,-1)<<"  JECdown: "<<mt_ele_JECdownRebin->Integral(0,-1)<<endl;

  //now muon

  TH2F* mt_mu_Vinvmass = (TH2F*)f.Get("ggMTvsInvarMass_Loose_1Mu_0_1Jets");
  TH1F* mt_mu= (TH1F*)mt_mu_Vinvmass->ProjectionY("mt_mu");
  //mt_mu->Rebin(50);
  mt_mu->SetTitle("");
  mt_mu->SetMarkerSize(0.5);

  TH1F* mt_mu_NoJEC = f.Get("ggMT_NoJEC_Loose_1Mu_0_1Jets");
  //mt_mu_NoJEC->Rebin(50);
  mt_mu_NoJEC->SetLineColor(kGreen);

  TH1F* mt_mu_JECup = f.Get("ggMT_JECup_Loose_1Mu_0_1Jets");
  //mt_mu_JECup->Rebin(50);
  mt_mu_JECup->SetLineColor(kRed);
  mt_mu_JECup->SetMarkerColor(kRed);
  mt_mu_JECup->SetMarkerSize(0.5);

  TH1F* mt_mu_JECdown = f.Get("ggMT_JECdown_Loose_1Mu_0_1Jets");
  //mt_mu_JECdown->Rebin(50);
  mt_mu_JECdown->SetLineColor(kBlue);
  mt_mu_JECdown->SetMarkerColor(kBlue);
  mt_mu_JECdown->SetMarkerSize(0.5);
  cout<<"Events:  raw: "<<mt_mu->Integral(0,-1)<<"  NoJEC: "<<mt_mu_NoJEC->Integral(0,-1)<<"  JECup: "<<mt_mu_JECup->Integral(0,-1)<<"  JECdown: "<<mt_mu_JECdown->Integral(0,-1)<<endl;

  TH1F* mt_mu_Rebin = (TH1F*)mt_mu->Rebin(nMTbins,"mt_mu_Rebin",MTbins);
  TH1F* mt_mu_NoJECRebin = (TH1F*)mt_mu_NoJEC->Rebin(nMTbins,"mt_mu_NoJECRebin",MTbins);
  TH1F* mt_mu_JECupRebin = (TH1F*)mt_mu_JECup->Rebin(nMTbins,"mt_mu_JECupRebin",MTbins);
  TH1F* mt_mu_JECdownRebin = (TH1F*)mt_mu_JECdown->Rebin(nMTbins,"mt_mu_JECdownRebin",MTbins);

  cout<<"Events Rebin:  raw: "<<mt_mu_Rebin->Integral(0,-1)<<"  NoJEC: "<<mt_mu_NoJECRebin->Integral(0,-1)<<"  JECup: "<<mt_mu_JECupRebin->Integral(0,-1)<<"  JECdown: "<<mt_mu_JECdownRebin->Integral(0,-1)<<endl;


  AddOverflowToLastBin(mt_mu_Rebin);
  AddOverflowToLastBin(mt_mu_NoJECRebin);
  AddOverflowToLastBin(mt_mu_JECupRebin);
  AddOverflowToLastBin(mt_mu_JECdownRebin);

  mt_mu_Rebin->GetXaxis()->SetTitle("M_{T} [GeV]"); mt_mu_Rebin->GetXaxis()->SetTitleOffset(1.);
  //mt_mu_Rebin->GetYaxis()->SetRangeUser(450,590);
  //mt_mu_Rebin->GetXaxis()->SetRangeUser(0,80);
  mt_mu_Rebin->Draw("pe");
  //mt_mu_NoJECRebin->Draw("pesames");
  mt_mu_JECupRebin->Draw("pesames");
  mt_mu_JECdownRebin->Draw("pesames");
  
  c1->Print("Plots/Higgs/JEC_up_down_mt_mu.png");
  c1->Print("Plots/Higgs/JEC_up_down_mt_mu.pdf");

  for(int i=1;i<=mt_mu_Rebin->GetNbinsX();i++){
    if(i==1)cout<<"mu mt"<<endl<<"No JEC  "<<endl;
    cout<<"bin:"<<i<<" events: "<<mt_mu_Rebin->GetBinContent(i)<<endl;
  }
  for(int i=1;i<=mt_mu_Rebin->GetNbinsX();i++){
    if(i==1)cout<<"JECup   "<<endl;
    cout<<"bin:"<<i<<" events: "<<mt_mu_JECupRebin->GetBinContent(i)<<endl;
  }
  for(int i=1;i<=mt_mu_Rebin->GetNbinsX();i++){
    if(i==1)cout<<"JEC down"<<endl;
    cout<<"bin:"<<i<<" events: "<<mt_mu_JECdownRebin->GetBinContent(i)<<endl;
  }
  
  cout<<"Events Rebin after Overflow add:  raw: "<<mt_mu_Rebin->Integral(0,-1)<<"  NoJEC: "<<mt_mu_NoJECRebin->Integral(0,-1)<<"  JECup: "<<mt_mu_JECupRebin->Integral(0,-1)<<"  JECdown: "<<mt_mu_JECdownRebin->Integral(0,-1)<<endl;

  //Now MET
  cout<<"Now do MET"<<endl;
  TH2F* met_ele_Vinvmass = (TH2F*)f.Get("ggMetVsInvarMass_Loose_1Ele_0_1Jets");
  TH1F* met_ele= (TH1F*)met_ele_Vinvmass->ProjectionY("met_ele");
  //met_ele->Rebin(50);
  met_ele->SetTitle("");
  met_ele->SetMarkerSize(0.5);

  TH1F* met_ele_NoJEC = f.Get("ggMet_NoJEC_Loose_1Ele_0_1Jets");
  //met_ele_NoJEC->Rebin(50);
  met_ele_NoJEC->SetLineColor(kGreen);

  TH1F* met_ele_JECup = f.Get("ggMet_JECup_Loose_1Ele_0_1Jets");
  //met_ele_JECup->Rebin(50);
  met_ele_JECup->SetLineColor(kRed);
  met_ele_JECup->SetMarkerColor(kRed);
  met_ele_JECup->SetMarkerSize(0.5);

  TH1F* met_ele_JECdown = f.Get("ggMet_JECdown_Loose_1Ele_0_1Jets");
  //met_ele_JECdown->Rebin(50);
  met_ele_JECdown->SetLineColor(kBlue);
  met_ele_JECdown->SetMarkerColor(kBlue);
  met_ele_JECdown->SetMarkerSize(0.5);
  cout<<"Events:  raw: "<<met_ele->Integral(0,-1)<<"  NoJEC: "<<met_ele_NoJEC->Integral(0,-1)<<"  JECup: "<<met_ele_JECup->Integral(0,-1)<<"  JECdown: "<<met_ele_JECdown->Integral(0,-1)<<endl;

  TH1F* met_ele_Rebin = (TH1F*)met_ele->Rebin(nMetbins,"met_ele_Rebin",Metbins);
  TH1F* met_ele_NoJECRebin = (TH1F*)met_ele_NoJEC->Rebin(nMetbins,"met_ele_NoJECRebin",Metbins);
  TH1F* met_ele_JECupRebin = (TH1F*)met_ele_JECup->Rebin(nMetbins,"met_ele_JECupRebin",Metbins);
  TH1F* met_ele_JECdownRebin = (TH1F*)met_ele_JECdown->Rebin(nMetbins,"met_ele_JECdownRebin",Metbins);

  cout<<"Events Rebin:  raw: "<<met_ele_Rebin->Integral(0,-1)<<"  NoJEC: "<<met_ele_NoJECRebin->Integral(0,-1)<<"  JECup: "<<met_ele_JECupRebin->Integral(0,-1)<<"  JECdown: "<<met_ele_JECdownRebin->Integral(0,-1)<<endl;


  AddOverflowToLastBin(met_ele_Rebin);
  AddOverflowToLastBin(met_ele_NoJECRebin);
  AddOverflowToLastBin(met_ele_JECupRebin);
  AddOverflowToLastBin(met_ele_JECdownRebin);

  met_ele_Rebin->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]"); met_ele_Rebin->GetXaxis()->SetTitleOffset(1.);
  //met_ele_Rebin->GetYaxis()->SetRangeUser(300,2380);
  //met_ele_Rebin->GetXaxis()->SetRangeUser(120,121);
  met_ele_Rebin->Draw("pe");
  //met_ele_NoJECRebin->Draw("pesames");
  met_ele_JECupRebin->Draw("pesames");
  met_ele_JECdownRebin->Draw("pesames");
  
  c1->Print("Plots/Higgs/JEC_up_down_met_ele.png");
  c1->Print("Plots/Higgs/JEC_up_down_met_ele.pdf");

  for(int i=1;i<=met_ele_Rebin->GetNbinsX();i++){
    if(i==1)cout<<"ele met"<<endl<<"No JEC  "<<endl;
    cout<<"bin:"<<i<<" events: "<<met_ele_Rebin->GetBinContent(i)<<endl;
  }
  for(int i=1;i<=met_ele_Rebin->GetNbinsX();i++){
    if(i==1)cout<<"JECup   "<<endl;
    cout<<"bin:"<<i<<" events: "<<met_ele_JECupRebin->GetBinContent(i)<<endl;
  }
  for(int i=1;i<=met_ele_Rebin->GetNbinsX();i++){
    if(i==1)cout<<"JEC down"<<endl;
    cout<<"bin:"<<i<<" events: "<<met_ele_JECdownRebin->GetBinContent(i)<<endl;
  }

  cout<<"Events Rebin after Overflow add:  raw: "<<met_ele_Rebin->Integral(0,-1)<<"  NoJEC: "<<met_ele_NoJECRebin->Integral(0,-1)<<"  JECup: "<<met_ele_JECupRebin->Integral(0,-1)<<"  JECdown: "<<met_ele_JECdownRebin->Integral(0,-1)<<endl;

  //now muon

  TH2F* met_mu_Vinvmass = (TH2F*)f.Get("ggMetVsInvarMass_Loose_1Mu_0_1Jets");
  TH1F* met_mu= (TH1F*)met_mu_Vinvmass->ProjectionY("met_mu");
  //met_mu->Rebin(50);
  met_mu->SetTitle("");
  met_mu->SetMarkerSize(0.5);

  TH1F* met_mu_NoJEC = f.Get("ggMet_NoJEC_Loose_1Mu_0_1Jets");
  //met_mu_NoJEC->Rebin(50);
  met_mu_NoJEC->SetLineColor(kGreen);

  TH1F* met_mu_JECup = f.Get("ggMet_JECup_Loose_1Mu_0_1Jets");
  //met_mu_JECup->Rebin(50);
  met_mu_JECup->SetLineColor(kRed);
  met_mu_JECup->SetMarkerColor(kRed);
  met_mu_JECup->SetMarkerSize(0.5);

  TH1F* met_mu_JECdown = f.Get("ggMet_JECdown_Loose_1Mu_0_1Jets");
  //met_mu_JECdown->Rebin(50);
  met_mu_JECdown->SetLineColor(kBlue);
  met_mu_JECdown->SetMarkerColor(kBlue);
  met_mu_JECdown->SetMarkerSize(0.5);
  cout<<"Events:  raw: "<<met_mu->Integral(0,-1)<<"  NoJEC: "<<met_mu_NoJEC->Integral(0,-1)<<"  JECup: "<<met_mu_JECup->Integral(0,-1)<<"  JECdown: "<<met_mu_JECdown->Integral(0,-1)<<endl;

  TH1F* met_mu_Rebin = (TH1F*)met_mu->Rebin(nMetbins,"met_mu_Rebin",Metbins);
  TH1F* met_mu_NoJECRebin = (TH1F*)met_mu_NoJEC->Rebin(nMetbins,"met_mu_NoJECRebin",Metbins);
  TH1F* met_mu_JECupRebin = (TH1F*)met_mu_JECup->Rebin(nMetbins,"met_mu_JECupRebin",Metbins);
  TH1F* met_mu_JECdownRebin = (TH1F*)met_mu_JECdown->Rebin(nMetbins,"met_mu_JECdownRebin",Metbins);

  cout<<"Events Rebin:  raw: "<<met_mu_Rebin->Integral(0,-1)<<"  NoJEC: "<<met_mu_NoJECRebin->Integral(0,-1)<<"  JECup: "<<met_mu_JECupRebin->Integral(0,-1)<<"  JECdown: "<<met_mu_JECdownRebin->Integral(0,-1)<<endl;


  AddOverflowToLastBin(met_mu_Rebin);
  AddOverflowToLastBin(met_mu_NoJECRebin);
  AddOverflowToLastBin(met_mu_JECupRebin);
  AddOverflowToLastBin(met_mu_JECdownRebin);

  met_mu_Rebin->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]"); met_mu_Rebin->GetXaxis()->SetTitleOffset(1.);
  //met_mu_Rebin->GetYaxis()->SetRangeUser(450,590);
  //met_mu_Rebin->GetXaxis()->SetRangeUser(0,80);
  met_mu_Rebin->Draw("pe");
  //met_mu_NoJECRebin->Draw("pesames");
  met_mu_JECupRebin->Draw("pesames");
  met_mu_JECdownRebin->Draw("pesames");
  
  c1->Print("Plots/Higgs/JEC_up_down_met_mu.png");
  c1->Print("Plots/Higgs/JEC_up_down_met_mu.pdf");

  for(int i=1;i<=met_mu_Rebin->GetNbinsX();i++){
    if(i==1)cout<<"mu met"<<endl<<"No JEC  "<<endl;
    cout<<"bin:"<<i<<" events: "<<met_mu_Rebin->GetBinContent(i)<<endl;
  }
  for(int i=1;i<=met_mu_Rebin->GetNbinsX();i++){
    if(i==1)cout<<"JECup   "<<endl;
    cout<<"bin:"<<i<<" events: "<<met_mu_JECupRebin->GetBinContent(i)<<endl;
  }
  for(int i=1;i<=met_mu_Rebin->GetNbinsX();i++){
    if(i==1)cout<<"JEC down"<<endl;
    cout<<"bin:"<<i<<" events: "<<met_mu_JECdownRebin->GetBinContent(i)<<endl;
  }

  cout<<"Events Rebin after Overflow add:  raw: "<<met_mu_Rebin->Integral(0,-1)<<"  NoJEC: "<<met_mu_NoJECRebin->Integral(0,-1)<<"  JECup: "<<met_mu_JECupRebin->Integral(0,-1)<<"  JECdown: "<<met_mu_JECdownRebin->Integral(0,-1)<<endl;

  return;
  
}

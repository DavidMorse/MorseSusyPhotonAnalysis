void MHT(){

  gStyle->SetOptStat(0);

  TFile f("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Photon_handMade_QCD_Pt-350_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug26.root","READ");

  f.cd();

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  c1->SetLogy(1);

  Double_t xbins[]={0,25,50,100,150,200,300,400,500,600,800,1000};///ICHEP binning
  int nBins=(sizeof(xbins)/sizeof(Double_t))-1;
  float PlotXmax = xbins[nBins];

  TH1F* ggMHTana = (TH1F*)f.Get("ggMHT");ggMHTana->SetTitle("");
  TH1F* ggMHTana_reweight = (TH1F*)ggMHTana->Clone();
  TH1F* ffMHTana = (TH1F*)f.Get("ffMHT");
  TH1F* gfMHTana = (TH1F*)f.Get("gammafakeMHT");

  float x = ggMHTana->Integral()/ffMHTana->Integral();
  ffMHTana->Scale(x);
  x = ggMHTana->Integral()/gfMHTana->Integral();
  gfMHTana->Scale(x);

  TH1F* gg = (TH1F*)ggMHTana->Rebin(nBins,"gg",xbins);
  TH1F* ff = (TH1F*)ffMHTana->Rebin(nBins,"ff",xbins);
  TH1F* gf = (TH1F*)gfMHTana->Rebin(nBins,"gf",xbins);

  for(int i=0;i<gg->GetNbinsX()+1;i++){
    float val = gg->GetBinContent(i);
    float err = gg->GetBinError(i);
    gg->SetBinContent(i,val/gg->GetBinWidth(i));gg->SetBinError(i,err/gg->GetBinWidth(i));
    val = ff->GetBinContent(i);
    err = ff->GetBinError(i);
    ff->SetBinContent(i,val/gg->GetBinWidth(i));ff->SetBinError(i,err/gg->GetBinWidth(i));
    val = gf->GetBinContent(i);
    err = gf->GetBinError(i);
    gf->SetBinContent(i,val/gg->GetBinWidth(i));gf->SetBinError(i,err/gg->GetBinWidth(i));
  }

  ff->SetLineColor(kRed);ff->SetMarkerColor(kRed);
  gf->SetLineColor(kBlue);gf->SetMarkerColor(kBlue);

  gg->GetYaxis()->SetRangeUser(1e-2,80);
  gg->Draw("pe");
  ff->Draw("peSAMES");
  gf->Draw("peSAMES");

  c1->Print("Plots/Closure/handmade/moreEvents/MHT.png");

  TH1F* ffMHTana_reweight = (TH1F*)f.Get("ffMHT_reweight");
  TH1F* gfMHTana_reweight = (TH1F*)f.Get("gammafakeMHT_reweight");

  float x = ggMHTana_reweight->Integral()/ffMHTana_reweight->Integral();
  ffMHTana_reweight->Scale(x);
  x = ggMHTana_reweight->Integral()/gfMHTana_reweight->Integral();
  gfMHTana_reweight->Scale(x);

  TH1F* gg_reweight = (TH1F*)ggMHTana_reweight->Rebin(nBins,"gg_reweight",xbins);
  TH1F* ff_reweight = (TH1F*)ffMHTana_reweight->Rebin(nBins,"ff_reweight",xbins);
  TH1F* gf_reweight = (TH1F*)gfMHTana_reweight->Rebin(nBins,"gf_reweight",xbins);

  for(int i=0;i<gg_reweight->GetNbinsX()+1;i++){
    float val = gg_reweight->GetBinContent(i);
    float err = gg_reweight->GetBinError(i);
    gg_reweight->SetBinContent(i,val/gg_reweight->GetBinWidth(i));gg_reweight->SetBinError(i,err/gg_reweight->GetBinWidth(i));
    val = ff_reweight->GetBinContent(i);
    err = ff_reweight->GetBinError(i);
    ff_reweight->SetBinContent(i,val/gg_reweight->GetBinWidth(i));ff_reweight->SetBinError(i,err/gg_reweight->GetBinWidth(i));
    val = gf_reweight->GetBinContent(i);
    err = gf_reweight->GetBinError(i);
    gf_reweight->SetBinContent(i,val/gg_reweight->GetBinWidth(i));gf_reweight->SetBinError(i,err/gg_reweight->GetBinWidth(i));
  }

  ff_reweight->SetLineColor(kRed);ff_reweight->SetMarkerColor(kRed);
  gf_reweight->SetLineColor(kBlue);gf_reweight->SetMarkerColor(kBlue);

  gg_reweight->GetYaxis()->SetRangeUser(1e-2,80);
  gg_reweight->Draw("pe");
  ff_reweight->Draw("peSAMES");
  gf_reweight->Draw("peSAMES");

  c1->Print("Plots/Closure/handmade/moreEvents/MHT_reweight.png");

  TH1F* ggMetana_reweight = (TH1F*)f.Get("ggMet");
  TH1F* ffMetana_reweight = (TH1F*)f.Get("ffMet_reweightJet_binned");
  TH1F* gfMetana_reweight = (TH1F*)f.Get("gammafakeMet_reweightJet_binned");

  float x = ggMetana_reweight->Integral()/ffMetana_reweight->Integral();
  ffMetana_reweight->Scale(x);
  x = ggMetana_reweight->Integral()/gfMetana_reweight->Integral();
  gfMetana_reweight->Scale(x);

  TH1F* gg_reweight = (TH1F*)ggMetana_reweight->Rebin(nBins,"gg_reweight",xbins);
  TH1F* ff_reweight = (TH1F*)ffMetana_reweight->Rebin(nBins,"ff_reweight",xbins);
  TH1F* gf_reweight = (TH1F*)gfMetana_reweight->Rebin(nBins,"gf_reweight",xbins);

  for(int i=0;i<gg_reweight->GetNbinsX()+1;i++){
    float val = gg_reweight->GetBinContent(i);
    float err = gg_reweight->GetBinError(i);
    gg_reweight->SetBinContent(i,val/gg_reweight->GetBinWidth(i));gg_reweight->SetBinError(i,err/gg_reweight->GetBinWidth(i));
    val = ff_reweight->GetBinContent(i);
    err = ff_reweight->GetBinError(i);
    ff_reweight->SetBinContent(i,val/gg_reweight->GetBinWidth(i));ff_reweight->SetBinError(i,err/gg_reweight->GetBinWidth(i));
    val = gf_reweight->GetBinContent(i);
    err = gf_reweight->GetBinError(i);
    gf_reweight->SetBinContent(i,val/gg_reweight->GetBinWidth(i));gf_reweight->SetBinError(i,err/gg_reweight->GetBinWidth(i));
  }

  ff_reweight->SetLineColor(kRed);ff_reweight->SetMarkerColor(kRed);
  gf_reweight->SetLineColor(kBlue);gf_reweight->SetMarkerColor(kBlue);

  gg_reweight->GetYaxis()->SetRangeUser(1e-4,400);
  gg_reweight->Draw("pe");
  ff_reweight->Draw("peSAMES");
  gf_reweight->Draw("peSAMES");

  c1->Print("Plots/Closure/handmade/moreEvents/MHT_MET_reweight.png");




}

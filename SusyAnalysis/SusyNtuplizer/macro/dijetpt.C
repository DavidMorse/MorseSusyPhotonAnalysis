{

  TFile f("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_15GevCHfakeCut_PixelVetoOnFakes_NewDiEMPtBins_19499pb_NewJetMatching_2JetReq_May7.root","READ");

  Double_t bins[19]={0,50,55,60,65,70,75,80,85,90,95,100,110,130,150,175,200,300,800};int nbins=(sizeof(bins)/sizeof(Double_t))-1;


  TH1F* gg = (((TH1F*)f.Get("ggDiJetPt_scalar"))->Rebin(nbins,"gg",bins));
  TH1F* ee = (((TH1F*)f.Get("eeDiJetPt_scalar"))->Rebin(nbins,"ee",bins));
  TH1F* ff = (((TH1F*)f.Get("ffDiJetPt_scalar"))->Rebin(nbins,"ff",bins));
  TH1F* gf = (((TH1F*)f.Get("gammafakeDiJetPt_scalar"))->Rebin(nbins,"gf",bins));
  TH1F* gg0j = (((TH1F*)f.Get("ggDiJetPt_scalar_0Jet"))->Rebin(nbins,"gg0j",bins));
  TH1F* gg1j = (((TH1F*)f.Get("ggDiJetPt_scalar_1Jet"))->Rebin(nbins,"gg1j",bins));
  TH1F* gg2j = (((TH1F*)f.Get("ggDiJetPt_scalar_2Jet"))->Rebin(nbins,"gg2j",bins));
  TH1F* ee0j = (((TH1F*)f.Get("eeDiJetPt_scalar_0Jet"))->Rebin(nbins,"ee0j",bins));
  TH1F* ee1j = (((TH1F*)f.Get("eeDiJetPt_scalar_1Jet"))->Rebin(nbins,"ee1j",bins));
  TH1F* ee2j = (((TH1F*)f.Get("eeDiJetPt_scalar_2Jet"))->Rebin(nbins,"ee2j",bins));
  TH1F* ff0j = (((TH1F*)f.Get("ffDiJetPt_scalar_0Jet"))->Rebin(nbins,"ff0j",bins));
  TH1F* ff1j = (((TH1F*)f.Get("ffDiJetPt_scalar_1Jet"))->Rebin(nbins,"ff1j",bins));
  TH1F* ff2j = (((TH1F*)f.Get("ffDiJetPt_scalar_2Jet"))->Rebin(nbins,"ff2j",bins));
  TH1F* gf0j = (((TH1F*)f.Get("gammafakeDiJetPt_scalar_0Jet"))->Rebin(nbins,"gf0j",bins));
  TH1F* gf1j = (((TH1F*)f.Get("gammafakeDiJetPt_scalar_1Jet"))->Rebin(nbins,"gf1j",bins));
  TH1F* gf2j = (((TH1F*)f.Get("gammafakeDiJetPt_scalar_2Jet"))->Rebin(nbins,"gf2j",bins));

  /*gg->Rebin(2);ee->Rebin(2);ff->Rebin(2);gf->Rebin(2);
  gg0j->Rebin(2);ee0j->Rebin(2);ff0j->Rebin(2);gf0j->Rebin(2);
  gg1j->Rebin(2);ee1j->Rebin(2);ff1j->Rebin(2);gf1j->Rebin(2);
  gg2j->Rebin(2);ee2j->Rebin(2);ff2j->Rebin(2);gf2j->Rebin(2);*/

  TCanvas c1("c1","c1",800,600);c1.cd();

  float ggint = gg->Integral(),eeint = ee->Integral(),ffint = ff->Integral(),gfint = gf->Integral();

  ee->Scale(ggint/eeint);ff->Scale(ggint/ffint);gf->Scale(ggint/gfint);
  ee0j->Scale(ggint/eeint);ff0j->Scale(ggint/ffint);gf0j->Scale(ggint/gfint);
  ee1j->Scale(ggint/eeint);ff1j->Scale(ggint/ffint);gf1j->Scale(ggint/gfint);
  ee2j->Scale(ggint/eeint);ff2j->Scale(ggint/ffint);gf2j->Scale(ggint/gfint);

  TH1F* ggee = (TH1F*)gg->Clone();ggee->Divide(ee);
  TH1F* ggff = (TH1F*)gg->Clone();ggff->Divide(ee);
  TH1F* gggf = (TH1F*)gg->Clone();gggf->Divide(ee);

  TH1F* ggee0j = (TH1F*)gg0j->Clone();ggee0j->Divide(ee0j);
  TH1F* ggff0j = (TH1F*)gg0j->Clone();ggff0j->Divide(ff0j);
  TH1F* gggf0j = (TH1F*)gg0j->Clone();gggf0j->Divide(gf0j);

  TH1F* ggee1j = (TH1F*)gg1j->Clone();ggee1j->Divide(ee1j);
  TH1F* ggff1j = (TH1F*)gg1j->Clone();ggff1j->Divide(ff1j);
  TH1F* gggf1j = (TH1F*)gg1j->Clone();gggf1j->Divide(gf1j);

  TH1F* ggee2j = (TH1F*)gg2j->Clone();ggee2j->Divide(ee2j);
  TH1F* ggff2j = (TH1F*)gg2j->Clone();ggff2j->Divide(ff2j);
  TH1F* gggf2j = (TH1F*)gg2j->Clone();gggf2j->Divide(gf2j);

  ee0j->SetLineColor(kRed);ee0j->SetMarkerColor(kRed);ee0j->Draw();
  ff0j->SetLineColor(kBlue);ff0j->SetMarkerColor(kBlue);ff0j->Draw("SAMES");
  gg0j->Draw("SAMES");

  //ggee->Draw();
  ggee2j->Draw();ggff0j->SetLineColor(kBlue);//ggff0j->SetMarkerColor(kBlue);ggff0j->Draw("sames");


}

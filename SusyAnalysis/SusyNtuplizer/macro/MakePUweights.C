void MakePUweights(){

  TFile f_data("pileupData.root","READ");
  //TFile f_ggHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Jun14.root","READ");
  //TFile f_WZHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Jun14.root","READ");
  //TFile f_TTHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_TTH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Jun14.root","READ");
  //TFile f_VBFHgg("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Higgs_cms533v1_VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Jun14.root","READ");
  //TFile f_DiPhoBox("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Photon_handMade_DiPho_Box_Pt20_doubleEMEnriched_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July10_2.root","READ");  
  //TFile f_DiPhoBorn("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Photon_handMade_DiPho_Born_Pt20_doubleEMEnriched_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July10.root","READ");
  //TFile f_QCD30_40("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Photon_handMade_QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July12.root","READ");  
  //TFile f_QCD40_inf("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Photon_handMade_QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July12.root","READ");  
  //TFile f_GJet20_40("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Photon_handMade_GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July10.root","READ");  
  //TFile f_GJet40_inf("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Photon_handMade_GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July8.root","READ"); 
  //TFile f_aaw("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Photon_cms533_WinoNLSP_chargino130_bino1_5_10_15_hw_aaw_full_July4.root","READ");  
  TFile f_SMS_WH("hist_Photon_SMS_TChiWH_WincHgg_2J.root");
  //TFile f_SMS_ZH("hist_Photon_SMS_TChiZH_ZincHgg_2J.root");
  //TFile f_SMS_HH_2W2g("hist_Photon_SMS_TChiHH_2W2g_2J.root");
  //TFile f_SMS_HH_2tau2g("hist_Photon_SMS_TChiHH_2tau2g_2J.root");
  //TFile f_SMS_HH_2Z2g("hist_Photon_SMS_TChiHH_2Z2g_2J.root");
  //TFile f_SMS_HH_2b2g("hist_Photon_SMS_TChiHH_2b2g_2J.root");
  //TFile f_ZG("hist_Photon_ZGToLLG.root");


  //TH1F* h_numTrueInt = (TH1F*)f_ggHgg.Get("numTrueInt");h_numTrueInt->Sumw2();
  //TH1F* h_numTrueInt = (TH1F*)f_WZHgg.Get("numTrueInt");h_numTrueInt->Sumw2();
  //TH1F* h_numTrueInt = (TH1F*)f_TTHgg.Get("numTrueInt");h_numTrueInt->Sumw2();
  //TH1F* h_numTrueInt = (TH1F*)f_VBFHgg.Get("numTrueInt");h_numTrueInt->Sumw2();
  //TH1F* h_numTrueInt = (TH1F*)f_DiPhoBorn.Get("numTrueInt");h_numTrueInt->Sumw2();
  //TH1F* h_numTrueInt = (TH1F*)f_DiPhoBox.Get("numTrueInt");h_numTrueInt->Sumw2();
  //TH1F* h_numTrueInt = (TH1F*)f_QCD30_40.Get("numTrueInt");h_numTrueInt->Sumw2();
  //TH1F* h_numTrueInt = (TH1F*)f_QCD40_inf.Get("numTrueInt");h_numTrueInt->Sumw2();
  //TH1F* h_numTrueInt = (TH1F*)f_GJet20_40.Get("numTrueInt");h_numTrueInt->Sumw2();
  //TH1F* h_numTrueInt = (TH1F*)f_GJet40_inf.Get("numTrueInt");h_numTrueInt->Sumw2();
  //TH1F* h_numTrueInt = (TH1F*)f_aaw.Get("numTrueInt");h_numTrueInt->Sumw2();
  TH1F* h_numTrueInt = (TH1F*)f_SMS_WH.Get("numTrueInt");h_numTrueInt->Sumw2();
  //TH1F* h_numTrueInt = (TH1F*)f_SMS_ZH.Get("numTrueInt");h_numTrueInt->Sumw2();
  //TH1F* h_numTrueInt = (TH1F*)f_SMS_HH_2W2g.Get("numTrueInt");h_numTrueInt->Sumw2();
  //TH1F* h_numTrueInt = (TH1F*)f_SMS_HH_2Z2g.Get("numTrueInt");h_numTrueInt->Sumw2();
  //TH1F* h_numTrueInt = (TH1F*)f_SMS_HH_2b2g.Get("numTrueInt");h_numTrueInt->Sumw2();
  //TH1F* h_numTrueInt = (TH1F*)f_SMS_HH_2tau2g.Get("numTrueInt");h_numTrueInt->Sumw2();
  //TH1F* h_numTrueInt = (TH1F*)f_ZG.Get("numTrueInt");h_numTrueInt->Sumw2();

  TH1F* h_dataPileup = (TH1F*)f_data.Get("pileup");h_dataPileup->Sumw2();
  
  float scale = 1./h_dataPileup->Integral(0,-1);h_dataPileup->Scale(scale);
  scale = 1./h_numTrueInt->Integral(0,-1);h_numTrueInt->Scale(scale);
  
  /*
  TFile f_ggOut("PUweightsGluGluH.root","RECREATE");f_ggOut.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_ggOut.Write();
  f_ggOut.Close();
  */
  /*
  TFile f_WZOut("PUweightsWZH.root","RECREATE");f_WZOut.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_WZOut.Write();
  f_WZOut.Close();
  */
  /*
  TFile f_TTOut("PUweightsTTH.root","RECREATE");f_TTOut.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_TTOut.Write();
  f_TTOut.Close();
  */
  /*
  TFile f_VBFOut("PUweightsVBFH.root","RECREATE");f_VBFOut.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_VBFOut.Write();
  f_VBFOut.Close();
  */
  /*
  TFile f_DiPhoBornOut("PUweightsDiPhoBorn.root","RECREATE");f_DiPhoBornOut.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_DiPhoBornOut.Write();
  f_DiPhoBornOut.Close();
  */
  /*
  TFile f_DiPhoBoxOut("PUweightsDiPhoBox.root","RECREATE");f_DiPhoBoxOut.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_DiPhoBoxOut.Write();
  f_DiPhoBoxOut.Close();
  */
  /*
  TFile f_QCD30_40Out("PUweightsQCD30_40.root","RECREATE");f_QCD30_40Out.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_QCD30_40Out.Write();
  f_QCD30_40Out.Close();
  */
  /*
  TFile f_QCD40_infOut("PUweightsQCD40_inf.root","RECREATE");f_QCD40_infOut.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_QCD40_infOut.Write();
  f_QCD40_infOut.Close();
  */
  /*
  TFile f_GJet20_40Out("PUweightsGJet20_40.root","RECREATE");f_GJet20_40Out.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_GJet20_40Out.Write();
  f_GJet20_40Out.Close();
  */
  /*
  TFile f_GJet40_infOut("PUweightsGJet40_inf.root","RECREATE");f_GJet40_infOut.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_GJet40_infOut.Write();
  f_GJet40_infOut.Close();
  */
  /*
  TFile f_aawOut("PUweightsAAW.root","RECREATE");f_aawOut.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_aawOut.Write();
  f_aawOut.Close();
  */
  
  TFile f_SMS_WHOut("PUweightsSMS_WH_test.root","RECREATE");f_SMS_WHOut.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_SMS_WHOut.Write();
  f_SMS_WHOut.Close();
  
  /*
  TFile f_SMS_ZHOut("PUweightsSMS_ZH.root","RECREATE");f_SMS_ZHOut.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_SMS_ZHOut.Write();
  f_SMS_ZHOut.Close();
  */
  /*
  TFile f_SMS_HH_2W2gOut("PUweightsSMS_HH_2W2g.root","RECREATE");f_SMS_HH_2W2gOut.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_SMS_HH_2W2gOut.Write();
  f_SMS_HH_2W2gOut.Close();
  */
  /*
  TFile f_SMS_HH_2Z2gOut("PUweightsSMS_HH_2Z2g.root","RECREATE");f_SMS_HH_2Z2gOut.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_SMS_HH_2Z2gOut.Write();
  f_SMS_HH_2Z2gOut.Close();
  */
  /*
  TFile f_SMS_HH_2b2gOut("PUweightsSMS_HH_2b2g.root","RECREATE");f_SMS_HH_2b2gOut.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_SMS_HH_2b2gOut.Write();
  f_SMS_HH_2b2gOut.Close();
  */
  /*
  TFile f_SMS_HH_2tau2gOut("PUweightsSMS_HH_2tau2g.root","RECREATE");f_SMS_HH_2tau2gOut.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_SMS_HH_2tau2gOut.Write();
  f_SMS_HH_2tau2gOut.Close();
  */
  /*
  TFile f_ZGOut("PUweightsZG.root","RECREATE");f_ZGOut.cd();
  TH1F* PUweights = (TH1F*)h_dataPileup->Clone("PUweights");
  PUweights->Divide(h_numTrueInt);
  f_ZGOut.Write();
  f_ZGOut.Close();
  */
}

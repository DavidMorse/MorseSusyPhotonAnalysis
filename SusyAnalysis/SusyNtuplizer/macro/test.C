void test(){


  TFile *eos = TFile::Open("root://eoscms//eos/cms/store/user/dmorse/AnalysisOutput/cms533v1/newLepDZ/hist_HiggsAna_Photon_SMS_TChiHH_2W2g_2J.root","READ"); 
  TFile *hnew = TFile::Open("hist_HiggsAna_Photon_SMS_TChiHH_2W2g_2J.root","READ");

  TH3F* muOld = (TH3F*)eos->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT");
  TH3F* muOldNoSF = (TH3F*)eos->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor");
  TH1F* muOldProj = (TH1F*)muOld->ProjectionZ("muOldProj",1,1,1,1,"o");
  TH1F* muOldNoSFProj = (TH1F*)muOldNoSF->ProjectionZ("muOldNoSFProj",1,1,1,1,"o");

  cout<<"muOld: "<<muOld->Integral()<<"  muOld: "<<muOld->GetEntries()<<endl;
  cout<<"muOldNoSF: "<<muOldNoSF->Integral()<<"  muOldNoSF: "<<muOldNoSF->GetEntries()<<endl;
  cout<<"muOld: "<<muOld->Integral(1,1,1,1,0,-1)<<"  muOldProj: "<<muOldProj->Integral()<<endl;
  cout<<"muOldProj: "<<muOldProj->Integral(0,-1)<<"  muOldProj: "<<muOldProj->GetEntries()<<endl;
  cout<<"muOldNoSFProj: "<<muOldNoSFProj->Integral(0,-1)<<"  muOldNoSFProj: "<<muOldNoSFProj->GetEntries()<<endl;

  TH3F* muNew = (TH3F*)hnew->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT");
  TH1F* muNewProj = (TH1F*)muNew->ProjectionZ("h",1,1,1,1,"o");
  TH3F* muNewNoSF = (TH3F*)hnew->Get("gg_SMS_Loose_1Ele_0_1Jets_mChi_mBino_MT_noScaleFactor");
  TH1F* muNewNoSFProj = (TH1F*)muNewNoSF->ProjectionZ("muNewNoSFProj",1,1,1,1,"o");

  cout<<"muNew: "<<muNew->Integral()<<"  muNew: "<<muNew->GetEntries()<<endl;
  cout<<"muNew: "<<muNew->Integral(1,1,1,1,0,-1)<<"  muNewProj: "<<muNewProj->Integral()<<endl;
  cout<<"muNewNoSFProj: "<<muNewNoSFProj->Integral(0,-1)<<"  muNewNoSFProj: "<<muNewNoSFProj->GetEntries()<<endl;

  TH1F* muPt = (TH1F*)hnew->Get("ggEle_ElePt");
  cout<<"muPt: "<<muPt->Integral()<<"  "<<muPt->Integral(0,-1);
}

void met_Phi(){

  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","",800,600);

  TFile F("hist_HiggsAna_Data2012A_13July_06Aug_recover_Data2012B_13July_Data2012C_Prompt_Runs198022-198903_Runs198941-203742_Data2012D_Prompt_Runs207920-209151_Runs203777-207905_Filter_Ana_HLT_JSON_Two40-25GeVbarrelPhotons_ALL-Higgsino-selectedEvents_temp.root","READ");

  TH1F* gg = (TH1F*)F.Get("met_Phi_ggLoose");
  TH1F* ggC = (TH1F*)gg->Clone("ggC");
  TH1F* ggEle = (TH1F*)F.Get("met_Phi_ggLoose_Ele");
  TH1F* ggEleC = (TH1F*)ggEle->Clone("ggEleC");
  TH1F* ggMu = (TH1F*)F.Get("met_Phi_ggLoose_Mu");
  TH1F* ggMuC = (TH1F*)ggMu->Clone("ggMuC");

  gg->Fit("pol5","","",0,6.4);
  ggC->Fit("pol0","","",0,6.4);

  gg->Draw("PE");ggC->Draw("PESAMES");
  c1->Print("Plots/Higgs/met_Phi_ggLoose_all_Fit.png");
  c1->Print("Plots/Higgs/met_Phi_ggLoose_all_Fit.pdf");

}

#include "TPaveStats.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH1F.h"

void pileup(){
  TFile f("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/Pileup/cms525v3_jec2012/hist_Data2012B_Filter_JsonNVertexOne25GeVBarrelPho_Photon_cms525v3_jec2012_Runs190456-195947_Pileup.root");
  f.cd();

  TCanvas *c1 = new TCanvas("c1","",800,600);

  gStyle->SetOptStat(0);
  gStyle->SetFitFormat("5.3f");

  TH1F *h_rho25 = (TH1F*)f.Get("rho25");
  h_rho25->GetXaxis()->SetRangeUser(0,50);
  h_rho25->Draw("histo");
  c1->Print("Plots/Pileup/rho25.png");
  c1->Print("Plots/Pileup/rho25.pdf");

  c1->SetLogz(1);

  TPaveText *Text;
  Text = new TPaveText(.25,.65,.55,.79,"NDC");
  Text->AddText("CMS Preliminary 2012");
  Text->SetFillStyle(4000);
  Text->SetFillColor(0);
  Text->SetBorderSize(0);
  Text->AddText("#sqrt{s} = 8 TeV, #intL = 4.04 fb^{-1}");

  TH2F* h_EcalDR03 = (TH2F*)f.Get("EcalIsoVsRho25DR03_ee");
  h_EcalDR03->GetYaxis()->SetRangeUser(-25,70);
  h_EcalDR03->GetXaxis()->SetRangeUser(0,40);
  h_EcalDR03->SetTitle("");
  h_EcalDR03->GetYaxis()->SetTitle("Ecal Isolation (DR03)");
  h_EcalDR03->GetXaxis()->SetTitle("Rho25 (E_{T} / unit area)");
  h_EcalDR03->Draw("colz");
  Text->Draw();
  c1->Print("Plots/Pileup/Figure5_EcalIsoVsRho25DR03_ee.png");
  c1->Print("Plots/Pileup/Figure5_EcalIsoVsRho25DR03_ee.pdf");

  TH2F* h_HcalDR03 = (TH2F*)f.Get("HcalIsoVsRho25DR03_ee");
  h_HcalDR03->GetYaxis()->SetRangeUser(-3,45);
  h_HcalDR03->GetXaxis()->SetRangeUser(0,40);
  h_HcalDR03->SetTitle("");
  h_HcalDR03->GetYaxis()->SetTitle("Hcal Isolation (DR03)");
  h_HcalDR03->GetXaxis()->SetTitle("Rho25 (E_{T} / unit area)");
  h_HcalDR03->Draw("colz");
  Text->Draw();
  c1->Print("Plots/Pileup/Figure5_HcalIsoVsRho25DR03_ee.png");
  c1->Print("Plots/Pileup/Figure5_HcalIsoVsRho25DR03_ee.pdf");

  TH2F* h_TrackDR03 = (TH2F*)f.Get("TrackIsoVsRho25DR03_ee");
  h_TrackDR03->GetYaxis()->SetRangeUser(-3,80);
  h_TrackDR03->GetXaxis()->SetRangeUser(0,40);
  h_TrackDR03->SetTitle("");
  h_TrackDR03->GetYaxis()->SetTitle("Track Isolation (DR03)");
  h_TrackDR03->GetXaxis()->SetTitle("Rho25 (E_{T} / unit area)");
  h_TrackDR03->Draw("colz");
  Text->Draw();
  c1->Print("Plots/Pileup/Figure5_TrackIsoVsRho25DR03_ee.png");
  c1->Print("Plots/Pileup/Figure5_TrackIsoVsRho25DR03_ee.pdf");

  TH2F* h_HOverE = (TH2F*)f.Get("HOverEVsRho25_ee");
  h_HOverE->GetYaxis()->SetRangeUser(0,.16);
  h_HOverE->GetXaxis()->SetRangeUser(0,40);
  h_HOverE->SetTitle("");
  h_HOverE->GetXaxis()->SetTitle("Rho25 (E_{T} / unit area)");
  h_HOverE->Draw("colz");
  Text->Draw();
  c1->Print("Plots/Pileup/Figure5_HOverEVsRho25_ee.png");
  c1->Print("Plots/Pileup/Figure5_HOverEVsRho25_ee.pdf");


  TH2F* h_HOverEAfterIso = (TH2F*)f.Get("HOverEAfterIsoVsRho25_ee");
  h_HOverEAfterIso->GetYaxis()->SetRangeUser(0,.16);
  h_HOverEAfterIso->GetXaxis()->SetRangeUser(0,40);
  h_HOverEAfterIso->SetTitle("");
  h_HOverEAfterIso->GetXaxis()->SetTitle("Rho25 (E_{T} / unit area)");
  h_HOverEAfterIso->Draw("colz");
  Text->Draw();
  c1->Print("Plots/Pileup/Figure5_HOverEAfterIsoVsRho25_ee.png");
  c1->Print("Plots/Pileup/Figure5_HOverEAfterIsoVsRho25_ee.pdf");

  h_EcalDR03->ProfileX("ecal_pfx",0,450);
  ecal_pfx->GetYaxis()->SetTitle("Average Ecal Isolation (DR03)");
  ecal_pfx->GetYaxis()->SetRangeUser(0,4.6);
  ecal_pfx->Fit("pol1","","",0,45);
  ecal_pfx->Draw();
  c1->Update();
  TPaveStats *st = (TPaveStats*)ecal_pfx->FindObject("stats");
  st->SetName("ecalStats");
  st->SetTextSize(.04);
  st->SetX1NDC(.25);st->SetX2NDC(.55);
  st->SetY1NDC(.52);st->SetY2NDC(.65);
  Text->Draw();
  c1->Print("Plots/Pileup/Figure6_EcalIsoVsRho25_ee_pfx.png");
  c1->Print("Plots/Pileup/Figure6_EcalIsoVsRho25_ee_pfx.pdf");

  h_HcalDR03->ProfileX("hcal_pfx",0,450);
  hcal_pfx->GetYaxis()->SetTitle("Average Hcal Isolation (DR03)");
  hcal_pfx->GetYaxis()->SetRangeUser(0.1,1.5);
  hcal_pfx->Fit("pol1","","",0,45);
  hcal_pfx->Draw();
  c1->Update();
  TPaveStats *st = (TPaveStats*)hcal_pfx->FindObject("stats");
  st->SetName("hcalStats");
  st->SetTextSize(.04);
  st->SetX1NDC(.25);st->SetX2NDC(.55);
  st->SetY1NDC(.52);st->SetY2NDC(.65);
  Text->Draw();
  c1->Print("Plots/Pileup/Figure6_HcalIsoVsRho25_ee_pfx.png");
  c1->Print("Plots/Pileup/Figure6_HcalIsoVsRho25_ee_pfx.pdf");

  h_TrackDR03->ProfileX("track_pfx",0,450);
  track_pfx->GetYaxis()->SetTitle("Average Track Isolation (DR03)");
  track_pfx->GetYaxis()->SetRangeUser(0.1,1.2);
  track_pfx->Fit("pol1","","",0,30);
  track_pfx->Draw();
  c1->Update();
  TPaveStats *st = (TPaveStats*)track_pfx->FindObject("stats");
  st->SetName("trackStats");
  st->SetTextSize(.04);
  st->SetX1NDC(.25);st->SetX2NDC(.55);
  st->SetY1NDC(.52);st->SetY2NDC(.65);
  Text->Draw();
  c1->Print("Plots/Pileup/Figure6_TrackIsoVsRho25_ee_pfx.png");
  c1->Print("Plots/Pileup/Figure6_TrackIsoVsRho25_ee_pfx.pdf");

  gStyle->SetFitFormat("5.6f");

  h_HOverE->ProfileX("hovere_pfx",0,450);
  hovere_pfx->GetYaxis()->SetTitle("Average H/E");
  hovere_pfx->Fit("pol1","","",0,45);
  hovere_pfx->Draw();
  Text->Draw();
  c1->Print("Plots/Pileup/Figure6_HOverEIsoVsRho25_ee_pfx.png");
  c1->Print("Plots/Pileup/Figure6_HOverEIsoVsRho25_ee_pfx.pdf");

  h_HOverEAfterIso->ProfileX("hovereafteriso_pfx",0,450);
  hovereafteriso_pfx->GetYaxis()->SetTitle("Average H/E After Combined Isolation (DR03) < 6.0");
  hovereafteriso_pfx->Fit("pol1","","",0,45);
  hovereafteriso_pfx->Draw();
  Text->Draw();
  c1->Print("Plots/Pileup/Figure6_HOverEAfterIsoVsRho25_ee_pfx.png");
  c1->Print("Plots/Pileup/Figure6_HOverEAfterIsoVsRho25_ee_pfx.pdf");

  TH2F* h_chargedHadron = (TH2F*)f.Get("chargedHadronIsoVsRho25_ee");
  h_chargedHadron->GetYaxis()->SetRangeUser(-3,120);
  h_chargedHadron->GetXaxis()->SetRangeUser(0,40);
  h_chargedHadron->SetTitle("");
  h_chargedHadron->Draw("colz");
  Text->Draw();
  c1->Print("Plots/Pileup/chargedHadronIso_ee.png");
  h_chargedHadron->ProfileX("chargedHad_pfx",0,450);
  chargedHad_pfx->GetYaxis()->SetTitle("Average chargedHadron Isolation");
  chargedHad_pfx->GetYaxis()->SetRangeUser(0.4,2);
  chargedHad_pfx->Fit("pol1","","",0,50);
  chargedHad_pfx->Draw();
  c1->Update();
  TPaveStats *st = (TPaveStats*)chargedHad_pfx->FindObject("stats");
  st->SetName("hcalStats");
  st->SetTextSize(.04);
  st->SetX1NDC(.21);st->SetX2NDC(.55);
  st->SetY1NDC(.52);st->SetY2NDC(.65);
  Text->Draw();
  c1->Print("Plots/Pileup/chargedHadronIso_ee_pfx.png");

  TH2F* h_neutralHadron = (TH2F*)f.Get("neutralHadronIsoVsRho25_ee");
  h_neutralHadron->GetYaxis()->SetRangeUser(-3,40);
  h_neutralHadron->GetXaxis()->SetRangeUser(0,40);
  h_neutralHadron->SetTitle("");
  h_neutralHadron->Draw("colz");
  Text->Draw();
  c1->Print("Plots/Pileup/neutralHadronIso_ee.png");
  h_neutralHadron->ProfileX("neutralHad_pfx",0,450);
  neutralHad_pfx->GetYaxis()->SetTitle("Average neutralHadron Isolation");
  neutralHad_pfx->GetYaxis()->SetRangeUser(0,1.2);
  neutralHad_pfx->Fit("pol1","","",0,50);
  neutralHad_pfx->Draw();
  c1->Update();
  TPaveStats *st = (TPaveStats*)neutralHad_pfx->FindObject("stats");
  st->SetName("hcalStats");
  st->SetTextSize(.04);
  st->SetX1NDC(.21);st->SetX2NDC(.55);
  st->SetY1NDC(.52);st->SetY2NDC(.65);
  Text->Draw();
  c1->Print("Plots/Pileup/neutralHadronIso_ee_pfx.png");

  TH2F* h_photon = (TH2F*)f.Get("photonIsoVsRho25_ee");
  h_photon->GetYaxis()->SetRangeUser(-3,100);
  h_photon->GetXaxis()->SetRangeUser(0,40);
  h_photon->SetTitle("");
  h_photon->Draw("colz");
  Text->Draw();
  c1->Print("Plots/Pileup/photonIso_ee.png");
  h_photon->ProfileX("photon_pfx",0,450);
  photon_pfx->GetYaxis()->SetTitle("Average photon Isolation");
  //photon_pfx->GetYaxis()->SetRangeUser(0,1.2);
  photon_pfx->Fit("pol1","","",0,50);
  photon_pfx->Draw();
  c1->Update();
  TPaveStats *st = (TPaveStats*)photon_pfx->FindObject("stats");
  st->SetName("hcalStats");
  st->SetTextSize(.04);
  st->SetX1NDC(.21);st->SetX2NDC(.55);
  st->SetY1NDC(.52);st->SetY2NDC(.65);
  Text->Draw();
  c1->Print("Plots/Pileup/photonIso_ee_pfx.png");

  f.Close();
}

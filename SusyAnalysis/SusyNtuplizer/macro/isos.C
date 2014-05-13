{

TFile f("Isos.root","RECREATE");
_file0.cd();
gStyle->SetOptStat(0);

_file0.cd();
float x = 1./EcalIsoDR03RhoCorr_gg->Integral();
EcalIsoDR03RhoCorr_gg->Scale(x);
EcalIsoDR03RhoCorr_gg->SetTitle("");
float x = 1./EcalIsoDR03RhoCorr_ee->Integral();
EcalIsoDR03RhoCorr_ee->Scale(x);
EcalIsoDR03RhoCorr_ee->SetLineColor(kRed);
EcalIsoDR03RhoCorr_ee->SetTitle("");
float x = 1./EcalIsoDR03RhoCorr_ff->Integral();
EcalIsoDR03RhoCorr_ff->Scale(x);
EcalIsoDR03RhoCorr_ff->SetLineColor(kBlue);
EcalIsoDR03RhoCorr_ff->SetTitle("");
EcalIsoDR03RhoCorr_ee->Draw();
EcalIsoDR03RhoCorr_gg->Draw("SAME");
EcalIsoDR03RhoCorr_ff->Draw("SAME");
f.cd();
c1->Write("EcalIsoRhoCorr");

_file0.cd();
float x = 1./HcalIsoDR03RhoCorr_gg->Integral();
HcalIsoDR03RhoCorr_gg->Scale(x);
HcalIsoDR03RhoCorr_gg->SetTitle("");
float x = 1./HcalIsoDR03RhoCorr_ee->Integral();
HcalIsoDR03RhoCorr_ee->Scale(x);
HcalIsoDR03RhoCorr_ee->SetLineColor(kRed);
HcalIsoDR03RhoCorr_ee->SetTitle("");
float x = 1./HcalIsoDR03RhoCorr_ff->Integral();
HcalIsoDR03RhoCorr_ff->Scale(x);
HcalIsoDR03RhoCorr_ff->SetLineColor(kBlue);
HcalIsoDR03RhoCorr_ff->SetTitle("");
HcalIsoDR03RhoCorr_ee->Draw();
HcalIsoDR03RhoCorr_gg->Draw("SAME");
HcalIsoDR03RhoCorr_ff->Draw("SAME");
f.cd();
c1->Write("HcalIsoRhoCorr");

_file0.cd();
float x = 1./TrackIsoDR03RhoCorr_gg->Integral();
TrackIsoDR03RhoCorr_gg->Scale(x);
TrackIsoDR03RhoCorr_gg->SetTitle("");
float x = 1./TrackIsoDR03RhoCorr_ee->Integral();
TrackIsoDR03RhoCorr_ee->Scale(x);
TrackIsoDR03RhoCorr_ee->SetLineColor(kRed);
TrackIsoDR03RhoCorr_ee->SetTitle("");
float x = 1./TrackIsoDR03RhoCorr_ff->Integral();
TrackIsoDR03RhoCorr_ff->Scale(x);
TrackIsoDR03RhoCorr_ff->SetLineColor(kBlue);
TrackIsoDR03RhoCorr_ff->SetTitle("");
TrackIsoDR03RhoCorr_ee->Draw();
TrackIsoDR03RhoCorr_gg->Draw("SAME");
TrackIsoDR03RhoCorr_ff->Draw("SAME");
f.cd();
c1->Write("TrackIsoRhoCorr");

_file0.cd();
float x = 1./CombIsoDR03RhoCorr_gg->Integral();
CombIsoDR03RhoCorr_gg->Scale(x);
CombIsoDR03RhoCorr_gg->SetTitle("");
float x = 1./CombIsoDR03RhoCorr_ee->Integral();
CombIsoDR03RhoCorr_ee->Scale(x);
CombIsoDR03RhoCorr_ee->SetLineColor(kRed);
CombIsoDR03RhoCorr_ee->SetTitle("");
float x = 1./CombIsoDR03RhoCorr_ff->Integral();
CombIsoDR03RhoCorr_ff->Scale(x);
CombIsoDR03RhoCorr_ff->SetLineColor(kBlue);
CombIsoDR03RhoCorr_ff->SetTitle("");
CombIsoDR03RhoCorr_ee->Draw();
CombIsoDR03RhoCorr_gg->Draw("SAME");
CombIsoDR03RhoCorr_ff->Draw("SAME");
f.cd();
c1->Write("CombIsoRhoCorr");

}

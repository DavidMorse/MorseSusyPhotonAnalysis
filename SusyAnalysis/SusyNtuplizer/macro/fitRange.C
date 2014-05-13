void fitRange(){

  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","",900,600);
  c1->cd();
 

  TH2F* h = new TH2F("h","",8,102.5,110.5,12,153.5,165.5);

  h->SetBinContent(1,1,9.67);h->SetBinError(1,1,1.9);
  h->SetBinContent(1,2,9.89);h->SetBinError(1,2,1.8);
  h->SetBinContent(1,3,9.71);h->SetBinError(1,3,1.75);
  h->SetBinContent(1,4,9.55);h->SetBinError(1,4,1.86);
  h->SetBinContent(1,5,9.40);h->SetBinError(1,5,1.64);
  h->SetBinContent(1,6,9.26);h->SetBinError(1,6,1.63);
  h->SetBinContent(1,7,9.26);h->SetBinError(1,7,1.63);
  h->SetBinContent(1,8,9.01);h->SetBinError(1,8,1.57);
  h->SetBinContent(1,9,8.91);h->SetBinError(1,9,1.6);
  h->SetBinContent(1,10,8.81);h->SetBinError(1,10,1.59);
  h->SetBinContent(1,11,8.72);h->SetBinError(1,11,1.52);
  h->SetBinContent(1,12,8.81);h->SetBinError(1,12,1.59);

  h->SetBinContent(2,1,9.57);h->SetBinError(2,1,1.86);
  h->SetBinContent(2,2,9.78);h->SetBinError(2,2,1.79);
  h->SetBinContent(2,3,9.60);h->SetBinError(2,3,1.75);
  h->SetBinContent(2,4,9.45);h->SetBinError(2,4,1.72);
  h->SetBinContent(2,5,9.30);h->SetBinError(2,5,1.69);
  h->SetBinContent(2,6,9.17);h->SetBinError(2,6,1.65);
  h->SetBinContent(2,7,9.05);h->SetBinError(2,7,1.63);
  h->SetBinContent(2,8,8.94);h->SetBinError(2,8,1.61);
  h->SetBinContent(2,9,8.84);h->SetBinError(2,9,1.6);
  h->SetBinContent(2,10,8.74);h->SetBinError(2,10,1.58);
  h->SetBinContent(2,11,8.66);h->SetBinError(2,11,1.56);
  h->SetBinContent(2,12,9.04);h->SetBinError(2,12,1.62);


  h->SetBinContent(3,1,10.0);h->SetBinError(3,1,1.85);
  h->SetBinContent(3,2,10.22);h->SetBinError(3,2,1.88);
  h->SetBinContent(3,3,10.05);h->SetBinError(3,3,1.83);
  h->SetBinContent(3,4,9.9);h->SetBinError(3,4,1.78);
  h->SetBinContent(3,5,9.75);h->SetBinError(3,5,1.76);
  h->SetBinContent(3,6,9.62);h->SetBinError(3,6,1.74);
  h->SetBinContent(3,7,9.51);h->SetBinError(3,7,1.71);
  h->SetBinContent(3,8,9.4);h->SetBinError(3,8,1.69);
  h->SetBinContent(3,9,9.3);h->SetBinError(3,9,1.69);
  h->SetBinContent(3,10,9.2);h->SetBinError(3,10,1.66);
  h->SetBinContent(3,11,9.11);h->SetBinError(3,11,1.65);
  h->SetBinContent(3,12,9.04);h->SetBinError(3,12,1.62);

  
  h->SetBinContent(4,1,10.49);h->SetBinError(4,1,1.93);
  h->SetBinContent(4,2,10.74);h->SetBinError(4,2,1.94);
  h->SetBinContent(4,3,10.57);h->SetBinError(4,3,1.92);
  h->SetBinContent(4,4,10.42);h->SetBinError(4,4,1.87);
  h->SetBinContent(4,5,10.28);h->SetBinError(4,5,1.87);
  h->SetBinContent(4,6,10.15);h->SetBinError(4,6,1.83);
  h->SetBinContent(4,7,10.04);h->SetBinError(4,7,1.81);
  h->SetBinContent(4,8,9.93);h->SetBinError(4,8,1.80);
  h->SetBinContent(4,9,9.83);h->SetBinError(4,9,1.78);
  h->SetBinContent(4,10,9.74);h->SetBinError(4,10,1.76);
  h->SetBinContent(4,11,9.65);h->SetBinError(4,11,1.74);
  h->SetBinContent(4,12,9.58);h->SetBinError(4,12,1.73);
  
  h->SetBinContent(5,1,11.08);h->SetBinError(5,1,2.04);
  h->SetBinContent(5,2,11.35);h->SetBinError(5,2,2.05);
  h->SetBinContent(5,3,11.19);h->SetBinError(5,3,2.02);
  h->SetBinContent(5,4,11.04);h->SetBinError(5,4,1.99);
  h->SetBinContent(5,5,10.9);h->SetBinError(5,5,1.96);
  h->SetBinContent(5,6,10.77);h->SetBinError(5,6,1.94);
  h->SetBinContent(5,7,10.66);h->SetBinError(5,7,1.92);
  h->SetBinContent(5,8,10.56);h->SetBinError(5,8,1.9);
  h->SetBinContent(5,9,10.46);h->SetBinError(5,9,1.88);
  h->SetBinContent(5,10,10.37);h->SetBinError(5,10,1.87);
  h->SetBinContent(5,11,10.29);h->SetBinError(5,11,1.85);
  h->SetBinContent(5,12,10.21);h->SetBinError(5,12,1.85);
  
  h->SetBinContent(6,1,11.78);h->SetBinError(6,1,2.26);
  h->SetBinContent(6,2,12.09);h->SetBinError(6,2,2.18);
  h->SetBinContent(6,3,11.93);h->SetBinError(6,3,2.14);
  h->SetBinContent(6,4,11.78);h->SetBinError(6,4,2.12);
  h->SetBinContent(6,5,11.64);h->SetBinError(6,5,2.16);
  h->SetBinContent(6,6,11.52);h->SetBinError(6,6,2.09);
  h->SetBinContent(6,7,11.41);h->SetBinError(6,7,2.06);
  h->SetBinContent(6,8,11.30);h->SetBinError(6,8,2.04);
  h->SetBinContent(6,9,11.21);h->SetBinError(6,9,2.04);
  h->SetBinContent(6,10,11.12);h->SetBinError(6,10,2.0);
  h->SetBinContent(6,11,11.04);h->SetBinError(6,11,1.99);
  h->SetBinContent(6,12,10.97);h->SetBinError(6,12,1.99);
  
  h->SetBinContent(7,1,12.28);h->SetBinError(7,1,2.33);
  h->SetBinContent(7,2,12.61);h->SetBinError(7,2,2.3);
  h->SetBinContent(7,3,12.45);h->SetBinError(7,3,2.28);
  h->SetBinContent(7,4,12.3);h->SetBinError(7,4,2.26);
  h->SetBinContent(7,5,12.17);h->SetBinError(7,5,2.23);
  h->SetBinContent(7,6,12.06);h->SetBinError(7,6,2.21);
  h->SetBinContent(7,7,11.95);h->SetBinError(7,7,2.18);
  h->SetBinContent(7,8,11.85);h->SetBinError(7,8,2.16);
  h->SetBinContent(7,9,11.76);h->SetBinError(7,9,2.15);
  h->SetBinContent(7,10,11.68);h->SetBinError(7,10,2.14);
  h->SetBinContent(7,11,11.6);h->SetBinError(7,11,2.12);
  h->SetBinContent(7,12,11.53);h->SetBinError(7,12,2.11);
  
  h->SetBinContent(8,1,12.03);h->SetBinError(8,1,2.36);
  h->SetBinContent(8,2,12.32);h->SetBinError(8,2,2.38);
  h->SetBinContent(8,3,12.18);h->SetBinError(8,3,2.34);
  h->SetBinContent(8,4,12.05);h->SetBinError(8,4,2.39);
  h->SetBinContent(8,5,11.93);h->SetBinError(8,5,2.29);
  h->SetBinContent(8,6,11.83);h->SetBinError(8,6,2.27);
  h->SetBinContent(8,7,11.73);h->SetBinError(8,7,2.26);
  h->SetBinContent(8,8,11.64);h->SetBinError(8,8,2.24);
  h->SetBinContent(8,9,11.57);h->SetBinError(8,9,2.30);
  h->SetBinContent(8,10,11.49);h->SetBinError(8,10,2.21);
  h->SetBinContent(8,11,11.42);h->SetBinError(8,11,2.21);
  h->SetBinContent(8,12,11.36);h->SetBinError(8,12,2.21);

  h->GetXaxis()->SetTitle("Fit Low Edge");
  h->GetYaxis()->SetTitle("Fit High Edge");
  h->GetZaxis()->SetRangeUser(8.65,12.4);
  h->Draw("colz");
  c1->Print("Plots/Higgs/fitYieldsWe.png");
  c1->Print("Plots/Higgs/fitYieldsWe.pdf");


}

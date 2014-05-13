// To be able to compile the above code, insert the following preprocessor directives:
// #include "Math/QuantFuncMathCore.h"
// #include "TMath.h"
// #include "TGraphAsymmErrors.h"
// - Thanks Guillelmo Gomez-Ceballos for these.

{
   const double alpha = 1 - 0.6827;
   TH1D * h1 = new TH1D("h1","h1",50,-4,4);
   h1->FillRandom("gaus",100);
  
    TGraphAsymmErrors * g = new TGraphAsymmErrors(h1);
    g->SetMarkerSize(0.5);
    g->SetMarkerStyle (20);

    for (int i = 0; i < g->GetN(); ++i) {
       int N = g->GetY()[i];
       double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
       double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
//      double U =  (N==0) ?  ( ROOT::Math::gamma_quantile_c(alpha,N+1,1) ) :
//         ( ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) );
       g->SetPointEYlow(i, N-L);
       g->SetPointEYhigh(i, U-N);
   }
   g->Draw("AP");
}

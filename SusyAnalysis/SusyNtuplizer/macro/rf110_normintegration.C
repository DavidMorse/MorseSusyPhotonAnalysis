/////////////////////////////////////////////////////////////////////////
//
// 'BASIC FUNCTIONALITY' RooFit tutorial macro #110
// 
// Examples on normalization of p.d.f.s,
// integration of p.d.fs, construction
// of cumulative distribution functions from p.d.f.s
// in one dimension
//
// 07/2008 - Wouter Verkerke 
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooAbsReal.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
using namespace RooFit ;


void rf110_normintegration()
{
  // S e t u p   m o d e l 
  // ---------------------

  // Create observables x,y
  RooRealVar x("x","x",-10,10) ;

  // Create p.d.f. gaussx(x,-2,3) 
  RooRealVar mean("mean","mean",0);
  RooRealVar width("width","width",2,0,5);
  RooGaussian gx("gx","gx",x,mean,width) ;
  //RooPolynomial gx("gx","gx",x,RooArgList(RooConst(-.3),RooConst(1)));

  RooDataSet *data = gx.generate(x,20000);
  gx.fitTo(*data,Save());

  // Plot cdf of gx versus x
  RooPlot* frame = x.frame(Title("Gaussian p.d.f")) ;
  data->plotOn(frame);
  gx.plotOn(frame) ;



  // Draw plot on canvas
  new TCanvas("rf110_normintegration","rf110_normintegration",1100,720) ;
  gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.6) ;
  frame->Draw() ;


  // R e t r i e v e   r a w  &   n o r m a l i z e d   v a l u e s   o f   R o o F i t   p . d . f . s
  // --------------------------------------------------------------------------------------------------

  // Return 'raw' unnormalized value of gx
  cout << "gx = " << gx.getVal() << endl ;
  
  // Return value of gx normalized over x in range [-10,10]
  RooArgSet nset(x) ;
  cout << "gx_Norm[x] = " << gx.getVal(&nset) << endl ;

  // Create object representing integral over gx
  // which is used to calculate  gx_Norm[x] == gx / gx_Int[x]
  RooAbsReal* igx = gx.createIntegral(x) ;
  cout << "gx_Int[x] = " << igx->getVal() << endl ;


  // I n t e g r a t e   n o r m a l i z e d   p d f   o v e r   s u b r a n g e
  // ----------------------------------------------------------------------------

  // Define a range named "signal" in x from -5,5
  x.setRange("signal",0,10) ;
  x.setRange("full",-10,10) ;
  
  // Create an integral of gx_Norm[x] over x in range "signal"
  // This is the fraction of of p.d.f. gx_Norm[x] which is in the
  // range named "signal"
  RooAbsReal* igx_sig = gx.createIntegral(x,NormSet(x),Range("signal")) ;
  RooAbsReal* igx_full = gx.createIntegral(x,NormSet(x),Range("full")) ;
  RooAbsReal* igx_fullNoNorm = gx.createIntegral(x,Range("full")) ;
  RooAbsReal* igx_sig2 = gx.createIntegral(x,Range("signal")) ;
  cout << "gx_Int[x|signal]_Norm[x]  = " << igx_sig->getVal() << endl ;
  cout << "gx_Int[x|full]_Norm[x]    = " << igx_full->getVal() << endl ;
  cout << "gx_Int[x|full]_NoNorm[x]  = " << igx_fullNoNorm->getVal() << endl ;
  cout << "gx_Int[x|signal]_Norm[x]2 = " << igx_sig2->getVal() << endl ;
  cout << "gx_Int[x|signal]_Norm[x] full  = " << igx_sig->getVal()*igx->getVal() << endl ;



  // C o n s t r u c t   c u m u l a t i v e   d i s t r i b u t i o n   f u n c t i o n   f r o m   p d f
  // -----------------------------------------------------------------------------------------------------

  // Create the cumulative distribution function of gx
  // i.e. calculate Int[-10,x] gx(x') dx'
  //RooAbsReal* gx_cdf = gx.createCdf(x) ;
  

}

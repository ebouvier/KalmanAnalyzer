#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooHist.h"

using namespace RooFit;

//---------------------------------------------------------------
void doTheFit()
//---------------------------------------------------------------
{
  gROOT->SetStyle("Plain");

  TFile *_file0 = TFile::Open("analyzeDmeson.root");
  gDirectory->cd("demo");
  TH1D* histo = (TH1D*)gDirectory->Get("HMASS4");

  RooRealVar x("x","D0 mass (GeV)",1.75,1.95);
  RooDataHist dh("dh","dh",x,histo);
  
  RooRealVar mean("mean","mean",1.86,1.8,2.0);
  RooRealVar sigma("sigma","sigma",0.005,0.001,0.200);

  RooRealVar lambda("lambda", "slope", -0.1, -5., 0.);

  RooGaussian sig("sig","sig",x,mean,sigma);
  RooExponential bck("bck", "exponential PDF", x, lambda);

  RooRealVar fsig("fsig","signal fraction 1",0.10,0.,1.);
  //  fsig.setConstant(kTRUE);

  // Signal+Background pdf
  RooAddPdf model("model","model",RooArgList(sig,bck),RooArgList(fsig)) ;
  // Fit and Plot

  RooPlot* frame = x.frame();
  
  model.fitTo(dh);
  dh.plotOn(frame);
  model.plotOn(frame);

  mass     = mean.getVal();
  masserr  = mean.getError();
  width    = sigma.getVal();
  widtherr = sigma.getError();

  mean.Print();
  sigma.Print();
  
  // Overlay the background component of model with a dashed line
  //  model.plotOn(frame,Components(bck),LineStyle(kDashed)) ;
  model.plotOn(frame,Components(bck),LineColor(kRed)) ;

  // Overlay the background+sig2 components of model with a dotted line
  //  model.plotOn(frame,Components(RooArgSet(sig)),LineStyle(kDotted)) ;
  model.plotOn(frame,Components(RooArgSet(sig)),LineColor(kGreen)) ;

  // Overlay the background+sig2 components of model with a dotted line
  //  model.plotOn(frame,Components(RooArgSet(sig2)),LineStyle(kDotted)) ;

  TCanvas* c = new TCanvas("D0 meson candidate","D0 meson candidate");

  gPad->SetLeftMargin(0.15); 
  frame->GetYaxis()->SetTitleOffset(1.4);
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->Draw();

  c->Update();

}


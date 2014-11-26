#include <iomanip>
#include <string>
#include <sstream>
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooLandau.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"

#pragma once

#define TITLE_FONTSIZE 26
#define LABEL_FONTSIZE 18

#define LEFT_MARGIN 0.17
#define RIGHT_MARGIN 0.03
#define TOP_MARGIN 0.05
#define BOTTOM_MARGIN 0.13

TStyle* createStyle() {
  TStyle *style = new TStyle("style", "style");

  // For the canvas:
  style->SetCanvasBorderMode(0);
  style->SetCanvasColor(kWhite);
  style->SetCanvasDefH(800); //Height of canvas
  style->SetCanvasDefW(800); //Width of canvas
  style->SetCanvasDefX(0);   //POsition on screen
  style->SetCanvasDefY(0);

  // For the Pad:
  style->SetPadBorderMode(0);
  style->SetPadColor(kWhite);
  style->SetPadGridX(false);
  style->SetPadGridY(false);
  style->SetGridColor(0);
  style->SetGridStyle(3);
  style->SetGridWidth(1);

  // For the frame:
  style->SetFrameBorderMode(0);
  style->SetFrameBorderSize(1);
  style->SetFrameFillColor(0);
  style->SetFrameFillStyle(0);
  style->SetFrameLineColor(1);
  style->SetFrameLineStyle(1);
  style->SetFrameLineWidth(1);

  // For the histo:
  style->SetHistLineStyle(1);
  style->SetHistLineWidth(2);
  style->SetEndErrorSize(2);
  style->SetMarkerStyle(20);

  //For the fit/function:
  style->SetFitFormat("5.4g");
  style->SetFuncColor(2);
  style->SetFuncStyle(1);
  style->SetFuncWidth(1);

  // For the statistics box:
  style->SetOptFile(0);
  style->SetStatColor(kWhite);
  //style->SetStatFont(43);
  //style->SetStatFontSize(0.025);
  style->SetStatTextColor(1);
  style->SetStatFormat("6.4g");
  style->SetStatBorderSize(1);
  //style->SetStatH(0.1);
  //style->SetStatW(0.15);

  //For the date:
  style->SetOptDate(0);

  // Margins:
  style->SetPadTopMargin(TOP_MARGIN);
  style->SetPadBottomMargin(BOTTOM_MARGIN);
  style->SetPadLeftMargin(LEFT_MARGIN);
  style->SetPadRightMargin(RIGHT_MARGIN);

  // For the Global title:
  style->SetOptTitle(0);
  style->SetTitleFont(63);
  style->SetTitleColor(1);
  style->SetTitleTextColor(1);
  style->SetTitleFillColor(10);
  style->SetTitleBorderSize(0);
  style->SetTitleAlign(33); 
  style->SetTitleX(1);
  style->SetTitleFontSize(TITLE_FONTSIZE);

  // For the axis titles:

  style->SetTitleColor(1, "XYZ");
  style->SetTitleFont(43, "XYZ");
  style->SetTitleSize(TITLE_FONTSIZE, "XYZ");
  style->SetTitleYOffset(2.5); 
  style->SetTitleXOffset(1.5);

  style->SetLabelColor(1, "XYZ");
  style->SetLabelFont(43, "XYZ");
  style->SetLabelOffset(0.01, "YZ");
  style->SetLabelOffset(0.015, "X");
  style->SetLabelSize(LABEL_FONTSIZE, "XYZ");

  style->SetAxisColor(1, "XYZ");
  style->SetStripDecimals(kTRUE);
  style->SetTickLength(0.03, "XYZ");
  style->SetNdivisions(510, "XYZ");
  style->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  style->SetPadTickY(1);

  style->SetOptLogx(0);
  style->SetOptLogy(0);
  style->SetOptLogz(0);

  style->SetHatchesSpacing(1.3);
  style->SetHatchesLineWidth(1);

  style->cd();

  return style;
}

void leg_style(TLegend *leg,
    //int text_size = 0.035,
    //int text_font = 43,
    int text_align = 22,
    int fill_style = 1,
    int fill_color = 10,
    int line_color = 0,
    int line_width = 0,
    int border_size = 1) {
  //leg->SetTextSize(text_size);
  //leg->SetTextFont(text_font);
  leg->SetTextAlign(text_align);
  leg->SetFillStyle(fill_style);
  leg->SetFillColor(fill_color);
  leg->SetLineColor(line_color);
  leg->SetLineWidth(line_width);
  leg->SetBorderSize(border_size);
}

void h1_style(TH1 *h,
    int line_color=1,
    int fill_color=50,
    int fill_style=1001,
    float y_min=-1111.,
    float y_max=-1111.,
    int ndivx=510,
    int ndivy=510,
    int marker_color=1,
    float marker_size=1.2,
    int optstat=0,
    const char* xtitle="") {

  h->SetLineColor(line_color);
  h->SetFillColor(fill_color);
  h->SetFillStyle(fill_style);
  h->SetMaximum(y_max);
  h->SetMinimum(y_min);
  h->GetXaxis()->SetNdivisions(ndivx);
  h->GetYaxis()->SetNdivisions(ndivy);
  h->GetYaxis()->SetTitleOffset(2.5);

  h->SetMarkerColor(marker_color);
  h->SetMarkerSize(marker_size);
  h->SetStats(optstat);

  h->GetXaxis()->SetTitle(xtitle);
  float binSize = h->GetXaxis()->GetBinWidth(1);
  std::stringstream ss;
  ss << "Events / " << std::fixed << std::setprecision(2) << binSize;
  h->GetYaxis()->SetTitle(ss.str().c_str());
}

void grapherrors_style(TGraphErrors *gr,
    const char* name="",
    int line_width=2,
    int line_color=1,
    int line_style=1,
    int fill_color=50,
    int fill_style=1001,
    float y_min=-1111.,
    float y_max=-1111.,
    int ndivx=510,
    int ndivy=510,
    int marker_style=20,
    int marker_color=1,
    float marker_size=1.2,
    const char* title="",
    const char* xtitle="",
    const char* ytitle="") {

  gr->SetLineWidth(line_width);
  gr->SetLineColor(line_color);
  gr->SetLineStyle(line_style);
  gr->SetFillColor(fill_color);
  gr->SetFillStyle(fill_style);
  gr->SetMaximum(y_max);
  gr->SetMinimum(y_min);
  gr->GetXaxis()->SetNdivisions(ndivx);
  gr->GetYaxis()->SetNdivisions(ndivy);

  gr->SetMarkerStyle(marker_style);
  gr->SetMarkerColor(marker_color);
  gr->SetMarkerSize(marker_size);

  gr->GetXaxis()->SetTitle(xtitle);
  gr->GetYaxis()->SetTitle(ytitle);
  gr->SetTitle(title);
  gr->SetName(name);
}

void cms_style(bool isData = false){
  std::string status = "Simulation preliminary";
  if (isData) status = "Preliminary";
  TPaveText* pt_exp = new TPaveText(LEFT_MARGIN, 1 - 0.5 * TOP_MARGIN, 1 - RIGHT_MARGIN, 1, "brNDC");
  pt_exp->SetFillStyle(0);
  pt_exp->SetBorderSize(0);
  pt_exp->SetMargin(0);
  pt_exp->SetTextFont(62);
  pt_exp->SetTextSize(0.75 * TOP_MARGIN);
  pt_exp->SetTextAlign(13);
  TString d = TString::Format("CMS #font[52]{#scale[0.76]{%s}}", status.c_str());
  pt_exp->AddText(d);
  pt_exp->Draw();

  TString lumi_s = "19.8 fb^{-1} (8 TeV)";
  TPaveText* pt_lumi = new TPaveText(LEFT_MARGIN, 1 - 0.5 * TOP_MARGIN, 1 - RIGHT_MARGIN, 1, "brNDC");
  pt_lumi->SetFillStyle(0);
  pt_lumi->SetBorderSize(0);
  pt_lumi->SetMargin(0);
  pt_lumi->SetTextFont(42);
  pt_lumi->SetTextSize(0.6 * TOP_MARGIN);
  pt_lumi->SetTextAlign(33);
  pt_lumi->AddText(lumi_s);
  pt_lumi->Draw();

}

//---------------------------------------------------------------
double *step(bool inBatch, TFile* fi, TString name, TString date)
//---------------------------------------------------------------
{
  using namespace RooFit;

  TH1D* histo = (TH1D*)fi->Get("ana/h_"+name);

  RooRealVar x("mass","D^{0} mass",1.7,2.,"GeV/c^{2}");
  RooDataHist dh("datahist","datahist",RooArgList(x),histo,1.);
  
  // Signal+Background pdf
  RooRealVar mean("mean","mean",1.86484,1.84554,1.88414);
  RooRealVar sigma("sigma","sigma",0.0193,0.01,0.03);
  RooRealVar lambda("lambda", "slope", -1.8, -2.2, -0.7);
  RooGaussian sig("sig","gaussian signal",x,mean,sigma);
  RooExponential bck("bck", "exponential background", x, lambda);
  RooRealVar nsig("nsig","number of signal events",0.1*histo->Integral(),0.,histo->Integral());
  RooRealVar nbck("nbck","number of background events",0.9*histo->Integral(),0.,histo->Integral());
  RooAddPdf model("model","model",RooArgList(sig,bck),RooArgList(nsig,nbck)) ;

  // Fit
  model.fitTo(dh);
  double mass     = mean.getVal();
  double masserr  = mean.getError();
  double width    = sigma.getVal();
  double widtherr = sigma.getError();
  double Nsig  = nsig.getVal();
  double dNsig = nsig.getError();
  double Nbck  = nbck.getVal();
  double dNbck = nbck.getError();
  double *SNR = new double[2];
  SNR[0] = Nsig*pow(Nbck,-0.5);
  SNR[1]  = dNsig*pow(Nbck,-0.5) + 0.5*SNR[0]*dNbck/Nbck;
  std::cout << "\nSignal mean = " << mass << " +/- " << masserr << std::endl;
  std::cout << "Signal width = " << width << " +/- " << widtherr << "\n" << std::endl;
  std::cout << "Nsig = " << Nsig << " +/- " << dNsig << std::endl;
  std::cout << "Nbck = " << Nbck << " +/- " << dNbck << std::endl;
  std::cout << "SNR  = " << SNR[0] << " +/- " << SNR[1] << "\n" << std::endl;

  TPaveText* fit_tex = new TPaveText(0.21,0.21,0.51,0.46,"BRNDC");
  fit_tex->AddText("Gaussian parameters :");
  fit_tex->AddText(TString::Format("#mu = (%2.4f #pm %2.4f) GeV/c^{2}",mass,masserr));
  fit_tex->AddText(TString::Format("#sigma = (%2.4f #pm %2.4f) GeV/c^{2}",width,widtherr));
  fit_tex->AddText("");
  fit_tex->AddText(TString::Format("N_{sig} = %5.0f #pm %3.0f", Nsig, dNsig));
  fit_tex->AddText(TString::Format("N_{bck} = %5.0f #pm %3.0f", Nbck, dNbck));
  fit_tex->AddText(TString::Format("N_{sig}/#sqrt{N_{bck}} = %2.1f #pm %2.1f", SNR[0], SNR[1]));
  fit_tex->SetTextFont(43);
  fit_tex->SetFillColor(0);
  fit_tex->SetTextSize(TITLE_FONTSIZE - 6);

  // Plot
  RooPlot* frame = x.frame();
  dh.plotOn(frame);
  model.plotOn(frame,LineColor(9));
  model.plotOn(frame,Components(bck),LineColor(kBlue));
  model.plotOn(frame,Components(RooArgSet(sig)),LineColor(kRed));

  TCanvas* cn = new TCanvas("cn_"+name,"cn_"+name,800,800);
  frame->Draw();
  if (Nsig >= 1) fit_tex->Draw("same");
  cms_style(); 
  cn->SaveAs("Plots"+date+"/fit_"+name+".C");
  cn->SaveAs("Plots"+date+"/fit_"+name+".pdf");
  cn->SaveAs("Plots"+date+"/fit_"+name+".eps");
  if (!inBatch) getchar();

  delete cn; delete frame; delete fit_tex; delete histo;
  return SNR;
}

//---------------------------------------------------------------
double *BRstep(bool inBatch, TFile* fi, TString name, TString date)
//---------------------------------------------------------------
{
  using namespace RooFit;

  TH1D* histo = (TH1D*)fi->Get("ana/h_"+name);

  RooRealVar x("mass","D^{0} mass",0.8,2.5,"GeV/c^{2}");
  RooDataHist dh("datahist","datahist",RooArgList(x),histo,1.);
  
  // Signal+Background pdf
  RooRealVar mean("mean","mean",1.86484,1.84554,1.88414);
  RooRealVar widthbr("width","width",0.03,0.,0.1);
  RooRealVar sigma("sigma","sigma",0.,-0.1,0.1);
  RooRealVar lambda("lambda", "slope", -1e-10, -10, 0.);
  RooVoigtian sig("sig","voigtian signal",x,mean,widthbr,sigma);
  RooPolynomial bck("bck", "linear background", x, RooArgList(lambda), 1);
  RooRealVar nsig("nsig","number of signal events",0.1*histo->Integral(),0.,histo->Integral());
  RooRealVar nbck("nbck","number of background events",0.9*histo->Integral(),0.,histo->Integral());
  RooAddPdf model("model","model",RooArgList(sig,bck),RooArgList(nsig,nbck)) ;

  // Fit
  model.fitTo(dh);
  double mass     = mean.getVal();
  double masserr  = mean.getError();
  double width    = sigma.getVal();
  double widtherr = sigma.getError();
  double Nsig  = nsig.getVal();
  double dNsig = nsig.getError();
  double Nbck  = nbck.getVal();
  double dNbck = nbck.getError();
  double *SNR = new double[2];
  SNR[0] = Nsig*pow(Nbck,-0.5);
  SNR[1]  = dNsig*pow(Nbck,-0.5) + 0.5*SNR[0]*dNbck/Nbck;
  std::cout << "\nSignal mean = " << mass << " +/- " << masserr << std::endl;
  std::cout << "Signal width = " << width << " +/- " << widtherr << "\n" << std::endl;
  std::cout << "Nsig = " << Nsig << " +/- " << dNsig << std::endl;
  std::cout << "Nbck = " << Nbck << " +/- " << dNbck << std::endl;
  std::cout << "SNR  = " << SNR[0] << " +/- " << SNR[1] << "\n" << std::endl;

  TPaveText* fit_tex = new TPaveText(0.51,0.61,0.81,0.86,"BRNDC");
  fit_tex->AddText("Gaussian parameters :");
  fit_tex->AddText(TString::Format("#mu = (%2.4f #pm %2.4f) GeV/c^{2}",mass,masserr));
  fit_tex->AddText(TString::Format("#sigma = (%2.4f #pm %2.4f) GeV/c^{2}",width,widtherr));
  fit_tex->AddText("");
  fit_tex->AddText(TString::Format("N_{sig} = %5.0f #pm %3.0f", Nsig, dNsig));
  fit_tex->AddText(TString::Format("N_{bck} = %5.0f #pm %3.0f", Nbck, dNbck));
  fit_tex->AddText(TString::Format("N_{sig}/#sqrt{N_{bck}} = %2.1f #pm %2.1f", SNR[0], SNR[1]));
  fit_tex->SetTextFont(43);
  fit_tex->SetFillColor(0);
  fit_tex->SetTextSize(TITLE_FONTSIZE - 6);

  // Plot
  RooPlot* frame = x.frame();
  dh.plotOn(frame);
  model.plotOn(frame,LineColor(9));
  model.plotOn(frame,Components(bck),LineColor(kBlue));
  model.plotOn(frame,Components(RooArgSet(sig)),LineColor(kRed));

  TCanvas* cn = new TCanvas("cn_"+name,"cn_"+name,800,800);
  frame->Draw();
  if (Nsig >= 1) fit_tex->Draw("same");
  cms_style(); 
  cn->SaveAs("Plots"+date+"/fit_"+name+".C");
  cn->SaveAs("Plots"+date+"/fit_"+name+".pdf");
  cn->SaveAs("Plots"+date+"/fit_"+name+".eps");
  if (!inBatch) getchar();

  delete cn; delete frame; delete fit_tex; delete histo;
  return SNR;
}

//---------------------------------------------------------------
int doTheFit(bool inBatch = true, TString date = "")
//---------------------------------------------------------------
{
  using namespace RooFit;
  TStyle* m_style = createStyle();
  m_style->cd();
  if (inBatch) gROOT->SetBatch(true);
  if (date.Length() > 0) date = "_" + date;
  gROOT->ProcessLine(".! mkdir Plots"+date);

  //TFile *fi = TFile::Open("../test/test/kalmanAnalyzed_1p7to2.root"); //For tests
  //TFile *fi = TFile::Open("../test/test/kalmanAnalyzed_3sigma.root"); //For tests
  //TFile *fi = TFile::Open("../test/test/kalmanAnalyzed_noskim.root"); //For tests
  //TFile *fi = TFile::Open("../test/test/kalmanAnalyzed_141124.root"); //For tests
  //TFile *fi = TFile::Open("../test/test/kalmanAnalyzed_141125_pTcut.root"); //For tests
  TFile *fi = TFile::Open("../test/kalmanAnalyzed.root"); 

  //=============================================================================================
  //                  Simple Kalman Vertex Fitter for the D0 -> K Pi reconstruction
  //=============================================================================================

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // optimize chi2 cut on D0 -> K Pi reconstruction
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  const unsigned int numberOfPoints = 15;
  double *SNR_D0Cand_MassChi2Inf1 = step(inBatch, fi, "D0Cand_MassChi2Inf1", date); 
  double *SNR_D0Cand_MassChi2Inf1p5 = step(inBatch, fi, "D0Cand_MassChi2Inf1p5", date); 
  double *SNR_D0Cand_MassChi2Inf2 = step(inBatch, fi, "D0Cand_MassChi2Inf2", date);
  double *SNR_D0Cand_MassChi2Inf2p5 = step(inBatch, fi, "D0Cand_MassChi2Inf2p5", date); 
  double *SNR_D0Cand_MassChi2Inf3 = step(inBatch, fi, "D0Cand_MassChi2Inf3", date); 
  double *SNR_D0Cand_MassChi2Inf3p5 = step(inBatch, fi, "D0Cand_MassChi2Inf3p5", date); 
  double *SNR_D0Cand_MassChi2Inf4 = step(inBatch, fi, "D0Cand_MassChi2Inf4", date); 
  double *SNR_D0Cand_MassChi2Inf4p5 = step(inBatch, fi, "D0Cand_MassChi2Inf4p5", date); 
  double *SNR_D0Cand_MassChi2Inf5 = step(inBatch, fi, "D0Cand_MassChi2Inf5", date); 
  double *SNR_D0Cand_MassChi2Inf5p5 = step(inBatch, fi, "D0Cand_MassChi2Inf5p5", date); 
  double *SNR_D0Cand_MassChi2Inf6 = step(inBatch, fi, "D0Cand_MassChi2Inf6", date); 
  double *SNR_D0Cand_MassChi2Inf6p5 = step(inBatch, fi, "D0Cand_MassChi2Inf6p5", date); 
  double *SNR_D0Cand_MassChi2Inf7 = step(inBatch, fi, "D0Cand_MassChi2Inf7", date); 
  double *SNR_D0Cand_MassChi2Inf7p5 = step(inBatch, fi, "D0Cand_MassChi2Inf7p5", date); 
  double *SNR_D0Cand_MassChi2Inf8 = step(inBatch, fi, "D0Cand_MassChi2Inf8", date); 

  double x[numberOfPoints] = {1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.};
  double ex[numberOfPoints] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double y[numberOfPoints] = {SNR_D0Cand_MassChi2Inf1[0], SNR_D0Cand_MassChi2Inf1p5[0], SNR_D0Cand_MassChi2Inf2[0], SNR_D0Cand_MassChi2Inf2p5[0], SNR_D0Cand_MassChi2Inf3[0], SNR_D0Cand_MassChi2Inf3p5[0], SNR_D0Cand_MassChi2Inf4[0], SNR_D0Cand_MassChi2Inf4p5[0], SNR_D0Cand_MassChi2Inf5[0], SNR_D0Cand_MassChi2Inf5p5[0], SNR_D0Cand_MassChi2Inf6[0], SNR_D0Cand_MassChi2Inf6p5[0], SNR_D0Cand_MassChi2Inf7[0], SNR_D0Cand_MassChi2Inf7p5[0], SNR_D0Cand_MassChi2Inf8[0]};
  double ey[numberOfPoints] = {SNR_D0Cand_MassChi2Inf1[1], SNR_D0Cand_MassChi2Inf1p5[1], SNR_D0Cand_MassChi2Inf2[1], SNR_D0Cand_MassChi2Inf2p5[1], SNR_D0Cand_MassChi2Inf3[1], SNR_D0Cand_MassChi2Inf3p5[1], SNR_D0Cand_MassChi2Inf4[1], SNR_D0Cand_MassChi2Inf4p5[1], SNR_D0Cand_MassChi2Inf5[1], SNR_D0Cand_MassChi2Inf5p5[1], SNR_D0Cand_MassChi2Inf6[1], SNR_D0Cand_MassChi2Inf6p5[1], SNR_D0Cand_MassChi2Inf7[1], SNR_D0Cand_MassChi2Inf7p5[1], SNR_D0Cand_MassChi2Inf8[1]};

  TGraphErrors *gr_SNR_D0Cand = new TGraphErrors(numberOfPoints,x,y,ex,ey);
  grapherrors_style(gr_SNR_D0Cand, "SNR_D0Cand", 2, 30, 1, 30, 1001, -1111., -1111., 510, 510, 21, 36, 1., "D^{0} candidates", "#chi^{2}/NDOF (D^{0}#rightarrow#kappa#pi)", "N_{sig}/#sqrt{N_{bck}} (D^{0}#rightarrow#kappa#pi)");
  TCanvas *cn_SNR_D0Cand = new TCanvas("cn_SNR_D0Cand", "cn_SNR_D0Cand", 800, 800);
  cn_SNR_D0Cand->cd();
  gr_SNR_D0Cand->Draw("AP");
  cms_style();
  cn_SNR_D0Cand->SaveAs("Plots"+date+"/SNR_D0Cand.C");
  cn_SNR_D0Cand->SaveAs("Plots"+date+"/SNR_D0Cand.eps");
  cn_SNR_D0Cand->SaveAs("Plots"+date+"/SNR_D0Cand.pdf");

  double *SNR_D0_Mass = step(inBatch, fi, "D0_Mass", date);
  double *SNR_B_D0Mass = step(inBatch, fi, "B_D0Mass", date);

  // Print chi2 otpimization results 
  std::cout << "\n\nD0 candidate mass SNR for \n" << std::endl;
  std::cout << "Chi2 < 1   : " <<  SNR_D0Cand_MassChi2Inf1[0] << " +/- " << SNR_D0Cand_MassChi2Inf1[1] << std::endl;
  std::cout << "Chi2 < 1.5 : " <<  SNR_D0Cand_MassChi2Inf1p5[0] << " +/- " << SNR_D0Cand_MassChi2Inf1p5[1] << std::endl;
  std::cout << "Chi2 < 2   : " <<  SNR_D0Cand_MassChi2Inf2[0] << " +/- " << SNR_D0Cand_MassChi2Inf2[1] << std::endl;
  std::cout << "Chi2 < 2.5 : " <<  SNR_D0Cand_MassChi2Inf2p5[0] << " +/- " << SNR_D0Cand_MassChi2Inf2p5[1] << std::endl;
  std::cout << "Chi2 < 3   : " <<  SNR_D0Cand_MassChi2Inf3[0] << " +/- " << SNR_D0Cand_MassChi2Inf3[1] << std::endl;
  std::cout << "Chi2 < 3.5 : " <<  SNR_D0Cand_MassChi2Inf3p5[0] << " +/- " << SNR_D0Cand_MassChi2Inf3p5[1] << std::endl;
  std::cout << "Chi2 < 4   : " <<  SNR_D0Cand_MassChi2Inf4[0] << " +/- " << SNR_D0Cand_MassChi2Inf4[1] << std::endl;
  std::cout << "... and ctau/sigma(ctau) > 50 : " << SNR_D0_Mass[0] << " +/- " << SNR_D0_Mass[1] << std::endl;
  std::cout << "... and non iso mu : " << SNR_B_D0Mass[0] << " +/- " << SNR_B_D0Mass[1] << std::endl;
  std::cout << "Chi2 < 4.5 : " <<  SNR_D0Cand_MassChi2Inf4p5[0] << " +/- " << SNR_D0Cand_MassChi2Inf4p5[1] << std::endl;
  std::cout << "Chi2 < 5   : " <<  SNR_D0Cand_MassChi2Inf5[0] << " +/- " << SNR_D0Cand_MassChi2Inf5[1] << std::endl;
  std::cout << "Chi2 < 5.5 : " <<  SNR_D0Cand_MassChi2Inf5p5[0] << " +/- " << SNR_D0Cand_MassChi2Inf5p5[1] << std::endl;
  std::cout << "Chi2 < 6   : " <<  SNR_D0Cand_MassChi2Inf6[0] << " +/- " << SNR_D0Cand_MassChi2Inf6[1] << std::endl;
  std::cout << "Chi2 < 6.5 : " <<  SNR_D0Cand_MassChi2Inf6p5[0] << " +/- " << SNR_D0Cand_MassChi2Inf6p5[1] << std::endl;
  std::cout << "Chi2 < 7   : " <<  SNR_D0Cand_MassChi2Inf7[0] << " +/- " << SNR_D0Cand_MassChi2Inf7[1] << std::endl;
  std::cout << "Chi2 < 7.5 : " <<  SNR_D0Cand_MassChi2Inf7p5[0] << " +/- " << SNR_D0Cand_MassChi2Inf7p5[1] << std::endl;
  std::cout << "Chi2 < 8   : " <<  SNR_D0Cand_MassChi2Inf8[0] << " +/- " << SNR_D0Cand_MassChi2Inf8[1] << std::endl;

  if (!inBatch) getchar();
  delete cn_SNR_D0Cand; delete gr_SNR_D0Cand;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // have a look at other distributions 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  TH1D* h_nJets = (TH1D*)fi->Get("ana/h_nJets");
  h1_style(h_nJets, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "Number of jets");
  TCanvas* cn_nJets = new TCanvas("cn_nJets","cn_nJets",800,800);
  cn_nJets->cd();
  h_nJets->Draw("hist");
  cms_style();
  cn_nJets->SaveAs("Plots"+date+"/nJets.C");
  cn_nJets->SaveAs("Plots"+date+"/nJets.eps");
  cn_nJets->SaveAs("Plots"+date+"/nJets.pdf");

  TH1D* h_CSV = (TH1D*)fi->Get("ana/h_CSV");
  h1_style(h_CSV, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "CSV discriminant");
  TCanvas* cn_CSV = new TCanvas("cn_CSV","cn_CSV",800,800);
  cn_CSV->cd();
  h_CSV->Draw("hist");
  cms_style();
  cn_CSV->SaveAs("Plots"+date+"/CSV.C");
  cn_CSV->SaveAs("Plots"+date+"/CSV.eps");
  cn_CSV->SaveAs("Plots"+date+"/CSV.pdf");

  TH1D* h_B_cuts = (TH1D*)fi->Get("ana/h_B_cuts");
  h_B_cuts->GetXaxis()->SetRangeUser(0,11);
  h1_style(h_B_cuts, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "");
  TCanvas* cn_B_cuts = new TCanvas("cn_B_cuts","cn_B_cuts",800,800);
  cn_B_cuts->cd();
  cn_B_cuts->SetLogy(1);
  h_B_cuts->Draw("hist");
  cms_style();
  cn_B_cuts->SaveAs("Plots"+date+"/B_cuts.C");
  cn_B_cuts->SaveAs("Plots"+date+"/B_cuts.eps");
  cn_B_cuts->SaveAs("Plots"+date+"/B_cuts.pdf");

  TH1D* h_D0_dRJet = (TH1D*)fi->Get("ana/h_D0_dRJet");
  h_D0_dRJet->GetXaxis()->SetRangeUser(0,0.5);
  h1_style(h_D0_dRJet, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "#DeltaR(D^{0},jet)");
  TCanvas* cn_D0_dRJet = new TCanvas("cn_D0_dRJet","cn_D0_dRJet",800,800);
  cn_D0_dRJet->cd();
  h_D0_dRJet->Draw("hist");
  cms_style();
  cn_D0_dRJet->SaveAs("Plots"+date+"/D0_dRJet.C");
  cn_D0_dRJet->SaveAs("Plots"+date+"/D0_dRJet.eps");
  cn_D0_dRJet->SaveAs("Plots"+date+"/D0_dRJet.pdf");

  TH1D* h_D0_L = (TH1D*)fi->Get("ana/h_D0_L");
  h_D0_L->GetXaxis()->SetRangeUser(0,0.25);
  h1_style(h_D0_L, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "c#tau (D^{0}#rightarrow#kappa#pi) (cm)");
  TCanvas* cn_D0_L = new TCanvas("cn_D0_L","cn_D0_L",800,800);
  cn_D0_L->cd();
  h_D0_L->Draw("hist");
  cms_style();
  cn_D0_L->SaveAs("Plots"+date+"/D0_L.C");
  cn_D0_L->SaveAs("Plots"+date+"/D0_L.eps");
  cn_D0_L->SaveAs("Plots"+date+"/D0_L.pdf");

  TH1D* h_D0_SigmaL = (TH1D*)fi->Get("ana/h_D0_SigmaL");
  h_D0_SigmaL->GetXaxis()->SetRangeUser(0,0.001);
  h1_style(h_D0_SigmaL, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "#Delta(c#tau) (D^{0}#rightarrow#kappa#pi) (cm)");
  TCanvas* cn_D0_SigmaL = new TCanvas("cn_D0_SigmaL","cn_D0_SigmaL",800,800);
  cn_D0_SigmaL->cd();
  h_D0_SigmaL->Draw("hist");
  cms_style();
  cn_D0_SigmaL->SaveAs("Plots"+date+"/D0_SigmaL.C");
  cn_D0_SigmaL->SaveAs("Plots"+date+"/D0_SigmaL.eps");
  cn_D0_SigmaL->SaveAs("Plots"+date+"/D0_SigmaL.pdf");

  TH1D* h_D0_p = (TH1D*)fi->Get("ana/h_D0_p");
  h_D0_p->GetXaxis()->SetRangeUser(12.,250.);
  h1_style(h_D0_p, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "p(D^{0}#rightarrow#kappa#pi) (GeV/c)");
  TCanvas* cn_D0_p = new TCanvas("cn_D0_p","cn_D0_p",800,800);
  cn_D0_p->cd();
  h_D0_p->Draw("hist");
  cms_style();
  cn_D0_p->SaveAs("Plots"+date+"/D0_p.C");
  cn_D0_p->SaveAs("Plots"+date+"/D0_p.eps");
  cn_D0_p->SaveAs("Plots"+date+"/D0_p.pdf");

  TH1D* h_D0_pT = (TH1D*)fi->Get("ana/h_D0_pT");
  h_D0_pT->GetXaxis()->SetRangeUser(12.,120.);
  h1_style(h_D0_pT, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "p_{T}(D^{0}#rightarrow#kappa#pi) (GeV/c)");
  TCanvas* cn_D0_pT = new TCanvas("cn_D0_pT","cn_D0_pT",800,800);
  cn_D0_pT->cd();
  h_D0_pT->Draw("hist");
  cms_style();
  cn_D0_pT->SaveAs("Plots"+date+"/D0_pT.C");
  cn_D0_pT->SaveAs("Plots"+date+"/D0_pT.eps");
  cn_D0_pT->SaveAs("Plots"+date+"/D0_pT.pdf");

  TH1D* h_D0_eta = (TH1D*)fi->Get("ana/h_D0_eta");
  h1_style(h_D0_eta, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "#eta(D^{0}#rightarrow#kappa#pi)");
  TCanvas* cn_D0_eta = new TCanvas("cn_D0_eta","cn_D0_eta",800,800);
  cn_D0_eta->cd();
  h_D0_eta->Draw("hist");
  cms_style();
  cn_D0_eta->SaveAs("Plots"+date+"/D0_eta.C");
  cn_D0_eta->SaveAs("Plots"+date+"/D0_eta.eps");
  cn_D0_eta->SaveAs("Plots"+date+"/D0_eta.pdf");

  TH1D* h_D0_phi = (TH1D*)fi->Get("ana/h_D0_phi");
  h1_style(h_D0_phi, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "#phi(D^{0}#rightarrow#kappa#pi)");
  TCanvas* cn_D0_phi = new TCanvas("cn_D0_phi","cn_D0_phi",800,800);
  cn_D0_phi->cd();
  h_D0_phi->Draw("hist");
  cms_style();
  cn_D0_phi->SaveAs("Plots"+date+"/D0_phi.C");
  cn_D0_phi->SaveAs("Plots"+date+"/D0_phi.eps");
  cn_D0_phi->SaveAs("Plots"+date+"/D0_phi.pdf");

  TH1D* h_D0_Mass = (TH1D*)fi->Get("ana/h_D0_Mass");
  h_D0_Mass->GetXaxis()->SetRangeUser(1.7,2.);
  h1_style(h_D0_Mass, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "m(D^{0}#rightarrow#kappa#pi) (GeV/c^{2})");
  TCanvas* cn_D0_Mass = new TCanvas("cn_D0_Mass","cn_D0_Mass",800,800);
  cn_D0_Mass->cd();
  h_D0_Mass->Draw("hist");
  cms_style();
  cn_D0_Mass->SaveAs("Plots"+date+"/D0_Mass.C");
  cn_D0_Mass->SaveAs("Plots"+date+"/D0_Mass.eps");
  cn_D0_Mass->SaveAs("Plots"+date+"/D0_Mass.pdf");

  TH1D* h_B_D0Mass = (TH1D*)fi->Get("ana/h_B_D0Mass");
  h_B_D0Mass->GetXaxis()->SetRangeUser(1.7,2.);
  h1_style(h_B_D0Mass, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "m(D^{0}#rightarrow#kappa#pi) (GeV/c^{2})");
  TCanvas* cn_B_D0Mass = new TCanvas("cn_B_D0Mass","cn_B_D0Mass",800,800);
  cn_B_D0Mass->cd();
  h_B_D0Mass->Draw("hist");
  cms_style();
  cn_B_D0Mass->SaveAs("Plots"+date+"/B_D0Mass.C");
  cn_B_D0Mass->SaveAs("Plots"+date+"/B_D0Mass.eps");
  cn_B_D0Mass->SaveAs("Plots"+date+"/B_D0Mass.pdf");

  TH1D* h_B_DeltaRD0Mu = (TH1D*)fi->Get("ana/h_B_DeltaRD0Mu");
  h1_style(h_B_DeltaRD0Mu, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "#DeltaR(D^{0},#mu^{#pm})");
  TCanvas* cn_B_DeltaRD0Mu = new TCanvas("cn_B_DeltaRD0Mu","cn_B_DeltaRD0Mu",800,800);
  cn_B_DeltaRD0Mu->cd();
  h_B_DeltaRD0Mu->Draw("hist");
  cms_style();
  cn_B_DeltaRD0Mu->SaveAs("Plots"+date+"/B_DeltaRD0Mu.C");
  cn_B_DeltaRD0Mu->SaveAs("Plots"+date+"/B_DeltaRD0Mu.eps");
  cn_B_DeltaRD0Mu->SaveAs("Plots"+date+"/B_DeltaRD0Mu.pdf");

  TH1D* h_B_pD0pMu = (TH1D*)fi->Get("ana/h_B_pD0pMu");
  h_B_pD0pMu->GetXaxis()->SetRangeUser(0.1,10.);
  h1_style(h_B_pD0pMu, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "p(D^{0})/p(#mu^{#pm})");
  TCanvas* cn_B_pD0pMu = new TCanvas("cn_B_pD0pMu","cn_B_pD0pMu",800,800);
  cn_B_pD0pMu->cd();
  cn_B_pD0pMu->SetLogx();
  h_B_pD0pMu->Draw("hist");
  cms_style();
  cn_B_pD0pMu->SaveAs("Plots"+date+"/B_pD0pMu.C");
  cn_B_pD0pMu->SaveAs("Plots"+date+"/B_pD0pMu.eps");
  cn_B_pD0pMu->SaveAs("Plots"+date+"/B_pD0pMu.pdf");

  TH1D* h_B_pTD0pTMu = (TH1D*)fi->Get("ana/h_B_pTD0pTMu");
  h_B_pTD0pTMu->GetXaxis()->SetRangeUser(0.1,10.);
  h1_style(h_B_pTD0pTMu, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "p_{T}(D^{0})/p_{T}(#mu^{#pm})");
  TCanvas* cn_B_pTD0pTMu = new TCanvas("cn_B_pTD0pTMu","cn_B_pTD0pTMu",800,800);
  cn_B_pTD0pTMu->cd();
  cn_B_pTD0pTMu->SetLogx();
  h_B_pTD0pTMu->Draw("hist");
  cms_style();
  cn_B_pTD0pTMu->SaveAs("Plots"+date+"/B_pTD0pTMu.C");
  cn_B_pTD0pTMu->SaveAs("Plots"+date+"/B_pTD0pTMu.eps");
  cn_B_pTD0pTMu->SaveAs("Plots"+date+"/B_pTD0pTMu.pdf");

  TH1D* h_B_Mass = (TH1D*)fi->Get("ana/h_B_Mass");
  h_B_Mass->Rebin(5);
  h_B_Mass->GetXaxis()->SetRangeUser(1.5,6.);
  h1_style(h_B_Mass, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "m(D^{0}+#mu^{#pm}) (GeV/c^{2})");
  TCanvas* cn_B_Mass = new TCanvas("cn_B_Mass","cn_B_Mass",800,800);
  cn_B_Mass->cd();
  h_B_Mass->Draw("hist");
  cms_style();
  cn_B_Mass->SaveAs("Plots"+date+"/B_Mass.C");
  cn_B_Mass->SaveAs("Plots"+date+"/B_Mass.eps");
  cn_B_Mass->SaveAs("Plots"+date+"/B_Mass.pdf");
/*
  TH1D* h_B_DiracMass = (TH1D*)fi->Get("ana/h_B_DiracMass");
  h_B_DiracMass->Rebin(5);
  //h_B_DiracMass->GetXaxis()->SetRangeUser(FIXME,FIXME);
  h1_style(h_B_DiracMass, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "m(D^{0}+W^{#pm}) (GeV/c^{2})");
  TCanvas* cn_B_DiracMass = new TCanvas("cn_B_DiracMass","cn_B_DiracMass",800,800);
  cn_B_DiracMass->cd();
  h_B_DiracMass->Draw("hist");
  cms_style();
  cn_B_DiracMass->SaveAs("Plots"+date+"/B_DiracMass.C");
  cn_B_DiracMass->SaveAs("Plots"+date+"/B_DiracMass.eps");
  cn_B_DiracMass->SaveAs("Plots"+date+"/B_DiracMass.pdf");

  TH1D* h_B_GausMass = (TH1D*)fi->Get("ana/h_B_GausMass");
  h_B_GausMass->Rebin(5);
  //h_B_GausMass->GetXaxis()->SetRangeUser(FIXME,FIXME);
  h1_style(h_B_GausMass, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "m(D^{0}+W^{#pm}) (GeV/c^{2})");
  TCanvas* cn_B_GausMass = new TCanvas("cn_B_GausMass","cn_B_GausMass",800,800);
  cn_B_GausMass->cd();
  h_B_GausMass->Draw("hist");
  cms_style();
  cn_B_GausMass->SaveAs("Plots"+date+"/B_GausMass.C");
  cn_B_GausMass->SaveAs("Plots"+date+"/B_GausMass.eps");
  cn_B_GausMass->SaveAs("Plots"+date+"/B_GausMass.pdf");
*/
  if (!inBatch) getchar();
  delete cn_nJets; delete cn_CSV; delete cn_B_cuts; delete cn_D0_L; delete cn_D0_SigmaL; delete cn_D0_Mass; 
  delete cn_D0_dRJet; delete cn_D0_p; delete cn_D0_pT; delete cn_D0_eta; delete cn_D0_phi;  delete h_B_DeltaRD0Mu; delete cn_B_pD0pMu; delete cn_B_pTD0pTMu;
  delete cn_B_D0Mass; delete cn_B_Mass; //delete cn_B_DiracMass; delete cn_B_GausMass;
    
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Fit Bpm mass distributions 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  //~~~ D0+mu mass ~~~
  
  RooRealVar B_Mass("B_Mass","D^{0}#mu^{#pm} mass",1.5,6.,"GeV/c^{2}");
  RooDataHist dh_B_Mass("datahist_B_Mass","datahist_B_Mass",RooArgList(B_Mass),h_B_Mass,1.);
  
  // Signal+Background pdf
  //RooRealVar gausmean_B_Mass("gausmean_B_Mass","gausmean_B_Mass",2.6,2.55,2.65); //FIXME
  RooRealVar gausmean_B_Mass("gausmean_B_Mass","gausmean_B_Mass",3.75,2.5,5.);
  RooRealVar gaussig_B_Mass("gaussig_B_Mass","gaussig_B_Mass",0.0193,0.,0.05);
  //RooRealVar landmean_B_Mass("landmean_B_Mass","landmean_B_Mass",2.4,2.2,2.6); //FIXME
  //RooRealVar landsig_B_Mass("landsig_B_Mass","landsig_B_Mass",0.5,0.,1.); //FIXME
  //RooLandau bck_B_Mass("bck_B_Mass", "landau background", B_Mass, landmean_B_Mass, landsig_B_Mass); //FIXME
  RooRealVar alpha_bck("alpha_bck","alpha_bck",-1.5,-2.5,-0.5);
  RooExponential bck_B_Mass("bck_B_Mass","exp background",B_Mass,alpha_bck);
  RooGaussian sig_B_Mass("sig_B_Mass","gaussian signal",B_Mass,gausmean_B_Mass,gaussig_B_Mass);
  RooRealVar nsig_B_Mass("nsig_B_Mass","number of signal events",0.1*h_B_Mass->Integral(),0.,h_B_Mass->Integral());
  RooRealVar nbck_B_Mass("nbck_B_Mass","number of background events",0.9*h_B_Mass->Integral(),0.,h_B_Mass->Integral());
  RooAddPdf model_B_Mass("model_B_Mass","model_B_Mass",RooArgList(sig_B_Mass,bck_B_Mass),RooArgList(nsig_B_Mass,nbck_B_Mass)) ;

  // Fit
  model_B_Mass.fitTo(dh_B_Mass,Range(2.8,4.8));
  TPaveText* fit_tex_B = new TPaveText(0.5,0.41,0.8,0.61,"BRNDC");
  fit_tex_B->AddText("Gaussian parameters :");
  fit_tex_B->AddText(TString::Format("#mu = (%2.4f #pm %2.4f) GeV/c^{2}",gausmean_B_Mass.getVal(),gausmean_B_Mass.getError()));
  fit_tex_B->AddText(TString::Format("#sigma = (%2.4f #pm %2.4f) GeV/c^{2}",gaussig_B_Mass.getVal(),gaussig_B_Mass.getError()));
  fit_tex_B->AddText("");
  fit_tex_B->AddText(TString::Format("N_{sig} = %5.0f #pm %3.0f", nsig_B_Mass.getVal(), nsig_B_Mass.getError()));
  fit_tex_B->AddText(TString::Format("N_{bck} = %5.0f #pm %3.0f", nbck_B_Mass.getVal(), nbck_B_Mass.getError()));
  fit_tex_B->SetTextFont(43);
  fit_tex_B->SetFillColor(0);
  fit_tex_B->SetTextSize(TITLE_FONTSIZE - 6);

  // Plot
  RooPlot* frame_B = B_Mass.frame();
  dh_B_Mass.plotOn(frame_B);
  model_B_Mass.plotOn(frame_B,LineColor(9));
  model_B_Mass.plotOn(frame_B,Components(bck_B_Mass),LineColor(kBlue));
  model_B_Mass.plotOn(frame_B,Components(sig_B_Mass),LineColor(kRed));

  TCanvas* cn_fit_B_Mass = new TCanvas("cn_fit_B_Mass","cn_fit_B_Mass",800,800);
  frame_B->Draw();
  fit_tex_B->Draw("same");
  cms_style(); 
  cn_fit_B_Mass->SaveAs("Plots"+date+"/fit_B_Mass.C");
  cn_fit_B_Mass->SaveAs("Plots"+date+"/fit_B_Mass.pdf");
  cn_fit_B_Mass->SaveAs("Plots"+date+"/fit_B_Mass.eps");

  //~~~ D0+W mass ~~~
/*  
  RooRealVar B_WMass("B_WMass","D^{0}W^{#pm} mass",0.,10.,"GeV/c^{2}");
  RooDataHist dh_B_DiracMass("datahist_B_DiracMass","datahist_B_DiracMass",RooArgList(B_WMass),h_B_DiracMass,1.);
  RooDataHist dh_B_GausMass("datahist_B_GausMass","datahist_B_GausMass",RooArgList(B_WMass),h_B_GausMass,1.);

  // Plot
  RooPlot* frame_BDirac = B_WMass.frame();
  dh_B_DiracMass.plotOn(frame_BDirac);
  RooPlot* frame_BGaus = B_WMass.frame();
  dh_B_GausMass.plotOn(frame_BGaus);

  TCanvas* cn_fit_B_DiracMass = new TCanvas("cn_fit_B_DiracMass","cn_fit_B_DiracMass",800,800);
  frame_BDirac->Draw();
  cms_style(); 
  cn_fit_B_DiracMass->SaveAs("Plots"+date+"/fit_B_DiracMass.C");
  cn_fit_B_DiracMass->SaveAs("Plots"+date+"/fit_B_DiracMass.pdf");
  cn_fit_B_DiracMass->SaveAs("Plots"+date+"/fit_B_DiracMass.eps");
  TCanvas* cn_fit_B_GausMass = new TCanvas("cn_fit_B_GausMass","cn_fit_B_GausMass",800,800);
  frame_BGaus->Draw();
  cms_style(); 
  cn_fit_B_GausMass->SaveAs("Plots"+date+"/fit_B_GausMass.C");
  cn_fit_B_GausMass->SaveAs("Plots"+date+"/fit_B_GausMass.pdf");
  cn_fit_B_GausMass->SaveAs("Plots"+date+"/fit_B_GausMass.eps");
*/
  if (!inBatch) getchar();
  delete cn_fit_B_Mass; //delete cn_fit_B_DiracMass; delete cn_fit_B_GausMass;
  
  //=============================================================================================
  //                      Combinatorial reconstruction of the D0 -> K Pi
  //=============================================================================================

  double *SNR_D0combi_Mass = step(inBatch, fi, "D0combi_Mass", date);
  double *SNR_B_D0combiMass = step(inBatch, fi, "B_D0combiMass", date);
  std::cout << "D0 candidate : " << SNR_D0combi_Mass[0] << " +/- " << SNR_D0combi_Mass[1] << std::endl;
  std::cout << "D0 candidate when there is a non iso mu : " << SNR_B_D0combiMass[0] << " +/- " << SNR_B_D0combiMass[1] << std::endl;
  if (!inBatch) getchar();

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // have a look at other distributions 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  TH1D* h_D0combi_dRJet = (TH1D*)fi->Get("ana/h_D0combi_dRJet");
  h_D0combi_dRJet->GetXaxis()->SetRangeUser(0,0.5);
  h1_style(h_D0combi_dRJet, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "#DeltaR(D^{0},jet)");
  TCanvas* cn_D0combi_dRJet = new TCanvas("cn_D0combi_dRJet","cn_D0combi_dRJet",800,800);
  cn_D0combi_dRJet->cd();
  h_D0combi_dRJet->Draw("hist");
  cms_style();
  cn_D0combi_dRJet->SaveAs("Plots"+date+"/D0combi_dRJet.C");
  cn_D0combi_dRJet->SaveAs("Plots"+date+"/D0combi_dRJet.eps");
  cn_D0combi_dRJet->SaveAs("Plots"+date+"/D0combi_dRJet.pdf");

  TH1D* h_D0combi_p = (TH1D*)fi->Get("ana/h_D0combi_p");
  h_D0combi_p->GetXaxis()->SetRangeUser(12.,250.);
  h1_style(h_D0combi_p, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "p(D^{0}#rightarrow#kappa#pi) (GeV/c)");
  TCanvas* cn_D0combi_p = new TCanvas("cn_D0combi_p","cn_D0combi_p",800,800);
  cn_D0combi_p->cd();
  h_D0combi_p->Draw("hist");
  cms_style();
  cn_D0combi_p->SaveAs("Plots"+date+"/D0combi_p.C");
  cn_D0combi_p->SaveAs("Plots"+date+"/D0combi_p.eps");
  cn_D0combi_p->SaveAs("Plots"+date+"/D0combi_p.pdf");

  TH1D* h_D0combi_pT = (TH1D*)fi->Get("ana/h_D0combi_pT");
  h_D0combi_pT->GetXaxis()->SetRangeUser(12.,120.);
  h1_style(h_D0combi_pT, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "p_{T}(D^{0}#rightarrow#kappa#pi) (GeV/c)");
  TCanvas* cn_D0combi_pT = new TCanvas("cn_D0combi_pT","cn_D0combi_pT",800,800);
  cn_D0combi_pT->cd();
  h_D0combi_pT->Draw("hist");
  cms_style();
  cn_D0combi_pT->SaveAs("Plots"+date+"/D0combi_pT.C");
  cn_D0combi_pT->SaveAs("Plots"+date+"/D0combi_pT.eps");
  cn_D0combi_pT->SaveAs("Plots"+date+"/D0combi_pT.pdf");

  TH1D* h_D0combi_eta = (TH1D*)fi->Get("ana/h_D0combi_eta");
  h1_style(h_D0combi_eta, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "#eta(D^{0}#rightarrow#kappa#pi)");
  TCanvas* cn_D0combi_eta = new TCanvas("cn_D0combi_eta","cn_D0combi_eta",800,800);
  cn_D0combi_eta->cd();
  h_D0combi_eta->Draw("hist");
  cms_style();
  cn_D0combi_eta->SaveAs("Plots"+date+"/D0combi_eta.C");
  cn_D0combi_eta->SaveAs("Plots"+date+"/D0combi_eta.eps");
  cn_D0combi_eta->SaveAs("Plots"+date+"/D0combi_eta.pdf");

  TH1D* h_D0combi_phi = (TH1D*)fi->Get("ana/h_D0combi_phi");
  h1_style(h_D0combi_phi, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "#phi(D^{0}#rightarrow#kappa#pi)");
  TCanvas* cn_D0combi_phi = new TCanvas("cn_D0combi_phi","cn_D0combi_phi",800,800);
  cn_D0combi_phi->cd();
  h_D0combi_phi->Draw("hist");
  cms_style();
  cn_D0combi_phi->SaveAs("Plots"+date+"/D0combi_phi.C");
  cn_D0combi_phi->SaveAs("Plots"+date+"/D0combi_phi.eps");
  cn_D0combi_phi->SaveAs("Plots"+date+"/D0combi_phi.pdf");

  TH1D* h_D0combi_Mass = (TH1D*)fi->Get("ana/h_D0combi_Mass");
  h_D0combi_Mass->Rebin(2);
  h_D0combi_Mass->GetXaxis()->SetRangeUser(1.7,2.);
  h1_style(h_D0combi_Mass, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "m(D^{0}#rightarrow#kappa#pi) (GeV/c^{2})");
  TCanvas* cn_D0combi_Mass = new TCanvas("cn_D0combi_Mass","cn_D0combi_Mass",800,800);
  cn_D0combi_Mass->cd();
  h_D0combi_Mass->Draw("hist");
  cms_style();
  cn_D0combi_Mass->SaveAs("Plots"+date+"/D0combi_Mass.C");
  cn_D0combi_Mass->SaveAs("Plots"+date+"/D0combi_Mass.eps");
  cn_D0combi_Mass->SaveAs("Plots"+date+"/D0combi_Mass.pdf");

  TH1D* h_B_D0combiMass = (TH1D*)fi->Get("ana/h_B_D0combiMass");
  h_B_D0combiMass->GetXaxis()->SetRangeUser(1.7,2.);
  h1_style(h_B_D0combiMass, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "m(D^{0}#rightarrow#kappa#pi) (GeV/c^{2})");
  TCanvas* cn_B_D0combiMass = new TCanvas("cn_B_D0combiMass","cn_B_D0combiMass",800,800);
  cn_B_D0combiMass->cd();
  h_B_D0combiMass->Draw("hist");
  cms_style();
  cn_B_D0combiMass->SaveAs("Plots"+date+"/B_D0combiMass.C");
  cn_B_D0combiMass->SaveAs("Plots"+date+"/B_D0combiMass.eps");
  cn_B_D0combiMass->SaveAs("Plots"+date+"/B_D0combiMass.pdf");

  TH1D* h_B_DeltaRD0combiMu = (TH1D*)fi->Get("ana/h_B_DeltaRD0combiMu");
  h1_style(h_B_DeltaRD0combiMu, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "#DeltaR(D^{0},#mu^{#pm})");
  TCanvas* cn_B_DeltaRD0combiMu = new TCanvas("cn_B_DeltaRD0combiMu","cn_B_DeltaRD0combiMu",800,800);
  cn_B_DeltaRD0combiMu->cd();
  h_B_DeltaRD0combiMu->Draw("hist");
  cms_style();
  cn_B_DeltaRD0combiMu->SaveAs("Plots"+date+"/B_DeltaRD0combiMu.C");
  cn_B_DeltaRD0combiMu->SaveAs("Plots"+date+"/B_DeltaRD0combiMu.eps");
  cn_B_DeltaRD0combiMu->SaveAs("Plots"+date+"/B_DeltaRD0combiMu.pdf");

  TH1D* h_B_pD0combipMu = (TH1D*)fi->Get("ana/h_B_pD0combipMu");
  h_B_pD0combipMu->Rebin(2);
  h_B_pD0combipMu->GetXaxis()->SetRangeUser(0.1,10.);
  h1_style(h_B_pD0combipMu, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "p(D^{0})/p(#mu^{#pm})");
  TCanvas* cn_B_pD0combipMu = new TCanvas("cn_B_pD0combipMu","cn_B_pD0combipMu",800,800);
  cn_B_pD0combipMu->cd();
  cn_B_pD0combipMu->SetLogx();
  h_B_pD0combipMu->Draw("hist");
  cms_style();
  cn_B_pD0combipMu->SaveAs("Plots"+date+"/B_pD0combipMu.C");
  cn_B_pD0combipMu->SaveAs("Plots"+date+"/B_pD0combipMu.eps");
  cn_B_pD0combipMu->SaveAs("Plots"+date+"/B_pD0combipMu.pdf");

  TH1D* h_B_pTD0combipTMu = (TH1D*)fi->Get("ana/h_B_pTD0combipTMu");
  h_B_pTD0combipTMu->Rebin(2);
  h_B_pTD0combipTMu->GetXaxis()->SetRangeUser(0.1,10.);
  h1_style(h_B_pTD0combipTMu, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "p_{T}(D^{0})/p_{T}(#mu^{#pm})");
  TCanvas* cn_B_pTD0combipTMu = new TCanvas("cn_B_pTD0combipTMu","cn_B_pTD0combipTMu",800,800);
  cn_B_pTD0combipTMu->cd();
  cn_B_pTD0combipTMu->SetLogx();
  h_B_pTD0combipTMu->Draw("hist");
  cms_style();
  cn_B_pTD0combipTMu->SaveAs("Plots"+date+"/B_pTD0combipTMu.C");
  cn_B_pTD0combipTMu->SaveAs("Plots"+date+"/B_pTD0combipTMu.eps");
  cn_B_pTD0combipTMu->SaveAs("Plots"+date+"/B_pTD0combipTMu.pdf");

  TH1D* h_B_combiMass = (TH1D*)fi->Get("ana/h_B_combiMass");
  h_B_combiMass->Rebin(20);
  h_B_combiMass->GetXaxis()->SetRangeUser(1.5,6.);
  h1_style(h_B_combiMass, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "m(D^{0}+#mu^{#pm}) (GeV/c^{2})");
  TCanvas* cn_B_combiMass = new TCanvas("cn_B_combiMass","cn_B_combiMass",800,800);
  cn_B_combiMass->cd();
  h_B_combiMass->Draw("hist");
  cms_style();
  cn_B_combiMass->SaveAs("Plots"+date+"/B_combiMass.C");
  cn_B_combiMass->SaveAs("Plots"+date+"/B_combiMass.eps");
  cn_B_combiMass->SaveAs("Plots"+date+"/B_combiMass.pdf");

  if (!inBatch) getchar();
  delete cn_B_D0combiMass; delete cn_B_combiMass; 
  delete cn_D0combi_dRJet; delete cn_D0combi_p; delete cn_D0combi_pT; delete cn_D0combi_eta; delete cn_D0combi_phi;  delete h_B_DeltaRD0combiMu; delete cn_B_pD0combipMu; delete cn_B_pTD0combipTMu;
    
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Fit Bpm mass distributions 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  //~~~ D0+mu mass ~~~
  
  RooRealVar B_combiMass("B_combiMass","D^{0}#mu^{#pm} mass",1.5,6.,"GeV/c^{2}");
  RooDataHist dh_B_combiMass("datahist_B_combiMass","datahist_B_combiMass",RooArgList(B_combiMass),h_B_combiMass,1.);
  
  // Signal+Background pdf
  //RooRealVar gausmean_B_combiMass("gausmean_B_combiMass","gausmean_B_combiMass",2.6,2.55,2.65); //FIXME
  RooRealVar gausmean_B_combiMass("gausmean_B_combiMass","gausmean_B_combiMass",3.4,3.2,3.6);
  RooRealVar gaussig_B_combiMass("gaussig_B_combiMass","gaussig_B_combiMass",0.5,0.1,1.);
  //RooRealVar landmean_B_combiMass("landmean_B_combiMass","landmean_B_combiMass",2.4,2.2,2.6); //FIXME
  //RooRealVar landsig_B_combiMass("landsig_B_combiMass","landsig_B_combiMass",0.5,0.,1.); //FIXME
  //RooLandau bck_B_combiMass("bck_B_combiMass", "landau background", B_combiMass, landmean_B_combiMass, landsig_B_combiMass); //FIXME
  RooRealVar alpha_bck2("alpha_bck2","alpha_bck2",-1.,-1.5,-0.5);
  RooExponential bck_B_combiMass("bck_B_combiMass","exp background",B_combiMass,alpha_bck2);
  RooGaussian sig_B_combiMass("sig_B_combiMass","gaussian signal",B_combiMass,gausmean_B_combiMass,gaussig_B_combiMass);
  RooRealVar nsig_B_combiMass("nsig_B_combiMass","number of signal events",0.1*h_B_combiMass->Integral(),0.,h_B_combiMass->Integral());
  RooRealVar nbck_B_combiMass("nbck_B_combiMass","number of background events",0.9*h_B_combiMass->Integral(),0.,h_B_combiMass->Integral());
  RooAddPdf model_B_combiMass("model_B_combiMass","model_B_combiMass",RooArgList(sig_B_combiMass,bck_B_combiMass),RooArgList(nsig_B_combiMass,nbck_B_combiMass)) ;

  // Fit
  model_B_combiMass.fitTo(dh_B_combiMass,Range(2.2,5.5));
  TPaveText* fit_tex_Bcombi = new TPaveText(0.65,0.61,0.95,0.81,"BRNDC");
  fit_tex_Bcombi->AddText("Gaussian parameters :");
  fit_tex_Bcombi->AddText(TString::Format("#mu = (%2.4f #pm %2.4f) GeV/c^{2}",gausmean_B_combiMass.getVal(),gausmean_B_combiMass.getError()));
  fit_tex_Bcombi->AddText(TString::Format("#sigma = (%2.4f #pm %2.4f) GeV/c^{2}",gaussig_B_combiMass.getVal(),gaussig_B_combiMass.getError()));
  fit_tex_Bcombi->AddText("");
  fit_tex_Bcombi->AddText(TString::Format("N_{sig} = %5.0f #pm %3.0f", nsig_B_combiMass.getVal(), nsig_B_combiMass.getError()));
  fit_tex_Bcombi->AddText(TString::Format("N_{bck} = %5.0f #pm %3.0f", nbck_B_combiMass.getVal(), nbck_B_combiMass.getError()));
  fit_tex_Bcombi->SetTextFont(43);
  fit_tex_Bcombi->SetFillColor(0);
  fit_tex_Bcombi->SetTextSize(TITLE_FONTSIZE - 6);

  // Plot
  RooPlot* frame_Bcombi = B_combiMass.frame();
  dh_B_combiMass.plotOn(frame_Bcombi); 
  model_B_combiMass.plotOn(frame_Bcombi,LineColor(9));
  model_B_combiMass.plotOn(frame_Bcombi,Components(bck_B_combiMass),LineColor(kBlue));
  model_B_combiMass.plotOn(frame_Bcombi,Components(sig_B_combiMass),LineColor(kRed));

  TCanvas* cn_fit_B_combiMass = new TCanvas("cn_fit_B_combiMass","cn_fit_B_combiMass",800,800);
  frame_Bcombi->Draw();
  fit_tex_Bcombi->Draw("same");
  cms_style(); 
  cn_fit_B_combiMass->SaveAs("Plots"+date+"/fit_B_combiMass.C");
  cn_fit_B_combiMass->SaveAs("Plots"+date+"/fit_B_combiMass.pdf");
  cn_fit_B_combiMass->SaveAs("Plots"+date+"/fit_B_combiMass.eps");

  if (!inBatch) getchar();
  delete cn_fit_B_combiMass;

  //=============================================================================================
  //                  Constrained Kalman Vertex Fitter for the D0 -> K Pi reconstruction
  //=============================================================================================

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // optimize chi2 cut on D0 -> K Pi reconstruction
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  double *SNR_D0cons_Mass = BRstep(inBatch, fi, "D0cons_Mass", date);
  double *SNR_B_D0consMass = BRstep(inBatch, fi, "B_D0consMass", date);

  // Print chi2 otpimization results 
  std::cout << "\n\nD0 candidate mass SNR for \n" << std::endl;
  std::cout << "... ctau/sigma(ctau) > 50 : " << SNR_D0cons_Mass[0] << " +/- " << SNR_D0cons_Mass[1] << std::endl;
  std::cout << "... and non iso mu : " << SNR_B_D0consMass[0] << " +/- " << SNR_B_D0consMass[1] << std::endl;

  if (!inBatch) getchar();

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // have a look at other distributions 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  TH1D* h_D0cons_dRJet = (TH1D*)fi->Get("ana/h_D0cons_dRJet");
  h_D0cons_dRJet->GetXaxis()->SetRangeUser(0,0.5);
  h1_style(h_D0cons_dRJet, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "#DeltaR(D^{0},jet)");
  TCanvas* cn_D0cons_dRJet = new TCanvas("cn_D0cons_dRJet","cn_D0cons_dRJet",800,800);
  cn_D0cons_dRJet->cd();
  h_D0cons_dRJet->Draw("hist");
  cms_style();
  cn_D0cons_dRJet->SaveAs("Plots"+date+"/D0cons_dRJet.C");
  cn_D0cons_dRJet->SaveAs("Plots"+date+"/D0cons_dRJet.eps");
  cn_D0cons_dRJet->SaveAs("Plots"+date+"/D0cons_dRJet.pdf");

  TH1D* h_D0cons_Chi2NDOF = (TH1D*)fi->Get("ana/h_D0cons_Chi2NDOF");
  h_D0cons_Chi2NDOF->GetXaxis()->SetRangeUser(0,0.25);
  h1_style(h_D0cons_Chi2NDOF, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "#chi^{2}/NDOF");
  TCanvas* cn_D0cons_Chi2NDOF = new TCanvas("cn_D0cons_Chi2NDOF","cn_D0cons_Chi2NDOF",800,800);
  cn_D0cons_Chi2NDOF->cd();
  h_D0cons_Chi2NDOF->Draw("hist");
  cms_style();
  cn_D0cons_Chi2NDOF->SaveAs("Plots"+date+"/D0cons_Chi2NDOF.C");
  cn_D0cons_Chi2NDOF->SaveAs("Plots"+date+"/D0cons_Chi2NDOF.eps");
  cn_D0cons_Chi2NDOF->SaveAs("Plots"+date+"/D0cons_Chi2NDOF.pdf");

  TH1D* h_D0cons_L = (TH1D*)fi->Get("ana/h_D0cons_L");
  h_D0cons_L->GetXaxis()->SetRangeUser(0,0.25);
  h1_style(h_D0cons_L, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "c#tau (D^{0}#rightarrow#kappa#pi) (cm)");
  TCanvas* cn_D0cons_L = new TCanvas("cn_D0cons_L","cn_D0cons_L",800,800);
  cn_D0cons_L->cd();
  h_D0cons_L->Draw("hist");
  cms_style();
  cn_D0cons_L->SaveAs("Plots"+date+"/D0cons_L.C");
  cn_D0cons_L->SaveAs("Plots"+date+"/D0cons_L.eps");
  cn_D0cons_L->SaveAs("Plots"+date+"/D0cons_L.pdf");

  TH1D* h_D0cons_SigmaL = (TH1D*)fi->Get("ana/h_D0cons_SigmaL");
  h_D0cons_SigmaL->GetXaxis()->SetRangeUser(0,0.001);
  h1_style(h_D0cons_SigmaL, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "#Delta(c#tau) (D^{0}#rightarrow#kappa#pi) (cm)");
  TCanvas* cn_D0cons_SigmaL = new TCanvas("cn_D0cons_SigmaL","cn_D0cons_SigmaL",800,800);
  cn_D0cons_SigmaL->cd();
  h_D0cons_SigmaL->Draw("hist");
  cms_style();
  cn_D0cons_SigmaL->SaveAs("Plots"+date+"/D0cons_SigmaL.C");
  cn_D0cons_SigmaL->SaveAs("Plots"+date+"/D0cons_SigmaL.eps");
  cn_D0cons_SigmaL->SaveAs("Plots"+date+"/D0cons_SigmaL.pdf");

  TH1D* h_D0cons_p = (TH1D*)fi->Get("ana/h_D0cons_p");
  h_D0cons_p->Rebin(2);
  h_D0cons_p->GetXaxis()->SetRangeUser(12.,250.);
  h1_style(h_D0cons_p, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "p(D^{0}#rightarrow#kappa#pi) (GeV/c)");
  TCanvas* cn_D0cons_p = new TCanvas("cn_D0cons_p","cn_D0cons_p",800,800);
  cn_D0cons_p->cd();
  h_D0cons_p->Draw("hist");
  cms_style();
  cn_D0cons_p->SaveAs("Plots"+date+"/D0cons_p.C");
  cn_D0cons_p->SaveAs("Plots"+date+"/D0cons_p.eps");
  cn_D0cons_p->SaveAs("Plots"+date+"/D0cons_p.pdf");

  TH1D* h_D0cons_pT = (TH1D*)fi->Get("ana/h_D0cons_pT");
  h_D0cons_pT->GetXaxis()->SetRangeUser(12.,120.);
  h1_style(h_D0cons_pT, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "p_{T}(D^{0}#rightarrow#kappa#pi) (GeV/c)");
  TCanvas* cn_D0cons_pT = new TCanvas("cn_D0cons_pT","cn_D0cons_pT",800,800);
  cn_D0cons_pT->cd();
  h_D0cons_pT->Draw("hist");
  cms_style();
  cn_D0cons_pT->SaveAs("Plots"+date+"/D0cons_pT.C");
  cn_D0cons_pT->SaveAs("Plots"+date+"/D0cons_pT.eps");
  cn_D0cons_pT->SaveAs("Plots"+date+"/D0cons_pT.pdf");

  TH1D* h_D0cons_eta = (TH1D*)fi->Get("ana/h_D0cons_eta");
  h1_style(h_D0cons_eta, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "#eta(D^{0}#rightarrow#kappa#pi)");
  TCanvas* cn_D0cons_eta = new TCanvas("cn_D0cons_eta","cn_D0cons_eta",800,800);
  cn_D0cons_eta->cd();
  h_D0cons_eta->Draw("hist");
  cms_style();
  cn_D0cons_eta->SaveAs("Plots"+date+"/D0cons_eta.C");
  cn_D0cons_eta->SaveAs("Plots"+date+"/D0cons_eta.eps");
  cn_D0cons_eta->SaveAs("Plots"+date+"/D0cons_eta.pdf");

  TH1D* h_D0cons_phi = (TH1D*)fi->Get("ana/h_D0cons_phi");
  h1_style(h_D0cons_phi, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "#phi(D^{0}#rightarrow#kappa#pi)");
  TCanvas* cn_D0cons_phi = new TCanvas("cn_D0cons_phi","cn_D0cons_phi",800,800);
  cn_D0cons_phi->cd();
  h_D0cons_phi->Draw("hist");
  cms_style();
  cn_D0cons_phi->SaveAs("Plots"+date+"/D0cons_phi.C");
  cn_D0cons_phi->SaveAs("Plots"+date+"/D0cons_phi.eps");
  cn_D0cons_phi->SaveAs("Plots"+date+"/D0cons_phi.pdf");

  TH1D* h_D0cons_Mass = (TH1D*)fi->Get("ana/h_D0cons_Mass");
  h_D0cons_Mass->GetXaxis()->SetRangeUser(1.7,2.);
  h1_style(h_D0cons_Mass, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "m(D^{0}#rightarrow#kappa#pi) (GeV/c^{2})");
  TCanvas* cn_D0cons_Mass = new TCanvas("cn_D0cons_Mass","cn_D0cons_Mass",800,800);
  cn_D0cons_Mass->cd();
  h_D0cons_Mass->Draw("hist");
  cms_style();
  cn_D0cons_Mass->SaveAs("Plots"+date+"/D0cons_Mass.C");
  cn_D0cons_Mass->SaveAs("Plots"+date+"/D0cons_Mass.eps");
  cn_D0cons_Mass->SaveAs("Plots"+date+"/D0cons_Mass.pdf");

  TH1D* h_B_D0consMass = (TH1D*)fi->Get("ana/h_B_D0consMass");
  h_B_D0consMass->GetXaxis()->SetRangeUser(1.7,2.);
  h1_style(h_B_D0consMass, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "m(D^{0}#rightarrow#kappa#pi) (GeV/c^{2})");
  TCanvas* cn_B_D0consMass = new TCanvas("cn_B_D0consMass","cn_B_D0consMass",800,800);
  cn_B_D0consMass->cd();
  h_B_D0consMass->Draw("hist");
  cms_style();
  cn_B_D0consMass->SaveAs("Plots"+date+"/B_D0consMass.C");
  cn_B_D0consMass->SaveAs("Plots"+date+"/B_D0consMass.eps");
  cn_B_D0consMass->SaveAs("Plots"+date+"/B_D0consMass.pdf");

  TH1D* h_B_DeltaRD0consMu = (TH1D*)fi->Get("ana/h_B_DeltaRD0consMu");
  h1_style(h_B_DeltaRD0consMu, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "#DeltaR(D^{0},#mu^{#pm})");
  TCanvas* cn_B_DeltaRD0consMu = new TCanvas("cn_B_DeltaRD0consMu","cn_B_DeltaRD0consMu",800,800);
  cn_B_DeltaRD0consMu->cd();
  h_B_DeltaRD0consMu->Draw("hist");
  cms_style();
  cn_B_DeltaRD0consMu->SaveAs("Plots"+date+"/B_DeltaRD0consMu.C");
  cn_B_DeltaRD0consMu->SaveAs("Plots"+date+"/B_DeltaRD0consMu.eps");
  cn_B_DeltaRD0consMu->SaveAs("Plots"+date+"/B_DeltaRD0consMu.pdf");

  TH1D* h_B_pD0conspMu = (TH1D*)fi->Get("ana/h_B_pD0conspMu");
  h_B_pD0conspMu->Rebin(5);
  h_B_pD0conspMu->GetXaxis()->SetRangeUser(0.1,10.);
  h1_style(h_B_pD0conspMu, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "p(D^{0})/p(#mu^{#pm})");
  TCanvas* cn_B_pD0conspMu = new TCanvas("cn_B_pD0conspMu","cn_B_pD0conspMu",800,800);
  cn_B_pD0conspMu->cd();
  cn_B_pD0conspMu->SetLogx();
  h_B_pD0conspMu->Draw("hist");
  cms_style();
  cn_B_pD0conspMu->SaveAs("Plots"+date+"/B_pD0conspMu.C");
  cn_B_pD0conspMu->SaveAs("Plots"+date+"/B_pD0conspMu.eps");
  cn_B_pD0conspMu->SaveAs("Plots"+date+"/B_pD0conspMu.pdf");

  TH1D* h_B_pTD0conspTMu = (TH1D*)fi->Get("ana/h_B_pTD0conspTMu");
  h_B_pTD0conspTMu->Rebin(5);
  h_B_pTD0conspTMu->GetXaxis()->SetRangeUser(0.1,10.);
  h1_style(h_B_pTD0conspTMu, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "p_{T}(D^{0})/p_{T}(#mu^{#pm})");
  TCanvas* cn_B_pTD0conspTMu = new TCanvas("cn_B_pTD0conspTMu","cn_B_pTD0conspTMu",800,800);
  cn_B_pTD0conspTMu->cd();
  cn_B_pTD0conspTMu->SetLogx();
  h_B_pTD0conspTMu->Draw("hist");
  cms_style();
  cn_B_pTD0conspTMu->SaveAs("Plots"+date+"/B_pTD0conspTMu.C");
  cn_B_pTD0conspTMu->SaveAs("Plots"+date+"/B_pTD0conspTMu.eps");
  cn_B_pTD0conspTMu->SaveAs("Plots"+date+"/B_pTD0conspTMu.pdf");

  TH1D* h_B_consMass = (TH1D*)fi->Get("ana/h_B_consMass");
  h_B_consMass->Rebin(5);
  h_B_consMass->GetXaxis()->SetRangeUser(1.5,6.);
  h1_style(h_B_consMass, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "m(D^{0}+#mu^{#pm}) (GeV/c^{2})");
  TCanvas* cn_B_consMass = new TCanvas("cn_B_consMass","cn_B_consMass",800,800);
  cn_B_consMass->cd();
  h_B_consMass->Draw("hist");
  cms_style();
  cn_B_consMass->SaveAs("Plots"+date+"/B_consMass.C");
  cn_B_consMass->SaveAs("Plots"+date+"/B_consMass.eps");
  cn_B_consMass->SaveAs("Plots"+date+"/B_consMass.pdf");

  if (!inBatch) getchar();
  delete cn_D0cons_Chi2NDOF; delete cn_D0cons_L; delete cn_D0cons_SigmaL; delete cn_D0cons_Mass; 
  delete cn_D0cons_dRJet; delete cn_D0cons_p; delete cn_D0cons_pT; delete cn_D0cons_eta; delete cn_D0cons_phi;  delete h_B_DeltaRD0consMu; delete cn_B_pD0conspMu; delete cn_B_pTD0conspTMu;
  delete cn_B_D0consMass; delete cn_B_consMass; 
    
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Fit Bpm mass distributions 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  //~~~ D0+mu mass ~~~
  
  RooRealVar B_consMass("B_consMass","D^{0}#mu^{#pm} mass",1.5,6.,"GeV/c^{2}");
  RooDataHist dh_B_consMass("datahist_B_consMass","datahist_B_consMass",RooArgList(B_consMass),h_B_consMass,1.);
  
  // Signal+Background pdf
  //RooRealVar gausmean_B_consMass("gausmean_B_consMass","gausmean_B_consMass",2.6,2.55,2.65); //FIXME
  RooRealVar gausmean_B_consMass("gausmean_B_consMass","gausmean_B_consMass",3.75,2.5,5.);
  RooRealVar gaussig_B_consMass("gaussig_B_consMass","gaussig_B_consMass",0.0193,0.,0.05);
  //RooRealVar landmean_B_consMass("landmean_B_consMass","landmean_B_consMass",2.4,2.2,2.6); //FIXME
  //RooRealVar landsig_B_consMass("landsig_B_consMass","landsig_B_consMass",0.5,0.,1.); //FIXME
  //RooLandau bck_B_consMass("bck_B_consMass", "landau background", B_consMass, landmean_B_consMass, landsig_B_consMass); //FIXME
  RooRealVar alpha_consbck("alpha_consbck","alpha_consbck",-1.5,-2.5,-0.5);
  RooExponential bck_B_consMass("bck_B_consMass","exp background",B_consMass,alpha_consbck);
  RooGaussian sig_B_consMass("sig_B_consMass","gaussian signal",B_consMass,gausmean_B_consMass,gaussig_B_consMass);
  RooRealVar nsig_B_consMass("nsig_B_consMass","number of signal events",0.1*h_B_consMass->Integral(),0.,h_B_consMass->Integral());
  RooRealVar nbck_B_consMass("nbck_B_consMass","number of background events",0.9*h_B_consMass->Integral(),0.,h_B_consMass->Integral());
  RooAddPdf model_B_consMass("model_B_consMass","model_B_consMass",RooArgList(sig_B_consMass,bck_B_consMass),RooArgList(nsig_B_consMass,nbck_B_consMass)) ;

  // Fit
  model_B_consMass.fitTo(dh_B_consMass,Range(2.8,4.8));
  TPaveText* fit_tex_Bcons = new TPaveText(0.5,0.41,0.8,0.61,"BRNDC");
  fit_tex_Bcons->AddText("Gaussian parameters :");
  fit_tex_Bcons->AddText(TString::Format("#mu = (%2.4f #pm %2.4f) GeV/c^{2}",gausmean_B_consMass.getVal(),gausmean_B_consMass.getError()));
  fit_tex_Bcons->AddText(TString::Format("#sigma = (%2.4f #pm %2.4f) GeV/c^{2}",gaussig_B_consMass.getVal(),gaussig_B_consMass.getError()));
  fit_tex_Bcons->AddText("");
  fit_tex_Bcons->AddText(TString::Format("N_{sig} = %5.0f #pm %3.0f", nsig_B_consMass.getVal(), nsig_B_consMass.getError()));
  fit_tex_Bcons->AddText(TString::Format("N_{bck} = %5.0f #pm %3.0f", nbck_B_consMass.getVal(), nbck_B_consMass.getError()));
  fit_tex_Bcons->SetTextFont(43);
  fit_tex_Bcons->SetFillColor(0);
  fit_tex_Bcons->SetTextSize(TITLE_FONTSIZE - 6);

  // Plot
  RooPlot* frame_Bcons = B_consMass.frame();
  dh_B_consMass.plotOn(frame_Bcons);
  model_B_consMass.plotOn(frame_Bcons,LineColor(9));
  model_B_consMass.plotOn(frame_Bcons,Components(bck_B_consMass),LineColor(kBlue));
  model_B_consMass.plotOn(frame_Bcons,Components(sig_B_consMass),LineColor(kRed));

  TCanvas* cn_fit_B_consMass = new TCanvas("cn_fit_B_consMass","cn_fit_B_consMass",800,800);
  frame_Bcons->Draw();
  fit_tex_Bcons->Draw("same");
  cms_style(); 
  cn_fit_B_consMass->SaveAs("Plots"+date+"/fit_B_consMass.C");
  cn_fit_B_consMass->SaveAs("Plots"+date+"/fit_B_consMass.pdf");
  cn_fit_B_consMass->SaveAs("Plots"+date+"/fit_B_consMass.eps");

  if (!inBatch) getchar();
  delete cn_fit_B_consMass; 

  //close
  fi->Close();
  return 0;
}

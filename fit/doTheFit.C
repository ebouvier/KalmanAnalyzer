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
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooGaussian.h"
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
int doTheFit(bool inBatch = true)
//---------------------------------------------------------------
{
  TStyle* m_style = createStyle();
  m_style->cd();
  using namespace RooFit;
  if (inBatch) gROOT->SetBatch(true);

  TFile *fi = TFile::Open("../test/kalmanAnalyzed.root");
  TH1D* histo = (TH1D*)fi->Get("ana/h_D0_Mass");

  RooRealVar x("mass","D^{0} mass",1.7,2.,"GeV/c^{2}");
  RooDataHist dh("datahist","datahist",RooArgList(x),histo,1.);
  
  // Signal+Background pdf
  RooRealVar mean("mean","mean",1.86,1.8,2.0);
  RooRealVar sigma("sigma","sigma",0.005,0.001,0.200);
  RooRealVar lambda("lambda", "slope", -0.1, -5., 0.);
  RooGaussian sig("sig","gaussian signal",x,mean,sigma);
  RooExponential bck("bck", "exponential background", x, lambda);
  RooRealVar fsig("fsig","signal fraction 1",0.10,0.,1.);
  RooAddPdf model("model","model",RooArgList(sig,bck),RooArgList(fsig)) ;

  // Fit
  model.fitTo(dh);
  double mass     = mean.getVal();
  double masserr  = mean.getError();
  double width    = sigma.getVal();
  double widtherr = sigma.getError();
  std::cout << "\nSignal mean = " << mass << " +/- " << masserr << std::endl;
  std::cout << "Signal width = " << width << " +/- " << widtherr << "\n" << std::endl;
  TPaveText* fit_tex = new TPaveText(0.22,0.47,0.52,0.6,"BRNDC");
  fit_tex->AddText("Gaussian parameters :");
  fit_tex->AddText(TString::Format("#mu = (%2.4f #pm %2.4f) GeV/c^{2}",mass,masserr));
  fit_tex->AddText(TString::Format("#sigma = (%2.4f #pm %2.4f) GeV/c^{2}",width,widtherr));
  fit_tex->SetTextFont(43);
  fit_tex->SetFillColor(0);
  fit_tex->SetTextSize(TITLE_FONTSIZE - 6);

  // Plot
  RooPlot* frame = x.frame();
  dh.plotOn(frame);
  model.plotOn(frame,LineColor(9));
  model.plotOn(frame,Components(bck),LineColor(kBlue));
  model.plotOn(frame,Components(RooArgSet(sig)),LineColor(kRed));

  TCanvas* cn = new TCanvas("cn","D0 meson candidate",800,800);
  frame->Draw();
  fit_tex->Draw("same");
  cms_style(); 
  cn->SaveAs("D0_Mass.png");
  cn->SaveAs("D0_Mass.pdf");
  cn->SaveAs("D0_Mass.eps");
  if (!inBatch) getchar();

  fi->Close();
  return 0;
}


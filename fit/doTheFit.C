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
  RooRealVar sigma("sigma","sigma",0.0193,0.,0.05);
  RooRealVar lambda("lambda", "slope", -0.1, -5., 0.);
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

  TPaveText* fit_tex = new TPaveText(0.22,0.36,0.52,0.61,"BRNDC");
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
  fit_tex->Draw("same");
  cms_style(); 
  cn->SaveAs("Plots"+date+"/"+name+".png");
  cn->SaveAs("Plots"+date+"/"+name+".pdf");
  cn->SaveAs("Plots"+date+"/"+name+".eps");
  if (!inBatch) getchar();

  delete cn; delete frame; delete fit_tex; delete histo;
  return SNR;
}

//---------------------------------------------------------------
int doTheFit(bool inBatch = true, TString date = "")
//---------------------------------------------------------------
{
  TStyle* m_style = createStyle();
  m_style->cd();
  if (inBatch) gROOT->SetBatch(true);
  gROOT->ProcessLine(".! mkdir Plots"+date);

  TFile *fi = TFile::Open("../test/kalmanAnalyzed.root"); 
  const unsigned int numberOfPoints = 7;
  double *SNR_D0Cand_MassChi2Inf1 = step(inBatch, fi, "D0Cand_MassChi2Inf1", date); 
  double *SNR_D0Cand_MassChi2Inf1p5 = step(inBatch, fi, "D0Cand_MassChi2Inf1p5", date); 
  double *SNR_D0Cand_MassChi2Inf2 = step(inBatch, fi, "D0Cand_MassChi2Inf2", date);
  double *SNR_D0Cand_MassChi2Inf2p5 = step(inBatch, fi, "D0Cand_MassChi2Inf2p5", date); 
  double *SNR_D0Cand_MassChi2Inf3 = step(inBatch, fi, "D0Cand_MassChi2Inf3", date); 
  double *SNR_D0Cand_MassChi2Inf3p5 = step(inBatch, fi, "D0Cand_MassChi2Inf3p5", date); 
  double *SNR_D0Cand_MassChi2Inf4 = step(inBatch, fi, "D0Cand_MassChi2Inf4", date); 

  double x[numberOfPoints] = {1., 1.5, 2., 2.5, 3., 3.5, 4.};
  double ex[numberOfPoints] = {0., 0., 0., 0., 0., 0., 0.};
  double y[numberOfPoints] = {SNR_D0Cand_MassChi2Inf1[0], SNR_D0Cand_MassChi2Inf1p5[0], SNR_D0Cand_MassChi2Inf2[0], SNR_D0Cand_MassChi2Inf2p5[0], SNR_D0Cand_MassChi2Inf3[0], SNR_D0Cand_MassChi2Inf3p5[0], SNR_D0Cand_MassChi2Inf4[0]};
  double ey[numberOfPoints] = {SNR_D0Cand_MassChi2Inf1[1], SNR_D0Cand_MassChi2Inf1p5[1], SNR_D0Cand_MassChi2Inf2[1], SNR_D0Cand_MassChi2Inf2p5[1], SNR_D0Cand_MassChi2Inf3[1], SNR_D0Cand_MassChi2Inf3p5[1], SNR_D0Cand_MassChi2Inf4[1]};

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

  std::cout << "\n\nD0 candidate mass SNR for \n" << std::endl;
  std::cout << "Chi2 < 1   : " <<  SNR_D0Cand_MassChi2Inf1[0] << " +/- " << SNR_D0Cand_MassChi2Inf1[1] << std::endl;
  std::cout << "Chi2 < 1.5 : " <<  SNR_D0Cand_MassChi2Inf1p5[0] << " +/- " << SNR_D0Cand_MassChi2Inf1p5[1] << std::endl;
  std::cout << "Chi2 < 2   : " <<  SNR_D0Cand_MassChi2Inf2[0] << " +/- " << SNR_D0Cand_MassChi2Inf2[1] << std::endl;
  std::cout << "Chi2 < 2.5 : " <<  SNR_D0Cand_MassChi2Inf2p5[0] << " +/- " << SNR_D0Cand_MassChi2Inf2p5[1] << std::endl;
  std::cout << "Chi2 < 3   : " <<  SNR_D0Cand_MassChi2Inf3[0] << " +/- " << SNR_D0Cand_MassChi2Inf3[1] << std::endl;
  std::cout << "Chi2 < 3.5 : " <<  SNR_D0Cand_MassChi2Inf3p5[0] << " +/- " << SNR_D0Cand_MassChi2Inf3p5[1] << std::endl;
  std::cout << "Chi2 < 4   : " <<  SNR_D0Cand_MassChi2Inf4[0] << " +/- " << SNR_D0Cand_MassChi2Inf4[1] << "\n" << std::endl;
  std::cout << "... and ctau/sigma(ctau) > 50 : " << SNR_D0_Mass[0] << " +/- " << SNR_D0_Mass[1] << std::endl;

  fi->Close();
  return 0;
}

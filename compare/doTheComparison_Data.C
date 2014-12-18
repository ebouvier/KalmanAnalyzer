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
#include "RooCBShape.h"
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
  style->SetStatH(0.12);
  style->SetStatW(0.3);
  style->SetStatY(0.92);
  style->SetStatX(0.94);

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

void cms_style(bool isData = true){
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

  TPaveText* fit_tex = new TPaveText(0.7,0.8,0.91,0.91,"BRNDC");
//  TPaveText* fit_tex = new TPaveText(0.21,0.21,0.51,0.46,"BRNDC");
//  fit_tex->AddText("Gaussian parameters :");
//  fit_tex->AddText(TString::Format("#mu = (%2.4f #pm %2.4f) GeV/c^{2}",mass,masserr));
//  fit_tex->AddText(TString::Format("#sigma = (%2.4f #pm %2.4f) GeV/c^{2}",width,widtherr));
//  fit_tex->AddText("");
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
void doTheComparison_step(bool inBatch = true, TString date = "", TString file = "")
//---------------------------------------------------------------
{
  using namespace RooFit;
  TStyle* m_style = createStyle();
  m_style->cd();
  if (inBatch) gROOT->SetBatch(true);
  if (date.Length() > 0) date = "_" + date;
  gROOT->ProcessLine(".! mkdir Plots"+date);

  TFile *fi = TFile::Open(file); 

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // B_D0 mass fit signal+background
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  double *SNR_B_D0Mass = step(inBatch, fi, "B_D0Mass", date); 
  double *SNR_B_D0consMass = step(inBatch, fi, "B_D0consMass", date); 
  double *SNR_B_D0combiMass = step(inBatch, fi, "B_D0combiMass", date); 
  double *SNR_B_D0optcombiMass = step(inBatch, fi, "B_D0optcombiMass", date); 

  // Print chi2 otpimization results 
  std::cout << "\n\nD0 candidate mass SNR with \n" << std::endl;
  std::cout << "a simple KVF           : " <<  SNR_B_D0Mass[0] << " +/- " << SNR_B_D0Mass[1] << std::endl;
  std::cout << "a constrained KVF      : " <<  SNR_B_D0consMass[0] << " +/- " << SNR_B_D0consMass[1] << std::endl;
  std::cout << "a simple combination   : " <<  SNR_B_D0combiMass[0] << " +/- " << SNR_B_D0combiMass[1] << std::endl;
  std::cout << "a biased combination   : " <<  SNR_B_D0optcombiMass[0] << " +/- " << SNR_B_D0optcombiMass[1] << std::endl;

  if (!inBatch) getchar();

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // B_D0 mass single distributions
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  TH1D* h_B_D0Mass = (TH1D*)fi->Get("ana/h_B_D0Mass");
  h_B_D0Mass->Rebin(2);
//  h_B_D0Mass->GetXaxis()->SetRangeUser(0,5);
  h_B_D0Mass->SetName("Simple Kalman Vertex Fit");
  h1_style(h_B_D0Mass, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 1, "D^{0} mass (GeV/c^{2})");
  TCanvas* cn_B_D0Mass = new TCanvas("cn_B_D0Mass","cn_B_D0Mass",800,800);
  cn_B_D0Mass->cd();
  h_B_D0Mass->Draw("hist");
  cms_style();
  cn_B_D0Mass->SaveAs("Plots"+date+"/B_D0Mass.C");
  cn_B_D0Mass->SaveAs("Plots"+date+"/B_D0Mass.eps");
  cn_B_D0Mass->SaveAs("Plots"+date+"/B_D0Mass.pdf");

  TH1D* h_B_D0consMass = (TH1D*)fi->Get("ana/h_B_D0consMass");
  h_B_D0consMass->Rebin(2);
//  h_B_D0consMass->GetXaxis()->SetRangeUser(0,5);
  h_B_D0consMass->SetName("Constrained Kalman Vertex Fit");
  h1_style(h_B_D0consMass, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 1, "D^{0} mass (GeV/c^{2})");
  TCanvas* cn_B_D0consMass = new TCanvas("cn_B_D0consMass","cn_B_D0consMass",800,800);
  cn_B_D0consMass->cd();
  h_B_D0consMass->Draw("hist");
  cms_style();
  cn_B_D0consMass->SaveAs("Plots"+date+"/B_D0consMass.C");
  cn_B_D0consMass->SaveAs("Plots"+date+"/B_D0consMass.eps");
  cn_B_D0consMass->SaveAs("Plots"+date+"/B_D0consMass.pdf");

  TH1D* h_B_D0combiMass = (TH1D*)fi->Get("ana/h_B_D0combiMass");
  h_B_D0combiMass->Rebin(2);
//  h_B_D0combiMass->GetXaxis()->SetRangeUser(0,5);
  h_B_D0combiMass->SetName("Simple PF combination");
  h1_style(h_B_D0combiMass, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 1, "D^{0} mass (GeV/c^{2})");
  TCanvas* cn_B_D0combiMass = new TCanvas("cn_B_D0combiMass","cn_B_D0combiMass",800,800);
  cn_B_D0combiMass->cd();
  h_B_D0combiMass->Draw("hist");
  cms_style();
  cn_B_D0combiMass->SaveAs("Plots"+date+"/B_D0combiMass.C");
  cn_B_D0combiMass->SaveAs("Plots"+date+"/B_D0combiMass.eps");
  cn_B_D0combiMass->SaveAs("Plots"+date+"/B_D0combiMass.pdf");

  TH1D* h_B_D0optcombiMass = (TH1D*)fi->Get("ana/h_B_D0optcombiMass");
  h_B_D0optcombiMass->Rebin(2);
//  h_B_D0optcombiMass->GetXaxis()->SetRangeUser(0,5);
  h_B_D0optcombiMass->SetName("Biased PF combination");
  h1_style(h_B_D0optcombiMass, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 1, "D^{0} mass (GeV/c^{2})");
  TCanvas* cn_B_D0optcombiMass = new TCanvas("cn_B_D0optcombiMass","cn_B_D0optcombiMass",800,800);
  cn_B_D0optcombiMass->cd();
  h_B_D0optcombiMass->Draw("hist");
  cms_style();
  cn_B_D0optcombiMass->SaveAs("Plots"+date+"/B_D0optcombiMass.C");
  cn_B_D0optcombiMass->SaveAs("Plots"+date+"/B_D0optcombiMass.eps");
  cn_B_D0optcombiMass->SaveAs("Plots"+date+"/B_D0optcombiMass.pdf");

  if (!inBatch) getchar();
  delete cn_B_D0Mass; delete cn_B_D0consMass; delete cn_B_D0combiMass; delete cn_B_D0optcombiMass; 
  delete h_B_D0Mass; delete h_B_D0consMass; delete h_B_D0combiMass; delete h_B_D0optcombiMass;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // B_D0 mass overlayed distributions
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  TH1D* h_B_D0Mass_zoom = (TH1D*)fi->Get("ana/h_B_D0Mass");
  h_B_D0Mass_zoom->SetName("Simple Kalman Vertex Fit");
  //h_B_D0Mass_zoom->Scale(1./h_B_D0Mass_zoom->Integral());

  RooRealVar x_B_D0Mass_zoom("mass","D^{0} mass",1.7,2.,"GeV/c^{2}");
  RooDataHist dh_B_D0Mass_zoom("datahist","datahist",RooArgList(x_B_D0Mass_zoom),h_B_D0Mass_zoom,1.);
  
  // Signal+Background pdf
  RooRealVar mean_B_D0Mass_zoom("mean","mean",1.86484,1.84554,1.88414);
  RooRealVar sigma_B_D0Mass_zoom("sigma","sigma",0.0193,0.01,0.03);
  RooRealVar lambda_B_D0Mass_zoom("lambda", "slope", -2., -5., -1.5);
  RooGaussian sig_B_D0Mass_zoom("sig","gaussian signal",x_B_D0Mass_zoom,mean_B_D0Mass_zoom,sigma_B_D0Mass_zoom);
  RooExponential bck_B_D0Mass_zoom("bck", "exponential background", x_B_D0Mass_zoom, lambda_B_D0Mass_zoom);
  RooRealVar nsig_B_D0Mass_zoom("nsig","number of signal events",0.1*h_B_D0Mass_zoom->Integral(),0.,h_B_D0Mass_zoom->Integral());
  RooRealVar nbck_B_D0Mass_zoom("nbck","number of background events",0.9*h_B_D0Mass_zoom->Integral(),0.,h_B_D0Mass_zoom->Integral());
  RooAddPdf model_B_D0Mass_zoom("model","model",RooArgList(sig_B_D0Mass_zoom,bck_B_D0Mass_zoom),RooArgList(nsig_B_D0Mass_zoom,nbck_B_D0Mass_zoom)) ;

  // Fit
  model_B_D0Mass_zoom.fitTo(dh_B_D0Mass_zoom);
  double mass_B_D0Mass_zoom     = mean_B_D0Mass_zoom.getVal();
  double masserr_B_D0Mass_zoom  = mean_B_D0Mass_zoom.getError();
  double width_B_D0Mass_zoom    = sigma_B_D0Mass_zoom.getVal();
  double widtherr_B_D0Mass_zoom = sigma_B_D0Mass_zoom.getError();
  double Nsig_B_D0Mass_zoom  = nsig_B_D0Mass_zoom.getVal();
  double dNsig_B_D0Mass_zoom = nsig_B_D0Mass_zoom.getError();
  double Nbck_B_D0Mass_zoom  = nbck_B_D0Mass_zoom.getVal();
  double dNbck_B_D0Mass_zoom = nbck_B_D0Mass_zoom.getError();
  double *SNR_B_D0Mass_zoom = new double[2];
  SNR_B_D0Mass_zoom[0] = Nsig_B_D0Mass_zoom*pow(Nbck_B_D0Mass_zoom,-0.5);
  SNR_B_D0Mass_zoom[1]  = dNsig_B_D0Mass_zoom*pow(Nbck_B_D0Mass_zoom,-0.5) + 0.5*SNR_B_D0Mass_zoom[0]*dNbck_B_D0Mass_zoom/Nbck_B_D0Mass_zoom;
  std::cout << "\nSignal mean = " << mass_B_D0Mass_zoom << " +/- " << masserr_B_D0Mass_zoom << std::endl;
  std::cout << "Signal width = " << width_B_D0Mass_zoom << " +/- " << widtherr_B_D0Mass_zoom << "\n" << std::endl;
  std::cout << "Nsig = " << Nsig_B_D0Mass_zoom << " +/- " << dNsig_B_D0Mass_zoom << std::endl;
  std::cout << "Nbck = " << Nbck_B_D0Mass_zoom << " +/- " << dNbck_B_D0Mass_zoom << std::endl;
  std::cout << "SNR  = " << SNR_B_D0Mass_zoom[0] << " +/- " << SNR_B_D0Mass_zoom[1] << "\n" << std::endl;

  TH1D* h_B_D0combiMass_zoom = (TH1D*)fi->Get("ana/h_B_D0combiMass");
  h_B_D0combiMass_zoom->SetName("Simple PF combination");
  h_B_D0combiMass_zoom->Scale(h_B_D0Mass_zoom->Integral()/h_B_D0combiMass_zoom->Integral());

  RooRealVar x_B_D0combiMass_zoom("mass","D^{0} mass",1.7,2.,"GeV/c^{2}");
  RooDataHist dh_B_D0combiMass_zoom("datahist","datahist",RooArgList(x_B_D0combiMass_zoom),h_B_D0combiMass_zoom,1.);
  
  // Signal+Background pdf
  RooRealVar mean_B_D0combiMass_zoom("mean","mean",1.86484,1.84554,1.88414);
  RooRealVar sigma_B_D0combiMass_zoom("sigma","sigma",0.0193,0.01,0.03);
  RooRealVar lambda_B_D0combiMass_zoom("lambda", "slope", -1.7, -5., -0.7);
  RooGaussian sig_B_D0combiMass_zoom("sig","gaussian signal",x_B_D0combiMass_zoom,mean_B_D0combiMass_zoom,sigma_B_D0combiMass_zoom);
  RooExponential bck_B_D0combiMass_zoom("bck", "exponential background", x_B_D0combiMass_zoom, lambda_B_D0combiMass_zoom);
  RooRealVar nsig_B_D0combiMass_zoom("nsig","number of signal events",0.1*h_B_D0combiMass_zoom->Integral(),0.,h_B_D0combiMass_zoom->Integral());
  RooRealVar nbck_B_D0combiMass_zoom("nbck","number of background events",0.9*h_B_D0combiMass_zoom->Integral(),0.,h_B_D0combiMass_zoom->Integral());
  RooAddPdf model_B_D0combiMass_zoom("model","model",RooArgList(sig_B_D0combiMass_zoom,bck_B_D0combiMass_zoom),RooArgList(nsig_B_D0combiMass_zoom,nbck_B_D0combiMass_zoom)) ;

  // Fit
  model_B_D0combiMass_zoom.fitTo(dh_B_D0combiMass_zoom);
  double mass_B_D0combiMass_zoom     = mean_B_D0combiMass_zoom.getVal();
  double masserr_B_D0combiMass_zoom  = mean_B_D0combiMass_zoom.getError();
  double width_B_D0combiMass_zoom    = sigma_B_D0combiMass_zoom.getVal();
  double widtherr_B_D0combiMass_zoom = sigma_B_D0combiMass_zoom.getError();
  double Nsig_B_D0combiMass_zoom  = nsig_B_D0combiMass_zoom.getVal();
  double dNsig_B_D0combiMass_zoom = nsig_B_D0combiMass_zoom.getError();
  double Nbck_B_D0combiMass_zoom  = nbck_B_D0combiMass_zoom.getVal();
  double dNbck_B_D0combiMass_zoom = nbck_B_D0combiMass_zoom.getError();
  double *SNR_B_D0combiMass_zoom = new double[2];
  SNR_B_D0combiMass_zoom[0] = Nsig_B_D0combiMass_zoom*pow(Nbck_B_D0combiMass_zoom,-0.5);
  SNR_B_D0combiMass_zoom[1]  = dNsig_B_D0combiMass_zoom*pow(Nbck_B_D0combiMass_zoom,-0.5) + 0.5*SNR_B_D0combiMass_zoom[0]*dNbck_B_D0combiMass_zoom/Nbck_B_D0combiMass_zoom;
  std::cout << "\nSignal mean = " << mass_B_D0combiMass_zoom << " +/- " << masserr_B_D0combiMass_zoom << std::endl;
  std::cout << "Signal width = " << width_B_D0combiMass_zoom << " +/- " << widtherr_B_D0combiMass_zoom << "\n" << std::endl;
  std::cout << "Nsig = " << Nsig_B_D0combiMass_zoom << " +/- " << dNsig_B_D0combiMass_zoom << std::endl;
  std::cout << "Nbck = " << Nbck_B_D0combiMass_zoom << " +/- " << dNbck_B_D0combiMass_zoom << std::endl;
  std::cout << "SNR  = " << SNR_B_D0combiMass_zoom[0] << " +/- " << SNR_B_D0combiMass_zoom[1] << "\n" << std::endl;

  RooRealVar x_B_D0unbias_zoom("mass","D^{0} mass",1.7,2.,"GeV/c^{2}");
  RooPlot* frame_B_D0unbias_zoom = x_B_D0unbias_zoom.frame();
  dh_B_D0Mass_zoom.plotOn(frame_B_D0unbias_zoom,MarkerColor(50),LineColor(50),Name("model_B_D0Mass_zoom"));
  model_B_D0Mass_zoom.plotOn(frame_B_D0unbias_zoom,LineColor(50));
  dh_B_D0combiMass_zoom.plotOn(frame_B_D0unbias_zoom,MarkerColor(38),LineColor(38),Name("model_B_D0combiMass_zoom"));
  model_B_D0combiMass_zoom.plotOn(frame_B_D0unbias_zoom,LineColor(38));
  TCanvas* cn_B_D0unbias_zoom = new TCanvas("cn_B_D0unbias_zoom","cn_B_D0unbias_zoom",800,800);
  cn_B_D0unbias_zoom->cd();
  frame_B_D0unbias_zoom->Draw();
  TLegend *leg_B_D0unbias_zoom = new TLegend(0.7,0.72,0.94,0.92);
  leg_style(leg_B_D0unbias_zoom,12);
  leg_B_D0unbias_zoom->AddEntry("model_B_D0combiMass_zoom","#splitline{Simple PF}{combination}","LP");
  leg_B_D0unbias_zoom->AddEntry("model_B_D0Mass_zoom","#splitline{Simple Kalman}{Vertex Fit}","LP");
  leg_B_D0unbias_zoom->Draw();
  cms_style();
  cn_B_D0unbias_zoom->SaveAs("Plots"+date+"/B_D0unbias_zoom.C");
  cn_B_D0unbias_zoom->SaveAs("Plots"+date+"/B_D0unbias_zoom.eps");
  cn_B_D0unbias_zoom->SaveAs("Plots"+date+"/B_D0unbias_zoom.pdf");

  if (!inBatch) getchar();
  delete cn_B_D0unbias_zoom;

  TH1D* h_B_D0consMass_zoom = (TH1D*)fi->Get("ana/h_B_D0consMass");
  h_B_D0consMass_zoom->SetName("Constrained Kalman Vertex Fit");

  RooRealVar x_B_D0consMass_zoom("mass","D^{0} mass",1.7,2.,"GeV/c^{2}");
  RooDataHist dh_B_D0consMass_zoom("datahist","datahist",RooArgList(x_B_D0consMass_zoom),h_B_D0consMass_zoom,1.);
  
  // Signal+Background pdf
  RooRealVar mean_B_D0consMass_zoom("mean","mean",1.86484,1.84554,1.88414);
  RooRealVar sigma_B_D0consMass_zoom("sigma","sigma",0.0193,0.01,0.03);
  RooRealVar lambda_B_D0consMass_zoom("lambda", "slope", -2., -5., -1.5);
  RooGaussian sig_B_D0consMass_zoom("sig","gaussian signal",x_B_D0consMass_zoom,mean_B_D0consMass_zoom,sigma_B_D0consMass_zoom);
  RooExponential bck_B_D0consMass_zoom("bck", "exponential background", x_B_D0consMass_zoom, lambda_B_D0consMass_zoom);
  RooRealVar nsig_B_D0consMass_zoom("nsig","number of signal events",0.1*h_B_D0consMass_zoom->Integral(),0.,h_B_D0consMass_zoom->Integral());
  RooRealVar nbck_B_D0consMass_zoom("nbck","number of background events",0.9*h_B_D0consMass_zoom->Integral(),0.,h_B_D0consMass_zoom->Integral());
  RooAddPdf model_B_D0consMass_zoom("model","model",RooArgList(sig_B_D0consMass_zoom,bck_B_D0consMass_zoom),RooArgList(nsig_B_D0consMass_zoom,nbck_B_D0consMass_zoom)) ;

  // Fit
  model_B_D0consMass_zoom.fitTo(dh_B_D0consMass_zoom);
  double mass_B_D0consMass_zoom     = mean_B_D0consMass_zoom.getVal();
  double masserr_B_D0consMass_zoom  = mean_B_D0consMass_zoom.getError();
  double width_B_D0consMass_zoom    = sigma_B_D0consMass_zoom.getVal();
  double widtherr_B_D0consMass_zoom = sigma_B_D0consMass_zoom.getError();
  double Nsig_B_D0consMass_zoom  = nsig_B_D0consMass_zoom.getVal();
  double dNsig_B_D0consMass_zoom = nsig_B_D0consMass_zoom.getError();
  double Nbck_B_D0consMass_zoom  = nbck_B_D0consMass_zoom.getVal();
  double dNbck_B_D0consMass_zoom = nbck_B_D0consMass_zoom.getError();
  double *SNR_B_D0consMass_zoom = new double[2];
  SNR_B_D0consMass_zoom[0] = Nsig_B_D0consMass_zoom*pow(Nbck_B_D0consMass_zoom,-0.5);
  SNR_B_D0consMass_zoom[1]  = dNsig_B_D0consMass_zoom*pow(Nbck_B_D0consMass_zoom,-0.5) + 0.5*SNR_B_D0consMass_zoom[0]*dNbck_B_D0consMass_zoom/Nbck_B_D0consMass_zoom;
  std::cout << "\nSignal mean = " << mass_B_D0consMass_zoom << " +/- " << masserr_B_D0consMass_zoom << std::endl;
  std::cout << "Signal width = " << width_B_D0consMass_zoom << " +/- " << widtherr_B_D0consMass_zoom << "\n" << std::endl;
  std::cout << "Nsig = " << Nsig_B_D0consMass_zoom << " +/- " << dNsig_B_D0consMass_zoom << std::endl;
  std::cout << "Nbck = " << Nbck_B_D0consMass_zoom << " +/- " << dNbck_B_D0consMass_zoom << std::endl;
  std::cout << "SNR  = " << SNR_B_D0consMass_zoom[0] << " +/- " << SNR_B_D0consMass_zoom[1] << "\n" << std::endl;

  TH1D* h_B_D0optcombiMass_zoom = (TH1D*)fi->Get("ana/h_B_D0optcombiMass");
  h_B_D0optcombiMass_zoom->SetName("Biased PF combination");

  RooRealVar x_B_D0optcombiMass_zoom("mass","D^{0} mass",1.7,2.,"GeV/c^{2}");
  RooDataHist dh_B_D0optcombiMass_zoom("datahist","datahist",RooArgList(x_B_D0optcombiMass_zoom),h_B_D0optcombiMass_zoom,1.);
  
  // Signal+Background pdf
  RooRealVar mean_B_D0optcombiMass_zoom("mean","mean",1.86484,1.84554,1.88414);
  RooRealVar sigma_B_D0optcombiMass_zoom("sigma","sigma",0.0193,0.01,0.03);
  RooRealVar lambda_B_D0optcombiMass_zoom("lambda", "slope", -1.7, -5., -0.7);
  RooGaussian sig_B_D0optcombiMass_zoom("sig","gaussian signal",x_B_D0optcombiMass_zoom,mean_B_D0optcombiMass_zoom,sigma_B_D0optcombiMass_zoom);
  RooExponential bck_B_D0optcombiMass_zoom("bck", "exponential background", x_B_D0optcombiMass_zoom, lambda_B_D0optcombiMass_zoom);
  RooRealVar nsig_B_D0optcombiMass_zoom("nsig","number of signal events",0.1*h_B_D0optcombiMass_zoom->Integral(),0.,h_B_D0optcombiMass_zoom->Integral());
  RooRealVar nbck_B_D0optcombiMass_zoom("nbck","number of background events",0.9*h_B_D0optcombiMass_zoom->Integral(),0.,h_B_D0optcombiMass_zoom->Integral());
  RooAddPdf model_B_D0optcombiMass_zoom("model","model",RooArgList(sig_B_D0optcombiMass_zoom,bck_B_D0optcombiMass_zoom),RooArgList(nsig_B_D0optcombiMass_zoom,nbck_B_D0optcombiMass_zoom)) ;

  // Fit
  model_B_D0optcombiMass_zoom.fitTo(dh_B_D0optcombiMass_zoom);
  double mass_B_D0optcombiMass_zoom     = mean_B_D0optcombiMass_zoom.getVal();
  double masserr_B_D0optcombiMass_zoom  = mean_B_D0optcombiMass_zoom.getError();
  double width_B_D0optcombiMass_zoom    = sigma_B_D0optcombiMass_zoom.getVal();
  double widtherr_B_D0optcombiMass_zoom = sigma_B_D0optcombiMass_zoom.getError();
  double Nsig_B_D0optcombiMass_zoom  = nsig_B_D0optcombiMass_zoom.getVal();
  double dNsig_B_D0optcombiMass_zoom = nsig_B_D0optcombiMass_zoom.getError();
  double Nbck_B_D0optcombiMass_zoom  = nbck_B_D0optcombiMass_zoom.getVal();
  double dNbck_B_D0optcombiMass_zoom = nbck_B_D0optcombiMass_zoom.getError();
  double *SNR_B_D0optcombiMass_zoom = new double[2];
  SNR_B_D0optcombiMass_zoom[0] = Nsig_B_D0optcombiMass_zoom*pow(Nbck_B_D0optcombiMass_zoom,-0.5);
  SNR_B_D0optcombiMass_zoom[1]  = dNsig_B_D0optcombiMass_zoom*pow(Nbck_B_D0optcombiMass_zoom,-0.5) + 0.5*SNR_B_D0optcombiMass_zoom[0]*dNbck_B_D0optcombiMass_zoom/Nbck_B_D0optcombiMass_zoom;
  std::cout << "\nSignal mean = " << mass_B_D0optcombiMass_zoom << " +/- " << masserr_B_D0optcombiMass_zoom << std::endl;
  std::cout << "Signal width = " << width_B_D0optcombiMass_zoom << " +/- " << widtherr_B_D0optcombiMass_zoom << "\n" << std::endl;
  std::cout << "Nsig = " << Nsig_B_D0optcombiMass_zoom << " +/- " << dNsig_B_D0optcombiMass_zoom << std::endl;
  std::cout << "Nbck = " << Nbck_B_D0optcombiMass_zoom << " +/- " << dNbck_B_D0optcombiMass_zoom << std::endl;
  std::cout << "SNR  = " << SNR_B_D0optcombiMass_zoom[0] << " +/- " << SNR_B_D0optcombiMass_zoom[1] << "\n" << std::endl;

  RooRealVar x_B_D0bias_zoom("mass","D^{0} mass",1.7,2.,"GeV/c^{2}");
  RooPlot* frame_B_D0bias_zoom = x_B_D0bias_zoom.frame();
  dh_B_D0consMass_zoom.plotOn(frame_B_D0bias_zoom,MarkerColor(47),LineColor(47),Name("model_B_D0consMass_zoom"));
  model_B_D0consMass_zoom.plotOn(frame_B_D0bias_zoom,LineColor(47));
  dh_B_D0optcombiMass_zoom.plotOn(frame_B_D0bias_zoom,MarkerColor(30),LineColor(30),Name("model_B_D0optcombiMass_zoom"));
  model_B_D0optcombiMass_zoom.plotOn(frame_B_D0bias_zoom,LineColor(30));
  TCanvas* cn_B_D0bias_zoom = new TCanvas("cn_B_D0bias_zoom","cn_B_D0bias_zoom",800,800);
  cn_B_D0bias_zoom->cd();
  frame_B_D0bias_zoom->Draw();
  TLegend *leg_B_D0bias_zoom = new TLegend(0.7,0.72,0.94,0.92);
  leg_style(leg_B_D0bias_zoom,12);
  leg_B_D0bias_zoom->AddEntry("model_B_D0optcombiMass_zoom","#splitline{Biased PF}{combination}","LP");
  leg_B_D0bias_zoom->AddEntry("model_B_D0consMass_zoom","#splitline{Constrained Kalman}{Vertex Fit}","LP");
  leg_B_D0bias_zoom->Draw();
  cms_style();
  cn_B_D0bias_zoom->SaveAs("Plots"+date+"/B_D0bias_zoom.C");
  cn_B_D0bias_zoom->SaveAs("Plots"+date+"/B_D0bias_zoom.eps");
  cn_B_D0bias_zoom->SaveAs("Plots"+date+"/B_D0bias_zoom.pdf");

  if (!inBatch) getchar();
  delete cn_B_D0bias_zoom;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // have a look at other distributions 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  TH1D* h_nJets = (TH1D*)fi->Get("ana/h_nJets");
  h1_style(h_nJets, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "Number of jets");
  h_nJets->GetYaxis()->SetNoExponent();
  TCanvas* cn_nJets = new TCanvas("cn_nJets","cn_nJets",800,800);
  cn_nJets->cd();
  h_nJets->Draw("hist");
  cms_style();
  cn_nJets->SaveAs("Plots"+date+"/nJets.C");
  cn_nJets->SaveAs("Plots"+date+"/nJets.eps");
  cn_nJets->SaveAs("Plots"+date+"/nJets.pdf");

  TH1D* h_CSV = (TH1D*)fi->Get("ana/h_CSV");
  h1_style(h_CSV, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "CSV discriminant");
  h_CSV->GetYaxis()->SetNoExponent();
  TCanvas* cn_CSV = new TCanvas("cn_CSV","cn_CSV",800,800);
  cn_CSV->cd();
  h_CSV->Draw("hist");
  cms_style();
  cn_CSV->SaveAs("Plots"+date+"/CSV.C");
  cn_CSV->SaveAs("Plots"+date+"/CSV.eps");
  cn_CSV->SaveAs("Plots"+date+"/CSV.pdf");

  TH1D* h_B_cuts = (TH1D*)fi->Get("ana/h_B_cuts");
  h1_style(h_B_cuts, 38, 38, 3003, -1111., -1111., 510, 510, 38, 1.2, 0, "");
  h_B_cuts->GetXaxis()->SetRangeUser(0,11);
  h_B_cuts->GetYaxis()->SetNoExponent();
  h_B_cuts->GetYaxis()->SetTitle();
  TCanvas* cn_B_cuts = new TCanvas("cn_B_cuts","cn_B_cuts",800,800);
  cn_B_cuts->cd();
  cn_B_cuts->SetLogy(1);
  h_B_cuts->Draw("hist");
  cms_style();
  cn_B_cuts->SaveAs("Plots"+date+"/B_cuts.C");
  cn_B_cuts->SaveAs("Plots"+date+"/B_cuts.eps");
  cn_B_cuts->SaveAs("Plots"+date+"/B_cuts.pdf");

  if (!inBatch) getchar();
  delete cn_nJets; delete cn_CSV; delete cn_B_cuts;
  delete h_nJets; delete h_CSV; delete h_B_cuts;

  //close
  fi->Close();
}

//---------------------------------------------------------------
int doTheComparison_Data(bool inBatch = true, TString date = "")
//---------------------------------------------------------------
{

  doTheComparison_step(inBatch, "El_" + date, "../test/crab_results/17Dec14/kalmanAnalyzed_El_merged.root");
  doTheComparison_step(inBatch, "Mu_" + date, "../test/crab_results/17Dec14/kalmanAnalyzed_Mu_merged.root");
  return 0;
}

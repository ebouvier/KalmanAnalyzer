#include <iomanip>
#include <string>
#include <iostream>
#include <sstream>
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TFile.h"
#include "THStack.h"
#include "TH1D.h"
#include "TCanvas.h"

#pragma once

#define TITLE_FONTSIZE 26
#define LABEL_FONTSIZE 18

#define LEFT_MARGIN 0.17
#define RIGHT_MARGIN 0.03
#define TOP_MARGIN 0.05
#define BOTTOM_MARGIN 0.13

TStyle* createMyStyle() {
  TStyle *myStyle = new TStyle("myStyle", "myStyle");

  // For the canvas:
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetCanvasColor(kWhite);
  myStyle->SetCanvasDefH(800); //Height of canvas
  myStyle->SetCanvasDefW(800); //Width of canvas
  myStyle->SetCanvasDefX(0);   //POsition on screen
  myStyle->SetCanvasDefY(0);

  // For the Pad:
  myStyle->SetPadBorderMode(0);
  myStyle->SetPadColor(kWhite);
  myStyle->SetPadGridX(false);
  myStyle->SetPadGridY(false);
  myStyle->SetGridColor(0);
  myStyle->SetGridStyle(3);
  myStyle->SetGridWidth(1);

  // For the frame:
  myStyle->SetFrameBorderMode(0);
  myStyle->SetFrameBorderSize(1);
  myStyle->SetFrameFillColor(0);
  myStyle->SetFrameFillStyle(0);
  myStyle->SetFrameLineColor(1);
  myStyle->SetFrameLineStyle(1);
  myStyle->SetFrameLineWidth(1);

  // For the histo:
  myStyle->SetHistLineStyle(1);
  myStyle->SetHistLineWidth(2);
  myStyle->SetEndErrorSize(2);

  //For the fit/function:
  myStyle->SetFitFormat("5.4g");
  myStyle->SetFuncColor(2);
  myStyle->SetFuncStyle(1);
  myStyle->SetFuncWidth(1);

  // For the statistics box:
  myStyle->SetOptFile(0);
  myStyle->SetStatColor(kWhite);
  //myStyle->SetStatFont(43);
  //myStyle->SetStatFontSize(0.025);
  myStyle->SetStatTextColor(1);
  myStyle->SetStatFormat("6.4g");
  myStyle->SetStatBorderSize(1);
  myStyle->SetStatH(0.12);
  myStyle->SetStatW(0.3);
  myStyle->SetStatY(0.92);
  myStyle->SetStatX(0.94);

  //For the date:
  myStyle->SetOptDate(0);

  // Margins:
  myStyle->SetPadTopMargin(TOP_MARGIN);
  myStyle->SetPadBottomMargin(BOTTOM_MARGIN);
  myStyle->SetPadLeftMargin(LEFT_MARGIN);
  myStyle->SetPadRightMargin(RIGHT_MARGIN);

  // For the Global title:
  myStyle->SetOptTitle(0);
  myStyle->SetTitleFont(63);
  myStyle->SetTitleColor(1);
  myStyle->SetTitleTextColor(1);
  myStyle->SetTitleFillColor(10);
  myStyle->SetTitleBorderSize(0);
  myStyle->SetTitleAlign(33); 
  myStyle->SetTitleX(1);
  myStyle->SetTitleFontSize(TITLE_FONTSIZE);

  // For the axis titles:

  myStyle->SetTitleColor(1, "XYZ");
  myStyle->SetTitleFont(43, "XYZ");
  myStyle->SetTitleSize(TITLE_FONTSIZE, "XYZ");
  myStyle->SetTitleYOffset(2.5); 
  myStyle->SetTitleXOffset(1.5);

  myStyle->SetLabelColor(1, "XYZ");
  myStyle->SetLabelFont(43, "XYZ");
  myStyle->SetLabelOffset(0.01, "YZ");
  myStyle->SetLabelOffset(0.015, "X");
  myStyle->SetLabelSize(LABEL_FONTSIZE, "XYZ");

  myStyle->SetAxisColor(1, "XYZ");
  myStyle->SetStripDecimals(kTRUE);
  myStyle->SetTickLength(0.03, "XYZ");
  myStyle->SetNdivisions(510, "XYZ");
  myStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  myStyle->SetPadTickY(1);

  myStyle->SetOptLogx(0);
  myStyle->SetOptLogy(0);
  myStyle->SetOptLogz(0);

  myStyle->SetHatchesSpacing(1.3);
  myStyle->SetHatchesLineWidth(1);

  myStyle->cd();

  return myStyle;
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

void h_style(TH1 *h,
    int line_color=1,
    int fill_color=50,
    int fill_style=1001,
    float y_min=-1111.,
    float y_max=-1111.,
    int ndivx=510,
    int ndivy=510,
    int marker_style=20,
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

  h->SetMarkerStyle(marker_style);
  h->SetMarkerColor(marker_color);
  h->SetMarkerSize(marker_size);
  h->SetStats(optstat);

  h->GetXaxis()->SetTitle(xtitle);
  float binSize = h->GetXaxis()->GetBinWidth(1);
  std::stringstream ss;
  ss << "Events / " << std::fixed << std::setprecision(2) << binSize;
  h->GetYaxis()->SetTitle(ss.str().c_str());
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
void plotHisto_1(bool inBatch, TFile* fi_pT, TFile* fi_muTag, TString h_name, TString xtitle, double xmin, double xmax)
//---------------------------------------------------------------
{
  TStyle* my_style = createMyStyle();
  my_style->cd();
  if (inBatch) gROOT->SetBatch(true);
  
  TLatex* channel_tex = new TLatex(0.8, 0.9, "X-check");
  channel_tex->SetNDC(true);
  channel_tex->SetTextFont(43);
  channel_tex->SetTextSize(TITLE_FONTSIZE - 6);  
  
  TH1D* h_pT = (TH1D*)fi_pT->Get("pTbasedSelection/"+h_name+"-b-jets");
  TH1D* h_muTag = (TH1D*)fi_muTag->Get("muTagBasedSelection/"+h_name+"-b-jets");
  h_pT->Scale(h_muTag->Integral()/h_pT->Integral());
  h_style(h_pT, kBlue, kBlue, 4000, -1111., -1111., 510, 510, 24, kBlue, 1.2, 0, xtitle);
  h_style(h_muTag, kRed, kRed, 4000, -1111., -1111., 510, 510, 20, kRed, 1.2, 0, xtitle);
  
  TCanvas* cn = new TCanvas("cn","cn",800,800);
  cn->cd();
  if (h_muTag->GetMaximum() < h_pT->GetMaximum()) {
    h_pT->Draw("ep");
    if (xmin < xmax) h_pT->GetXaxis()->SetRangeUser(xmin, xmax);
    h_muTag->Draw("epsame");
  }
  else {
    h_muTag->Draw("ep");
    if (xmin < xmax) h_muTag->GetXaxis()->SetRangeUser(xmin, xmax);
    h_pT->Draw("epsame");
    h_muTag->Draw("epsame");
  }
  TLegend *leg = new TLegend(0.64,0.65,0.93,0.84);
  leg_style(leg, 12);
  leg->AddEntry(h_pT,"2 jets of highest p_{T}","LP");
  leg->AddEntry(h_muTag,"all jets with a #mu^{#pm}","LP");
  leg->SetHeader("Data - Run 2012 A,B,C,D");
  leg->Draw();
  
  channel_tex->Draw("same");  
  cms_style();
  cn->SaveAs("11Feb15/muTagVSpT/"+h_name+".C");
  cn->SaveAs("11Feb15/muTagVSpT/"+h_name+".eps");
  cn->SaveAs("11Feb15/muTagVSpT/"+h_name+".pdf");
  
  if (!inBatch) getchar();
  
  delete cn;
  delete h_pT;
  delete h_muTag;
  delete channel_tex;
  delete my_style;
}

//---------------------------------------------------------------
void plotHisto_2(bool inBatch, TFile* fi_pT, TFile* fi_muTag, TString h_name, TString xtitle, double xmin, double xmax)
//---------------------------------------------------------------
{
  TStyle* my_style = createMyStyle();
  my_style->cd();
  if (inBatch) gROOT->SetBatch(true);
  
  TLatex* channel_tex = new TLatex(0.8, 0.9, "X-check");
  channel_tex->SetNDC(true);
  channel_tex->SetTextFont(43);
  channel_tex->SetTextSize(TITLE_FONTSIZE - 6);  
  
  TH1D* h_pT = (TH1D*)fi_pT->Get("pTbasedSelection/"+h_name+"-b-jet1");
  TH1D* h_pT_2 = (TH1D*)fi_pT->Get("pTbasedSelection/"+h_name+"-b-jet2");
  TH1D* h_muTag = (TH1D*)fi_muTag->Get("muTagBasedSelection/"+h_name+"-b-jets");
  h_pT->Add(h_pT_2); 
  h_pT->Scale(h_muTag->Integral()/h_pT->Integral());
  h_style(h_pT, kBlue, kBlue, 4000, -1111., -1111., 510, 510, 24, kBlue, 1.2, 0, xtitle);
  h_style(h_muTag, kRed, kRed, 4000, -1111., -1111., 510, 510, 20, kRed, 1.2, 0, xtitle);

  TCanvas* cn = new TCanvas("cn","cn",800,800);
  cn->cd();
  if (h_muTag->GetMaximum() < h_pT->GetMaximum()) {
    h_pT->Draw("ep");
    if (xmin < xmax) h_pT->GetXaxis()->SetRangeUser(xmin, xmax);
    h_muTag->Draw("epsame");
  }
  else {
    h_muTag->Draw("ep");
    if (xmin < xmax) h_muTag->GetXaxis()->SetRangeUser(xmin, xmax);
    h_pT->Draw("epsame");
    h_muTag->Draw("epsame");
  }
  TLegend *leg = new TLegend(0.64,0.65,0.93,0.84);
  leg_style(leg, 12);
  leg->AddEntry(h_pT,"2 jets of highest p_{T}","LP");
  leg->AddEntry(h_muTag,"all jets with a #mu^{#pm}","LP");
  leg->SetHeader("Data - Run 2012 A,B,C,D");
  leg->Draw();
  
  channel_tex->Draw("same");  
  cms_style();
  cn->SaveAs("11Feb15/muTagVSpT/"+h_name+".C");
  cn->SaveAs("11Feb15/muTagVSpT/"+h_name+".eps");
  cn->SaveAs("11Feb15/muTagVSpT/"+h_name+".pdf");
  
  if (!inBatch) getchar();

  delete cn;
  delete h_pT_2;
  delete h_pT;
  delete h_muTag;
  delete channel_tex;
  delete my_style;
}

//---------------------------------------------------------------
void xcheckForRoberto_d0(bool inBatch)
//---------------------------------------------------------------
{
  TFile* fi_pT = TFile::Open("../test/crab_results/09Feb15/D0ForRivet_pT_Data_merged.root");
  TFile* fi_muTag = TFile::Open("../test/crab_results/11Feb15/MuTagForRivet_Data_merged.root");
  
  plotHisto_1(inBatch, fi_pT, fi_muTag, "CSV", "CSV discriminant", 0., 0.);
  plotHisto_1(inBatch, fi_pT, fi_muTag, "TransverseMomentum", "p(jets) (GeV/c)", 0., 300.);
  plotHisto_2(inBatch, fi_pT, fi_muTag, "Nch", "Tracks multiplicity", 0., 20.);
  plotHisto_2(inBatch, fi_pT, fi_muTag, "Sump", "Scalar sum of tracks momenta (GeV/c)", 0., 400.);
  plotHisto_2(inBatch, fi_pT, fi_muTag, "VectorialSump", "Vectorial sum of tracks momenta (GeV/c)", 0., 400.);
  plotHisto_2(inBatch, fi_pT, fi_muTag, "Highestp", "Highest track momentum (GeV/c)", 0., 200.);
  plotHisto_2(inBatch, fi_pT, fi_muTag, "Sum3p", "Scalar sum of the 3 highest tracks momenta (GeV/c)", 0., 0.);
  plotHisto_2(inBatch, fi_pT, fi_muTag, "R1", "R_{1}", 0., 0.);
  plotHisto_2(inBatch, fi_pT, fi_muTag, "R3", "R_{3}", 0., 0.);
  plotHisto_2(inBatch, fi_pT, fi_muTag, "D0Mass", "m(D^{0}#rightarrow#kappa^{+}#pi^{-}) (GeV/c^{2})", 0., 0.);
  plotHisto_2(inBatch, fi_pT, fi_muTag, "D0p", "p(D^{0}#rightarrow#kappa^{+}#pi^{-}) (GeV/c)", 0., 0.);
  plotHisto_2(inBatch, fi_pT, fi_muTag, "BMomentum", "p(#kappa^{+}#pi^{-}+#mu^{-}) (GeV/c)", 0., 350.);
  
  fi_pT->Close();
  fi_muTag->Close();
  delete fi_pT;
  delete fi_muTag;
}
  
//---------------------------------------------------------------
int xcheckForRoberto(bool inBatch = true)
//---------------------------------------------------------------
{  
  gROOT->ProcessLine(".! mkdir 11Feb15/muTagVSpT");
  xcheckForRoberto_d0(inBatch);
  return 0;
}

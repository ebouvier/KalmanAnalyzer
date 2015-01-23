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
void plotHisto(bool inBatch, TString date, bool isCSVbased, TFile* fi_data, TFile* fi_sl, TFile* fi_dl, TString channel, TString h_name, TString xtitle, double xmin, double xmax)
//---------------------------------------------------------------
{
  TStyle* my_style = createMyStyle();
  my_style->cd();
  if (inBatch) gROOT->SetBatch(true);
  
  TLatex* channel_tex = new TLatex(0.87, 0.9, channel);
  channel_tex->SetNDC(true);
  channel_tex->SetTextFont(43);
  channel_tex->SetTextSize(TITLE_FONTSIZE - 6);  

  TString dir_name = "basedSelection/";
  TString rep_name = date;
  if (isCSVbased) {
    dir_name = "CSV" + dir_name;
    rep_name = rep_name + "/csv/";
  }
  else { 
    dir_name = "pT" + dir_name;
    rep_name = rep_name + "/pT/";
  }
  
  TH1D* h_data = (TH1D*)fi_data->Get(dir_name+h_name);
  TH1D* h_sl = (TH1D*)fi_sl->Get(dir_name+h_name);
  TH1D* h_dl = (TH1D*)fi_dl->Get(dir_name+h_name);
  h_sl->Scale(2.*25.8031*19769./12119013.);
  h_dl->Scale(2.*107.6722*19769./25424818.);
  TH1D* h_mc = (TH1D*)h_sl->Clone();
  h_mc->Add(h_dl);
  
  // because we miss some data...
  h_sl->Scale(h_data->Integral()/h_mc->Integral());
  h_dl->Scale(h_data->Integral()/h_mc->Integral());
  h_mc->Scale(h_data->Integral()/h_mc->Integral());
 
  h_style(h_data, 1, 1, 4000, -1111., -1111., 510, 510, 21, 1, 1.2, 0, xtitle);
  h_style(h_sl, 2, 2, 1001, -1111., -1111., 510, 510, 21, 2, 1.2, 0, xtitle);
  h_style(h_dl, 3, 3, 1001, -1111., -1111., 510, 510, 21, 3, 1.2, 0, xtitle);
  h_style(h_mc, 0, 12, 3004, -1111., -1111., 510, 510, 1, 12, 1, 0, xtitle);

  THStack* hs_mc = new THStack("hs", "hs");
  hs_mc->Add(h_dl);
  hs_mc->Add(h_sl);

  TCanvas* cn = new TCanvas("cn","cn",800,800);
  cn->cd();
  if (h_data->GetMaximum() < h_mc->GetMaximum()) {
    h_mc->Draw("e2");
    if (xmin < xmax) h_mc->GetXaxis()->SetRangeUser(xmin, xmax);
  }
  else {
    h_data->Draw("ep");
    if (xmin < xmax) h_data->GetXaxis()->SetRangeUser(xmin, xmax);
    h_mc->Draw("e2same");
  }
  hs_mc->Draw("histsame");
  h_data->Draw("epsame");
  TLegend *leg = new TLegend(0.64,0.62,0.93,0.84);
  leg_style(leg, 12);
  leg->AddEntry(h_data,"Data - Run 2012 A,B,C,D","LP");
  leg->AddEntry(h_sl,"MG+PY6 Z2* Semilept. t#bar{t}","F");
  leg->AddEntry(h_dl,"MG+PY6 Z2* Dilept. t#bar{t}","F");
  leg->AddEntry(h_mc,"Stat. uncertainty","F");
  leg->Draw();
  
  channel_tex->Draw("same");  
  cms_style();
  cn->SaveAs(rep_name+h_name+".C");
  cn->SaveAs(rep_name+h_name+".eps");
  cn->SaveAs(rep_name+h_name+".pdf");
  
  if (!inBatch) getchar();

  delete cn;
  delete hs_mc;
  delete h_data;
  delete h_sl;
  delete h_dl;
  delete channel_tex;
  delete my_style;
}

//---------------------------------------------------------------
void doRivetJob_cut(bool inBatch, TString date, bool isCSVbased)
//---------------------------------------------------------------
{
  TString fi_data_name = "../test/crab_results/"+date;
  TString fi_sl_name = "../test/crab_results/"+date;
  TString fi_dl_name = "../test/crab_results/"+date;
  std::cout << "\t\t /!\\ Don't forget in " << fi_data_name << " /!\\ " << std::endl;
  if (isCSVbased)
    std::cout << "hadd D0ForRivet_csv_Data_merged.root D0ForRivet_csv_El_merged.root D0ForRivet_csv_Mu_merged.root" << std::endl;
  else 
    std::cout << "hadd D0ForRivet_pT_Data_merged.root D0ForRivet_pT_El_merged.root D0ForRivet_pT_Mu_merged.root" << std::endl;
  getchar();
  if (isCSVbased) {
    fi_data_name = fi_data_name + "/D0ForRivet_csv_Data_merged.root";
    fi_sl_name = fi_sl_name + "/D0ForRivet_csv_TTJets_SemiLeptMGDecays.root";
    fi_dl_name = fi_dl_name + "/D0ForRivet_csv_TTJets_FullLeptMGDecays.root";
  }
  else {
    fi_data_name = fi_data_name + "/D0ForRivet_pT_Data_merged.root";
    fi_sl_name = fi_sl_name + "/D0ForRivet_pT_TTJets_SemiLeptMGDecays.root";
    fi_dl_name = fi_dl_name + "/D0ForRivet_pT_TTJets_FullLeptMGDecays.root";
  }
  
  TFile* fi_data = TFile::Open(fi_data_name);
  TFile* fi_sl = TFile::Open(fi_sl_name);
  TFile* fi_dl = TFile::Open(fi_dl_name);

  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "All", "NPrimaryVtx", "Primary vertices multiplicity", 0., 0.);
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "All", "NJets", "Jets multiplicity", 0., 0.);
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "All", "TransverseMomentumJets", "p(jets) (GeV/c)", 0., 0.);
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "All", "NCSVJets", "b-tagged jets multiplicity", 0., 0.);

  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "Both jets", "CSV-b-jets", "CSV discriminant", 0., 0.);
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "Both jets", "TransverseMomentum-b-jets", "p(jets) (GeV/c)", 0., 0.);

  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "1^{st} jet", "Nch-b-jet1", "Tracks multiplicity", 0., 20.);
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "2^{nd} jet", "Nch-b-jet2", "Tracks multiplicity", 0., 20.);

  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "1^{st} jet", "Sump-b-jet1", "Scalar sum of tracks momenta (GeV/c)", 0., 0.);
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "2^{nd} jet", "Sump-b-jet2", "Scalar sum of tracks momenta (GeV/c)", 0., 0.);
  
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "1^{st} jet", "VectorialSump-b-jet1", "Vectorial sum of tracks momenta (GeV/c)", 0., 0.);
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "2^{nd} jet", "VectorialSump-b-jet2", "Vectorial sum of tracks momenta (GeV/c)", 0., 0.);
  
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "1^{st} jet", "Highestp-b-jet1", "Highest track momentum (GeV/c)", 0., 0.);
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "2^{nd} jet", "Highestp-b-jet2", "Highest track momentum (GeV/c)", 0., 0.);
  
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "1^{st} jet", "Highestp-b-jet1", "Scalar sum of the 3 highest tracks momenta (GeV/c)", 0., 0.);
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "2^{nd} jet", "Highestp-b-jet2", "Scalar sum of the 3 highest tracks momenta (GeV/c)", 0., 0.);
  
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "1^{st} jet", "R1-b-jet1", "#frac{Highest track momentum}{Scalar sum of tracks momenta}", 0., 0.);
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "2^{nd} jet", "R1-b-jet2", "#frac{Highest track momentum}{Scalar sum of tracks momenta}", 0., 0.);
  
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "1^{st} jet", "R3-b-jet1", "#frac{Scalar sum of the 3 highest tracks momenta}{Scalar sum of tracks momenta}", 0., 0.);
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "2^{nd} jet", "R3-b-jet2", "#frac{Scalar sum of the 3 highest tracks momenta}{Scalar sum of tracks momenta}", 0., 0.);
  
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "1^{st} jet", "D0Mass-b-jet1", "m(D^{0}#rightarrow#kappa^{#pm}#pi^{#mp}) (GeV/c^{2})", 0., 0.);
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "2^{nd} jet", "D0Mass-b-jet2", "m(D^{0}#rightarrow#kappa^{#pm}#pi^{#mp}) (GeV/c^{2})", 0., 0.);
  
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "1^{st} jet", "D0p-b-jet1", "p(D^{0}#rightarrow#kappa^{#pm}#pi^{#mp}) (GeV/c)", 0., 0.);
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "2^{nd} jet", "D0p-b-jet2", "p(D^{0}#rightarrow#kappa^{#pm}#pi^{#mp}) (GeV/c)", 0., 0.);

  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "1^{st} jet", "BMomentum-b-nobias-jet1", "p(#kappa^{#pm}#pi^{#mp}+#mu^{#mp}) (GeV/c)", 0., 0.);
  plotHisto(inBatch, date, isCSVbased, fi_data, fi_sl, fi_dl, "2^{nd} jet", "BMomentum-b-nobias-jet2", "p(#kappa^{#pm}#pi^{#mp}+#mu^{#mp}) (GeV/c)", 0., 0.);

  fi_data->Close();
  fi_sl->Close();
  fi_dl->Close();
  delete fi_data;
  delete fi_sl;
  delete fi_dl;
}
  
//---------------------------------------------------------------
int doRivetJob(bool inBatch = true, TString date = "")
//---------------------------------------------------------------
{  
  if (date.Length() > 0)  {
    gROOT->ProcessLine(".! mkdir "+date);
    gROOT->ProcessLine(".! mkdir "+date+"/csv");
    gROOT->ProcessLine(".! mkdir "+date+"/pT");

    doRivetJob_cut(inBatch, date, true);
    doRivetJob_cut(inBatch, date, false);
    return 0;
  }
  else 
    return 1;
}

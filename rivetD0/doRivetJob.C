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
void plotHisto(bool inBatch, TString date, TString type, TFile* fi_data, TFile* fi_sl, TFile* fi_dl, TFile* fi_ts, TFile* fi_tt, TFile* fi_tw, TFile* fi_Ts, TFile* fi_Tt, TFile* fi_Tw, TFile* fi_w1, TFile* fi_w2, TFile* fi_w3, TFile* fi_w4, TFile* fi_ww, TFile* fi_wz, TFile* fi_zz, TString channel, TString h_name, TString xtitle, double xmin, double xmax, TString position)
//---------------------------------------------------------------
{
  TStyle* my_style = createMyStyle();
  my_style->cd();
  if (inBatch) gROOT->SetBatch(true);
  
  //TLatex* channel_tex = new TLatex(0.87, 0.9, channel);
  TLatex* channel_tex = new TLatex(0.8, 0.9, channel);
  channel_tex->SetNDC(true);
  channel_tex->SetTextFont(43);
  channel_tex->SetTextSize(TITLE_FONTSIZE - 6);  

  TString dir_name = "basedSelection/";
  TString rep_name = date;
  if (type.Contains("csv", TString::kIgnoreCase)) {
    dir_name = "CSV" + dir_name;
    rep_name = rep_name + "/csv/";
  }
  else if (type.Contains("pT", TString::kIgnoreCase)) { 
    dir_name = "pT" + dir_name;
    rep_name = rep_name + "/pT/";
  }
  else if (type.Contains("muTag", TString::kIgnoreCase)) {
    dir_name = "muTagBasedSelection/";
    rep_name = rep_name + "/muTag/";
  }
  
  TH1D* h_data = (TH1D*)fi_data->Get(dir_name+h_name);
  TH1D* h_sl = (TH1D*)fi_sl->Get(dir_name+h_name);
  TH1D* h_dl = (TH1D*)fi_dl->Get(dir_name+h_name);
  TH1D* h_ts = (TH1D*)fi_ts->Get(dir_name+h_name);
  TH1D* h_tt = (TH1D*)fi_tt->Get(dir_name+h_name);
  TH1D* h_tw = (TH1D*)fi_tw->Get(dir_name+h_name);
  TH1D* h_Ts = (TH1D*)fi_Ts->Get(dir_name+h_name);
  TH1D* h_Tt = (TH1D*)fi_Tt->Get(dir_name+h_name);
  TH1D* h_Tw = (TH1D*)fi_Tw->Get(dir_name+h_name);
  TH1D* h_w1 = (TH1D*)fi_w1->Get(dir_name+h_name);
  TH1D* h_w2 = (TH1D*)fi_w2->Get(dir_name+h_name);
  TH1D* h_w3 = (TH1D*)fi_w3->Get(dir_name+h_name);
  TH1D* h_w4 = (TH1D*)fi_w4->Get(dir_name+h_name);
  TH1D* h_ww = (TH1D*)fi_ww->Get(dir_name+h_name);
  TH1D* h_wz = (TH1D*)fi_wz->Get(dir_name+h_name);
  TH1D* h_zz = (TH1D*)fi_zz->Get(dir_name+h_name);
  h_sl->Scale(2.*25.8031*19769./12031276.);
  h_dl->Scale(2.*107.6722*19769./25339818.);
  h_ts->Scale(2.*3.79*19769./259961.);
  h_tt->Scale(2.*56.4*19769./3758227.);
  h_tw->Scale(2.*11.1*19769./497658.);
  h_Ts->Scale(2.*1.76*19769./139974.);
  h_Tt->Scale(2.*30.7*19769./1935072.);
  h_Tw->Scale(2.*11.1*19769./493460.);
  h_w1->Scale(2.*5562.36*19769./23138212.);
  h_w2->Scale(2.*1802.62*19769./33892560.);
  h_w3->Scale(2.*534.604*19769./15464503.);
  h_w4->Scale(2.*220.434*19769./13382803.);
  h_ww->Scale(2.*54.838*19769./10000431.);
  h_wz->Scale(2.*33.21*19769./10000283.);
  h_zz->Scale(2.*17.654*19769./9799908.);
  TH1D* h_mc = (TH1D*)h_sl->Clone();
  h_mc->Add(h_dl);
  h_mc->Add(h_ts);
  h_mc->Add(h_tt);
  h_mc->Add(h_tw);
  h_mc->Add(h_Ts);
  h_mc->Add(h_Tt);
  h_mc->Add(h_Tw);
  h_mc->Add(h_ww);
  h_mc->Add(h_wz);
  h_mc->Add(h_zz);
  h_mc->Add(h_w1);
  h_mc->Add(h_w2);
  h_mc->Add(h_w3);
  h_mc->Add(h_w4);
  
  // because we miss some data...
  h_sl->Scale(h_data->Integral()/h_mc->Integral());
  h_dl->Scale(h_data->Integral()/h_mc->Integral());
  h_ts->Scale(h_data->Integral()/h_mc->Integral());
  h_tt->Scale(h_data->Integral()/h_mc->Integral());
  h_tw->Scale(h_data->Integral()/h_mc->Integral());
  h_Ts->Scale(h_data->Integral()/h_mc->Integral());
  h_Tt->Scale(h_data->Integral()/h_mc->Integral());
  h_Tw->Scale(h_data->Integral()/h_mc->Integral());
  h_w1->Scale(h_data->Integral()/h_mc->Integral());
  h_w2->Scale(h_data->Integral()/h_mc->Integral());
  h_w3->Scale(h_data->Integral()/h_mc->Integral());
  h_w4->Scale(h_data->Integral()/h_mc->Integral());
  h_ww->Scale(h_data->Integral()/h_mc->Integral());
  h_wz->Scale(h_data->Integral()/h_mc->Integral());
  h_zz->Scale(h_data->Integral()/h_mc->Integral());
  h_mc->Scale(h_data->Integral()/h_mc->Integral());
 
  h_style(h_data, 1, 1, 4000, -1111., -1111., 510, 510, 21, 1, 1.2, 0, xtitle);
  h_style(h_sl, 2, 2, 1001, -1111., -1111., 510, 510, 21, 2, 1.2, 0, xtitle);
  h_style(h_dl, 3, 3, 1001, -1111., -1111., 510, 510, 21, 3, 1.2, 0, xtitle);
  h_style(h_ts, 28, 28, 1001, -1111., -1111., 510, 510, 21, 28, 1.2, 0, xtitle);
  h_style(h_tt, 28, 28, 1001, -1111., -1111., 510, 510, 21, 28, 1.2, 0, xtitle);
  h_style(h_tw, 28, 28, 1001, -1111., -1111., 510, 510, 21, 28, 1.2, 0, xtitle);
  h_style(h_Ts, 28, 28, 1001, -1111., -1111., 510, 510, 21, 28, 1.2, 0, xtitle);
  h_style(h_Tt, 28, 28, 1001, -1111., -1111., 510, 510, 21, 28, 1.2, 0, xtitle);
  h_style(h_Tw, 28, 28, 1001, -1111., -1111., 510, 510, 21, 28, 1.2, 0, xtitle);
  h_style(h_w1, 33, 33, 1001, -1111., -1111., 510, 510, 21, 33, 1.2, 0, xtitle);
  h_style(h_w2, 33, 33, 1001, -1111., -1111., 510, 510, 21, 33, 1.2, 0, xtitle);
  h_style(h_w3, 33, 33, 1001, -1111., -1111., 510, 510, 21, 33, 1.2, 0, xtitle);
  h_style(h_w4, 33, 33, 1001, -1111., -1111., 510, 510, 21, 33, 1.2, 0, xtitle);
  h_style(h_ww, 5, 5, 1001, -1111., -1111., 510, 510, 21, 5, 1.2, 0, xtitle);
  h_style(h_wz, 5, 5, 1001, -1111., -1111., 510, 510, 21, 5, 1.2, 0, xtitle);
  h_style(h_zz, 5, 5, 1001, -1111., -1111., 510, 510, 21, 5, 1.2, 0, xtitle);
  h_style(h_mc, 0, 12, 3004, -1111., -1111., 510, 510, 1, 12, 1, 0, xtitle);

  THStack* hs_mc = new THStack("hs", "hs");
  hs_mc->Add(h_w1);
  hs_mc->Add(h_w2);
  hs_mc->Add(h_w3);
  hs_mc->Add(h_w4);
  hs_mc->Add(h_ww);
  hs_mc->Add(h_wz);
  hs_mc->Add(h_zz);
  hs_mc->Add(h_ts);
  hs_mc->Add(h_tt);
  hs_mc->Add(h_tw);
  hs_mc->Add(h_Ts);
  hs_mc->Add(h_Tt);
  hs_mc->Add(h_Tw);
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
  if (position.Contains("right", TString::kIgnoreCase)) {
    TLegend *leg = new TLegend(0.64,0.52,0.93,0.89);
    leg_style(leg, 12);
    leg->AddEntry(h_data,"Data - Run 2012 A,B,C,D","LP");
    leg->AddEntry(h_sl,"MG+PY6 Z2* Semilept. t#bar{t}","F");
    leg->AddEntry(h_dl,"MG+PY6 Z2* Dilept. t#bar{t}","F");
    leg->AddEntry(h_ts,"Single-top","F");
    leg->AddEntry(h_ww,"WW,WZ,ZZ","F");
    leg->AddEntry(h_w1,"W#rightarrow l#nu + jets","F");
    leg->AddEntry(h_mc,"Stat. uncertainty","F");
    leg->Draw();
  }
  else if (position.Contains("short", TString::kIgnoreCase)) {
    TLegend *leg = new TLegend(0.78,0.52,0.93,0.89);
    leg_style(leg, 12);
    leg->AddEntry(h_data,"Data","LP");
    leg->AddEntry(h_sl,"Semilept. t#bar{t}","F");
    leg->AddEntry(h_dl,"Dilept. t#bar{t}","F");
    leg->AddEntry(h_ts,"Single-top","F");
    leg->AddEntry(h_ww,"WW,WZ,ZZ","F");
    leg->AddEntry(h_w1,"W#rightarrow l#nu + jets","F");
    leg->AddEntry(h_mc,"Stat. uncertainty","F");
    leg->Draw();
  }
  else if (position.Contains("left", TString::kIgnoreCase)) {
    TLegend *leg = new TLegend(0.19,0.52,0.47,0.89);
    leg_style(leg, 12);
    leg->AddEntry(h_data,"Data - Run 2012 A,B,C,D","LP");
    leg->AddEntry(h_sl,"MG+PY6 Z2* Semilept. t#bar{t}","F");
    leg->AddEntry(h_dl,"MG+PY6 Z2* Dilept. t#bar{t}","F");
    leg->AddEntry(h_ts,"Single-top","F");
    leg->AddEntry(h_ww,"WW,WZ,ZZ","F");
    leg->AddEntry(h_w1,"W#rightarrow l#nu + jets","F");
    leg->AddEntry(h_mc,"Stat. uncertainty","F");
    leg->Draw();
  }  
  
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
  delete h_ts;
  delete h_tt;
  delete h_tw;
  delete h_Ts;
  delete h_Tt;
  delete h_Tw;
  delete h_w1;
  delete h_w2;
  delete h_w3;
  delete h_w4;
  delete h_ww;
  delete h_wz;
  delete h_zz;
  delete channel_tex;
  delete my_style;
}

//---------------------------------------------------------------
void doRivetJob_d0(bool inBatch, TString date, bool isCSVbased)
//---------------------------------------------------------------
{
  TString fi_data_name = "../test/crab_results/"+date;
  TString fi_sl_name = "../test/crab_results/"+date;
  TString fi_dl_name = "../test/crab_results/"+date;
  TString fi_ts_name = "../test/crab_results/"+date;
  TString fi_tt_name = "../test/crab_results/"+date;
  TString fi_tw_name = "../test/crab_results/"+date;
  TString fi_Ts_name = "../test/crab_results/"+date;
  TString fi_Tt_name = "../test/crab_results/"+date;
  TString fi_Tw_name = "../test/crab_results/"+date;
  TString fi_w1_name = "../test/crab_results/"+date;
  TString fi_w2_name = "../test/crab_results/"+date;
  TString fi_w3_name = "../test/crab_results/"+date;
  TString fi_w4_name = "../test/crab_results/"+date;
  TString fi_ww_name = "../test/crab_results/"+date;
  TString fi_wz_name = "../test/crab_results/"+date;
  TString fi_zz_name = "../test/crab_results/"+date;
  TString type;
  /*
  std::cout << "\t\t /!\\ Don't forget in " << fi_data_name << " /!\\ " << std::endl;
  if (isCSVbased)
    std::cout << "hadd D0ForRivet_csv_Data_merged.root D0ForRivet_csv_El_merged.root D0ForRivet_csv_Mu_merged.root" << std::endl;
  else 
    std::cout << "hadd D0ForRivet_pT_Data_merged.root D0ForRivet_pT_El_merged.root D0ForRivet_pT_Mu_merged.root" << std::endl;
  getchar();
  */
  if (isCSVbased) {
    fi_data_name = fi_data_name + "/D0ForRivet_csv_Data_merged.root";
    fi_sl_name = fi_sl_name + "/D0ForRivet_csv_TTJets_SemiLeptMGDecays.root";
    fi_dl_name = fi_dl_name + "/D0ForRivet_csv_TTJets_FullLeptMGDecays.root";
    fi_ts_name = fi_ts_name + "/D0ForRivet_csv_T_s-channel.root";
    fi_tt_name = fi_tt_name + "/D0ForRivet_csv_T_t-channel.root";
    fi_tw_name = fi_tw_name + "/D0ForRivet_csv_T_tW-channel.root";
    fi_Ts_name = fi_Ts_name + "/D0ForRivet_csv_Tbar_s-channel.root";
    fi_Tt_name = fi_Tt_name + "/D0ForRivet_csv_Tbar_t-channel.root";
    fi_Tw_name = fi_Tw_name + "/D0ForRivet_csv_Tbar_tW-channel.root";
    fi_w1_name = fi_w1_name + "/D0ForRivet_csv_W1JetsToLNu.root";
    fi_w2_name = fi_w2_name + "/D0ForRivet_csv_W2JetsToLNu.root";
    fi_w3_name = fi_w3_name + "/D0ForRivet_csv_W3JetsToLNu.root";
    fi_w4_name = fi_w4_name + "/D0ForRivet_csv_W4JetsToLNu.root";
    fi_ww_name = fi_ww_name + "/D0ForRivet_csv_WW-incl.root";
    fi_wz_name = fi_wz_name + "/D0ForRivet_csv_WZ-incl.root";
    fi_zz_name = fi_zz_name + "/D0ForRivet_csv_ZZ-incl.root";
    type = "csv";
  }
  else {
    fi_data_name = fi_data_name + "/D0ForRivet_pT_Data_merged.root";
    fi_sl_name = fi_sl_name + "/D0ForRivet_pT_TTJets_SemiLeptMGDecays.root";
    fi_dl_name = fi_dl_name + "/D0ForRivet_pT_TTJets_FullLeptMGDecays.root";
    fi_ts_name = fi_ts_name + "/D0ForRivet_pT_T_s-channel.root";
    fi_tt_name = fi_tt_name + "/D0ForRivet_pT_T_t-channel.root";
    fi_tw_name = fi_tw_name + "/D0ForRivet_pT_T_tW-channel.root";
    fi_Ts_name = fi_Ts_name + "/D0ForRivet_pT_Tbar_s-channel.root";
    fi_Tt_name = fi_Tt_name + "/D0ForRivet_pT_Tbar_t-channel.root";
    fi_Tw_name = fi_Tw_name + "/D0ForRivet_pT_Tbar_tW-channel.root";
    fi_w1_name = fi_w1_name + "/D0ForRivet_pT_W1JetsToLNu.root";
    fi_w2_name = fi_w2_name + "/D0ForRivet_pT_W2JetsToLNu.root";
    fi_w3_name = fi_w3_name + "/D0ForRivet_pT_W3JetsToLNu.root";
    fi_w4_name = fi_w4_name + "/D0ForRivet_pT_W4JetsToLNu.root";
    fi_ww_name = fi_ww_name + "/D0ForRivet_pT_WW-incl.root";
    fi_wz_name = fi_wz_name + "/D0ForRivet_pT_WZ-incl.root";
    fi_zz_name = fi_zz_name + "/D0ForRivet_pT_ZZ-incl.root";
    type = "pT";
  }
  
  TFile* fi_data = TFile::Open(fi_data_name);
  TFile* fi_sl = TFile::Open(fi_sl_name);
  TFile* fi_dl = TFile::Open(fi_dl_name);
  TFile* fi_ts = TFile::Open(fi_ts_name);
  TFile* fi_tt = TFile::Open(fi_tt_name);
  TFile* fi_tw = TFile::Open(fi_tw_name);
  TFile* fi_Ts = TFile::Open(fi_Ts_name);
  TFile* fi_Tt = TFile::Open(fi_Tt_name);
  TFile* fi_Tw = TFile::Open(fi_Tw_name);
  TFile* fi_w1 = TFile::Open(fi_w1_name);
  TFile* fi_w2 = TFile::Open(fi_w2_name);
  TFile* fi_w3 = TFile::Open(fi_w3_name);
  TFile* fi_w4 = TFile::Open(fi_w4_name);
  TFile* fi_ww = TFile::Open(fi_ww_name);
  TFile* fi_wz = TFile::Open(fi_wz_name);
  TFile* fi_zz = TFile::Open(fi_zz_name);

  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "All", "NPrimaryVtx", "Primary vertices multiplicity", 0., 0., "right");
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "All", "NJets", "Jets multiplicity", 0., 0., "right");
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "All", "TransverseMomentumJets", "p(jets) (GeV/c)", 0., 250., "right");
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "All", "NCSVJets", "b-tagged jets multiplicity", 0., 0., "right");

  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Both jets", "CSV-b-jets", "CSV discriminant", 0., 0., "right");
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Both jets", "TransverseMomentum-b-jets", "p(jets) (GeV/c)", 0., 300., "right");

  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "1^{st} jet", "Nch-b-jet1", "Tracks multiplicity", 0., 20., "right");
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "2^{nd} jet", "Nch-b-jet2", "Tracks multiplicity", 0., 20., "right");

  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "1^{st} jet", "Sump-b-jet1", "Scalar sum of tracks momenta (GeV/c)", 0., 400., "right");
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "2^{nd} jet", "Sump-b-jet2", "Scalar sum of tracks momenta (GeV/c)", 0., 400., "right");
  
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "1^{st} jet", "VectorialSump-b-jet1", "Vectorial sum of tracks momenta (GeV/c)", 0., 400., "right");
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "2^{nd} jet", "VectorialSump-b-jet2", "Vectorial sum of tracks momenta (GeV/c)", 0., 400., "right");
  
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "1^{st} jet", "Highestp-b-jet1", "Highest track momentum (GeV/c)", 0., 200., "right");
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "2^{nd} jet", "Highestp-b-jet2", "Highest track momentum (GeV/c)", 0., 200., "right");
  
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "1^{st} jet", "Sum3p-b-jet1", "Scalar sum of the 3 highest tracks momenta (GeV/c)", 0., 0., "right");
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "2^{nd} jet", "Sum3p-b-jet2", "Scalar sum of the 3 highest tracks momenta (GeV/c)", 0., 0., "right");
  
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "1^{st} jet", "R1-b-jet1", "R_{1}", 0., 0., "right");
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "2^{nd} jet", "R1-b-jet2", "R_{1}", 0., 0., "right");
  
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "1^{st} jet", "R3-b-jet1", "R_{3}", 0., 0., "right");
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "2^{nd} jet", "R3-b-jet2", "R_{3}", 0., 0., "right");
  
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "1^{st} jet", "D0Mass-b-jet1", "m(D^{0}#rightarrow#kappa^{+}#pi^{-}) (GeV/c^{2})", 0., 0., "right");
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "2^{nd} jet", "D0Mass-b-jet2", "m(D^{0}#rightarrow#kappa^{+}#pi^{-}) (GeV/c^{2})", 0., 0., "right");
  
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "1^{st} jet", "D0p-b-jet1", "p(D^{0}#rightarrow#kappa^{+}#pi^{-}) (GeV/c)", 0., 0., "right");
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "2^{nd} jet", "D0p-b-jet2", "p(D^{0}#rightarrow#kappa^{+}#pi^{-}) (GeV/c)", 0., 0., "right");

  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "1^{st} jet", "BMomentum-nobias-b-jet1", "p(#kappa^{+}#pi^{-}+#mu^{-}) (GeV/c)", 0., 350., "right");
  plotHisto(inBatch, date, type, fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "2^{nd} jet", "BMomentum-nobias-b-jet2", "p(#kappa^{+}#pi^{-}+#mu^{-}) (GeV/c)", 0., 350., "right");

  fi_data->Close();
  fi_sl->Close();
  fi_dl->Close();
  fi_ts->Close();
  fi_tt->Close();
  fi_tw->Close();
  fi_Ts->Close();
  fi_Tt->Close();
  fi_Tw->Close();
  fi_w1->Close();
  fi_w2->Close();
  fi_w3->Close();
  fi_w4->Close();
  fi_ww->Close();
  fi_wz->Close();
  fi_zz->Close();
  delete fi_data;
  delete fi_sl;
  delete fi_dl;
  delete fi_ts;
  delete fi_tt;
  delete fi_tw;
  delete fi_Ts;
  delete fi_Tt;
  delete fi_Tw;
  delete fi_w1;
  delete fi_w2;
  delete fi_w3;
  delete fi_w4;
  delete fi_ww;
  delete fi_wz;
  delete fi_zz;
}

//---------------------------------------------------------------
void doRivetJob_mutag(bool inBatch, TString date)
//---------------------------------------------------------------
{
  TString fi_data_name = "../test/crab_results/"+date+"/MuTagForRivet_Data_merged.root";
  TString fi_sl_name = "../test/crab_results/"+date+"/MuTagForRivet_TTJets_SemiLeptMGDecays.root";
  TString fi_dl_name = "../test/crab_results/"+date+"/MuTagForRivet_TTJets_FullLeptMGDecays.root";
  TString fi_ts_name = "../test/crab_results/"+date+"/MuTagForRivet_T_s-channel.root";
  TString fi_tt_name = "../test/crab_results/"+date+"/MuTagForRivet_T_t-channel.root";
  TString fi_tw_name = "../test/crab_results/"+date+"/MuTagForRivet_T_tW-channel.root";
  TString fi_Ts_name = "../test/crab_results/"+date+"/MuTagForRivet_Tbar_s-channel.root";
  TString fi_Tt_name = "../test/crab_results/"+date+"/MuTagForRivet_Tbar_t-channel.root";
  TString fi_Tw_name = "../test/crab_results/"+date+"/MuTagForRivet_Tbar_tW-channel.root";
  TString fi_w1_name = "../test/crab_results/"+date+"/MuTagForRivet_W1JetsToLNu.root";
  TString fi_w2_name = "../test/crab_results/"+date+"/MuTagForRivet_W2JetsToLNu.root";
  TString fi_w3_name = "../test/crab_results/"+date+"/MuTagForRivet_W3JetsToLNu.root";
  TString fi_w4_name = "../test/crab_results/"+date+"/MuTagForRivet_W4JetsToLNu.root";
  TString fi_ww_name = "../test/crab_results/"+date+"/MuTagForRivet_WW-incl.root";
  TString fi_wz_name = "../test/crab_results/"+date+"/MuTagForRivet_WZ-incl.root";
  TString fi_zz_name = "../test/crab_results/"+date+"/MuTagForRivet_ZZ-incl.root";
  
  TFile* fi_data = TFile::Open(fi_data_name);
  TFile* fi_sl = TFile::Open(fi_sl_name);
  TFile* fi_dl = TFile::Open(fi_dl_name);
  TFile* fi_ts = TFile::Open(fi_ts_name);
  TFile* fi_tt = TFile::Open(fi_tt_name);
  TFile* fi_tw = TFile::Open(fi_tw_name);
  TFile* fi_Ts = TFile::Open(fi_Ts_name);
  TFile* fi_Tt = TFile::Open(fi_Tt_name);
  TFile* fi_Tw = TFile::Open(fi_Tw_name);
  TFile* fi_w1 = TFile::Open(fi_w1_name);
  TFile* fi_w2 = TFile::Open(fi_w2_name);
  TFile* fi_w3 = TFile::Open(fi_w3_name);
  TFile* fi_w4 = TFile::Open(fi_w4_name);
  TFile* fi_ww = TFile::Open(fi_ww_name);
  TFile* fi_wz = TFile::Open(fi_wz_name);
  TFile* fi_zz = TFile::Open(fi_zz_name);

  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "All", "NPrimaryVtx", "Primary vertices multiplicity", 0., 0., "right");

  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Tagged jets", "CSV-b-jets", "CSV discriminant", 0., 0., "right");
  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Tagged jets", "TransverseMomentum-b-jets", "p(jets) (GeV/c)", 0., 300., "right");
  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Tagged jets", "Nch-b-jets", "Tracks multiplicity", 0., 0., "right");
  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Tagged jets", "Sump-b-jets", "Scalar sum of tracks momenta (GeV/c)", 0., 400., "right");
  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Tagged jets", "VectorialSump-b-jets", "Vectorial sum of tracks momenta (GeV/c)", 0., 400., "right");
  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Tagged jets", "Highestp-b-jets", "Highest track momentum (GeV/c)", 0., 200., "right");
  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Tagged jets", "Sum2p-b-jets", "Scalar sum of the 2 highest tracks momenta (GeV/c)", 0., 0., "right");
  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Tagged jets", "Sum3p-b-jets", "Scalar sum of the 3 highest tracks momenta (GeV/c)", 0., 0., "right");
  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Tagged jets", "R1-b-jets", "R_{1}", 0., 0., "right");
  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Tagged jets", "R2-b-jets", "R_{2}", 1.5, 1.5, "short");
  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Tagged jets", "R3-b-jets", "R_{3}", 1.5, 1.5, "left");
  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Tagged jets", "D0Mass-b-jets", "m(D^{0}#rightarrow#kappa^{+}#pi^{-}) (GeV/c^{2})", 0., 0., "right");
  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Tagged jets", "D0p-b-jets", "p(D^{0}#rightarrow#kappa^{+}#pi^{-}) (GeV/c)", 0., 0., "right");
  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Tagged jets", "BMomentum-nobias-b-jets", "p(#kappa^{+}#pi^{-}+#mu^{-}) (GeV/c)", 0., 350., "right");
  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Tagged jets", "BMass-nobias-b-jets", "M(#kappa^{+}#pi^{-}+#mu^{-}) (GeV/c^{2})", 0., 0, "right");
  plotHisto(inBatch, date, "muTag", fi_data, fi_sl, fi_dl, fi_ts, fi_tt, fi_tw, fi_Ts, fi_Tt, fi_Tw, fi_w1, fi_w2, fi_w3, fi_w4, fi_ww, fi_wz, fi_zz, "Tagged jets", "Muonp-nobias-b-jets", "p(#mu) (GeV/c)", 0., 0, "right");

  fi_data->Close();
  fi_sl->Close();
  fi_dl->Close();
  fi_ts->Close();
  fi_tt->Close();
  fi_tw->Close();
  fi_Ts->Close();
  fi_Tt->Close();
  fi_Tw->Close();
  fi_w1->Close();
  fi_w2->Close();
  fi_w3->Close();
  fi_w4->Close();
  fi_ww->Close();
  fi_wz->Close();
  fi_zz->Close();
  delete fi_data;
  delete fi_sl;
  delete fi_dl;
  delete fi_ts;
  delete fi_tt;
  delete fi_tw;
  delete fi_Ts;
  delete fi_Tt;
  delete fi_Tw;
  delete fi_w1;
  delete fi_w2;
  delete fi_w3;
  delete fi_w4;
  delete fi_ww;
  delete fi_wz;
  delete fi_zz;
}
  
//---------------------------------------------------------------
int doRivetJob(TString date = "", TString type = "muTag", bool inBatch = true)
//---------------------------------------------------------------
{  
  if (date.Length() > 0)  {
    gROOT->ProcessLine(".! mkdir "+date);

    if (type.Contains("D0", TString::kIgnoreCase)) {
      gROOT->ProcessLine(".! mkdir "+date+"/csv");
      gROOT->ProcessLine(".! mkdir "+date+"/pT");
      doRivetJob_d0(inBatch, date, true);
      doRivetJob_d0(inBatch, date, false);
    }
    else if (type.Contains("MuTag", TString::kIgnoreCase)) {
      gROOT->ProcessLine(".! mkdir "+date+"/muTag");
      doRivetJob_mutag(inBatch, date);
    }
    return 0;
  }
  else 
    return 1;
}

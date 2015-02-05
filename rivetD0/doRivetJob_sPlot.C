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
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TTree.h"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooStats/SPlot.h"

#pragma once

#define TITLE_FONTSIZE 26
#define LABEL_FONTSIZE 18

#define LEFT_MARGIN 0.17
#define RIGHT_MARGIN 0.03
#define TOP_MARGIN 0.05
#define BOTTOM_MARGIN 0.13

using namespace RooFit;

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

void graphasymerror_mystyle(TGraphAsymmErrors *gr,
                          TString name="",
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
                          float marker_size=1.2) {
  
  gr->SetName(name);
  gr->SetTitle("");
  
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

  gr->GetYaxis()->SetTitleOffset(2.5);
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
TGraphAsymmErrors **treatHisto(bool inBatch, TStyle* my_style, TString date, bool isCSVbased, const int NsPlots, TFile* fi, TString name, TString channel)
//---------------------------------------------------------------
{
  my_style->cd();
  if (inBatch) gROOT->SetBatch(true);
  
  TLatex* channel_tex_r = new TLatex(0.66, 0.9, channel);
  channel_tex_r->SetNDC(true);
  channel_tex_r->SetTextFont(43);
  channel_tex_r->SetTextSize(TITLE_FONTSIZE - 6);  
  
  TLatex* channel_tex_l = new TLatex(0.22, 0.9, channel);
  channel_tex_l->SetNDC(true);
  channel_tex_l->SetTextFont(43);
  channel_tex_l->SetTextSize(TITLE_FONTSIZE - 6);  
  
  TString dir_name = "basedSelection/";
  TString rep_name = date;
  if (isCSVbased) {
    dir_name = "CSV" + dir_name;
    rep_name = rep_name + "/csv/sPlot/";
  }
  else { 
    dir_name = "pT" + dir_name;
    rep_name = rep_name + "/pT/sPlot/";
  }
  
  // Import the tree
  TTree *tree = (TTree*) fi->Get(dir_name+"/D0window-b-jets");
  RooRealVar weight("Weight", "Weight", 0., 2.);
  RooRealVar csvdisc("CSVdisc", "CSV discriminant", 0., 1.);
  RooRealVar d0mass("D0mass", "D^{0} mass", 1.7, 2., "GeV/c^{2}");
  RooRealVar bmomentum("Bmomentum", "p(#kappa^{+}#pi^{-}+#mu^{-})", 0., 400., "GeV/c");
  RooRealVar r1("R1", "R1", 0., 1.);
  RooRealVar r3("R3", "R3", 0., 1.);
  RooRealVar nch("Nch", "Tracks multiplicity", 0., 30.);
  RooRealVar sumpT("SumpT", "#sum p_{T}", 0., 1000., "GeV/c");
  RooRealVar averpT("AveragepT", "#LT pT #GT", 0., 100., "GeV/c");
  RooRealVar r1_nomu("R1_nomu", "R1 (no #mu)", 0., 1.);
  RooRealVar r3_nomu("R3_nomu", "R3 (no #mu)", 0., 1.);
  RooRealVar genId("GenId", "True jets flavour", 0., 22.);
  RooArgSet variables(weight, csvdisc, d0mass, bmomentum, r1, r3, nch, sumpT, averpT);
  variables.add(r1_nomu);
  variables.add(r3_nomu);
  if (!name.Contains("Data")) 
    variables.add(genId);
  RooDataSet *data = new RooDataSet("data", "data", variables, Import(*tree), WeightVar(weight));
  
  // Fit the D0 mass peak
  RooRealVar nsig("nsig", "number of signal events", 0.1*tree->GetEntries(), 0., tree->GetEntries());
  RooRealVar nbck("nbck", "number of background events", 0.9*tree->GetEntries(), 0., tree->GetEntries());
  RooRealVar mean("mean", "mean", 1.86484, 1.84554, 1.88414);
  RooRealVar sigma("sigma", "sigma", 0.0193, 0.01, 0.03);
  RooRealVar lambda("lambda", "slope", -1.8, -2.2, -0.7);
  RooGaussian sig("sig", "gaussian signal", d0mass, mean, sigma);
  RooExponential bck("bck", "exponential background", d0mass, lambda);
  RooAddPdf *model = new RooAddPdf("model", "model", RooArgList(sig,bck), RooArgList(nsig,nbck));
  
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(variables);
  ws->import(*data);
  ws->import(*model);
  
  model->fitTo(*data);
  RooPlot* fr_discvar = d0mass.frame();
  data->plotOn(fr_discvar, Name("Data"), Binning(30));
  model->plotOn(fr_discvar, FillStyle(0), MoveToBack(), Name ("total fit")); // LineColor(9)
  model->plotOn(fr_discvar, Components(bck), LineColor(1), LineWidth(2), FillStyle(1001), FillColor(920), DrawOption("LF"), MoveToBack(), Name("bck fit"));
  // model->plotOn(fr_discvar, Components(bck), LineColor(kBlue), Name("bck fit"));
  // model->plotOn(fr_discvar, Components(RooArgSet(sig)), LineColor(kRed), Name("signal fit"));
  
  // Apply the sPlot technique
  double nsigVal = nsig.getVal();
  double nsigErr = nsig.getError();
  double nbckVal = nbck.getVal();
  double nbckErr = nbck.getError();
  mean.setConstant();
  sigma.setConstant();
  lambda.setConstant();
  
  RooStats::SPlot* sData = new RooStats::SPlot("sData", "SPlot for m(D0)", *data, model, RooArgList(nsig, nbck));
  ws->import(*data, Rename("dataWithSWeights"));
  data = (RooDataSet*)ws->data("dataWithSWeights");
  RooDataSet *sigData = new RooDataSet(data->GetName(), data->GetTitle(), data, *data->get(), 0, "nsig_sw");
  RooDataSet *bckData = new RooDataSet(data->GetName(), data->GetTitle(), data, *data->get(), 0, "nbck_sw");
  
  RooPlot* fr_Bmomentum = bmomentum.frame();
  sigData->plotOn(fr_Bmomentum, Name("sigData_Bmomentum"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kRed), LineWidth(2), MarkerColor(kRed), MarkerStyle(20));
  bckData->plotOn(fr_Bmomentum, Name("bckData_Bmomentum"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kBlue), LineWidth(2), MarkerColor(kBlue), MarkerStyle(24));
  
  RooPlot* fr_CSVdisc = csvdisc.frame();
  sigData->plotOn(fr_CSVdisc, Name("sigData_CSVdisc"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kRed), LineWidth(2), MarkerColor(kRed), MarkerStyle(20));
  bckData->plotOn(fr_CSVdisc, Name("bckData_CSVdisc"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kBlue), LineWidth(2), MarkerColor(kBlue), MarkerStyle(24));
  
  RooPlot* fr_R1 = r1.frame();
  sigData->plotOn(fr_R1, Name("sigData_R1"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kRed), LineWidth(2), MarkerColor(kRed), MarkerStyle(20));
  bckData->plotOn(fr_R1, Name("bckData_R1"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kBlue), LineWidth(2), MarkerColor(kBlue), MarkerStyle(24));
  
  RooPlot* fr_R3 = r3.frame();
  sigData->plotOn(fr_R3, Name("sigData_R3"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kRed), LineWidth(2), MarkerColor(kRed), MarkerStyle(20));
  bckData->plotOn(fr_R3, Name("bckData_R3"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kBlue), LineWidth(2), MarkerColor(kBlue), MarkerStyle(24));
  
  RooPlot* fr_Nch = nch.frame();
  sigData->plotOn(fr_Nch, Name("sigData_Nch"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kRed), LineWidth(2), MarkerColor(kRed), MarkerStyle(20));
  bckData->plotOn(fr_Nch, Name("bckData_Nch"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kBlue), LineWidth(2), MarkerColor(kBlue), MarkerStyle(24));
  
  RooPlot* fr_SumpT = sumpT.frame();
  sigData->plotOn(fr_SumpT, Name("sigData_SumpT"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kRed), LineWidth(2), MarkerColor(kRed), MarkerStyle(20));
  bckData->plotOn(fr_SumpT, Name("bckData_SumpT"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kBlue), LineWidth(2), MarkerColor(kBlue), MarkerStyle(24));
  
  RooPlot* fr_AveragepT = averpT.frame();
  sigData->plotOn(fr_AveragepT, Name("sigData_AveragepT"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kRed), LineWidth(2), MarkerColor(kRed), MarkerStyle(20));
  bckData->plotOn(fr_AveragepT, Name("bckData_AveragepT"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kBlue), LineWidth(2), MarkerColor(kBlue), MarkerStyle(24));
  
  RooPlot* fr_R1_nomu = r1_nomu.frame();
  sigData->plotOn(fr_R1_nomu, Name("sigData_R1_nomu"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kRed), LineWidth(2), MarkerColor(kRed), MarkerStyle(20));
  bckData->plotOn(fr_R1_nomu, Name("bckData_R1_nomu"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kBlue), LineWidth(2), MarkerColor(kBlue), MarkerStyle(24));
  
  RooPlot* fr_R3_nomu = r3_nomu.frame();
  sigData->plotOn(fr_R3_nomu, Name("sigData_R3_nomu"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kRed), LineWidth(2), MarkerColor(kRed), MarkerStyle(20));
  bckData->plotOn(fr_R3_nomu, Name("bckData_R3_nomu"), Binning(20), DataError(RooAbsData::SumW2), LineColor(kBlue), LineWidth(2), MarkerColor(kBlue), MarkerStyle(24));
  
  TCanvas* cn_discvar = new TCanvas("cn_discvar","cn_discvar",800,800);
  cn_discvar->cd();
  fr_discvar->Draw();
  TLegend *leg_discvar = new TLegend(0.69,0.74,0.94,0.92);
  leg_style(leg_discvar, 12);
  leg_discvar->SetHeader("Discriminating variable");
  if (name.Contains("Data"))
    leg_discvar->AddEntry("Data","Data - Run 2012 A,B,C,D","P");
  else if (name.Contains("Semi"))
    leg_discvar->AddEntry("Data","MG+PY6 Z2* Semilept. t#bar{t}","P");
  else 
    leg_discvar->AddEntry("Data","MG+PY6 Z2* Dilept. t#bar{t}","P");  
  leg_discvar->AddEntry("total fit",TString::Format("#splitline{Signal:}{N_{sig} = %.0f #pm %.0f}", nsigVal, nsigErr),"F");
  leg_discvar->AddEntry("bck fit",TString::Format("#splitline{Background:}{N_{bck} = %.0f #pm %.0f}", nbckVal, nbckErr),"F");
  // leg_discvar->AddEntry("signal fit","Gaussian component","L");    
  // leg_discvar->AddEntry("bck fit","Exponential component","F");
  leg_discvar->Draw();
  if (name.Contains("Data"))
    cms_style();
  else 
    cms_style(false);
  cn_discvar->SaveAs(rep_name+"D0mass_Fit_"+name+".C");
  cn_discvar->SaveAs(rep_name+"D0mass_Fit_"+name+".eps");
  cn_discvar->SaveAs(rep_name+"D0mass_Fit_"+name+".pdf");
  
  TCanvas* cn_Bmomentum = new TCanvas("cn_Bmomentum","cn_Bmomentum",800,800);
  cn_Bmomentum->cd();
  fr_Bmomentum->Draw();
  TLegend *leg_Bmomentum = new TLegend(0.6,0.8,0.89,0.88);
  leg_style(leg_Bmomentum, 12);
  leg_Bmomentum->AddEntry("sigData_Bmomentum","Weighted signal","P");
  leg_Bmomentum->AddEntry("bckData_Bmomentum","Weighted background","P");
  leg_Bmomentum->Draw();
  channel_tex_r->Draw("same");  
  if (name.Contains("Data"))
    cms_style();
  else 
    cms_style(false);
  cn_Bmomentum->SaveAs(rep_name+"Bmomentum_sPlot_"+name+".C");
  cn_Bmomentum->SaveAs(rep_name+"Bmomentum_sPlot_"+name+".eps");
  cn_Bmomentum->SaveAs(rep_name+"Bmomentum_sPlot_"+name+".pdf");
  
  TCanvas* cn_CSVdisc = new TCanvas("cn_CSVdisc","cn_CSVdisc",800,800);
  cn_CSVdisc->cd();
  fr_CSVdisc->Draw();
  TLegend *leg_CSVdisc = new TLegend(0.2,0.8,0.49,0.88);
  leg_style(leg_CSVdisc, 12);
  leg_CSVdisc->AddEntry("sigData_CSVdisc","Weighted signal","P");
  leg_CSVdisc->AddEntry("bckData_CSVdisc","Weighted background","P");
  leg_CSVdisc->Draw();
  channel_tex_l->Draw("same");  
  if (name.Contains("Data"))
    cms_style();
  else 
    cms_style(false);
  cn_CSVdisc->SaveAs(rep_name+"CSVdisc_sPlot_"+name+".C");
  cn_CSVdisc->SaveAs(rep_name+"CSVdisc_sPlot_"+name+".eps");
  cn_CSVdisc->SaveAs(rep_name+"CSVdisc_sPlot_"+name+".pdf");
  
  TCanvas* cn_R1 = new TCanvas("cn_R1","cn_R1",800,800);
  cn_R1->cd();
  fr_R1->Draw();
  TLegend *leg_R1 = new TLegend(0.6,0.8,0.89,0.88);
  leg_style(leg_R1, 12);
  leg_R1->AddEntry("sigData_R1","Weighted signal","P");
  leg_R1->AddEntry("bckData_R1","Weighted background","P");  
  leg_R1->Draw();
  channel_tex_r->Draw("same");  
  if (name.Contains("Data"))
    cms_style();
  else 
    cms_style(false);
  cn_R1->SaveAs(rep_name+"R1_sPlot_"+name+".C");
  cn_R1->SaveAs(rep_name+"R1_sPlot_"+name+".eps");
  cn_R1->SaveAs(rep_name+"R1_sPlot_"+name+".pdf");

  TCanvas* cn_R3 = new TCanvas("cn_R3","cn_R3",800,800);
  cn_R3->cd();
  fr_R3->Draw();
  TLegend *leg_R3 = new TLegend(0.2,0.8,0.49,0.88);
  leg_style(leg_R3, 12);
  leg_R3->AddEntry("sigData_R3","Weighted signal","P");
  leg_R3->AddEntry("bckData_R3","Weighted background","P");
  leg_R3->Draw();  
  channel_tex_l->Draw("same");  
  if (name.Contains("Data"))
    cms_style();
  else 
    cms_style(false);
  cn_R3->SaveAs(rep_name+"R3_sPlot_"+name+".C");
  cn_R3->SaveAs(rep_name+"R3_sPlot_"+name+".eps");
  cn_R3->SaveAs(rep_name+"R3_sPlot_"+name+".pdf");
  
  TCanvas* cn_Nch = new TCanvas("cn_Nch","cn_Nch",800,800);
  cn_Nch->cd();
  fr_Nch->Draw();
  TLegend *leg_Nch = new TLegend(0.6,0.8,0.89,0.88);
  leg_style(leg_Nch, 12);
  leg_Nch->AddEntry("sigData_Nch","Weighted signal","P");
  leg_Nch->AddEntry("bckData_Nch","Weighted background","P");
  leg_Nch->Draw();
  channel_tex_r->Draw("same");  
  if (name.Contains("Data"))
    cms_style();
  else 
    cms_style(false);
  cn_Nch->SaveAs(rep_name+"Nch_sPlot_"+name+".C");
  cn_Nch->SaveAs(rep_name+"Nch_sPlot_"+name+".eps");
  cn_Nch->SaveAs(rep_name+"Nch_sPlot_"+name+".pdf");
  
  TCanvas* cn_SumpT = new TCanvas("cn_SumpT","cn_SumpT",800,800);
  cn_SumpT->cd();
  fr_SumpT->Draw();
  TLegend *leg_SumpT = new TLegend(0.6,0.8,0.89,0.88);
  leg_style(leg_SumpT, 12);
  leg_SumpT->AddEntry("sigData_SumpT","Weighted signal","P");
  leg_SumpT->AddEntry("bckData_SumpT","Weighted background","P");
  leg_SumpT->Draw();
  channel_tex_r->Draw("same");  
  if (name.Contains("Data"))
    cms_style();
  else 
    cms_style(false);
  cn_SumpT->SaveAs(rep_name+"SumpT_sPlot_"+name+".C");
  cn_SumpT->SaveAs(rep_name+"SumpT_sPlot_"+name+".eps");
  cn_SumpT->SaveAs(rep_name+"SumpT_sPlot_"+name+".pdf");
  
  TCanvas* cn_AveragepT = new TCanvas("cn_AveragepT","cn_AveragepT",800,800);
  cn_AveragepT->cd();
  fr_AveragepT->Draw();
  TLegend *leg_AveragepT = new TLegend(0.6,0.8,0.89,0.88);
  leg_style(leg_AveragepT, 12);
  leg_AveragepT->AddEntry("sigData_AveragepT","Weighted signal","P");
  leg_AveragepT->AddEntry("bckData_AveragepT","Weighted background","P");
  leg_AveragepT->Draw();
  channel_tex_r->Draw("same");  
  if (name.Contains("Data"))
    cms_style();
  else 
    cms_style(false);
  cn_AveragepT->SaveAs(rep_name+"AveragepT_sPlot_"+name+".C");
  cn_AveragepT->SaveAs(rep_name+"AveragepT_sPlot_"+name+".eps");
  cn_AveragepT->SaveAs(rep_name+"AveragepT_sPlot_"+name+".pdf");
  
  TCanvas* cn_R1_nomu = new TCanvas("cn_R1_nomu","cn_R1_nomu",800,800);
  cn_R1_nomu->cd();
  fr_R1_nomu->Draw();
  TLegend *leg_R1_nomu = new TLegend(0.6,0.8,0.89,0.88);
  leg_style(leg_R1_nomu, 12);
  leg_R1_nomu->AddEntry("sigData_R1_nomu","Weighted signal","P");
  leg_R1_nomu->AddEntry("bckData_R1_nomu","Weighted background","P");  
  leg_R1_nomu->Draw();
  channel_tex_r->Draw("same");  
  if (name.Contains("Data"))
    cms_style();
  else 
    cms_style(false);
  cn_R1_nomu->SaveAs(rep_name+"R1_nomu_sPlot_"+name+".C");
  cn_R1_nomu->SaveAs(rep_name+"R1_nomu_sPlot_"+name+".eps");
  cn_R1_nomu->SaveAs(rep_name+"R1_nomu_sPlot_"+name+".pdf");

  TCanvas* cn_R3_nomu = new TCanvas("cn_R3_nomu","cn_R3_nomu",800,800);
  cn_R3_nomu->cd();
  fr_R3_nomu->Draw();
  TLegend *leg_R3_nomu = new TLegend(0.2,0.8,0.49,0.88);
  leg_style(leg_R3_nomu, 12);
  leg_R3_nomu->AddEntry("sigData_R3_nomu","Weighted signal","P");
  leg_R3_nomu->AddEntry("bckData_R3_nomu","Weighted background","P");
  leg_R3_nomu->Draw();  
  channel_tex_l->Draw("same");  
  if (name.Contains("Data"))
    cms_style();
  else 
    cms_style(false);
  cn_R3_nomu->SaveAs(rep_name+"R3_nomu_sPlot_"+name+".C");
  cn_R3_nomu->SaveAs(rep_name+"R3_nomu_sPlot_"+name+".eps");
  cn_R3_nomu->SaveAs(rep_name+"R3_nomu_sPlot_"+name+".pdf");
  
  if (!inBatch) getchar();
  
  TGraphAsymmErrors **gr_all = new TGraphAsymmErrors*[NsPlots];
  int ig = 0;
  
  assert (ig < NsPlots);
  gr_all[ig] = (TGraphAsymmErrors*)cn_CSVdisc->GetPrimitive("sigData_CSVdisc"); 
  gr_all[ig]->GetXaxis()->SetRangeUser(0.,1.);
  if (name.Contains("Data"))
    graphasymerror_mystyle(gr_all[ig], "Data_Sig_CSVdisc", 2, 1, 1, 0, 0, -1111., -1111., 510, 510, 20, 1, 1.2);  
  else 
    graphasymerror_mystyle(gr_all[ig], "MC_Sig_CSVdisc", 2, 862, 1, 862, 1001, -1111., -1111., 510, 510, 1, 862, 1.2);
  ig++;

  assert (ig < NsPlots);
  gr_all[ig] = (TGraphAsymmErrors*)cn_Bmomentum->GetPrimitive("sigData_Bmomentum");
  gr_all[ig]->GetXaxis()->SetRangeUser(0.,400.);
  if (name.Contains("Data"))
    graphasymerror_mystyle(gr_all[ig], "Data_Sig_Bmomentum", 2, 1, 1, 0, 0, -1111., -1111., 510, 510, 20, 1, 1.2);  
  else 
    graphasymerror_mystyle(gr_all[ig], "MC_Sig_Bmomentum", 2, 862, 1, 862, 1001, -1111., -1111., 510, 510, 1, 862, 1.2);
  ig++;
  
  assert (ig < NsPlots);
  gr_all[ig] = (TGraphAsymmErrors*)cn_R1->GetPrimitive("sigData_R1"); 
  gr_all[ig]->GetXaxis()->SetRangeUser(0.,1.);
  if (name.Contains("Data")) 
    graphasymerror_mystyle(gr_all[ig], "Data_Sig_R1", 2, 1, 1, 0, 0, -1111., -1111., 510, 510, 20, 1, 1.2);  
  else 
    graphasymerror_mystyle(gr_all[ig], "MC_Sig_R1", 2, 862, 1, 862, 1001, -1111., -1111., 510, 510, 1, 862, 1.2);
  ig++;
  
  assert (ig < NsPlots);
  gr_all[ig] = (TGraphAsymmErrors*)cn_R3->GetPrimitive("sigData_R3"); 
  gr_all[ig]->GetXaxis()->SetRangeUser(0.,1.);
  if (name.Contains("Data"))
    graphasymerror_mystyle(gr_all[ig], "Data_Sig_R3", 2, 1, 1, 0, 0, -1111., -1111., 510, 510, 20, 1, 1.2);  
  else 
    graphasymerror_mystyle(gr_all[ig], "MC_Sig_R3", 2, 862, 1, 862, 1001, -1111., -1111., 510, 510, 1, 862, 1.2);
  ig++;
  
  assert (ig < NsPlots);
  gr_all[ig] = (TGraphAsymmErrors*)cn_Nch->GetPrimitive("sigData_Nch"); 
  gr_all[ig]->GetXaxis()->SetRangeUser(0.,30.);
  if (name.Contains("Data"))
    graphasymerror_mystyle(gr_all[ig], "Data_Sig_Nch", 2, 1, 1, 0, 0, -1111., -1111., 510, 510, 20, 1, 1.2);  
  else 
    graphasymerror_mystyle(gr_all[ig], "MC_Sig_Nch", 2, 862, 1, 862, 1001, -1111., -1111., 510, 510, 1, 862, 1.2);
  ig++;
  
  assert (ig < NsPlots);
  gr_all[ig] = (TGraphAsymmErrors*)cn_SumpT->GetPrimitive("sigData_SumpT"); 
  gr_all[ig]->GetXaxis()->SetRangeUser(0.,1000.);
  if (name.Contains("Data"))
    graphasymerror_mystyle(gr_all[ig], "Data_Sig_SumpT", 2, 1, 1, 0, 0, -1111., -1111., 510, 510, 20, 1, 1.2);  
  else 
    graphasymerror_mystyle(gr_all[ig], "MC_Sig_SumpT", 2, 862, 1, 862, 1001, -1111., -1111., 510, 510, 1, 862, 1.2);
  ig++;
  
  assert (ig < NsPlots);
  gr_all[ig] = (TGraphAsymmErrors*)cn_AveragepT->GetPrimitive("sigData_AveragepT"); 
  gr_all[ig]->GetXaxis()->SetRangeUser(0.,100.);
  if (name.Contains("Data"))
    graphasymerror_mystyle(gr_all[ig], "Data_Sig_AveragepT", 2, 1, 1, 0, 0, -1111., -1111., 510, 510, 20, 1, 1.2);  
  else 
    graphasymerror_mystyle(gr_all[ig], "MC_Sig_AveragepT", 2, 862, 1, 862, 1001, -1111., -1111., 510, 510, 1, 862, 1.2);
  ig++;
  
  assert (ig < NsPlots);
  gr_all[ig] = (TGraphAsymmErrors*)cn_R1_nomu->GetPrimitive("sigData_R1_nomu"); 
  gr_all[ig]->GetXaxis()->SetRangeUser(0.,1.);
  if (name.Contains("Data")) 
    graphasymerror_mystyle(gr_all[ig], "Data_Sig_R1_nomu", 2, 1, 1, 0, 0, -1111., -1111., 510, 510, 20, 1, 1.2);  
  else 
    graphasymerror_mystyle(gr_all[ig], "MC_Sig_R1_nomu", 2, 862, 1, 862, 1001, -1111., -1111., 510, 510, 1, 862, 1.2);
  ig++;
  
  assert (ig < NsPlots);
  gr_all[ig] = (TGraphAsymmErrors*)cn_R3_nomu->GetPrimitive("sigData_R3_nomu"); 
  gr_all[ig]->GetXaxis()->SetRangeUser(0.,1.);
  if (name.Contains("Data"))
    graphasymerror_mystyle(gr_all[ig], "Data_Sig_R3_nomu", 2, 1, 1, 0, 0, -1111., -1111., 510, 510, 20, 1, 1.2);  
  else 
    graphasymerror_mystyle(gr_all[ig], "MC_Sig_R3_nomu", 2, 862, 1, 862, 1001, -1111., -1111., 510, 510, 1, 862, 1.2);
  ig++;

  if (!name.Contains("Data")) {
    RooPlot* fr_GenId = genId.frame();
    sigData->plotOn(fr_GenId, Name("sigData_GenId"), Binning(22), DataError(RooAbsData::SumW2), LineColor(kRed), LineWidth(2), MarkerColor(kRed), MarkerStyle(20));
    bckData->plotOn(fr_GenId, Name("bckData_GenId"), Binning(22), DataError(RooAbsData::SumW2), LineColor(kBlue), LineWidth(2), MarkerColor(kBlue), MarkerStyle(24));

    TCanvas* cn_GenId = new TCanvas("cn_GenId","cn_GenId",800,800);
    cn_GenId->cd();
    fr_GenId->Draw();
    TLegend *leg_GenId = new TLegend(0.6,0.8,0.89,0.88);
    leg_style(leg_GenId, 12);
    leg_GenId->AddEntry("sigData_GenId","Weighted signal","P");
    leg_GenId->AddEntry("bckData_GenId","Weighted background","P");
    leg_GenId->Draw();  
    channel_tex_r->Draw("same");  
    if (name.Contains("Data"))
      cms_style();
    else 
      cms_style(false);
    cn_GenId->SaveAs(rep_name+"GenId_sPlot_"+name+".C");
    cn_GenId->SaveAs(rep_name+"GenId_sPlot_"+name+".eps");
    cn_GenId->SaveAs(rep_name+"GenId_sPlot_"+name+".pdf");

    if (!inBatch) getchar();
    
    delete cn_GenId;
    delete fr_GenId;
    delete leg_GenId;
  }

  delete cn_discvar;
/*  delete cn_Bmomentum;
  delete cn_CSVdisc;
  delete cn_R1;
  delete cn_R3;
  delete cn_Nch;
  delete cn_SumpT;
  delete cn_AveragepT;
  delete cn_R1_nomu;
  delete cn_R3_nomu;
*/  
  delete leg_discvar;
  delete leg_Bmomentum;
  delete leg_CSVdisc;
  delete leg_R1;
  delete leg_R3;
  delete leg_Nch;
  delete leg_SumpT;
  delete leg_AveragepT;
  delete leg_R1_nomu;
  delete leg_R3_nomu;
  
  delete fr_discvar;
/*  delete fr_Bmomentum;
  delete fr_CSVdisc;
  delete fr_R1;
  delete fr_R3;
  delete fr_Nch;
  delete fr_SumpT;
  delete fr_AveragepT;
  delete fr_R1_nomu;
  delte fr_R3_nomu;
*/  
//  delete bckData;
//  delete sigData;
//  delete sData;
  delete model;
//  delete data;
//  delete tree;
  delete channel_tex_r;
  delete channel_tex_l;
  
  return gr_all;
}

//---------------------------------------------------------------
void doRivetJob_file(bool inBatch, TString date, bool isCSVbased)
//---------------------------------------------------------------
{
  TStyle* my_style = createMyStyle();
  my_style->cd();
  if (inBatch) gROOT->SetBatch(true);

  TString rep_name = date;
  TString fi_data_name = "../test/crab_results/"+date;
  TString fi_sl_name = "../test/crab_results/"+date;
  TString fi_dl_name = "../test/crab_results/"+date;
  if (isCSVbased) {
    rep_name = rep_name + "/csv/sPlot/";
    fi_data_name = fi_data_name + "/D0ForRivet_csv_Data_merged.root";
    fi_sl_name = fi_sl_name + "/D0ForRivet_csv_TTJets_SemiLeptMGDecays.root";
    fi_dl_name = fi_dl_name + "/D0ForRivet_csv_TTJets_FullLeptMGDecays.root";
  }
  else {
    rep_name = rep_name + "/pT/sPlot/";
    fi_data_name = fi_data_name + "/D0ForRivet_pT_Data_merged.root";
    fi_sl_name = fi_sl_name + "/D0ForRivet_pT_TTJets_SemiLeptMGDecays.root";
    fi_dl_name = fi_dl_name + "/D0ForRivet_pT_TTJets_FullLeptMGDecays.root";
  }
    
  TFile* fi_data = TFile::Open(fi_data_name);
  TFile* fi_sl = TFile::Open(fi_sl_name);
  TFile* fi_dl = TFile::Open(fi_dl_name);
  
  const double norm_sl = 2.*25.8031*19769./12031276.;
  const double norm_dl = 2.*107.6722*19769./25339818.;
  const int NsPlots = 9;
  vector<TString> xtitle; xtitle.push_back("CSV discriminant"); xtitle.push_back("p(#kappa^{+}#pi^{-}+#mu^{-}) (GeV/c)"); xtitle.push_back("R1"); xtitle.push_back("R3");
  xtitle.push_back("Tracks multiplicity"); xtitle.push_back("#sum p_{T} (GeV/c)"); xtitle.push_back("#LT p_{T} #GT = #frac{#sum p_{T}}{N_{tracks}} (GeV/c)"); xtitle.push_back("R1 (no #mu)"); xtitle.push_back("R3 (no #mu)");
  vector<TString> grtitle; grtitle.push_back("CSVdisc_Sig_Data2MC"); grtitle.push_back("Bmomentum_Sig_Data2MC"); grtitle.push_back("R1_Sig_Data2MC"); grtitle.push_back("R3_Sig_Data2MC");
  grtitle.push_back("Nch_Sig_Data2MC"); grtitle.push_back("SumpT_Sig_Data2MC"); grtitle.push_back("AveragepT_Sig_Data2MC"); grtitle.push_back("R1_nomu_Sig_Data2MC"); grtitle.push_back("R3_nomu_Sig_Data2MC");
  assert (xtitle.size() == grtitle.size());
  assert (xtitle.size() == NsPlots);
  
  TGraphAsymmErrors **gr_all_data = treatHisto(inBatch, my_style, date, isCSVbased, NsPlots, fi_data, "Data_merged", "Data - Run 2012 A,B,C,D");
  TGraphAsymmErrors **gr_all_sl = treatHisto(inBatch, my_style, date, isCSVbased, NsPlots, fi_sl, "TTJets_SemiLeptMGDecays", "MG+PY6 Z2* Semilept. t#bar{t}");
  TGraphAsymmErrors **gr_all_dl = treatHisto(inBatch, my_style, date, isCSVbased, NsPlots, fi_dl, "TTJets_FullLeptMGDecays", "MG+PY6 Z2* Dilept. t#bar{t}");

  for (int ig = 0; ig < NsPlots; ig++) {
    TGraphAsymmErrors *gr_alldata = (TGraphAsymmErrors*)gr_all_data[ig]->Clone(); 
    TGraphAsymmErrors *gr_allmc = (TGraphAsymmErrors*)gr_all_sl[ig]->Clone(); 
    double x_data = 0.;
    double y_data = 0.;
    double x_sl = 0.;
    double y_sl = 0.;
    double x_dl = 0.;
    double y_dl = 0.;
    double Y_data = 0.;
    double Y_mc = 0.;
    for (int ip = 0; ip < gr_all_data[ig]->GetN(); ip++) {
      gr_all_data[ig]->GetPoint(ip, x_data, y_data);
      Y_data += y_data;
    }
    for (int ip = 0; ip < gr_allmc->GetN(); ip++) {
      gr_all_sl[ig]->GetPoint(ip, x_sl, y_sl);
      gr_all_dl[ig]->GetPoint(ip, x_dl, y_dl);
      assert (fabs(x_sl-x_dl) < 1e-6);
      assert (fabs(gr_all_sl[ig]->GetErrorXlow(ip)-gr_all_dl[ig]->GetErrorXlow(ip)) < 1e-6);
      assert (fabs(gr_all_sl[ig]->GetErrorXhigh(ip)-gr_all_dl[ig]->GetErrorXhigh(ip)) < 1e-6);
      Y_mc += y_sl*norm_sl + y_dl*norm_dl;
    }
    for (int ip = 0; ip < gr_alldata->GetN(); ip++) {
      gr_all_data[ig]->GetPoint(ip, x_data, y_data);
      if (Y_data > 0) {
        gr_alldata->SetPoint(ip, x_data, y_data);
        gr_alldata->SetPointError(ip, gr_all_data[ig]->GetErrorXlow(ip), gr_all_data[ig]->GetErrorXhigh(ip), gr_all_data[ig]->GetErrorYlow(ip), gr_all_data[ig]->GetErrorYhigh(ip));
      }
      else {
        gr_alldata->SetPoint(ip, x_data, 0.);
        gr_alldata->SetPointError(ip, gr_all_data[ig]->GetErrorXlow(ip), gr_all_data[ig]->GetErrorXhigh(ip), 0., 0.);
      }
    }
    for (int ip = 0; ip < gr_allmc->GetN(); ip++) {
      gr_all_sl[ig]->GetPoint(ip, x_sl, y_sl);
      gr_all_dl[ig]->GetPoint(ip, x_dl, y_dl);
      if (Y_mc > 0) {
        gr_allmc->SetPoint(ip, x_sl, (y_sl*norm_sl+y_dl*norm_dl)*Y_data/Y_mc);
        gr_allmc->SetPointError(ip, gr_all_sl[ig]->GetErrorXlow(ip), 
                                      gr_all_sl[ig]->GetErrorXhigh(ip), 
                                     (gr_all_sl[ig]->GetErrorYlow(ip)*norm_sl+gr_all_dl[ig]->GetErrorYlow(ip)*norm_dl)*Y_data/Y_mc, 
                                     (gr_all_sl[ig]->GetErrorYhigh(ip)*norm_sl+gr_all_dl[ig]->GetErrorYhigh(ip)*norm_dl)*Y_data/Y_mc
                                 );
      }
      else {
        gr_allmc->SetPoint(ip, x_sl, 0.);
        gr_allmc->SetPointError(ip, gr_all_sl[ig]->GetErrorXlow(ip), gr_all_sl[ig]->GetErrorXhigh(ip), 0., 0.);
      }
    }
    double x1_data; double x2_data;
    double y1_data; double y2_data;
    gr_all_data[ig]->GetPoint(0, x1_data, y1_data);
    gr_all_data[ig]->GetPoint(1, x2_data, y2_data);
    gr_alldata->GetXaxis()->SetTitle(xtitle[ig]);
    if (xtitle[ig].Contains("GeV/c"))
      gr_alldata->GetYaxis()->SetTitle(TString::Format("Events / (%.2f GeV/c)", x2_data-x1_data));
    else
      gr_alldata->GetYaxis()->SetTitle(TString::Format("Events / %.2f", x2_data-x1_data));
        
    TCanvas* cn_all = new TCanvas("cn_all","cn_all",800,800);
    cn_all->cd();
    gr_alldata->Draw("AP");
    gr_allmc->Draw("2P");
    TLegend *leg_all = new TLegend(0.56,0.84,0.85,0.92); 
    leg_style(leg_all, 12);
    leg_all->AddEntry(gr_alldata,"Data - Run 2012 A,B,C,D","P");
    leg_all->AddEntry(gr_allmc,"MG+PY6 Z2 t#bar{t}","F");
    leg_all->Draw();
    gr_alldata->Draw("P");
    cms_style();
    cn_all->SaveAs(rep_name+grtitle[ig]+".C"); 
    cn_all->SaveAs(rep_name+grtitle[ig]+".eps"); 
    cn_all->SaveAs(rep_name+grtitle[ig]+".pdf"); 
    if (!inBatch) getchar();

    delete cn_all;
    delete leg_all;
    delete gr_alldata;
    delete gr_allmc;
  }
  
  delete gr_all_data;
  delete gr_all_sl;
  delete gr_all_dl;
  
  fi_data->Close();
  fi_sl->Close();
  fi_dl->Close();
  delete fi_data;
  delete fi_sl;
  delete fi_dl;
    
  delete my_style;
}
  
//---------------------------------------------------------------
int doRivetJob_sPlot(bool inBatch = true, TString date = "")
//---------------------------------------------------------------
{  
  if (date.Length() > 0)  {
    gROOT->ProcessLine(".! mkdir "+date);
    gROOT->ProcessLine(".! mkdir "+date+"/csv");
    gROOT->ProcessLine(".! mkdir "+date+"/csv/sPlot");
    gROOT->ProcessLine(".! mkdir "+date+"/pT");
    gROOT->ProcessLine(".! mkdir "+date+"/pT/sPlot");

    doRivetJob_file(inBatch, date, true);
    doRivetJob_file(inBatch, date, false);
    return 0;
  }
  else 
    return 1;
}

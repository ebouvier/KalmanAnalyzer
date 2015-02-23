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
TGraphAsymmErrors **treatHisto(bool inBatch, TStyle* my_style, TString date, TString type, bool isKVF, const int NsPlots, vector<TString> varName, vector<TString> varTitle, vector<double> varMin, vector<double> varMax, vector<TString> varUnit, vector<bool> atRight, TFile* fi, TString name, TString channel)
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
  if (type.Contains("csv", TString::kIgnoreCase)) {
    dir_name = "CSV" + dir_name;
    rep_name = rep_name + "/csv/sPlot/";
  }
  else if (type.Contains("pT", TString::kIgnoreCase)) { 
    dir_name = "pT" + dir_name;
    rep_name = rep_name + "/pT/sPlot/";
  }
  else if (type.Contains("muTag", TString::kIgnoreCase)) { 
    dir_name = "muTagBasedSelection/";
    rep_name = rep_name + "/muTag/sPlot/";
  }
  if (isKVF)
    rep_name = rep_name + "KVF/";
  else
    rep_name = rep_name + "Combi/";
  
  // Import the tree
  TString tree_name = dir_name + "/D0window-b-jets";
  if (isKVF)
    tree_name = dir_name + "/D0KVFwindow-b-jets";
  TTree *tree = (TTree*) fi->Get(tree_name);
  RooRealVar weight("Weight", "Weight", 0., 2.);
  RooRealVar d0mass("D0mass", "D^{0} mass", 1.7, 2., "GeV/c^{2}");  
  RooArgSet variables(weight, d0mass);
  RooRealVar *var[NsPlots];
  for (int iVar = 0; iVar < NsPlots; iVar++) {
    var[iVar] = new RooRealVar(varName[iVar], varTitle[iVar], varMin[iVar], varMax[iVar], varUnit[iVar]);
    variables.add(*var[iVar]);
  }
  RooRealVar genId("GenId", "True jets flavour", 0., 22.);
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
  double nsigVal = nsig.getVal();
  double nsigErr = nsig.getError();
  double nbckVal = nbck.getVal();
  double nbckErr = nbck.getError();

  RooPlot* fr_discvar = d0mass.frame();
  data->plotOn(fr_discvar, Name("Data"), Binning(30));
  model->plotOn(fr_discvar, FillStyle(0), MoveToBack(), Name ("total fit")); // LineColor(9)
  model->plotOn(fr_discvar, Components(bck), LineColor(1), LineWidth(2), FillStyle(1001), FillColor(920), DrawOption("LF"), MoveToBack(), Name("bck fit"));
  // model->plotOn(fr_discvar, Components(bck), LineColor(kBlue), Name("bck fit"));
  // model->plotOn(fr_discvar, Components(RooArgSet(sig)), LineColor(kRed), Name("signal fit"));

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

  if (!inBatch) getchar();

  delete cn_discvar;
  delete leg_discvar;
  delete fr_discvar;
  
  // Apply the sPlot technique
  mean.setConstant();
  sigma.setConstant();
  lambda.setConstant();
  
  RooStats::SPlot* sData = new RooStats::SPlot("sData", "SPlot for m(D0)", *data, model, RooArgList(nsig, nbck));
  ws->import(*data, Rename("dataWithSWeights"));
  data = (RooDataSet*)ws->data("dataWithSWeights");
  RooDataSet *sigData = new RooDataSet(data->GetName(), data->GetTitle(), data, *data->get(), 0, "nsig_sw");
  RooDataSet *bckData = new RooDataSet(data->GetName(), data->GetTitle(), data, *data->get(), 0, "nbck_sw");
  
  RooPlot *fr[NsPlots];
  TCanvas *cn[NsPlots];
  TGraphAsymmErrors **gr_all = new TGraphAsymmErrors*[NsPlots];
  for (int iVar = 0; iVar < NsPlots; iVar++) {
    fr[iVar] = var[iVar]->frame();
    sigData->plotOn(fr[iVar], Name("sigData_"+varName[iVar]), Binning(20), DataError(RooAbsData::SumW2), LineColor(kRed), LineWidth(2), MarkerColor(kRed), MarkerStyle(20));
    bckData->plotOn(fr[iVar], Name("bckData_"+varName[iVar]), Binning(20), DataError(RooAbsData::SumW2), LineColor(kBlue), LineWidth(2), MarkerColor(kBlue), MarkerStyle(24));
    
    cn[iVar] = new TCanvas("cn_"+varName[iVar], "cn_"+varName[iVar], 800, 800);
    cn[iVar]->cd();
    fr[iVar]->Draw();
    double xpos[2] = {0.2, 0.49};
    if (atRight[iVar]) {
      xpos[0] = 0.6;
      xpos[1] = 0.89;
    }
    TLegend *leg = new TLegend(xpos[0], 0.8, xpos[1], 0.88);
    leg_style(leg, 12);
    leg->AddEntry("sigData_"+varName[iVar],"Weighted signal","P");
    leg->AddEntry("bckData_"+varName[iVar],"Weighted background","P");
    leg->Draw();
    if (atRight[iVar]) 
      channel_tex_r->Draw("same");  
    else
      channel_tex_l->Draw("same");  
    if (name.Contains("Data"))
      cms_style();
    else 
      cms_style(false);
    cn[iVar]->SaveAs(rep_name+varName[iVar]+"_sPlot_"+name+".C");
    cn[iVar]->SaveAs(rep_name+varName[iVar]+"_sPlot_"+name+".eps");
    cn[iVar]->SaveAs(rep_name+varName[iVar]+"_sPlot_"+name+".pdf");

    if (!inBatch) getchar();
    delete leg;

    gr_all[iVar] = (TGraphAsymmErrors*)cn[iVar]->GetPrimitive("sigData_"+varName[iVar]); 
    gr_all[iVar]->GetXaxis()->SetRangeUser(varMin[iVar],varMax[iVar]);
    if (name.Contains("Data"))
      graphasymerror_mystyle(gr_all[iVar], "Data_Sig_"+varName[iVar], 2, 1, 1, 0, 0, -1111., -1111., 510, 510, 20, 1, 1.2);  
    else 
      graphasymerror_mystyle(gr_all[iVar], "MC_Sig_"+varName[iVar], 2, 862, 1, 862, 1001, -1111., -1111., 510, 510, 1, 862, 1.2);
  }
  
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

//  delete sData;
  delete model;
//  delete data;
//  delete tree;
  delete channel_tex_r;
  delete channel_tex_l;
  
  return gr_all;
}

//---------------------------------------------------------------
void doRivetJob_file(bool inBatch, TString date, TString type, bool isKVF)
//---------------------------------------------------------------
{
  TStyle* my_style = createMyStyle();
  my_style->cd();
  if (inBatch) gROOT->SetBatch(true);

  TString rep_name = date;
  TString fi_data_name = "../test/crab_results/"+date;
  TString fi_sl_name = "../test/crab_results/"+date;
  TString fi_dl_name = "../test/crab_results/"+date;
  if (type.Contains("csv", TString::kIgnoreCase)) {
    rep_name = rep_name + "/csv/sPlot/";
    fi_data_name = fi_data_name + "/D0ForRivet_csv_Data_merged.root";
    fi_sl_name = fi_sl_name + "/D0ForRivet_csv_TTJets_SemiLeptMGDecays.root";
    fi_dl_name = fi_dl_name + "/D0ForRivet_csv_TTJets_FullLeptMGDecays.root";
  }
  else if (type.Contains("pT", TString::kIgnoreCase))  {
    rep_name = rep_name + "/pT/sPlot/";
    fi_data_name = fi_data_name + "/D0ForRivet_pT_Data_merged.root";
    fi_sl_name = fi_sl_name + "/D0ForRivet_pT_TTJets_SemiLeptMGDecays.root";
    fi_dl_name = fi_dl_name + "/D0ForRivet_pT_TTJets_FullLeptMGDecays.root";
  }
  else if (type.Contains("muTag", TString::kIgnoreCase))  {
    rep_name = rep_name + "/muTag/sPlot/";
    fi_data_name = fi_data_name + "/MuTagForRivet_Data_merged.root";
    fi_sl_name = fi_sl_name + "/MuTagForRivet_TTJets_SemiLeptMGDecays.root";
    fi_dl_name = fi_dl_name + "/MuTagForRivet_TTJets_FullLeptMGDecays.root";
  }
  if (isKVF)
    rep_name = rep_name + "KVF/";
  else
    rep_name = rep_name + "Combi/";
    
  TFile* fi_data = TFile::Open(fi_data_name);
  TFile* fi_sl = TFile::Open(fi_sl_name);
  TFile* fi_dl = TFile::Open(fi_dl_name);
  TFile* fi_out = TFile::Open(rep_name+"/MuTagForRivet_Data_sPlot.root", "recreate");
  
  const double norm_sl = 2.*25.8031*19769./12031276.;
  const double norm_dl = 2.*107.6722*19769./25339818.;
  const int NsPlots = 13;
  vector<TString> xTitle; 
  xTitle.push_back("CSV discriminant"); xTitle.push_back("p(#kappa^{+}#pi^{-}+#mu^{-}) (GeV/c)"); xTitle.push_back("M(#kappa^{+}#pi^{-}+#mu^{-}) (GeV/c^{2})"); 
  xTitle.push_back("p(#mu^{#pm}) (GeV/c)"); xTitle.push_back("R1");  xTitle.push_back("R2");
  xTitle.push_back("R3");  xTitle.push_back("Tracks multiplicity"); xTitle.push_back("#sum p_{T} (GeV/c)"); 
  xTitle.push_back("#LT p_{T} #GT = #frac{#sum p_{T}}{N_{tracks}} (GeV/c)"); xTitle.push_back("R1 (no #mu)"); xTitle.push_back("R2 (no #mu)");
  xTitle.push_back("R3 (no #mu)");
  vector<TString> varName; 
  varName.push_back("CSVdisc"); varName.push_back("Bmomentum");  varName.push_back("Bmass");
   varName.push_back("Mup"); varName.push_back("R1"); varName.push_back("R2"); 
  varName.push_back("R3");  varName.push_back("Nch"); varName.push_back("SumpT"); 
  varName.push_back("AveragepT"); varName.push_back("R1_nomu"); varName.push_back("R2_nomu");
  varName.push_back("R3_nomu");
  vector<TString> varTitle; 
  varTitle.push_back("CSV discriminant"); varTitle.push_back("p(#kappa^{+}#pi^{-}+#mu^{-})"); varTitle.push_back("M(#kappa^{+}#pi^{-}+#mu^{-})"); 
  varTitle.push_back("p(#mu^{#pm})"); varTitle.push_back("R1"); varTitle.push_back("R2");
  varTitle.push_back("R3");  varTitle.push_back("Tracks multiplicity"); varTitle.push_back("#sum p_{T}"); 
  varTitle.push_back("#LT p_{T} #GT"); varTitle.push_back("R1 (no #mu)"); varTitle.push_back("R2 (no #mu)");
  varTitle.push_back("R3 (no #mu)");
  vector<double> varMin; 
  varMin.push_back(0.); varMin.push_back(0.); varMin.push_back(0.); 
  varMin.push_back(0.); varMin.push_back(0.); varMin.push_back(0.);
  varMin.push_back(0.); varMin.push_back(0.); varMin.push_back(0.);
  varMin.push_back(0.); varMin.push_back(0.); varMin.push_back(0.);
  varMin.push_back(0.);
  vector<double> varMax; 
  varMax.push_back(1.); varMax.push_back(400.); varMax.push_back(10.);
  varMax.push_back(300.); varMax.push_back(1.); varMax.push_back(1.); 
  varMax.push_back(1.); varMax.push_back(30.); varMax.push_back(400.); 
  varMax.push_back(80.); varMax.push_back(1.); varMax.push_back(1.);
  varMax.push_back(1.);
  vector<TString> varUnit; 
  varUnit.push_back(""); varUnit.push_back("GeV/c"); varUnit.push_back("GeV/c^{2}"); 
  varUnit.push_back("GeV/c"); varUnit.push_back(""); varUnit.push_back(""); 
  varUnit.push_back(""); varUnit.push_back(""); varUnit.push_back("GeV/c"); 
  varUnit.push_back("GeV/c"); varUnit.push_back(""); varUnit.push_back("");
  varUnit.push_back("");
  vector<bool> atRight; 
  atRight.push_back(false); atRight.push_back(true); atRight.push_back(true); 
  atRight.push_back(true); atRight.push_back(true); atRight.push_back(false); 
  atRight.push_back(false); atRight.push_back(true); atRight.push_back(true); 
  atRight.push_back(true); atRight.push_back(true); atRight.push_back(true);
  atRight.push_back(false);
  assert (xTitle.size() == NsPlots);
  assert (varName.size() == NsPlots);
  assert (varTitle.size() == NsPlots);
  assert (varMin.size() == NsPlots);
  assert (varMax.size() == NsPlots);
  assert (varUnit.size() == NsPlots);
  assert (atRight.size() == NsPlots);

  TH1F **graph_data_sig = new TH1F*[NsPlots]; 
  
  TGraphAsymmErrors **gr_all_data = treatHisto(inBatch, my_style, date, type, isKVF, NsPlots, varName, varTitle, varMin, varMax, varUnit, atRight, fi_data, "Data_merged", "Data - Run 2012 A,B,C,D");
  TGraphAsymmErrors **gr_all_sl = treatHisto(inBatch, my_style, date, type, isKVF, NsPlots, varName, varTitle, varMin, varMax, varUnit, atRight, fi_sl, "TTJets_SemiLeptMGDecays", "MG+PY6 Z2* Semilept. t#bar{t}");
  TGraphAsymmErrors **gr_all_dl = treatHisto(inBatch, my_style, date, type, isKVF, NsPlots, varName, varTitle, varMin, varMax, varUnit, atRight, fi_dl, "TTJets_FullLeptMGDecays", "MG+PY6 Z2* Dilept. t#bar{t}");

  for (int ig = 0; ig < NsPlots; ig++) {
    graph_data_sig[ig] = new TH1F(varName[ig], varTitle[ig], gr_all_data[ig]->GetN(), varMin[ig], varMax[ig]);
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
        graph_data_sig[ig]->SetBinContent(1+ip, y_data);
        graph_data_sig[ig]->SetBinError(1+ip, 0.5*(gr_all_data[ig]->GetErrorYlow(ip)+gr_all_data[ig]->GetErrorYhigh(ip)));
        gr_alldata->SetPoint(ip, x_data, y_data);
        gr_alldata->SetPointError(ip, gr_all_data[ig]->GetErrorXlow(ip), gr_all_data[ig]->GetErrorXhigh(ip), gr_all_data[ig]->GetErrorYlow(ip), gr_all_data[ig]->GetErrorYhigh(ip));
      }
      else {
        graph_data_sig[ig]->SetBinContent(1+ip, 0.);
        graph_data_sig[ig]->SetBinError(1+ip, 0.);
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
    graph_data_sig[ig]->GetXaxis()->SetTitle(xTitle[ig]);
    gr_alldata->GetXaxis()->SetTitle(xTitle[ig]);
    if (xTitle[ig].Contains("GeV/c")) {
      graph_data_sig[ig]->GetYaxis()->SetTitle(TString::Format("Events / (%.2f GeV/c)", x2_data-x1_data));
      gr_alldata->GetYaxis()->SetTitle(TString::Format("Events / (%.2f GeV/c)", x2_data-x1_data));
    }
    else {
      graph_data_sig[ig]->GetYaxis()->SetTitle(TString::Format("Events / %.2f", x2_data-x1_data));
      gr_alldata->GetYaxis()->SetTitle(TString::Format("Events / %.2f", x2_data-x1_data));
    }
        
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
    cn_all->SaveAs(rep_name+varName[ig]+"_Sig_Data2MC.C"); 
    cn_all->SaveAs(rep_name+varName[ig]+"_Sig_Data2MC.eps"); 
    cn_all->SaveAs(rep_name+varName[ig]+"_Sig_Data2MC.pdf"); 
    if (!inBatch) getchar();

    fi_out->cd();
    graph_data_sig[ig]->Write();

    delete cn_all;
    delete leg_all;
    delete gr_alldata;
    delete gr_allmc;
  }

  fi_out->Close();
  delete fi_out;
  delete graph_data_sig;
  
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
int doRivetJob_sPlot(TString date = "", TString type = "muTag", bool inBatch = true)
//---------------------------------------------------------------
{  
  if (date.Length() > 0)  {
    gROOT->ProcessLine(".! mkdir "+date);

    if (type.Contains("D0", TString::kIgnoreCase)) {    
      gROOT->ProcessLine(".! mkdir "+date+"/csv");
      gROOT->ProcessLine(".! mkdir "+date+"/csv/sPlot");
      gROOT->ProcessLine(".! mkdir "+date+"/csv/sPlot/Combi");
      gROOT->ProcessLine(".! mkdir "+date+"/csv/sPlot/KVF");
      gROOT->ProcessLine(".! mkdir "+date+"/pT");
      gROOT->ProcessLine(".! mkdir "+date+"/pT/sPlot");
      gROOT->ProcessLine(".! mkdir "+date+"/pT/sPlot/Combi");
      gROOT->ProcessLine(".! mkdir "+date+"/pT/sPlot/KVF");

      doRivetJob_file(inBatch, date, "csv", true);
      doRivetJob_file(inBatch, date, "csv", false);
      doRivetJob_file(inBatch, date, "pT", true);
      doRivetJob_file(inBatch, date, "pT", false);
    }
    else if (type.Contains("MuTag", TString::kIgnoreCase)) {
      gROOT->ProcessLine(".! mkdir "+date+"/muTag");
      gROOT->ProcessLine(".! mkdir "+date+"/muTag/sPlot");
      gROOT->ProcessLine(".! mkdir "+date+"/muTag/sPlot/Combi");
      gROOT->ProcessLine(".! mkdir "+date+"/muTag/sPlot/KVF");

      doRivetJob_file(inBatch, date, "muTag", true);
      doRivetJob_file(inBatch, date, "muTag", false);
    }
    return 0;
  }
  else 
    return 1;
}

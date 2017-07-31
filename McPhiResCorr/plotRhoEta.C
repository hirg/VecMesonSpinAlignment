#include <string>
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "../Utility/functions.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"

void plotRhoEta()
{
  string Order[6] = {"01","02","05","1","2","5"};
  int style[6] = {20,24,21,25,29,30};
  int color[6] = {kGray+2,kAzure-2,kGray+2,kAzure-2,kGray+2,kAzure-2};
  TFile *File_InPut[6];
  TGraphAsymmErrors *g_Kaon[6];
  for(int i_eta = 0; i_eta < 6; ++i_eta)
  {
    string InPutHist = Form("/Users/xusun/Data/SpinAlignment/AuAu200GeV/MonteCarlo/McRho_%s.root",Order[i_eta].c_str());
    File_InPut[i_eta] = TFile::Open(InPutHist.c_str());
    g_Kaon[i_eta] = (TGraphAsymmErrors*)File_InPut[i_eta]->Get("g_Kaon");
  }
  TH1F *h_play = (TH1F*)File_InPut[0]->Get("h_rho");

  TCanvas *c_eta = new TCanvas("c_eta","c_eta",10,10,800,800);
  c_eta->cd()->SetLeftMargin(0.15);
  c_eta->cd()->SetBottomMargin(0.15);
  c_eta->cd()->SetTicks(1,1);
  c_eta->cd()->SetGrid(0,0);
  h_play->Draw("pE");
  for(int i_eta = 0; i_eta < 6; ++i_eta)
  {
    g_Kaon[i_eta]->SetMarkerStyle(style[i_eta]);
    g_Kaon[i_eta]->SetMarkerColor(color[i_eta]);
    g_Kaon[i_eta]->SetMarkerSize(1.4);
    g_Kaon[i_eta]->Draw("pE same");
  }
  PlotLine(0.0,10.0,1.0/3.0,1.0/3.0,1,2,2);

  TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(g_Kaon[0],"0.1*p_{T}","p");
  leg->AddEntry(g_Kaon[1],"0.2*p_{T}","p");
  leg->AddEntry(g_Kaon[2],"0.5*p_{T}","p");
  leg->AddEntry(g_Kaon[3],"p_{T}","p");
  leg->AddEntry(g_Kaon[4],"2.0*p_{T}","p");
  leg->AddEntry(g_Kaon[5],"5.0*p_{T}","p");
  leg->Draw("same");

  c_eta->SaveAs("../figures/c_RhoEta.eps");
}

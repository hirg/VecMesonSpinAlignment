#include <string>
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TProfile.h"
#include "../Utility/functions.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"

void plotMcLambdaEta()
{
  string InPutHist = "/Users/xusun/Data/SpinAlignment/AuAu200GeV/McLambdaEta.root";
  TFile *File_InPut = TFile::Open(InPutHist.c_str());

  TProfile *p_cosRP = (TProfile*)File_InPut->Get("p_cosRP");

  TProfile *p_CosEtaKaon[20], *p_CosEtaPhi[20];
  TGraphAsymmErrors *g_Kaon = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_Phi = new TGraphAsymmErrors();
  for(int i_eta = 0; i_eta < 20; ++i_eta)
  {
    string HistName;
    HistName = Form("p_CosEtaKaon_%d",i_eta);
    p_CosEtaKaon[i_eta] = (TProfile*)File_InPut->Get(HistName.c_str());
    g_Kaon->SetPoint(i_eta,vmsa::McEtaBin[i_eta],p_CosEtaKaon[i_eta]->GetBinContent(1));
    g_Kaon->SetPointError(i_eta,0.0,0.0,p_CosEtaKaon[i_eta]->GetBinError(1),p_CosEtaKaon[i_eta]->GetBinError(1));

    HistName = Form("p_CosEtaPhi_%d",i_eta);
    p_CosEtaPhi[i_eta] = (TProfile*)File_InPut->Get(HistName.c_str());
    g_Phi->SetPoint(i_eta,vmsa::McEtaBin[i_eta],p_CosEtaPhi[i_eta]->GetBinContent(1));
    g_Phi->SetPointError(i_eta,0.0,0.0,p_CosEtaPhi[i_eta]->GetBinError(1),p_CosEtaPhi[i_eta]->GetBinError(1));
  }
  g_Kaon->SetPoint(20,8.0,p_cosRP->GetBinContent(1));
  g_Kaon->SetPointError(20,0.0,0.0,p_cosRP->GetBinError(1),p_cosRP->GetBinError(1));
  g_Phi->SetPoint(20,8.0,p_cosRP->GetBinContent(1));
  g_Phi->SetPointError(20,0.0,0.0,p_cosRP->GetBinError(1),p_cosRP->GetBinError(1));

  TCanvas *c_rhoEta = new TCanvas("c_rhoEta","c_rhoEta",10,10,800,800);
  c_rhoEta->cd()->SetLeftMargin(0.15);
  c_rhoEta->cd()->SetBottomMargin(0.15);
  c_rhoEta->cd()->SetTicks(1,1);
  c_rhoEta->cd()->SetGrid(0,0);

  TH1F *h_rho = new TH1F("h_rho","h_rho",100,0.0,10.0);
  for(int i_bin = 1; i_bin < 101; ++i_bin)
  {
    h_rho->SetBinContent(i_bin,-10.0);
    h_rho->SetBinError(i_bin,1.0);
  }
  h_rho->SetTitle("");
  h_rho->SetStats(0);
  h_rho->GetXaxis()->SetTitle("|#eta| <");
  h_rho->GetXaxis()->CenterTitle();
  h_rho->GetXaxis()->SetRangeUser(0.0,10.0);
  h_rho->GetYaxis()->SetTitle("<cos(#theta*)>");
  h_rho->GetYaxis()->SetTitleSize(0.04);
  h_rho->GetXaxis()->SetNdivisions(505);

  h_rho->GetYaxis()->CenterTitle();
  h_rho->GetYaxis()->SetRangeUser(-0.1,0.1);
  h_rho->GetYaxis()->SetNdivisions(505);
  h_rho->Draw("pE");
  PlotLine(0.0,10.0,0.0,0.0,1,2,2);

  g_Kaon->SetMarkerStyle(24);
  g_Kaon->SetMarkerColor(2);
  g_Kaon->SetMarkerSize(1.4);
  g_Kaon->Draw("pE same");

  g_Phi->SetMarkerStyle(24);
  g_Phi->SetMarkerColor(4);
  g_Phi->SetMarkerSize(1.4);
  g_Phi->Draw("pE same");

  TLegend *leg = new TLegend(0.4,0.6,0.8,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(g_Kaon,"cut on daughters and #Lambda","p");
  leg->AddEntry(g_Phi,"cut on #Lambda only","p");
  leg->Draw("same");
  c_rhoEta->SaveAs("../figures/c_LambdaRhoEta.eps");
}

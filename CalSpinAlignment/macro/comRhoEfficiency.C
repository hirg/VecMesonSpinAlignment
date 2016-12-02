#include "draw.h"
#include "TString.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TCanvas.h"

void Com_rho_Efficiency()
{
  TString InPut_200 = "/Users/xusun/STAR/SpinAlignment/OutPut/AuAu200GeV/Phi/rho00_pT.root";
  TFile *File_200 = TFile::Open(InPut_200.Data());
  TGraphAsymmErrors *g_200 = (TGraphAsymmErrors*)File_200->Get("Rho00_Phi_Centrality_0_EtaGap_0_BW");
  TH1F *h_play = (TH1F*)File_200->Get("h_play");

  TString InPut_200_Eff = "/Users/xusun/STAR/SpinAlignment/OutPut/AuAu200GeV/Phi/rho00_pT_Eff.root";
  TFile *File_200_Eff = TFile::Open(InPut_200_Eff.Data());
  TGraphAsymmErrors *g_200_Eff = (TGraphAsymmErrors*)File_200_Eff->Get("Rho00_Phi_Centrality_0_EtaGap_0_BW");

  TString InPut_200_Eff_EP = "/Users/xusun/STAR/SpinAlignment/OutPut/AuAu200GeV/Phi/rho00_pT_Eff_EP.root";
  TFile *File_200_Eff_EP = TFile::Open(InPut_200_Eff_EP.Data());
  TGraphAsymmErrors *g_200_Eff_EP = (TGraphAsymmErrors*)File_200_Eff_EP->Get("Rho00_Phi_Centrality_0_EtaGap_0_BW");

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->cd()->SetLeftMargin(0.15);
  c_play->cd()->SetBottomMargin(0.15);
  c_play->cd()->SetTicks(1,1);
  c_play->cd()->SetGrid(0,0);
  h_play->GetYaxis()->SetRangeUser(0.3,0.45);
  h_play->Draw("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_200,20,kGray+3,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.44,0.0,0.0,0.0,0.0,24,kGray+3,1.3);
  plotTopLegend("200 GeV w/o correction",0.6,0.4385,0.03,1,0.0,42,0);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_200_Eff,24,2,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.43,0.0,0.0,0.0,0.0,24,2,1.3);
  plotTopLegend("200 GeV with TPC efficiency correction (RP)",0.6,0.4285,0.03,1,0.0,42,0);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_200_Eff_EP,24,4,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.42,0.0,0.0,0.0,0.0,24,4,1.3);
  plotTopLegend("200 GeV with TPC efficiency correction (EP)",0.6,0.4185,0.03,1,0.0,42,0);
  c_play->SaveAs("./figures/Com_rho_Efficiency.eps");
}

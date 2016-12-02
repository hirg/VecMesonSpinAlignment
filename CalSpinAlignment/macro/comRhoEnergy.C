#include "draw.h"
#include "TString.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TCanvas.h"

void Com_rho_Energy()
{
  TString InPut_200 = "/Users/xusun/STAR/SpinAlignment/OutPut/AuAu200GeV/Phi/rho00_pT.root";
  TFile *File_200 = TFile::Open(InPut_200.Data());
  TGraphAsymmErrors *g_200 = (TGraphAsymmErrors*)File_200->Get("Rho00_Phi_Centrality_0_EtaGap_0_BW");
  TH1F *h_play = (TH1F*)File_200->Get("h_play");

  TString InPut_39 = "/Users/xusun/STAR/SpinAlignment/OutPut/AuAu39GeV/Phi/rho00_pT.root";
  TFile *File_39 = TFile::Open(InPut_39.Data());
  TGraphAsymmErrors *g_39 = (TGraphAsymmErrors*)File_39->Get("Rho00_Phi_Centrality_0_EtaGap_0_BW");

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->cd()->SetLeftMargin(0.15);
  c_play->cd()->SetBottomMargin(0.15);
  c_play->cd()->SetTicks(1,1);
  c_play->cd()->SetGrid(0,0);
  h_play->GetYaxis()->SetRangeUser(0.3,0.4);
  h_play->Draw("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_200,24,kGray+3,1.1);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_39,24,2,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.38,0.0,0.0,0.0,0.0,24,kGray+3,1.3);
  plotTopLegend("200 GeV",0.6,0.3785,0.03,1,0.0,42,0);
  Draw_TGAE_Point_new_Symbol(0.5,0.37,0.0,0.0,0.0,0.0,24,2,1.3);
  plotTopLegend("39 GeV",0.6,0.3685,0.03,1,0.0,42,0);
  
  plotTopLegend("w/o corrctions",2.6,0.3685,0.03,1,0.0,42,0);
  c_play->SaveAs("./figures/Com_rho_Energy.eps");
}

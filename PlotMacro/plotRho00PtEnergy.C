#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "../Utility/draw.h"
#include "TCanvas.h"
#include "TH1F.h"

using namespace std;

void plotRho00PtEnergy()
{
  TGraphAsymmErrors *g_rho[5];
  string input19 = "/Users/xusun/Data/SpinAlignment/AuAu19GeV/RawRhoPt.root.oldCut";
  TFile *File_19 = TFile::Open(input19.c_str());
  g_rho[0] = (TGraphAsymmErrors*)File_19->Get("Rho00_Phi_Centrality_0_EtaGap_0_Count");

  string input27 = "/Users/xusun/Data/SpinAlignment/AuAu27GeV/RawRhoPt.root.oldCut";
  TFile *File_27 = TFile::Open(input27.c_str());
  g_rho[1] = (TGraphAsymmErrors*)File_27->Get("Rho00_Phi_Centrality_0_EtaGap_0_Count");

  string input39 = "/Users/xusun/Data/SpinAlignment/AuAu39GeV/RawRhoPt.root.oldCut";
  TFile *File_39 = TFile::Open(input39.c_str());
  g_rho[2] = (TGraphAsymmErrors*)File_39->Get("Rho00_Phi_Centrality_0_EtaGap_0_Count");

  string input62 = "/Users/xusun/Data/SpinAlignment/AuAu62GeV/RawRhoPt.root.oldCut";
  TFile *File_62 = TFile::Open(input62.c_str());
  g_rho[3] = (TGraphAsymmErrors*)File_62->Get("Rho00_Phi_Centrality_0_EtaGap_0_Count");

  string input200 = "/Users/xusun/Data/SpinAlignment/AuAu200GeV/RawRhoPt.root.oldCut";
  TFile *File_200 = TFile::Open(input200.c_str());
  g_rho[4] = (TGraphAsymmErrors*)File_200->Get("Rho00_Phi_Centrality_0_EtaGap_0_Count");

  int const Color[5] = {kGray+2,kAzure,kRed,kCyan,kMagenta};
  int const Style[5] = {20,22,23,24,29};
  double const delta[5] = {-0.1,-0.05,0.0,0.05,0.1};
  float const pos_x[5] = {0.8,2.0,0.8,2.0,0.8};
  float const pos_y[5] = {0.47,0.47,0.45,0.45,0.43};
  string const BeamEnergy[5] = {"19GeV","27GeV","39GeV","62GeV","200GeV"};
  TH1F *h_frame = new TH1F("h_frame","h_frame",100,0,10);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,-10.0);
  }
  TCanvas *c_rho00 = new TCanvas("c_rho00","c_rho00",10,10,800,800);
  c_rho00->cd();
  c_rho00->cd()->SetLeftMargin(0.15);
  c_rho00->cd()->SetBottomMargin(0.15);
  c_rho00->cd()->SetTicks(1,1);
  c_rho00->cd()->SetGrid(0,0);
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,5.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.03);
  h_frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_frame->GetXaxis()->SetTitleSize(0.05);
  h_frame->GetXaxis()->SetTitleOffset(1.2);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.2,0.5);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00}^{raw}");
  h_frame->GetYaxis()->SetTitleSize(0.05);
  h_frame->GetYaxis()->SetLabelSize(0.03);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);
  for(int i_energy = 0; i_energy < 5; ++i_energy)
  {
    for(int i_point = 0; i_point < g_rho[4]->GetN(); ++i_point)
    {
      double pt, rho00;
      g_rho[i_energy]->GetPoint(i_point,pt,rho00);
      g_rho[i_energy]->SetPoint(i_point,pt+delta[i_energy],rho00);
    }
    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho[i_energy],Style[i_energy],Color[i_energy],1.1);
    Draw_TGAE_Point_new_Symbol(pos_x[i_energy],pos_y[i_energy],0.0,0.0,0.0,0.0,Style[i_energy],Color[i_energy],1.3);
    plotTopLegend((char*)BeamEnergy[i_energy].c_str(),pos_x[i_energy]+0.1,pos_y[i_energy]-0.003,0.03,1,0.0,42,0);
  }
  c_rho00->SaveAs("../figures/c_rho00PtEnergy.eps");
}

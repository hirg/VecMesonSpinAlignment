#include <string>
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "../Utility/draw.h"

const float BeamEnergy[7]    = {200.0,62.4,39.0,27.0,19.6,11.5,7.7};
const float rho00[7]         = {3.435e-01,3.472e-01,3.496e-01,3.462e-01,3.456e-01,3.548e-01,3.548e-01};
const float ErrStat[7]       = {1.168e-03,3.570e-03,2.596e-03,4.007e-03,6.035e-03,1.202e-02,2.716e-02};
const float ErrStat_BESII[7] = {0.000e-03,0.000e-03,0.000e-03,0.000e-03,6.035e-03/3.0,1.202e-02/4.5,2.716e-02/5.0};
const float ErrSys_Xu[7]     = {1.083e-02,1.145e-02,1.093e-02,1.665e-02,1.890e-02,0.0,0.0};
const float ErrSys_Chigh[7]  = {1.000e-03,0.000e-01,8.000e-03,3.800e-03,2.940e-02,2.020e-02,0.0};
const float ErrSys_Clow[7]   = {3.000e-03,0.000e-01,2.000e-03,1.020e-02,8.600e-03,4.480e-02,0.0};

void plotBESII()
{
  TGraphAsymmErrors *g_rho_stat = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_rho_sys  = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_rho_sysII  = new TGraphAsymmErrors();
  for(int i_point = 0; i_point < 7; ++i_point)
  {
    g_rho_stat->SetPoint(i_point,BeamEnergy[i_point],rho00[i_point]);
    g_rho_stat->SetPointError(i_point,0.0,0.0,ErrStat[i_point],ErrStat[i_point]);

    g_rho_sys->SetPoint(i_point,BeamEnergy[i_point],rho00[i_point]);
    g_rho_sys->SetPointError(i_point,0.0,0.0,ErrSys_Xu[i_point],ErrSys_Xu[i_point]);

    g_rho_sysII->SetPoint(i_point,BeamEnergy[i_point],rho00[i_point]);
    g_rho_sysII->SetPointError(i_point,0.0,0.0,ErrSys_Clow[i_point],ErrSys_Chigh[i_point]);
  }

  TGraphAsymmErrors *g_rho_stat_est = new TGraphAsymmErrors();
  for(int i_point = 0; i_point < 2; ++i_point)
  {
    g_rho_stat_est->SetPoint(i_point,BeamEnergy[i_point+5],rho00[i_point+5]);
    g_rho_stat_est->SetPointError(i_point,0.0,0.0,ErrStat[i_point+5],ErrStat[i_point+5]);
  }

  TGraphAsymmErrors *g_rho_statII = new TGraphAsymmErrors();
  for(int i_point = 0; i_point < 3; ++i_point)
  {
    g_rho_statII->SetPoint(i_point,BeamEnergy[i_point+4]*0.9,rho00[i_point+4]);
    g_rho_statII->SetPointError(i_point,0.0,0.0,ErrStat_BESII[i_point+4],ErrStat_BESII[i_point+4]);
  }

  TH1F *h_frame = new TH1F("h_frame","h_frame",1000,0,1000);
  for(int i_bin = 0; i_bin < 1000; ++i_bin)
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
  c_rho00->cd()->SetLogx();
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(5.00,500.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("#sqrt{s_{NN}}  (GeV)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.301,0.43);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.2);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(5.0,500.0,1.0/3.0,1.0/3.0,1,3,2);
  plotTopLegend((char*)"#rho_{00} = 1/3",100,0.328,0.04,1,0.0,42,0);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_stat,29,kRed,1.8);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_stat_est,30,kRed,1.8);
  for(int i_point = 0; i_point < 5; ++i_point)
  {
    double energy = BeamEnergy[i_point];
    double rho = rho00[i_point];
    double err = ErrSys_Xu[i_point];

    PlotLine(energy*0.95,energy*1.05,rho+err,rho+err,kGray+2,2,1);
    PlotLine(energy*0.95,energy*0.95,rho+err-0.001,rho+err,kGray+2,2,1);
    PlotLine(energy*1.05,energy*1.05,rho+err-0.001,rho+err,kGray+2,2,1);
    PlotLine(energy*0.95,energy*1.05,rho-err,rho-err,kGray+2,2,1);
    PlotLine(energy*0.95,energy*0.95,rho-err+0.001,rho-err,kGray+2,2,1);
    PlotLine(energy*1.05,energy*1.05,rho-err+0.001,rho-err,kGray+2,2,1);
  }
  g_rho_statII->SetLineColor(kAzure+2);
  g_rho_statII->SetLineStyle(1);
  g_rho_statII->SetLineWidth(6);
  g_rho_statII->Draw("pE same");

  plotTopLegend((char*)"Au+Au (20-60\%)",7,0.42,0.04,1,0.0,42,0);
  plotTopLegend((char*)"|#eta| < 1",200,0.42,0.04,1,0.0,42,0);
  Draw_TGAE_Point_new_Symbol(7,0.41,0.0,0.0,0.0,0.0,29,kRed,1.8);
  plotTopLegend((char*)"#phi-meson (0.4 < p_{T}< 3.0 GeV/c)",8,0.409,0.04,1,0.0,42,0);
  Draw_TGAE_Point_new_Symbol(7,0.40,0.0,0.0,0.0,0.0,30,kRed,1.8);
  plotTopLegend((char*)"estimated #rho_{00}",8,0.399,0.04,1,0.0,42,0);
  PlotLine(6.5,7.5,0.391,0.391,kAzure+2,6,1);
  plotTopLegend((char*)"BES-II statistic",8,0.389,0.04,1,0.0,42,0);
  c_rho00->SaveAs("../figures/c_rho00_Xu.eps");



  TCanvas *c_rho00_sys = new TCanvas("c_rho00_sys","c_rho00_sys",810,10,800,800);
  c_rho00_sys->cd();
  c_rho00_sys->cd()->SetLeftMargin(0.15);
  c_rho00_sys->cd()->SetBottomMargin(0.15);
  c_rho00_sys->cd()->SetTicks(1,1);
  c_rho00_sys->cd()->SetGrid(0,0);
  c_rho00_sys->cd()->SetLogx();
  h_frame->DrawCopy("pE");
  PlotLine(5.0,500.0,1.0/3.0,1.0/3.0,1,3,2);
  plotTopLegend((char*)"#rho_{00} = 1/3",100,0.328,0.04,1,0.0,42,0);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_stat,29,kRed,1.8);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_stat_est,30,kRed,1.8);
  for(int i_point = 0; i_point < 5; ++i_point)
  {
    if(i_point == 1) continue;
    double energy = BeamEnergy[i_point];
    double rho = rho00[i_point];
    double errl = ErrSys_Clow[i_point];
    double errh = ErrSys_Chigh[i_point];

    PlotLine(energy*0.95,energy*1.05,rho+errh,rho+errh,kGray+2,2,1);
    PlotLine(energy*0.95,energy*0.95,rho+errh-0.001,rho+errh,kGray+2,2,1);
    PlotLine(energy*1.05,energy*1.05,rho+errh-0.001,rho+errh,kGray+2,2,1);
    PlotLine(energy*0.95,energy*1.05,rho-errl,rho-errl,kGray+2,2,1);
    PlotLine(energy*0.95,energy*0.95,rho-errl+0.001,rho-errl,kGray+2,2,1);
    PlotLine(energy*1.05,energy*1.05,rho-errl+0.001,rho-errl,kGray+2,2,1);
  }
  g_rho_statII->SetLineColor(kAzure+2);
  g_rho_statII->SetLineStyle(1);
  g_rho_statII->SetLineWidth(6);
  g_rho_statII->Draw("pE same");

  plotTopLegend((char*)"Au+Au (20-60\%)",7,0.42,0.04,1,0.0,42,0);
  plotTopLegend((char*)"|#eta| < 1",200,0.42,0.04,1,0.0,42,0);
  Draw_TGAE_Point_new_Symbol(7,0.41,0.0,0.0,0.0,0.0,29,kRed,1.8);
  plotTopLegend((char*)"#phi-meson (0.4 < p_{T}< 3.0 GeV/c)",8,0.409,0.04,1,0.0,42,0);
  Draw_TGAE_Point_new_Symbol(7,0.40,0.0,0.0,0.0,0.0,30,kRed,1.8);
  plotTopLegend((char*)"estimated #rho_{00}",8,0.399,0.04,1,0.0,42,0);
  PlotLine(6.5,7.5,0.391,0.391,kAzure+2,6,1);
  plotTopLegend((char*)"BES-II statistic",8,0.389,0.04,1,0.0,42,0);
  c_rho00_sys->SaveAs("../figures/c_rho00_Chensheng.eps");

  TFile *File_OutPut = new TFile("/Users/xusun/Data/SpinAlignment/BESII_rho00.root","RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  g_rho_stat->SetName("g_rho_stat");
  g_rho_stat->Write();
  g_rho_sys->SetName("g_rho_sys_Xu");
  g_rho_sys->Write();
  g_rho_sysII->SetName("g_rho_sys_Chensheng");
  g_rho_sysII->Write();
  g_rho_statII->SetName("g_rho_stat_BESII");
  g_rho_statII->Write();
  File_OutPut->Close();
}

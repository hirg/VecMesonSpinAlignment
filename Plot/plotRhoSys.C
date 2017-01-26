#include "string.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "../Utility/draw.h"
#include "TCanvas.h"

using namespace std;

void plotRhoSys()
{
  string Energy[5] = {"19GeV","27GeV","39GeV","62GeV","200GeV"};
  int plot_style[5] = {20,29,34,22,24};
  int plot_color[5] = {kGray+3,kRed-4,kMagenta,kAzure+4,4};
  int pt_bin[5]  = {5,5,5,5,7};
  TFile *File_InPut[5];
  TGraphAsymmErrors *g_stat[5];
  TGraphAsymmErrors *g_sys[5];
  TH1F *h_frame;
  for(int i_energy = 0; i_energy < 5; ++i_energy)
  {
    string inputfile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/rho00/Rho_SysErrors.root",Energy[i_energy].c_str());
    File_InPut[i_energy] = TFile::Open(inputfile.c_str());
    if(i_energy == 0) h_frame = (TH1F*)File_InPut[i_energy]->Get("h_frame");

    string stat = Form("g_rho00_%s_Phi_StatError",Energy[i_energy].c_str());
    TGraphAsymmErrors *g_stat_temp = (TGraphAsymmErrors*)File_InPut[i_energy]->Get(stat.c_str());
    g_stat[i_energy] = new TGraphAsymmErrors();

    string sys = Form("g_rho00_%s_Phi_SysError",Energy[i_energy].c_str());
    TGraphAsymmErrors *g_sys_temp = (TGraphAsymmErrors*)File_InPut[i_energy]->Get(sys.c_str());
    g_sys[i_energy] = new TGraphAsymmErrors();

    for(int i_point = 0; i_point < pt_bin[i_energy]; ++i_point)
    {
      double pt_stat, rho_stat;
      g_stat_temp->GetPoint(i_point,pt_stat,rho_stat);
      double err_stat = g_stat_temp->GetErrorYhigh(i_point);
      g_stat[i_energy]->SetPoint(i_point,pt_stat-0.25+i_energy*0.05,rho_stat);
      g_stat[i_energy]->SetPointError(i_point,0.0,0.0,err_stat,err_stat);

      double pt_sys, rho_sys;
      g_sys_temp->GetPoint(i_point,pt_sys,rho_sys);
      double err_sys = g_sys_temp->GetErrorYhigh(i_point);
      g_sys[i_energy]->SetPoint(i_point,pt_sys-0.25+i_energy*0.05,rho_sys);
      g_sys[i_energy]->SetPointError(i_point,0.0,0.0,err_sys,err_sys);
    }
  }

  TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,800,800);
  c_rho->SetLeftMargin(0.15);
  c_rho->SetBottomMargin(0.15);
  c_rho->SetGrid(0,0);
  c_rho->SetTicks(1,1);
  c_rho->cd();
  h_frame->GetXaxis()->SetRangeUser(0.0,5.0);
  h_frame->GetYaxis()->SetRangeUser(0.3,0.4);
  h_frame->Draw("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);
  plotTopLegend("20\%-60\%",1.0,0.385,0.04,1,0.0,42,0);
  for(int i_energy = 0; i_energy < 5; ++i_energy)
  {
    g_sys[i_energy]->SetLineColor(kGray);
    g_sys[i_energy]->SetLineStyle(1);
    g_sys[i_energy]->SetLineWidth(6);
    g_sys[i_energy]->Draw("pE same");
    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_stat[i_energy],plot_style[i_energy],plot_color[i_energy],1.0);
    Draw_TGAE_Point_new_Symbol(3.5,0.385-i_energy*0.005,0.0,0.0,0.0,0.0,plot_style[i_energy],plot_color[i_energy],1.4);
    plotTopLegend((char*)Energy[i_energy].c_str(),3.65,0.385-i_energy*0.005-0.0014,0.04,1,0.0,42,0);
  }

  c_rho->SaveAs("./figures/c_rhoSys.png");
}

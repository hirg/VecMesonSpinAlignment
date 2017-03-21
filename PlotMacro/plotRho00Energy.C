#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TStyle.h"
#include "TF1.h"

using namespace std;

void plotRho00Energy()
{
  gStyle->SetOptDate(0);
  TFile *File_Input[5];
  TGraphAsymmErrors *g_rho_stat[5], *g_rho_sys[5];
  for(int i_energy = 0; i_energy < 5; ++i_energy)
  {
    string inputfile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/rho00/InteRhoSys.root",vmsa::mBeamEnergy[i_energy+2].c_str());
    cout << "Open InPut file: " << inputfile.c_str() << endl;
    File_Input[i_energy] = TFile::Open(inputfile.c_str());
    g_rho_stat[i_energy] = (TGraphAsymmErrors*)File_Input[i_energy]->Get("g_rho_stat");
    g_rho_sys[i_energy] = (TGraphAsymmErrors*)File_Input[i_energy]->Get("g_rho_sys");
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
  // c_rho00->cd()->SetLogx();
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,220.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("#sqrt{s_{NN}}  (GeV)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.301,0.40);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.2);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(0.0,220.0,1.0/3.0,1.0/3.0,1,3,2);
  plotTopLegend((char*)"#rho_{00} = 1/3",100,0.328,0.04,1,0.0,42,0);

  // TBox *sysbox[5];
  TGraphAsymmErrors *g_rho = new TGraphAsymmErrors();
  for(int i_energy = 0; i_energy < 5; ++i_energy)
  {
    double energy, rho;
    g_rho_stat[i_energy]->GetPoint(0,energy,rho);
    double err_stat = g_rho_stat[i_energy]->GetErrorYhigh(0);
    g_rho->SetPoint(i_energy,energy,rho);
    g_rho->SetPointError(i_energy,0.0,0.0,err_stat,err_stat);
  }
  TF1 *f_pol = new TF1("f_pol","pol1",0,300);
  f_pol->SetParameter(0,0.35);
  f_pol->SetParameter(1,-1);
  f_pol->SetRange(18,200);
  g_rho->Fit(f_pol,"MFNR");
  double mean_err = f_pol->GetParError(0);
  TGraphAsymmErrors *g_band = new TGraphAsymmErrors();
  g_band->SetPoint(0,19.6,f_pol->Eval(19.6));
  g_band->SetPointError(0,0.0,0.0,mean_err,mean_err);
  g_band->SetPoint(1,200.0,f_pol->Eval(200.0));
  g_band->SetPointError(1,0.0,0.0,mean_err,mean_err);
  g_band->SetMarkerColor(kYellow);
  g_band->SetLineColor(kYellow);
  g_band->SetFillColor(kYellow);
  g_band->SetFillStyle(3001);
  g_band->Draw("pE3 same");

  for(int i_energy = 0; i_energy < 5; ++i_energy)
  {
    double energy, rho;
    g_rho_sys[i_energy]->GetPoint(0,energy,rho);
    double err = g_rho_sys[i_energy]->GetErrorYhigh(0);

    // sysbox[i_energy] = new TBox(energy*0.95,rho-err,energy*1.05,rho+err);
    // sysbox[i_energy]->SetFillColor(kGray+1);
    // sysbox[i_energy]->SetFillStyle(3001);
    // sysbox[i_energy]->Draw("same");

    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_stat[i_energy],29,kRed,1.8);
    PlotLine(energy-2,energy+2,rho+err,rho+err,kGray+2,2,1);
    PlotLine(energy-2,energy-2,rho+err-0.001,rho+err,kGray+2,2,1);
    PlotLine(energy+2,energy+2,rho+err-0.001,rho+err,kGray+2,2,1);
    PlotLine(energy-2,energy+2,rho-err,rho-err,kGray+2,2,1);
    PlotLine(energy-2,energy-2,rho-err+0.001,rho-err,kGray+2,2,1);
    PlotLine(energy+2,energy+2,rho-err+0.001,rho-err,kGray+2,2,1);
  }
  plotTopLegend((char*)"Au+Au (20-60\%)",14,0.39,0.04,1,0.0,42,0);
  plotTopLegend((char*)"|#eta| < 1",150,0.39,0.04,1,0.0,42,0);
  Draw_TGAE_Point_new_Symbol(14,0.38,0.0,0.0,0.0,0.0,29,kRed,1.8);
  plotTopLegend((char*)"#phi-meson (0.4 < p_{T}< 3.0 GeV/c)",20,0.379,0.04,1,0.0,42,0);
  TBox *avebox = new TBox(6,0.370,18,0.373);
  avebox->SetFillColor(kYellow);
  avebox->SetFillStyle(3001);
  avebox->Draw("same");
  // plotTopLegend((char*)"averaged #rho_{00} (statistical errors)",20,0.37,0.04,1,0.0,42,0);
  plotTopLegend((char*)"pol1 fit (statistical errors)",20,0.37,0.04,1,0.0,42,0);

  f_pol->SetLineColor(kYellow);
  f_pol->SetLineWidth(3);
  f_pol->SetLineStyle(2);
  // f_pol->Draw("lf same");

  // plotTopLegend((char*)"STAR Preliminary",60,0.325,0.04,1,0.0,42,0);

  c_rho00->SaveAs("../figures/rhoSys_energy.eps");
}

#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include "../Utility/StSpinAlignmentCons.h"
// #include "../Utility/type.h"
#include "../Utility/draw.h"
#include "../Utility/functions.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

void calInteRho(int energy = 2, int pid = 0)
{
  string inputspec = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/MonteCarlo/Data/Phi_Spec.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  cout << "Open input spectra: " << inputspec.c_str() << endl;
  TFile *File_Spec = TFile::Open(inputspec.c_str());

  TGraphAsymmErrors *g_spec = (TGraphAsymmErrors*)File_Spec->Get("g_spec");
  TF1 *f_Levy = new TF1("f_Levy",Levy,vmsa::ptMin,vmsa::ptMax,3);
  f_Levy->SetParameter(0,1);
  f_Levy->SetParameter(1,10);
  f_Levy->SetParameter(2,0.1);
  f_Levy->SetLineStyle(2);
  f_Levy->SetLineColor(4);
  f_Levy->SetLineWidth(2);
  g_spec->Fit(f_Levy,"N");

  TF1 *f_spec = new TF1("f_spec",pTLevy,vmsa::ptMin,vmsa::ptMax,3);
  f_spec->SetParameter(0,f_Levy->GetParameter(0));
  f_spec->SetParameter(1,f_Levy->GetParameter(1));
  f_spec->SetParameter(2,f_Levy->GetParameter(2));
  f_spec->SetLineStyle(2);
  f_spec->SetLineColor(2);
  f_spec->SetLineWidth(2);

#if _PlotQA_
  TCanvas *c_spec = new TCanvas("c_spec","c_spec",10,10,800,800);
  c_spec->cd()->SetLeftMargin(0.15);
  c_spec->cd()->SetBottomMargin(0.15);
  c_spec->cd()->SetTicks(1,1);
  c_spec->cd()->SetGrid(0,0);
  c_spec->SetLogy();
  TH1F *h_spec = new TH1F("h_spec","h_spec",100,0.0,10.0);
  for(int i_bin = 1; i_bin < 101; ++i_bin)
  {
    h_spec->SetBinContent(i_bin,-10.0);
    h_spec->SetBinError(i_bin,1.0);
  }
  h_spec->SetTitle("");
  h_spec->SetStats(0);
  h_spec->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_spec->GetXaxis()->CenterTitle();
  h_spec->GetYaxis()->SetTitle("dN/p_{T}dp_{T}");
  h_spec->GetYaxis()->CenterTitle();
  h_spec->GetYaxis()->SetRangeUser(1E-6,10);
  h_spec->Draw("pE");
  g_spec->Draw("pE same");
  f_Levy->Draw("l same");
  f_spec->Draw("l same");
#endif

  string inputrho = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/Rho_SysErrors.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_Rho = TFile::Open(inputrho.c_str());
  TH1F *h_frame = (TH1F*)File_Rho->Get("h_frame");
  string StatErrorRho = Form("g_rho00_%s_%s_StatError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TGraphAsymmErrors *g_StatErrors = (TGraphAsymmErrors*)File_Rho->Get(StatErrorRho.c_str());
  string SysErrorRho = Form("g_rho00_%s_%s_SysError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TGraphAsymmErrors *g_SysErrors = (TGraphAsymmErrors*)File_Rho->Get(SysErrorRho.c_str());

  // integration range 0.4-3.0 GeV/c
  int Nmax = 5;
  double yields = f_spec->Integral(vmsa::pt_low[energy][0],vmsa::pt_low[energy][Nmax]);
  cout << "yields in [0.4,3.0]: " << yields << endl;
  double rho = 0.0;
  double err_stat = 0.0;
  double err_sys = 0.0;
  for(int i_point = 0; i_point < Nmax; ++i_point)
  {
    double pt, rho_pt;
    g_StatErrors->GetPoint(i_point,pt,rho_pt);
    double yields_pt = f_spec->Integral(vmsa::pt_low[energy][i_point],vmsa::pt_up[energy][i_point]);
    double weight = yields_pt/yields;
    rho += rho_pt*weight;
    cout << "pt_bin: [" << vmsa::pt_low[energy][i_point] << "," << vmsa::pt_up[energy][i_point] << "]: rho_pt = " << rho_pt << ", weight = " << 100.0*weight << "%"<< endl;

    double err_stat_pt = g_StatErrors->GetErrorYhigh(i_point);
    err_stat += weight*err_stat_pt*weight*err_stat_pt;

    double err_sys_pt = g_SysErrors->GetErrorYhigh(i_point);
    cout << "err_sys_pt = " << err_sys_pt << endl;
    // err_sys += weight*err_sys_pt*weight*err_sys_pt;
    err_sys += weight*err_sys_pt;
  }
  cout << "integrated rho = " << rho << " +/- " << TMath::Sqrt(err_stat) << " +/- " << err_sys << endl;

#if _PlotQA_
  TCanvas *c_rho_SysError = new TCanvas("c_rho_SysError","c_rho_SysError",600,10,800,800);
  c_rho_SysError->cd();
  c_rho_SysError->cd()->SetLeftMargin(0.15);
  c_rho_SysError->cd()->SetBottomMargin(0.15);
  c_rho_SysError->cd()->SetTicks(1,1);
  c_rho_SysError->cd()->SetGrid(0,0);
  h_frame->Draw("pE");
  g_StatErrors->Draw("pE same");
  g_SysErrors->Draw("pE same");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  TGraphAsymmErrors *g_temp_stat = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_temp_sys  = new TGraphAsymmErrors();

  g_temp_stat->SetPoint(0,3.0,rho);
  g_temp_stat->SetPointError(0,0.0,0.0,TMath::Sqrt(err_stat),TMath::Sqrt(err_stat));
  g_temp_stat->SetMarkerStyle(29);
  g_temp_stat->SetMarkerColor(2);
  g_temp_stat->SetMarkerSize(2);
  g_temp_stat->SetLineColor(2);
  g_temp_stat->Draw("pE same");

  g_temp_sys->SetPoint(0,3.0,rho);
  // g_temp_sys->SetPointError(0,0.0,0.0,TMath::Sqrt(err_sys),TMath::Sqrt(err_sys));
  g_temp_sys->SetPointError(0,0.0,0.0,err_sys,err_sys);
  g_temp_sys->SetMarkerStyle(29);
  g_temp_sys->SetMarkerColor(2);
  g_temp_sys->SetMarkerSize(2);
  g_temp_sys->SetLineColor(4);
  g_temp_sys->Draw("pE same");
#endif

  float val_energy[7] = {7.7,11.5,19.6,27.0,39.0,62.4,200.0};
  string outputfile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/InteRhoSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  cout << "OutPut file set to: " << outputfile.c_str() << endl;
  TFile *File_Out = new TFile(outputfile.c_str(),"RECREATE");
  File_Out->cd();
  TGraphAsymmErrors *g_rho_stat = new TGraphAsymmErrors();
  g_rho_stat->SetPoint(0,val_energy[energy],rho);
  g_rho_stat->SetPointError(0,0.0,0.0,TMath::Sqrt(err_stat),TMath::Sqrt(err_stat));
  g_rho_stat->SetName("g_rho_stat");
  g_rho_stat->SetMarkerStyle(29);
  g_rho_stat->SetMarkerColor(2);
  g_rho_stat->SetMarkerSize(2);
  g_rho_stat->SetLineColor(2);
  g_rho_stat->Write();

  TGraphAsymmErrors *g_rho_sys  = new TGraphAsymmErrors();
  g_rho_sys->SetPoint(0,val_energy[energy],rho);
  // g_rho_sys->SetPointError(0,0.0,0.0,TMath::Sqrt(err_sys),TMath::Sqrt(err_sys));
  g_rho_sys->SetPointError(0,0.0,0.0,err_sys,err_sys);
  g_rho_sys->SetName("g_rho_sys");
  g_rho_sys->SetMarkerStyle(29);
  g_rho_sys->SetMarkerColor(kGray+2);
  g_rho_sys->SetMarkerSize(2);
  g_rho_sys->SetLineColor(kGray+2);
  g_rho_sys->Write();
  File_Out->Close();
}

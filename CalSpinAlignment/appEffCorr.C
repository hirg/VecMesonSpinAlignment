#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TF1.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
#include "../Utility/draw.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

// used for 20%-60% and eta_gap = 0.05 only
void appEffCorr(int energy = 6, int pid = 0)
{
  string InPutRho = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/RawRhoPt.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_Rho = TFile::Open(InPutRho.c_str());

  TH1FMap h_mCounts;
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    string KEY_Count = Form("Count_pt_%d_Centrality_0_EtaGap_0_2nd_%s_SM",i_pt,vmsa::mPID[pid].c_str());
    h_mCounts[KEY_Count] = (TH1F*)File_Rho->Get(KEY_Count.c_str());
    string KEY_BW = Form("BW_pt_%d_Centrality_0_EtaGap_0_2nd_%s_SM",i_pt,vmsa::mPID[pid].c_str());
    h_mCounts[KEY_BW] = (TH1F*)File_Rho->Get(KEY_BW.c_str());
  }

  string InPutEff = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Embedding/%s/Efficiency/Eff_%s_SingleKaon.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Eff = TFile::Open(InPutEff.c_str());
  TH1DMap h_mEff;
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt)
  {
    string KEY_Eff = Form("h_mEffCos_Cent_9_Pt_%d",i_pt); // 20%-60%
    h_mEff[KEY_Eff] = (TH1D*)File_Eff->Get(KEY_Eff.c_str());
  }


#if _PlotQA_
  TCanvas *c_raw = new TCanvas("c_raw","c_raw",10,10,1000,500);
  c_raw->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_raw->cd(i_pad+1);
    c_raw->cd(i_pad+1)->SetLeftMargin(0.15);
    c_raw->cd(i_pad+1)->SetBottomMargin(0.15);
    c_raw->cd(i_pad+1)->SetGrid(0,0);
    c_raw->cd(i_pad+1)->SetTicks(1,1);
  }
  c_raw->cd(1);
  string QA_Count = Form("Count_pt_%d_Centrality_0_EtaGap_0_2nd_%s_SM",vmsa::pt_QA[energy],vmsa::mPID[pid].c_str());
  h_mCounts[QA_Count]->SetStats(0);
  h_mCounts[QA_Count]->SetTitle("20-60%");
  h_mCounts[QA_Count]->SetTitleSize(0.08);
  h_mCounts[QA_Count]->SetLineColor(4);
  h_mCounts[QA_Count]->SetMarkerColor(4);
  h_mCounts[QA_Count]->SetMarkerStyle(24);
  h_mCounts[QA_Count]->SetMarkerSize(0.8);
  h_mCounts[QA_Count]->GetXaxis()->SetTitle("cos(#theta^{*}) (w.r.t. 2^{nd} event plane)");
  h_mCounts[QA_Count]->GetXaxis()->SetTitleSize(0.05);
  h_mCounts[QA_Count]->GetXaxis()->CenterTitle();
  h_mCounts[QA_Count]->GetYaxis()->SetRangeUser(0.9*h_mCounts[QA_Count]->GetMinimum(),1.1*h_mCounts[QA_Count]->GetMaximum());
  h_mCounts[QA_Count]->DrawCopy("pE");
  string QA_BW = Form("BW_pt_%d_Centrality_0_EtaGap_0_2nd_%s_SM",vmsa::pt_QA[energy],vmsa::mPID[pid].c_str());
  h_mCounts[QA_BW]->SetLineColor(2);
  h_mCounts[QA_BW]->SetMarkerColor(2);
  h_mCounts[QA_BW]->SetMarkerStyle(24);
  h_mCounts[QA_BW]->SetMarkerSize(0.8);
  h_mCounts[QA_BW]->DrawCopy("pE same");

  TLegend *leg_temp = new TLegend(0.20,0.2,0.45,0.4);
  leg_temp->SetFillColor(10);
  leg_temp->SetBorderSize(0.0);
  leg_temp->AddEntry(h_mCounts[QA_Count],"bin counting","p");
  leg_temp->AddEntry(h_mCounts[QA_BW],"breit wigner","p");
  leg_temp->Draw("same");

  c_raw->cd(2);
  string QA_Eff = Form("h_mEffCos_Cent_9_Pt_%d",vmsa::pt_QA[energy]); // 20%-60%
  h_mEff[QA_Eff]->SetStats(0);
  h_mEff[QA_Eff]->SetTitle("TPC efficiency 20-60%");
  h_mEff[QA_Eff]->SetTitleSize(0.08);
  h_mEff[QA_Eff]->SetLineColor(kGray+2);
  h_mEff[QA_Eff]->SetMarkerColor(kGray+2);
  h_mEff[QA_Eff]->SetMarkerStyle(24);
  h_mEff[QA_Eff]->SetMarkerSize(0.8);
  h_mEff[QA_Eff]->GetXaxis()->SetTitle("cos(#theta^{*})");
  h_mEff[QA_Eff]->GetXaxis()->SetTitleSize(0.05);
  h_mEff[QA_Eff]->GetXaxis()->CenterTitle();
  h_mEff[QA_Eff]->GetYaxis()->SetRangeUser(0.9*h_mEff[QA_Eff]->GetMinimum(),1.1*h_mEff[QA_Eff]->GetMaximum());
  h_mEff[QA_Eff]->DrawCopy("pE");
  TF1 *f_polQA = new TF1("f_polQA","pol0",0,1);
  h_mEff[QA_Eff]->Fit(f_polQA,"N");
  f_polQA->SetLineColor(kGray+2);
  f_polQA->SetLineStyle(2);
  f_polQA->SetLineWidth(2);
  f_polQA->Draw("l same");
  float effQA = f_polQA->GetParameter(0);
  float err_effQA = f_polQA->GetParError(0);
  float chi2QA = f_polQA->GetChisquare();
  float NDFQA = f_polQA->GetNDF();
  string leg_effQA = Form("eff = %2.4f #pm %2.4f, #chi^{2}/NDF = %2.2f/%2.2f",effQA,err_effQA,chi2QA,NDFQA);
  plotTopLegend((char*)leg_effQA.c_str(),0.2,0.8,0.04,1,0.0,42,1);
  string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][vmsa::pt_QA[energy]],vmsa::pt_up[energy][vmsa::pt_QA[energy]]);
  plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.04,1,0.0,42,1);
#endif

}

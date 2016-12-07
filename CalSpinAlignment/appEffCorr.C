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
#include "../Utility/functions.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

using namespace std;

int const Color_Counts = 2;
int const Color_BW     = 4;

// TH1F* CalEffError(TH1F *h_counts, TH1D *h_eff, std::string HistName)
// {
//   TH1F* h_ratio = (TH1F*)h_Rc->Clone();
//   h_ratio->Divide(h_Mc);
//   for(int i_bin = 1; i_bin < h_ratio->GetNbinsX()+1; ++i_bin)
//   {
//     double n = h_Mc->GetBinContent(i_bin);
//     double k = h_Rc->GetBinContent(i_bin);
//     double variance = (k+1.0)*(k+2.0)/((n+2.0)*(n+3.0))-(k+1.0)*(k+1.0)/((n+2.0)*(n+2.0));
//     double sigma = TMath::Sqrt(variance);
//     if(n > 0.0 && k > 0.0) h_ratio->SetBinError(i_bin,sigma);
//   }
//   h_ratio->SetName(HistName.c_str());
//
//   return h_ratio;
// }

// used for 20%-60% and eta_gap = 0.05 only
void appEffCorr(int energy = 6, int pid = 0)
{
  string InPutRho = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/RawRhoPt.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_Rho = TFile::Open(InPutRho.c_str());

  string InPutEff = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Embedding/%s/Efficiency/Eff_%s_SingleKaon.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Eff = TFile::Open(InPutEff.c_str());

  TH1FMap h_mCounts, h_mCountEff;
  TH1DMap h_mEff;
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    string KEY_Count = Form("Count_pt_%d_Centrality_0_EtaGap_0_2nd_%s_SM",i_pt,vmsa::mPID[pid].c_str());
    h_mCounts[KEY_Count] = (TH1F*)File_Rho->Get(KEY_Count.c_str());
    h_mCountEff[KEY_Count] = (TH1F*)h_mCounts[KEY_Count]->Clone();

    string KEY_BW = Form("BW_pt_%d_Centrality_0_EtaGap_0_2nd_%s_SM",i_pt,vmsa::mPID[pid].c_str());
    h_mCounts[KEY_BW] = (TH1F*)File_Rho->Get(KEY_BW.c_str());
    h_mCountEff[KEY_BW] = (TH1F*)h_mCounts[KEY_BW]->Clone();

    string KEY_Eff = Form("h_mEffCos_Cent_9_Pt_%d",i_pt); // 20%-60%
    h_mEff[KEY_Eff] = (TH1D*)File_Eff->Get(KEY_Eff.c_str());

    h_mCountEff[KEY_Count]->Divide(h_mEff[KEY_Eff]);
    h_mCountEff[KEY_BW]->Divide(h_mEff[KEY_Eff]);
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
  h_mCounts[QA_Count]->SetLineColor(Color_Counts);
  h_mCounts[QA_Count]->SetMarkerColor(Color_Counts);
  h_mCounts[QA_Count]->SetMarkerStyle(24);
  h_mCounts[QA_Count]->SetMarkerSize(0.8);
  h_mCounts[QA_Count]->GetXaxis()->SetTitle("cos(#theta^{*}) (w.r.t. 2^{nd} event plane)");
  h_mCounts[QA_Count]->GetXaxis()->SetTitleSize(0.05);
  h_mCounts[QA_Count]->GetXaxis()->CenterTitle();
  h_mCounts[QA_Count]->GetXaxis()->SetNdivisions(505);
  h_mCounts[QA_Count]->GetXaxis()->SetLabelSize(0.03);
  h_mCounts[QA_Count]->GetYaxis()->SetTitle("Counts");
  h_mCounts[QA_Count]->GetYaxis()->SetTitleSize(0.05);
  h_mCounts[QA_Count]->GetYaxis()->CenterTitle();
  h_mCounts[QA_Count]->GetYaxis()->SetNdivisions(505);
  h_mCounts[QA_Count]->GetYaxis()->SetLabelSize(0.03);
  h_mCounts[QA_Count]->GetYaxis()->SetRangeUser(0.9*h_mCounts[QA_Count]->GetMinimum(),1.1*h_mCounts[QA_Count]->GetMaximum());
  h_mCounts[QA_Count]->DrawCopy("pE");
  TF1 *f_CountsQA = new TF1("f_CountsQA",SpinDensity,0,1.0,2);
  f_CountsQA->SetParameter(0,0.33);
  f_CountsQA->SetParameter(1,100);
  h_mCounts[QA_Count]->Fit(f_CountsQA,"N");
  f_CountsQA->SetLineColor(Color_Counts);
  f_CountsQA->SetLineWidth(2);
  f_CountsQA->SetLineStyle(2);
  f_CountsQA->Draw("l Same");
  float rhoC = f_CountsQA->GetParameter(0);
  float err_rhoC = f_CountsQA->GetParError(0);
  float chi2C = f_CountsQA->GetChisquare();
  float NDFC = f_CountsQA->GetNDF();
  string leg_rhoC = Form("#rho_{00} = %2.4f #pm %2.4f, #chi^{2}/NDF = %2.2f/%2.2f",rhoC,err_rhoC,chi2C,NDFC);
  plotTopLegend((char*)leg_rhoC.c_str(),0.2,0.8,0.04,Color_Counts,0.0,42,1);

  string QA_BW = Form("BW_pt_%d_Centrality_0_EtaGap_0_2nd_%s_SM",vmsa::pt_QA[energy],vmsa::mPID[pid].c_str());
  h_mCounts[QA_BW]->SetLineColor(Color_BW);
  h_mCounts[QA_BW]->SetMarkerColor(Color_BW);
  h_mCounts[QA_BW]->SetMarkerStyle(24);
  h_mCounts[QA_BW]->SetMarkerSize(0.8);
  h_mCounts[QA_BW]->DrawCopy("pE same");
  TF1 *f_BWQA= new TF1("f_BWQA",SpinDensity,0,1.0,2);
  f_BWQA->SetParameter(0,0.33);
  f_BWQA->SetParameter(1,100);
  h_mCounts[QA_BW]->Fit(f_BWQA,"N");
  f_BWQA->SetLineColor(Color_BW);
  f_BWQA->SetLineWidth(2);
  f_BWQA->SetLineStyle(2);
  f_BWQA->Draw("l Same");
  float rhoB = f_BWQA->GetParameter(0);
  float err_rhoB = f_BWQA->GetParError(0);
  float chi2B = f_BWQA->GetChisquare();
  float NDFB = f_BWQA->GetNDF();
  string leg_rhoB = Form("#rho_{00} = %2.4f #pm %2.4f, #chi^{2}/NDF = %2.2f/%2.2f",rhoB,err_rhoB,chi2B,NDFB);
  plotTopLegend((char*)leg_rhoB.c_str(),0.2,0.75,0.04,Color_BW,0.0,42,1);

  TLegend *leg_rhoQA = new TLegend(0.20,0.2,0.45,0.4);
  leg_rhoQA->SetFillColor(10);
  leg_rhoQA->SetBorderSize(0.0);
  leg_rhoQA->AddEntry(h_mCounts[QA_Count],"bin counting","p");
  leg_rhoQA->AddEntry(h_mCounts[QA_BW],"breit wigner","p");
  leg_rhoQA->Draw("same");

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

  TCanvas *c_Eff = new TCanvas("c_Eff","c_Eff",10,10,800,800);
  c_Eff->cd();
  c_Eff->SetLeftMargin(0.15);
  c_Eff->SetBottomMargin(0.15);
  c_Eff->SetTicks(1,1);
  c_Eff->SetGrid(0,0);
  h_mCountEff[QA_Count]->SetStats(0);
  h_mCountEff[QA_Count]->SetTitle("20-60%");
  h_mCountEff[QA_Count]->SetTitleSize(0.08);
  h_mCountEff[QA_Count]->SetLineColor(Color_Counts);
  h_mCountEff[QA_Count]->SetMarkerColor(Color_Counts);
  h_mCountEff[QA_Count]->SetMarkerStyle(24);
  h_mCountEff[QA_Count]->SetMarkerSize(0.8);
  h_mCountEff[QA_Count]->GetXaxis()->SetTitle("cos(#theta^{*}) (w.r.t. 2^{nd} event plane)");
  h_mCountEff[QA_Count]->GetXaxis()->SetTitleSize(0.05);
  h_mCountEff[QA_Count]->GetXaxis()->CenterTitle();
  h_mCountEff[QA_Count]->GetXaxis()->SetNdivisions(505);
  h_mCountEff[QA_Count]->GetXaxis()->SetLabelSize(0.03);
  h_mCountEff[QA_Count]->GetYaxis()->SetTitle("Counts");
  h_mCountEff[QA_Count]->GetYaxis()->SetTitleSize(0.05);
  h_mCountEff[QA_Count]->GetYaxis()->CenterTitle();
  h_mCountEff[QA_Count]->GetYaxis()->SetNdivisions(505);
  h_mCountEff[QA_Count]->GetYaxis()->SetLabelSize(0.03);
  h_mCountEff[QA_Count]->GetYaxis()->SetRangeUser(0.9*h_mCountEff[QA_Count]->GetMinimum(),1.1*h_mCountEff[QA_Count]->GetMaximum());
  h_mCountEff[QA_Count]->DrawCopy("pE");
  TF1 *f_CeffQA = new TF1("f_CeffQA",SpinDensity,0,1.0,2);
  f_CeffQA->SetParameter(0,0.33);
  f_CeffQA->SetParameter(1,100);
  h_mCountEff[QA_Count]->Fit(f_CeffQA,"N");
  f_CeffQA->SetLineColor(Color_Counts);
  f_CeffQA->SetLineWidth(2);
  f_CeffQA->SetLineStyle(2);
  f_CeffQA->Draw("l Same");
  float rhoCeff = f_CeffQA->GetParameter(0);
  float err_rhoCeff = f_CeffQA->GetParError(0);
  float chi2Ceff = f_CeffQA->GetChisquare();
  float NDFCeff = f_CeffQA->GetNDF();
  string leg_rhoCeff = Form("#rho_{00} = %2.4f #pm %2.4f, #chi^{2}/NDF = %2.2f/%2.2f",rhoCeff,err_rhoCeff,chi2Ceff,NDFCeff);
  plotTopLegend((char*)leg_rhoCeff.c_str(),0.2,0.8,0.04,Color_Counts,0.0,42,1);

  h_mCountEff[QA_BW]->SetLineColor(Color_BW);
  h_mCountEff[QA_BW]->SetMarkerColor(Color_BW);
  h_mCountEff[QA_BW]->SetMarkerStyle(24);
  h_mCountEff[QA_BW]->SetMarkerSize(0.8);
  h_mCountEff[QA_BW]->DrawCopy("pE same");
  TF1 *f_BeffQA= new TF1("f_BeffQA",SpinDensity,0,1.0,2);
  f_BeffQA->SetParameter(0,0.33);
  f_BeffQA->SetParameter(1,100);
  h_mCountEff[QA_BW]->Fit(f_BeffQA,"N");
  f_BeffQA->SetLineColor(Color_BW);
  f_BeffQA->SetLineWidth(2);
  f_BeffQA->SetLineStyle(2);
  f_BeffQA->Draw("l Same");
  float rhoBeff = f_BeffQA->GetParameter(0);
  float err_rhoBeff = f_BeffQA->GetParError(0);
  float chi2Beff = f_BeffQA->GetChisquare();
  float NDFBeff = f_BeffQA->GetNDF();
  string leg_rhoBeff = Form("#rho_{00} = %2.4f #pm %2.4f, #chi^{2}/NDF = %2.2f/%2.2f",rhoBeff,err_rhoBeff,chi2Beff,NDFBeff);
  plotTopLegend((char*)leg_rhoBeff.c_str(),0.2,0.75,0.04,Color_BW,0.0,42,1);

  TLegend *leg_eff = new TLegend(0.20,0.2,0.45,0.4);
  leg_eff->SetFillColor(10);
  leg_eff->SetBorderSize(0.0);
  leg_eff->AddEntry(h_mCountEff[QA_Count],"bin counting","p");
  leg_eff->AddEntry(h_mCountEff[QA_BW],"breit wigner","p");
  leg_eff->Draw("same");
#endif

}

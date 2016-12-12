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

int const Color_Counts = kGray+2;
int const Color_BW     = kAzure-4;

TH1F* CalEffCorr(TH1F *h_counts, TH1D *h_eff, std::string HistName)
{
  TH1F* h_effCorr = (TH1F*)h_counts->Clone();
  h_effCorr->Divide(h_eff);
  for(int i_bin = 1; i_bin < h_effCorr->GetNbinsX()+1; ++i_bin)
  {
    float counts = h_counts->GetBinContent(i_bin);
    float err_counts = h_counts->GetBinError(i_bin);
    float eff = h_eff->GetBinContent(i_bin);
    float err_eff = h_eff->GetBinError(i_bin);
    float err = ErrDiv(counts,eff,err_counts,err_eff);
    h_effCorr->SetBinError(i_bin,err);
  }
  h_effCorr->SetName(HistName.c_str());

  return h_effCorr;
}

TH1F* CalEffCorr(TH1F *h_counts, TF1 *f_eff, std::string HistName)
{
  TH1F* h_effCorr = (TH1F*)h_counts->Clone();
  h_effCorr->Divide(f_eff);
  float eff = f_eff->GetParameter(0);
  float err_eff = f_eff->GetParError(0);
  for(int i_bin = 1; i_bin < h_effCorr->GetNbinsX()+1; ++i_bin)
  {
    float counts = h_counts->GetBinContent(i_bin);
    float err_counts = h_counts->GetBinError(i_bin);
    float err = ErrDiv(counts,eff,err_counts,err_eff);
    h_effCorr->SetBinError(i_bin,err);
  }
  h_effCorr->SetName(HistName.c_str());
  h_effCorr->SetTitle(HistName.c_str());

  return h_effCorr;
}

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
    string KEY_Eff = Form("h_mEffCos_Cent_9_Pt_%d",i_pt); // 20%-60%
    h_mEff[KEY_Eff] = (TH1D*)File_Eff->Get(KEY_Eff.c_str());
    TF1 *f_pol = new TF1("f_pol","pol0",0.0,1.0);
    h_mEff[KEY_Eff]->Fit(f_pol,"N");

    string KEY_Count = Form("Count_pt_%d_Centrality_0_EtaGap_0_2nd_%s_SM",i_pt,vmsa::mPID[pid].c_str());
    h_mCounts[KEY_Count] = (TH1F*)File_Rho->Get(KEY_Count.c_str());
    string KEY_CPoint = Form("CPoint_pt_%d_Centrality_0_EtaGap_0_%s",i_pt,vmsa::mPID[pid].c_str());
    h_mCountEff[KEY_CPoint] = CalEffCorr(h_mCounts[KEY_Count],h_mEff[KEY_Eff],KEY_CPoint);
    string KEY_CPoly = Form("CPoly_pt_%d_Centrality_0_EtaGap_0_%s",i_pt,vmsa::mPID[pid].c_str());
    h_mCountEff[KEY_CPoly] = CalEffCorr(h_mCounts[KEY_Count],f_pol,KEY_CPoly);

    string KEY_BW = Form("BW_pt_%d_Centrality_0_EtaGap_0_2nd_%s_SM",i_pt,vmsa::mPID[pid].c_str());
    h_mCounts[KEY_BW] = (TH1F*)File_Rho->Get(KEY_BW.c_str());
    string KEY_BPoint = Form("BPoint_pt_%d_Centrality_0_EtaGap_0_%s",i_pt,vmsa::mPID[pid].c_str());
    h_mCountEff[KEY_BPoint] = CalEffCorr(h_mCounts[KEY_BW],h_mEff[KEY_Eff],KEY_BPoint);
    string KEY_BPoly = Form("BPoly_pt_%d_Centrality_0_EtaGap_0_%s",i_pt,vmsa::mPID[pid].c_str());
    h_mCountEff[KEY_BPoly] = CalEffCorr(h_mCounts[KEY_BW],f_pol,KEY_BPoly);
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
#endif

#if _PlotQA_
  TCanvas *c_Eff = new TCanvas("c_Eff","c_Eff",10,10,1000,500);
  c_Eff->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_Eff->cd(i_pad+1);
    c_Eff->cd(i_pad+1)->SetLeftMargin(0.15);
    c_Eff->cd(i_pad+1)->SetBottomMargin(0.15);
    c_Eff->cd(i_pad+1)->SetGrid(0,0);
    c_Eff->cd(i_pad+1)->SetTicks(1,1);
  }
  c_Eff->cd(1);
  string QA_CPoint = Form("CPoint_pt_%d_Centrality_0_EtaGap_0_%s",vmsa::pt_QA[energy],vmsa::mPID[pid].c_str());
  h_mCountEff[QA_CPoint]->SetStats(0);
  h_mCountEff[QA_CPoint]->SetTitle("20-60%");
  h_mCountEff[QA_CPoint]->SetTitleSize(0.08);
  h_mCountEff[QA_CPoint]->SetLineColor(Color_Counts);
  h_mCountEff[QA_CPoint]->SetMarkerColor(Color_Counts);
  h_mCountEff[QA_CPoint]->SetMarkerStyle(24);
  h_mCountEff[QA_CPoint]->SetMarkerSize(0.8);
  h_mCountEff[QA_CPoint]->GetXaxis()->SetTitle("cos(#theta^{*}) (w.r.t. 2^{nd} event plane)");
  h_mCountEff[QA_CPoint]->GetXaxis()->SetTitleSize(0.05);
  h_mCountEff[QA_CPoint]->GetXaxis()->CenterTitle();
  h_mCountEff[QA_CPoint]->GetXaxis()->SetNdivisions(505);
  h_mCountEff[QA_CPoint]->GetXaxis()->SetLabelSize(0.03);
  h_mCountEff[QA_CPoint]->GetYaxis()->SetTitle("Counts");
  h_mCountEff[QA_CPoint]->GetYaxis()->SetTitleSize(0.05);
  h_mCountEff[QA_CPoint]->GetYaxis()->CenterTitle();
  h_mCountEff[QA_CPoint]->GetYaxis()->SetNdivisions(505);
  h_mCountEff[QA_CPoint]->GetYaxis()->SetLabelSize(0.03);
  h_mCountEff[QA_CPoint]->GetYaxis()->SetRangeUser(0.9*h_mCountEff[QA_CPoint]->GetMinimum(),1.1*h_mCountEff[QA_CPoint]->GetMaximum());
  h_mCountEff[QA_CPoint]->DrawCopy("pE");
  TF1 *f_CPointQA = new TF1("f_CPointQA",SpinDensity,0,1.0,2);
  f_CPointQA->SetParameter(0,0.33);
  f_CPointQA->SetParameter(1,100);
  h_mCountEff[QA_CPoint]->Fit(f_CPointQA,"N");
  f_CPointQA->SetLineColor(Color_Counts);
  f_CPointQA->SetLineWidth(2);
  f_CPointQA->SetLineStyle(2);
  f_CPointQA->Draw("l Same");
  float rhoCPoint = f_CPointQA->GetParameter(0);
  float err_rhoCPoint = f_CPointQA->GetParError(0);
  float chi2CPoint = f_CPointQA->GetChisquare();
  float NDFCPoint = f_CPointQA->GetNDF();
  string leg_rhoCPoint = Form("#rho_{00} = %2.4f #pm %2.4f, #chi^{2}/NDF = %2.2f/%2.2f",rhoCPoint,err_rhoCPoint,chi2CPoint,NDFCPoint);
  plotTopLegend((char*)leg_rhoCPoint.c_str(),0.2,0.8,0.04,Color_Counts,0.0,42,1);

  string QA_BPoint = Form("BPoint_pt_%d_Centrality_0_EtaGap_0_%s",vmsa::pt_QA[energy],vmsa::mPID[pid].c_str());
  h_mCountEff[QA_BPoint]->SetLineColor(Color_BW);
  h_mCountEff[QA_BPoint]->SetMarkerColor(Color_BW);
  h_mCountEff[QA_BPoint]->SetMarkerStyle(24);
  h_mCountEff[QA_BPoint]->SetMarkerSize(0.8);
  h_mCountEff[QA_BPoint]->DrawCopy("pE same");
  TF1 *f_BPointQA = new TF1("f_BPointQA",SpinDensity,0,1.0,2);
  f_BPointQA->SetParameter(0,0.33);
  f_BPointQA->SetParameter(1,100);
  h_mCountEff[QA_BPoint]->Fit(f_BPointQA,"N");
  f_BPointQA->SetLineColor(Color_BW);
  f_BPointQA->SetLineWidth(2);
  f_BPointQA->SetLineStyle(2);
  f_BPointQA->Draw("l Same");
  float rhoBPoint = f_BPointQA->GetParameter(0);
  float err_rhoBPoint = f_BPointQA->GetParError(0);
  float chi2BPoint = f_BPointQA->GetChisquare();
  float NDFBPoint = f_BPointQA->GetNDF();
  string leg_rhoBPoint = Form("#rho_{00} = %2.4f #pm %2.4f, #chi^{2}/NDF = %2.2f/%2.2f",rhoBPoint,err_rhoBPoint,chi2BPoint,NDFBPoint);
  plotTopLegend((char*)leg_rhoBPoint.c_str(),0.2,0.75,0.04,Color_BW,0.0,42,1);

  TLegend *leg_Point = new TLegend(0.20,0.2,0.45,0.4);
  leg_Point->SetFillColor(10);
  leg_Point->SetBorderSize(0.0);
  leg_Point->AddEntry(h_mCountEff[QA_CPoint],"bin counting","p");
  leg_Point->AddEntry(h_mCountEff[QA_BPoint],"breit wigner","p");
  leg_Point->Draw("same");

  c_Eff->cd(2);
  string QA_CPoly = Form("CPoly_pt_%d_Centrality_0_EtaGap_0_%s",vmsa::pt_QA[energy],vmsa::mPID[pid].c_str());
  h_mCountEff[QA_CPoly]->SetStats(0);
  h_mCountEff[QA_CPoly]->SetTitle("20-60%");
  h_mCountEff[QA_CPoly]->SetTitleSize(0.08);
  h_mCountEff[QA_CPoly]->SetLineColor(Color_Counts);
  h_mCountEff[QA_CPoly]->SetMarkerColor(Color_Counts);
  h_mCountEff[QA_CPoly]->SetMarkerStyle(24);
  h_mCountEff[QA_CPoly]->SetMarkerSize(0.8);
  h_mCountEff[QA_CPoly]->GetXaxis()->SetTitle("cos(#theta^{*}) (w.r.t. 2^{nd} event plane)");
  h_mCountEff[QA_CPoly]->GetXaxis()->SetTitleSize(0.05);
  h_mCountEff[QA_CPoly]->GetXaxis()->CenterTitle();
  h_mCountEff[QA_CPoly]->GetXaxis()->SetNdivisions(505);
  h_mCountEff[QA_CPoly]->GetXaxis()->SetLabelSize(0.03);
  h_mCountEff[QA_CPoly]->GetYaxis()->SetTitle("Counts");
  h_mCountEff[QA_CPoly]->GetYaxis()->SetTitleSize(0.05);
  h_mCountEff[QA_CPoly]->GetYaxis()->CenterTitle();
  h_mCountEff[QA_CPoly]->GetYaxis()->SetNdivisions(505);
  h_mCountEff[QA_CPoly]->GetYaxis()->SetLabelSize(0.03);
  h_mCountEff[QA_CPoly]->GetYaxis()->SetRangeUser(0.9*h_mCountEff[QA_CPoly]->GetMinimum(),1.1*h_mCountEff[QA_CPoly]->GetMaximum());
  h_mCountEff[QA_CPoly]->DrawCopy("pE");
  TF1 *f_CPolyQA = new TF1("f_CPolyQA",SpinDensity,0,1.0,2);
  f_CPolyQA->SetParameter(0,0.33);
  f_CPolyQA->SetParameter(1,100);
  h_mCountEff[QA_CPoly]->Fit(f_CPolyQA,"N");
  f_CPolyQA->SetLineColor(Color_Counts);
  f_CPolyQA->SetLineWidth(2);
  f_CPolyQA->SetLineStyle(2);
  f_CPolyQA->Draw("l Same");
  float rhoCPoly = f_CPolyQA->GetParameter(0);
  float err_rhoCPoly = f_CPolyQA->GetParError(0);
  float chi2CPoly = f_CPolyQA->GetChisquare();
  float NDFCPoly = f_CPolyQA->GetNDF();
  string leg_rhoCPoly = Form("#rho_{00} = %2.4f #pm %2.4f, #chi^{2}/NDF = %2.2f/%2.2f",rhoCPoly,err_rhoCPoly,chi2CPoly,NDFCPoly);
  plotTopLegend((char*)leg_rhoCPoly.c_str(),0.2,0.8,0.04,Color_Counts,0.0,42,1);

  string QA_BPoly = Form("BPoly_pt_%d_Centrality_0_EtaGap_0_%s",vmsa::pt_QA[energy],vmsa::mPID[pid].c_str());
  h_mCountEff[QA_BPoly]->SetLineColor(Color_BW);
  h_mCountEff[QA_BPoly]->SetMarkerColor(Color_BW);
  h_mCountEff[QA_BPoly]->SetMarkerStyle(24);
  h_mCountEff[QA_BPoly]->SetMarkerSize(0.8);
  h_mCountEff[QA_BPoly]->DrawCopy("pE same");
  TF1 *f_BPolyQA = new TF1("f_BPolyQA",SpinDensity,0,1.0,2);
  f_BPolyQA->SetParameter(0,0.33);
  f_BPolyQA->SetParameter(1,100);
  h_mCountEff[QA_BPoly]->Fit(f_BPolyQA,"N");
  f_BPolyQA->SetLineColor(Color_BW);
  f_BPolyQA->SetLineWidth(2);
  f_BPolyQA->SetLineStyle(2);
  f_BPolyQA->Draw("l Same");
  float rhoBPoly = f_BPolyQA->GetParameter(0);
  float err_rhoBPoly = f_BPolyQA->GetParError(0);
  float chi2BPoly = f_BPolyQA->GetChisquare();
  float NDFBPoly = f_BPolyQA->GetNDF();
  string leg_rhoBPoly = Form("#rho_{00} = %2.4f #pm %2.4f, #chi^{2}/NDF = %2.2f/%2.2f",rhoBPoly,err_rhoBPoly,chi2BPoly,NDFBPoly);
  plotTopLegend((char*)leg_rhoBPoly.c_str(),0.2,0.75,0.04,Color_BW,0.0,42,1);

  TLegend *leg_Poly = new TLegend(0.20,0.2,0.45,0.4);
  leg_Poly->SetFillColor(10);
  leg_Poly->SetBorderSize(0.0);
  leg_Poly->AddEntry(h_mCountEff[QA_CPoly],"bin counting","p");
  leg_Poly->AddEntry(h_mCountEff[QA_BPoly],"breit wigner","p");
  leg_Poly->Draw("same");
#endif

  TGraMap g_mRho00;
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    {
      string KEY_Rho00_CPoint = Form("Rho00_CPoint_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00[KEY_Rho00_CPoint] = new TGraphAsymmErrors();

      string KEY_Rho00_CPoly = Form("Rho00_CPoly_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00[KEY_Rho00_CPoly] = new TGraphAsymmErrors();

      string KEY_Rho00_BPoint = Form("Rho00_BPoint_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00[KEY_Rho00_BPoint] = new TGraphAsymmErrors();

      string KEY_Rho00_BPoly = Form("Rho00_BPoly_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00[KEY_Rho00_BPoly] = new TGraphAsymmErrors();

      for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++)
      {
	float pt_mean = (vmsa::pt_low[energy][i_pt]+vmsa::pt_up[energy][i_pt])/2.0;

	TF1 *f_CPoint = new TF1("f_CPoint",SpinDensity,0,1.0,2);
	f_CPoint->SetParameter(0,0.33);
	f_CPoint->SetParameter(1,100);
	string KEY_CPoint = Form("CPoint_pt_%d_Centrality_%d_EtaGap_%d_%s",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str());
	h_mCountEff[KEY_CPoint]->Fit(f_CPoint,"N");
	g_mRho00[KEY_Rho00_CPoint]->SetPoint(i_pt,pt_mean,f_CPoint->GetParameter(0));
	g_mRho00[KEY_Rho00_CPoint]->SetPointError(i_pt,0.0,0.0,f_CPoint->GetParError(0),f_CPoint->GetParError(0));
	g_mRho00[KEY_Rho00_CPoint]->SetName(KEY_Rho00_CPoint.c_str());

	TF1 *f_CPoly = new TF1("f_CPoly",SpinDensity,0,1.0,2);
	f_CPoly->SetParameter(0,0.33);
	f_CPoly->SetParameter(1,100);
	string KEY_CPoly = Form("CPoly_pt_%d_Centrality_%d_EtaGap_%d_%s",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str());
	h_mCountEff[KEY_CPoly]->Fit(f_CPoly,"N");
	g_mRho00[KEY_Rho00_CPoly]->SetPoint(i_pt,pt_mean,f_CPoly->GetParameter(0));
	g_mRho00[KEY_Rho00_CPoly]->SetPointError(i_pt,0.0,0.0,f_CPoly->GetParError(0),f_CPoly->GetParError(0));
	g_mRho00[KEY_Rho00_CPoly]->SetName(KEY_Rho00_CPoly.c_str());

	TF1 *f_BPoint = new TF1("f_BPoint",SpinDensity,0,1.0,2);
	f_BPoint->SetParameter(0,0.33);
	f_BPoint->SetParameter(1,100);
	string KEY_BPoint = Form("BPoint_pt_%d_Centrality_%d_EtaGap_%d_%s",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str());
	h_mCountEff[KEY_BPoint]->Fit(f_BPoint,"N");
	g_mRho00[KEY_Rho00_BPoint]->SetPoint(i_pt,pt_mean,f_BPoint->GetParameter(0));
	g_mRho00[KEY_Rho00_BPoint]->SetPointError(i_pt,0.0,0.0,f_BPoint->GetParError(0),f_BPoint->GetParError(0));
	g_mRho00[KEY_Rho00_BPoint]->SetName(KEY_Rho00_BPoint.c_str());

	TF1 *f_BPoly = new TF1("f_BPoly",SpinDensity,0,1.0,2);
	f_BPoly->SetParameter(0,0.33);
	f_BPoly->SetParameter(1,100);
	string KEY_BPoly = Form("BPoly_pt_%d_Centrality_%d_EtaGap_%d_%s",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str());
	h_mCountEff[KEY_BPoly]->Fit(f_BPoly,"N");
	g_mRho00[KEY_Rho00_BPoly]->SetPoint(i_pt,pt_mean,f_BPoly->GetParameter(0));
	g_mRho00[KEY_Rho00_BPoly]->SetPointError(i_pt,0.0,0.0,f_BPoly->GetParError(0),f_BPoly->GetParError(0));
	g_mRho00[KEY_Rho00_BPoly]->SetName(KEY_Rho00_BPoly.c_str());
      }
    }
  }

  TH1F *h_play = new TH1F("h_play","h_play",100,0.0,10.0);
  for(int i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetRangeUser(0.0,5.0);
  h_play->GetXaxis()->SetNdivisions(505,'N');
  h_play->GetXaxis()->SetLabelSize(0.03);
  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetXaxis()->SetTitleSize(0.05);
  h_play->GetXaxis()->SetTitleOffset(1.2);
  h_play->GetXaxis()->CenterTitle();

  h_play->GetYaxis()->SetRangeUser(0.2,0.5);
  h_play->GetYaxis()->SetNdivisions(505,'N');
  h_play->GetYaxis()->SetTitle("#rho_{00}");
  h_play->GetYaxis()->SetTitleSize(0.05);
  h_play->GetYaxis()->SetLabelSize(0.03);
  h_play->GetYaxis()->CenterTitle();

  TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,1000,500);
  c_rho->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_rho->cd(i_pad+1);
    c_rho->cd(i_pad+1)->SetLeftMargin(0.15);
    c_rho->cd(i_pad+1)->SetBottomMargin(0.15);
    c_rho->cd(i_pad+1)->SetGrid(0,0);
    c_rho->cd(i_pad+1)->SetTicks(1,1);
  }
  c_rho->cd(1);
  h_play->DrawCopy("pE");
  PlotLine(0,vmsa::ptMax,1.0/3.0,1.0/3.0,1,2,2);
  plotTopLegend((char*)"bin counting",0.6,0.47,0.04,1,0.0,42,0);

  string KEY_Rho00_CPointQA = Form("Rho00_CPoint_Centrality_%d_EtaGap_%d_%s",vmsa::Cent_QA,vmsa::Eta_QA,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00[KEY_Rho00_CPointQA] ,24,Color_Counts,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.45,0.0,0.0,0.0,0.0,24,Color_Counts,1.3);
  plotTopLegend((char*)"#rho_{00} with P2P efficiency correction",0.7,0.447,0.04,1,0.0,42,0);

  string KEY_Rho00_CPolyQA = Form("Rho00_CPoly_Centrality_%d_EtaGap_%d_%s",vmsa::Cent_QA,vmsa::Eta_QA,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00[KEY_Rho00_CPolyQA] ,20,Color_Counts,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.43,0.0,0.0,0.0,0.0,20,Color_Counts,1.3);
  plotTopLegend((char*)"#rho_{00} with pol0 efficiency correction",0.7,0.427,0.04,1,0.0,42,0);



  c_rho->cd(2);
  h_play->DrawCopy("pE");
  PlotLine(0,vmsa::ptMax,1.0/3.0,1.0/3.0,1,2,2);
  plotTopLegend((char*)"breit wigner fitting",0.6,0.47,0.04,1,0.0,42,0);

  string KEY_Rho00_BPointQA = Form("Rho00_BPoint_Centrality_%d_EtaGap_%d_%s",vmsa::Cent_QA,vmsa::Eta_QA,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00[KEY_Rho00_BPointQA] ,24,Color_BW,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.45,0.0,0.0,0.0,0.0,24,Color_BW,1.3);
  plotTopLegend((char*)"#rho_{00} with P2P efficiency correction",0.7,0.447,0.04,1,0.0,42,0);


  string KEY_Rho00_BPolyQA = Form("Rho00_BPoly_Centrality_%d_EtaGap_%d_%s",vmsa::Cent_QA,vmsa::Eta_QA,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00[KEY_Rho00_BPolyQA] ,20,Color_BW,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.43,0.0,0.0,0.0,0.0,20,Color_BW,1.3);
  plotTopLegend((char*)"#rho_{00} with pol0 efficiency correction",0.7,0.427,0.04,1,0.0,42,0);

  string OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/EffRhoPt.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_play->Write();
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    {
      string KEY_Rho00_CPoint = Form("Rho00_CPoint_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00[KEY_Rho00_CPoint]->SetMarkerSize(1.1);
      g_mRho00[KEY_Rho00_CPoint]->SetMarkerStyle(24);
      g_mRho00[KEY_Rho00_CPoint]->SetMarkerColor(Color_Counts);
      g_mRho00[KEY_Rho00_CPoint]->Write();
      string KEY_Rho00_CPoly = Form("Rho00_CPoly_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00[KEY_Rho00_CPoly]->SetMarkerSize(1.1);
      g_mRho00[KEY_Rho00_CPoly]->SetMarkerStyle(20);
      g_mRho00[KEY_Rho00_CPoly]->SetMarkerColor(Color_Counts);
      g_mRho00[KEY_Rho00_CPoly]->Write();
      string KEY_Rho00_BPoint = Form("Rho00_BPoint_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00[KEY_Rho00_BPoint]->SetMarkerSize(1.1);
      g_mRho00[KEY_Rho00_BPoint]->SetMarkerStyle(24);
      g_mRho00[KEY_Rho00_BPoint]->SetMarkerColor(Color_BW);
      g_mRho00[KEY_Rho00_BPoint]->Write();
      string KEY_Rho00_BPoly = Form("Rho00_BPoly_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00[KEY_Rho00_BPoly]->SetMarkerSize(1.1);
      g_mRho00[KEY_Rho00_BPoly]->SetMarkerStyle(20);
      g_mRho00[KEY_Rho00_BPoly]->SetMarkerColor(Color_BW);
      g_mRho00[KEY_Rho00_BPoly]->Write();
    }
  }
  File_OutPut->Close();
}

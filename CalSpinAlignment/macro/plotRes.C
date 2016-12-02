#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "draw.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

static const TString Energy[2] = {"200GeV","39GeV"};
static const TString PID[2] = {"Phi","KStar"};

static const Int_t Centrality_total = 10;
static const Int_t Centrality_start = 0;
static const Int_t Centrality_stop  = 10;
static const TString mCentrality[10] = {"2060","0005","0510","1020","2030","3040","4050","5060","6070","7080"};
static const TString mCentLabel[10] = {"20-60%","0-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%"};

static const TString mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
static const TString mMode_AMPT[2] = {"Default","StringMelting"};
static const TString mScreenMass_AMPT[3] = {"1mb","3mb","6mb"};

void PlotRes(Int_t mEnergy = 0, Int_t mPID = 0)
{
  const Int_t mColor[Centrality_total] = {kGray+2,kBlack,kRed,kOrange,kMagenta,kAzure,kCyan,kGreen,kViolet,kBlue};
  const Int_t mStyle[Centrality_total] = {20,21,22,23,24,25,26,27,28,29};

  TString input_AMPT = Form("./Data/AuAu%s/%s/Resolution_200GeV.root",Energy[mEnergy].Data(),PID[mPID].Data());
  TFile *File_InPutAMPT = TFile::Open(input_AMPT.Data());
  TProfile *p_mRes_before = (TProfile*)File_InPutAMPT->Get("p_mRes_before");
  TProfile *p_mRes_after = (TProfile*)File_InPutAMPT->Get("p_mRes_after");

  TString input_Data = Form("./Data/AuAu%s/%s/file_200GeV_Resolution.root",Energy[mEnergy].Data(),PID[mPID].Data());
  TFile *File_InPutData = TFile::Open(input_Data.Data());
  TProfile *p_mRes_data = (TProfile*)File_InPutData->Get("Res2_EtaGap_0_EP");

  Float_t Centrality_start[9] = {0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05, 0.0};
  Float_t Centrality_stop[9]  = {0.8,0.7,0.6,0.5,0.4,0.3,0.2, 0.1,0.05};

  TGraphAsymmErrors *g_mRes_before = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_mRes_after = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_mRes_data = new TGraphAsymmErrors();
  for(Int_t i_cent = 1; i_cent < 10; i_cent++)
  {
    if(p_mRes_data->GetBinContent(i_cent) > 0.0)
    {
      Float_t res2 = TMath::Sqrt(p_mRes_data->GetBinContent(i_cent));
      Float_t err_res2 = p_mRes_data->GetBinError(i_cent)/(2*TMath::Sqrt(p_mRes_data->GetBinContent(i_cent)));
      g_mRes_data->SetPoint(i_cent-1,50.0*(Centrality_start[i_cent-1]+Centrality_stop[i_cent-1]),res2*100.0);
      g_mRes_data->SetPointError(i_cent-1,0.0,0.0,err_res2*100.0,err_res2*100.0);
    }
    if(p_mRes_before->GetBinContent(i_cent) > 0.0)
    {
      Float_t res2 = TMath::Sqrt(p_mRes_before->GetBinContent(i_cent));
      Float_t err_res2 = p_mRes_before->GetBinError(i_cent)/(2*TMath::Sqrt(p_mRes_before->GetBinContent(i_cent)));
      g_mRes_before->SetPoint(i_cent-1,50.0*(Centrality_start[i_cent-1]+Centrality_stop[i_cent-1]),res2*100.0);
      g_mRes_before->SetPointError(i_cent-1,0.0,0.0,err_res2*100.0,err_res2*100.0);
    }
    if(p_mRes_after->GetBinContent(i_cent) > 0.0)
    {
      Float_t res2 = TMath::Sqrt(p_mRes_after->GetBinContent(i_cent));
      Float_t err_res2 = p_mRes_after->GetBinError(i_cent)/(2*TMath::Sqrt(p_mRes_after->GetBinContent(i_cent)));
      g_mRes_after->SetPoint(i_cent-1,50.0*(Centrality_start[i_cent-1]+Centrality_stop[i_cent-1]),res2*100.0);
      g_mRes_after->SetPointError(i_cent-1,0.0,0.0,err_res2*100.0,err_res2*100.0);
    }
  }

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->cd();
  c_play->cd()->SetLeftMargin(0.15);
  c_play->cd()->SetBottomMargin(0.15);
  c_play->cd()->SetTicks(1,1);
  c_play->cd()->SetGrid(0,0);
  TH1F *h_play = new TH1F("h_play","h_play",100,0,100);
  for(Int_t i_bin = 1; i_bin < 1001; i_bin++)
  {
    h_play->SetBinContent(i_bin,-10.0);
    h_play->SetBinError(i_bin,0.1);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetNdivisions(505,'X');
  h_play->GetYaxis()->SetNdivisions(505,'Y');
  h_play->GetXaxis()->SetTitle("Centrality(%)");
  h_play->GetYaxis()->SetTitle("Resolution");
  h_play->GetXaxis()->SetTitleSize(0.04);
  h_play->GetYaxis()->SetTitleSize(0.04);
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->CenterTitle();
  h_play->GetXaxis()->SetRangeUser(0,85);
  h_play->GetYaxis()->SetRangeUser(0,100);
  h_play->SetLineColor(1);
  h_play->SetMarkerStyle(24);
  h_play->SetMarkerColor(1);
  h_play->SetMarkerSize(0.8);
  h_play->Draw("pE");
  Draw_TGAE_new_Symbol(g_mRes_before,mStyle[0],mColor[0],1.2);
  Draw_TGAE_Point_new_Symbol(42,90,0.0,0.0,0.0,0.0,mStyle[0],mColor[0],1.0);
  plotTopLegend("AMPT w/o Efficiency Correction",43.5,89.0,0.025,1,0.0,42,0);

  // Draw_TGAE_new_Symbol(g_mRes_after,24,4,1.2);
  // Draw_TGAE_Point_new_Symbol(42,86,0.0,0.0,0.0,0.0,24,4,1.0);
  // plotTopLegend("AMPT with Efficiency Correction",43.5,85.0,0.025,1,0.0,42,0);

  Draw_TGAE_new_Symbol(g_mRes_data,29,2,1.2);
  Draw_TGAE_Point_new_Symbol(42,82,0.0,0.0,0.0,0.0,29,2,1.0);
  plotTopLegend("STAR data",43.5,81.0,0.025,1,0.0,42,0);

  plotTopLegend("Au+Au 200 GeV",20,20,0.03,1,0.0,42,0);
  c_play->SaveAs("./figures/Res_Com.eps");

  TString input_SE = Form("./Data/AuAu%s/Phi/SpinAlignment_%s_SE.root",Energy[mEnergy].Data(),Energy[mEnergy].Data());
  TFile *File_InPutSE = TFile::Open(input_SE.Data());
  TH1F *h_RP = (TH1F*)File_InPutSE->Get("h_SpinAlignment_Phi_2060_CosThetaStar_3_14_RP");
  TH1F *h_TPC = (TH1F*)File_InPutSE->Get("h_SpinAlignment_Phi_2060_CosThetaStar_3_14_TPC");

  TCanvas *c_Yield = new TCanvas("c_Yield","c_Yield",10,10,800,800);
  c_Yield->cd();
  c_Yield->cd()->SetLeftMargin(0.15);
  c_Yield->cd()->SetBottomMargin(0.15);
  c_Yield->cd()->SetTicks(1,1);
  c_Yield->cd()->SetGrid(0,0);
  h_RP->SetTitle("");
  h_RP->SetStats(0);
  h_RP->GetXaxis()->SetNdivisions(505,'X');
  h_RP->GetYaxis()->SetNdivisions(505,'Y');
  h_RP->GetXaxis()->SetTitle("M_{K^{+},K^{-}} (GeV/c^{2})");
  h_RP->GetYaxis()->SetTitle("Yields");
  h_RP->GetXaxis()->SetTitleSize(0.04);
  h_RP->GetYaxis()->SetTitleSize(0.04);
  h_RP->GetYaxis()->SetLabelSize(0.02);
  h_RP->GetXaxis()->CenterTitle();
  h_RP->GetYaxis()->CenterTitle();
  // h_RP->GetXaxis()->SetRangeUser(0,85);
  // h_RP->GetYaxis()->SetRangeUser(0,100);
  h_RP->SetLineColor(1);
  h_RP->SetMarkerStyle(24);
  h_RP->SetMarkerColor(1);
  h_RP->SetMarkerSize(0.8);
  h_RP->Draw("PE");

  h_TPC->SetLineColor(1);
  h_TPC->SetMarkerStyle(24);
  h_TPC->SetMarkerColor(2);
  h_TPC->SetMarkerSize(0.8);
  h_TPC->Draw("PE same");

  TLegend *leg = new TLegend(0.2,0.7,0.5,0.85);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_RP,"w/o TPC Efficiency","p");
  leg->AddEntry(h_TPC,"with TPC Efficiency","p");
  leg->Draw("same");

  c_Yield->SaveAs("./figures/Yield_Eff.eps");
}

#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
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

TGraphAsymmErrors* calResCorr(TGraphAsymmErrors *g_rho, TF1 *f_res, std::string GraName)
{
  TGraphAsymmErrors *g_ResCorr = new TGraphAsymmErrors();
  float resCorr = f_res->GetParameter(1);
  float err_resCorr = f_res->GetParError(1);
  for(int i_point = 0; i_point < g_rho->GetN(); ++i_point)
  {
    double pt, rho;
    g_rho->GetPoint(i_point,pt,rho);
    double err_rho = g_rho->GetErrorYhigh(i_point);
    double err = ErrDiv(rho,resCorr,err_rho,err_resCorr);
    double rhoCorr = f_res->GetX(rho);
    g_ResCorr->SetPoint(i_point,pt,rhoCorr);
    g_ResCorr->SetPointError(i_point,0.0,0.0,err,err);
  }
  g_ResCorr->SetName(GraName.c_str());

  return g_ResCorr;
}

void appResCorr(int energy = 6, int pid = 0)
{
  gStyle->SetOptDate(0);
  string InPutRes  = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/MonteCarlo/McResCorr/Mc%sResCorrFactor.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_Res = TFile::Open(InPutRes.c_str());
  TGraMap g_Res;
  TF1Map f_Res;
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; ++i_cent)
  {
    string KEY_Res = Form("ResEP_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // use Res EP only
    g_Res[KEY_Res] = (TGraphAsymmErrors*)File_Res->Get(KEY_Res.c_str());
    f_Res[KEY_Res] = new TF1(KEY_Res.c_str(),PolyRes,0.0,0.6,2);
    f_Res[KEY_Res]->SetParameter(0,0.01);
    f_Res[KEY_Res]->SetParameter(1,0.8);
    g_Res[KEY_Res]->Fit(f_Res[KEY_Res],"N");
  }

#if _PlotQA_
  TH1F *h_res = (TH1F*)File_Res->Get("h_res");
  TCanvas *c_res = new TCanvas("c_res","c_res",10,10,800,800);
  c_res->cd();
  c_res->cd()->SetLeftMargin(0.15);
  c_res->cd()->SetBottomMargin(0.15);
  c_res->cd()->SetGrid(0,0);
  c_res->cd()->SetTicks(1,1);
  h_res->Draw("pE");
  PlotLine(1.0/3.0,1.0/3.0,0.0,0.6,1,2,2);
  PlotLine(0.0,0.6,1.0/3.0,1.0/3.0,1,2,2);
  string KEY_ResQA = Form("ResEP_Centrality_%d_%s",vmsa::Cent_QA,vmsa::mPID[pid].c_str()); // use Res EP only
  g_Res[KEY_ResQA]->Draw("pE same");
  f_Res[KEY_ResQA]->SetLineStyle(2);
  f_Res[KEY_ResQA]->SetLineColor(2);
  f_Res[KEY_ResQA]->SetLineWidth(4);
  f_Res[KEY_ResQA]->Draw("l same");
  string leg_p0   = Form("p_{0}^{Sergei} = %2.4f #pm %2.4f",f_Res[KEY_ResQA]->GetParameter(0),f_Res[KEY_ResQA]->GetParError(0));
  string leg_p1   = Form("p_{1}^{Sergei} = %2.4f #pm %2.4f",f_Res[KEY_ResQA]->GetParameter(1),f_Res[KEY_ResQA]->GetParError(1));
  string leg_chi2 = Form("#chi^{2}/NDF = %2.2f/%d",f_Res[KEY_ResQA]->GetChisquare(),f_Res[KEY_ResQA]->GetNDF());
  plotTopLegend((char*)leg_p0.c_str(),0.58,0.40,0.03,kRed,0.0,42,1);
  plotTopLegend((char*)leg_p1.c_str(),0.58,0.34,0.03,kRed,0.0,42,1);
  plotTopLegend((char*)leg_chi2.c_str(),0.58,0.28,0.03,kRed,0.0,42,1);
#endif

  string InPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/EffRhoPt.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_InPut = TFile::Open(InPutFile.c_str());
  TGraMap g_mRho00Eff, g_mRho00;
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    {
      string KEY_Res = Form("ResEP_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // use Res EP only
      string KEY_Rho00_CPoint = Form("Rho00_CPoint_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00Eff[KEY_Rho00_CPoint] = (TGraphAsymmErrors*)File_InPut->Get(KEY_Rho00_CPoint.c_str());
      g_mRho00[KEY_Rho00_CPoint] = calResCorr(g_mRho00Eff[KEY_Rho00_CPoint],f_Res[KEY_Res],KEY_Rho00_CPoint.c_str());

      string KEY_Rho00_CPoly = Form("Rho00_CPoly_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00Eff[KEY_Rho00_CPoly] = (TGraphAsymmErrors*)File_InPut->Get(KEY_Rho00_CPoly.c_str());
      g_mRho00[KEY_Rho00_CPoly] = calResCorr(g_mRho00Eff[KEY_Rho00_CPoly],f_Res[KEY_Res],KEY_Rho00_CPoly.c_str());

      string KEY_Rho00_BPoint = Form("Rho00_BPoint_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00Eff[KEY_Rho00_BPoint] = (TGraphAsymmErrors*)File_InPut->Get(KEY_Rho00_BPoint.c_str());
      g_mRho00[KEY_Rho00_BPoint] = calResCorr(g_mRho00Eff[KEY_Rho00_BPoint],f_Res[KEY_Res],KEY_Rho00_BPoint.c_str());

      string KEY_Rho00_BPoly = Form("Rho00_BPoly_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00Eff[KEY_Rho00_BPoly] = (TGraphAsymmErrors*)File_InPut->Get(KEY_Rho00_BPoly.c_str());
      g_mRho00[KEY_Rho00_BPoly] = calResCorr(g_mRho00Eff[KEY_Rho00_BPoly],f_Res[KEY_Res],KEY_Rho00_BPoly.c_str());
    }
  }

#if _PlotQA_
  TH1F *h_play = (TH1F*)File_InPut->Get("h_play");
  TCanvas *c_resCorr =  new TCanvas("c_resCorr","c_resCorr",10,10,1000,500);
  c_resCorr->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_resCorr->cd(i_pad+1);
    c_resCorr->cd(i_pad+1)->SetLeftMargin(0.15);
    c_resCorr->cd(i_pad+1)->SetBottomMargin(0.15);
    c_resCorr->cd(i_pad+1)->SetGrid(0,0);
    c_resCorr->cd(i_pad+1)->SetTicks(1,1);
    h_play->Draw("PE");
    PlotLine(0.0,vmsa::ptMax,1.0/3.0,1.0/3.0,1,2,2);
  }

  c_resCorr->cd(1);
  plotTopLegend((char*)"bin counting",0.7,0.47,0.04,1,0.0,42,0);
  string KEY_Rho00_CPointQA = Form("Rho00_CPoint_Centrality_%d_EtaGap_%d_%s",vmsa::Cent_QA,vmsa::Eta_QA,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00Eff[KEY_Rho00_CPointQA],24,Color_Counts,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.45,0.0,0.0,0.0,0.0,24,Color_Counts,1.3);
  plotTopLegend((char*)"#rho_{00} with P2P efficiency correction",0.7,0.447,0.04,1,0.0,42,0);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00[KEY_Rho00_CPointQA],29,kRed,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.43,0.0,0.0,0.0,0.0,29,kRed,1.3);
  plotTopLegend((char*)"#rho_{00} with resolution correction",0.7,0.427,0.04,1,0.0,42,0);

  // c_resCorr->cd(2);
  // string KEY_Rho00_CPolyQA = Form("Rho00_CPoly_Centrality_%d_EtaGap_%d_%s",vmsa::Cent_QA,vmsa::Eta_QA,vmsa::mPID[pid].c_str());
  // Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00Eff[KEY_Rho00_CPolyQA],20,Color_Counts,1.1);
  // Draw_TGAE_Point_new_Symbol(0.5,0.45,0.0,0.0,0.0,0.0,20,Color_Counts,1.3);
  // plotTopLegend((char*)"#rho_{00} with pol0 efficiency correction",0.7,0.447,0.04,1,0.0,42,0);
  // Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00[KEY_Rho00_CPolyQA],29,kRed,1.1);
  // Draw_TGAE_Point_new_Symbol(0.5,0.43,0.0,0.0,0.0,0.0,29,kRed,1.3);
  // plotTopLegend((char*)"#rho_{00} with resolution correction",0.7,0.427,0.04,1,0.0,42,0);

  c_resCorr->cd(2);
  plotTopLegend((char*)"breit wigner fitting",0.7,0.47,0.04,1,0.0,42,0);
  string KEY_Rho00_BPointQA = Form("Rho00_BPoint_Centrality_%d_EtaGap_%d_%s",vmsa::Cent_QA,vmsa::Eta_QA,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00Eff[KEY_Rho00_BPointQA],24,Color_BW,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.45,0.0,0.0,0.0,0.0,24,Color_BW,1.3);
  plotTopLegend((char*)"#rho_{00} with P2P efficiency correction",0.7,0.447,0.04,1,0.0,42,0);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00[KEY_Rho00_BPointQA],29,kRed,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.43,0.0,0.0,0.0,0.0,29,kRed,1.3);
  plotTopLegend((char*)"#rho_{00} with resolution correction",0.7,0.427,0.04,1,0.0,42,0);

  // c_resCorr->cd(4);
  // string KEY_Rho00_BPolyQA = Form("Rho00_BPoly_Centrality_%d_EtaGap_%d_%s",vmsa::Cent_QA,vmsa::Eta_QA,vmsa::mPID[pid].c_str());
  // Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00Eff[KEY_Rho00_BPolyQA],20,Color_BW,1.1);
  // Draw_TGAE_Point_new_Symbol(0.5,0.45,0.0,0.0,0.0,0.0,20,Color_BW,1.3);
  // plotTopLegend((char*)"#rho_{00} with pol0 efficiency correction",0.7,0.447,0.04,1,0.0,42,0);
  // Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00[KEY_Rho00_BPolyQA],29,kRed,1.1);
  // Draw_TGAE_Point_new_Symbol(0.5,0.43,0.0,0.0,0.0,0.0,29,kRed,1.3);
  // plotTopLegend((char*)"#rho_{00} with resolution correction",0.7,0.427,0.04,1,0.0,42,0);

  TCanvas *c_resPlot = new TCanvas("c_resPlot","c_resPlot",10,10,800,800);
  c_resPlot->cd();
  c_resPlot->cd()->SetLeftMargin(0.15);
  c_resPlot->cd()->SetBottomMargin(0.15);
  c_resPlot->cd()->SetGrid(0,0);
  c_resPlot->cd()->SetTicks(1,1);
  h_play->Draw("PE");
  PlotLine(0.0,vmsa::ptMax,1.0/3.0,1.0/3.0,1,2,2);
  g_mRho00[KEY_Rho00_CPointQA]->SetMarkerStyle(30);
  g_mRho00[KEY_Rho00_CPointQA]->SetMarkerColor(2);
  g_mRho00[KEY_Rho00_CPointQA]->SetMarkerSize(1.5);
  g_mRho00[KEY_Rho00_CPointQA]->SetLineColor(2);
  g_mRho00[KEY_Rho00_CPointQA]->Draw("pE same");
  plotTopLegend((char*)"AuAu 200 GeV, 20%-60%",0.7,0.45,0.04,1,0.0,42,0);
  Draw_TGAE_Point_new_Symbol(0.4,0.43,0.0,0.0,0.0,0.0,30,2,1.5);
  plotTopLegend((char*)"#phi-meson #rho_{00} with efficiency and resolution correction",0.6,0.427,0.03,1,0.0,42,0);
#endif

  string OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/ResRhoPt.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    {
      string KEY_Rho00_CPoint = Form("Rho00_CPoint_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00[KEY_Rho00_CPoint]->Write();

      string KEY_Rho00_CPoly = Form("Rho00_CPoly_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00[KEY_Rho00_CPoly]->Write();

      string KEY_Rho00_BPoint = Form("Rho00_BPoint_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00[KEY_Rho00_BPoint]->Write();

      string KEY_Rho00_BPoly = Form("Rho00_BPoly_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00[KEY_Rho00_BPoly]->Write();
    }
  }
  File_OutPut->Close();

  c_resCorr->SaveAs("../figures/resRhoPt.eps");
}

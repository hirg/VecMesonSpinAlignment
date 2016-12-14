#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "../../Utility/type.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/draw.h"
#include "../../Utility/functions.h"

TGraphAsymmErrors* calRatio(TGraphAsymmErrors *g_rho, float norm, std::string GraName)
{
  TGraphAsymmErrors *g_ResCorr = new TGraphAsymmErrors();
  for(int i_point = 0; i_point < g_rho->GetN(); ++i_point)
  {
    double pt, rho;
    g_rho->GetPoint(i_point,pt,rho);
    double err_rho = g_rho->GetErrorYhigh(i_point);
    double err = ErrDiv(rho,norm,err_rho,0.0);
    double ratio = rho/norm;
    g_ResCorr->SetPoint(i_point,pt,ratio);
    g_ResCorr->SetPointError(i_point,0.0,0.0,err,err);
  }
  g_ResCorr->SetName(GraName.c_str());

  return g_ResCorr;
}

void plotRhoRatio(int energy = 6, int pid = 0)
{
  gStyle->SetOptDate(0);
  string InPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/ResRhoPt.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_InPut = TFile::Open(InPutFile.c_str());
  TGraMap g_mRho00;
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    {
      string KEY_Rho00_CPoint = Form("Rho00_CPoint_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str());
      g_mRho00[KEY_Rho00_CPoint] = (TGraphAsymmErrors*)File_InPut->Get(KEY_Rho00_CPoint.c_str());
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------
  TCanvas* c_v2_Sim = new TCanvas("c_v2_Sim","c_v2_Sim",600,20,620,700);
  c_v2_Sim->SetFillColor(10);
  c_v2_Sim->SetTopMargin(0.03);
  c_v2_Sim->SetBottomMargin(0.2);
  c_v2_Sim->SetRightMargin(0.25);
  c_v2_Sim->SetLeftMargin(0.23);
  c_v2_Sim->Divide(1,2,0.0,0.0,10);
  //c_v2_Sim->Divide(1,2);
  TH1F* h_test_sim_B2 = new TH1F("h_test_sim_B2","h_test_sim_B2",1000,-0.5,6);
  TH1F* h_test_sim_B3 = new TH1F("h_test_sim_B3","h_test_sim_B3",1000,-0.5,6);

  for(Int_t i = 0; i < h_test_sim_B2->GetNbinsX(); i++)
  {
    h_test_sim_B2->SetBinContent(i,-10);
    h_test_sim_B3->SetBinContent(i,-10);
  }

  for(Int_t i = 0; i < 2; i++)
  {
    c_v2_Sim->cd(i+1);
    c_v2_Sim->cd(i+1);
    c_v2_Sim->cd(i+1)->SetTicks(1,1);
    c_v2_Sim->cd(i+1)->SetGrid(0,0);
    c_v2_Sim->cd(i+1)->SetFillColor(10);
    if(i == 0 || i == 1) c_v2_Sim->cd(i+1)->SetLeftMargin(0.25);
    Float_t x1_array[1] = {0.0};
    Float_t x2_array[1] = {0.92};
    if(i < 1)
    {
      c_v2_Sim->cd(i+1)->SetBottomMargin(0.0);
      c_v2_Sim->cd(i+1)->SetTopMargin(0.1);
      c_v2_Sim->cd(i+1)->SetPad(x1_array[i],0.4,x2_array[i],1.0); // x1, y1, x2, y2

      h_test_sim_B2->SetStats(0);
      h_test_sim_B2->SetTitle("");
      h_test_sim_B2->GetXaxis()->SetTitleOffset(1.2);
      h_test_sim_B2->GetYaxis()->SetTitleOffset(0.93);
      h_test_sim_B2->GetYaxis()->SetLabelOffset(0.009);
      h_test_sim_B2->GetXaxis()->SetLabelSize(0.04);
      h_test_sim_B2->GetYaxis()->SetLabelSize(0.04);
      h_test_sim_B2->GetXaxis()->SetTitleSize(0.11);
      h_test_sim_B2->GetYaxis()->SetTitleSize(0.11);
      h_test_sim_B2->GetXaxis()->SetNdivisions(505,'N');
      h_test_sim_B2->GetYaxis()->SetNdivisions(505,'N');
      h_test_sim_B2->GetXaxis()->CenterTitle();
      h_test_sim_B2->GetYaxis()->CenterTitle();
      h_test_sim_B2->GetYaxis()->SetTitle("#rho_{00}");
      h_test_sim_B2->GetXaxis()->SetTitle("");
      h_test_sim_B2->GetXaxis()->SetRangeUser(0.0,vmsa::ptMax);
      h_test_sim_B2->GetYaxis()->SetRangeUser(0.27,0.5);
      h_test_sim_B2->DrawCopy("h");

      PlotLine(0.0,vmsa::ptMax,1./3.,1./3.,1,2,2); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle
      string KEY_Rho00_CPointQA = Form("Rho00_CPoint_Centrality_%d_EtaGap_%d_%s",vmsa::Cent_QA,vmsa::Eta_QA,vmsa::mPID[pid].c_str());
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00[KEY_Rho00_CPointQA],30,kRed,1.3);
      plotTopLegend((char*)"AuAu 200 GeV, 20%-60%",0.7,0.45,0.04,1,0.0,42,0);
      Draw_TGAE_Point_new_Symbol(0.4,0.43,0.0,0.0,0.0,0.0,30,2,1.5);
      plotTopLegend((char*)"#phi-meson #rho_{00} with efficiency and resolution correction",0.6,0.427,0.04,1,0.0,42,0);
    }
    if(i >= 1)
    {
      Float_t pad_factor = 1.5;
      c_v2_Sim->cd(i+1)->SetTopMargin(0.0);
      c_v2_Sim->cd(i+1)->SetBottomMargin(0.5);
      c_v2_Sim->cd(i+1)->SetPad(x1_array[i-1],0,x2_array[i-1],0.4); // x1, y1, x2, y2

      h_test_sim_B3->SetStats(0);
      h_test_sim_B3->SetTitle("");
      h_test_sim_B3->GetXaxis()->SetTitleOffset(1.2);
      h_test_sim_B3->GetYaxis()->SetTitleOffset(0.65);
      h_test_sim_B3->GetXaxis()->SetLabelOffset(0.02);
      h_test_sim_B3->GetYaxis()->SetLabelOffset(0.009);
      h_test_sim_B3->GetXaxis()->SetLabelSize(0.04*pad_factor);
      h_test_sim_B3->GetYaxis()->SetLabelSize(0.04*pad_factor);
      h_test_sim_B3->GetXaxis()->SetTitleSize(0.11*pad_factor);
      h_test_sim_B3->GetYaxis()->SetTitleSize(0.11*pad_factor);
      h_test_sim_B3->GetXaxis()->SetNdivisions(505,'N');
      h_test_sim_B3->GetYaxis()->SetNdivisions(303,'N');
      h_test_sim_B3->GetXaxis()->CenterTitle();
      h_test_sim_B3->GetYaxis()->CenterTitle();
      h_test_sim_B3->GetYaxis()->SetTitle("ratio");
      h_test_sim_B3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h_test_sim_B3->GetXaxis()->SetRangeUser(0.0,vmsa::ptMax);
      h_test_sim_B3->GetYaxis()->SetRangeUser(0.75,1.25);
      h_test_sim_B3->DrawCopy("h");

      PlotLine(0.0,vmsa::ptMax,1.,1.,1,2,2); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle
      string KEY_Rho00_CPointQA = Form("Rho00_CPoint_Centrality_%d_EtaGap_%d_%s",vmsa::Cent_QA,vmsa::Eta_QA,vmsa::mPID[pid].c_str());
      TGraphAsymmErrors *g_ratio = calRatio(g_mRho00[KEY_Rho00_CPointQA],1./3.,"g_ratio");
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ratio,30,kRed,1.3);
    }
  }
  //------------------------------------------------------------------------------------------------------------------------------------

  c_v2_Sim->SaveAs("../../figures/rhoRatioPt.eps");
}

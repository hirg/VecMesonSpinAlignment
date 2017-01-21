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

void appResCorrSys(int energy = 6, int pid = 0)
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

  // string InPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/EffRhoPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string InPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/RawRhoPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_InPut = TFile::Open(InPutFile.c_str());
  TGraMap g_mRhoEff, g_mRho;
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; ++i_eta) // EtaGap loop
    {
      for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
      {
	for(int i_func = vmsa::Func_start; i_func < vmsa::Func_stop; ++i_func)
	{
	  for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	  {
	    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	    {
	      string KEY_rho = Form("rhoRaw_Centrality_%d_EtaGap_%d_2nd_%s_Norm_%d_Func_%d_Sigma_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str(),i_norm,i_func,i_sigma,vmsa::mInteMethod[i_method].c_str());
	      g_mRhoEff[KEY_rho] = (TGraphAsymmErrors*)File_InPut->Get(KEY_rho.c_str());
	      string KEY_Res = Form("ResEP_Centrality_%d_%s",i_cent,vmsa::mPID[pid].c_str()); // use Res EP only
	      g_mRho[KEY_rho] = calResCorr(g_mRhoEff[KEY_rho],f_Res[KEY_Res],KEY_rho.c_str());
	    }
	  }
	}
      }
    }
  }

  TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,800,800);
  c_rho->cd();
  c_rho->cd()->SetLeftMargin(0.15);
  c_rho->cd()->SetBottomMargin(0.15);
  c_rho->cd()->SetTicks(1,1);
  c_rho->cd()->SetGrid(0,0);
  TH1F *h_frame = new TH1F("h_frame","h_frame",100,-0.05,9.95);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,1.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,5.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.03);
  h_frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_frame->GetXaxis()->SetTitleSize(0.05);
  h_frame->GetXaxis()->SetTitleOffset(1.2);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.2,0.5);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.05);
  h_frame->GetYaxis()->SetLabelSize(0.03);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; ++i_eta) // EtaGap loop
    {
      for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
      {
	for(int i_func = vmsa::Func_start; i_func < vmsa::Func_stop; ++i_func)
	{
	  for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	  {
	    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	    {
	      string KEY_rho = Form("rhoRaw_Centrality_%d_EtaGap_%d_2nd_%s_Norm_%d_Func_%d_Sigma_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str(),i_norm,i_func,i_sigma,vmsa::mInteMethod[i_method].c_str());
	      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho[KEY_rho],24,i_sigma+10*i_method+1,1.1);
	    }
	  }
	}
      }
    }
  }


  string OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/ResRhoPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  h_frame->Write();
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; ++i_eta) // EtaGap loop
    {
      for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
      {
	for(int i_func = vmsa::Func_start; i_func < vmsa::Func_stop; ++i_func)
	{
	  for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	  {
	    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	    {
	      string KEY_rho = Form("rhoRaw_Centrality_%d_EtaGap_%d_2nd_%s_Norm_%d_Func_%d_Sigma_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str(),i_norm,i_func,i_sigma,vmsa::mInteMethod[i_method].c_str());
	      g_mRho[KEY_rho]->SetName(KEY_rho.c_str());
	      g_mRho[KEY_rho]->Write();
	    }
	  }
	}
      }
    }
  }
  File_OutPut->Close();
}

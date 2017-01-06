#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TProfile2D.h"
#include "../Utility/functions.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

void calSigBgRatio(int energy = 3, int pid = 0)
{
  TGaxis::SetMaxDigits(4);

  string InPutFile_SE = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/Yields/merged_file/Yields_SE_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  
  string InPutFile_ME = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/Yields/merged_file/Yields_ME_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_SE = TFile::Open(InPutFile_SE.c_str());
  TFile *File_ME = TFile::Open(InPutFile_ME.c_str());

  // read in histogram for same event and mixed event
  // merged cos(theta*) bins
  TH1FMap h_mMass_SE, h_mMass_ME;
  TH1FMap h_mSumMass_SE, h_mSumMass_ME;
  for(int i_pt = vmsa::pt_start; i_pt < vmsa::pt_stop; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
    {
      for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
      {
	string KEY_sumSE = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_sumSE",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str());
	string KEY_sumME = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_sumME",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str());
	for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	{
	  string KEY_SE = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_SE",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str());
	  h_mMass_SE[KEY_SE] = (TH1F*)File_SE->Get(KEY_SE.c_str())->Clone(); 

	  string KEY_ME = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_ME",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str());
	  h_mMass_ME[KEY_ME] = (TH1F*)File_ME->Get(KEY_ME.c_str())->Clone(); 

	  if(i_theta == vmsa::CTS_start)
	  {
	    h_mSumMass_SE[KEY_sumSE] = (TH1F*)h_mMass_SE[KEY_SE]->Clone();
	    h_mSumMass_ME[KEY_sumME] = (TH1F*)h_mMass_ME[KEY_ME]->Clone();
	  }
	  else
	  {
	    h_mSumMass_SE[KEY_sumSE]->Add(h_mMass_SE[KEY_SE],1.0);
	    h_mSumMass_ME[KEY_sumME]->Add(h_mMass_ME[KEY_ME],1.0);
	  }
	}
      }
    }
  }

  // pT rebin
  TH1FMap h_mMassSE, h_mMassME, h_mMassSM; // rebinned InvMass distribution, SE-ME

  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    {
      for(int pt_bin = vmsa::pt_rebin_first[energy]; pt_bin < vmsa::pt_rebin_last[energy]; pt_bin++) // pt loop
      {
	string KeySE = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_sumSE",pt_bin,i_cent,i_eta,vmsa::mPID[pid].c_str());
	string KeyME = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_sumME",pt_bin,i_cent,i_eta,vmsa::mPID[pid].c_str());
	for(int i_pt = vmsa::pt_rebin_start[energy][pt_bin]; i_pt <= vmsa::pt_rebin_stop[energy][pt_bin]; i_pt++)
	{
	  string KEY_sumSE = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_sumSE",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str());
	  string KEY_sumME = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_sumME",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str());
	  if(i_pt == vmsa::pt_rebin_start[energy][pt_bin])
	  {
	    h_mMassSE[KeySE] = (TH1F*)h_mSumMass_SE[KEY_sumSE]->Clone();
	    h_mMassME[KeyME] = (TH1F*)h_mSumMass_ME[KEY_sumME]->Clone();
	  }
	  else
	  {
	    h_mMassSE[KeySE]->Add(h_mSumMass_SE[KEY_sumSE],1.0);
	    h_mMassME[KeyME]->Add(h_mSumMass_ME[KEY_sumME],1.0);
	  }
	}
	int Norm_bin_start = h_mMassSE[KeySE]->FindBin(vmsa::Norm_Start[pid]);
	int Norm_bin_stop  = h_mMassSE[KeySE]->FindBin(vmsa::Norm_Stop[pid]);
	float InteSE = h_mMassSE[KeySE]->Integral(Norm_bin_start,Norm_bin_stop);
	float InteME = h_mMassME[KeyME]->Integral(Norm_bin_start,Norm_bin_stop);
	h_mMassME[KeyME]->Scale(InteSE/InteME);

	string KeySM = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_sumSM",pt_bin,i_cent,i_eta,vmsa::mPID[pid].c_str());
	h_mMassSM[KeySM] = (TH1F*)h_mMassSE[KeySE]->Clone();
	h_mMassSM[KeySM]->Add(h_mMassME[KeyME],-1.0);
      }
    }
  }

#if _PlotQA_
  TCanvas *c_SigBg = new TCanvas("c_SigBg","c_SigBg",10,10,900,900);
  c_SigBg->Divide(3,3);
  for(int i_pad = vmsa::pt_rebin_first[energy]; i_pad < vmsa::pt_rebin_last[energy]; ++i_pad)
  {
    c_SigBg->cd(i_pad+1);
    c_SigBg->cd(i_pad+1)->SetLeftMargin(0.15);
    c_SigBg->cd(i_pad+1)->SetBottomMargin(0.15);
    c_SigBg->cd(i_pad+1)->SetTicks(1,1);
    c_SigBg->cd(i_pad+1)->SetGrid(0,0);

    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
    {
      for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
      {
	string KeySE = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_sumSE",i_pad,i_cent,i_eta,vmsa::mPID[pid].c_str());
	h_mMassSE[KeySE]->SetMarkerColor(1);
	h_mMassSE[KeySE]->SetMarkerStyle(24);
	h_mMassSE[KeySE]->SetMarkerSize(0.8);
	h_mMassSE[KeySE]->DrawCopy("pE");

	string KeyME = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_sumME",i_pad,i_cent,i_eta,vmsa::mPID[pid].c_str());
	h_mMassME[KeyME]->SetLineColor(2);
	h_mMassME[KeyME]->SetFillColor(2);
	h_mMassME[KeyME]->SetFillStyle(3002);
	h_mMassME[KeyME]->DrawCopy("h same");

	string KeySM = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_sumSM",i_pad,i_cent,i_eta,vmsa::mPID[pid].c_str());
	h_mMassSM[KeySM]->SetLineColor(4);
	h_mMassSM[KeySM]->SetFillColor(4);
	h_mMassSM[KeySM]->SetFillStyle(3004);
	h_mMassSM[KeySM]->DrawCopy("h same");
      }
    }
  }
#endif

  TGraphAsymmErrors *g_ratio = new TGraphAsymmErrors();
  float inteStart = vmsa::InvMass[pid]+vmsa::nSigVec*vmsa::Width[pid];
  float inteStop  = vmsa::InvMass[pid]-vmsa::nSigVec*vmsa::Width[pid];
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    {
      for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt)
      {
	string KeySM = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_sumSM",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str());
	int inteBinStart = h_mMassSM[KeySM]->FindBin(inteStart);
	int inteBinStop = h_mMassSM[KeySM]->FindBin(inteStop);
	float Sig = h_mMassSM[KeySM]->Integral(inteBinStart,inteBinStop);

	string KeyME = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_sumME",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str());
	float Bg = h_mMassME[KeyME]->Integral(inteBinStart,inteBinStop);
	float raito = Sig/Bg;
	float pt = 0.5*(vmsa::pt_low[energy][i_pt]+vmsa::pt_up[energy][i_pt]);
	float ptWidth = 0.5*(vmsa::pt_up[energy][i_pt]-vmsa::pt_low[energy][i_pt]);
	g_ratio->SetPoint(i_pt,pt,raito);
	g_ratio->SetPointError(i_pt,ptWidth,ptWidth,0.0,0.0);
      }
    }
  }

  TCanvas *c_ratio = new TCanvas("c_ratio","c_ratio",10,10,800,800);
  c_ratio->cd();
  c_ratio->cd()->SetLeftMargin(0.15);
  c_ratio->cd()->SetBottomMargin(0.15);
  c_ratio->cd()->SetTicks(1,1);
  c_ratio->cd()->SetGrid(0,0);
  TH1F *h_frame = new TH1F("h_frame","h_frame",100,0,10);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,-10.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,7.2);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.03);
  h_frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_frame->GetXaxis()->SetTitleSize(0.05);
  h_frame->GetXaxis()->SetTitleOffset(1.2);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(-0.01,0.2);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("Sig/Bg");
  h_frame->GetYaxis()->SetTitleSize(0.05);
  h_frame->GetYaxis()->SetLabelSize(0.03);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ratio,24,2,1.1);
  g_ratio->SetName("g_SigBgRatio");
  g_ratio->SetMarkerStyle(24);
  g_ratio->SetMarkerSize(1.4);
  g_ratio->SetMarkerColor(2);

  string OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/SigBgRatio.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  g_ratio->Write();
  File_OutPut->Close();
}

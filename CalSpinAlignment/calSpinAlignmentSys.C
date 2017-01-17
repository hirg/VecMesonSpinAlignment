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


void calSpinAlignmentSys(int energy = 3, int pid = 0, int year = 0)
{
  string inputfile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  File_InPut->cd();
  TH1FMap h_mMass, h_mMass_InteTheta;
  vecFMap Par_InteTheta;
  // read in histograms
  // integrated over cos(theta*) and do breit wiger fit to extract common fit parameter
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
    {
      for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
      {
	for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	{
	  for(int i_func = vmsa::Func_start; i_func < vmsa::Func_stop; ++i_func)
	  {
	    string KEY_InteTheta = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_Norm_%d_Func_%d",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str(),i_norm,i_func);
	    for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	    {
	      string KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_Func_%d",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str(),i_norm,i_func);
	      h_mMass[KEY] = (TH1F*)File_InPut->Get(KEY.c_str());

	      if(i_theta == vmsa::CTS_start) h_mMass_InteTheta[KEY_InteTheta] = (TH1F*)h_mMass[KEY]->Clone(KEY_InteTheta.c_str());
	      else h_mMass_InteTheta[KEY_InteTheta]->Add(h_mMass[KEY],1.0);
	    }
	    TF1 *f_bw = new TF1("f_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
	    f_bw->SetParameter(0,vmsa::InvMass[pid]);
	    f_bw->SetParLimits(0,vmsa::InvMass[pid]-0.001,vmsa::InvMass[pid]+0.001);
	    f_bw->SetParameter(1,vmsa::Width[pid]);
	    f_bw->SetParameter(2,1.0);
	    float norm = h_mMass_InteTheta[KEY_InteTheta]->GetMaximum()/f_bw->GetMaximum();
	    f_bw->SetParameter(2,norm);
	    f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	    h_mMass_InteTheta[KEY_InteTheta]->Fit(f_bw,"MQNR");
	    Par_InteTheta[KEY_InteTheta].clear();
	    Par_InteTheta[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(0)));
	    Par_InteTheta[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(1)));
	    Par_InteTheta[KEY_InteTheta].push_back(static_cast<float>(f_bw->GetParameter(2)));
	  }
	}
      }
    }
  }

#if _PlotQA_
  TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,900,900);
  c_diff->Divide(3,3);
  for(int i_theta = 0; i_theta < 9; ++i_theta)
  {
    c_diff->cd(i_theta+1);
    c_diff->cd(i_theta+1)->SetLeftMargin(0.15);
    c_diff->cd(i_theta+1)->SetBottomMargin(0.15);
    c_diff->cd(i_theta+1)->SetTicks(1,1);
    c_diff->cd(i_theta+1)->SetGrid(0,0);
    if(i_theta < vmsa::CTS_stop)
    {
      string KEY_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_Func_%d",vmsa::pt_QA[energy],vmsa::Cent_start,vmsa::Eta_start,i_theta,vmsa::mPID[pid].c_str(),vmsa::Norm_QA,vmsa::Func_QA);
      h_mMass[KEY_QA]->SetTitle("");
      h_mMass[KEY_QA]->GetXaxis()->SetNdivisions(505,'N');
      h_mMass[KEY_QA]->GetXaxis()->SetLabelSize(0.03);
      h_mMass[KEY_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
      h_mMass[KEY_QA]->GetXaxis()->SetTitleSize(0.05);
      h_mMass[KEY_QA]->GetXaxis()->SetTitleOffset(1.2);
      h_mMass[KEY_QA]->GetXaxis()->CenterTitle();

      h_mMass[KEY_QA]->GetYaxis()->SetRangeUser(h_mMass[KEY_QA]->GetMinimum(),1.1*h_mMass[KEY_QA]->GetMaximum());
      h_mMass[KEY_QA]->GetYaxis()->SetNdivisions(505,'N');
      h_mMass[KEY_QA]->GetYaxis()->SetTitle("Yields");
      h_mMass[KEY_QA]->GetYaxis()->SetTitleSize(0.05);
      h_mMass[KEY_QA]->GetYaxis()->SetLabelSize(0.03);
      h_mMass[KEY_QA]->GetYaxis()->CenterTitle();

      h_mMass[KEY_QA]->SetMarkerStyle(24);
      h_mMass[KEY_QA]->SetMarkerColor(kGray+2);
      h_mMass[KEY_QA]->SetMarkerSize(1.2);
      h_mMass[KEY_QA]->Draw("pE");
      PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
    }
    if(i_theta == vmsa::CTS_stop)
    {
      string KEY_InteTheta_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_Norm_%d_Func_%d",vmsa::pt_QA[energy],vmsa::Cent_start,vmsa::Eta_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA,vmsa::Func_QA);
      h_mMass_InteTheta[KEY_InteTheta_QA]->SetTitle("");
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetNdivisions(505,'N');
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetLabelSize(0.03);
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleSize(0.05);
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleOffset(1.2);
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->CenterTitle();

      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetRangeUser(h_mMass_InteTheta[KEY_InteTheta_QA]->GetMinimum(),1.1*h_mMass_InteTheta[KEY_InteTheta_QA]->GetMaximum());
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetNdivisions(505,'N');
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetTitle("Yields");
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetTitleSize(0.05);
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetLabelSize(0.03);
      h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->CenterTitle();

      h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerStyle(24);
      h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerColor(kGray+2);
      h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerSize(1.2);
      h_mMass_InteTheta[KEY_InteTheta_QA]->Draw("pE");
      PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
      TF1 *f_bw = new TF1("f_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
      f_bw->SetParameter(0,Par_InteTheta[KEY_InteTheta_QA][0]);
      f_bw->SetParameter(1,Par_InteTheta[KEY_InteTheta_QA][1]);
      f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta_QA][2]);
      f_bw->SetLineColor(2);
      f_bw->SetLineStyle(2);
      f_bw->SetLineWidth(2);
      f_bw->Draw("l same");
    }
  }
#endif

  // fit theta-differential bin to extract fit parameters for integration
  vecFMap Par;
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
    {
      for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
      {
	for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	{
	  for(int i_func = vmsa::Func_start; i_func < vmsa::Func_stop; ++i_func)
	  {
	    for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	    {
	      string KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_Func_%d",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str(),i_norm,i_func);
	    }
	  }
	}
      }
    }
  }

  // extract counts vs. pT with diffenretial integration ranges and methods
  TH1FMap h_mCounts;
  vecFMap Par_rhoFit;
  TGraMap g_mRho;
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
	      g_mRho[KEY_rho] = new TGraphAsymmErrors();
	      for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt) // pt loop
	      {
		string KEY_counts = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_Norm_%d_Func_%d_Sigma_%d_%s",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str(),i_norm,i_func,i_sigma,vmsa::mInteMethod[i_method].c_str());
		h_mCounts[KEY_counts] = new TH1F(KEY_counts.c_str(),KEY_counts.c_str(),7,0.0,1.0);
		string KEY_InteTheta = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_Norm_%d_Func_%d",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str(),i_norm,i_func);
		for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; ++i_theta) // cos(theta*) loop
		{
		  string KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_Func_%d",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str(),i_norm,i_func);
		  TF1 *f_bw = new TF1("f_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
		  f_bw->FixParameter(0,Par_InteTheta[KEY_InteTheta][0]);
		  f_bw->FixParameter(1,Par_InteTheta[KEY_InteTheta][1]);
		  f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta][2]);
		  f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
		  h_mMass[KEY]->Fit(f_bw,"MQNR");
		  float bin_center = 1/14.0+i_theta/7.0;
		  if(i_method == 0)
		  {
		    int bin_start = h_mMass[KEY]->FindBin(Par_InteTheta[KEY_InteTheta][0]-vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]);
		    int bin_stop  = h_mMass[KEY]->FindBin(Par_InteTheta[KEY_InteTheta][0]+vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]);
		    float counts = 0.0;
		    float errors = 0.0;
		    for(int i_bin = bin_start; i_bin <= bin_stop; i_bin++)
		    {
		      counts += h_mMass[KEY]->GetBinContent(i_bin);
		      errors += h_mMass[KEY]->GetBinError(i_bin)*h_mMass[KEY]->GetBinError(i_bin);
		    }
		    h_mCounts[KEY_counts]->SetBinContent(h_mCounts[KEY_counts]->FindBin(bin_center),counts);
		    h_mCounts[KEY_counts]->SetBinError(h_mCounts[KEY_counts]->FindBin(bin_center),TMath::Sqrt(errors));
		  }
		  if(i_method == 1)
		  {
		    float bin_width = h_mMass[KEY]->GetBinWidth(1);
		    float Inte_start = Par_InteTheta[KEY_InteTheta][0]-vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]-0.5*bin_width;
		    float Inte_stop  = Par_InteTheta[KEY_InteTheta][0]+vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta][1]+0.5*bin_width;
		    float counts_bw = f_bw->Integral(Inte_start,Inte_stop)/bin_width;
		    float errors_bw = f_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
		    h_mCounts[KEY_counts]->SetBinContent(h_mCounts[KEY_counts]->FindBin(bin_center),counts_bw);
		    h_mCounts[KEY_counts]->SetBinError(h_mCounts[KEY_counts]->FindBin(bin_center),errors_bw);
		  }
		  Par[KEY].clear();
		  Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(0)));
		  Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(1)));
		  Par[KEY].push_back(static_cast<float>(f_bw->GetParameter(2)));
		}
		float pt_mean = (vmsa::pt_low[energy][i_pt]+vmsa::pt_up[energy][i_pt])/2.0;

		TF1 *f_rho = new TF1("f_rho",SpinDensity,0.0,1.0,2);
		f_rho->SetParameter(0,0.33);
		f_rho->SetParameter(1,h_mCounts[KEY_counts]->GetMaximum());
		h_mCounts[KEY_counts]->Fit(f_rho,"NMRI");
		Par_rhoFit[KEY_counts].clear();
		Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(0)));
		Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(1)));
		g_mRho[KEY_rho]->SetPoint(i_pt,pt_mean,f_rho->GetParameter(0));
		g_mRho[KEY_rho]->SetPointError(i_pt,0.0,0.0,f_rho->GetParError(0),f_rho->GetParError(0));
	      }
	    }
	  }
	}
      }
    }
  }
  
#if _PlotQA_
  for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; ++i_theta)
  {
    c_diff->cd(i_theta+1);
    string KEY_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_Func_%d",vmsa::pt_QA[energy],vmsa::Cent_start,vmsa::Eta_start,i_theta,vmsa::mPID[pid].c_str(),vmsa::Norm_QA,vmsa::Func_QA);
    TF1 *f_bw = new TF1("f_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
    f_bw->SetParameter(0,Par[KEY_QA][0]);
    f_bw->SetParameter(1,Par[KEY_QA][1]);
    f_bw->SetParameter(2,Par[KEY_QA][2]);
    f_bw->SetLineColor(4);
    f_bw->SetLineStyle(2);
    f_bw->SetLineWidth(2);
    f_bw->Draw("l same");
  }

  c_diff->cd(9);
  for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
  {
    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
    {
      string KEY_counts_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_Norm_%d_Func_%d_Sigma_%d_%s",vmsa::pt_QA[energy],vmsa::Cent_start,vmsa::Eta_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA,vmsa::Func_QA,i_sigma,vmsa::mInteMethod[i_method].c_str());
      if(i_sigma == vmsa::Sig_start && i_method == vmsa::Method_start)
      {
	h_mCounts[KEY_counts_QA];
	h_mCounts[KEY_counts_QA]->SetTitle("");
	h_mCounts[KEY_counts_QA]->GetXaxis()->SetNdivisions(505,'N');
	h_mCounts[KEY_counts_QA]->GetXaxis()->SetLabelSize(0.03);
	h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitle("cos(#theta*)");
	h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleSize(0.05);
	h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleOffset(1.2);
	h_mCounts[KEY_counts_QA]->GetXaxis()->CenterTitle();

	h_mCounts[KEY_counts_QA]->GetYaxis()->SetRangeUser(0.8*h_mCounts[KEY_counts_QA]->GetMinimum(),1.2*h_mCounts[KEY_counts_QA]->GetMaximum());
	h_mCounts[KEY_counts_QA]->GetYaxis()->SetNdivisions(505,'N');
	h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitle("Counts");
	h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitleSize(0.05);
	h_mCounts[KEY_counts_QA]->GetYaxis()->SetLabelSize(0.03);
	h_mCounts[KEY_counts_QA]->GetYaxis()->CenterTitle();

	h_mCounts[KEY_counts_QA]->SetMarkerStyle(20);
	h_mCounts[KEY_counts_QA]->SetMarkerColor(1);
	h_mCounts[KEY_counts_QA]->SetMarkerSize(1.2);
	h_mCounts[KEY_counts_QA]->Draw("pE");
      }
      else
      {
	h_mCounts[KEY_counts_QA]->SetMarkerStyle(24);
	h_mCounts[KEY_counts_QA]->SetMarkerColor(i_sigma+10*i_method+1);
	h_mCounts[KEY_counts_QA]->SetMarkerSize(1.2);
	h_mCounts[KEY_counts_QA]->Draw("pE same");
      }
      TF1 *f_rho = new TF1("f_rho",SpinDensity,0.0,1.0,2);
      f_rho->SetParameter(0,Par_rhoFit[KEY_counts_QA][0]);
      f_rho->SetParameter(1,Par_rhoFit[KEY_counts_QA][1]);
      f_rho->SetLineColor(i_sigma+10*i_method+1);
      f_rho->SetLineWidth(2);
      f_rho->SetLineStyle(2);
      f_rho->Draw("l same");
    }
  }
#endif

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

  string outputfile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/RawRhoPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
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
	      for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt) // pt loop
	      {
		string KEY_counts = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_Norm_%d_Func_%d_Sigma_%d_%s",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str(),i_norm,i_func,i_sigma,vmsa::mInteMethod[i_method].c_str());
		h_mCounts[KEY_counts]->Write();
	      }
	    }
	  }
	}
      }
    }
  }
  File_OutPut->Close();
}

#include <iostream>
#include <map>
#include <vector>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
#include "../Utility/functions.h"
#include "../Utility/draw.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooBreitWigner.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooCurve.h"
#include "RooHist.h"
#include "RooAbsReal.h"

using namespace RooFit;

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

void RooFitSpinAlignment(int energy = 6, int pid = 0, int year = 0)
{
  TGaxis::SetMaxDigits(4);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  double mMassKaon = 0.49368;
  // RooRealVar x("x","x",vmsa::InvMass_low[pid],vmsa::InvMass_high[pid]); // variable for inv. mass
  RooRealVar x("x","x",vmsa::BW_Start[pid],vmsa::BW_Stop[pid]); // variable for inv. mass
  // x.setRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
  // x.setRange("plotRange",vmsa::InvMass_low[pid],vmsa::InvMass_high[pid]);
  x.setRange("plotRange",vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);

  // string InPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/Yields/merged_file/Yields_SE_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string InPutFile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/Yields_SE_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(InPutFile.c_str());

  TH1FMap h_mInPut;
  RooHistMap r_mInput;
  for(int i_pt = vmsa::pt_start; i_pt < vmsa::pt_stop; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
    {
      for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
      {
	for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	{
	  string KEY_InPut = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_SE",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str());
	  h_mInPut[KEY_InPut] = (TH1F*)File_InPut->Get(KEY_InPut.c_str())->Clone(); 
	  h_mInPut[KEY_InPut]->SetMarkerStyle(24);
	  h_mInPut[KEY_InPut]->SetMarkerSize(0.8);
	  h_mInPut[KEY_InPut]->SetMarkerColor(kGray+2);
	  // h_mInPut[KEY_InPut]->GetXaxis()->SetRangeUser(x.getMin(),x.getMax());
	  h_mInPut[KEY_InPut]->GetYaxis()->SetRangeUser(-0.1*h_mInPut[KEY_InPut]->GetMaximum(),1.1*h_mInPut[KEY_InPut]->GetMaximum());
	  r_mInput[KEY_InPut] = new RooDataHist(KEY_InPut.c_str(),KEY_InPut.c_str(),x,h_mInPut[KEY_InPut]);
	}
      }
    }
  }

#if _PlotQA_
  // QA Plots for SE vs. ME
  TCanvas *c_peak = new TCanvas("c_peak","c_peak",10,10,1200,600);
  c_peak->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_peak->cd(i_pad+1);
    c_peak->cd(i_pad+1)->SetLeftMargin(0.15);
    c_peak->cd(i_pad+1)->SetBottomMargin(0.15);
    c_peak->cd(i_pad+1)->SetTicks(1,1);
    c_peak->cd(i_pad+1)->SetGrid(0,0);
  }
  c_peak->cd(1);
  string KEY_InPut_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_SE",vmsa::pt_RawQA[energy],vmsa::Cent_start,vmsa::Eta_start,vmsa::CTS_start,vmsa::mPID[pid].c_str());
  h_mInPut[KEY_InPut_QA]->GetYaxis()->SetRangeUser(-0.1*h_mInPut[KEY_InPut_QA]->GetMaximum(),1.1*h_mInPut[KEY_InPut_QA]->GetMaximum());
  h_mInPut[KEY_InPut_QA]->SetTitle("invMass from TH1F");

  h_mInPut[KEY_InPut_QA]->GetXaxis()->SetLabelSize(0.03);
  h_mInPut[KEY_InPut_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
  h_mInPut[KEY_InPut_QA]->GetXaxis()->SetTitleSize(0.05);
  h_mInPut[KEY_InPut_QA]->GetXaxis()->SetNdivisions(505);

  h_mInPut[KEY_InPut_QA]->GetYaxis()->SetLabelSize(0.03);
  h_mInPut[KEY_InPut_QA]->GetYaxis()->SetTitle("Yileds (arb. units)");
  h_mInPut[KEY_InPut_QA]->GetYaxis()->SetTitleSize(0.05);
  h_mInPut[KEY_InPut_QA]->GetYaxis()->SetTitleOffset(1.1);
  h_mInPut[KEY_InPut_QA]->GetYaxis()->SetNdivisions(505);

  h_mInPut[KEY_InPut_QA]->DrawCopy("PE");

  c_peak->cd(2);
  RooPlot *frame_peak = x.frame(Name("frame_peak"),Title("invMass from RooDataHist"),Range("x"));
  r_mInput[KEY_InPut_QA]->plotOn(frame_peak,MarkerStyle(24),MarkerSize(0.8),MarkerColor(kGray+2));
  frame_peak->SetLabelSize(0.03,"X");
  frame_peak->SetXTitle("M(K^{+},K^{-}) (GeV/c^{2})");
  frame_peak->SetTitleSize(0.05,"X");

  frame_peak->SetLabelSize(0.03,"Y");
  frame_peak->SetYTitle("Yileds (arb. units)");
  frame_peak->SetTitleSize(0.05,"Y");
  frame_peak->SetTitleOffset(1.1,"Y");
  frame_peak->SetAxisRange(h_mInPut[KEY_InPut_QA]->GetMinimum(),h_mInPut[KEY_InPut_QA]->GetMaximum(),"Y");
  frame_peak->Draw();

  // QA Plots for pT bins
  TCanvas *c_Pt = new TCanvas("c_Pt","c_Pt",10,10,1400,1400);
  c_Pt->Divide(5,5);
  for(int i_pt = vmsa::pt_start; i_pt < vmsa::pt_stop; i_pt++) // pt loop
  {
    c_Pt->cd(i_pt+1);
    c_Pt->cd(i_pt+1)->SetLeftMargin(0.15);
    c_Pt->cd(i_pt+1)->SetBottomMargin(0.15);
    c_Pt->cd(i_pt+1)->SetTicks(1,1);
    c_Pt->cd(i_pt+1)->SetGrid(0,0);
    string KEY_InPut_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_SE",i_pt,vmsa::Cent_start,vmsa::Eta_start,vmsa::CTS_start,vmsa::mPID[pid].c_str());
    h_mInPut[KEY_InPut_QA]->GetXaxis()->SetLabelSize(0.03);
    h_mInPut[KEY_InPut_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
    h_mInPut[KEY_InPut_QA]->GetXaxis()->SetTitleSize(0.05);
    h_mInPut[KEY_InPut_QA]->GetXaxis()->SetNdivisions(505);

    h_mInPut[KEY_InPut_QA]->GetYaxis()->SetLabelSize(0.03);
    h_mInPut[KEY_InPut_QA]->GetYaxis()->SetTitle("Yileds (arb. units)");
    h_mInPut[KEY_InPut_QA]->GetYaxis()->SetTitleSize(0.05);
    h_mInPut[KEY_InPut_QA]->GetYaxis()->SetTitleOffset(1.1);
    h_mInPut[KEY_InPut_QA]->GetYaxis()->SetNdivisions(505);
    h_mInPut[KEY_InPut_QA]->DrawCopy("PE");

    string pT_range = Form("[%.2f,%.2f]",vmsa::ptRawStart[i_pt],vmsa::ptRawStop[i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
  }

  TCanvas *c_RooPt = new TCanvas("c_RooPt","c_RooPt",10,10,1400,1400);
  c_RooPt->Divide(5,5);
  for(int i_pt = vmsa::pt_start; i_pt < vmsa::pt_stop; i_pt++) // pt loop
  {
    c_RooPt->cd(i_pt+1);
    c_RooPt->cd(i_pt+1)->SetLeftMargin(0.15);
    c_RooPt->cd(i_pt+1)->SetBottomMargin(0.15);
    c_RooPt->cd(i_pt+1)->SetTicks(1,1);
    c_RooPt->cd(i_pt+1)->SetGrid(0,0);
    string KEY_InPut_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_SE",i_pt,vmsa::Cent_start,vmsa::Eta_start,vmsa::CTS_start,vmsa::mPID[pid].c_str());
    RooPlot *frame_Pt = x.frame(Name(KEY_InPut_QA.c_str()),Title(KEY_InPut_QA.c_str()),Range("x"));
    r_mInput[KEY_InPut_QA]->plotOn(frame_Pt,MarkerStyle(24),MarkerSize(0.8),MarkerColor(kGray+2));
    frame_Pt->SetLabelSize(0.03,"X");
    frame_Pt->SetXTitle("M(K^{+},K^{-}) (GeV/c^{2})");
    frame_Pt->SetTitleSize(0.05,"X");

    frame_Pt->SetLabelSize(0.03,"Y");
    frame_Pt->SetYTitle("Yileds (arb. units)");
    frame_Pt->SetTitleSize(0.05,"Y");
    frame_Pt->SetTitleOffset(1.1,"Y");
    frame_Pt->SetAxisRange(h_mInPut[KEY_InPut_QA]->GetMinimum(),h_mInPut[KEY_InPut_QA]->GetMaximum(),"Y");
    frame_Pt->Draw();

    string pT_range = Form("[%.2f,%.2f]",vmsa::ptRawStart[i_pt],vmsa::ptRawStop[i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
  }
#endif

  // pT rebin
  TH1FMap h_mMass; // rebinned SE InvMass distribution
  RooHistMap r_mMass;
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    {
      for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
      {
	for(int pt_bin = vmsa::pt_rebin_first[energy]; pt_bin < vmsa::pt_rebin_last[energy]; pt_bin++) // pt loop
	{
	  string KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s",pt_bin,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str());
	  for(int i_pt = vmsa::pt_rebin_start[energy][pt_bin]; i_pt <= vmsa::pt_rebin_stop[energy][pt_bin]; i_pt++)
	  {
	    string KEY_InPut = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_SE",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str());
	    // cout << "KEY = " << KEY.c_str() << ", KEY_InPut = " << KEY_InPut.c_str() << endl;
	    if(i_pt == vmsa::pt_rebin_start[energy][pt_bin])
	    {
	      h_mMass[KEY] = (TH1F*)h_mInPut[KEY_InPut]->Clone();
	    }
	    else
	    {
	      h_mMass[KEY]->Add(h_mInPut[KEY_InPut],1.0);
	    }
	  }
	  r_mMass[KEY] = new RooDataHist(KEY.c_str(),KEY.c_str(),x,h_mMass[KEY]);
	}
      }
    }
  }

#if _PlotQA_
  // QA Plots for pT rebins
  TCanvas *c_pT_rebin = new TCanvas("c_pT_rebin","c_pT_rebin",10,10,1400,1400);
  c_pT_rebin->Divide(5,5);
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);
    string KEY_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s",i_pt,vmsa::Cent_start,vmsa::Eta_start,vmsa::CTS_start,vmsa::mPID[pid].c_str());
    h_mMass[KEY_QA]->GetYaxis()->SetRangeUser(-0.1*h_mMass[KEY_QA]->GetMaximum(),1.1*h_mMass[KEY_QA]->GetMaximum());
    h_mMass[KEY_QA]->DrawCopy("pE");
    // cout << "KEY_QA = " << KEY_QA.c_str() << endl;

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
    // PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
    PlotLine(x.getMin(),x.getMax(),0.0,0.0,1,2,2);
  }

  TCanvas *c_pT_rebinRoo = new TCanvas("c_pT_rebinRoo","c_pT_rebinRoo",10,10,1400,1400);
  c_pT_rebinRoo->Divide(5,5);
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    c_pT_rebinRoo->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
    c_pT_rebinRoo->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
    c_pT_rebinRoo->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
    c_pT_rebinRoo->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
    c_pT_rebinRoo->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);
    string KEY_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s",i_pt,vmsa::Cent_start,vmsa::Eta_start,vmsa::CTS_start,vmsa::mPID[pid].c_str());
    RooPlot *frame_pT_rebin = x.frame(Name(KEY_QA.c_str()),Title(KEY_QA.c_str()),Range("x"));
    r_mMass[KEY_QA]->plotOn(frame_pT_rebin,MarkerStyle(24),MarkerSize(0.8),MarkerColor(kGray+2));
    frame_pT_rebin->SetLabelSize(0.03,"X");
    frame_pT_rebin->SetXTitle("M(K^{+},K^{-}) (GeV/c^{2})");
    frame_pT_rebin->SetTitleSize(0.05,"X");

    frame_pT_rebin->SetLabelSize(0.03,"Y");
    frame_pT_rebin->SetYTitle("Yileds (arb. units)");
    frame_pT_rebin->SetTitleSize(0.05,"Y");
    frame_pT_rebin->SetTitleOffset(1.1,"Y");
    frame_pT_rebin->SetAxisRange(h_mMass[KEY_QA]->GetMinimum(),h_mMass[KEY_QA]->GetMaximum(),"Y");
    frame_pT_rebin->Draw();

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
    // PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
    PlotLine(x.getMin(),x.getMax(),0.0,0.0,1,2,2);
  }
#endif

  TH1FMap h_mInteMass; // integrated over cos(theta*)
  RooHistMap r_mInteMass; // integrated over cos(theta*)
  RooVarMap Invmass, Width, Nsig, C0, C1, C2, C3, C4, Nbkg;
  RooBWMap rf_mInteBW;
  RooGePdfMap rf_mInteDELPHI;
  RooPdfMap rf_mInteModel;
  TF1Map f_mInteDELPHI;
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
    {
      for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
      {
	string KEY_theta = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str());
	for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	{
	  string KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str());
	  if(i_theta == vmsa::CTS_start) h_mInteMass[KEY_theta] = (TH1F*)h_mMass[KEY]->Clone();
	  else h_mInteMass[KEY_theta]->Add(h_mMass[KEY],1.0);
	}
	r_mInteMass[KEY_theta] = new RooDataHist(KEY_theta.c_str(),KEY_theta.c_str(),x,h_mInteMass[KEY_theta]);

	int bin_start = h_mInteMass[KEY_theta]->FindBin(vmsa::BW_Start[pid]);
	int bin_stop  = h_mInteMass[KEY_theta]->FindBin(vmsa::BW_Stop[pid]);
	double yields = h_mInteMass[KEY_theta]->Integral(bin_start,bin_stop);
	// double yields = h_mInteMass[KEY_theta]->Integral();
	// breit wigner for signal
	Invmass[KEY_theta] = new RooRealVar("Invmass","Invmass",1.019,1.015,1.024); 
	Width[KEY_theta]   = new RooRealVar("Width","Width",0.0055,0.004,0.007);
	Nsig[KEY_theta]    = new RooRealVar("Nsig","signal fraction",0.05*yields,0.,yields);
	rf_mInteBW[KEY_theta] = new RooBreitWigner("sig","signal p.d.f.",x,*Invmass[KEY_theta],*Width[KEY_theta]);

	// DELPHI for background
	C0[KEY_theta]   = new RooRealVar("C0","C0", 0.5,-100.,100.);
	C1[KEY_theta]   = new RooRealVar("C1","C1",-0.5,-100.,100.);
	C2[KEY_theta]   = new RooRealVar("C2","C2",-0.5,-100.,100.);
	C3[KEY_theta]   = new RooRealVar("C3","C3",-0.5,-100.,100.);
	C4[KEY_theta]   = new RooRealVar("C4","C4", 0.5,-100.,100.);
	if(i_pt == 0) Nbkg[KEY_theta] = new RooRealVar("Nbkg","background fraction",0.01*yields,0.,yields);
	else Nbkg[KEY_theta] = new RooRealVar("Nbkg","background fraction",0.95*yields,0.,yields);
	rf_mInteDELPHI[KEY_theta] = new RooGenericPdf("bkg","bkg","pow((x-2.0*0.49368),C0)*exp(C1*x+C2*x*x+C3*x*x*x+C4*x*x*x*x)",RooArgSet(x,*C0[KEY_theta],*C1[KEY_theta],*C2[KEY_theta],*C3[KEY_theta],*C4[KEY_theta]));

	// bulid model: bw sig + DELPHI bkg
	rf_mInteModel[KEY_theta] = new RooAddPdf(KEY_theta.c_str(),KEY_theta.c_str(),RooArgSet(*rf_mInteBW[KEY_theta],*rf_mInteDELPHI[KEY_theta]),RooArgSet(*Nsig[KEY_theta],*Nbkg[KEY_theta]));
	RooFitResult* rooResult = rf_mInteModel[KEY_theta]->fitTo(*r_mInteMass[KEY_theta],Extended(kTRUE),Range("x"),Save(kTRUE));
	rooResult->Print();
	cout << "pt = " << i_pt << ", RooFitResult::covQual() = " << rooResult->covQual() << endl;

	// ROOT Fit
	f_mInteDELPHI[KEY_theta] = new TF1(KEY_theta.c_str(),DELPHIBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],9);
	f_mInteDELPHI[KEY_theta]->FixParameter(0,Invmass[KEY_theta]->getVal());
	f_mInteDELPHI[KEY_theta]->FixParameter(1,Width[KEY_theta]->getVal());
	// f_mInteDELPHI[KEY_theta]->SetParameter(2,0.002*h_mInteMass[KEY_theta]->GetMaximum());
	f_mInteDELPHI[KEY_theta]->SetParameter(2,1.);

	f_mInteDELPHI[KEY_theta]->SetParameter(3, 0.5);
	f_mInteDELPHI[KEY_theta]->SetParLimits(3,0.0,1.0);
	f_mInteDELPHI[KEY_theta]->SetParameter(4, 0.1);
	f_mInteDELPHI[KEY_theta]->SetParameter(5,-1.0);
	f_mInteDELPHI[KEY_theta]->SetParameter(6,-0.1);
	f_mInteDELPHI[KEY_theta]->SetParameter(7, 0.1);
	// f_mInteDELPHI[KEY_theta]->SetParameter(8,5.*h_mInteMass[KEY_theta]->GetMaximum());
	f_mInteDELPHI[KEY_theta]->SetParameter(8,1.);

	f_mInteDELPHI[KEY_theta]->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	h_mInteMass[KEY_theta]->Fit(f_mInteDELPHI[KEY_theta],"NR");
	double chi2 = f_mInteDELPHI[KEY_theta]->GetChisquare();
	double ndf = f_mInteDELPHI[KEY_theta]->GetNDF();
	cout << "pt = " << i_pt << ", chi2_ndf = " << chi2 << "/" << ndf << " = " << chi2/ndf << endl;
      }
    }
  }

#if _PlotQA_
  // QA: breit wigner + DELPHI fits to cos(theta*) integrated InvMass
  TCanvas *c_mInteMass = new TCanvas("c_mInteMass","c_mInteMass",10,10,1400,1400);
  c_mInteMass->Divide(5,5);
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    c_mInteMass->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
    c_mInteMass->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
    c_mInteMass->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
    c_mInteMass->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
    c_mInteMass->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);
    string KEY_theta_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s",i_pt,vmsa::Cent_start,vmsa::Eta_start,vmsa::mPID[pid].c_str());
    h_mInteMass[KEY_theta_QA]->GetYaxis()->SetRangeUser(-0.1*h_mInteMass[KEY_theta_QA]->GetMaximum(),1.1*h_mInteMass[KEY_theta_QA]->GetMaximum());
    h_mInteMass[KEY_theta_QA]->DrawCopy("PE");
    f_mInteDELPHI[KEY_theta_QA]->SetLineColor(2);
    f_mInteDELPHI[KEY_theta_QA]->SetLineStyle(1);
    f_mInteDELPHI[KEY_theta_QA]->SetLineWidth(2);
    f_mInteDELPHI[KEY_theta_QA]->Draw("l same");

    TF1 *f_bkg = new TF1("f_bkg",DELPHI,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
    for(int i_par = 0; i_par < 6; ++i_par)
    {
      f_bkg->FixParameter(i_par,f_mInteDELPHI[KEY_theta_QA]->GetParameter(i_par+3));
    }
    f_bkg->SetLineColor(4);
    f_bkg->SetLineStyle(2);
    f_bkg->SetLineWidth(2);
    f_bkg->Draw("l same");

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
    // PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
    PlotLine(x.getMin(),x.getMax(),0.0,0.0,1,2,2);
  }

  TCanvas *c_mInteMassRoo = new TCanvas("c_mInteMassRoo","c_mInteMassRoo",10,10,1400,1400);
  c_mInteMassRoo->Divide(5,5);
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    c_mInteMassRoo->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
    c_mInteMassRoo->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
    c_mInteMassRoo->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
    c_mInteMassRoo->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
    c_mInteMassRoo->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);
    string KEY_theta_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s",i_pt,vmsa::Cent_start,vmsa::Eta_start,vmsa::mPID[pid].c_str());
    RooPlot *frame_InteMass = x.frame(Name(KEY_theta_QA.c_str()),Title(KEY_theta_QA.c_str()),Range("plotRange"));
    r_mInteMass[KEY_theta_QA]->plotOn(frame_InteMass,Range("plotRange"),MarkerStyle(24),MarkerSize(0.8),MarkerColor(kGray+1));
    rf_mInteModel[KEY_theta_QA]->plotOn(frame_InteMass,Range("x"),LineColor(kRed));
    rf_mInteModel[KEY_theta_QA]->plotOn(frame_InteMass,Range("x"),Components(*rf_mInteDELPHI[KEY_theta_QA]),LineColor(kAzure+2),LineStyle(kDashed));
    rf_mInteModel[KEY_theta_QA]->plotOn(frame_InteMass,Range("x"),Components(*rf_mInteBW[KEY_theta_QA]),LineColor(kRed+2),LineStyle(kDashed));

    frame_InteMass->SetLabelSize(0.03,"X");
    frame_InteMass->SetXTitle("M(K^{+},K^{-}) (GeV/c^{2})");
    frame_InteMass->SetTitleSize(0.05,"X");

    frame_InteMass->SetLabelSize(0.03,"Y");
    frame_InteMass->SetYTitle("Yileds (arb. units)");
    frame_InteMass->SetTitleSize(0.05,"Y");
    frame_InteMass->SetTitleOffset(1.1,"Y");
    frame_InteMass->SetAxisRange(h_mInteMass[KEY_theta_QA]->GetMinimum(),h_mInteMass[KEY_theta_QA]->GetMaximum(),"Y");

    frame_InteMass->Draw();

    RooCurve *curve = (RooCurve*)frame_InteMass->getObject(1);
    RooHist *histo = (RooHist*)frame_InteMass->getHist(("h_"+KEY_theta_QA).c_str());
    double chi2_ndf = curve->chiSquare(*histo,9);
    cout << "pt = " << i_pt << ", chi2_ndf = " << chi2_ndf << endl;

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
    // PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
    PlotLine(x.getMin(),x.getMax(),0.0,0.0,1,2,2);
  }
#endif

  // fit invmass vs. cos(theta*) bin with different fitting
  RooVarMap invmass, width, nsig, c0, c1, c2, c3, c4, nbkg;
  RooBWMap rf_mBW;
  RooGePdfMap rf_mDELPHI;
  RooPdfMap rf_mModel;
  TF1Map f_mDELPHI;
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    {
      for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
      {
	string KEY_theta = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str());
	for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	{
	  string KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str());

	  // RooFit
	  invmass[KEY] = new RooRealVar("invmass","invmass",Invmass[KEY_theta]->getVal(),1.015,1.024); 
	  width[KEY]   = new RooRealVar("width","width",Width[KEY_theta]->getVal(),0.004,0.007);
	  nsig[KEY]    = new RooRealVar("nsig","signal fraction",0.15*Nsig[KEY_theta]->getVal(),0.,0.5*Nsig[KEY_theta]->getVal());
	  rf_mBW[KEY] = new RooBreitWigner("sig","signal p.d.f.",x,*invmass[KEY],*width[KEY]);
	  invmass[KEY]->setConstant(kTRUE);
	  width[KEY]->setConstant(kTRUE);

	  // DELPHI for background
	  c0[KEY]   = new RooRealVar("c0","c0",C0[KEY_theta]->getVal(),0.,1.);
	  c1[KEY]   = new RooRealVar("c1","c1",C1[KEY_theta]->getVal(),-100.,100.);
	  c2[KEY]   = new RooRealVar("c2","c2",C2[KEY_theta]->getVal(),-100.,100.);
	  c3[KEY]   = new RooRealVar("c3","c3",C3[KEY_theta]->getVal(),-100.,100.);
	  c4[KEY]   = new RooRealVar("c4","c4",C4[KEY_theta]->getVal(),-100.,100.);
	  nbkg[KEY] = new RooRealVar("nbkg","background fraction",0.15*Nbkg[KEY_theta]->getVal(),0.,0.5*Nbkg[KEY_theta]->getVal());
	  rf_mDELPHI[KEY] = new RooGenericPdf("bkg","bkg","pow((x-2.0*0.49368),c0)*exp(c1*x+c2*x*x+c3*x*x*x+c4*x*x*x*x)",RooArgSet(x,*c0[KEY],*c1[KEY],*c2[KEY],*c3[KEY],*c4[KEY]));

	  // bulid model: bw sig + DELPHI bkg
	  rf_mModel[KEY] = new RooAddPdf(KEY.c_str(),KEY.c_str(),RooArgSet(*rf_mBW[KEY],*rf_mDELPHI[KEY]),RooArgSet(*nsig[KEY],*nbkg[KEY]));
	  RooFitResult* rooResult = rf_mModel[KEY]->fitTo(*r_mMass[KEY],Extended(kTRUE),Range("x"),Save(kTRUE));
	  rooResult->Print();
	  cout << "RooFit: pt = " << i_pt << ", i_theta = " << i_theta << ", RooFitResult::covQual() = " << rooResult->covQual() << endl;

	  // ROOT Fit
	  f_mDELPHI[KEY] = new TF1(KEY.c_str(),DELPHIBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],9);
	  f_mDELPHI[KEY]->FixParameter(0,f_mInteDELPHI[KEY_theta]->GetParameter(0));
	  f_mDELPHI[KEY]->FixParameter(1,f_mInteDELPHI[KEY_theta]->GetParameter(1));
	  f_mDELPHI[KEY]->SetParameter(2,0.1*f_mInteDELPHI[KEY_theta]->GetParameter(2));

	  f_mDELPHI[KEY]->SetParameter(3,f_mInteDELPHI[KEY_theta]->GetParameter(3));
	  f_mDELPHI[KEY]->SetParLimits(3,0.0,1.0);
	  f_mDELPHI[KEY]->SetParameter(4,f_mInteDELPHI[KEY_theta]->GetParameter(4));
	  f_mDELPHI[KEY]->SetParameter(5,f_mInteDELPHI[KEY_theta]->GetParameter(5));
	  f_mDELPHI[KEY]->SetParameter(6,f_mInteDELPHI[KEY_theta]->GetParameter(6));
	  f_mDELPHI[KEY]->SetParameter(7,f_mInteDELPHI[KEY_theta]->GetParameter(7));
	  f_mDELPHI[KEY]->SetParameter(8,0.1*f_mInteDELPHI[KEY_theta]->GetParameter(8));
	  f_mDELPHI[KEY]->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	  h_mMass[KEY]->Fit(f_mDELPHI[KEY],"NR");
	  double chi2 = f_mDELPHI[KEY]->GetChisquare();
	  double ndf = f_mDELPHI[KEY]->GetNDF();
	  cout << "ROOT: pt = " << i_pt << ", i_theta = " << i_theta << ", chi2_ndf = " << chi2 << "/" << ndf << " = " << chi2/ndf << endl;
	}
      }
    }
  }

#if _PlotQA_
  // QA InvMass vs. phi for gaussian and breit wigner fits
  TCanvas *c_mMass[vmsa::pt_rebin_last[energy]]; 
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt)
  {
    string CanName = Form("c_mMass_%d",i_pt);
    c_mMass[i_pt] = new TCanvas(CanName.c_str(),CanName.c_str(),10,10,900,900);
    c_mMass[i_pt]->Divide(3,3);
    for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
    {
      c_mMass[i_pt]->cd(i_theta+1)->SetLeftMargin(0.15);
      c_mMass[i_pt]->cd(i_theta+1)->SetBottomMargin(0.15);
      c_mMass[i_pt]->cd(i_theta+1)->SetTicks(1,1);
      c_mMass[i_pt]->cd(i_theta+1)->SetGrid(0,0);
      c_mMass[i_pt]->cd(i_theta+1);
      string KEY_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s",i_pt,vmsa::Cent_start,vmsa::Eta_start,i_theta,vmsa::mPID[pid].c_str());
      h_mMass[KEY_QA]->SetMarkerColor(1);
      h_mMass[KEY_QA]->SetMarkerStyle(24);
      h_mMass[KEY_QA]->SetMarkerSize(0.8);
      h_mMass[KEY_QA]->DrawCopy("pE");

      f_mDELPHI[KEY_QA]->SetLineColor(2);
      f_mDELPHI[KEY_QA]->SetLineWidth(2);
      f_mDELPHI[KEY_QA]->Draw("l same");

      TF1 *f_bkg = new TF1("f_bkg",DELPHI,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
      for(int i_par = 0; i_par < 6; ++i_par)
      {
	f_bkg->FixParameter(i_par,f_mDELPHI[KEY_QA]->GetParameter(i_par+3));
      }
      f_bkg->SetLineColor(4);
      f_bkg->SetLineStyle(2);
      f_bkg->SetLineWidth(2);
      f_bkg->Draw("l same");
    }
    c_mMass[i_pt]->SaveAs((CanName+".eps").c_str());
  }

  TCanvas *c_mMassRoo[vmsa::pt_rebin_last[energy]];
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt)
  {
    string CanName = Form("c_mMassRoo_%d",i_pt);
    c_mMassRoo[i_pt] = new TCanvas(CanName.c_str(),CanName.c_str(),10,10,900,900);
    c_mMassRoo[i_pt]->Divide(3,3);
    for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
    {
      c_mMassRoo[i_pt]->cd(i_theta+1)->SetLeftMargin(0.15);
      c_mMassRoo[i_pt]->cd(i_theta+1)->SetBottomMargin(0.15);
      c_mMassRoo[i_pt]->cd(i_theta+1)->SetTicks(1,1);
      c_mMassRoo[i_pt]->cd(i_theta+1)->SetGrid(0,0);
      c_mMassRoo[i_pt]->cd(i_theta+1);
      string KEY_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s",i_pt,vmsa::Cent_start,vmsa::Eta_start,i_theta,vmsa::mPID[pid].c_str());
      RooPlot *frame_Mass = x.frame(Name(KEY_QA.c_str()),Title(KEY_QA.c_str()),Range("plotRange"));
      r_mMass[KEY_QA]->plotOn(frame_Mass,Range("plotRange"),MarkerStyle(24),MarkerSize(0.8),MarkerColor(kGray+2));
      rf_mModel[KEY_QA]->plotOn(frame_Mass,Range("x"),LineColor(kRed));
      rf_mModel[KEY_QA]->plotOn(frame_Mass,Range("x"),Components(*rf_mDELPHI[KEY_QA]),LineColor(kAzure+2),LineStyle(kDashed));
      rf_mModel[KEY_QA]->plotOn(frame_Mass,Range("x"),Components(*rf_mBW[KEY_QA]),LineColor(kRed+2),LineStyle(kDashed));

      frame_Mass->SetLabelSize(0.03,"X");
      frame_Mass->SetXTitle("M(K^{+},K^{-}) (GeV/c^{2})");
      frame_Mass->SetTitleSize(0.05,"X");

      frame_Mass->SetLabelSize(0.03,"Y");
      frame_Mass->SetYTitle("Yileds (arb. units)");
      frame_Mass->SetTitleSize(0.05,"Y");
      frame_Mass->SetTitleOffset(1.1,"Y");
      // frame_Mass->SetAxisRange(h_mMass[KEY_QA]->GetMinimum(),h_mMass[KEY_QA]->GetMaximum(),"Y");
      frame_Mass->SetAxisRange(-0.1*h_mMass[KEY_QA]->GetMaximum(),1.1*h_mMass[KEY_QA]->GetMaximum(),"Y");

      frame_Mass->Draw();
      PlotLine(x.getMin(),x.getMax(),0.0,0.0,1,2,2);

      RooCurve *curve = (RooCurve*)frame_Mass->getObject(1);
      RooHist *histo = (RooHist*)frame_Mass->getHist(("h_"+KEY_QA).c_str());
      double chi2_ndf = curve->chiSquare(*histo,9);
      cout << "i_pt = " << i_pt << ", theta = " << i_theta << ", chi2_ndf = " << chi2_ndf << endl;
    }
    c_mMassRoo[i_pt]->SaveAs((CanName+".eps").c_str());
  }
#endif

  double nSig = 2.0;
  TH1FMap h_mCounts;
  vecFMap ParSpin_Count, ParSpin_RooBW;
  TGraMap g_mRho00;
  // calculate counts from different fitting method and extract raw rho_00
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    {
      string KEY_Rho00_RooBW = Form("Rho00_RooBW_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str()); // breit wigner fits
      g_mRho00[KEY_Rho00_RooBW] = new TGraphAsymmErrors();
      string KEY_Rho00_Count = Form("Rho00_Count_Centrality_%d_EtaGap_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str()); // bin counting
      g_mRho00[KEY_Rho00_Count] = new TGraphAsymmErrors();
      for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
      {
	string KEY_RooBW = Form("RooBW_pt_%d_Centrality_%d_EtaGap_%d_2nd_%s",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str());
	h_mCounts[KEY_RooBW] = new TH1F(KEY_RooBW.c_str(),KEY_RooBW.c_str(),7,0.0,1.0);
	string KEY_Count = Form("Count_pt_%d_Centrality_%d_EtaGap_%d_2nd_%s",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str());
	h_mCounts[KEY_Count] = new TH1F(KEY_Count.c_str(),KEY_Count.c_str(),7,0.0,1.0);

	string KEY_theta = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str());
	double mean = Invmass[KEY_theta]->getVal();
	double sigma = Width[KEY_theta]->getVal();
	x.setRange(KEY_theta.c_str(),mean-nSig*sigma,mean+nSig*sigma);
	for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	{
	  string KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str());
	  RooAbsReal *inteSig = rf_mBW[KEY]->createIntegral(x,NormSet(x),Range(KEY_theta.c_str()));
	  double counts_Roo = inteSig->getVal()*nsig[KEY]->getVal();
	  double errors_Roo = inteSig->getVal()*nsig[KEY]->getError();
	  // cout << "i_pt = " << i_pt << ", i_theta = " << i_theta << ", inteSig = " << counts_Roo << " +/- " << errors_Roo << endl;
	  float bin_center = 1/14.0+i_theta/7.0;
	  h_mCounts[KEY_RooBW]->SetBinContent(h_mCounts[KEY_RooBW]->FindBin(bin_center),counts_Roo);
	  h_mCounts[KEY_RooBW]->SetBinError(h_mCounts[KEY_RooBW]->FindBin(bin_center),errors_Roo);

	  TF1 *f_bkg = new TF1("f_bkg",DELPHI,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
	  for(int i_par = 0; i_par < 6; ++i_par)
	  {
	    f_bkg->FixParameter(i_par,f_mDELPHI[KEY]->GetParameter(i_par+3));
	  }
	  h_mMass[KEY]->Add(f_bkg,-1.0); // subtract background

	  float counts = 0.0;
	  float errors = 0.0;
	  int bin_start = h_mMass[KEY]->FindBin(mean-nSig*sigma);
	  int bin_stop  = h_mMass[KEY]->FindBin(mean+nSig*sigma);
	  for(int i_bin = bin_start; i_bin <= bin_stop; i_bin++)
	  {
	    counts += h_mMass[KEY]->GetBinContent(i_bin);
	    errors += h_mMass[KEY]->GetBinError(i_bin)*h_mMass[KEY]->GetBinError(i_bin);
	  }
	  h_mCounts[KEY_Count]->SetBinContent(h_mCounts[KEY_Count]->FindBin(bin_center),counts);
	  h_mCounts[KEY_Count]->SetBinError(h_mCounts[KEY_Count]->FindBin(bin_center),TMath::Sqrt(errors));
	  // cout << "i_pt = " << i_pt << ", i_theta = " << i_theta << ", countSig = " << counts << " +/- " << TMath::Sqrt(errors) << endl; cout << endl;
	}

	float pt_mean = (vmsa::pt_low[energy][i_pt]+vmsa::pt_up[energy][i_pt])/2.0;

	TF1 *f_cos_roobw = new TF1("f_cos_roobw",SpinDensity,0.0,1.0,2);
	f_cos_roobw->SetParameter(0,0.33);
	f_cos_roobw->SetParameter(1,1000);
	h_mCounts[KEY_RooBW]->Fit(f_cos_roobw,"NMQI");
	ParSpin_RooBW[KEY_RooBW].clear();
	ParSpin_RooBW[KEY_RooBW].push_back(static_cast<float>(f_cos_roobw->GetParameter(0)));
	ParSpin_RooBW[KEY_RooBW].push_back(static_cast<float>(f_cos_roobw->GetParameter(1)));
	g_mRho00[KEY_Rho00_RooBW]->SetPoint(i_pt,pt_mean,f_cos_roobw->GetParameter(0));
	g_mRho00[KEY_Rho00_RooBW]->SetPointError(i_pt,0.0,0.0,f_cos_roobw->GetParError(0),f_cos_roobw->GetParError(0));

	TF1 *f_cos_count = new TF1("f_cos_count",SpinDensity,0.0,1.0,2);
	f_cos_count->SetParameter(0,0.33);
	f_cos_count->SetParameter(1,1000);
	h_mCounts[KEY_Count]->Fit(f_cos_count,"NMQI");
	ParSpin_Count[KEY_Count].clear();
	ParSpin_Count[KEY_Count].push_back(static_cast<float>(f_cos_count->GetParameter(0)));
	ParSpin_Count[KEY_Count].push_back(static_cast<float>(f_cos_count->GetParameter(1)));
	g_mRho00[KEY_Rho00_Count]->SetPoint(i_pt,pt_mean,f_cos_count->GetParameter(0));
	g_mRho00[KEY_Rho00_Count]->SetPointError(i_pt,0.0,0.0,f_cos_count->GetParError(0),f_cos_count->GetParError(0));
      }
    }
  }

#if _PlotQA_
  // QA InvMass vs. phi for gaussian and breit wigner fits
  TCanvas *c_mMassSig[vmsa::pt_rebin_last[energy]]; 
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt)
  {
    string CanName = Form("c_mMassSig_%d",i_pt);
    c_mMassSig[i_pt] = new TCanvas(CanName.c_str(),CanName.c_str(),10,10,900,900);
    c_mMassSig[i_pt]->Divide(3,3);
    string KEY_theta_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s",i_pt,0,0,vmsa::mPID[pid].c_str());
    double mean = Invmass[KEY_theta_QA]->getVal();
    double sigma = Width[KEY_theta_QA]->getVal();
    for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
    {
      c_mMassSig[i_pt]->cd(i_theta+1)->SetLeftMargin(0.15);
      c_mMassSig[i_pt]->cd(i_theta+1)->SetBottomMargin(0.15);
      c_mMassSig[i_pt]->cd(i_theta+1)->SetTicks(1,1);
      c_mMassSig[i_pt]->cd(i_theta+1)->SetGrid(0,0);
      c_mMassSig[i_pt]->cd(i_theta+1);
      string KEY_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s",i_pt,vmsa::Cent_start,vmsa::Eta_start,i_theta,vmsa::mPID[pid].c_str());
      h_mMass[KEY_QA]->GetXaxis()->SetRangeUser(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
      h_mMass[KEY_QA]->SetMarkerColor(1);
      h_mMass[KEY_QA]->SetMarkerStyle(24);
      h_mMass[KEY_QA]->SetMarkerSize(0.8);
      h_mMass[KEY_QA]->DrawCopy("pE");
      PlotLine(vmsa::BW_Start[pid],vmsa::BW_Stop[pid],0.0,0.0,1,2,2);
      PlotLine(mean-nSig*sigma,mean-nSig*sigma,0.0,0.5*h_mMass[KEY_QA]->GetMaximum(),4,2,2);
      PlotLine(mean+nSig*sigma,mean+nSig*sigma,0.0,0.5*h_mMass[KEY_QA]->GetMaximum(),4,2,2);
    }

    c_mMassSig[i_pt]->cd(8)->SetLeftMargin(0.15);
    c_mMassSig[i_pt]->cd(8)->SetBottomMargin(0.15);
    c_mMassSig[i_pt]->cd(8)->SetTicks(1,1);
    c_mMassSig[i_pt]->cd(8)->SetGrid(0,0);
    c_mMassSig[i_pt]->cd(8);

    string KEY_RooBW_QA = Form("RooBW_pt_%d_Centrality_%d_EtaGap_%d_2nd_%s",i_pt,0,0,vmsa::mPID[pid].c_str());
    h_mCounts[KEY_RooBW_QA]->SetStats(0);
    h_mCounts[KEY_RooBW_QA]->SetTitle("");
    h_mCounts[KEY_RooBW_QA]->SetTitleSize(0.08);
    h_mCounts[KEY_RooBW_QA]->GetXaxis()->SetTitle("cos(#theta*) (w.r.t. 2^{nd} event plane)");
    h_mCounts[KEY_RooBW_QA]->GetXaxis()->SetTitleSize(0.06);
    h_mCounts[KEY_RooBW_QA]->GetXaxis()->CenterTitle();
    h_mCounts[KEY_RooBW_QA]->GetXaxis()->SetNdivisions(505);
    h_mCounts[KEY_RooBW_QA]->GetXaxis()->SetLabelSize(0.05);

    h_mCounts[KEY_RooBW_QA]->GetYaxis()->SetTitle("Yields");
    h_mCounts[KEY_RooBW_QA]->GetYaxis()->SetTitleSize(0.06);
    h_mCounts[KEY_RooBW_QA]->GetYaxis()->CenterTitle();
    h_mCounts[KEY_RooBW_QA]->GetYaxis()->SetNdivisions(505);
    h_mCounts[KEY_RooBW_QA]->GetYaxis()->SetLabelSize(0.05);
    h_mCounts[KEY_RooBW_QA]->SetLineColor(2);
    h_mCounts[KEY_RooBW_QA]->SetMarkerColor(2);
    h_mCounts[KEY_RooBW_QA]->SetMarkerStyle(24);
    h_mCounts[KEY_RooBW_QA]->SetMarkerSize(1.2);
    h_mCounts[KEY_RooBW_QA]->DrawCopy("pE");
    TF1 *f_cos_roobw_QA = new TF1("f_cos_roobw_QA",SpinDensity,0.0,1.0,2);
    f_cos_roobw_QA->FixParameter(0,ParSpin_RooBW[KEY_RooBW_QA][0]);
    f_cos_roobw_QA->FixParameter(1,ParSpin_RooBW[KEY_RooBW_QA][1]);
    f_cos_roobw_QA->SetLineStyle(2);
    f_cos_roobw_QA->SetLineColor(2);
    f_cos_roobw_QA->DrawCopy("l same");


    string KEY_Count_QA = Form("Count_pt_%d_Centrality_%d_EtaGap_%d_2nd_%s",i_pt,0,0,vmsa::mPID[pid].c_str());
    h_mCounts[KEY_Count_QA]->SetLineColor(4);
    h_mCounts[KEY_Count_QA]->SetMarkerColor(4);
    h_mCounts[KEY_Count_QA]->SetMarkerStyle(24);
    h_mCounts[KEY_Count_QA]->SetMarkerSize(1.2);
    h_mCounts[KEY_Count_QA]->DrawCopy("pE same");
    TF1 *f_cos_count_QA = new TF1("f_cos_count_QA",SpinDensity,0.0,1.0,2);
    f_cos_count_QA->FixParameter(0,ParSpin_Count[KEY_Count_QA][0]);
    f_cos_count_QA->FixParameter(1,ParSpin_Count[KEY_Count_QA][1]);
    f_cos_count_QA->SetLineStyle(2);
    f_cos_count_QA->SetLineColor(4);
    f_cos_count_QA->DrawCopy("l same");
    c_mMassSig[i_pt]->SaveAs((CanName+".eps").c_str());
  }
#endif

  // QA Plots for rho00 vs. pt of bin counting and breit wigner integrating
  TCanvas *c_Rho00_pT = new TCanvas("c_Rho00_pT","c_Rho00_pT",10,10,800,800);
  c_Rho00_pT->cd();
  c_Rho00_pT->cd()->SetLeftMargin(0.15);
  c_Rho00_pT->cd()->SetBottomMargin(0.15);
  c_Rho00_pT->cd()->SetTicks(1,1);
  c_Rho00_pT->cd()->SetGrid(0,0);
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
  h_play->DrawCopy("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  string KEY_Rho00_RooBW_QA = Form("Rho00_RooBW_Centrality_%d_EtaGap_%d_%s",vmsa::Cent_start,vmsa::Eta_start,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00[KEY_Rho00_RooBW_QA],24,2,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.43,0.0,0.0,0.0,0.0,24,2,1.3);
  string leg_bw = "#phi (RooFit bw integrating)";
  plotTopLegend((char*)leg_bw.c_str(),0.6,0.427,0.03,1,0.0,42,0);

  string KEY_Rho00_Count_QA = Form("Rho00_Count_Centrality_%d_EtaGap_%d_%s",vmsa::Cent_start,vmsa::Eta_start,vmsa::mPID[pid].c_str());
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho00[KEY_Rho00_Count_QA],24,4,1.1);
  Draw_TGAE_Point_new_Symbol(0.5,0.45,0.0,0.0,0.0,0.0,24,4,1.3);
  string leg_count = "#phi (bin counting)";
  plotTopLegend((char*)leg_count.c_str(),0.6,0.447,0.03,1,0.0,42,0);

  PlotLine(3.0,3.5,0.45,0.45,1,2,2);
  string leg_line = "#rho_{00} = 1/3";
  plotTopLegend((char*)leg_line.c_str(),3.6,0.447,0.03,1,0.0,42,0);

  string leg_energy = Form("AuAu %s (%s)", vmsa::mBeamEnergy[energy].c_str(), vmsa::mYear[year].c_str());
  plotTopLegend((char*)leg_energy.c_str(),1.4,0.27,0.04,1,0.0,42,0);
  string leg_centrality = "20%-60%";
  plotTopLegend((char*)leg_centrality.c_str(),2.0,0.24,0.04,1,0.0,42,0);

  string figure_name = Form("./figures/rho00_pT_%s.eps",vmsa::mBeamEnergy[energy].c_str());
  c_Rho00_pT->SaveAs(figure_name.c_str());
}

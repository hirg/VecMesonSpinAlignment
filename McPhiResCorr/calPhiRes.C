#include <iostream>
#include <string>
#include <map>
#include <vector>
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TProfile.h"
#include "../Utility/functions.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"

using namespace std;

// typedef std::map<std::string,TH1F*> TH1FMap;
// typedef std::map<TString,TProfile*> TProMap;
// typedef std::map<string,std::vector<float> > vecFMap;

void calPhiRes(int energy = 6, int pid = 0)
{
  TGaxis::SetMaxDigits(4);

  // use v3 histogram right now
  string InPutFile_SE = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/%s/flow_%s/merged_file/Yields_SE_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  cout << "InPut Same Event: " << InPutFile_SE.c_str() << endl;
  string InPutFile_ME = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/%s/flow_%s/merged_file/Yields_ME_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  cout << "InPut Mixed Event: " << InPutFile_ME.c_str() << endl;
  
  TFile *File_SE = TFile::Open(InPutFile_SE.c_str());
  TFile *File_ME = TFile::Open(InPutFile_ME.c_str());

  TH1FMap h_mYield_SE, h_mYield_ME, h_mYield;
  for(int i_cent = 0; i_cent < 9; ++i_cent) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; ++i_eta) // EtaGap loop
    {
      for(int i_sys = vmsa::Sys_start; i_sys < vmsa::Sys_stop; ++i_sys) // Systematic loop
      {
	string KEY_Yield_SE = Form("Yields_Centrality_%d_EtaGap_%d_%s_SE_SysErrors_%d",i_cent,i_eta,vmsa::mPID[pid].c_str(),i_sys);
	h_mYield_SE[KEY_Yield_SE] = (TH1F*)File_SE->Get(KEY_Yield_SE.c_str())->Clone(); 
	int Norm_bin_start = h_mYield_SE[KEY_Yield_SE]->FindBin(vmsa::Norm_Start[pid]);
	int Norm_bin_stop  = h_mYield_SE[KEY_Yield_SE]->FindBin(vmsa::Norm_Stop[pid]);
	float Inte_SE = h_mYield_SE[KEY_Yield_SE]->Integral(Norm_bin_start,Norm_bin_stop);

	string KEY_Yield_ME = Form("Yields_Centrality_%d_EtaGap_%d_%s_ME_SysErrors_%d",i_cent,i_eta,vmsa::mPID[pid].c_str(),i_sys);
	h_mYield_ME[KEY_Yield_ME] = (TH1F*)File_ME->Get(KEY_Yield_ME.c_str())->Clone(); 
	float Inte_ME = h_mYield_ME[KEY_Yield_ME]->Integral(Norm_bin_start,Norm_bin_stop);
	h_mYield_ME[KEY_Yield_ME]->Scale(Inte_SE/Inte_ME);

	string KEY_Yield = Form("Yields_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,i_eta,vmsa::mPID[pid].c_str(),i_sys);
	h_mYield[KEY_Yield] = (TH1F*)h_mYield_SE[KEY_Yield_SE]->Clone();
	h_mYield[KEY_Yield]->Add(h_mYield_ME[KEY_Yield_ME],-1.0);
      }
    }
  }

  /*
  TCanvas *c_RawYields = new TCanvas("c_RawYields","c_RawYields",10,10,900,900);
  c_RawYields->Divide(3,3);
  for(int i_pad = 0; i_pad < 9; ++i_pad)
  {
    c_RawYields->cd(i_pad+1);
    c_RawYields->cd(i_pad+1)->SetLeftMargin(0.15);
    c_RawYields->cd(i_pad+1)->SetBottomMargin(0.15);
    c_RawYields->cd(i_pad+1)->SetTicks(1,1);
    c_RawYields->cd(i_pad+1)->SetGrid(0,0);

    string KEY_Yield_SE = Form("Yields_Centrality_%d_EtaGap_%d_%s_SE_SysErrors_%d",i_pad,vmsa::Eta_QA,vmsa::mPID[pid].c_str(),vmsa::Sys_QA);
    h_mYield_SE[KEY_Yield_SE]->SetTitle("");
    h_mYield_SE[KEY_Yield_SE]->SetStats(0);
    h_mYield_SE[KEY_Yield_SE]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
    h_mYield_SE[KEY_Yield_SE]->GetXaxis()->SetTitleOffset(0.9);
    h_mYield_SE[KEY_Yield_SE]->GetXaxis()->CenterTitle();
    h_mYield_SE[KEY_Yield_SE]->GetXaxis()->SetLabelSize(0.04);
    h_mYield_SE[KEY_Yield_SE]->GetXaxis()->SetNdivisions(505);

    h_mYield_SE[KEY_Yield_SE]->GetYaxis()->SetTitle("Counts");
    h_mYield_SE[KEY_Yield_SE]->GetYaxis()->CenterTitle();
    h_mYield_SE[KEY_Yield_SE]->GetYaxis()->SetRangeUser(-0.05*h_mYield_SE[KEY_Yield_SE]->GetMaximum(),1.2*h_mYield_SE[KEY_Yield_SE]->GetMaximum());
    h_mYield_SE[KEY_Yield_SE]->GetYaxis()->SetLabelSize(0.04);
    h_mYield_SE[KEY_Yield_SE]->GetYaxis()->SetNdivisions(505);
    h_mYield_SE[KEY_Yield_SE]->SetMarkerStyle(24);
    h_mYield_SE[KEY_Yield_SE]->SetMarkerColor(kGray+2);
    h_mYield_SE[KEY_Yield_SE]->SetMarkerSize(1.0);
    h_mYield_SE[KEY_Yield_SE]->Draw("pE");

    string KEY_Yield_ME = Form("Yields_Centrality_%d_EtaGap_%d_%s_ME_SysErrors_%d",i_pad,vmsa::Eta_QA,vmsa::mPID[pid].c_str(),vmsa::Sys_QA);
    h_mYield_ME[KEY_Yield_ME]->SetFillColor(2);
    h_mYield_ME[KEY_Yield_ME]->SetFillStyle(3003);
    h_mYield_ME[KEY_Yield_ME]->Draw("h same");

    string KEY_Yield = Form("Yields_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_pad,vmsa::Eta_QA,vmsa::mPID[pid].c_str(),vmsa::Sys_QA);
    h_mYield[KEY_Yield]->SetFillColor(4);
    h_mYield[KEY_Yield]->SetFillStyle(3003);
    h_mYield[KEY_Yield]->DrawCopy("h same");
    plotTopLegend((char*)vmsa::Centrality[i_pad].c_str(),0.65,0.80,0.06,1,0.0,42,1);
  }
  */

  // Poly + Breit Wignar fit to Yields
  TH1FMap h_mYield_SM; // for QA plot only
  vecFMap ParYield_SM;

  for(int i_cent = 0; i_cent < 9; ++i_cent) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; ++i_eta) // EtaGap loop
    {
      for(int i_sys = vmsa::Sys_start; i_sys < vmsa::Sys_stop; ++i_sys) // Systematic loop
      {
	string KEY_Yield = Form("Yields_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,i_eta,vmsa::mPID[pid].c_str(),i_sys);
	h_mYield_SM[KEY_Yield] = (TH1F*)h_mYield[KEY_Yield]->Clone();
	TF1 *f_bw = new TF1("f_bw",PolyBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
	for(int i_par = 0; i_par < 5; i_par++)
	{
	  f_bw->ReleaseParameter(i_par);
	}
	f_bw->SetParameter(0,1.019);
	f_bw->SetParLimits(0,1.014,1.024);
	f_bw->SetParameter(1,0.0055);
	f_bw->SetParameter(2,10000);
	f_bw->SetParameter(3,-6000);
	f_bw->SetParameter(4,0.5);
	f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	ParYield_SM[KEY_Yield].clear();
	h_mYield[KEY_Yield]->Fit(f_bw,"NQR");
	for(int n_par = 0; n_par < 5; n_par++)
	{
	  ParYield_SM[KEY_Yield].push_back(static_cast<float>(f_bw->GetParameter(n_par)));
	}

	TF1 *f_poly = new TF1("f_poly",Poly,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],2);
	f_poly->FixParameter(0,f_bw->GetParameter(3));
	f_poly->FixParameter(1,f_bw->GetParameter(4));
	h_mYield[KEY_Yield]->Add(f_poly,-1.0); // subtract linear background for phi differential InvMass
      }
    }
  }

  /*
  //QA: Yields subtract linear background
  TCanvas *c_Yields = new TCanvas("c_Yields","c_Yields",10,10,900,900);
  c_Yields->Divide(3,3);
  for(int i_cent = 0; i_cent < 9; i_cent++)
  {
    c_Yields->cd(i_cent+1);
    c_Yields->cd(i_cent+1)->SetLeftMargin(0.15);
    c_Yields->cd(i_cent+1)->SetBottomMargin(0.15);
    c_Yields->cd(i_cent+1)->SetTicks(1,1);
    c_Yields->cd(i_cent+1)->SetGrid(0,0);

    string KEY_Yield = Form("Yields_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,vmsa::Eta_QA,vmsa::mPID[pid].c_str(),vmsa::Sys_QA);
    h_mYield[KEY_Yield]->SetTitle("");
    h_mYield[KEY_Yield]->SetStats(0);
    h_mYield[KEY_Yield]->SetMarkerStyle(24);
    h_mYield[KEY_Yield]->SetMarkerColor(kGray+3);
    h_mYield[KEY_Yield]->SetMarkerSize(0.8);
    h_mYield[KEY_Yield]->GetXaxis()->SetNdivisions(505,'N');
    h_mYield[KEY_Yield]->GetXaxis()->SetTitle("InvMass (GeV/c^{2})");
    h_mYield[KEY_Yield]->GetXaxis()->SetTitleOffset(0.9);
    h_mYield[KEY_Yield]->GetXaxis()->CenterTitle();
    h_mYield[KEY_Yield]->GetXaxis()->SetLabelSize(0.04);

    h_mYield[KEY_Yield]->GetYaxis()->SetTitle("Counts");
    h_mYield[KEY_Yield]->GetYaxis()->CenterTitle();
    h_mYield[KEY_Yield]->GetYaxis()->SetLabelSize(0.04);
    h_mYield[KEY_Yield]->DrawCopy("pE");
    PlotLine(0.98,1.05,0.0,0.0,1,2,2);

    h_mYield_SM[KEY_Yield]->SetMarkerStyle(24);
    h_mYield_SM[KEY_Yield]->SetMarkerColor(4);
    h_mYield_SM[KEY_Yield]->SetMarkerSize(0.8);
    h_mYield_SM[KEY_Yield]->DrawCopy("pE Same");

    TF1 *f_bw = new TF1("f_bw",PolyBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
    for(int i_par = 0; i_par < 5; i_par++)
    {
      f_bw->ReleaseParameter(i_par);
    }
    f_bw->FixParameter(0,ParYield_SM[KEY_Yield][0]);
    f_bw->FixParameter(1,ParYield_SM[KEY_Yield][1]);
    f_bw->FixParameter(2,ParYield_SM[KEY_Yield][2]);
    f_bw->FixParameter(3,ParYield_SM[KEY_Yield][3]);
    f_bw->FixParameter(4,ParYield_SM[KEY_Yield][4]);
    f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
    f_bw->SetLineColor(2);
    f_bw->SetLineWidth(2);
    f_bw->SetLineStyle(1);
    f_bw->DrawCopy("l same");

    TF1 *f_poly = new TF1("f_poly",Poly,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],2);
    f_poly->FixParameter(0,ParYield_SM[KEY_Yield][3]);
    f_poly->FixParameter(1,ParYield_SM[KEY_Yield][4]);
    f_poly->SetLineColor(2);
    f_poly->SetLineWidth(2);
    f_poly->SetLineStyle(2);
    f_poly->DrawCopy("l same");
    plotTopLegend((char*)vmsa::Centrality[i_cent].c_str(),0.65,0.80,0.06,1,0.0,42,1);
  }
  */

  // calculate total yields for each centrality bin via gaussian and breit wigner fits
  vecFMap ParYield_BW;
  vecFMap yields_Gaus, yields_BW;
  for(int i_cent = 0; i_cent < 9; ++i_cent) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; ++i_eta) // EtaGap loop
    {
      for(int i_sys = vmsa::Sys_start; i_sys < vmsa::Sys_stop; ++i_sys) // Systematic loop
      {
	string KEY_Yield = Form("Yields_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,i_eta,vmsa::mPID[pid].c_str(),i_sys);

	// integrating for breit wigner
	TF1 *f_yields_bw = new TF1("f_yields_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
	f_yields_bw->SetParameter(0,vmsa::InvMass[pid]);
	f_yields_bw->SetParLimits(0,vmsa::InvMass[pid]-0.001,vmsa::InvMass[pid]+0.001);
	f_yields_bw->SetParameter(1,vmsa::Width[pid]);
	f_yields_bw->SetParameter(2,1000);
	f_yields_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	h_mYield[KEY_Yield]->Fit(f_yields_bw,"MQNR");
	ParYield_BW[KEY_Yield].clear();
	ParYield_BW[KEY_Yield].push_back(static_cast<float>(f_yields_bw->GetParameter(0)));
	ParYield_BW[KEY_Yield].push_back(static_cast<float>(f_yields_bw->GetParameter(1)));
	ParYield_BW[KEY_Yield].push_back(static_cast<float>(f_yields_bw->GetParameter(2)));

	float bin_width = h_mYield[KEY_Yield]->GetBinWidth(1);
	float Inte_start = ParYield_BW[KEY_Yield][0]-vmsa::nSigVec*ParYield_BW[KEY_Yield][1]-0.5*bin_width;
	float Inte_stop  = ParYield_BW[KEY_Yield][0]+vmsa::nSigVec*ParYield_BW[KEY_Yield][1]+0.5*bin_width;
	float counts_bw = f_yields_bw->Integral(Inte_start,Inte_stop)/bin_width;
	float errors_bw = f_yields_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
	yields_BW[KEY_Yield].clear();
	yields_BW[KEY_Yield].push_back(static_cast<float>(counts_bw));
	yields_BW[KEY_Yield].push_back(static_cast<float>(errors_bw));

	// counting for guassian
	float counts_gaus = 0.0;
	float errors_gaus = 0.0;
	int bin_start = h_mYield[KEY_Yield]->FindBin(ParYield_BW[KEY_Yield][0]-vmsa::nSigVec*ParYield_BW[KEY_Yield][1]);
	int bin_stop  = h_mYield[KEY_Yield]->FindBin(ParYield_BW[KEY_Yield][0]+vmsa::nSigVec*ParYield_BW[KEY_Yield][1]);
	for(int i_bin = bin_start; i_bin <= bin_stop; ++i_bin)
	{
	  counts_gaus += h_mYield[KEY_Yield]->GetBinContent(i_bin);
	  errors_gaus += h_mYield[KEY_Yield]->GetBinError(i_bin)*h_mYield[KEY_Yield]->GetBinError(i_bin);
	}
	yields_Gaus[KEY_Yield].clear();
	yields_Gaus[KEY_Yield].push_back(static_cast<float>(counts_gaus));
	yields_Gaus[KEY_Yield].push_back(static_cast<float>(TMath::Sqrt(errors_gaus)));
      }
    }
  }

  /*
  // QA: different counting method: bin counting vs breit wigner integrating
  TCanvas *c_Yields_counts = new TCanvas("c_Yields_counts","c_Yields_counts",10,10,900,900);
  c_Yields_counts->Divide(3,3);
  for(int i_cent = 0; i_cent < 9; i_cent++)
  {
    c_Yields_counts->cd(i_cent+1);
    c_Yields_counts->cd(i_cent+1)->SetLeftMargin(0.20);
    c_Yields_counts->cd(i_cent+1)->SetBottomMargin(0.20);
    c_Yields_counts->cd(i_cent+1)->SetTicks(1,1);
    c_Yields_counts->cd(i_cent+1)->SetGrid(0,0);

    string KEY_Yield = Form("Yields_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,vmsa::Eta_QA,vmsa::mPID[pid].c_str(),vmsa::Sys_QA);
    h_mYield[KEY_Yield]->SetTitle("");
    h_mYield[KEY_Yield]->SetStats(0);
    h_mYield[KEY_Yield]->SetMarkerStyle(24);
    h_mYield[KEY_Yield]->SetMarkerColor(kGray+3);
    h_mYield[KEY_Yield]->SetMarkerSize(0.8);
    h_mYield[KEY_Yield]->GetXaxis()->SetNdivisions(505,'N');
    h_mYield[KEY_Yield]->GetXaxis()->SetTitle("InvMass (GeV/c^{2})");
    h_mYield[KEY_Yield]->GetXaxis()->CenterTitle();
    h_mYield[KEY_Yield]->GetYaxis()->SetTitle("Counts");
    h_mYield[KEY_Yield]->GetYaxis()->CenterTitle();
    h_mYield[KEY_Yield]->GetYaxis()->SetTitleOffset(1.2);
    h_mYield[KEY_Yield]->DrawCopy("pE");
    PlotLine(0.98,1.05,0.0,0.0,1,2,2);

    TF1 *f_yields_bw = new TF1("f_yields_bw",BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
    f_yields_bw->SetParameter(0,ParYield_BW[KEY_Yield][0]);
    f_yields_bw->SetParameter(1,ParYield_BW[KEY_Yield][1]);
    f_yields_bw->SetParameter(2,ParYield_BW[KEY_Yield][2]);
    f_yields_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
    f_yields_bw->SetLineColor(2);
    f_yields_bw->SetLineStyle(1);
    f_yields_bw->SetLineWidth(2);
    f_yields_bw->DrawCopy("l same");

    float x1 = ParYield_BW[KEY_Yield][0] - vmsa::nSigVec*ParYield_BW[KEY_Yield][1];
    float x2 = ParYield_BW[KEY_Yield][0] + vmsa::nSigVec*ParYield_BW[KEY_Yield][1];
    float y = 0.7*h_mYield[KEY_Yield]->GetBinContent(h_mYield[KEY_Yield]->FindBin(ParYield_BW[KEY_Yield][0]));
    PlotLine(x1,x1,0,y,4,2,2);
    PlotLine(x2,x2,0,y,4,2,2);
    plotTopLegend((char*)vmsa::Centrality[i_cent].c_str(),0.65,0.80,0.06,1,0.0,42,1);
  }
  */

  // calculate final resolution correction factors and correct flow
  string InPutFile_Res = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Resolution/file_%s_Resolution.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Res = TFile::Open(InPutFile_Res.c_str());
  TProMap p_mRes;
  vecFMap ResValue;
  TH1FMap h_mRes;
  for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; ++i_eta) // eta_gap loop
  {
    string KEY_eta = Form("Res2_EtaGap_%d_EP",i_eta);
    p_mRes[KEY_eta] = (TProfile*)File_Res->Get(KEY_eta.c_str()); // read in resolution file
    for(int i_sys = vmsa::Sys_start; i_sys < vmsa::Sys_stop; i_sys++) // Systematic errors loop
    {
      for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin loop
      {
	string KEY_Res = Form("h_mRes_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,i_eta,vmsa::mPID[pid].c_str(),i_sys);
	h_mRes[KEY_Res] = new TH1F(KEY_Res.c_str(),KEY_Res.c_str(),2,-0.5,1.5);
	float yields_total_gaus = 0.0;
	float yields_total_bw   = 0.0;
	for(int cent = vmsa::cent_low[i_cent]; cent <= vmsa::cent_up[i_cent]; cent++) // calculate resolution and total yields in selected centrality bin
	{
	  if(p_mRes[KEY_eta]->GetBinContent(cent+1) > 0) 
	  {
	    string KEY_Yield = Form("Yields_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",cent,i_eta,vmsa::mPID[pid].c_str(),i_sys);
	    ResValue[KEY_Yield].push_back(static_cast<float>(TMath::Sqrt(p_mRes[KEY_eta]->GetBinContent(cent+1))));
	    yields_total_gaus += yields_Gaus[KEY_Yield][0];
	    yields_total_bw += yields_BW[KEY_Yield][0];
	  }
	}

	string KEY_ResCorr = Form("Res_2nd_Centrality_%d_EtaGap_%d_SysErrors_%d",i_cent,i_eta,i_sys); // KEY for final resolution correction factor
	float mean_res_gaus = 0.0;
	float mean_res_bw = 0.0;
	for(int cent = vmsa::cent_low[i_cent]; cent <= vmsa::cent_up[i_cent]; cent++) // calculate final resolution correction factor <R(centrality)>
	{
	  string KEY_Yield = Form("Yields_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",cent,i_eta,vmsa::mPID[pid].c_str(),i_sys);
	  mean_res_gaus += ResValue[KEY_Yield][0]*yields_Gaus[KEY_Yield][0]/yields_total_gaus;
	  mean_res_bw   += ResValue[KEY_Yield][0]*yields_BW[KEY_Yield][0]/yields_total_bw;
	}
	cout << "i_eta = " << i_eta << ", i_sys = " << i_sys << ", centrality_bin = " << i_cent << ", mean_res_gaus = " << mean_res_gaus << endl;
	cout << "i_eta = " << i_eta << ", i_sys = " << i_sys << ", centrality_bin = " << i_cent << ", mean_res_bw = " << mean_res_bw << endl;
	h_mRes[KEY_Res]->SetBinContent(1,mean_res_gaus);
	h_mRes[KEY_Res]->SetBinContent(2,mean_res_bw);
      }
    }
  }

  string OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/MonteCarlo/Data/Resolution.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; ++i_eta) // eta_gap loop
  {
    for(int i_sys = vmsa::Sys_start; i_sys < vmsa::Sys_stop; i_sys++) // Systematic errors loop
    {
      for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin loop
      {
	string KEY_Res = Form("h_mRes_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,i_eta,vmsa::mPID[pid].c_str(),i_sys);
	h_mRes[KEY_Res]->Write();
      }
    }
  }
  File_OutPut->Close();
}

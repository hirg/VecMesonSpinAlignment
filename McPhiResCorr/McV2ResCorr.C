#include <string>
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/functions.h"

#ifndef _PlotQA_
#define _PlotQA_ 0
#endif

using namespace std;
float const Resolution = 0.20;

double getChi(double *var, double *par)
{
  double chi = var[0];
  double arg = chi*chi/4.0;
  double norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  double y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return y;
}

double EventPlane(double *var, double *par)
{
  double DeltaPsi = var[0];
  double chi = par[0];
  double arg = chi/TMath::Sqrt(2.0);
  double arg2 = -0.5*chi*chi;
  double pi = TMath::Pi();
  double norm = 0.5/pi;

  double cos = TMath::Cos(2.0*DeltaPsi);
  double sin2 = TMath::Sin(2.0*DeltaPsi)*TMath::Sin(2.0*DeltaPsi);
  double y = norm*(TMath::Exp(arg2)+TMath::Sqrt(pi)*arg*cos*TMath::Exp(arg2*sin2)*(1.0+TMath::Erf(arg*cos)));

  return y;
}

// histograms
TH3F *h_Tracks;
TH2F *h_phiRP, *h_phiEP;
TH1F *h_Psi2;
TProfile *p_cosRP, *p_sinRP;
TProfile *p_cosEP, *p_sinEP;

// sampling functions
TF1 *f_v2, *f_spec, *f_flow, *f_gaus;
TF1 *f_com, f_res;

TF1* readv2(int energy, int pid, int centrality)
{
  string InPutV2 = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Phi_v2_1040.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_v2 = TFile::Open(InPutV2.c_str());
  TGraphAsymmErrors *g_v2 = (TGraphAsymmErrors*)File_v2->Get("g_v2");
  TF1 *f_v2 = new TF1("f_v2",v2_pT_FitFunc,vmsa::ptMin,vmsa::ptMax,5);
  f_v2->FixParameter(0,2);
  f_v2->SetParameter(1,0.1);
  f_v2->SetParameter(2,0.1);
  f_v2->SetParameter(3,0.1);
  f_v2->SetParameter(4,0.1);
  f_v2->SetLineColor(2);
  f_v2->SetLineWidth(2);
  f_v2->SetLineStyle(2);
  g_v2->Fit(f_v2,"N");

#if _PlotQA_
  TCanvas *c_v2 = new TCanvas("c_v2","c_v2",10,10,800,800);
  c_v2->cd()->SetLeftMargin(0.15);
  c_v2->cd()->SetBottomMargin(0.15);
  c_v2->cd()->SetTicks(1,1);
  c_v2->cd()->SetGrid(0,0);
  TH1F *h_v2 = new TH1F("h_v2","h_v2",100,0.0,10.0);
  for(int i_bin = 1; i_bin < 101; ++i_bin)
  {
    h_v2->SetBinContent(i_bin,-10.0);
    h_v2->SetBinError(i_bin,1.0);
  }
  h_v2->SetTitle("");
  h_v2->SetStats(0);
  h_v2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_v2->GetXaxis()->CenterTitle();
  h_v2->GetYaxis()->SetTitle("v_{2}");
  h_v2->GetYaxis()->CenterTitle();
  h_v2->GetYaxis()->SetRangeUser(0.0,0.2);
  h_v2->Draw("pE");
  g_v2->Draw("pE same");
  f_v2->Draw("l same");
#endif

  return f_v2;
}

TF1* readspec(int energy, int pid, int centrality)
{
  string InPutSpec = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Phi_Spec.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Spec = TFile::Open(InPutSpec.c_str());
  TGraphAsymmErrors *g_spec = (TGraphAsymmErrors*)File_Spec->Get("g_spec");
  TF1 *f_Levy = new TF1("f_Levy",Levy,vmsa::ptMin,vmsa::ptMax,3);
  f_Levy->SetParameter(0,1);
  f_Levy->SetParameter(1,10);
  f_Levy->SetParameter(2,0.1);
  f_Levy->SetLineStyle(2);
  f_Levy->SetLineColor(4);
  f_Levy->SetLineWidth(2);
  g_spec->Fit(f_Levy,"N");

  TF1 *f_spec = new TF1("f_spec",pTLevy,vmsa::ptMin,vmsa::ptMax,3);
  f_spec->SetParameter(0,f_Levy->GetParameter(0));
  f_spec->SetParameter(1,f_Levy->GetParameter(1));
  f_spec->SetParameter(2,f_Levy->GetParameter(2));
  f_spec->SetLineStyle(2);
  f_spec->SetLineColor(2);
  f_spec->SetLineWidth(2);

#if _PlotQA_
  TCanvas *c_spec = new TCanvas("c_spec","c_spec",10,10,800,800);
  c_spec->cd()->SetLeftMargin(0.15);
  c_spec->cd()->SetBottomMargin(0.15);
  c_spec->cd()->SetTicks(1,1);
  c_spec->cd()->SetGrid(0,0);
  c_spec->SetLogy();
  TH1F *h_spec = new TH1F("h_spec","h_spec",100,0.0,10.0);
  for(int i_bin = 1; i_bin < 101; ++i_bin)
  {
    h_spec->SetBinContent(i_bin,-10.0);
    h_spec->SetBinError(i_bin,1.0);
  }
  h_spec->SetTitle("");
  h_spec->SetStats(0);
  h_spec->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_spec->GetXaxis()->CenterTitle();
  h_spec->GetYaxis()->SetTitle("dN/p_{T}dp_{T}");
  h_spec->GetYaxis()->CenterTitle();
  h_spec->GetYaxis()->SetRangeUser(1E-6,10);
  h_spec->Draw("pE");
  g_spec->Draw("pE same");
  f_Levy->Draw("l same");
  f_spec->Draw("l same");
#endif

  return f_spec;
}

void getKinematics(TLorentzVector& lPhi, double const mass)
{
  f_flow->ReleaseParameter(0);
  // double const pt = f_spec->GetRandom(vmsa::ptMin, vmsa::ptMax);
  double const pt = gRandom->Uniform(vmsa::ptMin, vmsa::ptMax);
  double const y = gRandom->Uniform(-vmsa::acceptanceRapidity, vmsa::acceptanceRapidity);
  f_flow->SetParameter(0,f_v2->Eval(pt));
  double const phi = f_flow->GetRandom();

  double const mT = sqrt(mass * mass + pt * pt);
  double const pz = mT * sinh(y);
  double const E = mT * cosh(y);

  lPhi.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
}

float GausSmearing(TF1 *f_gaus)
{
  // float Psi2 = f_gaus->GetRandom(-TMath::TwoPi(),TMath::TwoPi());
  float Psi2 = f_gaus->GetRandom();
  return Psi2;
}

void fill(TLorentzVector* lPhi)
{
  h_phiRP->Fill(lPhi->Pt(),lPhi->Phi());

  float cosRP = cos(2.0*lPhi->Phi());
  p_cosRP->Fill(lPhi->Pt(),cosRP);
  float sinRP = sin(2.0*lPhi->Phi());
  p_sinRP->Fill(lPhi->Pt(),sinRP);

  float Psi2 = GausSmearing(f_gaus);
  h_Psi2->Fill(Psi2);

  float phiSmear = lPhi->Phi()-Psi2;
  if(phiSmear > TMath::Pi()) phiSmear -= TMath::TwoPi();
  if(phiSmear < -TMath::Pi()) phiSmear += TMath::TwoPi();
  h_phiEP->Fill(lPhi->Pt(),phiSmear);
  float cosEP = cos(2.0*lPhi->Phi())*cos(2.0*Psi2);
  p_cosEP->Fill(lPhi->Pt(),cosEP);
  float sinEP = sin(2.0*lPhi->Phi())*sin(2.0*Psi2);
  p_sinEP->Fill(lPhi->Pt(),sinEP);
}

void write(int energy)
{
  string OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/McV2.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_Tracks->Write();
  h_phiRP->Write();
  h_phiEP->Write();
  h_Psi2->Write();
  p_cosRP->Write();
  p_cosEP->Write();
  p_sinRP->Write();
  p_sinEP->Write();
  File_OutPut->Close();
}

void McV2ResCorr(int energy = 6, int pid = 0, int cent = 0, int NMax = 100000)
{
  int   const BinPt    = vmsa::BinPt;
  int   const BinY     = vmsa::BinY;
  int   const BinPhi   = vmsa::BinPhi;

  f_v2 = readv2(energy,pid,cent);
  f_spec = readspec(energy,pid,cent);
  f_flow = new TF1("f_flow",flowSample,-TMath::Pi(),TMath::Pi(),1);
  f_gaus = new TF1("f_gaus",GausSmearing,-TMath::TwoPi(),TMath::TwoPi(),1);
  f_gaus->FixParameter(0,Resolution);

  string HistName;
  h_Tracks = new TH3F("h_Tracks","h_Tracks",BinPt,vmsa::ptMin,vmsa::ptMax,BinY,-1.0,1.0,BinPhi,-TMath::Pi(),TMath::Pi());
  h_phiRP = new TH2F("h_phiRP","h_phiRP",BinPt,vmsa::ptMin,vmsa::ptMax,BinPhi,-TMath::Pi(),TMath::Pi());
  h_phiEP = new TH2F("h_phiEP","h_phiEP",BinPt,vmsa::ptMin,vmsa::ptMax,BinPhi,-TMath::Pi(),TMath::Pi());
  h_Psi2 = new TH1F("h_Psi2","h_Psi2",BinPhi*10,-TMath::TwoPi(),TMath::TwoPi());
  
  p_cosRP = new TProfile("p_cosRP","p_cosRP",BinPt,vmsa::ptMin,vmsa::ptMax);
  p_sinRP = new TProfile("p_sinRP","p_sinRP",BinPt,vmsa::ptMin,vmsa::ptMax);
  p_cosEP = new TProfile("p_cosEP","p_cosEP",BinPt,vmsa::ptMin,vmsa::ptMax);
  p_sinEP = new TProfile("p_sinEP","p_sinEP",BinPt,vmsa::ptMin,vmsa::ptMax);

  TStopwatch* stopWatch = new TStopwatch();
  stopWatch->Start();
  if(gRandom) delete gRandom;
  gRandom = new TRandom3();
  gRandom->SetSeed();

  TLorentzVector *lPhi = new TLorentzVector();
  for(int i_ran = 0; i_ran < NMax; ++i_ran)
  {
    if (floor(10.0*i_ran/ static_cast<float>(NMax)) > floor(10.0*(i_ran-1)/ static_cast<float>(NMax)))
    cout << "=> processing data: " << 100.0*i_ran/ static_cast<float>(NMax) << "%" << endl;

    getKinematics(*lPhi,vmsa::InvMass[pid]);
    fill(lPhi);
  }
  cout << "=> processing data: 100%" << endl;
  cout << "work done!" << endl;

  write(energy);

  stopWatch->Stop();   
  stopWatch->Print();
}

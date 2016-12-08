#include <iostream>
#include <string> 
#include <map>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TPythia6.h"
#include "TPythia6Decayer.h"
#include "TParticle.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TVector3.h"
#include "/global/homes/x/xusun/STAR/VecMesonSpinAlignment/Utility/functions.h"
#include "/global/homes/x/xusun/STAR/VecMesonSpinAlignment/Utility/StSpinAlignmentCons.h"

using namespace std;

float readRes(int energy, int pid, int centrality);
TF1* readv2(int energy, int pid, int centrality);
TF1* readspec(int energy, int pid, int centrality);
void getKinematics(TLorentzVector& lPhi, double const mass);
void setDecayChannels(int const pid);
void decayAndFill(int const kf, TLorentzVector* lPhi, TClonesArray& daughters);
void fill(TLorentzVector* lPhi, TLorentzVector const& lKplus, TLorentzVector const& lKminus);
void write(int energy,int Nrho);
TVector3 CalBoostedVector(TLorentzVector const lMcDau, TLorentzVector *lMcVec);
bool Sampling(TF1 *f_rhoPhy,float CosThetaStar);
float GausSmearing(TF1 *f_gaus);
float getChi(float Resolution);
float EventPlaneSmearing(TF1 *f_gaus);

// histograms
TH3F *h_Tracks;
TH2F *h_phiQA, *h_phiRP, *h_phiGaus, *h_phiEP;
TH2F *h_cosQA, *h_cosRP, *h_cosGaus, *h_cosEP;
TH1F *h_PsiGaus, *h_PsiEP;

// sampling functions
TF1 *f_v2, *f_spec, *f_flow, *f_rhoPhy, *f_gaus, *f_EP;

TPythia6Decayer* pydecay;

void McPhiResCorr(int energy = 6, int pid = 0, int cent = 0, int Nrho = 40, int const NMax = 100000)
{
  int   const BinPt    = vmsa::BinPt;
  int   const BinY     = vmsa::BinY;
  int   const BinPhi   = vmsa::BinPhi;
  float const rhoDelta = vmsa::rhoDelta;

  string Info = Form("sampling rhophy = %.2f with %d tracks!!!!",rhoDelta*Nrho,NMax);
  cout << Info.c_str() << endl;

  string HistName;
  HistName = Form("h_Tracks_%d",Nrho);
  h_Tracks = new TH3F(HistName.c_str(),HistName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax,BinY,-1.0,1.0,BinPhi,-TMath::Pi(),TMath::Pi());

  HistName = Form("h_phiQA_%d",Nrho);
  h_phiQA = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax,BinPhi,-TMath::Pi(),TMath::Pi());
  HistName = Form("h_phiRP_%d",Nrho);
  h_phiRP = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax,BinPhi,-TMath::Pi(),TMath::Pi());
  HistName = Form("h_phiGaus_%d",Nrho);
  h_phiGaus = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax,BinPhi,-TMath::Pi(),TMath::Pi());
  HistName = Form("h_phiEP_%d",Nrho);
  h_phiEP = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax,BinPhi,-TMath::Pi(),TMath::Pi());

  HistName = Form("h_cosQA_%d",Nrho);
  h_cosQA = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax,BinY,-1.0,1.0);
  HistName = Form("h_cosRP_%d",Nrho);
  h_cosRP = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax,BinY,-1.0,1.0);
  HistName = Form("h_cosGaus_%d",Nrho);
  h_cosGaus = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax,BinY,-1.0,1.0);
  HistName = Form("h_cosEP_%d",Nrho);
  h_cosEP = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax,BinY,-1.0,1.0);

  HistName = Form("h_PsiGaus_%d",Nrho);
  h_PsiGaus  = new TH1F(HistName.c_str(),HistName.c_str(),BinPhi*10,-TMath::PiOver2(),TMath::PiOver2());
  HistName = Form("h_PsiEP_%d",Nrho);
  h_PsiEP = new TH1F(HistName.c_str(),HistName.c_str(),BinPhi*10,-TMath::PiOver2(),TMath::PiOver2());

  f_v2   = readv2(energy,pid,cent);
  f_spec = readspec(energy,pid,cent);

  f_flow = new TF1("f_flow",flowSample,-TMath::Pi(),TMath::Pi(),1);

  float rhoPhy = Nrho*rhoDelta;
  f_rhoPhy = new TF1("f_rhoPhy",SpinDensity,-1.0,1.0,2);
  f_rhoPhy->FixParameter(0,rhoPhy);
  f_rhoPhy->FixParameter(1,0.75);

  float const Resolution = 0.81; // xin's resolution
  // float const Resolution = readRes(energy,pid,cent);
  cout << "InPut Resolution = " << Resolution << endl;
  f_gaus = new TF1("f_gaus",EventPlaneGaus,-TMath::PiOver2(),TMath::PiOver2(),1);
  f_gaus->FixParameter(0,Resolution);

  float const chi = getChi(Resolution);
  f_EP = new TF1("f_EP",EventPlaneDist,-TMath::PiOver2(),TMath::PiOver2(),1);
  f_EP->FixParameter(0,chi);

  TStopwatch* stopWatch = new TStopwatch();
  stopWatch->Start();
  if(gRandom) delete gRandom;
  gRandom = new TRandom3();
  gRandom->SetSeed();

  pydecay = TPythia6Decayer::Instance();
  pydecay->Init();
  setDecayChannels(pid); // phi--> K+K-

  TClonesArray ptl("TParticle", 10);
  TLorentzVector *lPhi = new TLorentzVector();
  for(int i_ran = 0; i_ran < NMax; ++i_ran)
  {
    if (floor(10.0*i_ran/ static_cast<float>(NMax)) > floor(10.0*(i_ran-1)/ static_cast<float>(NMax)))
    cout << "=> processing data: " << 100.0*i_ran/ static_cast<float>(NMax) << "%" << endl;

    getKinematics(*lPhi,vmsa::InvMass[pid]);
    decayAndFill(vmsa::decayMother[pid],lPhi,ptl);
  }
  cout << "=> processing data: 100%" << endl;
  cout << "work done!" << endl;

  write(energy,Nrho);

  stopWatch->Stop();   
  stopWatch->Print();
}

float readRes(int energy, int pid, int centrality)
{
  string InPutRes = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Resolution.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Res = TFile::Open(InPutRes.c_str());
  string HistName = Form("h_mRes_Centrality_%d_EtaGap_0_Phi_SysErrors_0",centrality);
  TH1F *h_mRes = (TH1F*)File_Res->Get(HistName.c_str());
  float resGaus = h_mRes->GetBinContent(1);
  float resBW   = h_mRes->GetBinContent(2);
  cout << "resGaus = " << resGaus << ", resBW = " << resBW << endl;

  return 0.5*(resGaus+resBW);
}

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

  /*
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
  */

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

  /*
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
  */

  return f_spec;
}

void getKinematics(TLorentzVector& lPhi, double const mass)
{
  f_flow->ReleaseParameter(0);
  double const pt = f_spec->GetRandom(vmsa::ptMin, vmsa::ptMax);
  double const y = gRandom->Uniform(-vmsa::acceptanceRapidity, vmsa::acceptanceRapidity);
  f_flow->SetParameter(0,f_v2->Eval(pt));
  double const phi = f_flow->GetRandom();

  // double const pt = gRandom->Uniform(vmsa::ptMin, vmsa::ptMax);
  // double const y = gRandom->Uniform(-vmsa::acceptanceRapidity, vmsa::acceptanceRapidity);
  // double const phi = TMath::TwoPi() * gRandom->Rndm();

  double const mT = sqrt(mass * mass + pt * pt);
  double const pz = mT * sinh(y);
  double const E = mT * cosh(y);

  lPhi.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
}

void setDecayChannels(int const pid)
{
  int const mdme = vmsa::decayChannels[pid];
  cout << "mdme = " << mdme << endl;
  for (int idc = vmsa::decayChannelsFirst[pid]; idc < vmsa::decayChannelsSecond[pid] + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0); // close all decay channel
  TPythia6::Instance()->SetMDME(mdme, 1, 1); // open the one we need
  int *PYSeed = new int;
  TPythia6::Instance()->SetMRPY(1,(int)PYSeed); // Random seed
}

void decayAndFill(int const kf, TLorentzVector* lPhi, TClonesArray& daughters)
{
  pydecay->Decay(kf, lPhi);
  pydecay->ImportParticles(&daughters);

  TLorentzVector lKplus;
  TLorentzVector lKminus;

  int nTrk = daughters.GetEntriesFast();
  for (int iTrk = 0; iTrk < nTrk; ++iTrk)
  {
    TParticle* ptl0 = (TParticle*)daughters.At(iTrk);

    switch (ptl0->GetPdgCode())
    {
      case 321:
	ptl0->Momentum(lKplus);
	break;
      case -321:
	ptl0->Momentum(lKminus);
	break;
      default:
	break;
    }
  }
  daughters.Clear("C");

  fill(lPhi,lKplus,lKminus);
}

void fill(TLorentzVector* lPhi, TLorentzVector const& lKplus, TLorentzVector const& lKminus)
{
  TVector3 vMcKpBoosted = CalBoostedVector(lKplus,lPhi); // boost Kplus back to phi-meson rest frame

  TVector3 nQ(0.0,-1.0,0.0); // direction of angular momentum with un-smeared EP
  float CosThetaStarRP = vMcKpBoosted.Dot(nQ);
  h_phiQA->Fill(lPhi->Pt(),lPhi->Phi());
  h_cosQA->Fill(lPhi->Pt(),CosThetaStarRP);
  h_Tracks->Fill(lPhi->Pt(),lPhi->Rapidity(),lPhi->Phi());

  if(!Sampling(f_rhoPhy,CosThetaStarRP)) return;
  h_phiRP->Fill(lPhi->Pt(),lPhi->Phi());
  h_cosRP->Fill(lPhi->Pt(),CosThetaStarRP);

  float PsiGaus = GausSmearing(f_gaus);
  TVector3 nQGausSmear(sin(PsiGaus),-cos(PsiGaus),0.0); // direction of angular momentum with Gaussian Smearing
  float CosThetaStarGaus = vMcKpBoosted.Dot(nQGausSmear);
  float phiGausSmear = lPhi->Phi()-PsiGaus;
  if(phiGausSmear > TMath::Pi())  phiGausSmear -= TMath::TwoPi();
  if(phiGausSmear < -TMath::Pi()) phiGausSmear += TMath::TwoPi();
  h_phiGaus->Fill(lPhi->Pt(),phiGausSmear);
  h_cosGaus->Fill(lPhi->Pt(),CosThetaStarGaus);
  h_PsiGaus->Fill(PsiGaus);

  float PsiEP = EventPlaneSmearing(f_EP);
  TVector3 nQSmear(sin(PsiEP),-cos(PsiEP),0.0); // direction of angular momentum with EP Smearing
  float CosThetaStarEP = vMcKpBoosted.Dot(nQSmear);
  float phiSmear = lPhi->Phi()-PsiEP;
  if(phiSmear > TMath::Pi())  phiSmear -= TMath::TwoPi();
  if(phiSmear < -TMath::Pi()) phiSmear += TMath::TwoPi();
  h_phiEP->Fill(lPhi->Pt(),phiSmear);
  h_cosEP->Fill(lPhi->Pt(),CosThetaStarEP);
  h_PsiEP->Fill(PsiEP);
}

float GausSmearing(TF1 *f_gaus)
{
  float Psi2 = f_gaus->GetRandom(-TMath::PiOver2(),TMath::PiOver2());
  return Psi2;
}

float getChi(float Resolution)
{
  TF1 *f_res = new TF1("f_res",EventPlaneResolution,0,10,0);
  double chi = f_res->GetX(Resolution);

  return chi;
}

float EventPlaneSmearing(TF1 *f_EP)
{
  float Psi2 = f_EP->GetRandom(-TMath::PiOver2(),TMath::PiOver2());
  return Psi2;
}

TVector3 CalBoostedVector(TLorentzVector const lMcDau, TLorentzVector *lMcVec)
{
  TVector3 vMcBeta = -1.0*lMcVec->BoostVector(); // boost vector

  TLorentzVector lKaon = lMcDau;
  lKaon.Boost(vMcBeta); // boost Kplus back to phi-meson rest frame
  TVector3 vMcDauStar = lKaon.Vect().Unit(); // momentum direction of Kplus in phi-meson rest frame

  return vMcDauStar;
}

bool Sampling(TF1 *f_rhoPhy, float CosThetaStar)
{
  float wMax;
  if(f_rhoPhy->GetParameter(0) <= 1.0/3.0) wMax = f_rhoPhy->Eval(0.0);
  else wMax = f_rhoPhy->Eval(1.0);
  return !(gRandom->Rndm() > f_rhoPhy->Eval(CosThetaStar)/wMax);
}

void write(int energy, int Nrho)
{
  string OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/McPhiResCorr_%d.root",vmsa::mBeamEnergy[energy].c_str(),Nrho);
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_Tracks->Write();
  h_phiQA->Write();
  h_phiRP->Write();
  h_phiGaus->Write();
  h_phiEP->Write();
  h_cosQA->Write();
  h_cosRP->Write();
  h_cosGaus->Write();
  h_cosEP->Write();
  h_PsiGaus->Write();
  h_PsiEP->Write();
  File_OutPut->Close();
}


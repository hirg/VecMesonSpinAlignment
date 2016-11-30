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

using namespace std;

std::pair<int, int> const decayChannels(656,666); // phi decay channel
std::pair<double, double> const momentumRange(0.2,5.0);
double const acceptanceRapidity = 1.0;
double const invMass = 1.01940;
std::string const  mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
int const BinPt = 20;
int const BinY = 20;
int const BinPhi = 36;
// float const Resolution = 0.81;
float const Resolution = 0.999;
float const rhoDelta = 0.01;

double flow(double *x_val, double *par)
{
  double x, y, v2;
  x = x_val[0];
  v2 = par[0];
  y = 1.0 + 2.0*v2*cos(2.0*x);

  return y;
}

double v2_pT_FitFunc(double *x_val, double *par)
{
  // Fit function for v2 vs. pT
  // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
  double v2, pT, a, b, c, d, n;
  pT = x_val[0];
  n  = par[0]; // number-of-constituent quarks
  a  = par[1];
  b  = par[2];
  c  = par[3];
  d  = par[4];

  if(c != 0.0)
  {
    v2 = a*n/(1.0 + exp(-(pT/n - b)/c)) - d*n;
  }
  else v2 = 0.0;

  return v2;
}

double Levy(double *var, double *par)
{
  double const m0 = 1.01940; // phi-meson mass
  double pT   = var[0];
  double mT   = sqrt(pT*pT+m0*m0);
  double dNdy = par[0];
  double n    = par[1];
  double T    = par[2];

  double numer = dNdy*(n-1)*(n-2);
  double denom = n*T*(n*T+m0*(n-2));
  double power = pow(1+(mT-m0)/(n*T),-1.0*n);

  double y = numer*power/denom;

  return y;
}


double pTLevy(double *var, double *par)
{
  double const m0 = 1.01940; // phi-meson mass
  double pT   = var[0];
  double mT   = sqrt(pT*pT+m0*m0);
  double dNdy = par[0];
  double n    = par[1];
  double T    = par[2];

  double numer = dNdy*(n-1)*(n-2);
  double denom = n*T*(n*T+m0*(n-2));
  double power = pow(1+(mT-m0)/(n*T),-1.0*n);

  double y = pT*numer*power/denom;

  return y;
}

double SpinDist(double *var, double *par)
{
  double CosThetaStar = var[0];
  double rho = par[0];
  double y = 0.75*((1.0-rho)+(3.0*rho-1)*CosThetaStar*CosThetaStar);

  return y;
}

double GausSmearing(double *var, double *par)
{
  double Psi2 = var[0];
  double res = par[0]; // res = <cos(2*(Psi2-Psi_RP))>

  double sigma = acos(res)/2.0;
  double sigmaSquare = sigma*sigma;
  double norm = 1.0/sqrt(2.0*sigmaSquare*TMath::Pi());
  double power = -1.0*Psi2*Psi2/(2.0*sigmaSquare);

  double y = norm*exp(power);

  return y;
}

TF1* readv2(int energy);
TF1* readspec(int energy);
void getKinematics(TLorentzVector& lPhi, double const mass);
void setDecayChannels(int const mdme);
void decayAndFill(int const kf, TLorentzVector* lPhi, TClonesArray& daughters);
void fill(TLorentzVector* lPhi, TLorentzVector const& lKplus, TLorentzVector const& lKminus);
void write(int energy,int Nrho);
TVector3 CalBoostedVector(TLorentzVector const lMcDau, TLorentzVector *lMcVec);
bool Sampling(TF1 *f_rhoPhy,float CosThetaStar);
float EventPlaneSmearing(TF1 *f_gaus);

// histograms
TH3F *h_Tracks;
TH2F *h_phiQA, *h_phiRP, *h_phiEP;
TH2F *h_cosQA, *h_cosRP, *h_cosEP;
TH1F *h_Psi2;

// sampling functions
TF1 *f_v2, *f_spec, *f_flow, *f_rhoPhy, *f_gaus;

TPythia6Decayer* pydecay;

void McPhiResCorr(int energy = 6, int Nrho = 40, int const NMax = 1000000)
{
  string Info = Form("sampling rhophy = %.2f with %d tracks!!!!",0.01*Nrho,NMax);
  cout << Info.c_str() << endl;
  string HistName;
  HistName = Form("h_Tracks_%d",Nrho);
  h_Tracks = new TH3F(HistName.c_str(),HistName.c_str(),BinPt,momentumRange.first,momentumRange.second,BinY,-1.0,1.0,BinPhi,-TMath::Pi(),TMath::Pi());
  HistName = Form("h_phiQA_%d",Nrho);
  h_phiQA = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,momentumRange.first,momentumRange.second,BinPhi,-TMath::Pi(),TMath::Pi());
  HistName = Form("h_phiRP_%d",Nrho);
  h_phiRP = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,momentumRange.first,momentumRange.second,BinPhi,-TMath::Pi(),TMath::Pi());
  HistName = Form("h_phiEP_%d",Nrho);
  h_phiEP = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,momentumRange.first,momentumRange.second,BinPhi,-TMath::Pi(),TMath::Pi());
  HistName = Form("h_cosQA_%d",Nrho);
  h_cosQA = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,momentumRange.first,momentumRange.second,BinY,-1.0,1.0);
  HistName = Form("h_cosRP_%d",Nrho);
  h_cosRP = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,momentumRange.first,momentumRange.second,BinY,-1.0,1.0);
  HistName = Form("h_cosEP_%d",Nrho);
  h_cosEP = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,momentumRange.first,momentumRange.second,BinY,-1.0,1.0);
  HistName = Form("h_Psi2_%d",Nrho);
  h_Psi2 = new TH1F(HistName.c_str(),HistName.c_str(),BinPhi*10,-TMath::PiOver2(),TMath::PiOver2());

  f_v2   = readv2(energy);
  f_spec = readspec(energy);

  f_flow = new TF1("f_flow",flow,-TMath::Pi(),TMath::Pi(),1);
  float rhoPhy = Nrho*rhoDelta;
  f_rhoPhy = new TF1("f_rhoPhy",SpinDist,-1.0,1.0,1);
  f_rhoPhy->FixParameter(0,rhoPhy);
  f_gaus = new TF1("f_gaus",GausSmearing,-TMath::PiOver2(),TMath::PiOver2(),1);
  f_gaus->FixParameter(0,Resolution);

  TStopwatch* stopWatch = new TStopwatch();
  stopWatch->Start();
  if(gRandom) delete gRandom;
  gRandom = new TRandom3();
  gRandom->SetSeed();

  pydecay = TPythia6Decayer::Instance();
  pydecay->Init();
  setDecayChannels(656); // phi--> K+K-

  TClonesArray ptl("TParticle", 10);
  TLorentzVector *lPhi = new TLorentzVector();
  for(int i_ran = 0; i_ran < NMax; ++i_ran)
  {
    if (floor(10.0*i_ran/ static_cast<float>(NMax)) > floor(10.0*(i_ran-1)/ static_cast<float>(NMax)))
    cout << "=> processing data: " << 100.0*i_ran/ static_cast<float>(NMax) << "%" << endl;

    getKinematics(*lPhi,invMass);
    decayAndFill(333,lPhi,ptl);
  }
  cout << "=> processing data: 100%" << endl;
  cout << "work done!" << endl;

  write(energy,Nrho);

  stopWatch->Stop();   
  stopWatch->Print();
}

TF1* readv2(int energy)
{
  string InPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Phi_v2_1040.root",mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(InPutFile.c_str());
  TGraphAsymmErrors *g_v2 = (TGraphAsymmErrors*)File_InPut->Get("g_v2");
  TF1 *f_v2 = new TF1("f_v2",v2_pT_FitFunc,momentumRange.first,momentumRange.second,5);
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

TF1* readspec(int energy)
{
  string InPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Phi_Spec.root",mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(InPutFile.c_str());
  TGraphAsymmErrors *g_spec = (TGraphAsymmErrors*)File_InPut->Get("g_spec");
  TF1 *f_Levy = new TF1("f_Levy",Levy,momentumRange.first,momentumRange.second,3);
  f_Levy->SetParameter(0,1);
  f_Levy->SetParameter(1,10);
  f_Levy->SetParameter(2,0.1);
  f_Levy->SetLineStyle(2);
  f_Levy->SetLineColor(4);
  f_Levy->SetLineWidth(2);
  g_spec->Fit(f_Levy,"N");

  TF1 *f_spec = new TF1("f_spec",pTLevy,momentumRange.first,momentumRange.second,3);
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
  double const pt = f_spec->GetRandom(momentumRange.first, momentumRange.second);
  double const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
  f_flow->SetParameter(0,f_v2->Eval(pt));
  double const phi = f_flow->GetRandom();

  // double const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
  // double const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
  // double const phi = TMath::TwoPi() * gRandom->Rndm();

  double const mT = sqrt(mass * mass + pt * pt);
  double const pz = mT * sinh(y);
  double const E = mT * cosh(y);

  lPhi.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
}

void setDecayChannels(int const mdme)
{
  for (int idc = decayChannels.first; idc < decayChannels.second + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0); // close all decay channel
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

  float Psi2 = EventPlaneSmearing(f_gaus);
  h_Psi2->Fill(Psi2);

  TVector3 nQSmear(sin(Psi2),-cos(Psi2),0.0); // direction of angular momentum, Psi2 = 0
  float CosThetaStarEP = vMcKpBoosted.Dot(nQSmear);
  float phiSmear = lPhi->Phi()-Psi2;
  if(phiSmear > TMath::Pi()) phiSmear -= TMath::TwoPi();
  if(phiSmear < -TMath::Pi()) phiSmear += TMath::TwoPi();
  h_phiEP->Fill(lPhi->Pt(),phiSmear);
  h_cosEP->Fill(lPhi->Pt(),CosThetaStarEP);
}

float EventPlaneSmearing(TF1 *f_gaus)
{
  float Psi2 = f_gaus->GetRandom(-TMath::PiOver2(),TMath::PiOver2());
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

bool Sampling(TF1 *f_rhoPhy,float CosThetaStar)
{
  float wMax;
  if(f_rhoPhy->GetParameter(0) <= 1.0/3.0) wMax = f_rhoPhy->Eval(0.0);
  else wMax = f_rhoPhy->Eval(1.0);
  return !(gRandom->Rndm() > f_rhoPhy->Eval(CosThetaStar)/wMax);
}

void write(int energy,int Nrho)
{
  string OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/McPhiV2_%d.root",mBeamEnergy[energy].c_str(),Nrho);
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_Tracks->Write();
  h_phiQA->Write();
  h_phiRP->Write();
  h_phiEP->Write();
  h_cosQA->Write();
  h_cosRP->Write();
  h_cosEP->Write();
  h_Psi2->Write();
  File_OutPut->Close();
}


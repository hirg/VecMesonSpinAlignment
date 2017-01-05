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
#include "TH1D.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TNtuple.h"
#include "../../Utility/StSpinAlignmentCons.h"

using namespace std;

typedef std::map<std::string,TH1D*> TH1DMap;

void readEfficiency(int energy, int year, int cut, int jobID);
void getKinematics(TLorentzVector& lPhi, double const mass);
void setDecayChannels(int const pid);
void decayAndFill(int const kf, TLorentzVector* lPhi, TClonesArray& daughters);
void fill(int const kf, TLorentzVector* lPhi, TLorentzVector const& lKplus, TLorentzVector const& lKminus);
bool tpcReconstructed(int iParticleIndex, int cent, TLorentzVector const& lKaon);
void findHist(TLorentzVector const& lKaon, int iParticleIndex, int& EtaBin, int& PhiBin); // iParticleIndex = 0 => K+
void write();

TPythia6Decayer* pydecay;
TNtuple* McPhiMeson;
TH1DMap h_EffKplus;
TH1DMap h_EffKminus;

TH1D *h_FrameEta[2];
TH1D *h_FramePhi[2];

TFile *File_OutPut;

void toyMcPhiDecay(const int energy = 4, const int pid = 0, const int year = 1, const int cut = 0, const int NMax = 100000, const int jobID = 3)
{
  TStopwatch* stopWatch = new TStopwatch();
  stopWatch->Start();
  gRandom->SetSeed();
  readEfficiency(energy,year,cut,jobID);

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
    // if (fabs(lPhi->Phi()) > TMath::Pi()) continue;
    decayAndFill(vmsa::decayMother[pid],lPhi,ptl);

    if (i_ran % 1000 == 1) McPhiMeson->AutoSave("SaveSelf");
  }
  cout << "=> processing data: 100%" << endl;
  cout << "work done!" << endl;

  write();

  stopWatch->Stop();   
  stopWatch->Print();
}

void getKinematics(TLorentzVector& lPhi, double const mass)
{
   double const pt = gRandom->Uniform(vmsa::ptMin, vmsa::ptEffMax);
   double const y = gRandom->Uniform(-vmsa::acceptanceRapidity, vmsa::acceptanceRapidity);
   double const phi = TMath::TwoPi() * gRandom->Rndm();

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

  fill(kf,lPhi,lKplus,lKminus);
}

void fill(int const kf, TLorentzVector* lPhi, TLorentzVector const& lKplus, TLorentzVector const& lKminus)
{
   int const centrality = floor(vmsa::NCentMax * gRandom->Rndm());
   TLorentzVector lRcPhi = lKplus + lKminus; // phi meson reconstruction
   // cout << "lPhi.pt = " << lPhi->Pt() << ", lPhi.eta = " << lPhi->Eta() << ", lPhi.phi = " << lPhi->Phi() << ", lPhi.m = " << lPhi->M() << endl;
   // cout << "lRcPhi.pt = " << lRcPhi.Pt() << ", lRcPhi.eta = " << lRcPhi.Eta() << ", lRcPhi.phi = " << lRcPhi.Phi() << ", lRcPhi.m = " << lRcPhi.M() << endl;
   // cout << endl;

   float arr[110];
   int iArr = 0;
   arr[iArr++] = centrality; // McPhi
   arr[iArr++] = lPhi->Pt();
   arr[iArr++] = lPhi->P();
   arr[iArr++] = lPhi->PseudoRapidity();
   arr[iArr++] = lPhi->Rapidity();
   arr[iArr++] = lPhi->Phi();
   arr[iArr++] = lPhi->M();
   arr[iArr++] = kf;

   arr[iArr++] = lKplus.Pt();
   arr[iArr++] = lKplus.PseudoRapidity();
   arr[iArr++] = lKplus.Rapidity();
   arr[iArr++] = lKplus.Phi();
   arr[iArr++] = lKplus.M();
   arr[iArr++] = 321;

   arr[iArr++] = lKplus.Pt();
   arr[iArr++] = lKplus.PseudoRapidity();
   arr[iArr++] = lKplus.Rapidity();
   arr[iArr++] = lKplus.Phi();
   arr[iArr++] = lKplus.M();
   arr[iArr++] = tpcReconstructed(0,centrality,lKplus);

   arr[iArr++] = lKminus.Pt();
   arr[iArr++] = lKminus.PseudoRapidity();
   arr[iArr++] = lKminus.Rapidity();
   arr[iArr++] = lKminus.Phi();
   arr[iArr++] = lKminus.M();
   arr[iArr++] = -321;

   arr[iArr++] = lKminus.Pt();
   arr[iArr++] = lKminus.PseudoRapidity();
   arr[iArr++] = lKminus.Rapidity();
   arr[iArr++] = lKminus.Phi();
   arr[iArr++] = lKminus.M();
   arr[iArr++] = tpcReconstructed(1,centrality,lKminus);

   arr[iArr++] = lRcPhi.Pt();
   arr[iArr++] = lRcPhi.P();
   arr[iArr++] = lRcPhi.PseudoRapidity();
   arr[iArr++] = lRcPhi.Rapidity();
   arr[iArr++] = lRcPhi.Phi();
   arr[iArr++] = lRcPhi.M();

   McPhiMeson->Fill(arr);
   // if(lRcPhi.Pt() < 10e-4) cout << "lPhi->Pt = " << lPhi->Pt() << ", lRcPhi.Pt = " << lRcPhi.Pt() << endl;
}

bool tpcReconstructed(int iParticleIndex, int cent, TLorentzVector const& lKaon)
{
   TH1D* h = NULL;
   int EtaBin = -1;
   int PhiBin = -1;

   if(fabs(lKaon.Eta()) > 1.0) return false;
   findHist(lKaon,iParticleIndex,EtaBin,PhiBin);

   // string KEY = Form("h_mRcEffPt_Cent_%d",cent);
   string KEY = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin,PhiBin);
   // cout << KEY.c_str() << endl;
   if (iParticleIndex == 0)
   {
     h = h_EffKplus[KEY];
   }
   else
   {
     h = h_EffKminus[KEY];
   }

   int const bin = h->FindBin(lKaon.Perp());

   return gRandom->Rndm() < h->GetBinContent(bin);
}

void findHist(TLorentzVector const& lKaon, int iParticleIndex, int& EtaBin, int& PhiBin)
{
  float eta = lKaon.Eta();
  EtaBin = h_FrameEta[iParticleIndex]->FindBin(eta)-1;
  // cout << "eta = " << eta << ", EtaBin = " << EtaBin << endl;
  float phi = lKaon.Phi();
  PhiBin = h_FramePhi[iParticleIndex]->FindBin(phi)-1;
}

void readEfficiency(int energy, int year, int cut, int jobID)
{
  string inputKplus = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[0].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str());
  TFile *File_Kplus = TFile::Open(inputKplus.c_str());
  cout << "OPEN Efficiency File for K+: " << inputKplus.c_str() << endl;

  string inputKminus = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[1].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str());
  TFile *File_Kminus = TFile::Open(inputKminus.c_str());
  cout << "OPEN Efficiency File for K-: " << inputKminus.c_str() << endl;

  h_FrameEta[0] = (TH1D*)File_Kplus->Get("h_FrameEta");
  h_FramePhi[0] = (TH1D*)File_Kplus->Get("h_FramePhi");
  h_FrameEta[1] = (TH1D*)File_Kminus->Get("h_FrameEta");
  h_FramePhi[1] = (TH1D*)File_Kminus->Get("h_FramePhi");

  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta)
    {
      for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
      {
	string KEY = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	// cout << "KEY = " << KEY.c_str() << endl;
	h_EffKplus[KEY] = (TH1D*)File_Kplus->Get(KEY.c_str());
	h_EffKminus[KEY] = (TH1D*)File_Kminus->Get(KEY.c_str());
      }	
    }
  }

  string outputfile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/Efficiency/Eff_%s_SingleKaon_%s_%s_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str(),jobID);
  cout << "OutPut File set to: " << outputfile.c_str() << endl;
  File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();

  int BufSize = (int)pow(2., 16.);
  // int Split = 1;

  const char* varlist = "Centrality:McPt:McP:McEta:McY:McPhi:McInvMass:McPid:" // MC phi 
                        "KpMcPt:KpMcEta:KpMcY:KpMcPhi:KpMcM:KpMcPid:" // MC K+ 
                        "KpRcPt:KpRcEta:KpRcY:KpRcPhi:KpRcM:KpRcTpc:" // RC K+
                        "KmMcPt:KmMcEta:KmMcY:KmMcPhi:KmMcM:KmMcPid:" // MC K-
                        "KmRcPt:KmRcEta:KmRcY:KmRcPhi:KmRcM:KmRcTpc:" // RC K-
                        "RcPt:RcP:RcEta:RcY:RcPhi:RcInvMass"; // reconstructed phi

  McPhiMeson = new TNtuple("McPhiMeson", "McPhiMeson", varlist, BufSize);
}

void write()
{
  File_OutPut->cd();
  McPhiMeson->Write();
  File_OutPut->Close();
}

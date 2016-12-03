#include "StEffHistManger.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"
#include <iostream>
#include "TLorentzVector.h"
// #include "TVector3.h"

ClassImp(StEffHistManger)
//
StEffHistManger::StEffHistManger()
{
  /* */
}

StEffHistManger::~StEffHistManger()
{
  /* */
}

void StEffHistManger::InitHist()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    std::string HistName = Form("h_mMcTracks_%d",i_cent);
    h_mMcTracks[i_cent] = new TH3D(HistName.c_str(),HistName.c_str(),Efficiency::BinPt,0.0,Efficiency::ptMax,Efficiency::BinEta,-1.0,1.0,Efficiency::BinPhi,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mRcTracks_%d",i_cent);
    h_mRcTracks[i_cent] = new TH3D(HistName.c_str(),HistName.c_str(),Efficiency::BinPt,0.0,Efficiency::ptMax,Efficiency::BinEta,-1.0,1.0,Efficiency::BinPhi,-1.0*TMath::Pi(),TMath::Pi());
    for(int i_pt = 0; i_pt < Efficiency::BinPt; ++i_pt)
    {
      HistName = Form("h_mMcEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mMcEffCos[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),Efficiency::BinCos,-1.0,1.0);
      HistName = Form("h_mRcEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mRcEffCos[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),Efficiency::BinCos,-1.0,1.0);
    }
  }
  h_FrameEta = new TH1D("h_FrameEta","h_FrameEta",Efficiency::BinEta,-1.0,1.0);
  h_FramePhi = new TH1D("h_FramePhi","h_FramePhi",Efficiency::BinPhi,-TMath::Pi(),TMath::Pi());

  flag_eff = 0;
  flag_eff_PtEtaPhi = 0;
  flag_eff_Cos = 0;
  flag_QA_Cos = 0;
  // InitQACos();
  InitQAEP();
}

void StEffHistManger::InitQACos()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_pt = 0; i_pt < Efficiency::BinPt; ++i_pt)
    {
      std::string HistName;

      HistName = Form("h_mMcKpCosStar_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mMcKpCosStar[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),10.0*Efficiency::BinCos,-1.0,1.0);
      HistName = Form("h_mMcKpCos_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mMcKpCos[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),10.0*Efficiency::BinCos,-1.0,1.0);
      HistName = Form("h_mMcKpCosTheta_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mMcKpCosTheta[i_cent][i_pt] = new TH2D(HistName.c_str(),HistName.c_str(),10.0*Efficiency::BinCos,-1.0,1.0,10.0*Efficiency::BinCos,-1.0,1.0);

      HistName = Form("h_mMcKmCosStar_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mMcKmCosStar[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),10.0*Efficiency::BinCos,-1.0,1.0);
      HistName = Form("h_mMcKmCos_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mMcKmCos[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),10.0*Efficiency::BinCos,-1.0,1.0);
      HistName = Form("h_mMcKmCosTheta_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mMcKmCosTheta[i_cent][i_pt] = new TH2D(HistName.c_str(),HistName.c_str(),10.0*Efficiency::BinCos,-1.0,1.0,10.0*Efficiency::BinCos,-1.0,1.0);

      HistName = Form("h_mMcKpKmCos_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mMcKpKmCos[i_cent][i_pt] = new TH2D(HistName.c_str(),HistName.c_str(),10.0*Efficiency::BinCos,-1.0,1.0,10.0*Efficiency::BinCos,-1.0,1.0);

      HistName = Form("h_mRcKpCosStar_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mRcKpCosStar[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),10.0*Efficiency::BinCos,-1.0,1.0);
      HistName = Form("h_mRcKpCos_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mRcKpCos[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),10.0*Efficiency::BinCos,-1.0,1.0);
      HistName = Form("h_mRcKpCosTheta_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mRcKpCosTheta[i_cent][i_pt] = new TH2D(HistName.c_str(),HistName.c_str(),10.0*Efficiency::BinCos,-1.0,1.0,10.0*Efficiency::BinCos,-1.0,1.0);

      HistName = Form("h_mRcKmCosStar_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mRcKmCosStar[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),10.0*Efficiency::BinCos,-1.0,1.0);
      HistName = Form("h_mRcKmCos_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mRcKmCos[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),10.0*Efficiency::BinCos,-1.0,1.0);
      HistName = Form("h_mRcKmCosTheta_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mRcKmCosTheta[i_cent][i_pt] = new TH2D(HistName.c_str(),HistName.c_str(),10.0*Efficiency::BinCos,-1.0,1.0,10.0*Efficiency::BinCos,-1.0,1.0);

      HistName = Form("h_mRcKpKmCos_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mRcKpKmCos[i_cent][i_pt] = new TH2D(HistName.c_str(),HistName.c_str(),10.0*Efficiency::BinCos,-1.0,1.0,10.0*Efficiency::BinCos,-1.0,1.0);
    }
  }
  flag_QA_Cos = 1;
}

void StEffHistManger::InitQAEP()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    std::string HistName = Form("h_mMcEP_%d",i_cent);
    h_mMcEP[i_cent] = new TH2D(HistName.c_str(),HistName.c_str(),Efficiency::BinPt,0.0,Efficiency::ptMax,7,0.0,0.5*TMath::Pi());
    HistName = Form("h_mRcEP_%d",i_cent);
    h_mRcEP[i_cent] = new TH2D(HistName.c_str(),HistName.c_str(),Efficiency::BinPt,0.0,Efficiency::ptMax,7,0.0,0.5*TMath::Pi());
  }
  flag_QA_EP = 0;
}

void StEffHistManger::FillHistMc(int cent, float pt, float eta, float phi, float cos)
{
  h_mMcTracks[cent]->Fill(pt,eta,phi);
  h_mMcTracks[9]->Fill(pt,eta,phi);
  float delta_pt = Efficiency::ptMax/((float)Efficiency::BinPt);
  for(int i_pt = 0; i_pt < Efficiency::BinPt; ++i_pt)
  {
    if(!(pt > i_pt*delta_pt && pt < (i_pt+1)*delta_pt)) continue;
    h_mMcEffCos[cent][i_pt]->Fill(cos);
    h_mMcEffCos[9][i_pt]->Fill(cos);
  }
}

void StEffHistManger::FillHistRc(int cent, float pt, float eta, float phi, float cos)
{
  h_mRcTracks[cent]->Fill(pt,eta,phi);
  h_mRcTracks[9]->Fill(pt,eta,phi);
  float delta_pt = Efficiency::ptMax/((float)Efficiency::BinPt);
  for(int i_pt = 0; i_pt < Efficiency::BinPt; ++i_pt)
  {
    if(!(pt > i_pt*delta_pt && pt < (i_pt+1)*delta_pt)) continue;
    h_mRcEffCos[cent][i_pt]->Fill(cos);
    h_mRcEffCos[9][i_pt]->Fill(cos);
  }
}

void StEffHistManger::FillQACosMc(int cent, McVecMeson McPhi, McDecayDau McKP, McDecayDau McKM)
{
  TLorentzVector lMcPhi;
  lMcPhi.SetPtEtaPhiM(McPhi.pt,McPhi.eta,McPhi.phi,Efficiency::mMassPhi);
  TVector3 vMcPhi = lMcPhi.Vect().Unit(); // momentum direction of phi

  TLorentzVector lMcKP;
  lMcKP.SetPtEtaPhiM(McKP.pt,McKP.eta,McKP.phi,Efficiency::mMassKaon);
  TVector3 vMcKP = lMcKP.Vect().Unit(); // momentum direction of Kplus in lab frame
  float CosThetaKp = vMcKP.Dot(vMcPhi);
  TVector3 vMcKPStar = CalBoostedVector(McKP,McPhi); // momentum direction of Kplus in phi-meson rest frame
  float CosThetaKpStar = vMcKPStar.Dot(vMcPhi);

  TLorentzVector lMcKM;
  lMcKM.SetPtEtaPhiM(McKM.pt,McKM.eta,McKM.phi,Efficiency::mMassKaon);
  TVector3 vMcKM = lMcKM.Vect().Unit();
  float CosThetaKm = vMcKM.Dot(vMcPhi);
  TVector3 vMcKMStar = CalBoostedVector(McKM,McPhi);
  float CosThetaKmStar = vMcKMStar.Dot(vMcPhi);

  float CosThetaKpKm = vMcKP.Dot(vMcKM);

  float delta_pt = Efficiency::ptMax/((float)Efficiency::BinPt);
  for(int i_pt = 0; i_pt < Efficiency::BinPt; ++i_pt)
  {
    if(!(McPhi.pt > i_pt*delta_pt && McPhi.pt < (i_pt+1)*delta_pt)) continue;

    h_mMcKpCosStar[cent][i_pt]->Fill(CosThetaKpStar);
    h_mMcKpCos[cent][i_pt]->Fill(CosThetaKp);
    h_mMcKpCosTheta[cent][i_pt]->Fill(CosThetaKpStar,CosThetaKp);

    h_mMcKmCosStar[cent][i_pt]->Fill(CosThetaKmStar);
    h_mMcKmCos[cent][i_pt]->Fill(CosThetaKm);
    h_mMcKmCosTheta[cent][i_pt]->Fill(CosThetaKmStar,CosThetaKm);

    h_mMcKpKmCos[cent][i_pt]->Fill(CosThetaKp,CosThetaKm);

    h_mMcKpCosStar[9][i_pt]->Fill(CosThetaKpStar);
    h_mMcKpCos[9][i_pt]->Fill(CosThetaKp);
    h_mMcKpCosTheta[9][i_pt]->Fill(CosThetaKpStar,CosThetaKp);

    h_mMcKmCosStar[9][i_pt]->Fill(CosThetaKmStar);
    h_mMcKmCos[9][i_pt]->Fill(CosThetaKm);
    h_mMcKmCosTheta[9][i_pt]->Fill(CosThetaKmStar,CosThetaKm);

    h_mMcKpKmCos[9][i_pt]->Fill(CosThetaKp,CosThetaKm);
  }
}

void StEffHistManger::FillQACosRc(int cent, McVecMeson McPhi, McDecayDau McKP, McDecayDau McKM)
{
  TLorentzVector lMcPhi;
  lMcPhi.SetPtEtaPhiM(McPhi.pt,McPhi.eta,McPhi.phi,Efficiency::mMassPhi);
  TVector3 vMcPhi = lMcPhi.Vect().Unit(); // momentum direction of phi

  TLorentzVector lMcKP;
  lMcKP.SetPtEtaPhiM(McKP.pt,McKP.eta,McKP.phi,Efficiency::mMassKaon);
  TVector3 vMcKP = lMcKP.Vect().Unit(); // momentum direction of Kplus in lab frame
  float CosThetaKp = vMcKP.Dot(vMcPhi);
  TVector3 vMcKPStar = CalBoostedVector(McKP,McPhi); // momentum direction of Kplus in phi-meson rest frame
  float CosThetaKpStar = vMcKPStar.Dot(vMcPhi);

  TLorentzVector lMcKM;
  lMcKM.SetPtEtaPhiM(McKM.pt,McKM.eta,McKM.phi,Efficiency::mMassKaon);
  TVector3 vMcKM = lMcKM.Vect().Unit();
  float CosThetaKm = vMcKM.Dot(vMcPhi);
  TVector3 vMcKMStar = CalBoostedVector(McKM,McPhi);
  float CosThetaKmStar = vMcKMStar.Dot(vMcPhi);

  float CosThetaKpKm = vMcKP.Dot(vMcKM);

  float delta_pt = Efficiency::ptMax/((float)Efficiency::BinPt);
  for(int i_pt = 0; i_pt < Efficiency::BinPt; ++i_pt)
  {
    if(!(McPhi.pt > i_pt*delta_pt && McPhi.pt < (i_pt+1)*delta_pt)) continue;

    h_mRcKpCosStar[cent][i_pt]->Fill(CosThetaKpStar);
    h_mRcKpCos[cent][i_pt]->Fill(CosThetaKp);
    h_mRcKpCosTheta[cent][i_pt]->Fill(CosThetaKpStar,CosThetaKp);

    h_mRcKmCosStar[cent][i_pt]->Fill(CosThetaKmStar);
    h_mRcKmCos[cent][i_pt]->Fill(CosThetaKm);
    h_mRcKmCosTheta[cent][i_pt]->Fill(CosThetaKmStar,CosThetaKm);

    h_mRcKpKmCos[cent][i_pt]->Fill(CosThetaKp,CosThetaKm);

    h_mRcKpCosStar[9][i_pt]->Fill(CosThetaKpStar);
    h_mRcKpCos[9][i_pt]->Fill(CosThetaKp);
    h_mRcKpCosTheta[9][i_pt]->Fill(CosThetaKpStar,CosThetaKp);

    h_mRcKmCosStar[9][i_pt]->Fill(CosThetaKmStar);
    h_mRcKmCos[9][i_pt]->Fill(CosThetaKm);
    h_mRcKmCosTheta[9][i_pt]->Fill(CosThetaKmStar,CosThetaKm);

    h_mRcKpKmCos[9][i_pt]->Fill(CosThetaKp,CosThetaKm);
  }
}

void StEffHistManger::FillQAEPMc(int cent, float pt, float phi)
{
  float phi_final = -999.0;
  for(int psi_bin = 0; psi_bin < 3; ++psi_bin)
  {
    if(phi >= Efficiency::Psi2_low[psi_bin] && phi < Efficiency::Psi2_up[psi_bin])
    {
      phi_final = phi - (psi_bin-1)*2.0*TMath::Pi()/2.0;
    }
  }
  h_mMcEP[cent]->Fill(pt,TMath::Abs(phi_final));
  if(cent >= Efficiency::cent_low[0] && cent <= Efficiency::cent_high[0]) h_mMcEP[9]->Fill(pt,TMath::Abs(phi_final));
}

void StEffHistManger::FillQAEPRc(int cent, float pt, float phi)
{
  float phi_final = -999.0;
  for(int psi_bin = 0; psi_bin < 3; ++psi_bin)
  {
    if(phi >= Efficiency::Psi2_low[psi_bin] && phi < Efficiency::Psi2_up[psi_bin])
    {
      phi_final = phi - (psi_bin-1)*2.0*TMath::Pi()/2.0;
    }
  }
  h_mRcEP[cent]->Fill(pt,TMath::Abs(phi_final));
  if(cent >= Efficiency::cent_low[0] && cent <= Efficiency::cent_high[0]) h_mRcEP[9]->Fill(pt,TMath::Abs(phi_final));
}

TH1D* StEffHistManger::CalEffError(TH1D *h_Mc, TH1D *h_Rc, std::string HistName)
{
  TH1D* h_ratio = (TH1D*)h_Rc->Clone();
  h_ratio->Divide(h_Mc);
  for(int i_bin = 1; i_bin < h_ratio->GetNbinsX()+1; ++i_bin)
  {
    double n = h_Mc->GetBinContent(i_bin);
    double k = h_Rc->GetBinContent(i_bin);
    double variance = (k+1.0)*(k+2.0)/((n+2.0)*(n+3.0))-(k+1.0)*(k+1.0)/((n+2.0)*(n+2.0));
    double sigma = TMath::Sqrt(variance);
    if(n > 0.0 && k > 0.0) h_ratio->SetBinError(i_bin,sigma);
  }
  h_ratio->SetName(HistName.c_str());
  
  return h_ratio;
}

void StEffHistManger::CalEfficiency()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    std::string HistName;

    HistName = Form("h_mMcEffPt_Cent_%d",i_cent);
    h_mMcEffPt[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("x")->Clone(HistName.c_str());
    HistName = Form("h_mRcEffPt_Cent_%d",i_cent);
    h_mRcEffPt[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("x")->Clone(HistName.c_str());
    HistName = Form("h_mEffPt_Cent_%d",i_cent);
    h_mEffPt[i_cent] = CalEffError(h_mMcEffPt[i_cent],h_mRcEffPt[i_cent],HistName.c_str());

    HistName = Form("h_mMcEffEta_Cent_%d",i_cent);
    h_mMcEffEta[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("y")->Clone(HistName.c_str());
    HistName = Form("h_mRcEffEta_Cent_%d",i_cent);
    h_mRcEffEta[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("y")->Clone(HistName.c_str());
    HistName = Form("h_mEffEta_Cent_%d",i_cent);
    h_mEffEta[i_cent] = CalEffError(h_mMcEffEta[i_cent],h_mRcEffEta[i_cent],HistName.c_str());

    HistName = Form("h_mMcEffPhi_Cent_%d",i_cent);
    h_mMcEffPhi[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("z")->Clone(HistName.c_str());
    HistName = Form("h_mRcEffPhi_Cent_%d",i_cent);
    h_mRcEffPhi[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("z")->Clone(HistName.c_str());
    HistName = Form("h_mEffPhi_Cent_%d",i_cent);
    h_mEffPhi[i_cent] = CalEffError(h_mMcEffPhi[i_cent],h_mRcEffPhi[i_cent],HistName.c_str());
  }
  flag_eff = 1;
}

void StEffHistManger::CalEffPtEtaPhi()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_eta = 0; i_eta < Efficiency::BinEta; ++i_eta)
    {
      for(int i_phi = 0; i_phi < Efficiency::BinPhi; ++i_phi)
      {
	std::string HistNameMc = Form("h_mMcEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_mMcEffPEP[HistNameMc] = (TH1D*)h_mMcTracks[i_cent]->ProjectionX(HistNameMc.c_str(),i_eta+1,i_eta+1,i_phi+1,i_phi+1);

	std::string HistNameRc = Form("h_mRcEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_mRcEffPEP[HistNameRc] = (TH1D*)h_mRcTracks[i_cent]->ProjectionX(HistNameRc.c_str(),i_eta+1,i_eta+1,i_phi+1,i_phi+1);

	std::string HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_mEffPEP[HistNameEff] = CalEffError(h_mMcEffPEP[HistNameMc],h_mRcEffPEP[HistNameRc],HistNameEff.c_str());
      }
    }
  }
  flag_eff_PtEtaPhi = 1;
}

void StEffHistManger::CalEffCosThetaStar()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_pt = 0; i_pt < Efficiency::BinPt; ++i_pt)
    {
      std::string HistName = Form("h_mEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mEffCos[i_cent][i_pt] = CalEffError(h_mMcEffCos[i_cent][i_pt],h_mRcEffCos[i_cent][i_pt],HistName.c_str());
    }
  }
  flag_eff_Cos = 1;
}

void StEffHistManger::CalEffEP()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_pt = 0; i_pt < Efficiency::BinPt; ++i_pt)
    {
      std::string HistNameMc = Form("h_mMcEP_Cent_%d_Pt_%d",i_cent,i_pt);
      TH1D *h_mMc = (TH1D*)h_mMcEP[i_cent]->ProjectionY(HistNameMc.c_str(),i_pt+1,i_pt+1);
      std::string HistNameRc = Form("h_mRcEP_Cent_%d_Pt_%d",i_cent,i_pt);
      TH1D *h_mRc = (TH1D*)h_mRcEP[i_cent]->ProjectionY(HistNameRc.c_str(),i_pt+1,i_pt+1);

      std::string HistName = Form("h_mEffEP_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mEffEP[i_cent][i_pt] = CalEffError(h_mMc,h_mRc,HistName.c_str());
    }
  }
  flag_QA_EP = 1;
}

TVector3 StEffHistManger::CalBoostedVector(McDecayDau dau, McVecMeson vec)
{
  TLorentzVector lMcVec;
  lMcVec.SetPtEtaPhiM(vec.pt,vec.eta,vec.phi,Efficiency::mMassPhi);
  TVector3 vMcBeta = -1.0*lMcVec.BoostVector(); // boost vector

  TLorentzVector lMcDau;
  lMcDau.SetPtEtaPhiM(dau.pt,dau.eta,dau.phi,Efficiency::mMassKaon);
  lMcDau.Boost(vMcBeta); // boost Kplus back to phi-meson rest frame
  TVector3 vMcDauStar = lMcDau.Vect().Unit(); // momentum direction of Kplus in phi-meson rest frame

  return vMcDauStar;
}

void StEffHistManger::WriteHist()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    h_mMcTracks[i_cent]->Write();
    h_mRcTracks[i_cent]->Write();
    if(flag_eff > 0)
    {
      h_mEffPt[i_cent]->Write();
      h_mEffEta[i_cent]->Write();
      h_mEffPhi[i_cent]->Write();
    }

    if(flag_eff_PtEtaPhi > 0)
    {
      for(int i_eta = 0; i_eta < Efficiency::BinEta; ++i_eta)
      {
	for(int i_phi = 0; i_phi < Efficiency::BinPhi; ++i_phi)
	{
	  std::string HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	  h_mEffPEP[HistNameEff]->Write();
	}
      }
    }

    for(int i_pt = 0; i_pt < Efficiency::BinPt; ++i_pt)
    {
      if(flag_eff_Cos > 0)
      {
	// h_mMcEffCos[i_cent][i_pt]->Write();
	// h_mRcEffCos[i_cent][i_pt]->Write();
	h_mEffCos[i_cent][i_pt]->Write();
      }

      if(flag_QA_Cos > 0)
      {
	h_mMcKpCosStar[i_cent][i_pt]->Write();
	h_mMcKpCos[i_cent][i_pt]->Write();
	h_mMcKpCosTheta[i_cent][i_pt]->Write();

	h_mMcKmCosStar[i_cent][i_pt]->Write();
	h_mMcKmCos[i_cent][i_pt]->Write();
	h_mMcKmCosTheta[i_cent][i_pt]->Write();

	h_mMcKpKmCos[i_cent][i_pt]->Write();

	h_mRcKpCosStar[i_cent][i_pt]->Write();
	h_mRcKpCos[i_cent][i_pt]->Write();
	h_mRcKpCosTheta[i_cent][i_pt]->Write();

	h_mRcKmCosStar[i_cent][i_pt]->Write();
	h_mRcKmCos[i_cent][i_pt]->Write();
	h_mRcKmCosTheta[i_cent][i_pt]->Write();

	h_mRcKpKmCos[i_cent][i_pt]->Write();
      }

      if(flag_QA_EP > 0)
      {
	h_mEffEP[i_cent][i_pt]->Write();
      }
    }
    if(flag_QA_EP > 0)
    {
      h_mMcEP[i_cent]->Write();
      h_mRcEP[i_cent]->Write();
    }
  }
  h_FrameEta->Write();
  h_FramePhi->Write();
}

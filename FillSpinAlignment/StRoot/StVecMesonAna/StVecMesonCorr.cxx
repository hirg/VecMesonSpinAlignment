#include "StRoot/StVecMesonAna/StVecMesonCorr.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include "StMessMgr.h"

Double_t Resolution_Full(Double_t *x_val, Double_t *par)
{
  Double_t y;
  Double_t chi = x_val[0];
  Double_t arg = chi*chi/4.0;
  Double_t norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return y;
}

ClassImp(StVecMesonCorr)

TString StVecMesonCorr::mVStr[2] = {"pos","neg"};
TString StVecMesonCorr::mOrder = "2nd";
//--------------------------------------------------------
StVecMesonCorr::StVecMesonCorr(Int_t energy)
{
  mEnergy = energy;
}

StVecMesonCorr::~StVecMesonCorr()
{
}
//--------------------------------------------------------
// ReCenter Correction
void StVecMesonCorr::InitReCenterCorrection()
{
  TString InPutFile_ReCenter = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ReCenterParameter/file_%s_ReCenterPar.root",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mBeamEnergy[mEnergy].c_str());

  mInPutFile_ReCenter = TFile::Open(InPutFile_ReCenter.Data());
}

// get ReCenter Parameter
TVector2 StVecMesonCorr::getReCenterPar_East(Int_t Cent9, Int_t RunIndex, Int_t vz_sign)
{
  Float_t mean_qx, mean_qy;
  TVector2 qVector(0.0,0.0);

  TString ProName_x = Form("qx_%s_Vertex_%s_East",mOrder.Data(),mVStr[vz_sign].Data());
  TString ProName_y = Form("qy_%s_Vertex_%s_East",mOrder.Data(),mVStr[vz_sign].Data());

  TProfile2D *p_x = (TProfile2D*)mInPutFile_ReCenter->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile_ReCenter->Get(ProName_y.Data());

  mean_qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex,(Double_t)Cent9));
  mean_qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex,(Double_t)Cent9));

  qVector.Set(mean_qx,mean_qy);

  return qVector;
}

TVector2 StVecMesonCorr::getReCenterPar_West(Int_t Cent9, Int_t RunIndex, Int_t vz_sign)
{
  Float_t mean_qx, mean_qy;
  TVector2 qVector(0.0,0.0);

  TString ProName_x = Form("qx_%s_Vertex_%s_West",mOrder.Data(),mVStr[vz_sign].Data());
  TString ProName_y = Form("qy_%s_Vertex_%s_West",mOrder.Data(),mVStr[vz_sign].Data());

  TProfile2D *p_x = (TProfile2D*)mInPutFile_ReCenter->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile_ReCenter->Get(ProName_y.Data());

  mean_qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex,(Double_t)Cent9));
  mean_qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex,(Double_t)Cent9));

  qVector.Set(mean_qx,mean_qy);

  return qVector;
}

TVector2 StVecMesonCorr::getReCenterPar_Full(Int_t Cent9, Int_t RunIndex, Int_t vz_sign)
{
  Float_t mean_qx, mean_qy;
  TVector2 qVector(0.0,0.0);

  TString ProName_x = Form("qx_%s_Vertex_%s_Full",mOrder.Data(),mVStr[vz_sign].Data());
  TString ProName_y = Form("qy_%s_Vertex_%s_Full",mOrder.Data(),mVStr[vz_sign].Data());

  TProfile2D *p_x = (TProfile2D*)mInPutFile_ReCenter->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile_ReCenter->Get(ProName_y.Data());

  mean_qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex,(Double_t)Cent9));
  mean_qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex,(Double_t)Cent9));

  qVector.Set(mean_qx,mean_qy);

  return qVector;
}

// calculate qVector of Track
TVector2 StVecMesonCorr::calq2Vector(TLorentzVector lTrack)
{
  const Float_t phi = lTrack.Phi();
  TVector2 q2Vector(0.0,0.0);

  const Float_t q2x = TMath::Cos(2.0*phi);
  const Float_t q2y = TMath::Sin(2.0*phi);
  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

Float_t StVecMesonCorr::getWeight(TLorentzVector lTrack)
{
  Float_t pt = lTrack.Perp();
  Float_t w;
  if(pt <= vmsa::mPrimPtWeight)
  {
    w = pt;
  }
  if(pt > vmsa::mPrimPtWeight)
  {
    w = vmsa::mPrimPtWeight;
  }

  return w;
}
//--------------------------------------------------------
// Shift Correction
void StVecMesonCorr::InitShiftCorrection()
{
  TString InPutFile_Shift = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ShiftParameter/file_%s_ShiftPar.root",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mBeamEnergy[mEnergy].c_str());
  mInPutFile_Shift = TFile::Open(InPutFile_Shift.Data());
}

bool StVecMesonCorr::passTrackEtaNumCut(Int_t NumTrackEast, Int_t NumTrackWest)
{
  if(!(NumTrackEast > vmsa::mTrackMin && NumTrackWest > vmsa::mTrackMin))
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StVecMesonCorr::passTrackFullNumCut(Int_t NumTrackFull, Int_t NumTrackFull_East, Int_t NumTrackFull_West)
{
  if(!(NumTrackFull > vmsa::mTrackMin_Full && NumTrackFull_East > 0 && NumTrackFull_West > 0))
  {
    return kFALSE;
  }

  return kTRUE;
}

Float_t StVecMesonCorr::AngleShift(Float_t Psi_raw)
{
  Float_t Psi_Corr = Psi_raw;
  if(Psi_raw > 0.5*TMath::Pi())
  {
    Psi_Corr = Psi_raw - TMath::Pi();
  }
  if(Psi_raw < -0.5*TMath::Pi())
  {
    Psi_Corr = Psi_raw + TMath::Pi();
  }

  return Psi_Corr;
}

// calculate EP angle after Shift Correction
// 2nd
Float_t StVecMesonCorr::calShiftAngle2East_EP(TVector2 Q2Vector_East, Int_t runIndex, Int_t Cent9, Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(Q2Vector_East.Y(),Q2Vector_East.X())/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_East",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_East",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(vmsa::mShiftOrder[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(vmsa::mShiftOrder[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  return Psi_Shift;
}

Float_t StVecMesonCorr::calShiftAngle2West_EP(TVector2 Q2Vector_West, Int_t runIndex, Int_t Cent9, Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(Q2Vector_West.Y(),Q2Vector_West.X())/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_West",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_West",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(vmsa::mShiftOrder[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(vmsa::mShiftOrder[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  return Psi_Shift;
}

Float_t StVecMesonCorr::calShiftAngle2Full_EP(TVector2 Q2Vector_Full, Int_t runIndex, Int_t Cent9, Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(Q2Vector_Full.Y(),Q2Vector_Full.X())/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_Full",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_Full",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(vmsa::mShiftOrder[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(vmsa::mShiftOrder[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  return Psi_Shift;
}

Float_t StVecMesonCorr::calShiftAngle2Full_EP(TVector2 Q2Vector_Full, Int_t runIndex, Int_t Cent9, Int_t vz_sign, TLorentzVector lTrack, bool UsedEP)
{
  TVector2 QVector_sub = Q2Vector_Full;
  if(UsedEP)
  {
    Float_t w = getWeight(lTrack);
    QVector_sub -= w*(calq2Vector(lTrack) - getReCenterPar_Full(Cent9,runIndex,vz_sign));
  }

  Float_t Psi_ReCenter = TMath::ATan2(QVector_sub.Y(),QVector_sub.X())/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(vmsa::mShiftOrder[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(vmsa::mShiftOrder[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  return Psi_Shift;
}

//--------------------------------------------------------
// Resolution Correction
void StVecMesonCorr::InitResolutionCorr()
{
  TString InPutFile_Res = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Resolution/file_%s_Resolution.root",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mBeamEnergy[mEnergy].c_str());
  mInPutFile_Res = TFile::Open(InPutFile_Res.Data());
}

Float_t StVecMesonCorr::getResolution2_EP(Int_t Cent9)
{
  TProfile *p_res2 = (TProfile*)mInPutFile_Res->Get("p_mRes2_Sub");
  Float_t Res_raw = p_res2->GetBinContent(p_res2->FindBin(Cent9));
  if(Res_raw <= 0)
  {
    return -999.9;
  }
  else
  {
    Float_t Res = TMath::Sqrt(Res_raw);
    return Res;
  }
}
Float_t StVecMesonCorr::getResolution2_Full_EP(Int_t Cent9)
{
  TProfile *p_res2 = (TProfile*)mInPutFile_Res->Get("p_mRes2_Ran");
  Float_t Res_raw = p_res2->GetBinContent(p_res2->FindBin(Cent9));
  if(Res_raw <= 0)
  {
    return -999.9;
  }
  else
  {
    Float_t Res_sub = TMath::Sqrt(Res_raw);
    TF1 *f_res = new TF1("f_res",Resolution_Full,0,10,0);
    Float_t chi_sub = f_res->GetX(Res_sub);
    Float_t chi_full = chi_sub*TMath::Sqrt(2.0);
    Float_t Res_full = f_res->Eval(chi_full);
    return Res_full;
  }
}


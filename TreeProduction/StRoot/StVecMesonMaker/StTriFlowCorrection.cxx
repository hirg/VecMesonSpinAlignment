#include "StTriFlowCorrection.h"
#include "StTriFlowConstants.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StMessMgr.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"

Double_t Resolution_Full(Double_t *x_val, Double_t *par)
{
  Double_t y;
  Double_t chi = x_val[0];
  Double_t arg = chi*chi/4.0;
  Double_t norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return y;
}

ClassImp(StTriFlowCorrection)

TString StTriFlowCorrection::mVStr[2] = {"pos","neg"};
TString StTriFlowCorrection::mOrder[2] = {"2nd","3rd"};
TString StTriFlowCorrection::mMethod[2] = {"EP","SP"};
//---------------------------------------------------------------------------------

StTriFlowCorrection::StTriFlowCorrection(Int_t energy)
{
  mEnergy = energy;
}

//---------------------------------------------------------------------------------

StTriFlowCorrection::~StTriFlowCorrection()
{
  /* */
}

//---------------------------------------------------------------------------------

void StTriFlowCorrection::InitReCenterCorrection(Int_t mEnergy)
{
  TString InPutFile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/RecenterParameter/file_%s_ReCenterPar.root",TriFlow::Energy[mEnergy].Data(),TriFlow::Energy[mEnergy].Data());

  mInPutFile = TFile::Open(InPutFile.Data());

  for(Int_t i = 0; i < 4; i++)
  {
    mQ2Vector_East_EP[i].Set(0.0,0.0);
    mQ3Vector_East_EP[i].Set(0.0,0.0);
    mQ2Vector_East_SP[i].Set(0.0,0.0);
    mQ3Vector_East_SP[i].Set(0.0,0.0);
    mQCounter_East[i] = 0;

    mQ2Vector_West_EP[i].Set(0.0,0.0);
    mQ3Vector_West_EP[i].Set(0.0,0.0);
    mQ2Vector_West_SP[i].Set(0.0,0.0);
    mQ3Vector_West_SP[i].Set(0.0,0.0);
    mQCounter_West[i] = 0;
  }

  mQ2Vector_Full_EP.Set(0.0,0.0);
  mQ3Vector_Full_EP.Set(0.0,0.0);
  mQ2Vector_Full_SP.Set(0.0,0.0);
  mQ3Vector_Full_SP.Set(0.0,0.0);
  mQCounter_Full = 0;
  mQCounter_Full_East = 0;
  mQCounter_Full_West = 0;

  mQ2Vector_A_EP.Set(0.0,0.0);
  mQ3Vector_A_EP.Set(0.0,0.0);
  mQ2Vector_A_SP.Set(0.0,0.0);
  mQ3Vector_A_SP.Set(0.0,0.0);
  mQCounter_A = 0;

  mQ2Vector_B_EP.Set(0.0,0.0);
  mQ3Vector_B_EP.Set(0.0,0.0);
  mQ2Vector_B_SP.Set(0.0,0.0);
  mQ3Vector_B_SP.Set(0.0,0.0);
  mQCounter_B = 0;
}

//---------------------------------------------------------------------------------

void StTriFlowCorrection::InitShiftCorrection(Int_t mEnergy)
{
  TString InPutFile_Shift = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Correction/Shift/file_%s_Corr_Shift.root",TriFlow::Energy[mEnergy].Data(),TriFlow::Energy[mEnergy].Data());
  mInPutFile_Shift = TFile::Open(InPutFile_Shift.Data());

  TString InPutFile_Res = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Resolution/file_%s_Resolution.root",TriFlow::Energy[mEnergy].Data(),TriFlow::Energy[mEnergy].Data());
  mInPutFile_Res = TFile::Open(InPutFile_Res.Data());
}

//---------------------------------------------------------------------------------

bool StTriFlowCorrection::passTrackEtaEast(StPicoTrack *track, Int_t i, Int_t Mode) // neg || i = different eta_gap || Mode = 0 Event Plane Mode, Mode = 1 Flow Mode
{
  Float_t eta = track->pMom().pseudoRapidity();
  
  if(Mode == 0) // Event Plane Mode
  {
    // eta cut
    // eta_gap between two sub event plane is 2*mEta_Gap[i]
    if(!(eta > -1.0*TriFlow::mEtaMax && eta < -1.0*TriFlow::mEta_Gap[i]))
    {
      return kFALSE;
    }

    return kTRUE;
  }
  if(Mode == 1) // Flow Mode
  {
    // eta cut
    // eta_gap between two sub event plane is mEta_Gap[i]
    if(!(eta > -1.0*TriFlow::mEtaMax && eta < 0.0))
    {
      return kFALSE;
    }

    return kTRUE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------

bool StTriFlowCorrection::passTrackEtaWest(StPicoTrack *track, Int_t i, Int_t Mode) // pos || i = different eta_gap || Mode = 0 Event Plane Mode, Mode = 1 Flow Mode
{
  Float_t eta = track->pMom().pseudoRapidity();

  if(Mode == 0) // Event Plane Mode
  {
    // eta cut
    // eta_gap between two sub event plane is 2*mEta_Gap[i]
    if(!(eta > TriFlow::mEta_Gap[i] && eta < TriFlow::mEtaMax))
    {
      return kFALSE;
    }

    return kTRUE;
  }
  if(Mode == 1) // Flow Mode
  {
    // eta cut
    // eta_gap between two sub event plane is mEta_Gap[i]
    if(!(eta > 0.0 && eta < TriFlow::mEtaMax))
    {
      return kFALSE;
    }

    return kTRUE;
  }

  return kTRUE;
}
//---------------------------------------------------------------------------------

bool StTriFlowCorrection::passTrackFull(StPicoTrack *track) // Full Event Plane pt Cut
{
  // pt cut 0.2 - 2.0 GeV/c
  Float_t pt = track->pMom().perp();
  if(!(pt > TriFlow::mPrimPtMin[mEnergy] && pt < TriFlow::mPrimPtMax))
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------
// Event Plane method and Scalar Product method
TVector2 StTriFlowCorrection::calq2Vector(StPicoTrack *track)
{
  const Float_t phi = track->pMom().phi();
  TVector2 q2Vector(0.0,0.0);

  const Float_t q2x = TMath::Cos(2.0*phi);
  const Float_t q2y = TMath::Sin(2.0*phi);
  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

TVector2 StTriFlowCorrection::calq3Vector(StPicoTrack *track)
{
  const Float_t phi = track->pMom().phi();
  TVector2 q3Vector(0.0,0.0);

  const Float_t q3x = TMath::Cos(3.0*phi);
  const Float_t q3y = TMath::Sin(3.0*phi);
  q3Vector.Set(q3x,q3y);

  return q3Vector;
}

Float_t StTriFlowCorrection::getWeight(StPicoTrack *track)
{
  Float_t pt = track->pMom().perp();
  Float_t w;
  if(pt <= TriFlow::mPrimPtWeight)
  {
    w = pt;
  }
  if(pt > TriFlow::mPrimPtWeight)
  {
    w = TriFlow::mPrimPtWeight;
  }

  return w;
}
//---------------------------------------------------------------------------------

TVector2 StTriFlowCorrection::getReCenterPar_East(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t eta_gap, Int_t method)
{
  Float_t mean_qx, mean_qy;
  TVector2 qVector(0.0,0.0);

  TString ProName_x = Form("qx_%s_Vertex_%s_EtaGap_%d_East_%s",mOrder[order].Data(),mVStr[vz_sign].Data(),eta_gap,mMethod[method].Data());
  TString ProName_y = Form("qy_%s_Vertex_%s_EtaGap_%d_East_%s",mOrder[order].Data(),mVStr[vz_sign].Data(),eta_gap,mMethod[method].Data());

  TProfile2D *p_x = (TProfile2D*)mInPutFile->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile->Get(ProName_y.Data());

  mean_qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex,(Double_t)Cent9));
  mean_qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex,(Double_t)Cent9));

  qVector.Set(mean_qx,mean_qy);

  return qVector;
}

//---------------------------------------------------------------------------------

TVector2 StTriFlowCorrection::getReCenterPar_West(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t eta_gap, Int_t method)
{
  Float_t mean_qx, mean_qy;
  TVector2 qVector(0.0,0.0);

  TString ProName_x = Form("qx_%s_Vertex_%s_EtaGap_%d_West_%s",mOrder[order].Data(),mVStr[vz_sign].Data(),eta_gap,mMethod[method].Data());
  TString ProName_y = Form("qy_%s_Vertex_%s_EtaGap_%d_West_%s",mOrder[order].Data(),mVStr[vz_sign].Data(),eta_gap,mMethod[method].Data());

  TProfile2D *p_x = (TProfile2D*)mInPutFile->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile->Get(ProName_y.Data());

  mean_qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex,(Double_t)Cent9));
  mean_qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex,(Double_t)Cent9));

  qVector.Set(mean_qx,mean_qy);

  return qVector;
}

//---------------------------------------------------------------------------------

TVector2 StTriFlowCorrection::getReCenterPar_Full(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t method)
{
  Float_t mean_qx, mean_qy;
  TVector2 qVector(0.0,0.0);

  TString ProName_x = Form("qx_%s_Vertex_%s_Full_%s",mOrder[order].Data(),mVStr[vz_sign].Data(),mMethod[method].Data());
  TString ProName_y = Form("qy_%s_Vertex_%s_Full_%s",mOrder[order].Data(),mVStr[vz_sign].Data(),mMethod[method].Data());

  TProfile2D *p_x = (TProfile2D*)mInPutFile->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile->Get(ProName_y.Data());

  mean_qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex,(Double_t)Cent9));
  mean_qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex,(Double_t)Cent9));

  qVector.Set(mean_qx,mean_qy);

  return qVector;
}

//---------------------------------------------------------------------------------

void StTriFlowCorrection::addTrack_East(StPicoTrack *track, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j)
{
  // Event Plane method
  Float_t w = getWeight(track);
  mQ2Vector_East_EP[j] += w*(calq2Vector(track) - getReCenterPar_East(0,Cent9,RunIndex,i,j,0));
  mQ3Vector_East_EP[j] += w*(calq3Vector(track) - getReCenterPar_East(1,Cent9,RunIndex,i,j,0));

  // Scalor Product method
  mQ2Vector_East_SP[j] += calq2Vector(track) - getReCenterPar_East(0,Cent9,RunIndex,i,j,1);
  mQ3Vector_East_SP[j] += calq3Vector(track) - getReCenterPar_East(1,Cent9,RunIndex,i,j,1);

  mQCounter_East[j]++;
}

//---------------------------------------------------------------------------------

void StTriFlowCorrection::addTrack_West(StPicoTrack *track, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j)
{
  // Event Plane method
  Float_t w = getWeight(track);
  mQ2Vector_West_EP[j] += w*(calq2Vector(track) - getReCenterPar_West(0,Cent9,RunIndex,i,j,0));
  mQ3Vector_West_EP[j] += w*(calq3Vector(track) - getReCenterPar_West(1,Cent9,RunIndex,i,j,0));

  // Scalor Product method
  mQ2Vector_West_SP[j] += calq2Vector(track) - getReCenterPar_West(0,Cent9,RunIndex,i,j,1);
  mQ3Vector_West_SP[j] += calq3Vector(track) - getReCenterPar_West(1,Cent9,RunIndex,i,j,1);

  mQCounter_West[j]++;
}

//---------------------------------------------------------------------------------

void StTriFlowCorrection::addTrack_Full(StPicoTrack *track, Int_t Cent9, Int_t RunIndex, Int_t i)
{
  // Event Plane method
  Float_t w = getWeight(track);
  mQ2Vector_Full_EP += w*(calq2Vector(track) - getReCenterPar_Full(0,Cent9,RunIndex,i,0));
  mQ3Vector_Full_EP += w*(calq3Vector(track) - getReCenterPar_Full(1,Cent9,RunIndex,i,0));

  // Scalor Product method
  mQ2Vector_Full_SP += calq2Vector(track) - getReCenterPar_Full(0,Cent9,RunIndex,i,1);
  mQ3Vector_Full_SP += calq3Vector(track) - getReCenterPar_Full(1,Cent9,RunIndex,i,1);

  mQCounter_Full++;

  Float_t eta = track->pMom().pseudoRapidity();
  if(eta >= 0.0)
  {
    mQCounter_Full_West++;
  }
  if(eta < 0.0)
  {
    mQCounter_Full_East++;
  }
}

void StTriFlowCorrection::addTrack_A(StPicoTrack *track, Int_t Cent9, Int_t RunIndex, Int_t i)
{
  // Event Plane method
  Float_t w = getWeight(track);
  mQ2Vector_A_EP += w*(calq2Vector(track) - getReCenterPar_Full(0,Cent9,RunIndex,i,0));
  mQ3Vector_A_EP += w*(calq3Vector(track) - getReCenterPar_Full(1,Cent9,RunIndex,i,0));

  // Scalor Product method
  mQ2Vector_A_SP += calq2Vector(track) - getReCenterPar_Full(0,Cent9,RunIndex,i,1);
  mQ3Vector_A_SP += calq3Vector(track) - getReCenterPar_Full(1,Cent9,RunIndex,i,1);

  mQCounter_A++;
}

void StTriFlowCorrection::addTrack_B(StPicoTrack *track, Int_t Cent9, Int_t RunIndex, Int_t i)
{
  // Event Plane method
  Float_t w = getWeight(track);
  mQ2Vector_B_EP += w*(calq2Vector(track) - getReCenterPar_Full(0,Cent9,RunIndex,i,0));
  mQ3Vector_B_EP += w*(calq3Vector(track) - getReCenterPar_Full(1,Cent9,RunIndex,i,0));

  // Scalor Product method
  mQ2Vector_B_SP += calq2Vector(track) - getReCenterPar_Full(0,Cent9,RunIndex,i,1);
  mQ3Vector_B_SP += calq3Vector(track) - getReCenterPar_Full(1,Cent9,RunIndex,i,1);

  mQCounter_B++;
}

void StTriFlowCorrection::Randomization()
{
  TRandom3 Ran;
  TVector2 Q2Switch_EP,Q3Switch_EP;
  TVector2 Q2Switch_SP,Q3Switch_SP;
  Int_t CSwitch;
  Ran.SetSeed();
  Float_t ran = Ran.Rndm(); // random number between [0,1]
  if(ran < 0.5)
  {
    // switch Event Plane Q Vector
    Q2Switch_EP = mQ2Vector_A_EP;
    mQ2Vector_A_EP = mQ2Vector_B_EP;
    mQ2Vector_B_EP = Q2Switch_EP;

    Q3Switch_EP = mQ3Vector_A_EP;
    mQ3Vector_A_EP = mQ3Vector_B_EP;
    mQ3Vector_B_EP = Q3Switch_EP;

    // switch Scalar Product Q Vector
    Q2Switch_SP = mQ2Vector_A_SP;
    mQ2Vector_A_SP = mQ2Vector_B_SP;
    mQ2Vector_B_SP = Q2Switch_SP;

    Q3Switch_SP = mQ3Vector_A_SP;
    mQ3Vector_A_SP = mQ3Vector_B_SP;
    mQ3Vector_B_SP = Q3Switch_SP;

    // switch Counter
    CSwitch = mQCounter_A;
    mQCounter_A = mQCounter_B;
    mQCounter_B = CSwitch;
  }
//  cout << "random = " << ran << endl;
//  cout << "mQCounter_A = " << mQCounter_A << endl;
//  cout << "mQCounter_B = " << mQCounter_B << endl;
//  cout << "mQCounter_Full = " << mQCounter_Full << endl;
}

//---------------------------------------------------------------------------------

void StTriFlowCorrection::print(TVector2 vector)
{
  cout << "qx = " << vector.X() << endl;
  cout << "qy = " << vector.Y() << endl;
  cout << endl;
}

//---------------------------------------------------------------------------------

void StTriFlowCorrection::clear()
{
  for(Int_t i = 0; i < 4; i++)
  {
    mQ2Vector_East_EP[i].Set(0.0,0.0);
    mQ3Vector_East_EP[i].Set(0.0,0.0);
    mQ2Vector_East_SP[i].Set(0.0,0.0);
    mQ3Vector_East_SP[i].Set(0.0,0.0);
    mQCounter_East[i] = 0;

    mQ2Vector_West_EP[i].Set(0.0,0.0);
    mQ3Vector_West_EP[i].Set(0.0,0.0);
    mQ2Vector_West_SP[i].Set(0.0,0.0);
    mQ3Vector_West_SP[i].Set(0.0,0.0);
    mQCounter_West[i] = 0;
  }
  
  mQ2Vector_Full_EP.Set(0.0,0.0);
  mQ3Vector_Full_EP.Set(0.0,0.0);
  mQ2Vector_Full_SP.Set(0.0,0.0);
  mQ3Vector_Full_SP.Set(0.0,0.0);
  mQCounter_Full = 0;
  mQCounter_Full_East = 0;
  mQCounter_Full_West = 0;

  mQ2Vector_A_EP.Set(0.0,0.0);
  mQ3Vector_A_EP.Set(0.0,0.0);
  mQ2Vector_A_SP.Set(0.0,0.0);
  mQ3Vector_A_SP.Set(0.0,0.0);
  mQCounter_A = 0;

  mQ2Vector_B_EP.Set(0.0,0.0);
  mQ3Vector_B_EP.Set(0.0,0.0);
  mQ2Vector_B_SP.Set(0.0,0.0);
  mQ3Vector_B_SP.Set(0.0,0.0);
  mQCounter_B = 0;
}

//---------------------------------------------------------------------------------

void StTriFlowCorrection::InitNtuple()
{
  mNtuple = new TNtuple("Ntuple","Ntuple","runId:eventId:RefMult:ZDCx:BBCx:vzVpd:Centrality9:vx:vy:vz:nToFMatched:mQ2X_East_EP_0:mQ2X_East_EP_1:mQ2X_East_EP_2:mQ2X_East_EP_3:mQ2Y_East_EP_0:mQ2Y_East_EP_1:mQ2Y_East_EP_2:mQ2Y_East_EP_3:mQ3X_East_EP_0:mQ3X_East_EP_1:mQ3X_East_EP_2:mQ3X_East_EP_3:mQ3Y_East_EP_0:mQ3Y_East_EP_1:mQ3Y_East_EP_2:mQ3Y_East_EP_3:mQ2X_West_EP_0:mQ2X_West_EP_1:mQ2X_West_EP_2:mQ2X_West_EP_3:mQ2Y_West_EP_0:mQ2Y_West_EP_1:mQ2Y_West_EP_2:mQ2Y_West_EP_3:mQ3X_West_EP_0:mQ3X_West_EP_1:mQ3X_West_EP_2:mQ3X_West_EP_3:mQ3Y_West_EP_0:mQ3Y_West_EP_1:mQ3Y_West_EP_2:mQ3Y_West_EP_3:mQ2X_Full_EP:mQ2Y_Full_EP:mQ3X_Full_EP:mQ3Y_Full_EP:mQ2X_East_SP_0:mQ2X_East_SP_1:mQ2X_East_SP_2:mQ2X_East_SP_3:mQ2Y_East_SP_0:mQ2Y_East_SP_1:mQ2Y_East_SP_2:mQ2Y_East_SP_3:mQ3X_East_SP_0:mQ3X_East_SP_1:mQ3X_East_SP_2:mQ3X_East_SP_3:mQ3Y_East_SP_0:mQ3Y_East_SP_1:mQ3Y_East_SP_2:mQ3Y_East_SP_3:mQ2X_West_SP_0:mQ2X_West_SP_1:mQ2X_West_SP_2:mQ2X_West_SP_3:mQ2Y_West_SP_0:mQ2Y_West_SP_1:mQ2Y_West_SP_2:mQ2Y_West_SP_3:mQ3X_West_SP_0:mQ3X_West_SP_1:mQ3X_West_SP_2:mQ3X_West_SP_3:mQ3Y_West_SP_0:mQ3Y_West_SP_1:mQ3Y_West_SP_2:mQ3Y_West_SP_3:mQ2X_Full_SP:mQ2Y_Full_SP:mQ3X_Full_SP:mQ3Y_Full_SP:mQCounter_East_0:mQCounter_East_1:mQCounter_East_2:mQCounter_East_3:mQCounter_West_0:mQCounter_West_1:mQCounter_West_2:mQCounter_West_3:mQCounter_Full:mQCounter_Full_East:mQCounter_Full_West:runIndex:mQ2X_A_EP:mQ2Y_A_EP:mQ3X_A_EP:mQ3Y_A_EP:mQ2X_A_SP:mQ2Y_A_SP:mQ3X_A_SP:mQ3Y_A_SP:mQCounter_A:mQ2X_B_EP:mQ2Y_B_EP:mQ3X_B_EP:mQ3Y_B_EP:mQ2X_B_SP:mQ2Y_B_SP:mQ3X_B_SP:mQ3Y_B_SP:mQCounter_B");
  mNtuple->SetAutoSave(50000000);
}

//---------------------------------------------------------------------------------

Int_t StTriFlowCorrection::fillNtuple(StPicoDst *pico, Int_t Cent9, Int_t nToFMatched, Int_t runIndex)
{
  StPicoEvent *event = pico->event();
  if(!event)
  {
    return kFALSE;
  }

  Randomization();

  // event information
  mFillNtuple[0]  = (Float_t)event->runId();
  mFillNtuple[1]  = (Float_t)event->eventId();
  mFillNtuple[2]  = (Float_t)event->refMult();
  mFillNtuple[3]  = (Float_t)event->ZDCx();
  mFillNtuple[4]  = (Float_t)event->BBCx();
  mFillNtuple[5]  = (Float_t)event->vzVpd();
  mFillNtuple[6]  = (Float_t)Cent9;
  mFillNtuple[7]  = (Float_t)event->primaryVertex().x();
  mFillNtuple[8]  = (Float_t)event->primaryVertex().y();
  mFillNtuple[9]  = (Float_t)event->primaryVertex().z();
  mFillNtuple[10] = (Float_t)nToFMatched;

  // Q Vector Event Plane method
  // East
  mFillNtuple[11] = (Float_t)mQ2Vector_East_EP[0].X();
  mFillNtuple[12] = (Float_t)mQ2Vector_East_EP[1].X();
  mFillNtuple[13] = (Float_t)mQ2Vector_East_EP[2].X();
  mFillNtuple[14] = (Float_t)mQ2Vector_East_EP[3].X();

  mFillNtuple[15] = (Float_t)mQ2Vector_East_EP[0].Y();
  mFillNtuple[16] = (Float_t)mQ2Vector_East_EP[1].Y();
  mFillNtuple[17] = (Float_t)mQ2Vector_East_EP[2].Y();
  mFillNtuple[18] = (Float_t)mQ2Vector_East_EP[3].Y();

  mFillNtuple[19] = (Float_t)mQ3Vector_East_EP[0].X();
  mFillNtuple[20] = (Float_t)mQ3Vector_East_EP[1].X();
  mFillNtuple[21] = (Float_t)mQ3Vector_East_EP[2].X();
  mFillNtuple[22] = (Float_t)mQ3Vector_East_EP[3].X();

  mFillNtuple[23] = (Float_t)mQ3Vector_East_EP[0].Y();
  mFillNtuple[24] = (Float_t)mQ3Vector_East_EP[1].Y();
  mFillNtuple[25] = (Float_t)mQ3Vector_East_EP[2].Y();
  mFillNtuple[26] = (Float_t)mQ3Vector_East_EP[3].Y();
  // West
  mFillNtuple[27] = (Float_t)mQ2Vector_West_EP[0].X();
  mFillNtuple[28] = (Float_t)mQ2Vector_West_EP[1].X();
  mFillNtuple[29] = (Float_t)mQ2Vector_West_EP[2].X();
  mFillNtuple[30] = (Float_t)mQ2Vector_West_EP[3].X();

  mFillNtuple[31] = (Float_t)mQ2Vector_West_EP[0].Y();
  mFillNtuple[32] = (Float_t)mQ2Vector_West_EP[1].Y();
  mFillNtuple[33] = (Float_t)mQ2Vector_West_EP[2].Y();
  mFillNtuple[34] = (Float_t)mQ2Vector_West_EP[3].Y();

  mFillNtuple[35] = (Float_t)mQ3Vector_West_EP[0].X();
  mFillNtuple[36] = (Float_t)mQ3Vector_West_EP[1].X();
  mFillNtuple[37] = (Float_t)mQ3Vector_West_EP[2].X();
  mFillNtuple[38] = (Float_t)mQ3Vector_West_EP[3].X();

  mFillNtuple[39] = (Float_t)mQ3Vector_West_EP[0].Y();
  mFillNtuple[40] = (Float_t)mQ3Vector_West_EP[1].Y();
  mFillNtuple[41] = (Float_t)mQ3Vector_West_EP[2].Y();
  mFillNtuple[42] = (Float_t)mQ3Vector_West_EP[3].Y();
  // Full
  mFillNtuple[43] = (Float_t)mQ2Vector_Full_EP.X();
  mFillNtuple[44] = (Float_t)mQ2Vector_Full_EP.Y();
  mFillNtuple[45] = (Float_t)mQ3Vector_Full_EP.X();
  mFillNtuple[46] = (Float_t)mQ3Vector_Full_EP.Y();

  // Q Vector Scalor Product method
  // East
  mFillNtuple[47] = (Float_t)mQ2Vector_East_SP[0].X();
  mFillNtuple[48] = (Float_t)mQ2Vector_East_SP[1].X();
  mFillNtuple[49] = (Float_t)mQ2Vector_East_SP[2].X();
  mFillNtuple[50] = (Float_t)mQ2Vector_East_SP[3].X();

  mFillNtuple[51] = (Float_t)mQ2Vector_East_SP[0].Y();
  mFillNtuple[52] = (Float_t)mQ2Vector_East_SP[1].Y();
  mFillNtuple[53] = (Float_t)mQ2Vector_East_SP[2].Y();
  mFillNtuple[54] = (Float_t)mQ2Vector_East_SP[3].Y();

  mFillNtuple[55] = (Float_t)mQ3Vector_East_SP[0].X();
  mFillNtuple[56] = (Float_t)mQ3Vector_East_SP[1].X();
  mFillNtuple[57] = (Float_t)mQ3Vector_East_SP[2].X();
  mFillNtuple[58] = (Float_t)mQ3Vector_East_SP[3].X();

  mFillNtuple[59] = (Float_t)mQ3Vector_East_SP[0].Y();
  mFillNtuple[60] = (Float_t)mQ3Vector_East_SP[1].Y();
  mFillNtuple[61] = (Float_t)mQ3Vector_East_SP[2].Y();
  mFillNtuple[62] = (Float_t)mQ3Vector_East_SP[3].Y();
  // West
  mFillNtuple[63] = (Float_t)mQ2Vector_West_SP[0].X();
  mFillNtuple[64] = (Float_t)mQ2Vector_West_SP[1].X();
  mFillNtuple[65] = (Float_t)mQ2Vector_West_SP[2].X();
  mFillNtuple[66] = (Float_t)mQ2Vector_West_SP[3].X();

  mFillNtuple[67] = (Float_t)mQ2Vector_West_SP[0].Y();
  mFillNtuple[68] = (Float_t)mQ2Vector_West_SP[1].Y();
  mFillNtuple[69] = (Float_t)mQ2Vector_West_SP[2].Y();
  mFillNtuple[70] = (Float_t)mQ2Vector_West_SP[3].Y();

  mFillNtuple[71] = (Float_t)mQ3Vector_West_SP[0].X();
  mFillNtuple[72] = (Float_t)mQ3Vector_West_SP[1].X();
  mFillNtuple[73] = (Float_t)mQ3Vector_West_SP[2].X();
  mFillNtuple[74] = (Float_t)mQ3Vector_West_SP[3].X();

  mFillNtuple[75] = (Float_t)mQ3Vector_West_SP[0].Y();
  mFillNtuple[76] = (Float_t)mQ3Vector_West_SP[1].Y();
  mFillNtuple[77] = (Float_t)mQ3Vector_West_SP[2].Y();
  mFillNtuple[78] = (Float_t)mQ3Vector_West_SP[3].Y();
  // Full
  mFillNtuple[79] = (Float_t)mQ2Vector_Full_SP.X();
  mFillNtuple[80] = (Float_t)mQ2Vector_Full_SP.Y();
  mFillNtuple[81] = (Float_t)mQ3Vector_Full_SP.X();
  mFillNtuple[82] = (Float_t)mQ3Vector_Full_SP.Y();

  // Counter
  // East
  mFillNtuple[83] = (Float_t)mQCounter_East[0];
  mFillNtuple[84] = (Float_t)mQCounter_East[1];
  mFillNtuple[85] = (Float_t)mQCounter_East[2];
  mFillNtuple[86] = (Float_t)mQCounter_East[3];
  // West
  mFillNtuple[87] = (Float_t)mQCounter_West[0];
  mFillNtuple[88] = (Float_t)mQCounter_West[1];
  mFillNtuple[89] = (Float_t)mQCounter_West[2];
  mFillNtuple[90] = (Float_t)mQCounter_West[3];
  // Full
  mFillNtuple[91] = (Float_t)mQCounter_Full;
  mFillNtuple[92] = (Float_t)mQCounter_Full_East;
  mFillNtuple[93] = (Float_t)mQCounter_Full_West;

  // runIndex
  mFillNtuple[94] = (Float_t)runIndex;

  // random sub
  mFillNtuple[95]  = (Float_t)mQ2Vector_A_EP.X();
  mFillNtuple[96]  = (Float_t)mQ2Vector_A_EP.Y();
  mFillNtuple[97]  = (Float_t)mQ3Vector_A_EP.X();
  mFillNtuple[98]  = (Float_t)mQ3Vector_A_EP.Y();
  mFillNtuple[99]  = (Float_t)mQ2Vector_A_SP.X();
  mFillNtuple[100] = (Float_t)mQ2Vector_A_SP.Y();
  mFillNtuple[101] = (Float_t)mQ3Vector_A_SP.X();
  mFillNtuple[102] = (Float_t)mQ3Vector_A_SP.Y();
  mFillNtuple[103] = (Float_t)mQCounter_A;

  mFillNtuple[104] = (Float_t)mQ2Vector_B_EP.X();
  mFillNtuple[105] = (Float_t)mQ2Vector_B_EP.Y();
  mFillNtuple[106] = (Float_t)mQ3Vector_B_EP.X();
  mFillNtuple[107] = (Float_t)mQ3Vector_B_EP.Y();
  mFillNtuple[108] = (Float_t)mQ2Vector_B_SP.X();
  mFillNtuple[109] = (Float_t)mQ2Vector_B_SP.Y();
  mFillNtuple[110] = (Float_t)mQ3Vector_B_SP.X();
  mFillNtuple[111] = (Float_t)mQ3Vector_B_SP.Y();
  mFillNtuple[112] = (Float_t)mQCounter_B;

  mNtuple->Fill(mFillNtuple);

  return kTRUE;
}

//---------------------------------------------------------------------------------

void StTriFlowCorrection::writeNtuple()
{
  mNtuple->Write();
}

//---------------------------------------------------------------------------------

bool StTriFlowCorrection::passTrackEtaNumCut(Int_t j)
{
  if(!(mQCounter_East[j] > TriFlow::mTrackMin && mQCounter_West[j] > TriFlow::mTrackMin))
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------

bool StTriFlowCorrection::passTrackFullNumCut()
{
  if(!(mQCounter_Full > TriFlow::mTrackMin_Full && mQCounter_Full_East > 0 && mQCounter_Full_West > 0))
  {
    return kFALSE;
  }
  
  return kTRUE;
}

//---------------------------------------------------------------------------------
// Event Plane method
// 2nd
TVector2 StTriFlowCorrection::calPsi2_East_EP(Int_t j, Int_t k) // j = eta_gap, k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ2Vector_East_EP[j].X();
  Float_t Qy = mQ2Vector_East_EP[j].Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/2.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder2[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder2[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

TVector2 StTriFlowCorrection::calPsi2_West_EP(Int_t j, Int_t k) // j = eta_gap, k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ2Vector_West_EP[j].X();
  Float_t Qy = mQ2Vector_West_EP[j].Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/2.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder2[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder2[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

TVector2 StTriFlowCorrection::calPsi2_Full_EP(Int_t k) // k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ2Vector_Full_EP.X();
  Float_t Qy = mQ2Vector_Full_EP.Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/2.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder2[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder2[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

//---------------------------------------------------------------------------------
// 3rd
TVector2 StTriFlowCorrection::calPsi3_East_EP(Int_t j, Int_t k) // j = eta_gap, k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ3Vector_East_EP[j].X();
  Float_t Qy = mQ3Vector_East_EP[j].Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/3.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder3[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder3[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

TVector2 StTriFlowCorrection::calPsi3_West_EP(Int_t j, Int_t k) // j = eta_gap, k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ3Vector_West_EP[j].X();
  Float_t Qy = mQ3Vector_West_EP[j].Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/3.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder3[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder3[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

TVector2 StTriFlowCorrection::calPsi3_Full_EP(Int_t k) // k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ3Vector_Full_EP.X();
  Float_t Qy = mQ3Vector_Full_EP.Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/3.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder3[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder3[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

//---------------------------------------------------------------------------------
// Scalor Product method
// 2nd
TVector2 StTriFlowCorrection::calPsi2_East_SP(Int_t j, Int_t k) // j = eta_gap, k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ2Vector_East_SP[j].X();
  Float_t Qy = mQ2Vector_East_SP[j].Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/2.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder2[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder2[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

TVector2 StTriFlowCorrection::calPsi2_West_SP(Int_t j, Int_t k) // j = eta_gap, k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ2Vector_West_SP[j].X();
  Float_t Qy = mQ2Vector_West_SP[j].Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/2.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder2[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder2[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

TVector2 StTriFlowCorrection::calPsi2_Full_SP(Int_t k) // k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ2Vector_Full_SP.X();
  Float_t Qy = mQ2Vector_Full_SP.Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/2.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder2[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder2[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

//---------------------------------------------------------------------------------
// 3rd
TVector2 StTriFlowCorrection::calPsi3_East_SP(Int_t j, Int_t k) // j = eta_gap, k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ3Vector_East_SP[j].X();
  Float_t Qy = mQ3Vector_East_SP[j].Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/3.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder3[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder3[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

TVector2 StTriFlowCorrection::calPsi3_West_SP(Int_t j, Int_t k) // j = eta_gap, k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ3Vector_West_SP[j].X();
  Float_t Qy = mQ3Vector_West_SP[j].Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/3.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder3[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder3[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

TVector2 StTriFlowCorrection::calPsi3_Full_SP(Int_t k) // k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ3Vector_Full_SP.X();
  Float_t Qy = mQ3Vector_Full_SP.Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/3.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder3[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder3[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

//---------------------------------------------------------------------------------

Float_t StTriFlowCorrection::AngleShift(Float_t Psi_raw, Float_t order)
{
  Float_t Psi_Corr = Psi_raw;
  if(Psi_raw > TMath::Pi()/order)
  {
    Psi_Corr = Psi_raw - 2.0*TMath::Pi()/order;
  }
  if(Psi_raw < -1.0*TMath::Pi()/order)
  {
    Psi_Corr = Psi_raw + 2.0*TMath::Pi()/order;
  }

  return Psi_Corr;
}

//---------------------------------------------------------------------------------

// Event Plane method
// 2nd
Float_t StTriFlowCorrection::calShiftAngle2East_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Vector_East_EP[eta_gap].Y(),mQ2Vector_East_EP[eta_gap].X())/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[vz_sign].Data(),eta_gap,k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[vz_sign].Data(),eta_gap,k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,2.0);

  return Psi_Shift;
}

Float_t StTriFlowCorrection::calShiftAngle2West_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Vector_West_EP[eta_gap].Y(),mQ2Vector_West_EP[eta_gap].X())/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,2.0);

  return Psi_Shift;
}

Float_t StTriFlowCorrection::calShiftAngle2Full_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, StPicoTrack *track)
{
  TVector2 QVector_sub(0.0,0.0);
  if(passTrackFull(track))
  {
    Float_t w = getWeight(track);
    QVector_sub = mQ2Vector_Full_EP - w*(calq2Vector(track) - getReCenterPar_Full(0,Cent9,runIndex,vz_sign,0));
//    QVector_sub = mQ2Vector_Full_EP; 
  }
  else
  {
    QVector_sub = mQ2Vector_Full_EP;
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

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,2.0);

  return Psi_Shift;
}

// 3rd 
Float_t StTriFlowCorrection::calShiftAngle3East_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ3Vector_East_EP[eta_gap].Y(),mQ3Vector_East_EP[eta_gap].X())/3.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[vz_sign].Data(),eta_gap,k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[vz_sign].Data(),eta_gap,k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder3[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,3.0);

  return Psi_Shift;
}

Float_t StTriFlowCorrection::calShiftAngle3West_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ3Vector_West_EP[eta_gap].Y(),mQ3Vector_West_EP[eta_gap].X())/3.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder3[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,3.0);

  return Psi_Shift;
}

Float_t StTriFlowCorrection::calShiftAngle3Full_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, StPicoTrack *track)
{
  TVector2 QVector_sub(0.0,0.0);
  if(passTrackFull(track))
  {
    Float_t w = getWeight(track);
    QVector_sub = mQ3Vector_Full_EP - w*(calq3Vector(track) - getReCenterPar_Full(1,Cent9,runIndex,vz_sign,0));
//    QVector_sub = mQ3Vector_Full_EP;
  }
  else
  {
    QVector_sub = mQ3Vector_Full_EP;
  }
  Float_t Psi_ReCenter = TMath::ATan2(QVector_sub.Y(),QVector_sub.X())/3.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi3_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi3_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder3[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,3.0);

  return Psi_Shift;
}

//---------------------------------------------------------------------------------
// Scalor Product method
// 2nd
TVector2 StTriFlowCorrection::calQVector2East_SP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Vector_East_SP[eta_gap].Y(),mQ2Vector_East_SP[eta_gap].X())/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;
  TVector2 QVector(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_SP",mVStr[vz_sign].Data(),eta_gap,k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_SP",mVStr[vz_sign].Data(),eta_gap,k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,2.0);
  Float_t Qx = mQ2Vector_East_SP[eta_gap].Mod()*TMath::Cos(2.0*Psi_Shift);
  Float_t Qy = mQ2Vector_East_SP[eta_gap].Mod()*TMath::Sin(2.0*Psi_Shift);
  QVector.Set(Qx,Qy);

  return QVector;
}

TVector2 StTriFlowCorrection::calQVector2West_SP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Vector_West_SP[eta_gap].Y(),mQ2Vector_West_SP[eta_gap].X())/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;
  TVector2 QVector(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_SP",mVStr[vz_sign].Data(),eta_gap,k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_SP",mVStr[vz_sign].Data(),eta_gap,k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,2.0);
  Float_t Qx = mQ2Vector_West_SP[eta_gap].Mod()*TMath::Cos(2.0*Psi_Shift);
  Float_t Qy = mQ2Vector_West_SP[eta_gap].Mod()*TMath::Sin(2.0*Psi_Shift);
  QVector.Set(Qx,Qy);

  return QVector;
}

TVector2 StTriFlowCorrection::calQVector2Full_SP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, StPicoTrack *track)
{
  TVector2 QVector_sub(0.0,0.0);
  if(passTrackFull(track))
  {
    QVector_sub = mQ2Vector_Full_SP - (calq2Vector(track) - getReCenterPar_Full(0,Cent9,runIndex,vz_sign,1));
//    QVector_sub = mQ2Vector_Full_SP;
  }
  else
  {
    QVector_sub = mQ2Vector_Full_SP;
  }
  Float_t Psi_ReCenter = TMath::ATan2(QVector_sub.Y(),QVector_sub.X())/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;
  TVector2 QVector(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_Full_SP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_Full_SP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,2.0);
  Float_t Qx = QVector_sub.Mod()*TMath::Cos(2.0*Psi_Shift);
  Float_t Qy = QVector_sub.Mod()*TMath::Sin(2.0*Psi_Shift);
  QVector.Set(Qx,Qy);

  return QVector;
}

// 3rd 
TVector2 StTriFlowCorrection::calQVector3East_SP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ3Vector_East_SP[eta_gap].Y(),mQ3Vector_East_SP[eta_gap].X())/3.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;
  TVector2 QVector(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_SP",mVStr[vz_sign].Data(),eta_gap,k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_SP",mVStr[vz_sign].Data(),eta_gap,k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder3[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,3.0);
  Float_t Qx = mQ3Vector_East_SP[eta_gap].Mod()*TMath::Cos(3.0*Psi_Shift);
  Float_t Qy = mQ3Vector_East_SP[eta_gap].Mod()*TMath::Sin(3.0*Psi_Shift);
  QVector.Set(Qx,Qy);

  return QVector;
}

TVector2 StTriFlowCorrection::calQVector3West_SP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ3Vector_West_SP[eta_gap].Y(),mQ3Vector_West_SP[eta_gap].X())/3.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;
  TVector2 QVector(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_SP",mVStr[vz_sign].Data(),eta_gap,k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_SP",mVStr[vz_sign].Data(),eta_gap,k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder3[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,3.0);
  Float_t Qx = mQ3Vector_West_SP[eta_gap].Mod()*TMath::Cos(3.0*Psi_Shift);
  Float_t Qy = mQ3Vector_West_SP[eta_gap].Mod()*TMath::Sin(3.0*Psi_Shift);
  QVector.Set(Qx,Qy);

  return QVector;
}

TVector2 StTriFlowCorrection::calQVector3Full_SP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, StPicoTrack *track)
{
  TVector2 QVector_sub(0.0,0.0);
  if(passTrackFull(track))
  {
    QVector_sub = mQ3Vector_Full_SP - (calq3Vector(track) - getReCenterPar_Full(1,Cent9,runIndex,vz_sign,1));
//    QVector_sub = mQ3Vector_Full_SP;
  }
  else
  {
    QVector_sub = mQ3Vector_Full_SP;
  }
  Float_t Psi_ReCenter = TMath::ATan2(QVector_sub.Y(),QVector_sub.X())/3.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;
  TVector2 QVector(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi3_Vertex_%s_Order_%d_Full_SP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi3_Vertex_%s_Order_%d_Full_SP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder3[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,3.0);
  Float_t Qx = QVector_sub.Mod()*TMath::Cos(3.0*Psi_Shift);
  Float_t Qy = QVector_sub.Mod()*TMath::Sin(3.0*Psi_Shift);
  QVector.Set(Qx,Qy);

  return QVector;
}

//---------------------------------------------------------------------------------
//Event Plane method
Float_t StTriFlowCorrection::getResolution2_EP(Int_t Cent9, Int_t eta_gap)
{
  TString ProName = Form("Res2_EtaGap_%d_EP",eta_gap);
  TProfile *p_res2 = (TProfile*)mInPutFile_Res->Get(ProName.Data());
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

Float_t StTriFlowCorrection::getResolution3_EP(Int_t Cent9, Int_t eta_gap)
{
  TString ProName = Form("Res3_EtaGap_%d_EP",eta_gap);
  TProfile *p_res3 = (TProfile*)mInPutFile_Res->Get(ProName.Data());
  Float_t Res_raw = p_res3->GetBinContent(p_res3->FindBin(Cent9));
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

//Scalar Product method
Float_t StTriFlowCorrection::getResolution2_SP(Int_t Cent9, Int_t eta_gap)
{
  TString ProName = Form("Res2_EtaGap_%d_SP",eta_gap);
  TProfile *p_res2 = (TProfile*)mInPutFile_Res->Get(ProName.Data());
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

Float_t StTriFlowCorrection::getResolution3_SP(Int_t Cent9, Int_t eta_gap)
{
  TString ProName = Form("Res3_EtaGap_%d_SP",eta_gap);
  TProfile *p_res3 = (TProfile*)mInPutFile_Res->Get(ProName.Data());
  Float_t Res_raw = p_res3->GetBinContent(p_res3->FindBin(Cent9));
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

//---------------------------------------------------------------------------------
//Event Plane method
Float_t StTriFlowCorrection::getResolution2_Full_EP(Int_t Cent9)
{
  TString ProName = "Res2_Ran_EP";
  TProfile *p_res2 = (TProfile*)mInPutFile_Res->Get(ProName.Data());
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

Float_t StTriFlowCorrection::getResolution3_Full_EP(Int_t Cent9)
{
  TString ProName = "Res3_Ran_EP";
  TProfile *p_res3 = (TProfile*)mInPutFile_Res->Get(ProName.Data());
  Float_t Res_raw = p_res3->GetBinContent(p_res3->FindBin(Cent9));
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

// Scalar Product method
Float_t StTriFlowCorrection::getResolution2_Full_SP(Int_t Cent9)
{
  TString ProName = "Res2_Ran_SP";
  TProfile *p_res2 = (TProfile*)mInPutFile_Res->Get(ProName.Data());
  Float_t Res_raw = p_res2->GetBinContent(p_res2->FindBin(Cent9));
  if(Res_raw <= 0)
  {
    return -999.9;
  }
  else
  {
    Float_t Res = 2.0*TMath::Sqrt(Res_raw);
    return Res;
  }
}

Float_t StTriFlowCorrection::getResolution3_Full_SP(Int_t Cent9)
{
  TString ProName = "Res3_Ran_SP";
  TProfile *p_res3 = (TProfile*)mInPutFile_Res->Get(ProName.Data());
  Float_t Res_raw = p_res3->GetBinContent(p_res3->FindBin(Cent9));
  if(Res_raw <= 0)
  {
    return -999.9;
  }
  else
  {
    Float_t Res = 2.0*TMath::Sqrt(Res_raw);
    return Res;
  }
}
//---------------------------------------------------------------------------------
TVector2 StTriFlowCorrection::getQVector(Int_t j, Int_t k, Int_t l) // 0 = eta_gap, 1 = flow type, 2 = east/west
{
  if(k == 0 && l == 0) return mQ2Vector_East_EP[j];
  if(k == 0 && l == 1) return mQ2Vector_West_EP[j];
  if(k == 1 && l == 0) return mQ3Vector_East_EP[j];
  if(k == 1 && l == 1) return mQ3Vector_West_EP[j];
}

Int_t StTriFlowCorrection::getNumTrack(Int_t j, Int_t l) // 0 = eta_gap, 1 = east/west
{
  if(l == 0) return mQCounter_East[j];
  if(l == 1) return mQCounter_West[j];
}
//---------------------------------------------------------------------------------

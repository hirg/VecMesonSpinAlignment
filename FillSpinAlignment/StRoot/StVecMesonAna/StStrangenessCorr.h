#ifndef StStrangenessCorr_h
#define StStrangenessCorr_h

#include "TObject.h"
#include "TVector2.h"
#include "TString.h"

class TLorentzVector;
class TProfile2D;
class TFile;

class StStrangenessCorr : public TObject
{
  public:
    StStrangenessCorr();
    ~StStrangenessCorr();

    // ReCenter Correction
    void InitReCenterCorrection(Int_t);

    TVector2 getReCenterPar_East(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t eta_gap, Int_t method); // order: 0 = 2nd, 1 = 3rd || method: 0 = EP, 1 = SP
    TVector2 getReCenterPar_West(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t eta_gap, Int_t method); // order: 0 = 2nd, 1 = 3rd || method: 0 = EP, 1 = SP
 
    TVector2 calq2Vector(TLorentzVector);
    TVector2 calq3Vector(TLorentzVector);
    Float_t getWeight(TLorentzVector);

    // Shift Correction
    void InitShiftCorrection(Int_t);
    bool passTrackNumCut(Int_t, Int_t); // Num of Tracks East, Num of Tracks West

    Float_t AngleShift(Float_t Psi_raw, Float_t order);
    Float_t calShiftAngle2East_EP(TVector2, Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap);
    Float_t calShiftAngle2West_EP(TVector2, Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap);

    Float_t calShiftAngle3East_EP(TVector2, Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap);
    Float_t calShiftAngle3West_EP(TVector2, Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap);

    // Resolution Correction
    Float_t getResolution2_EP(Int_t, Int_t); // centrality, eta_gap
    Float_t getResolution3_EP(Int_t, Int_t); // centrality, eta_gap

  private:
    TFile *mInPutFile_ReCenter; // input file for ReCenter Correction
    TFile *mInPutFile_Shift; // input file for ReCenter Correction
    TFile *mInPutFile_Res; // input file for ReCenter Correction

    static TString mVStr[2];
    static TString mOrder[2];
    static TString mMethod[2];

  ClassDef(StStrangenessCorr,1)
};
#endif

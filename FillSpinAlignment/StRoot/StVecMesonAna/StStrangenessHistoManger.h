#ifndef StStrangenessHistoManger_h
#define StStrangenessHistoManger_h

#include "StMessMgr.h"
#include <map>

class TH1F;

typedef std::map<TString,TH1F*> TH1FMap;

class StStrangenessHistoManger
{
  public:
    StStrangenessHistoManger();
    ~StStrangenessHistoManger();

    void Init(Int_t X_flag, Int_t mode);
    void Fill(Float_t pt, Int_t Cent9, Int_t eta_gap, Float_t CosThetaStar, Float_t Res2, Float_t Mass2, Double_t reweight, Int_t X_flag, Int_t mode);
    void Write(Int_t X_flag, Int_t mode);

  private:
    // flow analysis
    // 0 = pt bin
    // 1 = centrality: 0 = 20-60%(0-70%), 1 = 0-10%, 2 = 10-40%, 3 = 40-80%(40-70%)
    // 2 = eta_gap
    // 3 = phi - Psi
    TH1FMap h_mMass2;

    // raw pt spectra
    // 0 = pt bin
    // 1 = pt bin finer
    // 2 = centrality: 0 = 20-60%(0-70%), 1 = 0-10%, 2 = 10-40%, 3 = 40-80%(40-70%)
    // 3 = eta_gap
    TH1FMap h_mMass_Spec;

    // event plane resolution correction
    // 0 = centrality
    // 1 = eta_gap
    TH1FMap h_mMass_Yields;

  ClassDef(StStrangenessHistoManger,1)
};
#endif

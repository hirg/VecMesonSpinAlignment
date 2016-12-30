#ifndef StVecMesonHistoManger_h
#define StVecMesonHistoManger_h

#include "StMessMgr.h"
#include <map>

class TH1F;

typedef std::map<TString,TH1F*> TH1FMap;

class StVecMesonHistoManger
{
  public:
    StVecMesonHistoManger();
    ~StVecMesonHistoManger();

    void Init(Int_t X_flag, Int_t mode);
    void Fill(Float_t pt, Int_t Cent9, Float_t CosThetaStar, Float_t Res2, Float_t Mass2, Double_t reweight, Int_t X_flag, Int_t mode);
    void Write(Int_t X_flag, Int_t mode);

  private:
    // spin alignment analysis
    // 0 = pt bin
    // 1 = centrality: 0 = 20%-60%, 1 = 0-80%, 2 = 0-10%, 3 = 10-40%, 4 = 40-80%
    // 2 = cos(theta*)
    TH1FMap h_mMass2;

    // raw pt spectra
    // 0 = pt bin
    // 1 = pt bin finer
    // 2 = centrality: 0 = 20%-60%, 1 = 0-80%, 2 = 0-10%, 3 = 10-40%, 4 = 40-80%
    TH1FMap h_mMass_Spec;

    // event plane resolution correction
    // centrality: 0 = 20%-60%, 1 = 0-80%, 2 = 0-10%, 3 = 10-40%, 4 = 40-80%
    TH1FMap h_mMass_Yields;

  ClassDef(StVecMesonHistoManger,1)
};
#endif

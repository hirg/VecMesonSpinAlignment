#ifndef StZdcSmdHistoManger_h
#define StZdcSmdHistoManger_h

#include "StMessMgr.h"
#include <map>

class TH1F;

typedef std::map<TString,TH1F*> TH1FMap;

class StZdcSmdHistoManger
{
  public:
    StZdcSmdHistoManger();
    ~StZdcSmdHistoManger();

    void Init(int X_flag, int mode);
    void Fill(float pt, int Cent9, float CosThetaStar, float Res, float Mass2, double reweight, int X_flag, int mode);
    void Write(int X_flag, int mode);

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

  ClassDef(StZdcSmdHistoManger,1)
};
#endif

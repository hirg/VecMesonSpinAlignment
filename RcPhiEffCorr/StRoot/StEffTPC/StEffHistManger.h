#ifndef StEffHistManger_h
#define StEffHistManger_h
#include "TObject.h"
#include "StEffCons.h"
#include <map>
#include <string>
#include "TVector3.h"

class TH1D;
class TH2D;
class TH3D;

typedef std::map<std::string,TH1D*> TH1DMap;

class StEffHistManger : public TObject
{
  public:
    StEffHistManger();
    virtual ~StEffHistManger();
    void InitHist();
    void InitQACos();
    void InitQAEP();
    void FillHistMc(int,float,float,float,float);
    void FillHistRc(int,float,float,float,float);
    void FillQACosMc(int,McVecMeson,McDecayDau,McDecayDau);
    void FillQACosRc(int,McVecMeson,McDecayDau,McDecayDau);
    void FillQAEPMc(int,float,float);
    void FillQAEPRc(int,float,float);
    void CalEfficiency();
    void CalEffPtEtaPhi();
    void CalEffCosThetaStar();
    void CalEffEP();
    TH1D* CalEffError(TH1D*,TH1D*,std::string);
    TVector3 CalBoostedVector(McDecayDau,McVecMeson);
    void WriteHist();

  private:
    TH3D *h_mMcTracks[10]; // pt, eta, phi distribution as a function of centrality, centrality = 9 is for miniBias
    TH3D *h_mRcTracks[10];

    TH1D *h_mMcEffPt[10]; // pt distritbution as a function of centrality
    TH1D *h_mRcEffPt[10];
    TH1D *h_mEffPt[10];

    TH1D *h_mMcEffEta[10]; // pt distritbution as a function of centrality and eta
    TH1D *h_mRcEffEta[10];
    TH1D *h_mEffEta[10];

    TH1D *h_mMcEffPhi[10]; // pt distritbution as a function of centrality and phi
    TH1D *h_mRcEffPhi[10];
    TH1D *h_mEffPhi[10];

    TH1D *h_FrameEta;
    TH1D *h_FramePhi;

    TH1DMap h_mMcEffPEP;
    TH1DMap h_mRcEffPEP; // efficiency as a fucntion of centrality, pt, eta and phi
    TH1DMap h_mEffPEP; // efficiency as a fucntion of centrality, pt, eta and phi

    TH1D *h_mMcEffCos[10][Efficiency::BinPt]; // efficiency vs CosThetaStar as a function of centrality and pt
    TH1D *h_mRcEffCos[10][Efficiency::BinPt];
    TH1D *h_mEffCos[10][Efficiency::BinPt];

    TH1D *h_mMcKpCosStar[10][Efficiency::BinPt];
    TH1D *h_mMcKpCos[10][Efficiency::BinPt];
    TH2D *h_mMcKpCosTheta[10][Efficiency::BinPt];

    TH1D *h_mMcKmCosStar[10][Efficiency::BinPt];
    TH1D *h_mMcKmCos[10][Efficiency::BinPt];
    TH2D *h_mMcKmCosTheta[10][Efficiency::BinPt];

    TH2D *h_mMcKpKmCos[10][Efficiency::BinPt];

    TH1D *h_mRcKpCosStar[10][Efficiency::BinPt];
    TH1D *h_mRcKpCos[10][Efficiency::BinPt];
    TH2D *h_mRcKpCosTheta[10][Efficiency::BinPt];

    TH1D *h_mRcKmCosStar[10][Efficiency::BinPt];
    TH1D *h_mRcKmCos[10][Efficiency::BinPt];
    TH2D *h_mRcKmCosTheta[10][Efficiency::BinPt];

    TH2D *h_mRcKpKmCos[10][Efficiency::BinPt];

    TH2D *h_mMcEP[10]; // pt, phi-Psi2 distribution as a function of centrality, centrality = 9 is for miniBias
    TH2D *h_mRcEP[10];
    TH1D *h_mEffEP[10][Efficiency::BinPt];

    int flag_eff;
    int flag_eff_PtEtaPhi;
    int flag_eff_Cos;
    int flag_QA_Cos;
    int flag_QA_EP;

  ClassDef(StEffHistManger,1)
};

#endif

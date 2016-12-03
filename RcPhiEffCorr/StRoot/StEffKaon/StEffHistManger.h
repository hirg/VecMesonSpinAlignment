#ifndef StEffHistManger_h
#define StEffHistManger_h
#include "TObject.h"
#include "../../../Utility/StSpinAlignmentCons.h"
#include <map>
#include <string>

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
    void FillHistMc(int,float,float,float);
    void FillHistRc(int,float,float,float);
    void FillHistPt(int,float,float,float);
    void CalEfficiency();
    void CalEffPtEtaPhi();
    TH1D* CalEffError(TH1D*,TH1D*,std::string);
    void WriteHist();

  private:
    TH3D *h_mMcTracks[10]; // pt, eta, phi distribution as a function of centrality, centrality = 9 is for miniBias
    TH3D *h_mRcTracks[10];
    TH2D *h_mPtGl[10]; // McPt vs. gRcPt
    TH2D *h_mPtPr[10]; // McPt vs. pRcPt

    TH1D *h_mMcEffPt[10]; // pt distritbution as a function of centrality
    TH1D *h_mRcEffPt[10];
    TH1D *h_mEffPt[10];

    TH1D *h_mMcEffEta[10]; // eta distritbution as a function of centrality
    TH1D *h_mRcEffEta[10];
    TH1D *h_mEffEta[10];

    TH1D *h_mMcEffPhi[10]; // phi distritbution as a function of centrality
    TH1D *h_mRcEffPhi[10];
    TH1D *h_mEffPhi[10];

    TH1D *h_FrameEta;
    TH1D *h_FramePhi;

    TH1DMap h_mMcEffPEP;  
    TH1DMap h_mRcEffPEP; // efficiency as a fucntion of centrality, pt, eta and phi
    TH1DMap h_mEffPEP; // efficiency as a fucntion of centrality, pt, eta and phi

    int flag_eff;
    int flag_eff_PtEtaPhi;

  ClassDef(StEffHistManger,1)
};

#endif

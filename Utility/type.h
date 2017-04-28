#include <string>
#include <map>
#include <vector>
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"

#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooBreitWigner.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"

typedef std::map<std::string,TH1F*> TH1FMap;
typedef std::map<std::string,TH1D*> TH1DMap;
typedef std::map<std::string,TProfile*> TProMap;
typedef std::map<std::string,std::vector<float> > vecFMap;
typedef std::map<std::string,std::vector<double> > vecDMap;
typedef std::map<std::string,TGraphAsymmErrors*> TGraMap;
typedef std::map<std::string,TF1*> TF1Map;

typedef std::map<std::string,RooDataHist*> RooHistMap;
typedef std::map<std::string,RooRealVar*> RooVarMap;
typedef std::map<std::string,RooBreitWigner*> RooBWMap;
typedef std::map<std::string,RooGenericPdf*> RooGePdfMap;
typedef std::map<std::string,RooAddPdf*> RooPdfMap;

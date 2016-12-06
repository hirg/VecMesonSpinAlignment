#include <string>
#include <map>
#include <vector>
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"

typedef std::map<std::string,TH1F*> TH1FMap;
typedef std::map<std::string,TH1D*> TH1DMap;
typedef std::map<std::string,TProfile*> TProMap;
typedef std::map<std::string,std::vector<float> > vecFMap;
typedef std::map<std::string,TGraphAsymmErrors*> TGraMap;


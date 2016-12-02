#include <string>
#include <map>
#include <vector>
#include "TH1F.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"

typedef std::map<std::string,TH1F*> TH1FMap;
typedef std::map<string,TProfile*> TProMap;
typedef std::map<string,std::vector<float> > vecFMap;
typedef std::map<TString,TGraphAsymmErrors*> TGraMap;


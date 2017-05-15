#include "StRoot/StZdcSmdAna/StZdcSmdCorr.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/funtions.h"
#include "TMath.h"
#include "TFile.h"
#include "TProfile.h"
#include "TF1.h"
#include "StMessMgr.h"

ClassImp(StZdcSmdCorr)

//--------------------------------------------------------

StZdcSmdCorr::StZdcSmdCorr(Int_t energy)
{
  mEnergy = energy;
  for(int i_cent = 0; i_cent < 9; ++i_cent) mResolution[i_cent] = 0.0;
}

StZdcSmdCorr::~StZdcSmdCorr()
{
}

//--------------------------------------------------------

void StZdcSmdCorr::ReadResolution()
{
  string InPutFile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ZDCSMD/Resolution/merged_file/file_%s_Resolution.root",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mBeamEnergy[mEnergy].c_str());
  mFile_Resolution = TFile::Open(InPutFile.c_str());
  p_mResolution = (TProfile*)mFile_Resolution->Get("p_mResolution");
}

void StZdcSmdCorr::CalResolution()
{
  TF1 *f_Res_Full = new TF1("f_Res_Full",Resolution_Full,0,10,0);
  for(Int_t i_cent = 0; i_cent < 9; i_cent++)
  {
    float val_res_full, err_res_full;
    float val_res_raw = p_mResolution->GetBinContent(i_cent+1);
    float err_res_raw = p_mResolution->GetBinError(i_cent+1);
    //    cout << "val_res_raw = " << val_res_raw << ", err_res_raw = " << err_res_raw << endl;
    if(val_res_raw <= 0)
    {
      val_res_full = -999.9;
      err_res_full = 1.0;
    }
    else
    {
      float val_res_sub = TMath::Sqrt(val_res_raw);
      float err_res_sub = err_res_raw/(2*val_res_sub);

      // calculate full event plane resolution
      float chi_sub = f_Res_Full->GetX(val_res_sub);
      float chi_full = chi_sub*TMath::Sqrt(2.0);
      val_res_full = f_Res_Full->Eval(chi_full);
      // error propagation
      float err_chi_sub = err_res_sub/f_Res_Full->Derivative(chi_sub);
      err_res_full = f_Res_Full->Derivative(chi_full)*err_chi_sub*TMath::Sqrt(2.0);
    }
    mResolution[i_cent] = val_res_full;
  }
}

float StZdcSmdCorr::GetResolution(int Cent9)
{
  return mResolution[Cent9];
}

//--------------------------------------------------------

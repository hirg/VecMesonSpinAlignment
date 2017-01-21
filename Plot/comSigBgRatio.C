#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "../Utility/draw.h"
#include "TCanvas.h"

using namespace std;

void comSigBgRatio(int const energy = 2)
{
  string const mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  string inputOld = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/SigBgRatio.root.oldCut",mBeamEnergy[energy].c_str());
  TFile *File_Old = TFile::Open(inputOld.c_str());
  cout << "inputOld = " << inputOld.c_str() << endl;
  cout << inputOld.c_str() << endl;
  TGraphAsymmErrors *g_ratioOld = (TGraphAsymmErrors*)File_Old->Get("g_SigBgRatio");
  TH1F *h_frame = (TH1F*)File_Old->Get("h_frame");

  string inputNew = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/SigBgRatio.root",mBeamEnergy[energy].c_str());
  TFile *File_New = TFile::Open(inputNew.c_str());
  cout << "inputNew = " << inputNew.c_str() << endl;
  TGraphAsymmErrors *g_ratioNew = (TGraphAsymmErrors*)File_New->Get("g_SigBgRatio");


  TCanvas *c_ratio = new TCanvas("c_ratio","c_ratio",10,10,800,800);
  c_ratio->cd();
  c_ratio->cd()->SetLeftMargin(0.15);
  c_ratio->cd()->SetBottomMargin(0.15);
  c_ratio->cd()->SetTicks(1,1);
  c_ratio->cd()->SetGrid(0,0);
  h_frame->GetXaxis()->SetRangeUser(0.0,5.0);
  h_frame->GetYaxis()->SetRangeUser(-0.01,0.1);
  h_frame->Draw("pE");
  plotTopLegend((char*)mBeamEnergy[energy].c_str(),0.9,0.09,0.03,1,0.0,42,0);
  plotTopLegend((char*)"20-60%",1.7,0.09,0.03,1,0.0,42,0);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ratioOld,24,2,1.1);
  Draw_TGAE_Point_new_Symbol(0.8,0.08,0.0,0.0,0.0,0.0,24,2,1.3);
  plotTopLegend((char*)"Old Cut: TPC + ToF(when available)",0.9,0.0785,0.03,1,0.0,42,0);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ratioNew,20,kGray+2,1.1);
  Draw_TGAE_Point_new_Symbol(0.8,0.07,0.0,0.0,0.0,0.0,20,kGray+2,1.3);
  plotTopLegend((char*)"New Cut: TPC + ToF(always)",0.9,0.0685,0.03,1,0.0,42,0);

  string outputFigure = Form("../figures/c_ratio_%s.eps",mBeamEnergy[energy].c_str());
  c_ratio->SaveAs(outputFigure.c_str());
}

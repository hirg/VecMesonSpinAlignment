#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "../Utility/draw.h"

using namespace std;

void comRho00PtRaw(int energy = 2)
{
  string const mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  string inputOld = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/RawRhoPt.root.oldCut",mBeamEnergy[energy].c_str());
  cout << "inputOld = " << inputOld.c_str() << endl;
  TFile *File_Old = TFile::Open(inputOld.c_str());
  TGraphAsymmErrors *g_rhoOld = (TGraphAsymmErrors*)File_Old->Get("Rho00_Phi_Centrality_0_EtaGap_0_Count");

  string inputNew = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/RawRhoPt.root",mBeamEnergy[energy].c_str());
  cout << "inputNew = " << inputNew.c_str() << endl;
  TFile *File_New = TFile::Open(inputNew.c_str());
  TGraphAsymmErrors *g_rhoNew = (TGraphAsymmErrors*)File_New->Get("Rho00_Phi_Centrality_0_EtaGap_0_Count");

  TH1F *h_frame = new TH1F("h_frame","h_frame",100,0,10);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,-10.0);
  }
  TCanvas *c_rho00 = new TCanvas("c_rho00","c_rho00",10,10,800,800);
  c_rho00->cd();
  c_rho00->cd()->SetLeftMargin(0.15);
  c_rho00->cd()->SetBottomMargin(0.15);
  c_rho00->cd()->SetTicks(1,1);
  c_rho00->cd()->SetGrid(0,0);
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,5.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.03);
  h_frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_frame->GetXaxis()->SetTitleSize(0.05);
  h_frame->GetXaxis()->SetTitleOffset(1.2);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.3,0.4);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00}^{raw}");
  h_frame->GetYaxis()->SetTitleSize(0.05);
  h_frame->GetYaxis()->SetLabelSize(0.03);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);
  plotTopLegend((char*)mBeamEnergy[energy].c_str(),0.9,0.39,0.03,1,0.0,42,0);
  plotTopLegend((char*)"20-60%",1.5,0.39,0.03,1,0.0,42,0);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoOld,24,2,1.1);
  Draw_TGAE_Point_new_Symbol(0.8,0.38,0.0,0.0,0.0,0.0,24,2,1.3);
  plotTopLegend((char*)"Old Cut: TPC + ToF(when available)",0.9,0.3786,0.03,1,0.0,42,0);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoNew,20,kGray+2,1.1);
  Draw_TGAE_Point_new_Symbol(0.8,0.37,0.0,0.0,0.0,0.0,20,kGray+2,1.3);
  plotTopLegend((char*)"New Cut: TPC + ToF(always)",0.9,0.3686,0.03,1,0.0,42,0);

  string outputFigure = Form("../figures/c_rho00_%s.eps",mBeamEnergy[energy].c_str());
  c_rho00->SaveAs(outputFigure.c_str());
}

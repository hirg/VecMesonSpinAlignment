#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "../../../Utility/StSpinAlignmentCons.h"
#include "../../../Utility/draw.h"

using namespace std;

float pos_x[8] = {0.10,0.55,0.10,0.55,0.10,0.55,0.10,0.55};
float pos_y[8] = {0.95,0.95,0.90,0.90,0.85,0.85,0.80,0.80};

void plotMcPhiEff(int energy = 6, int cent = 9)
{
  gStyle->SetOptDate(0);
  string InPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/Efficiency/Eff_%s_SingleKaon.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(InPutFile.c_str());
  TH1D *h_mEff[vmsa::pt_rebin];
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt)
  {
    string HistName = Form("h_mEffCos_Cent_%d_Pt_%d",cent,i_pt);
    h_mEff[i_pt] = (TH1D*)File_InPut->Get(HistName.c_str());
  }

  TCanvas *c_eff = new TCanvas("c_eff","c_eff",10,10,800,800);
  c_eff->SetLeftMargin(0.15);
  c_eff->SetBottomMargin(0.15);
  c_eff->SetGrid(0,0);
  c_eff->SetTicks(1,1);
  TH1D *h_frame = new TH1D("h_frame","h_frame",100,0.0,1.0);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10);
    h_frame->SetBinError(i_bin+1,1);
  }
  string legEnergy = Form("AuAu %s 20%%-60%%",vmsa::mBeamEnergy[energy].c_str());
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetTitle("cos(#theta^{*})");
  h_frame->GetXaxis()->CenterTitle();
  h_frame->GetXaxis()->SetLabelSize(0.03);
  h_frame->GetXaxis()->SetNdivisions(505);

  h_frame->GetYaxis()->SetTitle("Efficiency");
  h_frame->GetYaxis()->SetTitleSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->GetYaxis()->SetLabelSize(0.03);
  h_frame->GetYaxis()->SetNdivisions(505);
  h_frame->GetYaxis()->SetRangeUser(0.0,0.95);
  h_frame->Draw("pE");
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt)
  {
    h_mEff[i_pt]->SetMarkerStyle(vmsa::Style[i_pt]);
    h_mEff[i_pt]->SetMarkerColor(vmsa::Color[i_pt]);
    h_mEff[i_pt]->SetMarkerSize(1.4);
    h_mEff[i_pt]->Draw("pEX0 same");
    TF1 *f_poly = new TF1("f_poly","pol1",0.0,1.0);
    h_mEff[i_pt]->Fit(f_poly,"N");
    f_poly->SetLineColor(vmsa::Color[i_pt]);
    f_poly->SetLineStyle(2);
    f_poly->SetLineWidth(2);
    f_poly->Draw("l same");
    string pt_range = Form("p_{T} = %1.1f-%1.1f GeV/c",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
    Draw_TGAE_Point_new_Symbol(pos_x[i_pt],pos_y[i_pt]-0.05,0.0,0.0,0.0,0.0,vmsa::Style[i_pt],vmsa::Color[i_pt],1.0);
    plotTopLegend((char*)pt_range.c_str(),pos_x[i_pt]+0.03,pos_y[i_pt]-0.055,0.03,1,0.0,42,0,1);
  }
  plotTopLegend((char*)legEnergy.c_str(),0.2,0.95,0.06,1,0.0,42,0,1);

  c_eff->SaveAs("../../../figures/effPt.eps");
}

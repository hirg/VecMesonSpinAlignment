#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "draw.h"

static const TString Energy[2] = {"200GeV","39GeV"};
static const TString PID[2] = {"Phi","KStar"};

static const Int_t pt_total = 8; // pt bin
static const Float_t pt_low[pt_total] = {0.4,0.8,1.2,1.8,2.4,3.0,4.2,5.4};
static const Float_t pt_up[pt_total]  = {0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2};

static const Int_t Color[pt_total] = {kGray+2,kBlack,kRed,kCyan,kMagenta,kAzure,kViolet,kBlue};
static const Int_t Style[pt_total] = {20,21,22,23,24,25,26,28};

void PlotEffEvent(Int_t mEnergy = 0, Int_t mPID = 0)
{
  TString InPutFile = Form("/Users/xusun/STAR/SpinAlignment/Data/AuAu%s/%s/Eff_%s_StMcEvent.root",Energy[mEnergy].Data(),PID[mPID].Data(),Energy[mEnergy].Data());
  cout << InPutFile.Data() << endl;
  TFile *File_InPut = TFile::Open(InPutFile.Data());

  TH1F*h_eff[pt_total];
  for(Int_t i_pt = 0; i_pt < pt_total; i_pt++)
  {
    TString HistName = Form("h_mEff_CosTheta_%d",i_pt);
    h_eff[i_pt] = (TH1F*)File_InPut->Get(HistName.Data());
  }

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->cd()->SetLeftMargin(0.15);
  c_play->cd()->SetBottomMargin(0.15);
  c_play->cd()->SetTicks(1,1);
  c_play->cd()->SetGrid(0,0);

  TH1F *h_play = new TH1F("h_play","h_play",7,0.0,1.0);
  for(Int_t i_bin = 0; i_bin < 7; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetRangeUser(0.0,1.0);
  h_play->GetXaxis()->SetNdivisions(505,'N');
  h_play->GetXaxis()->SetLabelSize(0.03);
  h_play->GetXaxis()->SetTitle("cos(#theta*)");
  h_play->GetXaxis()->SetTitleSize(0.05);
  h_play->GetXaxis()->SetTitleOffset(1.2);
  h_play->GetXaxis()->CenterTitle();

  h_play->GetYaxis()->SetRangeUser(0.1,0.95);
  h_play->GetYaxis()->SetNdivisions(505,'N');
  h_play->GetYaxis()->SetTitle("efficiency");
  h_play->GetYaxis()->SetTitleSize(0.05);
  h_play->GetYaxis()->SetLabelSize(0.03);
  h_play->GetYaxis()->CenterTitle();
  h_play->DrawCopy("pE");

  Float_t pos_x[pt_total] = {0.10,0.55,0.10,0.55,0.10,0.55,0.10,0.55};
  Float_t pos_y[pt_total] = {0.95,0.95,0.90,0.90,0.85,0.85,0.80,0.80};
  for(Int_t i_pt = 0; i_pt < pt_total; i_pt++)
  {
    h_eff[i_pt]->SetMarkerStyle(Style[i_pt]);
    h_eff[i_pt]->SetMarkerColor(Color[i_pt]);
    h_eff[i_pt]->SetMarkerSize(1.1);
    h_eff[i_pt]->Draw("pX0 same");
    TString pt_range = Form("p_{T} = %1.1f-%1.1f GeV/c",pt_low[i_pt],pt_up[i_pt]);
    Draw_TGAE_Point_new_Symbol(pos_x[i_pt],pos_y[i_pt]-0.1,0.0,0.0,0.0,0.0,Style[i_pt],Color[i_pt],1.0);
    plotTopLegend(pt_range.Data(),pos_x[i_pt]+0.03,pos_y[i_pt]-0.105,0.03,1,0.0,42,0,1);
    // for(Int_t i_theta = 0; i_theta < 7; i_theta++)
    // {
    //   cout << "Eff = " << p_mEfficiency->GetBinContent(i_pt,i_theta+1) << "+/- " << p_mEfficiency->GetBinError(i_pt+1,i_theta+1) << endl;
    // }
  }

  // c_play->SaveAs("./figures/Eff_CosThetaStar.eps");
}

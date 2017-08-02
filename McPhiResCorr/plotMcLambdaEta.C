#include <string>
#include "TFile.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/draw.h"

void plotMcLambdaEta(int pid = 0)
{
  string InPutHist = Form("/Users/xusun/Data/SpinAlignment/AuAu200GeV/MonteCarlo/McLambdaEta_%d.root",pid);
  TFile *File_InPut = TFile::Open(InPutHist.c_str());

  TProfile *p_cosRP, *p_sinRP;
  p_cosRP = (TProfile*)File_InPut->Get("p_cosRP");
  p_sinRP = (TProfile*)File_InPut->Get("p_sinRP");

  TProfile *p_cosInteDau[20], *p_cosInteLambda[20];
  TProfile *p_sinInteDau[20], *p_sinInteLambda[20];
  for(int i_eta = 0; i_eta < 20; ++i_eta)
  {
    string ProName;
    ProName = Form("p_cosInteDau_%d",i_eta);
    p_cosInteDau[i_eta] = (TProfile*)File_InPut->Get(ProName.c_str());

    ProName = Form("p_cosInteLambda_%d",i_eta);
    p_cosInteLambda[i_eta] = (TProfile*)File_InPut->Get(ProName.c_str());

    ProName = Form("p_sinInteDau_%d",i_eta);
    p_sinInteDau[i_eta] = (TProfile*)File_InPut->Get(ProName.c_str());

    ProName = Form("p_sinInteLambda_%d",i_eta);
    p_sinInteLambda[i_eta] = (TProfile*)File_InPut->Get(ProName.c_str());
  }

  TCanvas *c_PolaPt = new TCanvas("c_PolaPt","c_PolaPt",10,10,800,800);
  c_PolaPt->cd()->SetLeftMargin(0.15);
  c_PolaPt->cd()->SetBottomMargin(0.15);
  c_PolaPt->cd()->SetTicks(1,1);
  c_PolaPt->cd()->SetGrid(0,0);

  p_cosRP->SetTitle("");
  p_cosRP->SetStats(0);
  p_cosRP->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  p_cosRP->GetXaxis()->CenterTitle();
  p_cosRP->GetXaxis()->SetLabelSize(0.04);
  p_cosRP->GetXaxis()->SetNdivisions(505);

  p_cosRP->GetYaxis()->SetTitle("P_{H}");
  p_cosRP->GetYaxis()->SetTitleSize(0.04);
  p_cosRP->GetYaxis()->CenterTitle();
  p_cosRP->GetYaxis()->SetLabelSize(0.04);
  p_cosRP->GetYaxis()->SetNdivisions(505);
  p_cosRP->SetMarkerStyle(20);
  p_cosRP->SetMarkerSize(1.4);
  p_cosRP->SetMarkerColor(kGray+2);
  p_cosRP->Draw("pE");

  p_sinRP->SetMarkerStyle(21);
  p_sinRP->SetMarkerSize(1.4);
  p_sinRP->SetMarkerColor(kAzure-2);
  p_sinRP->Draw("pE same");
  PlotLine(vmsa::ptMin,vmsa::ptMax,0.0,0.0,1,2,2);


  TCanvas *c_PolaEta = new TCanvas("c_PolaEta","c_PolaEta",10,10,800,800);
  c_PolaEta->cd()->SetLeftMargin(0.15);
  c_PolaEta->cd()->SetBottomMargin(0.15);
  c_PolaEta->cd()->SetTicks(1,1);
  c_PolaEta->cd()->SetGrid(0,0);
  TH1F *h_play = new TH1F("h_play","h_play",100,0.0,10.0);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("|#eta| <");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetXaxis()->SetNdivisions(505);

  h_play->GetYaxis()->SetTitle("P_{H}");
  h_play->GetYaxis()->SetTitleSize(0.04);
  h_play->GetYaxis()->CenterTitle();
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetNdivisions(505);
  h_play->GetYaxis()->SetRangeUser(-0.01,0.01);
  h_play->Draw("pE");
  PlotLine(0.0,10.0,0.0,0.0,1,2,2);

  for(int i_eta = 0; i_eta < 20; ++i_eta)
  {
    p_cosInteLambda[i_eta]->SetMarkerStyle(20);
    p_cosInteLambda[i_eta]->SetMarkerSize(1.4);
    p_cosInteLambda[i_eta]->SetMarkerColor(kGray+2);
    p_cosInteLambda[i_eta]->Draw("pE same");

    p_cosInteDau[i_eta]->SetMarkerStyle(24);
    p_cosInteDau[i_eta]->SetMarkerSize(1.4);
    p_cosInteDau[i_eta]->SetMarkerColor(kGray+2);
    p_cosInteDau[i_eta]->Draw("pE same");

    p_sinInteLambda[i_eta]->SetMarkerStyle(21);
    p_sinInteLambda[i_eta]->SetMarkerSize(1.4);
    p_sinInteLambda[i_eta]->SetMarkerColor(kAzure-2);
    // p_sinInteLambda[i_eta]->Draw("pE same");

    p_sinInteDau[i_eta]->SetMarkerStyle(25);
    p_sinInteDau[i_eta]->SetMarkerSize(1.4);
    p_sinInteDau[i_eta]->SetMarkerColor(kAzure-2);
    // p_sinInteDau[i_eta]->Draw("pE same");
  }
}

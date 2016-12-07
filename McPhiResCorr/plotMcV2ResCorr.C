#include <string>
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "../Utility/functions.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"

using namespace std;

std::pair<double, double> const momentumRange(0.2,5.0);
int const pT_low = 1;
int const pT_high = 14;
int const MarkerColorQA = 2;
int const MarkerStyleQA = 24;
int const MarkerColorRP = kGray+2;
int const MarkerStyleRP = 22;
int const MarkerColorEP = kAzure;
int const MarkerStyleEP = 23;

void plotMcV2ResCorr(int energy = 6)
{
  string InPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Phi_v2_1040.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(InPutFile.c_str());
  TGraphAsymmErrors *g_v2 = (TGraphAsymmErrors*)File_InPut->Get("g_v2");
  TF1 *f_v2 = new TF1("f_v2",v2_pT_FitFunc,momentumRange.first,momentumRange.second,5);
  f_v2->FixParameter(0,2);
  f_v2->SetParameter(1,0.1);
  f_v2->SetParameter(2,0.1);
  f_v2->SetParameter(3,0.1);
  f_v2->SetParameter(4,0.1);
  f_v2->SetLineColor(kGray+2);
  f_v2->SetLineWidth(2);
  f_v2->SetLineStyle(2);
  g_v2->Fit(f_v2,"N");

  string InPutHist = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/McV2.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Hist = TFile::Open(InPutHist.c_str());

  TF1 *f_gaus = new TF1("f_gaus","gaus",-TMath::TwoPi(),TMath::TwoPi());
  TH1F *h_Psi2 = (TH1F*)File_Hist->Get("h_Psi2");
  TCanvas *c_Psi2 = new TCanvas("c_Psi2","c_Psi2",100,10,800,800);
  c_Psi2->cd()->SetLeftMargin(0.15);
  c_Psi2->cd()->SetBottomMargin(0.15);
  c_Psi2->cd()->SetTicks(1,1);
  c_Psi2->cd()->SetGrid(0,0);
  h_Psi2->SetTitle("");
  h_Psi2->SetStats(0);
  h_Psi2->GetXaxis()->SetTitle("#Psi_{2}-#Psi_{RP}");
  h_Psi2->GetXaxis()->CenterTitle();
  h_Psi2->GetXaxis()->SetRangeUser(-TMath::TwoPi(),TMath::TwoPi());
  h_Psi2->GetYaxis()->SetTitle("Counts");
  // h_Psi2->GetXaxis()->SetNdivisions(505);

  h_Psi2->GetYaxis()->CenterTitle();
  h_Psi2->GetYaxis()->SetRangeUser(0,1.2*h_Psi2->GetMaximum());
  h_Psi2->GetYaxis()->SetNdivisions(505);
  h_Psi2->SetMarkerStyle(24);
  h_Psi2->SetMarkerColor(kGray+2);
  h_Psi2->SetMarkerSize(1.0);
  h_Psi2->Draw("pE");
  h_Psi2->Fit(f_gaus,"N");
  float sig = f_gaus->GetParameter(2);
  f_gaus->SetRange(-3.0*sig,3.0*sig);
  h_Psi2->Fit(f_gaus,"NR");
  // f_gaus->SetRange(-TMath::PiOver2(),TMath::PiOver2());

  f_gaus->SetLineColor(2);
  f_gaus->SetLineWidth(4);
  f_gaus->SetLineStyle(2);
  f_gaus->Draw("l same");

  float sigma = f_gaus->GetParameter(2);
  float resolution = cos(2.0*sigma);
  string leg_Psi = Form("#sigma = %2.2f, resolution = cos(2*#sigma) = %2.2f",sigma,resolution);
  plotTopLegend(leg_Psi.c_str(),0.2,0.8,0.04,1,0.0,42,1);

  TCanvas *c_v2fitRP = new TCanvas("c_v2fitRP","c_v2fitRP",10,10,1200,600);
  c_v2fitRP->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_v2fitRP->cd(i_pad+1)->SetLeftMargin(0.15);
    c_v2fitRP->cd(i_pad+1)->SetBottomMargin(0.15);
    c_v2fitRP->cd(i_pad+1)->SetTicks(1,1);
    c_v2fitRP->cd(i_pad+1)->SetGrid(0,0);
  }
  TH2F *h_phiRP = (TH2F*)File_Hist->Get("h_phiRP");
  TH1F *h_phiRPproj[20];
  TGraphAsymmErrors *g_v2RP = new TGraphAsymmErrors();
  float delta_pt = (momentumRange.second-momentumRange.first)/20.0;
  for(int i_pt = 0; i_pt < 20; ++i_pt)
  {
    float pt = momentumRange.first+(i_pt+0.5)*delta_pt;
    string phiRPproj = Form("h_phiRPproj_%d",i_pt);
    h_phiRPproj[i_pt] = (TH1F*)h_phiRP->ProjectionY(phiRPproj.c_str(),i_pt+1,i_pt+1);
    TF1 *f_flow = new TF1("f_flow",flow,-1.0*TMath::Pi(),1.0*TMath::Pi(),2);
    f_flow->SetParameter(0,0.1);
    f_flow->SetParameter(1,100);
    h_phiRPproj[i_pt]->Fit(f_flow,"N");
    float v2 = f_flow->GetParameter(0);
    float err_v2 = f_flow->GetParError(0);
    g_v2RP->SetPoint(i_pt,pt+0.02,v2);
    g_v2RP->SetPointError(i_pt,0.0,0.0,err_v2,err_v2);
    if(i_pt == pT_low) c_v2fitRP->cd(1);
    if(i_pt == pT_high) c_v2fitRP->cd(2);
    if(i_pt == pT_low || i_pt == pT_high)
    {
      h_phiRPproj[i_pt]->SetTitle("");
      h_phiRPproj[i_pt]->SetStats(0);
      h_phiRPproj[i_pt]->GetXaxis()->SetTitle("#phi-#Psi_{RP}");
      h_phiRPproj[i_pt]->GetXaxis()->CenterTitle();
      h_phiRPproj[i_pt]->GetXaxis()->SetLabelSize(0.04);
      h_phiRPproj[i_pt]->GetXaxis()->SetNdivisions(505);

      h_phiRPproj[i_pt]->GetYaxis()->SetTitle("Counts");
      h_phiRPproj[i_pt]->GetYaxis()->SetTitleSize(0.04);
      h_phiRPproj[i_pt]->GetYaxis()->CenterTitle();
      h_phiRPproj[i_pt]->GetYaxis()->SetLabelSize(0.04);
      h_phiRPproj[i_pt]->GetYaxis()->SetNdivisions(505);
      h_phiRPproj[i_pt]->GetYaxis()->SetRangeUser(0.7*h_phiRPproj[i_pt]->GetMinimum(),1.2*h_phiRPproj[i_pt]->GetMaximum());

      h_phiRPproj[i_pt]->SetMarkerStyle(MarkerStyleRP);
      h_phiRPproj[i_pt]->SetMarkerSize(1.1);
      h_phiRPproj[i_pt]->SetMarkerColor(MarkerColorRP);
      h_phiRPproj[i_pt]->SetLineColor(MarkerColorRP);
      h_phiRPproj[i_pt]->DrawCopy("pE");
      f_flow->SetLineStyle(2);
      f_flow->SetLineWidth(2);
      f_flow->SetLineColor(2);
      f_flow->Draw("l same");
      string legPt = Form("p_{T} = %2.2f GeV/c",pt);
      string legV2 = Form("v_{2} = %2.3f #pm %0.3f",v2,err_v2);
      plotTopLegend(legPt.c_str(),0.4,0.35,0.04,1,0.0,42,1);
      plotTopLegend(legV2.c_str(),0.4,0.25,0.04,1,0.0,42,1);
    } 
  }

  TCanvas *c_v2fitEP = new TCanvas("c_v2fitEP","c_v2fitEP",10,10,1200,600);
  c_v2fitEP->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_v2fitEP->cd(i_pad+1)->SetLeftMargin(0.15);
    c_v2fitEP->cd(i_pad+1)->SetBottomMargin(0.15);
    c_v2fitEP->cd(i_pad+1)->SetTicks(1,1);
    c_v2fitEP->cd(i_pad+1)->SetGrid(0,0);
  }
  TH2F *h_phiEP = (TH2F*)File_Hist->Get("h_phiEP");
  TH1F *h_phiEPproj[20];
  TGraphAsymmErrors *g_v2EP = new TGraphAsymmErrors();
  float delta_pt = (momentumRange.second-momentumRange.first)/20.0;
  for(int i_pt = 0; i_pt < 20; ++i_pt)
  {
    float pt = momentumRange.first+(i_pt+0.5)*delta_pt;
    string phiEPproj = Form("h_phiEPproj_%d",i_pt);
    h_phiEPproj[i_pt] = (TH1F*)h_phiEP->ProjectionY(phiEPproj.c_str(),i_pt+1,i_pt+1);
    TF1 *f_flow = new TF1("f_flow",flow,-1.0*TMath::Pi(),1.0*TMath::Pi(),2);
    f_flow->SetParameter(0,0.1);
    f_flow->SetParameter(1,100);
    h_phiEPproj[i_pt]->Fit(f_flow,"N");
    float v2 = f_flow->GetParameter(0);
    float err_v2 = f_flow->GetParError(0);
    g_v2EP->SetPoint(i_pt,pt+0.05,v2/resolution);
    g_v2EP->SetPointError(i_pt,0.0,0.0,err_v2/resolution,err_v2/resolution);
    if(i_pt == pT_low) c_v2fitEP->cd(1);
    if(i_pt == pT_high) c_v2fitEP->cd(2);
    if(i_pt == pT_low || i_pt == pT_high)
    {
      h_phiEPproj[i_pt]->SetTitle("");
      h_phiEPproj[i_pt]->SetStats(0);
      h_phiEPproj[i_pt]->GetXaxis()->SetTitle("#phi-#Psi_{EP}");
      h_phiEPproj[i_pt]->GetXaxis()->CenterTitle();
      h_phiEPproj[i_pt]->GetXaxis()->SetLabelSize(0.04);
      h_phiEPproj[i_pt]->GetXaxis()->SetNdivisions(505);

      h_phiEPproj[i_pt]->GetYaxis()->SetTitle("Counts");
      h_phiEPproj[i_pt]->GetYaxis()->SetTitleSize(0.04);
      h_phiEPproj[i_pt]->GetYaxis()->CenterTitle();
      h_phiEPproj[i_pt]->GetYaxis()->SetLabelSize(0.04);
      h_phiEPproj[i_pt]->GetYaxis()->SetNdivisions(505);
      h_phiEPproj[i_pt]->GetYaxis()->SetRangeUser(0.7*h_phiEPproj[i_pt]->GetMinimum(),1.2*h_phiEPproj[i_pt]->GetMaximum());

      h_phiEPproj[i_pt]->SetMarkerStyle(MarkerStyleEP);
      h_phiEPproj[i_pt]->SetMarkerSize(1.1);
      h_phiEPproj[i_pt]->SetMarkerColor(MarkerColorEP);
      h_phiEPproj[i_pt]->SetLineColor(MarkerColorEP);
      h_phiEPproj[i_pt]->DrawCopy("pE");
      f_flow->SetLineStyle(2);
      f_flow->SetLineWidth(2);
      f_flow->SetLineColor(2);
      f_flow->Draw("l same");
      string legPt = Form("p_{T} = %2.2f GeV/c",pt);
      string legV2 = Form("v_{2} = %2.3f #pm %0.3f",v2,err_v2);
      plotTopLegend(legPt.c_str(),0.4,0.35,0.04,1,0.0,42,1);
      plotTopLegend(legV2.c_str(),0.4,0.25,0.04,1,0.0,42,1);
    } 
  }

  TCanvas *c_v2 = new TCanvas("c_v2","c_v2",10,10,800,800);
  c_v2->cd()->SetLeftMargin(0.15);
  c_v2->cd()->SetBottomMargin(0.15);
  c_v2->cd()->SetTicks(1,1);
  c_v2->cd()->SetGrid(0,0);
  TH1F *h_play = new TH1F("h_play","h_play",100,0.0,10.0);
  for(int i_bin = 1; i_bin < 101; ++i_bin)
  {
    h_play->SetBinContent(i_bin,-10.0);
    h_play->SetBinError(i_bin,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetXaxis()->SetRangeUser(momentumRange.first,momentumRange.second);
  h_play->GetXaxis()->SetNdivisions(505);

  h_play->GetYaxis()->SetTitle("v_{2}");
  h_play->GetYaxis()->CenterTitle();
  h_play->GetYaxis()->SetRangeUser(0.0,0.2);
  h_play->GetYaxis()->SetNdivisions(505);
  h_play->Draw("pE");
  f_v2->Draw("l same");

  g_v2RP->SetMarkerStyle(MarkerStyleRP);
  g_v2RP->SetMarkerColor(MarkerColorRP);
  g_v2RP->SetLineColor(MarkerColorRP);
  g_v2RP->SetMarkerSize(1.4);
  g_v2RP->Draw("pE same");

  g_v2EP->SetMarkerStyle(MarkerStyleEP);
  g_v2EP->SetMarkerColor(MarkerColorEP);
  g_v2EP->SetLineColor(MarkerColorEP);
  g_v2EP->SetMarkerSize(1.4);
  g_v2EP->Draw("pE same");

  g_v2->SetMarkerStyle(30);
  g_v2->SetMarkerColor(2);
  g_v2->SetLineColor(2);
  g_v2->SetMarkerSize(2.4);
  g_v2->Draw("pE same");

  plotTopLegend("AuAu 39 GeV 10%-40%",0.4,0.35,0.04,1,0.0,42,1);
  TLegend *legv2 = new TLegend(0.2,0.7,0.55,0.85);
  legv2->SetBorderSize(0.0);
  legv2->SetFillColor(10);
  legv2->AddEntry(g_v2,"#phi STAR PRC 93 014907","p");
  legv2->AddEntry(f_v2,"fit","l");
  legv2->AddEntry(g_v2RP,"v_{2}^{RP}","p");
  legv2->AddEntry(g_v2EP,"v_{2}^{EP}","p");
  legv2->Draw("same");


  TProfile *p_cosRP = (TProfile*)File_Hist->Get("p_cosRP");
  TProfile *p_sinRP = (TProfile*)File_Hist->Get("p_sinRP");
  TProfile *p_cosEP = (TProfile*)File_Hist->Get("p_cosEP");
  TProfile *p_sinEP = (TProfile*)File_Hist->Get("p_sinEP");

  TCanvas *c_cos = new TCanvas("c_cos","c_cos",10,10,800,800);
  c_cos->cd()->SetLeftMargin(0.15);
  c_cos->cd()->SetBottomMargin(0.15);
  c_cos->cd()->SetTicks(1,1);
  c_cos->cd()->SetGrid(0,0);
  p_cosRP->SetTitle("");
  p_cosRP->SetStats(0);
  p_cosRP->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  p_cosRP->GetXaxis()->CenterTitle();
  p_cosRP->GetXaxis()->SetRangeUser(momentumRange.first,momentumRange.second);
  p_cosRP->GetXaxis()->SetNdivisions(505);

  p_cosRP->GetYaxis()->SetTitle("<cos(2#phi)cos(2#Delta#Psi)>");
  p_cosRP->GetYaxis()->CenterTitle();
  p_cosRP->GetYaxis()->SetRangeUser(0.0,0.2);
  p_cosRP->GetYaxis()->SetNdivisions(505);
  p_cosRP->SetMarkerStyle(MarkerStyleRP);
  p_cosRP->SetMarkerColor(MarkerColorRP);
  p_cosRP->SetLineColor(MarkerColorRP);
  p_cosRP->SetMarkerSize(1.4);
  p_cosRP->Draw("pE");

  p_cosEP->SetMarkerStyle(MarkerStyleEP);
  p_cosEP->SetMarkerColor(MarkerColorEP);
  p_cosEP->SetLineColor(MarkerColorEP);
  p_cosEP->SetMarkerSize(1.4);
  p_cosEP->Draw("pE same");

  TCanvas *c_sin = new TCanvas("c_sin","c_sin",10,10,800,800);
  c_sin->cd()->SetLeftMargin(0.15);
  c_sin->cd()->SetBottomMargin(0.15);
  c_sin->cd()->SetTicks(1,1);
  c_sin->cd()->SetGrid(0,0);
  p_sinRP->SetTitle("");
  p_sinRP->SetStats(0);
  p_sinRP->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  p_sinRP->GetXaxis()->CenterTitle();
  p_sinRP->GetXaxis()->SetRangeUser(momentumRange.first,momentumRange.second);
  p_sinRP->GetXaxis()->SetNdivisions(505);

  p_sinRP->GetYaxis()->SetTitle("<sin(2#phi)sin(2#Delta#Psi)>");
  p_sinRP->GetYaxis()->CenterTitle();
  p_sinRP->GetYaxis()->SetRangeUser(-0.1,0.1);
  p_sinRP->GetYaxis()->SetNdivisions(505);
  p_sinRP->SetMarkerStyle(MarkerStyleRP);
  p_sinRP->SetMarkerColor(MarkerColorRP);
  p_sinRP->SetLineColor(MarkerColorRP);
  p_sinRP->SetMarkerSize(1.4);
  p_sinRP->Draw("pE");

  p_sinEP->SetMarkerStyle(MarkerStyleEP);
  p_sinEP->SetMarkerColor(MarkerColorEP);
  p_sinEP->SetLineColor(MarkerColorEP);
  p_sinEP->SetMarkerSize(1.4);
  p_sinEP->Draw("pE same");
  PlotLine(momentumRange.first,momentumRange.second,0.0,0.0,1,2,2);

}

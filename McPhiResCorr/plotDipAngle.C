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

void plotDipAngle(int energy = 6)
{
  string InPutFile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/Phi_v2_1040.root",vmsa::mBeamEnergy[energy].c_str());
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

  string InPutHist = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/McPhiDipAngle.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Hist = TFile::Open(InPutHist.c_str());

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
  TH1F *h_phiRPproj[vmsa::BinPt];
  TGraphAsymmErrors *g_v2RP = new TGraphAsymmErrors();
  float delta_pt = (momentumRange.second-momentumRange.first)/vmsa::BinPt;
  for(int i_pt = 0; i_pt < vmsa::BinPt; ++i_pt)
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

  TH1F *h_PsiEP = (TH1F*)File_Hist->Get("h_PsiEP");
  TCanvas *c_Psi2 = new TCanvas("c_Psi2","c_Psi2",100,10,800,800);
  c_Psi2->cd()->SetLeftMargin(0.15);
  c_Psi2->cd()->SetBottomMargin(0.15);
  c_Psi2->cd()->SetTicks(1,1);
  c_Psi2->cd()->SetGrid(0,0);
  h_PsiEP->SetTitle("");
  h_PsiEP->SetStats(0);
  h_PsiEP->GetXaxis()->SetLabelSize(0.05);
  h_PsiEP->GetXaxis()->SetTitle("#Psi_{2}-#Psi_{RP}");
  h_PsiEP->GetXaxis()->CenterTitle();
  h_PsiEP->GetXaxis()->SetTitleSize(0.05);
  h_PsiEP->GetXaxis()->SetRangeUser(-TMath::PiOver2(),TMath::PiOver2());
  h_PsiEP->GetXaxis()->SetNdivisions(505);

  h_PsiEP->GetYaxis()->SetLabelSize(0.05);
  h_PsiEP->GetYaxis()->SetTitle("Counts");
  h_PsiEP->GetYaxis()->CenterTitle();
  h_PsiEP->GetYaxis()->SetTitleSize(0.05);
  h_PsiEP->GetYaxis()->SetNdivisions(505);
  h_PsiEP->GetYaxis()->SetRangeUser(0,1.2*h_PsiEP->GetMaximum());
  h_PsiEP->SetMarkerStyle(24);
  h_PsiEP->SetMarkerColor(2);
  h_PsiEP->SetMarkerSize(1.0);
  h_PsiEP->Draw("pE");
  TF1 *f_EP = new TF1("f_EP",EventPlaneDist,-TMath::PiOver2(),TMath::PiOver2(),2);
  f_EP->SetParameter(0,1.5);
  f_EP->SetParameter(1,1000);
  h_PsiEP->Fit(f_EP,"N");
  f_EP->SetLineColor(2);
  f_EP->SetLineWidth(4);
  f_EP->SetLineStyle(2);
  f_EP->Draw("l same");
  float chi = f_EP->GetParameter(0);
  TF1 *f_res = new TF1("f_res",EventPlaneResolution,0,10,0);
  float resolution = f_res->Eval(chi);
  string leg_PsiEP = Form("Event Plane: #chi = %2.2f, resolution = %2.2f",chi,resolution);
  plotTopLegend((char*)leg_PsiEP.c_str(),0.2,0.85,0.04,2,0.0,42,1);

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
  TH1F *h_phiEPproj[vmsa::BinPt];
  TGraphAsymmErrors *g_v2EP = new TGraphAsymmErrors();
  for(int i_pt = 0; i_pt < vmsa::BinPt; ++i_pt)
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

  TCanvas *c_v2DA = new TCanvas("c_v2DA","c_v2DA",10,10,1200,600);
  c_v2DA->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_v2DA->cd(i_pad+1)->SetLeftMargin(0.15);
    c_v2DA->cd(i_pad+1)->SetBottomMargin(0.15);
    c_v2DA->cd(i_pad+1)->SetTicks(1,1);
    c_v2DA->cd(i_pad+1)->SetGrid(0,0);
  }
  TH2F *h_phiDA = (TH2F*)File_Hist->Get("h_phiDA");
  TH1F *h_phiDAproj[vmsa::BinPt];
  TGraphAsymmErrors *g_v2DA = new TGraphAsymmErrors();
  for(int i_pt = 0; i_pt < vmsa::BinPt; ++i_pt)
  {
    float pt = momentumRange.first+(i_pt+0.5)*delta_pt;
    string phiDAproj = Form("h_phiDAproj_%d",i_pt);
    h_phiDAproj[i_pt] = (TH1F*)h_phiDA->ProjectionY(phiDAproj.c_str(),i_pt+1,i_pt+1);
    TF1 *f_flow = new TF1("f_flow",flow,-1.0*TMath::Pi(),1.0*TMath::Pi(),2);
    f_flow->SetParameter(0,0.1);
    f_flow->SetParameter(1,100);
    h_phiDAproj[i_pt]->Fit(f_flow,"N");
    float v2 = f_flow->GetParameter(0);
    float err_v2 = f_flow->GetParError(0);
    g_v2DA->SetPoint(i_pt,pt+0.05,v2/resolution);
    g_v2DA->SetPointError(i_pt,0.0,0.0,err_v2/resolution,err_v2/resolution);
    if(i_pt == pT_low) c_v2DA->cd(1);
    if(i_pt == pT_high) c_v2DA->cd(2);
    if(i_pt == pT_low || i_pt == pT_high)
    {
      h_phiDAproj[i_pt]->SetTitle("");
      h_phiDAproj[i_pt]->SetStats(0);
      h_phiDAproj[i_pt]->GetXaxis()->SetTitle("#phi-#Psi_{EP}");
      h_phiDAproj[i_pt]->GetXaxis()->CenterTitle();
      h_phiDAproj[i_pt]->GetXaxis()->SetLabelSize(0.04);
      h_phiDAproj[i_pt]->GetXaxis()->SetNdivisions(505);

      h_phiDAproj[i_pt]->GetYaxis()->SetTitle("Counts");
      h_phiDAproj[i_pt]->GetYaxis()->SetTitleSize(0.04);
      h_phiDAproj[i_pt]->GetYaxis()->CenterTitle();
      h_phiDAproj[i_pt]->GetYaxis()->SetLabelSize(0.04);
      h_phiDAproj[i_pt]->GetYaxis()->SetNdivisions(505);
      h_phiDAproj[i_pt]->GetYaxis()->SetRangeUser(0.7*h_phiDAproj[i_pt]->GetMinimum(),1.2*h_phiDAproj[i_pt]->GetMaximum());

      h_phiDAproj[i_pt]->SetMarkerStyle(MarkerStyleEP);
      h_phiDAproj[i_pt]->SetMarkerSize(1.1);
      h_phiDAproj[i_pt]->SetMarkerColor(MarkerColorEP);
      h_phiDAproj[i_pt]->SetLineColor(MarkerColorEP);
      h_phiDAproj[i_pt]->DrawCopy("pE");
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

  g_v2DA->SetMarkerStyle(MarkerStyleQA);
  g_v2DA->SetMarkerColor(MarkerColorQA);
  g_v2DA->SetLineColor(MarkerColorQA);
  g_v2DA->SetMarkerSize(1.4);
  g_v2DA->Draw("pE same");

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
}

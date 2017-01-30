#include <string>
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "../Utility/functions.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"

using namespace std;

void calPhiResCorr(int energy = 6, int pid = 0, int centrality = 0)
{
  gStyle->SetOptDate(0);
  string InPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/McPhiResCorr.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(InPutFile.c_str());
  TH2F *h_cosRP[60], *h_cosGaus[60], *h_cosEP[60];
  TGraphAsymmErrors *g_rhoRP[60], *g_rhoGaus[60], *g_rhoEP[60];
  TF1 *f_polRP[60], *f_polGaus[60], *f_polEP[60]; 
  TGraphAsymmErrors *g_ResGaus = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_ResEP = new TGraphAsymmErrors();
  float delta_pt = (vmsa::ptMax-vmsa::ptMin)/(float)vmsa::BinPt;
  for(int i_rho = 0; i_rho < 60; ++i_rho)
  {
    string HistCosRP = Form("h_cosRP_%d",i_rho);
    h_cosRP[i_rho] = (TH2F*)File_InPut->Get(HistCosRP.c_str());
    g_rhoRP[i_rho] = new TGraphAsymmErrors();
    string FuncCosRP = Form("f_cosRP_%d",i_rho);
    f_polRP[i_rho] = new TF1(FuncCosRP.c_str(),"pol0",vmsa::ptMin,vmsa::ptMax);
    f_polRP[i_rho]->SetParameter(0,0.01*(i_rho));

    string HistCosGaus = Form("h_cosGaus_%d",i_rho);
    h_cosGaus[i_rho] = (TH2F*)File_InPut->Get(HistCosGaus.c_str());
    g_rhoGaus[i_rho] = new TGraphAsymmErrors();
    string FuncCosGaus = Form("f_cosGaus_%d",i_rho);
    f_polGaus[i_rho] = new TF1(FuncCosGaus.c_str(),"pol0",vmsa::ptMin,vmsa::ptMax);
    f_polGaus[i_rho]->SetParameter(0,0.01*(i_rho));

    string HistCosEP = Form("h_cosEP_%d",i_rho);
    h_cosEP[i_rho] = (TH2F*)File_InPut->Get(HistCosEP.c_str());
    g_rhoEP[i_rho] = new TGraphAsymmErrors();
    string FuncCosEP = Form("f_cosEP_%d",i_rho);
    f_polEP[i_rho] = new TF1(FuncCosEP.c_str(),"pol0",vmsa::ptMin,vmsa::ptMax);
    f_polEP[i_rho]->SetParameter(0,0.01*(i_rho));

    for(int i_pt = 0; i_pt < vmsa::BinPt; ++i_pt)
    {
      float pt = vmsa::ptMin+(i_pt+0.5)*delta_pt;
      string HistCosProjRP = Form("h_rhoRP_%d_proj_%d",i_rho,i_pt);
      TH1F *h_cosRPproj = (TH1F*)h_cosRP[i_rho]->ProjectionY(HistCosProjRP.c_str(),i_pt+1,i_pt+1);
      TF1 *f_rhoRP = new TF1("f_rhoRP",SpinDensity,-1.0,1.0,2);
      f_rhoRP->SetParameter(0,0.01*(i_rho));
      f_rhoRP->SetParameter(1,100);
      h_cosRPproj->Fit(f_rhoRP,"NQ");
      float rhoRP = f_rhoRP->GetParameter(0);
      float err_rhoRP = f_rhoRP->GetParError(0);
      g_rhoRP[i_rho]->SetPoint(i_pt,pt,rhoRP);
      g_rhoRP[i_rho]->SetPointError(i_pt,0.0,0.0,err_rhoRP,err_rhoRP);

      string HistCosProjGaus = Form("h_rhoGaus_%d_proj_%d",i_rho,i_pt);
      TH1F *h_cosGausproj = (TH1F*)h_cosGaus[i_rho]->ProjectionY(HistCosProjGaus.c_str(),i_pt+1,i_pt+1);
      TF1 *f_rhoGaus = new TF1("f_rhoGaus",SpinDensity,-1.0,1.0,2);
      f_rhoGaus->SetParameter(0,0.01*(i_rho));
      f_rhoGaus->SetParameter(1,100);
      h_cosGausproj->Fit(f_rhoGaus,"NQ");
      float rhoGaus = f_rhoGaus->GetParameter(0);
      float err_rhoGaus = f_rhoGaus->GetParError(0);
      g_rhoGaus[i_rho]->SetPoint(i_pt,pt,rhoGaus);
      g_rhoGaus[i_rho]->SetPointError(i_pt,0.0,0.0,err_rhoGaus,err_rhoGaus);

      string HistCosProjEP = Form("h_rhoEP_%d_proj_%d",i_rho,i_pt);
      TH1F *h_cosEPproj = (TH1F*)h_cosEP[i_rho]->ProjectionY(HistCosProjEP.c_str(),i_pt+1,i_pt+1);
      TF1 *f_rhoEP = new TF1("f_rhoEP",SpinDensity,-1.0,1.0,2);
      f_rhoEP->SetParameter(0,0.01*(i_rho));
      f_rhoEP->SetParameter(1,100);
      h_cosEPproj->Fit(f_rhoEP,"NQ");
      float rhoEP = f_rhoEP->GetParameter(0);
      float err_rhoEP = f_rhoEP->GetParError(0);
      g_rhoEP[i_rho]->SetPoint(i_pt,pt,rhoEP);
      g_rhoEP[i_rho]->SetPointError(i_pt,0.0,0.0,err_rhoEP,err_rhoEP);
    }
    g_rhoRP[i_rho]->Fit(f_polRP[i_rho],"N");
    g_rhoGaus[i_rho]->Fit(f_polGaus[i_rho],"N");
    g_rhoEP[i_rho]->Fit(f_polEP[i_rho],"N");
    g_ResGaus->SetPoint(i_rho,f_polRP[i_rho]->GetParameter(0),f_polGaus[i_rho]->GetParameter(0));
    g_ResGaus->SetPointError(i_rho,f_polRP[i_rho]->GetParError(0),f_polRP[i_rho]->GetParError(0),f_polGaus[i_rho]->GetParError(0),f_polGaus[i_rho]->GetParError(0));
    g_ResEP->SetPoint(i_rho,f_polRP[i_rho]->GetParameter(0),f_polEP[i_rho]->GetParameter(0));
    g_ResEP->SetPointError(i_rho,f_polRP[i_rho]->GetParError(0),f_polRP[i_rho]->GetParError(0),f_polEP[i_rho]->GetParError(0),f_polEP[i_rho]->GetParError(0));
  }

  TCanvas *c_res = new TCanvas("c_res","c_res",10,10,800,800);
  c_res->SetLeftMargin(0.15);
  c_res->SetBottomMargin(0.15);
  c_res->SetTicks(1,1);
  c_res->SetGrid(0,0);
  TH1F *h_res = new TH1F("h_res","h_res",100,0.0,0.6);
  h_res->SetTitle("");
  h_res->SetStats(0);
  h_res->GetXaxis()->SetLabelSize(0.05);
  h_res->GetXaxis()->SetTitle("#rho_{00}^{phy}");
  h_res->GetXaxis()->CenterTitle();
  h_res->GetXaxis()->SetTitleSize(0.05);
  h_res->GetXaxis()->SetRangeUser(0.0,0.6);
  h_res->GetXaxis()->SetNdivisions(505);

  h_res->GetYaxis()->SetLabelSize(0.05);
  h_res->GetYaxis()->SetTitle("#rho_{00}^{obs}");
  h_res->GetYaxis()->SetTitleSize(0.05);
  h_res->GetYaxis()->CenterTitle();
  h_res->GetYaxis()->SetRangeUser(0.0,0.6);
  h_res->GetYaxis()->SetNdivisions(505);
  h_res->Draw("pE");
  PlotLine(1.0/3.0,1.0/3.0,0.0,0.6,1,2,2);
  PlotLine(0.0,0.6,1.0/3.0,1.0/3.0,1,2,2);

  g_ResGaus->SetMarkerStyle(20);
  g_ResGaus->SetMarkerColor(1);
  g_ResGaus->SetMarkerSize(1.4);
  g_ResGaus->Draw("PE same");
  TF1 *f_ResGaus = new TF1("f_ResGaus",PolyRes,0.0,0.6,2);
  f_ResGaus->SetParameter(0,0.01);
  f_ResGaus->SetParameter(1,0.8);
  g_ResGaus->Fit(f_ResGaus,"N");
  f_ResGaus->SetLineColor(kGray+2);
  f_ResGaus->SetLineWidth(2);
  f_ResGaus->SetLineStyle(2);
  f_ResGaus->Draw("l same");
  // string leg_p0Gaus   = Form("p_{0}^{Gaus} = %2.4f #pm %2.4f",f_ResGaus->GetParameter(0),f_ResGaus->GetParError(0));
  // string leg_p1Gaus   = Form("p_{1}^{Gaus} = %2.4f #pm %2.4f",f_ResGaus->GetParameter(1),f_ResGaus->GetParError(1));
  string leg_p1Gaus   = Form("p_{1}^{Gaus} = %2.2f",f_ResGaus->GetParameter(1));
  string leg_chi2Gaus = Form("#chi^{2}/NDF = %2.2f/%d",f_ResGaus->GetChisquare(),f_ResGaus->GetNDF());
  // plotTopLegend((char*)leg_p0Gaus.c_str(),0.2,0.80,0.03,1,0.0,42,1);
  plotTopLegend((char*)leg_p1Gaus.c_str(),0.2,0.74,0.04,1,0.0,42,1);
  plotTopLegend((char*)leg_chi2Gaus.c_str(),0.2,0.68,0.04,1,0.0,42,1);

  g_ResEP->SetMarkerStyle(24);
  g_ResEP->SetMarkerColor(kRed);
  g_ResEP->SetMarkerSize(1.4);
  g_ResEP->Draw("PE same");
  TF1 *f_ResEP = new TF1("f_ResEP",PolyRes,0.0,0.6,2);
  f_ResEP->SetParameter(0,0.01);
  f_ResEP->SetParameter(1,0.8);
  g_ResEP->Fit(f_ResEP,"N");
  f_ResEP->SetLineColor(2);
  f_ResEP->SetLineWidth(2);
  f_ResEP->SetLineStyle(2);
  f_ResEP->Draw("l same");
  // string leg_p0EP   = Form("p_{0}^{EP} = %2.4f #pm %2.4f",f_ResEP->GetParameter(0),f_ResEP->GetParError(0));
  // string leg_p1EP   = Form("p_{1}^{EP} = %2.4f #pm %2.4f",f_ResEP->GetParameter(1),f_ResEP->GetParError(1));
  string leg_p1EP   = Form("p_{1}^{EP} = %2.2f",f_ResEP->GetParameter(1));
  string leg_chi2EP = Form("#chi^{2}/NDF = %2.2f/%d",f_ResEP->GetChisquare(),f_ResEP->GetNDF());
  // plotTopLegend((char*)leg_p0EP.c_str(),0.58,0.40,0.03,kRed,0.0,42,1);
  plotTopLegend((char*)leg_p1EP.c_str(),0.58,0.34,0.04,kRed,0.0,42,1);
  plotTopLegend((char*)leg_chi2EP.c_str(),0.58,0.28,0.04,kRed,0.0,42,1);

  string OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/MonteCarlo/McResCorr/Mc%sResCorrFactor.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_res->Write();
  string KEY_ResGaus = Form("ResGaus_Centrality_%d_%s",centrality,vmsa::mPID[pid].c_str());
  g_ResGaus->SetName(KEY_ResGaus.c_str());
  g_ResGaus->Write();
  string KEY_ResEP = Form("ResEP_Centrality_%d_%s",centrality,vmsa::mPID[pid].c_str());
  g_ResEP->SetName(KEY_ResEP.c_str());
  g_ResEP->Write();
  File_OutPut->Close();

  string FigureName = Form("../figures/resCorrFactor%s.png",vmsa::mBeamEnergy[energy].c_str());
  c_res->SaveAs(FigureName.c_str());
}


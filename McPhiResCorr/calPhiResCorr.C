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

using namespace std;

double SpinDensity(double *x_val, double *par)
{
  double x = x_val[0];
  double rho00 = par[0];
  double Norm = par[1];

  double dNdCosThetaStar = Norm*((1.0-rho00)+(3.0*rho00-1)*x*x);

  return dNdCosThetaStar;
}

double Poly(double *x_val, double *par)
{
  Double_t x = x_val[0];
  Double_t y = par[0] + par[1]*(x-1.0/3.0) + 1.0/3.0;

  return y;
}

void PlotLine(double x1_val, double x2_val, double y1_val, double y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
{
  TLine* Zero_line = new TLine();
  Zero_line -> SetX1(x1_val);
  Zero_line -> SetX2(x2_val);
  Zero_line -> SetY1(y1_val);
  Zero_line -> SetY2(y2_val);
  Zero_line -> SetLineWidth(LineWidth);
  Zero_line -> SetLineStyle(LineStyle);
  Zero_line -> SetLineColor(Line_Col);
  Zero_line -> Draw();
  //delete Zero_line;
}

TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1)
{
  // coordinates in NDC!
  // plots the string label in position x and y in NDC coordinates
  // size is the text size
  // color is the text color

  //    if(x<0||y<0)
  //    {   // defaults
  //      x=gPad->GetLeftMargin()*1.15;
  //      y=(1-gPad->GetTopMargin())*1.04;
  //    }
  TLatex* text=new TLatex(x,y,label);
  text->SetTextFont(font);
  text->SetTextSize(size);
  if(NDC == 1) text->SetNDC();
  text->SetTextColor(color);
  text->SetTextAngle(angle);
  text->Draw();
  return text;
}

string const  mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
pair<double, double> const momentumRange(0.2,5.0);

void calPhiResCorr(int energy = 6)
{
  string InPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/McPhiV2.root",mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(InPutFile.c_str());
  TH2F *h_cosRP[60], *h_cosEP[60];
  TGraphAsymmErrors *g_rhoRP[60], *g_rhoEP[60];
  TF1 *f_polRP[60], *f_polEP[60]; 
  TGraphAsymmErrors *g_Res = new TGraphAsymmErrors();
  float delta_pt = (momentumRange.second-momentumRange.first)/20.0;
  for(int i_rho = 0; i_rho < 60; ++i_rho)
  {
    string HistCosRP = Form("h_cosRP_%d",i_rho+1);
    h_cosRP[i_rho] = (TH2F*)File_InPut->Get(HistCosRP.c_str());
    g_rhoRP[i_rho] = new TGraphAsymmErrors();
    string FuncCosRP = Form("f_cosRP_%d",i_rho+1);
    f_polRP[i_rho] = new TF1(FuncCosRP.c_str(),"pol0",momentumRange.first,momentumRange.second);
    f_polRP[i_rho]->SetParameter(0,0.01*(i_rho+1));

    string HistCosEP = Form("h_cosEP_%d",i_rho+1);
    h_cosEP[i_rho] = (TH2F*)File_InPut->Get(HistCosEP.c_str());
    g_rhoEP[i_rho] = new TGraphAsymmErrors();
    string FuncCosEP = Form("f_cosEP_%d",i_rho+1);
    f_polEP[i_rho] = new TF1(FuncCosEP.c_str(),"pol0",momentumRange.first,momentumRange.second);
    f_polEP[i_rho]->SetParameter(0,0.01*(i_rho+1));

    for(int i_pt = 0; i_pt < 20; ++i_pt)
    {
      float pt = momentumRange.first+(i_pt+0.5)*delta_pt;
      string HistCosProjRP = Form("h_rhoRP_%d_proj_%d",i_rho+1,i_pt);
      TH1F *h_cosRPproj = (TH1F*)h_cosRP[i_rho]->ProjectionY(HistCosProjRP.c_str(),i_pt+1,i_pt+1);
      TF1 *f_rhoRP = new TF1("f_rhoRP",SpinDensity,-1.0,1.0,2);
      f_rhoRP->SetParameter(0,0.01*(i_rho+1));
      f_rhoRP->SetParameter(1,100);
      h_cosRPproj->Fit(f_rhoRP,"NQ");
      float rhoRP = f_rhoRP->GetParameter(0);
      float err_rhoRP = f_rhoRP->GetParError(0);
      g_rhoRP[i_rho]->SetPoint(i_pt,pt,rhoRP);
      g_rhoRP[i_rho]->SetPointError(i_pt,0.0,0.0,err_rhoRP,err_rhoRP);

      string HistCosProjEP = Form("h_rhoEP_%d_proj_%d",i_rho+1,i_pt);
      TH1F *h_cosEPproj = (TH1F*)h_cosEP[i_rho]->ProjectionY(HistCosProjEP.c_str(),i_pt+1,i_pt+1);
      TF1 *f_rhoEP = new TF1("f_rhoEP",SpinDensity,-1.0,1.0,2);
      f_rhoEP->SetParameter(0,0.01*(i_rho+1));
      f_rhoEP->SetParameter(1,100);
      h_cosEPproj->Fit(f_rhoEP,"NQ");
      float rhoEP = f_rhoEP->GetParameter(0);
      float err_rhoEP = f_rhoEP->GetParError(0);
      g_rhoEP[i_rho]->SetPoint(i_pt,pt,rhoEP);
      g_rhoEP[i_rho]->SetPointError(i_pt,0.0,0.0,err_rhoEP,err_rhoEP);
    }
    g_rhoRP[i_rho]->Fit(f_polRP[i_rho],"N");
    g_rhoEP[i_rho]->Fit(f_polEP[i_rho],"N");
    g_Res->SetPoint(i_rho,f_polRP[i_rho]->GetParameter(0),f_polEP[i_rho]->GetParameter(0));
    g_Res->SetPointError(i_rho,f_polRP[i_rho]->GetParError(0),f_polRP[i_rho]->GetParError(0),f_polEP[i_rho]->GetParError(0),f_polEP[i_rho]->GetParError(0));
  }

  TCanvas *c_res = new TCanvas("c_res","c_res",10,10,800,800);
  c_res->SetLeftMargin(0.15);
  c_res->SetBottomMargin(0.15);
  c_res->SetTicks(1,1);
  c_res->SetGrid(0,0);
  TH1F *h_res = new TH1F("h_res","h_res",100,0.0,0.6);
  h_res->SetTitle("");
  h_res->SetStats(0);
  h_res->GetXaxis()->SetTitle("#rho_{00}^{phy}");
  h_res->GetXaxis()->CenterTitle();
  h_res->GetXaxis()->SetRangeUser(0.0,0.6);
  h_res->GetXaxis()->SetNdivisions(505);

  h_res->GetYaxis()->SetTitle("#rho_{00}^{obs}");
  h_res->GetYaxis()->CenterTitle();
  h_res->GetYaxis()->SetRangeUser(0.0,0.6);
  h_res->GetYaxis()->SetNdivisions(505);
  h_res->Draw("pE");
  PlotLine(1.0/3.0,1.0/3.0,0.0,0.6,1,2,2);
  PlotLine(0.0,0.6,1.0/3.0,1.0/3.0,1,2,2);
  g_Res->SetMarkerStyle(24);
  g_Res->SetMarkerColor(kGray+2);
  g_Res->SetMarkerSize(1.2);
  g_Res->Draw("PE same");
  
  TF1 *f_Res = new TF1("f_Res",Poly,0.0,0.6,2);
  f_Res->SetParameter(0,0.01);
  f_Res->SetParameter(1,0.8);
  g_Res->Fit(f_Res,"N");
  f_Res->SetLineColor(2);
  f_Res->SetLineWidth(2);
  f_Res->SetLineStyle(2);
  f_Res->Draw("l same");
  string leg_p0 = Form("p_{0} = %2.4f #pm %2.4f",f_Res->GetParameter(0),f_Res->GetParError(0));
  string leg_p1 = Form("p_{1} = %2.4f #pm %2.4f",f_Res->GetParameter(1),f_Res->GetParError(1));
  string leg_chi2 = Form("#chi^{2}/NDF = %2.2f/%d",f_Res->GetChisquare(),f_Res->GetNDF());
  plotTopLegend(leg_p0.c_str(),0.2,0.80,0.04,1,0.0,42,1);
  plotTopLegend(leg_p1.c_str(),0.2,0.74,0.04,1,0.0,42,1);
  plotTopLegend(leg_chi2.c_str(),0.2,0.68,0.04,1,0.0,42,1);

  // c_res->SaveAs("c_res.eps");
}


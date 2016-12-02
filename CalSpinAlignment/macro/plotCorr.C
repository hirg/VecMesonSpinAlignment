#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "draw.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TMath.h"

Double_t Poly(Double_t *x_val, Double_t *par) 
{
  Double_t x = x_val[0];
  Double_t y = par[0] + par[1]*(x-1.0/3.0) + 1.0/3.0;

  return y;
}

static const TString Energy[2] = {"200GeV","39GeV"};
static const TString PID[2] = {"Phi","KStar"};

static const Int_t Centrality_total = 10;
static const Int_t Centrality_start = 0;
static const Int_t Centrality_stop  = 10;
static const TString mCentrality[10] = {"2060","0005","0510","1020","2030","3040","4050","5060","6070","7080"};
static const TString mCentLabel[10] = {"20-60%","0-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%"};

static const TString mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
static const TString mMode_AMPT[2] = {"Default","StringMelting"};
static const TString mScreenMass_AMPT[3] = {"1mb","3mb","6mb"};

void PlotCorr(Int_t mEnergy = 0, Int_t mPID = 0)
{
  const Int_t mColor[Centrality_total] = {kGray+2,kBlack,kRed,kOrange,kMagenta,kAzure,kCyan,kGreen,kViolet,kBlue};
  const Int_t mStyle[Centrality_total] = {20,21,22,23,24,25,26,27,28,29};

  TString inputfile = Form("./Data/AuAu%s/%s/ResCorr_AMPT.root",Energy[mEnergy].Data(),PID[mPID].Data());
  TFile *File_InPut = TFile::Open(inputfile.Data());

  TGraphAsymmErrors *g_rho[Centrality_total];
  TGraphAsymmErrors *g_rho_Eff[Centrality_total];
  Double_t val_p0[10], err_p0[10], val_p1[10], err_p1[10];
  for(Int_t i_cent = Centrality_start; i_cent < Centrality_stop; i_cent++)
  {
    TString GraName = Form("g_rho00_Cent_%d",i_cent);
    g_rho[i_cent] = (TGraphAsymmErrors*)File_InPut->Get(GraName.Data());
    TF1 *f_poly = new TF1("f_poly",Poly,0.0,0.6,2);
    f_poly->SetParameter(0,0.2);
    f_poly->SetParameter(1,1.0);
    g_rho[i_cent]->Fit(f_poly,"NQ");
    val_p0[i_cent] = f_poly->GetParameter(0);
    err_p0[i_cent] = f_poly->GetParError(0);
    val_p1[i_cent] = f_poly->GetParameter(1);
    err_p1[i_cent] = f_poly->GetParError(1);
    TString GraName_Eff = Form("g_rho00_Cent_%d_Eff",i_cent);
    g_rho_Eff[i_cent] = (TGraphAsymmErrors*)File_InPut->Get(GraName_Eff.Data());
  }

  Float_t x_start = 0.0;
  Float_t x_stop  = 0.6;
  Float_t y_start = 0.0;
  Float_t y_stop  = 0.6;
  TCanvas *c_rho00 = new TCanvas("c_rho00","c_rho00",10,10,800,800);
  c_rho00->cd();
  c_rho00->cd()->SetLeftMargin(0.15);
  c_rho00->cd()->SetBottomMargin(0.15);
  c_rho00->cd()->SetTicks(1,1);
  c_rho00->cd()->SetGrid(0,0);
  TH1F *h_play = new TH1F("h_play","h_play",1000,-0.001,0.999);
  for(Int_t i_bin = 1; i_bin < 1001; i_bin++)
  {
    h_play->SetBinContent(i_bin,-10.0);
    h_play->SetBinError(i_bin,0.1);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetNdivisions(505,'X');
  h_play->GetYaxis()->SetNdivisions(505,'Y');
  h_play->GetXaxis()->SetTitle("#rho_{00}^{phy}");
  h_play->GetYaxis()->SetTitle("#rho_{00}^{obs}");
  h_play->GetXaxis()->SetTitleSize(0.04);
  h_play->GetYaxis()->SetTitleSize(0.04);
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->CenterTitle();
  h_play->GetXaxis()->SetRangeUser(x_start,x_stop);
  h_play->GetYaxis()->SetRangeUser(y_start,y_stop);
  h_play->SetLineColor(1);
  h_play->SetMarkerStyle(24);
  h_play->SetMarkerColor(1);
  h_play->SetMarkerSize(0.8);
  h_play->DrawCopy("pE");
  PlotLine(x_start,x_stop,1.0/3.0,1.0/3.0,1,2,2);
  PlotLine(1.0/3.0,1.0/3.0,y_start,y_stop,1,2,2);

  /*
  for(Int_t i_cent = 0; i_cent < Centrality_total; i_cent++)
  {
    Draw_TGAE_new_Symbol(g_rho[i_cent],mStyle[i_cent],mColor[i_cent],1.2);
    Float_t pos_x = x_start+0.1;
    Float_t pos_y = y_stop-0.04-i_cent*0.02;
    Draw_TGAE_Point_new_Symbol(pos_x,pos_y,0.0,0.0,0.0,0.0,mStyle[i_cent],mColor[i_cent],1.0);
    plotTopLegend((char*)mCentLabel[i_cent].Data(),pos_x+0.01,pos_y-0.006,0.025,1,0.0,42,0);

    TF1 *f_poly = new TF1("f_poly",Poly,0.0,0.6,2);
    f_poly->SetParameter(0,val_p0[i_cent]);
    f_poly->SetParameter(1,val_p1[i_cent]);
    f_poly->SetLineColor(mColor[i_cent]);
    f_poly->Draw("same");
  }
  */
  
  Draw_TGAE_new_Symbol(g_rho[2],20,kGray,1.2);
  Float_t pos_x = x_start+0.1;
  Float_t pos_y = y_stop-0.06;
  Draw_TGAE_Point_new_Symbol(pos_x,pos_y,0.0,0.0,0.0,0.0,20,kGray,1.0);
  plotTopLegend((char*)mCentLabel[2].Data(),pos_x+0.01,pos_y-0.006,0.025,1,0.0,42,0);

  Draw_TGAE_new_Symbol(g_rho[4],22,8,1.2);
  Float_t pos_x = x_start+0.1;
  Float_t pos_y = y_stop-0.08;
  Draw_TGAE_Point_new_Symbol(pos_x,pos_y,0.0,0.0,0.0,0.0,22,8,1.0);
  plotTopLegend((char*)mCentLabel[4].Data(),pos_x+0.01,pos_y-0.006,0.025,1,0.0,42,0);

  Draw_TGAE_new_Symbol(g_rho[7],24,kAzure,1.2);
  Float_t pos_x = x_start+0.1;
  Float_t pos_y = y_stop-0.10;
  Draw_TGAE_Point_new_Symbol(pos_x,pos_y,0.0,0.0,0.0,0.0,24,kAzure,1.0);
  plotTopLegend((char*)mCentLabel[7].Data(),pos_x+0.01,pos_y-0.006,0.025,1,0.0,42,0);

  Draw_TGAE_new_Symbol(g_rho[8],34,38,1.2);
  Float_t pos_x = x_start+0.1;
  Float_t pos_y = y_stop-0.12;
  Draw_TGAE_Point_new_Symbol(pos_x,pos_y,0.0,0.0,0.0,0.0,34,38,1.0);
  plotTopLegend((char*)mCentLabel[9].Data(),pos_x+0.01,pos_y-0.006,0.025,1,0.0,42,0);

  Draw_TGAE_new_Symbol(g_rho[0],29,2,1.4);
  Float_t pos_x = x_start+0.1;
  Float_t pos_y = y_stop-0.04;
  Draw_TGAE_Point_new_Symbol(pos_x,pos_y,0.0,0.0,0.0,0.0,29,2,1.0);
  plotTopLegend((char*)mCentLabel[0].Data(),pos_x+0.01,pos_y-0.006,0.025,1,0.0,42,0);
  TF1 *f_poly = new TF1("f_poly",Poly,0.0,0.6,2);
  f_poly->SetParameter(0,val_p0[0]);
  f_poly->SetParameter(1,val_p1[0]);
  f_poly->SetLineColor(2);
  f_poly->Draw("same");

  plotTopLegend("Au+Au 200 GeV",x_stop-0.2,y_start+0.15,0.025,1,0.0,42,0);
  TString AMPT_version = Form("AMPT %s %s",mMode_AMPT[1].Data(),mScreenMass_AMPT[1].Data());
  plotTopLegend((char*)AMPT_version.Data(),x_stop-0.23,y_start+0.13,0.025,1,0.0,42,0);
  plotTopLegend("w/o TPC Efficiency Correction",x_stop-0.26,y_stop-0.066,0.025,1,0.0,42,0);
  c_rho00->SaveAs("./figures/rho00_Centrality.eps");

  /*
  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->cd();
  c_play->cd()->SetLeftMargin(0.15);
  c_play->cd()->SetBottomMargin(0.15);
  c_play->cd()->SetTicks(1,1);
  c_play->cd()->SetGrid(0,0);
  h_play->Draw("pE");
  PlotLine(x_start,x_stop,1.0/3.0,1.0/3.0,1,2,2);
  PlotLine(1.0/3.0,1.0/3.0,y_start,y_stop,1,2,2);
  Draw_TGAE_new_Symbol(g_rho[0],mStyle[0],mColor[0],1.2);
  Draw_TGAE_Point_new_Symbol(x_start+0.05,y_stop-0.04,0.0,0.0,0.0,0.0,mStyle[0],mColor[0],1.0);
  plotTopLegend("w/o TPC Efficiency Correction",x_start+0.06,y_stop-0.046,0.025,1,0.0,42,0);

  Draw_TGAE_new_Symbol(g_rho_Eff[0],24,2,1.2);
  Draw_TGAE_Point_new_Symbol(x_start+0.05,y_stop-0.06,0.0,0.0,0.0,0.0,24,2,1.0);
  plotTopLegend("with TPC Efficiency Correction",x_start+0.06,y_stop-0.066,0.025,1,0.0,42,0);

  plotTopLegend("Au+Au 200 GeV",x_stop-0.2,y_start+0.15,0.025,1,0.0,42,0);
  TString AMPT_version = Form("AMPT %s %s",mMode_AMPT[1].Data(),mScreenMass_AMPT[1].Data());
  plotTopLegend((char*)AMPT_version.Data(),x_stop-0.23,y_start+0.13,0.025,1,0.0,42,0);
  plotTopLegend((char*)mCentLabel[0].Data(),x_stop-0.16,y_stop-0.066,0.025,1,0.0,42,0);
  c_play->SaveAs("./figures/rho00_Com.eps");
  */

  TString input_AMPT = Form("./Data/AuAu%s/%s/Resolution_200GeV.root",Energy[mEnergy].Data(),PID[mPID].Data());
  TFile *File_InPutAMPT = TFile::Open(input_AMPT.Data());
  TProfile *p_mRes_before = (TProfile*)File_InPutAMPT->Get("p_mRes_before");
  Float_t res2[10];
  res2[0] = 0.0;
  for(Int_t i_cent = 1; i_cent < 10; i_cent++)
  {
    res2[i_cent] = TMath::Sqrt(p_mRes_before->GetBinContent(10-i_cent));
  }

  TGraphAsymmErrors *g_p0_Cent = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_p1_Cent = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_p0_Res = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_p1_Res = new TGraphAsymmErrors();
  Float_t Centrality_start[10] = {20,0, 5,10,20,30,40,50,60,70};
  Float_t Centrality_stop[10]  = {60,5,10,20,30,40,50,60,70,80};
  for(Int_t i_cent = 0; i_cent < 10; i_cent++)
  {
    g_p0_Cent->SetPoint(i_cent,0.5*(Centrality_start[i_cent]+Centrality_stop[i_cent]),val_p0[i_cent]);
    g_p0_Cent->SetPointError(i_cent,0.0,0.0,err_p0[i_cent],err_p0[i_cent]);
    g_p1_Cent->SetPoint(i_cent,0.5*(Centrality_start[i_cent]+Centrality_stop[i_cent]),val_p1[i_cent]);
    g_p1_Cent->SetPointError(i_cent,0.0,0.0,err_p1[i_cent],err_p1[i_cent]);

    g_p0_Res->SetPoint(i_cent,res2[i_cent],val_p0[i_cent]);
    g_p0_Res->SetPointError(i_cent,0.0,0.0,err_p0[i_cent],err_p0[i_cent]);
    g_p1_Res->SetPoint(i_cent,res2[i_cent],val_p1[i_cent]);
    g_p1_Res->SetPointError(i_cent,0.0,0.0,err_p1[i_cent],err_p1[i_cent]);
  }
  g_p0_Res->RemovePoint(0);
  g_p1_Res->RemovePoint(0);

  TCanvas *c_p0_Cent = new TCanvas("c_p0_Cent","c_p0_Cent",10,10,800,800);
  c_p0_Cent->cd();
  c_p0_Cent->cd()->SetLeftMargin(0.15);
  c_p0_Cent->cd()->SetBottomMargin(0.15);
  c_p0_Cent->cd()->SetTicks(1,1);
  c_p0_Cent->cd()->SetGrid(0,0);
  TH1F *h_poly = new TH1F("h_poly","h_poly",100,0,100);
  for(Int_t i_bin = 1; i_bin < 101; i_bin++)
  {
    h_poly->SetBinContent(i_bin,-10.0);
    h_poly->SetBinError(i_bin,0.1);
  }
  h_poly->SetTitle("");
  h_poly->SetStats(0);
  h_poly->GetXaxis()->SetNdivisions(505,'X');
  h_poly->GetYaxis()->SetNdivisions(505,'Y');
  h_poly->GetXaxis()->SetTitle("Centrality(%)");
  h_poly->GetYaxis()->SetTitle("parameters");
  h_poly->GetXaxis()->SetTitleSize(0.04);
  h_poly->GetYaxis()->SetTitleSize(0.04);
  h_poly->GetYaxis()->SetTitleOffset(1.2);
  h_poly->GetXaxis()->CenterTitle();
  h_poly->GetYaxis()->CenterTitle();
  h_poly->GetXaxis()->SetRangeUser(0,80);
  h_poly->GetYaxis()->SetRangeUser(-0.1,0.8);
  h_poly->Draw("pE");
  PlotLine(0.0,80,0.0,0.0,1,2,2);
  Draw_TGAE_new_Symbol(g_p0_Cent,24,kGray,1.2);
  Draw_TGAE_Point_new_Symbol(65,0.75,0.0,0.0,0.0,0.0,24,kGray,1.4);
  plotTopLegend("p_{0}",67,0.749,0.03,1,0.0,42,0);
  Draw_TGAE_new_Symbol(g_p1_Cent,24,2,1.2);
  Draw_TGAE_Point_new_Symbol(65,0.71,0.0,0.0,0.0,0.0,24,2,1.4);
  plotTopLegend("p_{1}",67,0.699,0.03,1,0.0,42,0);
  c_p0_Cent->SaveAs("./figures/c_p0_Cent.eps");

  TCanvas *c_p0_Res = new TCanvas("c_p0_Res","c_p0_Res",10,10,800,800);
  c_p0_Res->cd();
  c_p0_Res->cd()->SetLeftMargin(0.15);
  c_p0_Res->cd()->SetBottomMargin(0.15);
  c_p0_Res->cd()->SetTicks(1,1);
  c_p0_Res->cd()->SetGrid(0,0);
  TH1F *h_Res = new TH1F("h_Res","h_Res",100,0,1.0);
  for(Int_t i_bin = 1; i_bin < 101; i_bin++)
  {
    h_Res->SetBinContent(i_bin,-10.0);
    h_Res->SetBinError(i_bin,0.1);
  }
  h_Res->SetTitle("");
  h_Res->SetStats(0);
  h_Res->GetXaxis()->SetNdivisions(505,'X');
  h_Res->GetYaxis()->SetNdivisions(505,'Y');
  h_Res->GetXaxis()->SetTitle("2^{nd} Event Plane Resolution");
  h_Res->GetYaxis()->SetTitle("parameters");
  h_Res->GetXaxis()->SetTitleSize(0.04);
  h_Res->GetYaxis()->SetTitleSize(0.04);
  h_Res->GetYaxis()->SetTitleOffset(1.2);
  h_Res->GetXaxis()->CenterTitle();
  h_Res->GetYaxis()->CenterTitle();
  h_Res->GetXaxis()->SetRangeUser(0.2,0.8);
  h_Res->GetYaxis()->SetRangeUser(-0.1,0.8);
  h_Res->Draw("pE");
  PlotLine(0.2,0.8,0.0,0.0,1,2,2);
  Draw_TGAE_new_Symbol(g_p0_Res,24,kGray,1.2);
  Draw_TGAE_Point_new_Symbol(0.35,0.75,0.0,0.0,0.0,0.0,24,kGray,1.4);
  plotTopLegend("p_{0}",0.37,0.749,0.03,1,0.0,42,0);
  Draw_TGAE_new_Symbol(g_p1_Res,24,2,1.2);
  Draw_TGAE_Point_new_Symbol(0.35,0.71,0.0,0.0,0.0,0.0,24,2,1.4);
  plotTopLegend("p_{1}",0.37,0.699,0.03,1,0.0,42,0);
  c_p0_Res->SaveAs("./figures/c_p0_Res.eps");
}

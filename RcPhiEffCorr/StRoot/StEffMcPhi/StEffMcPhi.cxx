#include "StEffMcPhi.h"
#include "StEffHistManger.h"
#include "StEffCut.h"
#include <string>
#include "TFile.h"
#include "TNtuple.h"
#include "TBranch.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "../../../Utility/StSpinAlignmentCons.h"

ClassImp(StEffMcPhi)

int StEffMcPhi::mInput_flag = 1;

StEffMcPhi::StEffMcPhi(int Energy, long StartEvent, long StopEvent, int PID, int year, int cut)
{
  energy = Energy;
  pid = PID;

  string InPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/Efficiency/Eff_%s_SingleKaon_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str());

  SetInPutFile(InPutFile); // set input list

  SetStartEvent(StartEvent); // set start event
  SetStopEvent(StopEvent); // set stop event

  string OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/Efficiency/Eff_%s_SingleKaon.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  SetOutPutFile(OutPutFile); // set output file

  mEffCut = new StEffCut();
  mEffHistManger = new StEffHistManger(energy);
}

StEffMcPhi::~StEffMcPhi()
{
}

//------------------------------------------------------------
void StEffMcPhi::SetInPutFile(const string inputfile)
{
  mInPutFile = inputfile;
  cout << "Input file was set to: " << mInPutFile.c_str() << endl;
}

void StEffMcPhi::SetOutPutFile(const string outputfile)
{
  mOutPutFile = outputfile;
  cout << "Output file was set to: " << mOutPutFile.c_str() << endl;
}

void StEffMcPhi::SetStartEvent(const long StartEvent)
{
  mStartEvent = StartEvent;
  cout << "nStartEvent = " << mStartEvent << endl;
}

void StEffMcPhi::SetStopEvent(const long StopEvent)
{
  mStopEvent = StopEvent;
  cout << "nStopEvent = " << mStopEvent << endl;
}
//------------------------------------------------------------

void StEffMcPhi::Init()
{
  mEffHistManger->InitHist();

  // initialize the TNtuple
  mFile_InPut = TFile::Open(mInPutFile.c_str());
  cout << "OPEN InPut File: " << mInPutFile.c_str() << endl;
  mNtuple = (TNtuple*)mFile_InPut->Get("McPhiMeson");

  // initialize Ntuple
  mNtuple->SetBranchAddress("Centrality",&mCentrality);
  mNtuple->SetBranchAddress("McPt",&mMcPt);
  mNtuple->SetBranchAddress("McP",&mMcP);
  mNtuple->SetBranchAddress("McEta",&mMcEta);
  mNtuple->SetBranchAddress("McY",&mMcY);
  mNtuple->SetBranchAddress("McPhi",&mMcPhi);
  mNtuple->SetBranchAddress("McInvMass",&mMcInvMass);
  mNtuple->SetBranchAddress("McPid",&mMcPid);

  mNtuple->SetBranchAddress("KpMcPt",&mKpMcPt);
  mNtuple->SetBranchAddress("KpMcEta",&mKpMcEta);
  mNtuple->SetBranchAddress("KpMcY",&mKpMcY);
  mNtuple->SetBranchAddress("KpMcPhi",&mKpMcPhi);
  mNtuple->SetBranchAddress("KpMcM",&mKpMcM);
  mNtuple->SetBranchAddress("KpMcPid",&mKpMcPid);

  mNtuple->SetBranchAddress("KpRcPt",&mKpRcPt);
  mNtuple->SetBranchAddress("KpRcEta",&mKpRcEta);
  mNtuple->SetBranchAddress("KpRcY",&mKpRcY);
  mNtuple->SetBranchAddress("KpRcPhi",&mKpRcPhi);
  mNtuple->SetBranchAddress("KpRcM",&mKpRcM);
  mNtuple->SetBranchAddress("KpRcTpc",&mKpRcTpc);
  
  mNtuple->SetBranchAddress("KmMcPt",&mKmMcPt);
  mNtuple->SetBranchAddress("KmMcEta",&mKmMcEta);
  mNtuple->SetBranchAddress("KmMcY",&mKmMcY);
  mNtuple->SetBranchAddress("KmMcPhi",&mKmMcPhi);
  mNtuple->SetBranchAddress("KmMcM",&mKmMcM);
  mNtuple->SetBranchAddress("KmMcPid",&mKmMcPid);

  mNtuple->SetBranchAddress("KmRcPt",&mKmRcPt);
  mNtuple->SetBranchAddress("KmRcEta",&mKmRcEta);
  mNtuple->SetBranchAddress("KmRcY",&mKmRcY);
  mNtuple->SetBranchAddress("KmRcPhi",&mKmRcPhi);
  mNtuple->SetBranchAddress("KmRcM",&mKmRcM);
  mNtuple->SetBranchAddress("KmRcTpc",&mKmRcTpc);

  mNtuple->SetBranchAddress("RcPt",&mRcPt);
  mNtuple->SetBranchAddress("RcP",&mRcP);
  mNtuple->SetBranchAddress("RcEta",&mRcEta);
  mNtuple->SetBranchAddress("RcY",&mRcY);
  mNtuple->SetBranchAddress("RcPhi",&mRcPhi);
  mNtuple->SetBranchAddress("RcInvMass",&mRcInvMass);

  int num_tracks = mNtuple->GetEntriesFast();
  cout << "Number of tracks in McPhiMeson = " << num_tracks<< endl;

  if(mStartEvent > num_tracks) mStartEvent = num_tracks;
  if(mStopEvent  > num_tracks) mStopEvent  = num_tracks;
  cout << "New nStartEvent = " << mStartEvent << ", new nStopEvent = " << mStopEvent << endl;

  mFile_OutPut = new TFile(mOutPutFile.c_str(),"RECREATE");
}

void StEffMcPhi::Make()
{
  long start_event_use = mStartEvent;
  long stop_event_use  = mStopEvent;

  gRandom->SetSeed();
  mNtuple->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

  for(long i_track = start_event_use; i_track < stop_event_use; ++i_track)
  {
    if (!mNtuple->GetEntry(i_track)) // take track information
      break;  // end of data chunk
    if (floor(10.0*i_track/ static_cast<float>(stop_event_use)) > floor(10.0*(i_track-1)/ static_cast<float>(stop_event_use)))
      cout << "=> processing data: " << 100.0*i_track/ static_cast<float>(stop_event_use) << "%" << endl;

    McVecMeson McPhi; // initialize McPhi
    McPhi.Centrality = mCentrality;
    McPhi.McPt       = mMcPt;
    McPhi.McP        = mMcP;
    McPhi.McEta      = mMcEta;
    McPhi.McY        = mMcY;
    McPhi.McPhi      = mMcPhi;
    McPhi.McInvMass  = mMcInvMass;
    McPhi.McPid      = mMcPid;

    McDecayDau McKP; // initialize McKplus
    McKP.McPt  = mKpMcPt;
    McKP.McEta = mKpMcEta;
    McKP.McY   = mKpMcY;
    McKP.McPhi = mKpMcPhi;
    McKP.McM   = mKpMcM;
    McKP.McPid = mKpMcPid;

    RcDecayDau RcKP; // initialize RcKplus
    RcKP.RcPt  = mKpRcPt;
    RcKP.RcEta = mKpRcEta;
    RcKP.RcY   = mKpRcY;
    RcKP.RcPhi = mKpRcPhi;
    RcKP.RcM   = mKpRcM;
    RcKP.RcTpc = mKpRcTpc;

    McDecayDau McKM; // initialize McKminus
    McKM.McPt  = mKmMcPt;
    McKM.McEta = mKmMcEta;
    McKM.McY   = mKmMcY;
    McKM.McPhi = mKmMcPhi;
    McKM.McM   = mKmMcM;
    McKM.McPid = mKmMcPid;

    RcDecayDau RcKM; // initialize RcKminus
    RcKM.RcPt  = mKmRcPt;
    RcKM.RcEta = mKmRcEta;
    RcKM.RcY   = mKmRcY;
    RcKM.RcPhi = mKmRcPhi;
    RcKM.RcM   = mKmRcM;
    RcKM.RcTpc = mKmRcTpc;

    RcVecMeson RcPhi; // initialize RcPhi
    RcPhi.RcPt      = mRcPt;
    RcPhi.RcP       = mRcP;
    RcPhi.RcEta     = mRcEta;
    RcPhi.RcY       = mRcY;
    RcPhi.RcPhi     = mRcPhi;
    RcPhi.RcInvMass = mRcInvMass;

    if( !mEffCut->passTrackCutPhi(McPhi) ) continue; // eta cuts for McPhi 
    if( !mEffCut->passTrackCut(McKP) ) continue; // eta cut for McKplus
    if( !mEffCut->passTrackCut(McKM) ) continue; // eta cut for McKminus

    TLorentzVector lMcPhi;
    lMcPhi.SetPtEtaPhiM(McPhi.McPt,McPhi.McEta,McPhi.McPhi,vmsa::InvMass[pid]);
    TVector3 vMcPhiBeta = -1.0*lMcPhi.BoostVector();

    TLorentzVector lMcKP;
    lMcKP.SetPtEtaPhiM(McKP.McPt,McKP.McEta,McKP.McPhi,vmsa::mMassKaon);
    lMcKP.Boost(vMcPhiBeta);
    TVector3 vMcKP = lMcKP.Vect().Unit(); // direction of K+ momentum in phi-meson rest frame
    float Psi2 = gRandom->Uniform(-0.5*TMath::Pi(),0.5*TMath::Pi()); // random event plane angle
    // float Psi2 = gRandom->Uniform(-TMath::Pi(),TMath::Pi()); // random event plane angle
    TVector3 QVector(TMath::Sin(Psi2),-1.0*TMath::Cos(Psi2),0.0);
    TVector3 nQ = QVector.Unit(); // direction of QVector
    float McCosThetaStar = vMcKP.Dot(nQ);

    mEffHistManger->FillHistMc(McPhi.Centrality,McPhi.McPt,McPhi.McEta,McPhi.McPhi,McCosThetaStar);
    // mEffHistManger->FillHistMc(McPhi.Centrality,McKP.McPt,McKP.McEta,McKP.McPhi,McCosThetaStar);
    // mEffHistManger->FillHistMc(McPhi.Centrality,McKM.McPt,McKM.McEta,McKM.McPhi,McCosThetaStar);

    if( !mEffCut->passTrackCut(RcKP) ) continue; // eta and TPC cuts for RcKplus
    if( !mEffCut->passTrackCut(RcKM) ) continue; // eta and TPC cuts for RcKminus
    if( !mEffCut->passTrackCutPhi(RcPhi) ) continue;  // eta cuts for RcPhi 
    mEffHistManger->FillHistRc(McPhi.Centrality,McPhi.McPt,McPhi.McEta,McPhi.McPhi,McCosThetaStar);
    // mEffHistManger->FillHistRc(McPhi.Centrality,McKP.McPt,McKP.McEta,McKP.McPhi,McCosThetaStar);
    // mEffHistManger->FillHistRc(McPhi.Centrality,McKM.McPt,McKM.McEta,McKM.McPhi,McCosThetaStar);
  }
  cout << "=> processing data: 100%" << endl;
  cout << "work done!" << endl;

  mEffHistManger->CalEfficiency();
  // mEffHistManger->CalEffPtEtaPhi();
  mEffHistManger->CalEffCosThetaStar();
}

void StEffMcPhi::Finish()
{
  mFile_OutPut->cd();
  mEffHistManger->WriteHist();
  mFile_OutPut->Close();
}

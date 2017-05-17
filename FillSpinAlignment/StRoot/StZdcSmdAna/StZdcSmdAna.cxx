#include "StRoot/StZdcSmdAna/StZdcSmdAna.h"
// #include "StRoot/StZdcSmdAna/StZdcSmdCut.h"
#include "StRoot/StZdcSmdAna/StZdcSmdCorr.h"
#include "StRoot/StZdcSmdAna/StZdcSmdHistoManger.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StAlexPhiMesonEvent/StAlexPhiMesonEvent.h"
#include "StRoot/StRunIdEventsDb/StRunIdEventsDb.h"
#include "StThreeVectorF.hh"
#include "StMessMgr.h"
#include "TFile.h"
#include "TChain.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <fstream>
#include "TStopwatch.h"

ClassImp(StZdcSmdAna)

StRefMultCorr* StZdcSmdAna::mRefMultCorr = NULL;
int StZdcSmdAna::mInPut_flag = 1;
char* StZdcSmdAna::VM_EVENT_TREE = NULL;
char* StZdcSmdAna::VM_EVENT_BRANCH = NULL;

//----------------------------------------------------
StZdcSmdAna::StZdcSmdAna(int energy, int X_flag, int List, Long64_t start_event, Long64_t stop_event, int mode)
{
  mEnergy = energy;
  mX_flag = X_flag;
  mList = List;
  mStart_Event = start_event;
  mStop_Event = stop_event;
  mMode = mode;
  if(!mRefMultCorr)
  {
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
  mZdcSmdCorr = new StZdcSmdCorr(mEnergy);
  // mZdcSmdCut = new StZdcSmdCut();
  mZdcSmdHistoManger = new StZdcSmdHistoManger();
  mStopWatch = new TStopwatch();
}

StZdcSmdAna::~StZdcSmdAna()
{
}
//----------------------------------------------------
// set Input/Output
void StZdcSmdAna::setInputDir(const TString inputdir)
{
  mInputdir = inputdir.Copy();
  cout << "Input directory was set to: " << mInputdir.Data() << endl;
}
void StZdcSmdAna::setOutputfile(const TString outputfile)
{
  mOutputfile = outputfile.Copy();
  cout << "Output file was set to: " << mOutputfile.Data() << endl;
}
void StZdcSmdAna::setInPutList(const TString InPutList)
{
  mInPutList = InPutList.Copy();
  TString InFo_InPutList = Form("InPut %s list was set to: %s",vmsa::MixEvent[mX_flag].Data(),mInPutList.Data());
  cout << InFo_InPutList.Data() << endl;
}
void StZdcSmdAna::setStopEvent(const Long64_t StopEvent)
{
    mStopEvent = StopEvent;
    cout << "nStopEvent = " << mStopEvent << endl;
}
void StZdcSmdAna::setStartEvent(const Long64_t StartEvent)
{
    mStartEvent = StartEvent;
    cout << "nStartEvent = " << mStartEvent << endl;
}
//----------------------------------------------------
// initial functions
void StZdcSmdAna::Init()
{
  mStopWatch->Start();
  mZdcSmdCorr->ReadResolution();
  mZdcSmdCorr->CalResolution();
  mZdcSmdHistoManger->Init(mX_flag,mMode);

  TString inputdir = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ZDCSMD/%s/Forest/",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mPID[mMode].c_str());
  setInputDir(inputdir);

  const int list_start = vmsa::mList_Delta*mList + 1; // start list
  const int list_stop  = vmsa::mList_Delta*(mList+1); // stop list

  TString InPutList = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ZDCSMD/%s/List/Split_%s_%s_%d_%d.list",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mPID[mMode].c_str(),vmsa::MixEvent[mX_flag].Data(),vmsa::mBeamEnergy[mEnergy].c_str(),list_start,list_stop);
  setInPutList(InPutList);

  TString outputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ZDCSMD/%s/Yields/Yields_%s_%s_%d.root",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mPID[mMode].c_str(),vmsa::MixEvent[mX_flag].Data(),vmsa::mBeamEnergy[mEnergy].c_str(),mList);
  setOutputfile(outputfile);

  setStartEvent(Long64_t(mStart_Event));
  setStopEvent(Long64_t(mStop_Event));
  //----------------------------------------------------------------------------------------------------

  TString Notification = Form("Initializing parameters and input/output for %s %s",vmsa::mPID[mMode].c_str(),vmsa::MixEvent[mX_flag].Data());
  cout << Notification.Data() << endl;
  mFile_OutPut = new TFile(mOutputfile.Data(),"RECREATE");

  VM_EVENT_TREE       = (char*)vmsa::vm_tree[mMode].Data();
  VM_EVENT_BRANCH     = (char*)vmsa::vm_branch[mMode].Data();

  //----------------------------------------------------------------------------------------------------
  // input
  if (!mInPutList.IsNull())   // if input file is ok
  {
    TString InFo_List = Form("Open %s file list ",vmsa::MixEvent[mX_flag].Data());
    cout << InFo_List.Data() << endl;
    ifstream in(mInPutList);  // input stream
    if(in)
    {
      cout << "input file list is ok" << endl;
      mInPut = new TChain( VM_EVENT_TREE, VM_EVENT_TREE );
      char str[255];       // char array for each file name
      Long64_t entries_save = 0;
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  TString addfile;
	  addfile = str;
	  addfile = mInputdir+addfile;
	  mInPut->AddFile(addfile.Data(),-1, VM_EVENT_TREE );
	  Long64_t file_entries = mInPut->GetEntries();
	  cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
	  entries_save = file_entries;
	}
      }
    }
    else
    {
      TString InFo_Warning = Form("WARNING: %s file input is problemtic",vmsa::MixEvent[mX_flag].Data());
      cout << InFo_Warning.Data() << endl;
      mInPut_flag = 0;
    }
  }

  // Set the input tree
  if (mInPut_flag == 1 && !mInPut->GetBranch( VM_EVENT_BRANCH ))
  {
    cerr << "ERROR: Could not find branch '"
      << VM_EVENT_BRANCH << "'in tree!" << endl;
  }

  if(mMode == 0) mPhiMeson_event = new StAlexPhiMesonEvent();

  if(mInPut_flag == 1)
  {
    if(mMode == 0) mInPut->SetBranchAddress( VM_EVENT_BRANCH, &mPhiMeson_event );

    int num_events = mInPut->GetEntriesFast();
    cout << "Number of events in file(s) = " << num_events << endl;
    if(mStartEvent > num_events) mStartEvent = num_events;
    if(mStopEvent > num_events) mStopEvent   = num_events;

    cout << "New nStartEvent = " << mStartEvent << ", new nStopEvent = " << mStopEvent << endl;
  }
}

void StZdcSmdAna::Make()
{
  if(mMode == 0) MakePhi();
}


// loop phi meson Same Event
void StZdcSmdAna::MakePhi()
{
  Long64_t start_event_use;
  Long64_t stop_event_use;

  start_event_use = mStartEvent;
  stop_event_use  = mStopEvent;
  mInPut->SetBranchAddress( VM_EVENT_BRANCH, &mPhiMeson_event);
  mInPut->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

  // Initialise Event Head
  StThreeVectorF PrimaryVertex(0.0,0.0,0.0);
  int          RunId = 0;
  int          EventId = 0;
  int          RefMult = 0;
  int          Centrality = 0;
  int          N_prim = 0;
  int          N_non_prim = 0;
  int          N_Tof_match = 0;
  float        ZDCx = 0.0; 
  float        BBCx = 0.0; 
  float        VzVpd = 0.0;
  int          NumTrackUsed = 0;
  // ---------------------------------------QVector---------------------------------------------
  TVector2 QEast(0.0,0.0);
  TVector2 QWest(0.0,0.0);
  TVector2 QFull(0.0,0.0);
  // -----------------------------------Number of Tracks----------------------------------------
  int   NumTrackEast = 0;
  int   NumTrackWest = 0;
  int   NumTrackFull = 0;
  int   NumTrackFullEast = 0;
  int   NumTrackFullWest = 0;

  for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
  {
    if (!mInPut->GetEntry( counter )) // take the event -> information is stored in event
      break;  // end of data chunk

    // get Event Header
    PrimaryVertex    = mPhiMeson_event->getPrimaryVertex();
    RunId            = mPhiMeson_event->getRunId();
    EventId          = mPhiMeson_event->getEventId();
    RefMult          = mPhiMeson_event->getRefMult();
    Centrality       = mPhiMeson_event->getCentrality();
    N_prim           = mPhiMeson_event->getN_prim();
    N_non_prim       = mPhiMeson_event->getN_non_prim();
    N_Tof_match      = mPhiMeson_event->getN_Tof_match();
    ZDCx             = mPhiMeson_event->getZDCx(); 
    BBCx             = mPhiMeson_event->getBBCx(); 
    VzVpd            = mPhiMeson_event->getVzVpd();
    NumTrackUsed     = mPhiMeson_event->getNumTracks();
    QEast           = mPhiMeson_event->getQ2East();
    QWest           = mPhiMeson_event->getQ2West();
    QFull           = mPhiMeson_event->getQ2Full();
    NumTrackEast     = mPhiMeson_event->getNumTrackEast();
    NumTrackWest     = mPhiMeson_event->getNumTrackWest();
    NumTrackFull     = mPhiMeson_event->getNumTrackFull();
    NumTrackFullEast = mPhiMeson_event->getNumTrackFullEast();
    NumTrackFullWest = mPhiMeson_event->getNumTrackFullWest();

    // Initialise Track 
    Float_t m2A = -999.9;
    Float_t m2B = -999.9;
    Float_t nsA = -999.9;
    Float_t nsB = -999.9;
    Float_t dcaA = -999.9;
    Float_t dcaB = -999.9;
    TLorentzVector lTrackA(0.0,0.0,0.0,0.0);
    TLorentzVector lTrackB(0.0,0.0,0.0,0.0);
    int flagA = -1;
    int flagB = -1;

    // Centrality
    mRefMultCorr->init(RunId);
    if(mEnergy == 6) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z(),ZDCx); // 200 GeV
    if(mEnergy != 6) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z());
    const int cent9 = Centrality;
    const Double_t reweight = mRefMultCorr->getWeight();

    if (counter != 0  &&  counter % 1000 == 0)
      cout << "." << flush;
    if (counter != 0  &&  counter % 10000 == 0)
    {
      if((stop_event_use-start_event_use) > 0)
      {
	Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
	cout << " " << counter-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (ZdcSmdSpinAlignment) " << flush;
      }
    }

    // get Track Information
    for(UShort_t nTracks = 0; nTracks < NumTrackUsed; nTracks++) // loop over all tracks of the actual event
    {
      mPhiMeson_track = mPhiMeson_event->getTrack(nTracks);
      m2A = mPhiMeson_track->getMass2A();
      m2B = mPhiMeson_track->getMass2B();
      nsA = mPhiMeson_track->getNSigKaonA();
      nsB = mPhiMeson_track->getNSigKaonB();
      dcaA = mPhiMeson_track->getDcaA();
      dcaB = mPhiMeson_track->getDcaB();
      lTrackA = mPhiMeson_track->getTrackA();
      lTrackB = mPhiMeson_track->getTrackB();
      flagA = mPhiMeson_track->getFlagA();
      flagB = mPhiMeson_track->getFlagB();

      Float_t pA = lTrackA.P();
      Float_t pB = lTrackB.P();
      TLorentzVector lTrack = lTrackA + lTrackB; // phi-meson
      Float_t pt_lTrack = lTrack.Perp();

      if( // with TPC+ToF (if possible) at low momentum and TPC+ToF (always) at high momentum 
	  ((fabs(pA) <= 0.65 && m2A < -10) || (m2A > 0 && ((fabs(pA) < 1.5 && m2A > 0.16 && m2A < 0.36) || (fabs(pA) >= 1.5 && m2A > 0.125 && m2A < 0.36)) )) &&
	  ((fabs(pB) <= 0.65 && m2B < -10) || (m2B > 0 && ((fabs(pB) < 1.5 && m2B > 0.16 && m2B < 0.36) || (fabs(pB) >= 1.5 && m2B > 0.125 && m2B < 0.36)) )) &&
	  (pt_lTrack < 0.8 || (pt_lTrack >= 0.8 && ( (m2A > 0.16 && m2A < 0.36) || (m2B > 0.16 && m2B < 0.36)))) &&
	  (
	   ((m2A < -10 && nsA < 2.5 && nsA > -1.5) || (m2A > 0.16 && m2A < 0.36)) &&
	   ((m2B < -10 && nsB < 2.5 && nsB > -1.5) || (m2B > 0.16 && m2B < 0.36))
	  )
	)
      {
	Float_t eta_lTrack = lTrack.Eta();
	if(TMath::Abs(eta_lTrack) > 1.0) continue;

	Float_t InvMass_lTrack = lTrack.M();
	TVector3 vBetaPhi = -1.0*lTrack.BoostVector(); // get phi beta
	TLorentzVector lKpRest = lTrackA;
	lKpRest.Boost(vBetaPhi); // boost K+ back to phi rest frame
	TVector3 vKpRest = lKpRest.Vect().Unit(); // K+ momentum direction in phi rest frame

	float resolution = mZdcSmdCorr->GetResolution(cent9);
	float Psi = TMath::ATan2(QFull.Y(),QFull.X());
	TVector3 nQ_Full(TMath::Sin(Psi),-1.0*TMath::Cos(Psi),0.0); // normal vector of 1st Event Plane
	TVector3 nQ = nQ_Full.Unit();
	float CosThetaStar = vKpRest.Dot(nQ);
	mZdcSmdHistoManger->Fill(pt_lTrack,cent9,CosThetaStar,resolution,InvMass_lTrack,reweight,mX_flag,mMode);
      }
    }
  }

  cout << "." << flush;
  cout << " " << stop_event_use-start_event_use << "(" << 100 << "%)";
  cout << endl;
}

//-------------------------------------------------------------------
void StZdcSmdAna::Finish()
{
  mFile_OutPut->cd();
  mZdcSmdHistoManger->Write(mX_flag,mMode);
  mFile_OutPut->Close();
  mStopWatch->Stop();
  mStopWatch->Print();
}

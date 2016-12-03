#include "StEffKaon.h"
#include "TFile.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TBranch.h"
#include "TMath.h"
#include "TH1F.h"
#include "StThreeVectorF.hh"
#include <fstream>
#include <vector>
#include "../../../Utility/StSpinAlignmentCons.h"
#include "StEffCut.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "StEffHistManger.h"
#include <string>

ClassImp(StEffKaon)

int StEffKaon::mInput_flag = 1;

StEffKaon::StEffKaon(int Energy, long StartEvent, long StopEvent, int PID)
{
  mEnergy = Energy;
  mPID = PID;

  // string InPutList = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/List/Kaon_list/run14/%s_StMcEvent.list",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mParType[mPID].c_str());
  string InPutList = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/List/Kaon_list/run11/%s_StMcEvent.list",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mParType[mPID].c_str());

  SetInPutList(InPutList); // set input list

  SetStartEvent(StartEvent); // set start event
  SetStopEvent(StopEvent); // set stop event

  string OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Embedding/%s/Efficiency/Eff_%s_StMcEvent.root",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mParType[mPID].c_str(),vmsa::mBeamEnergy[mEnergy].c_str());
  SetOutPutFile(OutPutFile); // set output file

  mEffCut = new StEffCut(mEnergy);
  mEffHistManger = new StEffHistManger();
}

StEffKaon::~StEffKaon()
{
}

//------------------------------------------------------------
void StEffKaon::SetInPutList(const string inputlist)
{
  mInPutList = inputlist;
  cout << "Input list was set to: " << mInPutList.c_str() << endl;
}

void StEffKaon::SetOutPutFile(const string outputfile)
{
  mOutPutFile = outputfile;
  cout << "Output file was set to: " << mOutPutFile.c_str() << endl;
}

void StEffKaon::SetStartEvent(const long StartEvent)
{
  mStartEvent = StartEvent;
  cout << "nStartEvent = " << mStartEvent << endl;
}

void StEffKaon::SetStopEvent(const long StopEvent)
{
  mStopEvent = StopEvent;
  cout << "nStopEvent = " << mStopEvent << endl;
}
//------------------------------------------------------------

void StEffKaon::Init()
{
  mEffHistManger->InitHist();

  // initialize the TChain
  if (!mInPutList.empty())   // if input file is ok
  {
    string COUT = "Open Embedding file list ";
    cout << COUT.c_str() << mInPutList.c_str() << endl;
    ifstream in(mInPutList);  // input stream
    if(in)
    {
      cout << "file list is ok" << endl;
      mChain_Event = new TChain("eventCount");
      mChain_Track = new TChain("phi");
      char str[255];       // char array for each file name
      long entries_save_Evnet = 0;
      // long entries_save_Track = 0;
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  string addfile;
	  addfile = str;

	  mChain_Event->AddFile(addfile.c_str(),-1,"eventCount"); // event header
	  long file_entries_Evnet = mChain_Event->GetEntries();
	  // cout << "File added to data chain_Event: " << addfile.c_str() << " with " << (file_entries_Evnet-entries_save_Evnet) << " entries" << endl;
	  entries_save_Evnet = file_entries_Evnet;

	  mChain_Track->AddFile(addfile.c_str(),-1,"Kaon"); // tracks
	  // long file_entries_Track = mChain_Track->GetEntries();
	  // cout << "File added to data chain_Track: " << addfile.c_str() << " with " << (file_entries_Track-entries_save_Track) << " entries" << endl;
	  // entries_save_Track = file_entries_Track;
	}
      }
    }
    else
    {
      cout << "WARNING: file input is problemtic" << endl;
      mInput_flag = 0;
    }
  }

  // set input event header 
  if(mInput_flag == 1 && !mChain_Event->GetBranch("eventId"))
  {
    cerr << "ERROR: Could not find branch 'eventId' in tree!" << endl;
  }
  if(mInput_flag == 1 && !mChain_Track->GetBranch("McGeantId"))
  {
    cerr << "ERROR: Could not find branch 'geantId' in tree!" << endl;
  }

  if(mInput_flag == 1)
  {
    mChain_Event->SetBranchAddress("runId",&mRunId);
    mChain_Event->SetBranchAddress("eventId",&mEventId);
    mChain_Event->SetBranchAddress("mcVx",&mMcVx);
    mChain_Event->SetBranchAddress("mcVy",&mMcVy);
    mChain_Event->SetBranchAddress("mcVz",&mMcVz);
    mChain_Event->SetBranchAddress("vx",&mVx);
    mChain_Event->SetBranchAddress("vy",&mVy);
    mChain_Event->SetBranchAddress("vz",&mVz);
    mChain_Event->SetBranchAddress("vzVpd",&mVzVpd);
    mChain_Event->SetBranchAddress("centrality",&mCentrality);
    mChain_Event->SetBranchAddress("gRefMult",&mGRefMult);
    mChain_Event->SetBranchAddress("RefMult",&mRefMult);
    mChain_Event->SetBranchAddress("posRefMult",&mPosRefMult);
    mChain_Event->SetBranchAddress("negRefMult",&mNegRefMult);
    mChain_Event->SetBranchAddress("zdc",&mZdc);
    mChain_Event->SetBranchAddress("bbc",&mBbc);
    mChain_Event->SetBranchAddress("nMcTracks",&mNMcTracks);
    mChain_Event->SetBranchAddress("nRTracks",&mNRcTracks);
    mChain_Event->SetBranchAddress("magField",&mMagField);
    mChain_Event->SetBranchAddress("t0",&mT0);
    mChain_Event->SetBranchAddress("t1",&mT1);
    mChain_Event->SetBranchAddress("t2",&mT2);
    mChain_Event->SetBranchAddress("t3",&mT3);
    mChain_Event->SetBranchAddress("t4",&mT4);
    mChain_Event->SetBranchAddress("t5",&mT5);

    int num_events = mChain_Event->GetEntriesFast();
    cout << "Number of events in mChain_Event = " << num_events << endl;

    if(mStartEvent > num_events) mStartEvent = num_events;
    if(mStopEvent  > num_events) mStopEvent  = num_events;

    cout << "New nStartEvent = " << mStartEvent << ", new nStopEvent = " << mStopEvent << endl;
    mChain_Track->SetBranchAddress("McPt",&mMcPt);
    mChain_Track->SetBranchAddress("McP",&mMcP);
    mChain_Track->SetBranchAddress("McEta",&mMcEta);
    mChain_Track->SetBranchAddress("McY",&mMcY);
    mChain_Track->SetBranchAddress("McPhi",&mMcPhi);
    mChain_Track->SetBranchAddress("McGeantId",&mMcGeantId);
    mChain_Track->SetBranchAddress("eventGenLabel",&mMcEventGenLabel);
    mChain_Track->SetBranchAddress("StartVtxX",&mMcStartVtxX);
    mChain_Track->SetBranchAddress("StartVtxY",&mMcStartVtxY);
    mChain_Track->SetBranchAddress("StartVtxZ",&mMcStartVtxZ);
    mChain_Track->SetBranchAddress("StopVtxX",&mMcStopVtxX);
    mChain_Track->SetBranchAddress("StopVtxY",&mMcStopVtxY);
    mChain_Track->SetBranchAddress("StopVtxZ",&mMcStopVtxZ);

    mChain_Track->SetBranchAddress("gRcPt",&mGRcPt);
    mChain_Track->SetBranchAddress("gRcEta",&mGRcEta);
    mChain_Track->SetBranchAddress("gRcPhi",&mGRcPhi);
    mChain_Track->SetBranchAddress("gRcNfit",&mGRcNfit);
    mChain_Track->SetBranchAddress("gRcNmax",&mGRcNmax);
    mChain_Track->SetBranchAddress("gRcNcom",&mGRcNcom);
    mChain_Track->SetBranchAddress("gRcNdedx",&mGRcNdedx);
    mChain_Track->SetBranchAddress("gRcDedx",&mGRcDedx);
    mChain_Track->SetBranchAddress("gRcNsigKP",&mGRcNsigKP);
    mChain_Track->SetBranchAddress("gRcNsigKM",&mGRcNsigKM);
    mChain_Track->SetBranchAddress("gRcDca",&mGRcDca);
    mChain_Track->SetBranchAddress("gRcDcaXY",&mGRcDcaXY);
    mChain_Track->SetBranchAddress("gRcDcaZ",&mGRcDcaZ);

    mChain_Track->SetBranchAddress("pRcPt",&mPRcPt);
    mChain_Track->SetBranchAddress("pRcEta",&mPRcEta);
    mChain_Track->SetBranchAddress("pRcPhi",&mPRcPhi);
    mChain_Track->SetBranchAddress("pRcNfit",&mPRcNfit);
    mChain_Track->SetBranchAddress("pRcNmax",&mPRcNmax);
    mChain_Track->SetBranchAddress("pRcNcom",&mPRcNcom);
    mChain_Track->SetBranchAddress("pRcNdedx",&mPRcNdedx);
    mChain_Track->SetBranchAddress("pRcDedx",&mPRcDedx);
    mChain_Track->SetBranchAddress("pRcNsigKP",&mPRcNsigKP);
    mChain_Track->SetBranchAddress("pRcNsigKM",&mPRcNsigKM);
    mChain_Track->SetBranchAddress("pRcDca",&mPRcDca);
    mChain_Track->SetBranchAddress("pRcDcaXY",&mPRcDcaXY);
    mChain_Track->SetBranchAddress("pRcDcaZ",&mPRcDcaZ);

    int num_tracks = mChain_Track->GetEntriesFast();
    cout << "Number of tracks in mChain_Track = " << num_tracks<< endl;
  }

  mFile_OutPut = new TFile(mOutPutFile.c_str(),"RECREATE");
}

void StEffKaon::Make()
{
  long start_event_use = mStartEvent;
  long stop_event_use  = mStopEvent;

  mChain_Event->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry
  long NMcTraks = 0;
  long McTraks_Start = 0;
  long McTraks_Stop  = 0;
  // long NEvent = 0;

  for(long i_event = start_event_use; i_event < stop_event_use; ++i_event)
  {
    if (!mChain_Event->GetEntry(i_event)) // take the event -> information is stored in event
      break;  // end of data chunk

    if (floor(10.0*i_event/ static_cast<float>(stop_event_use)) > floor(10.0*(i_event-1)/ static_cast<float>(stop_event_use)))
      cout << "=> processing data: " << 100.0*i_event/ static_cast<float>(stop_event_use) << "%" << endl;

    if(mNMcTracks < 0) continue;
    McTraks_Start = NMcTraks;
    NMcTraks += mNMcTracks;
    McTraks_Stop = NMcTraks;
    // cout << "McTraks_Start = " << McTraks_Start << ", McTraks_Stop = " << McTraks_Stop << ", Num of McTraks = " << mNMcTracks << endl;

    McEvent EventHeader;
    EventHeader.runId = mRunId;
    EventHeader.eventId = mEventId;
    EventHeader.McVx = mMcVx;
    EventHeader.McVy = mMcVy;
    EventHeader.McVz = mMcVz;
    EventHeader.vx = mVx;
    EventHeader.vy = mVy;
    EventHeader.vz = mVz;
    EventHeader.vzVpd = mVzVpd;
    EventHeader.Centrality = mCentrality;
    EventHeader.gRefMult = mGRefMult;
    EventHeader.RefMult = mRefMult;
    EventHeader.posRefMult = mPosRefMult;
    EventHeader.negRefMult = mNegRefMult;
    EventHeader.zdc = mZdc;
    EventHeader.bbc = mBbc;
    EventHeader.nMcTracks = mNMcTracks;
    EventHeader.nRcTracks = mNRcTracks;
    EventHeader.magField = mMagField;
    EventHeader.t0 = mT0;
    EventHeader.t1 = mT1;
    EventHeader.t2 = mT2;
    EventHeader.t3 = mT3;
    EventHeader.t4 = mT4;
    EventHeader.t5 = mT5;
    if(!mEffCut->passEventCut(EventHeader)) continue;
    // NEvent++;
    for(int i_track = McTraks_Start; i_track < McTraks_Stop; ++i_track)
    {
      mChain_Track->GetEntry(i_track);

      McDecayDau McKaon; // initialize McKaon
      McKaon.McPt          = mMcPt;
      McKaon.McP           = mMcP;
      McKaon.McEta         = mMcEta;
      McKaon.McY           = mMcY;
      McKaon.McPhi         = mMcPhi;
      McKaon.McGeantId     = mMcGeantId;
      McKaon.eventGenLabel = mMcEventGenLabel;
      McKaon.StartVtxX     = mMcStartVtxX;
      McKaon.StartVtxY     = mMcStartVtxY;
      McKaon.StartVtxZ     = mMcStartVtxZ;
      McKaon.StopVtxX      = mMcStopVtxX;
      McKaon.StopVtxY      = mMcStopVtxY;
      McKaon.StopVtxZ      = mMcStopVtxZ;

      RcDecayDau RcKaon; // initialize RcKaon
      RcKaon.gRcPt     = mGRcPt;
      RcKaon.gRcEta    = mGRcEta;
      RcKaon.gRcPhi    = mGRcPhi;
      RcKaon.gRcNfit   = mGRcNfit;
      RcKaon.gRcNmax   = mGRcNmax;
      RcKaon.gRcNcom   = mGRcNcom;
      RcKaon.gRcNdedx  = mGRcNdedx;
      RcKaon.gRcDedx   = mGRcDedx;
      RcKaon.gRcNsigKP = mGRcNsigKP;
      RcKaon.gRcNsigKM = mGRcNsigKM;
      RcKaon.gRcDca    = mGRcDca;
      RcKaon.gRcDcaXY  = mGRcDcaXY;
      RcKaon.gRcDcaZ   = mGRcDcaZ;

      RcKaon.pRcPt     = mPRcPt;
      RcKaon.pRcEta    = mPRcEta;
      RcKaon.pRcPhi    = mPRcPhi;
      RcKaon.pRcNfit   = mPRcNfit;
      RcKaon.pRcNmax   = mPRcNmax;
      RcKaon.pRcNcom   = mPRcNcom;
      RcKaon.pRcNdedx  = mPRcNdedx;
      RcKaon.pRcDedx   = mPRcDedx;
      RcKaon.pRcNsigKP = mPRcNsigKP;
      RcKaon.pRcNsigKM = mPRcNsigKM;
      RcKaon.pRcDca    = mPRcDca;
      RcKaon.pRcDcaXY  = mPRcDcaXY;
      RcKaon.pRcDcaZ   = mPRcDcaZ;

      RcKaon.McPt      = mMcPt;
      RcKaon.McEta     = mMcEta;
      RcKaon.McPhi     = mMcPhi;

      //-------------------------McKaon-----------------------------------------------------
      if( !mEffCut->passTrackCut(McKaon) ) continue; // eta and momentum cuts for McKplus
      mEffHistManger->FillHistMc(EventHeader.Centrality,McKaon.McPt,McKaon.McEta,McKaon.McPhi);
      //-------------------------McKaon-----------------------------------------------------

      //-------------------------RcKaon-----------------------------------------------------
      if( !mEffCut->passTrackCut(RcKaon) ) continue; // track cuts for RcKplus
      // mEffHistManger->FillHistRc(RcKaon.RcPt,RcKaon.RcEta,RcKaon.RcPhi);
      mEffHistManger->FillHistRc(EventHeader.Centrality,McKaon.McPt,McKaon.McEta,McKaon.McPhi);
      mEffHistManger->FillHistPt(EventHeader.Centrality,McKaon.McPt,RcKaon.gRcPt,RcKaon.pRcPt);
      //-------------------------RcKaon-----------------------------------------------------
    }
  }
  mEffHistManger->CalEfficiency();
  mEffHistManger->CalEffPtEtaPhi();
  cout << "." << flush;
  cout << " " << stop_event_use << "(" << 100 << "%)";
  cout << endl;
  // cout << "Num of Events from Chain = " << mChain_Event->GetEntriesFast() << ", Num of Events: " << NEvent << ", Num of Tracks = " << NMcTraks << ", Num of Tracks from Chain = " << mChain_Track->GetEntries() << endl;
}

void StEffKaon::Finish()
{
  mFile_OutPut->cd();
  mEffHistManger->WriteHist();
  mFile_OutPut->Close();
}

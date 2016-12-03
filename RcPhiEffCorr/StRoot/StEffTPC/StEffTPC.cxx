#include "StEffTPC.h"
#include "TFile.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TBranch.h"
#include "TMath.h"
#include "TH1F.h"
#include "StThreeVectorF.hh"
#include <fstream>
#include <vector>
#include "StEffCons.h"
#include "StEffCut.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "StEffHistManger.h"
#include "TRandom3.h"

ClassImp(StEffTPC)

Int_t StEffTPC::mInput_flag = 1;

StEffTPC::StEffTPC(Int_t Energy, Long64_t StartEvent, Long64_t StopEvent, Int_t PID)
{
  mEnergy = Energy;
  mPID = PID;

  TString InPutList = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/List/%s_list/embedding_List/eff_List/%s_StMcEvent.list",Efficiency::mBeamEnergy[mEnergy].c_str(),Efficiency::mParType[mPID].c_str(),Efficiency::mParType[mPID].c_str());

  SetInPutList(InPutList); // set input list

  SetStartEvent(StartEvent); // set start event
  SetStopEvent(StopEvent); // set stop event

  TString OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Embedding/%s/Efficiency/Eff_%s_StMcEvent.root",Efficiency::mBeamEnergy[mEnergy].c_str(),Efficiency::mParType[mPID].c_str(),Efficiency::mBeamEnergy[mEnergy].c_str());
  SetOutPutFile(OutPutFile); // set output file

  mEffCut = new StEffCut(mEnergy);
  mEffHistManger = new StEffHistManger();
}

StEffTPC::~StEffTPC()
{
}

//------------------------------------------------------------
void StEffTPC::SetInPutList(const TString inputlist)
{
  mInPutList = inputlist.Copy();
  cout << "Input list was set to: " << mInPutList.Data() << endl;
}

void StEffTPC::SetOutPutFile(const TString outputfile)
{
  mOutPutFile = outputfile.Copy();
  cout << "Output file was set to: " << mOutPutFile.Data() << endl;
}

void StEffTPC::SetStartEvent(const Long64_t StartEvent)
{
  mStartEvent = StartEvent;
  cout << "nStartEvent = " << mStartEvent << endl;
}

void StEffTPC::SetStopEvent(const Long64_t StopEvent)
{
  mStopEvent = StopEvent;
  cout << "nStopEvent = " << mStopEvent << endl;
}
//------------------------------------------------------------

void StEffTPC::InitMap()
{
  TString InPutFile = "/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu200GeV/Embedding/Phi/embedding_map_pico.txt"; // embedding map pico
  cout << "Input mapping was set to: " << InPutFile.Data() << endl;
  FILE *fp = fopen(InPutFile.Data(),"r");
  if(fp == NULL)
  {
    perror("Error opening mapping file");
  }
  else
  {
    Int_t eventId;
    Float_t psiEast, psiWest;
    char line[80];

    Int_t line_counter = 0;
    while(fgets(line,80,fp))
    {
      sscanf(&line[0],"%d %f %f", &eventId, &psiEast, &psiWest);
      line_counter++;
      // cout << "eventId  = " << eventId << ", psiEast = " << psiEast << ", psiWest = " << psiWest << endl;
      mEventID[eventId] = eventId;
      mPsiEast[eventId] = psiEast;
      mPsiWest[eventId] = psiWest;
      // cout << "eventId  = " << mEventID[eventId] << ", psiEast = " << mPsiEast[eventId] << ", psiWest = " << mPsiWest[eventId] << endl;
    }
  }
}

void StEffTPC::Init()
{
  InitMap();
  mEffHistManger->InitHist();

  // initialize the TChain
  if (!mInPutList.IsNull())   // if input file is ok
  {
    TString COUT = "Open Embedding file list ";
    cout << COUT.Data() << mInPutList.Data() << endl;
    ifstream in(mInPutList);  // input stream
    if(in)
    {
      cout << "file list is ok" << endl;
      mChain_Event = new TChain("eventCount");
      mChain_Track = new TChain("phi");
      char str[255];       // char array for each file name
      Long64_t entries_save_Evnet = 0;
      // Long64_t entries_save_Track = 0;
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  TString addfile;
	  addfile = str;

	  mChain_Event->AddFile(addfile.Data(),-1,"eventCount"); // event header
	  Long64_t file_entries_Evnet = mChain_Event->GetEntries();
	  // cout << "File added to data chain_Event: " << addfile.Data() << " with " << (file_entries_Evnet-entries_save_Evnet) << " entries" << endl;
	  entries_save_Evnet = file_entries_Evnet;

	  mChain_Track->AddFile(addfile.Data(),-1,"vecMeson"); // tracks
	  // Long64_t file_entries_Track = mChain_Track->GetEntries();
	  // cout << "File added to data chain_Track: " << addfile.Data() << " with " << (file_entries_Track-entries_save_Track) << " entries" << endl;
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
  if(mInput_flag == 1 && !mChain_Track->GetBranch("geantId"))
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
    mChain_Event->SetBranchAddress("nRTracks",&mNRTracks);
    mChain_Event->SetBranchAddress("magField",&mMagField);
    mChain_Event->SetBranchAddress("t0",&mT0);
    mChain_Event->SetBranchAddress("t1",&mT1);
    mChain_Event->SetBranchAddress("t2",&mT2);
    mChain_Event->SetBranchAddress("t3",&mT3);
    mChain_Event->SetBranchAddress("t4",&mT4);
    mChain_Event->SetBranchAddress("t5",&mT5);

    Int_t num_events = mChain_Event->GetEntriesFast();
    cout << "Number of events in mChain_Event = " << num_events << endl;

    if(mStartEvent > num_events) mStartEvent = num_events;
    if(mStopEvent  > num_events) mStopEvent  = num_events;

    cout << "New nStartEvent = " << mStartEvent << ", new nStopEvent = " << mStopEvent << endl;
    mChain_Track->SetBranchAddress("pt",&mPt);
    mChain_Track->SetBranchAddress("p",&mP);
    mChain_Track->SetBranchAddress("eta",&mEta);
    mChain_Track->SetBranchAddress("y",&mY);
    mChain_Track->SetBranchAddress("phi",&mPhi);
    mChain_Track->SetBranchAddress("label",&mLabel);
    mChain_Track->SetBranchAddress("geantId",&mGeantId);
    mChain_Track->SetBranchAddress("kpPt",&mKpPt);
    mChain_Track->SetBranchAddress("kpEta",&mKpEta);
    mChain_Track->SetBranchAddress("kpPhi",&mKpPhi);
    mChain_Track->SetBranchAddress("kpGeantId",&mKpGeantId);
    mChain_Track->SetBranchAddress("kpStartVtxX",&mKpStartVtxX);
    mChain_Track->SetBranchAddress("kpStartVtxY",&mKpStartVtxY);
    mChain_Track->SetBranchAddress("kpStartVtxZ",&mKpStartVtxZ);
    mChain_Track->SetBranchAddress("kpStopVtxX",&mKpStopVtxX);
    mChain_Track->SetBranchAddress("kpStopVtxY",&mKpStopVtxY);
    mChain_Track->SetBranchAddress("kpStopVtxZ",&mKpStopVtxZ);
    mChain_Track->SetBranchAddress("kpRpt",&mKpRpt);
    mChain_Track->SetBranchAddress("kpReta",&mKpReta);
    mChain_Track->SetBranchAddress("kpRphi",&mKpRphi);
    mChain_Track->SetBranchAddress("kpNfit",&mKpNfit);
    mChain_Track->SetBranchAddress("kpNmax",&mKpNmax);
    mChain_Track->SetBranchAddress("kpNcom",&mKpNcom);
    mChain_Track->SetBranchAddress("kpNdedx",&mKpNdedx);
    mChain_Track->SetBranchAddress("kpDedx",&mKpDedx);
    mChain_Track->SetBranchAddress("kpNsigKP",&mKpNsigKP);
    mChain_Track->SetBranchAddress("kpNsigKM",&mKpNsigKM);
    mChain_Track->SetBranchAddress("kpDca",&mKpDca);
    mChain_Track->SetBranchAddress("kpDcaXY",&mKpDcaXY);
    mChain_Track->SetBranchAddress("kpDcaZ",&mKpDcaZ);
    mChain_Track->SetBranchAddress("kmPt",&mKmPt);
    mChain_Track->SetBranchAddress("kmEta",&mKmEta);
    mChain_Track->SetBranchAddress("kmPhi",&mKmPhi);
    mChain_Track->SetBranchAddress("kmGeantId",&mKmGeantId);
    mChain_Track->SetBranchAddress("kmStartVtxX",&mKmStartVtxX);
    mChain_Track->SetBranchAddress("kmStartVtxY",&mKmStartVtxY);
    mChain_Track->SetBranchAddress("kmStartVtxZ",&mKmStartVtxZ);
    mChain_Track->SetBranchAddress("kmStopVtxX",&mKmStopVtxX);
    mChain_Track->SetBranchAddress("kmStopVtxY",&mKmStopVtxY);
    mChain_Track->SetBranchAddress("kmStopVtxZ",&mKmStopVtxZ);
    mChain_Track->SetBranchAddress("kmRpt",&mKmRpt);
    mChain_Track->SetBranchAddress("kmReta",&mKmReta);
    mChain_Track->SetBranchAddress("kmRphi",&mKmRphi);
    mChain_Track->SetBranchAddress("kmNfit",&mKmNfit);
    mChain_Track->SetBranchAddress("kmNmax",&mKmNmax);
    mChain_Track->SetBranchAddress("kmNcom",&mKmNcom);
    mChain_Track->SetBranchAddress("kmNdedx",&mKmNdedx);
    mChain_Track->SetBranchAddress("kmDedx",&mKmDedx);
    mChain_Track->SetBranchAddress("kmNsigKP",&mKmNsigKP);
    mChain_Track->SetBranchAddress("kmNsigKM",&mKmNsigKM);
    mChain_Track->SetBranchAddress("kmDca",&mKmDca);
    mChain_Track->SetBranchAddress("kmDcaXY",&mKmDcaXY);
    mChain_Track->SetBranchAddress("kmDcaZ",&mKmDcaZ);
    mChain_Track->SetBranchAddress("invMass",&mInvMass);
    mChain_Track->SetBranchAddress("rPt",&mRPt);
    mChain_Track->SetBranchAddress("rEta",&mREta);
    mChain_Track->SetBranchAddress("rY",&mRY);
    mChain_Track->SetBranchAddress("rphi",&mRphi);

    Int_t num_tracks = mChain_Track->GetEntriesFast();
    cout << "Number of tracks in mChain_Track = " << num_tracks<< endl;
  }

  mFile_OutPut = new TFile(mOutPutFile.Data(),"RECREATE");
}

void StEffTPC::Make()
{
  Long64_t start_event_use = mStartEvent;
  Long64_t stop_event_use  = mStopEvent;

  mChain_Event->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry
  Long64_t NMcTraks = 0;
  Long64_t McTraks_Start = 0;
  Long64_t McTraks_Stop  = 0;
  // Long64_t NEvent = 0;

  for(Long64_t i_event = start_event_use; i_event < stop_event_use; ++i_event)
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
    EventHeader.mcVx = mMcVx;
    EventHeader.mcVy = mMcVy;
    EventHeader.mcVz = mMcVz;
    EventHeader.vx = mVx;
    EventHeader.vy = mVy;
    EventHeader.vz = mVz;
    EventHeader.vzVpd = mVzVpd;
    EventHeader.centrality = mCentrality;
    EventHeader.gRefMult = mGRefMult;
    EventHeader.RefMult = mRefMult;
    EventHeader.posRefMult = mPosRefMult;
    EventHeader.negRefMult = mNegRefMult;
    EventHeader.zdc = mZdc;
    EventHeader.bbc = mBbc;
    EventHeader.nMcTracks = mNMcTracks;
    EventHeader.nRTracks = mNRTracks;
    EventHeader.magField = mMagField;
    EventHeader.t0 = mT0;
    EventHeader.t1 = mT1;
    EventHeader.t2 = mT2;
    EventHeader.t3 = mT3;
    EventHeader.t4 = mT4;
    EventHeader.t5 = mT5;
    if(!mEventID[EventHeader.eventId]) continue;
    if(!mEffCut->passEventCut(EventHeader)) continue;
    // cout << "eventId = " << EventHeader.eventId << ", mEventID = " << mEventID[EventHeader.eventId] << endl;
    // NEvent++;
    float Psi2East = mPsiEast[EventHeader.eventId];
    float Psi2West = mPsiWest[EventHeader.eventId];
    for(Int_t i_track = McTraks_Start; i_track < McTraks_Stop; ++i_track)
    {
      mChain_Track->GetEntry(i_track);

      McVecMeson McPhi; // initialize McPhi
      McPhi.pt = mPt;
      McPhi.p = mP;
      McPhi.eta = mEta;
      McPhi.y = mY;
      McPhi.phi = mPhi;
      McPhi.label = mLabel;
      McPhi.geantId = mGeantId;

      McDecayDau McKP; // initialize McKplus
      McKP.pt = mKpPt;
      McKP.eta = mKpEta;
      McKP.phi = mKpPhi;
      McKP.geantId = mKpGeantId;
      McKP.startVtxX = mKpStartVtxX;
      McKP.startVtxY = mKpStartVtxY;
      McKP.startVtxZ = mKpStartVtxZ;
      McKP.stopVtxX = mKpStopVtxX;
      McKP.stopVtxY = mKpStopVtxY;
      McKP.stopVtxZ = mKpStopVtxZ;

      RcDecayDau RcKP; // initialize RcKplus
      RcKP.pt = mKpRpt;
      RcKP.eta = mKpReta;
      RcKP.phi = mKpRphi;
      RcKP.nFit = mKpNfit;
      RcKP.nMax = mKpNmax;
      RcKP.nCom = mKpNcom;
      RcKP.nDedx = mKpNdedx;
      RcKP.Dedx = mKpDedx;
      RcKP.nSigKP = mKpNsigKP;
      RcKP.nSigKM = mKpNsigKM;
      RcKP.dca = mKpDca;
      RcKP.dcaXY = mKpDcaXY;
      RcKP.dcaZ = mKpDcaZ;
      RcKP.McPt = mKpPt;
      RcKP.McEta = mKpEta;
      RcKP.McPhi = mKpPhi;

      McDecayDau McKM; // initialize McKminus
      McKM.pt = mKmPt;
      McKM.eta = mKmEta;
      McKM.phi = mKmPhi;
      McKM.geantId = mKmGeantId;
      McKM.startVtxX = mKmStartVtxX;
      McKM.startVtxY = mKmStartVtxY;
      McKM.startVtxZ = mKmStartVtxZ;
      McKM.stopVtxX = mKmStopVtxX;
      McKM.stopVtxY = mKmStopVtxY;
      McKM.stopVtxZ = mKmStopVtxZ;

      RcDecayDau RcKM; // initialize RcKminus
      RcKM.pt = mKmRpt;
      RcKM.eta = mKmReta;
      RcKM.phi = mKmRphi;
      RcKM.nFit = mKmNfit;
      RcKM.nMax = mKmNmax;
      RcKM.nCom = mKmNcom;
      RcKM.nDedx = mKmNdedx;
      RcKM.Dedx = mKmDedx;
      RcKM.nSigKP = mKmNsigKP;
      RcKM.nSigKM = mKmNsigKM;
      RcKM.dca = mKmDca;
      RcKM.dcaXY = mKmDcaXY;
      RcKM.dcaZ = mKmDcaZ;
      RcKM.McPt = mKmPt;
      RcKM.McEta = mKmEta;
      RcKM.McPhi = mKmPhi;

      RcVecMeson RcPhi; // initialize RcPhi
      RcPhi.InvMass = mInvMass;
      RcPhi.pt = mRPt;
      RcPhi.eta = mREta;
      RcPhi.y = mRY;
      RcPhi.phi = mRphi;
      RcPhi.McPt = mPt;
      RcPhi.McEta = mEta;
      RcPhi.McPhi = mPhi;

      //-------------------------McPhi-----------------------------------------------------
      if( !mEffCut->passTrackCutPhi(McPhi) ) continue; // eta cuts for McPhi [-1.0,1.0]

      if( !mEffCut->passTrackCut(McKP) ) continue; // eta and momentum cuts for McKplus
      if( !mEffCut->passTrackCut(McKM) ) continue; // eta and momentum cuts for McKminus

      TVector3 vMcKPBoost = mEffHistManger->CalBoostedVector(McKP,McPhi); // K+ momentum in phi rest frame

      // calculate cos(theta*)
      TVector3 nQ;
      float Psi2 = -999.0;
      if(mEffCut->passPhiEtaEast(McKP)) // K+ neg eta(east)
      { // Below is West Only
	TVector3 nQ_West(TMath::Sin(Psi2West),-1.0*TMath::Cos(Psi2West),0.0); // direction of angular momentum
	nQ = nQ_West.Unit();
	Psi2 = Psi2West;
      }
      if(mEffCut->passPhiEtaWest(McKP)) // K+ pos eta (west)
      { // Below is East Only
	TVector3 nQ_East(TMath::Sin(Psi2East),-1.0*TMath::Cos(Psi2East),0.0); // direction of angular momentum
	nQ = nQ_East.Unit();
	Psi2 = Psi2East;
      }
      if(Psi2 < -100.0) continue;
      float McCosThetaStar = vMcKPBoost.Dot(nQ);
      // float McCosThetaStar = TMath::Sin(vMcKPBoost.Theta())*TMath::Sin(Psi2-vMcKPBoost.Phi());
      /*
      float Psi2 = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
      TVector3 QVector(TMath::Sin(Psi2),-1.0*TMath::Cos(Psi2),0.0); // direction of angular momentum
      TVector3 nQ = QVector.Unit();
      float McCosThetaStar = vMcKPBoost.Dot(nQ);
      */
      mEffHistManger->FillHistMc(EventHeader.centrality,McPhi.pt,McPhi.eta,McPhi.phi,McCosThetaStar);
      mEffHistManger->FillQAEPMc(EventHeader.centrality,McPhi.pt,McPhi.phi-Psi2);
      // mEffHistManger->FillHistMc(EventHeader.centrality,McKP.pt,McKP.eta,McKP.phi,McCosThetaStar);
      // mEffHistManger->FillHistMc(EventHeader.centrality,McKM.pt,McKM.eta,McKM.phi,McCosThetaStar);
      // mEffHistManger->FillQACosMc(EventHeader.centrality,McPhi,McKP,McKM);
      //-------------------------McPhi-----------------------------------------------------

      //-------------------------RcPhi-----------------------------------------------------
      if( !mEffCut->passTrackCut(RcKP) ) continue; // track cuts for RcKplus
      if( !mEffCut->passTrackCut(RcKM) ) continue; // track cuts for RcKminus
      if( !mEffCut->passTrackCutPhi(RcPhi) ) continue;  // eta cuts for RcPhi [-1.0,1.0]

      // mEffHistManger->FillHistRc(RcPhi.pt,RcPhi.eta,RcPhi.phi,RcCosThetaStar);
      mEffHistManger->FillHistRc(EventHeader.centrality,McPhi.pt,McPhi.eta,McPhi.phi,McCosThetaStar);
      mEffHistManger->FillQAEPRc(EventHeader.centrality,McPhi.pt,McPhi.phi-Psi2);
      // mEffHistManger->FillHistRc(EventHeader.centrality,McKP.pt,McKP.eta,McKP.phi,McCosThetaStar);
      // mEffHistManger->FillHistRc(EventHeader.centrality,McKM.pt,McKM.eta,McKM.phi,McCosThetaStar);
      // mEffHistManger->FillQACosRc(EventHeader.centrality,McPhi,McKP,McKM);
      //-------------------------RcPhi-----------------------------------------------------
    }
  }
  // mEffHistManger->CalEfficiency();
  // mEffHistManger->CalEffPtEtaPhi();
  mEffHistManger->CalEffCosThetaStar();
  mEffHistManger->CalEffEP();
  cout << "." << flush;
  cout << " " << stop_event_use << "(" << 100 << "%)";
  cout << endl;
  // cout << "Num of Events from Chain = " << mChain_Event->GetEntriesFast() << ", Num of Events: " << NEvent << ", Num of Tracks = " << NMcTraks << ", Num of Tracks from Chain = " << mChain_Track->GetEntries() << endl;
}

void StEffTPC::Finish()
{
  mFile_OutPut->cd();
  mEffHistManger->WriteHist();
  mFile_OutPut->Close();
}

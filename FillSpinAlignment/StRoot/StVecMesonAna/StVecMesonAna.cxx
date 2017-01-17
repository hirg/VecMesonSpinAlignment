#include "StRoot/StVecMesonAna/StVecMesonAna.h"
#include "StRoot/StVecMesonAna/StVecMesonCut.h"
#include "StRoot/StVecMesonAna/StVecMesonCorr.h"
#include "StRoot/StVecMesonAna/StVecMesonHistoManger.h"
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

ClassImp(StVecMesonAna)

StRefMultCorr* StVecMesonAna::mRefMultCorr = NULL;
Int_t StVecMesonAna::mInPut_flag = 1;
char* StVecMesonAna::VM_EVENT_TREE = NULL;
char* StVecMesonAna::VM_EVENT_BRANCH = NULL;

//----------------------------------------------------
StVecMesonAna::StVecMesonAna(Int_t energy, Int_t X_flag, Int_t List, Long64_t start_event, Long64_t stop_event, Int_t mode)
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
  mVecMesonCorr = new StVecMesonCorr(mEnergy);
  mVecMesonCut = new StVecMesonCut(mEnergy);
  mVecMesonHistoManger = new StVecMesonHistoManger();
}

StVecMesonAna::~StVecMesonAna()
{
}
//----------------------------------------------------
// set Input/Output
void StVecMesonAna::setInputDir(const TString inputdir)
{
  mInputdir = inputdir.Copy();
  cout << "Input directory was set to: " << mInputdir.Data() << endl;
}
void StVecMesonAna::setOutputfile(const TString outputfile)
{
  mOutputfile = outputfile.Copy();
  cout << "Output file was set to: " << mOutputfile.Data() << endl;
}
void StVecMesonAna::setInPutList(const TString iInPutList)
{
  mInPutList = iInPutList.Copy();
  TString InFo_InPutList = Form("InPut %s list was set to: %s",vmsa::MixEvent[mX_flag].Data(),mInPutList.Data());
  cout << InFo_InPutList.Data() << endl;
}
void StVecMesonAna::setStopEvent(const Long64_t StopEvent)
{
    mStopEvent = StopEvent;
    cout << "nStopEvent = " << mStopEvent << endl;
}
void StVecMesonAna::setStartEvent(const Long64_t StartEvent)
{
    mStartEvent = StartEvent;
    cout << "nStartEvent = " << mStartEvent << endl;
}
//----------------------------------------------------
// initial functions
void StVecMesonAna::Init()
{
  mVecMesonCorr->InitReCenterCorrection();
  mVecMesonCorr->InitShiftCorrection();
  mVecMesonCorr->InitResolutionCorr();
  mVecMesonHistoManger->Init(mX_flag,mMode);

  TString inputdir = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/Forest/",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mPID[mMode].c_str());
  setInputDir(inputdir);

  const Int_t list_start = vmsa::mList_Delta*mList + 1; // start list
  const Int_t list_stop  = vmsa::mList_Delta*(mList+1); // stop list

  TString InPutList = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/List/Split_%s_%s_%d_%d.list",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mPID[mMode].c_str(),vmsa::MixEvent[mX_flag].Data(),vmsa::mBeamEnergy[mEnergy].c_str(),list_start,list_stop);
  setInPutList(InPutList);

  TString outputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/Yields/Yields_%s_%s_%d.root",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mPID[mMode].c_str(),vmsa::MixEvent[mX_flag].Data(),vmsa::mBeamEnergy[mEnergy].c_str(),mList);
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

    Int_t num_events = mInPut->GetEntriesFast();
    cout << "Number of events in file(s) = " << num_events << endl;
    if(mStartEvent > num_events) mStartEvent = num_events;
    if(mStopEvent > num_events) mStopEvent   = num_events;

    cout << "New nStartEvent = " << mStartEvent << ", new nStopEvent = " << mStopEvent << endl;
  }
}

void StVecMesonAna::Make()
{
  if(mMode == 0) MakePhi();
}


// loop phi meson Same Event
void StVecMesonAna::MakePhi()
{
  Long64_t start_event_use;
  Long64_t stop_event_use;

  start_event_use = mStartEvent;
  stop_event_use  = mStopEvent;
  mInPut->SetBranchAddress( VM_EVENT_BRANCH, &mPhiMeson_event);
  mInPut->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

  // Initialise Event Head
  StThreeVectorF PrimaryVertex(0.0,0.0,0.0);
  Int_t          RunId = 0;
  Int_t          EventId = 0;
  Int_t          RefMult = 0;
  Int_t          Centrality = 0;
  Int_t          N_prim = 0;
  Int_t          N_non_prim = 0;
  Int_t          N_Tof_match = 0;
  Float_t        ZDCx = 0.0; 
  Float_t        BBCx = 0.0; 
  Float_t        VzVpd = 0.0;
  Int_t          NumTrackUsed = 0;
  // ---------------------------------------QVector---------------------------------------------
  TVector2 Q2East(0.0,0.0);
  TVector2 Q2West(0.0,0.0);
  TVector2 Q2Full(0.0,0.0);
  // -----------------------------------Number of Tracks----------------------------------------
  Int_t   NumTrackEast = 0;
  Int_t   NumTrackWest = 0;
  Int_t   NumTrackFull = 0;
  Int_t   NumTrackFullEast = 0;
  Int_t   NumTrackFullWest = 0;

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
    Q2East           = mPhiMeson_event->getQ2East();
    Q2West           = mPhiMeson_event->getQ2West();
    Q2Full           = mPhiMeson_event->getQ2Full();
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
    Int_t flagA = -1;
    Int_t flagB = -1;

    // vz sign
    Int_t vz_sign;
    if(PrimaryVertex.z() > 0.0)
    {
      vz_sign = 0;
    }
    else
    {
      vz_sign = 1;
    }

    // Centrality
    mRefMultCorr->init(RunId);
    if(mEnergy == 6) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z(),ZDCx); // 200 GeV
    if(mEnergy != 6) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z());
    const Int_t cent9 = Centrality;
    const Double_t reweight = mRefMultCorr->getWeight();

    // runIndex
    mRunIdEventsDb = StRunIdEventsDb::Instance(vmsa::mEnergyValue[mEnergy],vmsa::mBeamYear[mEnergy]);
    const Int_t runIndex = mRunIdEventsDb->getRunIdIndex(RunId); // expensive
    // cout << runIndex << endl;

    if (counter != 0  &&  counter % 1000 == 0)
      cout << "." << flush;
    if (counter != 0  &&  counter % 10000 == 0)
    {
      if((stop_event_use-start_event_use) > 0)
      {
	Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
	cout << " " << counter-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (VecMesonSpinAlignment) " << flush;
      }
    }

    // get Track Information
    if(mVecMesonCorr->passTrackEtaNumCut(NumTrackEast,NumTrackWest))
    {
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

	// if( // always require ToF for daughter particles
	//     ( m2A > 0 && ((fabs(pA) < 1.5 && m2A > 0.16 && m2A < 0.36) || (fabs(pA) >= 1.5 && m2A > 0.125 && m2A < 0.36)) ) &&
	//     ( m2B > 0 && ((fabs(pB) < 1.5 && m2B > 0.16 && m2B < 0.36) || (fabs(pB) >= 1.5 && m2B > 0.125 && m2B < 0.36)) )
	//   )
	if( // previous cut with TPC+ToF (if possible) at low momentum and TPC+ToF (always) at high momentum 
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

	  if(mVecMesonCut->passPhiEtaEast(lTrackA)) // K+ neg eta(east)
	  { // Below is West Only
	    TVector2 Q2Vector = Q2West;
	    // subtract auto-correlation from pos eta(west) event plane
	    if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaWest(lTrackB)) // trackB
	    {
	      Float_t  w = mVecMesonCorr->getWeight(lTrackB);
	      TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
	      TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_West(cent9,runIndex,vz_sign);
	      Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
	    }
	    Float_t Res2 = mVecMesonCorr->getResolution2_EP(cent9);
	    Float_t Psi2_west = mVecMesonCorr->calShiftAngle2West_EP(Q2Vector,runIndex,cent9,vz_sign);
	    TVector3 nQ_West(TMath::Sin(Psi2_west),-1.0*TMath::Cos(Psi2_west),0.0); // normal vector of 2nd Event Plane
	    TVector3 nQ = nQ_West.Unit();
	    Double_t CosThetaStar = vKpRest.Dot(nQ);

	    mVecMesonHistoManger->Fill(pt_lTrack,cent9,CosThetaStar,Res2,InvMass_lTrack,reweight,mX_flag,mMode);
	  }

	  if(mVecMesonCut->passPhiEtaWest(lTrackA)) // K+ pos eta (west)
	  { // Below is East Only
	    TVector2 Q2Vector = Q2East;
	    // subtract auto-correlation from pos eta(west) event plane
	    if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaEast(lTrackB)) // trackB
	    {
	      Float_t  w = mVecMesonCorr->getWeight(lTrackB);
	      TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
	      TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_East(cent9,runIndex,vz_sign);
	      Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
	    }
	    Float_t Res2 = mVecMesonCorr->getResolution2_EP(cent9);
	    Float_t Psi2_east = mVecMesonCorr->calShiftAngle2East_EP(Q2Vector,runIndex,cent9,vz_sign);
	    TVector3 nQ_East(TMath::Sin(Psi2_east),-1.0*TMath::Cos(Psi2_east),0.0); // normal vector of 2nd Event Plane
	    TVector3 nQ = nQ_East.Unit();
	    Double_t CosThetaStar = vKpRest.Dot(nQ);

	    mVecMesonHistoManger->Fill(pt_lTrack,cent9,CosThetaStar,Res2,InvMass_lTrack,reweight,mX_flag,mMode);
	  }
	}
      }
    }
  }

  cout << "." << flush;
  cout << " " << stop_event_use-start_event_use << "(" << 100 << "%)";
  cout << endl;
}

//-------------------------------------------------------------------
void StVecMesonAna::Finish()
{
  mFile_OutPut->cd();
  mVecMesonHistoManger->Write(mX_flag,mMode);
  mFile_OutPut->Close();
}

#include "StStrangenessAna.h"
#include "StStrangenessCons.h"
#include "StStrangenessCut.h"
#include "StStrangenessCorr.h"
#include "StStrangenessHistoManger.h"
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

ClassImp(StStrangenessAna)

StRefMultCorr* StStrangenessAna::mRefMultCorr = NULL;
Int_t StStrangenessAna::mInPut_flag = 1;
char* StStrangenessAna::XUV0_EVENT_TREE = NULL;
char* StStrangenessAna::XUV0_EVENT_BRANCH = NULL;

//----------------------------------------------------
StStrangenessAna::StStrangenessAna(Int_t energy, Int_t X_flag, Int_t List, Long64_t start_event, Long64_t stop_event, Int_t mode, Int_t flag_Embedding)
{
  mEnergy = energy;
  mX_flag = X_flag;
  mList = List;
  mStart_Event = start_event;
  mStop_Event = stop_event;
  mMode = mode;
  mFlag_Embedding = flag_Embedding;
  if(!mRefMultCorr)
  {
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
  mStrangenessCorr = new StStrangenessCorr();
  mStrangenessCut = new StStrangenessCut(mEnergy);
  mStrangenessHistoManger = new StStrangenessHistoManger();
}

StStrangenessAna::~StStrangenessAna()
{
}
//----------------------------------------------------
// set Input/Output
void StStrangenessAna::setInputDir(const TString inputdir)
{
  mInputdir = inputdir.Copy();
  cout << "Input directory was set to: " << mInputdir.Data() << endl;
}
void StStrangenessAna::setOutputfile(const TString outputfile)
{
  mOutputfile = outputfile.Copy();
  cout << "Output file was set to: " << mOutputfile.Data() << endl;
}
void StStrangenessAna::setInPutList(const TString iInPutList)
{
  mInPutList = iInPutList.Copy();
  TString InFo_InPutList = Form("InPut %s list was set to: %s",Strangeness::XFlag[mX_flag].Data(),mInPutList.Data());
  cout << InFo_InPutList.Data() << endl;
}
void StStrangenessAna::setStopEvent(const Long64_t StopEvent)
{
    mStopEvent = StopEvent;
    cout << "nStopEvent = " << mStopEvent << endl;
}
void StStrangenessAna::setStartEvent(const Long64_t StartEvent)
{
    mStartEvent = StartEvent;
    cout << "nStartEvent = " << mStartEvent << endl;
}
//----------------------------------------------------
// initial functions
void StStrangenessAna::Init()
{
  mStrangenessCorr->InitReCenterCorrection(mEnergy);
  mStrangenessCorr->InitShiftCorrection(mEnergy);
  mStrangenessHistoManger->Init(mX_flag,mMode);

  TString inputdir;
  // if(mFlag_Embedding == 0) inputdir = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/",Strangeness::Energy[mEnergy].Data(),Strangeness::Partype[mMode].Data());
  if(mFlag_Embedding == 0) inputdir = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/%s/",Strangeness::Energy[mEnergy].Data(),Strangeness::Partype[mMode].Data());
  if(mFlag_Embedding == 1) inputdir = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Embedding/%s/",Strangeness::Energy[mEnergy].Data(),Strangeness::Partype[mMode].Data());
  setInputDir(inputdir);

  const Int_t list_start = Strangeness::mList_Delta*mList + 1; // start list
  const Int_t list_stop  = Strangeness::mList_Delta*(mList+1); // stop list

  TString InPutList;
  if(mFlag_Embedding == 0) InPutList = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/List/%s_list/spin_List/Split_%s_%s_%d_%d.list",Strangeness::Energy[mEnergy].Data(),Strangeness::Partype[mMode].Data(),Strangeness::XFlag[mX_flag].Data(),Strangeness::Energy[mEnergy].Data(),list_start,list_stop);
  if(mFlag_Embedding == 1) InPutList = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/List/%s_list/embedding_List/embedding_map_pico.list",Strangeness::Energy[mEnergy].Data(),Strangeness::Partype[mMode].Data());
  setInPutList(InPutList);

  TString outputfile;
  if(mFlag_Embedding == 0) outputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/Yields/Yields_%s_%s_%d.root",Strangeness::Energy[mEnergy].Data(),Strangeness::Partype[mMode].Data(),Strangeness::XFlag[mX_flag].Data(),Strangeness::Energy[mEnergy].Data(),mList);
  if(mFlag_Embedding == 1) outputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Embedding/%s/Yields/Yields_%s_%s_%d.root",Strangeness::Energy[mEnergy].Data(),Strangeness::Partype[mMode].Data(),Strangeness::XFlag[mX_flag].Data(),Strangeness::Energy[mEnergy].Data(),mList);
  setOutputfile(outputfile);

  setStartEvent(Long64_t(mStart_Event));
  setStopEvent(Long64_t(mStop_Event));
  //----------------------------------------------------------------------------------------------------

  TString Notification = Form("Initializing parameters and input/output for %s %s",Strangeness::Partype[mMode].Data(),Strangeness::XFlag[mX_flag].Data());
  cout << Notification.Data() << endl;
  mFile_OutPut = new TFile(mOutputfile.Data(),"RECREATE");

  XUV0_EVENT_TREE       = (char*)Strangeness::v0_tree[mMode].Data();
  XUV0_EVENT_BRANCH     = (char*)Strangeness::v0_branch[mMode].Data();

  //----------------------------------------------------------------------------------------------------
  // input
  if (!mInPutList.IsNull())   // if input file is ok
  {
    TString InFo_List = Form("Open %s file list ",Strangeness::XFlag[mX_flag].Data());
    cout << InFo_List.Data() << endl;
    ifstream in(mInPutList);  // input stream
    if(in)
    {
      cout << "input file list is ok" << endl;
      mInPut = new TChain( XUV0_EVENT_TREE, XUV0_EVENT_TREE );
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
	  mInPut->AddFile(addfile.Data(),-1, XUV0_EVENT_TREE );
	  Long64_t file_entries = mInPut->GetEntries();
	  cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
	  entries_save = file_entries;
	}
      }
    }
    else
    {
      TString InFo_Warning = Form("WARNING: %s file input is problemtic",Strangeness::XFlag[mX_flag].Data());
      cout << InFo_Warning.Data() << endl;
      mInPut_flag = 0;
    }
  }

  // Set the input tree
  if (mInPut_flag == 1 && !mInPut->GetBranch( XUV0_EVENT_BRANCH ))
  {
    cerr << "ERROR: Could not find branch '"
      << XUV0_EVENT_BRANCH << "'in tree!" << endl;
  }

  if(mMode == 0) mXuPhiMeson_event = new StAlexPhiMesonEvent();

  if(mInPut_flag == 1)
  {
    if(mMode == 0) mInPut->SetBranchAddress( XUV0_EVENT_BRANCH, &mXuPhiMeson_event );

    Int_t num_events = mInPut->GetEntriesFast();
    cout << "Number of events in file(s) = " << num_events << endl;
    if(mStartEvent > num_events) mStartEvent = num_events;
    if(mStopEvent > num_events) mStopEvent   = num_events;

    cout << "New nStartEvent = " << mStartEvent << ", new nStopEvent = " << mStopEvent << endl;
  }
}

void StStrangenessAna::Make()
{
  if(mMode == 0) MakePhi();
}


// loop phi meson Same Event
void StStrangenessAna::MakePhi()
{
  ofstream out_map;
  if(mFlag_Embedding == 1)
  {
    TString OutMap = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Embedding/%s/embedding_map_pico.txt",Strangeness::Energy[mEnergy].Data(),Strangeness::Partype[mMode].Data());
    out_map.open(OutMap);
  }

  Long64_t start_event_use;
  Long64_t stop_event_use;

  start_event_use = mStartEvent;
  stop_event_use  = mStopEvent;
  mInPut->SetBranchAddress( XUV0_EVENT_BRANCH, &mXuPhiMeson_event);
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
  TVector2       Q2East[4];
  TVector2       Q2West[4];
  TVector2       Q3East[4];
  TVector2       Q3West[4];
  // -----------------------------------Number of Tracks----------------------------------------
  Int_t          NumTrackEast[4];
  Int_t          NumTrackWest[4];
  for(Int_t j = 0; j < 4; j++)
  {
    Q2East[j].Set(0.0,0.0);
    Q2West[j].Set(0.0,0.0);
    Q3East[j].Set(0.0,0.0);
    Q3West[j].Set(0.0,0.0);
    NumTrackEast[j] = 0;
    NumTrackWest[j] = 0;
  }

  for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
  {
    if (!mInPut->GetEntry( counter )) // take the event -> information is stored in event
      break;  // end of data chunk

    // get Event Header
    PrimaryVertex     = mXuPhiMeson_event->getPrimaryVertex();
    RunId             = mXuPhiMeson_event->getRunId();
    EventId           = mXuPhiMeson_event->getEventId();
    RefMult           = mXuPhiMeson_event->getRefMult();
    Centrality        = mXuPhiMeson_event->getCentrality();
    N_prim            = mXuPhiMeson_event->getN_prim();
    N_non_prim        = mXuPhiMeson_event->getN_non_prim();
    N_Tof_match       = mXuPhiMeson_event->getN_Tof_match();
    ZDCx              = mXuPhiMeson_event->getZDCx(); 
    BBCx              = mXuPhiMeson_event->getBBCx(); 
    VzVpd             = mXuPhiMeson_event->getVzVpd();
    NumTrackUsed      = mXuPhiMeson_event->getNumTracks();

    for(Int_t j = 0; j < 4; j++)
    {
      Q2East[j]       = mXuPhiMeson_event->getQ2East(j);
      Q2West[j]       = mXuPhiMeson_event->getQ2West(j);
      Q3East[j]       = mXuPhiMeson_event->getQ3East(j);
      Q3West[j]       = mXuPhiMeson_event->getQ3West(j);
      NumTrackEast[j] = mXuPhiMeson_event->getNumTrackEast(j);
      NumTrackWest[j] = mXuPhiMeson_event->getNumTrackWest(j);
    }
 
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
    if(mEnergy == 0) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z(),ZDCx);
    if(mEnergy != 0) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z());
    const Int_t cent9 = Centrality;
    const Double_t reweight = mRefMultCorr->getWeight();

    // runIndex
    mRunIdEventsDb = StRunIdEventsDb::Instance(Strangeness::mBeamEnergy[mEnergy],Strangeness::mBeamYear[mEnergy]);
    const Int_t runIndex = mRunIdEventsDb->getRunIdIndex(RunId); // expensive
    // cout << runIndex << endl;

    if (counter != 0  &&  counter % 1000 == 0)
      cout << "." << flush;
    if (counter != 0  &&  counter % 10000 == 0)
    {
      if((stop_event_use-start_event_use) > 0)
      {
	Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
	cout << " " << counter-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (strangeness_SpinAlignment) " << flush;
      }
    }

    if(mFlag_Embedding == 1)
    {
      Float_t EP_east = mStrangenessCorr->calShiftAngle2East_EP(Q2East[0],runIndex,cent9,vz_sign,0); // eta_gap = +/-0.05
      Float_t EP_west = mStrangenessCorr->calShiftAngle2West_EP(Q2West[0],runIndex,cent9,vz_sign,0);
      // cout << "eventId = " << EventId << ", Psi2_east = "<< Psi2_east << endl;
      out_map << EventId << "  " << EP_east << "  " << EP_west << endl;
    }

    // get Track Information
    for(Int_t j = Strangeness::EtaGap_start; j < Strangeness::EtaGap_stop; j++)
    {
      if(mStrangenessCorr->passTrackNumCut(NumTrackEast[j],NumTrackWest[j]))
      {
	for(UShort_t nTracks = 0; nTracks < NumTrackUsed; nTracks++) // loop over all tracks of the actual event
	{
	  mXuPhiMeson_track = mXuPhiMeson_event->getTrack(nTracks);
	  m2A = mXuPhiMeson_track->getMass2A();
	  m2B = mXuPhiMeson_track->getMass2B();
	  nsA = mXuPhiMeson_track->getNSigKaonA();
	  nsB = mXuPhiMeson_track->getNSigKaonB();
	  dcaA = mXuPhiMeson_track->getDcaA();
	  dcaB = mXuPhiMeson_track->getDcaB();
	  lTrackA = mXuPhiMeson_track->getTrackA();
	  lTrackB = mXuPhiMeson_track->getTrackB();
	  flagA = mXuPhiMeson_track->getFlagA();
	  flagB = mXuPhiMeson_track->getFlagB();

	  Float_t pA = lTrackA.P();
	  Float_t pB = lTrackB.P();
	  TLorentzVector lTrack = lTrackA + lTrackB;
	  Float_t pt_lTrack = lTrack.Perp();

	  if(
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
	    TVector3 beta_lTrack = -1.0*lTrack.BoostVector(); // get phi beta
	    TLorentzVector lTrackA_Boosted = lTrackA;
//	    lTrackA.Boost(beta_lTrack); // boost K+ back to phi rest frame
	    lTrackA_Boosted.Boost(beta_lTrack); // boost K+ back to phi rest frame
	    TVector3 momentumA_Boosted = lTrackA_Boosted.Vect(); // K+ momentum in phi rest frame

	    if(mStrangenessCut->passPhiEtaEast(lTrackA)) // K+ neg eta(east)
	    { // Below is West Only
	      TVector2 Q2Vector = Q2West[j];
	      // subtract auto-correlation from pos eta(west) event plane
	      if(flagB == 0 && mStrangenessCut->passTrackEP(lTrackB,dcaB) && mStrangenessCut->passTrackEtaWest(lTrackB,j)) // trackB
	      {
		Float_t  w = mStrangenessCorr->getWeight(lTrackB);
		TVector2 q2VectorB = mStrangenessCorr->calq2Vector(lTrackB);
		TVector2 q2CorrB   = mStrangenessCorr->getReCenterPar_West(0,cent9,runIndex,vz_sign,j,0); // 2nd
		Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
	      }
	      Float_t Res2 = mStrangenessCorr->getResolution2_EP(cent9,j);
	      Float_t Psi2_west = mStrangenessCorr->calShiftAngle2West_EP(Q2Vector,runIndex,cent9,vz_sign,j);
	      TVector3 nQ_West(TMath::Sin(Psi2_west),-1.0*TMath::Cos(Psi2_west),0.0); // normal vector of 2nd Event Plane
	      Double_t CosThetaStar = momentumA_Boosted.Dot(nQ_West)/(momentumA_Boosted.Mag()*nQ_West.Mag());

	      mStrangenessHistoManger->Fill(pt_lTrack,cent9,j,CosThetaStar,Res2,InvMass_lTrack,reweight,mX_flag,mMode);
	    }

	    if(mStrangenessCut->passPhiEtaWest(lTrackA)) // K+ pos eta (west)
	    { // Below is East Only
	      TVector2 Q2Vector = Q2East[j];
	      // subtract auto-correlation from pos eta(west) event plane
	      if(flagB == 0 && mStrangenessCut->passTrackEP(lTrackB,dcaB) && mStrangenessCut->passTrackEtaEast(lTrackB,j)) // trackB
	      {
		Float_t  w = mStrangenessCorr->getWeight(lTrackB);
		TVector2 q2VectorB = mStrangenessCorr->calq2Vector(lTrackB);
		TVector2 q2CorrB   = mStrangenessCorr->getReCenterPar_East(0,cent9,runIndex,vz_sign,j,0); // 2nd
		Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
	      }
	      Float_t Res2 = mStrangenessCorr->getResolution2_EP(cent9,j);
	      Float_t Psi2_east = mStrangenessCorr->calShiftAngle2East_EP(Q2Vector,runIndex,cent9,vz_sign,j);
	      TVector3 nQ_East(TMath::Sin(Psi2_east),-1.0*TMath::Cos(Psi2_east),0.0); // normal vector of 2nd Event Plane
	      Double_t CosThetaStar = momentumA_Boosted.Dot(nQ_East)/(momentumA_Boosted.Mag()*nQ_East.Mag());

	      mStrangenessHistoManger->Fill(pt_lTrack,cent9,j,CosThetaStar,Res2,InvMass_lTrack,reweight,mX_flag,mMode);
	    }
	  }
	}
      }
    }
  }

  cout << "." << flush;
  cout << " " << stop_event_use-start_event_use << "(" << 100 << "%)";
  cout << endl;
  if(mFlag_Embedding == 1)
  {
    out_map.close();
  }
}

//-------------------------------------------------------------------
void StStrangenessAna::Finish()
{
  mFile_OutPut->cd();
  mStrangenessHistoManger->Write(mX_flag,mMode);
  mFile_OutPut->Close();
}

#include <iostream>
#include <string>

void Load();

using namespace std;

// void run_StMcPhiMaker(const char* file = "/projecta/projectdirs/starprod/embedding/AuAu200_production_2011/Phi_114_20140915/P11id.SL11d_embed/2011/140/12140027/st_physics_adc_12140027_raw_4510001.event.root", std::string outFile = "test")
void run_StMcPhiMaker(const char* file = "/projecta/projectdirs/starprod/embedding/AuAu200_production_2011/Phi_108_20140915/P11id.SL11d_embed/2011/160/12160019/st_physics_adc_12160019_raw_4490002.event.root", std::string outFile = "test")
{
   //Check STAR Library. Please set SL_version to the original star library used
   // in the production from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl

   string SL_version = "SL11d";
   string env_SL = getenv("STAR");
   // cout << env_SL << endl;

   if (env_SL.find(SL_version) == string::npos)
   {
      cout << SL_version << " " << env_SL << endl;
      cout << "Environment Star Library does not match the requested library in run_st_etree.C. Exiting..." << endl;
      return;
   }

   // Load shared libraries
   Load();

   // Create chain
   StChain* chain = new StChain;

   // I/O maker
   StIOMaker* ioMaker = new StIOMaker;
   ioMaker->SetFile(file);
   ioMaker->SetIOMode("r");
   ioMaker->SetBranch("*", 0, "0");
   //ioMaker->SetBranch("McEventBranch",0,"r");
   ioMaker->SetBranch("geantBranch", 0, "r");
   ioMaker->SetBranch("eventBranch", 0, "r");

   TString mudstfile = file;

   if(mudstfile.First("$") != -1)
   {
     mudstfile.ReplaceAll("$","");
     mudstfile = getenv(mudstfile.Data());
   }

   mudstfile.ReplaceAll(".event.root", ".MuDst.root");
   mudstfile.ReplaceAll(".geant.root", ".MuDst.root");
   cout << "Reading MuDst file " << mudstfile << endl;
   StMuDstMaker* muDstMaker = new StMuDstMaker(0, 0, "", mudstfile.Data(), "", 100000, "MuDst");

   StMcEventMaker *mcEventMaker = new StMcEventMaker();
   mcEventMaker->doPrintEventInfo = false;
   mcEventMaker->doPrintMemoryInfo = false;

   StAssociationMaker* assoc = new StAssociationMaker;
   assoc->useInTracker();
   assoc->SetDebug();

   //.. see example in CVS: StRoot/macros/mudst/exampleEmc.C
   // Need St_db_Maker for Emc calibration
   // St_db_Maker* dbMk = new St_db_Maker("StarDb", "MySQL:StarDb");
   //dbMk->SetMaxEntryTime(20100301,0);
   //dbMk->SetDateTime(20080101,000001);
   int refMultCorrLoad = gSystem->Load("StRefMultCorr");
   StRefMultCorr* grefmultCorrUtil = NULL;

   if (refMultCorrLoad == -1)
   {
      cout << "StRefMultCorr library is not available" << endl;
   }
   else
   {
      grefmultCorrUtil = CentralityMaker::instance()->getRefMultCorr();
   }

   // Monte Carlo event maker
   StMcPhiMaker* analysis = new StMcPhiMaker;
   analysis->setOutFileName(outFile);
   analysis->setRefMultCorr(grefmultCorrUtil);

   // Initialize chain
   chain->Init();
   chain->EventLoop(1e6);
   // chain->EventLoop(4000);
   chain->Finish();

   //delete chain;
}

void Load()
{
   gROOT->Macro("loadMuDst.C");
   gROOT->Macro("LoadLogger.C");
   gSystem->Load("StMcEvent");
   gSystem->Load("StMcEventMaker");
   gSystem->Load("StAssociationMaker");
   gSystem->Load("StDbLib.so");
   gSystem->Load("StDbBroker.so");
   gSystem->Load("libglobal_Tables.so");
   gSystem->Load("St_db_Maker.so");
   gSystem->Load("StDetectorDbMaker");
   gSystem->Load("StTpcDb");
   gSystem->Load("StDbUtilities");
   gSystem->Load("StMcEvent");
   gSystem->Load("StMcEventMaker");
   gSystem->Load("StAssociationMaker");
   gSystem->Load("StEmcRawMaker");
   gSystem->Load("StEmcADCtoEMaker");
   gSystem->Load("StPreEclMaker");
   gSystem->Load("StEpcMaker");
   gSystem->Load("StEmcSimulatorMaker");
   gSystem->Load("StEmcUtil");
   gSystem->Load("StEEmcUtil");
   gSystem->Load("StEEmcDbMaker");
   gSystem->Load("StEmcTriggerMaker");
   gSystem->Load("StDaqLib");
   gSystem->Load("StMcPhiMaker");
}

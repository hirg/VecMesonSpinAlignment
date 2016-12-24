
#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;


StChain *chain;
void VecMesonTree(const Char_t *inputFile="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu19GeV/List/run_list/19GeV_140.list", const Int_t jobCounter = 140, const Int_t Mode = 0, const Int_t energy = 2, const Int_t flag_ME = 0)
{
  // mBeamEnergy[NumBeamEnergy] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  // Mode: 0 for re-center correction, 1 for shift correction, 3 for phi meson
  // flag_ME: 0 for Same Event, 1 for Mixed Event

  Int_t nEvents = 10000000;
  // Int_t nEvents = 5000;

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StRefMultCorr");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StAlexPhiMesonEvent");
  gSystem->Load("StVecMesonMaker");
  gSystem->Load("StRunIdEventsDb");

  chain = new StChain();

  StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst");

  StVecMesonMaker *VecMesonMaker = new StVecMesonMaker("VecMeson",picoMaker,jobCounter,Mode,energy,flag_ME);

  chain->Init();
  cout<<"chain->Init();"<<endl;
  int total = picoMaker->chain()->GetEntries();
  cout << " Total entries = " << total << endl;
  if(nEvents>total) nEvents = total;
  for (Int_t i=0; i<nEvents; i++)
  {
    if(i != 0 && i%50 == 0)
    {
      cout << "." << flush;
    }
    if(i != 0 && i%500 == 0)
    {
      Float_t event_percent = 100.0*i/nEvents;
      cout << " " << i << "(" << event_percent << "%)" << "\n" << " ==> Processing data " << flush;
    }

    chain->Clear();
    int iret = chain->Make(i);

    if (iret)
    { 
      cout << "Bad return code!" << iret << endl;
      break;
    }

    total++;

  }

  cout << "." << flush;
  cout << " " << nEvents << "(" << 100 << "%)";
  cout << endl;
  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!"<< endl;
  cout << "****************************************** " << endl;
  chain->Finish();
  cout << "****************************************** " << endl;
  cout << "total number of events  " << nEvents << endl;
  cout << "****************************************** " << endl;

  delete chain;
}


#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;


StChain *chain;
void ZdcSmdTree(const Char_t *inputFile="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu200GeV/List/run_list/200GeV_140.list", const Int_t jobCounter = 140, const Int_t Mode = 4, const Int_t energy = 6, const Int_t flag_ME = 0)
{
  // mBeamEnergy[NumBeamEnergy] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  // Mode: 
  //       0 = gain correction for East/West;
  //       1 = re-center correction for East/West;
  //       2 = shift correction for East/West;
  //       3 = shift correction for Full;
  //       4 = event plane resolution for East/West
  //       5 = charge particle v1 w.r.t. full ZDCSMD event plane (QA)
  //       6 = phi meson tree production w.r.t. full ZDCSMD event plane
  // flag_ME: 0 for Same Event, 1 for Mixed Event

  Int_t nEvents = 10000000;
  // Int_t nEvents = 50000;

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StRefMultCorr");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StAlexPhiMesonEvent");
  gSystem->Load("StZdcSmdMaker");
  gSystem->Load("StRunIdEventsDb");

  chain = new StChain();

  StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst");

  StZdcSmdMaker *ZdcSmdMaker = new StZdcSmdMaker("ZdcSmd",picoMaker,jobCounter,Mode,energy,flag_ME);

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

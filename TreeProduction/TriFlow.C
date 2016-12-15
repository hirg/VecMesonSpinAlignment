
#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;


StChain *chain;
void TriFlow(const Char_t *inputFile="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu200GeV/List/run_list/200GeV_511.list", const Int_t jobCounter = 511, const Int_t Mode = 2, const Int_t energy = 0, const Int_t flag_ME = 0, const Int_t flag_Embedding = 0)
// void TriFlow(const Char_t *inputFile="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu200GeV/List/Phi_list/embedding_List/run_list/200GeV_3.list", const Int_t jobCounter = 3, const Int_t Mode = 2, const Int_t energy = 0, const Int_t flag_ME = 0, const Int_t flag_Embedding = 1)
{
// Mode: 0 for re-center correction, 1 for shift correction, 2 for phi meson
// energy: 0 for 200 GeV, 1 for 39 GeV, 2 for 27 GeV (on progress)
// flag_ME: 0 for Same Event, 1 for Mixed Event
// flag_Embedding: 0 for production, 1 for embedding

        // Int_t nEvents = 10000000;
	Int_t nEvents = 5000;
	
        gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();

	gSystem->Load("StRefMultCorr");
	gSystem->Load("StPicoDstMaker");
	gSystem->Load("StAlexPhiMesonEvent");
        gSystem->Load("StTriFlowMaker");
	gSystem->Load("StRunIdEventsDb");

	chain = new StChain();

	StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst");

        StTriFlowMaker *TriFlowMaker = new StTriFlowMaker("TriFlow",picoMaker,jobCounter,Mode,energy,flag_ME,flag_Embedding);

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

#include <TSystem>

void FillSpinAlignmentZdcSmd(const Int_t energy = 6, const Int_t X_flag = 0, const Int_t List = 0, const Long64_t start_event = 0, const Long64_t stop_event = 1024, const Int_t mode = 0)
{
  // mBeamEnergy[NumBeamEnergy] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  // X_flag: 0 for Same Event, 1 for Mixed Event
  // List: different number for different TTree list
  // mode: 0 for phi meson, 1 for K*, 2 for K0S

  gSystem->Load("StRefMultCorr");
  gSystem->Load("StAlexPhiMesonEvent");
  gSystem->Load("StZdcSmdAna");
  gSystem->Load("StRunIdEventsDb");

  cout << "All libraries are loaded!!!!" << endl;

  cout << "Start to Read Trees!" << endl;

  StZdcSmdAna *mZdcSmdAna = new StZdcSmdAna(energy,X_flag,List,start_event,stop_event,mode);
  mZdcSmdAna->Init();
  mZdcSmdAna->Make();
  mZdcSmdAna->Finish();

//  cout << "End of the Calculation!!" << endl;
}

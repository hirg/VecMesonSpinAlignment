#include <TSystem>

// root4star -b -q Efficiency.C\(1,0,1,1024,0\)
void Efficiency(const Int_t Energy = 6, const Long64_t StartEvent = 0, const Long64_t StopEvent = 1000024, const Int_t PID = 0) // Energy = 0: 7GeV, 1: 11GeV, 2: 19GeV, 3: 27GeV, 4: 39GeV, 5: 62GeV, 6: 200 GeV | PID = 0: phi, 1 Kstar
{
  gSystem->Load("StEffTPC");

  cout << "All libraries are loaded!!!!" << endl;
  cout << "Start to Calculate Efficiency!!" << endl;

  StEffTPC *mEfficiency = new StEffTPC(Energy,StartEvent,StopEvent,PID);
  mEfficiency->Init();
  mEfficiency->Make();
  mEfficiency->Finish();
  delete mEfficiency;

  cout << "End of the Calculation!!" << endl;
}

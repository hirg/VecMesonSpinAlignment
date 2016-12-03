#include <TSystem>

// root4star -b -q EffMcPhi.C\(6,1,1024,0\)
void EffMcPhi(const int Energy = 6, const long StartEvent = 0, const long StopEvent = 10000024, const int PID = 0) // Energy = 0: 7GeV, 1: 11GeV, 2: 19GeV, 3: 27GeV, 4: 39GeV, 5: 62GeV, 6: 200 GeV | PID = 0: phi, 1 Kstar
{
  gSystem->Load("../.sl64_gcc482/lib/libStEffMcPhi.so");

  cout << "All libraries are loaded!!!!" << endl;
  cout << "Start to Calculate Efficiency!!" << endl;

  StEffMcPhi *mEffMcPhi = new StEffMcPhi(Energy,StartEvent,StopEvent,PID);
  mEffMcPhi->Init();
  mEffMcPhi->Make();
  mEffMcPhi->Finish();
  delete mEffMcPhi;

  cout << "End of the Calculation!!" << endl;
}

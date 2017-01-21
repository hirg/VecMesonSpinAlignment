#include <TSystem>

// root4star -b -q EffMcPhi.C\(2,0,100000024,0,0,0\)
void EffMcPhi(const int Energy = 6, const long StartEvent = 0, const long StopEvent = 100000024, const int PID = 0, const int Year = 0, const int Cut = 0) // Energy = 0: 7GeV, 1: 11GeV, 2: 19GeV, 3: 27GeV, 4: 39GeV, 5: 62GeV, 6: 200 GeV | PID = 0: phi, 1 Kstar | Year = 0: run11, 1 run10 | Cut = 0: pr, 1: gl
{
  gSystem->Load("../.sl64_gcc482/lib/libStEffMcPhi.so");

  cout << "All libraries are loaded!!!!" << endl;
  cout << "Start to Calculate Efficiency!!" << endl;

  StEffMcPhi *mEffMcPhi = new StEffMcPhi(Energy,StartEvent,StopEvent,PID,Year,Cut);
  mEffMcPhi->Init();
  mEffMcPhi->Make();
  mEffMcPhi->Finish();
  delete mEffMcPhi;

  cout << "End of the Calculation!!" << endl;
}

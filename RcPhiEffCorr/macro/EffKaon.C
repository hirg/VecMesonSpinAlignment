#include <TSystem>

// root4star -b -q Efficiency.C\(1,0,1,1024,0\)
void EffKaon(const int Energy = 4, const long StartEvent = 0, const long StopEvent = 10000024, const int PID = 1, const int Year = 1, const int Cut = 0) // Energy = 0: 7GeV, 1: 11GeV, 2: 19GeV, 3: 27GeV, 4: 39GeV, 5: 62GeV, 6: 200 GeV | PID = 0: K+, 1 K- | Year = 0: run11, 1 run10 | Cut = 0: pr, 1: gl
{
  gSystem->Load("../.sl64_gcc482/lib/libStEffKaon.so");

  cout << "All libraries are loaded!!!!" << endl;
  cout << "Start to Calculate Efficiency!!" << endl;

  StEffKaon *mEffKaon = new StEffKaon(Energy,StartEvent,StopEvent,PID,Year,Cut);
  mEffKaon->Init();
  mEffKaon->Make();
  mEffKaon->Finish();
  delete mEffKaon;

  cout << "End of the Calculation!!" << endl;
}

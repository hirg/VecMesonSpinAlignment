#include <string>
#include <map>
using namespace std;

void testStRunIdEventsDb(Int_t energy, Int_t year) {

  map<string, Int_t> test_runids;
  test_runids["7_2010"]   = 11122112;
  test_runids["11_2010"]  = 11151074;
  test_runids["19_2011"]  = 12119049;
  test_runids["27_2011"]  = 12172050;
  test_runids["39_2010"]  = 11108040;
  test_runids["62_2010"]  = 11092087;
  test_runids["200_2010"] = 11021031;
  test_runids["200_2011"] = 12158041;

  gSystem->Load("src/StRunIdEventsDb/StRunIdEventsDb_cxx.so");
  StRunIdEventsDb* runevdb = new StRunIdEventsDb(energy,year);
  runevdb->readDb();

  cout << runevdb->getTotalNrRunIds() << endl;
  TString strEnYr = Form("%d_%d",energy,year);
  Int_t runid = test_runids[strEnYr.Data()];
  cout << "runid = " << runid << endl;
  Int_t index = runevdb->getRunIdIndex(runid); 
  cout << index << "  :  " << runevdb->getNrEvents("total",runid) << endl;
  runevdb->printAllTrigNames();
  string tn = runevdb->getTrigName(0);
  cout << index << "  :  " << runevdb->getNrEvents(tn,runid) << endl;
  Int_t indexNew = 3;
  Int_t runidNew = runevdb->getRunId(indexNew);
  cout << indexNew << "  :  " << runevdb->getNrEvents(tn,runidNew) << endl;

}

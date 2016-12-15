#include <Riostream.h>
#include "StRunIdEventsDb.h"
#include <TSystem.h>
#include <limits>
#include <cstdlib>
using namespace std;

ClassImp(StRunIdEventsDb)

StRunIdEventsDb* StRunIdEventsDb::mRunIdEventsDb = NULL;
map<Int_t,Float_t> StRunIdEventsDb::mEnergyMap = StRunIdEventsDb::createEnergyMap();

StRunIdEventsDb::StRunIdEventsDb(Float_t energy, Int_t year) {

  mNrRunIds = 0;
  SetEnergyYear(energy,year);

  // define input filenames
  TString srcdir = "StRoot/";
  TString indir = "StRunIdEventsDb/database/";
  TString instr = Form("dbfiles/database_%.1f_%d",mEnergy,mYear);
  mDbFile = srcdir + indir + instr;
  mCfgFile = srcdir + indir + "config";

  readDb();
}

StRunIdEventsDb* StRunIdEventsDb::Instance(Float_t energy, Int_t year, Bool_t kForce) {
  if ( !mRunIdEventsDb || kForce ) {
    mRunIdEventsDb = new StRunIdEventsDb(energy,year);
  }
  return mRunIdEventsDb;
}

StRunIdEventsDb* StRunIdEventsDb::Instance(Int_t energy, Int_t year, Bool_t kForce) {
  if ( !mRunIdEventsDb || kForce ) {
    mRunIdEventsDb = new StRunIdEventsDb(mEnergyMap[energy],year);
  }
  return mRunIdEventsDb;
}

void StRunIdEventsDb::readConfigFile() {

  ifstream in; 
  in.open(mCfgFile.Data());
  if ( in.fail() ) { exitProg(mCfgFile); }

  Int_t year, ncols;
  Float_t energy;
  while ( in >> energy >> year >> ncols ) {
    if ( mTrigNames[0].empty() ) {
      if ( energy == mEnergy && year == mYear ) {
	mNTrigs = ncols;
	for ( Int_t n = 0; n < mNTrigs; ++n ) {
	  in >> mTrigNames[n];
	  cout << "  " << mTrigNames[n] << endl;
	}
      }
      else {
	in.ignore(numeric_limits<streamsize>::max(),'\n');
      }
    }
    else break;
  }

  in.close();
  cout << "StRunIdEventsDb(): " << mNTrigs << " triggers loaded." << endl;
}

void StRunIdEventsDb::readDb() { 

  // read the config file first
  readConfigFile();

  // read the database file second
  ifstream in; 
  in.open(mDbFile.Data());
  if ( in.fail() ) { exitProg(mDbFile); }

  Int_t runid;
  while ( in >> runid ) {
    for ( Int_t n = 0; n < mNTrigs; ++n ) {
      string key = getKey(mTrigNames[n],runid);
      mMapRunIdsEvts[key].index = mNrRunIds;
      in >> mMapRunIdsEvts[key].nevts;
      if ( mNrRunIds < 5 ) {
	cout << key << "  " << mMapRunIdsEvts[key].index << "  " << mMapRunIdsEvts[key].nevts << endl;
      }
    }
    mNrRunIds++;
  }

  in.close();
  printLoadMessage();
}

string StRunIdEventsDb::getKey(string tn, Int_t runid) {
  return Form("%s_%08d",tn.c_str(),runid); 
}

Int_t StRunIdEventsDb::getRunIdFromKey(string key) {
  size_t found = key.find_last_of("_");
  return atoi(key.substr(found+1).c_str());
}

Int_t StRunIdEventsDb::getRunId(Int_t ind) {
  if ( ind > mNrRunIds-1 ) {
    printf("StRunIdEventsDb(): index for mNrRunIds out of bounds (max = %d)!",mNrRunIds-1);
    return -1;
  }
  Int_t runid = -1;
  map<string,data>::iterator it;
  for ( it = mMapRunIdsEvts.begin(); it != mMapRunIdsEvts.end(); ++it ) {
    data dt = (*it).second;
    if ( dt.index == ind ) {
      string key = (*it).first;
      runid = getRunIdFromKey(key);
      break;
    }
  }
  return runid;
}

Int_t StRunIdEventsDb::getRunIdIndex(Int_t runid) {
  string key = getKey(mTrigNames[0],runid);
  if ( mMapRunIdsEvts.find(key) == mMapRunIdsEvts.end() ) {
    printMissingKeyMessage(key); 
    return -1;
  }
  return mMapRunIdsEvts[key].index;
}

Int_t StRunIdEventsDb::getNrEvents(string tn, Int_t runid) {
  string key = getKey(tn,runid);
  if ( ! mMapRunIdsEvts.count(key) ) {
    printMissingKeyMessage(key); 
    return -1;
  }
  return mMapRunIdsEvts[key].nevts;
}

string StRunIdEventsDb::getTrigName(Int_t n) {
  if ( n > mNTrigs-1 ) {
    printf("StRunIdEventsDb(): index for mTrigNames out of bounds (max = %d)!",mNTrigs-1);
    return "";
  }
  return mTrigNames[n];
}

void StRunIdEventsDb::printAllTrigNames() {
  for ( Int_t n = 0; n < mNTrigs; ++n ) {
    cout << "n = " << n << " : " << mTrigNames[n].c_str() << endl;
  }
}

void StRunIdEventsDb::printMissingKeyMessage(string key) {
  printf("StRunIdEventsDb(): key %s does not exist! returning -1!\n",key.c_str());
}

void StRunIdEventsDb::exitProg(TString fstr) {
  printf("StRunIdEventsDb(): failed to open file: \n   %s!\n",fstr.Data());
  gSystem->Exit(0); 
}

void StRunIdEventsDb::printLoadMessage() {
  printf("--> loaded StRunIdEventsDb's \n %s \n with %d RunIds\n",mDbFile.Data(),mNrRunIds);
}


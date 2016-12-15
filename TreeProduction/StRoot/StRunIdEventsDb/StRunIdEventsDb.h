#ifndef StRunIdEventsDb_h
#define StRunIdEventsDb_h

#include <TROOT.h>
#include <string>
#include <map>
#include <TString.h>
using namespace std;

const Int_t MaxNrTrigNames = 10;

class StRunIdEventsDb {

 private: 
   StRunIdEventsDb(Float_t energy, Int_t year);
   static StRunIdEventsDb* mRunIdEventsDb;

   Int_t mNrRunIds;
   Float_t mEnergy;
   Int_t mYear;
   Int_t mNTrigs;
   
   string mTrigNames[MaxNrTrigNames];
   TString mDbFile;
   TString mCfgFile;

   struct data {
     Int_t index;
     Int_t nevts;
   };
   map<string,data> mMapRunIdsEvts;

   void SetEnergyYear(Float_t e, Int_t y) { mEnergy = e; mYear = y; }
   void readConfigFile();
   string getKey(string tn, Int_t runid);
   Int_t getRunIdFromKey(string key);
   void printLoadMessage();
   void exitProg(TString fstr);
   void printMissingKeyMessage(string key);
   void  readDb();
   static map<Int_t,Float_t> createEnergyMap() {
     map<Int_t,Float_t> m;
     m[7] = 7.7;
     m[11] = 11.5;
     m[19] = 19.6;
     m[27] = 27.0;
     m[39] = 39.0;
     m[62] = 62.4;
     m[200] = 200.0;
     return m;
   }

 public:
   static StRunIdEventsDb* Instance(Float_t energy, Int_t year, Bool_t kForce = 0);
   static StRunIdEventsDb* Instance(Int_t energy, Int_t year, Bool_t kForce = 0);
   virtual ~StRunIdEventsDb() {};

   static map<Int_t,Float_t> mEnergyMap;
   Int_t getTotalNrRunIds() { return mNrRunIds; }
   Int_t getRunIdIndex(Int_t runid);
   Int_t getRunId(Int_t ind);
   Int_t getNrEvents(string tn, Int_t runid);
   void printAllTrigNames();
   string getTrigName(Int_t n);

  ClassDef(StRunIdEventsDb,0)
};
#endif


#ifndef StZdcSmdMEKey_h
#define StZdcSmdMEKey_h
#include "StPhysicalHelixD.hh"
#include "Rtypes.h"
#include <map>
#include <vector>
#include "TLorentzVector.h"
#include "StThreeVectorF.hh"

class MEKey
{
  public:
    int Centrality;
    int vertex;
    int EventPlane;
    int Event_depth;
    int charge;

    MEKey() {}
    MEKey(int i, int j, int k, int l, int m) : Centrality(i), vertex(j), EventPlane(k), Event_depth(l), charge(m) {}
    virtual ~MEKey() {/* */};

    bool operator<(const MEKey &r) const
    { // "r == right"
      if(Centrality == r.Centrality)
      {
        if(vertex == r.vertex)
	{
	  if(EventPlane == r.EventPlane)
	  {
	    if(Event_depth == r.Event_depth)
	    {
	      return charge < r.charge;
	    }
	    else
	    {
	      return Event_depth < r.Event_depth;
	    }
	  }
	  else
	  {
	    return EventPlane < r.EventPlane;
	  }
	}
	else
	{
	  return vertex < r.vertex;
	}
      }
      else
      {
        return Centrality < r.Centrality;
      }
    }

    void print()
    {
      cout << "Centrality = " << MEKey::Centrality << endl;
      cout << "vertex = " << MEKey::vertex << endl;
      cout << "EventPlane = " << MEKey::EventPlane << endl;
      cout << "Event_depth = " << MEKey::Event_depth << endl;
      cout << "charge = " << MEKey::charge << endl;
    }

    ClassDef(MEKey, 1)
};

typedef std::map<MEKey, std::vector<int> > vectorIntMap;
typedef std::map<MEKey, std::vector<float> > vectorFloatMap;
typedef std::map<MEKey, std::vector<StPhysicalHelixD> > vectorHelixMap;
typedef std::map<MEKey, std::vector<TLorentzVector> > vectorLorentzMap;
typedef std::map<MEKey, std::vector<StThreeVectorF> > vectorStThreeFMap;

#endif

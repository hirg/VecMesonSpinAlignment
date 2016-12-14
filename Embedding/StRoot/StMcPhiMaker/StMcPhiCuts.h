#ifndef StMcPhiCuts_H
#define StMcPhiCuts_H

/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */
#include <vector>

#include "Rtypes.h"
#include "StEvent/StEnumerations.h"

namespace McPhiCuts
{
  std::vector<unsigned int> getAllTriggers()
  {
    std::vector<unsigned int> t;
    t.push_back(350003); // minBias triger
    t.push_back(350013); 
    t.push_back(350023); 
    t.push_back(350033); 
    t.push_back(350043);

    return t;
  }

  std::vector<unsigned int> const interesting_triggers = getAllTriggers();

  float const mcTrackStartVtxR = 1.0; // maximum
  int const geantId_phi = 10151;
  int const geantId_KP  = 11;
  int const geantId_KM  = 12;

  StDedxMethod dedxMethod = kLikelihoodFitId;

  int const maxNumberOfTriggers = 6;
}
#endif

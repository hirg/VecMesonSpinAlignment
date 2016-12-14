#ifndef StMcAnaCuts_H
#define StMcAnaCuts_H

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

namespace McAnaCuts
{
  std::vector<unsigned int> getAllTriggers()
  {
    std::vector<unsigned int> t;
    // t.push_back(350003); // miniBias triger @ 200 GeV
    // t.push_back(350013); 
    // t.push_back(350023); 
    // t.push_back(350033); 
    // t.push_back(350043);

    // t.push_back(270001); // miniBias triger @ 62.4 GeV
    // t.push_back(270011); 
    // t.push_back(270021); 
    // t.push_back(270005); 

    t.push_back(280001); // miniBias triger @ 39 GeV
    t.push_back(280002); 

    // t.push_back(360001); // miniBias triger @ 27 GeV
    // t.push_back(360002); 

    // t.push_back(340001); // miniBias triger @ 19.6 GeV
    // t.push_back(340011); 
    // t.push_back(340021); 
    // t.push_back(340002); 
    // t.push_back(340012);
    // t.push_back(340022);

    return t;
  }

  std::vector<unsigned int> const interesting_triggers = getAllTriggers();

  float const mcTrackStartVtxR = 1.0; // maximum
  int const geantId = 11; // K+
  // int const geantId = 12; // K-

  StDedxMethod dedxMethod = kLikelihoodFitId;

  int const maxNumberOfTriggers = 6;
}
#endif

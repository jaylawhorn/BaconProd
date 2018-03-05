#ifndef BACONPROD_UTILS_TRIGGERTOOLS_HH
#define BACONPROD_UTILS_TRIGGERTOOLS_HH

#include "FWCore/Framework/interface/Event.h"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/Utils/interface/TriggerRecord.hh"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include <vector>

namespace baconhep {

class TriggerTools
{
public:
  static TriggerObjects  matchHLT(const double eta, const double phi, 
				  const std::vector<TriggerRecord> &triggerRecords,
				  const trigger::TriggerEvent &triggerEvent);

  static TriggerObjects matchHLT(const double eta, const double phi, 
				 const edm::Event &iEvent,
				 const std::vector<TriggerRecord> &triggerRecords,
				 const edm::TriggerNames &triggerNames,
				 const edm::TriggerResults &triggerResults,
				 pat::TriggerObjectStandAloneCollection &triggerObjects);
  
};
  
}
#endif

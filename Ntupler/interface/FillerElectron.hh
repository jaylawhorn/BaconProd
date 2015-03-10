#ifndef BACONPROD_NTUPLER_FILLERELECTRON_HH
#define BACONPROD_NTUPLER_FILLERELECTRON_HH

#include "BaconProd/Utils/interface/TriggerTools.hh"
#include <vector>
#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
class TClonesArray;
namespace trigger {
  class TriggerEvent;
}


namespace baconhep
{
  class FillerElectron
  {
    public:
      FillerElectron(const edm::ParameterSet &iConfig);
      ~FillerElectron();
      
      void fill(TClonesArray                                *array,           // output array to be filled
                const edm::Event                            &iEvent,          // event info
		const edm::EventSetup                       &iSetup,          // event setup info
		const reco::Vertex                          &pv,              // event primary vertex
		const std::vector<TriggerRecord>            &triggerRecords,  // list of trigger names and objects
		const trigger::TriggerEvent                 &triggerEvent);   // event trigger objects
  
    protected:
      double dEtaInSeed(const reco::GsfElectron& ele);
      
      
      // Electron cuts
      double fMinPt;
      
      // EDM object collection names
      std::string fEleName;
      std::string fPFCandName;
      std::string fTrackName;
      std::string fConvName;
      std::string fSCName;
  };
}
#endif

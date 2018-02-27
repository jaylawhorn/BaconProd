#ifndef BACONPROD_NTUPLER_FILLERELECTRON_HH
#define BACONPROD_NTUPLER_FILLERELECTRON_HH

#include "BaconProd/Utils/interface/TriggerTools.hh"
#include <vector>
#include <string>

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

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
    FillerElectron(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC);
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
      edm::EDGetTokenT<reco::GsfElectronCollection> fEleName_token;
      std::string fPFCandName;
      edm::EDGetTokenT<reco::PFCandidateCollection> fPFCandName_token;
      std::string fTrackName;
      edm::EDGetTokenT<reco::TrackCollection> fTrackName_token;
      std::string fBeamspotName;
      edm::EDGetTokenT<reco::BeamSpot> fBeamspotName_token;
      std::string fConvName;
      edm::EDGetTokenT<reco::ConversionCollection> fConvName_token;
      std::string fSCName;
      edm::EDGetTokenT<reco::SuperClusterCollection> fSCName_token;
  };
}
#endif

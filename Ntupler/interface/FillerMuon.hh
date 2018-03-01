#ifndef BACONPROD_NTUPLER_FILLERMUON_HH
#define BACONPROD_NTUPLER_FILLERMUON_HH

#include "BaconProd/Utils/interface/TriggerTools.hh"
#include <vector>
#include <string>

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class TClonesArray;
namespace trigger {
  class TriggerEvent;
}


namespace baconhep
{
  class FillerMuon
  {
  public:
    FillerMuon(const edm::ParameterSet &iConfig,const bool useAOD,edm::ConsumesCollector && iC);
    ~FillerMuon();
    
    // AOD
    void fill(TClonesArray			 *array,	   // output array to be filled
	      const edm::Event		 &iEvent,	   // event info
	      const edm::EventSetup		 &iSetup,	   // event setup info
	      const reco::Vertex		 &pv,	           // event primary vertex
	      const std::vector<TriggerRecord> &triggerRecords,  // list of trigger names and objects
	      const trigger::TriggerEvent	 &triggerEvent);   // event trigger objects
    
    // miniAOD ===
    void fill(TClonesArray                                 *array,            // output array to be filled
	      const edm::Event                             &iEvent,           // event info
	      const edm::EventSetup                        &iSetup,           // event setup info
	      const reco::Vertex                           &pv,               // event primary vertex
	      const std::vector<TriggerRecord>             &triggerRecords,   // list of trigger names and objects
	      const pat::TriggerObjectStandAloneCollection &triggerObjects); // event trigger objects
    
  protected:

    void computeIso(double &iEta,double &iPhi, const double extRadius,
		    const reco::PFCandidateCollection    &puppi,
		    float &out_chHadIso, float &out_gammaIso, float &out_neuHadIso) const;
    
    void computeIso(double &iEta,double &iPhi, const double extRadius,
		    const pat::PackedCandidateCollection    &puppi,
		    float &out_chHadIso, float &out_gammaIso, float &out_neuHadIso) const;
    
    // Muon cuts
    double fMinPt;
    
    // EDM object collection names
    std::string fMuonName;
    std::string fPFCandName;
    std::string fTrackName;

    // general tracks cuts
    bool   fSaveTracks;
    double fTrackMinPt;

    bool fUseAOD;

    edm::EDGetTokenT<reco::MuonCollection> fMuonName_token;
    edm::EDGetTokenT<reco::PFCandidateCollection> fPFCandName_token;
    edm::EDGetTokenT<reco::TrackCollection> fTrackName_token;

    edm::EDGetTokenT<pat::MuonCollection> fMuonPATName_token;
    
  };
}
#endif

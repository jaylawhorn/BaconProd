#ifndef BACONPROD_NTUPLER_FILLERELECTRON_HH
#define BACONPROD_NTUPLER_FILLERELECTRON_HH

#include "BaconProd/Utils/interface/TriggerTools.hh"
#include <vector>
#include <string>

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

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
namespace pat {
  class Electron;
}

namespace baconhep
{
  class FillerElectron
  {
  public:
    FillerElectron(const edm::ParameterSet &iConfig, const bool useAOD, edm::ConsumesCollector && iC);
    ~FillerElectron();
    
    // AOD
    void fill(TClonesArray                                *array,           // output array to be filled
	      const edm::Event                            &iEvent,          // event info
	      const edm::EventSetup                       &iSetup,          // event setup info
	      const reco::Vertex                          &pv,              // event primary vertex
	      const std::vector<TriggerRecord>            &triggerRecords,  // list of trigger names and objects
	      const trigger::TriggerEvent                 &triggerEvent);   // event trigger objects

    //miniAOD
    void fill(TClonesArray                                *array,           // output array to be filled
	      const edm::Event                            &iEvent,          // event info
	      const edm::EventSetup                       &iSetup,          // event setup info
	      const reco::Vertex                          &pv,              // event primary vertex
	      //const trigger::Results &triggerResults,
	      const std::vector<TriggerRecord>            &triggerRecords,  // list of trigger names and objects
	      const pat::TriggerObjectStandAloneCollection &triggerObjects);   // event trigger objects
    
  protected:
    double dEtaInSeed(const reco::GsfElectron& ele);
    double dEtaInSeed(const pat::Electron& ele); 
    void computeIso(double &iEta,double &iPhi, const double extRadius,
		    const reco::PFCandidateCollection    &puppi,
		    float &out_chHadIso, float &out_gammaIso, float &out_neuHadIso) const;    

    void computeIso(double &iEta,double &iPhi, const double extRadius,
		    const pat::PackedCandidateCollection    &puppi,
		    float &out_chHadIso, float &out_gammaIso, float &out_neuHadIso) const;
      
    // Electron cuts
    double fMinPt;
    
    // EDM object collection names
    std::string fEleName;
    std::string fBeamspotName;
    std::string fPFCandName;
    std::string fTrackName;
    std::string fConvName;

    edm::InputTag fSCName;
    // PF cluster isolation info (not in AOD)
    edm::InputTag fEcalPFClusterIsoMapTag;
    edm::InputTag fHcalPFClusterIsoMapTag;

    bool fUseAOD;
    edm::EDGetTokenT<reco::GsfElectronCollection> fEleName_token;
    edm::EDGetTokenT<pat::ElectronCollection> fElePATName_token;
    edm::EDGetTokenT<reco::BeamSpot> fBeamspotName_token;
    edm::EDGetTokenT<reco::PFCandidateCollection> fPFCandName_token;
    edm::EDGetTokenT<reco::TrackCollection> fTrackName_token;
    edm::EDGetTokenT<reco::ConversionCollection> fConvName_token;
    edm::EDGetTokenT<reco::SuperClusterCollection> fSCName_token;
    edm::EDGetTokenT<edm::ValueMap<float> > fEcalPFClusterIsoMap_token;
    edm::EDGetTokenT<edm::ValueMap<float> > fHcalPFClusterIsoMap_token;
  };
}
#endif

#ifndef BACONPROD_NTUPLER_FILLEREVENTINFO_HH
#define BACONPROD_NTUPLER_FILLEREVENTINFO_HH

#include <string>
#include "FWCore/Framework/interface/Frameworkfwd.h"     // declaration of EDM types
#include "FWCore/Framework/interface/EDAnalyzer.h"       // EDAnalyzer class
#include "FWCore/ParameterSet/interface/ParameterSet.h"  // Parameters
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/Framework/interface/MakerMacros.h"      // definitions for declaring plug-in modules
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
namespace trigger {
  class TriggerEvent;
}


namespace baconhep
{
  class TEventInfo;  // foward declaration
  class FillerEventInfo
  {
    public:
    FillerEventInfo(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC);
      ~FillerEventInfo();
      
      void fill(TEventInfo         *evtInfo,       // output object to be filled
                const edm::Event   &iEvent,        // EDM event info
		const reco::Vertex &pv,            // event primary vertex
		const bool          hasGoodPV,     // flag for if PV passing cuts is found
		const TriggerBits   triggerBits);//,   // bits for corresponding fired triggers
	       
    protected:
      void computeTrackMET(const reco::Vertex &pv, 
                           const reco::PFCandidateCollection *pfCandCol,
                           float &out_met, float &out_metphi);
    
    
      // EDM object collection names
      std::string fPFCandName;
      std::string fPUInfoName;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > fPUInfoName_token;
      std::string fBSName;
      edm::EDGetTokenT<reco::BeamSpot> fBSName_token;
      std::string fPFMETName;
      edm::EDGetTokenT<reco::PFMETCollection> fPFMETName_token;
      std::string fPFMETCName;
      edm::EDGetTokenT<reco::PFMETCollection> fPFMETCName_token;
      std::string fPuppETName;
      edm::EDGetTokenT<reco::PFMETCollection> fPuppETName_token;
      std::string fMVAMETName;
      std::string fCHMETName;
      edm::EDGetTokenT<reco::PFMETCollection> fCHMETName_token;
//      std::string fMVAMET0Name;
      std::string fRhoIsoName;
      edm::EDGetTokenT<double> rhoIsoTag_token;
      std::string fRhoJetName;
      edm::EDGetTokenT<double> rhoJetTag_token;
      bool        fFillMET;
      bool        fFillMETFilters;
//      bool        fAddSusyGen;
  };
}
#endif

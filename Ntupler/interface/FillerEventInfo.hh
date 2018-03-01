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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/MET.h"


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
//  class TSusyGen;    // ditto
  class FillerEventInfo
  {
  public:
    FillerEventInfo(const edm::ParameterSet &iConfig,const bool useAOD, edm::ConsumesCollector && iC);
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

    void computeTrackMET(const pat::PackedCandidateCollection *pfCandCol,
			 float &out_met, float &out_metphi);
    
    // EDM object collection names
    std::string fPFCandName;
    std::string fPUInfoName;
    std::string fPVName;
    std::string fBSName;

    std::string fMETName;
    std::string fPFMETName;
    std::string fPFMETCName;
    std::string fPuppETName;
    std::string fPuppETCName;

    std::string fRhoIsoName;
    std::string fRhoJetName;
    
    bool        fFillMET;
    bool        fFillMETFilters;
    bool        fUseAOD;

    edm::EDGetTokenT<reco::PFMETCollection> fPFMETName_token;
    edm::EDGetTokenT<reco::PFMETCollection> fPFMETCName_token;
    edm::EDGetTokenT<reco::PFMETCollection> fPuppETName_token;
    edm::EDGetTokenT<reco::PFMETCollection> fPuppETCName_token;

    edm::EDGetTokenT<pat::METCollection> fPFMETPATName_token;
    //edm::EDGetTokenT<pat::METCollection> fPFMETCPATName_token;
    edm::EDGetTokenT<pat::METCollection> fPuppETPATName_token;
    //edm::EDGetTokenT<pat::METCollection> fPuppETCPATName_token;

    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > fPUInfoName_token;
    edm::EDGetTokenT<reco::BeamSpot> fBSName_token;

    edm::EDGetTokenT<reco::PFCandidateCollection> fPFCandName_token;
    edm::EDGetTokenT<pat::PackedCandidateCollection> fPackCandName_token;

    edm::EDGetTokenT<double> rhoIsoTag_token;
    edm::EDGetTokenT<double> rhoJetTag_token;
    
  };
}
#endif

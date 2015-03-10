#ifndef BACONPROD_NTUPLER_FILLERPHOTON_HH
#define BACONPROD_NTUPLER_FILLERPHOTON_HH

#include "BaconProd/Utils/interface/TriggerTools.hh"
//#include "BaconProd/Utils/interface/PhotonMVACalculator.hh"
#include <vector>
#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class TClonesArray;
namespace trigger {
  class TriggerEvent;
}


namespace baconhep
{
  class FillerPhoton
  {
    public:
      FillerPhoton(const edm::ParameterSet &iConfig);
      ~FillerPhoton();
      
      void fill(TClonesArray			            *array,	      // output array to be filled
                const edm::Event		            &iEvent,	      // event info
	        const edm::EventSetup		            &iSetup,	      // event setup info
                const reco::Vertex                          &pv,              // event primary vertex
	        const std::vector<TriggerRecord>            &triggerRecords,  // list of trigger names and objects
	        const trigger::TriggerEvent	            &triggerEvent);   // event trigger objects
/*            
      void computeVtxIso    (const reco::Photon &photon,
			     const std::vector<reco::PFCandidate>   &pf,
			     const std::vector<reco::Vertex>        &iVetex,
			     float &out_chHadIsoWvtx,float &out_chHadIsoFirstVtx) const;
      std::vector<float>  getESHits(double X, double Y, double Z, std::map<DetId, EcalRecHit> rechits_map, const CaloGeometry& geometry, CaloSubdetectorTopology *topology_p, int row);
      std::vector<float>  getESShape(std::vector<float> ESHits0);      
*/

      // Photon cuts
      double fMinPt;
      
      // EDM object collection names
      std::string fPhotonName;
      std::string fPFCandName;
      std::string fEleName;
      std::string fConvName;
      std::string fSCName;
//      std::string fEBRecHitName;
//      std::string fEERecHitName;
      std::string fPVName;  
  
//      PhotonMVACalculator *fPhotonMVA; 
  };
}
#endif

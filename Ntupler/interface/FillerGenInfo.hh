#ifndef BACONPROD_NTUPLER_FILLERGENINFO_HH
#define BACONPROD_NTUPLER_FILLERGENINFO_HH

#include <string>
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
class TClonesArray;


namespace baconhep
{
  class TGenEventInfo;  // foward declaration

  class FillerGenInfo
  {
    public:
    FillerGenInfo(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC);
      ~FillerGenInfo();
      
      void fill(TGenEventInfo    *genEvtInfo,     // output object to be filled
                TClonesArray     *particlesArr,   // output array of particles to be filled
		const edm::Event &iEvent);        // EDM event info
  
    protected:  
      
      // EDM object collection names
      std::string fGenEvtInfoName;
      edm::EDGetTokenT<GenEventInfoProduct> fGenEvtInfoName_token;
      std::string fLHEEventName;
      edm::EDGetTokenT<LHEEventProduct> fLHEEventName_token;
      std::string fGenParName;
      edm::EDGetTokenT<reco::GenParticleCollection> fGenParName_token;
      bool        fFillAll;
  };
}
#endif

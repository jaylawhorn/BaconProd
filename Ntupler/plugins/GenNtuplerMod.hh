#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"      // definitions for declaring plug-in modules
#include "FWCore/Framework/interface/Frameworkfwd.h"     // declaration of EDM types
#include "FWCore/Framework/interface/EDAnalyzer.h"       // EDAnalyzer class
#include "FWCore/ParameterSet/interface/ParameterSet.h"  // Parameters
#include "FWCore/Utilities/interface/InputTag.h"
#include <string>                                        // string class


// forward class declarations
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


class TFile;
class TH1D;
class TTree;
class TClonesArray;

namespace baconhep {
  class TGenEventInfo;
  class FillerGenInfo;
}

//
class GenNtuplerMod : public edm::EDAnalyzer {
  public:
    explicit GenNtuplerMod(const edm::ParameterSet &iConfig);
    ~GenNtuplerMod();

  private:
    virtual void beginJob();
    virtual void endJob();
    virtual void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup);

    virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup);
    virtual void endRun  (const edm::Run &iRun, const edm::EventSetup &iSetup);
    virtual void beginLuminosityBlock(const edm::LuminosityBlock &iLumi, const edm::EventSetup &iSetup);
    virtual void endLuminosityBlock  (const edm::LuminosityBlock &iLumi, const edm::EventSetup &iSetup);
  
    // bacon fillers
    baconhep::FillerGenInfo   *fFillerGenInfo;
   
    bool fIsActiveGenInfo;   
    // Objects and arrays for output file
    std::string              fOutputName;
    TFile                   *fOutputFile;
    TH1D                    *fTotalEvents;
    TTree                   *fEventTree;
    baconhep::TGenEventInfo *fGenEvtInfo;
    TClonesArray            *fGenParArr;
  //TClonesArray            *fGenJetArr;
};

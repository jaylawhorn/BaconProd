#include "FWCore/Framework/interface/MakerMacros.h"      // definitions for declaring plug-in modules
#include "FWCore/Framework/interface/Frameworkfwd.h"     // declaration of EDM types
#include "FWCore/Framework/interface/EDAnalyzer.h"       // EDAnalyzer class
#include "FWCore/ParameterSet/interface/ParameterSet.h"  // Parameters
#include "FWCore/Utilities/interface/InputTag.h"
#include <string>                                        // string class

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

// forward class declarations
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
class TFile;
class TH1D;
class TTree;
class TClonesArray;
namespace edm {
  class TriggerResults;
  class TriggerNames;
}
namespace baconhep {
  class TEventInfo;
  class TGenEventInfo;
  class TTrigger;
  class FillerEventInfo;
  class FillerGenInfo;
  class FillerVertex;
  class FillerElectron;
  class FillerMuon;
  class FillerPhoton;
  class FillerPF;
  class FillerJet;
}

//
class NtuplerMod : public edm::EDAnalyzer {
  public:
    explicit NtuplerMod(const edm::ParameterSet &iConfig);
    ~NtuplerMod();

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  private:
    virtual void beginJob();
    virtual void endJob();
    virtual void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup);

    virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup);
    virtual void endRun  (const edm::Run &iRun, const edm::EventSetup &iSetup);
    virtual void beginLuminosityBlock(const edm::LuminosityBlock &iLumi, const edm::EventSetup &iSetup);
    virtual void endLuminosityBlock  (const edm::LuminosityBlock &iLumi, const edm::EventSetup &iSetup);

    // specify trigger paths of interest
    void setTriggers();
    
    // initialization from HLT menu; needs to be called on every change in HLT menu
    void initHLT(const edm::TriggerResults&, const edm::TriggerNames&);
    
    //
    void separatePileUp(const edm::Event &iEvent, const reco::Vertex &pv);


    //--------------------------------------------------------------------------------------------------
    //  data members
    //==================================================================================================   
    bool fSkipOnHLTFail;
    
    // variables to handle triggers
    edm::ParameterSetID fTriggerNamesID;
    edm::InputTag	fHLTTag;
    edm::InputTag       fHLTObjTag;
    edm::EDGetTokenT<edm::TriggerResults> fHLTTag_token;
    edm::EDGetTokenT<trigger::TriggerEvent> fHLTObjTag_token;
    edm::EDGetTokenT<reco::PFCandidateCollection> fPFCandName_token;
    edm::EDGetTokenT<reco::VertexCollection> fPVName_token;
    std::string         fHLTFile;

    std::vector<const reco::PFCandidate*> fPFNoPU;
    std::vector<const reco::PFCandidate*> fPFPU;
    
    // AOD collection names
    std::string fPVName;
    std::string fPFCandName;
      
    // bacon fillers
    baconhep::FillerEventInfo *fFillerEvtInfo;
    baconhep::FillerGenInfo   *fFillerGenInfo;
    baconhep::FillerVertex    *fFillerPV;
    baconhep::FillerElectron  *fFillerEle;
    baconhep::FillerMuon      *fFillerMuon;
    baconhep::FillerPhoton    *fFillerPhoton;
    baconhep::FillerPF        *fFillerPF;    
    baconhep::FillerJet       **fFillerJet;
    
    baconhep::TTrigger        *fTrigger;
    
    bool fIsActiveEvtInfo;
    bool fIsActiveGenInfo;
    bool fIsActivePV;
    bool fIsActiveEle;
    bool fIsActiveMuon;
    bool fIsActivePhoton;
    bool fIsActiveJet;
    bool fIsActivePF;
    
    // Objects and arrays for output file
    std::string              fOutputName;
    TFile                   *fOutputFile;
    TH1D                    *fTotalEvents;
    TTree                   *fEventTree;
    baconhep::TEventInfo    *fEvtInfo;
    baconhep::TGenEventInfo *fGenEvtInfo;
    TClonesArray            *fGenParArr;
    TClonesArray	    *fEleArr;
    TClonesArray	    *fMuonArr;
    TClonesArray	    **fJetArr;
    TClonesArray	    *fPhotonArr;
    TClonesArray	    *fPVArr;
    TClonesArray	    *fPFParArr;
};

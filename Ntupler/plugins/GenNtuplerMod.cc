#include "GenNtuplerMod.hh"

// bacon classes and constants
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
//#include "BaconAna/DataFormats/interface/TGenJet.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconProd/Ntupler/interface/FillerGenInfo.hh"

// tools to parse HLT name patterns
#include <boost/foreach.hpp>
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

// data format classes
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

// ROOT classes
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMath.h>


//--------------------------------------------------------------------------------------------------
GenNtuplerMod::GenNtuplerMod(const edm::ParameterSet &iConfig):
  fFillerGenInfo     (0),
  fIsActiveGenInfo   (false),
  fOutputName        (iConfig.getUntrackedParameter<std::string>("outputName", "ntuple.root")),
  fOutputFile        (0),
  fGenEvtInfo        (0),
  fGenParArr         (0)//,
  //fGenJetArr         (0)
{
  // Don't write TObject part of the objects
  baconhep::TGenEventInfo::Class()->IgnoreTObjectStreamer();
  baconhep::TGenParticle::Class()->IgnoreTObjectStreamer();
  //baconhep::TGenJet::Class()->IgnoreTObjectStreamer();
  //
  // Set up bacon objects and configure fillers
  //   
  if(iConfig.existsAs<edm::ParameterSet>("GenInfo",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("GenInfo"));
    fIsActiveGenInfo = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActiveGenInfo) {
      fGenEvtInfo    = new baconhep::TGenEventInfo();                   assert(fGenEvtInfo);
      fGenParArr     = new TClonesArray("baconhep::TGenParticle",5000); assert(fGenParArr);
      //fGenJetArr     = new TClonesArray("baconhep::TGenJet",5000); assert(fGenJetArr);
      fFillerGenInfo = new baconhep::FillerGenInfo(cfg,consumesCollector());                assert(fFillerGenInfo);
    }
  }
}

//--------------------------------------------------------------------------------------------------
GenNtuplerMod::~GenNtuplerMod()
{
  delete fFillerGenInfo;
  delete fGenEvtInfo;
  delete fGenParArr;
  //delete fGenJetArr;
}
//--------------------------------------------------------------------------------------------------
void GenNtuplerMod::beginJob()
{  
  //
  // Create output file, trees, and histograms
  //
  fOutputFile  = new TFile(fOutputName.c_str(), "RECREATE");
  fTotalEvents = new TH1D("TotalEvents","TotalEvents",1,-10,10);
  fEventTree   = new TTree("Events","Events");
  
  if(fIsActiveGenInfo) {
    fEventTree->Branch("GenEvtInfo",fGenEvtInfo);
    fEventTree->Branch("GenParticle",&fGenParArr);
    //fEventTree->Branch("GenJet",&fGenJetArr);
  }
 
}

//--------------------------------------------------------------------------------------------------
void GenNtuplerMod::endJob() 
{
  //
  // Save to ROOT file
  //
  //fEventTree->Print();
  fOutputFile->cd();
  fTotalEvents->Write();
  fOutputFile->Write();
  fOutputFile->Close();
}
//--------------------------------------------------------------------------------------------------
void GenNtuplerMod::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  fTotalEvents->Fill(1);

  if(fIsActiveGenInfo) {
    fGenParArr->Clear();
    //fGenJetArr->Clear();
    //fFillerGenInfo->fill(fGenEvtInfo, fGenParArr,fGenJetArr, iEvent);
    fFillerGenInfo->fill(fGenEvtInfo, fGenParArr, iEvent);
  }
     
  fEventTree->Fill();
}

void GenNtuplerMod::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){}
void GenNtuplerMod::endRun  (const edm::Run& iRun, const edm::EventSetup& iSetup){}
void GenNtuplerMod::beginLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){}
void GenNtuplerMod::endLuminosityBlock  (const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){}


//define this as a plug-in
DEFINE_FWK_MODULE(GenNtuplerMod);

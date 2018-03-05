#include "NtuplerMod.hh"

// bacon classes and constants
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

#include "BaconProd/Ntupler/interface/FillerEventInfo.hh"
#include "BaconProd/Ntupler/interface/FillerGenInfo.hh"
#include "BaconProd/Ntupler/interface/FillerVertex.hh"
#include "BaconProd/Ntupler/interface/FillerMuon.hh"
#include "BaconProd/Ntupler/interface/FillerElectron.hh"
#include "BaconProd/Ntupler/interface/FillerPhoton.hh"
#include "BaconProd/Ntupler/interface/FillerJet.hh"
#include "BaconProd/Ntupler/interface/FillerPF.hh"

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
NtuplerMod::NtuplerMod(const edm::ParameterSet &iConfig):
  fSkipOnHLTFail     (iConfig.getUntrackedParameter<bool>("skipOnHLTFail",false)),
  fUseAOD            (iConfig.getUntrackedParameter<bool>("useAOD",false)),
  fHLTTag            ("TriggerResults","","HLT"),
  fHLTObjTag         ("slimmedPatTrigger","","RECO"),
  fHLTFile           (iConfig.getUntrackedParameter<std::string>("TriggerFile","HLT")),
  fPVName            (iConfig.getUntrackedParameter<std::string>("edmPVName","offlinePrimaryVertices")),
  fPFCandName        (iConfig.getUntrackedParameter<std::string>("edmPFCandName","particleFlow")),
  fFillerEvtInfo     (0),
  fFillerGenInfo     (0),
  fFillerPV          (0),
  fFillerEle         (0),
  fFillerMuon        (0),
  fFillerPhoton      (0),
  fFillerPF          (0),
  fTrigger           (0),
  fIsActiveEvtInfo   (false),
  fIsActiveGenInfo   (false),
  fIsActivePV        (false),
  fIsActiveEle       (false),
  fIsActiveMuon      (false),
  fIsActivePhoton    (false),
  fIsActivePF        (false),
  fOutputName        (iConfig.getUntrackedParameter<std::string>("outputName", "ntuple.root")),
  fOutputFile        (0),
  fTotalEvents       (0),
  fEventTree         (0),
  fEvtInfo           (0),
  fGenEvtInfo        (0),
  fGenParArr         (0),
  fEleArr            (0),
  fMuonArr           (0),
  fJetArr            (0),
  fPhotonArr         (0),
  fPVArr             (0),
  fPFParArr          (0)
{
  fHLTTag_token = consumes<edm::TriggerResults>(edm::InputTag(fHLTTag));
  fHLTObjTag_token = consumes<trigger::TriggerEvent>(edm::InputTag(fHLTObjTag));
  fTrgObj_token = consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag(fHLTObjTag)); 
  fPFCandName_token = consumes<reco::PFCandidateCollection>(fPFCandName);
  fPVName_token = consumes<reco::VertexCollection>(fPVName);
  // Don't write TObject part of the objects
  baconhep::TEventInfo::Class()->IgnoreTObjectStreamer();
  baconhep::TGenEventInfo::Class()->IgnoreTObjectStreamer();
  baconhep::TGenParticle::Class()->IgnoreTObjectStreamer();
  baconhep::TMuon::Class()->IgnoreTObjectStreamer();
  baconhep::TElectron::Class()->IgnoreTObjectStreamer();
  baconhep::TJet::Class()->IgnoreTObjectStreamer();
  baconhep::TPhoton::Class()->IgnoreTObjectStreamer();
  baconhep::TVertex::Class()->IgnoreTObjectStreamer();
  
  //
  // Set up bacon objects and configure fillers
  // 
  if(iConfig.existsAs<edm::ParameterSet>("Info",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Info"));
    fIsActiveEvtInfo = cfg.getUntrackedParameter<bool>("isActive");
    
    if(fIsActiveEvtInfo) {
      fEvtInfo       = new baconhep::TEventInfo();            assert(fEvtInfo);
      fFillerEvtInfo = new baconhep::FillerEventInfo(cfg, fUseAOD, consumesCollector());    assert(fFillerEvtInfo);
    }
  }

  if(iConfig.existsAs<edm::ParameterSet>("GenInfo",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("GenInfo"));
    fIsActiveGenInfo = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActiveGenInfo) {
      fGenEvtInfo    = new baconhep::TGenEventInfo();                   assert(fGenEvtInfo);
      fGenParArr     = new TClonesArray("baconhep::TGenParticle",5000); assert(fGenParArr);
      fFillerGenInfo = new baconhep::FillerGenInfo(cfg,consumesCollector());                assert(fFillerGenInfo);
    }
  }

  if(iConfig.existsAs<edm::ParameterSet>("PV",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("PV"));
    fIsActivePV = cfg.getUntrackedParameter<bool>("isActive");
    
    // create array and filler even if vertices won't be saved to output (i.e. fIsActivePV == false),
    // because FillerVertex::fill(...) is used to find the event primary vertex
    // (not elegant, but I suppose a dedicated PV finding function can be implemented somewhere...)
    fPVArr    = new TClonesArray("baconhep::TVertex"); assert(fPVArr);
    fFillerPV = new baconhep::FillerVertex(cfg,fUseAOD,consumesCollector());       assert(fFillerPV);
  }
    
  if(iConfig.existsAs<edm::ParameterSet>("Electron",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Electron"));
    fIsActiveEle = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActiveEle) {
      fEleArr    = new TClonesArray("baconhep::TElectron"); assert(fEleArr);
      fFillerEle = new baconhep::FillerElectron(cfg,fUseAOD,consumesCollector());       assert(fFillerEle);
    }
  }  

  if(iConfig.existsAs<edm::ParameterSet>("Muon",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Muon"));
    fIsActiveMuon = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActiveMuon) {
      fMuonArr    = new TClonesArray("baconhep::TMuon"); assert(fMuonArr);
      fFillerMuon = new baconhep::FillerMuon(cfg,fUseAOD,consumesCollector());       assert(fFillerMuon);
    }
  }  

  if(iConfig.existsAs<edm::ParameterSet>("Photon",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Photon"));
    fIsActivePhoton = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActivePhoton) {
      fPhotonArr    = new TClonesArray("baconhep::TPhoton"); assert(fPhotonArr);
      fFillerPhoton = new baconhep::FillerPhoton(cfg);       assert(fFillerPhoton);
    }
  } 
  if(iConfig.existsAs<edm::ParameterSet>("PFCand",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("PFCand"));
    fIsActivePF = cfg.getUntrackedParameter<bool>("isActive");
    if(fIsActivePF) {
      fPFParArr = new TClonesArray("baconhep::TPFPart",5000); assert(fPFParArr);
      fFillerPF = new baconhep::FillerPF(cfg,consumesCollector());                assert(fFillerPF);
    }
  } 
  if(iConfig.existsAs<edm::ParameterSet>("Jet",false)) {
    edm::ParameterSet cfg(iConfig.getUntrackedParameter<edm::ParameterSet>("Jet"));
    fIsActiveJet = cfg.getUntrackedParameter<bool>("isActive");

    if(fIsActiveJet) {
      fFillerJet = new baconhep::FillerJet*[1];
      fJetArr    = new TClonesArray*[1];
      fFillerJet[0] = new baconhep::FillerJet(cfg,consumesCollector());
      assert(fFillerJet[0]);
      fJetArr[0] = new TClonesArray("baconhep::TJet");
      assert(fJetArr[0]);
    }
  }
}

//--------------------------------------------------------------------------------------------------
NtuplerMod::~NtuplerMod()
{
  delete fFillerEvtInfo;
  delete fFillerGenInfo;
  delete fFillerPV;
  delete fFillerEle;
  delete fFillerMuon;
  delete fFillerPhoton;
  
  delete fTrigger;
  delete fEvtInfo;
  delete fGenEvtInfo;
  delete fGenParArr;
  delete fEleArr;
  delete fMuonArr;
  delete fPhotonArr;
  delete fPVArr;
  delete fPFParArr;
  
  if(fIsActiveJet) {
    delete fFillerJet[0];
    delete [] fFillerJet;

    delete fJetArr[0];
    delete [] fJetArr;
  } 
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::beginJob()
{  
  //
  // Create output file, trees, and histograms
  //
  fOutputFile  = new TFile(fOutputName.c_str(), "RECREATE");
  fTotalEvents = new TH1D("TotalEvents","TotalEvents",1,-10,10);
  fEventTree   = new TTree("Events","Events");
  
  if(fIsActiveEvtInfo) { 
    fEventTree->Branch("Info",fEvtInfo); 
  }
  if(fIsActiveGenInfo) {
    fEventTree->Branch("GenEvtInfo",fGenEvtInfo);
    fEventTree->Branch("GenParticle",&fGenParArr);
  }
  if(fIsActiveEle)    { fEventTree->Branch("Electron", &fEleArr); }
  if(fIsActiveMuon)   { fEventTree->Branch("Muon",     &fMuonArr); }
  if(fIsActivePhoton) { fEventTree->Branch("Photon",   &fPhotonArr); }
  if(fIsActivePV)     { fEventTree->Branch("PV",       &fPVArr); }
  if(fIsActiveJet) {
    fEventTree->Branch("AK4CHS", &fJetArr[0]);
  }
  if(fIsActivePF) { fEventTree->Branch("PFPart", &fPFParArr); }
  //
  // Triggers
  //
  setTriggers();
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::endJob() 
{
  //
  // Save to ROOT file
  //
//  fEventTree->Print();
  fOutputFile->cd();
  fTotalEvents->Write();
  fOutputFile->Write();
  fOutputFile->Close();
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::setTriggers()
{
  std::string cmssw_base_src = getenv("CMSSW_BASE"); cmssw_base_src+="/src/";
  fTrigger = new baconhep::TTrigger(cmssw_base_src + fHLTFile);

  std::cout << fTrigger->fRecords.size() << std::endl;
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  fTotalEvents->Fill(1);
  
  edm::Handle<edm::TriggerResults> hTrgRes;
  iEvent.getByToken(fHLTTag_token,hTrgRes);
  assert(hTrgRes.isValid());  
  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*hTrgRes);
  Bool_t config_changed = false;
  if(fTriggerNamesID != triggerNames.parameterSetID()) {
    fTriggerNamesID = triggerNames.parameterSetID();
    config_changed  = true;
  }
  if(config_changed) {
    initHLT(*hTrgRes, triggerNames);
  }
  
  TriggerBits triggerBits;
  for(unsigned int irec=0; irec<fTrigger->fRecords.size(); irec++) {
    if(fTrigger->fRecords[irec].hltPathIndex == (unsigned int)-1) continue;
    if(hTrgRes->accept(fTrigger->fRecords[irec].hltPathIndex)) {
      triggerBits [fTrigger->fRecords[irec].baconTrigBit] = 1;
    }
  }
  if(fSkipOnHLTFail && triggerBits == 0) return;

  if(fIsActiveGenInfo) {
    fGenParArr->Clear();
    fFillerGenInfo->fill(fGenEvtInfo, fGenParArr, iEvent);
  }
   
  fPVArr->Clear();
  int nvertices = 0;
  const reco::Vertex *pv = fFillerPV->fill(fPVArr, nvertices, iEvent);
  assert(pv);
  
  //separatePileUp(iEvent, *pv);
  
  if(fIsActiveEvtInfo) {
    fFillerEvtInfo->fill(fEvtInfo, iEvent, *pv, (nvertices>0), triggerBits);//,fSusyGen);
  }
  
  edm::Handle<trigger::TriggerEvent> hTrgEvt;
  edm::Handle<pat::TriggerObjectStandAloneCollection> hTrgObjs;
  //edm::Handle  hTrgObjs;

  //const trigger::TriggerEvent*                  hTrgEvtDummy  = 0; 
  pat::TriggerObjectStandAloneCollection* hTrgObjsDummy = 0; 

  iEvent.getByToken(fHLTObjTag_token,hTrgEvt);  
  iEvent.getByToken(fTrgObj_token,hTrgObjs);

  if (!fUseAOD) {
    hTrgObjsDummy = new pat::TriggerObjectStandAloneCollection(*hTrgObjs);
  }

  if(fIsActiveEle) {
    fEleArr->Clear();
    if (fUseAOD) { fFillerEle->fill(fEleArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgEvt); }
    else         { fFillerEle->fill(fEleArr, iEvent, iSetup, *pv, triggerNames, *hTrgRes, fTrigger->fRecords, *hTrgObjsDummy); }
  }

  if(fIsActiveMuon) {
    fMuonArr->Clear();  
    if (fUseAOD) { fFillerMuon->fill(fMuonArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgEvt); }
    else         { fFillerMuon->fill(fMuonArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgObjs); }    
  }

  if(fIsActivePhoton) {
    fPhotonArr->Clear();  
    if (fUseAOD) { fFillerPhoton->fill(fPhotonArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgEvt); }
    //else         { fFillerPhoton->fill(fPhotonArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgObjs); }
  }
  if(fIsActivePF) { 
    fPFParArr->Clear();
    if (fUseAOD) { fFillerPF->fill(fPFParArr,fPVArr,iEvent); }
    else         { fFillerPF->fillMiniAOD(fPFParArr, fPVArr, iEvent); }
  }
  if(fIsActiveJet) {
    fJetArr[0]->Clear();
    fFillerJet[0]->fill(fJetArr[0], iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgEvt);
  }

  fEventTree->Fill();
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::initHLT(const edm::TriggerResults& result, const edm::TriggerNames& triggerNames)
{
  for(unsigned int irec=0; irec<fTrigger->fRecords.size(); irec++) {
    fTrigger->fRecords[irec].hltPathName  = "";
    fTrigger->fRecords[irec].hltPathIndex = (unsigned int)-1;
    const std::string pattern = fTrigger->fRecords[irec].hltPattern;
    if(edm::is_glob(pattern)) {  // handle pattern with wildcards (*,?)
      std::vector<std::vector<std::string>::const_iterator> matches = edm::regexMatch(triggerNames.triggerNames(), pattern);
      if(matches.empty()) {
        std::cout << "requested pattern [" << pattern << "] does not match any HLT paths" << std::endl;
      } else {
        BOOST_FOREACH(std::vector<std::string>::const_iterator match, matches) {
          fTrigger->fRecords[irec].hltPathName = *match;
        }
      }
    } else {  // take full HLT path name given
      fTrigger->fRecords[irec].hltPathName = pattern;
    }
    // Retrieve index in trigger menu corresponding to HLT path
    unsigned int index = triggerNames.triggerIndex(fTrigger->fRecords[irec].hltPathName);
    if(index < result.size()) {  // check for valid index
      fTrigger->fRecords[irec].hltPathIndex = index;
    }
  }
}

//--------------------------------------------------------------------------------------------------
/*void NtuplerMod::separatePileUp(const edm::Event &iEvent, const reco::Vertex &pv)
{
  // recipe from Matthew Chan

  // Get PF-candidates collection
  edm::Handle<reco::PFCandidateCollection> hPFCandProduct;
  iEvent.getByToken(fPFCandName_token,hPFCandProduct);
  assert(hPFCandProduct.isValid());
  const reco::PFCandidateCollection *pfCandCol = hPFCandProduct.product();  
  
  // Get vertex collection
  edm::Handle<reco::VertexCollection> hVertexProduct;
  iEvent.getByToken(fPVName_token,hVertexProduct);
  assert(hVertexProduct.isValid());
  const reco::VertexCollection *pvCol = hVertexProduct.product();
  
  fPFNoPU.clear();
  fPFPU.clear();
  
  for(reco::PFCandidateCollection::const_iterator iP = pfCandCol->begin(); iP!=pfCandCol->end(); ++iP) {
    if(iP->particleId() == reco::PFCandidate::h) {  // charged hadrons
      if(iP->trackRef().isNonnull() && pv.trackWeight(iP->trackRef())>0) {
        // charged hadrons with track used to compute PV
	fPFNoPU.push_back(&(*iP)); 
      
      } else {
        // Find closest vertex to charged hadron's vertex source
	bool vertexFound = false;
	const reco::Vertex *closestVtx = 0;
	double dzmin = 10000;
	
	for(reco::VertexCollection::const_iterator iV = pvCol->begin(); iV!=pvCol->end(); ++iV) {
	  if(iP->trackRef().isNonnull() && iV->trackWeight(iP->trackRef())>0) {
	    vertexFound = true;
	    closestVtx  = &(*iV);
	    break;
	  }
	  
	  double dz = fabs(iP->vertex().z() - iV->z());
	  if(dz < dzmin) {
	    closestVtx = &(*iV);
	    dzmin      = dz;
	  }
	}
	
	if(vertexFound || closestVtx != &pv) {
	  fPFPU.push_back(&(*iP));
	} else {
	  fPFNoPU.push_back(&(*iP));  // Note: when no associated vertex found, assume to come from PV
	}
      }
      
    } else {  // all non-charged-hadron PFCandidates are considered to be from PV
      fPFNoPU.push_back(&(*iP));
    }
  }
}
*/
//--------------------------------------------------------------------------------------------------
void NtuplerMod::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void NtuplerMod::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){}
void NtuplerMod::endRun  (const edm::Run& iRun, const edm::EventSetup& iSetup){}
void NtuplerMod::beginLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){}
void NtuplerMod::endLuminosityBlock  (const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){}


//define this as a plug-in
DEFINE_FWK_MODULE(NtuplerMod);

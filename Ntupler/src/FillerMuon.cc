#include "BaconProd/Ntupler/interface/FillerMuon.hh"
#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMath.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerMuon::FillerMuon(const edm::ParameterSet &iConfig):
  fMinPt     (iConfig.getUntrackedParameter<double>("minPt",0)),
  fMuonName  (iConfig.getUntrackedParameter<std::string>("edmName","muons")),
  fPFCandName(iConfig.getUntrackedParameter<std::string>("edmPFCandName","particleFlow")),
  fTrackName (iConfig.getUntrackedParameter<std::string>("edmTrackName","generalTracks")),
  fSaveTracks(iConfig.getUntrackedParameter<bool>("doSaveTracks",false)),
  fTrackMinPt(iConfig.getUntrackedParameter<double>("minTrackPt",20))
{}

//--------------------------------------------------------------------------------------------------
FillerMuon::~FillerMuon(){}

//--------------------------------------------------------------------------------------------------
void FillerMuon::fill(TClonesArray *array,
                      const edm::Event &iEvent, const edm::EventSetup &iSetup, const reco::Vertex &pv, 
		      const std::vector<TriggerRecord> &triggerRecords,
		      const trigger::TriggerEvent &triggerEvent)
{
  assert(array);
  
  // Get muon collection
  edm::Handle<reco::MuonCollection> hMuonProduct;
  iEvent.getByLabel(fMuonName,hMuonProduct);
  assert(hMuonProduct.isValid());
  const reco::MuonCollection *muonCol = hMuonProduct.product();
  
  // Get PF-candidates collection
  edm::Handle<reco::PFCandidateCollection> hPFCandProduct;
  iEvent.getByLabel(fPFCandName,hPFCandProduct);
  assert(hPFCandProduct.isValid());
  const reco::PFCandidateCollection *pfCandCol = hPFCandProduct.product();
  
  // Get track collection
  edm::Handle<reco::TrackCollection> hTrackProduct;
  iEvent.getByLabel(fTrackName,hTrackProduct);
  assert(hTrackProduct.isValid());
  const reco::TrackCollection *trackCol = hTrackProduct.product();
  
  // Track builder for computing 3D impact parameter
  edm::ESHandle<TransientTrackBuilder> hTransientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",hTransientTrackBuilder);
  const TransientTrackBuilder *transientTrackBuilder = hTransientTrackBuilder.product();  

  
  for(reco::MuonCollection::const_iterator itMu = muonCol->begin(); itMu!=muonCol->end(); ++itMu) {    
    
    // muon pT cut
    if(itMu->pt() < fMinPt) continue;
    
    // construct object and place in array
    TClonesArray &rArray = *array;
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TMuon();
    baconhep::TMuon *pMuon = (baconhep::TMuon*)rArray[index];
    
    //
    // Kinematics
    //==============================
    pMuon->pt      = itMu->muonBestTrack()->pt();
    pMuon->eta     = itMu->muonBestTrack()->eta();
    pMuon->phi     = itMu->muonBestTrack()->phi();
    pMuon->ptErr   = itMu->muonBestTrack()->ptError();
    pMuon->q       = itMu->muonBestTrack()->charge();
    pMuon->staPt   = itMu->standAloneMuon().isNonnull() ? itMu->standAloneMuon()->pt()  : 0;
    pMuon->staEta  = itMu->standAloneMuon().isNonnull() ? itMu->standAloneMuon()->eta() : 0;
    pMuon->staPhi  = itMu->standAloneMuon().isNonnull() ? itMu->standAloneMuon()->phi() : 0;   
    
    pMuon->pfPt  = itMu->pfP4().pt();
    pMuon->pfEta = itMu->pfP4().eta();
    pMuon->pfPhi = itMu->pfP4().phi();
    
    //
    // Isolation
    //==============================
    pMuon->trkIso  = itMu->isolationR03().sumPt;
    pMuon->ecalIso = itMu->isolationR03().emEt;
    pMuon->hcalIso = itMu->isolationR03().hadEt;

    pMuon->chHadIso  = itMu->pfIsolationR04().sumChargedHadronPt;
    pMuon->gammaIso  = itMu->pfIsolationR04().sumPhotonEt;
    pMuon->neuHadIso = itMu->pfIsolationR04().sumNeutralHadronEt;
    pMuon->puIso     = itMu->pfIsolationR04().sumPUPt;

    // PF-Weighted

    // PUPPI
    
    //
    // Impact Parameter
    //==============================
    pMuon->d0 = (-1)*(itMu->muonBestTrack()->dxy(pv.position()));  // note: d0 = -dxy
    pMuon->dz = itMu->muonBestTrack()->dz(pv.position());
    
    const reco::TransientTrack &tt = transientTrackBuilder->build(itMu->muonBestTrack());
    const double thesign = (pMuon->d0 >= 0) ? 1. : -1.;
    const std::pair<bool,Measurement1D> &ip3d = IPTools::absoluteImpactParameter3D(tt,pv);
    pMuon->sip3d = ip3d.first ? thesign*ip3d.second.value() / ip3d.second.error() : -999.;
    
    //
    // Identification
    //==============================
    pMuon->tkNchi2  = itMu->innerTrack().isNonnull() ? itMu->innerTrack()->normalizedChi2()  : -999.;
    pMuon->muNchi2  = itMu->isGlobalMuon()           ? itMu->globalTrack()->normalizedChi2() : -999.;    
    pMuon->trkKink  = itMu->combinedQuality().trkKink;
    pMuon->glbKink  = itMu->combinedQuality().glbKink;        
    pMuon->typeBits = itMu->type();
        
    pMuon->selectorBits=0;
    if(muon::isGoodMuon(*itMu,muon::All))                                    pMuon->selectorBits |= baconhep::kAll;
    if(muon::isGoodMuon(*itMu,muon::AllGlobalMuons))                         pMuon->selectorBits |= baconhep::kAllGlobalMuons;
    if(muon::isGoodMuon(*itMu,muon::AllStandAloneMuons))                     pMuon->selectorBits |= baconhep::kAllStandAloneMuons;
    if(muon::isGoodMuon(*itMu,muon::AllTrackerMuons))                        pMuon->selectorBits |= baconhep::kAllTrackerMuons;
    if(muon::isGoodMuon(*itMu,muon::TrackerMuonArbitrated))                  pMuon->selectorBits |= baconhep::kTrackerMuonArbitrated;
    if(muon::isGoodMuon(*itMu,muon::AllArbitrated))                          pMuon->selectorBits |= baconhep::kAllArbitrated;
    if(muon::isGoodMuon(*itMu,muon::GlobalMuonPromptTight))                  pMuon->selectorBits |= baconhep::kGlobalMuonPromptTight;
    if(muon::isGoodMuon(*itMu,muon::TMLastStationLoose))                     pMuon->selectorBits |= baconhep::kTMLastStationLoose;
    if(muon::isGoodMuon(*itMu,muon::TMLastStationTight))                     pMuon->selectorBits |= baconhep::kTMLastStationTight;
    if(muon::isGoodMuon(*itMu,muon::TM2DCompatibilityLoose))                 pMuon->selectorBits |= baconhep::kTM2DCompatibilityLoose;
    if(muon::isGoodMuon(*itMu,muon::TM2DCompatibilityTight))                 pMuon->selectorBits |= baconhep::kTM2DCompatibilityTight;
    if(muon::isGoodMuon(*itMu,muon::TMOneStationLoose))                      pMuon->selectorBits |= baconhep::kTMOneStationLoose;
    if(muon::isGoodMuon(*itMu,muon::TMOneStationTight))                      pMuon->selectorBits |= baconhep::kTMOneStationTight;
    if(muon::isGoodMuon(*itMu,muon::TMLastStationOptimizedLowPtLoose))       pMuon->selectorBits |= baconhep::kTMLastStationOptimizedLowPtLoose;
    if(muon::isGoodMuon(*itMu,muon::TMLastStationOptimizedLowPtTight))       pMuon->selectorBits |= baconhep::kTMLastStationOptimizedLowPtTight;
    if(muon::isGoodMuon(*itMu,muon::GMTkChiCompatibility))                   pMuon->selectorBits |= baconhep::kGMTkChiCompatibility;
    if(muon::isGoodMuon(*itMu,muon::GMStaChiCompatibility))                  pMuon->selectorBits |= baconhep::kGMStaChiCompatibility;
    if(muon::isGoodMuon(*itMu,muon::GMTkKinkTight))                          pMuon->selectorBits |= baconhep::kGMTkKinkTight;
    if(muon::isGoodMuon(*itMu,muon::TMLastStationAngLoose))                  pMuon->selectorBits |= baconhep::kTMLastStationAngLoose;
    if(muon::isGoodMuon(*itMu,muon::TMLastStationAngTight))                  pMuon->selectorBits |= baconhep::kTMLastStationAngTight;
    if(muon::isGoodMuon(*itMu,muon::TMOneStationAngLoose))                   pMuon->selectorBits |= baconhep::kTMOneStationAngLoose;
    if(muon::isGoodMuon(*itMu,muon::TMOneStationAngTight))                   pMuon->selectorBits |= baconhep::kTMOneStationAngTight;
    if(muon::isGoodMuon(*itMu,muon::TMLastStationOptimizedBarrelLowPtLoose)) pMuon->selectorBits |= baconhep::kTMLastStationOptimizedBarrelLowPtLoose;
    if(muon::isGoodMuon(*itMu,muon::TMLastStationOptimizedBarrelLowPtTight)) pMuon->selectorBits |= baconhep::kTMLastStationOptimizedBarrelLowPtTight;
    if(muon::isGoodMuon(*itMu,muon::RPCMuLoose))                             pMuon->selectorBits |= baconhep::kRPCMuLoose;

    pMuon->pogIDBits=0;
    if(muon::isLooseMuon(*itMu))      pMuon->pogIDBits |= baconhep::kPOGLooseMuon;
    if(muon::isTightMuon(*itMu, pv))  pMuon->pogIDBits |= baconhep::kPOGTightMuon;
    if(muon::isSoftMuon(*itMu, pv))   pMuon->pogIDBits |= baconhep::kPOGSoftMuon;
    if(muon::isHighPtMuon(*itMu, pv)) pMuon->pogIDBits |= baconhep::kPOGHighPtMuon;

    pMuon->nValidHits = itMu->isGlobalMuon()           ? itMu->globalTrack()->hitPattern().numberOfValidMuonHits()       : 0;        
    pMuon->nTkHits    = itMu->innerTrack().isNonnull() ? itMu->innerTrack()->hitPattern().numberOfValidTrackerHits()     : 0;
    pMuon->nPixHits   = itMu->innerTrack().isNonnull() ? itMu->innerTrack()->hitPattern().numberOfValidPixelHits()       : 0;
    pMuon->nTkLayers  = itMu->innerTrack().isNonnull() ? itMu->innerTrack()->hitPattern().trackerLayersWithMeasurement() : 0;
    pMuon->nPixLayers = itMu->innerTrack().isNonnull() ? itMu->innerTrack()->hitPattern().pixelLayersWithMeasurement()   : 0;
    pMuon->nMatchStn  = itMu->numberOfMatchedStations();
    
    pMuon->trkID = -1;
    if(itMu->innerTrack().isNonnull()) {
      int trkIndex = -1;
      for(reco::TrackCollection::const_iterator itTrk = trackCol->begin(); itTrk!=trackCol->end(); ++itTrk) {
        trkIndex++;
	if(itMu->innerTrack().get() == &(*itTrk)) {
	  pMuon->trkID = index;
	  break;
	}
      }
    }
    
    pMuon->hltMatchBits = TriggerTools::matchHLT(pMuon->eta, pMuon->phi, triggerRecords, triggerEvent);
  }
  
  
  if(fSaveTracks) {    
    int trkIndex = -1;
    for(reco::TrackCollection::const_iterator itTrk = trackCol->begin(); itTrk!=trackCol->end(); ++itTrk) {
      trkIndex++;
      
      // check track is not a muon
      bool isMuon = false;
      for(reco::MuonCollection::const_iterator itMu = muonCol->begin(); itMu!=muonCol->end(); ++itMu) {
        if(itMu->innerTrack().isNonnull() && itMu->innerTrack().get() == &(*itTrk)) {
	  isMuon = true;
	  break;
	}
      }
      if(isMuon) continue;

      // track pT cut
      if(itTrk->pt() < fTrackMinPt) continue;    
      
      
      TClonesArray &rArray = *array;
      assert(rArray.GetEntries() < rArray.GetSize());
      const int index = rArray.GetEntries();
      new(rArray[index]) baconhep::TMuon();
      baconhep::TMuon *pMuon = (baconhep::TMuon*)rArray[index];
    
      //
      // Kinematics
      //==============================
      pMuon->pt      = itTrk->pt();
      pMuon->eta     = itTrk->eta();
      pMuon->phi     = itTrk->phi();
      pMuon->ptErr   = itTrk->ptError();
      pMuon->q       = itTrk->charge();
      pMuon->staPt   = 0;
      pMuon->staEta  = 0;
      pMuon->staPhi  = 0;   
    
      pMuon->pfPt  = 0;
      pMuon->pfEta = 0;
      pMuon->pfPhi = 0;
      for(reco::PFCandidateCollection::const_iterator itPF = pfCandCol->begin(); itPF!=pfCandCol->end(); ++itPF) {
        if(itPF->trackRef().isNonnull() && &(*itTrk) == itPF->trackRef().get()) {
	  pMuon->pfPt  = itPF->pt();
	  pMuon->pfEta = itPF->eta();
	  pMuon->pfPhi = itPF->phi();
        }
      }
    
      //
      // Isolation
      //==============================
      pMuon->trkIso  = -1;
      pMuon->ecalIso = -1;
      pMuon->hcalIso = -1;
/*    
      computeIso(*itTrk, 0.3, pfNoPU, pfPU,
                 pMuon->chHadIso,
                 pMuon->gammaIso,
                 pMuon->neuHadIso,
                 pMuon->puIso);
*/    
      //
      // Impact Parameter
      //==============================
      pMuon->d0 = (-1)*(itTrk->dxy(pv.position()));  // note: d0 = -dxy
      pMuon->dz =  itTrk->dz(pv.position());

      const reco::TransientTrack &tt = transientTrackBuilder->build(&(*itTrk));
      const double thesign = (pMuon->d0 >= 0) ? 1. : -1.;
      const std::pair<bool,Measurement1D> &ip3d = IPTools::absoluteImpactParameter3D(tt,pv);
      pMuon->sip3d = ip3d.first ? thesign*ip3d.second.value() / ip3d.second.error() : -999.;      

      //
      // Identification
      //==============================
      pMuon->tkNchi2      = itTrk->normalizedChi2();
      pMuon->muNchi2      = -999.;    
      pMuon->trkKink      = 0;
      pMuon->glbKink      = 0;        
      pMuon->typeBits     = 0;
      pMuon->selectorBits = 0;
      pMuon->nValidHits   = 0;        
      pMuon->nTkHits      = itTrk->hitPattern().numberOfValidTrackerHits();
      pMuon->nPixHits     = itTrk->hitPattern().numberOfValidPixelHits();
      pMuon->nTkLayers    = itTrk->hitPattern().trackerLayersWithMeasurement();
      pMuon->nPixLayers   = itTrk->hitPattern().pixelLayersWithMeasurement();
      pMuon->nMatchStn    = 0;
      pMuon->trkID        = 0;//trkIndex;
      pMuon->hltMatchBits = 0;//TriggerTools::matchHLT(pMuon->eta, pMuon->phi, triggerRecords, triggerEvent);
    }    
  } 
}

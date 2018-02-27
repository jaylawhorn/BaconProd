#include "BaconProd/Ntupler/interface/FillerElectron.hh"
#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <utility>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerElectron::FillerElectron(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC):
  fMinPt       (iConfig.getUntrackedParameter<double>("minPt",7)),
  fEleName     (iConfig.getUntrackedParameter<std::string>("edmName","gedGsfElectrons")),
  fPFCandName  (iConfig.getUntrackedParameter<std::string>("edmPFCandName","particleFlow")),
  fTrackName   (iConfig.getUntrackedParameter<std::string>("edmTrackName","generalTracks")),
  fBeamspotName (iConfig.getUntrackedParameter<std::string>("edmBeamspotName","beamspot")),
  fConvName    (iConfig.getUntrackedParameter<std::string>("edmConversionName","allConversions")),
  fSCName      (iConfig.getUntrackedParameter<std::string>("edmSCName","particleFlowEGamma"))
{
  fEleName_token = iC.consumes<reco::GsfElectronCollection>(fEleName);
  fPFCandName_token = iC.consumes<reco::PFCandidateCollection>(fPFCandName);
  fTrackName_token = iC.consumes<reco::TrackCollection>(fTrackName);
  fBeamspotName_token = iC.consumes<reco::BeamSpot>(fBeamspotName);
  fConvName_token = iC.consumes<reco::ConversionCollection>(fConvName);
  fSCName_token = iC.consumes<reco::SuperClusterCollection>(fSCName);
  std::string cmssw_base_src = getenv("CMSSW_BASE");
  cmssw_base_src += "/src/";
}

//--------------------------------------------------------------------------------------------------
FillerElectron::~FillerElectron(){}

//--------------------------------------------------------------------------------------------------
void FillerElectron::fill(TClonesArray *array,	    
	                  const edm::Event &iEvent, const edm::EventSetup &iSetup,      
	                  const reco::Vertex &pv,
			  const std::vector<TriggerRecord> &triggerRecords,
			  const trigger::TriggerEvent &triggerEvent)
{
  assert(array);
  
  // Get electron collection
  edm::Handle<reco::GsfElectronCollection> hEleProduct;
  iEvent.getByToken(fEleName_token,hEleProduct);
  assert(hEleProduct.isValid());
  const reco::GsfElectronCollection *eleCol = hEleProduct.product();

  // Get PF-candidates collection
  edm::Handle<reco::PFCandidateCollection> hPFCandProduct;
  iEvent.getByToken(fPFCandName_token,hPFCandProduct);
  assert(hPFCandProduct.isValid());
  const reco::PFCandidateCollection *pfCandCol = hPFCandProduct.product();
  
  // Get track collection
  edm::Handle<reco::TrackCollection> hTrackProduct;
  iEvent.getByToken(fTrackName_token,hTrackProduct);
  assert(hTrackProduct.isValid());
  const reco::TrackCollection *trackCol = hTrackProduct.product();
  
  // Get beam spot
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(fBeamspotName_token,theBeamSpot);
 
  // Get conversions collection
  edm::Handle<reco::ConversionCollection> hConvProduct;
  iEvent.getByToken(fConvName_token,hConvProduct);
  assert(hConvProduct.isValid());

  // Get SuperCluster collection
  edm::Handle<reco::SuperClusterCollection> hSCProduct;
  iEvent.getByToken(fSCName_token,hSCProduct);
  assert(hSCProduct.isValid());
  const reco::SuperClusterCollection *scCol = hSCProduct.product();

  // Track builder for computing 3D impact parameter    
  edm::ESHandle<TransientTrackBuilder> hTransientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",hTransientTrackBuilder);
  const TransientTrackBuilder *transientTrackBuilder = hTransientTrackBuilder.product();


  for(reco::GsfElectronCollection::const_iterator itEle = eleCol->begin(); itEle!=eleCol->end(); ++itEle) {

    const reco::GsfTrackRef gsfTrack = itEle->gsfTrack();
    const reco::SuperClusterRef sc   = itEle->superCluster();

    // electron pT cut
    if(itEle->pt() < fMinPt) continue;

    // construct object and place in array    
    TClonesArray &rElectronArr = *array;
    assert(rElectronArr.GetEntries() < rElectronArr.GetSize());
    const int index = rElectronArr.GetEntries();  
    new(rElectronArr[index]) baconhep::TElectron();
    baconhep::TElectron *pElectron = (baconhep::TElectron*)rElectronArr[index];
    
    //
    // Kinematics
    //==============================    
    pElectron->pt         = itEle->pt();
    pElectron->eta        = itEle->eta();
    pElectron->phi        = itEle->phi();
    pElectron->q          = itEle->charge();
    pElectron->ecalEnergy = itEle->ecalEnergy();
    pElectron->scEt       = (sc->energy())*(sc->position().Rho())/(sc->position().R());
    pElectron->scEta      = sc->eta();
    pElectron->scPhi      = sc->phi();

/** double check this...seems like too many non-matches **/
    pElectron->pfPt  = 0;
    pElectron->pfEta = 0;
    pElectron->pfPhi = 0;
    for(reco::PFCandidateCollection::const_iterator itPF = pfCandCol->begin(); itPF!=pfCandCol->end(); ++itPF) {
      if( (itEle->track().isNonnull() && itPF->trackRef().isNonnull() && itEle->track() == itPF->trackRef()) ||
          (itPF->gsfTrackRef().isNonnull() && itPF->gsfTrackRef() == gsfTrack) ) {
        
	pElectron->pfPt  = itPF->pt();
	pElectron->pfEta = itPF->eta();
	pElectron->pfPhi = itPF->phi();  
      }
    }

    //
    // Isolation
    //==============================
    pElectron->trkIso  = itEle->dr03TkSumPt();
    pElectron->ecalIso = itEle->dr03EcalRecHitSumEt();
    pElectron->hcalIso = itEle->dr03HcalTowerSumEt();
    pElectron->hcalDepth1Iso = itEle->dr03HcalDepth1TowerSumEt();

    pElectron->chHadIso  = itEle->pfIsolationVariables().sumChargedHadronPt;
    pElectron->gammaIso  = itEle->pfIsolationVariables().sumPhotonEt;
    pElectron->neuHadIso = itEle->pfIsolationVariables().sumNeutralHadronEt;
    pElectron->puIso     = itEle->pfIsolationVariables().sumPUPt;

    //
    // Impact Parameter
    //==============================
    if(gsfTrack.isNonnull()) { 
      pElectron->d0 = (-1)*(gsfTrack->dxy(pv.position()));
      pElectron->dz = gsfTrack->dz(pv.position());

      pElectron->chi2 = gsfTrack->chi2();
      pElectron->ndof = gsfTrack->ndof();
      pElectron->npixmatch = -9;

/** double check recipe **/
      const reco::TransientTrack &tt = transientTrackBuilder->build(gsfTrack);
      const double gsfsign = (pElectron->d0 >= 0) ? 1. : -1.;
      const std::pair<bool,Measurement1D> &ip3d = IPTools::absoluteImpactParameter3D(tt,pv);
      pElectron->sip3d = ip3d.first ? gsfsign*ip3d.second.value() / ip3d.second.error() : -999.;
    }
    
    //
    // Identification
    //==============================
    pElectron->sieie    = itEle->full5x5_sigmaIetaIeta();
    pElectron->e1x5     = itEle->full5x5_e1x5();
    pElectron->e2x5     = itEle->full5x5_e2x5Max();
    pElectron->e5x5     = itEle->full5x5_e5x5();
    pElectron->r9       = itEle->full5x5_r9();
    pElectron->hovere   = itEle->hcalOverEcal();
    pElectron->eoverp   = itEle->eSuperClusterOverP();
    pElectron->fbrem    = itEle->fbrem();
    pElectron->dEtaInSeed = dEtaInSeed(*itEle);
    pElectron->dEtaIn   = itEle->deltaEtaSuperClusterTrackAtVtx();
    pElectron->dPhiIn   = itEle->deltaPhiSuperClusterTrackAtVtx();
    
    pElectron->mva = -999;
/** check conversion stuff **/    
    pElectron->isConv = ConversionTools::hasMatchedConversion(*itEle, hConvProduct, theBeamSpot->position());
    
    if(gsfTrack.isNonnull()) {
      pElectron->nMissingHits = gsfTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
    }

    pElectron->typeBits=0;
    if(itEle->ecalDrivenSeed())    pElectron->typeBits |= baconhep::kEcalDriven;
    if(itEle->trackerDrivenSeed()) pElectron->typeBits |= baconhep::kTrackerDriven;
    
    pElectron->fiducialBits=0;
    if(itEle->isEB())        pElectron->fiducialBits |= kIsEB;
    if(itEle->isEE())        pElectron->fiducialBits |= kIsEE;
    if(itEle->isGap())       pElectron->fiducialBits |= kIsGap;
    if(itEle->isEBEEGap())   pElectron->fiducialBits |= kIsEBEEGap;
    if(itEle->isEBGap())     pElectron->fiducialBits |= kIsEBGap;
    if(itEle->isEBEtaGap())  pElectron->fiducialBits |= kIsEBEtaGap;
    if(itEle->isEBPhiGap())  pElectron->fiducialBits |= kIsEBPhiGap;
    if(itEle->isEEGap())     pElectron->fiducialBits |= kIsEEGap;
    if(itEle->isEEDeeGap())  pElectron->fiducialBits |= kIsEEDeeGap;
    if(itEle->isEERingGap()) pElectron->fiducialBits |= kIsEERingGap;

    pElectron->classification = itEle->classification();

    // Obtain a supercluster ID, unique per event. The SC ID is the index in the SC collection.
    pElectron->scID = -1;
    int scIndex = -1;
    for(reco::SuperClusterCollection::const_iterator iS = scCol->begin(); iS!=scCol->end(); ++iS) {
      scIndex++;
      if(itEle->superCluster().get() == &(*iS)) {
        pElectron->scID = scIndex;
	break;
      }
    }

    // Obtain a track ID, unique per event. The track ID is the index in the general tracks collection
    pElectron->trkID = -1;
    if(itEle->closestTrack().isNonnull()) {
      int trkIndex = -1;    
      for(reco::TrackCollection::const_iterator itTrk = trackCol->begin(); itTrk!=trackCol->end(); ++itTrk) {
        trkIndex++;
        if(itEle->closestTrack().get() == &(*itTrk)) {
          pElectron->trkID = trkIndex;
	  break;
        }
      }
    }
    
    pElectron->hltMatchBits = TriggerTools::matchHLT(pElectron->eta, pElectron->phi, triggerRecords, triggerEvent);
  }
}

//--------------------------------------------------------------------------------------------------
double FillerElectron::dEtaInSeed(const reco::GsfElectron& ele) {
  return (ele.superCluster().isNonnull() && ele.superCluster()->seed().isNonnull()) ? 
         ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta() : 
         std::numeric_limits<float>::max();
}

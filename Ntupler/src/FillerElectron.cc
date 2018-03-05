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
FillerElectron::FillerElectron(const edm::ParameterSet &iConfig,const bool useAOD,edm::ConsumesCollector && iC):
  fMinPt       (iConfig.getUntrackedParameter<double>("minPt",7)),
  fEleName     (iConfig.getUntrackedParameter<std::string>("edmName","gedGsfElectrons")),
  fBeamspotName (iConfig.getUntrackedParameter<std::string>("edmBeamspotName","beamspot")),
  fPFCandName  (iConfig.getUntrackedParameter<std::string>("edmPFCandName","particleFlow")),
  fTrackName   (iConfig.getUntrackedParameter<std::string>("edmTrackName","generalTracks")),
  fConvName    (iConfig.getUntrackedParameter<std::string>("edmConversionName","allConversions")),
  fSCName      (iConfig.getUntrackedParameter<edm::InputTag>("edmSCName")),
  fEcalPFClusterIsoMapTag(iConfig.getUntrackedParameter<edm::InputTag>("edmEcalPFClusterIsoMapTag")),
  fHcalPFClusterIsoMapTag(iConfig.getUntrackedParameter<edm::InputTag>("edmHcalPFClusterIsoMapTag")),
  fUseAOD      (useAOD)
{
  if (fUseAOD) fEleName_token = iC.consumes<reco::GsfElectronCollection>(fEleName);
  if (!fUseAOD) fElePATName_token = iC.consumes<pat::ElectronCollection>(fEleName);

  fBeamspotName_token = iC.consumes<reco::BeamSpot>(fBeamspotName);
  fPFCandName_token = iC.consumes<reco::PFCandidateCollection>(fPFCandName);
  fTrackName_token = iC.consumes<reco::TrackCollection>(fTrackName);
  fConvName_token = iC.consumes<reco::ConversionCollection>(fConvName);
  fSCName_token = iC.consumes<reco::SuperClusterCollection>(fSCName);
  fEcalPFClusterIsoMap_token = iC.consumes<edm::ValueMap<float> >(fEcalPFClusterIsoMapTag);
  fHcalPFClusterIsoMap_token = iC.consumes<edm::ValueMap<float> >(fHcalPFClusterIsoMapTag);
  std::string cmssw_base_src = getenv("CMSSW_BASE");
  cmssw_base_src += "/src/";
}

//--------------------------------------------------------------------------------------------------
FillerElectron::~FillerElectron(){}

//--------------------------------------------------------------------------------------------------
//AOD
void FillerElectron::fill(TClonesArray *array,	    
	                  const edm::Event &iEvent, const edm::EventSetup &iSetup,      
	                  const reco::Vertex &pv,
			  //const trigger::Results &triggerResults,
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

// === filler for MINIAOD ===
void FillerElectron::fill(TClonesArray *array,
                          const edm::Event &iEvent, const edm::EventSetup &iSetup,
                          const reco::Vertex &pv,
			  const edm::TriggerNames &triggerNames,
			  const edm::TriggerResults &triggerResults,
                          const std::vector<TriggerRecord> &triggerRecords,
                          pat::TriggerObjectStandAloneCollection &triggerObjects)
{
  assert(array);

  // Get electron collection
  edm::Handle<pat::ElectronCollection> hEleProduct;
  iEvent.getByToken(fElePATName_token,hEleProduct);
  assert(hEleProduct.isValid());
  const pat::ElectronCollection *eleCol = hEleProduct.product();

  // Get supercluster collection
  edm::Handle<reco::SuperClusterCollection> hSCProduct;
  iEvent.getByToken(fSCName_token,hSCProduct);
  assert(hSCProduct.isValid());
  const reco::SuperClusterCollection *scCol = hSCProduct.product();

  for(pat::ElectronCollection::const_iterator itEle = eleCol->begin(); itEle!=eleCol->end(); ++itEle) {

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

    pElectron->pfPt  = 0;
    pElectron->pfEta = 0;
    pElectron->pfPhi = 0;
    double pfDR = 0.05;
    for(unsigned int ipf=0; ipf<itEle->numberOfSourceCandidatePtrs(); ipf++) {
      if(ipf==0 && itEle->pfCandidateRef().isNonnull()) continue;  // PF-candidate collection not in MINIAOD, asking for reference will crash
      const reco::CandidatePtr pfcand = itEle->sourceCandidatePtr(ipf);
      double dR = reco::deltaR(itEle->eta(), itEle->phi(), pfcand->eta(), pfcand->phi());
      if(dR < pfDR) {
        pfDR = dR;
        pElectron->pfPt  = pfcand->pt();
        pElectron->pfEta = pfcand->eta();
        pElectron->pfPhi = pfcand->phi();
      }
    }
    

    //
    // Isolation
    //==============================
    pElectron->trkIso        = itEle->dr03TkSumPt();
    pElectron->ecalIso       = itEle->dr03EcalRecHitSumEt();
    pElectron->hcalIso       = itEle->dr03HcalTowerSumEt();
    pElectron->hcalDepth1Iso = itEle->dr03HcalDepth1TowerSumEt();

    pElectron->chHadIso  = itEle->pfIsolationVariables().sumChargedHadronPt;
    pElectron->gammaIso  = itEle->pfIsolationVariables().sumPhotonEt;
    pElectron->neuHadIso = itEle->pfIsolationVariables().sumNeutralHadronEt;
    pElectron->puIso     = itEle->pfIsolationVariables().sumPUPt;
    
    //pElectron->ecalPFClusIso = itEle->ecalPFClusterIso();
    //pElectron->hcalPFClusIso = itEle->hcalPFClusterIso();


    //
    // Impact Parameter
    //==============================
    if(gsfTrack.isNonnull()) {
      pElectron->d0    = (-1)*(gsfTrack->dxy(pv.position()));  // note: d0 = -dxy
      pElectron->dz    = gsfTrack->dz(pv.position());
      pElectron->sip3d = (itEle->edB(pat::Electron::PV3D) > 0) ? itEle->dB(pat::Electron::PV3D)/itEle->edB(pat::Electron::PV3D) : -999;
      }
    

    //
    // Identification
    //==============================
    pElectron->sieie      = itEle->full5x5_sigmaIetaIeta();
    pElectron->e1x5       = itEle->full5x5_e1x5();
    pElectron->e2x5       = itEle->full5x5_e2x5Max();
    pElectron->e5x5       = itEle->full5x5_e5x5();
    pElectron->r9         = itEle->full5x5_r9();
    pElectron->hovere     = itEle->hcalOverEcal();
    pElectron->eoverp     = itEle->eSuperClusterOverP();
    pElectron->fbrem      = itEle->fbrem();
    pElectron->dEtaInSeed = dEtaInSeed(*itEle);
    pElectron->dEtaIn     = itEle->deltaEtaSuperClusterTrackAtVtx();
    pElectron->dPhiIn     = itEle->deltaPhiSuperClusterTrackAtVtx();

    pElectron->mva = -999;
    pElectron->isConv = !itEle->passConversionVeto();


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
    pElectron->trkID = -1;  // general tracks not in MINIAOD

    pElectron->hltMatchBits = TriggerTools::matchHLT(pElectron->eta, pElectron->phi, iEvent, triggerRecords, triggerNames, triggerResults, triggerObjects);
  }
}


//--------------------------------------------------------------------------------------------------
double FillerElectron::dEtaInSeed(const reco::GsfElectron& ele) {
  return (ele.superCluster().isNonnull() && ele.superCluster()->seed().isNonnull()) ? 
         ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta() : 
         std::numeric_limits<float>::max();
}
double FillerElectron::dEtaInSeed(const pat::Electron& ele) {
  return (ele.superCluster().isNonnull() && ele.superCluster()->seed().isNonnull()) ?
    ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta() :
    std::numeric_limits<float>::max();
}
void FillerElectron::computeIso(double &iEta,double &iPhi, const double extRadius,
				const reco::PFCandidateCollection    &puppi,
				float &out_chHadIso, float &out_gammaIso, float &out_neuHadIso) const
{
  // Muon PF isolation with delta-beta PU correction:
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation
  
  double chHadIso=0, gammaIso=0, neuHadIso=0;
  
  const double ptMin           = 0.1;  
  double intRadius = 0;
  
  for(unsigned int ipf=0; ipf<puppi.size(); ipf++) {
    const reco::PFCandidate pfcand = puppi.at(ipf);    
    bool pPass = true;
    double dr = reco::deltaR(pfcand.eta(), pfcand.phi(), iEta, iPhi);
    if(dr < 0.0001) pPass = false; //Use this to avoid float/double bullshit
    if(pPass) { 
      if     (pfcand.particleId() == reco::PFCandidate::h)     { intRadius = 0;}//intRadiusChHad; }
      else if(pfcand.particleId() == reco::PFCandidate::gamma) { intRadius = 0;}//intRadiusGamma;  }
      else if(pfcand.particleId() == reco::PFCandidate::h0)    { intRadius = 0;}//intRadiusNeuHad; }
            
      if(dr>=extRadius || dr<intRadius) continue;
            
      if     (pfcand.particleId() == reco::PFCandidate::h)                             { chHadIso  += pfcand.pt(); }
      else if(pfcand.particleId() == reco::PFCandidate::gamma && pfcand.pt() > ptMin) { gammaIso  += pfcand.pt(); }
      else if(pfcand.particleId() == reco::PFCandidate::h0    && pfcand.pt() > ptMin) { neuHadIso += pfcand.pt(); }
    }
  }
  // compute PU iso
  out_chHadIso  = chHadIso;
  out_gammaIso  = gammaIso;
  out_neuHadIso = neuHadIso;
}
void FillerElectron::computeIso(double &iEta,double &iPhi, const double extRadius,
				const pat::PackedCandidateCollection    &puppi,
				float &out_chHadIso, float &out_gammaIso, float &out_neuHadIso) const
{
  // Muon PF isolation with delta-beta PU correction:
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation
  
  double chHadIso=0, gammaIso=0, neuHadIso=0;
  
  const double ptMin           = 0.1;  
  //const double intRadiusChHad  = 0.0001;
  //const double intRadiusGamma  = 0.01;
  //const double intRadiusNeuHad = 0.01;
  double intRadius = 0;
  
  for(unsigned int ipf=0; ipf<puppi.size(); ipf++) {
    const pat::PackedCandidate pfcand = puppi.at(ipf);    
    bool pPass = true;
    double dr = reco::deltaR(pfcand.eta(), pfcand.phi(), iEta, iPhi);
    if(dr < 0.0001) pPass = false; //Use this to avoid float/double bullshit
    if(pPass) { 
      if     (abs(pfcand.pdgId()) == 211)     { intRadius = 0;}//intRadiusChHad; }
      else if(abs(pfcand.pdgId()) == 22)      { intRadius = 0;}//intRadiusGamma;  }
      else if(abs(pfcand.pdgId()) == 130)     { intRadius = 0;}//intRadiusNeuHad; }
            
      if(dr>=extRadius || dr<intRadius) continue;
      if     (abs(pfcand.pdgId()) == 211)                        { chHadIso  += pfcand.pt(); }
      else if(abs(pfcand.pdgId()) == 22  && pfcand.pt() > ptMin) { gammaIso  += pfcand.pt(); }
      else if(abs(pfcand.pdgId()) == 130 && pfcand.pt() > ptMin) { neuHadIso += pfcand.pt(); }      
    }
  }
  // compute PU iso
  out_chHadIso  = chHadIso;
  out_gammaIso  = gammaIso;
  out_neuHadIso = neuHadIso;
}

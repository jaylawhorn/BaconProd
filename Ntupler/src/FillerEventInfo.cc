#include "BaconProd/Ntupler/interface/FillerEventInfo.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"
#include <TLorentzVector.h>
#include <string>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerEventInfo::FillerEventInfo(const edm::ParameterSet &iConfig, const bool useAOD, edm::ConsumesCollector && iC):
  fPFCandName (iConfig.getUntrackedParameter<std::string>("edmPFCandName","particleFlow")),
  fPUInfoName (iConfig.getUntrackedParameter<std::string>("edmPileupInfoName","addPileupInfo")),
  fPVName     (iConfig.getUntrackedParameter<std::string>("edmPVName","offlinePrimaryVertices")),
  fBSName     (iConfig.getUntrackedParameter<std::string>("edmBeamspotName","offlineBeamSpot")),
  fMETName     (iConfig.getUntrackedParameter<std::string>("edmMETName","slimmedMETs")),
  fPFMETName   (iConfig.getUntrackedParameter<std::string>("edmPFMETName","slimmedMETs")),
  fPFMETCName  (iConfig.getUntrackedParameter<std::string>("edmPFMETCorrName","pfType1CorrectedMet")),
  fPuppETName  (iConfig.getUntrackedParameter<std::string>("edmPuppETName")),
  fPuppETCName (iConfig.getUntrackedParameter<std::string>("edmPuppETCorrName","pfType1CorrectedMetPuppi")),
  fRhoIsoName (iConfig.getUntrackedParameter<std::string>("edmRhoForIsoName","fixedGridRhoFastjetAll")),
  fRhoJetName (iConfig.getUntrackedParameter<std::string>("edmRhoForJetEnergy","fixedGridRhoFastjetAll")),
  fFillMET    (iConfig.getUntrackedParameter<bool>("doFillMET",true)),
  fFillMETFilters(iConfig.getUntrackedParameter<bool>("doFillMETFilters",true)),
  fUseAOD      (useAOD)
{
  fPUInfoName_token = iC.consumes< std::vector<PileupSummaryInfo> >(fPUInfoName);
  fBSName_token   = iC.consumes<reco::BeamSpot>(fBSName);

  fPFMETName_token   = iC.consumes<reco::PFMETCollection>(fPFMETName);
  fPFMETCName_token  = iC.consumes<reco::PFMETCollection> (fPFMETCName);
  fPuppETName_token  = iC.consumes<reco::PFMETCollection> (fPuppETName);
  fPuppETCName_token = iC.consumes<reco::PFMETCollection> (fPuppETCName);
  fPFMETPATName_token   = iC.consumes<pat::METCollection>(fPFMETName);
  //fPFMETCPATName_token  = iC.consumes<pat::METCollection>(fPFMETCName);
  fPuppETPATName_token  = iC.consumes<pat::METCollection>(fPuppETName);
  //fPuppETCPATName_token = iC.consumes<pat::METCollection>(fPuppETName);

  if (fUseAOD) fPFCandName_token = iC.consumes<reco::PFCandidateCollection> (fPFCandName);
  if(!fUseAOD) fPackCandName_token = iC.consumes<pat::PackedCandidateCollection> (fPFCandName);

  rhoIsoTag_token = iC.consumes<double>(fRhoIsoName);
  rhoJetTag_token = iC.consumes<double>(fRhoJetName);
}

//--------------------------------------------------------------------------------------------------
FillerEventInfo::~FillerEventInfo(){}

//--------------------------------------------------------------------------------------------------
void FillerEventInfo::fill(TEventInfo *evtInfo,
                           const edm::Event &iEvent, const reco::Vertex &pv,
                           const bool hasGoodPV,
			   const TriggerBits triggerBits)
{
  assert(evtInfo);
  
  evtInfo->runNum  = iEvent.id().run();
  evtInfo->lumiSec = iEvent.luminosityBlock();
  evtInfo->evtNum  = iEvent.id().event();

  
  //
  // Pile-up info
  //==============================
  if(!iEvent.isRealData()) {
    edm::Handle< std::vector<PileupSummaryInfo> > hPileupInfoProduct;
    iEvent.getByToken(fPUInfoName_token,hPileupInfoProduct);
    assert(hPileupInfoProduct.isValid());
    const std::vector<PileupSummaryInfo> *inPUInfos = hPileupInfoProduct.product();
    for (std::vector<PileupSummaryInfo>::const_iterator itPUInfo = inPUInfos->begin(); itPUInfo!=inPUInfos->end(); ++itPUInfo) {
      if(itPUInfo->getBunchCrossing()==0) {
        evtInfo->nPU      = itPUInfo->getPU_NumInteractions();
        evtInfo->nPUmean  = itPUInfo->getTrueNumInteractions();
      } else if(itPUInfo->getBunchCrossing()==-1) { 
        evtInfo->nPUm     = itPUInfo->getPU_NumInteractions();
        evtInfo->nPUmeanm = itPUInfo->getTrueNumInteractions();
      } else if(itPUInfo->getBunchCrossing()==1) {
        evtInfo->nPUp     = itPUInfo->getPU_NumInteractions();
        evtInfo->nPUmeanp = itPUInfo->getTrueNumInteractions();
      }
    }
  }

  
  //
  // primary vertex info
  //==============================
  evtInfo->pvx = pv.x();
  evtInfo->pvy = pv.y();
  evtInfo->pvz = pv.z();
  evtInfo->hasGoodPV  = hasGoodPV;
 
  
  //
  // beam spot info
  //==============================
  edm::Handle<reco::BeamSpot> hBeamSpotProduct;
  iEvent.getByToken(fBSName_token,hBeamSpotProduct);
  assert(hBeamSpotProduct.isValid());
  const reco::BeamSpot *bs = hBeamSpotProduct.product();
  evtInfo->bsx = bs->x0();
  evtInfo->bsy = bs->y0();
  evtInfo->bsz = bs->z0();
  	

  //
  // MET filter
  //==============================
  evtInfo->metFilterFailBits=0;
  if(fFillMETFilters) { 
    if (fUseAOD) {
      // beam halo filter using CSCs
      edm::Handle<reco::BeamHaloSummary> hBeamHaloSummary;
      iEvent.getByLabel("BeamHaloSummary",hBeamHaloSummary);
      assert(hBeamHaloSummary.isValid());
      const reco::BeamHaloSummary *beamHaloSummary = hBeamHaloSummary.product();
      if(beamHaloSummary->CSCTightHaloId()) {  // if true, then event has identified beam halo
	evtInfo->metFilterFailBits |= kCSCTightHaloFilter;
      }
      
      // HB,HE anomalous noise filter
      edm::Handle<bool> hHBHENoiseFilterResult;
      iEvent.getByLabel("HBHENoiseFilterResultProducer","HBHENoiseFilterResult",hHBHENoiseFilterResult);
      assert(hHBHENoiseFilterResult.isValid());
      if(!(*hHBHENoiseFilterResult)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kHBHENoiseFilter;
      }
      
      // HCAL laser filter
      edm::Handle<bool> hHCALLaserEventFilter;
      iEvent.getByLabel("hcalLaserEventFilter",hHCALLaserEventFilter);
      assert(hHCALLaserEventFilter.isValid());
      if(!(*hHCALLaserEventFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kHCALLaserEventFilter;
      }
      
      // bad EE SuperCrystal filter
      edm::Handle<bool> hEEBadScFilter;
      iEvent.getByLabel("eeBadScFilter",hEEBadScFilter);
      assert(hEEBadScFilter.isValid());
      if(!(*hEEBadScFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kEEBadScFilter;
      }
      
      // ECAL dead cell filter using trigger primitives
      edm::Handle<bool> hECALDeadCellTriggerPrimitiveFilter;
      iEvent.getByLabel("EcalDeadCellTriggerPrimitiveFilter",hECALDeadCellTriggerPrimitiveFilter);
      assert(hECALDeadCellTriggerPrimitiveFilter.isValid());
      if(!(*hECALDeadCellTriggerPrimitiveFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kECALDeadCellTriggerPrimitiveFilter;
      }
      
      // ECAL bad laser correction filter
      edm::Handle<bool> hECALLaserCorrFilter;
      iEvent.getByLabel("ecalLaserCorrFilter",hECALLaserCorrFilter);
      assert(hECALLaserCorrFilter.isValid());
      if(!(*hECALLaserCorrFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kECALLaserCorrFilter;
      }
      
      // tracking failure filter
      edm::Handle<bool> hTrackingFailureFilter;
      iEvent.getByLabel("trackingFailureFilter",hTrackingFailureFilter);
      assert(hTrackingFailureFilter.isValid());
      if(!(*hTrackingFailureFilter)) {  // if result is "false", then event is flagged as bad
	evtInfo->metFilterFailBits |= kTrackingFailureFilter;
      } 
    }
  }
  //else {//MiniAOD
  //}

  //
  // MET info
  //==============================
  if(fUseAOD) { 

    // PF MET
    edm::Handle<reco::PFMETCollection> hPFMETProduct;
    iEvent.getByToken(fPFMETName_token,hPFMETProduct);
    assert(hPFMETProduct.isValid());
    const reco::PFMET &inPFMET = hPFMETProduct.product()->front();
    evtInfo->pfMET      = inPFMET.pt();
    evtInfo->pfMETphi   = inPFMET.phi();

    // Corrected PF MET
    edm::Handle<reco::PFMETCollection> hPFMETCProduct;
    iEvent.getByToken(fPFMETCName_token,hPFMETCProduct);
    assert(hPFMETCProduct.isValid());
    const reco::PFMET &inPFMETC = hPFMETCProduct.product()->front();
    evtInfo->pfMETC      = inPFMETC.pt();
    evtInfo->pfMETCphi   = inPFMETC.phi();
    
    //  =============== MVA MET ======================
    evtInfo->mvaMET      = 0.0;//inMVAMET.pt();
    evtInfo->mvaMETphi   = 0.0;//inMVAMET.phi();
   
    // ============ Puppi Party ===================
    edm::Handle<reco::PFMETCollection> hPuppET;
    iEvent.getByToken(fPuppETName_token,hPuppET);
    assert(hPuppET.isValid());
    const reco::PFMET &inPuppET = hPuppET.product()->front();
    evtInfo->puppET      = inPuppET.pt();
    evtInfo->puppETphi   = inPuppET.phi();

    // Track MET
    evtInfo->trkMET      = 0.0;
    evtInfo->trkMETphi   = 0.0;
  }
  else { //MiniAOD

    edm::Handle<pat::METCollection> hMETProduct;
    iEvent.getByToken(fPFMETPATName_token,hMETProduct);
    assert(hMETProduct.isValid());
    const pat::MET &inMET = hMETProduct->front();
    // Raw PF MET
    evtInfo->pfMET      = inMET.uncorPt();
    evtInfo->pfMETphi   = inMET.uncorPhi();
    //evtInfo->pfMETCov00 = inMET.getSignificanceMatrix()(0,0);
    //evtInfo->pfMETCov01 = inMET.getSignificanceMatrix()(0,1);
    //evtInfo->pfMETCov11 = inMET.getSignificanceMatrix()(1,1);

    // Corrected PF MET
    evtInfo->pfMETC      = inMET.pt();
    evtInfo->pfMETCphi   = inMET.phi();
    //evtInfo->pfMETCCov00 = inMET.getSignificanceMatrix()(0,0);
    //evtInfo->pfMETCCov01 = inMET.getSignificanceMatrix()(0,1);
    //evtInfo->pfMETCCov11 = inMET.getSignificanceMatrix()(1,1);

    // PUPPI MET
    edm::Handle<pat::METCollection> hPuppET;
    iEvent.getByToken(fPuppETPATName_token,hPuppET);
    assert(hPuppET.isValid());
    const pat::MET &inPuppET = hPuppET.product()->front();
    evtInfo->puppET      = inPuppET.uncorPt();
    evtInfo->puppETphi   = inPuppET.uncorPhi();

    //Type1 PUPPI MET
    //evtInfo->puppETC      = inPuppET.pt();
    //evtInfo->puppETCphi   = inPuppET.phi();
    //evtInfo->puppETCov00  = inPuppET.getSignificanceMatrix()(0,0);
    //evtInfo->puppETCov01  = inPuppET.getSignificanceMatrix()(0,1);
    //evtInfo->puppETCov11 = inPuppET.getSignificanceMatrix()(1,1);

  }
  
  
  //
  // event energy density
  //==============================
  
  // Rho for isolation correction
  edm::Handle<double> hRhoIso;
  edm::InputTag rhoIsoTag(fRhoIsoName,"");
  iEvent.getByToken(rhoIsoTag_token,hRhoIso);
  assert(hRhoIso.isValid());
  evtInfo->rhoIso = *hRhoIso;
  
  // Rho for jet energy correction
  edm::Handle<double> hRhoJet;
  edm::InputTag rhoJetTag(fRhoJetName,"");
  iEvent.getByToken(rhoJetTag_token,hRhoJet);
  assert(hRhoJet.isValid());
  evtInfo->rhoJet = *hRhoJet;


  //
  // fired triggers
  //==============================
  evtInfo->triggerBits = triggerBits;
}


//--------------------------------------------------------------------------------------------------
void FillerEventInfo::computeTrackMET(const reco::Vertex &pv, const reco::PFCandidateCollection *pfCandCol,
                                      float &out_met, float &out_metphi)
{  
  out_met    = 0;
  out_metphi = 0;
  
  double metx=0, mety=0;
  for(reco::PFCandidateCollection::const_iterator itPF = pfCandCol->begin(); itPF!=pfCandCol->end(); ++itPF) {
    if(itPF->trackRef().isNonnull() && pv.trackWeight(itPF->trackRef())>0 && fabs(itPF->vz()-pv.z())<0.2) {
      metx  -= itPF->px();
      mety  -= itPF->py();
     }
  }
  
  TLorentzVector met;
  met.SetPxPyPzE(metx,mety,0,0);
  out_met    = met.Pt();
  out_metphi = met.Phi();
}


void FillerEventInfo::computeTrackMET(const pat::PackedCandidateCollection *pfCandCol,
                                      float &out_met, float &out_metphi)
{
  out_met    = 0;
  out_metphi = 0;

  double metx=0, mety=0;
  for(pat::PackedCandidateCollection::const_iterator itPF = pfCandCol->begin(); itPF!=pfCandCol->end(); ++itPF) {
    // track is:
    // 1) used in the fit of the PV (status=3)
    // 2) not used in fit of any PV but closest in z to the PV (status=2)
    if(itPF->bestTrack()!=0 && itPF->fromPV()>1) {  // (!) MINIAOD: with respect to PV[0]
      metx  -= itPF->px();
      mety  -= itPF->py();
    }
  }

  TLorentzVector met;
  met.SetPxPyPzE(metx,mety,0,0);
  out_met    = met.Pt();
  out_metphi = met.Phi();
}

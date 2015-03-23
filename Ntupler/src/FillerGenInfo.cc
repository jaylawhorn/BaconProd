#include "BaconProd/Ntupler/interface/FillerGenInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include <TClonesArray.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerGenInfo::FillerGenInfo(const edm::ParameterSet &iConfig):
  fGenEvtInfoName(iConfig.getUntrackedParameter<std::string>("edmGenEventInfoName","generator")),
  fLHEEventName  (iConfig.getUntrackedParameter<std::string>("edLHEEventName","source")),
  fGenParName    (iConfig.getUntrackedParameter<std::string>("edmGenParticlesName","genParticles")),
  fFillAll       (iConfig.getUntrackedParameter<bool>("fillAllGen",true))
{}

//--------------------------------------------------------------------------------------------------
FillerGenInfo::~FillerGenInfo(){}

//--------------------------------------------------------------------------------------------------
void FillerGenInfo::fill(TGenEventInfo *genEvtInfo, TClonesArray *array,      
                         const edm::Event &iEvent)
{
  assert(array);
  
  // Get generator event information

  edm::Handle<GenEventInfoProduct> hGenEvtInfoProduct;
  iEvent.getByLabel(fGenEvtInfoName,hGenEvtInfoProduct);
  assert(hGenEvtInfoProduct.isValid());
  /*
  const gen::PdfInfo *pdfInfo = (hGenEvtInfoProduct->hasPDF()) ? hGenEvtInfoProduct->pdf() : 0;
  genEvtInfo->id_1     = (hGenEvtInfoProduct->hasPDF()) ? pdfInfo->id.first    : 0;
  genEvtInfo->id_2     = (hGenEvtInfoProduct->hasPDF()) ? pdfInfo->id.second   : 0;
  genEvtInfo->x_1      = (hGenEvtInfoProduct->hasPDF()) ? pdfInfo->x.first     : 0;
  genEvtInfo->x_2      = (hGenEvtInfoProduct->hasPDF()) ? pdfInfo->x.second    : 0;
  genEvtInfo->weight   = hGenEvtInfoProduct->weight();
  genEvtInfo->scalePDF = hGenEvtInfoProduct->qScale();
  */

  edm::Handle<LHEEventProduct> hLHEEventProduct;
  iEvent.getByLabel(fLHEEventName,hLHEEventProduct);
  assert(hLHEEventProduct.isValid());

  const gen::PdfInfo *pdfInfo = (hLHEEventProduct->pdf()!=0) ? hLHEEventProduct->pdf() : 0;
  genEvtInfo->id_1     = (hLHEEventProduct->pdf()!=0) ? pdfInfo->id.first    : 0;
  genEvtInfo->id_2     = (hLHEEventProduct->pdf()!=0) ? pdfInfo->id.second   : 0;
  genEvtInfo->x_1      = (hLHEEventProduct->pdf()!=0) ? pdfInfo->x.first     : 0;
  genEvtInfo->x_2      = (hLHEEventProduct->pdf()!=0) ? pdfInfo->x.second    : 0;
  genEvtInfo->xPDF_1   = (hLHEEventProduct->pdf()!=0) ? pdfInfo->xPDF.first  : 0;
  genEvtInfo->xPDF_2   = (hLHEEventProduct->pdf()!=0) ? pdfInfo->xPDF.second : 0;
  genEvtInfo->scalePDF = (hLHEEventProduct->pdf()!=0) ? pdfInfo->scalePDF    : 0;

  genEvtInfo->weight   = hGenEvtInfoProduct->weight();

  
  // Get generator particles collection
  edm::Handle<reco::GenParticleCollection> hGenParProduct;
  iEvent.getByLabel(fGenParName,hGenParProduct);
  assert(hGenParProduct.isValid());  
  const reco::GenParticleCollection genParticles = *(hGenParProduct.product());  

  // loop over GEN particles
  std::vector<edm::Ptr<reco::GenParticle>> lMothers;
  TClonesArray &rArray = *array;
  for (reco::GenParticleCollection::const_iterator itGenP = genParticles.begin(); itGenP!=genParticles.end(); ++itGenP) {

    // if not storing all gen particles, then do selective storing
    // Note: assuming Pythia8 status codes
    // Reference: http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
    if(!fFillAll) {
      bool skip=true;

      if(itGenP->status()>20 && itGenP->status()<30)           { skip = false; }  // keep particles from hard scatter process
      if(abs(itGenP->pdgId())>= 5 && abs(itGenP->pdgId())<= 8) { skip = false; }  // keep b, t, b', t'
      if(abs(itGenP->pdgId())>=11 && abs(itGenP->pdgId())<=18) { skip = false; }  // keep leptons
      if(abs(itGenP->pdgId())>=23 && abs(itGenP->pdgId())<=39) { skip = false; }  // keep bosons except photons and gluons
      if(abs(itGenP->pdgId())>10000)                           { skip = false; }  // keep exotic particles

      ///// FSR/ISR photons?

      if(skip) continue;
    }

    // construct object and place in array
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TGenParticle();
    baconhep::TGenParticle *pGenPart = (baconhep::TGenParticle*)rArray[index];
    pGenPart->pdgId  = itGenP->pdgId();
    pGenPart->status = itGenP->status();
    pGenPart->pt     = itGenP->pt();
    pGenPart->eta    = itGenP->eta();
    pGenPart->phi    = itGenP->phi();
    pGenPart->y      = itGenP->rapidity();
    pGenPart->mass   = itGenP->mass();
    if(itGenP->numberOfMothers() >  0 ) {
      int lId = -2;
      edm::Ptr<reco::GenParticle> lMomPtr = edm::refToPtr(itGenP->motherRef()); 
      for(unsigned int im=0; im < lMothers.size(); ++im) { 
	if(lMothers[im] == lMomPtr) {
	  lId = im;  
	  break;
	}
      }
      pGenPart->parent = lId;
    }
    if(itGenP->numberOfMothers() == 0 ) pGenPart->parent = -1;
    edm::Ptr<reco::GenParticle> thePtr(hGenParProduct, itGenP - genParticles.begin());
    lMothers.push_back(thePtr);
  }
}

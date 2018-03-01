#ifndef BACONPROD_NTUPLER_FILLERJET_HH
#define BACONPROD_NTUPLER_FILLERJET_HH

#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "TRandom2.h"
#include <vector>
#include <string>

#include "FWCore/Framework/interface/ConsumesCollector.h"

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
class TClonesArray;
class FactorizedJetCorrector;
class JetCorrectionUncertainty;
namespace trigger {
  class TriggerEvent;
}

namespace baconhep
{
  class FillerJet
  {
    public:
    FillerJet(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC);		
      ~FillerJet();
      
      
      void fill(TClonesArray                     *array,           // output array to be filled
//		TClonesArray                     *iExtraArray,     // Extra Array to be filled
//		TClonesArray                     *iTopArray,       // Top Jet Array to be filled
                const edm::Event                 &iEvent,          // event info
		const edm::EventSetup            &iSetup,          // event setup info
	        const reco::Vertex		 &pv,	           // event primary vertex
		const std::vector<TriggerRecord> &triggerRecords,  // list of trigger names and objects
		const trigger::TriggerEvent      &triggerEvent);   // event trigger objects
            
    protected:
      void initJetCorr(const std::vector<std::string> &jecFiles, 
                       const std::vector<std::string> &jecUncFiles);
      
//      double correction(fastjet::PseudoJet &iJet,double iRho);      
//      void   addJet(TAddJet *pPFJet,const reco::PFJet &itJet,double iRho);
//      void   topJet(TTopJet *pPFJet,const reco::PFJet &itJet,double iRho);
//      fastjet::PseudoJet CACluster(fastjet::PseudoJet &iJet, fastjet::ClusterSequenceArea &iCAClustering); 
      //float              getTau( fastjet::PseudoJet &iJet,int iN, float iKappa );
//      const reco::BasicJet*    match( const reco::PFJet *jet,const reco::BasicJetCollection  *jets );
//      const reco::GenJet*      match( const reco::PFJet *jet,const reco::GenJetCollection    *jets );
      
      // Jet cuts
      double fMinPt;
 
      // Do matching to GenJets?
      bool fUseGen;
      
      // EDM object collection names
      std::string fPVName;
      std::string fRhoName;
      edm::EDGetTokenT<double> rhoTag_token;
      std::string fJetName;
      edm::EDGetTokenT<reco::PFJetCollection> fJetName_token;
      std::string fGenJetName;
      std::string fJetFlavorName;
    //edm::EDGetTokenT<reco::JetFlavourMatchingCollection> fJetFlavorName_token;
      std::string fJetFlavorPhysName;
    //edm::EDGetTokenT<reco::JetFlavourMatchingCollection> fJetFlavorPhysName_token;
//      std::string fPruneJetName;
//      std::string fSubJetName;
      std::string fCSVbtagName;
      edm::EDGetTokenT<reco::JetTagCollection> fCSVbtagName_token;
//      std::string fCSVbtagSubJetName;
//      std::string fJettinessName;
//      std::string fQGLikelihood;
//      std::string fQGLikelihoodSubJets;
      double      fConeSize;
      bool        fComputeFullJetInfo;
      
      // Jet ID MVA
 //     JetPUIDMVACalculator fJetPUIDMVACalc;


//      fastjet::JetDefinition*       fJetDef;
//      fastjet::JetDefinition*       fGenJetDef;
//      fastjet::JetDefinition*       fCAJetDef;
    
//      fastjet::ActiveAreaSpec*      fActiveArea;
//      fastjet::AreaDefinition*      fAreaDefinition;
//      fastjet::ClusterSequenceArea* fClustering;
      
//      fastjet::Pruner* fPruner1;
//      fastjet::Pruner* fPruner2;

//      fastjet::Filter* fFilter1;
//      fastjet::Filter* fFilter2;

//      fastjet::contrib::SoftDropTagger *fSoftDrop1;
//      fastjet::contrib::SoftDropTagger *fSoftDrop2;
//      fastjet::contrib::SoftDropTagger *fSoftDrop3;

//      fastjet::Filter* fTrimmer1;
//      fastjet::Filter* fTrimmer2;
//      fastjet::Filter* fTrimmer3;
//      fastjet::Filter* fTrimmer4;

//      fastjet::CMSTopTagger* fCMSTopTagger;
//      fastjet::HEPTopTagger* fHEPTopTagger;

      // Random number generator for Q-jet volatility
      TRandom2* fRand;

      // function to check if jet passes loose PFJet ID
      bool passPFLooseID();
      
      // JEC corrector
      FactorizedJetCorrector   *fJetCorr;
      JetCorrectionUncertainty *fJetUnc;
  };
}
#endif

import FWCore.ParameterSet.Config as cms

# load jet energy correction parameters
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

#--------------------------------------------------------------------------------
# produce Type 1 MET corrections for PFJets
pfJetMETcorr = cms.EDProducer("PFJetMETcorrInputProducer",
    src = cms.InputTag('ak4PFJets'),
    offsetCorrLabel = cms.string("ak4PFL1Fastjet"),
    jetCorrLabel = cms.string("ak4PFL1FastL2L3"), # NOTE: use "ak4PFL1FastL2L3" for MC / "ak4PFL1FastL2L3Residual" for Data
    jetCorrEtaMax = cms.double(9.9),
    type1JetPtThreshold = cms.double(10.0),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.90),
    skipMuons = cms.bool(True),
    skipMuonSelection = cms.string("isGlobalMuon | isStandAloneMuon")
)					  
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# use MET corrections to produce Type 1 corrected PFMET objects
pfType1CorrectedMet = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag('pfMet'),
    applyType0Corrections = cms.bool(False),
    applyType1Corrections = cms.bool(True),
    srcType1Corrections = cms.VInputTag(
	cms.InputTag('pfJetMETcorr', 'type1')
    ),
    applyType2Corrections = cms.bool(False)
)   
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define sequence to run all modules
producePFMETCorrections = cms.Sequence(
     pfJetMETcorr
   * pfType1CorrectedMet
)
#--------------------------------------------------------------------------------

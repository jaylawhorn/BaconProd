import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.JetCorrectors_cff import *

corrPfMetType1 = cms.EDProducer(
    "PFJetMETcorrInputProducer",
    src = cms.InputTag('ak4PFJets'),
    offsetCorrLabel = cms.InputTag("ak4PFL1FastjetCorrector"),
    jetCorrLabel = cms.InputTag("ak4PFL1FastL2L3Corrector"), # NOTE: use "ak4PFL1FastL2L3Corrector" for MC / "ak4PFL1FastL2L3ResidualCorrector" for Data
    jetCorrEtaMax = cms.double(9.9),
    type1JetPtThreshold = cms.double(10.0),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.90),
    skipMuons = cms.bool(True),
    skipMuonSelection = cms.string("isGlobalMuon | isStandAloneMuon")
)

pfMetT1 = cms.EDProducer(
    "AddCorrectionsToPFMET",
    src = cms.InputTag('pfMet'),
    srcCorrections = cms.VInputTag(
        cms.InputTag('corrPfMetType1', 'type1')
    ),
)   

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define sequence to run all modules
producePFMETCorrections = cms.Sequence(
    ak4PFL1FastL2L3CorrectorChain
    * corrPfMetType1
    * pfMetT1
)
#--------------------------------------------------------------------------------

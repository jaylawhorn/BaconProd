import FWCore.ParameterSet.Config as cms
from JetMETCorrections.Configuration.JetCorrectorsAllAlgos_cff  import *
 
from CondCore.DBCommon.CondDBSetup_cfi import *
puppijec =  cms.ESSource("PoolDBESSource",
                         DBParameters = cms.PSet(messageLevel = cms.untracked.int32(0)),
                         timetype = cms.string('runnumber'),
                         toGet = cms.VPSet(
                           cms.PSet(record  = cms.string('JetCorrectionsRecord'),
                                    tag     = cms.string('JetCorrectorParametersCollection_Summer15_50nsV4_DATA_AK4PFPuppi'),
                                    label   = cms.untracked.string('AK4PFPuppi')
                                    )
                           ),
                         connect = cms.string('sqlite_file:Summer15_50nsV4_DATA.db'),
                         )                                
#Puppi Sequence AK4
ak4PuppiL1FastjetCorrector  = ak4PFCHSL1FastjetCorrector.clone (algorithm   = cms.string('AK4PFPuppi'))
ak4PuppiL2RelativeCorrector = ak4PFCHSL2RelativeCorrector.clone(algorithm   = cms.string('AK4PFPuppi'))
ak4PuppiL3AbsoluteCorrector = ak4PFCHSL3AbsoluteCorrector.clone(algorithm   = cms.string('AK4PFPuppi'))
ak4PuppiResidualCorrector   = ak4PFCHSResidualCorrector.clone  (algorithm   = cms.string('AK4PFPuppi'))
ak4PuppiL1FastL2L3Corrector = ak4PFL1FastL2L3Corrector.clone(
    correctors = cms.VInputTag("ak4PuppiL1FastjetCorrector", "ak4PuppiL2RelativeCorrector", "ak4PuppiL3AbsoluteCorrector")
    )
ak4PuppiL1FastL2L3ResidualCorrector = ak4PFL1FastL2L3Corrector.clone(
    correctors = cms.VInputTag("ak4PuppiL1FastjetCorrector", "ak4PuppiL2RelativeCorrector", "ak4PuppiL3AbsoluteCorrector",'ak4PuppiResidualCorrector')
    )
ak4PuppiL1FastL2L3Chain = cms.Sequence(
    ak4PuppiL1FastjetCorrector * ak4PuppiL2RelativeCorrector * ak4PuppiL3AbsoluteCorrector * ak4PuppiL1FastL2L3Corrector
)
ak4PuppiL1FastL2L3ResidualChain = cms.Sequence(
    ak4PuppiL1FastjetCorrector * ak4PuppiL2RelativeCorrector * ak4PuppiL3AbsoluteCorrector * ak4PuppiResidualCorrector * ak4PuppiL1FastL2L3ResidualCorrector
)

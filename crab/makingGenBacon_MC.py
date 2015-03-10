import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.RecoGenJets_cff     import ak5GenJets

process = cms.Process('MakingGenBacon')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.GlobalTag.globaltag = 'START53_V7G::ALL'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/j/jlawhorn/PYTHIA-GEN-SIM/Zee.root')
)
process.source.inputCommands = cms.untracked.vstring("keep *",
                                                     "drop *_MEtoEDMConverter_*_*")

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(False),
  Rethrow     = cms.untracked.vstring('ProductNotFound'),
  fileMode    = cms.untracked.string('NOMERGE')
)

process.genParticlesForJets = cms.EDProducer("InputGenJetsParticleSelector",
    src = cms.InputTag("genParticles"),
    ignoreParticleIDs = cms.vuint32(
         1000022,
         1000012, 1000014, 1000016,
         2000012, 2000014, 2000016,
         1000039, 5100039,
         4000012, 4000014, 4000016,
         9900012, 9900014, 9900016,
         39,12,14,16),
    partonicFinalState = cms.bool(False),
    excludeResonances = cms.bool(True),
    excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
    tausAsJets = cms.bool(False)
)


process.partons = cms.EDProducer("PartonSelector",
withLeptons = cms.bool(False),
src = cms.InputTag("genParticles")
)

#process.AK4GenJets = ak5GenJets.clone(
#    rParam = cms.double(0.4)
#    )

#process.AK4byRef = cms.EDProducer("JetPartonMatcher",
#jets = cms.InputTag("AK4GenJets"),
#coneSizeToAssociate = cms.double(0.3),
#partons = cms.InputTag("partons")
#)

#process.AK4byValPhys = cms.EDProducer("JetFlavourIdentifier",
#srcByReference = cms.InputTag("AK4byRef"),
#physicsDefinition = cms.bool(True),
#leptonInfo = cms.bool(True)
#)

#process.AK4byValAlgo = cms.EDProducer("JetFlavourIdentifier",
#srcByReference = cms.InputTag("AK4byRef"),
#physicsDefinition = cms.bool(False),
#leptonInfo = cms.bool(True))

#process.AK4jetFlavor = cms.Sequence(AK4byRef*AK4byValPhys*AK4byValAlgo)

process.ntupler = cms.EDAnalyzer('GenNtuplerMod',
  skipOnHLTFail = cms.untracked.bool(False),
  outputName    = cms.untracked.string('/afs/cern.ch/work/j/jlawhorn/PYTHIA-GEN-SIM/Zee_bacon.root'),
  TriggerFile   = cms.untracked.string("tt"),
  edmPVName     = cms.untracked.string('offlinePrimaryVertices'),
  edmPFCandName = cms.untracked.string('particleFlow'),
  
  GenInfo = cms.untracked.PSet(
    isActive            = cms.untracked.bool(True),
    edmGenEventInfoName = cms.untracked.string('generator'),
    edmGenMETName = cms.untracked.string('genMetTrue'),
    edmGenParticlesName = cms.untracked.string('genParticles'),
#    edmGenJetName       = cms.untracked.string('AK4GenJets'),
    fillAllGen          = cms.untracked.bool(True)
  )
)
 
#process.baconSequence = cms.Sequence(process.genParticlesForJets+process.AK4GenJets+process.partons+process.AK4byRef+process.AK4byValPhys+process.AK4byValAlgo+process.ntupler)
process.baconSequence = cms.Sequence(process.partons+process.ntupler)
process.p = cms.Path(process.baconSequence)

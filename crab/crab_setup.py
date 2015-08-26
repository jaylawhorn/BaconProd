from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'data7'
config.General.workArea = 'w_mc2'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'selectWe_cfg.py'
config.JobType.psetName = 'makingBaconPuppiMVAMets_MC.py'
#config.JobType.psetName = 'makingBaconPuppiMVAMets_DATA.py'
config.JobType.outputFiles = ['Output.root']
config.section_("Data")
config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/AODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/AODSIM'
#config.Data.inputDataset = '/ExpressPhysics/Run2015A-Express-v1/FEVT'
#config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/AODSIM'
#config.Data.inputDataset = '/SingleElectron/Run2015B-PromptReco-v1/AOD'
#config.Data.inputDataset = '/SingleMuon/Run2015B-PromptReco-v1/AOD'
#config.Data.inputDataset = '/WWTo2L2Nu_13TeV-powheg/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/AODSIM'
#config.Data.inputDataset = '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/AODSIM'
#config.Data.inputDataset = '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/AODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/AODSIM'
#config.JobType.inputFiles  = ['Summer15_50nsV4_DATA.db']
config.JobType.inputFiles  = ['Summer15_50nsV4_MC.db']
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 7 
config.Data.outLFNDirBase = '/store/user/arapyan/Run2' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = 'Samples'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'

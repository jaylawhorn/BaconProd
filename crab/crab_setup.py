from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'MC_analysis_wm_new_pu20'
config.General.workArea = 'crab_projects7'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'selectWe_cfg.py'
config.JobType.psetName = 'selectWm_cfg.py'
config.JobType.outputFiles = ['output.root']
config.section_("Data")
#config.Data.inputDataset = '/DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU4bx50_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDataset = '/WJetsToLNu_13TeV-madgraph-pythia8-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset = '/WJetsToLNu_13TeV-madgraph-pythia8-tauola/Phys14DR-PU4bx50_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset='/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset='/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_castor_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset='/QCD_Pt-80to170_EMEnriched_Tune4C_13TeV_pythia8/Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset='/QCD_Pt-30to80_EMEnriched_Tune4C_13TeV_pythia8/Phys14DR-PU4bx50_castor_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 
config.Data.outLFN = '/store/user/arapyan/Phys14_final_new/muons/PU20bx25' # or '/store/group/<subdir>'
#config.Data.outLFN = '/store/user/arapyan/phys14/electrons/PU20bx25' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = 'Phys14_samples'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'

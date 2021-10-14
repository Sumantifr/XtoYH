from CRABClient.UserUtilities import config
config = config()

config.section_("General")
config.General.requestName = "tmp"
config.General.requestName       = 'MINIAOD_WH-1j-SMEFTsim-topU3l-v1-MINIAODSIM'
config.General.workArea = 'crab_MINIAOD_WH-1j-SMEFTsim-topU3l-v1-MINIAODSIM' 
config.General.transferLogs = True
config.General.transferOutputs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'miniAOD-prod_PAT.py'
config.JobType.disableAutomaticOutputCollection = False
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 2700
config.JobType.maxMemoryMB = 5000
config.JobType.numCores = 2

config.section_("Data")
config.Data.inputDataset = '/WH-1j-SMEFTsim-topU3l-v1/chatterj-WH-1j-SMEFTsim-topU3l-v1-eaaa4c5d1276b0cb8322b27174d184b2/USER' 
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased' #'EventBased'
config.Data.unitsPerJob = 20
#config.Data.totalUnits  = options.totalUnits
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'WH-1j-SMEFTsim-topU3l-v1-MINIAODSIM'
config.Data.outLFNDirBase = '/store/user/chatterj/'

config.section_("Site")
config.Site.storageSite = 'T2_AT_Vienna'
config.Site.whitelist = ['T2_AT_Vienna']

#config.section_("User")

#config.Data.outputDatasetTag     = 'MINIAOD_' + options.production_label
#config.General.requestName       = 'MINIAOD_WH-1j-SMEFTsim-v1-MINIAODSIM'
#config.Data.outputPrimaryDataset = config.General.requestName # dataset name


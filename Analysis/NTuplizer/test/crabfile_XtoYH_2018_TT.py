from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'crab_XtoYH_2018_TT'
config.General.workArea = 'crab_XtoYH_2018_TT'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'JEC_MC_MINIAOD_cfg.py'
config.JobType.inputFiles = ['Summer19UL18_V5_MC','Summer19UL18_JRV2_MC']
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = ['hist_jerc_l5.root','rootuple_jerc_l5.root']
config.JobType.maxJobRuntimeMin = 2700
#config.JobType.maxMemoryMB = 2480
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.publishDataName = 'XtoYH_2018_TT'
config.Data.publishDBS = 'phys03'

config.Site.storageSite = 'T2_IN_TIFR'
#config.Site.whitelist = ["T2_IN_TIFR"]

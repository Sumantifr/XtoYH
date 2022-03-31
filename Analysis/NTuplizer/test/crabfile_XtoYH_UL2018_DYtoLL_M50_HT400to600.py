from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'crab_XtoYH_UL2018_DYtoLL_M50_HT400to600'
config.General.workArea = 'crab_XtoYH_UL2018_DYtoLL_M50_HT400to600'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'RunJets_MC_MINIAOD_cfg.py'
config.JobType.inputFiles = ['Summer19UL18_V5_MC','Summer19UL18_JRV2_MC','BtagRecommendation106XUL18']
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = ['hist.root','rootuple.root']
config.JobType.maxJobRuntimeMin = 2700
#config.JobType.maxMemoryMB = 2200
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/chatterj/'
config.Data.publication = True
config.Data.outputDatasetTag = 'XtoYH_UL2018_DYtoLL_M50_HT400to600'
config.Data.publishDBS = 'phys03'

config.Site.storageSite = 'T2_IN_TIFR'
#config.Site.whitelist = ["T2_IN_TIFR"]

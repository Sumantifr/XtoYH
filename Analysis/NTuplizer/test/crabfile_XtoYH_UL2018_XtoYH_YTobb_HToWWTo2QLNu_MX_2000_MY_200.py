from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'crab_XtoYH_UL2018_XtoYH_YTobb_HToWWTo2QLNu_MX_2000_MY_200'
config.General.workArea = 'crab_XtoYH_UL2018_XtoYH_YTobb_HToWWTo2QLNu_MX_2000_MY_200'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'RunJets_FastSIM_MINIAOD_cfg.py'
config.JobType.inputFiles = ['Summer19UL18_V5_MC','Summer19UL18_JRV2_MC','BtagRecommendation106XUL18']
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = ['hist.root','rootuple.root']
config.JobType.maxJobRuntimeMin = 2700
#config.JobType.maxMemoryMB = 2200
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2000_MY_200-v1/mukherje-NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2000_MY_200-v1_MINIAODSIM-58e8664f142fdf477807c95cd1ce2e2a/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/chatterj/'
config.Data.publication = True
config.Data.outputDatasetTag = 'XtoYH_UL2018_XtoYH_YTobb_HToWWTo2QLNu_MX_2000_MY_200'
config.Data.publishDBS = 'phys03'

config.Site.storageSite = 'T2_IN_TIFR'
#config.Site.whitelist = ["T2_IN_TIFR"]

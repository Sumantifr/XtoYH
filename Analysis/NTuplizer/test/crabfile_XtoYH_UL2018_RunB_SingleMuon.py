from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'crab_XtoYH_UL2018_RunB_SingleMuon'
config.General.workArea = 'crab_XtoYH_UL2018_RunB_SingleMuon'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'RunJets_Data_2018B_MINIAOD_cfg.py'
config.JobType.inputFiles = ['Summer19UL18_RunB_V5_DATA','Summer19UL18_JRV2_MC','BtagRecommendation106XUL18']
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = ['hist.root','rootuple.root']
config.JobType.maxJobRuntimeMin = 2700
#config.JobType.maxMemoryMB = 2200
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/SingleMuon/Run2018B-UL2018_MiniAODv2-v2/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 5
config.Data.lumiMask ='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
config.Data.outLFNDirBase = '/store/user/chatterj/'
config.Data.publication = True
config.Data.outputDatasetTag = 'XtoYH_UL2018_RunB_SingleMuon'
config.Data.publishDBS = 'phys03'

config.Site.storageSite = 'T2_IN_TIFR'
#config.Site.whitelist = ["T2_IN_TIFR"]

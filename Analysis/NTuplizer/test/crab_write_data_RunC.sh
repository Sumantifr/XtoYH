#!/bin/sh
production_tag=$1
config=$2
Dataset=$3
publication=$4
site=$5
DBS=$6

temp=crabfile_${1}.py

truncate -s 0 $temp

echo "from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'crab_${production_tag}'
config.General.workArea = 'crab_${production_tag}'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '${config}'
config.JobType.inputFiles = ['Summer19UL18_RunC_V5_DATA','Summer19UL18_JRV2_MC','BtagRecommendation106XUL18']
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = ['hist.root','rootuple.root']
config.JobType.maxJobRuntimeMin = 2700
#config.JobType.maxMemoryMB = 2200
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '$Dataset'
config.Data.inputDBS = '$DBS'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 5
config.Data.lumiMask ='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
config.Data.outLFNDirBase = '/store/user/chatterj/'
config.Data.publication = $publication
config.Data.outputDatasetTag = '${production_tag}'
config.Data.publishDBS = 'phys03'

config.Site.storageSite = '$site'
#config.Site.whitelist = [\"T2_IN_TIFR\"]" | cat >>$temp

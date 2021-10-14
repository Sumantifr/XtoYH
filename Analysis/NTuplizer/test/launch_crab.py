from CRABAPI.RawCommand import crabCommand
#from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config

import imp, os, sys
from optparse import OptionParser
import re

cfgPath    = os.path.expandvars( "$CMSSW_BASE/src/Analysis/NTuplizer/test" )
allConfigs = [ x.split(".")[0] for x in os.listdir( cfgPath ) if x.endswith(".py") ]

parser = OptionParser(usage="python launch.py [options] component1 [ component2 ...]", \
                          description="Launch heppy jobs with CRAB3. Components correspond to the variables defined in heppy_samples.py (their name attributes)")
parser.add_option("--production_label", dest="production_label",                                  default="heppy", help="production label")
parser.add_option("--remoteDir",        dest="remoteDir",                                         default="",      help="remote subdirectory")
parser.add_option("--unitsPerJob",      dest="unitsPerJob",      type=int,                        default=1,       help="Nr. of units (files) / crab job")
parser.add_option("--totalUnits",       dest="totalUnits",       type=int,                        default=None,    help="Total nr. of units (files)")
parser.add_option("--config",           dest="config",                     choices = allConfigs,                   help="Which config?")
parser.add_option("--publish",          action='store_true',                                      default=False,   help="Publish on dbs?")
parser.add_option("--dryrun",           action='store_true',                                      default=False,   help="Test script?")
parser.add_option("--Dataset",         dest="Dataset",                                          default=None,    help="Dataset?")
( options, args ) = parser.parse_args()

print "## Starting submission to crab for sample %s ##"%(options.Dataset)

cfgFile      = os.path.join( cfgPath, "%s.py" % options.config )

config = config()

config.section_("General")
config.General.requestName       = options.production_label
config.General.workArea = 'crab_' + options.production_label
config.General.transferLogs = True
config.General.transferOutputs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = cfgFile
#config.JobType.disableAutomaticOutputCollection = False
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 2700
config.JobType.maxMemoryMB = 5000
config.JobType.numCores = 2
config.JobType.inputFiles        = ['Summer19UL18_V5_MC', 'Summer19UL18_JRV2_MC', 'jetToolbox_cff.py']
config.JobType.outputFiles = ['hist.root','rootuple.root']

config.section_("Data")
config.Data.inputDataset = options.Dataset
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = options.unitsPerJob
#config.Data.totalUnits  = options.totalUnits
config.Data.publication = options.publish
config.Data.publishDBS = 'phys03'

config.Data.outLFNDirBase = '/store/user/chatterj/'

config.section_("Site")
config.Site.storageSite = 'T2_AT_Vienna'
#config.Site.whitelist = ['T2_DE_RWTH','T2_US_Nebraska']

#config.section_("User")

if options.dryrun:
    print "Processing %s" %( options.Dataset )
    print "## Dryrun, continue..."
    sys.exit(0)

crabCommand('submit', config=config)

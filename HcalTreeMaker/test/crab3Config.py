from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
#config.General.requestName = 'HcalTupleMaker_RelValTTbar_PU25ns_10_2_0_pre2_MCfull_v03'
config.General.requestName = 'HcalTupleMaker_2018_ZeroBias_DataFull_v6'
#config.General.requestName = 'HcalTupleMaker_2018_HLTPhysics_DataFull_v3'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'run_HcalTupleMaker_2018_MCfull.py'
config.JobType.psetName = 'run_HcalTupleMaker_2018_DataFull.py'
config.JobType.allowUndistributedCMSSW = False
#config.JobType.outputFiles=['trees_relval_ttbar_10_2_0_pre2_2018_MCfull.root']
config.JobType.outputFiles=['trees_2018_ZeroBias_DataFull.root']
#config.JobType.outputFiles=['trees_2018_HLTPhysics_DataFull.root']

config.section_("Data")
# MC example
#config.Data.inputDataset = '/RelValMinBias_13/CMSSW_10_2_0_pre3-101X_upgrade2018_realistic_v7-v1/GEN-SIM-RECO'
#config.Data.inputDataset = '/RelValTTbar_13/CMSSW_10_2_0_pre2-PU25ns_101X_upgrade2018_realistic_v7-v1/GEN-SIM-RECO'
#config.Data.useParent = True
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 1
#KH config.Data.totalUnits = 100
# MC example ends
# Data example
config.Data.inputDataset = '/ZeroBias/Run2018A-v1/RAW'
#config.Data.inputDataset = '/HLTPhysics/Run2018A-v1/RAW'
config.Data.runRange = '315361-315690'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
# Data example ends

config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.publication = False
config.Data.ignoreLocality = False

config.section_("Site")
config.Site.storageSite = 'T3_US_Baylor'
#KH (this whitelisting below is not really necessary. we can use any T2/T3 for running jobs. we can still send output to Baylor)
#KH config.Site.whitelist = ['T3_US_Baylor']

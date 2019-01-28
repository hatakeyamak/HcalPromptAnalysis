from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = ''

config.section_("JobType")
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'QCDForPF_13TeV_TuneCUETP8M1_cfi_GEN_SIM.py' # <== Check!!!
#config.JobType.allowNonProductionCMSSW = False 
config.JobType.allowUndistributedCMSSW = False # Parameter JobType.allowNonProductionCMSSW has been renamed to JobType.allowUndistributedCMSSW

config.section_("Data")

config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 100
#config.Data.unitsPerJob = 25
NJOBS = 1000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS

config.Data.publication = True
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/' # Parameter Data.publishDbsUrl has been renamed to Data.publishDBS
config.Data.outputDatasetTag = 'CMSSW_10_4_0_pre3_v2' # <== Check!!!

config.Data.outLFNDirBase = '/store/user/hatake/crab_outputs'  # Data.outLFN has been renamed to Data.outLFNDirBase
config.Data.outputPrimaryDataset = 'QCDForPF'
config.Data.ignoreLocality = True

config.section_("Site")
config.Site.storageSite = 'T3_US_Baylor'
#KH (this whitelisting below is not really necessary. we can use any T2/T3 for running jobs. we can still send output to Baylor)
#config.Site.whitelist = ['T3_US_Baylor']
config.Site.whitelist = ['T3_US_Baylor','T2_US_*']
config.Site.blacklist = ['T3_US_UCR','T3_US_UMiss']

config.General.transferLogs=True 

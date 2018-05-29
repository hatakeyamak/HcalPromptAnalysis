
## HCAL tuple maker 

### Normal mode
```
cmsRun run_HcalTupleMaker_2018.py
```

### Full(?) content mode (including digis, simhits)
```
cmsRun run_HcalTupleMaker_2018_MCfull.py
## for a short test
cmsRun run_HcalTupleMaker_2018_MCfull.py maxEvents=10 skipEvents=0
```

### Submit jobs to condor (for test) or to CMSConnect

```
#
# prepare for submission to condor
voms-proxy-init -valid 192:0 -voms cms
mkdir -p /cms/data/store/user/${USER}/condor/tarballs
mkdir -p /cms/data/store/user/${USER}/condor/outputs
mkdir -p log
../make_tarball.sh
# also check condor.jdl, cmsRun.csh, run_HcalTupleMaker_2018_MCfull.py

#
# submit to condor or CMSconnect
# for local condor on n147
condor_submit submit.jdl
#
# or from another sessoon where cmsenv is not run yet, set up CMSconnect tool
# see: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CMSConnectHATSatLPC2017
# in particular source <some path>/cmsconnect_client.sh
# then
connect submit submit.jdl
# 
```

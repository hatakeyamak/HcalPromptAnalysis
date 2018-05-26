#!/bin/bash

export MYCMSSW=`echo $CMSSW_BASE | rev | cut -d'/' -f-1 | rev`

cd ${CMSSW_BASE}/..
tar -zcvf /cms/data/store/user/${USER}/condor/${MYCMSSW}_condor.tgz ${MYCMSSW} --exclude=tmp --exclude=test --exclude='*.root'
cd -


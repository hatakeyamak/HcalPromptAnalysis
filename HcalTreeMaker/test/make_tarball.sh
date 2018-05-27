#!/bin/bash

export MYCMSSW=`echo $CMSSW_BASE | rev | cut -d'/' -f-1 | rev`

cd ${CMSSW_BASE}/..
mkdir -p /cms/data/store/user/${USER}/condor/tarballs
tar -zcvf /cms/data/store/user/${USER}/condor/tarballs/${MYCMSSW}_condor.tgz ${MYCMSSW} --exclude=tmp --exclude=test --exclude='*.root'
cd -


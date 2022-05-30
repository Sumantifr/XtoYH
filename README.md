# XtoYH

Framework for ntuple production in X->Y(->bb)H(->WW) analysis

- Log in to your lxplus account

- cd work/private

- mkdir XtoYH

- cd XtoYH

- cmsrel CMSSW_10_6_29 <br/>
#For CMSSW_10_6_X (X<=25), some changes need to be made for ECAL prefiring weights

- cd CMSSW_10_6_29/src

- git clone https://github.com/Sumantifr/XtoYH.git . <br/>
  *(Don't forget '.')*

- scram b -j10 

## For a test run: 

cd CMSSSW_10_6_29/src

- cmsenv

- cd $CMSSW_BASE/src/Analysis/NTuplizer/test/

- voms-proxy-init -rfc -voms cms -valid 48:00

- cmsRun PNet_MC_MINIAOD_cfg.py

Enjoy!

## For submitting crab jobs (MC):

- cd CMSSSW_10_6_29/src

- cmsenv

- cd $CMSSW_BASE/src/Analysis/NTuplizer/test/

- voms-proxy-init -rfc -voms cms -valid 48:00

- ./generate_crab_miniaod_UL18.sh <br/>
  *This will create all the files necessary to submit jobs, but the command will not submit jobs!!*

./crab_submit.sh

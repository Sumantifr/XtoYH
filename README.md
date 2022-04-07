# XtoYH

Log in to your lxplus account

cd work/private

mkdir XtoYH

cd XtoYH

cmsrel CMSSW_10_6_29 

#For CMSSW_10_6_X (X<=25), some changes need to be made for ECAL prefiring weights

cd CMSSW_10_6_29/src

git clone https://github.com/Sumantifr/XtoYH.git .

scram b -j10

Enjoy! 

For a test run: 
==============
cd CMSSSW_10_6_29/src

cmsenv

cd $CMSSW_BASE/src/Analysis/NTuplizer/test/

voms-proxy-init -rfc -voms cms -valid 48:00

cmsRun PNet_MC_MINIAOD_cfg.py

For submitting crab jobs (MC):
=========================

cd CMSSSW_10_6_29/src

cmsenv

cd $CMSSW_BASE/src/Analysis/NTuplizer/test/

voms-proxy-init -rfc -voms cms -valid 48:00

./generate_crab_miniaod_UL18.sh 
  (This will create all the files necessary to submit jobs, but the command will not submit jobs!!)

./crab_submit.sh

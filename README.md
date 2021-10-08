# XtoYH

Log in to your lxplus account

cd work/private
mkdir XtoYH
cd XtoYH
cmsrel CMSSSW_10_6_26  
#For CMSSW_10_6_X (X<=25), some changes need to be made for ECAL prefiring weights
cd CMSSW_10_6_26/src
git clone https://github.com/Sumantifr/XtoYH.git .
scram b -j9
Enjoy! 

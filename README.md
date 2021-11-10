# XtoYH

Log in to your lxplus account

cd work/private

mkdir XtoYH

cd XtoYH

cmsrel CMSSSW_10_6_26  

#For CMSSW_ version should not matter much here

cd CMSSW_10_6_26/src

git clone https://github.com/Sumantifr/XtoYH.git -b Offline_Analysis

Compile the code with:
./Makefile Anal_XtoYH

Run the executable:
./Anal_XtoYH.exe

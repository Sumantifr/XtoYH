#!/bin/bash

config=RunJets_FastSIM_MINIAOD_cfg.py
publish=True
site=T2_IN_TIFR
DBS=phys03

sample_names=( 
XtoYH_YTobb_HToWWTo2QLNu_MX_1500_MY_200
XtoYH_YTobb_HToWWTo2QLNu_MX_2000_MY_200
)
sample_data=(
/NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_1500_MY_200-v1/mukherje-NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_1500_MY_200-v1_MINIAODSIM-58e8664f142fdf477807c95cd1ce2e2a/USER
/NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2000_MY_200-v1/mukherje-NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2000_MY_200-v1_MINIAODSIM-58e8664f142fdf477807c95cd1ce2e2a/USER
)

nsamples=${#sample_data[*]}
if [ $nsamples != ${#sample_names[*]} ]; 
then
	echo "No of names & samples are not same!! please check! (samples $nsamples names ${#sample_names[*]}"
	exit
fi

fil_list=crab_submit_fastsim
mon_list=crab_monitor_fastsim
truncate -s 0 ${fil_list}.sh
echo "#!/bin/bash" | cat >>${fil_list}.sh
truncate -s 0 ${mon_list}.sh

i=1
while [[ $i -le $nsamples ]]
do
	echo ${sample_data[i-1]} ${sample_names[i-1]}
	label=XtoYH_UL2018_${sample_names[i-1]}
	./crab_write.sh $label $config  ${sample_data[i-1]} $publish $site $DBS
	echo "crab submit -c crabfile_${label}.py" | cat >>${fil_list}.sh
	echo "crab status -d crab_${label}/crab_crab_${label}/" | cat >>${mon_list}.sh
	((i = i + 1))
done

chmod 744 ${fil_list}.sh
chmod 744 ${mon_list}.sh

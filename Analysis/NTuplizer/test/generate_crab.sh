#!/bin/bash

config=PNet_MC_MINIAOD_cfg.py
publish=False
site=T2_IN_TIFR

sample_names=( 
TT_Dilep TT_Had TT_Semilep
)
#ST_tch_top ST_tch_atop ST_sch_thad ST_sch_tlep ST_tW_top ST_tW_atop 
#Wj)
sample_data=(
/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM
/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM
/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM
)

nsamples=${#sample_data[*]}
if [ $nsamples != ${#sample_names[*]} ]; 
then
	echo "No of names & samples are not same!! please check! (samples $nsamples names ${#sample_names[*]}"
	exit
fi

fil_list=crab_submit
mon_list=crab_monitor
truncate -s 0 ${fil_list}.sh
echo "#!/bin/bash" | cat >>${fil_list}.sh
truncate -s 0 ${mon_list}.sh

i=1
while [[ $i -le $nsamples ]]
do
	echo ${sample_data[i-1]} ${sample_names[i-1]}
	label=XtoYH_2018_${sample_names[i-1]}
	./crab_write.sh $label $config  ${sample_data[i-1]} $publish $site
	#echo "crab submit -c crabfile_${label}.py" | cat >>${fil_list}.sh
	echo "crab status -d crab_${label}_${sample_names[i-1]}/crab_crab_${label}_${sample_names[i-1]}/" | cat >>${mon_list}.sh
	((i = i + 1))
done

chmod 744 ${fil_list}.sh
chmod 744 ${mon_list}.sh

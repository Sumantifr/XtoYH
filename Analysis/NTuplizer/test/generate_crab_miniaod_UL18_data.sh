#!/bin/bash

config=RunJets_MC_MINIAOD_cfg.py
publish=True
site=T2_IN_TIFR
DBS=global

sample_names=(
RunA_SingleMuon
RunB_SingleMuon
RunC_SingleMuon
RunD_SingleMuon
RunA_EGamma
RunB_EGamma
RunC_EGamma
RunD_EGamma
)

sample_data=(
/SingleMuon/Run2018A-UL2018_MiniAODv2-v3/MINIAOD
/SingleMuon/Run2018B-UL2018_MiniAODv2-v2/MINIAOD
/SingleMuon/Run2018C-UL2018_MiniAODv2-v2/MINIAOD
/SingleMuon/Run2018D-UL2018_MiniAODv2-v3/MINIAOD
/EGamma/Run2018A-UL2018_MiniAODv2-v1/MINIAOD
/EGamma/Run2018B-UL2018_MiniAODv2-v1/MINIAOD
/EGamma/Run2018C-UL2018_MiniAODv2-v1/MINIAOD
/EGamma/Run2018D-UL2018_MiniAODv2-v2/MINIAOD
)

nsamples=${#sample_data[*]}
if [ $nsamples != ${#sample_names[*]} ]; 
then
	echo "No of names & samples are not same!! please check! (samples $nsamples names ${#sample_names[*]}"
	exit
fi

fil_list=crab_submit_2018_data
mon_list=crab_monitor_2018_data
truncate -s 0 ${fil_list}.sh
echo "#!/bin/bash" | cat >>${fil_list}.sh
truncate -s 0 ${mon_list}.sh

i=1
while [[ $i -le $nsamples ]]
do
	echo ${sample_data[i-1]} ${sample_names[i-1]}
	label=XtoYH_UL2018_${sample_names[i-1]}
	runtag=`expr $i % 4`
	echo "runtag" $runtag
	if [ $runtag = 1 ];
	then
		./crab_write_data_RunA.sh $label RunJets_Data_2018A_MINIAOD_cfg.py ${sample_data[i-1]} $publish $site $DBS
	fi
	if [ $runtag = 2 ];
	then
		./crab_write_data_RunB.sh $label RunJets_Data_2018B_MINIAOD_cfg.py ${sample_data[i-1]} $publish $site $DBS
	fi
	if [ $runtag = 3 ];
	then
		./crab_write_data_RunC.sh $label RunJets_Data_2018C_MINIAOD_cfg.py ${sample_data[i-1]} $publish $site $DBS
	fi
	if [ $runtag = 0 ];
	then
		./crab_write_data_RunD.sh $label RunJets_Data_2018D_MINIAOD_cfg.py ${sample_data[i-1]} $publish $site $DBS
	fi
	#./crab_write.sh $label $config  ${sample_data[i-1]} $publish $site $DBS
	echo "crab submit -c crabfile_${label}.py" | cat >>${fil_list}.sh
	echo "crab status -d crab_${label}/crab_crab_${label}/" | cat >>${mon_list}.sh
	((i = i + 1))
done

chmod 744 ${fil_list}.sh
chmod 744 ${mon_list}.sh

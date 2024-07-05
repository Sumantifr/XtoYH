#!/bin/bash

config=RunJets_MC_MINIAOD_cfg.py
publish=True
site=T2_AT_Vienna
DBS=global

sample_names=(
RunF_SingleMuon
RunG_SingleMuon
RunH_SingleMuon
RunF_DoubleMuon
RunG_DoubleMuon
RunH_DoubleMuon
RunF_DoubleEG
RunG_DoubleEG
RunH_DoubleEG
RunF_SingleElectron
RunG_SingleElectron
RunH_SingleElectron
RunF_JetHT
RunG_JetHT
RunH_JetHT
RunF_MuonEG
RunG_MuonEG
RunH_MuonEG
)

sample_data=(
/SingleMuon/Run2016F-UL2016_MiniAODv2-v2/MINIAOD
/SingleMuon/Run2016G-UL2016_MiniAODv2-v2/MINIAOD
/SingleMuon/Run2016H-UL2016_MiniAODv2-v2/MINIAOD
/DoubleMuon/Run2016F-UL2016_MiniAODv2-v1/MINIAOD
/DoubleMuon/Run2016G-UL2016_MiniAODv2-v1/MINIAOD
/DoubleMuon/Run2016H-UL2016_MiniAODv2-v2/MINIAOD
/DoubleEG/Run2016F-UL2016_MiniAODv2-v1/MINIAOD
/DoubleEG/Run2016G-UL2016_MiniAODv2-v1/MINIAOD
/DoubleEG/Run2016H-UL2016_MiniAODv2-v1/MINIAOD
/SingleElectron/Run2016F-UL2016_MiniAODv2-v2/MINIAOD
/SingleElectron/Run2016G-UL2016_MiniAODv2-v2/MINIAOD
/SingleElectron/Run2016H-UL2016_MiniAODv2-v2/MINIAOD
/JetHT/Run2016F-UL2016_MiniAODv2-v2/MINIAOD
/JetHT/Run2016G-UL2016_MiniAODv2-v2/MINIAOD
/JetHT/Run2016H-UL2016_MiniAODv2-v2/MINIAOD
/MuonEG/Run2016F-UL2016_MiniAODv2-v2/MINIAOD
/MuonEG/Run2016G-UL2016_MiniAODv2-v2/MINIAOD
/MuonEG/Run2016H-UL2016_MiniAODv2-v2/MINIAOD
)

nsamples=${#sample_data[*]}
if [ $nsamples != ${#sample_names[*]} ]; 
then
	echo "No of names & samples are not same!! please check! (samples $nsamples names ${#sample_names[*]}"
	exit
fi

fil_list=crab_submit_2016_postVFP_data
mon_list=crab_monitor_2016_postVFP_data
truncate -s 0 ${fil_list}.sh
echo "#!/bin/bash" | cat >>${fil_list}.sh
truncate -s 0 ${mon_list}.sh

i=1
while [[ $i -le $nsamples ]]
do
	echo ${sample_data[i-1]} ${sample_names[i-1]}
	label=XtoYH_UL2016postVFP_${sample_names[i-1]}
	runtag=`expr $i % 3`
	echo "runtag" $runtag
        if [ $runtag = 1 ];
        then
                ./crab_write_data_2016postVFP.sh $label RunJets_Data_2016postVFP_MINIAOD_cfg.py ${sample_data[i-1]} $publish $site $DBS
        fi
        if [ $runtag = 2 ];
        then
                ./crab_write_data_2016postVFP.sh $label RunJets_Data_2016postVFP_MINIAOD_cfg.py ${sample_data[i-1]} $publish $site $DBS
        fi
        if [ $runtag = 0 ];
        then
                ./crab_write_data_2016postVFP.sh $label RunJets_Data_2016postVFP_MINIAOD_cfg.py ${sample_data[i-1]} $publish $site $DBS
        fi
	echo "crab submit -c crabfile_${label}.py" | cat >>${fil_list}.sh
	echo "crab status -d crab_${label}/crab_crab_${label}/" | cat >>${mon_list}.sh
	((i = i + 1))
done

chmod 744 ${fil_list}.sh
chmod 744 ${mon_list}.sh

#!/bin/bash

config=RunJets_MC_MINIAOD_cfg.py
publish=True
site=T2_IN_TIFR
DBS=global

sample_names=(
RunB_SingleMuon
RunC_SingleMuon
RunD_SingleMuon
RunE_SingleMuon
RunF_SingleMuon
RunB_DoubleMuon
RunC_DoubleMuon
RunD_DoubleMuon
RunE_DoubleMuon
RunF_DoubleMuon
RunB_DoubleEG
RunC_DoubleEG
RunD_DoubleEG
RunE_DoubleEG
RunF_DoubleEG
RunB_SingleElectron
RunC_SingleElectron
RunD_SingleElectron
RunE_SingleElectron
RunF_SingleElectron
RunB_JetHT
RunC_JetHT
RunD_JetHT
RunE_JetHT
RunF_JetHT
RunB_MuonEG
RunC_MuonEG
RunD_MuonEG
RunE_MuonEG
RunF_MuonEG
)

sample_data=(
/SingleMuon/Run2017B-UL2017_MiniAODv2_GT36-v2/MINIAOD
/SingleMuon/Run2017C-UL2017_MiniAODv2_GT36-v2/MINIAOD
/SingleMuon/Run2017D-UL2017_MiniAODv2_GT36-v2/MINIAOD
/SingleMuon/Run2017E-UL2017_MiniAODv2_GT36-v2/MINIAOD
/SingleMuon/Run2017F-UL2017_MiniAODv2_GT36-v2/MINIAOD
/DoubleMuon/Run2017B-UL2017_MiniAODv2-v1/MINIAOD
/DoubleMuon/Run2017C-UL2017_MiniAODv2-v1/MINIAOD
/DoubleMuon/Run2017D-UL2017_MiniAODv2-v1/MINIAOD
/DoubleMuon/Run2017E-UL2017_MiniAODv2-v2/MINIAOD
/DoubleMuon/Run2017F-UL2017_MiniAODv2-v1/MINIAOD
/DoubleEG/Run2017B-UL2017_MiniAODv2-v1/MINIAOD
/DoubleEG/Run2017C-UL2017_MiniAODv2-v2/MINIAOD
/DoubleEG/Run2017D-UL2017_MiniAODv2-v1/MINIAOD
/DoubleEG/Run2017E-UL2017_MiniAODv2-v1/MINIAOD
/DoubleEG/Run2017F-UL2017_MiniAODv2-v2/MINIAOD
/SingleElectron/Run2017B-UL2017_MiniAODv2-v1/MINIAOD
/SingleElectron/Run2017C-UL2017_MiniAODv2-v1/MINIAOD
/SingleElectron/Run2017D-UL2017_MiniAODv2-v1/MINIAOD
/SingleElectron/Run2017E-UL2017_MiniAODv2-v1/MINIAOD
/SingleElectron/Run2017F-UL2017_MiniAODv2-v1/MINIAOD
/JetHT/Run2017B-UL2017_MiniAODv2-v1/MINIAOD
/JetHT/Run2017C-UL2017_MiniAODv2-v1/MINIAOD
/JetHT/Run2017D-UL2017_MiniAODv2-v1/MINIAOD
/JetHT/Run2017E-UL2017_MiniAODv2-v1/MINIAOD
/JetHT/Run2017F-UL2017_MiniAODv2-v1/MINIAOD
/MuonEG/Run2017B-UL2017_MiniAODv2-v1/MINIAOD
/MuonEG/Run2017C-UL2017_MiniAODv2-v1/MINIAOD
/MuonEG/Run2017D-UL2017_MiniAODv2-v1/MINIAOD
/MuonEG/Run2017E-UL2017_MiniAODv2-v1/MINIAOD
/MuonEG/Run2017F-UL2017_MiniAODv2-v1/MINIAOD
)

nsamples=${#sample_data[*]}
if [ $nsamples != ${#sample_names[*]} ]; 
then
	echo "No of names & samples are not same!! please check! (samples $nsamples names ${#sample_names[*]}"
	exit
fi

fil_list=crab_submit_2017_data
mon_list=crab_monitor_2017_data
truncate -s 0 ${fil_list}.sh
echo "#!/bin/bash" | cat >>${fil_list}.sh
truncate -s 0 ${mon_list}.sh

i=1
while [[ $i -le $nsamples ]]
do
	echo ${sample_data[i-1]} ${sample_names[i-1]}
	label=XtoYH_UL2017_${sample_names[i-1]}
	runtag=`expr $i % 5`
	echo "runtag" $runtag
        if [ $runtag = 1 ];
        then
                ./crab_write_data_2017B.sh $label RunJets_Data_2017B_MINIAOD_cfg.py ${sample_data[i-1]} $publish $site $DBS
        fi
        if [ $runtag = 2 ];
        then
                ./crab_write_data_2017C.sh $label RunJets_Data_2017C_MINIAOD_cfg.py ${sample_data[i-1]} $publish $site $DBS
        fi
        if [ $runtag = 3 ];
        then
                ./crab_write_data_2017D.sh $label RunJets_Data_2017D_MINIAOD_cfg.py ${sample_data[i-1]} $publish $site $DBS
        fi
        if [ $runtag = 4 ];
        then
                ./crab_write_data_2017E.sh $label RunJets_Data_2017E_MINIAOD_cfg.py ${sample_data[i-1]} $publish $site $DBS
        fi
        if [ $runtag = 0 ];
        then
                ./crab_write_data_2017F.sh $label RunJets_Data_2017F_MINIAOD_cfg.py ${sample_data[i-1]} $publish $site $DBS
        fi
	echo "crab submit -c crabfile_${label}.py" | cat >>${fil_list}.sh
	echo "crab status -d crab_${label}/crab_crab_${label}/" | cat >>${mon_list}.sh
	((i = i + 1))
done

chmod 744 ${fil_list}.sh
chmod 744 ${mon_list}.sh

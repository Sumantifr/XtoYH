#include<stdio.h>
#include <bits/stdc++.h>
using namespace std;

void jobfilemaker(bool isDATA=0, bool isDL=0)
{

string Filenames_MC[] = {
"DYJetsToLL_M-10to50_XtoYH",
"DYJetsToLL_M-50_HT-100To200_XtoYH",
"DYJetsToLL_M-50_HT-1200To2500_XtoYH",
"DYJetsToLL_M-50_HT-200To400_XtoYH",
"DYJetsToLL_M-50_HT-2500ToInf_XtoYH",
"DYJetsToLL_M-50_HT-400To600_XtoYH",
"DYJetsToLL_M-50_HT-600To800_XtoYH",
"DYJetsToLL_M-50_HT-70to100_XtoYH",
"DYJetsToLL_M-50_HT-800To1200_XtoYH",
"ST_s-channel_XtoYH",
"ST_t-channel_antitop_XtoYH",
"ST_t-channel_top_XtoYH",
"ST_tW_antitop_XtoYH",
"ST_tW_top_XtoYH",
"TTTo2L2Nu_XtoYH",
"TTToHadronic_XtoYH",
"TTToSemiLeptonic_XtoYH",
"WJetsToLNu_HT-100To200_XtoYH",
"WJetsToLNu_HT-1200To2500_XtoYH",
"WJetsToLNu_HT-200To400_XtoYH",
"WJetsToLNu_HT-2500ToInf_XtoYH",
"WJetsToLNu_HT-400To600_XtoYH",
"WJetsToLNu_HT-600To800_XtoYH",
"WJetsToLNu_HT-70To100_XtoYH",
"WJetsToLNu_HT-800To1200_XtoYH",
"WWTo1L1Nu2Q_XtoYH",
"WWTo2L2Nu_XtoYH",
"WWTo4Q_XtoYH",
"WZTo1L1Nu2Q_XtoYH",
"WZTo2Q2L_XtoYH",
"WZTo2Q2Nu_XtoYH",
"WZTo3LNu_XtoYH",
"ZZTo2L2Nu_XtoYH",
"ZZTo2Q2L_XtoYH",
"ZZTo2Q2Nu_XtoYH",
"ZZTo4L_XtoYH",
"ZZTo4Q_XtoYH",
"QCD_HT300to500_XtoYH",
"QCD_HT500to700_XtoYH",
"QCD_HT1000to1500_XtoYH",
"QCD_HT1500to2000_XtoYH",
"QCD_HT2000toInf_XtoYH",
"QCD_HT700to1000_XtoYH"
};

string Filenames_Data_DL[] = {
"EGamma_UL2018A_XtoYH_Nov_2021",
"EGamma_UL2018B_XtoYH_Nov_2021",
"EGamma_UL2018C_XtoYH_Nov_2021",
"EGamma_UL2018D_XtoYH_Nov_2021",
"DoubleMuon_UL2018A_XtoYH_Nov_2021",
"DoubleMuon_UL2018B_XtoYH_Nov_2021",
"DoubleMuon_UL2018C_XtoYH_Nov_2021",
"DoubleMuon_UL2018D_XtoYH_Nov_2021",
"MuonEG_UL2018A_XtoYH_Nov_2021",
"MuonEG_UL2018B_XtoYH_Nov_2021",
"MuonEG_UL2018C_XtoYH_Nov_2021",
"MuonEG_UL2018D_XtoYH_Nov_2021",
};

string Filenames_Data_SL[] = {
"EGamma_UL2018A_XtoYH_Nov_2021",
"EGamma_UL2018B_XtoYH_Nov_2021",
"EGamma_UL2018C_XtoYH_Nov_2021",
"EGamma_UL2018D_XtoYH_Nov_2021",
"SingleMuon_UL2018A_XtoYH_Nov_2021",
"SingleMuon_UL2018B_XtoYH_Nov_2021",
"SingleMuon_UL2018C_XtoYH_Nov_2021",
"SingleMuon_UL2018D_XtoYH_Nov_2021",
"JetHT_UL2018A_XtoYH_Nov_2021",
"JetHT_UL2018B_XtoYH_Nov_2021",
"JetHT_UL2018C_XtoYH_Nov_2021",
"JetHT_UL2018D_XtoYH_Nov_2021",
};


string path = "/afs/cern.ch/user/c/chatterj/work/private/XToYH/CMSSW_10_6_27/src/XtoYH/HistoMaker/";
string proxy = "/tmp/x509up_u81649";

fstream file_sub;
char name_submit[100];
if(isDATA){
sprintf(name_submit,"condor_submit_data.sh");
}
else{
sprintf(name_submit,"condor_submit.sh");
}
file_sub.open(name_submit,ios::out);

int nfile_data_DL = sizeof(Filenames_Data_DL)/sizeof(Filenames_Data_DL[0]);
int nfile_data_SL = sizeof(Filenames_Data_SL)/sizeof(Filenames_Data_SL[0]);
int nfile_mc = sizeof(Filenames_MC)/sizeof(Filenames_MC[0]);

vector<string> files;
if(isDATA){
if(isDL){
for(int fg=0; fg<nfile_data_DL; fg++){
files.push_back((Filenames_Data_DL[fg]));
}
}
else{
for(int fg=0; fg<nfile_data_SL; fg++){
files.push_back((Filenames_Data_SL[fg]));
}
}
}
else{
for(int fg=0; fg<nfile_mc; fg++){
files.push_back((Filenames_MC[fg]));
}
}

int nfile = int(files.size());

for (int ii=0;ii<nfile;ii++)
{

fstream file;
char name_buffer[512];
sprintf(name_buffer,"execute_%s.csh",files[ii].c_str());
file.open(name_buffer,ios::out);
if(!file)
   {
       cout<<"Error in creating file!!!";
       //return 0;
   }

   file <<"#!/bin/bash\n";
   file << "source /cvmfs/cms.cern.ch/cmsset_default.sh\n";
   file <<"cd "<<path<<" \n";
   file<<"export X509_USER_PROXY="<<proxy<<"\n";
   file<<"eval `scramv1 runtime -sh`\n";
   file<<"./histomaker_comb.exe  "<<isDL<<" "<<isDATA<<" "<<files[ii] << ".root 0" ;

   cout<<"execute created successfully." << endl;
   file.close();

fstream file1;
char name_buffer1[512];
sprintf(name_buffer1,"submit_%s.sh",files[ii].c_str());
file1.open(name_buffer1,ios::out);
   if(!file1)
   {
       cout<<"Error in creating file!!!";
       //return 0;
   }
file1 << "universe = vanilla\n";
file1 << "executable = execute_"<< files[ii] << ".csh\n";
file1 << "getenv = TRUE\n";
file1 << "log =" << files[ii] << ".log\n";
file1 << "output ="<< files[ii] <<".out\n";
file1 << "error = "<< files[ii] <<".error\n";
file1 << "notification = never\n";
file1 << "should_transfer_files = YES\n";
file1 << "when_to_transfer_output = ON_EXIT\n";
file1 << "queue\n";
//file1 << "+JobFlavour = \"workday\"\n";
file1<< "+MaxRuntime = 30000\n";
cout<<"sub created successfully." << endl;
file1.close();

file_sub<<"condor_submit "<<name_buffer1<<"\n";
}

file_sub.close();
  
//return 0;
}

#include<stdio.h>
#include <bits/stdc++.h>
using namespace std;

void jobfilemaker_dilep()
{

TString Filenames[] = {
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


int nfile = sizeof(Filenames)/sizeof(Filenames[0]);
for (int ii=0;ii<nfile;ii++)
{

fstream file;
char name_buffer[512];
sprintf(name_buffer,"execute_"+Filenames[ii]+".csh");
file.open(name_buffer,ios::out);
if(!file)
   {
       cout<<"Error in creating file!!!";
       return 0;
   }

   file << "#!/bin/bash\nsource /cvmfs/cms.cern.ch/cmsset_default.sh\n cd /afs/cern.ch/work/m/mukherje/public/For_Suman_da/HISTMAKER/CMSSW_10_2_10/src/ \nexport X509_USER_PROXY=/afs/cern.ch/work/m/mukherje/public/For_Suman_da/HISTMAKER/CMSSW_10_2_10/src/x509up_u56596\neval `scramv1 runtime -sh`\n./histomaker_comb.exe  1 0 "<< Filenames[ii] << ".root 0" ;
cout<<"execute created successfully." << endl;
file.close();

fstream file1;
char name_buffer1[512];
sprintf(name_buffer1,"submit_"+Filenames[ii]+".sh");
file1.open(name_buffer1,ios::out);
   if(!file1)
   {
       cout<<"Error in creating file!!!";
       return 0;
   }
file1 << "universe = vanilla\nexecutable = execute_"<< Filenames[ii] << ".csh\ngetenv = TRUE\nlog =" << Filenames[ii] << ".log\noutput ="<< Filenames[ii] <<".out\nerror = "<< Filenames[ii] <<".error\nnotification = never\nshould_transfer_files = YES\nwhen_to_transfer_output = ON_EXIT\nqueue";
cout<<"sub created successfully." << endl;
file1.close();
}
  return 0;
}

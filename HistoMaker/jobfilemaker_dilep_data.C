#include<stdio.h>
#include <bits/stdc++.h>
using namespace std;

void jobfilemaker_dilep_data()
{

TString Filenames[] = {
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


   file << "#!/bin/bash\nsource /cvmfs/cms.cern.ch/cmsset_default.sh\n cd /afs/cern.ch/work/m/mukherje/public/For_Suman_da/HISTMAKER/CMSSW_10_2_10/src/ \nexport X509_USER_PROXY=/afs/cern.ch/work/m/mukherje/public/For_Suman_da/HISTMAKER/CMSSW_10_2_10/src/x509up_u56596\neval `scramv1 runtime -sh`\n./histomaker_comb.exe  1 1 "<< Filenames[ii] << ".root 0" ;

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

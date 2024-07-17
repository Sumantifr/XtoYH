#include<stdio.h>
#include <bits/stdc++.h>
using namespace std;

void signaljobmaker(bool isDL=0)
{

const int nymass = 19;

int xymass_DL[][1+nymass]={
{500,  100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{600,  100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{650,  100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{700,  100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{800,  100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{900,  100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{1000, 100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{1200, 100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{1400, 100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{1600, 100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{1800, 100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{2000, 100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{2200, 100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{2400, 100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{2500, 100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{2600, 100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{2800, 100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{3000, 100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{3500, 100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000},
{4000, 100, 125, 150, 190, 250, 300, 350, 400, 450, 500, 600, 60, 700, 70, 800, 80, 900, 90, 1000}
};

int nxmass = sizeof(xymass_DL)/sizeof(xymass_DL[0]);

string path = "/afs/cern.ch/user/c/chatterj/work/private/XToYH/CMSSW_10_6_27/src/XtoYH/HistoMaker/";
string proxy = "/tmp/x509up_u81649";

fstream file_sub;
fstream file_local;
char name_submit[100], name_local[100];
sprintf(name_submit,"condor_submit_signal.sh");
sprintf(name_local,"run_local_signal.sh");

file_sub.open(name_submit,ios::out);
file_local.open(name_local,ios::out);

vector<string> files;
for(int ix=0; ix<nxmass; ix++){
  for(int iy=0; iy<nymass; iy++){

	if(xymass_DL[ix][1+iy]>xymass_DL[ix][0]) continue;
	if(xymass_DL[ix][0]==1000 && xymass_DL[ix][1+iy]>800) continue;
	if(xymass_DL[ix][0]==500 &&  (xymass_DL[ix][1+iy]>300||xymass_DL[ix][1+iy]==190)) continue;
	if(xymass_DL[ix][0]==600 &&  (xymass_DL[ix][1+iy]>400||xymass_DL[ix][1+iy]==350||xymass_DL[ix][1+iy]==190)) continue;
	if(xymass_DL[ix][0]==650 &&  (xymass_DL[ix][1+iy]>500)) continue;
	if(xymass_DL[ix][0]==700 &&  (xymass_DL[ix][1+iy]>500||xymass_DL[ix][1+iy]==350||xymass_DL[ix][1+iy]==450||xymass_DL[ix][1+iy]==190)) continue;
	if(xymass_DL[ix][0]==800 &&  (xymass_DL[ix][1+iy]>600||xymass_DL[ix][1+iy]==350||xymass_DL[ix][1+iy]==450||xymass_DL[ix][1+iy]==190)) continue;
	if(xymass_DL[ix][0]==900 &&  (xymass_DL[ix][1+iy]>700||xymass_DL[ix][1+iy]==350||xymass_DL[ix][1+iy]==450||xymass_DL[ix][1+iy]==190)) continue;
	
	char signal_name[500];
	if(isDL){
	sprintf(signal_name,"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-%i_MY-%i_TuneCP5_13TeV-madgraph-pythia8",xymass_DL[ix][0],xymass_DL[ix][1+iy]);
	}
	else{
	sprintf(signal_name,"NMSSM_XToYHTo2B2WTo2B2Q1L1Nu_MX-%i_MY-%i_TuneCP5_13TeV-madgraph-pythia8",xymass_DL[ix][0],xymass_DL[ix][1+iy]);
	}
	cout<<"SIG: "<<signal_name<<endl;
	files.push_back(string(signal_name));
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
   file <<"cmssw-el7\n";
   file << "source /cvmfs/cms.cern.ch/cmsset_default.sh\n";
   file <<"cd "<<path<<" \n";
   file<<"export X509_USER_PROXY="<<proxy<<"\n";
   file<<"eval `scramv1 runtime -sh`\n";
   file<<"./histomaker_comb.exe  "<<isDL<<" "<<" 0 "<<" "<<files[ii] << ".root 1" ;

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
file_local<<"sh "<<name_buffer<<"\n";
}

file_sub.close();
  
//return 0;
}

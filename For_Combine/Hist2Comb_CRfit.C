#include <math.h>

void check_zero_bin(TH1D *hin)
{
	for(int bn=0; bn<=hin->GetNbinsX(); bn++){
	  if((hin->GetBinContent(bn)) < 0){ hin->SetBinContent(bn,0);  }
	}
	//hin->Rebin(2);
}

void make_positive_bins(TH1D *hin){
	for(int bn=0; bn<hin->GetNbinsX(); bn++){
		hin->SetBinContent(bn+1,fabs(hin->GetBinContent(bn+1)));
	}
}

void Make_Average(TH1D *hin){
	for(int bn=0; bn<hin->GetNbinsX(); bn++){
		hin->SetBinContent(bn+1,0.5*hin->GetBinContent(bn+1));
		hin->SetBinError(bn+1,0.5*hin->GetBinError(bn+1));
	}
}


//const char* bkg_names[] = {"TT","ST","DYj","Wj","Diboson","QCD"};

vector<pair<const char *,const char *>> samples;

const char* sysnames[] = {
	 "JES_AbsoluteStat", "JES_AbsoluteScale","JES_AbsoluteMPFBias", 
	 "JES_FlavorQCD", "JES_Fragmentation", 
	 "JES_PileUpDataMC",  "JES_PileUpPtBB", "JES_PileUpPtEC1", "JES_PileUpPtEC2", 
	 "JES_PileUpPtRef",
	 "JES_RelativeFSR", "JES_RelativeJEREC1", "JES_RelativeJEREC2", 
	 "JES_RelativePtBB", "JES_RelativePtEC1", "JES_RelativePtEC2", 
	 "JES_RelativeBal", "JES_RelativeSample", "JES_RelativeStatEC", "JES_RelativeStatFSR", 
	 "JES_SinglePionECAL", "JES_SinglePionHCAL","JES_TimePtEta",
	 "JES_Total",
	 "JER",
	 "PU","LeptonSF","LeptonSF2","Prefire","PNbbSF","PNWSF","BTG","TrigSF1","TrigSF2"
}; 

struct sample_info
{
string name;
float norm;	
};

void Hist2Comb_CRfit(int year=2018, bool isDL=false)
{

samples.push_back(make_pair("Top","_Top_fullymerged"));
samples.push_back(make_pair("Top","_Top_semimerged"));
samples.push_back(make_pair("Top","_Top_unmerged"));
samples.push_back(make_pair("DYj",""));
samples.push_back(make_pair("Wj",""));
samples.push_back(make_pair("Diboson",""));
samples.push_back(make_pair("QCD",""));

//int nbkg = sizeof(bkg_names)/sizeof(bkg_names[0]);
int nbkg = samples.size();
int nsys = sizeof(sysnames)/sizeof(sysnames[0]);

float data_lumi = 139.;

if(year==2016) { data_lumi = 35.9; }
if(year==2017) { data_lumi = 41.5; }
if(year==2018) { data_lumi = 59.7; }
cout<<"Lumi:"<<data_lumi<<endl;

string filepath, output_filepath, data_file;
if(isDL){
	filepath = "/groups/hephy/cms/suman.chatterjee/XtoYH/Histograms/Dileptonic/October2022/";
	data_file = "MuEGamma";
}
else{
	//filepath = "/groups/hephy/cms/suman.chatterjee/XtoYH/Histograms/Semileptonic/Feb2023/";   
	filepath = "/groups/hephy/cms/suman.chatterjee/XtoYH/Histograms/Semileptonic/April2023_2/";
	//"/users/suman.chatterjee/XtoYH/CMSSW_10_6_0/src/Histograms/March_2022/Singlelep/October2022/";
	data_file = "MuEGammaJetHT"; 
}
output_filepath = "";

vector<string> histnames;
//histnames.push_back("ST");
//histnames.push_back("unrolled_HTlep_pt");
histnames.push_back("Y_msoftdrop");
histnames.push_back("Y_msoftdrop_xbin");
histnames.push_back("unrolled_ST");
if(!isDL){
	histnames.push_back("unrolled_bin1_X_mass");
	histnames.push_back("unrolled_bin2_X_mass");
}

const char *catnames[] = {"opt2"};
const char *regnames[] = {"CR2_nb1","CR4_nb1","CR6_nb1","CR3_nb0"};
const char *lepids[] = {"","_Mu","_El","_EMu"};

int nvar = int(histnames.size()); //sizeof(histnames)/sizeof(histnames[0]);
int ncat = sizeof(catnames)/sizeof(catnames[0]);
int nreg = sizeof(regnames)/sizeof(regnames[0]);
int nlid = sizeof(lepids)/sizeof(lepids[0]);

TH1D *hist_b[ncat][nreg][nlid][nvar][nbkg];
TH1D *hist_b_up[ncat][nreg][nlid][nvar][nbkg][nsys];
TH1D *hist_b_dn[ncat][nreg][nlid][nvar][nbkg][nsys];

TH1D *hist_data[ncat][nreg][nlid][nvar];

char name[1000];

for(int fg=0; fg<nbkg; fg++){
	
	sprintf(name,"%s/Output_%s.root",(filepath).c_str(),samples[fg].first);
	cout<<name<<endl;
	TFile *_file0 = new TFile(name,"read");
	
	for(int icat=0; icat<ncat; icat++){
		for(int ireg=0; ireg<nreg; ireg++){
			for(int lc=0; lc<nlid; lc++){
					
				for(int ivar=0; ivar<(nvar); ivar++){
			
					sprintf(name,"h_Y_T_W_T_%s_%s_%s%s",histnames[ivar].c_str(),catnames[icat],regnames[ireg],lepids[lc]);
					
					char histname[200];
					sprintf(histname,"%s%s",name,samples[fg].second);
					
					hist_b[icat][ireg][lc][ivar][fg] = (TH1D*)_file0->Get(histname);
					hist_b[icat][ireg][lc][ivar][fg]->Scale(data_lumi);
					check_zero_bin(hist_b[icat][ireg][lc][ivar][fg]);
					
					for(int isys=0; isys<nsys; isys++){
						char name_sys[300];
						sprintf(name_sys,"%s_%s_up",histname,sysnames[isys]);
						hist_b_up[icat][ireg][lc][ivar][fg][isys] = (TH1D*)_file0->Get(name_sys);
						hist_b_up[icat][ireg][lc][ivar][fg][isys]->Scale(data_lumi);
						check_zero_bin(hist_b_up[icat][ireg][lc][ivar][fg][isys]);
						sprintf(name_sys,"%s_%s_dn",histname,sysnames[isys]);
						hist_b_dn[icat][ireg][lc][ivar][fg][isys] = (TH1D*)_file0->Get(name_sys);
						hist_b_dn[icat][ireg][lc][ivar][fg][isys]->Scale(data_lumi);
						check_zero_bin(hist_b_dn[icat][ireg][lc][ivar][fg][isys]);
					}
						
					// end //
	
				}//var
			}//lc
		}//reg
	}//cat
	
//	_file0->Close();
	
}//files

sprintf(name,"%s/Output_%s.root",(filepath).c_str(),(data_file).c_str());
TFile *_file_data = new TFile(name,"read");

for(int icat=0; icat<ncat; icat++){
	for(int ireg=0; ireg<nreg; ireg++){
		for(int lc=0; lc<nlid; lc++){
					
			for(int ivar=0; ivar<(nvar); ivar++){
			
				sprintf(name,"h_Y_T_W_T_%s_%s_%s%s",histnames[ivar].c_str(),catnames[icat],regnames[ireg],lepids[lc]);
			
				hist_data[icat][ireg][lc][ivar] = (TH1D*)_file_data->Get(name);
				check_zero_bin(hist_data[icat][ireg][lc][ivar]);
	
				}//var
			}//lc
		}//reg
	}//cat

// Now write all these to a root file using convention of CombineTool //

TFile *fileout;

if(year!=2016 && year!=2017 && year!=2018){
	if(isDL){ sprintf(name,"Combine_input_XYH_CRFit_Run2_DL.root");	}
	else { sprintf(name,"Combine_input_XYH_CRFit_Run2_SL.root");	}
}
else{
	if(isDL){ sprintf(name,"Combine_input_XYH_CRFit_%i_DL.root",year); }
	else { sprintf(name,"Combine_input_XYH_CRFit_%i_SL.root",year); }
}

cout<<"Producing "<<name<<endl;
fileout = new TFile(name,"RECREATE");

for(int icat=0; icat<ncat; icat++){
		for(int ireg=0; ireg<nreg; ireg++){
			for(int lc=0; lc<nlid; lc++){
				for(int ivar=0; ivar<(nvar); ivar++){

					TDirectory *dir;
					fileout->cd();

					sprintf(name,"%s_%s%s_%s",catnames[icat],regnames[ireg],lepids[lc],histnames[ivar].c_str());
					dir = fileout->mkdir(name);
					dir->cd();

					for(int fg=0; fg<nbkg; fg++){

						TH1D *h_bkg = (TH1D*)hist_b[icat][ireg][lc][ivar][fg]->Clone();
						sprintf(name,"%s%s",samples[fg].first,samples[fg].second);
						h_bkg->SetName(name);
						//h_bkg->Write();
				
						for(int isys=0; isys<nsys; isys++){
					
							h_bkg = (TH1D*)hist_b_up[icat][ireg][lc][ivar][fg][isys]->Clone();
							sprintf(name,"%s%s_%sUp",samples[fg].first,samples[fg].second,sysnames[isys]);
							h_bkg->SetName(name);
							//h_bkg->Write();
					
							h_bkg = (TH1D*)hist_b_dn[icat][ireg][lc][ivar][fg][isys]->Clone();
							sprintf(name,"%s%s_%sDown",samples[fg].first,samples[fg].second,sysnames[isys]);
							h_bkg->SetName(name);
							//h_bkg->Write();	
					
						}
					}
			
					TH1D *h_data = (TH1D*)hist_data[icat][ireg][lc][ivar]->Clone();						
					h_data->SetName("data_obs");
					//h_data->Write();	
						
			}//ivar
		}//lc
	}//ireg
}//icat

fileout->Write();
fileout->Close();

}

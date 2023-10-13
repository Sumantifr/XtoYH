#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCut.h"
#include "TPaveText.h"
#include "TGaxis.h"
#include "TProfile.h"
#include "TCut.h"
#include <iostream>
#include <iostream>
#include <cstring>
#include <string>
#include <vector>

#include "Functions.h"

using namespace std;

int get_Y_id(float  Y_msoftdrop ) 
{
	int Y_id;
    if ( Y_msoftdrop < 50.0  )                               Y_id = 0;
    else if ( Y_msoftdrop >= 50.0  &&  Y_msoftdrop < 100.0 ) Y_id = 1;
    else if ( Y_msoftdrop >= 100.0 &&  Y_msoftdrop < 150.0 ) Y_id = 2;
    else if ( Y_msoftdrop >= 150.0 &&  Y_msoftdrop < 200.0 ) Y_id = 3;
    else if ( Y_msoftdrop >= 200.0 &&  Y_msoftdrop < 250.0 ) Y_id = 4;
    else if ( Y_msoftdrop >= 250.0 &&  Y_msoftdrop < 300.0 ) Y_id = 5;
    else if ( Y_msoftdrop >= 300.0 &&  Y_msoftdrop < 350.0 ) Y_id = 6;
    else if ( Y_msoftdrop >= 350.0 &&  Y_msoftdrop < 400.0 ) Y_id = 7;
    else if ( Y_msoftdrop >= 400.0 &&  Y_msoftdrop < 450.0 ) Y_id = 8;
    else if ( Y_msoftdrop >= 450.0 &&  Y_msoftdrop < 500.0 ) Y_id = 9;
    else if ( Y_msoftdrop >= 500.0 &&  Y_msoftdrop < 550.0 ) Y_id = 10;
    else if ( Y_msoftdrop >= 550.0 &&  Y_msoftdrop < 600.0 ) Y_id = 11;
    return Y_id;
}

int get_region(vector<bool> reg_tags){

  int freg = -1;
  for(unsigned ireg=0; ireg<reg_tags.size(); ireg++){
	  if(reg_tags[ireg]){
		  freg = int(ireg);
		  break;
		  }
	  }
  return freg;
  
}

vector<int> get_regions(vector<bool> reg_tags){

  vector<int> regs;

  int freg = -1;
  for(unsigned ireg=0; ireg<reg_tags.size(); ireg++){
	  if(reg_tags[ireg]){
		  freg = int(ireg);
		  regs.push_back(freg);
		  }
	  }
  return regs;
  
}

TH1F* get_histo_symbin(TString Y_op_name, TString W_op_name, TString reg_name, TString bcat_name, TString Wops_name, TString lep_name, string var, string addText, int nbins, float low_edge, float up_edge)
{
char name[200];
sprintf(name,"h_Y_%s_W_%s_%s_%s_%s%s%s%s",Y_op_name.Data(),W_op_name.Data(),var.c_str(),Wops_name.Data(),reg_name.Data(),bcat_name.Data(),lep_name.Data(),addText.c_str());
TH1F *hout = new TH1F(name, "",nbins,low_edge,up_edge);
hout->Sumw2();
return hout;
}

TH1F* get_histo_symbin_II(TString Y_op_name, TString W_op_name, TString reg_name, TString bcat_name, TString Wops_name, TString lep_name, TString top_name, string var, string addText, int nbins, float low_edge, float up_edge)
{
char name[500];
sprintf(name,"h_Y_%s_W_%s_%s_%s_%s%s%s%s%s",Y_op_name.Data(),W_op_name.Data(),var.c_str(),Wops_name.Data(),reg_name.Data(),bcat_name.Data(),lep_name.Data(),top_name.Data(),addText.c_str());
TH1F *hout = new TH1F(name, "",nbins,low_edge,up_edge);
hout->Sumw2();
return hout;
}

TH1F* get_histo_symbin_onetopreg(TString Y_op_name, TString W_op_name, TString Wops_name, TString lep_name, TString top_name, string var, string addText, int nbins, float low_edge, float up_edge)
{
char name[500];
sprintf(name,"h_Y_%s_W_%s_%s_%s_%s%s%s",Y_op_name.Data(),W_op_name.Data(),var.c_str(),Wops_name.Data(),lep_name.Data(),top_name.Data(),addText.c_str());
TH1F *hout = new TH1F(name, "",nbins,low_edge,up_edge);
hout->Sumw2();
return hout;
}


TH1F* get_histo_asymbin(TString Y_op_name, TString W_op_name, TString reg_name, TString bcat_name, TString Wops_name, TString lep_name, string var, string addText, int nbins, float *bins)
{
char name[200];
sprintf(name,"h_Y_%s_W_%s_%s_%s_%s%s%s%s",Y_op_name.Data(),W_op_name.Data(),var.c_str(),Wops_name.Data(),reg_name.Data(),bcat_name.Data(),lep_name.Data(),addText.c_str());
TH1F *hout = new TH1F(name, "",nbins,bins);
hout->Sumw2();
return hout;
}

TH1F* get_histo_asymbin_II(TString Y_op_name, TString W_op_name, TString reg_name, TString bcat_name, TString Wops_name, TString lep_name, TString top_name, string var, string addText, int nbins, float *bins)
{
char name[500];
sprintf(name,"h_Y_%s_W_%s_%s_%s_%s%s%s%s%s",Y_op_name.Data(),W_op_name.Data(),var.c_str(),Wops_name.Data(),reg_name.Data(),bcat_name.Data(),lep_name.Data(),top_name.Data(),addText.c_str());
TH1F *hout = new TH1F(name, "",nbins,bins);
hout->Sumw2();
return hout;
}

TH1F* get_histo_asymbin_onetopreg(TString Y_op_name, TString W_op_name, TString Wops_name, TString lep_name, TString top_name, string var, string addText, int nbins, float *bins)
{
char name[500];
sprintf(name,"h_Y_%s_W_%s_%s_%s_%s%s%s",Y_op_name.Data(),W_op_name.Data(),var.c_str(),Wops_name.Data(),lep_name.Data(),top_name.Data(),addText.c_str());
TH1F *hout = new TH1F(name, "",nbins,bins);
hout->Sumw2();
return hout;
}



TH1F* getHisto1F(const char *name, const char *title, int nbins, float low_edge, float up_edge)
{
TH1F *hout = new TH1F(name,title,nbins,low_edge,up_edge);
hout->Sumw2();
return hout;
}

TH1F* getHisto1F(const char *name, const char *title, int nbins, float *bins)
{
TH1F *hout = new TH1F(name,title,nbins,bins);
hout->Sumw2();
return hout;
}

TH1D* getHisto1D(const char *name, const char *title, int nbins, float low_edge, float up_edge)
{
TH1D *hout = new TH1D(name,title,nbins,low_edge,up_edge);
hout->Sumw2();
return hout;
}

TH1D* getHisto1D(const char *name, const char *title, int nbins, float *bins)
{
TH1D *hout = new TH1D(name,title,nbins,bins);
hout->Sumw2();
return hout;
}

static const int PNbb_SF_nptbins = 4;
float PNbb_SF_ptbins[PNbb_SF_nptbins+1] = {300,400,500,600,10000};
double PNbb_SF_HP[PNbb_SF_nptbins] = {1.192,1.137,1.211,1.350};
double PNbb_SF_HP_up[PNbb_SF_nptbins] = {1.192-0.101,1.137-0.081,1.211-0.151,1.350-0.237};
double PNbb_SF_HP_dn[PNbb_SF_nptbins] = {1.192+0.101,1.137+0.082,1.211+0.151,1.350+0.238};
// SFs taken from: https://indico.cern.ch/event/1011640/contributions/4460748/attachments/2285317/3885997/21.07.21_BTV_ParticleNet%20SFs%20for%20UL1718%20v2.pdf

static const int PNW_SF_nptbins = 3;
float PNW_SF_ptbins[PNW_SF_nptbins+1] = {200,300,400,10000};
double PNW_SF_T[PNW_SF_nptbins] = {0.81,0.81,0.77};
double PNW_SF_T_up[PNW_SF_nptbins] = {0.81+0.02,0.81+0.02,0.77+0.04};
double PNW_SF_T_dn[PNW_SF_nptbins] = {0.81-0.02,0.81-0.02,0.77-0.04};
// SFs taken from: https://indico.cern.ch/event/1103765/contributions/4647556/attachments/2364610/4037250/ParticleNet_2018_ULNanoV9_JMAR_14Dec2021_PK.pdf

static const int PNTop_SF_nptbins = 4;
float PNTop_SF_ptbins[PNTop_SF_nptbins+1] = {300,400,480,600,100000};
double PNTop_SF_M[PNTop_SF_nptbins] = {1.12,0.98,0.97,0.99};
double PNTop_SF_M_up[PNTop_SF_nptbins] = {1.12+0.12,0.98+0.04,0.97+0.03,0.99+0.05};
double PNTop_SF_M_dn[PNTop_SF_nptbins] = {1.12-0.07,0.98-0.03,0.97-0.03,0.99-0.05};
// SFs taken from: https://indico.cern.ch/event/1103765/contributions/4647556/attachments/2364610/4037250/ParticleNet_2018_ULNanoV9_JMAR_14Dec2021_PK.pdf

float top_pt_cor_MC = 1.;

float PN_Top_med = 0.8;
float deep_btag_cut = 0.2783; 

//masscut on AK8 jets //
float msd_cut = 30;

// SR2 cut on bb tagging score //

//Z veto window //
float Z_mass_min = 75;
float Z_mass_max = 105;//120;

//mini-isolation cut for lepton//
float miniso_cut = 0.1;

//PNbb score cut for SR2 //
float PNbb_cut_SR2 = 0.4;

//trigger cuts //

float mu_trig_pt_SL = 50;
float el_trig_pt_SL = 32;
float jet_trig_pt_SL = 550;
float mu_trig_pt_DL_0 = 37;
float mu_trig_pt_DL_1 = 27;
float el_trig_pt_DL   = 25;
float emu_trig_pt_DL_0 = 37;
float emu_trig_pt_DL_1 = 27;

double SF_Trig, SF_Trig_stat, SF_Trig_syst;
double SF_Trig_1_up, SF_Trig_1_dn;
double SF_Trig_2_up, SF_Trig_2_dn;


double b_SF, b_SF_up, b_SF_dn;
double bb_SF, bb_SF_up, bb_SF_dn;
double W_SF, W_SF_up, W_SF_dn;
double Top_SF, Top_SF_up, Top_SF_dn;

double analysis_SF, analysis_SF_up, analysis_SF_dn;

TString proc_Name[] = {
//"TTTo2L2Nu_XtoYH.root"
/*
"DYJetsToLL_M-10to50_XtoYH.root",
"DYJetsToLL_M-50_HT-100To200_XtoYH.root",
"DYJetsToLL_M-50_HT-1200To2500_XtoYH.root",
"DYJetsToLL_M-50_HT-200To400_XtoYH.root",
"DYJetsToLL_M-50_HT-2500ToInf_XtoYH.root",
"DYJetsToLL_M-50_HT-400To600_XtoYH.root",
"DYJetsToLL_M-50_HT-600To800_XtoYH.root",
"DYJetsToLL_M-50_HT-70to100_XtoYH.root",
"DYJetsToLL_M-50_HT-800To1200_XtoYH.root",
"QCD_HT1000to1500_XtoYH.root",
"QCD_HT1500to2000_XtoYH.root",
"QCD_HT2000toInf_XtoYH.root",
"QCD_HT700to1000_XtoYH.root",
"QCD_HT300to500_XtoYH.root",
"QCD_HT500to700_XtoYH.root",
"ST_s-channel_XtoYH.root",
"ST_t-channel_antitop_XtoYH.root",
"ST_t-channel_top_XtoYH.root",
"ST_tW_antitop_XtoYH.root",
"ST_tW_top_XtoYH.root",
"TTTo2L2Nu_XtoYH.root",
"TTToHadronic_XtoYH.root",
"TTToSemiLeptonic_XtoYH.root",
"WJetsToLNu_HT-100To200_XtoYH.root",
"WJetsToLNu_HT-1200To2500_XtoYH.root",
"WJetsToLNu_HT-200To400_XtoYH.root",
"WJetsToLNu_HT-2500ToInf_XtoYH.root",
"WJetsToLNu_HT-400To600_XtoYH.root",
"WJetsToLNu_HT-600To800_XtoYH.root",
"WJetsToLNu_HT-70To100_XtoYH.root",
"WJetsToLNu_HT-800To1200_XtoYH.root",
"WWTo1L1Nu2Q_XtoYH.root",
"WWTo2L2Nu_XtoYH.root",
"WWTo4Q_XtoYH.root",
"WZTo1L1Nu2Q_XtoYH.root",
"WZTo2Q2L_XtoYH.root",
"WZTo2Q2Nu_XtoYH.root",
"WZTo3LNu_XtoYH.root",
"ZZTo2L2Nu_XtoYH.root",
"ZZTo2Q2L_XtoYH.root",
"ZZTo2Q2Nu_XtoYH.root",
"ZZTo4L_XtoYH.root",
"ZZTo4Q_XtoYH.root",
*/
"WJetsToLNu_HT-800To1200_XtoYH.root"
/*
"NMSSM_XYH_YTobb_HToWWTo2L2Nu_MX_2000_MY_200_v3.root",
"NMSSM_XYH_YTobb_HToWWTo2L2Nu_MX_3000_MY_100_v3.root",
"NMSSM_XYH_YTobb_HToWWTo2L2Nu_MX_1500_MY_200_v3.root",
"NMSSM_XYH_YTobb_HToWWTo2L2Nu_MX_3000_MY_500_v3.root"
*/
/*
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_3000_MY_500.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_3000_MY_100.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2400_MY_300.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2000_MY_200.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_1500_MY_200.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_1000_MY_100.root",
*/
/*
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_1000_MY_100_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_1500_MY_200_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_1500_MY_250_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2000_MY_125_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2000_MY_200_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2400_MY_300_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2600_MY_125_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_3000_MY_100_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_3000_MY_300_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_3000_MY_500_TuneCP5_13TeV-madgraph-pythia8.root"
*/
//Running on the data samples
/*
"EGamma_UL2018A_XtoYH_Nov_2021.root",
"EGamma_UL2018B_XtoYH_Nov_2021.root",
"EGamma_UL2018C_XtoYH_Nov_2021.root",
"EGamma_UL2018D_XtoYH_Nov_2021.root",
"SingleMuon_UL2018A_XtoYH_Nov_2021.root",
"SingleMuon_UL2018B_XtoYH_Nov_2021.root",
"SingleMuon_UL2018C_XtoYH_Nov_2021.root",
"SingleMuon_UL2018D_XtoYH_Nov_2021.root",
"JetHT_UL2018A_XtoYH_Nov_2021.root",
"JetHT_UL2018B_XtoYH_Nov_2021.root",
"JetHT_UL2018C_XtoYH_Nov_2021.root",
"JetHT_UL2018D_XtoYH_Nov_2021.root"
*/
//Datasets for the Dilepton analysis part
/*
"MuonEG_UL2018A_XtoYH_Nov_2021.root",
"MuonEG_UL2018B_XtoYH_Nov_2021.root",
"MuonEG_UL2018C_XtoYH_Nov_2021.root",
"MuonEG_UL2018D_XtoYH_Nov_2021.root",
"EGamma_UL2018A_XtoYH_Nov_2021.root",
"EGamma_UL2018B_XtoYH_Nov_2021.root",
"EGamma_UL2018C_XtoYH_Nov_2021.root",
"EGamma_UL2018D_XtoYH_Nov_2021.root",
"DoubleMuon_UL2018A_XtoYH_Nov_2021.root",
"DoubleMuon_UL2018B_XtoYH_Nov_2021.root",
"DoubleMuon_UL2018C_XtoYH_Nov_2021.root",
"DoubleMuon_UL2018D_XtoYH_Nov_2021.root"
*/
/*
"NMSSM_XYH_YTobb_HToWWTo2L2Nu_MX_1500_MY_200_v3.root",
"NMSSM_XYH_YTobb_HToWWTo2L2Nu_MX_2000_MY_200_v3.root",
"NMSSM_XYH_YTobb_HToWWTo2L2Nu_MX_3000_MY_100_v3.root",
"NMSSM_XYH_YTobb_HToWWTo2L2Nu_MX_3000_MY_500_v3.root"
*/
/*
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1000_MY-100_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1000_MY-125_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1000_MY-150_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1000_MY-250_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1000_MY-300_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1000_MY-400_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1000_MY-500_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1000_MY-60_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1000_MY-70_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1000_MY-80_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1000_MY-90_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1200_MY-100_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1200_MY-125_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1200_MY-150_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1200_MY-190_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1200_MY-250_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1200_MY-300_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1200_MY-350_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1200_MY-400_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1200_MY-500_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1200_MY-600_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1200_MY-60_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1200_MY-70_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1200_MY-80_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1200_MY-90_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1400_MY-100_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1400_MY-125_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1400_MY-150_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1400_MY-190_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1400_MY-250_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1400_MY-300_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1400_MY-400_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1400_MY-500_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1400_MY-600_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1400_MY-60_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1400_MY-70_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1400_MY-80_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1600_MY-100_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1600_MY-125_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1600_MY-150_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1600_MY-250_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1600_MY-300_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1600_MY-350_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1600_MY-400_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1600_MY-500_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1600_MY-600_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1600_MY-60_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1600_MY-70_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1600_MY-80_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1600_MY-90_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1800_MY-100_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1800_MY-125_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1800_MY-150_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1800_MY-300_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1800_MY-400_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1800_MY-500_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1800_MY-600_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1800_MY-60_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1800_MY-70_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1800_MY-80_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-1800_MY-90_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2000_MY-125_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2000_MY-250_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2000_MY-300_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2000_MY-400_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2000_MY-500_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2000_MY-600_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2000_MY-60_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2000_MY-70_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2000_MY-80_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2000_MY-90_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2500_MY-100_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2500_MY-125_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2500_MY-150_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2500_MY-190_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2500_MY-250_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2500_MY-300_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2500_MY-350_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2500_MY-400_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2500_MY-500_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2500_MY-600_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2500_MY-60_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2500_MY-70_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2500_MY-80_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-2500_MY-90_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-3000_MY-100_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-3000_MY-125_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-3000_MY-150_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-3000_MY-250_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-3000_MY-300_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-3000_MY-400_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-3000_MY-500_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-3000_MY-600_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-3000_MY-60_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-3000_MY-70_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-3000_MY-80_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-3500_MY-125_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-3500_MY-450_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-3500_MY-600_TuneCP5_13TeV-madgraph-pythia8.root",
"NMSSM_XToYHTo2B2WTo2B2L2Nu_MX-3500_MY-60_TuneCP5_13TeV-madgraph-pythia8.root"
*/
};
////

  static const int narray = 20;
  static const int njetmx = 6;
  static const int ncutmax = 7;
  static const int nlhescalemax = 9;
  static const int nlhepdfmax = 101;
  static const int nalpsmax = 3;
  static const int nlhepsmax = 8;
  
  Int_t           irun;
  Int_t           ilumi;
  UInt_t          ievt;
  Int_t			  npvert;
  Int_t 		  PV_npvsGood;
   
  Int_t           nleptons;
  Int_t           nfatjets;
  Int_t           ncuts;
  Bool_t          Flag_event_cuts[ncutmax];  
  Bool_t	  	  Flag_pass_baseline;  
  Bool_t	  	  Flag_pass_baseline_no_LJet;
 
  Bool_t          hlt_IsoMu24;
  Bool_t          hlt_Mu50;
  Bool_t          hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;
  Bool_t          hlt_Ele115_CaloIdVT_GsfTrkIdT;
  Bool_t          hlt_Ele40_WPTight_Gsf;
  Bool_t          hlt_Ele32_WPTight_Gsf;
  Bool_t          hlt_Ele28_eta2p1_WPTight_Gsf_HT150;
  Bool_t          hlt_Mu37_Ele27_CaloIdL_MW;
  Bool_t          hlt_Mu27_Ele37_CaloIdL_MW;
  Bool_t          hlt_Mu37_TkMu27;
  Bool_t          hlt_DoubleEle25_CaloIdL_MW;
  Bool_t          hlt_AK8PFJet500;
  Bool_t          hlt_PFJet500;
  Bool_t          hlt_HT1050;
  Bool_t          hlt_AK8PFJet400_TrimMass30;
  Bool_t          hlt_AK8PFHT800_TrimMass50;
  //Bool_t          hlt_Photon200;
  //Bool_t          hlt_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60; //exist?
  //Bool_t          hlt_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60; //exist?
  //Bool_t          hlt_PFMETNoMu140_PFMHTNoMu140_IDTight; //exist?
  //Bool_t          hlt_PFMETTypeOne140_PFMHT140_IDTight; //exist?
  
  Bool_t          Muon_trig_pass;
  Bool_t          Electron_trig_pass;
  Bool_t          MuonElectron_trig_pass;
  
  Double_t        LHE_weight;
  Double_t        Generator_weight;
  Double_t        Event_weight;
  
  Double_t        prefiringweight;
  Double_t        prefiringweightup;
  Double_t        prefiringweightdown;
  Float_t         puWeight;
  Float_t         puWeightup;
  Float_t         puWeightdown;
  Float_t         leptonsf_weight;
  Float_t         leptonsf_weight_up;
  Float_t         leptonsf_weight_dn;
  Float_t         leptonsf_weight_stat;
  Float_t         leptonsf_weight_syst;
  
  Float_t         l1_pt;
  Float_t         l1_eta;
  Float_t         l1_phi;
  Float_t         l1_mass;
  Int_t           l1_pdgId;
  Float_t         l1_minisoch;
  Float_t         l1_minisonh;
  Float_t         l1_minisoph;
  Float_t         l1_minisoall;
  Int_t           l1_genindex;
  
  Float_t         l2_pt;
  Float_t         l2_eta;
  Float_t         l2_phi;
  Float_t         l2_mass;
  Int_t           l2_pdgId;
  Float_t         l2_minisoch;
  Float_t         l2_minisonh;
  Float_t         l2_minisoph;
  Float_t         l2_minisoall;
  Int_t           l2_genindex;
  
  Float_t         MET_pt;
  Float_t         MET_phi;
  Float_t         MET_sig;
  Float_t         MET_sumEt;
  Float_t         MET_pt_JESup;
  Float_t         MET_pt_JESdn;
  Float_t         MET_pt_JERup;
  Float_t         MET_pt_JERdn;
  Float_t         MET_pt_UnclusEup;
  Float_t         MET_pt_UnclusEdn;
  vector<float>   *MET_pt_JESup_split;
  vector<float>   *MET_pt_JESdn_split;
  Float_t         MET_phi_JESup;
  Float_t         MET_phi_JESdn;
  Float_t         MET_phi_JERup;
  Float_t         MET_phi_JERdn;
  Float_t         MET_phi_UnclusEup;
  Float_t         MET_phi_UnclusEdn;
  vector<float>   *MET_phi_JESup_split;
  vector<float>   *MET_phi_JESdn_split;
  
  Float_t         Y_pt;
  Float_t         Y_y;
  Float_t         Y_eta;
  Float_t         Y_phi;
  Float_t         Y_mass;
  Float_t         Y_msoftdrop;
  Float_t         Y_tau21;
  Float_t         Y_tau32;
  Float_t         Y_DeepTag_DAK8MD_TvsQCD;
  Float_t         Y_DeepTag_DAK8MD_WvsQCD;
  Float_t         Y_DeepTag_DAK8MD_ZvsQCD;
  Float_t         Y_DeepTag_DAK8MD_HvsQCD;
  Float_t         Y_DeepTag_DAK8MD_bbvsQCD;
  Float_t         Y_DeepTag_PNet_TvsQCD;
  Float_t         Y_DeepTag_PNet_WvsQCD;
  Float_t         Y_DeepTag_PNet_ZvsQCD;
  Float_t         Y_DeepTag_PNetMD_XbbvsQCD;
  Float_t         Y_DeepTag_PNetMD_XccvsQCD;
  Float_t         Y_DeepTag_PNetMD_XqqvsQCD;
  Float_t         Y_DeepTag_PNetMD_WvsQCD;
  Float_t         Y_DeepTag_PNetMD_QCD;
  Float_t         Y_PN_bb;
  Bool_t          Y_label_Top_bq;
  Bool_t          Y_label_Top_bc;
  Bool_t          Y_label_Top_bcq;
  Bool_t          Y_label_Top_bqq;
  Bool_t          Y_label_W_qq;
  Bool_t          Y_label_W_cq;
  Float_t         Y_sub1_pt;
  Float_t         Y_sub1_eta;
  Float_t         Y_sub1_phi;
  Float_t         Y_sub1_mass;
  Float_t         Y_sub1_btag;
  Float_t         Y_sub2_pt;
  Float_t         Y_sub2_eta;
  Float_t         Y_sub2_phi;
  Float_t         Y_sub2_mass;
  Float_t         Y_sub2_btag;
  Int_t           Y_genindex;
  Int_t           Y_genbindex[2];
  Float_t         Y_JESup;
  Float_t         Y_JESdn;
  Float_t         Y_JERup;
  Float_t         Y_JERdn;
  vector<float>   *Y_JESup_split;
  vector<float>   *Y_JESdn_split;
  
  Float_t         W_pt_opt1;
  Float_t         W_y_opt1;
  Float_t         W_eta_opt1;
  Float_t         W_phi_opt1;
  Float_t         W_mass_opt1;
  Float_t         W_msoftdrop_opt1;
  Float_t         W_tau21_opt1;
  Float_t         W_tau32_opt1;
  Float_t         W_DeepTag_DAK8MD_TvsQCD_opt1;
  Float_t         W_DeepTag_DAK8MD_WvsQCD_opt1;
  Float_t         W_DeepTag_DAK8MD_ZvsQCD_opt1;
  Float_t         W_DeepTag_DAK8MD_HvsQCD_opt1;
  Float_t         W_DeepTag_DAK8MD_bbvsQCD_opt1;
  Float_t         W_DeepTag_PNet_TvsQCD_opt1;
  Float_t         W_DeepTag_PNet_WvsQCD_opt1;
  Float_t         W_DeepTag_PNet_ZvsQCD_opt1;
  Float_t         W_DeepTag_PNetMD_XbbvsQCD_opt1;
  Float_t         W_DeepTag_PNetMD_XccvsQCD_opt1;
  Float_t         W_DeepTag_PNetMD_XqqvsQCD_opt1;
  Float_t         W_DeepTag_PNetMD_WvsQCD_opt1;
  Float_t         W_DeepTag_PNetMD_QCD_opt1;
  Float_t         W_DAK8_W_opt1;
  Float_t         W_PN_W_opt1;
  Bool_t          W_label_W_qq_opt1;
  Bool_t          W_label_W_cq_opt1;
  Float_t         W_sub1_pt_opt1;
  Float_t         W_sub1_eta_opt1;
  Float_t         W_sub1_phi_opt1;
  Float_t         W_sub1_mass_opt1;
  Float_t         W_sub1_btag_opt1;
  Float_t         W_sub2_pt_opt1;
  Float_t         W_sub2_eta_opt1;
  Float_t         W_sub2_phi_opt1;
  Float_t         W_sub2_mass_opt1;
  Float_t         W_sub2_btag_opt1;
  Int_t           W_genindex_opt1;
  Float_t         W_JESup_opt1;
  Float_t         W_JESdn_opt1;
  Float_t         W_JERup_opt1;
  Float_t         W_JERdn_opt1;
  vector<float>   *W_JESup_split_opt1;
  vector<float>   *W_JESdn_split_opt1;
  Float_t         H_pt_opt1;
  Float_t         H_y_opt1;
  Float_t         H_eta_opt1;
  Float_t         H_phi_opt1;
  Float_t         H_mass_opt1;
  Int_t           H_genindex_opt1;
  Float_t         H_JESup_opt1;
  Float_t         H_JESdn_opt1;
  Float_t         H_JERup_opt1;
  Float_t         H_JERdn_opt1;
  vector<float>   *H_JESup_split_opt1;
  vector<float>   *H_JESdn_split_opt1;
  Float_t         X_mass_opt1;
  vector<float>   *X_mass_JESup_split_opt1;
  vector<float>   *X_mass_JESdn_split_opt1;
  Float_t         dR_lW_opt1;
  Float_t         dy_lW_opt1;
  Float_t         dphi_lW_opt1;
 
  Float_t         W_pt_opt2;
  Float_t         W_y_opt2;
  Float_t         W_eta_opt2;
  Float_t         W_phi_opt2;
  Float_t         W_mass_opt2;
  Float_t         W_msoftdrop_opt2;
  Float_t         W_tau21_opt2;
  Float_t         W_tau32_opt2;
  Float_t         W_DeepTag_DAK8MD_TvsQCD_opt2;
  Float_t         W_DeepTag_DAK8MD_WvsQCD_opt2;
  Float_t         W_DeepTag_DAK8MD_ZvsQCD_opt2;
  Float_t         W_DeepTag_DAK8MD_HvsQCD_opt2;
  Float_t         W_DeepTag_DAK8MD_bbvsQCD_opt2;
  Float_t         W_DeepTag_PNet_TvsQCD_opt2;
  Float_t         W_DeepTag_PNet_WvsQCD_opt2;
  Float_t         W_DeepTag_PNet_ZvsQCD_opt2;
  Float_t         W_DeepTag_PNetMD_XbbvsQCD_opt2;
  Float_t         W_DeepTag_PNetMD_XccvsQCD_opt2;
  Float_t         W_DeepTag_PNetMD_XqqvsQCD_opt2;
  Float_t         W_DeepTag_PNetMD_WvsQCD_opt2;
  Float_t         W_DeepTag_PNetMD_QCD_opt2;
  Float_t         W_DAK8_W_opt2;
  Float_t         W_PN_W_opt2;
  Bool_t          W_label_W_qq_opt2;
  Bool_t          W_label_W_cq_opt2;
  Float_t         W_sub1_pt_opt2;
  Float_t         W_sub1_eta_opt2;
  Float_t         W_sub1_phi_opt2;
  Float_t         W_sub1_mass_opt2;
  Float_t         W_sub1_btag_opt2;
  Float_t         W_sub2_pt_opt2;
  Float_t         W_sub2_eta_opt2;
  Float_t         W_sub2_phi_opt2;
  Float_t         W_sub2_mass_opt2;
  Float_t         W_sub2_btag_opt2;
  Int_t           W_genindex_opt2;
  Float_t         W_JESup_opt2;
  Float_t         W_JESdn_opt2;
  Float_t         W_JERup_opt2;
  Float_t         W_JERdn_opt2;
  vector<float>   *W_JESup_split_opt2;
  vector<float>   *W_JESdn_split_opt2;
  Float_t         H_pt_opt2;
  Float_t         H_y_opt2;
  Float_t         H_eta_opt2;
  Float_t         H_phi_opt2;
  Float_t         H_mass_opt2;
  Int_t           H_genindex_opt2;
  Float_t         H_JESup_opt2;
  Float_t         H_JESdn_opt2;
  Float_t         H_JERup_opt2;
  Float_t         H_JERdn_opt2;
  vector<float>   *H_JESup_split_opt2;
  vector<float>   *H_JESdn_split_opt2;
  Float_t         X_mass_opt2;
  vector<float>   *X_mass_JESup_split_opt2;
  vector<float>   *X_mass_JESdn_split_opt2;
  Float_t         dR_lW_opt2;
  Float_t         dy_lW_opt2;
  Float_t         dphi_lW_opt2;
   
  Float_t         dR_l1Y;
  Float_t         dy_l1Y;
  Float_t         dphi_l1Y;
  Float_t         dR_l2Y;
  Float_t         dy_l2Y;
  Float_t         dphi_l2Y;
  Float_t         l1l2_mass;
  Float_t         l1l2_deta;
  Float_t         l1l2_dphi;
  Float_t         l1l2_dR;
  Float_t         dphi_MET_l1l2;
  //additionally added  //
  Float_t         l1l2_pt;
  
  Float_t         HTlep_pt;
  Float_t         HTlep_pt_JESup;
  Float_t         HTlep_pt_JESdn;
  Float_t         HTlep_pt_JERup;
  Float_t         HTlep_pt_JERdn;
  vector<float>   *HTlep_pt_JESup_split;
  vector<float>   *HTlep_pt_JESdn_split;
  Float_t         ST;
  Float_t         ST_JESup;
  Float_t         ST_JESdn;
  Float_t         ST_JERup;
  Float_t         ST_JERdn;
  vector<float>   *ST_JESup_split;
  vector<float>   *ST_JESdn_split;
  Int_t           nbjets_other;
  Int_t           nbjets_outY;
  Int_t           nbjets_outY_L;
  Int_t           nbjets;
  Int_t           nbjets_L;
  
  Bool_t          Flag_Y_bb_pass_T;
  Bool_t          Flag_Y_bb_pass_M;
  Bool_t          Flag_Y_bb_pass_L;
  Bool_t          Flag_H_W_pass_T_opt1;
  Bool_t          Flag_H_W_pass_M_opt1;
  Bool_t          Flag_H_W_pass_L_opt1;
  Bool_t          Flag_H_m_pass_opt1;
  Bool_t          Flag_H_m_pass_opt2;
  Bool_t          Flag_dR_lW_pass_opt1;
  Bool_t          Flag_H_W_pass_T_opt2;
  Bool_t          Flag_H_W_pass_M_opt2;
  Bool_t          Flag_H_W_pass_L_opt2;
  Bool_t          Flag_dR_lW_pass_opt2;
  Bool_t          Flag_MET_pass;
  
  Int_t           nPFJetAK8;
  Float_t         PFJetAK8_pt[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_eta[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_phi[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_mass[njetmx];   //[_s_nPFJetAK8]
  Bool_t          PFJetAK8_jetID[njetmx];   //[_s_nPFJetAK8]
  Bool_t          PFJetAK8_jetID_tightlepveto[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_msoftdrop[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_tau21[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_tau32[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_DeepTag_PNetMD_XbbvsQCD[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_DeepTag_PNetMD_WvsQCD[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_DeepTag_PNet_TvsQCD[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_DeepTag_PNet_WvsQCD[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_DeepTag_DAK8MD_TvsQCD[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_DeepTag_DAK8MD_WvsQCD[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_DeepTag_DAK8MD_bbvsQCD[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_JESup[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_JESdn[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_JERup[njetmx];   //[_s_nPFJetAK8]
  Float_t         PFJetAK8_JERdn[njetmx];   //[_s_nPFJetAK8]
  Int_t           PFJetAK8_Y_index;
  Int_t           PFJetAK8_W_index_opt1;
  Int_t           PFJetAK8_W_index_opt2;
  Int_t           nJetAK4;
  Float_t         JetAK4_pt[njetmx];   //[_s_nJetAK4]
  Float_t         JetAK4_eta[njetmx];   //[_s_nJetAK4]
  Float_t         JetAK4_phi[njetmx];   //[_s_nJetAK4]
  Float_t         JetAK4_mass[njetmx];   //[_s_nJetAK4]
  Float_t         JetAK4_btag_DeepCSV[njetmx];   //[_s_nJetAK4]
  Float_t         JetAK4_btag_DeepFlav[njetmx];   //[_s_nJetAK4]
  Int_t           JetAK4_hadronflav[njetmx];   //[_s_nJetAK4]
  Int_t           JetAK4_partonflav[njetmx];   //[_s_nJetAK4]
  Float_t         JetAK4_qgl[njetmx];   //[_s_nJetAK4]
  Float_t         JetAK4_PUID[njetmx];   //[_s_nJetAK4]
  Float_t         JetAK4_JESup[njetmx];   //[_s_nJetAK4]
  Float_t         JetAK4_JESdn[njetmx];   //[_s_nJetAK4]
  Float_t         JetAK4_JERup[njetmx];   //[_s_nJetAK4]
  Float_t         JetAK4_JERdn[njetmx];   //[_s_nJetAK4]
  Float_t         JetAK4_btag_DeepFlav_SF[njetmx];   //[_s_nJetAK4]
  Float_t         JetAK4_btag_DeepFlav_SF_up[njetmx];   //[_s_nJetAK4]
  Float_t         JetAK4_btag_DeepFlav_SF_dn[njetmx];   //[_s_nJetAK4]
  Int_t           nGenLep;
  Float_t         GenLep_pt[narray];   //[nGenLep]
  Float_t         GenLep_eta[narray];   //[nGenLep]
  Float_t         GenLep_phi[narray];   //[nGenLep]
  Float_t         GenLep_mass[narray];   //[nGenLep]
  Int_t           GenLep_pdgId[narray];   //[nGenLep]
  Int_t           GenLep_mompdgId[narray];   //[nGenLep]
  Int_t           GenLep_grmompdgId[narray];   //[nGenLep]
  Int_t           nGenNu;
  Float_t         GenNu_pt[narray];   //[nGenNu]
  Float_t         GenNu_eta[narray];   //[nGenNu]
  Float_t         GenNu_phi[narray];   //[nGenNu]
  Float_t         GenNu_mass[narray];   //[nGenNu]
  Int_t           GenNu_pdgId[narray];   //[nGenNu]
  Int_t           GenNu_mompdgId[narray];   //[nGenNu]
  Int_t           GenNu_grmompdgId[narray];   //[nGenNu]
  Int_t           nGenBPart;
  Float_t         GenBPart_pt[narray];   //[nGenBPart]
  Float_t         GenBPart_eta[narray];   //[nGenBPart]
  Float_t         GenBPart_phi[narray];   //[nGenBPart]
  Float_t         GenBPart_mass[narray];   //[nGenBPart]
  Int_t           GenBPart_pdgId[narray];   //[nGenBPart]
  Int_t           GenBPart_mompdgId[narray];   //[nGenBPart]
  Int_t           GenBPart_grmompdgId[narray];   //[nGenNu]
  Int_t           nGenV;
  Float_t         GenV_pt[narray];   //[nGenV]
  Float_t         GenV_eta[narray];   //[nGenV]
  Float_t         GenV_phi[narray];   //[nGenV]
  Float_t         GenV_mass[narray];   //[nGenV]
  Int_t           GenV_pdgId[narray];   //[nGenV]
  Int_t           GenV_mompdgId[narray];   //[nGenV]
  Int_t           GenV_grmompdgId[narray];   //[nGenNu]
  Int_t 		  nLHETop;
  Float_t 		  LHETop_pt[narray];
  Float_t		  LHETop_eta[narray];  
  Float_t		  LHETop_phi[narray];
  Float_t		  LHETop_mass[narray];  
  Int_t 		  nGenTop;
  Float_t 		  GenTop_pt[narray];
  Float_t		  GenTop_eta[narray];  
  Float_t		  GenTop_phi[narray];  
  Float_t		  GenTop_mass[narray];  
  
  Int_t           nLHEScaleWeights;
  Float_t         LHEScaleWeights[9];   //[nLHEScaleWeights]
  Int_t           nLHEPDFWeights;
  Float_t         LHEPDFWeights[103];   //[nLHEPDFWeights]
  Int_t           nLHEAlpsWeights;
  Float_t         LHEAlpsWeights[3];   //[nLHEAlpsWeights]
  Int_t           nLHEPSWeights;
  Float_t         LHEPSWeights[8];   //[nLHEPSWeights]


  float EWK_cor;
  float QCD_cor;

  float ptedges[] = {20, 25, 30, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 175, 200, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1101, 1248, 1410, 1588, 1784, 2000};//, 2366, 2787, 3450};
  const int nptbins = sizeof(ptedges)/sizeof(ptedges[0])-1;

  //float msdbins[] = {30,45,65,90,120,160,205,255,310,370,430,500,600}; //roughly 3 sigma bin width
  //float msdbins[] = {30,45,60,75,90,105,120,135,150,165,185,210,250,300,600}; //my bin
  float msdbins_semilep[]  = {30,50,70,90,120,150,170,190,210,260,600};
  float msdbins_dilep[]    = {30,50,70,90,110,130,160,200,600};
  float invmassbins_semilep[] =  {500,650,800,950,1150,1300,1500,1800,2100,2500,4000};
  float invmassbins_dilep[]   =  {500,550,600,650,725,850,1000,1250,1500,4000};
  //const int nunrollbins = 8*ninvmassbins;
  
  float HEM_weight;
  
  float yptbins[] = {200, 300, 450, 600, 3000};
  const int nyptbins  = sizeof(yptbins)/sizeof(yptbins[0])-1;
  
  int njecmax = 0;

  // change between SL & DL //
  TString rgn[] = {"SR1","SR2","CR2","CR3","CR4","CR5","CR6","CR7","CR8","QCDCR1","QCDCR2","QCDVR1","QCDVR2","QCDVR3","QCDVR4","VRTT"};
  //TString rgn[] = {"SR1","SR2","CR2","CR3","CR4","CR6","CR8","CRVjL","ARVjL","CRVjM","ARVjM","CRVjT","ARVjT","VRVj","VRTT"};

  int nrgn = sizeof(rgn)/sizeof(rgn[0]);
  
  TString rgn_CR[] = {"CR3_nb0","CR2_nb0","CR4_nb0","CR6_nb1"};
  int nCR = sizeof(rgn_CR)/sizeof(rgn_CR[0]);

  TString Ytype[] = {"T", "M" , "L"};
  //int nYtype = sizeof(Ytype)/sizeof(Ytype[0]);
  int y_wp = 0;

  TString Wtype[] = {"L", "M" , "T"};
  //int nWtype = sizeof(Wtype)/sizeof(Wtype[0]);
  int w_wp = 2;
  
  TString bcats[] = {"","_nb0","_nb1"};
  int nbcat = sizeof(bcats)/sizeof(bcats[0]);
  
  TString Wops[] = {"opt2"};  //,"opt1"};
  int nWop = sizeof(Wops)/sizeof(Wops[0]);
 
  TString tops[] = {"","_Top_fullymerged","_Top_semimerged","_Top_unmerged"};
  int ntop = sizeof(tops)/sizeof(tops[0]);
 
  TString lepids[] = {"","_Mu","_El","_EMu"};
  int nlid = sizeof(lepids)/sizeof(lepids[0]);
  
  TString sysnames[] = {
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
	 "PU","LeptonSF","LeptonSF2","Prefire","PNbbSF","PNWSF","BTG","TrigSF1","TrigSF2",
	 "CR_SF"
	 }; 
	 
  int nsys = sizeof(sysnames)/sizeof(sysnames[0]);
  
  bool isDL = false;
  bool useBout = true;
  bool isDATA = true;
  bool isSignal = false;
  
  int YEAR = 2018;
  
  // Scale factors derived in the analysis //
  
  //float DY_SF_DL[][nyptbins]          = {{1.0338326,0.8479128,0.8218863,0.7610330}, {1.00842,0.821746,0.798714,0.715957}, {1.06785,0.883104,0.861406,0.844771}, {1.0338326,0.8479128,0.8218863,0.7610330}};
  //float DY_SF_DL_uncs[][nyptbins]     = {{0.0270255,0.0259449,0.0336987,0.0421685}, {0.0255998,0.0239224,0.0347324,0.0462003}, {0.0388015,0.0403951,0.0526211,0.0702999}, {0.0270255,0.0259449,0.0336987,0.0421685}};
  float DY_SF_DL[][nyptbins]          = {{0.868614,0.772042,0.708655,0.616403}, {0.846313,0.74796,0.690209,0.586556}, {0.901945,0.80394,0.740487,0.670203}, {0.868614,0.772042,0.708655,0.616403}};
  float DY_SF_DL_uncs[][nyptbins]     = {{0.023041,0.0232153,0.0290738,0.0342703}, {0.0214107,0.0215712,0.0300266,0.0379569}, {0.0324458,0.0353412,0.0452637,0.0557093}, {0.023041,0.0232153,0.0290738,0.0342703}};
  
  // fitting tt+st per pt bin in CR2 //
  float Top_um_SF_DL[][nyptbins]      = {{0.9630,0.9850,0.9676,0.7616}, {0.9977,0.9932,0.9436,0.7799}, {0.8949,0.9833,1.0133,0.816}, {0.9744,0.9796,0.9637,0.671}};
  float Top_um_SF_DL_uncs[][nyptbins] = {{0.0245,0.0299,0.0499,0.0688}, {0.0291,0.0399,0.0715,0.0966}, {0.0356,0.0511,0.0963,0.137}, {0.0283,0.0380,0.0808,0.128}};
  // fitting tt+st inclusive in pT in CR2 //
  float Top_um_SF_DL_inc[] = {0.965,0.987,0.921,0.972};
  // fitting tt+st inclusive in pT in CR2+CR6 //
  //float Top_um_SF_DL_inc[] = {0.966,0.989,0.921,0.972};
  
  // fitting tt+st (separately) inclusive in pT in CR2 //
  float Top_um_SF_incl_TT[] = {0.987,1.091,0.876,0.829};
  float Top_um_SF_incl_ST[] = {0.457,0.000,0.000,2.660};
  // fitting tt+st (separately) inclusive in pT in CR2+CR6 //
  //float Top_um_SF_incl_TT[] = {0.926,1.093,0.873,0.829};
  //float Top_um_SF_incl_ST[] = {1.152,0.000,0.002,2.660};
  
  // fitting only tt per pt bin in CR2 //
  float TT_um_SF_DL_pt[][nyptbins] =       {{0.960,0.984,0.964,0.721},{0.998,0.993,0.936,0.741},{0.888,0.983,1.017,0.779},{0.972,0.978,0.959,0.631}};
  // fitting only tt per pt bin in CR2+CR6 //
  //float TT_um_SF_DL[][nyptbins] =       {{0.960,0.984,0.964,0.721},{0.998,0.993,0.936,0.741},{0.888,0.983,1.017,0.779},{0.972,0.978,0.959,0.631}};
  // fitting tt inclusive in pT in CR2 //
  float TT_um_SF_DL_inc[] = {0.962,0.987,0.915,0.970};
  // fitting tt inclusive in pT in CR2+CR6 //
  //float TT_um_SF_DL[] = {0.963,0.988,0.916,0.970};
  
  void read_branches(TTree *tree, bool isDL=false)
  {
	  
   tree->SetBranchAddress("irun", &irun);	
   tree->SetBranchAddress("ilumi", &ilumi);	
   //tree->SetBranchAddress("ievt", &ievt);	
   tree->SetBranchAddress("npvert", &npvert);	
   tree->SetBranchAddress("PV_npvsGood", &PV_npvsGood);	
	
   tree->SetBranchAddress("nleptons", &nleptons);
   tree->SetBranchAddress("nfatjets", &nfatjets);
   tree->SetBranchAddress("ncuts", &ncuts);
   tree->SetBranchAddress("Flag_event_cuts", Flag_event_cuts);

   tree->SetBranchAddress("Flag_pass_baseline",&Flag_pass_baseline);
   tree->SetBranchAddress("Flag_pass_baseline_no_LJet",&Flag_pass_baseline_no_LJet);

   tree->SetBranchAddress("hlt_IsoMu24", &hlt_IsoMu24);
   tree->SetBranchAddress("hlt_Mu50", &hlt_Mu50);
   tree->SetBranchAddress("hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", &hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
   tree->SetBranchAddress("hlt_Ele115_CaloIdVT_GsfTrkIdT", &hlt_Ele115_CaloIdVT_GsfTrkIdT);
   tree->SetBranchAddress("hlt_Ele40_WPTight_Gsf", &hlt_Ele40_WPTight_Gsf);
   tree->SetBranchAddress("hlt_Ele32_WPTight_Gsf", &hlt_Ele32_WPTight_Gsf);
   tree->SetBranchAddress("hlt_Ele28_eta2p1_WPTight_Gsf_HT150", &hlt_Ele28_eta2p1_WPTight_Gsf_HT150);
   tree->SetBranchAddress("hlt_Mu37_Ele27_CaloIdL_MW", &hlt_Mu37_Ele27_CaloIdL_MW);
   tree->SetBranchAddress("hlt_Mu27_Ele37_CaloIdL_MW", &hlt_Mu27_Ele37_CaloIdL_MW);
   tree->SetBranchAddress("hlt_Mu37_TkMu27", &hlt_Mu37_TkMu27);
   tree->SetBranchAddress("hlt_DoubleEle25_CaloIdL_MW", &hlt_DoubleEle25_CaloIdL_MW);
   tree->SetBranchAddress("hlt_AK8PFJet500", &hlt_AK8PFJet500);
   tree->SetBranchAddress("hlt_PFJet500", &hlt_PFJet500);
   tree->SetBranchAddress("hlt_HT1050", &hlt_HT1050);
   tree->SetBranchAddress("hlt_AK8PFJet400_TrimMass30", &hlt_AK8PFJet400_TrimMass30);
   tree->SetBranchAddress("hlt_AK8PFHT800_TrimMass50", &hlt_AK8PFHT800_TrimMass50);
   tree->SetBranchAddress("Muon_trig_pass", &Muon_trig_pass);
   tree->SetBranchAddress("Electron_trig_pass", &Electron_trig_pass);
   tree->SetBranchAddress("MuonElectron_trig_pass", &MuonElectron_trig_pass);

   tree->SetBranchAddress("l1_pt", &l1_pt);
   tree->SetBranchAddress("l1_eta", &l1_eta);
   tree->SetBranchAddress("l1_phi", &l1_phi);
   tree->SetBranchAddress("l1_mass", &l1_mass);
   tree->SetBranchAddress("l1_pdgId", &l1_pdgId);
   tree->SetBranchAddress("l1_minisoch", &l1_minisoch);
   tree->SetBranchAddress("l1_minisonh", &l1_minisonh);
   tree->SetBranchAddress("l1_minisoph", &l1_minisoph);
   tree->SetBranchAddress("l1_minisoall", &l1_minisoall);
   tree->SetBranchAddress("l1_genindex", &l1_genindex);
   if(isDL){
   tree->SetBranchAddress("l2_pt", &l2_pt);
   tree->SetBranchAddress("l2_eta", &l2_eta);
   tree->SetBranchAddress("l2_phi", &l2_phi);
   tree->SetBranchAddress("l2_mass", &l2_mass);
   tree->SetBranchAddress("l2_pdgId", &l2_pdgId);
   tree->SetBranchAddress("l2_minisoch", &l2_minisoch);
   tree->SetBranchAddress("l2_minisonh", &l2_minisonh);
   tree->SetBranchAddress("l2_minisoph", &l2_minisoph);
   tree->SetBranchAddress("l2_minisoall", &l2_minisoall);
   tree->SetBranchAddress("l2_genindex", &l2_genindex);
   }
   
   tree->SetBranchAddress("MET_pt", &MET_pt);
   tree->SetBranchAddress("MET_phi", &MET_phi);
   tree->SetBranchAddress("MET_sig", &MET_sig);
   tree->SetBranchAddress("MET_sumEt", &MET_sumEt);
   tree->SetBranchAddress("MET_pt_JESup", &MET_pt_JESup);
   tree->SetBranchAddress("MET_pt_JESdn", &MET_pt_JESdn);
   tree->SetBranchAddress("MET_pt_JERup", &MET_pt_JERup);
   tree->SetBranchAddress("MET_pt_JERdn", &MET_pt_JERdn);
   tree->SetBranchAddress("MET_pt_UnclusEup", &MET_pt_UnclusEup);
   tree->SetBranchAddress("MET_pt_UnclusEdn", &MET_pt_UnclusEdn);
   tree->SetBranchAddress("MET_phi_JESup", &MET_phi_JESup);
   tree->SetBranchAddress("MET_phi_JESdn", &MET_phi_JESdn);
   tree->SetBranchAddress("MET_phi_JERup", &MET_phi_JERup);
   tree->SetBranchAddress("MET_phi_JERdn", &MET_phi_JERdn);
   tree->SetBranchAddress("MET_phi_UnclusEup", &MET_phi_UnclusEup);
   tree->SetBranchAddress("MET_phi_UnclusEdn", &MET_phi_UnclusEdn);
   tree->SetBranchAddress("MET_pt_JESup_split", &MET_pt_JESup_split);
   tree->SetBranchAddress("MET_pt_JESdn_split", &MET_pt_JESdn_split);
   tree->SetBranchAddress("MET_phi_JESup_split", &MET_phi_JESup_split);
   tree->SetBranchAddress("MET_phi_JESdn_split", &MET_phi_JESdn_split);
    
 
   tree->SetBranchAddress("Y_pt", &Y_pt);
   tree->SetBranchAddress("Y_y", &Y_y);
   tree->SetBranchAddress("Y_eta", &Y_eta);
   tree->SetBranchAddress("Y_phi", &Y_phi);
   tree->SetBranchAddress("Y_mass", &Y_mass);
   tree->SetBranchAddress("Y_msoftdrop", &Y_msoftdrop);
   tree->SetBranchAddress("Y_tau21", &Y_tau21);
   tree->SetBranchAddress("Y_tau32", &Y_tau32);
   tree->SetBranchAddress("Y_DeepTag_DAK8MD_TvsQCD", &Y_DeepTag_DAK8MD_TvsQCD);
   tree->SetBranchAddress("Y_DeepTag_DAK8MD_WvsQCD", &Y_DeepTag_DAK8MD_WvsQCD);
   tree->SetBranchAddress("Y_DeepTag_DAK8MD_ZvsQCD", &Y_DeepTag_DAK8MD_ZvsQCD);
   tree->SetBranchAddress("Y_DeepTag_DAK8MD_HvsQCD", &Y_DeepTag_DAK8MD_HvsQCD);
   tree->SetBranchAddress("Y_DeepTag_DAK8MD_bbvsQCD", &Y_DeepTag_DAK8MD_bbvsQCD);
   tree->SetBranchAddress("Y_DeepTag_PNet_TvsQCD", &Y_DeepTag_PNet_TvsQCD);
   tree->SetBranchAddress("Y_DeepTag_PNet_WvsQCD", &Y_DeepTag_PNet_WvsQCD);
   tree->SetBranchAddress("Y_DeepTag_PNet_ZvsQCD", &Y_DeepTag_PNet_ZvsQCD);
   tree->SetBranchAddress("Y_DeepTag_PNetMD_XbbvsQCD", &Y_DeepTag_PNetMD_XbbvsQCD);
   tree->SetBranchAddress("Y_DeepTag_PNetMD_XccvsQCD", &Y_DeepTag_PNetMD_XccvsQCD);
   tree->SetBranchAddress("Y_DeepTag_PNetMD_XqqvsQCD", &Y_DeepTag_PNetMD_XqqvsQCD);
   tree->SetBranchAddress("Y_DeepTag_PNetMD_WvsQCD", &Y_DeepTag_PNetMD_WvsQCD);
   tree->SetBranchAddress("Y_DeepTag_PNetMD_QCD", &Y_DeepTag_PNetMD_QCD);
   tree->SetBranchAddress("Y_PN_bb", &Y_PN_bb);
   tree->SetBranchAddress("Y_label_Top_bq", &Y_label_Top_bq);
   tree->SetBranchAddress("Y_label_Top_bc", &Y_label_Top_bc);
   tree->SetBranchAddress("Y_label_Top_bcq", &Y_label_Top_bcq);
   tree->SetBranchAddress("Y_label_Top_bqq", &Y_label_Top_bqq);
   tree->SetBranchAddress("Y_label_W_qq", &Y_label_W_qq);
   tree->SetBranchAddress("Y_label_W_cq", &Y_label_W_cq);
   tree->SetBranchAddress("Y_sub1_pt", &Y_sub1_pt);
   tree->SetBranchAddress("Y_sub1_eta", &Y_sub1_eta);
   tree->SetBranchAddress("Y_sub1_phi", &Y_sub1_phi);
   tree->SetBranchAddress("Y_sub1_mass", &Y_sub1_mass);
   tree->SetBranchAddress("Y_sub1_btag", &Y_sub1_btag);
   tree->SetBranchAddress("Y_sub2_pt", &Y_sub2_pt);
   tree->SetBranchAddress("Y_sub2_eta", &Y_sub2_eta);
   tree->SetBranchAddress("Y_sub2_phi", &Y_sub2_phi);
   tree->SetBranchAddress("Y_sub2_mass", &Y_sub2_mass);
   tree->SetBranchAddress("Y_sub2_btag", &Y_sub2_btag);
   tree->SetBranchAddress("Y_genindex", &Y_genindex);
   tree->SetBranchAddress("Y_genbindex", Y_genbindex);
   tree->SetBranchAddress("Y_JESup", &Y_JESup);
   tree->SetBranchAddress("Y_JESdn", &Y_JESdn);
   tree->SetBranchAddress("Y_JERup", &Y_JERup);
   tree->SetBranchAddress("Y_JERdn", &Y_JERdn);
   tree->SetBranchAddress("Y_JESup_split", &Y_JESup_split);
   tree->SetBranchAddress("Y_JESdn_split", &Y_JESdn_split);
   if(!isDL){
   
   tree->SetBranchAddress("W_pt_opt1", &W_pt_opt1);
   tree->SetBranchAddress("W_y_opt1", &W_y_opt1);
   tree->SetBranchAddress("W_eta_opt1", &W_eta_opt1);
   tree->SetBranchAddress("W_phi_opt1", &W_phi_opt1);
   tree->SetBranchAddress("W_mass_opt1", &W_mass_opt1);
   tree->SetBranchAddress("W_msoftdrop_opt1", &W_msoftdrop_opt1);
   tree->SetBranchAddress("W_tau21_opt1", &W_tau21_opt1);
   tree->SetBranchAddress("W_tau32_opt1", &W_tau32_opt1);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_TvsQCD_opt1", &W_DeepTag_DAK8MD_TvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_WvsQCD_opt1", &W_DeepTag_DAK8MD_WvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_ZvsQCD_opt1", &W_DeepTag_DAK8MD_ZvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_HvsQCD_opt1", &W_DeepTag_DAK8MD_HvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_bbvsQCD_opt1", &W_DeepTag_DAK8MD_bbvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_PNet_TvsQCD_opt1", &W_DeepTag_PNet_TvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_PNet_WvsQCD_opt1", &W_DeepTag_PNet_WvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_PNet_ZvsQCD_opt1", &W_DeepTag_PNet_ZvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_PNetMD_XbbvsQCD_opt1", &W_DeepTag_PNetMD_XbbvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_PNetMD_XccvsQCD_opt1", &W_DeepTag_PNetMD_XccvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_PNetMD_XqqvsQCD_opt1", &W_DeepTag_PNetMD_XqqvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_PNetMD_WvsQCD_opt1", &W_DeepTag_PNetMD_WvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_PNetMD_QCD_opt1", &W_DeepTag_PNetMD_QCD_opt1);
   tree->SetBranchAddress("W_DAK8_W_opt1", &W_DAK8_W_opt1);
   tree->SetBranchAddress("W_PN_W_opt1", &W_PN_W_opt1);
   tree->SetBranchAddress("W_label_W_qq_opt1", &W_label_W_qq_opt1);
   tree->SetBranchAddress("W_label_W_cq_opt1", &W_label_W_cq_opt1);
   tree->SetBranchAddress("W_sub1_pt_opt1", &W_sub1_pt_opt1);
   tree->SetBranchAddress("W_sub1_eta_opt1", &W_sub1_eta_opt1);
   tree->SetBranchAddress("W_sub1_phi_opt1", &W_sub1_phi_opt1);
   tree->SetBranchAddress("W_sub1_mass_opt1", &W_sub1_mass_opt1);
   tree->SetBranchAddress("W_sub1_btag_opt1", &W_sub1_btag_opt1);
   tree->SetBranchAddress("W_sub2_pt_opt1", &W_sub2_pt_opt1);
   tree->SetBranchAddress("W_sub2_eta_opt1", &W_sub2_eta_opt1);
   tree->SetBranchAddress("W_sub2_phi_opt1", &W_sub2_phi_opt1);
   tree->SetBranchAddress("W_sub2_mass_opt1", &W_sub2_mass_opt1);
   tree->SetBranchAddress("W_sub2_btag_opt1", &W_sub2_btag_opt1);
   tree->SetBranchAddress("W_genindex_opt1", &W_genindex_opt1);
   tree->SetBranchAddress("W_JESup_opt1", &W_JESup_opt1);
   tree->SetBranchAddress("W_JESdn_opt1", &W_JESdn_opt1);
   tree->SetBranchAddress("W_JERup_opt1", &W_JERup_opt1);
   tree->SetBranchAddress("W_JERdn_opt1", &W_JERdn_opt1);
   tree->SetBranchAddress("W_JESup_split_opt1", &W_JESup_split_opt1);
   tree->SetBranchAddress("W_JESdn_split_opt1", &W_JESdn_split_opt1);
   tree->SetBranchAddress("H_pt_opt1", &H_pt_opt1);
   tree->SetBranchAddress("H_y_opt1", &H_y_opt1);
   tree->SetBranchAddress("H_eta_opt1", &H_eta_opt1);
   tree->SetBranchAddress("H_phi_opt1", &H_phi_opt1);
   tree->SetBranchAddress("H_mass_opt1", &H_mass_opt1);
   tree->SetBranchAddress("H_genindex_opt1", &H_genindex_opt1);
   tree->SetBranchAddress("H_JESup_opt1", &H_JESup_opt1);
   tree->SetBranchAddress("H_JESdn_opt1", &H_JESdn_opt1);
   tree->SetBranchAddress("H_JERup_opt1", &H_JERup_opt1);
   tree->SetBranchAddress("H_JERdn_opt1", &H_JERdn_opt1);
   tree->SetBranchAddress("H_JESup_split_opt1", &H_JESup_split_opt1);
   tree->SetBranchAddress("H_JESdn_split_opt1", &H_JESdn_split_opt1);
   tree->SetBranchAddress("X_mass_opt1", &X_mass_opt1);
   //tree->SetBranchAddress("X_mass_JESup_split_opt1", &X_mass_JESup_split_opt1);
   //tree->SetBranchAddress("X_mass_JESdn_split_opt1", &X_mass_JESdn_split_opt1);
   tree->SetBranchAddress("dR_lW_opt1", &dR_lW_opt1);
   tree->SetBranchAddress("dy_lW_opt1", &dy_lW_opt1);
   tree->SetBranchAddress("dphi_lW_opt1", &dphi_lW_opt1);
   
   tree->SetBranchAddress("W_pt_opt2", &W_pt_opt2);
   tree->SetBranchAddress("W_y_opt2", &W_y_opt2);
   tree->SetBranchAddress("W_eta_opt2", &W_eta_opt2);
   tree->SetBranchAddress("W_phi_opt2", &W_phi_opt2);
   tree->SetBranchAddress("W_mass_opt2", &W_mass_opt2);
   tree->SetBranchAddress("W_msoftdrop_opt2", &W_msoftdrop_opt2);
   tree->SetBranchAddress("W_tau21_opt2", &W_tau21_opt2);
   tree->SetBranchAddress("W_tau32_opt2", &W_tau32_opt2);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_TvsQCD_opt2", &W_DeepTag_DAK8MD_TvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_WvsQCD_opt2", &W_DeepTag_DAK8MD_WvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_ZvsQCD_opt2", &W_DeepTag_DAK8MD_ZvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_HvsQCD_opt2", &W_DeepTag_DAK8MD_HvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_bbvsQCD_opt2", &W_DeepTag_DAK8MD_bbvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_PNet_TvsQCD_opt2", &W_DeepTag_PNet_TvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_PNet_WvsQCD_opt2", &W_DeepTag_PNet_WvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_PNet_ZvsQCD_opt2", &W_DeepTag_PNet_ZvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_PNetMD_XbbvsQCD_opt2", &W_DeepTag_PNetMD_XbbvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_PNetMD_XccvsQCD_opt2", &W_DeepTag_PNetMD_XccvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_PNetMD_XqqvsQCD_opt2", &W_DeepTag_PNetMD_XqqvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_PNetMD_WvsQCD_opt2", &W_DeepTag_PNetMD_WvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_PNetMD_QCD_opt2", &W_DeepTag_PNetMD_QCD_opt2);
   tree->SetBranchAddress("W_DAK8_W_opt2", &W_DAK8_W_opt2);
   tree->SetBranchAddress("W_PN_W_opt2", &W_PN_W_opt2);
   tree->SetBranchAddress("W_label_W_qq_opt2", &W_label_W_qq_opt2);
   tree->SetBranchAddress("W_label_W_cq_opt2", &W_label_W_cq_opt2);
   tree->SetBranchAddress("W_sub1_pt_opt2", &W_sub1_pt_opt2);
   tree->SetBranchAddress("W_sub1_eta_opt2", &W_sub1_eta_opt2);
   tree->SetBranchAddress("W_sub1_phi_opt2", &W_sub1_phi_opt2);
   tree->SetBranchAddress("W_sub1_mass_opt2", &W_sub1_mass_opt2);
   tree->SetBranchAddress("W_sub1_btag_opt2", &W_sub1_btag_opt2);
   tree->SetBranchAddress("W_sub2_pt_opt2", &W_sub2_pt_opt2);
   tree->SetBranchAddress("W_sub2_eta_opt2", &W_sub2_eta_opt2);
   tree->SetBranchAddress("W_sub2_phi_opt2", &W_sub2_phi_opt2);
   tree->SetBranchAddress("W_sub2_mass_opt2", &W_sub2_mass_opt2);
   tree->SetBranchAddress("W_sub2_btag_opt2", &W_sub2_btag_opt2);
   tree->SetBranchAddress("W_genindex_opt2", &W_genindex_opt2);
   tree->SetBranchAddress("W_JESup_opt2", &W_JESup_opt2);
   tree->SetBranchAddress("W_JESdn_opt2", &W_JESdn_opt2);
   tree->SetBranchAddress("W_JERup_opt2", &W_JERup_opt2);
   tree->SetBranchAddress("W_JERdn_opt2", &W_JERdn_opt2);
   tree->SetBranchAddress("W_JESup_split_opt2", &W_JESup_split_opt2);
   tree->SetBranchAddress("W_JESdn_split_opt2", &W_JESdn_split_opt2);
   tree->SetBranchAddress("H_pt_opt2", &H_pt_opt2);
   tree->SetBranchAddress("H_y_opt2", &H_y_opt2);
   tree->SetBranchAddress("H_eta_opt2", &H_eta_opt2);
   tree->SetBranchAddress("H_phi_opt2", &H_phi_opt2);
   tree->SetBranchAddress("H_mass_opt2", &H_mass_opt2);
   tree->SetBranchAddress("H_genindex_opt2", &H_genindex_opt2);
   tree->SetBranchAddress("H_JESup_opt2", &H_JESup_opt2);
   tree->SetBranchAddress("H_JESdn_opt2", &H_JESdn_opt2);
   tree->SetBranchAddress("H_JERup_opt2", &H_JERup_opt2);
   tree->SetBranchAddress("H_JERdn_opt2", &H_JERdn_opt2);
   tree->SetBranchAddress("H_JESup_split_opt2", &H_JESup_split_opt2);
   tree->SetBranchAddress("H_JESdn_split_opt2", &H_JESdn_split_opt2);
   tree->SetBranchAddress("X_mass_JESup_split_opt2", &X_mass_JESup_split_opt2);
   tree->SetBranchAddress("X_mass_JESdn_split_opt2", &X_mass_JESdn_split_opt2);
   tree->SetBranchAddress("X_mass_opt2", &X_mass_opt2);
   tree->SetBranchAddress("dR_lW_opt2", &dR_lW_opt2);
   tree->SetBranchAddress("dy_lW_opt2", &dy_lW_opt2);
   tree->SetBranchAddress("dphi_lW_opt2", &dphi_lW_opt2);
   
   }
  
   tree->SetBranchAddress("dR_l1Y", &dR_l1Y);
   tree->SetBranchAddress("dy_l1Y", &dy_l1Y);
   tree->SetBranchAddress("dphi_l1Y", &dphi_l1Y);
   
   if(isDL){
   
   tree->SetBranchAddress("dR_l2Y", &dR_l2Y);
   tree->SetBranchAddress("dy_l2Y", &dy_l2Y);
   tree->SetBranchAddress("dphi_l2Y", &dphi_l2Y);
   
   tree->SetBranchAddress("l1l2_mass", &l1l2_mass);
   tree->SetBranchAddress("l1l2_deta", &l1l2_deta);
   tree->SetBranchAddress("l1l2_dphi", &l1l2_dphi);
   tree->SetBranchAddress("l1l2_dR", &l1l2_dR);
   tree->SetBranchAddress("dphi_MET_l1l2", &dphi_MET_l1l2);
   
   }
   
   tree->SetBranchAddress("HTlep_pt", &HTlep_pt);
   tree->SetBranchAddress("HTlep_pt_JESup", &HTlep_pt_JESup);
   tree->SetBranchAddress("HTlep_pt_JESdn", &HTlep_pt_JESdn);
   tree->SetBranchAddress("HTlep_pt_JERup", &HTlep_pt_JERup);
   tree->SetBranchAddress("HTlep_pt_JERdn", &HTlep_pt_JERdn);
   tree->SetBranchAddress("HTlep_pt_JESup_split", &HTlep_pt_JESup_split);
   tree->SetBranchAddress("HTlep_pt_JESdn_split", &HTlep_pt_JESdn_split);
   tree->SetBranchAddress("ST", &ST);
   tree->SetBranchAddress("ST_JESup", &ST_JESup);
   tree->SetBranchAddress("ST_JESdn", &ST_JESdn);
   tree->SetBranchAddress("ST_JERup", &ST_JERup);
   tree->SetBranchAddress("ST_JERdn", &ST_JERdn);
   tree->SetBranchAddress("ST_JESup_split", &ST_JESup_split);
   tree->SetBranchAddress("ST_JESdn_split", &ST_JESdn_split);
   
   tree->SetBranchAddress("nbjets_other", &nbjets_other);
   tree->SetBranchAddress("nbjets_outY", &nbjets_outY);
   tree->SetBranchAddress("nbjets_outY_L", &nbjets_outY_L);
   tree->SetBranchAddress("nbjets", &nbjets);
   tree->SetBranchAddress("nbjets_L", &nbjets_L);
   
   tree->SetBranchAddress("Flag_Y_bb_pass_T", &Flag_Y_bb_pass_T);
   tree->SetBranchAddress("Flag_Y_bb_pass_M", &Flag_Y_bb_pass_M);
   tree->SetBranchAddress("Flag_Y_bb_pass_L", &Flag_Y_bb_pass_L);
   if(!isDL){
   tree->SetBranchAddress("Flag_H_W_pass_T_opt1", &Flag_H_W_pass_T_opt1);
   tree->SetBranchAddress("Flag_H_W_pass_M_opt1", &Flag_H_W_pass_M_opt1);
   tree->SetBranchAddress("Flag_H_W_pass_L_opt1", &Flag_H_W_pass_L_opt1);
   tree->SetBranchAddress("Flag_H_m_pass_opt1", &Flag_H_m_pass_opt1);
   tree->SetBranchAddress("Flag_H_m_pass_opt2", &Flag_H_m_pass_opt2);
   tree->SetBranchAddress("Flag_dR_lW_pass_opt1", &Flag_dR_lW_pass_opt1);
   tree->SetBranchAddress("Flag_H_W_pass_T_opt2", &Flag_H_W_pass_T_opt2);
   tree->SetBranchAddress("Flag_H_W_pass_M_opt2", &Flag_H_W_pass_M_opt2);
   tree->SetBranchAddress("Flag_H_W_pass_L_opt2", &Flag_H_W_pass_L_opt2);
   tree->SetBranchAddress("Flag_dR_lW_pass_opt2", &Flag_dR_lW_pass_opt2);
   }
   tree->SetBranchAddress("Flag_MET_pass", &Flag_MET_pass);
   
   tree->SetBranchAddress("nPFJetAK8", &nPFJetAK8);
   tree->SetBranchAddress("PFJetAK8_pt", PFJetAK8_pt);
   tree->SetBranchAddress("PFJetAK8_eta", PFJetAK8_eta);
   tree->SetBranchAddress("PFJetAK8_phi", PFJetAK8_phi);
   tree->SetBranchAddress("PFJetAK8_mass", PFJetAK8_mass);
   tree->SetBranchAddress("PFJetAK8_jetID", PFJetAK8_jetID);
   tree->SetBranchAddress("PFJetAK8_jetID_tightlepveto", PFJetAK8_jetID_tightlepveto);
   tree->SetBranchAddress("PFJetAK8_msoftdrop", PFJetAK8_msoftdrop);
   tree->SetBranchAddress("PFJetAK8_tau21", PFJetAK8_tau21);
   tree->SetBranchAddress("PFJetAK8_tau32", PFJetAK8_tau32);
   tree->SetBranchAddress("PFJetAK8_DeepTag_PNetMD_XbbvsQCD", PFJetAK8_DeepTag_PNetMD_XbbvsQCD);
   tree->SetBranchAddress("PFJetAK8_DeepTag_PNetMD_WvsQCD", PFJetAK8_DeepTag_PNetMD_WvsQCD);
   tree->SetBranchAddress("PFJetAK8_DeepTag_PNet_TvsQCD", PFJetAK8_DeepTag_PNet_TvsQCD);
   tree->SetBranchAddress("PFJetAK8_DeepTag_PNet_WvsQCD", PFJetAK8_DeepTag_PNet_WvsQCD);
   tree->SetBranchAddress("PFJetAK8_DeepTag_DAK8MD_TvsQCD", PFJetAK8_DeepTag_DAK8MD_TvsQCD);
   tree->SetBranchAddress("PFJetAK8_DeepTag_DAK8MD_WvsQCD", PFJetAK8_DeepTag_DAK8MD_WvsQCD);
   tree->SetBranchAddress("PFJetAK8_DeepTag_DAK8MD_bbvsQCD", PFJetAK8_DeepTag_DAK8MD_bbvsQCD);
   tree->SetBranchAddress("PFJetAK8_JESup", PFJetAK8_JESup);
   tree->SetBranchAddress("PFJetAK8_JESdn", PFJetAK8_JESdn);
   //tree->SetBranchAddress("PFJetAK8_JERup", PFJetAK8_JERup);
   //tree->SetBranchAddress("PFJetAK8_JERdn", PFJetAK8_JERdn);
   tree->SetBranchAddress("PFJetAK8_Y_index", &PFJetAK8_Y_index);
   tree->SetBranchAddress("PFJetAK8_W_index_opt1", &PFJetAK8_W_index_opt1);
   tree->SetBranchAddress("PFJetAK8_W_index_opt2", &PFJetAK8_W_index_opt2);
   tree->SetBranchAddress("nJetAK4", &nJetAK4);
   tree->SetBranchAddress("JetAK4_pt", JetAK4_pt);
   tree->SetBranchAddress("JetAK4_eta", JetAK4_eta);
   tree->SetBranchAddress("JetAK4_phi", JetAK4_phi);
   tree->SetBranchAddress("JetAK4_mass", JetAK4_mass);
   tree->SetBranchAddress("JetAK4_btag_DeepCSV", JetAK4_btag_DeepCSV);
   tree->SetBranchAddress("JetAK4_btag_DeepFlav", JetAK4_btag_DeepFlav);
   tree->SetBranchAddress("JetAK4_hadronflav", JetAK4_hadronflav);
   tree->SetBranchAddress("JetAK4_partonflav", JetAK4_partonflav);
   tree->SetBranchAddress("JetAK4_qgl", JetAK4_qgl);
   tree->SetBranchAddress("JetAK4_PUID", JetAK4_PUID);
   tree->SetBranchAddress("JetAK4_JESup", JetAK4_JESup);
   tree->SetBranchAddress("JetAK4_JESdn", JetAK4_JESdn);
   //tree->SetBranchAddress("JetAK4_JERup", JetAK4_JERup);
   //tree->SetBranchAddress("JetAK4_JERdn", JetAK4_JERdn);
   tree->SetBranchAddress("JetAK4_btag_DeepFlav_SF", JetAK4_btag_DeepFlav_SF);
   tree->SetBranchAddress("JetAK4_btag_DeepFlav_SF_up", JetAK4_btag_DeepFlav_SF_up);
   tree->SetBranchAddress("JetAK4_btag_DeepFlav_SF_dn", JetAK4_btag_DeepFlav_SF_dn);
   
   if(!isDATA){
	
	tree->SetBranchAddress("LHE_weight", &LHE_weight);
	tree->SetBranchAddress("Generator_weight", &Generator_weight);
	tree->SetBranchAddress("Event_weight", &Event_weight);
	
	tree->SetBranchAddress("puWeight", &puWeight);
	tree->SetBranchAddress("puWeightup", &puWeightup);
	tree->SetBranchAddress("puWeightdown", &puWeightdown);
	tree->SetBranchAddress("leptonsf_weight", &leptonsf_weight);
	tree->SetBranchAddress("leptonsf_weight_up", &leptonsf_weight_up);
	tree->SetBranchAddress("leptonsf_weight_dn", &leptonsf_weight_dn);
	tree->SetBranchAddress("leptonsf_weight_stat", &leptonsf_weight_stat);
	tree->SetBranchAddress("leptonsf_weight_syst", &leptonsf_weight_syst);
	tree->SetBranchAddress("prefiringweight", &prefiringweight);
	tree->SetBranchAddress("prefiringweightup", &prefiringweightup);
	tree->SetBranchAddress("prefiringweightdown", &prefiringweightdown);
   
	tree->SetBranchAddress("nGenLep", &nGenLep);
	tree->SetBranchAddress("GenLep_pt", GenLep_pt);
	tree->SetBranchAddress("GenLep_eta", GenLep_eta);
	tree->SetBranchAddress("GenLep_phi", GenLep_phi);
	tree->SetBranchAddress("GenLep_mass", GenLep_mass);
	tree->SetBranchAddress("GenLep_pdgId", GenLep_pdgId);
	tree->SetBranchAddress("GenLep_mompdgId", GenLep_mompdgId);
	tree->SetBranchAddress("GenLep_grmompdgId", GenLep_grmompdgId);
	tree->SetBranchAddress("nGenNu", &nGenNu);
	tree->SetBranchAddress("GenNu_pt", GenNu_pt);
	tree->SetBranchAddress("GenNu_eta", GenNu_eta);
	tree->SetBranchAddress("GenNu_phi", GenNu_phi);
	tree->SetBranchAddress("GenNu_mass", GenNu_mass);
	tree->SetBranchAddress("GenNu_pdgId", GenNu_pdgId);
	tree->SetBranchAddress("GenNu_mompdgId", GenNu_mompdgId);
	tree->SetBranchAddress("GenNu_grmompdgId", GenNu_grmompdgId);
	tree->SetBranchAddress("nGenBPart", &nGenBPart);
	tree->SetBranchAddress("GenBPart_pt", &GenBPart_pt);
	tree->SetBranchAddress("GenBPart_eta", &GenBPart_eta);
	tree->SetBranchAddress("GenBPart_phi", &GenBPart_phi);
	tree->SetBranchAddress("GenBPart_mass", &GenBPart_mass);
	tree->SetBranchAddress("GenBPart_pdgId", &GenBPart_pdgId);
	tree->SetBranchAddress("GenBPart_mompdgId", &GenBPart_mompdgId);
	tree->SetBranchAddress("GenBPart_grmompdgId", GenBPart_grmompdgId);
	tree->SetBranchAddress("nGenV", &nGenV);
	tree->SetBranchAddress("GenV_pt", GenV_pt);
	tree->SetBranchAddress("GenV_eta", GenV_eta);
	tree->SetBranchAddress("GenV_phi", GenV_phi);
	tree->SetBranchAddress("GenV_mass", GenV_mass);
	tree->SetBranchAddress("GenV_pdgId", GenV_pdgId);
	tree->SetBranchAddress("GenV_mompdgId", GenV_mompdgId);
	tree->SetBranchAddress("GenV_grmompdgId", GenV_grmompdgId);
	tree->SetBranchAddress("nLHETop", &nLHETop);
	tree->SetBranchAddress("LHETop_pt", LHETop_pt);
	tree->SetBranchAddress("LHETop_eta", LHETop_eta);
	tree->SetBranchAddress("LHETop_phi", LHETop_phi);
	tree->SetBranchAddress("LHETop_mass", LHETop_mass);
	tree->SetBranchAddress("nGenTop", &nGenTop);
	tree->SetBranchAddress("GenTop_pt", GenTop_pt);
	tree->SetBranchAddress("GenTop_eta", GenTop_eta);
	tree->SetBranchAddress("GenTop_phi", GenTop_phi);
	tree->SetBranchAddress("GenTop_mass", GenTop_mass);
	
	tree->SetBranchAddress("nLHEScaleWeights", &nLHEScaleWeights);
    tree->SetBranchAddress("LHEScaleWeights", LHEScaleWeights);
    tree->SetBranchAddress("nLHEPDFWeights", &nLHEPDFWeights);
    tree->SetBranchAddress("LHEPDFWeights", LHEPDFWeights);
    //tree->SetBranchAddress("nLHEAlpsWeights", &nLHEAlpsWeights);
    //tree->SetBranchAddress("LHEAlpsWeights", LHEAlpsWeights);
    tree->SetBranchAddress("nLHEPSWeights", &nLHEPSWeights);
    tree->SetBranchAddress("LHEPSWeights", LHEPSWeights);

    }
       
	  
  }
  /*
  float fit_val_pt_L = 0.0902;
  float fit_val_pt_M = 0.0743;
  float fit_val_pt_T = 0.0030;

  float fit_val_msd_L_p0 = 0.0765;
  float fit_val_msd_M_p0 = 0.0611;
  float fit_val_msd_T_p0 = 0.0031;
  
  float fit_val_msd_L_p1 = 0.000205;
  float fit_val_msd_M_p1 = 0.000207;
  float fit_val_msd_T_p1 = 0;

  float fit_val_msd_L_high = 0.06953;
  float fit_val_msd_M_high = 0.07375;
  float fit_val_msd_T_high = 0.018;
  */
  
  float fit_val_pt_L = 0.0897;
  float fit_val_pt_M = 0.0711;
  float fit_val_pt_T = 0.0022;

  float fit_val_msd_L_p0 = 0.0786;
  float fit_val_msd_M_p0 = 0.0607;
  float fit_val_msd_T_p0 = 0.0025;
  
  float fit_val_msd_L_p1 = 0.000173;
  float fit_val_msd_M_p1 = 0.000176;
  float fit_val_msd_T_p1 = 0;

  float fit_val_msd_L_high = 0.0802486;
  float fit_val_msd_M_high = 0.0706991;
  float fit_val_msd_T_high = 0.0260286;

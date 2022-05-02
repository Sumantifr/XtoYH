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

int get_Y_id(float  Y_msoftdrop ) {
                        int Y_id;
                        if ( Y_msoftdrop < 50.0  )                               Y_id = 0;
                        else if ( Y_msoftdrop >= 50.0 &&  Y_msoftdrop < 100.0 )  Y_id = 1;
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

TH1F* get_histo_symbin(TString Y_op_name, TString W_op_name, TString reg_name, TString bcat_name, TString Wops_name, TString lep_name, string var, string addText, int nbins, float low_edge, float up_edge)
{
char name[200];
sprintf(name,"h_Y_%s_W_%s_%s_%s_%s%s%s%s",Y_op_name.Data(),W_op_name.Data(),var.c_str(),Wops_name.Data(),reg_name.Data(),bcat_name.Data(),lep_name.Data(),addText.c_str());
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

float PN_Top_med = 0.8;
float deep_btag_cut = 0.2783; 

//int main()
void histomaker_VI()
{
Bool_t isDATA = false;
TString proc_Name[] = {
"DYBJetsToLL_M-50_Zpt-100to200_XtoYH.root",
"DYBJetsToLL_M-50_Zpt-200toInf_XtoYH.root",
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
"ST_s-channel_XtoYH.root",
"ST_t-channel_antitop_XtoYH.root",
"ST_t-channel_top_XtoYH.root",
"ST_tW_antitop_XtoYH.root",
"ST_tW_top_XtoYH.root",
"TTTo2L2Nu_XtoYH.root",
"TTToHadronic_XtoYH.root",
"TTToSemiLeptonic_XtoYH.root",
"WBJetsToLNu_Wpt-100to200_XtoYH.root",
"WBJetsToLNu_Wpt-200toInf_XtoYH.root",
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
"QCD_HT300to500_XtoYH.root",
"QCD_HT500to700_XtoYH.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_1000_MY_100_XtoYH_Nov_2021_v2.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_1500_MY_200_XtoYH_Nov_2021_v2.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2000_MY_200_XtoYH_Nov_2021_v2.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2400_MY_300_XtoYH_Nov_2021_v2.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_3000_MY_100_XtoYH_Nov_2021_v2.root",
"NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_3000_MY_500_XtoYH_Nov_2021_v2.root"
};
 int nproc = sizeof(proc_Name)/sizeof(proc_Name[0]);
 for (int ii=0;ii<nproc;ii++)
  {
  //Calling efficiency files 
   TFile *FEff = TFile::Open("Efficiency/Efficiency_"+ proc_Name[ii] );
   TH2F *h_AK4M_flv_b_eff = (TH2F*)FEff->Get("Efficiency_h_Ak4_b_flv_pass_M");
   TH2F *h_AK4M_flv_c_eff = (TH2F*)FEff->Get("Efficiency_h_Ak4_c_flv_pass_M");
   TH2F *h_AK4M_flv_l_eff = (TH2F*)FEff->Get("Efficiency_h_Ak4_l_flv_pass_M");
   TH2F *h_YtagT_eff = (TH2F*)FEff->Get("Efficiency_h_Ak8_DeepTag_PNetMD_XbbvsQCD_pass_T");
   TH2F *h_WtagT_eff = (TH2F*)FEff->Get("Efficiency_h_Ak8_DeepTag_PNetMD_WvsQCD_pass_T");



   std::cout << proc_Name[ii] << std::endl;
   TFile* final_file = TFile::Open("OUTPUTS/Histogram_"+proc_Name[ii], "RECREATE");  

   TFile *file = TFile::Open(proc_Name[ii]);
   TTree *tree = (TTree*)file->Get("Tout");

   int narray = 20;
   Int_t           nleptons;
   Int_t           nfatjets;
   Bool_t          Flag_event_cuts;
   Float_t         puWeight;
   Float_t         puWeightup;
   Float_t         puWeightdown;
   Float_t         leptonsf_weight;
   Float_t         leptonsf_weight_up;
   Float_t         leptonsf_weight_dn;
   Float_t         leptonsf_weight_stat;
   Float_t         leptonsf_weight_syst;
   Float_t         l_pt;
   Float_t         l_eta;
   Float_t         l_phi;
   Float_t         l_mass;
   Int_t           l_pdgId;
   Float_t         l_minisoch;
   Float_t         l_minisonh;
   Float_t         l_minisoph;
   Float_t         l_minisoall;
   Int_t           l_genindex;
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
   Float_t         MET_phi_JESup;
   Float_t         MET_phi_JESdn;
   Float_t         MET_phi_JERup;
   Float_t         MET_phi_JERdn;
   Float_t         MET_phi_UnclusEup;
   Float_t         MET_phi_UnclusEdn;
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
   Float_t         X_mass_opt1;
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
   Float_t         X_mass_opt2;
   Float_t         dR_lW_opt2;
   Float_t         dy_lW_opt2;
   Float_t         dphi_lW_opt2;
   Float_t         dR_lY;
   Float_t         dy_lY;
   Float_t         dphi_lY;
   Int_t           nbjets_other;
   Bool_t          Flag_Y_bb_pass_T;
   Bool_t          Flag_Y_bb_pass_M;
   Bool_t          Flag_Y_bb_pass_L;
   Bool_t          Flag_H_W_pass_T_opt1;
   Bool_t          Flag_H_W_pass_M_opt1;
   Bool_t          Flag_H_W_pass_L_opt1;
   Bool_t          Flag_H_W_pass_T_opt2;
   Bool_t          Flag_H_W_pass_M_opt2;
   Bool_t          Flag_H_W_pass_L_opt2;
   Bool_t          Flag_H_m_pass_opt1;
   Bool_t          Flag_H_m_pass_opt2;
   Bool_t          Flag_dR_lW_pass_opt1;
   Bool_t          Flag_dR_lW_pass_opt2;
   Bool_t          Flag_MET_pass;
   Bool_t          Reg_SR_opt1;
   Bool_t          Reg_Wj_CR_opt1;
   Bool_t          Reg_SR_opt2;
   Bool_t          Reg_Wj_CR_opt2;
   Double_t        LHE_weight;
   Double_t        Generator_weight;
   Double_t        Event_weight;
   Double_t        prefiringweight;
   Double_t        prefiringweightup;
   Double_t        prefiringweightdown;
   Int_t           nPFJetAK8;
   Float_t         PFJetAK8_pt[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_eta[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_phi[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_mass[narray];   //[_s_nPFJetAK8]
   Bool_t          PFJetAK8_jetID[narray];   //[_s_nPFJetAK8]
   Bool_t          PFJetAK8_jetID_tightlepveto[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_msoftdrop[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_tau21[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_tau32[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_PNetMD_XbbvsQCD[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_PNetMD_WvsQCD[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_PNet_TvsQCD[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_PNet_WvsQCD[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_DAK8MD_TvsQCD[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_DAK8MD_WvsQCD[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_DAK8MD_bbvsQCD[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_JESup[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_JESdn[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_JERup[narray];   //[_s_nPFJetAK8]
   Float_t         PFJetAK8_JERdn[narray];   //[_s_nPFJetAK8]
   Int_t           PFJetAK8_Y_index;
   Int_t           PFJetAK8_W_index_opt1;
   Int_t           PFJetAK8_W_index_opt2;
   Int_t           nJetAK4;
   Float_t         JetAK4_pt[6];   //[_s_nJetAK4]
   Float_t         JetAK4_eta[6];   //[_s_nJetAK4]
   Float_t         JetAK4_phi[6];   //[_s_nJetAK4]
   Float_t         JetAK4_mass[6];   //[_s_nJetAK4]
   Float_t         JetAK4_btag_DeepCSV[6];   //[_s_nJetAK4]
   Float_t         JetAK4_btag_DeepFlav[6];   //[_s_nJetAK4]
   Int_t           JetAK4_hadronflav[6];   //[_s_nJetAK4]
   Int_t           JetAK4_partonflav[6];   //[_s_nJetAK4]
   Float_t         JetAK4_qgl[6];   //[_s_nJetAK4]
   Float_t         JetAK4_PUID[6];   //[_s_nJetAK4]
   Float_t         JetAK4_JESup[6];   //[_s_nJetAK4]
   Float_t         JetAK4_JESdn[6];   //[_s_nJetAK4]
   Float_t         JetAK4_JERup[6];   //[_s_nJetAK4]
   Float_t         JetAK4_JERdn[6];   //[_s_nJetAK4]
   Float_t         JetAK4_btag_DeepFlav_SF[6];   //[_s_nJetAK4]
   Float_t         JetAK4_btag_DeepFlav_SF_up[6];   //[_s_nJetAK4]
   Float_t         JetAK4_btag_DeepFlav_SF_dn[6];   //[_s_nJetAK4]
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


  // List of branches
   TBranch        *b_nleptons;   //!
   TBranch        *b_nfatjets;   //!
   TBranch        *b_Flag_event_cuts;   //!
   TBranch        *b_puWeight;   //!
   TBranch        *b_puWeightup;   //!
   TBranch        *b_puWeightdown;   //!
   TBranch        *b_leptonsf_weight;   //!
   TBranch        *b_leptonsf_weight_up;   //!
   TBranch        *b_leptonsf_weight_dn;   //!
   TBranch        *b_leptonsf_weight_stat;   //!
   TBranch        *b_leptonsf_weight_syst;   //!
   TBranch        *b_l_pt;   //!
   TBranch        *b_l_eta;   //!
   TBranch        *b_l_phi;   //!
   TBranch        *b_l_mass;   //!
   TBranch        *b_l_pdgId;   //!
   TBranch        *b_l_minisoch;   //!
   TBranch        *b_l_minisonh;   //!
   TBranch        *b_l_minisoph;   //!
   TBranch        *b_l_minisoall;   //!
   TBranch        *b_l_genindex;   //!
   TBranch        *b_MET_pt;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_sig;   //!
   TBranch        *b_MET_sumEt;   //!
   TBranch        *b_MET_pt_JESup;   //!
   TBranch        *b_MET_pt_JESdn;   //!
   TBranch        *b_MET_pt_JERup;   //!
   TBranch        *b_MET_pt_JERdn;   //!
   TBranch        *b_MET_pt_UnclusEup;   //!
   TBranch        *b_MET_pt_UnclusEdn;   //!
   TBranch        *b_MET_phi_JESup;   //!
   TBranch        *b_MET_phi_JESdn;   //!
   TBranch        *b_MET_phi_JERup;   //!
   TBranch        *b_MET_phi_JERdn;   //!
   TBranch        *b_MET_phi_UnclusEup;   //!
   TBranch        *b_MET_phi_UnclusEdn;   //!
   TBranch        *b_Y_pt;   //!
   TBranch        *b_Y_y;   //!
   TBranch        *b_Y_eta;   //!
   TBranch        *b_Y_phi;   //!
   TBranch        *b_Y_mass;   //!
   TBranch        *b_Y_msoftdrop;   //!
   TBranch        *b_Y_tau21;   //!
   TBranch        *b_Y_tau32;   //!
   TBranch        *b_Y_DeepTag_DAK8MD_TvsQCD;   //!
   TBranch        *b_Y_DeepTag_DAK8MD_WvsQCD;   //!
   TBranch        *b_Y_DeepTag_DAK8MD_ZvsQCD;   //!
   TBranch        *b_Y_DeepTag_DAK8MD_HvsQCD;   //!
   TBranch        *b_Y_DeepTag_DAK8MD_bbvsQCD;   //!
   TBranch        *b_Y_DeepTag_PNet_TvsQCD;   //!
   TBranch        *b_Y_DeepTag_PNet_WvsQCD;   //!
   TBranch        *b_Y_DeepTag_PNet_ZvsQCD;   //!
   TBranch        *b_Y_DeepTag_PNetMD_XbbvsQCD;   //!
   TBranch        *b_Y_DeepTag_PNetMD_XccvsQCD;   //!
   TBranch        *b_Y_DeepTag_PNetMD_XqqvsQCD;   //!
   TBranch        *b_Y_DeepTag_PNetMD_WvsQCD;   //!
   TBranch        *b_Y_DeepTag_PNetMD_QCD;   //!
   TBranch        *b_Y_PN_bb;   //!
   TBranch        *b_Y_sub1_pt;   //!
   TBranch        *b_Y_sub1_eta;   //!
   TBranch        *b_Y_sub1_phi;   //!
   TBranch        *b_Y_sub1_mass;   //!
   TBranch        *b_Y_sub1_btag;   //!
   TBranch        *b_Y_sub2_pt;   //!
   TBranch        *b_Y_sub2_eta;   //!
   TBranch        *b_Y_sub2_phi;   //!
   TBranch        *b_Y_sub2_mass;   //!
   TBranch        *b_Y_sub2_btag;   //!
   TBranch        *b_Y_genindex;   //!
   TBranch        *b_Y_genbindex;   //!
   TBranch        *b_Y_JESup;   //!
   TBranch        *b_Y_JESdn;   //!
   TBranch        *b_Y_JERup;   //!
   TBranch        *b_Y_JERdn;   //!
   TBranch        *b_W_pt_opt1;   //!
   TBranch        *b_W_y_opt1;   //!
   TBranch        *b_W_eta_opt1;   //!
   TBranch        *b_W_phi_opt1;   //!
   TBranch        *b_W_mass_opt1;   //!
   TBranch        *b_W_msoftdrop_opt1;   //!
   TBranch        *b_W_tau21_opt1;   //!
   TBranch        *b_W_tau32_opt1;   //!
   TBranch        *b_W_DeepTag_DAK8MD_TvsQCD_opt1;   //!
   TBranch        *b_W_DeepTag_DAK8MD_WvsQCD_opt1;   //!
   TBranch        *b_W_DeepTag_DAK8MD_ZvsQCD_opt1;   //!
   TBranch        *b_W_DeepTag_DAK8MD_HvsQCD_opt1;   //!
   TBranch        *b_W_DeepTag_DAK8MD_bbvsQCD_opt1;   //!
   TBranch        *b_W_DeepTag_PNet_TvsQCD_opt1;   //!
   TBranch        *b_W_DeepTag_PNet_WvsQCD_opt1;   //!
   TBranch        *b_W_DeepTag_PNet_ZvsQCD_opt1;   //!
   TBranch        *b_W_DeepTag_PNetMD_XbbvsQCD_opt1;   //!
   TBranch        *b_W_DeepTag_PNetMD_XccvsQCD_opt1;   //!
   TBranch        *b_W_DeepTag_PNetMD_XqqvsQCD_opt1;   //!
   TBranch        *b_W_DeepTag_PNetMD_WvsQCD_opt1;   //!
   TBranch        *b_W_DeepTag_PNetMD_QCD_opt1;   //!
   TBranch        *b_W_DAK8_W_opt1;   //!
   TBranch        *b_W_PN_W_opt1;   //!
   TBranch        *b_W_sub1_pt_opt1;   //!
   TBranch        *b_W_sub1_eta_opt1;   //!
   TBranch        *b_W_sub1_phi_opt1;   //!
   TBranch        *b_W_sub1_mass_opt1;   //!
   TBranch        *b_W_sub1_btag_opt1;   //!
   TBranch        *b_W_sub2_pt_opt1;   //!
   TBranch        *b_W_sub2_eta_opt1;   //!
   TBranch        *b_W_sub2_phi_opt1;   //!
   TBranch        *b_W_sub2_mass_opt1;   //!
   TBranch        *b_W_sub2_btag_opt1;   //!
   TBranch        *b_W_genindex_opt1;   //!
   TBranch        *b_W_JESup_opt1;   //!
   TBranch        *b_W_JESdn_opt1;   //!
   TBranch        *b_W_JERup_opt1;   //!
   TBranch        *b_W_JERdn_opt1;   //!
   TBranch        *b_H_pt_opt1;   //!
   TBranch        *b_H_y_opt1;   //!
   TBranch        *b_H_eta_opt1;   //!
   TBranch        *b_H_phi_opt1;   //!
   TBranch        *b_H_mass_opt1;   //!
   TBranch        *b_H_genindex_opt1;   //!
   TBranch        *b_H_JESup_opt1;   //!
   TBranch        *b_H_JESdn_opt1;   //!
   TBranch        *b_H_JERup_opt1;   //!
   TBranch        *b_H_JERdn_opt1;   //!
   TBranch        *b_X_mass_opt1;   //!
   TBranch        *b_dR_lW_opt1;   //!
   TBranch        *b_dy_lW_opt1;   //!
   TBranch        *b_dphi_lW_opt1;   //!
   TBranch        *b_W_pt_opt2;   //!
   TBranch        *b_W_y_opt2;   //!
   TBranch        *b_W_eta_opt2;   //!
   TBranch        *b_W_phi_opt2;   //!
   TBranch        *b_W_mass_opt2;   //!
   TBranch        *b_W_msoftdrop_opt2;   //!
   TBranch        *b_W_tau21_opt2;   //!
   TBranch        *b_W_tau32_opt2;   //!
   TBranch        *b_W_DeepTag_DAK8MD_TvsQCD_opt2;   //!
   TBranch        *b_W_DeepTag_DAK8MD_WvsQCD_opt2;   //!
   TBranch        *b_W_DeepTag_DAK8MD_ZvsQCD_opt2;   //!
   TBranch        *b_W_DeepTag_DAK8MD_HvsQCD_opt2;   //!
   TBranch        *b_W_DeepTag_DAK8MD_bbvsQCD_opt2;   //!
   TBranch        *b_W_DeepTag_PNet_TvsQCD_opt2;   //!
   TBranch        *b_W_DeepTag_PNet_WvsQCD_opt2;   //!
   TBranch        *b_W_DeepTag_PNet_ZvsQCD_opt2;   //!
   TBranch        *b_W_DeepTag_PNetMD_XbbvsQCD_opt2;   //!
   TBranch        *b_W_DeepTag_PNetMD_XccvsQCD_opt2;   //!
   TBranch        *b_W_DeepTag_PNetMD_XqqvsQCD_opt2;   //!
   TBranch        *b_W_DeepTag_PNetMD_WvsQCD_opt2;   //!
   TBranch        *b_W_DeepTag_PNetMD_QCD_opt2;   //!
   TBranch        *b_W_DAK8_W_opt2;   //!
   TBranch        *b_W_PN_W_opt2;   //!
   TBranch        *b_W_sub1_pt_opt2;   //!
   TBranch        *b_W_sub1_eta_opt2;   //!
   TBranch        *b_W_sub1_phi_opt2;   //!
   TBranch        *b_W_sub1_mass_opt2;   //!
   TBranch        *b_W_sub1_btag_opt2;   //!
   TBranch        *b_W_sub2_pt_opt2;   //!
   TBranch        *b_W_sub2_eta_opt2;   //!
   TBranch        *b_W_sub2_phi_opt2;   //!
   TBranch        *b_W_sub2_mass_opt2;   //!
   TBranch        *b_W_sub2_btag_opt2;   //!
   TBranch        *b_W_genindex_opt2;   //!
   TBranch        *b_W_JESup_opt2;   //!
   TBranch        *b_W_JESdn_opt2;   //!
   TBranch        *b_W_JERup_opt2;   //!
   TBranch        *b_W_JERdn_opt2;   //!
   TBranch        *b_H_pt_opt2;   //!
   TBranch        *b_H_y_opt2;   //!
   TBranch        *b_H_eta_opt2;   //!
   TBranch        *b_H_phi_opt2;   //!
   TBranch        *b_H_mass_opt2;   //!
   TBranch        *b_H_genindex_opt2;   //!
   TBranch        *b_H_JESup_opt2;   //!
   TBranch        *b_H_JESdn_opt2;   //!
   TBranch        *b_H_JERup_opt2;   //!
   TBranch        *b_H_JERdn_opt2;   //!
   TBranch        *b_X_mass_opt2;   //!
   TBranch        *b_dR_lW_opt2;   //!
   TBranch        *b_dy_lW_opt2;   //!
   TBranch        *b_dphi_lW_opt2;   //!
   TBranch        *b_dR_lY;   //!
   TBranch        *b_dy_lY;   //!
   TBranch        *b_dphi_lY;   //!
   TBranch        *b_nbjets_other;   //!
   TBranch        *b_Flag_Y_bb_pass_T;   //!
   TBranch        *b_Flag_Y_bb_pass_M;   //!
   TBranch        *b_Flag_Y_bb_pass_L;   //!
   TBranch        *b_Flag_H_W_pass_T_opt1;   //!
   TBranch        *b_Flag_H_W_pass_M_opt1;   //!
   TBranch        *b_Flag_H_W_pass_L_opt1;   //!
   TBranch        *b_Flag_H_W_pass_T_opt2;   //!
   TBranch        *b_Flag_H_W_pass_M_opt2;   //!
   TBranch        *b_Flag_H_W_pass_L_opt2;   //!
   TBranch        *b_Flag_H_m_pass_opt1;   //!
   TBranch        *b_Flag_H_m_pass_opt2;   //!
   TBranch        *b_Flag_dR_lW_pass_opt1;   //!
   TBranch        *b_Flag_dR_lW_pass_opt2;   //!
   TBranch        *b_Flag_MET_pass;   //!
   TBranch        *b_Reg_SR_opt1;   //!
   TBranch        *b_Reg_Wj_CR_opt1;   //!
   TBranch        *b_Reg_SR_opt2;   //!
   TBranch        *b_Reg_Wj_CR_opt2;   //!
   TBranch        *b_LHE_weight;   //!
   TBranch        *b_Generator_weight;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_prefiringweight;   //!
   TBranch        *b_prefiringweightup;   //!
   TBranch        *b_prefiringweightdown;   //!
   TBranch        *b__s_nPFJetAK8;   //!
   TBranch        *b_PFJetAK8_pt;   //!
   TBranch        *b_PFJetAK8_eta;   //!
   TBranch        *b_PFJetAK8_phi;   //!
   TBranch        *b_PFJetAK8_mass;   //!
   TBranch        *b_PFJetAK8_jetID;   //!
   TBranch        *b_PFJetAK8_jetID_tightlepveto;   //!
   TBranch        *b_PFJetAK8_msoftdrop;   //!
   TBranch        *b_PFJetAK8_tau21;   //!
   TBranch        *b_PFJetAK8_tau32;   //!
   TBranch        *b_PFJetAK8_DeepTag_PNetMD_XbbvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_PNetMD_WvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_PNet_TvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_PNet_WvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_DAK8MD_TvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_DAK8MD_WvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_DAK8MD_bbvsQCD;   //!
   TBranch        *b_PFJetAK8_JESup;   //!
   TBranch        *b_PFJetAK8_JESdn;   //!
   TBranch        *b_PFJetAK8_JERup;   //!
   TBranch        *b_PFJetAK8_JERdn;   //!
   TBranch        *b__s_PFJetAK8_Y_index;   //!
   TBranch        *b__s_PFJetAK8_W_index_opt1;   //!
   TBranch        *b__s_PFJetAK8_W_index_opt2;   //!
   TBranch        *b__s_nJetAK4;   //!
   TBranch        *b_JetAK4_pt;   //!
   TBranch        *b_JetAK4_eta;   //!
   TBranch        *b_JetAK4_phi;   //!
   TBranch        *b_JetAK4_mass;   //!
   TBranch        *b_JetAK4_btag_DeepCSV;   //!
   TBranch        *b_JetAK4_btag_DeepFlav;   //!
   TBranch        *b_JetAK4_hadronflav;   //!
   TBranch        *b_JetAK4_partonflav;   //!
   TBranch        *b_JetAK4_qgl;   //!
   TBranch        *b_JetAK4_PUID;   //!
   TBranch        *b_JetAK4_JESup;   //!
   TBranch        *b_JetAK4_JESdn;   //!
   TBranch        *b_JetAK4_JERup;   //!
   TBranch        *b_JetAK4_JERdn;   //!
   TBranch        *b_JetAK4_btag_DeepFlav_SF;   //!
   TBranch        *b_JetAK4_btag_DeepFlav_SF_up;   //!
   TBranch        *b_JetAK4_btag_DeepFlav_SF_dn;   //!
   TBranch        *b_nGenLep;   //!
   TBranch        *b_GenLep_pt;   //!
   TBranch        *b_GenLep_eta;   //!
   TBranch        *b_GenLep_phi;   //!
   TBranch        *b_GenLep_mass;   //!
   TBranch        *b_GenLep_pdgId;   //!
   TBranch        *b_GenLep_mompdgId;   //!
   TBranch        *b_GenLep_grmompdgId;   //!
   TBranch        *b_nGenNu;   //!
   TBranch        *b_GenNu_pt;   //!
   TBranch        *b_GenNu_eta;   //!
   TBranch        *b_GenNu_phi;   //!
   TBranch        *b_GenNu_mass;   //!
   TBranch        *b_GenNu_pdgId;   //!
   TBranch        *b_GenNu_mompdgId;   //!
   TBranch        *b_GenNu_grmompdgId;   //!
   TBranch        *b_nGenBPart;   //!
   TBranch        *b_GenBPart_pt;   //!
   TBranch        *b_GenBPart_eta;   //!
   TBranch        *b_GenBPart_phi;   //!
   TBranch        *b_GenBPart_mass;   //!
   TBranch        *b_GenBPart_pdgId;   //!
   TBranch        *b_GenBPart_mompdgId;   //!
   TBranch        *b_GenBPart_grmompdgId;   //!
   TBranch        *b_nGenV;   //!
   TBranch        *b_GenV_pt;   //!
   TBranch        *b_GenV_eta;   //!
   TBranch        *b_GenV_phi;   //!
   TBranch        *b_GenV_mass;   //!
   TBranch        *b_GenV_pdgId;   //!
   TBranch        *b_GenV_mompdgId;   //!
   TBranch        *b_GenV_grmompdgId;   //!



   Double_t        event_weight_LHE = 1.0;

   tree->SetBranchAddress("nleptons", &nleptons);
   tree->SetBranchAddress("nfatjets", &nfatjets);
   tree->SetBranchAddress("Flag_event_cuts", &Flag_event_cuts);
   tree->SetBranchAddress("puWeight", &puWeight);
   tree->SetBranchAddress("puWeightup", &puWeightup);
   tree->SetBranchAddress("puWeightdown", &puWeightdown);
   tree->SetBranchAddress("leptonsf_weight", &leptonsf_weight);
   tree->SetBranchAddress("leptonsf_weight_up", &leptonsf_weight_up, &b_leptonsf_weight_up);
   tree->SetBranchAddress("leptonsf_weight_dn", &leptonsf_weight_dn, &b_leptonsf_weight_dn);
   tree->SetBranchAddress("leptonsf_weight_stat", &leptonsf_weight_stat, &b_leptonsf_weight_stat);
   tree->SetBranchAddress("leptonsf_weight_syst", &leptonsf_weight_syst, &b_leptonsf_weight_syst);
   tree->SetBranchAddress("l_pt", &l_pt, &b_l_pt);
   tree->SetBranchAddress("l_eta", &l_eta, &b_l_eta);
   tree->SetBranchAddress("l_phi", &l_phi, &b_l_phi);
   tree->SetBranchAddress("l_mass", &l_mass, &b_l_mass);
   tree->SetBranchAddress("l_pdgId", &l_pdgId, &b_l_pdgId);
   tree->SetBranchAddress("l_minisoch", &l_minisoch, &b_l_minisoch);
   tree->SetBranchAddress("l_minisonh", &l_minisonh, &b_l_minisonh);
   tree->SetBranchAddress("l_minisoph", &l_minisoph, &b_l_minisoph);
   tree->SetBranchAddress("l_minisoall", &l_minisoall, &b_l_minisoall);
   tree->SetBranchAddress("l_genindex", &l_genindex, &b_l_genindex);
   tree->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
   tree->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   tree->SetBranchAddress("MET_sig", &MET_sig, &b_MET_sig);
   tree->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   tree->SetBranchAddress("MET_pt_JESup", &MET_pt_JESup, &b_MET_pt_JESup);
   tree->SetBranchAddress("MET_pt_JESdn", &MET_pt_JESdn, &b_MET_pt_JESdn);
   tree->SetBranchAddress("MET_pt_JERup", &MET_pt_JERup, &b_MET_pt_JERup);
   tree->SetBranchAddress("MET_pt_JERdn", &MET_pt_JERdn, &b_MET_pt_JERdn);
   tree->SetBranchAddress("MET_pt_UnclusEup", &MET_pt_UnclusEup, &b_MET_pt_UnclusEup);
   tree->SetBranchAddress("MET_pt_UnclusEdn", &MET_pt_UnclusEdn, &b_MET_pt_UnclusEdn);
   tree->SetBranchAddress("MET_phi_JESup", &MET_phi_JESup, &b_MET_phi_JESup);
   tree->SetBranchAddress("MET_phi_JESdn", &MET_phi_JESdn, &b_MET_phi_JESdn);
   tree->SetBranchAddress("MET_phi_JERup", &MET_phi_JERup, &b_MET_phi_JERup);
   tree->SetBranchAddress("MET_phi_JERdn", &MET_phi_JERdn, &b_MET_phi_JERdn);
   tree->SetBranchAddress("MET_phi_UnclusEup", &MET_phi_UnclusEup, &b_MET_phi_UnclusEup);
   tree->SetBranchAddress("MET_phi_UnclusEdn", &MET_phi_UnclusEdn, &b_MET_phi_UnclusEdn);
   tree->SetBranchAddress("Y_pt", &Y_pt, &b_Y_pt);
   tree->SetBranchAddress("Y_y", &Y_y, &b_Y_y);
   tree->SetBranchAddress("Y_eta", &Y_eta, &b_Y_eta);
   tree->SetBranchAddress("Y_phi", &Y_phi, &b_Y_phi);
   tree->SetBranchAddress("Y_mass", &Y_mass, &b_Y_mass);
   tree->SetBranchAddress("Y_msoftdrop", &Y_msoftdrop, &b_Y_msoftdrop);
   tree->SetBranchAddress("Y_tau21", &Y_tau21, &b_Y_tau21);
   tree->SetBranchAddress("Y_tau32", &Y_tau32, &b_Y_tau32);
   tree->SetBranchAddress("Y_DeepTag_DAK8MD_TvsQCD", &Y_DeepTag_DAK8MD_TvsQCD, &b_Y_DeepTag_DAK8MD_TvsQCD);
   tree->SetBranchAddress("Y_DeepTag_DAK8MD_WvsQCD", &Y_DeepTag_DAK8MD_WvsQCD, &b_Y_DeepTag_DAK8MD_WvsQCD);
   tree->SetBranchAddress("Y_DeepTag_DAK8MD_ZvsQCD", &Y_DeepTag_DAK8MD_ZvsQCD, &b_Y_DeepTag_DAK8MD_ZvsQCD);
   tree->SetBranchAddress("Y_DeepTag_DAK8MD_HvsQCD", &Y_DeepTag_DAK8MD_HvsQCD, &b_Y_DeepTag_DAK8MD_HvsQCD);
   tree->SetBranchAddress("Y_DeepTag_DAK8MD_bbvsQCD", &Y_DeepTag_DAK8MD_bbvsQCD, &b_Y_DeepTag_DAK8MD_bbvsQCD);
   tree->SetBranchAddress("Y_DeepTag_PNet_TvsQCD", &Y_DeepTag_PNet_TvsQCD, &b_Y_DeepTag_PNet_TvsQCD);
   tree->SetBranchAddress("Y_DeepTag_PNet_WvsQCD", &Y_DeepTag_PNet_WvsQCD, &b_Y_DeepTag_PNet_WvsQCD);
   tree->SetBranchAddress("Y_DeepTag_PNet_ZvsQCD", &Y_DeepTag_PNet_ZvsQCD, &b_Y_DeepTag_PNet_ZvsQCD);
   tree->SetBranchAddress("Y_DeepTag_PNetMD_XbbvsQCD", &Y_DeepTag_PNetMD_XbbvsQCD, &b_Y_DeepTag_PNetMD_XbbvsQCD);
   tree->SetBranchAddress("Y_DeepTag_PNetMD_XccvsQCD", &Y_DeepTag_PNetMD_XccvsQCD, &b_Y_DeepTag_PNetMD_XccvsQCD);
   tree->SetBranchAddress("Y_DeepTag_PNetMD_XqqvsQCD", &Y_DeepTag_PNetMD_XqqvsQCD, &b_Y_DeepTag_PNetMD_XqqvsQCD);
   tree->SetBranchAddress("Y_DeepTag_PNetMD_WvsQCD", &Y_DeepTag_PNetMD_WvsQCD, &b_Y_DeepTag_PNetMD_WvsQCD);
   tree->SetBranchAddress("Y_DeepTag_PNetMD_QCD", &Y_DeepTag_PNetMD_QCD, &b_Y_DeepTag_PNetMD_QCD);
   tree->SetBranchAddress("Y_PN_bb", &Y_PN_bb, &b_Y_PN_bb);
   tree->SetBranchAddress("Y_sub1_pt", &Y_sub1_pt, &b_Y_sub1_pt);
   tree->SetBranchAddress("Y_sub1_eta", &Y_sub1_eta, &b_Y_sub1_eta);
   tree->SetBranchAddress("Y_sub1_phi", &Y_sub1_phi, &b_Y_sub1_phi);
   tree->SetBranchAddress("Y_sub1_mass", &Y_sub1_mass, &b_Y_sub1_mass);
   tree->SetBranchAddress("Y_sub1_btag", &Y_sub1_btag, &b_Y_sub1_btag);
   tree->SetBranchAddress("Y_sub2_pt", &Y_sub2_pt, &b_Y_sub2_pt);
   tree->SetBranchAddress("Y_sub2_eta", &Y_sub2_eta, &b_Y_sub2_eta);
   tree->SetBranchAddress("Y_sub2_phi", &Y_sub2_phi, &b_Y_sub2_phi);
   tree->SetBranchAddress("Y_sub2_mass", &Y_sub2_mass, &b_Y_sub2_mass);
   tree->SetBranchAddress("Y_sub2_btag", &Y_sub2_btag, &b_Y_sub2_btag);
   tree->SetBranchAddress("Y_genindex", &Y_genindex, &b_Y_genindex);
   tree->SetBranchAddress("Y_genbindex", Y_genbindex, &b_Y_genbindex);
   tree->SetBranchAddress("Y_JESup", &Y_JESup, &b_Y_JESup);
   tree->SetBranchAddress("Y_JESdn", &Y_JESdn, &b_Y_JESdn);
   tree->SetBranchAddress("Y_JERup", &Y_JERup, &b_Y_JERup);
   tree->SetBranchAddress("Y_JERdn", &Y_JERdn, &b_Y_JERdn);
   tree->SetBranchAddress("W_pt_opt1", &W_pt_opt1, &b_W_pt_opt1);
   tree->SetBranchAddress("W_y_opt1", &W_y_opt1, &b_W_y_opt1);
   tree->SetBranchAddress("W_eta_opt1", &W_eta_opt1, &b_W_eta_opt1);
   tree->SetBranchAddress("W_phi_opt1", &W_phi_opt1, &b_W_phi_opt1);
   tree->SetBranchAddress("W_mass_opt1", &W_mass_opt1, &b_W_mass_opt1);
   tree->SetBranchAddress("W_msoftdrop_opt1", &W_msoftdrop_opt1, &b_W_msoftdrop_opt1);
   tree->SetBranchAddress("W_tau21_opt1", &W_tau21_opt1, &b_W_tau21_opt1);
   tree->SetBranchAddress("W_tau32_opt1", &W_tau32_opt1, &b_W_tau32_opt1);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_TvsQCD_opt1", &W_DeepTag_DAK8MD_TvsQCD_opt1, &b_W_DeepTag_DAK8MD_TvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_WvsQCD_opt1", &W_DeepTag_DAK8MD_WvsQCD_opt1, &b_W_DeepTag_DAK8MD_WvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_ZvsQCD_opt1", &W_DeepTag_DAK8MD_ZvsQCD_opt1, &b_W_DeepTag_DAK8MD_ZvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_HvsQCD_opt1", &W_DeepTag_DAK8MD_HvsQCD_opt1, &b_W_DeepTag_DAK8MD_HvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_bbvsQCD_opt1", &W_DeepTag_DAK8MD_bbvsQCD_opt1, &b_W_DeepTag_DAK8MD_bbvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_PNet_TvsQCD_opt1", &W_DeepTag_PNet_TvsQCD_opt1, &b_W_DeepTag_PNet_TvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_PNet_WvsQCD_opt1", &W_DeepTag_PNet_WvsQCD_opt1, &b_W_DeepTag_PNet_WvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_PNet_ZvsQCD_opt1", &W_DeepTag_PNet_ZvsQCD_opt1, &b_W_DeepTag_PNet_ZvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_PNetMD_XbbvsQCD_opt1", &W_DeepTag_PNetMD_XbbvsQCD_opt1, &b_W_DeepTag_PNetMD_XbbvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_PNetMD_XccvsQCD_opt1", &W_DeepTag_PNetMD_XccvsQCD_opt1, &b_W_DeepTag_PNetMD_XccvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_PNetMD_XqqvsQCD_opt1", &W_DeepTag_PNetMD_XqqvsQCD_opt1, &b_W_DeepTag_PNetMD_XqqvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_PNetMD_WvsQCD_opt1", &W_DeepTag_PNetMD_WvsQCD_opt1, &b_W_DeepTag_PNetMD_WvsQCD_opt1);
   tree->SetBranchAddress("W_DeepTag_PNetMD_QCD_opt1", &W_DeepTag_PNetMD_QCD_opt1, &b_W_DeepTag_PNetMD_QCD_opt1);
   tree->SetBranchAddress("W_DAK8_W_opt1", &W_DAK8_W_opt1, &b_W_DAK8_W_opt1);
   tree->SetBranchAddress("W_PN_W_opt1", &W_PN_W_opt1, &b_W_PN_W_opt1);
   tree->SetBranchAddress("W_sub1_pt_opt1", &W_sub1_pt_opt1, &b_W_sub1_pt_opt1);
   tree->SetBranchAddress("W_sub1_eta_opt1", &W_sub1_eta_opt1, &b_W_sub1_eta_opt1);
   tree->SetBranchAddress("W_sub1_phi_opt1", &W_sub1_phi_opt1, &b_W_sub1_phi_opt1);
   tree->SetBranchAddress("W_sub1_mass_opt1", &W_sub1_mass_opt1, &b_W_sub1_mass_opt1);
   tree->SetBranchAddress("W_sub1_btag_opt1", &W_sub1_btag_opt1, &b_W_sub1_btag_opt1);
   tree->SetBranchAddress("W_sub2_pt_opt1", &W_sub2_pt_opt1, &b_W_sub2_pt_opt1);
   tree->SetBranchAddress("W_sub2_eta_opt1", &W_sub2_eta_opt1, &b_W_sub2_eta_opt1);
   tree->SetBranchAddress("W_sub2_phi_opt1", &W_sub2_phi_opt1, &b_W_sub2_phi_opt1);
   tree->SetBranchAddress("W_sub2_mass_opt1", &W_sub2_mass_opt1, &b_W_sub2_mass_opt1);
   tree->SetBranchAddress("W_sub2_btag_opt1", &W_sub2_btag_opt1, &b_W_sub2_btag_opt1);
   tree->SetBranchAddress("W_genindex_opt1", &W_genindex_opt1, &b_W_genindex_opt1);
   tree->SetBranchAddress("W_JESup_opt1", &W_JESup_opt1, &b_W_JESup_opt1);
   tree->SetBranchAddress("W_JESdn_opt1", &W_JESdn_opt1, &b_W_JESdn_opt1);
   tree->SetBranchAddress("W_JERup_opt1", &W_JERup_opt1, &b_W_JERup_opt1);
   tree->SetBranchAddress("W_JERdn_opt1", &W_JERdn_opt1, &b_W_JERdn_opt1);
   tree->SetBranchAddress("H_pt_opt1", &H_pt_opt1, &b_H_pt_opt1);
   tree->SetBranchAddress("H_y_opt1", &H_y_opt1, &b_H_y_opt1);
   tree->SetBranchAddress("H_eta_opt1", &H_eta_opt1, &b_H_eta_opt1);
   tree->SetBranchAddress("H_phi_opt1", &H_phi_opt1, &b_H_phi_opt1);
   tree->SetBranchAddress("H_mass_opt1", &H_mass_opt1, &b_H_mass_opt1);
   tree->SetBranchAddress("H_genindex_opt1", &H_genindex_opt1, &b_H_genindex_opt1);
   tree->SetBranchAddress("H_JESup_opt1", &H_JESup_opt1, &b_H_JESup_opt1);
   tree->SetBranchAddress("H_JESdn_opt1", &H_JESdn_opt1, &b_H_JESdn_opt1);
   tree->SetBranchAddress("H_JERup_opt1", &H_JERup_opt1, &b_H_JERup_opt1);
   tree->SetBranchAddress("H_JERdn_opt1", &H_JERdn_opt1, &b_H_JERdn_opt1);
   tree->SetBranchAddress("X_mass_opt1", &X_mass_opt1, &b_X_mass_opt1);
   tree->SetBranchAddress("dR_lW_opt1", &dR_lW_opt1, &b_dR_lW_opt1);
   tree->SetBranchAddress("dy_lW_opt1", &dy_lW_opt1, &b_dy_lW_opt1);
   tree->SetBranchAddress("dphi_lW_opt1", &dphi_lW_opt1, &b_dphi_lW_opt1);
   tree->SetBranchAddress("W_pt_opt2", &W_pt_opt2, &b_W_pt_opt2);
   tree->SetBranchAddress("W_y_opt2", &W_y_opt2, &b_W_y_opt2);
   tree->SetBranchAddress("W_eta_opt2", &W_eta_opt2, &b_W_eta_opt2);
   tree->SetBranchAddress("W_phi_opt2", &W_phi_opt2, &b_W_phi_opt2);
   tree->SetBranchAddress("W_mass_opt2", &W_mass_opt2, &b_W_mass_opt2);
   tree->SetBranchAddress("W_msoftdrop_opt2", &W_msoftdrop_opt2, &b_W_msoftdrop_opt2);
   tree->SetBranchAddress("W_tau21_opt2", &W_tau21_opt2, &b_W_tau21_opt2);
   tree->SetBranchAddress("W_tau32_opt2", &W_tau32_opt2, &b_W_tau32_opt2);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_TvsQCD_opt2", &W_DeepTag_DAK8MD_TvsQCD_opt2, &b_W_DeepTag_DAK8MD_TvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_WvsQCD_opt2", &W_DeepTag_DAK8MD_WvsQCD_opt2, &b_W_DeepTag_DAK8MD_WvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_ZvsQCD_opt2", &W_DeepTag_DAK8MD_ZvsQCD_opt2, &b_W_DeepTag_DAK8MD_ZvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_HvsQCD_opt2", &W_DeepTag_DAK8MD_HvsQCD_opt2, &b_W_DeepTag_DAK8MD_HvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_DAK8MD_bbvsQCD_opt2", &W_DeepTag_DAK8MD_bbvsQCD_opt2, &b_W_DeepTag_DAK8MD_bbvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_PNet_TvsQCD_opt2", &W_DeepTag_PNet_TvsQCD_opt2, &b_W_DeepTag_PNet_TvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_PNet_WvsQCD_opt2", &W_DeepTag_PNet_WvsQCD_opt2, &b_W_DeepTag_PNet_WvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_PNet_ZvsQCD_opt2", &W_DeepTag_PNet_ZvsQCD_opt2, &b_W_DeepTag_PNet_ZvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_PNetMD_XbbvsQCD_opt2", &W_DeepTag_PNetMD_XbbvsQCD_opt2, &b_W_DeepTag_PNetMD_XbbvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_PNetMD_XccvsQCD_opt2", &W_DeepTag_PNetMD_XccvsQCD_opt2, &b_W_DeepTag_PNetMD_XccvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_PNetMD_XqqvsQCD_opt2", &W_DeepTag_PNetMD_XqqvsQCD_opt2, &b_W_DeepTag_PNetMD_XqqvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_PNetMD_WvsQCD_opt2", &W_DeepTag_PNetMD_WvsQCD_opt2, &b_W_DeepTag_PNetMD_WvsQCD_opt2);
   tree->SetBranchAddress("W_DeepTag_PNetMD_QCD_opt2", &W_DeepTag_PNetMD_QCD_opt2, &b_W_DeepTag_PNetMD_QCD_opt2);
   tree->SetBranchAddress("W_DAK8_W_opt2", &W_DAK8_W_opt2, &b_W_DAK8_W_opt2);
   tree->SetBranchAddress("W_PN_W_opt2", &W_PN_W_opt2, &b_W_PN_W_opt2);
   tree->SetBranchAddress("W_sub1_pt_opt2", &W_sub1_pt_opt2, &b_W_sub1_pt_opt2);
   tree->SetBranchAddress("W_sub1_eta_opt2", &W_sub1_eta_opt2, &b_W_sub1_eta_opt2);
   tree->SetBranchAddress("W_sub1_phi_opt2", &W_sub1_phi_opt2, &b_W_sub1_phi_opt2);
   tree->SetBranchAddress("W_sub1_mass_opt2", &W_sub1_mass_opt2, &b_W_sub1_mass_opt2);
   tree->SetBranchAddress("W_sub1_btag_opt2", &W_sub1_btag_opt2, &b_W_sub1_btag_opt2);
   tree->SetBranchAddress("W_sub2_pt_opt2", &W_sub2_pt_opt2, &b_W_sub2_pt_opt2);
   tree->SetBranchAddress("W_sub2_eta_opt2", &W_sub2_eta_opt2, &b_W_sub2_eta_opt2);
   tree->SetBranchAddress("W_sub2_phi_opt2", &W_sub2_phi_opt2, &b_W_sub2_phi_opt2);
   tree->SetBranchAddress("W_sub2_mass_opt2", &W_sub2_mass_opt2, &b_W_sub2_mass_opt2);
   tree->SetBranchAddress("W_sub2_btag_opt2", &W_sub2_btag_opt2, &b_W_sub2_btag_opt2);
   tree->SetBranchAddress("W_genindex_opt2", &W_genindex_opt2, &b_W_genindex_opt2);
   tree->SetBranchAddress("W_JESup_opt2", &W_JESup_opt2, &b_W_JESup_opt2);
   tree->SetBranchAddress("W_JESdn_opt2", &W_JESdn_opt2, &b_W_JESdn_opt2);
   tree->SetBranchAddress("W_JERup_opt2", &W_JERup_opt2, &b_W_JERup_opt2);
   tree->SetBranchAddress("W_JERdn_opt2", &W_JERdn_opt2, &b_W_JERdn_opt2);
   tree->SetBranchAddress("H_pt_opt2", &H_pt_opt2, &b_H_pt_opt2);
   tree->SetBranchAddress("H_y_opt2", &H_y_opt2, &b_H_y_opt2);
   tree->SetBranchAddress("H_eta_opt2", &H_eta_opt2, &b_H_eta_opt2);
   tree->SetBranchAddress("H_phi_opt2", &H_phi_opt2, &b_H_phi_opt2);
   tree->SetBranchAddress("H_mass_opt2", &H_mass_opt2, &b_H_mass_opt2);
   tree->SetBranchAddress("H_genindex_opt2", &H_genindex_opt2, &b_H_genindex_opt2);
   tree->SetBranchAddress("H_JESup_opt2", &H_JESup_opt2, &b_H_JESup_opt2);
   tree->SetBranchAddress("H_JESdn_opt2", &H_JESdn_opt2, &b_H_JESdn_opt2);
   tree->SetBranchAddress("H_JERup_opt2", &H_JERup_opt2, &b_H_JERup_opt2);
   tree->SetBranchAddress("H_JERdn_opt2", &H_JERdn_opt2, &b_H_JERdn_opt2);
   tree->SetBranchAddress("X_mass_opt2", &X_mass_opt2, &b_X_mass_opt2);
   tree->SetBranchAddress("dR_lW_opt2", &dR_lW_opt2, &b_dR_lW_opt2);
   tree->SetBranchAddress("dy_lW_opt2", &dy_lW_opt2, &b_dy_lW_opt2);
   tree->SetBranchAddress("dphi_lW_opt2", &dphi_lW_opt2, &b_dphi_lW_opt2);
   tree->SetBranchAddress("dR_lY", &dR_lY, &b_dR_lY);
   tree->SetBranchAddress("dy_lY", &dy_lY, &b_dy_lY);
   tree->SetBranchAddress("dphi_lY", &dphi_lY, &b_dphi_lY);
   tree->SetBranchAddress("nbjets_other", &nbjets_other, &b_nbjets_other);
   tree->SetBranchAddress("Flag_Y_bb_pass_T", &Flag_Y_bb_pass_T, &b_Flag_Y_bb_pass_T);
   tree->SetBranchAddress("Flag_Y_bb_pass_M", &Flag_Y_bb_pass_M, &b_Flag_Y_bb_pass_M);
   tree->SetBranchAddress("Flag_Y_bb_pass_L", &Flag_Y_bb_pass_L, &b_Flag_Y_bb_pass_L);
   tree->SetBranchAddress("Flag_H_W_pass_T_opt1", &Flag_H_W_pass_T_opt1, &b_Flag_H_W_pass_T_opt1);
   tree->SetBranchAddress("Flag_H_W_pass_M_opt1", &Flag_H_W_pass_M_opt1, &b_Flag_H_W_pass_M_opt1);
   tree->SetBranchAddress("Flag_H_W_pass_L_opt1", &Flag_H_W_pass_L_opt1, &b_Flag_H_W_pass_L_opt1);
   tree->SetBranchAddress("Flag_H_W_pass_T_opt2", &Flag_H_W_pass_T_opt2, &b_Flag_H_W_pass_T_opt2);
   tree->SetBranchAddress("Flag_H_W_pass_M_opt2", &Flag_H_W_pass_M_opt2, &b_Flag_H_W_pass_M_opt2);
   tree->SetBranchAddress("Flag_H_W_pass_L_opt2", &Flag_H_W_pass_L_opt2, &b_Flag_H_W_pass_L_opt2);
   tree->SetBranchAddress("Flag_H_m_pass_opt1", &Flag_H_m_pass_opt1, &b_Flag_H_m_pass_opt1);
   tree->SetBranchAddress("Flag_H_m_pass_opt2", &Flag_H_m_pass_opt2, &b_Flag_H_m_pass_opt2);
   tree->SetBranchAddress("Flag_dR_lW_pass_opt1", &Flag_dR_lW_pass_opt1, &b_Flag_dR_lW_pass_opt1);
   tree->SetBranchAddress("Flag_dR_lW_pass_opt2", &Flag_dR_lW_pass_opt2, &b_Flag_dR_lW_pass_opt2);
   tree->SetBranchAddress("Flag_MET_pass", &Flag_MET_pass, &b_Flag_MET_pass);
   tree->SetBranchAddress("Reg_SR_opt1", &Reg_SR_opt1, &b_Reg_SR_opt1);
   tree->SetBranchAddress("Reg_Wj_CR_opt1", &Reg_Wj_CR_opt1, &b_Reg_Wj_CR_opt1);
   tree->SetBranchAddress("Reg_SR_opt2", &Reg_SR_opt2, &b_Reg_SR_opt2);
   tree->SetBranchAddress("Reg_Wj_CR_opt2", &Reg_Wj_CR_opt2, &b_Reg_Wj_CR_opt2);
   tree->SetBranchAddress("LHE_weight", &LHE_weight, &b_LHE_weight);
   tree->SetBranchAddress("Generator_weight", &Generator_weight, &b_Generator_weight);
   tree->SetBranchAddress("Event_weight", &Event_weight, &b_weight);
   tree->SetBranchAddress("prefiringweight", &prefiringweight, &b_prefiringweight);
   tree->SetBranchAddress("prefiringweightup", &prefiringweightup, &b_prefiringweightup);
   tree->SetBranchAddress("prefiringweightdown", &prefiringweightdown, &b_prefiringweightdown);
   tree->SetBranchAddress("nPFJetAK8", &nPFJetAK8, &b__s_nPFJetAK8);
   tree->SetBranchAddress("PFJetAK8_pt", PFJetAK8_pt, &b_PFJetAK8_pt);
   tree->SetBranchAddress("PFJetAK8_eta", PFJetAK8_eta, &b_PFJetAK8_eta);
   tree->SetBranchAddress("PFJetAK8_phi", PFJetAK8_phi, &b_PFJetAK8_phi);
   tree->SetBranchAddress("PFJetAK8_mass", PFJetAK8_mass, &b_PFJetAK8_mass);
   tree->SetBranchAddress("PFJetAK8_jetID", PFJetAK8_jetID, &b_PFJetAK8_jetID);
   tree->SetBranchAddress("PFJetAK8_jetID_tightlepveto", PFJetAK8_jetID_tightlepveto, &b_PFJetAK8_jetID_tightlepveto);
   tree->SetBranchAddress("PFJetAK8_msoftdrop", PFJetAK8_msoftdrop, &b_PFJetAK8_msoftdrop);
   tree->SetBranchAddress("PFJetAK8_tau21", PFJetAK8_tau21, &b_PFJetAK8_tau21);
   tree->SetBranchAddress("PFJetAK8_tau32", PFJetAK8_tau32, &b_PFJetAK8_tau32);
   tree->SetBranchAddress("PFJetAK8_DeepTag_PNetMD_XbbvsQCD", PFJetAK8_DeepTag_PNetMD_XbbvsQCD, &b_PFJetAK8_DeepTag_PNetMD_XbbvsQCD);
   tree->SetBranchAddress("PFJetAK8_DeepTag_PNetMD_WvsQCD", PFJetAK8_DeepTag_PNetMD_WvsQCD, &b_PFJetAK8_DeepTag_PNetMD_WvsQCD);
   tree->SetBranchAddress("PFJetAK8_DeepTag_PNet_TvsQCD", PFJetAK8_DeepTag_PNet_TvsQCD, &b_PFJetAK8_DeepTag_PNet_TvsQCD);
   tree->SetBranchAddress("PFJetAK8_DeepTag_PNet_WvsQCD", PFJetAK8_DeepTag_PNet_WvsQCD, &b_PFJetAK8_DeepTag_PNet_WvsQCD);
   tree->SetBranchAddress("PFJetAK8_DeepTag_DAK8MD_TvsQCD", PFJetAK8_DeepTag_DAK8MD_TvsQCD, &b_PFJetAK8_DeepTag_DAK8MD_TvsQCD);
   tree->SetBranchAddress("PFJetAK8_DeepTag_DAK8MD_WvsQCD", PFJetAK8_DeepTag_DAK8MD_WvsQCD, &b_PFJetAK8_DeepTag_DAK8MD_WvsQCD);
   tree->SetBranchAddress("PFJetAK8_DeepTag_DAK8MD_bbvsQCD", PFJetAK8_DeepTag_DAK8MD_bbvsQCD, &b_PFJetAK8_DeepTag_DAK8MD_bbvsQCD);
   tree->SetBranchAddress("PFJetAK8_JESup", PFJetAK8_JESup, &b_PFJetAK8_JESup);
   tree->SetBranchAddress("PFJetAK8_JESdn", PFJetAK8_JESdn, &b_PFJetAK8_JESdn);
   tree->SetBranchAddress("PFJetAK8_JERup", PFJetAK8_JERup, &b_PFJetAK8_JERup);
   tree->SetBranchAddress("PFJetAK8_JERdn", PFJetAK8_JERdn, &b_PFJetAK8_JERdn);
   tree->SetBranchAddress("PFJetAK8_Y_index", &PFJetAK8_Y_index, &b__s_PFJetAK8_Y_index);
   tree->SetBranchAddress("PFJetAK8_W_index_opt1", &PFJetAK8_W_index_opt1, &b__s_PFJetAK8_W_index_opt1);
   tree->SetBranchAddress("PFJetAK8_W_index_opt2", &PFJetAK8_W_index_opt2, &b__s_PFJetAK8_W_index_opt2);
   tree->SetBranchAddress("nJetAK4", &nJetAK4, &b__s_nJetAK4);
   tree->SetBranchAddress("JetAK4_pt", JetAK4_pt, &b_JetAK4_pt);
   tree->SetBranchAddress("JetAK4_eta", JetAK4_eta, &b_JetAK4_eta);
   tree->SetBranchAddress("JetAK4_phi", JetAK4_phi, &b_JetAK4_phi);
   tree->SetBranchAddress("JetAK4_mass", JetAK4_mass, &b_JetAK4_mass);
   tree->SetBranchAddress("JetAK4_btag_DeepCSV", JetAK4_btag_DeepCSV, &b_JetAK4_btag_DeepCSV);
   tree->SetBranchAddress("JetAK4_btag_DeepFlav", JetAK4_btag_DeepFlav, &b_JetAK4_btag_DeepFlav);
   tree->SetBranchAddress("JetAK4_hadronflav", JetAK4_hadronflav, &b_JetAK4_hadronflav);
   tree->SetBranchAddress("JetAK4_partonflav", JetAK4_partonflav, &b_JetAK4_partonflav);
   tree->SetBranchAddress("JetAK4_qgl", JetAK4_qgl, &b_JetAK4_qgl);
   tree->SetBranchAddress("JetAK4_PUID", JetAK4_PUID, &b_JetAK4_PUID);
   tree->SetBranchAddress("JetAK4_JESup", JetAK4_JESup, &b_JetAK4_JESup);
   tree->SetBranchAddress("JetAK4_JESdn", JetAK4_JESdn, &b_JetAK4_JESdn);
   tree->SetBranchAddress("JetAK4_JERup", JetAK4_JERup, &b_JetAK4_JERup);
   tree->SetBranchAddress("JetAK4_JERdn", JetAK4_JERdn, &b_JetAK4_JERdn);
   tree->SetBranchAddress("JetAK4_btag_DeepFlav_SF", JetAK4_btag_DeepFlav_SF, &b_JetAK4_btag_DeepFlav_SF);
   tree->SetBranchAddress("JetAK4_btag_DeepFlav_SF_up", JetAK4_btag_DeepFlav_SF_up, &b_JetAK4_btag_DeepFlav_SF_up);
   tree->SetBranchAddress("JetAK4_btag_DeepFlav_SF_dn", JetAK4_btag_DeepFlav_SF_dn, &b_JetAK4_btag_DeepFlav_SF_dn);
   tree->SetBranchAddress("nGenLep", &nGenLep, &b_nGenLep);
   tree->SetBranchAddress("GenLep_pt", GenLep_pt, &b_GenLep_pt);
   tree->SetBranchAddress("GenLep_eta", GenLep_eta, &b_GenLep_eta);
   tree->SetBranchAddress("GenLep_phi", GenLep_phi, &b_GenLep_phi);
   tree->SetBranchAddress("GenLep_mass", GenLep_mass, &b_GenLep_mass);
   tree->SetBranchAddress("GenLep_pdgId", GenLep_pdgId, &b_GenLep_pdgId);
   tree->SetBranchAddress("GenLep_mompdgId", GenLep_mompdgId, &b_GenLep_mompdgId);
   tree->SetBranchAddress("GenLep_grmompdgId", GenLep_grmompdgId, &b_GenLep_grmompdgId);
   tree->SetBranchAddress("nGenNu", &nGenNu, &b_nGenNu);
   tree->SetBranchAddress("GenNu_pt", GenNu_pt, &b_GenNu_pt);
   tree->SetBranchAddress("GenNu_eta", GenNu_eta, &b_GenNu_eta);
   tree->SetBranchAddress("GenNu_phi", GenNu_phi, &b_GenNu_phi);
   tree->SetBranchAddress("GenNu_mass", GenNu_mass, &b_GenNu_mass);
   tree->SetBranchAddress("GenNu_pdgId", GenNu_pdgId, &b_GenNu_pdgId);
   tree->SetBranchAddress("GenNu_mompdgId", GenNu_mompdgId, &b_GenNu_mompdgId);
   tree->SetBranchAddress("GenNu_grmompdgId", GenNu_grmompdgId, &b_GenNu_grmompdgId);
   tree->SetBranchAddress("nGenBPart", &nGenBPart, &b_nGenBPart);
   tree->SetBranchAddress("GenBPart_pt", &GenBPart_pt, &b_GenBPart_pt);
   tree->SetBranchAddress("GenBPart_eta", &GenBPart_eta, &b_GenBPart_eta);
   tree->SetBranchAddress("GenBPart_phi", &GenBPart_phi, &b_GenBPart_phi);
   tree->SetBranchAddress("GenBPart_mass", &GenBPart_mass, &b_GenBPart_mass);
   tree->SetBranchAddress("GenBPart_pdgId", &GenBPart_pdgId, &b_GenBPart_pdgId);
   tree->SetBranchAddress("GenBPart_mompdgId", &GenBPart_mompdgId, &b_GenBPart_mompdgId);
   tree->SetBranchAddress("GenBPart_grmompdgId", GenBPart_grmompdgId, &b_GenBPart_grmompdgId);
   tree->SetBranchAddress("nGenV", &nGenV, &b_nGenV);
   tree->SetBranchAddress("GenV_pt", GenV_pt, &b_GenV_pt);
   tree->SetBranchAddress("GenV_eta", GenV_eta, &b_GenV_eta);
   tree->SetBranchAddress("GenV_phi", GenV_phi, &b_GenV_phi);
   tree->SetBranchAddress("GenV_mass", GenV_mass, &b_GenV_mass);
   tree->SetBranchAddress("GenV_pdgId", GenV_pdgId, &b_GenV_pdgId);
   tree->SetBranchAddress("GenV_mompdgId", GenV_mompdgId, &b_GenV_mompdgId);
   tree->SetBranchAddress("GenV_grmompdgId", GenV_grmompdgId, &b_GenV_grmompdgId);
   //tree->SetBranchAddress("event_weight_LHE", &event_weight_LHE);

  float ptedges[] = {20, 25, 30, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1101, 1248, 1410, 1588, 1784, 2000, 2366, 2787, 3450};
  const int nptbins = sizeof(ptedges)/sizeof(ptedges[0])-1;

  TString rgn[] = {"SR","CR1","CR2","CR3","CR4","CR5","CR6","CR7","CR8"};
  int nrgn = sizeof(rgn)/sizeof(rgn[0]);

  TString Ytype[] = {"T", "M" , "L"};
  //int nYtype = sizeof(Ytype)/sizeof(Ytype[0]);
  int y_wp = 0;

  TString Wtype[] = {"L", "M" , "T"};
  //int nWtype = sizeof(Wtype)/sizeof(Wtype[0]);
  int w_wp = 2;
  
  TString bcats[] = {"","_nb0","_nb1"};
  int nbcat = sizeof(bcats)/sizeof(bcats[0]);
  
  TString Wops[] = {"opt1","opt2"};
  int nWop = sizeof(Wops)/sizeof(Wops[0]);
  
  TString lepids[] = {"","_Mu","_El"};
  int nlid = sizeof(lepids)/sizeof(lepids[0]);
  
  TString sysnames[] = {"JES","JER","PU","LeptonSF","LeptonSF2","Prefire","PNbbSF","PNWSF","BTG"}; 
  int nsys = sizeof(sysnames)/sizeof(sysnames[0]);
  
  
  // Declaration of histograms //
  
  final_file->cd();
  
  TH1D* h_nom;
  TH1D* h_nom_reg[nrgn];  
  TH1D* h_reg = new TH1D("h_reg_id","h_reg_id",10,0.0,10.0);

  TH1F* h_MET_pt_ref[nrgn][nbcat][nWop][nlid];
  TH1F *h_Y_PNetMD_XbbvsQCD_ref[nrgn][nbcat][nWop][nlid];
  TH1F* h_W_msoftdrop_ref[nrgn][nbcat][nWop][nlid]; 
  
  TH1F* h_l_pt[nrgn][nbcat][nWop][nlid];
  TH1F* h_l_eta[nrgn][nbcat][nWop][nlid];
  TH1F* h_l_minisoall[nrgn][nbcat][nWop][nlid];
  
  TH1F* h_MET_pt[nrgn][nbcat][nWop][nlid];
  TH1F* h_MET_sig[nrgn][nbcat][nWop][nlid]; 
  TH1F* h_MT[nrgn][nbcat][nWop][nlid];  
  
  TH1F* h_leadbjet_pt[nrgn][nbcat][nWop][nlid];  
  TH1F* h_leadbjet_btag_DeepFlav[nrgn][nbcat][nWop][nlid];  
  
  TH1F *h_Y_pt[nrgn][nbcat][nWop][nlid];
  TH1F *h_Y_msoftdrop[nrgn][nbcat][nWop][nlid];
  TH1F *h_Y_PNetMD_XbbvsQCD[nrgn][nbcat][nWop][nlid];
  TH1F *h_Y_PNetMD_WvsQCD[nrgn][nbcat][nWop][nlid];
  TH1F *h_Y_PNet_TvsQCD[nrgn][nbcat][nWop][nlid];
  TH1F *h_Y_sub1_mass[nrgn][nbcat][nWop][nlid];
  TH1F *h_Y_sub2_mass[nrgn][nbcat][nWop][nlid];
  TH1F *h_Y_sub1_btag[nrgn][nbcat][nWop][nlid];
  TH1F *h_Y_sub2_btag[nrgn][nbcat][nWop][nlid];
		   
  TH1F* h_W_pt[nrgn][nbcat][nWop][nlid];
  TH1F* h_W_msoftdrop[nrgn][nbcat][nWop][nlid]; 
  TH1F* h_W_PNetMD_XbbvsQCD[nrgn][nbcat][nWop][nlid];
  TH1F* h_W_PNetMD_WvsQCD[nrgn][nbcat][nWop][nlid];
  TH1F* h_W_PNet_TvsQCD[nrgn][nbcat][nWop][nlid];
  TH1F* h_W_DAK8MD_WvsQCD[nrgn][nbcat][nWop][nlid];
  TH1F* h_W_sub1_mass[nrgn][nbcat][nWop][nlid];
  TH1F* h_W_sub2_mass[nrgn][nbcat][nWop][nlid];
  TH1F* h_W_sub1_btag[nrgn][nbcat][nWop][nlid];
  TH1F* h_W_sub2_btag[nrgn][nbcat][nWop][nlid];
  
  TH1F* h_dR_lW[nrgn][nbcat][nWop][nlid];
  TH1F* h_dy_lW[nrgn][nbcat][nWop][nlid];
  TH1F* h_dphi_lW[nrgn][nbcat][nWop][nlid];
				
  TH1F* h_H_mass[nrgn][nbcat][nWop][nlid]; 
  TH1F* h_X_mass[nrgn][nbcat][nWop][nlid]; 
  TH1F* h_nbjets_other[nrgn][nbcat][nWop][nlid]; 
  
  TH2F* h_X_Y_mass[nrgn][nbcat][nWop][nlid];   
  
  TH1F *h_Y_msoftdrop_sys[nrgn][nbcat][nWop][nlid][1+2*nsys];
  TH1F* h_X_mass_sys[nrgn][nbcat][nWop][nlid][1+2*nsys];

  TH1F *h_for_limit_X_mass[nrgn][nbcat][nWop][nlid];
  TH1F *h_for_limit_X_mass_sys[nrgn][nbcat][nWop][nlid][1+2*nsys];

  TH1F *h_for_limit_X_mass_v2[nrgn][nbcat][nWop][nlid];
  TH1F *h_for_limit_X_mass_sys_v2[nrgn][nbcat][nWop][nlid][1+2*nsys];

  // End of declaration //
  
  // Definition of histograms //
  
  h_nom = new TH1D("h_nom","h_nom",7,0.0,7.0);
  h_nom->Sumw2();
  
  for(int ij=0; ij<nrgn; ij++){
	char name[50];
	sprintf(name,"h_nom_%s",rgn[ij].Data());
	h_nom_reg[ij] = new TH1D(name,name,7,0.,7.);
	h_nom_reg[ij]->Sumw2();
  }
  
  for (int ij=0 ; ij< nrgn ; ij++)
  {
	  for(int jk=0; jk<nbcat; jk++)
	  {
		  for(int kl=0; kl<nWop; kl++)
		  {
			  for(int lm=0; lm<nlid; lm++)
			  {
				  
				//ref hist //
				
				h_MET_pt_ref[ij][jk][kl][lm] 			= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"MET_pt_ref","",40,0,1000);
				h_Y_PNetMD_XbbvsQCD_ref[ij][jk][kl][lm] = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_PNetMD_XbbvsQCD_ref","", 100, 0.0, 1.0 );
				h_W_msoftdrop_ref[ij][jk][kl][lm] 		 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"W_msoftdrop_ref","",50,0,350);
				
				// all hist//
                             
				h_l_pt[ij][jk][kl][lm] 				= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"l_pt","",40,0,1000);
				h_l_eta[ij][jk][kl][lm]				= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"l_eta","",50,-2.5,2.5);
				h_l_minisoall[ij][jk][kl][lm]		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"l_minisoall","",150,0,1.5);
	       	
				h_MET_pt[ij][jk][kl][lm] 			= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"MET_pt","",40,0,1000);
				h_MET_sig[ij][jk][kl][lm] 			= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"MET_sig","",50,0,300);				
				h_MT[ij][jk][kl][lm] 				= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"MT","",25,0,300);
			
				//h_leadbjet_pt[ij][jk][kl][lm] 		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"LeadBJet_pt","",40,0,1000);
				//h_leadbjet_btag_DeepFlav[ij][jk][kl][lm] = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"LeadBJet_btag_DeepFlav","",50,0,1);
				
				h_Y_pt[ij][jk][kl][lm] 				= get_histo_asymbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_pt","",nptbins, ptedges);
				h_Y_msoftdrop[ij][jk][kl][lm] 		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_msoftdrop","",40,0,600);
				h_Y_PNetMD_XbbvsQCD[ij][jk][kl][lm] = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_PNetMD_XbbvsQCD","", 100, 0.0, 1.0 );
				h_Y_PNetMD_WvsQCD[ij][jk][kl][lm] 	= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_PNetMD_WvsQCD","", 100, 0.0, 1.0 );
				h_Y_PNet_TvsQCD[ij][jk][kl][lm] 	= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_PNet_TvsQCD","", 100, 0.0, 1.0 );
				h_Y_sub1_mass[ij][jk][kl][lm] 		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_sub1_mass","", 40, 0.0, 300 );
				h_Y_sub2_mass[ij][jk][kl][lm] 		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_sub2_mass","", 40, 0.0, 300 );
				h_Y_sub1_btag[ij][jk][kl][lm] 		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_sub1_btag","", 50, 0.0, 1.0 );
				h_Y_sub2_btag[ij][jk][kl][lm] 		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_sub2_btag","", 50, 0.0, 1.0 );
				
				h_W_pt[ij][jk][kl][lm]	 			 = get_histo_asymbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"W_pt","",nptbins, ptedges);
				h_W_msoftdrop[ij][jk][kl][lm] 		 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"W_msoftdrop","",50,0,350);
				h_W_PNetMD_XbbvsQCD[ij][jk][kl][lm] = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"W_PNetMD_XbbvsQCD","",100,0,1);
				h_W_PNetMD_WvsQCD[ij][jk][kl][lm] 	 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"W_PNetMD_WvsQCD","",100,0,1);
				h_W_PNet_TvsQCD[ij][jk][kl][lm] 	 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"W_PNet_TvsQCD","",100,0,1);
				h_W_DAK8MD_WvsQCD[ij][jk][kl][lm]	 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"W_DAK8MD_WvsQCD","",100,0,1);
				h_W_sub1_mass[ij][jk][kl][lm] 		 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"W_sub1_mass","",40, 0.0, 300);
				h_W_sub2_mass[ij][jk][kl][lm] 		 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"W_sub2_mass","",40, 0.0, 300);
				h_W_sub1_btag[ij][jk][kl][lm] 		 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"W_sub1_btag","",50, 0.0, 1);
				h_W_sub2_btag[ij][jk][kl][lm] 		 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"W_sub2_btag","",50, 0.0, 1);
				
				h_dR_lW[ij][jk][kl][lm] 			 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"dR_lW","",100, 0.0, 5.0);
				h_dy_lW[ij][jk][kl][lm] 			 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"dy_lW","",100, -5.0, 5.0);
				h_dphi_lW[ij][jk][kl][lm]			 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"dphi_lW","",90, -M_PI, +M_PI);
				
				h_H_mass[ij][jk][kl][lm] 			 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"H_mass","",35, 0.0, 350.0);
				
				h_X_mass[ij][jk][kl][lm] 			 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"X_mass","",40, 0.0, 4000.0);
				
				h_nbjets_other[ij][jk][kl][lm]		 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"nbjets_other","",5 ,0.0, 5.0 );
               
				h_X_Y_mass[ij][jk][kl][lm] 		     = new TH2F("h_Y_"+Ytype[y_wp]+"_W_"+Wtype[w_wp]+"_X_Y_mass_"+Wops[kl]+"_"+rgn[ij]+bcats[jk]+lepids[lm], "", 40, 0.0, 4000.0, 40, 0.0, 600.0);
				
				h_Y_msoftdrop_sys[ij][jk][kl][lm][0] 		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_msoftdrop","_nom",40,0,600);
				h_X_mass_sys[ij][jk][kl][lm][0] 			 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"X_mass","_nom",40, 0.0, 4000.0);

                h_for_limit_X_mass[ij][jk][kl][lm]                 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"unrolled_bin1_X_mass","",260,0,52000);
                h_for_limit_X_mass_v2[ij][jk][kl][lm]              = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"unrolled_bin2_X_mass","",216,400,43600);

				h_for_limit_X_mass_sys[ij][jk][kl][lm][0] = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"unrolled_bin1_X_mass","_nom",260,0,52000);
                h_for_limit_X_mass_sys_v2[ij][jk][kl][lm][0] = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"unrolled_bin2_X_mass","_nom",216,400,43600);
	
				for(int isys=0; isys<nsys; isys++){
					
					char name[50];
					//up systematics
					sprintf(name,"_%s_up",sysnames[isys].Data());
					h_Y_msoftdrop_sys[ij][jk][kl][lm][2*(isys+1)-1] 		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_msoftdrop",name,40,0,600);
					h_X_mass_sys[ij][jk][kl][lm][2*(isys+1)-1]			 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"X_mass",name,40, 0.0, 4000.0);
					h_for_limit_X_mass_sys[ij][jk][kl][lm][2*(isys+1)-1] = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"unrolled_bin1_X_mass",name,260,0,52000);
                    h_for_limit_X_mass_sys_v2[ij][jk][kl][lm][2*(isys+1)-1] = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"unrolled_bin2_X_mass",name,216,400,43600);

                    //dn systematics
					sprintf(name,"_%s_dn",sysnames[isys].Data());
					h_Y_msoftdrop_sys[ij][jk][kl][lm][2*(isys+1)] 		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_msoftdrop",name,40,0,600);
					h_X_mass_sys[ij][jk][kl][lm][2*(isys+1)]			 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"X_mass",name,40, 0.0, 4000.0);
                    h_for_limit_X_mass_sys[ij][jk][kl][lm][2*(isys+1)] = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"unrolled_bin1_X_mass",name,260,0,52000);
                    h_for_limit_X_mass_sys_v2[ij][jk][kl][lm][2*(isys+1)] = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"unrolled_bin2_X_mass",name,216,400,43600);

				}
			}
	      }
       }
	}
 
   file->cd();
 
   Long64_t nn = tree->GetEntries();
   for(Long64_t j =0; j < nn ; j++)
   {
   tree->GetEntry(j);
   if( j % 10000 == 0) { std::cout <<j<<" events processed" << std::endl;}

   // Condition to avoid double counting in data //
   if((string(proc_Name[ii].Data()).find("SingleMuon")!=string::npos) && (abs(l_pdgId)!=13)) continue;
   if((string(proc_Name[ii].Data()).find("EGamma")!=string::npos) && (abs(l_pdgId)!=11)) continue;
   // end of condition

   // Tagger scale factors // (simplified version for time being)
   double b_SF, b_SF_up, b_SF_dn;
   double bb_SF, bb_SF_up, bb_SF_dn;
   double W_SF, W_SF_up, W_SF_dn;
   double Top_SF, Top_SF_up, Top_SF_dn;
   // Read tagger SFs //
   
   b_SF = b_SF_up = b_SF_dn = 1.;
   bb_SF = bb_SF_up = bb_SF_dn = 1.;
   W_SF = W_SF_up = W_SF_dn = 1.;
   Top_SF = Top_SF_up = Top_SF_dn = 1.;
   
   if(!isDATA){		 
			 
	for( int jj = 0; jj < nJetAK4 ; jj++)
	{
                       
		if(delta2R(JetAK4_eta[jj],JetAK4_phi[jj],W_eta_opt2,W_phi_opt2)>0.6 && delta2R(JetAK4_eta[jj],JetAK4_phi[jj],Y_eta,Y_phi)>0.6){
					
			int etabin = std::max(1, std::min(h_AK4M_flv_b_eff->GetNbinsY(), h_AK4M_flv_b_eff->GetYaxis()->FindBin(fabs(JetAK4_eta[jj]))));
			int ptbin  = std::max(1, std::min(h_AK4M_flv_b_eff->GetNbinsX(), h_AK4M_flv_b_eff->GetXaxis()->FindBin(JetAK4_pt[jj])));
					
			if(JetAK4_btag_DeepFlav[jj] > deep_btag_cut )
			{
				b_SF *= JetAK4_btag_DeepFlav_SF[jj];
				b_SF_up *= JetAK4_btag_DeepFlav_SF_up[jj]; 
				b_SF_dn *= JetAK4_btag_DeepFlav_SF_dn[jj];
			}
			else
			{
				if ( abs(JetAK4_hadronflav[jj]) == 5 ) {
					b_SF *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF[jj]*h_AK4M_flv_b_eff->GetBinContent(ptbin, etabin))/ ( 1.- h_AK4M_flv_b_eff->GetBinContent(ptbin, etabin)));
                    b_SF_up *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF_up[jj]*h_AK4M_flv_b_eff->GetBinContent(ptbin, etabin))/ ( 1.- h_AK4M_flv_b_eff->GetBinContent(ptbin, etabin)));
					b_SF_dn *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF_dn[jj]*h_AK4M_flv_b_eff->GetBinContent(ptbin, etabin))/ ( 1.- h_AK4M_flv_b_eff->GetBinContent(ptbin, etabin)));
                }
                else if ( abs(JetAK4_hadronflav[jj]) == 4 ) {
					b_SF *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF[jj]*h_AK4M_flv_c_eff->GetBinContent(ptbin, etabin))/ ( 1.- h_AK4M_flv_c_eff->GetBinContent(ptbin, etabin)));
                    b_SF_up *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF_up[jj]*h_AK4M_flv_c_eff->GetBinContent(ptbin, etabin))/ ( 1.- h_AK4M_flv_c_eff->GetBinContent(ptbin, etabin)));
					b_SF_dn *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF_dn[jj]*h_AK4M_flv_c_eff->GetBinContent(ptbin, etabin))/ ( 1.- h_AK4M_flv_c_eff->GetBinContent(ptbin, etabin)));
                }
				else {
                    b_SF *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF[jj]*h_AK4M_flv_l_eff->GetBinContent(ptbin, etabin))/ ( 1.- h_AK4M_flv_l_eff->GetBinContent(ptbin, etabin)));
                    b_SF_up *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF_up[jj]*h_AK4M_flv_l_eff->GetBinContent(ptbin, etabin))/ ( 1.- h_AK4M_flv_l_eff->GetBinContent(ptbin, etabin)));
                    b_SF_dn *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF_dn[jj]*h_AK4M_flv_l_eff->GetBinContent(ptbin, etabin))/ ( 1.- h_AK4M_flv_l_eff->GetBinContent(ptbin, etabin)));
                }
            }
		}
	}
                                   	    	
	int Y_pt_bin = getbinid(Y_pt,PNbb_SF_nptbins,PNbb_SF_ptbins);
	if(Y_pt_bin>=0 && Y_pt_bin<PNbb_SF_nptbins) { 
		if(Flag_Y_bb_pass_T){
			bb_SF = PNbb_SF_HP[Y_pt_bin]; 
			bb_SF_up = PNbb_SF_HP_up[Y_pt_bin]; 
			bb_SF_dn = PNbb_SF_HP_dn[Y_pt_bin]; 
		}
		/*
		else{
            int etabin = std::max(1, std::min(h_YtagT_eff->GetNbinsY(), h_YtagT_eff->GetYaxis()->FindBin(fabs(Y_eta))));
		   	int ptbin  = std::max(1, std::min(h_YtagT_eff->GetNbinsX(), h_YtagT_eff->GetXaxis()->FindBin(Y_pt)));	
		    //std::cout << etabin << "	" << ptbin << "		" << Y_DeepTag_PNetMD_XbbvsQCD << std::endl;
		    bb_SF     = std::max(1.e-6,(1.0 - PNbb_SF_HP[Y_pt_bin] * h_YtagT_eff->GetBinContent(ptbin, etabin) )/(1. - h_YtagT_eff->GetBinContent(ptbin, etabin)));    
			bb_SF_up  = std::max(1.e-6,(1.0 - PNbb_SF_HP_up[Y_pt_bin] * h_YtagT_eff->GetBinContent(ptbin, etabin) )/(1. - h_YtagT_eff->GetBinContent(ptbin, etabin)));
            bb_SF_dn  = std::max(1.e-6,(1.0 - PNbb_SF_HP_dn[Y_pt_bin] * h_YtagT_eff->GetBinContent(ptbin, etabin) )/(1. - h_YtagT_eff->GetBinContent(ptbin, etabin)));
		}
		*/ 
	}

                 
	int W_pt_bin = getbinid(W_pt_opt2,PNW_SF_nptbins,PNW_SF_ptbins);
	if(W_pt_bin>=0 && W_pt_bin<PNW_SF_nptbins) { 
		if(Flag_H_W_pass_T_opt2 && (W_msoftdrop_opt2>=65. && W_msoftdrop_opt2<=105.)){
			W_SF = PNW_SF_T[W_pt_bin]; 
			W_SF_up = PNW_SF_T_up[W_pt_bin]; 
			W_SF_dn = PNW_SF_T_dn[W_pt_bin]; 
		}
		/*
		else {
			int wetabin = std::max(1, std::min(h_WtagT_eff->GetNbinsY(), h_WtagT_eff->GetYaxis()->FindBin(fabs(W_eta_opt2))));
			int wptbin  = std::max(1, std::min(h_WtagT_eff->GetNbinsX(), h_WtagT_eff->GetXaxis()->FindBin(W_pt_opt2)));
			//std::cout << wetabin << "      " << wptbin << "         " << W_DeepTag_PNetMD_WvsQCD_opt2 << std::endl;
			W_SF     = std::max(1.e-6,(1.0 - PNW_SF_T[W_pt_bin]  * h_WtagT_eff->GetBinContent(wptbin, wetabin) )/(1. - h_WtagT_eff->GetBinContent(wptbin, wetabin)));
			W_SF_up  = std::max(1.e-6,(1.0 - PNW_SF_T_up[W_pt_bin] * h_WtagT_eff->GetBinContent(wptbin, wetabin) )/(1. - h_WtagT_eff->GetBinContent(wptbin, wetabin)));
			W_SF_dn  = std::max(1.e-6,(1.0 - PNW_SF_T_dn[W_pt_bin] * h_WtagT_eff->GetBinContent(wptbin, wetabin) )/(1. - h_WtagT_eff->GetBinContent(wptbin, wetabin)));
		}
		*/ 
	}
    //std::cout << Flag_H_W_pass_T_opt2 << "	" << W_SF << "	" << W_SF_up  << "	" << W_SF_dn << std::endl; 
   
	int Y_top_pt_bin = getbinid(Y_pt,PNTop_SF_nptbins,PNTop_SF_ptbins);
	if(Y_top_pt_bin>=0 && Y_top_pt_bin<PNTop_SF_nptbins) { 
		if(Y_DeepTag_PNet_TvsQCD>=PN_Top_med){
			Top_SF = PNTop_SF_M[Y_top_pt_bin]; 
			Top_SF_up = PNTop_SF_M_up[Y_top_pt_bin]; 
			Top_SF_dn = PNTop_SF_M_dn[Y_top_pt_bin]; 
		}
	}
	
   }//isDATA

   h_nom->Fill(0.0,1.0);
   h_nom->Fill(1.0,b_SF);
   h_nom->Fill(2.0,bb_SF);
   h_nom->Fill(3.0,W_SF);
   h_nom->Fill(4.0,bb_SF*W_SF);
   h_nom->Fill(5.0,b_SF*W_SF*bb_SF);
   h_nom->Fill(6.0,Top_SF);

   // end of tagger SF //
   
   float weight_nom;
   if(isDATA) {weight_nom = 1.0;}
   else 
   {
   if ( proc_Name[ii] == "NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_1000_MY_100_XtoYH_Nov_2021_v2.root" || proc_Name[ii] == "NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_1500_MY_200_XtoYH_Nov_2021_v2.root" || proc_Name[ii] == "NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2000_MY_200_XtoYH_Nov_2021_v2.root" || proc_Name[ii] == "NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2400_MY_300_XtoYH_Nov_2021_v2.root" || proc_Name[ii] == "NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_3000_MY_100_XtoYH_Nov_2021_v2.root" || proc_Name[ii] == "NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_3000_MY_500_XtoYH_Nov_2021_v2.root")
   {
     weight_nom = 1.0;
   }
   else
   {
     weight_nom = Generator_weight;
   }
     weight_nom *= prefiringweight;
     weight_nom *= puWeight;
     weight_nom *= leptonsf_weight;
     //weight_nom *= b_SF;           (apply b tagging scale factor later only while using b tagging based selection condition)
     weight_nom *= bb_SF;
     weight_nom *= W_SF;
   }
   
   for(int kl=0; kl<nWop; kl++)
   {
   
	vector<bool> reg_tags;
   
	bool isSR(false);
	bool isCR1(false);
	bool isCR2(false);
	bool isCR3(false);
	bool isCR4(false);
	bool isCR5(false);
	bool isCR6(false);
	bool isCR7(false);
	bool isCR8(false);
   
	bool lep_miniso = (l_minisoall<0.1)?true:false;
   
	if(kl==0)
	{   
		isSR = (Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass && lep_miniso);
		isCR1 = (!Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass && lep_miniso);
		isCR2 = (!Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass && lep_miniso);
		isCR3 = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass && lep_miniso);
		isCR4 = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && !Flag_MET_pass && lep_miniso && (Y_DeepTag_PNet_TvsQCD>=PN_Top_med));
		isCR5 = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && !Flag_MET_pass && !lep_miniso);
		isCR6 = (Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass && lep_miniso);
		isCR7 = (Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && !Flag_MET_pass && lep_miniso);
		isCR8 = (Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass && !lep_miniso);
	}
	else
	{
		isSR = (Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass && lep_miniso);
		isCR1 = (!Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass && lep_miniso);
		isCR2 = (!Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass && lep_miniso);
		isCR3 = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass && lep_miniso);
		isCR4 = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && !Flag_MET_pass && lep_miniso && (Y_DeepTag_PNet_TvsQCD>=PN_Top_med));
		isCR5 = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && !Flag_MET_pass && !lep_miniso);
		isCR6 = (Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass && lep_miniso); 
		isCR7 = (Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && !Flag_MET_pass && lep_miniso);
		isCR8 = (Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass && !lep_miniso);  
	}
   
	reg_tags.push_back(isSR);
	reg_tags.push_back(isCR1);
	reg_tags.push_back(isCR2);
	reg_tags.push_back(isCR3);
	reg_tags.push_back(isCR4);
	reg_tags.push_back(isCR5);
	reg_tags.push_back(isCR6);
	reg_tags.push_back(isCR7);
	reg_tags.push_back(isCR8);
	
	int ireg = get_region(reg_tags);
   
	if(ireg<0||ireg>=nrgn) continue;
	
	// conditions to purify control regions //
	
	bool pure_pass = true;
	
	if(isCR1){	pure_pass = (Y_DeepTag_PNetMD_XbbvsQCD>=0.4);	}
	if(isCR2){	pure_pass = (Y_DeepTag_PNetMD_XbbvsQCD>=0.1);	}
	if(isCR3){										
		if(kl==0) { pure_pass = (W_msoftdrop_opt1<=60.); }
		if(kl==1) { pure_pass = (W_msoftdrop_opt2<=60.); }
	}
	if(isCR4){	pure_pass = (Y_DeepTag_PNetMD_XbbvsQCD>=0.1);	}
	if(isCR6){	pure_pass = (MET_pt>=75.);	}
	
	// end of purification conditions //		
	
	if(kl==1){
		h_nom_reg[ireg]->Fill(0.0,1.0);
		h_nom_reg[ireg]->Fill(1.0,b_SF);
		h_nom_reg[ireg]->Fill(2.0,bb_SF);
		h_nom_reg[ireg]->Fill(3.0,W_SF);
		h_nom_reg[ireg]->Fill(4.0,bb_SF*W_SF);
		h_nom_reg[ireg]->Fill(5.0,b_SF*W_SF*bb_SF);
		h_nom_reg[ireg]->Fill(6.0,Top_SF);
	}
   
	int jk_b = (nbjets_other==0)?1:2;
   
	float MT = sqrt(2*l_pt*MET_pt*(1-cos(PhiInRange(l_phi-MET_phi))));
   
	bool lepid_info[nlid];
	lepid_info[0] = true;
	if(abs(l_pdgId)==13) { lepid_info[1] = true;  }
	if(abs(l_pdgId)==11) { lepid_info[2] = true;  }
   
	// first fill few general histograms //
   
	h_reg->Fill(ireg,weight_nom);
   
	// now fill histograms binned in analysis categories //
   
    float weight = weight_nom;
    
    if(isCR4) { weight = weight_nom*Top_SF; } // since top tagging condition is used only in CR4
   
	for(int jk=0; jk<nbcat; jk++){
	   
		if(!(jk==0 || jk==jk_b)) continue;
		
		if(jk!=0) { weight = weight_nom*b_SF; } // applying b tagging SF only if any condition on number of b-tagged jets is used
		
		for(int lm=0; lm<nlid; lm++){
			
			if(!lepid_info[lm]) continue;
			
			// reference histograms corresponding to purity conditions //
			
			h_Y_PNetMD_XbbvsQCD_ref[ireg][jk][kl][lm]->Fill(Y_DeepTag_PNetMD_XbbvsQCD,weight);
			h_MET_pt_ref[ireg][jk][kl][lm]->Fill(MET_pt,weight); 
			if(kl==0){  
				h_W_msoftdrop_ref[ireg][jk][kl][lm]->Fill(W_msoftdrop_opt1,weight); 
			}
			if(kl==1){  
				h_W_msoftdrop_ref[ireg][jk][kl][lm]->Fill(W_msoftdrop_opt2,weight); 
			}
			
			// ref end //
			
			// Now apply the purity condition & fill all histograms (including sys)
			
			if(pure_pass)
			{
   
				h_l_pt[ireg][jk][kl][lm]->Fill(l_pt,weight); 
				h_l_eta[ireg][jk][kl][lm]->Fill(l_eta,weight); 
				h_l_minisoall[ireg][jk][kl][lm]->Fill(l_minisoall,weight); 
					
				h_MET_pt[ireg][jk][kl][lm]->Fill(MET_pt,weight); 
				h_MET_sig[ireg][jk][kl][lm]->Fill(MET_sig,weight); 
				h_MT[ireg][jk][kl][lm]->Fill(MT,weight); 
	        
				h_Y_pt[ireg][jk][kl][lm]->Fill(Y_pt,weight); 
				h_Y_msoftdrop[ireg][jk][kl][lm]->Fill(Y_msoftdrop,weight);
				h_Y_PNetMD_XbbvsQCD[ireg][jk][kl][lm]->Fill(Y_DeepTag_PNetMD_XbbvsQCD,weight);
				h_Y_PNetMD_WvsQCD[ireg][jk][kl][lm]->Fill(Y_DeepTag_PNetMD_WvsQCD,weight);
				h_Y_PNet_TvsQCD[ireg][jk][kl][lm]->Fill(Y_DeepTag_PNet_TvsQCD,weight);
				h_Y_sub1_mass[ireg][jk][kl][lm]->Fill(Y_sub1_mass,weight);
				h_Y_sub2_mass[ireg][jk][kl][lm]->Fill(Y_sub2_mass,weight);
				h_Y_sub1_btag[ireg][jk][kl][lm]->Fill(Y_sub1_btag,weight);
				h_Y_sub2_btag[ireg][jk][kl][lm]->Fill(Y_sub2_btag,weight);
                        
				float X_conv_mass;

				if(kl==0){           
			
					h_W_pt[ireg][jk][kl][lm]->Fill(W_pt_opt1,weight); 
					h_W_msoftdrop[ireg][jk][kl][lm]->Fill(W_msoftdrop_opt1,weight); 
					h_W_PNetMD_XbbvsQCD[ireg][jk][kl][lm]->Fill(W_DeepTag_PNetMD_XbbvsQCD_opt1,weight); 
					h_W_PNetMD_WvsQCD[ireg][jk][kl][lm]->Fill(W_DeepTag_PNetMD_WvsQCD_opt1,weight); 
					h_W_PNet_TvsQCD[ireg][jk][kl][lm]->Fill(W_DeepTag_PNet_TvsQCD_opt1,weight); 
					h_W_DAK8MD_WvsQCD[ireg][jk][kl][lm]->Fill(W_DeepTag_DAK8MD_WvsQCD_opt1,weight); 
					h_W_sub1_mass[ireg][jk][kl][lm]->Fill(W_sub1_mass_opt1,weight); 
					h_W_sub2_mass[ireg][jk][kl][lm]->Fill(W_sub2_mass_opt1,weight); 
					h_W_sub1_btag[ireg][jk][kl][lm]->Fill(W_sub1_btag_opt1,weight); 
					h_W_sub2_btag[ireg][jk][kl][lm]->Fill(W_sub2_btag_opt1,weight); 
			
					h_dR_lW[ireg][jk][kl][lm]->Fill(dR_lW_opt1,weight); 
					h_dy_lW[ireg][jk][kl][lm]->Fill(dy_lW_opt1,weight); 
					h_dphi_lW[ireg][jk][kl][lm]->Fill(dphi_lW_opt1,weight);
			
					h_H_mass[ireg][jk][kl][lm]->Fill(H_mass_opt1,weight); 
					h_X_mass[ireg][jk][kl][lm]->Fill(X_mass_opt1,weight); 
					
					h_X_Y_mass[ireg][jk][kl][lm]->Fill(X_mass_opt1, Y_msoftdrop, weight);
			        
					X_conv_mass = X_mass_opt1 + 4000.0 * get_Y_id(Y_msoftdrop);
					h_for_limit_X_mass[ireg][jk][kl][lm]->Fill(X_conv_mass,weight);
					if(X_mass_opt1>400) {h_for_limit_X_mass_v2[ireg][jk][kl][lm]->Fill(X_mass_opt1 + (4000.0 - 400.0) * get_Y_id(Y_msoftdrop)   ,weight);}

				}
				else
				{
			
					h_W_pt[ireg][jk][kl][lm]->Fill(W_pt_opt2,weight); 
					h_W_msoftdrop[ireg][jk][kl][lm]->Fill(W_msoftdrop_opt2,weight); 
					h_W_PNetMD_XbbvsQCD[ireg][jk][kl][lm]->Fill(W_DeepTag_PNetMD_XbbvsQCD_opt2,weight); 
					h_W_PNetMD_WvsQCD[ireg][jk][kl][lm]->Fill(W_DeepTag_PNetMD_WvsQCD_opt2,weight); 
					h_W_PNet_TvsQCD[ireg][jk][kl][lm]->Fill(W_DeepTag_PNet_TvsQCD_opt2,weight); 
					h_W_DAK8MD_WvsQCD[ireg][jk][kl][lm]->Fill(W_DeepTag_DAK8MD_WvsQCD_opt2,weight); 
					h_W_sub1_mass[ireg][jk][kl][lm]->Fill(W_sub1_mass_opt2,weight); 
					h_W_sub2_mass[ireg][jk][kl][lm]->Fill(W_sub2_mass_opt2,weight); 
					h_W_sub1_btag[ireg][jk][kl][lm]->Fill(W_sub1_btag_opt2,weight); 
					h_W_sub2_btag[ireg][jk][kl][lm]->Fill(W_sub2_btag_opt2,weight); 
			
					h_dR_lW[ireg][jk][kl][lm]->Fill(dR_lW_opt2,weight); 
					h_dy_lW[ireg][jk][kl][lm]->Fill(dy_lW_opt2,weight); 
					h_dphi_lW[ireg][jk][kl][lm]->Fill(dphi_lW_opt2,weight);
			
					h_H_mass[ireg][jk][kl][lm]->Fill(H_mass_opt2,weight); 
					h_X_mass[ireg][jk][kl][lm]->Fill(X_mass_opt2,weight); 
					
					h_X_Y_mass[ireg][jk][kl][lm]->Fill(X_mass_opt2, Y_msoftdrop, weight);
					X_conv_mass = X_mass_opt2 + 4000.0 * get_Y_id(Y_msoftdrop);
					h_for_limit_X_mass[ireg][jk][kl][lm]->Fill(X_conv_mass,weight);
					if(X_mass_opt2>400){h_for_limit_X_mass_v2[ireg][jk][kl][lm]->Fill(X_mass_opt2 + (4000.0 - 400.0) * get_Y_id(Y_msoftdrop)   ,weight);}

				}
		
				h_nbjets_other[ireg][jk][kl][lm]->Fill(nbjets_other,weight);
		
				float X_mass;
				if(kl==0) { X_mass = X_mass_opt1; }
				else { X_mass = X_mass_opt2; }
			
				h_X_mass_sys[ireg][jk][kl][lm][0]->Fill(X_mass,weight); 
				h_Y_msoftdrop_sys[ireg][jk][kl][lm][0]->Fill(Y_msoftdrop,weight);
				h_for_limit_X_mass_sys[ireg][jk][kl][lm][0]->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight);
				if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][0]->Fill(X_mass + (4000.0 -400.0) * get_Y_id(Y_msoftdrop),weight);}
		
				float X_JESup, X_JESdn, X_JERup, X_JERdn;
		
				TLorentzVector Y_cand, H_cand;
		
				if(kl==0){
					H_cand.SetPtEtaPhiM(H_pt_opt1*H_JESup_opt1,H_eta_opt1,H_phi_opt1,H_mass_opt1*H_JESup_opt1);
					Y_cand.SetPtEtaPhiM(Y_pt*Y_JESup,Y_eta,Y_phi,Y_mass*Y_JESup);
					X_JESup = (Y_cand+H_cand).M()/X_mass;
					H_cand.SetPtEtaPhiM(H_pt_opt1*H_JESdn_opt1,H_eta_opt1,H_phi_opt1,H_mass_opt1*H_JESdn_opt1);
					Y_cand.SetPtEtaPhiM(Y_pt*Y_JESdn,Y_eta,Y_phi,Y_mass*Y_JESdn);
					X_JESdn = (Y_cand+H_cand).M()/X_mass;
					H_cand.SetPtEtaPhiM(H_pt_opt1*H_JERup_opt1,H_eta_opt1,H_phi_opt1,H_mass_opt1);
					Y_cand.SetPtEtaPhiM(Y_pt*Y_JERup,Y_eta,Y_phi,Y_mass);
					X_JERup = (Y_cand+H_cand).M()/X_mass;
					H_cand.SetPtEtaPhiM(H_pt_opt1*H_JERdn_opt1,H_eta_opt1,H_phi_opt1,H_mass_opt1);
					Y_cand.SetPtEtaPhiM(Y_pt*Y_JERdn,Y_eta,Y_phi,Y_mass);
					X_JERdn = (Y_cand+H_cand).M()/X_mass;
				}
				else
				{
					H_cand.SetPtEtaPhiM(H_pt_opt2*H_JESup_opt2,H_eta_opt2,H_phi_opt2,H_mass_opt2*H_JESup_opt2);
					Y_cand.SetPtEtaPhiM(Y_pt*Y_JESup,Y_eta,Y_phi,Y_mass*Y_JESup);
					X_JESup = (Y_cand+H_cand).M()/X_mass;
					H_cand.SetPtEtaPhiM(H_pt_opt2*H_JESdn_opt2,H_eta_opt2,H_phi_opt2,H_mass_opt2*H_JESdn_opt2);
					Y_cand.SetPtEtaPhiM(Y_pt*Y_JESdn,Y_eta,Y_phi,Y_mass*Y_JESdn);
					X_JESdn = (Y_cand+H_cand).M()/X_mass;
					H_cand.SetPtEtaPhiM(H_pt_opt2*H_JERup_opt2,H_eta_opt2,H_phi_opt2,H_mass_opt2);
					Y_cand.SetPtEtaPhiM(Y_pt*Y_JERup,Y_eta,Y_phi,Y_mass);
					X_JERup = (Y_cand+H_cand).M()/X_mass;
					H_cand.SetPtEtaPhiM(H_pt_opt2*H_JERdn_opt2,H_eta_opt2,H_phi_opt2,H_mass_opt2);
					Y_cand.SetPtEtaPhiM(Y_pt*Y_JERdn,Y_eta,Y_phi,Y_mass);
					X_JERdn = (Y_cand+H_cand).M()/X_mass;
				}

				for(int isys=0; isys<nsys; isys++){
				
					if(isys==0) // JES
					{
						h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass*X_JESup,weight); 
						h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]	->Fill(X_mass*X_JESdn,weight); 
				
						h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(Y_msoftdrop*Y_JESup,weight);
						h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)]	 ->Fill(Y_msoftdrop*Y_JESdn,weight);

						h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass*X_JESup + 4000.0 * get_Y_id(Y_msoftdrop*Y_JESup),weight);
						h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass*X_JESdn + 4000.0 * get_Y_id(Y_msoftdrop*Y_JESdn),weight);
						if(X_mass*X_JESup>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass*X_JESup + (4000.0 -400.0) * get_Y_id(Y_msoftdrop*Y_JESup),weight);}
						if(X_mass*X_JESdn>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)]	->Fill(X_mass*X_JESdn + (4000.0 -400.0) * get_Y_id(Y_msoftdrop*Y_JESdn),weight);}
					}
					if(isys==1) // JER
					{
						h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass*X_JERup,weight); 
						h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]	->Fill(X_mass*X_JERdn,weight); 
				
						h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(Y_msoftdrop,weight);
						h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)]	 ->Fill(Y_msoftdrop,weight);

						h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass*X_JERup + 4000.0 * get_Y_id(Y_msoftdrop),weight);
						h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass*X_JERdn + 4000.0 * get_Y_id(Y_msoftdrop),weight);
						if(X_mass*X_JERup>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass*X_JERup + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight);}
						if(X_mass*X_JERdn>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass*X_JERdn + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight);}
					}
					if(isys==2)  // PU reweighting 
					{
						h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass,weight*puWeightup/TMath::Max(float(1.e-6),puWeight)); 
						h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]	->Fill(X_mass,weight*puWeightdown/TMath::Max(float(1.e-6),puWeight)); 
				
						h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(Y_msoftdrop,weight*puWeightup/TMath::Max(float(1.e-6),puWeight));
						h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)]	 ->Fill(Y_msoftdrop,weight*puWeightdown/TMath::Max(float(1.e-6),puWeight));

						h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight*puWeightup/TMath::Max(float(1.e-6),puWeight));
						h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight*puWeightdown/TMath::Max(float(1.e-6),puWeight));
						if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight*puWeightup/TMath::Max(float(1.e-6),puWeight));}
						if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight*puWeightdown/TMath::Max(float(1.e-6),puWeight));}
					}
					if(isys==3) // LeptonSF stat & syst
					{
						h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass,weight*leptonsf_weight_stat/TMath::Max(float(1.e-6),leptonsf_weight)); 
						h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]	->Fill(X_mass,weight*leptonsf_weight_syst/TMath::Max(float(1.e-6),leptonsf_weight)); 
				
						h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(Y_msoftdrop,weight*leptonsf_weight_stat/TMath::Max(float(1.e-6),leptonsf_weight));
						h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)]	 ->Fill(Y_msoftdrop,weight*leptonsf_weight_syst/TMath::Max(float(1.e-6),leptonsf_weight));

						h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight*leptonsf_weight_stat/TMath::Max(float(1.e-6),leptonsf_weight));
						h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight*leptonsf_weight_syst/TMath::Max(float(1.e-6),leptonsf_weight));
						if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight*leptonsf_weight_stat/TMath::Max(float(1.e-6),leptonsf_weight));}
						if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight*leptonsf_weight_syst/TMath::Max(float(1.e-6),leptonsf_weight));}
					}
					if(isys==4) // Lepton SF statistical component only
					{
						h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass,weight*leptonsf_weight_up/TMath::Max(float(1.e-6),leptonsf_weight));
						h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]	->Fill(X_mass,weight*leptonsf_weight_dn/TMath::Max(float(1.e-6),leptonsf_weight));

						h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(Y_msoftdrop,weight*leptonsf_weight_up/TMath::Max(float(1.e-6),leptonsf_weight));
						h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)]	 ->Fill(Y_msoftdrop,weight*leptonsf_weight_dn/TMath::Max(float(1.e-6),leptonsf_weight));

						h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight*leptonsf_weight_up/TMath::Max(float(1.e-6),leptonsf_weight));
						h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight*leptonsf_weight_dn/TMath::Max(float(1.e-6),leptonsf_weight));
						if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight*leptonsf_weight_up/TMath::Max(float(1.e-6),leptonsf_weight));}
						if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight*leptonsf_weight_dn/TMath::Max(float(1.e-6),leptonsf_weight));}
					}
					if(isys==5) // Prefiring weight 
					{
						h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass,weight*prefiringweightup/TMath::Max(double(1.e-6),prefiringweight)); 
						h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]	->Fill(X_mass,weight*prefiringweightdown/TMath::Max(double(1.e-6),prefiringweight)); 
				
						h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(Y_msoftdrop,weight*prefiringweightup/TMath::Max(double(1.e-6),prefiringweight));
						h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)]	 ->Fill(Y_msoftdrop,weight*prefiringweightdown/TMath::Max(double(1.e-6),prefiringweight));

						h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight*prefiringweightup/TMath::Max(double(1.e-6),prefiringweight));
						h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight*prefiringweightdown/TMath::Max(double(1.e-6),prefiringweight));
						if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight*prefiringweightup/TMath::Max(double(1.e-6),prefiringweight));}
						if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight*prefiringweightdown/TMath::Max(double(1.e-6),prefiringweight));}
					}
					if(isys==6) // Xbb tagging scale factors
					{
						h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass,weight*bb_SF_up/TMath::Max(double(1.e-6),bb_SF)); 
						h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]	->Fill(X_mass,weight*bb_SF_dn/TMath::Max(double(1.e-6),bb_SF)); 
				
						h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(Y_msoftdrop,weight*bb_SF_up/TMath::Max(double(1.e-6),bb_SF));
						h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)]	 ->Fill(Y_msoftdrop,weight*bb_SF_dn/TMath::Max(double(1.e-6),bb_SF));

						h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight*bb_SF_up/TMath::Max(double(1.e-6),bb_SF));
						h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight*bb_SF_dn/TMath::Max(double(1.e-6),bb_SF));
						if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight*bb_SF_up/TMath::Max(double(1.e-6),bb_SF));}
						if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight*bb_SF_dn/TMath::Max(double(1.e-6),bb_SF));}
					}
					if(isys==7) // W tagging scale factors
					{
						h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass,weight*W_SF_up/TMath::Max(double(1.e-6),W_SF)); 
						h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]	->Fill(X_mass,weight*W_SF_dn/TMath::Max(double(1.e-6),W_SF)); 
				
						h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(Y_msoftdrop,weight*W_SF_up/TMath::Max(double(1.e-6),W_SF));
						h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)]	 ->Fill(Y_msoftdrop,weight*W_SF_dn/TMath::Max(double(1.e-6),W_SF));

						h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight*W_SF_up/TMath::Max(double(1.e-6),W_SF));
						h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight*W_SF_dn/TMath::Max(double(1.e-6),W_SF));
						if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight*W_SF_up/TMath::Max(double(1.e-6),W_SF));}
						if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight*W_SF_dn/TMath::Max(double(1.e-6),W_SF));}
					}
					if(isys==8) // B tagging scale factors
					{
						if(jk!=0){
							h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass,weight*b_SF_up/TMath::Max(double(1.e-6),b_SF));
							h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]	->Fill(X_mass,weight*b_SF_dn/TMath::Max(double(1.e-6),b_SF));

							h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(Y_msoftdrop,weight*b_SF_up/TMath::Max(double(1.e-6),b_SF));
							h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)]	 ->Fill(Y_msoftdrop,weight*b_SF_dn/TMath::Max(double(1.e-6),b_SF));

							h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight*b_SF_up/TMath::Max(double(1.e-6),b_SF));
							h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight*b_SF_dn/TMath::Max(double(1.e-6),b_SF));
							if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight*b_SF_up/TMath::Max(double(1.e-6),b_SF));}
							if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight*b_SF_dn/TMath::Max(double(1.e-6),b_SF));}
						}
						else{
							h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass,weight);
							h_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]	->Fill(X_mass,weight);

							h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(Y_msoftdrop,weight);
							h_Y_msoftdrop_sys[ireg][jk][kl][lm][2*(isys+1)]	 ->Fill(Y_msoftdrop,weight);

							h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight);
							h_for_limit_X_mass_sys[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight);
							if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)-1]->Fill(X_mass + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight);}
							if(X_mass>400){h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][2*(isys+1)]  ->Fill(X_mass + (4000.0-400.0) * get_Y_id(Y_msoftdrop),weight);}
						}
					}	
                                
				
				}//sys loop
		
			}//purity condition
			
		}// lepid (lm)
		
	  }// b cat (jk)
	}// W opt (kl)
	
   }// end of event loop
   
    final_file->Write();
    final_file->cd();
   
    final_file->Close();
    FEff->Close();
  }
}

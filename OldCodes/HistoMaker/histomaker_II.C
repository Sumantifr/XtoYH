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
using namespace std;
double PhiInRange(const double& phi) {
  double phiout = phi;

  if( phiout > 2*M_PI || phiout < -2*M_PI) {
    phiout = fmod( phiout, 2*M_PI);
  }
  if (phiout <= -M_PI) phiout += 2*M_PI;
  else if (phiout >  M_PI) phiout -= 2*M_PI;

  return phiout;
}

void histomaker_II()
{

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
"QCD_HT300to500_XtoYH_Nov_2021.root",
"QCD_HT500to700_XtoYH_Nov_2021.root",
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

   std::cout << proc_Name[ii] << std::endl;
   TFile* final_file = TFile::Open("Histograms_II_l_iso/Histogram_"+proc_Name[ii], "RECREATE");  

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
   Int_t           nBJetAK4;
   Float_t         BJetAK4_pt[narray];   //[_s_nBJetAK4]
   Float_t         BJetAK4_eta[narray];   //[_s_nBJetAK4]
   Float_t         BJetAK4_phi[narray];   //[_s_nBJetAK4]
   Float_t         BJetAK4_mass[narray];   //[_s_nBJetAK4]
   Float_t         BJetAK4_btag_DeepCSV[narray];   //[_s_nBJetAK4]
   Float_t         BJetAK4_btag_DeepFlav[narray];   //[_s_nBJetAK4]
   Int_t           BJetAK4_hadronflav[narray];   //[_s_nBJetAK4]
   Int_t           BJetAK4_partonflav[narray];   //[_s_nBJetAK4]
   Float_t         BJetAK4_qgl[narray];   //[_s_nBJetAK4]
   Float_t         BJetAK4_PUID[narray];   //[_s_nBJetAK4]
   Float_t         BJetAK4_JESup[narray];   //[_s_nBJetAK4]
   Float_t         BJetAK4_JESdn[narray];   //[_s_nBJetAK4]
   Float_t         BJetAK4_JERup[narray];   //[_s_nBJetAK4]
   Float_t         BJetAK4_JERdn[narray];   //[_s_nBJetAK4]
   Float_t         BJetAK4_btag_DeepFlav_SF[narray];   //[_s_nBJetAK4]
   Float_t         BJetAK4_btag_DeepFlav_SF_up[narray];   //[_s_nBJetAK4]
   Float_t         BJetAK4_btag_DeepFlav_SF_dn[narray];   //[_s_nBJetAK4]
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
   TBranch        *b__s_nBJetAK4;   //!
   TBranch        *b_BJetAK4_pt;   //!
   TBranch        *b_BJetAK4_eta;   //!
   TBranch        *b_BJetAK4_phi;   //!
   TBranch        *b_BJetAK4_mass;   //!
   TBranch        *b_BJetAK4_btag_DeepCSV;   //!
   TBranch        *b_BJetAK4_btag_DeepFlav;   //!
   TBranch        *b_BJetAK4_hadronflav;   //!
   TBranch        *b_BJetAK4_partonflav;   //!
   TBranch        *b_BJetAK4_qgl;   //!
   TBranch        *b_BJetAK4_PUID;   //!
   TBranch        *b_BJetAK4_JESup;   //!
   TBranch        *b_BJetAK4_JESdn;   //!
   TBranch        *b_BJetAK4_JERup;   //!
   TBranch        *b_BJetAK4_JERdn;   //!
   TBranch        *b_BJetAK4_btag_DeepFlav_SF;   //!
   TBranch        *b_BJetAK4_btag_DeepFlav_SF_up;   //!
   TBranch        *b_BJetAK4_btag_DeepFlav_SF_dn;   //!
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
   tree->SetBranchAddress("nBJetAK4", &nBJetAK4, &b__s_nBJetAK4);
   tree->SetBranchAddress("BJetAK4_pt", BJetAK4_pt, &b_BJetAK4_pt);
   tree->SetBranchAddress("BJetAK4_eta", BJetAK4_eta, &b_BJetAK4_eta);
   tree->SetBranchAddress("BJetAK4_phi", BJetAK4_phi, &b_BJetAK4_phi);
   tree->SetBranchAddress("BJetAK4_mass", BJetAK4_mass, &b_BJetAK4_mass);
   tree->SetBranchAddress("BJetAK4_btag_DeepCSV", BJetAK4_btag_DeepCSV, &b_BJetAK4_btag_DeepCSV);
   tree->SetBranchAddress("BJetAK4_btag_DeepFlav", BJetAK4_btag_DeepFlav, &b_BJetAK4_btag_DeepFlav);
   tree->SetBranchAddress("BJetAK4_hadronflav", BJetAK4_hadronflav, &b_BJetAK4_hadronflav);
   tree->SetBranchAddress("BJetAK4_partonflav", BJetAK4_partonflav, &b_BJetAK4_partonflav);
   tree->SetBranchAddress("BJetAK4_qgl", BJetAK4_qgl, &b_BJetAK4_qgl);
   tree->SetBranchAddress("BJetAK4_PUID", BJetAK4_PUID, &b_BJetAK4_PUID);
   tree->SetBranchAddress("BJetAK4_JESup", BJetAK4_JESup, &b_BJetAK4_JESup);
   tree->SetBranchAddress("BJetAK4_JESdn", BJetAK4_JESdn, &b_BJetAK4_JESdn);
   tree->SetBranchAddress("BJetAK4_JERup", BJetAK4_JERup, &b_BJetAK4_JERup);
   tree->SetBranchAddress("BJetAK4_JERdn", BJetAK4_JERdn, &b_BJetAK4_JERdn);
   tree->SetBranchAddress("BJetAK4_btag_DeepFlav_SF", BJetAK4_btag_DeepFlav_SF, &b_BJetAK4_btag_DeepFlav_SF);
   tree->SetBranchAddress("BJetAK4_btag_DeepFlav_SF_up", BJetAK4_btag_DeepFlav_SF_up, &b_BJetAK4_btag_DeepFlav_SF_up);
   tree->SetBranchAddress("BJetAK4_btag_DeepFlav_SF_dn", BJetAK4_btag_DeepFlav_SF_dn, &b_BJetAK4_btag_DeepFlav_SF_dn);
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

 
 
   // tree->SetBranchAddress("event_weight_LHE", &event_weight_LHE);

  double XEDGES[] = {20, 25, 30, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1101, 2500};
  const int XBINS = sizeof(XEDGES)/sizeof(XEDGES[0])-1;

  TString rgn[] = {"SR","CR1","CR2","CR3","CR4","CR5","CR6","CR7"};
  int nrgn = sizeof(rgn)/sizeof(rgn[0]);

  TString Ytype[] = {"T", "M" , "L"};
  int nYtype = sizeof(Ytype)/sizeof(Ytype[0]);

  TString Wtype[] = {"L", "M" , "T"};
  int nWtype = sizeof(Wtype)/sizeof(Wtype[0]);
  
  TH1F* h_opt1_l_pt[nrgn*9*2];
  TH1F* h_opt1_l_eta[nrgn*9*2];
  TH1F* h_opt1_MET_pt[nrgn*9*2];
  TH1F* h_opt1_MET_sig[nrgn*9*2]; 
  TH1F* h_opt1_MT[nrgn*9*2];  
  TH1F* h_opt1_leadbjet_pt[nrgn*9*2];  
  TH1F* h_opt1_leadbjet_btag_DeepFlav[nrgn*9*2];  
  
  TH1F *h_opt1_Y_pt[nrgn*9*2];
  TH1F *h_opt1_Y_msoftdrop[nrgn*9*2];
  TH1F *h_opt1_Y_PNetMD_XbbvsQCD[nrgn*9*2];
  TH1F *h_opt1_Y_PNetMD_WvsQCD[nrgn*9*2];
  TH1F *h_opt1_Y_PNet_TvsQCD[nrgn*9*2];
  TH1F *h_opt1_Y_sub1_mass[nrgn*9*2];
  TH1F *h_opt1_Y_sub2_mass[nrgn*9*2];
  TH1F *h_opt1_Y_sub1_btag[nrgn*9*2];
  TH1F *h_opt1_Y_sub2_btag[nrgn*9*2];
		   
  TH1F* h_opt1_W_pt_opt1[nrgn*9*2];
  TH1F* h_opt1_W_msoftdrop_opt1[nrgn*9*2]; 
  TH1F* h_opt1_W_PNetMD_XbbvsQCD_opt1[nrgn*9*2];
  TH1F* h_opt1_W_PNetMD_WvsQCD_opt1[nrgn*9*2];
  TH1F* h_opt1_W_PNet_TvsQCD_opt1[nrgn*9*2];
  TH1F* h_opt1_W_DAK8MD_WvsQCD_opt1[nrgn*9*2];
  TH1F* h_opt1_W_sub1_mass_opt1[nrgn*9*2];
  TH1F* h_opt1_W_sub2_mass_opt1[nrgn*9*2];
  TH1F* h_opt1_W_sub1_btag_opt1[nrgn*9*2];
  TH1F* h_opt1_W_sub2_btag_opt1[nrgn*9*2];
  TH1F* h_opt1_dR_lW_opt1[nrgn*9*2];
  TH1F* h_opt1_dy_lW_opt1[nrgn*9*2];
  TH1F* h_opt1_dphi_lW_opt1[nrgn*9*2];
				
  TH1F* h_opt1_H_mass_opt1[nrgn*9*2]; 
  TH1F* h_opt1_X_mass_opt1[nrgn*9*2]; 
  TH1F* h_opt1_nbjets_other[nrgn*9*2]; 
  TH2F* h_opt1_X_Y_mass[nrgn*9*2];   
  
  TH1F* h_opt2_l_pt[nrgn*9*2];
  TH1F* h_opt2_l_eta[nrgn*9*2];
  TH1F* h_opt2_MET_pt[nrgn*9*2];
  TH1F* h_opt2_MET_sig[nrgn*9*2]; 
  TH1F* h_opt2_MT[nrgn*9*2];  
  TH1F *h_opt2_leadbjet_pt[nrgn*9*2];  
  TH1F *h_opt2_leadbjet_btag_DeepFlav[nrgn*9*2];  
  
  TH1F *h_opt2_Y_pt[nrgn*9*2];
  TH1F *h_opt2_Y_msoftdrop[nrgn*9*2];
  TH1F *h_opt2_Y_PNetMD_XbbvsQCD[nrgn*9*2];
  TH1F *h_opt2_Y_PNetMD_WvsQCD[nrgn*9*2];
  TH1F *h_opt2_Y_PNet_TvsQCD[nrgn*9*2];
  TH1F *h_opt2_Y_sub1_mass[nrgn*9*2];
  TH1F *h_opt2_Y_sub2_mass[nrgn*9*2];
  TH1F *h_opt2_Y_sub1_btag[nrgn*9*2];
  TH1F *h_opt2_Y_sub2_btag[nrgn*9*2];
  
  TH1F *h_opt2_W_pt_opt2[nrgn*9*2];
  TH1F* h_opt2_W_msoftdrop_opt2[nrgn*9*2];
  TH1F* h_opt2_W_PNetMD_XbbvsQCD_opt2[nrgn*9*2];
  TH1F* h_opt2_W_PNetMD_WvsQCD_opt2[nrgn*9*2];
  TH1F* h_opt2_W_PNet_TvsQCD_opt2[nrgn*9*2];
  TH1F* h_opt2_W_DAK8MD_WvsQCD_opt2[nrgn*9*2];
  TH1F* h_opt2_W_sub1_mass_opt2[nrgn*9*2];
  TH1F* h_opt2_W_sub2_mass_opt2[nrgn*9*2];
  TH1F* h_opt2_W_sub1_btag_opt2[nrgn*9*2];
  TH1F* h_opt2_W_sub2_btag_opt2[nrgn*9*2];
  TH1F* h_opt2_dR_lW_opt2[nrgn*9*2];
  TH1F* h_opt2_dy_lW_opt2[nrgn*9*2];
  TH1F* h_opt2_dphi_lW_opt2[nrgn*9*2];
  
  TH1F* h_opt2_H_mass_opt2[nrgn*9*2];
  TH1F* h_opt2_X_mass_opt2[nrgn*9*2];
  TH1F* h_opt2_nbjets_other[nrgn*9*2];
  TH2F* h_opt2_X_Y_mass[nrgn*9*2];
 

  int pc = 0;
  for (int ij=0 ; ij< nrgn ; ij++)
  {
	  
      for (int yt=0 ; yt< 3 ; yt++)
       {
            for (int wt=0 ; wt< 3 ; wt++)
                {

                                //*******************************************************************************//
                                //*******************************************************************************//
				h_opt1_l_pt[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_l_pt_"+rgn[ij], "",40,0,1000);
				h_opt1_l_eta[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_l_eta_"+rgn[ij], "",50,-2.5,2.5);
	       	                h_opt1_MET_pt[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_MET_pt_"+rgn[ij], "",40,0,1000);
				h_opt1_MET_sig[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_MET_sig_"+rgn[ij], "",20,0,300);
				h_opt1_MT[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_MT_"+rgn[ij], "",20,0,300);
				h_opt1_leadbjet_pt[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_LeadBJet_pt_"+rgn[ij], "",40,0,1000);
				h_opt1_leadbjet_btag_DeepFlav[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_LeadBJet_btag_DeepFlav_"+rgn[ij], "",50,0,1);
				h_opt1_Y_pt[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_Y_pt_"+rgn[ij], "",XBINS, XEDGES );
				h_opt1_Y_msoftdrop[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_Y_msoftdrop_"+rgn[ij], "", 20, 0.0, 500.0 );
				h_opt1_Y_PNetMD_XbbvsQCD[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_Y_PNetMD_XbbvsQCD_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt1_Y_PNetMD_WvsQCD[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_Y_PNetMD_WvsQCD_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt1_Y_PNet_TvsQCD[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_Y_PNet_TvsQCD_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt1_Y_sub1_mass[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_Y_sub1_mass_"+rgn[ij], "", 20, 0.0, 300 );
				h_opt1_Y_sub2_mass[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_Y_sub2_mass_"+rgn[ij], "", 20, 0.0, 300 );
				h_opt1_Y_sub1_btag[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_Y_sub1_btag_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt1_Y_sub2_btag[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_Y_sub2_btag_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt1_W_pt_opt1[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_W_pt_opt1_"+rgn[ij], "",XBINS, XEDGES );
				h_opt1_W_msoftdrop_opt1[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_W_msoftdrop_opt1_"+rgn[ij], "", 35, 0.0, 350.0 );
				h_opt1_W_PNetMD_XbbvsQCD_opt1[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_W_PNetMD_XbbvsQCD_opt1_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt1_W_PNetMD_WvsQCD_opt1[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_W_PNetMD_WvsQCD_opt1_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt1_W_PNet_TvsQCD_opt1[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_W_PNet_TvsQCD_opt1_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt1_W_DAK8MD_WvsQCD_opt1[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_W_DAK8MD_WvsQCD_opt1_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt1_W_sub1_mass_opt1[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_W_sub1_mass_opt1_"+rgn[ij], "", 20, 0.0, 300 );
				h_opt1_W_sub2_mass_opt1[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_W_sub2_mass_opt1_"+rgn[ij], "", 20, 0.0, 300 );
				h_opt1_W_sub1_btag_opt1[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_W_sub1_btag_opt1_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt1_W_sub2_btag_opt1[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_W_sub2_btag_opt1_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt1_dR_lW_opt1[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_dR_lW_opt1_"+rgn[ij], "", 50, 0.0, 5.0 );
				h_opt1_dy_lW_opt1[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_dy_lW_opt1_"+rgn[ij], "", 50, -5.0, 5.0 );
				h_opt1_dphi_lW_opt1[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_dphi_lW_opt1_"+rgn[ij], "", 60, -M_PI, +M_PI );
				h_opt1_H_mass_opt1[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_H_mass_opt1_"+rgn[ij], "", 35, 0.0, 350.0 );
				h_opt1_X_mass_opt1[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_X_mass_opt1_"+rgn[ij], "", 40 ,0.0, 4000.0 );
				h_opt1_nbjets_other[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_nbjets_other_"+rgn[ij], "", 5 ,0.0, 5.0 );
				
                                //*******************************************************************************//
				//*******************************************************************************//

				h_opt2_l_pt[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_l_pt_"+rgn[ij], "",40,0,1000);
				h_opt2_l_eta[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_l_eta_"+rgn[ij], "",50,-2.5,2.5);
	       	                h_opt2_MET_pt[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_MET_pt_"+rgn[ij], "",40,0,1000);
				h_opt2_MET_sig[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_MET_sig_"+rgn[ij], "",20,0,300);
				h_opt2_MT[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_MT_"+rgn[ij], "",20,0,300);
				h_opt2_leadbjet_pt[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_LeadBJet_pt_"+rgn[ij], "",40,0,1000);
				h_opt2_leadbjet_btag_DeepFlav[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_LeadBJet_btag_DeepFlav_"+rgn[ij], "",50,0,1);
				h_opt2_Y_pt[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_Y_pt_"+rgn[ij], "",XBINS, XEDGES );
				h_opt2_Y_msoftdrop[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_Y_msoftdrop_"+rgn[ij], "", 20, 0.0, 500.0 );
				h_opt2_Y_PNetMD_XbbvsQCD[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_Y_PNetMD_XbbvsQCD_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt2_Y_PNetMD_WvsQCD[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_Y_PNetMD_WvsQCD_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt2_Y_PNet_TvsQCD[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_Y_PNet_TvsQCD_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt2_Y_sub1_mass[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_Y_sub1_mass_"+rgn[ij], "", 20, 0.0, 300 );
				h_opt2_Y_sub2_mass[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_Y_sub2_mass_"+rgn[ij], "", 20, 0.0, 300 );
				h_opt2_Y_sub1_btag[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_Y_sub1_btag_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt2_Y_sub2_btag[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_Y_sub2_btag_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt2_W_pt_opt2[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_W_pt_opt2_"+rgn[ij], "",XBINS, XEDGES );
				h_opt2_W_msoftdrop_opt2[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_W_msoftdrop_opt2_"+rgn[ij], "", 35, 0.0, 350.0 );
				h_opt2_W_PNetMD_XbbvsQCD_opt2[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_W_PNetMD_XbbvsQCD_opt2_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt2_W_PNetMD_WvsQCD_opt2[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_W_PNetMD_WvsQCD_opt2_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt2_W_PNet_TvsQCD_opt2[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_W_PNet_TvsQCD_opt2_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt2_W_DAK8MD_WvsQCD_opt2[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_W_DAK8MD_WvsQCD_opt2_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt2_W_sub1_mass_opt2[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_W_sub1_mass_opt2_"+rgn[ij], "", 20, 0.0, 300 );
				h_opt2_W_sub2_mass_opt2[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_W_sub2_mass_opt2_"+rgn[ij], "", 20, 0.0, 300 );
				h_opt2_W_sub1_btag_opt2[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_W_sub1_btag_opt2_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt2_W_sub2_btag_opt2[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_W_sub2_btag_opt2_"+rgn[ij], "", 50, 0.0, 1.0 );
				h_opt2_dR_lW_opt2[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_dR_lW_opt2_"+rgn[ij], "", 50, 0.0, 5.0 );
				h_opt2_dy_lW_opt2[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_dy_lW_opt2_"+rgn[ij], "", 50, -5.0, 5.0 );
				h_opt2_dphi_lW_opt2[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_dphi_lW_opt2_"+rgn[ij], "", 60, -M_PI, +M_PI );
				h_opt2_H_mass_opt2[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_H_mass_opt2_"+rgn[ij], "", 35, 0.0, 350.0 );
				h_opt2_X_mass_opt2[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_X_mass_opt2_"+rgn[ij], "", 40 ,0.0, 4000.0 );
				h_opt2_nbjets_other[pc] = new TH1F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_nbjets_other_"+rgn[ij], "", 5 ,0.0, 5.0 );
				h_opt1_X_Y_mass[pc] = new TH2F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt1_X_Y_mass_"+rgn[ij], "", 40, 0.0, 4000.0, 20, 0.0, 500.0);
				h_opt2_X_Y_mass[pc] = new TH2F("h_Y_"+Ytype[yt]+"_W_"+Wtype[wt]+"_opt2_X_Y_mass_"+rgn[ij], "", 40, 0.0, 4000.0, 20, 0.0, 500.0);
                   
                pc++;
	      	}
       }
  }
   cout << pc << endl;
 
   Long64_t nn = tree->GetEntries();
   for(Long64_t j =0; j < nn ; j++)
   {
   tree->GetEntry(j);
   if( j % 10000 == 0) { std::cout << " 10000 events processed" << std::endl;}

   bool logic_opt1[nrgn*9];
   bool logic_opt2[nrgn*9];
   //if (nbjets_other > 0)
   if(l_minisoall < 0.1)
   {
   logic_opt1[0] = (Flag_Y_bb_pass_T && Flag_H_W_pass_L_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[1] = (Flag_Y_bb_pass_T && Flag_H_W_pass_M_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[2] = (Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[3] = (Flag_Y_bb_pass_M && Flag_H_W_pass_L_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[4] = (Flag_Y_bb_pass_M && Flag_H_W_pass_M_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[5] = (Flag_Y_bb_pass_M && Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);   
   logic_opt1[6] = (Flag_Y_bb_pass_L && Flag_H_W_pass_L_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[7] = (Flag_Y_bb_pass_L && Flag_H_W_pass_M_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[8] = (Flag_Y_bb_pass_L && Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);

   logic_opt1[9] =  (!Flag_Y_bb_pass_T && Flag_H_W_pass_L_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[10] = (!Flag_Y_bb_pass_T && Flag_H_W_pass_M_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[11] = (!Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[12] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_L_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[13] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_M_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[14] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[15] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_L_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[16] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_M_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[17] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass);

   logic_opt1[18] = (!Flag_Y_bb_pass_T && Flag_H_W_pass_L_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[19] = (!Flag_Y_bb_pass_T && Flag_H_W_pass_M_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[20] = (!Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[21] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_L_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[22] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_M_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[23] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[24] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_L_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[25] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_M_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[26] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);

   logic_opt1[27] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_L_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[28] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_M_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[29] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[30] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_L_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[31] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_M_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[32] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[33] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_L_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[34] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_M_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[35] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);

   logic_opt1[36] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_L_opt1 && Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   logic_opt1[37] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_M_opt1 && Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   logic_opt1[38] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   logic_opt1[39] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_L_opt1 && Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   logic_opt1[40] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_M_opt1 && Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   logic_opt1[41] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   logic_opt1[42] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_L_opt1 && Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   logic_opt1[43] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_M_opt1 && Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   logic_opt1[44] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && !Flag_MET_pass);

   logic_opt1[45] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_L_opt1 && !Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   logic_opt1[46] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_M_opt1 && !Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   logic_opt1[47] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   logic_opt1[48] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_L_opt1 && !Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   logic_opt1[49] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_M_opt1 && !Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   logic_opt1[50] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   logic_opt1[51] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_L_opt1 && !Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   logic_opt1[52] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_M_opt1 && !Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   logic_opt1[53] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && !Flag_MET_pass);
   
   logic_opt1[54] = (Flag_Y_bb_pass_T && !Flag_H_W_pass_L_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[55] = (Flag_Y_bb_pass_T && !Flag_H_W_pass_M_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[56] = (Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[57] = (Flag_Y_bb_pass_M && !Flag_H_W_pass_L_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[58] = (Flag_Y_bb_pass_M && !Flag_H_W_pass_M_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[59] = (Flag_Y_bb_pass_M && !Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[60] = (Flag_Y_bb_pass_L && !Flag_H_W_pass_L_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[61] = (Flag_Y_bb_pass_L && !Flag_H_W_pass_M_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[62] = (Flag_Y_bb_pass_L && !Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);

   logic_opt1[63] = (!Flag_Y_bb_pass_T && Flag_H_W_pass_L_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[64] = (!Flag_Y_bb_pass_T && Flag_H_W_pass_M_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[65] = (!Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[66] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_L_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[67] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_M_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[68] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[69] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_L_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[70] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_M_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);
   logic_opt1[71] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass);

   logic_opt2[0] = (Flag_Y_bb_pass_T && Flag_H_W_pass_L_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[1] = (Flag_Y_bb_pass_T && Flag_H_W_pass_M_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[2] = (Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[3] = (Flag_Y_bb_pass_M && Flag_H_W_pass_L_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[4] = (Flag_Y_bb_pass_M && Flag_H_W_pass_M_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[5] = (Flag_Y_bb_pass_M && Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);   
   logic_opt2[6] = (Flag_Y_bb_pass_L && Flag_H_W_pass_L_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[7] = (Flag_Y_bb_pass_L && Flag_H_W_pass_M_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[8] = (Flag_Y_bb_pass_L && Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);

   logic_opt2[9] =  (!Flag_Y_bb_pass_T && Flag_H_W_pass_L_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[10] = (!Flag_Y_bb_pass_T && Flag_H_W_pass_M_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[11] = (!Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[12] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_L_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[13] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_M_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[14] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[15] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_L_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[16] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_M_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[17] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass);

   logic_opt2[18] = (!Flag_Y_bb_pass_T && Flag_H_W_pass_L_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[19] = (!Flag_Y_bb_pass_T && Flag_H_W_pass_M_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[20] = (!Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[21] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_L_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[22] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_M_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[23] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[24] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_L_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[25] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_M_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[26] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);

   logic_opt2[27] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_L_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[28] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_M_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[29] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[30] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_L_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[31] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_M_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[32] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[33] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_L_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[34] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_M_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[35] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);

   logic_opt2[36] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_L_opt2 && Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   logic_opt2[37] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_M_opt2 && Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   logic_opt2[38] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   logic_opt2[39] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_L_opt2 && Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   logic_opt2[40] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_M_opt2 && Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   logic_opt2[41] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   logic_opt2[42] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_L_opt2 && Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   logic_opt2[43] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_M_opt2 && Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   logic_opt2[44] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && !Flag_MET_pass);

   logic_opt2[45] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_L_opt2 && !Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   logic_opt2[46] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_M_opt2 && !Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   logic_opt2[47] = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   logic_opt2[48] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_L_opt2 && !Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   logic_opt2[49] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_M_opt2 && !Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   logic_opt2[50] = (!Flag_Y_bb_pass_M && !Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   logic_opt2[51] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_L_opt2 && !Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   logic_opt2[52] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_M_opt2 && !Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   logic_opt2[53] = (!Flag_Y_bb_pass_L && !Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && !Flag_MET_pass);
   
   logic_opt2[54] = (Flag_Y_bb_pass_T && !Flag_H_W_pass_L_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[55] = (Flag_Y_bb_pass_T && !Flag_H_W_pass_M_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[56] = (Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[57] = (Flag_Y_bb_pass_M && !Flag_H_W_pass_L_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[58] = (Flag_Y_bb_pass_M && !Flag_H_W_pass_M_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[59] = (Flag_Y_bb_pass_M && !Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[60] = (Flag_Y_bb_pass_L && !Flag_H_W_pass_L_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[61] = (Flag_Y_bb_pass_L && !Flag_H_W_pass_M_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[62] = (Flag_Y_bb_pass_L && !Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);

   logic_opt2[63] = (!Flag_Y_bb_pass_T && Flag_H_W_pass_L_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[64] = (!Flag_Y_bb_pass_T && Flag_H_W_pass_M_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[65] = (!Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[66] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_L_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[67] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_M_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[68] = (!Flag_Y_bb_pass_M && Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[69] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_L_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[70] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_M_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   logic_opt2[71] = (!Flag_Y_bb_pass_L && Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass);
   
   int t=0;
   float weight_nom;
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
   
   float MT = sqrt(2*l_pt*MET_pt*(1-cos(PhiInRange(l_phi-MET_phi))));
   
   for (int i=0 ; i< nrgn ; i++)
   {
	   
      for (int yt=0 ; yt< 3 ; yt++)
       {
           for (int wt=0 ; wt< 3 ; wt++)
               {
                  if(logic_opt1[t]) {  // probably, it is logic_opt1
					
					h_opt1_l_pt[t]->Fill(l_pt,weight_nom); 
					h_opt1_l_eta[t]->Fill(l_eta,weight_nom); 
					
 			        	h_opt1_MET_pt[t]->Fill(MET_pt,weight_nom); 
					h_opt1_MET_sig[t]->Fill(MET_sig,weight_nom); 
					h_opt1_MT[t]->Fill(MT,weight_nom); 
	  
					if(nBJetAK4>0){
						h_opt1_leadbjet_pt[t]->Fill(BJetAK4_pt[0],weight_nom); 
						h_opt1_leadbjet_btag_DeepFlav[t]->Fill(BJetAK4_btag_DeepFlav[0],weight_nom); 
					}	
					  
					h_opt1_Y_pt[t]->Fill(Y_pt,weight_nom); 
					h_opt1_Y_msoftdrop[t]->Fill(Y_msoftdrop,weight_nom);
					h_opt1_Y_PNetMD_XbbvsQCD[t]->Fill(Y_DeepTag_PNetMD_XbbvsQCD,weight_nom);
					h_opt1_Y_PNetMD_WvsQCD[t]->Fill(Y_DeepTag_PNetMD_WvsQCD,weight_nom);
					h_opt1_Y_PNet_TvsQCD[t]->Fill(Y_DeepTag_PNet_TvsQCD,weight_nom);
					h_opt1_Y_sub1_mass[t]->Fill(Y_sub1_mass,weight_nom);
					h_opt1_Y_sub2_mass[t]->Fill(Y_sub2_mass,weight_nom);
					h_opt1_Y_sub1_btag[t]->Fill(Y_sub1_btag,weight_nom);
					h_opt1_Y_sub2_btag[t]->Fill(Y_sub2_btag,weight_nom);
                   
					h_opt1_W_pt_opt1[t]->Fill(W_pt_opt1,weight_nom); 
					h_opt1_W_msoftdrop_opt1[t]->Fill(W_msoftdrop_opt1,weight_nom); 
					h_opt1_W_PNetMD_XbbvsQCD_opt1[t]->Fill(W_DeepTag_PNetMD_XbbvsQCD_opt1,weight_nom); 
					h_opt1_W_PNetMD_WvsQCD_opt1[t]->Fill(W_DeepTag_PNetMD_WvsQCD_opt1,weight_nom); 
					h_opt1_W_PNet_TvsQCD_opt1[t]->Fill(W_DeepTag_PNet_TvsQCD_opt1,weight_nom); 
					h_opt1_W_DAK8MD_WvsQCD_opt1[t]->Fill(W_DeepTag_DAK8MD_WvsQCD_opt1,weight_nom); 
					h_opt1_W_sub1_mass_opt1[t]->Fill(W_sub1_mass_opt1,weight_nom); 
					h_opt1_W_sub2_mass_opt1[t]->Fill(W_sub2_mass_opt1,weight_nom); 
					h_opt1_W_sub1_btag_opt1[t]->Fill(W_sub1_btag_opt1,weight_nom); 
					h_opt1_W_sub2_btag_opt1[t]->Fill(W_sub2_btag_opt1,weight_nom); 
					h_opt1_dR_lW_opt1[t]->Fill(dR_lW_opt1,weight_nom); 
					h_opt1_dy_lW_opt1[t]->Fill(dy_lW_opt1,weight_nom); 
					h_opt1_dphi_lW_opt1[t]->Fill(dphi_lW_opt1,weight_nom); 
				
					h_opt1_H_mass_opt1[t]->Fill(H_mass_opt1,weight_nom); 
					h_opt1_X_mass_opt1[t]->Fill(X_mass_opt1,weight_nom); 
					
					h_opt1_nbjets_other[t]->Fill(nbjets_other,weight_nom);
					
					h_opt1_X_Y_mass[t]->Fill(X_mass_opt1, Y_msoftdrop, weight_nom);
				   
		                    }
                   if(logic_opt2[t]) {
                   
                                        h_opt2_l_pt[t]->Fill(l_pt,weight_nom); 
					h_opt2_l_eta[t]->Fill(l_eta,weight_nom); 
                    
					h_opt2_MET_pt[t]->Fill(MET_pt,weight_nom); 
					h_opt2_MET_sig[t]->Fill(MET_sig,weight_nom); 
					h_opt2_MT[t]->Fill(MT,weight_nom); 
	  
					if(nBJetAK4>0){
						h_opt2_leadbjet_pt[t]->Fill(BJetAK4_pt[0],weight_nom); 
						h_opt2_leadbjet_btag_DeepFlav[t]->Fill(BJetAK4_btag_DeepFlav[0],weight_nom); 
					}	
                   
					h_opt2_Y_pt[t]->Fill(Y_pt,weight_nom);
					h_opt2_Y_msoftdrop[t]->Fill(Y_msoftdrop,weight_nom);
					h_opt2_Y_PNetMD_XbbvsQCD[t]->Fill(Y_DeepTag_PNetMD_XbbvsQCD,weight_nom);
					h_opt2_Y_PNetMD_WvsQCD[t]->Fill(Y_DeepTag_PNetMD_WvsQCD,weight_nom);
					h_opt2_Y_PNet_TvsQCD[t]->Fill(Y_DeepTag_PNet_TvsQCD,weight_nom);
					h_opt2_Y_sub1_mass[t]->Fill(Y_sub1_mass,weight_nom);
					h_opt2_Y_sub2_mass[t]->Fill(Y_sub2_mass,weight_nom);
					h_opt2_Y_sub1_btag[t]->Fill(Y_sub1_btag,weight_nom);
					h_opt2_Y_sub2_btag[t]->Fill(Y_sub2_btag,weight_nom);
                   
					h_opt2_W_pt_opt2[t]->Fill(W_pt_opt2,weight_nom);
					h_opt2_W_msoftdrop_opt2[t]->Fill(W_msoftdrop_opt2,weight_nom);
					h_opt2_W_PNetMD_XbbvsQCD_opt2[t]->Fill(W_DeepTag_PNetMD_XbbvsQCD_opt2,weight_nom); 
					h_opt2_W_PNetMD_WvsQCD_opt2[t]->Fill(W_DeepTag_PNetMD_WvsQCD_opt2,weight_nom); 
					h_opt2_W_PNet_TvsQCD_opt2[t]->Fill(W_DeepTag_PNet_TvsQCD_opt2,weight_nom); 
					h_opt2_W_DAK8MD_WvsQCD_opt2[t]->Fill(W_DeepTag_DAK8MD_WvsQCD_opt2,weight_nom); 
					h_opt2_W_sub1_mass_opt2[t]->Fill(W_sub1_mass_opt2,weight_nom); 
					h_opt2_W_sub2_mass_opt2[t]->Fill(W_sub2_mass_opt2,weight_nom); 
					h_opt2_W_sub1_btag_opt2[t]->Fill(W_sub1_btag_opt2,weight_nom); 
					h_opt2_W_sub2_btag_opt2[t]->Fill(W_sub2_btag_opt2,weight_nom); 
					h_opt2_dR_lW_opt2[t]->Fill(dR_lW_opt2,weight_nom); 
					h_opt2_dy_lW_opt2[t]->Fill(dy_lW_opt2,weight_nom); 
					h_opt2_dphi_lW_opt2[t]->Fill(dphi_lW_opt2,weight_nom); 
                   
					h_opt2_H_mass_opt2[t]->Fill(H_mass_opt2,weight_nom);
					h_opt2_X_mass_opt2[t]->Fill(X_mass_opt2,weight_nom);
					h_opt2_nbjets_other[t]->Fill(nbjets_other,weight_nom);
					h_opt2_X_Y_mass[t]->Fill(X_mass_opt2, Y_msoftdrop, weight_nom);
                   
                                     }
                   t++;
      	       }	  
	    }
       }
   }
   }
    final_file->Write();
    final_file->cd();
    int z=0;
    for (int i=0 ; i< nrgn ; i++)
    {
      for (int yt=0 ; yt< 3 ; yt++)
       {
           for (int wt=0 ; wt< 3 ; wt++)
               {
                    h_opt1_l_pt[z]->Write();
                    h_opt1_l_eta[z]->Write();
                    h_opt1_MET_pt[z]->Write();
                    h_opt1_MET_pt[z]->Write();
                    h_opt1_MET_sig[z]->Write();
                    h_opt1_MT[z]->Write();
                    h_opt1_leadbjet_pt[z]->Write();
                    h_opt1_leadbjet_btag_DeepFlav[z]->Write();

                    h_opt1_Y_pt[z]->Write();
                    h_opt1_Y_msoftdrop[z]->Write();
                    h_opt1_Y_PNetMD_XbbvsQCD[z]->Write();
                    h_opt1_Y_PNetMD_WvsQCD[z]->Write();
                    h_opt1_Y_PNet_TvsQCD[z]->Write();
                    h_opt1_Y_sub1_mass[z]->Write();
                    h_opt1_Y_sub2_mass[z]->Write();
                    h_opt1_Y_sub1_btag[z]->Write();
                    h_opt1_Y_sub2_btag[z]->Write();

                    h_opt1_W_pt_opt1[z]->Write();
                    h_opt1_W_msoftdrop_opt1[z]->Write();
                    h_opt1_W_PNetMD_XbbvsQCD_opt1[z]->Write();
                    h_opt1_W_PNetMD_WvsQCD_opt1[z]->Write();
                    h_opt1_W_PNet_TvsQCD_opt1[z]->Write();
                    h_opt1_W_DAK8MD_WvsQCD_opt1[z]->Write();
                    h_opt1_W_sub1_mass_opt1[z]->Write();
                    h_opt1_W_sub2_mass_opt1[z]->Write();
                    h_opt1_W_sub1_btag_opt1[z]->Write();
                    h_opt1_W_sub2_btag_opt1[z]->Write();
                    h_opt1_dR_lW_opt1[z]->Write();
                    h_opt1_dy_lW_opt1[z]->Write();
                    h_opt1_dphi_lW_opt1[z]->Write();

                    h_opt1_H_mass_opt1[z]->Write();
                    h_opt1_X_mass_opt1[z]->Write();
                    h_opt1_nbjets_other[z]->Write();
                    h_opt1_X_Y_mass[z]->Write();

                    h_opt2_l_pt[z]->Write();
                    h_opt2_l_eta[z]->Write();
                    h_opt2_MET_pt[z]->Write();
                    h_opt2_MET_sig[z]->Write();
                    h_opt2_MT[z]->Write();
                    h_opt2_leadbjet_pt[z]->Write();
                    h_opt2_leadbjet_btag_DeepFlav[z]->Write();

                    h_opt2_Y_pt[z]->Write();
                    h_opt2_Y_msoftdrop[z]->Write();
                    h_opt2_Y_PNetMD_XbbvsQCD[z]->Write();
                    h_opt2_Y_PNetMD_WvsQCD[z]->Write();
                    h_opt2_Y_PNet_TvsQCD[z]->Write();
                    h_opt2_Y_sub1_mass[z]->Write();
                    h_opt2_Y_sub2_mass[z]->Write();
                    h_opt2_Y_sub1_btag[z]->Write();
                    h_opt2_Y_sub2_btag[z]->Write();

                    h_opt2_W_pt_opt2[z]->Write();
                    h_opt2_W_msoftdrop_opt2[z]->Write();
                    h_opt2_W_PNetMD_XbbvsQCD_opt2[z]->Write();
                    h_opt2_W_PNetMD_WvsQCD_opt2[z]->Write();
                    h_opt2_W_PNet_TvsQCD_opt2[z]->Write();
                    h_opt2_W_DAK8MD_WvsQCD_opt2[z]->Write();
                    h_opt2_W_sub1_mass_opt2[z]->Write();
                    h_opt2_W_sub2_mass_opt2[z]->Write();
                    h_opt2_W_sub1_btag_opt2[z]->Write();
                    h_opt2_W_sub2_btag_opt2[z]->Write();
                    h_opt2_dR_lW_opt2[z]->Write();
                    h_opt2_dy_lW_opt2[z]->Write();
                    h_opt2_dphi_lW_opt2[z]->Write();

                    h_opt2_H_mass_opt2[z]->Write();
                    h_opt2_X_mass_opt2[z]->Write();
                    h_opt2_nbjets_other[z]->Write();
                    h_opt2_X_Y_mass[z]->Write();
		    z++;
               }
           }
     }
    final_file->Close();
  }
}

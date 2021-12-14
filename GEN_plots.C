//#include "Anal_XtoYH.h"
#include "getobjects.h"
#include<iostream>
#include<vector>
#include <bits/stdc++.h>
#include<fstream>
#include<string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

float* Muon_SF(TFile *file_mu_sf, string id, float pt, float eta){

	char name[100];
	
	sprintf(name,"NUM_%sID_DEN_TrackerMuons_abseta_pt",id.c_str());
	TH2F *h_SF = (TH2F*)file_mu_sf->Get(name);
	sprintf(name,"NUM_%sID_DEN_TrackerMuons_abseta_pt_stat",id.c_str());
	TH2F *h_SF_stat = (TH2F*)file_mu_sf->Get(name);
	sprintf(name,"NUM_%sID_DEN_TrackerMuons_abseta_pt_syst",id.c_str());
	TH2F *h_SF_sys = (TH2F*)file_mu_sf->Get(name);
	
	int eta_bin_id = h_SF->GetXaxis()->FindBin(fabs(eta));
	int pt_bin_id = h_SF->GetYaxis()->FindBin(pt);
	
	float sf, sf_stat, sf_sys;
		
	if(eta_bin_id>0 && eta_bin_id<=(h_SF->GetNbinsX()) && pt_bin_id>0 && pt_bin_id<=(h_SF->GetNbinsY())){
		sf = h_SF->GetBinContent(eta_bin_id,pt_bin_id);
		sf_stat = h_SF_stat->GetBinContent(eta_bin_id,pt_bin_id);
		sf_sys = h_SF_sys->GetBinContent(eta_bin_id,pt_bin_id);
	}else{
		sf = 1;
		sf_stat = sf_sys = 0;
		}
	
	static float sfvalues[3];
	sfvalues[0] = sf;
	sfvalues[1] = sf_stat;
	sfvalues[2] = sf_stat;
		
	return sfvalues;
}

float* Electron_SF(TFile *file_el_sf, float pt, float eta){

	char name[100];

	TH2F *h_SF = (TH2F*)file_el_sf->Get("EGamma_SF2D");
	TH2F *h_SF_statData = (TH2F*)file_el_sf->Get("statData");
	TH2F *h_SF_statMC = (TH2F*)file_el_sf->Get("statMC");
	TH2F *h_SF_altBkgModel = (TH2F*)file_el_sf->Get("altBkgModel");
	TH2F *h_SF_altSignalModel = (TH2F*)file_el_sf->Get("altSignalModel");
	TH2F *h_SF_altMCEff = (TH2F*)file_el_sf->Get("altMCEff");
	TH2F *h_SF_altTagSelection = (TH2F*)file_el_sf->Get("altTagSelection");
	
	int eta_bin_id = h_SF->GetXaxis()->FindBin(fabs(eta));
	int pt_bin_id = h_SF->GetYaxis()->FindBin(pt);
	
	float sf, sf_stat, sf_sys;
	
	static float sfvalues[7] = {-100,-100,-100,-100,-100,-100,-100};
	
	if(eta_bin_id>0 && eta_bin_id<=(h_SF->GetNbinsX()) && pt_bin_id>0 && pt_bin_id<=(h_SF->GetNbinsY())){
		sfvalues[0] = h_SF->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[1] = h_SF_statData->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[2] = h_SF_statMC->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[3] = h_SF_altBkgModel->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[4] = h_SF_altSignalModel->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[5] = h_SF_altMCEff->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[6] = h_SF_altTagSelection->GetBinContent(eta_bin_id,pt_bin_id);
	}
	else{
		sfvalues[0] = 1;
		sfvalues[1] = sfvalues[2] = sfvalues[3] = sfvalues[4] = sfvalues[5] = sfvalues[6] = 0;
		}
			
	return sfvalues;
}


//int main(int argc, char *argv[])
int main()
{
	
   

   char name[1000];
   
   char rootfiles[100];
   char infile[1000];
   char datafile[1000];
   char outfilx[100];

   cout <<"Give the input file name"<<endl;
   cin>> rootfiles;
   
   // declare output file & histograms //
   
   int len = strlen(rootfiles);
   strncpy(outfilx, rootfiles, len-4);
   outfilx[len-4]='\0';
   sprintf(name,"Output_GEN_%s.root",outfilx);
   TFile *fileout = new TFile(name,"recreate");
   
   // Histograms
   
   TH1F *h_LHE_b1b2_m = new TH1F("LHE_b1b2_M", "", 100, 195, 205);
   TH1F *h_LHE_j1j2_m = new TH1F("LHE_j1j2_M", "", 100, 0, 150);
   TH1F *h_LHE_H_m = new TH1F("LHE_H_M", "", 50, 124, 126);
   TH1F *h_LHE_X_m = new TH1F("LHE_X_M", "", 350, 1400, 2100);
   
   TH1F *h_LHE_b1b2_pt = new TH1F("LHE_b1b2_pt", "", 40, 0, 1000);
   TH1F *h_LHE_j1j2_pt = new TH1F("LHE_j1j2_pt", "", 40, 0, 1000);
   TH1F *h_LHE_l_pt = new TH1F("LHE_l_pt", "", 40, 0, 1000);
   TH1F *h_LHE_nu_pt = new TH1F("LHE_nu_pt", "", 40, 0, 1000);
   TH1F *h_LHE_H_pt = new TH1F("LHE_H_pt", "", 40, 0, 1000);
   
   TH1F *h_LHE_b1b2_dR = new TH1F("LHE_b1b2_dR", "", 100, 0, 5);
   TH1F *h_LHE_j1j2_dR = new TH1F("LHE_j1j2_dR", "", 100, 0, 5);
   TH1F *h_LHE_j1j2l_dR = new TH1F("LHE_j1j2l_dR", "", 100, 0, 5);
   TH1F *h_LHE_j1j2l_b1b2_dR = new TH1F("LHE_j1j2l_b1b2_dR", "", 100, 0, 5);
   
   TH1F *h_GEN_b1b2_m = new TH1F("GEN_b1b2_M", "", 100, 195, 205);
   TH1F *h_GEN_j1j2_m = new TH1F("GEN_j1j2_M", "", 100, 0, 150);
   TH1F *h_GEN_H_m = new TH1F("GEN_H_M", "", 50, 124, 126);
   TH1F *h_GEN_X_m = new TH1F("GEN_X_M", "", 350, 1400, 2100);
   
   TH1F *h_GEN_b1b2_pt = new TH1F("GEN_b1b2_pt", "", 40, 0, 1000);
   TH1F *h_GEN_j1j2_pt = new TH1F("GEN_j1j2_pt", "", 40, 0, 1000);
   TH1F *h_GEN_l_pt = new TH1F("GEN_l_pt", "", 40, 0, 1000);
   TH1F *h_GEN_nu_pt = new TH1F("GEN_nu_pt", "", 40, 0, 1000);
   TH1F *h_GEN_H_pt = new TH1F("GEN_H_pt", "", 40, 0, 1000);
   
   TH1F *h_GEN_b1b2_dR = new TH1F("GEN_b1b2_dR", "", 100, 0, 5);
   TH1F *h_GEN_j1j2_dR = new TH1F("GEN_j1j2_dR", "", 100, 0, 5);
   TH1F *h_GEN_j1j2l_dR = new TH1F("GEN_j1j2l_dR", "", 100, 0, 5);
   TH1F *h_GEN_j1j2l_b1b2_dR = new TH1F("GEN_j1j2l_b1b2_dR", "", 100, 0, 5);
   
   TH1F *h_GEN_Y_m = new TH1F("GEN_Y_M", "", 100, 0, 300);
   TH1F *h_GEN_Wh_m = new TH1F("GEN_Wh_M", "", 100, 0, 150);
   TH1F *h_GEN_HJ_m = new TH1F("GEN_HJ_M", "", 60, 50, 350);
   TH1F *h_GEN_XJ_m = new TH1F("GEN_XJ_M", "", 160, 1000, 2600);
   TH2F *h_GEN_Y_XJ_m = new TH2F("GEN_XJ_Y_M", "", 160, 1000, 2600, 100, 0, 300);
   TH1F *h_GEN_HJ_m_MET = new TH1F("GEN_MET_HJ_M", "", 60, 50, 350);
   TH1F *h_GEN_XJ_m_MET = new TH1F("GEN_MET_XJ_M", "", 160, 1000, 2600);
   TH2F *h_GEN_Y_XJ_m_MET = new TH2F("GEN_MET_XJ_Y_M", "", 160, 1000, 2600, 100, 0, 300);
   
   // end //

   ifstream file_db;
   file_db.open(rootfiles);

    while(!(file_db.eof())){

	file_db >> datafile;
    if(file_db.eof()) break;

    sprintf(infile, "%s", datafile);
    cout<<"infile "<<infile<<endl;

    TFile *fileIn = new TFile(infile,"read");

   TTree *fChain;
   fChain = (TTree*)fileIn->Get("Events");
   fChain->SetBranchAddress("irun", &irun, &b_irun);
   fChain->SetBranchAddress("ilumi", &ilumi, &b_ilumi);
   fChain->SetBranchAddress("ievt", &ievt, &b_ievt);
   fChain->SetBranchAddress("nprim", &nprim, &b_nprim);
   fChain->SetBranchAddress("npvert", &npvert, &b_npvert);
   fChain->SetBranchAddress("Rho", &Rho, &b_Rho);
   fChain->SetBranchAddress("trig_value", &trig_value, &b_trig_value);
   fChain->SetBranchAddress("hlt_IsoMu24", &hlt_IsoMu24, &b_hlt_IsoMu24);
   fChain->SetBranchAddress("hlt_Mu50", &hlt_Mu50, &b_hlt_Mu50);
   fChain->SetBranchAddress("hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", &hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165, &b_hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
   fChain->SetBranchAddress("hlt_Ele115_CaloIdVT_GsfTrkIdT", &hlt_Ele115_CaloIdVT_GsfTrkIdT, &b_hlt_Ele115_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("hlt_Ele40_WPTight_Gsf", &hlt_Ele40_WPTight_Gsf, &b_hlt_Ele40_WPTight_Gsf);
   fChain->SetBranchAddress("hlt_Ele32_WPTight_Gsf", &hlt_Ele32_WPTight_Gsf, &b_hlt_Ele32_WPTight_Gsf);
   fChain->SetBranchAddress("hlt_Ele28_eta2p1_WPTight_Gsf_HT150", &hlt_Ele28_eta2p1_WPTight_Gsf_HT150, &b_hlt_Ele28_eta2p1_WPTight_Gsf_HT150);
   fChain->SetBranchAddress("hlt_Mu37_Ele27_CaloIdL_MW", &hlt_Mu37_Ele27_CaloIdL_MW, &b_hlt_Mu37_Ele27_CaloIdL_MW);
   fChain->SetBranchAddress("hlt_Mu27_Ele37_CaloIdL_MW", &hlt_Mu27_Ele37_CaloIdL_MW, &b_hlt_Mu27_Ele37_CaloIdL_MW);
   fChain->SetBranchAddress("hlt_Mu37_TkMu27", &hlt_Mu37_TkMu27, &b_hlt_Mu37_TkMu27);
   fChain->SetBranchAddress("hlt_DoubleEle25_CaloIdL_MW", &hlt_DoubleEle25_CaloIdL_MW, &b_hlt_DoubleEle25_CaloIdL_MW);
   fChain->SetBranchAddress("hlt_AK8PFJet500", &hlt_AK8PFJet500, &b_hlt_AK8PFJet500);
   fChain->SetBranchAddress("hlt_PFJet500", &hlt_PFJet500, &b_hlt_PFJet500);
   fChain->SetBranchAddress("hlt_HT1050", &hlt_HT1050, &b_hlt_HT1050);
   fChain->SetBranchAddress("hlt_AK8PFJet400_TrimMass30", &hlt_AK8PFJet400_TrimMass30, &b_hlt_AK8PFJet400_TrimMass30);
   fChain->SetBranchAddress("hlt_AK8PFHT800_TrimMass50", &hlt_AK8PFHT800_TrimMass50, &b_hlt_AK8PFHT800_TrimMass50);
   fChain->SetBranchAddress("hlt_Photon200", &hlt_Photon200, &b_hlt_Photon200);
   fChain->SetBranchAddress("hlt_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &hlt_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, &b_hlt_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
   fChain->SetBranchAddress("hlt_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60", &hlt_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60, &b_hlt_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60);
   fChain->SetBranchAddress("hlt_PFMETNoMu140_PFMHTNoMu140_IDTight", &hlt_PFMETNoMu140_PFMHTNoMu140_IDTight, &b_hlt_PFMETNoMu140_PFMHTNoMu140_IDTight);
   fChain->SetBranchAddress("hlt_PFMETTypeOne140_PFMHT140_IDTight", &hlt_PFMETTypeOne140_PFMHT140_IDTight, &b_hlt_PFMETTypeOne140_PFMHT140_IDTight);
   fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   fChain->SetBranchAddress("TrigObj_pt", TrigObj_pt, &b_TrigObj_pt);
   fChain->SetBranchAddress("TrigObj_eta", TrigObj_eta, &b_TrigObj_eta);
   fChain->SetBranchAddress("TrigObj_phi", TrigObj_phi, &b_TrigObj_phi);
   fChain->SetBranchAddress("TrigObj_mass", TrigObj_mass, &b_TrigObj_mass);
   fChain->SetBranchAddress("TrigObj_HLT", TrigObj_HLT, &b_TrigObj_HLT);
   fChain->SetBranchAddress("TrigObj_L1", TrigObj_L1, &b_TrigObj_L1);
   fChain->SetBranchAddress("TrigObj_Ihlt", TrigObj_Ihlt, &b_TrigObj_Ihlt);
   fChain->SetBranchAddress("TrigObj_pdgId", TrigObj_pdgId, &b_TrigObj_pdgId);
   fChain->SetBranchAddress("TrigObj_type", TrigObj_type, &b_TrigObj_type);
   fChain->SetBranchAddress("CHSMET_pt", &CHSMET_pt, &b_miset);
   fChain->SetBranchAddress("CHSMET_phi", &CHSMET_phi, &b_misphi);
   fChain->SetBranchAddress("CHSMET_sig", &CHSMET_sig, &b_misetsig);
   fChain->SetBranchAddress("CHSMET_sumEt", &CHSMET_sumEt, &b_sumEt);
   fChain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_miset_PUPPI);
   fChain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_misphi_PUPPI);
   fChain->SetBranchAddress("PuppiMET_sig", &PuppiMET_sig, &b_misetsig_PUPPI);
   fChain->SetBranchAddress("PuppiMET_sumEt", &PuppiMET_sumEt, &b_sumEt_PUPPI);
   fChain->SetBranchAddress("PuppiMET_pt_JESup", &PuppiMET_pt_JESup, &b_miset_PUPPI_JESup);
   fChain->SetBranchAddress("PuppiMET_pt_JESdn", &PuppiMET_pt_JESdn, &b_miset_PUPPI_JESdn);
   fChain->SetBranchAddress("PuppiMET_pt_JERup", &PuppiMET_pt_JERup, &b_miset_PUPPI_JERup);
   fChain->SetBranchAddress("PuppiMET_pt_JERdn", &PuppiMET_pt_JERdn, &b_miset_PUPPI_JERdn);
   fChain->SetBranchAddress("PuppiMET_pt_UnclusEup", &PuppiMET_pt_UnclusEup, &b_miset_PUPPI_UnclusEup);
   fChain->SetBranchAddress("PuppiMET_pt_UnclusEdn", &PuppiMET_pt_UnclusEdn, &b_miset_PUPPI_UnclusEdn);
   fChain->SetBranchAddress("PuppiMET_phi_JESup", &PuppiMET_phi_JESup, &b_misphi_PUPPI_JESup);
   fChain->SetBranchAddress("PuppiMET_phi_JESdn", &PuppiMET_phi_JESdn, &b_misphi_PUPPI_JESdn);
   fChain->SetBranchAddress("PuppiMET_phi_JERup", &PuppiMET_phi_JERup, &b_misphi_PUPPI_JERup);
   fChain->SetBranchAddress("PuppiMET_phi_JERdn", &PuppiMET_phi_JERdn, &b_misphi_PUPPI_JERdn);
   fChain->SetBranchAddress("PuppiMET_phi_UnclusEup", &PuppiMET_phi_UnclusEup, &b_misphi_PUPPI_UnclusEup);
   fChain->SetBranchAddress("PuppiMET_phi_UnclusEdn", &PuppiMET_phi_UnclusEdn, &b_misphi_PUPPI_UnclusEdn);
   fChain->SetBranchAddress("nPFJetAK8", &nPFJetAK8, &b_nPFJetAK8);
   fChain->SetBranchAddress("PFJetAK8_pt", PFJetAK8_pt, &b_PFJetAK8_pt);
   fChain->SetBranchAddress("PFJetAK8_y", PFJetAK8_y, &b_PFJetAK8_y);
   fChain->SetBranchAddress("PFJetAK8_eta", PFJetAK8_eta, &b_PFJetAK8_eta);
   fChain->SetBranchAddress("PFJetAK8_phi", PFJetAK8_phi, &b_PFJetAK8_phi);
   fChain->SetBranchAddress("PFJetAK8_mass", PFJetAK8_mass, &b_PFJetAK8_mass);
   fChain->SetBranchAddress("PFJetAK8_jetID_tightlepveto", PFJetAK8_jetID_tightlepveto, &b_PFJetAK8_jetID_tightlepveto);
   fChain->SetBranchAddress("PFJetAK8_jetID", PFJetAK8_jetID, &b_PFJetAK8_jetID);
   fChain->SetBranchAddress("PFJetAK8_JEC", PFJetAK8_JEC, &b_PFJetAK8_JEC);
   fChain->SetBranchAddress("PFJetAK8_CHF", PFJetAK8_CHF, &b_PFJetAK8_CHF);
   fChain->SetBranchAddress("PFJetAK8_NHF", PFJetAK8_NHF, &b_PFJetAK8_NHF);
   fChain->SetBranchAddress("PFJetAK8_CEMF", PFJetAK8_CEMF, &b_PFJetAK8_CEMF);
   fChain->SetBranchAddress("PFJetAK8_NEMF", PFJetAK8_NEMF, &b_PFJetAK8_NEMF);
   fChain->SetBranchAddress("PFJetAK8_MUF", PFJetAK8_MUF, &b_PFJetAK8_MUF);
   fChain->SetBranchAddress("PFJetAK8_PHF", PFJetAK8_PHF, &b_PFJetAK8_PHF);
   fChain->SetBranchAddress("PFJetAK8_EEF", PFJetAK8_EEF, &b_PFJetAK8_EEF);
   fChain->SetBranchAddress("PFJetAK8_HFHF", PFJetAK8_HFHF, &b_PFJetAK8_HFHF);
   fChain->SetBranchAddress("PFJetAK8_CHM", PFJetAK8_CHM, &b_PFJetAK8_CHM);
   fChain->SetBranchAddress("PFJetAK8_NHM", PFJetAK8_NHM, &b_PFJetAK8_NHM);
   fChain->SetBranchAddress("PFJetAK8_MUM", PFJetAK8_MUM, &b_PFJetAK8_MUM);
   fChain->SetBranchAddress("PFJetAK8_PHM", PFJetAK8_PHM, &b_PFJetAK8_PHM);
   fChain->SetBranchAddress("PFJetAK8_EEM", PFJetAK8_EEM, &b_PFJetAK8_EEM);
   fChain->SetBranchAddress("PFJetAK8_HFHM", PFJetAK8_HFHM, &b_PFJetAK8_HFHM);
   fChain->SetBranchAddress("PFJetAK8_Neucons", PFJetAK8_Neucons, &b_PFJetAK8_Neucons);
   fChain->SetBranchAddress("PFJetAK8_Chcons", PFJetAK8_Chcons, &b_PFJetAK8_Chcons);
   fChain->SetBranchAddress("PFJetAK8_reso", PFJetAK8_reso, &b_PFJetAK8_reso);
   fChain->SetBranchAddress("PFJetAK8_resoup", PFJetAK8_resoup, &b_PFJetAK8_resoup);
   fChain->SetBranchAddress("PFJetAK8_resodn", PFJetAK8_resodn, &b_PFJetAK8_resodn);
   fChain->SetBranchAddress("PFJetAK8_sdmass", PFJetAK8_sdmass, &b_PFJetAK8_sdmass);
   fChain->SetBranchAddress("PFJetAK8_tau1", PFJetAK8_tau1, &b_PFJetAK8_tau1);
   fChain->SetBranchAddress("PFJetAK8_tau2", PFJetAK8_tau2, &b_PFJetAK8_tau2);
   fChain->SetBranchAddress("PFJetAK8_tau3", PFJetAK8_tau3, &b_PFJetAK8_tau3);
   fChain->SetBranchAddress("PFJetAK8_btag_DeepCSV", PFJetAK8_btag_DeepCSV, &b_PFJetAK8_btag_DeepCSV);
   fChain->SetBranchAddress("PFJetAK8_DeepTag_DAK8MD_TvsQCD", PFJetAK8_DeepTag_DAK8MD_TvsQCD, &b_PFJetAK8_DeepTag_DAK8MD_TvsQCD);
   fChain->SetBranchAddress("PFJetAK8_DeepTag_DAK8MD_WvsQCD", PFJetAK8_DeepTag_DAK8MD_WvsQCD, &b_PFJetAK8_DeepTag_DAK8MD_WvsQCD);
   fChain->SetBranchAddress("PFJetAK8_DeepTag_DAK8MD_ZvsQCD", PFJetAK8_DeepTag_DAK8MD_ZvsQCD, &b_PFJetAK8_DeepTag_DAK8MD_ZvsQCD);
   fChain->SetBranchAddress("PFJetAK8_DeepTag_DAK8MD_HvsQCD", PFJetAK8_DeepTag_DAK8MD_HvsQCD, &b_PFJetAK8_DeepTag_DAK8MD_HvsQCD);
   fChain->SetBranchAddress("PFJetAK8_DeepTag_DAK8MD_bbvsQCD", PFJetAK8_DeepTag_DAK8MD_bbvsQCD, &b_PFJetAK8_DeepTag_DAK8MD_bbvsQCD);
   fChain->SetBranchAddress("PFJetAK8_DeepTag_PNet_TvsQCD", PFJetAK8_DeepTag_PNet_TvsQCD, &b_PFJetAK8_DeepTag_PNet_TvsQCD);
   fChain->SetBranchAddress("PFJetAK8_DeepTag_PNet_WvsQCD", PFJetAK8_DeepTag_PNet_WvsQCD, &b_PFJetAK8_DeepTag_PNet_WvsQCD);
   fChain->SetBranchAddress("PFJetAK8_DeepTag_PNet_ZvsQCD", PFJetAK8_DeepTag_PNet_ZvsQCD, &b_PFJetAK8_DeepTag_PNet_ZvsQCD);
   fChain->SetBranchAddress("PFJetAK8_DeepTag_PNetMD_XbbvsQCD", PFJetAK8_DeepTag_PNetMD_XbbvsQCD, &b_PFJetAK8_DeepTag_PNetMD_XbbvsQCD);
   fChain->SetBranchAddress("PFJetAK8_DeepTag_PNetMD_XccvsQCD", PFJetAK8_DeepTag_PNetMD_XccvsQCD, &b_PFJetAK8_DeepTag_PNetMD_XccvsQCD);
   fChain->SetBranchAddress("PFJetAK8_DeepTag_PNetMD_XqqvsQCD", PFJetAK8_DeepTag_PNetMD_XqqvsQCD, &b_PFJetAK8_DeepTag_PNetMD_XqqvsQCD);
   fChain->SetBranchAddress("PFJetAK8_DeepTag_PNetMD_QCD", PFJetAK8_DeepTag_PNetMD_QCD, &b_PFJetAK8_DeepTag_PNetMD_QCD);
   fChain->SetBranchAddress("PFJetAK8_sub1pt", PFJetAK8_sub1pt, &b_PFJetAK8_sub1pt);
   fChain->SetBranchAddress("PFJetAK8_sub1eta", PFJetAK8_sub1eta, &b_PFJetAK8_sub1eta);
   fChain->SetBranchAddress("PFJetAK8_sub1phi", PFJetAK8_sub1phi, &b_PFJetAK8_sub1phi);
   fChain->SetBranchAddress("PFJetAK8_sub1mass", PFJetAK8_sub1mass, &b_PFJetAK8_sub1mass);
   fChain->SetBranchAddress("PFJetAK8_sub1btag", PFJetAK8_sub1btag, &b_PFJetAK8_sub1btag);
   fChain->SetBranchAddress("PFJetAK8_sub2pt", PFJetAK8_sub2pt, &b_PFJetAK8_sub2pt);
   fChain->SetBranchAddress("PFJetAK8_sub2eta", PFJetAK8_sub2eta, &b_PFJetAK8_sub2eta);
   fChain->SetBranchAddress("PFJetAK8_sub2phi", PFJetAK8_sub2phi, &b_PFJetAK8_sub2phi);
   fChain->SetBranchAddress("PFJetAK8_sub2mass", PFJetAK8_sub2mass, &b_PFJetAK8_sub2mass);
   fChain->SetBranchAddress("PFJetAK8_sub2btag", PFJetAK8_sub2btag, &b_PFJetAK8_sub2btag);
   fChain->SetBranchAddress("PFJetAK8_jesup_AbsoluteStat", PFJetAK8_jesup_AbsoluteStat, &b_PFJetAK8_jesup_AbsoluteStat);
   fChain->SetBranchAddress("PFJetAK8_jesup_AbsoluteScale", PFJetAK8_jesup_AbsoluteScale, &b_PFJetAK8_jesup_AbsoluteScale);
   fChain->SetBranchAddress("PFJetAK8_jesup_AbsoluteMPFBias", PFJetAK8_jesup_AbsoluteMPFBias, &b_PFJetAK8_jesup_AbsoluteMPFBias);
   fChain->SetBranchAddress("PFJetAK8_jesup_FlavorQCD", PFJetAK8_jesup_FlavorQCD, &b_PFJetAK8_jesup_FlavorQCD);
   fChain->SetBranchAddress("PFJetAK8_jesup_Fragmentation", PFJetAK8_jesup_Fragmentation, &b_PFJetAK8_jesup_Fragmentation);
   fChain->SetBranchAddress("PFJetAK8_jesup_PileUpDataMC", PFJetAK8_jesup_PileUpDataMC, &b_PFJetAK8_jesup_PileUpDataMC);
   fChain->SetBranchAddress("PFJetAK8_jesup_PileUpPtBB", PFJetAK8_jesup_PileUpPtBB, &b_PFJetAK8_jesup_PileUpPtBB);
   fChain->SetBranchAddress("PFJetAK8_jesup_PileUpPtEC1", PFJetAK8_jesup_PileUpPtEC1, &b_PFJetAK8_jesup_PileUpPtEC1);
   fChain->SetBranchAddress("PFJetAK8_jesup_PileUpPtEC2", PFJetAK8_jesup_PileUpPtEC2, &b_PFJetAK8_jesup_PileUpPtEC2);
   fChain->SetBranchAddress("PFJetAK8_jesup_PileUpPtRef", PFJetAK8_jesup_PileUpPtRef, &b_PFJetAK8_jesup_PileUpPtRef);
   fChain->SetBranchAddress("PFJetAK8_jesup_RelativeFSR", PFJetAK8_jesup_RelativeFSR, &b_PFJetAK8_jesup_RelativeFSR);
   fChain->SetBranchAddress("PFJetAK8_jesup_RelativeJEREC1", PFJetAK8_jesup_RelativeJEREC1, &b_PFJetAK8_jesup_RelativeJEREC1);
   fChain->SetBranchAddress("PFJetAK8_jesup_RelativeJEREC2", PFJetAK8_jesup_RelativeJEREC2, &b_PFJetAK8_jesup_RelativeJEREC2);
   fChain->SetBranchAddress("PFJetAK8_jesup_RelativePtBB", PFJetAK8_jesup_RelativePtBB, &b_PFJetAK8_jesup_RelativePtBB);
   fChain->SetBranchAddress("PFJetAK8_jesup_RelativePtEC1", PFJetAK8_jesup_RelativePtEC1, &b_PFJetAK8_jesup_RelativePtEC1);
   fChain->SetBranchAddress("PFJetAK8_jesup_RelativePtEC2", PFJetAK8_jesup_RelativePtEC2, &b_PFJetAK8_jesup_RelativePtEC2);
   fChain->SetBranchAddress("PFJetAK8_jesup_RelativeBal", PFJetAK8_jesup_RelativeBal, &b_PFJetAK8_jesup_RelativeBal);
   fChain->SetBranchAddress("PFJetAK8_jesup_RelativeSample", PFJetAK8_jesup_RelativeSample, &b_PFJetAK8_jesup_RelativeSample);
   fChain->SetBranchAddress("PFJetAK8_jesup_RelativeStatEC", PFJetAK8_jesup_RelativeStatEC, &b_PFJetAK8_jesup_RelativeStatEC);
   fChain->SetBranchAddress("PFJetAK8_jesup_RelativeStatFSR", PFJetAK8_jesup_RelativeStatFSR, &b_PFJetAK8_jesup_RelativeStatFSR);
   fChain->SetBranchAddress("PFJetAK8_jesup_SinglePionECAL", PFJetAK8_jesup_SinglePionECAL, &b_PFJetAK8_jesup_SinglePionECAL);
   fChain->SetBranchAddress("PFJetAK8_jesup_SinglePionHCAL", PFJetAK8_jesup_SinglePionHCAL, &b_PFJetAK8_jesup_SinglePionHCAL);
   fChain->SetBranchAddress("PFJetAK8_jesup_TimePtEta", PFJetAK8_jesup_TimePtEta, &b_PFJetAK8_jesup_TimePtEta);
   fChain->SetBranchAddress("PFJetAK8_jesup_Total", PFJetAK8_jesup_Total, &b_PFJetAK8_jesup_Total);
   fChain->SetBranchAddress("PFJetAK8_jesdn_AbsoluteStat", PFJetAK8_jesdn_AbsoluteStat, &b_PFJetAK8_jesdn_AbsoluteStat);
   fChain->SetBranchAddress("PFJetAK8_jesdn_AbsoluteScale", PFJetAK8_jesdn_AbsoluteScale, &b_PFJetAK8_jesdn_AbsoluteScale);
   fChain->SetBranchAddress("PFJetAK8_jesdn_AbsoluteMPFBias", PFJetAK8_jesdn_AbsoluteMPFBias, &b_PFJetAK8_jesdn_AbsoluteMPFBias);
   fChain->SetBranchAddress("PFJetAK8_jesdn_FlavorQCD", PFJetAK8_jesdn_FlavorQCD, &b_PFJetAK8_jesdn_FlavorQCD);
   fChain->SetBranchAddress("PFJetAK8_jesdn_Fragmentation", PFJetAK8_jesdn_Fragmentation, &b_PFJetAK8_jesdn_Fragmentation);
   fChain->SetBranchAddress("PFJetAK8_jesdn_PileUpDataMC", PFJetAK8_jesdn_PileUpDataMC, &b_PFJetAK8_jesdn_PileUpDataMC);
   fChain->SetBranchAddress("PFJetAK8_jesdn_PileUpPtBB", PFJetAK8_jesdn_PileUpPtBB, &b_PFJetAK8_jesdn_PileUpPtBB);
   fChain->SetBranchAddress("PFJetAK8_jesdn_PileUpPtEC1", PFJetAK8_jesdn_PileUpPtEC1, &b_PFJetAK8_jesdn_PileUpPtEC1);
   fChain->SetBranchAddress("PFJetAK8_jesdn_PileUpPtEC2", PFJetAK8_jesdn_PileUpPtEC2, &b_PFJetAK8_jesdn_PileUpPtEC2);
   fChain->SetBranchAddress("PFJetAK8_jesdn_PileUpPtRef", PFJetAK8_jesdn_PileUpPtRef, &b_PFJetAK8_jesdn_PileUpPtRef);
   fChain->SetBranchAddress("PFJetAK8_jesdn_RelativeFSR", PFJetAK8_jesdn_RelativeFSR, &b_PFJetAK8_jesdn_RelativeFSR);
   fChain->SetBranchAddress("PFJetAK8_jesdn_RelativeJEREC1", PFJetAK8_jesdn_RelativeJEREC1, &b_PFJetAK8_jesdn_RelativeJEREC1);
   fChain->SetBranchAddress("PFJetAK8_jesdn_RelativeJEREC2", PFJetAK8_jesdn_RelativeJEREC2, &b_PFJetAK8_jesdn_RelativeJEREC2);
   fChain->SetBranchAddress("PFJetAK8_jesdn_RelativePtBB", PFJetAK8_jesdn_RelativePtBB, &b_PFJetAK8_jesdn_RelativePtBB);
   fChain->SetBranchAddress("PFJetAK8_jesdn_RelativePtEC1", PFJetAK8_jesdn_RelativePtEC1, &b_PFJetAK8_jesdn_RelativePtEC1);
   fChain->SetBranchAddress("PFJetAK8_jesdn_RelativePtEC2", PFJetAK8_jesdn_RelativePtEC2, &b_PFJetAK8_jesdn_RelativePtEC2);
   fChain->SetBranchAddress("PFJetAK8_jesdn_RelativeBal", PFJetAK8_jesdn_RelativeBal, &b_PFJetAK8_jesdn_RelativeBal);
   fChain->SetBranchAddress("PFJetAK8_jesdn_RelativeSample", PFJetAK8_jesdn_RelativeSample, &b_PFJetAK8_jesdn_RelativeSample);
   fChain->SetBranchAddress("PFJetAK8_jesdn_RelativeStatEC", PFJetAK8_jesdn_RelativeStatEC, &b_PFJetAK8_jesdn_RelativeStatEC);
   fChain->SetBranchAddress("PFJetAK8_jesdn_RelativeStatFSR", PFJetAK8_jesdn_RelativeStatFSR, &b_PFJetAK8_jesdn_RelativeStatFSR);
   fChain->SetBranchAddress("PFJetAK8_jesdn_SinglePionECAL", PFJetAK8_jesdn_SinglePionECAL, &b_PFJetAK8_jesdn_SinglePionECAL);
   fChain->SetBranchAddress("PFJetAK8_jesdn_SinglePionHCAL", PFJetAK8_jesdn_SinglePionHCAL, &b_PFJetAK8_jesdn_SinglePionHCAL);
   fChain->SetBranchAddress("PFJetAK8_jesdn_TimePtEta", PFJetAK8_jesdn_TimePtEta, &b_PFJetAK8_jesdn_TimePtEta);
   fChain->SetBranchAddress("PFJetAK8_jesdn_Total", PFJetAK8_jesdn_Total, &b_PFJetAK8_jesdn_Total);
   fChain->SetBranchAddress("nPFJetAK4", &nPFJetAK4, &b_nPFJetAK4);
   fChain->SetBranchAddress("PFJetAK4_jetID", PFJetAK4_jetID, &b_PFJetAK4_jetID);
   fChain->SetBranchAddress("PFJetAK4_jetID_tightlepveto", PFJetAK4_jetID_tightlepveto, &b_PFJetAK4_jetID_tightlepveto);
   fChain->SetBranchAddress("PFJetAK4_pt", PFJetAK4_pt, &b_PFJetAK4_pt);
   fChain->SetBranchAddress("PFJetAK4_eta", PFJetAK4_eta, &b_PFJetAK4_eta);
   fChain->SetBranchAddress("PFJetAK4_y", PFJetAK4_y, &b_PFJetAK4_y);
   fChain->SetBranchAddress("PFJetAK4_phi", PFJetAK4_phi, &b_PFJetAK4_phi);
   fChain->SetBranchAddress("PFJetAK4_mass", PFJetAK4_mass, &b_PFJetAK4_mass);
   fChain->SetBranchAddress("PFJetAK4_JEC", PFJetAK4_JEC, &b_PFJetAK4_JEC);
   fChain->SetBranchAddress("PFJetAK4_btag_DeepCSV", PFJetAK4_btag_DeepCSV, &b_PFJetAK4_btag_DeepCSV);
   fChain->SetBranchAddress("PFJetAK4_btag_DeepFlav", PFJetAK4_btag_DeepFlav, &b_PFJetAK4_btag_DeepFlav);
   fChain->SetBranchAddress("PFJetAK4_reso", PFJetAK4_reso, &b_PFJetAK4_reso);
   fChain->SetBranchAddress("PFJetAK4_resoup", PFJetAK4_resoup, &b_PFJetAK4_resoup);
   fChain->SetBranchAddress("PFJetAK4_resodn", PFJetAK4_resodn, &b_PFJetAK4_resodn);
   fChain->SetBranchAddress("PFJetAK4_hadronflav", PFJetAK4_hadronflav, &b_PFJetAK4_hadronflav);
   fChain->SetBranchAddress("PFJetAK4_partonflav", PFJetAK4_partonflav, &b_PFJetAK4_partonflav);
   fChain->SetBranchAddress("PFJetAK4_qgl", PFJetAK4_qgl, &b_PFJetAK4_qgl);
   fChain->SetBranchAddress("PFJetAK4_PUID", PFJetAK4_PUID, &b_PFJetAK4_PUID);
   fChain->SetBranchAddress("PFJetAK4_jesup_AbsoluteStat", PFJetAK4_jesup_AbsoluteStat, &b_PFJetAK4_jesup_AbsoluteStat);
   fChain->SetBranchAddress("PFJetAK4_jesup_AbsoluteScale", PFJetAK4_jesup_AbsoluteScale, &b_PFJetAK4_jesup_AbsoluteScale);
   fChain->SetBranchAddress("PFJetAK4_jesup_AbsoluteMPFBias", PFJetAK4_jesup_AbsoluteMPFBias, &b_PFJetAK4_jesup_AbsoluteMPFBias);
   fChain->SetBranchAddress("PFJetAK4_jesup_FlavorQCD", PFJetAK4_jesup_FlavorQCD, &b_PFJetAK4_jesup_FlavorQCD);
   fChain->SetBranchAddress("PFJetAK4_jesup_Fragmentation", PFJetAK4_jesup_Fragmentation, &b_PFJetAK4_jesup_Fragmentation);
   fChain->SetBranchAddress("PFJetAK4_jesup_PileUpDataMC", PFJetAK4_jesup_PileUpDataMC, &b_PFJetAK4_jesup_PileUpDataMC);
   fChain->SetBranchAddress("PFJetAK4_jesup_PileUpPtBB", PFJetAK4_jesup_PileUpPtBB, &b_PFJetAK4_jesup_PileUpPtBB);
   fChain->SetBranchAddress("PFJetAK4_jesup_PileUpPtEC1", PFJetAK4_jesup_PileUpPtEC1, &b_PFJetAK4_jesup_PileUpPtEC1);
   fChain->SetBranchAddress("PFJetAK4_jesup_PileUpPtEC2", PFJetAK4_jesup_PileUpPtEC2, &b_PFJetAK4_jesup_PileUpPtEC2);
   fChain->SetBranchAddress("PFJetAK4_jesup_PileUpPtRef", PFJetAK4_jesup_PileUpPtRef, &b_PFJetAK4_jesup_PileUpPtRef);
   fChain->SetBranchAddress("PFJetAK4_jesup_RelativeFSR", PFJetAK4_jesup_RelativeFSR, &b_PFJetAK4_jesup_RelativeFSR);
   fChain->SetBranchAddress("PFJetAK4_jesup_RelativeJEREC1", PFJetAK4_jesup_RelativeJEREC1, &b_PFJetAK4_jesup_RelativeJEREC1);
   fChain->SetBranchAddress("PFJetAK4_jesup_RelativeJEREC2", PFJetAK4_jesup_RelativeJEREC2, &b_PFJetAK4_jesup_RelativeJEREC2);
   fChain->SetBranchAddress("PFJetAK4_jesup_RelativePtBB", PFJetAK4_jesup_RelativePtBB, &b_PFJetAK4_jesup_RelativePtBB);
   fChain->SetBranchAddress("PFJetAK4_jesup_RelativePtEC1", PFJetAK4_jesup_RelativePtEC1, &b_PFJetAK4_jesup_RelativePtEC1);
   fChain->SetBranchAddress("PFJetAK4_jesup_RelativePtEC2", PFJetAK4_jesup_RelativePtEC2, &b_PFJetAK4_jesup_RelativePtEC2);
   fChain->SetBranchAddress("PFJetAK4_jesup_RelativeBal", PFJetAK4_jesup_RelativeBal, &b_PFJetAK4_jesup_RelativeBal);
   fChain->SetBranchAddress("PFJetAK4_jesup_RelativeSample", PFJetAK4_jesup_RelativeSample, &b_PFJetAK4_jesup_RelativeSample);
   fChain->SetBranchAddress("PFJetAK4_jesup_RelativeStatEC", PFJetAK4_jesup_RelativeStatEC, &b_PFJetAK4_jesup_RelativeStatEC);
   fChain->SetBranchAddress("PFJetAK4_jesup_RelativeStatFSR", PFJetAK4_jesup_RelativeStatFSR, &b_PFJetAK4_jesup_RelativeStatFSR);
   fChain->SetBranchAddress("PFJetAK4_jesup_SinglePionECAL", PFJetAK4_jesup_SinglePionECAL, &b_PFJetAK4_jesup_SinglePionECAL);
   fChain->SetBranchAddress("PFJetAK4_jesup_SinglePionHCAL", PFJetAK4_jesup_SinglePionHCAL, &b_PFJetAK4_jesup_SinglePionHCAL);
   fChain->SetBranchAddress("PFJetAK4_jesup_TimePtEta", PFJetAK4_jesup_TimePtEta, &b_PFJetAK4_jesup_TimePtEta);
   fChain->SetBranchAddress("PFJetAK4_jesup_Total", PFJetAK4_jesup_Total, &b_PFJetAK4_jesup_Total);
   fChain->SetBranchAddress("PFJetAK4_jesdn_AbsoluteStat", PFJetAK4_jesdn_AbsoluteStat, &b_PFJetAK4_jesdn_AbsoluteStat);
   fChain->SetBranchAddress("PFJetAK4_jesdn_AbsoluteScale", PFJetAK4_jesdn_AbsoluteScale, &b_PFJetAK4_jesdn_AbsoluteScale);
   fChain->SetBranchAddress("PFJetAK4_jesdn_AbsoluteMPFBias", PFJetAK4_jesdn_AbsoluteMPFBias, &b_PFJetAK4_jesdn_AbsoluteMPFBias);
   fChain->SetBranchAddress("PFJetAK4_jesdn_FlavorQCD", PFJetAK4_jesdn_FlavorQCD, &b_PFJetAK4_jesdn_FlavorQCD);
   fChain->SetBranchAddress("PFJetAK4_jesdn_Fragmentation", PFJetAK4_jesdn_Fragmentation, &b_PFJetAK4_jesdn_Fragmentation);
   fChain->SetBranchAddress("PFJetAK4_jesdn_PileUpDataMC", PFJetAK4_jesdn_PileUpDataMC, &b_PFJetAK4_jesdn_PileUpDataMC);
   fChain->SetBranchAddress("PFJetAK4_jesdn_PileUpPtBB", PFJetAK4_jesdn_PileUpPtBB, &b_PFJetAK4_jesdn_PileUpPtBB);
   fChain->SetBranchAddress("PFJetAK4_jesdn_PileUpPtEC1", PFJetAK4_jesdn_PileUpPtEC1, &b_PFJetAK4_jesdn_PileUpPtEC1);
   fChain->SetBranchAddress("PFJetAK4_jesdn_PileUpPtEC2", PFJetAK4_jesdn_PileUpPtEC2, &b_PFJetAK4_jesdn_PileUpPtEC2);
   fChain->SetBranchAddress("PFJetAK4_jesdn_PileUpPtRef", PFJetAK4_jesdn_PileUpPtRef, &b_PFJetAK4_jesdn_PileUpPtRef);
   fChain->SetBranchAddress("PFJetAK4_jesdn_RelativeFSR", PFJetAK4_jesdn_RelativeFSR, &b_PFJetAK4_jesdn_RelativeFSR);
   fChain->SetBranchAddress("PFJetAK4_jesdn_RelativeJEREC1", PFJetAK4_jesdn_RelativeJEREC1, &b_PFJetAK4_jesdn_RelativeJEREC1);
   fChain->SetBranchAddress("PFJetAK4_jesdn_RelativeJEREC2", PFJetAK4_jesdn_RelativeJEREC2, &b_PFJetAK4_jesdn_RelativeJEREC2);
   fChain->SetBranchAddress("PFJetAK4_jesdn_RelativePtBB", PFJetAK4_jesdn_RelativePtBB, &b_PFJetAK4_jesdn_RelativePtBB);
   fChain->SetBranchAddress("PFJetAK4_jesdn_RelativePtEC1", PFJetAK4_jesdn_RelativePtEC1, &b_PFJetAK4_jesdn_RelativePtEC1);
   fChain->SetBranchAddress("PFJetAK4_jesdn_RelativePtEC2", PFJetAK4_jesdn_RelativePtEC2, &b_PFJetAK4_jesdn_RelativePtEC2);
   fChain->SetBranchAddress("PFJetAK4_jesdn_RelativeBal", PFJetAK4_jesdn_RelativeBal, &b_PFJetAK4_jesdn_RelativeBal);
   fChain->SetBranchAddress("PFJetAK4_jesdn_RelativeSample", PFJetAK4_jesdn_RelativeSample, &b_PFJetAK4_jesdn_RelativeSample);
   fChain->SetBranchAddress("PFJetAK4_jesdn_RelativeStatEC", PFJetAK4_jesdn_RelativeStatEC, &b_PFJetAK4_jesdn_RelativeStatEC);
   fChain->SetBranchAddress("PFJetAK4_jesdn_RelativeStatFSR", PFJetAK4_jesdn_RelativeStatFSR, &b_PFJetAK4_jesdn_RelativeStatFSR);
   fChain->SetBranchAddress("PFJetAK4_jesdn_SinglePionECAL", PFJetAK4_jesdn_SinglePionECAL, &b_PFJetAK4_jesdn_SinglePionECAL);
   fChain->SetBranchAddress("PFJetAK4_jesdn_SinglePionHCAL", PFJetAK4_jesdn_SinglePionHCAL, &b_PFJetAK4_jesdn_SinglePionHCAL);
   fChain->SetBranchAddress("PFJetAK4_jesdn_TimePtEta", PFJetAK4_jesdn_TimePtEta, &b_PFJetAK4_jesdn_TimePtEta);
   fChain->SetBranchAddress("PFJetAK4_jesdn_Total", PFJetAK4_jesdn_Total, &b_PFJetAK4_jesdn_Total);
   fChain->SetBranchAddress("PFJetAK4_btag_DeepCSV_SF", PFJetAK4_btag_DeepCSV_SF, &b_PFJetAK4_btag_DeepCSV_SF);
   fChain->SetBranchAddress("PFJetAK4_btag_DeepCSV_SF_up", PFJetAK4_btag_DeepCSV_SF_up, &b_PFJetAK4_btag_DeepCSV_SF_up);
   fChain->SetBranchAddress("PFJetAK4_btag_DeepCSV_SF_dn", PFJetAK4_btag_DeepCSV_SF_dn, &b_PFJetAK4_btag_DeepCSV_SF_dn);
   fChain->SetBranchAddress("PFJetAK4_btag_DeepFlav_SF", PFJetAK4_btag_DeepFlav_SF, &b_PFJetAK4_btag_DeepFlav_SF);
   fChain->SetBranchAddress("PFJetAK4_btag_DeepFlav_SF_up", PFJetAK4_btag_DeepFlav_SF_up, &b_PFJetAK4_btag_DeepFlav_SF_up);
   fChain->SetBranchAddress("PFJetAK4_btag_DeepFlav_SF_dn", PFJetAK4_btag_DeepFlav_SF_dn, &b_PFJetAK4_btag_DeepFlav_SF_dn);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("Muon_isPF", Muon_isPF, &b_Muon_isPF);
   fChain->SetBranchAddress("Muon_isGL", Muon_isGL, &b_Muon_isGL);
   fChain->SetBranchAddress("Muon_isTRK", Muon_isTRK, &b_Muon_isTRK);
   fChain->SetBranchAddress("Muon_isLoose", Muon_isLoose, &b_Muon_isLoose);
   fChain->SetBranchAddress("Muon_isGoodGL", Muon_isGoodGL, &b_Muon_isGoodGL);
   fChain->SetBranchAddress("Muon_isMed", Muon_isMed, &b_Muon_isMed);
   fChain->SetBranchAddress("Muon_isMedPr", Muon_isMedPr, &b_Muon_isMedPr);
   fChain->SetBranchAddress("Muon_isTight", Muon_isTight, &b_Muon_isTight);
   fChain->SetBranchAddress("Muon_isHighPt", Muon_isHighPt, &b_Muon_isHighPt);
   fChain->SetBranchAddress("Muon_isHighPttrk", Muon_isHighPttrk, &b_Muon_isHighPttrk);
   fChain->SetBranchAddress("Muon_MediumID", Muon_MediumID, &b_Muon_MediumID);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_p", Muon_p, &b_Muon_p);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_minchiso", Muon_minchiso, &b_Muon_minchiso);
   fChain->SetBranchAddress("Muon_minnhiso", Muon_minnhiso, &b_Muon_minnhiso);
   fChain->SetBranchAddress("Muon_minphiso", Muon_minphiso, &b_Muon_minphiso);
   fChain->SetBranchAddress("Muon_minisoall", Muon_minisoall, &b_Muon_minisoall);
   fChain->SetBranchAddress("Muon_dxy", Muon_dxy, &b_Muon_dxy);
   fChain->SetBranchAddress("Muon_dz", Muon_dz, &b_Muon_dz);
   fChain->SetBranchAddress("Muon_dxyErr", Muon_dxyErr, &b_Muon_dxyErr);
   fChain->SetBranchAddress("Muon_ip3d", Muon_ip3d, &b_Muon_ip3d);
   fChain->SetBranchAddress("Muon_ptErr", Muon_ptErr, &b_Muon_ptErr);
   fChain->SetBranchAddress("Muon_chi", Muon_chi, &b_Muon_chi);
   fChain->SetBranchAddress("Muon_ndf", Muon_ndf, &b_Muon_ndf);
   fChain->SetBranchAddress("Muon_ecal", Muon_ecal, &b_Muon_ecal);
   fChain->SetBranchAddress("Muon_hcal", Muon_hcal, &b_Muon_hcal);
   fChain->SetBranchAddress("Muon_pfiso", Muon_pfiso, &b_Muon_pfiso);
   fChain->SetBranchAddress("Muon_posmatch", Muon_posmatch, &b_Muon_posmatch);
   fChain->SetBranchAddress("Muon_trkink", Muon_trkink, &b_Muon_trkink);
   fChain->SetBranchAddress("Muon_segcom", Muon_segcom, &b_Muon_segcom);
   fChain->SetBranchAddress("Muon_hit", Muon_hit, &b_Muon_hit);
   fChain->SetBranchAddress("Muon_pixhit", Muon_pixhit, &b_Muon_pixhit);
   fChain->SetBranchAddress("Muon_mst", Muon_mst, &b_Muon_mst);
   fChain->SetBranchAddress("Muon_trklay", Muon_trklay, &b_Muon_trklay);
   fChain->SetBranchAddress("Muon_valfrac", Muon_valfrac, &b_Muon_valfrac);
   fChain->SetBranchAddress("Muon_dxy_sv", Muon_dxy_sv, &b_Muon_dxy_sv);
   fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_p", Electron_p, &b_Electron_p);
   fChain->SetBranchAddress("Electron_e", Electron_e, &b_Electron_e);
   fChain->SetBranchAddress("Electron_e_ECAL", Electron_e_ECAL, &b_Electron_e_ECAL);
   fChain->SetBranchAddress("Electron_mvaid_Fallv2WP90", Electron_mvaid_Fallv2WP90, &b_Electron_mvaid_Fallv2WP90);
   fChain->SetBranchAddress("Electron_mvaid_Fallv2WP90_noIso", Electron_mvaid_Fallv2WP90_noIso, &b_Electron_mvaid_Fallv2WP90_noIso);
   fChain->SetBranchAddress("Electron_mvaid_Fallv2WP80", Electron_mvaid_Fallv2WP80, &b_Electron_mvaid_Fallv2WP80);
   fChain->SetBranchAddress("Electron_mvaid_Fallv2WP80_noIso", Electron_mvaid_Fallv2WP80_noIso, &b_Electron_mvaid_Fallv2WP80_noIso);
   fChain->SetBranchAddress("Electron_dxy", Electron_dxy, &b_Electron_dxy);
   fChain->SetBranchAddress("Electron_dxyErr", Electron_dxyErr, &b_Electron_dxyErr);
   fChain->SetBranchAddress("Electron_dz", Electron_dz, &b_Electron_dz);
   fChain->SetBranchAddress("Electron_dzErr", Electron_dzErr, &b_Electron_dzErr);
   fChain->SetBranchAddress("Electron_ip3d", Electron_ip3d, &b_Electron_ip3d);
   fChain->SetBranchAddress("Electron_dxy_sv", Electron_dxy_sv, &b_Electron_dxy_sv);
   fChain->SetBranchAddress("Electron_hovere", Electron_hovere, &b_Electron_hovere);
   fChain->SetBranchAddress("Electron_chi", Electron_chi, &b_Electron_chi);
   fChain->SetBranchAddress("Electron_ndf", Electron_ndf, &b_Electron_ndf);
   fChain->SetBranchAddress("Electron_eoverp", Electron_eoverp, &b_Electron_eoverp);
   fChain->SetBranchAddress("Electron_ietaieta", Electron_ietaieta, &b_Electron_ietaieta);
   fChain->SetBranchAddress("Electron_misshits", Electron_misshits, &b_Electron_misshits);
   fChain->SetBranchAddress("Electron_pfiso_drcor", Electron_pfiso_drcor, &b_Electron_pfiso_drcor);
   fChain->SetBranchAddress("Electron_pfiso_eacor", Electron_pfiso_eacor, &b_Electron_pfiso_eacor);
   fChain->SetBranchAddress("Electron_pfiso04_eacor", Electron_pfiso04_eacor, &b_Electron_pfiso04_eacor);
   fChain->SetBranchAddress("Electron_eccalTrkEnergyPostCorr", Electron_eccalTrkEnergyPostCorr, &b_Electron_eccalTrkEnergyPostCorr);
   fChain->SetBranchAddress("Electron_energyScaleValue", Electron_energyScaleValue, &b_Electron_energyScaleValue);
   fChain->SetBranchAddress("Electron_energyScaleUp", Electron_energyScaleUp, &b_Electron_energyScaleUp);
   fChain->SetBranchAddress("Electron_energyScaleDown", Electron_energyScaleDown, &b_Electron_energyScaleDown);
   fChain->SetBranchAddress("Electron_energySigmaValue", Electron_energySigmaValue, &b_Electron_energySigmaValue);
   fChain->SetBranchAddress("Electron_energySigmaUp", Electron_energySigmaUp, &b_Electron_energySigmaUp);
   fChain->SetBranchAddress("Electron_energySigmaDown", Electron_energySigmaDown, &b_Electron_energySigmaDown);
   fChain->SetBranchAddress("Electron_supcl_eta", Electron_supcl_eta, &b_Electron_supcl_eta);
   fChain->SetBranchAddress("Electron_supcl_phi", Electron_supcl_phi, &b_Electron_supcl_phi);
   fChain->SetBranchAddress("Electron_supcl_e", Electron_supcl_e, &b_Electron_supcl_e);
   fChain->SetBranchAddress("Electron_supcl_rawE", Electron_supcl_rawE, &b_Electron_supcl_rawE);
   fChain->SetBranchAddress("Electron_sigmaieta", Electron_sigmaieta, &b_Electron_sigmaieta);
   fChain->SetBranchAddress("Electron_sigmaiphi", Electron_sigmaiphi, &b_Electron_sigmaiphi);
   fChain->SetBranchAddress("Electron_r9full", Electron_r9full, &b_Electron_r9full);
   fChain->SetBranchAddress("Electron_hcaloverecal", Electron_hcaloverecal, &b_Electron_hcaloverecal);
   fChain->SetBranchAddress("Electron_hitsmiss", Electron_hitsmiss, &b_Electron_hitsmiss);
   fChain->SetBranchAddress("Electron_ecloverpout", Electron_ecloverpout, &b_Electron_ecloverpout);
   fChain->SetBranchAddress("Electron_convVeto", Electron_convVeto, &b_Electron_convVeto);
   fChain->SetBranchAddress("Electron_pfisolsumphet", Electron_pfisolsumphet, &b_Electron_pfisolsumphet);
   fChain->SetBranchAddress("Electron_pfisolsumchhadpt", Electron_pfisolsumchhadpt, &b_Electron_pfisolsumchhadpt);
   fChain->SetBranchAddress("Electron_pfsiolsumneuhadet", Electron_pfsiolsumneuhadet, &b_Electron_pfsiolsumneuhadet);
   fChain->SetBranchAddress("Electron_minchiso", Electron_minchiso, &b_Electron_minchiso);
   fChain->SetBranchAddress("Electron_minnhiso", Electron_minnhiso, &b_Electron_minnhiso);
   fChain->SetBranchAddress("Electron_minphiso", Electron_minphiso, &b_Electron_minphiso);
   fChain->SetBranchAddress("Electron_minisoall", Electron_minisoall, &b_Electron_minisoall);
   fChain->SetBranchAddress("nPhoton", &nPhoton, &b_nPhoton);
   fChain->SetBranchAddress("Photon_e", Photon_e, &b_Photon_e);
   fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_mvaid_Fall17V2_raw", Photon_mvaid_Fall17V2_raw, &b_Photon_mvaid_Fall17V2_raw);
   fChain->SetBranchAddress("Photon_mvaid_Fall17V2_WP90", Photon_mvaid_Fall17V2_WP90, &b_Photon_mvaid_Fall17V2_WP90);
   fChain->SetBranchAddress("Photon_mvaid_Fall17V2_WP80", Photon_mvaid_Fall17V2_WP80, &b_Photon_mvaid_Fall17V2_WP80);
   fChain->SetBranchAddress("Photon_mvaid_Spring16V1_WP90", Photon_mvaid_Spring16V1_WP90, &b_Photon_mvaid_Spring16V1_WP90);
   fChain->SetBranchAddress("Photon_mvaid_Spring16V1_WP80", Photon_mvaid_Spring16V1_WP80, &b_Photon_mvaid_Spring16V1_WP80);
   fChain->SetBranchAddress("Photon_e1by9", Photon_e1by9, &b_Photon_e1by9);
   fChain->SetBranchAddress("Photon_e9by25", Photon_e9by25, &b_Photon_e9by25);
   fChain->SetBranchAddress("Photon_trkiso", Photon_trkiso, &b_Photon_trkiso);
   fChain->SetBranchAddress("Photon_emiso", Photon_emiso, &b_Photon_emiso);
   fChain->SetBranchAddress("Photon_hadiso", Photon_hadiso, &b_Photon_hadiso);
   fChain->SetBranchAddress("Photon_chhadiso", Photon_chhadiso, &b_Photon_chhadiso);
   fChain->SetBranchAddress("Photon_neuhadiso", Photon_neuhadiso, &b_Photon_neuhadiso);
   fChain->SetBranchAddress("Photon_phoiso", Photon_phoiso, &b_Photon_phoiso);
   fChain->SetBranchAddress("Photon_PUiso", Photon_PUiso, &b_Photon_PUiso);
   fChain->SetBranchAddress("Photon_hadbyem", Photon_hadbyem, &b_Photon_hadbyem);
   fChain->SetBranchAddress("Photon_ietaieta", Photon_ietaieta, &b_Photon_ietaieta);
   fChain->SetBranchAddress("nTau", &nTau, &b_nTau);
   fChain->SetBranchAddress("Tau_isPF", Tau_isPF, &b_Tau_isPF);
   fChain->SetBranchAddress("Tau_pt", Tau_pt, &b_Tau_pt);
   fChain->SetBranchAddress("Tau_eta", Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau_phi", Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau_e", Tau_e, &b_Tau_e);
   fChain->SetBranchAddress("Tau_charge", Tau_charge, &b_Tau_charge);
   fChain->SetBranchAddress("Tau_dxy", Tau_dxy, &b_Tau_dxy);
   fChain->SetBranchAddress("Tau_leadtrkdxy", Tau_leadtrkdxy, &b_Tau_leadtrkdxy);
   fChain->SetBranchAddress("Tau_leadtrkdz", Tau_leadtrkdz, &b_Tau_leadtrkdz);
   fChain->SetBranchAddress("Tau_leadtrkpt", Tau_leadtrkpt, &b_Tau_leadtrkpt);
   fChain->SetBranchAddress("Tau_leadtrketa", Tau_leadtrketa, &b_Tau_leadtrketa);
   fChain->SetBranchAddress("Tau_leadtrkphi", Tau_leadtrkphi, &b_Tau_leadtrkphi);
   fChain->SetBranchAddress("Tau_decayMode", Tau_decayMode, &b_Tau_decayMode);
   fChain->SetBranchAddress("Tau_decayModeinding", Tau_decayModeinding, &b_Tau_decayModeinding);
   fChain->SetBranchAddress("Tau_decayModeindingNewDMs", Tau_decayModeindingNewDMs, &b_Tau_decayModeindingNewDMs);
   fChain->SetBranchAddress("Tau_eiso2018_raw", Tau_eiso2018_raw, &b_Tau_eiso2018_raw);
   fChain->SetBranchAddress("Tau_eiso2018", Tau_eiso2018, &b_Tau_eiso2018);
   fChain->SetBranchAddress("Tau_jetiso_deeptau2017v2p1_raw", Tau_jetiso_deeptau2017v2p1_raw, &b_Tau_jetiso_deeptau2017v2p1_raw);
   fChain->SetBranchAddress("Tau_jetiso_deeptau2017v2p1", Tau_jetiso_deeptau2017v2p1, &b_Tau_jetiso_deeptau2017v2p1);
   fChain->SetBranchAddress("Tau_eiso_deeptau2017v2p1_raw", Tau_eiso_deeptau2017v2p1_raw, &b_Tau_eiso_deeptau2017v2p1_raw);
   fChain->SetBranchAddress("Tau_eiso_deeptau2017v2p1", Tau_eiso_deeptau2017v2p1, &b_Tau_eiso_deeptau2017v2p1);
   fChain->SetBranchAddress("Tau_muiso_deeptau2017v2p1_raw", Tau_muiso_deeptau2017v2p1_raw, &b_Tau_muiso_deeptau2017v2p1_raw);
   fChain->SetBranchAddress("Tau_muiso_deeptau2017v2p1", Tau_muiso_deeptau2017v2p1, &b_Tau_muiso_deeptau2017v2p1);
   fChain->SetBranchAddress("Tau_rawiso", Tau_rawiso, &b_Tau_rawiso);
   fChain->SetBranchAddress("Tau_rawisodR03", Tau_rawisodR03, &b_Tau_rawisodR03);
   fChain->SetBranchAddress("Tau_puCorr", Tau_puCorr, &b_Tau_puCorr);
   fChain->SetBranchAddress("event_weight", &event_weight, &b_event_weight);
   fChain->SetBranchAddress("qscale", &qscale, &b_qscale);
   fChain->SetBranchAddress("Generator_x1", &Generator_x1, &b_Generator_x1);
   fChain->SetBranchAddress("Generator_x2", &Generator_x2, &b_Generator_x2);
   fChain->SetBranchAddress("Generator_xpdf1", &Generator_xpdf1, &b_Generator_xpdf1);
   fChain->SetBranchAddress("Generator_xpdf2", &Generator_xpdf2, &b_Generator_xpdf2);
   fChain->SetBranchAddress("Generator_id1", &Generator_id1, &b_Generator_id1);
   fChain->SetBranchAddress("Generator_id2", &Generator_id2, &b_Generator_id2);
   fChain->SetBranchAddress("Generator_scalePDF", &Generator_scalePDF, &b_Generator_scalePDF);
   fChain->SetBranchAddress("npu_vert", &npu_vert, &b_npu_vert);
   fChain->SetBranchAddress("npu_vert_true", &npu_vert_true, &b_npu_vert_true);
   fChain->SetBranchAddress("GENMET_pt", &GENMET_pt, &b_genmiset);
   fChain->SetBranchAddress("GENMET_phi", &GENMET_phi, &b_genmisphi);
   fChain->SetBranchAddress("nGenJetAK8", &nGenJetAK8, &b_nGenJetAK8);
   fChain->SetBranchAddress("GenJetAK8_pt", GenJetAK8_pt, &b_GenJetAK8_pt);
   fChain->SetBranchAddress("GenJetAK8_eta", GenJetAK8_eta, &b_GenJetAK8_eta);
   fChain->SetBranchAddress("GenJetAK8_phi", GenJetAK8_phi, &b_GenJetAK8_phi);
   fChain->SetBranchAddress("GenJetAK8_mass", GenJetAK8_mass, &b_GenJetAK8_mass);
   fChain->SetBranchAddress("GenJetAK8_sdmass", GenJetAK8_sdmass, &b_GenJetAK8_sdmass);
   fChain->SetBranchAddress("GenJetAK8_hadronflav", GenJetAK8_hadronflav, &b_GenJetAK8_hadronflav);
   fChain->SetBranchAddress("GenJetAK8_partonflav", GenJetAK8_partonflav, &b_GenJetAK8_partonflav);
   fChain->SetBranchAddress("nGenJetAK4", &nGenJetAK4, &b_nGenJetAK4);
   fChain->SetBranchAddress("GenJetAK4_pt", GenJetAK4_pt, &b_GenJetAK4_pt);
   fChain->SetBranchAddress("GenJetAK4_eta", GenJetAK4_eta, &b_GenJetAK4_eta);
   fChain->SetBranchAddress("GenJetAK4_phi", GenJetAK4_phi, &b_GenJetAK4_phi);
   fChain->SetBranchAddress("GenJetAK4_mass", GenJetAK4_mass, &b_GenJetAK4_mass);
   fChain->SetBranchAddress("GenJetAK4_hadronflav", GenJetAK4_hadronflav, &b_GenJetAK4_hadronflav);
   fChain->SetBranchAddress("GenJetAK4_partonflav", GenJetAK4_partonflav, &b_GenJetAK4_partonflav);
   fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   fChain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
   fChain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
   fChain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
   fChain->SetBranchAddress("GenPart_m", GenPart_m, &b_GenPart_m);
   fChain->SetBranchAddress("GenPart_status", GenPart_status, &b_GenPart_status);
   fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
   fChain->SetBranchAddress("GenPart_mompdgId", GenPart_mompdgId, &b_GenPart_mompdgId);
   fChain->SetBranchAddress("GenPart_grmompdgId", GenPart_grmompdgId, &b_GenPart_grmompdgId);
   fChain->SetBranchAddress("GenPart_daugno", GenPart_daugno, &b_GenPart_daugno);
   fChain->SetBranchAddress("GenPart_fromhard", GenPart_fromhard, &b_GenPart_fromhard);
   fChain->SetBranchAddress("GenPart_fromhardbFSR", GenPart_fromhardbFSR, &b_GenPart_fromhardbFSR);
   fChain->SetBranchAddress("GenPart_isPromptFinalState", GenPart_isPromptFinalState, &b_GenPart_isPromptFinalState);
   fChain->SetBranchAddress("GenPart_isLastCopyBeforeFSR", GenPart_isLastCopyBeforeFSR, &b_GenPart_isLastCopyBeforeFSR);
   fChain->SetBranchAddress("nLHEPart", &nLHEPart, &b_nLHEPart);
   fChain->SetBranchAddress("LHEPart_pdg", LHEPart_pdg, &b_LHEPart_pdg);
   fChain->SetBranchAddress("LHEPart_pt", LHEPart_pt, &b_LHEPart_pt);
   fChain->SetBranchAddress("LHEPart_eta", LHEPart_eta, &b_LHEPart_eta);
   fChain->SetBranchAddress("LHEPart_phi", LHEPart_phi, &b_LHEPart_phi);
   fChain->SetBranchAddress("LHEPart_m", LHEPart_m, &b_LHEPart_m);
   fChain->SetBranchAddress("event_weight_LHE", &event_weight_LHE, &b_event_weight_LHE);
   fChain->SetBranchAddress("nLHEScaleWeights", &nLHEScaleWeights, &b_nLHEScaleWeights);
   fChain->SetBranchAddress("LHEScaleWeights", LHEScaleWeights, &b_LHEScaleWeights);
   fChain->SetBranchAddress("nLHEPDFWeights", &nLHEPDFWeights, &b_nLHEPDFWeights);
   fChain->SetBranchAddress("LHEPDFWeights", LHEPDFWeights, &b_LHEPDFWeights);
   fChain->SetBranchAddress("nLHEAlpsWeights", &nLHEAlpsWeights, &b_nLHEAlpsWeights);
   fChain->SetBranchAddress("LHEAlpsWeights", LHEAlpsWeights, &b_LHEAlpsWeights);
   fChain->SetBranchAddress("nLHEPSWeights", &nLHEPSWeights, &b_nLHEPSWeights);
   fChain->SetBranchAddress("LHEPSWeights", LHEPSWeights, &b_LHEPSWeights);
   
   gRandom = new TRandom3();
   
   int nentries = fChain->GetEntries();
   cout<<"nentries "<<nentries<<endl;
   
   isMC = true;
   isFastSIM = false;
   
   for (int ij=0; ij<nentries; ij++) {
   
	//if(ij%100==0){ cout<<"event "<<ij+1<<endl; }
	fileIn->cd();

	fChain->GetEntry(ij);

	//cout<<"Before: "<<nPFJetAK8<<" "<<nPFJetAK4<<" "<<nMuon<<" "<<nElectron<<" "<<nPhoton<<endl;
	
	//Here you get gen particles
	vector<GenParton> genpartons;
	getPartons(genpartons);
	
	vector<GenParton> lheparts;
	getLHEParts(lheparts);
	
	vector <AK8GenJet> genLJets;
    getAK8Genjets(genLJets,AK8GenJet_pt_cut,absetacut,isMC);

	
	weight = 1;
	if(!isFastSIM){
		weight *= event_weight;
	}
	
	// LHE story //
	
	//cout<<"lheparts.size() "<<lheparts.size()<<endl;
	unsigned LHE_b1, LHE_b2, LHE_j1, LHE_j2, LHE_l, LHE_nu;
	for(unsigned ilhe=0; ilhe<lheparts.size(); ++ilhe){
		//	cout<<ilhe+1<<" ID "<<lheparts[ilhe].pdgId<<" pt "<<lheparts[ilhe].p4.Pt()<<endl;	
			if(lheparts[ilhe].pdgId==5) { LHE_b1 = ilhe; }
			if(lheparts[ilhe].pdgId==-5) { LHE_b2 = ilhe; }
			if(abs(lheparts[ilhe].pdgId)==1||abs(lheparts[ilhe].pdgId)==3) { LHE_j1 = ilhe; }
			if(abs(lheparts[ilhe].pdgId)==2||abs(lheparts[ilhe].pdgId)==4) { LHE_j2 = ilhe; }
			if(abs(lheparts[ilhe].pdgId)==11||abs(lheparts[ilhe].pdgId)==13/*||abs(lheparts[ilhe].pdgId)==15*/) { LHE_l = ilhe; }
			if(abs(lheparts[ilhe].pdgId)==12||abs(lheparts[ilhe].pdgId)==14/*||abs(lheparts[ilhe].pdgId)==16*/) { LHE_nu = ilhe; }
		}
	if(LHE_b1>=0 && LHE_b2>=0 && LHE_j1>=0 && LHE_j2>=0 && LHE_l>=0 && LHE_nu>=0){
		
		TLorentzVector p4_b1, p4_b2, p4_j1, p4_j2, p4_l, p4_nu;
		TLorentzVector p4_Wh, p4_Wl, p4_H, p4_Y;
		
		p4_b1 = lheparts[LHE_b1].p4;
		p4_b2 = lheparts[LHE_b2].p4;
		p4_j1 = lheparts[LHE_j1].p4;
		p4_j2 = lheparts[LHE_j2].p4;
		p4_l = lheparts[LHE_l].p4;
		p4_nu = lheparts[LHE_nu].p4;
		
		h_LHE_b1b2_pt->Fill((p4_b1+p4_b2).Pt());
		h_LHE_j1j2_pt->Fill((p4_j1+p4_j2).Pt());
		h_LHE_l_pt->Fill((p4_l).Pt());
		h_LHE_nu_pt->Fill((p4_nu).Pt());
		h_LHE_H_pt->Fill((p4_nu+p4_l+p4_j1+p4_j2).Pt());
		
		h_LHE_b1b2_m->Fill((p4_b1+p4_b2).M());
		h_LHE_j1j2_m->Fill((p4_j1+p4_j2).M());
		h_LHE_H_m->Fill((p4_j1+p4_j2+p4_l+p4_nu).M());
		h_LHE_X_m->Fill((p4_b1+p4_b2+p4_j1+p4_j2+p4_l+p4_nu).M());
		
		h_LHE_b1b2_dR->Fill(delta2R_vec(p4_b1,p4_b2));
		h_LHE_j1j2_dR->Fill(delta2R_vec(p4_j1,p4_j2));
		h_LHE_j1j2l_dR->Fill(delta2R_vec(p4_j1+p4_j2,p4_l));
		h_LHE_j1j2l_b1b2_dR->Fill(delta2R_vec(p4_b1+p4_b2,p4_j1+p4_j2+p4_l));
		
		}
	
	// GEN story //
	
		unsigned GEN_b1, GEN_b2, GEN_Y, GEN_W, GEN_j1, GEN_j2, GEN_l, GEN_nu;
		GEN_Y = GEN_b1 = GEN_b2 = GEN_W = GEN_j1 = GEN_j2 = GEN_l = GEN_nu = 1000000;
		
		for(unsigned igen=0; igen<genpartons.size(); ++igen){
			if((genpartons[igen].status==23 || genpartons[igen].status==1)&&(genpartons[igen].fromhard)){
				if(genpartons[igen].pdgId==5 &&  abs(genpartons[igen].mompdgId)==35 && abs(genpartons[igen].grmompdgId)==45) { GEN_b1 = igen;  }
				if(genpartons[igen].pdgId==-5 &&  abs(genpartons[igen].mompdgId)==35 && abs(genpartons[igen].grmompdgId)==45) { GEN_b2 = igen;  }
				if((abs(genpartons[igen].pdgId)==1 || abs(genpartons[igen].pdgId)==3) &&  ((abs(genpartons[igen].mompdgId)==24 && abs(genpartons[igen].grmompdgId)==25)||(abs(genpartons[igen].mompdgId)==25 && abs(genpartons[igen].grmompdgId)==45))) { GEN_j1 = igen;  }
				if((abs(genpartons[igen].pdgId)==2 || abs(genpartons[igen].pdgId)==4) &&  ((abs(genpartons[igen].mompdgId)==24 && abs(genpartons[igen].grmompdgId)==25)||(abs(genpartons[igen].mompdgId)==25 && abs(genpartons[igen].grmompdgId)==45))) { GEN_j2 = igen;  }
				if((abs(genpartons[igen].pdgId)==11 || abs(genpartons[igen].pdgId)==13) &&  ((abs(genpartons[igen].mompdgId)==24 && abs(genpartons[igen].grmompdgId)==25)||(abs(genpartons[igen].mompdgId)==25 && abs(genpartons[igen].grmompdgId)==45))) { GEN_l = igen;  }
				if((abs(genpartons[igen].pdgId)==12 || abs(genpartons[igen].pdgId)==14) &&  ((abs(genpartons[igen].mompdgId)==24 && abs(genpartons[igen].grmompdgId)==25)||(abs(genpartons[igen].mompdgId)==25 && abs(genpartons[igen].grmompdgId)==45))) { GEN_nu = igen;  }
		
		////		cout<<"GEN:"<<igen+1<<" "<<genpartons[igen].pdgId<<" "<<genpartons[igen].mompdgId<<" "<<genpartons[igen].grmompdgId<<endl; 
			}
		}
		
		//cout<<GEN_b1<<" "<<GEN_b2<<" "<<GEN_j1<<" "<<GEN_j2<<" "<<GEN_l<<" "<<GEN_nu<<endl;
		
		if(LHE_b1>=0 && LHE_b1<(lheparts.size()) && LHE_b2>=0 && LHE_b2<(lheparts.size())){
			
			float dRmin = 0.8;
			for(unsigned ig=0; ig<(genLJets.size()); ig++){
				if(delta2R_vec(genLJets[ig].p4,(lheparts[LHE_b1].p4+lheparts[LHE_b2].p4))<dRmin){
					dRmin = delta2R_vec(genLJets[ig].p4,(lheparts[LHE_b1].p4+lheparts[LHE_b2].p4));
					GEN_Y = ig;
					}
				}
			}
			
		if(LHE_j1>=0 && LHE_j1<(lheparts.size()) && LHE_j2>=0 && LHE_j2<(lheparts.size())){
			
			float dRmin = 0.8;
			for(unsigned ig=0; ig<(genLJets.size()); ig++){
				if(delta2R_vec(genLJets[ig].p4,(lheparts[LHE_j1].p4+lheparts[LHE_j2].p4))<dRmin){
					dRmin = delta2R_vec(genLJets[ig].p4,(lheparts[LHE_j1].p4+lheparts[LHE_j2].p4));
					GEN_W = ig;
					}
				}
			}
		
		if(GEN_b1>=0 && GEN_b1<(genpartons.size()) && GEN_b2>=0 && GEN_b2<(genpartons.size()) && GEN_j1>=0 && GEN_j1<(genpartons.size()) && GEN_j2>=0 && GEN_j2<(genpartons.size()) && 
		GEN_l>=0 && GEN_l<(genpartons.size()) && GEN_nu>=0 && GEN_nu<(genpartons.size())){
		
			TLorentzVector p4_b1, p4_b2, p4_j1, p4_j2, p4_l, p4_nu;
			
			p4_b1 = genpartons[GEN_b1].p4;
			p4_b2 = genpartons[GEN_b2].p4;
			p4_j1 = genpartons[GEN_j1].p4;
			p4_j2 = genpartons[GEN_j2].p4;
			p4_l = genpartons[GEN_l].p4;
			p4_nu = genpartons[GEN_nu].p4;
		
			h_GEN_b1b2_pt->Fill((p4_b1+p4_b2).Pt());
			h_GEN_j1j2_pt->Fill((p4_j1+p4_j2).Pt());
			h_GEN_l_pt->Fill((p4_l).Pt());
			h_GEN_nu_pt->Fill((p4_nu).Pt());
			h_GEN_H_pt->Fill((p4_nu+p4_l+p4_j1+p4_j2).Pt());
		
			h_GEN_b1b2_m->Fill((p4_b1+p4_b2).M());
			h_GEN_j1j2_m->Fill((p4_j1+p4_j2).M());
			h_GEN_H_m->Fill((p4_j1+p4_j2+p4_l+p4_nu).M());
			h_GEN_X_m->Fill((p4_b1+p4_b2+p4_j1+p4_j2+p4_l+p4_nu).M());
		
			h_GEN_b1b2_dR->Fill(delta2R_vec(p4_b1,p4_b2));
			h_GEN_j1j2_dR->Fill(delta2R_vec(p4_j1,p4_j2));
			h_GEN_j1j2l_dR->Fill(delta2R_vec(p4_j1+p4_j2,p4_l));
			h_GEN_j1j2l_b1b2_dR->Fill(delta2R_vec(p4_b1+p4_b2,p4_j1+p4_j2+p4_l));
			
			
		}
		
		if(GEN_l>=0 && GEN_l<(genpartons.size()) && GEN_nu>=0 && GEN_nu<(genpartons.size())){
		
			TLorentzVector p4_Wh, p4_Wl, p4_H, p4_Y, p4_GMet, p4_l, p4_nu;
			p4_l = genpartons[GEN_l].p4;
			p4_nu = genpartons[GEN_nu].p4;
			
			if(GEN_W>=0 && GEN_W<(genLJets.size()) && GEN_Y>=0 && GEN_Y<(genLJets.size()) && GEN_Y!=GEN_W){
				
				p4_Wh = genLJets[GEN_W].p4;
				p4_Y = genLJets[GEN_Y].p4;
				double random_no = gRandom->Uniform(0,1);
				p4_GMet = neutrino_mom_fromH((p4_l+p4_Wh), GENMET_pt, GENMET_phi, random_no);
						
				h_GEN_Y_m->Fill((p4_Y).M());
				h_GEN_Wh_m->Fill((p4_Wh).M());
				h_GEN_HJ_m->Fill((p4_Wh+p4_l+p4_nu).M());
				h_GEN_XJ_m->Fill((p4_Y+p4_Wh+p4_l+p4_nu).M());
				h_GEN_HJ_m_MET->Fill((p4_Wh+p4_l+p4_GMet).M());
				h_GEN_XJ_m_MET->Fill((p4_Y+p4_Wh+p4_l+p4_GMet).M());
				
				h_GEN_Y_XJ_m->Fill((p4_Y+p4_Wh+p4_l+p4_nu).M(),(p4_Y).M());
				h_GEN_Y_XJ_m_MET->Fill((p4_Y+p4_Wh+p4_l+p4_GMet).M(),(p4_Y).M());
			}
		}
		
	
	}// entry
     
     fileIn->cd();
     delete fChain;
   
      }
      
    file_db.close();    
    fileout->cd();
	fileout->Write();
	fileout->Close();
}

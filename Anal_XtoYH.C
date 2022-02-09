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


int main(int argc, char *argv[])
{
  cout<<"Program started"<<endl;
  char fOut[50], fout1[50],fout2[50];
  string inputFile=argv[3];
  if(inputFile=="FILELIST_2018UL/XtoYH_2018_TT_Had.log"){
      sprintf(fOut,"HIST_2018UL/JetHT_EraB_2017_UL_HFNoise_%s_%s.root",argv[1],argv[2]);
   }
  else if(inputFile=="FILELIST_2018UL/TTToSemiLeptonic_XtoYH_Nov_2021.log"){
     sprintf(fOut,"HIST_2018UL/TTToSemiLeptonic_XtoYH_Nov_2021_%s_%s.root",argv[1],argv[2]);
   }
  else if(inputFile=="FILELIST_2018UL/ZZ_XtoYH_Nov_2021.log"){
     sprintf(fOut,"HIST_2018UL/ZZ_XtoYH_Nov_2021_%s_%s.root",argv[1],argv[2]);
   }
  else if(inputFile=="FILELIST_2018UL/WZ_XtoYH_Nov_2021.log"){
     sprintf(fOut,"HIST_2018UL/WZ_XtoYH_Nov_2021_%s_%s.root",argv[1],argv[2]);
   }

  else if(inputFile=="FILELIST_2018UL/WW_XtoYH_Nov_2021.log"){
     sprintf(fOut,"HIST_2018UL/WW_XtoYH_Nov_2021_%s_%s.root",argv[1],argv[2]);
   }

  else if(inputFile=="FILELIST_2018UL/TTToHadronic_XtoYH_Nov_2021.log"){
     sprintf(fOut,"HIST_2018UL/TTToHadronic_XtoYH_Nov_2021_%s_%s.root",argv[1],argv[2]);
   }

  else if(inputFile=="FILELIST_2018UL/TTTo2L2Nu_XtoYH_Nov_2021.log"){
     sprintf(fOut,"HIST_2018UL/TTTo2L2Nu_XtoYH_Nov_2021_%s_%s.root",argv[1],argv[2]);
   }

  else if(inputFile=="FILELIST_2018UL/ST_tW_top_XtoYH_Nov_2021.log"){
     sprintf(fOut,"HIST_2018UL/ST_tW_top_XtoYH_Nov_2021_%s_%s.root",argv[1],argv[2]);
   }

  else if(inputFile=="FILELIST_2018UL/ST_tW_antitop_XtoYH_Nov_2021.log"){
     sprintf(fOut,"HIST_2018UL/ST_tW_antitop_XtoYH_Nov_2021_%s_%s.root",argv[1],argv[2]);
   }

  else if(inputFile=="FILELIST_2018UL/ST_s-channel_XtoYH_Nov_2021.log"){
     sprintf(fOut,"HIST_2018UL/ST_s-channel_XtoYH_Nov_2021_%s_%s.root",argv[1],argv[2]);
   }

  else if(inputFile=="FILELIST_2018UL/DYJetsToLL_M-50_XtoYH_Nov_2021.log"){
     sprintf(fOut,"HIST_2018UL/DYJetsToLL_M-50_XtoYH_Nov_2021_%s_%s.root",argv[1],argv[2]);
   }

  else if(inputFile=="FILELIST_2018UL/DYJetsToLL_M-10to50_XtoYH_Nov_2021.log"){
     sprintf(fOut,"HIST_2018UL/DYJetsToLL_M-10to50_XtoYH_Nov_2021_%s_%s.root",argv[1],argv[2]);
   }
  else if(inputFile=="FILELIST_2018UL/WJetsToLNu_XtoYH_Nov_2021.log"){
     sprintf(fOut,"HIST_2018UL/WJetsToLNu_XtoYH_Nov_2021_%s_%s.root",argv[1],argv[2]);
   }
  else if(inputFile=="FILELIST_2018UL/NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_1500_MY_200_XtoYH_Nov_2021.log"){
     sprintf(fOut,"HIST_2018UL/NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_1500_MY_200_XtoYH_Nov_2021_%s_%s.root",argv[1],argv[2]);
   }
  else if(inputFile=="FILELIST_2018UL/NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2000_MY_200_XtoYH_Nov_2021.log"){
     sprintf(fOut,"HIST_2018UL/NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2000_MY_200_XtoYH_Nov_2021_%s_%s.root",argv[1],argv[2]);
   }

   TFile *fileout = new TFile(fOut,"recreate");
   
   Tout = new TTree("Tout", "Results");
   // Histogram definition for PU, NV 
   TH1F *h_PU = new TH1F("h_PU", "h_PU", 100, 0, 100);
   TH1F *h_TPU = new TH1F("h_TPU", "h_TPU", 100, 0, 100);
   TH1F *h_PV = new TH1F("h_PV", "h_PV", 100, 0, 100);

   //Binning definition
   Double_t xEdges0[] =  {200,300,400,800,7000};
   Double_t yEdges0[] =  {0,0.6,1.2,2.4}; 
   const int XBINS0 = sizeof(xEdges0)/sizeof(xEdges0[0])-1; 
   const int YBINS0 = sizeof(yEdges0)/sizeof(yEdges0[0])-1;
   // Histogram definition for AK4 
   TH2F* h_Ak4_b_flv = new TH2F("h_Ak4_b_flv", "h_Ak4_b_flv", XBINS0, xEdges0, YBINS0, yEdges0);
   TH2F* h_Ak4_c_flv = new TH2F("h_Ak4_c_flv", "h_Ak4_c_flv", XBINS0, xEdges0, YBINS0, yEdges0);
   TH2F* h_Ak4_l_flv = new TH2F("h_Ak4_l_flv", "h_Ak4_l_flv", XBINS0, xEdges0, YBINS0, yEdges0);

   // Histogram definition for AK4 passed L DeepFlv
   TH2F* h_Ak4_b_flv_pass_L = new TH2F("h_Ak4_b_flv_pass_L", "h_Ak4_b_flv_pass_L", XBINS0, xEdges0, YBINS0, yEdges0);
   TH2F* h_Ak4_c_flv_pass_L = new TH2F("h_Ak4_c_flv_pass_L", "h_Ak4_c_flv_pass_L", XBINS0, xEdges0, YBINS0, yEdges0);
   TH2F* h_Ak4_l_flv_pass_L = new TH2F("h_Ak4_l_flv_pass_L", "h_Ak4_l_flv_pass_L", XBINS0, xEdges0, YBINS0, yEdges0);
   
   // Histogram definition for AK4 passed M DeepFlv
   TH2F* h_Ak4_b_flv_pass_M = new TH2F("h_Ak4_b_flv_pass_M", "h_Ak4_b_flv_pass_M", XBINS0, xEdges0, YBINS0, yEdges0);
   TH2F* h_Ak4_c_flv_pass_M = new TH2F("h_Ak4_c_flv_pass_M", "h_Ak4_c_flv_pass_M", XBINS0, xEdges0, YBINS0, yEdges0);
   TH2F* h_Ak4_l_flv_pass_M = new TH2F("h_Ak4_l_flv_pass_M", "h_Ak4_l_flv_pass_M", XBINS0, xEdges0, YBINS0, yEdges0);

   // Histogram definition for AK4 passed T DeepFlv
   TH2F* h_Ak4_b_flv_pass_T = new TH2F("h_Ak4_b_flv_pass_T", "h_Ak4_b_flv_pass_T", XBINS0, xEdges0, YBINS0, yEdges0);
   TH2F* h_Ak4_c_flv_pass_T = new TH2F("h_Ak4_c_flv_pass_T", "h_Ak4_c_flv_pass_T", XBINS0, xEdges0, YBINS0, yEdges0);
   TH2F* h_Ak4_l_flv_pass_T = new TH2F("h_Ak4_l_flv_pass_T", "h_Ak4_l_flv_pass_T", XBINS0, xEdges0, YBINS0, yEdges0);

   //Binning definition   
   Double_t xEdges1[] =  {200,300,400,800,7000};
   Double_t yEdges1[] =  {0,0.6,1.2,2.4};
   const int XBINS1 = sizeof(xEdges1)/sizeof(xEdges1[0])-1;
   const int YBINS1 = sizeof(yEdges1)/sizeof(yEdges1[0])-1;
   
   TH2F* h_Ak8_DeepTag_PNetMD_WvsQCD = new TH2F("h_Ak8_DeepTag_PNetMD_WvsQCD", "h_Ak8_DeepTag_PNetMD_WvsQCD", XBINS1, xEdges1, YBINS1, yEdges1);
   TH2F* h_Ak8_DeepTag_PNetMD_WvsQCD_pass_L = new TH2F("h_Ak8_DeepTag_PNetMD_WvsQCD_pass_L", "h_Ak8_DeepTag_PNetMD_WvsQCD_pass_L", XBINS1, xEdges1, YBINS1, yEdges1);
   TH2F* h_Ak8_DeepTag_PNetMD_WvsQCD_pass_M = new TH2F("h_Ak8_DeepTag_PNetMD_WvsQCD_pass_M", "h_Ak8_DeepTag_PNetMD_WvsQCD_pass_M", XBINS1, xEdges1, YBINS1, yEdges1);
   TH2F* h_Ak8_DeepTag_PNetMD_WvsQCD_pass_T = new TH2F("h_Ak8_DeepTag_PNetMD_WvsQCD_pass_T", "h_Ak8_DeepTag_PNetMD_WvsQCD_pass_T", XBINS1, xEdges1, YBINS1, yEdges1);

   Double_t xEdges2[] =  {200,300,400,500,600,800,7000};
   Double_t yEdges2[] =  {0,0.6,1.2,2.4};
   const int XBINS2 = sizeof(xEdges2)/sizeof(xEdges2[0])-1;
   const int YBINS2 = sizeof(yEdges2)/sizeof(yEdges2[0])-1;

   TH2F* h_Ak8_DeepTag_PNetMD_XbbvsQCD = new TH2F("h_Ak8_DeepTag_PNetMD_XbbvsQCD", "h_Ak8_DeepTag_PNetMD_XbbvsQCD", XBINS2, xEdges2, YBINS2, yEdges2);
   TH2F* h_Ak8_DeepTag_PNetMD_XbbvsQCD_pass_L = new TH2F("h_Ak8_DeepTag_PNetMD_XbbvsQCD_pass_L", "h_Ak8_DeepTag_PNetMD_XbbvsQCD_pass_L", XBINS2, xEdges2, YBINS2, yEdges2);
   TH2F* h_Ak8_DeepTag_PNetMD_XbbvsQCD_pass_M = new TH2F("h_Ak8_DeepTag_PNetMD_XbbvsQCD_pass_M", "h_Ak8_DeepTag_PNetMD_XbbvsQCD_pass_M", XBINS2, xEdges2, YBINS2, yEdges2);
   TH2F* h_Ak8_DeepTag_PNetMD_XbbvsQCD_pass_T = new TH2F("h_Ak8_DeepTag_PNetMD_XbbvsQCD_pass_T", "h_Ak8_DeepTag_PNetMD_XbbvsQCD_pass_T", XBINS2, xEdges2, YBINS2, yEdges2);


   Tout->Branch("nleptons", &nleptons, "nleptons/I");
   Tout->Branch("nfatjets", &nfatjets, "nfatjets/I");	
   
   Tout->Branch("Flag_event_cuts", Flag_event_cuts, "Flag_event_cuts/O");	
   
   Tout->Branch("puWeight", &puWeight, "puWeight/F");
   Tout->Branch("puWeightup", &puWeightup, "puWeightup/F");	
   Tout->Branch("puWeightdown", &puWeightdown, "puWeightdown/F");	
   
   Tout->Branch("leptonsf_weight", &leptonsf_weight, "leptonsf_weight/F");
   Tout->Branch("leptonsf_weight_stat", &leptonsf_weight_stat, "leptonsf_weight_stat/F");	
   Tout->Branch("leptonsf_weight_syst", &leptonsf_weight_syst, "leptonsf_weight_syst/F");	
   
   Tout->Branch("l_pt", &l_pt, "l_pt/F");	
   Tout->Branch("l_eta", &l_eta, "l_eta/F");	
   Tout->Branch("l_phi", &l_phi, "l_phi/F");	
   Tout->Branch("l_mass", &l_mass, "l_mass/F");	
   Tout->Branch("l_genindex", &l_genindex, "l_genindex/I");	

   Tout->Branch("Y_pt", &Y_pt, "Y_pt/F");	
   Tout->Branch("Y_y", &Y_y, "Y_y/F");	
   Tout->Branch("Y_eta", &Y_eta, "Y_eta/F");
   Tout->Branch("Y_phi", &Y_phi, "Y_phi/F");	
   Tout->Branch("Y_mass", &Y_mass, "Y_mass/F");	
   Tout->Branch("Y_genindex", &Y_genindex, "Y_genindex/I");	
   Tout->Branch("Y_genbindex", Y_genbindex, "Y_genbindex[2]/I");	
   Tout->Branch("Y_sdmass", &Y_sdmass, "Y_sdmass/F");	
   Tout->Branch("Y_PN_bb", &Y_PN_bb, "Y_PN_bb/F");	
   Tout->Branch("Y_JESup", &Y_JESup, "Y_JESup/F");	
   Tout->Branch("Y_JESdn", &Y_JESdn, "Y_JESdn/F");	
   Tout->Branch("Y_JERup", &Y_JERup, "Y_JERup/F");	
   Tout->Branch("Y_JERdn", &Y_JERdn, "Y_JERdn/F");	
   
   // W boson related branches based on option -1 
   Tout->Branch("W_pt_opt1", &W_pt_opt1, "W_pt_opt1/F");	
   Tout->Branch("W_y_opt1", &W_y_opt1, "W_y_opt1/F");	
   Tout->Branch("W_eta_opt1", &W_eta_opt1, "W_eta_opt1/F");
   Tout->Branch("W_phi_opt1", &W_phi_opt1, "W_phi_opt1/F");	
   Tout->Branch("W_mass_opt1", &W_mass_opt1, "W_mass_opt1/F");	
   Tout->Branch("W_sdmass_opt1", &W_sdmass_opt1, "W_sdmass_opt1/F");	
   Tout->Branch("W_DAK8_W_opt1", &W_DAK8_W_opt1, "W_DAK8_W_opt1/F");
   Tout->Branch("W_PN_W_opt1", &W_PN_W_opt1, "W_PN_W_opt1/F");	
   Tout->Branch("W_genindex_opt1", &W_genindex_opt1, "W_genindex_opt1/I");	
   Tout->Branch("W_JESup_opt1", &W_JESup_opt1, "W_JESup_opt1/F");	
   Tout->Branch("W_JESdn_opt1", &W_JESdn_opt1, "W_JESdn_opt1/F");	
   Tout->Branch("W_JERup_opt1", &W_JERup_opt1, "W_JERup_opt1/F");	
   Tout->Branch("W_JERdn_opt1", &W_JERdn_opt1, "W_JERdn_opt1/F");	
   
   Tout->Branch("H_pt_opt1", &H_pt_opt1, "H_pt_opt1/F");	
   Tout->Branch("H_y_opt1", &H_y_opt1, "H_y_opt1/F");
   Tout->Branch("H_eta_opt1", &H_eta_opt1, "H_eta_opt1/F");
   Tout->Branch("H_phi_opt1", &H_phi_opt1, "H_phi_opt1/F");	
   Tout->Branch("H_mass_opt1", &H_mass_opt1, "H_mass_opt1/F");	
   Tout->Branch("H_genindex_opt1", &H_genindex_opt1, "H_genindex_opt1/I");	
   Tout->Branch("H_JESup_opt1", &H_JESup_opt1, "H_JESup_opt1/F");	
   Tout->Branch("H_JESdn_opt1", &H_JESdn_opt1, "H_JESdn_opt1/F");	
   Tout->Branch("H_JERup_opt1", &H_JERup_opt1, "H_JERup_opt1/F");	
   Tout->Branch("H_JERdn_opt1", &H_JERdn_opt1, "H_JERdn_opt1/F");	
   
   Tout->Branch("X_mass_opt1", &X_mass_opt1, "X_mass_opt1/F");	
   
   Tout->Branch("dR_lW_opt1", &dR_lW_opt1, "dR_lW_opt1/F");	
   Tout->Branch("dy_lW_opt1", &dy_lW_opt1, "dy_lW_opt1/F");	
   Tout->Branch("dphi_lW_opt1", &dphi_lW_opt1, "dphi_lW_opt1/F");	
   
   // W boson related branches based on option 2
   Tout->Branch("W_pt_opt2", &W_pt_opt2, "W_pt_opt2/F");
   Tout->Branch("W_y_opt2", &W_y_opt2, "W_y_opt2/F");
   Tout->Branch("W_eta_opt2", &W_eta_opt2, "W_eta_opt2/F");
   Tout->Branch("W_phi_opt2", &W_phi_opt2, "W_phi_opt2/F");
   Tout->Branch("W_mass_opt2", &W_mass_opt2, "W_mass_opt2/F");
   Tout->Branch("W_sdmass_opt2", &W_sdmass_opt2, "W_sdmass_opt2/F");
   Tout->Branch("W_DAK8_W_opt2", &W_DAK8_W_opt2, "W_DAK8_W_opt2/F");
   Tout->Branch("W_PN_W_opt2", &W_PN_W_opt2, "W_PN_W_opt2/F");
   Tout->Branch("W_genindex_opt2", &W_genindex_opt2, "W_genindex_opt2/I");	
   Tout->Branch("W_JESup_opt2", &W_JESup_opt2, "W_JESup_opt2/F");	
   Tout->Branch("W_JESdn_opt2", &W_JESdn_opt2, "W_JESdn_opt2/F");	
   Tout->Branch("W_JERup_opt2", &W_JERup_opt2, "W_JERup_opt2/F");	
   Tout->Branch("W_JERdn_opt2", &W_JERdn_opt2, "W_JERdn_opt2/F");	

   Tout->Branch("H_pt_opt2", &H_pt_opt2, "H_pt_opt2/F");
   Tout->Branch("H_y_opt2", &H_y_opt2, "H_y_opt2/F");
   Tout->Branch("H_eta_opt2", &H_eta_opt2, "H_eta_opt2/F");
   Tout->Branch("H_phi_opt2", &H_phi_opt2, "H_phi_opt2/F");
   Tout->Branch("H_mass_opt2", &H_mass_opt2, "H_mass_opt2/F");
   Tout->Branch("H_genindex_opt2", &H_genindex_opt2, "H_genindex_opt2/I");	
   Tout->Branch("H_JESup_opt2", &H_JESup_opt2, "H_JESup_opt2/F");	
   Tout->Branch("H_JESdn_opt2", &H_JESdn_opt2, "H_JESdn_opt2/F");	
   Tout->Branch("H_JERup_opt2", &H_JERup_opt2, "H_JERup_opt2/F");	
   Tout->Branch("H_JERdn_opt2", &H_JERdn_opt2, "H_JERdn_opt2/F");	

   Tout->Branch("X_mass_opt2", &X_mass_opt2, "X_mass_opt2/F");

   Tout->Branch("dR_lW_opt2", &dR_lW_opt2, "dR_lW_opt2/F");
   Tout->Branch("dy_lW_opt2", &dy_lW_opt2, "dy_lW_opt2/F");
   Tout->Branch("dphi_lW_opt2", &dphi_lW_opt2, "dphi_lW_opt2/F");

   Tout->Branch("dR_lY", &dR_lY, "dR_lY/F");	
   Tout->Branch("dy_lY", &dy_lY, "dy_lY/F");	
   Tout->Branch("dphi_lY", &dphi_lY, "dphi_lY/F");
   
   Tout->Branch("nbjets_other", &nbjets_other, "nbjets_other/I");	
   Tout->Branch("MET", &MET, "MET/F");	
   
   // Different flags and control regions   

   Tout->Branch("Flag_Y_bb_pass_T", &Flag_Y_bb_pass_T, "Flag_Y_bb_pass_T/O");	
   Tout->Branch("Flag_Y_bb_pass_M", &Flag_Y_bb_pass_M, "Flag_Y_bb_pass_M/O");
   Tout->Branch("Flag_Y_bb_pass_L", &Flag_Y_bb_pass_L, "Flag_Y_bb_pass_L/O");	
   Tout->Branch("Flag_H_W_pass_T_opt1", &Flag_H_W_pass_T_opt1, "Flag_H_W_pass_T_opt1/O");	
   Tout->Branch("Flag_H_W_pass_M_opt1", &Flag_H_W_pass_M_opt1, "Flag_H_W_pass_M_opt1/O");
   Tout->Branch("Flag_H_W_pass_L_opt1", &Flag_H_W_pass_L_opt1, "Flag_H_W_pass_L_opt1/O");

   Tout->Branch("Flag_H_W_pass_T_opt2", &Flag_H_W_pass_T_opt2, "Flag_H_W_pass_T_opt2/O");
   Tout->Branch("Flag_H_W_pass_M_opt2", &Flag_H_W_pass_M_opt2, "Flag_H_W_pass_M_opt2/O");
   Tout->Branch("Flag_H_W_pass_L_opt2", &Flag_H_W_pass_L_opt2, "Flag_H_W_pass_L_opt2/O");
   
   Tout->Branch("Flag_H_m_pass_opt1", &Flag_H_m_pass_opt1, "Flag_H_m_pass_opt1/O");	
   Tout->Branch("Flag_H_m_pass_opt2", &Flag_H_m_pass_opt2, "Flag_H_m_pass_opt2/O");
   
   Tout->Branch("Flag_dR_lW_pass_opt1", &Flag_dR_lW_pass_opt1, "Flag_dR_lW_pass_opt1/O");	
   Tout->Branch("Flag_dR_lW_pass_opt2", &Flag_dR_lW_pass_opt2, "Flag_dR_lW_pass_opt2/O");
  
   Tout->Branch("Flag_MET_pass", &Flag_MET_pass, "Flag_MET_pass/O");	
   Tout->Branch("Reg_SR_opt1", &Reg_SR_opt1, "Reg_SR_opt1/O");	
   Tout->Branch("Reg_Wj_CR_opt1", &Reg_Wj_CR_opt1, "Reg_Wj_CR_opt1/O");
   Tout->Branch("Reg_SR_opt2", &Reg_SR_opt2, "Reg_SR_opt2/O");
   Tout->Branch("Reg_Wj_CR_opt2", &Reg_Wj_CR_opt2, "Reg_Wj_CR_opt2/O");
   
   Tout->Branch("LHE_weight", &LHE_weight, "LHE_weight/D");	
   Tout->Branch("Generator_weight", &Generator_weight, "Generator_weight/D");	
   Tout->Branch("Event_weight", &weight, "weight/D");
   
   Tout->Branch("prefiringweight", &prefiringweight, "prefiringweight/D");	
   Tout->Branch("prefiringweightup", &prefiringweightup, "prefiringweightup/D");	
   Tout->Branch("prefiringweightdown", &prefiringweightdown, "prefiringweightdown/D");	
   
   // GEN particles //
   
   Tout->Branch("nGenLep",&nGenLep, "nGenLep/I");
   Tout->Branch("GenLep_pt",GenLep_pt,"GenLep_pt[nGenLep]/F");
   Tout->Branch("GenLep_eta",GenLep_eta,"GenLep_eta[nGenLep]/F");
   Tout->Branch("GenLep_phi",GenLep_phi,"GenLep_phi[nGenLep]/F");
   Tout->Branch("GenLep_mass",GenLep_mass,"GenLep_mass[nGenLep]/F");
   Tout->Branch("GenLep_pdgId",GenLep_pdgId,"GenLep_pdgId[nGenLep]/I");
   Tout->Branch("GenLep_mompdgId",GenLep_mompdgId,"GenLep_mompdgId[nGenLep]/I");
   Tout->Branch("GenLep_grmompdgId",GenLep_grmompdgId,"GenLep_grmompdgId[nGenLep]/I");
   
   Tout->Branch("nGenNu",&nGenNu, "nGenNu/I");
   Tout->Branch("GenNu_pt",GenNu_pt,"GenNu_pt[nGenNu]/F");
   Tout->Branch("GenNu_eta",GenNu_eta,"GenNu_eta[nGenNu]/F");
   Tout->Branch("GenNu_phi",GenNu_phi,"GenNu_phi[nGenNu]/F");
   Tout->Branch("GenNu_mass",GenNu_mass,"GenNu_mass[nGenNu]/F");
   Tout->Branch("GenNu_pdgId",GenNu_pdgId,"GenNu_pdgId[nGenNu]/I");
   Tout->Branch("GenNu_mompdgId",GenNu_mompdgId,"GenNu_mompdgId[nGenNu]/I");
   Tout->Branch("GenNu_grmompdgId",GenNu_grmompdgId,"GenNu_grmompdgId[nGenNu]/I");
   
   Tout->Branch("nGenBPart",&nGenBPart, "nGenBPart/I");
   Tout->Branch("GenBPart_pt",GenBPart_pt,"GenBPart_pt[nGenBPart]/F");
   Tout->Branch("GenBPart_eta",GenBPart_eta,"GenBPart_eta[nGenBPart]/F");
   Tout->Branch("GenBPart_phi",GenBPart_phi,"GenBPart_phi[nGenBPart]/F");
   Tout->Branch("GenBPart_mass",GenBPart_mass,"GenBPart_mass[nGenBPart]/F");
   Tout->Branch("GenBPart_pdgId",GenBPart_pdgId,"GenBPart_pdgId[nGenBPart]/I");
   Tout->Branch("GenBPart_mompdgId",GenBPart_mompdgId,"GenBPart_mompdgId[nGenBPart]/I");
   Tout->Branch("GenBPart_grmompdgId",GenBPart_grmompdgId,"GenBPart_grmompdgId[nGenNu]/I");
   
   Tout->Branch("nGenV",&nGenV, "nGenV/I");
   Tout->Branch("GenV_pt",GenV_pt,"GenV_pt[nGenV]/F");
   Tout->Branch("GenV_eta",GenV_eta,"GenV_eta[nGenV]/F");
   Tout->Branch("GenV_phi",GenV_phi,"GenV_phi[nGenV]/F");
   Tout->Branch("GenV_mass",GenV_mass,"GenV_mass[nGenV]/F");
   Tout->Branch("GenV_pdgId",GenV_pdgId,"GenV_pdgId[nGenV]/I");
   Tout->Branch("GenV_mompdgId",GenV_mompdgId,"GenV_mompdgId[nGenV]/I");
   Tout->Branch("GenV_grmompdgId",GenV_grmompdgId,"GenV_grmompdgId[nGenNu]/I");
   
   calib_deepflav = BTagCalibration("DeepJet", "BtagRecommendation106XUL18/DeepJet_106XUL18SF_WPonly_V1p1.csv");
   reader_deepflav = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"}); 
   reader_deepflav.load(calib_deepflav, BTagEntry::FLAV_B, "comb");
   reader_deepflav.load(calib_deepflav, BTagEntry::FLAV_C, "comb");
   reader_deepflav.load(calib_deepflav, BTagEntry::FLAV_UDSG, "incl");
   
   char name[1000];
   
   TFile *file_mu_sf;
   sprintf(name,"data/Efficiencies_muon_generalTracks_Z_Run%i_UL_ID.root",year);
   file_mu_sf = new TFile(name,"read");
   
   TFile *file_el_sf;
   sprintf(name,"data/egammaEffi.txt_Ele_%s_EGM2D_UL%i.root",electron_id_name.c_str(),year);
   file_el_sf = new TFile(name,"read");
   
   int count =0;
   string fileName;
   ifstream infile;
   infile.open(argv[3]);
   while(!infile.eof()){
   count = count+1;
   getline(infile,fileName);

   int L_lim = stof(argv[1]);
   int H_lim = stof(argv[2]);
   if(count<=L_lim)continue;
   if(count>H_lim)continue;
   TFile *f = TFile::Open(fileName.data());
   if(f==0) continue;

   cout<<fileName<<endl;

   TTree *fChain;
   fChain = (TTree*)f->Get("Events");
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
   fChain->SetBranchAddress("prefiringweight", &prefiringweight, &b_prefiringweight);
   fChain->SetBranchAddress("prefiringweightup", &prefiringweightup, &b_prefiringweightup);
   fChain->SetBranchAddress("prefiringweightdown", &prefiringweightdown, &b_prefiringweightdown);
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
   fChain->SetBranchAddress("PFJetAK8_JER", PFJetAK8_JER, &b_PFJetAK8_JER);
   fChain->SetBranchAddress("PFJetAK8_JERup", PFJetAK8_JERup, &b_PFJetAK8_JERup);
   fChain->SetBranchAddress("PFJetAK8_JERdn", PFJetAK8_JERdn, &b_PFJetAK8_JERdn);
   fChain->SetBranchAddress("PFJetAK8_msoftdrop", PFJetAK8_msoftdrop, &b_PFJetAK8_msoftdrop);
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
   fChain->SetBranchAddress("PFJetAK8_sub1JEC", PFJetAK8_sub1JEC, &b_PFJetAK8_sub1JEC);
   fChain->SetBranchAddress("PFJetAK8_sub2pt", PFJetAK8_sub2pt, &b_PFJetAK8_sub2pt);
   fChain->SetBranchAddress("PFJetAK8_sub2eta", PFJetAK8_sub2eta, &b_PFJetAK8_sub2eta);
   fChain->SetBranchAddress("PFJetAK8_sub2phi", PFJetAK8_sub2phi, &b_PFJetAK8_sub2phi);
   fChain->SetBranchAddress("PFJetAK8_sub2mass", PFJetAK8_sub2mass, &b_PFJetAK8_sub2mass);
   fChain->SetBranchAddress("PFJetAK8_sub2btag", PFJetAK8_sub2btag, &b_PFJetAK8_sub2btag);
   fChain->SetBranchAddress("PFJetAK8_sub2JEC", PFJetAK8_sub2JEC, &b_PFJetAK8_sub2JEC);
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
   fChain->SetBranchAddress("PFJetAK4_JER", PFJetAK4_JER, &b_PFJetAK4_JER);
   fChain->SetBranchAddress("PFJetAK4_JERup", PFJetAK4_JERup, &b_PFJetAK4_JERup);
   fChain->SetBranchAddress("PFJetAK4_JERdn", PFJetAK4_JERdn, &b_PFJetAK4_JERdn);
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
   fChain->SetBranchAddress("Muon_TightID", Muon_TightID, &b_Muon_Muon_TightID);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_p", Muon_p, &b_Muon_p);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_minisoch", Muon_minisoch, &b_Muon_minisoch);
   fChain->SetBranchAddress("Muon_minisonh", Muon_minisonh, &b_Muon_minisonh);
   fChain->SetBranchAddress("Muon_minisoph", Muon_minisoph, &b_Muon_minisoph);
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
   fChain->SetBranchAddress("Electron_minisoch", Electron_minisoch, &b_Electron_minisoch);
   fChain->SetBranchAddress("Electron_minisonh", Electron_minisonh, &b_Electron_minisonh);
   fChain->SetBranchAddress("Electron_minisoph", Electron_minisoph, &b_Electron_minisoph);
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
   fChain->SetBranchAddress("Generator_weight", &Generator_weight, &b_Generator_weight);
   fChain->SetBranchAddress("Generator_qscale", &Generator_qscale, &b_Generator_qscale);
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
   fChain->SetBranchAddress("LHE_weight", &LHE_weight, &b_LHE_weight);
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
   isFastSIM = true;
   
   for (int ij=0; ij<nentries; ij++) {
   
	//if(ij%100==0){ cout<<"event "<<ij+1<<endl; }
	f->cd();

	double event_weight;
	if(fabs(LHE_weight)>1.e-12) { event_weight = LHE_weight; }
	else { event_weight = Generator_weight;	}

	fChain->GetEntry(ij);
        h_PU->Fill(npu_vert,event_weight);
        h_TPU->Fill(npu_vert_true,event_weight);
        h_PV->Fill(npvert,event_weight);
	
	//cout<<"Before: "<<nPFJetAK8<<" "<<nPFJetAK4<<" "<<nMuon<<" "<<nElectron<<" "<<nPhoton<<endl;
	
	//Here you get gen particles
	vector<GenParton> genpartons;
	getPartons(genpartons);
	
	vector<GenParton> genleps;
	for(unsigned ig=0; ig<(genpartons).size(); ig++){
		if((abs(genpartons[ig].pdgId)==11||abs(genpartons[ig].pdgId)==13) && ((genpartons[ig].status)==1||(genpartons[ig].status)==23) && (abs(genpartons[ig].mompdgId)==24||abs(genpartons[ig].mompdgId)==25)){
			genleps.push_back(genpartons[ig]);
		}
	}
	
	vector<GenParton> gennus;
	for(unsigned ig=0; ig<(genpartons).size(); ig++){
		if((abs(genpartons[ig].pdgId)==12||abs(genpartons[ig].pdgId)==14||abs(genpartons[ig].pdgId)==16) && ((genpartons[ig].status)==1||(genpartons[ig].status)==23)){
			gennus.push_back(genpartons[ig]);
		}
	}
	
	vector<GenParton> genbs;
	for(unsigned ig=0; ig<(genpartons).size(); ig++){
		if((abs(genpartons[ig].pdgId)==5) && (genpartons[ig].status==23) && (abs(genpartons[ig].mompdgId)==25||abs(genpartons[ig].mompdgId)==6||abs(genpartons[ig].mompdgId)==35)){
			genbs.push_back(genpartons[ig]);
		}
	}
	
	vector<GenParton> genVs;
	for(unsigned ig=0; ig<(genpartons).size(); ig++){
		if(abs(genpartons[ig].pdgId)==23||abs(genpartons[ig].pdgId)==24||abs(genpartons[ig].pdgId)==25||abs(genpartons[ig].pdgId)==35){
			genVs.push_back(genpartons[ig]);
		}
	}
	
	
	vector<GenParton> lheparts;
	getLHEParts(lheparts);
	
	//Here you get electrons with your criteria
	vector <Electron> velectrons;
	getelectrons(velectrons,electron_pt_cut,absetacut,3); // index after eta cut is id: 0=WP80_Iso, 1=WP80_noIso, 2=WP90_Iso, 3=WP90_noIso
	
	//Here you get muons with your criteria (no iso used by default)
    vector <Muon> vmuons;
    getmuons(vmuons,muon_pt_cut,absetacut,2); // index after eta cut is id: 0=Loose, 1=Medium, 2=Tight
    
    //Make lepton collection from electrons & muons (using only common variables)
    vector <Lepton> vleptons;
    getLeptons(vleptons,vmuons,velectrons,lepton_pt_cut);    
    
    //Here you get AK4 jets with your criteria
    vector <AK4Jet> Jets;
    getAK4jets(Jets,AK4jet_pt_cut,absetacut,isMC);
    //LeptonJet_cleaning(Jets,vleptons,AK4jet_pt_cut,absetacut);
    
    // Add b tag SF (if not in ntuple)//
    for(auto & jet: Jets){
		
		BTagEntry::JetFlavor btv_flav;
		if(abs(jet.hadronFlavour)==5){ btv_flav = BTagEntry::FLAV_B; }
		else if (abs(jet.hadronFlavour)==4){ btv_flav = BTagEntry::FLAV_C; }
		else { btv_flav = BTagEntry::FLAV_UDSG; }
		
		jet.btag_DeepFlav_SF = reader_deepflav.eval_auto_bounds("central",btv_flav,fabs(jet.eta),jet.pt); 
		jet.btag_DeepFlav_SF_up = reader_deepflav.eval_auto_bounds("up",btv_flav,fabs(jet.eta),jet.pt);
		jet.btag_DeepFlav_SF_dn = reader_deepflav.eval_auto_bounds("down",btv_flav,fabs(jet.eta),jet.pt);
	
	}
	    
    //Get b-tagged jets from AK4 jets
    vector <AK4Jet> BJets;
    for(auto & jet: Jets){
		if(isBJet(jet,deep_btag_cut)){
			BJets.push_back(jet);
		}
    }
    
	//Here you get AK8 jets with your criteria
    vector <AK8Jet> LJets;
    getAK8jets(LJets,AK8jet_pt_cut,absetacut,isMC);
    //LeptonJet_cleaning(LJets,vleptons,AK8jet_pt_cut,absetacut);

    nleptons = (vleptons.size());
	nfatjets = (LJets.size());
    
    ////trigger object along with pdgid////
    vector<TrigObj> trigobjects;
    getTrigObjs(trigobjects);
    
    for (unsigned tr=0; tr<trigobjects.size(); tr++) {
		
        int trig_id = -1;
        if(*(decToBinary(trigobjects[tr].type)) == 1) { trig_id = 11; }
        else if (*(decToBinary(trigobjects[tr].type)+1)==1) { trig_id = 11; }
        else if (*(decToBinary(trigobjects[tr].type)+2)==1 || *(decToBinary(trigobjects[tr].type)+3)==1) { trig_id = 13; }
        else if (*(decToBinary(trigobjects[tr].type)+4)==1 || *(decToBinary(trigobjects[tr].type)+5)==1) { trig_id = 15; }
        else if (*(decToBinary(trigobjects[tr].type)+6)==1) { trig_id = 0; }
        else if (*(decToBinary(trigobjects[tr].type)+7)==1 || *(decToBinary(trigobjects[tr].type)+8)==1) { trig_id = 1; }
        //cout<<"pdg "<<TrigObj_pdgId[tr]<<" id "<<trig_id<<endl;
        trigobjects[tr].ID = trig_id;
    
    }

	// trigger decisions //
	vector<bool> double_hlts; vector<vector<float>> double_pt_cuts; vector<vector<int>> double_pids;
    vector<bool> single_hlts; vector<float> single_pt_cuts; vector<int> single_pids; vector<float> single_other_pt_cuts; vector<int> single_other_pids;
	vector<bool> jet_hlts; vector<float> jet_pt_cuts; vector<int> jet_pids;
	
	// For now offline cuts for leptons is kept 3 GeV higher than at trigger-level, these numbers need to be finalized from trigger efficiency plots //
		
	double_hlts.push_back(hlt_Mu37_Ele27_CaloIdL_MW);
    double_pt_cuts.push_back({37+3,27+3});
    double_pids.push_back({13,11});
    double_hlts.push_back(hlt_Mu27_Ele37_CaloIdL_MW);
    double_pt_cuts.push_back({37+3,27+3});
    double_pids.push_back({11,13});
	double_hlts.push_back(hlt_DoubleEle25_CaloIdL_MW);
    double_pt_cuts.push_back({25+15,25+5});
    double_pids.push_back({11,11});
	double_hlts.push_back(hlt_Mu37_TkMu27);
    double_pt_cuts.push_back({37+3,27+3});
    double_pids.push_back({13,13});
    
    single_hlts.push_back(hlt_Mu50);
    single_pt_cuts.push_back(50+3);
    single_pids.push_back(13);
    single_other_pt_cuts.push_back(-100);
    single_other_pids.push_back(0);
    single_hlts.push_back(hlt_Ele32_WPTight_Gsf);
    single_pt_cuts.push_back(32+3);
    single_pids.push_back(11);
    single_other_pt_cuts.push_back(-100);
    single_other_pids.push_back(0);
    single_hlts.push_back(hlt_Ele40_WPTight_Gsf);
    single_pt_cuts.push_back(40+3);
    single_pids.push_back(11);
    single_other_pt_cuts.push_back(-100);
    single_other_pids.push_back(0);
    single_hlts.push_back(hlt_Ele115_CaloIdVT_GsfTrkIdT);
    single_pt_cuts.push_back(115+3);
    single_pids.push_back(11);
    single_other_pt_cuts.push_back(-100);
    single_other_pids.push_back(0);
    single_hlts.push_back(hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
    single_pt_cuts.push_back(50+3);
    single_pids.push_back(11);
    single_other_pt_cuts.push_back(165+15);
    single_other_pids.push_back(1);

	jet_hlts.push_back(hlt_AK8PFJet500);
	jet_pt_cuts.push_back(500+50);
	jet_pids.push_back(2);
	jet_hlts.push_back(hlt_PFJet500);
	jet_pt_cuts.push_back(500+50);
	jet_pids.push_back(1);
	
	bool anytrig_pass(false);
	bool trig_threshold_pass(false); 
	bool trig_matching_pass(false); 
	bool muon_trig_pass(false);
	bool electron_trig_pass(false);
	
	if(isFastSIM){
				
		anytrig_pass = true; 
		if(vmuons.size()>0 && vmuons[0].pt>53) { trig_threshold_pass = true; trig_matching_pass = true; muon_trig_pass = true; }
		else if(velectrons.size()>0 && velectrons[0].pt>35) { trig_threshold_pass = true; trig_matching_pass = true; electron_trig_pass = true; }
		else if(Jets.size()>0 && Jets[0].pt>550) { trig_threshold_pass = true; trig_matching_pass = true; }
		else if(LJets.size()>0 && LJets[0].pt>550) { trig_threshold_pass = true; trig_matching_pass = true; }
		
		}
	else{
		Match_trigger(double_hlts, double_pt_cuts, double_pids,
					  single_hlts, single_pt_cuts, single_pids, single_other_pt_cuts, single_other_pids,
					  jet_hlts, jet_pt_cuts, jet_pids,
					  trigobjects,
					  vmuons, velectrons, vleptons, Jets, LJets,
					  anytrig_pass, trig_threshold_pass, trig_matching_pass, muon_trig_pass, electron_trig_pass
					);
		}
    //*****************************************************************************************
    //                            AK4 histogram filling for btag SF                          //
    //*****************************************************************************************                            
    if(anytrig_pass && trig_threshold_pass && trig_matching_pass && !(muon_trig_pass && electron_trig_pass))
    {
    vector <AK4GenJet> genJets;
    getAK4Genjets(genJets,AK4GenJet_pt_cut,absetacut,isMC);

    for(unsigned ijet=0; ijet<Jets.size(); ijet++){
           double dR = 9999.9;
           for(unsigned gjet=0; gjet<genJets.size(); gjet++)
	     {
                 double temp_dR = delta2R(Jets[ijet].y,Jets[ijet].phi,genJets[gjet].p4.Rapidity(),genJets[gjet].phi) ;
		 if (temp_dR < dR )
		 {
			 dR = temp_dR;
		 }
	     }
            if(dR < 0.4)
	              {
                           if( fabs(Jets[ijet].hadronFlavour) == 5 )  {  h_Ak4_b_flv->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),event_weight); 
				                                   if (Jets[ijet].btag_DeepFlav   > DAK4_T  )  { h_Ak4_b_flv_pass_T->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),event_weight); } 
                                                                   if (Jets[ijet].btag_DeepFlav   > DAK4_M ) { h_Ak4_b_flv_pass_M->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),event_weight); } 
			                                           if (Jets[ijet].btag_DeepFlav   > DAK4_L ) { h_Ak4_b_flv_pass_L->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),event_weight); }
			                                        }    
                            else if( fabs(Jets[ijet].hadronFlavour) == 4 )  {  h_Ak4_c_flv->Fill(Jets[ijet].pt,Jets[ijet].eta,LHE_weight);  
                                                                   if (Jets[ijet].btag_DeepFlav   > DAK4_T  )  { h_Ak4_c_flv_pass_T->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),event_weight); }
                                                                   if (Jets[ijet].btag_DeepFlav   > DAK4_M ) { h_Ak4_c_flv_pass_M->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),event_weight); }
                                                                   if (Jets[ijet].btag_DeepFlav   > DAK4_L ) { h_Ak4_c_flv_pass_L->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),event_weight); }
                                                                }
                           if( Jets[ijet].hadronFlavour == 0 )  {  h_Ak4_l_flv->Fill(Jets[ijet].pt,Jets[ijet].eta,LHE_weight);  
                                                                   if (Jets[ijet].btag_DeepFlav   > DAK4_T  )  { h_Ak4_l_flv_pass_T->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),event_weight); }
                                                                   if (Jets[ijet].btag_DeepFlav   > DAK4_M ) { h_Ak4_l_flv_pass_M->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),event_weight); }
                                                                   if (Jets[ijet].btag_DeepFlav   > DAK4_L ) { h_Ak4_l_flv_pass_L->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),event_weight); }
                                                                }
		      }		                                   
             }
    //*****************************************************************************************
    //                                           AK8 histogram filling                       //
    //*****************************************************************************************                            
    vector <AK8GenJet> genLJets;
    getAK8Genjets(genLJets,AK8GenJet_pt_cut,absetacut,isMC);
    for(unsigned ijet=0; ijet<LJets.size(); ijet++){
         double dR = 9999.9;
         for(unsigned gjet=0; gjet<genLJets.size(); gjet++)
             {
		     double temp_dR = delta2R(LJets[ijet].y,LJets[ijet].phi,genLJets[gjet].p4.Rapidity(),genLJets[gjet].phi) ;
		     if(temp_dR < dR )
                     {
                         dR = temp_dR;
                     }
            }
	 if(dR < 0.4) 
	   { 
              h_Ak8_DeepTag_PNetMD_WvsQCD->Fill(LJets[ijet].pt,LJets[ijet].eta,event_weight);
			  h_Ak8_DeepTag_PNetMD_XbbvsQCD->Fill(LJets[ijet].pt,LJets[ijet].eta,LHE_weight);
		      if(LJets[ijet].DeepTag_PNetMD_WvsQCD>PNetW_cut_T) {h_Ak8_DeepTag_PNetMD_WvsQCD_pass_T->Fill(LJets[ijet].pt,fabs(LJets[ijet].eta),event_weight);  }
	          if(LJets[ijet].DeepTag_PNetMD_WvsQCD>PNetW_cut_M) {h_Ak8_DeepTag_PNetMD_WvsQCD_pass_M->Fill(LJets[ijet].pt,fabs(LJets[ijet].eta),event_weight);  }
              if(LJets[ijet].DeepTag_PNetMD_WvsQCD>PNetW_cut_L) {h_Ak8_DeepTag_PNetMD_WvsQCD_pass_L->Fill(LJets[ijet].pt,fabs(LJets[ijet].eta),event_weight);  }
              if(LJets[ijet].DeepTag_PNetMD_XbbvsQCD>PNetbb_cut_T){h_Ak8_DeepTag_PNetMD_XbbvsQCD_pass_T->Fill(LJets[ijet].pt,fabs(LJets[ijet].eta),event_weight);}
              if(LJets[ijet].DeepTag_PNetMD_XbbvsQCD>PNetbb_cut_M){h_Ak8_DeepTag_PNetMD_XbbvsQCD_pass_M->Fill(LJets[ijet].pt,fabs(LJets[ijet].eta),event_weight);}
              if(LJets[ijet].DeepTag_PNetMD_XbbvsQCD>PNetbb_cut_L){h_Ak8_DeepTag_PNetMD_XbbvsQCD_pass_L->Fill(LJets[ijet].pt,fabs(LJets[ijet].eta),event_weight);}
           }
    }
    }
    // Event selection cuts //

    vector<bool> event_cuts;
    
    event_cuts.push_back(vleptons.size()==1); // for single lepton channel 
    event_cuts.push_back(LJets.size()>=2);	// At least two AK8 jets 
    event_cuts.push_back(PuppiMET_pt>=20);	// Minimum MET cut
	event_cuts.push_back(anytrig_pass);		// At least one trigger should be fired 
	event_cuts.push_back(trig_threshold_pass);	// Offline objects should pass trigger threshold 
	event_cuts.push_back(trig_matching_pass);	// Offline objects should match to trigger object
	event_cuts.push_back(!(muon_trig_pass && electron_trig_pass));   // reject if both mu & el triggers are fired (for single lepton channel)
	
	for(unsigned icut=0; icut<event_cuts.size(); icut++){
		Flag_event_cuts[icut] = event_cuts[icut];
	}
		
    bool event_pass = true;
	for(unsigned icut=0; icut<event_cuts.size(); icut++){
		event_pass *= event_cuts[icut];
		if(!event_pass) break;
    }

	if(!event_pass) continue;
	
    //cout<<"After: "<<LJets.size()<<" "<<Jets.size()<<" "<<vmuons.size()<<" "<<velectrons.size()<<" "<<nPhoton<<endl;
    
    // store gen particles first //
    
    nGenLep = int(genleps.size());
    for(unsigned ig=0; ig<(genleps.size()); ig++){
		GenLep_pt[ig] = genleps[ig].pt;
		GenLep_eta[ig] = genleps[ig].eta;
		GenLep_phi[ig] = genleps[ig].phi;
		GenLep_mass[ig] = genleps[ig].mass;
		GenLep_pdgId[ig] = genleps[ig].pdgId;
		GenLep_mompdgId[ig] = genleps[ig].mompdgId;
		GenLep_grmompdgId[ig] = genleps[ig].grmompdgId;
		}
    
    nGenNu = int(gennus.size());
    for(unsigned ig=0; ig<(gennus.size()); ig++){
		GenNu_pt[ig] = gennus[ig].pt;
		GenNu_eta[ig] = gennus[ig].eta;
		GenNu_phi[ig] = gennus[ig].phi;
		GenNu_mass[ig] = gennus[ig].mass;
		GenNu_pdgId[ig] = gennus[ig].pdgId;
		GenNu_mompdgId[ig] = gennus[ig].mompdgId;
		GenNu_grmompdgId[ig] = gennus[ig].grmompdgId;
		}
    
    nGenBPart = int(genbs.size());
    for(unsigned ig=0; ig<(genbs.size()); ig++){
		GenBPart_pt[ig] = genbs[ig].pt;
		GenBPart_eta[ig] = genbs[ig].eta;
		GenBPart_phi[ig] = genbs[ig].phi;
		GenBPart_mass[ig] = genbs[ig].mass;
		GenBPart_pdgId[ig] = genbs[ig].pdgId;
		GenBPart_mompdgId[ig] = genbs[ig].mompdgId;
		GenBPart_grmompdgId[ig] = genbs[ig].grmompdgId;
		}
    
    nGenV = int(genVs.size());
    for(unsigned ig=0; ig<(genVs.size()); ig++){
		GenV_pt[ig] = genVs[ig].pt;
		GenV_eta[ig] = genVs[ig].eta;
		GenV_phi[ig] = genVs[ig].phi;
		GenV_mass[ig] = genVs[ig].mass;
		GenV_pdgId[ig] = genVs[ig].pdgId;
		GenV_mompdgId[ig] = genVs[ig].mompdgId;
		GenV_grmompdgId[ig] = genVs[ig].grmompdgId;
		}
    
    // check for reco objects //
		
	int Y_cand = -1;
	int W_cand_opt1 = -1;
        int W_cand_opt2 = -1;

	float max_PNet_bb = -100;
	float max_PNet_W = -100;
	// Y candidate
	for(unsigned ijet=0; ijet<LJets.size(); ijet++){
		
		if(LJets[ijet].DeepTag_PNetMD_XbbvsQCD > max_PNet_bb){
			
			max_PNet_bb = LJets[ijet].DeepTag_PNetMD_XbbvsQCD;
			Y_cand = int(ijet);
			
		}
		
	}
        //W candidate option 1
	for(unsigned ijet=0; ijet<LJets.size(); ijet++){
		
		if(LJets[ijet].DeepTag_PNetMD_WvsQCD > max_PNet_W && int(ijet)!=Y_cand){
			
			max_PNet_W = LJets[ijet].DeepTag_PNetMD_WvsQCD;
			W_cand_opt1 = int(ijet);
			
		}
	}
	if(Y_cand==-1 || W_cand_opt1==-1) continue;

	// W candidate option 2
	double dR_LJet_lep =  9999.9;
	for(unsigned ijet=0; ijet<LJets.size(); ijet++){
	        if (int(ijet) == Y_cand) continue;
                double temp_dR = delta2R(LJets[ijet].y,LJets[ijet].phi,vleptons[0].eta,vleptons[0].phi);
		if( temp_dR < dR_LJet_lep )
                  {
      		      dR_LJet_lep =  temp_dR;
	              W_cand_opt2 =  int(ijet);
	          }
         }  		
         if(Y_cand==-1 || W_cand_opt2==-1) continue;



	TLorentzVector pnu_opt1, pnu_opt2;
	double random_no = gRandom->Uniform(0,1);
	
	// regions can be constructed using: 
	// 1. MET
	// 2. bb score of Y candidate
	// 3. W score of W candidate
	// 4. Mass of H candidate (l+nu+W)
	// 5. dR of lepton & W candidate
	// 6. dphi & dR between lepton & Y candidate
	
	bool Y_bb_pass_T = (LJets[Y_cand].DeepTag_PNetMD_XbbvsQCD > PNetbb_cut_T);
    bool Y_bb_pass_M = (LJets[Y_cand].DeepTag_PNetMD_XbbvsQCD > PNetbb_cut_M);
	bool Y_bb_pass_L = (LJets[Y_cand].DeepTag_PNetMD_XbbvsQCD > PNetbb_cut_L);

	bool H_W_pass_T_opt1 = (LJets[W_cand_opt1].DeepTag_PNetMD_WvsQCD > PNetW_cut_T);
    bool H_W_pass_M_opt1 = (LJets[W_cand_opt1].DeepTag_PNetMD_WvsQCD > PNetW_cut_M);
	bool H_W_pass_L_opt1 = (LJets[W_cand_opt1].DeepTag_PNetMD_WvsQCD > PNetW_cut_L);

    bool H_W_pass_T_opt2 = (LJets[W_cand_opt2].DeepTag_PNetMD_WvsQCD > PNetW_cut_T);
    bool H_W_pass_M_opt2 = (LJets[W_cand_opt2].DeepTag_PNetMD_WvsQCD > PNetW_cut_M);
    bool H_W_pass_L_opt2 = (LJets[W_cand_opt2].DeepTag_PNetMD_WvsQCD > PNetW_cut_L);

	pnu_opt1 = neutrino_mom_fromH(vleptons[0].p4+LJets[W_cand_opt1].p4, PuppiMET_pt, PuppiMET_phi, random_no);

	bool H_m_pass_opt1 = ((vleptons[0].p4+LJets[W_cand_opt1].p4+pnu_opt1).M()>90. && (vleptons[0].p4+LJets[W_cand_opt1].p4+pnu_opt1).M()<150.);
	bool dR_lW_pass_opt1 = (delta2R(LJets[W_cand_opt1].y,LJets[W_cand_opt1].phi,vleptons[0].eta,vleptons[0].phi) < 1.2);

	pnu_opt2 = neutrino_mom_fromH(vleptons[0].p4+LJets[W_cand_opt2].p4, PuppiMET_pt, PuppiMET_phi, random_no);

    bool H_m_pass_opt2 = ((vleptons[0].p4+LJets[W_cand_opt2].p4+pnu_opt2).M()>90. && (vleptons[0].p4+LJets[W_cand_opt2].p4+pnu_opt2).M()<150.);
    bool dR_lW_pass_opt2 = (delta2R(LJets[W_cand_opt2].y,LJets[W_cand_opt2].phi,vleptons[0].eta,vleptons[0].phi) < 1.2);
	
    bool MET_pass = (PuppiMET_pt > 50);
    
    bool SR_opt1 = (Y_bb_pass_T && H_W_pass_T_opt1 && H_m_pass_opt1 && dR_lW_pass_opt1 && MET_pass);
    bool Wj_CR_opt1 = (!Y_bb_pass_T && H_W_pass_T_opt1 && H_m_pass_opt1 && !dR_lW_pass_opt1 && MET_pass);

	bool SR_opt2 = (Y_bb_pass_T && H_W_pass_T_opt2 && H_m_pass_opt2 && dR_lW_pass_opt2 && MET_pass);
    bool Wj_CR_opt2 = (!Y_bb_pass_T && H_W_pass_T_opt2 && H_m_pass_opt2 && !dR_lW_pass_opt2 && MET_pass);
    
    if(vleptons.size()>0){
		
		l_pt = vleptons[0].pt;
		l_eta = vleptons[0].eta;
		l_phi = vleptons[0].phi;
		l_mass = vleptons[0].mass;
		
		l_genindex = get_nearest_Parton(genleps,vleptons[0].p4,0.4);

	}
    
    if(Y_cand>=0) {
		Y_pt = LJets[Y_cand].pt;
		Y_y = LJets[Y_cand].y;
        Y_eta = LJets[Y_cand].eta;
		Y_phi = LJets[Y_cand].phi;
		Y_mass = LJets[Y_cand].mass;
		Y_sdmass = LJets[Y_cand].sdmass;
		Y_PN_bb = LJets[Y_cand].DeepTag_PNetMD_XbbvsQCD;
		
		Y_genbindex[0] = get_nearest_Parton(genbs,LJets[Y_cand].p4,0.8);
		if(Y_genbindex[0]>=0){
			genbs.erase(genbs.begin()+Y_genbindex[0]);
			Y_genbindex[1] = get_nearest_Parton(genbs,LJets[Y_cand].p4,0.8);
		}
		else{
			Y_genbindex[1] = -1;
			}
		
		int gen_match = get_nearest_Parton(genVs,LJets[Y_cand].p4,0.8);
		if(gen_match>=0 && abs(genVs[gen_match].pdgId)==35){
			Y_genindex = gen_match; 
		}else{
			Y_genindex = -1;
			}
		
		Y_JESup = LJets[Y_cand].jesup_Total;
		Y_JESdn = LJets[Y_cand].jesdn_Total;
		Y_JERup = LJets[Y_cand].JERup/LJets[Y_cand].JER;
		Y_JERdn = LJets[Y_cand].JERdn/LJets[Y_cand].JER;
				
		if(vleptons.size()>0){
			
			dR_lY = delta2R(LJets[Y_cand].eta,LJets[Y_cand].phi,vleptons[0].eta,vleptons[0].phi);
			dy_lY = (LJets[Y_cand].eta - vleptons[0].eta);
			dphi_lY = PhiInRange(LJets[Y_cand].phi - vleptons[0].p4.Phi());
			
		}
	}
	
	if(W_cand_opt1>=0) {
		
		W_pt_opt1 = LJets[W_cand_opt1].pt;
		W_y_opt1 = LJets[W_cand_opt1].y;
        W_eta_opt1 = LJets[W_cand_opt1].eta;
	    W_phi_opt1 = LJets[W_cand_opt1].phi;
		W_mass_opt1 = LJets[W_cand_opt1].mass;
		W_sdmass_opt1 = LJets[W_cand_opt1].sdmass;
		W_DAK8_W_opt1 = LJets[W_cand_opt1].DeepTag_DAK8MD_WvsQCD;
		W_PN_W_opt1 = LJets[W_cand_opt1].DeepTag_PNetMD_WvsQCD;
		
		int gen_match = get_nearest_Parton(genVs,LJets[W_cand_opt1].p4,0.8);
		if(gen_match>=0 && abs(genVs[gen_match].pdgId)==24){
			W_genindex_opt1 = gen_match; 
		}else{
			W_genindex_opt1 = -1;
			}
		
		W_JESup_opt1 = LJets[W_cand_opt1].jesup_Total;
		W_JESdn_opt1 = LJets[W_cand_opt1].jesdn_Total;
		W_JERup_opt1 = LJets[W_cand_opt1].JERup/LJets[W_cand_opt1].JER;
		W_JERdn_opt1 = LJets[W_cand_opt1].JERdn/LJets[W_cand_opt1].JER;
		
		if(vleptons.size()>0 && pnu_opt1.Eta()>-100){
			
			TLorentzVector W_mom = LJets[W_cand_opt1].p4;
			TLorentzVector H_mom = (W_mom + vleptons[0].p4 + pnu_opt1);
			
			H_pt_opt1 = H_mom.Pt();
			H_y_opt1 = H_mom.Rapidity();
            H_eta_opt1 = H_mom.Eta();
			H_phi_opt1 = H_mom.Phi();
			H_mass_opt1 = H_mom.M();
			
			int gen_match = get_nearest_Parton(genVs,H_mom,0.8);
			if(gen_match>=0 && abs(genVs[gen_match].pdgId)==25){
				H_genindex_opt1 = gen_match; 
			}else{
				H_genindex_opt1 = -1;
			}
			
			W_mom.SetPtEtaPhiM(LJets[W_cand_opt1].jesup_Total*LJets[W_cand_opt1].p4.Pt(),LJets[W_cand_opt1].p4.Eta(),LJets[W_cand_opt1].p4.Phi(),LJets[W_cand_opt1].jesup_Total*LJets[W_cand_opt1].p4.M());
			H_JESup_opt1 = (W_mom + vleptons[0].p4 + pnu_opt1).Pt()/H_mom.Pt();
			
			W_mom.SetPtEtaPhiM(LJets[W_cand_opt1].jesdn_Total*LJets[W_cand_opt1].p4.Pt(),LJets[W_cand_opt1].p4.Eta(),LJets[W_cand_opt1].p4.Phi(),LJets[W_cand_opt1].jesdn_Total*LJets[W_cand_opt1].p4.M());
			H_JESdn_opt1 = (W_mom + vleptons[0].p4 + pnu_opt1).Pt()/H_mom.Pt();
			
			W_mom.SetPtEtaPhiM(LJets[W_cand_opt1].JERup*LJets[W_cand_opt1].p4.Pt(),LJets[W_cand_opt1].p4.Eta(),LJets[W_cand_opt1].p4.Phi(),LJets[W_cand_opt1].p4.M());
			H_JERup_opt1 = (W_mom + vleptons[0].p4 + pnu_opt1).Pt()/H_mom.Pt();
			
			W_mom.SetPtEtaPhiM(LJets[W_cand_opt1].JERdn*LJets[W_cand_opt1].p4.Pt(),LJets[W_cand_opt1].p4.Eta(),LJets[W_cand_opt1].p4.Phi(),LJets[W_cand_opt1].p4.M());
			H_JERdn_opt1 = (W_mom + vleptons[0].p4 + pnu_opt1).Pt()/H_mom.Pt();
			
			if(Y_cand>=0) {
				X_mass_opt1 = (LJets[Y_cand].p4 + LJets[W_cand_opt1].p4 + vleptons[0].p4 + pnu_opt1).M();
				}
		
			}
	}
        if(W_cand_opt2>=0) {

                W_pt_opt2 = LJets[W_cand_opt2].pt;
                W_y_opt2 = LJets[W_cand_opt2].y;
                W_eta_opt2 = LJets[W_cand_opt2].eta;
                W_phi_opt2 = LJets[W_cand_opt2].phi;
                W_mass_opt2 = LJets[W_cand_opt2].mass;
                W_sdmass_opt2 = LJets[W_cand_opt2].sdmass;
                W_DAK8_W_opt2 = LJets[W_cand_opt2].DeepTag_DAK8MD_WvsQCD;
                W_PN_W_opt2 = LJets[W_cand_opt2].DeepTag_PNetMD_WvsQCD;
                
                int gen_match = get_nearest_Parton(genVs,LJets[W_cand_opt2].p4,0.8);
				if(gen_match>=0 && abs(genVs[gen_match].pdgId)==24){
					W_genindex_opt2 = gen_match; 
				}else{
					W_genindex_opt2 = -1;
				}
                
                W_JESup_opt2 = LJets[W_cand_opt2].jesup_Total;
				W_JESdn_opt2 = LJets[W_cand_opt2].jesdn_Total;
				W_JERup_opt2 = LJets[W_cand_opt2].JERup/LJets[W_cand_opt2].JER;
				W_JERdn_opt2 = LJets[W_cand_opt2].JERdn/LJets[W_cand_opt2].JER;

                if(vleptons.size()>0 && pnu_opt2.Eta()>-100){

						TLorentzVector W_mom = LJets[W_cand_opt2].p4;
						TLorentzVector H_mom = (W_mom + vleptons[0].p4 + pnu_opt2);

                        H_pt_opt2 = H_mom.Pt();
                        H_y_opt2 = H_mom.Rapidity();
                        H_eta_opt2 = H_mom.Eta();
                        H_phi_opt2 = H_mom.Phi();
                        H_mass_opt2 = H_mom.M();
                        
                        int gen_match = get_nearest_Parton(genVs,H_mom,0.8);
						if(gen_match>=0 && abs(genVs[gen_match].pdgId)==25){
							H_genindex_opt2 = gen_match; 
						}else{
							H_genindex_opt2 = -1;
						}
                        
                        W_mom.SetPtEtaPhiM(LJets[W_cand_opt2].jesup_Total*LJets[W_cand_opt2].p4.Pt(),LJets[W_cand_opt2].p4.Eta(),LJets[W_cand_opt2].p4.Phi(),LJets[W_cand_opt2].jesup_Total*LJets[W_cand_opt2].p4.M());
						H_JESup_opt2 = (W_mom + vleptons[0].p4 + pnu_opt2).Pt()/H_mom.Pt();
						
						W_mom.SetPtEtaPhiM(LJets[W_cand_opt2].jesdn_Total*LJets[W_cand_opt2].p4.Pt(),LJets[W_cand_opt2].p4.Eta(),LJets[W_cand_opt2].p4.Phi(),LJets[W_cand_opt2].jesdn_Total*LJets[W_cand_opt2].p4.M());
						H_JESdn_opt2 = (W_mom + vleptons[0].p4 + pnu_opt2).Pt()/H_mom.Pt();
						
						W_mom.SetPtEtaPhiM(LJets[W_cand_opt2].JERup*LJets[W_cand_opt2].p4.Pt(),LJets[W_cand_opt2].p4.Eta(),LJets[W_cand_opt2].p4.Phi(),LJets[W_cand_opt2].p4.M());
						H_JERup_opt2 = (W_mom + vleptons[0].p4 + pnu_opt2).Pt()/H_mom.Pt();
						
						W_mom.SetPtEtaPhiM(LJets[W_cand_opt2].JERdn*LJets[W_cand_opt2].p4.Pt(),LJets[W_cand_opt2].p4.Eta(),LJets[W_cand_opt2].p4.Phi(),LJets[W_cand_opt2].p4.M());
						H_JERdn_opt2 = (W_mom + vleptons[0].p4 + pnu_opt2).Pt()/H_mom.Pt();

                        if(Y_cand>=0) {
                                X_mass_opt2 = (LJets[Y_cand].p4 + LJets[W_cand_opt2].p4 + vleptons[0].p4 + pnu_opt2).M();
                                }

                        }
        }

	
	dR_lW_opt1 = delta2R(LJets[W_cand_opt1].eta,LJets[W_cand_opt1].phi,vleptons[0].eta,vleptons[0].phi);
	dy_lW_opt1 = (LJets[W_cand_opt1].eta - vleptons[0].eta);
	dphi_lW_opt1 = PhiInRange(LJets[W_cand_opt1].phi - vleptons[0].p4.Phi());

        dR_lW_opt2 = delta2R(LJets[W_cand_opt2].eta,LJets[W_cand_opt2].phi,vleptons[0].eta,vleptons[0].phi);
        dy_lW_opt2 = (LJets[W_cand_opt2].eta - vleptons[0].eta);
        dphi_lW_opt2 = PhiInRange(LJets[W_cand_opt2].phi - vleptons[0].p4.Phi());

	nbjets_other = 0;
	for(auto & bjet: BJets){
		if(delta2R(bjet.y,bjet.phi,Jets[W_cand_opt1].y,LJets[W_cand_opt1].phi)>0.6 && delta2R(bjet.y,bjet.phi,Jets[Y_cand].y,LJets[Y_cand].phi)>0.6){
			nbjets_other++;
			}
	}

	MET = pnu_opt1.Pt();

	Flag_Y_bb_pass_T = Y_bb_pass_T;
        Flag_Y_bb_pass_M = Y_bb_pass_M;
	Flag_Y_bb_pass_L = Y_bb_pass_L;
	Flag_H_W_pass_T_opt1 = H_W_pass_T_opt1;
        Flag_H_W_pass_M_opt1 = H_W_pass_M_opt1;
	Flag_H_W_pass_L_opt1 = H_W_pass_L_opt1;

	Flag_H_W_pass_T_opt2 = H_W_pass_T_opt2;
        Flag_H_W_pass_M_opt2 = H_W_pass_M_opt2;
        Flag_H_W_pass_L_opt2 = H_W_pass_L_opt2;

	Flag_H_m_pass_opt1 = H_m_pass_opt1;
	Flag_dR_lW_pass_opt1 = dR_lW_pass_opt1;

        Flag_H_m_pass_opt2 = H_m_pass_opt2;
        Flag_dR_lW_pass_opt2 = dR_lW_pass_opt2;	

	Flag_MET_pass = MET_pass;
	
	Reg_SR_opt1 = SR_opt1;
	Reg_Wj_CR_opt1 = Wj_CR_opt1;
        Reg_SR_opt2 = SR_opt2;
        Reg_Wj_CR_opt2 = Wj_CR_opt2;	
	
    if(npu_vert_true>=0 && npu_vert_true<100){
		puWeight = pu_rat18[npu_vert_true];
		puWeightup = pu_rat18_up[npu_vert_true];
		puWeightdown = pu_rat18_dn[npu_vert_true];
	}
	
	leptonsf_weight = 1.0;
	leptonsf_weight_stat = 1.0;
	leptonsf_weight_syst = 1.0;
	
	for(unsigned lep=0; lep<vleptons.size(); lep++){
		if(abs(vleptons[lep].pdgId)==11) { 
			float *sfvalues = Electron_SF(file_el_sf, vleptons[lep].pt, vleptons[lep].eta);
			leptonsf_weight *= sfvalues[0];
			leptonsf_weight_stat *= (sfvalues[0] + sqrt(sfvalues[1]*sfvalues[1] + sfvalues[2]*sfvalues[2]));  // like this for time being 
			leptonsf_weight_syst *= (sfvalues[0] + sqrt(sfvalues[3]*sfvalues[3] + sfvalues[4]*sfvalues[4] + sfvalues[5]*sfvalues[5] + sfvalues[6]*sfvalues[6]));  // like this for time being 
		}
		if(abs(vleptons[lep].pdgId)==13) { 
			float *sfvalues;
			sfvalues = Muon_SF(file_mu_sf, muon_id_name, vleptons[lep].pt, vleptons[lep].eta);
			leptonsf_weight *= *(sfvalues+0);
			leptonsf_weight_stat *= *(sfvalues+1);
			leptonsf_weight_syst *= *(sfvalues+2);
		}
	}
	
	weight = 1;
	if(!isFastSIM){
		weight *= Generator_weight;
	}
	weight *= puWeight;
	weight *= leptonsf_weight;
	weight *= prefiringweight;
	
	// LHE story //
	/*
	for(unsigned ilhe=0; ilhe<lheparts.size(); ++ilhe){
		if(abs(lheparts[ilhe].pdgId)>10){
			cout<<ilhe+1<<" ID "<<lheparts[ilhe].pdgId<<" pt "<<lheparts[ilhe].p4.Pt()<<endl;	
			}	
		}
	cout<<"Lep "<<vleptons[0].pt<<" Nu "<<pnu.Pt()<<endl;
	*/
	Tout->Fill();
	}// entry
      f->Close();
      }
        infile.close();
        fileout->cd();
	fileout->Write();
	fileout->Close();
}

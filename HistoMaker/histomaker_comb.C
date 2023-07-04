#include "histomaker_comb.h"

using namespace std;

string input_path = "";

//int main()
//void histomaker_comb()
int main(int argc, char *argv[])
{
  if((argc-1)!=4){
	cout<<"Need exactly 4 arguments. Exiting!"<<endl;
	return 0;
  }
   
 std::istringstream(argv[1]) >> isDL; 
 std::istringstream(argv[2]) >> isDATA; 
 std::istringstream(argv[4]) >> isSignal; 
 cout<<"Running with options: isDileptonic? "<<isDL<<" isDATA? "<<isDATA<<" Signal? "<<isSignal<<endl;
 cout<<"Running on file : " << argv[3] << std::endl;	

 if(isDL) { input_path = "/eos/user/m/mukherje/XToYH/DL/"; }
 else { input_path = "/eos/user/m/mukherje/XToYH/SL/"; }

 //int nproc = sizeof(inputFile)/sizeof(inputFile[0]);
 
 //for (int ii=0;ii<nproc;ii++)
 // {

  //xbin histogram binning definitions
   int nmsdbins;
   if(isDL)
          {
              const int nmsdbins_dilep = sizeof(msdbins_dilep)/sizeof(msdbins_dilep[0])-1;
              nmsdbins  = nmsdbins_dilep;
              std::cout << "MY for dilep " << nmsdbins  << " :-->  " ;
          }
   else
         {
              const int nmsdbins_semilep = sizeof(msdbins_semilep)/sizeof(msdbins_semilep[0])-1;
              nmsdbins  = nmsdbins_semilep;
              std::cout << "MY for semilep " << nmsdbins  << " :-->  " ;
         }
   float msdbins[nmsdbins+1];
   if(isDL)
          {
              for(int aa = 0; aa < nmsdbins+1; aa++)
              {
                     std::cout << msdbins_dilep[aa] << " | " ;
                     msdbins[aa] = msdbins_dilep[aa];
              }
          }
   else
         {
              for(int aa = 0; aa < nmsdbins+1; aa++)
              {
                     std::cout << msdbins_semilep[aa] << " | " ;
                     msdbins[aa] = msdbins_semilep[aa];
              }
         }
   std::cout << std::endl;
   std::cout << "-----------------------------------------" << std::endl;
   int ninvmassbins;
   if(isDL)
          {
              const int ninvmassbins_dilep = sizeof(invmassbins_dilep)/sizeof(invmassbins_dilep[0])-1;
              ninvmassbins  = ninvmassbins_dilep;
              std::cout << "ST for dilep " << ninvmassbins  << " :-->  " ;
          }
   else
         {
              const int ninvmassbins_semilep = sizeof(invmassbins_semilep)/sizeof(invmassbins_semilep[0])-1;
              ninvmassbins  = ninvmassbins_semilep;
              std::cout << "MX for semilep " << ninvmassbins  << " :-->  " ;
         }
   float invmassbins[ninvmassbins+1];
   if(isDL)
          {
              for(int bb = 0; bb < ninvmassbins+1; bb++)
              {
                     std::cout << invmassbins_dilep[bb] << " | " ;
                     invmassbins[bb] = invmassbins_dilep[bb];
              }
          }
   else
         {
              for(int bb = 0; bb < ninvmassbins+1; bb++)
              {
                     std::cout << invmassbins_semilep[bb] << " | " ;
                     invmassbins[bb] = invmassbins_semilep[bb];
              }
         }
   std::cout << std::endl;
   std::cout << "-----------------------------------------" << std::endl;
   const int nunrollbins = nmsdbins*ninvmassbins;
   
    //Calling efficiency files 
  
   TFile *FEff; TH2F *h_AK4M_flv_b_eff; TH2F *h_AK4M_flv_c_eff; TH2F *h_AK4M_flv_l_eff;
   TH2F *h_YtagT_eff; TH2F *h_WtagT_eff;

   TFile *fSF_singlemuon;
   TH2F *h_singlemuon; TH2F *h_singlemuon_syst_up;  TH2F *h_singlemuon_syst_dn; 

   TFile *fSF_singleelectron;
   TH2F  *h_singleele; TH2F  *h_singleele_stat; TH2F  *h_singleele_syst;

   TFile *fSF_jet;
   TH2F  *h_jet; TH2F  *h_jet_stat; TH2F  *h_jet_syst;


   //https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopTrigger#Dilepton_triggers (miniAODv2)
   TFile *fSF_dilepton;
   TH2F  *h_trigmumu; TH2F  *h_trigemu; TH2F  *h_trigee;


   if(!isDATA)
   {
         if(isDL)
             {
		FEff = TFile::Open("Efficiency/Efficiency_TTTo2L2Nu_XtoYH.root" );
                std::cout << "We are taking the Efficiency/Efficiency_TTTo2L2Nu_XtoYH.root" << std::endl; 
             }
          else 
             {
                FEff = TFile::Open("Efficiency/Efficiency_TTToSemiLeptonic_XtoYH.root" );
                std::cout << "We are taking the Efficiency/Efficiency_TTToSemiLeptonic_XtoYH.root" << std::endl;
             }
 		h_AK4M_flv_b_eff = (TH2F*)FEff->Get("Efficiency_h_Ak4_b_flv_pass_M");
		h_AK4M_flv_c_eff = (TH2F*)FEff->Get("Efficiency_h_Ak4_c_flv_pass_M");
		h_AK4M_flv_l_eff = (TH2F*)FEff->Get("Efficiency_h_Ak4_l_flv_pass_M");
		h_YtagT_eff = (TH2F*)FEff->Get("Efficiency_h_Ak8_DeepTag_PNetMD_XbbvsQCD_pass_T");
		h_WtagT_eff = (TH2F*)FEff->Get("Efficiency_h_Ak8_DeepTag_PNetMD_WvsQCD_pass_T");
   }

   if(!isDATA)
   {
      if(!isDL)
      {
		//single muon trigger
                fSF_singlemuon = TFile::Open("Trigger_SF/OutFile-v20190510-Combined-Run2016BtoH_Run2017BtoF_Run2018AtoD-M120to10000.root");
		h_singlemuon      = (TH2F*)fSF_singlemuon->Get("SF_2018_var"); 
		h_singlemuon_syst_up = (TH2F*)fSF_singlemuon->Get("SF_2018_errorUpper");
		h_singlemuon_syst_dn = (TH2F*)fSF_singlemuon->Get("SF_2018_errorLower");

		//single electron trigger
		fSF_singleelectron  = TFile::Open("Trigger_SF/Trigger_electron_SF_Unc_MC_XtoYH_Pt_Eta.root");
		h_singleele    = (TH2F*)fSF_singleelectron->Get("h_electron_sf_trig");
                h_singleele_stat    = (TH2F*)fSF_singleelectron->Get("h_electron_sf_trig_stat_unc");
                h_singleele_syst    = (TH2F*)fSF_singleelectron->Get("h_electron_sf_trig_syst_unc");
               
                //ak8 jet trigger
                fSF_jet = TFile::Open("Trigger_SF/Trigger_LJet_SF_Unc_MC_XtoYH_Pt_Eta.root");
                h_jet          = (TH2F*)fSF_jet->Get("h_LJet_sf_trig");
                h_jet_stat          = (TH2F*)fSF_jet->Get("h_LJet_sf_trig_stat_unc");
                h_jet_syst          = (TH2F*)fSF_jet->Get("h_LJet_sf_trig_syst_unc");
      }
      else 
      {
        //Dilepton trigger
		fSF_dilepton = TFile::Open("Trigger_SF/DileptonTriggerSF_UL_miniAODv2/TriggerSF_2018_ULv2.root");
		h_trigmumu = (TH2F*)fSF_dilepton->Get("h2D_SF_mumu_lepABpt_FullError");
                h_trigee   = (TH2F*)fSF_dilepton->Get("h2D_SF_ee_lepABpt_FullError");
                h_trigemu  = (TH2F*)fSF_dilepton->Get("h2D_SF_emu_lepABpt_FullError");
      }
   }

   TString inputFile= input_path+argv[3];
   std::cout << inputFile << std::endl;
   TFile* final_file = TFile::Open("OUTPUTS/Histogram_"+inputFile, "RECREATE");  

   TFile *file = TFile::Open(inputFile);
   TTree *tree = (TTree*)file->Get("Tout");
      
   // read branches //
  read_branches(tree, isDL);
   
  // Declaration of histograms //
  
  final_file->cd();
 
  TH1D* h_nom;
  TH1D* h_nom_reg[nrgn];  
  TH1D* h_reg = new TH1D("h_reg_id","h_reg_id",12,0.0,12.0);

  TH1F* h_MET_pt_allreg;
  TH1F* h_Y_PNetMD_XbbvsQCD_allreg;
  TH1F* h_W_PNetMD_WvsQCD_opt1_allreg;
  TH1F* h_W_PNetMD_WvsQCD_opt2_allreg;
  TH1F* h_dR_lW_opt1_allreg;
  TH1F* h_dR_lW_opt2_allreg;
  TH1F* h_l1l2_mass_allreg;
  TH1F* h_l1l2_dR_allreg;
  TH1F* h_dphi_MET_l1l2_allreg;

  TH1F* h_MET_pt_ref[nrgn][nbcat][nWop][nlid];
  TH1F *h_Y_PNetMD_XbbvsQCD_ref[nrgn][nbcat][nWop][nlid];
  TH1F* h_W_msoftdrop_ref[nrgn][nbcat][nWop][nlid]; 
  
  TH1F* h_l1_pt[nrgn][nbcat][nWop][nlid];
  TH1F* h_l1_eta[nrgn][nbcat][nWop][nlid];
  TH1F* h_l1_minisoall[nrgn][nbcat][nWop][nlid];
  
  TH1F* h_l2_pt[nrgn][nbcat][nWop][nlid];
  TH1F* h_l2_eta[nrgn][nbcat][nWop][nlid];
  TH1F* h_l2_minisoall[nrgn][nbcat][nWop][nlid];
  
  TH1F* h_l1l2_mass[nrgn][nbcat][nWop][nlid];
  TH1F* h_l1l2_dR[nrgn][nbcat][nWop][nlid];
  TH1F* h_l1l2_deta[nrgn][nbcat][nWop][nlid];
  TH1F* h_l1l2_dphi[nrgn][nbcat][nWop][nlid];
  TH1F* h_dphi_MET_l1l2[nrgn][nbcat][nWop][nlid];
  
  TH1F* h_MET_pt[nrgn][nbcat][nWop][nlid];
  TH1F* h_MET_sig[nrgn][nbcat][nWop][nlid]; 
  
  TH1F* h_leadbjet_pt[nrgn][nbcat][nWop][nlid];  
  TH1F* h_leadbjet_btag_DeepFlav[nrgn][nbcat][nWop][nlid];  
  
  TH1F *h_Y_pt[nrgn][nbcat][nWop][nlid];
  
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

  TH1F* h_X_mass[nrgn][nbcat][nWop][nlid][ntop]; 
  TH1F* h_X_mass_xbin[nrgn][nbcat][nWop][nlid][ntop]; 

  TH1F* h_nbjets_other[nrgn][nbcat][nWop][nlid]; 
  TH1F* h_nbjets_outY[nrgn][nbcat][nWop][nlid]; 
  TH1F* h_nbjets_outY_L[nrgn][nbcat][nWop][nlid]; 
  TH1F* h_nbjets[nrgn][nbcat][nWop][nlid]; 
  TH1F* h_nbjets_L[nrgn][nbcat][nWop][nlid]; 
  
  TH1F* h_HTlep_pt[nrgn][nbcat][nWop][nlid][ntop]; 
  TH1F* h_ST[nrgn][nbcat][nWop][nlid][ntop];
  TH1F* h_ST_xbin[nrgn][nbcat][nWop][nlid][ntop]; // xbin ST distribution
  TH1F* h_MT[nrgn][nbcat][nWop][nlid][ntop]; 
 
  TH2F* h_X_Y_mass[nrgn][nbcat][nWop][nlid];   
  TH2F* h_X_Y_mass_xbin[nrgn][nbcat][nWop][nlid];   

  TH2F* h_X_Y_mass_sys[nrgn][nbcat][nWop][nlid][ntop][1+2*nsys];  
  TH2F* h_ST_Y_mass_sys[nrgn][nbcat][nWop][nlid][ntop][1+2*nsys];  
  
  
  TH1F *h_Y_msoftdrop[nrgn][nbcat][nWop][nlid][ntop];
  TH1F *h_Y_mass[nrgn][nbcat][nWop][nlid][ntop];
  TH1F *h_Y_msoftdrop_xbin[nrgn][nbcat][nWop][nlid][ntop];
  TH1F *h_Y_mass_xbin[nrgn][nbcat][nWop][nlid][ntop];
  TH1F *h_Y_PNetMD_XbbvsQCD[nrgn][nbcat][nWop][nlid][ntop];
  
  TH1F *h_Y_msoftdrop_sys[nrgn][nbcat][nWop][nlid][ntop][1+2*nsys];
  TH1F *h_Y_msoftdrop_xbin_sys[nrgn][nbcat][nWop][nlid][ntop][1+2*nsys];
  
  TH1F* h_X_mass_sys[nrgn][nbcat][nWop][nlid][ntop][1+2*nsys];
  TH1F* h_X_mass_xbin_sys[nrgn][nbcat][nWop][nlid][ntop][1+2*nsys];
  

  TH1F *h_for_limit_X_mass[nrgn][nbcat][nWop][nlid][ntop];
  TH1F *h_for_limit_X_mass_sys[nrgn][nbcat][nWop][nlid][ntop][1+2*nsys];

  TH1F *h_for_limit_X_mass_v2[nrgn][nbcat][nWop][nlid][ntop];
  TH1F *h_for_limit_X_mass_sys_v2[nrgn][nbcat][nWop][nlid][ntop][1+2*nsys];
  
  TH2F* h_HTlep_pt_Y_mass[nrgn][nbcat][nWop][nlid];   
  TH2F* h_ST_Y_mass[nrgn][nbcat][nWop][nlid];   
  
  TH1F* h_HTlep_pt_sys[nrgn][nbcat][nWop][nlid][ntop][1+2*nsys];
  TH1F* h_ST_sys[nrgn][nbcat][nWop][nlid][ntop][1+2*nsys];
  TH1F* h_ST_xbin_sys[nrgn][nbcat][nWop][nlid][ntop][1+2*nsys]; // xbin ST distribution 
  TH1F* h_MT_sys[nrgn][nbcat][nWop][nlid][ntop][1+2*nsys];
  //TH1F* h_ST_sys_v2[nrgn][nbcat][nWop][nlid][ntop][1+2*nsys];

  TH1F *h_for_limit_HTlep_pt[nrgn][nbcat][nWop][nlid][ntop];
  TH1F *h_for_limit_HTlep_pt_sys[nrgn][nbcat][nWop][nlid][ntop][1+2*nsys];

  TH1F *h_for_limit_ST[nrgn][nbcat][nWop][nlid][ntop];
  TH1F *h_for_limit_ST_sys[nrgn][nbcat][nWop][nlid][ntop][1+2*nsys];
  
  TH1F *h_for_limit_ST_v2[nrgn][nbcat][nWop][nlid][ntop];
  TH1F *h_for_limit_ST_sys_v2[nrgn][nbcat][nWop][nlid][ntop][1+2*nsys];
  
 
  // End of declaration //
    
  // Definition of histograms //
  
  h_nom = getHisto1D("h_nom","h_nom",10,0.0,10.0);
 
  for(int ij=0; ij<nrgn; ij++){
	char name[50];
	sprintf(name,"h_nom_%s",rgn[ij].Data());
	h_nom_reg[ij] = new TH1D(name,name,7,0.,7.);
	h_nom_reg[ij]->Sumw2();
  }
 
  h_MET_pt_allreg = getHisto1F("h_MET_pt_allreg","",40,0,1000);  
  h_Y_PNetMD_XbbvsQCD_allreg = getHisto1F("h_Y_PNetMD_XbbvsQCD_allreg","",100, 0.0, 1.0 );
  h_W_PNetMD_WvsQCD_opt1_allreg = getHisto1F("h_W_PNetMD_WvsQCD_opt1_allreg","",100, 0.0, 1.0 );
  h_W_PNetMD_WvsQCD_opt2_allreg = getHisto1F("h_W_PNetMD_WvsQCD_opt2_allreg","",100, 0.0, 1.0 );
  h_dR_lW_opt1_allreg = getHisto1F("h_dR_lW_opt1_allreg","",100, 0.0, 1.0 );
  h_dR_lW_opt2_allreg = getHisto1F("h_dR_lW_opt2_allreg","",100, 0.0, 1.0 );
  if(isDL){
	h_l1l2_mass_allreg = getHisto1F("h_l1l2_mass_allreg","",40,0,250);
	h_l1l2_dR_allreg = getHisto1F("h_l1l2_dR_allreg","",120,0,6);
	h_dphi_MET_l1l2_allreg = getHisto1F("h_dphi_MET_l1l2_allreg","",65,-M_PI,M_PI);
  }
  for (int ij=0 ; ij< nrgn ; ij++)
  {
          for(int jk=0; jk<nbcat; jk++)
          {
                  for(int kl=0; kl<nWop; kl++)
                  {
                         if(isDL && kl<0) { continue;  }
                         for(int lm=0; lm<nlid; lm++)
                          {
                                 for(int mn=0; mn<ntop; mn++)
                                 {
									 
									h_Y_msoftdrop[ij][jk][kl][lm][mn] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"Y_msoftdrop","",38,30,600); 
									h_X_mass[ij][jk][kl][lm][mn] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"X_mass","",40, 0.0, 4000.0);
									h_HTlep_pt[ij][jk][kl][lm][mn] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"HTlep_pt","",40, 0.0, 4000.0);
                                                                        h_ST[ij][jk][kl][lm][mn] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"ST","",40, 0.0, 4000.0);    
                                                                        h_MT[ij][jk][kl][lm][mn] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"MT","",25,0,300); 
                                                                        h_Y_msoftdrop_xbin[ij][jk][kl][lm][mn]   = get_histo_asymbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"Y_msoftdrop_xbin","",nmsdbins, msdbins);
                                                                        h_Y_mass_xbin[ij][jk][kl][lm][mn]        = get_histo_asymbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"Y_mass_xbin","",nmsdbins, msdbins);
                                                                        h_X_mass_xbin[ij][jk][kl][lm][mn]        = get_histo_asymbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"X_mass_xbin","",ninvmassbins,invmassbins);
                                                                        h_ST_xbin[ij][jk][kl][lm][mn]            = get_histo_asymbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"ST_xbin","",ninvmassbins,invmassbins);

									if(!isSignal){
									 
                                      h_Y_PNetMD_XbbvsQCD[ij][jk][kl][lm][mn] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"Y_PNetMD_XbbvsQCD","", 100, 0.0, 1.0 );
									  //h_Y_msoftdrop[ij][jk][kl][lm][mn] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"Y_msoftdrop","",38,30,600);
									  h_Y_mass[ij][jk][kl][lm][mn] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"Y_mass","",38,30,600);
                                      //h_X_mass[ij][jk][kl][lm][mn] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"X_mass","",40, 0.0, 4000.0);
                                      //h_HTlep_pt[ij][jk][kl][lm][mn] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"HTlep_pt","",40, 0.0, 4000.0);
                                      //h_ST[ij][jk][kl][lm][mn] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"ST","",40, 0.0, 4000.0);    

				      h_Y_msoftdrop_sys[ij][jk][kl][lm][mn][0]           = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"Y_msoftdrop","_nom",38,30,600);
                                      h_Y_msoftdrop_xbin_sys[ij][jk][kl][lm][mn][0]  = get_histo_asymbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"Y_msoftdrop_xbin","_nom",nmsdbins, msdbins);

                                      h_X_mass_sys[ij][jk][kl][lm][mn][0]                          = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"X_mass","_nom",40, 0.0, 4000.0);
                                      h_X_mass_xbin_sys[ij][jk][kl][lm][mn][0]         = get_histo_asymbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"X_mass_xbin","_nom",ninvmassbins,invmassbins);
                                      h_ST_xbin_sys[ij][jk][kl][lm][mn][0]             = get_histo_asymbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"ST_xbin","_nom",ninvmassbins,invmassbins);
                                      h_HTlep_pt_sys[ij][jk][kl][lm][mn][0]               = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"HTlep_pt","_nom",40, 0.0, 4000.0);
                                      h_ST_sys[ij][jk][kl][lm][mn][0]                     = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"ST","_nom",40, 0.0, 4000.0);
                                      h_MT_sys[ij][jk][kl][lm][mn][0]                     = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"MT","_nom",25,0,300);
									}
									
								      h_for_limit_X_mass[ij][jk][kl][lm][mn]                 = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_bin1_X_mass","",240,0,48000);
                                      h_for_limit_X_mass_v2[ij][jk][kl][lm][mn]              = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_bin2_X_mass","",nunrollbins,float(0.0),float(nunrollbins));

                                      h_for_limit_X_mass_sys[ij][jk][kl][lm][mn][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_bin1_X_mass","_nom",240,0,48000);
                                      h_for_limit_X_mass_sys_v2[ij][jk][kl][lm][mn][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_bin2_X_mass","_nom",nunrollbins,float(0.0),float(nunrollbins));

                                      h_for_limit_HTlep_pt[ij][jk][kl][lm][mn]                 = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_HTlep_pt","",240,0,48000);
                                      h_for_limit_ST[ij][jk][kl][lm][mn]                 = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_ST","",240,0,48000);
				      h_for_limit_ST_v2[ij][jk][kl][lm][mn]              = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_ST_v2","",nunrollbins,float(0.0),float(nunrollbins));

                                      h_for_limit_HTlep_pt_sys[ij][jk][kl][lm][mn][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_HTlep_pt","_nom",240,0,48000);
                                      h_for_limit_ST_sys[ij][jk][kl][lm][mn][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_ST","_nom",240,0,48000);                                         
                                      h_for_limit_ST_sys_v2[ij][jk][kl][lm][mn][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_ST_v2","_nom",nunrollbins,float(0.0),float(nunrollbins));                                         
                                          
                                      h_X_Y_mass_sys[ij][jk][kl][lm][mn][0]              = new TH2F("h_Y_"+Ytype[y_wp]+"_W_"+Wtype[w_wp]+"_X_Y_mass_"+Wops[kl]+"_"+rgn[ij]+bcats[jk]+lepids[lm]+tops[mn]+"_nom", "",40, 0.0, 4000.0, 40, 0.0, 600.0);
                                      h_ST_Y_mass_sys[ij][jk][kl][lm][mn][0]             = new TH2F("h_Y_"+Ytype[y_wp]+"_W_"+Wtype[w_wp]+"_ST_Y_mass_"+Wops[kl]+"_"+rgn[ij]+bcats[jk]+lepids[lm]+tops[mn]+"_nom","", 40, 0.0, 4000.0, 40, 0.0, 600.0);
                                        
                                      //Systematic uncertainties
									  for(int isys=0; isys<nsys; isys++){
                                        
                                        char name[100];
                                        //up systematics
                                        sprintf(name,"_%s_up",sysnames[isys].Data());
                                        if(!isSignal){
											h_Y_msoftdrop_sys[ij][jk][kl][lm][mn][2*(isys+1)-1]          = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"Y_msoftdrop",name,38,30,600);
											h_X_mass_sys[ij][jk][kl][lm][mn][2*(isys+1)-1]               = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"X_mass",name,40, 0.0, 4000.0);
										}
                                        h_for_limit_X_mass_sys[ij][jk][kl][lm][mn][2*(isys+1)-1]     = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_bin1_X_mass",name,240,0,48000);
										h_for_limit_X_mass_sys_v2[ij][jk][kl][lm][mn][2*(isys+1)-1]  = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_bin2_X_mass",name,nunrollbins,float(0.0),float(nunrollbins));
                                        if(!isSignal){
											h_Y_msoftdrop_xbin_sys[ij][jk][kl][lm][mn][2*(isys+1)-1]     = get_histo_asymbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"Y_msoftdrop_xbin",name,nmsdbins, msdbins);
											h_X_mass_xbin_sys[ij][jk][kl][lm][mn][2*(isys+1)-1]          = get_histo_asymbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"X_mass_xbin",name,ninvmassbins,invmassbins);
											h_ST_xbin_sys[ij][jk][kl][lm][mn][2*(isys+1)-1]              = get_histo_asymbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"ST_xbin",name,ninvmassbins,invmassbins);
											h_HTlep_pt_sys[ij][jk][kl][lm][mn][2*(isys+1)-1]             = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"HTlep_pt",name,40, 0.0, 4000.0);
											h_ST_sys[ij][jk][kl][lm][mn][2*(isys+1)-1]                   = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"ST",name,40, 0.0, 4000.0);
                                                                                        h_MT_sys[ij][jk][kl][lm][mn][2*(isys+1)-1]                   = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"MT",name,25,0,300);
										}
                                        h_for_limit_HTlep_pt_sys[ij][jk][kl][lm][mn][2*(isys+1)-1]   = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_HTlep_pt",name,240,0,48000);
                                        h_for_limit_ST_sys[ij][jk][kl][lm][mn][2*(isys+1)-1]         = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_ST",name,240,0,48000);
										h_for_limit_ST_sys_v2[ij][jk][kl][lm][mn][2*(isys+1)-1]         = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_ST_v2",name,nunrollbins,float(0.0),float(nunrollbins));
										
										h_X_Y_mass_sys[ij][jk][kl][lm][mn][2*(isys+1)-1]             = new TH2F("h_Y_"+Ytype[y_wp]+"_W_"+Wtype[w_wp]+"_X_Y_mass_"+Wops[kl]+"_"+rgn[ij]+bcats[jk]+lepids[lm]+tops[mn]+name, "", 40, 0.0, 4000.0, 40, 0.0, 600.0);
                                        h_ST_Y_mass_sys[ij][jk][kl][lm][mn][2*(isys+1)-1]            = new TH2F("h_Y_"+Ytype[y_wp]+"_W_"+Wtype[w_wp]+"_ST_Y_mass_"+Wops[kl]+"_"+rgn[ij]+bcats[jk]+lepids[lm]+tops[mn]+name, "", 40, 0.0, 4000.0, 40, 0.0, 600.0);
                                        
										//dn systematics
                                        sprintf(name,"_%s_dn",sysnames[isys].Data());
                                        if(!isSignal){
											h_Y_msoftdrop_sys[ij][jk][kl][lm][mn][2*(isys+1)]            = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"Y_msoftdrop",name,38,30,600);
											h_X_mass_sys[ij][jk][kl][lm][mn][2*(isys+1)]                 = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"X_mass",name,40, 0.0, 4000.0);
										}
										h_for_limit_X_mass_sys[ij][jk][kl][lm][mn][2*(isys+1)]       = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_bin1_X_mass",name,240,0,48000);
										h_for_limit_X_mass_sys_v2[ij][jk][kl][lm][mn][2*(isys+1)]    = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_bin2_X_mass",name,nunrollbins,float(0.0),float(nunrollbins));
                                        if(!isSignal){
											h_Y_msoftdrop_xbin_sys[ij][jk][kl][lm][mn][2*(isys+1)]       = get_histo_asymbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"Y_msoftdrop_xbin",name,nmsdbins, msdbins);
											h_X_mass_xbin_sys[ij][jk][kl][lm][mn][2*(isys+1)]            = get_histo_asymbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"X_mass_xbin",name,ninvmassbins,invmassbins);
										        h_ST_xbin_sys[ij][jk][kl][lm][mn][2*(isys+1)]                = get_histo_asymbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"ST_xbin",name,ninvmassbins,invmassbins);
											h_HTlep_pt_sys[ij][jk][kl][lm][mn][2*(isys+1)]               = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"HTlep_pt",name,40, 0.0, 4000.0);
											h_ST_sys[ij][jk][kl][lm][mn][2*(isys+1)]                     = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"ST",name,40, 0.0, 4000.0);
                                                                                        h_MT_sys[ij][jk][kl][lm][mn][2*(isys+1)]                     = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"MT",name,25,0,300);
                                                                                 }
                                        h_for_limit_HTlep_pt_sys[ij][jk][kl][lm][mn][2*(isys+1)]     = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_HTlep_pt",name,240,0,48000);
                                        h_for_limit_ST_sys[ij][jk][kl][lm][mn][2*(isys+1)]           = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_ST",name,240,0,48000);
										h_for_limit_ST_sys_v2[ij][jk][kl][lm][mn][2*(isys+1)]        = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],tops[mn],"unrolled_ST_v2",name,nunrollbins,float(0.0),float(nunrollbins));
			
										h_X_Y_mass_sys[ij][jk][kl][lm][mn][2*(isys+1)]               = new TH2F("h_Y_"+Ytype[y_wp]+"_W_"+Wtype[w_wp]+"_X_Y_mass_"+Wops[kl]+"_"+rgn[ij]+bcats[jk]+lepids[lm]+tops[mn]+"_"+name, "", 40, 0.0, 4000.0, 40, 0.0, 600.0);
                                        h_ST_Y_mass_sys[ij][jk][kl][lm][mn][2*(isys+1)]              = new TH2F("h_Y_"+Ytype[y_wp]+"_W_"+Wtype[w_wp]+"_ST_Y_mass_"+Wops[kl]+"_"+rgn[ij]+bcats[jk]+lepids[lm]+tops[mn]+"_"+name, "", 40, 0.0, 4000.0, 40, 0.0, 600.0);
                                        
                                     } //sys_unc               
				 }//mn
			  }//lm
		  }//kl
	  }//jk
  }//ij
  for (int ij=0 ; ij< nrgn ; ij++)
  {
	  for(int jk=0; jk<nbcat; jk++)
	  {
		  for(int kl=0; kl<nWop; kl++)
		  {
			  if(isDL && kl<0) { continue;  } // no need to loop over W candidates for dileptonic channel
			  
			  for(int lm=0; lm<nlid; lm++)
			  {
				  
				//ref hist //
				
				if(!isSignal){
				
				h_MET_pt_ref[ij][jk][kl][lm] 			= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"MET_pt_ref","",40,0,1000);
				h_Y_PNetMD_XbbvsQCD_ref[ij][jk][kl][lm] = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_PNetMD_XbbvsQCD_ref","", 100, 0.0, 1.0 );
				h_W_msoftdrop_ref[ij][jk][kl][lm] 		 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"W_msoftdrop_ref","",50,0,350);
				
				}
				
				// all hist//
                
                if(!isSignal){
                             
				h_l1_pt[ij][jk][kl][lm] 				= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"l1_pt","",40,0,1000);
				h_l1_eta[ij][jk][kl][lm]				= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"l1_eta","",50,-2.5,2.5);
				h_l1_minisoall[ij][jk][kl][lm]		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"l1_minisoall","",150,0,1.5);
				
				h_l2_pt[ij][jk][kl][lm] 				= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"l2_pt","",40,0,1000);
				h_l2_eta[ij][jk][kl][lm]				= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"l2_eta","",50,-2.5,2.5);
				h_l2_minisoall[ij][jk][kl][lm]		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"l2_minisoall","",150,0,1.5);
				
				h_l1l2_mass[ij][jk][kl][lm] 				= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"l1l2_mass","",40,0,250);
				h_l1l2_dR[ij][jk][kl][lm] 				= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"l1l2_dR","",120,0,6);
				h_l1l2_deta[ij][jk][kl][lm] 				= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"l1l2_deta","",50,-5,5);
				h_l1l2_dphi[ij][jk][kl][lm] 				= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"l1l2_dphi","",65,-M_PI,M_PI);
				h_dphi_MET_l1l2[ij][jk][kl][lm] 				= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"dphi_MET_l1l2","",65,-M_PI,M_PI);
	       	
				h_MET_pt[ij][jk][kl][lm] 			= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"MET_pt","",40,0,1000);
				h_MET_sig[ij][jk][kl][lm] 			= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"MET_sig","",50,0,300);				
			
				}
				//h_leadbjet_pt[ij][jk][kl][lm] 		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"LeadBJet_pt","",40,0,1000);
				//h_leadbjet_btag_DeepFlav[ij][jk][kl][lm] = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"LeadBJet_btag_DeepFlav","",50,0,1);
				
				h_Y_pt[ij][jk][kl][lm] 				= get_histo_asymbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_pt","",nptbins, ptedges);
				//h_Y_msoftdrop[ij][jk][kl][lm] 		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_msoftdrop","",38,30,600);
				//h_Y_msoftdrop_xbin[ij][jk][kl][lm]  = get_histo_asymbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_msoftdrop_xbin","",nmsdbins, msdbins);
				//h_Y_PNetMD_XbbvsQCD[ij][jk][kl][lm] = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_PNetMD_XbbvsQCD","", 100, 0.0, 1.0 );
				
				if(!isSignal){
				
				h_Y_PNetMD_WvsQCD[ij][jk][kl][lm] 	= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_PNetMD_WvsQCD","", 100, 0.0, 1.0 );
				h_Y_PNet_TvsQCD[ij][jk][kl][lm] 	= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_PNet_TvsQCD","", 100, 0.0, 1.0 );
				h_Y_sub1_mass[ij][jk][kl][lm] 		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_sub1_mass","", 40, 0.0, 300 );
				h_Y_sub2_mass[ij][jk][kl][lm] 		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_sub2_mass","", 40, 0.0, 300 );
				h_Y_sub1_btag[ij][jk][kl][lm] 		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_sub1_btag","", 50, 0.0, 1.0 );
				h_Y_sub2_btag[ij][jk][kl][lm] 		= get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"Y_sub2_btag","", 50, 0.0, 1.0 );
				
				}
				
				h_W_msoftdrop[ij][jk][kl][lm] 		 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"W_msoftdrop","",50,0,350);
				
				if(!isSignal){
				
				h_W_pt[ij][jk][kl][lm]	 			 = get_histo_asymbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"W_pt","",nptbins, ptedges);
				//h_W_msoftdrop[ij][jk][kl][lm] 		 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"W_msoftdrop","",50,0,350);
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
				
				}
				
				h_nbjets_other[ij][jk][kl][lm]		 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"nbjets_other","",5 ,0.0, 5.0 );
				h_nbjets_outY[ij][jk][kl][lm]		 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"nbjets_outY","",5 ,0.0, 5.0 );
				h_nbjets_outY_L[ij][jk][kl][lm]		 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"nbjets_outY_L","",5 ,0.0, 5.0 );
				h_nbjets[ij][jk][kl][lm]		 	 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"nbjets","",5 ,0.0, 5.0 );
				h_nbjets_L[ij][jk][kl][lm]		 	 = get_histo_symbin(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],lepids[lm],"nbjets_L","",5 ,0.0, 5.0 );
               
				if(!isSignal){
               
					h_X_Y_mass[ij][jk][kl][lm] 		    = new TH2F("h_Y_"+Ytype[y_wp]+"_W_"+Wtype[w_wp]+"_X_Y_mass_"+Wops[kl]+"_"+rgn[ij]+bcats[jk]+lepids[lm], "", 40, 0.0, 4000.0, 40, 0.0, 600.0);
					h_X_Y_mass_xbin[ij][jk][kl][lm] 	= new TH2F("h_Y_"+Ytype[y_wp]+"_W_"+Wtype[w_wp]+"_X_Y_mass_"+Wops[kl]+"_"+rgn[ij]+bcats[jk]+lepids[lm]+"_xbin", "", ninvmassbins, invmassbins, nmsdbins, msdbins);
							
					h_HTlep_pt_Y_mass[ij][jk][kl][lm] 	= new TH2F("h_Y_"+Ytype[y_wp]+"_W_"+Wtype[w_wp]+"_HTlep_pt_Y_mass_"+Wops[kl]+"_"+rgn[ij]+bcats[jk]+lepids[lm], "", 40, 0.0, 4000.0, 40, 0.0, 600.0);
					h_ST_Y_mass[ij][jk][kl][lm] 		= new TH2F("h_Y_"+Ytype[y_wp]+"_W_"+Wtype[w_wp]+"_ST_Y_mass_"+Wops[kl]+"_"+rgn[ij]+bcats[jk]+lepids[lm], "", 40, 0.0, 4000.0, 40, 0.0, 600.0);
	
				}
	
				}
	         }
		}
    }
    
   // end of histogram defitions // 
    	 
   file->cd();
 
   Long64_t nn = tree->GetEntries();
   
   //// Event loop ////
   std::cout << nn << std::endl;  
   for(Long64_t jentry =0; jentry < nn ; jentry++)
   {
	      
	tree->GetEntry(jentry);
	if( jentry % 10000 == 0) { std::cout <<jentry<<" events processed" << std::endl;}
   
	// read number of JES+JER uncs. //
	njecmax = (*Y_JESup_split).size();
	// Condition to avoid double counting in data //
        bool mu_trig = Muon_trig_pass;
        bool el_trig = Electron_trig_pass;	
		bool jet_trig = hlt_AK8PFJet500 || hlt_PFJet500;
        bool emu_trig = hlt_Mu37_Ele27_CaloIdL_MW || hlt_Mu27_Ele37_CaloIdL_MW;
		if(isDATA){
			if (string(inputFile.Data()).find("SingleMuon")!=string::npos)
			{                        
				if(!mu_trig) continue;
			}
            else if (string(inputFile.Data()).find("DoubleMuon")!=string::npos)
            {
                if(!mu_trig) continue;
            }
		else if (string(inputFile.Data()).find("EGamma")!=string::npos)
		{
			if(mu_trig || !el_trig) continue;
		}
		else if (string(inputFile.Data()).find("JetHT")!=string::npos)
		{
			if(mu_trig || el_trig || !jet_trig) continue;
		}
        else if (string(inputFile.Data()).find("MuonEG")!=string::npos)
        {
            if(mu_trig || el_trig || !emu_trig) continue;
        }
		else{
			continue;
		}
	}
        
    //Trigger_pt Scale factor determination and the corresponding unc
	double SF_Trig, SF_Trig_stat, SF_Trig_syst;
	SF_Trig = 1; SF_Trig_stat = 1; SF_Trig_syst = 1;
    if(!isDATA) {
		if(!isDL) {
	        if(mu_trig && abs(l1_pdgId)==13) {
		      if(l1_pt > 50)
		          {
				int etabin = std::max(1, std::min(h_singlemuon->GetNbinsY(), h_singlemuon->GetYaxis()->FindBin(fabs(l1_eta))));
                                int ptbin  = std::max(1, std::min(h_singlemuon->GetNbinsX(), h_singlemuon->GetXaxis()->FindBin(l1_pt)));
                                SF_Trig = TMath::Max(float(1.e-3),float(h_singlemuon->GetBinContent(ptbin,etabin))) ;
		                SF_Trig_stat = 0.0 ;	
		                SF_Trig_syst = std::max(TMath::Max(float(1.e-3),float(h_singlemuon_syst_up->GetBinContent(ptbin,etabin))), TMath::Max(float(1.e-3),float(h_singlemuon_syst_dn->GetBinContent(ptbin,etabin)))) ;
				//std::cout << fabs(l1_pdgId) << "\t" << l1_pt <<"\t" << l1_eta << "\t" << ptbin << "\t" << etabin << "\t" << SF_Trig << "\t" << SF_Trig_stat << "\t" << SF_Trig_syst << std::endl;
                          } // l1_pt threshold
                     } //muon trigger
	        else if ( !mu_trig  && (hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || hlt_Ele115_CaloIdVT_GsfTrkIdT || hlt_Ele40_WPTight_Gsf || hlt_Ele28_eta2p1_WPTight_Gsf_HT150 ||  hlt_Ele32_WPTight_Gsf) &&  abs(l1_pdgId)==11) {
			  if(l1_pt > 32 )
			  {
				int etabin = std::max(1, std::min(h_singleele->GetNbinsY(), h_singleele->GetYaxis()->FindBin(fabs(l1_eta))));
                                int ptbin  = std::max(1, std::min(h_singleele->GetNbinsX(), h_singleele->GetXaxis()->FindBin(l1_pt)));
                                SF_Trig      = TMath::Max(float(1.e-3),float(h_singleele->GetBinContent(ptbin, etabin))) ;
				SF_Trig_stat = TMath::Max(float(1.e-3),float(h_singleele_stat->GetBinContent(ptbin, etabin))) ;
				SF_Trig_syst = TMath::Max(float(1.e-3),float(h_singleele_syst->GetBinContent(ptbin, etabin))) ;
                        //std::cout << fabs(l1_pdgId) << "\t" << l1_pt <<"\t" << l1_eta << "\t" << ptbin << "\t" << etabin << "\t" << SF_Trig << "\t" << SF_Trig_stat << "\t" << SF_Trig_syst << std::endl;
                          } //l1_pt threshold
		     } //electron trigger
                else if ( !mu_trig || !el_trig && jet_trig ){
                           if( PFJetAK8_pt[0] > 550  )
                           {
                                int etabin = std::max(1, std::min(h_jet->GetNbinsY(), h_jet->GetYaxis()->FindBin(fabs(PFJetAK8_eta[0]))));
                                int ptbin  = std::max(1, std::min(h_jet->GetNbinsX(), h_jet->GetXaxis()->FindBin(PFJetAK8_pt[0])));
                                SF_Trig      = TMath::Max(float(1.e-3),float(h_jet->GetBinContent( ptbin, etabin))) ;
                                SF_Trig_stat = TMath::Max(float(1.e-3),float(h_jet_stat->GetBinContent(ptbin, etabin))) ;
                                SF_Trig_syst = TMath::Max(float(1.e-3),float(h_jet_syst->GetBinContent(ptbin, etabin))) ;
                           }//jet pt thresholds
                     } //jet trig
        } //Semi-lepton
    else{
		if( mu_trig && abs(l1_pdgId)==13 && abs(l2_pdgId)==13)
		{
			if(l1_pt > 37 && l2_pt > 27  )
		        {
				int pt1bin = std::max(1, std::min(h_trigmumu->GetNbinsX(), h_trigmumu->GetXaxis()->FindBin(l1_pt)));
				int pt2bin = std::max(1, std::min(h_trigmumu->GetNbinsY(), h_trigmumu->GetYaxis()->FindBin(l2_pt)));
				SF_Trig         = TMath::Max(float(1.e-3),float(h_trigmumu->GetBinContent(pt1bin, pt2bin))) ;	
				SF_Trig_stat    = 0.0;
				SF_Trig_syst    = TMath::Max(float(1.e-3),float(h_trigmumu->GetBinError(pt1bin, pt2bin))) ;
				//std::cout << fabs(l1_pdgId) << "\t" << fabs(l2_pdgId) << "\t" << l1_pt <<"\t" << l2_pt << "\t" << pt1bin << "\t" << pt2bin << "\t" << SF_Trig << std::endl;
                        }
   		 } //mu-mu trigger
         else if( !mu_trig && el_trig && abs(l1_pdgId)==11 && abs(l2_pdgId)==11)
         {
		        if(l1_pt > 25 && l2_pt > 25 )
		        {
				int pt1bin = std::max(1, std::min(h_trigmumu->GetNbinsX(), h_trigmumu->GetXaxis()->FindBin(l1_pt)));
                                int pt2bin = std::max(1, std::min(h_trigmumu->GetNbinsY(), h_trigmumu->GetYaxis()->FindBin(l2_pt)));
                                SF_Trig         = TMath::Max(float(1.e-3),float(h_trigmumu->GetBinContent(pt1bin, pt2bin))) ;
				SF_Trig_stat    = 0.0 ;
				SF_Trig_syst    = TMath::Max(float(1.e-3),float(h_trigmumu->GetBinError(pt1bin, pt2bin))) ;
                                //std::cout << fabs(l1_pdgId) << "\t" << fabs(l2_pdgId) << "\t" << l1_pt <<"\t" << l2_pt << "\t" << pt1bin << "\t" << pt2bin << "\t" << SF_Trig << std::endl;
                       }
   		 } // e-e trigger
         else if(!mu_trig && !el_trig &&  emu_trig && ((abs(l1_pdgId)==11 && abs(l2_pdgId)==13) || (abs(l1_pdgId)==13 && abs(l2_pdgId)==11)) )
		 {
			if(l1_pt > 37 && l2_pt > 27)
            {
				int ptelebin, ptmubin;
                if(abs(l1_pdgId)==11 && abs(l2_pdgId)==13) 
				{
					ptelebin = std::max(1, std::min(h_trigemu->GetNbinsX(), h_trigemu->GetXaxis()->FindBin(l1_pt)));
					ptmubin  = std::max(1, std::min(h_trigemu->GetNbinsY(), h_trigemu->GetYaxis()->FindBin(l2_pt)));
				}
                else if(abs(l1_pdgId)==13 && abs(l2_pdgId)==11)
				{
					ptmubin  = std::max(1, std::min(h_trigemu->GetNbinsX(), h_trigemu->GetXaxis()->FindBin(l1_pt)));
                                        ptelebin = std::max(1, std::min(h_trigemu->GetNbinsY(), h_trigemu->GetYaxis()->FindBin(l2_pt)));
				}
				SF_Trig         = TMath::Max(float(1.e-3),float(h_trigemu->GetBinContent(ptelebin, ptmubin))) ;
				SF_Trig_stat    = 0.0 ;
				SF_Trig_syst    = TMath::Max(float(1.e-3),float(h_trigemu->GetBinError(ptelebin, ptmubin))) ;
		    }
		 }  //e-mu trigger
	    } //di-lepton
	} //Only MC	
       
    //std::cout << fabs(l1_pdgId) << "\t" << l1_pt <<"\t" << l1_eta << "\t"  << SF_Trig << "\t" << SF_Trig_stat << "\t" << SF_Trig_syst << std::endl;

    double SF_Trig_1_up = SF_Trig + abs(SF_Trig_stat);  
    double SF_Trig_1_dn = SF_Trig - abs(SF_Trig_stat);

    double SF_Trig_2_up = SF_Trig + abs(SF_Trig_syst);
    double SF_Trig_2_dn = SF_Trig - abs(SF_Trig_syst);

    //cout << SF_Trig << "\t" << SF_Trig_stat << "\t" << SF_Trig_syst << "\t" << SF_Trig_1_up << "\t" << SF_Trig_1_dn << "\t" << SF_Trig_2_up << "\t" << SF_Trig_2_dn << endl;

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
		if(delta2R(JetAK4_eta[jj],JetAK4_phi[jj],Y_eta,Y_phi)>1.2){	
					
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
					b_SF *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF[jj]*TMath::Max(float(1.e-3),float(h_AK4M_flv_b_eff->GetBinContent(ptbin, etabin)) ))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4M_flv_b_eff->GetBinContent(ptbin, etabin))) ));
                                        b_SF_up *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF_up[jj]*TMath::Max(float(1.e-3),float(h_AK4M_flv_b_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4M_flv_b_eff->GetBinContent(ptbin, etabin))) ));
					b_SF_dn *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF_dn[jj]*TMath::Max(float(1.e-3),float(h_AK4M_flv_b_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4M_flv_b_eff->GetBinContent(ptbin, etabin))) ));
                }
                else if ( abs(JetAK4_hadronflav[jj]) == 4 ) {
					b_SF *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF[jj]*TMath::Max(float(1.e-3),float(h_AK4M_flv_c_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4M_flv_c_eff->GetBinContent(ptbin, etabin))) ));
                                        b_SF_up *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF_up[jj]*TMath::Max(float(1.e-3),float(h_AK4M_flv_c_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4M_flv_c_eff->GetBinContent(ptbin, etabin))) ));
					b_SF_dn *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF_dn[jj]*TMath::Max(float(1.e-3),float(h_AK4M_flv_c_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4M_flv_c_eff->GetBinContent(ptbin, etabin))) ));
                }
				else {
                    b_SF *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF[jj]*TMath::Max(float(1.e-3),float(h_AK4M_flv_l_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4M_flv_l_eff->GetBinContent(ptbin, etabin)))));
                    b_SF_up *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF_up[jj]*TMath::Max(float(1.e-3),float(h_AK4M_flv_l_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4M_flv_l_eff->GetBinContent(ptbin, etabin)))));
                    b_SF_dn *= std::max(1.e-6,(1.-JetAK4_btag_DeepFlav_SF_dn[jj]*TMath::Max(float(1.e-3),float(h_AK4M_flv_l_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4M_flv_l_eff->GetBinContent(ptbin, etabin)))));
                }
            }
		}
	}//jj

        //Y-cand gen bb matching
	TLorentzVector Y_cand;
	Y_cand.SetPtEtaPhiM(Y_pt,Y_eta,Y_phi,Y_mass);
	Bool_t Y_cand_gen_bb_matched = false;
        for(int pp = 0 ; pp  < (nGenBPart-1); pp++)
        {
                for(int qq = pp+1; qq < nGenBPart; qq++)
                         {
                                 TLorentzVector gen_b1, gen_b2;
                                 gen_b1.SetPtEtaPhiM(GenBPart_pt[pp],GenBPart_eta[pp],GenBPart_phi[pp],GenBPart_mass[pp]);
                                 gen_b2.SetPtEtaPhiM(GenBPart_pt[qq],GenBPart_eta[qq],GenBPart_phi[qq],GenBPart_mass[qq]);
                                 Float_t DR_Yb1, DR_Yb2;
                                 DR_Yb1 = Y_cand.DeltaR(gen_b1);
                                 DR_Yb2 = Y_cand.DeltaR(gen_b2);
                                 if(DR_Yb1 < 0.8 && DR_Yb2 < 0.8 && GenBPart_pdgId[pp]*GenBPart_pdgId[qq] < 0 && GenBPart_mompdgId[pp] == GenBPart_mompdgId[qq] )
                                 {
                                       Y_cand_gen_bb_matched = true;
                                       break;
                                 }
                         } //qq
        }//pp               

	int Y_pt_bin = getbinid(Y_pt,PNbb_SF_nptbins,PNbb_SF_ptbins);
        int Y_etabin_YtagT = std::max(1, std::min(h_YtagT_eff->GetNbinsY(), h_YtagT_eff->GetYaxis()->FindBin(fabs(Y_eta))));
        int Y_ptbin_YtagT  = std::max(1, std::min(h_YtagT_eff->GetNbinsX(), h_YtagT_eff->GetXaxis()->FindBin(Y_pt)));       
	if(Y_pt_bin>=0 && Y_pt_bin<PNbb_SF_nptbins) {
                if(Y_cand_gen_bb_matched){
                       if(Flag_Y_bb_pass_T)
                        {
			  bb_SF = PNbb_SF_HP[Y_pt_bin]; 
			  bb_SF_up = PNbb_SF_HP_up[Y_pt_bin]; 
			  bb_SF_dn = PNbb_SF_HP_dn[Y_pt_bin]; 
                        }
                       else 
                        {
                          bb_SF = std::max(1.e-6,(1.0 - PNbb_SF_HP[Y_pt_bin]*TMath::Max(float(1.e-3),float(h_YtagT_eff->GetBinContent(Y_ptbin_YtagT, Y_etabin_YtagT))) )/ (1.0 - TMath::Max(float(1.e-3),float(h_YtagT_eff->GetBinContent(Y_ptbin_YtagT, Y_etabin_YtagT)))) );
                          bb_SF_up = std::max(1.e-6,(1.0 - PNbb_SF_HP_up[Y_pt_bin]*TMath::Max(float(1.e-3),float(h_YtagT_eff->GetBinContent(Y_ptbin_YtagT, Y_etabin_YtagT))) )/ (1.0 - TMath::Max(float(1.e-3),float(h_YtagT_eff->GetBinContent(Y_ptbin_YtagT, Y_etabin_YtagT)))) );
                          bb_SF_dn = std::max(1.e-6,(1.0 - PNbb_SF_HP_dn[Y_pt_bin]*TMath::Max(float(1.e-3),float(h_YtagT_eff->GetBinContent(Y_ptbin_YtagT, Y_etabin_YtagT))) )/ (1.0 - TMath::Max(float(1.e-3),float(h_YtagT_eff->GetBinContent(Y_ptbin_YtagT, Y_etabin_YtagT)))) );
                          //std::cout << bb_SF << "	" << bb_SF_up << "	" << bb_SF_dn << " " << Y_ptbin_YtagT << " " << Y_etabin_YtagT << " " <<  TMath::Max(float(1.e-3),float(h_YtagT_eff->GetBinContent(Y_ptbin_YtagT, Y_etabin_YtagT)))  <<  std::endl; 
                        }
		}
	}
       
	int W_pt_bin = getbinid(W_pt_opt2,PNW_SF_nptbins,PNW_SF_ptbins);
	if(W_pt_bin>=0 && W_pt_bin<PNW_SF_nptbins) { 
		if(W_label_W_qq_opt2||W_label_W_cq_opt2){// use only GEN-matched W for applying SFs 
			if(Flag_H_W_pass_T_opt2){
				W_SF = PNW_SF_T[W_pt_bin]; 
				W_SF_up = PNW_SF_T_up[W_pt_bin]; 
				W_SF_dn = PNW_SF_T_dn[W_pt_bin]; 
			}
		}
	}
   
	int Y_top_pt_bin = getbinid(Y_pt,PNTop_SF_nptbins,PNTop_SF_ptbins);
	if(Y_top_pt_bin>=0 && Y_top_pt_bin<PNTop_SF_nptbins) { 
		if(Y_label_Top_bqq||Y_label_Top_bcq){
			if(Y_DeepTag_PNet_TvsQCD>=PN_Top_med){
				Top_SF = PNTop_SF_M[Y_top_pt_bin]; 
				Top_SF_up = PNTop_SF_M_up[Y_top_pt_bin]; 
				Top_SF_dn = PNTop_SF_M_dn[Y_top_pt_bin]; 
			}
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
   h_nom->Fill(7.0,SF_Trig);
   h_nom->Fill(8.0,SF_Trig*b_SF*W_SF*bb_SF); 
   // end of tagger SF //
   
   float weight_nom;
   if(isDATA) {weight_nom = 1.0;}
   else 
   {
   if (isSignal)
   {
     weight_nom = 1.0; //std::cout << "signal event" << std::endl;
   }
   else
   {
     weight_nom = Generator_weight; //std::cout << "Background event" << std::endl;
   }
     weight_nom *= prefiringweight;
     weight_nom *= puWeight;
     weight_nom *= leptonsf_weight;
     //weight_nom *= b_SF;           (apply b tagging scale factor later only while using b tagging based selection condition)
     weight_nom *= bb_SF;
     weight_nom *= SF_Trig; //Trigger scale factors 
     if(!isDL){
		weight_nom *= W_SF;
	 }
   }
   
   // Selections so that regions where ParticleNet bb tagger is not trianed or SFs are not derived are not used //

   if(Y_msoftdrop < msd_cut) continue;
   
   // inclusive histograms //
   
   h_MET_pt_allreg->Fill(MET_pt,weight_nom); 
   h_Y_PNetMD_XbbvsQCD_allreg->Fill(Y_DeepTag_PNetMD_XbbvsQCD,weight_nom); 
   if(isDL){
	h_l1l2_mass_allreg->Fill(l1l2_mass,weight_nom); 
	h_l1l2_dR_allreg->Fill(l1l2_dR,weight_nom); 
	h_dphi_MET_l1l2_allreg->Fill(dphi_MET_l1l2,weight_nom); 
   }
   else{
	h_W_PNetMD_WvsQCD_opt1_allreg->Fill(W_DeepTag_PNetMD_WvsQCD_opt1,weight_nom);
	h_W_PNetMD_WvsQCD_opt2_allreg->Fill(W_DeepTag_PNetMD_WvsQCD_opt2,weight_nom);
	h_dR_lW_opt1_allreg->Fill(dR_lW_opt1,weight_nom);
	h_dR_lW_opt2_allreg->Fill(dR_lW_opt2,weight_nom);
   }
 
   // looping over choices of W candidate
   
   for(int kl=0; kl<nWop; kl++)
   {
   
   // only run once for dileptonic channel
   if(isDL && kl>0) break; 
   
   // Opposite-charge condition for dileptonic channel
   //std::cout << l1_pdgId << "\t" << l2_pdgId << "\t" << l1_pdgId*l2_pdgId << std::endl;
   if(isDL && l1_pdgId*l2_pdgId>0) continue;
   
	// Defining booleans for signal & control regions //
	
	vector<bool> reg_tags;
   
	bool isSR1(false);
	bool isSR2(false);
	bool isCR2(false);
	bool isCR3(false);
	bool isCR4(false);
	bool isCR5(false);
	bool isCR6(false);
	bool isCR7(false);
	bool isCR8(false);
        bool isCR9(false);
        bool isCR10(false);  
        //ABCD method vkg estimation
        bool isQCDVR1(false);
        bool isQCDVR2(false);
        bool isQCDVR3(false);

 
	bool lep_miniso;
	if(isDL) {  lep_miniso = (l1_minisoall<miniso_cut && l2_minisoall<miniso_cut)?true:false; }
	else {  (lep_miniso = l1_minisoall<miniso_cut)?true:false; }
	 
    bool Z_veto = ((l1l2_mass>10. && l1l2_mass<Z_mass_min) || (l1l2_mass>Z_mass_max));
    bool Z_pass = (l1l2_mass>=Z_mass_min && l1l2_mass<=Z_mass_max);
  
	// Apply 30 GeV cut on W candidate soft-drop mass in SL channel (PNet is trained only after mSD > 30 GeV) //
	if(!isDL){	
		if(kl==0) { if( W_msoftdrop_opt2 < msd_cut) continue; }
		if(kl==1) { if( W_msoftdrop_opt1 < msd_cut) continue; }
	}
	// end of W candidate soft-drop mass cut //
  
	if(isDL){
		//signal regions//
		isSR1 = (Flag_Y_bb_pass_T && Z_veto && (l1l2_dR<0.8) && Flag_MET_pass && lep_miniso && abs(dphi_MET_l1l2)<0.5*M_PI);
		isSR2 = (!Flag_Y_bb_pass_T && Z_veto && (l1l2_dR<0.8) && Flag_MET_pass && lep_miniso && abs(dphi_MET_l1l2)<0.5*M_PI);
		//TT CRs //
		isCR2 = (!Flag_Y_bb_pass_T && Z_veto && (l1l2_dR>1.0) && Flag_MET_pass && lep_miniso);// && abs(dphi_MET_l1l2)<0.5*M_PI);
		isCR6 = ( Flag_Y_bb_pass_T && Z_veto && (l1l2_dR>1.0) && Flag_MET_pass && lep_miniso);// && abs(dphi_MET_l1l2)<0.5*M_PI);
		isCR8 = (					  Z_veto && (l1l2_dR>1.0) && Flag_MET_pass && lep_miniso && abs(dphi_MET_l1l2)<0.5*M_PI);
		//DY+j CRs //
		isCR3 = (!Flag_Y_bb_pass_T && Z_pass && (l1l2_dR<1.0) && !Flag_MET_pass && lep_miniso);// && abs(dphi_MET_l1l2)>0.5*M_PI);
		isCR4 = (!Flag_Y_bb_pass_T && Z_pass && (l1l2_dR<1.0) && Flag_MET_pass && lep_miniso);// && abs(dphi_MET_l1l2)<0.5*M_PI); //dphi_MET_l1l2 cut does not have any impact (for CR3 & CR4)
		isCR7 = ( Flag_Y_bb_pass_T && Z_pass && (l1l2_dR<1.0) && lep_miniso);// && abs(dphi_MET_l1l2)<0.5*M_PI);
		//QCD CR //
		isCR5 = (!Flag_Y_bb_pass_T && Z_veto && (l1l2_dR<0.8) && !Flag_MET_pass && !lep_miniso && abs(dphi_MET_l1l2)<0.5*M_PI);
	}
	
	else{
  
		if(kl==0)
		{   
			isSR1 = (Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass && lep_miniso);
			isSR2 = (!Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass && lep_miniso);
			isCR2 = (!Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass && lep_miniso);
			isCR3 = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass && lep_miniso);
			isCR4 = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && lep_miniso && (Y_DeepTag_PNet_TvsQCD>=PN_Top_med));
			isCR5 = (!Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && !Flag_MET_pass && MET_pt > 30  && !lep_miniso &&  l1_minisoall < 0.25);
			isCR6 = (Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass && lep_miniso);
			isCR7 = (Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && !Flag_MET_pass && lep_miniso);
			isCR8  = (Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && Flag_MET_pass && !lep_miniso);
                        isCR9  = (Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && !Flag_MET_pass && lep_miniso);
                        isCR10 = (Flag_H_W_pass_T_opt2 && Flag_dR_lW_pass_opt2 && !Flag_MET_pass && !lep_miniso);
                        isQCDVR1 = (!Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass && lep_miniso);
                        isQCDVR2 = (!Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && Flag_MET_pass && !lep_miniso &&  l1_minisoall < 0.25 );
                        isQCDVR3 = (!Flag_H_W_pass_T_opt2 && !Flag_dR_lW_pass_opt2 && !Flag_MET_pass && MET_pt > 30  && lep_miniso);
		}
		else
		{
			isSR1 = (Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass && lep_miniso);
			isSR2 = (!Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass && lep_miniso);
			isCR2 = (!Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass && lep_miniso);
			isCR3 = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass && lep_miniso);
			isCR4 = (!Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && lep_miniso && (Y_DeepTag_PNet_TvsQCD>=PN_Top_med));
			isCR5 = (!Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && !Flag_MET_pass && !lep_miniso);
			isCR6 = (Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass && lep_miniso); 
			isCR7 = (Flag_Y_bb_pass_T && !Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && !Flag_MET_pass && lep_miniso);
			isCR8 = (Flag_Y_bb_pass_T && Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass && !lep_miniso);  
                        isCR8  = (Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && Flag_MET_pass && !lep_miniso);
                        isCR9  = (Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && !Flag_MET_pass && lep_miniso);
                        isCR10 = (Flag_H_W_pass_T_opt1 && Flag_dR_lW_pass_opt1 && !Flag_MET_pass && !lep_miniso);
                        isQCDVR1 = (!Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass && lep_miniso);
                        isQCDVR2 = (!Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && Flag_MET_pass && !lep_miniso);
                        isQCDVR3 = (!Flag_H_W_pass_T_opt1 && !Flag_dR_lW_pass_opt1 && !Flag_MET_pass && lep_miniso);

		}
	}
   
	reg_tags.push_back(isSR1);
	reg_tags.push_back(isSR2);
	reg_tags.push_back(isCR2);
	reg_tags.push_back(isCR3);
	reg_tags.push_back(isCR4);
	reg_tags.push_back(isCR5);
	reg_tags.push_back(isCR6);
	reg_tags.push_back(isCR7);
	reg_tags.push_back(isCR8);
	reg_tags.push_back(isCR9);
        reg_tags.push_back(isCR10);
        reg_tags.push_back(isQCDVR1);
        reg_tags.push_back(isQCDVR2);
        reg_tags.push_back(isQCDVR3);
	//int ireg = get_region(reg_tags);
	
	vector<int> regs;
	regs = get_regions(reg_tags);
	if(regs.size()<1) continue;
	
	for(unsigned ir=0; ir<regs.size(); ir++){
   
		int ireg = regs[ir];
   
		if(ireg<0||ireg>=nrgn) continue;
	
		// conditions to purify control regions //
	
		bool pure_pass = true;
	
                if( isSR2 || isCR8 || isCR9 || isCR10 ) { pure_pass = (Y_DeepTag_PNetMD_XbbvsQCD>=0.4); }	
		if(!isDL){
		
			if(isCR2){	pure_pass = (Y_DeepTag_PNetMD_XbbvsQCD>=0.1);	}
			if(isCR3){										
				if(kl==0) { pure_pass = (W_msoftdrop_opt2<=60.); }
				if(kl==1) { pure_pass = (W_msoftdrop_opt1<=60.); }
			}
			if(isCR4){	pure_pass = (Y_DeepTag_PNetMD_XbbvsQCD>=0.1);	}
			if(isCR6){	pure_pass = (MET_pt>=75.);	}
	
		}

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
   
		int jk_b = -1;
		if(useBout){
			jk_b = (nbjets_outY==0)?1:2;
		}
		else{
			jk_b = (nbjets_other==0)?1:2;
		}
	   
		float MT = sqrt(2*l1_pt*MET_pt*(1-cos(PhiInRange(l1_phi-MET_phi))));
	   
		bool lepid_info[nlid];
		lepid_info[0] = true;
		if(isDL){
			if(abs(l1_pdgId)==13 && abs(l2_pdgId)==13) { lepid_info[1] = true;  }
			if(abs(l1_pdgId)==11 && abs(l2_pdgId)==11) { lepid_info[2] = true;  }
			if((abs(l1_pdgId)==11 && abs(l2_pdgId)==13) || (abs(l1_pdgId)==13 && abs(l2_pdgId)==11)) { lepid_info[3] = true;  }
		}
		else{
			if(abs(l1_pdgId)==13) { lepid_info[1] = true;  }
			if(abs(l1_pdgId)==11) { lepid_info[2] = true;  }
		}
	        	
		// first fill few general histograms //
	   
		h_reg->Fill(ireg,weight_nom);
	   
		// now fill histograms binned in analysis categories //
	   
		float weight = weight_nom;
	        //std::cout << weight << std::endl;	
		if(isCR4 && !isDL) { weight = weight_nom*Top_SF; } // since top tagging condition is used only in CR4
		for(int jk=0; jk<nbcat; jk++){
		   
			if(!(jk==0 || jk==jk_b)) continue;
			
			if(jk!=0) { weight = weight_nom*b_SF; } // applying b tagging SF only if any condition on number of b-tagged jets is used
			
			for(int lm=0; lm<nlid; lm++){
				
				if(!lepid_info[lm]) continue;
				
				// reference histograms corresponding to purity conditions //
				
				if(!isSignal){
					h_Y_PNetMD_XbbvsQCD_ref[ireg][jk][kl][lm]->Fill(Y_DeepTag_PNetMD_XbbvsQCD,weight);
					h_MET_pt_ref[ireg][jk][kl][lm]->Fill(MET_pt,weight); 
					if(kl==0){  
						h_W_msoftdrop_ref[ireg][jk][kl][lm]->Fill(W_msoftdrop_opt2,weight); 
					}
					if(kl==1){  
						h_W_msoftdrop_ref[ireg][jk][kl][lm]->Fill(W_msoftdrop_opt1,weight); 
					}
				}
				// ref end //
				
				// Now apply the purity condition & fill all histograms (including sys)
				
				if(pure_pass)
				{
	   
					if(!isSignal){
					
					h_l1_pt[ireg][jk][kl][lm]->Fill(l1_pt,weight); 
					h_l1_eta[ireg][jk][kl][lm]->Fill(l1_eta,weight); 
					h_l1_minisoall[ireg][jk][kl][lm]->Fill(l1_minisoall,weight); 
					
					if(isDL){
					
						h_l2_pt[ireg][jk][kl][lm]->Fill(l2_pt,weight); 
						h_l2_eta[ireg][jk][kl][lm]->Fill(l2_eta,weight); 
						h_l2_minisoall[ireg][jk][kl][lm]->Fill(l2_minisoall,weight); 
					
						h_l1l2_mass[ireg][jk][kl][lm]->Fill(l1l2_mass,weight); 
						h_l1l2_dR[ireg][jk][kl][lm]->Fill(l1l2_dR,weight); 
						h_l1l2_deta[ireg][jk][kl][lm]->Fill(l1l2_deta,weight); 
						h_l1l2_dphi[ireg][jk][kl][lm]->Fill(l1l2_dphi,weight); 
						
						h_dphi_MET_l1l2[ireg][jk][kl][lm]->Fill(dphi_MET_l1l2,weight); 
					
					}
						
					h_MET_pt[ireg][jk][kl][lm]->Fill(MET_pt,weight); 
					h_MET_sig[ireg][jk][kl][lm]->Fill(MET_sig,weight); 
					
					}
					
					h_Y_pt[ireg][jk][kl][lm]->Fill(Y_pt,weight); 
					//h_Y_msoftdrop[ireg][jk][kl][lm]->Fill(Y_msoftdrop,weight);
					//h_Y_msoftdrop_xbin[ireg][jk][kl][lm]->Fill(Y_msoftdrop,weight);
					//h_Y_PNetMD_XbbvsQCD[ireg][jk][kl][lm]->Fill(Y_DeepTag_PNetMD_XbbvsQCD,weight);
					
					if(!isSignal){
					
					h_Y_PNetMD_WvsQCD[ireg][jk][kl][lm]->Fill(Y_DeepTag_PNetMD_WvsQCD,weight);
					h_Y_PNet_TvsQCD[ireg][jk][kl][lm]->Fill(Y_DeepTag_PNet_TvsQCD,weight);
					h_Y_sub1_mass[ireg][jk][kl][lm]->Fill(Y_sub1_mass,weight);
					h_Y_sub2_mass[ireg][jk][kl][lm]->Fill(Y_sub2_mass,weight);
					h_Y_sub1_btag[ireg][jk][kl][lm]->Fill(Y_sub1_btag,weight);
					h_Y_sub2_btag[ireg][jk][kl][lm]->Fill(Y_sub2_btag,weight);
					
					}
					
					//h_HTlep_pt[ireg][jk][kl][lm]->Fill(HTlep_pt,weight); 
					//h_ST[ireg][jk][kl][lm]->Fill(ST,weight); 
							
					float X_conv_mass;
					
					float X_mass;
					if(kl==0) { X_mass = X_mass_opt2; }
					else { X_mass = X_mass_opt1; }
					float unrol_mass = -1.0;
					if(X_mass>=invmassbins[0] && X_mass<invmassbins[ninvmassbins] && Y_msoftdrop>=msdbins[0] && Y_msoftdrop<msdbins[nmsdbins]){
						//unrol_mass = float(1.0*getbinid(Y_msoftdrop,nmsdbins,msdbins)*getbinid(X_mass,ninvmassbins,invmassbins));
						unrol_mass = float(getbinid(X_mass,ninvmassbins,invmassbins) + getbinid(Y_msoftdrop,nmsdbins,msdbins)* ninvmassbins);
					}
					float unrol_ST = -1.0;
					if(ST>=invmassbins[0] && ST<invmassbins[ninvmassbins] && Y_msoftdrop>=msdbins[0] && Y_msoftdrop<msdbins[nmsdbins]){
						//unrol_ST = float(1.0*getbinid(Y_msoftdrop,nmsdbins,msdbins)*getbinid(ST,ninvmassbins,invmassbins));
						unrol_ST = float(getbinid(ST,ninvmassbins,invmassbins) + getbinid(Y_msoftdrop,nmsdbins,msdbins)* ninvmassbins);
					}
                    bool top_semimerged = false; 
                    if(!isDATA) { top_semimerged  = Y_label_Top_bq || Y_label_Top_bc || Y_label_W_qq || Y_label_W_cq  ;}
                    bool top_fullymerged = false; 
                    if(!isDATA ) { top_fullymerged = Y_label_Top_bcq || Y_label_Top_bqq;}

					if(kl==0){           
				
						h_W_msoftdrop[ireg][jk][kl][lm]->Fill(W_msoftdrop_opt2,weight); 
				
						if(!isSignal){
				
						h_W_pt[ireg][jk][kl][lm]->Fill(W_pt_opt2,weight); 
						//h_W_msoftdrop[ireg][jk][kl][lm]->Fill(W_msoftdrop_opt2,weight); 
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
					
						}
					}
					else
					{
				
						h_W_msoftdrop[ireg][jk][kl][lm]->Fill(W_msoftdrop_opt1,weight); 
				
						if(!isSignal){
				
						h_W_pt[ireg][jk][kl][lm]->Fill(W_pt_opt1,weight); 
						//h_W_msoftdrop[ireg][jk][kl][lm]->Fill(W_msoftdrop_opt1,weight); 
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
				
						}
				
					}
					
					if(!isSignal){
					
					h_X_Y_mass[ireg][jk][kl][lm]->Fill(X_mass, Y_msoftdrop, weight);
					h_X_Y_mass_xbin[ireg][jk][kl][lm]->Fill(X_mass, Y_msoftdrop, weight);
					   
					}
					   
					X_conv_mass = X_mass + 4000.0 * get_Y_id(Y_msoftdrop);
			
					h_nbjets_other[ireg][jk][kl][lm]->Fill(nbjets_other,weight);
					h_nbjets_outY[ireg][jk][kl][lm]->Fill(nbjets_outY,weight);
					h_nbjets_outY_L[ireg][jk][kl][lm]->Fill(nbjets_outY_L,weight);
					h_nbjets[ireg][jk][kl][lm]->Fill(nbjets,weight);
					h_nbjets_L[ireg][jk][kl][lm]->Fill(nbjets_L,weight);
				
	            //cout << "check! dude: 1" << endl;			
                    for(int mn = 0 ; mn < ntop ; mn++)
                    {
						if( mn == 0 || ( mn == 1 && top_fullymerged) || ( mn == 2 && !top_fullymerged && top_semimerged) || ( mn == 3 && !top_fullymerged && !top_semimerged)  )
                        {
							h_Y_msoftdrop[ireg][jk][kl][lm][mn]->Fill(Y_msoftdrop,weight);
							if(!isSignal){
								h_Y_mass[ireg][jk][kl][lm][mn]->Fill(Y_mass,weight);
								h_Y_PNetMD_XbbvsQCD[ireg][jk][kl][lm][mn]->Fill(Y_DeepTag_PNetMD_XbbvsQCD,weight);
							}
                            h_HTlep_pt[ireg][jk][kl][lm][mn]->Fill(HTlep_pt,weight);
                            h_ST[ireg][jk][kl][lm][mn]->Fill(ST,weight);
                            h_MT[ireg][jk][kl][lm][mn]->Fill(MT,weight);					
                            h_Y_msoftdrop_xbin[ireg][jk][kl][lm][mn]->Fill(Y_msoftdrop,weight);
                            h_Y_mass_xbin[ireg][jk][kl][lm][mn]->Fill(Y_msoftdrop,weight);
                            h_ST_xbin[ireg][jk][kl][lm][mn]->Fill(ST,weight);		
							if(!isSignal){
								h_Y_msoftdrop_sys[ireg][jk][kl][lm][mn][0]->Fill(Y_msoftdrop,weight);
								h_Y_msoftdrop_xbin_sys[ireg][jk][kl][lm][mn][0]->Fill(Y_msoftdrop,weight);
								h_HTlep_pt_sys[ireg][jk][kl][lm][mn][0]->Fill(HTlep_pt,weight);
								h_ST_sys[ireg][jk][kl][lm][mn][0]->Fill(ST,weight);
                                                                h_ST_xbin_sys[ireg][jk][kl][lm][mn][0]->Fill(ST,weight);
                                                                h_MT_sys[ireg][jk][kl][lm][mn][0]->Fill(MT,weight);
							}
							
                            h_for_limit_HTlep_pt_sys[ireg][jk][kl][lm][mn][0]->Fill(HTlep_pt + 4000.0 * get_Y_id(Y_msoftdrop),weight);
                            h_for_limit_ST_sys[ireg][jk][kl][lm][mn][0]->Fill(ST + 4000.0 * get_Y_id(Y_msoftdrop),weight);
                            h_for_limit_ST_sys_v2[ireg][jk][kl][lm][mn][0]->Fill(unrol_ST,weight);

                            h_for_limit_HTlep_pt[ireg][jk][kl][lm][mn]->Fill(HTlep_pt + 4000.0 * get_Y_id(Y_msoftdrop),weight);
                            h_for_limit_ST[ireg][jk][kl][lm][mn]->Fill(ST + 4000.0 * get_Y_id(Y_msoftdrop),weight);
                            h_for_limit_ST_v2[ireg][jk][kl][lm][mn]->Fill(unrol_ST,weight);
                            
                            h_ST_Y_mass_sys[ireg][jk][kl][lm][mn][0]->Fill(ST, Y_msoftdrop, weight);
                            
                            if(!isDL){
                            
								h_X_mass[ireg][jk][kl][lm][mn]->Fill(X_mass,weight);
								h_for_limit_X_mass[ireg][jk][kl][lm][mn]->Fill(X_conv_mass,weight);
								h_for_limit_X_mass_v2[ireg][jk][kl][lm][mn]->Fill(unrol_mass,weight);
							        h_X_mass_xbin[ireg][jk][kl][lm][mn]->Fill(X_mass,weight);	
								if(!isSignal){
									h_X_mass_sys[ireg][jk][kl][lm][mn][0]->Fill(X_mass,weight);
									h_X_mass_xbin_sys[ireg][jk][kl][lm][mn][0]->Fill(X_mass,weight);
								}
								
								h_for_limit_X_mass_sys[ireg][jk][kl][lm][mn][0]->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight);
								h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][mn][0]->Fill(unrol_mass,weight);
                            
								h_X_Y_mass_sys[ireg][jk][kl][lm][mn][0]->Fill(X_mass, Y_msoftdrop, weight);
								
				    }
			   }
                       }//mn
                       //std::cout << "check! dude: 2" << std::endl; 
                       if(!isDATA) 
                       {
					vector<float> X_JESup_split, X_JESdn_split;
					
					/*
					TLorentzVector Y_cand, H_cand;
		            if (!isDL){		        	
						for(int ijes=0; ijes<njecmax; ijes++){
							if(kl==0){
								H_cand.SetPtEtaPhiM(H_pt_opt2*(*H_JESup_split_opt2)[ijes],H_eta_opt2,H_phi_opt2,H_mass_opt2*(*H_JESup_split_opt2)[ijes]);
								Y_cand.SetPtEtaPhiM(Y_pt*(*Y_JESup_split)[ijes],Y_eta,Y_phi,Y_mass*(*Y_JESup_split)[ijes]);
								X_JESup_split.push_back((Y_cand+H_cand).M()/X_mass);
								H_cand.SetPtEtaPhiM(H_pt_opt2*(*H_JESdn_split_opt2)[ijes],H_eta_opt2,H_phi_opt2,H_mass_opt2*(*H_JESdn_split_opt2)[ijes]);
								Y_cand.SetPtEtaPhiM(Y_pt*(*Y_JESdn_split)[ijes],Y_eta,Y_phi,Y_mass*(*Y_JESdn_split)[ijes]);
								X_JESdn_split.push_back((Y_cand+H_cand).M()/X_mass);
							}
							else{
								H_cand.SetPtEtaPhiM(H_pt_opt1*(*H_JESup_split_opt1)[ijes],H_eta_opt1,H_phi_opt1,H_mass_opt1*(*H_JESup_split_opt1)[ijes]);
								Y_cand.SetPtEtaPhiM(Y_pt*(*Y_JESup_split)[ijes],Y_eta,Y_phi,Y_mass*(*Y_JESup_split)[ijes]);
								X_JESup_split.push_back((Y_cand+H_cand).M()/X_mass);
								H_cand.SetPtEtaPhiM(H_pt_opt1*(*H_JESdn_split_opt1)[ijes],H_eta_opt1,H_phi_opt1,H_mass_opt1*(*H_JESdn_split_opt1)[ijes]);
								Y_cand.SetPtEtaPhiM(Y_pt*(*Y_JESdn_split)[ijes],Y_eta,Y_phi,Y_mass*(*Y_JESdn_split)[ijes]);
								X_JESdn_split.push_back((Y_cand+H_cand).M()/X_mass);
							}
						}
                    }
				    */
				    
				    if (!isDL && !isDATA){	
						for(int ijes=0; ijes<njecmax; ijes++){
						//	if(kl==0){
								X_JESup_split.push_back((*X_mass_JESup_split_opt2)[ijes]);
								X_JESdn_split.push_back((*X_mass_JESdn_split_opt2)[ijes]);
						//	}
						//	else{
						//		X_JESup_split.push_back((*X_mass_JESup_split_opt1)[ijes]);
						//		X_JESdn_split.push_back((*X_mass_JESdn_split_opt1)[ijes]);
						//	}
						}
					}
				      
				        //cout << "check! dude: 3" << endl;	
					vector<float> shape_weight_up, shape_weight_dn, shape_weight_nom;
					
					shape_weight_up.push_back(puWeightup);				shape_weight_dn.push_back(puWeightdown);  			shape_weight_nom.push_back(puWeight);
					shape_weight_up.push_back(leptonsf_weight_stat);	shape_weight_dn.push_back(leptonsf_weight_syst);	shape_weight_nom.push_back(leptonsf_weight);
					shape_weight_up.push_back(leptonsf_weight_up);		shape_weight_dn.push_back(leptonsf_weight_dn);		shape_weight_nom.push_back(leptonsf_weight);
					shape_weight_up.push_back(prefiringweightup);		shape_weight_dn.push_back(prefiringweightdown);		shape_weight_nom.push_back(prefiringweight);
					shape_weight_up.push_back(bb_SF_up);				shape_weight_dn.push_back(bb_SF_dn);				shape_weight_nom.push_back(bb_SF);
					if(!isDL){
					shape_weight_up.push_back(W_SF_up);					shape_weight_dn.push_back(W_SF_dn); 				shape_weight_nom.push_back(W_SF); // for single-lepton channel only
					}
					else{
						shape_weight_up.push_back(1);					shape_weight_dn.push_back(1); 				shape_weight_nom.push_back(1); 
					}
					if(jk!=0){
						shape_weight_up.push_back(b_SF_up);					shape_weight_dn.push_back(b_SF_dn); 				shape_weight_nom.push_back(b_SF); // only with condition on # of b jets
					}
					else{
						shape_weight_up.push_back(1);					shape_weight_dn.push_back(1); 				shape_weight_nom.push_back(1);
					}
					shape_weight_up.push_back(SF_Trig_1_up);			shape_weight_dn.push_back(SF_Trig_1_dn);			shape_weight_nom.push_back(SF_Trig);
					shape_weight_up.push_back(SF_Trig_2_up);			shape_weight_dn.push_back(SF_Trig_2_dn);			shape_weight_nom.push_back(SF_Trig);
					
		                        //cout << "check! dude: 4 " << endl;			
   				    for(int isys=0; isys<nsys; isys++){
						for(int mn = 0 ; mn < ntop ; mn++)
						{
							bool top_cond = false;
							if( mn == 0 || ( mn == 1 && top_fullymerged) || ( mn == 2 && !top_fullymerged && top_semimerged) || ( mn == 3 && !top_fullymerged && !top_semimerged)  ) { top_cond = true; }
							if(!top_cond) continue; 
				            //cout << "check! dude" << endl;			
     					    if(isys==0 || isys<njecmax) // JES + JER (splitted in sources)
		 				    {
						
								if(!isSignal){
								
									h_Y_msoftdrop_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(Y_msoftdrop*(*Y_JESup_split)[isys],weight);
									h_Y_msoftdrop_sys[ireg][jk][kl][lm][mn][2*(isys+1)]	 ->Fill(Y_msoftdrop*(*Y_JESdn_split)[isys],weight);
							
									h_Y_msoftdrop_xbin_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(Y_msoftdrop*(*Y_JESup_split)[isys],weight);
									h_Y_msoftdrop_xbin_sys[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(Y_msoftdrop*(*Y_JESdn_split)[isys],weight);
							
									h_HTlep_pt_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(HTlep_pt*(*HTlep_pt_JESup_split)[isys],weight); 
									h_HTlep_pt_sys[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(HTlep_pt*(*HTlep_pt_JESdn_split)[isys],weight); 
							
									h_ST_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(ST*(*ST_JESup_split)[isys],weight); 
									h_ST_sys[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(ST*(*ST_JESdn_split)[isys],weight);
                                                                        
									//h_ST_xbin_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(ST*(*ST_JESup_split)[isys],weight);
                                                                        //h_ST_xbin_sys[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(ST*(*ST_JESdn_split)[isys],weight);

				                                        float MT_JESup_split = sqrt(2*l1_pt*(*MET_pt_JESup_split)[isys]*(1-cos(PhiInRange(l1_phi- (*MET_phi_JESup_split)[isys]))));
                                                                        float MT_JESdn_split = sqrt(2*l1_pt*(*MET_pt_JESdn_split)[isys]*(1-cos(PhiInRange(l1_phi- (*MET_phi_JESdn_split)[isys]))));
                                                                        h_MT_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(MT_JESup_split,weight);
                                                                        h_MT_sys[ireg][jk][kl][lm][mn][2*(isys+1)]->Fill(MT_JESdn_split,weight);	
								}
				                                //cout << "check! dude" << endl;	
								h_for_limit_HTlep_pt_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(HTlep_pt*(*HTlep_pt_JESup_split)[isys] + 4000.0 * get_Y_id(Y_msoftdrop*(*Y_JESup_split)[isys]),weight);
								h_for_limit_HTlep_pt_sys[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(HTlep_pt*(*HTlep_pt_JESdn_split)[isys] + 4000.0 * get_Y_id(Y_msoftdrop*(*Y_JESdn_split)[isys]),weight);
							
								h_for_limit_ST_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(ST*(*ST_JESup_split)[isys] + 4000.0 * get_Y_id(Y_msoftdrop*(*Y_JESup_split)[isys]),weight);
								h_for_limit_ST_sys[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(ST*(*ST_JESdn_split)[isys] + 4000.0 * get_Y_id(Y_msoftdrop*(*Y_JESdn_split)[isys]),weight);
							
								if((ST*(*ST_JESup_split)[isys])>=invmassbins[0] && (ST*(*ST_JESup_split)[isys])<invmassbins[ninvmassbins] && Y_msoftdrop*(*Y_JESup_split)[isys]>=msdbins[0] && Y_msoftdrop*(*Y_JESup_split)[isys]<msdbins[nmsdbins]){
								         h_for_limit_ST_sys_v2[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(float(getbinid((ST*(*ST_JESup_split)[isys]),ninvmassbins,invmassbins) + getbinid(Y_msoftdrop*(*Y_JESup_split)[isys],nmsdbins,msdbins)* ninvmassbins),weight);
								}
								if((ST*(*ST_JESdn_split)[isys])>=invmassbins[0] && (ST*(*ST_JESdn_split)[isys])<invmassbins[ninvmassbins] && Y_msoftdrop*(*Y_JESdn_split)[isys]>=msdbins[0] && Y_msoftdrop*(*Y_JESdn_split)[isys]<msdbins[nmsdbins]){
									h_for_limit_ST_sys_v2[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(float(getbinid((ST*(*ST_JESdn_split)[isys]),ninvmassbins,invmassbins) + getbinid(Y_msoftdrop*(*Y_JESdn_split)[isys],nmsdbins,msdbins)* ninvmassbins),weight);
								}
							
								h_ST_Y_mass_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(ST*(*ST_JESup_split)[isys], Y_msoftdrop*(*Y_JESup_split)[isys], weight);
								h_ST_Y_mass_sys[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(ST*(*ST_JESdn_split)[isys], Y_msoftdrop*(*Y_JESdn_split)[isys], weight);
							
								if(!isDL){
							
									if(!isSignal){
							
										h_X_mass_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(X_mass*X_JESup_split[isys],weight); 
										h_X_mass_sys[ireg][jk][kl][lm][mn][2*(isys+1)]	->Fill(X_mass*X_JESdn_split[isys],weight); 
							
										h_X_mass_xbin_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(X_mass*X_JESup_split[isys],weight); 
										h_X_mass_xbin_sys[ireg][jk][kl][lm][mn][2*(isys+1)]	->Fill(X_mass*X_JESdn_split[isys],weight); 
							
									}
							
									h_for_limit_X_mass_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(X_mass*X_JESup_split[isys] + 4000.0 * get_Y_id(Y_msoftdrop*(*Y_JESup_split)[isys]),weight);
									h_for_limit_X_mass_sys[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(X_mass*X_JESdn_split[isys] + 4000.0 * get_Y_id(Y_msoftdrop*(*Y_JESdn_split)[isys]),weight);
									if(X_mass*X_JESup_split[isys]>=invmassbins[0] && X_mass*X_JESup_split[isys]<invmassbins[ninvmassbins] && Y_msoftdrop*(*Y_JESup_split)[isys]>=msdbins[0] && Y_msoftdrop*(*Y_JESup_split)[isys]<msdbins[nmsdbins])
								    {h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(float(getbinid(X_mass*X_JESup_split[isys],ninvmassbins,invmassbins) + getbinid(Y_msoftdrop*(*Y_JESup_split)[isys],nmsdbins,msdbins)* ninvmassbins),weight);}
									if(X_mass*X_JESdn_split[isys]>=invmassbins[0] && X_mass*X_JESdn_split[isys]<invmassbins[ninvmassbins] && Y_msoftdrop*(*Y_JESdn_split)[isys]>=msdbins[0] && Y_msoftdrop*(*Y_JESdn_split)[isys]<msdbins[nmsdbins])
								    {h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(float(getbinid(X_mass*X_JESdn_split[isys],ninvmassbins,invmassbins) + getbinid(Y_msoftdrop*(*Y_JESdn_split)[isys],nmsdbins,msdbins)* ninvmassbins),weight);}
								
									h_X_Y_mass_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(X_mass*X_JESup_split[isys], Y_msoftdrop*(*Y_JESup_split)[isys], weight);
									h_X_Y_mass_sys[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(X_mass*X_JESdn_split[isys], Y_msoftdrop*(*Y_JESdn_split)[isys], weight);
					
								}
							
							}
                                                        
                            X_JESup_split.clear();
                            X_JESdn_split.clear();  

							//cout << "check! dude-5" << endl;
							if(isys>=njecmax) // Weight-based systematics (SFs)
							{
							
								int isys1 = isys-(njecmax);
							
								float weight_up = weight*shape_weight_up[isys1]/TMath::Max(float(1.e-6),shape_weight_nom[isys1]);
								float weight_dn = weight*shape_weight_dn[isys1]/TMath::Max(float(1.e-6),shape_weight_nom[isys1]);
							
								if(jk==0 && isys1==6) {  weight_up = 1.; weight_dn = 1.; } // b tagging SF
							
								if(!isSignal){
							
									h_Y_msoftdrop_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(Y_msoftdrop,weight_up);
									h_Y_msoftdrop_sys[ireg][jk][kl][lm][mn][2*(isys+1)]	 ->Fill(Y_msoftdrop,weight_dn);
							
									h_Y_msoftdrop_xbin_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(Y_msoftdrop,weight_up);
									h_Y_msoftdrop_xbin_sys[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(Y_msoftdrop,weight_dn);

									h_HTlep_pt_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(HTlep_pt,weight_up); 
									h_HTlep_pt_sys[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(HTlep_pt,weight_dn); 
							
									h_ST_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(ST,weight_up); 
									h_ST_sys[ireg][jk][kl][lm][mn][2*(isys+1)]	->Fill(ST,weight_dn);

                                                                        h_ST_xbin_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(ST,weight_up);
                                                                        h_ST_xbin_sys[ireg][jk][kl][lm][mn][2*(isys+1)]      ->Fill(ST,weight_dn);
								 	
							                h_MT_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(MT,weight_up);
                                                                        h_MT_sys[ireg][jk][kl][lm][mn][2*(isys+1)]      ->Fill(MT,weight_dn);	
								}
								
								h_for_limit_HTlep_pt_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(HTlep_pt + 4000.0 * get_Y_id(Y_msoftdrop),weight_up);
								h_for_limit_HTlep_pt_sys[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(HTlep_pt + 4000.0 * get_Y_id(Y_msoftdrop),weight_dn);
							
								h_for_limit_ST_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(ST + 4000.0 * get_Y_id(Y_msoftdrop),weight_up);
								h_for_limit_ST_sys[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(ST + 4000.0 * get_Y_id(Y_msoftdrop),weight_dn);
								
								h_for_limit_ST_sys_v2[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(unrol_ST,weight_up);
								h_for_limit_ST_sys_v2[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(unrol_ST,weight_dn);
							
								h_ST_Y_mass_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(ST, Y_msoftdrop, weight_up);
								h_ST_Y_mass_sys[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(ST, Y_msoftdrop, weight_dn);
							
								if(!isDL){
							
									if(!isSignal){
							
										h_X_mass_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(X_mass,weight_up); 
										h_X_mass_sys[ireg][jk][kl][lm][mn][2*(isys+1)]	->Fill(X_mass,weight_dn); 
					
										h_X_mass_xbin_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(X_mass,weight_up); 
										h_X_mass_xbin_sys[ireg][jk][kl][lm][mn][2*(isys+1)]	 ->Fill(X_mass,weight_dn); 
				
									}
				
									h_for_limit_X_mass_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight_up);
									h_for_limit_X_mass_sys[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(X_mass + 4000.0 * get_Y_id(Y_msoftdrop),weight_dn);
									h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(unrol_mass,weight_up);
									h_for_limit_X_mass_sys_v2[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(unrol_mass,weight_dn);
														
									h_X_Y_mass_sys[ireg][jk][kl][lm][mn][2*(isys+1)-1]->Fill(X_mass, Y_msoftdrop, weight_up);
									h_X_Y_mass_sys[ireg][jk][kl][lm][mn][2*(isys+1)]  ->Fill(X_mass, Y_msoftdrop, weight_dn);
								
								}
							}
						        //std::cout << "Dude-6" << std::endl;	
							shape_weight_up.clear();
							shape_weight_dn.clear();
							shape_weight_nom.clear();
                            
                                            }//top conditions 

					}//sys loop
		                    } //No Data condition	
				}//purity condition
				
			}// lepid (lm)
			
		  }// b cat (jk)
	  
		}//ir (loop over regions)
	  
	}// W opt (kl)
	
	
   }// end of event loop
   
    final_file->Write();
    final_file->cd();
   
    final_file->Close();
    if(!isDATA)
    {
    	FEff->Close();
    }
  //}//ii

}

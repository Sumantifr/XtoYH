#include "Histo_PNetSF.h"

using namespace std;

//int main()
//void Histo_PNetSF()
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
   float msdbins[nmsdbins];
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
   float invmassbins[ninvmassbins];
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

   TString inputFile=argv[3];
   std::cout << inputFile << std::endl;
   TFile* final_file = TFile::Open("PFNetSF_OUTPUTS/Histogram_"+inputFile, "RECREATE");  

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
  
  float ak8ptbins[] = {200,300,400,600,10000};
  int nak8ptbins =  (sizeof(ak8ptbins)/sizeof(ak8ptbins[0]))-1;
  
  char name[1000];

  TH1F* h_AK8_msoftdrop_sys_Top_fullymerged_pass[nbcat][nptbins][1+2*nsys];
  TH1F* h_AK8_msoftdrop_sys_Top_fullymerged_v2_pass[nbcat][nptbins][1+2*nsys];
  TH1F* h_AK8_msoftdrop_sys_Top_semimerged_pass[nbcat][nptbins][1+2*nsys];
  TH1F* h_AK8_msoftdrop_sys_W_hadronic_pass[nbcat][nptbins][1+2*nsys];
  TH1F* h_AK8_msoftdrop_sys_rest_pass[nbcat][nptbins][1+2*nsys];
  TH1F* h_AK8_msoftdrop_sys_rest_v2_pass[nbcat][nptbins][1+2*nsys]; 

  TH1F* h_AK8_msoftdrop_sys_Top_fullymerged_fail[nbcat][nptbins][1+2*nsys];
  TH1F* h_AK8_msoftdrop_sys_Top_fullymerged_v2_fail[nbcat][nptbins][1+2*nsys];
  TH1F* h_AK8_msoftdrop_sys_Top_semimerged_fail[nbcat][nptbins][1+2*nsys];
  TH1F* h_AK8_msoftdrop_sys_W_hadronic_fail[nbcat][nptbins][1+2*nsys];
  TH1F* h_AK8_msoftdrop_sys_rest_fail[nbcat][nptbins][1+2*nsys];
  TH1F* h_AK8_msoftdrop_sys_rest_v2_fail[nbcat][nptbins][1+2*nsys];
  
  TH1F* h_AK8_msoftdrop_pass[nbcat][nptbins];
  TH1F* h_AK8_msoftdrop_fail[nbcat][nptbins];

  for(int jk=0; jk<nbcat; jk++)
  {
	  
	for(int kl=0; kl<nak8ptbins; kl++){ 
	
		sprintf(name,"_ptbin%i",kl+1);
		  
		h_AK8_msoftdrop_pass[jk][kl] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_pass",name,38,30,600);
		h_AK8_msoftdrop_fail[jk][kl] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_fail",name,38,30,600);

		sprintf(name,"_ptbin%i_nom",kl+1);

		//nom systematics  
		h_AK8_msoftdrop_sys_Top_fullymerged_pass[jk][kl][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_fullymerged_pass",name,38,30,600);
		h_AK8_msoftdrop_sys_Top_fullymerged_v2_pass[jk][kl][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_fullymerged_v2_pass",name,38,30,600);
		h_AK8_msoftdrop_sys_Top_semimerged_pass[jk][kl][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_semimerged_pass",name,38,30,600);
		h_AK8_msoftdrop_sys_W_hadronic_pass[jk][kl][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_W_hadronic_pass",name,38,30,600);
		h_AK8_msoftdrop_sys_rest_pass[jk][kl][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_rest_pass",name,38,30,600);
		h_AK8_msoftdrop_sys_rest_v2_pass[jk][kl][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_rest_v2_pass",name,38,30,600);

		h_AK8_msoftdrop_sys_Top_fullymerged_fail[jk][kl][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_fullymerged_fail",name,38,30,600);
		h_AK8_msoftdrop_sys_Top_fullymerged_v2_fail[jk][kl][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_fullymerged_v2_fail",name,38,30,600);
		h_AK8_msoftdrop_sys_Top_semimerged_fail[jk][kl][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_semimerged_fail",name,38,30,600);
		h_AK8_msoftdrop_sys_W_hadronic_fail[jk][kl][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_W_hadronic_fail",name,38,30,600);
		h_AK8_msoftdrop_sys_rest_fail[jk][kl][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_rest_fail",name,38,30,600);
		h_AK8_msoftdrop_sys_rest_v2_fail[jk][kl][0] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_rest_v2_fail",name,38,30,600);

		for(int isys=0; isys<nsys; isys++)
		{               
			//up systematics
			sprintf(name,"_ptbin%i_%s_up",kl+1,sysnames[isys].Data());
			h_AK8_msoftdrop_sys_Top_fullymerged_pass[jk][kl][2*(isys+1)-1]  = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_fullymerged_pass",name,38,30,600);
			h_AK8_msoftdrop_sys_Top_fullymerged_v2_pass[jk][kl][2*(isys+1)-1]  = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_fullymerged_v2_pass",name,38,30,600);
			h_AK8_msoftdrop_sys_Top_semimerged_pass[jk][kl][2*(isys+1)-1]   = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_semimerged_pass",name,38,30,600);
			h_AK8_msoftdrop_sys_W_hadronic_pass[jk][kl][2*(isys+1)-1]  = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_W_hadronic_pass",name,38,30,600);
			h_AK8_msoftdrop_sys_rest_pass[jk][kl][2*(isys+1)-1]        = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_rest_pass",name,38,30,600);
			h_AK8_msoftdrop_sys_rest_v2_pass[jk][kl][2*(isys+1)-1]  = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_rest_v2_pass",name,38,30,600);
						   
			h_AK8_msoftdrop_sys_Top_fullymerged_fail[jk][kl][2*(isys+1)-1] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_fullymerged_fail",name,38,30,600);
			h_AK8_msoftdrop_sys_Top_fullymerged_v2_fail[jk][kl][2*(isys+1)-1] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_fullymerged_v2_fail",name,38,30,600);
			h_AK8_msoftdrop_sys_Top_semimerged_fail[jk][kl][2*(isys+1)-1]  = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_semimerged_fail",name,38,30,600);
			h_AK8_msoftdrop_sys_W_hadronic_fail[jk][kl][2*(isys+1)-1]  = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_W_hadronic_fail",name,38,30,600);
			h_AK8_msoftdrop_sys_rest_fail[jk][kl][2*(isys+1)-1]        = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_rest_fail",name,38,30,600);
			h_AK8_msoftdrop_sys_rest_v2_fail[jk][kl][2*(isys+1)-1]  = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_rest_v2_fail",name,38,30,600); 

			//dn systmatics
			sprintf(name,"_ptbin%i_%s_dn",kl+1,sysnames[isys].Data());
			h_AK8_msoftdrop_sys_Top_fullymerged_pass[jk][kl][2*(isys+1)]  = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_fullymerged_pass",name,38,30,600);
			h_AK8_msoftdrop_sys_Top_fullymerged_v2_pass[jk][kl][2*(isys+1)]  = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_fullymerged_v2_pass",name,38,30,600);
			h_AK8_msoftdrop_sys_Top_semimerged_pass[jk][kl][2*(isys+1)]   = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_semimerged_pass",name,38,30,600);
			h_AK8_msoftdrop_sys_W_hadronic_pass[jk][kl][2*(isys+1)]  = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_W_hadronic_pass",name,38,30,600);
			h_AK8_msoftdrop_sys_rest_pass[jk][kl][2*(isys+1)]        = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_rest_pass",name,38,30,600);
			h_AK8_msoftdrop_sys_rest_v2_pass[jk][kl][2*(isys+1)]  = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_rest_v2_pass",name,38,30,600);

			h_AK8_msoftdrop_sys_Top_fullymerged_fail[jk][kl][2*(isys+1)] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_fullymerged_fail",name,38,30,600);
			h_AK8_msoftdrop_sys_Top_fullymerged_v2_fail[jk][kl][2*(isys+1)] = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_fullymerged_v2_fail",name,38,30,600);
			h_AK8_msoftdrop_sys_Top_semimerged_fail[jk][kl][2*(isys+1)]  = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_Top_semimerged_fail",name,38,30,600);
			h_AK8_msoftdrop_sys_W_hadronic_fail[jk][kl][2*(isys+1)]  = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_W_hadronic_fail",name,38,30,600);
			h_AK8_msoftdrop_sys_rest_fail[jk][kl][2*(isys+1)]        = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_rest_fail",name,38,30,600);
			h_AK8_msoftdrop_sys_rest_v2_fail[jk][kl][2*(isys+1)]  = get_histo_symbin_II(Ytype[y_wp],Wtype[w_wp],rgn[ij],bcats[jk],Wops[kl],"AK8_msoftdrop_rest_v2_fail",name,38,30,600);	      
		  
			} 
		}
	}
	 

  // End of declaration //
    
  // Definition of histograms //
  
  h_nom = getHisto1D("h_nom","h_nom",10,0.0,10.0);
 
 
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

    	 
   file->cd();
 
   Long64_t nn = tree->GetEntries();
   
   //// Event loop ////
   std::cout << nn << std::endl;  
   for(Long64_t jentry =0; jentry < nn ; jentry++)
   {
	      
	tree->GetEntry(jentry);
	if( jentry % 10000 == 0) { std::cout <<jentry<<" events processed" << std::endl;}
   
	
	// Condition to avoid double counting in data //
    bool mu_trig = Muon_trig_pass;
    bool el_trig = Electron_trig_pass;	
	bool jet_trig = hlt_AK8PFJet500 || hlt_PFJet500;
    bool emu_trig = hlt_Mu37_Ele27_CaloIdL_MW || hlt_Mu27_Ele37_CaloIdL_MW;
    
    //Muon trigger selection
	if(!mu_trig) {continue;}

	//Offline muon selection
	if(fabs(l1_pdgId) != 13 ) {continue;} //only muon selection
	if(l1_minisoall > miniso_cut) {continue;} //muon isolation
	
	//At least one AK8 jet
    if(nPFJetAK8<1) continue;
    
    TLorentzVector muon_cand;
    muon_cand.SetPtEtaPhiM(l1_pt,l1_eta,l1_phi,l1_mass);	
     	
    TLorentzVector J_cand;
    J_cand.SetPtEtaPhiM(PFJetAK8_pt[0],PFJetAK8_eta[0],PFJetAK8_phi[0],PFJetAK8_mass[0]);
    
    //DR(muon_cand,Y_cand) > 1.2 application to be opposite direction
	Float_t DR_muon_J;
	DR_muon_J = muon_cand.DeltaR(J_cand);
	if(DR_muon_J < 1.2) {continue;}
        	
    // Selections so that regions where ParticleNet bb tagger is not trianed or SFs are not derived are not used //
	if(PFJetAK8_msoftdrop[0] < msd_cut) continue;
		
	// MET condition
	if(MET_pt < 30) continue;
	
	// read number of JES+JER uncs. //
	njecmax = ((*PFJetAK8_JESup_split)[0]).size();

    //Trigger_pt Scale factor determination and the corresponding unc 
    double SF_Trig, SF_Trig_stat, SF_Trig_syst;
    SF_Trig = 1; SF_Trig_stat = 1; SF_Trig_syst = 1;
    if(!isDATA) {
		if(mu_trig && abs(l1_pdgId)==13) {
			if(l1_pt > 50)
		    {
				int etabin = std::max(1, std::min(h_singlemuon->GetNbinsY(), h_singlemuon->GetYaxis()->FindBin(fabs(l1_eta))));
                int ptbin  = std::max(1, std::min(h_singlemuon->GetNbinsX(), h_singlemuon->GetXaxis()->FindBin(l1_pt)));
                SF_Trig = TMath::Max(float(1.e-3),float(h_singlemuon->GetBinContent(ptbin,etabin))) ;
				SF_Trig_stat = 0.0 ;	
				SF_Trig_syst = std::max(TMath::Max(float(1.e-3),float(h_singlemuon_syst_up->GetBinContent(ptbin,etabin))), TMath::Max(float(1.e-3),float(h_singlemuon_syst_dn->GetBinContent(ptbin,etabin)))) ;
            } // l1_pt threshold
		} //muon trigger
	} //Only MC	

	// # of b-tagged jets outside leading AK8 jet
	
	int nbjets_outJ = 0;
		
	for( int jj = 0; jj < nJetAK4 ; jj++)
	{
		if(delta2R(JetAK4_eta[jj],JetAK4_phi[jj],J_cand.Eta(),J_cand.Phi())>1.2){
			if(JetAK4_btag_DeepFlav[jj] > deep_btag_cut ){
				nbjets_outJ++;
			}
		}
	}
        
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
		//b-tagging SF for the Ak4 jets	 
		for( int jj = 0; jj < nJetAK4 ; jj++)
		{
		if(delta2R(JetAK4_eta[jj],JetAK4_phi[jj],J_cand.Eta(),J_cand.Phi())>1.2){	
					
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
	
	int jk_b = -1;
	jk_b = (nbjets_outJ==0)?1:2;
	
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
		weight_nom *= SF_Trig; //Trigger scale factors 
	 }
  
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
 
	bool Flag_AK8_bb_pass = (PFJetAK8_DeepTag_PNetMD_XbbvsQCD[0]>=PNetbb_cut_T)?true:false;
 
	int iptbin = getbinid(PFJetAK8_pt[0],nak8ptbins,ak8ptbins);
	
	if(iptbin<0||iptbin>=nak8ptbins) continue;
     
    float weight = weight_nom; 
    
     //Looping over nbjet
    
	for(int jk=0; jk<nbcat; jk++){
		 
		if(!(jk==0 || jk==jk_b)) continue;
		
		if(jk!=0) { weight = weight_nom*b_SF; } // applying b tagging SF only if any condition on number of b-tagged jets is used
		
		if(Flag_AK8_bb_pass)
        {
			h_AK8_msoftdrop_pass[jk][iptbin]->Fill(PFJetAK8_msoftdrop[0],weight);
	  	}
		else
		{
            h_AK8_msoftdrop_fail[jk][iptbin]->Fill(PFJetAK8_msoftdrop[0],weight);
		}
            
        if(!isDATA) 
		{
			vector<float> shape_weight_up, shape_weight_dn, shape_weight_nom;
					
			shape_weight_up.push_back(puWeightup);		        shape_weight_dn.push_back(puWeightdown);  		shape_weight_nom.push_back(puWeight);
			shape_weight_up.push_back(leptonsf_weight_stat);	shape_weight_dn.push_back(leptonsf_weight_syst);	shape_weight_nom.push_back(leptonsf_weight);
			shape_weight_up.push_back(leptonsf_weight_up);		shape_weight_dn.push_back(leptonsf_weight_dn);		shape_weight_nom.push_back(leptonsf_weight);
			shape_weight_up.push_back(prefiringweightup);		shape_weight_dn.push_back(prefiringweightdown);		shape_weight_nom.push_back(prefiringweight);
			shape_weight_up.push_back(bb_SF_up);			shape_weight_dn.push_back(bb_SF_dn);			shape_weight_nom.push_back(bb_SF);
			if(!isDL){
				shape_weight_up.push_back(W_SF_up);			shape_weight_dn.push_back(W_SF_dn); 			shape_weight_nom.push_back(W_SF); // for single-lepton channel only
			}
			else{
				shape_weight_up.push_back(1);				shape_weight_dn.push_back(1); 				shape_weight_nom.push_back(1); 
			}
				
			if(jk!=0){
				shape_weight_up.push_back(b_SF_up);			shape_weight_dn.push_back(b_SF_dn); 			shape_weight_nom.push_back(b_SF); // only with condition on # of b jets
			}
			else{
				shape_weight_up.push_back(1);					shape_weight_dn.push_back(1); 				shape_weight_nom.push_back(1);
			}
				
			shape_weight_up.push_back(SF_Trig_1_up);		shape_weight_dn.push_back(SF_Trig_1_dn);			shape_weight_nom.push_back(SF_Trig);
			shape_weight_up.push_back(SF_Trig_2_up);		shape_weight_dn.push_back(SF_Trig_2_dn);			shape_weight_nom.push_back(SF_Trig);
					
			bool top_fullymerged = (PFJetAK8_label_Top_bqq[0] || PFJetAK8_label_Top_bcq[0]);
			bool top_semimerged = (PFJetAK8_label_Top_bq[0] || PFJetAK8_label_Top_bc[0] || PFJetAK8_label_W_qq[0] || PFJetAK8_label_W_cq[0]);
			bool W_merged = (PFJetAK8_label_W_qq[0] || PFJetAK8_label_W_cq[0]);

            //Y_DeepTag_PNetMD_XbbvsQCD pass & fail histogram
            if(Flag_AK8_bb_pass)
			{
				//SET-A
				if(top_fullymerged)      {h_AK8_msoftdrop_sys_Top_fullymerged_pass[jk][iptbin][0]->Fill(PFJetAK8_msoftdrop[0],weight);} 
				else if (top_semimerged) {h_AK8_msoftdrop_sys_Top_semimerged_pass[jk][iptbin][0]->Fill(PFJetAK8_msoftdrop[0],weight);}
                else                     {h_AK8_msoftdrop_sys_rest_pass[jk][iptbin][0]->Fill(PFJetAK8_msoftdrop[0],weight);}                                             
                //SET-B
                if(top_fullymerged) {h_AK8_msoftdrop_sys_Top_fullymerged_v2_pass[jk][iptbin][0]->Fill(PFJetAK8_msoftdrop[0],weight);}
				else if(W_merged)   {h_AK8_msoftdrop_sys_W_hadronic_pass[jk][iptbin][0]->Fill(PFJetAK8_msoftdrop[0],weight);}
				else                {h_AK8_msoftdrop_sys_rest_v2_pass[jk][iptbin][0]->Fill(PFJetAK8_msoftdrop[0],weight);}
			}
            else
			{
				//SET-A
                if(top_fullymerged)     {h_AK8_msoftdrop_sys_Top_fullymerged_fail[jk][iptbin][0]->Fill(PFJetAK8_msoftdrop[0],weight);}
                else if(top_semimerged) {h_AK8_msoftdrop_sys_Top_semimerged_fail[jk][iptbin][0]->Fill(PFJetAK8_msoftdrop[0],weight);}
                else                    {h_AK8_msoftdrop_sys_rest_fail[jk][iptbin][0]->Fill(PFJetAK8_msoftdrop[0],weight);}
                //SET-B
                if(top_fullymerged) 	{h_AK8_msoftdrop_sys_Top_fullymerged_v2_fail[jk][iptbin][0]->Fill(PFJetAK8_msoftdrop[0],weight);}
                else if(W_merged)  	{h_AK8_msoftdrop_sys_W_hadronic_fail[jk][iptbin][0]->Fill(PFJetAK8_msoftdrop[0],weight);}
                else                 {h_AK8_msoftdrop_sys_rest_v2_fail[jk][iptbin][0]->Fill(PFJetAK8_msoftdrop[0],weight);}
			}
   				
   			for(int isys=0; isys<nsys; isys++){
				if(isys==0 || isys<njecmax) // JES + JER (splitted in sources)
		 		{
					
					int iptbin_up = getbinid(PFJetAK8_pt[0]*(*PFJetAK8_JESup_split)[0][isys],nak8ptbins,ak8ptbins);
					int iptbin_dn = getbinid(PFJetAK8_pt[0]*(*PFJetAK8_JESdn_split)[0][isys],nak8ptbins,ak8ptbins);
					
					if(iptbin_up<0||iptbin_dn<0) continue;
					
					if(Flag_AK8_bb_pass)
                    {
						//SET-A
						if(top_fullymerged) 
						{ 
							h_AK8_msoftdrop_sys_Top_fullymerged_pass[jk][iptbin_up][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESup_split)[isys],weight);
							h_AK8_msoftdrop_sys_Top_fullymerged_pass[jk][iptbin_dn][2*(isys+1)]  ->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESdn_split)[isys],weight); 
						}	 
                        else if(top_semimerged) 
						{ 
							h_AK8_msoftdrop_sys_Top_semimerged_pass[jk][iptbin_up][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESup_split)[isys],weight);
							h_AK8_msoftdrop_sys_Top_semimerged_pass[jk][iptbin_dn][2*(isys+1)]  ->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESdn_split)[isys],weight); 
						}	
                        else            
						{ 
							h_AK8_msoftdrop_sys_rest_pass[jk][iptbin_up][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESup_split)[isys],weight);
							h_AK8_msoftdrop_sys_rest_pass[jk][iptbin_dn][2*(isys+1)]  ->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESdn_split)[isys],weight); 
						}
                        //SET-B
                        if(top_fullymerged) 
						{ 
							h_AK8_msoftdrop_sys_Top_fullymerged_v2_pass[jk][iptbin_up][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESup_split)[isys],weight);
                            h_AK8_msoftdrop_sys_Top_fullymerged_v2_pass[jk][iptbin_dn][2*(isys+1)]  ->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESdn_split)[isys],weight); 
                        }
                        else if(W_merged)
						{ 
							h_AK8_msoftdrop_sys_W_hadronic_pass[jk][iptbin_up][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESup_split)[isys],weight);
                            h_AK8_msoftdrop_sys_W_hadronic_pass[jk][iptbin_dn][2*(isys+1)]  ->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESdn_split)[isys],weight); 
                        }
                        else 
                        { 
							 h_AK8_msoftdrop_sys_rest_v2_pass[jk][iptbin_up][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESup_split)[isys],weight);
                             h_AK8_msoftdrop_sys_rest_v2_pass[jk][iptbin_dn][2*(isys+1)]  ->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESdn_split)[isys],weight);  
                        }
					}
                    else
                    {
						//SET-A
                        if(top_fullymerged)
                        { 
							h_AK8_msoftdrop_sys_Top_fullymerged_fail[jk][iptbin_up][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESup_split)[isys],weight);
                            h_AK8_msoftdrop_sys_Top_fullymerged_fail[jk][iptbin_dn][2*(isys+1)]  ->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESdn_split)[isys],weight); 
                        }
                        else if(top_semimerged)
                        { 
							h_AK8_msoftdrop_sys_Top_semimerged_fail[jk][iptbin_up][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESup_split)[isys],weight);
                            h_AK8_msoftdrop_sys_Top_semimerged_fail[jk][iptbin_dn][2*(isys+1)]  ->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESdn_split)[isys],weight); 
                        }
                        else
                        { 
							h_AK8_msoftdrop_sys_rest_fail[jk][iptbin_up][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESup_split)[isys],weight);
                            h_AK8_msoftdrop_sys_rest_fail[jk][iptbin_dn][2*(isys+1)]  ->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESdn_split)[isys],weight); 
                        }
                        //SET-B
                        if(top_fullymerged) 
						{  
							h_AK8_msoftdrop_sys_Top_fullymerged_v2_fail[jk][iptbin_up][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESup_split)[isys],weight);
                            h_AK8_msoftdrop_sys_Top_fullymerged_v2_fail[jk][iptbin_dn][2*(isys+1)]  ->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESdn_split)[isys],weight); 
                        }
                        else if(W_merged)
                        { 
							h_AK8_msoftdrop_sys_W_hadronic_fail[jk][iptbin_up][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESup_split)[isys],weight);
                            h_AK8_msoftdrop_sys_W_hadronic_fail[jk][iptbin_dn][2*(isys+1)]  ->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESdn_split)[isys],weight); 
                        }
                        else
                        { 
							h_AK8_msoftdrop_sys_rest_v2_fail[jk][iptbin_up][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESup_split)[isys],weight);
                            h_AK8_msoftdrop_sys_rest_v2_fail[jk][iptbin_dn][2*(isys+1)]->Fill(PFJetAK8_msoftdrop[0],weight);//*(*PFJetAK8_JESdn_split)[isys],weight);  
                        }
                            
					} 
				}

				//cout << "check! dude-5" << endl;
				if(isys>=njecmax) // Weight-based systematics (SFs)
				{
							
					int isys1 = isys-(njecmax);
							
					float weight_up = weight*shape_weight_up[isys1]/TMath::Max(float(1.e-6),shape_weight_nom[isys1]);
					float weight_dn = weight*shape_weight_dn[isys1]/TMath::Max(float(1.e-6),shape_weight_nom[isys1]);
						
					if(jk==0 && isys1==4) {  weight_up = 1.; weight_dn = 1.; } // b tagging SF
					if(Flag_AK8_bb_pass)
                    {
						//SET-A
						if(top_fullymerged)
                        { 
							h_AK8_msoftdrop_sys_Top_fullymerged_pass[jk][iptbin][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight_up);
                            h_AK8_msoftdrop_sys_Top_fullymerged_pass[jk][iptbin][2*(isys+1)]  ->Fill(PFJetAK8_msoftdrop[0],weight_dn); 
                        }
                        else if(top_semimerged)
                        { 
							h_AK8_msoftdrop_sys_Top_semimerged_pass[jk][iptbin][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight_up);
                            h_AK8_msoftdrop_sys_Top_semimerged_pass[jk][iptbin][2*(isys+1)]   ->Fill(PFJetAK8_msoftdrop[0],weight_dn); 
                        }
                        else
                        { 
							h_AK8_msoftdrop_sys_rest_pass[jk][iptbin][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight_up);
                            h_AK8_msoftdrop_sys_rest_pass[jk][iptbin][2*(isys+1)]->Fill(PFJetAK8_msoftdrop[0],weight_dn); 
                        }
                        //SET-B
                        if(top_fullymerged)
						{
							h_AK8_msoftdrop_sys_Top_fullymerged_v2_pass[jk][iptbin][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight_up);
                            h_AK8_msoftdrop_sys_Top_fullymerged_v2_pass[jk][iptbin][2*(isys+1)]  ->Fill(PFJetAK8_msoftdrop[0],weight_dn); 
                        }
                        else if(W_merged)
                        { 
							h_AK8_msoftdrop_sys_W_hadronic_pass[jk][iptbin][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight_up);
                            h_AK8_msoftdrop_sys_W_hadronic_pass[jk][iptbin][2*(isys+1)]   ->Fill(PFJetAK8_msoftdrop[0],weight_dn); 
                        }
                        else
                        { 
							h_AK8_msoftdrop_sys_rest_v2_pass[jk][iptbin][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight_up);
                            h_AK8_msoftdrop_sys_rest_v2_pass[jk][iptbin][2*(isys+1)]->Fill(PFJetAK8_msoftdrop[0],weight_dn);  
                        }
                    }
                    else
                    {
						//SET-A
                        if(top_fullymerged)
                        {
							h_AK8_msoftdrop_sys_Top_fullymerged_fail[jk][iptbin][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight_up);
                            h_AK8_msoftdrop_sys_Top_fullymerged_fail[jk][iptbin][2*(isys+1)]  ->Fill(PFJetAK8_msoftdrop[0],weight_dn); 
                        }
                        else if(top_semimerged)
                        { 
							h_AK8_msoftdrop_sys_Top_semimerged_fail[jk][iptbin][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight_up);
                            h_AK8_msoftdrop_sys_Top_semimerged_fail[jk][iptbin][2*(isys+1)]   ->Fill(PFJetAK8_msoftdrop[0],weight_dn); 
                        }
                        else
                        { 
							h_AK8_msoftdrop_sys_rest_fail[jk][iptbin][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight_up);
                             h_AK8_msoftdrop_sys_rest_fail[jk][iptbin][2*(isys+1)]->Fill(PFJetAK8_msoftdrop[0],weight_dn); 
                        }
                        //SET-B
                        if(top_fullymerged) 
						{ 
							h_AK8_msoftdrop_sys_Top_fullymerged_v2_fail[jk][iptbin][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight_up);
                            h_AK8_msoftdrop_sys_Top_fullymerged_v2_fail[jk][iptbin][2*(isys+1)]  ->Fill(PFJetAK8_msoftdrop[0],weight_dn); 
                        }
                        else if(W_merged)
                        { 
							h_AK8_msoftdrop_sys_W_hadronic_fail[jk][iptbin][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight_up);
                            h_AK8_msoftdrop_sys_W_hadronic_fail[jk][iptbin][2*(isys+1)]   ->Fill(PFJetAK8_msoftdrop[0],weight_dn); 
                        }
                        else
                        { 
							h_AK8_msoftdrop_sys_rest_v2_fail[jk][iptbin][2*(isys+1)-1]->Fill(PFJetAK8_msoftdrop[0],weight_up);
                            h_AK8_msoftdrop_sys_rest_v2_fail[jk][iptbin][2*(isys+1)]->Fill(PFJetAK8_msoftdrop[0],weight_dn);  
                        }

                    }
				}
					   //std::cout << "Dude-6" << std::endl;	
				shape_weight_up.clear();
				shape_weight_dn.clear();
				shape_weight_nom.clear();
				
				}//sys loop
			} //No Data condition	
				
		}// b cat (jk)
	   
	
	
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

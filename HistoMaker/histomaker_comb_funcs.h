#include "histomaker_comb.h"

bool check_double_count(TString inputFile, bool mu_trig, bool el_trig, bool jet_trig, bool emu_trig)
{

	bool pass_doublecount = true;

	if (string(inputFile.Data()).find("SingleMuon")!=string::npos)
	{     
		if(!mu_trig)  { pass_doublecount = false; }               
	}
    else if (string(inputFile.Data()).find("DoubleMuon")!=string::npos)
    {
		if(!mu_trig)  { pass_doublecount = false; }     
    }
	else if (string(inputFile.Data()).find("EGamma")!=string::npos)
	{
		if(mu_trig || !el_trig)  { pass_doublecount = false; }     
	}
	else if (string(inputFile.Data()).find("JetHT")!=string::npos)
	{
		if(mu_trig || el_trig || !jet_trig)  { pass_doublecount = false; }    
	}
    else if (string(inputFile.Data()).find("MuonEG")!=string::npos)
    {
		if(mu_trig || el_trig || !emu_trig)  { pass_doublecount = false; }    
    }
	else{
		pass_doublecount = false;
	}
		
	return pass_doublecount;
	
}

float* get_Top_SF(float pt, bool gen_match, bool score_pass)
{
	
	static float SFs[3];
	SFs[0] = 1.;
	SFs[1] = 1.;
	SFs[2] = 1.;
	
	int pt_bin = getbinid(pt,PNTop_SF_nptbins,PNTop_SF_ptbins);
	if(pt_bin>=0 && pt_bin<PNTop_SF_nptbins) { 
		// use only GEN-matched top for applying SFs 
		if(gen_match){
 			if(score_pass){
				SFs[0] = PNTop_SF_M[pt_bin]; 
				SFs[1] = PNTop_SF_M_up[pt_bin]; 
				SFs[2] = PNTop_SF_M_dn[pt_bin]; 
			}
		}
	}
	
    return 	SFs;
}

float* get_W_SF(float pt, bool gen_match, bool score_pass)
{
	static float SFs[3];
	SFs[0] = 1.;
	SFs[1] = 1.;
	SFs[2] = 1.;
	
	int pt_bin = getbinid(pt,PNW_SF_nptbins,PNW_SF_ptbins);
	if(pt_bin>=0 && pt_bin<PNW_SF_nptbins) { 
		if(gen_match){
			if(score_pass){
				// we are using SFs for tight WP
				SFs[0] = PNW_SF_T[pt_bin]; 
				SFs[1] = PNW_SF_T_up[pt_bin]; 
				SFs[2] = PNW_SF_T_dn[pt_bin]; 
			}
		}
	}
	
	return SFs;
}

float *get_bb_SF(float pt, bool gen_match, bool score_pass)
{
	static float SFs[3];
	SFs[0] = 1.;
	SFs[1] = 1.;
	SFs[2] = 1.;
	
	int pt_bin = getbinid(pt,PNbb_SF_nptbins,PNbb_SF_ptbins);
	if(pt_bin>=0 && pt_bin<PNbb_SF_nptbins) {
		if(gen_match){
			if(score_pass)
            {
				SFs[0] = PNbb_SF_HP[pt_bin]; 
				SFs[1] = PNbb_SF_HP_up[pt_bin]; 
				SFs[2] = PNbb_SF_HP_dn[pt_bin]; 
            }
		}
	}
	
	return SFs;
}

float* get_AK8_tagging_SF(float pt, bool gen_match, bool score_pass, string tagging)
{
	static float SFs[3];
	SFs[0] = 1.;
	SFs[1] = 1.;
	SFs[2] = 1.;
	
	if(gen_match){
		if(score_pass){
	
			if(tagging.find("W_PN")!=string::npos){
	
				int pt_bin = getbinid(pt,PNW_SF_nptbins,PNW_SF_ptbins);
				if(pt_bin>=0 && pt_bin<PNW_SF_nptbins) { 
					// we are using SFs for tight WP
					SFs[0] = PNW_SF_T[pt_bin]; 
					SFs[1] = PNW_SF_T_up[pt_bin]; 
					SFs[2] = PNW_SF_T_dn[pt_bin]; 
				}
			}
	
			else if(tagging.find("Xbb_PN")!=string::npos){
				
				int pt_bin = getbinid(pt,PNbb_SF_nptbins,PNbb_SF_ptbins);
				if(pt_bin>=0 && pt_bin<PNbb_SF_nptbins) {
					SFs[0] = PNbb_SF_HP[pt_bin]; 
					SFs[1] = PNbb_SF_HP_up[pt_bin]; 
					SFs[2] = PNbb_SF_HP_dn[pt_bin]; 
				}
			}
			
			else if(tagging.find("Top_PN")!=string::npos){
				
				int pt_bin = getbinid(pt,PNTop_SF_nptbins,PNTop_SF_ptbins);
				if(pt_bin>=0 && pt_bin<PNTop_SF_nptbins) { 
					SFs[0] = PNTop_SF_M[pt_bin]; 
					SFs[1] = PNTop_SF_M_up[pt_bin]; 
					SFs[2] = PNTop_SF_M_dn[pt_bin]; 
				}
			}
	
		}
	}
	
	return SFs;
}

float *get_btag_SF( float pt, float eta, 
					int hadronflav, bool btag_pass, 
					float btag_SF, float btag_SF_up, float btag_SF_dn, 
				    TH2F *h_AK4_flv_b_eff, TH2F *h_AK4_flv_c_eff, TH2F *h_AK4_flv_l_eff)
{
	static float SFs[3];
	SFs[0] = 1.;
	SFs[1] = 1.;
	SFs[2] = 1.;
	
	int etabin = std::max(1, std::min(h_AK4_flv_b_eff->GetNbinsY(), h_AK4_flv_b_eff->GetYaxis()->FindBin(fabs(eta))));
	int ptbin  = std::max(1, std::min(h_AK4_flv_b_eff->GetNbinsX(), h_AK4_flv_b_eff->GetXaxis()->FindBin(pt)));
	
	if(btag_pass)
	{
		SFs[0] = btag_SF;
		SFs[1] = btag_SF_up; 
		SFs[2] = btag_SF_dn;
	}
	else
	{
		if ( abs(hadronflav) == 5 ) {
			SFs[0] = std::max(1.e-6,(1.-btag_SF*TMath::Max(float(1.e-3),float(h_AK4_flv_b_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4_flv_b_eff->GetBinContent(ptbin, etabin))) ));
            SFs[1] = std::max(1.e-6,(1.-btag_SF_up*TMath::Max(float(1.e-3),float(h_AK4_flv_b_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4_flv_b_eff->GetBinContent(ptbin, etabin))) ));
			SFs[2] = std::max(1.e-6,(1.-btag_SF_dn*TMath::Max(float(1.e-3),float(h_AK4_flv_b_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4_flv_b_eff->GetBinContent(ptbin, etabin))) ));
        }
        else if ( abs(hadronflav) == 4 ) {
			SFs[0] = std::max(1.e-6,(1.-btag_SF*TMath::Max(float(1.e-3),float(h_AK4_flv_c_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4_flv_c_eff->GetBinContent(ptbin, etabin))) ));
            SFs[1] = std::max(1.e-6,(1.-btag_SF_up*TMath::Max(float(1.e-3),float(h_AK4_flv_c_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4_flv_c_eff->GetBinContent(ptbin, etabin))) ));
			SFs[2] = std::max(1.e-6,(1.-btag_SF_dn*TMath::Max(float(1.e-3),float(h_AK4_flv_c_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4_flv_c_eff->GetBinContent(ptbin, etabin))) ));
        }
		else {
            SFs[0] = std::max(1.e-6,(1.-btag_SF*TMath::Max(float(1.e-3),float(h_AK4_flv_l_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4_flv_l_eff->GetBinContent(ptbin, etabin)))));
            SFs[1] = std::max(1.e-6,(1.-btag_SF_up*TMath::Max(float(1.e-3),float(h_AK4_flv_l_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4_flv_l_eff->GetBinContent(ptbin, etabin)))));
            SFs[2] = std::max(1.e-6,(1.-btag_SF_dn*TMath::Max(float(1.e-3),float(h_AK4_flv_l_eff->GetBinContent(ptbin, etabin))))/ ( 1.- TMath::Max(float(1.e-3),float(h_AK4_flv_l_eff->GetBinContent(ptbin, etabin)))));
        }
    }
    
    return SFs;
}

float *get_trigger_SF(bool isDL, 
					  bool mu_trig, bool el_trig, bool jet_trig, bool emu_trig,
					  float l1_pt, float l1_eta, int l1_pdgId, float l2_pt, float l2_eta, int l2_pdgId, float jet_pt, float jet_eta,
					  TH2F *h_singlemuon, TH2F *h_singlemuon_syst_up, TH2F *h_singlemuon_syst_dn,
					  TH2F *h_singleele, TH2F *h_singleele_stat, TH2F *h_singleele_syst,
					  TH2F *h_jet, TH2F *h_jet_stat, TH2F *h_jet_syst, 
					  TH2F *h_trigmumu, TH2F *h_trigee, TH2F *h_trigemu
					)
{

float SF = 1;
float SF_stat = 0;
float SF_syst = 0;

if(!isDL) {
	if(mu_trig && abs(l1_pdgId)==13) {
		if(l1_pt > mu_trig_pt_SL)
		{
			int etabin = std::max(1, std::min(h_singlemuon->GetNbinsY(), h_singlemuon->GetYaxis()->FindBin(fabs(l1_eta))));
			int ptbin  = std::max(1, std::min(h_singlemuon->GetNbinsX(), h_singlemuon->GetXaxis()->FindBin(l1_pt)));
			SF = TMath::Max(float(1.e-3),float(h_singlemuon->GetBinContent(ptbin,etabin))) ;
			SF_stat = 0.0 ;	
			SF_syst = std::max(TMath::Max(float(1.e-3),float(h_singlemuon_syst_up->GetBinContent(ptbin,etabin))), TMath::Max(float(1.e-3),float(h_singlemuon_syst_dn->GetBinContent(ptbin,etabin)))) ;
		} // l1_pt threshold
    } //muon trigger
	else if ( !mu_trig  && el_trig &&  abs(l1_pdgId)==11) {
		if(l1_pt > el_trig_pt_SL )
		{
			int etabin = std::max(1, std::min(h_singleele->GetNbinsY(), h_singleele->GetYaxis()->FindBin(fabs(l1_eta))));
			int ptbin  = std::max(1, std::min(h_singleele->GetNbinsX(), h_singleele->GetXaxis()->FindBin(l1_pt)));
			SF      = TMath::Max(float(1.e-3),float(h_singleele->GetBinContent(ptbin, etabin))) ;
			SF_stat = TMath::Max(float(1.e-3),float(h_singleele_stat->GetBinContent(ptbin, etabin))) ;
			SF_syst = TMath::Max(float(1.e-3),float(h_singleele_syst->GetBinContent(ptbin, etabin))) ;
		} //l1_pt threshold
	} //electron trigger
    else if ( !mu_trig && !el_trig && jet_trig ){
		if( jet_pt > jet_trig_pt_SL  )
        {
			int etabin = std::max(1, std::min(h_jet->GetNbinsY(), h_jet->GetYaxis()->FindBin(fabs(jet_eta))));
            int ptbin  = std::max(1, std::min(h_jet->GetNbinsX(), h_jet->GetXaxis()->FindBin(jet_pt)));
            SF      = TMath::Max(float(1.e-3),float(h_jet->GetBinContent( ptbin, etabin))) ;
            SF_stat = TMath::Max(float(1.e-3),float(h_jet_stat->GetBinContent(ptbin, etabin))) ;
            SF_syst = TMath::Max(float(1.e-3),float(h_jet_syst->GetBinContent(ptbin, etabin))) ;
        }//jet pt thresholds
    } //jet trig
} //Semi-lepton
else{
	if( mu_trig && abs(l1_pdgId)==13 && abs(l2_pdgId)==13)
	{
		if(l1_pt > mu_trig_pt_DL_0 && l2_pt > mu_trig_pt_DL_1  )
		{
			int pt1bin = std::max(1, std::min(h_trigmumu->GetNbinsX(), h_trigmumu->GetXaxis()->FindBin(l1_pt)));
			int pt2bin = std::max(1, std::min(h_trigmumu->GetNbinsY(), h_trigmumu->GetYaxis()->FindBin(l2_pt)));
			SF        = TMath::Max(float(1.e-3),float(h_trigmumu->GetBinContent(pt1bin, pt2bin))) ;	
			SF_stat    = 0.0;
			SF_syst    = TMath::Max(float(1.e-3),float(h_trigmumu->GetBinError(pt1bin, pt2bin))) ;
        }
   	} //mu-mu trigger
    else if( !mu_trig && el_trig && abs(l1_pdgId)==11 && abs(l2_pdgId)==11)
    {
		if(l1_pt > el_trig_pt_DL && l2_pt > el_trig_pt_DL )
		{
			int pt1bin = std::max(1, std::min(h_trigee->GetNbinsX(), h_trigee->GetXaxis()->FindBin(l1_pt)));
            int pt2bin = std::max(1, std::min(h_trigee->GetNbinsY(), h_trigee->GetYaxis()->FindBin(l2_pt)));
            SF         = TMath::Max(float(1.e-3),float(h_trigee->GetBinContent(pt1bin, pt2bin))) ;
			SF_stat    = 0.0 ;
			SF_syst    = TMath::Max(float(1.e-3),float(h_trigee->GetBinError(pt1bin, pt2bin))) ;
        }
   	} // e-e trigger
    else if(!mu_trig && !el_trig &&  emu_trig && ((abs(l1_pdgId)==11 && abs(l2_pdgId)==13) || (abs(l1_pdgId)==13 && abs(l2_pdgId)==11)) )
	{
		if(l1_pt > emu_trig_pt_DL_0 && l2_pt > emu_trig_pt_DL_1)
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
			SF         = TMath::Max(float(1.e-3),float(h_trigemu->GetBinContent(ptelebin, ptmubin))) ;
			SF_stat    = 0.0 ;
			SF_syst    = TMath::Max(float(1.e-3),float(h_trigemu->GetBinError(ptelebin, ptmubin))) ;
		 }
	}  //e-mu trigger
} //di-lepton	
	
	static float sfvalues[3];
	sfvalues[0] = SF;
	sfvalues[1] = SF_stat;
	sfvalues[2] = SF_syst;
	
	return sfvalues;
}

bool check_Y_genmatching(TLorentzVector Y_cand, vector<TLorentzVector> genbs, vector<int> genpdgids, vector<int> genmompdgids)
{

bool isMatched = false;

for(int pp = 0 ; pp  < (int(genbs.size())-1); pp++)
    {
		for(int qq = pp+1; qq < int(genbs.size()); qq++)
        {
			TLorentzVector gen_b1, gen_b2;
			gen_b1 = genbs[pp];
			gen_b2 = genbs[qq];
            Float_t DR_Yb1, DR_Yb2;
            DR_Yb1 = Y_cand.DeltaR(gen_b1);
            DR_Yb2 = Y_cand.DeltaR(gen_b2);
            if(DR_Yb1 < 0.8 && DR_Yb2 < 0.8 && genpdgids[pp]*genpdgids[qq] < 0 && genmompdgids[pp] == genmompdgids[qq] )
            //GenBPart_pdgId[pp]*GenBPart_pdgId[qq] < 0 && GenBPart_mompdgId[pp] == GenBPart_mompdgId[qq] )
            {
				isMatched = true;
                break;
            }
         }//qq
    }//pp      	
	
}

enum TheRunEra{
  y2016B,y2016C,y2016D,y2016E,y2016F,y2016G,y2016H,
  y2017B,y2017C,y2017D,y2017E,y2017F,
  y2018A,y2018B,y2018C,y2018D,
  y2016MC,
  y2017MC,
  y2018MC,
  yUL2016B,yUL2016C,yUL2016D,yUL2016E,yUL2016F,yUL2016Flate,yUL2016G,yUL2016H,
  yUL2017B,yUL2017C,yUL2017D,yUL2017E,yUL2017F,
  yUL2018A,yUL2018B,yUL2018C,yUL2018D,
  yUL2016MCAPV,
  yUL2016MCnonAPV,
  yUL2017MC,
  yUL2018MC
};


std::pair<double,double> METXYCorr_Met_MetPhi(double uncormet, double uncormet_phi, int runnb, TString year, bool isMC, int npv, bool isUL =true, bool ispuppi=true, TString eraID=""){

  std::pair<double,double>  TheXYCorr_Met_MetPhi(uncormet,uncormet_phi);
  
  if(npv>100) npv=100;
  int runera =-1;
  bool usemetv2 =false;
  
  if(isMC && year == "2016APV" && isUL) runera = yUL2016MCAPV;
  else if(isMC && year == "2016nonAPV" && isUL) runera = yUL2016MCnonAPV;
  else if(isMC && year == "2017" && isUL) runera = yUL2017MC;
  else if(isMC && year == "2018" && isUL) runera = yUL2018MC;
  
  else if(!isMC && runnb >=315252 && runnb <=316995 && isUL) runera = yUL2018A;
  else if(!isMC && runnb >=316998 && runnb <=319312 && isUL) runera = yUL2018B;
  else if(!isMC && runnb >=319313 && runnb <=320393 && isUL) runera = yUL2018C;
  else if(!isMC && runnb >=320394 && runnb <=325273 && isUL) runera = yUL2018D;
  
  else if(!isMC && runnb >=297020 && runnb <=299329 && isUL) runera = yUL2017B;
  else if(!isMC && runnb >=299337 && runnb <=302029 && isUL) runera = yUL2017C; 
  else if(!isMC && runnb >=302030 && runnb <=303434 && isUL) runera = yUL2017D; 
  else if(!isMC && runnb >=303435 && runnb <=304826 && isUL) runera = yUL2017E; 
  else if(!isMC && runnb >=304911 && runnb <=306462 && isUL) runera = yUL2017F; 
  
  else if(!isMC && runnb >=272007 && runnb <=275376 && isUL) runera = yUL2016B;
  else if(!isMC && runnb >=275657 && runnb <=276283 && isUL) runera = yUL2016C;
  else if(!isMC && runnb >=276315 && runnb <=276811 && isUL) runera = yUL2016D;
  else if(!isMC && runnb >=276831 && runnb <=277420 && isUL) runera = yUL2016E;
  else if(!isMC && ((runnb >=277772 && runnb <=278768) || runnb==278770) && isUL) runera = yUL2016F;
  else if(!isMC && ((runnb >=278801 && runnb <=278808) || runnb==278769) && isUL) runera = yUL2016Flate;
  else if(!isMC && runnb >=278820 && runnb <=280385 && isUL) runera = yUL2016G;
  else if(!isMC && runnb >=280919 && runnb <=284044 && isUL) runera = yUL2016H;
  
  /*
  else if(!isMC && year == "2018" && eraID=="A" && isUL) runera = yUL2018A;
  else if(!isMC && year == "2018" && eraID=="B" && isUL) runera = yUL2018B;
  else if(!isMC && year == "2018" && eraID=="C" && isUL) runera = yUL2018C;
  else if(!isMC && year == "2018" && eraID=="D" && isUL) runera = yUL2018D;
  
  else if(!isMC && year == "2017" && eraID=="B" && isUL){ runera = yUL2017B; }
  else if(!isMC && year == "2017" && eraID=="C" && isUL){ runera = yUL2017C; }
  else if(!isMC && year == "2017" && eraID=="D" && isUL){ runera = yUL2017D; }
  else if(!isMC && year == "2017" && eraID=="E" && isUL){ runera = yUL2017E; }
  else if(!isMC && year == "2017" && eraID=="F" && isUL){ runera = yUL2017F; }
  
  else if(!isMC && year == "2016" && eraID=="B" && isUL) runera = yUL2016B;
  else if(!isMC && year == "2016" && eraID=="C" && isUL) runera = yUL2016C;
  else if(!isMC && year == "2016" && eraID=="D" && isUL) runera = yUL2016D;
  else if(!isMC && year == "2016" && eraID=="E" && isUL) runera = yUL2016E;
  else if(!isMC && year == "2016" && eraID=="F" && isUL) runera = yUL2016F;
  else if(!isMC && year == "2016" && eraID=="Flate" && isUL) runera = yUL2016Flate;
  else if(!isMC && year == "2016" && eraID=="G" && isUL) runera = yUL2016G;
  else if(!isMC && year == "2016" && eraID=="H" && isUL) runera = yUL2016H;
  */ 

  else {
    //Couldn't find data/MC era => no correction applied
    //cout<<"could not find era"<<endl;
    return TheXYCorr_Met_MetPhi;
  }
  
  //cout<<"runera "<<runera<<endl;
  
  double METxcorr(0.),METycorr(0.);
  
  if(ispuppi){
   
   //UL2017Puppi
   if(runera==yUL2017B) METxcorr = -(-0.00382117*npv +-0.666228);
   if(runera==yUL2017B) METycorr = -(0.0109034*npv +0.172188);
   if(runera==yUL2017C) METxcorr = -(-0.00110699*npv +-0.747643);
   if(runera==yUL2017C) METycorr = -(-0.0012184*npv +0.303817);
   if(runera==yUL2017D) METxcorr = -(-0.00141442*npv +-0.721382);
   if(runera==yUL2017D) METycorr = -(-0.0011873*npv +0.21646);
   if(runera==yUL2017E) METxcorr = -(0.00593859*npv +-0.851999);
   if(runera==yUL2017E) METycorr = -(-0.00754254*npv +0.245956);
   if(runera==yUL2017F) METxcorr = -(0.00765682*npv +-0.945001);
   if(runera==yUL2017F) METycorr = -(-0.0154974*npv +0.804176);
   if(runera==yUL2017MC) METxcorr = -(-0.0102265*npv +-0.446416);
   if(runera==yUL2017MC) METycorr = -(0.0198663*npv +0.243182);

   //UL2018Puppi
   if(runera==yUL2018A) METxcorr = -(-0.0073377*npv +0.0250294);
   if(runera==yUL2018A) METycorr = -(-0.000406059*npv +0.0417346);
   if(runera==yUL2018B) METxcorr = -(0.00434261*npv +0.00892927);
   if(runera==yUL2018B) METycorr = -(0.00234695*npv +0.20381);
   if(runera==yUL2018C) METxcorr = -(0.00198311*npv +0.37026);
   if(runera==yUL2018C) METycorr = -(-0.016127*npv +0.402029);
   if(runera==yUL2018D) METxcorr = -(0.00220647*npv +0.378141);
   if(runera==yUL2018D) METycorr = -(-0.0160244*npv +0.471053);
   if(runera==yUL2018MC) METxcorr = -(-0.0214557*npv +0.969428);
   if(runera==yUL2018MC) METycorr = -(0.0167134*npv +0.199296);

   //UL2016Puppi
   if(runera==yUL2016B) METxcorr = -(-0.00109025*npv +-0.338093);
   if(runera==yUL2016B) METycorr = -(-0.00356058*npv +0.128407);
   if(runera==yUL2016C) METxcorr = -(-0.00271913*npv +-0.342268);
   if(runera==yUL2016C) METycorr = -(0.00187386*npv +0.104);
   if(runera==yUL2016D) METxcorr = -(-0.00254194*npv +-0.305264);
   if(runera==yUL2016D) METycorr = -(-0.00177408*npv +0.164639);
   if(runera==yUL2016E) METxcorr = -(-0.00358835*npv +-0.225435);
   if(runera==yUL2016E) METycorr = -(-0.000444268*npv +0.180479);
   if(runera==yUL2016F) METxcorr = -(0.0056759*npv +-0.454101);
   if(runera==yUL2016F) METycorr = -(-0.00962707*npv +0.35731);
   if(runera==yUL2016Flate) METxcorr = -(0.0234421*npv +-0.371298);
   if(runera==yUL2016Flate) METycorr = -(-0.00997438*npv +0.0809178);
   if(runera==yUL2016G) METxcorr = -(0.0182134*npv +-0.335786);
   if(runera==yUL2016G) METycorr = -(-0.0063338*npv +0.093349);
   if(runera==yUL2016H) METxcorr = -(0.015702*npv +-0.340832);
   if(runera==yUL2016H) METycorr = -(-0.00544957*npv +0.199093);
   if(runera==yUL2016MCnonAPV) METxcorr = -(-0.0058341*npv +-0.395049);
   if(runera==yUL2016MCnonAPV) METycorr = -(0.00971595*npv +-0.101288);
   if(runera==yUL2016MCAPV) METxcorr = -(-0.0060447*npv +-0.4183);
   if(runera==yUL2016MCAPV) METycorr = -(0.008331*npv +-0.0990046);

  }
  
  double CorrectedMET_x = uncormet *cos( uncormet_phi)+METxcorr;
  double CorrectedMET_y = uncormet *sin( uncormet_phi)+METycorr;

  double CorrectedMET = sqrt(CorrectedMET_x*CorrectedMET_x+CorrectedMET_y*CorrectedMET_y);
  double CorrectedMETPhi;
  if(CorrectedMET_x==0 && CorrectedMET_y>0) CorrectedMETPhi = TMath::Pi();
  else if(CorrectedMET_x==0 && CorrectedMET_y<0 )CorrectedMETPhi = -TMath::Pi();
  else if(CorrectedMET_x >0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x);
  else if(CorrectedMET_x <0&& CorrectedMET_y>0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) + TMath::Pi();
  else if(CorrectedMET_x <0&& CorrectedMET_y<0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) - TMath::Pi();
  else CorrectedMETPhi =0;

  TheXYCorr_Met_MetPhi.first= CorrectedMET;
  TheXYCorr_Met_MetPhi.second= CorrectedMETPhi;
  return TheXYCorr_Met_MetPhi;

}

float get_EWK_cor(TFile *file_SF, float pt, int pdgId, int choice = 2, float pt_min=100)
{
	
	// choice = 2-> taking cor from dark matter paper  
	// any other option-> taking cor from UHH files 
		
	if(choice==0) { pt_min = 150; }
	
	char name[100];
		
	if(choice==2){
		if      (pdgId==23){ sprintf(name,"eej_pTV_kappa_EW"); }
		else if (pdgId==24){ sprintf(name,"evj_pTV_kappa_EW"); }
		else    { sprintf(name,"ewcorr"); }
	}
	else { sprintf(name,"ewcorr"); }
	
	//cout<<name<<endl;
	
	TH1F *h_ewk = (TH1F*)file_SF->Get(name);
	
	float cor = 1;
	if(pt>=pt_min){
		int pt_bin_id = h_ewk->GetXaxis()->FindBin(pt);
		cor = h_ewk->GetBinContent(pt_bin_id);
	}
	
	if(choice==2) { cor += 1; }
		
	return cor;
}

float get_QCD_cor(TFile *file_SF, float pt, int pdgId, int choice = 1, float pt_min=100)
{
	// choice = 2-> taking cor from dark matter paper  
	// choice = 1-> taking cor from CMS NLO samples  
	// any other option-> taking cor from UHH files 
	
	if(choice==0) { pt_min = 150; }
	
	char name[100];
	
	if(choice==1) { sprintf(name,"kFactor"); }
	else if(choice==2) { 
		if      (pdgId==23){ sprintf(name,"eej_pTV_K_NLO"); }
		else if (pdgId==24){ sprintf(name,"evj_pTV_K_NLO"); }
		else { sprintf(name,"kfactor"); }
	}
	else { sprintf(name,"kfactor"); }
	
	//cout<<name<<endl;
	
	TH1F *h_qcd = (TH1F*)file_SF->Get(name);
	float cor = 1;
	if(pt>=pt_min){
		int pt_bin_id = h_qcd->GetXaxis()->FindBin(pt);
		cor = h_qcd->GetBinContent(pt_bin_id);
	}
		
	return cor;
	
}

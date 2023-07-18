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

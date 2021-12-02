#include "Anal_XtoYH.h"
//#include "Functions.h"

void getmuons(std::vector<Muon> &vmuons, float ptcut=25, float etacut=2.5, int id=-1, int maxsize=njetmx, float dxy_cut=0.2, float dz_cut=0.5)
{

for(int mu=0; mu<(nMuon); mu++){
 
    if(Muon_pt[mu]<ptcut) continue; 
    if(fabs(Muon_eta[mu])>etacut)  continue; 

	/*
    bool mu_id = Muon_TightID(Muon_isGL[mu],Muon_isPF[mu],
			      Muon_chi[mu],Muon_hit[mu],Muon_mst[mu],
			      Muon_dxy[mu],Muon_dz[mu],Muon_pixhit[mu],Muon_trklay[mu]);
    if(!mu_id) continue;
    */
    if(id==0 && !Muon_isLoose[mu]) 	continue;
    if(id==1 && !Muon_isMed[mu]) 	continue;
    if(id==2 && !Muon_MediumID[mu]) continue;  //Muon_MediumID variable is actually tight id
    
    //bool mu_iso = Muon_Iso_ID(Muon_pfiso[mu]);
    //if(!mu_iso) continue;
    
    if(fabs(Muon_dxy[mu])>dxy_cut || fabs(Muon_dz[mu])>dz_cut) continue;
    
    Muon vmuon;
    
    vmuon.pt = Muon_pt[mu];
    vmuon.eta = Muon_eta[mu];
    vmuon.phi = Muon_phi[mu];
    vmuon.mass = 0.105658;
    vmuon.charge = Sign(Muon_p[mu]);
    vmuon.p = Muon_p[mu];
    
    vmuon.dxy = Muon_dxy[mu];
    vmuon.dxyErr = Muon_dxyErr[mu];
    vmuon.dz = Muon_dz[mu];
    vmuon.ip3d = Muon_ip3d[mu];
    
    vmuon.isTRK = Muon_isTRK[mu];
    vmuon.isGL = Muon_isGL[mu];
    vmuon.isPF = Muon_isPF[mu];
    vmuon.isLoose = Muon_isLoose[mu];
    vmuon.isGoodGL = Muon_isGoodGL[mu];
    vmuon.isMed = Muon_isMed[mu];
    vmuon.isMedPr = Muon_isMedPr[mu];
    vmuon.isTight = Muon_isTight[mu];
    vmuon.isHighPt = Muon_isHighPt[mu];
    vmuon.isHighPttrk = Muon_isHighPttrk[mu];
    vmuon.MediumID = Muon_MediumID[mu];
    
    vmuon.chi = Muon_chi[mu];
    vmuon.posmatch = Muon_posmatch[mu];
    vmuon.trkink = Muon_trkink[mu];
    vmuon.segcom = Muon_segcom[mu];
    vmuon.hit = Muon_hit[mu];
    vmuon.mst = Muon_mst[mu];
    vmuon.pixhit = Muon_pixhit[mu];
    vmuon.trklay = Muon_trklay[mu];
    vmuon.valfrac = Muon_valfrac[mu];
    vmuon.pfiso = Muon_pfiso[mu];   
   
    vmuon.minisoall = Muon_minisoall[mu];
    vmuon.minchiso = Muon_minchiso[mu];
    vmuon.minnhiso = Muon_minnhiso[mu];
    vmuon.minphiso = Muon_minphiso[mu];
    
    vmuon.p4.SetPtEtaPhiM(vmuon.pt,vmuon.eta,vmuon.phi,vmuon.mass);
    
    vmuons.push_back(vmuon);
    
    if(int(vmuons.size())>=maxsize) break;
  }
  
  sorted_by_pt(vmuons);
	
}

void getelectrons(std::vector<Electron> &velectrons, float ptcut=25, float etacut=2.5, int id=-1, int maxsize=njetmx, float dxy_cut=0.05, float dz_cut=0.1)
{

 for(int ie=0; ie<(nElectron); ie++) {
		
    if (Electron_pt[ie]<ptcut) continue;
    if(fabs(Electron_eta[ie])>etacut)  continue; 
    
    if(id==0 && !Electron_mvaid_Fallv2WP80[ie]) continue;
    if(id==1 && !Electron_mvaid_Fallv2WP80_noIso[ie]) continue;
    if(id==2 && !Electron_mvaid_Fallv2WP90[ie]) continue;
    if(id==3 && !Electron_mvaid_Fallv2WP90_noIso[ie]) continue;
    
    bool impact_pass = 	((fabs(Electron_supcl_eta[ie])<1.4442 && fabs(Electron_dxy[ie])<dxy_cut && fabs(Electron_dz[ie])<dz_cut)
					   ||(fabs(Electron_supcl_eta[ie])>1.5660 && fabs(Electron_dxy[ie])<(2*dxy_cut) && fabs(Electron_dz[ie])<(2*dz_cut)));

	if(!impact_pass) continue;
        
    Electron velectron;

    velectron.pt = Electron_pt[ie];
    velectron.eta = Electron_eta[ie];
    velectron.phi = Electron_phi[ie];
    velectron.mass = 0.000511;
    velectron.p = Electron_p[ie];
    velectron.e_ECAL = Electron_e_ECAL[ie];
    velectron.charge = Sign(Electron_p[ie]);
    
    velectron.Fallv2WP80 = Electron_mvaid_Fallv2WP80[ie];
    velectron.Fallv2WP80_noIso = Electron_mvaid_Fallv2WP80_noIso[ie];
    velectron.Fallv2WP90 = Electron_mvaid_Fallv2WP90[ie];
    velectron.Fallv2WP90_noIso = Electron_mvaid_Fallv2WP90_noIso[ie];
    
    velectron.dxy = Electron_dxy[ie];
    velectron.dxyErr = Electron_dxyErr[ie];
    velectron.dz = Electron_dz[ie];
    velectron.ip3d = Electron_ip3d[ie];
    velectron.dxy_sv = Electron_dxy_sv[ie];
    
    velectron.eccalTrkEnergyPostCorr = Electron_eccalTrkEnergyPostCorr[ie];
    velectron.energyScaleValue = Electron_energyScaleValue[ie];
    velectron.energyScaleUp = Electron_energyScaleUp[ie];
    velectron.energyScaleDown = Electron_energyScaleDown[ie];
    velectron.energySigmaValue = Electron_energySigmaValue[ie]; 
    velectron.energySigmaUp = Electron_energySigmaUp[ie];
    velectron.energySigmaDown = Electron_energySigmaDown[ie];
    
    velectron.supcl_eta = Electron_supcl_eta[ie];
    velectron.supcl_phi = Electron_supcl_phi[ie];
    velectron.supcl_e = Electron_supcl_e[ie];
    velectron.supcl_rawE = Electron_supcl_rawE[ie];
    velectron.sigmaieta = Electron_sigmaieta[ie];
    velectron.sigmaiphi = Electron_sigmaiphi[ie];
    velectron.r9full = Electron_r9full[ie];
    
    velectron.hcaloverecal = Electron_hcaloverecal[ie];  
    velectron.ecloverpout = Electron_ecloverpout[ie];
    velectron.eoverp = Electron_eoverp[ie];
    velectron.hovere = Electron_hovere[ie];
  
    velectron.pfiso_drcor = Electron_pfiso_drcor[ie];
    velectron.pfiso_eacor = Electron_pfiso_eacor[ie];
    velectron.pfiso04_eacor = Electron_pfiso04_eacor[ie];
    velectron.pfisolsumphet = Electron_pfisolsumphet[ie];
    velectron.pfisolsumchhadpt = Electron_pfisolsumchhadpt[ie];
    velectron.pfsiolsumneuhadet = Electron_pfsiolsumneuhadet[ie];
    velectron.minchiso  = Electron_minchiso[ie];
    velectron.minnhiso  = Electron_minnhiso[ie];
    velectron.minphiso  = Electron_minphiso[ie];
    velectron.minisoall  = Electron_minisoall[ie];
    
    //velectron.fbrem = Electron_fbrem[ie];
    //velectron.normchi2 = Electron_normchi2[ie];
    //velectron.hitsmiss = Electron_hitsmiss[ie];
    //velectron.trkmeasure = Electron_trkmeasure[ie];
    //velectron.ecaletrkmomentum = Electron_ecaletrkmomentum[ie];
    //velectron.deltaetacltrkcalo = Electron_deltaetacltrkcalo[ie];
    //velectron.supcl_preshvsrawe = Electron_supcl_preshvsrawe[ie];
    //velectron.supcl_etaw = Electron_supcl_etaw[ie];
    //velectron.supcl_phiw = Electron_supcl_phiw[ie];
    //velectron.cloctftrkn = Electron_cloctftrkn[ie];
    //velectron.cloctftrkchi2 = Electron_cloctftrkchi2[ie];
    //velectron.e1x5bye5x5 = Electron_e1x5bye5x5[ie];
    // velectron.etain = Electron_etain[ie];
    //velectron.phiin = Electron_phiin[ie];
    
    velectron.p4.SetPtEtaPhiM(velectron.pt,velectron.eta,velectron.phi,velectron.mass);
     
    velectrons.push_back(velectron);
    
    if(int(velectrons.size()) >= maxsize) break;
    
  }
  
  sorted_by_pt(velectrons);	
	
}

void getLeptons(std::vector<Lepton> &vleptons, std::vector<Muon> vmuons, std::vector<Electron> velectrons, float pt_cut=30)
{ 
  for(unsigned imu=0; imu<vmuons.size(); imu++){
	if(vmuons[imu].pt < pt_cut) continue;
    Lepton vlepton;
    vlepton.pt = vmuons[imu].pt;
    vlepton.eta = vmuons[imu].eta;
    vlepton.phi = vmuons[imu].phi;
    vlepton.mass = vmuons[imu].mass;
    vlepton.charge = vmuons[imu].charge;
    vlepton.lepton_id = 1;
    vlepton.pdgId = 13;
    vlepton.p4 = vmuons[imu].p4;
	vlepton.indexemu = imu; 
    vleptons.push_back(vlepton);
  }
  for(unsigned ie=0; ie<velectrons.size(); ie++){
	if(velectrons[ie].pt < pt_cut) continue;
    Lepton vlepton;
    vlepton.pt = velectrons[ie].pt;
    vlepton.eta = velectrons[ie].eta;
    vlepton.phi = velectrons[ie].phi;
    vlepton.mass = velectrons[ie].mass;
    vlepton.charge = velectrons[ie].charge;
    vlepton.lepton_id = 2;
    vlepton.pdgId = 11;
    vlepton.p4 = velectrons[ie].p4;
    vlepton.indexemu=ie;
    vleptons.push_back(vlepton);
  }
  sorted_by_pt(vleptons);
}

void getAK4Genjets(std::vector<AK4GenJet> &genJets, float ptcut = 8.0, float etacut=2.5, bool isMC=false, int maxsize=njetmx)
{
	for(int ijet=0; ijet<(nGenJetAK4); ijet++){
		AK4GenJet gJet;
		if (fabs(GenJetAK4_eta[ijet]) >etacut) continue;
                if (GenJetAK4_pt[ijet] < ptcut) continue;

		gJet.pt = GenJetAK4_pt[ijet];
		gJet.eta = GenJetAK4_eta[ijet];
                gJet.phi = GenJetAK4_phi[ijet];
                gJet.mass = GenJetAK4_mass[ijet];
		gJet.hadronFlavor = GenJetAK4_hadronflav[ijet];
		gJet.partonFlavor = GenJetAK4_partonflav[ijet];
                gJet.p4.SetPtEtaPhiM(gJet.pt,gJet.eta,gJet.phi,gJet.mass);

   genJets.push_back(gJet);
   if(int(genJets.size())>=maxsize) break;
  }
  sorted_by_pt(genJets);
}

void getAK8Genjets(std::vector<AK8GenJet> &Jets, float ptcut=50, float etacut=2.5, bool isMC=false, int maxsize=njetmx)
{
    for(int ijet=0; ijet<(nGenJetAK8); ijet++){
    AK8GenJet lJet; 
	  
    if(fabs(GenJetAK8_eta[ijet])>etacut) continue;
    if(GenJetAK8_pt[ijet]<ptcut) continue;
    
    lJet.pt = GenJetAK8_pt[ijet];
    lJet.mass = GenJetAK8_mass[ijet];
    lJet.eta = GenJetAK8_eta[ijet];
    lJet.phi = GenJetAK8_phi[ijet];
    lJet.p4.SetPtEtaPhiM(lJet.pt,lJet.eta,lJet.phi,lJet.mass);
    
    Jets.push_back(lJet);
    
    if(int(Jets.size())>=maxsize) break;
    
  }
}

void getAK4jets(std::vector<AK4Jet> &Jets, float ptcut=30, float etacut=2.5, bool isMC=false, int maxsize=njetmx)
{
  
  for(int ijet=0; ijet<(nPFJetAK4); ijet++){
	 
	AK4Jet sJet; 
	  
    if(!PFJetAK4_jetID[ijet]) continue;
    
    PFJetAK4_pt[ijet] *= PFJetAK4_JEC[ijet] ;
    PFJetAK4_mass[ijet] *= PFJetAK4_JEC[ijet];
    
    if(isMC){
      PFJetAK4_pt[ijet] *= (1+PFJetAK4_reso[ijet]) ;
      PFJetAK4_mass[ijet] *= (1+PFJetAK4_reso[ijet]) ;
    }
	  
    if(fabs(PFJetAK4_eta[ijet])>etacut) continue;
    if(PFJetAK4_pt[ijet]<ptcut) continue;
    
    sJet.pt = PFJetAK4_pt[ijet];
    sJet.mass = PFJetAK4_mass[ijet];
    sJet.eta = PFJetAK4_eta[ijet];
    sJet.y = PFJetAK4_y[ijet];
    sJet.phi = PFJetAK4_phi[ijet];
    sJet.p4.SetPtEtaPhiM(sJet.pt,sJet.eta,sJet.phi,sJet.mass);
  
    sJet.jetID = PFJetAK4_jetID[ijet];
    sJet.jetID_tightlepveto = PFJetAK4_jetID_tightlepveto[ijet];
    sJet.hadronFlavour = PFJetAK4_hadronflav[ijet];
    sJet.partonFlavour = PFJetAK4_partonflav[ijet];
    sJet.btag_DeepFlav = PFJetAK4_btag_DeepFlav[ijet];
    sJet.btag_DeepCSV = PFJetAK4_btag_DeepCSV[ijet];
    sJet.PUID = PFJetAK4_PUID[ijet];
    sJet.qgl = PFJetAK4_qgl[ijet];
    
    sJet.reso = PFJetAK4_reso[ijet];
    sJet.resoup = PFJetAK4_resoup[ijet];
    sJet.resodn = PFJetAK4_resodn[ijet];
    
    sJet.JEC = PFJetAK4_JEC[ijet];
    sJet.jesup_AbsoluteStat = PFJetAK4_jesup_AbsoluteStat[ijet];
    sJet.jesup_AbsoluteScale = PFJetAK4_jesup_AbsoluteScale[ijet];
    sJet.jesup_AbsoluteMPFBias = PFJetAK4_jesup_AbsoluteMPFBias[ijet];
    sJet.jesup_FlavorQCD = PFJetAK4_jesup_FlavorQCD[ijet];
    sJet.jesup_Fragmentation = PFJetAK4_jesup_Fragmentation[ijet];
    sJet.jesup_PileUpDataMC = PFJetAK4_jesup_PileUpDataMC[ijet];
    sJet.jesup_PileUpPtBB = PFJetAK4_jesup_PileUpPtBB[ijet];
    sJet.jesup_PileUpPtEC1 = PFJetAK4_jesup_PileUpPtEC1[ijet];
    sJet.jesup_PileUpPtEC2 = PFJetAK4_jesup_PileUpPtEC2[ijet];
    sJet.jesup_PileUpPtRef = PFJetAK4_jesup_PileUpPtRef[ijet];
    sJet.jesup_RelativeFSR = PFJetAK4_jesup_RelativeFSR[ijet];
    sJet.jesup_RelativeJEREC1 = PFJetAK4_jesup_RelativeJEREC1[ijet];
    sJet.jesup_RelativeJEREC2 = PFJetAK4_jesup_RelativeJEREC2[ijet];
    sJet.jesup_RelativePtBB = PFJetAK4_jesup_RelativePtBB[ijet];
    sJet.jesup_RelativePtEC1 = PFJetAK4_jesup_RelativePtEC1[ijet];
    sJet.jesup_RelativePtEC2 = PFJetAK4_jesup_RelativePtEC2[ijet];
    sJet.jesup_RelativeBal = PFJetAK4_jesup_RelativeBal[ijet];
    sJet.jesup_RelativeSample = PFJetAK4_jesup_RelativeSample[ijet];
    sJet.jesup_RelativeStatEC = PFJetAK4_jesup_RelativeStatEC[ijet];
    sJet.jesup_RelativeStatFSR = PFJetAK4_jesup_RelativeStatFSR[ijet];
    sJet.jesup_SinglePionECAL = PFJetAK4_jesup_SinglePionECAL[ijet];
    sJet.jesup_SinglePionHCAL = PFJetAK4_jesup_SinglePionHCAL[ijet];
    sJet.jesup_TimePtEta = PFJetAK4_jesup_TimePtEta[ijet];
    sJet.jesup_Total = PFJetAK4_jesup_Total[ijet];
    sJet.jesdn_AbsoluteStat = PFJetAK4_jesdn_AbsoluteStat[ijet];
    sJet.jesdn_AbsoluteScale = PFJetAK4_jesdn_AbsoluteScale[ijet];
    sJet.jesdn_AbsoluteMPFBias = PFJetAK4_jesdn_AbsoluteMPFBias[ijet];
    sJet.jesdn_FlavorQCD = PFJetAK4_jesdn_FlavorQCD[ijet];
    sJet.jesdn_Fragmentation = PFJetAK4_jesdn_Fragmentation[ijet];
    sJet.jesdn_PileUpDataMC = PFJetAK4_jesdn_PileUpDataMC[ijet];
    sJet.jesdn_PileUpPtBB = PFJetAK4_jesdn_PileUpPtBB[ijet];
    sJet.jesdn_PileUpPtEC1 = PFJetAK4_jesdn_PileUpPtEC1[ijet];
    sJet.jesdn_PileUpPtEC2 = PFJetAK4_jesdn_PileUpPtEC2[ijet];
    sJet.jesdn_PileUpPtRef = PFJetAK4_jesdn_PileUpPtRef[ijet];
    sJet.jesdn_RelativeFSR = PFJetAK4_jesdn_RelativeFSR[ijet];
    sJet.jesdn_RelativeJEREC1 = PFJetAK4_jesdn_RelativeJEREC1[ijet];
    sJet.jesdn_RelativeJEREC2 = PFJetAK4_jesdn_RelativeJEREC2[ijet];
    sJet.jesdn_RelativePtBB = PFJetAK4_jesdn_RelativePtBB[ijet];
    sJet.jesdn_RelativePtEC1 = PFJetAK4_jesdn_RelativePtEC1[ijet];
    sJet.jesdn_RelativePtEC2 = PFJetAK4_jesdn_RelativePtEC2[ijet];
    sJet.jesdn_RelativeBal = PFJetAK4_jesdn_RelativeBal[ijet];
    sJet.jesdn_RelativeSample = PFJetAK4_jesdn_RelativeSample[ijet];
    sJet.jesdn_RelativeStatEC = PFJetAK4_jesdn_RelativeStatEC[ijet];
    sJet.jesdn_RelativeStatFSR = PFJetAK4_jesdn_RelativeStatFSR[ijet];
    sJet.jesdn_SinglePionECAL = PFJetAK4_jesdn_SinglePionECAL[ijet];
    sJet.jesdn_SinglePionHCAL = PFJetAK4_jesdn_SinglePionHCAL[ijet];
    sJet.jesdn_TimePtEta = PFJetAK4_jesdn_TimePtEta[ijet];
    sJet.jesdn_Total = PFJetAK4_jesdn_Total[ijet];
    
    Jets.push_back(sJet);
    
    if(int(Jets.size())>=maxsize) break;
    
  }

  sorted_by_pt(Jets);
			
}

void getAK8jets(std::vector<AK8Jet> &LJets, float ptcut=200, float etacut=2.5, bool isMC=false, int maxsize=njetmxAK8)
{

  for(int ijet=0; ijet<(nPFJetAK8); ijet++){
  
    if(!PFJetAK8_jetID[ijet]) continue;
    
    PFJetAK8_pt[ijet] *= PFJetAK8_JEC[ijet] ;
    PFJetAK8_mass[ijet] *= PFJetAK8_JEC[ijet];
    
    if(isMC){
      PFJetAK8_pt[ijet] *= (1+PFJetAK8_reso[ijet]) ;
      PFJetAK8_mass[ijet] *= (1+PFJetAK8_reso[ijet]) ;
    }
				
    if(fabs(PFJetAK8_eta[ijet])>etacut) continue;
    if(PFJetAK8_pt[ijet] < ptcut) continue;
    
    AK8Jet LJet;
    
    LJet.jetID = PFJetAK8_jetID[ijet];
    LJet.jetID_tightlepveto = PFJetAK8_jetID_tightlepveto[ijet];
    
    LJet.pt = PFJetAK8_pt[ijet];
    LJet.eta = PFJetAK8_eta[ijet];
    LJet.mass = PFJetAK8_mass[ijet];
    LJet.phi = PFJetAK8_phi[ijet];
    LJet.y = PFJetAK8_y[ijet];
    LJet.p4.SetPtEtaPhiM(LJet.pt,LJet.eta,LJet.phi,LJet.mass);
    
    LJet.sdmass = PFJetAK8_sdmass[ijet];
    LJet.tau21 = PFJetAK8_tau2[ijet]*1./max(float(1.e-6),PFJetAK8_tau1[ijet]);
    LJet.tau32 = PFJetAK8_tau3[ijet]*1./max(float(1.e-6),PFJetAK8_tau2[ijet]);
   
    LJet.btag_DeepCSV = PFJetAK8_btag_DeepCSV[ijet];
    LJet.DeepTag_DAK8MD_TvsQCD = PFJetAK8_DeepTag_DAK8MD_TvsQCD[ijet];
    LJet.DeepTag_DAK8MD_WvsQCD = PFJetAK8_DeepTag_DAK8MD_WvsQCD[ijet];
    LJet.DeepTag_DAK8MD_ZvsQCD = PFJetAK8_DeepTag_DAK8MD_ZvsQCD[ijet];
    LJet.DeepTag_DAK8MD_HvsQCD = PFJetAK8_DeepTag_DAK8MD_HvsQCD[ijet];
    LJet.DeepTag_DAK8MD_bbvsQCD = PFJetAK8_DeepTag_DAK8MD_bbvsQCD[ijet];
    LJet.DeepTag_PNet_TvsQCD = PFJetAK8_DeepTag_PNet_TvsQCD[ijet];
    LJet.DeepTag_PNet_WvsQCD = PFJetAK8_DeepTag_PNet_WvsQCD[ijet];
    LJet.DeepTag_PNet_ZvsQCD = PFJetAK8_DeepTag_PNet_ZvsQCD[ijet];
    LJet.DeepTag_PNetMD_XbbvsQCD = PFJetAK8_DeepTag_PNetMD_XbbvsQCD[ijet];
    LJet.DeepTag_PNetMD_XccvsQCD = PFJetAK8_DeepTag_PNetMD_XccvsQCD[ijet];
    LJet.DeepTag_PNetMD_XqqvsQCD = PFJetAK8_DeepTag_PNetMD_XqqvsQCD[ijet];
    LJet.DeepTag_PNetMD_QCD = PFJetAK8_DeepTag_PNetMD_QCD[ijet];
    
    LJet.CHF = PFJetAK8_CHF[ijet];
    LJet.NHF = PFJetAK8_NHF[ijet];
    LJet.CEMF = PFJetAK8_CEMF[ijet];
    LJet.NEMF = PFJetAK8_NEMF[ijet];
    LJet.MUF = PFJetAK8_MUF[ijet];
    LJet.PHF = PFJetAK8_PHF[ijet];
    LJet.EEF = PFJetAK8_EEF[ijet];
    LJet.HFHF = PFJetAK8_HFHF[ijet];
    LJet.CHM = PFJetAK8_CHM[ijet];
    LJet.NHM = PFJetAK8_NHM[ijet];
    LJet.MUM = PFJetAK8_MUM[ijet];
    LJet.EEM = PFJetAK8_EEM[ijet];
    LJet.PHM = PFJetAK8_PHM[ijet];
    LJet.HFHM = PFJetAK8_HFHM[ijet];
    LJet.Neucons = PFJetAK8_Neucons[ijet];
    LJet.Chcons = PFJetAK8_Chcons[ijet];
    
    LJet.sub1pt = PFJetAK8_sub1pt[ijet];
    LJet.sub1eta = PFJetAK8_sub1eta[ijet];
    LJet.sub1phi = PFJetAK8_sub1phi[ijet];
    LJet.sub1mass = PFJetAK8_sub1mass[ijet];
    LJet.sub1btag = PFJetAK8_sub1btag[ijet];
   
    LJet.sub2pt = PFJetAK8_sub2pt[ijet];
    LJet.sub2eta = PFJetAK8_sub2eta[ijet];
    LJet.sub2phi = PFJetAK8_sub2phi[ijet];
    LJet.sub2mass = PFJetAK8_sub2mass[ijet];
    LJet.sub2btag = PFJetAK8_sub2btag[ijet];
        
    LJet.reso = PFJetAK8_reso[ijet];
    LJet.resoup = PFJetAK8_resoup[ijet];
    LJet.resodn = PFJetAK8_resodn[ijet];
    
    LJet.JEC = PFJetAK8_JEC[ijet];
    LJet.jesup_AbsoluteStat = PFJetAK8_jesup_AbsoluteStat[ijet];
    LJet.jesup_AbsoluteScale = PFJetAK8_jesup_AbsoluteScale[ijet];
    LJet.jesup_AbsoluteMPFBias = PFJetAK8_jesup_AbsoluteMPFBias[ijet];
    LJet.jesup_FlavorQCD = PFJetAK8_jesup_FlavorQCD[ijet];
    LJet.jesup_Fragmentation = PFJetAK8_jesup_Fragmentation[ijet];
    LJet.jesup_PileUpDataMC = PFJetAK8_jesup_PileUpDataMC[ijet];
    LJet.jesup_PileUpPtBB = PFJetAK8_jesup_PileUpPtBB[ijet];
    LJet.jesup_PileUpPtEC1 = PFJetAK8_jesup_PileUpPtEC1[ijet];
    LJet.jesup_PileUpPtEC2 = PFJetAK8_jesup_PileUpPtEC2[ijet];
    LJet.jesup_PileUpPtRef = PFJetAK8_jesup_PileUpPtRef[ijet];
    LJet.jesup_RelativeFSR = PFJetAK8_jesup_RelativeFSR[ijet];
    LJet.jesup_RelativeJEREC1 = PFJetAK8_jesup_RelativeJEREC1[ijet];
    LJet.jesup_RelativeJEREC2 = PFJetAK8_jesup_RelativeJEREC2[ijet];
    LJet.jesup_RelativePtBB = PFJetAK8_jesup_RelativePtBB[ijet];
    LJet.jesup_RelativePtEC1 = PFJetAK8_jesup_RelativePtEC1[ijet];
    LJet.jesup_RelativePtEC2 = PFJetAK8_jesup_RelativePtEC2[ijet];
    LJet.jesup_RelativeBal = PFJetAK8_jesup_RelativeBal[ijet];
    LJet.jesup_RelativeSample = PFJetAK8_jesup_RelativeSample[ijet];
    LJet.jesup_RelativeStatEC = PFJetAK8_jesup_RelativeStatEC[ijet];
    LJet.jesup_RelativeStatFSR = PFJetAK8_jesup_RelativeStatFSR[ijet];
    LJet.jesup_SinglePionECAL = PFJetAK8_jesup_SinglePionECAL[ijet];
    LJet.jesup_SinglePionHCAL = PFJetAK8_jesup_SinglePionHCAL[ijet];
    LJet.jesup_TimePtEta = PFJetAK8_jesup_TimePtEta[ijet];
    LJet.jesup_Total = PFJetAK8_jesup_Total[ijet];
    LJet.jesdn_AbsoluteStat = PFJetAK8_jesdn_AbsoluteStat[ijet];
    LJet.jesdn_AbsoluteScale = PFJetAK8_jesdn_AbsoluteScale[ijet];
    LJet.jesdn_AbsoluteMPFBias = PFJetAK8_jesdn_AbsoluteMPFBias[ijet];
    LJet.jesdn_FlavorQCD = PFJetAK8_jesdn_FlavorQCD[ijet];
    LJet.jesdn_Fragmentation = PFJetAK8_jesdn_Fragmentation[ijet];
    LJet.jesdn_PileUpDataMC = PFJetAK8_jesdn_PileUpDataMC[ijet];
    LJet.jesdn_PileUpPtBB = PFJetAK8_jesdn_PileUpPtBB[ijet];
    LJet.jesdn_PileUpPtEC1 = PFJetAK8_jesdn_PileUpPtEC1[ijet];
    LJet.jesdn_PileUpPtEC2 = PFJetAK8_jesdn_PileUpPtEC2[ijet];
    LJet.jesdn_PileUpPtRef = PFJetAK8_jesdn_PileUpPtRef[ijet];
    LJet.jesdn_RelativeFSR = PFJetAK8_jesdn_RelativeFSR[ijet];
    LJet.jesdn_RelativeJEREC1 = PFJetAK8_jesdn_RelativeJEREC1[ijet];
    LJet.jesdn_RelativeJEREC2 = PFJetAK8_jesdn_RelativeJEREC2[ijet];
    LJet.jesdn_RelativePtBB = PFJetAK8_jesdn_RelativePtBB[ijet];
    LJet.jesdn_RelativePtEC1 = PFJetAK8_jesdn_RelativePtEC1[ijet];
    LJet.jesdn_RelativePtEC2 = PFJetAK8_jesdn_RelativePtEC2[ijet];
    LJet.jesdn_RelativeBal = PFJetAK8_jesdn_RelativeBal[ijet];
    LJet.jesdn_RelativeSample = PFJetAK8_jesdn_RelativeSample[ijet];
    LJet.jesdn_RelativeStatEC = PFJetAK8_jesdn_RelativeStatEC[ijet];
    LJet.jesdn_RelativeStatFSR = PFJetAK8_jesdn_RelativeStatFSR[ijet];
    LJet.jesdn_SinglePionECAL = PFJetAK8_jesdn_SinglePionECAL[ijet];
    LJet.jesdn_SinglePionHCAL = PFJetAK8_jesdn_SinglePionHCAL[ijet];
    LJet.jesdn_TimePtEta = PFJetAK8_jesdn_TimePtEta[ijet];
    LJet.jesdn_Total = PFJetAK8_jesdn_Total[ijet];
    
    LJet.haselectron = LJet.hasmuon = LJet.hastau = LJet.hasqg = LJet.hasb = LJet.hasleptop = LJet.hashadtop = LJet.hastop = LJet.hasmatchmu = LJet.hasmatche = false;
    LJet.hasleptop_alldecay = LJet.hashadtop_alldecay = false;
    LJet.matchAK4deepb = -100;
    
    LJets.push_back(LJet);
    
    if(int(LJets.size())>=maxsize) break;
  }

  sorted_by_pt(LJets);
  	
}

void LeptonJet_cleaning(std::vector<AK4Jet> &Jets, std::vector<Lepton> Leptons, float dR_cut=0.4, float ptcut=30, float etacut=2.5)
{
  
   auto jet = Jets.begin();
   while (jet != Jets.end()){
	
	for(auto & lepton: Leptons){
		
		if(delta2R(jet->eta, jet->phi, lepton.eta, lepton.phi)<dR_cut){
			jet->p4 -= lepton.p4;
			jet->pt = (jet->p4).Pt();
			jet->eta = (jet->p4).Eta();
			jet->y = (jet->p4).Rapidity();
			jet->phi = (jet->p4).Phi();
			jet->mass = (jet->p4).M();
		}
	}
	
	if(jet->pt<ptcut || fabs(jet->eta)>etacut){
		Jets.erase(jet);
		}
	else{
		++jet;
		}
  }	

  sorted_by_pt(Jets);

}

void LeptonJet_cleaning(std::vector<AK8Jet> &Jets, std::vector<Lepton> Leptons, float dR_cut=0.4, float ptcut=30, float etacut=2.5)
{
  
   auto jet = Jets.begin();
   while (jet != Jets.end()){
	
	for(auto & lepton: Leptons){
		
		if(delta2R(jet->eta, jet->phi, lepton.eta, lepton.phi)<dR_cut){
			jet->p4 -= lepton.p4;
			jet->pt = (jet->p4).Pt();
			jet->eta = (jet->p4).Eta();
			jet->y = (jet->p4).Rapidity();
			jet->phi = (jet->p4).Phi();
			jet->mass = (jet->p4).M();
		}
	}
	
	if(jet->pt<ptcut || fabs(jet->eta)>etacut){
		Jets.erase(jet);
		}
	else{
		++jet;
		}
  }	

  sorted_by_pt(Jets);

}

void getLHEParts(std::vector<GenParton> &LHEParts, int maxsize=njetmx)
{
  
  for(int igen=0; igen<(nLHEPart); igen++){
	 
	GenParton parton; 
	
	parton.pt = LHEPart_pt[igen];
	parton.eta = LHEPart_eta[igen];
	parton.phi = LHEPart_phi[igen];
	parton.mass = LHEPart_m[igen];  
	parton.p4.SetPtEtaPhiM(parton.pt,parton.eta,parton.phi,parton.mass);
    
	parton.pdgId = LHEPart_pdg[igen];
		
    LHEParts.push_back(parton);
    
    if(int(LHEParts.size())>=maxsize) break;
    
  }

//	GenPartons.sorted_by_pt();	
	
}

void getPartons(std::vector<GenParton> &GenPartons, int maxsize=npartmx)
{
  
  for(int igen=0; igen<(nGenPart); igen++){
	 
	GenParton parton; 
	
	parton.pt = GenPart_pt[igen];
	parton.eta = GenPart_eta[igen];
	parton.phi = GenPart_phi[igen];
	parton.mass = GenPart_m[igen];  
	parton.p4.SetPtEtaPhiM(parton.pt,parton.eta,parton.phi,parton.mass);
    
    parton.status = GenPart_status[igen];
	parton.pdgId = GenPart_pdgId[igen];
	parton.mompdgId = GenPart_pdgId[igen];
	parton.grmompdgId = GenPart_pdgId[igen]; 
	
	parton.fromhard = GenPart_fromhard[igen];
	parton.fromhardbFSR = GenPart_fromhardbFSR[igen];
	parton.isPromptFinalState = GenPart_isPromptFinalState[igen];
	parton.isLastCopyBeforeFSR = GenPart_isLastCopyBeforeFSR[igen]; 
		
    GenPartons.push_back(parton);
    
    if(int(GenPartons.size())>=maxsize) break;
    
  }

//	GenPartons.sorted_by_pt();	
	
}

void getTrigObjs(std::vector<TrigObj> &trigobjects, int maxsize=njetmx)
{
  
  for(int iobj=0; iobj<(nTrigObj); iobj++){
	 
	TrigObj trigobj; 
	
	trigobj.pt = TrigObj_pt[iobj];
	trigobj.eta = TrigObj_eta[iobj];
	trigobj.phi = TrigObj_phi[iobj];
	trigobj.mass = TrigObj_mass[iobj];  
	trigobj.p4.SetPtEtaPhiM(trigobj.pt,trigobj.eta,trigobj.phi,trigobj.mass);
    
	trigobj.pdgId = TrigObj_pdgId[iobj];
	trigobj.type = TrigObj_type[iobj];
	trigobj.L1 = TrigObj_L1[iobj];
	trigobj.HLT = TrigObj_HLT[iobj];
		
    trigobjects.push_back(trigobj);
    
    if(int(trigobjects.size())>=maxsize) break;
    
  }
	
}

void getLHETops(std::vector<GenParton> &LHETops, std::vector<GenParton> GenPartons)
{
  
  for(auto & part: GenPartons){
	  
	  if(abs(part.status)!=22) continue;
      if(!(part.fromhard)) continue;
      if(abs(part.pdgId)!=6) continue;
	  
	  LHETops.push_back(part);
	  
	  }
}

void getGENTops(vector<TopQuark> &gentops, vector<GenParton> genpartons)  // with daughters after shower
{     
    vector<GenParton> W_dau;
    vector<GenParton> t_bp;
  
    for(unsigned igen=0; igen<genpartons.size(); igen++){
      
      if(!(genpartons[igen].status==23 || genpartons[igen].status==1)) continue;
      if(!(genpartons[igen].fromhard)) continue;
      
      if(!((abs(genpartons[igen].pdgId)>=1 && abs(genpartons[igen].pdgId)<=5)||(abs(genpartons[igen].pdgId)>=11 && abs(genpartons[igen].pdgId)<=16))) continue;
      if(!(abs(genpartons[igen].mompdgId)==6 || abs(genpartons[igen].mompdgId)==24)) continue;
      if(abs(genpartons[igen].mompdgId)==24 && abs(genpartons[igen].grmompdgId)!=6) continue;
            
      if(abs(genpartons[igen].pdgId)>=1 && abs(genpartons[igen].pdgId)<5 && abs(genpartons[igen].mompdgId)==24 && abs(genpartons[igen].grmompdgId)==6)   {  W_dau.push_back(genpartons[igen]); }
      if(abs(genpartons[igen].pdgId)>=11 && abs(genpartons[igen].pdgId)<=16 && abs(genpartons[igen].mompdgId)==24 && abs(genpartons[igen].grmompdgId)==6) {  W_dau.push_back(genpartons[igen]); }
      if(abs(genpartons[igen].pdgId)==5 && abs(genpartons[igen].mompdgId)==6) {  t_bp.push_back(genpartons[igen]); }
    }
    
	
	for(unsigned ipart=0; ipart<W_dau.size(); ipart++){
		  
	  if(!((abs(W_dau[ipart].pdgId)==11 || abs(W_dau[ipart].pdgId)==13 || abs(W_dau[ipart].pdgId)==15 || abs(W_dau[ipart].pdgId)==1 || abs(W_dau[ipart].pdgId)==3))) continue;
		  
		unsigned partner = -1;
		unsigned match_b = -1;
		  
		for(unsigned jpart=(ipart+1); jpart<W_dau.size(); jpart++){
			if((W_dau[ipart].mompdgId==W_dau[jpart].mompdgId) && (W_dau[ipart].grmompdgId==W_dau[jpart].grmompdgId) && (W_dau[ipart].pdgId*W_dau[jpart].pdgId<0)){
				 partner = jpart;
				 break;
			 }
		 }
		  
		 for(unsigned ib=0; ib<t_bp.size(); ib++){
		   if(t_bp[ib].mompdgId==W_dau[ipart].grmompdgId){
		    	 match_b = ib;
				 break;
			}
		 }
		  
		 GenParton q1, q2, b;
		  
		  q1 = W_dau[ipart];
		 
		  if(int(partner)>=0 && partner<W_dau.size()){
			  q2 = W_dau[partner];
		  }
		  if(int(match_b)>=0 && match_b<t_bp.size()){
			  b = t_bp[match_b];
		  }
		  
		  if(int(partner)>=0 && partner<W_dau.size() && int(match_b)>=0 && match_b<t_bp.size()){
			
			TopQuark topQ;
			  
			topQ.p4 = (b.p4 + q1.p4 + q2.p4);
			topQ.daughter.push_back(q1);
			topQ.daughter.push_back(q2);
			topQ.daughter.push_back(b);
			
			gentops.push_back(topQ);	
		  }
    }
  
}



void TopAssignment_toJet(std::vector<AK8Jet> &LJets, std::vector<GenParton> lhetops, std::vector<TopQuark> gentops)
{

  for(unsigned ijet=0; ijet<LJets.size(); ijet++){
     
	for(unsigned itop=0; itop<lhetops.size(); itop++)
	{
	  if(delta2R(LJets[ijet].y,LJets[ijet].phi,lhetops[itop].eta,lhetops[itop].phi)<0.6){
	    LJets[ijet].hastop = true;
	    break;
	  }
	}
    
    for(unsigned itop=0; itop<gentops.size(); itop++){
		
		if(abs(gentops[itop].daughter[0].pdgId)==11 || abs(gentops[itop].daughter[0].pdgId)==13 || abs(gentops[itop].daughter[0].pdgId)==15){
			
			if(delta2R(LJets[ijet].y,LJets[ijet].phi,gentops[itop].p4.Rapidity(),gentops[itop].p4.Phi())<0.8){
				LJets[ijet].hasleptop = true;	
			}
			
			if((delta2R(LJets[ijet].y,LJets[ijet].phi,gentops[itop].daughter[0].eta,gentops[itop].daughter[0].phi)<0.8)
			&& (delta2R(LJets[ijet].y,LJets[ijet].phi,gentops[itop].daughter[1].eta,gentops[itop].daughter[1].phi)<0.8)
			&& (delta2R(LJets[ijet].y,LJets[ijet].phi,gentops[itop].daughter[2].eta,gentops[itop].daughter[2].phi)<0.8))
			{
				LJets[ijet].hasleptop_alldecay = true;
			}
		
		}
		
		else if(abs(gentops[itop].daughter[0].pdgId)==1 || abs(gentops[itop].daughter[0].pdgId)==3){
			
			if(delta2R(LJets[ijet].y,LJets[ijet].phi,gentops[itop].p4.Rapidity(),gentops[itop].p4.Phi())<0.8){
				LJets[ijet].hashadtop = true;	
			}
			
			if((delta2R(LJets[ijet].y,LJets[ijet].phi,gentops[itop].daughter[0].eta,gentops[itop].daughter[0].phi)<0.8)
			&& (delta2R(LJets[ijet].y,LJets[ijet].phi,gentops[itop].daughter[1].eta,gentops[itop].daughter[1].phi)<0.8)
			&& (delta2R(LJets[ijet].y,LJets[ijet].phi,gentops[itop].daughter[2].eta,gentops[itop].daughter[2].phi)<0.8))
			{
				LJets[ijet].hashadtop_alldecay = true;
			}
		
		}
	}
      
  }//ijet
	
	
}

void AssignGen(std::vector<AK8Jet> &LJets, std::vector<GenParton> GenPartons){

	for(auto & ljet: LJets){
		
		for(auto & part: GenPartons){
			
			if(abs(part.status)!=23 && part.status!=1) continue;
			if(!(part.fromhard)) continue;
			
			if(abs(part.pdgId)==11 && delta2R(ljet.y,ljet.phi,part.eta,part.phi)<0.8){
				ljet.haselectron  = true;
				break;
			 }
		
			if(abs(part.pdgId)==13 && delta2R(ljet.y,ljet.phi,part.eta,part.phi)<0.8){
				ljet.hasmuon  = true;
				break;
			 }
			
			if(abs(part.pdgId)==15 && delta2R(ljet.y,ljet.phi,part.eta,part.phi)<0.8){
				ljet.hastau  = true;
				break;
			 }
		
			if(((abs(part.pdgId)>=1 && abs(part.pdgId)<5) || abs(part.pdgId)==21) && delta2R(ljet.y,ljet.phi,part.eta,part.phi)<0.8){
				ljet.hasqg  = true;
				break;
			 }
			 
			if(abs(part.pdgId)==5 && delta2R(ljet.y,ljet.phi,part.eta,part.phi)<0.8){
				ljet.hasb  = true;
				break;
			 } 
		}
	}
}

bool isBJet(AK4Jet jet, float btag_cut)
{
if(jet.btag_DeepFlav >= btag_cut){
	return true;
}
else{
	return false;
	}	
}

/*
void Match_trigger(vector<bool> double_hlts, vector<bool> single_hlts, 
									  vector<vector<float>> double_pt_cuts, vector<float> single_pt_cuts, 
									  vector<vector<int>> double_pids, vector<int> single_pids, 
									  vector<float> single_other_pt_cuts, vector<int> single_other_pids,
									  vector<std::pair<int,TLorentzVector> > TrigRefObj,
									  Lepton lepcand_1, Lepton lepcand_2, vector<AK4Jet> Jets,
									  bool &trig_threshold_pass,
									  bool &trig_matching_pass,
									  vector<TH1D*> &hist_init
									  )
{
	
  if(double_hlts.size()<1 && single_hlts.size()<1) { 
	  return; 
	  }
  
  else{

	// checking if offline objects passed trigger thresholds //
	
	bool double_trig_pass(false), single_trig_pass(false);
	bool any_double_hlt_pass = false;
	
	if(double_hlts.size()>0){
		for(unsigned ihlt=0; ihlt<double_hlts.size(); ihlt++){
			if (double_hlts[ihlt]){ any_double_hlt_pass = true; }
			if (double_hlts[ihlt] && ((lepcand_1.pt>double_pt_cuts[ihlt][0] && lepcand_1.pdgId==double_pids[ihlt][0] && lepcand_2.pt>double_pt_cuts[ihlt][1] && lepcand_2.pdgId==double_pids[ihlt][1])
			||  (lepcand_1.pt>double_pt_cuts[ihlt][1] && lepcand_1.pdgId==double_pids[ihlt][1] && lepcand_2.pt>double_pt_cuts[ihlt][0] && lepcand_2.pdgId==double_pids[ihlt][0]))
			) { 
				trig_threshold_pass = true; 
				double_trig_pass = true;
				break; }
			}
		}
	
	int fired_single_trig = -1;
	
	if(!any_double_hlt_pass){
		if(single_hlts.size()>0){
			for(unsigned ihlt=0; ihlt<single_hlts.size(); ihlt++){
				if (single_hlts[ihlt] && ((lepcand_1.pt>single_pt_cuts[ihlt] && lepcand_1.pdgId==single_pids[ihlt]) || (lepcand_2.pt>single_pt_cuts[ihlt] && lepcand_2.pdgId==single_pids[ihlt])) && single_other_pt_cuts[ihlt]>0 && Jets.size()>0 && Jets[0].pt>single_other_pt_cuts[ihlt])  { 
					trig_threshold_pass = true; 
					single_trig_pass = true;
					fired_single_trig = ihlt;
					break; 
					}
				else if (single_hlts[ihlt] && ((lepcand_1.pt>single_pt_cuts[ihlt] && lepcand_1.pdgId==single_pids[ihlt]) || (lepcand_2.pt>single_pt_cuts[ihlt] && lepcand_2.pdgId==single_pids[ihlt])) && single_other_pt_cuts[ihlt]<0)
					{
					trig_threshold_pass = true; 
					single_trig_pass = true;
					fired_single_trig = ihlt;
					break; 
					}
				}
			}
		}
  
	// check if offline objects match to trigger objects //
  
	bool lep1_match = false;
	bool lep2_match = false;
	bool jet_match = false;
  
    float ptratmin_lep1(0.3), ptratmin_lep2(0.3), ptratmin_jet(0.3);
  
	for (uint trv=0; trv<TrigRefObj.size(); trv++) {
		
		bool mutrobj(false), eltrobj(false), jettrobj(false);
		if (abs(TrigRefObj[trv].first)==13) mutrobj=true;
		else if (abs(TrigRefObj[trv].first)==0 && TrigRefObj[trv].second.M() < 1.e-3) eltrobj=true;
		else if (abs(TrigRefObj[trv].first)==0 && TrigRefObj[trv].second.M() > 1.e-3) jettrobj=true;

		TVector3 Trv = TrigRefObj[trv].second.Vect();
		
		if (mutrobj||eltrobj) {
			
			TVector3 flep1v = lepcand_1.p4.Vect();
			TVector3 flep2v = lepcand_2.p4.Vect();

			double tmprat1 = fabs((flep1v-Trv).Mag()/max(1.e-6,flep1v.Mag()));
			double tmprat2 = fabs((flep2v-Trv).Mag()/max(1.e-6,flep2v.Mag()));
			 
			//hist_init[0]->Fill(tmprat1, weight);
			//hist_init[1]->Fill(tmprat2, weight);
      
			if (tmprat1 < ptratmin_lep1) {
				if((mutrobj && lepcand_1.pdgId==13)||(eltrobj && lepcand_1.pdgId==11)){
					lep1_match = true;
					ptratmin_lep1 = tmprat1;
				}
			}
			
			if (tmprat2 < ptratmin_lep2) {
				if((mutrobj && lepcand_2.pdgId==13)||(eltrobj && lepcand_2.pdgId==11)){
					lep2_match = true;
					ptratmin_lep2 = tmprat2;
				}
			}
		}
		
		if(jettrobj){
			
			double tmprat = fabs((Jets[0].p4.Vect()-Trv).Mag()/max(1.e-6,(Jets[0].p4.Vect().Mag())));
			if (tmprat < ptratmin_jet) {
				jet_match = true;
				ptratmin_jet = tmprat;
			}
			
		}
	}
	
	//if (((double_hlts.size()>0 && double_trig_pass && lep1_match && lep2_match) || (single_hlts.size()>0 && single_trig_pass && lep1_match))) { trig_matching_pass = true; }
	if (double_hlts.size()>0 && double_trig_pass && lep1_match && lep2_match)  { trig_matching_pass = true; }
	else if (single_hlts.size()>0 && single_trig_pass && (lep1_match||lep2_match) && single_other_pt_cuts[fired_single_trig]>0 && jet_match) { trig_matching_pass = true; }
	else if (single_hlts.size()>0 && single_trig_pass && (lep1_match||lep2_match) && single_other_pt_cuts[fired_single_trig]<0) { trig_matching_pass = true; }
	
  }
   
}
*/
void Match_trigger(vector<bool> double_hlts, vector<vector<float>> double_pt_cuts, vector<vector<int>> double_pids,
				   vector<bool> single_hlts, vector<float> single_pt_cuts, vector<int> single_pids, vector<float> single_other_pt_cuts, vector<int> single_other_pids,
				   vector<bool> jet_hlts, vector<float> jet_pt_cuts, vector<int> jet_pids,
				   vector<TrigObj> trigobj,
				   vector<Muon> vmuons, vector<Electron> velectrons, vector<Lepton> vleptons, vector<AK4Jet> Jets, vector<AK8Jet> LJets,
				   bool &anytrig_pass,
				   bool &trig_threshold_pass,
  				   bool &trig_matching_pass,
  				   bool &muon_trig_pass, bool &electron_trig_pass)
{
	
	anytrig_pass = false;
	trig_threshold_pass = false;
	trig_matching_pass = false;
	
	/*
	if(channel_DL){
		if(double_hlts.size()>0){
			for(unsigned ihlt=0; ihlt<double_hlts.size(); ihlt++){
				if(double_hlts[ihlt]){ anytrig_pass = true; break; }
			}
		}
	}
    */
    
    vector<tuple<int, unsigned, bool> > pass_hlt;
    
    if(!anytrig_pass && single_hlts.size()>0){
        
        for(unsigned ihlt=0; ihlt<single_hlts.size(); ihlt++){
			
			if(single_hlts[ihlt]){ 
				
				anytrig_pass = true; 
				
				bool trig_cut = false;
				
				for(unsigned ilep=0; ilep<vleptons.size(); ilep++){
					if(vleptons[ilep].pt>single_pt_cuts[ihlt] && abs(vleptons[ilep].pdgId)==single_pids[ihlt]){
						if(single_other_pids[ihlt]==0){
							trig_threshold_pass = true;
							trig_cut = true;
							}
						else if(single_other_pids[ihlt]==1){
							for(unsigned ijet=0; ijet<Jets.size(); ijet++){
								if(Jets[ilep].pt>single_other_pt_cuts[ihlt]){
									trig_threshold_pass = true;
									trig_cut = true;
									break;
									}
								}
							}
						}
					}
				
				pass_hlt.push_back(make_tuple(0,ihlt,trig_cut));
				
				}
								
			}
		}

	if(!anytrig_pass && jet_hlts.size()>0){
		
        for(unsigned ihlt=0; ihlt<jet_hlts.size(); ihlt++){
			
			if(jet_hlts[ihlt]){ 
				
				anytrig_pass = true; 
				
				bool trig_cut = false;
				
				if(jet_pids[ihlt]==2){
					for(unsigned ijet=0; ijet<LJets.size(); ijet++){
						if(LJets[ijet].pt > jet_pt_cuts[ihlt]) {
							trig_threshold_pass = true;
							trig_cut = true;
							break;
						}
					}
					pass_hlt.push_back(make_tuple(2,ihlt,trig_cut));
				}
				
				if(jet_pids[ihlt]==1){
					for(unsigned ijet=0; ijet<Jets.size(); ijet++){
						if(Jets[ijet].pt > jet_pt_cuts[ihlt]) {
							trig_threshold_pass = true;
							trig_cut = true;
							break;
						}
					}
					pass_hlt.push_back(make_tuple(1,ihlt,trig_cut));
				}
			}
		}
	}
     
    // Matching with trigger object // 
     
	float ptvar_min(0.3);
	float dR_min(0.4);
	
	if(trig_threshold_pass && pass_hlt.size()>0){
		
		if(get<0>(pass_hlt[0])==0){
				
			if(!trig_matching_pass && single_pids[get<1>(pass_hlt[0])]==11 && velectrons.size()>0){
				for(unsigned iobj=0; iobj<trigobj.size(); iobj++){
				
					float ptvar = fabs(velectrons[0].pt - trigobj[iobj].pt)/velectrons[0].pt;
					if(ptvar<ptvar_min && delta2R(velectrons[0].eta,velectrons[0].phi,trigobj[iobj].eta,trigobj[iobj].phi)<dR_min && abs(trigobj[iobj].pdgId)!=13) { trig_matching_pass = true; }

				}				
			}
				
			if(!trig_matching_pass && single_pids[get<1>(pass_hlt[0])]==13 && vmuons.size()>0){
				for(unsigned iobj=0; iobj<trigobj.size(); iobj++){
				
					float ptvar = fabs(vmuons[0].pt - trigobj[iobj].pt)/vmuons[0].pt;
					if(ptvar<ptvar_min && delta2R(vmuons[0].eta,vmuons[0].phi,trigobj[iobj].eta,trigobj[iobj].phi)<dR_min && abs(trigobj[iobj].pdgId)==13) { trig_matching_pass = true; }

				}				
			}
		}
		
		else if(!trig_matching_pass && get<0>(pass_hlt[0])==1){
			
			if(Jets.size()>0){
				for(unsigned iobj=0; iobj<trigobj.size(); iobj++){
				
					float ptvar = fabs(Jets[0].pt - trigobj[iobj].pt)/Jets[0].pt;
					if(ptvar<ptvar_min && delta2R(Jets[0].eta,Jets[0].phi,trigobj[iobj].eta,trigobj[iobj].phi)<dR_min && abs(trigobj[iobj].pdgId)!=13) { trig_matching_pass = true; }

				}				
			}
		}
			
		else if(!trig_matching_pass && get<0>(pass_hlt[0])==2){
			
			if( LJets.size()>0){
				for(unsigned iobj=0; iobj<trigobj.size(); iobj++){
				
					float ptvar = fabs(LJets[0].pt - trigobj[iobj].pt)/LJets[0].pt;
					if(ptvar<ptvar_min && delta2R(LJets[0].eta,LJets[0].phi,trigobj[iobj].eta,trigobj[iobj].phi)<dR_min && abs(trigobj[iobj].pdgId)!=13) { trig_matching_pass = true; }

				}				
			}
		}
		
	}
	
	// end of trigger object matching //
	
	for(unsigned ihlt=0; ihlt<(pass_hlt.size()); ihlt++){
		if(get<0>(pass_hlt[ihlt])==0){
			if(single_pids[get<1>(pass_hlt[ihlt])]==13) {  muon_trig_pass = true;  }
			if(single_pids[get<1>(pass_hlt[ihlt])]==11) {  electron_trig_pass = true;  }
		}
	}
}

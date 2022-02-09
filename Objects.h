#include <string>
#include <cmath>
#include <cassert>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <TLorentzVector.h>

#include "Functions.h"

using namespace std;

class AK4Jet {

 public:
  
  int  jetID;
  int  jetID_tightlepveto;
  
  float  pt;
  float  eta;
  float  mass;
  float  phi;
  float  y;
  
  int    hadronFlavour;                                                                        
  int    partonFlavour;
  float  PUID;
  float  qgl;
  
  bool   closebymu;
  bool   closebyel;
  
  float  btag_DeepFlav;
  float  btag_DeepFlav_SF;
  float  btag_DeepFlav_SF_up;
  float  btag_DeepFlav_SF_dn;
  
  float  btag_DeepCSV;
  float  btag_DeepCSV_SF;
  float  btag_DeepCSV_SF_up;
  float  btag_DeepCSV_SF_dn;
 
  float  JER;
  float  JERup;
  float  JERdn;
  
  float JEC;
  float jesup_AbsoluteStat;
  float jesup_AbsoluteScale;
  float jesup_AbsoluteMPFBias;
  float jesup_FlavorQCD;
  float jesup_Fragmentation;
  float jesup_PileUpDataMC;
  float jesup_PileUpPtBB;
  float jesup_PileUpPtEC1;
  float jesup_PileUpPtEC2;
  float jesup_PileUpPtRef;
  float jesup_RelativeFSR;
  float jesup_RelativeJEREC1;
  float jesup_RelativeJEREC2;
  float jesup_RelativePtBB;
  float jesup_RelativePtEC1;
  float jesup_RelativePtEC2;
  float jesup_RelativeBal;
  float jesup_RelativeSample;
  float jesup_RelativeStatEC;
  float jesup_RelativeStatFSR;
  float jesup_SinglePionECAL;
  float jesup_SinglePionHCAL;
  float jesup_TimePtEta;
  float jesup_Total;
  float jesdn_AbsoluteStat;
  float jesdn_AbsoluteScale;
  float jesdn_AbsoluteMPFBias;
  float jesdn_FlavorQCD;
  float jesdn_Fragmentation;
  float jesdn_PileUpDataMC;
  float jesdn_PileUpPtBB;
  float jesdn_PileUpPtEC1;
  float jesdn_PileUpPtEC2;
  float jesdn_PileUpPtRef;
  float jesdn_RelativeFSR;
  float jesdn_RelativeJEREC1;
  float jesdn_RelativeJEREC2;
  float jesdn_RelativePtBB;
  float jesdn_RelativePtEC1;
  float jesdn_RelativePtEC2;
  float jesdn_RelativeBal;
  float jesdn_RelativeSample;
  float jesdn_RelativeStatEC;
  float jesdn_RelativeStatFSR;
  float jesdn_SinglePionECAL;
  float jesdn_SinglePionHCAL;
  float jesdn_TimePtEta;
  float jesdn_Total;
  
  TLorentzVector p4;
};

class AK8Jet {

 public:
  
  int  jetID;
  int  jetID_tightlepveto;
  
  float  pt;
  float  eta;
  float  mass;
  float  phi;
  float  y;
  
  float CHF;
  float NHF;
  float CEMF;
  float NEMF;
  float MUF;
  float PHF;
  float EEF;
  float HFHF;
  int CHM;
  int NHM;
  int MUM;
  int PHM;
  int EEM;
  int HFHM;
  int Neucons;
  int Chcons;

  float sdmass;
  float tau21;
  float tau32;
  
  float btag_DeepCSV;
  float DeepTag_DAK8MD_TvsQCD;
  float DeepTag_DAK8MD_WvsQCD;
  float DeepTag_DAK8MD_ZvsQCD;
  float DeepTag_DAK8MD_HvsQCD;
  float DeepTag_DAK8MD_bbvsQCD;
  float DeepTag_PNet_TvsQCD;
  float DeepTag_PNet_WvsQCD;
  float DeepTag_PNet_ZvsQCD;
  float DeepTag_PNetMD_XbbvsQCD;
  float DeepTag_PNetMD_XccvsQCD;
  float DeepTag_PNetMD_XqqvsQCD;
  float DeepTag_PNetMD_QCD;
  float DeepTag_PNetMD_WvsQCD;
  
  float sub1pt;
  float sub1eta;
  float sub1phi;
  float sub1mass;
  float sub1btag;
  float sub1JEC;
 
  float sub2pt;
  float sub2eta;
  float sub2phi;
  float sub2mass;
  float sub2btag;
  float sub2JEC;
  
  int match_muon_index;
  int match_electron_index;
  int match_lepton_index;
  int match_AK4_index;
  float matchAK4deepb;

  bool haselectron;
  bool hasmuon;
  bool hastau;
  bool hasqg; 
  bool hasb;
  bool hasleptop; 
  bool hashadtop;
  bool hastop;
  bool hasleptop_alldecay;
  bool hashadtop_alldecay;
  bool haspfelectron;
  bool haspfmuon;
  bool hasmatchmu;
  bool hasmatche;
  
  float  JER;
  float  JERup;
  float  JERdn;
  
  float JEC;
  float jesup_AbsoluteStat;
  float jesup_AbsoluteScale;
  float jesup_AbsoluteMPFBias;
  float jesup_FlavorQCD;
  float jesup_Fragmentation;
  float jesup_PileUpDataMC;
  float jesup_PileUpPtBB;
  float jesup_PileUpPtEC1;
  float jesup_PileUpPtEC2;
  float jesup_PileUpPtRef;
  float jesup_RelativeFSR;
  float jesup_RelativeJEREC1;
  float jesup_RelativeJEREC2;
  float jesup_RelativePtBB;
  float jesup_RelativePtEC1;
  float jesup_RelativePtEC2;
  float jesup_RelativeBal;
  float jesup_RelativeSample;
  float jesup_RelativeStatEC;
  float jesup_RelativeStatFSR;
  float jesup_SinglePionECAL;
  float jesup_SinglePionHCAL;
  float jesup_TimePtEta;
  float jesup_Total;
  float jesdn_AbsoluteStat;
  float jesdn_AbsoluteScale;
  float jesdn_AbsoluteMPFBias;
  float jesdn_FlavorQCD;
  float jesdn_Fragmentation;
  float jesdn_PileUpDataMC;
  float jesdn_PileUpPtBB;
  float jesdn_PileUpPtEC1;
  float jesdn_PileUpPtEC2;
  float jesdn_PileUpPtRef;
  float jesdn_RelativeFSR;
  float jesdn_RelativeJEREC1;
  float jesdn_RelativeJEREC2;
  float jesdn_RelativePtBB;
  float jesdn_RelativePtEC1;
  float jesdn_RelativePtEC2;
  float jesdn_RelativeBal;
  float jesdn_RelativeSample;
  float jesdn_RelativeStatEC;
  float jesdn_RelativeStatFSR;
  float jesdn_SinglePionECAL;
  float jesdn_SinglePionHCAL;
  float jesdn_TimePtEta;
  float jesdn_Total;
  
  TLorentzVector p4;
};


class Muon {

 public:

  float pt;
  float ptErr;
  float eta;
  float phi;
  float  mass;
  float charge;
  float p;
  
  float dxy;
  float dxyErr;
  float dz;  
  float ip3d;
  float dxy_sv;
 
  bool ip;
  bool isTRK;
  bool isGL;
  bool isPF;
  bool isLoose;
  bool isGoodGL;
  bool isMed;
  bool isMedPr;
  bool isTight;
  bool isHighPt;
  bool isHighPttrk;
  bool TightID;
  
  float minisoall;
  float minisoch;
  float minisonh;
  float minisoph;
  
  float chi;
  int ndf;
  float ecal;
  float hcal;
  float posmatch;
  float trkink;
  float segcom;
  float hit;
  float mst;
  float pixhit;
  float trklay;
  float valfrac;
  float pfiso;
  
  TLorentzVector p4;
};

class Electron {

 public:

  float pt;
  float ptErr;
  float eta;
  float phi;
  float  mass;
  float e_ECAL;
  float charge;
  float p;
    
  bool Fallv2WP80;
  bool Fallv2WP80_noIso;
  bool Fallv2WP90;
  bool Fallv2WP90_noIso;
  
  float dxy;
  float dxyErr;
  float dz;
  float dzErr;
  float ip3d;
  float dxy_sv;
  
  bool ip;
 
  float hcaloverecal;
  float cloctftrkn;
  float cloctftrkchi2;
  float e1x5bye5x5;
  float etain;
  float phiin;
  float fbrem;
  float eoverp;
  float hovere;
  float chi;
  int npdf;
  float ietaieta;
  float pfiso_drcor;
  float pfiso_eacor;
  float pfiso04_eacor;
  bool convVeto;
  float ecloverpout;
  float misshits;
  
  //float normchi2;
  // float trkmeasure;
  //float ecaletrkmomentum;
  //float deltaetacltrkcalo;
  //float supcl_preshvsrawe;
 
  float eccalTrkEnergyPostCorr;
  float energyScaleValue;
  float energyScaleUp;
  float energyScaleDown;
  float energySigmaValue;
  float energySigmaUp;
  float energySigmaDown;
   
  float supcl_eta;
  float supcl_phi;
  float supcl_e;
  float supcl_rawE;
  float sigmaieta;
  float sigmaiphi;
  float r9full;
  float supcl_etaw;
  float supcl_phiw;
 
  float pfisolsumphet;
  float pfisolsumchhadpt;
  float pfsiolsumneuhadet;
  float minisoall;
  float minisoch;
  float minisonh;
  float minisoph;
  
  TLorentzVector p4;

};

class Lepton {
 
 public:

  float pt;
  float eta;
  float phi;
  float mass;
  float charge;
  int lepton_id;
  int indexemu;
  int pdgId;
  int AK8_neighbor_index;
  
  TLorentzVector p4;
};

class AK4GenJet {

 public:

  float  eta;
  float  mass;
  float  phi;
  float  pt;
  int hadronFlavor;
  int partonFlavor;

  TLorentzVector p4;

};

class AK8GenJet {

 public:

  float  eta;
  float  mass;
  float  phi;
  float  pt;
  int hadronFlavor;
  int partonFlavor;

  TLorentzVector p4;

} ;

class GenParton{

 public:

  float  eta;
  float  mass;
  float  phi;
  float  pt;

  int status;
  int pdgId;
  int mompdgId;
  int grmompdgId;
  
  bool fromhard;
  bool fromhardbFSR;
  bool isPromptFinalState;
  bool isLastCopyBeforeFSR;
  
  TLorentzVector p4;

} ;

class TopQuark{
 // gives 4-momentum of top quark and a vector of its daughters 
 //(length of the vector of daughters should be 3)
 // each daughter is of GenParton-type
 public:

  TLorentzVector p4;
  vector<GenParton> daughter;

} ;

bool AK4Jet_sort_by_pt(AK4Jet i1, AK4Jet i2)
{
  return (i1.pt > i2.pt);
}
void sorted_by_pt(vector<AK4Jet> & objs) {
  sort(objs.begin(), objs.end(), AK4Jet_sort_by_pt);
}
bool AK8Jet_sort_by_pt(AK8Jet i1, AK8Jet i2)                                                        
{                                                                                                   
  return (i1.pt > i2.pt);                                                                           
}                                                                                                   
void sorted_by_pt(vector<AK8Jet> & objs) {                                                          
  sort(objs.begin(), objs.end(), AK8Jet_sort_by_pt);                                                
}

bool AK4GenJet_sort_by_pt(AK4GenJet i1, AK4GenJet i2)
{                        
  return (i1.pt > i2.pt);
}
void sorted_by_pt(vector<AK4GenJet> & objs) {
  sort(objs.begin(), objs.end(), AK4GenJet_sort_by_pt);
}

bool Muon_sort_by_pt(Muon i1, Muon i2)                                                                           
{                                                                                                                
  return (i1.pt > i2.pt);                                                                                        
}                                                                                                                
void sorted_by_pt(vector<Muon> & objs) {                                                                         
  sort(objs.begin(), objs.end(), Muon_sort_by_pt);                                                               
}
bool Electron_sort_by_pt(Electron i1, Electron i2)
{
  return (i1.pt > i2.pt);
}
void sorted_by_pt(vector<Electron> & objs) {
  sort(objs.begin(), objs.end(), Electron_sort_by_pt);
}
bool Lepton_sort_by_pt(Lepton i1, Lepton i2)
{
  return (i1.pt > i2.pt);
}
void sorted_by_pt(vector<Lepton> & objs) {
  sort(objs.begin(), objs.end(), Lepton_sort_by_pt);
}

bool Parton_sort_by_pt(GenParton i1, GenParton i2)
{
  return (i1.pt > i2.pt);
}
void sorted_by_pt(vector<GenParton> & objs) {
  sort(objs.begin(), objs.end(), Parton_sort_by_pt);
}

bool Top_sort_by_pt(TopQuark i1, TopQuark i2)
{
  return (i1.p4.Pt() > i2.p4.Pt());
}
void sorted_by_pt(vector<TopQuark> & objs) {
  sort(objs.begin(), objs.end(), Top_sort_by_pt);
}


float compute_HT(vector<AK4Jet>  & objs, float ptcut, float etacut){
  float HT = 0;
  for(unsigned iobs=0; iobs<objs.size(); iobs++){
    if(objs[iobs].pt > ptcut && abs(objs[iobs].eta)<=etacut){
  HT += objs[iobs].pt;
    }
  }
  return HT;
}


int get_nearest_AK4(vector<AK4Jet>  & objs, TLorentzVector tmp_vec, float minR = 0.6) {
	// gives the index of AK4 jet nearest to the tmp_vec vector
    int nearest = -1;

    for(unsigned iobs=0; iobs<objs.size(); iobs++){

		if(delta2R(objs[iobs].eta,objs[iobs].phi,tmp_vec.Eta(),tmp_vec.Phi()) < minR){

			minR = delta2R(objs[iobs].eta,objs[iobs].phi,tmp_vec.Eta(),tmp_vec.Phi()) ;
            nearest = iobs;

           }
	}

    return  nearest;
}

int get_nearest_lepton(vector<Lepton>  & objs, TLorentzVector tmp_vec, int lepton_pdgid=-1, float minR = 0.8) {
    // gives the index of lepton nearest to the tmp_vec vector
	int nearest = -1;

    for(unsigned iobs=0; iobs<objs.size(); iobs++){

		if(lepton_pdgid>=0 && objs[iobs].pdgId!=lepton_pdgid) continue;

		if(delta2R(objs[iobs].eta,objs[iobs].phi,tmp_vec.Eta(),tmp_vec.Phi()) < minR){

			minR = delta2R(objs[iobs].eta,objs[iobs].phi,tmp_vec.Eta(),tmp_vec.Phi()) ;
            nearest = iobs;

           }
	}

    return  nearest;
}

int get_nearest_AK8Jet(vector<AK8Jet>  & objs, TLorentzVector tmp_vec, float minR = 0.8) {
    // gives the index of AK8 jet nearest to the tmp_vec vector
	int nearest = -1;

    for(unsigned iobs=0; iobs<objs.size(); iobs++){

		if(delta2R(objs[iobs].y,objs[iobs].phi,tmp_vec.Rapidity(),tmp_vec.Phi()) < minR){

			minR = delta2R(objs[iobs].y,objs[iobs].phi,tmp_vec.Rapidity(),tmp_vec.Phi()) ;
            nearest = iobs;

           }
	}

    return  nearest;
}

int get_nearest_Parton(vector<GenParton>  & objs, TLorentzVector tmp_vec, float minR = 0.4) {
    // gives the index of parton to the tmp_vec vector
	int nearest = -1;

    for(unsigned iobs=0; iobs<objs.size(); iobs++){

		if(delta2R(objs[iobs].p4.Rapidity(),objs[iobs].phi,tmp_vec.Rapidity(),tmp_vec.Phi()) < minR){

			minR = delta2R(objs[iobs].p4.Rapidity(),objs[iobs].phi,tmp_vec.Rapidity(),tmp_vec.Phi()) ;
            nearest = iobs;

           }
	}

    return  nearest;
}


/*
  bool AK4Jet_sort_by_DeepFlav(AK4Jet i1, AK4Jet i2)
  {
  return (i1.btagDeepFlavB > i2.btagDeepFlavB);
  }
  void sorted_by_DeepFlav(vector<AK4Jet> & objs) {
  sort(objs.begin(), objs.end(), AK4Jet_sort_by_DeepFlav);
  }
  bool AK8Jet_sort_by_DeepAK8_Htag(AK8Jet i1, AK8Jet i2)
  {
  return (i1.deepTagMD_bbvsLight > i2.deepTagMD_bbvsLight);
  }
  void sorted_by_DeepAK8_Htag(vector<AK8Jet> & objs) {
  sort(objs.begin(), objs.end(), AK8Jet_sort_by_DeepAK8_Htag);
  }
  bool GenAK4Jet_sort_by_pt(AK4GenJet i1, AK4GenJet i2)
  {
  return (i1.pt > i2.pt);
  }
  void sorted_by_pt(vector<AK4GenJet> & objs) {
  sort(objs.begin(), objs.end(), GenAK4Jet_sort_by_pt);
  }
  bool GenAK8Jet_sort_by_pt(AK8GenJet i1, AK8GenJet i2)
  {
  return (i1.pt > i2.pt);
  }
  void sorted_by_pt(vector<AK8GenJet> & objs) {
  sort(objs.begin(), objs.end(), GenAK8Jet_sort_by_pt);
  }
*/

class TrigObj {
 
 public:

  float pt;
  float eta;
  float phi;
  float mass;
  bool HLT;
  bool L1;
  int pdgId;
  int ID;
  int type;
  
  TLorentzVector p4;
};

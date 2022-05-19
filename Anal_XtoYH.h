//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct 29 13:04:45 2021 by ROOT version 6.20/08
// from TTree Events/XtoYH
// found on file: rootuple_eg.root
//////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <math.h>

#include <TRandom3.h>

#include "BTagCalibrationStandalone.h"
#include "BTagCalibrationStandalone.cpp"

#include "Objects.h"

   static const int njetmx = 20; 
   static const int njetmxAK8 =10;
   static const int npartmx = 50; 
   static const int nlhemax = 10;
   
   static const int nlhescalemax = 9;
   static const int nlhepdfmax = 103;
   static const int nalpsmax = 3;
   static const int nlhepsmax = 8;

   // Declaration of leaf types
   Int_t           irun;
   Int_t           ilumi;
   UInt_t          ievt;
   Int_t           nprim;
   Int_t           npvert;
   Double_t        Rho;
   Int_t           trig_value;
   Bool_t          hlt_IsoMu24;
   Bool_t          hlt_Mu50;
   Bool_t          hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;
   Bool_t          hlt_Ele115_CaloIdVT_GsfTrkIdT;
   Bool_t          hlt_Ele40_WPTight_Gsf;
   Bool_t          hlt_Ele32_WPTight_Gsf;
   Bool_t          hlt_Ele28_eta2p1_WPTight_Gsf_HT150;
   Bool_t          hlt_Mu37_Ele27_CaloIdL_MW;
   Bool_t          hlt_Mu27_Ele37_CaloIdL_MW;
   Bool_t          hlt_Mu37_TkMu27;
   Bool_t          hlt_DoubleEle25_CaloIdL_MW;
   Bool_t          hlt_AK8PFJet500;
   Bool_t          hlt_PFJet500;
   Bool_t          hlt_HT1050;
   Bool_t          hlt_AK8PFJet400_TrimMass30;
   Bool_t          hlt_AK8PFHT800_TrimMass50;
   Bool_t          hlt_Photon200;
   Bool_t          hlt_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
   Bool_t          hlt_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;
   Bool_t          hlt_PFMETNoMu140_PFMHTNoMu140_IDTight;
   Bool_t          hlt_PFMETTypeOne140_PFMHT140_IDTight;
   Int_t           nTrigObj;
   int ncuts;
   Bool_t          Flag_event_cuts[njetmxAK8];
   Float_t         TrigObj_pt[njetmx];   //[nTrigObj]
   Float_t         TrigObj_eta[njetmx];   //[nTrigObj]
   Float_t         TrigObj_phi[njetmx];   //[nTrigObj]
   Float_t         TrigObj_mass[njetmx];   //[nTrigObj]
   Bool_t          TrigObj_HLT[njetmx];   //[nTrigObj]
   Bool_t          TrigObj_L1[njetmx];   //[nTrigObj]
   Int_t           TrigObj_Ihlt[njetmx];   //[nTrigObj]
   Int_t           TrigObj_pdgId[njetmx];   //[nTrigObj]
   Int_t           TrigObj_type[njetmx];   //[nTrigObj]
   Double_t        prefiringweight;
   Double_t        prefiringweightup;
   Double_t        prefiringweightdown;
   Float_t         CHSMET_pt;
   Float_t         CHSMET_phi;
   Float_t         CHSMET_sig;
   Float_t         CHSMET_sumEt;
   Float_t         PuppiMET_pt;
   Float_t         PuppiMET_phi;
   Float_t         PuppiMET_sig;
   Float_t         PuppiMET_sumEt;
   Float_t         PuppiMET_pt_JESup;
   Float_t         PuppiMET_pt_JESdn;
   Float_t         PuppiMET_pt_JERup;
   Float_t         PuppiMET_pt_JERdn;
   Float_t         PuppiMET_pt_UnclusEup;
   Float_t         PuppiMET_pt_UnclusEdn;
   Float_t         PuppiMET_phi_JESup;
   Float_t         PuppiMET_phi_JESdn;
   Float_t         PuppiMET_phi_JERup;
   Float_t         PuppiMET_phi_JERdn;
   Float_t         PuppiMET_phi_UnclusEup;
   Float_t         PuppiMET_phi_UnclusEdn;
   Int_t           nPFJetAK8;
   Float_t         PFJetAK8_pt[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_y[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_eta[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_phi[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_mass[njetmxAK8];   //[nPFJetAK8]
   Bool_t          PFJetAK8_jetID_tightlepveto[njetmxAK8];   //[nPFJetAK8]
   Bool_t          PFJetAK8_jetID[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_JEC[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_CHF[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_NHF[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_CEMF[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_NEMF[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_MUF[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_PHF[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_EEF[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_HFHF[njetmxAK8];   //[nPFJetAK8]
   Int_t           PFJetAK8_CHM[njetmxAK8];   //[nPFJetAK8]
   Int_t           PFJetAK8_NHM[njetmxAK8];   //[nPFJetAK8]
   Int_t           PFJetAK8_MUM[njetmxAK8];   //[nPFJetAK8]
   Int_t           PFJetAK8_PHM[njetmxAK8];   //[nPFJetAK8]
   Int_t           PFJetAK8_EEM[njetmxAK8];   //[nPFJetAK8]
   Int_t           PFJetAK8_HFHM[njetmxAK8];   //[nPFJetAK8]
   Int_t           PFJetAK8_Neucons[njetmxAK8];   //[nPFJetAK8]
   Int_t           PFJetAK8_Chcons[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_JER[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_JERup[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_JERdn[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_msoftdrop[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_tau1[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_tau2[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_tau3[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_btag_DeepCSV[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_DAK8MD_TvsQCD[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_DAK8MD_WvsQCD[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_DAK8MD_ZvsQCD[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_DAK8MD_HvsQCD[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_DAK8MD_bbvsQCD[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_PNet_TvsQCD[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_PNet_WvsQCD[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_PNet_ZvsQCD[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_PNetMD_XbbvsQCD[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_PNetMD_XccvsQCD[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_PNetMD_XqqvsQCD[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_DeepTag_PNetMD_QCD[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_sub1pt[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_sub1eta[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_sub1phi[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_sub1mass[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_sub1btag[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_sub1JEC[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_sub2pt[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_sub2eta[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_sub2phi[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_sub2mass[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_sub2btag[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_sub2JEC[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_AbsoluteStat[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_AbsoluteScale[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_AbsoluteMPFBias[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_FlavorQCD[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_Fragmentation[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_PileUpDataMC[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_PileUpPtBB[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_PileUpPtEC1[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_PileUpPtEC2[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_PileUpPtRef[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_RelativeFSR[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_RelativeJEREC1[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_RelativeJEREC2[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_RelativePtBB[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_RelativePtEC1[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_RelativePtEC2[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_RelativeBal[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_RelativeSample[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_RelativeStatEC[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_RelativeStatFSR[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_SinglePionECAL[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_SinglePionHCAL[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_TimePtEta[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesup_Total[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_AbsoluteStat[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_AbsoluteScale[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_AbsoluteMPFBias[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_FlavorQCD[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_Fragmentation[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_PileUpDataMC[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_PileUpPtBB[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_PileUpPtEC1[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_PileUpPtEC2[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_PileUpPtRef[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_RelativeFSR[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_RelativeJEREC1[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_RelativeJEREC2[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_RelativePtBB[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_RelativePtEC1[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_RelativePtEC2[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_RelativeBal[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_RelativeSample[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_RelativeStatEC[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_RelativeStatFSR[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_SinglePionECAL[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_SinglePionHCAL[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_TimePtEta[njetmxAK8];   //[nPFJetAK8]
   Float_t         PFJetAK8_jesdn_Total[njetmxAK8];   //[nPFJetAK8]
   Int_t           nPFJetAK4;
   Bool_t          PFJetAK4_jetID[njetmx];   //[nPFJetAK4]
   Bool_t          PFJetAK4_jetID_tightlepveto[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_pt[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_eta[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_y[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_phi[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_mass[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_JEC[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_btag_DeepCSV[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_btag_DeepFlav[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_JER[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_JERup[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_JERdn[njetmx];   //[nPFJetAK4]
   Int_t           PFJetAK4_hadronflav[njetmx];   //[nPFJetAK4]
   Int_t           PFJetAK4_partonflav[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_qgl[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_PUID[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_AbsoluteStat[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_AbsoluteScale[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_AbsoluteMPFBias[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_FlavorQCD[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_Fragmentation[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_PileUpDataMC[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_PileUpPtBB[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_PileUpPtEC1[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_PileUpPtEC2[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_PileUpPtRef[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_RelativeFSR[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_RelativeJEREC1[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_RelativeJEREC2[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_RelativePtBB[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_RelativePtEC1[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_RelativePtEC2[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_RelativeBal[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_RelativeSample[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_RelativeStatEC[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_RelativeStatFSR[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_SinglePionECAL[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_SinglePionHCAL[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_TimePtEta[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesup_Total[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_AbsoluteStat[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_AbsoluteScale[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_AbsoluteMPFBias[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_FlavorQCD[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_Fragmentation[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_PileUpDataMC[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_PileUpPtBB[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_PileUpPtEC1[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_PileUpPtEC2[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_PileUpPtRef[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_RelativeFSR[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_RelativeJEREC1[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_RelativeJEREC2[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_RelativePtBB[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_RelativePtEC1[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_RelativePtEC2[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_RelativeBal[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_RelativeSample[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_RelativeStatEC[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_RelativeStatFSR[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_SinglePionECAL[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_SinglePionHCAL[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_TimePtEta[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_jesdn_Total[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_btag_DeepCSV_SF[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_btag_DeepCSV_SF_up[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_btag_DeepCSV_SF_dn[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_btag_DeepFlav_SF[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_btag_DeepFlav_SF_up[njetmx];   //[nPFJetAK4]
   Float_t         PFJetAK4_btag_DeepFlav_SF_dn[njetmx];   //[nPFJetAK4]
   Int_t           nMuon;
   Bool_t          Muon_isPF[njetmx];   //[nMuon]
   Bool_t          Muon_isGL[njetmx];   //[nMuon]
   Bool_t          Muon_isTRK[njetmx];   //[nMuon]
   Bool_t          Muon_isLoose[njetmx];   //[nMuon]
   Bool_t          Muon_isGoodGL[njetmx];   //[nMuon]
   Bool_t          Muon_isMed[njetmx];   //[nMuon]
   Bool_t          Muon_isMedPr[njetmx];   //[nMuon]
   Bool_t          Muon_isTight[njetmx];   //[nMuon]
   Bool_t          Muon_isHighPt[njetmx];   //[nMuon]
   Bool_t          Muon_isHighPttrk[njetmx];   //[nMuon]
   Bool_t          Muon_TightID[njetmx];   //[nMuon]
   Float_t         Muon_pt[njetmx];   //[nMuon]
   Float_t         Muon_p[njetmx];   //[nMuon]
   Float_t         Muon_eta[njetmx];   //[nMuon]
   Float_t         Muon_phi[njetmx];   //[nMuon]
   Float_t         Muon_minisoch[njetmx];   //[nMuon]
   Float_t         Muon_minisonh[njetmx];   //[nMuon]
   Float_t         Muon_minisoph[njetmx];   //[nMuon]
   Float_t         Muon_minisoall[njetmx];   //[nMuon]
   Float_t         Muon_dxy[njetmx];   //[nMuon]
   Float_t         Muon_dz[njetmx];   //[nMuon]
   Float_t         Muon_dxyErr[njetmx];   //[nMuon]
   Float_t         Muon_ip3d[njetmx];   //[nMuon]
   Float_t         Muon_ptErr[njetmx];   //[nMuon]
   Float_t         Muon_chi[njetmx];   //[nMuon]
   Int_t           Muon_ndf[njetmx];   //[nMuon]
   Float_t         Muon_ecal[njetmx];   //[nMuon]
   Float_t         Muon_hcal[njetmx];   //[nMuon]
   Float_t         Muon_pfiso[njetmx];   //[nMuon]
   Float_t         Muon_posmatch[njetmx];   //[nMuon]
   Float_t         Muon_trkink[njetmx];   //[nMuon]
   Float_t         Muon_segcom[njetmx];   //[nMuon]
   Float_t         Muon_hit[njetmx];   //[nMuon]
   Float_t         Muon_pixhit[njetmx];   //[nMuon]
   Float_t         Muon_mst[njetmx];   //[nMuon]
   Float_t         Muon_trklay[njetmx];   //[nMuon]
   Float_t         Muon_valfrac[njetmx];   //[nMuon]
   Float_t         Muon_dxy_sv[njetmx];   //[nMuon]
   Int_t           nElectron;
   Float_t         Electron_pt[njetmx];   //[nElectron]
   Float_t         Electron_eta[njetmx];   //[nElectron]
   Float_t         Electron_phi[njetmx];   //[nElectron]
   Float_t         Electron_p[njetmx];   //[nElectron]
   Float_t         Electron_e[njetmx];   //[nElectron]
   Float_t         Electron_e_ECAL[njetmx];   //[nElectron]
   Bool_t          Electron_mvaid_Fallv2WP90[njetmx];   //[nElectron]
   Bool_t          Electron_mvaid_Fallv2WP90_noIso[njetmx];   //[nElectron]
   Bool_t          Electron_mvaid_Fallv2WP80[njetmx];   //[nElectron]
   Bool_t          Electron_mvaid_Fallv2WP80_noIso[njetmx];   //[nElectron]
   Float_t         Electron_dxy[njetmx];   //[nElectron]
   Float_t         Electron_dxyErr[njetmx];   //[nElectron]
   Float_t         Electron_dz[njetmx];   //[nElectron]
   Float_t         Electron_dzErr[njetmx];   //[nElectron]
   Float_t         Electron_ip3d[njetmx];   //[nElectron]
   Float_t         Electron_dxy_sv[njetmx];   //[nElectron]
   Float_t         Electron_hovere[njetmx];   //[nElectron]
   Float_t         Electron_chi[njetmx];   //[nElectron]
   Int_t           Electron_ndf[njetmx];   //[nElectron]
   Float_t         Electron_eoverp[njetmx];   //[nElectron]
   Float_t         Electron_ietaieta[njetmx];   //[nElectron]
   Float_t         Electron_misshits[njetmx];   //[nElectron]
   Float_t         Electron_pfiso_drcor[njetmx];   //[nElectron]
   Float_t         Electron_pfiso_eacor[njetmx];   //[nElectron]
   Float_t         Electron_pfiso04_eacor[njetmx];   //[nElectron]
   Float_t         Electron_eccalTrkEnergyPostCorr[njetmx];   //[nElectron]
   Float_t         Electron_energyScaleValue[njetmx];   //[nElectron]
   Float_t         Electron_energyScaleUp[njetmx];   //[nElectron]
   Float_t         Electron_energyScaleDown[njetmx];   //[nElectron]
   Float_t         Electron_energySigmaValue[njetmx];   //[nElectron]
   Float_t         Electron_energySigmaUp[njetmx];   //[nElectron]
   Float_t         Electron_energySigmaDown[njetmx];   //[nElectron]
   Float_t         Electron_supcl_eta[njetmx];   //[nElectron]
   Float_t         Electron_supcl_phi[njetmx];   //[nElectron]
   Float_t         Electron_supcl_e[njetmx];   //[nElectron]
   Float_t         Electron_supcl_rawE[njetmx];   //[nElectron]
   Float_t         Electron_sigmaieta[njetmx];   //[nElectron]
   Float_t         Electron_sigmaiphi[njetmx];   //[nElectron]
   Float_t         Electron_r9full[njetmx];   //[nElectron]
   Float_t         Electron_hcaloverecal[njetmx];   //[nElectron]
   Float_t         Electron_hitsmiss[njetmx];   //[nElectron]
   Float_t         Electron_ecloverpout[njetmx];   //[nElectron]
   Bool_t          Electron_convVeto[njetmx];   //[nElectron]
   Float_t         Electron_pfisolsumphet[njetmx];   //[nElectron]
   Float_t         Electron_pfisolsumchhadpt[njetmx];   //[nElectron]
   Float_t         Electron_pfsiolsumneuhadet[njetmx];   //[nElectron]
   Float_t         Electron_minisoch[njetmx];   //[nElectron]
   Float_t         Electron_minisonh[njetmx];   //[nElectron]
   Float_t         Electron_minisoph[njetmx];   //[nElectron]
   Float_t         Electron_minisoall[njetmx];   //[nElectron]
   Int_t           nPhoton;
   Float_t         Photon_e[njetmx];   //[nPhoton]
   Float_t         Photon_eta[njetmx];   //[nPhoton]
   Float_t         Photon_phi[njetmx];   //[nPhoton]
   Float_t         Photon_mvaid_Fall17V2_raw[njetmx];   //[nPhoton]
   Bool_t          Photon_mvaid_Fall17V2_WP90[njetmx];   //[nPhoton]
   Bool_t          Photon_mvaid_Fall17V2_WP80[njetmx];   //[nPhoton]
   Bool_t          Photon_mvaid_Spring16V1_WP90[njetmx];   //[nPhoton]
   Bool_t          Photon_mvaid_Spring16V1_WP80[njetmx];   //[nPhoton]
   Float_t         Photon_e1by9[njetmx];   //[nPhoton]
   Float_t         Photon_e9by25[njetmx];   //[nPhoton]
   Float_t         Photon_trkiso[njetmx];   //[nPhoton]
   Float_t         Photon_emiso[njetmx];   //[nPhoton]
   Float_t         Photon_hadiso[njetmx];   //[nPhoton]
   Float_t         Photon_chhadiso[njetmx];   //[nPhoton]
   Float_t         Photon_neuhadiso[njetmx];   //[nPhoton]
   Float_t         Photon_phoiso[njetmx];   //[nPhoton]
   Float_t         Photon_PUiso[njetmx];   //[nPhoton]
   Float_t         Photon_hadbyem[njetmx];   //[nPhoton]
   Float_t         Photon_ietaieta[njetmx];   //[nPhoton]
   Int_t           nTau;
   Bool_t          Tau_isPF[njetmx];   //[nTau]
   Float_t         Tau_pt[njetmx];   //[nTau]
   Float_t         Tau_eta[njetmx];   //[nTau]
   Float_t         Tau_phi[njetmx];   //[nTau]
   Float_t         Tau_e[njetmx];   //[nTau]
   Int_t           Tau_charge[njetmx];   //[nTau]
   Float_t         Tau_dxy[njetmx];   //[nTau]
   Float_t         Tau_leadtrkdxy[njetmx];   //[nTau]
   Float_t         Tau_leadtrkdz[njetmx];   //[nTau]
   Float_t         Tau_leadtrkpt[njetmx];   //[nTau]
   Float_t         Tau_leadtrketa[njetmx];   //[nTau]
   Float_t         Tau_leadtrkphi[njetmx];   //[nTau]
   Int_t           Tau_decayMode[njetmx];   //[nTau]
   Bool_t          Tau_decayModeinding[njetmx];   //[nTau]
   Bool_t          Tau_decayModeindingNewDMs[njetmx];   //[nTau]
   Float_t         Tau_eiso2018_raw[njetmx];   //[nTau]
   Int_t           Tau_eiso2018[njetmx];   //[nTau]
   Float_t         Tau_jetiso_deeptau2017v2p1_raw[njetmx];   //[nTau]
   Int_t           Tau_jetiso_deeptau2017v2p1[njetmx];   //[nTau]
   Float_t         Tau_eiso_deeptau2017v2p1_raw[njetmx];   //[nTau]
   Int_t           Tau_eiso_deeptau2017v2p1[njetmx];   //[nTau]
   Float_t         Tau_muiso_deeptau2017v2p1_raw[njetmx];   //[nTau]
   Int_t           Tau_muiso_deeptau2017v2p1[njetmx];   //[nTau]
   Float_t         Tau_rawiso[njetmx];   //[nTau]
   Float_t         Tau_rawisodR03[njetmx];   //[nTau]
   Float_t         Tau_puCorr[njetmx];   //[nTau]
   Double_t        Generator_weight;
   Float_t         Generator_qscale;
   Float_t         Generator_x1;
   Float_t         Generator_x2;
   Float_t         Generator_xpdf1;
   Float_t         Generator_xpdf2;
   Int_t           Generator_id1;
   Int_t           Generator_id2;
   Float_t         Generator_scalePDF;
   Int_t           npu_vert;
   Int_t           npu_vert_true;
   Float_t         GENMET_pt;
   Float_t         GENMET_phi;
   Int_t           nGenJetAK8;
   Float_t         GenJetAK8_pt[njetmxAK8];   //[nGenJetAK8]
   Float_t         GenJetAK8_eta[njetmxAK8];   //[nGenJetAK8]
   Float_t         GenJetAK8_phi[njetmxAK8];   //[nGenJetAK8]
   Float_t         GenJetAK8_mass[njetmxAK8];   //[nGenJetAK8]
   Float_t         GenJetAK8_sdmass[njetmxAK8];   //[nGenJetAK8]
   Int_t           GenJetAK8_hadronflav[njetmxAK8];   //[nGenJetAK8]
   Int_t           GenJetAK8_partonflav[njetmxAK8];   //[nGenJetAK8]
   Int_t           nGenJetAK4;
   Float_t         GenJetAK4_pt[njetmx];   //[nGenJetAK4]
   Float_t         GenJetAK4_eta[njetmx];   //[nGenJetAK4]
   Float_t         GenJetAK4_phi[njetmx];   //[nGenJetAK4]
   Float_t         GenJetAK4_mass[njetmx];   //[nGenJetAK4]
   Int_t           GenJetAK4_hadronflav[njetmx];   //[nGenJetAK4]
   Int_t           GenJetAK4_partonflav[njetmx];   //[nGenJetAK4]
   Int_t           nGenPart;
   Float_t         GenPart_pt[npartmx];   //[nGenPart]
   Float_t         GenPart_eta[npartmx];   //[nGenPart]
   Float_t         GenPart_phi[npartmx];   //[nGenPart]
   Float_t         GenPart_m[npartmx];   //[nGenPart]
   Int_t           GenPart_status[npartmx];   //[nGenPart]
   Int_t           GenPart_pdgId[npartmx];   //[nGenPart]
   Int_t           GenPart_mompdgId[npartmx];   //[nGenPart]
   Int_t           GenPart_grmompdgId[npartmx];   //[nGenPart]
   Int_t           GenPart_daugno[npartmx];   //[nGenPart]
   Bool_t          GenPart_fromhard[npartmx];   //[nGenPart]
   Bool_t          GenPart_fromhardbFSR[npartmx];   //[nGenPart]
   Bool_t          GenPart_isPromptFinalState[npartmx];   //[nGenPart]
   Bool_t          GenPart_isLastCopyBeforeFSR[npartmx];   //[nGenPart]
   Int_t           nLHEPart;
   Int_t           LHEPart_pdg[nlhemax];   //[nLHEPart]
   Float_t         LHEPart_pt[nlhemax];   //[nLHEPart]
   Float_t         LHEPart_eta[nlhemax];   //[nLHEPart]
   Float_t         LHEPart_phi[nlhemax];   //[nLHEPart]
   Float_t         LHEPart_m[nlhemax];   //[nLHEPart]
   Double_t        LHE_weight;
   Int_t           nLHEScaleWeights;
   Float_t         LHEScaleWeights[nlhescalemax];   //[nLHEScaleWeights]
   Int_t           nLHEPDFWeights;
   Float_t         LHEPDFWeights[nlhepdfmax];   //[nLHEPDFWeights]
   Int_t           nLHEAlpsWeights;
   Float_t         LHEAlpsWeights[nalpsmax];   //[nLHEAlpsWeights]
   Int_t           nLHEPSWeights;
   Float_t         LHEPSWeights[nlhepsmax];   //[nLHEPSWeights]

   // List of branches
   TBranch        *b_irun;   //!
   TBranch        *b_ilumi;   //!
   TBranch        *b_ievt;   //!
   TBranch        *b_nprim;   //!
   TBranch        *b_npvert;   //!
   TBranch        *b_Rho;   //!
   TBranch        *b_trig_value;   //!
   TBranch        *b_hlt_IsoMu24;   //!
   TBranch        *b_hlt_Mu50;   //!
   TBranch        *b_hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;   //!
   TBranch        *b_hlt_Ele115_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_hlt_Ele40_WPTight_Gsf;   //!
   TBranch        *b_hlt_Ele32_WPTight_Gsf;   //!
   TBranch        *b_hlt_Ele28_eta2p1_WPTight_Gsf_HT150;   //!
   TBranch        *b_hlt_Mu37_Ele27_CaloIdL_MW;   //!
   TBranch        *b_hlt_Mu27_Ele37_CaloIdL_MW;   //!
   TBranch        *b_hlt_Mu37_TkMu27;   //!
   TBranch        *b_hlt_DoubleEle25_CaloIdL_MW;   //!
   TBranch        *b_hlt_AK8PFJet500;   //!
   TBranch        *b_hlt_PFJet500;   //!
   TBranch        *b_hlt_HT1050;   //!
   TBranch        *b_hlt_AK8PFJet400_TrimMass30;   //!
   TBranch        *b_hlt_AK8PFHT800_TrimMass50;   //!
   TBranch        *b_hlt_Photon200;   //!
   TBranch        *b_hlt_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;   //!
   TBranch        *b_hlt_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;   //!
   TBranch        *b_hlt_PFMETNoMu140_PFMHTNoMu140_IDTight;   //!
   TBranch        *b_hlt_PFMETTypeOne140_PFMHT140_IDTight;   //!
   TBranch        *b_nTrigObj;   //!
   TBranch        *b_TrigObj_pt;   //!
   TBranch        *b_TrigObj_eta;   //!
   TBranch        *b_TrigObj_phi;   //!
   TBranch        *b_TrigObj_mass;   //!
   TBranch        *b_TrigObj_HLT;   //!
   TBranch        *b_TrigObj_L1;   //!
   TBranch        *b_TrigObj_Ihlt;   //!
   TBranch        *b_TrigObj_pdgId;   //!
   TBranch        *b_TrigObj_type;   //!
   TBranch        *b_prefiringweight;   //!
   TBranch        *b_prefiringweightup;   //!
   TBranch        *b_prefiringweightdown;   //!
   TBranch        *b_miset;   //!
   TBranch        *b_misphi;   //!
   TBranch        *b_misetsig;   //!
   TBranch        *b_sumEt;   //!
   TBranch        *b_miset_PUPPI;   //!
   TBranch        *b_misphi_PUPPI;   //!
   TBranch        *b_misetsig_PUPPI;   //!
   TBranch        *b_sumEt_PUPPI;   //!
   TBranch        *b_miset_PUPPI_JESup;   //!
   TBranch        *b_miset_PUPPI_JESdn;   //!
   TBranch        *b_miset_PUPPI_JERup;   //!
   TBranch        *b_miset_PUPPI_JERdn;   //!
   TBranch        *b_miset_PUPPI_UnclusEup;   //!
   TBranch        *b_miset_PUPPI_UnclusEdn;   //!
   TBranch        *b_misphi_PUPPI_JESup;   //!
   TBranch        *b_misphi_PUPPI_JESdn;   //!
   TBranch        *b_misphi_PUPPI_JERup;   //!
   TBranch        *b_misphi_PUPPI_JERdn;   //!
   TBranch        *b_misphi_PUPPI_UnclusEup;   //!
   TBranch        *b_misphi_PUPPI_UnclusEdn;   //!
   TBranch        *b_nPFJetAK8;   //!
   TBranch        *b_PFJetAK8_pt;   //!
   TBranch        *b_PFJetAK8_y;   //!
   TBranch        *b_PFJetAK8_eta;   //!
   TBranch        *b_PFJetAK8_phi;   //!
   TBranch        *b_PFJetAK8_mass;   //!
   TBranch        *b_PFJetAK8_jetID_tightlepveto;   //!
   TBranch        *b_PFJetAK8_jetID;   //!
   TBranch        *b_PFJetAK8_JEC;   //!
   TBranch        *b_PFJetAK8_CHF;   //!
   TBranch        *b_PFJetAK8_NHF;   //!
   TBranch        *b_PFJetAK8_CEMF;   //!
   TBranch        *b_PFJetAK8_NEMF;   //!
   TBranch        *b_PFJetAK8_MUF;   //!
   TBranch        *b_PFJetAK8_PHF;   //!
   TBranch        *b_PFJetAK8_EEF;   //!
   TBranch        *b_PFJetAK8_HFHF;   //!
   TBranch        *b_PFJetAK8_CHM;   //!
   TBranch        *b_PFJetAK8_NHM;   //!
   TBranch        *b_PFJetAK8_MUM;   //!
   TBranch        *b_PFJetAK8_PHM;   //!
   TBranch        *b_PFJetAK8_EEM;   //!
   TBranch        *b_PFJetAK8_HFHM;   //!
   TBranch        *b_PFJetAK8_Neucons;   //!
   TBranch        *b_PFJetAK8_Chcons;   //!
   TBranch        *b_PFJetAK8_JER;   //!
   TBranch        *b_PFJetAK8_JERup;   //!
   TBranch        *b_PFJetAK8_JERdn;   //!
   TBranch        *b_PFJetAK8_msoftdrop;   //!
   TBranch        *b_PFJetAK8_tau1;   //!
   TBranch        *b_PFJetAK8_tau2;   //!
   TBranch        *b_PFJetAK8_tau3;   //!
   TBranch        *b_PFJetAK8_btag_DeepCSV;   //!
   TBranch        *b_PFJetAK8_DeepTag_DAK8MD_TvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_DAK8MD_WvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_DAK8MD_ZvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_DAK8MD_HvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_DAK8MD_bbvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_PNet_TvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_PNet_WvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_PNet_ZvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_PNetMD_XbbvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_PNetMD_XccvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_PNetMD_XqqvsQCD;   //!
   TBranch        *b_PFJetAK8_DeepTag_PNetMD_QCD;   //!
   TBranch        *b_PFJetAK8_sub1pt;   //!
   TBranch        *b_PFJetAK8_sub1eta;   //!
   TBranch        *b_PFJetAK8_sub1phi;   //!
   TBranch        *b_PFJetAK8_sub1mass;   //!
   TBranch        *b_PFJetAK8_sub1btag;   //!
   TBranch        *b_PFJetAK8_sub1JEC;   //!
   TBranch        *b_PFJetAK8_sub2pt;   //!
   TBranch        *b_PFJetAK8_sub2eta;   //!
   TBranch        *b_PFJetAK8_sub2phi;   //!
   TBranch        *b_PFJetAK8_sub2mass;   //!
   TBranch        *b_PFJetAK8_sub2btag;   //!
   TBranch        *b_PFJetAK8_sub2JEC;   //!
   TBranch        *b_PFJetAK8_jesup_AbsoluteStat;   //!
   TBranch        *b_PFJetAK8_jesup_AbsoluteScale;   //!
   TBranch        *b_PFJetAK8_jesup_AbsoluteMPFBias;   //!
   TBranch        *b_PFJetAK8_jesup_FlavorQCD;   //!
   TBranch        *b_PFJetAK8_jesup_Fragmentation;   //!
   TBranch        *b_PFJetAK8_jesup_PileUpDataMC;   //!
   TBranch        *b_PFJetAK8_jesup_PileUpPtBB;   //!
   TBranch        *b_PFJetAK8_jesup_PileUpPtEC1;   //!
   TBranch        *b_PFJetAK8_jesup_PileUpPtEC2;   //!
   TBranch        *b_PFJetAK8_jesup_PileUpPtRef;   //!
   TBranch        *b_PFJetAK8_jesup_RelativeFSR;   //!
   TBranch        *b_PFJetAK8_jesup_RelativeJEREC1;   //!
   TBranch        *b_PFJetAK8_jesup_RelativeJEREC2;   //!
   TBranch        *b_PFJetAK8_jesup_RelativePtBB;   //!
   TBranch        *b_PFJetAK8_jesup_RelativePtEC1;   //!
   TBranch        *b_PFJetAK8_jesup_RelativePtEC2;   //!
   TBranch        *b_PFJetAK8_jesup_RelativeBal;   //!
   TBranch        *b_PFJetAK8_jesup_RelativeSample;   //!
   TBranch        *b_PFJetAK8_jesup_RelativeStatEC;   //!
   TBranch        *b_PFJetAK8_jesup_RelativeStatFSR;   //!
   TBranch        *b_PFJetAK8_jesup_SinglePionECAL;   //!
   TBranch        *b_PFJetAK8_jesup_SinglePionHCAL;   //!
   TBranch        *b_PFJetAK8_jesup_TimePtEta;   //!
   TBranch        *b_PFJetAK8_jesup_Total;   //!
   TBranch        *b_PFJetAK8_jesdn_AbsoluteStat;   //!
   TBranch        *b_PFJetAK8_jesdn_AbsoluteScale;   //!
   TBranch        *b_PFJetAK8_jesdn_AbsoluteMPFBias;   //!
   TBranch        *b_PFJetAK8_jesdn_FlavorQCD;   //!
   TBranch        *b_PFJetAK8_jesdn_Fragmentation;   //!
   TBranch        *b_PFJetAK8_jesdn_PileUpDataMC;   //!
   TBranch        *b_PFJetAK8_jesdn_PileUpPtBB;   //!
   TBranch        *b_PFJetAK8_jesdn_PileUpPtEC1;   //!
   TBranch        *b_PFJetAK8_jesdn_PileUpPtEC2;   //!
   TBranch        *b_PFJetAK8_jesdn_PileUpPtRef;   //!
   TBranch        *b_PFJetAK8_jesdn_RelativeFSR;   //!
   TBranch        *b_PFJetAK8_jesdn_RelativeJEREC1;   //!
   TBranch        *b_PFJetAK8_jesdn_RelativeJEREC2;   //!
   TBranch        *b_PFJetAK8_jesdn_RelativePtBB;   //!
   TBranch        *b_PFJetAK8_jesdn_RelativePtEC1;   //!
   TBranch        *b_PFJetAK8_jesdn_RelativePtEC2;   //!
   TBranch        *b_PFJetAK8_jesdn_RelativeBal;   //!
   TBranch        *b_PFJetAK8_jesdn_RelativeSample;   //!
   TBranch        *b_PFJetAK8_jesdn_RelativeStatEC;   //!
   TBranch        *b_PFJetAK8_jesdn_RelativeStatFSR;   //!
   TBranch        *b_PFJetAK8_jesdn_SinglePionECAL;   //!
   TBranch        *b_PFJetAK8_jesdn_SinglePionHCAL;   //!
   TBranch        *b_PFJetAK8_jesdn_TimePtEta;   //!
   TBranch        *b_PFJetAK8_jesdn_Total;   //!
   TBranch        *b_nPFJetAK4;   //!
   TBranch        *b_PFJetAK4_jetID;   //!
   TBranch        *b_PFJetAK4_jetID_tightlepveto;   //!
   TBranch        *b_PFJetAK4_pt;   //!
   TBranch        *b_PFJetAK4_eta;   //!
   TBranch        *b_PFJetAK4_y;   //!
   TBranch        *b_PFJetAK4_phi;   //!
   TBranch        *b_PFJetAK4_mass;   //!
   TBranch        *b_PFJetAK4_JEC;   //!
   TBranch        *b_PFJetAK4_btag_DeepCSV;   //!
   TBranch        *b_PFJetAK4_btag_DeepFlav;   //!
   TBranch        *b_PFJetAK4_JER;   //!
   TBranch        *b_PFJetAK4_JERup;   //!
   TBranch        *b_PFJetAK4_JERdn;   //!
   TBranch        *b_PFJetAK4_hadronflav;   //!
   TBranch        *b_PFJetAK4_partonflav;   //!
   TBranch        *b_PFJetAK4_qgl;   //!
   TBranch        *b_PFJetAK4_PUID;   //!
   TBranch        *b_PFJetAK4_jesup_AbsoluteStat;   //!
   TBranch        *b_PFJetAK4_jesup_AbsoluteScale;   //!
   TBranch        *b_PFJetAK4_jesup_AbsoluteMPFBias;   //!
   TBranch        *b_PFJetAK4_jesup_FlavorQCD;   //!
   TBranch        *b_PFJetAK4_jesup_Fragmentation;   //!
   TBranch        *b_PFJetAK4_jesup_PileUpDataMC;   //!
   TBranch        *b_PFJetAK4_jesup_PileUpPtBB;   //!
   TBranch        *b_PFJetAK4_jesup_PileUpPtEC1;   //!
   TBranch        *b_PFJetAK4_jesup_PileUpPtEC2;   //!
   TBranch        *b_PFJetAK4_jesup_PileUpPtRef;   //!
   TBranch        *b_PFJetAK4_jesup_RelativeFSR;   //!
   TBranch        *b_PFJetAK4_jesup_RelativeJEREC1;   //!
   TBranch        *b_PFJetAK4_jesup_RelativeJEREC2;   //!
   TBranch        *b_PFJetAK4_jesup_RelativePtBB;   //!
   TBranch        *b_PFJetAK4_jesup_RelativePtEC1;   //!
   TBranch        *b_PFJetAK4_jesup_RelativePtEC2;   //!
   TBranch        *b_PFJetAK4_jesup_RelativeBal;   //!
   TBranch        *b_PFJetAK4_jesup_RelativeSample;   //!
   TBranch        *b_PFJetAK4_jesup_RelativeStatEC;   //!
   TBranch        *b_PFJetAK4_jesup_RelativeStatFSR;   //!
   TBranch        *b_PFJetAK4_jesup_SinglePionECAL;   //!
   TBranch        *b_PFJetAK4_jesup_SinglePionHCAL;   //!
   TBranch        *b_PFJetAK4_jesup_TimePtEta;   //!
   TBranch        *b_PFJetAK4_jesup_Total;   //!
   TBranch        *b_PFJetAK4_jesdn_AbsoluteStat;   //!
   TBranch        *b_PFJetAK4_jesdn_AbsoluteScale;   //!
   TBranch        *b_PFJetAK4_jesdn_AbsoluteMPFBias;   //!
   TBranch        *b_PFJetAK4_jesdn_FlavorQCD;   //!
   TBranch        *b_PFJetAK4_jesdn_Fragmentation;   //!
   TBranch        *b_PFJetAK4_jesdn_PileUpDataMC;   //!
   TBranch        *b_PFJetAK4_jesdn_PileUpPtBB;   //!
   TBranch        *b_PFJetAK4_jesdn_PileUpPtEC1;   //!
   TBranch        *b_PFJetAK4_jesdn_PileUpPtEC2;   //!
   TBranch        *b_PFJetAK4_jesdn_PileUpPtRef;   //!
   TBranch        *b_PFJetAK4_jesdn_RelativeFSR;   //!
   TBranch        *b_PFJetAK4_jesdn_RelativeJEREC1;   //!
   TBranch        *b_PFJetAK4_jesdn_RelativeJEREC2;   //!
   TBranch        *b_PFJetAK4_jesdn_RelativePtBB;   //!
   TBranch        *b_PFJetAK4_jesdn_RelativePtEC1;   //!
   TBranch        *b_PFJetAK4_jesdn_RelativePtEC2;   //!
   TBranch        *b_PFJetAK4_jesdn_RelativeBal;   //!
   TBranch        *b_PFJetAK4_jesdn_RelativeSample;   //!
   TBranch        *b_PFJetAK4_jesdn_RelativeStatEC;   //!
   TBranch        *b_PFJetAK4_jesdn_RelativeStatFSR;   //!
   TBranch        *b_PFJetAK4_jesdn_SinglePionECAL;   //!
   TBranch        *b_PFJetAK4_jesdn_SinglePionHCAL;   //!
   TBranch        *b_PFJetAK4_jesdn_TimePtEta;   //!
   TBranch        *b_PFJetAK4_jesdn_Total;   //!
   TBranch        *b_PFJetAK4_btag_DeepCSV_SF;   //!
   TBranch        *b_PFJetAK4_btag_DeepCSV_SF_up;   //!
   TBranch        *b_PFJetAK4_btag_DeepCSV_SF_dn;   //!
   TBranch        *b_PFJetAK4_btag_DeepFlav_SF;   //!
   TBranch        *b_PFJetAK4_btag_DeepFlav_SF_up;   //!
   TBranch        *b_PFJetAK4_btag_DeepFlav_SF_dn;   //!
   TBranch        *b_Muon_Muon_TightID;
   TBranch        *b_Muon_minisoch;
   TBranch        *b_Muon_minisonh;
   TBranch        *b_Muon_minisoph;
   TBranch        *b_nMuon;   //!
   TBranch        *b_Muon_isPF;   //!
   TBranch        *b_Muon_isGL;   //!
   TBranch        *b_Muon_isTRK;   //!
   TBranch        *b_Muon_isLoose;   //!
   TBranch        *b_Muon_isGoodGL;   //!
   TBranch        *b_Muon_isMed;   //!
   TBranch        *b_Muon_isMedPr;   //!
   TBranch        *b_Muon_isTight;   //!
   TBranch        *b_Muon_isHighPt;   //!
   TBranch        *b_Muon_isHighPttrk;   //!
   TBranch        *b_Muon_MediumID;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_p;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_minchiso;   //!
   TBranch        *b_Muon_minnhiso;   //!
   TBranch        *b_Muon_minphiso;   //!
   TBranch        *b_Muon_minisoall;   //!
   TBranch        *b_Muon_dxy;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_dxyErr;   //!
   TBranch        *b_Muon_ip3d;   //!
   TBranch        *b_Muon_ptErr;   //!
   TBranch        *b_Muon_chi;   //!
   TBranch        *b_Muon_ndf;   //!
   TBranch        *b_Muon_ecal;   //!
   TBranch        *b_Muon_hcal;   //!
   TBranch        *b_Muon_pfiso;   //!
   TBranch        *b_Muon_posmatch;   //!
   TBranch        *b_Muon_trkink;   //!
   TBranch        *b_Muon_segcom;   //!
   TBranch        *b_Muon_hit;   //!
   TBranch        *b_Muon_pixhit;   //!
   TBranch        *b_Muon_mst;   //!
   TBranch        *b_Muon_trklay;   //!
   TBranch        *b_Muon_valfrac;   //!
   TBranch        *b_Muon_dxy_sv;   //!
   TBranch        *b_nElectron;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_p;   //!
   TBranch        *b_Electron_e;   //!
   TBranch        *b_Electron_e_ECAL;   //!
   TBranch        *b_Electron_mvaid_Fallv2WP90;   //!
   TBranch        *b_Electron_mvaid_Fallv2WP90_noIso;   //!
   TBranch        *b_Electron_mvaid_Fallv2WP80;   //!
   TBranch        *b_Electron_mvaid_Fallv2WP80_noIso;   //!
   TBranch        *b_Electron_dxy;   //!
   TBranch        *b_Electron_dxyErr;   //!
   TBranch        *b_Electron_dz;   //!
   TBranch        *b_Electron_dzErr;   //!
   TBranch        *b_Electron_ip3d;   //!
   TBranch        *b_Electron_dxy_sv;   //!
   TBranch        *b_Electron_hovere;   //!
   TBranch        *b_Electron_chi;   //!
   TBranch        *b_Electron_ndf;   //!
   TBranch        *b_Electron_eoverp;   //!
   TBranch        *b_Electron_ietaieta;   //!
   TBranch        *b_Electron_misshits;   //!
   TBranch        *b_Electron_pfiso_drcor;   //!
   TBranch        *b_Electron_pfiso_eacor;   //!
   TBranch        *b_Electron_pfiso04_eacor;   //!
   TBranch        *b_Electron_eccalTrkEnergyPostCorr;   //!
   TBranch        *b_Electron_energyScaleValue;   //!
   TBranch        *b_Electron_energyScaleUp;   //!
   TBranch        *b_Electron_energyScaleDown;   //!
   TBranch        *b_Electron_energySigmaValue;   //!
   TBranch        *b_Electron_energySigmaUp;   //!
   TBranch        *b_Electron_energySigmaDown;   //!
   TBranch        *b_Electron_supcl_eta;   //!
   TBranch        *b_Electron_supcl_phi;   //!
   TBranch        *b_Electron_supcl_e;   //!
   TBranch        *b_Electron_supcl_rawE;   //!
   TBranch        *b_Electron_sigmaieta;   //!
   TBranch        *b_Electron_sigmaiphi;   //!
   TBranch        *b_Electron_r9full;   //!
   TBranch        *b_Electron_hcaloverecal;   //!
   TBranch        *b_Electron_hitsmiss;   //!
   TBranch        *b_Electron_ecloverpout;   //!
   TBranch        *b_Electron_convVeto;   //!
   TBranch        *b_Electron_pfisolsumphet;   //!
   TBranch        *b_Electron_pfisolsumchhadpt;   //!
   TBranch        *b_Electron_pfsiolsumneuhadet;   //!
   TBranch        *b_Electron_minchiso;   //!
   TBranch        *b_Electron_minnhiso;   //!
   TBranch        *b_Electron_minphiso;   //!
   TBranch        *b_Electron_minisoall;   //!
   TBranch        *b_Electron_minisoch;
   TBranch        *b_Electron_minisonh;
   TBranch        *b_Electron_minisoph;
   TBranch        *b_Generator_weight;
   TBranch        *b_Generator_qscale;
   TBranch        *b_LHE_weight;
   TBranch        *b_nPhoton;   //!
   TBranch        *b_Photon_e;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_mvaid_Fall17V2_raw;   //!
   TBranch        *b_Photon_mvaid_Fall17V2_WP90;   //!
   TBranch        *b_Photon_mvaid_Fall17V2_WP80;   //!
   TBranch        *b_Photon_mvaid_Spring16V1_WP90;   //!
   TBranch        *b_Photon_mvaid_Spring16V1_WP80;   //!
   TBranch        *b_Photon_e1by9;   //!
   TBranch        *b_Photon_e9by25;   //!
   TBranch        *b_Photon_trkiso;   //!
   TBranch        *b_Photon_emiso;   //!
   TBranch        *b_Photon_hadiso;   //!
   TBranch        *b_Photon_chhadiso;   //!
   TBranch        *b_Photon_neuhadiso;   //!
   TBranch        *b_Photon_phoiso;   //!
   TBranch        *b_Photon_PUiso;   //!
   TBranch        *b_Photon_hadbyem;   //!
   TBranch        *b_Photon_ietaieta;   //!
   TBranch        *b_nTau;   //!
   TBranch        *b_Tau_isPF;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_e;   //!
   TBranch        *b_Tau_charge;   //!
   TBranch        *b_Tau_dxy;   //!
   TBranch        *b_Tau_leadtrkdxy;   //!
   TBranch        *b_Tau_leadtrkdz;   //!
   TBranch        *b_Tau_leadtrkpt;   //!
   TBranch        *b_Tau_leadtrketa;   //!
   TBranch        *b_Tau_leadtrkphi;   //!
   TBranch        *b_Tau_decayMode;   //!
   TBranch        *b_Tau_decayModeinding;   //!
   TBranch        *b_Tau_decayModeindingNewDMs;   //!
   TBranch        *b_Tau_eiso2018_raw;   //!
   TBranch        *b_Tau_eiso2018;   //!
   TBranch        *b_Tau_jetiso_deeptau2017v2p1_raw;   //!
   TBranch        *b_Tau_jetiso_deeptau2017v2p1;   //!
   TBranch        *b_Tau_eiso_deeptau2017v2p1_raw;   //!
   TBranch        *b_Tau_eiso_deeptau2017v2p1;   //!
   TBranch        *b_Tau_muiso_deeptau2017v2p1_raw;   //!
   TBranch        *b_Tau_muiso_deeptau2017v2p1;   //!
   TBranch        *b_Tau_rawiso;   //!
   TBranch        *b_Tau_rawisodR03;   //!
   TBranch        *b_Tau_puCorr;   //!
   TBranch        *b_event_weight;   //!
   TBranch        *b_qscale;   //!
   TBranch        *b_Generator_x1;   //!
   TBranch        *b_Generator_x2;   //!
   TBranch        *b_Generator_xpdf1;   //!
   TBranch        *b_Generator_xpdf2;   //!
   TBranch        *b_Generator_id1;   //!
   TBranch        *b_Generator_id2;   //!
   TBranch        *b_Generator_scalePDF;   //!
   TBranch        *b_npu_vert;   //!
   TBranch        *b_npu_vert_true;   //!
   TBranch        *b_genmiset;   //!
   TBranch        *b_genmisphi;   //!
   TBranch        *b_nGenJetAK8;   //!
   TBranch        *b_GenJetAK8_pt;   //!
   TBranch        *b_GenJetAK8_eta;   //!
   TBranch        *b_GenJetAK8_phi;   //!
   TBranch        *b_GenJetAK8_mass;   //!
   TBranch        *b_GenJetAK8_sdmass;   //!
   TBranch        *b_GenJetAK8_hadronflav;   //!
   TBranch        *b_GenJetAK8_partonflav;   //!
   TBranch        *b_nGenJetAK4;   //!
   TBranch        *b_GenJetAK4_pt;   //!
   TBranch        *b_GenJetAK4_eta;   //!
   TBranch        *b_GenJetAK4_phi;   //!
   TBranch        *b_GenJetAK4_mass;   //!
   TBranch        *b_GenJetAK4_hadronflav;   //!
   TBranch        *b_GenJetAK4_partonflav;   //!
   TBranch        *b_nGenPart;   //!
   TBranch        *b_GenPart_pt;   //!
   TBranch        *b_GenPart_eta;   //!
   TBranch        *b_GenPart_phi;   //!
   TBranch        *b_GenPart_m;   //!
   TBranch        *b_GenPart_status;   //!
   TBranch        *b_GenPart_pdgId;   //!
   TBranch        *b_GenPart_mompdgId;   //!
   TBranch        *b_GenPart_grmompdgId;   //!
   TBranch        *b_GenPart_daugno;   //!
   TBranch        *b_GenPart_fromhard;   //!
   TBranch        *b_GenPart_fromhardbFSR;   //!
   TBranch        *b_GenPart_isPromptFinalState;   //!
   TBranch        *b_GenPart_isLastCopyBeforeFSR;   //!
   TBranch        *b_nLHEPart;   //!
   TBranch        *b_LHEPart_pdg;   //!
   TBranch        *b_LHEPart_pt;   //!
   TBranch        *b_LHEPart_eta;   //!
   TBranch        *b_LHEPart_phi;   //!
   TBranch        *b_LHEPart_m;   //!
   TBranch        *b_event_weight_LHE;   //!
   TBranch        *b_nLHEScaleWeights;   //!
   TBranch        *b_LHEScaleWeights;   //!
   TBranch        *b_nLHEPDFWeights;   //!
   TBranch        *b_LHEPDFWeights;   //!
   TBranch        *b_nLHEAlpsWeights;   //!
   TBranch        *b_LHEAlpsWeights;   //!
   TBranch        *b_nLHEPSWeights;   //!
   TBranch        *b_LHEPSWeights;   //!
   
   double pu_rat18[100] =    {15.0025,9.25619,7.25497,5.06682,3.81133,3.00109,2.47446,2.15125,1.91921,1.76245,1.65799,1.5848,1.53433,1.49373,1.46908,1.44313,1.42046,1.40154,1.37988,1.35679,1.3361,1.31139,1.28604,1.26065,1.23868,1.21481,1.19261,1.17143,1.15268,1.13366,1.11664,1.10102,1.08379,1.06829,1.05251,1.03584,1.01745,1.00004,0.980328,0.958154,0.935215,0.910953,0.885202,0.856835,0.827565,0.797445,0.765605,0.733814,0.701484,0.668361,0.634682,0.60224,0.570207,0.537462,0.505992,0.475353,0.445231,0.417069,0.390368,0.363391,0.339587,0.316478,0.293616,0.272703,0.253076,0.23398,0.216635,0.200789,0.185823,0.171907,0.159236,0.148382,0.13732,0.127015,0.11736,0.108589,0.0994979,0.0923745,0.0861297,0.0803695,0.0741731,0.0689201,0.0636846,0.0574831,0.053418,0.0491854,0.0468157,0.0431244,0.0407697,0.0352331,0.0327207,0.0311201,0.0284833,0.0257933,0.0234082,0.0200422,0.0199332,0.0189272,0.020733,0.0166171};
   double pu_rat18_up[100] = {0,11.3701,49.1593,16.3978,10.4484,7.79227,5.70396,4.15872,3.02768,2.28549,1.82582,1.52983,1.3595,1.2554,1.19605,1.1684,1.16115,1.17185,1.18964,1.20936,1.22873,1.23491,1.23159,1.21107,1.18259,1.14644,1.11133,1.08136,1.05384,1.03331,1.01987,1.01367,1.01107,1.01298,1.01865,1.02593,1.03512,1.0447,1.05099,1.0554,1.05447,1.04466,1.02824,1.00332,0.965566,0.923431,0.871249,0.814665,0.752156,0.689408,0.624858,0.564,0.505617,0.452167,0.402,0.359344,0.321227,0.285921,0.258403,0.233682,0.210464,0.192413,0.174424,0.159861,0.146181,0.131623,0.119227,0.10899,0.0963316,0.086803,0.0773651,0.0712667,0.0629173,0.0552031,0.0481823,0.0455058,0.0376989,0.0339163,0.0298286,0.0264131,0.0255965,0.0179475,0.0169746,0.0136435,0.0117583,0.00988318,0.00674005,0.00661599,0.00316237,0.00149674,0.0010104,0.00106782,0.000384941,0.000591271,0.000423128,0.000165822,7.60044e-05,4.96232e-05,7.51979e-05,1.05862e-05};
   double pu_rat18_dn[100] = {0,15.0557,67.8751,22.3278,14.1211,10.4821,7.88069,5.86513,4.31762,3.35551,2.78627,2.40097,2.16428,2.00485,1.9056,1.85092,1.82051,1.80608,1.78719,1.75544,1.71117,1.64481,1.57234,1.49261,1.42092,1.35612,1.3043,1.26517,1.23118,1.20443,1.18302,1.16596,1.14834,1.13047,1.11055,1.08517,1.05388,1.01479,0.96502,0.907499,0.841466,0.767187,0.68971,0.610695,0.530471,0.45611,0.385995,0.32355,0.268127,0.221267,0.181416,0.149012,0.122387,0.100955,0.0832931,0.0694147,0.0579993,0.0482614,0.0406839,0.0341693,0.0284128,0.0238208,0.0196651,0.0163071,0.0134164,0.0108213,0.00875349,0.00713274,0.00561523,0.00450669,0.00357902,0.00293888,0.00231295,0.00180802,0.00140385,0.00117654,0.000861839,0.000682485,0.000525487,0.000404909,0.00033922,0.000204219,0.000164688,0.000112084,8.12391e-05,5.70485e-05,3.2298e-05,2.61592e-05,1.02574e-05,3.96059e-06,2.16985e-06,1.85204e-06,5.36884e-07,6.60936e-07,3.78607e-07,1.19189e-07,4.4536e-08,2.4673e-08,3.47283e-08,5.35281e-09};

   float muon_pt_cut = 25;
   float electron_pt_cut = 25;
   float lepton_pt_cut = 30;
   float AK4GenJet_pt_cut = 8.0;
   float AK4jet_pt_cut = 30;
   float AK8jet_pt_cut = 200;
   float AK8GenJet_pt_cut = 50;
   float absetacut = 2.5;
   
   float deep_btag_cut = 0.2783; 
   //for UL18 => 0.0490: loose, 0.2783: medium, 0.7100: tight 

   //DeepTag_PNetMD_XbbvsQCD
   float PNetbb_cut_T = 0.98;
   float PNetbb_cut_M = 0.94;
   float PNetbb_cut_L = 0.90;
   //DeepTag_DAK8MD_WvsQCD cut values
   float DAK8W_cut_T = 0.806;
   float DAK8W_cut_M = 0.704;
   float DAK8W_cut_L = 0.479;
   //DeepTag_PNetMD_WvsQCD cut values
   float PNetW_cut_T = 0.90;
   float PNetW_cut_M = 0.82;
   float PNetW_cut_L = 0.59;
   // Deep Ak4 Flv
   float DAK4_T = 0.71;
   float DAK4_M = 0.2783;
   float DAK4_L = 0.0490;

   bool isMC;
   bool isFastSIM;
   
   bool channel_SL = true;
   bool channel_DL = false;
   
   TRandom3* gxRandom;
   
   TTree *Tout ;
   TTree *Tout_presel; 
   int nleptons, nfatjets;
   
   float l_pt, l_eta, l_phi, l_mass;
   int l_pdgId;
   float l_minisoch, l_minisonh, l_minisoph, l_minisoall; 
   int l_genindex;
   
   float MET_pt, MET_phi, MET_sig, MET_sumEt;
   float MET_pt_JESup, MET_pt_JESdn, MET_pt_JERup, MET_pt_JERdn, MET_pt_UnclusEup, MET_pt_UnclusEdn;
   float MET_phi_JESup, MET_phi_JESdn, MET_phi_JERup, MET_phi_JERdn, MET_phi_UnclusEup, MET_phi_UnclusEdn;
   
   float Y_pt, Y_y, Y_eta, Y_phi, Y_mass;
   float Y_msoftdrop, Y_tau21, Y_tau32;
   float Y_DeepTag_DAK8MD_TvsQCD, Y_DeepTag_DAK8MD_WvsQCD, Y_DeepTag_DAK8MD_ZvsQCD, Y_DeepTag_DAK8MD_HvsQCD, Y_DeepTag_DAK8MD_bbvsQCD; 
   float Y_DeepTag_PNet_TvsQCD, Y_DeepTag_PNet_WvsQCD, Y_DeepTag_PNet_ZvsQCD, Y_DeepTag_PNetMD_XbbvsQCD, Y_DeepTag_PNetMD_XccvsQCD, Y_DeepTag_PNetMD_XqqvsQCD, Y_DeepTag_PNetMD_QCD, Y_DeepTag_PNetMD_WvsQCD; 
   float Y_PN_bb;
   bool Y_label_Top_bq, Y_label_Top_bc, Y_label_Top_bcq, Y_label_Top_bqq, Y_label_W_qq, Y_label_W_cq;
   float Y_sub1_pt, Y_sub1_eta, Y_sub1_phi, Y_sub1_mass, Y_sub1_btag;
   float Y_sub2_pt, Y_sub2_eta, Y_sub2_phi, Y_sub2_mass, Y_sub2_btag;
   int Y_genindex, Y_genbindex[2];
   float Y_JESup, Y_JESdn, Y_JERup, Y_JERdn;
   
   float W_pt_opt1, W_y_opt1, W_eta_opt1, W_phi_opt1, W_mass_opt1; 
   float W_msoftdrop_opt1, W_tau21_opt1, W_tau32_opt1;
   float W_DeepTag_DAK8MD_TvsQCD_opt1, W_DeepTag_DAK8MD_WvsQCD_opt1, W_DeepTag_DAK8MD_ZvsQCD_opt1, W_DeepTag_DAK8MD_HvsQCD_opt1, W_DeepTag_DAK8MD_bbvsQCD_opt1; 
   float W_DeepTag_PNet_TvsQCD_opt1, W_DeepTag_PNet_WvsQCD_opt1, W_DeepTag_PNet_ZvsQCD_opt1, W_DeepTag_PNetMD_XbbvsQCD_opt1, W_DeepTag_PNetMD_XccvsQCD_opt1, W_DeepTag_PNetMD_XqqvsQCD_opt1, W_DeepTag_PNetMD_QCD_opt1, W_DeepTag_PNetMD_WvsQCD_opt1; 
   float W_DAK8_W_opt1, W_PN_W_opt1;
   bool W_label_W_qq_opt1, W_label_W_cq_opt1;
   float W_sub1_pt_opt1, W_sub1_eta_opt1, W_sub1_phi_opt1, W_sub1_mass_opt1, W_sub1_btag_opt1;
   float W_sub2_pt_opt1, W_sub2_eta_opt1, W_sub2_phi_opt1, W_sub2_mass_opt1, W_sub2_btag_opt1;
   int W_genindex_opt1;
   float W_JESup_opt1, W_JESdn_opt1, W_JERup_opt1, W_JERdn_opt1;
   
   float W_pt_opt2, W_y_opt2, W_eta_opt2, W_phi_opt2, W_mass_opt2;
   float W_msoftdrop_opt2, W_tau21_opt2, W_tau32_opt2;
   float W_DeepTag_DAK8MD_TvsQCD_opt2, W_DeepTag_DAK8MD_WvsQCD_opt2, W_DeepTag_DAK8MD_ZvsQCD_opt2, W_DeepTag_DAK8MD_HvsQCD_opt2, W_DeepTag_DAK8MD_bbvsQCD_opt2; 
   float W_DeepTag_PNet_TvsQCD_opt2, W_DeepTag_PNet_WvsQCD_opt2, W_DeepTag_PNet_ZvsQCD_opt2, W_DeepTag_PNetMD_XbbvsQCD_opt2, W_DeepTag_PNetMD_XccvsQCD_opt2, W_DeepTag_PNetMD_XqqvsQCD_opt2, W_DeepTag_PNetMD_QCD_opt2, W_DeepTag_PNetMD_WvsQCD_opt2; 
   float W_DAK8_W_opt2, W_PN_W_opt2;
   bool W_label_W_qq_opt2, W_label_W_cq_opt2;
   float W_sub1_pt_opt2, W_sub1_eta_opt2, W_sub1_phi_opt2, W_sub1_mass_opt2, W_sub1_btag_opt2;
   float W_sub2_pt_opt2, W_sub2_eta_opt2, W_sub2_phi_opt2, W_sub2_mass_opt2, W_sub2_btag_opt2;
   int W_genindex_opt2;
   float W_JESup_opt2, W_JESdn_opt2, W_JERup_opt2, W_JERdn_opt2;
   
   float H_pt_opt1, H_y_opt1, H_eta_opt1, H_phi_opt1, H_mass_opt1;
   int H_genindex_opt1;
   float H_JESup_opt1, H_JESdn_opt1, H_JERup_opt1, H_JERdn_opt1;
   
   float H_pt_opt2, H_y_opt2, H_eta_opt2, H_phi_opt2, H_mass_opt2;
   int H_genindex_opt2;
   float H_JESup_opt2, H_JESdn_opt2, H_JERup_opt2, H_JERdn_opt2;

   float X_mass_opt1, X_mass_opt2; 
   
   float dR_lW_opt1, dphi_lW_opt1, dy_lW_opt1;
   float dR_lW_opt2, dphi_lW_opt2, dy_lW_opt2;

   float dR_lY, dphi_lY, dy_lY;
   int nbjets_other, nbjets_outY;
   
   bool Flag_Y_bb_pass_T, Flag_Y_bb_pass_M, Flag_Y_bb_pass_L, Flag_H_W_pass_T_opt1, Flag_H_W_pass_M_opt1, Flag_H_W_pass_L_opt1, Flag_H_m_pass_opt1, Flag_dR_lW_pass_opt1, Flag_MET_pass;
   bool Flag_H_W_pass_T_opt2, Flag_H_W_pass_M_opt2, Flag_H_W_pass_L_opt2, Flag_H_m_pass_opt2, Flag_dR_lW_pass_opt2;
   bool Reg_SR_opt1, Reg_Wj_CR_opt1;
   bool Reg_SR_opt2, Reg_Wj_CR_opt2;
   
   int _s_nPFJetAK8; 
   float _s_PFJetAK8_pt[njetmxAK8], _s_PFJetAK8_eta[njetmxAK8], _s_PFJetAK8_phi[njetmxAK8], _s_PFJetAK8_mass[njetmxAK8];
   float _s_PFJetAK8_msoftdrop[njetmxAK8], _s_PFJetAK8_tau21[njetmxAK8], _s_PFJetAK8_tau32[njetmxAK8];
   bool _s_PFJetAK8_jetID[njetmxAK8], _s_PFJetAK8_jetID_tightlepveto[njetmxAK8];
   float _s_PFJetAK8_DeepTag_PNetMD_XbbvsQCD[njetmxAK8],  _s_PFJetAK8_DeepTag_PNetMD_WvsQCD[njetmxAK8], _s_PFJetAK8_DeepTag_PNet_TvsQCD[njetmxAK8], _s_PFJetAK8_DeepTag_PNet_WvsQCD[njetmxAK8];
   float _s_PFJetAK8_DeepTag_DAK8MD_TvsQCD[njetmxAK8], _s_PFJetAK8_DeepTag_DAK8MD_WvsQCD[njetmxAK8], _s_PFJetAK8_DeepTag_DAK8MD_bbvsQCD[njetmxAK8]; 
   float _s_PFJetAK8_JESup[njetmxAK8], _s_PFJetAK8_JESdn[njetmxAK8], _s_PFJetAK8_JERup[njetmxAK8], _s_PFJetAK8_JERdn[njetmxAK8];
   int _s_PFJetAK8_Y_index, _s_PFJetAK8_W_index_opt1, _s_PFJetAK8_W_index_opt2;
   
   int _s_nBJetAK4;
   float _s_BJetAK4_pt[njetmx], _s_BJetAK4_eta[njetmx], _s_BJetAK4_phi[njetmx], _s_BJetAK4_mass[njetmx];
   float _s_BJetAK4_btag_DeepFlav[njetmx], _s_BJetAK4_btag_DeepCSV[njetmx];
   int _s_BJetAK4_hadronflav[njetmx], _s_BJetAK4_partonflav[njetmx];
   float _s_BJetAK4_qgl[njetmx], _s_BJetAK4_PUID[njetmx]; 
   float _s_BJetAK4_JESup[njetmx], _s_BJetAK4_JESdn[njetmx], _s_BJetAK4_JERup[njetmx], _s_BJetAK4_JERdn[njetmx];
   float _s_BJetAK4_btag_DeepFlav_SF[njetmx], _s_BJetAK4_btag_DeepFlav_SF_up[njetmx], _s_BJetAK4_btag_DeepFlav_SF_dn[njetmx];

   int _s_nJetAK4;
   float _s_JetAK4_pt[njetmx], _s_JetAK4_eta[njetmx], _s_JetAK4_phi[njetmx], _s_JetAK4_mass[njetmx];
   float _s_JetAK4_btag_DeepFlav[njetmx], _s_JetAK4_btag_DeepCSV[njetmx];
   int _s_JetAK4_hadronflav[njetmx], _s_JetAK4_partonflav[njetmx];
   float _s_JetAK4_qgl[njetmx], _s_JetAK4_PUID[njetmx];
   float _s_JetAK4_JESup[njetmx], _s_JetAK4_JESdn[njetmx], _s_JetAK4_JERup[njetmx], _s_JetAK4_JERdn[njetmx];
   float _s_JetAK4_btag_DeepFlav_SF[njetmx], _s_JetAK4_btag_DeepFlav_SF_up[njetmx], _s_JetAK4_btag_DeepFlav_SF_dn[njetmx];

   int nGenLep;
   float GenLep_pt[njetmx], GenLep_eta[njetmx], GenLep_phi[njetmx], GenLep_mass[njetmx];
   int GenLep_pdgId[njetmx], GenLep_mompdgId[njetmx], GenLep_grmompdgId[njetmx];
   
   int nGenNu;
   float GenNu_pt[njetmx], GenNu_eta[njetmx], GenNu_phi[njetmx], GenNu_mass[njetmx];
   int GenNu_pdgId[njetmx], GenNu_mompdgId[njetmx], GenNu_grmompdgId[njetmx];
   
   int nGenBPart;
   float GenBPart_pt[njetmx], GenBPart_eta[njetmx], GenBPart_phi[njetmx], GenBPart_mass[njetmx];
   int GenBPart_pdgId[njetmx], GenBPart_mompdgId[njetmx], GenBPart_grmompdgId[njetmx];
   
   int nGenV;
   float GenV_pt[njetmx], GenV_eta[njetmx], GenV_phi[njetmx], GenV_mass[njetmx];
   int GenV_pdgId[njetmx], GenV_mompdgId[njetmx], GenV_grmompdgId[njetmx];
   
   double weight;
   
   float puWeight, puWeightup, puWeightdown;
   float leptonsf_weight, leptonsf_weight_up, leptonsf_weight_dn, leptonsf_weight_stat, leptonsf_weight_syst;
   
   BTagCalibration calib_deepcsv, calib_deepflav;
   BTagCalibrationReader reader_deepcsv, reader_deepflav;
   
   int year = 2018; 
   string muon_id_name = "Tight";
   string electron_id_name = "wp90noiso";
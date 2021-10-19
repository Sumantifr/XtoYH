// -*- C++ -*-
//
// Package:    Run2_2016/TopplusB
// Class:      TopplusB
// 
/**\class NTuplizer_XYH NTuplizer_XYH.cc 
   
   Description: [one line class summary]
   
   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Suman Chatterjee
//         Created:  Fri, 1 Oct 2021 16:22:44 GMT
//

//Twikis used:
/*
Prefiring weights: https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1PrefiringWeightRecipe#2018_UL 
Electron MVA ID: https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
Photon MVA ID: https://twiki.cern.ch/twiki/bin/view/CMS/MultivariatePhotonIdentificationRun2 
Main EGamma: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2#MVA_based_electron_Identificatio
JEC: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECUncertaintySources
DeepAKX: https://twiki.cern.ch/twiki/bin/viewauth/CMS/DeepAKXTagging
*/
// system include files-
#include <memory>
// user include files
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TAxis.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "TRandom.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "GeneratorInterface/Pythia8Interface/plugins/ReweightUserHooks.h"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include <string>
#include <iostream>
#include <fstream>
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "FWCore/Utilities/interface/typelookup.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include  "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "RecoEgamma/ElectronIdentification/interface/ElectronMVAEstimatorRun2.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtTrigReportEntry.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtTrigReport.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include "FWCore/Utilities/interface/typelookup.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include <fastjet/GhostedAreaSpec.hh>
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"

using namespace std;
using namespace edm;
using namespace reco;  
using namespace CLHEP;
using namespace trigger;
using namespace math;
using namespace fastjet;
using namespace fastjet::contrib;

const float mu_mass = 0.105658;
const float el_mass = 0.000511;
const float pival = acos(-1.);

static const int nsrc = 24;
const char* jecsrcnames[nsrc] = {
	 "AbsoluteStat", "AbsoluteScale","AbsoluteMPFBias", 
	 "FlavorQCD", "Fragmentation", 
	 "PileUpDataMC",  "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", //"PileUpPtHF",
	 "PileUpPtRef",
	 "RelativeFSR", "RelativeJEREC1", "RelativeJEREC2", //"RelativeJERHF",
	 "RelativePtBB", "RelativePtEC1", "RelativePtEC2", //"RelativePtHF", 
	 "RelativeBal", "RelativeSample", "RelativeStatEC", "RelativeStatFSR", //"RelativeStatHF", 
	 "SinglePionECAL", "SinglePionHCAL","TimePtEta",
	 "Total"
	};
const int njecmcmx = 2*nsrc + 1 ;

struct triggervar{
  TLorentzVector  trg4v;
  bool		  both;
  bool            level1;
  bool            highl;
  int             ihlt;
  int             prescl;
  int             pdgId;
};

int getbinid(double val, int nbmx, double* array) {
  if (val<array[0]) return -2;
  for (int ix=0; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return -3;
}

double theta_to_eta(double theta) { return -log(tan(theta/2.)); }

double PhiInRange(const double& phi) {
  double phiout = phi;
  if( phiout > 2*M_PI || phiout < -2*M_PI) {
    phiout = fmod( phiout, 2*M_PI);
  }
  if (phiout <= -M_PI) phiout += 2*M_PI;
  else if (phiout >  M_PI) phiout -= 2*M_PI;
  return phiout;
}

double delta2R(double eta1, double phi1, double eta2, double phi2) {
  return sqrt(pow(eta1 - eta2,2) +pow(PhiInRange(phi1 - phi2),2));
}

double diff_func(double f1, double f2){
  double ff = pow(f1-f2,2)*1./pow(f1+f2,2);
  return ff;
}


TLorentzVector productX(TLorentzVector X, TLorentzVector Y, float pro1, float pro2)
{
  float b1, b2, b3;
  float c1, c2, c3;
  
  b1 = X.Px();
  b2 = X.Py();
  b3 = X.Pz();
  
  c1 = Y.Px();
  c2 = Y.Py();
  c3 = Y.Pz();
  
  float d1, d2, e1, e2, X1, X2;
  
  X1 = pro1;
  X2 = pro2;
  
  d1 = (c2*X1 - b2*X2)*1./(b1*c2 - b2*c1);
  d2 = (c1*X1 - b1*X2)*1./(b2*c1 - b1*c2);
  e1 = (b2*c3 - b3*c2)*1./(b1*c2 - b2*c1);
  e2 = (b1*c3 - b3*c1)*1./(b2*c1 - b1*c2);
  
  float A, B, C;
  A = (e1*e1 + e2*e2+ 1);
  B = 2*(d1*e1 + d2*e2);
  C = d1*d1 + d2*d2 - 1;
  
  float sol;
  
  if((pow(B,2) - (4*A*C)) < 0){
    sol = -1*B/(2*A);
    
    float A1, A2, A3;
    A3 = sol;
    A1 = d1 + e1*A3;
    A2 = d2 + e2*A3;
    
    TLorentzVector vec4;
    vec4.SetPxPyPzE(A1,A2,A3,0);
    return vec4;
  }
  else{
    float sol1 = (-1*B+sqrt((pow(B,2) - (4*A*C))))*1./(2*A);
    float sol2 =  (-1*B-sqrt((pow(B,2) - (4*A*C))))*1./(2*A);
    (sol1>sol2)?sol=sol1:sol=sol2;
    
    float A1, A2, A3;
    A3 = sol;
    A1 = d1 + e1*A3;
    A2 = d2 + e2*A3;
    
    TLorentzVector vec4;
    vec4.SetPxPyPzE(A1,A2,A3,0);
    return vec4;;
  }
}

struct JetIDVars
{
  float NHF, NEMF, MUF, CHF, CEMF;
  int NumConst, NumNeutralParticle, CHM;
};

bool getJetID(JetIDVars vars, string jettype="CHS", int year=2018, double eta=0, bool tightLepVeto=true, bool UltraLegacy=false){
  
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
  
  if (jettype!="CHS" && jettype!="PUPPI"){
    cout<<"Don't know your jet type! I know only CHS & PUPPI :D"<<endl;
    return false;
  }
  
  float NHF, NEMF, MUF, CHF, CEMF;
  int NumConst, NumNeutralParticle, CHM;
  
  NHF = vars.NHF; 
  NEMF = vars.NEMF;
  MUF = vars.MUF;
  CHF = vars.CHF;
  CEMF = vars.CEMF;
  NumConst = vars.NumConst;
  NumNeutralParticle = vars.NumNeutralParticle;
  CHM = vars.CHM;
  
  bool JetID = false;
  
  if(!UltraLegacy){
    
    if(year==2018 && jettype=="CHS"){
      
      JetID = ( (fabs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.6 && fabs(eta)<=2.7 && CEMF<0.8 && CHM>0 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.6 && fabs(eta)<=2.7 && CHM>0 && NEMF<0.99 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.7 && fabs(eta)<=3.0 && NEMF>0.02 && NEMF<0.99 && NumNeutralParticle>2) || (fabs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10));
    }
    
    if(year==2018 && jettype=="PUPPI"){
      
      JetID = ( (fabs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.6 && fabs(eta)<=2.7 && CEMF<0.8 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.6 && fabs(eta)<=2.7 && NEMF<0.99 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)<=3.0 && NHF<0.99) || (fabs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticle>2 && NumNeutralParticle<15));
    }
    
    if(year==2017 && jettype=="CHS"){
      
      JetID = ( (fabs(eta)<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.4 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 &&  NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 &&  NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.7 && fabs(eta)<=3.0 && NEMF>0.02 && NEMF<0.99 && NumNeutralParticle>2) || (fabs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticle>10));
    }
    
    if(year==2017 && jettype=="PUPPI"){
      
      JetID = ( (fabs(eta)<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.4 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) ||
 (fabs(eta)>2.4 && fabs(eta)<=2.7 &&  NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 &&  NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.7 && fabs(eta)<=3.0 && NHF<0.99) || (fabs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticle>2 && NumNeutralParticle<15));
    }

    if(year==2016 && jettype=="CHS"){
      
      JetID = ( (fabs(eta)<=2.4 && CEMF<0.90 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.4 && CEMF<0.99 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9  && !tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.7 && fabs(eta)<=3.0 && NEMF>0.01 && NHF<0.98 && NumNeutralParticle>2) || (fabs(eta)>3.0 && NEMF<0.90 && NumNeutralParticle>10));
	}
    
    if(year==2016 && jettype=="PUPPI"){
      
      JetID = ( (fabs(eta)<=2.4 && CEMF<0.9 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.4 && CEMF<0.99 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ));
      if(fabs(eta)>2.7) { JetID = false; }
	}
  }
  
  else {
    
    if(year==2017||year==2018){
      
      if(jettype=="CHS"){
	
	JetID = ( fabs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || ( fabs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || ( fabs(eta)>2.6 && fabs(eta)<=2.7 && CEMF<0.8 && CHM>0 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 ) || ( fabs(eta)>2.6 && fabs(eta)<=2.7 && CHM>0 && NEMF<0.99 && NHF < 0.9 ) || ( fabs(eta)>2.7 && fabs(eta)<=3.0 && NEMF>0.01 && NEMF<0.99 && NumNeutralParticle>1 ) || ( fabs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10) ;
      }
      
      if(jettype=="PUPPI"){
	
	JetID =  ( fabs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || ( fabs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || ( fabs(eta)>2.6 && fabs(eta)<=2.7 && CEMF<0.8 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || ( fabs(eta)>2.6 && fabs(eta)<=2.7 && NEMF<0.99 && NHF < 0.9 && !tightLepVeto ) || ( fabs(eta)>2.7 && fabs(eta)<=3.0 && NHF<0.9999 ) ||( fabs(eta)>3.0 && NEMF<0.90 && NumNeutralParticle>2 ) ;
      }
      // there is a inconsistency between table & lines in https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL
      // table is chosen as it is consistent with the slides https://indico.cern.ch/event/937597/contributions/3940302/attachments/2073315/3481068/ULJetID_UL17_UL18_AK4PUPPI.pdf 
    }
    
    if(year==2016){
      
      if(jettype=="CHS"){
	
	JetID =  ( fabs(eta)<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || ( fabs(eta)<=2.4 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || ( fabs(eta)>2.4 && fabs(eta)<=2.7 && NEMF<0.99 && NHF < 0.9 ) || ( fabs(eta)>2.7 && fabs(eta)<=3.0 && NEMF>0.0 && NEMF<0.99 && NHF<0.9 && NumNeutralParticle>1 ) || ( fabs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10) ;
	
      }
      
      if(jettype=="PUPPI"){
	
	JetID = ( fabs(eta)<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || ( fabs(eta)<=2.4 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || ( fabs(eta)>2.4 && fabs(eta)<=2.7 && NEMF<0.99 && NHF < 0.98 ) || ( fabs(eta)>2.7 && fabs(eta)<=3.0 && NumNeutralParticle>=1 ) || ( fabs(eta)>3.0 && NEMF<0.90 && NumNeutralParticle>2  ) ;
      }
    }	
  }
  
  return JetID;
  
}

bool Muon_TightID(bool muonisGL,bool muonisPF, float muonchi, float muonhit, float muonmst, float muontrkvtx, float muondz, float muonpixhit, float muontrklay){
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Tight_Muon
	bool tightid = false;
	if(muonisGL && muonisPF){
		if(muonchi<10 && muonhit>0 && muonmst>1){
			if(fabs(muontrkvtx)<0.2 && fabs(muondz)<0.5){
				if(muonpixhit>0 && muontrklay>5){
					tightid = true;
				}
			}
		}
	}
	return tightid;
}

bool StoreMuon(pat::Muon muon1, float ptcut, float etacut){
	
	if (((muon1.isTrackerMuon() || muon1.isGlobalMuon()) && (muon1.isPFMuon())) && (muon1.pt()>=ptcut) && (fabs(muon1.eta())<=etacut)) {                                                                
			return true;
	}
	else{
			return false;
		}
}

bool StoreElectron(pat::Electron electron1, float ptcut, float etacut){
	
	GsfTrackRef gsftrk1 = electron1.gsfTrack();                                                                                                      
    if ((!gsftrk1.isNull()) && (electron1.pt()>=ptcut) && (fabs(electron1.eta())<=etacut) && (gsftrk1->ndof()>=9)) {
			return true;
		}
    else{
			return false;
		}
}

TLorentzVector LeptonJet_subtraction(vector<auto> leps, pat::Jet jet, TLorentzVector jet4v){
	
	TLorentzVector newjet4v;
	newjet4v = jet4v;
	
	if (leps.size()>0) {                                                                                           
		for (unsigned int ilep = 0; ilep<leps.size(); ilep++) {
	  
			bool lepmember = false;
			
			for(unsigned int jd = 0 ; jd < leps[ilep].numberOfSourceCandidatePtrs() ; ++jd) {
				
				if(leps[ilep].sourceCandidatePtr(jd).isNonnull() && leps[ilep].sourceCandidatePtr(jd).isAvailable()){
					const reco::Candidate* jcand = leps[ilep].sourceCandidatePtr(jd).get();
				
					for(unsigned int ic = 0 ; ic < jet.numberOfSourceCandidatePtrs() ; ++ic) {  
					
						if(jet.sourceCandidatePtr(ic).isNonnull() && jet.sourceCandidatePtr(ic).isAvailable()){
							const reco::Candidate* icand = jet.sourceCandidatePtr(ic).get();
							if (delta2R(jcand->eta(),jcand->phi(),icand->eta(),icand->phi()) < 0.00001)    
							{
								TLorentzVector tmpvec(jcand->px(),jcand->py(),jcand->pz(),jcand->energy());
								newjet4v = jet4v - tmpvec;
								lepmember = true; 
								break;
							}
						}
					}		
				
				}
								
			}//jd
			
			if(lepmember) break;
			
		}//ilep
	}
    
    return newjet4v;
}

void Read_JEC(double &total_JEC,  double &tmprecpt, 
			  double jeteta, double Rho, bool isData,
			  pat::Jet jet,
			  FactorizedJetCorrector *jecL1Fast, FactorizedJetCorrector *jecL2Relative, FactorizedJetCorrector *jecL3Absolute, FactorizedJetCorrector*jecL2L3Residual)
{
	
    double total_cor =1;
      
    jecL1Fast->setJetPt(tmprecpt); jecL1Fast->setJetA(jet.jetArea()); jecL1Fast->setRho(Rho);jecL1Fast->setJetEta(jeteta);
    double corFactorL1Fast = jecL1Fast->getCorrection();
    total_cor *= corFactorL1Fast;
    tmprecpt = tmprecpt * corFactorL1Fast;
      
    jecL2Relative->setJetPt(tmprecpt); jecL2Relative->setJetEta(jeteta);
    double corFactorL2Relative = jecL2Relative->getCorrection();
    total_cor *= corFactorL2Relative ;
    tmprecpt = tmprecpt * corFactorL2Relative;
      
    jecL3Absolute->setJetPt(tmprecpt); jecL3Absolute->setJetEta(jeteta);
    double corFactorL3Absolute = jecL3Absolute->getCorrection();
    total_cor *= corFactorL3Absolute ;
    tmprecpt = tmprecpt * corFactorL3Absolute;
      
    double corFactorL2L3Residual=1.;
      
    if(isData){
		jecL2L3Residual->setJetPt(tmprecpt); jecL2L3Residual->setJetEta(jeteta);
		corFactorL2L3Residual = jecL2L3Residual->getCorrection();
		total_cor*= corFactorL2L3Residual;
		tmprecpt *=corFactorL2L3Residual;
	}
	
	total_JEC = total_cor;
	
	return;     
}

void Read_JER(std::string mPtResoFile, std::string mPtSFFile, double tmprecpt, TLorentzVector pfjet4v, double Rho, edm::Handle<reco::GenJetCollection>  genjets, vector<double> &SFs)
{
 
	JME::JetResolution resolution;
	resolution = JME::JetResolution(mPtResoFile.c_str());
	JME::JetResolutionScaleFactor res_sf;
	res_sf = JME::JetResolutionScaleFactor(mPtSFFile.c_str());
	
	JME::JetParameters parameters_5 = {{JME::Binning::JetPt, tmprecpt}, {JME::Binning::JetEta, pfjet4v.Eta()}, {JME::Binning::Rho, Rho}};
	double rp = resolution.getResolution(parameters_5);
	double gaus_rp = gRandom->Gaus(0.,rp);
	double sf = res_sf.getScaleFactor(parameters_5, Variation::NOMINAL);
	double sf_up = res_sf.getScaleFactor(parameters_5, Variation::UP);
	double sf_dn = res_sf.getScaleFactor(parameters_5, Variation::DOWN);
	
	bool match = false;
	int match_gen = -1;
		
	for (unsigned get = 0; get<(genjets->size()); get++) {
		TLorentzVector genjet4v((*genjets)[get].px(),(*genjets)[get].py(),(*genjets)[get].pz(), (*genjets)[get].energy());
		if((delta2R(pfjet4v.Rapidity(),pfjet4v.Phi(),genjet4v.Rapidity(),genjet4v.Phi()) < (0.5*0.8)) &&(fabs(tmprecpt-genjet4v.Pt())<(3*fabs(rp)*tmprecpt))){
			match = true;
			match_gen = get;
			break;
		}
	}
		
	if(match && (match_gen>=0)){
	  
		SFs.push_back((sf-1.)*(tmprecpt-(*genjets)[match_gen].pt())*1./tmprecpt);
		SFs.push_back((sf_up-1.)*(tmprecpt-(*genjets)[match_gen].pt())*1./tmprecpt);
		SFs.push_back((sf_dn-1.)*(tmprecpt-(*genjets)[match_gen].pt())*1./tmprecpt);
	  
	}else{
	  
		SFs.push_back(sqrt(max(0.,(sf*sf-1))) * gaus_rp);
		SFs.push_back(sqrt(max(0.,(sf_up*sf_up-1))) * gaus_rp);
		SFs.push_back(sqrt(max(0.,(sf_dn*sf_dn-1))) * gaus_rp);
	}
      	
}

float getEtaForEA(auto obj){
	float eta;
	if(abs(obj->pdgId())==11||abs(obj->pdgId())==22) { eta = obj->superCluster()->eta(); }     
	else { eta = obj->eta(); }
	return eta;    
}

std::unique_ptr<EffectiveAreas> ea_mu_miniiso_, ea_el_miniiso_;

void Read_MiniIsolation(auto obj, double Rho, vector<float> &isovalues)
{
	pat::PFIsolation iso = obj->miniPFIsolation();                                                                                                                                                                                                   
	float chg = iso.chargedHadronIso();                                                                                                                     
	float neu = iso.neutralHadronIso();                                                                                                                     
	float pho = iso.photonIso();                                                                                       
	                                                                                   
	float ea;
	if(abs(obj->pdgId())==13) { ea = ea_mu_miniiso_->getEffectiveArea(fabs(getEtaForEA(obj))); }
	else { ea = ea_el_miniiso_->getEffectiveArea(fabs(getEtaForEA(obj))); }  
	                                                                                    
	float R = 10.0/std::min(std::max(obj->pt(), 50.0),200.0);                                                                      
	ea *= std::pow(R / 0.3, 2);                                                                                                                  	
	float tot = (chg+std::max(0.0,neu+pho-(Rho)*ea));
	
	isovalues.push_back(tot);
	isovalues.push_back(chg);
	isovalues.push_back(neu);
	isovalues.push_back(pho);	
	
	for(unsigned ij=0; ij<isovalues.size(); ij++){
		isovalues[ij] *= 1./obj->pt();
	}
}

std::unique_ptr<EffectiveAreas> ea_el_pfiso_;

void Read_ElePFIsolation(auto obj, double Rho, vector<float> &isovalues)
{
	auto iso = obj->pfIsolationVariables();   
	auto  ea = ea_el_pfiso_->getEffectiveArea(fabs(getEtaForEA(obj)));                                                    
    float val = iso.sumChargedHadronPt + max(0., iso.sumNeutralHadronEt + iso.sumPhotonEt - (Rho)*ea); 
    float val04 = (obj->chargedHadronIso()+std::max(0.0,obj->neutralHadronIso()+obj->photonIso()-(Rho)*ea*16./9.));
    isovalues.push_back(val);
    isovalues.push_back(val04);
    
    for(unsigned ij=0; ij<isovalues.size(); ij++){
		isovalues[ij] *= 1./obj->pt();
	}    
}

//class declaration
//
class Leptop : public edm::EDAnalyzer {
public:
  explicit Leptop(const edm::ParameterSet&);
  ~Leptop();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  void fillmetarray();
  // ----------member data ---------------------------
  int Nevt;
  bool isData;
  bool isMC;
  int year;
  bool isUltraLegacy;
  bool isSoftDrop;
  bool add_prefireweights;
  
  uint nPDFsets;
  
  std::string theRootFileName;
  std::string theHLTTag;
  std::string softdropmass;
  std::string tau1;
  std::string tau2;
  std::string tau3;
  std::string subjets;
  std::string toptagger_DAK8;
  std::string Wtagger_DAK8;
  std::string Ztagger_DAK8;
  std::string Htagger_DAK8;
  std::string bbtagger_DAK8;
  std::string toptagger_PNet;
  std::string Wtagger_PNet;
  std::string Ztagger_PNet;
  std::string Xbbtagger_PNet;
  std::string Xcctagger_PNet;
  std::string Xqqtagger_PNet;
  
  edm::EDGetTokenT<double> tok_Rho_;
  edm::EDGetTokenT<reco::BeamSpot> tok_beamspot_;
  edm::EDGetTokenT<reco::VertexCollection> tok_primaryVertices_;
  edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> tok_sv;
  edm::EDGetTokenT<pat::METCollection>tok_mets_, tok_mets_PUPPI_;
  edm::EDGetTokenT<pat::PackedCandidateCollection>tok_pfcands_;
  edm::EDGetTokenT<edm::View<pat::Jet>>tok_pfjetAK8s_;
  edm::EDGetTokenT<edm::View<pat::Jet>>tok_pfjetAK4s_;
  edm::EDGetTokenT<edm::View<pat::Muon>> tok_muons_;
  edm::EDGetTokenT<edm::View<pat::Electron>> tok_electrons_;
  edm::EDGetTokenT<edm::View<pat::Photon>>tok_photons_;
  edm::EDGetTokenT<edm::View<pat::Tau>>tok_taus_;
  
  //std::unique_ptr<EffectiveAreas> ea_miniiso_;
  
  edm::EDGetTokenT<reco::GenMETCollection>tok_genmets_;
  edm::EDGetTokenT<reco::GenJetCollection>tok_genjetAK8s_;
  edm::EDGetTokenT<reco::GenJetCollection>tok_genjetAK4s_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>>tok_genparticles_;
  edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;
  
  edm::EDGetTokenT<HepMCProduct> tok_HepMC ;
  edm::EDGetTokenT<GenEventInfoProduct> tok_wt_;
  edm::EDGetTokenT<LHEEventProduct> lheEventProductToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileup_;
  
  //edm::EDGetTokenT <edm::ValueMap <bool> > tok_mvaPhoID_FallV2_WP90;
  //edm::EDGetTokenT <edm::ValueMap <bool> > tok_mvaPhoID_FallV2_WP80;
  edm::EDGetTokenT <edm::ValueMap <float> > tok_mvaPhoID_FallV2_raw;

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  
  edm::EDGetTokenT< double > prefweight_token;
  edm::EDGetTokenT< double > prefweightup_token;
  edm::EDGetTokenT< double > prefweightdown_token;
  
  // object cuts //
  int iTag;
  int iTagMET;
  
  double minjPt;
  double minGenPt;
  double maxEta;
  double maxgenEta;
  double AK8PtCut;
  double AK8GenPtCut;
  double minmuPt;
  double minePt;
  double mingmPt;
  double mintauPt;
  
  double beta ;
  double z_cut;
  
  // Root file & tree //
  
  TFile* theFile;
  
  TTree* T1;
  
  // HLTConfigProvider hltConfig_;
  
  unsigned ievt;
  
  static const int njetmx = 20; 
  static const int njetmxAK8 =10;
  static const int npartmx = 50; 
  
  int irunold;
  int irun, ilumi, ifltr, nprim, npvert, ibrnch;
  double event_weight;
  double weights[njetmx];
  
  double Rho ;
  
  double prefiringweight, prefiringweightup, prefiringweightdown;
  
  int npfjetAK8;
  float pfjetAK8pt[njetmxAK8], pfjetAK8y[njetmxAK8], pfjetAK8eta[njetmxAK8], pfjetAK8phi[njetmxAK8], pfjetAK8mass[njetmxAK8];
  
  bool pfjetAK8jetID_tightlepveto[njetmxAK8], pfjetAK8jetID[njetmxAK8];
  
  float pfjetAK8btag_DeepCSV[njetmxAK8];
  float pfjetAK8DeepTag_DAK8_TvsQCD[njetmxAK8], pfjetAK8DeepTag_DAK8_WvsQCD[njetmxAK8], pfjetAK8DeepTag_DAK8_ZvsQCD[njetmxAK8], pfjetAK8DeepTag_DAK8_HvsQCD[njetmxAK8], pfjetAK8DeepTag_DAK8_bbvsQCD[njetmxAK8]; 
  float pfjetAK8DeepTag_PNet_TvsQCD[njetmxAK8], pfjetAK8DeepTag_PNet_WvsQCD[njetmxAK8], pfjetAK8DeepTag_PNet_ZvsQCD[njetmxAK8], pfjetAK8DeepTag_PNet_XbbvsQCD[njetmxAK8], pfjetAK8DeepTag_PNet_XccvsQCD[njetmxAK8], pfjetAK8DeepTag_PNet_XqqvsQCD[njetmxAK8]; 
  
  float pfjetAK8CHF[njetmxAK8], pfjetAK8NHF[njetmxAK8], pfjetAK8MUF[njetmxAK8], pfjetAK8PHF[njetmxAK8], pfjetAK8CEMF[njetmxAK8], pfjetAK8NEMF[njetmxAK8], pfjetAK8EEF[njetmxAK8], pfjetAK8HFHF[njetmxAK8], /*pfjetAK8HFEMF[njetmxAK8],*/ pfjetAK8HOF[njetmxAK8];
  int pfjetAK8CHM[njetmxAK8], pfjetAK8NHM[njetmxAK8], pfjetAK8MUM[njetmxAK8], pfjetAK8PHM[njetmxAK8], pfjetAK8Neucons[njetmxAK8], pfjetAK8Chcons[njetmxAK8], pfjetAK8EEM[njetmxAK8], pfjetAK8HFHM[njetmxAK8];// pfjetAK8HFEMM[njetmxAK8];
  
  float pfjetAK8chrad[njetmxAK8], pfjetAK8pTD[njetmxAK8]; 
  float pfjetAK8sdmass[njetmxAK8], pfjetAK8tau1[njetmxAK8], pfjetAK8tau2[njetmxAK8], pfjetAK8tau3[njetmxAK8];
  
  float pfjetAK8sub1pt[njetmxAK8], pfjetAK8sub1eta[njetmxAK8], pfjetAK8sub1phi[njetmxAK8], pfjetAK8sub1mass[njetmxAK8], pfjetAK8sub1btag[njetmxAK8]; 
  float pfjetAK8sub2pt[njetmxAK8], pfjetAK8sub2eta[njetmxAK8], pfjetAK8sub2phi[njetmxAK8], pfjetAK8sub2mass[njetmxAK8], pfjetAK8sub2btag[njetmxAK8];
  
  float pfjetAK8JEC[njetmxAK8];
  float pfjetAK8reso[njetmxAK8], pfjetAK8resoup[njetmxAK8], pfjetAK8resodn[njetmxAK8];
  float pfjetAK8jesup_pu[njetmx], pfjetAK8jesup_rel[njetmx], pfjetAK8jesup_scale[njetmx], pfjetAK8jesup_total[njetmx], pfjetAK8jesdn_pu[njetmx], pfjetAK8jesdn_rel[njetmx], pfjetAK8jesdn_scale[njetmx], pfjetAK8jesdn_total[njetmx];

  float pfjetAK8jesup_AbsoluteStat[njetmxAK8], pfjetAK8jesdn_AbsoluteStat[njetmxAK8];
  float pfjetAK8jesup_AbsoluteScale[njetmxAK8], pfjetAK8jesdn_AbsoluteScale[njetmxAK8];
  float pfjetAK8jesup_AbsoluteMPFBias[njetmxAK8], pfjetAK8jesdn_AbsoluteMPFBias[njetmxAK8];
  float pfjetAK8jesup_FlavorQCD[njetmxAK8], pfjetAK8jesdn_FlavorQCD[njetmxAK8];
  float pfjetAK8jesup_Fragmentation[njetmxAK8], pfjetAK8jesdn_Fragmentation[njetmxAK8];
  float pfjetAK8jesup_PileUpDataMC[njetmxAK8], pfjetAK8jesdn_PileUpDataMC[njetmxAK8];
  float pfjetAK8jesup_PileUpPtBB[njetmxAK8], pfjetAK8jesdn_PileUpPtBB[njetmxAK8];
  float pfjetAK8jesup_PileUpPtEC1[njetmxAK8], pfjetAK8jesdn_PileUpPtEC1[njetmxAK8];
  float pfjetAK8jesup_PileUpPtEC2[njetmxAK8], pfjetAK8jesdn_PileUpPtEC2[njetmxAK8];
  float pfjetAK8jesup_PileUpPtRef[njetmxAK8], pfjetAK8jesdn_PileUpPtRef[njetmxAK8];
  float pfjetAK8jesup_RelativeFSR[njetmxAK8], pfjetAK8jesdn_RelativeFSR[njetmxAK8];
  float pfjetAK8jesup_RelativeJEREC1[njetmxAK8], pfjetAK8jesdn_RelativeJEREC1[njetmxAK8];
  float pfjetAK8jesup_RelativeJEREC2[njetmxAK8], pfjetAK8jesdn_RelativeJEREC2[njetmxAK8];
  float pfjetAK8jesup_RelativePtBB[njetmxAK8], pfjetAK8jesdn_RelativePtBB[njetmxAK8];
  float pfjetAK8jesup_RelativePtEC1[njetmxAK8], pfjetAK8jesdn_RelativePtEC1[njetmxAK8];
  float pfjetAK8jesup_RelativePtEC2[njetmxAK8], pfjetAK8jesdn_RelativePtEC2[njetmxAK8];
  float pfjetAK8jesup_RelativeBal[njetmxAK8], pfjetAK8jesdn_RelativeBal[njetmxAK8];
  float pfjetAK8jesup_RelativeSample[njetmxAK8], pfjetAK8jesdn_RelativeSample[njetmxAK8];
  float pfjetAK8jesup_RelativeStatEC[njetmxAK8], pfjetAK8jesdn_RelativeStatEC[njetmxAK8];
  float pfjetAK8jesup_RelativeStatFSR[njetmxAK8], pfjetAK8jesdn_RelativeStatFSR[njetmxAK8];
  float pfjetAK8jesup_SinglePionECAL[njetmxAK8], pfjetAK8jesdn_SinglePionECAL[njetmxAK8];
  float pfjetAK8jesup_SinglePionHCAL[njetmxAK8], pfjetAK8jesdn_SinglePionHCAL[njetmxAK8];
  float pfjetAK8jesup_TimePtEta[njetmxAK8], pfjetAK8jesdn_TimePtEta[njetmxAK8];
  float pfjetAK8jesup_Total[njetmxAK8], pfjetAK8jesdn_Total[njetmxAK8];
  
  int npfjetAK4;
  float pfjetAK4pt[njetmx], pfjetAK4eta[njetmx], pfjetAK4y[njetmx], pfjetAK4phi[njetmx], pfjetAK4mass[njetmx];
  float pfjetAK4btag_DeepCSV[njetmx], pfjetAK4btag_DeepFlav[njetmx]; 
  bool pfjetAK4jetID[njetmx], pfjetAK4jetID_tightlepveto[njetmx];
  
  float pfjetAK4reso[njetmx], pfjetAK4resoup[njetmx], pfjetAK4resodn[njetmx];
  
  float pfjetAK4JEC[njetmx];
  
  int pfjetAK4hadronflav[njetmx], pfjetAK4partonflav[njetmx];
  int pfjetAK4Ncons[njetmx];
  float pfjetAK4qgl[njetmx], pfjetAK4PUID[njetmx];
  
  float pfjetAK4jesup_AbsoluteStat[njetmx], pfjetAK4jesdn_AbsoluteStat[njetmx];
  float pfjetAK4jesup_AbsoluteScale[njetmx], pfjetAK4jesdn_AbsoluteScale[njetmx];
  float pfjetAK4jesup_AbsoluteMPFBias[njetmx], pfjetAK4jesdn_AbsoluteMPFBias[njetmx];
  float pfjetAK4jesup_FlavorQCD[njetmx], pfjetAK4jesdn_FlavorQCD[njetmx];
  float pfjetAK4jesup_Fragmentation[njetmx], pfjetAK4jesdn_Fragmentation[njetmx];
  float pfjetAK4jesup_PileUpDataMC[njetmx], pfjetAK4jesdn_PileUpDataMC[njetmx];
  float pfjetAK4jesup_PileUpPtBB[njetmx], pfjetAK4jesdn_PileUpPtBB[njetmx];
  float pfjetAK4jesup_PileUpPtEC1[njetmx], pfjetAK4jesdn_PileUpPtEC1[njetmx];
  float pfjetAK4jesup_PileUpPtEC2[njetmx], pfjetAK4jesdn_PileUpPtEC2[njetmx];
  float pfjetAK4jesup_PileUpPtRef[njetmx], pfjetAK4jesdn_PileUpPtRef[njetmx];
  float pfjetAK4jesup_RelativeFSR[njetmx], pfjetAK4jesdn_RelativeFSR[njetmx];
  float pfjetAK4jesup_RelativeJEREC1[njetmx], pfjetAK4jesdn_RelativeJEREC1[njetmx];
  float pfjetAK4jesup_RelativeJEREC2[njetmx], pfjetAK4jesdn_RelativeJEREC2[njetmx];
  float pfjetAK4jesup_RelativePtBB[njetmx], pfjetAK4jesdn_RelativePtBB[njetmx];
  float pfjetAK4jesup_RelativePtEC1[njetmx], pfjetAK4jesdn_RelativePtEC1[njetmx];
  float pfjetAK4jesup_RelativePtEC2[njetmx], pfjetAK4jesdn_RelativePtEC2[njetmx];
  float pfjetAK4jesup_RelativeBal[njetmx], pfjetAK4jesdn_RelativeBal[njetmx];
  float pfjetAK4jesup_RelativeSample[njetmx], pfjetAK4jesdn_RelativeSample[njetmx];
  float pfjetAK4jesup_RelativeStatEC[njetmx], pfjetAK4jesdn_RelativeStatEC[njetmx];
  float pfjetAK4jesup_RelativeStatFSR[njetmx], pfjetAK4jesdn_RelativeStatFSR[njetmx];
  float pfjetAK4jesup_SinglePionECAL[njetmx], pfjetAK4jesdn_SinglePionECAL[njetmx];
  float pfjetAK4jesup_SinglePionHCAL[njetmx], pfjetAK4jesdn_SinglePionHCAL[njetmx];
  float pfjetAK4jesup_TimePtEta[njetmx], pfjetAK4jesdn_TimePtEta[njetmx];
  float pfjetAK4jesup_Total[njetmx], pfjetAK4jesdn_Total[njetmx];
  
  static const int ngenjetAK8mx =10;
  
  int ngenjetAK8;
  float genjetAK8pt[njetmxAK8], genjetAK8eta[njetmxAK8], genjetAK8phi[njetmxAK8], genjetAK8mass[njetmxAK8], genjetAK8sdmass[njetmxAK8]; 
  int genjetAK8hadronflav[njetmxAK8], genjetAK8partonflav[njetmxAK8];
  
  int ngenjetAK4;
  float genjetAK4pt[njetmx], genjetAK4eta[njetmx], genjetAK4phi[njetmx], genjetAK4mass[njetmx];
  int genjetAK4hadronflav[njetmx], genjetAK4partonflav[njetmx];
  
  int ngenparticles;
  int genpartstatus[npartmx], genpartpdg[npartmx], genpartmompdg[npartmx], genpartgrmompdg[npartmx], genpartmomid[npartmx], genpartdaugno[npartmx];
  float genpartpt[npartmx], genparteta[npartmx], genpartphi[npartmx], genpartm[npartmx]; //genpartq[npartmx];
  bool genpartfromhard[npartmx], genpartfromhardbFSR[npartmx], genpartisPromptFinalState[npartmx], genpartisLastCopyBeforeFSR[npartmx];
  
  static const int nlhemax = 10;
  int nLHEparticles;
  float LHEpartpt[nlhemax], LHEparteta[nlhemax], LHEpartphi[nlhemax], LHEpartm[nlhemax];
  int LHEpartpdg[nlhemax];
  
  static const int nlhescalemax = 9;
  int nLHEScaleWeights;
  float LHEScaleWeights[nlhescalemax];
  
  static const int nlhepdfmax = 103; // be consistent with nPDFsets (nlhepdfmax should be >= nPDFsets)
  int nLHEPDFWeights;
  float LHEPDFWeights[nlhepdfmax];
  
  static const int nalpsmax = 3;
  int nLHEAlpsWeights;
  float LHEAlpsWeights[nalpsmax];
  
  static const int nlhepsmax = 8;
  int nLHEPSWeights;
  float LHEPSWeights[nlhepsmax];
  
  double event_weight_LHE;
  
  float miset , misphi , sumEt, misetsig;
  float miset_PUPPI , misphi_PUPPI , sumEt_PUPPI, misetsig_PUPPI;
  float genmiset, genmisphi, genmisetsig;
  
  int nmuons;
  
  float muonminchiso[njetmx], muonminnhiso[njetmx], muonminphiso[njetmx], muonminisoall[njetmx]; 
  float muoncharge[njetmx], muonp[njetmx], muonpt[njetmx], muoneta[njetmx], muonphi[njetmx], muondz[njetmx], muonpter[njetmx], muonchi[njetmx], muonecal[njetmx], muonhcal[njetmx]; //muonemiso[njetmx], muonhadiso[njetmx], muontkpt03[njetmx], muontkpt05[njetmx];
  
  float muonposmatch[njetmx], muontrkink[njetmx], muonsegcom[njetmx], muonpfiso[njetmx], muontrkvtx[njetmx], muonhit[njetmx], muonpixhit[njetmx], muonmst[njetmx], muontrklay[njetmx], muonvalfrac[njetmx],mudxy_sv[njetmx];
  int muonndf[njetmx];
  
  bool muonisPF[njetmx], muonisGL[njetmx], muonisTRK[njetmx];
  bool muonisGoodGL[njetmx], muonisTight[njetmx], muonisHighPt[njetmx], muonisHighPttrk[njetmx], muonisMed[njetmx], muonisMedPr[njetmx], muonisLoose[njetmx];
  
  int nelecs;
  bool elmvaid_Fallv2WP90[njetmx], elmvaid_Fallv2WP90_noIso[njetmx];
  bool elmvaid_Fallv2WP80[njetmx], elmvaid_Fallv2WP80_noIso[njetmx];

  float elcharge[njetmx], elpt[njetmx], eleta[njetmx], elphi[njetmx], ele[njetmx], elp[njetmx], eldxytrk[njetmx], eldxy_sv[njetmx], eldztrk[njetmx],elhovere[njetmx], elqovrper[njetmx], elchi[njetmx]; //elemiso03[njetmx], elhadiso03[njetmx], elemiso04[njetmx], elhadiso04[njetmx];
  float eleoverp[njetmx], elietaieta[njetmx], eletain[njetmx], elphiin[njetmx], elfbrem[njetmx]; 
  float elnohits[njetmx], elmisshits[njetmx];
  float elpfiso_drcor[njetmx];
  float elpfiso_eacor[njetmx];
  float elpfiso04_eacor[njetmx];
  int elndf[njetmx];
  
  float elsupcl_eta[njetmx]; 
  float elsupcl_phi[njetmx]; 
  float elsupcl_rawE[njetmx]; 
  float elsigmaieta[njetmx];
  float elsigmaiphi[njetmx];
  float elr9full[njetmx];
  float elsupcl_etaw[njetmx];
  float elsupcl_phiw[njetmx];
  float elhcaloverecal[njetmx];
  float elcloctftrkn[njetmx];
  float elcloctftrkchi2[njetmx];
  float ele1x5bye5x5[njetmx];
  float elnormchi2[njetmx];
  float elhitsmiss[njetmx];
  float eltrkmeasure[njetmx];
  float elconvtxprob[njetmx];
  float elecloverpout[njetmx];
  float elecaletrkmomentum[njetmx];
  float eldeltaetacltrkcalo[njetmx];
  float elsupcl_preshvsrawe[njetmx];
  float elpfisolsumphet[njetmx];
  float elpfisolsumchhadpt[njetmx];
  float elpfsiolsumneuhadet[njetmx];
  float elminchiso[njetmx];
  float elminnhiso[njetmx];
  float elminphiso[njetmx];
  float elminisoall[njetmx]; 
  
  int nphotons;
  bool phomvaid_Fall17V2_WP90[njetmx];
  bool phomvaid_Fall17V2_WP80[njetmx];
  float phomvaid_Fall17V2_raw[njetmx];
  bool phomvaid_Spring16V1_WP90[njetmx];
  bool phomvaid_Spring16V1_WP80[njetmx];
  float phoe[njetmx];
  float phoeta[njetmx];
  float phophi[njetmx];
  float phoe1by9[njetmx];
  float phoe9by25[njetmx];
  float phohadbyem[njetmx];
  float photrkiso[njetmx];
  float phoemiso[njetmx];
  float phohadiso[njetmx];
  float phochhadiso[njetmx];
  float phoneuhadiso[njetmx];
  float phoPUiso[njetmx];
  float phophoiso[njetmx];
  float phoietaieta[njetmx];
  
  int ntaus;
  float taupt[njetmx];
  float taueta[njetmx];
  float tauphi[njetmx];
  float taue[njetmx];
  bool tauisPF[njetmx];
  float taudxy[njetmx];
  float taudz[njetmx];
  int taucharge[njetmx];
  int taudecayMode[njetmx];
  bool taudecayModeinding[njetmx];
  bool taudecayModeindingNewDMs[njetmx];

  int taumuiso[njetmx];
  float taueiso_raw[njetmx];
  int taueiso[njetmx];
  float taueiso2018_raw[njetmx];
  int taueiso2018[njetmx];

  float taujetiso_deeptau2017v2p1_raw[njetmx];
  int taujetiso_deeptau2017v2p1[njetmx];
  float taueiso_deeptau2017v2p1_raw[njetmx];
  int taueiso_deeptau2017v2p1[njetmx];
  float taumuiso_deeptau2017v2p1_raw[njetmx];
  int taumuiso_deeptau2017v2p1[njetmx];

  float taurawiso[njetmx];
  float taurawisodR03[njetmx];
  float taupuCorr[njetmx];
  float tauleadtrkpt[njetmx];
  float tauleadtrketa[njetmx];
  float tauleadtrkphi[njetmx];
  float tauleadtrkdxy[njetmx];
  float tauleadtrkdz[njetmx];

  // Trigger Info //
  
  int ntrigobjs;
  float trigobjpt[njetmx], trigobjeta[njetmx],trigobjphi[njetmx], trigobjmass[njetmx];
  bool trigobjHLT[njetmx], trigobjL1[njetmx],  trigobjBoth[njetmx];
  int  trigobjIhlt[njetmx], trigobjpdgId[njetmx];
  
  // Collision Info //
  
  float qscale, Generator_x1, Generator_x2, Generator_xpdf1, Generator_xpdf2, Generator_scalePDF;
  int Generator_id1, Generator_id2;
  
  float wtfact;
  int npu_vert;
  int npu_vert_true;
  
  // HL triggers //
  
  static const int nHLTmx = 13;
  const char *hlt_name[nHLTmx] = {"HLT_IsoMu24_v","HLT_Mu50_v","HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v","HLT_Ele115_CaloIdVT_GsfTrkIdT_v", "HLT_Ele40_WPTight_Gsf_v",  // single-lepton triggers
								  "HLT_Mu37_Ele27_CaloIdL_MW_v", "HLT_Mu27_Ele37_CaloIdL_MW_v", "HLT_Mu37_TkMu27_v", "HLT_DoubleEle25_CaloIdL_MW_v", // double-lepton triggers
								  "HLT_AK8PFJet500_v", "HLT_PFJet500_v", "HLT_PFHT1050_v", "HLT_Photon200_v"}; //reference & backup triggers
  
  bool hlt_IsoMu24, hlt_Mu50, hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165, hlt_Ele115_CaloIdVT_GsfTrkIdT, hlt_Ele40_WPTight_Gsf, 
       hlt_Mu37_Ele27_CaloIdL_MW, hlt_Mu27_Ele37_CaloIdL_MW, hlt_Mu37_TkMu27, hlt_DoubleEle25_CaloIdL_MW, 
       hlt_AK8PFJet500, hlt_PFJet500, hlt_HT1050, hlt_Photon200;
  
  int trig_value;
  
  HLTPrescaleProvider hltPrescaleProvider_;
    
  // ---- Jet Corrector Parameter ---- //
  
  JetCorrectorParameters *L1FastAK4, *L2RelativeAK4, *L3AbsoluteAK4, *L2L3ResidualAK4;
  vector<JetCorrectorParameters> vecL1FastAK4, vecL2RelativeAK4, vecL3AbsoluteAK4, vecL2L3ResidualAK4;
  FactorizedJetCorrector *jecL1FastAK4, *jecL2RelativeAK4, *jecL3AbsoluteAK4, *jecL2L3ResidualAK4;
  
  JetCorrectorParameters *L1FastAK8, *L2RelativeAK8, *L3AbsoluteAK8, *L2L3ResidualAK8;
  vector<JetCorrectorParameters> vecL1FastAK8, vecL2RelativeAK8, vecL3AbsoluteAK8, vecL2L3ResidualAK8;
  FactorizedJetCorrector *jecL1FastAK8, *jecL2RelativeAK8, *jecL3AbsoluteAK8, *jecL2L3ResidualAK8;
  
  // ---- Jet Corrector Parameter End---- //

  // ---- Jet Resolution Parameter ---- //
  
  std::string mJECL1FastFileAK4, mJECL2RelativeFileAK4, mJECL3AbsoluteFileAK4, mJECL2L3ResidualFileAK4, mJECL1FastFileAK8, mJECL2RelativeFileAK8, mJECL3AbsoluteFileAK8, mJECL2L3ResidualFileAK8;
  std::string mPtResoFileAK4, mPtResoFileAK8, mPtSFFileAK4, mPtSFFileAK8;
  
  std::string mJECUncFileAK4;
  std::vector<JetCorrectionUncertainty*> vsrc ;
  
  std::string mJECUncFileAK8;
  std::vector<JetCorrectionUncertainty*> vsrcAK8 ;
  
  // ---- Jet Resolution Parameter End---- //
  
  // Electron MVA ID //
  
  std::string melectronID_isowp90, melectronID_noisowp90;
  std::string melectronID_isowp80, melectronID_noisowp80;
  
  // Photon MVA ID //
  
  std::string mPhoID_FallV2_WP90, mPhoID_FallV2_WP80;
  std::string mPhoID_SpringV1_WP90, mPhoID_SpringV1_WP80;
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

Leptop::Leptop(const edm::ParameterSet& pset):
  hltPrescaleProvider_(pset, consumesCollector(), *this)  
{
  //now do what ever initialization is needed
  
  edm::Service<TFileService> fs;
  /*
  isMC      = pset.getUntrackedParameter<bool>("MonteCarlo", false);
  isData = !isMC;
  year		= pset.getUntrackedParameter<int>("YEAR", 2018);
  isUltraLegacy = pset.getUntrackedParameter<bool>("UltraLegacy", false);
  isSoftDrop      = pset.getUntrackedParameter<bool>("SoftDrop_ON",false);
  theHLTTag = pset.getUntrackedParameter<string>("HLTTag", "HLT");
  
  theRootFileName = pset.getUntrackedParameter<string>("RootFileName");
  
  triggerBits_ = consumes<edm::TriggerResults> ( pset.getParameter<edm::InputTag>("bits"));
  triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(pset.getParameter<edm::InputTag>("objects"));
  triggerPrescales_ = consumes<pat::PackedTriggerPrescales>(pset.getParameter<edm::InputTag>("prescales"));
  
  minPt = pset.getUntrackedParameter<double>("minPt",25.);
  minGenPt = pset.getUntrackedParameter<double>("minGenPt",15.);
  maxEta = pset.getUntrackedParameter<double>("maxEta",3.);
  maxgenEta = pset.getUntrackedParameter<double>("maxgenEta",3.);
  AK8PtCut = pset.getUntrackedParameter<double>("AK8PtCut",180.);
  AK8GenPtCut = pset.getUntrackedParameter<double>("AK8GenPtCut",150.);
  
  softdropmass = pset.getUntrackedParameter<string>("softdropmass");
  tau1 = pset.getUntrackedParameter<string>("tau1");
  tau2 = pset.getUntrackedParameter<string>("tau2");
  tau3 = pset.getUntrackedParameter<string>("tau3");
  subjets = pset.getUntrackedParameter<string>("subjets");
  toptagger = pset.getUntrackedParameter<string>("toptagger");
  Wtagger = pset.getUntrackedParameter<string>("Wtagger");
  Ztagger = pset.getUntrackedParameter<string>("Ztagger");
  Htagger = pset.getUntrackedParameter<string>("Htagger");
  bbtagger = pset.getUntrackedParameter<string>("bbtagger");
 
  tok_beamspot_ = consumes<reco::BeamSpot> (pset.getParameter<edm::InputTag>("Beamspot"));
  tok_primaryVertices_ =consumes<reco::VertexCollection>( pset.getParameter<edm::InputTag>("PrimaryVertices"));
  tok_sv =consumes<reco::VertexCompositePtrCandidateCollection>( pset.getParameter<edm::InputTag>("SecondaryVertices"));
  tok_Rho_ = consumes<double>(pset.getParameter<edm::InputTag>("PFRho"));
  ea_miniiso_.reset(new EffectiveAreas((pset.getParameter<edm::FileInPath>("EAFile_MiniIso")).fullPath()));

  relative_ = pset.getParameter<bool>("relative");
  
  tok_mets_= consumes<pat::METCollection> ( pset.getParameter<edm::InputTag>("PFMet"));
  tok_pfcands_ = consumes<pat::PackedCandidateCollection>( pset.getParameter<edm::InputTag>("pfCands"));
  tok_muons_ = consumes<edm::View<pat::Muon>> ( pset.getParameter<edm::InputTag>("Muons"));
  tok_electrons_ = consumes<edm::View<pat::Electron>> ( pset.getParameter<edm::InputTag>("Electrons"));
  //tok_photons_ = consumes<edm::View<pat::Photon>>  ( pset.getParameter<edm::InputTag>("Photons"));
  tok_pfjetAK8s_= consumes<edm::View<pat::Jet>>( pset.getParameter<edm::InputTag>("PFJetsAK8"));
  tok_pfjetAK4s_= consumes<edm::View<pat::Jet>>( pset.getParameter<edm::InputTag>("PFJetsAK4"));
  
  if(isMC){
	tok_genmets_= consumes<reco::GenMETCollection> ( pset.getParameter<edm::InputTag>("GENMet"));
    tok_genjetAK8s_= consumes<reco::GenJetCollection>( pset.getParameter<edm::InputTag>("GENJetAK8"));
    tok_genjetAK4s_= consumes<reco::GenJetCollection>( pset.getParameter<edm::InputTag>("GENJetAK4"));
    tok_genparticles_ = consumes<std::vector<reco::GenParticle>>( pset.getParameter<edm::InputTag>("GenParticles"));
    jetFlavourInfosToken_ = consumes<reco::JetFlavourInfoMatchingCollection>(pset.getParameter<edm::InputTag>("jetFlavourInfos"));
  }
  
  melectronID_isowp90       = pset.getParameter<std::string>("electronID_isowp90");
  melectronID_noisowp90     = pset.getParameter<std::string>("electronID_noisowp90");
  melectronID_isowp80       = pset.getParameter<std::string>("electronID_isowp80");
  melectronID_noisowp80     = pset.getParameter<std::string>("electronID_noisowp80");

  mJECL1FastFileAK4         = pset.getParameter<std::string>("jecL1FastFileAK4");
  mJECL1FastFileAK8         = pset.getParameter<std::string>("jecL1FastFileAK8");
  mJECL2RelativeFileAK4     = pset.getParameter<std::string>("jecL2RelativeFileAK4");
  mJECL2RelativeFileAK8     = pset.getParameter<std::string>("jecL2RelativeFileAK8");
  mJECL3AbsoluteFileAK4     = pset.getParameter<std::string>("jecL3AbsoluteFileAK4");
  mJECL3AbsoluteFileAK8     = pset.getParameter<std::string> ("jecL3AbsoluteFileAK8");
  mJECL2L3ResidualFileAK4   = pset.getParameter<std::string> ("jecL2L3ResidualFileAK4");
  mJECL2L3ResidualFileAK8   = pset.getParameter<std::string> ("jecL2L3ResidualFileAK8");
  
  mPtResoFileAK4  = pset.getParameter<std::string>("PtResoFileAK4");
  mPtResoFileAK8  = pset.getParameter<std::string>("PtResoFileAK8");
  mPtSFFileAK4  = pset.getParameter<std::string>("PtSFFileAK4");
  mPtSFFileAK8  = pset.getParameter<std::string>("PtSFFileAK8");
  
  mJECUncFileAK4 = pset.getParameter<std::string>("JECUncFileAK4");
  mJECUncFileAK8 = pset.getParameter<std::string>("JECUncFileAK8");
  
  if(isMC){    
    tok_HepMC = consumes<HepMCProduct>(pset.getParameter<edm::InputTag>("Generator"));
    tok_wt_ = consumes<GenEventInfoProduct>(pset.getParameter<edm::InputTag>("Generator")) ;
    lheEventProductToken_ = consumes<LHEEventProduct>(pset.getParameter<edm::InputTag>("LHEEventProductInputTag")) ;
    pileup_ = consumes<std::vector<PileupSummaryInfo> >(pset.getParameter<edm::InputTag>("slimmedAddPileupInfo"));
  } 
  
  beta = pset.getUntrackedParameter<double>("beta",0.);
  z_cut = pset.getUntrackedParameter<double>("z_cut",0.1); 
  */
  // old begin //
  isData    = pset.getUntrackedParameter<bool>("Data",false);
  isMC      = pset.getUntrackedParameter<bool>("MonteCarlo", false);
  year		= pset.getUntrackedParameter<int>("YEAR", 2018);
  isUltraLegacy = pset.getUntrackedParameter<bool>("UltraLegacy", false);
  isSoftDrop      = pset.getUntrackedParameter<bool>("SoftDrop_ON",false);
  theRootFileName = pset.getUntrackedParameter<string>("RootFileName");
  theHLTTag = pset.getUntrackedParameter<string>("HLTTag", "HLT");
  add_prefireweights = pset.getUntrackedParameter<bool>("add_prefireweights", false);
  
  minjPt = pset.getUntrackedParameter<double>("minjPt",25.);
  minGenPt = pset.getUntrackedParameter<double>("minGenPt",15.);
  AK8PtCut = pset.getUntrackedParameter<double>("AK8PtCut",180.);
  AK8GenPtCut = pset.getUntrackedParameter<double>("AK8GenPtCut",150.);
  minmuPt = pset.getUntrackedParameter<double>("minmuPt",10.);
  minePt = pset.getUntrackedParameter<double>("minePt",10.);
  mingmPt = pset.getUntrackedParameter<double>("mingmPt",10.);
  mintauPt = pset.getUntrackedParameter<double>("mintauPt",10.);
  
  maxEta = pset.getUntrackedParameter<double>("maxEta",3.);
  maxgenEta = pset.getUntrackedParameter<double>("maxgenEta",3.);
 
  softdropmass = pset.getUntrackedParameter<string>("softdropmass");
  beta = pset.getUntrackedParameter<double>("beta",0);
  z_cut = pset.getUntrackedParameter<double>("z_cut",0.1); 
  tau1 = pset.getUntrackedParameter<string>("tau1");
  tau2 = pset.getUntrackedParameter<string>("tau2");
  tau3 = pset.getUntrackedParameter<string>("tau3");
  subjets = pset.getUntrackedParameter<string>("subjets");
  toptagger_DAK8 = pset.getUntrackedParameter<string>("toptagger_DAK8");
  Wtagger_DAK8 = pset.getUntrackedParameter<string>("Wtagger_DAK8");
  Ztagger_DAK8 = pset.getUntrackedParameter<string>("Ztagger_DAK8");
  Htagger_DAK8 = pset.getUntrackedParameter<string>("Htagger_DAK8");
  bbtagger_DAK8 = pset.getUntrackedParameter<string>("bbtagger_DAK8");  
  toptagger_PNet = pset.getUntrackedParameter<string>("toptagger_PNet");
  Wtagger_PNet = pset.getUntrackedParameter<string>("Wtagger_PNet");
  Ztagger_PNet = pset.getUntrackedParameter<string>("Ztagger_PNet");
  Xbbtagger_PNet = pset.getUntrackedParameter<string>("Xbbtagger_PNet");
  Xcctagger_PNet = pset.getUntrackedParameter<string>("Xcctagger_PNet");  
  Xqqtagger_PNet = pset.getUntrackedParameter<string>("Xqqtagger_PNet");  
  
  ea_mu_miniiso_.reset(new EffectiveAreas((pset.getParameter<edm::FileInPath>("EAFile_MuonMiniIso")).fullPath()));
  ea_el_miniiso_.reset(new EffectiveAreas((pset.getParameter<edm::FileInPath>("EAFile_EleMiniIso")).fullPath()));
  ea_el_pfiso_.reset(new EffectiveAreas((pset.getParameter<edm::FileInPath>("EAFile_ElePFIso")).fullPath()));
 
  tok_beamspot_ = consumes<reco::BeamSpot> (pset.getParameter<edm::InputTag>("Beamspot"));
  tok_primaryVertices_ =consumes<reco::VertexCollection>( pset.getParameter<edm::InputTag>("PrimaryVertices"));
  tok_sv =consumes<reco::VertexCompositePtrCandidateCollection>( pset.getParameter<edm::InputTag>("SecondaryVertices"));
  tok_Rho_ = consumes<double>(pset.getParameter<edm::InputTag>("PFRho"));
     
  tok_mets_= consumes<pat::METCollection> ( pset.getParameter<edm::InputTag>("PFMet"));
  tok_mets_PUPPI_ = consumes<pat::METCollection> ( pset.getParameter<edm::InputTag>("PuppiMet"));
  tok_genmets_= consumes<reco::GenMETCollection> ( pset.getParameter<edm::InputTag>("GENMet"));
  
  tok_pfcands_ = consumes<pat::PackedCandidateCollection>( pset.getParameter<edm::InputTag>("pfCands"));
  
  tok_muons_ = consumes<edm::View<pat::Muon>> ( pset.getParameter<edm::InputTag>("Muons"));
  tok_electrons_ = consumes<edm::View<pat::Electron>> ( pset.getParameter<edm::InputTag>("Electrons"));
  tok_photons_ = consumes<edm::View<pat::Photon>>  ( pset.getParameter<edm::InputTag>("Photons"));
  tok_taus_ = consumes<edm::View<pat::Tau>>  ( pset.getParameter<edm::InputTag>("Taus"));
  
  tok_pfjetAK8s_= consumes<edm::View<pat::Jet>>( pset.getParameter<edm::InputTag>("PFJetsAK8"));
  tok_pfjetAK4s_= consumes<edm::View<pat::Jet>>( pset.getParameter<edm::InputTag>("PFJetsAK4"));
  
  if(isMC){
	  
    tok_genjetAK8s_= consumes<reco::GenJetCollection>( pset.getParameter<edm::InputTag>("GENJetAK8"));
    tok_genjetAK4s_= consumes<reco::GenJetCollection>( pset.getParameter<edm::InputTag>("GENJetAK4"));
    tok_genparticles_ = consumes<std::vector<reco::GenParticle>>( pset.getParameter<edm::InputTag>("GenParticles"));
    jetFlavourInfosToken_ = consumes<reco::JetFlavourInfoMatchingCollection>(pset.getParameter<edm::InputTag>("jetFlavourInfos"));
    
    tok_HepMC = consumes<HepMCProduct>(pset.getParameter<edm::InputTag>("Generator"));
    tok_wt_ = consumes<GenEventInfoProduct>(pset.getParameter<edm::InputTag>("Generator")) ;
    lheEventProductToken_ = consumes<LHEEventProduct>(pset.getParameter<edm::InputTag>("LHEEventProductInputTag")) ;
    pileup_ = consumes<std::vector<PileupSummaryInfo> >(pset.getParameter<edm::InputTag>("slimmedAddPileupInfo"));
    nPDFsets      = pset.getUntrackedParameter<uint>("nPDFsets", 103);
    
  }
  
  melectronID_isowp90       = pset.getParameter<std::string>("electronID_isowp90");
  melectronID_noisowp90     = pset.getParameter<std::string>("electronID_noisowp90");
  melectronID_isowp80       = pset.getParameter<std::string>("electronID_isowp80");
  melectronID_noisowp80     = pset.getParameter<std::string>("electronID_noisowp80");
  
  mPhoID_FallV2_WP90       = pset.getParameter<std::string>("PhoID_FallV2_WP90");
  mPhoID_FallV2_WP80       = pset.getParameter<std::string>("PhoID_FallV2_WP80");
  mPhoID_SpringV1_WP90       = pset.getParameter<std::string>("PhoID_SpringV1_WP90");
  mPhoID_SpringV1_WP80       = pset.getParameter<std::string>("PhoID_SpringV1_WP80");
  
  //tok_mvaPhoID_FallV2_WP90 = consumes<edm::ValueMap <bool> >(pset.getParameter<edm::InputTag>("label_mvaPhoID_FallV2_WP90"));
  //tok_mvaPhoID_FallV2_WP80 = consumes<edm::ValueMap <bool> >(pset.getParameter<edm::InputTag>("label_mvaPhoID_FallV2_WP80"));
  tok_mvaPhoID_FallV2_raw = consumes<edm::ValueMap <float> >(pset.getParameter<edm::InputTag>("label_mvaPhoID_FallV2_Value"));
  
  mJECL1FastFileAK4         = pset.getParameter<std::string>("jecL1FastFileAK4");
  mJECL1FastFileAK8         = pset.getParameter<std::string>("jecL1FastFileAK8");
  mJECL2RelativeFileAK4     = pset.getParameter<std::string>("jecL2RelativeFileAK4");
  mJECL2RelativeFileAK8     = pset.getParameter<std::string>("jecL2RelativeFileAK8");
  mJECL3AbsoluteFileAK4     = pset.getParameter<std::string>("jecL3AbsoluteFileAK4");
  mJECL3AbsoluteFileAK8     = pset.getParameter<std::string> ("jecL3AbsoluteFileAK8");
  mJECL2L3ResidualFileAK4   = pset.getParameter<std::string> ("jecL2L3ResidualFileAK4");
  mJECL2L3ResidualFileAK8   = pset.getParameter<std::string> ("jecL2L3ResidualFileAK8");
  
  mPtResoFileAK4  = pset.getParameter<std::string>("PtResoFileAK4");
  mPtResoFileAK8  = pset.getParameter<std::string>("PtResoFileAK8");
  mPtSFFileAK4  = pset.getParameter<std::string>("PtSFFileAK4");
  mPtSFFileAK8  = pset.getParameter<std::string>("PtSFFileAK8");
  
  mJECUncFileAK4 = pset.getParameter<std::string>("JECUncFileAK4");
  mJECUncFileAK8 = pset.getParameter<std::string>("JECUncFileAK8");
  
  triggerBits_ = consumes<edm::TriggerResults> ( pset.getParameter<edm::InputTag>("bits"));
  triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(pset.getParameter<edm::InputTag>("objects"));
  triggerPrescales_ = consumes<pat::PackedTriggerPrescales>(pset.getParameter<edm::InputTag>("prescales"));
  if(add_prefireweights){
	prefweight_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"));
	prefweightup_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbUp"));
	prefweightdown_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbDown"));
	}
  
  //old end //
  
  theFile = new TFile(theRootFileName.c_str(), "RECREATE");
  theFile->cd();
  
  T1 = new TTree("Events", "XtoYH");
 
  T1->Branch("irun", &irun, "irun/I");  
  T1->Branch("ilumi", &ilumi, "ilumi/I");  
  
  // primary vertices //
  
  T1->Branch("ievt", &ievt, "ievt/i");
  T1->Branch("nprim", &nprim, "nprim/I");
  T1->Branch("npvert", &npvert, "npvert/I");
  
  // energy density //
  
  T1->Branch("Rho", &Rho, "Rho/D") ;
  
  // trigger info //
  
  T1->Branch("trig_value",&trig_value,"trig_value/I");  
  
  T1->Branch("hlt_IsoMu24",&hlt_IsoMu24,"hlt_IsoMu24/O");
  T1->Branch("hlt_Mu50",&hlt_Mu50,"hlt_Mu50/O");
  T1->Branch("hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165",&hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165,"hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165/O");
  T1->Branch("hlt_Ele115_CaloIdVT_GsfTrkIdT",&hlt_Ele115_CaloIdVT_GsfTrkIdT,"hlt_Ele115_CaloIdVT_GsfTrkIdT/O");
  T1->Branch("hlt_Ele40_WPTight_Gsf",&hlt_Ele40_WPTight_Gsf,"hlt_Ele40_WPTight_Gsf/O");
  T1->Branch("hlt_Mu37_Ele27_CaloIdL_MW",&hlt_Mu37_Ele27_CaloIdL_MW,"hlt_Mu37_Ele27_CaloIdL_MW/O");
  T1->Branch("hlt_Mu27_Ele37_CaloIdL_MW",&hlt_Mu27_Ele37_CaloIdL_MW,"hlt_Mu27_Ele37_CaloIdL_MW/O");
  T1->Branch("hlt_Mu37_TkMu27",&hlt_Mu37_TkMu27,"hlt_Mu37_TkMu27/O");
  T1->Branch("hlt_DoubleEle25_CaloIdL_MW",&hlt_DoubleEle25_CaloIdL_MW,"hlt_DoubleEle25_CaloIdL_MW/O");
  T1->Branch("hlt_AK8PFJet500",&hlt_AK8PFJet500,"hlt_AK8PFJet500/O");
  T1->Branch("hlt_PFJet500",&hlt_PFJet500,"hlt_PFJet500/O");
  T1->Branch("hlt_HT1050",&hlt_HT1050,"hlt_HT1050/O");
  T1->Branch("hlt_Photon200",&hlt_Photon200,"hlt_Photon200/O");
  
  
  T1->Branch("ntrigobjs",&ntrigobjs,"ntrigobjs/I");
  T1->Branch("trigobjpt",trigobjpt,"trigobjpt[ntrigobjs]/F");
  T1->Branch("trigobjeta",trigobjeta,"trigobjeta[ntrigobjs]/F");
  T1->Branch("trigobjphi",trigobjphi,"trigobjphi[ntrigobjs]/F");
  T1->Branch("trigobjmass",trigobjmass,"trigobjmass[ntrigobjs]/F");
  T1->Branch("trigobjHLT",trigobjHLT,"trigobjHLT[ntrigobjs]/O");
  T1->Branch("trigobjL1",trigobjL1,"trigobjL1[ntrigobjs]/O");
  //T1->Branch("trigobjBoth",trigobjBoth,"trigobjBoth[ntrigobjs]/O");
  T1->Branch("trigobjIhlt",trigobjIhlt,"trigobjIhlt[ntrigobjs]/I");
  T1->Branch("trigobjpdgId",trigobjpdgId,"trigobjpdgId[ntrigobjs]/I");
  
  // MET info //
  
  T1->Branch("CHS_MET",&miset,"miset/F") ;
  T1->Branch("CHS_METPhi",&misphi,"misphi/F") ;
  T1->Branch("CHS_METSig",&misetsig,"misetsig/F");
  T1->Branch("CHS_sumEt",&sumEt,"sumEt/F");
  
  T1->Branch("PUPPI_MET",&miset_PUPPI,"miset_PUPPI/F") ;
  T1->Branch("PUPPI_METPhi",&misphi_PUPPI,"misphi_PUPPI/F") ;
  T1->Branch("PUPPI_METSig",&misetsig_PUPPI,"misetsig_PUPPI/F");
  T1->Branch("PUPPI_sumEt",&sumEt_PUPPI,"sumEt_PUPPI/F");
  
  // AK8 jet info //
  
  T1->Branch("npfjetAK8",&npfjetAK8, "npfjetAK8/I"); 
  T1->Branch("pfjetAK8pt",pfjetAK8pt,"pfjetAK8pt[npfjetAK8]/F");
  T1->Branch("pfjetAK8y",pfjetAK8y,"pfjetAK8y[npfjetAK8]/F");
  T1->Branch("pfjetAK8eta",pfjetAK8eta,"pfjetAK8eta[npfjetAK8]/F");
  T1->Branch("pfjetAK8phi",pfjetAK8phi,"pfjetAK8phi[npfjetAK8]/F");
  T1->Branch("pfjetAK8mass",pfjetAK8mass,"pfjetAK8mass[npfjetAK8]/F");
  T1->Branch("pfjetAK8jetID_tightlepveto",pfjetAK8jetID_tightlepveto,"pfjetAK8jetID_tightlepveto[npfjetAK8]/O");
  T1->Branch("pfjetAK8jetID",pfjetAK8jetID,"pfjetAK8jetID[npfjetAK8]/O");
  T1->Branch("pfjetAK8JEC",pfjetAK8JEC,"pfjetAK8JEC[npfjetAK8]/F");
  T1->Branch("pfjetAK8CHF",pfjetAK8CHF,"pfjetAK8CHF[npfjetAK8]/F");
  T1->Branch("pfjetAK8NHF",pfjetAK8NHF,"pfjetAK8NHF[npfjetAK8]/F");
  T1->Branch("pfjetAK8CEMF",pfjetAK8CEMF,"pfjetAK8CEMF[npfjetAK8]/F");
  T1->Branch("pfjetAK8NEMF",pfjetAK8NEMF,"pfjetAK8NEMF[npfjetAK8]/F");
  T1->Branch("pfjetAK8MUF",pfjetAK8MUF,"pfjetAK8MUF[npfjetAK8]/F");
  T1->Branch("pfjetAK8PHF",pfjetAK8PHF,"pfjetAK8PHF[npfjetAK8]/F");
  T1->Branch("pfjetAK8EEF",pfjetAK8EEF,"pfjetAK8EEF[npfjetAK8]/F");
  T1->Branch("pfjetAK8HFHF",pfjetAK8HFHF,"pfjetAK8HFHF[npfjetAK8]/F");
  T1->Branch("pfjetAK8CHM",pfjetAK8CHM,"pfjetAK8CHM[npfjetAK8]/I");
  T1->Branch("pfjetAK8NHM",pfjetAK8NHM,"pfjetAK8NHM[npfjetAK8]/I");
  T1->Branch("pfjetAK8MUM",pfjetAK8MUM,"pfjetAK8MUM[npfjetAK8]/I");
  T1->Branch("pfjetAK8PHM",pfjetAK8PHM,"pfjetAK8PHM[npfjetAK8]/I");
  T1->Branch("pfjetAK8EEM",pfjetAK8EEM,"pfjetAK8EEM[npfjetAK8]/I");
  T1->Branch("pfjetAK8HFHM",pfjetAK8HFHM,"pfjetAK8HFHM[npfjetAK8]/I");
  T1->Branch("pfjetAK8Neucons",pfjetAK8Neucons,"pfjetAK8Neucons[npfjetAK8]/I");
  T1->Branch("pfjetAK8Chcons",pfjetAK8Chcons,"pfjetAK8Chcons[npfjetAK8]/I");
  
  T1->Branch("pfjetAK8reso",pfjetAK8reso,"pfjetAK8reso[npfjetAK8]/F");
  T1->Branch("pfjetAK8resoup",pfjetAK8resoup,"pfjetAK8resoup[npfjetAK8]/F");
  T1->Branch("pfjetAK8resodn",pfjetAK8resodn,"pfjetAK8resodn[npfjetAK8]/F");
  
  T1->Branch("pfjetAK8sdmass",pfjetAK8sdmass,"pfjetAK8sdmass[npfjetAK8]/F");
  T1->Branch("pfjetAK8tau1",pfjetAK8tau1,"pfjetAK8tau1[npfjetAK8]/F");
  T1->Branch("pfjetAK8tau2",pfjetAK8tau2,"pfjetAK8tau2[npfjetAK8]/F");
  T1->Branch("pfjetAK8tau3",pfjetAK8tau3,"pfjetAK8tau3[npfjetAK8]/F");
  
  T1->Branch("pfjetAK8btag_DeepCSV",pfjetAK8btag_DeepCSV,"pfjetAK8btag_DeepCSV[npfjetAK8]/F");
  T1->Branch("pfjetAK8DeepTag_DAK8_TvsQCD",pfjetAK8DeepTag_DAK8_TvsQCD,"pfjetAK8DeepTag_DAK8_TvsQCD[npfjetAK8]/F");
  T1->Branch("pfjetAK8DeepTag_DAK8_WvsQCD",pfjetAK8DeepTag_DAK8_WvsQCD,"pfjetAK8DeepTag_DAK8_WvsQCD[npfjetAK8]/F");
  T1->Branch("pfjetAK8DeepTag_DAK8_ZvsQCD",pfjetAK8DeepTag_DAK8_ZvsQCD,"pfjetAK8DeepTag_DAK8_ZvsQCD[npfjetAK8]/F");
  T1->Branch("pfjetAK8DeepTag_DAK8_HvsQCD",pfjetAK8DeepTag_DAK8_HvsQCD,"pfjetAK8DeepTag_DAK8_HvsQCD[npfjetAK8]/F");
  T1->Branch("pfjetAK8DeepTag_DAK8_bbvsQCD",pfjetAK8DeepTag_DAK8_bbvsQCD,"pfjetAK8DeepTag_DAK8_bbvsQCD[npfjetAK8]/F");
  T1->Branch("pfjetAK8DeepTag_PNet_TvsQCD",pfjetAK8DeepTag_PNet_TvsQCD,"pfjetAK8DeepTag_PNet_TvsQCD[npfjetAK8]/F");
  T1->Branch("pfjetAK8DeepTag_PNet_WvsQCD",pfjetAK8DeepTag_PNet_WvsQCD,"pfjetAK8DeepTag_PNet_WvsQCD[npfjetAK8]/F");
  T1->Branch("pfjetAK8DeepTag_PNet_ZvsQCD",pfjetAK8DeepTag_PNet_ZvsQCD,"pfjetAK8DeepTag_PNet_ZvsQCD[npfjetAK8]/F");
  T1->Branch("pfjetAK8DeepTag_PNet_XbbvsQCD",pfjetAK8DeepTag_PNet_XbbvsQCD,"pfjetAK8DeepTag_PNet_XbbvsQCD[npfjetAK8]/F");
  T1->Branch("pfjetAK8DeepTag_PNet_XccvsQCD",pfjetAK8DeepTag_PNet_XccvsQCD,"pfjetAK8DeepTag_PNet_XccvsQCD[npfjetAK8]/F");
  T1->Branch("pfjetAK8DeepTag_PNet_XqqvsQCD",pfjetAK8DeepTag_PNet_XqqvsQCD,"pfjetAK8DeepTag_PNet_XqqvsQCD[npfjetAK8]/F");
  
  T1->Branch("pfjetAK8sub1pt",pfjetAK8sub1pt,"pfjetAK8sub1pt[npfjetAK8]/F");
  T1->Branch("pfjetAK8sub1eta",pfjetAK8sub1eta,"pfjetAK8sub1eta[npfjetAK8]/F");
  T1->Branch("pfjetAK8sub1phi",pfjetAK8sub1phi,"pfjetAK8sub1phi[npfjetAK8]/F");
  T1->Branch("pfjetAK8sub1mass",pfjetAK8sub1mass,"pfjetAK8sub1mass[npfjetAK8]/F");
  T1->Branch("pfjetAK8sub1btag",pfjetAK8sub1btag,"pfjetAK8sub1btag[npfjetAK8]/F");
  
  T1->Branch("pfjetAK8sub2pt",pfjetAK8sub2pt,"pfjetAK8sub2pt[npfjetAK8]/F");
  T1->Branch("pfjetAK8sub2eta",pfjetAK8sub2eta,"pfjetAK8sub2eta[npfjetAK8]/F");
  T1->Branch("pfjetAK8sub2phi",pfjetAK8sub2phi,"pfjetAK8sub2phi[npfjetAK8]/F");
  T1->Branch("pfjetAK8sub2mass",pfjetAK8sub2mass,"pfjetAK8sub2mass[npfjetAK8]/F");
  T1->Branch("pfjetAK8sub2btag",pfjetAK8sub2btag,"pfjetAK8sub2btag[npfjetAK8]/F");
  
  T1->Branch("pfjetAK8jesup_AbsoluteStat",pfjetAK8jesup_AbsoluteStat,"pfjetAK8jesup_AbsoluteStat[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_AbsoluteScale",pfjetAK8jesup_AbsoluteScale,"pfjetAK8jesup_AbsoluteScale[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_AbsoluteMPFBias",pfjetAK8jesup_AbsoluteMPFBias,"pfjetAK8jesup_AbsoluteMPFBias[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_FlavorQCD",pfjetAK8jesup_FlavorQCD,"pfjetAK8jesup_FlavorQCD[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_Fragmentation",pfjetAK8jesup_Fragmentation,"pfjetAK8jesup_Fragmentation[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_PileUpDataMC",pfjetAK8jesup_PileUpDataMC,"pfjetAK8jesup_PileUpDataMC[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_PileUpPtBB",pfjetAK8jesup_PileUpPtBB,"pfjetAK8jesup_PileUpPtBB[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_PileUpPtEC1",pfjetAK8jesup_PileUpPtEC1,"pfjetAK8jesup_PileUpPtEC1[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_PileUpPtEC2",pfjetAK8jesup_PileUpPtEC2,"pfjetAK8jesup_PileUpPtEC2[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_PileUpPtRef",pfjetAK8jesup_PileUpPtRef,"pfjetAK8jesup_PileUpPtRef[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_RelativeFSR",pfjetAK8jesup_RelativeFSR,"pfjetAK8jesup_RelativeFSR[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_RelativeJEREC1",pfjetAK8jesup_RelativeJEREC1,"pfjetAK8jesup_RelativeJEREC1[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_RelativeJEREC2",pfjetAK8jesup_RelativeJEREC2,"pfjetAK8jesup_RelativeJEREC2[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_RelativePtBB",pfjetAK8jesup_RelativePtBB,"pfjetAK8jesup_RelativePtBB[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_RelativePtEC1",pfjetAK8jesup_RelativePtEC1,"pfjetAK8jesup_RelativePtEC1[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_RelativeJEREC2",pfjetAK8jesup_RelativeJEREC2,"pfjetAK8jesup_RelativeJEREC2[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_RelativePtBB",pfjetAK8jesup_RelativePtBB,"pfjetAK8jesup_RelativePtBB[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_RelativePtEC1",pfjetAK8jesup_RelativePtEC1,"pfjetAK8jesup_RelativePtEC1[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_RelativePtEC2",pfjetAK8jesup_RelativePtEC2,"pfjetAK8jesup_RelativePtEC2[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_RelativeBal",pfjetAK8jesup_RelativeBal,"pfjetAK8jesup_RelativeBal[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_RelativeSample",pfjetAK8jesup_RelativeSample,"pfjetAK8jesup_RelativeSample[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_RelativeStatEC",pfjetAK8jesup_RelativeStatEC,"pfjetAK8jesup_RelativeStatEC[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_RelativeStatFSR",pfjetAK8jesup_RelativeStatFSR,"pfjetAK8jesup_RelativeStatFSR[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_SinglePionECAL",pfjetAK8jesup_SinglePionECAL,"pfjetAK8jesup_SinglePionECAL[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_SinglePionHCAL",pfjetAK8jesup_SinglePionHCAL,"pfjetAK8jesup_SinglePionHCAL[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_TimePtEta",pfjetAK8jesup_TimePtEta,"pfjetAK8jesup_TimePtEta[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesup_Total",pfjetAK8jesup_Total,"pfjetAK8jesup_Total[npfjetAK8]/F");
  
  T1->Branch("pfjetAK8jesdn_AbsoluteStat",pfjetAK8jesdn_AbsoluteStat,"pfjetAK8jesdn_AbsoluteStat[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_AbsoluteScale",pfjetAK8jesdn_AbsoluteScale,"pfjetAK8jesdn_AbsoluteScale[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_AbsoluteMPFBias",pfjetAK8jesdn_AbsoluteMPFBias,"pfjetAK8jesdn_AbsoluteMPFBias[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_FlavorQCD",pfjetAK8jesdn_FlavorQCD,"pfjetAK8jesdn_FlavorQCD[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_Fragmentation",pfjetAK8jesdn_Fragmentation,"pfjetAK8jesdn_Fragmentation[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_PileUpDataMC",pfjetAK8jesdn_PileUpDataMC,"pfjetAK8jesdn_PileUpDataMC[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_PileUpPtBB",pfjetAK8jesdn_PileUpPtBB,"pfjetAK8jesdn_PileUpPtBB[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_PileUpPtEC1",pfjetAK8jesdn_PileUpPtEC1,"pfjetAK8jesdn_PileUpPtEC1[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_PileUpPtEC2",pfjetAK8jesdn_PileUpPtEC2,"pfjetAK8jesdn_PileUpPtEC2[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_PileUpPtRef",pfjetAK8jesdn_PileUpPtRef,"pfjetAK8jesdn_PileUpPtRef[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_RelativeFSR",pfjetAK8jesdn_RelativeFSR,"pfjetAK8jesdn_RelativeFSR[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_RelativeJEREC1",pfjetAK8jesdn_RelativeJEREC1,"pfjetAK8jesdn_RelativeJEREC1[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_RelativeJEREC2",pfjetAK8jesdn_RelativeJEREC2,"pfjetAK8jesdn_RelativeJEREC2[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_RelativePtBB",pfjetAK8jesdn_RelativePtBB,"pfjetAK8jesdn_RelativePtBB[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_RelativePtEC1",pfjetAK8jesdn_RelativePtEC1,"pfjetAK8jesdn_RelativePtEC1[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_RelativeJEREC2",pfjetAK8jesdn_RelativeJEREC2,"pfjetAK8jesdn_RelativeJEREC2[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_RelativePtBB",pfjetAK8jesdn_RelativePtBB,"pfjetAK8jesdn_RelativePtBB[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_RelativePtEC1",pfjetAK8jesdn_RelativePtEC1,"pfjetAK8jesdn_RelativePtEC1[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_RelativePtEC2",pfjetAK8jesdn_RelativePtEC2,"pfjetAK8jesdn_RelativePtEC2[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_RelativeBal",pfjetAK8jesdn_RelativeBal,"pfjetAK8jesdn_RelativeBal[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_RelativeSample",pfjetAK8jesdn_RelativeSample,"pfjetAK8jesdn_RelativeSample[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_RelativeStatEC",pfjetAK8jesdn_RelativeStatEC,"pfjetAK8jesdn_RelativeStatEC[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_RelativeStatFSR",pfjetAK8jesdn_RelativeStatFSR,"pfjetAK8jesdn_RelativeStatFSR[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_SinglePionECAL",pfjetAK8jesdn_SinglePionECAL,"pfjetAK8jesdn_SinglePionECAL[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_SinglePionHCAL",pfjetAK8jesdn_SinglePionHCAL,"pfjetAK8jesdn_SinglePionHCAL[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_TimePtEta",pfjetAK8jesdn_TimePtEta,"pfjetAK8jesdn_TimePtEta[npfjetAK8]/F");
  T1->Branch("pfjetAK8jesdn_Total",pfjetAK8jesdn_Total,"pfjetAK8jesdn_Total[npfjetAK8]/F");
  
  // AK4 jet info //
 
  T1->Branch("npfjetAK4",&npfjetAK4,"npfjetAK4/I"); 

  T1->Branch("pfjetAK4jetID",pfjetAK4jetID,"pfjetAK4jetID[npfjetAK4]/O");
  T1->Branch("pfjetAK4jetID_tightlepveto",pfjetAK4jetID_tightlepveto,"pfjetAK4jetID_tightlepveto[npfjetAK4]/O");
  
  T1->Branch("pfjetAK4pt",pfjetAK4pt,"pfjetAK4pt[npfjetAK4]/F");
  T1->Branch("pfjetAK4eta",pfjetAK4eta,"pfjetAK4eta[npfjetAK4]/F");
  T1->Branch("pfjetAK4y",pfjetAK4y,"pfjetAK4y[npfjetAK4]/F");
  T1->Branch("pfjetAK4phi",pfjetAK4phi,"pfjetAK4phi[npfjetAK4]/F");
  T1->Branch("pfjetAK4mass",pfjetAK4mass,"pfjetAK4mass[npfjetAK4]/F");
  T1->Branch("pfjetAK4JEC",pfjetAK4JEC,"pfjetAK4JEC[npfjetAK4]/F");
  T1->Branch("pfjetAK4btag_DeepCSV",pfjetAK4btag_DeepCSV,"pfjetAK4btag_DeepCSV[npfjetAK4]/F");
  T1->Branch("pfjetAK4btag_DeepFlav",pfjetAK4btag_DeepFlav,"pfjetAK4btag_DeepFlav[npfjetAK4]/F");
 
  T1->Branch("pfjetAK4reso",pfjetAK4reso,"pfjetAK4reso[npfjetAK4]/F");
  T1->Branch("pfjetAK4resoup",pfjetAK4resoup,"pfjetAK4resoup[npfjetAK4]/F");
  T1->Branch("pfjetAK4resodn",pfjetAK4resodn,"pfjetAK4resodn[npfjetAK4]/F"); 
  
  T1->Branch("pfjetAK4hadronflav",pfjetAK4hadronflav,"pfjetAK4hadronflav[npfjetAK4]/I");
  T1->Branch("pfjetAK4partonflav",pfjetAK4partonflav,"pfjetAK4partonflav[npfjetAK4]/I");
  T1->Branch("pfjetAK4qgl",pfjetAK4qgl,"pfjetAK4qgl[npfjetAK4]/F");
  T1->Branch("pfjetAK4PUID",pfjetAK4PUID,"pfjetAK4PUID[npfjetAK4]/F");
  
  T1->Branch("pfjetAK4jesup_AbsoluteStat",pfjetAK4jesup_AbsoluteStat,"pfjetAK4jesup_AbsoluteStat[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_AbsoluteScale",pfjetAK4jesup_AbsoluteScale,"pfjetAK4jesup_AbsoluteScale[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_AbsoluteMPFBias",pfjetAK4jesup_AbsoluteMPFBias,"pfjetAK4jesup_AbsoluteMPFBias[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_FlavorQCD",pfjetAK4jesup_FlavorQCD,"pfjetAK4jesup_FlavorQCD[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_Fragmentation",pfjetAK4jesup_Fragmentation,"pfjetAK4jesup_Fragmentation[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_PileUpDataMC",pfjetAK4jesup_PileUpDataMC,"pfjetAK4jesup_PileUpDataMC[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_PileUpPtBB",pfjetAK4jesup_PileUpPtBB,"pfjetAK4jesup_PileUpPtBB[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_PileUpPtEC1",pfjetAK4jesup_PileUpPtEC1,"pfjetAK4jesup_PileUpPtEC1[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_PileUpPtEC2",pfjetAK4jesup_PileUpPtEC2,"pfjetAK4jesup_PileUpPtEC2[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_PileUpPtRef",pfjetAK4jesup_PileUpPtRef,"pfjetAK4jesup_PileUpPtRef[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_RelativeFSR",pfjetAK4jesup_RelativeFSR,"pfjetAK4jesup_RelativeFSR[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_RelativeJEREC1",pfjetAK4jesup_RelativeJEREC1,"pfjetAK4jesup_RelativeJEREC1[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_RelativeJEREC2",pfjetAK4jesup_RelativeJEREC2,"pfjetAK4jesup_RelativeJEREC2[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_RelativePtBB",pfjetAK4jesup_RelativePtBB,"pfjetAK4jesup_RelativePtBB[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_RelativePtEC1",pfjetAK4jesup_RelativePtEC1,"pfjetAK4jesup_RelativePtEC1[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_RelativeJEREC2",pfjetAK4jesup_RelativeJEREC2,"pfjetAK4jesup_RelativeJEREC2[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_RelativePtBB",pfjetAK4jesup_RelativePtBB,"pfjetAK4jesup_RelativePtBB[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_RelativePtEC1",pfjetAK4jesup_RelativePtEC1,"pfjetAK4jesup_RelativePtEC1[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_RelativePtEC2",pfjetAK4jesup_RelativePtEC2,"pfjetAK4jesup_RelativePtEC2[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_RelativeBal",pfjetAK4jesup_RelativeBal,"pfjetAK4jesup_RelativeBal[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_RelativeSample",pfjetAK4jesup_RelativeSample,"pfjetAK4jesup_RelativeSample[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_RelativeStatEC",pfjetAK4jesup_RelativeStatEC,"pfjetAK4jesup_RelativeStatEC[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_RelativeStatFSR",pfjetAK4jesup_RelativeStatFSR,"pfjetAK4jesup_RelativeStatFSR[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_SinglePionECAL",pfjetAK4jesup_SinglePionECAL,"pfjetAK4jesup_SinglePionECAL[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_SinglePionHCAL",pfjetAK4jesup_SinglePionHCAL,"pfjetAK4jesup_SinglePionHCAL[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_TimePtEta",pfjetAK4jesup_TimePtEta,"pfjetAK4jesup_TimePtEta[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesup_Total",pfjetAK4jesup_Total,"pfjetAK4jesup_Total[npfjetAK4]/F");
  
  T1->Branch("pfjetAK4jesdn_AbsoluteStat",pfjetAK4jesdn_AbsoluteStat,"pfjetAK4jesdn_AbsoluteStat[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_AbsoluteScale",pfjetAK4jesdn_AbsoluteScale,"pfjetAK4jesdn_AbsoluteScale[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_AbsoluteMPFBias",pfjetAK4jesdn_AbsoluteMPFBias,"pfjetAK4jesdn_AbsoluteMPFBias[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_FlavorQCD",pfjetAK4jesdn_FlavorQCD,"pfjetAK4jesdn_FlavorQCD[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_Fragmentation",pfjetAK4jesdn_Fragmentation,"pfjetAK4jesdn_Fragmentation[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_PileUpDataMC",pfjetAK4jesdn_PileUpDataMC,"pfjetAK4jesdn_PileUpDataMC[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_PileUpPtBB",pfjetAK4jesdn_PileUpPtBB,"pfjetAK4jesdn_PileUpPtBB[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_PileUpPtEC1",pfjetAK4jesdn_PileUpPtEC1,"pfjetAK4jesdn_PileUpPtEC1[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_PileUpPtEC2",pfjetAK4jesdn_PileUpPtEC2,"pfjetAK4jesdn_PileUpPtEC2[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_PileUpPtRef",pfjetAK4jesdn_PileUpPtRef,"pfjetAK4jesdn_PileUpPtRef[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_RelativeFSR",pfjetAK4jesdn_RelativeFSR,"pfjetAK4jesdn_RelativeFSR[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_RelativeJEREC1",pfjetAK4jesdn_RelativeJEREC1,"pfjetAK4jesdn_RelativeJEREC1[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_RelativeJEREC2",pfjetAK4jesdn_RelativeJEREC2,"pfjetAK4jesdn_RelativeJEREC2[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_RelativePtBB",pfjetAK4jesdn_RelativePtBB,"pfjetAK4jesdn_RelativePtBB[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_RelativePtEC1",pfjetAK4jesdn_RelativePtEC1,"pfjetAK4jesdn_RelativePtEC1[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_RelativeJEREC2",pfjetAK4jesdn_RelativeJEREC2,"pfjetAK4jesdn_RelativeJEREC2[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_RelativePtBB",pfjetAK4jesdn_RelativePtBB,"pfjetAK4jesdn_RelativePtBB[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_RelativePtEC1",pfjetAK4jesdn_RelativePtEC1,"pfjetAK4jesdn_RelativePtEC1[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_RelativePtEC2",pfjetAK4jesdn_RelativePtEC2,"pfjetAK4jesdn_RelativePtEC2[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_RelativeBal",pfjetAK4jesdn_RelativeBal,"pfjetAK4jesdn_RelativeBal[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_RelativeSample",pfjetAK4jesdn_RelativeSample,"pfjetAK4jesdn_RelativeSample[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_RelativeStatEC",pfjetAK4jesdn_RelativeStatEC,"pfjetAK4jesdn_RelativeStatEC[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_RelativeStatFSR",pfjetAK4jesdn_RelativeStatFSR,"pfjetAK4jesdn_RelativeStatFSR[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_SinglePionECAL",pfjetAK4jesdn_SinglePionECAL,"pfjetAK4jesdn_SinglePionECAL[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_SinglePionHCAL",pfjetAK4jesdn_SinglePionHCAL,"pfjetAK4jesdn_SinglePionHCAL[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_TimePtEta",pfjetAK4jesdn_TimePtEta,"pfjetAK4jesdn_TimePtEta[npfjetAK4]/F");
  T1->Branch("pfjetAK4jesdn_Total",pfjetAK4jesdn_Total,"pfjetAK4jesdn_Total[npfjetAK4]/F");
  
  // Muon info //
  
  T1->Branch("nmuons",&nmuons,"nmuons/I");
  T1->Branch("muonisPF",muonisPF,"muonisPF[nmuons]/O");
  T1->Branch("muonisGL",muonisGL,"muonisGL[nmuons]/O");
  T1->Branch("muonisTRK",muonisTRK,"muonisTRK[nmuons]/O");
  T1->Branch("muonisLoose",muonisLoose,"muonisLoose[nmuons]/O");
  T1->Branch("muonisGoodGL",muonisGoodGL,"muonisGoodGL[nmuons]/O");
  T1->Branch("muonisMed",muonisMed,"muonisMed[nmuons]/O");
  T1->Branch("muonisMedPr",muonisMedPr,"muonisMedPr[nmuons]/O");
  T1->Branch("muonisTight",muonisTight,"muonisTight[nmuons]/O");
  T1->Branch("muonisHighPt",muonisHighPt,"muonisHighPt[nmuons]/O"); 
  T1->Branch("muonisHighPttrk",muonisHighPttrk,"muonisHighPttrk[nmuons]/O");
  
  T1->Branch("muonpt",muonpt,"muonpt[nmuons]/F");
  T1->Branch("muonp",muonp,"muonp[nmuons]/F");
  T1->Branch("muoneta",muoneta,"muoneta[nmuons]/F");
  T1->Branch("muonphi",muonphi,"muonphi[nmuons]/F");
  
  T1->Branch("muonminchiso", muonminchiso, "muonminchiso[nmuons]/F");
  T1->Branch("muonminnhiso", muonminnhiso, "muonminnhiso[nmuons]/F");
  T1->Branch("muonminphiso", muonminphiso, "muonminphiso[nmuons]/F");
  T1->Branch("muonminisoall", muonminisoall, "muonminisoall[nmuons]/F");
  T1->Branch("muontrkvtx",muontrkvtx,"muontrkvtx[nmuons]/F");
  T1->Branch("muondz",muondz,"muondz[nmuons]/F");
  T1->Branch("muonpter",muonpter,"muonpter[nmuons]/F");
  T1->Branch("muonchi",muonchi,"muonchi[nmuons]/F");
  T1->Branch("muonndf",muonndf,"muonndf[nmuons]/I");
  T1->Branch("muonecal",muonecal,"muonecal[nmuons]/F");
  T1->Branch("muonhcal",muonhcal,"muonhcal[nmuons]/F");
  T1->Branch("muonpfiso",muonpfiso,"muonpfiso[nmuons]/F");
  T1->Branch("muonposmatch",muonposmatch,"muonposmatch[nmuons]/F");
  T1->Branch("muontrkink",muontrkink,"muontrkink[nmuons]/F");
  T1->Branch("muonsegcom",muonsegcom,"muonsegcom[nmuons]/F");
  T1->Branch("muonhit",muonhit,"muonhit[nmuons]/F");
  T1->Branch("muonpixhit",muonpixhit,"muonpixhit[nmuons]/F");
  T1->Branch("muonmst",muonmst,"muonmst[nmuons]/F");
  T1->Branch("muontrklay",muontrklay,"muontrklay[nmuons]/F"); 
  T1->Branch("muonvalfrac",muonvalfrac,"muonvalfrac[nmuons]/F"); 
  T1->Branch("mudxy_sv",mudxy_sv,"mudxy_sv[nmuons]/F");
  
  // Electron info //
  
  T1->Branch("nelecs",&nelecs,"nelecs/I");
  T1->Branch("elpt",elpt,"elpt[nelecs]/F");
  T1->Branch("eleta",eleta,"eleta[nelecs]/F");
  T1->Branch("elphi",elphi,"elphi[nelecs]/F");
  T1->Branch("elp",elp,"elp[nelecs]/F");
  T1->Branch("ele",ele,"ele[nelecs]/F");
  
  T1->Branch("elmvaid_Fallv2WP90",elmvaid_Fallv2WP90,"elmvaid_Fallv2WP90[nelecs]/O");
  T1->Branch("elmvaid_Fallv2WP90_noIso",elmvaid_Fallv2WP90_noIso,"elmvaid_Fallv2WP90_noIso[nelecs]/O");
  T1->Branch("elmvaid_Fallv2WP80",elmvaid_Fallv2WP80,"elmvaid_Fallv2WP80[nelecs]/O");
  T1->Branch("elmvaid_Fallv2WP80_noIso",elmvaid_Fallv2WP80_noIso,"elmvaid_Fallv2WP80_noIso[nelecs]/O");
  
  T1->Branch("eldxytrk",eldxytrk,"eldxytrk[nelecs]/F");
  T1->Branch("eldxy_sv",eldxy_sv,"eldxy_sv[nelecs]/F");
  T1->Branch("eldztrk",eldztrk,"eldztrk[nelecs]/F");
  T1->Branch("elhovere",elhovere,"elhovere[nelecs]/F");
  T1->Branch("elchi",elchi,"elchi[nelecs]/F");
  T1->Branch("elndf",elndf,"elndf[nelecs]/I");
  T1->Branch("eletain",eletain,"eletain[nelecs]/F");
  T1->Branch("elphiin",elphiin,"elphiin[nelecs]/F");
  T1->Branch("elfbrem",elfbrem,"elfbrem[nelecs]/F");
  T1->Branch("eleoverp",eleoverp,"eleoverp[nelecs]/F");
  T1->Branch("elietaieta",elietaieta,"elietaieta[nelecs]/F");
  T1->Branch("elmisshits",elmisshits,"elmisshits[nelecs]/F");
  T1->Branch("elpfiso_drcor",elpfiso_drcor,"elpfiso_drcor[nelecs]/F");
  T1->Branch("elpfiso_eacor",elpfiso_eacor,"elpfiso_eacor[nelecs]/F");
  T1->Branch("elpfiso04_eacor",elpfiso04_eacor,"elpfiso04_eacor[nelecs]/F");
  
  T1->Branch("elsupcl_eta",elsupcl_eta,"elsupcl_eta[nelecs]/F");
  T1->Branch("elsupcl_phi",elsupcl_phi,"elsupcl_phi[nelecs]/F");
  T1->Branch("elsupcl_rawE",elsupcl_rawE,"elsupcl_rawE[nelecs]/F");
  T1->Branch("elsigmaieta", elsigmaieta, "elsigmaieta[nelecs]/F");
  T1->Branch("elsigmaiphi", elsigmaiphi, "elsigmaiphi[nelecs]/F");
  T1->Branch("elr9full", elr9full, "elr9full[nelecs]/F");
  T1->Branch("elsupcl_etaw", elsupcl_etaw, "elsupcl_etaw[nelecs]/F");
  T1->Branch("elsupcl_phiw", elsupcl_phiw, "elsupcl_phiw[nelecs]/F");
  T1->Branch("elhcaloverecal", elhcaloverecal, "elhcaloverecal[nelecs]/F");
  T1->Branch("elcloctftrkn", elcloctftrkn, "elcloctftrkn[nelecs]/F");
  T1->Branch("elcloctftrkchi2", elcloctftrkchi2, "elcloctftrkchi2[nelecs]/F");
  T1->Branch("ele1x5bye5x5", ele1x5bye5x5, "ele1x5bye5x5[nelecs]/F");
  T1->Branch("elnormchi2", elnormchi2, "elnormchi2[nelecs]/F");
  T1->Branch("elhitsmiss", elhitsmiss, "elhitsmiss[nelecs]/F");
  T1->Branch("eltrkmeasure", eltrkmeasure, "eltrkmeasure[nelecs]/F");
  T1->Branch("elconvtxprob", elconvtxprob, "elconvtxprob[nelecs]/F");
  T1->Branch("elecloverpout", elecloverpout, "elecloverpout[nelecs]/F");
  T1->Branch("elecaletrkmomentum", elecaletrkmomentum, "elecaletrkmomentum[nelecs]/F");
  T1->Branch("eldeltaetacltrkcalo", eldeltaetacltrkcalo, "eldeltaetacltrkcalo[nelecs]/F");
  T1->Branch("elsupcl_preshvsrawe", elsupcl_preshvsrawe, "elsupcl_preshvsrawe[nelecs]/F");
  T1->Branch("elpfisolsumphet", elpfisolsumphet, "elpfisolsumphet[nelecs]/F");
  T1->Branch("elpfisolsumchhadpt", elpfisolsumchhadpt, "elpfisolsumchhadpt[nelecs]/F");
  T1->Branch("elpfsiolsumneuhadet", elpfsiolsumneuhadet, "elpfsiolsumneuhadet[nelecs]/F");
  
  T1->Branch("elminchiso", elminchiso, "elminchiso[nelecs]/F");
  T1->Branch("elminnhiso", elminnhiso, "elminnhiso[nelecs]/F");
  T1->Branch("elminphiso", elminphiso, "elminphiso[nelecs]/F");
  T1->Branch("elminisoall", elminisoall, "elminisoall[nelecs]/F");
  
  // Photon Info //
  
  T1->Branch("nphotons",&nphotons,"nphotons/I");
  T1->Branch("phoe",phoe,"phoe[nphotons]/F");
  T1->Branch("phoeta",phoeta,"phoeta[nphotons]/F");
  T1->Branch("phophi",phophi,"phophi[nphotons]/F");
  T1->Branch("phomvaid_Fall17V2_raw",phomvaid_Fall17V2_raw,"phomvaid_Fall17V2_raw[nphotons]/F");
  T1->Branch("phomvaid_Fall17V2_WP90",phomvaid_Fall17V2_WP90,"phomvaid_Fall17V2_WP90[nphotons]/O");
  T1->Branch("phomvaid_Fall17V2_WP80",phomvaid_Fall17V2_WP80,"phomvaid_Fall17V2_WP80[nphotons]/O");
  T1->Branch("phomvaid_Spring16V1_WP90",phomvaid_Spring16V1_WP90,"phomvaid_Spring16V1_WP90[nphotons]/O");
  T1->Branch("phomvaid_Spring16V1_WP80",phomvaid_Spring16V1_WP80,"phomvaid_Spring16V1_WP80[nphotons]/O");
  T1->Branch("phoe1by9",phoe1by9,"phoe1by9[nphotons]/F");
  T1->Branch("phoe9by25",phoe9by25,"phoe9by25[nphotons]/F");
  T1->Branch("photrkiso",photrkiso,"photrkiso[nphotons]/F");
  T1->Branch("phoemiso",phoemiso,"phoemiso[nphotons]/F");
  T1->Branch("phohadiso",phohadiso,"phohadiso[nphotons]/F");
  T1->Branch("phochhadiso",phochhadiso,"phochhadiso[nphotons]/F");
  T1->Branch("phoneuhadiso",phoneuhadiso,"phoneuhadiso[nphotons]/F");
  T1->Branch("phophoiso",phophoiso,"phophoiso[nphotons]/F");
  T1->Branch("phoPUiso",phoPUiso,"phoPUiso[nphotons]/F");
  T1->Branch("phohadbyem",phohadbyem,"phohadbyem[nphotons]/F");
  T1->Branch("phoietaieta",phoietaieta,"phoietaieta[nphotons]/F");

  // Tau Info //
  
  T1->Branch("ntaus",&ntaus,"ntaus/I");
  T1->Branch("tauisPF",tauisPF,"tauisPF[ntaus]/O");
  T1->Branch("taupt",taupt,"taupt[ntaus]/F");
  T1->Branch("taueta",taueta,"taueta[ntaus]/F");
  T1->Branch("tauphi",tauphi,"tauphi[ntaus]/F");
  T1->Branch("taue",taue,"taue[ntaus]/F");
  T1->Branch("taucharge",taucharge,"taucharge[ntaus]/I");
  T1->Branch("taudxy",taudxy,"taudxy[ntaus]/F");
  T1->Branch("tauleadtrkdxy",tauleadtrkdxy,"tauleadtrkdxy[ntaus]/F");
  T1->Branch("tauleadtrkdz",tauleadtrkdz,"tauleadtrkdz[ntaus]/F");
  T1->Branch("tauleadtrkpt",tauleadtrkpt,"tauleadtrkpt[ntaus]/F");
  T1->Branch("tauleadtrketa",tauleadtrketa,"tauleadtrketa[ntaus]/F");
  T1->Branch("tauleadtrkphi",tauleadtrkphi,"tauleadtrkphi[ntaus]/F");
  T1->Branch("taudecayMode",taudecayMode,"taudecayMode[ntaus]/I");
  T1->Branch("taudecayModeinding",taudecayModeinding,"taudecayModeinding[ntaus]/O");
  T1->Branch("taudecayModeindingNewDMs",taudecayModeindingNewDMs,"taudecayModeindingNewDMs[ntaus]/O");
  T1->Branch("taueiso2018_raw",taueiso2018_raw,"taueiso2018_raw[ntaus]/F");
  T1->Branch("taueiso2018",taueiso2018,"taueiso2018[ntaus]/I");
  T1->Branch("taujetiso_deeptau2017v2p1_raw",taujetiso_deeptau2017v2p1_raw,"taujetiso_deeptau2017v2p1_raw[ntaus]/F");
  T1->Branch("taujetiso_deeptau2017v2p1",taujetiso_deeptau2017v2p1,"taujetiso_deeptau2017v2p1[ntaus]/I");
  T1->Branch("taueiso_deeptau2017v2p1_raw",taueiso_deeptau2017v2p1_raw,"taueiso_deeptau2017v2p1_raw[ntaus]/F");
  T1->Branch("taueiso_deeptau2017v2p1",taueiso_deeptau2017v2p1,"taueiso_deeptau2017v2p1[ntaus]/I");
  T1->Branch("taumuiso_deeptau2017v2p1_raw",taumuiso_deeptau2017v2p1_raw,"taumuiso_deeptau2017v2p1_raw[ntaus]/F");
  T1->Branch("taumuiso_deeptau2017v2p1",taumuiso_deeptau2017v2p1,"taumuiso_deeptau2017v2p1[ntaus]/I");
  T1->Branch("taurawiso",taurawiso,"taurawiso[ntaus]/F");
  T1->Branch("taurawisodR03",taurawisodR03,"taurawisodR03[ntaus]/F");
  T1->Branch("taupuCorr",taupuCorr,"taupuCorr[ntaus]/F");
  
  // MC Info //
  
  if(isMC){
	  
  // generator-related info //
  
  T1->Branch("event_weight", &event_weight, "event_weight/D") ;
  T1->Branch("qscale",&qscale,"qscale/F");
  T1->Branch("Generator_x1",&Generator_x1,"Generator_x1/F");
  T1->Branch("Generator_x2",&Generator_x2,"Generator_x2/F");
  T1->Branch("Generator_xpdf1",&Generator_xpdf1,"Generator_xpdf1/F");
  T1->Branch("Generator_xpdf2",&Generator_xpdf2,"Generator_xpdf2/F");
  T1->Branch("Generator_id1",&Generator_id1,"Generator_id1/I");
  T1->Branch("Generator_id2",&Generator_id2,"Generator_id2/I");
  T1->Branch("Generator_scalePDF",&Generator_scalePDF,"Generator_scalePDF/F");
  
  T1->Branch("npu_vert",&npu_vert,"npu_vert/I");
  T1->Branch("npu_vert_true",&npu_vert_true,"npu_vert_true/I");
	  
  // GEN MET info //    
  
  T1->Branch("GENMET",&genmiset,"genmiset/F") ;
  T1->Branch("GENMETPhi",&genmisphi,"genmisphi/F") ;
  
  // GEN AK8 jet info //  
  
  T1->Branch("ngenjetAK8",&ngenjetAK8, "ngenjetAK8/I");
  T1->Branch("genjetAK8pt",genjetAK8pt,"genjetAK8pt[ngenjetAK8]/F");
  T1->Branch("genjetAK8eta",genjetAK8eta,"genjetAK8eta[ngenjetAK8]/F");
  T1->Branch("genjetAK8phi",genjetAK8phi,"genjetAK8phi[ngenjetAK8]/F");
  T1->Branch("genjetAK8mass",genjetAK8mass,"genjetAK8mass[ngenjetAK8]/F"); 
  T1->Branch("genjetAK8sdmass",genjetAK8sdmass,"genjetAK8sdmass[ngenjetAK8]/F");
  T1->Branch("genjetAK8hadronflav",genjetAK8hadronflav,"genjetAK8hadronflav[ngenjetAK8]/I");
  T1->Branch("genjetAK8partonflav",genjetAK8partonflav,"genjetAK8partonflav[ngenjetAK8]/I");

  // GEN AK4 jet info //  
 
  T1->Branch("ngenjetAK4",&ngenjetAK4, "ngenjetAK4/I");
  T1->Branch("genjetAK4pt",genjetAK4pt,"genjetAK4pt[ngenjetAK4]/F");
  T1->Branch("genjetAK4eta",genjetAK4eta,"genjetAK4eta[ngenjetAK4]/F");
  T1->Branch("genjetAK4phi",genjetAK4phi,"genjetAK4phi[ngenjetAK4]/F");
  T1->Branch("genjetAK4mass",genjetAK4mass,"genjetAK4mass[ngenjetAK4]/F");
  T1->Branch("genjetAK4hadronflav",genjetAK4hadronflav,"genjetAK4hadronflav[ngenjetAK4]/I");
  T1->Branch("genjetAK4partonflav",genjetAK4partonflav,"genjetAK4partonflav[ngenjetAK4]/I");
  
  // GEN particles info //  
  
  T1->Branch("ngenparticles",&ngenparticles, "ngenparticles/I");
  T1->Branch("genpartstatus",genpartstatus,"genpartstatus[ngenparticles]/I");
  T1->Branch("genpartpdg",genpartpdg,"genpartpdg[ngenparticles]/I");
  T1->Branch("genpartmompdg",genpartmompdg,"genpartmompdg[ngenparticles]/I");
  T1->Branch("genpartgrmompdg",genpartgrmompdg,"genpartgrmompdg[ngenparticles]/I");
  T1->Branch("genpartdaugno",genpartdaugno,"genpartdaugno[ngenparticles]/I");
  T1->Branch("genpartfromhard",genpartfromhard,"genpartfromhard[ngenparticles]/O");
  T1->Branch("genpartfromhardbFSR",genpartfromhardbFSR,"genpartfromhardbFSR[ngenparticles]/O");
  T1->Branch("genpartisPromptFinalState",genpartisPromptFinalState,"genpartisPromptFinalState[ngenparticles]/O");
  T1->Branch("genpartisLastCopyBeforeFSR",genpartisLastCopyBeforeFSR,"genpartisLastCopyBeforeFSR[ngenparticles]/O");
  T1->Branch("genpartpt",genpartpt,"genpartpt[ngenparticles]/F");
  T1->Branch("genparteta",genparteta,"genparteta[ngenparticles]/F");
  T1->Branch("genpartphi",genpartphi,"genpartphi[ngenparticles]/F");
  T1->Branch("genpartm",genpartm,"genpartm[ngenparticles]/F");
  
  // LHE Info //
  
  T1->Branch("nLHEparticles",&nLHEparticles, "nLHEparticles/I");
  T1->Branch("LHEpartpdg",LHEpartpdg,"LHEpartpdg[nLHEparticles]/I");
  T1->Branch("LHEpartpt",LHEpartpt,"LHEpartpt[nLHEparticles]/F");
  T1->Branch("LHEparteta",LHEparteta,"LHEparteta[nLHEparticles]/F");
  T1->Branch("LHEpartphi",LHEpartphi,"LHEpartphi[nLHEparticles]/F");
  T1->Branch("LHEpartm",LHEpartm,"LHEpartm[nLHEparticles]/F");
  
  T1->Branch("event_weight_LHE",&event_weight_LHE, "event_weight_LHE/D");
  T1->Branch("nLHEScaleWeights",&nLHEScaleWeights, "nLHEScaleWeights/I");
  T1->Branch("LHEScaleWeights",LHEScaleWeights,"LHEScaleWeights[nLHEScaleWeights]/F");
  T1->Branch("nLHEPDFWeights",&nLHEPDFWeights, "nLHEPDFWeights/I");
  T1->Branch("LHEPDFWeights",LHEPDFWeights,"LHEPDFWeights[nLHEPDFWeights]/F");
  T1->Branch("nLHEAlpsWeights",&nLHEAlpsWeights, "nLHEAlpsWeights/I");
  T1->Branch("LHEAlpsWeights",LHEAlpsWeights,"LHEAlpsWeights[nLHEAlpsWeights]/F");
  T1->Branch("nLHEPSWeights",&nLHEPSWeights, "nLHEPSWeights/I");
  T1->Branch("LHEPSWeights",LHEPSWeights,"LHEPSWeights[nLHEPSWeights]/F");
  
  } //isMC
  
  Nevt=0;
}


Leptop::~Leptop()
{
 
  // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Leptop::analyze(const edm::Event& iEvent, const edm::EventSetup& pset) {
  
  using namespace edm;
  Nevt++;
  
  irun = iEvent.id().run();
  ilumi = iEvent.luminosityBlock();
  ievt = iEvent.id().event();
  
  if (Nevt%100==1)cout <<"Leptop::analyze "<<Nevt<<" "<<iEvent.id().run()<<" "<<iEvent.id().event()<<endl;
  
  wtfact = 1.;
  
  // First store all MC information //
  
  edm::Handle<reco::GenJetCollection> genjetAK8s;
  edm::Handle<reco::GenJetCollection> genjetAK4s;
  
  iEvent.getByToken(tok_genjetAK4s_, genjetAK4s);
  
  if(isMC){
	
	// MC weights //
	  
    edm::Handle<GenEventInfoProduct>eventinfo ;  
    iEvent.getByToken(tok_wt_,eventinfo) ;
    
    if (eventinfo.isValid()){
		
		event_weight = eventinfo->weight();
		wtfact *= event_weight;
		
		qscale = eventinfo->qScale();
		Generator_x1 = (*eventinfo->pdf()).x.first;
        Generator_x2 = (*eventinfo->pdf()).x.second;
        Generator_id1 = (*eventinfo->pdf()).id.first;
        Generator_id2 = (*eventinfo->pdf()).id.second;
        Generator_xpdf1 = (*eventinfo->pdf()).xPDF.first;
        Generator_xpdf2 = (*eventinfo->pdf()).xPDF.second;
        Generator_scalePDF = (*eventinfo->pdf()).scalePDF;
        
        //cout<<"eventinfo->weights().size() "<<eventinfo->weights().size()<<" GEN weight "<<event_weight<<endl;
        // Parton shower weights //
        nLHEPSWeights = 0;
        if(eventinfo->weights().size()>2){
			for(unsigned int i=2; i<eventinfo->weights().size(); ++i){
				LHEPSWeights[nLHEPSWeights] = eventinfo->weights()[i]/eventinfo->weights()[1];
				nLHEPSWeights++;
				if(nLHEPSWeights >= nlhepsmax) break;
			}
		}

    }
    else
      {
		event_weight = qscale = wtfact = -10000;
      }
    
    // LHE-level particles //
      
    edm::Handle<LHEEventProduct>lheeventinfo ;
    iEvent.getByToken(lheEventProductToken_,lheeventinfo) ;
    
    if(lheeventinfo.isValid()){
      
      const auto & hepeup = lheeventinfo->hepeup();
      const auto & pup = hepeup.PUP;
      
      nLHEparticles = 0;
      
      for (unsigned int i = 0; i  < pup.size(); ++i) {
		if(hepeup.ISTUP[i]==1){// status==1 --> particles stay up to final state                          
			TLorentzVector p4(pup[i][0], pup[i][1], pup[i][2], pup[i][3]);
			LHEpartpt[nLHEparticles] = p4.Pt();
			LHEparteta[nLHEparticles] = p4.Eta();
			LHEpartphi[nLHEparticles] = p4.Phi();
			LHEpartm[nLHEparticles] = p4.M();
			LHEpartpdg[nLHEparticles] = (hepeup.IDUP[i]);
			nLHEparticles++;
			if(nLHEparticles>=nlhemax) break;
			}
		}
	
	 // LHE-level weights //
	
	  event_weight_LHE = lheeventinfo->originalXWGTUP();
	  
	  nLHEScaleWeights = 0;
	  nLHEPDFWeights = 0;
	  nLHEAlpsWeights = 0;
	  	  
	  for ( unsigned int index = 0; index < lheeventinfo->weights().size(); ++index ) {	
		//cout<<"Index "<<index+1<<" Id "<<lheeventinfo->weights()[index].id<<" weight "<<lheeventinfo->weights()[index].wgt/lheeventinfo->originalXWGTUP()<<endl;//" muR "<<lheeventinfo->weights()[index].MUR<<" muF "<<lheeventinfo->weights()[index].MUF<<" DYN Scale "<<lheeventinfo->weights()[index].DYN_SCALE<<endl;
		if(index<nlhescalemax && nLHEScaleWeights<nlhescalemax){
			LHEScaleWeights[nLHEScaleWeights] = lheeventinfo->weights()[index].wgt/lheeventinfo->originalXWGTUP();
			nLHEScaleWeights++;
		}
		if(index>=nlhescalemax && index<(nlhescalemax+nPDFsets)  && nLHEPDFWeights<nlhepdfmax){
			LHEPDFWeights[nLHEPDFWeights] = lheeventinfo->weights()[index].wgt/lheeventinfo->originalXWGTUP();
			nLHEPDFWeights++;
		}
		if(index>=(nlhescalemax+nPDFsets) && index<(nlhescalemax+nPDFsets+nalpsmax) && nLHEAlpsWeights<nalpsmax){
			LHEAlpsWeights[nLHEAlpsWeights] = lheeventinfo->weights()[index].wgt/lheeventinfo->originalXWGTUP();
			nLHEAlpsWeights++;
			}
	  }
	  		
	}

    // Flavor tagging of GEN jets using ghost-matching //                                                             

	edm::Handle<reco::JetFlavourInfoMatchingCollection> jetFlavourInfos;
    iEvent.getByToken(jetFlavourInfosToken_, jetFlavourInfos);

    ngenjetAK8 = 0;
    iEvent.getByToken(tok_genjetAK8s_, genjetAK8s);
    
    JetDefinition pfjetAK8Def(antikt_algorithm,0.8,E_scheme);
	SoftDrop sd(beta,z_cut,0.8);
    
    if(genjetAK8s.isValid()){
		
		std::vector<int> partonFlavour_AK8;
		std::vector<uint8_t> hadronFlavour_AK8;

		for (const reco::GenJet & jet : *genjetAK8s) {
			
			bool matched = false;
			for (const reco::JetFlavourInfoMatching & jetFlavourInfoMatching : *jetFlavourInfos) {
				if (deltaR(jet.p4(), jetFlavourInfoMatching.first->p4()) < 0.1) {
					partonFlavour_AK8.push_back(jetFlavourInfoMatching.second.getPartonFlavour());
					hadronFlavour_AK8.push_back(jetFlavourInfoMatching.second.getHadronFlavour());
					matched = true;
					break;
				}
			}
		
			if (!matched) {
				partonFlavour_AK8.push_back(-100);
				hadronFlavour_AK8.push_back(-100);
			}
		}
      
		for(unsigned gjet = 0; gjet<genjetAK8s->size(); gjet++)	{
	
			TLorentzVector genjetAK84v((*genjetAK8s)[gjet].px(),(*genjetAK8s)[gjet].py(),(*genjetAK8s)[gjet].pz(), (*genjetAK8s)[gjet].energy());
			if(genjetAK84v.Pt()<AK8GenPtCut) continue;
			if(abs(genjetAK84v.Eta())>maxgenEta) continue;
	
			genjetAK8pt[ngenjetAK8] = genjetAK84v.Pt();
			genjetAK8eta[ngenjetAK8] = genjetAK84v.Eta();
			genjetAK8phi[ngenjetAK8] = genjetAK84v.Phi();
			genjetAK8mass[ngenjetAK8] = (*genjetAK8s)[gjet].mass();
			genjetAK8hadronflav[ngenjetAK8] = (int)hadronFlavour_AK8[gjet];
			genjetAK8partonflav[ngenjetAK8] = partonFlavour_AK8[gjet];
	
			std::vector<reco::CandidatePtr> daught((*genjetAK8s)[gjet].daughterPtrVector());
	
			vector <fastjet::PseudoJet> fjInputs;
			fjInputs.resize(0);
			for (unsigned int i2 = 0; i2< daught.size(); ++i2) {
				PseudoJet psjet = PseudoJet( (*daught[i2]).px(),(*daught[i2]).py(),(*daught[i2]).pz(),(*daught[i2]).energy() );
				psjet.set_user_index(i2);
				fjInputs.push_back(psjet); 
			} //i2
	
			vector <fastjet::PseudoJet> sortedJets;
			fastjet::ClusterSequence clustSeq(fjInputs, pfjetAK8Def);
			fjInputs.clear();
			sortedJets    = sorted_by_pt(clustSeq.inclusive_jets());
	
			if(sortedJets.size()>0){
				genjetAK8sdmass[ngenjetAK8] = (sd(sortedJets[0])).m();
			}
	
			if (++ngenjetAK8>=njetmx) break;
	
		}
	}
      
      
	ngenjetAK4 = 0;
	
	if(genjetAK4s.isValid()){
		
		std::vector<int> partonFlavour_AK4;
		std::vector<uint8_t> hadronFlavour_AK4;
      
		for (const reco::GenJet & jet : *genjetAK4s) {
		  
			bool matched = false;
			for (const reco::JetFlavourInfoMatching & jetFlavourInfoMatching : *jetFlavourInfos) {
				if (deltaR(jet.p4(), jetFlavourInfoMatching.first->p4()) < 0.1) {
					partonFlavour_AK4.push_back(jetFlavourInfoMatching.second.getPartonFlavour());
					hadronFlavour_AK4.push_back(jetFlavourInfoMatching.second.getHadronFlavour());
					matched = true;
					break;
				}
			}
		
			if (!matched) {
				partonFlavour_AK4.push_back(-100);
				hadronFlavour_AK4.push_back(-100);
			}	
		}
	
		for(unsigned gjet = 0; gjet<genjetAK4s->size(); gjet++)	{
	
			TLorentzVector genjetAK44v((*genjetAK4s)[gjet].px(),(*genjetAK4s)[gjet].py(),(*genjetAK4s)[gjet].pz(), (*genjetAK4s)[gjet].energy());
			if(genjetAK44v.Pt()<minGenPt) continue;
			if(abs(genjetAK44v.Eta())>maxgenEta) continue;
	
			genjetAK4pt[ngenjetAK4] = genjetAK44v.Pt();
			genjetAK4eta[ngenjetAK4] = genjetAK44v.Eta();
			genjetAK4phi[ngenjetAK4] = genjetAK44v.Phi();
			genjetAK4mass[ngenjetAK4] = (*genjetAK4s)[gjet].mass();
			genjetAK4hadronflav[ngenjetAK4] = (int)hadronFlavour_AK4[gjet];
			genjetAK4partonflav[ngenjetAK4] = partonFlavour_AK4[gjet];

			if (++ngenjetAK4>=njetmx) break;
      
		}
		
    }
    
	ngenparticles = 0;
    edm::Handle<std::vector<reco::GenParticle>> genparticles;  
	iEvent.getByToken(tok_genparticles_,genparticles);
    
	if(genparticles.isValid()){
	
		for(unsigned ig=0; ig<(genparticles->size()); ig++){
			
			if(!(((*genparticles)[ig].status()==1)||((*genparticles)[ig].status()==22)||((*genparticles)[ig].status()==23))) continue;
			if(!((*genparticles)[ig].isHardProcess())) continue;
	  
			if(!((abs((*genparticles)[ig].pdgId())>=1&&abs((*genparticles)[ig].pdgId())<=6)||(abs((*genparticles)[ig].pdgId())>=11&&abs((*genparticles)[ig].pdgId())<=16)||(abs((*genparticles)[ig].pdgId())==24))) continue;
	  
			const Candidate * mom = (*genparticles)[ig].mother();
	  
			genpartstatus[ngenparticles] = (*genparticles)[ig].status();
			genpartpdg[ngenparticles] = (*genparticles)[ig].pdgId();
			genpartmompdg[ngenparticles] = mom->pdgId();
			const Candidate * momtmp = (*genparticles)[ig].mother();
			while(genpartmompdg[ngenparticles] == genpartpdg[ngenparticles])
			{
				genpartmompdg[ngenparticles] = momtmp->mother()->pdgId();
				momtmp = momtmp->mother();
			}
			genpartdaugno[ngenparticles] = (*genparticles)[ig].numberOfDaughters();
			genpartfromhard[ngenparticles] = (*genparticles)[ig].isHardProcess();
			genpartfromhardbFSR[ngenparticles] = (*genparticles)[ig].fromHardProcessBeforeFSR();
			genpartisLastCopyBeforeFSR[ngenparticles] = (*genparticles)[ig].isLastCopyBeforeFSR();
			genpartisPromptFinalState[ngenparticles] = (*genparticles)[ig].isPromptFinalState();
			genpartpt[ngenparticles] = (*genparticles)[ig].pt();
			genparteta[ngenparticles] = (*genparticles)[ig].eta();
			genpartphi[ngenparticles] = (*genparticles)[ig].phi();
			genpartm[ngenparticles] = (*genparticles)[ig].mass();
	  
			bool found_grmom = false;
			if(mom->numberOfMothers()>0){
				const Candidate * grmom  = mom->mother();
				for(int iter=0; iter<10; iter++){
					if(grmom->pdgId() != mom->pdgId()){
						genpartgrmompdg[ngenparticles]  = grmom->pdgId();
						found_grmom = true;
						break;
					}else{
						if(grmom->numberOfMothers()>0){
							grmom = grmom->mother();
						}
						else{ break; }
					}
				}
			}
			if(!found_grmom){
				genpartgrmompdg[ngenparticles]  = -10000000;
			}
	  
			ngenparticles++;
			if(ngenparticles>=npartmx) break;
		}
    }
    
    // pileup information //
    
    npu_vert = 0;
	npu_vert_true = 0;
    
    edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByToken(pileup_, PupInfo);
    if (PupInfo.isValid()) {
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
		if (PVI->getBunchCrossing()==0) {
			npu_vert = PVI->getPU_NumInteractions();
			npu_vert_true = PVI->getTrueNumInteractions();
			break;
			}
		}
    }

  }//isMC
  
  Handle<VertexCollection> primaryVertices;
  iEvent.getByToken(tok_primaryVertices_, primaryVertices);
  reco::Vertex vertex;
  
  if (primaryVertices.isValid()) {
	  
    vertex = primaryVertices->at(0);  
    int ndofct_org=0;
    int nchict_org=0;
    int nvert_org = 0;
    int nprimi_org = 0;
    
    for (reco::VertexCollection::const_iterator vert=primaryVertices->begin(); vert<primaryVertices->end(); vert++) {
      nvert_org++;
      if (vert->isValid() && !vert->isFake()) {
		if (vert->ndof() > 4 && fabs(vert->position().z()) <= 24 && fabs(vert->position().Rho()) <= 2) {
			nprimi_org++;
			}
		if (vert->ndof()>7) {
			ndofct_org++;
			if (vert->normalizedChi2()<5) nchict_org++;
			}
		}
    }
    
    nprim = min(999,nvert_org) + 1000*min(999,ndofct_org) + 1000000*min(999,nchict_org);
    npvert = nchict_org;
    
  } else { 
    nprim = -100;
    npvert = -100;
  }
  
  reco::TrackBase::Point beamPoint(0,0, 0);
  edm::Handle<reco::BeamSpot> beamSpotH;
  
  iEvent.getByToken(tok_beamspot_, beamSpotH);  //Label("offlineBeamSpot",beamSpotH);
  if (beamSpotH.isValid()){
    beamPoint = beamSpotH->position();
  }

  edm::Handle<reco::VertexCompositePtrCandidateCollection> svin;
  iEvent.getByToken(tok_sv,svin);
  
  edm::Handle<double> Rho_PF;
  
  iEvent.getByToken(tok_Rho_,Rho_PF);
  Rho = *Rho_PF;
  
  // Store trigger information //
  
  const char* variab1;
  
  edm::Handle<edm::TriggerResults> trigRes;
  iEvent.getByToken(triggerBits_, trigRes);
  
  const edm::TriggerNames &names = iEvent.triggerNames(*trigRes);
  
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  
  int ihlttrg[nHLTmx+1]= {0};
  bool booltrg[nHLTmx]= {false};
  
  for (int jk=0; jk<nHLTmx; jk++) {
    for(unsigned ij = 0; ij<trigRes->size(); ++ij) {
      std::string name = names.triggerName(ij);
      variab1 = name.c_str(); 
      if (strstr(variab1,hlt_name[jk]) && ((strlen(variab1)-strlen(hlt_name[jk]))<5))
		{
			if ((trigRes->accept(ij))){   //||(isMC)) {
				ihlttrg[jk] = ihlttrg[nHLTmx] = 1;
				booltrg[jk] = true;
				break;
			}
		}
    }//ij     
  }//jk
  
  trig_value = 1; 
  
  for (int jk=1; jk<(nHLTmx+1); jk++) {
    if(ihlttrg[nHLTmx-jk]>0) {
      trig_value+=(1<<jk);
    }
  }
  
  // Trigger objects //
  
  vector<triggervar> alltrgobj;
  if (trigRes.isValid()) { 
    
    const char* variab2 ;
    
    alltrgobj.clear(); 
    
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
      
      obj.unpackPathNames(names);
      std::vector<std::string> pathNamesAll  = obj.pathNames(false);
      
      for (unsigned ih = 0, n = pathNamesAll.size(); ih < n; ih++) {
	
		variab2 = pathNamesAll[ih].c_str(); 
	
		for (int jk=0; jk<nHLTmx; jk++) {
			if (strstr(variab2,hlt_name[jk]) && (strlen(variab2)-strlen(hlt_name[jk])<5)) {
	    
				if(obj.pt()>20 && fabs(obj.eta())<3.0) {
	      
					triggervar tmpvec1;
	      
					tmpvec1.both = obj.hasPathName( pathNamesAll[ih], true, true );
					tmpvec1.highl  = obj.hasPathName( pathNamesAll[ih], false, true );
					tmpvec1.level1 = obj.hasPathName( pathNamesAll[ih], true, false );
					tmpvec1.trg4v = TLorentzVector(obj.px(), obj.py(), obj.pz(), obj.energy());
					tmpvec1.pdgId = obj.pdgId();
					tmpvec1.prescl = 1;    //triggerPrescales->getPrescaleForIndex(ih);
					tmpvec1.ihlt = jk;
					alltrgobj.push_back(tmpvec1);
					break;
				}
			}
		}//jk 
      }//ih
    }
  }
  
  int xht=0;
  ntrigobjs = alltrgobj.size();
  if(ntrigobjs>njetmx) { ntrigobjs = njetmx; }
  if(ntrigobjs > 0){
    for(unsigned int iht=0; iht<(unsigned int)ntrigobjs; iht++){
      if(alltrgobj[iht].trg4v.Pt()>20 && fabs(alltrgobj[iht].trg4v.Eta())<3.0) {
		trigobjpt[xht] = alltrgobj[iht].trg4v.Pt();
		trigobjeta[xht] = alltrgobj[iht].trg4v.Eta();
		trigobjphi[xht] = alltrgobj[iht].trg4v.Phi();
		trigobjmass[xht] = alltrgobj[iht].trg4v.M();
		trigobjHLT[xht] = alltrgobj[iht].highl;
		trigobjL1[xht] = alltrgobj[iht].level1;
		trigobjBoth[xht] = alltrgobj[iht].both;
		trigobjIhlt[xht] = alltrgobj[iht].ihlt;
		trigobjpdgId[xht] = alltrgobj[iht].pdgId;
		xht++;
      }
      if(iht == (njetmx-1)) break;
    }
  }
  
  // End of trigger info //
  
  // Prefire weights //
  
  if(add_prefireweights){
  
	edm::Handle< double > theprefweight;
	iEvent.getByToken(prefweight_token, theprefweight ) ;	
	prefiringweight =(*theprefweight);

	edm::Handle< double > theprefweightup;
	iEvent.getByToken(prefweightup_token, theprefweightup ) ;
	prefiringweightup =(*theprefweightup);
    
	edm::Handle< double > theprefweightdown;
	iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;   
	prefiringweightdown =(*theprefweightdown);
  
  }
  
  // End of prefire weights //
  
  // ====== RECO-objects now  ==========//
  
  // MET //
  
  miset = misphi = misetsig = sumEt = genmiset = genmisphi = genmisetsig = -1000 ;
  
  edm::Handle<pat::METCollection> pfmet_ ;
  iEvent.getByToken(tok_mets_,pfmet_) ;
  
  edm::Handle<pat::METCollection> pfmet_PUPPI_ ;
  iEvent.getByToken(tok_mets_PUPPI_,pfmet_PUPPI_) ;
  
  if(pfmet_.isValid()){
	  
    const pat::MET &met = pfmet_->front();
    
    miset = met.corPt(); //met.pt();
    misphi = met.corPhi();//met.phi();
    misetsig = met.significance();
    sumEt = met.corSumEt();//sumEt();
    
    if(isMC){
      genmiset = met.genMET()->pt();
      genmisphi = met.genMET()->phi();
      genmisetsig = met.genMET()->significance();
    }
  }
  
  if(pfmet_PUPPI_.isValid()){
	
	const pat::MET &met = pfmet_PUPPI_->front();
	  
	miset_PUPPI = met.corPt(); 
    misphi_PUPPI = met.corPhi();
    misetsig_PUPPI = met.significance();
    sumEt_PUPPI = met.corSumEt();
	  
  }
  
  // Muons //
  
  nmuons = 0;                                                                                                                                        
  std::vector<pat::Muon> tlvmu;
  edm::Handle<edm::View<pat::Muon>> muons;                                                                                                          
  iEvent.getByToken(tok_muons_, muons);                                                                                                             
    
  if(muons.isValid() && muons->size()>0) {                                                                                                           
    
	edm::View<pat::Muon>::const_iterator muon1;                                                                                                      

    for( muon1 = muons->begin(); muon1 < muons->end(); muon1++ ) {                                                                                   

		if (StoreMuon(*muon1,minmuPt,maxEta)) {                                                                
			
			muonpt[nmuons] = muon1->pt();                                                                         
			TrackRef trktrk = muon1->innerTrack();                                                                                                       
			muonp[nmuons] = trktrk->p()*muon1->charge();                                                                                                                 
			muoneta[nmuons] = muon1->eta();                                                                                                              
			muonphi[nmuons] = muon1->phi();                                                                                                                                                                                                                                
                                                                                                                                                       
			//MiniIsolation: begin//                                                                                      
			vector<float> isovalues;
			Read_MiniIsolation(muon1,Rho,isovalues);
			muonminisoall[nmuons] = isovalues[0];
			muonminchiso[nmuons] = isovalues[1];
			muonminnhiso[nmuons] = isovalues[2];
			muonminphiso[nmuons] = isovalues[3];
			//MiniIsolation: end//  
			                                                              
			muonisPF[nmuons] = muon1->isPFMuon();                                                                                                        
			muonisGL[nmuons] = muon1->isGlobalMuon();                                                                                                    
			muonisTRK[nmuons] = muon1->isTrackerMuon();                                                                                                  
			muonisLoose[nmuons] = (muon::isLooseMuon(*muon1));                                                                                           
			muonisMed[nmuons] = (muon::isMediumMuon(*muon1));                                                                                            
			muonisMedPr[nmuons] = false;                                                                          
			if(muon::isMediumMuon(*muon1)) {                                                                                                             
				if ((std::abs(muon1->muonBestTrack()->dz(vertex.position())) < 0.1) && (std::abs(muon1->muonBestTrack()->dxy(vertex.position())) < 0.02)){                                                                                                                  
					muonisMedPr[nmuons] = true;                                                                                                              
				}                                                                                                                                          
			}                                                                                                                                      
			muonisGoodGL[nmuons] = (muon1->isGlobalMuon() && muon1->globalTrack()->normalizedChi2() < 3 && muon1->combinedQuality().chi2LocalPosition < 12 && muon1->combinedQuality().trkKink < 20 && (muon::segmentCompatibility(*muon1)) > 0.303);                     
			muonisTight[nmuons] = (muon::isTightMuon(*muon1,vertex));                                                                                    
			muonisHighPt[nmuons] = (muon::isHighPtMuon(*muon1,vertex));                                                                                  
			muonisHighPttrk[nmuons] = (muon::isTrackerHighPtMuon(*muon1,vertex));   
			                                                                     
			muonecal[nmuons] = (muon1->calEnergy()).em;                                                                                                  
			muonhcal[nmuons] = (muon1->calEnergy()).had;     
			                                                                            
			muonposmatch[nmuons] = muon1->combinedQuality().chi2LocalPosition;                                                                           
			muontrkink[nmuons] = muon1->combinedQuality().trkKink;                                                                                       
			muonsegcom[nmuons] = muon::segmentCompatibility(*muon1);                                                                                     
			muontrkvtx[nmuons] = muon1->muonBestTrack()->dxy(vertex.position());                                                                         
			muondz[nmuons] = muon1->muonBestTrack()->dz(vertex.position());          
			                                                                    
			float dzmumin = 1000;                                                                                                                        
			float dxymumin = 1000;                                                                                                                       
			if(svin.isValid()){                                                                                                                          
				for(unsigned int isv=0; isv<(svin->size()); isv++){                                                                                        
					const auto &sv = (*svin)[isv];                                                                                                           
					reco::TrackBase::Point svpoint(sv.vx(),sv.vy(),sv.vz());
					float dztmp = fabs(muon1->muonBestTrack()->dz(svpoint));
					if(dztmp < dzmumin){
						dzmumin = dztmp;                                                                                   
						dxymumin = muon1->muonBestTrack()->dxy(svpoint);                                                                                       
					}                                                                                                                                        
				}                                                                                                                                          
			}                                                                                                                                            
			mudxy_sv[nmuons] = dxymumin;                                                                                                                 
			muonpter[nmuons] = trktrk->ptError();                                                                                                        
			TrackRef trkglb =muon1->globalTrack();                                                                                                       
			if ((!muon1->isGlobalMuon())) {                                                                                                              
				if (muon1->isTrackerMuon()) {                                                                                                              
					trkglb =muon1->innerTrack();                                                                                                             
				} else {                                                                                                                                   
					trkglb =muon1->outerTrack();                                                                                                             
				}                                                                                                                                          
			}
			                                                                                                                                            
			muonchi[nmuons] = trkglb->normalizedChi2();                                                                                                  
			muonndf[nmuons] = (int)trkglb->ndof();                                                                                                       
			muonhit[nmuons] = trkglb->hitPattern().numberOfValidMuonHits();                                                                              
			muonmst[nmuons] = muon1->numberOfMatchedStations();                                                                                          
			muonpixhit[nmuons] = trktrk->hitPattern().numberOfValidPixelHits();                                                                          
			muontrklay[nmuons] = trktrk->hitPattern().trackerLayersWithMeasurement();                                                                    
			muonvalfrac[nmuons] = trktrk->validFraction();                                                        
			muonpfiso[nmuons] = (muon1->pfIsolationR04().sumChargedHadronPt + max(0., muon1->pfIsolationR04().sumNeutralHadronEt + muon1->pfIsolationR04().sumPhotonEt - 0.5*muon1->pfIsolationR04().sumPUPt))/muon1->pt();                                               
			
			TLorentzVector tlmu;
			bool mu_id = Muon_TightID(muonisGL[nmuons],muonisPF[nmuons],
				    muonchi[nmuons],muonhit[nmuons],muonmst[nmuons],
				    muontrkvtx[nmuons],muondz[nmuons],muonpixhit[nmuons],muontrklay[nmuons]);
			if (muonpt[nmuons]>15 && fabs(muoneta[nmuons])<2.5 && mu_id && muontrkvtx[nmuons]<0.2 && muondz[nmuons]<0.5) {
				tlvmu.push_back(*muon1);
			}
	  
			if (++nmuons>=njetmx) break;                                                                                                                 
		
		}                                                                                                                                              
      }                                                                                                                                               
  }// muon loop 
  
  // Electrons //
    
  nelecs = 0;             
  std::vector<pat::Electron> tlvel;
  int iE1 = 0;                                                                                                                                       
    
  for(const auto& electron1 : iEvent.get(tok_electrons_) ) {                                                                                          
 
    if (!StoreElectron(electron1,minePt,maxEta)) continue;
    iE1++;       
    
	elmvaid_Fallv2WP90[nelecs] = electron1.electronID(melectronID_isowp90);                                                                                 
    elmvaid_Fallv2WP90_noIso[nelecs] = electron1.electronID(melectronID_noisowp90);                                                                             
    elmvaid_Fallv2WP80[nelecs] = electron1.electronID(melectronID_isowp80);                                                                                 
    elmvaid_Fallv2WP80_noIso[nelecs] = electron1.electronID(melectronID_noisowp80);                                                                             
                                                   
	GsfTrackRef gsftrk1 = electron1.gsfTrack();   																														
    TrackRef ctftrk = electron1.closestCtfTrackRef();    
    
    elpt[nelecs] = electron1.pt();                                                                                                                   
    eleta[nelecs] = electron1.eta();                                                                                                                 
    elphi[nelecs] = electron1.phi();                                                                                                                 
    ele[nelecs] = electron1.ecalEnergy();                                                                                                            
    elp[nelecs] = electron1.trackMomentumAtVtx().R()*electron1.charge();                                                                                                
    elsupcl_eta[nelecs] = electron1.superCluster()->eta();                                                                                           
    elsupcl_phi[nelecs] = electron1.superCluster()->phi();                                                                                           
    elsupcl_rawE[nelecs] = electron1.superCluster()->rawEnergy();                                                                                    
    eldxytrk[nelecs] = gsftrk1->dxy(vertex.position());                                                                                              
    eldztrk[nelecs] = gsftrk1->dz(vertex.position());                                                                                                   
                                                                                                                                                                                                                                              
	elsigmaieta[nelecs] = electron1.full5x5_sigmaIetaIeta();                                                                                         
    elsigmaiphi[nelecs] = electron1.full5x5_sigmaIphiIphi();                                                                                         
    elr9full[nelecs] = electron1.full5x5_r9();                                                                                                       
    elsupcl_etaw[nelecs] = electron1.superCluster()->etaWidth();                                                                                     
    elsupcl_phiw[nelecs] = electron1.superCluster()->phiWidth();                                                                                     
    elhcaloverecal[nelecs] = electron1.full5x5_hcalOverEcal();                                                                                       
    elcloctftrkn[nelecs] = electron1.closestCtfTrackNLayers();                                                                                       
    elcloctftrkchi2[nelecs] = electron1.closestCtfTrackNormChi2();                                                                                   
    ele1x5bye5x5[nelecs] = 1.-electron1.full5x5_e1x5()/electron1.full5x5_e5x5();                                                                     
    elnormchi2[nelecs] =  electron1.gsfTrack()->normalizedChi2();                                                                                    
    elhitsmiss[nelecs] =  electron1.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);                                 
    eltrkmeasure[nelecs] = electron1.gsfTrack()->hitPattern().trackerLayersWithMeasurement();                                                        
    elconvtxprob[nelecs] = electron1.convVtxFitProb();                                                                                               
    elecloverpout[nelecs] = electron1.eEleClusterOverPout();                                                                                         
    elecaletrkmomentum[nelecs] = 1.0/(electron1.ecalEnergy())-1.0/(electron1.trackMomentumAtVtx().R());                                              
    eldeltaetacltrkcalo[nelecs] = electron1.deltaEtaSeedClusterTrackAtCalo();                                                                        
    elsupcl_preshvsrawe[nelecs] = electron1.superCluster()->preshowerEnergy()/electron1.superCluster()->rawEnergy();                                 
    elpfisolsumphet[nelecs] = electron1.pfIsolationVariables().sumPhotonEt;                                                                          
    elpfisolsumchhadpt[nelecs] = electron1.pfIsolationVariables().sumChargedHadronPt;                                                                
    elpfsiolsumneuhadet[nelecs] = electron1.pfIsolationVariables().sumNeutralHadronEt;                                                               
    eleoverp[nelecs] = electron1.eSuperClusterOverP();                                                                                               
    elhovere[nelecs] = electron1.hadronicOverEm();                                                                                                   
    elietaieta[nelecs] = electron1.sigmaIetaIeta();                                                                                                  
    eletain[nelecs] = electron1.deltaEtaSuperClusterTrackAtVtx();                                                                                    
    elphiin[nelecs] = electron1.deltaPhiSuperClusterTrackAtVtx();                                                                                    
    elfbrem[nelecs] = electron1.fbrem();   
                                                                                                                
    const reco::GsfElectron::PflowIsolationVariables& pfIso = electron1.pfIsolationVariables();                                                      
    elpfiso_drcor[nelecs] = (pfIso.sumChargedHadronPt + max(0., pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5*pfIso.sumPUPt))*1./electron1.pt();      
    vector<float> pfisovalues;                                                                                     
    Read_ElePFIsolation(&electron1,Rho,pfisovalues);
    elpfiso_eacor[nelecs] = pfisovalues[0];
    elpfiso04_eacor[nelecs] = pfisovalues[1];
    
    float dzmin = 1000;                                                                                                                              
    float dxymin = 1000;
    if(svin.isValid()){                                                                                                                              
		for(unsigned int isv=0; isv<(svin->size()); isv++){                                                                                            
		const auto &sv = (*svin)[isv];                                                                                                               
     	  reco::TrackBase::Point svpoint(sv.vx(),sv.vy(),sv.vz());
		  float dztmp =fabs(gsftrk1->dz(svpoint));
		  if(dztmp < dzmin){                                                                                                      
			dzmin = dztmp;                                                                                                        
			dxymin = gsftrk1->dxy(svpoint);                                                                                                            
			}                                                                                                                                            
		}                                                                                                                                              
    }     
                                                                                                                                                 
    eldxy_sv[nelecs] = dxymin;                                                                                                                       
    elchi[nelecs] = gsftrk1->chi2();                                                                                                                 
    elndf[nelecs] = (int)gsftrk1->ndof();                                                                                                            
    elmisshits[nelecs] = (int)gsftrk1->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    
    //MiniIsolation: begin//                                                                                      
	vector<float> isovalues;
	Read_MiniIsolation(&electron1,Rho,isovalues);
	elminisoall[nmuons] = isovalues[0];
	elminchiso[nmuons] = isovalues[1];
	elminnhiso[nmuons] = isovalues[2];
	elminphiso[nmuons] = isovalues[3];
	//MiniIsolation: end//  
	
	bool impact_pass = 	((fabs(elsupcl_eta[nelecs])<1.4442 && fabs(eldxytrk[nelecs])<0.05 && fabs(eldztrk[nelecs])<0.1)
					   ||(fabs(elsupcl_eta[nelecs])>1.5660 && fabs(eldxytrk[nelecs])<(2*0.05) && fabs(eldztrk[nelecs])<(2*0.1)));

	
	if(elpt[nelecs]>15 && fabs(eleta[nelecs])<2.5 && elmvaid_Fallv2WP90_noIso[nelecs] && impact_pass){
		tlvel.push_back(electron1);
	}
  
    if(++nelecs>=njetmx) break;                                                                                                                      
  }
    
  // PF candidates //
  
  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(tok_pfcands_, pfs);  
  
  // AK8 jets //
  
  npfjetAK8 = 0;
  edm::Handle<edm::View<pat::Jet>> pfjetAK8s;
  iEvent.getByToken(tok_pfjetAK8s_, pfjetAK8s);	
  
  if(pfjetAK8s.isValid()){
    
    for (unsigned jet = 0; jet< pfjetAK8s->size(); jet++) {
      
      const auto &ak8jet = (*pfjetAK8s)[jet];

      TLorentzVector pfjetAK84v(ak8jet.correctedP4("Uncorrected").px(),ak8jet.correctedP4("Uncorrected").py(),ak8jet.correctedP4("Uncorrected").pz(), ak8jet.correctedP4("Uncorrected").energy());
     
      pfjetAK84v = LeptonJet_subtraction(tlvmu,ak8jet,pfjetAK84v);
	  pfjetAK84v = LeptonJet_subtraction(tlvel,ak8jet,pfjetAK84v);
	
      double tmprecpt = pfjetAK84v.Pt();
      
      double total_cor =1;
      Read_JEC(total_cor,tmprecpt,pfjetAK84v.Eta(),Rho,isData,ak8jet,jecL1FastAK8,jecL2RelativeAK8,jecL3AbsoluteAK8,jecL2L3ResidualAK8);
      pfjetAK8JEC[npfjetAK8] = total_cor;
            
      if(tmprecpt <AK8PtCut) continue;
      if(abs(pfjetAK84v.Rapidity())>maxEta) continue;
      
      pfjetAK8pt[npfjetAK8] = 	pfjetAK84v.Pt();
      pfjetAK8y[npfjetAK8] = pfjetAK84v.Rapidity();
      pfjetAK8eta[npfjetAK8] = pfjetAK84v.Eta();
      pfjetAK8phi[npfjetAK8] = pfjetAK84v.Phi();
      pfjetAK8mass[npfjetAK8] = ak8jet.correctedP4("Uncorrected").mass();
      pfjetAK8btag_DeepCSV[npfjetAK8] = ak8jet.bDiscriminator("pfDeepCSVJetTags:probb")+ak8jet.bDiscriminator("pfDeepCSVJetTags:probbb");
      
      pfjetAK8DeepTag_DAK8_TvsQCD[npfjetAK8] = ak8jet.bDiscriminator(toptagger_DAK8);
      pfjetAK8DeepTag_DAK8_WvsQCD[npfjetAK8] = ak8jet.bDiscriminator(Wtagger_DAK8);
      pfjetAK8DeepTag_DAK8_ZvsQCD[npfjetAK8] = ak8jet.bDiscriminator(Ztagger_DAK8);
      pfjetAK8DeepTag_DAK8_HvsQCD[npfjetAK8] = ak8jet.bDiscriminator(Htagger_DAK8);
      pfjetAK8DeepTag_DAK8_bbvsQCD[npfjetAK8] = ak8jet.bDiscriminator(bbtagger_DAK8);
      
      pfjetAK8DeepTag_PNet_TvsQCD[npfjetAK8] = ak8jet.bDiscriminator(toptagger_PNet);
      pfjetAK8DeepTag_PNet_WvsQCD[npfjetAK8] = ak8jet.bDiscriminator(Wtagger_PNet);
      pfjetAK8DeepTag_PNet_ZvsQCD[npfjetAK8] = ak8jet.bDiscriminator(Ztagger_PNet);
      pfjetAK8DeepTag_PNet_XbbvsQCD[npfjetAK8] = ak8jet.bDiscriminator(Xbbtagger_PNet);
      pfjetAK8DeepTag_PNet_XccvsQCD[npfjetAK8] = ak8jet.bDiscriminator(Xcctagger_PNet);
      pfjetAK8DeepTag_PNet_XqqvsQCD[npfjetAK8] = ak8jet.bDiscriminator(Xqqtagger_PNet);
       
      if(isMC){
	
      	vector<double> SFs;
      	Read_JER(mPtResoFileAK8, mPtSFFileAK8, tmprecpt, pfjetAK84v, Rho, genjetAK8s, SFs);
      	
      	pfjetAK8reso[npfjetAK8] = SFs[0];
      	pfjetAK8resoup[npfjetAK8] = SFs[1];
      	pfjetAK8resodn[npfjetAK8] = SFs[2];
      	      	
      }//isMC
      
      for(int isrc =0 ; isrc<njecmcmx; isrc++){
	
        double sup = 1.0 ;
	
		if((isrc>0)&&(isrc<=nsrc)){
		  
		  JetCorrectionUncertainty *jecUnc = vsrcAK8[isrc-1];
		  jecUnc->setJetEta(ak8jet.eta());
		  jecUnc->setJetPt(tmprecpt);
		  
		  sup += jecUnc->getUncertainty(true);         
		  if(isrc==1){ pfjetAK8jesup_AbsoluteStat[npfjetAK8] = sup; }
		  if(isrc==2){ pfjetAK8jesup_AbsoluteScale[npfjetAK8] = sup; }
		  if(isrc==3){ pfjetAK8jesup_AbsoluteMPFBias[npfjetAK8] = sup; }
		  if(isrc==4){ pfjetAK8jesup_FlavorQCD[npfjetAK8] = sup; }
		  if(isrc==5){ pfjetAK8jesup_Fragmentation[npfjetAK8] = sup; }
		  if(isrc==6){ pfjetAK8jesup_PileUpDataMC[npfjetAK8] = sup; }
		  if(isrc==7){ pfjetAK8jesup_PileUpPtBB[npfjetAK8] = sup; }
		  if(isrc==8){ pfjetAK8jesup_PileUpPtEC1[npfjetAK8] = sup; }
		  if(isrc==9){ pfjetAK8jesup_PileUpPtEC2[npfjetAK8] = sup; }
		  if(isrc==10){ pfjetAK8jesup_PileUpPtRef[npfjetAK8] = sup; }
		  if(isrc==11){ pfjetAK8jesup_RelativeFSR[npfjetAK8] = sup; }
		  if(isrc==12){ pfjetAK8jesup_RelativeJEREC1[npfjetAK8] = sup; }
		  if(isrc==13){ pfjetAK8jesup_RelativeJEREC2[npfjetAK8] = sup; }
		  if(isrc==14){ pfjetAK8jesup_RelativePtBB[npfjetAK8] = sup; }
		  if(isrc==15){ pfjetAK8jesup_RelativePtEC1[npfjetAK8] = sup; }
		  if(isrc==16){ pfjetAK8jesup_RelativePtEC2[npfjetAK8] = sup; }
		  if(isrc==17){ pfjetAK8jesup_RelativeBal[npfjetAK8] = sup; }
		  if(isrc==18){ pfjetAK8jesup_RelativeSample[npfjetAK8] = sup; }
		  if(isrc==19){ pfjetAK8jesup_RelativeStatEC[npfjetAK8] = sup; }
		  if(isrc==20){ pfjetAK8jesup_RelativeStatFSR[npfjetAK8] = sup; }
		  if(isrc==21){ pfjetAK8jesup_SinglePionECAL[npfjetAK8] = sup; }
		  if(isrc==22){ pfjetAK8jesup_SinglePionHCAL[npfjetAK8] = sup; }
		  if(isrc==23){ pfjetAK8jesup_TimePtEta[npfjetAK8] = sup; }
		  if(isrc==24){ pfjetAK8jesup_Total[npfjetAK8] = sup; }
		
		}
		else if(isrc>nsrc){
		  
		  JetCorrectionUncertainty *jecUnc = vsrcAK8[isrc-1-nsrc];
		  jecUnc->setJetEta(ak8jet.eta());
		  jecUnc->setJetPt(tmprecpt);
		  
		  sup -= jecUnc->getUncertainty(false);
		  if(isrc==(nsrc+1)){ pfjetAK8jesdn_AbsoluteStat[npfjetAK8] = sup; }
		  if(isrc==(nsrc+2)){ pfjetAK8jesdn_AbsoluteScale[npfjetAK8] = sup; }
		  if(isrc==(nsrc+3)){ pfjetAK8jesdn_AbsoluteMPFBias[npfjetAK8] = sup; }
		  if(isrc==(nsrc+4)){ pfjetAK8jesdn_FlavorQCD[npfjetAK8] = sup; }
		  if(isrc==(nsrc+5)){ pfjetAK8jesdn_Fragmentation[npfjetAK8] = sup; }
		  if(isrc==(nsrc+6)){ pfjetAK8jesdn_PileUpDataMC[npfjetAK8] = sup; }
		  if(isrc==(nsrc+7)){ pfjetAK8jesdn_PileUpPtBB[npfjetAK8] = sup; }
		  if(isrc==(nsrc+8)){ pfjetAK8jesdn_PileUpPtEC1[npfjetAK8] = sup; }
		  if(isrc==(nsrc+9)){ pfjetAK8jesdn_PileUpPtEC2[npfjetAK8] = sup; }
		  if(isrc==(nsrc+10)){ pfjetAK8jesdn_PileUpPtRef[npfjetAK8] = sup; }
		  if(isrc==(nsrc+11)){ pfjetAK8jesdn_RelativeFSR[npfjetAK8] = sup; }
		  if(isrc==(nsrc+12)){ pfjetAK8jesdn_RelativeJEREC1[npfjetAK8] = sup; }
		  if(isrc==(nsrc+13)){ pfjetAK8jesdn_RelativeJEREC2[npfjetAK8] = sup; }
		  if(isrc==(nsrc+14)){ pfjetAK8jesdn_RelativePtBB[npfjetAK8] = sup; }
		  if(isrc==(nsrc+15)){ pfjetAK8jesdn_RelativePtEC1[npfjetAK8] = sup; }
		  if(isrc==(nsrc+16)){ pfjetAK8jesdn_RelativePtEC2[npfjetAK8] = sup; }
		  if(isrc==(nsrc+17)){ pfjetAK8jesdn_RelativeBal[npfjetAK8] = sup; }
		  if(isrc==(nsrc+18)){ pfjetAK8jesdn_RelativeSample[npfjetAK8] = sup; }
		  if(isrc==(nsrc+19)){ pfjetAK8jesdn_RelativeStatEC[npfjetAK8] = sup; }
		  if(isrc==(nsrc+20)){ pfjetAK8jesdn_RelativeStatFSR[npfjetAK8] = sup; }
		  if(isrc==(nsrc+21)){ pfjetAK8jesdn_SinglePionECAL[npfjetAK8] = sup; }
		  if(isrc==(nsrc+22)){ pfjetAK8jesdn_SinglePionHCAL[npfjetAK8] = sup; }
		  if(isrc==(nsrc+23)){ pfjetAK8jesdn_TimePtEta[npfjetAK8] = sup; }
		  if(isrc==(nsrc+24)){ pfjetAK8jesdn_Total[npfjetAK8] = sup; }
		}
		
      }
       
      pfjetAK8CHF[npfjetAK8] = ak8jet.chargedHadronEnergyFraction();
      pfjetAK8NHF[npfjetAK8] = ak8jet.neutralHadronEnergyFraction();
      pfjetAK8CEMF[npfjetAK8] = ak8jet.chargedEmEnergyFraction();
      pfjetAK8NEMF[npfjetAK8] = ak8jet.neutralEmEnergyFraction();
      pfjetAK8MUF[npfjetAK8] = ak8jet.muonEnergyFraction();
      pfjetAK8PHF[npfjetAK8] = ak8jet.photonEnergyFraction();
      pfjetAK8EEF[npfjetAK8] = ak8jet.electronEnergyFraction();
      pfjetAK8HFHF[npfjetAK8] = ak8jet.HFHadronEnergyFraction();
      
      pfjetAK8CHM[npfjetAK8] = ak8jet.chargedHadronMultiplicity();
      pfjetAK8NHM[npfjetAK8] = ak8jet.neutralHadronMultiplicity();
      pfjetAK8MUM[npfjetAK8] = ak8jet.muonMultiplicity();
      pfjetAK8PHM[npfjetAK8] = ak8jet.photonMultiplicity();
      pfjetAK8EEM[npfjetAK8] = ak8jet.electronMultiplicity();
      pfjetAK8HFHM[npfjetAK8] = ak8jet.HFHadronMultiplicity();
      
      pfjetAK8Chcons[npfjetAK8] = ak8jet.chargedMultiplicity();
      pfjetAK8Neucons[npfjetAK8] = ak8jet.neutralMultiplicity();
      
      JetIDVars idvars; 
      idvars.NHF = pfjetAK8NHF[npfjetAK8];
      idvars.NEMF = pfjetAK8NEMF[npfjetAK8];
      idvars.MUF = pfjetAK8MUF[npfjetAK8];
      idvars.CHF = pfjetAK8CHF[npfjetAK8];
      idvars.CEMF = pfjetAK8CEMF[npfjetAK8];
      idvars.NumConst = (pfjetAK8Chcons[npfjetAK8]+pfjetAK8Neucons[npfjetAK8]);
      idvars.NumNeutralParticle = pfjetAK8Neucons[npfjetAK8];
      idvars.CHM = pfjetAK8CHM[npfjetAK8];
      
      pfjetAK8jetID[npfjetAK8] = getJetID(idvars,"PUPPI",year,pfjetAK8eta[npfjetAK8],false,isUltraLegacy);
      pfjetAK8jetID_tightlepveto[npfjetAK8] = getJetID(idvars,"PUPPI",year,pfjetAK8eta[npfjetAK8],true,isUltraLegacy);  
     
      pfjetAK8sub1pt[npfjetAK8] = pfjetAK8sub1eta[npfjetAK8] = pfjetAK8sub1phi[npfjetAK8] = pfjetAK8sub1mass[npfjetAK8] = pfjetAK8sub1btag[npfjetAK8] = -100;              
      pfjetAK8sub2pt[npfjetAK8] = pfjetAK8sub2eta[npfjetAK8] = pfjetAK8sub2phi[npfjetAK8] = pfjetAK8sub2mass[npfjetAK8] = pfjetAK8sub2btag[npfjetAK8] = -100;                                                        
      pfjetAK8sdmass[npfjetAK8] = -100;                                                                      
      
      if(isSoftDrop){
	
		pfjetAK8tau1[npfjetAK8] = ak8jet.userFloat(tau1);
		pfjetAK8tau2[npfjetAK8] = ak8jet.userFloat(tau2);
		pfjetAK8tau3[npfjetAK8] = ak8jet.userFloat(tau3);
		
		pfjetAK8sdmass[npfjetAK8] = (ak8jet.groomedMass(subjets) > 0)? ak8jet.groomedMass(subjets) : 0;
		
		/*
		for(unsigned int ic = 0 ; ic < ak8jet.numberOfDaughters() ; ++ic) {                                            
			const pat::PackedCandidate* con = dynamic_cast<const pat::PackedCandidate*>(ak8jet.daughter(ic));    
			cout<<"part "<<ic+1<<" pt "<<con->pt()<<endl;
		}
		*/
		for(unsigned int isub=0; isub<((ak8jet.subjets(subjets)).size()); isub++){
		
			const auto ak8subjet = (ak8jet.subjets(subjets))[isub];
	    
			if(isub==0){
				pfjetAK8sub1pt[npfjetAK8] = ak8subjet->correctedP4("Uncorrected").pt();
				pfjetAK8sub1eta[npfjetAK8] = ak8subjet->eta();
				pfjetAK8sub1phi[npfjetAK8] = ak8subjet->phi();
				pfjetAK8sub1mass[npfjetAK8] = ak8subjet->correctedP4("Uncorrected").mass();	 
				pfjetAK8sub1btag[npfjetAK8] = ak8subjet->bDiscriminator("pfDeepCSVJetTags:probb")+ak8subjet->bDiscriminator("pfDeepCSVJetTags:probbb");
			}
			else if(isub==1){
				pfjetAK8sub2pt[npfjetAK8] = ak8subjet->correctedP4("Uncorrected").pt();
				pfjetAK8sub2eta[npfjetAK8] = ak8subjet->eta();
				pfjetAK8sub2phi[npfjetAK8] = ak8subjet->phi();
				pfjetAK8sub2mass[npfjetAK8] = ak8subjet->correctedP4("Uncorrected").mass();	 
				pfjetAK8sub2btag[npfjetAK8] = ak8subjet->bDiscriminator("pfDeepCSVJetTags:probb")+ak8subjet->bDiscriminator("pfDeepCSVJetTags:probbb");
			}
		}
		
	
      }//isSoftDrop
      
      npfjetAK8++;	
      if(npfjetAK8 >= njetmxAK8) { break;}
      
    }
  }
 
  npfjetAK4 = 0;
  edm::Handle<edm::View<pat::Jet>> pfjetAK4s;
  iEvent.getByToken(tok_pfjetAK4s_, pfjetAK4s);
    
  for (unsigned jet = 0; jet< pfjetAK4s->size(); jet++) {
      
	const auto &ak4jet = (*pfjetAK4s)[jet];
    TLorentzVector pfjetAK44v(ak4jet.correctedP4("Uncorrected").px(),ak4jet.correctedP4("Uncorrected").py(),ak4jet.correctedP4("Uncorrected").pz(), ak4jet.correctedP4("Uncorrected").energy());
    
    pfjetAK44v = LeptonJet_subtraction(tlvmu,ak4jet,pfjetAK44v);
    pfjetAK44v = LeptonJet_subtraction(tlvel,ak4jet,pfjetAK44v);
    
    double tmprecpt = pfjetAK44v.Pt();
   
    double total_cor =1;
    Read_JEC(total_cor,tmprecpt,pfjetAK44v.Eta(),Rho,isData,ak4jet,jecL1FastAK4,jecL2RelativeAK4,jecL3AbsoluteAK4,jecL2L3ResidualAK4);  
    pfjetAK4JEC[npfjetAK4] = total_cor;
    
    tmprecpt = pfjetAK44v.Pt();
    if(tmprecpt<minjPt) continue;
    if(abs(pfjetAK44v.Rapidity())>maxEta) continue;
      
    pfjetAK4pt[npfjetAK4] = 	tmprecpt;
    pfjetAK4eta[npfjetAK4] = 	pfjetAK44v.Eta();
    pfjetAK4y[npfjetAK4] = pfjetAK44v.Rapidity();
    pfjetAK4phi[npfjetAK4] = pfjetAK44v.Phi();
    pfjetAK4mass[npfjetAK4] = pfjetAK44v.M(); 
     
    pfjetAK4btag_DeepCSV[npfjetAK4] = ak4jet.bDiscriminator("pfDeepCSVJetTags:probb")+ak4jet.bDiscriminator("pfDeepCSVJetTags:probbb");
    pfjetAK4btag_DeepFlav[npfjetAK4] = ak4jet.bDiscriminator("pfDeepFlavourJetTags:probb") + ak4jet.bDiscriminator("pfDeepFlavourJetTags:probbb")+ak4jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
     
    if(isMC){
		
		vector<double> SFs;
      	Read_JER(mPtResoFileAK4, mPtSFFileAK4, tmprecpt, pfjetAK44v, Rho, genjetAK4s, SFs);
      	
      	pfjetAK4reso[npfjetAK4] = SFs[0];
      	pfjetAK4resoup[npfjetAK4] = SFs[1];
      	pfjetAK4resodn[npfjetAK4] = SFs[2];
		
    }//isMC
      
     // JES uncertainty //
      
    for(int isrc =0 ; isrc<njecmcmx; isrc++){
	
		double sup = 1.0 ;
	
		if((isrc>0)&&(isrc<=nsrc)){
	  
			JetCorrectionUncertainty *jecUnc = vsrc[isrc-1];
			jecUnc->setJetEta(ak4jet.eta());
			jecUnc->setJetPt(tmprecpt);
	  
			sup += jecUnc->getUncertainty(true);         
			if(isrc==1){ pfjetAK4jesup_AbsoluteStat[npfjetAK4] = sup; }
			if(isrc==2){ pfjetAK4jesup_AbsoluteScale[npfjetAK4] = sup; }
			if(isrc==3){ pfjetAK4jesup_AbsoluteMPFBias[npfjetAK4] = sup; }
			if(isrc==4){ pfjetAK4jesup_FlavorQCD[npfjetAK4] = sup; }
			if(isrc==5){ pfjetAK4jesup_Fragmentation[npfjetAK4] = sup; }
			if(isrc==6){ pfjetAK4jesup_PileUpDataMC[npfjetAK4] = sup; }
			if(isrc==7){ pfjetAK4jesup_PileUpPtBB[npfjetAK4] = sup; }
			if(isrc==8){ pfjetAK4jesup_PileUpPtEC1[npfjetAK4] = sup; }
			if(isrc==9){ pfjetAK4jesup_PileUpPtEC2[npfjetAK4] = sup; }
			if(isrc==10){ pfjetAK4jesup_PileUpPtRef[npfjetAK4] = sup; }
			if(isrc==11){ pfjetAK4jesup_RelativeFSR[npfjetAK4] = sup; }
			if(isrc==12){ pfjetAK4jesup_RelativeJEREC1[npfjetAK4] = sup; }
			if(isrc==13){ pfjetAK4jesup_RelativeJEREC2[npfjetAK4] = sup; }
			if(isrc==14){ pfjetAK4jesup_RelativePtBB[npfjetAK4] = sup; }
			if(isrc==15){ pfjetAK4jesup_RelativePtEC1[npfjetAK4] = sup; }
			if(isrc==16){ pfjetAK4jesup_RelativePtEC2[npfjetAK4] = sup; }
			if(isrc==17){ pfjetAK4jesup_RelativeBal[npfjetAK4] = sup; }
			if(isrc==18){ pfjetAK4jesup_RelativeSample[npfjetAK4] = sup; }
			if(isrc==19){ pfjetAK4jesup_RelativeStatEC[npfjetAK4] = sup; }
			if(isrc==20){ pfjetAK4jesup_RelativeStatFSR[npfjetAK4] = sup; }
			if(isrc==21){ pfjetAK4jesup_SinglePionECAL[npfjetAK4] = sup; }
			if(isrc==22){ pfjetAK4jesup_SinglePionHCAL[npfjetAK4] = sup; }
			if(isrc==23){ pfjetAK4jesup_TimePtEta[npfjetAK4] = sup; }
			if(isrc==24){ pfjetAK4jesup_Total[npfjetAK4] = sup; }
		}
	
		else if(isrc>nsrc){
	  
			JetCorrectionUncertainty *jecUnc = vsrc[isrc-1-nsrc];
		    jecUnc->setJetEta(ak4jet.eta());
		    jecUnc->setJetPt(tmprecpt);
	  
			sup -= jecUnc->getUncertainty(false);
			if(isrc==(nsrc+1)){ pfjetAK4jesdn_AbsoluteStat[npfjetAK4] = sup; }
			if(isrc==(nsrc+2)){ pfjetAK4jesdn_AbsoluteScale[npfjetAK4] = sup; }
			if(isrc==(nsrc+3)){ pfjetAK4jesdn_AbsoluteMPFBias[npfjetAK4] = sup; }
			if(isrc==(nsrc+4)){ pfjetAK4jesdn_FlavorQCD[npfjetAK4] = sup; }
			if(isrc==(nsrc+5)){ pfjetAK4jesdn_Fragmentation[npfjetAK4] = sup; }
			if(isrc==(nsrc+6)){ pfjetAK4jesdn_PileUpDataMC[npfjetAK4] = sup; }
			if(isrc==(nsrc+7)){ pfjetAK4jesdn_PileUpPtBB[npfjetAK4] = sup; }
			if(isrc==(nsrc+8)){ pfjetAK4jesdn_PileUpPtEC1[npfjetAK4] = sup; }
			if(isrc==(nsrc+9)){ pfjetAK4jesdn_PileUpPtEC2[npfjetAK4] = sup; }
			if(isrc==(nsrc+10)){ pfjetAK4jesdn_PileUpPtRef[npfjetAK4] = sup; }
			if(isrc==(nsrc+11)){ pfjetAK4jesdn_RelativeFSR[npfjetAK4] = sup; }
			if(isrc==(nsrc+12)){ pfjetAK4jesdn_RelativeJEREC1[npfjetAK4] = sup; }
			if(isrc==(nsrc+13)){ pfjetAK4jesdn_RelativeJEREC2[npfjetAK4] = sup; }
			if(isrc==(nsrc+14)){ pfjetAK4jesdn_RelativePtBB[npfjetAK4] = sup; }
			if(isrc==(nsrc+15)){ pfjetAK4jesdn_RelativePtEC1[npfjetAK4] = sup; }
			if(isrc==(nsrc+16)){ pfjetAK4jesdn_RelativePtEC2[npfjetAK4] = sup; }
			if(isrc==(nsrc+17)){ pfjetAK4jesdn_RelativeBal[npfjetAK4] = sup; }
			if(isrc==(nsrc+18)){ pfjetAK4jesdn_RelativeSample[npfjetAK4] = sup; }
			if(isrc==(nsrc+19)){ pfjetAK4jesdn_RelativeStatEC[npfjetAK4] = sup; }
			if(isrc==(nsrc+20)){ pfjetAK4jesdn_RelativeStatFSR[npfjetAK4] = sup; }
			if(isrc==(nsrc+21)){ pfjetAK4jesdn_SinglePionECAL[npfjetAK4] = sup; }
			if(isrc==(nsrc+22)){ pfjetAK4jesdn_SinglePionHCAL[npfjetAK4] = sup; }
			if(isrc==(nsrc+23)){ pfjetAK4jesdn_TimePtEta[npfjetAK4] = sup; }
			if(isrc==(nsrc+24)){ pfjetAK4jesdn_Total[npfjetAK4] = sup; }
			
		}
	
    }
      
    // JES uncertainty Ends //
      
    JetIDVars AK4idvars;
      
    AK4idvars.NHF = ak4jet.neutralHadronEnergyFraction();
    AK4idvars.NEMF = ak4jet.neutralEmEnergyFraction();
    AK4idvars.MUF = ak4jet.muonEnergyFraction();
    AK4idvars.CHF = ak4jet.chargedHadronEnergyFraction();
    AK4idvars.CEMF = ak4jet.chargedEmEnergyFraction();
    AK4idvars.NumConst = (ak4jet.chargedMultiplicity()+ak4jet.neutralMultiplicity());
    AK4idvars.NumNeutralParticle = ak4jet.neutralMultiplicity();
    AK4idvars.CHM = ak4jet.chargedHadronMultiplicity();
     
    pfjetAK4jetID[npfjetAK4] = getJetID(AK4idvars,"CHS",year,pfjetAK4eta[npfjetAK4],false,isUltraLegacy);
    pfjetAK4jetID_tightlepveto[npfjetAK4] = getJetID(AK4idvars,"CHS",year,pfjetAK4eta[npfjetAK4],true,isUltraLegacy);
      
    pfjetAK4hadronflav[npfjetAK4] = ak4jet.hadronFlavour();
    pfjetAK4partonflav[npfjetAK4] = ak4jet.partonFlavour();
      
    pfjetAK4qgl[npfjetAK4] = ak4jet.userFloat("QGTagger:qgLikelihood");
    pfjetAK4PUID[npfjetAK4] = ak4jet.userFloat("pileupJetId:fullDiscriminant");
     
    npfjetAK4++;	
    if(npfjetAK4 >= njetmx) { break;}
    
  }
  
  
  ntaus = 0;

  for(const auto& tau1 : iEvent.get(tok_taus_) ) {

	if (tau1.pt()<mintauPt) continue;
    if(fabs(tau1.eta())>2.3) continue;

    taupt[ntaus] = tau1.pt();
    taueta[ntaus] = tau1.eta();
    tauphi[ntaus] = tau1.phi();
    taue[ntaus] = tau1.energy();
    taucharge[ntaus] = tau1.charge();
    tauisPF[ntaus] = tau1.isPFTau();
    taudxy[ntaus] = tau1.dxy();

    if(!tau1.leadTrack().isNull()){
		tauleadtrkdxy[ntaus] = tau1.leadTrack()->dxy(vertex.position());
        tauleadtrkdz[ntaus] = tau1.leadTrack()->dz(vertex.position());
    }

    if(!tau1.leadChargedHadrCand().isNull()){

	//  taudxy[ntaus] = tau1.leadChargedHadrCand()->dxy(vertex.position());
	//  taudz[ntaus] = tau1.leadChargedHadrCand()->dz(vertex.position());

		tauleadtrkpt[ntaus] = tau1.leadChargedHadrCand()->pt();
        tauleadtrketa[ntaus] = tau1.leadChargedHadrCand()->eta();
        tauleadtrkphi[ntaus] = tau1.leadChargedHadrCand()->phi();
    }
    
    taueiso2018_raw[ntaus] = tau1.tauID("againstElectronMVA6Raw2018");
    taueiso2018[ntaus] = (0 + (int(tau1.tauID("againstElectronVLooseMVA62018"))) + (1<<(1*int(tau1.tauID("againstElectronLooseMVA62018")))) + (1<<(2*int(tau1.tauID("againstElectronMediumMVA62018")))) + (1<<(3*int(tau1.tauID("againstElectronTightMVA62018")))) + (1<<(4*int(tau1.tauID("againstElectronVTightMVA62018")))));
    
    taujetiso_deeptau2017v2p1_raw[ntaus] = tau1.tauID("byDeepTau2017v2p1VSjetraw");
    taujetiso_deeptau2017v2p1[ntaus] = (0 + (int(tau1.tauID("byVVVLooseDeepTau2017v2p1VSjet"))) + (1<<(1*int(tau1.tauID("byVVLooseDeepTau2017v2p1VSjet")))) + (1<<(2*int(tau1.tauID("byVLooseDeepTau2017v2p1VSjet")))) + (1<<(3*int(tau1.tauID("byLooseDeepTau2017v2p1VSjet")))) + (1<<(4*int(tau1.tauID("byMediumDeepTau2017v2p1VSjet")))) + (1<<(5*int(tau1.tauID("byTightDeepTau2017v2p1VSjet")))) + (1<<(6*int(tau1.tauID("byVTightDeepTau2017v2p1VSjet")))) + (1<<(7*int(tau1.tauID("byVVTightDeepTau2017v2p1VSjet")))) );

    taueiso_deeptau2017v2p1_raw[ntaus] = tau1.tauID("byDeepTau2017v2p1VSeraw");
    taueiso_deeptau2017v2p1[ntaus] = (0 + (int(tau1.tauID("byVVVLooseDeepTau2017v2p1VSe"))) + (1<<(1*int(tau1.tauID("byVVLooseDeepTau2017v2p1VSe")))) + (1<<(2*int(tau1.tauID("byVLooseDeepTau2017v2p1VSe")))) + (1<<(3*int(tau1.tauID("byLooseDeepTau2017v2p1VSe")))) + (1<<(4*int(tau1.tauID("byMediumDeepTau2017v2p1VSe")))) + (1<<(5*int(tau1.tauID("byTightDeepTau2017v2p1VSe")))) + (1<<(6*int(tau1.tauID("byVTightDeepTau2017v2p1VSe")))) + (1<<(7*int(tau1.tauID("byVVTightDeepTau2017v2p1VSe")))) );

    taumuiso_deeptau2017v2p1_raw[ntaus] = tau1.tauID("byDeepTau2017v2p1VSmuraw");
    taumuiso_deeptau2017v2p1[ntaus] = (0 + (int(tau1.tauID("byVLooseDeepTau2017v2p1VSmu"))) + (1<<(1*int(tau1.tauID("byLooseDeepTau2017v2p1VSmu")))) + (1<<(2*int(tau1.tauID("byMediumDeepTau2017v2p1VSmu")))) + (1<<(3*int(tau1.tauID("byTightDeepTau2017v2p1VSmu")))) );
    
    taudecayMode[ntaus] = tau1.decayMode();
    taudecayModeinding[ntaus] = tau1.tauID("decayModeFinding");
    taudecayModeindingNewDMs[ntaus] = tau1.tauID("decayModeFindingNewDMs");

    taurawiso[ntaus] = tau1.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    taurawisodR03[ntaus] = (tau1.tauID("chargedIsoPtSumdR03")+TMath::Max(0.,tau1.tauID("neutralIsoPtSumdR03")-0.072*tau1.tauID("puCorrPtSum")));
    taupuCorr[ntaus] = tau1.tauID("puCorrPtSum");

    if (++ntaus>=njetmx) break;

  }
  
  // Photons //
  
  nphotons = 0;

  edm::Handle<edm::View<pat::Photon>> photons;

  edm::Handle <edm::ValueMap <float> > mvaPhoID_FallV2_raw;
  iEvent.getByToken(tok_mvaPhoID_FallV2_raw, mvaPhoID_FallV2_raw);

  iEvent.getByToken(tok_photons_, photons);

  for(const auto& gamma1 : photons->ptrs() ) {

	if(gamma1->pt() < mingmPt) continue;

    phoe[nphotons] = gamma1->energy();
    phoeta[nphotons] = gamma1->eta();
    phophi[nphotons] = gamma1->phi();
    phoe1by9[nphotons] = gamma1->maxEnergyXtal()/max(float(1),gamma1->e3x3());
    if (gamma1->hasConversionTracks()) { phoe1by9[nphotons] *= -1; }
	phoe9by25[nphotons] = gamma1->r9();
    phohadbyem[nphotons] = gamma1->hadronicOverEm();

    phomvaid_Fall17V2_WP90[nphotons] =  gamma1->photonID(mPhoID_FallV2_WP90);
    phomvaid_Fall17V2_WP80[nphotons] = gamma1->photonID(mPhoID_FallV2_WP80);
    phomvaid_Fall17V2_raw[nphotons] = (*mvaPhoID_FallV2_raw)[gamma1];
	
	phomvaid_Spring16V1_WP90[nphotons] = gamma1->photonID(mPhoID_SpringV1_WP90);
    phomvaid_Spring16V1_WP80[nphotons] = gamma1->photonID(mPhoID_SpringV1_WP80);

    photrkiso[nphotons] = gamma1->trkSumPtSolidConeDR04();
    phoemiso[nphotons] = gamma1->ecalRecHitSumEtConeDR04();
    phohadiso[nphotons] = gamma1->hcalTowerSumEtConeDR04();
    phophoiso[nphotons] = gamma1->photonIso() ;
    phochhadiso[nphotons] = gamma1->chargedHadronIso();
    phoneuhadiso[nphotons] = gamma1->neutralHadronIso();
    phoietaieta[nphotons] = gamma1->sigmaIetaIeta();

    if (++nphotons>=njetmx) break;

  }
   
 
  //booltrg
 
  for(int jk=0; jk<nHLTmx; jk++) {
	  if(jk==0) {  hlt_IsoMu24 = booltrg[jk]; }
	  else if(jk==1) {  hlt_Mu50 = booltrg[jk]; }
	  else if(jk==2) {  hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 = booltrg[jk]; }
	  else if(jk==3) {  hlt_Ele115_CaloIdVT_GsfTrkIdT = booltrg[jk]; }
	  else if(jk==4) {  hlt_Ele40_WPTight_Gsf = booltrg[jk]; }
	  else if(jk==5) {  hlt_Mu37_Ele27_CaloIdL_MW = booltrg[jk]; }
	  else if(jk==6) {  hlt_Mu27_Ele37_CaloIdL_MW = booltrg[jk]; }
	  else if(jk==7) {  hlt_Mu37_TkMu27 = booltrg[jk]; }
	  else if(jk==8) {  hlt_DoubleEle25_CaloIdL_MW = booltrg[jk]; }
	  else if(jk==9) {  hlt_AK8PFJet500 = booltrg[jk]; }
	  else if(jk==10) {  hlt_PFJet500 = booltrg[jk]; }
	  else if(jk==11) {  hlt_HT1050 = booltrg[jk]; }
	  else if(jk==12) {  hlt_Photon200 = booltrg[jk]; }
  }
    
  //cout<<"done!"<<endl;
  T1->Fill();
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
Leptop::beginJob()
{
  
  Nevt = 0;
  
  ////JEC /////
  
  L1FastAK4       = new JetCorrectorParameters(mJECL1FastFileAK4.c_str());
  L2RelativeAK4   = new JetCorrectorParameters(mJECL2RelativeFileAK4.c_str());
  L3AbsoluteAK4   = new JetCorrectorParameters(mJECL3AbsoluteFileAK4.c_str());
  L2L3ResidualAK4 = new JetCorrectorParameters(mJECL2L3ResidualFileAK4.c_str());
  
  vecL1FastAK4.push_back(*L1FastAK4);
  vecL2RelativeAK4.push_back(*L2RelativeAK4);
  vecL3AbsoluteAK4.push_back(*L3AbsoluteAK4);
  vecL2L3ResidualAK4.push_back(*L2L3ResidualAK4);
  
  jecL1FastAK4       = new FactorizedJetCorrector(vecL1FastAK4);
  jecL2RelativeAK4   = new FactorizedJetCorrector(vecL2RelativeAK4);
  jecL3AbsoluteAK4   = new FactorizedJetCorrector(vecL3AbsoluteAK4);
  jecL2L3ResidualAK4 = new FactorizedJetCorrector(vecL2L3ResidualAK4);
  
  L1FastAK8       = new JetCorrectorParameters(mJECL1FastFileAK8.c_str());
  L2RelativeAK8   = new JetCorrectorParameters(mJECL2RelativeFileAK8.c_str());
  L3AbsoluteAK8   = new JetCorrectorParameters(mJECL3AbsoluteFileAK8.c_str());
  L2L3ResidualAK8 = new JetCorrectorParameters(mJECL2L3ResidualFileAK8.c_str());
  
  vecL1FastAK8.push_back(*L1FastAK8);
  vecL2RelativeAK8.push_back(*L2RelativeAK8);
  vecL3AbsoluteAK8.push_back(*L3AbsoluteAK8);
  vecL2L3ResidualAK8.push_back(*L2L3ResidualAK8);
  
  jecL1FastAK8       = new FactorizedJetCorrector(vecL1FastAK8);
  jecL2RelativeAK8   = new FactorizedJetCorrector(vecL2RelativeAK8);
  jecL3AbsoluteAK8   = new FactorizedJetCorrector(vecL3AbsoluteAK8);
  jecL2L3ResidualAK8 = new FactorizedJetCorrector(vecL2L3ResidualAK8);
  
  for (int isrc = 0; isrc < nsrc; isrc++) {
    const char *name = jecsrcnames[isrc];
    JetCorrectorParameters *pAK4 = new JetCorrectorParameters(mJECUncFileAK4.c_str(), name) ;
    JetCorrectionUncertainty *uncAK4 = new JetCorrectionUncertainty(*pAK4);
    vsrc.push_back(uncAK4);
    JetCorrectorParameters *pAK8 = new JetCorrectorParameters(mJECUncFileAK8.c_str(), name) ;
    JetCorrectionUncertainty *uncAK8 = new JetCorrectionUncertainty(*pAK8);
    vsrcAK8.push_back(uncAK8);
  }
  
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Leptop::endJob() 
{
  theFile->cd();
  theFile->Write();
  theFile->Close();
}

// ------------ method called when starting to processes a run  ------------
void 
Leptop::beginRun(edm::Run const& iRun, edm::EventSetup const& pset)
{	
  bool changed(true);
  
  hltPrescaleProvider_.init(iRun,pset,theHLTTag,changed);
  HLTConfigProvider const&  hltConfig_ = hltPrescaleProvider_.hltConfigProvider();
  //hltConfig_.dump("Triggers");
  //hltConfig_.dump("PrescaleTable");
 
}

// ------------ method called when ending the processing of a run  ------------
void 
Leptop::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Leptop::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Leptop::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Leptop::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Leptop);

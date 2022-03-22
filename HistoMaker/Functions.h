#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <math.h>
#include "TLorentzVector.h"
#include <TProofOutputFile.h>
#include <TProofServ.h>

#include <string>
#include <cmath>
#include <cassert>
#include <sstream>
#include <string>
#include "TFileCollection.h"
#include "THashList.h"
#include "TBenchmark.h"

#include <iostream>
#include <fstream>

#include "TSystem.h"

//#include "boost/config.hpp"
//#include "boost/lexical_cast.hpp"

#include "TMatrixDBase.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"

#include <TRandom3.h>

#include <iostream>
#include <fstream>

int Sign(float xx)
{
if(xx<0) { return -1; }
else { return +1; }	
}

int* decToBinary(int n)
{
	const int length = 15;
    static int binaryNum[length];
    int i = 0;
    while (n > 0) {
        binaryNum[i] = n % 2;
        n = n / 2;
        i++;
    }
    int sum=0;
    
    return binaryNum;
    /*
    if(i>10) { i = 10; }
    // printing binary array in reverse order
    for (int j = i - 1; j >= 0; j--) { sum += binaryNum[j]*pow(10,j); }
    return sum;
    */ 
} 

int getbinid(double val, int nbmx, float* array) {
  if (val<array[0]) return -2;
  for (int ix=0; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return -3;
}

int getbinid(double val, int nbmx, double* array) {
  if (val<array[0]) return -2;
  for (int ix=0; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return -3;
}

double theta_to_eta(double theta) { return -log(tan(theta/2.)); }

double eta_to_theta(double eta){
return(2*atan(exp(-2*eta)));
}

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

double delta2R_vec(TLorentzVector vec1, TLorentzVector vec2) {
  return sqrt(pow(vec1.Eta() - vec2.Eta(),2) +pow(PhiInRange(vec1.Phi() - vec2.Phi()),2));
}

double EW_toppt_cor(double pt){
return (exp(-1.08872-(pt*0.011998)) + 0.895139);
}
/*
double SF_TOP(double alpha, double beta, double pt0, double pt1)
{
	double sfwt = sqrt(exp(alpha-beta*pt0) * exp(alpha-beta*pt1));
	return sfwt;
}
*/

float Calc_MT(const TLorentzVector t1, const TLorentzVector t2)
{

float mT2 = 0;

TLorentzVector vec1, vec2;
vec1 = t1;
vec2 = t2;

vec1.SetPz(0);
vec2.SetPz(0);

mT2 = (vec1+vec2).M();

return mT2;
	
}

TLorentzVector neutrino_mom(TLorentzVector vec_lep, float MET_pt, float MET_phi, double seed){

	float W_mass = 80.4;

	TLorentzVector pnu;
	
	if(vec_lep.E()<1.e-6) {pnu.SetPtEtaPhiM(0,-100,-100,0);}
	
	else{
	
		float mT2 = 2*vec_lep.Pt()*MET_pt*(1-cos(PhiInRange(vec_lep.Phi()-MET_phi)));
	
		float Delta2 = (W_mass*W_mass - mT2)*1./(2*MET_pt*vec_lep.Pt());
	
		if(Delta2>=0){
	
			float nueta;
			nueta = (seed>=0.5)?(vec_lep.Eta() + fabs(acosh(1+(Delta2)))):(vec_lep.Eta() - fabs(acosh(1+(Delta2))));
	
			pnu.SetPtEtaPhiM(MET_pt,nueta,MET_phi,0);
	
		}
		else{
			//pnu.SetPtEtaPhiM(0,-100,-100,0);
			pnu.SetPtEtaPhiM(MET_pt,vec_lep.Eta(),MET_phi,0);
			}
	}
	
	return pnu;
}


TLorentzVector neutrino_mom_fromH(TLorentzVector vec_X, float MET_pt, float MET_phi, double seed){

	float H_mass = 125;

	TLorentzVector pnu;
	
	if(vec_X.E()<1.e-6) {pnu.SetPtEtaPhiM(0,-100,-100,0);}
	
	else{
	
		float Delta2 = (H_mass*H_mass - vec_X.M()*vec_X.M() + 2*vec_X.Pt()*MET_pt*cos(PhiInRange(vec_X.Phi()-MET_phi)))*1./(2*MET_pt*vec_X.Mt());
		if(Delta2>=1){
			float nueta;
			nueta = (seed>=0.5)?(vec_X.Rapidity() + acosh(Delta2)):(vec_X.Rapidity() - acosh(Delta2));
			pnu.SetPtEtaPhiM(MET_pt,nueta,MET_phi,0);
		}
		else{
			//pnu.SetPtEtaPhiM(0,-100,-100,0);
			float nueta;
			nueta = (vec_X.Rapidity());
			pnu.SetPtEtaPhiM(MET_pt,nueta,MET_phi,0);
			}
	}
	
	return pnu;
}

bool Muon_Tight_ID(bool muonisGL,bool muonisPF, float muonchi, float muonhit, float muonmst,
                                  float muontrkvtx, float muondz, float muonpixhit, float muontrklay){
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
//https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Tight_Muon

bool Muon_Iso_ID(float muonpfiso)
{
bool isoid = false;
if(muonpfiso<0.15) { isoid = true; } //SR
//if(muonpfiso>0.15) { isoid = true; }  // CR
return isoid;
}

void Normalize_h(TH1D *hin,bool normalize=false, bool dividebywidth=false){

if(dividebywidth){
	for(int bn=0; bn<hin->GetNbinsX(); bn++){
        hin->SetBinContent(bn+1,hin->GetBinContent(bn+1)*1./hin->GetBinWidth(bn+1));
	}
}

if(normalize==1){
    hin->Scale(1./hin->Integral());
}

}

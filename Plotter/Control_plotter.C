#include "My_Style_Suman.C"
#include <math.h>
#include "Functions.h"

struct variable_plot
{
	string name;
	string title;
	bool plot_sys;
	bool normalize;
};

void SetLimits(TH1D *hin, int maxbin, float max, float min, int Logy){
	
if(Logy==0){
	if(min<0){
		hin->SetMaximum(2.*max);
		hin->SetMinimum(1.5*min);
		}
	else{
		hin->SetMaximum(2.*max);
		hin->SetMinimum(2.5*min);
		hin->SetMinimum(0.0);
		} 
	}
else{
	if(maxbin<int(0.3*hin->GetNbinsX())){
		hin->SetMaximum(1000*max);
	}else{
		hin->SetMaximum(1000*max);
		}
	hin->SetMinimum(TMath::Max(float(1.e-2),min));
	}
}

float leg_x1 = 0.55;//75;
float leg_x2 = 0.885;
float leg_y1 = 0.75;
float leg_y2 = 0.915;

vector<pair<const char *,const char *>> samples;
vector <pair<const char *,int>> bkgs;
vector <pair<const char *,int>> sigs;

const char *sysnames[] = {
	 "JES_AbsoluteStat", "JES_AbsoluteScale","JES_AbsoluteMPFBias", 
	 "JES_FlavorQCD", "JES_Fragmentation", 
	 "JES_PileUpDataMC",  "JES_PileUpPtBB", "JES_PileUpPtEC1", "JES_PileUpPtEC2", 
	 "JES_PileUpPtRef",
	 "JES_RelativeFSR", "JES_RelativeJEREC1", "JES_RelativeJEREC2", 
	 "JES_RelativePtBB", "JES_RelativePtEC1", "JES_RelativePtEC2", 
	 "JES_RelativeBal", "JES_RelativeSample", "JES_RelativeStatEC", "JES_RelativeStatFSR", 
	 "JES_SinglePionECAL", "JES_SinglePionHCAL","JES_TimePtEta",
	 "JER",
	 "PU","LeptonSF","LeptonSF2","PNbbSF"};
	 //,"BTG"};//,"Prefire"};

int nsys = sizeof(sysnames)/sizeof(sysnames[0]);	

void output_err(double b0, double berr, double b1, double b2, double* c1, double* c2){

double a1,a2;

a1 = (b1-b0)*1./b0;
a2 = (b2-b0)*1./b0;

berr = 0.*berr*1./b0;
/*	
if(a1*a2>0) {
	if(a1>0){ a1 = max(a1,a2); a2 = 0; }
	else{ a1 = 0; a2 = abs(min(a1,a2));	}
}
*/
a1 = abs(a1);
a2 = abs(a2);

a1 = (abs(a1)>abs(berr))?sqrt(a1*a1-berr*berr):a1;
a2 = (abs(a2)>abs(berr))?sqrt(a2*a2-berr*berr):a2;

if(abs(a1)>1) { a1 = 1; }
if(abs(a2)>1) { a2 = 1; }
if(!(a1>=0)) { a1 = 0; }
if(!(a2>=0)) { a2 = 0; }
if(isnan(a1)) { a1 = 0; }
if(isnan(a2)) { a2 = 0; }

*c1 = a1;
*c2 = a2;
}

string signal_process = " ";			

void Plot(vector<TH1D*> hists_in, vector<vector<TH1D*>> hists_sysup, vector<vector<TH1D*>> hists_sysdn, 
		  string canvas_name, 
		  int runtag, 
		  vector<pair<const char*,int>> bkgs, vector<pair<const char*,int>> sigs, 
		  string signal_process, 
		  string output_filepath, 
		  bool plot_sys, 
		  bool show_data, 
		  bool signal_purity_plot, 
		  bool show_significance, 
		  float data_lumi,
		  bool isDL){
			  
	int nbkg = bkgs.size();
	int nsig = sigs.size();
	
	if(hists_in.size()<2) {
		cout<<"Need at least two samples!";
		return false;
		}
	
	vector<TH1D*> hists;
	for(unsigned ih=0; ih<hists_in.size(); ih++){
		hists.push_back(hists_in[ih]);
	}
	
	unsigned imax = 0;
	
	char name[1000];
	sprintf(name,"Plot_%s",(canvas_name).c_str());
	TCanvas *canv;
	if(signal_purity_plot){
		canv = tdrCanvas(name,hists[0],runtag,0,true,signal_process);
	}else{
		for(unsigned fg=0; fg<(hists.size()); fg++){
			if(hists[fg]->GetMaximum()>hists[imax]->GetMaximum()){
				imax = fg;
				}
		}
		canv = tdrDiCanvas(name,hists[0],hists[imax],runtag,0,true,signal_process);
	}
	
	TLegend *leg;
	if(signal_purity_plot){
		leg = tdrLeg(leg_x1,leg_y1-0.085,leg_x2+0.075,leg_y2-0.06,42,0.04);
	}else{
		leg = tdrLeg(leg_x1,leg_y1,leg_x2,leg_y2,42,0.04);
	}
	
	if((nbkg)>8){ leg->SetNColumns(3); }
	else if((nbkg)>4){ leg->SetNColumns(2); }
	else{ leg->SetNColumns(2); }
	
	TLegend *leg_sig;
	if(signal_purity_plot){ leg_sig = tdrLeg(leg_x1-0.05,leg_y1-0.085,leg_x2+0.075,leg_y2-0.06,42,0.025); }
	else{ leg_sig = tdrLeg(leg_x1,leg_y1-0.1,leg_x2,leg_y1,42,0.025); }
	leg_sig->SetNColumns(2);
	
	TLegend *leg_data = tdrLeg(leg_x1,leg_y1-0.15,leg_x1+0.1,leg_y1-0.1,42,0.04);
	
	TH1D *hists_mcsum;	
	TH1D *hists_mcsum_sysup[nsys];
	TH1D *hists_mcsum_sysdn[nsys];
	
	canv->cd(1);
//	gPad->SetRightMargin(0.225);
	
	char hist_name[500];
	sprintf(hist_name,"%s",hists[0]->GetName());
	if(string(hist_name).find("miniso")!=string::npos || string(hist_name).find("_pt")!=string::npos)
	{
		gPad->SetLogy(1);
	}
	
	sprintf(name,"Stack_%s",hists[0]->GetName());
	THStack *h_var_stack = new THStack(name,"");
	
	for(int fg=0; fg<(hists.size()); fg++){
		
		if(fg==0){
			hists_mcsum = (TH1D*)hists[fg]->Clone();
			if(plot_sys){
				for(int isys=0; isys<nsys; isys++){
					hists_mcsum_sysup[isys] = (TH1D*)hists_sysup[fg][isys]->Clone();
					hists_mcsum_sysdn[isys] = (TH1D*)hists_sysdn[fg][isys]->Clone();
					}
				}
		}			
		else if(fg>0 && fg<nbkg){
			if(isnan(hists[fg]->Integral())) continue;
			hists_mcsum->Add(hists[fg]);
			if(plot_sys){
				for(int isys=0; isys<nsys; isys++){
					hists_mcsum_sysup[isys]->Add(hists_sysup[fg][isys]);
					hists_mcsum_sysdn[isys]->Add(hists_sysdn[fg][isys]);
				}
			}
		}
		
	}//fg
	
	float maxval=-1000, minval=1000;
	
	vector <float> err_up, err_dn;
	
	maxval = max(hists_mcsum->GetMaximum(),hists[imax]->GetMaximum());
	
	for(int bn=0; bn<hists_mcsum->GetNbinsX(); bn++){
		
		if(hists_mcsum->GetBinContent(bn+1)<minval && fabs(hists_mcsum->GetBinContent(bn+1))>1.e-12){
			minval = hists_mcsum->GetBinContent(bn+1);
			}
			
		float err_up_proxy = 0;
		float err_dn_proxy = 0;
		
		if(plot_sys){
			for(int isys=0; isys<nsys; isys++){
				double err_u, err_d;
				output_err(hists_mcsum->GetBinContent(bn+1),double(0.),hists_mcsum_sysup[isys]->GetBinContent(bn+1),hists_mcsum_sysdn[isys]->GetBinContent(bn+1),&err_u,&err_d);
				err_up_proxy += err_u*err_u;
				err_dn_proxy += err_d*err_d;
			}
		}else{
			err_up_proxy = 0.0000025;
			err_dn_proxy = 0.0000025;
			}
			
		err_up_proxy = sqrt(err_up_proxy);
		err_dn_proxy = sqrt(err_dn_proxy);

		err_up.push_back(err_up_proxy);
		err_dn.push_back(err_dn_proxy);
	}
	
	SetLimits(hists[0],hists_mcsum->GetMaximumBin(),maxval,minval,gPad->GetLogy());
	
	for(int fg=0; fg<int(hists.size()-1); fg++){
	
		SetAxisStyle(hists[fg]);
		hists[fg]->GetXaxis()->SetNdivisions(506);
		if(signal_purity_plot){
			hists[fg]->GetYaxis()->SetTitleOffset(1.25);
		}
		if(fg<nbkg){
			SetPlotStyle(hists[fg],kSolid,bkgs[fg].second,1,0,bkgs[fg].second,0,1001,bkgs[fg].second,1);
			if(!isnan(hists[fg]->Integral())){
				h_var_stack->Add(hists[fg]);
				leg->AddEntry(hists[fg],bkgs[fg].first,"lf");
			}
		}
		else{
			SetPlotStyle(hists[fg],kSolid,sigs[fg-nbkg].second,2,0,sigs[fg-nbkg].second,0,0,0,0);
			leg_sig->AddEntry(hists[fg],sigs[fg-nbkg].first,"l");
		}	  
	}//fg
    
	h_var_stack->Draw("hist:SAME");			// bkg
  
	for(int fg=nbkg; fg<int(hists.size()-1); fg++){
		hists[fg]->Draw("hist:SAME");		// signal
	}
	
	TH1D *h_data;
	
	if(show_data){   h_data = (TH1D*)hists[hists.size()-1]->Clone();  }
	else{
		h_data = (TH1D*)hists_mcsum->Clone();
		for(int bn=0; bn<h_data->GetNbinsX(); bn++){
			h_data->SetBinError(bn+1,0);
			}
	}
	
	if(show_data){
	//	gROOT->ForceStyle();
		SetPlotStyle(hists_mcsum,1,1,1,20,1,1,0,0,0);
		SetPlotStyle(h_data,1,1,1,20,1,1,0,0,0);
		
		h_data->Draw("PE SAME");	// data
		leg_data->AddEntry(h_data,"Data","pe");
	}
	
	TH1D *hrat;
	hrat = (TH1D*)h_data->Clone();
	hrat->Divide(hists_mcsum);
	
	TLine *line = new TLine(hrat->GetBinLowEdge(1),1,hrat->GetBinLowEdge(hrat->GetNbinsX()+1),1);
	line->SetLineColor(kBlack);
	line->SetLineStyle(kDashed);
	
	TGraphAsymmErrors *gr_err = new TGraphAsymmErrors();
	TGraphAsymmErrors *gr_err_rat = new TGraphAsymmErrors();
	
	for(int bn=0; bn<(hrat->GetNbinsX()); bn++){
	
		if(isnan(err_up[bn])) { err_up[bn] = 1.e-2; }
		if(isnan(err_dn[bn])) { err_dn[bn] = 1.e-2; }
		
		float binvalue = 1.e-20;
		if(!isnan(hists_mcsum->GetBinContent(bn+1)) and fabs(hists_mcsum->GetBinContent(bn+1))>1.e-15){
			binvalue = hists_mcsum->GetBinContent(bn+1);
			}
		
		gr_err->SetPoint(bn,hists_mcsum->GetBinCenter(bn+1),binvalue);
		gr_err->SetPointEYhigh(bn,err_up[bn]*binvalue);
		gr_err->SetPointEYlow(bn,err_dn[bn]*binvalue);
		gr_err->SetPointEXhigh(bn,0.5*hists_mcsum->GetBinWidth(bn+1));
		gr_err->SetPointEXlow(bn,0.5*hists_mcsum->GetBinWidth(bn+1));
		
		if(!isnan(hrat->GetBinContent(bn+1)) and fabs(hrat->GetBinContent(bn+1))>1.e-15){
			binvalue = 1;
		}
		
		gr_err_rat->SetPoint(bn,hrat->GetBinCenter(bn+1),1);
		gr_err_rat->SetPointEYhigh(bn,err_up[bn]*binvalue);
		gr_err_rat->SetPointEYlow(bn,err_dn[bn]*binvalue);
		gr_err_rat->SetPointEXhigh(bn,0.5*hrat->GetBinWidth(bn+1));
		gr_err_rat->SetPointEXlow(bn,0.5*hrat->GetBinWidth(bn+1));
	
	}
	
	gr_err->SetFillColor(kGray+3);
	gr_err->SetFillStyle(3004);
	
	gr_err_rat->SetFillColor(kGray+3);
	gr_err_rat->SetFillStyle(3004);
	
	if(show_data){
	//	gr_err->Draw("E2:Z");
	}
	
	TLatex latex;
	latex.SetNDC();	
	latex.SetTextAngle(0);
	latex.SetTextColor(kBlack);
	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.05);

	float lab_x, lab_y;
	if(signal_purity_plot) { lab_x = 0.69; lab_y = 0.85; }
	else { lab_x = 0.73; lab_y = 0.85; } 

	latex.SetTextSize(0.055);
	latex.SetTextFont(42);

	if(string(hists[0]->GetName()).find("_SR1")!=string::npos){
		latex.DrawLatex(0.265,0.85,"SR1");
	}else if(string(hists[0]->GetName()).find("_SR2")!=string::npos){
		latex.DrawLatex(0.265,0.85,"SR2");
	}else if(string(hists[0]->GetName()).find("_CR2")!=string::npos){
		latex.DrawLatex(0.265,0.85,"CR2");
	}else if(string(hists[0]->GetName()).find("_CR3")!=string::npos){
		latex.DrawLatex(0.265,0.85,"CR3");
	}else if(string(hists[0]->GetName()).find("_CR4")!=string::npos){
		latex.DrawLatex(0.265,0.85,"CR4");
	}else if(string(hists[0]->GetName()).find("_CR5")!=string::npos){
		latex.DrawLatex(0.265,0.85,"CR5");
	}else if(string(hists[0]->GetName()).find("_CR6")!=string::npos){
		latex.DrawLatex(0.265,0.85,"CR6");
	}else if(string(hists[0]->GetName()).find("_CR7")!=string::npos){
		latex.DrawLatex(0.265,0.85,"CR7");
	}else if(string(hists[0]->GetName()).find("_CR8")!=string::npos){
		latex.DrawLatex(0.265,0.85,"CR8");
	}
	
	
	latex.SetTextSize(0.045);
	
	if(string(hists[0]->GetName()).find("_Y_M_")!=string::npos){
		latex.DrawLatex(0.325,0.775,"Y: Medium");
	}else if(string(hists[0]->GetName()).find("_Y_L_")!=string::npos){
		latex.DrawLatex(0.325,0.775,"Y: Loose");
	}else if(string(hists[0]->GetName()).find("_Y_T_")!=string::npos){
		latex.DrawLatex(0.325,0.775,"Y: Tight");
	}
	
	if(string(hists[0]->GetName()).find("_W_M_")!=string::npos){
		latex.DrawLatex(0.475,0.775,"W: Medium");
	}else if(string(hists[0]->GetName()).find("_W_L_")!=string::npos){
		latex.DrawLatex(0.475,0.775,"W: Loose");
	}else if(string(hists[0]->GetName()).find("_W_T_")!=string::npos){
		latex.DrawLatex(0.475,0.775,"W: Tight");
	}
	
	if(string(hists[0]->GetName()).find("_Mu")!=string::npos){
		latex.DrawLatex(0.495,0.675,"#mu-channel");
	}else if(string(hists[0]->GetName()).find("_El")!=string::npos){
		latex.DrawLatex(0.495,0.675,"e-channel");
	}
	
	if(string(hists[0]->GetName()).find("_nb0")!=string::npos){
		latex.DrawLatex(0.305,0.675,"n_{b jet}^{out Y}=0");
	}else if(string(hists[0]->GetName()).find("_nb1")!=string::npos){
		latex.DrawLatex(0.305,0.675,"n_{b jet}^{out Y}>0");
	}
	
	latex.SetTextSize(0.0275);
	if(isDL){
		latex.DrawLatex(0.525,0.6,"#sigma(pp#rightarrowX#rightarrowY(#rightarrowbb)H(#rightarrowl#nul#nu)) = 1 fb");
	}
	else{
		latex.DrawLatex(0.525,0.6,"#sigma(pp#rightarrowX#rightarrowY(#rightarrowbb)H(#rightarrowl#nujj)) = 1 fb");
	}

	gPad->RedrawAxis();

	canv->cd(2);
//	gPad->SetRightMargin(0.225);
	gPad->SetLogy(0);
	
	vector <TH1D*> hsvsb;
	vector <float> pun_sig;
	
	for(int fg=nbkg; fg<int(hists.size()-1); fg++){
		
		TH1D *hcopy = (TH1D*)hists[fg]->Clone();
		sprintf(name,"Rat_%s",hcopy->GetName());
		hcopy->SetName(name);
		hsvsb.push_back(hcopy);
		
		pun_sig.push_back((hists[fg]->Integral()*1./(data_lumi*100))*1./(1+sqrt(hists_mcsum->GetSumOfWeights())));
	}

	float maxrat = -100;

	for(unsigned isig=0; isig<(hsvsb.size()); isig++){
				
		for(int bn=0; bn<hsvsb[isig]->GetNbinsX(); bn++){
			hsvsb[isig]->SetBinContent(bn+1,hsvsb[isig]->GetBinContent(bn+1)*1./sqrt(hists_mcsum->GetBinContent(bn+1)));
			if(isnan(hsvsb[isig]->GetBinContent(bn+1))||(hsvsb[isig]->GetBinContent(bn+1)>1.e+6)) { hsvsb[isig]->SetBinContent(bn+1,1.e-6); }
			hsvsb[isig]->SetBinError(bn+1,0);
		}
		
		if((hsvsb[isig]->GetMaximum())>maxrat) { maxrat = float(hsvsb[isig]->GetMaximum()); }
	}
	
	hrat->SetMinimum(0);
	hrat->SetMaximum(2);
	
	if(show_data)
	{	
		SetAxisStyle(hrat,true);
		hrat->GetYaxis()->SetTitle("Data / Bkg.");
		
		for(int bn=0; bn<hrat->GetNbinsX(); bn++){
			if(isnan(hrat->GetBinContent(bn+1))) {
				hrat->SetBinContent(bn+1,-100);
			}
		}
    
		hrat->Draw("pe");
		gr_err_rat->Draw("E2");
		line->Draw("sames");
        hrat->Draw("SAME:pe");
	}
	else{
	
		if(!show_significance && plot_sys){
			
			SetAxisStyle(hrat,true);
			hrat->GetYaxis()->SetTitle("Relative sys. unc.");
			
			for(int bn=0; bn<hrat->GetNbinsX(); bn++){
				hrat->SetBinContent(bn+1,1);
			}
    
			hrat->Draw("hist");
			gr_err_rat->Draw("E2");
			line->Draw("sames");
			//  hrat->Draw("SAME:p");
		}
		else
		{
			for(unsigned isig=0; isig<(hsvsb.size()); isig++){
		
				hsvsb[isig]->SetMaximum(1.25*maxrat);
		
				SetAxisStyle(hsvsb[isig],true);
				hsvsb[isig]->GetYaxis()->SetTitle("S/#sqrt{B}");
				hsvsb[isig]->SetMinimum(0);
				if(isig==0){
					hsvsb[isig]->Draw("hist");
				}else{
					hsvsb[isig]->Draw("hist:SAME");
				}
		
				latex.SetTextSize(0.075);
				latex.SetTextColor(sigs[isig].second);
				if((string(hists[0]->GetName()).find("_SR")!=string::npos || string(hists[0]->GetName()).find("_CR1")!=string::npos) && string(hists[0]->GetName()).find("_X_mass")!=string::npos)// && string(hists[0]->GetName()).find("unrolled")!=string::npos)
				{
					sprintf(name,"Punzi sig: %6.2e\n",pun_sig[isig]);
					latex.DrawLatex(0.425+0.3*(isig/3),0.775-(0.1*(isig%3)),name);
				}
			}
		}
	}
	
	sprintf(name,"%s/%s.png",(output_filepath).c_str(),canv->GetName());
	canv->SaveAs(name);
	
	hists.clear();
	hists.shrink_to_fit();
	
	canv->Close();
	gSystem->ProcessEvents();
	delete canv;
}

void Plot2D(TH2D *hist, 
		  string canvas_name, string sample_name,
		  int runtag, string signal_process, 
		  string output_filepath, float data_lumi){
	
	TString cmsText     = "CMS";
	float cmsTextFont   = 61;

	TString extraText   = "Simulation  Preliminary";
	float extraTextFont = 52;  // default is helvetica-italics

	float lumiTextSize     = 0.06;
	float lumiTextOffset   = 0.2;
	float cmsTextSize      = 0.065;
	float cmsTextOffset    = 0.1;

	TLatex cms_latex;
	TLatex cms_ex_latex;
	TLatex cms_lumi;

	cms_latex.SetNDC();
	cms_latex.SetTextFont(cmsTextFont);
	cms_latex.SetTextAlign(31);
	cms_latex.SetTextSize(cmsTextSize);
	cms_latex.SetTextAngle(0);

	cms_ex_latex.SetNDC();
	cms_ex_latex.SetTextFont(extraTextFont);
	cms_ex_latex.SetTextAlign(31);
	cms_ex_latex.SetTextSize(cmsTextSize);
	cms_ex_latex.SetTextAngle(0);
	
	cms_lumi.SetNDC();
	cms_lumi.SetTextFont(42);
	cms_lumi.SetTextAlign(31);
	cms_lumi.SetTextSize(lumiTextSize);
	cms_lumi.SetTextAngle(0);

	float lp = 0.2;
	float tp = 0.055;
	float rp = 0.1;
	
	char name[1000];
	sprintf(name,"Plot_%s",(canvas_name).c_str());
	TCanvas *canv;
	//canv = tdrCanvas(name,hist,runtag,0,true,signal_process);
	canv = new TCanvas(name,name,50,50,800,600);
	
	TLegend *leg;
	leg = tdrLeg(leg_x1,leg_y1-0.085,leg_x2+0.075,leg_y2-0.06,42,0.04);
	
	canv->SetFillColor(0);
    canv->SetBorderMode(0);
    canv->SetFrameFillStyle(0);
    canv->SetFrameBorderMode(0);
    canv->SetLeftMargin( 0.12);
    canv->SetRightMargin( 0.12 );
    canv->SetTopMargin( 0.08 );
    canv->SetBottomMargin( 0.12);

	canv->cd(1);
	gStyle->SetPalette(kInvertedDarkBodyRadiator);

	hist->SetBit(TH1::kNoStats);
	SetAxisStyle(hist);
	
	hist->Draw("COLZ");	
	leg->AddEntry(hist,sample_name.c_str(),"b");
	
	TLatex latex;
	latex.SetNDC();	
	latex.SetTextAngle(0);
	latex.SetTextColor(kBlack);
	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.05);

	float lab_x, lab_y;
	lab_x = 0.69; lab_y = 0.85; 

	latex.SetTextSize(0.055);
	latex.SetTextFont(42);

	if(string(hist->GetName()).find("_SR1")!=string::npos){
		latex.DrawLatex(0.265,0.85,"SR1");
	}else if(string(hist->GetName()).find("_SR2")!=string::npos){
		latex.DrawLatex(0.265,0.85,"SR2");
	}else if(string(hist->GetName()).find("_CR2")!=string::npos){
		latex.DrawLatex(0.265,0.85,"CR2");
	}else if(string(hist->GetName()).find("_CR3")!=string::npos){
		latex.DrawLatex(0.265,0.85,"CR3");
	}else if(string(hist->GetName()).find("_CR4")!=string::npos){
		latex.DrawLatex(0.265,0.85,"CR4");
	}else if(string(hist->GetName()).find("_CR5")!=string::npos){
		latex.DrawLatex(0.265,0.85,"CR5");
	}else if(string(hist->GetName()).find("_CR6")!=string::npos){
		latex.DrawLatex(0.265,0.85,"CR6");
	}else if(string(hist->GetName()).find("_CR7")!=string::npos){
		latex.DrawLatex(0.265,0.85,"CR7");
	}else if(string(hist->GetName()).find("_CR8")!=string::npos){
		latex.DrawLatex(0.265,0.85,"CR8");
	}
	
	latex.SetTextSize(0.045);
	
	if(string(hist->GetName()).find("_Y_M_")!=string::npos){
		latex.DrawLatex(0.325,0.775,"Y: Medium");
	}else if(string(hist->GetName()).find("_Y_L_")!=string::npos){
		latex.DrawLatex(0.325,0.775,"Y: Loose");
	}else if(string(hist->GetName()).find("_Y_T_")!=string::npos){
		latex.DrawLatex(0.325,0.775,"Y: Tight");
	}
	
	if(string(hist->GetName()).find("_W_M_")!=string::npos){
		latex.DrawLatex(0.475,0.775,"W: Medium");
	}else if(string(hist->GetName()).find("_W_L_")!=string::npos){
		latex.DrawLatex(0.475,0.775,"W: Loose");
	}else if(string(hist->GetName()).find("_W_T_")!=string::npos){
		latex.DrawLatex(0.475,0.775,"W: Tight");
	}
	
	if(string(hist->GetName()).find("_Mu")!=string::npos){
		latex.DrawLatex(0.495,0.675,"#mu-channel");
	}else if(string(hist->GetName()).find("_El")!=string::npos){
		latex.DrawLatex(0.495,0.675,"e-channel");
	}
	
	if(string(hist->GetName()).find("_nb0")!=string::npos){
		latex.DrawLatex(0.305,0.675,"n_{b jet}^{other}=0");
	}else if(string(hist->GetName()).find("_nb1")!=string::npos){
		latex.DrawLatex(0.305,0.675,"n_{b jet}^{other}>0");
	}
	
	latex.DrawLatex(0.5,0.85,sample_name.c_str());
	
	cms_latex.DrawLatex(lp+0.025,1-tp,cmsText);                                                                                                                                                           
    cms_ex_latex.DrawLatex(lp+0.5,1-tp,extraText);
    cms_lumi.DrawLatex(1-rp,1-tp,"(13 TeV)");

	gPad->RedrawAxis();
	
	canv->Update();
    TPaletteAxis *palette = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.875);
    palette->SetX2NDC(0.925);
    palette->SetY1NDC(0.12);
    palette->SetY2NDC(0.92);

    canv->Modified();
    canv->Update();

	sprintf(name,"%s/%s.png",(output_filepath).c_str(),canv->GetName());
	canv->SaveAs(name);
	
	canv->Close();
	gSystem->ProcessEvents();
	delete canv;
}


void Control_plotter(int year=2018, bool isDL=false, bool istopmerging=false, bool signal_purity_plot=false, bool test_run=false, bool run_significance=false)
{
if(istopmerging){	
	samples.push_back(make_pair("Top","_Top_fullymerged"));
	samples.push_back(make_pair("Top","_Top_semimerged"));
	samples.push_back(make_pair("Top","_Top_unmerged"));
}
else{
	samples.push_back(make_pair("TT",""));
	samples.push_back(make_pair("ST",""));
}
samples.push_back(make_pair("DYj",""));
samples.push_back(make_pair("Wj",""));
samples.push_back(make_pair("Diboson",""));
samples.push_back(make_pair("QCD",""));

if(!istopmerging && (1<0)){ //keep it off for now
if(isDL){
	samples.push_back(make_pair("NMSSM_MX_1500_MY_200_DL",""));
	samples.push_back(make_pair("NMSSM_MX_2000_MY_200_DL",""));
	samples.push_back(make_pair("NMSSM_MX_3000_MY_100_DL",""));
	samples.push_back(make_pair("NMSSM_MX_3000_MY_500_DL",""));
}
else{
	samples.push_back(make_pair("NMSSM_MX_1000_MY_100_FullSIM",""));
	samples.push_back(make_pair("NMSSM_MX_1500_MY_200_FullSIM",""));
	samples.push_back(make_pair("NMSSM_MX_2000_MY_200_FullSIM",""));
	samples.push_back(make_pair("NMSSM_MX_2400_MY_300_FullSIM",""));
	samples.push_back(make_pair("NMSSM_MX_3000_MY_100_FullSIM",""));
	samples.push_back(make_pair("NMSSM_MX_3000_MY_500_FullSIM",""));
}
}

if(istopmerging){
	bkgs.push_back(make_pair("t (bqq)",kViolet-8));
	bkgs.push_back(make_pair("t (bq/qq)",kYellow-6));
	bkgs.push_back(make_pair("t (q/b)",kMagenta-7));
}
else{
	bkgs.push_back(make_pair("t #bar{t}",kViolet-8));
	bkgs.push_back(make_pair("Single t",kYellow-6));
}
bkgs.push_back(make_pair("DY+j",kRed));
bkgs.push_back(make_pair("W+j",kGreen));
bkgs.push_back(make_pair("Diboson",kAzure+1));
bkgs.push_back(make_pair("QCD",kBlue));
if(!istopmerging && (1<0)){ // keep it off for now
if(isDL){
	sigs.push_back(make_pair("(1500, 200)",kBlack));
	sigs.push_back(make_pair("(2000, 200)",kRed+2));
	sigs.push_back(make_pair("(3000, 100)",kOrange+7));
	sigs.push_back(make_pair("(3000, 500)",kBlue-2));
}
else{
	sigs.push_back(make_pair("(1000, 100)",kBlack));
	sigs.push_back(make_pair("(1500, 200)",kRed+2));
	sigs.push_back(make_pair("(2000, 200)",kOrange+7));
	sigs.push_back(make_pair("(2400, 300)",kCyan-1));
	sigs.push_back(make_pair("(3000, 100)",kMagenta-2));
	sigs.push_back(make_pair("(3000, 500)",kBlue-2));
}
}

vector<variable_plot> vars;

if(istopmerging){
	vars.push_back({"Y_msoftdrop","m_{SD}^{Y} (GeV)",1});
	vars.push_back({"HTlep_pt","#H_{T}^{jets+leptons}"});
	vars.push_back({"ST","#S_{T}^{jets+leptons+MET}"});
	vars.push_back({"MT","#M_{T} (GeV)"});
}
else{
	if(isDL){
		vars.push_back({"l1_pt","p_{T}^{l1} (GeV)",false});
		vars.push_back({"l1_eta","#eta^{l1}"});
		vars.push_back({"l1_minisoall","Iso_{min}^{l1}"});
		vars.push_back({"l2_pt","p_{T}^{l2} (GeV)"});
		vars.push_back({"l2_eta","#eta^{l2}"});
		vars.push_back({"l2_minisoall","Iso_{min}^{l2}"});
		vars.push_back({"l1l2_mass","M_{l1l2} (GeV)"});
		vars.push_back({"l1l2_dR","#DeltaR_{l1l2}"});
		vars.push_back({"l1l2_deta","#Delta#eta_{l1l2}"});
		vars.push_back({"l1l2_dphi","#Delta#phi_{l1l2}"});
		vars.push_back({"dphi_MET_l1l2","#Delta#phi_{MET,(l1+l2)}"});
		vars.push_back({"MET_pt", "MET p_{T} (GeV)"});
		vars.push_back({"MT","M_{T} (GeV)"});
		vars.push_back({"Y_pt","p_{T}^{Y} (GeV)"});
		vars.push_back({"Y_msoftdrop","m_{SD}^{Y} (GeV)",1});
		vars.push_back({"Y_PNetMD_XbbvsQCD","ParticleNet (MD) XbbvsQCD (Y)"});
		vars.push_back({"Y_PNetMD_WvsQCD","ParticleNet (MD) WvsQCD (Y)"});
		vars.push_back({"Y_PNet_TvsQCD","ParticleNet TvsQCD (Y)"});
		vars.push_back({"Y_sub1_mass","m^{Y}_{subjet 1} (GeV)"});
		vars.push_back({"Y_sub2_mass","m^{Y}_{subjet 2} (GeV)"});
		vars.push_back({"Y_sub1_btag","DeepCSV (subjet 1 of Y)"});
		vars.push_back({"Y_sub2_btag","DeepCSV (subjet 2 of Y)"});
		vars.push_back({"nbjets_other","# of additional medium b-tagged jets"});
		vars.push_back({"nbjets_outY","# of medium b-tagged jets outside Y"});
		vars.push_back({"nbjets_outY_L","# of loose b-tagged jets outside Y"});
		vars.push_back({"HTlep_pt","#H_{T}^{jets+leptons} (GeV)"});
		vars.push_back({"ST","#S_{T}^{jets+leptons+MET} (GeV)"});
		vars.push_back({"unrolled_HTlep_pt","(Unrolled) #H_{T}^{jets+leptons} (GeV)",1});
		vars.push_back({"unrolled_ST","(Unrolled) #S_{T}^{jets+leptons+MET} (GeV)",1});
	}
	else{		
		vars.push_back({"l1_pt","p_{T}^{l1} (GeV)",false});
//		vars.push_back({"l1_eta","#eta^{l1}"});
		vars.push_back({"l1_minisoall","Iso_{min}^{l1}"});
		vars.push_back({"MET_pt", "MET p_{T} (GeV)"});
		vars.push_back({"MT","M_{T} (GeV)"});
		vars.push_back({"dR_lW","#Delta R (l1,W)"});
//		vars.push_back({"dphi_lW","#Delta #phi (l1,W)"});
////		vars.push_back({"Y_pt","p_{T}^{Y} (GeV)"});
		vars.push_back({"Y_msoftdrop","m_{SD}^{Y} (GeV)",1});
		vars.push_back({"Y_PNetMD_XbbvsQCD","ParticleNet (MD) XbbvsQCD (Y)"});
//		vars.push_back({"Y_PNetMD_WvsQCD","ParticleNet (MD) WvsQCD (Y)"});
////		vars.push_back({"Y_PNet_TvsQCD","ParticleNet TvsQCD (Y)"});
//		vars.push_back({"Y_sub1_mass","m^{Y}_{subjet 1} (GeV)"});
//		vars.push_back({"Y_sub2_mass","m^{Y}_{subjet 2} (GeV)"});
//		vars.push_back({"Y_sub1_btag","DeepCSV (subjet 1 of Y)"});
//		vars.push_back({"Y_sub2_btag","DeepCSV (subjet 2 of Y)"});
////		vars.push_back({"nbjets_other","# of additional medium b-tagged jets"});
		vars.push_back({"nbjets_outY","# of medium b-tagged jets outside Y"});
//		vars.push_back({"nbjets_outY_L","# of loose b-tagged jets outside Y"});
////		vars.push_back({"HTlep_pt","#H_{T}^{jets+leptons} (GeV)"});
////		vars.push_back({"ST","#S_{T}^{jets+leptons+MET} (GeV)"});
////		vars.push_back({"unrolled_HTlep_pt","(Unrolled) #H_{T}^{jets+leptons} (GeV)",1});
////		vars.push_back({"unrolled_ST","(Unrolled) #S_{T}^{jets+leptons+MET} (GeV)",1});
////		vars.push_back({"W_pt","p_{T}^{W} (GeV)"});
		vars.push_back({"W_msoftdrop","m_{SD}^{W} (GeV)"});
//		vars.push_back({"W_PNetMD_XbbvsQCD","ParticleNet (MD) XbbvsQCD (W)"});
////		vars.push_back({"W_PNetMD_WvsQCD","ParticleNet (MD) WvsQCD (W)"});
//		vars.push_back({"W_PNet_TvsQCD","ParticleNet TvsQCD (W)"});
////		vars.push_back({"W_DAK8MD_WvsQCD","DeepA8 (MD) WvsQCD (W)"});
//		vars.push_back({"W_sub1_mass","m^{W}_{subjet 1} (GeV)"});
//		vars.push_back({"W_sub2_mass","m^{W}_{subjet 2} (GeV)"});
//		vars.push_back({"W_sub1_btag","DeepCSV (subjet 1 of W)"});
//		vars.push_back({"W_sub2_btag","DeepCSV (subjet 2 of W)"});
//		vars.push_back({"H_mass","m^{H} (GeV)"});
		vars.push_back({"X_mass","m^{X} (GeV)",1});
////		vars.push_back({"unrolled_bin1_X_mass","(Unrolled) m^{X} (GeV)",1}); 
////		vars.push_back({"unrolled_bin2_X_mass","(Unrolled) m^{X} (GeV)",1});
	}
}

int nsamples = samples.size();
int nbkg = bkgs.size();
int nsig = sigs.size();
int nvar = vars.size();			

int runtag=12;
float data_lumi = 300.;
float data_scale = 1;

if(year==2016) { runtag = 9; data_lumi = 35.9; }
if(year==2017) { runtag = 8; data_lumi = 41.5; }
if(year==2018) { runtag = 10; data_lumi = 59.7; }

gStyle->SetTitleFont(42,"XYZ");
gStyle->SetLabelFont(42,"XYZ");

string filepath, output_filepath, data_file;

if(isDL){
	filepath = "/groups/hephy/cms/suman.chatterjee/XtoYH/Histograms/Dileptonic/October2022/";
	if(istopmerging){
		output_filepath = "Control_plots/DL/TopMerging";
	}
	else{
		output_filepath = "Control_plots/DL";
	}
	data_file = "MuEGamma"; 
}
else{
	//filepath = "/groups/hephy/cms/suman.chatterjee/XtoYH/Histograms/Semileptonic/November2022/";
	filepath = "/groups/hephy/cms/suman.chatterjee/XtoYH/Histograms/Semileptonic/Feb2023_3/";
	if(istopmerging){
		output_filepath = "Control_plots/SL/TopMerging";
	}
	else{
		output_filepath = "Control_plots/SL/QCD";
	}
	data_file = "MuEGammaJetHT"; 
}

float scale_factor[nsamples];

for(int is=0; is<nsamples; is++){
	scale_factor[is] = 1;
}

if(istopmerging){
	// fitting Y msoftdrop excluding CR3 //
	/*
	scale_factor[0] = 0.955;
	scale_factor[1] = 0.478;
	scale_factor[2] = 1.082;
	*/
	// fitting Y msoftdrop including CR3 //
	/*
	scale_factor[0] = 0.900;
	scale_factor[1] = 0.678;
	scale_factor[2] = 1.050;
	scale_factor[4] = 0.914;
	*/
	// fitting unrolled X mass (i.e. mX-mY) including CR3 //
	scale_factor[0] = 0.944;
	scale_factor[1] = 0.503;
	scale_factor[2] = 1.107;
	scale_factor[4] = 0.964;
}

cout<<"Starting ..."<<endl;

float ptrange = 1000;
			
const char *catnames[] = {"opt2"}; //;{"opt1","opt2"};
int ncat = sizeof(catnames)/sizeof(catnames[0]);

//const char *regnames[] = {"SR1","SR2","CR2","CR3","CR4","CR5","CR6","CR7","CR8"};
const char *regnames[] = {"CR3","CR5","QCDVR1","QCDVR2","QCDVR3"};
int nreg = sizeof(regnames)/sizeof(regnames[0]);

const char *bcats[] = {"","_nb0","_nb1"};
int nbcat = sizeof(bcats)/sizeof(bcats[0]);

const char *lepids[] = {"","_Mu","_El","_EMu"};
int nlid = sizeof(lepids)/sizeof(lepids[0]);

const char *cutnames[] = {"L","M","T"};
int ncut = sizeof(cutnames)/sizeof(cutnames[0]);

char name[1000], filename[1000];

if((nbkg+nsig)!=nsamples) {cout<<"Check number of files"; }

TFile *_file[nsamples];

for(int fg=0; fg<nsamples; fg++){
	sprintf(filename,"%s/Output_%s.root",(filepath).c_str(),samples[fg].first);
	cout<<filename<<endl;
	_file[fg] = new TFile(filename,"read");
}

sprintf(filename,"%s/Output_%s.root",(filepath).c_str(),(data_file).c_str());
TFile *_file_data = new TFile(filename,"read");

for(int icat=0; icat<ncat; icat++){
	for(int ireg=0; ireg<nreg; ireg++){
		for(int jb=0; jb<nbcat; jb++){
			for(int lc=0; lc<nlid; lc++){
				
				if(!isDL && nlid>1) { if(lc>=3) continue; }
					
				for(int ivar=0; ivar<(nvar); ivar++){
						
					vector<TH1D*> hists;
					vector< vector<TH1D*> > hists_sys_up;
					vector< vector<TH1D*> > hists_sys_dn;
				
					sprintf(name,"h_Y_T_W_T_%s_%s_%s%s%s",(vars[ivar].name).c_str(),catnames[icat],regnames[ireg],bcats[jb],lepids[lc]);

					for(int fg=0; fg<nsamples; fg++){
						
						hists_sys_up.push_back(std::vector<TH1D*>());
						hists_sys_dn.push_back(std::vector<TH1D*>());
						
						char histname[200];
						sprintf(histname,"%s%s",name,samples[fg].second);

						hists.push_back((TH1D*)_file[fg]->Get(histname));
												
						if(vars[ivar].plot_sys){
							for(int isys=0; isys<nsys; isys++){
								char name_sys[200];
								sprintf(name_sys,"%s_%s_up",histname,sysnames[isys]);
								hists_sys_up[fg].push_back((TH1D*)_file[fg]->Get(name_sys));
								sprintf(name_sys,"%s_%s_dn",histname,sysnames[isys]);
								hists_sys_dn[fg].push_back((TH1D*)_file[fg]->Get(name_sys));
							}
						}
					
						
					}//fg
										
					for(unsigned fg=0; fg<hists.size(); fg++){
						hists[fg]->Scale(scale_factor[fg]*data_lumi);
						if(vars[ivar].normalize){
							Normalize_h(hists[fg],false,true);
						}
						if(vars[ivar].plot_sys){
							for(int isys=0; isys<nsys; isys++){
								hists_sys_up[fg][isys]->Scale(scale_factor[fg]*data_lumi);
								hists_sys_dn[fg][isys]->Scale(scale_factor[fg]*data_lumi);
								if(vars[ivar].normalize){
									Normalize_h(hists_sys_up[fg][isys],false,true);
									Normalize_h(hists_sys_dn[fg][isys],false,true);
								}
							}
						}
						
					}//fg
					
					hists.push_back((TH1D*)_file_data->Get(name));
					hists[hists.size()-1]->Scale(data_scale);
					
					for(unsigned fg=0; fg<hists.size(); fg++){
						hists[fg]->GetXaxis()->SetTitle((vars[ivar].title).c_str());
						hists[fg]->GetYaxis()->SetTitle("Number of events");
					}
												
					sprintf(name,"%s_%s%s%s_%s",catnames[icat],regnames[ireg],bcats[jb],lepids[lc],(vars[ivar].name).c_str());
					
					if(((string(regnames[ireg]).find("SR1")!=string::npos) || (string(regnames[ireg]).find("SR2")!=string::npos))){
					
							Plot(hists,hists_sys_up,hists_sys_dn,
							string(name),
							runtag,
							bkgs, sigs,
							signal_process,
							output_filepath,
							vars[ivar].plot_sys,
							false,signal_purity_plot,false,
							data_lumi,isDL);
					
							if(run_significance && string(name).find("X_mass")!=string::npos){
								sprintf(name,"%s_%s%s%s_%s_Sig",catnames[icat],regnames[ireg],bcats[jb],lepids[lc],(vars[ivar].name).c_str());
								Plot(hists,hists_sys_up,hists_sys_dn,
								string(name),
								runtag,
								bkgs, sigs, 
								signal_process,
								output_filepath,
								vars[ivar].plot_sys,
								false,signal_purity_plot,true,
								data_lumi,isDL);
							}
					}
					else{
							Plot(hists,hists_sys_up,hists_sys_dn,
							string(name),
							runtag,
							bkgs, sigs,
							signal_process,
							output_filepath,
							vars[ivar].plot_sys,
							true,signal_purity_plot,false,
							data_lumi,isDL);
					}
					
					hists.clear();
					hists_sys_up.clear();
					hists_sys_dn.clear();					
					
				}//var
					
			}//lc
		}//jb
	}//reg
}//cat

samples.clear();
bkgs.clear();
sigs.clear();
vars.clear();

for(int fg=0; fg<nsamples; fg++){
	_file[fg]->Close();
}
_file_data->Close();
	
}

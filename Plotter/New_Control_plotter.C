#include "My_Style_Suman.C"
#include <math.h>

void Normalize_h(TH1D *hin,int normalize){

for(int bn=0; bn<hin->GetNbinsX(); bn++){
	hin->SetBinContent(bn+1,hin->GetBinContent(bn+1)*1./hin->GetBinWidth(bn+1));
}

if(normalize==1){
	hin->Scale(1./hin->Integral());
}
	
}

void check_zero_bin(TH1D *hin)
{
        for(int bn=0; bn<hin->GetNbinsX(); bn++){
          if((hin->GetBinContent(bn+1)) < 1.e-12){ hin->SetBinContent(bn+1,1.e-12);  }
        }
}

void check_zero_bin_2D(TH2D *hin)
{
        for(int bx=0; bx<hin->GetNbinsX(); bx++){
			for(int by=0; by<hin->GetNbinsX(); by++){
				if((hin->GetBinContent(bx+1,by+1)) < 1.e-12){ hin->SetBinContent(bx+1,by+1,1.e-12);  }
			}
        }
}

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
		/*
		if(string(hin->GetName()).find("_H_m")!=string::npos || string(hin->GetName()).find("_H_msd")!=string::npos){
			hin->SetMinimum(0.0);
			}
		*/ 
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
float leg_y1 = 0.725;
float leg_y2 = 0.885;

const int nsamples = 12;
const char *samples[nsamples] = {"TT","ST","DYj","Wj","Diboson","QCD",
								 "NMSSM_MX_1000_MY_100","NMSSM_MX_1500_MY_200","NMSSM_MX_2000_MY_200","NMSSM_MX_2400_MY_300","NMSSM_MX_3000_MY_100","NMSSM_MX_3000_MY_500"};

const int nsig = 6;
static const int col_sig[nsig] = {kBlack,kRed+2,kOrange+7,kCyan-1,kMagenta-2,kBlue-2};
const char *sig_tites[nsig] = {"1000, 100","1500, 200","2000, 200", "2400, 300", "3000, 100", "3000, 500"};

/*
const int nsamples = 7;
const char *samples[nsamples] = {"TT","ST","DYj","Wj","Diboson","QCD",
								 "NMSSM_MX_2000_MY_200_DL"};

const int nsig = 1;
static const int col_sig[nsig] = {kBlack};
const char *sig_tites[nsig] = {"NMSSM_MX_2000_MY_200"};
*/	
const int nbkg = 6;
static const int col_bkg[nbkg] = {kViolet-8,kYellow-6,kRed,kGreen,kAzure+1,kBlue};
const char *bkg_tites[nbkg] = {"t#bar{t}","Single t","DY+j","W+j","Diboson","QCD"};

const int nsys = 9;
//const char *sysnames[nsys] = {"JER"};
const char *sysnames[nsys] = {"JES","JER","PU","Prefire","LeptonSF","LeptonSF2","PNbbSF","PNWSF","BTG"};

void output_err(double b0, double berr, double b1, double b2, double* c1, double* c2){

double a1,a2;

a1 = (b1-b0)*1./b0;
a2 = (b2-b0)*1./b0;

berr = 0.*berr*1./b0;
/*	
if(a1*a2<0){
	if(a1>0){
		a1 = a1;
		a2 = abs(a2);			
	}else{
		a1 = a2;
		a2 = abs(a1);
	}
}else{
	if(a1*a2>0) {
		if(a1>0){
			a1 = max(a1,a2);
			a2 = 0;
		}
		else{
			a1 = 0;
			a2 = abs(min(a1,a2));	
		}
	}
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

string signal_process = " ";//"pp #rightarrow X #rightarrow Y (#rightarrow b #bar{b}) h (#rightarrow j j l #nu)";				

void Plot(vector<TH1D*> hists_in, vector<vector<TH1D*>> hists_sysup, vector<vector<TH1D*>> hists_sysdn, 
		  string canvas_name, 
		  int runtag, string signal_process, vector<string> sig_titles, 
		  string output_filepath, bool plot_sys, bool show_data, bool signal_purity_plot, bool show_significance, float data_lumi,
		  bool isDL=false){
	
	if(hists_in.size()<2) {
		cout<<"Need at least two samples!";
		return false;
		}
	
	vector<TH1D*> hists;
	for(unsigned ih=0; ih<hists_in.size(); ih++){
		hists.push_back(hists_in[ih]);
	}
	
	unsigned imax = 0;
	float maxcon = -100;
	
	char name[1000];
	sprintf(name,"Plot_%s",(canvas_name).c_str());
	TCanvas *canv;
	if(signal_purity_plot){
		canv = tdrCanvas(name,hists[0],runtag,0,true,signal_process);
	}else{
		for(unsigned fg=0; fg<(hists.size()); fg++){
			if(hists[fg]->GetMaximum()>maxcon){
				imax = fg;
				maxcon = hists[fg]->GetMaximum();
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
	leg->SetNColumns(2);
	
	TLegend *leg_sig;
	if(signal_purity_plot){
		leg_sig = tdrLeg(leg_x1-0.05,leg_y1-0.085,leg_x2+0.075,leg_y2-0.06,42,0.04);
	}else{
		leg_sig = tdrLeg(leg_x1-0.05,leg_y1-0.1,leg_x2,leg_y1,42,0.04);
	}
	leg_sig->SetNColumns(2);
	leg_sig->SetTextSize(0.025);
	
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
			
		if(fg>0 && fg<nbkg){
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
	
	maxval = max(float(hists_mcsum->GetMaximum()),maxcon);
	
	for(int bn=0; bn<hists_mcsum->GetNbinsX(); bn++){
		
		if(hists_mcsum->GetBinContent(bn+1)<minval && fabs(hists_mcsum->GetBinContent(bn+1))>1.e-12){
			minval = hists_mcsum->GetBinContent(bn+1);
			}
			
		float err_up_proxy = 0;
		float err_dn_proxy = 0;
		
		if(plot_sys){
			for(int isys=0; isys<nsys; isys++){
				/*
				err_up_proxy += pow((hists_mcsum_sysup[isys]->GetBinContent(bn+1) - hists_mcsum->GetBinContent(bn+1))*1./hists_mcsum->GetBinContent(bn+1),2);
				err_dn_proxy += pow((hists_mcsum_sysdn[isys]->GetBinContent(bn+1) - hists_mcsum->GetBinContent(bn+1))*1./hists_mcsum->GetBinContent(bn+1),2);
				*/
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
		if(signal_purity_plot){
			hists[fg]->GetYaxis()->SetTitleOffset(1.25);
		}
		if(fg<nbkg){
			SetPlotStyle(hists[fg],kSolid,col_bkg[fg],1,0,col_bkg[fg],0,1001,col_bkg[fg],1);
			if(!isnan(hists[fg]->Integral())){
				h_var_stack->Add(hists[fg]);
			}
			//cout<<"Fill style "<<hists[fg]->GetFillStyle()<<" color "<<col_bkg[fg]<<endl;
			leg->AddEntry(hists[fg],bkg_tites[fg],"lf");
		}
		else{
			//if(int(nbkg%2)>0 && fg==nbkg)  { leg_sig->AddEntry(hists[fg],"",""); }
			SetPlotStyle(hists[fg],kSolid,col_sig[fg-nbkg],2,0,col_sig[fg-nbkg],0,0,0,0);
			leg_sig->AddEntry(hists[fg],sig_titles[fg-nbkg].c_str(),"l");
			//if(fg!=nsamples)  { leg_sig->AddEntry(hists[fg],"",""); }
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
	/*hrat = (TH1D*)hists_mcsum->Clone();
	hrat->Divide(h_data);*/
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
		/*
		if(!isnan(h_data->GetBinContent(bn+1)) && h_data->GetBinContent(bn+1)>1.e-9){
			hrat->SetBinError(bn+1,h_data->GetBinError(bn+1)*1./h_data->GetBinContent(bn+1)); // setting error bar on ratio by hand (same as data error)
		}
		*/ 
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
		
		//hsvsb[isig]->Divide(hists_mcsum);
		
		for(int bn=0; bn<hsvsb[isig]->GetNbinsX(); bn++){
			hsvsb[isig]->SetBinContent(bn+1,hsvsb[isig]->GetBinContent(bn+1)*1./sqrt(hists_mcsum->GetBinContent(bn+1)));
			if(isnan(hsvsb[isig]->GetBinContent(bn+1))||(hsvsb[isig]->GetBinContent(bn+1)>1.e+6)) { hsvsb[isig]->SetBinContent(bn+1,1.e-6); }
			hsvsb[isig]->SetBinError(bn+1,0);
		}
		
		if((hsvsb[isig]->GetMaximum())>maxrat) { maxrat = float(hsvsb[isig]->GetMaximum()); }
	}
	
	if(!show_data)
	{
		if(!show_significance && plot_sys){
			
			SetAxisStyle(hrat,true);
			hrat->GetYaxis()->SetTitle("Relative sys. unc.");
			hrat->SetMinimum(0);
			hrat->SetMaximum(2);

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
				latex.SetTextColor(col_sig[isig]);
				if((string(hists[0]->GetName()).find("_SR")!=string::npos || string(hists[0]->GetName()).find("_CR1")!=string::npos) && string(hists[0]->GetName()).find("_X_mass")!=string::npos)// && string(hists[0]->GetName()).find("unrolled")!=string::npos)
				{
					sprintf(name,"Punzi sig: %6.2e\n",pun_sig[isig]);
					latex.DrawLatex(0.425+0.3*(isig/3),0.775-(0.1*(isig%3)),name);
				}
			}
		}
	}
	
	if(show_data)
	{	
		SetAxisStyle(hrat,true);
		hrat->GetYaxis()->SetTitle("Data / Bkg.");
		hrat->SetMinimum(0);
		hrat->SetMaximum(2);

		for(int bn=0; bn<hrat->GetNbinsX(); bn++){
			if(isnan(hrat->GetBinContent(bn+1))) {
				hrat->SetBinContent(bn+1,-100);
			}
		}
    
		hrat->Draw("p");
		gr_err_rat->Draw("E2");
		line->Draw("sames");
        hrat->Draw("SAME:p");
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


void New_Control_plotter(int year=2018, bool isDL=false, bool signal_purity_plot=false, bool test_run=false)
{

cout<<"Starting ..."<<endl;

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
	filepath = "/users/suman.chatterjee/XtoYH/CMSSW_10_6_0/src/Histograms/March_2022/Dilep/July2022/";
	output_filepath = "Control_plots/DL";
	data_file = "MuEGamma"; 
}
else{
	filepath = "/users/suman.chatterjee/XtoYH/CMSSW_10_6_0/src/Histograms/March_2022/Singlelep/July2022/";// MC_CR/IsoCut_unrolled/YMsdcut/";
	output_filepath = "Control_plots/SL";
	data_file = "MuEGammaJetHT"; 
}
	
vector<string> sig_titles;
if(isDL){
	sig_titles.push_back("(2000, 200)");
}
else{
	sig_titles.push_back("(1000, 100)");
	sig_titles.push_back("(1500, 200)");
	sig_titles.push_back("(2000, 200)");
	sig_titles.push_back("(2400, 300)");
	sig_titles.push_back("(3000, 100)");
	sig_titles.push_back("(3000, 500)");
}

float scale_factor[nsamples];

for(int is=0; is<nsamples; is++){
	scale_factor[is] = 1;
}

float ptrange = 1000;
/*
variable_plot vars[] = 
							{
							    {"l1_pt","p_{T}^{l1} (GeV)"},
							    {"l1_eta","#eta^{l1}"},
							    {"l1_minisoall","Iso_{min}^{l1}"},
							    {"l2_pt","p_{T}^{l2} (GeV)"},
							    {"l2_eta","#eta^{l2}"},
							    {"l2_minisoall","Iso_{min}^{l2}"},
							    {"l1l2_mass","M_{l1l2} (GeV)"},
							    {"l1l2_dR","#DeltaR_{l1l2}"},
							    {"l1l2_deta","#Delta#eta_{l1l2}"},
							    {"l1l2_dphi","#Delta#phi_{l1l2}"},
							    {"dphi_MET_l1l2","#Delta#phi_{MET,(l1+l2)}"},
							    {"MET_pt", "MET p_{T} (GeV)"},
							    {"MT","M_{T} (GeV)"},
								{"Y_pt","p_{T}^{Y} (GeV)"},
								{"Y_msoftdrop","m_{SD}^{Y} (GeV)",1},
								{"Y_PNetMD_XbbvsQCD","ParticleNet (MD) XbbvsQCD (Y)"},
								{"Y_PNetMD_WvsQCD","ParticleNet (MD) WvsQCD (Y)"},
								{"Y_PNet_TvsQCD","ParticleNet TvsQCD (Y)"}, 
								{"Y_sub1_mass","m^{Y}_{subjet 1} (GeV)"},
								{"Y_sub2_mass","m^{Y}_{subjet 2} (GeV)"},
								{"Y_sub1_btag","DeepCSV (subjet 1 of Y)"},
								{"Y_sub2_btag","DeepCSV (subjet 2 of Y)"},
								{"nbjets_other","# of additional medium b-tagged jets"},
								{"nbjets_outY","# of medium b-tagged jets outside Y"},
								{"nbjets_outY_L","# of loose b-tagged jets outside Y"},
								{"HTlep_pt","#H_{T}^{jets+leptons}"},
								{"ST","#S_{T}^{jets+leptons+MET}"}
							};
*/						
variable_plot vars[] = 
							{
							    {"l1_pt","p_{T}^{l} (GeV)"},
							    {"l1_eta","#eta^{l}"},
							    {"l1_minisoall","Iso_{min}^{l}"},
							    {"MET_pt", "MET p_{T} (GeV)"},
//							    {"MET_sig","MET significance"},
							    {"MT","M_{T} (GeV)"},
								{"Y_pt","p_{T}^{Y} (GeV)"},
								{"Y_msoftdrop","m_{SD}^{Y} (GeV)",1},
								{"Y_PNetMD_XbbvsQCD","ParticleNet (MD) XbbvsQCD (Y)"},
								{"Y_PNetMD_WvsQCD","ParticleNet (MD) WvsQCD (Y)"},
								{"Y_PNet_TvsQCD","ParticleNet TvsQCD (Y)"}, 
								{"Y_sub1_mass","m^{Y}_{subjet 1} (GeV)"},
								{"Y_sub2_mass","m^{Y}_{subjet 2} (GeV)"},
								{"Y_sub1_btag","DeepCSV (subjet 1 of Y)"},
								{"Y_sub2_btag","DeepCSV (subjet 2 of Y)"},
								{"nbjets_other","# of additional medium b-tagged jets"},
								{"nbjets_outY","# of medium b-tagged jets outside Y"},
								{"nbjets_outY_L","# of loose b-tagged jets outside Y"},
								{"W_pt","p_{T}^{W} (GeV)"},
								{"W_msoftdrop","m_{SD}^{W} (GeV)"},
//								{"W_PNetMD_XbbvsQCD","ParticleNet (MD) XbbvsQCD (W)"},
								{"W_PNetMD_WvsQCD","ParticleNet (MD) WvsQCD (W)"},
//								{"W_PNet_TvsQCD","ParticleNet TvsQCD (W)"},  
//								{"W_DAK8MD_WvsQCD","DeepA8 (MD) WvsQCD (W)"},
//								{"W_sub1_mass","m^{W}_{subjet 1} (GeV)"},
//								{"W_sub2_mass","m^{W}_{subjet 2} (GeV)"},
//								{"W_sub1_btag","DeepCSV (subjet 1 of W)"},
//								{"W_sub2_btag","DeepCSV (subjet 2 of W)"},
								{"H_mass","m^{H} (GeV)"},
								{"X_mass","m^{X} (GeV)",1},
								{"unrolled_bin1_X_mass","",1}, 
								{"unrolled_bin2_X_mass","",1},
								{"HTlep_pt","#H_{T}^{jets+leptons}"},
								{"ST","#S_{T}^{jets+leptons+MET}"}
//								{"MET_pt_ref",""},
//								{"Y_PNetMD_XbbvsQCD_ref",""},
//								{"W_msoftdrop_ref",""}
//								{"Y_msoftdrop_xbin","m_{SD}^{Y} (GeV)",1},
//								{"X_mass_xbin","m^{X} (GeV)",1},
//								{"unrolled_bin2_X_mass","",1}
							  } ;

int nvar = sizeof(vars)/sizeof(vars[0]);			
			
const int ncat = 1;
const char *catnames[ncat] = {"opt2"};//;{"opt1","opt2"};
const int nreg = 8;
const char *regnames[nreg] = {"SR1","SR2","CR2","CR3","CR4","CR5","CR6","CR7"};
const int nbcat = 3;
const char *bcats[nbcat] = {"","_nb0","_nb1"};
const int nlid = 4;
const char *lepids[nlid] = {"","_Mu","_El","_EMu"};

const int ncut = 3;
const char *cutnames[ncut] = {"L","M","T"};

char name[1000], filename[1000];

if((nbkg+nsig)!=nsamples) {cout<<"Check number of files"; }

for(int icat=0; icat<ncat; icat++){
	for(int ireg=0; ireg<nreg; ireg++){
		for(int jb=0; jb<nbcat; jb++){
			for(int lc=0; lc<nlid; lc++){
				
				if(!isDL && nlid>1) { if(lc>=3) continue; }
					
				for(int ivar=0; ivar<(nvar); ivar++){
						
					vector<TH1D*> hists;
					vector< vector<TH1D*> > hists_sys_up;
					vector< vector<TH1D*> > hists_sys_dn;
					
					TH1D *hist_d[nsamples];
					TH1D *hist_d_up[nsamples][nsys];
					TH1D *hist_d_dn[nsamples][nsys];
	
					sprintf(name,"h_Y_T_W_T_%s_%s_%s%s%s",(vars[ivar].name).c_str(),catnames[icat],regnames[ireg],bcats[jb],lepids[lc]);
					
					TFile *_file[nsamples];
						
					for(int fg=0; fg<nsamples; fg++){
						
						sprintf(filename,"%s/Output_%s.root",(filepath).c_str(),samples[fg]);
						_file[fg] = new TFile(filename,"read");
						hist_d[fg] = (TH1D*)_file[fg]->Get(name);
						hist_d[fg]->Scale(scale_factor[fg]*data_lumi);
						
						if(vars[ivar].plot_sys){
							for(int isys=0; isys<nsys; isys++){
								char name_sys[200];
								sprintf(name_sys,"%s_%s_up",name,sysnames[isys]);
								hist_d_up[fg][isys] = (TH1D*)_file[fg]->Get(name_sys);
								hist_d_up[fg][isys]->Scale(scale_factor[fg]*data_lumi);
								sprintf(name_sys,"%s_%s_dn",name,sysnames[isys]);
								hist_d_dn[fg][isys] = (TH1D*)_file[fg]->Get(name_sys);
								hist_d_dn[fg][isys]->Scale(scale_factor[fg]*data_lumi);
							}
						}
					
						// Set axes titles //
	
						hist_d[fg]->GetXaxis()->SetTitle((vars[ivar].title).c_str());
						hist_d[fg]->GetXaxis()->SetNdivisions(506);
						hist_d[fg]->GetYaxis()->SetTitle("Number of events");
						
						if(vars[ivar].normalize){
							Normalize_h(hist_d[fg],0);
							if(vars[ivar].plot_sys){
								for(int isys=0; isys<nsys; isys++){
									Normalize_h(hist_d_up[fg][isys],0);
									Normalize_h(hist_d_dn[fg][isys],0);
								}
							}
						}
							
						hists.push_back(hist_d[fg]);
						
						hists_sys_up.push_back(std::vector<TH1D*>());
						hists_sys_dn.push_back(std::vector<TH1D*>());
						//if(plot_sys[ivar]){
						if(vars[ivar].plot_sys){
							for(int isys=0; isys<nsys; isys++){
								hists_sys_up[fg].push_back(hist_d_up[fg][isys]);
								hists_sys_dn[fg].push_back(hist_d_dn[fg][isys]);
							}
						}
												
					}//fg
					
					TH1D *hist_data;
						
					sprintf(filename,"%s/Output_%s.root",(filepath).c_str(),(data_file).c_str());
					TFile *_file_data = new TFile(filename,"read");
					hist_data = (TH1D*)_file_data->Get(name);
					hist_data->GetXaxis()->SetTitle((vars[ivar].title).c_str());
					hist_data->GetXaxis()->SetNdivisions(506);
					hist_data->GetYaxis()->SetTitle("Number of events");
			
					hist_data->Scale(data_scale);
						
					hists.push_back(hist_data);
						
					//_file_data->Close();
						
					sprintf(name,"%s_%s%s%s_%s",catnames[icat],regnames[ireg],bcats[jb],lepids[lc],(vars[ivar].name).c_str());
					if(((string(regnames[ireg]).find("SR1")!=string::npos) || (string(regnames[ireg]).find("SR2")!=string::npos))){
							Plot(hists,hists_sys_up,hists_sys_dn,string(name),runtag,signal_process,sig_titles,output_filepath,vars[ivar].plot_sys,false,signal_purity_plot,false,data_lumi,isDL);
							if(string(name).find("X_mass")!=string::npos){
							//	sprintf(name,"%s_%s%s%s_%s_Sig",catnames[icat],regnames[ireg],bcats[jb],lepids[lc],(vars[ivar].name).c_str());
							//	Plot(hists,hists_sys_up,hists_sys_dn,string(name),runtag,signal_process,sig_titles,output_filepath,vars[ivar].plot_sys,false,signal_purity_plot,true,data_lumi);
							}
						
					}
					else{
							Plot(hists,hists_sys_up,hists_sys_dn,string(name),runtag,signal_process,sig_titles,output_filepath,vars[ivar].plot_sys,true,signal_purity_plot,false,data_lumi,isDL);
					}
					
					hists.clear();
					hists_sys_up.clear();
					hists_sys_dn.clear();
					
					for(int fg=0; fg<nsamples; fg++){
						_file[fg]->Close();
					}
					_file_data->Close();
					
				}//var
					
			}//lc
		}//jb
	}//reg
}//cat
	

}

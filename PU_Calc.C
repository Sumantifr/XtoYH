#include <iostream>
#include <TH1D.h>
#include <TFile.h>

using namespace std;
void PU_Calc(int year){

char name[1000];

sprintf(name,"data/pileup/PileupHistogram-UL%i-100bins_withVar.root",year);
TFile *file_data = new TFile(name,"read");
TH1F *h_data = (TH1F*)file_data->Get("pileup");
TH1F *h_data_plus = (TH1F*)file_data->Get("pileup_plus");
TH1F *h_data_minus = (TH1F*)file_data->Get("pileup_minus");

sprintf(name,"data/pileup/mcPileupUL%i.root",year);
TFile *file_mc = new TFile(name,"read");
TH1F *h_mc = (TH1F*)file_mc->Get("pu_mc");

h_data->Scale(1./h_data->Integral());
h_data_plus->Scale(1./h_data_plus->Integral());
h_data_minus->Scale(1./h_data_minus->Integral());
h_mc->Scale(1./h_mc->Integral());

h_data->Divide(h_mc);
h_data_plus->Divide(h_mc);
h_data_minus->Divide(h_mc);

sprintf(name,"data/pileup/RatioPileup-UL%i-100bins.root",year);
TFile *fileout = new TFile(name,"recreate");
fileout->cd();
h_data->SetName("pileup_weight");
h_data_plus->SetName("pileup_plus_weight");
h_data_minus->SetName("pileup_minus_weight");
h_data->Write();
h_data_plus->Write();
h_data_minus->Write();
fileout->Close();

/*
int bn = h_data->FindBin(npu);
cout<<"puweight "<<h_data->GetBinContent(bn+1)<<'\n';

cout<<"NBins "<<h_mc->GetNbinsX()<<endl;
cout<<"{";
for(int bn=0; bn< h_mc->GetNbinsX(); bn++){
cout<<h_data->GetBinContent(bn+1);
if(bn!=(h_mc->GetNbinsX() - 1)){cout<<",";}
}
cout<<"}\n";

cout<<"{";
for(int bn=0; bn< h_mc->GetNbinsX(); bn++){
cout<<h_data_plus->GetBinContent(bn+1);
if(bn!=(h_mc->GetNbinsX() - 1)){cout<<",";}
}
cout<<"}\n";

cout<<"{";
for(int bn=0; bn< h_mc->GetNbinsX(); bn++){
cout<<h_data_minus->GetBinContent(bn+1);
if(bn!=(h_mc->GetNbinsX() - 1)){cout<<",";}
}
cout<<"}\n";
*/
}

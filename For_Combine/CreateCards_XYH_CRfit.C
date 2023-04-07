#include <string>
#include <map>
#include <set>
#include <iostream>
#include <utility>
#include <vector>
#include <cstdlib>
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/BinByBin.h"

using namespace std;

int main(int argc, char **argv) {
	
  bool isDL = true;	
  std::istringstream(argv[1]) >> isDL;

  int lepchannel = 0;
  lepchannel = atoi(argv[2]);

  //! [part1]
  // First define the location of the "auxiliaries" directory where we can
  // source the input files containing the datacard shapes
//  string aux_shapes = string(getenv("CMSSW_BASE")) + "/src/auxiliaries/shapes/";
   string aux_shapes = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/CombineTools/bin/";

  // Create an empty CombineHarvester instance that will hold all of the
  // datacard configuration and histograms etc.
  ch::CombineHarvester cb;
  // Uncomment this next line to see a *lot* of debug information
  // cb.SetVerbosity(3);

  // Here we will just define two categories for an 8TeV analysis. Each entry in
  // the vector below specifies a bin name and corresponding bin_id.
  ch::Categories cats = {};
  //string variable = "Y_msoftdrop";
  string variable = "unrolled_bin2_X_mass";
  //string variable = "unrolled_ST";
  /* 
  cats.push_back({1, "opt2_CR2_nb1_unrolled_bin1_X_mass"});
  cats.push_back({2, "opt2_CR4_nb1_unrolled_bin1_X_mass"});
  cats.push_back({3, "opt2_CR6_nb1_unrolled_bin1_X_mass"});
  cats.push_back({4, "opt2_CR3_nb0_unrolled_bin1_X_mass"});
  */

  /*
  cats.push_back({1, "opt2_CR2_nb1_Y_msoftdrop"});
  if(!isDL){
  cats.push_back({2, "opt2_CR4_nb1_Y_msoftdrop"});
  cats.push_back({3, "opt2_CR6_nb1_Y_msoftdrop"});
  }
  cats.push_back({4, "opt2_CR3_nb0_Y_msoftdrop"});
  */
  char var_name[100];
  if(lepchannel==1){
  	cats.push_back({1, "opt2_CR2_nb1_Mu_"+variable});
  	if(!isDL){
  		cats.push_back({2, "opt2_CR4_nb1_Mu_"+variable});
		cats.push_back({3, "opt2_CR6_nb1_Mu_"+variable});
  	}
  	cats.push_back({4, "opt2_CR3_nb0_Mu_"+variable});
  }
  else if (lepchannel==2){
  	cats.push_back({1, "opt2_CR2_nb1_El_"+variable});
  	if(!isDL){
        	cats.push_back({2, "opt2_CR4_nb1_El_"+variable});
        	cats.push_back({3, "opt2_CR6_nb1_El_"+variable});
	}
  	cats.push_back({4, "opt2_CR3_nb0_El_"+variable});
  }
  else if (lepchannel==3){
  	cats.push_back({1, "opt2_CR2_nb1_EMu_"+variable});
  	if(!isDL){
        	cats.push_back({2, "opt2_CR4_nb1_EMu_"+variable});
        	cats.push_back({3, "opt2_CR6_nb1_EMu_"+variable});
  	}
  	cats.push_back({4, "opt2_CR3_nb0_EMu_"+variable});
  }
  else{
  	cats.push_back({1, "opt2_CR2_nb1_"+variable});
	 if(!isDL){
		cats.push_back({2, "opt2_CR4_nb1_"+variable});
		cats.push_back({3, "opt2_CR6_nb1_"+variable});
	 }
	cats.push_back({4, "opt2_CR3_nb0_"+variable});
  }


  // ch::Categories is just a typedef of vector<pair<int, string>>
  //! [part1]

  //! [part2]
//  vector<string> masses = ch::MassesFromRange("120-135:5");
  // Or equivalently, specify the mass points explicitly:
  vector<string> masses;
  masses.push_back("");
  //! [part2]

  //! [part3]
//  cb.AddObservations({"*"}, {"WH"}, {"13TeV_2018"}, {"bb_cpq3i"}, cats);
  //! [part3]
  //! [part4]
  //backgrounds
  vector<string> bkg_procs;
  //bkg_procs.push_back("DYj");
  bkg_procs.push_back("Diboson");
  if(!isDL){
  	bkg_procs.push_back("Wj");
  	bkg_procs.push_back("DYj");
	bkg_procs.push_back("QCD");
  }
  //signal
  vector<string> sig_procs;
  sig_procs.push_back("Top_Top_fullymerged");
  sig_procs.push_back("Top_Top_semimerged");
  sig_procs.push_back("Top_Top_unmerged");
  if(isDL){
  	sig_procs.push_back("DYj");
  }
  else{
 // 	sig_procs.push_back("Wj");
  }

  if(isDL){
	cb.AddObservations({"*"}, {"XYH"}, {"13TeV_2018"}, {"bbllnunu_TopBkg"}, cats);
	cb.AddProcesses({"*"}, {"XYH"}, {"13TeV_2018"},  {"bbllnunu_TopBkg"}, bkg_procs, cats, false);
	cb.AddProcesses(masses, {"XYH"}, {"13TeV_2018"}, {"bbllnunu_TopBkg"}, sig_procs, cats, true);
  }
  else{
	cb.AddObservations({"*"}, {"XYH"}, {"13TeV_2018"}, {"bblnujj_TopBkg"}, cats);
	cb.AddProcesses({"*"}, {"XYH"}, {"13TeV_2018"},  {"bblnujj_TopBkg"}, bkg_procs, cats, false);
	cb.AddProcesses(masses, {"XYH"}, {"13TeV_2018"}, {"bblnujj_TopBkg"}, sig_procs, cats, true);
  }
  
  //! [part4]

  //Some of the code for this is in a nested namespace, so
  // we'll make some using declarations first to simplify things a bit.
  using ch::syst::SystMap;
  using ch::syst::era;
  using ch::syst::bin_id;
  using ch::syst::process;

  //! [part5]
//  cb.cp().signals()
    cb.cp()
      .AddSyst(cb, "lumi_13TeV", "lnN", SystMap<era>::init
      ({"13TeV_2018"}, 1.018));
  //! [part5]

  //! [part6]
 
  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_AbsoluteStat", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_AbsoluteScale", "shape", SystMap<>::init(1.00));
 
  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_AbsoluteMPFBias", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_FlavorQCD", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_Fragmentation", "shape", SystMap<>::init(1.00));
  
  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_PileUpDataMC", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_PileUpPtBB", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_PileUpPtEC1", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_PileUpPtEC2", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_PileUpPtRef", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_RelativeFSR", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_RelativeJEREC1", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_RelativeJEREC2", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_RelativePtBB", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_RelativePtEC1", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_RelativePtEC2", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_RelativeBal", "shape", SystMap<>::init(1.00)); 

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_RelativeSample", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_RelativeStatEC", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_RelativeStatFSR", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_SinglePionECAL", "shape", SystMap<>::init(1.00));
 
  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_SinglePionHCAL", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JES_TimePtEta", "shape", SystMap<>::init(1.00));

//  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
//      .AddSyst(cb, "JES_Total", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "JER", "shape", SystMap<>::init(1.00));
  
  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "PU", "shape", SystMap<>::init(1.00));
      
  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "LeptonSF", "shape", SystMap<>::init(1.00));

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      .AddSyst(cb, "LeptonSF2", "shape", SystMap<>::init(1.00));
/*
  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      	.AddSyst(cb, "Prefire", "shape", SystMap<>::init(1.00));
*/
  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
      	.AddSyst(cb, "PNbbSF", "shape", SystMap<>::init(1.00));
  
  if(!isDL){
  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
        .AddSyst(cb, "PNWSF", "shape", SystMap<>::init(1.00));
  }
    /*
  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
	.AddSyst(cb, "BTG", "shape", SystMap<>::init(1.00));
  */

  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
        .AddSyst(cb, "TrigSF1", "shape", SystMap<>::init(1.00));
 
  cb.cp().process(ch::JoinStr({sig_procs, bkg_procs}))
        .AddSyst(cb, "TrigSF2", "shape", SystMap<>::init(1.00));
  //! [part6]

  //! [part7]
  string input_filename; 
  if(isDL){ input_filename = "Combine_input_XYH_CRFit_2018_DL.root"; }
  else{ input_filename = "Combine_input_XYH_CRFit_2018_SL.root"; }

  cb.cp().backgrounds().ExtractShapes(
      aux_shapes + input_filename,
      "$BIN/$PROCESS",
      "$BIN/$PROCESS_$SYSTEMATIC");
  cb.cp().signals().ExtractShapes(
      aux_shapes + input_filename,
      "$BIN/$PROCESS$MASS",
      "$BIN/$PROCESS$MASS_$SYSTEMATIC");
  //! [part7]

  //! [part8]
  
  auto bbb = ch::BinByBinFactory()
    .SetAddThreshold(1.)
    .SetFixNorm(true);
  //bbb.AddBinByBin(cb.cp().backgrounds(), cb);
  bbb.AddBinByBin(cb.cp().signals(), cb);
 
  //! [part8]
  //
  // This function modifies every entry to have a standardised bin name of
  // the form: {analysis}_{channel}_{bin_id}_{era}
  // which is commonly used in the htt analyses
  ch::SetStandardBinNames(cb);
  //! [part8]

  //! [part9]
  // First we generate a set of bin names:
  set<string> bins = cb.bin_set();
  // This method will produce a set of unique bin names by considering all
  // Observation, Process and Systematic entries in the CombineHarvester
  // instance.

  // We create the output root file that will contain all the shapes.
  // Finally we iterate through each bin,mass combination and write a datacard.
  char filename[100];
  /*
  sprintf(filename,"NMSSM_XYH_2018.input.root" ); 
  TFile output(filename, "RECREATE");
  */
  for (auto m : masses) {
	if(isDL){ sprintf(filename,"NMSSM_XYH_CRfit_2018_%s_DL.input.root",m.c_str()); }
	else { sprintf(filename,"NMSSM_XYH_CRfit_2018_%s_SL.input.root",m.c_str()); }
	TFile output(filename, "RECREATE");

	for (auto b : bins) {
		cout << ">> Writing datacard for bin: " << b << " and mass: " << m<< "\n";
		cb.cp().bin({b}).mass({m, "*"}).WriteDatacard(b + "_" + m + ".txt", output);
	}
  }
  //! [part9]

}

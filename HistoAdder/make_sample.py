import os

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--isDL',          action='store',      default=False, type=bool,             help = "WC value")
args = argParser.parse_args()

import decimal
ctx = decimal.Context()
ctx.prec = 20

def float_to_str(f):
    d1 = ctx.create_decimal(repr(f))
    return format(d1, 'f')

samples=[]

with open('mc_info.txt') as f:
    for line in f:
        x = line.split()
        samples.append(dict({"name":x[0], "sumofevents":float(x[1]), "sumofweights":float(x[2]), "xsec":float(x[3]), "kfactor":float(x[4]), "weight":float(1.)}))

for sm in samples:
    sm['weight'] = sm['xsec']*sm['kfactor']*1000./sm['sumofweights']                   

weights = [sm['weight'] for sm in samples]
names = [sm['name'] for sm in samples]

#isDL = bool(True)

pass_fraction_EGamma_RunA = 0.96
pass_fraction_EGamma_RunB = 1.0
pass_fraction_EGamma_RunC = 1.0
pass_fraction_EGamma_RunD = 1.0

pass_fraction_DoubleMuon_RunA = 0.96
pass_fraction_DoubleMuon_RunB = 1.0
pass_fraction_DoubleMuon_RunC = 1.0
pass_fraction_DoubleMuon_RunD = 1.0

pass_fraction_MuEGamma_RunA = 0.96
pass_fraction_MuEGamma_RunB = 1.0
pass_fraction_MuEGamma_RunC = 1.0
pass_fraction_MuEGamma_RunD = 1.0

pass_fraction_SingleMuon_RunA = 0.96
pass_fraction_SingleMuon_RunB = 1.0
pass_fraction_SingleMuon_RunC = 1.0
pass_fraction_SingleMuon_RunD = 1.0

pass_fraction_JetHT_RunA = 0.96
pass_fraction_JetHT_RunB = 1.0
pass_fraction_JetHT_RunC = 1.0
pass_fraction_JetHT_RunD = 1.0


sample_directory="" 
if args.isDL:
	sample_directory="/groups/hephy/cms/suman.chatterjee/XtoYH/Histograms/Dileptonic/October2022/" 
else:
	sample_directory="/groups/hephy/cms/suman.chatterjee/XtoYH/Histograms/Semileptonic/November2022/" 

os.system("./haddws.exe  "+sample_directory+"/Histogram_DYJetsToLL_M-10to50_XtoYH.root  "+sample_directory+"/Histogram_DYJetsToLL_M-50_HT-70to100_XtoYH.root "+sample_directory+"/Histogram_DYJetsToLL_M-50_HT-100To200_XtoYH.root "+sample_directory+"/Histogram_DYJetsToLL_M-50_HT-200To400_XtoYH.root "+sample_directory+"/Histogram_DYJetsToLL_M-50_HT-400To600_XtoYH.root "+sample_directory+"/Histogram_DYJetsToLL_M-50_HT-600To800_XtoYH.root "+sample_directory+"/Histogram_DYJetsToLL_M-50_HT-800To1200_XtoYH.root "+sample_directory+"/Histogram_DYJetsToLL_M-50_HT-1200To2500_XtoYH.root "+sample_directory+"/Histogram_DYJetsToLL_M-50_HT-2500ToInf_XtoYH.root "+float_to_str(weights[names.index('DYJetsToLL_M-10to50')])+" "+float_to_str(weights[names.index('DYJetsToLL_M-50_HT-70to100')])+" "+float_to_str(weights[names.index('DYJetsToLL_M-50_HT-100To200')])+" "+float_to_str(weights[names.index('DYJetsToLL_M-50_HT-200To400')])+" "+float_to_str(weights[names.index('DYJetsToLL_M-50_HT-400To600')])+" "+float_to_str(weights[names.index('DYJetsToLL_M-50_HT-600To800')])+" "+float_to_str(weights[names.index('DYJetsToLL_M-50_HT-800To1200')])+" "+float_to_str(weights[names.index('DYJetsToLL_M-50_HT-1200To2500')])+" "+float_to_str(weights[names.index('DYJetsToLL_M-50_HT-2500ToInf')]))
os.system("mv result.root "+sample_directory+"/Output_DYj.root")
print("Output_DYj.root done!")

os.system("./haddws.exe  "+sample_directory+"/Histogram_QCD_HT300to500_XtoYH.root "+sample_directory+"/Histogram_QCD_HT500to700_XtoYH.root  "+sample_directory+"/Histogram_QCD_HT700to1000_XtoYH.root "+sample_directory+"/Histogram_QCD_HT1000to1500_XtoYH.root "+sample_directory+"/Histogram_QCD_HT1500to2000_XtoYH.root "+sample_directory+"/Histogram_QCD_HT2000toInf_XtoYH.root "+float_to_str(weights[names.index('QCD_HT300to500')])+" "+float_to_str(weights[names.index('QCD_HT500to700')])+" "+float_to_str(weights[names.index('QCD_HT700to1000')])+" "+float_to_str(weights[names.index('QCD_HT1000to1500')])+" "+float_to_str(weights[names.index('QCD_HT1500to2000')])+" "+float_to_str(weights[names.index('QCD_HT2000toInf')]))
os.system("mv result.root "+sample_directory+"/Output_QCD.root")
print("Output_QCD.root done!")

os.system("./haddws.exe  "+sample_directory+"/Histogram_ST_s-channel_XtoYH.root "+sample_directory+"/Histogram_ST_t-channel_antitop_XtoYH.root "+sample_directory+"/Histogram_ST_t-channel_top_XtoYH.root "+sample_directory+"/Histogram_ST_tW_antitop_XtoYH.root "+sample_directory+"/Histogram_ST_tW_top_XtoYH.root "+float_to_str(weights[names.index('ST_s-channel')])+" "+float_to_str(weights[names.index('ST_t-channel_antitop')])+" "+float_to_str(weights[names.index('ST_t-channel_top')])+" "+float_to_str(weights[names.index('ST_tW_antitop')])+" "+float_to_str(weights[names.index('ST_tW_top')]))
os.system("mv result.root "+sample_directory+"/Output_ST.root")
print("Output_ST.root done!")

os.system("./haddws.exe  "+sample_directory+"/Histogram_TTTo2L2Nu_XtoYH.root "+sample_directory+"/Histogram_TTToHadronic_XtoYH.root "+sample_directory+"/Histogram_TTToSemiLeptonic_XtoYH.root "+float_to_str(weights[names.index('TTTo2L2Nu')])+" "+float_to_str(weights[names.index('TTToHadronic')])+" "+float_to_str(weights[names.index('TTToSemiLeptonic')]))
os.system("mv result.root "+sample_directory+"/Output_TT.root")
print("Output_TT.root done!")

os.system("hadd -fk "+sample_directory+"/Output_Top.root "+sample_directory+"/Output_TT.root "+sample_directory+"/Output_ST.root")
print("Output_Top.root done!")

os.system("./haddws.exe  "+sample_directory+"/Histogram_WJetsToLNu_HT-70To100_XtoYH.root "+sample_directory+"/Histogram_WJetsToLNu_HT-100To200_XtoYH.root "+sample_directory+"/Histogram_WJetsToLNu_HT-200To400_XtoYH.root "+sample_directory+"/Histogram_WJetsToLNu_HT-400To600_XtoYH.root "+sample_directory+"/Histogram_WJetsToLNu_HT-600To800_XtoYH.root "+sample_directory+"/Histogram_WJetsToLNu_HT-800To1200_XtoYH.root "+sample_directory+"/Histogram_WJetsToLNu_HT-1200To2500_XtoYH.root "+sample_directory+"/Histogram_WJetsToLNu_HT-2500ToInf_XtoYH.root "+float_to_str(weights[names.index('WJetsToLNu_HT-70To100')])+" "+float_to_str(weights[names.index('WJetsToLNu_HT-100To200')])+" "+float_to_str(weights[names.index('WJetsToLNu_HT-200To400')])+" "+float_to_str(weights[names.index('WJetsToLNu_HT-400To600')])+" "+float_to_str(weights[names.index('WJetsToLNu_HT-600To800')])+" "+float_to_str(weights[names.index('WJetsToLNu_HT-800To1200')])+" "+float_to_str(weights[names.index('WJetsToLNu_HT-1200To2500')])+" "+float_to_str(weights[names.index('WJetsToLNu_HT-2500ToInf')]))
os.system("mv result.root "+sample_directory+"/Output_Wj.root")
print("Output_Wj.root done!")

#os.system("./haddws.exe "+sample_directory+"/Histogram_WWTo1L1Nu2Q_XtoYH.root "+sample_directory+"/Histogram_WWTo2L2Nu_XtoYH.root "+sample_directory+"/Histogram_WWTo4Q_XtoYH.root "+float_to_str(WWTo1L1Nu2Q)+" "+float_to_str(WWTo2L2Nu)+" "+float_to_str(WWTo4Q))
if args.isDL:
	os.system("./haddws.exe "+sample_directory+"/Histogram_WWTo2L2Nu_XtoYH.root "+float_to_str(weights[names.index('WWTo2L2Nu')]))
else:
	os.system("./haddws.exe "+sample_directory+"/Histogram_WWTo1L1Nu2Q_XtoYH.root "+float_to_str(weights[names.index('WWTo1L1Nu2Q')]))
os.system("mv result.root "+sample_directory+"/Output_WW.root")
print("Output_WW.root done!")

os.system("./haddws.exe "+sample_directory+"/Histogram_WZTo1L1Nu2Q_XtoYH.root "+sample_directory+"/Histogram_WZTo2Q2L_XtoYH.root "+sample_directory+"/Histogram_WZTo2Q2Nu_XtoYH.root "+sample_directory+"/Histogram_WZTo3LNu_XtoYH.root "+float_to_str(weights[names.index('WZTo1L1Nu2Q')])+" "+float_to_str(weights[names.index('WZTo2Q2L')])+" "+float_to_str(weights[names.index('WZTo2Q2Nu')])+" "+float_to_str(weights[names.index('WZTo3LNu')]))
os.system("mv result.root "+sample_directory+"/Output_WZ.root")
print("Output_WZ.root done!")

os.system("./haddws.exe "+sample_directory+"/Histogram_ZZTo2L2Nu_XtoYH.root "+sample_directory+"/Histogram_ZZTo2Q2L_XtoYH.root "+sample_directory+"/Histogram_ZZTo2Q2Nu_XtoYH.root "+sample_directory+"/Histogram_ZZTo4L_XtoYH.root "+sample_directory+"/Histogram_ZZTo4Q_XtoYH.root "+float_to_str(weights[names.index('ZZTo2L2Nu')])+" "+float_to_str(weights[names.index('ZZTo2Q2L')])+" "+float_to_str(weights[names.index('ZZTo2Q2Nu')])+" "+float_to_str(weights[names.index('ZZTo4L')])+" "+float_to_str(weights[names.index('ZZTo4Q')]))
os.system("mv result.root "+sample_directory+"/Output_ZZ.root")
print("Output_ZZ.root done!")

os.system("hadd -fk "+sample_directory+"/Output_Diboson.root "+sample_directory+"/Output_WW.root "+sample_directory+"/Output_WZ.root "+sample_directory+"/Output_ZZ.root")
print("Output_Diboson.root done!")

'''
if not args.isDL:
	
	# FullSIM #

	os.system("./haddws.exe  "+sample_directory+"/Histogram_NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_1000_MY_100_FullSIM.root "+float_to_str(weights[names.index('NMSSM_MX_1000_MY_100_FullSIM')]))
	os.system("mv result.root "+sample_directory+"/Output_NMSSM_MX_1000_MY_100_FullSIM.root")

	os.system("./haddws.exe  "+sample_directory+"/Histogram_NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_1500_MY_200_FullSIM.root "+float_to_str(weights[names.index('NMSSM_MX_1500_MY_200_FullSIM')]))
	os.system("mv result.root "+sample_directory+"/Output_NMSSM_MX_1500_MY_200_FullSIM.root")

	os.system("./haddws.exe  "+sample_directory+"/Histogram_NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_1500_MY_250_FullSIM.root "+float_to_str(weights[names.index('NMSSM_MX_1500_MY_250_FullSIM')]))
	os.system("mv result.root "+sample_directory+"/Output_NMSSM_MX_1500_MY_250_FullSIM.root")

	os.system("./haddws.exe  "+sample_directory+"/Histogram_NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2000_MY_125_FullSIM.root "+float_to_str(weights[names.index('NMSSM_MX_2000_MY_125_FullSIM')]))
	os.system("mv result.root "+sample_directory+"/Output_NMSSM_MX_2000_MY_125_FullSIM.root")

	os.system("./haddws.exe  "+sample_directory+"/Histogram_NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2000_MY_200_FullSIM.root "+float_to_str(weights[names.index('NMSSM_MX_2000_MY_200_FullSIM')]))
	os.system("mv result.root "+sample_directory+"/Output_NMSSM_MX_2000_MY_200_FullSIM.root")

	os.system("./haddws.exe  "+sample_directory+"/Histogram_NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2400_MY_300_FullSIM.root "+float_to_str(weights[names.index('NMSSM_MX_2400_MY_300_FullSIM')]))
	os.system("mv result.root "+sample_directory+"/Output_NMSSM_MX_2400_MY_300_FullSIM.root")

	os.system("./haddws.exe  "+sample_directory+"/Histogram_NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_2600_MY_125_FullSIM.root "+float_to_str(weights[names.index('NMSSM_MX_2600_MY_125_FullSIM')]))
	os.system("mv result.root "+sample_directory+"/Output_NMSSM_MX_2600_MY_125_FullSIM.root")

	os.system("./haddws.exe  "+sample_directory+"/Histogram_NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_3000_MY_100_FullSIM.root "+float_to_str(weights[names.index('NMSSM_MX_3000_MY_100_FullSIM')]))
	os.system("mv result.root "+sample_directory+"/Output_NMSSM_MX_3000_MY_100_FullSIM.root")

	os.system("./haddws.exe  "+sample_directory+"/Histogram_NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_3000_MY_300_FullSIM.root "+float_to_str(weights[names.index('NMSSM_MX_3000_MY_300_FullSIM')]))
	os.system("mv result.root "+sample_directory+"/Output_NMSSM_MX_3000_MY_300_FullSIM.root")

	os.system("./haddws.exe  "+sample_directory+"/Histogram_NMSSM_XYH_YTobb_HToWWTo2QLNu_MX_3000_MY_500_FullSIM.root "+float_to_str(weights[names.index('NMSSM_MX_3000_MY_500_FullSIM')]))
	os.system("mv result.root "+sample_directory+"/Output_NMSSM_MX_3000_MY_500_FullSIM.root")

if args.isDL:

	os.system("./haddws.exe  "+sample_directory+"/Histogram_NMSSM_XYH_YTobb_HToWWTo2L2Nu_MX_2000_MY_200_v3.root "+float_to_str(weights[names.index('NMSSM_MX_2000_MY_200_DL')]))
	os.system("mv result.root "+sample_directory+"/Output_NMSSM_MX_2000_MY_200_DL.root")

	os.system("./haddws.exe  "+sample_directory+"/Histogram_NMSSM_XYH_YTobb_HToWWTo2L2Nu_MX_1500_MY_200_v3.root "+float_to_str(weights[names.index('NMSSM_MX_1500_MY_200_DL')]))
	os.system("mv result.root "+sample_directory+"/Output_NMSSM_MX_1500_MY_200_DL.root")

	os.system("./haddws.exe  "+sample_directory+"/Histogram_NMSSM_XYH_YTobb_HToWWTo2L2Nu_MX_3000_MY_100_v3.root "+float_to_str(weights[names.index('NMSSM_MX_3000_MY_100_DL')]))
	os.system("mv result.root "+sample_directory+"/Output_NMSSM_MX_3000_MY_100_DL.root")

	os.system("./haddws.exe  "+sample_directory+"/Histogram_NMSSM_XYH_YTobb_HToWWTo2L2Nu_MX_3000_MY_500_v3.root "+float_to_str(weights[names.index('NMSSM_MX_3000_MY_500_DL')]))
	os.system("mv result.root "+sample_directory+"/Output_NMSSM_MX_3000_MY_500_DL.root")

print("Signal done!")

'''

os.system("./haddws.exe "+sample_directory+"/Histogram_EGamma_UL2018.root "+sample_directory+"/Histogram_EGamma_UL2018A_XtoYH_Nov_2021.root "+sample_directory+"/Histogram_EGamma_UL2018B_XtoYH_Nov_2021.root " +sample_directory+"/Histogram_EGamma_UL2018C_XtoYH_Nov_2021.root " +sample_directory+"/Histogram_EGamma_UL2018D_XtoYH_Nov_2021.root "+float_to_str(pass_fraction_EGamma_RunA)+" "+float_to_str(pass_fraction_EGamma_RunB)+" "+float_to_str(pass_fraction_EGamma_RunC)+" "+float_to_str(pass_fraction_EGamma_RunD))

if args.isDL:
	os.system("hadd -fk "+sample_directory+"/Histogram_DoubleMuon_UL2018.root "+sample_directory+"/Histogram_DoubleMuon_UL2018A_XtoYH_Nov_2021.root "+sample_directory+"/Histogram_DoubleMuon_UL2018B_XtoYH_Nov_2021.root " +sample_directory+"/Histogram_DoubleMuon_UL2018C_XtoYH_Nov_2021.root " +sample_directory+"/Histogram_DoubleMuon_UL2018D_XtoYH_Nov_2021.root "+float_to_str(pass_fraction_DoubleMuon_RunA)+" "+float_to_str(pass_fraction_DoubleMuon_RunB)+" "+float_to_str(pass_fraction_DoubleMuon_RunC)+" "+float_to_str(pass_fraction_DoubleMuon_RunD))
	os.system("hadd -fk "+sample_directory+"/Histogram_MuonEG_UL2018.root "+sample_directory+"/Histogram_MuonEG_UL2018A_XtoYH_Nov_2021.root "+sample_directory+"/Histogram_MuonEG_UL2018B_XtoYH_Nov_2021.root " +sample_directory+"/Histogram_MuonEG_UL2018C_XtoYH_Nov_2021.root " +sample_directory+"/Histogram_MuonEG_UL2018D_XtoYH_Nov_2021.root "+float_to_str(pass_fraction_MuonEG_RunA)+" "+float_to_str(pass_fraction_MuonEG_RunB)+" "+float_to_str(pass_fraction_MuonEG_RunC)+" "+float_to_str(pass_fraction_MuonEG_RunD))
else:
	os.system("hadd -fk "+sample_directory+"/Histogram_SingleMuon_UL2018.root "+sample_directory+"/Histogram_SingleMuon_UL2018A_XtoYH_Nov_2021.root "+sample_directory+"/Histogram_SingleMuon_UL2018B_XtoYH_Nov_2021.root " +sample_directory+"/Histogram_SingleMuon_UL2018C_XtoYH_Nov_2021.root " +sample_directory+"/Histogram_SingleMuon_UL2018D_XtoYH_Nov_2021.root "+float_to_str(pass_fraction_SingleMuon_RunA)+" "+float_to_str(pass_fraction_SingleMuon_RunB)+" "+float_to_str(pass_fraction_SingleMuon_RunC)+" "+float_to_str(pass_fraction_SingleMuon_RunD))
	os.system("hadd -fk "+sample_directory+"/Histogram_JetHT_UL2018.root "+sample_directory+"/Histogram_JetHT_UL2018A_XtoYH_Nov_2021.root "+sample_directory+"/Histogram_JetHT_UL2018B_XtoYH_Nov_2021.root " +sample_directory+"/Histogram_JetHT_UL2018C_XtoYH_Nov_2021.root " +sample_directory+"/Histogram_JetHT_UL2018D_XtoYH_Nov_2021.root "+float_to_str(pass_fraction_JetHT_RunA)+" "+float_to_str(pass_fraction_JetHT_RunB)+" "+float_to_str(pass_fraction_JetHT_RunC)+" "+float_to_str(pass_fraction_JetHT_RunD))

if args.isDL:
	os.system("hadd -fk "+sample_directory+"/Output_MuEGamma.root "+sample_directory+"/Histogram_EGamma_UL2018.root "+sample_directory+"/Histogram_DoubleMuon_UL2018.root "+sample_directory+"/Histogram_MuonEG_UL2018.root")
	os.system("hadd -fk "+sample_directory+"/Output_MuEGammaJetHT.root "+sample_directory+"/Histogram_EGamma_UL2018.root "+sample_directory+"/Histogram_DoubleMuon_UL2018.root "+sample_directory+"/Histogram_MuonEG_UL2018.root "+sample_directory+"/Histogram_JetHT_UL2018.root")
else:
	os.system("hadd -fk "+sample_directory+"/Output_MuEGamma.root "+sample_directory+"/Histogram_EGamma_UL2018.root "+sample_directory+"/Histogram_SingleMuon_UL2018.root")
	os.system("hadd -fk "+sample_directory+"/Output_MuEGammaJetHT.root "+sample_directory+"/Histogram_EGamma_UL2018.root "+sample_directory+"/Histogram_SingleMuon_UL2018.root "+sample_directory+"/Histogram_JetHT_UL2018.root")
print("Data done!")


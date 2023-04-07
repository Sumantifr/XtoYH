combineTool.py -M Impacts -d XYH_bblnujj_TopBkg_13TeV_2018.root -m 125 --doInitialFit --robustFit 1 --rMax 1 --rMin -1
combineTool.py -M Impacts -d XYH_bblnujj_TopBkg_13TeV_2018.root -m 125 --robustFit 1 --doFits --rMax 1 --rMin -1
combineTool.py -M Impacts -d XYH_bblnujj_TopBkg_13TeV_2018.root -m 125  -o impacts_topCRfit.json
plotImpacts.py -i impacts_topCRfit.json -o impacts_topCRfit_r_fm --POI r_fm
plotImpacts.py -i impacts_topCRfit.json -o impacts_topCRfit_r_sm --POI r_sm
plotImpacts.py -i impacts_topCRfit.json -o impacts_topCRfit_r_um --POI r_um

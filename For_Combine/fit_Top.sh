lepid=$1                                                                                                                                                                                                    
echo "Creating cards"
#lepid=0
CreateCards_XYH_CRfit 0 $lepid
lepname=""
if [ $lepid == 1 ]
then
	lepname="_Mu"
elif [ $lepid == 2 ]
then
	lepname="_El"
elif [ $lepid == 3 ]
then
	lepname="_EMu"
fi
echo $lepname

echo "Combining different regions"
#combineCards.py CR2=XYH_bblnujj_TopBkg_1_13TeV_2018_.txt CR4=XYH_bblnujj_TopBkg_2_13TeV_2018_.txt CR6=XYH_bblnujj_TopBkg_3_13TeV_2018_.txt CR3=XYH_bblnujj_TopBkg_4_13TeV_2018_.txt > XYH_bblnujj_TopWBkg_13TeV_2018.txt
combineCards.py CR2=XYH_bblnujj_TopBkg_1_13TeV_2018_.txt CR4=XYH_bblnujj_TopBkg_2_13TeV_2018_.txt CR6=XYH_bblnujj_TopBkg_3_13TeV_2018_.txt > XYH_bblnujj_TopBkg_13TeV_2018${lepname}.txt
echo "Creating workspace"
#text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/Top_Top_fullymerged:r_fm[1,0,10]' --PO 'map=.*/Top_Top_semimerged:r_sm[1,0,10]' --PO 'map=.*/Top_Top_unmerged:r_um[1,0,10]'  --PO 'map=CR3/Wj:r_wj[1,0,10]' XYH_bblnujj_TopWBkg_13TeV_2018.txt -o XYH_bblnujj_TopWBkg_13TeV_2018.root
text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/Top_Top_fullymerged:r_fm[1,0,10]' --PO 'map=.*/Top_Top_semimerged:r_sm[1,0,10]' --PO 'map=.*/Top_Top_unmerged:r_um[1,0,10]'  --PO 'map=.*/Wj:1'   XYH_bblnujj_TopBkg_13TeV_2018${lepname}.txt -o XYH_bblnujj_TopBkg_13TeV_2018${lepname}.root
echo "performing best fit"
#combine -M MultiDimFit -m 125 XYH_bblnujj_TopWBkg_13TeV_2018.root --saveWorkspace --redefineSignalPOIs=r_fm,r_sm,r_um,r_wj  --setParameters r=1,r_fm=1,r_sm=1,r_um=1,r_wj=1  -n .bestfit_TopWFit 
combine -M MultiDimFit -m 125 XYH_bblnujj_TopBkg_13TeV_2018${lepname}.root --saveWorkspace --redefineSignalPOIs=r_fm,r_sm,r_um  --setParameters r=1,r_fm=1,r_sm=1,r_um=1  -n .bestfit_TopFit${lepname}
exit
echo "Generating post-fit distributions"
combine -M FitDiagnostics -m 125 higgsCombine.bestfit_TopFit${lepname}.MultiDimFit.mH125.root --verbose -1 --saveShapes --saveWithUncertainties --saveOverall  --redefineSignalPOIs=r_fm --setParameters r_fm=1 --setParameterRanges r_fm=[0,10] -n _TopFit_fm${lepname}
PostFitShapesFromWorkspace -d XYH_bblnujj_TopBkg_13TeV_2018${lepname}.txt -w XYH_bblnujj_TopBkg_13TeV_2018${lepname}.root -o fittedshapes_TopFit_fm${lepname}.root -f fitDiagnostics_TopFit_fm${lepname}.root:fit_s --postfit --sampling
combine -M FitDiagnostics -m 125 higgsCombine.bestfit_TopFit${lepname}.MultiDimFit.mH125.root --verbose -1 --saveShapes --saveWithUncertainties --saveOverall  --redefineSignalPOIs=r_sm --setParameters r_sm=1 --setParameterRanges r_sm=[0,10] -n _TopFit_sm${lepname}
PostFitShapesFromWorkspace -d XYH_bblnujj_TopBkg_13TeV_2018${lepname}.txt -w XYH_bblnujj_TopBkg_13TeV_2018${lepname}.root -o fittedshapes_TopFit_sm${lepname}.root -f fitDiagnostics_TopFit_sm${lepname}.root:fit_s --postfit --sampling
combine -M FitDiagnostics -m 125 higgsCombine.bestfit_TopFit${lepname}.MultiDimFit.mH125.root --verbose -1 --saveShapes --saveWithUncertainties --saveOverall  --redefineSignalPOIs=r_um --setParameters r_um=1 --setParameterRanges r_um=[0,10] -n _TopFit_um${lepname}
PostFitShapesFromWorkspace -d XYH_bblnujj_TopBkg_13TeV_2018${lepname}.txt -w XYH_bblnujj_TopBkg_13TeV_2018${lepname}.root -o fittedshapes_TopFit_um${lepname}.root -f fitDiagnostics_TopFit_um${lepname}.root:fit_s --postfit --sampling
#combine -M FitDiagnostics -m 125 higgsCombine.bestfit_TopWFit.MultiDimFit.mH125.root --verbose -1 --saveShapes --saveWithUncertainties --saveOverall  --redefineSignalPOIs=r_wj --setParameters r_wj=1 --setParameterRanges r_wj=[0,10] -n _TopFit_wj
#PostFitShapesFromWorkspace -d XYH_bblnujj_TopWBkg_13TeV_2018.txt -w XYH_bblnujj_TopWBkg_13TeV_2018.root -o fittedshapes_TopWFit_wj.root -f fitDiagnostics_TopWFit_wj.root:fit_s --postfit --sampling
exit

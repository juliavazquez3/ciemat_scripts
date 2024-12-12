Para hacer un fit con la datacard basta con hacer

source $VO_CMS_SW_DIR/cmsset_default.sh
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/tutorials/longexercise/

combineCards.py onebin_finalfit/datacard_final_unblinded_v2.txt > onebin_finalfit/datacard_final_unblinded_v2.text

text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose --PO 'map=.*/ttWcq:rWqq[1,0,10]' --PO 'map=.*/ttWuq:rWqq[1,0,10]'  --PO 'map=.*/stWcq:rWqq[1,0,10]' --PO 'map=.*/stWuq:rWqq[1,0,10]' onebin_finalfit/datacard_final_unblinded_v2.text -o workspace_test.root;

combine -M MultiDimFit workspace_test.root -n test --algo=singles --robustFit=1 --redefineSignalPOIs rWqq,r 


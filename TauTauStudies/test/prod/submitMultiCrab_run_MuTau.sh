#!/bin/sh

#back up
cp ../runMuTauStreamAnalyzerFullAnalysis_Recoil_cfg.py ../runMuTauStreamAnalyzerFullAnalysis_Recoil_cfg.py.mine

rm ../runMuTauStreamAnalyzerFullAnalysis_Recoil_DATA_cfg.py
sed 's/runOnMC     = True/runOnMC     = False/g' ../runMuTauStreamAnalyzerFullAnalysis_Recoil_cfg.py >> ../runMuTauStreamAnalyzerFullAnalysis_Recoil_DATA_cfg.py
echo 'multicrab -create -submit -cfg multicrab_run_MuTau_16Nov2011_DATA.cfg'
multicrab -create -submit -cfg multicrab_run_MuTau_16Nov2011_DATA.cfg

exit

rm ../runMuTauStreamAnalyzerFullAnalysis_Recoil_MC_cfg.py
sed 's/runOnMC     = True/runOnMC     = True/g' ../runMuTauStreamAnalyzerFullAnalysis_Recoil_cfg.py >> ../runMuTauStreamAnalyzerFullAnalysis_Recoil_MC_cfg.py
echo 'multicrab -create -submit -cfg multicrab_run_MuTau_16Nov2011_MC.cfg'
multicrab -create -submit -cfg multicrab_run_MuTau_16Nov2011_MC.cfg


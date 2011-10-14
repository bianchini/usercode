#!/bin/sh

#back up
cp ../runElecTauStreamAnalyzerFullAnalysis_A_cfg.py ../runElecTauStreamAnalyzerFullAnalysis_A_cfg.py.mine

sed 's/runOnMC     = True/runOnMC     = False/g' ../runElecTauStreamAnalyzerFullAnalysis_A_cfg.py
echo 'multicrab -create -submit -cfg multicrab_run_ElecTau_13Oct2011_DATA.cfg'
#multicrab -create -submit -cfg multicrab_run_ElecTau_13Oct2011_DATA.cfg
sed 's/runOnMC     = False/runOnMC     = True/g' ../runElecTauStreamAnalyzerFullAnalysis_A_cfg.py
echo 'multicrab -create -submit -cfg multicrab_run_ElecTau_13Oct2011_MC.cfg'
#multicrab -create -submit -cfg multicrab_run_ElecTau_13Oct2011_MC.cfg

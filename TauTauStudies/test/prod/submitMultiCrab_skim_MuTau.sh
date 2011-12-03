#!/bin/sh

#back up
cp ../patTuple_PAT_SkimMuTauStream_cfg.py ../patTuple_PAT_SkimMuTauStream_cfg.py.mine

#rm ../patTuple_PAT_SkimMuTauStream_DATA_cfg.py
#sed 's/runOnMC     = True/runOnMC     = False/g' ../patTuple_PAT_SkimMuTauStream_cfg.py >> ../patTuple_PAT_SkimMuTauStream_DATA_cfg.py 
#echo 'multicrab -create -submit -cfg multicrab_skim_MuTau_16Nov2011_DATA.cfg'
#multicrab -create -submit -cfg multicrab_skim_MuTau_16Nov2011_DATA.cfg

rm ../patTuple_PAT_SkimMuTauStream_MC_cfg.py
sed 's/runOnMC     = True/runOnMC     = True/g' ../patTuple_PAT_SkimMuTauStream_cfg.py >> ../patTuple_PAT_SkimMuTauStream_MC_cfg.py
echo 'multicrab -create -submit -cfg multicrab_skim_MuTau_16Nov2011_MC.cfg'
multicrab -create -submit -cfg multicrab_skim_MuTau_16Nov2011_MC.cfg

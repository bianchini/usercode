#! /bin/sh

# SVfit
cvs co -r b4_2_x_2012Apr03 AnalysisDataFormats/TauAnalysis
cvs co -r b4_2_x_2012Apr03 TauAnalysis/CandidateTools
cvs up -r joseNov14      TauAnalysis/SVFitStandAlone 

# PAT
addpkg PhysicsTools/PatAlgos      V08-06-55
addpkg DataFormats/PatCandidates  V06-04-19-04
addpkg CommonTools/ParticleFlow   B4_2_X_V00-03-04

# Tau
cvs co -r V01-02-01 RecoTauTag/TauTagTools
cvs co -r V01-02-16 RecoTauTag/RecoTau
cvs co -r V01-02-12 RecoTauTag/Configuration 
cvs up -r 1.47 PhysicsTools/PatAlgos/python/tools/tauTools.py

# MET
cvs co -r CMSSW_4_2_8_patch7 RecoMET/METAlgorithms                            
cvs up -r 1.2 RecoMET/METAlgorithms/interface/SigInputObj.h

# Elec and muons
cvs co UserCode/MitPhysics/data/ElectronMVAWeights
cvs co -r For2011-October-21st-reload HiggsAnalysis/HiggsToWW2Leptons
cvs co -r V00-00-05 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools
cvs co -r V00-00-04 -d Muon/MuonAnalysisTools UserCode/sixie/Muon/MuonAnalysisTools
cvs co -r V00-00-00 UserCode/sixie/EGamma/EGammaAnalysisTools/data/

# Jet
cvs co -r V00-00-09 -d CMGTools/External UserCode/CMG/CMGTools/External
cvs co -r V00-00    -d  pharris/MVAMet UserCode/pharris/MVAMet
cvs co -r V00-04-01 CondFormats/EgammaObjects

# limits
addpkg HiggsAnalysis/CombinedLimit V01-13-02
#! /bin/sh

# SVfit
#cvs co -r b4_2_x_2012Apr03 AnalysisDataFormats/TauAnalysis
#cvs co -r b4_2_x_2012Apr03 TauAnalysis/CandidateTools
cvs up -r b4_2_x_2012May04 TauAnalysis/CandidateTools
cvs up -r b4_2_x_2012May04 AnalysisDataFormats/TauAnalysis
cd TauAnalysis/CandidateTools/bin/
rm test*
rm train*
cd -
cvs co -r joseNov14  UserCode/TauAnalysis/SVFitStandAlone

# PAT
addpkg PhysicsTools/PatAlgos      V08-06-55
addpkg DataFormats/PatCandidates  V06-04-19-04
addpkg CommonTools/ParticleFlow   B4_2_X_V00-03-04
cvs co -r V03-09-18 PhysicsTools/PatUtils 
cvs co -r V01-05-07      PhysicsTools/IsolationAlgos
cvs co -r V08-06-55      PhysicsTools/PatAlgos


# Tau
cvs co -r V01-04-16 RecoTauTag/RecoTau
cvs co -r V01-04-01 RecoTauTag/Configuration 
cvs up -r 1.47 PhysicsTools/PatAlgos/python/tools/tauTools.py
addpkg DataFormats/TauReco CMSSW_5_2_4
addpkg RecoTauTag/TauTagTools CMSSW_5_2_4

# Elec and muons
cvs co UserCode/MitPhysics/data/ElectronMVAWeights
cvs co -r For2011-October-21st-reload HiggsAnalysis/HiggsToWW2Leptons
cvs co -r V00-00-05 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools
cvs co -r V00-00-09 -d Muon/MuonAnalysisTools UserCode/sixie/Muon/MuonAnalysisTools
cvs co -r V00-00-00 UserCode/sixie/EGamma/EGammaAnalysisTools/data/

# Jet
cvs co -r V00-00-09 -d CMGTools/External UserCode/CMG/CMGTools/External
cvs co -r V00-01root -d  pharris/MVAMet UserCode/pharris/MVAMet
#cvs co -r CMSSW_4_2_8_patch7 RecoMET/METAlgorithms 
#cvs co -r CMSSW_4_2_8_patch7 RecoMET/METProducers
cvs co -r V00-04-01 CondFormats/EgammaObjects 
cvs co -r CMSSW_4_2_8_patch7 PhysicsTools/SelectorUtils
cvs up -r 1.22 PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h

cvs co -r CMSSW_4_2_8_patch7 RecoMET/METAlgorithms
cvs up -r b5_2_X_cvMEtCorr_2012May04 RecoMET/METAlgorithms/interface/PFMETAlgorithmMVA.h
cvs up -r b5_2_X_cvMEtCorr_2012May04 RecoMET/METAlgorithms/interface/mvaMEtUtilities.h
cvs up -r b5_2_X_cvMEtCorr_2012May04 RecoMET/METAlgorithms/src/PFMETAlgorithmMVA.cc
cvs up -r b5_2_X_cvMEtCorr_2012May04 RecoMET/METAlgorithms/src/mvaMEtUtilities.cc
cvs up -r b5_2_X_cvMEtCorr_2012May04 RecoMET/METAlgorithms/BuildFile.xml
cvs co -r CMSSW_4_2_8_patch7 RecoMET/METProducers
cvs up -r b5_2_X_cvMEtCorr_2012May04 RecoMET/METProducers/interface/PFMETProducerMVA.h
cvs up -r b5_2_X_cvMEtCorr_2012May04 RecoMET/METProducers/src/PFMETProducerMVA.cc
cvs up -r b5_2_X_cvMEtCorr_2012May04 RecoMET/METProducers/src/SealModule.cc
cvs up -r b5_2_X_cvMEtCorr_2012May04 RecoMET/METProducers/python/mvaPFMET_cff.py
cvs up -r b5_2_X_cvMEtCorr_2012May04 RecoMET/METProducers/BuildFile.xml
cd RecoMET/METProducers/src/
cp /afs/cern.ch/user/b/bianchi/public/SealModule.cc ./
cd -

# limits
addpkg HiggsAnalysis/CombinedLimit V01-13-02

# lumi
cvs co -r V03-05-05      RecoLuminosity/LumiDB

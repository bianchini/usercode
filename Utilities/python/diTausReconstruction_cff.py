import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.svFitAlgorithm_cfi import *

svFitLikelihoodElecTauPairKinematicsPhaseSpace = copy.deepcopy(svFitLikelihoodDiTauKinematicsPhaseSpace)
svFitLikelihoodElecTauPairKinematicsPhaseSpace.pluginType = "SVfitLikelihoodElecTauPairKinematics"
svFitLikelihoodElecTauPairKinematicsPhaseSpace.leg1.pluginType = "SVfitElectronLikelihoodPhaseSpace"
svFitLikelihoodElecTauPairKinematicsPhaseSpace.leg2.pluginType = "SVfitTauLikelihoodPhaseSpace"

svFitLikelihoodMuTauPairKinematicsPhaseSpace = copy.deepcopy(svFitLikelihoodDiTauKinematicsPhaseSpace)
svFitLikelihoodMuTauPairKinematicsPhaseSpace.pluginType = "SVfitLikelihoodMuTauPairKinematics"
svFitLikelihoodMuTauPairKinematicsPhaseSpace.leg1.pluginType = "SVfitMuonLikelihoodPhaseSpace"
svFitLikelihoodMuTauPairKinematicsPhaseSpace.leg2.pluginType = "SVfitTauLikelihoodPhaseSpace"

#svFitLikelihoodDiCandidatePairKinematicsPolarized = copy.deepcopy(svFitLikelihoodDiTauKinematicsPolarized)
#svFitLikelihoodDiCandidatePairKinematicsPolarized.pluginType = "SVfitLikelihoodDiCandidatePairKinematics"
#svFitLikelihoodDiCandidatePairKinematicsPolarized.leg2.pluginType = "SVfitCandidateLikelihoodPolarization"
#svFitLikelihoodDiCandidatePairKinematicsPolarized.leg1 = svFitLikelihoodDiCandidatePairKinematicsPolarized.leg2

svFitLikelihoodElecTauPairMEt = copy.deepcopy(svFitLikelihoodDiTauMEt)
svFitLikelihoodElecTauPairMEt.pluginType = cms.string("SVfitLikelihoodElecTauPairMEt")

svFitLikelihoodMuTauPairMEt = copy.deepcopy(svFitLikelihoodDiTauMEt)
svFitLikelihoodMuTauPairMEt.pluginType = cms.string("SVfitLikelihoodMuTauPairMEt")

#svFitLikelihoodDiCandidatePairMEt = copy.deepcopy(svFitLikelihoodMEt)
#svFitLikelihoodDiCandidatePairMEt.pluginType = cms.string("SVfitLikelihoodMEtDiCandidate")

svFitLikelihoodElecTauPairPtBalance = copy.deepcopy(svFitLikelihoodDiTauPtBalance)
svFitLikelihoodElecTauPairPtBalance.pluginType = cms.string("SVfitLikelihoodElecTauPairPtBalance")

svFitLikelihoodMuTauPairPtBalance = copy.deepcopy(svFitLikelihoodDiTauPtBalance)
svFitLikelihoodMuTauPairPtBalance.pluginType = cms.string("SVfitLikelihoodMuTauPairPtBalance")

#svFitLikelihoodDiCandidatePairPtBalance = copy.deepcopy(svFitLikelihoodDiTauPtBalance)
#svFitLikelihoodDiCandidatePairPtBalance.pluginType = cms.string("SVfitLikelihoodDiCandidatePairPtBalance")

#svFitLikelihoodDiCandidatePairZprod = copy.deepcopy(svFitLikelihoodDiTauProdZ0)
#svFitLikelihoodDiCandidatePairZprod.pluginType = cms.string("SVfitLikelihoodDiCandidatePairProd")
#svFitLikelihoodDiCandidatePairZprod.process = cms.string("Z0")

#--------------------------------------------------------------------------------
# produce combinations of tau-jet + tau-jet pairs
#--------------------------------------------------------------------------------

elecTauPairs = cms.EDProducer("PATElecTauPairProducer",
                              useLeadingTausOnly = cms.bool(False),
                              srcLeg1 = cms.InputTag('patElectrons'),
                              srcLeg2 = cms.InputTag('patTaus'),
                              dRmin12 = cms.double(0.3),
                              srcMET = cms.InputTag('patMETs'),
                              srcPrimaryVertex = cms.InputTag("offlinePrimaryVerticesWithBS"),
                              srcBeamSpot = cms.InputTag("offlineBeamSpot"),
                              srcGenParticles = cms.InputTag('genParticles'),                  
                              recoMode = cms.string(""),
                              doSVreco = cms.bool(True),
                              
                              svFit = cms.PSet(
    psKine = cms.PSet(
    likelihoodFunctions = cms.VPSet(
    #svFitLikelihoodDiTauPairKinematicsPhaseSpace
    svFitLikelihoodElecTauPairKinematicsPhaseSpace
    ),
    estUncertainties = cms.PSet(
    numSamplings = cms.int32(-1)
    )
    ),
    psKine_MEt = cms.PSet(
    likelihoodFunctions = cms.VPSet(
    #svFitLikelihoodDiCandidatePairKinematicsPhaseSpace,
    #svFitLikelihoodDiCandidatePairMEt
    svFitLikelihoodElecTauPairKinematicsPhaseSpace,
    svFitLikelihoodElecTauPairMEt
    ),
    estUncertainties = cms.PSet(
    numSamplings = cms.int32(-1)
    )
    ),
    psKine_MEt_ptBalance = cms.PSet(
    likelihoodFunctions = cms.VPSet(
    #svFitLikelihoodDiCandidatePairKinematicsPhaseSpace,
    #svFitLikelihoodDiCandidatePairMEt,
    #svFitLikelihoodDiCandidatePairPtBalance
    svFitLikelihoodElecTauPairKinematicsPhaseSpace,
    svFitLikelihoodElecTauPairMEt,
    svFitLikelihoodElecTauPairPtBalance
    ),
    estUncertainties = cms.PSet(
    numSamplings = cms.int32(-1)
    )
    )
    ),
                               
                              
                              scaleFuncImprovedCollinearApprox = cms.string('1'),                           
                              verbosity = cms.untracked.int32(0)
                              )



muTauPairs = cms.EDProducer("PATMuTauPairProducer",
                            useLeadingTausOnly = cms.bool(False),
                            srcLeg1 = cms.InputTag('patMuons'),
                            srcLeg2 = cms.InputTag('patTaus'),
                            dRmin12 = cms.double(0.3),
                            srcMET = cms.InputTag('patMETs'),
                            srcPrimaryVertex = cms.InputTag("offlinePrimaryVerticesWithBS"),
                            srcBeamSpot = cms.InputTag("offlineBeamSpot"),
                            srcGenParticles = cms.InputTag('genParticles'),                  
                            recoMode = cms.string(""),
                            doSVreco = cms.bool(True),
                            
                            svFit = cms.PSet(
    psKine = cms.PSet(
    likelihoodFunctions = cms.VPSet(
    #svFitLikelihoodDiTauPairKinematicsPhaseSpace
    svFitLikelihoodMuTauPairKinematicsPhaseSpace
    ),
    estUncertainties = cms.PSet(
    numSamplings = cms.int32(-1)
    )
    ),
    psKine_MEt = cms.PSet(
    likelihoodFunctions = cms.VPSet(
    #svFitLikelihoodDiCandidatePairKinematicsPhaseSpace,
    #svFitLikelihoodDiCandidatePairMEt
    svFitLikelihoodMuTauPairKinematicsPhaseSpace,
    svFitLikelihoodMuTauPairMEt
    ),
    estUncertainties = cms.PSet(
    numSamplings = cms.int32(-1)
    )
    ),
    psKine_MEt_ptBalance = cms.PSet(
    likelihoodFunctions = cms.VPSet(
    #svFitLikelihoodDiCandidatePairKinematicsPhaseSpace,
    #svFitLikelihoodDiCandidatePairMEt,
    #svFitLikelihoodDiCandidatePairPtBalance
    svFitLikelihoodMuTauPairKinematicsPhaseSpace,
    svFitLikelihoodMuTauPairMEt,
    svFitLikelihoodMuTauPairPtBalance
    ),
    estUncertainties = cms.PSet(
    numSamplings = cms.int32(-1)
    )
    )
    ),
                               
                               
                               scaleFuncImprovedCollinearApprox = cms.string('1'),                           
                               verbosity = cms.untracked.int32(0)
                              )

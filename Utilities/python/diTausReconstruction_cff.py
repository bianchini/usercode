import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.svFitAlgorithm_cfi import *

svFitLikelihoodElecTauPairKinematicsPhaseSpace = copy.deepcopy(svFitLikelihoodDiTauKinematicsPhaseSpace)
svFitLikelihoodElecTauPairKinematicsPhaseSpace.pluginType = "SVfitLikelihoodElecTauPairKinematics"
svFitLikelihoodElecTauPairKinematicsPhaseSpace.leg1.pluginType = "SVfitElectronLikelihoodPhaseSpace"
svFitLikelihoodElecTauPairKinematicsPhaseSpace.leg2.pluginType = "SVfitTauLikelihoodPhaseSpace"

svFitLikelihoodDiCandidatePairKinematicsPolarized = copy.deepcopy(svFitLikelihoodDiTauKinematicsPolarized)
svFitLikelihoodDiCandidatePairKinematicsPolarized.pluginType = "SVfitLikelihoodDiCandidatePairKinematics"
svFitLikelihoodDiCandidatePairKinematicsPolarized.leg2.pluginType = "SVfitCandidateLikelihoodPolarization"
svFitLikelihoodDiCandidatePairKinematicsPolarized.leg1 = svFitLikelihoodDiCandidatePairKinematicsPolarized.leg2

svFitLikelihoodDiCandidatePairMEt = copy.deepcopy(svFitLikelihoodMEt)
svFitLikelihoodDiCandidatePairMEt.pluginType = cms.string("SVfitLikelihoodMEtDiCandidate")

svFitLikelihoodDiCandidatePairPtBalance = copy.deepcopy(svFitLikelihoodDiTauPtBalance)
svFitLikelihoodDiCandidatePairPtBalance.pluginType = cms.string("SVfitLikelihoodDiCandidatePairPtBalance")

svFitLikelihoodDiCandidatePairZprod = copy.deepcopy(svFitLikelihoodDiTauProdZ0)
svFitLikelihoodDiCandidatePairZprod.pluginType = cms.string("SVfitLikelihoodDiCandidatePairProd")
svFitLikelihoodDiCandidatePairZprod.process = cms.string("Z0")

#--------------------------------------------------------------------------------
# produce combinations of tau-jet + tau-jet pairs
#--------------------------------------------------------------------------------

allDiTauPairs = cms.EDProducer("PATElecTauPairProducer",
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
    #psKine_MEt = cms.PSet(
    #likelihoodFunctions = cms.VPSet(
    #svFitLikelihoodDiCandidatePairKinematicsPhaseSpace,
    #svFitLikelihoodDiCandidatePairMEt
    #),
    #estUncertainties = cms.PSet(
    #numSamplings = cms.int32(-1)
    #)
    #),
    #psKine_MEt_ptBalance = cms.PSet(
    #likelihoodFunctions = cms.VPSet(
    #svFitLikelihoodDiCandidatePairKinematicsPhaseSpace,
    #svFitLikelihoodDiCandidatePairMEt,
    #svFitLikelihoodDiCandidatePairPtBalance
    #),
    #estUncertainties = cms.PSet(
    #numSamplings = cms.int32(-1)
    #)
    #)
    ),
                               
                               
                               scaleFuncImprovedCollinearApprox = cms.string('1'),                           
                               verbosity = cms.untracked.int32(0)
                               )

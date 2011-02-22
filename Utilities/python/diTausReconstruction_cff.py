import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.svFitAlgorithm_cfi import *

svFitLikelihoodDiTauPairKinematicsPhaseSpace = copy.deepcopy(svFitLikelihoodDiTauKinematicsPhaseSpace)
svFitLikelihoodDiTauPairKinematicsPhaseSpace.pluginType = "SVfitLikelihoodDiTauPairKinematics"
svFitLikelihoodDiTauPairKinematicsPhaseSpace.leg1.pluginType = "SVfitTauLikelihoodPhaseSpace"
svFitLikelihoodDiTauPairKinematicsPhaseSpace.leg2.pluginType = "SVfitTauLikelihoodPhaseSpace"

svFitLikelihoodDiTauPairKinematicsPolarized = copy.deepcopy(svFitLikelihoodDiTauKinematicsPolarized)
svFitLikelihoodDiTauPairKinematicsPolarized.pluginType = "SVfitLikelihoodDiTauPairKinematics"
svFitLikelihoodDiTauPairKinematicsPolarized.leg2.pluginType = "SVfitTauLikelihoodPolarization"
svFitLikelihoodDiTauPairKinematicsPolarized.leg1 = svFitLikelihoodDiTauPairKinematicsPolarized.leg2

svFitLikelihoodDiTauPairMEt = copy.deepcopy(svFitLikelihoodMEt)
svFitLikelihoodDiTauPairMEt.pluginType = cms.string("SVfitLikelihoodMEtDiTau")

svFitLikelihoodDiTauPairPtBalance = copy.deepcopy(svFitLikelihoodDiTauPtBalance)
svFitLikelihoodDiTauPairPtBalance.pluginType = cms.string("SVfitLikelihoodDiTauPairPtBalance")

svFitLikelihoodDiTauPairZprod = copy.deepcopy(svFitLikelihoodDiTauProdZ0)
svFitLikelihoodDiTauPairZprod.pluginType = cms.string("SVfitLikelihoodDiTauPairProd")
svFitLikelihoodDiTauPairZprod.process = cms.string("Z0")

#--------------------------------------------------------------------------------
# produce combinations of tau-jet + tau-jet pairs
#--------------------------------------------------------------------------------

allDiTauPairs = cms.EDProducer("DiCandidatePairProducer",
                               useLeadingTausOnly = cms.bool(False),
                               srcLeg1 = cms.InputTag('patTaus'),
                               srcLeg2 = cms.InputTag('patElectrons'),
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
    svFitLikelihoodDiTauPairKinematicsPhaseSpace         
    ),
    estUncertainties = cms.PSet(
    numSamplings = cms.int32(-1)
    )
    ),
    psKine_MEt = cms.PSet(
    likelihoodFunctions = cms.VPSet(
    svFitLikelihoodDiTauPairKinematicsPhaseSpace,
    svFitLikelihoodDiTauPairMEt
    ),
    estUncertainties = cms.PSet(
    numSamplings = cms.int32(-1)
    )
    ),
    psKine_MEt_ptBalance = cms.PSet(
    likelihoodFunctions = cms.VPSet(
    svFitLikelihoodDiTauPairKinematicsPhaseSpace,
    svFitLikelihoodDiTauPairMEt,
    svFitLikelihoodDiTauPairPtBalance
    ),
    estUncertainties = cms.PSet(
    numSamplings = cms.int32(-1)
    )
    )
    ),
                               
                               
                               scaleFuncImprovedCollinearApprox = cms.string('1'),                           
                               verbosity = cms.untracked.int32(0)
                               )

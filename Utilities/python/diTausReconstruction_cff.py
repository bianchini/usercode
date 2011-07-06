import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.tools.objProdConfigurator import *
from TauAnalysis.CandidateTools.resolutions_cfi import *
from TauAnalysis.CandidateTools.nSVfitAlgorithmDiTau_cfi import *
from TauAnalysis.CandidateTools.nSVfitAlgorithmTauDecayKineMC_cfi import *
from RecoMET.METProducers.METSigParams_cfi import *

#--------------------------------------------------------------------------------
# produce combinations of electron + tau-jet pairs
#--------------------------------------------------------------------------------

allElecTauPairs = cms.EDProducer(
    "PATElecTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedPatElectrons'),
    srcLeg2 = cms.InputTag('selectedPatTaus'),
    dRmin12 = cms.double(0.3),
    srcMET = cms.InputTag('patMETs'),
    srcPrimaryVertex = cms.InputTag("offlinePrimaryVertices"),
    srcBeamSpot = cms.InputTag("offlineBeamSpot"),
    srcGenParticles = cms.InputTag('genParticles'),
    recoMode = cms.string(""),
    doSVreco = cms.bool(True),
    nSVfit = cms.PSet(),
    scaleFuncImprovedCollinearApprox = cms.string('1'),
    verbosity = cms.untracked.int32(0),
    )


allElecTauPairs.nSVfit.psKine_MEt_logM_fit = cms.PSet()
allElecTauPairs.nSVfit.psKine_MEt_logM_fit.config = copy.deepcopy(nSVfitConfig_template)
allElecTauPairs.nSVfit.psKine_MEt_logM_fit.config.event.resonances.A.daughters.leg1 = cms.PSet(
        src = allElecTauPairs.srcLeg1,
            likelihoodFunctions = cms.VPSet(nSVfitElectronLikelihoodPhaseSpace),
            builder = nSVfitTauToElecBuilder
        )
allElecTauPairs.nSVfit.psKine_MEt_logM_fit.config.event.resonances.A.daughters.leg2 = cms.PSet(
        src = allElecTauPairs.srcLeg2,
            likelihoodFunctions = cms.VPSet(nSVfitTauLikelihoodPhaseSpace),
            builder = nSVfitTauToHadBuilder
        )
allElecTauPairs.nSVfit.psKine_MEt_logM_fit.algorithm = cms.PSet(
        pluginName = cms.string("nSVfitAlgorithmByLikelihoodMaximization"),
            pluginType = cms.string("NSVfitAlgorithmByLikelihoodMaximization"),
            minimizer  = cms.vstring("Minuit2", "Migrad"),
            maxObjFunctionCalls = cms.uint32(5000),
            verbosity = cms.int32(0)
        )


allElecTauPairs.nSVfit.psKine_MEt_logM_int = cms.PSet()
allElecTauPairs.nSVfit.psKine_MEt_logM_int.config = allElecTauPairs.nSVfit.psKine_MEt_logM_fit.config
allElecTauPairs.nSVfit.psKine_MEt_logM_int.algorithm = nSVfitProducerByIntegration.algorithm
'''
cms.PSet(
        pluginName = cms.string("nSVfitAlgorithmByIntegration"),
            pluginType = cms.string("NSVfitAlgorithmByIntegration"),
            parameters = cms.PSet(
            mass_A = cms.PSet(
                min = cms.double(20.),
                            max = cms.double(750.),
                            stepSize = cms.double(5.),
                            replace = cms.string("leg1.x"),
                            by = cms.string("(A.p4.mass/mass_A)*(A.p4.mass/mass_A)/leg2.x")
                        )
                ),
            vegasOptions = cms.PSet(
            numCalls = cms.uint32(10000)
                )
        )
'''
#--------------------------------------------------------------------------------

elecTauPairProdConfigurator = objProdConfigurator(
            allElecTauPairs,
            pyModuleName = __name__
            )

produceElecTauPairs = elecTauPairProdConfigurator.configure(pyNameSpace = locals())

##################################################################################

allMuTauPairs = cms.EDProducer(
    "PATMuTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedPatMuons'),
    srcLeg2 = cms.InputTag('selectedPatTaus'),
    dRmin12 = cms.double(0.3),
    srcMET = cms.InputTag('patMETs'),
    srcPrimaryVertex = cms.InputTag("offlinePrimaryVertices"),
    srcBeamSpot = cms.InputTag("offlineBeamSpot"),
    srcGenParticles = cms.InputTag('genParticles'),
    recoMode = cms.string(""),
    doSVreco = cms.bool(True),
    nSVfit = cms.PSet(),
    scaleFuncImprovedCollinearApprox = cms.string('1'),
    doPFMEtSign = cms.bool(True),
    pfMEtSign = cms.PSet(
    srcPFJets = cms.InputTag('ak5PFJets'),
    srcPFCandidates = cms.InputTag('particleFlow'),
    resolution = METSignificance_params,
    dRoverlapPFJet = cms.double(0.3),
    dRoverlapPFCandidate = cms.double(0.1)
    ),
    verbosity = cms.untracked.int32(0)
    )

#--------------------------------------------------------------------------------
# configure (new) SVfit algorithm
# (using combination of PS + MET likelihoods + logM regularization term
#  to reconstruct mass of tau lepton pair, as described in CMS AN-11-165)
allMuTauPairs.nSVfit.psKine_MEt_logM_fit = cms.PSet()
allMuTauPairs.nSVfit.psKine_MEt_logM_fit.config = copy.deepcopy(nSVfitConfig_template)
allMuTauPairs.nSVfit.psKine_MEt_logM_fit.config.event.resonances.A.daughters.leg1 = cms.PSet(
        src = allMuTauPairs.srcLeg1,
            likelihoodFunctions = cms.VPSet(nSVfitMuonLikelihoodPhaseSpace),
            builder = nSVfitTauToMuBuilder
        )
allMuTauPairs.nSVfit.psKine_MEt_logM_fit.config.event.resonances.A.daughters.leg2 = cms.PSet(
        src = allMuTauPairs.srcLeg2,
            likelihoodFunctions = cms.VPSet(nSVfitTauLikelihoodPhaseSpace),
            builder = nSVfitTauToHadBuilder
        )
allMuTauPairs.nSVfit.psKine_MEt_logM_fit.algorithm = cms.PSet(
        pluginName = cms.string("nSVfitAlgorithmByLikelihoodMaximization"),
            pluginType = cms.string("NSVfitAlgorithmByLikelihoodMaximization"),
            minimizer  = cms.vstring("Minuit2", "Migrad"),
            maxObjFunctionCalls = cms.uint32(5000),
            verbosity = cms.int32(0)
        )


allMuTauPairs.nSVfit.psKine_MEt_logM_int = cms.PSet()
allMuTauPairs.nSVfit.psKine_MEt_logM_int.config = allMuTauPairs.nSVfit.psKine_MEt_logM_fit.config
allMuTauPairs.nSVfit.psKine_MEt_logM_int.algorithm = nSVfitProducerByIntegration.algorithm
'''
cms.PSet(
    pluginName = cms.string("nSVfitAlgorithmByIntegration"),
    pluginType = cms.string("NSVfitAlgorithmByIntegration"),
    parameters = cms.PSet(
    mass_A = cms.PSet(
    min = cms.double(20.),
    max = cms.double(750.),
    #stepSize = cms.double(5.),
    stepSizeFactor = cms.double(1.025), # nextM = max(stepSizeFactor*currentM, minStepSize)
    minStepSize = cms.double(2.5),      
    replace = cms.string("leg1.x"),
    by = cms.string("(A.p4.mass/mass_A)*(A.p4.mass/mass_A)/leg2.x")
    )
    ),
    vegasOptions = cms.PSet(
    #numCalls = cms.uint32(10000)
    numCallsGridOpt = cms.uint32(1000),
    numCallsIntEval = cms.uint32(10000),
    maxChi2         = cms.double(2.),
    maxIntEvalIter  = cms.uint32(5),                                          
    precision       = cms.double(0.00001)
    )
    )
'''

muTauPairProdConfigurator = objProdConfigurator(
    allMuTauPairs,
    pyModuleName = __name__
    )

produceMuTauPairs = muTauPairProdConfigurator.configure(pyNameSpace = locals())


'''
import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.nSVfitAlgorithmDiTau_cfi import *

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
'''

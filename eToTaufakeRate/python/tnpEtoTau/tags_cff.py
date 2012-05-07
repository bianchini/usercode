import FWCore.ParameterSet.Config as cms

#hltMatching = "triggerObjectMatchesByPath('HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v').size() != 0"

elePreID = "(userFloat('antiConv') > 0.5 && userFloat('nHits') <= 0 && userFloat('dxyWrtPV') < 0.02 && userFloat('dzWrtPV') < 0.1) && ((abs(superClusterPosition.Eta) < 1.479 && userFloat('sihih') < 0.01 && userFloat('dEta') < 0.007 && userFloat('dPhi') < 0.15 && userFloat('HoE') < 0.12) || (abs(superClusterPosition.Eta) > 1.479 && userFloat('sihih') < 0.03 && userFloat('dEta') < 0.009 && userFloat('dPhi') < 0.10 && userFloat('HoE') < 0.10))"

electronMVA = "((pt < 20 && abs(superClusterPosition.Eta) >= 0.0 && abs(superClusterPosition.Eta) < 1.0 && userFloat('mva') > 0.133) || (pt < 20 && abs(superClusterPosition.Eta) >= 1.0 && abs(superClusterPosition.Eta) < 1.5 && userFloat('mva') > 0.465) || (pt < 20 && abs(superClusterPosition.Eta) >= 1.5 && abs(superClusterPosition.Eta) < 2.5 && userFloat('mva') > 0.518) || (pt > 20 && abs(superClusterPosition.Eta) >= 0.0 && abs(superClusterPosition.Eta) < 1.0 && userFloat('mva') > 0.942) || (pt > 20 && abs(superClusterPosition.Eta) >= 1.0 && abs(superClusterPosition.Eta) < 1.5 && userFloat('mva') > 0.947) || (pt > 20 && abs(superClusterPosition.Eta) >= 1.5 && abs(superClusterPosition.Eta) < 2.5 && userFloat('mva') > 0.878))"

simpleCutsWP95 = "(userFloat('nHits')<=1"+ \
                 " && (" + \
                 " (isEB && userFloat('sihih')<0.010 && userFloat('dPhi')<0.80 && "+ \
                 "          userFloat('dEta') <0.007 && userFloat('HoE') <0.15)"   + \
                 " || "  + \
                 " (isEE && userFloat('sihih')<0.030 && userFloat('dPhi')<0.70 && "+ \
                 "          userFloat('dEta') <0.010 && userFloat('HoE') <0.07)"   + \
                 "     )"+ \
                 ")"

simpleCutsWP80 = "(userFloat('nHits')==0 && userInt('antiConv')>0.5 "+ \
                 " && ("   + \
                 " (isEB && userFloat('sihih')<0.010 && userFloat('dPhi')<0.06 && "  + \
                 "          userFloat('dEta')< 0.004 && userFloat('HoE') <0.04)"     + \
                 " ||"+ \
                 " (isEE && userFloat('sihih')<0.030 && userFloat('dPhi')<0.030 && " + \
                 "          userFloat('dEta') <0.007 && userFloat('HoE') <0.025) "   + \
                 "    )"   + \
                 ")"

simpleCutsWP70 = "(userFloat('nHits')==0 && userInt('antiConv')>0.5 "+ \
                 " && ("   + \
                 " (isEB && userFloat('sihih')<0.010 && userFloat('dPhi')<0.03 && "  + \
                 "          userFloat('dEta')< 0.004 && userFloat('HoE') <0.025)"     + \
                 " ||"+ \
                 " (isEE && userFloat('sihih')<0.030 && userFloat('dPhi')<0.020 && " + \
                 "          userFloat('dEta') <0.005 && userFloat('HoE') <0.025) "   + \
                 "    )"   + \
                 ")"

#WW ID (HoE nell'endcap sarebbe 0.025 ma in WW e' 0.1)
eleWWID = "(userFloat('nHits')==0 && userInt('antiConv')>0.5 "+ \
                 " && ("   + \
                 " (pt>=20 && ("    + \
                 "               (isEB && userFloat('sihih')<0.010 && userFloat('dPhi')<0.06 && "  + \
                 "                        userFloat('dEta')< 0.004 && userFloat('HoE') <0.04)"     + \
                 "               ||"+ \
                 "               (isEE && userFloat('sihih')<0.030 && userFloat('dPhi')<0.030 && " + \
                 "                        userFloat('dEta') <0.007 && userFloat('HoE') <0.1) )) "+ \
                 "     || "+ \
                 " (pt<20 && (fbrem>0.15 || (abs(superClusterPosition.Eta)<1. && eSuperClusterOverP>0.95) ) && ( "+ \
                 "               (isEB && userFloat('sihih')<0.010 && userFloat('dPhi')<0.030 && " + \
                 "                        userFloat('dEta') <0.004 && userFloat('HoE') <0.025) "   + \
                 "               ||"+ \
                 "               (isEE && userFloat('sihih')<0.030 && userFloat('dPhi')<0.020 &&"  + \
                 "                        userFloat('dEta') <0.005 && userFloat('HoE') <0.1) ))" + \
                 "    )"   + \
                 ")"

tagAnyEle = cms.EDFilter("PATElectronRefSelector",
                         src = cms.InputTag("electronVariables"),
                         #cut = cms.string("pt > 25.0 && abs(eta) < 2.5 && !(1.4442 < abs(eta) < 1.566) && electronID('simpleEleId70relIso') > 6.5"),
                         #cut = cms.string("pt > 25.0 && abs(eta) < 2.5 && !(1.4442 < abs(eta) < 1.566) && " + elePreID + electronMVA + "&& userFloat('PFRelIsoDB04v2') < 0.1"),
			 cut = cms.string("pt > 25.0 && abs(eta) < 2.5 && !isEBEEGap &&" + simpleCutsWP80 + "&& userFloat('PFRelIsoDB04v2') < 0.1"),
                         filter = cms.bool(False)
                         )

#eleWP95 = cms.EDFilter("PATElectronRefSelector",
#                         src = cms.InputTag("electronVariables"),
#			  cut = cms.string("pt > 25.0 && abs(eta) < 2.5 && !isEBEEGap &&" + simpleCutsWP95 + "&& userFloat('PFRelIsoDB04v2') < 0.1"),
#                         filter = cms.bool(False)
#                         )

#OneTaggedEle = cms.EDFilter("PATCandViewCountFilter",
#                     	src = cms.InputTag("tagAnyEle"),
#                     	maxNumber = cms.uint32(1),
#                     	minNumber = cms.uint32(1),
#                       filter = cms.bool(True)
#                     	)

vecMC = (0.0020804,0.0109929,0.0228113,0.0429041,0.0566695,0.0721029,0.0729181,0.0711529,0.0713194,0.0605273,0.0531286,0.0434069,0.0446853,0.0439974,0.0422232,0.0381212,0.0286534,0.0336355,0.0270326,0.0256726,0.0236528,0.0217785,0.0180944,0.0171807,0.0145082,0.0112167,0.00852673,0.00395336,0.00458083,0.00435781,0.00354159,0.00167699,0.00103091,0.000206458,0.000618599,0.000206121,0.000412579,8.83568e-06,0,0.000412545,0,0,0,0,0,0,0,0,0,0)
vecDataA = (0,0.000115513,0.00262453,0.0231253,0.122659,0.234608,0.23849,0.171807,0.110434,0.0654693,0.0246198,0.00540239,0.000590312,5.10257e-05,2.9901e-06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
vecDataB = (0,2.0743e-05,1.69768e-05,4.20265e-05,0.000413519,0.00233585,0.0219252,0.0625536,0.0892844,0.112663,0.130604,0.135136,0.130953,0.116884,0.0911061,0.0578308,0.0295598,0.0122804,0.00433753,0.00146337,0.000449148,0.000114735,2.54458e-05,9.34371e-07,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

addUserVariables = cms.EDProducer("UserDefinedVariables",
    objects = cms.InputTag("electronVariables"),
    objects2 = cms.InputTag("patTaus"),
    triggerResults = cms.InputTag("patTriggerEvent"),
    met = cms.InputTag("patMETsPF"),
    isMC = cms.bool(False), #####
    TrueDist2011A = cms.vdouble(vecDataA),
    TrueDist2011B = cms.vdouble(vecDataB),
    MCDist = cms.vdouble(vecMC),
    )

electronVariables = cms.EDProducer('ElectronsUserEmbedder',
    electronTag = cms.InputTag("selectedPatElectronsTriggerMatch"),
    vertexTag = cms.InputTag("offlinePrimaryVerticesWithBS"),
    isMC = cms.bool(True),
    doMVA = cms.bool(True),
    inputFileName0 = cms.FileInPath("UserCode/MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_NoIPInfo_BDTG.weights.xml"),
    inputFileName1 = cms.FileInPath("UserCode/MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_NoIPInfo_BDTG.weights.xml"),
    inputFileName2 = cms.FileInPath("UserCode/MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_NoIPInfo_BDTG.weights.xml"),
    inputFileName3 = cms.FileInPath("UserCode/MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_NoIPInfo_BDTG.weights.xml"),
    inputFileName4 = cms.FileInPath("UserCode/MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_NoIPInfo_BDTG.weights.xml"),
    inputFileName5 = cms.FileInPath("UserCode/MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_NoIPInfo_BDTG.weights.xml"),
)


from PhysicsTools.PatAlgos.patTemplate_cfg import *



process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

postfix     = "PF"
runOnMC     = False

#from Configuration.PyReleaseValidation.autoCond import autoCond
#process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

if runOnMC:
    process.GlobalTag.globaltag = cms.string('START52_V7::All')

else:
    process.GlobalTag.globaltag = cms.string('GR_R_52_V7::All')
    
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source.fileNames = cms.untracked.vstring(
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Test_52X/439eae0c5ab5b9ee314b9645ffd1dade/DYJets_1_1_MAU.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/TauPlusX/Test_52X/439eae0c5ab5b9ee314b9645ffd1dade/Run2012_9_1_Vc9.root'
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Test_51X/439eae0c5ab5b9ee314b9645ffd1dade/out_10_1_8Te.root'
    )

#process.source.eventsToProcess = cms.untracked.VEventRange(
#    '1:1080'
#    )
################### event content ##################

process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

################### filters log  ####################

process.allEventsFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.primaryVertexFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.atLeastOneElecTauFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.elecPtEtaFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.elecPtEtaIDFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.tauPtEtaFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.tauPtEtaIDFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.tauPtEtaIDAgMuFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.tauPtEtaIDAgMuAgElecFilter = cms.EDFilter(
    "AllEventsFilter"
    )

################### gen listing  ####################

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree1 = cms.EDAnalyzer(
    "ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint  = cms.untracked.int32(1)
    )

################### jet sequence ####################

#process.load('RecoJets.Configuration.RecoPFJets_cff')

#process.kt6PFJets.doRhoFastjet  = True
#process.kt6PFJets.Rho_EtaMax    = cms.double(4.4)
#process.kt6PFJets.Ghost_EtaMax  = cms.double(5.0)
#process.ak5PFJets.doAreaFastjet = True

## re-run kt4PFJets within lepton acceptance to compute rho
#process.load('RecoJets.JetProducers.kt4PFJets_cfi')

#process.kt6PFJetsCentral = process.kt4PFJets.clone(
#    rParam = 0.6,
#    doRhoFastjet = True )
#process.kt6PFJetsCentral.Rho_EtaMax   = cms.double(1.9)
#process.kt6PFJetsCentral.Ghost_EtaMax = cms.double(2.5)

#process.kt6PFJetsNeutral = process.kt4PFJets.clone(
#    rParam = 0.6,
#    doRhoFastjet = True,
#    src="pfAllNeutral" )
#process.kt6PFJetsNeutral.Rho_EtaMax   = cms.double(1.9)
#process.kt6PFJetsNeutral.Ghost_EtaMax = cms.double(2.5)

#process.fjSequence = cms.Sequence(process.kt6PFJets+process.ak5PFJets+process.kt6PFJetsCentral)


################### vertex sequence ####################

process.selectedPrimaryVertices = cms.EDFilter(
    "VertexSelector",
    src = cms.InputTag('offlinePrimaryVertices'),
    cut = cms.string("isValid & ndof >= 4 & z > -24 & z < +24 & position.Rho < 2."),
    filter = cms.bool(False)                                          
)

process.primaryVertexCounter = cms.EDFilter(
    "VertexCountFilter",
    src = cms.InputTag('selectedPrimaryVertices'),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

################### pat specific ####################

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag

from PhysicsTools.PatAlgos.tools.coreTools import *
if not runOnMC:
    removeMCMatching(process,["All"])
    
#removeSpecificPATObjects(process, ['Photons'],
#                         outputInProcess=False)
#removeCleaning(process,
#               outputInProcess=False)

restrictInputToAOD(process, ['All'])

from Bianchi.Utilities.customizePAT  import *
addSelectedPFlowParticle(process)

from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')


from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process,cms.InputTag('ak5PFJets'),
                    doJTA        = True,
                    doBTagging   = True,
                    jetCorrLabel = ('AK5PF', ['L1FastJet','L2Relative', 'L3Absolute']),
                    doType1MET   = True,
                    genJetCollection=cms.InputTag("ak5GenJets"),
                    doJetID      = True,
                    )

JEClevels = cms.vstring(['L2Relative', 'L3Absolute'])
if runOnMC:
    JEClevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
else:
    JEClevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']

process.patJetCorrFactors.levels = JEClevels
process.patJetCorrFactors.rho = cms.InputTag('kt6PFJets','rho')
process.patJetCorrFactors.useRho = True

#process.patJetCorrFactorsL1Offset = process.patJetCorrFactors.clone(
#    levels = cms.vstring('L1Offset',
#                         'L2Relative',
#                         'L3Absolute')
#    )
#
#if runOnMC:
#    process.patJetCorrFactorsL1Offset.levels = ['L1Offset', 'L2Relative', 'L3Absolute']
#else:
#    process.patJetCorrFactorsL1Offset.levels = ['L1Offset', 'L2Relative', 'L3Absolute', 'L2L3Residual']
#    
#process.patJets.jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactors"),cms.InputTag("patJetCorrFactorsL1Offset"))
#process.patDefaultSequence.replace(process.patJetCorrFactors,
#                                   process.patJetCorrFactors+process.patJetCorrFactorsL1Offset)



if runOnMC:
    process.load("RecoJets.Configuration.GenJetParticles_cff")
    process.load("RecoJets.Configuration.RecoGenJets_cff")
    process.genJetsNoNu = cms.Sequence(process.genParticlesForJetsNoNu*
                                       process.ak5GenJetsNoNu)
    process.patDefaultSequence.replace(process.patJetGenJetMatch,
                                       process.genJetsNoNu*
                                       process.patJetGenJetMatch)
    process.patJetGenJetMatch.matched = cms.InputTag("ak5GenJetsNoNu")


#################### tau sequence #######################

#process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
#from PhysicsTools.PatAlgos.tools.tauTools import *
#switchToPFTauHPS(process, 
#                 pfTauLabelOld = 'shrinkingConePFTauProducer',
#                 pfTauLabelNew = 'hpsPFTauProducer'
#                 )

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")


getattr(process,"patTaus").embedIsolationTracks             = cms.bool(True)
getattr(process,"patTaus").embedSignalTracks                = cms.bool(True)
getattr(process,"patTaus").embedGenMatch                    = cms.bool(True)
getattr(process,"patTaus").embedLeadTrack                   = cms.bool(True)
getattr(process,"patTaus").embedLeadPFCand                  = cms.bool(True)
getattr(process,"patTaus").embedLeadPFChargedHadrCand       = cms.bool(True)
getattr(process,"patTaus").embedLeadPFNeutralCand           = cms.bool(True)
getattr(process,"patTaus").embedSignalPFCands               = cms.bool(True)
getattr(process,"patTaus").embedSignalPFChargedHadrCands    = cms.bool(True)
getattr(process,"patTaus").embedSignalPFNeutralHadrCands    = cms.bool(True)
getattr(process,"patTaus").embedSignalPFGammaCands          = cms.bool(True)
getattr(process,"patTaus").embedIsolationPFCands            = cms.bool(True)
getattr(process,"patTaus").embedIsolationPFChargedHadrCands = cms.bool(True)
getattr(process,"patTaus").embedIsolationPFNeutralHadrCands = cms.bool(True)
getattr(process,"patTaus").embedIsolationPFGammaCands       = cms.bool(True)
getattr(process,"patTaus").embedGenJetMatch                 = cms.bool(True)

process.tauMatch.maxDeltaR                = 0.15
process.tauMatch.resolveAmbiguities       = cms.bool(False)
process.tauGenJetMatch.resolveAmbiguities = cms.bool(False)
process.tauGenJetMatch.maxDeltaR          = 0.15
process.tauGenJetMatch.maxDPtRel          = 999

##################################################################

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFMuonIso, setupPFElectronIso
process.muIsoSequence       = setupPFMuonIso(process,'muons')
process.electronIsoSequence = setupPFElectronIso(process,'gsfElectrons')
from CommonTools.ParticleFlow.pfParticleSelection_cff import pfParticleSelectionSequence
process.pfParticleSelectionSequence = pfParticleSelectionSequence

process.patMuons.isoDeposits = cms.PSet(
    pfAllParticles   = cms.InputTag("muPFIsoDepositPUPFIso"),      # all PU   CH+MU+E
    pfChargedHadrons = cms.InputTag("muPFIsoDepositChargedPFIso"), # all noPU CH
    pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutralPFIso"), # all NH
    pfPhotons        = cms.InputTag("muPFIsoDepositGammaPFIso"),   # all PH
    user = cms.VInputTag(
    cms.InputTag("muPFIsoDepositChargedAllPFIso"),                 # all noPU CH+MU+E
    )
    )
process.patMuons.isolationValues = cms.PSet(
    pfAllParticles   = cms.InputTag("muPFIsoValuePU04PFIso"),
    pfChargedHadrons = cms.InputTag("muPFIsoValueCharged04PFIso"),
    pfNeutralHadrons = cms.InputTag("muPFIsoValueGamma04PFIso"),
    pfPhotons        = cms.InputTag("muPFIsoValueNeutral04PFIso"),
    user = cms.VInputTag(
    cms.InputTag("muPFIsoValueChargedAll04PFIso"),
    )
    )

process.patElectrons.isoDeposits = cms.PSet(
    pfAllParticles   = cms.InputTag("elPFIsoDepositPUPFIso"),      # all PU   CH+MU+E
    pfChargedHadrons = cms.InputTag("elPFIsoDepositChargedPFIso"), # all noPU CH
    pfNeutralHadrons = cms.InputTag("elPFIsoDepositNeutralPFIso"), # all NH
    pfPhotons        = cms.InputTag("elPFIsoDepositGammaPFIso"),   # all PH
    user = cms.VInputTag(
    cms.InputTag("elPFIsoDepositChargedAllPFIso"),                 # all noPU CH+MU+E
    )
    )
process.patElectrons.isolationValues = cms.PSet(
    pfAllParticles   = cms.InputTag("elPFIsoValuePU04PFIdPFIso"),
    pfChargedHadrons = cms.InputTag("elPFIsoValueCharged04PFIdPFIso"),
    pfNeutralHadrons = cms.InputTag("elPFIsoValueGamma04PFIdPFIso"),
    pfPhotons        = cms.InputTag("elPFIsoValueNeutral04PFIdPFIso"),
    user = cms.VInputTag(
    cms.InputTag("elPFIsoValueChargedAll04PFIdPFIso"),
    cms.InputTag("elPFIsoValueChargedAll04NoPFIdPFIso"),
    cms.InputTag("elPFIsoValuePU04NoPFIdPFIso"),
    cms.InputTag("elPFIsoValueCharged04NoPFIdPFIso"),
    cms.InputTag("elPFIsoValueGamma04NoPFIdPFIso"),
    cms.InputTag("elPFIsoValueNeutral04NoPFIdPFIso")
    )
    )

########################  pat::muon  #############################

getattr(process,"patMuons").embedTrack = True

######################## pat::electron ###########################

getattr(process,"patElectrons").embedTrack = True
getattr(process,"patElectrons").embedGsfTrack = True

import RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi as cic
process.eidCiCVeryLoose = cic.eidVeryLoose.clone()
process.eidCiCLoose     = cic.eidLoose.clone()
process.eidCiCMedium    = cic.eidMedium.clone()
process.eidCiCTight     = cic.eidTight.clone() 

import RecoEgamma.ElectronIdentification.cutsInCategoriesHZZElectronIdentificationV06_cfi as cicHZZ
process.eidCiCHZZVeryLoose   = cicHZZ.eidHZZVeryLoose.clone()
process.eidCiCHZZLoose       = cicHZZ.eidHZZLoose.clone()
process.eidCiCHZZMedium      = cicHZZ.eidHZZMedium.clone()
process.eidCiCHZZTight       = cicHZZ.eidHZZTight.clone()
process.eidCiCHZZSuperTight  = cicHZZ.eidHZZSuperTight.clone()
process.eidCiCHZZHyperTight1 = cicHZZ.eidHZZHyperTight1.clone()

process.patElectrons.addElectronID = True
process.patElectrons.electronIDSources.eidCiCVeryLoose      = cms.InputTag("eidCiCVeryLoose")
process.patElectrons.electronIDSources.eidCiCLoose          = cms.InputTag("eidCiCLoose")
process.patElectrons.electronIDSources.eidCiCMedium         = cms.InputTag("eidCiCMedium")
process.patElectrons.electronIDSources.eidCiCTight          = cms.InputTag("eidCiCTight")
process.patElectrons.electronIDSources.eidCiCHZZVeryLoose   = cms.InputTag("eidCiCHZZVeryLoose")
process.patElectrons.electronIDSources.eidCiCHZZLoose       = cms.InputTag("eidCiCHZZLoose")
process.patElectrons.electronIDSources.eidCiCHZZMedium      = cms.InputTag("eidCiCHZZMedium")
process.patElectrons.electronIDSources.eidCiCHZZTight       = cms.InputTag("eidCiCHZZTight")
process.patElectrons.electronIDSources.eidCiCHZZSuperTight  = cms.InputTag("eidCiCHZZSuperTight")
process.patElectrons.electronIDSources.eidCiCHZZHyperTight1 = cms.InputTag("eidCiCHZZHyperTight1")

process.cicIDSequence = cms.Sequence(
    process.eidCiCVeryLoose +
    process.eidCiCLoose +
    process.eidCiCMedium +
    process.eidCiCTight + 
    process.eidCiCHZZVeryLoose +
    process.eidCiCHZZLoose +
    process.eidCiCHZZMedium +
    process.eidCiCHZZTight +
    process.eidCiCHZZSuperTight +
    process.eidCiCHZZHyperTight1
)
process.patDefaultSequence.replace(process.patElectrons,
                                   process.cicIDSequence+process.patElectrons)

######################## pat::tau ################################


######################## pat::jet ################################

getattr(process,"selectedPatJets").cut = cms.string('pt>10 && abs(eta)<5.0')

######################## pat::trigger ############################

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )
process.patTriggerEvent.processName = '*'

if hasattr(process,"patTrigger"):
    process.patTrigger.processName  = '*'


######################## embedding ###############################

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
                 " (pt>=20 && ("    + \
                 "               (isEB && userFloat('sihih')<0.010 && userFloat('dPhi')<0.06 && "  + \
                 "                        userFloat('dEta')< 0.004 && userFloat('HoE') <0.04)"     + \
                 "               ||"+ \
                 "               (isEE && userFloat('sihih')<0.030 && userFloat('dPhi')<0.030 && " + \
                 "                        userFloat('dEta') <0.007 && userFloat('HoE') <0.025) )) "+ \
                 "     || "+ \
                 " (pt<20 && (fbrem>0.15 || (abs(superClusterPosition.Eta)<1. && eSuperClusterOverP>0.95) ) && ( "+ \
                 "               (isEB && userFloat('sihih')<0.010 && userFloat('dPhi')<0.030 && " + \
                 "                        userFloat('dEta') <0.004 && userFloat('HoE') <0.025) "   + \
                 "               ||"+ \
                 "               (isEE && userFloat('sihih')<0.030 && userFloat('dPhi')<0.020 &&"  + \
                 "                        userFloat('dEta') <0.005 && userFloat('HoE') <0.025) ))" + \
                 "    )"   + \
                 ")"

CiCTight = "((electronID('eidCiCTight')==1  ||" + \
           "  electronID('eidCiCTight')==3  ||" + \
           "  electronID('eidCiCTight')==5  ||" + \
           "  electronID('eidCiCTight')==7  ||" + \
           "  electronID('eidCiCTight')==9  ||" + \
           "  electronID('eidCiCTight')==11 ||" + \
           "  electronID('eidCiCTight')==13 ||" + \
           "  electronID('eidCiCTight')==15)" + \
           " || " + \
           " (electronID('eidCiCHZZTight')==1  ||" + \
           "  electronID('eidCiCHZZTight')==3  ||" + \
           "  electronID('eidCiCHZZTight')==5  ||" + \
           "  electronID('eidCiCHZZTight')==7  ||" + \
           "  electronID('eidCiCHZZTight')==9  ||" + \
           "  electronID('eidCiCHZZTight')==11 ||" + \
           "  electronID('eidCiCHZZTight')==13 ||" + \
           "  electronID('eidCiCHZZTight')==15)"+ \
           ")"

MVA = "((pt<=20 && abs(superClusterPosition.Eta)>=0.0 && abs(superClusterPosition.Eta)<1.0 && userFloat('mva')>0.133) ||" + \
      "(pt<=20 && abs(superClusterPosition.Eta)>=1.0 && abs(superClusterPosition.Eta)<1.5 && userFloat('mva')>0.465) ||" + \
      "(pt<=20 && abs(superClusterPosition.Eta)>=1.5 && abs(superClusterPosition.Eta)<2.5 && userFloat('mva')>0.518) ||" + \
      "(pt>20  && abs(superClusterPosition.Eta)>=0.0 && abs(superClusterPosition.Eta)<1.0 && userFloat('mva')>0.942) ||" + \
      "(pt>20  && abs(superClusterPosition.Eta)>=1.0 && abs(superClusterPosition.Eta)<1.5 && userFloat('mva')>0.947) ||" + \
      "(pt>20  && abs(superClusterPosition.Eta)>=1.5 && abs(superClusterPosition.Eta)<2.5 && userFloat('mva')>0.878))"

process.selectedPatMuonsUserEmbedded = cms.EDProducer(
    "MuonsUserEmbedded",
    muonTag           = cms.InputTag("selectedPatMuons"),
    vertexTag         = cms.InputTag("offlinePrimaryVertices"),
    fitUnbiasedVertex = cms.bool(False)
    )

process.selectedPatElectronsUserEmbedded = cms.EDProducer(
    "ElectronsUserEmbedded",
    electronTag = cms.InputTag("selectedPatElectrons"),
    vertexTag   = cms.InputTag("offlinePrimaryVertices"),
    isMC        = cms.bool(runOnMC),
    doMVAMIT    = cms.bool(True),
    doMVADaniele= cms.bool(True),
    inputFileName0 = cms.FileInPath('UserCode/MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_NoIPInfo_BDTG.weights.xml'),
    inputFileName1 = cms.FileInPath('UserCode/MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_NoIPInfo_BDTG.weights.xml'),
    inputFileName2 = cms.FileInPath('UserCode/MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_NoIPInfo_BDTG.weights.xml'),
    inputFileName3 = cms.FileInPath('UserCode/MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_NoIPInfo_BDTG.weights.xml'),
    inputFileName4 = cms.FileInPath('UserCode/MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_NoIPInfo_BDTG.weights.xml'),
    inputFileName5 = cms.FileInPath('UserCode/MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_NoIPInfo_BDTG.weights.xml'),
    inputFileName0v2 = cms.FileInPath('Bianchi/Utilities/data/mvaEleId/Electrons_BDTG_TrigV0_Cat1.weights.xml'),
    inputFileName1v2 = cms.FileInPath('Bianchi/Utilities/data/mvaEleId/Electrons_BDTG_TrigV0_Cat2.weights.xml'),
    inputFileName2v2 = cms.FileInPath('Bianchi/Utilities/data/mvaEleId/Electrons_BDTG_TrigV0_Cat3.weights.xml'),
    inputFileName3v2 = cms.FileInPath('Bianchi/Utilities/data/mvaEleId/Electrons_BDTG_TrigV0_Cat4.weights.xml'),
    inputFileName4v2 = cms.FileInPath('Bianchi/Utilities/data/mvaEleId/Electrons_BDTG_TrigV0_Cat5.weights.xml'),
    inputFileName5v2 = cms.FileInPath('Bianchi/Utilities/data/mvaEleId/Electrons_BDTG_TrigV0_Cat6.weights.xml'),
    #inputFileNameMVADaniele = cms.FileInPath('Bianchi/Utilities/data/mvaEleId/TMVA_BDTSimpleCat.weights.xml')
    )

process.selectedPatTausUserEmbedded = cms.EDProducer(
    "TausUserEmbedded",
    tauTag = cms.InputTag("selectedPatTaus"),
    vertexTag = cms.InputTag("offlinePrimaryVertices"),
    )


####################### pairing ##################################

process.atLeastOneElecTau = cms.EDProducer(
    "CandViewShallowCloneCombiner",
    decay = cms.string("selectedPatElectronsUserEmbedded selectedPatTausUserEmbedded"),
    cut = cms.string("sqrt((daughter(0).eta-daughter(1).eta)*(daughter(0).eta-daughter(1).eta)+  min( abs(daughter(0).phi-daughter(1).phi), 2*3.1415926 - abs(daughter(0).phi-daughter(1).phi)  ) *  min( abs(daughter(0).phi-daughter(1).phi), 2*3.1415926 - abs(daughter(0).phi-daughter(1).phi)  )  )>0.5"),
    checkCharge = cms.bool(False)
    )

process.atLeastOneElecTauCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneElecTau"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

process.elecPtEta = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("selectedPatElectronsUserEmbedded"),
    cut = cms.string("pt>14 && abs(eta)<2.5"),
    filter = cms.bool(False)
    )
process.atLeastOneElecTauelecPtEta = process.atLeastOneElecTau.clone(
    decay=cms.string("elecPtEta selectedPatTausUserEmbedded")
    )
process.elecPtEtaCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneElecTauelecPtEta"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

process.elecPtEtaID = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("selectedPatElectronsUserEmbedded"),
    cut = cms.string(process.elecPtEta.cut.value()+
                     " && abs(userFloat('dxyWrtPV'))<0.045 && abs(userFloat('dzWrtPV'))<0.2 &&"
                     #+"("+simpleCutsWP80+" || "+CiCTight+" || "+MVA+")"),
                     #+simpleCutsWP95
                     +" (userInt('isTriggerElectron')>0.5 || userInt('mvaPreselection')>0.5)"),
    filter = cms.bool(False)
    )
process.atLeastOneElecTauelecPtEtaID = process.atLeastOneElecTau.clone(
    decay=cms.string("elecPtEtaID selectedPatTausUserEmbedded")
    )
process.elecPtEtaIDCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneElecTauelecPtEtaID"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

process.elecPtEtaRelID = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("selectedPatElectronsUserEmbedded"),
    cut = cms.string("pt>14 && abs(eta)<2.5 && "+
                     simpleCutsWP95),
    filter = cms.bool(False)
    )

process.tauPtEta  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("selectedPatTausUserEmbedded"),
    cut = cms.string("pt>19 && abs(eta)<2.3"),
    filter = cms.bool(False)
    )
process.atLeastOneElecTautauPtEta = process.atLeastOneElecTau.clone(
    decay=cms.string("elecPtEtaID tauPtEta")
    )
process.tauPtEtaCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneElecTautauPtEta"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )


process.tauPtEtaID  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("selectedPatTausUserEmbedded"),
    cut = cms.string(process.tauPtEta.cut.value()+
                     " && tauID('decayModeFinding')>0.5"+
                     " && userFloat('dzWrtPV')<0.2"),
    filter = cms.bool(False)
    )
process.atLeastOneElecTautauPtEtaID = process.atLeastOneElecTau.clone(
    decay=cms.string("elecPtEtaID tauPtEtaID")
    )
process.tauPtEtaIDCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneElecTautauPtEtaID"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

process.tauPtEtaIDAgMu  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("selectedPatTausUserEmbedded"),
    cut = cms.string(process.tauPtEtaID.cut.value()+
                     " && tauID('againstMuonLoose')>0.5"),
    filter = cms.bool(False)
    )
process.atLeastOneElecTautauPtEtaIDAgMu = process.atLeastOneElecTau.clone(
    decay=cms.string("elecPtEtaID tauPtEtaIDAgMu")
    )
process.tauPtEtaIDAgMuCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneElecTautauPtEtaIDAgMu"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

process.tauPtEtaIDAgMuAgElec  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("selectedPatTausUserEmbedded"),
    cut = cms.string(process.tauPtEtaIDAgMu.cut.value()+
                     " && tauID('againstElectronMedium')>0.5"),
    filter = cms.bool(False)
    )
process.atLeastOneElecTautauPtEtaIDAgMuAgElec = process.atLeastOneElecTau.clone(
    decay=cms.string("elecPtEtaID tauPtEtaIDAgMuAgElec")
    )
process.tauPtEtaIDAgMuAgElecCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneElecTautauPtEtaIDAgMuAgElec"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

###################### muons ####################################

process.muPtEtaRelID = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsUserEmbedded"),
    cut = cms.string("pt>15 && abs(eta)<2.4 && isGlobalMuon"),
    filter = cms.bool(False)
    )


###################### final sequences ##############################

process.atLeastOneGoodVertexSequence = cms.Sequence(
    process.selectedPrimaryVertices*
    process.primaryVertexCounter*
    process.primaryVertexFilter
    )

process.alLeastOneElecTauSequence = cms.Sequence(
    process.atLeastOneElecTau*process.atLeastOneElecTauCounter*process.atLeastOneElecTauFilter
    )

process.elecLegSequence = cms.Sequence(
    (process.elecPtEta*process.atLeastOneElecTauelecPtEta*process.elecPtEtaCounter*process.elecPtEtaFilter) *
    (process.elecPtEtaID*process.atLeastOneElecTauelecPtEtaID*process.elecPtEtaIDCounter*process.elecPtEtaIDFilter) *
    process.elecPtEtaRelID
    )
process.tauLegSequence = cms.Sequence(
    (process.tauPtEta*process.atLeastOneElecTautauPtEta*process.tauPtEtaCounter*process.tauPtEtaFilter) *
    (process.tauPtEtaID*process.atLeastOneElecTautauPtEtaID*process.tauPtEtaIDCounter*process.tauPtEtaIDFilter) *
    (process.tauPtEtaIDAgMu*process.atLeastOneElecTautauPtEtaIDAgMu*process.tauPtEtaIDAgMuCounter*process.tauPtEtaIDAgMuFilter)*
    (process.tauPtEtaIDAgMuAgElec*process.atLeastOneElecTautauPtEtaIDAgMuAgElec*process.tauPtEtaIDAgMuAgElecCounter*process.tauPtEtaIDAgMuAgElecFilter)
    )

####################### x-cleaning of jets #########################

process.deltaRJetElectrons = cms.EDProducer(
    "DeltaRNearestElectronComputer",
    probes = cms.InputTag("selectedPatJets"),
    objects = cms.InputTag("elecPtEtaID"),
    )
process.selectedPatJetsNoElectrons = cms.EDProducer(
    "JetsCleaner",
    jets =  cms.InputTag("selectedPatJets"),
    valueMap = cms.InputTag("deltaRJetElectrons"),
    minDeltaR = cms.double(0.5)
    )

process.deltaRJetTaus = cms.EDProducer(
    "DeltaRNearestTauComputer",
    probes = cms.InputTag("selectedPatJetsNoElectrons"),
    objects = cms.InputTag("tauPtEtaIDAgMuAgElec"),
    )
process.selectedPatJetsNoElectronsNoTaus = cms.EDProducer(
    "JetsCleaner",
    jets =  cms.InputTag("selectedPatJetsNoElectrons"),
    valueMap = cms.InputTag("deltaRJetTaus"),
    minDeltaR = cms.double(0.5)
    )

process.jetCleaningSequence = cms.Sequence(
    process.deltaRJetElectrons*process.selectedPatJetsNoElectrons*
    process.deltaRJetTaus*process.selectedPatJetsNoElectronsNoTaus
    )


########################## path ###############################

process.skim = cms.Sequence(
    process.allEventsFilter+
    process.atLeastOneGoodVertexSequence*
    process.PFTau*
    #process.fjSequence*
    process.pfParticleSelectionSequence*
    process.muIsoSequence*
    process.electronIsoSequence*
    process.patDefaultSequence*
    #process.kt6PFJetsNeutral*
    process.selectedPatMuonsUserEmbedded*
    process.muPtEtaRelID*
    process.selectedPatElectronsUserEmbedded*
    process.selectedPatTausUserEmbedded*
    process.alLeastOneElecTauSequence*
    process.elecLegSequence*
    process.tauLegSequence#*
    #process.jetCleaningSequence*
    #process.printTree1
    )

massSearchReplaceAnyInputTag(process.skim,
                             "offlinePrimaryVertices",
                             "selectedPrimaryVertices",
                             verbose=False)
process.selectedPrimaryVertices.src = cms.InputTag('offlinePrimaryVertices')

if not runOnMC:
    process.skim.remove(process.printTree1)

process.p = cms.Path(process.skim)

########################## output ###############################

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
process.out.outputCommands = cms.untracked.vstring('drop *',
                                                   *patEventContentNoCleaning ) 
process.out.outputCommands.extend( cms.vstring(
    'keep *_TriggerResults_*_*',
    'keep *_hltTriggerSummaryAOD_*_*',
    'keep recoGenParticles_genParticles*_*_*',
    'keep *_patTriggerEvent_*_*',
    'keep *_patTrigger_*_*',
    'keep *_selectedPatJets_*_*',
    'keep *_ak5PFJets_*_*',
    'keep *_particleFlow__*',
    'keep *_genParticles_*_*',
    'keep *_offlinePrimaryVertices_*_*',
    'keep *_selectedPrimaryVertices_*_*',
    'keep *_offlinePrimaryVerticesWithBS_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_patMETsPFlow_*_*',
    'keep *_tauGenJetsSelectorAllHadrons_*_*',
    'keep *_kt6PFJetsCentral_rho_*',
    'keep *_kt6PFJetsNeutral_rho_*',
    'keep *_elecPtEtaID_*_*',
    'keep *_elecPtEtaRelID_*_*',
    'keep *_muPtEtaRelID_*_*',
    'keep *_muons_*_*',
    'keep *_addPileupInfo_*_*',
    'keep *_generalTracks_*_*',
    'keep *_tmfTracks_*_*',
    'keep *_electronGsfTracks_*_*',
    'keep recoTrackExtras_*_*_*',
    'keep recoGsfTrackExtras_*_*_*',
    'keep *_tauPtEtaIDAgMuAgElec_*_*',
    'keep *_generator_*_*',
    'drop *_TriggerResults_*_HLT',
    'drop *_TriggerResults_*_RECO',
    'drop *_selectedPatElectrons_*_*',
    'drop *_selectedPatMuons_*_*',
    'drop *_selectedPatTaus_*_*',
    'drop *_patMETs_*_*',
    'drop *_selectedPatElectrons_*_*',
    'drop *_selectedPatMuons_*_*',
    'drop *_selectedPatTaus_*_*',
    'drop *_selectedPatElectronsUserEmbedded_*_*',
    'drop *_selectedPatTausUserEmbedded_*_*',
    )
                                   )

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("skimElecTauStream.root")
    )

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

process.out.fileName = cms.untracked.string('patTuples_ElecTauStream.root')

process.outpath = cms.EndPath(process.out)


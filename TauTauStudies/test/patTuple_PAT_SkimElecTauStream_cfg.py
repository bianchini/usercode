from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load("JetMETCorrections.Configuration.JetCorrectionServices_cff")

postfix     = "PFlow"
runOnMC     = True
runOnEmbed  = False

from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

if runOnMC:
    process.GlobalTag.globaltag = cms.string('START42_V17::All')

else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V23::All')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source.fileNames = cms.untracked.vstring(
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat//store/mc/Summer11/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/AODSIM/PU_S4_START42_V11-v1/0000/0E47FBF8-0295-E011-818F-0030487E3026.root'
    #'root://polgrid4.in2p3.fr//dpm/in2p3.fr/home/cms/trivcat//store/mc/Fall11/VBF_HToTauTau_M-115_7TeV-powheg-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/F4ACA82D-FDF8-E011-A31A-E0CB4E29C51E.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/results/higgs/DoubleMu/StoreResults-DoubleMu_2011A_Aug05thRR_v1_embedded_trans1_tau115_ptelec1_17had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/DoubleMu/USER/StoreResults-DoubleMu_2011A_Aug05thRR_v1_embedded_trans1_tau115_ptelec1_17had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/0000/82A9FCE3-A8F9-E011-8974-00266CFCC618.root'
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/results/higgs/DoubleMu/StoreResults-DoubleMu_2011B_PR_v1_embedded_trans1_tau115_ptelec1_17had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/DoubleMu/USER/StoreResults-DoubleMu_2011B_PR_v1_embedded_trans1_tau115_ptelec1_17had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/0000/0011BCB7-8800-E111-85F3-0023AEFDE9AC.root'
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/embeddedTest/embedded_1_1_Zo4.root'
    #'file:pickevents.root'
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/results/higgs/DoubleMu/StoreResults-DoubleMu_2011A_PR_v6_embedded_trans1_tau115_ptelec1_17had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/DoubleMu/USER/StoreResults-DoubleMu_2011A_PR_v6_embedded_trans1_tau115_ptelec1_17had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/0000/F8AFEED6-5FFA-E011-A8F8-0023AEFDEEE0.root'
    #'file:./root/GluGluToHToTauTau_M-125_7TeV-powheg-pythia6_12628F24-31FB-E011-883A-90E6BA19A248.root',
    'file:./root/VBF_HToTauTau_M-125_7TeV-powheg-pythia6-tauola_668A54D7-53F8-E011-9D81-E0CB4E29C502.root',
    )

#process.source.eventsToProcess = cms.untracked.VEventRange(
#    '1:1080'
#    )


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

################### HLT trigger  ####################

process.HLTFilter = cms.EDFilter(
    "HLTHighLevel",
    TriggerResultsTag  = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths           = cms.vstring("HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1"),                                     
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True),
    throw              = cms.bool(False)
    )

################### gen listing  ####################

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree1 = cms.EDAnalyzer(
    "ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint  = cms.untracked.int32(1)
    )

################### jet sequence ####################

process.load('RecoJets.Configuration.RecoPFJets_cff')

process.kt6PFJetsForRhoComputationVoronoi = process.kt6PFJets.clone(
    doRhoFastjet = True,
    voronoiRfact = 0.9
    )

process.kt6PFJets.doRhoFastjet  = True
process.kt6PFJets.doAreaFastjet = True
process.kt6PFJets.Rho_EtaMax    = cms.double(4.4)
process.kt6PFJets.Ghost_EtaMax  = cms.double(5.0)
process.ak5PFJets.doAreaFastjet = True

## re-run kt4PFJets within lepton acceptance to compute rho
process.load('RecoJets.JetProducers.kt4PFJets_cfi')

process.kt6PFJetsCentral = process.kt6PFJets.clone(
    rParam = 0.6,
    doRhoFastjet = True )
process.kt6PFJetsCentral.Rho_EtaMax   = cms.double(1.9)
process.kt6PFJetsCentral.Ghost_EtaMax = cms.double(2.5)

#process.kt6PFJetsNeutral = process.kt4PFJets.clone(
#    rParam = 0.6,
#    doRhoFastjet = True,
#    src="pfAllNeutral" )
#process.kt6PFJetsNeutral.Rho_EtaMax   = cms.double(1.9)
#process.kt6PFJetsNeutral.Ghost_EtaMax = cms.double(2.5)

process.fjSequence = cms.Sequence(process.kt6PFJets+
                                  process.ak5PFJets+
                                  process.kt6PFJetsCentral
                                  #+process.kt6PFJetsForRhoComputationVoronoi
                                  )

# load the PU JetID sequence
process.load("CMGTools.External.pujetidsequence_cff")

################### met ################################

process.load("RecoMET.METProducers.mvaPFMET_cff")
if runOnMC:
    process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3")
else:
    process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3Residual") 

process.pfMEtMVA.srcLeptons = cms.VInputTag( cms.InputTag('elecPtEtaRelIDRelIso'), cms.InputTag('tauPtEtaIDAgMuAgElecRelIso') )

process.patPFMetByMVA = process.patMETs.clone(
    metSource = cms.InputTag('pfMEtMVA'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue')
    )

process.patPFMetByMVA.addGenMET = cms.bool(runOnMC)

################### bTag ##############################

if runOnEmbed:
    process.load('RecoBTag/Configuration/RecoBTag_cff')
    process.load('RecoJets/JetAssociationProducers/ak5JTA_cff')
    process.ak5JetTracksAssociatorAtVertex.jets   = cms.InputTag("ak5PFJets")
    process.ak5JetTracksAssociatorAtVertex.tracks = cms.InputTag("tmfTracks")

## Plus, add this to your path:
#process.ak5JetTracksAssociatorAtVertex*process.btagging

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

from PhysicsTools.PatAlgos.tools.coreTools import *

if not runOnMC:
    removeMCMatching(process,["All"])
    
removeSpecificPATObjects(process, ['Photons'],
                         outputInProcess=None)
removeCleaning(process,
               outputInProcess=None)

restrictInputToAOD(process, ['All'])

from Bianchi.Utilities.customizePAT  import *
addSelectedPFlowParticle(process)

from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, postfix)

from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process,cms.InputTag('ak5PFJets'),
                    doJTA        = True,
                    doBTagging   = True,
                    jetCorrLabel = ('AK5PF', ['L2Relative', 'L3Absolute']),
                    doType1MET   = False,
                    genJetCollection=cms.InputTag("ak5GenJets"),
                    doJetID      = True,
                    jetIdLabel   = 'ak5'
                    )

JEClevels = cms.vstring(['L2Relative', 'L3Absolute'])
if runOnMC:
    JEClevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
else:
    JEClevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']

process.patJetCorrFactors.levels = JEClevels
process.patJetCorrFactors.rho = cms.InputTag('kt6PFJets','rho')
process.patJetCorrFactors.useRho = True

process.patJetCorrFactorsL1Offset = process.patJetCorrFactors.clone(
    levels = cms.vstring('L1Offset',
                         'L2Relative',
                         'L3Absolute')
    )

if runOnMC:
    process.patJetCorrFactorsL1Offset.levels = ['L1Offset', 'L2Relative', 'L3Absolute']
else:
    process.patJetCorrFactorsL1Offset.levels = ['L1Offset', 'L2Relative', 'L3Absolute', 'L2L3Residual']
    
process.patJets.jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactors"),cms.InputTag("patJetCorrFactorsL1Offset"))
process.patDefaultSequence.replace(process.patJetCorrFactors,
                                   process.patJetCorrFactors+process.patJetCorrFactorsL1Offset)



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

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process, 
                 pfTauLabelOld = 'shrinkingConePFTauProducer',
                 pfTauLabelNew = 'hpsPFTauProducer'
                 )

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

#from RecoTauTag.RecoTau.TauDiscriminatorTools import requireLeadTrack
#from RecoTauTag.RecoTau.PFRecoTauDiscriminationByMVAIsolation_cfi import *
#process.hpsPFTauDiscriminationByMVAIsolation = pfRecoTauDiscriminationByMVAIsolation.clone(
#    PFTauProducer = cms.InputTag('hpsPFTauProducer'),
#    Prediscriminants = cms.PSet(
#    BooleanOperator = cms.string('and'),
#    decayMode = cms.PSet(
#    cut = cms.double(0.5),
#    Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding")
#    )
#    ),
#    )
# add to the PFTau sequence
#process.recoTauClassicHPSSequence.replace( process.hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr,
#                                           process.hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr+process.hpsPFTauDiscriminationByMVAIsolation )


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
    pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral04PFIso"),
    pfPhotons        = cms.InputTag("muPFIsoValueGamma04PFIso"),
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
    pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral04PFIdPFIso"),
    pfPhotons        = cms.InputTag("elPFIsoValueGamma04PFIdPFIso"),
    user = cms.VInputTag(
    cms.InputTag("elPFIsoValueChargedAll04PFIdPFIso"),
    cms.InputTag("elPFIsoValueChargedAll04NoPFIdPFIso"),
    cms.InputTag("elPFIsoValuePU04NoPFIdPFIso"),
    cms.InputTag("elPFIsoValueCharged04NoPFIdPFIso"),
    cms.InputTag("elPFIsoValueGamma04NoPFIdPFIso"),
    cms.InputTag("elPFIsoValueNeutral04NoPFIdPFIso")
    )
    )

process.patElectrons.isolationValues = cms.PSet()

########################  pat::muon  #############################

#addPFMuonIsolation(process,process.patMuons)
#addTriggerMatchingMuon(process,isMC=runOnMC)
getattr(process,"patMuons").embedTrack = True


######################## pat::electron ###########################

#addPFElectronIsolation(process,process.patElectrons)
getattr(process,"patElectrons").embedTrack = True
getattr(process,"patElectrons").embedGsfTrack = True
#addTriggerMatchingElectron(process,isMC=runOnMC)

#import RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi as cic
#process.eidCiCVeryLoose = cic.eidVeryLoose.clone()
#process.eidCiCLoose     = cic.eidLoose.clone()
#process.eidCiCMedium    = cic.eidMedium.clone()
#process.eidCiCTight     = cic.eidTight.clone() 

#import RecoEgamma.ElectronIdentification.cutsInCategoriesHZZElectronIdentificationV06_cfi as cicHZZ
#process.eidCiCHZZVeryLoose   = cicHZZ.eidHZZVeryLoose.clone()
#process.eidCiCHZZLoose       = cicHZZ.eidHZZLoose.clone()
#process.eidCiCHZZMedium      = cicHZZ.eidHZZMedium.clone()
#process.eidCiCHZZTight       = cicHZZ.eidHZZTight.clone()
#process.eidCiCHZZSuperTight  = cicHZZ.eidHZZSuperTight.clone()
#process.eidCiCHZZHyperTight1 = cicHZZ.eidHZZHyperTight1.clone()

#process.patElectrons.addElectronID = True
#process.patElectrons.electronIDSources.eidCiCVeryLoose      = cms.InputTag("eidCiCVeryLoose")
#process.patElectrons.electronIDSources.eidCiCLoose          = cms.InputTag("eidCiCLoose")
#process.patElectrons.electronIDSources.eidCiCMedium         = cms.InputTag("eidCiCMedium")
#process.patElectrons.electronIDSources.eidCiCTight          = cms.InputTag("eidCiCTight")
#process.patElectrons.electronIDSources.eidCiCHZZVeryLoose   = cms.InputTag("eidCiCHZZVeryLoose")
#process.patElectrons.electronIDSources.eidCiCHZZLoose       = cms.InputTag("eidCiCHZZLoose")
#process.patElectrons.electronIDSources.eidCiCHZZMedium      = cms.InputTag("eidCiCHZZMedium")
#process.patElectrons.electronIDSources.eidCiCHZZTight       = cms.InputTag("eidCiCHZZTight")
#process.patElectrons.electronIDSources.eidCiCHZZSuperTight  = cms.InputTag("eidCiCHZZSuperTight")
#process.patElectrons.electronIDSources.eidCiCHZZHyperTight1 = cms.InputTag("eidCiCHZZHyperTight1")

#process.cicIDSequence = cms.Sequence(
#    process.eidCiCVeryLoose +
#    process.eidCiCLoose +
#    process.eidCiCMedium +
#    process.eidCiCTight + 
#    process.eidCiCHZZVeryLoose +
#    process.eidCiCHZZLoose +
#    process.eidCiCHZZMedium +
#    process.eidCiCHZZTight +
#    process.eidCiCHZZSuperTight +
#    process.eidCiCHZZHyperTight1
#)
#process.patDefaultSequence.replace(process.patElectrons,
#                                   process.cicIDSequence+process.patElectrons)

######################## pat::tau ################################

addTriggerMatchingTau(process,isMC=runOnMC,postfix="",XtriggerMu=False)


######################## pat::jet ################################

getattr(process,"selectedPatJets").cut = cms.string('pt>10 && abs(eta)<5.0')

######################## pat::trigger ############################

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )
process.patTriggerEvent.processName = '*'

if hasattr(process,"patTrigger"):
    process.patTrigger.processName  = '*'


######################## embedding ###############################

simpleCutsVeto = "(userFloat('nHits')<=999"+ \
                 " && (" + \
                 " (isEB && userFloat('sihih')<0.010 && userFloat('dPhi')<0.80 && "+ \
                 "          userFloat('dEta') <0.007 && userFloat('HoE') <0.15)"   + \
                 " || "  + \
                 " (isEE && userFloat('sihih')<0.030 && userFloat('dPhi')<0.70 && "+ \
                 "          userFloat('dEta') <0.010 && userFloat('HoE') <999)"   + \
                 "     )"+ \
                 ")"
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
    muonTag            = cms.InputTag("selectedPatMuons"),
    vertexTag          = cms.InputTag("offlinePrimaryVertices"),
    fitUnbiasedVertex  = cms.bool(False)
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
    inputFileName0v3 = cms.FileInPath('Bianchi/Utilities/data/mvaEleId/Electrons_BDTG_NonTrigV0_Cat1.weights.xml'),
    inputFileName1v3 = cms.FileInPath('Bianchi/Utilities/data/mvaEleId/Electrons_BDTG_NonTrigV0_Cat2.weights.xml'),
    inputFileName2v3 = cms.FileInPath('Bianchi/Utilities/data/mvaEleId/Electrons_BDTG_NonTrigV0_Cat3.weights.xml'),
    inputFileName3v3 = cms.FileInPath('Bianchi/Utilities/data/mvaEleId/Electrons_BDTG_NonTrigV0_Cat4.weights.xml'),
    inputFileName4v3 = cms.FileInPath('Bianchi/Utilities/data/mvaEleId/Electrons_BDTG_NonTrigV0_Cat5.weights.xml'),
    inputFileName5v3 = cms.FileInPath('Bianchi/Utilities/data/mvaEleId/Electrons_BDTG_NonTrigV0_Cat6.weights.xml'),
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
                     " && abs(userFloat('dxyWrtPV'))<0.045 && abs(userFloat('dzWrtPV'))<0.2"+
                     " && dr03TkSumPt/pt<0.30 &&"+
                     simpleCutsVeto
                     #" && userInt('isTriggerElectron')>0.5"
                     #+"("+simpleCutsWP80+" || "+CiCTight+" || "+MVA+")"),
                     ),
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
    cut = cms.string("pt>14 && abs(eta)<2.5 &&"+
                     "abs(userFloat('dxyWrtPV'))<0.045 && abs(userFloat('dzWrtPV'))<0.2 &&"+
                     simpleCutsVeto
                     ),
    filter = cms.bool(False)
    )

process.elecPtEtaRelIDRelIso = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("selectedPatElectronsUserEmbedded"),
    cut = cms.string(process.elecPtEtaRelID.cut.value()+
                     " && userFloat('PFRelIsoDB04')<0.20"
                     ),
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
                     " && abs(userFloat('dzWrtPV'))<0.2"),
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
                     " && tauID('againstElectronLoose')>0.5"),
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

process.tauPtEtaIDAgMuAgElecRelIso  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("selectedPatTausUserEmbedded"),
    cut = cms.string(process.tauPtEtaIDAgMuAgElec.cut.value()+
                     " && tauID('byVLooseCombinedIsolationDeltaBetaCorr')>0.5"),
    filter = cms.bool(False)
    )

###################### muons ####################################

process.muPtEtaRelID = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsUserEmbedded"),
    cut = cms.string("pt>15 && abs(eta)<2.4 && isGlobalMuon"),
    filter = cms.bool(False)
    )
process.muPtEtaRelIDRelIso = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsUserEmbedded"),
    cut = cms.string(process.muPtEtaRelID.cut.value()+
                     " && userFloat('PFRelIsoDB04')<0.20"),
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
    process.elecPtEtaRelID*process.elecPtEtaRelIDRelIso
    )
process.tauLegSequence = cms.Sequence(
    (process.tauPtEta*process.atLeastOneElecTautauPtEta*process.tauPtEtaCounter*process.tauPtEtaFilter) *
    (process.tauPtEtaID*process.atLeastOneElecTautauPtEtaID*process.tauPtEtaIDCounter*process.tauPtEtaIDFilter) *
    (process.tauPtEtaIDAgMu*process.atLeastOneElecTautauPtEtaIDAgMu*process.tauPtEtaIDAgMuCounter*process.tauPtEtaIDAgMuFilter)*
    (process.tauPtEtaIDAgMuAgElec*process.atLeastOneElecTautauPtEtaIDAgMuAgElec*process.tauPtEtaIDAgMuAgElecCounter*process.tauPtEtaIDAgMuAgElecFilter)*
    process.tauPtEtaIDAgMuAgElecRelIso
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
    process.HLTFilter*
    process.atLeastOneGoodVertexSequence*
    process.fjSequence*
    process.PFTau*
    process.pfParticleSelectionSequence*
    process.muIsoSequence*
    process.electronIsoSequence*
    (process.ak5JetTracksAssociatorAtVertex*process.btagging)*
    process.patDefaultSequence*
    process.puJetIdSqeuence *
    ##process.kt6PFJetsNeutral*
    process.selectedPatMuonsUserEmbedded*
    (process.muPtEtaRelID*process.muPtEtaRelIDRelIso)*
    process.selectedPatElectronsUserEmbedded*
    process.selectedPatTausUserEmbedded*
    process.alLeastOneElecTauSequence*
    process.elecLegSequence*
    process.tauLegSequence*
    (process.pfMEtMVAsequence*process.patPFMetByMVA) +
    ##process.jetCleaningSequence*
    process.printTree1
    )

massSearchReplaceAnyInputTag(process.skim,
                             "offlinePrimaryVertices",
                             "selectedPrimaryVertices",
                             verbose=False)
process.selectedPrimaryVertices.src = cms.InputTag('offlinePrimaryVertices')

if not runOnMC:
    process.skim.remove(process.printTree1)
    process.skim.remove(process.HLTFilter)

if not runOnEmbed:
     process.skim.remove(process.ak5JetTracksAssociatorAtVertex)
     process.skim.remove(process.btagging)

if runOnMC and runOnEmbed:
    process.skim.remove(process.HLTFilter)
    
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
    'keep *_patPFMetByMVA_*_*',
    'keep *_tauGenJetsSelectorAllHadrons_*_*',
    'keep *_kt6PFJets_rho_*',
    'keep *_kt6PFJetsCentral_rho_*',
    'keep *_kt6PFJetsNeutral_rho_*',
    'keep *_kt6PFJetsForRhoComputationVoronoi_rho_*',
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
    'keep *_reducedEcalRecHitsEB_*_*',
    'keep *_reducedEcalRecHitsEE_*_*',
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
    'keep *_puJetId_*_*',
    'keep *_puJetMva_*_*',
    'keep *_hfEMClusters_*_*',
    'keep *_hybridSuperClusters_*_*',
    'keep *_multi5x5BasicClusters_*_*',
    'keep *_pfElectronTranslator_*_*',
    'keep *_pfPhotonTranslator_*_*',
    'keep *_hybridSuperClusters_*_*',
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


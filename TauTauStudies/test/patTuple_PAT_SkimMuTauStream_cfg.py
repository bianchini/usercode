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
    process.GlobalTag.globaltag = cms.string('START42_V17::All') #START42_V14B

else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V23::All') #GR_R_42_V19

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source.fileNames = cms.untracked.vstring(
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/results/higgs/DoubleMu/StoreResults-DoubleMu_2011B_PR_v1_embedded_trans1_tau116_ptmu1_13had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/DoubleMu/USER/StoreResults-DoubleMu_2011B_PR_v1_embedded_trans1_tau116_ptmu1_13had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/0000/FCAE02CE-7800-E111-A2CB-0022198904D4.root'
    #'root://polgrid4.in2p3.fr//dpm/in2p3.fr/home/cms/trivcat//store/mc/Fall11/VBF_HToTauTau_M-115_7TeV-powheg-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/F4ACA82D-FDF8-E011-A31A-E0CB4E29C51E.root',
    #'root://polgrid4.in2p3.fr//dpm/in2p3.fr/home/cms/trivcat/store/user/akalinow/TauPlusX/428_mutau_skim_Run2011A-05Aug2011-v1_v4/b8ede77eca865a3526029ca11820f552/tautauSkimmAOD_9_1_coF.root'
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/akalinow/TauPlusX/428_mutau_skim_Run2011A-05Aug2011-v1_v4/b8ede77eca865a3526029ca11820f552/tautauSkimmAOD_9_1_coF.root'
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/results/higgs/DoubleMu/StoreResults-DoubleMu_2011B_PR_v1_embedded_trans1_tau116_ptmu1_13had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/DoubleMu/USER/StoreResults-DoubleMu_2011B_PR_v1_embedded_trans1_tau116_ptmu1_13had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/0000/FCAE02CE-7800-E111-A2CB-0022198904D4.root'
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/rbonieck/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/428_mutau_Fall11_skimNoTauIsol_v1/c69db8fd1d6ad63518141f89915d27aa/tautauSkimmAOD_87_1_poE.root'
    #'file:./root/pickevents.root',
    #'file:./root/GluGluToHToTauTau_M-125_7TeV-powheg-pythia6_12628F24-31FB-E011-883A-90E6BA19A248.root',
    #'file:./root/VBF_HToTauTau_M-125_7TeV-powheg-pythia6-tauola_668A54D7-53F8-E011-9D81-E0CB4E29C502.root',
    'file:./root/pfembTauTau_TTJets_Fall11_PU_S6_START42_V14B_v2_1_0_pt_0_398_embedded.root'
    )

#process.source.eventsToProcess = cms.untracked.VEventRange(
#    '1:15460'
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
process.atLeastOneMuTauFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.muPtEtaFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.muPtEtaIDFilter = cms.EDFilter(
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
process.atLeast1selectedDiTauFilter = cms.EDFilter(
    "AllEventsFilter"
    )


################### HLT trigger  ####################

process.HLTFilter = cms.EDFilter(
    "HLTHighLevel",
    TriggerResultsTag  = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths           = cms.vstring("HLT_IsoMu15_LooseIsoPFTau15_v9",
                                     "HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v1"),
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
    rParam       = 0.6,
    doRhoFastjet = True )
process.kt6PFJetsCentral.Rho_EtaMax   = cms.double(1.9)
process.kt6PFJetsCentral.Ghost_EtaMax = cms.double(2.5)

#process.kt6PFJetsNeutral = process.kt4PFJets.clone(
#    rParam       = 0.6,
#    doRhoFastjet = True,
#    src          = "pfAllNeutral" )
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

process.pfMEtMVA.srcLeptons = cms.VInputTag( cms.InputTag('muPtEtaRelIDRelIso'), cms.InputTag('tauPtEtaIDAgMuAgElecRelIso') )

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

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag

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
                    jetCorrLabel = ('AK5PF', ['L2Relative', 'L3Absolute',]),
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
process.patJetCorrFactors.rho    = cms.InputTag('kt6PFJets','rho')
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

process.patJets.jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactors"),
                                                     cms.InputTag("patJetCorrFactorsL1Offset"))
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

#switchToPFTauHPS(process)


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

#process.patElectrons.isoDeposits = cms.PSet(
#    pfAllParticles   = cms.InputTag("elPFIsoDepositPUPFIso"),
#    )
process.patElectrons.isolationValues = cms.PSet(
    #pfAllParticles   = cms.InputTag("elPFIsoValuePU04PFIdPFIso"),
    )

########################  pat::muon  #############################

#addPFMuonIsolation(process,process.patMuons)
#addTriggerMatchingMuon(process,isMC=runOnMC)
getattr(process,"patMuons").embedTrack = True


######################## pat::electron ###########################

#addPFElectronIsolation(process,process.patElectrons)
getattr(process,"patElectrons").embedTrack    = True
getattr(process,"patElectrons").embedGsfTrack = True
#addTriggerMatchingElectron(process,isMC=runOnMC)


######################## pat::tau ################################

#addTriggerMatchingTau(process,isMC=runOnMC,postfix="",XtriggerMu=True)


######################## pat::jet ################################

getattr(process,"selectedPatJets").cut = cms.string('pt>10 && abs(eta)<5.0')


######################## pat::trigger ############################

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )
process.patTriggerEvent.processName = '*'

if hasattr(process,"patTrigger"):
    process.patTrigger.processName = '*'
    

######################## embedding ###############################

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
    tauTag    = cms.InputTag("selectedPatTaus"),
    vertexTag = cms.InputTag("offlinePrimaryVertices"),
    )

####################### pairing ##################################

process.atLeastOneMuTau = cms.EDProducer(
    "CandViewShallowCloneCombiner",
    decay = cms.string("selectedPatMuonsUserEmbedded selectedPatTausUserEmbedded"),
    cut = cms.string("sqrt((daughter(0).eta-daughter(1).eta)*(daughter(0).eta-daughter(1).eta)+  min( abs(daughter(0).phi-daughter(1).phi), 2*3.1415926 - abs(daughter(0).phi-daughter(1).phi)  ) *  min( abs(daughter(0).phi-daughter(1).phi), 2*3.1415926 - abs(daughter(0).phi-daughter(1).phi)  )  )>0.5"),
    checkCharge = cms.bool(False)
    )

process.atLeastOneMuTauCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneMuTau"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

process.muPtEta = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsUserEmbedded"),
    cut = cms.string("pt>14 && abs(eta)<2.4"),
    filter = cms.bool(False)
    )
process.atLeastOneMuTaumuPtEta = process.atLeastOneMuTau.clone(
    decay=cms.string("muPtEta selectedPatTausUserEmbedded")
    )
process.muPtEtaCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneMuTaumuPtEta"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

process.muPtEtaRelID = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsUserEmbedded"),
    cut = cms.string("pt>14 && abs(eta)<2.4 && isGlobalMuon"+
                     "&& abs(userFloat('dxyWrtPV'))<0.045 && abs(userFloat('dzWrtPV'))<0.2"
                     ),
    filter = cms.bool(False)
    )

process.muPtEtaRelIDRelIso = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsUserEmbedded"),
    cut = cms.string(process.muPtEtaRelID.cut.value()+
                     " && userFloat('PFRelIsoDB04')<0.20"),
    filter = cms.bool(False)
    )

process.muPtEtaID = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsUserEmbedded"),
    cut = cms.string(process.muPtEta.cut.value()+
                     " && abs(userFloat('dxyWrtPV'))<0.045 && abs(userFloat('dzWrtPV'))<0.2"+
                     " && ("+
                     "(   isGlobalMuon"+
                     " && globalTrack.isNonnull "+
                     " && globalTrack.normalizedChi2<10"+
                     " && globalTrack.hitPattern.numberOfValidMuonHits>0"+                     
                     " && numberOfMatchedStations>1"+                     
                     " && innerTrack.hitPattern.numberOfValidPixelHits>0"+
                     " && track.hitPattern.trackerLayersWithMeasurement > 5)"+
                     " || userInt('isPFMuon')>0.5)"
                     ),
    filter = cms.bool(False)
    )
process.atLeastOneMuTaumuPtEtaID = process.atLeastOneMuTau.clone(
    decay=cms.string("muPtEtaID selectedPatTausUserEmbedded")
    )
process.muPtEtaIDCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneMuTaumuPtEtaID"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

process.tauPtEta  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("selectedPatTausUserEmbedded"),
    cut = cms.string("pt>19 && abs(eta)<2.3"),
    filter = cms.bool(False)
    )
process.atLeastOneMuTautauPtEta = process.atLeastOneMuTau.clone(
    decay=cms.string("muPtEtaID tauPtEta")
    )
process.tauPtEtaCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneMuTautauPtEta"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

process.tauPtEtaID  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("selectedPatTausUserEmbedded"),
    cut = cms.string(process.tauPtEta.cut.value()+
                     " && tauID('decayModeFinding')>0.5"+
                     " && abs(userFloat('dzWrtPV'))<0.2"
                     ),
    filter = cms.bool(False)
    )
process.atLeastOneMuTautauPtEtaID = process.atLeastOneMuTau.clone(
    decay=cms.string("muPtEtaID tauPtEtaID")
    )
process.tauPtEtaIDCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneMuTautauPtEtaID"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

process.tauPtEtaIDAgMu  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("selectedPatTausUserEmbedded"),
    cut = cms.string(process.tauPtEtaID.cut.value()+
                     " && tauID('againstMuonTight')>0.5"),
    filter = cms.bool(False)
    )
process.atLeastOneMuTautauPtEtaIDAgMu = process.atLeastOneMuTau.clone(
    decay=cms.string("muPtEtaID tauPtEtaIDAgMu")
    )
process.tauPtEtaIDAgMuCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneMuTautauPtEtaIDAgMu"),
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
process.atLeastOneMuTautauPtEtaIDAgMuAgElec = process.atLeastOneMuTau.clone(
    decay=cms.string("muPtEtaID tauPtEtaIDAgMuAgElec")
    )
process.tauPtEtaIDAgMuAgElecCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneMuTautauPtEtaIDAgMuAgElec"),
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
###################### electrons ####################################

#https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification

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

process.elecPtEtaRelID = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("selectedPatElectronsUserEmbedded"),
    cut = cms.string("pt>15 && abs(eta)<2.4" +
                     " && abs(userFloat('dxyWrtPV'))<0.045 && abs(userFloat('dzWrtPV'))<0.2 &&"+
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

###################### final sequences ##############################

process.atLeastOneGoodVertexSequence = cms.Sequence(
    process.selectedPrimaryVertices*
    process.primaryVertexCounter*
    process.primaryVertexFilter
    )

process.alLeastOneMuTauSequence = cms.Sequence(
    process.atLeastOneMuTau*
    process.atLeastOneMuTauCounter*
    process.atLeastOneMuTauFilter
    )

process.muLegSequence = cms.Sequence(
    (process.muPtEta*process.atLeastOneMuTaumuPtEta*process.muPtEtaCounter*process.muPtEtaFilter) *
    (process.muPtEtaID*process.atLeastOneMuTaumuPtEtaID*process.muPtEtaIDCounter*process.muPtEtaIDFilter) +
    process.muPtEtaRelID*process.muPtEtaRelIDRelIso
    )

process.tauLegSequence = cms.Sequence(
    (process.tauPtEta*process.atLeastOneMuTautauPtEta*process.tauPtEtaCounter*process.tauPtEtaFilter) *
    (process.tauPtEtaID*process.atLeastOneMuTautauPtEtaID*process.tauPtEtaIDCounter*process.tauPtEtaIDFilter) *
    (process.tauPtEtaIDAgMu*process.atLeastOneMuTautauPtEtaIDAgMu*process.tauPtEtaIDAgMuCounter*process.tauPtEtaIDAgMuFilter)*
    (process.tauPtEtaIDAgMuAgElec*process.atLeastOneMuTautauPtEtaIDAgMuAgElec*process.tauPtEtaIDAgMuAgElecCounter*process.tauPtEtaIDAgMuAgElecFilter)*
    process.tauPtEtaIDAgMuAgElecRelIso
    )

####################### x-cleaning of jets #########################

process.deltaRJetMuons = cms.EDProducer(
    "DeltaRNearestMuonComputer",
    probes = cms.InputTag("selectedPatJets"),
    objects = cms.InputTag("muPtEtaID"),
    )
process.selectedPatJetsNoMuons = cms.EDProducer(
    "JetsCleaner",
    jets =  cms.InputTag("selectedPatJets"),
    valueMap = cms.InputTag("deltaRJetMuons"),
    minDeltaR = cms.double(0.5)
    )

process.deltaRJetTaus = cms.EDProducer(
    "DeltaRNearestTauComputer",
    probes = cms.InputTag("selectedPatJetsNoMuons"),
    objects = cms.InputTag("tauPtEtaIDAgMuAgElec"),
    )
process.selectedPatJetsNoMuonsNoTaus = cms.EDProducer(
    "JetsCleaner",
    jets =  cms.InputTag("selectedPatJetsNoMuons"),
    valueMap = cms.InputTag("deltaRJetTaus"),
    minDeltaR = cms.double(0.3)
    )

process.jetCleaningSequence = cms.Sequence(
    process.deltaRJetMuons*process.selectedPatJetsNoMuons*
    process.deltaRJetTaus*process.selectedPatJetsNoMuonsNoTaus
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
    process.selectedPatElectronsUserEmbedded*
    (process.elecPtEtaRelID+process.elecPtEtaRelIDRelIso)*
    process.selectedPatTausUserEmbedded*
    process.alLeastOneMuTauSequence*
    process.muLegSequence*
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

#process.p = cms.Path(process.printEventContent+process.skim)
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
    'keep *_genParticles_*_*',
    'keep *_particleFlow__*',
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
    'keep *_muPtEtaID_*_*',
    'keep *_muPtEtaRelID_*_*',
    'keep *_muons_*_*',
    'keep *_elecPtEtaRelID_*_*',
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
    'drop *_selectedPatMuons_*_*',
    'drop *_selectedPatElectrons_*_*',
    'drop *_selectedPatTaus_*_*',
    'drop *_selectedPatMuonsUserEmbedded_*_*',
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
    fileName = cms.string("skimMuTauStream.root")
    )

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

process.out.fileName = cms.untracked.string('patTuples_MuTauStream.root')

process.outpath = cms.EndPath(process.out)


processDumpFile = open('patTuplePATSkimMuTauStream.dump', 'w')
print >> processDumpFile, process.dumpPython()

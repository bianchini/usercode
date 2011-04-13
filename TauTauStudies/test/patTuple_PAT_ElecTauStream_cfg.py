from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.Rho_EtaMax = cms.double(4.4)
process.kt6PFJets.Ghost_EtaMax = cms.double(5.0)
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.Rho_EtaMax = cms.double(4.4)
process.ak5PFJets.Ghost_EtaMax = cms.double(5.0)

process.ak5PFL1Fastjet.useCondDB = False

process.fjSequence = cms.Sequence(process.kt6PFJets+process.ak5PFJets)

process.source.fileNames = cms.untracked.vstring(
    #'file:/data_CMS/cms/lbianchini/ZTT_RelVal386_1.root',
    #'file:/data_CMS/cms/lbianchini/ZMuMu_RelVal386.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat//store/mc/Fall10/VBF_HToTauTau_M-115_7TeV-powheg-pythia6-tauola/GEN-SIM-RECO/START38_V12-v1/0000/044E940A-55EC-DF11-89D6-0023AEFDEE60.root',
    #'file:/data_CMS/cms/lbianchini/F41A3437-7AED-DF11-A50D-002618943894.root',
    #'file:/data_CMS/cms/lbianchini/ZElEl_RelVal386_1.root',
    )

#process.source.eventsToProcess = cms.untracked.VEventRange(
#    '1:1080','1:2003','1:2028','1:6867','1:7016'
#    )

postfix           = "PFlow"
sample            = ""
runOnMC           = True


process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree1 = cms.EDAnalyzer("ParticleListDrawer",
                                    src = cms.InputTag("genParticles"),
                                    maxEventsToPrint  = cms.untracked.int32(1)
                                    )

process.primaryVertexFilter = cms.EDFilter(
    "GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
    )

process.scrapping = cms.EDFilter(
    "FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
    )

process.allEventsFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.vertexScrapingFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.oneElectronFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.noMuonFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.electronLegFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.tauLegFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.atLeastOneDiTauFilter = cms.EDFilter(
    "AllEventsFilter"
    )


from PhysicsTools.PatAlgos.tools.coreTools import *

if not runOnMC:
    removeMCMatching(process,["All"])
    
removeSpecificPATObjects(process, ['Photons'],
                         outputInProcess=False)
removeCleaning(process,
               outputInProcess=False)

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
    JEClevels = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']

process.patJetCorrFactors.levels = JEClevels
process.patJetCorrFactors.rho = cms.InputTag('kt6PFJets','rho')


## <tau part>
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process, 
                 pfTauLabelOld = 'shrinkingConePFTauProducer',
                 pfTauLabelNew = 'hpsPFTauProducer'
                 )

getattr(process,"patTaus").embedIsolationTracks = cms.bool(True)
getattr(process,"patTaus").embedSignalTracks = cms.bool(True)
getattr(process,"patTaus").embedGenMatch = cms.bool(True)
getattr(process,"patTaus").embedLeadTrack = cms.bool(True)
getattr(process,"patTaus").embedLeadPFCand = True
getattr(process,"patTaus").embedLeadPFChargedHadrCand = True
getattr(process,"patTaus").embedLeadPFNeutralCand = True
getattr(process,"patTaus").embedSignalPFCands = True
getattr(process,"patTaus").embedSignalPFChargedHadrCands = True
getattr(process,"patTaus").embedSignalPFNeutralHadrCands = True
getattr(process,"patTaus").embedSignalPFGammaCands = True
getattr(process,"patTaus").embedIsolationPFCands = True
getattr(process,"patTaus").embedIsolationPFChargedHadrCands = True
getattr(process,"patTaus").embedIsolationPFNeutralHadrCands = True
getattr(process,"patTaus").embedIsolationPFGammaCands = True
getattr(process,"patTaus").embedGenJetMatch = cms.bool(True)


setattr(process,"hpsPFTauDiscriminationAgainstElectronCrackRem",
        getattr(process,"hpsPFTauDiscriminationAgainstElectron").clone(
    ApplyCut_EcalCrackCut = cms.bool(True),
    ApplyCut_PFElectronMVA =  cms.bool(False)
    )
        )

process.patHPSPFTauDiscrimination += process.hpsPFTauDiscriminationAgainstElectronCrackRem

getattr(process,"makePatTaus").replace(
    getattr(process,"patTaus"),
    process.patHPSPFTauDiscrimination + getattr(process,"patTaus")
    )

getattr(process,"patTaus").tauIDSources = cms.PSet(
    leadingTrackFinding = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
    byLooseIsolation = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"),
    byMediumIsolation = cms.InputTag("hpsPFTauDiscriminationByMediumIsolation"),
    byTightIsolation = cms.InputTag("hpsPFTauDiscriminationByTightIsolation"),
    againstElectron = cms.InputTag("hpsPFTauDiscriminationAgainstElectron"),
    againstElectronCrackRem = cms.InputTag("hpsPFTauDiscriminationAgainstElectronCrackRem"),
    againstMuon = cms.InputTag("hpsPFTauDiscriminationAgainstMuon")
    )

process.tauMatch.maxDeltaR = 0.5
process.tauMatch.resolveAmbiguities = cms.bool(False)
process.tauGenJetMatch.resolveAmbiguities = cms.bool(False)
process.tauGenJetMatch.maxDeltaR = 0.15
process.tauGenJetMatch.maxDPtRel = 999
## <\tau part>

addPFMuonIsolation(process,process.patMuons)
#addPFMuon(process,postfix)
addTriggerMatchingMuon(process)
getattr(process,"patMuons").embedTrack = True
#getattr(process,"patMuons"+postfix).embedTrack = True
#addTriggerMatchingMuon(process,postfix)

from Bianchi.Utilities.electrons import *
addCutBasedID(process)
addPFElectronIsolation(process,process.patElectrons)
#addPFElectron(process,postfix)
getattr(process,"patElectrons").embedTrack = True
#getattr(process,"patElectrons"+postfix).embedTrack = True
getattr(process,"patElectrons").embedGsfTrack = True
#getattr(process,"patElectrons"+postfix).embedGsfTrack = True
addTriggerMatchingElectron(process)
#addTriggerMatchingElectron(process,postfix)

if hasattr(process,"patTrigger"):
    process.patTrigger.processName = '*'

# SC sequence
process.mergedSuperClusters = cms.EDProducer(
    "SuperClusterMerger",
    src = cms.VInputTag(cms.InputTag("correctedHybridSuperClusters"),
                        cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"))
    )
process.selectedSuperClusters = cms.EDFilter(
    "SuperClusterSelector",
    src = cms.InputTag("mergedSuperClusters"),
    cut = cms.string("energy()*sin(position().theta()) > 4"),
    filter = cms.bool(False)
    )
process.makeSCs = cms.Sequence(process.mergedSuperClusters*process.selectedSuperClusters)

    
#########
## PAT

process.selectedPatMuonsTriggerMatchUserEmbedded = cms.EDProducer(
    "MuonsUserEmbedded",
    muonTag = cms.InputTag("selectedPatMuonsTriggerMatch"),
    vertexTag = cms.InputTag("offlinePrimaryVertices")
    )

process.looseMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsTriggerMatchUserEmbedded"),
    cut = cms.string("pt>15 && (eta<2.4&&eta>-2.4) && isGlobalMuon && globalTrack.isNonnull"),
    filter = cms.bool(False)
    )

process.selectedPatElectronsTriggerMatchUserEmbedded = cms.EDProducer(
    "ElectronsUserEmbedded",
    electronTag = cms.InputTag("selectedPatElectronsTriggerMatch"),
    vertexTag = cms.InputTag("offlinePrimaryVertices")
    )

process.looseElectrons = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("selectedPatElectronsTriggerMatchUserEmbedded"),
    cut = cms.string("pt>15 && (eta<2.4&&eta>-2.4)  && !isEBEEGap  && electronID('simpleEleId95relIso')>6.5"),
    filter = cms.bool(False)
    )

process.electronLeg = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("looseElectrons"),
    #cut = cms.string("pt>15 && (eta<2.1&&eta>-2.1) && isTrackerMuon && numberOfMatches>=2 && globalTrack.isNonnull && globalTrack.hitPattern.numberOfValidMuonHits>=1 && globalTrack.hitPattern.numberOfValidPixelHits>=1 && globalTrack.normalizedChi2<=10 && userFloat('dxyWrtPV')<0.2 && userFloat('PFRelIso04')<0.1 && (triggerObjectMatchesByPath('HLT_Mu11').size()!=0 || triggerObjectMatchesByPath('HLT_Mu9').size()!=0 || triggerObjectMatchesByPath('HLT_Mu15').size()!=0 || triggerObjectMatchesByPath('HLT_Mu15_v1').size()!=0)"),
    #cut = cms.string("( (electronID('simpleEleId80relIso')>4.5 && electronID('simpleEleId80relIso')<5.5) || electronID('simpleEleId80relIso')>6.5 ) && userFloat('dxyWrtPV')<0.2 && userFloat('PFRelIso04')<999 && (triggerObjectMatchesByPath('HLT_Ele10_LW_L1R').size()!=0 || triggerObjectMatchesByPath('HLT_Ele15_SW_L1R').size()!=0 || triggerObjectMatchesByPath('HLT_Ele15_SW_CaloEleId_L1R').size()!=0 || triggerObjectMatchesByPath('HLT_Ele17_SW_CaloEleId_L1R').size()!=0)"),
    cut = cms.string("( (electronID('simpleEleId80relIso')>4.5 && electronID('simpleEleId80relIso')<5.5) || electronID('simpleEleId80relIso')>6.5 ) && userFloat('dxyWrtPV')<0.2 && userFloat('PFRelIso04')<999"),
    filter = cms.bool(False)
    )

process.tauLeg = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("selectedPatTaus"),
    #cut = cms.string("pt>20 && (eta<2.3&&eta>-2.3) && tauID('leadingTrackFinding')>0.5 && tauID('byLooseIsolation')>0.5 && tauID('againstMuon') && leadPFChargedHadrCand.mva_e_pi < 0.6 && tauID('againstElectronCrackRem')>0.5 "),
    cut = cms.string("pt>20 && (eta<2.3&&eta>-2.3) && tauID('leadingTrackFinding')>0.5 && tauID('byLooseIsolation')>-1 && tauID('againstMuon') && leadPFChargedHadrCand.mva_e_pi < -0.1 && tauID('againstElectronCrackRem')>0.5 "),
    filter = cms.bool(False)
    )


process.secondLeptonVeto = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("looseElectrons"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1),
    )
process.noMuonVeto = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("looseMuons"),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(0),
    )
process.oneElectronLeg = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("electronLeg"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )
process.oneTauLeg = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("tauLeg"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

process.load("Bianchi.Utilities.diTausReconstruction_cff")
process.diTau = process.elecTauPairs.clone()
process.diTau.srcLeg1 = cms.InputTag("electronLeg")
process.diTau.srcLeg2 = cms.InputTag("tauLeg")
if not runOnMC:
    process.diTau.srcGenParticles = ""


process.selectedDiTau = cms.EDFilter(
    "ElecTauPairSelector",
    src = cms.InputTag("diTau"),
    #cut = cms.string("charge==0 && mt1MET<40")
    cut = cms.string("charge>-99")
    )

process.atLeast1selectedDiTau = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("selectedDiTau"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

getattr(process,"selectedPatJets").cut = cms.string('pt>5 && abs(eta)<4.5')

process.deltaRJetElectrons = cms.EDProducer(
    "DeltaRNearestElectronComputer",
    probes = cms.InputTag("selectedPatJets"),
    objects = cms.InputTag("electronLeg"),
    )
process.selectedPatJetsNoElectrons = cms.EDProducer(
    "JetsCleaner",
    jets =  cms.InputTag("selectedPatJets"),
    valueMap = cms.InputTag("deltaRJetElectrons"),
    minDeltaR = cms.double(0.3)
    )

process.deltaRJetTaus = cms.EDProducer(
    "DeltaRNearestTauComputer",
    probes = cms.InputTag("selectedPatJetsNoElectrons"),
    objects = cms.InputTag("tauLeg"),
    )
process.selectedPatJetsNoElectronsNoTaus = cms.EDProducer(
    "JetsCleaner",
    jets =  cms.InputTag("selectedPatJetsNoElectrons"),
    valueMap = cms.InputTag("deltaRJetTaus"),
    minDeltaR = cms.double(0.3)
    )

process.elecTauStreamAnalyzer = cms.EDAnalyzer(
    "ElecTauStreamAnalyzer",
    diTaus =  cms.InputTag("selectedDiTau"),
    jets =  cms.InputTag("selectedPatJetsNoElectrons"),
    isMC = cms.bool(runOnMC),
    deltaRLegJet  = cms.untracked.double(0.3),
    minCorrPt = cms.untracked.double(5.),
    minJetID  = cms.untracked.double(0.5), # 1=loose,2=medium,3=tight
    applyTauSignalSel =  cms.bool( True ),
    verbose =  cms.untracked.bool( False ),
    )

process.pat = cms.Sequence(
    process.allEventsFilter+
    (process.primaryVertexFilter+process.scrapping)*
    process.vertexScrapingFilter +
    process.makeSCs +
    process.fjSequence*
    process.patDefaultSequence*
    process.selectedPatMuonsTriggerMatchUserEmbedded*
    process.selectedPatElectronsTriggerMatchUserEmbedded*
    process.looseMuons*
    process.looseElectrons*
    (process.secondLeptonVeto*process.oneElectronFilter) +
    (process.noMuonVeto*process.noMuonFilter) +
    (process.electronLeg * process.oneElectronLeg * process.electronLegFilter) +
    (process.tauLeg * process.oneTauLeg * process.tauLegFilter) +
    (process.diTau*process.selectedDiTau*process.atLeast1selectedDiTau*process.atLeastOneDiTauFilter) +
    (process.deltaRJetElectrons*process.selectedPatJetsNoElectrons)*
    (process.deltaRJetTaus*process.selectedPatJetsNoElectronsNoTaus)*
    process.elecTauStreamAnalyzer+
    process.printTree1
    )

if not runOnMC:
    process.pat.remove(process.printTree1)

process.p = cms.Path(process.pat)

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
process.out.outputCommands = cms.untracked.vstring('drop *',
                                                   *patEventContentNoCleaning ) 
process.out.outputCommands.extend( cms.vstring(
    'keep *_TriggerResults_*_*',
    'keep *_hltTriggerSummaryAOD_*_*',
    'keep recoGenParticles_genParticles*_*_*',
    'keep *_generalTracks_*_*',
    'keep *_electronGsfTracks_*_*',
    'keep recoTrackExtras_*_*_*',
    'keep recoGsfTrackExtras_*_*_*',
    #'keep TrackingRecHitsOwned_*_*_*',
    'keep *_selectedSuperClusters_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_offlinePrimaryVertices*_*_*',
    'keep *_particleFlow_*_*',
    'keep *_selectedPatJetsNoMuonsNoTaus_*_*',
    'keep *_selectedDiTau_*_*',
    )
                                   )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("treeElecTauStream.root")
                                   )

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

process.out.fileName = cms.untracked.string('patTuples_'+sample+'.root')

#process.outpath = cms.EndPath(process.out)
process.outpath = cms.EndPath()


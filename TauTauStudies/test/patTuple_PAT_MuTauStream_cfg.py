from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
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
     'rfio:/dpm/in2p3.fr/home/cms/trivcat//store/mc/Spring11/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/AODSIM/PU_S1_START311_V1G1-v2/0000/FA5943AB-A756-E011-A6C8-002618FDA208.root'
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat//store/mc/Fall10/VBF_HToTauTau_M-115_7TeV-powheg-pythia6-tauola/GEN-SIM-RECO/START38_V12-v1/0000/044E940A-55EC-DF11-89D6-0023AEFDEE60.root',
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
process.oneMuonFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.noElecFilter = cms.EDFilter(
    "AllEventsFilter"
    )
process.muonLegFilter = cms.EDFilter(
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
#addSelectedPFlowParticle(process)

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

from RecoTauTag.RecoTau.PFRecoTauDiscriminationAgainstElectron_cfi import pfRecoTauDiscriminationAgainstElectron
process.hpsPFTauDiscriminationAgainstElectronCrackRem = pfRecoTauDiscriminationAgainstElectron.clone(
    PFTauProducer = cms.InputTag('hpsPFTauProducer'),
    Prediscriminants = cms.PSet(
    BooleanOperator = cms.string("and"),
    leadTrack = cms.PSet(
    Producer = cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding'),
    cut = cms.double(0.5)
    )
    ),
    ApplyCut_EcalCrackCut = cms.bool(True),
    ApplyCut_PFElectronMVA =  cms.bool(False)
    )

# add the crack-rem sequence before running patTau
process.patDefaultSequence.replace(process.patTaus,
                                   process.hpsPFTauDiscriminationAgainstElectronCrackRem+process.patTaus)


getattr(process,"patTaus").tauIDSources = cms.PSet(
    leadingTrackFinding = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
    byLooseIsolation = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"),
    byMediumIsolation = cms.InputTag("hpsPFTauDiscriminationByMediumIsolation"),
    byTightIsolation = cms.InputTag("hpsPFTauDiscriminationByTightIsolation"),
    againstElectronLoose = cms.InputTag("hpsPFTauDiscriminationByLooseElectronRejection"),
    againstElectronMedium = cms.InputTag("hpsPFTauDiscriminationByMediumElectronRejection"),
    againstElectronTight = cms.InputTag("hpsPFTauDiscriminationByTightElectronRejection"),
    againstElectronCrackRem = cms.InputTag("hpsPFTauDiscriminationAgainstElectronCrackRem"),
    againstMuonLoose = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection"),
    againstMuonTight = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection")
    )


process.tauMatch.maxDeltaR = 0.5
process.tauMatch.resolveAmbiguities = cms.bool(False)
process.tauGenJetMatch.resolveAmbiguities = cms.bool(False)
process.tauGenJetMatch.maxDeltaR = 0.15
process.tauGenJetMatch.maxDPtRel = 999
## <\tau part>

addPFMuonIsolation(process,process.patMuons)
addTriggerMatchingMuon(process)
getattr(process,"patMuons").embedTrack = True

from Bianchi.Utilities.electrons import *
addCutBasedID(process)
addPFElectronIsolation(process,process.patElectrons)

getattr(process,"patElectrons").embedTrack = True
getattr(process,"patElectrons").embedGsfTrack = True
addTriggerMatchingElectron(process)

addTriggerMatchingTau(process)


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

process.muonLeg = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsTriggerMatchUserEmbedded"),
    #cut = cms.string("pt>15 && (eta<2.1&&eta>-2.1) && isTrackerMuon && numberOfMatches>=2 && globalTrack.isNonnull && globalTrack.hitPattern.numberOfValidMuonHits>=1 && globalTrack.hitPattern.numberOfValidPixelHits>=1 && globalTrack.normalizedChi2<=10 && userFloat('dxyWrtPV')<0.2 && userFloat('PFRelIso04')<0.1 && (triggerObjectMatchesByPath('HLT_Mu11').size()!=0 || triggerObjectMatchesByPath('HLT_Mu9').size()!=0 || triggerObjectMatchesByPath('HLT_Mu15').size()!=0 || triggerObjectMatchesByPath('HLT_Mu15_v1').size()!=0)"),
    cut = cms.string("pt>15 && (eta<2.1&&eta>-2.1) && isTrackerMuon && numberOfMatches>=2 && globalTrack.isNonnull && globalTrack.hitPattern.numberOfValidMuonHits>=1 && globalTrack.hitPattern.numberOfValidPixelHits>=1 && globalTrack.normalizedChi2<=10 && userFloat('dxyWrtPV')<0.2 && userFloat('PFRelIso04')<999"),
    filter = cms.bool(False)
    )

process.tauLeg = cms.EDFilter(
    "PATTauSelector",
    #src = cms.InputTag("selectedPatTaus"),
    src = cms.InputTag("selectedPatTausTriggerMatch"),
    #cut = cms.string("pt>20 && (eta<2.3&&eta>-2.3) && tauID('leadingTrackFinding')>0.5 && tauID('byLooseIsolation')>0.5 && tauID('againstMuon') && leadPFChargedHadrCand.mva_e_pi < 0.6 && tauID('againstElectronCrackRem')>0.5 "),
    cut = cms.string("pt>20 && (eta<2.3&&eta>-2.3) && tauID('leadingTrackFinding')>0.5 && tauID('byLooseIsolation')>-1 && tauID('againstMuonTight')>0.5 &&  tauID('againstElectronLoose')>0.5 && tauID('againstElectronCrackRem')>0.5 "),
    filter = cms.bool(False)
    )

getattr(process,"selectedPatElectronsTriggerMatch").cut = cms.string("(eta<2.4&&eta>-2.4) && !isEBEEGap && et>15 && electronID('simpleEleId95relIso')>6.5")

process.secondLeptonVeto = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("looseMuons"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1),
    )
process.noElectronVeto = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("selectedPatElectronsTriggerMatch"),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(0),
    )
process.oneMuonLeg = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("muonLeg"),
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
process.diTau = process.muTauPairs.clone()
process.diTau.srcLeg1 = cms.InputTag("muonLeg")
process.diTau.srcLeg2 = cms.InputTag("tauLeg")
if not runOnMC:
    process.diTau.srcGenParticles = ""


process.selectedDiTau = cms.EDFilter(
    "MuTauPairSelector",
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

process.deltaRJetMuons = cms.EDProducer(
    "DeltaRNearestMuonComputer",
    probes = cms.InputTag("selectedPatJets"),
    objects = cms.InputTag("muonLeg"),
    )
process.selectedPatJetsNoMuons = cms.EDProducer(
    "JetsCleaner",
    jets =  cms.InputTag("selectedPatJets"),
    valueMap = cms.InputTag("deltaRJetMuons"),
    minDeltaR = cms.double(0.3)
    )

process.deltaRJetTaus = cms.EDProducer(
    "DeltaRNearestTauComputer",
    probes = cms.InputTag("selectedPatJetsNoMuons"),
    objects = cms.InputTag("tauLeg"),
    )
process.selectedPatJetsNoMuonsNoTaus = cms.EDProducer(
    "JetsCleaner",
    jets =  cms.InputTag("selectedPatJetsNoMuons"),
    valueMap = cms.InputTag("deltaRJetTaus"),
    minDeltaR = cms.double(0.3)
    )

process.muTauStreamAnalyzer = cms.EDAnalyzer(
    "MuTauStreamAnalyzer",
    diTaus =  cms.InputTag("selectedDiTau"),
    jets =  cms.InputTag("selectedPatJetsNoMuons"),
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
    process.PFTau*
    process.patDefaultSequence*
    process.selectedPatMuonsTriggerMatchUserEmbedded*
    process.looseMuons*
    (process.secondLeptonVeto*process.oneMuonFilter) +
    (process.noElectronVeto*process.noElecFilter) +
    (process.muonLeg * process.oneMuonLeg * process.muonLegFilter) +
    (process.tauLeg * process.oneTauLeg * process.tauLegFilter) +
    (process.diTau*process.selectedDiTau*process.atLeast1selectedDiTau*process.atLeastOneDiTauFilter) +
    (process.deltaRJetMuons*process.selectedPatJetsNoMuons)*
    (process.deltaRJetTaus*process.selectedPatJetsNoMuonsNoTaus)*
    process.muTauStreamAnalyzer+
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
                                   fileName = cms.string("treeMuTauStream.root")
                                   )

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

process.out.fileName = cms.untracked.string('patTuples_'+sample+'.root')

#process.outpath = cms.EndPath(process.out)
process.outpath = cms.EndPath()


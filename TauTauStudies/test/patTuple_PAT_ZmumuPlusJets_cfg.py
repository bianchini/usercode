from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

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
    'file:/data_CMS/cms/lbianchini/F41A3437-7AED-DF11-A50D-002618943894.root',
    #'file:/data_CMS/cms/lbianchini/ZElEl_RelVal386_1.root',
    )

#process.source.eventsToProcess = cms.untracked.VEventRange(
#    '1:1080','1:2003','1:2028','1:6867','1:7016'
#    )

postfix           = "PFlow"
sample            = ""
runOnMC           = False


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

addPFMuonIsolation(process,process.patMuons)
addPFMuon(process,postfix)
addTriggerMatchingMuon(process)
getattr(process,"patMuons").embedTrack = True
getattr(process,"patMuons"+postfix).embedTrack = True
addTriggerMatchingMuon(process,postfix)

from Bianchi.Utilities.electrons import *
addCutBasedID(process)
addPFElectronIsolation(process,process.patElectrons)
addPFElectron(process,postfix)
getattr(process,"patElectrons").embedTrack = True
getattr(process,"patElectrons"+postfix).embedTrack = True
getattr(process,"patElectrons").embedGsfTrack = True
getattr(process,"patElectrons"+postfix).embedGsfTrack = True
addTriggerMatchingElectron(process)
addTriggerMatchingElectron(process,postfix)

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

process.tightMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsTriggerMatchUserEmbedded"),
    cut = cms.string("pt>20 && (eta<2.1&&eta>-2.1) && isTrackerMuon && numberOfMatches>=2 && globalTrack.isNonnull && globalTrack.hitPattern.numberOfValidMuonHits>=1 && globalTrack.hitPattern.numberOfValidPixelHits>=1 && globalTrack.normalizedChi2<=10 && userFloat('dxyWrtPV')<0.2 && (ecalIso+hcalIso+trackIso)/pt<0.15 && (triggerObjectMatchesByPath('HLT_Mu11').size()!=0 || triggerObjectMatchesByPath('HLT_Mu9').size()!=0 || triggerObjectMatchesByPath('HLT_Mu15').size()!=0 || triggerObjectMatchesByPath('HLT_Mu15_v1').size()!=0)"),
    filter = cms.bool(False)
    )

process.looseMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsTriggerMatchUserEmbedded"),
    cut = cms.string("pt>10 && (eta<2.4&&eta>-2.4) && isGlobalMuon && globalTrack.isNonnull"),
    filter = cms.bool(False)
    )

process.diMuons = cms.EDProducer(
    "CandViewShallowCloneCombiner",
    decay = cms.string("tightMuons@- looseMuons@+"),
    cut   = cms.string("60 < mass < 120"),
    )


process.atLeast1diMuon = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("diMuons"),
    minNumber = cms.uint32(1),
    )

getattr(process,"selectedPatElectronsTriggerMatch"+postfix).cut = cms.string('(eta<2.4&&eta>-2.4 && !isEBEEGap) && et>10')
getattr(process,"selectedPatElectronsTriggerMatch").cut = cms.string('(eta<2.4&&eta>-2.4 && !isEBEEGap) && et>10')

getattr(process,"selectedPatJets").cut = cms.string('pt>5')

process.deltaRJetMuons = cms.EDProducer(
    "DeltaRNearestMuonComputer",
    probes = cms.InputTag("selectedPatJets"),
    objects = cms.InputTag("looseMuons"),
    )

process.selectedPatJetsNoMuons = cms.EDProducer(
    "JetsCleaner",
    jets =  cms.InputTag("selectedPatJets"),
    valueMap = cms.InputTag("deltaRJetMuons"),
    minDeltaR = cms.double(0.3)
    )

process.zPlusJetsAnalyzer = cms.EDAnalyzer(
    "ZmumuPlusJetsAnalyzer",
    diMuons =  cms.InputTag("diMuons"),
    jets =  cms.InputTag("selectedPatJetsNoMuons"),
    isMC = cms.bool(runOnMC),
    minCorrPt = cms.untracked.double(5.),
    minJetID  = cms.untracked.double(0.5), # 1=loose,2=medium,3=tight
    verbose =  cms.untracked.bool(False),
    )

process.pat = cms.Sequence(
    process.allEventsFilter+
    process.primaryVertexFilter+
    process.scrapping +
    process.makeSCs +
    process.fjSequence*
    process.patDefaultSequence*
    process.selectedPatMuonsTriggerMatchUserEmbedded*
    (process.looseMuons+ process.tightMuons)*
    process.diMuons*
    process.atLeast1diMuon*
    process.deltaRJetMuons*
    process.selectedPatJetsNoMuons*
    process.zPlusJetsAnalyzer+
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
    'keep *_selectedPatJetsNoMuons_*_*',
    'keep *_diMuons_*_*',
    )
                                   )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("treeZmumuPlusJets.root")
                                   )

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

process.out.fileName = cms.untracked.string('patTuples_'+sample+'.root')

#process.outpath = cms.EndPath(process.out)
process.outpath = cms.EndPath()


from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START38_V13::All'

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')


process.source.fileNames = cms.untracked.vstring(
    #'file:/data_CMS/cms/lbianchini/DY_ZTauTau_Fall10_RECO.root'
    #'file:/data_CMS/cms/lbianchini/ZTT_RelVal386_1.root',
    #'file:/data_CMS/cms/lbianchini/ZTT_RelVal386_2.root',
    #'file:/data_CMS/cms/lbianchini/ZTT_RelVal386_3.root',
    #'file:/data_CMS/cms/lbianchini/ZTT_RelVal386_4.root',
    #'file:/data_CMS/cms/lbianchini/F41A3437-7AED-DF11-A50D-002618943894.root',
    'file:/data_CMS/cms/lbianchini/ZElEl_RelVal386_1.root',
    #'file:/data_CMS/cms/lbianchini/ZElEl_RelVal386_2.root',
    #'file:/data_CMS/cms/lbianchini/ZElEl_RelVal386_3.root',
    #'file:/data_CMS/cms/lbianchini/ZElEl_RelVal386_4.root',
    )

#process.source.eventsToProcess = cms.untracked.VEventRange(
#    '1:1080','1:2003','1:2028','1:6867','1:7016'
#    )

postfix           = "PFlow"
sample            = ""
runOnMC           = False

### TnP
##########################
makeEtoTauFakeRate= True
###########       ########
makeMCtrees       = False
makeUnbiased      = False
removeWenuFilters = False
##########################

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree1 = cms.EDAnalyzer("ParticleListDrawer",
                                    src = cms.InputTag("genParticles"),
                                    maxEventsToPrint  = cms.untracked.int32(1)
                                    )

process.load("Bianchi.eToTaufakeRate.filters_ZmumuPlusJets_cff")

from PhysicsTools.PatAlgos.tools.coreTools import *

if not runOnMC:
    removeMCMatching(process,["All"])
    
removeSpecificPATObjects(process, ['Photons'],
                         outputInProcess=False)
removeCleaning(process,
               outputInProcess=False)

restrictInputToAOD(process, ['All'])

process.load("PhysicsTools.PFCandProducer.ParticleSelectors.pfSortByType_cff")
process.load("PhysicsTools.PFCandProducer.pfNoPileUp_cff")

process.pfCandidateSelectionByType = cms.Sequence(
    process.pfNoPileUpSequence *
    ( process.pfAllNeutralHadrons +
      process.pfAllChargedHadrons +
      process.pfAllPhotons
      )  +
    process.pfAllMuons +
    process.pfAllElectrons
    )
process.pfPileUp.Enable    = True # enable pile-up filtering
process.pfPileUp.Vertices  = "offlinePrimaryVertices" # use vertices w/o BS
process.pfAllMuons.src     = "particleFlow"
process.pfAllElectrons.src = "particleFlow"

process.patDefaultSequence.replace(process.patCandidates,
                                   process.pfCandidateSelectionByType+
                                   process.patCandidates)

process.makePatTaus.remove(process.patPFCandidateIsoDepositSelection)
from Bianchi.eToTaufakeRate.customizePAT import *
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

from PhysicsTools.PatAlgos.tools.tauTools import *

switchToPFTauHPS(process, 
                 pfTauLabelOld = 'shrinkingConePFTauProducer',
                 pfTauLabelNew = 'hpsPFTauProducer'
                 )

cloneProcessingSnippet(process, process.makePatTaus, postfix)


switchToPFTauShrinkingCone(process,
                           pfTauLabelOld = 'hpsPFTauProducer',
                           pfTauLabelNew = 'shrinkingConePFTauProducer'
                           )

process.patCandidates.replace(process.makePatTaus,
                              process.makePatTaus+
                              getattr(process,"makePatTaus"+postfix)
                              )
process.patCandidateSummary.candidates.append(cms.InputTag("patTaus"+postfix))

setattr(process,"selectedPatTaus"+postfix,process.selectedPatTaus.clone())
getattr(process,"selectedPatTaus"+postfix).src = 'patTaus'+postfix
process.selectedPatCandidates.replace(process.selectedPatTaus,
                                      process.selectedPatTaus+
                                      getattr(process,"selectedPatTaus"+postfix)
                                      )
process.selectedPatCandidateSummary.candidates.append(cms.InputTag("selectedPatTaus"+postfix))
# and add counter
setattr(process,"countPatTaus"+postfix,process.countPatTaus.clone())
getattr(process,"countPatTaus"+postfix).src = 'selectedPatTaus'+postfix
process.countPatCandidates.replace(process.countPatTaus,
                                   process.countPatTaus+
                                   getattr(process,"countPatTaus"+postfix)
                                   )

getattr(process,"patTaus"+postfix).embedIsolationTracks = cms.bool(True)
getattr(process,"patTaus"+postfix).embedSignalTracks = cms.bool(True)
getattr(process,"patTaus"+postfix).embedGenMatch = cms.bool(True)
getattr(process,"patTaus"+postfix).embedLeadTrack = cms.bool(True)
getattr(process,"patTaus"+postfix).embedLeadPFCand = True
getattr(process,"patTaus"+postfix).embedLeadPFChargedHadrCand = True
getattr(process,"patTaus"+postfix).embedLeadPFNeutralCand = True
getattr(process,"patTaus"+postfix).embedSignalPFCands = True
getattr(process,"patTaus"+postfix).embedSignalPFChargedHadrCands = True
getattr(process,"patTaus"+postfix).embedSignalPFNeutralHadrCands = True
getattr(process,"patTaus"+postfix).embedSignalPFGammaCands = True
getattr(process,"patTaus"+postfix).embedIsolationPFCands = True
getattr(process,"patTaus"+postfix).embedIsolationPFChargedHadrCands = True
getattr(process,"patTaus"+postfix).embedIsolationPFNeutralHadrCands = True
getattr(process,"patTaus"+postfix).embedIsolationPFGammaCands = True
getattr(process,"patTaus"+postfix).embedGenJetMatch = cms.bool(True)
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

setattr(process,"hpsPFTauDiscriminationAgainstElectron2D",
        getattr(process,"hpsPFTauDiscriminationAgainstElectron").clone(
    ApplyCut_ElectronPreID_2D = cms.bool(True),
    ApplyCut_PFElectronMVA =  cms.bool(False)
    )
            )
setattr(process,"hpsPFTauDiscriminationAgainstElectronCrackRem",
        getattr(process,"hpsPFTauDiscriminationAgainstElectron").clone(
    ApplyCut_EcalCrackCut = cms.bool(True),
    ApplyCut_PFElectronMVA =  cms.bool(False)
    )
        )
    
setattr(process,"shrinkingConePFTauDiscriminationAgainstElectron2D",
        getattr(process,"shrinkingConePFTauDiscriminationAgainstElectron").clone(
    ApplyCut_ElectronPreID_2D = cms.bool(True),
    ApplyCut_PFElectronMVA =  cms.bool(False)
    )
        )
setattr(process,"shrinkingConePFTauDiscriminationAgainstElectronCrackRem",
        getattr(process,"shrinkingConePFTauDiscriminationAgainstElectron").clone(
    ApplyCut_EcalCrackCut = cms.bool(True),
    ApplyCut_PFElectronMVA =  cms.bool(False)
    )
        )

process.patHPSPFTauDiscrimination += process.hpsPFTauDiscriminationAgainstElectron2D
process.patHPSPFTauDiscrimination += process.hpsPFTauDiscriminationAgainstElectronCrackRem
process.patShrinkingConePFTauDiscrimination += process.shrinkingConePFTauDiscriminationAgainstElectron2D
process.patShrinkingConePFTauDiscrimination += process.shrinkingConePFTauDiscriminationAgainstElectronCrackRem

getattr(process,"makePatTaus"+postfix).replace(
    getattr(process,"patTaus"+postfix),
    process.patHPSPFTauDiscrimination + getattr(process,"patTaus"+postfix)
    )
getattr(process,"makePatTaus").replace(
    getattr(process,"patTaus"),
    process.patShrinkingConePFTauDiscrimination + getattr(process,"patTaus")
    )


getattr(process,"patTaus"+postfix).tauIDSources = cms.PSet(
    leadingTrackFinding = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
    byLooseIsolation = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"),
    byMediumIsolation = cms.InputTag("hpsPFTauDiscriminationByMediumIsolation"),
    byTightIsolation = cms.InputTag("hpsPFTauDiscriminationByTightIsolation"),
    againstElectron = cms.InputTag("hpsPFTauDiscriminationAgainstElectron"),
    againstElectron2D = cms.InputTag("hpsPFTauDiscriminationAgainstElectron2D"),
    againstElectronCrackRem = cms.InputTag("hpsPFTauDiscriminationAgainstElectronCrackRem"),
    againstMuon = cms.InputTag("hpsPFTauDiscriminationAgainstMuon")
    )

getattr(process,"patTaus").tauIDSources = cms.PSet(
    leadingTrackFinding = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding"),
    leadingTrackPtCut = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackPtCut"),
    leadingPionPtCut = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingPionPtCut"),
    trackIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolation"),
    trackIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion"),
    ecalIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolation"),
    ecalIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion"),
    byIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByIsolation"),
    byIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByIsolationUsingLeadingPion"),
    againstElectron = cms.InputTag("shrinkingConePFTauDiscriminationAgainstElectron"),
    againstElectron2D = cms.InputTag("shrinkingConePFTauDiscriminationAgainstElectron2D"),
    againstElectronCrackRem = cms.InputTag("shrinkingConePFTauDiscriminationAgainstElectronCrackRem"),
    againstMuon = cms.InputTag("shrinkingConePFTauDiscriminationAgainstMuon"),
    byTaNC = cms.InputTag("shrinkingConePFTauDiscriminationByTaNC"),
    byTaNCfrOnePercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrOnePercent"),
    byTaNCfrHalfPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrHalfPercent"),
    byTaNCfrQuarterPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrQuarterPercent"),
    byTaNCfrTenthPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrTenthPercent")
    )

process.patMuons.embedTrack = True
addPFMuonIsolation(process,process.patMuons)
addPFMuon(process,postfix)
addTriggerMatchingMuon(process)
getattr(process,"muonTriggerMatchHLTMuons").pathNames=cms.vstring('*')
addTriggerMatchingMuon(process,postfix)
getattr(process,"muonTriggerMatchHLTMuons"+postfix).pathNames=cms.vstring('*')

process.patElectrons.embedGsfTrack = True
from Bianchi.eToTaufakeRate.electrons import *
addCutBasedID(process)
addPFElectronIsolation(process,process.patElectrons)
addPFElectron(process,postfix)
process.patElectronsPFlow.embedGsfTrack = True
addTriggerMatchingElectron(process)
addTriggerMatchingElectron(process,postfix)

if hasattr(process,"patTrigger"):
    process.patTrigger.processName = '*'

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


process.genTauDecaysToHadrons = cms.EDFilter("TauGenJetDecayModeSelector",
                                             src = cms.InputTag("tauGenJets"),
                                             select = cms.vstring('oneProng0Pi0', 'oneProng1Pi0', 'oneProng2Pi0', 'oneProngOther','threeProng0Pi0', 'threeProng1Pi0', 'threeProngOther', 'rare'),
                                             filter = cms.bool(False)
                                             )
process.genTauDecaysToElectron = cms.EDFilter("TauGenJetDecayModeSelector",
                                              src = cms.InputTag("tauGenJets"),
                                              select = cms.vstring('electron'),
                                              filter = cms.bool(False)
                                              )

process.tauFakeRateAnalyzerHPS = cms.EDAnalyzer("TauFakeRateAnalyzer",
                                                isMC = cms.bool(runOnMC),
                                                matchTo = cms.string("tau"),
                                                tauTag = cms.InputTag("selectedPatTausPFlow"),
                                                electronsTag = cms.InputTag("patElectrons"),
                                                )
process.tauFakeRateAnalyzerSHC = cms.EDAnalyzer("TauFakeRateAnalyzer",
                                                isMC = cms.bool(runOnMC),
                                                matchTo = cms.string("tau"),
                                                tauTag = cms.InputTag("selectedPatTaus"),
                                                electronsTag = cms.InputTag("patElectrons"),
                                                )
process.makeEtoTauEff = cms.Sequence(process.genTauDecaysToHadrons*
                                     (process.tauFakeRateAnalyzerHPS+process.tauFakeRateAnalyzerSHC))
    
#########
## PAT

getattr(process,"selectedPatJets").cut = cms.string('pt>10')
getattr(process,"selectedPatTaus").cut = cms.string('pt>15 && (eta<2.4&&eta>-2.4) && tauID("leadingTrackFinding") > 0.5 && tauID("leadingPionPtCut") > 0.5 && tauID("byIsolation")>0.5 && tauID("againstElectronCrackRem") > 0.5')
getattr(process,"selectedPatTaus"+postfix).cut = cms.string('pt>15 && (eta<2.4&&eta>-2.4) && tauID("leadingTrackFinding") > 0.5 && tauID("byLooseIsolation")>0.5 && tauID("againstElectronCrackRem") > 0.5')
getattr(process,"selectedPatMuonsTriggerMatch"+postfix).cut = cms.string('pt>10 && (eta<2.1&&eta>-2.1) && isTrackerMuon && numberOfMatches>=2 && globalTrack.isNonnull  && globalTrack.hitPattern.numberOfValidMuonHits>=1 && globalTrack.hitPattern.numberOfValidPixelHits>=1 && globalTrack.normalizedChi2<=10')
getattr(process,"selectedPatMuonsTriggerMatch").cut = cms.string('pt>10 && (eta<2.1&&eta>-2.1) && isTrackerMuon && numberOfMatches>=2 && globalTrack.isNonnull && globalTrack.hitPattern.numberOfValidMuonHits>=1 && globalTrack.hitPattern.numberOfValidPixelHits>=1 && globalTrack.normalizedChi2<=10')

getattr(process,"patElectrons").embedTrack = cms.bool(True)
getattr(process,"patElectrons"+postfix).embedTrack = cms.bool(True)

getattr(process,"selectedPatElectronsTriggerMatch"+postfix).cut = cms.string('(eta<2.4&&eta>-2.4 && !isEBEEGap) && et>10')
getattr(process,"selectedPatElectronsTriggerMatch").cut = cms.string('(eta<2.4&&eta>-2.4 && !isEBEEGap) && et>10')


process.pat = cms.Sequence(
    process.allEventsFilter+
    process.primaryVertexFilter+
    process.scrapping +
    process.makeSCs +
    process.patDefaultSequence
    #*process.makeEtoTauEff+
    #process.printTree1
    )

if not runOnMC:
    process.pat.remove(process.printTree1)

from Bianchi.eToTaufakeRate.tnpEtoTau.addTnPSequences import addTnPSequences
if makeEtoTauFakeRate:
    addTnPSequences(process,"pat",makeMCtrees,makeUnbiased,removeWenuFilters)
else:
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
    'keep TrackingRecHitsOwned_*_*_*',
    'keep *_selectedSuperClusters_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_offlinePrimaryVertices*_*_*',
    'keep *_particleFlow_*_*',
    'keep *_pfAllMuons_*_*',
    'keep *_pfAllElectrons_*_*',
    'keep *_tauGenJets_*_*',
    'keep *_genTauDecaysToElectron_*_*',
    'keep *_genTauDecaysToMuon_*_*',
    'keep *_genTauDecaysToHadrons_*_*',
    'keep *_selectedSuperClusters_*_*',
    'keep *_WenuPair90_*_*'
    )
                                   )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("testNewWriteFromPAT.root")
                                   )

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

process.out.fileName = cms.untracked.string('patTuples_'+sample+'.root')

if makeEtoTauFakeRate:
    process.outpath = cms.EndPath()
else:
    #process.outpath = cms.EndPath(process.out)
    process.outpath = cms.EndPath()


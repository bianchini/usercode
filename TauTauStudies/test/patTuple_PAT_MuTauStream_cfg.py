from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

## temporary JEC
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.jec = cms.ESSource(
    "PoolDBESSource",
    DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0)
    ),
    timetype = cms.string('runnumber'),
    toGet = cms.VPSet(
    cms.PSet(
    record = cms.string('JetCorrectionsRecord'),
    tag    = cms.string('JetCorrectorParametersCollection_Jec10V3_AK5PF'),
    label  = cms.untracked.string('AK5PF')
    )
    ),
    ## here you add as many jet types as you need (AK5Calo, AK5JPT, AK7PF, AK7Calo, KT4PF, KT4Calo, KT6PF, KT6Calo)
    connect = cms.string('sqlite_file:Jec10V3.db')
    )

process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')


process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.Rho_EtaMax = cms.double(4.4)
process.kt6PFJets.Ghost_EtaMax = cms.double(5.0)
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.Rho_EtaMax = cms.double(4.4)
process.ak5PFJets.Ghost_EtaMax = cms.double(5.0)

process.ak5PFL1Fastjet.useCondDB = False

## re-run kt4PFJets within lepton acceptance to compute rho
process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJetsCentral = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsCentral.Rho_EtaMax = cms.double(2.5)

process.fjSequence = cms.Sequence(process.kt6PFJets+process.ak5PFJets+process.kt6PFJetsCentral)

process.source.fileNames = cms.untracked.vstring(
    #'file:/data_CMS/cms/lbianchini/ZTT_RelVal386_1.root',
    #'file:/data_CMS/cms/lbianchini/ZMuMu_RelVal386.root',
    # 'rfio:/dpm/in2p3.fr/home/cms/trivcat//store/mc/Spring11/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/AODSIM/PU_S1_START311_V1G1-v2/0000/FA5943AB-A756-E011-A6C8-002618FDA208.root',
    # 'rfio:/dpm/in2p3.fr/home/cms/trivcat//store/mc/Spring11/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/AODSIM/PU_S1_START311_V1G1-v2/0000/02700D19-A756-E011-B4D5-0030486792B6.root'
    #'file:goodDataEvents_84_1_yHS.root'
    'file:/data_CMS/cms/akalinow/VBF_HToTauTau_M-115_7TeV-powheg-pythia6-tauola/PU_S1_START311_V1G1-v1/AOD/8EC598C7-3453-E011-AC82-002481E14F8C.root'
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat//store/mc/Fall10/VBF_HToTauTau_M-115_7TeV-powheg-pythia6-tauola/GEN-SIM-RECO/START38_V12-v1/0000/044E940A-55EC-DF11-89D6-0023AEFDEE60.root',
    #'file:/data_CMS/cms/lbianchini/F41A3437-7AED-DF11-A50D-002618943894.root',
    #'file:/data_CMS/cms/lbianchini/ZElEl_RelVal386_1.root',
    )

#process.source.eventsToProcess = cms.untracked.VEventRange(
#    '1:751063','1:751088','1:751100','1:751117','1:751141','1:751237','1:751325','1:751392','1:751466','1:977809','1:977870','1:977898',
#    '1:977915','1:978085','1:978214','1:811102','1:811190'
#    )

postfix           = "PFlow"
sample            = ""
runOnMC           = True

if runOnMC:
    process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
else:
    #process.GlobalTag.globaltag = cms.string(autoCond[ 'com10' ])
     process.GlobalTag.globaltag = cms.string('GR_R_311_V4::All')

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
    decayModeFinding = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
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


process.tauMatch.maxDeltaR = 0.15
process.tauMatch.resolveAmbiguities = cms.bool(False)
process.tauGenJetMatch.resolveAmbiguities = cms.bool(False)
process.tauGenJetMatch.maxDeltaR = 0.15
process.tauGenJetMatch.maxDPtRel = 999
## <\tau part>

addPFMuonIsolation(process,process.patMuons)
addTriggerMatchingMuon(process,isMC=runOnMC)
getattr(process,"patMuons").embedTrack = True

from Bianchi.Utilities.electrons import *
addCutBasedID(process)
addPFElectronIsolation(process,process.patElectrons)

getattr(process,"patElectrons").embedTrack = True
getattr(process,"patElectrons").embedGsfTrack = True
addTriggerMatchingElectron(process,isMC=runOnMC)

addTriggerMatchingTau(process,isMC=runOnMC,postfix="",XtriggerMu=True)

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )
process.patTriggerEvent.processName = '*'

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

    
########### PAT

process.selectedPatMuonsTriggerMatchUserEmbedded = cms.EDProducer(
    "MuonsUserEmbedded",
    muonTag = cms.InputTag("selectedPatMuonsTriggerMatch"),
    vertexTag = cms.InputTag("offlinePrimaryVertices")
    )

process.atLeastOneMuTau = cms.EDProducer("CandViewShallowCloneCombiner",
                                         decay = cms.string("selectedPatMuonsTriggerMatch selectedPatTausTriggerMatch"),
                                         cut   = cms.string("sqrt((daughter(0).eta-daughter(1).eta)*(daughter(0).eta-daughter(1).eta)+(daughter(0).phi-daughter(1).phi)*(daughter(0).phi-daughter(1).phi))>0.3"),
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
    src = cms.InputTag("selectedPatMuonsTriggerMatchUserEmbedded"),
    cut = cms.string("pt>15 && abs(eta)<2.1"),
    filter = cms.bool(False)
    )
process.atLeastOneMuTaumuPtEta = process.atLeastOneMuTau.clone(
    decay=cms.string("muPtEta selectedPatTausTriggerMatch")
    )
process.muPtEtaCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneMuTaumuPtEta"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

process.muPtEtaRelID = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsTriggerMatchUserEmbedded"),
    cut = cms.string("pt>15 && abs(eta)<2.4 &&  isTrackerMuon && numberOfMatches>=2 && globalTrack.isNonnull && globalTrack.hitPattern.numberOfValidMuonHits>=1 && globalTrack.hitPattern.numberOfValidPixelHits>=1 && globalTrack.hitPattern.numberOfValidTrackerHits>10 && globalTrack.normalizedChi2<10 && globalTrack.ptError/globalTrack.pt<0.1 && abs(userFloat('dxyWrtPV'))<0.045 && abs(userFloat('dzWrtPV'))<0.2"),
    filter = cms.bool(False)
    )

process.muPtEtaID = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsTriggerMatchUserEmbedded"),
    cut = cms.string(process.muPtEta.cut.value()+" &&  isTrackerMuon && numberOfMatches>=2 && globalTrack.isNonnull && globalTrack.hitPattern.numberOfValidMuonHits>=1 && globalTrack.hitPattern.numberOfValidPixelHits>=1 && globalTrack.hitPattern.numberOfValidTrackerHits>10 && globalTrack.normalizedChi2<10 && globalTrack.ptError/globalTrack.pt<0.1 && abs(userFloat('dxyWrtPV'))<0.045 && abs(userFloat('dzWrtPV'))<0.2"),
    filter = cms.bool(False)
    )
process.atLeastOneMuTaumuPtEtaID = process.atLeastOneMuTau.clone(
    decay=cms.string("muPtEtaID selectedPatTausTriggerMatch")
    )
process.muPtEtaIDCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneMuTaumuPtEtaID"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

process.tauPtEta  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("selectedPatTausTriggerMatch"),
    cut = cms.string("pt>20 && abs(eta)<2.3"),
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
    src = cms.InputTag("selectedPatTausTriggerMatch"),
    cut = cms.string(process.tauPtEta.cut.value()+" && tauID('decayModeFinding')>0.5"),
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
    src = cms.InputTag("selectedPatTausTriggerMatch"),
    cut = cms.string(process.tauPtEtaID.cut.value()+" && tauID('againstMuonTight')>0.5"),
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
    src = cms.InputTag("selectedPatTausTriggerMatch"),
    cut = cms.string(process.tauPtEtaIDAgMu.cut.value()+" && tauID('againstElectronLoose')>0.5"),
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

process.alLeastOneMuTauSequence = cms.Sequence(
    process.atLeastOneMuTau*process.atLeastOneMuTauCounter*process.atLeastOneMuTauFilter
    )

process.muLegSequence = cms.Sequence(
    (process.muPtEta*process.atLeastOneMuTaumuPtEta*process.muPtEtaCounter*process.muPtEtaFilter) *
    (process.muPtEtaID*process.atLeastOneMuTaumuPtEtaID*process.muPtEtaIDCounter*process.muPtEtaIDFilter) *
    process.muPtEtaRelID
    )
process.tauLegSequence = cms.Sequence(
    (process.tauPtEta*process.atLeastOneMuTautauPtEta*process.tauPtEtaCounter*process.tauPtEtaFilter) *
    (process.tauPtEtaID*process.atLeastOneMuTautauPtEtaID*process.tauPtEtaIDCounter*process.tauPtEtaIDFilter) *
    (process.tauPtEtaIDAgMu*process.atLeastOneMuTautauPtEtaIDAgMu*process.tauPtEtaIDAgMuCounter*process.tauPtEtaIDAgMuFilter)*
    (process.tauPtEtaIDAgMuAgElec*process.atLeastOneMuTautauPtEtaIDAgMuAgElec*process.tauPtEtaIDAgMuAgElecCounter*process.tauPtEtaIDAgMuAgElecFilter)
    )

process.load("Bianchi.Utilities.diTausReconstruction_cff")
process.diTau = process.muTauPairs.clone()
process.diTau.srcLeg1 = cms.InputTag("muPtEtaID")
process.diTau.srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElec")
process.diTau.srcMET  = cms.InputTag("patMETsPFlow")

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

process.diTauSequence = cms.Sequence(
    process.diTau*
    process.selectedDiTau*
    process.atLeast1selectedDiTau*
    process.atLeast1selectedDiTauFilter
    )


getattr(process,"selectedPatJets").cut = cms.string('et>10 && abs(eta)<5.0')

process.deltaRJetMuons = cms.EDProducer(
    "DeltaRNearestMuonComputer",
    probes = cms.InputTag("selectedPatJets"),
    #objects = cms.InputTag("muonLeg"),
    objects = cms.InputTag("muPtEtaID"),
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
    #objects = cms.InputTag("tauLeg"),
    objects = cms.InputTag("tauPtEtaIDAgMuAgElec"),
    )
process.selectedPatJetsNoMuonsNoTaus = cms.EDProducer(
    "JetsCleaner",
    jets =  cms.InputTag("selectedPatJetsNoMuons"),
    valueMap = cms.InputTag("deltaRJetTaus"),
    minDeltaR = cms.double(0.3)
    )

process.jetCleaningSequence = cms.Sequence(
    process.deltaRJetMuons*process.selectedPatJetsNoMuons*process.deltaRJetTaus*process.selectedPatJetsNoMuonsNoTaus
    )

process.muTauStreamAnalyzer = cms.EDAnalyzer(
    "MuTauStreamAnalyzer",
    diTaus =  cms.InputTag("selectedDiTau"),
    #jets =  cms.InputTag("selectedPatJetsNoMuons"),
    jets =  cms.InputTag("selectedPatJets"),
    triggerResults = cms.InputTag("patTriggerEvent"),
    isMC = cms.bool(runOnMC),
    deltaRLegJet  = cms.untracked.double(0.3),
    minCorrPt = cms.untracked.double(10.), # is Et !!
    minJetID  = cms.untracked.double(0.5), # 1=loose,2=medium,3=tight
    applyTauSignalSel =  cms.bool( True ),
    verbose =  cms.untracked.bool( False ),
    )

#process.load("Bianchi.TauTauStudies.preselectionSeq_cff")

process.pat = cms.Sequence(
    process.allEventsFilter+
    (process.primaryVertexFilter+process.scrapping)*
    process.vertexScrapingFilter +
    process.makeSCs +
    process.fjSequence*
    process.PFTau*
    process.patDefaultSequence*
    process.selectedPatMuonsTriggerMatchUserEmbedded*
    #process.arturPreselection*
    process.alLeastOneMuTauSequence*
    process.muLegSequence*
    process.tauLegSequence*
    process.diTauSequence*
    process.jetCleaningSequence*
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


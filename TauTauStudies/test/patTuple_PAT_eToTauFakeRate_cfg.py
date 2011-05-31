from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

process.source.fileNames = cms.untracked.vstring(
    #'file:/data_CMS/cms/lbianchini/ZTT_RelVal386_1.root',
    #'file:/data_CMS/cms/lbianchini/ZMuMu_RelVal386.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat//store/mc/Spring11/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/AODSIM/PU_S1_START311_V1G1-v2/0000/FA5943AB-A756-E011-A6C8-002618FDA208.root',
    #'file:goodDataEvents_84_1_yHS.root'
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat//store/mc/Spring11/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/AODSIM/PU_S1_START311_V1G1-v2/0000/FA5943AB-A756-E011-A6C8-002618FDA208.root'
    'rfio:/dpm/in2p3.fr/home/cms/trivcat//store/mc/Spring11/VBF_ToHToZZTo4L_M-550_7TeV-powheg-pythia6/AODSIM/PU_S1_START311_V1G1-v1/0034/126992AB-6554-E011-BD51-003048D476B8.root'
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat//store/mc/Spring11/DYToEE_M-20_TuneZ2_7TeV-pythia6/GEN-SIM-RECODEBUG/E7TeV_FlatDist10_2011EarlyData_50ns_START311_V1G1-v1/0000/0053C2AC-423C-E011-976F-00215E21DB3A.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat//store/mc/Spring11/DYToEE_M-20_TuneZ2_7TeV-pythia6/GEN-SIM-RECODEBUG/E7TeV_FlatDist10_2011EarlyData_50ns_START311_V1G1-v1/0000/00DCEC58-433B-E011-8FF8-E41F13181A70.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat//store/mc/Spring11/DYToEE_M-20_TuneZ2_7TeV-pythia6/GEN-SIM-RECODEBUG/E7TeV_FlatDist10_2011EarlyData_50ns_START311_V1G1-v1/0000/02300FF5-423B-E011-B949-00215E2221E4.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat//store/mc/Spring11/DYToEE_M-20_TuneZ2_7TeV-pythia6/GEN-SIM-RECODEBUG/E7TeV_FlatDist10_2011EarlyData_50ns_START311_V1G1-v1/0000/023A43AC-433B-E011-A582-00215E21DD26.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat//store/mc/Fall10/VBF_HToTauTau_M-115_7TeV-powheg-pythia6-tauola/GEN-SIM-RECO/START38_V12-v1/0000/044E940A-55EC-DF11-89D6-0023AEFDEE60.root',
    #'file:/data_CMS/cms/lbianchini/F41A3437-7AED-DF11-A50D-002618943894.root',
    #'file:/data_CMS/cms/lbianchini/ZElEl_RelVal386_1.root',
    )

#process.source.eventsToProcess = cms.untracked.VEventRange(
#    '1:1080','1:2003','1:2028','1:6867','1:7016'
#    )

postfix           = "PFlow"
runOnMC           = True
doEfficiency      = False

if not doEfficiency:
    FileName = "treeEtoTauTnP.root"
else:
    FileName = "treeEtoTauEff.root"

if runOnMC:
    process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
else:
    #process.GlobalTag.globaltag = cms.string(autoCond[ 'com10' ])
    process.GlobalTag.globaltag = cms.string('GR_R_41_V0::All')

process.primaryVertexFilter = cms.EDFilter(
    "GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVerticesDA'),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
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

addPFElectronIsolation(process,process.patElectrons)

process.pfPileUp.Vertices = "offlinePrimaryVerticesDA"

getattr(process,"patElectrons").embedTrack = True
getattr(process,"patElectrons").embedGsfTrack = True
addTriggerMatchingElectron(process,isMC=runOnMC)
process.eleTriggerMatchHLTElectrons.matchedCuts =  cms.string('(path("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*") || path("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*") ) && type("TriggerElectron")')


if hasattr(process,"patTrigger"):
    process.patTrigger.processName = '*'


process.selectedPatElectronsTriggerMatchUserEmbedded = cms.EDProducer(
    "ElectronsUserEmbedded",
    electronTag = cms.InputTag("selectedPatElectronsTriggerMatch"),
    vertexTag = cms.InputTag("offlinePrimaryVerticesDA"),
    isMC = cms.bool(runOnMC)
    )

#####################################################################################
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
process.genTauDecaysToMuon = cms.EDFilter("TauGenJetDecayModeSelector",
                                          src = cms.InputTag("tauGenJets"),
                                          select = cms.vstring('muon'),
                                          filter = cms.bool(False)
                                          )
process.taus = cms.EDFilter("PATTauSelector",
                            src = cms.InputTag("selectedPatTaus"),
                            cut = cms.string("pt>15.0 && abs(eta)<2.5 && tauID('decayModeFinding')>0.5"),
                            filter = cms.bool(False)
                            )
process.tauFakeRateAnalyzer = cms.EDAnalyzer("TauFakeRateAnalyzer",
                                             matchTo = cms.string("muon"),
                                             tauTag = cms.InputTag("taus"),
                                             electronsTag = cms.InputTag("selectedPatElectronsTriggerMatchUserEmbedded"),
                                             )

process.makeEtoTauEff = cms.Sequence(process.genTauDecaysToHadrons+
                                     process.taus*
                                     process.tauFakeRateAnalyzer)
#####################################################################################
if not doEfficiency:
    process.load("Bianchi.TauTauStudies.eToTauFakeRate_cff")

if (not doEfficiency) and (not runOnMC):
    process.sequence70.remove(process.tag70McMatch)
    process.sequence70.remove(process.probeIDLooseMcMatch)
    process.sequence70.remove(process.probeIDMediumMcMatch)
    process.sequence70.remove(process.probeIDTightMcMatch)
    process.etoTau70IDLoose.isMC = cms.bool(False)
    process.etoTau70IDMedium.isMC = cms.bool(False)
    process.etoTau70IDTight.isMC = cms.bool(False)
    process.etoTau70IDLoose.makeMCUnbiasTree = cms.bool(False)
    process.etoTau70IDMedium.makeMCUnbiasTree= cms.bool(False)
    process.etoTau70IDTight.makeMCUnbiasTree = cms.bool(False)
    process.addUserVariables.isMC = cms.bool(False)

process.pat = cms.Sequence(
    process.allEventsFilter+
    process.PFTau*
    process.primaryVertexFilter*
    process.patDefaultSequence*
    process.selectedPatElectronsTriggerMatchUserEmbedded
    )

if not doEfficiency:
    process.pat.replace(process.selectedPatElectronsTriggerMatchUserEmbedded,
                        process.selectedPatElectronsTriggerMatchUserEmbedded*process.sequence70)
else:
    process.pat.replace(process.selectedPatElectronsTriggerMatchUserEmbedded,
                        process.selectedPatElectronsTriggerMatchUserEmbedded*process.makeEtoTauEff)

process.p = cms.Path(process.pat)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(FileName)
                                   )

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

process.out.fileName = cms.untracked.string('patTuples_eToTauFakeRate.root')
process.outpath = cms.EndPath()

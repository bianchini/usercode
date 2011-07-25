import FWCore.ParameterSet.Config as cms

process = cms.Process("MUTAUANA")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

runOnMC = True

if runOnMC:
    process.GlobalTag.globaltag = cms.string('START42_V12::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V14::All')
    
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 2000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    'file:./patTuples_MuTauStream.root'
    )
    )

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    )

process.load("Bianchi.Utilities.diTausReconstruction_cff")
process.diTau = process.allMuTauPairs.clone()
process.diTau.srcLeg1 = cms.InputTag("muPtEtaID")
process.diTau.srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElec")
process.diTau.srcMET  = cms.InputTag("patMETsPFlow")
process.diTau.dRmin12  = cms.double(0.5)
process.diTau.doSVreco = cms.bool(False)

if not runOnMC:
    process.diTau.srcGenParticles = ""
        
process.selectedDiTau = cms.EDFilter(
    "MuTauPairSelector",
    src = cms.InputTag("diTau"),
    cut = cms.string("dR12>0.5")
    )

#######################################################################3

process.diTau2 = process.diTau.clone(doSVreco = cms.bool(True))
process.selectedDiTau2 = process.selectedDiTau.clone(src = cms.InputTag("diTau2") )
process.tauPtEtaIDAgMuAgElecIso  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("tauPtEtaIDAgMuAgElec"),
    cut = cms.string("tauID('byLooseCombinedIsolationDeltaBetaCorr')>0.5"),
    filter = cms.bool(False)
    )
process.tauPtEtaIDAgMuAgElecIsoCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )
process.muPtEtaIDIso  = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("muPtEtaID"),
    cut = cms.string("userFloat('PFRelIsoDB04v2')<0.2"),
    filter = cms.bool(False)
    )
process.muPtEtaIDIsoCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("muPtEtaIDIso"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )
process.selectedPatJets25  = cms.EDFilter(
    "PATJetSelector",
    src = cms.InputTag("selectedPatJets"),
    cut = cms.string("pt>25"),
    filter = cms.bool(True)
    )
process.deltaRJetMuons = cms.EDProducer(
    "DeltaRNearestMuonComputer",
    probes = cms.InputTag("selectedPatJets25"),
    objects = cms.InputTag("muPtEtaIDIso"),
    )
process.selectedPatJetsNoMuons = cms.EDProducer(
    "JetsCleaner",
    jets =  cms.InputTag("selectedPatJets25"),
    valueMap = cms.InputTag("deltaRJetMuons"),
    minDeltaR = cms.double(0.3)
    )
process.deltaRJetTaus = cms.EDProducer(
    "DeltaRNearestTauComputer",
    probes = cms.InputTag("selectedPatJetsNoMuons"),
    objects = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
    )
process.selectedPatJetsNoMuonsNoTaus = cms.EDProducer(
    "JetsCleaner",
    jets =  cms.InputTag("selectedPatJetsNoMuons"),
    valueMap = cms.InputTag("deltaRJetTaus"),
    minDeltaR = cms.double(0.3)
    )
process.selectedPatJetsNoMuonsNoTausCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("selectedPatJetsNoMuonsNoTaus"),
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(999),
    )
process.diJet = cms.EDProducer(
    "CandViewShallowCloneCombiner",
    decay = cms.string("selectedPatJetsNoMuonsNoTaus selectedPatJetsNoMuonsNoTaus"),
    cut   = cms.string("mass > 300 && daughter(0).pt>25 && daughter(1).pt>25 && abs(daughter(0).eta-daughter(1).eta)>2.5 && daughter(0).eta*daughter(1).eta<0"),
    checkCharge = cms.bool(False),   
    )
process.diJetCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("diJet"),
    minNumber = cms.uint32(1)
    )
process.twoJetSequence = cms.Sequence(
    process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter*
    process.muPtEtaIDIso*process.muPtEtaIDIsoCounter*
    process.selectedPatJets25*process.deltaRJetMuons*process.selectedPatJetsNoMuons*
    process.deltaRJetTaus*process.selectedPatJetsNoMuonsNoTaus*
    process.selectedPatJetsNoMuonsNoTausCounter*
    process.diJet*process.diJetCounter
    )

#######################################################################


process.muTauStreamAnalyzer = cms.EDAnalyzer(
    "MuTauStreamAnalyzer",
    diTaus =  cms.InputTag("selectedDiTau"),
    jets =  cms.InputTag("selectedPatJets"),
    triggerResults = cms.InputTag("patTriggerEvent"),
    isMC = cms.bool(runOnMC),
    deltaRLegJet  = cms.untracked.double(0.5),
    minCorrPt = cms.untracked.double(15.),
    minJetID  = cms.untracked.double(0.5), # 1=loose,2=medium,3=tight
    applyTauSignalSel =  cms.bool( True ),
    verbose =  cms.untracked.bool( False ),
    )
process.muTauStreamAnalyzer2 = process.muTauStreamAnalyzer.clone(diTaus =  cms.InputTag("selectedDiTau2"))


process.analysis  = cms.Sequence(
    process.diTau*process.selectedDiTau* process.muTauStreamAnalyzer
    )
process.analysis2 = cms.Sequence(
    process.twoJetSequence*
    process.diTau2*process.selectedDiTau2*process.muTauStreamAnalyzer2
    )

process.p = cms.Path(
    process.analysis
    )
process.p2 = cms.Path(
    process.analysis2
    )

process.out = cms.OutputModule(
    "PoolOutputModule",
    outputCommands = cms.untracked.vstring( 'keep *'),
    fileName = cms.untracked.string('patTuplesSkimmed_MuTauStream.root'),
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("treeMuTauStream.root")
                                   )

process.outpath = cms.EndPath()

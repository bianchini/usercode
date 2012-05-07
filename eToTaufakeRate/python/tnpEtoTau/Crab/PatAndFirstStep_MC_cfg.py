import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
from PhysicsTools.TagAndProbe.Crab.customizePAT import *
from PhysicsTools.PatAlgos.tools.coreTools import *
process.load("PhysicsTools.TagAndProbe.tags_cff")
process.load("PhysicsTools.TagAndProbe.probes_cff")
process.load("PhysicsTools.TagAndProbe.tagandprobes_cff")
process.load("PhysicsTools.TagAndProbe.onetp_cff")
process.load("PhysicsTools.TagAndProbe.mcMatch_cff")
process.load("PhysicsTools.TagAndProbe.treeProducers_cff")

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
process.load("Configuration.StandardSequences.MagneticField_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/5E6BB841-A7F8-E011-95AE-485B39800C1B.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/905B0FB0-A1F8-E011-9293-E0CB4E19F9BD.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/1EB7FB80-B4F8-E011-89B1-001EC9D81A4A.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/BC4089C3-A8F8-E011-921A-00261834B548.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/0870F5FC-B9F8-E011-B97A-E0CB4EA0A8EA.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/12D19FBC-9EF8-E011-8FA3-90E6BAE8CC1C.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/167A26E4-A7F8-E011-B653-00261834B5B1.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/0E9F1594-96F8-E011-9F0A-E0CB4E1A118D.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/58FFDB31-97F8-E011-B89E-485B398971EA.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/4A36FE0C-ADF8-E011-8A08-20CF300E9ECF.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/12A2EBEA-ADF8-E011-B8F7-485B39800BDA.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/BAC098B6-AAF8-E011-949F-90E6BA0D099E.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/C6FBA935-B0F8-E011-9567-00261834B529.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/B864105B-B5F8-E011-B346-E0CB4E19F982.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/364EE1D2-9BF8-E011-AE0D-90E6BA0D0989.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/6C4A5277-AEF8-E011-BAEE-0022198F5AF7.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/50056D35-AAF8-E011-9176-90E6BA0D0996.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/FE743448-A3F8-E011-8BE4-E0CB4E29C513.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/9C8E989B-A7F8-E011-B90C-485B39800C17.root',
	'/store/mc/Fall11/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/72C689D4-BAF8-E011-B8DD-90E6BA442F0C.root',

    )
    
)

## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

process.load('RecoJets.Configuration.RecoPFJets_cff')

#--------------------------------------------------------------------------------

from PhysicsTools.PatAlgos.tools.jetTools import *

jec = [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]
#if not isMC:
#        jec.extend([ 'L2L3Residual' ])
addJetCollection(process, cms.InputTag('ak5PFJets'),
     'AK5', 'PF',
     doJTA            = True,
     doBTagging       = True,
     jetCorrLabel     = ('AK5PF', cms.vstring(jec)),
     doType1MET       = False,
     doL1Cleaning     = True,
     doL1Counters     = False,
     genJetCollection = cms.InputTag("ak5GenJets"),
     doJetID          = True,
     jetIdLabel       = "ak5",
     outputModule     = ''
)

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# configure Jet Energy Corrections
#
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.jec = cms.ESSource("PoolDBESSource",
     DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
     ),
     timetype = cms.string('runnumber'),
     toGet = cms.VPSet(
       cms.PSet(
           record = cms.string('JetCorrectionsRecord'),
           tag    = cms.string('JetCorrectorParametersCollection_Jec11V2_AK5PF'),
           label  = cms.untracked.string('AK5PF')
       ),
       cms.PSet(
           record = cms.string('JetCorrectionsRecord'),
           tag    = cms.string('JetCorrectorParametersCollection_Jec11V2_AK5Calo'),
           label  = cms.untracked.string('AK5Calo')
       )
    ),
    connect = cms.string('sqlite_fip:TauAnalysis/Configuration/data/Jec11V2.db')
)
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')

#--------------------------------------------------------------------------------

process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.Rho_EtaMax = cms.double(4.4)
#process.kt6PFJets.Ghost_EtaMax = cms.double(5.0)
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.Rho_EtaMax = cms.double(4.4)
#process.ak5PFJets.Ghost_EtaMax = cms.double(5.0)

## re-run kt4PFJets within lepton acceptance to compute rho
process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJetsCentral = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsCentral.Rho_EtaMax = cms.double(2.5)

process.fjSequence = cms.Sequence(process.kt6PFJets+process.ak5PFJets+process.kt6PFJetsCentral)

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process) # For HPS Taus
#switchToPFTauHPSpTaNC(process) # For HPS TaNC Taus

addSelectedPFlowParticle(process)

process.electronVariables.isMC = cms.bool(True)

addTriggerMatchingElectron(process,isMC=True)
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import *
patTrigger.processName    = cms.string( "HLT" )
process.eleTriggerMatchHLTElectrons.matchedCuts = cms.string('path("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v*")')

process.etoTau.isMC = cms.bool(True)
process.etoTau.makeMCUnbiasTree = cms.bool(True)
process.etoTau.checkMotherInUnbiasEff = cms.bool(True)
process.addUserVariables.isMC = cms.bool(True)

# require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter = cms.untracked.bool(True),
                                    debugOn = cms.untracked.bool(False),
                                    numtrack = cms.untracked.uint32(10),
                                    thresh = cms.untracked.double(0.2)
                                    )

process.p = cms.Path(process.scrapingVeto *
    		     process.PFTau *
    		     process.fjSequence *
		     process.patDefaultSequence *
    		     process.electronVariables *
		     process.addUserVariables *
		     (process.tagAnyEle + process.passingProbes) *
        	     process.tnpAnyEleAnyTau *
        	     process.oneTp *
        	     (process.tagMcMatch + process.probeMcMatch) *
        	     process.etoTau 
)

################################################################################################
###    P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s   f r o m    A O D  ###
################################################################################################

## pat sequences to be loaded:
#process.load("PhysicsTools.PFCandProducer.PF2PAT_cff")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
#process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
                       
# load the coreTools of PAT
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## modify the final pat sequence: keep only electrons + METS (muons are needed for met corrections)
process.load("RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequence_cff")
#process.patElectronIsolation = cms.Sequence(process.egammaIsolationSequence)

process.patElectrons.isoDeposits = cms.PSet()
process.patElectrons.userIsolation = cms.PSet()
process.patElectrons.addElectronID = cms.bool(True)
process.patElectrons.electronIDSources = cms.PSet(
    simpleEleId95relIso= cms.InputTag("simpleEleId95relIso"),
    simpleEleId90relIso= cms.InputTag("simpleEleId90relIso"),
    simpleEleId85relIso= cms.InputTag("simpleEleId85relIso"),
    simpleEleId80relIso= cms.InputTag("simpleEleId80relIso"),
    simpleEleId70relIso= cms.InputTag("simpleEleId70relIso"),
    simpleEleId60relIso= cms.InputTag("simpleEleId60relIso"),
    simpleEleId95cIso= cms.InputTag("simpleEleId95cIso"),
    simpleEleId90cIso= cms.InputTag("simpleEleId90cIso"),
    simpleEleId85cIso= cms.InputTag("simpleEleId85cIso"),
    simpleEleId80cIso= cms.InputTag("simpleEleId80cIso"),
    simpleEleId70cIso= cms.InputTag("simpleEleId70cIso"),
    simpleEleId60cIso= cms.InputTag("simpleEleId60cIso"),    
    )
##
process.patElectrons.addGenMatch = cms.bool(False)
process.patElectrons.embedGenMatch = cms.bool(False)
process.patElectrons.usePV = cms.bool(False)
##
process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
# you have to tell the ID that it is data
process.simpleEleId95relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId90relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId85relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId80relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId70relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId60relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId95cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId90cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId85cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId80cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId70cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId60cIso.dataMagneticFieldSetUp = cms.bool(True)
#
process.patElectronIDs = cms.Sequence(process.simpleEleIdSequence)
process.makePatElectrons = cms.Sequence(process.patElectronIDs*process.patElectrons)
# process.makePatMuons may be needed depending on how you calculate the MET
#process.makePatCandidates = cms.Sequence(process.makePatElectrons+process.makePatMETs)
#process.patDefaultSequence = cms.Sequence(process.makePatCandidates)
##
##  ################################################################################

addPFMuonIsolation(process,process.patMuons)
addPFElectronIsolation(process,process.patElectrons)

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple.root'),
                               # save only events passing the full path
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               # save PAT Layer 1 output; you need a '*' to
                               # unpack the list of commands 'patEventContent'
                               outputCommands = cms.untracked.vstring('keep *', *patEventContent )
                               )

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )
process.patTriggerEvent.processName = '*'

if hasattr(process,"patTrigger"):
    process.patTrigger.processName = '*'

## remove MC matching from the default sequence
#removeMCMatching(process, ['All'])
#removeMCMatching(process, ['METs'], "TC")
#removeMCMatching(process, ['METs'], "PF")
#runOnData(process)

process.TFileService = cms.Service("TFileService", fileName = cms.string("testTagAndProbe.root") )

#process.outpath = cms.EndPath(process.out)


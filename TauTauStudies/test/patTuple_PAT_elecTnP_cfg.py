from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 2000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

process.source.fileNames = cms.untracked.vstring(
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat//store/mc/Summer11/WH_ZH_TTH_HToWW_M-120_7TeV-pythia6//AODSIM/PU_S3_START42_V11-v1//0000/E430E306-347C-E011-A75A-00261834B5C6.root'
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/mc/Summer11/DYToEE_M-120_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S3_START42_V11-v2/0000/74B7A773-8F88-E011-84A7-1CC1DE1D03DE.root'
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/data/Run2011A/SingleElectron/AOD/May10ReReco-v1/0000/FEFE0831-A27B-E011-B146-00266CF32EAC.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/data/Run2011A/SingleElectron/AOD/May10ReReco-v1/0000/FE6D4370-AA7B-E011-9B9D-0025901D4C74.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/data/Run2011A/SingleElectron/AOD/May10ReReco-v1/0000/FE321C2B-A27B-E011-8CFB-00266CF327E0.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/data/Run2011A/SingleElectron/AOD/May10ReReco-v1/0000/FCD648BB-AC7B-E011-ADF8-0025901D4764.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/data/Run2011A/SingleElectron/AOD/May10ReReco-v1/0000/FC922AA3-957B-E011-8C6C-00266CF33318.root',
    )

postfix           = "PFlow"
runOnMC           =  True
FileName          = "treeElecTnP.root"

if runOnMC:
    process.GlobalTag.globaltag = cms.string('START42_V13::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V19::All')


process.primaryVertexFilter = cms.EDFilter(
    "GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
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

addPFElectronIsolation(process,process.patElectrons)

getattr(process,"patElectrons").embedTrack = True
getattr(process,"patElectrons").embedGsfTrack = True
addTriggerMatchingElectron(process,isMC=runOnMC)
process.eleTriggerMatchHLTElectrons.matchedCuts = cms.string('( (path("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",0,0) || path("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",0,0) ) || filter("hltEle15CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter") || filter("hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter"))  && type("TriggerElectron")')

if hasattr(process,"patTrigger"):
    process.patTrigger.processName = '*'


process.selectedPatElectronsTriggerMatchUserEmbedded = cms.EDProducer(
    "ElectronsUserEmbedded",
    electronTag = cms.InputTag("selectedPatElectronsTriggerMatch"),
    vertexTag = cms.InputTag("offlinePrimaryVertices"),
    isMC = cms.bool(runOnMC)
    )

#####################################################################################
#####################################################################################
process.load("Bianchi.TauTauStudies.elecTnP_cff")

if not runOnMC:
    process.sequence.remove(process.tagMcMatch)
    process.sequence.remove(process.probeMcMatch)
    process.elecTnP.isMC = cms.bool(False)
    process.elecTnP.makeMCUnbiasTree = cms.bool(False)
    process.addUserVariables.isMC = cms.bool(False)


process.pat = cms.Sequence(
    process.allEventsFilter+
    process.primaryVertexFilter*
    process.patDefaultSequence*
    process.selectedPatElectronsTriggerMatchUserEmbedded*
    process.sequence
    )

process.p = cms.Path(process.pat)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(FileName)
                                   )

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

process.out.fileName = cms.untracked.string('patTuples_elecTnP.root')
process.outpath = cms.EndPath()

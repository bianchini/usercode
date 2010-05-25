# This is an example PAT configuration showing the usage of PF2PAT+PAT

from PhysicsTools.PatAlgos.patTemplate_cfg import *

#process.GlobalTag.globaltag = cms.string('GR10_P_V2::All')
process.GlobalTag.globaltag = cms.string('GR10_P_V4::All')

# utilities to print the generator information
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree1 = cms.EDAnalyzer("ParticleListDrawer",
                                    src = cms.InputTag("genParticles"),
                                    maxEventsToPrint  = cms.untracked.int32(20)
                                    )

process.printTree2 = cms.EDAnalyzer("ParticleTreeDrawer",
                                    src = cms.InputTag('genParticles'),
                                    printP4 = cms.untracked.bool(False),
                                    printPtEtaPhi =  cms.untracked.bool(True),
                                    printVertex = cms.untracked.bool(False),
                                    printStatus = cms.untracked.bool(False),
                                    printIndex = cms.untracked.bool(False),
                                    status = cms.untracked.vint32( 3, 2 )
                                    )



readFiles = cms.untracked.vstring()
readFilesReReco = cms.untracked.vstring()


readFiles.extend([
    'rfio:/castor/cern.ch/user/b/bianchi/Zee_2/ZeeSummer09_1.root',
    'rfio:/castor/cern.ch/user/b/bianchi/Zee_2/ZeeSummer09_2.root',
    'rfio:/castor/cern.ch/user/b/bianchi/Zee_2/ZeeSummer09_3.root',
    'rfio:/castor/cern.ch/user/b/bianchi/Zee_2/ZeeSummer09_4.root',
    'rfio:/castor/cern.ch/user/b/bianchi/Zee_2/ZeeSummer09_5.root']
                 )

readFilesReReco.extend([
    #'rfio:/castor/cern.ch/user/b/bianchi/CMSSW357/reco/recoFileReco.root'
    'rfio:/castor/cern.ch/user/b/bianchi/CMSSW357/rereco/recoFileReReco.root'
    ]
                       )

process.source = cms.Source("PoolSource", 
                            #fileNames = readFiles
                            fileNames = readFilesReReco
)


# process.load("PhysicsTools.PFCandProducer.Sources.source_ZtoMus_DBS_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.out.fileName = cms.untracked.string('patLayer_PF2PAT_ZCAND_RERECO.root')

# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.pfTools import *

usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=False) 
removeMCMatching(process, ['All'])
switchJECSet( process, "Summer09_7TeV_ReReco332") # Data

process.patElectrons.embedTrack = True

# remove pt threshold from pfIsolatedElectrons
process.pfElectronsPtGt5.ptMin = 0.0

# remove isolation from PAT
process.pfIsolatedElectrons.combinedIsolationCut = cms.double(9999)
process.pfIsolatedElectrons.isolationCuts = cms.vdouble(9999,9999)

# remove vertex constraint on isoDepsoits
process.isoDepElectronWithCharged.ExtractorPSet.Diff_r = 9999.99
process.isoDepElectronWithCharged.ExtractorPSet.Diff_z = 9999.99

# https://twiki.cern.ch/twiki/bin/view/CMS/L1TechnicalTriggerBits
# 0: BPTX+ and BPTX- ;
# 40: BSC minbias: hits +z side >= 1 and hits -z side >= 1 (all segments: inner rings, outer petals)
# 41: BSC minbias: hits +z side >= 2 and hits -z side >= 2 (all segments: inner rings, outer petals)
# 36/37/38/39: BSC halo: beam 2/1 inner/outer
# 42/43: BSC splash trigger: beam 1/2, threshold = 2; inner ring hits on -z side >= 2 
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
#process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('((40 OR 41) AND NOT (36 OR 37 OR 38 OR 39))')
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39) AND NOT ((42 AND NOT 43) OR (43 AND NOT 42))')


#DPGAnalysis/Skims/python/skim_noscrape_cfg.py
process.skim_noscrape = cms.EDFilter("FilterOutScraping",
                                 applyfilter = cms.untracked.bool(True),
                                 debugOn = cms.untracked.bool(False),
                                 numtrack = cms.untracked.uint32(10),
                                 thresh = cms.untracked.double(0.2)
                                 )

# DPGAnalysis/Skims/python/GoodVertex_cfg.py
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,                                           
                                           maxAbsZ = cms.double(15), 
                                           maxd0 = cms.double(2) )
                                           

# DPGAnalysis/Skims/python/skim_physdecl_cfg.py
process.skim_physdecl = cms.EDFilter("PhysDecl",
     applyfilter = cms.untracked.bool(True)
)


# filter based on a FWLite analyzer
process.ZeeFilterAnalyzer = cms.EDFilter("EDZeeFilterAnalyzer",
                                         cfgFileName=cms.untracked.string("Demo_cfg.py")
                                         )

process.electronsCountFilter = cms.EDFilter("PATCandViewCountFilter",
                                            minNumber = cms.uint32(2),
                                            maxNumber = cms.uint32(999999),
                                            src = cms.InputTag("selectedPatElectrons")
                                            )

# Let it run
process.p = cms.Path(
    process.hltLevel1GTSeed +
    process.skim_noscrape + 
    process.skim_physdecl +
    process.primaryVertexFilter +
    process.patDefaultSequence +
    process.electronsCountFilter +
    process.ZeeFilterAnalyzer
    #process.printTree1 +
    #process.printTree2
)


# Add PF2PAT output to the created file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
from Configuration.EventContent.EventContent_cff import AODEventContent
#process.load("PhysicsTools.PFCandProducer.PF2PAT_EventContent_cff")

process.out.outputCommands = AODEventContent.outputCommands

additionalOutput = cms.vstring(
    'drop *_*_*_PAT',
    *patEventContentNoCleaning
    )

process.out.outputCommands.extend( additionalOutput )

#'keep *_offlinePrimaryVerticesWithBS_*_*',                                                                                             
#'keep *_offlinePrimaryVertices_*_*',     
#'keep *_particleFlow_*_*',               
#'keep *_generalTracks_*_*',              
#'keep *_TriggerResults__HLT',
#*patEventContentNoCleaning


process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )


process.MessageLogger.cerr.FwkReport.reportEvery = 10

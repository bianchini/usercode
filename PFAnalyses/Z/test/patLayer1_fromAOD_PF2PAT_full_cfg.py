# This is an example PAT configuration showing the usage of PF2PAT+PAT

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *


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
    #'rfio:/castor/cern.ch/user/b/bianchi/Zee_2/ZeeSummer09_2.root',
    #'rfio:/castor/cern.ch/user/b/bianchi/Zee_2/ZeeSummer09_3.root',
    #'rfio:/castor/cern.ch/user/b/bianchi/Zee_2/ZeeSummer09_4.root',
    #'rfio:/castor/cern.ch/user/b/bianchi/Zee_2/ZeeSummer09_5.root'
    ]
                 )

readFilesReReco.extend([
    '/store/relval/CMSSW_3_5_7/RelValZEE/GEN-SIM-RECO/START3X_V26-v1/0012/020A72FB-4749-DF11-A27E-003048679076.root'
    ]
                       )

process.source = cms.Source("PoolSource", 
                            #fileNames = readFiles
                            fileNames = readFilesReReco
)


process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.out.fileName = cms.untracked.string('/tmp/bianchi/patLayer1_fromAOD_PF2PAT_full.root')

process.GlobalTag.globaltag = 'MC_3XY_V26::All'

#process.out.fileName = cms.untracked.string('/tmp/bianchi/patLayer/patLayer1_fromAOD_PF2PAT_full_ReReco.root')

# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.pfTools import *


usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=True) 
switchJECSet( process, "Summer09_7TeV_ReReco332") # Data

process.pfElectronsPtGt5.ptMin = 0.0

# remove isolation from PAT
process.pfIsolatedElectrons.combinedIsolationCut = cms.double(9999)
process.pfIsolatedElectrons.isolationCuts = cms.vdouble(9999,9999)


# remove vertex constraints
process.isoDepElectronWithCharged.ExtractorPSet.Diff_r = 9999.99
process.isoDepElectronWithCharged.ExtractorPSet.Diff_z = 9999.99

# turn to false when running on data
process.patElectrons.embedGenMatch = True
process.patElectrons.embedGsfTrack = True 
process.patElectrons.embedTrack    = True

# customize matching
process.electronMatch.maxDeltaR = 0.15 # def 0.5
process.electronMatch.maxDPtRel = 0.5 # def 0.5

# filter based on a FWLite analyzer
process.ZeeFilterAnalyzer = cms.EDFilter("EDZeeFilterAnalyzer",
                                         cfgFileName=cms.untracked.string("Demo_cfg.py")
                                         )


process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
#process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('((40 OR 41) AND NOT (36 OR 37 OR 38 OR 39))')
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39) AND NOT ((42 AND NOT 43) OR (43 AND NOT 42))')

process.electronsCountFilter = cms.EDFilter("PATCandViewCountFilter",
                                            minNumber = cms.uint32(2),
                                            maxNumber = cms.uint32(999999),
                                            src = cms.InputTag("selectedPatElectrons")
                                            )


# Let it run

process.p = cms.Path(
    #process.hltLevel1GTSeed +
    process.patDefaultSequence +
    process.electronsCountFilter +
    process.ZeeFilterAnalyzer
    # + process.printTree1
    # + process.printTree2
)


# Add PF2PAT output to the created file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
#process.load("PhysicsTools.PFCandProducer.PF2PAT_EventContent_cff")

'''
process.out.outputCommands = cms.untracked.vstring('drop *',
                                                   'keep recoGenParticles_genParticles_*_*',
                                                   'keep *_offlinePrimaryVerticesWithBS_*_*',
                                                   'keep *_offlinePrimaryVertices_*_*',
                                                   'keep *_particleFlow_*_*',
                                                   'keep *_generalTracks_*_*',
                                                   'keep *_TriggerResults__HLT',
                                                   *patEventContentNoCleaning
                                                   )
'''
from Configuration.EventContent.EventContent_cff import AODEventContent
process.out.outputCommands = AODEventContent.outputCommands

additionalOutput = cms.vstring(
    'drop *_*_*_PAT',
    *patEventContentNoCleaning
    )
process.out.outputCommands.extend( additionalOutput )


process.MessageLogger.cerr.FwkReport.reportEvery = 10



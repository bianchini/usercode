import FWCore.ParameterSet.Config as cms

process = cms.Process("Zee")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

readFiles1 = cms.untracked.vstring()
readFiles2 = cms.untracked.vstring()
readFiles3 = cms.untracked.vstring()

readFiles1.extend([
    'file:./ZCandidate/rereco/patLayer_PF2PAT_ZCAND_RERECO.root'
    #'file:./ZCandidate/reco/patLayer_PF2PAT_ZCAND_RECO.root'
    ]
                  )

readFiles2.extend([
    'file:./Data/patLayer_PF2PAT_SD_EG_1_1.root',
    'file:./Data/patLayer_PF2PAT_SD_EG_2_1.root',
    'file:./Data/patLayer_PF2PAT_SD_EG_19_1.root',
    'file:./Data/patLayer_PF2PAT_SD_EG_21_1.root',
    'file:./Data/patLayer_PF2PAT_SD_EG_16_1.root'
    ]
                  )

readFiles3.extend([
    'file:./MC/patLayer_PF2PAT_SD_EG_10_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_11_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_12_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_16_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_17_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_18_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_19_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_1_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_20_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_21_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_24_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_25_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_26_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_28_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_29_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_2_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_30_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_31_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_33_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_35_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_36_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_38_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_39_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_42_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_44_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_45_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_48_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_49_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_4_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_51_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_52_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_54_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_55_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_57_1.root', 
    'file:./MC/patLayer_PF2PAT_SD_EG_9_1.root'
    ]
                  )

# data from 05 May on GOODCOLL_v9
from PFAnalyses.Z.data05May_cff import * 

process.source = cms.Source("PoolSource",
                            fileNames = readFiles1
                            #fileNames = readFiles2
                            #fileNames = readFiles3
                            #fileNames = readThisData
                            )


process.configurationMetadata = cms.untracked.PSet(
     version = cms.untracked.string('1.0'),
     name = cms.untracked.string('ReReco'),
     annotation = cms.untracked.string('None'),
     preselectionEff = cms.untracked.double(1.0)
 )

process.Summary =  cms.PSet(
    selectionFlavours = cms.untracked.vstring("")
    )

process.PatZeeAnalyzer = cms.PSet(

    patElectronsLabel = cms.InputTag('selectedPatElectrons'),
    offlinePrimaryVertecesLabel = cms.InputTag("offlinePrimaryVertices"),
    triggerResultsLabel = cms.InputTag("TriggerResults","","HLT"),

    useTrigger = cms.bool( True ),
    triggerItemNames = cms.vstring('HLT_Photon10_L1R',
                                   'HLT_Photon15_L1R',
                                   'HLT_Photon15_LooseEcalIso_L1R',
                                   'HLT_Photon20_L1R',
                                   'HLT_Photon30_L1R_8E29',
                                   'HLT_DoublePhoton4_Jpsi_L1R',
                                   'HLT_DoublePhoton4_Upsilon_L1R',
                                   'HLT_DoublePhoton4_eeRes_L1R',
                                   'HLT_DoublePhoton5_Jpsi_L1R',
                                   'HLT_DoublePhoton5_Upsilon_L1R',
                                   'HLT_DoublePhoton5_L1R',
                                   'HLT_DoublePhoton10_L1R',
                                   'HLT_DoubleEle5_SW_L1R',
                                   'HLT_Ele20_LW_L1R',
                                   'HLT_Ele15_SiStrip_L1R',
                                   'HLT_Ele15_SC10_LW_L1R',
                                   'HLT_Ele15_LW_L1R',
                                   'HLT_Ele10_LW_EleId_L1R',
                                   'HLT_Ele10_LW_L1R',
                                   'HLT_Photon15_TrackIso_L1R'),
    triggerItemAND = cms.bool( False ),

    verbose = cms.bool( False ),
    triggerVerbose = cms.bool( False ),
    scanEvsByM = cms.bool( True ),
    minMass = cms.double(60),
    maxMass = cms.double(120),

    minPtSubSub = cms.double(5.),

    elePhaseSpaceSelector = cms.PSet(
       pt = cms.VPSet(
            cms.PSet( min = cms.double(3.0), max = cms.double(9999.) ), #3.0
            ),
       eta = cms.VPSet(
            cms.PSet( min = cms.double(-2.4), max = cms.double(2.4) ), #2.4
            ),
       phi = cms.VPSet(
            cms.PSet( min = cms.double(-10), max = cms.double(10) ) 
            ),
       verbose = cms.untracked.bool( False )
       ),
    
    eleSelector = cms.PSet(
       version = cms.int32(0),
       mvaCut = cms.double(-0.3), #-0.3
       maxMissingHitCut = cms.double(999),
       verbose = cms.untracked.bool( False )
       ),

    eleIsoSelector = cms.PSet(
       verbose = cms.untracked.bool( False ),
       isolationCuts = cms.VPSet(
        cms.PSet(
          name = cms.string("combined"), 
          absoluteCut = cms.double(999.0), #3.0
          relativeCut = cms.double(0.15),
          chargedWeight = cms.double(1.0),
          photonWeight = cms.double(1.0),
          neutralWeight = cms.double(0.0)
          ),
        cms.PSet(
          name = cms.string("charged"), 
          absoluteCut = cms.double(5),
          relativeCut = cms.double(0.14),
          chargedWeight = cms.double(1),
          photonWeight = cms.double(0),
          neutralWeight = cms.double(0)
          )
       ),
      chargedIsoPar = cms.PSet(
        coneSize = cms.double(0.3),
        vetoes = cms.VPSet(
          cms.PSet(
           type = cms.string("ThresholdVeto"),
           threshold = cms.double(0.0)  # threshold on deposit
          )
        )
       ),
       neutralIsoPar = cms.PSet(
        coneSize = cms.double(0.3),
        vetoes = cms.VPSet(
          cms.PSet(
           type = cms.string("ThresholdVeto"),
           threshold = cms.double(0.5)
          ),
          cms.PSet(
           type = cms.string("ConeVeto"),
           deltaR = cms.double(0.00)  # radius of cone veto
          )
        )
       ),
       photonsIsoPar = cms.PSet(
        coneSize = cms.double(0.3), 
        vetoes = cms.VPSet(
          cms.PSet(
           type = cms.string("ThresholdVeto"),
           threshold = cms.double(0.5)
          ),
          cms.PSet(
           type = cms.string("RectangularEtaPhiVeto"),
           dEta = cms.double(0.00), # eta width of the strip to be vetoed
           dPhi = cms.double(0.00)   # phi width of the strip to be vetoed
           )
          )
       )
       
       )
    
    )



process.PatElectronAnalyzer = cms.PSet(
    genParticlesLabel = cms.InputTag('genParticles'),
    patElectronsLabel = cms.InputTag('selectedPatElectrons'),

    patPfCandidateLabel = cms.InputTag('selectedPatPFParticles'),
    tracksLabel = cms.InputTag('generalTracks'),

    offlinePrimaryVertecesLabel = cms.InputTag("offlinePrimaryVertices"),

    minElePt = cms.double(15.),
    maxEleEta = cms.double(3.),
    minMVA = cms.untracked.double(-0.3),
    maxDxy = cms.untracked.double(0.2),
    maxDz = cms.untracked.double(1.0),
    minNumValidHits = cms.untracked.int32(10),

    runOnMC = cms.bool(True),

    ElectronSelector = cms.PSet(
       verbose = cms.untracked.bool( False ),

       isolationCuts = cms.VPSet(
        cms.PSet(
          name = cms.string("combinedFlorent"), 
          absoluteCut = cms.double(5),
          relativeCut = cms.double(0.04),
          chargedWeight = cms.double(1),
          photonWeight = cms.double(1),
          neutralWeight = cms.double(0.33)
          ),
        cms.PSet(
          name = cms.string("charged"), 
          absoluteCut = cms.double(5),
          relativeCut = cms.double(0.14),
          chargedWeight = cms.double(1),
          photonWeight = cms.double(0),
          neutralWeight = cms.double(0)
          )
       ),
       
       chargedIsoPar = cms.PSet(
        coneSize = cms.double(0.4),
        vetoes = cms.VPSet(
         cms.PSet(
          type = cms.string("ThresholdVeto"),
          threshold = cms.double(0.5)
         )
        )
       ),
       
       neutralIsoPar = cms.PSet(
        coneSize = cms.double(0.4),
        vetoes = cms.VPSet(
         cms.PSet(
          type = cms.string("ThresholdVeto"),
          threshold = cms.double(0.5)
          ),
         cms.PSet(
          type = cms.string("ConeVeto"),
          deltaR = cms.double(0.07)
          )
         )
        ),

       photonsIsoPar = cms.PSet(
        coneSize = cms.double(0.4), 
        vetoes = cms.VPSet(
         cms.PSet(
          type = cms.string("ThresholdVeto"),
          threshold = cms.double(0.5)
          ),
         cms.PSet(
          type = cms.string("RectangularEtaPhiVeto"),
          dEta = cms.double(0.02),
          dPhi = cms.double(0.3)
          )
         )
        ) 
       )
    
    )



process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printTree = cms.EDAnalyzer("ParticleListDrawer",
     src = cms.InputTag("genParticles"),
     maxEventsToPrint  = cms.untracked.int32(-1)
)

process.ZeeFilterAnalyzer = cms.EDFilter("EDZeeFilterAnalyzer",
                                         cfgFileName=cms.untracked.string("Demo_cfg.py"))

process.EleSelectorFilter = cms.EDFilter(
    "PatElectronSelectorFilter",
    src = cms.InputTag("selectedPatElectrons"),
    selector = process.PatZeeAnalyzer.eleSelector
    )


process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    )

process.p = cms.Path(
    #process.ZeeFilterAnalyzer+
    #process.EleSelectorFilter
    process.printTree
    )

process.out = cms.OutputModule(
    "PoolOutputModule",
    outputCommands = cms.untracked.vstring( 'keep *'),
    fileName = cms.untracked.string('filteredPatLayerZee.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring("p")
    )

)

process.outpath = cms.EndPath(process.out )




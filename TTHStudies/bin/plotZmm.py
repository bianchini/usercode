import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms


process = cms.Process("FWLitePlots")

process.fwliteInput = cms.PSet(

    pathToFile    = cms.string("/store/HBB_EDMNtuple/V42/Oct22/env/sys/MVAout/"),
    ordering      = cms.string("ZllH.DiJetPt.Oct22."),
    lumi          = cms.double(12.1),
    

    samples       = cms.VPSet(

    cms.PSet(
    name     = cms.string('DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball'),
    nickName = cms.string('DYJets'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71*3)
    ),

    cms.PSet(
    name     = cms.string('TTJets_Merged'),
    nickName = cms.string('TTJets'),
    color    = cms.int32(5),
    xSec     = cms.double(234)
    ),

    cms.PSet(
    name     = cms.string('T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('TtW'),
    color    = cms.int32(6),
    xSec     = cms.double(11.1)
    ),

    cms.PSet(
    name     = cms.string('T_t-channel_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('Tt'),
    color    = cms.int32(6),
    xSec     = cms.double(56.4)
    ),

    cms.PSet(
    name     = cms.string('T_s-channel_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('Ts'),
    color    = cms.int32(6),
    xSec     = cms.double(3.79)
    ),

    cms.PSet(
    name     = cms.string('Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('TbartW'),
    color    = cms.int32(6),
    xSec     = cms.double(11.1)
    ),

    cms.PSet(
    name     = cms.string('Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('Tbart'),
    color    = cms.int32(6),
    xSec     = cms.double(30.7)
    ),

    cms.PSet(
    name     = cms.string('Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('Tbars'),
    color    = cms.int32(6),
    xSec     = cms.double(1.76)
    ),

    cms.PSet(
    name     = cms.string('WW_TuneZ2star_8TeV_pythia6_tauola'),
    nickName = cms.string('WW'),
    color    = cms.int32(4),
    xSec     = cms.double(56.75)
    ),

    cms.PSet(
    name     = cms.string('WZ_TuneZ2star_8TeV_pythia6_tauola'),
    nickName = cms.string('WZ'),
    color    = cms.int32(3),
    xSec     = cms.double(33.85)
    ),

    cms.PSet(
    name     = cms.string('ZZ_TuneZ2star_8TeV_pythia6_tauola'),
    nickName = cms.string('ZZ'),
    color    = cms.int32(8),
    xSec     = cms.double(8.297)
    ),

    cms.PSet(
    name     = cms.string('ZH_ZToLL_HToBB_M-110_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH110'),
    color    = cms.int32(2),
    xSec     = cms.double(0.04414)
    ),

    cms.PSet(
    name     = cms.string('ZH_ZToLL_HToBB_M-115_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH115'),
    color    = cms.int32(2),
    xSec     = cms.double(0.036375)
    ),

    cms.PSet(
    name     = cms.string('ZH_ZToLL_HToBB_M-120_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH120'),
    color    = cms.int32(2),
    xSec     = cms.double(0.0293327854)
    ),

    cms.PSet(
    name     = cms.string('ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH125'),
    color    = cms.int32(2),
    xSec     = cms.double(0.0229727058)
    ),

    cms.PSet(
    name     = cms.string('ZH_ZToLL_HToBB_M-130_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH130'),
    color    = cms.int32(2),
    xSec     = cms.double(0.017288657)
    ),

    cms.PSet(
    name     = cms.string('ZH_ZToLL_HToBB_M-135_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH135'),
    color    = cms.int32(2),
    xSec     = cms.double(0.01250888)
    ),

    cms.PSet(
    name     = cms.string('DataZmm'),
    nickName = cms.string('DataZmm'),
    color    = cms.int32(1),
    xSec     = cms.double(-1)
    ),

    cms.PSet(
    name     = cms.string('DataZee'),
    nickName = cms.string('DataZee'),
    color    = cms.int32(1),
    xSec     = cms.double(-1)
    ),
    
    
    ),



    plots      = cms.VPSet(

    cms.PSet(
    xLow      = cms.double(10),
    xHigh     = cms.double(190),
    nBins     = cms.int32(90),
    variable  = cms.string("vLepton_pt[0]"),
    xTitle    = cms.string("p_{T} lead #mu"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_leadMuPt"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    xLow      = cms.double(10),
    xHigh     = cms.double(190),
    nBins     = cms.int32(90),
    variable  = cms.string("vLepton_pt[1]"),
    xTitle    = cms.string("p_{T} trail #mu"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_trailMuPt"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    xLow      = cms.double(0),
    xHigh     = cms.double(30),
    nBins     = cms.int32(30),
    variable  = cms.string("nPVs"),
    xTitle    = cms.string("vertices"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_vertex"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    xLow      = cms.double(0),
    xHigh     = cms.double(200),
    nBins     = cms.int32(50),
    variable  = cms.string("METtype1p2corr"),
    xTitle    = cms.string("MET"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_met"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    xLow      = cms.double(60),
    xHigh     = cms.double(200),
    nBins     = cms.int32(70),
    variable  = cms.string("V.mass"),
    xTitle    = cms.string("di#mu mass"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_diMuMass"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    xLow      = cms.double(-3),
    xHigh     = cms.double(3),
    nBins     = cms.int32(60),
    variable  = cms.string("vLepton_eta[0]"),
    xTitle    = cms.string("#eta lead #mu"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_leadMuEta"),
    logy      = cms.int32(0),
    ),

    
    cms.PSet(
    xLow      = cms.double(-3),
    xHigh     = cms.double(3),
    nBins     = cms.int32(60),
    variable  = cms.string("vLepton_eta[1]"),
    xTitle    = cms.string("#eta trail #mu"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_trailMuEta"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    xLow      = cms.double(20),
    xHigh     = cms.double(420),
    nBins     = cms.int32(100),
    variable  = cms.string("hJet_pt[0]"),
    xTitle    = cms.string("p_{T} lead jet"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_leadJetPt"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    xLow      = cms.double(20),
    xHigh     = cms.double(420),
    nBins     = cms.int32(100),
    variable  = cms.string("hJet_pt[1]"),
    xTitle    = cms.string("p_{T} trail jet"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_trailJetPt"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    xLow      = cms.double(0),
    xHigh     = cms.double(10),
    nBins     = cms.int32(10),
    variable  = cms.string("numJets"),
    xTitle    = cms.string("jet mult"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_numJet"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    xLow      = cms.double(0),
    xHigh     = cms.double(10),
    nBins     = cms.int32(10),
    variable  = cms.string("numBJets"),
    xTitle    = cms.string("b-jet mult"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_numBJet"),
    logy      = cms.int32(1),
    ),

    # cms.PSet(
    #xLow      = cms.double(),
    #xHigh     = cms.double(),
    #nBins     = cms.int32(),
    #xTitle    = cms.string(),
    #yTitle    = cms.string(),
    #histoName = cms.string(),
    #logy      = cms.int32(),
    #),

    ),


    fileList   = cms.vstring(
    'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball',
    'DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball',
    'DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball',
    'DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph',
    'DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph',
    'DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph',
    'DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph',
    'DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph',
    'DY3JetsToLL_M-50_TuneZ2Star_8TeV-madgraph',
    'DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph',
    'ZJetsToLL_Pt-100_8TeV-herwigpp',
    'TTJets_Merged',    
    'TToTlepWhad_tW-channel-DR_8TeV-powheg-tauola',
    'TToThadWlep_tW-channel-DR_8TeV-powheg-tauola',
    'TToDilepton_tW-channel-DR_8TeV-powheg-tauola',
    'TToLeptons_t-channel_8TeV-powheg-tauola',
    'TToLeptons_s-channel_8TeV-powheg-tauola',
    'T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola',
    'T_t-channel_TuneZ2star_8TeV-powheg-tauola',
    'T_s-channel_TuneZ2star_8TeV-powheg-tauola',
    'TBarToDilepton_tW-channel-DR_8TeV-powheg-tauola',
    'TBarToLeptons_t-channel_8TeV-powheg-tauola',
    'TBarToLeptons_s-channel_8TeV-powheg-tauola',   
    'Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola',
    'Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola',
    'Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola',
    'WW_TuneZ2star_8TeV_pythia6_tauola',
    'WZ_TuneZ2star_8TeV_pythia6_tauola',
    'WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola',
    'ZZ_TuneZ2star_8TeV_pythia6_tauola',
    'ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola',
    'ZH_ZToLL_HToBB_M-110_8TeV-powheg-herwigpp',
    'ZH_ZToLL_HToBB_M-115_8TeV-powheg-herwigpp',
    'ZH_ZToLL_HToBB_M-120_8TeV-powheg-herwigpp',
    'ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp',
    'ZH_ZToLL_HToBB_M-130_8TeV-powheg-herwigpp',
    'ZH_ZToLL_HToBB_M-135_8TeV-powheg-herwigpp',
    'DataZmm',
    'DataZee',
    ),
    
    xsection      = cms.vdouble()

    )

import FWCore.ParameterSet.Types  as CfgTypes
import FWCore.ParameterSet.Config as cms

VType     = "_VType2"
xsecTT_SL = 103.0

process = cms.Process("MEValidator")

process.fwliteInput = cms.PSet(

    outFileName   = cms.string("./root/MEValidator_SLW1j.root"),
    pathToTF      = cms.string("./root/transferFunctions_partonE.root"),
    pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt"+VType+"/v2/"),
    ordering      = cms.string("DiJetPt_"),
    lumi          = cms.double(12.1),

    samples       = cms.VPSet(


    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTH_HToBB_M-110_8TeV-pythia6'+VType),
    nickName = cms.string('TTH110'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1887*0.744)
    ),

    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTH_HToBB_M-115_8TeV-pythia6'+VType),
    nickName = cms.string('TTH115'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1663*0.703)
    ),


    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTH_HToBB_M-120_8TeV-pythia6'+VType),
    nickName = cms.string('TTH120'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1470*0.648)
    ),


    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTH_HToBB_M-125_8TeV-pythia6'+VType),
    nickName = cms.string('TTH125'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1302*0.569)
    ),


    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTH_HToBB_M-130_8TeV-pythia6'+VType),
    nickName = cms.string('TTH130'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1157*0.494)
    ),

    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTH_HToBB_M-135_8TeV-pythia6'+VType),
    nickName = cms.string('TTH135'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1031*0.404)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTWJets_8TeV-madgraph'+VType),
    nickName = cms.string('TTW'),
    color    = cms.int32(18),
    xSec     = cms.double(0.232),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTZJets_8TeV-madgraph'+VType),
    nickName = cms.string('TTZ'),
    color    = cms.int32(18),
    xSec     = cms.double(0.2057),
    ),


    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsSemiLept'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_SL),
    ),
    
    
    ),


    #SLNoBLep 10000
    #SLNoBHad 10000
    #SL2wj     2000
    #SL1wj     4000
    #DL       10000
    
    vegasPoints   = cms.int32(1000),

    mode          = cms.untracked.int32(0),

    norm          = cms.untracked.int32(2),
    functions     = cms.vstring('1.39e+23*x^(-2.81e+00)',
                                '7.84e+17*TMath::Landau(x,6.17e+01,1.61e+01)',
                                'x>=12 ? x^(-2.010e-01)*exp((-1.5785e-02)*x) : 4.184e-02*x'
                                ),

    switchoffOL   = cms.untracked.int32(0),

    useME         = cms.int32(1),
    useJac        = cms.int32(1),
    useMET        = cms.int32(1),
    useTF         = cms.int32(1),
    usePDF        = cms.int32(1),

    doParton        = cms.int32(1),
    doSmear         = cms.int32(0),   
    doMassScan      = cms.int32(1),
    doPermutations  = cms.int32(0),
    
    printP4         = cms.int32(0),
    
    verbose       = cms.bool(False),

    met           = cms.double(125.),
    masses        = cms.vdouble(60,  65,  70,  75, 80 , 85,  90, 95, 100, 105, 110, 
                                115, 120, 125, 130, 135, 140, 145, 150, 155, 
                                160, 165, 170, 175, 180 ,185, 190, 195, 200, 205, 210, 215, 220, 225, 230),
    #masses        = cms.vdouble(125),
    evLimits      = cms.vint32(1,1),

    scaleH        = cms.double(1.0),
    scaleL        = cms.double(1.0),
    scaleMET      = cms.double(1.0),

    )

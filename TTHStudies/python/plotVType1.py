import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

stitchDY = False

# trigger matching? ID? ISO?

mmBasic      = "(Vtype==1 && H.HiggsFlag==1 && vLepton_charge[0]*vLepton_charge[1]<0 && V.mass>60)"
mmKin        = "(vLepton_pt[0]>30 && vLepton_pt[1]>20 && abs(vLepton_eta[0])<2.1 && abs(vLepton_eta[1])<2.4)"
mmId         = "(vLepton_pfCorrIso[0]<0.10 && vLepton_pfCorrIso[1]<0.20)"
mmJetBasic   = "(numJets30>=2 && pt1>30 && pt2>30)"
common = mmBasic+" && "+mmKin+" && "+mmId+" && "+mmJetBasic

dataCut      = "((EVENT.json == 1 || EVENT.run < 196532) && (triggerFlags[5]>0 || triggerFlags[6]>0))"


#ttjetsLF = \
#         "(lheNj - 2 - (abs(genTop.wdau1id)<6 + abs(genTop.wdau2id)<6 + abs(genTbar.wdau1id)<6 + abs(genTbar.wdau2id)<6 ) )==0 "+ \
#         " || ((lheNj - 2 - (abs(genTop.wdau1id)<6 + abs(genTop.wdau2id)<6 + abs(genTbar.wdau1id)<6 + abs(genTbar.wdau2id)<6 ) )>0 && "+ \
#         " ( nC-nCTop == 0 && nB-nBTop == 0)" + \
#         ")"
ttjetsLF = \
         "nSimBs==2 && nC-nCTop==0"


#ttjetsC = \
#        "((lheNj - 2 - (abs(genTop.wdau1id)<6 + abs(genTop.wdau2id)<6 + abs(genTbar.wdau1id)<6 + abs(genTbar.wdau2id)<6 ) )>0 && "+ \
#        " ( nC-nCTop > 0 && nB-nBTop == 0)" + \
#        ")"
ttjetsC  = \
        "nSimBs==2 && nC-nCTop>0"


#ttjetsB = \
#        "((lheNj - 2 - (abs(genTop.wdau1id)<6 + abs(genTop.wdau2id)<6 + abs(genTbar.wdau1id)<6 + abs(genTbar.wdau2id)<6 ) )>0 && "+ \
#        " ( nC-nCTop >= 0 && nB-nBTop > 0)" + \
#        ")"
ttjetsB  = \
        "nSimBs>2"

#zjetsLF = \
#        "lheNj==0 || (lheNj>0 && nC==0 && nB==0)"
zjetsLF = \
        "nSimBs==0 && nC==0"

#zjetsC  = \
#       "lheNj>0 && nC>0 && nB==0"
zjetsC  = \
       "nSimBs==0 && nC>0"

#zjetsB  = \
#       "lheNj>0 && nC>=0 && nB>0"
zjetsB  = \
       "nSimBs>0"



process = cms.Process("FWLitePlots")

process.fwliteInput = cms.PSet(

    #pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/HBB_EDMNtuple/V42/Oct22/env/sys/MVAout/"),
    pathToFile    = cms.string("/scratch/bianchi/HBB_EDMNtuple/Zll.H.DiJetPt/"),
    ordering      = cms.string("ZllH.DiJetPt.Oct22."),
    lumi          = cms.double(12.1),


    samples       = cms.VPSet(

    cms.PSet(
    skip     = cms.bool(False),
    name     = cms.string('DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball'),
    nickName = cms.string('DYJets'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(not stitchDY),
    name     = cms.string('DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJets1'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(not stitchDY),
    name     = cms.string('DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJets2'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(not stitchDY),  
    name     = cms.string('DY3JetsToLL_M-50_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJets3'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(not stitchDY),  
    name     = cms.string('DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJets4'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(not stitchDY),  
    name     = cms.string('DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJetsHT1'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(not stitchDY),  
    name     = cms.string('DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJetsHT2'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    ##cms.PSet(
    ##name     = cms.string('DYJetsToLL_PtZ-180_TuneZ2star_8TeV-madgraph'),
    ##nickName = cms.string('DYJetsPt4'),
    ##color    = cms.int32(18),
    ##xSec     = cms.double(3503.71)
    ##),
    
    cms.PSet(
    skip     = cms.bool(not stitchDY),  
    name     = cms.string('DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph'),
    nickName = cms.string('DYJetsPt3'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),
    
    cms.PSet(
    skip     = cms.bool(not stitchDY),  
    name     = cms.string('DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball'),
    nickName = cms.string('DYJetsPt2'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(not stitchDY),  
    name     = cms.string('DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball'),
    nickName = cms.string('DYJetsPt1'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_Merged'),
    nickName = cms.string('TTJets'),
    color    = cms.int32(5),
    xSec     = cms.double(234)
    ),


   cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_HadronicMGDecays_8TeV-madgraph-part'),
    nickName = cms.string('TTJetsFullHad'),
    color    = cms.int32(41),
    xSec     = cms.double(133.62),
    ),
    
    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_HadronicMGDecays_8TeV-madgraph-part'),
    nickName = cms.string('TTJetsFullHad_LF'),
    color    = cms.int32(41),
    xSec     = cms.double(133.62),
    cut      = cms.string(ttjetsLF),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_HadronicMGDecays_8TeV-madgraph-part'),
    nickName = cms.string('TTJetsFullHad_C'),
    color    = cms.int32(44),
    xSec     = cms.double(133.62),
    cut      = cms.string(ttjetsC),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_HadronicMGDecays_8TeV-madgraph-part'),
    nickName = cms.string('TTJetsFullHad_B'),
    color    = cms.int32(46),
    xSec     = cms.double(133.62),
    cut      = cms.string(ttjetsB),
    ),
    

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_FullLeptMGDecays_8TeV-madgraph-part'),
    nickName = cms.string('TTJetsFullLept'),
    color    = cms.int32(41),
    xSec     = cms.double(24.56),
    ),
    
    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_FullLeptMGDecays_8TeV-madgraph-part'),
    nickName = cms.string('TTJetsFullLept_LF'),
    color    = cms.int32(41),
    xSec     = cms.double(24.56),
    cut      = cms.string(ttjetsLF),
    ),
    
    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_FullLeptMGDecays_8TeV-madgraph-part'),
    nickName = cms.string('TTJetsFullLept_C'),
    color    = cms.int32(44),
    xSec     = cms.double(24.56),
    cut      = cms.string(ttjetsC)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_FullLeptMGDecays_8TeV-madgraph-part'),
    nickName = cms.string('TTJetsFullLept_B'),
    color    = cms.int32(46),
    xSec     = cms.double(24.56),
    cut      = cms.string(ttjetsB)
    ),
    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'),
    nickName = cms.string('TTJetsSemiLept'),
    color    = cms.int32(41),
    xSec     = cms.double(75.82),
    ),
    
    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'),
    nickName = cms.string('TTJetsSemiLept_LF'),
    color    = cms.int32(41),
    xSec     = cms.double(75.82),
    cut      = cms.string(ttjetsLF),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'),
    nickName = cms.string('TTJetsSemiLept_C'),
    color    = cms.int32(44),
    xSec     = cms.double(75.82),
    cut      = cms.string(ttjetsC),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'),
    nickName = cms.string('TTJetsSemiLept_B'),
    color    = cms.int32(46),
    xSec     = cms.double(75.82),
    cut      = cms.string(ttjetsB),
    ),


    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('TtW'),
    color    = cms.int32(6),
    xSec     = cms.double(11.1)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('T_t-channel_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('Tt'),
    color    = cms.int32(6),
    xSec     = cms.double(56.4)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('T_s-channel_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('Ts'),
    color    = cms.int32(6),
    xSec     = cms.double(3.79)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('TbartW'),
    color    = cms.int32(6),
    xSec     = cms.double(11.1)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('Tbart'),
    color    = cms.int32(6),
    xSec     = cms.double(30.7)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('Tbars'),
    color    = cms.int32(6),
    xSec     = cms.double(1.76)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('WW_TuneZ2star_8TeV_pythia6_tauola'),
    nickName = cms.string('WW'),
    color    = cms.int32(4),
    xSec     = cms.double(56.75)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('WZ_TuneZ2star_8TeV_pythia6_tauola'),
    nickName = cms.string('WZ'),
    color    = cms.int32(4),
    xSec     = cms.double(33.85)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('ZZ_TuneZ2star_8TeV_pythia6_tauola'),
    nickName = cms.string('ZZ'),
    color    = cms.int32(4),
    xSec     = cms.double(8.297)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-110_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH110'),
    color    = cms.int32(2),
    xSec     = cms.double(0.04414)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-115_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH115'),
    color    = cms.int32(2),
    xSec     = cms.double(0.036375)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-120_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH120'),
    color    = cms.int32(2),
    xSec     = cms.double(0.0293327854)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH125'),
    color    = cms.int32(2),
    xSec     = cms.double(0.0229727058)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-130_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH130'),
    color    = cms.int32(2),
    xSec     = cms.double(0.017288657)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-135_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH135'),
    color    = cms.int32(2),
    xSec     = cms.double(0.01250888)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('DataZmm'),
    nickName = cms.string('DataZmm'),
    color    = cms.int32(1),
    xSec     = cms.double(-1)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('DataZee'),
    nickName = cms.string('DataZee'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut)
    ),
    
    
    ),



    plots      = cms.VPSet(

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(10),
    xHigh     = cms.double(250),
    nBins     = cms.int32(120),
    variable  = cms.string("vLepton_pt[0]"),
    xTitle    = cms.string("p_{T} lead e"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_leadElePt"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(10),
    xHigh     = cms.double(250),
    nBins     = cms.int32(120),
    variable  = cms.string("vLepton_pt[1]"),
    xTitle    = cms.string("p_{T} trail e"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_trailElePt"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(30),
    nBins     = cms.int32(30),
    variable  = cms.string("nPVs"),
    xTitle    = cms.string("vertices"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_vertex"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(200),
    nBins     = cms.int32(50),
    variable  = cms.string("METtype1p2corr.et"),
    xTitle    = cms.string("MET"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_met"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(400),
    nBins     = cms.int32(100),
    variable  = cms.string("V.mass"),
    xTitle    = cms.string("di-e mass"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_diEleMass"),
    cut       = cms.string(common),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(-3),
    xHigh     = cms.double(3),
    nBins     = cms.int32(60),
    variable  = cms.string("vLepton_eta[0]"),
    xTitle    = cms.string("#eta lead e"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_leadEleEta"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    
    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(-3),
    xHigh     = cms.double(3),
    nBins     = cms.int32(60),
    variable  = cms.string("vLepton_eta[1]"),
    xTitle    = cms.string("#eta trail e"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_trailEleEta"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),


    ########################## JETS ###############################
    
    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(10),
    nBins     = cms.int32(10),
    variable  = cms.string("numJets30"),
    xTitle    = cms.string("jet mult p_{T}>30 GeV"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_numJet30"),
    cut       = cms.string(common),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(10),
    nBins     = cms.int32(10),
    variable  = cms.string("numJets30bTag"),
    xTitle    = cms.string("b-jet mult"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_numBJet30"),
    cut       = cms.string(common),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(1.1),
    nBins     = cms.int32(55),
    variable  = cms.string("csv1"),
    xTitle    = cms.string("CSV jet 1"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_csv1"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(1.1),
    nBins     = cms.int32(55),
    variable  = cms.string("csv2"),
    xTitle    = cms.string("CSV jet 2"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_csv2"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(-5.5),
    xHigh     = cms.double(5.5),
    nBins     = cms.int32(55),
    variable  = cms.string("eta1"),
    xTitle    = cms.string("#eta jet 1"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_eta1"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(-5.5),
    xHigh     = cms.double(5.5),
    nBins     = cms.int32(55),
    variable  = cms.string("eta2"),
    xTitle    = cms.string("#eta jet 2"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_eta2"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(30),
    xHigh     = cms.double(300),
    nBins     = cms.int32(90),
    variable  = cms.string("pt1"),
    xTitle    = cms.string("p_{T} jet 1"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_pt1"),
    cut       = cms.string(common),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(30),
    xHigh     = cms.double(300),
    nBins     = cms.int32(90),
    variable  = cms.string("pt2"),
    xTitle    = cms.string("p_{T} jet 2"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_pt2"),
    cut       = cms.string(common),
    logy      = cms.int32(1),
    ),


    ########################## 2j2b ###########################


    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(10),
    nBins     = cms.int32(10),
    variable  = cms.string("numJets30"),
    xTitle    = cms.string("jet mult p_{T}>30 GeV"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_numJet30_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(10),
    nBins     = cms.int32(10),
    variable  = cms.string("numJets30bTag"),
    xTitle    = cms.string("b-jet mult"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_numBJet30_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(1.1),
    nBins     = cms.int32(55),
    variable  = cms.string("csv1"),
    xTitle    = cms.string("CSV jet 1"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_csv1_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(1.1),
    nBins     = cms.int32(55),
    variable  = cms.string("csv2"),
    xTitle    = cms.string("CSV jet 2"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_csv3_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(-5.5),
    xHigh     = cms.double(5.5),
    nBins     = cms.int32(55),
    variable  = cms.string("eta1"),
    xTitle    = cms.string("#eta jet 1"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_eta1_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(-5.5),
    xHigh     = cms.double(5.5),
    nBins     = cms.int32(55),
    variable  = cms.string("eta2"),
    xTitle    = cms.string("#eta jet 2"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_eta2_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(30),
    xHigh     = cms.double(390),
    nBins     = cms.int32(60),
    variable  = cms.string("pt1"),
    xTitle    = cms.string("p_{T} jet 1"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_pt1_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(30),
    xHigh     = cms.double(390),
    nBins     = cms.int32(60),
    variable  = cms.string("pt2"),
    xTitle    = cms.string("p_{T} jet 2"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_pt2_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(8),
    nBins     = cms.int32(8),
    variable  = cms.string("(lheNj - 2 - (abs(genTop.wdau1id)<6 + abs(genTop.wdau2id)<6 + abs(genTbar.wdau1id)<6 + abs(genTbar.wdau2id)<6 ) )"),
    xTitle    = cms.string("LHE N_{j}"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_lheNj_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(3),
    nBins     = cms.int32(3),
    variable  = cms.string("nLF"),
    xTitle    = cms.string("LF jets"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_nLF_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(3),
    nBins     = cms.int32(3),
    variable  = cms.string("nLFTop"),
    xTitle    = cms.string("LF jets from W"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_nLFTop_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(3),
    nBins     = cms.int32(3),
    variable  = cms.string("nC"),
    xTitle    = cms.string("charm jets"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_nC_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),
    
    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(3),
    nBins     = cms.int32(3),
    variable  = cms.string("nCTop"),
    xTitle    = cms.string("charm jets from W"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_nCTop_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),
    
    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(3),
    nBins     = cms.int32(3),
    variable  = cms.string("nB"),
    xTitle    = cms.string("bottom jets"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_nB_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    logy      = cms.int32(0),
    normalize = cms.int32(1),
    ),
    
    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(3),
    nBins     = cms.int32(3),
    variable  = cms.string("nBTop"),
    xTitle    = cms.string("bottom jets from top"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_nBTop_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(-10),
    xHigh     = cms.double(25),
    nBins     = cms.int32(35),
    variable  = cms.string("flavor1"),
    xTitle    = cms.string("flavor jet 1"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_flavor1_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(-10),
    xHigh     = cms.double(25),
    nBins     = cms.int32(35),
    variable  = cms.string("flavor2"),
    xTitle    = cms.string("flavor jet 2"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_flavor2_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(2),
    nBins     = cms.int32(2),
    variable  = cms.string("(topB1 || topW1)"),
    xTitle    = cms.string("jet 1 from top decay?"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_top1_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(2),
    nBins     = cms.int32(2),
    variable  = cms.string("(topB2 || topW2)"),
    xTitle    = cms.string("jet 2 from top decay?"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_top2_2j2b"),
    cut       = cms.string(common+" && numJets30bTag==2"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),


    ########################## 3j3b ###########################


    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(10),
    nBins     = cms.int32(10),
    variable  = cms.string("numJets30"),
    xTitle    = cms.string("jet mult p_{T}>30 GeV"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_numJet30_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(10),
    nBins     = cms.int32(10),
    variable  = cms.string("numJets30bTag"),
    xTitle    = cms.string("b-jet mult"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_numBJet30_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(1.1),
    nBins     = cms.int32(22),
    variable  = cms.string("csv1"),
    xTitle    = cms.string("CSV jet 1"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_csv1_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(1.1),
    nBins     = cms.int32(22),
    variable  = cms.string("csv2"),
    xTitle    = cms.string("CSV jet 2"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_csv2_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(1.1),
    nBins     = cms.int32(22),
    variable  = cms.string("csv3"),
    xTitle    = cms.string("CSV jet 3"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_csv3_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-5.5),
    xHigh     = cms.double(5.5),
    nBins     = cms.int32(22),
    variable  = cms.string("eta1"),
    xTitle    = cms.string("#eta jet 1"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_eta1_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-5.5),
    xHigh     = cms.double(5.5),
    nBins     = cms.int32(22),
    variable  = cms.string("eta2"),
    xTitle    = cms.string("#eta jet 2"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_eta2_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-5.5),
    xHigh     = cms.double(5.5),
    nBins     = cms.int32(22),
    variable  = cms.string("eta3"),
    xTitle    = cms.string("#eta jet 3"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_eta2_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(30),
    xHigh     = cms.double(390),
    nBins     = cms.int32(60),
    variable  = cms.string("pt1"),
    xTitle    = cms.string("p_{T} jet 1"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_pt1_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(30),
    xHigh     = cms.double(390),
    nBins     = cms.int32(60),
    variable  = cms.string("pt2"),
    xTitle    = cms.string("p_{T} jet 2"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_pt2_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(30),
    xHigh     = cms.double(390),
    nBins     = cms.int32(60),
    variable  = cms.string("pt3"),
    xTitle    = cms.string("p_{T} jet 3"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_pt3_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(8),
    nBins     = cms.int32(8),
    variable  = cms.string("(lheNj - 2 - (genTop.bpt>0.01)*(abs(genTop.wdau1id)<6 + abs(genTop.wdau2id)<6) - (genTbar.bpt>0.01)*(abs(genTbar.wdau1id)<6 + abs(genTbar.wdau2id)<6 ) )"),
    xTitle    = cms.string("LHE N_{j}"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_lheNj_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),
    
    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(6),
    nBins     = cms.int32(6),
    variable  = cms.string("nLF"),
    xTitle    = cms.string("LF jets"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_nLF_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(6),
    nBins     = cms.int32(6),
    variable  = cms.string("nLFTop"),
    xTitle    = cms.string("LF jets from W"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_nLFTop_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(6),
    nBins     = cms.int32(6),
    variable  = cms.string("nC"),
    xTitle    = cms.string("charm jets"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_nC_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),
    
    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(6),
    nBins     = cms.int32(6),
    variable  = cms.string("nCTop"),
    xTitle    = cms.string("charm jets from W"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_nCTop_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),
    
    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(6),
    nBins     = cms.int32(6),
    variable  = cms.string("nB"),
    xTitle    = cms.string("bottom jets"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_nB_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),
    
    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(6),
    nBins     = cms.int32(6),
    variable  = cms.string("nBTop"),
    xTitle    = cms.string("bottom jets from top"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_nBTop_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-10),
    xHigh     = cms.double(25),
    nBins     = cms.int32(35),
    variable  = cms.string("flavor1"),
    xTitle    = cms.string("flavor jet 1"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_flavor1_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-10),
    xHigh     = cms.double(25),
    nBins     = cms.int32(35),
    variable  = cms.string("flavor2"),
    xTitle    = cms.string("flavor jet 2"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_flavor2_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-10),
    xHigh     = cms.double(25),
    nBins     = cms.int32(35),
    variable  = cms.string("flavor3"),
    xTitle    = cms.string("flavor jet 3"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_flavor3_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(2),
    nBins     = cms.int32(2),
    variable  = cms.string("(topB1 || topW1)"),
    xTitle    = cms.string("jet 1 from top decay?"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_top1_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(2),
    nBins     = cms.int32(2),
    variable  = cms.string("(topB2 || topW2)"),
    xTitle    = cms.string("jet 2 from top decay?"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_top2_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),


    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(2),
    nBins     = cms.int32(2),
    variable  = cms.string("(topB3 || topW3)"),
    xTitle    = cms.string("jet 3 from top decay?"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_top3_3j3b"),
    cut       = cms.string(common+" && numJets30bTag>=3 && pt3>30"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
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

    )

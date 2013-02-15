import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms


VType = "_VType2"

mmBasic      = "(V.pt>0 && V.Mt>40 && H.pt>0 && Vtype==2 && H.HiggsFlag==1)"
mmKin        = "(vLepton_pt[0]>30 && abs(vLepton_eta[0])<2.1)"
mmId         = "(vLepton_pfCorrIso[0]<0.10)"
mmJetBasic   = "(numJets30>=3 && pt1>30 && pt2>30 && pt4>30 && numJets30bTag>=2)"
common       =  mmJetBasic + " && " + mmBasic+" && "+mmKin+" && "+mmId 



dataCut      = "((EVENT.json == 1 || EVENT.run < 196532) && (triggerFlags[14] || triggerFlags[15]  || triggerFlags[21]  || triggerFlags[22] || triggerFlags[23]))"

ttjetsLF = \
         "nSimBs==2 && nC-nCTop==0"
ttjetsC  = \
        "nSimBs==2 && nC-nCTop>0"
ttjetsB  = \
        "nSimBs>2"

zjetsLF = \
        "nSimBs==0 && nC==0"
zjetsC  = \
       "nSimBs==0 && nC>0"
zjetsB  = \
       "nSimBs>0"

xsecTT_FH = 106.9
xsecTT_SL = 103.0
xsecTT_FL = 24.8

process = cms.Process("FWLitePlots")

process.fwliteInput = cms.PSet(

    pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt"+VType),
    #pathToFile    = cms.string("/scratch/bianchi/HBB_EDMNtuple/All.H.DiJetPt/"),
    ordering      = cms.string("DiJetPt_"),
    lumi          = cms.double(12.1),
    debug         = cms.bool(False),

    samples       = cms.VPSet(

    cms.PSet(
    skip     = cms.bool(False),
    name     = cms.string('DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball'+VType),
    nickName = cms.string('DYJets'),
    color    = cms.int32(19),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('WJetsToLNu_PtW-70To100_TuneZ2star_8TeV-madgraph'+VType),
    nickName = cms.string('WJetsPt70100'),
    color    = cms.int32(29),
    xSec     = cms.double(37509.0*0.01793),  ## USE xsection from DY+jets!!!!
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('WJetsToLNu_PtW-100_TuneZ2star_8TeV-madgraph'+VType),
    nickName = cms.string('WJetsPt100'),
    color    = cms.int32(29),
    xSec     = cms.double(37509.0*0.0118),   ## USE xsection from DY+jets!!!!
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_HadronicMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullHad'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_FH),
    ),
    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_HadronicMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullHad_LF'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_FH),
    cut      = cms.string(ttjetsLF),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_HadronicMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullHad_C'),
    color    = cms.int32(44),
    xSec     = cms.double(xsecTT_FH),
    cut      = cms.string(ttjetsC),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_HadronicMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullHad_B'),
    color    = cms.int32(46),
    xSec     = cms.double(xsecTT_FH),
    cut      = cms.string(ttjetsB),
    ),
    

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_FullLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullLept'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_FL),
    ),
    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_FullLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullLept_LF'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_FL),
    cut      = cms.string(ttjetsLF),
    ),
    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_FullLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullLept_C'),
    color    = cms.int32(44),
    xSec     = cms.double(xsecTT_FL),
    cut      = cms.string(ttjetsC)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_FullLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullLept_B'),
    color    = cms.int32(46),
    xSec     = cms.double(xsecTT_FL),
    cut      = cms.string(ttjetsB)
    ),
    
    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsSemiLept'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_SL),
    ),
    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsSemiLept_LF'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_SL),
    cut      = cms.string(ttjetsLF),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsSemiLept_C'),
    color    = cms.int32(44),
    xSec     = cms.double(xsecTT_SL),
    cut      = cms.string(ttjetsC),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsSemiLept_B'),
    color    = cms.int32(46),
    xSec     = cms.double(xsecTT_SL),
    cut      = cms.string(ttjetsB),
    ),

    
    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('TtW'),
    color    = cms.int32(6),
    xSec     = cms.double(11.1)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('T_t-channel_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('Tt'),
    color    = cms.int32(6),
    xSec     = cms.double(56.4)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('T_s-channel_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('Ts'),
    color    = cms.int32(6),
    xSec     = cms.double(3.79)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('TbartW'),
    color    = cms.int32(6),
    xSec     = cms.double(11.1)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('Tbart'),
    color    = cms.int32(6),
    xSec     = cms.double(30.7)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('Tbars'),
    color    = cms.int32(6),
    xSec     = cms.double(1.76)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('WW_TuneZ2star_8TeV_pythia6_tauola'+VType),
    nickName = cms.string('WW'),
    color    = cms.int32(4),
    xSec     = cms.double(56.75)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('WZ_TuneZ2star_8TeV_pythia6_tauola'+VType),
    nickName = cms.string('WZ'),
    color    = cms.int32(4),
    xSec     = cms.double(33.85)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('ZZ_TuneZ2star_8TeV_pythia6_tauola'+VType),
    nickName = cms.string('ZZ'),
    color    = cms.int32(4),
    xSec     = cms.double(8.297)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-110_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('ZH110'),
    color    = cms.int32(2),
    xSec     = cms.double(0.04414)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-115_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('ZH115'),
    color    = cms.int32(2),
    xSec     = cms.double(0.036375)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-120_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('ZH120'),
    color    = cms.int32(2),
    xSec     = cms.double(0.0293327854)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('ZH125'),
    color    = cms.int32(2),
    xSec     = cms.double(0.0229727058)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-130_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('ZH130'),
    color    = cms.int32(2),
    xSec     = cms.double(0.017288657)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-135_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('ZH135'),
    color    = cms.int32(2),
    xSec     = cms.double(0.01250888)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WH_WToLNu_HToBB_M-110_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('WH110'),
    color    = cms.int32(2),
    xSec     = cms.double(0.78864),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WH_WToLNu_HToBB_M-115_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('WH115'),
    color    = cms.int32(2),
    xSec     = cms.double(0.644299499),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WH_WToLNu_HToBB_M-120_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('WH120'),
    color    = cms.int32(2),
    xSec     = cms.double(0.5161968)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('WH125'),
    color    = cms.int32(2),
    xSec     = cms.double(0.4019381),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WH_WToLNu_HToBB_M-130_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('WH130'),
    color    = cms.int32(2),
    xSec     = cms.double(0.3010930)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WH_WToLNu_HToBB_M-135_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('WH135'),
    color    = cms.int32(2),
    xSec     = cms.double(0.230628)
    ),


    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('SingleMuRun2012AAug06EdmV42'+VType),
    nickName = cms.string('DataSingleMu_1'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('SingleMuRun2012AJul13EdmV42'+VType),
    nickName = cms.string('DataSingleMu_2'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('SingleMuRun2012BJul13EdmV42'+VType),
    nickName = cms.string('DataSingleMu_3'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

     cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('SingleMuRun2012CAug24RerecoEdmV42'+VType),
    nickName = cms.string('DataSingleMu_4'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('SingleMuRun2012CPromptv2EdmV42'+VType),
    nickName = cms.string('DataSingleMu_5'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('SingleMuRun2012CPromptV2TopUpEdmV42'+VType),
    nickName = cms.string('DataSingleMu_6'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('SingleMuRun2012'),
    nickName = cms.string('DataSingleMu'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    
    ),



    plots      = cms.VPSet(

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(10),
    xHigh     = cms.double(250),
    nBins     = cms.int32(120),
    variable  = cms.string("vLepton_pt[0]"),
    xTitle    = cms.string("p_{T} lead #mu"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_leadMuPt"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(0.1),
    nBins     = cms.int32(40),
    variable  = cms.string("vLepton_pfCorrIso[0]"),
    xTitle    = cms.string("iso"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_iso"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),
    
    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(30),
    nBins     = cms.int32(30),
    variable  = cms.string("nPVs"),
    xTitle    = cms.string("vertices"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_vertex"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(200),
    nBins     = cms.int32(50),
    variable  = cms.string("METtype1p2corr.et"),
    xTitle    = cms.string("MET"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_met"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-3),
    xHigh     = cms.double(3),
    nBins     = cms.int32(60),
    variable  = cms.string("vLepton_eta[0]"),
    xTitle    = cms.string("#eta lead #mu"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_leadMuEta"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),


    ########################## JETS ###############################
    
    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(10),
    nBins     = cms.int32(10),
    variable  = cms.string("numJets30"),
    xTitle    = cms.string("jet mult p_{T}>30 GeV"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_numJet30"),
    cut       = cms.string(common),
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
    histoName = cms.string("VType2_numBJet30"),
    cut       = cms.string(common),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(1.1),
    nBins     = cms.int32(55),
    variable  = cms.string("csv1"),
    xTitle    = cms.string("CSV jet 1"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_csv1"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(1.1),
    nBins     = cms.int32(55),
    variable  = cms.string("csv2"),
    xTitle    = cms.string("CSV jet 2"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_csv2"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(1.1),
    nBins     = cms.int32(55),
    variable  = cms.string("csv3"),
    xTitle    = cms.string("CSV jet 3"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_csv3"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-5.5),
    xHigh     = cms.double(5.5),
    nBins     = cms.int32(55),
    variable  = cms.string("eta1"),
    xTitle    = cms.string("#eta jet 1"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_eta1"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-5.5),
    xHigh     = cms.double(5.5),
    nBins     = cms.int32(55),
    variable  = cms.string("eta2"),
    xTitle    = cms.string("#eta jet 2"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_eta2"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-5.5),
    xHigh     = cms.double(5.5),
    nBins     = cms.int32(55),
    variable  = cms.string("eta3"),
    xTitle    = cms.string("#eta jet 3"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_eta3"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(30),
    xHigh     = cms.double(300),
    nBins     = cms.int32(90),
    variable  = cms.string("pt1"),
    xTitle    = cms.string("p_{T} jet 1"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_pt1"),
    cut       = cms.string(common),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(30),
    xHigh     = cms.double(300),
    nBins     = cms.int32(90),
    variable  = cms.string("pt2"),
    xTitle    = cms.string("p_{T} jet 2"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_pt2"),
    cut       = cms.string(common),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(30),
    xHigh     = cms.double(300),
    nBins     = cms.int32(90),
    variable  = cms.string("pt3"),
    xTitle    = cms.string("p_{T} jet 3"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_pt3"),
    cut       = cms.string(common),
    logy      = cms.int32(1),
    ),


    ########################## 4j3b ###########################


    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(10),
    nBins     = cms.int32(10),
    variable  = cms.string("numJets30"),
    xTitle    = cms.string("jet mult p_{T}>30 GeV"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_numJet30_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_numBJet30_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_csv1_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_csv2_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(1.1),
    nBins     = cms.int32(55),
    variable  = cms.string("csv3"),
    xTitle    = cms.string("CSV jet 3"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_csv3_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(1.1),
    nBins     = cms.int32(55),
    variable  = cms.string("csv4"),
    xTitle    = cms.string("CSV jet 4"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_csv4_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_eta1_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_eta2_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(-5.5),
    xHigh     = cms.double(5.5),
    nBins     = cms.int32(55),
    variable  = cms.string("eta3"),
    xTitle    = cms.string("#eta jet 3"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_eta3_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(-5.5),
    xHigh     = cms.double(5.5),
    nBins     = cms.int32(55),
    variable  = cms.string("eta4"),
    xTitle    = cms.string("#eta jet 4"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_eta4_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_pt1_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_pt2_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(30),
    xHigh     = cms.double(390),
    nBins     = cms.int32(60),
    variable  = cms.string("pt3"),
    xTitle    = cms.string("p_{T} jet 3"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_pt3_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(30),
    xHigh     = cms.double(390),
    nBins     = cms.int32(60),
    variable  = cms.string("pt4"),
    xTitle    = cms.string("p_{T} jet 4"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_pt4_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_lheNj_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_nLF_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_nLFTop_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_nC_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_nCTop_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_nB_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_nBTop_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_flavor1_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_flavor2_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(-10),
    xHigh     = cms.double(25),
    nBins     = cms.int32(35),
    variable  = cms.string("flavor3"),
    xTitle    = cms.string("flavor jet 3"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_flavor3_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(-10),
    xHigh     = cms.double(25),
    nBins     = cms.int32(35),
    variable  = cms.string("flavor4"),
    xTitle    = cms.string("flavor jet 4"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_flavor4_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_top1_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
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
    histoName = cms.string("VType2_top2_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(2),
    nBins     = cms.int32(2),
    variable  = cms.string("(topB3 || topW3)"),
    xTitle    = cms.string("jet 3 from top decay?"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_top3_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(2),
    nBins     = cms.int32(2),
    variable  = cms.string("(topB4 || topW4)"),
    xTitle    = cms.string("jet 4 from top decay?"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_top4_4j3b"),
    cut       = cms.string(common+" && numJets30==4 && pt4>30 && numJets30bTag==3"),
    normalize = cms.int32(1),
    logy      = cms.int32(0),
    ),


    ########################## 4j4b ###########################

    ########################## 5j3b ###########################

    ########################## 5j4b ###########################

    ########################## 5j5b ###########################

    ########################## 6j2b ###########################

    ########################## 6j3b ###########################

    ########################## 6j4b ###########################




    
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

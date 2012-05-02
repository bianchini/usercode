import FWCore.ParameterSet.Config as cms

process = cms.Process("MUTAUANA")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

#from Configuration.PyReleaseValidation.autoCond import autoCond
#process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

runOnMC     = True
doSVFitReco = False

if runOnMC:
    print "Running on MC"
else:
    print "Running on Data"
        
if runOnMC:
    process.GlobalTag.globaltag = cms.string('START52_V7::All')

else:
    process.GlobalTag.globaltag = cms.string('GR_R_52_V7::All')
    
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 50
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    #'file:./patTuples_MuTauStream.root'

    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012-Summer12-50_DYJets-MuTau-madgraph-50_skim/91f7bad7a7ab94b1978a2d9f0d8f649d/patTuples_MuTauStream_100_1_I8e.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012-Summer12-50_DYJets-MuTau-madgraph-50_skim/91f7bad7a7ab94b1978a2d9f0d8f649d/patTuples_MuTauStream_101_1_P8u.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012-Summer12-50_DYJets-MuTau-madgraph-50_skim/91f7bad7a7ab94b1978a2d9f0d8f649d/patTuples_MuTauStream_102_1_sG6.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012-Summer12-50_DYJets-MuTau-madgraph-50_skim/91f7bad7a7ab94b1978a2d9f0d8f649d/patTuples_MuTauStream_103_1_Li8.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012-Summer12-50_DYJets-MuTau-madgraph-50_skim/91f7bad7a7ab94b1978a2d9f0d8f649d/patTuples_MuTauStream_104_1_Fgr.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012-Summer12-50_DYJets-MuTau-madgraph-50_skim/91f7bad7a7ab94b1978a2d9f0d8f649d/patTuples_MuTauStream_105_1_zFG.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012-Summer12-50_DYJets-MuTau-madgraph-50_skim/91f7bad7a7ab94b1978a2d9f0d8f649d/patTuples_MuTauStream_106_1_zvC.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012-Summer12-50_DYJets-MuTau-madgraph-50_skim/91f7bad7a7ab94b1978a2d9f0d8f649d/patTuples_MuTauStream_107_1_kMf.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012-Summer12-50_DYJets-MuTau-madgraph-50_skim/91f7bad7a7ab94b1978a2d9f0d8f649d/patTuples_MuTauStream_108_1_Z4o.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012-Summer12-50_DYJets-MuTau-madgraph-50_skim/91f7bad7a7ab94b1978a2d9f0d8f649d/patTuples_MuTauStream_10_1_Yo7.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012-Summer12-50_DYJets-MuTau-madgraph-50_skim/91f7bad7a7ab94b1978a2d9f0d8f649d/patTuples_MuTauStream_110_1_ggd.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012-Summer12-50_DYJets-MuTau-madgraph-50_skim/91f7bad7a7ab94b1978a2d9f0d8f649d/patTuples_MuTauStream_111_1_SeT.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012-Summer12-50_DYJets-MuTau-madgraph-50_skim/91f7bad7a7ab94b1978a2d9f0d8f649d/patTuples_MuTauStream_112_1_1lY.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012-Summer12-50_DYJets-MuTau-madgraph-50_skim/91f7bad7a7ab94b1978a2d9f0d8f649d/patTuples_MuTauStream_113_1_6vR.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012-Summer12-50_DYJets-MuTau-madgraph-50_skim/91f7bad7a7ab94b1978a2d9f0d8f649d/patTuples_MuTauStream_114_1_NC0.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012-Summer12-50_DYJets-MuTau-madgraph-50_skim/91f7bad7a7ab94b1978a2d9f0d8f649d/patTuples_MuTauStream_115_1_b9F.root',

    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_100_1_V3f.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_101_1_UFq.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_10_1_Ath.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_11_1_fGK.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_12_1_veW.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_13_1_vN3.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_14_1_jeL.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_15_1_4gS.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_16_1_3VO.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_17_1_vIm.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_18_1_U3T.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_19_1_cgu.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_1_1_Y0i.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_20_1_bt8.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_21_1_VpR.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_22_1_8OE.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_23_1_nr4.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_24_1_u2T.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_25_1_whT.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_26_1_Ion.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_27_1_6oe.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_28_1_isB.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_29_1_Y73.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_2_1_De8.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_30_1_1YF.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_31_1_Vbj.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_32_1_Iih.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_33_1_49D.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_34_1_23M.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_35_1_Tvg.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_36_1_U65.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_37_1_iMB.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_38_1_rmw.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_39_1_Snz.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBFHToTauTau_M-120_8TeV-pythia6/MuTauStream-30Mar2012-Summer12-50_VBFH120-MuTau-pythia-50_skim/867500062c684278db589e0b621afc8f/patTuples_MuTauStream_3_1_Euz.root',


    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_100_2_YSa.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_101_2_tgJ.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_102_2_9zL.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_103_1_FvJ.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_104_2_enx.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_105_2_9A1.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_107_2_WbE.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_108_1_8Zu.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_109_1_I76.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_10_2_QML.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_11_2_fwd.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_12_1_soI.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_13_1_jtq.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_15_2_MLo.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_16_2_1lH.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_17_2_s3G.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_18_1_BXJ.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_19_1_e67.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_1_1_IS7.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_20_1_zIb.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_21_2_ksw.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_22_2_5yr.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_23_1_zFn.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_24_1_f2w.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_25_2_U3z.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_26_1_Sjw.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_29_2_uvn.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_2_1_ers.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_30_1_Ydo.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_31_2_ZYR.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_32_1_SaM.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_33_2_2rm.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_34_1_QTQ.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_35_2_eMb.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_36_2_9XZ.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_37_2_Axu.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_38_2_4wf.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_39_1_JAj.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_3_1_JDL.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_40_2_lIR.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_41_1_N7C.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_42_1_EAr.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_43_2_KiP.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_44_2_31H.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_45_1_NG3.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_46_1_uUg.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_47_2_cYC.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_48_2_kB9.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_49_2_wam.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_4_1_hJO.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_50_2_B5m.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_51_2_iPo.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_52_2_VXq.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_53_1_kRx.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_54_1_D4p.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_55_2_HAr.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_56_2_wFL.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_57_1_Ptb.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_58_1_MMR.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_59_1_HNg.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_5_2_apl.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_60_2_dRE.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_61_2_HVf.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_62_2_8NP.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_63_1_f8j.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_64_1_S3w.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_65_1_ojW.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_66_1_dvk.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_68_1_sQj.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_69_1_jab.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_6_1_fHc.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_70_1_Pks.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_71_1_cdo.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_72_1_ynK.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_73_2_BCr.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_74_2_Rsu.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_75_1_EUu.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_76_1_YQm.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_78_2_gxo.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_7_1_Min.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_80_1_XGF.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_81_2_MR7.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_82_1_yQE.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_83_1_F9e.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_84_2_RpW.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_85_2_Pom.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_87_1_HRi.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_88_2_1lf.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_89_2_kaS.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_8_1_IWq.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_90_1_eBo.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_92_2_OGR.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_93_1_d2F.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_94_1_DZ5.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_95_1_AWJ.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_96_2_XEj.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_97_1_JOL.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_98_1_bHV.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_99_1_SAu.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-madgraph-Tarball-v2_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_9_1_azg.root'

    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-Tarball_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_10_1_gOR.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-Tarball_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_11_1_b5k.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-Tarball_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_12_1_zot.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-Tarball_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_13_1_t63.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-Tarball_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_1_1_3xN.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-Tarball_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_2_1_o8H.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-Tarball_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_3_1_QyJ.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-Tarball_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_4_1_iPl.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-Tarball_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_5_1_Wkg.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-Tarball_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_6_1_SzB.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-Tarball_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_7_1_zqZ.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-Tarball_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_8_1_7AW.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/MuTauStream-30Mar2012_DYJets-MuTau-Tarball_skim/788ef0e87cfaf0843a8d9198e4ce2a38/patTuples_MuTauStream_9_1_FHW.root',
    )
    )

process.allEventsFilter = cms.EDFilter(
    "AllEventsFilter"
    )

###################################################################################
process.rescaledMET = cms.EDProducer(
    "MEtRescalerProducer",
    metTag         = cms.InputTag("metRecoilCorrector",  "N"),
    jetTag         = cms.InputTag("selectedPatJets"),
    electronTag    = cms.InputTag(""),
    muonTag        = cms.InputTag("muPtEtaIDIso"),
    tauTag         = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
    unClusterShift = cms.double(0.10),
    tauShift       = cms.vdouble(0.03,0.03),
    muonShift      = cms.vdouble(0.01,0.01),
    electronShift  = cms.vdouble(0.01,0.025),
    jetThreshold   = cms.double(10),
    numOfSigmas    = cms.double(1.0),
    verbose = cms.bool(False)
    )

process.rescaledMETjet = process.rescaledMET.clone(
    unClusterShift = cms.double(0.10),
    tauShift       = cms.vdouble(0.0,0.0),
    muonShift      = cms.vdouble(0.0,0.0),
    electronShift  = cms.vdouble(0.0,0.0),
    )
process.rescaledMETtau = process.rescaledMET.clone(
    unClusterShift = cms.double(0.0),
    tauShift       = cms.vdouble(0.03,0.03),
    muonShift      = cms.vdouble(0.0,0.0),
    electronShift  = cms.vdouble(0.0,0.0),
    )
process.rescaledMETmuon = process.rescaledMET.clone(
    unClusterShift = cms.double(0.0),
    tauShift       = cms.vdouble(0.0,0.0),
    muonShift      = cms.vdouble(0.01,0.01),
    electronShift  = cms.vdouble(0.0,0.0),
    )

process.rescaledTaus = cms.EDProducer(
    "TauRescalerProducer",
    inputCollection = cms.InputTag("tauPtEtaIDAgMuAgElecIsoPtRel"),
    shift           = cms.vdouble(0.03,0.03),
    numOfSigmas     = cms.double(1.0),
    #verbose         = cms.bool(True)
    )
process.rescaledMuons = cms.EDProducer(
    "MuonRescalerProducer",
    inputCollection = cms.InputTag("muPtEtaIDIsoPtRel"),
    shift           = cms.vdouble(0.01,0.01),
    numOfSigmas     = cms.double(1.0),
    #verbose         = cms.bool(True)
    )
process.rescaledMuonsRel = cms.EDProducer(
    "MuonRescalerProducer",
    inputCollection = cms.InputTag("muPtEtaRelID"),
    shift           = cms.vdouble(0.01,0.01),
    numOfSigmas     = cms.double(1.0),
    #verbose         = cms.bool(True)
    )

process.rescaledObjects = cms.Sequence(
    process.rescaledMETjet+
    process.rescaledMETtau+
    process.rescaledMETmuon+
    process.rescaledTaus+
    process.rescaledMuons+
    process.rescaledMuonsRel
    )

###################################################################################

process.metRecoilCorrector = cms.EDProducer(
    "MEtRecoilCorrectorProducer",
    genParticleTag      = cms.InputTag("genParticles"),
    jetTag              = cms.InputTag("selectedPatJets"),
    metTag              = cms.InputTag("patMETsPF"),
    electronTag         = cms.InputTag(""),
    muonTag             = cms.InputTag("muPtEtaIDIso"),
    tauTag              = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
    inputFileNamezmm42X = cms.FileInPath("Bianchi/Utilities/data/recoilv4/RecoilCorrector_v4/recoilfits/recoilfit_zmm42X_njet.root"),
    inputFileNamedatamm = cms.FileInPath("Bianchi/Utilities/data/recoilv4/RecoilCorrector_v4/recoilfits/recoilfit_datamm_njet.root"),
    inputFileNamewjets  = cms.FileInPath("Bianchi/Utilities/data/recoilv4/RecoilCorrector_v4/recoilfits/recoilfit_wjets_njet.root"),
    inputFileNamezjets  = cms.FileInPath("Bianchi/Utilities/data/recoilv4/RecoilCorrector_v4/recoilfits/recoilfit_zjets_ltau_njet.root"),
    inputFileNamehiggs  = cms.FileInPath("Bianchi/Utilities/data/recoilv4/RecoilCorrector_v4/recoilfits/recoilfit_higgs_njet.root"),
    numOfSigmas         = cms.double(1.0),
    minJetPt            = cms.double(30.0),
    verbose             = cms.bool(False),
    isMC                = cms.bool(runOnMC),
    )


###################################################################################

process.load("Bianchi.Utilities.diTausReconstruction_cff")
process.diTau = process.allMuTauPairs.clone()
process.diTau.srcLeg1  = cms.InputTag("muPtEtaIDIso")
process.diTau.srcLeg2  = cms.InputTag("tauPtEtaIDAgMuAgElecIso")
process.diTau.srcMET   = cms.InputTag("metRecoilCorrector",  "N")
process.diTau.dRmin12  = cms.double(0.5)
process.diTau.doSVreco = cms.bool(doSVFitReco)

if not runOnMC:
    process.diTau.srcGenParticles = ""
        
process.selectedDiTau = cms.EDFilter(
    "MuTauPairSelector",
    src = cms.InputTag("diTau"),
    cut = cms.string("dR12>0.5")
    )
process.selectedDiTauCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("selectedDiTau"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

#######################################################################

process.diTauJetUp =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                          srcLeg1 = cms.InputTag("muPtEtaIDIso"),
                                          srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                          srcMET  = cms.InputTag("rescaledMETjet",  "UNNNU")
                                          )
process.selectedDiTauJetUp = process.selectedDiTau.clone(src = cms.InputTag("diTauJetUp") )
process.selectedDiTauJetUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauJetUp"))

process.diTauJetDown =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                            srcLeg1 = cms.InputTag("muPtEtaIDIso"),
                                            srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                            srcMET  = cms.InputTag("rescaledMETjet",  "DNNND")
                                            )
process.selectedDiTauJetDown = process.selectedDiTau.clone(src = cms.InputTag("diTauJetDown") )
process.selectedDiTauJetDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauJetDown"))


process.diTauMEtResponseUp =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                          srcLeg1 = cms.InputTag("muPtEtaIDIso"),
                                          srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                          srcMET  = cms.InputTag("metRecoilCorrector",  "ResponseU")
                                          )
process.selectedDiTauMEtResponseUp = process.selectedDiTau.clone(src = cms.InputTag("diTauMEtResponseUp") )
process.selectedDiTauMEtResponseUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauMEtResponseUp"))

process.diTauMEtResponseDown =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                            srcLeg1 = cms.InputTag("muPtEtaIDIso"),
                                            srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                            srcMET  = cms.InputTag("metRecoilCorrector",  "ResponseD")
                                            )
process.selectedDiTauMEtResponseDown = process.selectedDiTau.clone(src = cms.InputTag("diTauMEtResponseDown") )
process.selectedDiTauMEtResponseDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauMEtResponseDown"))


process.diTauMEtResolutionUp =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                          srcLeg1 = cms.InputTag("muPtEtaIDIso"),
                                          srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                          srcMET  = cms.InputTag("metRecoilCorrector",  "ResolutionU")
                                          )
process.selectedDiTauMEtResolutionUp = process.selectedDiTau.clone(src = cms.InputTag("diTauMEtResolutionUp") )
process.selectedDiTauMEtResolutionUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauMEtResolutionUp"))

process.diTauMEtResolutionDown =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                            srcLeg1 = cms.InputTag("muPtEtaIDIso"),
                                            srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                            srcMET  = cms.InputTag("metRecoilCorrector",  "ResolutionD")
                                            )
process.selectedDiTauMEtResolutionDown = process.selectedDiTau.clone(src = cms.InputTag("diTauMEtResolutionDown") )
process.selectedDiTauMEtResolutionDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauMEtResolutionDown"))


process.diTauMuUp = process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                        srcLeg1 = cms.InputTag("rescaledMuons","U"),
                                        srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                        srcMET  = cms.InputTag("rescaledMETmuon","NNUNN")
                                        )
process.selectedDiTauMuUp = process.selectedDiTau.clone(src = cms.InputTag("diTauMuUp") )
process.selectedDiTauMuUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauMuUp"))

process.diTauMuDown = process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                          srcLeg1 = cms.InputTag("rescaledMuons","D"),
                                          srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                          srcMET  = cms.InputTag("rescaledMETmuon","NNDNN")
                                          )
process.selectedDiTauMuDown = process.selectedDiTau.clone(src = cms.InputTag("diTauMuDown") )
process.selectedDiTauMuDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauMuDown"))



process.diTauTauUp = process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                         srcLeg1 = cms.InputTag("muPtEtaIDIso"),
                                         srcLeg2 = cms.InputTag("rescaledTaus", "U"),
                                         srcMET  = cms.InputTag("rescaledMETtau","NNNUN")
                                         )
process.selectedDiTauTauUp = process.selectedDiTau.clone(src = cms.InputTag("diTauTauUp") )
process.selectedDiTauTauUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauTauUp"))

process.diTauTauDown = process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                           srcLeg1 = cms.InputTag("muPtEtaIDIso"),
                                           srcLeg2 = cms.InputTag("rescaledTaus", "D"),
                                           srcMET  = cms.InputTag("rescaledMETtau","NNNDN")
                                           )
process.selectedDiTauTauDown = process.selectedDiTau.clone(src = cms.InputTag("diTauTauDown") )
process.selectedDiTauTauDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauTauDown"))

process.allDiTau = cms.Sequence(
    (process.diTau*process.selectedDiTau*process.selectedDiTauCounter)+
    (process.diTauJetUp*process.selectedDiTauJetUp*process.selectedDiTauJetUpCounter +
     process.diTauJetDown*process.selectedDiTauJetDown*process.selectedDiTauJetDownCounter) +
    (process.diTauMEtResolutionUp*process.selectedDiTauMEtResolutionUp*process.selectedDiTauMEtResolutionUpCounter +
     process.diTauMEtResolutionDown*process.selectedDiTauMEtResolutionDown*process.selectedDiTauMEtResolutionDownCounter) +
    (process.diTauMEtResponseUp*process.selectedDiTauMEtResponseUp*process.selectedDiTauMEtResponseUpCounter +
     process.diTauMEtResponseDown*process.selectedDiTauMEtResponseDown*process.selectedDiTauMEtResponseDownCounter) +
    (process.diTauMuUp*process.selectedDiTauMuUp*process.selectedDiTauMuUpCounter +
     process.diTauMuDown*process.selectedDiTauMuDown*process.selectedDiTauMuDownCounter) +
    (process.diTauTauUp*process.selectedDiTauTauUp*process.selectedDiTauTauUpCounter +
     process.diTauTauDown*process.selectedDiTauTauDown*process.selectedDiTauTauDownCounter)
    )
#######################################################################

process.tauPtEtaIDAgMuAgElecIso  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("tauPtEtaIDAgMuAgElec"),
    cut = cms.string("tauID('byLooseCombinedIsolationDeltaBetaCorr')>0.5 && pt>20 && abs(eta)<2.3"),
    filter = cms.bool(False)
    )

process.tauPtEtaIDAgMuAgElecIsoPtRel = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("tauPtEtaIDAgMuAgElec"),
    cut = cms.string("tauID('byLooseCombinedIsolationDeltaBetaCorr')>0.5 && pt>19 && abs(eta)<2.3"),
    filter = cms.bool(False)
    )

process.tauPtEtaIDAgMuAgElecIsoCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )
process.tauPtEtaIDAgMuAgElecIsoTauUp  =  process.tauPtEtaIDAgMuAgElecIso.clone(
    src = cms.InputTag("rescaledTaus", "U")
    )
process.tauPtEtaIDAgMuAgElecIsoTauUpCounter = process.tauPtEtaIDAgMuAgElecIsoCounter.clone(
    src = cms.InputTag("tauPtEtaIDAgMuAgElecIsoTauUp"),
    )
process.tauPtEtaIDAgMuAgElecIsoTauDown  =  process.tauPtEtaIDAgMuAgElecIso.clone(
    src = cms.InputTag("rescaledTaus", "D")
    )
process.tauPtEtaIDAgMuAgElecIsoTauDownCounter = process.tauPtEtaIDAgMuAgElecIsoCounter.clone(
    src = cms.InputTag("tauPtEtaIDAgMuAgElecIsoTauDown"),
    )


process.muPtEtaIDIso  = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("muPtEtaID"),
    cut = cms.string("userFloat('PFRelIsoDB04')<0.50 && pt>15 && abs(eta)<2.4"),
    filter = cms.bool(False)
    )
process.muPtEtaIDIsoPtRel  = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("muPtEtaID"),
    cut = cms.string("userFloat('PFRelIsoDB04')<0.50 && pt>14 && abs(eta)<2.4"),
    filter = cms.bool(False)
    )

process.muPtEtaIDIsoCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("muPtEtaIDIso"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )
process.muPtEtaIDIsoMuUp = process.muPtEtaIDIso.clone(
    src = cms.InputTag("rescaledMuons","U")
    )
process.muPtEtaIDIsoMuUpCounter = process.muPtEtaIDIsoCounter.clone(
    src = cms.InputTag("muPtEtaIDIsoMuUp"),
    )
process.muPtEtaIDIsoMuDown = process.muPtEtaIDIso.clone(
    src = cms.InputTag("rescaledMuons","D")
    )
process.muPtEtaIDIsoMuDownCounter = process.muPtEtaIDIsoCounter.clone(
    src = cms.InputTag("muPtEtaIDIsoMuDown"),
    )
process.muPtEtaRelID = process.muPtEtaIDIso.clone(
    src = cms.InputTag("muPtEtaRelID"),
    cut = cms.string("pt>15")
    )
process.muPtEtaRelIDMuUp   = process.muPtEtaRelID.clone(
    src = cms.InputTag("rescaledMuonsRel","U")
    )
process.muPtEtaRelIDMuDown = process.muPtEtaRelID.clone(
    src = cms.InputTag("rescaledMuonsRel","D")
    )

process.filterSequence = cms.Sequence(
    (process.tauPtEtaIDAgMuAgElecIso       * process.tauPtEtaIDAgMuAgElecIsoCounter) +
    (process.tauPtEtaIDAgMuAgElecIsoTauUp  * process.tauPtEtaIDAgMuAgElecIsoTauUpCounter) +
    (process.tauPtEtaIDAgMuAgElecIsoTauDown* process.tauPtEtaIDAgMuAgElecIsoTauDownCounter) +
    (process.muPtEtaIDIso      * process.muPtEtaIDIsoCounter) +
    (process.muPtEtaIDIsoMuUp  * process.muPtEtaIDIsoMuUpCounter) +
    (process.muPtEtaIDIsoMuDown* process.muPtEtaIDIsoMuDownCounter) +
    (process.muPtEtaRelID+process.muPtEtaRelIDMuUp+process.muPtEtaRelIDMuDown)
    )

#######################################################################


process.muTauStreamAnalyzer = cms.EDAnalyzer(
    "MuTauStreamAnalyzer",
    diTaus         = cms.InputTag("selectedDiTau"),
    jets           = cms.InputTag("selectedPatJets"),
    newJets        = cms.InputTag(""),
    met            = cms.InputTag("metRecoilCorrector",  "N"),
    rawMet         = cms.InputTag("patMETsPF"),
    muons          = cms.InputTag("muPtEtaIDIso"),
    muonsRel       = cms.InputTag("muPtEtaRelID"),
    vertices       = cms.InputTag("selectedPrimaryVertices"),
    triggerResults = cms.InputTag("patTriggerEvent"),
    genParticles   = cms.InputTag("genParticles"),
    genTaus        = cms.InputTag("tauGenJetsSelectorAllHadrons"),
    isMC           = cms.bool(runOnMC),
    deltaRLegJet   = cms.untracked.double(0.5),
    minCorrPt      = cms.untracked.double(15.),
    minJetID       = cms.untracked.double(0.5), # 1=loose,2=medium,3=tight
    verbose        = cms.untracked.bool( False ),
    )
process.muTauStreamAnalyzerJetUp   = process.muTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauJetUp"),
    met    =  cms.InputTag("rescaledMETjet",  "UNNNU"),
    )
process.muTauStreamAnalyzerJetDown = process.muTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauJetDown"),
    met    =  cms.InputTag("rescaledMETjet",  "DNNND"),
    )
process.muTauStreamAnalyzerMEtResolutionUp   = process.muTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauMEtResolutionUp"),
    met    =  cms.InputTag("metRecoilCorrector",  "ResolutionU"),
    )
process.muTauStreamAnalyzerMEtResolutionDown = process.muTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauMEtResolutionDown"),
    met    =  cms.InputTag("metRecoilCorrector",  "ResolutionD"),
    )
process.muTauStreamAnalyzerMEtResponseUp   = process.muTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauMEtResponseUp"),
    met    =  cms.InputTag("metRecoilCorrector",  "ResponseU"),
    )
process.muTauStreamAnalyzerMEtResponseDown = process.muTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauMEtResponseDown"),
    met    =  cms.InputTag("metRecoilCorrector",  "ResponseD"),
    )
process.muTauStreamAnalyzerMuUp    = process.muTauStreamAnalyzer.clone(
    diTaus   =  cms.InputTag("selectedDiTauMuUp"),
    met      =  cms.InputTag("rescaledMETmuon","NNUNN"),
    muons    =  cms.InputTag("muPtEtaIDIsoMuUp"),
    muonsRel =  cms.InputTag("muPtEtaRelIDMuUp"),
    )
process.muTauStreamAnalyzerMuDown  = process.muTauStreamAnalyzer.clone(
    diTaus   =  cms.InputTag("selectedDiTauMuDown"),
    met      =  cms.InputTag("rescaledMETmuon","NNDNN"),
    muons    =  cms.InputTag("muPtEtaIDIsoMuDown"),
    muonsRel =  cms.InputTag("muPtEtaRelIDMuDown"),
    )
process.muTauStreamAnalyzerTauUp   = process.muTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauTauUp"),
    met    =  cms.InputTag("rescaledMETtau","NNNUN"),
    )
process.muTauStreamAnalyzerTauDown = process.muTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauTauDown"),
    met    =  cms.InputTag("rescaledMETtau","NNNDN")
    )

process.allAnalyzers = cms.Sequence(
    process.muTauStreamAnalyzer+
    process.muTauStreamAnalyzerJetUp+
    process.muTauStreamAnalyzerJetDown+
    process.muTauStreamAnalyzerMEtResponseUp+
    process.muTauStreamAnalyzerMEtResponseDown+
    process.muTauStreamAnalyzerMEtResolutionUp+
    process.muTauStreamAnalyzerMEtResolutionDown+
    process.muTauStreamAnalyzerMuUp+
    process.muTauStreamAnalyzerMuDown+
    process.muTauStreamAnalyzerTauUp+
    process.muTauStreamAnalyzerTauDown
    )
#######################################################################

process.analysis = cms.Sequence(
    process.allEventsFilter*
    process.filterSequence*
    process.rescaledObjects*
    process.allDiTau*
    process.allAnalyzers
    )

#######################################################################

if runOnMC:
    
    process.pNominal = cms.Path(
        process.allEventsFilter*
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        (process.muPtEtaIDIso *process.muPtEtaIDIsoCounter) *
        process.muPtEtaRelID *
        process.metRecoilCorrector*
        process.diTau*process.selectedDiTau*process.selectedDiTauCounter*
        process.muTauStreamAnalyzer
        )
    '''
    process.pJetUp = cms.Path(
        process.allEventsFilter*
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        (process.muPtEtaIDIso *process.muPtEtaIDIsoCounter) *
        process.muPtEtaRelID *
        process.metRecoilCorrector*
        process.rescaledMETjet *
        process.diTauJetUp*process.selectedDiTauJetUp*process.selectedDiTauJetUpCounter*
        process.muTauStreamAnalyzerJetUp
        )
    process.pJetDown = cms.Path(
        process.allEventsFilter*
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        (process.muPtEtaIDIso *process.muPtEtaIDIsoCounter) *
        process.muPtEtaRelID *
        process.metRecoilCorrector*
        process.rescaledMETjet *
        process.diTauJetDown*process.selectedDiTauJetDown*process.selectedDiTauJetDownCounter*
        process.muTauStreamAnalyzerJetDown
        )
    '''
    '''
    process.pMEtResolutionUp = cms.Path(
        process.allEventsFilter*
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        (process.muPtEtaIDIso *process.muPtEtaIDIsoCounter) *
        process.muPtEtaRelID *
        process.metRecoilCorrector*
        process.diTauMEtResolutionUp*process.selectedDiTauMEtResolutionUp*process.selectedDiTauMEtResolutionUpCounter*
        process.muTauStreamAnalyzerMEtResolutionUp
        )
    process.pMEtResolutionDown = cms.Path(
        process.allEventsFilter*
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        (process.muPtEtaIDIso *process.muPtEtaIDIsoCounter) *
        process.muPtEtaRelID *
        process.metRecoilCorrector*
        process.diTauMEtResolutionDown*process.selectedDiTauMEtResolutionDown*process.selectedDiTauMEtResolutionDownCounter*
        process.muTauStreamAnalyzerMEtResolutionDown
        )

    process.pMEtResponseUp = cms.Path(
        process.allEventsFilter*
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        (process.muPtEtaIDIso *process.muPtEtaIDIsoCounter) *
        process.muPtEtaRelID *
        process.metRecoilCorrector*
        process.diTauMEtResponseUp*process.selectedDiTauMEtResponseUp*process.selectedDiTauMEtResponseUpCounter*
        process.muTauStreamAnalyzerMEtResponseUp
        )
    process.pMEtResponseDown = cms.Path(
        process.allEventsFilter*
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        (process.muPtEtaIDIso *process.muPtEtaIDIsoCounter) *
        process.muPtEtaRelID *
        process.metRecoilCorrector*
        process.diTauMEtResponseDown*process.selectedDiTauMEtResponseDown*process.selectedDiTauMEtResponseDownCounter*
        process.muTauStreamAnalyzerMEtResponseDown
        )
    process.pMuUp = cms.Path(
    process.allEventsFilter*
    process.muPtEtaIDIsoPtRel *
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.metRecoilCorrector*
    (process.rescaledMETmuon+process.rescaledMuons+process.rescaledMuonsRel)*
    (process.muPtEtaIDIsoMuUp*process.muPtEtaIDIsoMuUpCounter) *
    process.muPtEtaRelIDMuUp *
    process.diTauMuUp*process.selectedDiTauMuUp*process.selectedDiTauMuUpCounter*
    process.muTauStreamAnalyzerMuUp
    )
    process.pMuDown = cms.Path(
    process.allEventsFilter*
    process.muPtEtaIDIsoPtRel *
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.metRecoilCorrector*
    (process.muPtEtaIDIsoMuUp*process.muPtEtaIDIsoMuUpCounter) *
    (process.muPtEtaIDIsoMuDown*process.muPtEtaIDIsoMuDownCounter) *
    process.muPtEtaRelIDMuDown *
    process.diTauMuDown*process.selectedDiTauMuDown*process.selectedDiTauMuDownCounter*
    process.muTauStreamAnalyzerMuDown
    )
    '''

    '''
    process.pTauUp = cms.Path(
        process.allEventsFilter*
        (process.muPtEtaIDIso*process.muPtEtaIDIsoCounter) *
        process.tauPtEtaIDAgMuAgElecIsoPtRel*
        process.muPtEtaRelID *
        process.metRecoilCorrector*
        (process.rescaledMETtau+process.rescaledTaus)*
        (process.tauPtEtaIDAgMuAgElecIsoTauUp*process.tauPtEtaIDAgMuAgElecIsoTauUpCounter)*
        process.diTauTauUp*process.selectedDiTauTauUp*process.selectedDiTauTauUpCounter*
        process.muTauStreamAnalyzerTauUp
    )
    process.pTauDown = cms.Path(
        process.allEventsFilter*
        (process.muPtEtaIDIso*process.muPtEtaIDIsoCounter) *
        process.tauPtEtaIDAgMuAgElecIsoPtRel*
        process.muPtEtaRelID *
        process.metRecoilCorrector*
        (process.rescaledMETtau+process.rescaledTaus)*
        (process.tauPtEtaIDAgMuAgElecIsoTauDown*process.tauPtEtaIDAgMuAgElecIsoTauDownCounter)*
        process.diTauTauDown*process.selectedDiTauTauDown*process.selectedDiTauTauDownCounter*
        process.muTauStreamAnalyzerTauDown
        )
    '''
else:
    process.pNominal = cms.Path(
        process.allEventsFilter*
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        (process.muPtEtaIDIso *process.muPtEtaIDIsoCounter) *
        process.muPtEtaRelID *
        process.metRecoilCorrector*
        process.diTau*process.selectedDiTau*process.selectedDiTauCounter*
        process.muTauStreamAnalyzer
        )
    '''
    process.pTauUp = cms.Path(
        process.allEventsFilter*
        (process.muPtEtaIDIso*process.muPtEtaIDIsoCounter) *
        process.tauPtEtaIDAgMuAgElecIsoPtRel*
        process.muPtEtaRelID *
        process.metRecoilCorrector*
        (process.rescaledMETtau+process.rescaledTaus)*
        (process.tauPtEtaIDAgMuAgElecIsoTauUp*process.tauPtEtaIDAgMuAgElecIsoTauUpCounter)*
        process.diTauTauUp*process.selectedDiTauTauUp*process.selectedDiTauTauUpCounter*
        process.muTauStreamAnalyzerTauUp
    )
    process.pTauDown = cms.Path(
        process.allEventsFilter*
        (process.muPtEtaIDIso*process.muPtEtaIDIsoCounter) *
        process.tauPtEtaIDAgMuAgElecIsoPtRel*
        process.muPtEtaRelID *
        process.metRecoilCorrector*
        (process.rescaledMETtau+process.rescaledTaus)*
        (process.tauPtEtaIDAgMuAgElecIsoTauDown*process.tauPtEtaIDAgMuAgElecIsoTauDownCounter)*
        process.diTauTauDown*process.selectedDiTauTauDown*process.selectedDiTauTauDownCounter*
        process.muTauStreamAnalyzerTauDown
        )
        
    process.pMuUp = cms.Path(
    process.allEventsFilter*
    process.muPtEtaIDIsoPtRel *
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.metRecoilCorrector*
    (process.rescaledMETmuon+process.rescaledMuons+process.rescaledMuonsRel)*
    (process.muPtEtaIDIsoMuUp*process.muPtEtaIDIsoMuUpCounter) *
    process.muPtEtaRelIDMuUp *
    process.diTauMuUp*process.selectedDiTauMuUp*process.selectedDiTauMuUpCounter*
    process.muTauStreamAnalyzerMuUp
    )
    process.pMuDown = cms.Path(
    process.allEventsFilter*
    process.muPtEtaIDIsoPtRel *
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.metRecoilCorrector*
    (process.muPtEtaIDIsoMuUp*process.muPtEtaIDIsoMuUpCounter) *
    (process.muPtEtaIDIsoMuDown*process.muPtEtaIDIsoMuDownCounter) *
    process.muPtEtaRelIDMuDown *
    process.diTauMuDown*process.selectedDiTauMuDown*process.selectedDiTauMuDownCounter*
    process.muTauStreamAnalyzerMuDown
    )
    '''

process.out = cms.OutputModule(
    "PoolOutputModule",
    outputCommands = cms.untracked.vstring( 'drop *',
                                            'keep *_metRecoilCorrector_*_*'),
    fileName = cms.untracked.string('patTuplesSkimmed_MuTauStream.root'),
    )

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("treeMuTauStream.root")
    )


## To work on Artur's skim
#from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
#massSearchReplaceAnyInputTag(process.pNominal,
#                             "muPtEtaID",
#                             "selectedPatMuonsTriggerMatch",
#                             verbose=False)
#massSearchReplaceAnyInputTag(process.pNominal,
#                             "muPtEtaRelID",
#                             "selectedPatMuonsTriggerMatch",
#                             verbose=False)
#massSearchReplaceAnyInputTag(process.pNominal,
#                             "selectedPrimaryVertices",
#                             "offlinePrimaryVertices",
#                             verbose=False)
#massSearchReplaceAnyInputTag(process.pNominal,
#                            "tauPtEtaIDAgMuAgElec",
#                             "selectedPatTausTriggerMatch",
#                             verbose=False)
#massSearchReplaceAnyInputTag(process.pNominal,
#                             "genParticles",
#                             "prunedGenParticles",
#                             verbose=False)
#massSearchReplaceAnyInputTag(process.pNominal,
#                             "tauGenJetsSelectorAllHadrons",
#                             "genTauDecaysToHadrons",
#                             verbose=False)


process.outpath = cms.EndPath()

processDumpFile = open('runMuTauStreamAnalyzerFullAnalysis_Recoil.dump', 'w')
print >> processDumpFile, process.dumpPython()

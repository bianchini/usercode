import FWCore.ParameterSet.Config as cms

process = cms.Process("ELECTAUANA")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

#from Configuration.PyReleaseValidation.autoCond import autoCond
#process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

runOnMC     = False
doSVFitReco = True

if runOnMC:
    print "Running on MC"
else:
    print "Running on Data"


if runOnMC:
    process.GlobalTag.globaltag = cms.string('START52_V7::All')

else:
    process.GlobalTag.globaltag = cms.string('GR_R_52_V7::All')
    
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    #'file:./patTuples_ElecTauStream.root'

    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_101_1_Q8V.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_102_1_fyV.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_103_1_xEA.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_107_2_Y1L.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_109_1_5pG.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_11_2_Wtu.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_12_1_jIz.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_13_1_mCk.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_16_1_WeO.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_17_2_BI3.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_18_1_k1t.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_19_2_sty.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_1_1_UQS.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_20_1_oEP.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_22_2_UDX.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_23_1_PWw.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_24_1_ya3.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_26_1_Bx1.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_2_1_GCq.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_30_1_0RW.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_32_1_f6Q.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_34_1_QJs.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_36_1_nKS.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_38_2_TNX.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_39_1_psz.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_3_1_wfg.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_41_1_7Im.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_42_1_EoO.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_45_2_uTH.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_46_1_r84.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_47_2_PYf.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_48_1_WJA.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_49_2_8A2.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_4_1_lTn.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_51_2_CZ6.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_52_1_Pjo.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_53_2_j8d.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_54_1_ogs.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_56_1_ArY.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_58_1_LjP.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_59_1_R7Z.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_5_1_3Yc.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_60_1_FWU.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_62_2_5OG.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_63_1_esB.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_64_1_e2W.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_65_1_IDs.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_66_1_zed.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_69_2_XWA.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_6_1_PQz.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_70_1_oBJ.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_71_1_OuU.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_72_1_CtV.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_75_1_MzX.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_76_1_fXJ.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_7_1_3GB.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_80_1_gB1.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_82_1_J23.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_83_1_RHl.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_84_1_3Fz.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_86_1_6eb.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_87_2_4Er.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_8_1_Hmt.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_90_1_kr2.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_93_1_b2e.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_94_1_fnZ.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_95_2_Cdi.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_97_1_ZM9.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_98_1_juL.root',
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-madgraph-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_99_1_DRT.root',


    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_10_2_9pt.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_11_2_GDL.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_12_2_4In.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_13_2_Cbc.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_1_2_6iw.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_2_2_8Oi.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_3_2_qlY.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_4_2_3B3.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_5_2_Ntc.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_6_2_fT2.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_7_2_OBt.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_8_2_AfI.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/ElecTauStream-30Mar2012_DYJets-ElecTau-Tarball-v2_skim/c8840777bcf9a6298dca78985d99645d/patTuples_ElecTauStream_9_2_pdJ.root',
    )
    )

#process.source.eventsToProcess = cms.untracked.VEventRange(
#    '163659:555:404811294','163659:556:405921039','163659:557:406814110'
#    )

process.allEventsFilter = cms.EDFilter(
    "AllEventsFilter"
    )

###################################################################################
process.rescaledMET = cms.EDProducer(
    "MEtRescalerProducer",
    metTag          = cms.InputTag("metRecoilCorrector",  "N"),
    jetTag          = cms.InputTag("selectedPatJets"),
    electronTag     = cms.InputTag("elecPtEtaIDIso"),
    muonTag         = cms.InputTag(""),
    tauTag          = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
    unClusterShift  = cms.double(0.10),
    tauShift        = cms.vdouble(0.03,0.03),
    muonShift       = cms.vdouble(0.01,0.01),
    electronShift   = cms.vdouble(0.01,0.025),
    jetThreshold    = cms.double(10),
    numOfSigmas     = cms.double(1.0),
    verbose         = cms.bool(False)
    )

process.rescaledMETjet = process.rescaledMET.clone(
    unClusterShift = cms.double(0.10),
    tauShift       = cms.vdouble(0.0),
    muonShift      = cms.vdouble(0.0),
    electronShift  = cms.vdouble(0.0),
    )
process.rescaledMETtau = process.rescaledMET.clone(
    unClusterShift = cms.double(0.0),
    tauShift       = cms.vdouble(0.03,0.03),
    muonShift      = cms.vdouble(0.0),
    electronShift  = cms.vdouble(0.0),
    )
process.rescaledMETelectron = process.rescaledMET.clone(
    unClusterShift = cms.double(0.0),
    tauShift       = cms.vdouble(0.0),
    muonShift      = cms.vdouble(0.0),
    electronShift  = cms.vdouble(0.01,0.025),
    )

process.rescaledTaus = cms.EDProducer(
    "TauRescalerProducer",
    inputCollection = cms.InputTag("tauPtEtaIDAgMuAgElecIsoPtRel"),
    shift           = cms.vdouble(0.03,0.03),
    numOfSigmas     = cms.double(1.0)
    )
process.rescaledElectrons = cms.EDProducer(
    "ElectronRescalerProducer",
    inputCollection = cms.InputTag("elecPtEtaIDIsoPtRel"),
    shift           = cms.vdouble(0.01,0.025),
    numOfSigmas     = cms.double(1.0),
    )
process.rescaledElectronsRel = cms.EDProducer(
    "ElectronRescalerProducer",
    inputCollection = cms.InputTag("elecPtEtaRelID"),
    shift           = cms.vdouble(0.01,0.025),
    numOfSigmas     = cms.double(1.0),
    )

process.rescaledObjects = cms.Sequence(
    process.rescaledMETjet+
    process.rescaledMETtau+
    process.rescaledMETelectron+
    process.rescaledTaus+
    process.rescaledElectrons+
    process.rescaledElectronsRel
    )

###################################################################################

process.metRecoilCorrector = cms.EDProducer(
    "MEtRecoilCorrectorProducer",
    genParticleTag      = cms.InputTag("genParticles"),
    jetTag              = cms.InputTag("selectedPatJets"),
    metTag              = cms.InputTag("patMETsPF"),
    electronTag         = cms.InputTag("elecPtEtaIDIso"),
    muonTag             = cms.InputTag(""),
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
process.diTau = process.allElecTauPairs.clone()
process.diTau.srcLeg1  = cms.InputTag("elecPtEtaIDIso")
process.diTau.srcLeg2  = cms.InputTag("tauPtEtaIDAgMuAgElecIso")
process.diTau.srcMET   = cms.InputTag("metRecoilCorrector",  "N")
process.diTau.dRmin12  = cms.double(0.5)
process.diTau.doSVreco = cms.bool(doSVFitReco)

if not runOnMC:
    process.diTau.srcGenParticles = ""
        
process.selectedDiTau = cms.EDFilter(
    "ElecTauPairSelector",
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
                                          srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                          srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                          srcMET  = cms.InputTag("rescaledMETjet",  "UNNNU")
                                          )
process.selectedDiTauJetUp = process.selectedDiTau.clone(src = cms.InputTag("diTauJetUp") )
process.selectedDiTauJetUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauJetUp"))

process.diTauJetDown =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                            srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                            srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                            srcMET  = cms.InputTag("rescaledMETjet",  "DNNND")
                                            )
process.selectedDiTauJetDown = process.selectedDiTau.clone(src = cms.InputTag("diTauJetDown") )
process.selectedDiTauJetDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauJetDown"))

process.diTauMEtResponseUp =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                          srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                          srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                          srcMET  = cms.InputTag("metRecoilCorrector",  "ResponseU")
                                          )
process.selectedDiTauMEtResponseUp = process.selectedDiTau.clone(src = cms.InputTag("diTauMEtResponseUp") )
process.selectedDiTauMEtResponseUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauMEtResponseUp"))

process.diTauMEtResponseDown =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                            srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                            srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                            srcMET  = cms.InputTag("metRecoilCorrector",  "ResponseD")
                                            )
process.selectedDiTauMEtResponseDown = process.selectedDiTau.clone(src = cms.InputTag("diTauMEtResponseDown") )
process.selectedDiTauMEtResponseDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauMEtResponseDown"))


process.diTauMEtResolutionUp =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                          srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                          srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                          srcMET  = cms.InputTag("metRecoilCorrector",  "ResolutionU")
                                          )
process.selectedDiTauMEtResolutionUp = process.selectedDiTau.clone(src = cms.InputTag("diTauMEtResolutionUp") )
process.selectedDiTauMEtResolutionUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauMEtResolutionUp"))

process.diTauMEtResolutionDown =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                            srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                            srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                            srcMET  = cms.InputTag("metRecoilCorrector",  "ResolutionD")
                                            )
process.selectedDiTauMEtResolutionDown = process.selectedDiTau.clone(src = cms.InputTag("diTauMEtResolutionDown") )
process.selectedDiTauMEtResolutionDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauMEtResolutionDown"))

process.diTauElecUp = process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                          srcLeg1 = cms.InputTag("rescaledElectrons","U"),
                                          srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                          srcMET  = cms.InputTag("rescaledMETelectron","NUNNN")
                                          )
process.selectedDiTauElecUp = process.selectedDiTau.clone(src = cms.InputTag("diTauElecUp") )
process.selectedDiTauElecUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauElecUp"))

process.diTauElecDown = process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                            srcLeg1 = cms.InputTag("rescaledElectrons","D"),
                                            srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                            srcMET  = cms.InputTag("rescaledMETelectron","NDNNN")
                                            )
process.selectedDiTauElecDown = process.selectedDiTau.clone(src = cms.InputTag("diTauElecDown") )
process.selectedDiTauElecDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauElecDown"))


process.diTauTauUp = process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                         srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                         srcLeg2 = cms.InputTag("rescaledTaus", "U"),
                                         srcMET  = cms.InputTag("rescaledMETtau","NNNUN")
                                         )
process.selectedDiTauTauUp = process.selectedDiTau.clone(src = cms.InputTag("diTauTauUp") )
process.selectedDiTauTauUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauTauUp"))

process.diTauTauDown = process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                           srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
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
    (process.diTauElecUp*process.selectedDiTauElecUp*process.selectedDiTauElecUpCounter +
     process.diTauElecDown*process.selectedDiTauElecDown*process.selectedDiTauElecDownCounter) +
    (process.diTauTauUp*process.selectedDiTauTauUp*process.selectedDiTauTauUpCounter +
     process.diTauTauDown*process.selectedDiTauTauDown*process.selectedDiTauTauDownCounter)
    )
#######################################################################

MVA = "((pt<=20 && abs(superClusterPosition.Eta)>=0.0 && abs(superClusterPosition.Eta)<1.0 && userFloat('mva')>0.133) ||" + \
      " (pt<=20 && abs(superClusterPosition.Eta)>=1.0 && abs(superClusterPosition.Eta)<1.5 && userFloat('mva')>0.465) ||" + \
      " (pt<=20 && abs(superClusterPosition.Eta)>=1.5 && abs(superClusterPosition.Eta)<2.5 && userFloat('mva')>0.518) ||" + \
      " (pt>20  && abs(superClusterPosition.Eta)>=0.0 && abs(superClusterPosition.Eta)<1.0 && userFloat('mva')>0.942) ||" + \
      " (pt>20  && abs(superClusterPosition.Eta)>=1.0 && abs(superClusterPosition.Eta)<1.5 && userFloat('mva')>0.947) ||" + \
      " (pt>20  && abs(superClusterPosition.Eta)>=1.5 && abs(superClusterPosition.Eta)<2.5 && userFloat('mva')>0.878) )"

simpleCutsWP95 = "(userFloat('nHits')<=1"+ \
                 " && (" + \
                 " (isEB && userFloat('sihih')<0.010 && userFloat('dPhi')<0.80 && "+ \
                 "          userFloat('dEta') <0.007 && userFloat('HoE') <0.15)"   + \
                 " || "  + \
                 " (isEE && userFloat('sihih')<0.030 && userFloat('dPhi')<0.70 && "+ \
                 "          userFloat('dEta') <0.010 && userFloat('HoE') <0.07)"   + \
                 "     )"+ \
                 ")"
simpleCutsWP80 = "(userFloat('nHits')==0 && userInt('antiConv')>0.5 "+ \
                 " && ("   + \
                 " (pt>=20 && ("    + \
                 "               (isEB && userFloat('sihih')<0.010 && userFloat('dPhi')<0.06 && "  + \
                 "                        userFloat('dEta')< 0.004 && userFloat('HoE') <0.04)"     + \
                 "               ||"+ \
                 "               (isEE && userFloat('sihih')<0.030 && userFloat('dPhi')<0.030 && " + \
                 "                        userFloat('dEta') <0.007 && userFloat('HoE') <0.025) )) "+ \
                 "     || "+ \
                 " (pt<20 && (fbrem>0.15 || (abs(superClusterPosition.Eta)<1. && eSuperClusterOverP>0.95) ) && ( "+ \
                 "               (isEB && userFloat('sihih')<0.010 && userFloat('dPhi')<0.030 && " + \
                 "                        userFloat('dEta') <0.004 && userFloat('HoE') <0.025) "   + \
                 "               ||"+ \
                 "               (isEE && userFloat('sihih')<0.030 && userFloat('dPhi')<0.020 &&"  + \
                 "                        userFloat('dEta') <0.005 && userFloat('HoE') <0.025) ))" + \
                 "    )"   + \
                 ")"


process.tauPtEtaIDAgMuAgElecIso  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("tauPtEtaIDAgMuAgElec"),
    cut = cms.string("tauID('byLooseCombinedIsolationDeltaBetaCorr')>0.5 && pt>20 && abs(eta)<2.3"
                     "&& tauID('againstElectronMVA')>0.5"
                     ),
    filter = cms.bool(False)
    )
process.tauPtEtaIDAgMuAgElecIsoPtRel  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("tauPtEtaIDAgMuAgElec"),
    cut = cms.string("tauID('byLooseCombinedIsolationDeltaBetaCorr')>0.5 && pt>19 && abs(eta)<2.3"
                     "&& tauID('againstElectronMVA')>0.5"
                     ),
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


process.elecPtEtaIDIso  = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("elecPtEtaID"),
    cut = cms.string("userFloat('PFRelIsoDB04')<0.50 && pt>20 && abs(eta)<2.1 && "+simpleCutsWP95),
    filter = cms.bool(False)
    )
process.elecPtEtaIDIsoPtRel  = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("elecPtEtaID"),
    cut = cms.string("userFloat('PFRelIsoDB04')<0.50 && pt>19 && abs(eta)<2.1 && "+simpleCutsWP95),
    filter = cms.bool(False)
    )

process.elecPtEtaIDIsoCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("elecPtEtaIDIso"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )
process.elecPtEtaIDIsoElecUp = process.elecPtEtaIDIso.clone(
    src = cms.InputTag("rescaledElectrons","U")
    )
process.elecPtEtaIDIsoElecUpCounter = process.elecPtEtaIDIsoCounter.clone(
    src = cms.InputTag("elecPtEtaIDIsoElecUp"),
    )
process.elecPtEtaIDIsoElecDown = process.elecPtEtaIDIso.clone(
    src = cms.InputTag("rescaledElectrons","D")
    )
process.elecPtEtaIDIsoElecDownCounter = process.elecPtEtaIDIsoCounter.clone(
    src = cms.InputTag("elecPtEtaIDIsoElecDown"),
    )
process.elecPtEtaRelID = process.elecPtEtaIDIso.clone(
    src = cms.InputTag("elecPtEtaRelID"),
    cut = cms.string("pt>15")
    )
process.elecPtEtaRelIDElecUp   = process.elecPtEtaRelID.clone(
    src = cms.InputTag("rescaledElectronsRel","U")
    )
process.elecPtEtaRelIDElecDown = process.elecPtEtaRelID.clone(
    src = cms.InputTag("rescaledElectronsRel","D")
    )

process.filterSequence = cms.Sequence(
    (process.tauPtEtaIDAgMuAgElecIso       * process.tauPtEtaIDAgMuAgElecIsoCounter) +
    (process.tauPtEtaIDAgMuAgElecIsoTauUp  * process.tauPtEtaIDAgMuAgElecIsoTauUpCounter) +
    (process.tauPtEtaIDAgMuAgElecIsoTauDown* process.tauPtEtaIDAgMuAgElecIsoTauDownCounter) +
    (process.elecPtEtaIDIso                * process.elecPtEtaIDIsoCounter) +
    (process.elecPtEtaIDIsoElecUp          * process.elecPtEtaIDIsoElecUpCounter) +
    (process.elecPtEtaIDIsoElecDown        * process.elecPtEtaIDIsoElecDownCounter) +
    (process.elecPtEtaRelID+process.elecPtEtaRelIDElecUp+process.elecPtEtaRelIDElecDown)
    )

#######################################################################


process.elecTauStreamAnalyzer = cms.EDAnalyzer(
    "ElecTauStreamAnalyzer",
    diTaus             = cms.InputTag("selectedDiTau"),
    jets               = cms.InputTag("selectedPatJets"),
    newJets            = cms.InputTag(""),
    met                = cms.InputTag("metRecoilCorrector",  "N"),
    rawMet             = cms.InputTag("patMETsPF"),
    electrons          = cms.InputTag("elecPtEtaID"),
    electronsRel       = cms.InputTag("elecPtEtaRelID"),
    vertices           = cms.InputTag("selectedPrimaryVertices"),
    triggerResults     = cms.InputTag("patTriggerEvent"),
    genParticles       = cms.InputTag("genParticles"),
    genTaus            = cms.InputTag("tauGenJetsSelectorAllHadrons"),
    isMC               = cms.bool(runOnMC),
    deltaRLegJet       = cms.untracked.double(0.5),
    minCorrPt          = cms.untracked.double(15.),
    minJetID           = cms.untracked.double(0.5), # 1=loose,2=medium,3=tight
    #inputFileNameX0BL  = cms.FileInPath("Bianchi/Utilities/data/antiE_v4/TMVAClassification_v2_X_0BL_BDT.weights.xml"),
    #inputFileName11BL  = cms.FileInPath("Bianchi/Utilities/data/antiE_v4/TMVAClassification_v2_1_1BL_BDT.weights.xml"),
    #inputFileName01BL  = cms.FileInPath("Bianchi/Utilities/data/antiE_v4/TMVAClassification_v2_0_1BL_BDT.weights.xml"),
    #inputFileNameX0EC  = cms.FileInPath("Bianchi/Utilities/data/antiE_v4/TMVAClassification_v2_X_0EC_BDT.weights.xml"),
    #inputFileName11EC  = cms.FileInPath("Bianchi/Utilities/data/antiE_v4/TMVAClassification_v2_1_1EC_BDT.weights.xml"),
    #inputFileName01EC  = cms.FileInPath("Bianchi/Utilities/data/antiE_v4/TMVAClassification_v2_0_1EC_BDT.weights.xml"),
    verbose            = cms.untracked.bool( False ),
    )
process.elecTauStreamAnalyzerJetUp     = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauJetUp"),
    met    =  cms.InputTag("rescaledMETjet",  "UNNNU"),
    )
process.elecTauStreamAnalyzerJetDown   = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauJetDown"),
    met    =  cms.InputTag("rescaledMETjet",  "DNNND"),
    )
process.elecTauStreamAnalyzerMEtResponseUp   = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauMEtResponseUp"),
    met    =  cms.InputTag("metRecoilCorrector",  "ResponseU"),
    )
process.elecTauStreamAnalyzerMEtResponseDown = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauMEtResponseDown"),
    met    =  cms.InputTag("metRecoilCorrector",  "ResponseD"),
    )
process.elecTauStreamAnalyzerMEtResolutionUp  = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauMEtResolutionUp"),
    met    =  cms.InputTag("metRecoilCorrector",  "ResolutionU"),
    )
process.elecTauStreamAnalyzerMEtResolutionDown = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauMEtResolutionDown"),
    met    =  cms.InputTag("metRecoilCorrector",  "ResolutionD"),
    )
process.elecTauStreamAnalyzerElecUp    = process.elecTauStreamAnalyzer.clone(
    diTaus       =  cms.InputTag("selectedDiTauElecUp"),
    met          =  cms.InputTag("rescaledMETelectron","NUNNN"),
    electrons    =  cms.InputTag("elecPtEtaIDIsoElecUp"),
    electronsRel =  cms.InputTag("elecPtEtaRelIDElecUp"),
    )
process.elecTauStreamAnalyzerElecDown  = process.elecTauStreamAnalyzer.clone(
    diTaus       =  cms.InputTag("selectedDiTauElecDown"),
    met          =  cms.InputTag("rescaledMETelectron","NDNNN"),
    electrons    =  cms.InputTag("elecPtEtaIDIsoElecDown"),
    electronsRel =  cms.InputTag("elecPtEtaRelIDElecDown"),
    )
process.elecTauStreamAnalyzerTauUp     = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauTauUp"),
    met    =  cms.InputTag("rescaledMETtau","NNNUN")
    )
process.elecTauStreamAnalyzerTauDown   = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauTauDown"),
    met    =  cms.InputTag("rescaledMETtau","NNNDN")
    )

process.allAnalyzers = cms.Sequence(
    process.elecTauStreamAnalyzer+
    process.elecTauStreamAnalyzerJetUp+
    process.elecTauStreamAnalyzerJetDown+
    process.elecTauStreamAnalyzerMEtResponseUp+
    process.elecTauStreamAnalyzerMEtResponseDown+
    process.elecTauStreamAnalyzerMEtResolutionUp+
    process.elecTauStreamAnalyzerMEtResolutionDown+
    process.elecTauStreamAnalyzerElecUp+
    process.elecTauStreamAnalyzerElecDown+
    process.elecTauStreamAnalyzerTauUp+
    process.elecTauStreamAnalyzerTauDown
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
        (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
        process.elecPtEtaRelID *
        process.metRecoilCorrector*
        process.diTau*process.selectedDiTau*process.selectedDiTauCounter*
        process.elecTauStreamAnalyzer
        )
    '''
    process.pJetUp = cms.Path(
        process.allEventsFilter*
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
        process.elecPtEtaRelID *
        process.metRecoilCorrector*
        process.rescaledMETjet *
        process.diTauJetUp*process.selectedDiTauJetUp*process.selectedDiTauJetUpCounter*
        process.elecTauStreamAnalyzerJetUp
        )
    process.pJetDown = cms.Path(
        process.allEventsFilter*
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
        process.elecPtEtaRelID *
        process.metRecoilCorrector*
        process.rescaledMETjet *
        process.diTauJetDown*process.selectedDiTauJetDown*process.selectedDiTauJetDownCounter*
        process.elecTauStreamAnalyzerJetDown
        )
    process.pMEtResolutionUp = cms.Path(
        process.allEventsFilter*
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
        process.elecPtEtaRelID *
        process.metRecoilCorrector*
        process.diTauMEtResolutionUp*process.selectedDiTauMEtResolutionUp*process.selectedDiTauMEtResolutionUpCounter*
        process.elecTauStreamAnalyzerMEtResolutionUp
        )
    process.pMEtResolutionDown = cms.Path(
        process.allEventsFilter*
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
        process.elecPtEtaRelID *
        process.metRecoilCorrector*
        process.diTauMEtResolutionDown*process.selectedDiTauMEtResolutionDown*process.selectedDiTauMEtResolutionDownCounter*
        process.elecTauStreamAnalyzerMEtResolutionDown
        )
    
    process.pMEtResponseUp = cms.Path(
        process.allEventsFilter*
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
        process.elecPtEtaRelID *
        process.metRecoilCorrector*
        process.diTauMEtResponseUp*process.selectedDiTauMEtResponseUp*process.selectedDiTauMEtResponseUpCounter*
        process.elecTauStreamAnalyzerMEtResponseUp
        )
    process.pMEtResponseDown = cms.Path(
        process.allEventsFilter*
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
        process.elecPtEtaRelID *
        process.metRecoilCorrector*
        process.diTauMEtResponseDown*process.selectedDiTauMEtResponseDown*process.selectedDiTauMEtResponseDownCounter*
        process.elecTauStreamAnalyzerMEtResponseDown
        )

    
    process.pElecUp = cms.Path(
        process.allEventsFilter*
        process.elecPtEtaIDIsoPtRel *
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        process.metRecoilCorrector*
        (process.rescaledMETelectron+process.rescaledElectrons+process.rescaledElectronsRel)*
        (process.elecPtEtaIDIsoElecUp*process.elecPtEtaIDIsoElecUpCounter) *
        process.elecPtEtaRelIDElecUp *
        process.diTauElecUp*process.selectedDiTauElecUp*process.selectedDiTauElecUpCounter*
        process.elecTauStreamAnalyzerElecUp
        )
    process.pElecDown = cms.Path(
        process.allEventsFilter*
        process.elecPtEtaIDIsoPtRel *
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        process.metRecoilCorrector*
        (process.rescaledMETelectron+process.rescaledElectrons+process.rescaledElectronsRel)*
        (process.elecPtEtaIDIsoElecDown*process.elecPtEtaIDIsoElecDownCounter) *
        process.elecPtEtaRelIDElecDown *
        process.diTauElecDown*process.selectedDiTauElecDown*process.selectedDiTauElecDownCounter*
        process.elecTauStreamAnalyzerElecDown
        )
    process.pTauUp = cms.Path(
        process.allEventsFilter*
        (process.elecPtEtaIDIso*process.elecPtEtaIDIsoCounter) *
        process.tauPtEtaIDAgMuAgElecIsoPtRel*
        process.elecPtEtaRelID *
        process.metRecoilCorrector*
        (process.rescaledMETtau+process.rescaledTaus)*
        (process.tauPtEtaIDAgMuAgElecIsoTauUp*process.tauPtEtaIDAgMuAgElecIsoTauUpCounter)*
        process.diTauTauUp*process.selectedDiTauTauUp*process.selectedDiTauTauUpCounter*
        process.elecTauStreamAnalyzerTauUp
        )
    process.pTauDown = cms.Path(
        process.allEventsFilter*
        (process.elecPtEtaIDIso*process.elecPtEtaIDIsoCounter) *
        process.tauPtEtaIDAgMuAgElecIsoPtRel*
        process.elecPtEtaRelID *
        process.metRecoilCorrector*
        (process.rescaledMETtau+process.rescaledTaus)*
        (process.tauPtEtaIDAgMuAgElecIsoTauDown*process.tauPtEtaIDAgMuAgElecIsoTauDownCounter)*
        process.diTauTauDown*process.selectedDiTauTauDown*process.selectedDiTauTauDownCounter*
        process.elecTauStreamAnalyzerTauDown
        )
    '''
    
else:
    
    process.pNominal = cms.Path(
        process.allEventsFilter*
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
        process.elecPtEtaRelID *
        process.metRecoilCorrector*
        process.diTau*process.selectedDiTau*process.selectedDiTauCounter*
        process.elecTauStreamAnalyzer
        )
    '''
    process.pTauUp = cms.Path(
        process.allEventsFilter*
        (process.elecPtEtaIDIso*process.elecPtEtaIDIsoCounter) *
        process.tauPtEtaIDAgMuAgElecIsoPtRel*
        process.elecPtEtaRelID *
        process.metRecoilCorrector*
        (process.rescaledMETtau+process.rescaledTaus)*
        (process.tauPtEtaIDAgMuAgElecIsoTauUp*process.tauPtEtaIDAgMuAgElecIsoTauUpCounter)*
        process.diTauTauUp*process.selectedDiTauTauUp*process.selectedDiTauTauUpCounter*
        process.elecTauStreamAnalyzerTauUp
        )
    process.pTauDown = cms.Path(
        process.allEventsFilter*
        (process.elecPtEtaIDIso*process.elecPtEtaIDIsoCounter) *
        process.tauPtEtaIDAgMuAgElecIsoPtRel*
        process.elecPtEtaRelID *
        process.metRecoilCorrector*
        (process.rescaledMETtau+process.rescaledTaus)*
        (process.tauPtEtaIDAgMuAgElecIsoTauDown*process.tauPtEtaIDAgMuAgElecIsoTauDownCounter)*
        process.diTauTauDown*process.selectedDiTauTauDown*process.selectedDiTauTauDownCounter*
        process.elecTauStreamAnalyzerTauDown
        )
       
    process.pElecUp = cms.Path(
        process.allEventsFilter*
        process.elecPtEtaIDIsoPtRel *
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        process.metRecoilCorrector*
        (process.rescaledMETelectron+process.rescaledElectrons+process.rescaledElectronsRel)*
        (process.elecPtEtaIDIsoElecUp*process.elecPtEtaIDIsoElecUpCounter) *
        process.elecPtEtaRelIDElecUp *
        process.diTauElecUp*process.selectedDiTauElecUp*process.selectedDiTauElecUpCounter*
        process.elecTauStreamAnalyzerElecUp
        )
    process.pElecDown = cms.Path(
        process.allEventsFilter*
        process.elecPtEtaIDIsoPtRel *
        (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
        process.metRecoilCorrector*
        (process.rescaledMETelectron+process.rescaledElectrons+process.rescaledElectronsRel)*
        (process.elecPtEtaIDIsoElecDown*process.elecPtEtaIDIsoElecDownCounter) *
        process.elecPtEtaRelIDElecDown *
        process.diTauElecDown*process.selectedDiTauElecDown*process.selectedDiTauElecDownCounter*
        process.elecTauStreamAnalyzerElecDown
        )
    '''

process.out = cms.OutputModule(
    "PoolOutputModule",
    outputCommands = cms.untracked.vstring( 'drop *'),
    fileName = cms.untracked.string('patTuplesSkimmed_ElecTauStream.root'),
    )

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("treeElecTauStream.root")
    )

process.outpath = cms.EndPath()

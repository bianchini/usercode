import FWCore.ParameterSet.Config as cms

process = cms.Process("MUTAUANA")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

runOnMC     = True
doSVFitReco = False

if runOnMC:
    print "Running on MC"
else:
    print "Running on Data"
        

if runOnMC:
    process.GlobalTag.globaltag = cms.string('START42_V17::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V23::All')
    
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    'file:./patTuples_MuTauStream.root'
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/MuTauStream-16Nov2011/700226eee9a93cb10580e91b3d6e5c18/patTuples_MuTauStream_10_1_w41.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/GluGluToHToTauTau_M-120_7TeV-powheg-pythia6/MuTauStream-16Nov2011/cea16f7ce381f249953e69976c5a9109/patTuples_MuTauStream_10_1_yRu.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/SUSYGluGluToHToTauTau_M-500_7TeV-pythia6-tauola/MuTauStream-16Nov2011/d67d72c0e0080ba7d4d6c9945ad37260/patTuples_MuTauStream_10_1_8Pf.root'
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi//TauPlusX/MuTauStream-13Oct2011-05AugReReco/c37c208594d74fa447903aef959eea7d/patTuples_MuTauStream_14_1_MHD.root'

    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_10_1_rrD.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_11_1_lXv.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_12_1_ttU.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_13_1_Lfe.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_14_1_lVu.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_15_1_EJp.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_16_1_4Si.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_17_1_MaR.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_18_1_YDx.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_19_1_d86.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_1_1_VIM.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_20_1_zSB.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_21_1_osQ.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_22_1_XVu.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_23_1_IiZ.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_24_1_2Cj.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_25_1_bSv.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_26_1_UFe.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_27_1_dXQ.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_28_1_msA.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_29_1_uAN.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_2_1_V4V.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_30_1_6Oc.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_31_1_i7T.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_32_1_HB5.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_33_1_siv.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_34_1_YZA.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_35_1_LZD.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_36_1_ypv.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_37_1_MED.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_38_1_k2p.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_39_1_W4J.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_3_1_PtW.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_40_1_FRE.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_4_1_zIn.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_5_1_RzE.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_6_1_WxC.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_7_1_HCE.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_8_1_363.root',
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/MuTauStream-02May2012-Test_VBFH120-MuTau-powheg-PUS6_skim/0eb1f9a73e1c7ada0f6173ee6aa7d17e/patTuples_MuTauStream_9_1_129.root',    
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
    metTag              = cms.InputTag("patMETsPFlow"),
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
    cut = cms.string("pt>20 && abs(eta)<2.3"+
                     " && (tauID('byLooseCombinedIsolationDeltaBetaCorr')>0.5 || tauID('byLooseIsolationMVA')>0.5)"
                     ),
    filter = cms.bool(False)
    )
process.tauPtEtaIDAgMuAgElecIsoPtRel  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("tauPtEtaIDAgMuAgElec"),
    cut = cms.string("pt>19 && abs(eta)<2.3"+
                     " && (tauID('byLooseCombinedIsolationDeltaBetaCorr')>0.5 || tauID('byLooseIsolationMVA')>0.5)"
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


process.muPtEtaIDIso  = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("muPtEtaID"),
    cut = cms.string("userFloat('PFRelIsoDB04')<0.50 && pt>15 && abs(eta)<2.1"),
    filter = cms.bool(False)
    )
process.muPtEtaIDIsoPtRel  = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("muPtEtaID"),
    cut = cms.string("userFloat('PFRelIsoDB04')<0.50 && pt>14 && abs(eta)<2.1"),
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
    rawMet         = cms.InputTag("patMETsPFlow"),
    mvaMet         = cms.InputTag("patPFMetByMVA"),
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
    doMuIsoMVA     = cms.bool( True ),
    inputFileName0 = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt5To10_V1_BDTG.weights.xml"),
    inputFileName1 = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt5To10_V1_BDTG.weights.xml"),
    inputFileName2 = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt10ToInf_V1_BDTG.weights.xml"),
    inputFileName3 = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt10ToInf_V1_BDTG.weights.xml"),                 
    inputFileName4 = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Tracker_V1_BDTG.weights.xml"),
    inputFileName5 = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Global_V1_BDTG.weights.xml"),
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

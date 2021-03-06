import FWCore.ParameterSet.Config as cms

process = cms.Process("ELECTAUANA")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

runOnMC     = True
doSVFitReco = True
usePFMEtMVA = False

if runOnMC:
    print "Running on MC"
else:
    print "Running on Data"


if runOnMC:
    process.GlobalTag.globaltag = cms.string('START42_V14B::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V19::All')
    
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/ElecTauStream-04May2012-Reload_DYJets-ElecTau-50-madgraph-PUS6_skim/4badcc5695438d7f3df80162f5ad7ed7/patTuples_ElecTauStream_9_1_Ug5.root'
    'file:./root/patTuples_ElecTauStream_VBFH125.root'
    #'file:./patTuples_ElecTauStream.root'
    #'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/TauPlusX/ElecTauStream-04May2012-Reload-05AugReReco/396c4fb61647929194f9a223b98504bc/patTuples_ElecTauStream_9_1_kgg.root'
    )
    )

#process.source.eventsToProcess = cms.untracked.VEventRange(
#    '1:108429','1:157957','1:157990'
#    )

process.allEventsFilter = cms.EDFilter(
    "AllEventsFilter"
    )

###################################################################################

process.load("RecoMET.METProducers.mvaPFMET_cff")
if runOnMC:
    process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3")
else:
    process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3Residual")
    

process.pfMEtMVA.srcLeptons = cms.VInputTag( cms.InputTag('elecPtEtaIDIso'), cms.InputTag('tauPtEtaIDAgMuAgElecIso') )

process.load("PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi")
process.patPFMetByMVA = process.patMETs.clone(
    metSource = cms.InputTag('pfMEtMVA'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue')
    )

process.pfMEtMVA.srcVertices = cms.InputTag("selectedPrimaryVertices")
process.patPFMetByMVA.addGenMET = cms.bool(False)


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

if usePFMEtMVA:
    process.rescaledMET.metTag = cms.InputTag("patPFMetByMVA")
process.rescaledMETRaw = process.rescaledMET.clone( metTag = cms.InputTag("metRecoilCorrector",  "N") )

process.rescaledMETjet = process.rescaledMET.clone(
    unClusterShift = cms.double(0.10),
    tauShift       = cms.vdouble(0.0),
    muonShift      = cms.vdouble(0.0),
    electronShift  = cms.vdouble(0.0),
    )
process.rescaledMETRawjet = process.rescaledMETjet.clone( metTag = cms.InputTag("metRecoilCorrector",  "N") )

process.rescaledMETtau = process.rescaledMET.clone(
    unClusterShift = cms.double(0.0),
    tauShift       = cms.vdouble(0.03,0.03),
    muonShift      = cms.vdouble(0.0),
    electronShift  = cms.vdouble(0.0),
    )
process.rescaledMETRawtau = process.rescaledMETtau.clone( metTag = cms.InputTag("metRecoilCorrector",  "N") )

process.rescaledMETelectron = process.rescaledMET.clone(
    unClusterShift = cms.double(0.0),
    tauShift       = cms.vdouble(0.0),
    muonShift      = cms.vdouble(0.0),
    electronShift  = cms.vdouble(0.01,0.025),
    )
process.rescaledMETRawelectron = process.rescaledMETelectron.clone( metTag = cms.InputTag("metRecoilCorrector",  "N") )


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
    process.rescaledMETRawjet+
    process.rescaledMETRawtau+
    process.rescaledMETRawelectron+
    process.rescaledTaus+
    process.rescaledElectrons+
    process.rescaledElectronsRel
    )

###################################################################################

process.metRecoilCorrector = cms.EDProducer(
    "MEtRecoilCorrectorProducer",
    genParticleTag      = cms.InputTag("genParticles"),
    jetTag              = cms.InputTag("selectedPatJets"),
    metTag              = cms.InputTag("patMETsPFlow"),
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
if usePFMEtMVA:
    process.diTau.srcMET = cms.InputTag("patPFMetByMVA")
    
if not runOnMC:
    process.diTau.srcGenParticles = ""

process.pfMEtMVACov = cms.EDProducer(
    "PFMEtSignCovMatrixUnembedder",
    src = cms.InputTag("patPFMetByMVA")
    )

if usePFMEtMVA:
    process.diTau.nSVfit.psKine_MEt_logM_fit.config.event.srcMEt = cms.InputTag("patPFMetByMVA")
    process.diTau.nSVfit.psKine_MEt_logM_fit.config.event.likelihoodFunctions[0].srcMEtCovMatrix = cms.InputTag("pfMEtMVACov")
    process.diTau.nSVfit.psKine_MEt_logM_int.config.event.srcMEt = cms.InputTag("patPFMetByMVA")
    process.diTau.nSVfit.psKine_MEt_logM_int.config.event.likelihoodFunctions[0].srcMEtCovMatrix = cms.InputTag("pfMEtMVACov")

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


process.diTauRaw = process.allElecTauPairs.clone()
process.diTauRaw.srcLeg1  = cms.InputTag("elecPtEtaIDIso")
process.diTauRaw.srcLeg2  = cms.InputTag("tauPtEtaIDAgMuAgElecIso")
process.diTauRaw.srcMET   = cms.InputTag("metRecoilCorrector",  "N")
process.diTauRaw.dRmin12  = cms.double(0.5)
process.diTauRaw.doSVreco = cms.bool(doSVFitReco)

if not runOnMC:
    process.diTauRaw.srcGenParticles = ""
        
process.diTauRaw.nSVfit.psKine_MEt_logM_fit.config.event.srcMEt = cms.InputTag("metRecoilCorrector",  "N")
process.diTauRaw.nSVfit.psKine_MEt_logM_int.config.event.srcMEt = cms.InputTag("metRecoilCorrector",  "N")

process.selectedDiTauRaw =  process.selectedDiTau.clone(   src = cms.InputTag("diTauRaw") )
process.selectedDiTauRawCounter =  process.selectedDiTauCounter.clone(   src = cms.InputTag("selectedDiTauRaw") )

#######################################################################
#######################################################################

process.diTauJetUp =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                          srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                          srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                          srcMET  = cms.InputTag("rescaledMETjet",  "UNNNU")
                                          )
process.selectedDiTauJetUp = process.selectedDiTau.clone(src = cms.InputTag("diTauJetUp") )
process.selectedDiTauJetUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauJetUp"))

process.diTauRawJetUp                = process.diTauJetUp.clone(  srcMET  = cms.InputTag("rescaledMETRawjet",  "UNNNU") )
process.selectedDiTauRawJetUp        = process.selectedDiTauJetUp.clone( src = cms.InputTag("diTauRawJetUp") )
process.selectedDiTauRawJetUpCounter = process.selectedDiTauRawCounter.clone(src =  cms.InputTag("selectedDiTauRawJetUp"))

#######################################################################

process.diTauJetDown =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                            srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                            srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                            srcMET  = cms.InputTag("rescaledMETjet",  "DNNND")
                                            )
process.selectedDiTauJetDown = process.selectedDiTau.clone(src = cms.InputTag("diTauJetDown") )
process.selectedDiTauJetDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauJetDown"))

process.diTauRawJetDown                = process.diTauJetDown.clone(  srcMET  = cms.InputTag("rescaledMETRawjet",  "DNNND") )
process.selectedDiTauRawJetDown        = process.selectedDiTauJetDown.clone( src = cms.InputTag("diTauRawJetDown") )
process.selectedDiTauRawJetDownCounter = process.selectedDiTauRawCounter.clone(src =  cms.InputTag("selectedDiTauRawJetDown"))

#######################################################################

process.diTauMEtResponseUp =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                          srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                          srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                          srcMET  = cms.InputTag("metRecoilCorrector",  "ResponseU")
                                          )
process.selectedDiTauMEtResponseUp = process.selectedDiTau.clone(src = cms.InputTag("diTauMEtResponseUp") )
process.selectedDiTauMEtResponseUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauMEtResponseUp"))

process.diTauRawMEtResponseUp                = process.diTauMEtResponseUp.clone(  srcMET  = cms.InputTag("metRecoilCorrector",  "ResponseU") )
process.selectedDiTauRawMEtResponseUp        = process.selectedDiTauMEtResponseUp.clone( src = cms.InputTag("diTauRawMEtResponseUp") )
process.selectedDiTauRawMEtResponseUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauRawMEtResponseUp"))

#######################################################################

process.diTauMEtResponseDown =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                            srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                            srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                            srcMET  = cms.InputTag("metRecoilCorrector",  "ResponseD")
                                            )
process.selectedDiTauMEtResponseDown = process.selectedDiTau.clone(src = cms.InputTag("diTauMEtResponseDown") )
process.selectedDiTauMEtResponseDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauMEtResponseDown"))

process.diTauRawMEtResponseDown                = process.diTauMEtResponseDown.clone(  srcMET  = cms.InputTag("metRecoilCorrector",  "ResponseU") )
process.selectedDiTauRawMEtResponseDown        = process.selectedDiTauMEtResponseDown.clone( src = cms.InputTag("diTauRawMEtResponseDown") )
process.selectedDiTauRawMEtResponseDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauRawMEtResponseDown"))

#######################################################################


process.diTauMEtResolutionUp =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                          srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                          srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                          srcMET  = cms.InputTag("metRecoilCorrector",  "ResolutionU")
                                          )
process.selectedDiTauMEtResolutionUp = process.selectedDiTau.clone(src = cms.InputTag("diTauMEtResolutionUp") )
process.selectedDiTauMEtResolutionUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauMEtResolutionUp"))

process.diTauRawMEtResolutionUp                = process.diTauMEtResolutionUp.clone(  srcMET  = cms.InputTag("metRecoilCorrector",  "ResolutionU") )
process.selectedDiTauRawMEtResolutionUp        = process.selectedDiTauMEtResolutionUp.clone( src = cms.InputTag("diTauRawMEtResolutionUp") )
process.selectedDiTauRawMEtResolutionUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauRawMEtResolutionUp"))

#######################################################################

process.diTauMEtResolutionDown =  process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                            srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                            srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                            srcMET  = cms.InputTag("metRecoilCorrector",  "ResolutionD")
                                            )
process.selectedDiTauMEtResolutionDown = process.selectedDiTau.clone(src = cms.InputTag("diTauMEtResolutionDown") )
process.selectedDiTauMEtResolutionDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauMEtResolutionDown"))

process.diTauRawMEtResolutionDown                = process.diTauMEtResolutionDown.clone(  srcMET  = cms.InputTag("metRecoilCorrector",  "ResolutionU") )
process.selectedDiTauRawMEtResolutionDown        = process.selectedDiTauMEtResolutionDown.clone( src = cms.InputTag("diTauRawMEtResolutionDown") )
process.selectedDiTauRawMEtResolutionDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauRawMEtResolutionDown"))


#######################################################################


process.diTauElecUp = process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                          srcLeg1 = cms.InputTag("rescaledElectrons","U"),
                                          srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                          srcMET  = cms.InputTag("rescaledMETelectron","NUNNN")
                                          )
process.selectedDiTauElecUp = process.selectedDiTau.clone(src = cms.InputTag("diTauElecUp") )
process.selectedDiTauElecUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauElecUp"))

process.diTauRawElecUp                = process.diTauElecUp.clone(  srcMET  = cms.InputTag("rescaledMETRawelectron",  "NUNNN") )
process.selectedDiTauRawElecUp        = process.selectedDiTauElecUp.clone( src = cms.InputTag("diTauRawElecUp") )
process.selectedDiTauRawElecUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauRawElecUp"))

#######################################################################

process.diTauElecDown = process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                            srcLeg1 = cms.InputTag("rescaledElectrons","D"),
                                            srcLeg2 = cms.InputTag("tauPtEtaIDAgMuAgElecIso"),
                                            srcMET  = cms.InputTag("rescaledMETelectron","NDNNN")
                                            )
process.selectedDiTauElecDown = process.selectedDiTau.clone(src = cms.InputTag("diTauElecDown") )
process.selectedDiTauElecDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauElecDown"))

process.diTauRawElecDown                = process.diTauElecDown.clone(  srcMET  = cms.InputTag("rescaledMETRawelectron",  "NDNNN") )
process.selectedDiTauRawElecDown        = process.selectedDiTauElecDown.clone( src = cms.InputTag("diTauRawElecDown") )
process.selectedDiTauRawElecDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauRawElecDown"))

#######################################################################


process.diTauTauUp = process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                         srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                         srcLeg2 = cms.InputTag("rescaledTaus", "U"),
                                         srcMET  = cms.InputTag("rescaledMETtau","NNNUN")
                                         )
process.selectedDiTauTauUp = process.selectedDiTau.clone(src = cms.InputTag("diTauTauUp") )
process.selectedDiTauTauUpCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauTauUp"))

process.diTauRawTauUp                = process.diTauRaw.clone(doSVreco = cms.bool(doSVFitReco),
                                                              srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                                              srcLeg2 = cms.InputTag("rescaledTaus", "U"),
                                                              srcMET  = cms.InputTag("rescaledMETRawtau","NNNUN")
                                                              )
process.selectedDiTauRawTauUp        = process.selectedDiTauRaw.clone( src = cms.InputTag("diTauRawTauUp") )
process.selectedDiTauRawTauUpCounter = process.selectedDiTauRawCounter.clone(src =  cms.InputTag("selectedDiTauRawTauUp"))

#######################################################################


process.diTauTauDown = process.diTau.clone(doSVreco = cms.bool(doSVFitReco),
                                           srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                           srcLeg2 = cms.InputTag("rescaledTaus", "D"),
                                           srcMET  = cms.InputTag("rescaledMETtau","NNNDN")
                                           )
process.selectedDiTauTauDown = process.selectedDiTau.clone(src = cms.InputTag("diTauTauDown") )
process.selectedDiTauTauDownCounter = process.selectedDiTauCounter.clone(src =  cms.InputTag("selectedDiTauTauDown"))

process.diTauRawTauDown                = process.diTauRaw.clone(doSVreco = cms.bool(doSVFitReco),
                                                                srcLeg1 = cms.InputTag("elecPtEtaIDIso"),
                                                                srcLeg2 = cms.InputTag("rescaledTaus", "U"),
                                                                srcMET  = cms.InputTag("rescaledMETRawtau","NNNDN")
                                                                )
process.selectedDiTauRawTauDown        = process.selectedDiTauRaw.clone( src = cms.InputTag("diTauRawTauDown") )
process.selectedDiTauRawTauDownCounter = process.selectedDiTauRawCounter.clone(src =  cms.InputTag("selectedDiTauRawTauDown"))

#######################################################################


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
     process.diTauTauDown*process.selectedDiTauTauDown*process.selectedDiTauTauDownCounter) +
    (process.diTauRaw*process.selectedDiTauRaw*process.selectedDiTauRawCounter)+
    (process.diTauRawJetUp*process.selectedDiTauRawJetUp*process.selectedDiTauRawJetUpCounter +
     process.diTauRawJetDown*process.selectedDiTauRawJetDown*process.selectedDiTauRawJetDownCounter) +
    (process.diTauRawMEtResolutionUp*process.selectedDiTauRawMEtResolutionUp*process.selectedDiTauRawMEtResolutionUpCounter +
     process.diTauRawMEtResolutionDown*process.selectedDiTauRawMEtResolutionDown*process.selectedDiTauRawMEtResolutionDownCounter) +
    (process.diTauRawMEtResponseUp*process.selectedDiTauRawMEtResponseUp*process.selectedDiTauRawMEtResponseUpCounter +
     process.diTauRawMEtResponseDown*process.selectedDiTauRawMEtResponseDown*process.selectedDiTauRawMEtResponseDownCounter) +
    (process.diTauRawElecUp*process.selectedDiTauRawElecUp*process.selectedDiTauRawElecUpCounter +
     process.diTauRawElecDown*process.selectedDiTauRawElecDown*process.selectedDiTauRawElecDownCounter) +
    (process.diTauRawTauUp*process.selectedDiTauRawTauUp*process.selectedDiTauRawTauUpCounter +
     process.diTauRawTauDown*process.selectedDiTauRawTauDown*process.selectedDiTauRawTauDownCounter)

    )

if usePFMEtMVA:
    process.allDiTau.replace(process.diTau,process.pfMEtMVACov*process.diTau)
#    process.allDiTau.replace(process.diTauJetUp,process.pfMEtMVACovJetUp*process.diTauJetUp)
#    process.allDiTau.replace(process.diTauJetDown,process.pfMEtMVACovJetDown*process.diTauJetDown)
#    process.allDiTau.replace(process.diTauElecUp,process.pfMEtMVACovElecUp*process.diTauElecUp)
#    process.allDiTau.replace(process.diTauElecDown,process.pfMEtMVACovElecDown*process.diTauElecDown)
#    process.allDiTau.replace(process.diTauTauUp,process.pfMEtMVACovTauUp*process.diTauTauUp)
#    process.allDiTau.replace(process.diTauTauDown,process.pfMEtMVACovTauDown*process.diTauTauDown)

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
    cut = cms.string("pt>20 && abs(eta)<2.3"+
                     " && tauID('byLooseIsolationMVA')>-0.5"+
                     " && tauID('againstElectronMVA')>-0.5"
                     ),
    filter = cms.bool(False)
    )
process.tauPtEtaIDAgMuAgElecIsoPtRel  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("tauPtEtaIDAgMuAgElec"),
    cut = cms.string("pt>19 && abs(eta)<2.3"+
                     " && tauID('byLooseIsolationMVA')>-0.5"+
                     " && tauID('againstElectronMVA')>-0.5"
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

#################################################
process.elecPtEtaID = cms.EDProducer(
    "ElectronsUserEmbeddedIso",
    electronTag = cms.InputTag("elecPtEtaID"),
    )
#################################################

process.elecPtEtaIDIso  = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("elecPtEtaID"),
    cut = cms.string("userFloat('PFRelIsoDB04v2')<0.5 && pt>20 && abs(eta)<2.1"+
                     "&& userInt('antiConv')>0.5 && userInt('nHits')<1"),
    filter = cms.bool(False)
    )
process.elecPtEtaIDIsoPtRel  = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("elecPtEtaID"),
    cut = cms.string("userFloat('PFRelIsoDB04v2')<0.5 && pt>19 && abs(eta)<2.1"+
                     "&& userInt('antiConv')>0.5 && userInt('nHits')<1"),
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
    rawMet             = cms.InputTag("patMETsPFlow"),
    mvaMet             = cms.InputTag("patPFMetByMVA"),
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
    verbose            = cms.untracked.bool( False ),
    doElecIsoMVA       = cms.bool( True ),
    inputFileName0     = cms.FileInPath("UserCode/sixie/EGamma/EGammaAnalysisTools/data/ElectronIso_BDTG_V0_BarrelPt5To10.weights.xml"),
    inputFileName1     = cms.FileInPath("UserCode/sixie/EGamma/EGammaAnalysisTools/data/ElectronIso_BDTG_V0_EndcapPt5To10.weights.xml"),
    inputFileName2     = cms.FileInPath("UserCode/sixie/EGamma/EGammaAnalysisTools/data/ElectronIso_BDTG_V0_BarrelPt10ToInf.weights.xml"),
    inputFileName3     = cms.FileInPath("UserCode/sixie/EGamma/EGammaAnalysisTools/data/ElectronIso_BDTG_V0_EndcapPt10ToInf.weights.xml"),
    )

if usePFMEtMVA:
    process.elecTauStreamAnalyzer.met = cms.InputTag("patPFMetByMVA")

process.elecTauStreamAnalyzerRaw   = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauRaw"),
    met    = cms.InputTag("metRecoilCorrector",  "N"),
    )


process.elecTauStreamAnalyzerJetUp     = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauJetUp"),
    met    =  cms.InputTag("rescaledMETjet",  "UNNNU"),
    )
process.elecTauStreamAnalyzerJetDown   = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauJetDown"),
    met    =  cms.InputTag("rescaledMETjet",  "DNNND"),
    )
process.elecTauStreamAnalyzerRawJetUp   = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauRawJetUp"),
    met    =  cms.InputTag("rescaledMETRawjet",  "UNNNU"),
    )
process.elecTauStreamAnalyzerRawJetDown = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauRawJetDown"),
    met    =  cms.InputTag("rescaledMETRawjet",  "DNNND"),
    )

process.elecTauStreamAnalyzerMEtResponseUp   = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauMEtResponseUp"),
    met    =  cms.InputTag("metRecoilCorrector",  "ResponseU"),
    )
process.elecTauStreamAnalyzerMEtResponseDown = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauMEtResponseDown"),
    met    =  cms.InputTag("metRecoilCorrector",  "ResponseD"),
    )
process.elecTauStreamAnalyzerRawMEtResponseUp   = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauRawMEtResponseUp"),
    met    =  cms.InputTag("metRecoilCorrector",  "ResponseU"),
    )
process.elecTauStreamAnalyzerRawMEtResponseDown = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauRawMEtResponseDown"),
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
process.elecTauStreamAnalyzerRawMEtResolutionUp   = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauRawMEtResolutionUp"),
    met    =  cms.InputTag("metRecoilCorrector",  "ResolutionU"),
    )
process.elecTauStreamAnalyzerRawMEtResolutionDown = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauRawMEtResolutionDown"),
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
process.elecTauStreamAnalyzerRawElecUp    = process.elecTauStreamAnalyzer.clone(
    diTaus   =  cms.InputTag("selectedDiRawTauElecUp"),
    met      =  cms.InputTag("rescaledMETRawelectron","NUNNN"),
    muons    =  cms.InputTag("elecPtEtaIDIsoElecUp"),
    muonsRel =  cms.InputTag("elecPtEtaRelIDElecUp"),
    )
process.elecTauStreamAnalyzerRawElecDown  = process.elecTauStreamAnalyzer.clone(
    diTaus   =  cms.InputTag("selectedDiTauRawElecDown"),
    met      =  cms.InputTag("rescaledMETRawelectron","NDNNN"),
    muons    =  cms.InputTag("elecPtEtaIDIsoElecDown"),
    muonsRel =  cms.InputTag("elecPtEtaRelIDElecDown"),
    )


process.elecTauStreamAnalyzerTauUp     = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauTauUp"),
    met    =  cms.InputTag("rescaledMETtau","NNNUN")
    )
process.elecTauStreamAnalyzerTauDown   = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauTauDown"),
    met    =  cms.InputTag("rescaledMETtau","NNNDN")
    )
process.elecTauStreamAnalyzerRawTauUp   = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauRawTauUp"),
    met    =  cms.InputTag("rescaledMETRawtau","NNNUN"),
    )
process.elecTauStreamAnalyzerRawTauDown = process.elecTauStreamAnalyzer.clone(
    diTaus =  cms.InputTag("selectedDiTauRawTauDown"),
    met    =  cms.InputTag("rescaledMETRawtau","NNNDN")
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
    process.elecTauStreamAnalyzerTauDown+
    process.elecTauStreamAnalyzerRaw+
    process.elecTauStreamAnalyzerRawJetUp+
    process.elecTauStreamAnalyzerRawJetDown+
    process.elecTauStreamAnalyzerRawMEtResponseUp+
    process.elecTauStreamAnalyzerRawMEtResponseDown+
    process.elecTauStreamAnalyzerRawMEtResolutionUp+
    process.elecTauStreamAnalyzerRawMEtResolutionDown+
    process.elecTauStreamAnalyzerRawElecUp+
    process.elecTauStreamAnalyzerRawElecDown+
    process.elecTauStreamAnalyzerRawTauUp+
    process.elecTauStreamAnalyzerRawTauDown
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


process.seqNominal = cms.Sequence(
    process.allEventsFilter*
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
    process.elecPtEtaRelID *
    (process.pfMEtMVAsequence*process.patPFMetByMVA)*
    process.metRecoilCorrector*
    process.pfMEtMVACov*
    process.diTau*process.selectedDiTau*process.selectedDiTauCounter*
    process.elecTauStreamAnalyzer
    )
process.seqJetUp = cms.Sequence(
    process.allEventsFilter*
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
    process.elecPtEtaRelID *
    (process.pfMEtMVAsequence*process.patPFMetByMVA)*
    process.metRecoilCorrector*
    process.rescaledMETjet *
    process.pfMEtMVACov*
    process.diTauJetUp*process.selectedDiTauJetUp*process.selectedDiTauJetUpCounter*
    process.elecTauStreamAnalyzerJetUp
    )
process.seqJetDown = cms.Sequence(
    process.allEventsFilter*
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
    process.elecPtEtaRelID *
    (process.pfMEtMVAsequence*process.patPFMetByMVA)*
    process.metRecoilCorrector*
    process.rescaledMETjet *
    process.pfMEtMVACov*
    process.diTauJetDown*process.selectedDiTauJetDown*process.selectedDiTauJetDownCounter*
    process.elecTauStreamAnalyzerJetDown
    )

process.seqMEtResolutionUp = cms.Sequence(
    process.allEventsFilter*
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
    process.elecPtEtaRelID *
    process.metRecoilCorrector*
    process.pfMEtMVACov*
    process.diTauMEtResolutionUp*process.selectedDiTauMEtResolutionUp*process.selectedDiTauMEtResolutionUpCounter*
    process.elecTauStreamAnalyzerMEtResolutionUp
    )
process.seqMEtResolutionDown = cms.Sequence(
    process.allEventsFilter*
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
    process.elecPtEtaRelID *
    process.metRecoilCorrector*
    process.pfMEtMVACov*
    process.diTauMEtResolutionDown*process.selectedDiTauMEtResolutionDown*process.selectedDiTauMEtResolutionDownCounter*
    process.elecTauStreamAnalyzerMEtResolutionDown
    )
    
process.seqMEtResponseUp = cms.Sequence(
    process.allEventsFilter*
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
    process.elecPtEtaRelID *
    process.metRecoilCorrector*
    process.pfMEtMVACov*
    process.diTauMEtResponseUp*process.selectedDiTauMEtResponseUp*process.selectedDiTauMEtResponseUpCounter*
    process.elecTauStreamAnalyzerMEtResponseUp
    )
process.seqMEtResponseDown = cms.Sequence(
    process.allEventsFilter*
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
    process.elecPtEtaRelID *
    process.metRecoilCorrector*
    process.pfMEtMVACov*
    process.diTauMEtResponseDown*process.selectedDiTauMEtResponseDown*process.selectedDiTauMEtResponseDownCounter*
    process.elecTauStreamAnalyzerMEtResponseDown
    )

process.seqElecUp = cms.Sequence(
    process.allEventsFilter*
    process.elecPtEtaID*
    process.elecPtEtaIDIsoPtRel *
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    (process.pfMEtMVAsequence*process.patPFMetByMVA)*
    process.metRecoilCorrector*
    (process.rescaledMETelectron+process.rescaledElectrons+process.rescaledElectronsRel)*
    (process.elecPtEtaIDIsoElecUp*process.elecPtEtaIDIsoElecUpCounter) *
    process.elecPtEtaRelIDElecUp *
    process.pfMEtMVACov*
    process.diTauElecUp*process.selectedDiTauElecUp*process.selectedDiTauElecUpCounter*
    process.elecTauStreamAnalyzerElecUp
    )
process.seqElecDown = cms.Sequence(
    process.allEventsFilter*
    process.elecPtEtaID*
    process.elecPtEtaIDIsoPtRel *
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    (process.pfMEtMVAsequence*process.patPFMetByMVA)*
    process.metRecoilCorrector*
    (process.rescaledMETelectron+process.rescaledElectrons+process.rescaledElectronsRel)*
    (process.elecPtEtaIDIsoElecDown*process.elecPtEtaIDIsoElecDownCounter) *
    process.elecPtEtaRelIDElecDown *
    process.pfMEtMVACov*
    process.diTauElecDown*process.selectedDiTauElecDown*process.selectedDiTauElecDownCounter*
    process.elecTauStreamAnalyzerElecDown
    )

process.seqTauUp = cms.Sequence(
    process.allEventsFilter*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso*process.elecPtEtaIDIsoCounter) *
    process.tauPtEtaIDAgMuAgElecIsoPtRel*
    process.elecPtEtaRelID *
    (process.pfMEtMVAsequence*process.patPFMetByMVA)*
    process.metRecoilCorrector*
    (process.rescaledMETtau+process.rescaledTaus)*
    (process.tauPtEtaIDAgMuAgElecIsoTauUp*process.tauPtEtaIDAgMuAgElecIsoTauUpCounter)*
    process.pfMEtMVACov*
    process.diTauTauUp*process.selectedDiTauTauUp*process.selectedDiTauTauUpCounter*
    process.elecTauStreamAnalyzerTauUp
    )
process.seqTauDown = cms.Sequence(
    process.allEventsFilter*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso*process.elecPtEtaIDIsoCounter) *
    process.tauPtEtaIDAgMuAgElecIsoPtRel*
    process.elecPtEtaRelID *
    (process.pfMEtMVAsequence*process.patPFMetByMVA)*
    process.metRecoilCorrector*
    (process.rescaledMETtau+process.rescaledTaus)*
    (process.tauPtEtaIDAgMuAgElecIsoTauDown*process.tauPtEtaIDAgMuAgElecIsoTauDownCounter)*
    process.pfMEtMVACov*
    process.diTauTauDown*process.selectedDiTauTauDown*process.selectedDiTauTauDownCounter*
    process.elecTauStreamAnalyzerTauDown
    )


if not usePFMEtMVA:
    process.seqNominal.remove(process.pfMEtMVACov)
    process.seqJetUp.remove(process.pfMEtMVACov)
    process.seqJetDown.remove(process.pfMEtMVACov)
    process.seqMEtResolutionUp.remove(process.pfMEtMVACov)
    process.seqMEtResolutionDown.remove(process.pfMEtMVACov)
    process.seqMEtResponseUp.remove(process.pfMEtMVACov)
    process.seqMEtResponseDown.remove(process.pfMEtMVACov)
    process.seqElecUp.remove(process.pfMEtMVACov)
    process.seqElecDown.remove(process.pfMEtMVACov)
    process.seqTauUp.remove(process.pfMEtMVACov)
    process.seqTauDown.remove(process.pfMEtMVACov)


process.seqRawNominal = cms.Sequence(
    process.allEventsFilter*
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
    process.elecPtEtaRelID *
    process.metRecoilCorrector*
    #process.pfMEtMVACov*
    process.diTauRaw*process.selectedDiTauRaw*process.selectedDiTauRawCounter*
    process.elecTauStreamAnalyzerRaw
    )

process.seqRawJetUp = cms.Sequence(
    process.allEventsFilter*
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
    process.elecPtEtaRelID *
    process.metRecoilCorrector*
    process.rescaledMETRawjet *
    #process.pfMEtMVACov*
    process.diTauRawJetUp*process.selectedDiTauRawJetUp*process.selectedDiTauRawJetUpCounter*
    process.elecTauStreamAnalyzerRawJetUp
    )

process.seqRawJetDown = cms.Sequence(
    process.allEventsFilter*
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
    process.elecPtEtaRelID *
    process.metRecoilCorrector*
    process.rescaledMETRawjet *
    #process.pfMEtMVACov*
    process.diTauRawJetDown*process.selectedDiTauRawJetDown*process.selectedDiTauRawJetDownCounter*
    process.elecTauStreamAnalyzerRawJetDown
    )

process.seqRawMEtResolutionUp = cms.Sequence(
    process.allEventsFilter*
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
    process.elecPtEtaRelID *
    process.metRecoilCorrector*
    #process.pfMEtMVACov*
    process.diTauRawMEtResolutionUp*process.selectedDiTauRawMEtResolutionUp*process.selectedDiTauRawMEtResolutionUpCounter*
    process.elecTauStreamAnalyzerRawMEtResolutionUp
    )

process.seqRawMEtResolutionDown = cms.Sequence(
    process.allEventsFilter*
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
    process.elecPtEtaRelID *
    process.metRecoilCorrector*
    #process.pfMEtMVACov*
    process.diTauRawMEtResolutionDown*process.selectedDiTauRawMEtResolutionDown*process.selectedDiTauRawMEtResolutionDownCounter*
    process.elecTauStreamAnalyzerRawMEtResolutionDown
    )

process.seqRawMEtResponseUp = cms.Sequence(
    process.allEventsFilter*
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
    process.elecPtEtaRelID *
    process.metRecoilCorrector*
    #process.pfMEtMVACov*
    process.diTauRawMEtResponseUp*process.selectedDiTauRawMEtResponseUp*process.selectedDiTauRawMEtResponseUpCounter*
    process.elecTauStreamAnalyzerRawMEtResponseUp
    )

process.seqRawMEtResponseDown = cms.Sequence(
    process.allEventsFilter*
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso *process.elecPtEtaIDIsoCounter) *
    process.elecPtEtaRelID *
    process.metRecoilCorrector*
    #process.pfMEtMVACov*
    process.diTauRawMEtResponseDown*process.selectedDiTauRawMEtResponseDown*process.selectedDiTauRawMEtResponseDownCounter*
    process.elecTauStreamAnalyzerRawMEtResponseDown
    )

process.seqRawElecUp = cms.Sequence(
    process.allEventsFilter*
    process.elecPtEtaID*
    process.elecPtEtaIDIsoPtRel *
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.metRecoilCorrector*
    (process.rescaledMETRawelectron+process.rescaledElectrons+process.rescaledElectronsRel)*
    (process.elecPtEtaIDIsoElecUp*process.elecPtEtaIDIsoElecUpCounter) *
    process.elecPtEtaRelIDElecUp *
    #process.pfMEtMVACov*
    process.diTauRawElecUp*process.selectedDiTauRawElecUp*process.selectedDiTauRawElecUpCounter*
    process.elecTauStreamAnalyzerRawElecUp
    )

process.seqRawElecDown = cms.Sequence(
    process.allEventsFilter*
    process.elecPtEtaID*
    process.elecPtEtaIDIsoPtRel *
    (process.tauPtEtaIDAgMuAgElecIso*process.tauPtEtaIDAgMuAgElecIsoCounter)*
    process.metRecoilCorrector*
    (process.rescaledMETRawelectron+process.rescaledElectrons+process.rescaledElectronsRel)*
    (process.elecPtEtaIDIsoElecDown*process.elecPtEtaIDIsoElecDownCounter) *
    process.elecPtEtaRelIDElecDown *
    #process.pfMEtMVACov*
    process.diTauRawElecDown*process.selectedDiTauRawElecDown*process.selectedDiTauRawElecDownCounter*
    process.elecTauStreamAnalyzerRawElecDown
    )

process.seqRawTauUp = cms.Sequence(
    process.allEventsFilter*
    process.elecPtEtaID*
    (process.elecPtEtaIDIso*process.elecPtEtaIDIsoCounter) *
    process.tauPtEtaIDAgMuAgElecIsoPtRel*
    process.elecPtEtaRelID *
    process.metRecoilCorrector*
    (process.rescaledMETRawtau+process.rescaledTaus)*
    (process.tauPtEtaIDAgMuAgElecIsoTauUp*process.tauPtEtaIDAgMuAgElecIsoTauUpCounter)*
    #process.pfMEtMVACov*
    process.diTauRawTauUp*process.selectedDiTauRawTauUp*process.selectedDiTauRawTauUpCounter*
    process.elecTauStreamAnalyzerRawTauUp
    )

process.seqRawTauDown = cms.Sequence(
    process.allEventsFilter*
    (process.elecPtEtaIDIso*process.elecPtEtaIDIsoCounter) *
    process.tauPtEtaIDAgMuAgElecIsoPtRel*
    process.elecPtEtaID*
    process.elecPtEtaRelID *
    process.metRecoilCorrector*
    (process.rescaledMETRawtau+process.rescaledTaus)*
    (process.tauPtEtaIDAgMuAgElecIsoTauDown*process.tauPtEtaIDAgMuAgElecIsoTauDownCounter)*
    #process.pfMEtMVACov*
    process.diTauRawTauDown*process.selectedDiTauRawTauDown*process.selectedDiTauRawTauDownCounter*
    process.elecTauStreamAnalyzerRawTauDown
    )

#######################################################################

if runOnMC:
    #process.pNominal            = cms.Path( process.seqNominal )
    #process.pJetUp                 = cms.Path( process.seqJetUp   )
    #process.pJetDown               = cms.Path( process.seqJetDown )
    #process.pMEtResolutionUp       = cms.Path( process.seqMEtResolutionUp )
    #process.pMEtResolutionDown     = cms.Path( process.seqMEtResolutionDown )
    #process.pMEtResponseUp         = cms.Path( process.seqMEtResponseUp)
    #process.pMEtResponseDown       = cms.Path( process.seqMEtResponseDown)
    #process.pElecUp                  = cms.Path( process.seqElecUp)
    #process.pElecDown                = cms.Path( process.seqElecDown)
    #process.pTauUp              = cms.Path( process.seqTauUp)
    #process.pTauDown            = cms.Path( process.seqTauDown )
    process.pRawNominal         = cms.Path( process.seqRawNominal )
    #process.pRawJetUp              = cms.Path( process.seqRawJetUp   )
    #process.pRawJetDown            = cms.Path( process.seqRawJetDown )
    #process.pRawMEtResolutionUp    = cms.Path( process.seqRawMEtResolutionUp )
    #process.pRawMEtResolutionDown  = cms.Path( process.seqRawMEtResolutionDown )
    #process.pRawMEtResponseUp      = cms.Path( process.seqRawMEtResponseUp)
    #process.pRawMEtResponseDown    = cms.Path( process.seqRawMEtResponseDown)
    #process.pRawElecUp               = cms.Path( process.seqRawElecUp)
    #process.pRawElecDown             = cms.Path( process.seqRawElecDown)
    process.pRawTauUp           = cms.Path( process.seqRawTauUp )
    process.pRawTauDown         = cms.Path( process.seqRawTauDown )

else:
    #process.pNominal            = cms.Path( process.seqNominal )
    #process.pTauUp              = cms.Path( process.seqTauUp)
    #process.pTauDown            = cms.Path( process.seqTauDown )
    process.pRawNominal         = cms.Path( process.seqRawNominal )
    process.pRawTauUp           = cms.Path( process.seqRawTauUp )
    process.pRawTauDown         = cms.Path( process.seqRawTauDown )

#######################################################################

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

processDumpFile = open('runElecTauStreamAnalyzerFullAnalysis_Recoil.dump', 'w')
print >> processDumpFile, process.dumpPython()

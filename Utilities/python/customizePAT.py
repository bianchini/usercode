import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PFCandProducer.Isolation.tools_cfi import *

###################a#################################################
def addSelectedPFlowParticle(process,verbose=False):
    if verbose:
        print "[Info] Adding pf-particles (for pf-isolation and pf-seed pat-leptons)"
    process.load("PhysicsTools.PFCandProducer.ParticleSelectors.pfSortByType_cff")
    process.load("PhysicsTools.PFCandProducer.pfNoPileUp_cff")
    process.pfCandidateSelectionByType = cms.Sequence(
        process.pfNoPileUpSequence *
        ( process.pfAllNeutralHadrons +
          process.pfAllChargedHadrons +
          process.pfAllPhotons
          )  +
        process.pfAllMuons +
        process.pfAllElectrons
        )
    process.pfPileUp.Enable = True # enable pile-up filtering
    process.pfPileUp.Vertices = "offlinePrimaryVertices" # use vertices w/o BS
    process.pfAllMuons.src = "particleFlow"
    process.pfAllElectrons.src = "particleFlow"
    
    #process.patDefaultSequence.replace(process.patCandidates,
    #                                   process.pfCandidateSelectionByType+
    #                                   process.patCandidates)
    process.patDefaultSequence.replace(process.electronMatch,
                                       process.pfCandidateSelectionByType+
                                       process.electronMatch)
    

###################a#################################################
def addPFMuonIsolation(process,module,postfix="",verbose=False):
    if verbose:
        print "[Info] Adding particle isolation to muon with postfix '"+postfix+"'"

    if not hasattr(process, "pfCandidateSelectionByType"):
        addSelectedPFlowParticle(process,verbose=verbose)
        
    #setup correct src of isolated object
    setattr(process,"isoDepMuonWithCharged"+postfix,
            isoDepositReplace(module.muonSource,
                              'pfAllChargedHadrons'))
    setattr(process,"isoDepMuonWithNeutral"+postfix,
            isoDepositReplace(module.muonSource,
                              'pfAllNeutralHadrons'))
    setattr(process,"isoDepMuonWithPhotons"+postfix,
            isoDepositReplace(module.muonSource,
                              'pfAllPhotons'))

    #compute isolation values form deposits
    process.load("PhysicsTools.PFCandProducer.Isolation.pfMuonIsolationFromDeposits_cff")
    if postfix!="":
        setattr(process,"isoValMuonWithCharged"+postfix,
                process.isoValMuonWithCharged.clone())
        getattr(process,"isoValMuonWithCharged"+postfix).deposits.src="isoDepMuonWithCharged"+postfix
        setattr(process,"isoValMuonWithNeutral"+postfix,
                process.isoValMuonWithNeutral.clone())
        getattr(process,"isoValMuonWithNeutral"+postfix).deposits.src="isoDepMuonWithNeutral"+postfix
        setattr(process,"isoValMuonWithPhotons"+postfix,
                process.isoValMuonWithPhotons.clone())
        getattr(process,"isoValMuonWithPhotons"+postfix).deposits.src="isoDepMuonWithPhotons"+postfix
        
    setattr(process,"patMuonIsolationFromDepositsSequence"+postfix,
            cms.Sequence(getattr(process,"isoValMuonWithCharged"+postfix) +
                         getattr(process,"isoValMuonWithNeutral"+postfix) +
                         getattr(process,"isoValMuonWithPhotons"+postfix)
                         )
            )

    setattr(process,"patMuonIsoDepositsSequence"+postfix,
            cms.Sequence(getattr(process,"isoDepMuonWithCharged"+postfix) +
                         getattr(process,"isoDepMuonWithNeutral"+postfix) +
                         getattr(process,"isoDepMuonWithPhotons"+postfix)
                         )
            )
    setattr(process,"patMuonIsolationSequence"+postfix,
            cms.Sequence(getattr(process,"patMuonIsoDepositsSequence"+postfix) +
                         getattr(process,"patMuonIsolationFromDepositsSequence"+postfix)
                         )
            )
    
    getattr(process,"isoDepMuonWithCharged"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadrons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
        )
    getattr(process,"isoDepMuonWithNeutral"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadrons"),
        DR_Veto = cms.double(0.08),
        DepositLabel = cms.untracked.string('')
        )
    getattr(process,"isoDepMuonWithPhotons"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllPhotons"),
        DR_Veto = cms.double(0.05),
        DepositLabel = cms.untracked.string('')
        )
    
    getattr(process,"isoValMuonWithPhotons"+postfix).deposits = cms.VPSet(
        cms.PSet(
        src = cms.InputTag("isoDepMuonWithPhotons"+postfix),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(1.0)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
        )
        )
    getattr(process,"isoValMuonWithNeutral"+postfix).deposits = cms.VPSet(
        cms.PSet(
        src = cms.InputTag("isoDepMuonWithNeutral"+postfix),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(1.0)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
        )
        )
    getattr(process,"isoValMuonWithCharged"+postfix).deposits = cms.VPSet(
        cms.PSet(
        src = cms.InputTag("isoDepMuonWithCharged"+postfix),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(0.5)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
        )
        )

    

    module.isoDeposits = cms.PSet(
        pfChargedHadrons = cms.InputTag("isoDepMuonWithCharged"+postfix),
        pfNeutralHadrons = cms.InputTag("isoDepMuonWithNeutral"+postfix),
        pfPhotons = cms.InputTag("isoDepMuonWithPhotons"+postfix)
        )
    module.isolationValues = cms.PSet(
        pfChargedHadrons = cms.InputTag("isoValMuonWithCharged"+postfix),
        pfNeutralHadrons = cms.InputTag("isoValMuonWithNeutral"+postfix),
        pfPhotons = cms.InputTag("isoValMuonWithPhotons"+postfix)
        )
    
    process.patDefaultSequence.replace(module,
                                       getattr(process,"patMuonIsolationSequence"+postfix)+
                                       module
                                       )

####################################################################
def addPFMuon(process,postfix,verbose=False):
    if verbose:
        print "[Info] Adding pf-muon collection with postfix '"+postfix+"'"
    cloneProcessingSnippet(process, process.makePatMuons, postfix)
    getattr(process,"patMuons"+postfix).useParticleFlow = True
    pfSelectedMuons = cms.EDFilter( #dummy selector
        "PtMinPFCandidateSelector",
        src = cms.InputTag('pfAllMuons'),
        ptMin = cms.double(0.0)
        )
    setattr(process,"pfSelectedMuons"+postfix,pfSelectedMuons.clone())
    getattr(process,"patMuons"+postfix).pfMuonSource = "pfSelectedMuons"+postfix
    if  hasattr(process,"muonMatch"+postfix):
        getattr(process,"muonMatch"+postfix).src = getattr(process,"patMuons"+postfix).pfMuonSource
    ## Add pfMuons to sequence
    process.patCandidates.replace(process.makePatMuons,
                                  process.makePatMuons+
                                  getattr(process,"pfSelectedMuons"+postfix)*
                                  getattr(process,"makePatMuons"+postfix)
        )
    process.patCandidateSummary.candidates.append(cms.InputTag("patMuons"+postfix))
    # check if there is pf-isolation, if not add it
    if not hasattr(process, "patMuonIsolationFromDepositsSequence"+postfix):
        addPFMuonIsolation(process,getattr(process,"patMuons"+postfix),
                           postfix=postfix,verbose=verbose)
    #setup the isolation
    getattr(process,"isoDepMuonWithCharged"+postfix).src = getattr(process,"pfSelectedMuons"+postfix).src
    getattr(process,"isoDepMuonWithNeutral"+postfix).src = getattr(process,"pfSelectedMuons"+postfix).src
    getattr(process,"isoDepMuonWithPhotons"+postfix).src = getattr(process,"pfSelectedMuons"+postfix).src
    # and now selected Muons
    setattr(process,"selectedPatMuons"+postfix,process.selectedPatMuons.clone())
    getattr(process,"selectedPatMuons"+postfix).src = 'patMuons'+postfix
    process.selectedPatCandidates.replace(process.selectedPatMuons,
                                          process.selectedPatMuons+
                                          getattr(process,"selectedPatMuons"+postfix)
                                          )
    process.selectedPatCandidateSummary.candidates.append(cms.InputTag("selectedPatMuons"+postfix))
    # and add counter
    setattr(process,"countPatMuons"+postfix,process.countPatMuons.clone())
    getattr(process,"countPatMuons"+postfix).src = 'selectedPatMuons'+postfix
    process.countPatCandidates.replace(process.countPatMuons,
                                       process.countPatMuons+
                                       getattr(process,"countPatMuons"+postfix)
                                       )

###################a#################################################
def addPFElectronIsolation(process,module,postfix="",verbose=False):
    if verbose:
        print "[Info] Adding particle isolation to electron with postfix '"+postfix+"'"

    if not hasattr(process, "pfCandidateSelectionByType"):
        addSelectedPFlowParticle(process,verbose=verbose)
        
    #setup correct src of isolated object
    setattr(process,"isoDepElectronWithCharged"+postfix,
            isoDepositReplace(module.electronSource,
                              'pfAllChargedHadrons'))
    setattr(process,"isoDepElectronWithNeutral"+postfix,
            isoDepositReplace(module.electronSource,
                              'pfAllNeutralHadrons'))
    setattr(process,"isoDepElectronWithPhotons"+postfix,
            isoDepositReplace(module.electronSource,
                              'pfAllPhotons'))

    #compute isolation values form deposits
    process.load("PhysicsTools.PFCandProducer.Isolation.pfElectronIsolationFromDeposits_cff")
    if postfix!="":
        setattr(process,"isoValElectronWithCharged"+postfix,
                process.isoValElectronWithCharged.clone())
        getattr(process,"isoValElectronWithCharged"+postfix).deposits.src="isoDepElectronWithCharged"+postfix
        setattr(process,"isoValElectronWithNeutral"+postfix,
                process.isoValElectronWithNeutral.clone())
        getattr(process,"isoValElectronWithNeutral"+postfix).deposits.src="isoDepElectronWithNeutral"+postfix
        setattr(process,"isoValElectronWithPhotons"+postfix,
                process.isoValElectronWithPhotons.clone())
        getattr(process,"isoValElectronWithPhotons"+postfix).deposits.src="isoDepElectronWithPhotons"+postfix
        
    setattr(process,"patElectronIsolationFromDepositsSequence"+postfix,
            cms.Sequence(getattr(process,"isoValElectronWithCharged"+postfix) +
                         getattr(process,"isoValElectronWithNeutral"+postfix) +
                         getattr(process,"isoValElectronWithPhotons"+postfix)
                         )
            )

    setattr(process,"patElectronIsoDepositsSequence"+postfix,
            cms.Sequence(getattr(process,"isoDepElectronWithCharged"+postfix) +
                         getattr(process,"isoDepElectronWithNeutral"+postfix) +
                         getattr(process,"isoDepElectronWithPhotons"+postfix)
                         )
            )
    setattr(process,"patElectronIsolationSequence"+postfix,
            cms.Sequence(getattr(process,"patElectronIsoDepositsSequence"+postfix) +
                         getattr(process,"patElectronIsolationFromDepositsSequence"+postfix)
                         )
            )

    getattr(process,"isoDepElectronWithCharged"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadrons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
        )
    getattr(process,"isoDepElectronWithNeutral"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadrons"),
        DR_Veto = cms.double(0.08),
        DepositLabel = cms.untracked.string('')
        )
    getattr(process,"isoDepElectronWithPhotons"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllPhotons"),
        DR_Veto = cms.double(0.05),
        DepositLabel = cms.untracked.string('')
        )
    
    getattr(process,"isoValElectronWithPhotons"+postfix).deposits = cms.VPSet(
        cms.PSet(
        src = cms.InputTag("isoDepElectronWithPhotons"+postfix),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(1.0)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
        )
        )
    getattr(process,"isoValElectronWithNeutral"+postfix).deposits = cms.VPSet(
        cms.PSet(
        src = cms.InputTag("isoDepElectronWithNeutral"+postfix),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(1.0)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
        )
        )
    getattr(process,"isoValElectronWithCharged"+postfix).deposits = cms.VPSet(
        cms.PSet(
        src = cms.InputTag("isoDepElectronWithCharged"+postfix),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(0.5)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
        )
        )

    module.isoDeposits = cms.PSet(
        pfChargedHadrons = cms.InputTag("isoDepElectronWithCharged"+postfix),
        pfNeutralHadrons = cms.InputTag("isoDepElectronWithNeutral"+postfix),
        pfPhotons = cms.InputTag("isoDepElectronWithPhotons"+postfix)
        )
    module.isolationValues = cms.PSet(
        pfChargedHadrons = cms.InputTag("isoValElectronWithCharged"+postfix),
        pfNeutralHadrons = cms.InputTag("isoValElectronWithNeutral"+postfix),
        pfPhotons = cms.InputTag("isoValElectronWithPhotons"+postfix)
        )
    
    process.patDefaultSequence.replace(module,
                                       getattr(process,"patElectronIsolationSequence"+postfix)+
                                       module
                                       )

####################################################################
def addPFElectron(process,postfix,verbose=False):
    if verbose:
        print "[Info] Adding pf-electron collection with postfix '"+postfix+"'"
    cloneProcessingSnippet(process, process.makePatElectrons, postfix)
    getattr(process,"patElectrons"+postfix).useParticleFlow = True
    pfSelectedElectrons = cms.EDFilter( #dummy selector
        "PtMinPFCandidateSelector",
        src = cms.InputTag('pfAllElectrons'),
        ptMin = cms.double(0.0)
        )
    setattr(process,"pfSelectedElectrons"+postfix,pfSelectedElectrons.clone())
    getattr(process,"patElectrons"+postfix).pfElectronSource = "pfSelectedElectrons"+postfix
    ## Add pfElectrons to sequence
    process.patCandidates.replace(process.makePatElectrons,
                                  process.makePatElectrons+
                                  getattr(process,"pfSelectedElectrons"+postfix)*
                                  getattr(process,"makePatElectrons"+postfix)
        )
    process.patCandidateSummary.candidates.append(cms.InputTag("patElectrons"+postfix))
    # check if there is pf-isolation, if not add it
    if not hasattr(process, "patElectronIsolationFromDepositsSequence"+postfix):
        addPFElectronIsolation(process,getattr(process,"patElectrons"+postfix),
                           postfix=postfix,verbose=verbose)
    #setup the isolation
    getattr(process,"isoDepElectronWithCharged"+postfix).src = getattr(process,"pfSelectedElectrons"+postfix).src
    getattr(process,"isoDepElectronWithNeutral"+postfix).src = getattr(process,"pfSelectedElectrons"+postfix).src
    getattr(process,"isoDepElectronWithPhotons"+postfix).src = getattr(process,"pfSelectedElectrons"+postfix).src
    # and now selected Electrons
    setattr(process,"selectedPatElectrons"+postfix,process.selectedPatElectrons.clone())
    getattr(process,"selectedPatElectrons"+postfix).src = 'patElectrons'+postfix
    process.selectedPatCandidates.replace(process.selectedPatElectrons,
                                          process.selectedPatElectrons+
                                          getattr(process,"selectedPatElectrons"+postfix)
                                          )
    process.selectedPatCandidateSummary.candidates.append(cms.InputTag("selectedPatElectrons"+postfix))
    # and add counter
    setattr(process,"countPatElectrons"+postfix,process.countPatElectrons.clone())
    getattr(process,"countPatElectrons"+postfix).src = 'selectedPatElectrons'+postfix
    process.countPatCandidates.replace(process.countPatElectrons,
                                       process.countPatElectrons+
                                       getattr(process,"countPatElectrons"+postfix)
                                       )


    
####################################################################
def addTriggerMatchingMuon(process,postfix="",verbose=False):
    if verbose:
        print "[Info] Addting trigger matching to patMuon with postfix '"+postfix+"'"
    
    process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi" )
    muonMatch = cms.EDProducer(
        "PATTriggerMatcherDRLessByR"
        , src     = cms.InputTag( "selectedPatMuons"+postfix )
        , matched = cms.InputTag( "patTrigger" )
        #, matchedCuts    = cms.string('path("HLT_Mu9") || path("HLT_Mu11") || path("HLT_IsoMu9") || path("HLT_Mu15_v*")')
        , matchedCuts    = cms.string('path("HLT_Mu11_PFTau15_v*",0) && type("TriggerMuon")')
        , maxDPtRel = cms.double( 999 )
        , maxDeltaR = cms.double( 0.5 )
        , resolveAmbiguities    = cms.bool( True )
        , resolveByMatchQuality = cms.bool( True )
        )

    setattr(process,"muonTriggerMatchHLTMuons"+postfix,muonMatch.clone()) 
    if not hasattr(process,"patTriggerSequence"):
        setattr(process,"patTriggerSequence",
                cms.Sequence(process.patTrigger))
        process.patDefaultSequence += process.patTriggerSequence

    if not hasattr(process,"patTriggerMatcher"):
        setattr(process,"patTriggerMatcher",
                cms.Sequence(getattr(process,"muonTriggerMatchHLTMuons"+postfix)))
        process.patTriggerSequence += process.patTriggerMatcher
    else:
        process.patTriggerMatcher += getattr(process,"muonTriggerMatchHLTMuons"+postfix)

    setattr(process,"selectedPatMuonsTriggerMatch"+postfix,
            cms.EDProducer("PATTriggerMatchMuonEmbedder",
                           src = cms.InputTag( "selectedPatMuons"+postfix ),
                           matches = cms.VInputTag("muonTriggerMatchHLTMuons"+postfix)
                           ))
    if not hasattr(process,"patTriggerMatchEmbedder"):
        setattr(process,"patTriggerMatchEmbedder",
                cms.Sequence(getattr(process,"selectedPatMuonsTriggerMatch"+postfix)))
        process.patTriggerSequence += process.patTriggerMatchEmbedder
    else:
        process.patTriggerMatchEmbedder += getattr(process,"selectedPatMuonsTriggerMatch"+postfix)
    

    process.patTrigger.onlyStandAlone = True


####################################################################
def addTriggerMatchingElectron(process,postfix="",verbose=False):
    if verbose:
        print "[Info] Addting trigger matching to patElectron with postfix '"+postfix+"'"
    
    process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi" )
    eleMatch = cms.EDProducer(
        "PATTriggerMatcherDRLessByR"
        , src     = cms.InputTag( "selectedPatElectrons"+postfix )
        , matched = cms.InputTag( "patTrigger" )
        #, matchedCuts = cms.string( 'path("HLT_Ele10_SW_L1R_v*") || path("HLT_Ele15_SW_L1R") || path("HLT_Ele17_SW_L1R_v*") || path("HLT_Ele17_SW_Isol_L1R_v*") || path("HLT_Ele17_SW_TighterEleIdIsol_L1R_v*") || path("HLT_Ele22_SW_L1R_v*")' )
        , matchedCuts = cms.string( 'path("HLT_IsoEle12_PFTau15_v*",0) && type("TriggerElectron")' )
        , maxDPtRel = cms.double( 999 )
        , maxDeltaR = cms.double( 0.5)
        , resolveAmbiguities    = cms.bool( True )
        , resolveByMatchQuality = cms.bool( True )
        )

    setattr(process,"eleTriggerMatchHLTElectrons"+postfix,eleMatch.clone()) 
    if not hasattr(process,"patTriggerSequence"):
        setattr(process,"patTriggerSequence",
                cms.Sequence(process.patTrigger))
        process.patDefaultSequence += process.patTriggerSequence

    if not hasattr(process,"patTriggerMatcher"):
        setattr(process,"patTriggerMatcher",
                cms.Sequence(getattr(process,"eleTriggerMatchHLTElectrons"+postfix)))
        process.patTriggerSequence += process.patTriggerMatcher
    else:
        process.patTriggerMatcher += getattr(process,"eleTriggerMatchHLTElectrons"+postfix)

    setattr(process,"selectedPatElectronsTriggerMatch"+postfix,
            cms.EDProducer("PATTriggerMatchElectronEmbedder",
                           src = cms.InputTag( "selectedPatElectrons"+postfix ),
                           matches = cms.VInputTag("eleTriggerMatchHLTElectrons"+postfix)
                           ))
    if not hasattr(process,"patTriggerMatchEmbedder"):
        setattr(process,"patTriggerMatchEmbedder",
                cms.Sequence(getattr(process,"selectedPatElectronsTriggerMatch"+postfix)))
        process.patTriggerSequence += process.patTriggerMatchEmbedder
    else:
        process.patTriggerMatchEmbedder += getattr(process,"selectedPatElectronsTriggerMatch"+postfix)
    

    process.patTrigger.onlyStandAlone = True


####################################################################


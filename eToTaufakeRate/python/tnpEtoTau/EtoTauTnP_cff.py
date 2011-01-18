import FWCore.ParameterSet.Config as cms


# make the tag
tagAnyEle = cms.EDFilter("PATElectronRefSelector",
                         src = cms.InputTag("selectedPatElectrons"),
                         cut = cms.string('pt>25.0 && abs(eta)<2.5 && !(1.4442< abs(eta) <1.566)'),
                         filter = cms.bool(False)
                         )

tagRS = tagAnyEle.clone(
    cut = cms.string( tagAnyEle.cut.value()+' && (( isEB && ( (dr03TkSumPt + max(0., dr03EcalRecHitSumEt - 1.) + dr03HcalTowerSumEt)/(p4.Pt) < 0.15 )) || (isEE && ((dr03TkSumPt + dr03EcalRecHitSumEt + dr03HcalTowerSumEt)/(p4.Pt) < 0.1 )))')
    )
tag80 = tagAnyEle.clone(
    cut = cms.string(tagAnyEle.cut.value()+' && electronID("simpleEleId80relIso") > 6.5'),
    )
tag90 = tagAnyEle.clone(
    cut = cms.string(tagAnyEle.cut.value()+' && electronID("simpleEleId90relIso") > 6.5'),
    )


## make the probe

#HPS Tau
probeAnyTau = cms.EDFilter("PATTauRefSelector",
                                   src = cms.InputTag("selectedPatTausPFlow"),
                                   cut = cms.string('pt>15.0 && abs(eta)<2.5'),
                                   filter = cms.bool(False)
                                   )
probeTauID = probeAnyTau.clone(
    cut = cms.string(probeAnyTau.cut.value()+' && tauID("leadingTrackFinding") > 0.5')
    )
probeTauMediumIso = probeAnyTau.clone(
    cut = cms.string(probeAnyTau.cut.value()+' && tauID("byMediumIsolation") > 0.5')
    )
probeTauLooseIso = probeAnyTau.clone(
    cut = cms.string(probeAnyTau.cut.value()+' && tauID("byLooseIsolation") > 0.5')
    )
probeTauTightIso = probeAnyTau.clone(
    cut = cms.string(probeAnyTau.cut.value()+' && tauID("byTightIsolation") > 0.5')
    )
probeTauIDIsoMedium = probeAnyTau.clone(
    cut = cms.string(probeTauID.cut.value()+' && '+probeTauMediumIso.cut.value())
    )
probeTauIDIsoLoose = probeAnyTau.clone(
    cut = cms.string(probeTauID.cut.value()+' && '+probeTauLooseIso.cut.value())
    )
probeTauIDIsoTight = probeAnyTau.clone(
    cut = cms.string(probeTauID.cut.value()+' && '+probeTauTightIso.cut.value())
    )

## Shrinking Cone
probeAnyTauSC = probeAnyTau.clone(
    src = cms.InputTag("selectedPatTausPFlowShrinkingCone"),
    )
probeTauSCID = probeAnyTauSC.clone(    
    cut = cms.string(probeAnyTau.cut.value()+' && tauID("leadingTrackFinding") > 0.5 && tauID("leadingTrackPtCut") > 0.5')
    )
probeTauSCIso = probeAnyTauSC.clone(    
    cut = cms.string(probeAnyTau.cut.value()+' && tauID("byIsolation") > 0.5')
    )
probeTauSCIDIso = probeAnyTauSC.clone(    
    cut = cms.string(probeAnyTau.cut.value()+' && '+probeTauSCID.cut.value()+ ' && ' +probeTauSCIso.cut.value())
    )




tnpAnyEleAnyTau = cms.EDProducer("CandViewShallowCloneCombiner",
                                 decay = cms.string("tagAnyEle@+ probeAnyTau@-"),
                                 cut   = cms.string("60 < mass < 120"),
                                 #checkCharge = cms.bool(False),   
                                 )
tnpRSAnyTau = tnpAnyEleAnyTau.clone(decay = cms.string("tagRS@+ probeAnyTau@-"))
tnp80AnyTau = tnpAnyEleAnyTau.clone(decay = cms.string("tag80@+ probeAnyTau@-"))
tnp90AnyTau = tnpAnyEleAnyTau.clone(decay = cms.string("tag90@+ probeAnyTau@-"))
tnp80TauIDLoose = tnpAnyEleAnyTau.clone(decay = cms.string("tag80@+ probeTauIDIsoLoose@-"))
tnp80TauIDMedium = tnpAnyEleAnyTau.clone(decay = cms.string("tag80@+ probeTauIDIsoMedium@-"))
tnp80TauIDTight = tnpAnyEleAnyTau.clone(decay = cms.string("tag80@+ probeTauIDIsoTight@-"))
tnp90TauIDLoose = tnpAnyEleAnyTau.clone(decay = cms.string("tag90@+ probeTauIDIsoLoose@-"))
tnp90TauIDMedium = tnpAnyEleAnyTau.clone(decay = cms.string("tag90@+ probeTauIDIsoMedium@-"))
tnp90TauIDTight = tnpAnyEleAnyTau.clone(decay = cms.string("tag90@+ probeTauIDIsoTight@-"))

tnpAnyEleAnyTauSC = tnpAnyEleAnyTau.clone(
    decay = cms.string("tagAnyEle@+ probeAnyTauSC@-")
    )
tnp90AnyTauSC = tnpAnyEleAnyTauSC.clone(
    decay = cms.string("tag90@+ probeAnyTauSC@-")
    )
tnp80AnyTauSC = tnpAnyEleAnyTauSC.clone(
    decay = cms.string("tag80@+ probeAnyTauSC@-")
    )
tnp90TauSCIDIso = tnpAnyEleAnyTauSC.clone(
    decay = cms.string("tag90@+ probeTauSCIDIso@-")
    )
tnp80TauSCIDIso = tnpAnyEleAnyTauSC.clone(
    decay = cms.string("tag80@+ probeTauSCIDIso@-")
    )


eleMcMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
                            pdgId = cms.vint32(11,-11),
                            src = cms.InputTag("selectedPatElectrons"),
                            distMin = cms.double(0.3),
                            matched = cms.InputTag("genParticles")
                            )

eleMcMatchTag80   = eleMcMatch.clone(src = cms.InputTag("tag80"))
eleMcMatchTag90   = eleMcMatch.clone(src = cms.InputTag("tag90"))

eleMcMatchProbeAnyTau = eleMcMatch.clone(src = cms.InputTag("probeAnyTau"))
eleMcMatchProbeTauIDIsoMedium = eleMcMatch.clone(src = cms.InputTag("probeTauIDIsoMedium"))
eleMcMatchProbeTauIDIsoLoose = eleMcMatch.clone(src = cms.InputTag("probeTauIDIsoLoose"))
eleMcMatchProbeTauIDIsoTight = eleMcMatch.clone(src = cms.InputTag("probeTauIDIsoTight"))

eleMcMatchProbeAnyTauSC = eleMcMatch.clone(src = cms.InputTag("probeAnyTauSC"))
eleMcMatchProbeTauSCIDIso = eleMcMatch.clone(src = cms.InputTag("probeTauSCIDIso"))


etoTauRS = cms.EDAnalyzer("TagProbeFitTreeProducer",
                          tagProbePairs = cms.InputTag("tnpRSAnyTau"),
                          arbitration   = cms.string("BestMass"),
                          massForArbitration = cms.double(91),
                          variables = cms.PSet(
    pt = cms.string("pt"), 
    abseta = cms.string("abs(eta)"),
    electronPreIDOutput = cms.string("electronPreIDOutput"),
    hcal3x3OverPLead = cms.string("hcal3x3OverPLead"),
    ecalStripSumEOverPLead = cms.string("ecalStripSumEOverPLead"),
    ),
                          flags = cms.PSet(
    tauID = cms.string('tauID("leadingTrackFinding")>0.5'),
    tauLooseIso =  cms.string('tauID("byLooseIsolation")>0.5'),
    tauMediumIso =  cms.string('tauID("byMediumIsolation")>0.5'),
    tauTightIso =  cms.string('tauID("byTightIsolation")>0.5'),
    tauAntiEMVA =  cms.string('tauID("againstElectron")>0.5'),
    tauAntiECrackRem =  cms.string('tauID("againstElectronCrackRem")>0.5'),
    tauAntiE2D =  cms.string('tauID("againstElectron2D")>0.5'),
    tauAntiEMVAandCrackRem =  cms.string('tauID("againstElectron")>0.5 && tauID("againstElectronCrackRem")>0.5'),
    tauAntiE2DandCrackRem =  cms.string('tauID("againstElectron2D")>0.5 && tauID("againstElectronCrackRem")>0.5')
    ),
                          isMC = cms.bool( True ),
                          tagMatches = cms.InputTag("eleMcMatchTag"),
                          probeMatches  = cms.InputTag("eleMcMatchProbeAnyTau"),
                          motherPdgId = cms.vint32(22,23),
                          makeMCUnbiasTree = cms.bool(True),
                          checkMotherInUnbiasEff = cms.bool(True),
                          allProbes = cms.InputTag("probeAnyTau"),
                          addRunLumiInfo = cms.bool(True)
                          )

etoTauSC = etoTauRS.clone(
    tagProbePairs = cms.InputTag("tnpAnyEleAnyTauSC"),
    flags = cms.PSet(
    tauID = cms.string('tauID("leadingTrackFinding")>0.5 && tauID("leadingTrackPtCut") > 0.5'),
    tauIso =  cms.string('tauID("byIsolation")>0.5'),
    tauAntiEMVA =  cms.string('tauID("againstElectron")>0.5'),
    tauAntiECrackRem =  cms.string('tauID("againstElectronCrackRem")>0.5'),
    tauAntiE2D =  cms.string('tauID("againstElectron2D")>0.5'),
    tauAntiEMVAandCrackRem =  cms.string('tauID("againstElectron")>0.5 && tauID("againstElectronCrackRem")>0.5'),
    tauAntiE2DandCrackRem =  cms.string('tauID("againstElectron2D")>0.5 && tauID("againstElectronCrackRem")>0.5')
    ),
    probeMatches  = cms.InputTag("eleMcMatchProbeAnyTauSC"),
    allProbes = cms.InputTag("probeAnyTauSC")
    )

etoTau80 = etoTauRS.clone(
    tagProbePairs = cms.InputTag("tnp80AnyTau"),
    tagMatches = cms.InputTag("eleMcMatchTag80"),
    #probeMatches  = cms.InputTag("eleMcMatchProbeAnyTau"),
)
etoTau90 = etoTauRS.clone(
    tagProbePairs = cms.InputTag("tnp90AnyTau"),
    tagMatches = cms.InputTag("eleMcMatchTag90"),
    #probeMatches  = cms.InputTag("eleMcMatchProbeAnyTau"),
    )

etoTauSC80 = etoTauSC.clone(
    tagProbePairs = cms.InputTag("tnp80AnyTauSC"),
    tagMatches = cms.InputTag("eleMcMatchTag80"),
    #probeMatches  = cms.InputTag("eleMcMatchProbeAnyTauSC"),
)
etoTauSC90 = etoTauSC.clone(
    tagProbePairs = cms.InputTag("tnp90AnyTauSC"),
    tagMatches = cms.InputTag("eleMcMatchTag90"),
    #probeMatches  = cms.InputTag("eleMcMatchProbeAnyTauSC"),
)


etoTauMargMedium80 = etoTauRS.clone(
    tagProbePairs = cms.InputTag("tnp80TauIDMedium"),
    flags = cms.PSet(
    tauAntiEMVA =  cms.string('tauID("againstElectron")>0.5'),
    tauAntiECrackRem =  cms.string('tauID("againstElectronCrackRem")>0.5'),
    tauAntiE2D =  cms.string('tauID("againstElectron2D")>0.5'),
    tauAntiEMVAandCrackRem =  cms.string('tauID("againstElectron")>0.5 && tauID("againstElectronCrackRem")>0.5'),
    tauAntiE2DandCrackRem =  cms.string('tauID("againstElectron2D")>0.5 && tauID("againstElectronCrackRem")>0.5')
    ),
    tagMatches = cms.InputTag("eleMcMatchTag80"),
    probeMatches  = cms.InputTag("eleMcMatchProbeTauIDIsoMedium"),
    makeMCUnbiasTree = cms.bool(True),
    checkMotherInUnbiasEff = cms.bool(True),
    allProbes = cms.InputTag("probeTauIDIsoMedium")
    )

etoTauMargLoose80 = etoTauMargMedium80.clone(
    tagProbePairs = cms.InputTag("tnp80TauIDLoose"),
    tagMatches = cms.InputTag("eleMcMatchTag80"),
    probeMatches  = cms.InputTag("eleMcMatchProbeTauIDIsoLoose"),
    allProbes = cms.InputTag("probeTauIDIsoLoose")
    )
etoTauMargTight80 = etoTauMargMedium80.clone(
    tagProbePairs = cms.InputTag("tnp80TauIDTight"),
    tagMatches = cms.InputTag("eleMcMatchTag80"),
    probeMatches  = cms.InputTag("eleMcMatchProbeTauIDIsoTight"),
    allProbes = cms.InputTag("probeTauIDIsoTight")
    )
etoTauMargMedium90 = etoTauMargMedium80.clone(
    tagProbePairs = cms.InputTag("tnp90TauIDMedium"),
    tagMatches = cms.InputTag("eleMcMatchTag90"),
    probeMatches  = cms.InputTag("eleMcMatchProbeTauIDIsoMedium"),
    allProbes = cms.InputTag("probeTauIDIsoMedium")
    )
etoTauMargLoose90 = etoTauMargMedium80.clone(
    tagProbePairs = cms.InputTag("tnp90TauIDLoose"),
    tagMatches = cms.InputTag("eleMcMatchTag90"),
    probeMatches  = cms.InputTag("eleMcMatchProbeTauIDIsoLoose"),
    allProbes = cms.InputTag("probeTauIDIsoLoose")
    )
etoTauMargTight90 = etoTauMargMedium80.clone(
    tagProbePairs = cms.InputTag("tnp90TauIDTight"),
    tagMatches = cms.InputTag("eleMcMatchTag90"),
    probeMatches  = cms.InputTag("eleMcMatchProbeTauIDIsoTight"),
    allProbes = cms.InputTag("probeTauIDIsoTight")
    )

etoTauSCMarg80 = etoTauSC.clone(
    tagProbePairs = cms.InputTag("tnp80TauSCIDIso"),
    flags = cms.PSet(
    tauAntiEMVA =  cms.string('tauID("againstElectron")>0.5'),
    tauAntiECrackRem =  cms.string('tauID("againstElectronCrackRem")>0.5'),
    tauAntiE2D =  cms.string('tauID("againstElectron2D")>0.5'),
    tauAntiEMVAandCrackRem =  cms.string('tauID("againstElectron")>0.5 && tauID("againstElectronCrackRem")>0.5'),
    tauAntiE2DandCrackRem =  cms.string('tauID("againstElectron2D")>0.5 && tauID("againstElectronCrackRem")>0.5')
    ),
    tagMatches = cms.InputTag("eleMcMatchTag80"),
    probeMatches  = cms.InputTag("eleMcMatchProbeTauSCIDIso"),
    allProbes = cms.InputTag("probeTauSCIDIso")
    )

etoTauSCMarg90 = etoTauSCMarg80.clone(
    tagProbePairs = cms.InputTag("tnp90TauSCIDIso"),
    tagMatches = cms.InputTag("eleMcMatchTag90"),
    )


oneTp80 = cms.EDFilter("CandViewCountFilter",
                               src = cms.InputTag("tnp80AnyTau"),
                               minNumber = cms.uint32(1)
                               )

oneTp90 = oneTp80.clone( src = cms.InputTag("tnp90AnyTau") )

oneTpSC80 = oneTp80.clone( src = cms.InputTag("tnp80AnyTauSC"))
oneTpSC90 = oneTp80.clone( src = cms.InputTag("tnp90AnyTauSC"))
oneTp90Loose=oneTp80.clone(src=cms.InputTag("tnp90TauIDLoose"))
oneTp80Loose=oneTp80.clone(src=cms.InputTag("tnp80TauIDLoose"))
oneTpSC80Loose=oneTp80.clone(src=cms.InputTag("tnp80TauSCIDIso"))
oneTpSC90Loose=oneTp80.clone(src=cms.InputTag("tnp90TauSCIDIso"))

oneTp90 = cms.EDFilter("CandViewCountFilter",
                               src = cms.InputTag("tnp90AnyTau"),
                               minNumber = cms.uint32(1)
                               )
oneTpSC80 = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("tnp80AnyTauSC"),
                                 minNumber = cms.uint32(1)
                                 )
oneTpSC90 = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("tnp90AnyTauSC"),
                                 minNumber = cms.uint32(1)
                                 )


sequence80 = cms.Sequence((tag80+
                           probeAnyTau+
                           probeTauIDIsoMedium+
                           probeTauIDIsoLoose+
                           probeTauIDIsoTight)*
                          (eleMcMatchTag80+
                           eleMcMatchProbeAnyTau+
                           eleMcMatchProbeTauIDIsoMedium+
                           eleMcMatchProbeTauIDIsoLoose+
                           eleMcMatchProbeTauIDIsoTight)*
                          (tnp80AnyTau+
                           tnp80TauIDLoose+
                           tnp80TauIDMedium+
                           tnp80TauIDTight)*
                          oneTp80*
                          (etoTau80+
                           etoTauMargLoose80+
                           etoTauMargMedium80+
                           etoTauMargTight80)
                          )

sequence90 = cms.Sequence((tag90+
                           probeAnyTau+
                           probeTauIDIsoMedium+
                           probeTauIDIsoLoose+
                           probeTauIDIsoTight)*
                          (eleMcMatchTag90+
                           eleMcMatchProbeAnyTau+
                           eleMcMatchProbeTauIDIsoMedium+
                           eleMcMatchProbeTauIDIsoLoose+
                           eleMcMatchProbeTauIDIsoTight)*
                          (tnp90AnyTau+
                           tnp90TauIDLoose+
                           tnp90TauIDMedium+
                           tnp90TauIDTight)*
                          oneTp90*
                          (etoTau90+
                           etoTauMargLoose90+
                           etoTauMargMedium90+
                           etoTauMargTight90)
                          )

sequenceSC80 = cms.Sequence((tag80+
                             probeAnyTauSC+
                             probeTauSCIDIso)*
                            (eleMcMatchTag80+
                             eleMcMatchProbeAnyTauSC+
                             eleMcMatchProbeTauSCIDIso)*
                            (tnp80AnyTauSC+
                             tnp80TauSCIDIso)*
                            oneTpSC80*
                            (etoTauSC80+
                             etoTauSCMarg80)
                            )
sequenceSC90 = cms.Sequence((tag90+
                             probeAnyTauSC+
                             probeTauSCIDIso)*
                            (eleMcMatchTag90+
                             eleMcMatchProbeAnyTauSC+
                             eleMcMatchProbeTauSCIDIso)*
                            (tnp90AnyTauSC+
                             tnp90TauSCIDIso)*
                            oneTpSC90*
                            (etoTauSC90+
                             etoTauSCMarg90)
                            )

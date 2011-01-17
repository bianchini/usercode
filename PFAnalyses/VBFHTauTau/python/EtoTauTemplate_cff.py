import FWCore.ParameterSet.Config as cms

# make the tag
tagAnyEle = cms.EDFilter("PATElectronRefSelector",
                         src = cms.InputTag("patElectronsTriggerMatchEmbedder"),
                         cut = cms.string('pt>25.0 && abs(eta)<2.5 && !(1.4442< abs(eta) <1.566)'),
                         filter = cms.bool(False)
                         )

probeAnyTau = cms.EDFilter("PATTauRefSelector",
                           src = cms.InputTag("selectedPatTausPFlow"),
                           cut = cms.string('pt>15.0 && abs(eta)<2.5'),
                           filter = cms.bool(False)
                           )

tnpAnyEleAnyTau = cms.EDProducer("CandViewShallowCloneCombiner",
                                 decay = cms.string("tagAnyEle@+ probeAnyTau@-"),
                                 cut   = cms.string("40 < mass < 120"),
                                 #checkCharge = cms.bool(False),   
                                 )

oneTp = cms.EDFilter("CandViewCountFilter",
                     src = cms.InputTag("tnpAnyEleAnyTau"),
                     minNumber = cms.uint32(1)
                     )

tagMcMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
                            pdgId = cms.vint32(11,-11),
                            src = cms.InputTag("tagAnyEle"),
                            distMin = cms.double(0.3),
                            matched = cms.InputTag("genParticles")
                            )

probeMcMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
                              pdgId = cms.vint32(11,-11),
                              src = cms.InputTag("probeAnyTau"),
                              distMin = cms.double(0.3),
                              matched = cms.InputTag("genParticles")
                              )


etoTau = cms.EDAnalyzer("TagProbeFitTreeProducer",
                        tagProbePairs = cms.InputTag("tnpAnyEleAnyTau"),
                        arbitration   = cms.string("BestMass"),
                        massForArbitration = cms.double(91),
                        variables = cms.PSet(
    ),
                        flags = cms.PSet(
    ),
                        isMC = cms.bool( True ),
                        tagMatches = cms.InputTag("tagMcMatch") ,
                        probeMatches  = cms.InputTag("probeMcMatch"),
                        motherPdgId = cms.vint32(23),
                        makeMCUnbiasTree = cms.bool(True),
                        checkMotherInUnbiasEff = cms.bool(True),
                        allProbes = cms.InputTag("probeAnyTau"),
                        addRunLumiInfo = cms.bool(True)
                        )


sequence = cms.Sequence(
    (tagAnyEle+probeAnyTau)*
    tnpAnyEleAnyTau*
    oneTp*
    (tagMcMatch+probeMcMatch)*
    etoTau
    )


import FWCore.ParameterSet.Config as cms


eleMcMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
                            pdgId = cms.vint32(11,-11),
                            src = cms.InputTag("selectedPatElectronsTriggerMatch"),
                            distMin = cms.double(0.3),
                            matched = cms.InputTag("genParticles")
                            )

eleMcMatchTag90   = eleMcMatch.clone(src = cms.InputTag("tag90"))
eleMcMatchTag80   = eleMcMatch.clone(src = cms.InputTag("tag80"))
eleMcMatchTag70   = eleMcMatch.clone(src = cms.InputTag("tag70"))
eleMcMatchTag60   = eleMcMatch.clone(src = cms.InputTag("tag60"))

eleMcMatchProbeAnyTau = eleMcMatch.clone(src = cms.InputTag("probeAnyTau"))
eleMcMatchProbeTauIDIsoMedium = eleMcMatch.clone(src = cms.InputTag("probeTauIDIsoMedium"))
eleMcMatchProbeTauIDIsoLoose = eleMcMatch.clone(src = cms.InputTag("probeTauIDIsoLoose"))
eleMcMatchProbeTauIDIsoTight = eleMcMatch.clone(src = cms.InputTag("probeTauIDIsoTight"))
eleMcMatchProbeTauIDIsoMediumNoCracks = eleMcMatch.clone(src = cms.InputTag("probeTauIDIsoMediumNoCracks"))
eleMcMatchProbeTauIDIsoLooseNoCracks = eleMcMatch.clone(src = cms.InputTag("probeTauIDIsoLooseNoCracks"))
eleMcMatchProbeTauIDIsoTightNoCracks = eleMcMatch.clone(src = cms.InputTag("probeTauIDIsoTightNoCracks"))


eleMcMatchProbeAnyTauSC = eleMcMatch.clone(src = cms.InputTag("probeAnyTauSC"))
eleMcMatchProbeTauSCIDIso = eleMcMatch.clone(src = cms.InputTag("probeTauSCIDIso"))
eleMcMatchProbeTauSCIDIsoNoCracks = eleMcMatch.clone(src = cms.InputTag("probeTauSCIDIsoNoCracks"))

mcMatches90 = cms.Sequence(eleMcMatchTag90+eleMcMatchProbeTauIDIsoLoose+eleMcMatchProbeTauIDIsoMedium+eleMcMatchProbeTauIDIsoTight+eleMcMatchProbeTauIDIsoMediumNoCracks+eleMcMatchProbeTauIDIsoLooseNoCracks+eleMcMatchProbeTauIDIsoTightNoCracks)

mcMatchesSC90 = cms.Sequence(eleMcMatchTag90+eleMcMatchProbeTauSCIDIso+eleMcMatchProbeTauSCIDIsoNoCracks)

mcMatches80 = cms.Sequence(eleMcMatchTag80+eleMcMatchProbeTauIDIsoLoose+eleMcMatchProbeTauIDIsoMedium+eleMcMatchProbeTauIDIsoTight+eleMcMatchProbeTauIDIsoMediumNoCracks+eleMcMatchProbeTauIDIsoLooseNoCracks+eleMcMatchProbeTauIDIsoTightNoCracks)

mcMatchesSC80 = cms.Sequence(eleMcMatchTag80+eleMcMatchProbeTauSCIDIso+eleMcMatchProbeTauSCIDIsoNoCracks)

mcMatches70 = cms.Sequence(eleMcMatchTag70+eleMcMatchProbeTauIDIsoLoose+eleMcMatchProbeTauIDIsoMedium+eleMcMatchProbeTauIDIsoTight+eleMcMatchProbeTauIDIsoMediumNoCracks+eleMcMatchProbeTauIDIsoLooseNoCracks+eleMcMatchProbeTauIDIsoTightNoCracks)

mcMatchesSC70 = cms.Sequence(eleMcMatchTag70+eleMcMatchProbeTauSCIDIso+eleMcMatchProbeTauSCIDIsoNoCracks)

mcMatches60 = cms.Sequence(eleMcMatchTag60+eleMcMatchProbeTauIDIsoLoose+eleMcMatchProbeTauIDIsoMedium+eleMcMatchProbeTauIDIsoTight+eleMcMatchProbeTauIDIsoMediumNoCracks+eleMcMatchProbeTauIDIsoLooseNoCracks+eleMcMatchProbeTauIDIsoTightNoCracks)

mcMatchesSC60 = cms.Sequence(eleMcMatchTag60+eleMcMatchProbeTauSCIDIso+eleMcMatchProbeTauSCIDIsoNoCracks)

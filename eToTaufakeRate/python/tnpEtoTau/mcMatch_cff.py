import FWCore.ParameterSet.Config as cms

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

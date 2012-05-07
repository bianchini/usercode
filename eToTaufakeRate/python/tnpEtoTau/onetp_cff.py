import FWCore.ParameterSet.Config as cms

oneTp = cms.EDFilter("CandViewCountFilter",
                     	src = cms.InputTag("tnpAnyEleAnyTau"),
                     	minNumber = cms.uint32(1)
                     	)

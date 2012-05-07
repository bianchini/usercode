import FWCore.ParameterSet.Config as cms

tnpAnyEleAnyTau = cms.EDProducer("CandViewShallowCloneCombiner",
                          decay = cms.string("tagAnyEle@+ probeAnyTau@-"),
			  roles = cms.vstring('ele', 'tau'),
                          cut   = cms.string("30 < mass < 120"),
                          #checkCharge = cms.bool(False),   
                          )

import FWCore.ParameterSet.Config as cms

WenuMtCut = cms.EDFilter("WenuFilter",
  eleLabel =  cms.InputTag("selectedPatElectronsPFlowIso"),
  metLabel =  cms.InputTag('patMETsPFlow'),
  minMt = cms.double(60),
)

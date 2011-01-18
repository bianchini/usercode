import FWCore.ParameterSet.Config as cms


probeAnyTauc = cms.EDFilter("PATTauRefSelector",
                            src = cms.InputTag("patTausPFlow"),
                            cut = cms.string('pt>15.0 && abs(eta)<2.5'),
                            filter = cms.bool(False)
                            )

## HPS
probeTauID = probeAnyTauc.clone(
    cut = cms.string(probeAnyTauc.cut.value()+' && tauID("leadingTrackFinding") > 0.5')
    )
probeTauNoCracks = probeAnyTauc.clone(
    cut = cms.string(probeAnyTauc.cut.value()+' && tauID("againstElectronCrackRem") > 0.5')
    )
probeTauLooseIso = probeAnyTauc.clone(
    cut = cms.string(probeAnyTauc.cut.value()+' && tauID("byLooseIsolation") > 0.5')
    )
probeTauMediumIso = probeAnyTauc.clone(
    cut = cms.string(probeAnyTauc.cut.value()+' && tauID("byMediumIsolation") > 0.5')
    )
probeTauTightIso = probeAnyTauc.clone(
    cut = cms.string(probeAnyTauc.cut.value()+' && tauID("byTightIsolation") > 0.5')
    )
probeTauIDIsoLoose = probeAnyTauc.clone(
    cut = cms.string(probeTauID.cut.value()+' && '+probeTauLooseIso.cut.value())
    )
probeTauIDIsoMedium = probeAnyTauc.clone(
    cut = cms.string(probeTauID.cut.value()+' && '+probeTauMediumIso.cut.value())
    )
probeTauIDIsoTight = probeAnyTauc.clone(
    cut = cms.string(probeTauID.cut.value()+' && '+probeTauTightIso.cut.value())
    )

probeTauIDIsoLooseNoCracks = probeAnyTauc.clone(
    cut = cms.string(probeTauID.cut.value()+' && '+probeTauLooseIso.cut.value()+' && '+probeTauNoCracks.cut.value())
    )
probeTauIDIsoMediumNoCracks = probeAnyTauc.clone(
    cut = cms.string(probeTauID.cut.value()+' && '+probeTauMediumIso.cut.value()+' && '+probeTauNoCracks.cut.value())
    )
probeTauIDIsoTightNoCracks = probeAnyTauc.clone(
    cut = cms.string(probeTauID.cut.value()+' && '+probeTauTightIso.cut.value()+' && '+probeTauNoCracks.cut.value())
    )

## Shrinking Cone
probeAnyTaucSC = probeAnyTauc.clone(
    src = cms.InputTag("patTaus"),
    )
probeTauSCID = probeAnyTaucSC.clone(    
    cut = cms.string(probeAnyTauc.cut.value()+' && tauID("leadingTrackFinding") > 0.5 && tauID("leadingPionPtCut") > 0.5')
    )
probeTauSCIso = probeAnyTaucSC.clone(    
    cut = cms.string(probeAnyTauc.cut.value()+' && tauID("byIsolation") > 0.5')
    )
probeTauSCIDIso = probeAnyTaucSC.clone(    
    cut = cms.string(probeAnyTauc.cut.value()+' && '+probeTauSCID.cut.value()+ ' && ' +probeTauSCIso.cut.value())
    )
probeTauSCIDIsoNoCracks = probeAnyTaucSC.clone(    
    cut = cms.string(probeAnyTauc.cut.value()+' && '+probeTauSCID.cut.value()+ ' && ' +probeTauSCIso.cut.value()+' && '+probeTauNoCracks.cut.value())
    )

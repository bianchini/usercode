import FWCore.ParameterSet.Config as cms

probeAnyTau = cms.EDFilter("PATTauRefSelector",
                         src = cms.InputTag("patTaus"),
                         cut = cms.string('pfJetRef.pt > 20.0 && abs(eta) < 2.3'),
                         filter = cms.bool(False)
                         )

passingIsoLooseEleVetoLoose = probeAnyTau.clone(
                        cut = cms.string(probeAnyTau.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstElectronLoose") > 0.5')
                        )

passingIsoLooseEleVetoMedium = probeAnyTau.clone(
                        cut = cms.string(probeAnyTau.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstElectronMedium") > 0.5')
                        )

passingIsoLooseEleVetoTight = probeAnyTau.clone(
                        cut = cms.string(probeAnyTau.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstElectronTight") > 0.5')
                        )

passingIsoLooseEleVetoMVA = probeAnyTau.clone(
                        cut = cms.string(probeAnyTau.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstElectronMedium") > 0.5 && tauID("againstElectronMVA") > 0.5')
                        )

passingProbes = cms.Sequence(probeAnyTau* 
			    (passingIsoLooseEleVetoLoose +
			     passingIsoLooseEleVetoMedium +
			     passingIsoLooseEleVetoTight +
			     passingIsoLooseEleVetoMVA
			    )
			   )

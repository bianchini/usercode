import FWCore.ParameterSet.Config as cms

etoTau = cms.EDAnalyzer("TagProbeFitTreeProducer",
                        tagProbePairs = cms.InputTag("tnpAnyEleAnyTau"),
                        arbitration   = cms.string("OneProbe"),
                        #massForArbitration = cms.double(91),
                        variables = cms.PSet(
				pt = cms.string("pt"),
				abseta = cms.string("abs(eta)"),
				phi = cms.string("phi"),
				decayModeFinding = cms.string("tauID('decayModeFinding')"),
				hpsLooseIso = cms.string("tauID('byLooseCombinedIsolationDeltaBetaCorr')"),
				hpsMediumIso = cms.string("tauID('byMediumCombinedIsolationDeltaBetaCorr')"),
				hpsTightIso = cms.string("tauID('byTightCombinedIsolationDeltaBetaCorr')"),
				looseEleDis = cms.string("tauID('againstElectronLoose')"),
				mediumEleDis = cms.string("tauID('againstElectronMedium')"),
				tightEleDis = cms.string("tauID('againstElectronTight')"),
				mvaEleDis = cms.string("tauID('againstElectronMVA')"),
				tauTriggerMatchingTauXTrigger = cms.InputTag("addUserVariables","tauXTriggersLeg2"),
    			),
                        flags = cms.PSet(
				passingIsoLooseEleVetoLoose =  cms.InputTag("passingIsoLooseEleVetoLoose"),
				passingIsoLooseEleVetoMedium =  cms.InputTag("passingIsoLooseEleVetoMedium"),
				passingIsoLooseEleVetoTight =  cms.InputTag("passingIsoLooseEleVetoTight"),
				passingIsoLooseEleVetoMVA =  cms.InputTag("passingIsoLooseEleVetoMVA"),
    			),
			tagVariables = cms.PSet(
				Mt =  cms.InputTag("addUserVariables","Mt"),
    				puMCWeight2011A =  cms.InputTag("addUserVariables","puMCWeightRun2011A"),
    				puMCWeight2011B =  cms.InputTag("addUserVariables","puMCWeightRun2011B"),
				eleTriggerBitTauXTrigger = cms.InputTag("addUserVariables","triggerBit"),
				eleTriggerMatchingTauXTrigger = cms.InputTag("addUserVariables","tauXTriggersLeg1"),
				eleTriggerBitEleXTrigger = cms.InputTag("addUserVariables","triggerBitSingleEle"),
				eleTriggerMatchingEleXTrigger = cms.InputTag("addUserVariables","eleXTriggersLeg1"),
				eleIdWP70 = cms.InputTag("addUserVariables","eleIdWP70"),
				numEleWP95 = cms.InputTag("addUserVariables","numEleWP95"),
				numEleWP80 = cms.InputTag("addUserVariables","numEleWP80"),
			),
			tagFlags = cms.PSet(),
			pairVariables = cms.PSet(
				dr = cms.string("deltaR(daughter('ele').eta,daughter('ele').phi,daughter('tau').eta,daughter('tau').phi)"),
				charge = cms.string("charge")
			),
			pairFlags = cms.PSet(
			),
                        isMC = cms.bool(False), ######
                        tagMatches = cms.InputTag("tagMcMatch") ,
                        probeMatches  = cms.InputTag("probeMcMatch"),
                        motherPdgId = cms.vint32(23),
                        makeMCUnbiasTree = cms.bool(False), ######
                        checkMotherInUnbiasEff = cms.bool(False), ######
                        allProbes = cms.InputTag("probeAnyTau"),
                        addRunLumiInfo = cms.bool(True),
			addEventVariablesInfo = cms.bool(True)
                        )

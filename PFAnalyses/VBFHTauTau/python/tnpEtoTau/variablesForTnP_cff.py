import FWCore.ParameterSet.Config as cms


variables = cms.PSet(
    pt = cms.string("pt"), 
    abseta = cms.string("abs(eta)"),
    phi = cms.string("phi"),
    emFraction = cms.string("emFraction"),
    hcalTotOverPLead = cms.string("hcalTotOverPLead"),
    hcalMaxOverPLead = cms.string("hcalMaxOverPLead"),
    hcalEnergy =  cms.string("leadPFChargedHadrCand.hcalEnergy"),
    ecalEnergy =  cms.string("leadPFChargedHadrCand.ecalEnergy"),
    hcalEnergyLeadPFCand =  cms.string("leadPFCand.hcalEnergy"),
    ecalEnergyLeadPFCand =  cms.string("leadPFCand.ecalEnergy"),
    bremsRecoveryEOverPLead = cms.string("bremsRecoveryEOverPLead"),
    electronPreIDOutput = cms.string("electronPreIDOutput"),
    electronPreIDOutputLeadPFCand = cms.string("leadPFCand.mva_e_pi"),
    hcal3x3OverPLead = cms.string("hcal3x3OverPLead"),
    ecalStripSumEOverPLead = cms.string("ecalStripSumEOverPLead"),
    leadPFChargedHadrCandTrackPt = cms.string("leadPFChargedHadrCand.trackRef.pt"),
    leadPFChargedHadrCandTrackP = cms.string("leadPFChargedHadrCand.trackRef.p"),
    leadPFCandPt = cms.string("leadPFCand.pt"),
    leadPFCandP = cms.string("leadPFCand.p"),
    signalPFChargedHadrCands = cms.string("signalPFChargedHadrCands.size"),
    signalPFGammaCands = cms.string("signalPFGammaCands.size"),
    nHitsLeadTrack = cms.string("leadPFChargedHadrCand.trackRef.numberOfValidHits"),
    matchedID = cms.InputTag("tightestMatchedID"),
    )

variablesSC = cms.PSet(
    pt = cms.string("pt"), 
    abseta = cms.string("abs(eta)"),
    phi = cms.string("phi"),
    emFraction = cms.string("emFraction"),
    hcalTotOverPLead = cms.string("hcalTotOverPLead"),
    hcalMaxOverPLead = cms.string("hcalMaxOverPLead"),
    hcalEnergy =  cms.string("leadPFChargedHadrCand.hcalEnergy"),
    ecalEnergy =  cms.string("leadPFChargedHadrCand.ecalEnergy"),
    hcalEnergyLeadPFCand =  cms.string("leadPFCand.hcalEnergy"),
    ecalEnergyLeadPFCand =  cms.string("leadPFCand.ecalEnergy"),
    bremsRecoveryEOverPLead = cms.string("bremsRecoveryEOverPLead"),
    electronPreIDOutput = cms.string("electronPreIDOutput"),
    electronPreIDOutputLeadPFCand = cms.string("leadPFCand.mva_e_pi"),
    hcal3x3OverPLead = cms.string("hcal3x3OverPLead"),
    ecalStripSumEOverPLead = cms.string("ecalStripSumEOverPLead"),
    leadPFChargedHadrCandTrackPt = cms.string("leadPFChargedHadrCand.trackRef.pt"),
    leadPFChargedHadrCandTrackP = cms.string("leadPFChargedHadrCand.trackRef.p"),
    leadPFCandPt = cms.string("leadPFCand.pt"),
    leadPFCandP = cms.string("leadPFCand.p"),
    signalPFChargedHadrCands = cms.string("signalPFChargedHadrCands.size"),
    signalPFGammaCands = cms.string("signalPFGammaCands.size"),
    nHitsLeadTrack = cms.string("leadPFChargedHadrCand.trackRef.numberOfValidHits"),
    matchedID = cms.InputTag("tightestMatchedIDSC"),
    )

flags = cms.PSet(
    tauID = cms.string('tauID("leadingTrackFinding")>0.5'),
    tauLooseIso =  cms.string('tauID("byLooseIsolation")>0.5'),
    tauMediumIso =  cms.string('tauID("byMediumIsolation")>0.5'),
    tauTightIso =  cms.string('tauID("byTightIsolation")>0.5'),
    tauAntiEMVA =  cms.string('tauID("againstElectron")>0.5'),
    tauAntiECrackRem =  cms.string('tauID("againstElectronCrackRem")>0.5'),
    #tauAntiE2D =  cms.string('tauID("againstElectron2D")>0.5'),
    #tauAntiEMVAandCrackRem =  cms.string('tauID("againstElectron")>0.5 && tauID("againstElectronCrackRem")>0.5'),
    #tauAntiE2DandCrackRem =  cms.string('tauID("againstElectron2D")>0.5 && tauID("againstElectronCrackRem")>0.5')
    )

flagsForMarg = cms.PSet(
    tauAntiEMVA      =  cms.string('tauID("againstElectron")>0.5'),
    tauAntiECrackRem =  cms.string('tauID("againstElectronCrackRem")>0.5'),
    tauAntiEMVA_HoP  =  cms.string('tauID("againstElectron")>0.5 && leadPFChargedHadrCand.hcalEnergy/leadPFChargedHadrCand.trackRef.p>0.1'),
    #tauAntiE2D =  cms.string('tauID("againstElectron2D")>0.5'),
    #tauAntiEMVAandCrackRem =  cms.string('tauID("againstElectron")>0.5 && tauID("againstElectronCrackRem")>0.5'),
    #tauAntiE2DandCrackRem =  cms.string('tauID("againstElectron2D")>0.5 && tauID("againstElectronCrackRem")>0.5'),
    #tauAntiE2D_2 = cms.string('(electronPreIDDecision>0.5 && ((leadPFChargedHadrCand.hcalEnergy)/leadPFChargedHadrCand.trackRef.p >0.15 || (leadPFChargedHadrCand.ecalEnergy)/leadPFChargedHadrCand.trackRef.p<0.8)) || (electronPreIDDecision<0.5 && ((leadPFChargedHadrCand.hcalEnergy)/leadPFChargedHadrCand.trackRef.p>0.05 || (leadPFChargedHadrCand.ecalEnergy)/leadPFChargedHadrCand.trackRef.p<0.95))')
    )

flagsForMargNoCracks = cms.PSet(
    tauAntiEMVA =  cms.string('tauID("againstElectron")>0.5'),
    tauAntiEMVA_HoP  =  cms.string('tauID("againstElectron")>0.5 && leadPFChargedHadrCand.hcalEnergy/leadPFChargedHadrCand.trackRef.p>0.1'),
    #tauAntiE2D =  cms.string('tauID("againstElectron2D")>0.5'),
    #tauAntiE2D_2 = cms.string('(electronPreIDDecision>0.5 && ((leadPFChargedHadrCand.hcalEnergy)/leadPFChargedHadrCand.trackRef.p >0.15 || (leadPFChargedHadrCand.ecalEnergy)/leadPFChargedHadrCand.trackRef.p<0.8)) || (electronPreIDDecision<0.5 && ((leadPFChargedHadrCand.hcalEnergy)/leadPFChargedHadrCand.trackRef.p>0.05 || (leadPFChargedHadrCand.ecalEnergy)/leadPFChargedHadrCand.trackRef.p<0.95))')
    )

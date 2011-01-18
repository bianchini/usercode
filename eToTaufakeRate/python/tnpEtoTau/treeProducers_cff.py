import FWCore.ParameterSet.Config as cms

from Bianchi.eToTaufakeRate.tnpEtoTau.variablesForTnP_cff import *


## HPS

etoTau = cms.EDAnalyzer("TagProbeFitTreeProducer",
                        tagProbePairs = cms.InputTag("tnpAnyEleAnyTau"),
                        arbitration   = cms.string("OneProbe"),
                        #massForArbitration = cms.double(91.2),
                        variables = cms.PSet(
    pt = cms.string("pt"),
    ),
                        flags = cms.PSet(
    dummy =  cms.string("pt>0"),
    ),
                        isMC = cms.bool( True ),
                        tagMatches = cms.InputTag("eleMcMatch") ,
                        probeMatches  = cms.InputTag("probeMatch"),
                        motherPdgId = cms.vint32(11,-11),
                        makeMCUnbiasTree = cms.bool(True),
                        checkMotherInUnbiasEff = cms.bool(True),
                        allProbes = cms.InputTag("probeAnyTau"),
                        addRunLumiInfo = cms.bool(True)
                        )

etoTauMargLoose90 = etoTau.clone(
    tagProbePairs = cms.InputTag("tnp90TauIDLoose"),
    variables = variables,
    flags = flagsForMarg,
    tagMatches = cms.InputTag("eleMcMatchTag90"),
    probeMatches  = cms.InputTag("eleMcMatchProbeTauIDIsoLoose"),
    allProbes = cms.InputTag("probeTauIDIsoLoose")
    )
etoTauMargLoose80 = etoTauMargLoose90.clone(
    tagProbePairs =  cms.InputTag("tnp80TauIDLoose"),
    tagMatches = cms.InputTag("eleMcMatchTag80"),
    )
etoTauMargLoose70 = etoTauMargLoose90.clone(
    tagProbePairs =  cms.InputTag("tnp70TauIDLoose"),
    tagMatches = cms.InputTag("eleMcMatchTag70"),
    )
etoTauMargLoose60 = etoTauMargLoose90.clone(
    tagProbePairs =  cms.InputTag("tnp60TauIDLoose"),
    tagMatches = cms.InputTag("eleMcMatchTag60"),
    )

etoTauMargMedium90 = etoTau.clone(
    tagProbePairs = cms.InputTag("tnp90TauIDMedium"),
    variables = variables,
    flags = flagsForMarg,
    tagMatches = cms.InputTag("eleMcMatchTag90"),
    probeMatches  = cms.InputTag("eleMcMatchProbeTauIDIsoMedium"),
    allProbes = cms.InputTag("probeTauIDIsoMedium")
    )
etoTauMargMedium80 = etoTauMargMedium90.clone(
    tagProbePairs =  cms.InputTag("tnp80TauIDMedium"),
    tagMatches = cms.InputTag("eleMcMatchTag80"),
    )
etoTauMargMedium70 = etoTauMargMedium90.clone(
    tagProbePairs =  cms.InputTag("tnp70TauIDMedium"),
    tagMatches = cms.InputTag("eleMcMatchTag70"),
    )
etoTauMargMedium60 = etoTauMargMedium90.clone(
    tagProbePairs =  cms.InputTag("tnp60TauIDMedium"),
    tagMatches = cms.InputTag("eleMcMatchTag60"),
    )

etoTauMargTight90 = etoTau.clone(
    tagProbePairs = cms.InputTag("tnp90TauIDTight"),
    variables = variables,
    flags = flagsForMarg,
    tagMatches = cms.InputTag("eleMcMatchTag90"),
    probeMatches  = cms.InputTag("eleMcMatchProbeTauIDIsoTight"),
    allProbes = cms.InputTag("probeTauIDIsoTight")
    )
etoTauMargTight80 = etoTauMargTight90.clone(
    tagProbePairs =  cms.InputTag("tnp80TauIDTight"),
    tagMatches = cms.InputTag("eleMcMatchTag80"),
    )
etoTauMargTight70 = etoTauMargTight90.clone(
    tagProbePairs =  cms.InputTag("tnp70TauIDTight"),
    tagMatches = cms.InputTag("eleMcMatchTag70"),
    )
etoTauMargTight60 = etoTauMargTight90.clone(
    tagProbePairs =  cms.InputTag("tnp60TauIDTight"),
    tagMatches = cms.InputTag("eleMcMatchTag60"),
    )



etoTauMargLooseNoCracks90 = etoTau.clone(
    tagProbePairs = cms.InputTag("tnp90TauIDLooseNoCracks"),
    variables = variables,
    flags = flagsForMargNoCracks,
    tagMatches = cms.InputTag("eleMcMatchTag90"),
    probeMatches  = cms.InputTag("eleMcMatchProbeTauIDIsoLooseNoCracks"),
    allProbes = cms.InputTag("probeTauIDIsoLooseNoCracks")
    )
etoTauMargLooseNoCracks80 = etoTauMargLooseNoCracks90.clone(
    tagProbePairs =  cms.InputTag("tnp80TauIDLooseNoCracks"),
    tagMatches = cms.InputTag("eleMcMatchTag80"),
    )
etoTauMargLooseNoCracks70 = etoTauMargLooseNoCracks90.clone(
    tagProbePairs =  cms.InputTag("tnp70TauIDLooseNoCracks"),
    tagMatches = cms.InputTag("eleMcMatchTag70"),
    )
etoTauMargLooseNoCracks60 = etoTauMargLooseNoCracks90.clone(
    tagProbePairs =  cms.InputTag("tnp60TauIDLooseNoCracks"),
    tagMatches = cms.InputTag("eleMcMatchTag60"),
    )

etoTauMargMediumNoCracks90 = etoTau.clone(
    tagProbePairs = cms.InputTag("tnp90TauIDMediumNoCracks"),
    variables = variables,
    flags = flagsForMargNoCracks,
    tagMatches = cms.InputTag("eleMcMatchTag90"),
    probeMatches  = cms.InputTag("eleMcMatchProbeTauIDIsoMediumNoCracks"),
    allProbes = cms.InputTag("probeTauIDIsoMediumNoCracks")
    )
etoTauMargMediumNoCracks80 = etoTauMargMediumNoCracks90.clone(
    tagProbePairs =  cms.InputTag("tnp80TauIDMediumNoCracks"),
    tagMatches = cms.InputTag("eleMcMatchTag80"),
    )
etoTauMargMediumNoCracks70 = etoTauMargMediumNoCracks90.clone(
    tagProbePairs =  cms.InputTag("tnp70TauIDMediumNoCracks"),
    tagMatches = cms.InputTag("eleMcMatchTag70"),
    )
etoTauMargMediumNoCracks60 = etoTauMargMediumNoCracks90.clone(
    tagProbePairs =  cms.InputTag("tnp60TauIDMediumNoCracks"),
    tagMatches = cms.InputTag("eleMcMatchTag60"),
    )

etoTauMargTightNoCracks90 = etoTau.clone(
    tagProbePairs = cms.InputTag("tnp90TauIDTightNoCracks"),
    variables = variables,
    flags = flagsForMargNoCracks,
    tagMatches = cms.InputTag("eleMcMatchTag90"),
    probeMatches  = cms.InputTag("eleMcMatchProbeTauIDIsoTightNoCracks"),
    allProbes = cms.InputTag("probeTauIDIsoTightNoCracks")
    )
etoTauMargTightNoCracks80 = etoTauMargTightNoCracks90.clone(
    tagProbePairs =  cms.InputTag("tnp80TauIDTightNoCracks"),
    tagMatches = cms.InputTag("eleMcMatchTag80"),
    )
etoTauMargTightNoCracks70 = etoTauMargTightNoCracks90.clone(
    tagProbePairs =  cms.InputTag("tnp70TauIDTightNoCracks"),
    tagMatches = cms.InputTag("eleMcMatchTag70"),
    )
etoTauMargTightNoCracks60 = etoTauMargTightNoCracks90.clone(
    tagProbePairs =  cms.InputTag("tnp60TauIDTightNoCracks"),
    tagMatches = cms.InputTag("eleMcMatchTag60"),
    )



## Shrinking Cone

etoTauSCMarg90 = etoTau.clone(
    tagProbePairs = cms.InputTag("tnp90TauSCIDIso"),
    variables = variablesSC,
    flags = flagsForMarg,
    tagMatches = cms.InputTag("eleMcMatchTag90"),
    probeMatches  = cms.InputTag("eleMcMatchProbeTauSCIDIso"),
    allProbes = cms.InputTag("probeTauSCIDIso")
    )
etoTauSCMarg80 = etoTauSCMarg90.clone(
    tagProbePairs = cms.InputTag("tnp80TauSCIDIso"),
    tagMatches = cms.InputTag("eleMcMatchTag80"),
    )
etoTauSCMarg70 = etoTauSCMarg90.clone(
    tagProbePairs = cms.InputTag("tnp70TauSCIDIso"),
    tagMatches = cms.InputTag("eleMcMatchTag70"),
    )
etoTauSCMarg60 = etoTauSCMarg90.clone(
    tagProbePairs = cms.InputTag("tnp60TauSCIDIso"),
    tagMatches = cms.InputTag("eleMcMatchTag60"),
    )

etoTauSCMargNoCracks90 = etoTau.clone(
    tagProbePairs = cms.InputTag("tnp90TauSCIDIsoNoCracks"),
    variables = variablesSC,
    flags = flagsForMargNoCracks,
    tagMatches = cms.InputTag("eleMcMatchTag90"),
    probeMatches  = cms.InputTag("eleMcMatchProbeTauSCIDIsoNoCracks"),
    allProbes = cms.InputTag("probeTauSCIDIsoNoCracks")
    )
etoTauSCMargNoCracks80 = etoTauSCMargNoCracks90.clone(
    tagProbePairs = cms.InputTag("tnp80TauSCIDIsoNoCracks"),
    tagMatches = cms.InputTag("eleMcMatchTag80"),
    )
etoTauSCMargNoCracks70 = etoTauSCMargNoCracks90.clone(
    tagProbePairs = cms.InputTag("tnp70TauSCIDIsoNoCracks"),
    tagMatches = cms.InputTag("eleMcMatchTag70"),
    )
etoTauSCMargNoCracks60 = etoTauSCMargNoCracks90.clone(
    tagProbePairs = cms.InputTag("tnp60TauSCIDIsoNoCracks"),
    tagMatches = cms.InputTag("eleMcMatchTag60"),
    )

import FWCore.ParameterSet.Config as cms


## HPS

tnpAnyEleAnyTau = cms.EDProducer("CandViewShallowCloneCombiner",
                                 decay = cms.string("tagAnyEle@+ probeAnyTau@-"),
                                 cut   = cms.string("30 < mass < 120"),
                                 #checkCharge = cms.bool(False),   
                                 )


tnp90AnyTau = tnpAnyEleAnyTau.clone(decay = cms.string("tag90@+ probeAnyTau@-"))
tnp80AnyTau = tnpAnyEleAnyTau.clone(decay = cms.string("tag80@+ probeAnyTau@-"))
tnp70AnyTau = tnpAnyEleAnyTau.clone(decay = cms.string("tag70@+ probeAnyTau@-"))
tnp60AnyTau = tnpAnyEleAnyTau.clone(decay = cms.string("tag60@+ probeAnyTau@-"))

tnp90TauIDLoose = tnpAnyEleAnyTau.clone(decay = cms.string("tag90@+ probeTauIDIsoLoose@-"))
tnp90TauIDMedium = tnpAnyEleAnyTau.clone(decay = cms.string("tag90@+ probeTauIDIsoMedium@-"))
tnp90TauIDTight = tnpAnyEleAnyTau.clone(decay = cms.string("tag90@+ probeTauIDIsoTight@-"))

tnp80TauIDLoose = tnpAnyEleAnyTau.clone(decay = cms.string("tag80@+ probeTauIDIsoLoose@-"))
tnp80TauIDMedium = tnpAnyEleAnyTau.clone(decay = cms.string("tag80@+ probeTauIDIsoMedium@-"))
tnp80TauIDTight = tnpAnyEleAnyTau.clone(decay = cms.string("tag80@+ probeTauIDIsoTight@-"))

tnp70TauIDLoose = tnpAnyEleAnyTau.clone(decay = cms.string("tag70@+ probeTauIDIsoLoose@-"))
tnp70TauIDMedium = tnpAnyEleAnyTau.clone(decay = cms.string("tag70@+ probeTauIDIsoMedium@-"))
tnp70TauIDTight = tnpAnyEleAnyTau.clone(decay = cms.string("tag70@+ probeTauIDIsoTight@-"))

tnp60TauIDLoose = tnpAnyEleAnyTau.clone(decay = cms.string("tag60@+ probeTauIDIsoLoose@-"))
tnp60TauIDMedium = tnpAnyEleAnyTau.clone(decay = cms.string("tag60@+ probeTauIDIsoMedium@-"))
tnp60TauIDTight = tnpAnyEleAnyTau.clone(decay = cms.string("tag60@+ probeTauIDIsoTight@-"))

tnp90TauIDLooseNoCracks = tnpAnyEleAnyTau.clone(decay = cms.string("tag90@+ probeTauIDIsoLooseNoCracks@-"))
tnp90TauIDMediumNoCracks = tnpAnyEleAnyTau.clone(decay = cms.string("tag90@+ probeTauIDIsoMediumNoCracks@-"))
tnp90TauIDTightNoCracks = tnpAnyEleAnyTau.clone(decay = cms.string("tag90@+ probeTauIDIsoTightNoCracks@-"))

tnp80TauIDLooseNoCracks = tnpAnyEleAnyTau.clone(decay = cms.string("tag80@+ probeTauIDIsoLooseNoCracks@-"))
tnp80TauIDMediumNoCracks = tnpAnyEleAnyTau.clone(decay = cms.string("tag80@+ probeTauIDIsoMediumNoCracks@-"))
tnp80TauIDTightNoCracks = tnpAnyEleAnyTau.clone(decay = cms.string("tag80@+ probeTauIDIsoTightNoCracks@-"))

tnp70TauIDLooseNoCracks = tnpAnyEleAnyTau.clone(decay = cms.string("tag70@+ probeTauIDIsoLooseNoCracks@-"))
tnp70TauIDMediumNoCracks = tnpAnyEleAnyTau.clone(decay = cms.string("tag70@+ probeTauIDIsoMediumNoCracks@-"))
tnp70TauIDTightNoCracks = tnpAnyEleAnyTau.clone(decay = cms.string("tag70@+ probeTauIDIsoTightNoCracks@-"))

tnp60TauIDLooseNoCracks = tnpAnyEleAnyTau.clone(decay = cms.string("tag60@+ probeTauIDIsoLooseNoCracks@-"))
tnp60TauIDMediumNoCracks = tnpAnyEleAnyTau.clone(decay = cms.string("tag60@+ probeTauIDIsoMediumNoCracks@-"))
tnp60TauIDTightNoCracks = tnpAnyEleAnyTau.clone(decay = cms.string("tag60@+ probeTauIDIsoTightNoCracks@-"))


## Shrinking Cone

tnpAnyEleAnyTauSC = tnpAnyEleAnyTau.clone(
    decay = cms.string("tagAnyEle@+ probeAnyTauSC@-")
    )

tnp90AnyTauSC = tnpAnyEleAnyTauSC.clone(
    decay = cms.string("tag90@+ probeAnyTauSC@-")
    )
tnp80AnyTauSC = tnpAnyEleAnyTauSC.clone(
    decay = cms.string("tag80@+ probeAnyTauSC@-")
    )
tnp70AnyTauSC = tnpAnyEleAnyTauSC.clone(
    decay = cms.string("tag70@+ probeAnyTauSC@-")
    )
tnp60AnyTauSC = tnpAnyEleAnyTauSC.clone(
    decay = cms.string("tag60@+ probeAnyTauSC@-")
    )

tnp90TauSCIDIso = tnpAnyEleAnyTauSC.clone(
    decay = cms.string("tag90@+ probeTauSCIDIso@-")
    )
tnp80TauSCIDIso = tnpAnyEleAnyTauSC.clone(
    decay = cms.string("tag80@+ probeTauSCIDIso@-")
    )
tnp70TauSCIDIso = tnpAnyEleAnyTauSC.clone(
    decay = cms.string("tag70@+ probeTauSCIDIso@-")
    )
tnp60TauSCIDIso = tnpAnyEleAnyTauSC.clone(
    decay = cms.string("tag60@+ probeTauSCIDIso@-")
    )

tnp90TauSCIDIsoNoCracks = tnpAnyEleAnyTauSC.clone(
    decay = cms.string("tag90@+ probeTauSCIDIsoNoCracks@-")
    )
tnp80TauSCIDIsoNoCracks = tnpAnyEleAnyTauSC.clone(
    decay = cms.string("tag80@+ probeTauSCIDIsoNoCracks@-")
    )
tnp70TauSCIDIsoNoCracks = tnpAnyEleAnyTauSC.clone(
    decay = cms.string("tag70@+ probeTauSCIDIsoNoCracks@-")
    )
tnp60TauSCIDIsoNoCracks = tnpAnyEleAnyTauSC.clone(
    decay = cms.string("tag60@+ probeTauSCIDIsoNoCracks@-")
    )

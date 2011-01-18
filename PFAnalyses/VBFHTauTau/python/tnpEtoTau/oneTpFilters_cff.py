import FWCore.ParameterSet.Config as cms


oneTp = cms.EDFilter("CandViewCountFilter",
                     src = cms.InputTag("tnpAnyEleAnyTau"),
                     minNumber = cms.uint32(1)
                     )

oneTp90Loose = oneTp.clone( src = cms.InputTag("tnp90TauIDLoose") )
oneTp80Loose = oneTp.clone( src = cms.InputTag("tnp80TauIDLoose") )
oneTp70Loose = oneTp.clone( src = cms.InputTag("tnp70TauIDLoose") )
oneTp60Loose = oneTp.clone( src = cms.InputTag("tnp60TauIDLoose") )

oneTpSC90 = oneTp.clone( src = cms.InputTag("tnp90TauSCIDIso") )
oneTpSC80 = oneTp.clone( src = cms.InputTag("tnp80TauSCIDIso") )
oneTpSC70 = oneTp.clone( src = cms.InputTag("tnp70TauSCIDIso") )
oneTpSC60 = oneTp.clone( src = cms.InputTag("tnp60TauSCIDIso") )

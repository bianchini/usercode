import FWCore.ParameterSet.Config as cms


WenuPair = cms.EDProducer("CandViewShallowCloneCombiner",
                          decay = cms.string("tagAnyEle patMETsPFlow"),
                          cut   = cms.string('sqrt((daughter(0).pt+daughter(1).pt)*(daughter(0).pt+daughter(1).pt)-pt*pt)>50'),
                          checkCharge = cms.bool(False)
                          )

WenuPair90 = WenuPair.clone( decay =  cms.string("tag90 patMETsPFlow") )
noWenuFilter90 = cms.EDFilter("CandViewCountFilter",
                              src = cms.InputTag("WenuPair90"),
                              minNumber = cms.uint32(1),
                              )


WenuPair80 = WenuPair90.clone( decay = cms.string("tag80 patMETsPFlow") )
noWenuFilter80 = noWenuFilter90.clone( src =  cms.InputTag("WenuPair80"))

WenuPair70 = WenuPair90.clone( decay = cms.string("tag70 patMETsPFlow") )
noWenuFilter70 = noWenuFilter90.clone( src =  cms.InputTag("WenuPair70"))

WenuPair60 = WenuPair90.clone( decay = cms.string("tag60 patMETsPFlow") )
noWenuFilter60 = noWenuFilter90.clone( src =  cms.InputTag("WenuPair60"))

## sequences

noWenu90 = cms.Sequence( WenuPair90*(~noWenuFilter90))
noWenu80 = cms.Sequence( WenuPair80*(~noWenuFilter80))
noWenu70 = cms.Sequence( WenuPair70*(~noWenuFilter70))
noWenu60 = cms.Sequence( WenuPair60*(~noWenuFilter60))

noWenuSC90 = cms.Sequence( WenuPair90*(~noWenuFilter90))
noWenuSC80 = cms.Sequence( WenuPair80*(~noWenuFilter80))
noWenuSC70 = cms.Sequence( WenuPair70*(~noWenuFilter70))
noWenuSC60 = cms.Sequence( WenuPair60*(~noWenuFilter60))

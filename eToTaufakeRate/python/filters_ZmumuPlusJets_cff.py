import FWCore.ParameterSet.Config as cms


hltHighLevelMuons = cms.EDFilter("HLTHighLevel",
                                 TriggerResultsTag = cms.InputTag("TriggerResults","","REDIGI38X"),
                                 HLTPaths = cms.vstring('HLT_Mu9','HLT_Mu11','HLT_IsoMu9'),
                                 eventSetupPathsKey = cms.string(''),                 
                                 andOr = cms.bool(True), #True (OR) accept if ANY is true, False (AND) accept if ALL are true
                                 throw = cms.bool(True)                               
                                 )



primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                   vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                   minimumNDOF = cms.uint32(4) ,
                                   maxAbsZ = cms.double(24),	
                                   maxd0 = cms.double(2)	
                                   )

scrapping = cms.EDFilter("FilterOutScraping",
                         applyfilter = cms.untracked.bool(True),
                         debugOn = cms.untracked.bool(False),
                         numtrack = cms.untracked.uint32(10),
                         thresh = cms.untracked.double(0.25)
                         )

selectedPatJetsCountFilter = cms.EDFilter("PATCandViewCountFilter",
                                          minNumber = cms.uint32(0),
                                          maxNumber = cms.uint32(999999),
                                          src = cms.InputTag('selectedPatJetsPFlow')
                                          )

allEventsFilter = cms.EDFilter("AllEventsFilter")

noiseCleanedEventsFilter = allEventsFilter.clone()

import FWCore.ParameterSet.Config as cms


superClusters = cms.EDProducer("SuperClusterMerger",
   src = cms.VInputTag(cms.InputTag("hybridSuperClusters","", "RECO"),
                       cms.InputTag("multi5x5SuperClustersWithPreshower","", "RECO"))  
)

superClusterCands = cms.EDProducer("ConcreteEcalCandidateProducer",
   src = cms.InputTag("superClusters"),
   particleType = cms.int32(11),
)

goodSuperClusters = cms.EDFilter("CandViewSelector",
                                         src = cms.InputTag("superClusterCands"),
                                         cut = cms.string("et>10.0 && abs(eta)<2.5"),
                                         filter = cms.bool(False)
                                         )                                         
                                         

jetsToRemoveFromSuperCluster = cms.EDFilter("CaloJetSelector",   
    src = cms.InputTag("ak5CaloJets"),
    cut = cms.string('pt>10 && energyFractionHadronic > 0.15')
)

goodSuperClustersClean = cms.EDProducer("CandViewCleaner",
    srcCands = cms.InputTag("goodSuperClusters"),
    module_label = cms.string(''),
    srcObjects = cms.VInputTag(cms.InputTag("jetsToRemoveFromSuperCluster")),
    deltaRMin = cms.double(0.1)
)

makeSuperClusters = cms.Sequence(
    superClusters*
    superClusterCands*
    goodSuperClusters*
    jetsToRemoveFromSuperCluster*
    goodSuperClustersClean
    )

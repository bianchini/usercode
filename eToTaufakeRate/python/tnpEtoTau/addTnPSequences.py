import FWCore.ParameterSet.Config as cms



def addTnPSequences(process,makePAT="pat",makeMCtrees=True,makeUnbiased=True,removeWenuFilters=False):

    process.load("Bianchi.eToTaufakeRate.tnpEtoTau.tags_cff")
    process.load("Bianchi.eToTaufakeRate.tnpEtoTau.probes_cff")
    process.load("Bianchi.eToTaufakeRate.tnpEtoTau.tagAndProbes_cff")
    process.load("Bianchi.eToTaufakeRate.tnpEtoTau.mcMatches_cff")
    process.load("Bianchi.eToTaufakeRate.tnpEtoTau.oneTpFilters_cff")
    process.load("Bianchi.eToTaufakeRate.tnpEtoTau.filtersForTnP_cff")
    process.load("Bianchi.eToTaufakeRate.tnpEtoTau.treeProducers_cff")
    process.load("Bianchi.eToTaufakeRate.tnpEtoTau.EtoTauTnPSequences_cff")


    process.etoTauMargLoose90.isMC = cms.bool( makeMCtrees ) 
    process.etoTauMargMedium90.isMC = cms.bool( makeMCtrees ) 
    process.etoTauMargTight90.isMC = cms.bool( makeMCtrees )
    
    process.etoTauMargLoose80.isMC = cms.bool( makeMCtrees ) 
    process.etoTauMargMedium80.isMC = cms.bool( makeMCtrees ) 
    process.etoTauMargTight80.isMC = cms.bool( makeMCtrees )

    process.etoTauMargLoose70.isMC = cms.bool( makeMCtrees ) 
    process.etoTauMargMedium70.isMC = cms.bool( makeMCtrees ) 
    process.etoTauMargTight70.isMC = cms.bool( makeMCtrees )
    
    process.etoTauMargLoose60.isMC = cms.bool( makeMCtrees ) 
    process.etoTauMargMedium60.isMC = cms.bool( makeMCtrees ) 
    process.etoTauMargTight60.isMC = cms.bool( makeMCtrees )
    
    process.etoTauMargLooseNoCracks90.isMC = cms.bool( makeMCtrees ) 
    process.etoTauMargMediumNoCracks90.isMC = cms.bool( makeMCtrees ) 
    process.etoTauMargTightNoCracks90.isMC = cms.bool( makeMCtrees )
    
    process.etoTauMargLooseNoCracks80.isMC = cms.bool( makeMCtrees ) 
    process.etoTauMargMediumNoCracks80.isMC = cms.bool( makeMCtrees ) 
    process.etoTauMargTightNoCracks80.isMC = cms.bool( makeMCtrees )
    
    process.etoTauMargLooseNoCracks70.isMC = cms.bool( makeMCtrees ) 
    process.etoTauMargMediumNoCracks70.isMC = cms.bool( makeMCtrees ) 
    process.etoTauMargTightNoCracks70.isMC = cms.bool( makeMCtrees )
    
    process.etoTauMargLooseNoCracks60.isMC = cms.bool( makeMCtrees ) 
    process.etoTauMargMediumNoCracks60.isMC = cms.bool( makeMCtrees ) 
    process.etoTauMargTightNoCracks60.isMC = cms.bool( makeMCtrees )

    process.etoTauSCMarg90.isMC = cms.bool( makeMCtrees )
    process.etoTauSCMarg80.isMC = cms.bool( makeMCtrees )
    process.etoTauSCMarg70.isMC = cms.bool( makeMCtrees )
    process.etoTauSCMarg60.isMC = cms.bool( makeMCtrees )

    process.etoTauSCMargNoCracks90.isMC = cms.bool( makeMCtrees )
    process.etoTauSCMargNoCracks80.isMC = cms.bool( makeMCtrees )
    process.etoTauSCMargNoCracks70.isMC = cms.bool( makeMCtrees )
    process.etoTauSCMargNoCracks60.isMC = cms.bool( makeMCtrees )
    

    process.tagAndProbe90 = cms.Path(
        getattr(process,makePAT)*
        process.sequence90
        )
    process.tagAndProbe80 = cms.Path(
        getattr(process,makePAT)*
        process.sequence80
        )
    process.tagAndProbe70 = cms.Path(
        getattr(process,makePAT)*
        process.sequence70
        )
    process.tagAndProbe60 = cms.Path(
        getattr(process,makePAT)*
        process.sequence60
        )
    process.tagAndProbeSC90 = cms.Path(
        getattr(process,makePAT)*
        process.sequenceSC90
        )
    process.tagAndProbeSC80 = cms.Path(
        getattr(process,makePAT)*
        process.sequenceSC80
        )
    process.tagAndProbeSC70 = cms.Path(
        getattr(process,makePAT)*
        process.sequenceSC70
        )
    process.tagAndProbeSC60 = cms.Path(
        getattr(process,makePAT)*
        process.sequenceSC60
    )
    
    if removeWenuFilters:
        process.sequence90.remove(process.noWenu90)
        process.sequence80.remove(process.noWenu80)
        process.sequence70.remove(process.noWenu70)
        process.sequence60.remove(process.noWenu60)
        process.sequenceSC90.remove(process.noWenuSC90)
        process.sequenceSC80.remove(process.noWenuSC80)
        process.sequenceSC70.remove(process.noWenuSC70)
        process.sequenceSC60.remove(process.noWenuSC60)
        

    if not makeMCtrees:
        process.sequence90.remove(process.mcMatches90)
        process.sequence80.remove(process.mcMatches80)
        process.sequence70.remove(process.mcMatches70)
        process.sequence60.remove(process.mcMatches60)
        process.sequenceSC90.remove(process.mcMatchesSC90)
        process.sequenceSC80.remove(process.mcMatchesSC80)
        process.sequenceSC70.remove(process.mcMatchesSC70)
        process.sequenceSC60.remove(process.mcMatchesSC60)
            
      
    if makeUnbiased:
        process.oneTp90Loose.minNumber = 0
        process.oneTp80Loose.minNumber = 0
        process.oneTp70Loose.minNumber = 0
        process.oneTp60Loose.minNumber = 0
        process.oneTpSC90.minNumber = 0
        process.oneTpSC80.minNumber = 0
        process.oneTpSC70.minNumber = 0
        process.oneTpSC60.minNumber = 0

########################################################
########################################################

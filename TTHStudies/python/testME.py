import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

process = cms.Process("TestME")

process.fwliteInput = cms.PSet(

    outFileName   = cms.string("transferFunctions_SL.root"),
    pathToFile    = cms.string("./transferFunctions_SL.root"),
    vegasPoints   = cms.int32(400000),
    verbose       = cms.bool(False),
    
    
    )

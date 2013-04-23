import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

process = cms.Process("TFuncSL")

process.fwliteInput = cms.PSet(

    outFileName   = cms.string("transferFunctions_SL.root"),
    pathToFile    = cms.string("./topKinFitterSL_6j4b.root"),
    verbose       = cms.bool(False),
    
    
    )

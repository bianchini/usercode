import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

process = cms.Process("TestMENew")

process.fwliteInput = cms.PSet(

    outFileName   = cms.string("./root/TestMENew.root"),
    pathToFile    = cms.string("./root/transferFunctions_partonE.root"),
    vegasPoints   = cms.int32(500),

    verbose       = cms.bool(False),
    met           = cms.double(120.),

    pertBLep      = cms.double(1.0),

    pertW1        = cms.double(1.0),
    pertW2        = cms.double(1.0),
    pertBHad      = cms.double(1.0),

    enlargeE1     = cms.double(0.),
    enlargeEh1    = cms.double(0.),
    enlargePt     = cms.double(0.),

    )

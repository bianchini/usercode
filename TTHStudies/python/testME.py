import FWCore.ParameterSet.Types  as CfgTypes
import FWCore.ParameterSet.Config as cms

VType = "_VType2"

process = cms.Process("TestMENew")

process.fwliteInput = cms.PSet(

    outFileName   = cms.string("./root/TestMENew.root"),
    pathToTF      = cms.string("./root/transferFunctions_partonE.root"),
    pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt"+VType+"/v2/"),
    ordering      = cms.string("DiJetPt_"),
    lumi          = cms.double(12.1),

    samples       = cms.VPSet(
    
    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTH_HToBB_M-125_8TeV-pythia6'+VType),
    nickName = cms.string('TTH125'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1302*0.569)
    ),
    
    ),

    vegasPoints   = cms.int32(500),

    useME         = cms.int32(1),
    useJac        = cms.int32(1),
    useMET        = cms.int32(1),
    useTF         = cms.int32(1),
    usePDF        = cms.int32(1),
    
    verbose       = cms.bool(False),
    met           = cms.double(120.),
    evLimits      = cms.vint32(1,100),


    pertBLep      = cms.double(1.0),
    pertW1        = cms.double(1.0),
    pertW2        = cms.double(1.0),
    pertBHad      = cms.double(1.0),

    enlargeE1     = cms.double(0.),
    enlargeEh1    = cms.double(0.),
    enlargePt     = cms.double(0.),

    )

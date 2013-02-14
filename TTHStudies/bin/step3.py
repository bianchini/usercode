import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms



process = cms.Process("Step3")


process.fwliteInput = cms.PSet(

    #pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/ZllHDiJetPt/")
    #pathToFile    = cms.string("/scratch/bianchi/HBB_EDMNtuple/All.H.DiJetPt/"),
    pathToFile    = cms.string("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt/"),
    outPath       = cms.string("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_Step2/"),
    #ordering      = cms.string("ZllH.DiJetPt.Oct22."),
    ordering      = cms.string("DiJetPt_"),
    lumi          = cms.double(12.1),
    verbose       = cms.bool(False),
    

    samples       = cms.VPSet(

    cms.PSet(
    skip     = cms.bool(True),
    name     = cms.string('DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball'),
    nickName = cms.string('DYJets'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71),
    update   = cms.bool(True)
    ),

    
    ),


    )

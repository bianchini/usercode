import FWCore.ParameterSet.Config as cms


hltMatching = cms.string('(triggerObjectMatchesByPath("HLT_Ele10_SW_L1R").size()!=0 || triggerObjectMatchesByPath("HLT_Ele10_LW_L1R").size()!=0 || triggerObjectMatchesByPath("HLT_Ele10_SW_EleId_L1R").size()!=0 || triggerObjectMatchesByPath("HLT_Ele10_LW_EleId_L1R").size()!=0 ||  triggerObjectMatchesByPath("HLT_Ele12_SW_TightEleIdIsol_L1R").size()!=0 ||  triggerObjectMatchesByPath("HLT_Ele12_SW_TighterEleIdIsol_L1R").size()!=0 || triggerObjectMatchesByPath("HLT_Ele15_SW_L1R").size()!=0 || triggerObjectMatchesByPath("HLT_Ele15_SW_CaloEleId_L1R").size()!=0 || triggerObjectMatchesByPath("HLT_Ele17_SW_CaloEleId_L1R").size()!=0 || triggerObjectMatchesByPath("HLT_Ele17_SW_TightEleId_L1R").size()!=0 || triggerObjectMatchesByPath("HLT_Ele17_SW_TightEleId_L1R").size()!=0 || triggerObjectMatchesByPath("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2").size()!=0 || triggerObjectMatchesByPath("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3").size()!=0)')

'''
    * HLT_Ele10_LW_L1R for runs [136033,139980]           (~0.1/pb)
    * HLT_Ele15_SW_L1R for runs [140058,141882]           (~0.23/pb)
    * HLT_Ele15_SW_CaloEleId_L1R for runs [141956,144114] (~3/pb)
    * HLT_Ele17_SW_CaloEleId_L1R for runs [146428,147116] (~7/pb)
    * HLT_Ele17_SW_TightEleId_L1R for runs [147196,148058]
    * HLT_Ele17_SW_TighterEleIdIsol_L1R_v2 for runs [148819,149064]
    * HLT_Ele17_SW_TighterEleIdIsol_L1R_v3 for runs [149181,149442]
'''

#hltMatching = cms.string('pt>0 ')

tagAnyElec = cms.EDFilter("PATElectronRefSelector",
                          src = cms.InputTag("selectedPatElectronsTriggerMatch"),
                          cut = cms.string(hltMatching.value()+' && pt>25.0 && abs(eta)<2.5 && !(1.4442< abs(eta) <1.566)'),
                          filter = cms.bool(False)
                          )

tag90 = tagAnyElec.clone(
    cut = cms.string(tagAnyElec.cut.value()+' && electronID("simpleEleId90relIso") > 6.5'),
    )
tag80 = tagAnyElec.clone(
    cut = cms.string(tagAnyElec.cut.value()+' && electronID("simpleEleId80relIso") > 6.5'),
    )
tag70 = tagAnyElec.clone(
    cut = cms.string(tagAnyElec.cut.value()+' && electronID("simpleEleId70relIso") > 6.5'),
    )
tag60 = tagAnyElec.clone(
    cut = cms.string(tagAnyElec.cut.value()+' && electronID("simpleEleId60relIso") > 6.5'),
    )

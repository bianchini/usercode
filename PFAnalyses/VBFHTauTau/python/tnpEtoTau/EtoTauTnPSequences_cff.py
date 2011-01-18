import FWCore.ParameterSet.Config as cms


## main sequence
#from PFAnalyses.Z.EtoTauTemplate_cff import *

### tags
from PFAnalyses.VBFHTauTau.tnpEtoTau.tags_cff import *

### probes
from PFAnalyses.VBFHTauTau.tnpEtoTau.probes_cff import *

### tag&probes
from PFAnalyses.VBFHTauTau.tnpEtoTau.tagAndProbes_cff import *

### matches
from PFAnalyses.VBFHTauTau.tnpEtoTau.mcMatches_cff import *

### one tag&probe filters
from PFAnalyses.VBFHTauTau.tnpEtoTau.oneTpFilters_cff import *

### anti-Wenu filters
from PFAnalyses.VBFHTauTau.tnpEtoTau.filtersForTnP_cff import *

### tree producers
from PFAnalyses.VBFHTauTau.tnpEtoTau.treeProducers_cff import *


tightestMatchedID = cms.EDProducer("MatchedIDComputer",
                                   probes = cms.InputTag("patTausPFlow"),
                                   electrons = cms.InputTag("selectedPatElectronsTriggerMatch")
                                   )
tightestMatchedIDSC = cms.EDProducer("MatchedIDComputer",
                                     probes = cms.InputTag("patTaus"),
                                     electrons = cms.InputTag("selectedPatElectronsTriggerMatch")
                                     )

sequence90 = cms.Sequence((tag90+
                           probeTauIDIsoLoose+
                           probeTauIDIsoMedium+
                           probeTauIDIsoTight+
                           probeTauIDIsoLooseNoCracks+
                           probeTauIDIsoMediumNoCracks+
                           probeTauIDIsoTightNoCracks+
                           tightestMatchedID
                           )*
                          mcMatches90*
                          (tnp90TauIDLoose+
                           tnp90TauIDMedium+
                           tnp90TauIDTight+
                           tnp90TauIDLooseNoCracks+
                           tnp90TauIDMediumNoCracks+
                           tnp90TauIDTightNoCracks
                           )*
                          oneTp90Loose*
                          noWenu90*
                          (etoTauMargLoose90+
                           etoTauMargMedium90+
                           etoTauMargTight90+
                           etoTauMargLooseNoCracks90+
                           etoTauMargMediumNoCracks90+
                           etoTauMargTightNoCracks90
                           )
                          )

sequenceSC90 = cms.Sequence((tag90+
                             probeTauSCIDIso+
                             probeTauSCIDIsoNoCracks+
                             tightestMatchedIDSC
                             )*
                            mcMatchesSC90*
                            (tnp90TauSCIDIso+
                             tnp90TauSCIDIsoNoCracks
                             )*
                            oneTpSC90*
                            noWenuSC90*
                            (etoTauSCMarg90+
                             etoTauSCMargNoCracks90
                             )
                            )

sequence80 = cms.Sequence((tag80+
                           probeTauIDIsoLoose+
                           probeTauIDIsoMedium+
                           probeTauIDIsoTight+
                           probeTauIDIsoLooseNoCracks+
                           probeTauIDIsoMediumNoCracks+
                           probeTauIDIsoTightNoCracks+
                           tightestMatchedID
                           )*
                          mcMatches80*
                          (tnp80TauIDLoose+
                           tnp80TauIDMedium+
                           tnp80TauIDTight+
                           tnp80TauIDLooseNoCracks+
                           tnp80TauIDMediumNoCracks+
                           tnp80TauIDTightNoCracks
                           )*
                          oneTp80Loose*
                          noWenu80*
                          (etoTauMargLoose80+
                           etoTauMargMedium80+
                           etoTauMargTight80+
                           etoTauMargLooseNoCracks80+
                           etoTauMargMediumNoCracks80+
                           etoTauMargTightNoCracks80
                           )
                          )

sequenceSC80 = cms.Sequence((tag80+
                             probeTauSCIDIso+
                             probeTauSCIDIsoNoCracks+
                             tightestMatchedIDSC
                             )*
                            mcMatchesSC80*
                            (tnp80TauSCIDIso+
                             tnp80TauSCIDIsoNoCracks
                             )*
                            oneTpSC80*
                            noWenuSC80*
                            (etoTauSCMarg80+
                             etoTauSCMargNoCracks80
                             )
                            )

sequence70 = cms.Sequence((tag70+
                           probeTauIDIsoLoose+
                           probeTauIDIsoMedium+
                           probeTauIDIsoTight+
                           probeTauIDIsoLooseNoCracks+
                           probeTauIDIsoMediumNoCracks+
                           probeTauIDIsoTightNoCracks+
                           tightestMatchedID
                           )*
                          mcMatches70*
                          (tnp70TauIDLoose+
                           tnp70TauIDMedium+
                           tnp70TauIDTight+
                           tnp70TauIDLooseNoCracks+
                           tnp70TauIDMediumNoCracks+
                           tnp70TauIDTightNoCracks
                           )*
                          oneTp70Loose*
                          noWenu70*
                          (etoTauMargLoose70+
                           etoTauMargMedium70+
                           etoTauMargTight70+
                           etoTauMargLooseNoCracks70+
                           etoTauMargMediumNoCracks70+
                           etoTauMargTightNoCracks70
                           )
                          )

sequenceSC70 = cms.Sequence((tag70+
                             probeTauSCIDIso+
                             probeTauSCIDIsoNoCracks+
                             tightestMatchedIDSC
                             )*
                            mcMatchesSC70*
                            (tnp70TauSCIDIso+
                             tnp70TauSCIDIsoNoCracks
                             )*
                            oneTpSC70*
                            noWenuSC70*
                            (etoTauSCMarg70+
                             etoTauSCMargNoCracks70
                             )
                            )

sequence60 = cms.Sequence((tag60+
                           probeTauIDIsoLoose+
                           probeTauIDIsoMedium+
                           probeTauIDIsoTight+
                           probeTauIDIsoLooseNoCracks+
                           probeTauIDIsoMediumNoCracks+
                           probeTauIDIsoTightNoCracks+
                           tightestMatchedID
                           )*
                          mcMatches60*
                          (tnp60TauIDLoose+
                           tnp60TauIDMedium+
                           tnp60TauIDTight+
                           tnp60TauIDLooseNoCracks+
                           tnp60TauIDMediumNoCracks+
                           tnp60TauIDTightNoCracks
                           )*
                          oneTp60Loose*
                          noWenu60*
                          (etoTauMargLoose60+
                           etoTauMargMedium60+
                           etoTauMargTight60+
                           etoTauMargLooseNoCracks60+
                           etoTauMargMediumNoCracks60+
                           etoTauMargTightNoCracks60
                           )
                          )

sequenceSC60 = cms.Sequence((tag60+
                             probeTauSCIDIso+
                             probeTauSCIDIsoNoCracks+
                             tightestMatchedIDSC
                             )*
                            mcMatchesSC60*
                            (tnp60TauSCIDIso+
                             tnp60TauSCIDIsoNoCracks
                             )*
                            oneTpSC60*
                            noWenuSC60*
                            (etoTauSCMarg60+
                             etoTauSCMargNoCracks60
                             )
                            )



// system include files
#include <memory>
#include <algorithm>
#include <sstream>

// user include files
#include "PFAnalyses/VBFHTauTau/interface/JetVetoAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"


#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "AnalysisDataFormats/PFAnalyses/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/PFAnalyses/interface/CompositePtrCandidateT1T2MEtFwd.h"
#include "PFAnalyses/CommonTools/plugins/VBFEventConcreteProducers.h"

#include "PFAnalyses/CommonTools/interface/FWLiteTreeAnalyzer.h"
#include "PFAnalyses/VBFHTauTau/interface/FWLiteTriggerAnalyzer.h"
#include "PFAnalyses/VBFHTauTau/interface/FWLiteJetVetoAnalyzer.h"
#include "PFAnalyses/VBFHTauTau/interface/FWLiteDiTauAnalyzer.h"
#include "PFAnalyses/VBFHTauTau/interface/FWLiteDiTauAnalyzer2.h"
#include "PFAnalyses/VBFHTauTau/interface/FWLiteEventMixer.h"

#include <string>
//
// constructors and destructor
//
JetVetoAnalyzer::JetVetoAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  std::vector<FWLiteAnalyzer*> myAnalyzers;
  myAnalyzers.push_back(new FWLiteTriggerAnalyzer("TriggerAnalyzer"));
  //myAnalyzers.push_back(new FWLiteEventMixer("EventMixer"));
  myAnalyzers.push_back(new FWLiteDiTauAnalyzer2("DiTau"));

  std::string cfgFileName = iConfig.getUntrackedParameter<std::string>("cfgFileName");
  fwLiteTreeAnalyzer_ =  new FWLiteTreeAnalyzer("TreeAnalyzer",cfgFileName);
  fwLiteTreeAnalyzer_->init(myAnalyzers);

}


JetVetoAnalyzer::~JetVetoAnalyzer(){
 
  fwLiteTreeAnalyzer_->finalize();
  delete fwLiteTreeAnalyzer_;

}


//
// member functions
//

// ------------ method called to for each event  ------------
bool
JetVetoAnalyzer::filter(edm::Event& iEvent, edm::EventSetup const& iSetup)
{
   using namespace edm;
   using namespace reco;

   fwLiteTreeAnalyzer_->analyze(iEvent);
   const pat::strbitset & mySelections = fwLiteTreeAnalyzer_->getSelections();   

   bool decision = true;
   for(unsigned int i=0;i<mySelections.strings().size();++i){
     if(mySelections.strings()[i].find("PF")!=std::string::npos) 
       decision&=mySelections.test(mySelections.strings()[i]);
   }

   return decision;
}


// ------------ method called once each run just before starting event loop  ------------
void JetVetoAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){ }

// ------------ method called once each job just before starting event loop  ------------
void JetVetoAnalyzer::beginJob(){ }

// ------------ method called once each job just after ending the event loop  ------------
void JetVetoAnalyzer::endJob() { }

//define this as a plug-in
DEFINE_FWK_MODULE(JetVetoAnalyzer);

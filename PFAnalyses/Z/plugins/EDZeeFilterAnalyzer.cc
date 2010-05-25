// system include files
#include <memory>
#include <algorithm>
#include <sstream>

// user include files
#include "PFAnalyses/Z/plugins/EDZeeFilterAnalyzer.h"
#include "PFAnalyses/Z/interface/PatZeeAnalyzer.h"
#include "PFAnalyses/CommonTools/interface/FWLiteTreeAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>
//
// constructors and destructor
//
EDZeeFilterAnalyzer::EDZeeFilterAnalyzer(const edm::ParameterSet& iConfig){
   //now do what ever initialization is needed


 std::vector<FWLiteAnalyzer*> myAnalyzers;

  myAnalyzers.push_back(new PatZeeAnalyzer("PatZeeAnalyzer")); 
  
  std::string cfgFileName = iConfig.getUntrackedParameter<std::string>("cfgFileName");

  fwLiteTreeAnalyzer_ =  new FWLiteTreeAnalyzer("TreeAnalyzer",cfgFileName);
  fwLiteTreeAnalyzer_->init(myAnalyzers);

}


EDZeeFilterAnalyzer::~EDZeeFilterAnalyzer(){
 
  fwLiteTreeAnalyzer_->finalize();
  delete fwLiteTreeAnalyzer_;

}


//
// member functions
//

// ------------ method called to for each event  ------------
bool
EDZeeFilterAnalyzer::filter(edm::Event& iEvent, edm::EventSetup const& iSetup)
{
   using namespace edm;
   using namespace reco;

   fwLiteTreeAnalyzer_->analyze(iEvent);
   const std::strbitset & mySelections = fwLiteTreeAnalyzer_->getSelections();   

   bool decision = true;
   for(unsigned i=0;i<mySelections.strings().size();++i){
     if(mySelections.strings()[i].find("leading electron p_{T} cut")!=std::string::npos ||
	mySelections.strings()[i].find("leading electron #eta cut")!=std::string::npos ||
	mySelections.strings()[i].find("leading electron loose isolated")!=std::string::npos )  
       decision&=mySelections.test(mySelections.strings()[i]);
   }

   return decision;
}


// ------------ method called once each run just before starting event loop  ------------
void EDZeeFilterAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){ }

// ------------ method called once each job just before starting event loop  ------------
void EDZeeFilterAnalyzer::beginJob(){ }

// ------------ method called once each job just after ending the event loop  ------------
void EDZeeFilterAnalyzer::endJob() { }

//define this as a plug-in
DEFINE_FWK_MODULE(EDZeeFilterAnalyzer);

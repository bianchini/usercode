#ifndef PFAnalyses_CommonTools_JetVetoAnalyzer_h
#define PFAnalyses_CommonTools_JetVetoAnalyzer_h

//
//
// Original Author:  Artur Kalinowski
//         Created:  Mon Mar 22 13:25:58 CET 2010
//
//


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"


class FWLiteTreeAnalyzer;
//
// class decleration
//

class EDZeeFilterAnalyzer: public edm::EDFilter {
   public:
      explicit EDZeeFilterAnalyzer(const edm::ParameterSet&);
      ~EDZeeFilterAnalyzer();


   private:
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, edm::EventSetup const&);
      virtual void endJob() ;

      // ----------member data ---------------------------

      FWLiteTreeAnalyzer* fwLiteTreeAnalyzer_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//


#endif

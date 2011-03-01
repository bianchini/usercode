#ifndef Bianchi_TauTauStudies_MuonsUserEmbedded_h
#define Bianchi_TauTauStudies_MuonsUserEmbedded_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

class MuonsUserEmbedded : public edm::EDProducer{


 public: 

  explicit MuonsUserEmbedded(const edm::ParameterSet&);
  virtual ~MuonsUserEmbedded();

 private:

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  edm::InputTag muonTag_;
  edm::InputTag vertexTag_;

};


#endif

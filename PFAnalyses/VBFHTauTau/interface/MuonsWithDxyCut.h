#ifndef PFAnalyses_VBFHTauTau_MuonsWithDxyCut_h
#define PFAnalyses_VBFHTauTau_MuonsWithDxyCut_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

class MuonsWithDxyCut : public edm::EDProducer{


 public: 

  explicit MuonsWithDxyCut(const edm::ParameterSet&);
  ~MuonsWithDxyCut();

 private:

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  float dxyCut_;
  edm::InputTag muonTag_;
  edm::InputTag vertexTag_;

};


#endif

#include "PFAnalyses/VBFHTauTau/interface/MuonsWithDxyCut.h"


MuonsWithDxyCut::MuonsWithDxyCut(const edm::ParameterSet & iConfig){

  dxyCut_ = iConfig.getParameter<double>("dxyCut");
  muonTag_ = iConfig.getParameter<edm::InputTag>("muonTag");
  vertexTag_ = iConfig.getParameter<edm::InputTag>("vertexTag");
 
  produces<pat::MuonCollection>("");

}

MuonsWithDxyCut::~MuonsWithDxyCut(){
}

void MuonsWithDxyCut::produce(edm::Event & iEvent, const edm::EventSetup & iSetup){

  edm::Handle<pat::MuonCollection> muonsHandle;
  iEvent.getByLabel(muonTag_,muonsHandle);
  const pat::MuonCollection* muons = muonsHandle.product();

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByLabel(vertexTag_,vertexHandle);
  const reco::VertexCollection* vertexes = vertexHandle.product();

  std::auto_ptr< pat::MuonCollection > muonsWithDxyCut( new pat::MuonCollection() ) ;


  if( vertexes->size()!=0){
    for(unsigned int i = 0; i < muons->size(); i++){
      pat::Muon aMuon( (*muons)[i] );
      if( (aMuon.globalTrack())->dxy( (*vertexes)[0].position() ) < dxyCut_ ) 
	muonsWithDxyCut->push_back(aMuon);
    }
  } 
  else{
    for(unsigned int i = 0; i < muons->size(); i++){
      pat::Muon aMuon( (*muons)[i] );
      muonsWithDxyCut->push_back(aMuon);
    }

  }


  iEvent.put( muonsWithDxyCut );
  return;
}


#include "FWCore/Framework/interface/MakerMacros.h"
 
DEFINE_FWK_MODULE(MuonsWithDxyCut);



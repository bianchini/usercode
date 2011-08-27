#include "Bianchi/TauTauStudies/interface/MEtRescalerProducer.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


#include <cstdio>
#include <map>

using namespace std;

MEtRescalerProducer::MEtRescalerProducer(const edm::ParameterSet & iConfig){

  metTag_  = iConfig.getParameter<edm::InputTag>("metTag");
  jetTag_  = iConfig.getParameter<edm::InputTag>("jetTag");
  electronTag_ = iConfig.getParameter<edm::InputTag>("electronTag");
  muonTag_   = iConfig.getParameter<edm::InputTag>("muonTag");
  tauTag_  = iConfig.getParameter<edm::InputTag>("tauTag");

  unClusterShift_  = iConfig.getParameter<double>("unClusterShift");
  tauShift_   = iConfig.getParameter<double>("tauShift");
  muonShift_  = iConfig.getParameter<double>("muonShift");
  electronShift_  = iConfig.getParameter<double>("electronShift");
  jetThreshold_    = iConfig.getParameter<double>("jetThreshold");

  verbose_    = iConfig.getParameter<bool>("verbose");
  
 
  // convention: jet - elec - mu - tau - uncluster
  char shift[] = {'U','D'};

  for(int i = 0; i < 2 ; i++){
    for(int j = 0; j < 2 ; j++){
      for(int k = 0; k < 2 ; k++){
	for(int l = 0; l < 2 ; l++){
	  for(int m = 0; m < 2 ; m++){
	    char buffer[5];
	    sprintf (buffer, "%c%c%c%c%c", 
		     shift[i],shift[j],shift[k],shift[l],shift[m]);
	    produces<pat::METCollection>(buffer);     
	  }
	}
      }
    }
  }

}

MEtRescalerProducer::~MEtRescalerProducer(){
}

void MEtRescalerProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup){

  edm::Handle<pat::METCollection> metHandle;
  iEvent.getByLabel(metTag_,metHandle);
  const pat::METCollection* mets = metHandle.product();

  edm::Handle<pat::JetCollection> jetHandle;
  iEvent.getByLabel(jetTag_,jetHandle);
  const pat::JetCollection* jets = jetHandle.product();

  edm::Handle<pat::MuonCollection> muonHandle;
  const pat::MuonCollection* muons = 0;
  if(muonTag_.label()!=""){
    iEvent.getByLabel(muonTag_,muonHandle);
    muons = muonHandle.product();
  }

  edm::Handle<pat::ElectronCollection> electronHandle;
  const pat::ElectronCollection* electrons = 0;
  if(electronTag_.label()!=""){
    iEvent.getByLabel(electronTag_,electronHandle);
    electrons = electronHandle.product();
  }

  edm::Handle<pat::TauCollection> tauHandle;
  iEvent.getByLabel(tauTag_,tauHandle);
  const pat::TauCollection* taus = tauHandle.product();

  // convention: jet - elec - mu - tau - uncluster
  char shift[] = {'U','D'};

  for(int i = 0; i < 2 ; i++){
    for(int j = 0; j < 2 ; j++){
      for(int k = 0; k < 2 ; k++){
	for(int l = 0; l < 2 ; l++){
	  for(int m = 0; m < 2 ; m++){
	    char buffer[5];
	    sprintf (buffer, "%c%c%c%c%c", 
		     shift[i],shift[j],shift[k],shift[l],shift[m]);

	    if(verbose_) printf("Now doing shift: jet=%c ele=%c mu=%c tau=%c uncluster=%c\n",
				shift[i],shift[j],shift[k],shift[l],shift[m]); 
	    float dPx = 0; float dPy=0; float dSumEt=0;
	    
	    math::XYZTLorentzVectorD unClusterEnergy = -(*mets)[0].p4();
	    if(verbose_){
	      cout << "Unclustered energy (Et,phi) " << (-1)*unClusterEnergy.Et() << ", "
		   << unClusterEnergy.Phi() << endl;
	    }

	    // electron
	    for(unsigned int it = 0; electrons!=0 && it<electrons->size(); it++){
	      if(j==0) dPx += (*electrons)[it].px() * electronShift_;
	      else     dPx -= (*electrons)[it].px() * electronShift_;
	      if(j==0) dPy += (*electrons)[it].py() * electronShift_;
	      else     dPy -= (*electrons)[it].py() * electronShift_;
	      if(j==0) dSumEt += (*electrons)[it].pt() * electronShift_;
	      else     dSumEt -= (*electrons)[it].pt() * electronShift_;

	      unClusterEnergy -= (*electrons)[it].p4();
	      if(verbose_){
		cout << "Shifting electron #" << it << " => dPx= " << dPx
		     << " , dPy= " << dPy
		     << " , dSumEt= " << dSumEt
		     << endl;
	      }
	    }
	    // muon
	    for(unsigned int it = 0; muons!=0 && it<muons->size(); it++){
	      if(k==0) dPx += (*muons)[it].px() * muonShift_;
	      else     dPx -= (*muons)[it].px() * muonShift_;
	      if(k==0) dPy += (*muons)[it].py() * muonShift_;
	      else     dPy -= (*muons)[it].py() * muonShift_;
	      if(k==0) dSumEt += (*muons)[it].pt() * muonShift_;
	      else     dSumEt -= (*muons)[it].pt() * muonShift_;

	      unClusterEnergy -= (*muons)[it].p4();
	      if(verbose_){
		cout << "Shifting muon #" << it << " => dPx= " << dPx
		     << " , dPy= " << dPy
		     << " , dSumEt= " << dSumEt
		     << endl;
	      }
	    }
	    // tau
	    for(unsigned int it = 0; taus!=0 && it<taus->size(); it++){
	      if(l==0) dPx += (*taus)[it].px() * tauShift_;
	      else     dPx -= (*taus)[it].px() * tauShift_;
	      if(l==0) dPy += (*taus)[it].py() * tauShift_;
	      else     dPy -= (*taus)[it].py() * tauShift_;
	      if(l==0) dSumEt += (*taus)[it].pt() * tauShift_;
	      else     dSumEt -= (*taus)[it].pt() * tauShift_;

	      unClusterEnergy -= (*taus)[it].p4();
	      if(verbose_){
		cout << "Shifting tau #" << it << " => dPx= " << dPx
		     << " , dPy= " << dPy
		     << " , dSumEt= " << dSumEt
		     << endl;
	      }
	    }
	    
	    // jet
	    for(unsigned int it = 0; it<jets->size(); it++){

	      bool overlaps = (*jets)[it].correctedJet("Uncorrected","none","patJetCorrFactors").pt() < jetThreshold_ ;
	      for(unsigned int jt = 0; electrons!=0 && jt<electrons->size(); jt++){
		if(Geom::deltaR((*jets)[it].p4(),(*electrons)[jt].p4())<0.5)
		  overlaps = true;
	      }
	      for(unsigned int jt = 0; muons!=0 && jt<muons->size(); jt++){
		if(Geom::deltaR((*jets)[it].p4(),(*muons)[jt].p4())<0.5)
		  overlaps = true;
	      }
	      for(unsigned int jt = 0; taus!=0 && jt<taus->size(); jt++){
		if(Geom::deltaR((*jets)[it].p4(),(*taus)[jt].p4())<0.5)
		  overlaps = true;
	      }
	      if( overlaps && verbose_) cout << "Jet overlaps! Do not double-count it!" << endl;
	      if(!overlaps){

		//// try to use JES uncertainties
		edm::ESHandle<JetCorrectorParametersCollection> jetCorrParameters;
		// get the jet corrector parameters collection from the global tag
		iSetup.get<JetCorrectionsRecord>().get("AK5PF", jetCorrParameters);
		// get the uncertainty parameters from the collection
		JetCorrectorParameters const & param = (*jetCorrParameters)["Uncertainty"];
		// instantiate the jec uncertainty object
		JetCorrectionUncertainty* deltaJEC = new JetCorrectionUncertainty(param);
		deltaJEC->setJetEta((*jets)[it].eta()); deltaJEC->setJetPt((*jets)[it].pt());
		float jetShift  = deltaJEC->getUncertainty( true );
		if(i==0) dPx += (*jets)[it].px() * jetShift;
		else     dPx -= (*jets)[it].px() * jetShift;
		if(i==0) dPy += (*jets)[it].py() * jetShift;
		else     dPy -= (*jets)[it].py() * jetShift;
		if(i==0) dSumEt += (*jets)[it].pt() * jetShift;
		else     dSumEt -= (*jets)[it].pt() * jetShift;

		unClusterEnergy -= (*jets)[it].p4();
		if(verbose_){
		  cout << "Shifting jet #" << it << " => dPx= " << dPx
		       << " , dPy= " << dPy
		       << " , dSumEt= " << dSumEt
		       << endl;
		}
	      }
 
	    }

	    // uncluster energy
	    if(m==0) dPx += unClusterEnergy.Px() * unClusterShift_;
	    else     dPx -= unClusterEnergy.Px() * unClusterShift_;
	    if(m==0) dPy += unClusterEnergy.Py() * unClusterShift_;
	    else     dPy -= unClusterEnergy.Py() * unClusterShift_;
	    if(m==0) dSumEt += (-1)*unClusterEnergy.Pt() * unClusterShift_;
	    else     dSumEt -=  (-1)*unClusterEnergy.Pt() * unClusterShift_;
	    if(verbose_){
	      cout << "Shifting uncluster energy " << " => dPx= " << dPx
		   << " , dPy= " << dPy
		   << " , dSumEt= " << dSumEt
		   << endl;
	    }

	    double scaledMETPx = (*mets)[0].px() - dPx;
	    double scaledMETPy = (*mets)[0].py() - dPy;
	    pat::MET scaledMET( reco::MET(((*mets)[0]).sumEt()+dSumEt, 
				reco::MET::LorentzVector(scaledMETPx, scaledMETPy, 0, sqrt(scaledMETPx*scaledMETPx+scaledMETPy*scaledMETPy)), 
					  reco::MET::Point(0,0,0)) );

	    if(verbose_){
	      cout << "Original MET (Et,phi) " <<  (*mets)[0].pt() << ", " << (*mets)[0].phi()
		   << " --- Rescaled MET " << scaledMET.pt() << ", " << scaledMET.phi() << endl;
	    }
	    std::auto_ptr< pat::METCollection > MEtRescaledColl( new pat::METCollection() ) ;
	    MEtRescaledColl->push_back(scaledMET);

	    iEvent.put(MEtRescaledColl, buffer);
	  }
	}
      }
    }
  }
  
 return;
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(MEtRescalerProducer);



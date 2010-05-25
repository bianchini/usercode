#include "PFAnalyses/CommonTools/interface/FWLiteTreeAnalyzer.h"
#include "PFAnalyses/CommonTools/interface/PatLeptonSelector.h"

#include "PFAnalyses/CommonTools/interface/ElectronHistograms.h"
#include "PFAnalyses/CommonTools/interface/IsolationHistograms.h"

#include "PFAnalyses/Z/interface/PatZeeAnalyzer.h"
#include "PFAnalyses/Z/interface/PatZeeHistograms.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Utilities/interface/Algorithms.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"

#include <TMath.h>
#include <bitset>
#include <map>

using namespace std;
using namespace reco;

PatZeeAnalyzer::PatZeeAnalyzer(const std::string & aName):FWLiteAnalyzer(aName){}


PatZeeAnalyzer::~PatZeeAnalyzer()
{
  std::cout<<"PatZeeAnalyzer::~PatZeeAnalyzer()"<<std::endl;
  
  delete zeeHistos_;
  delete leadEleHistos_;
  delete subLeadEleHistos_;
  delete leadEleIsoHistos_;
  delete subLeadEleIsoHistos_;
}


void PatZeeAnalyzer::initialize(const edm::ParameterSet& ps, 
				     TFileDirectory& aDir,
				     std::strbitset *aSelections){

  // offline vertexes
  vtxLabel_ = ps.getParameter<edm::InputTag>("offlinePrimaryVertecesLabel");

  // pat electrons
  patElectronsLabel_ =  ps.getParameter<edm::InputTag>("patElectronsLabel");

  // trigger results
  triggerResultsLabel_ = ps.getParameter<edm::InputTag>("triggerResultsLabel");
  triggerItemNames_ = ps.getParameter< std::vector<std::string> >("triggerItemNames"); 
  useTrigger_ = ps.getParameter<bool>("useTrigger");
  triggerItemAND_ = ps.getParameter<bool>("triggerItemAND");

  // initilaize the object selectors
  elePhaseSpaceSelector_.initialize( ps.getParameter<edm::ParameterSet>("elePhaseSpaceSelector") ); 
  eleSelector_.initialize( ps.getParameter<edm::ParameterSet>("eleSelector") ); 
  eleIsoSelector_.initialize( ps.getParameter<edm::ParameterSet>("eleIsoSelector") ); 
  
  // initilaize the parameters
  minPtSubSub_ = ps.getParameter<double>("minPtSubSub");
  verbose_ = ps.getParameter<bool>("verbose");
  triggerVerbose_ = ps.getParameter<bool>("triggerVerbose");
  scanEvsByM_ = ps.getParameter<bool>("scanEvsByM");
  minMass_ = ps.getParameter<double>("minMass");
  maxMass_ = ps.getParameter<double>("maxMass");

  // class to store the Zee specific histos
  zeeHistos_ = new PatZeeHistograms(&aDir,"Zee");

  // classes to store the electrons specific histos:
  leadEleHistos_ = new ElectronHistograms(&aDir,"LeadingElectron");
  subLeadEleHistos_ = new ElectronHistograms(&aDir,"SubLeadingElectron");

  // classes to store the electrons isolation cuts histos:
  leadEleIsoHistos_ = new IsolationHistograms(&aDir,"LeadingIsolation");
  subLeadEleIsoHistos_ = new IsolationHistograms(&aDir,"SubLeadingIsolation");

  mySelections_ = aSelections;

  registerCuts();

  if(verbose_) mySelections_->print(std::cout);

}



void PatZeeAnalyzer::registerCuts(){

  std::string labels[1] = {""};
  for(int i =0;i<1;++i){  
    mySelections_->push_back("HLT passed"+labels[i]);
    mySelections_->push_back("> 1 pat electrons"+labels[i]);
    mySelections_->push_back("> 1 sorted electrons"+labels[i]);
    mySelections_->push_back("leading electron p_{T} cut"+labels[i]);
    mySelections_->push_back("leading electron #eta cut"+labels[i]);
    mySelections_->push_back("subleading electron p_{T} cut"+labels[i]);
    mySelections_->push_back("subleading electron #eta cut"+labels[i]);
    mySelections_->push_back("2 accepted electrons"+labels[i]);
    mySelections_->push_back("subsubleading electron !accepted"+labels[i]);
    mySelections_->push_back("opposite charge"+labels[i]);
    mySelections_->push_back("leading electron loose isolated"+labels[i]);
    //mySelections_->push_back("leading electron hard isolated"+labels[i]);
    mySelections_->push_back("subleading electron loose isolated"+labels[i]);
    //mySelections_->push_back("subleading electron hard isolated"+labels[i]);
    mySelections_->push_back("m(ee) constrain"+labels[i]);
  }
}



bool PatZeeAnalyzer::analyze(const edm::EventBase& iEvent){

  if(verbose_) std::cout<<iEvent.id()<<std::endl;

  using namespace reco;

  float eventWeight = 1.0;

  clear();

  edm::Handle<pat::ElectronCollection> patElectrons;
  edm::Handle<reco::VertexCollection> recVtxs;
  edm::Handle<edm::TriggerResults> triggerResults_;

  try{
    iEvent.getByLabel(patElectronsLabel_,patElectrons);
    //iEvent.getByLabel(vtxLabel_,recVtxs);
  }
  catch(...){
    std::cout<<"pat::Electron collection cannot be found!";
    return false;
  };

  mySelections_->set("> 1 pat electrons", (*patElectrons).size()>1 );

  // trigger selections
  if(useTrigger_){

    iEvent.getByLabel(triggerResultsLabel_,triggerResults_);
    if(!triggerResults_.isValid()) return false;    

    fwlite::Event& fwliteEvent = (fwlite::Event&) iEvent;
    edm::TriggerNames const&  triggerNames = fwliteEvent.triggerNames(*triggerResults_);

    
    if( triggerVerbose_ ){
      std::vector<std::string> names = triggerNames.triggerNames();
      for(unsigned i = 0; i < names.size();i++){
	int index =  triggerNames.triggerIndex(names[i]);
	if( index == triggerNames.size()) continue;
	bool iaccept = triggerResults_->accept( index );
	const int fired = iaccept ? 1 : 0;
	std::cout<< names[i] << ":   " << fired << std::endl;
      }
      std::cout << "**** end of the event ****" << std::endl;
    }
    

    bool triggerPassed = false;
    unsigned goodTriggers = 0;
    for(unsigned i = 0; i<triggerItemNames_.size();i++){
      int triggerItemIndex =  triggerNames.triggerIndex(triggerItemNames_[i]);
      if( triggerItemIndex == triggerNames.size() ){
	std::cout << "*** Trigger name not found ***" << std::endl;
	continue;
      }
      goodTriggers++;
      bool iaccept = triggerResults_->accept(triggerItemIndex);
      if(triggerVerbose_ ){
	const int fired = iaccept ? 1 : 0;
	std::cout << triggerItemNames_[i] 
		  << " (" <<triggerItemIndex 
		  << ")" << ":" << fired 
		  <<std::endl;
      }
      if(triggerItemAND_){
	if(i==0) triggerPassed = true;
	triggerPassed &= iaccept ;
      }
      else triggerPassed |= iaccept ;
    }

    mySelections_->set("HLT passed", triggerPassed || goodTriggers==0);    
  } else{
    mySelections_->set("HLT passed", true );
  }



  typedef pat::ElectronCollection::const_iterator EI; 
  typedef map< double, EI, greater<double> > PtMap;
  typedef PtMap::const_iterator IM;

  PtMap sortedEles; 

  if(verbose_) cout<<"loop on electrons:"<<endl;

  for(EI ei = patElectrons->begin(); ei!=patElectrons->end(); ++ei)  {    

    if(verbose_) cout<<"\t ele: "<<ei->pt()<<" "<<ei->mva()<<endl;

    std::strbitset eleSel = eleSelector_.getBitTemplate();
    bool passed = eleSelector_( *ei, eleSel ); 

    if(!passed) continue;

    sortedEles.insert( make_pair(ei->pt(), ei) );
  }

  eleMultiplicity_ =  sortedEles.size();
  if(sortedEles.size() < 2 ) return false;
  else mySelections_->set("> 1 sorted electrons", true );

  const pat::Electron& leading    = *(sortedEles.begin()->second); 
  const pat::Electron& subleading = *((++(sortedEles.begin()))->second);

  leadEleHistos_->fill(leading, eventWeight);
  leadEleHistos_->fillWithVertex(leading, math::XYZPoint(0.,0.,0.),eventWeight);
  subLeadEleHistos_->fill(subleading, eventWeight);
  subLeadEleHistos_->fillWithVertex(subleading, math::XYZPoint(0.,0.,0.),eventWeight);

  selectedElectrons_.push_back(leading);
  selectedElectrons_.push_back(subleading);
 
  if(verbose_) cout<<"leading electron: pt/eta/phi "
		   <<leading.pt()<<", "
		   <<leading.eta()<<", "
		   <<leading.phi()
		   <<", subleading electron: pt/eta/phi "
		   <<subleading.pt()<<", "
		   <<subleading.eta()<<", "
		   <<subleading.phi()
		   <<endl;

  // two different selectors could be in principle applied to the
  // leading and subleading: for the moment we use the same for both 

  std::strbitset leadElePhaseSpaceSel = elePhaseSpaceSelector_.getBitTemplate();
  bool leadAccepted = elePhaseSpaceSelector_( leading, leadElePhaseSpaceSel ); 

  mySelections_->set("leading electron p_{T} cut", leadElePhaseSpaceSel[std::string("pt")] );
  mySelections_->set("leading electron #eta cut", leadElePhaseSpaceSel[std::string("eta")] );
 
  std::strbitset subLeadElePhaseSpaceSel = elePhaseSpaceSelector_.getBitTemplate();
  bool subLeadAccepted = elePhaseSpaceSelector_( subleading, subLeadElePhaseSpaceSel ); 

  mySelections_->set("subleading electron p_{T} cut", subLeadElePhaseSpaceSel[std::string("pt")] );
  mySelections_->set("subleading electron #eta cut", subLeadElePhaseSpaceSel[std::string("eta")] );

  mySelections_->set("2 accepted electrons", leadAccepted && subLeadAccepted );
   
  bool subSubLeadDisc = (sortedEles.size()>2 && (*(++(++(sortedEles.begin()))->second)).p4().Pt() < minPtSubSub_) 
    || sortedEles.size()==2 ;  
  mySelections_->set("subsubleading electron !accepted", subSubLeadDisc);

  mySelections_->set("opposite charge", leading.charge()*subleading.charge()==-1);

  /// isolation of the leading electron
  std::strbitset leadEleIsoSel = eleIsoSelector_.getBitTemplate();
  eleIsoSelector_( leading, leadEleIsoSel );
  mySelections_->set("leading electron loose isolated", leadEleIsoSel[std::string("combinedAbs")] );

  CombinedIsolation isoLead = eleIsoSelector_.combinedIsolation("combined");
  leadEleIsoHistos_->fill( isoLead.isolation(), isoLead.isolationRel(), eventWeight);

  /// isolation of the subleading electron
  std::strbitset subLeadEleIsoSel = eleIsoSelector_.getBitTemplate();
  eleIsoSelector_( subleading, subLeadEleIsoSel );
  mySelections_->set("subleading electron loose isolated", subLeadEleIsoSel[std::string("combinedAbs")] );
  
  CombinedIsolation isoSubLead = eleIsoSelector_.combinedIsolation("combined");
  subLeadEleIsoHistos_->fill( isoSubLead.isolation(), isoSubLead.isolationRel(), eventWeight);

  /// both eles isolated
  bool looseIsolated = leadEleIsoSel[std::string("combinedAbs")] && subLeadEleIsoSel[std::string("combinedAbs")];

  Candidate::LorentzVector Zp4 = leading.p4() + subleading.p4();

  ZMass_ = std::sqrt(Zp4.M2());
  mySelections_->set("m(ee) constrain", std::sqrt(Zp4.M2()) > minMass_ && std::sqrt(Zp4.M2()) < maxMass_ );

  ZEta_ = Zp4.Eta();
  ZPt_ = Zp4.Pt();

  if( leadAccepted && subLeadAccepted && subSubLeadDisc 
      && mySelections_->test("opposite charge") 
      && looseIsolated 
      //&& mySelections_->test("m(ee) constrain") 
      )
    {
      if( scanEvsByM_ ) std::cout<<iEvent.id()
				 << " => m(ee) " << std::sqrt(Zp4.M2()) 
				 <<std::endl;
      
      zeeHistos_->fill1DHistogram("hDiEleMass", std::sqrt(Zp4.M2()) );
      zeeHistos_->fill1DHistogram("hDiEleEta", Zp4.Eta()  );
      zeeHistos_->fill1DHistogram("hDiElePt",  Zp4.Pt()   );
    }


  return true;
}



void  PatZeeAnalyzer::addBranch(TTree *tree){
  tree->Branch("eleMultiplicity",  &eleMultiplicity_);
  tree->Branch("ZMass",&ZMass_);
  tree->Branch("ZEta",&ZEta_);
  tree->Branch("ZPt",&ZPt_);
  tree->Branch("SelectedElectrons",&selectedElectrons_);
}



void  PatZeeAnalyzer::addCutHistos(TList *aList){
  aList->Add(new TH1F("heleMultiplicity","Electrons multiplicity ; multiplicity ;Events",21,-0.5,20.5));
  aList->Add(new TH1F("hZMass","Di-Electrons mass ; mass ; Events",150,0,150));
  aList->Add(new TH1F("hZPt","Di-Electron pT ; p_{T} [GeV]; Events",150,0,150));
  aList->Add(new TH1F("hZEta","Di-Electron #eta; #eta ; Events",80,-4,4) );
}


void PatZeeAnalyzer::clear(){
  eleMultiplicity_ = 0;
  ZMass_ = 0;
  ZEta_ = 0;
  ZPt_ = 0;
  selectedElectrons_.clear();
}





#ifndef PatZeeAnalyzer_H
#define PatZeeAnalyzer_H

// -*- C++ -*-
//
//
// Original Author:  Artur Kalinowski
//         Created:  Wed Jul 22 14:16:44 CEST 2009
// $Id: PatZeeAnalyzer.h,v 1.7 2010/05/07 14:08:36 bianchi Exp $
//
//



#include "PFAnalyses/CommonTools/interface/FWLiteAnalyzer.h"

#include "PFAnalyses/CommonTools/interface/PatLeptonSelector.h"
#include "PFAnalyses/CommonTools/interface/CandidateSelector.h"
#include "PFAnalyses/CommonTools/interface/PatElectronSelector.h"
#include "PFAnalyses/CommonTools/interface/ElectronHistograms.h"
#include "PFAnalyses/CommonTools/interface/IsolationHistograms.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include <string>
#include <vector>

class PatZeeHistograms;
class ElectronHistograms;

class PatZeeAnalyzer : public FWLiteAnalyzer{

 public:

  ///Default contructor
  PatZeeAnalyzer(const std::string & aName);

  ///Default destructor
  virtual ~PatZeeAnalyzer();

  ///Method to initialise the object.
  ///It is called by the top level FWLiteTreeAnalyzer
  ///class.
  virtual void initialize(const edm::ParameterSet&, 
			  TFileDirectory&,
			  std::strbitset *aSelections);
  
  ///Method where the analysis is done.
  virtual bool analyze(const edm::EventBase& iEvent);

  ///Method which defined the branches to be added
  ///to the analysis TTree
  virtual void addBranch(TTree *tree);

  //virtual bool checkSelections(const edm::EventBase& iEvent, const std::string & type);

  ///Here we define histogram to be filled after each cut.
  ///Notice: only 1D here.
  ///Notice: name of the histogram the same as the name of the TBranch registered
  ///        in the relevant Analyzer, EXCEPT "h" at the beggining
  ///Notice: PF suffix in the name is necessary
  virtual void addCutHistos(TList *aList);

 private:

  // typedef
  typedef pat::ElectronCollection::const_iterator EI;
  //typedef pat::PFParticleCollection::const_iterator CI;
  typedef math::PtEtaPhiMLorentzVector LorentzVector;

  ///Method for registering the selections for this analyzer
  void registerCuts();

  ///Method reseting the values of class data members
  void clear();

  ///Class holding all the histograms to be filled
  PatZeeHistograms   *zeeHistos_;
  ElectronHistograms *leadEleHistos_;
  ElectronHistograms *subLeadEleHistos_;
  IsolationHistograms *leadEleIsoHistos_;
  IsolationHistograms *subLeadEleIsoHistos_;

  ///This name will be part of the output ROOT file
  std::string sampleName_;

  // electrons and offline vertexes collection
  edm::InputTag patElectronsLabel_;
  edm::InputTag vtxLabel_;
  edm::InputTag triggerResultsLabel_;
  
  // vector to store the selected electrons
  std::vector<pat::Electron> selectedElectrons_;

  // object selectors
  CandidateSelector elePhaseSpaceSelector_;
  PatElectronIsolationSelector eleIsoSelector_;
  PatElectronSelector eleSelector_;
 
  // variables to be put int the tree
  float  eleMultiplicity_; // keep it as a float!
  float ZMass_;
  float ZEta_;
  float ZPt_;

  // private members initialized from the cfg
  double minPtSubSub_;
  double minMass_;
  double maxMass_;
  bool verbose_;
  bool triggerVerbose_;
  bool scanEvsByM_;
  bool useTrigger_;
  bool triggerItemAND_;
  std::vector<std::string> triggerItemNames_;
  
};

#endif

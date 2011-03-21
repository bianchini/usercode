#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"

//#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif


void TMVAClassification(std::string ordering_ = "Pt"){

  TMVA::Tools::Instance();

  TString outfileName( "TMVA"+ordering_+"Ord.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification"+ordering_+"", outputFile, 
					      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D" );
  factory->AddVariable( "pt1", "pT-tag1", "GeV/c"         , 'F'  );
  factory->AddVariable( "pt2", "pT-tag2", "GeV/c"         , 'F'  );
  factory->AddVariable( "Deta","|y-tag1 - y-tag2|",""     , 'F'  );
  //factory->AddVariable( "opposite:=abs(eta1*eta2)/eta1/eta2","sign1*sign2",""             , 'F'  );
  //factory->AddVariable( "Dphi", "#Delta#phi" ,""             , 'F'  );
  factory->AddVariable( "Mjj", "M(tag1,tag2)", "GeV/c^{2}"  , 'F'  );

  factory->AddSpectator( "eta1",  "#eta_{tag1}" , 'F' );
  factory->AddSpectator( "eta2",  "#eta_{tag2}" , 'F' );

  TString fSignalName        = "/data_CMS/cms/lbianchini/VbfJetsStudy/nTupleVbf.root";
  TString fBackgroundName    = "/data_CMS/cms/lbianchini/VbfJetsStudy/nTupleZJets.root";

  TFile *fSignal(0); 
  TFile *fBackground(0); 

  fSignal = TFile::Open( fSignalName ); 
  fBackground = TFile::Open( fBackgroundName ); 

  if(!fSignal || !fBackground){
    std::cout << "ERROR: could not open files" << std::endl;
    exit(1);
  }

  TString tree = "outTree"+ordering_+"Ord";

  TTree *signal         = (TTree*)fSignal->Get(tree);
  TTree *background     = (TTree*)fBackground->Get(tree);

  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;

  factory->AddSignalTree    ( signal,     signalWeight     );
  factory->AddBackgroundTree( background, backgroundWeight );

  TCut mycuts = "pt1>0 && abs(eta1*eta2)/eta1/eta2<0"; 
  //TCut mycuts = ""; 
  TCut mycutb = "pt1>0 && abs(eta1*eta2)/eta1/eta2<0"; 
  //TCut mycutb = ""; 

  factory->PrepareTrainingAndTestTree( mycuts, mycutb,
				       "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
  
  factory->BookMethod( TMVA::Types::kCuts, "Cuts", 
		       "!H:!V:FitMethod=GA:EffSel:CutRangeMin[0]=20.:CutRangeMax[0]=999;CutRangeMin[1]=15.:CutRangeMax[1]=999.:CutRangeMin[2]=1.0:CutRangeMax[2]=9.:CutRangeMin[3]=100:CutRangeMax[3]=7000:VarProp=FSmart" );
  
  /*
  factory->BookMethod( TMVA::Types::kBDT, "BDT", 
		       "!H:!V:NTrees=200:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
  */

  factory->TrainAllMethods();
  
  factory->TestAllMethods();

  factory->EvaluateAllMethods();  

  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;      
  
  delete factory;

  //if (!gROOT->IsBatch()) TMVAGui( outfileName );

}

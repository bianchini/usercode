#include <iostream>

#include "PFAnalyses/Z/interface/PatZeeHistograms.h"



PatZeeHistograms::PatZeeHistograms(std::string fileName, int opt){

  AnalysisHistograms::init(fileName);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
PatZeeHistograms::PatZeeHistograms( TFileDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
PatZeeHistograms::PatZeeHistograms( TFileDirectory *myDir, const std::string & fileName ){

  AnalysisHistograms::init(myDir,fileName);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
PatZeeHistograms::~PatZeeHistograms(){ 

  std::cout<<"PatZeeHistograms::~PatZeeHistograms()"<<std::endl;

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void PatZeeHistograms::defineHistograms(){

  using namespace std;

 if(!histosInitialized_){

   add1DHistogram("hDiEleMass","Di-Electron mass; mass [GeV];Events",480,0,120,file_);
   add1DHistogram("hDiElePt","Di-Electron pT; p_{T} [GeV];Events",150,0,150,file_);
   add1DHistogram("hDiEleEta","Di-Electron #eta; #eta ;Events",80,4,-4,file_);

   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

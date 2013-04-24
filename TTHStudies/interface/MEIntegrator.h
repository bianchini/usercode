#ifndef MEINTEGRATOR_H
#define MEINTEGRATOR_H

#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TF2.h"
#include "TLegend.h"
#include "TList.h"
#include "THStack.h"
#include "TCut.h"
#include "TArrayF.h"
#include "TObjArray.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TKey.h"
#include "TMultiGraph.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"
#include "RooWorkspace.h"
#include "TH1F.h"

#include <string>
#include <map>

using namespace RooFit;
using namespace std;


class MEIntegrator {

 public:

  MEIntegrator( string , int );
  ~MEIntegrator();

  double Eval(const double* ) const;  
  void   setInputLV( TLorentzVector , TLorentzVector );
  void   SetPar(int);
  void   setJets( vector<TLorentzVector>* );
  void   setBtag( std::vector<float>* );
  double probability(const double*) const;
  void   saveJetParam( string );
  void   cachePdf( string , string , int );
  TH1F*  getCachedPdf( string );
  void   debug();


 private:
  
  RooWorkspace *w_;
  std::map<string, double> jetParam_; 
  std::map<string, TH1F*> variables_;
  vector<TLorentzVector> jets_;
  vector<float> bTagging_;
  TLorentzVector higgsLV_;
  TLorentzVector topLepLV_;
  int par_;

};


MEIntegrator::MEIntegrator( string fileName , int param  ) {

  cout << "Begin constructor" << endl;

  par_      = param;
  higgsLV_. SetPxPyPzE(0.,0.,0.,0.);
  topLepLV_.SetPxPyPzE(0.,0.,0.,0.);
  jets_.    clear();
  bTagging_.clear();


  TFile* file = TFile::Open(fileName.c_str(),"READ");
  w_ = (RooWorkspace*)file->Get("transferFuntions");

  // jets
  RooArgSet allVars = w_->allVars();
  TIterator* iter = allVars.createIterator();
  RooRealVar* var = 0;
  while( (var = (RooRealVar*)(*iter)() ) ){
    jetParam_[ string(var->GetName()) ] = var->getVal();
  }

  cachePdf( "pdfCsvLight",  "csvReco", 100);
  cachePdf( "pdfCsvHeavy",  "csvReco", 100);
  cachePdf( "pdfV1",             "X1", 100);
  cachePdf( "pdfV2",             "X2", 100);
  cachePdf( "pdfPtTTH",       "PtTTH", 100);
  cachePdf( "pdfPtHStar",   "PtHStar", 100);
  cachePdf( "pdfBetaW",       "BetaW", 100);
  cachePdf( "pdfGammaW",     "GammaW", 100);


  cout << "End constructor" << endl;
  //file->Close();

}


MEIntegrator::~MEIntegrator(){

  for(std::map<string, TH1F*>::iterator it = variables_.begin(); it!=variables_.end(); it++){
    if(it->second) 
      delete (it->second);
  }


}


double MEIntegrator::Eval(const double* x) const {
  double prob = probability(x);      
  if ( TMath::IsNaN(prob) ) prob = 0.;
  return prob;
}

void MEIntegrator::setInputLV( TLorentzVector lv1, TLorentzVector lv2){ 
  higgsLV_  = lv1;
  topLepLV_ = lv2;
}

void MEIntegrator::SetPar(int p){ 
  par_ = p; 
}
 

void MEIntegrator::setJets( std::vector<TLorentzVector>* jets){
  jets_.clear();
  for(unsigned int k = 0 ; k<jets->size() ; k++)
    jets_.push_back( (*jets)[k] );
}

void MEIntegrator::setBtag( std::vector<float>* bTagging ){
  bTagging_.clear();
  for(unsigned int k = 0 ; k<bTagging->size() ; k++)
    bTagging_.push_back( (*bTagging)[k] );
}


void MEIntegrator::saveJetParam( string paramName){

  RooRealVar* var = w_->var(paramName.c_str());
  jetParam_[ paramName ] = var->getVal();

}

void MEIntegrator::cachePdf( string pdfName, string varName, int nBins){
  
  RooAbsPdf* pdf = w_->pdf( pdfName.c_str() );
  if(pdf){
    RooRealVar* var = w_->var( varName.c_str() );

    TH1F* hCache = new TH1F("hCache","", nBins, var->getMin(), var->getMax() );

    for(int j = 1; j < hCache->GetNbinsX(); j++){
      float lvalue = ( TMath::Abs(var->getMax() - var->getMin() )/nBins)*(j-0.5) +  var->getMin();
      var->setVal( lvalue );
      hCache->SetBinContent(j, pdf->getVal( RooArgSet(*var) ) );
    }

    variables_[ pdfName ] = hCache;

  }

}


TH1F* MEIntegrator::getCachedPdf( string pdfName ){

  return variables_[ pdfName ];

}

 



double MEIntegrator::probability(const double* x) const{

  return 2.0;
  
}



void MEIntegrator::debug(){

  cout << "*** debug start *** " << endl;

  for(std::map<string, double>::iterator it = jetParam_.begin(); it!=jetParam_.end(); it++){
    cout << it->first << " => " << it->second << endl;
  }
  for(std::map<string, TH1F*>::iterator it = variables_.begin(); it!=variables_.end(); it++){
    if(it->second) 
      cout << it->first << " => " << (it->second)->Integral() << endl;
    else
      cout << it->first << " => null pointer" << endl;
  }

  cout << "*** debug end *** " << endl;

}


#endif

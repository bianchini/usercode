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
#include <limits>  

using namespace RooFit;
using namespace std;

#define MTOP 175.
#define MW    81.
#define Mb    5.


class MEIntegrator {

 public:

  MEIntegrator( string , int );
  ~MEIntegrator();

  double Eval(const double* ) const;  
  void   setInputLV( TLorentzVector , TLorentzVector );
  void   SetPar(int);
  void   setJets( vector<TLorentzVector>* );
  void   setBtag( std::vector<float>* );
  void   createMash();
  double probability(const double*) const;
  void   saveJetParam( string );
  void   cachePdf( string , string , int );
  TH1F*  getCachedPdf( string );
  TH2F*  getMash( );
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
  float pStar_;
  float EbStar_;
  float EWStar_;
  float EuStar_;
  TH2F* mash_;

};


MEIntegrator::MEIntegrator( string fileName , int param  ) {

  cout << "Begin constructor" << endl;

  par_      = param;
  higgsLV_. SetPxPyPzE(0.,0.,0.,0.);
  topLepLV_.SetPxPyPzE(0.,0.,0.,0.);
  jets_.    clear();
  bTagging_.clear();

  pStar_  = TMath::Sqrt( (MTOP*MTOP-(MW+Mb)*(MW+Mb) )*( MTOP*MTOP-(MW-Mb)*(MW-Mb) ) )/2./MTOP;
  EbStar_ =  (MTOP*MTOP - MW*MW + Mb*Mb)/2./MTOP;
  EWStar_ =  (MTOP*MTOP + MW*MW - Mb*Mb)/2./MTOP;
  EuStar_ =  MW/2;

  mash_ = new TH2F("mash","",100,-2.5,2.5, 100, -TMath::Pi(), TMath::Pi());


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

  delete mash_;

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


void   MEIntegrator::createMash(){

  mash_->Reset();

  for(unsigned int j = 0; j < jets_.size() ; j++){

    float eta = jets_[j].Eta();
    float phi = jets_[j].Phi();

    for(int binX = 1; binX<mash_->GetNbinsX(); binX++ ){
      float binCenterX = mash_->GetXaxis()->GetBinCenter(binX);
      for(int binY = 1; binY<mash_->GetNbinsY(); binY++ ){
	float binCenterY = mash_->GetYaxis()->GetBinCenter(binY);

	float dist = TMath::Sqrt( (eta-binCenterX)*(eta-binCenterX) + (phi-binCenterY)*(phi-binCenterY) );
	if( j<4 ){ // jets to veto
	  if( dist<0.50 || mash_->GetBinContent(binX,binY)==-1 ) 
	    mash_->SetBinContent(binX,binY,  -1);
	  else
	    mash_->SetBinContent(binX,binY,   0);
	}
	else{ // jets candidate W had
	  if( mash_->GetBinContent(binX,binY)!=-1 && (dist<0.50 || mash_->GetBinContent(binX,binY)==+1) ) 
	    mash_->SetBinContent(binX,binY,  +1);
	  else if( mash_->GetBinContent(binX,binY)!=-1 )
	    mash_->SetBinContent(binX,binY,   0);
	  else{}
	}

      }
    }

  }


}


TH2F*  MEIntegrator::getMash( ){
  return mash_;
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

  double prob = 1.0;
  
  double v1       = x[0];
  double v2       = x[1];
  double PtTTH    = x[2];
  double PhiTTH   = x[3];
  double BetaW    = x[4];
  double DeltaW   = x[5];
  double GammaW   = x[6];
  double EpsilonW = x[7];


  TLorentzVector higgsLV  = higgsLV_;
  TLorentzVector topLepLV = topLepLV_;
  TLorentzVector PTOT;
  PTOT.SetPxPyPzE( PtTTH*TMath::Cos( PhiTTH ), PtTTH*TMath::Sin( PhiTTH ), 8000*v2, 
		   TMath::Sqrt( PtTTH*PtTTH + 8000*v2*8000*v2  + 8000*8000*v1*v1 )  );

  TVector3 boostToCMS(PTOT.Px()/PTOT.E(), PTOT.Py()/PTOT.E(), PTOT.Pz()/PTOT.E() );

  TLorentzVector topHadLV = PTOT - higgsLV - topLepLV ;
  if( topHadLV.E()<MTOP ){
    prob =   std::numeric_limits<double>::min();
    return prob; 
  }

  topHadLV.SetE( TMath::Sqrt( (PTOT - higgsLV - topLepLV).P()*(PTOT - higgsLV - topLepLV).P() + MTOP*MTOP ) );
  //cout << "before boost... " << (topHadLV).M() << endl;

  higgsLV. Boost( -boostToCMS );
  topLepLV.Boost( -boostToCMS );

  double PtHStar   = higgsLV.P()/PTOT.M();
  double DphiHStar = TMath::Cos((higgsLV.Vect()).Angle(topLepLV.Vect()));


  //////////////////////// 
  TVector3 boostToTopHadCMS(topHadLV.Px()/topHadLV.E(), topHadLV.Py()/topHadLV.E(), topHadLV.Pz()/topHadLV.E() );
  TLorentzVector topHadLVInCMS = topHadLV;
  topHadLVInCMS.Boost( -boostToTopHadCMS );

  TVector3 versor1 = boostToTopHadCMS.Unit();
  TVector3 versor2 = (versor1.Orthogonal()).Unit();
  TVector3 versor3 = versor1.Cross(versor2);

  TVector3 dirB = BetaW*versor1 + TMath::Sqrt(1-BetaW*BetaW)*TMath::Cos( DeltaW )*versor2 + TMath::Sqrt(1-BetaW*BetaW)*TMath::Sin( DeltaW )*versor3;
  TLorentzVector bLV, wLV;
  bLV.SetPxPyPzE(  dirB.X()*pStar_,  dirB.Y()*pStar_,  dirB.Z()*pStar_ , EbStar_ );
  wLV.SetPxPyPzE( -dirB.X()*pStar_, -dirB.Y()*pStar_, -dirB.Z()*pStar_ , EWStar_ );
  TLorentzVector bInLab = bLV;
  bInLab.Boost( boostToTopHadCMS );


  
  TVector3 boostToWHadCMS(wLV.Px()/wLV.E(), wLV.Py()/wLV.E(), wLV.Pz()/wLV.E() );
  TLorentzVector WHadLVInCMS = wLV;
  WHadLVInCMS.Boost( -boostToWHadCMS );

  versor1 = boostToWHadCMS.Unit();
  versor2 = (versor1.Orthogonal()).Unit();
  versor3 = versor1.Cross(versor2);

  TVector3 dirW = GammaW*versor1 + TMath::Sqrt(1-GammaW*GammaW)*TMath::Cos( EpsilonW )*versor2 + TMath::Sqrt(1-GammaW*GammaW)*TMath::Sin( EpsilonW )*versor3;
  TLorentzVector w1LV, w2LV;
  w1LV.SetPxPyPzE(  dirW.X()*EuStar_,  dirW.Y()*EuStar_,  dirW.Z()*EuStar_ , EuStar_ );
  w2LV.SetPxPyPzE( -dirW.X()*EuStar_, -dirW.Y()*EuStar_, -dirW.Z()*EuStar_ , EuStar_ );
  TLorentzVector w1InLab = w1LV;
  TLorentzVector w2InLab = w2LV;
  w1InLab.Boost( boostToWHadCMS );
  w2InLab.Boost( boostToWHadCMS );
  w1InLab.Boost( boostToTopHadCMS );
  w2InLab.Boost( boostToTopHadCMS );
  //////////////////////// 





  ////////////////////////
  return prob;
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

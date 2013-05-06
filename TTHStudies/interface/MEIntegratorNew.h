#ifndef MEINTEGRATORNEW_H
#define MEINTEGRATORNEW_H

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

#define DEBUG 0

using namespace RooFit;
using namespace std;

#define Mtop_ 175.
#define MW    81.
#define Mb    5.


class MEIntegratorNew {

 public:

  MEIntegratorNew( string , int , int);
  ~MEIntegratorNew();

  double Eval(const double* ) const;  
  void   SetPar(int);
  void   setJets( vector<TLorentzVector>* );
  void   setBtag( std::vector<float>* );
  void   createMash();
  double        probability(const double*) const;
  unsigned int  findMatch(double, double) const;
  void   saveJetParam( string );
  void   cachePdf( string , string , int );
  void   cachePdf( string , string , string, int,    int );
  void   cachePdf( string , string , string, string, int, int, int );
  void   cachePdf( string , string , string, string, TArrayF, TArrayF, TArrayF);
  void   setMass(double);
  void   initTF() ;
  TH1*   getCachedPdf( string ) const;
  TH2F*  getMash( );
  TH1F*  getDebugHisto( );
  void   debug();

  void initializeVersors(int);

  void topHadEnergies    (double, double&, double&, double&, double&, int&);
  void topLepEnergies    (double, double,  double&, double&, double&, double&, int&);
  void topHadLostEnergies(double, double,  double,  double&, double&, double&, double&, int&);
  void higgsEnergies     (double, double&, double&, int&);
  double topHadJakobi    (double, double, double) const;
  double topLepJakobi    (double, double ) const;
  double higgsJakobi     (double, double ) const;
  double topHadDensity   (double, double)  const;
  double topLepDensity   (double, double)  const;
  double higgsDensity    (double)  const;
  double tthDensity      (double, double, double, double)  const;

 private:
  
  RooWorkspace *w_;
  std::map<string, double> jetParam_; 
  std::map<string, TH1F*> variables1D_;
  std::map<string, TH2F*> variables2D_;
  std::map<string, TH3F*> variables3D_;
  std::map<string, TF1* > transferFunctions_;

  // 0- lep 
  // 1- MET 
  // 2- b from top lep
  // 3- j1 from W
  // 4- j2 from W
  // 5- b from top hadr
  // 6- b1 from H
  // 7- b2 from H

  vector<TLorentzVector> jets_; 
  vector<float> bTagging_;
  TVector3 eLep_;
  TVector3 eBLep_;
  TVector3 eBHad_;
  TVector3 eW1Had_;
  TVector3 eW2Had_;
  TVector3 eB1_;
  TVector3 eB2_;

  int par_;
  int verbose_;
  float M_;
  float pStar_;
  float EbStar_;
  float EWStar_;
  float EuStar_;
  float dM2_;
  float dMh2_;
  float Mtop_;
  float Mb_;
  float Mw_;
  TH2F* mash_;

  TH1F* debugHisto1_;
  TH1F* pdfBetaWHad_, pdfBetaWHad_, pdfBetaWLep_, pdfGammaWLep_, pdfGammaTTH_;
  TH3F* pdf3D_;
  TF1* tfWjet1_,tfWjet2_, tfbHad_, tfbLep_, tfMetPt_,tfMetPhi_, tfHiggs1_, tfHiggs2_;

};


MEIntegratorNew::MEIntegratorNew( string fileName , int param , int verbose ) {

  cout << "Begin constructor" << endl;

  par_     = param;
  verbose_ = verbose;

  jets_.clear();
  bTagging_.clear();
  initializeVersors(0);

  Mtop_   = 174.3;
  Mb_     = 4.8;
  Mw_     = 80.19;
  pStar_  =  TMath::Sqrt( (Mtop_*Mtop_-(Mw_+Mb_)*(Mw_+Mb_) )*( Mtop_*Mtop_-(Mw_-Mb_)*(Mw_-Mb_) ) )/2./Mtop_;
  EbStar_ =  (Mtop_*Mtop_ - Mw_*Mw_ + Mb_*Mb_)/2./Mtop_;
  EWStar_ =  (Mtop_*Mtop_ + Mw_*Mw_ - Mb_*Mb_)/2./Mtop_;
  EuStar_ =  Mw_/2;
  dM2_    =  (Mtop_*Mtop_-Mb_*Mb_-Mw_*Mw_)*0.5;
  dMh2_   =  (Mtop_*Mtop_-M_*M_)*0.5;
  M_      = -99;

  mash_        = new TH2F("mash","",500,-2.5,2.5, 628, -TMath::Pi(), TMath::Pi());
  debugHisto1_ = new TH1F("debugHisto1","w1 pt", 100,0,400);

  TFile* file = TFile::Open(fileName.c_str(),"READ");
  w_ = (RooWorkspace*)file->Get("transferFuntions");

  // jets
  RooArgSet allVars = w_->allVars();
  TIterator* iter = allVars.createIterator();
  RooRealVar* var = 0;
  while( (var = (RooRealVar*)(*iter)() ) ){
    jetParam_[ string(var->GetName()) ] = var->getVal();
  }

  int nBinsX = 27;
  TArrayF binsX(nBinsX+1);
  cout << "Making histograms with " << nBinsX << " bins:" << endl;
  binsX[0] = 0.; binsX[1] = 0.06;  
  for(int k = 2; k < 21; k++)
    binsX[k] = 0.06 + (k-1)*0.005;
  binsX[21] = 0.17;  binsX[22] = 0.18;  binsX[23] = 0.19;  binsX[24] = 0.20;
  binsX[25] = 0.225; binsX[26] = 0.25;  binsX[27] = 0.30; 

  int nBinsMassTT = 23;
  TArrayF binsMassTT(nBinsMassTT+1);
  cout << "Making histograms with " << nBinsMassTT << " bins:" << endl;
  binsMassTT[0] = 350;
  for(int k = 1; k < 19; k++)
    binsMassTT[k] = 350 + 25*k;
  binsMassTT[19] = 900;  binsMassTT[20] = 1000; binsMassTT[21] = 1200; binsMassTT[22] = 1400;
  binsMassTT[23] = 2000;

  int nBinsGammaTT = 3;
  TArrayF binsGammaTT(nBinsGammaTT+1);
  cout << "Making histograms with " << nBinsGammaTT << " bins:" << endl;
  binsGammaTT[0] = -1.01;
  binsGammaTT[1] = -0.5;
  binsGammaTT[2] =  0.5;
  binsGammaTT[3] =  1.01;

  cachePdf( "pdfGammaWHad",     "GammaW",    100);
  cachePdf( "pdfBetaWHad",      "BetaW",     100);
  cachePdf( "pdfGammaWLep",     "GammaWLep", 100);
  cachePdf( "pdfBetaWLep",      "BetaWLep",  100);
  cachePdf( "pdfGammaTTH",      "GammaTTH",  100);
  cachePdf( "pdf3D",            "X1", "MassTT","GammaTT", binsX, binsMassTT, binsGammaTT);

  pdfGammaWHad_ = (TH1F*)this->getCachedPdf("pdfGammaWHad");
  pdfBetaWHad_  = (TH1F*)this->getCachedPdf("pdfBetaWHad");
  pdfBetaWLep_  = (TH1F*)this->getCachedPdf("pdfBetaWLep");
  pdfGammaWLep_ = (TH1F*)this->getCachedPdf("pdfGammaWLep");
  pdfGammaTTH_  = (TH1F*)this->getCachedPdf("pdfGammaTTH");
  pdf3D_        = (TH3F*)this->getCachedPdf("pdf3D");

  
  cout << "End constructor" << endl;

}


MEIntegratorNew::~MEIntegratorNew(){

  for(std::map<string, TH1F*>::iterator it = variables1D_.begin(); it!=variables1D_.end(); it++){
    if(it->second) 
      delete (it->second);
  }
  for(std::map<string, TH2F*>::iterator it = variables2D_.begin(); it!=variables2D_.end(); it++){
    if(it->second) 
      delete (it->second);
  }
  for(std::map<string, TH3F*>::iterator it = variables3D_.begin(); it!=variables3D_.end(); it++){
    if(it->second) 
      delete (it->second);
  }
  for(std::map<string, TF1*>::iterator it = transferFunctions_.begin(); it!=transferFunctions_.end(); it++){
    if(it->second) 
      delete (it->second);
  }

  delete mash_; delete debugHisto1_;

}

void MEIntegratorNew::initializeVersors(int withJetList){
  
  if(withJetList==0){
    eLep_   = TVector3(0.,0.,1.);
    eBLep_  = TVector3(0.,0.,1.);
    eBHad_  = TVector3(0.,0.,1.);
    eW1Had_ = TVector3(0.,0.,1.);
    eW2Had_ = TVector3(0.,0.,1.);
    eB1_    = TVector3(0.,0.,1.);
    eB2_    = TVector3(0.,0.,1.);
    return;
  }
  else if( withJetList==1 && jets_.size() < 7 ){
    cout << "Jets are not properly initliazied!!" << endl;
    return;
  }
  else if(  withJetList==1 && jets_.size() == 8 ){
    eLep_  = (jets_[0].Vect()).Unit();
    eBLep_ = (jets_[2].Vect()).Unit();
    eBHad_ = (jets_[5].Vect()).Unit();
    eW1Had_= (jets_[3].Vect()).Unit();
    eW2Had_= (jets_[4].Vect()).Unit();
    eB1_   = (jets_[6].Vect()).Unit();
    eB2_   = (jets_[7].Vect()).Unit();
   }
  else{cout << "Problems in MEIntegratorNew::initializeVersors" << endl;}
  
  return;
}

 
void initTF(){
  
  string bin = "Bin0";
  if(  TMath::Abs( eW1Had_.Eta() )<1.0 ){
    bin = "Bin0";
  }
  else{
    bin = "Bin1";
  }

  TF1* tfWjet1 = new TF1("tfWjet1",Form("(TMath::Erf(%f*x+%f)+%f) * TMath::Gaus( x, x , x*TMath::Sqrt( %f/x + %f/x/x) ) ", 
					(jetParam_.find("param0AccLight"+bin))->second, 
					(jetParam_.find("param1AccLight"+bin))->second,
					(jetParam_.find("param2AccLight"+bin))->second,
					(jetParam_.find("param0resolLight"+bin))->second*(jetParam_.find("param0resolLight"+bin))->second,
					(jetParam_.find("param1resolLight"+bin))->second*(jetParam_.find("param1resolLight"+bin))->second
					),1,1000);
  tfWjet1->SetNpx(1000);
  transferFunctions_["tfWjet1"] = tfWjet1;
  tfWjet1_ = tfWjet1;

  if(  TMath::Abs( eW2Had_.Eta() )<1.0 ){
    bin = "Bin0";
  }
  else{
    bin = "Bin1";
  }

  TF1* tfWjet2 = new TF1("tfWjet2",Form("(TMath::Erf(%f*x+%f)+%f) * TMath::Gaus( x, x , x*TMath::Sqrt( %f/x + %f/x/x) ) ", 
					(jetParam_.find("param0AccLight"+bin))->second, 
					(jetParam_.find("param1AccLight"+bin))->second,
					(jetParam_.find("param2AccLight"+bin))->second,
					(jetParam_.find("param0resolLight"+bin))->second*(jetParam_.find("param0resolLight"+bin))->second,
					(jetParam_.find("param1resolLight"+bin))->second*(jetParam_.find("param1resolLight"+bin))->second
					),1,1000);
  tfWjet2->SetNpx(1000);
  transferFunctions_["tfWjet2"] = tfWjet2;
  tfWjet2_ = tfWjet2;

  if(  TMath::Abs( eBHad_.Eta() )<1.0 ){
    bin = "Bin0";
  }
  else{
    bin = "Bin1";
  }

  TF1* tfbHad = new TF1("tfbHad",Form("(TMath::Erf(%f*x+%f)+%f) * TMath::Gaus( x, x , x*TMath::Sqrt( %f/x + %f/x/x) ) ", 
				      (jetParam_.find("param0AccHeavy"+bin))->second, 
				      (jetParam_.find("param1AccHeavy"+bin))->second,
				      (jetParam_.find("param2AccHeavy"+bin))->second,
				      (jetParam_.find("param0resolHeavy"+bin))->second*(jetParam_.find("param0resolHeavy"+bin))->second,
				      (jetParam_.find("param1resolHeavy"+bin))->second*(jetParam_.find("param1resolHeavy"+bin))->second
				      ),1,1000);
  tfbHad->SetNpx(1000);
  transferFunctions_["tfbHad"] = tfbHad;
  tfbHad_ = tfbHad;

  if(  TMath::Abs( eBLep_.Eta() )<1.0 ){
    bin = "Bin0";
  }
  else{
    bin = "Bin1";
  }

  TF1* tfbLep = new TF1("tfbLep",Form("(TMath::Erf(%f*x+%f)+%f) * TMath::Gaus( x, x , x*TMath::Sqrt( %f/x + %f/x/x) ) ", 
				      (jetParam_.find("param0AccHeavy"+bin))->second, 
				      (jetParam_.find("param1AccHeavy"+bin))->second,
				      (jetParam_.find("param2AccHeavy"+bin))->second,
				      (jetParam_.find("param0resolHeavy"+bin))->second*(jetParam_.find("param0resolHeavy"+bin))->second,
				      (jetParam_.find("param1resolHeavy"+bin))->second*(jetParam_.find("param1resolHeavy"+bin))->second
				      ),1,1000);
  tfbLep->SetNpx(1000);
  transferFunctions_["tfbLep"] = tfbLep;
  tfbLep_ = tfbLep;

  if(  TMath::Abs( eB1_.Eta() )<1.0 ){
    bin = "Bin0";
  }
  else{
    bin = "Bin1";
  }

  TF1* tfHiggs1 = new TF1("tfHiggs1",Form("(TMath::Erf(%f*x+%f)+%f) * TMath::Gaus( x, x , x*TMath::Sqrt( %f/x + %f/x/x) ) ", 
				      (jetParam_.find("param0AccHeavy"+bin))->second, 
				      (jetParam_.find("param1AccHeavy"+bin))->second,
				      (jetParam_.find("param2AccHeavy"+bin))->second,
				      (jetParam_.find("param0resolHeavy"+bin))->second*(jetParam_.find("param0resolHeavy"+bin))->second,
				      (jetParam_.find("param1resolHeavy"+bin))->second*(jetParam_.find("param1resolHeavy"+bin))->second
				      ),1,1000);
  tfHiggs1->SetNpx(1000);
  transferFunctions_["tfHiggs1"] = tfHiggs1;
  tfHiggs1_ = tfHiggs1;

  if(  TMath::Abs( eB2_.Eta() )<1.0 ){
    bin = "Bin0";
  }
  else{
    bin = "Bin1";
  }

  TF1* tfHiggs2 = new TF1("tfHiggs2",Form("(TMath::Erf(%f*x+%f)+%f) * TMath::Gaus( x, x , x*TMath::Sqrt( %f/x + %f/x/x) ) ", 
				      (jetParam_.find("param0AccHeavy"+bin))->second, 
				      (jetParam_.find("param1AccHeavy"+bin))->second,
				      (jetParam_.find("param2AccHeavy"+bin))->second,
				      (jetParam_.find("param0resolHeavy"+bin))->second*(jetParam_.find("param0resolHeavy"+bin))->second,
				      (jetParam_.find("param1resolHeavy"+bin))->second*(jetParam_.find("param1resolHeavy"+bin))->second
				      ),1,1000);
  tfHiggs2->SetNpx(1000);
  transferFunctions_["tfHiggs2"] = tfHiggs2;
  tfHiggs2_ = tfHiggs2;


  /////// do something here for the MET



}



void MEIntegratorNew::setMass(double mass){
  M_ = mass;
}

void MEIntegratorNew::topHadEnergies(double E1, double& E2, double& E3, double& cos1, double& cos2, int& errFlag ){
  
  TVector3 e1 = eW1Had_;
  TVector3 e2 = eW2Had_;
  TVector3 e3 = eBHad_;

  double a12 = e1.Angle(e2);
  double a13 = e1.Angle(e3);
  double a23 = e2.Angle(e3);
  
  E2 = Mw_*Mw_/E1/(4*TMath::Sin(a12/2.)*TMath::Sin(a12/2.));

  double a = E1+E2;
  double b = E1*TMath::Cos(a13)+E2*TMath::Cos(a23);

  if( dM2*dM2 - (a*a - b*b)*Mb_*Mb_ < 0){
    errFlag = 1;
    continue;
  }
  
  double E3_1 =  (a*dM2 + b*TMath::Sqrt(dM2*dM2 - (a*a - b*b)*Mb_*Mb_))/(a*a - b*b) ;
  double E3_2 =  (a*dM2 - b*TMath::Sqrt(dM2*dM2 - (a*a - b*b)*Mb_*Mb_))/(a*a - b*b) ;
  
  E3 = E3_1;

  if(E3<Mb_){
    errFlag = 1;
    return;
  }

  TLorentzVector w1( (e1.Unit())*E1,   E1);
  TLorentzVector w2( (e2.Unit())*E2,   E2);
  TLorentzVector blv((e3.Unit())*(TMath::Sqrt(E3*E3 - Mb_*Mb_)), E3);


  TVector3 boost = (w1+w2+blv).BoostVector();
  w1.Boost(  -boost );
  w2.Boost(  -boost );
  blv.Boost( -boost );

  cos1 = TMath::Cos( blv.Angle(boost) );
    
  TVector3 boostW = (w1+w2).BoostVector();
  w1.Boost(  -boostW );
  w2.Boost(  -boostW );
    
  cos2 = TMath::Cos( w1.Angle(boostW) );
  
  return;
}


void MEIntegratorNew::topLepEnergies(double nuPhi, double nuTheta, double& Enu, double& Eb, double& cos1, double& cos2, int& errFlag ){

  double Elep =  jets_[0].E();
  TVector3 e1 =  eLep_;
  TVector3 e2 =  eBLep_;
  TVector3 e3(0.,0.,1.); // neutrino
  e3.SetTheta( TMath::ACos( nuTheta ) ); 
  e3.SetPhi  ( nuPhi);   
  e3.SetMag  ( 1.);  

  double a12 = e1.Angle(e2); // lep - b
  double a13 = e1.Angle(e3); // lep - nu
  double a23 = e2.Angle(e3); // b   - nu
  
  Enu = Mw_*Mw_/ Elep / (4*TMath::Sin(a13/2.)*TMath::Sin(a13/2.));

  double a = Elep+Enu;
  double b = Elep*TMath::Cos(a12)+Enu*TMath::Cos(a23);

  if( dM2*dM2 - (a*a - b*b)*Mb_*Mb_ < 0){
    errFlag = 1;
    continue;
  }

  double E3_1 =  (a*dM2 + b*TMath::Sqrt(dM2*dM2 - (a*a - b*b)*Mb_*Mb_))/(a*a - b*b) ;
  double E3_2 =  (a*dM2 - b*TMath::Sqrt(dM2*dM2 - (a*a - b*b)*Mb_*Mb_))/(a*a - b*b) ;

  Eb = E3_1;

  if(E3<Mb_){
    errFlag = 1;
    return;
  }

  TLorentzVector wLep( (e1.Unit())*Elep,   Elep);
  TLorentzVector blv(  (e2.Unit())*(TMath::Sqrt(Eb*Eb - Mb_*Mb_)), Eb);
  TLorentzVector wNu ( (e3.Unit())*Enu,    Enu);


  TVector3 boost = (wLep+wNu+blv).BoostVector();
  wLep.Boost ( -boost );
  wNu.Boost  ( -boost );
  blv.Boost  ( -boost );

  cos1 = TMath::Cos( blv.Angle(boost) );
    
  TVector3 boostW = (wLep+wNu).BoostVector();
  wLep.Boost( -boostW );
  wNu.Boost ( -boostW );
    
  cos2 = TMath::Cos( wLep.Angle(boostW) );
  
  return;
}



void MEIntegratorNew::topHadLostEnergies(double missPhi, double missTheta, double missE, double& E1, double& Eb, double& cos1, double& cos2, int& errFlag ){

}


void MEIntegratorNew::higgsEnergies(double E1, double& E2, double& cos1, int& errFlag){

  if(E1<Mb_){
    errFlag = 1;
    return;
  }

  TVector3 e1 = eB1_;
  TVector3 e2 = eB2_;

  double a12 = e1.Angle(e2); // b1 - b2

  double a = E1;
  double b = TMath::Sqrt(E1*E1 - Mb_*Mb_)*TMath::Cos(a12);

  if( dM2*dM2 - (a*a - b*b)*Mb_*Mb_ < 0){
    errFlag = 1;
    return;
  }

  double E2_1 = (a*dM2 + b*TMath::Sqrt(dM2*dM2 - (a*a - b*b)*Mb_*Mb_))/(a*a - b*b);
  double E2_2 = (a*dM2 - b*TMath::Sqrt(dM2*dM2 - (a*a - b*b)*Mb_*Mb_))/(a*a - b*b);

  E2 = E2_1;

  if(E2<Mb_){
    errFlag = 1;
    return;
  }
  
  TLorentzVector b1( (e1.Unit())*(TMath::Sqrt(E1*E1 - Mb_*Mb_)),   E1);
  TLorentzVector b2( (e2.Unit())*(TMath::Sqrt(E2*E2 - Mb_*Mb_)),   E2);


  TVector3 boost = (b1+b2).BoostVector();
  b1.Boost( -boost );
  b2.Boost( -boost );

  cos1 = TMath::Cos( b1.Angle(boost) );

  return;
}


double MEIntegratorNew::Eval(const double* x) const {
  double prob = probability(x);      
  if ( TMath::IsNaN(prob) ) prob = 0.;
  return prob;
}


void MEIntegratorNew::SetPar(int p){ 
  par_ = p; 
}
 

void MEIntegratorNew::setJets( std::vector<TLorentzVector>* jets){
  jets_.clear();
  for(unsigned int k = 0 ; k<jets->size() ; k++)
    jets_.push_back( (*jets)[k] );
}

void MEIntegratorNew::setBtag( std::vector<float>* bTagging ){
  bTagging_.clear();
  for(unsigned int k = 0 ; k<bTagging->size() ; k++)
    bTagging_.push_back( (*bTagging)[k] );
}


void   MEIntegratorNew::createMash(){

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
	  if( mash_->GetBinContent(binX,binY)==0 && dist<0.10 ) 
	    mash_->SetBinContent(binX,binY,   j);
	}

      }
    }

  }


}


TH2F*  MEIntegratorNew::getMash( ){
  return mash_;
}

TH1F* MEIntegratorNew::getDebugHisto(){
  return debugHisto1_;
}

void MEIntegratorNew::saveJetParam( string paramName){

  RooRealVar* var = w_->var(paramName.c_str());
  jetParam_[ paramName ] = var->getVal();

}

void MEIntegratorNew::cachePdf( string pdfName, string varName, int nBins){
  
  RooAbsPdf* pdf = w_->pdf( pdfName.c_str() );
  if(pdf){
    RooRealVar* var = w_->var( varName.c_str() );
    
    TH1F* hCache = new TH1F("hCache","", nBins, var->getMin(), var->getMax() );
    
    for(int j = 1; j < hCache->GetNbinsX(); j++){
      float lvalue = hCache->GetXaxis()->GetBinCenter(j);
      var->setVal( lvalue );
      hCache->SetBinContent(j, pdf->getVal( RooArgSet(*var) ) );
    }
    
    variables1D_[ pdfName ] = hCache;
    
  }
  
}


void MEIntegratorNew::cachePdf( string pdfName, string varName1, string varName2, int nBins1, int nBins2){
  
  RooAbsPdf* pdf = w_->pdf( pdfName.c_str() );
  if(pdf){
    RooRealVar* var1 = w_->var( varName1.c_str() );
    RooRealVar* var2 = w_->var( varName2.c_str() );
    
    TH2F* hCache = new TH2F("hCache","", nBins1, var1->getMin(), var1->getMax(), nBins2, var2->getMin(), var2->getMax() );
    
    for(int j = 1; j < hCache->GetNbinsX(); j++){
      float lvalue1 = hCache->GetXaxis()->GetBinCenter(j);
      var1->setVal( lvalue1 );
      for(int k = 1; k < hCache->GetNbinsY(); k++){
	float lvalue2 = hCache->GetYaxis()->GetBinCenter(k);
	var2->setVal( lvalue2 );
	hCache->SetBinContent(j,k, pdf->getVal( RooArgSet(*var1,*var2) ) );
      }
    }
    
    variables2D_[ pdfName ] = hCache;
    
  }
  
}


void MEIntegratorNew::cachePdf( string pdfName, string varName1, string varName2, string varName3, int nBins1, int nBins2, int nBins3){
  
  RooAbsPdf* pdf = w_->pdf( pdfName.c_str() );
  if(pdf){
    RooRealVar* var1 = w_->var( varName1.c_str() );
    RooRealVar* var2 = w_->var( varName2.c_str() );
    RooRealVar* var3 = w_->var( varName3.c_str() );

    TH3F* hCache = new TH3F("hCache","", nBins1, var1->getMin(), var1->getMax(), nBins2, var2->getMin(), var2->getMax(), nBins3, var3->getMin(), var3->getMax() );

    for(int j = 1; j < hCache->GetNbinsX(); j++){
      float lvalue1 = hCache->GetXaxis()->GetBinCenter(j);
      var1->setVal( lvalue1 );
      for(int k = 1; k < hCache->GetNbinsY(); k++){
	float lvalue2 =  hCache->GetYaxis()->GetBinCenter(k);
	var2->setVal( lvalue2 );
	for(int m = 1; m < hCache->GetNbinsZ(); m++){
	  float lvalue3 =  hCache->GetZaxis()->GetBinCenter(m);
	  var3->setVal( lvalue3 );
	  hCache->SetBinContent(j,k,m, pdf->getVal( RooArgSet(*var1,*var2,*var3) ) );
	}
      }
    }
    
    variables3D_[ pdfName ] = hCache;

  }

}



void MEIntegratorNew::cachePdf( string pdfName, string varName1, string varName2, string varName3 , TArrayF bins1, TArrayF bins2, TArrayF bins3){
  
  RooAbsPdf* pdf = w_->pdf( pdfName.c_str() );
  if(pdf){
    RooRealVar* var1 = w_->var( varName1.c_str() );
    RooRealVar* var2 = w_->var( varName2.c_str() );
    RooRealVar* var3 = w_->var( varName3.c_str() );

    TH3F* hCache = new TH3F("hCache","", 
			    bins1.GetSize()-1, bins1.GetArray() , 
			    bins2.GetSize()-1, bins2.GetArray(), 
			    bins3.GetSize()-1, bins3.GetArray());

    for(int j = 1; j < hCache->GetNbinsX(); j++){
      float lvalue1 = hCache->GetXaxis()->GetBinCenter(j);
      var1->setVal( lvalue1 );
      for(int k = 1; k < hCache->GetNbinsY(); k++){
	float lvalue2 =  hCache->GetYaxis()->GetBinCenter(k);
	var2->setVal( lvalue2 );
	for(int m = 1; m < hCache->GetNbinsZ(); m++){
	  float lvalue3 = hCache->GetZaxis()->GetBinCenter(m);;
	  var3->setVal( lvalue3 );
	  hCache->SetBinContent(j,k,m, pdf->getVal( RooArgSet(*var1,*var2,*var3) ) );
	}
      }
    }
    
    variables3D_[ pdfName ] = hCache;

  }

}




TH1* MEIntegratorNew::getCachedPdf( string pdfName ) const{

  if( variables1D_.find(pdfName) !=  variables1D_.end() ) 
    return (variables1D_.find(pdfName))->second;
  else if( variables2D_.find(pdfName) !=  variables2D_.end() ) 
    return (variables2D_.find(pdfName))->second;
  else if( variables3D_.find(pdfName) !=  variables3D_.end() ) 
    return (variables3D_.find(pdfName))->second;

  else return 0;
}

 
double topHadJakobi( double Eq1, double Eq2, double EbHad) const {
  return EbHad * Eq2 * Eq1 ;
}

double topLepJakobi( double Enu, double EbLep) const {
  return EbLep * Enu ;
}

double higgsJakobi ( double Eh1, double Eh2) const {
  return Eh1 * Eh2 ;
}


double topHadDensity ( double cos1, double cos2) const {

  if( !pdfBetaWHad_ || !pdfGammaWHad_ ){
    if(verbose_) cout << "topHadDensity tries to access non-existing histogram" << endl;
    return 0.0;
  }


  double val1 = pdfBetaWHad_ ->Interpolate(cos1);
  double val2 = pdfGammaWHad_->Interpolate(cos2);

  return val1*val2;

}

double topLepDensity ( double cos1, double cos2) const {

  if( !pdfBetaWLep_ || !pdfGammaWLep_){
    if(verbose_) cout << "topLepDensity tries to access non-existing histogram" << endl;
    return 0.0;
  }

  double val1 = pdfBetaWLep_ ->Interpolate(cos1);
  double val2 = pdfGammaWLep_->Interpolate(cos2);

  return val1*val2;

}

double higgsDensity ( double cos1 ) const {
  return 1.0;
}

double tthDensity ( double Q, double m12, double cos1Star, double cos3) const {

  if( !pdfGammaTTH || !pdf3D){
    if(verbose_) cout << "tthDensity tries to access non-existing histogram" << endl;
    return 0.0;
  }

  double p1Star = (m12*m12 - 4*Mtop_*Mtop_) >=0 ? 
    TMath::Sqrt((m12*m12 - 4*Mtop_*Mtop_)*(m12*m12))/2./m12 : -999.;
  double p3     = (Q*Q - (m12+M_)*(m12+M_)) >=0 && (Q*Q - (m12 - M_)*(m12 - M_)) >=0 ? 
    TMath::Sqrt((Q*Q - (m12+M_)*(m12+M_) )*( Q*Q - (m12 - M_)* (m12 - M_)))/2./Q : -999;

  if( TMath::IsNaN(p1Star) || TMath::IsNaN(p3) || p1Star<=0 || p3 <=0){
    if(verbose_) cout << "tthDensity evaluated for non-acceptable values of p1Star and/or p3" << endl;
    return 0.0;
  }

  double val1 = pdf3D_->Interpolate(Q, m12, cos1Star);
  double val2 = pdfGammaTTH_->Interpolate( cos3 );

  return (val1*val2)/p1Star/p3;

}





double MEIntegratorNew::probability(const double* x) const{

  double prob = 1.0;
  
  double Eq1         = x[0];
  double cosThetaNu  = x[1];
  double phiNu       = x[2];
  double Eh1         = x[3];

  TVector3 eNu(0.,0.,1.); // neutrino
  eNu.SetTheta( TMath::ACos( cosThetaNu ) ); 
  eNu.SetPhi  ( phiNu );   
  eNu.SetMag  ( 1.);  

  if(verbose){
    cout << "Eq = "         << Eq << endl; 
    cout << "cosThetaNu = " << cosThetaNu << endl; 
    cout << "phiNu = "      << phiNu  << endl; 
    cout << "Eh = "         << Eh  << endl; 
  }

  int erFlag = 0;

  double Eq2, EbHad, cos1Had, cos2Had;
  topHadEnergies( Eq1, Eq2, EbHad, cos1Had, cos2Had, errFlag );

  if(errFlag){
    if(verbose_) cout << "Problems with topHadEnergies" << endl;
    return 0.0;
  }

  double Enu, EbLep, cos1Lep, cos2Lep;
  topLepEnergies( phiNu, cosThetaNu, Enu, EbLep, cos1Lep, cos2Lep, erFlag  );

  if(errFlag){
    if(verbose_) cout << "Problems with topLepEnergies" << endl;
    return 0.0;
  }

  double Eh2, cos1Higgs;
  higgsEnergies( Eh1, Eh2, erFlag );
  
  if(errFlag){
    if(verbose_) cout << "Problems with higgsEnergies" << endl;
    return 0.0;
  }

  TLorentzVector W1Had( eW1Had_*Eq1, Eq1 );
  TLorentzVector W2Had( eW2Had_*Eq2, Eq2 );
  TLorentzVector WHad = W1Had + W2Had;
  TLorentzVector bHad( eBHad_*TMath::Sqrt(EbHad*EbHad - Mb_*Mb_), EbHad );
  TLorentzVector topHad = WHad + bHad;

  TLorentzVector WLepLep = jets_[0];
  TLorentzVector WLepNu( eNu*Enu, Enu);
  TLorentzVector WLep = WLepLep + WLepNu;
  TLorentzVector bLep( eBLep_*TMath::Sqrt(EbLep*EbLep - Mb_*Mb_), EbLep );
  TLorentzVector topLep = WLep + bLep;

  TLorentzVector higgs1(  eB1_*TMath::Sqrt(Eh1*Eh1 -  Mb_*Mb_), Eh1);
  TLorentzVector higgs2(  eB2_*TMath::Sqrt(Eh2*Eh2 -  Mb_*Mb_), Eh2);
  TLorentzVector higgs = higgs1 + higgs2;

  TLorentzVector tot = topHad+topLep+higgs;

  TVector3 boostToCMS  = tot.BoostVector();

  topHad.Boost( -boostToCMS );
  topLep.Boost( -boostToCMS );
  higgs.Boost(  -boostToCMS );

  double cos3 = TMath::Cos((higgs.Vect()).Angle( boostToCMS ));

  TVector3 boostToTTCMS = (topLep+topHad).BoostVector();
  topLep.Boost( -boostToTTCMS );
  topHad.Boost( -boostToTTCMS );
  
  double cos1Star = TMath::Cos( (topLep.Vect()).Angle( boostToTTCMS ) );
  double m12      = (topHad+topLep).M();
  double Q        = tot.M();

  double MEpart = 
    topHadDensity(cos1Had,cos2Had) * topLepDensity(cos1Lep,cos2Lep) * higgsDensity(cos1Higgs) *
    topHadJakobi( Eq1, Eq2, EbHad) * topLepJakobi(Enu, EbLep) * higgsJakobi( Eh1, Eh2 ) *
    tthDensity( Q , m12, cos1Star, cos3);

  double TFpart = 
    tfWjet1_->Eval( W1Had.Pt() ) * tfWjet2_->Eval( W2Had.Pt() ) * tfbHad_->Eval( bHad.Pt() ) *
    tfbLep_->Eval( bLep.Pt() ) * /* something for MET  */ 
    tfHiggs1_->Eval( higgs1.Pt() ) * tfHiggs12_->Eval( higgs2.Pt() ) ;

  
  prob *= MEpart;
  prob *= TFpart;


  ////////////////////////

  return prob;
}


unsigned int MEIntegratorNew::findMatch(double eta, double phi) const{

  unsigned int match = 99;
  int bin   = mash_->FindBin(eta, phi);
  float res = mash_->GetBinContent(bin);
  if( res>-0.5 && res<0.5 )
    match = 99; // no match
  else if(  res<-0.5 )
    match = 98; // veto
  else
    match = res;
 

  return match;
}




void MEIntegratorNew::debug(){

  cout << "*** debug start *** " << endl;

  for(std::map<string, double>::iterator it = jetParam_.begin(); it!=jetParam_.end(); it++){
    cout << it->first << " => " << it->second << endl;
  }
  for(std::map<string, TH1F*>::iterator it = variables1D_.begin(); it!=variables1D_.end(); it++){
    if(it->second) 
      cout << it->first << " => " << (it->second)->Integral() << endl;
    else
      cout << it->first << " => null pointer" << endl;
  }

  cout << "*** debug end *** " << endl;

}


#endif

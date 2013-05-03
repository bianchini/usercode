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

#define MTOP 175.
#define MW    81.
#define Mb    5.


class MEIntegratorNew {

 public:

  MEIntegratorNew( string , int );
  ~MEIntegratorNew();

  double Eval(const double* ) const;  
  void   setInputLV( TLorentzVector , TLorentzVector );
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
  void   setMass(double);
  TH1F*  getCachedPdf( string ) const;
  TH2F*  getMash( );
  TH1F*  getDebugHisto( );
  void   debug();

  void initializeVersors(int);

  void topHadEnergies    (double, double&, double&, double&, double&, int&);
  void topLepEnergies    (double, double,  double&, double&, double&, double&, int&);
  void topHadLostEnergies(double, double,  double,  double&, double&, double&, double&, int&);
  void higgsEnergies     (double, double&, double&, int& errFlag)


 private:
  
  RooWorkspace *w_;
  std::map<string, double> jetParam_; 
  std::map<string, TH1F*> variables1D_;
  std::map<string, TH2F*> variables2D_;
  std::map<string, TH3F*> variables3D_;

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
  float M_;
  float pStar_;
  float EbStar_;
  float EWStar_;
  float EuStar_;
  float dM2_;
  float dMh2_;
  TH2F* mash_;

  TH1F* debugHisto1_;

};


MEIntegratorNew::MEIntegratorNew( string fileName , int param  ) {

  cout << "Begin constructor" << endl;

  par_      = param;
  jets_.    clear();
  bTagging_.clear();
  initializeVersors(0);

  pStar_  =  TMath::Sqrt( (MTOP*MTOP-(MW+Mb)*(MW+Mb) )*( MTOP*MTOP-(MW-Mb)*(MW-Mb) ) )/2./MTOP;
  EbStar_ =  (MTOP*MTOP - MW*MW + Mb*Mb)/2./MTOP;
  EWStar_ =  (MTOP*MTOP + MW*MW - Mb*Mb)/2./MTOP;
  EuStar_ =  MW/2;
  dM2_    =  (MTOP*MTOP-Mb*Mb-MW*MW)*0.5;
  dMh2_   =  (MTOP*MTOP-M_*M_)*0.5;
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


  cachePdf( "pdfGammaWHad",     "GammaW",    100);
  cachePdf( "histBetaWHad",     "BetaW",     100);
  cachePdf( "pdfGammaWLep",     "GammaWLep", 100);
  cachePdf( "pdfBetaWLep",      "BetaWLep",  100);
  cachePdf( "pdfGammaTTH",      "GammaTTH",  100);
  cachePdf( "pdf3D",            "GammaTT","MassTT","X1", 50, 50, 50);

  cout << "End constructor" << endl;

}


MEIntegratorNew::~MEIntegratorNew(){

  for(std::map<string, TH1F*>::iterator it = variables1D_.begin(); it!=variables1D_.end(); it++){
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
  
  E2 = MW*MW/E1/(4*TMath::Sin(a12/2.)*TMath::Sin(a12/2.));

  double a = E1+E2;
  double b = E1*TMath::Cos(a13)+E2*TMath::Cos(a23);

  if( dM2*dM2 - (a*a - b*b)*Mb*Mb < 0){
    errFlag = 1;
    continue;
  }
  
  double E3_1 =  (a*dM2 + b*TMath::Sqrt(dM2*dM2 - (a*a - b*b)*Mb*Mb))/(a*a - b*b) ;
  double E3_2 =  (a*dM2 - b*TMath::Sqrt(dM2*dM2 - (a*a - b*b)*Mb*Mb))/(a*a - b*b) ;
  
  E3 = E3_1;

  if(E3<Mb){
    errFlag = 1;
    return;
  }

  TLorentzVector w1( (e1.Unit())*E1,   E1);
  TLorentzVector w2( (e2.Unit())*E2,   E2);
  TLorentzVector blv((e3.Unit())*(TMath::Sqrt(E3*E3 - Mb*Mb)), E3);


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

  double Elep =  jets_[0].P();
  TVector3 e1 =  eLep_;
  TVector3 e2 =  eBLep_;
  TVector3 e3(0.,0.,1.); // neutrino
  e3.SetTheta( TMath::ACos( nuTheta ) ); 
  e3.SetPhi  ( nuPhi);   
  e3.SetMag  ( 1.);  

  double a12 = e1.Angle(e2); // lep - b
  double a13 = e1.Angle(e3); // lep - nu
  double a23 = e2.Angle(e3); // b   - nu
  
  Enu = MW*MW/ Elep / (4*TMath::Sin(a13/2.)*TMath::Sin(a13/2.));

  double a = Elep+Enu;
  double b = Elep*TMath::Cos(a12)+Enu*TMath::Cos(a23);

  if( dM2*dM2 - (a*a - b*b)*Mb*Mb < 0){
    errFlag = 1;
    continue;
  }

  double E3_1 =  (a*dM2 + b*TMath::Sqrt(dM2*dM2 - (a*a - b*b)*Mb*Mb))/(a*a - b*b) ;
  double E3_2 =  (a*dM2 - b*TMath::Sqrt(dM2*dM2 - (a*a - b*b)*Mb*Mb))/(a*a - b*b) ;

  Eb = E3_1;

  if(E3<Mb){
    errFlag = 1;
    return;
  }

  TLorentzVector wLep( (e1.Unit())*Elep,   Elep);
  TLorentzVector blv(  (e2.Unit())*(TMath::Sqrt(Eb*Eb - Mb*Mb)), Eb);
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

  if(E1<Mb){
    errFlag = 1;
    return;
  }

  TVector3 e1 = eB1_;
  TVector3 e2 = eB2_;

  double a12 = e1.Angle(e2); // b1 - b2

  double a = E1;
  double b = TMath::Sqrt(E1*E1 - Mb*Mb)*TMath::Cos(a12);

  if( dM2*dM2 - (a*a - b*b)*Mb*Mb < 0){
    errFlag = 1;
    return;
  }

  double E2_1 = (a*dM2 + b*TMath::Sqrt(dM2*dM2 - (a*a - b*b)*Mb*Mb))/(a*a - b*b);
  double E2_2 = (a*dM2 - b*TMath::Sqrt(dM2*dM2 - (a*a - b*b)*Mb*Mb))/(a*a - b*b);

  E2 = E2_1;

  if(E2<Mb){
    errFlag = 1;
    return;
  }
  
  TLorentzVector b1( (e1.Unit())*(TMath::Sqrt(E1*E1 - Mb*Mb)),   E1);
  TLorentzVector b2( (e2.Unit())*(TMath::Sqrt(E2*E2 - Mb*Mb)),   E2);


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

void MEIntegratorNew::setInputLV( TLorentzVector lv1, TLorentzVector lv2){ 
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
      float lvalue = ( TMath::Abs(var->getMax() - var->getMin() )/nBins)*(j-0.5) +  var->getMin();
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
      float lvalue1 = ( TMath::Abs(var1->getMax() - var1->getMin() )/nBins1)*(j-0.5) +  var1->getMin();
      var1->setVal( lvalue1 );
      for(int k = 1; k < hCache->GetNbinsY(); k++){
	float lvalue2 = ( TMath::Abs(var2->getMax() - var2->getMin() )/nBins2)*(k-0.5) +  var2->getMin();
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
      float lvalue1 = ( TMath::Abs(var1->getMax() - var1->getMin() )/nBins1)*(j-0.5) +  var1->getMin();
      var1->setVal( lvalue1 );
      for(int k = 1; k < hCache->GetNbinsY(); k++){
	float lvalue2 = ( TMath::Abs(var2->getMax() - var2->getMin() )/nBins2)*(k-0.5) +  var2->getMin();
	var2->setVal( lvalue2 );
	for(int m = 1; m < hCache->GetNbinsZ(); m++){
	  float lvalue3 = ( TMath::Abs(var3->getMax() - var3->getMin() )/nBins3)*(m-0.5) +  var3->getMin();
	  var3->setVal( lvalue3 );
	  hCache->SetBinContent(j,k,m, pdf->getVal( RooArgSet(*var1,*var2,*var3) ) );
	}
      }
    }
    
    variables3D_[ pdfName ] = hCache;

  }

}


TH1F* MEIntegratorNew::getCachedPdf( string pdfName ) const{

  return (variables1D_.find(pdfName))->second;

}

 



double MEIntegratorNew::probability(const double* x) const{

  double prob = 1.0;
  
  double v1       = x[0];
  double v2       = x[1];
  double PtTTH    = x[2];
  double PhiTTH   = x[3];
  double BetaW    = x[4];
  double DeltaW   = x[5];
  double GammaW   = x[6];
  double EpsilonW = x[7];

  if(DEBUG){
    cout << "v1    = " << v1 << endl; 
    cout << "v2    = " << v2 << endl; 
    cout << "PtTTH = " << PtTTH << endl; 
    cout << "PhiTTH = " <<PhiTTH  << endl; 
    cout << "BetaW = " << BetaW  << endl; 
    cout << "DeltaW = " << DeltaW  << endl; 
    cout << "GammaW = " <<  GammaW << endl; 
    cout << "EpsilonW = " <<  EpsilonW << endl; 
  }








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
  if(DEBUG){
    cout << "Top had in lab " << topHadLV.Pt() << "," << topHadLV.Eta() << ", " << topHadLV.Phi() << endl;
  }

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
  //TLorentzVector bLV, wLV;
  //bLV.SetPxPyPzE(  dirB.X()*pStar_,  dirB.Y()*pStar_,  dirB.Z()*pStar_ , EbStar_ );
  //wLV.SetPxPyPzE( -dirB.X()*pStar_, -dirB.Y()*pStar_, -dirB.Z()*pStar_ , EWStar_ );
  TLorentzVector bLV( (dirB*pStar_),    EbStar_ );
  TLorentzVector wLV( (dirB*(-pStar_)), EWStar_ );

  TLorentzVector bInLab = bLV;
  bInLab.Boost( boostToTopHadCMS );

  if(DEBUG){
    cout << "B had in CMS " << bLV.Pt() << "," << bLV.Eta() << ", " << bLV.Phi() << endl;
    cout << "B had in Lab " << bInLab.Pt() << "," << bInLab.Eta() << ", " << bInLab.Phi() << endl;
  }

  
  TVector3 boostToWHadCMS(wLV.Px()/wLV.E(), wLV.Py()/wLV.E(), wLV.Pz()/wLV.E() );
  TLorentzVector WHadLVInCMS = wLV;
  WHadLVInCMS.Boost( -boostToWHadCMS );

  versor1 = boostToWHadCMS.Unit();
  versor2 = (versor1.Orthogonal()).Unit();
  versor3 = versor1.Cross(versor2);

  TVector3 dirW = GammaW*versor1 + TMath::Sqrt(1-GammaW*GammaW)*TMath::Cos( EpsilonW )*versor2 + TMath::Sqrt(1-GammaW*GammaW)*TMath::Sin( EpsilonW )*versor3;
  //TLorentzVector w1LV, w2LV;
  //w1LV.SetPxPyPzE(  dirW.X()*EuStar_,  dirW.Y()*EuStar_,  dirW.Z()*EuStar_ , EuStar_ );
  //w2LV.SetPxPyPzE( -dirW.X()*EuStar_, -dirW.Y()*EuStar_, -dirW.Z()*EuStar_ , EuStar_ );
  TLorentzVector w1LV( (dirW*EuStar_),    EuStar_);
  TLorentzVector w2LV( (dirW*(-EuStar_)), EuStar_);

  TLorentzVector w1InLab = w1LV;
  TLorentzVector w2InLab = w2LV;
  w1InLab.Boost( boostToWHadCMS );
  w2InLab.Boost( boostToWHadCMS );
  w1InLab.Boost( boostToTopHadCMS );
  w2InLab.Boost( boostToTopHadCMS );
  //////////////////////// 

  if(DEBUG){
    //debugHisto1_->Fill( w1InLab.Pt() );
    cout << "W1 in Lab " << w1InLab.Pt() << ", " << w1InLab.Eta() << ", " << w1InLab.Phi() << endl;
    cout << "W2 in Lab " << w2InLab.Pt() << ", " << w2InLab.Eta() << ", " << w2InLab.Phi() << endl;

    cout << "Check Pt: " <<   "Top had in lab " << topHadLV.Pt() << " from components: " << (w1InLab+w2InLab+bInLab).Pt()   << endl;
    cout << "Check Eta: " <<  "Top had in lab " << topHadLV.Eta() << " from components: " << (w1InLab+w2InLab+bInLab).Eta() << endl;
    cout << "Check BetaW..." << endl;
    TLorentzVector bHelp = bInLab;
    bHelp.Boost( -boostToTopHadCMS );
    cout << BetaW << " <=> " << TMath::Cos( (bHelp.Vect()).Angle(boostToTopHadCMS) ) << endl;
    cout << "Check GammaW..." << endl;
    TLorentzVector w1Help = w1InLab;
    TLorentzVector wHelp  = w1InLab+w2InLab;
    wHelp.Boost( -boostToTopHadCMS );
    w1Help.Boost( -boostToTopHadCMS );
    w1Help.Boost( -wHelp.BoostVector() );
    cout << GammaW << " <=> " << TMath::Cos( (w1Help.Vect()).Angle(wHelp.BoostVector()) ) << endl;
  }


  unsigned int j1 = findMatch(w1InLab.Eta(), w1InLab.Phi());
  unsigned int j2 = findMatch(w2InLab.Eta(), w2InLab.Phi());
  unsigned int bb = findMatch(bInLab.Eta(),  bInLab.Phi());


  double AccJ1 =  TMath::Abs(w1InLab.Eta())<1.0 ?  
    TMath::Erf( (jetParam_.find("param0AccLightBin0"))->second * w1InLab.Pt()+(jetParam_.find("param1AccLightBin0"))->second )+(jetParam_.find("param2AccLightBin0"))->second :
    TMath::Erf( (jetParam_.find("param0AccLightBin1"))->second * w1InLab.Pt()+(jetParam_.find("param1AccLightBin1"))->second )+(jetParam_.find("param2AccLightBin1"))->second ;
  double AccJ2 =  TMath::Abs(w2InLab.Eta())<1.0 ?  
    TMath::Erf( (jetParam_.find("param0AccLightBin0"))->second * w2InLab.Pt()+(jetParam_.find("param1AccLightBin0"))->second )+(jetParam_.find("param2AccLightBin0"))->second :
    TMath::Erf( (jetParam_.find("param0AccLightBin1"))->second * w2InLab.Pt()+(jetParam_.find("param1AccLightBin1"))->second )+(jetParam_.find("param2AccLightBin1"))->second ;
  double AccB =  TMath::Abs(bInLab.Eta())<1.0 ?  
    TMath::Erf( (jetParam_.find("param0AccHeavyBin0"))->second * bInLab.Pt()+(jetParam_.find("param1AccHeavyBin0"))->second )+(jetParam_.find("param2AccHeavyBin0"))->second :
    TMath::Erf( (jetParam_.find("param0AccHeavyBin1"))->second * bInLab.Pt()+(jetParam_.find("param1AccHeavyBin1"))->second )+(jetParam_.find("param2AccHeavyBin1"))->second ;
   
  double ResolJ1 =  999.;
  if( w1InLab.Pt()>0){
    ResolJ1 = TMath::Abs(w1InLab.Eta())<1.0 ?  
      w1InLab.Pt()*TMath::Sqrt( (jetParam_.find("param0resolLightBin0"))->second*(jetParam_.find("param0resolLightBin0"))->second/w1InLab.Pt() +  
				(jetParam_.find("param1resolLightBin0"))->second*(jetParam_.find("param1resolLightBin0"))->second/w1InLab.Pt()/w1InLab.Pt()) :
      w1InLab.Pt()*TMath::Sqrt( (jetParam_.find("param0resolLightBin1"))->second*(jetParam_.find("param0resolLightBin1"))->second/w1InLab.Pt() +  
				(jetParam_.find("param1resolLightBin1"))->second*(jetParam_.find("param1resolLightBin1"))->second/w1InLab.Pt()/w1InLab.Pt()) ;
  }
  double ResolJ2 =  999.;
  if( w2InLab.Pt()>0){
    ResolJ2 = TMath::Abs(w2InLab.Eta())<1.0 ?  
      w2InLab.Pt()*TMath::Sqrt( (jetParam_.find("param0resolLightBin0"))->second*(jetParam_.find("param0resolLightBin0"))->second/w2InLab.Pt() +  
				(jetParam_.find("param1resolLightBin0"))->second*(jetParam_.find("param1resolLightBin0"))->second/w2InLab.Pt()/w2InLab.Pt()) :
      w2InLab.Pt()*TMath::Sqrt( (jetParam_.find("param0resolLightBin1"))->second*(jetParam_.find("param0resolLightBin1"))->second/w2InLab.Pt() +  
				(jetParam_.find("param1resolLightBin1"))->second*(jetParam_.find("param1resolLightBin1"))->second/w2InLab.Pt()/w2InLab.Pt()) ;
  }
  double ResolB =  999.;
  if( bInLab.Pt()>0){
    ResolB = TMath::Abs(bInLab.Eta())<1.0 ?  
      bInLab.Pt()*TMath::Sqrt( (jetParam_.find("param0resolHeavyBin0"))->second*(jetParam_.find("param0resolHeavyBin0"))->second/bInLab.Pt() +  
			       (jetParam_.find("param1resolHeavyBin0"))->second*(jetParam_.find("param1resolHeavyBin0"))->second/bInLab.Pt()/bInLab.Pt()) :
      bInLab.Pt()*TMath::Sqrt( (jetParam_.find("param0resolHeavyBin1"))->second*(jetParam_.find("param0resolHeavyBin1"))->second/bInLab.Pt() +  
			       (jetParam_.find("param1resolHeavyBin1"))->second*(jetParam_.find("param1resolHeavyBin1"))->second/bInLab.Pt()/bInLab.Pt()) ;
  }

  if(DEBUG){
    cout << "AccJ1 " << AccJ1 << endl;
    cout << "AccJ2 " << AccJ2 << endl;
    cout << "AccB "  << AccB << endl;

    cout << "ResolJ1 " << ResolJ1 << endl;
    cout << "ResolJ2 " << ResolJ2 << endl;
    cout << "ResolB "  << ResolB << endl;

  }

  TH1F* bTagLight = this->getCachedPdf(string("pdfCsvLight"));
  TH1F* bTagHeavy = this->getCachedPdf(string("pdfCsvHeavy"));

  double TFJ1 = 1.0;
  if( j1==98) 
    TFJ1 = 0.0;
  else if(j1==99)
    TFJ1 = TMath::Max(1-AccJ1, 0.0);
  else
    TFJ1 = TMath::Min( AccJ1,  1.0)*TMath::Gaus( w1InLab.Pt(), ResolJ1 );
  TFJ1 *= bTagLight->Interpolate(  bTagging_[j1] );

 double TFJ2 = 1.0;
  if( j2==98) 
    TFJ2 = 0.0;
  else if(j2==99)
    TFJ2 = TMath::Max(1-AccJ2, 0.0);
  else
    TFJ2 = TMath::Min( AccJ2,  1.0)*TMath::Gaus( w2InLab.Pt(), ResolJ2 );
  TFJ2 *= bTagLight->Interpolate(  bTagging_[j2] );

  double TFB = 1.0;
  if( bb==98) 
    TFB = 0.0;
  else if(bb==99)
    TFB = TMath::Max(1-AccB, 0.0);
  else
    TFB = TMath::Min( AccB,  1.0)*TMath::Gaus( bInLab.Pt(), ResolB );
  TFB *= bTagHeavy->Interpolate(  bTagging_[bb] );

  prob *= TFJ1;
  prob *= TFJ2;
  prob *= TFB;

  TH1F* pdfV1      = this->getCachedPdf(string("pdfV1"));
  TH1F* pdfV2      = this->getCachedPdf(string("pdfV2"));
  TH1F* pdfPtHStar = this->getCachedPdf(string("pdfPtHStar"));
  TH1F* pdfPtTTH   = this->getCachedPdf(string("pdfPtTTH"));
  TH1F* pdfBetaW   = this->getCachedPdf(string("pdfBetaW"));
  TH1F* pdfGammaW  = this->getCachedPdf(string("pdfGammaW"));


  double pdf = 
    pdfV1->Interpolate( v1 ) *
    pdfV2->Interpolate( v2 ) *
    pdfPtTTH->Interpolate( PtTTH ) *
    (1./2./TMath::Pi())    ;

  double cms = 
    pdfPtHStar->Interpolate( PtHStar ) ;

  double decay = 
    pdfBetaW->Interpolate(  BetaW  ) * 
    (1./2./TMath::Pi())  * 
    pdfGammaW->Interpolate( GammaW ) * 
    (1./2./TMath::Pi())  ;

  prob *= pdf;
  prob *= cms;
  prob *= decay;

 if(DEBUG){
   cout << "TFJ1 " << TFJ1 << endl;
   cout << "TFJ2 " << TFJ2 << endl;
   cout << "TFB  " << TFB << endl;
   debugHisto1_->Fill( w1InLab.Pt() , decay*pdf*cms);
 }

 prob = pdf*cms*decay;

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

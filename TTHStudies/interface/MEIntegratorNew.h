#ifndef MEINTEGRATORNEW_H
#define MEINTEGRATORNEW_H

#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"

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
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

#include <string>
#include <map>
#include <limits>  

#define DEBUG 0

using namespace RooFit;
using namespace std;


class MEIntegratorNew {

 public:

  MEIntegratorNew( string , int , int);
  ~MEIntegratorNew();

  double Eval(const double* ) const;  
  void   SetPar(int);
  void   setJets( vector<TLorentzVector>* );
  void   setBtag( std::vector<float>* );
  void   createMash();
  double        probability(const double*, int) const;
  unsigned int  findMatch(double, double) const;
  void   saveJetParam( string );
  void   cachePdf( string , string , int );
  void   cachePdf( string , string , string, int,    int );
  void   cachePdf( string , string , string, string, int, int, int );
  void   cachePdf( string , string , string, string, TArrayF, TArrayF, TArrayF);
  void   setMass(double);
  void   setSumEt(double);
  void   setPtPhiParam (int);
  void   resetEvaluation();

  TH1*   getCachedPdf( string ) const;
  TH1*   getCachedTF ( string ) const;
  TH2F*  getMash( );
  TH1F*  getDebugHisto( );
  void   debug();

  void   initVersors(int);
  void   initTF() ;
  void   adaptRange(TF1*, float&, float&, float, float);
  void   createTFjet(string, float , float, string , float, float);
  void   createTFmet(string , float , float, float, float);

  void   topHadEnergies    (double, double&, double&, double&, double&, int&) const;
  void   topLepEnergies    (double, double,  double&, double&, double&, double&, int&) const;
  void   topLepEnergiesFromPtPhi    (int, double, double,  double&, double&, double&, double&, int&) const;
  void   topHadLostEnergies(double, double,  double,  double&, double&, double&, double&, int&) const;
  void   higgsEnergies     (double, double&, double&, int&) const;
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
  std::map<string, TH1* > transferFunctions_;

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
  TVector3 eMEt_;

  int par_;
  int verbose_;
  int usePtPhiParam_;
  int evaluation_;
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
  double sumEt_;
  TH2F* mash_;

  TH1F* debugHisto1_;
  TH1F pdfBetaWHad_, pdfGammaWHad_, pdfBetaWLep_, pdfGammaWLep_, pdfGammaTTH_;
  TH3F pdf3D_;
  TH1F tfWjet1_, tfWjet2_, tfbHad_, tfbLep_, tfMetPt_, tfHiggs1_, tfHiggs2_;
  TH2F tfMetPhi_;

};


MEIntegratorNew::MEIntegratorNew( string fileName , int param , int verbose ) {

  cout << "Begin constructor" << endl;

  par_     = param;
  verbose_ = verbose;
  sumEt_   = 1500.; //dummy
  usePtPhiParam_ = 0;
  evaluation_    = 0;

  jets_.clear();
  bTagging_.clear();
  initVersors(0);

  Mtop_   =  174.3;
  Mb_     =  4.8;
  Mw_     =  80.19;
  pStar_  =  TMath::Sqrt( (Mtop_*Mtop_-(Mw_+Mb_)*(Mw_+Mb_) )*( Mtop_*Mtop_-(Mw_-Mb_)*(Mw_-Mb_) ) )/2./Mtop_;
  EbStar_ =  (Mtop_*Mtop_ - Mw_*Mw_ + Mb_*Mb_)/2./Mtop_;
  EWStar_ =  (Mtop_*Mtop_ + Mw_*Mw_ - Mb_*Mb_)/2./Mtop_;
  EuStar_ =  Mw_/2;
  dM2_    =  (Mtop_*Mtop_-Mb_*Mb_-Mw_*Mw_)*0.5;
  M_      =  125.;
  dMh2_   =  (M_*M_-2*Mb_*Mb_)*0.5;

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

  pdfGammaWHad_ = *((TH1F*)this->getCachedPdf("pdfGammaWHad"));
  pdfBetaWHad_  = *((TH1F*)this->getCachedPdf("pdfBetaWHad"));
  pdfBetaWLep_  = *((TH1F*)this->getCachedPdf("pdfBetaWLep"));
  pdfGammaWLep_ = *((TH1F*)this->getCachedPdf("pdfGammaWLep"));
  pdfGammaTTH_  = *((TH1F*)this->getCachedPdf("pdfGammaTTH"));
  pdf3D_        = *((TH3F*)this->getCachedPdf("pdf3D"));

  
  cout << "End constructor" << endl;

}


MEIntegratorNew::~MEIntegratorNew(){

  cout << "Start destructor" << endl;

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
  for(std::map<string, TH1*>::iterator it = transferFunctions_.begin(); it!=transferFunctions_.end(); it++){
    if(it->second) 
      delete (it->second);
  }

  delete mash_; 
  delete debugHisto1_;

  cout << "End destructor" << endl;

}

void MEIntegratorNew::initVersors(int withJetList){
  
  if(withJetList==0){
    eLep_   = TVector3(0.,0.,1.);
    eBLep_  = TVector3(0.,0.,1.);
    eBHad_  = TVector3(0.,0.,1.);
    eW1Had_ = TVector3(0.,0.,1.);
    eW2Had_ = TVector3(0.,0.,1.);
    eB1_    = TVector3(0.,0.,1.);
    eB2_    = TVector3(0.,0.,1.);
    eMEt_   = TVector3(0.,0.,1.);
    return;
  }
  else if( withJetList==1 && jets_.size() < 7 ){
    cout << "Jets are not properly initliazied!!" << endl;
    return;
  }
  else if(  withJetList==1 && jets_.size() == 8 ){
    eLep_  = (jets_[0].Vect()).Unit();
    eMEt_  = (jets_[1].Vect()).Unit();
    eBLep_ = (jets_[2].Vect()).Unit();
    eW1Had_= (jets_[3].Vect()).Unit();
    eW2Had_= (jets_[4].Vect()).Unit();
    eBHad_ = (jets_[5].Vect()).Unit();
    eB1_   = (jets_[6].Vect()).Unit();
    eB2_   = (jets_[7].Vect()).Unit();
   }
  else{cout << "Problems in MEIntegratorNew::initVersors" << endl;}
  
  return;
}


void MEIntegratorNew::adaptRange(TF1* f, float& xLow, float& xHigh, float quantile, float margin){

  double probSum[2] = {quantile, 1-quantile};
  double q[2];
  int n = f->GetQuantiles(2,q,probSum);

  if(margin<=1){
    xLow  = q[0]*(1-margin);
    xHigh = q[1]*(1+margin);
  }
  else{
    xLow  = q[0];
    xHigh = q[1];
  }
}


void MEIntegratorNew::createTFjet(string tfName, float eta, float pt, string flavor, float quantile, float margin){

  float xLow, xHigh;
  float gevStep = 2.;

  string bin = "Bin0";
  if(  TMath::Abs( eta )<1.0 ) bin = "Bin0";
  else bin = "Bin1";
  
  TF1* tf = new TF1("tf",Form("TMath::Gaus( x, [0] , [0]*TMath::Sqrt( %f/[0] + %f/[0]/[0]) , 1) ", //(TMath::Erf(%f*x+%f)+%f) *
			      //(jetParam_.find("param0AccLight"+bin))->second, 
			      //(jetParam_.find("param1AccLight"+bin))->second,
			      //(jetParam_.find("param2AccLight"+bin))->second,
			      (jetParam_.find("param0resol"+flavor+bin))->second*(jetParam_.find("param0resol"+flavor+bin))->second,
			      (jetParam_.find("param1resol"+flavor+bin))->second*(jetParam_.find("param1resol"+flavor+bin))->second
			      ),1,1000);
  tf->SetNpx(1000);
  tf->SetParameter(0, pt );
  adaptRange(tf, xLow, xHigh, quantile, margin);

  TH1F* htf_ = new TH1F(("htf"+tfName).c_str(), "", int((xHigh-xLow)/gevStep), xLow, xHigh);
  for( int i = 1; i <= htf_->GetNbinsX(); i++){
    float binC = htf_->GetBinCenter(i);
    tf->SetParameter(0, binC);
    htf_->SetBinContent(i, tf->Eval( pt ) );
  }
  if( transferFunctions_.find("tf"+tfName)!=transferFunctions_.end()){
    delete (transferFunctions_.find("tf"+tfName)->second);
    transferFunctions_.erase( transferFunctions_.find("tf"+tfName) );
    transferFunctions_["tf"+tfName] = htf_;
  }
  else
    transferFunctions_["tf"+tfName] = htf_;

  if(tfName.find("Wjet1")!=string::npos)
    tfWjet1_ = *htf_;
  else if (tfName.find("Wjet2")!=string::npos)
    tfWjet2_ = *htf_;
  else if (tfName.find("bHad")!=string::npos)
    tfbHad_ = *htf_;
  else if (tfName.find("bLep")!=string::npos)
    tfbLep_ = *htf_;
  else if (tfName.find("Higgs1")!=string::npos)
    tfHiggs1_ = *htf_;
  else if (tfName.find("Higgs2")!=string::npos)
    tfHiggs2_ = *htf_;
  else{
    if(verbose_) cout << "Name in createTFjet is not valid" << endl;
  }

  delete tf;
}


void MEIntegratorNew::createTFmet(string tfName, float phi, float pt, float quantile, float margin){

  float xLowEt, xHighEt, xLowPhi, xHighPhi;
  float gevStep = 4.;
  float etaStep = 0.04;

  string bin = "Bin0";
  if (sumEt_ < 1200.) 
    bin =  "Bin0";
  else if ( sumEt_ > 1200. && sumEt_ < 1800.)
    bin =  "Bin1";
  else 
    bin =  "Bin2";

  TF1* tfMetPt   = new TF1(("tf"+tfName+"Pt").c_str(), Form("TMath::Gaus(x, (%f + %f*TMath::Exp(%f*[0]+%f) ),  [0]*TMath::Sqrt(%f/[0] + %f/[0]/[0])  , 1)",
					   (jetParam_.find("param0EtMean"+bin))->second, (jetParam_.find("param1EtMean"+bin))->second, (jetParam_.find("param2EtMean"+bin))->second,(jetParam_.find("param3EtMean"+bin))->second,
					   (jetParam_.find("param0EtWidth"+bin))->second*(jetParam_.find("param0EtWidth"+bin))->second,(jetParam_.find("param1EtWidth"+bin))->second*(jetParam_.find("param1EtWidth"+bin))->second
					   ),   // x = reco - gen 
			   -1000.,1000.
			   );
  tfMetPt->SetNpx(2000);
  tfMetPt->SetParameter(0, pt );
  adaptRange(tfMetPt, xLowEt, xHighEt, quantile, margin);
  if(xLowEt>0){
    xLowEt  += pt;
    xHighEt += pt;
  }
  else{
    xLowEt  =  0.;
    xHighEt += pt;
  }


  TH1F* htfMetPt_ = new TH1F(("htf"+tfName+"Pt").c_str(), "", int((xHighEt-xLowEt)/gevStep), xLowEt, xHighEt);
  for( int i = 1; i <= htfMetPt_->GetNbinsX(); i++){
    float binC = htfMetPt_->GetBinCenter(i);
    tfMetPt->SetParameter(0, binC);
    htfMetPt_->SetBinContent(i, tfMetPt->Eval( pt - binC) );
  }
  if( transferFunctions_.find("tf"+tfName+"Pt")!=transferFunctions_.end()){
    delete (transferFunctions_.find("tf"+tfName+"Pt")->second);
    transferFunctions_.erase( transferFunctions_.find("tf"+tfName+"Pt") );
    transferFunctions_["tf"+tfName+"Pt"] = htfMetPt_;
  }
  else
    transferFunctions_["tf"+tfName+"Pt"] = htfMetPt_;

  tfMetPt_ = *htfMetPt_;
  delete tfMetPt;


  // check norm...
  TF1* tfMetPhi   = new TF1(("tf"+tfName+"Phi").c_str(),Form("(2./(TMath::Erf(TMath::Pi()/ (%f/[0] + %f/[0]/[0]) )))*TMath::Gaus(x, 0.0, %f/[0] + %f/[0]/[0] ,1)", 
							     (jetParam_.find("param0PhiWidth"+bin))->second,(jetParam_.find("param1PhiWidth"+bin))->second,
							     (jetParam_.find("param0PhiWidth"+bin))->second,(jetParam_.find("param1PhiWidth"+bin))->second
							     ),   // x = |reco-gen| 
			    0,TMath::Pi()); 

  tfMetPhi->SetNpx(1000);
  tfMetPhi->SetParameter(0, pt );
  adaptRange(tfMetPhi, xLowPhi, xHighPhi, quantile, margin );
  if((TMath::Pi() - xHighPhi) < 0.2) xHighPhi =  TMath::Pi();

 
  TH2F* htfMetPhi_ = new TH2F(("htf"+tfName+"Phi").c_str(), "", int((xHighEt-xLowEt)/gevStep), xLowEt, xHighEt, int((xHighPhi-xLowPhi)/etaStep), xLowPhi, xHighPhi);
  for( int i = 1; i <= htfMetPhi_->GetNbinsX(); i++){
    float binCX = htfMetPhi_->GetXaxis()->GetBinCenter(i);
    tfMetPhi->SetParameter(0, binCX);
     for( int j = 1; j <= htfMetPhi_->GetNbinsY(); j++){
       float binCY = htfMetPhi_->GetYaxis()->GetBinCenter(j);
       //htfMetPhi_->SetBinContent(i,j, tfMetPhi->Eval( TMath::ACos(TMath::Cos( phi - binCY)) ) );
       htfMetPhi_->SetBinContent(i,j, tfMetPhi->Eval( binCY ) );
     }
  }
  if( transferFunctions_.find("tf"+tfName+"Phi")!=transferFunctions_.end()){
    delete (transferFunctions_.find("tf"+tfName+"Phi")->second);
    transferFunctions_.erase( transferFunctions_.find("tf"+tfName+"Phi") );
    transferFunctions_["tf"+tfName+"Phi"] = htfMetPhi_;
  }
  else
    transferFunctions_["tf"+tfName+"Phi"] = htfMetPhi_;

  tfMetPhi_ = *htfMetPhi_;
  delete tfMetPhi;



}




 
void MEIntegratorNew::initTF(){
  
  createTFjet("Wjet1",  eW1Had_.Eta(), jets_[3].E(), "Light", 0.025, 0.30);
  createTFjet("Wjet2",  eW2Had_.Eta(), jets_[4].E(), "Light", 0.025, 0.30);
  createTFjet("bHad",   eBHad_.Eta(),  jets_[5].E(), "Heavy", 0.025, 0.30);
  createTFjet("bLep",   eBLep_.Eta(),  jets_[2].E(), "Heavy", 0.025, 0.30);
  createTFjet("Higgs1", eB1_.Eta(),    jets_[6].E(), "Heavy", 0.025, 0.30);
  createTFjet("Higgs2", eB2_.Eta(),    jets_[7].E(), "Heavy", 0.025, 0.30);

  createTFmet("Met", jets_[1].Phi() , jets_[1].Pt() , 0.025, 0.50);

}



void MEIntegratorNew::setMass(double mass){
  M_    = mass;
  dMh2_ = (M_*M_-2*Mb_*Mb_)*0.5; 
}

void MEIntegratorNew::setSumEt(double sumEt){
  sumEt_ = sumEt;
}

void MEIntegratorNew::setPtPhiParam(int usePtPhiParam){
  usePtPhiParam_ =  usePtPhiParam;
}

void MEIntegratorNew::topHadEnergies(double E1, double& E2, double& E3, double& cos1, double& cos2, int& errFlag ) const {
  
  TVector3 e1 = eW1Had_;
  TVector3 e2 = eW2Had_;
  TVector3 e3 = eBHad_;

  double a12 = e1.Angle(e2);
  double a13 = e1.Angle(e3);
  double a23 = e2.Angle(e3);
  
  E2 = Mw_*Mw_/E1/(4*TMath::Sin(a12/2.)*TMath::Sin(a12/2.));

  double a = E1+E2;
  double b = E1*TMath::Cos(a13)+E2*TMath::Cos(a23);

  if( dM2_*dM2_ - (a*a - b*b)*Mb_*Mb_ < 0){
    errFlag = 1;
    return;
  }
  
  double E3_1 =  (a*dM2_ + b*TMath::Sqrt(dM2_*dM2_ - (a*a - b*b)*Mb_*Mb_))/(a*a - b*b) ;
  //double E3_2 =  (a*dM2_ - b*TMath::Sqrt(dM2_*dM2_ - (a*a - b*b)*Mb_*Mb_))/(a*a - b*b) ;
  
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


void MEIntegratorNew::topLepEnergiesFromPtPhi(int sign , double nuPhi, double nuPt, double& Enu, double& Eb, double& cos1, double& cos2, int& errFlag ) const{
  

  TVector3 e3T( nuPt*TMath::Cos(nuPhi), nuPt*TMath::Sin(nuPhi) , 0.);
  TVector3 e1T( jets_[0].Px(),  jets_[0].Py(), 0.);
    
  double rho = e1T.Dot(e3T) + Mw_*Mw_/2.;

  double a   = jets_[0].Pz()/jets_[0].P();
  double b   = rho/jets_[0].P();
  double c2 =  nuPt*nuPt - b*b;     

  if( (a*a*b*b - c2*(1-a*a))<0 ){
    errFlag = 1;
    return;
  }

  double pz1 = (a*b + sqrt(a*a*b*b - c2*(1-a*a)))/(1-a*a);
  double pz2 = (a*b - sqrt(a*a*b*b - c2*(1-a*a)))/(1-a*a);

  double pz = sign==1 ? pz1 : pz2;
  Enu = TMath::Sqrt(nuPt*nuPt + pz*pz);
  double nuTheta = pz/Enu;

  topLepEnergies( nuPhi, nuTheta, Enu, Eb, cos1, cos2, errFlag  );

  return;
}



void MEIntegratorNew::topLepEnergies(double nuPhi, double nuTheta, double& Enu, double& Eb, double& cos1, double& cos2, int& errFlag ) const{

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

  if( dM2_*dM2_ - (a*a - b*b)*Mb_*Mb_ < 0){
    errFlag = 1;
    return;
  }

  double E3_1 =  (a*dM2_ + b*TMath::Sqrt(dM2_*dM2_ - (a*a - b*b)*Mb_*Mb_))/(a*a - b*b) ;
  //double E3_2 =  (a*dM2_ - b*TMath::Sqrt(dM2_*dM2_ - (a*a - b*b)*Mb_*Mb_))/(a*a - b*b) ;

  Eb = E3_1;

  if(Eb<Mb_){
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



void MEIntegratorNew::topHadLostEnergies(double missPhi, double missTheta, double missE, double& E1, double& Eb, double& cos1, double& cos2, int& errFlag ) const{

}


void MEIntegratorNew::higgsEnergies(double E1, double& E2, double& cos1, int& errFlag) const{

  if(E1<Mb_){
    errFlag = 1;
    return;
  }

  TVector3 e1 = eB1_;
  TVector3 e2 = eB2_;

  double a12 = e1.Angle(e2); // b1 - b2

  double a = E1;
  double b = TMath::Sqrt(E1*E1 - Mb_*Mb_)*TMath::Cos(a12);

  if( dMh2_*dMh2_ - (a*a - b*b)*Mb_*Mb_ < 0){
    errFlag = 1;
    return;
  }

  double E2_1 = (a*dMh2_ + b*TMath::Sqrt(dMh2_*dMh2_ - (a*a - b*b)*Mb_*Mb_))/(a*a - b*b);
  //double E2_2 = (a*dMh2_ - b*TMath::Sqrt(dMh2_*dMh2_ - (a*a - b*b)*Mb_*Mb_))/(a*a - b*b);

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

  //evaluation_++;

  double prob = 0.0;
  if(usePtPhiParam_==0)
    prob = probability(x,0);     
  else{
    prob += probability(x,-1);
    prob += probability(x,+1);
  }

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
    
    TH1F* hCache = new TH1F(pdfName.c_str(),"", nBins, var->getMin(), var->getMax() );
    
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




TH1* MEIntegratorNew::getCachedPdf( string pdfName ) const {

  if( variables1D_.find(pdfName) !=  variables1D_.end() ) 
    return (variables1D_.find(pdfName))->second;
  else if( variables2D_.find(pdfName) !=  variables2D_.end() ) 
    return (variables2D_.find(pdfName))->second;
  else if( variables3D_.find(pdfName) !=  variables3D_.end() ) 
    return (variables3D_.find(pdfName))->second;

  else return 0;
}


TH1* MEIntegratorNew::getCachedTF( string tfName ) const {

  if( transferFunctions_.find(tfName) != transferFunctions_.end() ) 
    return (transferFunctions_.find(tfName))->second;

  else return 0;
}

 
double MEIntegratorNew::topHadJakobi( double Eq1, double Eq2, double EbHad)  const{
  return EbHad * Eq2 * Eq1 ;
}

double MEIntegratorNew::topLepJakobi( double Enu, double EbLep)  const{
  return EbLep * Enu ;
}

double MEIntegratorNew::higgsJakobi ( double Eh1, double Eh2)  const{
  return Eh1 * Eh2 ;
}


double MEIntegratorNew::topHadDensity ( double cos1, double cos2) const{

  double val1 = const_cast<TH1F*>(&pdfBetaWHad_)->Interpolate(cos1);
  double val2 = const_cast<TH1F*>(&pdfGammaWHad_)->Interpolate(cos2);
  double res  = val1*val2;

  return res;

}

double MEIntegratorNew::topLepDensity ( double cos1, double cos2)  const{

  double val1 = const_cast<TH1F*>(&pdfBetaWLep_)->Interpolate(cos1);
  double val2 = const_cast<TH1F*>(&pdfGammaWLep_)->Interpolate(cos2);

  return val1*val2;

}

double MEIntegratorNew::higgsDensity ( double cos1 )  const{
  return 1.0;
}

double MEIntegratorNew::tthDensity ( double Q, double m12, double cos1Star, double cos3)  const{

  double p1Star = (m12*m12 - 4*Mtop_*Mtop_) >=0 ? 
    TMath::Sqrt((m12*m12 - 4*Mtop_*Mtop_)*(m12*m12))/2./m12 : -999.;
  double p3     = (Q*Q - (m12+M_)*(m12+M_)) >=0 && (Q*Q - (m12 - M_)*(m12 - M_)) >=0 ? 
    TMath::Sqrt((Q*Q - (m12+M_)*(m12+M_) )*( Q*Q - (m12 - M_)* (m12 - M_)))/2./Q : -999;

  if( TMath::IsNaN(p1Star) || TMath::IsNaN(p3) || p1Star<=0 || p3 <=0){
    if(verbose_) cout << "tthDensity evaluated for non-acceptable values of p1Star and/or p3" << endl;
    return 0.0;
  }

  double val1 = const_cast<TH3F*>(&pdf3D_)->Interpolate(Q, m12, cos1Star);
  double val2 = const_cast<TH1F*>(&pdfGammaTTH_)->Interpolate( cos3 );

  return (val1*val2)/p1Star/p3;

}


void MEIntegratorNew::resetEvaluation(){
  evaluation_ = 0;
}


double MEIntegratorNew::probability(const double* x, int sign) const{

  int errFlag = 0;
  double prob = 1.0;
  
  double Eq1         = x[0];
  double cosThetaNu  = x[1]; // or nuPt
  double phiNu       = x[2];
  double Eh1         = x[3];

  // trasnform phiNu into an absolute phi:
  phiNu += eMEt_.Phi();
  if( phiNu < -TMath::Pi() ) phiNu += (2*TMath::Pi());
  if( phiNu >  TMath::Pi())  phiNu -= (2*TMath::Pi());


  if(verbose_){
    cout << "#" << evaluation_ << endl;
    cout << "Eq1 = "         << Eq1 << endl; 
    cout << "cosThetaNu = " << cosThetaNu << endl; 
    cout << "phiNu = "      << phiNu  << endl; 
    cout << "Eh1 = "         << Eh1  << endl; 
    cout << "********" << endl;
  }

  double Eq2, EbHad, cos1Had, cos2Had;
  topHadEnergies( Eq1, Eq2, EbHad, cos1Had, cos2Had, errFlag );

  if(errFlag){
    if(verbose_) cout << "Problems with topHadEnergies" << endl;
    return 0.0;
  }

  TVector3 eNu(0.,0.,1.); // neutrino
  double Enu, EbLep, cos1Lep, cos2Lep;
  if(usePtPhiParam_==0){
    topLepEnergies( phiNu, cosThetaNu, Enu, EbLep, cos1Lep, cos2Lep, errFlag  );
    eNu.SetTheta( TMath::ACos( cosThetaNu ) ); 
    eNu.SetPhi  ( phiNu );   
    eNu.SetMag  ( 1.); 
  }
  else{
    topLepEnergiesFromPtPhi(sign, phiNu, cosThetaNu, Enu, EbLep, cos1Lep, cos2Lep, errFlag );
    eNu.SetTheta( TMath::ASin(cosThetaNu/Enu) ); 
    eNu.SetPhi  ( phiNu );   
    eNu.SetMag  ( 1.); 
    prob *= (cosThetaNu/Enu/Enu); // Jacobian  <--------- CORRECT !!!!
  }

  if(errFlag){
    if(verbose_) cout << "Problems with topLepEnergies" << endl;
    return 0.0;
  }
 

  double Eh2, cos1Higgs;
  higgsEnergies( Eh1, Eh2, cos1Higgs, errFlag );
  
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

  ///double MEpart = 
  //topHadDensity(cos1Had,cos2Had) * topLepDensity(cos1Lep,cos2Lep) * higgsDensity(cos1Higgs) *
  //topHadJakobi( Eq1, Eq2, EbHad) * topLepJakobi(Enu, EbLep) * higgsJakobi( Eh1, Eh2 ) *
  //tthDensity( Q , m12, cos1Star, cos3);

  double TFpart = 
    const_cast<TH1F*>(&tfWjet1_)  ->Interpolate( W1Had.E() ) *
    const_cast<TH1F*>(&tfWjet2_)  ->Interpolate( W2Had.E() ) *
    const_cast<TH1F*>(&tfbHad_)   ->Interpolate( bHad.E()  ) //*
    //const_cast<TH1F*>(&tfbLep_)   ->Interpolate( bLep.E()  ) *
    //const_cast<TH1F*>(&tfHiggs1_) ->Interpolate( higgs1.E()) *
    //const_cast<TH1F*>(&tfHiggs2_) ->Interpolate( higgs2.E()) //*
    //const_cast<TH1F*>(&tfMetPt_)  ->Interpolate( WLepNu.Pt()) *
    //const_cast<TH2F*>(&tfMetPhi_) ->Interpolate( WLepNu.Pt(), TMath::ACos(TMath::Cos(WLepNu.Phi()-jets_[1].Phi())) ) 
    ;

  //prob *= MEpart;
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
  for(std::map<string, TH2F*>::iterator it = variables2D_.begin(); it!=variables2D_.end(); it++){
    if(it->second) 
      cout << it->first << " => " << (it->second)->Integral() << endl;
    else
      cout << it->first << " => null pointer" << endl;
  }
  for(std::map<string, TH3F*>::iterator it = variables3D_.begin(); it!=variables3D_.end(); it++){
    if(it->second) 
      cout << it->first << " => " << (it->second)->Integral() << endl;
    else
      cout << it->first << " => null pointer" << endl;
  }


  cout << "*** debug end *** " << endl;

}


#endif

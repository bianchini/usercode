#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>

#include "TSystem.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"
#include "TVectorD.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"


#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooLandau.h"
#include "RooUniform.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooKeysPdf.h"
#include "RooBifurGauss.h"
#include "RooVoigtian.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooGaussian.h"
#include "RooAbsCategory.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "RooTruthModel.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooBernstein.h"
#include "RooNDKeysPdf.h"
#include "RooChiSquarePdf.h"

#include "Bianchi/TTHStudies/interface/MEIntegratorNew.h"

#define GENJETDR 0.3
#define VERBOSE  false

using namespace std;
using namespace RooFit;



int main(int argc, const char* argv[])
{

  std::cout << "TestMENew" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");
  //gSystem->Load("BianchiTTHStudiesPlugins");

  AutoLibraryLoader::enable();

  TFile* fout = 0;

  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");
  std::string outFileName( in.getParameter<std::string>  ("outFileName" ) );
  std::string pathToFile(  in.getParameter<std::string> ("pathToFile"   ) );
  int vegasPoints(  in.getParameter<int> ("vegasPoints"   ) );
  bool verbose   (  in.getParameter<bool>("verbose") );

  float met(  in.getParameter<double> ("met") );


  float pertW1(    in.getParameter<double> ("pertW1") );
  float pertW2(    in.getParameter<double> ("pertW2") );
  float pertBHad(  in.getParameter<double> ("pertBHad") );
  float pertBLep(  in.getParameter<double> ("pertBLep") );
  float enlargeE1(  in.getParameter<double>  ("enlargeE1") );
  float enlargeEh1(  in.getParameter<double> ("enlargeEh1") );
  float enlargePt(  in.getParameter<double>  ("enlargePt") );


  int par = 4;
  MEIntegratorNew* meIntegrator = new MEIntegratorNew( pathToFile , par, int(verbose));
  //meIntegrator->debug();

  vector<TLorentzVector> jets;
  TLorentzVector jet1,jet2,jet3,jet4,jet5,jet6,jet7,jet8;

  jet1.SetPtEtaPhiM(19.8454,            0.822914,-2.3529,0.10566);     // lep
  jet2.SetPtEtaPhiM(91.934,             0.0,-0.0161266,-2.33602e-06);  // MET
  jet3.SetPtEtaPhiM(175.197 * pertBLep, 0.0936224,-0.969521,4.8);      // b from top lep
  jet4.SetPtEtaPhiM(122.914 * pertW1,   1.17562,2.17839,0.33);         // j1 from W
  jet5.SetPtEtaPhiM(34.9939 * pertW2,   1.06994,0.879248,0.33);        // j2 from W
  jet6.SetPtEtaPhiM(95.5899 * pertBHad, 0.829954,3.10669,4.8);         // b from top hadr
  jet7.SetPtEtaPhiM(32.805,             0.678347,-0.0951163,4.8);      // b1 from H
  jet8.SetPtEtaPhiM(75.4309,            2.08021,2.75503,4.8);          // b2 from H

  jets.push_back( jet1 );  
  jets.push_back( jet2 );  
  jets.push_back( jet3 );
  jets.push_back( jet4 );  
  jets.push_back( jet5 );  
  jets.push_back( jet6 );
  jets.push_back( jet7 );
  jets.push_back( jet8 );

  /////////////////////////////////////////////////////////////////////////////

  meIntegrator->SetPar(par);
  meIntegrator->setJets(&jets);

  meIntegrator->createMash();
  meIntegrator->setMass( met );
  meIntegrator->setSumEt(1500.);
  meIntegrator->initVersors(1);
  meIntegrator->initTF();
  meIntegrator->setPtPhiParam(0);
  
  /////////////////////////////////////////////////////////////////////////////

  double E1low   = (meIntegrator->getCachedTF("tfWjet1"))->GetXaxis()->GetXmin() * (1-enlargeE1);
  double E1high  = (meIntegrator->getCachedTF("tfWjet1"))->GetXaxis()->GetXmax() * (1+enlargeE1);

  double Ptlow   = (meIntegrator->getCachedTF("tfMetPt"))->GetXaxis()->GetXmin() * (1-enlargePt);
  double Pthigh  = (meIntegrator->getCachedTF("tfMetPt"))->GetXaxis()->GetXmax() * (1+enlargePt);

  Ptlow  = -1.;
  Pthigh = +1.;

  double Philow   = -(meIntegrator->getCachedTF("tfMetPhi"))->GetYaxis()->GetXmax();
  double Phihigh  =  (meIntegrator->getCachedTF("tfMetPhi"))->GetYaxis()->GetXmax();

  double Eh1low   = (meIntegrator->getCachedTF("tfHiggs1"))->GetXaxis()->GetXmin() * (1-enlargeEh1);
  double Eh1high  = (meIntegrator->getCachedTF("tfHiggs1"))->GetXaxis()->GetXmax() * (1+enlargeEh1);

  /////////////////////////////////////////////////////////////////////////////

  cout << "E1  : ["  << E1low  << "," << E1high  << "]" << endl;
  cout << "Pt  : ["  << Ptlow  << "," << Pthigh  << "]" << endl;
  cout << "Phi : ["  << Philow << "," << Phihigh << "]" << endl;
  cout << "Eh1 : ["  << Eh1low << "," << Eh1high << "]" << endl;

  /////////////////////////////////////////////////////////////////////////////

  double xL[4] = {  E1low,  Ptlow,   Philow,   Eh1low};
  double xU[4] = {  E1high, Pthigh , Phihigh,  Eh1high};

  ROOT::Math::GSLMCIntegrator ig2("vegas", 1.e-12, 1.e-5, vegasPoints);
  ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par); 
  ig2.SetFunction(toIntegrate);

  double mH[7] = {100, 110, 120, 125, 130, 140, 150};

  for(int m = 0; m<7 ; m++){
    meIntegrator->setMass( mH[m] );
    double p = ig2.Integral(xL, xU);
    cout << "Mass " << mH[m] << " => prob  = " << p << endl;
  }

  fout = new TFile(outFileName.c_str(),"RECREATE");
  fout->cd();
  (meIntegrator->getCachedPdf("pdfGammaWHad"))->Write("pdfGammaWHad",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfBetaWHad"))->Write("pdfBetaWHad",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfGammaWLep"))->Write("pdfGammaWLep",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfBetaWLep"))->Write("pdfBetaWLep",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfGammaTTH"))->Write("pdfGammaTTH",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdf3D"))->Write("pdf3D",TObject::kOverwrite);
  (meIntegrator->getCachedTF("tfWjet1"))->Write("",TObject::kOverwrite);
  (meIntegrator->getCachedTF("tfWjet2"))->Write("",TObject::kOverwrite);  
  (meIntegrator->getCachedTF("tfbHad"))->Write("",TObject::kOverwrite);
  (meIntegrator->getCachedTF("tfbLep"))->Write("",TObject::kOverwrite);
  (meIntegrator->getCachedTF("tfHiggs1"))->Write("",TObject::kOverwrite);
  (meIntegrator->getCachedTF("tfHiggs2"))->Write("",TObject::kOverwrite);
  (meIntegrator->getCachedTF("tfMetPt"))->Write("",TObject::kOverwrite);
  (meIntegrator->getCachedTF("tfMetPhi"))->Write("",TObject::kOverwrite);
  fout->Close();
  delete fout;


  delete meIntegrator;
  return 0;
}

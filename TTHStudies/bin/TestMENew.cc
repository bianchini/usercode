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

  float met(  in.getParameter<double> ("met") );

  int par = 4;
  MEIntegratorNew* meIntegrator = new MEIntegratorNew( pathToFile , par, 1);
  meIntegrator->debug();

  vector<TLorentzVector> jets;
  TLorentzVector jet1,jet2,jet3,jet4,jet5,jet6,jet7,jet8;
  jet1.SetPtEtaPhiM(100,   0.0,   0.0,    0. );    // lep
  jet2.SetPtEtaPhiM(met,   0.0,   2.0,    0. );    // MET
  jet3.SetPtEtaPhiM(70,    1.0,  -1.0,   10  );    // b from top lep
  jet4.SetPtEtaPhiM(100,   -1.0,  -2.0,   10  );    // j1 from W
  jet5.SetPtEtaPhiM(150,   -0.5,  -1.0,   10  );    // j2 from W
  jet6.SetPtEtaPhiM(49.2, -2.14,  0.421,  0.1);    // b from top hadr
  jet7.SetPtEtaPhiM(170.8, -3.16, -0.967,  0.1);    // b1 from H
  jet8.SetPtEtaPhiM(54.7, -3.4,   2.96,     5);    // b2 from H

  jets.push_back( jet1 );  
  jets.push_back( jet2 );  
  jets.push_back( jet3 );
  jets.push_back( jet4 );  
  jets.push_back( jet5 );  
  jets.push_back( jet6 );
  jets.push_back( jet7 );
  jets.push_back( jet8 );

  vector<float> bTagging;
  bTagging.push_back( 0.0 ) ; 
  bTagging.push_back( 0.0 ) ; 
  bTagging.push_back( 0.0 ) ; 
  bTagging.push_back( 0.0 ) ; 
  bTagging.push_back( 0.0 ) ; 
  bTagging.push_back( 0.0 ) ; 
  bTagging.push_back( 0.0 ) ; 
  bTagging.push_back( 0.0 ) ; 


  double xL[4] = {  0., -1., -TMath::Pi(),   0.};
  double xU[4] = {500.,  1.,  TMath::Pi(), 500.};

  ROOT::Math::GSLMCIntegrator ig2("vegas", 1.e-12, 1.e-5, vegasPoints);
  ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par); 
  meIntegrator->SetPar(par);
  meIntegrator->setJets(&jets);
  meIntegrator->setBtag(&bTagging);
  meIntegrator->createMash();
  meIntegrator->setMass(125);
  meIntegrator->setSumEt(1500.);
  meIntegrator->initVersors(1);
  meIntegrator->initTF();

  ig2.SetFunction(toIntegrate);


  //double p = ig2.Integral(xL, xU);
  //cout << "Prob  = " << p << endl;

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

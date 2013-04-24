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

#include "Bianchi/TTHStudies/interface/MEIntegrator.h"

#define GENJETDR 0.3
#define VERBOSE  false

using namespace std;
using namespace RooFit;



int main(int argc, const char* argv[])
{

  std::cout << "TestME" << std::endl;
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


  MEIntegrator* meIntegrator = new MEIntegrator( pathToFile );
  meIntegrator->debug();

  fout = new TFile("TestME.root","UPDATE");
  fout->cd();
  (meIntegrator->getCachedPdf("pdfPtTTH"))->Write("pdfPtTTH",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfV2"))->Write("pdfV2",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfPtHStar"))->Write("pdfPtHStar",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfBetaW"))->Write("pdfBetaW",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfGammaW"))->Write("pdfGammaW",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfCsvHeavy"))->Write("pdfCsvHeavy",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfCsvLight"))->Write("pdfCsvLight",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfV1"))->Write("pdfV1",TObject::kOverwrite);
  fout->Close();
  delete fout;


  delete meIntegrator;
  return 0;
}

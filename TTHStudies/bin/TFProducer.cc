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
#include "TH3F.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"
#include "TVectorD.h"
#include "TArrayF.h"

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

#define GENJETDR 0.3
#define VERBOSE  false

using namespace std;
using namespace RooFit;



int main(int argc, const char* argv[])
{

  std::cout << "TFProducer" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();

  TFile* fout = 0;

  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");
  std::string outFileName( in.getParameter<std::string>("outFileName" ) );
  std::string pathToFile(  in.getParameter<std::string> ("pathToFile"  ) );

  TFile* file = TFile::Open(pathToFile.c_str(),"READ");
  if(!file || file->IsZombie()){
    cout << "No input file..." << endl;
    return 0;
  }
  TTree* tree = (TTree*)file->Get("genTree");
  if(!tree){
    cout << "The tree called 'genTree' is not there..."<< endl;
    return 0;
  }

  TTree* treeEvent = (TTree*)file->Get("genEventTree");
  if(!treeEvent){
    cout << "The tree called 'genEventTree' is not there..."<< endl;
    return 0;
  }

  TTree* treeJetsLight = (TTree*)file->Get("genJetLightTree");
  if(!treeJetsLight){
    cout << "The tree called 'genJetLightTree' is not there..."<< endl;
    return 0;
  }
  TTree* treeJetsHeavy = (TTree*)file->Get("genJetHeavyTree");
  if(!treeJetsHeavy){
    cout << "The tree called 'genJetHeavyTree' is not there..."<< endl;
    return 0;
  }


  float etaBinning[] = {0.0, 1.0, 2.5};

  TH1F* accLightBin0  = new TH1F("accLightBin0","",100,0,400);
  TH1F* accLightBin1  = new TH1F("accLightBin1","",100,0,400);
  TH1F* accHeavyBin0  = new TH1F("accHeavyBin0","",100,0,400);
  TH1F* accHeavyBin1  = new TH1F("accHeavyBin1","",100,0,400);

  TH1F* resolLightBin0  = new TH1F("resolLightBin0","",20,0,200);
  TH1F* respLightBin0   = new TH1F("respLightBin0","", 20,0,200);
  TH1F* resolLightBin1  = new TH1F("resolLightBin1","",20,0,200);
  TH1F* respLightBin1   = new TH1F("respLightBin1","", 20,0,200);

  TH1F* resolHeavyBin0  = new TH1F("resolHeavyBin0","",20,0,200);
  TH1F* respHeavyBin0   = new TH1F("respHeavyBin0","", 20,0,200);
  TH1F* resolHeavyBin1  = new TH1F("resolHeavyBin1","",20,0,200);
  TH1F* respHeavyBin1   = new TH1F("respHeavyBin1","", 20,0,200);

  for(int k = 0; k < 2; k++){

    TH1F* accPass   = (TH1F*)accLightBin0->Clone(Form("accPass%d",k)); //accPass->Sumw2();
    TH1F* accTot    = (TH1F*)accLightBin0->Clone(Form("accTot%d",k));  //accTot->Sumw2();

    TH2F* hCorrLight = new TH2F("hCorrLight","", 20, 0,200, 40, 0,200);
    treeJetsLight->Draw("ptReco:pt>>hCorrLight", Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",etaBinning[k],etaBinning[k+1]) );

    treeJetsLight->Draw(Form("pt>>accTot%d",k),  Form("(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",            etaBinning[k],etaBinning[k+1]) );
    treeJetsLight->Draw(Form("pt>>accPass%d",k), Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",etaBinning[k],etaBinning[k+1]) );
    accPass->Sumw2();
    accTot->Sumw2();

    for(int j = 1; j<=hCorrLight->GetNbinsX(); j++){

      TH1F* hSlice = (TH1F*)hCorrLight->ProjectionY("_py",j,j);
      if(hSlice->GetEntries()<8) continue;
      cout << hSlice->GetMean() << ", " << hSlice->GetRMS() << endl;

      TF1* resol = new TF1("resol","gaus", 30,300);
      hSlice->Fit(resol);

      if(k==0){
	resolLightBin0->SetBinContent(j, resol->GetParameter(2));
	resolLightBin0->SetBinError  (j, resol->GetParameter(2)/TMath::Sqrt(hSlice->GetEntries()));
	respLightBin0->SetBinContent (j, resol->GetParameter(1));
	respLightBin0->SetBinError   (j, resol->GetParameter(1)/TMath::Sqrt(hSlice->GetEntries()));
      }
      if(k==1){
	resolLightBin1->SetBinContent(j, resol->GetParameter(2));
	resolLightBin1->SetBinError  (j, resol->GetParameter(2)/TMath::Sqrt(hSlice->GetEntries()));
	respLightBin1->SetBinContent (j, resol->GetParameter(1));
	respLightBin1->SetBinError   (j, resol->GetParameter(1)/TMath::Sqrt(hSlice->GetEntries()));
      }

      delete resol;
    }


    if(k==0)
      accLightBin0->Divide( accPass,accTot, 1.0,1.0);
    if(k==1)
      accLightBin1->Divide( accPass,accTot, 1.0,1.0);

    accPass->Reset();
    accTot->Reset();

    TH2F* hCorrHeavy = new TH2F("hCorrHeavy","", 20, 0,200, 40, 0,200);
    treeJetsHeavy->Draw("ptReco:pt>>hCorrHeavy", Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",etaBinning[k],etaBinning[k+1]) );

    treeJetsHeavy->Draw(Form("pt>>accTot%d",k),  Form("(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",            etaBinning[k],etaBinning[k+1]) );
    treeJetsHeavy->Draw(Form("pt>>accPass%d",k), Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",etaBinning[k],etaBinning[k+1]) );
    accPass->Sumw2();
    accTot->Sumw2();

    for(int j = 1; j<=hCorrHeavy->GetNbinsX(); j++){

      TH1F* hSlice = (TH1F*)hCorrHeavy ->ProjectionY("_py",j,j);
      if(hSlice->GetEntries()<8) continue;
      cout << hSlice->GetMean() << ", " << hSlice->GetRMS() << endl;

      TF1* resol = new TF1("resol","gaus", 30,300);
      hSlice->Fit(resol);

      if(k==0){
	resolHeavyBin0->SetBinContent(j, resol->GetParameter(2));
	resolHeavyBin0->SetBinError  (j, resol->GetParameter(2)/TMath::Sqrt(hSlice->GetEntries()));
	respHeavyBin0->SetBinContent (j, resol->GetParameter(1));
	respHeavyBin0->SetBinError   (j, resol->GetParameter(1)/TMath::Sqrt(hSlice->GetEntries()));
      }
      if(k==1){
	resolHeavyBin1->SetBinContent(j, resol->GetParameter(2));
	resolHeavyBin1->SetBinError  (j, resol->GetParameter(2)/TMath::Sqrt(hSlice->GetEntries()));
	respHeavyBin1->SetBinContent (j, resol->GetParameter(1));
	respHeavyBin1->SetBinError   (j, resol->GetParameter(1)/TMath::Sqrt(hSlice->GetEntries()));
      }

      delete resol;
    }


    if(k==0)
      accHeavyBin0->Divide( accPass,accTot, 1.0,1.0);
    if(k==1)
      accHeavyBin1->Divide( accPass,accTot, 1.0,1.0);
      
    delete hCorrLight;  delete hCorrHeavy;
  }

  TF1* turnOnLightBin0 = new TF1("turnOnBin0","TMath::Erf([0]*x+[1])+[2]", 0,400);
  TF1* turnOnLightBin1 = new TF1("turnOnBin1","TMath::Erf([0]*x+[1])+[2]", 0,400);
  TF1* turnOnHeavyBin0 = new TF1("turnOnBin0","TMath::Erf([0]*x+[1])+[2]", 0,400);
  TF1* turnOnHeavyBin1 = new TF1("turnOnBin1","TMath::Erf([0]*x+[1])+[2]", 0,400);
  accLightBin0->Fit(turnOnLightBin0,"","",25,400);  
  accLightBin1->Fit(turnOnLightBin1,"","",25,400);
  accHeavyBin0->Fit(turnOnHeavyBin0,"","",25,400);  
  accHeavyBin1->Fit(turnOnHeavyBin1,"","",25,400);

  TF1* sigmaLightBin0 = new TF1("sigmaLightBin0","x*sqrt([0]*[0]/x + [1]*[1]/x/x)",40,200);
  TF1* sigmaLightBin1 = new TF1("sigmaLightBin1","x*sqrt([0]*[0]/x + [1]*[1]/x/x)",40,200);
  TF1* sigmaHeavyBin0 = new TF1("sigmaHeavyBin0","x*sqrt([0]*[0]/x + [1]*[1]/x/x)",40,200);
  TF1* sigmaHeavyBin1 = new TF1("sigmaHeavyBin1","x*sqrt([0]*[0]/x + [1]*[1]/x/x)",40,200);
  TF1* meanLightBin0  = new TF1("meanLightBin0", "x*[0] + [1]",40,200);
  TF1* meanLightBin1  = new TF1("meanLightBin1", "x*[0] + [1]",40,200);
  TF1* meanHeavyBin0  = new TF1("meanHeavyBin0", "x*[0] + [1]",40,200);
  TF1* meanHeavyBin1  = new TF1("meanHeavyBin1", "x*[0] + [1]",40,200);

  resolLightBin0->Fit(sigmaLightBin0,"","",40,160);
  resolLightBin1->Fit(sigmaLightBin1,"","",40,160);
  resolHeavyBin0->Fit(sigmaHeavyBin0,"","",40,160);
  resolHeavyBin1->Fit(sigmaHeavyBin1,"","",40,160);
  respLightBin0->Fit (meanLightBin0, "","",40,160);
  respLightBin1->Fit (meanLightBin1, "","",40,160);
  respHeavyBin0->Fit (meanHeavyBin0, "","",40,160);
  respHeavyBin1->Fit (meanHeavyBin1, "","",40,160);

  gSystem->cd("./root/");
  RooWorkspace *w = new RooWorkspace("transferFuntions","all functions for SL ttH");

  RooRealVar param0AccLightBin0("param0AccLightBin0","", turnOnLightBin0->GetParameter(0));
  RooRealVar param1AccLightBin0("param1AccLightBin0","", turnOnLightBin0->GetParameter(1));
  RooRealVar param2AccLightBin0("param2AccLightBin0","", turnOnLightBin0->GetParameter(2));
  RooRealVar param0AccLightBin1("param0AccLightBin1","", turnOnLightBin1->GetParameter(0));
  RooRealVar param1AccLightBin1("param1AccLightBin1","", turnOnLightBin1->GetParameter(1));
  RooRealVar param2AccLightBin1("param2AccLightBin1","", turnOnLightBin1->GetParameter(2));
  RooRealVar param0AccHeavyBin0("param0AccHeavyBin0","", turnOnHeavyBin0->GetParameter(0));
  RooRealVar param1AccHeavyBin0("param1AccHeavyBin0","", turnOnHeavyBin0->GetParameter(1));
  RooRealVar param2AccHeavyBin0("param2AccHeavyBin0","", turnOnHeavyBin0->GetParameter(2));
  RooRealVar param0AccHeavyBin1("param0AccHeavyBin1","", turnOnHeavyBin1->GetParameter(0));
  RooRealVar param1AccHeavyBin1("param1AccHeavyBin1","", turnOnHeavyBin1->GetParameter(1));
  RooRealVar param2AccHeavyBin1("param2AccHeavyBin1","", turnOnHeavyBin1->GetParameter(2));

  RooRealVar param0resolLightBin0("param0resolLightBin0","", sigmaLightBin0->GetParameter(0));
  RooRealVar param1resolLightBin0("param1resolLightBin0","", sigmaLightBin0->GetParameter(1));
  RooRealVar param0resolLightBin1("param0resolLightBin1","", sigmaLightBin1->GetParameter(0));
  RooRealVar param1resolLightBin1("param1resolLightBin1","", sigmaLightBin1->GetParameter(1));
  RooRealVar param0respLightBin0 ("param0respLightBin0","",  meanLightBin0->GetParameter(0));
  RooRealVar param1respLightBin0 ("param1respLightBin0","",  meanLightBin0->GetParameter(1));
  RooRealVar param0respLightBin1 ("param0respLightBin1","",  meanLightBin1->GetParameter(0));
  RooRealVar param1respLightBin1 ("param1respLightBin1","",  meanLightBin1->GetParameter(1));

  RooRealVar param0resolHeavyBin0("param0resolHeavyBin0","", sigmaHeavyBin0->GetParameter(0));
  RooRealVar param1resolHeavyBin0("param1resolHeavyBin0","", sigmaHeavyBin0->GetParameter(1));
  RooRealVar param0resolHeavyBin1("param0resolHeavyBin1","", sigmaHeavyBin1->GetParameter(0));
  RooRealVar param1resolHeavyBin1("param1resolHeavyBin1","", sigmaHeavyBin1->GetParameter(1));
  RooRealVar param0respHeavyBin0 ("param0respHeavyBin0","",  meanHeavyBin0->GetParameter(0));
  RooRealVar param1respHeavyBin0 ("param1respHeavyBin0","",  meanHeavyBin0->GetParameter(1));
  RooRealVar param0respHeavyBin1 ("param0respHeavyBin1","",  meanHeavyBin1->GetParameter(0));
  RooRealVar param1respHeavyBin1 ("param1respHeavyBin1","",  meanHeavyBin1->GetParameter(1));


  w->import( param0AccLightBin0 ) ;
  w->import( param1AccLightBin0 ) ;
  w->import( param2AccLightBin0 ) ;
  w->import( param0AccLightBin1 ) ;
  w->import( param1AccLightBin1 ) ;
  w->import( param2AccLightBin1 ) ;
  w->import( param0AccHeavyBin0 ) ;
  w->import( param1AccHeavyBin0 ) ;
  w->import( param2AccHeavyBin0 ) ;
  w->import( param0AccHeavyBin1 ) ;
  w->import( param1AccHeavyBin1 ) ;
  w->import( param2AccHeavyBin1 ) ;

  w->import(  param0resolLightBin0 );
  w->import(  param1resolLightBin0 );
  w->import(  param0resolLightBin1 );
  w->import(  param1resolLightBin1 );
  w->import(  param0respLightBin0 );
  w->import(  param1respLightBin0 );
  w->import(  param0respLightBin1 );
  w->import(  param1respLightBin1 );
  w->import(  param0resolHeavyBin0 );
  w->import(  param1resolHeavyBin0 );
  w->import(  param0resolHeavyBin1 );
  w->import(  param1resolHeavyBin1 );
  w->import(  param0respHeavyBin0 );
  w->import(  param1respHeavyBin0 );
  w->import(  param0respHeavyBin1 );
  w->import(  param1respHeavyBin1 );

  //RooRealVar X1       ("X1","sqrt(X1*X2)", 0.0, 0.5);
  //RooRealVar X2       ("X2","X1-X2",      -0.6, 0.6);
  //RooRealVar PtTTH    ("PtTTH", "PtTTH",     0, 600);
  //RooRealVar PhiTTH   ("PhiTTH","PhiTTH",     -TMath::Pi(),  TMath::Pi() );
  //RooRealVar PtHStar  ("PtHStar","PtHStar",0,0.50);  
  //RooRealVar DphiHStar("DphiHStar","DphiHStar",-1,1);  
  RooRealVar BetaWHad ("BetaW", "BetaWHad",    -1,1);  
  RooRealVar GammaWHad("GammaW","GammaWHad",   -1,1);  
  RooRealVar BetaWLep ("BetaWLep", "BetaWLep", -1,1);  
  RooRealVar GammaWLep("GammaWLep","GammaWLep",-1,1);  
  RooRealVar MassTT   ("MassTT",   "MassTT",  340,2000);
  RooRealVar GammaTT  ("GammaTT",  "GammaTT",  -1,1);
  RooRealVar GammaTTH ("GammaTTH", "GammaTTH", -1,1);  
  RooRealVar X1       ("X1",           "X1",   0.04, 0.3);

  //RooRealVar PxT      ("PxT","PxT",-600,600);  
  //RooRealVar PyT      ("PyT","PyT",-600,600);  
  //RooRealVar PzT      ("PzT","PzT",-600,600);
  //RooRealVar PxH      ("PxH","PxH",-600,600);  
  //RooRealVar PyH      ("PyH","PyH",-600,600);  
  //RooRealVar PzH      ("PzH","PzH",-600,600);  

  RooDataSet dataset("dataset","dataset", RooArgSet(BetaWHad,GammaWHad,BetaWLep,GammaWLep,MassTT,GammaTT,GammaTTH,X1), Import( *tree ) );

  // v1)
  RooDataHist histGammaWHad("histGammaWHad","histGammaWHad", RooArgSet(GammaWHad), dataset, 1.0);
  RooRealVar p0v1("p0v1","",  1,  0,100);
  RooRealVar p1v1("p1v1","",  1,  0,100);
  RooRealVar p2v1("p2v1","",  1,  0,100);
  RooRealVar p3v1("p3v1","",  1,  0,100);
  RooBernstein pdfGammaWHad( "pdfGammaWHad", "pdfGammaWHad", GammaWHad, RooArgList(p0v1,p1v1,p2v1,p3v1));
  pdfGammaWHad.fitTo(histGammaWHad,  Minos(1), Save(1), Range(-1.,1.) ); 
  w->import(pdfGammaWHad);


  // v2)
  RooDataHist histBetaWHad("histBetaWHad","histBetaWHad", RooArgSet(BetaWHad), dataset, 1.0);
  RooRealVar p0v2("p0v2","",  1,  0,100);
  RooRealVar p1v2("p1v2","",  1,  0,100);
  RooBernstein pdfBetaWHad( "pdfBetaWHad", "pdfBetaWHad", BetaWHad, RooArgList(p0v2,p1v2));
  pdfBetaWHad.fitTo(histBetaWHad,  Minos(1), Save(1), Range(-1.,1.) ); 
  w->import(pdfBetaWHad);

  // v3)
  RooDataHist histGammaWLep("histGammaWLep","histGammaWLep", RooArgSet(GammaWLep), dataset, 1.0);
  RooRealVar p0v3("p0v3","",  1,  0,100);
  RooRealVar p1v3("p1v3","",  1,  0,100);
  RooRealVar p2v3("p2v3","",  1,  0,100);
  RooRealVar p2v4("p2v4","",  1,  0,100);
  RooBernstein pdfGammaWLep( "pdfGammaWLep", "pdfGammaWLep", GammaWLep, RooArgList(p0v3,p1v3,p2v3,p2v4));
  pdfGammaWLep.fitTo(histGammaWLep,  Minos(1), Save(1), Range(-1.,1.) ); 
  w->import(pdfGammaWLep);


  // v4)
  RooDataHist histBetaWLep("histBetaWLep","histBetaWLep", RooArgSet(BetaWLep), dataset, 1.0);
  RooRealVar p0v4("p0v4","",  1,  0,100);
  RooRealVar p1v4("p1v4","",  1,  0,100);
  RooBernstein pdfBetaWLep( "pdfBetaWLep", "pdfBetaWLep", BetaWLep, RooArgList(p0v4,p1v4));
  pdfBetaWLep.fitTo(histBetaWLep,  Minos(1), Save(1), Range(-1.,1.) ); 
  w->import(pdfBetaWLep);

  // v5)
  RooDataHist histGammaTTH("histGammaTTH","histGammaTTH", RooArgSet(GammaTTH), dataset, 1.0);
  RooRealVar p0v5("p0v5","",  1,  0,100);
  RooRealVar p1v5("p1v5","",  1,  0,100);
  RooRealVar p2v5("p2v5","",  1,  0,100);
  RooBernstein pdfGammaTTH( "pdfGammaTTH", "pdfGammaTTH", GammaTTH, RooArgList(p0v5,p1v5,p2v5));
  pdfGammaTTH.fitTo(histGammaTTH,  Minos(1), Save(1), Range(-1.,1.) ); 
  w->import(pdfGammaTTH);

  // v6)
  //GammaTT.setBins(  3 );
  //MassTT.setBins ( 40 );
  //X1.setBins     ( 40 );

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

  TH1F* hBinCont = new TH1F("hBinCont","",200,0,200);
  TH3F* h3 = new TH3F("h3","", nBinsX, binsX.GetArray(), nBinsMassTT, binsMassTT.GetArray(), nBinsGammaTT, binsGammaTT.GetArray());
  tree->Draw("GammaTT:MassTT:X1>>h3");
  
  RooDataHist hist3D("hist3D","hist3D", RooArgSet(X1,MassTT,GammaTT), Import( *h3, 1 ));
  RooHistPdf pdf3D  ("pdf3D","", RooArgSet(X1,MassTT,GammaTT), hist3D);
  w->import( pdf3D );

  for( int i = 1; i <= h3->GetNbinsX(); i++){
    for( int j = 1; j <= h3->GetNbinsY(); j++){
      for( int k = 1; k <= h3->GetNbinsZ(); k++){
	int bin = h3->GetBin(i,j,k);
	float binC  = h3->GetBinContent( bin );
	float width = h3->GetBinWidth( bin );
	h3->SetBinContent(bin, binC/width );
      }
    }
  }
  
  for(int i = 0 ; i < 1000 ; i++){
    double var1,var2,var3;
    h3->GetRandom3(var1,var2,var3);
    hBinCont->Fill(  h3->GetBinContent(h3->FindBin(var1,var2,var3)) );
  }



  ////////////////////////////////////////////////////////////

  int nBinsEt = 16;
  TArrayF binsEt(nBinsEt+1);
  cout << "Making histograms with " << nBinsEt << " bins:" << endl;
  for(int k = 0; k < 11; k++)
    binsEt[k] = 10*k;
  for(int k = 11; k < 16; k++)
    binsEt[k] = 100 + 50*(k-11);
  binsEt[16] = 1000;

  int nBinsPullEt = 25;
  TArrayF binsPullEt(nBinsPullEt+1);
  cout << "Making histograms with " << nBinsPullEt << " bins:" << endl;
  for(int k = 0; k < 21; k++)
    binsPullEt[k] = -1+k*(2./20);
  binsPullEt[21] = 1.5;
  binsPullEt[22] = 2.0;
  binsPullEt[23] = 3.0;
  binsPullEt[24] = 4.0;
  binsPullEt[25] = 10 ;

  TH2F* h2_0 = new TH2F("h2_0","", nBinsEt, binsEt.GetArray(), nBinsPullEt, binsPullEt.GetArray());
  treeEvent->Draw("(etReco-et)/et:et>>h2_0","sumEt<1200");

  TH2F* h2_1 = new TH2F("h2_1","", nBinsEt, binsEt.GetArray(), nBinsPullEt, binsPullEt.GetArray());
  treeEvent->Draw("(etReco-et)/et:et>>h2_1","sumEt>1200 && sumEt<1800");

  TH2F* h2_2 = new TH2F("h2_2","", nBinsEt, binsEt.GetArray(), nBinsPullEt, binsPullEt.GetArray());
  treeEvent->Draw("(etReco-et)/et:et>>h2_2","sumEt>1800");


  for(int b = 1; b <= h2_0->GetNbinsX(); b++){
    int counter = 0.;
    for(int bb = 1; bb <= h2_0->GetNbinsY(); bb++){
      counter += h2_0->GetBinContent(b,bb);
    }
    for(int bb = 1; bb <= h2_0->GetNbinsY(); bb++){
      int bin = h2_0->GetBin(b,bb);
      float width = h2_0->GetYaxis()->GetBinWidth( bb );
      float newbc = h2_0->GetBinContent( bin );
      if(counter>0)
	h2_0->SetBinContent(bin, newbc/width/counter );
    }
  }

  for(int b = 1; b <= h2_1->GetNbinsX(); b++){
    int counter = 0.;
    for(int bb = 1; bb <= h2_1->GetNbinsY(); bb++){
      counter += h2_1->GetBinContent(b,bb);
    }
    for(int bb = 1; bb <= h2_1->GetNbinsY(); bb++){
      int bin = h2_1->GetBin(b,bb);
      float width = h2_1->GetYaxis()->GetBinWidth( bb );
      float newbc = h2_1->GetBinContent( bin );
      if(counter>0)
	h2_1->SetBinContent(bin, newbc/width/counter );
    }
  }

  for(int b = 1; b <= h2_2->GetNbinsX(); b++){
    int counter = 0.;
    for(int bb = 1; bb <= h2_2->GetNbinsY(); bb++){
      counter += h2_2->GetBinContent(b,bb);
    }
    for(int bb = 1; bb <= h2_2->GetNbinsY(); bb++){
      int bin = h2_2->GetBin(b,bb);
      float width = h2_2->GetYaxis()->GetBinWidth( bb );
      float newbc = h2_2->GetBinContent( bin );
      if(counter>0)
	h2_2->SetBinContent(bin, newbc/width/counter );
    }
  }


  // v5)
  //RooDataHist histMassTT("histMassTT","histMassTT", RooArgSet(MassTT), dataset, 1.0);
  //RooRealVar p0v5("p0v5","p0v5",400,300,600);
  //RooRealVar p1v5("p1v5","p1v5",100,20, 500);
  //RooRealVar p2v5("p2v5","p2v5",-1,0,-20);
  //RooRealVar p3v5("p3v5","p3v5",1,0,50);
  //RooCBShape pdfMassTT("pdfMassTT","",MassTT,p0v5,p1v5,p2v5,p3v5);
  //pdfMassTT.fitTo(histMassTT,  Minos(1), Save(1) , Range(360.,1400.) ); 


  /*
  RooRealVar csvRecoLight("csvReco","csv", -0.1, 1.1);
  RooDataSet datasetBtagLight("datasetBtagLight","dataset", RooArgSet(csvRecoLight), Import( *treeJetsLight ) );
  //RooKeysPdf pdfCsvLight("pdfCsvLight","", csvRecoLight, datasetBtagLight);
  //w->import( pdfCsvLight );

  RooRealVar csvRecoHeavy("csvReco","csv", -0.1, 1.1);
  RooDataSet datasetBtagHeavy("datasetBtagHeavy","dataset", RooArgSet(csvRecoHeavy), Import( *treeJetsHeavy ) );
  //RooKeysPdf pdfCsvHeavy("pdfCsvHeavy","", csvRecoHeavy, datasetBtagHeavy);
  //w->import( pdfCsvHeavy );

  RooRealVar v1       ("X1","sqrt(X1*X2)", 0.0, 0.5);
  RooRealVar v2       ("X2","X1-X2",      -0.6, 0.6);
  RooRealVar PtTTH    ("PtTTH", "PtTTH",     0, 600);
  RooRealVar PhiTTH   ("PhiTTH","PhiTTH",     -TMath::Pi(),  TMath::Pi() );
  RooRealVar PtHStar  ("PtHStar","PtHStar",0,0.50);  
  RooRealVar DphiHStar("DphiHStar","DphiHStar",-1,1);  
  RooRealVar BetaW    ("BetaW","BetaW",-1,1);  
  RooRealVar GammaW   ("GammaW","GammaW",-1,1);  
  RooRealVar DeltaW   ("DeltaW","DeltaW",-1,1);  

  RooRealVar PxT      ("PxT","PxT",-600,600);  
  RooRealVar PyT      ("PyT","PyT",-600,600);  
  RooRealVar PzT      ("PzT","PzT",-600,600);
  RooRealVar PxH      ("PxH","PxH",-600,600);  
  RooRealVar PyH      ("PyH","PyH",-600,600);  
  RooRealVar PzH      ("PzH","PzH",-600,600);  

  RooDataSet dataset("dataset","dataset", RooArgSet(v1,v2, PtTTH, PhiTTH,PtHStar,DphiHStar,BetaW,GammaW,DeltaW), Import( *tree ) );
  RooDataSet datasetCMS("datasetCMS","dataset", RooArgSet(PxH,PyH,PzH,PxT,PyT,PzT), Import( *tree ) );

  // v1)
  ///////////////////////////////////////////////////////////////////////////////
  RooDataHist histV1("histV1","histV1", RooArgSet(v1), dataset, 1.0);

  RooRealVar p0v1("p0v1","",0.1,0,1);
  RooRealVar p1v1("p1v1","",0.2,0,1);
  RooLandau  pdfV1("pdfV1","pdfV1", v1 ,p0v1,p1v1);

  pdfV1.fitTo(histV1,  Minos(1), Save(1), Range(0,0.5) ); 
  w->import(pdfV1);
  ///////////////////////////////////////////////////////////////////////////////

  // v2)
  ///////////////////////////////////////////////////////////////////////////////
  RooDataHist histV2("histV2","histV2", RooArgSet(v2), dataset, 1.0);

  RooRealVar p0v2("p0v2","",0.,-1,1);
  RooRealVar p1v2("p1v2","",0.2,-2,2);
  RooGaussian pdfV2("pdfV2","pdfV2", v2 ,p0v2,p1v2);

  pdfV2.fitTo(histV2,  Minos(1), Save(1), Range(-0.6,0.6) ); 
  w->import(pdfV2);
  ///////////////////////////////////////////////////////////////////////////////


  // v3)
  ///////////////////////////////////////////////////////////////////////////////
  RooDataHist histPtTTH("histPtTTH","histPtTTH", RooArgSet(PtTTH), dataset, 1.0);
  
  RooRealVar p0v3("p0v3","",0.,0,100);
  RooRealVar p1v3("p1v3","",0.,-100,100);

  RooGenericPdf pdfPtTTH("pdfPtTTH","@0^(@1)*exp(@2*@0)", RooArgList(PtTTH,p0v3,p1v3 ));

  //pdfPtTTH.fitTo(histPtTTH,  Minos(1), Save(1), Range(0.,500.) ); 
  //RooFitResult*  fitResults = 
  //RooArgSet fitParam(fitResults->floatParsFinal());
  //RooRealVar* p0v3Fit = (RooRealVar*)(&fitParam["p0v3"]);
  //RooRealVar* p1v3Fit = (RooRealVar*)(&fitParam["p1v3"]);

  pdfPtTTH.fitTo(histPtTTH,  Minos(1), Save(1), Range(0.,500.) ); 
  w->import(pdfPtTTH);
  ///////////////////////////////////////////////////////////////////////////////


  // v4)
  ///////////////////////////////////////////////////////////////////////////////
  RooDataHist histPtHStar("histPtHStar","histPtHStar", RooArgSet(PtHStar), dataset, 1.0);

  RooRealVar p0v4("p0v4","",0.2,0,1);
  RooRealVar p1v4("p1v4","",0.2,-2,2);
  RooGaussian pdfPtHStar("pdfPtHStar","pdfPtHStar",  PtHStar,p0v4,p1v4);

  pdfPtHStar.fitTo(histPtHStar,  Minos(1), Save(1), Range(0.,0.5) ); 
  w->import(pdfPtHStar);
  ///////////////////////////////////////////////////////////////////////////////


  // v5)
  ///////////////////////////////////////////////////////////////////////////////
  RooDataHist histDphiHStar("histDphiHStar","histDphiHStar", RooArgSet(DphiHStar), dataset, 1.0);
  
  RooRealVar p0v5("p0v5","",0.,  0,1000);
  RooRealVar p1v5("p1v5","",0.,  0.,10.);
  RooRealVar p2v5("p2v5","",0.,  0,1000);

  TH1F* param0DphiHStar = new TH1F("param0DphiHStar","",3,0.05,0.35);
  TH1F* param1DphiHStar = new TH1F("param1DphiHStar","",3,0.05,0.35);
  TH1F* param2DphiHStar = new TH1F("param2DphiHStar","",3,0.05,0.35);
  TH1F* ratioDphiHStar  = new TH1F("ratioDphiHStar","", 3,0.05,0.35);


  float bounds[4] = {0.05, 0.15, 0.25, 0.35};
  for(int k = 0; k < 3 ; k++){
    
    RooDataSet dataTmp("dataTmp","", RooArgSet(DphiHStar,PtHStar), Import( *tree ), 
		       Cut(Form("PtHStar>=%f && PtHStar<%f", bounds[k], bounds[k+1])) );
    RooRealVar p0v5Tmp("p0v5Tmp","",0.,  0,1000);
    RooRealVar p1v5Tmp("p1v5Tmp","",0.,  0., 10.);
    RooRealVar p2v5Tmp("p2v5Tmp","",0.,  0,1000);
    RooGenericPdf pdfTmp("pdfTmp","@1*exp(-(@2*@0)) + @3", RooArgList(DphiHStar,p0v5Tmp,p1v5Tmp,p2v5Tmp ));
    RooFitResult*  fitResults  = pdfTmp.fitTo(dataTmp,  Minos(1), Save(1), Range(-1.,1.) ); 
    RooArgSet fitParam(fitResults->floatParsFinal());
    RooRealVar* p0v5Fit = (RooRealVar*)(&fitParam["p0v5Tmp"]);
    RooRealVar* p1v5Fit = (RooRealVar*)(&fitParam["p1v5Tmp"]);
    RooRealVar* p2v5Fit = (RooRealVar*)(&fitParam["p2v5Tmp"]);

    param0DphiHStar->SetBinContent( k+1, p0v5Fit->getVal()   );
    param1DphiHStar->SetBinContent( k+1, p1v5Fit->getVal()   );
    param1DphiHStar->SetBinError  ( k+1, p1v5Fit->getError()   );
    param2DphiHStar->SetBinContent( k+1, p2v5Fit->getVal()   );
    ratioDphiHStar->SetBinContent(  k+1, p0v5Fit->getVal()/(p0v5Fit->getVal()+ p2v5Fit->getVal()) );
    //param0DphiHStar->SetBinError  ( k+1, p0v5Fit->getError() );

  }

  TF1* pol1 = new TF1("pol1","[0] + [1]*x",0.,0.5);
  TF1* pol2 = new TF1("pol2","[0] + [1]*x",0.,0.5);
  ratioDphiHStar->Fit(pol1);
  param1DphiHStar->Fit(pol2);
  RooRealVar interceptRatiopol1DphiHStar("interceptRatiopol1DphiHStar", "", pol1->GetParameter(0));
  RooRealVar slopeRatiopol1DphiHStar    ("slopeRatiopol1DphiHStar", ""    , pol1->GetParameter(1));
  RooRealVar interceptTaupol1DphiHStar  ("interceptTaupol1DphiHStar", "",   pol2->GetParameter(0));
  RooRealVar slopeTaupol1DphiHStar      ("slopeTaupol1DphiHStar", ""    ,   pol2->GetParameter(1));
  w->import( interceptRatiopol1DphiHStar );
  w->import( slopeRatiopol1DphiHStar );
  w->import( interceptTaupol1DphiHStar );
  w->import( slopeTaupol1DphiHStar );

  RooGenericPdf pdfDphiHStar("pdfDphiHStar","@1*exp(-(@2*@0)) + @3", RooArgList(DphiHStar,p0v5,p1v5,p2v5 ));
  pdfDphiHStar.fitTo(histDphiHStar,  Minos(1), Save(1), Range(-1.,1.) ); 
  w->import(pdfDphiHStar);
  ///////////////////////////////////////////////////////////////////////////////


  // v1)
  ///////////////////////////////////////////////////////////////////////////////
  //RooNDKeysPdf pdf2DStar("histBetaW","histBetaW", RooArgList(PtHStar,DphiHStar), dataset);
  //w->import(pdf2DStar);
  ///////////////////////////////////////////////////////////////////////////////

  
  // v7)
  ///////////////////////////////////////////////////////////////////////////////
  RooDataHist histBetaW("histBetaW","histBetaW", RooArgSet(BetaW), dataset, 1.0);
  //RooHistPdf  pdfBetaW  ("pdfBetaW", "pdfBetaW",  RooArgSet(BetaW), dataset);

  RooRealVar p0v7("p0v7","",  1,  0,100);
  RooRealVar p1v7("p1v7","",  1,  0,100);
  RooRealVar p2v7("p2v7","",  1,  0,100);
  RooRealVar p3v7("p3v7","",  1,  0,100);

  RooBernstein pdfBetaW( "pdfBetaW", "pdfBetaW", BetaW, RooArgList(p0v7,p1v7,p2v7,p3v7));
  pdfBetaW.fitTo(histBetaW,  Minos(1), Save(1), Range(-1,1) ); 
  w->import(pdfBetaW);
  ///////////////////////////////////////////////////////////////////////////////


  // v8)
  ///////////////////////////////////////////////////////////////////////////////
  RooDataHist histGammaW("histGammaW","histGammaW", RooArgSet(GammaW), dataset, 1.0);
  //RooHistPdf  pdfGammaW  ("pdfGammaW", "pdfGammaW",  RooArgSet(GammaW), dataset);

  RooRealVar p0v8("p0v8","",  1,  0,100);
  RooRealVar p1v8("p1v8","",  1,  0,100);
  RooRealVar p2v8("p2v8","",  1,  0,100);
  RooRealVar p3v8("p3v8","",  1,  0,100);

  RooBernstein pdfGammaW( "pdfGammaW", "pdfGammaW", GammaW, RooArgList(p0v8,p1v8,p2v8,p3v8));
  pdfGammaW.fitTo(histGammaW,  Minos(1), Save(1), Range(-1,1) ); 
  w->import(pdfGammaW);
  ///////////////////////////////////////////////////////////////////////////////


  // v9)
  ///////////////////////////////////////////////////////////////////////////////
  RooNDKeysPdf pdfCMS("pdfCMS","pdfCMS", RooArgList(PxT,PyT,PzT,PxH,PyH,PzH), datasetCMS);

  PxT.setVal( 20. );
  PyT.setVal( -20. );
  PzT.setVal( 20. );
  PxH.setVal( -20. );
  PyH.setVal( 20. );
  PzH.setVal( -20. );
  
  cout << "Comb PDF: " <<  pdfCMS.getVal( RooArgSet(*(&PxT),*(&PyT),*(&PzT),*(&PxH),*(&PyH),*(&PzH)) ) << endl;

  PxT.setVal( 20. );
  PyT.setVal( -20. );
  PzT.setVal( 20. );
  PxH.setVal( -20. );
  PyH.setVal( 200. );
  PzH.setVal( -20. );
  
  cout << "Comb PDF: " <<  pdfCMS.getVal( RooArgSet(*(&PxT),*(&PyT),*(&PzT),*(&PxH),*(&PyH),*(&PzH)) ) << endl;

  PxT.setVal( 20. );
  PyT.setVal( -20. );
  PzT.setVal( 20. );
  PxH.setVal( -20. );
  PyH.setVal( 500. );
  PzH.setVal( -20. );
  
  cout << "Comb PDF: " <<  pdfCMS.getVal( RooArgSet(*(&PxT),*(&PyT),*(&PzT),*(&PxH),*(&PyH),*(&PzH)) ) << endl;

  //w->import(pdf2DStar);
  ///////////////////////////////////////////////////////////////////////////////

  */

  w->writeToFile(outFileName.c_str(),kTRUE);

  /*
  TCanvas *c1CsvLight = new TCanvas("c1CsvLight","canvas",10,30,650,600);
  RooPlot* plotCsvLight = csvRecoLight.frame(Bins(20),Title("cvs light"));
  datasetBtagLight.plotOn(plotCsvLight);
  //pdfCsvLight.plotOn(plotCsvLight);

  TCanvas *c1CsvHeavy = new TCanvas("c1CsvHeavy","canvas",10,30,650,600);
  RooPlot* plotCsvHeavy = csvRecoHeavy.frame(Bins(20),Title("cvs heavy"));
  datasetBtagHeavy.plotOn(plotCsvHeavy);
  //pdfCsvHeavy.plotOn(plotCsvHeavy);

  TCanvas *c1V1 = new TCanvas("c1V1","canvas",10,30,650,600);
  RooPlot* plotV1 = v1.frame(Bins(20),Title("v1"));
  histV1.plotOn(plotV1);
  pdfV1.plotOn(plotV1);

  TCanvas *c1V2 = new TCanvas("c1V2","canvas",10,30,650,600);
  RooPlot* plotV2 = v2.frame(Bins(20),Title("v2"));
  histV2.plotOn(plotV2);
  pdfV2.plotOn(plotV2);

  TCanvas *c1PtTTH = new TCanvas("c1PtTTH","canvas",10,30,650,600);
  RooPlot* plotPtTTH = PtTTH.frame(Bins(20),Title("PtTTH"));
  histPtTTH.plotOn(plotPtTTH);
  pdfPtTTH.plotOn(plotPtTTH);

  TCanvas *c1PtHStar = new TCanvas("c1PtHStar","canvas",10,30,650,600);
  RooPlot* plotPtHStar = PtHStar.frame(Bins(20),Title("PtHStar"));
  histPtHStar.plotOn(plotPtHStar);
  pdfPtHStar.plotOn(plotPtHStar);

  TCanvas *c1DphiHStar = new TCanvas("c1DphiHStar","canvas",10,30,650,600);
  RooPlot* plotDphiHStar = DphiHStar.frame(Bins(20),Title("DphiHStar"));
  histDphiHStar.plotOn(plotDphiHStar);
  pdfDphiHStar.plotOn(plotDphiHStar);

  TCanvas *c1BetaW = new TCanvas("c1BetaW","canvas",10,30,650,600);
  RooPlot* plotBetaW = BetaW.frame(Bins(20),Title("BetaW"));
  histBetaW.plotOn(plotBetaW);
  pdfBetaW.plotOn(plotBetaW);

  TCanvas *c1GammaW = new TCanvas("c1GammaW","canvas",10,30,650,600);
  RooPlot* plotGammaW = GammaW.frame(Bins(20),Title("GammaW"));
  histGammaW.plotOn(plotGammaW);
  pdfGammaW.plotOn(plotGammaW);
  */

  TCanvas *c1BetaWHad = new TCanvas("c1BetaWHad","canvas",10,30,650,600);
  RooPlot* plotBetaWHad = BetaWHad.frame(Bins(20),Title("BetaWHad"));
  histBetaWHad.plotOn(plotBetaWHad);
  pdfBetaWHad.plotOn(plotBetaWHad);

  TCanvas *c1GammaWHad = new TCanvas("c1GammaWHad","canvas",10,30,650,600);
  RooPlot* plotGammaWHad = GammaWHad.frame(Bins(20),Title("GammaWHad"));
  histGammaWHad.plotOn(plotGammaWHad);
  pdfGammaWHad.plotOn(plotGammaWHad);

  TCanvas *c1BetaWLep = new TCanvas("c1BetaWLep","canvas",10,30,650,600);
  RooPlot* plotBetaWLep = BetaWLep.frame(Bins(20),Title("BetaWLep"));
  histBetaWLep.plotOn(plotBetaWLep);
  pdfBetaWLep.plotOn(plotBetaWLep);

  TCanvas *c1GammaWLep = new TCanvas("c1GammaWLep","canvas",10,30,650,600);
  RooPlot* plotGammaWLep = GammaWLep.frame(Bins(20),Title("GammaWLep"));
  histGammaWLep.plotOn(plotGammaWLep);
  pdfGammaWLep.plotOn(plotGammaWLep);

  TCanvas *c1GammaTTH = new TCanvas("c1GammaTTH","canvas",10,30,650,600);
  RooPlot* plotGammaTTH = GammaTTH.frame(Bins(20),Title("GammaTTH"));
  histGammaTTH.plotOn(plotGammaTTH);
  pdfGammaTTH.plotOn(plotGammaTTH);

  //TCanvas *c1MassTT = new TCanvas("c1MassTT","canvas",10,30,650,600);
  //RooPlot* plotMassTT = MassTT.frame(Bins(20),Title("MassTT"));
  //histMassTT.plotOn(plotMassTT);
  //pdfMassTT.plotOn(plotMassTT);


  fout = new TFile("ControlPlots.root","RECREATE");
  fout->cd();

  h2_0->Write("",TObject::kOverwrite);
  h2_1->Write("",TObject::kOverwrite);
  h2_2->Write("",TObject::kOverwrite);

  h3->Write("",TObject::kOverwrite);
  hBinCont->Write("",TObject::kOverwrite);
  
  c1BetaWHad->cd();
  plotBetaWHad->Draw();
  c1BetaWHad->Write("",TObject::kOverwrite);

  c1GammaWHad->cd();
  plotGammaWHad->Draw();
  c1GammaWHad->Write("",TObject::kOverwrite);

  c1BetaWLep->cd();
  plotBetaWLep->Draw();
  c1BetaWLep->Write("",TObject::kOverwrite);

  c1GammaWLep->cd();
  plotGammaWLep->Draw();
  c1GammaWLep->Write("",TObject::kOverwrite);

  c1GammaTTH->cd();
  plotGammaTTH->Draw();
  c1GammaTTH->Write("",TObject::kOverwrite);

  //c1MassTT->cd();
  //plotMassTT->Draw();
  //c1MassTT->Write("",TObject::kOverwrite);

  /*
  c1CsvLight->cd();
  plotCsvLight->Draw();
  c1CsvLight->Write("",TObject::kOverwrite);

  c1CsvHeavy->cd();
  plotCsvHeavy->Draw();
  c1CsvHeavy->Write("",TObject::kOverwrite);

  c1V1->cd();
  plotV1->Draw();
  c1V1->Write("",TObject::kOverwrite);

  c1V2->cd();
  plotV2->Draw();
  c1V2->Write("",TObject::kOverwrite);

  c1PtTTH->cd();
  plotPtTTH->Draw();
  c1PtTTH->Write("",TObject::kOverwrite);

  c1PtHStar->cd();
  plotPtHStar->Draw();
  c1PtHStar->Write("",TObject::kOverwrite);

  c1DphiHStar->cd();
  plotDphiHStar->Draw();
  c1DphiHStar->Write("",TObject::kOverwrite);
  param0DphiHStar->Write("",TObject::kOverwrite);
  param1DphiHStar->Write("",TObject::kOverwrite);
  param2DphiHStar->Write("",TObject::kOverwrite);
  ratioDphiHStar->Write( "",TObject::kOverwrite);

  c1BetaW->cd();
  plotBetaW->Draw();
  c1BetaW->Write("",TObject::kOverwrite);

  c1GammaW->cd();
  plotGammaW->Draw();
  c1GammaW->Write("",TObject::kOverwrite);
  */


  accLightBin0->Write("",TObject::kOverwrite);
  accLightBin1->Write("",TObject::kOverwrite);
  accHeavyBin0->Write("",TObject::kOverwrite);
  accHeavyBin1->Write("",TObject::kOverwrite);
  resolLightBin0->Write("",TObject::kOverwrite);
  respLightBin0->Write("",TObject::kOverwrite);
  resolLightBin1->Write("",TObject::kOverwrite);
  respLightBin1->Write("",TObject::kOverwrite);
  resolHeavyBin0->Write("",TObject::kOverwrite);
  respHeavyBin0->Write("",TObject::kOverwrite);
  resolHeavyBin1->Write("",TObject::kOverwrite);
  respHeavyBin1->Write("",TObject::kOverwrite);

  fout->Close();
  delete fout;


  return 0;
}

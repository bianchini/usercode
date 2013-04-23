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

#define GENJETDR 0.3
#define VERBOSE  false

using namespace std;
using namespace RooFit;



int main(int argc, const char* argv[])
{

  std::cout << "TransferFunctions_SL" << std::endl;
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

  for(int k = 0; k < 2; k++){

    TH1F* accPass   = (TH1F*)accLightBin0->Clone(Form("accPass%d",k)); //accPass->Sumw2();
    TH1F* accTot    = (TH1F*)accLightBin0->Clone(Form("accTot%d",k));  //accTot->Sumw2();
    treeJetsLight->Draw(Form("pt>>accTot%d",k),  Form("(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",            etaBinning[k],etaBinning[k+1]) );
    treeJetsLight->Draw(Form("pt>>accPass%d",k), Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",etaBinning[k],etaBinning[k+1]) );
    accPass->Sumw2();
    accTot->Sumw2();

    if(k==0)
      accLightBin0->Divide( accPass,accTot, 1.0,1.0);
    if(k==1)
      accLightBin1->Divide( accPass,accTot, 1.0,1.0);

    accPass->Reset();
    accTot->Reset();

    treeJetsHeavy->Draw(Form("pt>>accTot%d",k),  Form("(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",            etaBinning[k],etaBinning[k+1]) );
    treeJetsHeavy->Draw(Form("pt>>accPass%d",k), Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",etaBinning[k],etaBinning[k+1]) );
    accPass->Sumw2();
    accTot->Sumw2();

    if(k==0)
      accHeavyBin0->Divide( accPass,accTot, 1.0,1.0);
    if(k==1)
      accHeavyBin1->Divide( accPass,accTot, 1.0,1.0);
      
  }

  TF1* turnOnLightBin0 = new TF1("turnOnBin0","TMath::Erf([0]*x+[1])+[2]", 0,400);
  TF1* turnOnLightBin1 = new TF1("turnOnBin1","TMath::Erf([0]*x+[1])+[2]", 0,400);
  TF1* turnOnHeavyBin0 = new TF1("turnOnBin0","TMath::Erf([0]*x+[1])+[2]", 0,400);
  TF1* turnOnHeavyBin1 = new TF1("turnOnBin1","TMath::Erf([0]*x+[1])+[2]", 0,400);
  accLightBin0->Fit(turnOnLightBin0);  
  accLightBin1->Fit(turnOnLightBin1);
  accHeavyBin0->Fit(turnOnHeavyBin0);  
  accHeavyBin1->Fit(turnOnHeavyBin1);




  fout = new TFile("ControlPlots.root","UPDATE");
  fout->cd();

  accLightBin0->Write();
  accLightBin1->Write();
  accHeavyBin0->Write();
  accHeavyBin1->Write();

  fout->Close();
  delete fout;

  return 0;
  




  RooWorkspace *w = new RooWorkspace("transferFuntions","all functions for SL ttH");

  

  RooRealVar v1       ("X1","sqrt(X1*X2)", 0.0, 0.5);
  RooRealVar v2       ("X2","X1-X2",      -0.6, 0.6);
  RooRealVar PtTTH    ("PtTTH", "PtTTH",     0, 600);
  RooRealVar PhiTTH   ("PhiTTH","PhiTTH",     -TMath::Pi(),  TMath::Pi() );
  RooRealVar PtHStar  ("PtHStar","PtHStar",0,0.50);  
  RooRealVar DphiHStar("DphiHStar","DphiHStar",-1,1);  
  RooRealVar BetaW    ("BetaW","BetaW",-1,1);  
  RooRealVar GammaW   ("GammaW","GammaW",-1,1);  
  RooRealVar DeltaW   ("DeltaW","DeltaW",-1,1);  


  RooDataSet dataset("dataset","dataset", RooArgSet(v1,v2, PtTTH, PhiTTH,PtHStar,DphiHStar,BetaW,GammaW,DeltaW), Import( *tree ) );

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





  //w->import(v1);
  //w->import(v2);
  //w->import(PtTTH);
  //w->import(PhiTTH);
  //w->import(PtHStar);
  //w->import(DphiHStar);
  //w->import(BetaW);
  //w->import(GammaW);
  //w->import(DeltaW);
  //w->import(dataset);

  w->writeToFile(outFileName.c_str(),kTRUE);

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


  
  fout = new TFile("ControlPlots.root","UPDATE");
  fout->cd();

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

  accLightBin0->Write();
  accLightBin1->Write();
  accHeavyBin0->Write();
  accHeavyBin1->Write();

  fout->Close();
  delete fout;


  return 0;
}

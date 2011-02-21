
/*
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooBreitWigner.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
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
#include "RooAbsCategory.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "RooTruthModel.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooKeysPdf.h"
#include "RooStepFunction.h"
#include "RooParametricStepFunction.h"
#include "RooProdPdf.h"
#include "RooChebychev.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooConstVar.h"
#include "RooIntegralMorph.h"
#include "RooNumIntConfig.h"
#include "RooLognormal.h"
#include "RooCurve.h"
#include "RooNLLVar.h"

#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TArrayD.h"
#include "TGraphAsymmErrors.h"

#include <vector>


using namespace std;
using namespace RooFit;


void test()
*/
{

  using namespace std;
  using namespace RooFit;

  TCanvas *c2 = new TCanvas("canvas","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  TH1F* h1 = new TH1F("h1","pull nSgn; (nSgn^{fit}-nSgn^{exp})/#DeltanSgn^{fit}",80,-4,4);
  TH1F* h2 = new TH1F("h2","pull nSgn constrained; (nSgn^{fit}-nSgn^{exp})/#DeltanSgn^{fit}",80,-4,4);

  RooRealVar mass("mass","mass",-5.,5.);
  mass.setBins( 20 );

  // exact pdf
  //RooGaussian signal("signal","exact signal",mass,RooConst(1.0),RooConst(0.75));
  RooCBShape signal("signal","exact signal",mass,RooConst(1.0),RooConst(0.75),RooConst(0.75),RooConst(1.5));
  RooExponential background("background","exact background",mass,RooConst(-0.25));

  int expectedSgn = 50;
  int expectedBkg = 50;


  // model pdf
  RooRealVar mean("mean","mean",1.0,-5,5);
  RooRealVar sigma("sigma","sigma",0.5,0,10);
  RooRealVar alpha("alpha","alpha",0.5,0,10);
  RooRealVar n("n","n",1.5,0,30);
  //RooGaussian signalModel("signalModel","signal model",mass,mean,sigma);
  RooCBShape signalModel("signalModel","exact signal",mass,mean,sigma,alpha,n);

  RooRealVar c("c","lifetime",0.,-10.,10);
  RooExponential backgroundModel("backgroundModel","background model",mass,c);

  RooRealVar nSgn("nSgn","number of signal events",100,0,10000);
  RooRealVar nBkg("nBkg","number of background events",100,0,10000);

  // contraints on the nuisance parameters
  RooLognormal sigmaConstraint("sigmaConstraint","constraint on the mean",sigma,RooConst(0.75),RooConst(1.5));
  RooLognormal alphaConstraint("alphaConstraint","constraint on alpha",alpha,RooConst(0.75),RooConst(1.5));
  RooLognormal nConstraint("nConstraint","constraint on n",n,RooConst(1.5),RooConst(1.5));
  RooLognormal nBkgConstraint("nBkgConstraint","constraint on the background yield",nBkg,RooConst(50.),RooConst(1.5));


  RooAddPdf model("model","signal+background",RooArgList(signalModel,backgroundModel),RooArgList(nSgn,nBkg));

  RooPlot* plot1 = mass.frame(Title("Example of mass distribution"));
  RooPlot* plot2 = mass.frame(Title("Example of mass distribution"));

  for(unsigned int iToy = 0; iToy < 250; iToy++){
    RooDataSet iSet("iSet","",mass);
    iSet.reset();

    iSet.append(*(signal.generate(mass,     expectedSgn, Extended() )));
    iSet.append(*(background.generate(mass, expectedBkg, Extended() )));

    RooDataHist iHisto("iHisto","",RooArgSet(mass),iSet, 1.0);

    //reset
    nSgn.setVal(100);
    nBkg.setVal(100);
    mean.setVal(1.0);
    sigma.setVal(0.5);
    alpha.setVal(0.5);
    n.setVal(1.5);
    c.setVal(0.);

    RooFitResult* fitRes = model.fitTo(iHisto, Extended(),  Minos(1), Save(1), NumCPU(4));

    if(iToy==1){
      iHisto.plotOn(plot1);
      model.plotOn(plot1);
      model.plotOn(plot1,Components("signalModel"),LineColor(kRed));
      model.plotOn(plot1,Components("backgroundModel"),LineColor(kGreen));
    }
    

    //reset
    nSgn.setVal(100);
    nBkg.setVal(100);
    mean.setVal(1.0);
    sigma.setVal(0.5);
    alpha.setVal(0.5);
    n.setVal(1.5);
    c.setVal(0.);

    RooFitResult* fitResConstrained = model.fitTo(iHisto, ExternalConstraints(RooArgSet(nBkgConstraint/*sigmaConstraint*/,alphaConstraint,nConstraint)), Extended(),  Minos(1), Save(1), NumCPU(4));

    RooArgSet fitParam(fitRes->floatParsFinal());
    RooArgSet fitParamConstrained(fitResConstrained->floatParsFinal());

    RooRealVar* nSgn_i = (RooRealVar*)(&fitParam["nSgn"]); 
    RooRealVar* nSgnConstrained_i = (RooRealVar*)(&fitParamConstrained["nSgn"]); 
    h1->Fill( (nSgn_i->getVal() - expectedSgn)/nSgn_i->getError() );    
    h2->Fill( (nSgnConstrained_i->getVal() - expectedSgn)/nSgnConstrained_i->getError() );  

    if(iToy==1){
      iHisto.plotOn(plot2);
      model.plotOn(plot2);
      model.plotOn(plot2,Components("signalModel"),LineColor(kRed));
      model.plotOn(plot2,Components("backgroundModel"),LineColor(kGreen));
    }

  }

  TF1* gaus1 = new TF1("gaus1","gaus",-2,2);
  TF1* gaus2 = new TF1("gaus2","gaus",-2,2);

  c2->Divide(2,2);
  c2->cd(1);
  h1->Sumw2();
  h1->Fit(gaus1);
  h1->Draw();
  c2->cd(2);
  h2->Sumw2();
  h2->Fit(gaus2);
  h2->Draw();
  c2->cd(3);
  plot1->Draw();
  c2->cd(4);
  plot2->Draw();

}

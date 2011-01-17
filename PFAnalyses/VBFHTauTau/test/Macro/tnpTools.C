#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
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
#include "RooAbsCategory.h"

#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"

#include <vector>

using namespace std;
using namespace RooFit;

pair< Double_t*,Double_t* > simFitToMC(const string tnp_ = "etoTauMargLooseNoCracks70", 
				       const string category_ = "tauAntiEMVA",
				       const string var_ = "abseta",
				       const string bin_ = "abseta < 1.5",
				       const float binCenter_ = 0.75,
				       const float binWidth_ = 0.75,
				       const float xLow_=60, 
				       const float xHigh_=120,
				       const float BKGfracInF_ = 0.0035,
				       bool doCombined_ = true,
				       bool SumW2_ = false,
				       bool verbose_ = true){
  

  // files
  TFile f("./test_all.root");
  TFile fbkg("./test_bkg.root");
  TFile fsgn("./test_sgn.root");

  //////////////////////////////////////////////////////////////
  TTree *fullTreeSgn = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());
  // Zee signal from MC:
  TH1F* hSall = new TH1F("hSall","",1,0,150);
  TH1F* hSPall = new TH1F("hSPall","",1,0,150);
  TH1F* hS = new TH1F("hS","",1,0,150);
  TH1F* hSP = new TH1F("hSP","",1,0,150);
  fullTreeSgn->Draw("mass>>hS",Form("weight*(%s && mass>%f && mass<%f && mcTrue)",bin_.c_str(),xLow_,xHigh_));
  fullTreeSgn->Draw("mass>>hSall",Form("weight*(%s && mass>%f && mass<%f)",bin_.c_str(),xLow_,xHigh_));
  float SGN = hS->Integral();
  float SGNall = hSall->Integral();
  fullTreeSgn->Draw("mass>>hSP",Form("weight*(%s && %s>0 && mass>%f && mass<%f && mcTrue)",bin_.c_str(),category_.c_str(),xLow_,xHigh_));
  fullTreeSgn->Draw("mass>>hSPall",Form("weight*(%s && %s>0 && mass>%f && mass<%f)",bin_.c_str(),category_.c_str(),xLow_,xHigh_));
  float SGNPass = hSP->Integral();
  float SGNPassall = hSPall->Integral();


  TTree *fullTreeBkg = (TTree*)fbkg.Get((tnp_+"/fitter_tree").c_str());
  // bkg
  TH1F* hB = new TH1F("hB","",1,0,150);
  TH1F* hBP = new TH1F("hBP","",1,0,150);
  fullTreeBkg->Draw("mass>>hB",Form("weight*(%s && mass>%f && mass<%f)",bin_.c_str(),xLow_,xHigh_));
  float BKG = hB->Integral();
  float BKGUnWeighted = hB->GetEntries();
  fullTreeBkg->Draw("mass>>hBP",Form("weight*(%s && %s>0 && mass>%f && mass<%f)",bin_.c_str(),category_.c_str(),xLow_,xHigh_));
  float BKGPass = hBP->Integral();
  float BKGPassUnWeighted = hBP->GetEntries();
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////
  TFile *test = new TFile( Form("testMC_%s_%s_%f.root",tnp_.c_str(),category_.c_str(),binCenter_),"RECREATE");

  RooWorkspace *w = new RooWorkspace("w","w");
  w->factory("mass[30,120]");
  w->factory("weight[0,10000]");
  w->factory("abseta[0,2.5]");
  w->factory("pt[0,200]");
  w->factory((category_+"[0,1]").c_str());
  //w->factory("RooGaussian::signalPdfP(mass,meanP[89,85,95.],sigmaP[5,0.5,10])");
  //w->factory("RooGaussian::signalPdfF(mass,meanF[91.2,85,95.],sigmaF[2,0.5,8])");
  /////////////!!! constrain n param to be > 0.5 !!!////////////////
  w->factory("RooCBShape::signalPdfF(mass,meanF[91.2,88,95.],sigmaF[3,0.5,8],alfaF[1.8,0.,10],nF[1.0,1e-06,10])");
  w->factory("RooCBShape::signalPdfP(mass,meanP[87,85,95.],sigmaP[5,3.0,8.5],alfaP[1.8,0.,10],nP[1.0,1e-06,10])");
  //w->factory("SUM::signalPdfP(fVBP[0.5,0,1]*RooBifurGauss::bifP(mass,meanP[91.2,88,95],sigmaLP[10,0.5,40],sigmaRP[0.]), RooVoigtian::voigP(mass, meanP, widthP[2.9], sigmaVoigP[5,0.1,10]) )");
  //w->factory("SUM::signalPdfF(fVBF[0.5,0,1]*RooBifurGauss::bifF(mass,meanF[91.2,88,95],sigmaLF[10,0.5,40],sigmaRF[0.]), RooVoigtian::voigF(mass, meanF, widthF[2.9], sigmaVoigF[5,0.1,10]) )");
  w->factory("RooExponential::backgroundPdfF(mass,cF[0,-10,10])");
  w->factory("RooExponential::backgroundPdfP(mass,cP[0,-10,10])");
  w->factory("efficiency[0.1,0,1]");
  w->factory("numSgn[0,100000]");
  w->factory("numBkgP[0,100000]");
  //w->factory("numBkgF[0,100000]");
  if(verbose_) cout << Form("numBkgF[%f,%f]",(0.5*BKGfracInF_)*(SGN-SGNPass),(1.5*BKGfracInF_)*(SGN-SGNPass)) << endl;
  w->factory( Form("numBkgF[%f,%f]",(0.5*BKGfracInF_)*(SGN-SGNPass),(1.5*BKGfracInF_)*(SGN-SGNPass)) );
  w->factory("expr::numSgnP('efficiency*numSgn',efficiency, numSgn)");
  w->factory("expr::numSgnF('(1-efficiency)*numSgn',efficiency, numSgn)");
  w->factory("passing[pass=1,fail=0]");
  w->factory("SUM::modelP(numSgnP*signalPdfP, numBkgP*backgroundPdfP)");
  w->factory("SUM::modelF(numSgnF*signalPdfF, numBkgF*backgroundPdfF)");
  w->factory("SIMUL::model(passing,pass=modelP,fail=modelF)");
  w->Print("V");
  w->saveSnapshot("clean", w->allVars());

  // soup tree
  TTree *fullTreeSoup = (TTree*)f.Get((tnp_+"/fitter_tree").c_str());

  // load the pdfs
  w->loadSnapshot("clean");
  RooRealVar  *weight = w->var("weight");
  RooRealVar  *abseta = w->var("abseta");
  RooRealVar  *pt = w->var("pt");
  RooRealVar  *mass = w->var("mass");
  mass->setRange(xLow_,xHigh_);
  RooRealVar  *cut = w->var( category_.c_str() );
  RooCategory *category = w->cat("passing");  
 
  //////////////////////////////////////////////////////////////
  RooRealVar nsigF("nsigF","num of signal events",100,0.,1000000) ;
  RooRealVar nbkgF("nbkgF","num of bkg events",(0.0*BKGfracInF_)*(SGN-SGNPass),(1.0*BKGfracInF_)*(SGN-SGNPass)) ;
  RooRealVar nsigP("nsigP","num of signal events",100,0.,1000000) ;
  RooRealVar nbkgP("nbkgP","num of bkg events",100,0.,1000000) ;
  //////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////
  RooDataSet dataP("dataP","dataset pass for the soup", RooArgSet(*mass,*weight,*abseta,*pt,*cut), Import( *fullTreeSoup ), /*WeightVar( *weight ),*/ Cut( Form("%s>0.5 && %s",category_.c_str(),bin_.c_str()) ) );
  RooDataSet dataF("dataF","dataset fail for the soup", RooArgSet(*mass,*weight,*abseta,*pt,*cut), Import( *fullTreeSoup ), /*WeightVar( *weight ),*/ Cut( Form("%s<0.5 && %s",category_.c_str(),bin_.c_str()) ) );

  RooDataSet combData("combData","combined data", RooArgSet(*mass,*weight,*abseta,*pt,*cut), Index(*category),Import("pass",dataP),Import("fail",dataF)) ;


  // model for passing
  RooAddPdf modelP("modelP","model for passing signal+bkg", RooArgList(*(w->pdf("signalPdfP")),*(w->pdf("backgroundPdfP"))),RooArgList(nsigP,nbkgP));
  RooAddPdf modelF("modelF","model for passing signal+bkg", RooArgList(*(w->pdf("signalPdfF")),*(w->pdf("backgroundPdfF"))),RooArgList(nsigF,nbkgF));

  RooFitResult *resCom = 0;
  RooFitResult *resIndP = 0;
  RooFitResult *resIndF = 0;

  if(doCombined_){
    // combined fit
    resCom = w->pdf("model")->fitTo(combData, Extended(1), Minos(1), Save(1),  SumW2Error( SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));
    resCom->Write("fitresults_combined");
  }
  else{
    // independent fit
    resIndP = modelP.fitTo(dataP, Extended(1), Minos(1), Save(1),  SumW2Error(  SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));
    resIndF = modelF.fitTo(dataF, Extended(1), Minos(1), Save(1),  SumW2Error(  SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));
    resIndP->Write("fitresults_independed_passed");
    resIndF->Write("fitresults_independent_failed");
  }


  RooRealVar *effFit = 0;
  RooRealVar *numSigFit = 0;
  RooRealVar *numBkgPFit = 0;
  RooRealVar *numBkgFFit = 0;

  if(doCombined_){
    RooArgSet fpf(resCom->floatParsFinal());
    effFit = (RooRealVar*)(&fpf["efficiency"]);
    numSigFit = (RooRealVar*)(&fpf["numSgn"]);
    numBkgPFit = (RooRealVar*)(&fpf["numBkgP"]);
    numBkgFFit = (RooRealVar*)(&fpf["numBkgF"]);  
  }
  else{
    RooArgSet fpfP(resIndP->floatParsFinal());
    RooArgSet fpfF(resIndF->floatParsFinal());
    RooRealVar *nsigP = (RooRealVar*)(&fpfP["nsigP"]);
    RooRealVar *nsigF = (RooRealVar*)(&fpfF["nsigF"]);
    numBkgPFit = (RooRealVar*)(&fpfP["nbkgP"]);
    numBkgFFit = (RooRealVar*)(&fpfF["nbkgF"]);
    numSigFit = (RooRealVar*)(&fpfF["nsigF"]);
    numSigFit->setVal( nsigP->getVal()+nsigF->getVal());
    numSigFit->setAsymError( TMath::Sqrt(nsigP->getErrorLo()*nsigP->getErrorLo() + nsigF->getErrorLo()*nsigF->getErrorLo())
			     , TMath::Sqrt(nsigP->getErrorHi()*nsigP->getErrorHi() + nsigF->getErrorHi()*nsigF->getErrorHi()));
    effFit = (RooRealVar*)(&fpfP["nsigP"]);
    effFit->setVal( nsigP->getVal()/numSigFit->getVal());
    //approx errors:
    effFit->setAsymError(nsigP->getErrorLo()/numSigFit->getVal(),nsigP->getErrorHi()/numSigFit->getVal() );
  }


  float errorLo = effFit->getErrorLo()<0 ? effFit->getErrorLo() : (-1)*effFit->getError();
  float errorHi = effFit->getErrorHi()>0 ? effFit->getErrorHi() : effFit->getError();
  float binError = TMath::Sqrt(SGNPass/SGN*(1-SGNPass/SGN)/SGN);
 
  Double_t* tnpMC = new Double_t[6];
  Double_t* truthMC = new Double_t[6];
  tnpMC[0] = binCenter_;
  tnpMC[1] = binWidth_;
  tnpMC[2] = binWidth_;
  tnpMC[3] = effFit->getVal();
  tnpMC[4] = (-1)*errorLo;
  tnpMC[5] = errorHi;
  truthMC[0] = binCenter_;
  truthMC[1] = binWidth_;
  truthMC[2] = binWidth_;
  truthMC[3] = SGNPass/SGN;
  truthMC[4] = binError;
  truthMC[5] = binError;

  pair< Double_t*,Double_t* > fitRes = make_pair( tnpMC, truthMC);

  if(verbose_) cout << "returning MC " << endl;
  return fitRes;

}

Double_t* simFitToData(const string tnp_ = "etoTauMargLooseNoCracks70", 
		       const string category_ = "tauAntiEMVA",
		       const string var_ = "abseta",
		       const string bin_ = "abseta < 1.5",
		       const float binCenter_ = 0.75,
		       const float binWidth_ = 0.75,
		       const float xLow_=60, 
		       const float xHigh_=120,
		       float BKGfracInF_=0.0035,
		       bool doCombined_ = true,
		       bool SumW2_ = false,
		       bool verbose_ = true){

  TFile fsgn("./test_sgn.root");
  TTree *fullTreeSgn = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());

  // file
  TFile f("../tagAndProbe/trees/38X/testNewWriteFromPAT_Data_runSplit.root");
  
  ///////////////////////////////////////////////////////////////
  TTree *fullTree = (TTree*)f.Get((tnp_+"/fitter_tree").c_str());
  // data tree
  TH1F* hS = new TH1F("hS","",1,0,150);
  TH1F* hSP = new TH1F("hSP","",1,0,150);
  fullTree->Draw("mass>>hS",Form("(%s && mass>%f && mass<%f)",bin_.c_str(),xLow_,xHigh_));
  float DATA = hS->Integral();
  fullTree->Draw("mass>>hSP",Form("(%s && %s>0 && mass>%f && mass<%f)",bin_.c_str(),category_.c_str(),xLow_,xHigh_));
  float DATAPass = hSP->Integral(); 
  ///////////////////////////////////////////////////////////////

  TFile *test = new TFile( Form("testData_%s_%s_%f.root",tnp_.c_str(),category_.c_str(),binCenter_),"RECREATE");

  RooWorkspace *w = new RooWorkspace("w","w");
  w->factory("mass[30,120]");
  w->factory("abseta[0,2.5]");
  w->factory("pt[0,200]");
  w->factory((category_+"[0,1]").c_str());
  //w->factory("RooGaussian::signalPdfP(mass,meanP[89,85,95.],sigmaP[5,0.5,10])");
  w->factory("RooGaussian::signalPdfF(mass,meanF[91.2,85,95.],sigmaF[2,0.5,8])");
  //w->factory("RooCBShape::signalPdfF(mass,meanF[91.2,88,95.],sigmaF[3,0.5,8],alfaF[1.8,0.,10],nF[1.0,0,10])");
  w->factory("RooCBShape::signalPdfP(mass,meanP[87,85,95.],sigmaP[5,3.0,8.5],alfaP[1.8,0.,10],nP[1.0,1e-06,10])");
  //w->factory("SUM::signalPdfP(fVBP[0.5,0,1]*RooBifurGauss::bifP(mass,meanP[91.2,88,95],sigmaLP[10,0.5,40],sigmaRP[0.]), RooVoigtian::voigP(mass, meanP, widthP[2.9], sigmaVoigP[5,0.1,10]) )");
  //w->factory("SUM::signalPdfF(fVBF[0.5,0,1]*RooBifurGauss::bifF(mass,meanF[91.2,88,95],sigmaLF[10,0.5,40],sigmaRF[0.]), RooVoigtian::voigF(mass, meanF, widthF[2.9], sigmaVoigF[5,0.1,10]) )");
  w->factory("SUM::signalPdfMC(fVBMC[0.5,0,1]*RooBifurGauss::bifMC(mass,meanMC[91.2,88,95],sigmaLMC[10,0.5,40],sigmaRMC[0.5,0,20]), RooVoigtian::voigMC(mass, meanMC, widthMC[2.49], sigmaVoigMC[5,0.1,10]) )");
  w->factory("RooExponential::backgroundPdfF(mass,cF[0,-10,10])");
  w->factory("RooExponential::backgroundPdfP(mass,cP[0,-10,10])");
  w->factory("efficiency[0.1,0,1]");
  w->factory("numSgn[0,100000]");
  w->factory("numBkgP[0,100000]");
  //w->factory("numBkgF[0,100000]");
  if(verbose_) cout << Form("numBkgF[%f,%f]",(BKGfracInF_*0.5)*(DATA-DATAPass),(BKGfracInF_*1.5)*(DATA-DATAPass)) << endl;
  w->factory( Form("numBkgF[%f,%f]",(BKGfracInF_*0.5)*(DATA-DATAPass),(BKGfracInF_*1.5)*(DATA-DATAPass)) );
  w->factory("expr::numSgnP('efficiency*numSgn',efficiency, numSgn)");
  w->factory("expr::numSgnF('(1-efficiency)*numSgn',efficiency, numSgn)");
  w->factory("passing[pass=1,fail=0]");
  w->factory("SUM::modelP(numSgnP*signalPdfP, numBkgP*backgroundPdfP)");
  w->factory("SUM::modelF(numSgnF*signalPdfF, numBkgF*backgroundPdfF)");
  w->factory("SIMUL::model(passing,pass=modelP,fail=modelF)");
  w->Print("V");
  w->saveSnapshot("clean", w->allVars());

  // load the pdfs
  w->loadSnapshot("clean");
  RooRealVar  *abseta = w->var("abseta");
  RooRealVar  *pt = w->var("pt");
  RooRealVar  *mass = w->var("mass");
  mass->setRange(xLow_,xHigh_);
  RooRealVar  *cut = w->var( category_.c_str() );
  RooCategory *category = w->cat("passing");  
 
  //////////////////////////////////////////////////////////////
  RooRealVar nsigF("nsigF","num of signal events",100,0.,1000000) ;
  RooRealVar nbkgF("nbkgF","num of bkg events",(0.5*BKGfracInF_)*(DATA-DATAPass),(1.5*BKGfracInF_)*(DATA-DATAPass)) ;
  RooRealVar nsigP("nsigP","num of signal events",100,0.,1000000) ;
  RooRealVar nbkgP("nbkgP","num of bkg events",100,0.,1000000) ;
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  RooDataSet dataSgn("dataSgn","dataset pass for the signal", RooArgSet(*mass,*abseta,*pt,*cut), Import( *fullTreeSgn ),Cut( Form("%s>0.5 && %s",category_.c_str(),bin_.c_str()) ) );
  RooFitResult* resSignal = w->pdf("signalPdfMC")->fitTo(dataSgn, Minos(1), Save(1),  SumW2Error( SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));
  RooArgSet fpfSgn(resSignal->floatParsFinal());
  RooRealVar* fVBMC_       = (RooRealVar*)(&fpfSgn["fVBMC"]);
  RooRealVar* meanMC_      = (RooRealVar*)(&fpfSgn["meanMC"]);
  RooRealVar* sigmaLMC_    = (RooRealVar*)(&fpfSgn["sigmaLMC"]);
  RooRealVar* sigmaRMC_    = (RooRealVar*)(&fpfSgn["sigmaRMC"]);
  RooRealVar* sigmaVoigMC_ = (RooRealVar*)(&fpfSgn["sigmaVoigMC"]);
  RooRealVar widthC("widthMC_","widthMC",2.49);
  // constrained voigtian
  RooRealVar meanC("meanC","mean",91,85,95/*meanMC_->getVal()*/);
  RooRealVar sigmaVoigC("sigmaVoigC","sigmaVoig",3,0.5,10/*sigmaVoigMC_->getVal()*/);
  RooVoigtian voigC("signalVoigC","signalVoig",*mass,meanC,widthC,sigmaVoigC);
  // constrained bifurcated gaussian
  RooRealVar sigmaRC("sigmaRC","sigmaR",sigmaRMC_->getVal());
  RooRealVar sigmaLC("sigmaLC","sigmaL",sigmaLMC_->getVal());
  RooBifurGauss bifurcC("bifurcC","bifurc",*mass,meanC,sigmaLC,sigmaRC); 
  RooRealVar fVBC("fVBC","FVB",fVBMC_->getVal());
  RooAddPdf voigPlusBifurcC("voigPlusBifurcC","voigPlusBifurc",RooArgList(voigC,bifurcC),fVBC);
  //////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////
  RooDataSet dataP("dataP","dataset pass for the soup", RooArgSet(*mass,*abseta,*pt,*cut), Import( *fullTree ),Cut( Form("%s>0.5 && %s",category_.c_str(),bin_.c_str()) ) );
  RooDataSet dataF("dataF","dataset fail for the soup", RooArgSet(*mass,*abseta,*pt,*cut), Import( *fullTree ), Cut( Form("%s<0.5 && %s",category_.c_str(),bin_.c_str()) ) );

  RooDataSet combData("combData","combined data", RooArgSet(*mass,*abseta,*pt,*cut), Index(*category),Import("pass",dataP),Import("fail",dataF)) ;


  // model for passing
  RooAddPdf modelP("modelP","model for passing signal+bkg", RooArgList(voigPlusBifurcC /**(w->pdf("signalPdfP"))*/,*(w->pdf("backgroundPdfP"))),RooArgList(nsigP,nbkgP));
  RooAddPdf modelF("modelF","model for passing signal+bkg", RooArgList(*(w->pdf("signalPdfF")),*(w->pdf("backgroundPdfF"))),RooArgList(nsigF,nbkgF));

  RooFitResult *resCom = 0;
  RooFitResult *resIndP = 0;
  RooFitResult *resIndF = 0;

  if(doCombined_){
    // combined fit
    resCom = w->pdf("model")->fitTo(combData, Extended(1), Minos(1), Save(1),  SumW2Error( SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));
    resCom->Write("fitresults_combined");
  }
  else{
    // independent fit
    resIndP = modelP.fitTo(dataP, Extended(1), Minos(1), Save(1),  SumW2Error(  SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));
    resIndF = modelF.fitTo(dataF, Extended(1), Minos(1), Save(1),  SumW2Error(  SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));
    resIndP->Write("fitresults_independed_passed");
    resIndF->Write("fitresults_independent_failed");
  }

  /////////////////////////////////////////////
  RooPlot* frame = mass->frame();
  dataP.plotOn(frame);
  modelP.plotOn(frame,LineColor(kRed));
  voigPlusBifurcC.plotOn(frame);
  frame->Draw();
  /////////////////////////////////////////////


  RooRealVar *effFit = 0;
  RooRealVar *numSigFit = 0;
  RooRealVar *numBkgPFit = 0;
  RooRealVar *numBkgFFit = 0;

  if(doCombined_){
    RooArgSet fpf(resCom->floatParsFinal());
    effFit = (RooRealVar*)(&fpf["efficiency"]);
    numSigFit = (RooRealVar*)(&fpf["numSgn"]);
    numBkgPFit = (RooRealVar*)(&fpf["numBkgP"]);
    numBkgFFit = (RooRealVar*)(&fpf["numBkgF"]);  
  }
  else{
    RooArgSet fpfP(resIndP->floatParsFinal());
    RooArgSet fpfF(resIndF->floatParsFinal());
    RooRealVar *nsigP = (RooRealVar*)(&fpfP["nsigP"]);
    RooRealVar *nsigF = (RooRealVar*)(&fpfF["nsigF"]);
    numBkgPFit = (RooRealVar*)(&fpfP["nbkgP"]);
    numBkgFFit = (RooRealVar*)(&fpfF["nbkgF"]);
    numSigFit = (RooRealVar*)(&fpfF["nsigF"]);
    numSigFit->setVal( nsigP->getVal()+nsigF->getVal());
    numSigFit->setAsymError( TMath::Sqrt(nsigP->getErrorLo()*nsigP->getErrorLo() + nsigF->getErrorLo()*nsigF->getErrorLo())
			     , TMath::Sqrt(nsigP->getErrorHi()*nsigP->getErrorHi() + nsigF->getErrorHi()*nsigF->getErrorHi()));
    effFit = (RooRealVar*)(&fpfP["nsigP"]);
    effFit->setVal( nsigP->getVal()/numSigFit->getVal());
    //approx errors:
    effFit->setAsymError(nsigP->getErrorLo()/numSigFit->getVal(),nsigP->getErrorHi()/numSigFit->getVal() );
    effFit->setError(nsigP->getError()/numSigFit->getVal());
  }

  ///////////////////////////////////////////////////////
 
  double binsEdges[3] = {0.0,1.5,2.5};
  TH1F* h1 = new TH1F("h1","Validator of tnp",2, binsEdges);
  h1->SetXTitle( var_.c_str() );
  h1->SetYTitle( ("Efficiency of passing "+category_).c_str() );

 
  float errorLo = effFit->getErrorLo()<0 ? effFit->getErrorLo() : (-1)*effFit->getError();
  float errorHi = effFit->getErrorHi()>0 ? effFit->getErrorHi() : effFit->getError();
 
  Double_t* tnpDATA = new Double_t[6];
  tnpDATA[0] = binCenter_;
  tnpDATA[1] = binWidth_;
  tnpDATA[2] = binWidth_;
  tnpDATA[3] = effFit->getVal();
  tnpDATA[4] = (-1)*errorLo;
  tnpDATA[5] = errorHi;

  if(verbose_) cout << "returning data: " << binCenter_ << ";" << effFit->getVal()  << endl;
  if(verbose_) cout << "data " << DATA << " , data pass " << DATAPass << endl;

  return tnpDATA;

}



void makePlot(const string tnp_ = "etoTauMargLooseNoCracks70",
	      const string category_ = "tauAntiEMVA",
	      const string var_ = "abseta",
	      const float xLow_=60, 
	      const float xHigh_=120,
	      const float BKGfracInF_ = 0.0035,
	      bool doCombined_ = false,
	      bool SumW2_ = true,
	      bool verbose_ = true,
	      double bin1_ = 0.0,
	      double bin2_ = 1.5,
	      double bin3_ = 2.5
	      ){

  TCanvas *c1 = new TCanvas("c1","Canvas",10,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.4,0.6,0.89,0.89,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);

  double binsEdges[3] = {bin1_,bin2_,bin3_};
  TH1F* h1 = new TH1F("h1","Validator of tnp",2, binsEdges);
  h1->SetXTitle( var_.c_str() );
  h1->SetYTitle( ("Efficiency of passing "+category_).c_str() );

  /////////////////////////////////////////////////////////////////
  //pair<Double_t*,Double_t*> MC_bin1 = simFitToMC(tnp_, category_,var_,string(Form("abseta < %f",bin2_)),(bin2_+bin1_)/2,(bin2_-bin1_)/2, xLow_, xHigh_, BKGfracInF_, doCombined_,SumW2_, verbose_);
  cout << "******************** MC_bin1" << endl;
  Double_t* DATA_bin1 = simFitToData(tnp_, category_,var_,string(Form("abseta < %f",bin2_)),(bin2_+bin1_)/2,(bin2_-bin1_)/2 ,xLow_, xHigh_, BKGfracInF_, doCombined_,SumW2_, verbose_);

  return;

  cout << "******************** DATA_bin1" << endl;
  //pair<Double_t*,Double_t*> MC_bin2 = simFitToMC(tnp_, category_,var_,string(Form("(abseta > %f && abseta<%f)",bin2_,bin3_)),(bin3_+bin2_)/2,(bin3_-bin2_)/2, xLow_, xHigh_, BKGfracInF_, doCombined_,SumW2_, verbose_);
  cout << "******************** MC_bin2" << endl;
  Double_t* DATA_bin2 = simFitToData(tnp_, category_,var_,string(Form("(abseta < %f && abseta>%f)",bin3_,bin2_)),(bin3_+bin2_)/2,(bin3_-bin2_)/2, xLow_, xHigh_, BKGfracInF_, doCombined_,SumW2_, verbose_);
  cout << "******************** Data_bin2" << endl;
  /////////////////////////////////////////////////////////////////

  /*
  Double_t tnpMC_x[2]  = {(MC_bin1.first)[0],(MC_bin2.first)[0]};
  Double_t tnpMC_xL[2] = {(MC_bin1.first)[1],(MC_bin2.first)[1]};
  Double_t tnpMC_xH[2] = {(MC_bin1.first)[2],(MC_bin2.first)[2]};
  Double_t tnpMC_y[2]  = {(MC_bin1.first)[3],(MC_bin2.first)[3]};
  Double_t tnpMC_yL[2] = {(MC_bin1.first)[4],(MC_bin2.first)[4]};
  Double_t tnpMC_yH[2] = {(MC_bin1.first)[5],(MC_bin2.first)[5]};
  //
  
  Double_t truthMC_x[2]  = {(MC_bin1.second)[0],(MC_bin2.second)[0]};
  Double_t truthMC_xL[2] = {(MC_bin1.second)[1],(MC_bin2.second)[1]};
  Double_t truthMC_xH[2] = {(MC_bin1.second)[2],(MC_bin2.second)[2]};
  Double_t truthMC_y[2]  = {(MC_bin1.second)[3],(MC_bin2.second)[3]};
  Double_t truthMC_yL[2] = {(MC_bin1.second)[4],(MC_bin2.second)[4]};
  Double_t truthMC_yH[2] = {(MC_bin1.second)[5],(MC_bin2.second)[5]};
  //
  */
  Double_t tnpDATA_x[2]  = {DATA_bin1[0],DATA_bin2[0]};
  Double_t tnpDATA_xL[2] = {DATA_bin1[1],DATA_bin2[1]};
  Double_t tnpDATA_xH[2] = {DATA_bin1[2],DATA_bin2[2]};
  Double_t tnpDATA_y[2]  = {DATA_bin1[3],DATA_bin2[3]};
  Double_t tnpDATA_yL[2] = {DATA_bin1[4],DATA_bin2[4]};
  Double_t tnpDATA_yH[2] = {DATA_bin1[5],DATA_bin2[5]};
  ///////////////////////////////////////////////////////
  /*
  TGraphAsymmErrors* graph_truthMC = new TGraphAsymmErrors(2,truthMC_x,truthMC_y, truthMC_xL,truthMC_xH,truthMC_yL,truthMC_yH);
  graph_truthMC->SetMarkerColor(kBlue);
  graph_truthMC->SetMarkerStyle(20);
  graph_truthMC->SetMarkerSize(1);
  ///////////////////////////////////////////////////////
  TGraphAsymmErrors* graph_tnpMC = new TGraphAsymmErrors(2,tnpMC_x,tnpMC_y,tnpMC_xL,tnpMC_xH,tnpMC_yL,tnpMC_yH);
  graph_tnpMC->SetMarkerColor(kRed);
  graph_tnpMC->SetMarkerStyle(21);
  graph_tnpMC->SetMarkerSize(1);
  */
  ///////////////////////////////////////////////////////
  TGraphAsymmErrors* graph_tnpDATA = new TGraphAsymmErrors(2,tnpDATA_x, tnpDATA_y,tnpDATA_xL,tnpDATA_xH,tnpDATA_yL,tnpDATA_yH);
  graph_tnpDATA->SetMarkerColor(kBlack);
  graph_tnpDATA->SetMarkerStyle(22);
  graph_tnpDATA->SetMarkerSize(1);

  h1->Draw("");
  //graph_truthMC->Draw("PSAME");
  //graph_tnpMC->Draw("PSAME");
  graph_tnpDATA->Draw("PSAME");

  c1->Draw();

}

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
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"

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

vector<Double_t*> simFit(bool makeSoupFit_ = false,
			 const string tnp_ = "etoTauMargLooseNoCracks70", 
			 const string category_ = "tauAntiEMVA",
			 const string bin_ = "abseta<1.5",
			 const float binCenter_ = 0.75,
			 const float binWidth_ = 0.75,
			 const float xLow_=60, 
			 const float xHigh_=120,
			 bool SumW2_ = false,
			 bool verbose_ = true){

  vector<Double_t*> out;
  //return out;

  //TFile *test = new TFile( outFile->GetName(),"UPDATE");
  // output file
  TFile *test = new TFile( Form("EtoTauPlotsFit_%s_%s_%f_v5.root",tnp_.c_str(),category_.c_str(),binCenter_),"RECREATE");
  test->mkdir(Form("bin%.2f",binCenter_));

  TCanvas *c = new TCanvas("fitCanvas",Form("fitCanvas_%s_%s",tnp_.c_str(),bin_.c_str()),10,30,650,600);
  c->SetGrid(0,0);
  c->SetFillStyle(4000);
  c->SetFillColor(10);
  c->SetTicky();
  c->SetObjectStat(0);
  
  TCanvas *c2 = new TCanvas("fitCanvasTemplate",Form("fitCanvasTemplate_%s_%s",tnp_.c_str(),bin_.c_str()),10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  // input files
  TFile fsup("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup.root");
  TFile fbkg("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_bkg.root");
  TFile fsgn("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_sgn.root");
  TFile fdat("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_Data.root");
  //TFile fdat("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_Run2010A.root");

  //********************** signal only tree *************************/

  TTree *fullTreeSgn = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());
  TH1F* hSall = new TH1F("hSall","",1,0,150);
  TH1F* hSPall = new TH1F("hSPall","",1,0,150);
  TH1F* hS = new TH1F("hS","",1,0,150);
  TH1F* hSP = new TH1F("hSP","",1,0,150);
  fullTreeSgn->Draw("mass>>hS",Form("weight*(%s && mass>%f && mass<%f && mcTrue && signalPFChargedHadrCands<1.5)",bin_.c_str(),xLow_,xHigh_));
  fullTreeSgn->Draw("mass>>hSall",Form("weight*(%s && mass>%f && mass<%f)",bin_.c_str(),xLow_,xHigh_));

  float SGNtrue = hS->Integral();
  float SGNall  = hSall->Integral();
 
  fullTreeSgn->Draw("mass>>hSP",Form("weight*(%s && %s>0 && mass>%f && mass<%f && mcTrue && signalPFChargedHadrCands<1.5 )",bin_.c_str(),category_.c_str(),xLow_,xHigh_));
  fullTreeSgn->Draw("mass>>hSPall",Form("weight*(%s && %s>0 && mass>%f && mass<%f && signalPFChargedHadrCands<1.5 )",bin_.c_str(),category_.c_str(),xLow_,xHigh_));

  float SGNtruePass = hSP->Integral();
  float SGNallPass  = hSPall->Integral();

  //********************** background only tree *************************//

  TTree *fullTreeBkg = (TTree*)fbkg.Get((tnp_+"/fitter_tree").c_str());
  TH1F* hB = new TH1F("hB","",1,0,150);
  TH1F* hBP = new TH1F("hBP","",1,0,150);
  fullTreeBkg->Draw("mass>>hB",Form("weight*(%s && mass>%f && mass<%f && signalPFChargedHadrCands<1.5 )",bin_.c_str(),xLow_,xHigh_));
 
  float BKG = hB->Integral();
  float BKGUnWeighted = hB->GetEntries();
  
  fullTreeBkg->Draw("mass>>hBP",Form("weight*(%s && %s>0 && mass>%f && mass<%f && signalPFChargedHadrCands<1.5 )",bin_.c_str(),category_.c_str(),xLow_,xHigh_));
  
  float BKGPass = hBP->Integral();
  float BKGUnWeightedPass = hBP->GetEntries();

  //********************** soup tree *************************//

  TTree *fullTreeSoup = (TTree*)fsup.Get((tnp_+"/fitter_tree").c_str());

  //********************** data tree *************************//

  TTree *fullTreeData = (TTree*)fdat.Get((tnp_+"/fitter_tree").c_str());

  //********************** workspace ***********************//

  RooWorkspace *w = new RooWorkspace("w","w");
  // tree variables to be imported
  w->factory("mass[30,120]");
  w->factory("weight[0,10000]");
  w->factory("abseta[0,2.5]");
  w->factory("pt[0,200]");
  w->factory("mcTrue[0,1]");
  w->factory("signalPFChargedHadrCands[0,10]");
  w->factory((category_+"[0,1]").c_str());
  // background pass pdf for MC
  w->factory("RooExponential::McBackgroundPdfP(mass,McCP[0,-10,10])");
  // background fail pdf for MC
  w->factory("RooExponential::McBackgroundPdfF(mass,McCF[0,-10,10])");
  // background pass pdf for Data
  w->factory("RooExponential::DataBackgroundPdfP(mass,DataCP[-8.8761e-02,-1,1])");
  // background fail pdf for Data
  w->factory("RooExponential::DataBackgroundPdfF(mass,DataCF[-5.3356e-02,-1,1])");
  // fit parameters for background
  w->factory("McEfficiency[0.1,0,1]");
  w->factory("McNumSgnP[0,100000]");
  w->factory("McNumSgnF[0,1000000]");
  w->factory("McNumBkgP[0,100000]");
  w->factory("McNumBkgF[0,100000]"); 
  //w->factory("expr::McNumSgnP('McEfficiency*McNumSgn',McEfficiency,McNumSgn)");
  //w->factory("expr::McNumSgnF('(1-McEfficiency)*McNumSgn',McEfficiency,McNumSgn)");
  w->factory("McPassing[pass=1,fail=0]");
  // fit parameters for data
  w->factory("DataEfficiency[0.05,0,1]");
  w->factory("DataNumSgnP[200,0,1000000]");
  w->factory("DataNumSgnF[10000,0,1000000]");
  w->factory("DataNumBkgP[200,0,1000000]");
  w->factory("DataNumBkgF[100,0,1000000]");
  //w->factory("expr::DataNumSgnP('DataEfficiency*DataNumSgn',DataEfficiency,DataNumSgn)");
  //w->factory("expr::DataNumSgnF('(1-DataEfficiency)*DataNumSgn',DataEfficiency,DataNumSgn)");
  w->factory("DataPassing[pass=1,fail=0]");

  RooRealVar  *weight = w->var("weight");
  RooRealVar  *abseta = w->var("abseta");
  RooRealVar  *pt = w->var("pt");
  RooRealVar  *mass = w->var("mass");
  mass->setRange(xLow_,xHigh_);
  RooRealVar  *mcTrue = w->var("mcTrue");
  RooRealVar  *cut = w->var( category_.c_str() );
  RooRealVar *signalPFChargedHadrCands = w->var("signalPFChargedHadrCands");
 
  // build the template for the signal pass sample:
  RooDataSet templateP("templateP","dataset for signal-pass template", RooArgSet(*mass,*weight,*abseta,*pt,*cut,*mcTrue,*signalPFChargedHadrCands), Import( *fullTreeSgn ), /*WeightVar( *weight ),*/ Cut( Form("(mcTrue && %s>0.5 && %s && signalPFChargedHadrCands<1.5)",category_.c_str(),bin_.c_str()) ) );
  // build the template for the signal fail sample:
  RooDataSet templateF("templateF","dataset for signal-fail template", RooArgSet(*mass,*weight,*abseta,*pt,*cut,*mcTrue,*signalPFChargedHadrCands), Import( *fullTreeSgn ), /*WeightVar( *weight ),*/ Cut( Form("(mcTrue && %s<0.5 && %s && signalPFChargedHadrCands<1.5)",category_.c_str(),bin_.c_str()) ) );
  

  mass->setMin(60); mass->setMax(120); mass->setBins(20);
  RooDataHist templateHistP("templateHistP","",RooArgSet(*mass),templateP,1.0);
  RooHistPdf TemplateSignalPdfP("TemplateSignalPdfP","",RooArgSet(*mass),templateHistP);
  w->import(TemplateSignalPdfP);

  mass->setBins(60);
  RooDataHist templateHistF("templateHistF","",RooArgSet(*mass),templateF,1.0);
  RooHistPdf TemplateSignalPdfF("TemplateSignalPdfF","",RooArgSet(*mass),templateHistF);
  w->import(TemplateSignalPdfF);

  mass->setBins(10000,"fft");

  RooPlot* TemplateFrameP = mass->frame(Title("Template passing"));
  templateP.plotOn(TemplateFrameP);
  w->pdf("TemplateSignalPdfP")->plotOn(TemplateFrameP);

  RooPlot* TemplateFrameF = mass->frame(Title("Template failing"));
  templateF.plotOn(TemplateFrameF);
  w->pdf("TemplateSignalPdfF")->plotOn(TemplateFrameF);


  w->factory("RooFFTConvPdf::McSignalPdfP(mass,TemplateSignalPdfP,RooGaussian::McResolModP(mass,McMeanResP[0.5,0.,10.],McSigmaResP[0.5,0.,10]))");
  w->factory("RooFFTConvPdf::McSignalPdfF(mass,TemplateSignalPdfF,RooGaussian::McResolModF(mass,McMeanResF[0.5,0.,10.],McSigmaResF[0.5,0.,10]))");
  w->factory("RooFFTConvPdf::DataSignalPdfP(mass,TemplateSignalPdfP,RooGaussian::DataResolModP(mass,DataMeanResP[0.5,0.,10.],DataSigmaResP[0.5,0.,10]))");
  //w->factory("RooCBShape::DataSignalPdfF(mass,DataMeanF[91.2,88,95.],DataSigmaF[3,0.5,8],DataAlfaF[1.8,0.,10],DataNF[1.0,1e-06,10])");
  //w->factory("RooFFTConvPdf::DataSignalPdfF(mass,TemplateSignalPdfF,RooGaussian::DataResolModF(mass,DataMeanResF[0.5,0.,10.],DataSigmaResF[0.5,0.,10]))");
  w->factory("RooFFTConvPdf::DataSignalPdfF(mass,RooVoigtian::DataVoigF(mass,DataMeanF[90.8,89,92],DataWidthF[2.49],DataSigmaF[1.4129e+00,0.5,3]),RooCBShape::DataResolModF(mass,DataMeanResF[1.3920e-01,0.,3],DataSigmaResF[2.5792e-01,0.,5],DataAlphaResF[2.7043e-01,0.,1],DataNResF[2.4538e+00,1,5]))");
 
  

  // composite model pass for MC
  w->factory("SUM::McModelP(McNumSgnP*McSignalPdfP,McNumBkgP*McBackgroundPdfP)");  
  w->factory("SUM::McModelF(McNumSgnF*McSignalPdfF,McNumBkgF*McBackgroundPdfF)"); 
  // composite model pass for data
  w->factory("SUM::DataModelP(DataNumSgnP*DataSignalPdfP,DataNumBkgP*DataBackgroundPdfP)");  
  w->factory("SUM::DataModelF(DataNumSgnF*DataSignalPdfF,DataNumBkgF*DataBackgroundPdfF)");  
  // simultaneous fir for MC
  w->factory("SIMUL::McModel(McPassing,pass=McModelP,fail=McModelF)");
  // simultaneous fir for data
  w->factory("SIMUL::DataModel(DataPassing,pass=DataModelP,fail=DataModelF)");
  w->Print("V");
  w->saveSnapshot("clean", w->allVars());

  w->loadSnapshot("clean");
  /****************** sim fit to soup **************************/

  ///////////////////////////////////////////////////////////////
  RooDataSet McDataP("McDataP","dataset pass for the soup", RooArgSet(*mass,*weight,*abseta,*pt,*cut,*signalPFChargedHadrCands), Import( *fullTreeSoup ), /*WeightVar( *weight ),*/ Cut( Form("(%s>0.5 && %s && signalPFChargedHadrCands<1.5)",category_.c_str(),bin_.c_str()) ) );
  RooDataSet McDataF("McDataF","dataset fail for the soup", RooArgSet(*mass,*weight,*abseta,*pt,*cut,*signalPFChargedHadrCands), Import( *fullTreeSoup ), /*WeightVar( *weight ),*/ Cut( Form("(%s<0.5 && %s && signalPFChargedHadrCands<1.5)",category_.c_str(),bin_.c_str()) ) );
  RooDataSet McCombData("McCombData","combined data for the soup", RooArgSet(*mass,*weight,*abseta,*pt,*cut,*signalPFChargedHadrCands), Index(*(w->cat("McPassing"))), Import("pass",McDataP),Import("fail",McDataF)) ;

  RooPlot* McFrameP = 0;
  RooPlot* McFrameF = 0;
  RooRealVar McEffFit("McEffFit","McEffFit",0.5);
  if(makeSoupFit_){

    RooFitResult* ResMcFitP = w->pdf("McModelP")->fitTo(McDataP, Extended(1), Minos(1), Save(1),  SumW2Error( SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));
    test->cd(Form("bin%.2f",binCenter_));
    ResMcFitP->Write("McFitResults_Passed");
    RooFitResult* ResMcFitF = w->pdf("McModelF")->fitTo(McDataF, Extended(1), Minos(1), Save(1),  SumW2Error( SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));
    test->cd(Form("bin%.2f",binCenter_));
    ResMcFitF->Write("McFitResults_Failed");

    RooArgSet McFitParamP(ResMcFitP->floatParsFinal());
    RooArgSet McFitParamF(ResMcFitF->floatParsFinal());    

    RooRealVar* McNumSigPFit = (RooRealVar*)(&McFitParamP["McNumSgnP"]);
    RooRealVar* McNumSigFFit = (RooRealVar*)(&McFitParamF["McNumSgnF"]);
    RooRealVar* McNumBkgPFit = (RooRealVar*)(&McFitParamP["McNumBkgP"]);
    RooRealVar* McNumBkgFFit = (RooRealVar*)(&McFitParamF["McNumBkgF"]);

    McEffFit.setVal(McNumSigPFit->getVal()/(McNumSigPFit->getVal()+McNumSigFFit->getVal()));
    McEffFit.setError( 1./(McNumSigPFit->getVal()+McNumSigFFit->getVal())*1./(McNumSigPFit->getVal()+McNumSigFFit->getVal())*TMath::Sqrt(McNumSigPFit->getVal()*McNumSigPFit->getVal()*McNumSigFFit->getError()+McNumSigFFit->getVal()*+McNumSigFFit->getVal()*McNumSigPFit->getError()) );
    McEffFit.setAsymError( (-1)*McEffFit.getError(), McEffFit.getError()  );

    McFrameP = mass->frame(Bins(30),Title("MC: passing sample"));
    McDataP.plotOn(McFrameP);
    w->pdf("McModelP")->plotOn(McFrameP,LineColor(kBlue),Range(60,120));
    w->pdf("McModelP")->plotOn(McFrameP, Components("McSignalPdfP"), LineColor(kRed),Range(60,120));
    w->pdf("McModelP")->plotOn(McFrameP, Components("McBackgroundPdfP"), LineColor(kGreen),LineStyle(kDashed),Range(60,120));


    McFrameF = mass->frame(Bins(30));
    McDataF.plotOn(McFrameF);
    w->pdf("McModelF")->plotOn(McFrameF,LineColor(kBlue),Range(60,120));
    w->pdf("McModelF")->plotOn(McFrameF,Components("McSignalPdfF"), LineColor(kRed),Range(60,120)); 
    w->pdf("McModelF")->plotOn(McFrameF,Components("McBackgroundPdfF"), LineColor(kGreen),LineStyle(kDashed),Range(60,120)); 
  }
  
  ///////////////////////////////////////////////////////////////

  /****************** sim fit to data **************************/

  ///////////////////////////////////////////////////////////////
  RooDataSet DataDataP("DataDataP","dataset pass for the soup", RooArgSet(*mass,*weight,*abseta,*pt,*cut,*signalPFChargedHadrCands), Import( *fullTreeData ), /*WeightVar( *weight ),*/ Cut( Form("(%s>0.5 && %s && signalPFChargedHadrCands<1.5)",category_.c_str(),bin_.c_str()) ) );
  RooDataSet DataDataF("DataDataF","dataset fail for the soup", RooArgSet(*mass,*weight,*abseta,*pt,*cut,*signalPFChargedHadrCands), Import( *fullTreeData ), /*WeightVar( *weight ),*/ Cut( Form("(%s<0.5 && %s && signalPFChargedHadrCands<1.5)",category_.c_str(),bin_.c_str()) ) );
  RooDataSet DataCombData("DataCombData","combined data for the soup", RooArgSet(*mass,*weight,*abseta,*pt,*cut,*signalPFChargedHadrCands), Index(*(w->cat("DataPassing"))), Import("pass",DataDataP),Import("fail",DataDataF)) ;

  RooPlot* DataFrameP = 0;
  RooPlot* DataFrameF = 0;
  RooRealVar DataEffFit("DataEffFit","DataEffFit",0.5);

  RooFitResult* ResDataFitP = w->pdf("DataModelP")->fitTo(DataDataP, Extended(1), Minos(1), Save(1),  SumW2Error( SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));
  test->cd(Form("bin%.2f",binCenter_));
  ResDataFitP->Write("DataFitResults_Passed");
  RooFitResult* ResDataFitF = w->pdf("DataModelF")->fitTo(DataDataF, Extended(1), Minos(1), Save(1),  SumW2Error( SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));
  test->cd(Form("bin%.2f",binCenter_));
  ResDataFitF->Write("DataFitResults_Failed");

  RooArgSet DataFitParamP(ResDataFitP->floatParsFinal());
  RooArgSet DataFitParamF(ResDataFitF->floatParsFinal());    
  
  RooRealVar* DataNumSigPFit = (RooRealVar*)(&DataFitParamP["DataNumSgnP"]);
  RooRealVar* DataNumSigFFit = (RooRealVar*)(&DataFitParamF["DataNumSgnF"]);
  RooRealVar* DataNumBkgPFit = (RooRealVar*)(&DataFitParamP["DataNumBkgP"]);
  RooRealVar* DataNumBkgFFit = (RooRealVar*)(&DataFitParamF["DataNumBkgF"]);
  
  DataEffFit.setError( 1./(DataNumSigPFit->getVal()+DataNumSigFFit->getVal())*1./(DataNumSigPFit->getVal()+DataNumSigFFit->getVal())*TMath::Sqrt(DataNumSigPFit->getVal()*DataNumSigPFit->getVal()*DataNumSigFFit->getError()+DataNumSigFFit->getVal()*+DataNumSigFFit->getVal()*DataNumSigPFit->getError()) );
  DataEffFit.setAsymError( (-1)*DataEffFit.getError(),  DataEffFit.getError() );

  
  DataFrameP = mass->frame(Bins(30),Title("MC: passing sample"));
  DataDataP.plotOn(DataFrameP);
  w->pdf("DataModelP")->plotOn(DataFrameP,LineColor(kBlue),Range(60,120));
  w->pdf("DataModelP")->plotOn(DataFrameP, Components("DataSignalPdfP"), LineColor(kRed),Range(60,120));
  w->pdf("DataModelP")->plotOn(DataFrameP, Components("DataBackgroundPdfP"), LineColor(kGreen),LineStyle(kDashed),Range(60,120));


  DataFrameF = mass->frame(Bins(30));
  DataDataF.plotOn(DataFrameF);
  w->pdf("DataModelF")->plotOn(DataFrameF,LineColor(kBlue),Range(60,120));
  w->pdf("DataModelF")->plotOn(DataFrameF,Components("DataSignalPdfF"), LineColor(kRed),Range(60,120)); 
  w->pdf("DataModelF")->plotOn(DataFrameF,Components("DataBackgroundPdfF"), LineColor(kGreen),LineStyle(kDashed),Range(60,120));
  
  ///////////////////////////////////////////////////////////////

 
  c->Divide(2,1);

 

  c->cd(1);
  DataFrameP->Draw();
  c->cd(2);
  DataFrameF->Draw();
  if(makeSoupFit_){
    c->Divide(2,2);
    c->cd(3);
    McFrameP->Draw();
    c->cd(4);
    McFrameF->Draw();
  }
 
  c->Draw();
 
  test->cd(Form("bin%f",binCenter_));
 
  c->Write();
 
  c2->Divide(2,1);
  c2->cd(1);
  TemplateFrameP->Draw();
  c2->cd(2);
  TemplateFrameF->Draw();
  c2->Draw();
 
  test->cd(Form("bin%f",binCenter_));
  c2->Write();


  // MINOS errors, otherwise HESSE quadratic errors
  float McErrorLo = 0;
  float McErrorHi = 0;
  if(makeSoupFit_){
    McErrorLo = McEffFit.getErrorLo()<0 ? McEffFit.getErrorLo() : (-1)*McEffFit.getError();
    McErrorHi = McEffFit.getErrorHi()>0 ? McEffFit.getErrorHi() : McEffFit.getError();
  }
  float DataErrorLo = DataEffFit.getErrorLo()<0 ? DataEffFit.getErrorLo() : (-1)*DataEffFit.getError();
  float DataErrorHi = DataEffFit.getErrorHi()>0 ? DataEffFit.getErrorHi() : DataEffFit.getError();
  float BinomialError = TMath::Sqrt(SGNtruePass/SGNtrue*(1-SGNtruePass/SGNtrue)/SGNtrue);
 

  Double_t* truthMC = new Double_t[6];
  Double_t* tnpMC   = new Double_t[6];
  Double_t* tnpData = new Double_t[6];

  truthMC[0] = binCenter_;
  truthMC[1] = binWidth_;
  truthMC[2] = binWidth_;
  truthMC[3] = SGNtruePass/SGNtrue;
  truthMC[4] = BinomialError;
  truthMC[5] = BinomialError;
  if(makeSoupFit_){
    tnpMC[0] = binCenter_;
    tnpMC[1] = binWidth_;
    tnpMC[2] = binWidth_;
    tnpMC[3] = McEffFit.getVal();
    tnpMC[4] = (-1)*McErrorLo;
    tnpMC[5] = McErrorHi;
  }
  tnpData[0] = binCenter_;
  tnpData[1] = binWidth_;
  tnpData[2] = binWidth_;
  tnpData[3] = DataEffFit.getVal();
  tnpData[4] = (-1)*DataErrorLo;
  tnpData[5] = DataErrorHi;

  out.push_back(truthMC);
  out.push_back(tnpData);
  if(makeSoupFit_) out.push_back(tnpMC);

  test->Close();

  //delete c; delete c2;

  if(verbose_) cout << "returning from bin " << bin_ << endl;
  return out;

}


void makePlot(const string tnp_ = "etoTauMargLooseNoCracks70",
	      const string category_ = "tauAntiEMVA",
	      const string var_ = "abseta",
	      const float xLow_=65, 
	      const float xHigh_=120,
	      bool makeSoupFit_ = false,
	      bool SumW2_ = true,
	      bool verbose_ = true,
	      double bin1_ = 0.0,
	      double bin2_ = 1.5,
	      double bin3_ = 2.5
	      ){

  // output file
  TFile *outFile = new TFile( Form("EtoTauPlots_%s_%s_%s_v5.root",tnp_.c_str(),category_.c_str(),var_.c_str()),"RECREATE");

  cout << "******************** bin1" << endl;
  vector<Double_t*> bin1Results = simFit(false, tnp_, category_,string(Form("%s>%f && %s<%f",var_.c_str(),bin1_,var_.c_str(),bin2_)),(bin2_+bin1_)/2,(bin2_-bin1_)/2 ,xLow_, xHigh_, SumW2_, verbose_);
  cout << "******************** bin2" << endl;
  vector<Double_t*> bin2Results = simFit(false, tnp_, category_,string(Form("(%s > %f && %s<%f)",var_.c_str(),bin2_,var_.c_str(),bin3_)),(bin3_+bin2_)/2,(bin3_-bin3_)/2 ,xLow_, xHigh_, SumW2_, verbose_);
  
  Double_t truthMC_x[2]  = {(bin1Results[0])[0],(bin2Results[0])[0]};
  Double_t truthMC_xL[2] = {(bin1Results[0])[1],(bin2Results[0])[1]};
  Double_t truthMC_xH[2] = {(bin1Results[0])[2],(bin2Results[0])[2]};
  Double_t truthMC_y[2]  = {(bin1Results[0])[3],(bin2Results[0])[3]};
  Double_t truthMC_yL[2] = {(bin1Results[0])[4],(bin2Results[0])[4]};
  Double_t truthMC_yH[2] = {(bin1Results[0])[5],(bin2Results[0])[5]};
  //
  Double_t tnpMC_x[2]  = {0,0};
  Double_t tnpMC_xL[2] = {0,0};
  Double_t tnpMC_xH[2] = {0,0};
  Double_t tnpMC_y[2]  = {0,0};
  Double_t tnpMC_yL[2] = {0,0};
  Double_t tnpMC_yH[2] = {0,0};
  if(makeSoupFit_){
    tnpMC_x[0]  = (bin1Results[2])[0]; tnpMC_x[1]  = (bin2Results[2])[0];
    tnpMC_xL[0] = (bin1Results[2])[1]; tnpMC_xL[1] = (bin2Results[2])[1];
    tnpMC_xH[0] = (bin1Results[2])[2]; tnpMC_xH[1] = (bin2Results[2])[2];
    tnpMC_y[0]  = (bin1Results[2])[3]; tnpMC_y[1]  = (bin2Results[2])[3];
    tnpMC_yL[0] = (bin1Results[2])[4]; tnpMC_yL[1] = (bin2Results[2])[4];
    tnpMC_yH[0] = (bin1Results[2])[5]; tnpMC_yH[1] = (bin2Results[2])[5];
  }
  //
  Double_t tnpDATA_x[2]  = {(bin1Results[1])[0],(bin2Results[1])[0]};
  Double_t tnpDATA_xL[2] = {(bin1Results[1])[1],(bin2Results[1])[1]};
  Double_t tnpDATA_xH[2] = {(bin1Results[1])[2],(bin2Results[1])[2]};
  Double_t tnpDATA_y[2]  = {(bin1Results[1])[3],(bin2Results[1])[3]};
  Double_t tnpDATA_yL[2] = {(bin1Results[1])[4],(bin2Results[1])[4]};
  Double_t tnpDATA_yH[2] = {(bin1Results[1])[5],(bin2Results[1])[5]};
  
  TCanvas *c1 = new TCanvas("effCanvas","Efficiency canvas",10,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.4,0.6,0.89,0.89,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);

  double binsEdges[3] = {bin1_,bin2_,bin3_};
  TH1F* h1 = new TH1F("h1","Passing efficiency plot",2, binsEdges);
  h1->SetXTitle( var_.c_str() );
  h1->SetYTitle( ("Efficiency of passing "+category_).c_str() );
  
  TGraphAsymmErrors* graph_truthMC = new TGraphAsymmErrors(2,truthMC_x,truthMC_y, truthMC_xL,truthMC_xH,truthMC_yL,truthMC_yH);
  graph_truthMC->SetMarkerColor(kBlue);
  graph_truthMC->SetMarkerStyle(20);
  graph_truthMC->SetMarkerSize(1);
  ///////////////////////////////////////////////////////
  TGraphAsymmErrors* graph_tnpMC = new TGraphAsymmErrors(2,tnpMC_x,tnpMC_y,tnpMC_xL,tnpMC_xH,tnpMC_yL,tnpMC_yH);
  graph_tnpMC->SetMarkerColor(kRed);
  graph_tnpMC->SetMarkerStyle(21);
  graph_tnpMC->SetMarkerSize(1);
  ///////////////////////////////////////////////////////
  TGraphAsymmErrors* graph_tnpDATA = new TGraphAsymmErrors(2,tnpDATA_x, tnpDATA_y,tnpDATA_xL,tnpDATA_xH,tnpDATA_yL,tnpDATA_yH);
  graph_tnpDATA->SetMarkerColor(kBlack);
  graph_tnpDATA->SetMarkerStyle(22);
  graph_tnpDATA->SetMarkerSize(1);

  c1->cd();
  gPad->SetLeftMargin(0.15); 
  h1->GetYaxis()->SetTitleOffset(1.4);
  h1->Draw("");
  graph_truthMC->Draw("PSAME");
  if(makeSoupFit_) graph_tnpMC->Draw("PSAME");
  graph_tnpDATA->Draw("PSAME");
  TH1F* h1mcTruth = new TH1F("h1mcTruth","",1,0,1);
  h1mcTruth->SetLineColor(kBlue);
  h1mcTruth->SetMarkerStyle(20);
  h1mcTruth->SetMarkerSize(1);
  TH1F* h1tnpMC = new TH1F("h1tnpMC","",1,0,1);
  h1tnpMC->SetLineColor(kRed);
  h1tnpMC->SetMarkerStyle(21);
  h1tnpMC->SetMarkerSize(1);
  TH1F* h1tnpDATA = new TH1F("h1tnpDATA","",1,0,1);
  h1tnpDATA->SetLineColor(kBlack);
  h1tnpDATA->SetMarkerStyle(22);
  h1tnpDATA->SetMarkerSize(1);
  leg->SetHeader(Form("#splitline{L=32 pb^{-1} Run2010A/B}{%s passing %s}",tnp_.c_str(),category_.c_str()));
  leg->AddEntry(h1mcTruth,"mc-truth");
  if(makeSoupFit_) leg->AddEntry(h1tnpMC,"t&p: simulation");
  leg->AddEntry(h1tnpDATA,"t&p: 7 TeV Data");
  leg->Draw();

  outFile->cd();
  c1->Write();
  outFile->Close();

}


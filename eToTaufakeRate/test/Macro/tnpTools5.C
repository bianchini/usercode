
///////// 
//source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.26.00/slc4_ia32_gcc34/root/bin/thisroot.(c)sh
////////

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
#include "RooTruthModel.h"
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
  TFile *test = new TFile( Form("EtoTauPlotsFit_%s_%s_%f.root",tnp_.c_str(),category_.c_str(),binCenter_),"RECREATE");
  test->mkdir(Form("bin%f",binCenter_));

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
  TFile fsup("/data_CMS/cms/lbianchini/tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup.root");
  TFile fbkg("/data_CMS/cms/lbianchini/tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_bkg.root");
  TFile fsgn("/data_CMS/cms/lbianchini/tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_sgn.root");
  TFile fdat("/data_CMS/cms/lbianchini/tagAndProbe/trees/38XWcut/testNewWriteFromPAT_Data.root");
  // data from 2iter:
  //TFile fdat("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Data.root");
  
  //********************** signal only tree *************************/

  TTree *fullTreeSgn = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());
  TH1F* hSall        = new TH1F("hSall","",1,0,150);
  TH1F* hSPall       = new TH1F("hSPall","",1,0,150);
  TH1F* hS           = new TH1F("hS","",1,0,150);
  TH1F* hSP          = new TH1F("hSP","",1,0,150);
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
 
  float BKG           = hB->Integral();
  float BKGUnWeighted = hB->GetEntries();
  
  fullTreeBkg->Draw("mass>>hBP",Form("weight*(%s && %s>0 && mass>%f && mass<%f && signalPFChargedHadrCands<1.5 )",bin_.c_str(),category_.c_str(),xLow_,xHigh_));
  
  float BKGPass           = hBP->Integral();
  float BKGUnWeightedPass = hBP->GetEntries();
  float BKGFail           = BKG-BKGPass;
  cout << "*********** BKGFail " << BKGFail << endl;

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
  w->factory("RooExponential::DataBackgroundPdfP(mass,DataCP[0,-10,10])");
  // background fail pdf for Data
  w->factory("RooExponential::DataBackgroundPdfF(mass,DataCF[0,-10,10])");
  // fit parameters for background
  w->factory("McEfficiency[0.04,0,1]");
  w->factory("McNumSgn[0,1000000]");
  w->factory("McNumBkgP[0,100000]");
  w->factory("McNumBkgF[0,100000]"); 
  w->factory("expr::McNumSgnP('McEfficiency*McNumSgn',McEfficiency,McNumSgn)");
  w->factory("expr::McNumSgnF('(1-McEfficiency)*McNumSgn',McEfficiency,McNumSgn)");
  w->factory("McPassing[pass=1,fail=0]");
  // fit parameters for data
  w->factory("DataEfficiency[0.1,0,1]");
  w->factory("DataNumSgn[0,1000000]");
  w->factory("DataNumBkgP[0,1000000]");
  w->factory("DataNumBkgF[0,10000]");
  w->factory("expr::DataNumSgnP('DataEfficiency*DataNumSgn',DataEfficiency,DataNumSgn)");
  w->factory("expr::DataNumSgnF('(1-DataEfficiency)*DataNumSgn',DataEfficiency,DataNumSgn)");
  w->factory("DataPassing[pass=1,fail=0]");

  RooRealVar  *weight = w->var("weight");
  RooRealVar  *abseta = w->var("abseta");
  RooRealVar  *pt     = w->var("pt");
  RooRealVar  *mass   = w->var("mass");
  mass->setRange(xLow_,xHigh_);
  RooRealVar  *mcTrue = w->var("mcTrue");
  RooRealVar  *cut    = w->var( category_.c_str() );
  RooRealVar  *signalPFChargedHadrCands = w->var("signalPFChargedHadrCands");
 
  // build the template for the signal pass sample:
  RooDataSet templateP("templateP","dataset for signal-pass template", RooArgSet(*mass,*weight,*abseta,*pt,*cut,*mcTrue,*signalPFChargedHadrCands), Import( *fullTreeSgn ), /*WeightVar( *weight ),*/ Cut( Form("(mcTrue && %s>0.5 && %s && signalPFChargedHadrCands<1.5)",category_.c_str(),bin_.c_str()) ) );
  // build the template for the signal fail sample:
  RooDataSet templateF("templateF","dataset for signal-fail template", RooArgSet(*mass,*weight,*abseta,*pt,*cut,*mcTrue,*signalPFChargedHadrCands), Import( *fullTreeSgn ), /*WeightVar( *weight ),*/ Cut( Form("(mcTrue && %s<0.5 && %s && signalPFChargedHadrCands<1.5)",category_.c_str(),bin_.c_str()) ) );
  

  mass->setBins(24);
  RooDataHist templateHistP("templateHistP","",RooArgSet(*mass), templateP, 1.0);
  RooHistPdf TemplateSignalPdfP("TemplateSignalPdfP","",RooArgSet(*mass),templateHistP);
  w->import(TemplateSignalPdfP);

  mass->setBins(24);
  RooDataHist templateHistF("templateHistF","",RooArgSet(*mass),templateF,1.0);
  RooHistPdf TemplateSignalPdfF("TemplateSignalPdfF","",RooArgSet(*mass),templateHistF);
  w->import(TemplateSignalPdfF);

  mass->setBins(10000,"fft");

  RooPlot* TemplateFrameP = mass->frame(Bins(24),Title("Template passing"));
  templateP.plotOn(TemplateFrameP);
  w->pdf("TemplateSignalPdfP")->plotOn(TemplateFrameP);
  
  RooPlot* TemplateFrameF = mass->frame(Bins(24),Title("Template failing"));
  templateF.plotOn(TemplateFrameF);
  w->pdf("TemplateSignalPdfF")->plotOn(TemplateFrameF);

  //w->factory("RooFFTConvPdf::McSignalPdfP(mass,TemplateSignalPdfP,RooTruthModel::McResolModP(mass))");
  //w->factory("RooFFTConvPdf::McSignalPdfF(mass,TemplateSignalPdfF,RooTruthModel::McResolModF(mass))");

  // FOR GREGORY: PROBLEM WHEN TRY TO USE THE PURE TEMPLATE =>
  RooHistPdf McSignalPdfP("McSignalPdfP","McSignalPdfP",RooArgSet(*mass),templateHistP);
  RooHistPdf McSignalPdfF("McSignalPdfF","McSignalPdfF",RooArgSet(*mass),templateHistF);
  w->import(McSignalPdfP);
  w->import(McSignalPdfF);
  // FOR GREGORY: FOR DATA, CONVOLUTION IS OK =>
  w->factory("RooFFTConvPdf::DataSignalPdfP(mass,TemplateSignalPdfP,RooGaussian::DataResolModP(mass,DataMeanResP[0.0,-5.,5.],DataSigmaResP[0.5,0.,10]))");
  w->factory("RooFFTConvPdf::DataSignalPdfF(mass,TemplateSignalPdfF,RooGaussian::DataResolModF(mass,DataMeanResF[-5.,-10.,10.],DataSigmaResF[0.5,0.,10]))");
  //w->factory("RooCBShape::DataSignalPdfF(mass,DataMeanF[91.2,88,95.],DataSigmaF[3,0.5,8],DataAlfaF[1.8,0.,10],DataNF[1.0,1e-06,10])");
  //w->factory("RooFFTConvPdf::DataSignalPdfF(mass,RooVoigtian::DataVoigF(mass,DataMeanF[85,80,95],DataWidthF[2.49],DataSigmaF[3,0.5,10]),RooCBShape::DataResolModF(mass,DataMeanResF[0.5,0.,10.],DataSigmaResF[0.5,0.,10],DataAlphaResF[0.5,0.,10],DataNResF[1.0,1e-06,10]))");
  //w->factory("SUM::DataSignalPdfF(fVBP[0.5,0,1]*RooBifurGauss::bifF(mass,DataMeanResF[91.2,80,95],sigmaLF[10,0.5,40],sigmaRF[0.]), RooVoigtian::voigF(mass, DataMeanResF, widthF[2.49], sigmaVoigF[5,0.1,10]) )" );
  
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
  TFile *f = new TFile("dummySoup.root","RECREATE");
  TTree* cutTreeSoupP = fullTreeSoup->CopyTree(Form("(%s>0.5 && %s && signalPFChargedHadrCands<1.5)",category_.c_str(),bin_.c_str()));
  TTree* cutTreeSoupF = fullTreeSoup->CopyTree(Form("(%s<0.5 && %s && signalPFChargedHadrCands<1.5)",category_.c_str(),bin_.c_str()));
 
  RooDataSet McDataP("McDataP","dataset pass for the soup", RooArgSet(*mass), Import( *cutTreeSoupP ) );
 
  RooDataSet McDataF("McDataF","dataset fail for the soup", RooArgSet(*mass), Import( *cutTreeSoupF ) );
 
  RooDataHist McCombData("McCombData","combined data for the soup", RooArgSet(*mass), Index(*(w->cat("McPassing"))), Import("pass", *(McDataP.createHistogram("histoP",*mass)) ), Import("fail",*(McDataF.createHistogram("histoF",*mass)) ) ) ;

  RooPlot* McFrameP    = 0;
  RooPlot* McFrameF    = 0;
  RooRealVar* McEffFit = 0;

  if(makeSoupFit_){

    cout << "**************** N bins in mass " << w->var("mass")->getBins() << endl;

    RooFitResult* ResMcCombinedFit = w->pdf("McModel")->fitTo(McCombData, Extended(1), Minos(1), Save(1),  SumW2Error( SumW2_ ), Range(xLow_,xHigh_), NumCPU(4) /*, ExternalConstraints( *(w->pdf("ConstrainMcNumBkgF")) )*/ );
    test->cd(Form("bin%f",binCenter_));
    ResMcCombinedFit->Write("McFitResults_Combined");

    RooArgSet McFitParam(ResMcCombinedFit->floatParsFinal());
    McEffFit     = (RooRealVar*)(&McFitParam["McEfficiency"]);
    RooRealVar* McNumSigFit  = (RooRealVar*)(&McFitParam["McNumSgn"]);
    RooRealVar* McNumBkgPFit = (RooRealVar*)(&McFitParam["McNumBkgP"]);
    RooRealVar* McNumBkgFFit = (RooRealVar*)(&McFitParam["McNumBkgF"]);

    McFrameP = mass->frame(Bins(24),Title("MC: passing sample"));
    McCombData.plotOn(McFrameP,Cut("McPassing==McPassing::pass"));
    w->pdf("McModel")->plotOn(McFrameP,Slice(*(w->cat("McPassing")),"pass"), ProjWData(*(w->cat("McPassing")),McCombData), LineColor(kBlue),Range(xLow_,xHigh_));
    w->pdf("McModel")->plotOn(McFrameP,Slice(*(w->cat("McPassing")),"pass"), ProjWData(*(w->cat("McPassing")),McCombData), Components("McSignalPdfP"), LineColor(kRed),Range(xLow_,xHigh_));
    w->pdf("McModel")->plotOn(McFrameP,Slice(*(w->cat("McPassing")),"pass"), ProjWData(*(w->cat("McPassing")),McCombData), Components("McBackgroundPdfP"), LineColor(kGreen),Range(xLow_,xHigh_));
    
    McFrameF = mass->frame(Bins(24),Title("MC: failing sample"));
    McCombData.plotOn(McFrameF,Cut("McPassing==McPassing::fail"));
    w->pdf("McModel")->plotOn(McFrameF,Slice(*(w->cat("McPassing")),"fail"), ProjWData(*(w->cat("McPassing")),McCombData), LineColor(kBlue),Range(xLow_,xHigh_));
    w->pdf("McModel")->plotOn(McFrameF,Slice(*(w->cat("McPassing")),"fail"), ProjWData(*(w->cat("McPassing")),McCombData), Components("McSignalPdfF"), LineColor(kRed),Range(xLow_,xHigh_)); 
    w->pdf("McModel")->plotOn(McFrameF,Slice(*(w->cat("McPassing")),"fail"), ProjWData(*(w->cat("McPassing")),McCombData), Components("McBackgroundPdfF"), LineColor(kGreen),Range(xLow_,xHigh_)); 
  }
  
  ///////////////////////////////////////////////////////////////

  /****************** sim fit to data **************************/

  ///////////////////////////////////////////////////////////////
  TFile *f2 = new TFile("dummyData.root","RECREATE");
  TTree* cutTreeDataP = fullTreeData->CopyTree(Form("(%s>0.5 && %s && signalPFChargedHadrCands<1.5)",category_.c_str(),bin_.c_str()));
  TTree* cutTreeDataF = fullTreeData->CopyTree(Form("(%s<0.5 && %s && signalPFChargedHadrCands<1.5)",category_.c_str(),bin_.c_str()));
 
  RooDataSet DataDataP("DataDataP","dataset pass for the soup", RooArgSet(*mass), Import( *cutTreeDataP ) );
  RooDataSet DataDataF("DataDataF","dataset fail for the soup", RooArgSet(*mass), Import( *cutTreeDataF ) );
  RooDataHist DataCombData("DataCombData","combined data for the soup", RooArgSet(*mass), Index(*(w->cat("DataPassing"))), Import("pass",*(DataDataP.createHistogram("histoDataP",*mass))),Import("fail",*(DataDataF.createHistogram("histoDataF",*mass)))) ;

  RooFitResult* ResDataCombinedFit = w->pdf("DataModel")->fitTo(DataCombData, Extended(1), Minos(1), Save(1),  SumW2Error( SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));
  test->cd(Form("bin%f",binCenter_));
  ResDataCombinedFit->Write("DataFitResults_Combined");

  RooArgSet DataFitParam(ResDataCombinedFit->floatParsFinal());
  RooRealVar* DataEffFit     = (RooRealVar*)(&DataFitParam["DataEfficiency"]);
  RooRealVar* DataNumSigFit  = (RooRealVar*)(&DataFitParam["DataNumSgn"]);
  RooRealVar* DataNumBkgPFit = (RooRealVar*)(&DataFitParam["DataNumBkgP"]);
  RooRealVar* DataNumBkgFFit = (RooRealVar*)(&DataFitParam["DataNumBkgF"]);

  RooPlot* DataFrameP = mass->frame(Bins(24),Title("Data: passing sample"));
  DataCombData.plotOn(DataFrameP,Cut("DataPassing==DataPassing::pass"));
  w->pdf("DataModel")->plotOn(DataFrameP,Slice(*(w->cat("DataPassing")),"pass"), ProjWData(*(w->cat("DataPassing")),DataCombData), LineColor(kBlue),Range(xLow_,xHigh_));
  w->pdf("DataModel")->plotOn(DataFrameP,Slice(*(w->cat("DataPassing")),"pass"), ProjWData(*(w->cat("DataPassing")),DataCombData), Components("DataSignalPdfP"), LineColor(kRed),Range(xLow_,xHigh_));
  w->pdf("DataModel")->plotOn(DataFrameP,Slice(*(w->cat("DataPassing")),"pass"), ProjWData(*(w->cat("DataPassing")),DataCombData), Components("DataBackgroundPdfP"), LineColor(kGreen),LineStyle(kDashed),Range(xLow_,xHigh_));
  
  RooPlot* DataFrameF = mass->frame(Bins(24),Title("Data: failing sample"));
  DataCombData.plotOn(DataFrameF,Cut("DataPassing==DataPassing::fail"));
  w->pdf("DataModel")->plotOn(DataFrameF,Slice(*(w->cat("DataPassing")),"fail"), ProjWData(*(w->cat("DataPassing")),DataCombData), LineColor(kBlue),Range(xLow_,xHigh_));
  w->pdf("DataModel")->plotOn(DataFrameF,Slice(*(w->cat("DataPassing")),"fail"), ProjWData(*(w->cat("DataPassing")),DataCombData), Components("DataSignalPdfF"), LineColor(kRed),Range(xLow_,xHigh_));
  w->pdf("DataModel")->plotOn(DataFrameF,Slice(*(w->cat("DataPassing")),"fail"), ProjWData(*(w->cat("DataPassing")),DataCombData), Components("DataBackgroundPdfF"), LineColor(kGreen),LineStyle(kDashed),Range(xLow_,xHigh_));
  ///////////////////////////////////////////////////////////////

 
  if(makeSoupFit_) c->Divide(2,2);
  else c->Divide(2,1);
 
  c->cd(1);
  DataFrameP->Draw();
  c->cd(2);
  DataFrameF->Draw();

  if(makeSoupFit_){
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
    McErrorLo = McEffFit->getErrorLo()<0 ? McEffFit->getErrorLo() : (-1)*McEffFit->getError();
    McErrorHi = McEffFit->getErrorHi()>0 ? McEffFit->getErrorHi() : McEffFit->getError();
  }
  float DataErrorLo = DataEffFit->getErrorLo()<0 ? DataEffFit->getErrorLo() : (-1)*DataEffFit->getError();
  float DataErrorHi = DataEffFit->getErrorHi()>0 ? DataEffFit->getErrorHi() : DataEffFit->getError();
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
    tnpMC[3] = McEffFit->getVal();
    tnpMC[4] = (-1)*McErrorLo;
    tnpMC[5] = McErrorHi;
  }
  tnpData[0] = binCenter_;
  tnpData[1] = binWidth_;
  tnpData[2] = binWidth_;
  tnpData[3] = DataEffFit->getVal();
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
	      const float xHigh_=113,
	      bool makeSoupFit_ = true,
	      bool SumW2_ = false,
	      bool verbose_ = true,
	      double bin1_ = 0.0,
	      double bin2_ = 1.5,
	      double bin3_ = 2.5
	      ){

  // output file
  TFile *outFile = new TFile( Form("EtoTauPlots_%s_%s_%s_v5.root",tnp_.c_str(),category_.c_str(),var_.c_str()),"RECREATE");
  cout << "******************** bin1" << endl;
  vector<Double_t*> bin1Results = simFit(makeSoupFit_, tnp_, category_,string(Form("%s>%f && %s<%f",var_.c_str(),bin1_,var_.c_str(),bin2_)),(bin2_+bin1_)/2,(bin2_-bin1_)/2 ,xLow_, xHigh_, SumW2_, verbose_);
  cout << "******************** bin2" << endl;
  vector<Double_t*> bin2Results = simFit(makeSoupFit_, tnp_, category_,string(Form("(%s > %f && %s<%f)",var_.c_str(),bin2_,var_.c_str(),bin3_)),(bin3_+bin2_)/2,(bin3_-bin3_)/2 ,xLow_, xHigh_, SumW2_, verbose_);
  
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


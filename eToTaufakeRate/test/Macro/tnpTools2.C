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
  TFile *test = new TFile( Form("EtoTauPlotsFit_%s_%s_%f_v3.root",tnp_.c_str(),category_.c_str(),binCenter_),"RECREATE");
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
  TFile fsup("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup.root");
  TFile fbkg("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_bkg.root");
  TFile fsgn("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_sgn.root");
  TFile fdat("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_Data.root");

  //********************** signal only tree *************************/

  TTree *fullTreeSgn = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());
  TH1F* hSall = new TH1F("hSall","",1,0,150);
  TH1F* hSPall = new TH1F("hSPall","",1,0,150);
  TH1F* hS = new TH1F("hS","",1,0,150);
  TH1F* hSP = new TH1F("hSP","",1,0,150);
  fullTreeSgn->Draw("mass>>hS",Form("weight*(%s && mass>%f && mass<%f && mcTrue)",bin_.c_str(),xLow_,xHigh_));
  fullTreeSgn->Draw("mass>>hSall",Form("weight*(%s && mass>%f && mass<%f)",bin_.c_str(),xLow_,xHigh_));

  float SGNtrue = hS->Integral();
  float SGNall  = hSall->Integral();
 
  fullTreeSgn->Draw("mass>>hSP",Form("weight*(%s && %s>0 && mass>%f && mass<%f && mcTrue)",bin_.c_str(),category_.c_str(),xLow_,xHigh_));
  fullTreeSgn->Draw("mass>>hSPall",Form("weight*(%s && %s>0 && mass>%f && mass<%f)",bin_.c_str(),category_.c_str(),xLow_,xHigh_));

  float SGNtruePass = hSP->Integral();
  float SGNallPass  = hSPall->Integral();

  //********************** background only tree *************************//

  TTree *fullTreeBkg = (TTree*)fbkg.Get((tnp_+"/fitter_tree").c_str());
  TH1F* hB = new TH1F("hB","",1,0,150);
  TH1F* hBP = new TH1F("hBP","",1,0,150);
  fullTreeBkg->Draw("mass>>hB",Form("weight*(%s && mass>%f && mass<%f)",bin_.c_str(),xLow_,xHigh_));
 
  float BKG = hB->Integral();
  float BKGUnWeighted = hB->GetEntries();
  
  fullTreeBkg->Draw("mass>>hBP",Form("weight*(%s && %s>0 && mass>%f && mass<%f)",bin_.c_str(),category_.c_str(),xLow_,xHigh_));
  
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
   w->factory((category_+"[0,1]").c_str());
  // template for passing signal pdf: to be extracted from simulation
  //w->factory("SUM::TemplateSignalPdfP(TemFVBP[0.5,0,1]*RooBifurGauss::TemBifP(mass,TemMeanP[91.2,88,95],TemSigmaLP[10,0.5,40],TemSigmaRP[5.,0.5,15.]), RooVoigtian::TemVoigP(mass, TemMeanP, TemWidthP[2.49], TemSigmaVoigP[5,0.1,10]) )");
  // signal fail pdf for simulation
  w->factory("RooCBShape::McSignalPdfF(mass,MCmeanF[91.2,88,95.],MCsigmaF[3,0.5,8],MCalfaF[1.8,0.,10],MCnF[1.0,1e-06,10])");
  //w->factory("RooGaussian::McSignalPdfF(mass,MCmeanF[91.2,88,95.],MCsigmaF[3,0.5,8])");
  // signal fail pdf for data
  w->factory("RooCBShape::DataSignalPdfF(mass,DataMeanF[91.2,88,95.],DataSigmaF[3,0.5,8],DataAlfaF[1.8,0.,10],DataNF[1.0,1e-06,10])");
  //w->factory("RooGaussian::DataSignalPdfF(mass,DataMeanF[91.2,88,95.],DataSigmaF[3,0.5,8])");
  // background pass pdf for MC
  w->factory("RooExponential::McBackgroundPdfP(mass,McCP[0,-10,10])");
  // background fail pdf for MC
  w->factory("RooExponential::McBackgroundPdfF(mass,McCF[0,-10,10])");
  // background pass pdf for Data
  w->factory("RooExponential::DataBackgroundPdfP(mass,DataCP[0,-10,10])");
  // background fail pdf for Data
  w->factory("RooExponential::DataBackgroundPdfF(mass,DataCF[0,-10,10])");
  // fit parameters for background
  w->factory("McEfficiency[0.1,0,1]");
  w->factory("McNumSgn[0,100000]");
  w->factory("McNumBkgP[0,100000]");
  w->factory("McNumBkgF[0,10]"); 
  w->factory("expr::McNumSgnP('McEfficiency*McNumSgn',McEfficiency,McNumSgn)");
  w->factory("expr::McNumSgnF('(1-McEfficiency)*McNumSgn',McEfficiency,McNumSgn)");
  w->factory("McPassing[pass=1,fail=0]");
  // fit parameters for data
  w->factory("DataEfficiency[0.1,0,1]");
  w->factory("DataNumSgn[0,1000000]");
  w->factory("DataNumBkgP[0,1000000]");
  w->factory("DataNumBkgF[0,10]");
  w->factory("expr::DataNumSgnP('DataEfficiency*DataNumSgn',DataEfficiency,DataNumSgn)");
  w->factory("expr::DataNumSgnF('(1-DataEfficiency)*DataNumSgn',DataEfficiency,DataNumSgn)");
  w->factory("DataPassing[pass=1,fail=0]");

  RooRealVar  *weight = w->var("weight");
  RooRealVar  *abseta = w->var("abseta");
  RooRealVar  *pt = w->var("pt");
  RooRealVar  *mass = w->var("mass");
  mass->setRange(xLow_,xHigh_);
  RooRealVar  *mcTrue = w->var("mcTrue");
  RooRealVar  *cut = w->var( category_.c_str() );
 
  // build the template for the signal pass sample:
  RooDataSet templateP("templateP","dataset for signal-pass template", RooArgSet(*mass,*weight,*abseta,*pt,*cut,*mcTrue), Import( *fullTreeSgn ), /*WeightVar( *weight ),*/ Cut( Form("(mcTrue && %s>0.5 && %s)",category_.c_str(),bin_.c_str()) ) );
  
  mass->setMin(60); mass->setMax(120); mass->setBins(30);
  RooDataHist templateHistP("templateHistP","",RooArgSet(*mass),templateP,1.0);
  RooHistPdf TemplateSignalPdfP("TemplateSignalPdfP","",RooArgSet(*mass),templateHistP);
  w->import(TemplateSignalPdfP);

  mass->setBins(10000,"fft");
  //RooFitResult* ResTemplateSignalP =  w->pdf("TemplateSignalPdfP")->fitTo(templateP,Minos(1), Save(1), SumW2Error( SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));

  RooPlot* TemplateFrameP = mass->frame(Title("Template passing"));
  templateP.plotOn(TemplateFrameP);
  w->pdf("TemplateSignalPdfP")->plotOn(TemplateFrameP);

  //RooArgSet templateFitParam(ResTemplateSignalP->floatParsFinal()); 
  //RooRealVar *TemFVBP       = (RooRealVar*)(&templateFitParam["TemFVBP"]);
  //RooRealVar *TemMeanP      = (RooRealVar*)(&templateFitParam["TemMeanP"]);
  //RooRealVar *TemSigmaLP    = (RooRealVar*)(&templateFitParam["TemSigmaLP"]);
  //RooRealVar *TemSigmaRP    = (RooRealVar*)(&templateFitParam["TemSigmaRP"]);
  //RooRealVar *TemSigmaVoigP = (RooRealVar*)(&templateFitParam["TemSigmaVoigP"]);

  // signal pass pdf for simulation
  //w->factory( Form("SUM::McSignalPdfP(McFVBP[%f]*RooBifurGauss::McBifP(mass,McMeanP[%f,85,95],McSigmaLP[%f],McSigmaRP[%f]), RooVoigtian::McVoigP(mass, McMeanP, McWidthP[2.49], McSigmaVoigP[%f]) )",TemFVBP->getVal(),TemMeanP->getVal(),TemSigmaLP->getVal(), TemSigmaRP->getVal(),TemSigmaVoigP->getVal()));
  //w->factory("RooGaussian::McSignalPdfP(mass,McMeanP[91.2,88,95.],McSigmaP[3,0.5,8])");

  w->factory("RooFFTConvPdf::McSignalPdfP(mass,TemplateSignalPdfP,RooGaussian::McResolMod(mass,McMeanRes[0.5,0.,10.],McSigmaRes[0.5,0.,10]))");
  w->factory("RooFFTConvPdf::DataSignalPdfP(mass,TemplateSignalPdfP,RooGaussian::DataResolMod(mass,DataMeanRes[0.5,0.,10.],DataSigmaRes[0.5,0.,10]))");

  // signal pass pdf for data
  //w->factory( Form("SUM::DataSignalPdfP(DataFVBP[%f]*RooBifurGauss::DataBifP(mass,DataMeanP[%f,85,95],DataSigmaLP[%f],DataSigmaRP[%f]), RooVoigtian::DataVoigP(mass, DataMeanP, DataWidthP[2.49], DataSigmaVoigP[%f,0.5,10]) )",TemFVBP->getVal(),TemMeanP->getVal(),TemSigmaLP->getVal(), TemSigmaRP->getVal(),TemSigmaVoigP->getVal()));
  //w->factory("RooGaussian::DataSignalPdfP(mass,DataMeanP[91.2,88,95.],DataSigmaP[3,0.5,8])");
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
  RooDataSet McDataP("McDataP","dataset pass for the soup", RooArgSet(*mass,*weight,*abseta,*pt,*cut), Import( *fullTreeSoup ), /*WeightVar( *weight ),*/ Cut( Form("(%s>0.5 && %s)",category_.c_str(),bin_.c_str()) ) );
  RooDataSet McDataF("McDataF","dataset fail for the soup", RooArgSet(*mass,*weight,*abseta,*pt,*cut), Import( *fullTreeSoup ), /*WeightVar( *weight ),*/ Cut( Form("(%s<0.5 && %s)",category_.c_str(),bin_.c_str()) ) );
  RooDataSet McCombData("McCombData","combined data for the soup", RooArgSet(*mass,*weight,*abseta,*pt,*cut), Index(*(w->cat("McPassing"))), Import("pass",McDataP),Import("fail",McDataF)) ;

  RooPlot* McFrameP = 0;
  RooPlot* McFrameF = 0;
  RooRealVar* McEffFit = 0;
  if(makeSoupFit_){
    RooFitResult* ResMcCombinedFit = w->pdf("McModel")->fitTo(McCombData, Extended(1), Minos(1), Save(1),  SumW2Error( SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));
    test->cd(Form("bin%.2f",binCenter_));
    ResMcCombinedFit->Write("McFitResults_Combined");

    RooArgSet McFitParam(ResMcCombinedFit->floatParsFinal());
    McEffFit     = (RooRealVar*)(&McFitParam["McEfficiency"]);
    RooRealVar* McNumSigFit  = (RooRealVar*)(&McFitParam["McNumSgn"]);
    RooRealVar* McNumBkgPFit = (RooRealVar*)(&McFitParam["McNumBkgP"]);
    RooRealVar* McNumBkgFFit = (RooRealVar*)(&McFitParam["McNumBkgF"]);

    McFrameP = mass->frame(Bins(30),Title("MC: passing sample"));
    McCombData.plotOn(McFrameP,Cut("McPassing==McPassing::pass"));
    w->pdf("McModel")->plotOn(McFrameP,Slice(*(w->cat("McPassing")),"pass"), ProjWData(*(w->cat("McPassing")),McCombData), LineColor(kBlue),Range(60,120));
    w->pdf("McModel")->plotOn(McFrameP,Slice(*(w->cat("McPassing")),"pass"), ProjWData(*(w->cat("McPassing")),McCombData), Components("McSignalPdfP"), LineColor(kRed),Range(60,120));
    
    McFrameF = mass->frame(Bins(30),Title("MC: failing sample"));
    McCombData.plotOn(McFrameF,Cut("McPassing==McPassing::fail"));
    w->pdf("McModel")->plotOn(McFrameF,Slice(*(w->cat("McPassing")),"fail"), ProjWData(*(w->cat("McPassing")),McCombData), LineColor(kBlue),Range(60,120));
    w->pdf("McModel")->plotOn(McFrameF,Slice(*(w->cat("McPassing")),"fail"), ProjWData(*(w->cat("McPassing")),McCombData), Components("McSignalPdfF"), LineColor(kRed),Range(60,120)); 
  }
  
  ///////////////////////////////////////////////////////////////

  /****************** sim fit to data **************************/

  ///////////////////////////////////////////////////////////////
  RooDataSet DataDataP("DataDataP","dataset pass for the soup", RooArgSet(*mass,*weight,*abseta,*pt,*cut), Import( *fullTreeData ), /*WeightVar( *weight ),*/ Cut( Form("(%s>0.5 && %s)",category_.c_str(),bin_.c_str()) ) );
  RooDataSet DataDataF("DataDataF","dataset fail for the soup", RooArgSet(*mass,*weight,*abseta,*pt,*cut), Import( *fullTreeData ), /*WeightVar( *weight ),*/ Cut( Form("(%s<0.5 && %s)",category_.c_str(),bin_.c_str()) ) );
  RooDataSet DataCombData("DataCombData","combined data for the soup", RooArgSet(*mass,*weight,*abseta,*pt,*cut), Index(*(w->cat("DataPassing"))), Import("pass",DataDataP),Import("fail",DataDataF)) ;

  RooFitResult* ResDataCombinedFit = w->pdf("DataModel")->fitTo(DataCombData, Extended(1), Minos(1), Save(1),  SumW2Error( SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));
  test->cd(Form("bin%f",binCenter_));
  ResDataCombinedFit->Write("DataFitResults_Combined");

  RooArgSet DataFitParam(ResDataCombinedFit->floatParsFinal());
  RooRealVar* DataEffFit     = (RooRealVar*)(&DataFitParam["DataEfficiency"]);
  RooRealVar* DataNumSigFit  = (RooRealVar*)(&DataFitParam["DataNumSgn"]);
  RooRealVar* DataNumBkgPFit = (RooRealVar*)(&DataFitParam["DataNumBkgP"]);
  RooRealVar* DataNumBkgFFit = (RooRealVar*)(&DataFitParam["DataNumBkgF"]);

  RooPlot* DataFrameP = mass->frame(Bins(30),Title("Data: passing sample"));
  DataCombData.plotOn(DataFrameP,Cut("DataPassing==DataPassing::pass"));
  w->pdf("DataModel")->plotOn(DataFrameP,Slice(*(w->cat("DataPassing")),"pass"), ProjWData(*(w->cat("DataPassing")),DataCombData), LineColor(kBlue),Range(60,120));
  w->pdf("DataModel")->plotOn(DataFrameP,Slice(*(w->cat("DataPassing")),"pass"), ProjWData(*(w->cat("DataPassing")),DataCombData), Components("DataSignalPdfP"), LineColor(kRed),Range(60,120));
  w->pdf("DataModel")->plotOn(DataFrameP,Slice(*(w->cat("DataPassing")),"pass"), ProjWData(*(w->cat("DataPassing")),DataCombData), Components("DataBackgroundPdfP"), LineColor(kGreen),LineStyle(kDashed),Range(60,120));
  
  RooPlot* DataFrameF = mass->frame(Bins(30),Title("Data: failing sample"));
  DataCombData.plotOn(DataFrameF,Cut("DataPassing==DataPassing::fail"));
  w->pdf("DataModel")->plotOn(DataFrameF,Slice(*(w->cat("DataPassing")),"fail"), ProjWData(*(w->cat("DataPassing")),DataCombData), LineColor(kBlue),Range(60,120));
  w->pdf("DataModel")->plotOn(DataFrameF,Slice(*(w->cat("DataPassing")),"fail"), ProjWData(*(w->cat("DataPassing")),DataCombData), Components("DataSignalPdfF"), LineColor(kRed),Range(60,120));
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
 
  c2->cd();
  TemplateFrameP->Draw();
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
	      const float xLow_=60, 
	      const float xHigh_=120,
	      bool makeSoupFit_ = false,
	      bool SumW2_ = true,
	      bool verbose_ = true,
	      double bin1_ = 0.0,
	      double bin2_ = 1.5,
	      double bin3_ = 2.5
	      ){

  // output file
  TFile *outFile = new TFile( Form("EtoTauPlots_%s_%s_%s_v3.root",tnp_.c_str(),category_.c_str(),var_.c_str()),"RECREATE");

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


void DrawFit(const string tnp_ = "etoTauMargLooseNoCracks70", 
	    const string category_ = "tauAntiEMVA",
	    const string bin_ = "abseta<1.5",
	    const float binCenter_ = 0.75,
	    const float binWidth_ = 0.75,
	    const float xLow_=60, 
	    const float xHigh_=120,
	    bool SumW2_ = false,
	    bool verbose_ = true){

  
  TCanvas *c = new TCanvas("c",Form("fitCanvas_%s_%s",tnp_.c_str(),bin_.c_str()),10,30,650,600);
  c->SetGrid(0,0);
  c->SetFillStyle(4000);
  c->SetFillColor(10);
  c->SetTicky();
  c->SetObjectStat(0);

  // output file
  TFile *test = new TFile("myTestFile.root","RECREATE");
  
  // input files
  TFile fsup("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_reduced.root");
  TFile fbkg("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_bkg.root");
  TFile fsgn("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_sgn.root");
  TFile fdat("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_Data.root");

  //********************** signal only tree *************************/

  TTree *fullTreeSgn = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());
  TH1F* hSall = new TH1F("hSall","",1,0,150);
  TH1F* hSPall = new TH1F("hSPall","",1,0,150);
  TH1F* hS = new TH1F("hS","",1,0,150);
  TH1F* hSP = new TH1F("hSP","",1,0,150);
  fullTreeSgn->Draw("mass>>hS",Form("weight*(%s && mass>%f && mass<%f && mcTrue)",bin_.c_str(),xLow_,xHigh_));
  fullTreeSgn->Draw("mass>>hSall",Form("weight*(%s && mass>%f && mass<%f)",bin_.c_str(),xLow_,xHigh_));

  float SGNtrue = hS->Integral();
  float SGNall  = hSall->Integral();
 
  fullTreeSgn->Draw("mass>>hSP",Form("weight*(%s && %s>0 && mass>%f && mass<%f && mcTrue)",bin_.c_str(),category_.c_str(),xLow_,xHigh_));
  fullTreeSgn->Draw("mass>>hSPall",Form("weight*(%s && %s>0 && mass>%f && mass<%f)",bin_.c_str(),category_.c_str(),xLow_,xHigh_));

  float SGNtruePass = hSP->Integral();
  float SGNallPass  = hSPall->Integral();

  //********************** background only tree *************************//

  TTree *fullTreeBkg = (TTree*)fbkg.Get((tnp_+"/fitter_tree").c_str());
  TH1F* hB = new TH1F("hB","",1,0,150);
  TH1F* hBP = new TH1F("hBP","",1,0,150);
  fullTreeBkg->Draw("mass>>hB",Form("weight*(%s && mass>%f && mass<%f)",bin_.c_str(),xLow_,xHigh_));
 
  float BKG = hB->Integral();
  float BKGUnWeighted = hB->GetEntries();
  
  fullTreeBkg->Draw("mass>>hBP",Form("weight*(%s && %s>0 && mass>%f && mass<%f)",bin_.c_str(),category_.c_str(),xLow_,xHigh_));
  
  float BKGPass = hBP->Integral();
  float BKGUnWeightedPass = hBP->GetEntries();

  test->cd();
  //********************** soup tree *************************//

  TTree *fullTreeSoup = (TTree*)fsup.Get((tnp_+"/fitter_tree").c_str());
  TTree* passTreeSoup = (TTree*)fullTreeSoup->CopyTree(Form("(%s>0.5 && %s)",category_.c_str(),bin_.c_str()));
  TTree* failTreeSoup = (TTree*)fullTreeSoup->CopyTree(Form("(%s<0.5 && %s)",category_.c_str(),bin_.c_str()));
  
  //********************** data tree *************************//

  TTree *fullTreeData = (TTree*)fdat.Get((tnp_+"/fitter_tree").c_str());
  TTree* passTreeData = (TTree*)fullTreeData->CopyTree(Form("(%s>0.5 && %s)",category_.c_str(),bin_.c_str()));
  TTree* failTreeData = (TTree*)fullTreeData->CopyTree(Form("(%s<0.5 && %s)",category_.c_str(),bin_.c_str()));

  //********************** workspace ***********************//

  RooWorkspace *w = new RooWorkspace("w","w");
  // tree variables to be imported
  w->factory("mass[30,120]");
  w->factory("weight[0,10000]");
  w->factory("abseta[0,2.5]");
  w->factory("pt[0,200]");
  w->factory("mcTrue[0,1]");
  w->factory((category_+"[0,1]").c_str());
  // template for passing signal pdf: to be extracted from simulation
  //w->factory("SUM::TemplateSignalPdfP(TemFVBP[0.5,0,1]*RooBifurGauss::TemBifP(mass,TemMeanP[91.2,88,95],TemSigmaLP[10,0.5,40],TemSigmaRP[5.,0.5,15.]), RooVoigtian::TemVoigP(mass, TemMeanP, TemWidthP[2.79], TemSigmaVoigP[5,0.1,10]) )");
  //w->factory("RooHistPdf::TemplateSignalPdfP()");
  // signal fail pdf for simulation
  ///w->factory("RooCBShape::McSignalPdfF(mass,MCmeanF[91.2,88,95.],MCsigmaF[3,0.5,8],MCalfaF[1.8,0.,10],MCnF[1.0,1e-06,10])");
  w->factory("RooGaussian::McSignalPdfF(mass,MCmeanF[91.2,88,95.],MCsigmaF[3,0.5,8])");
  // signal fail pdf for data
  //w->factory("RooCBShape::DataSignalPdfF(mass,DataMeanF[91.2,88,95.],DataSigmaF[3,0.5,8],DataAlfaF[1.8,0.,10],DataNF[1.0,1e-06,10])");
  w->factory("RooGaussian::DataSignalPdfF(mass,DataMeanF[91.2,88,95.],DataSigmaF[3,0.5,8])");
  // background pass pdf for MC
  w->factory("RooExponential::McBackgroundPdfP(mass,McCP[0,-10,10])");
  // background fail pdf for MC
  w->factory("RooExponential::McBackgroundPdfF(mass,McCF[0,-10,10])");
  // background pass pdf for Data
  w->factory("RooExponential::DataBackgroundPdfP(mass,DataCP[0,-10,10])");
  // background fail pdf for Data
  w->factory("RooExponential::DataBackgroundPdfF(mass,DataCF[0,-10,10])");
  // fit parameters for background
  w->factory("McEfficiency[0.1,0,1]");
  w->factory("McNumSgn[0,100000]");
  w->factory("McNumBkgP[0,100000]");
  w->factory("McNumBkgF[0,1e-06]"); 
  w->factory("expr::McNumSgnP('McEfficiency*McNumSgn',McEfficiency,McNumSgn)");
  w->factory("expr::McNumSgnF('(1-McEfficiency)*McNumSgn',McEfficiency,McNumSgn)");
  //w->factory("McPassing[pass=1,fail=0]");
  // fit parameters for data
  w->factory("DataEfficiency[0.1,0,1]");
  w->factory("DataNumSgn[0,1000000]");
  w->factory("DataNumBkgP[0,1000000]");
  w->factory("DataNumBkgF[0,1e-06]");
  w->factory("expr::DataNumSgnP('DataEfficiency*DataNumSgn',DataEfficiency,DataNumSgn)");
  w->factory("expr::DataNumSgnF('(1-DataEfficiency)*DataNumSgn',DataEfficiency,DataNumSgn)");
  //w->factory("DataPassing[pass=1,fail=0]");
  w->factory("passing[pass=1,fail=0]");

  RooRealVar  *weight = w->var("weight");
  RooRealVar  *abseta = w->var("abseta");
  RooRealVar  *pt = w->var("pt");
  RooRealVar  *mass = w->var("mass");
  mass->setRange(xLow_,xHigh_);
  RooRealVar  *mcTrue = w->var("mcTrue");
  RooRealVar  *cut = w->var( category_.c_str() );
  RooCategory *category = w->cat("passing");  
 
  // build the template for the signal pass sample:
  RooDataSet templateP("templateP","dataset for signal-pass template", RooArgSet(*mass,*weight,*abseta,*pt,*cut,*mcTrue), Import( *fullTreeSgn ), /*WeightVar( *weight ),*/ Cut( Form("(mcTrue && %s>0.5 && %s)",category_.c_str(),bin_.c_str()) ) );

  mass->setMin(60); mass->setMax(120); mass->setBins(30);
  RooDataHist templateHistP("templateHistP","",RooArgSet(*mass),templateP,1.0);
  RooHistPdf TemplateSignalPdfP("TemplateSignalPdfP","",RooArgSet(*mass),templateHistP);
  w->import(TemplateSignalPdfP);

  mass->setBins(10000,"fft");

  //RooFitResult* ResTemplateSignalP =  w->pdf("TemplateSignalPdfP")->fitTo(templateP,Minos(1), Save(1), SumW2Error( SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));

  //RooArgSet templateFitParam(ResTemplateSignalP->floatParsFinal()); 
  //RooRealVar *TemFVBP       = (RooRealVar*)(&templateFitParam["TemFVBP"]);
  //RooRealVar *TemMeanP      = (RooRealVar*)(&templateFitParam["TemMeanP"]);
  //RooRealVar *TemSigmaLP    = (RooRealVar*)(&templateFitParam["TemSigmaLP"]);
  //RooRealVar *TemSigmaRP    = (RooRealVar*)(&templateFitParam["TemSigmaRP"]);
  //RooRealVar *TemSigmaVoigP = (RooRealVar*)(&templateFitParam["TemSigmaVoigP"]);

  // signal pass pdf for simulation
  //w->factory( Form("SUM::McSignalPdfP(McFVBP[%f]*RooBifurGauss::McBifP(mass,McMeanP[%f,85,95],McSigmaLP[%f],McSigmaRP[%f]), RooVoigtian::McVoigP(mass, McMeanP, McWidthP[2.79], McSigmaVoigP[%f]) )",TemFVBP->getVal(),TemMeanP->getVal(),TemSigmaLP->getVal(), TemSigmaRP->getVal(),TemSigmaVoigP->getVal()));
  w->factory("RooGaussian::McSignalPdfP(mass,McMeanP[91.2,88,95.],McSigmaP[3,0.5,8])");
  // signal pass pdf for data
  //w->factory(Form("DataMeanP[%f]",TemMeanP->getVal()));
  //w->factory(Form("DataFVBP[%f]",TemFVBP->getVal()));
  //w->factory(Form("RooBifurGauss::DataBifP(mass,DataMeanP ,DataSigmaLP[%f],DataSigmaRP[%f])",TemSigmaLP->getVal(), TemSigmaRP->getVal())); 
  //w->factory(Form("RooVoigtian::DataVoigP(mass, DataMeanP, DataWidthP[2.79], DataSigmaVoigP[%f])",TemSigmaVoigP->getVal() ) );
  //w->factory( "SUM::SimulDataSignalPdfP( DataFVBP*DataBifP, DataVoigP");
  //w->factory("RooGaussian::SimulDataSignalPdfP(mass,DataMeanP,DataSigmaP[3,0.5,8])");

  //w->factory( Form("SUM::SimulDataSignalPdfP(DataFVBP[%f,0.0,1.0]*RooBifurGauss::DataBifP(mass,DataMeanP[%f],DataSigmaLP[%f],DataSigmaRP[%f]), RooVoigtian::DataVoigP(mass, DataMeanP, DataWidthP[2.79], DataSigmaVoigP[%f]) )",TemFVBP->getVal(),TemMeanP->getVal(),TemSigmaLP->getVal(), TemSigmaRP->getVal(),TemSigmaVoigP->getVal()));


  w->factory("RooFFTConvPdf::DataSignalPdfP(mass,TemplateSignalPdfP,RooGaussian::resolMod(mass,meanRes[0.,0.,10.],sigmaRes[0,0.,10]))");

  //w->factory("RooGaussian::DataSignalPdfP(mass,DataMeanP[91.2,88,95.],DataSigmaP[3,0.5,8])");
 // composite model pass for MC
  w->factory("SUM::McModelP(McNumSgnP*McSignalPdfP,McNumBkgP*McBackgroundPdfP)");  
  w->factory("SUM::McModelF(McNumSgnF*McSignalPdfF,McNumBkgF*McBackgroundPdfF)"); 
  // composite model pass for data
  w->factory("SUM::DataModelP(DataNumSgnP*DataSignalPdfP,DataNumBkgP*DataBackgroundPdfP)");  
  w->factory("SUM::DataModelF(DataNumSgnF*DataSignalPdfF,DataNumBkgF*DataBackgroundPdfF)");  
  // simultaneous fir for MC
  //w->factory("SIMUL::McModel(McPassing,pass=McModelP,fail=McModelF)");
  w->factory("SIMUL::McModel(passing,pass=McModelP,fail=McModelF)");
  // simultaneous fir for data
  //w->factory("SIMUL::DataModel(DataPassing,pass=DataModelP,fail=DataModelF)");
  w->factory("SIMUL::DataModel(passing,pass=DataModelP,fail=DataModelF)");
  w->Print("V");
  w->saveSnapshot("clean", w->allVars());

  /****************** sim fit to soup **************************/
  w->loadSnapshot("clean");

  ///////////////////////////////////////////////////////////////
  RooDataSet McDataP("McDataP","dataset pass for the soup", RooArgSet(*mass), Import( *passTreeSoup ) /*, WeightVar( *weight ), Cut( Form("(%s>0.5 && %s)",category_.c_str(),bin_.c_str()) )*/ );
  RooDataSet McDataF("McDataF","dataset fail for the soup", RooArgSet(*mass), Import( *failTreeSoup ) /*, WeightVar( *weight ), Cut( Form("(%s<0.5 && %s)",category_.c_str(),bin_.c_str()) )*/ );
  RooDataSet McCombData("McCombData","combined data for the soup", RooArgSet(*mass), Index(*category), Import("pass",McDataP),Import("fail",McDataF)) ;

  RooFitResult* ResMcCombinedFit = w->pdf("McModel")->fitTo(McCombData, Extended(1), Minos(1), Save(1),  SumW2Error( SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));
  //ResMcCombinedFit->Write("McFitResults_Combined");

  RooArgSet McFitParam(ResMcCombinedFit->floatParsFinal());
  RooRealVar* McEffFit     = (RooRealVar*)(&McFitParam["McEfficiency"]);
  RooRealVar* McNumSigFit  = (RooRealVar*)(&McFitParam["McNumSgn"]);
  RooRealVar* McNumBkgPFit = (RooRealVar*)(&McFitParam["McNumBkgP"]);
  RooRealVar* McNumBkgFFit = (RooRealVar*)(&McFitParam["McNumBkgF"]);
  ///////////////////////////////////////////////////////////////

  /****************** sim fit to data **************************/

  ///////////////////////////////////////////////////////////////
  RooDataSet DataDataP("DataDataP","dataset pass for the soup", RooArgSet(*mass), Import( *passTreeData ) /*, WeightVar( *weight ), Cut( Form("(%s>0.5 && %s)",category_.c_str(),bin_.c_str()) )*/ );
  RooDataSet DataDataF("DataDataF","dataset fail for the soup", RooArgSet(*mass), Import( *failTreeData ) /*,WeightVar( *weight ), Cut( Form("(%s<0.5 && %s)",category_.c_str(),bin_.c_str()) )*/ );
  RooDataSet DataCombData("DataCombData","combined data for the soup", RooArgSet(*mass,*weight,*abseta,*pt,*cut), Index(*category), Import("pass",DataDataP),Import("fail",DataDataF)) ;

  RooFitResult* ResDataCombinedFit = w->pdf("DataModel")->fitTo(DataCombData, Extended(1), Minos(1), Save(1),  SumW2Error( SumW2_ ), Range(xLow_,xHigh_), NumCPU(4));
  //ResDataCombinedFit->Write("DataFitResults_Combined");

  RooArgSet DataFitParam(ResDataCombinedFit->floatParsFinal());
  RooRealVar* DataEffFit     = (RooRealVar*)(&DataFitParam["DataEfficiency"]);
  RooRealVar* DataNumSigFit  = (RooRealVar*)(&DataFitParam["DataNumSgn"]);
  RooRealVar* DataNumBkgPFit = (RooRealVar*)(&DataFitParam["DataNumBkgP"]);
  RooRealVar* DataNumBkgFFit = (RooRealVar*)(&DataFitParam["DataNumBkgF"]);

  RooPlot* TemplateFrameP = mass->frame(Title("Template passing"));
  templateP.plotOn(TemplateFrameP);
  w->pdf("TemplateSignalPdfP")->plotOn(TemplateFrameP);

  RooPlot* McFrameP = mass->frame(Title("MC passing"));
  McCombData.plotOn(McFrameP,Cut("passing==passing::pass"));
  w->pdf("McModel")->plotOn(McFrameP,Slice(*category,"pass"), ProjWData(*category,McCombData), LineColor(kBlue),Range(60,120));
  w->pdf("McModel")->plotOn(McFrameP,Slice(*(w->cat("McPassing")),"pass"), ProjWData(*category,McCombData), Components("McSignalPdfP"), LineColor(kRed),Range(60,120));

  RooPlot* McFrameF = mass->frame(Title("MC failing"));
  McCombData.plotOn(McFrameF,Cut("passing==passing::fail"));
  w->pdf("McModel")->plotOn(McFrameF,Slice(*category,"fail"), ProjWData(*category,McCombData), LineColor(kBlue),Range(60,120));
  w->pdf("McModel")->plotOn(McFrameF,Slice(*category,"fail"), ProjWData(*category,McCombData), Components("McSignalPdfF"), LineColor(kRed),Range(60,120));

  RooPlot* DataFrameP = mass->frame(Title("Data passing"));
  DataCombData.plotOn(DataFrameP,Cut("passing==passing::pass"));
  w->pdf("DataModel")->plotOn(DataFrameP,Slice(*category,"pass"), ProjWData(*category,DataCombData), LineColor(kBlue),Range(60,120));
  w->pdf("DataModel")->plotOn(DataFrameP,Slice(*category,"pass"), ProjWData(*category,DataCombData), Components("DataSignalPdfP"), LineColor(kRed),Range(60,120));
  
  RooPlot* DataFrameF = mass->frame(Title("Data failing"));
  DataCombData.plotOn(DataFrameF,Cut("passing==passing::fail"));
  w->pdf("DataModel")->plotOn(DataFrameF,Slice(*category,"fail"), ProjWData(*category,DataCombData), LineColor(kBlue),Range(60,120));
  w->pdf("DataModel")->plotOn(DataFrameF,Slice(*category,"fail"), ProjWData(*category,DataCombData), Components("DataSignalPdfF"), LineColor(kRed),Range(60,120));
  ///////////////////////////////////////////////////////////////

  c->Divide(2,3);

  c->cd(1);
  McFrameP->Draw();
  c->cd(2);
  McFrameF->Draw();
  c->cd(3);
  DataFrameP->Draw();
  c->cd(4);
  DataFrameF->Draw();
  c->cd(5);
  TemplateFrameP->Draw();
  c->cd(6);
  TemplateFrameP->Draw();
  //c->Update();
  c->Draw();

  test->cd();
  c->Write();
  test->Write();
  test->Close();

  delete c;

  if(verbose_) cout << "********** SGNtruePass " << SGNtruePass << endl; 
  if(verbose_) cout << "********** BKGPass " << BKGPass << endl; 
  if(verbose_) cout << "********** fullTreeSoup " << fullTreeSoup->GetEntries() << endl;
  if(verbose_) cout << "********** fullTreeData " << fullTreeData->GetEntries() << endl;

  //test->cd();
  //c->Write();

  //return;

}



void testTemplate( const string tnp_ = "etoTauMargLooseNoCracks70", 
		   const string category_ = "tauAntiEMVA",
		   const string bin_ = "abseta<1.5",
		   const float xLow_=60, 
		   const float xHigh_=120,
		   bool verbose_ = true){

  TFile *test = new TFile( "testTemplate.root" ,"RECREATE");

  TCanvas *c = new TCanvas("fitCanvas",Form("fitCanvas_%s_%s",tnp_.c_str(),bin_.c_str()),10,30,650,600);
  c->SetGrid(0,0);
  c->SetFillStyle(4000);
  c->SetFillColor(10);
  c->SetTicky();
  c->SetObjectStat(0);

  // input files
  TFile fsup("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup.root");
  TFile fbkg("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_bkg.root");
  TFile fsgn("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_sgn.root");
  TFile fdat("../tagAndProbe/trees/38XWcut/testNewWriteFromPAT_Data.root");

  //********************** signal only tree *************************/

  TTree *fullTreeSgn_ = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());
  test->cd();
  TTree *fullTreeSgn = (TTree*)fullTreeSgn_->CopyTree("","",100);

  TH1F* hSall = new TH1F("hSall","",1,0,150);
  TH1F* hSPall = new TH1F("hSPall","",1,0,150);
  TH1F* hS = new TH1F("hS","",1,0,150);
  TH1F* hSP = new TH1F("hSP","",1,0,150);
  fullTreeSgn->Draw("mass>>hS",Form("weight*(%s && mass>%f && mass<%f && mcTrue)",bin_.c_str(),xLow_,xHigh_));
  fullTreeSgn->Draw("mass>>hSall",Form("weight*(%s && mass>%f && mass<%f)",bin_.c_str(),xLow_,xHigh_));

  float SGNtrue = hS->Integral();
  float SGNall  = hSall->Integral();
 
  fullTreeSgn->Draw("mass>>hSP",Form("weight*(%s && %s>0 && mass>%f && mass<%f && mcTrue)",bin_.c_str(),category_.c_str(),xLow_,xHigh_));
  fullTreeSgn->Draw("mass>>hSPall",Form("weight*(%s && %s>0 && mass>%f && mass<%f)",bin_.c_str(),category_.c_str(),xLow_,xHigh_));

  float SGNtruePass = hSP->Integral();
  float SGNallPass  = hSPall->Integral();

  //********************** background only tree *************************//

  TTree *fullTreeBkg = (TTree*)fbkg.Get((tnp_+"/fitter_tree").c_str());
  TH1F* hB = new TH1F("hB","",1,0,150);
  TH1F* hBP = new TH1F("hBP","",1,0,150);
  fullTreeBkg->Draw("mass>>hB",Form("weight*(%s && mass>%f && mass<%f)",bin_.c_str(),xLow_,xHigh_));
 
  float BKG = hB->Integral();
  float BKGUnWeighted = hB->GetEntries();
  
  fullTreeBkg->Draw("mass>>hBP",Form("weight*(%s && %s>0 && mass>%f && mass<%f)",bin_.c_str(),category_.c_str(),xLow_,xHigh_));
  
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
  w->factory((category_+"[0,1]").c_str());

  RooRealVar  *weight = w->var("weight");
  RooRealVar  *abseta = w->var("abseta");
  RooRealVar  *pt = w->var("pt");
  RooRealVar  *mass = w->var("mass");
  mass->setRange(xLow_,xHigh_);
  RooRealVar  *mcTrue = w->var("mcTrue");
  RooRealVar  *cut = w->var( category_.c_str() );

  //mass->setBins(10000,"fft");

  w->factory("RooExponential::BackgroundPdf(mass,McCP[0,-10,10])");
  w->factory("SUM::SignalPdf(FVBP[0.5,0,1]*RooBifurGauss::Bif(mass,Mean[91.2,88,95],SigmaL[10,0.5,40],SigmaR[5.,0.5,15.]), RooVoigtian::Voig(mass, Mean, Width[2.49], SigmaVoig[5,0.1,10]) )");
  w->factory("RooCBShape::resPdf(mass,resMean[0.,0.,10],resSigma[0.5,0.0,8],resAlfa[1.8,0.,10],resN[1.0,1e-06,10])");
  //w->factory("RooFFTConvPdf::SignalPdf(mass,RooVoigtian::Voig(mass, Mean[91.2,88,95], Width[2.49], SigmaVoig[0.5,0,15]), resPdf)");
  w->factory("SUM::Model(NumSgn[1000,0,100000]*SignalPdf,NumBkg[100,0,100000]*BackgroundPdf)"); 

  w->saveSnapshot("clean", w->allVars());
  w->loadSnapshot("clean");
  
  RooDataSet sgnF("sgnF","dataset", RooArgSet(*mass,*weight,*abseta,*pt,*cut), Import( *fullTreeSgn ), Cut( Form("(%s<0.5 && %s)",category_.c_str(),bin_.c_str()) ) );
  RooDataSet sgnP("sgnP","dataset", RooArgSet(*mass,*weight,*abseta,*pt,*cut), Import( *fullTreeSgn ), Cut( Form("(%s>0.5 && %s)",category_.c_str(),bin_.c_str()) ) );
  
  RooFitResult* fitRes =  w->pdf("Model")->fitTo( sgnF ,Minos(1), Save(1), Range(xLow_,xHigh_), NumCPU(4));
  
  RooArgSet FitParam(fitRes->floatParsFinal());

  RooRealVar* NumSgn = (RooRealVar*)(&FitParam["NumSgn"]);
  RooRealVar* NumBkg = (RooRealVar*)(&FitParam["NumBkg"]);

  RooPlot* frame = mass->frame();
  dataF.plotOn(frame,Bins(60));
  w->pdf("Model")->plotOn(frame, LineColor(kBlue));
  w->pdf("Model")->plotOn(frame, LineColor(kRed), Components("SignalPdf"));
  w->pdf("Model")->plotOn(frame, LineColor(kGreen), LineStyle(kDashed), Components("BackgroundPdf"));
  frame->Draw();

  cout << "SGN from fit " << NumSgn->getVal() << " ; SGN from MC " << SGNtrue-SGNtruePass << endl;
  cout << "BKG from fit " << NumBkg->getVal() << " ; BKG from MC " << BKG-BKGPass+( (SGNall-SGNallPass)-(SGNtrue-SGNtruePass) )<< endl;

}

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooBreitWigner.h"
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



vector<Double_t*> simFit(TFile *outFile_        = 0,
			 const string tnp_      = "etoTauMargLooseNoCracks70", 
			 const string category_ = "tauAntiEMVA",
			 double cutValue_       = 0.5,
			 const string bin_      = "abseta<1.5",
			 const float binCenter_ = 0.75,
			 const float binWidth_  = 0.75,
			 const float xLow_      = 65, 
			 const float xHigh_     = 120,
			 const float nBins_     = 24,
			 bool doBinned_         = true,
			 float deltaAlpha_      = 0.0,
			 float deltaN_          = 0.0
			 ){

  vector<Double_t*> out;

  if(!outFile_){
    outFile_ = new TFile("outPutFile.root","RECREATE");
  }
  outFile_->mkdir(Form("bin%.2f",binCenter_));
  
  TCanvas *c2 = new TCanvas("fitCanvasTemplate","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);


  //TFile fsgn("/data_CMS/cms/lbianchini/tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_sgn.root");
  TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");
  TTree *fullTreeSgn  = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());
  //TFile fsup("/data_CMS/cms/lbianchini/tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup.root");
  TFile fsup("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA.root");
  TTree *fullTreeSoup = (TTree*)fsup.Get((tnp_+"/fitter_tree").c_str());
  //TFile fdat("/data_CMS/cms/lbianchini/tagAndProbe/trees/38XWcut/testNewWriteFromPAT_Data.root");
  TFile fdat("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Data.root");
  TTree *fullTreeData = (TTree*)fdat.Get((tnp_+"/fitter_tree").c_str());

  /////////////////////////////////////////////////

  TH1F* hS           = new TH1F("hS","",1,0,150);
  TH1F* hSP          = new TH1F("hSP","",1,0,150);

  fullTreeSgn->Draw("mass>>hS",Form("%s && mass>%f && mass<%f && mcTrue",bin_.c_str(),xLow_,xHigh_));
  float SGNtrue = hS->Integral();

  fullTreeSgn->Draw("mass>>hSP",Form("%s && %s>=%f && mass>%f && mass<%f && mcTrue",bin_.c_str(),category_.c_str(),cutValue_,xLow_,xHigh_));
  float SGNtruePass = hSP->Integral();

  float McTruthEff    = SGNtruePass/SGNtrue;
  float BinomialError = TMath::Sqrt(SGNtruePass/SGNtrue*(1-SGNtruePass/SGNtrue)/SGNtrue);
  
  cout << bin_.c_str() << " ==> MCTRUTH: " << McTruthEff << " +/- " << BinomialError << endl;

  delete hS; delete hSP;

  /////////////////////////////////////////////////

  // variables from the tree
  RooRealVar mcTrue("mcTrue","",0,1);
  RooRealVar matchedID("matchedID","",0,1);
  RooRealVar tauAntiEMVA("tauAntiEMVA","",0,1);
  RooRealVar signalPFChargedHadrCands("signalPFChargedHadrCands","",0,10);
  RooRealVar leadPFChargedHadrCandTrackPt("leadPFChargedHadrCandTrackPt","",0,200);
  RooRealVar pt("pt","",0,200);
  RooRealVar abseta("abseta","",0,10);

  // variables to constrct the pdfs
  ////////////////////////////////////////////////
  RooRealVar mass("mass","",xLow_,xHigh_);
  // ranges to define the fraction of events expected in the range:
 
  mass.setBins(10000,"fft");
  //mass.setBins( nBins_ );

  RooRealVar McNumBkgP("McNumBkgP","",0,1000000);
  RooRealVar McNumBkgF("McNumBkgF","",0,1000000);
  RooRealVar McNumSgn("McNumSgn","",0,1000000);
  RooRealVar McEfficiency("McEfficiency","",0.04,0,1);

  RooFormulaVar McNumSgnP("McNumSgnP","McEfficiency*McNumSgn",RooArgSet(McEfficiency,McNumSgn));
  RooFormulaVar McNumSgnF("McNumSgnF","(1-McEfficiency)*McNumSgn",RooArgSet(McEfficiency,McNumSgn));

  RooRealVar McCP("McCP","",0,-10,10);
  RooRealVar McCF("McCF","",0,-10,10);
  RooExponential McBackgroundPdfP("McBackgroundPdfP","",mass,McCP);
  RooExponential McBackgroundPdfF("McBackgroundPdfF","",mass,McCF);

  ////////////////////////////////////////////////

  RooRealVar DataNumBkgP("DataNumBkgP","",0,1000000);
  RooRealVar DataNumBkgF("DataNumBkgF","",0,1000000);
  RooRealVar DataNumSgn("DataNumSgn","",0,1000000);
  RooRealVar DataEfficiency("DataEfficiency","",0.04,0,1);

  RooFormulaVar DataNumSgnP("DataNumSgnP","DataEfficiency*DataNumSgn",RooArgSet(DataEfficiency,DataNumSgn));
  RooFormulaVar DataNumSgnF("DataNumSgnF","(1-DataEfficiency)*DataNumSgn",RooArgSet(DataEfficiency,DataNumSgn));
  RooRealVar DataCP("DataCP","",0,-10,10);
  RooRealVar DataCF("DataCF","",0,-10,10);

  RooExponential DataBackgroundPdfP("DataBackgroundPdfP","",mass,DataCP);
  RooExponential DataBackgroundPdfF("DataBackgroundPdfF","",mass,DataCF);

  ////////////////////////////////////////////////

  RooDataSet templateP("templateP","dataset for signal-pass template", RooArgSet(mass,mcTrue,tauAntiEMVA,signalPFChargedHadrCands,abseta,matchedID,leadPFChargedHadrCandTrackPt,pt), Import( *fullTreeSgn ), Cut(Form("(mcTrue && %s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) ) );

  //RooDataHist templateHistP("templateHistP","",RooArgSet(mass), templateP, 1.0);
  //RooHistPdf TemplateSignalPdfP("TemplateSignalPdfP","",RooArgSet(mass),templateHistP);

  // define the functional form for the passing signal in the template
  
  // breit-wigner
  RooRealVar mean("mean","mean",91.19);
  RooRealVar width("width","width",2.49);
  RooBreitWigner bw("bw","bw",mass,mean,width);

  // crystall ball
  RooRealVar m1("m1","m1",0,-20,20);
  RooRealVar sigma("sigma","sigma",0.5,0,20);
  RooRealVar alfa("alfa","alfa", 0.5,0,20);
  RooRealVar n("n","n", 1,1e-06,50);
  RooCBShape cb("cb","",mass,m1,sigma,alfa,n);

  // convolute
  RooFFTConvPdf TemplateSignalPdfP("TemplateSignalPdfP","",mass,bw,cb);

  // fit to the signal (unbinned)
  RooFitResult* ResTemplateFit = TemplateSignalPdfP.fitTo(templateP, Minos(1), Save(1), NumCPU(4) );
  RooArgSet FitParamTempl(ResTemplateFit->floatParsFinal());
  //RooRealVar* meanFit    = (RooRealVar*)(&FitParamTempl["mean"]);
  //RooRealVar* widthFit   = (RooRealVar*)(&FitParamTempl["width"]);
  RooRealVar* m1Fit    = (RooRealVar*)(&FitParamTempl["m1"]);
  RooRealVar* sigmaFit = (RooRealVar*)(&FitParamTempl["sigma"]);
  RooRealVar* alfaFit  = (RooRealVar*)(&FitParamTempl["alfa"]);
  RooRealVar* nFit     = (RooRealVar*)(&FitParamTempl["n"]);

  //RooRealVar mean_C("mean_C","mean",meanFit->getVal());
  //RooRealVar width_C("width_C","width",widthFit->getVal());
  RooBreitWigner bw_C("bw_C","bw",mass,mean/*_C*/,width/*_C*/);

  // constrained crystall ball
  RooRealVar m1_C("m1_C","m1",m1Fit->getVal());
  RooRealVar sigma_C("sigma_C","sigma",sigmaFit->getVal());
  RooRealVar alfa_C("alfa_C","alfa",alfaFit->getVal()*(1+deltaAlpha_));
  RooRealVar n_C("n_C","n", nFit->getVal()*(1+deltaN_));
  RooCBShape cb_tailC("cb_tailC","",mass,m1, sigma,   alfa_C,n_C);
  RooCBShape cb_C("cb_C","",        mass,m1_C,sigma_C,alfa_C,n_C);

  // convolute
  RooFFTConvPdf DataSignalPdfP("DataSignalPdfP","",mass,bw_C,cb_tailC);
  RooFFTConvPdf TemplateSignalPdfP_C("TemplateSignalPdfP_C","",mass,bw_C,cb_C);


  RooDataSet templateF("templateF","dataset for signal-pass template", RooArgSet(mass,mcTrue,tauAntiEMVA,signalPFChargedHadrCands,abseta,matchedID,leadPFChargedHadrCandTrackPt,pt), Import( *fullTreeSgn ), Cut(Form("(mcTrue && %s<%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) ));
  RooDataHist templateHistF("templateHistF","",RooArgSet(mass), templateF, 1.0);
  RooHistPdf TemplateSignalPdfF("TemplateSignalPdfF","",RooArgSet(mass),templateHistF);

  RooPlot* TemplateFrameP = mass.frame(Bins(nBins_),Title("Template passing"));
  templateP.plotOn(TemplateFrameP);
  TemplateSignalPdfP.plotOn(TemplateFrameP);
  
  RooPlot* TemplateFrameF = mass.frame(Bins(nBins_),Title("Template failing"));
  templateF.plotOn(TemplateFrameF);
  TemplateSignalPdfF.plotOn(TemplateFrameF);

  c2->Divide(2,1);
  c2->cd(1);
  TemplateFrameP->Draw();
  c2->cd(2);
  TemplateFrameF->Draw();
  c2->Draw();

  outFile_->cd(Form("bin%.2f",binCenter_));
  c2->Write();
  ////////////////////////////////////////////////

  RooCategory category("category","category") ;
  category.defineType("pass") ;
  category.defineType("fail") ;

  ////////////////////////////////////////////////
  //                 Soup
  ////////////////////////////////////////////////

  mass.setBins( nBins_ );

  RooAddPdf McModelP("McModelP","",RooArgList(McBackgroundPdfP,TemplateSignalPdfP_C),RooArgList(McNumBkgP,McNumSgnP));
  RooAddPdf McModelF("McModelF","",RooArgList(McBackgroundPdfF,TemplateSignalPdfF),RooArgList(McNumBkgF,McNumSgnF));

  TFile *McF = new TFile("dummySoup.root","RECREATE");
  TTree* treeSoupP = fullTreeSoup->CopyTree( Form("(%s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) );
  TTree* treeSoupF = fullTreeSoup->CopyTree( Form("(%s<%f && %s)", category_.c_str(),cutValue_,bin_.c_str()) );

  RooDataSet McDataP("McDataP","dataset pass for the soup pass", RooArgSet(mass), Import( *treeSoupP ) );
  RooDataSet McDataF("McDataF","dataset pass for the soup fail", RooArgSet(mass), Import( *treeSoupF ) );

  RooDataHist McCombData("McCombData","combined data",mass,Index(category),Import("pass", *(McDataP.createHistogram("histoSoupP",mass))) ,Import("fail", *(McDataF.createHistogram("histoSoupF",mass))) ) ;

  RooSimultaneous McSimPdf("McSimPdf","simultaneous pdf",category) ;
  McSimPdf.addPdf(McModelP,"pass") ;
  McSimPdf.addPdf(McModelF,"fail") ;

  RooFitResult* ResMcCombinedFit = McSimPdf.fitTo(McCombData , Extended(1), Minos(1), Save(1), NumCPU(4) );
  outFile_->cd(Form("bin%.2f",binCenter_));
  ResMcCombinedFit->Write("McFitResults_Combined");

  
  RooArgSet McFitParam(ResMcCombinedFit->floatParsFinal());
  RooRealVar* McEffFit     = (RooRealVar*)(&McFitParam["McEfficiency"]);
  RooRealVar* McNumSigFit  = (RooRealVar*)(&McFitParam["McNumSgn"]);
  RooRealVar* McNumBkgPFit = (RooRealVar*)(&McFitParam["McNumBkgP"]);
  RooRealVar* McNumBkgFFit = (RooRealVar*)(&McFitParam["McNumBkgF"]);

  RooPlot* McFrameP = mass.frame(Bins(nBins_),Title("MC: passing sample"));
  McCombData.plotOn(McFrameP,Cut("category==category::pass"));
  McSimPdf.plotOn(McFrameP,Slice(category,"pass"), ProjWData(category,McCombData), LineColor(kBlue));
  McSimPdf.plotOn(McFrameP,Slice(category,"pass"), ProjWData(category,McCombData), Components("TemplateSignalPdfP_C"), LineColor(kRed));
  McSimPdf.plotOn(McFrameP,Slice(category,"pass"), ProjWData(category,McCombData), Components("McBackgroundPdfP"), LineColor(kGreen), LineStyle(kDashed));

  RooPlot* McFrameF = mass.frame(Bins(nBins_),Title("MC: failing sample"));
  McCombData.plotOn(McFrameF,Cut("category==category::fail"));
  McSimPdf.plotOn(McFrameF,Slice(category,"fail"), ProjWData(category,McCombData), LineColor(kBlue));
  McSimPdf.plotOn(McFrameF,Slice(category,"fail"), ProjWData(category,McCombData), Components("TemplateSignalPdfF"), LineColor(kRed));
  McSimPdf.plotOn(McFrameF,Slice(category,"fail"), ProjWData(category,McCombData), Components("McBackgroundPdfF"), LineColor(kGreen), LineStyle(kDashed));

  ////////////////////////////////////////////////
  //                 Data
  ////////////////////////////////////////////////

  mass.setBins( nBins_ );

  RooRealVar DataMeanResP("DataMeanResP","",0,-10,10);
  RooRealVar DataSigmaResP("DataSigmaResP","",0.5,0,10);
  RooRealVar DataMeanResF("DataMeanResF","",0,-10,10);
  RooRealVar DataSigmaResF("DataSigmaResF","",0.5,0,10);

  // convolve the template with a resolution model
  //RooGaussian DataResolModP("DataResolModP","",mass,DataMeanResP,DataSigmaResP);
  //RooFFTConvPdf DataSignalPdfP("DataSignalPdfP","",mass,TemplateSignalPdfP,DataResolModP);
  // convolve the template with a resolution model 
  RooGaussian DataResolModF("DataResolModF","",mass,DataMeanResF,DataSigmaResF);
  RooFFTConvPdf DataSignalPdfF("DataSignalPdfF","",mass,TemplateSignalPdfF,DataResolModF);

  RooAddPdf DataModelP("DataModelP","",RooArgList(DataBackgroundPdfP,DataSignalPdfP),RooArgList(DataNumBkgP,DataNumSgnP));
  RooAddPdf DataModelF("DataModelF","",RooArgList(DataBackgroundPdfF,DataSignalPdfF),RooArgList(DataNumBkgF,DataNumSgnF));

  TFile *DataF = new TFile("dummyData.root","RECREATE");
  TTree* treeDataP = fullTreeData->CopyTree( Form("(%s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()));
  TTree* treeDataF = fullTreeData->CopyTree( Form("(%s<%f && %s)",category_.c_str(),cutValue_,bin_.c_str()));

  RooDataSet DataDataP("DataDataP","dataset pass for the data pass", RooArgSet(mass), Import( *treeDataP ) );
  RooDataSet DataDataF("McDataF","dataset pass for the data fail", RooArgSet(mass), Import( *treeDataF ) );

  // binned combined dataset
  RooDataHist DataCombData("DataCombData","combined data",mass,Index(category),Import("pass", *(DataDataP.createHistogram("histoDataP",mass)) ) ,Import("fail", *(DataDataF.createHistogram("histoDataF",mass))) ) ;
  // unbinned combined dataset
  RooDataSet DataCombDataUnBinned("DataCombDataUnBinned","combined data",mass,Index(category),Import("pass", DataDataP ) ,Import("fail",DataDataF) ) ;

  RooSimultaneous DataSimPdf("DataSimPdf","simultaneous pdf",category) ;
  DataSimPdf.addPdf(DataModelP,"pass") ;
  DataSimPdf.addPdf(DataModelF,"fail") ;

  RooFitResult* ResDataCombinedFit =  0;
  if(doBinned_)  ResDataCombinedFit = DataSimPdf.fitTo(DataCombData , Extended(1), Minos(1), Save(1), NumCPU(4) );
  else ResDataCombinedFit = DataSimPdf.fitTo(DataCombDataUnBinned , Extended(1), Minos(1), Save(1), NumCPU(4) );
  outFile_->cd(Form("bin%.2f",binCenter_));
  ResDataCombinedFit->Write("DataFitResults_Combined");

  
  RooArgSet DataFitParam(ResDataCombinedFit->floatParsFinal());
  RooRealVar* DataEffFit      = (RooRealVar*)(&DataFitParam["DataEfficiency"]);
  RooRealVar* DataNumSigFit   = (RooRealVar*)(&DataFitParam["DataNumSgn"]);
  RooRealVar* DataNumBkgPFit  = (RooRealVar*)(&DataFitParam["DataNumBkgP"]);
  RooRealVar* DataNumBkgFFit  = (RooRealVar*)(&DataFitParam["DataNumBkgF"]);
  RooRealVar* DataMeanResFFit = (RooRealVar*)(&DataFitParam["DataMeanResF"]);
  RooRealVar* DataSigmaResFFit= (RooRealVar*)(&DataFitParam["DataSigmaResF"]); 


  RooPlot* DataFrameP = mass.frame(Bins(nBins_),Title("Data: passing sample"));
  DataCombData.plotOn(DataFrameP,Cut("category==category::pass"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), LineColor(kBlue));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("DataSignalPdfP"), LineColor(kRed));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("DataBackgroundPdfP"), LineColor(kGreen), LineStyle(kDashed));

  RooPlot* DataFrameF = mass.frame(Bins(nBins_),Title("Data: failing sample"));
  DataCombData.plotOn(DataFrameF,Cut("category==category::fail"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), LineColor(kBlue));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("DataSignalPdfF"), LineColor(kRed));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("DataBackgroundPdfF"), LineColor(kGreen), LineStyle(kDashed));

  TCanvas *c = new TCanvas("fitCanvas","canvas",10,30,650,600);
  c->SetGrid(0,0);
  c->SetFillStyle(4000);
  c->SetFillColor(10);
  c->SetTicky();
  c->SetObjectStat(0);

  c->Divide(2,2);
  c->cd(1);
  DataFrameP->Draw();
  c->cd(2);
  DataFrameF->Draw();
  c->cd(3);
  McFrameP->Draw();
  c->cd(4);
  McFrameF->Draw();
  c->Draw();

  outFile_->cd(Form("bin%.2f",binCenter_));
  c->Write();

  //outFile_->Close();


  // MINOS errors, otherwise HESSE quadratic errors
  float McErrorLo = McEffFit->getErrorLo()<0 ? McEffFit->getErrorLo() : (-1)*McEffFit->getError();
  float McErrorHi = McEffFit->getErrorHi()>0 ? McEffFit->getErrorHi() : McEffFit->getError();
  
  float DataErrorLo = DataEffFit->getErrorLo()<0 ? DataEffFit->getErrorLo() : (-1)*DataEffFit->getError();
  float DataErrorHi = DataEffFit->getErrorHi()>0 ? DataEffFit->getErrorHi() : DataEffFit->getError();
 
  Double_t* truthMC = new Double_t[6];
  Double_t* tnpMC   = new Double_t[6];
  Double_t* tnpData = new Double_t[6];

  truthMC[0] = binCenter_;
  truthMC[1] = binWidth_;
  truthMC[2] = binWidth_;
  truthMC[3] = McTruthEff;
  truthMC[4] = BinomialError;
  truthMC[5] = BinomialError;

  tnpMC[0] = binCenter_;
  tnpMC[1] = binWidth_;
  tnpMC[2] = binWidth_;
  tnpMC[3] = McEffFit->getVal();
  tnpMC[4] = (-1)*McErrorLo;
  tnpMC[5] = McErrorHi;
 
  tnpData[0] = binCenter_;
  tnpData[1] = binWidth_;
  tnpData[2] = binWidth_;
  tnpData[3] = DataEffFit->getVal();
  tnpData[4] = (-1)*DataErrorLo;
  tnpData[5] = DataErrorHi;

  out.push_back(truthMC);
  out.push_back(tnpData);
  out.push_back(tnpMC);

  return out;

 
}





void makePlot(const string tnp_      = "etoTauMargLooseNoCracks70",
	      const string category_ = "tauAntiEMVA",
	      double cutValue_       = 0.5,
	      const string var_      = "abseta",
	      const string otherCuts_= "abseta>-1",
	      const float xLow_      = 65, 
	      const float xHigh_     = 120,
	      const float nBins_     = 18,
	      bool doBinned_         = true,
	      double bin1_           = 0.0,
	      double bin2_           = 1.5,
	      double bin3_           = 2.5,
	      float deltaAlpha_      = 0.0,
	      float deltaN_          = 0.0
	      ){
  
    // output file
  TFile *outFile = new TFile( Form("EtoTauPlots_%s_%s_%s.root",tnp_.c_str(),category_.c_str(),var_.c_str()),"RECREATE");

  cout << "******************** bin1" << endl;
  vector<Double_t*> bin1Results = simFit(outFile,tnp_, category_,cutValue_,
					 string(Form("%s && %s>%f && %s<%f",otherCuts_.c_str(),var_.c_str(),bin1_,var_.c_str(),bin2_)),
					 (bin2_+bin1_)/2,(bin2_-bin1_)/2 
					 ,xLow_, xHigh_, nBins_,doBinned_,
					 deltaAlpha_,deltaN_);

  cout << "******************** bin2" << endl;
  vector<Double_t*> bin2Results = simFit(outFile,tnp_, category_, cutValue_,
					 string(Form("(%s && %s > %f && %s<%f)",otherCuts_.c_str(),var_.c_str(),bin2_,var_.c_str(),bin3_)),
					 (bin3_+bin2_)/2,(bin3_-bin2_)/2 ,
					 xLow_, xHigh_, nBins_,doBinned_,
					 deltaAlpha_,deltaN_);
  
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
  
  tnpMC_x[0]  = (bin1Results[2])[0]; tnpMC_x[1]  = (bin2Results[2])[0];
  tnpMC_xL[0] = (bin1Results[2])[1]; tnpMC_xL[1] = (bin2Results[2])[1];
  tnpMC_xH[0] = (bin1Results[2])[2]; tnpMC_xH[1] = (bin2Results[2])[2];
  tnpMC_y[0]  = (bin1Results[2])[3]; tnpMC_y[1]  = (bin2Results[2])[3];
  tnpMC_yL[0] = (bin1Results[2])[4]; tnpMC_yL[1] = (bin2Results[2])[4];
  tnpMC_yH[0] = (bin1Results[2])[5]; tnpMC_yH[1] = (bin2Results[2])[5];
  
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
  TH1F* h1 = new TH1F("h1","",2, binsEdges);
  h1->SetAxisRange(0,0.06,"Y");
  h1->SetXTitle( "|#eta|" );
  string YTitle="fake-rate";
  h1->SetYTitle( YTitle.c_str() );
  
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
  graph_tnpMC->Draw("PSAME");
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

  string discr = "" ;
  if(category_.find("tauAntiEMVA")!=string::npos) discr="discr. by mva";
  else discr="discr. by WP95";
  string tau = "" ;
  if(tnp_.find("SC")==string::npos) tau="HPS #tau-candidate";
  else tau="Shrinking-Cone #tau-candidate";
  leg->SetHeader(Form("#splitline{CMS Preliminary L=35 pb^{-1}}{%s passing %s}",tau.c_str(),discr.c_str()));
  leg->AddEntry(h1mcTruth,"mc-truth");
  leg->AddEntry(h1tnpMC,"t&p: simulation");
  leg->AddEntry(h1tnpDATA,"t&p: 7 TeV Data");
  leg->Draw();

  outFile->cd();
  c1->Write();
  outFile->Close();

}



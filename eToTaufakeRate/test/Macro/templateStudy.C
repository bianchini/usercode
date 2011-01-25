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


void templateStudy(){

  TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");
  TFile fsoup("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA.root");
  //TFile fsoup("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");

  TTree *fullTreeSgn   = (TTree*)fsgn.Get("etoTauMargLooseNoCracks70/fitter_tree");
  TTree *fullTreeSoup  = (TTree*)fsoup.Get("etoTauMargTightNoCracks60/fitter_tree");

  RooRealVar mcTrue("mcTrue","",0,1);
  RooRealVar matchedID("matchedID","",0,1);
  RooRealVar tauAntiEMVA("tauAntiEMVA","",0,1);
  RooRealVar signalPFChargedHadrCands("signalPFChargedHadrCands","",0,10);
  RooRealVar electronPreIDOutput("electronPreIDOutput","",-1000,1.0);
  RooRealVar abseta("abseta","",0,10);
  RooRealVar mass("mass","mass",0,140);

  mass.setRange("integralAll",0,140);
  mass.setRange("integral",65,120);

  mass.setBins(24);

  TFile * f = new TFile("dimmy.root","RECREATE");
  TTree *fullTreeSoupCut      = fullTreeSoup->CopyTree("matchedID<=0.95 && electronPreIDOutput<0.1 && abseta<1.5");
  TTree *fullTreeSoupCutRight = fullTreeSoup->CopyTree("tauAntiEMVA>0.5 && abseta<1.5");
  TTree *fullTreeSgnCut       = fullTreeSgn->CopyTree("tauAntiEMVA>0.5 && abseta<1.5 && mcTrue");

  RooDataSet templateDataSet("templateDataSet","", RooArgSet(mass),   Import( *fullTreeSoupCut ));
  RooDataSet passingDataSet("passingDataSet","",   RooArgSet(mass),   Import( *fullTreeSgnCut ) );
  RooDataSet passingDataSetR("passingDataSetR","", RooArgSet(mass),   Import( *fullTreeSoupCutRight ) );


  RooDataHist passingHist("passingHist","",RooArgSet(mass), passingDataSet , 1.0);
  RooHistPdf passingPdf("passingPdf","",RooArgSet(mass),passingHist);

  //RooDataHist passingHistR("passingHistR","",RooArgSet(mass), passingDataSetR , 1.0);
  //RooHistPdf passingPdfR("passingPdfR","",RooArgSet(mass),passingHistR);
 
  RooKeysPdf passingPdfR("passingPdfR","",mass,passingDataSetR);
  RooKeysPdf templateP("templateP","",mass,templateDataSet);

  RooPlot* frame = mass.frame( Bins(24) );
  //templateDataSet.plotOn(frame);
  templateP.plotOn(frame,LineColor(kBlue),LineStyle(kDashed));
  passingPdf.plotOn(frame,LineColor(kRed));

  RooRealVar height1("height1","",1);
  RooRealVar height0("height0","",0);
  TArrayD limits(3);
  limits[0]= 60;
  limits[1]= 85;
  limits[2]= 120;
 
  RooArgSet set(mass);

  RooParametricStepFunction stepLeft("stepLeft","",mass,RooArgList(height1,height0),limits,2);
  RooParametricStepFunction stepRight("stepRight","",mass,RooArgList(height0,height1),limits,2);

  RooProdPdf pdfLeft("pdfLeft","",RooArgList(stepLeft,templateP));
  RooProdPdf pdfRight("pdfRight","",RooArgList(stepRight,passingPdfR));

  mass.setVal(limits[1]-0.01);
  double patchLeft  =  pdfLeft.getVal(&set); 
  mass.setVal(limits[1]+0.01);
  double patchRight =  pdfRight.getVal(&set);

  cout << patchLeft << "   " << patchRight << endl;

  RooRealVar NLeft("NLeft","",1.0);
  RooRealVar NRight("NRight","",patchLeft/patchRight);
  
  RooAddPdf sum("sum","",RooArgList(pdfLeft,pdfRight),RooArgList(NLeft,NRight));

  sum.plotOn(frame,LineColor(kGreen));

  frame->Draw();

  TLegend* leg = new TLegend(0.6,0.6,0.89,0.89,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);

  TH1F* tail = new TH1F("tail","",1,60,120);
  TH1F* core = new TH1F("core","",1,60,120);
  TH1F* all  = new TH1F("all","",1,60,120);
  tail->SetLineColor(kBlue);
  core->SetLineColor(kRed);
  all->SetLineColor(kGreen);
  tail->SetLineStyle(kDashed);
  core->SetLineStyle(kSolid);
  all->SetLineStyle(kSolid);
  tail->SetLineWidth(0.4);
  core->SetLineWidth(0.4);
  all->SetLineWidth(0.4);

  leg->SetHeader("HPS, |#eta|<1.5");
  leg->AddEntry(core,"MC: passing probes (MC-truth)","l");
  leg->AddEntry(tail, "MC: WP95 && mva<0.1 from soup","l");
  leg->AddEntry(all,"MC: template from soup","l");

  leg->Draw();

  RooRealVar* integral = (RooRealVar*)passingPdf.createIntegral(mass,Range("integral"));
  RooRealVar* integralAll = (RooRealVar*)passingPdf.createIntegral(mass,Range("integralAll"));
  cout << integral->getVal() << " and all " <<   integralAll->getVal()<< endl;

}

void templateStudyOnData(){

  TFile fdata("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Data.root");

  TTree *fullTreeData   = (TTree*)fdata.Get("etoTauMargLooseNoCracks60/fitter_tree");

  RooRealVar matchedID("matchedID","",0,1);
  RooRealVar tauAntiEMVA("tauAntiEMVA","",0,1);
  RooRealVar electronPreIDOutput("electronPreIDOutput","",-1000,1.0);
  RooRealVar abseta("abseta","",0,10);
  RooRealVar mass("mass","",60,120);

  mass.setBins(18);

  TFile * f = new TFile("dimmy.root","RECREATE");
  TTree *fullTreeDataCutLeft  = fullTreeData->CopyTree("matchedID<=0.95 && electronPreIDOutput<0.1 && abseta<1.5");
  TTree *fullTreeDataCutRight = fullTreeData->CopyTree("tauAntiEMVA>0.5 && abseta<1.5");

  RooDataSet passingDataSetL("passingDataSetL","", RooArgSet(mass),   Import( *fullTreeDataCutLeft )  );
  RooDataSet passingDataSetR("passingDataSetR","", RooArgSet(mass),   Import( *fullTreeDataCutRight ) );

  RooDataHist passingHistL("passingHistL","",RooArgSet(mass), passingDataSetL , 1.0);
  RooHistPdf passingPdfL("passingPdfL","",RooArgSet(mass),passingHistL);
  RooDataHist passingHistR("passingHistR","",RooArgSet(mass), passingDataSetR , 1.0);
  RooHistPdf passingPdfR("passingPdfR","",RooArgSet(mass),passingHistR);
 
  //RooKeysPdf passingPdfR("passingPdfR","",mass,passingDataSetR);
  //RooKeysPdf passingPdfL("passingPdfL","",mass,passingDataSetL);

  RooPlot* frame = mass.frame( Bins(12) );
  passingDataSetR.plotOn(frame);

  RooRealVar height1("height1","",1);
  RooRealVar height0("height0","",0);
  TArrayD limits(3);
  limits[0]= 60;
  limits[1]= 85;
  limits[2]= 120;
 
  RooArgSet set(mass);

  RooParametricStepFunction stepLeft("stepLeft","",mass,RooArgList(height1,height0),limits,2);
  RooParametricStepFunction stepRight("stepRight","",mass,RooArgList(height0,height1),limits,2);

  RooProdPdf pdfLeft("pdfLeft","",RooArgList(stepLeft,passingPdfL));
  RooProdPdf pdfRight("pdfRight","",RooArgList(stepRight,passingPdfR));

  mass.setVal(84.99);
  double patchLeft  =  pdfLeft.getVal(&set); 
  mass.setVal(85.01);
  double patchRight =  pdfRight.getVal(&set);

  cout << patchLeft << "   " << patchRight << endl;

  RooRealVar NLeft("NLeft","",1.0);
  RooRealVar NRight("NRight","",patchLeft/patchRight);
  
  RooAddPdf sum("sum","",RooArgList(pdfLeft,pdfRight),RooArgList(NLeft,NRight));
  passingPdfR.plotOn(frame,LineColor(kRed));
  passingPdfL.plotOn(frame,LineColor(kBlue));
  sum.plotOn(frame,LineColor(kGreen));

  frame->Draw();
}


void templateStudyOnData2(){

  TFile fdata("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Data.root");

  TTree *fullTreeData   = (TTree*)fdata.Get("etoTauSCMargNoCracks70/fitter_tree");

  RooRealVar tauAntiEMVA("tauAntiEMVA","",0,1);
  RooRealVar electronPreIDOutput("electronPreIDOutput","",-1000,1.0);
  RooRealVar abseta("abseta","",0,10);
  RooRealVar mass("mass","",60,120);
  RooRealVar pt("pt","",15,100);

  mass.setBins(20);

  TFile * f = new TFile("dimmy.root","RECREATE");
  TTree *fullTreeDataCutLeft  = fullTreeData->CopyTree("pt>35 && electronPreIDOutput<0.1 && abseta>1.5");

  RooDataSet passingDataSetL("passingDataSetL","", RooArgSet(mass),   Import( *fullTreeDataCutLeft )  );

 
  RooKeysPdf passingPdfL("passingPdfL","",mass,passingDataSetL);

  RooPlot* frame = mass.frame( Bins(20) );
  passingDataSetL.plotOn(frame);


  passingPdfL.plotOn(frame,LineColor(kBlue));
 
  frame->Draw();
}


/*

TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");
TFile fsoup("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA.root");

  TTree *fullTreeSgn   = (TTree*)fsgn.Get("etoTauMargTightNoCracks60/fitter_tree");
  TTree *fullTreeSoup  = (TTree*)fsoup.Get("etoTauMargTightNoCracks60/fitter_tree");

  TH1F* hVBTF = new TH1F("hVBTF","",15,60,120);
  TH1F* hPF   = new TH1F("hPF","",15,60,120);
 
  hVBTF->SetLineColor(kBlue);
  hPF->SetLineColor(kRed);

  fullTreeSoup->Draw("mass>>hVBTF","(matchedID<=0.95) && electronPreIDOutput<0.1 && abseta<1.5");
  fullTreeSgn->Draw("mass>>hPF","mcTrue && tauAntiEMVA>0.5 && abseta<1.5");

  hVBTF->DrawNormalized("HIST");
  hPF->DrawNormalized("HISTSAME");

*/




void extractFromData( const string tnp_      = "etoTauSCMargNoCracks60",
		      bool doBinned_         = true,
		      const string category_ = "tauAntiEMVA",
		      double cutValue_       = 0.5,
		      const string bin_      = "abseta>1.5",
		      double nBins_          = 18,
		      double xLow_           = 65,
		      double xHigh_          = 120
		      ){

  TCanvas *c2 = new TCanvas("fitCanvasTemplate","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  // soup  
  //TFile fsup("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA.root");
  TFile fsup("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Data.root");
  TTree *fullTreeSoup = (TTree*)fsup.Get((tnp_+"/fitter_tree").c_str());

  // signal
  TFile fsgnHighStat("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");
  TTree *fullTreeSgnHighStat = (TTree*)fsgnHighStat.Get((tnp_+"/fitter_tree").c_str());

  // variables
  RooRealVar mcTrue("mcTrue","",0,1);
  RooRealVar matchedID("matchedID","",0,1);
  RooRealVar electronPreIDOutput("electronPreIDOutput","",-1000,1);
  RooRealVar tauAntiEMVA("tauAntiEMVA","",0,1);
  RooRealVar signalPFChargedHadrCands("signalPFChargedHadrCands","",0,10);
  RooRealVar leadPFChargedHadrCandTrackPt("leadPFChargedHadrCandTrackPt","",0,200);
  RooRealVar leadPFCandPt("leadPFCandPt","",0,200);
  RooRealVar pt("pt","",0,200);
  RooRealVar abseta("abseta","",0,10);
  RooRealVar mass("mass","",xLow_,xHigh_);

  mass.setBins(10000,"fft");
  mass.setBins( nBins_ );

  ///////////////// signal pdf from MC ////////////////////
  RooDataSet templateP("templateP","dataset for signal-pass template", RooArgSet(mass,mcTrue,tauAntiEMVA,signalPFChargedHadrCands,abseta,matchedID,leadPFChargedHadrCandTrackPt,pt), Import( *fullTreeSgnHighStat ), Cut(Form("(mcTrue && %s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) ) );
  RooDataHist templateHistP("templateHistP","",RooArgSet(mass), templateP, 1.0);
  RooHistPdf TemplateSignalPdfP("TemplateSignalPdfP","",RooArgSet(mass),templateHistP);

  RooRealVar McMeanResP("McMeanResP","",0,-10,10);
  RooRealVar McSigmaResP("McSigmaResP","",0.5,0,10);
  RooGaussian McResolModP("McResolModP","",mass,McMeanResP,McSigmaResP);
  RooFFTConvPdf McSignalPdfP("McSignalPdfP","",mass,TemplateSignalPdfP,McResolModP);

  /////////////////      bkg pdf       ////////////////////
  RooRealVar McCP("McCP","",0,-10,10);
  RooExponential McBackgroundPdfP("McBackgroundPdfP","",mass,McCP);


  RooRealVar McNumBkgP("McNumBkgP","",0,1000000);
  RooRealVar McNumSgnP("McNumSgnP","",0,1000000);
  
  // pdf for combined sample
  RooAddPdf McModelPResol("McModelPResol","",RooArgList(McBackgroundPdfP,McSignalPdfP),RooArgList(McNumBkgP,McNumSgnP));

  // passing dataset for the soup
  TFile *McP = new TFile("dummy1.root","RECREATE");
  TTree* treeSoupP = fullTreeSoup->CopyTree( Form("(%s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) );
  RooDataSet  McDataPSet("McDataPSet","dataset pass for the soup pass", RooArgSet(mass), Import( *treeSoupP ) );
  RooDataHist McDataP("McDataP","dataset pass for the soup pass", RooArgSet(mass), McDataPSet, 1.0 );

  // fit the passing suop to the MC(X)resol model
  RooFitResult* ResMcCombinedResolFit = 0;
  if(doBinned_) ResMcCombinedResolFit = McModelPResol.fitTo(McDataP , Extended(1), Minos(1), Save(1), NumCPU(4) );
  else  ResMcCombinedResolFit = McModelPResol.fitTo(McDataPSet , Extended(1), Minos(1), Save(1), NumCPU(4) );
  RooArgSet McFitParamResol(ResMcCombinedResolFit->floatParsFinal());
  RooRealVar* McNumSigPFitResol  = (RooRealVar*)(&McFitParamResol["McNumSgnP"]);
  RooRealVar* McNumBkgPFitResol  = (RooRealVar*)(&McFitParamResol["McNumBkgP"]);

  // now, let's try the template
  RooDataSet templatePfromData("templatePfromData","dataset for signal-pass template", RooArgSet(mass,tauAntiEMVA,electronPreIDOutput,signalPFChargedHadrCands,abseta,matchedID,leadPFChargedHadrCandTrackPt,leadPFCandPt,pt), Import( *fullTreeSoup ), Cut(Form("(electronPreIDOutput < 0.0 && leadPFChargedHadrCandTrackPt>15  && leadPFCandPt>15 && signalPFChargedHadrCands<1.5 && %s)",bin_.c_str()) ) );
  
  // define an expo to account for residual bkg
  RooRealVar McCPfromData("McCPfromData","",0,-10,10);
  RooExponential McBackgroundPdfPfromData("McBackgroundPdfPfromData","",mass,McCPfromData);
  RooRealVar McNumBkgPfromData("McNumBkgPfromData","",0,10000);
  RooRealVar McNumSgnPfromData("McNumSgnPfromData","",0,10000);

  // define the functional form for the passing signal in the template
  // breit-wigner
  RooRealVar mean("mean","mean",91,85,95);
  RooRealVar width("widthc","width",2.49);
  RooBreitWigner bw("bw","bw",mass,mean,width);

  // constrained crystall ball
  RooRealVar m1("m1","m1",0,-10,10);
  RooRealVar sigma("sigma","sigma",0.5,0,10);
  RooRealVar alfa("alfa","alfa", 0.5,0,10);
  RooRealVar n("n","n", 1,1e-06,10);
  RooCBShape cb("cb","",mass,m1,sigma,alfa,n);

  // convolute
  RooFFTConvPdf McSignalPdfPfromData("McSignalPdfPfromData","",mass,bw,cb);
  RooAddPdf McModelPfromData("McModelPResolfromData","",RooArgList(McBackgroundPdfPfromData, McSignalPdfPfromData ),RooArgList(McNumBkgPfromData,McNumSgnPfromData));
  // fit the signal in the data template
  RooFitResult* ResMcCombinedTemplateFit = McModelPfromData.fitTo(templatePfromData , Extended(1), Minos(1), Save(1), NumCPU(4) );

  // define the constrained signal:
  RooArgSet McFitParamTempl(ResMcCombinedTemplateFit->floatParsFinal());
  RooRealVar* meanFit  = (RooRealVar*)(&McFitParamTempl["mean"]);
  RooRealVar* m1Fit    = (RooRealVar*)(&McFitParamTempl["m1"]);
  RooRealVar* sigmaFit = (RooRealVar*)(&McFitParamTempl["sigma"]);
  RooRealVar* alfaFit  = (RooRealVar*)(&McFitParamTempl["alfa"]);
  RooRealVar* nFit     = (RooRealVar*)(&McFitParamTempl["n"]);

  RooRealVar mean_C("mean_C","mean",meanFit->getVal());
  RooBreitWigner bw_C("bw_C","bw",mass,mean_C,width);

  // constrained crystall ball
  RooRealVar m1_C("m1_C","m1",m1Fit->getVal());
  RooRealVar sigma_C("sigma_C","sigma",sigmaFit->getVal());
  RooRealVar alfa_C("alfa_C","alfa",alfaFit->getVal());
  RooRealVar n_C("n_C","n", nFit->getVal());
  RooCBShape cb_C("cb_C","",mass,m1_C,sigma_C,alfa_C,n_C);

  // convolute
  RooFFTConvPdf McSignalPdfPfromData_C("McSignalPdfPfromData_C","",mass,bw_C,cb_C);

  // pdf for combined sample
  RooAddPdf McModelPfromData_C("McModelPfromData_C","",RooArgList(McBackgroundPdfPfromData,McSignalPdfPfromData_C),RooArgList(McNumBkgP,McNumSgnP));

  // ... and fit
  RooFitResult* ResMcCombinedFromDataFit = 0;
  if(doBinned_) ResMcCombinedFromDataFit = McModelPfromData_C.fitTo(McDataP , Extended(1), Minos(1), Save(1), NumCPU(4) );
  else   ResMcCombinedFromDataFit = McModelPfromData_C.fitTo(McDataPSet , Extended(1), Minos(1), Save(1), NumCPU(4) );

  RooArgSet McFitParamFromData(ResMcCombinedFromDataFit->floatParsFinal());
  RooRealVar* McNumSigPFitFromData  = (RooRealVar*)(&McFitParamFromData["McNumSgnP"]);
  RooRealVar* McNumBkgPFitFromData = (RooRealVar*)(&McFitParamFromData["McNumBkgP"]);
  

  c2->Divide(1,2);

  c2->cd(1);
  RooPlot* McFrameP = mass.frame(Bins(nBins_),Title("MC: passing sample"));
  McDataP.plotOn(McFrameP);
  McModelPResol.plotOn(McFrameP,LineColor(kRed)); 
  McModelPResol.plotOn(McFrameP,LineColor(kRed),Components("McSignalPdfP"),LineStyle(kDashed)); 
  McModelPResol.plotOn(McFrameP,LineColor(kRed),Components("McBackgroundPdfP"),LineStyle(kDotted));
 
  McModelPfromData_C.plotOn(McFrameP,LineColor(kBlue)); 
  McModelPfromData_C.plotOn(McFrameP,LineColor(kBlue),Components("McSignalPdfPfromData_C"),LineStyle(kDashed)); 
  McModelPfromData_C.plotOn(McFrameP,LineColor(kBlue),Components("McBackgroundPdfPfromData"),LineStyle(kDotted)); 

  McFrameP->Draw();

  c2->cd(2);
  RooPlot* McFramePfromData = mass.frame(Bins(nBins_),Title("MC: passing sample from data"));
  templatePfromData.plotOn(McFramePfromData);
  McModelPfromData.plotOn(McFramePfromData,LineColor(kBlue));
  McModelPfromData.plotOn(McFramePfromData,LineColor(kBlue),Components("McSignalPdfPfromData"),LineStyle(kDashed));

  McFramePfromData->Draw();

  c2->cd();

  c2->Draw();


  cout << "Num sgn from MC(X)resol fit: "  << McNumSigPFitResol->getVal() << "  from bkg: " << McNumBkgPFitResol->getVal() << endl;
  ResMcCombinedResolFit->Print();
  cout << "Num sgn from data template fit: " << McNumSigPFitFromData->getVal() << "  from bkg: " << McNumBkgPFitFromData->getVal() << endl;
  ResMcCombinedTemplateFit->Print();
}







void extractFromData2( const string tnp_      = "etoTauSCMargNoCracks60",
		       bool doBinned_         = true,
		       const string category_ = "tauAntiEMVA",
		       double cutValue_       = 0.5,
		       const string bin_      = "abseta>1.5",
		       double nBins_          = 18,
		       double xLow_           = 65,
		       double xHigh_          = 120
		       ){
  
  TCanvas *c2 = new TCanvas("fitCanvasTemplate","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  // soup  
  //TFile fsup("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA.root");
  TFile fsup("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Data.root");
  TTree *fullTreeSoup = (TTree*)fsup.Get((tnp_+"/fitter_tree").c_str());

  // signal
  TFile fsgnHighStat("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");
  TTree *fullTreeSgnHighStat = (TTree*)fsgnHighStat.Get((tnp_+"/fitter_tree").c_str());

  // variables
  RooRealVar mcTrue("mcTrue","",0,1);
  RooRealVar matchedID("matchedID","",0,1);
  RooRealVar electronPreIDOutput("electronPreIDOutput","",-1000,1);
  RooRealVar tauAntiEMVA("tauAntiEMVA","",0,1);
  RooRealVar signalPFChargedHadrCands("signalPFChargedHadrCands","",0,10);
  RooRealVar leadPFChargedHadrCandTrackPt("leadPFChargedHadrCandTrackPt","",0,200);
  RooRealVar leadPFCandPt("leadPFCandPt","",0,200);
  RooRealVar pt("pt","",0,200);
  RooRealVar abseta("abseta","",0,10);
  RooRealVar mass("mass","",xLow_,xHigh_);

  mass.setBins(10000,"fft");
  mass.setBins( nBins_ );

  ///////////////// signal pdf from MC ////////////////////
  RooDataSet templateP("templateP","dataset for signal-pass template", RooArgSet(mass,mcTrue,tauAntiEMVA,signalPFChargedHadrCands,abseta,matchedID,leadPFChargedHadrCandTrackPt,pt), Import( *fullTreeSgnHighStat ), Cut(Form("(leadPFChargedHadrCandTrackPt>3 && mcTrue && %s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) ) );

  RooDataHist templateHistP("templateHistP","",RooArgSet(mass), templateP, 1.0);
  RooHistPdf TemplateSignalPdfPHist("TemplateSignalPdfPHist","",RooArgSet(mass),templateHistP);

  RooRealVar McMeanResP("McMeanResP","",0,-10,10);
  RooRealVar McSigmaResP("McSigmaResP","",0.5,0,10);
  RooGaussian McResolModP("McResolModP","",mass,McMeanResP,McSigmaResP);
  RooFFTConvPdf McSignalPdfPHist("McSignalPdfPHist","",mass,TemplateSignalPdfPHist,McResolModP);

  // define the functional form for the passing signal in the template
  // breit-wigner
  RooRealVar mean("mean","mean",91.19,85,98);
  RooRealVar width("width","width",2.49,0,20);
  RooBreitWigner bw("bw","bw",mass,mean,width);

  // constrained crystall ball
  RooRealVar m1("m1","m1",0,-20,20);
  RooRealVar sigma("sigma","sigma",0.5,0,10);
  RooRealVar alfa("alfa","alfa", 0.5,0,10);
  RooRealVar n("n","n", 1,1e-06,10);
  RooCBShape cb("cb","",mass,m1,sigma,alfa,n);

  // convolute
  RooFFTConvPdf TemplateSignalPdfP("TemplateSignalPdfP","",mass,bw,cb);
  RooFitResult* ResTemplateFit = TemplateSignalPdfP.fitTo(templateP, Minos(1), Save(1), NumCPU(4) );
  RooArgSet FitParamTempl(ResTemplateFit->floatParsFinal());
  RooRealVar* meanFit    = (RooRealVar*)(&FitParamTempl["mean"]);
  RooRealVar* widthFit   = (RooRealVar*)(&FitParamTempl["width"]);
  //RooRealVar* m1Fit    = (RooRealVar*)(&FitParamTempl["m1"]);
  //RooRealVar* sigmaFit = (RooRealVar*)(&FitParamTempl["sigma"]);
  RooRealVar* alfaFit  = (RooRealVar*)(&FitParamTempl["alfa"]);
  RooRealVar* nFit     = (RooRealVar*)(&FitParamTempl["n"]);

  //RooRealVar mean_C("mean_C","mean",meanFit->getVal());
  //RooRealVar width_C("width_C","width",widthFit->getVal());
  RooRealVar mean_C("mean_C","mean",91,80,110);
  RooRealVar width_C("width_C","width",2.49,0,20);
  RooBreitWigner bw_C("bw_C","bw",mass,mean_C,width_C);

  // constrained crystall ball
  RooRealVar m1_C("m1_C","m1",0.,-20,20);
  RooRealVar sigma_C("sigma_C","sigma",0.5,0,10);
  RooRealVar alfa_C("alfa_C","alfa",alfaFit->getVal());
  RooRealVar n_C("n_C","n", nFit->getVal());
  RooCBShape cb_C("cb_C","",mass,m1_C,sigma_C,alfa_C,n_C);

  // convolute
  RooFFTConvPdf McSignalPdfP("McSignalPdfP","",mass,bw_C,cb_C);


  // bkg pdf
  RooRealVar McCP("McCP","",0,-10,10);
  RooExponential McBackgroundPdfP("McBackgroundPdfP","",mass,McCP);

  RooRealVar McNumBkgP("McNumBkgP","",0,1000000);
  RooRealVar McNumSgnP("McNumSgnP","",0,1000000);
  
  // pdf for combined sample
  RooAddPdf McModelP("McModelP","",RooArgList(McBackgroundPdfP,McSignalPdfP),RooArgList(McNumBkgP,McNumSgnP));
  RooAddPdf McModelPHist("McModelPHist","",RooArgList(McBackgroundPdfP,McSignalPdfPHist),RooArgList(McNumBkgP,McNumSgnP));

  // passing dataset for the soup
  TFile *McP = new TFile("dummy1.root","RECREATE");
  TTree* treeSoupP = fullTreeSoup->CopyTree( Form("(leadPFChargedHadrCandTrackPt>3 && %s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) );
  RooDataSet  McDataPSet("McDataPSet","dataset pass for the soup pass", RooArgSet(mass), Import( *treeSoupP ) );
  RooDataHist McDataP("McDataP","dataset pass for the soup pass", RooArgSet(mass), McDataPSet, 1.0 );

  RooFitResult* ResMcCombinedFit = 0;
  if(doBinned_) ResMcCombinedFit = McModelP.fitTo(McDataP , Extended(1), Minos(1), Save(1), NumCPU(4) );
  else  ResMcCombinedFit = McModelP.fitTo(McDataPSet , Extended(1), Minos(1), Save(1), NumCPU(4) );
  RooArgSet McFitParam(ResMcCombinedFit->floatParsFinal());
  RooRealVar* McNumSigPFit  = (RooRealVar*)(&McFitParam["McNumSgnP"]);
  RooRealVar* McNumBkgPFit  = (RooRealVar*)(&McFitParam["McNumBkgP"]);

  c2->Divide(1,2);

  c2->cd(1);
  RooPlot* TemplateFrameP = mass.frame(Bins(nBins_),Title("MC: passing template"));
  templateP.plotOn(TemplateFrameP);
  TemplateSignalPdfP.plotOn(TemplateFrameP,LineColor(kRed)); 
 
  TemplateFrameP->Draw();

  c2->cd(2);
  RooPlot* McFrameP = mass.frame(Bins(nBins_),Title("MC: passing sample"));
  McDataP.plotOn(McFrameP);
  McModelP.plotOn(McFrameP,LineColor(kBlue));
  McModelP.plotOn(McFrameP,LineColor(kRed),Components("McSignalPdfP"),LineStyle(kSolid));
  McModelP.plotOn(McFrameP,LineColor(kGreen),Components("McBackgroundPdfP"),LineStyle(kDashed));

  McFrameP->Draw();

  c2->cd();

  c2->Draw();

  ResTemplateFit->Print();
  ResMcCombinedFit->Print();
}


void extractFromData3(const string tnp_      = "etoTauSCMargNoCracks60",
		      bool doBinned_         = true,
		      const string category_ = "tauAntiEMVA",
		      const string condition_ = ">=",
		      double cutValue_       = 0.5,
		      const string bin_      = "abseta>1.5",
		      const string additionalCut_ = "abseta>-1",
		      double nBins_          = 18,
		      double xLow_           = 65,
		      double xHigh_          = 120,
		      float deltaAlpha_      = 0.0,
		      float deltaN_          = 0.0
		      ){
  
  TCanvas *c2 = new TCanvas("fitCanvasTemplate","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  // soup  
  //TFile fsup("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA.root");
  TFile fsup("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Data.root");
  TTree *fullTreeSoup = (TTree*)fsup.Get((tnp_+"/fitter_tree").c_str());

  // signal
  TFile fsgnHighStat("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");
  TTree *fullTreeSgnHighStat = (TTree*)fsgnHighStat.Get((tnp_+"/fitter_tree").c_str());

  // variables
  RooRealVar mcTrue("mcTrue","",0,1);
  RooRealVar matchedID("matchedID","",0,1);
  RooRealVar electronPreIDOutput("electronPreIDOutput","",-1000,1);
  RooRealVar tauAntiEMVA("tauAntiEMVA","",0,1);
  RooRealVar signalPFChargedHadrCands("signalPFChargedHadrCands","",0,10);
  RooRealVar leadPFChargedHadrCandTrackPt("leadPFChargedHadrCandTrackPt","",0,200);
  RooRealVar leadPFCandPt("leadPFCandPt","",0,200);
  RooRealVar pt("pt","",0,200);
  RooRealVar abseta("abseta","",0,10);
  RooRealVar mass("mass","",xLow_,xHigh_);

  mass.setBins(10000,"fft");

  ///////////////// signal pdf from MC ////////////////////
  RooDataSet templateP("templateP","dataset for signal-pass template", RooArgSet(mass,mcTrue,tauAntiEMVA,signalPFChargedHadrCands,abseta,leadPFChargedHadrCandTrackPt,pt,matchedID,electronPreIDOutput), Import( *fullTreeSgnHighStat ), Cut(Form("(mcTrue && %s%s%f && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str()) ) );
  RooDataHist templateHistP("templateHistP","",RooArgSet(mass), templateP, 1.0);
  // define the functional form for the passing signal in the template
  // breit-wigner
  RooRealVar mean("mean","mean",91.19);
  RooRealVar width("width","width",2.49);
  RooBreitWigner bw("bw","bw",mass,mean,width);

  // constrained crystall ball
  RooRealVar m1("m1","m1",0,-20,20);
  RooRealVar sigma("sigma","sigma",0.5,0,20);
  RooRealVar alfa("alfa","alfa", 0.5,0,20);
  RooRealVar n("n","n", 1,1e-06,50);
  RooCBShape cb("cb","",mass,m1,sigma,alfa,n);

  RooRealVar sigmaR("sigmaR","sigmaR",0.0);
  RooRealVar sigmaL("sigmaL","sigmaL",0.5,0.,20.);
  RooBifurGauss bifurc("bifurc","bifurc",mass,m1,sigmaL,sigmaR); 
  RooRealVar f("f","",0.5,0,1);
  RooAddPdf cbPlusBifurc("cbPlusBifurc","",RooArgList(cb,bifurc),f);

  // convolute
  RooFFTConvPdf TemplateSignalPdfP("TemplateSignalPdfP","",mass,bw,cb/*PlusBifurc*/);
  RooFitResult* ResTemplateFit = TemplateSignalPdfP.fitTo(templateP, Minos(1), Save(1), NumCPU(4) );
  RooArgSet FitParamTempl(ResTemplateFit->floatParsFinal());
  //RooRealVar* meanFit    = (RooRealVar*)(&FitParamTempl["mean"]);
  //RooRealVar* widthFit   = (RooRealVar*)(&FitParamTempl["width"]);
  RooRealVar* m1Fit    = (RooRealVar*)(&FitParamTempl["m1"]);
  RooRealVar* sigmaFit = (RooRealVar*)(&FitParamTempl["sigma"]);
  RooRealVar* alfaFit  = (RooRealVar*)(&FitParamTempl["alfa"]);
  RooRealVar* nFit     = (RooRealVar*)(&FitParamTempl["n"]);
  //RooRealVar* sigmaLFit= (RooRealVar*)(&FitParamTempl["sigmaL"]);
  //RooRealVar* fFit     = (RooRealVar*)(&FitParamTempl["f"]);

  //RooRealVar mean_C("mean_C","mean",meanFit->getVal());
  //RooRealVar width_C("width_C","width",widthFit->getVal());
  RooBreitWigner bw_C("bw_C","bw",mass,mean/*_C*/,width/*_C*/);

  // constrained crystall ball
  RooRealVar m1_C("m1_C","m1",/*m1Fit->getVal()*/0,-10,10);
  RooRealVar sigma_C("sigma_C","sigma",/*sigmaFit->getVal()*/0.5,0,10);
  RooRealVar alfa_tailC("alfa_tailC","alfa",alfaFit->getVal()*(1+deltaAlpha_));
  RooRealVar n_tailC("n_tailC","n",nFit->getVal()*(1+deltaN_));
  RooRealVar alfa_C("alfa_C","alfa", 0.5,0,20);
  RooRealVar n_C("n_C","n",1,1e-06,50);

  RooCBShape cb_tailC("cb_tailC","",mass,m1_C,sigma_C,alfa_tailC,n_tailC);
  RooCBShape cb_C("cb_C","",mass,m1_C,sigma_C,alfa_C,n_C);
  
  //RooRealVar sigmaL_C("sigmaL_C","sigmaL",sigmaLFit->getVal());
  //RooBifurGauss bifurc_C("bifurc_C","bifurc",mass,m1_C,sigmaL_C,sigmaR); 
  //RooRealVar f_C("f_C","",fFit->getVal());
  //RooAddPdf cbPlusBifurc_C("cbPlusBifurc_C","",RooArgList(cb_C,bifurc_C),f_C);

  // convolute
  RooFFTConvPdf TemplateSignalPdfP_C("TemplateSignalPdfP_C","",mass,bw_C,cb_tailC /*PlusBifurc*/);
  
  RooRealVar McMeanResP("McMeanResP","",0,-10,10);
  RooRealVar McSigmaResP("McSigmaResP","",0.5,0,10);
  RooGaussian McResolModP("McResolModP","",mass,McMeanResP,McSigmaResP);
  RooFFTConvPdf McSignalPdfP("McSignalPdfP","",mass,TemplateSignalPdfP_C,McResolModP);

  // bkg pdf
  RooRealVar McCP("McCP","",-0.5,-10,0.);
  RooExponential McBackgroundPdfP("McBackgroundPdfP","",mass,McCP);

  RooRealVar McNumBkgP("McNumBkgP","",0,1000000);
  RooRealVar McNumSgnP("McNumSgnP","",0,1000000);
  
  // pdf for combined sample
  RooAddPdf McModelP("McModelP","",RooArgList(McBackgroundPdfP,TemplateSignalPdfP_C),RooArgList(McNumBkgP,McNumSgnP));
  
  // passing dataset for the soup
  TFile *McP = new TFile("dummy1.root","RECREATE");
  TTree* treeSoupP = fullTreeSoup->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  RooDataSet  McDataPSet("McDataPSet","dataset pass for the soup pass", RooArgSet(mass), Import( *treeSoupP ) );

  mass.setBins( nBins_ );

  RooDataHist McDataP("McDataP","dataset pass for the soup pass", RooArgSet(mass), McDataPSet, 1.0 );


  RooFitResult* ResMcCombinedFit = 0;
  if(doBinned_) ResMcCombinedFit = McModelP.fitTo(McDataP , Extended(1), Minos(1), Save(1), NumCPU(4) );
  else  ResMcCombinedFit = McModelP.fitTo(McDataPSet , Extended(1), Minos(1), Save(1), NumCPU(4) );
  RooArgSet McFitParam(ResMcCombinedFit->floatParsFinal());
  RooRealVar* McNumSigPFit  = (RooRealVar*)(&McFitParam["McNumSgnP"]);
  RooRealVar* McNumBkgPFit  = (RooRealVar*)(&McFitParam["McNumBkgP"]);

  c2->Divide(1,2);

  c2->cd(1);
  RooPlot* TemplateFrameP = mass.frame(Bins(nBins_),Title("MC: passing template"));
  templateP.plotOn(TemplateFrameP);
  TemplateSignalPdfP.plotOn(TemplateFrameP,LineColor(kRed)); 
  alfa.setVal(alfaFit->getVal()*(1+deltaAlpha_));
  n.setVal(nFit->getVal()*(1+deltaN_));
  TemplateSignalPdfP.plotOn(TemplateFrameP,LineColor(kBlue)); 

  TemplateFrameP->Draw();

  c2->cd(2);
  RooPlot* McFrameP = mass.frame(Bins(nBins_),Title("MC: passing sample"));
  McDataP.plotOn(McFrameP);
  McModelP.plotOn(McFrameP,LineColor(kBlue));
  McModelP.plotOn(McFrameP,LineColor(kRed),Components("TemplateSignalPdfP_C"),LineStyle(kSolid));
  McModelP.plotOn(McFrameP,LineColor(kGreen),Components("McBackgroundPdfP"),LineStyle(kDashed));

  McFrameP->Draw();

  c2->cd();

  c2->Draw();

  ResTemplateFit->Print();
  ResMcCombinedFit->Print();
}






void plotPerturbatedShape(const string tnp_      = "etoTauSCMargNoCracks60",
			  bool doBinned_         = true,
			  const string category_ = "tauAntiEMVA",
			  const string condition_ = ">=",
			  double cutValue_       = 0.5,
			  const string bin_      = "abseta>1.5",
			  const string additionalCut_ = "abseta>-1",
			  double nBins_          = 18,
			  double xLow_           = 65,
			  double xHigh_          = 120,
			  float deltaAlpha_      = 0.0,
			  float deltaN_          = 0.0
			  ){
  
  TCanvas *c2 = new TCanvas("fitCanvasTemplate","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  // signal
  TFile fsgnHighStat("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");
  TTree *fullTreeSgnHighStat = (TTree*)fsgnHighStat.Get((tnp_+"/fitter_tree").c_str());

  // variables
  RooRealVar mcTrue("mcTrue","",0,1);
  RooRealVar matchedID("matchedID","",0,1);
  RooRealVar electronPreIDOutput("electronPreIDOutput","",-1000,1);
  RooRealVar tauAntiEMVA("tauAntiEMVA","",0,1);
  RooRealVar signalPFChargedHadrCands("signalPFChargedHadrCands","",0,10);
  RooRealVar leadPFChargedHadrCandTrackPt("leadPFChargedHadrCandTrackPt","",0,200);
  RooRealVar leadPFCandPt("leadPFCandPt","",0,200);
  RooRealVar pt("pt","",0,200);
  RooRealVar abseta("abseta","",0,10);
  RooRealVar mass("mass","mass",xLow_,xHigh_);

  mass.setBins(10000,"fft");

  TFile *McP = new TFile("dummy1.root","RECREATE");
  TTree* fullTreeSgnHighStatCut = fullTreeSgnHighStat->CopyTree( Form("(mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );

  ///////////////// signal pdf from MC ////////////////////
  RooDataSet templateP("templateP","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeSgnHighStatCut ) );

  // define the functional form for the passing signal in the template
  // breit-wigner
  RooRealVar mean("mean","mean",91.19);
  RooRealVar width("width","width",2.49);
  RooBreitWigner bw("bw","bw",mass,mean,width);

  // constrained crystall ball
  RooRealVar m1("m1","m1",0,-20,20);
  RooRealVar sigma("sigma","sigma",0.5,0,20);
  RooRealVar alfa("alfa","alfa", 0.5,0,20);
  RooRealVar n("n","n", 1,1e-06,50);
  RooCBShape cb("cb","",mass,m1,sigma,alfa,n);

  RooRealVar sigmaR("sigmaR","sigmaR",0.0);
  RooRealVar sigmaL("sigmaL","sigmaL",0.5,0.,20.);
  RooBifurGauss bifurc("bifurc","bifurc",mass,m1,sigmaL,sigmaR); 
  RooRealVar f("f","",0.5,0,1);
  RooAddPdf cbPlusBifurc("cbPlusBifurc","",RooArgList(cb,bifurc),f);

  // convolute
  RooFFTConvPdf TemplateSignalPdfP("TemplateSignalPdfP","",mass,bw,cb);
  RooFitResult* ResTemplateFit = TemplateSignalPdfP.fitTo(templateP, Minos(1), Save(1), NumCPU(4) );
  RooArgSet FitParamTempl(ResTemplateFit->floatParsFinal());
  //RooRealVar* meanFit    = (RooRealVar*)(&FitParamTempl["mean"]);
  //RooRealVar* widthFit   = (RooRealVar*)(&FitParamTempl["width"]);
  RooRealVar* m1Fit    = (RooRealVar*)(&FitParamTempl["m1"]);
  RooRealVar* sigmaFit = (RooRealVar*)(&FitParamTempl["sigma"]);
  RooRealVar* alfaFit  = (RooRealVar*)(&FitParamTempl["alfa"]);
  RooRealVar* nFit     = (RooRealVar*)(&FitParamTempl["n"]);
  //RooRealVar* sigmaLFit= (RooRealVar*)(&FitParamTempl["sigmaL"]);
  //RooRealVar* fFit     = (RooRealVar*)(&FitParamTempl["f"]);

  TLegend* leg = new TLegend(0.6,0.6,0.89,0.89,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);

  TH1F* templ = new TH1F("templ","",1,60,120);
  TH1F* htail = new TH1F("htail","",1,60,120);
  TH1F* ltail = new TH1F("ltail","",1,60,120);
  templ->SetLineColor(kRed);
  htail->SetLineColor(kBlue);
  ltail->SetLineColor(kGreen);
  templ->SetLineStyle(kSolid);
  htail->SetLineStyle(kDashed);
  ltail->SetLineStyle(kDotted);
  templ->SetLineWidth(0.4);
  htail->SetLineWidth(0.4);
  ltail->SetLineWidth(0.4);

  RooPlot* TemplateFrameP = mass.frame(Bins(nBins_),Title("MC: passing template"));
  templateP.plotOn(TemplateFrameP);
  TemplateSignalPdfP.plotOn(TemplateFrameP,LineColor(kRed)); 
  alfa.setVal(alfaFit->getVal()*(1-0.5));
  n.setVal(nFit->getVal()*(1-0.5));
  TemplateSignalPdfP.plotOn(TemplateFrameP,LineColor(kBlue),LineStyle(kDashed)); 
  alfa.setVal(alfaFit->getVal()*(1+2));
  n.setVal(nFit->getVal()*(1+2));
  TemplateSignalPdfP.plotOn(TemplateFrameP,LineColor(kGreen),LineStyle(kDotted)); 

  TemplateFrameP->Draw();

  leg->AddEntry(templ,"best fit to MC","l");
  leg->AddEntry(htail,"#alpha#rightarrow #alpha(1-0.5), n#rightarrow n(1-0.5)","l");
  leg->AddEntry(ltail,"#alpha#rightarrow #alpha(1+2), n#rightarrow n(1+2)","l");

  leg->Draw();

  c2->cd();
  c2->Draw();

}


void bkgTemplateStudy(int nToys_             = 100, 
		      int nTot_              = 440,
		      float purity_          = 0.32,
		      int pdf                = 2, 
		      const string tnp_      = "etoTauSCMargNoCracks80",
		      bool doBinned_         = true,
		      const string category_ = "tauAntiEMVA",
		      const string condition_ = ">=",
		      double cutValue_       = 0.5,
		      const string bin_      = "abseta<1.5",
		      const string additionalCut_ = "abseta>-1",
		      double nBins_          = 18,
		      double xLow_           = 65,
		      double xHigh_          = 120
		      ){
  
  TFile* file = new TFile("bkgShapeStudies.root","UPDATE");

  TCanvas *c2 = new TCanvas("fitCanvasTemplate","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  TH1F* h1 = new TH1F("h1","",40,-4,4);
  // signal
  TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");
  TTree *fullTreeSgn = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());
  // bkg
  TFile fbkg("/data_CMS/cms/lbianchini/tagAndProbe/trees/38XWcut/testNewWriteFromPAT_soup_bkg.root");
  TTree *fullTreeBkg = (TTree*)fbkg.Get((tnp_+"/fitter_tree").c_str());

  RooRealVar mass("mass","mass",40,120);
  //mass.setRange(xLow_,xHigh_);
  mass.setBins( nBins_ );

  TFile *McP = new TFile("dummy1.root","RECREATE");
  TTree* fullTreeSgnCutFakes = fullTreeSgn->CopyTree( Form("(!mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeSgnCutMcTrue = fullTreeSgn->CopyTree( Form("(mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeBkgCut = fullTreeBkg->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  
  //////////////////////////////////////
  //      bkg
  //////////////////////////////////////

  RooDataSet bkgDataSet("bkgDataSet","dataset for bkg-pass template", RooArgSet(mass), Import( *fullTreeBkgCut ) );
  RooDataHist bkgDataHist("bkgDataHist","",RooArgSet(mass),bkgDataSet, 1.0);
  //RooHistPdf  bkgPdf("bkgPdf","",RooArgSet(mass),bkgDataHist);
  
  RooRealVar alpha("alpha","",0,-10,10);
  RooGenericPdf bkgPdf1("bkgPdf1","mass^alpha",RooArgList(mass,alpha));
  if( pdf == 1  ) bkgPdf1.fitTo(bkgDataSet, Minos(1), Save(1), NumCPU(4), Range(xLow_,xHigh_) );

  RooRealVar co("co","",0,-10,10);
  RooExponential bkgPdf2("bkgPdf2","",mass,co);
  if( pdf == 2  ) bkgPdf2.fitTo(bkgDataSet, Minos(1), Save(1), NumCPU(4), Range(xLow_,xHigh_) );
  
  RooKeysPdf bkgPdf3("bkgPdf3","",mass,bkgDataSet);

  RooRealVar c("c","",0,-10,10);
  RooExponential bkgPdfExp("bkgPdfExp","",mass,c);

  mass.setRange(xLow_,xHigh_);

  //////////////////////////////////////
  //      sgn
  //////////////////////////////////////

  RooDataSet sgnDataSet("sgnDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeSgnCutMcTrue ) );
  RooDataHist sgnDataHist("sgnDataHist","",RooArgSet(mass),sgnDataSet, 1.0);
  RooHistPdf  sgnPdf("sgnPdf","",RooArgSet(mass),sgnDataHist);

  RooDataSet sgnDataSetFakes("sgnDataSetFakes","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeSgnCutFakes) );
  RooDataHist sgnDataHistFakes("sgnDataHistFakes","",RooArgSet(mass),sgnDataSetFakes, 1.0);
  RooHistPdf  sgnPdfFakes("sgnPdfFakes","",RooArgSet(mass),sgnDataHistFakes);

  float fake = sgnDataSetFakes.numEntries()/sgnDataSet.numEntries();

  RooRealVar Nsgn("Nsgn","",100,0,1000);
  RooRealVar Nbkg("Nbkg","",100,0,1000);

  RooAddPdf sumExp("sumExp","",RooArgList(sgnPdf,bkgPdfExp),RooArgList(Nsgn,Nbkg));

  int nBkg = nTot_*(1-purity_); int nSgn = nTot_*purity_; int nSgnFake = nSgn*fake;

  RooPlot* frameBkg = mass.frame(  Bins(nBins_) );
  RooPlot* frameSgn = mass.frame(  Bins(nBins_) );
  RooPlot* frameSgnFake = mass.frame(  Bins(nBins_) );
  RooPlot* frameToy = mass.frame(  Bins(nBins_) );

  for(int iToy = 0; iToy< nToys_; iToy++){
  
    RooDataSet iSet("iSet","",mass);
    RooDataSet* iBkg = 0;
    if( pdf == 1  ) iBkg = bkgPdf1.generate( mass, nBkg, Extended(1) );
    else if(pdf==2) iBkg = bkgPdf2.generate( mass, nBkg, Extended(1) );
    else            iBkg = bkgPdf3.generate( mass, nBkg, Extended(1) );
    int iNBkg = iBkg->numEntries();
    RooDataSet* iSgn = sgnPdf.generate( mass, nSgn, Extended(1) );
    int iNSgn = iSgn->numEntries();

    RooDataSet* iSgnFake = sgnPdfFakes.generate( mass, nSgnFake, Extended(1) );
    int iNSgnFake = iSgnFake->numEntries();
    
    iSet.append(*iBkg);
    iSet.append(*iSgn);
    iSet.append(*iSgnFake);

    RooFitResult* fitRes = sumExp.fitTo(iSet, Minos(1), Save(1), NumCPU(4) );
    RooArgSet fitParam(fitRes->floatParsFinal());
    RooRealVar* NsgnFit = (RooRealVar*)(&fitParam["Nsgn"]);

    h1->Fill( (NsgnFit->getVal() - nSgn)/NsgnFit->getError() );

    if(iToy == 0){
      iSet.plotOn(frameToy);
      sumExp.plotOn(frameToy);
      sumExp.plotOn(frameToy,LineColor(kRed),Components("sgnPdf"));
      sumExp.plotOn(frameToy,LineColor(kGreen),Components("bkgPdf"));
    }

  }

  c2->Divide(2,2);

  c2->cd(1);
  sgnDataSet.plotOn(frameSgn);
  sgnPdf.plotOn(frameSgn);
  frameSgn->Draw();
  c2->cd(2);
  bkgDataSet.plotOn(frameBkg);
  if( pdf == 1  )  bkgPdf1.plotOn(frameBkg);
  else if(pdf==2)  bkgPdf2.plotOn(frameBkg);
  else             bkgPdf3.plotOn(frameBkg);
  frameBkg->Draw();
  c2->cd(3);
  sgnDataSetFakes.plotOn(frameSgnFake);
  frameSgnFake->Draw();
  //frameToy->Draw();
  c2->cd(4);
  TF1* gaus = new TF1("gaus","gaus",-4,4);
  h1->Fit(gaus);
  h1->Draw();
  c2->Update();
  c2->Draw();

  file->mkdir( Form("%s%s_%d",tnp_.c_str(),bin_.c_str(),pdf) );
  file->cd(  Form("%s%s_%d",tnp_.c_str(),bin_.c_str(),pdf) );

  c2->Write();

  file->Close();

}




void checkOddBins(int nToys_             = 100, 
		  const string tnp_      = "etoTauMargLooseNoCracks70",
		  const string category_ = "matchedID",
		  const string condition_ = ">=",
		  double cutValue_       = 0.975,
		  const string bin_      = "abseta>1.5",
		  const string additionalCut_ = "abseta>-1",
		  double nBins_          = 18,
		  double xLow_           = 65,
		  double xHigh_          = 120
		  ){
  
  TFile* file = new TFile("checkOddBins.root","RECREATE");

  TCanvas *c2 = new TCanvas("fitCanvasTemplate","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  TH1F* h1 = new TH1F("h1","",40,-10,10);
  // signal
  TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");
  TTree *fullTreeSgn = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());
  // bkg
  TFile fbkg("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA.root");
  TTree *fullTreeBkg = (TTree*)fbkg.Get((tnp_+"/fitter_tree").c_str());

  RooRealVar mass("mass","mass",40,120);
  mass.setRange(xLow_,xHigh_);
  mass.setBins( nBins_ );

  TFile *McP = new TFile("dummy1.root","RECREATE");
  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  float all     = fullTreeSgn->GetEntries(  Form("(mcTrue && %s && %s)",bin_.c_str(),additionalCut_.c_str() ) );
  float passing = fullTreeSgn->GetEntries(  Form("(mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_, bin_.c_str(),additionalCut_.c_str() ) );
  float eff   = passing/all;
  float error = TMath::Sqrt(eff*(1-eff)/all);

  TTree* fullTreeSgnCutMcTrue = fullTreeSgn->CopyTree( Form("(mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeBkgCut = fullTreeBkg->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  int failingSoup = fullTreeBkg->GetEntries(  Form("(!(%s%s%f) && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str() ) );
  
  //////////////////////////////////////
  //      bkg
  //////////////////////////////////////

  RooDataSet bkgDataSet("bkgDataSet","dataset for bkg-pass template", RooArgSet(mass), Import( *fullTreeBkgCut ) );
  RooDataHist bkgDataHist("bkgDataHist","",RooArgSet(mass),bkgDataSet, 1.0);
  RooHistPdf  bkgPdf("bkgPdf","",RooArgSet(mass),bkgDataHist);

  RooRealVar c("c","",0,-10,10);
  RooExponential bkgPdfExp("bkgPdfExp","",mass,c);

  //////////////////////////////////////
  //      sgn
  //////////////////////////////////////

  RooDataSet sgnDataSet("sgnDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeSgnCutMcTrue ) );
  RooDataHist sgnDataHist("sgnDataHist","",RooArgSet(mass),sgnDataSet, 1.0);
  RooHistPdf  sgnPdf("sgnPdf","",RooArgSet(mass),sgnDataHist);

  RooRealVar Nsgn("Nsgn","",100,0,10000);
  RooRealVar Nbkg("Nbkg","",100,0,10000);

  RooAddPdf sumExp("sumExp","",RooArgList(sgnPdf,bkgPdfExp),RooArgList(Nsgn,Nbkg));

  RooPlot* frameBkg = mass.frame(  Bins(nBins_) );
  RooPlot* frameSgn = mass.frame(  Bins(nBins_) );
  RooPlot* frameToy = mass.frame(  Bins(nBins_) );

  for(int iToy = 0; iToy< nToys_; iToy++){
  
    RooDataSet* iBkg = bkgPdf.generate( mass, bkgDataSet.numEntries() , Extended(1) );

    RooFitResult* fitRes = sumExp.fitTo(*iBkg, Minos(1), Save(1), NumCPU(4) );
    RooArgSet fitParam(fitRes->floatParsFinal());
    RooRealVar* NsgnFit = (RooRealVar*)(&fitParam["Nsgn"]);

    h1->Fill( (NsgnFit->getVal()/(failingSoup+NsgnFit->getVal()) - eff)/(NsgnFit->getError()/(failingSoup+NsgnFit->getVal())) );

    if(iToy == 0){
      iBkg->plotOn(frameToy);
      sumExp.plotOn(frameToy);
      sumExp.plotOn(frameToy,LineColor(kRed),Components("sgnPdf"));
      sumExp.plotOn(frameToy,LineColor(kGreen),Components("bkgPdfExp"));
    }

  }

  c2->Divide(2,2);

  c2->cd(1);
  sgnDataSet.plotOn(frameSgn);
  sgnPdf.plotOn(frameSgn);
  frameSgn->Draw();
  c2->cd(2);
  bkgDataSet.plotOn(frameBkg);
  bkgPdf.plotOn(frameBkg);
  frameBkg->Draw();
  c2->cd(3);
  frameToy->Draw();
  c2->cd(4);
  TF1* gaus = new TF1("gaus","gaus",-4,4);
  h1->Fit(gaus);
  h1->Draw();
  c2->Update();
  c2->Draw();

  //file->mkdir( Form("%s%s",tnp_.c_str(),bin_.c_str()) );
  //file->cd(  Form("%s%s",tnp_.c_str(),bin_.c_str()) );

  //c2->Write();

  //file->Close();

}



void bkgStudy(int nToys_             = 10,
	      const string tnp_      = "etoTauSCMargNoCracks80",
	      const string category_ = "tauAntiEMVA",
	      const string condition_ = ">=",
	      double cutValue_       = 0.5,
	      const string bin_      = "abseta<1.5",
	      const string additionalCut_ = "abseta>-1",
	      double nBins_          = 18,
	      double xLow_           = 55,
	      double xHigh_          = 120
	      ){
  
  TCanvas *c2 = new TCanvas("fitCanvasTemplate","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  TH1F* h1 = new TH1F("h1","",40,-4,4);

  // signal
  TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA-PILEUP-NOHLT_500pb.root");
  TTree *fullTreeSgn = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());
  // bkg
  TFile fbkg("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_bkg_tauAntiEMVA.root");
  TTree *fullTreeBkg = (TTree*)fbkg.Get((tnp_+"/fitter_tree").c_str());
  // mix
  TFile fmix("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA.root");
  TTree *fullTreeMix = (TTree*)fmix.Get((tnp_+"/fitter_tree").c_str());

  RooRealVar mass_ctrl("mass","mass",xLow_,xHigh_);
  //RooRealVar mass("mass","mass",50,120);
  RooRealVar mass("mass","mass",xLow_,xHigh_);
  mass.setBins( nBins_ );
  //mass.setRange(xLow_,xHigh_);

  TFile *McP = new TFile("dummy1.root","RECREATE");
  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  int nSgn = fullTreeSgnCut->GetEntries(Form("mass>%f && mass<%f",xLow_,xHigh_));
  TTree* fullTreeBkgCut = fullTreeBkg->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  int nBkg = fullTreeBkgCut->GetEntries(Form("mass>%f && mass<%f",xLow_,xHigh_));
  TTree* fullTreeMixCut = fullTreeMix->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  
  //////////////////////////////////////
  //      bkg
  //////////////////////////////////////

  RooDataSet bkgDataSet("bkgDataSet","dataset for bkg-pass template", RooArgSet(mass_ctrl), Import( *fullTreeBkgCut ) );
  RooDataHist bkgDataHist("bkgDataHist","",RooArgSet(mass_ctrl),bkgDataSet, 1.0);
  RooHistPdf  bkgPdf("bkgPdf","",RooArgSet(mass_ctrl),bkgDataHist);
  
  RooRealVar alpha("alpha","",0,-10,10);
  RooRealVar c0("c0","",10,0,1000);
  //RooGenericPdf bkgPdf1_tail("bkgPdf1_tail","mass^alpha",RooArgList(mass,alpha));
  RooPolynomial bkgPdf1_tail("bkgPdf1_tail","",mass_ctrl,RooArgList());
  RooRealVar mean("mean","",55,40,65);
  RooRealVar sigma("sigma","",10,5,50);
  RooGaussian bkgPdf1_core("bkgPdf1_core","",mass_ctrl,mean,sigma);
  RooRealVar f("f","",0.5,0,1);
  RooAddPdf bkgPdf1("bkgPdf1","",RooArgList( bkgPdf1_tail, bkgPdf1_core),f);
  bkgPdf1.fitTo(bkgDataSet, Minos(1), Save(1), NumCPU(4) );

  RooRealVar co("co","",0,-10,10);
  RooExponential bkgPdf2("bkgPdf2","",mass_ctrl,co);
  bkgPdf2.fitTo(bkgDataSet, Minos(1), Save(1), NumCPU(4) );

  RooRealVar m1_C("m1_C","m1",60,50,70);
  RooRealVar sigma1_C("sigma1_C","sigma",13,5,30);
  RooRealVar alfa1_C("alfa1_C","", -1.,-10.,0.);
  RooRealVar n_C("n_C","n", 1,1e-06,10);
  RooCBShape bkgPdf4_C("bkgPdf4_C","",mass_ctrl,m1_C,sigma1_C,alfa1_C,n_C);
  RooFitResult* fitRes_C = bkgPdf4_C.fitTo(bkgDataSet, Minos(1), Save(1), NumCPU(4) );
  RooArgSet fitParam_C(fitRes_C->floatParsFinal());
  RooRealVar* n_CFit      = (RooRealVar*)(&fitParam_C["n_C"]);
  RooRealVar* alfa1_CFit  = (RooRealVar*)(&fitParam_C["alfa1_C"]);
  
  RooRealVar m1("m1","m1",60,50,70);
  RooRealVar sigma1("sigma1","sigma",13,5,30);
  RooRealVar alfa1("alfa1","",alfa1_CFit->getVal()/*,-10.,0.*/);
  RooRealVar n("n","n", n_CFit->getVal(), 0,10);
  RooCBShape bkgPdf4("bkgPdf4","",mass,m1,sigma1,alfa1,n);
 
  //RooKeysPdf bkgPdf3("bkgPdf3","",mass,bkgDataSet);

  RooRealVar c("c","",0,-10,10);
  RooExponential bkgPdfExp("bkgPdfExp","",mass,c);

  //////////////////////////////////////
  //      sgn
  //////////////////////////////////////

  RooDataSet sgnDataSet("sgnDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeSgnCut ) );
  RooDataHist sgnDataHist("sgnDataHist","",RooArgSet(mass),sgnDataSet, 1.0);
  RooHistPdf  sgnPdf("sgnPdf","",RooArgSet(mass),sgnDataHist);

  RooRealVar Nsgn("Nsgn","",1000,0,100000);
  RooRealVar Nbkg("Nbkg","",1000,0,100000);
  RooAddPdf sumExp("sumExp","",RooArgList(sgnPdf,bkgPdfExp),RooArgList(Nsgn,Nbkg));

  RooPlot* frameBkg = mass_ctrl.frame(  Bins(nBins_) );
  RooPlot* frameSgn = mass.frame(  Bins(nBins_) );
  RooPlot* frameMix = mass.frame(  Bins(nBins_) );
  
  /////////////////////////////////////////
  //      mix
  /////////////////////////////////////////

  RooDataSet mixDataSet("sgnDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeMixCut ) );
  RooDataHist mixDataHist("sgnDataHist","",RooArgSet(mass),mixDataSet, 1.0);

  RooFitResult* fitRes = sumExp.fitTo( mixDataSet, Minos(1), Save(1), NumCPU(4) );
  RooArgSet fitParam(fitRes->floatParsFinal());
  RooRealVar* NsgnFit = (RooRealVar*)(&fitParam["Nsgn"]);


  for(int iToy = 0; iToy< nToys_; iToy++){
    
    RooDataSet iSet("iSet","",mass);
    
    RooDataSet* iBkg = 0;
    RooDataSet* iSgn = 0;
    
    iBkg = bkgPdf4_C.generate( mass,  nBkg , Extended(1) );
    iSgn = sgnPdf.generate(    mass,  nSgn , Extended(1) );

    iSet.append( *iBkg );
    iSet.append( *iSgn );

    //int iNSgn = iSgn->numEntries();
    //int iNSgnFake = iSgnFake->numEntries();
    
    RooFitResult* fitRes = sumExp.fitTo( iSet, Minos(1), Save(1), NumCPU(4) );
    RooArgSet fitParam(fitRes->floatParsFinal());
    RooRealVar* NsgnFit = (RooRealVar*)(&fitParam["Nsgn"]);
    h1->Fill( (NsgnFit->getVal()-iSgn->numEntries())/NsgnFit->getError() );

  } 

  c2->Divide(2,2);

  c2->cd(1);
  bkgDataSet.plotOn(frameBkg);
  bkgPdf1.plotOn(frameBkg,LineColor(kRed));
  bkgPdf4_C.plotOn(frameBkg);
  frameBkg->Draw();
  c2->cd(2);
  mixDataSet.plotOn(frameMix);
  sumExp.plotOn(frameMix,LineColor(kBlue));
  sumExp.plotOn(frameMix,LineColor(kRed),  Components("sgnPdf"));
  sumExp.plotOn(frameMix,LineColor(kGreen),Components("bkgPdfExp"));
  frameMix->Draw();
  c2->cd(3);
  h1->Draw();
  c2->cd(4);
  h1->Draw();

  c2->Update();
  c2->Draw();

  //cout << "Num signal generated = " << sgnDataSet.numEntries() << " ---- Num signl fit = " <<
  //NsgnFit->getVal() << "+/-" << NsgnFit->getError() << endl;

}




void fitStudy(const string tnp_      = "etoTauSCMargNoCracks80",
	      const string category_ = "tauAntiEMVA",
	      const string condition_ = ">=",
	      double cutValue_       = 0.5,
	      const string bin_      = "abseta<1.5",
	      const string additionalCut_ = "abseta>-1",
	      double nBins_          = 18,
	      double xLow_           = 55,
	      double xHigh_          = 120,
	      float deltaM_          = 0.0,
	      bool doBinned_         = true
	      ){
  
  TCanvas *c2 = new TCanvas("fitCanvasTemplate","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  // signal
  TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA-PILEUP-NOHLT.root");
  TTree *fullTreeSgn = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());
  // bkg
  TFile fbkg("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_bkg_tauAntiEMVA.root");
  TTree *fullTreeBkg = (TTree*)fbkg.Get((tnp_+"/fitter_tree").c_str());
  // mix
  TFile fmix("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Data.root");
  TTree *fullTreeMix = (TTree*)fmix.Get((tnp_+"/fitter_tree").c_str());

  RooRealVar mass("mass","mass",xLow_,xHigh_);
  mass.setBins( 10000, "fft" );
  mass.setBins( nBins_ );
 
  TFile *McP = new TFile("dummy1.root","RECREATE");
  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeBkgCut = fullTreeBkg->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeMixCut = fullTreeMix->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  
  //////////////////////////////////////
  //      bkg
  //////////////////////////////////////

  RooDataSet bkgDataSet("bkgDataSet","dataset for bkg-pass template", RooArgSet(mass), Import( *fullTreeBkgCut ) );
  RooDataHist bkgDataHist("bkgDataHist","",RooArgSet(mass),bkgDataSet, 1.0);
  
  // pdf0: RooDataHist
  RooHistPdf  bkgPdf("bkgPdf","",RooArgSet(mass),bkgDataHist);
  
  // pdf00: mass law
  RooRealVar alpha("alpha","",0,-10,10);
  RooGenericPdf("bkgPdf00","mass^alpha",RooArgList(mass,alpha));

  // pdf1: gaussian plus pedestal
  //RooRealVar alpha("alpha","",0,-10,10);
  //RooRealVar c0("c0","",10,0,1000);
  RooPolynomial bkgPdf1_tail_C("bkgPdf1_tail_C","",mass,RooArgList());
  RooRealVar mean_C("mean_C","",55,40,65);
  RooRealVar sigma_C("sigma_C","",10,5,50);
  RooGaussian bkgPdf1_core_C("bkgPdf1_core_C","",mass,mean_C,sigma_C);
  RooRealVar f_C("f_C","",0.5,0,1);
  RooAddPdf bkgPdf1_C("bkgPdf1_C","",RooArgList( bkgPdf1_tail_C, bkgPdf1_core_C),f_C);
  RooFitResult* fitRes1_C = bkgPdf1_C.fitTo(bkgDataSet, Minos(1), Save(1), NumCPU(4) );
  RooArgSet fitParam1_C(fitRes1_C->floatParsFinal());
  RooRealVar* mean_CFit      = (RooRealVar*)(&fitParam1_C["mean_C"]);

  RooPolynomial bkgPdf1_tail("bkgPdf1_tail","",mass,RooArgList());
  RooRealVar mean("mean","",mean_CFit->getVal());
  RooRealVar sigma("sigma","",10,5,50);
  RooGaussian bkgPdf1_core("bkgPdf1_core","",mass,mean,sigma);
  RooRealVar f("f","",0.5,0,1);
  RooAddPdf bkgPdf1("bkgPdf1","",RooArgList( bkgPdf1_tail, bkgPdf1_core),f);

  // pdf2: exponential
  RooRealVar c_C("c_C","",0,-10,10);
  RooExponential bkgPdf2_C("bkgPdf2_C","",mass,c_C);
  bkgPdf2_C.fitTo(bkgDataSet, Minos(1), Save(1), NumCPU(4) );
  RooRealVar c("c","",0,-10,10);
  RooExponential bkgPdf2("bkgPdf2","",mass,c);

  // pdf3: RooKeysPdf
  RooKeysPdf bkgPdf3("bkgPdf3","",mass,bkgDataSet);

  // pdf4: CB
  RooRealVar m1_C("m1_C","m1",60,50,70);
  RooRealVar sigma1_C("sigma1_C","sigma",13,5,30);
  RooRealVar alfa1_C("alfa1_C","", -1.,-10.,0.);
  RooRealVar n_C("n_C","n", 1,1e-06,10);
  RooCBShape bkgPdf4_C("bkgPdf4_C","",mass,m1_C,sigma1_C,alfa1_C,n_C);
  RooFitResult* fitRes4_C = bkgPdf4_C.fitTo(bkgDataSet, Minos(1), Save(1), NumCPU(4) );
  RooArgSet fitParam4_C(fitRes4_C->floatParsFinal());
  RooRealVar* m1_CFit      = (RooRealVar*)(&fitParam4_C["m1_C"]);
  RooRealVar* n_CFit      = (RooRealVar*)(&fitParam4_C["n_C"]);
  RooRealVar* alfa1_CFit  = (RooRealVar*)(&fitParam4_C["alfa1_C"]);
  
  RooRealVar m1("m1","m1",m1_CFit->getVal()*(1+deltaM_),45,65);
  RooRealVar sigma1("sigma1","sigma",13,5,30);
  RooRealVar alfa1("alfa1","",alfa1_CFit->getVal(),-10.,0.);
  RooRealVar n("n","n", n_CFit->getVal(), 0,10);
  RooCBShape bkgPdf4("bkgPdf4","",mass,m1,sigma1,alfa1,n);

  //////////////////////////////////////
  //      sgn
  //////////////////////////////////////

  RooDataSet sgnDataSet("sgnDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeSgnCut ) );
  RooDataHist sgnDataHist("sgnDataHist","",RooArgSet(mass),sgnDataSet, 1.0);
  RooHistPdf  sgnTemplatePdf("sgnTemplatePdf","",RooArgSet(mass),sgnDataHist);

  RooRealVar meanRes("meanRes","",0,-10,10);
  RooRealVar sigmaRes("sigmaRes","",0.5,0,10);

  RooGaussian resolMod("resolModP","",mass,meanRes,sigmaRes);
  RooFFTConvPdf sgnPdf("sgnPdf","",mass,sgnTemplatePdf,resolMod);

  /////////////////////////////////////////
  //      mix
  /////////////////////////////////////////

  RooDataSet mixDataSet("sgnDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeMixCut ) );
  mass.setBins( nBins_ );
  RooDataHist mixDataHist("sgnDataHist","",RooArgSet(mass),mixDataSet, 1.0);

  /////////////////////////////////////////
  //     sum
  /////////////////////////////////////////

  RooRealVar Nsgn("Nsgn","",1000,0,100000);
  RooRealVar Nbkg("Nbkg","",1000,0,100000);

  RooAddPdf sum("sum","",RooArgList(sgnPdf,bkgPdf4),RooArgList(Nsgn,Nbkg));

  RooPlot* frameBkg = mass.frame(  Bins(nBins_) );
  RooPlot* frameSgn = mass.frame(  Bins(nBins_) );
  RooPlot* frameMix = mass.frame(  Bins(nBins_) );

  RooFitResult* fitRes = 0;
  if(doBinned_){
    mass.setBins( 10000,"fft");   
    //mass.setRange(xLow_,xHigh_);
    fitRes = sum.fitTo( mixDataHist,/* Range(xLow_,xHigh_),*/ Minos(1), Save(1), NumCPU(4) );
  }
  else{
    mass.setBins( 10000,"fft");   
    //mass.setRange(xLow_,xHigh_);
    fitRes = sum.fitTo( mixDataSet,/* Range(xLow_,xHigh_),*/ Minos(1), Save(1), NumCPU(4) );
  }
  //RooArgSet fitParam(fitRes->floatParsFinal());

  c2->Divide(1,2);

  c2->cd(1);
  bkgDataSet.plotOn(frameBkg);
  bkgPdf4_C.plotOn(frameBkg,LineColor(kRed));
  //bkgPdf4_C.plotOn(frameBkg);
  frameBkg->Draw();
  c2->cd(2);
  mixDataSet.plotOn(frameMix);
  sum.plotOn(frameMix,LineColor(kBlue));
  sum.plotOn(frameMix,LineColor(kRed),  Components("sgnPdf"));
  sum.plotOn(frameMix,LineColor(kGreen),Components("bkgPdf4"));
  frameMix->Draw();

  c2->Update();
  c2->Draw();

  
}



void fitStudy2(const string tnp_      = "etoTauSCMargNoCracks80",
	      const string category_ = "tauAntiEMVA",
	      const string condition_ = ">=",
	      double cutValue_       = 0.5,
	      const string bin_      = "abseta<1.5",
	      const string additionalCut_ = "abseta>-1",
	      double nBins_          = 18,
	      double xLow_           = 55,
	      double xHigh_          = 120,
	      float deltaM_          = 0.0,
	      bool doBinned_         = true
	      ){
  
  TCanvas *c2 = new TCanvas("fitCanvasTemplate","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  // signal
  TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA-PILEUP-NOHLT.root");
  TTree *fullTreeSgn = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());
  // bkg
  TFile fbkg("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_bkg_tauAntiEMVA.root");
  TTree *fullTreeBkg = (TTree*)fbkg.Get((tnp_+"/fitter_tree").c_str());
  
  // QCD
  TFile fqcd("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_QCD.root");
  TTree *fullTreeQcd = (TTree*)fqcd.Get((tnp_+"/fitter_tree").c_str());
  // bkgQCD
  TFile fbkgQcd("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA_33pb_withQCD.root");
  TTree *fullTreeBkgQcd = (TTree*)fbkgQcd.Get((tnp_+"/fitter_tree").c_str());

  // mix
  TFile fmix("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Data.root");
  TTree *fullTreeMix = (TTree*)fmix.Get((tnp_+"/fitter_tree").c_str());



  RooRealVar mass("mass","mass",xLow_,xHigh_);
  mass.setBins( 10000, "fft" );
  mass.setBins( nBins_ );
 
  TFile *McP = new TFile("dummy1.root","RECREATE");
  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeBkgCut = fullTreeBkg->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeMixCut = fullTreeMix->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeQcdCut = fullTreeQcd->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeBkgQcdCut = fullTreeBkgQcd->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );

  // compute fraction of QCD
  TH1F* hBkgQcd = new TH1F("hBkgQcd","",1,40,140);
  TH1F* hQcd    = new TH1F("hQcd","",1,40,140);
  fullTreeBkgQcdCut->Draw("mass>>hBkgQcd","weight");
  fullTreeQcdCut->Draw("mass>>hQcd","weight");

  float fractionQCD = (float)hQcd->Integral()/(float)hBkgQcd->Integral();

  //////////////////////////////////////
  //      qcd
  //////////////////////////////////////
  mass.setBins( 5 );
  RooDataSet qcdDataSet("qcdDataSet","dataset for qcd", RooArgSet(mass), Import( *fullTreeQcdCut ) );
  RooDataHist qcdDataHist("bkgDataHist","",RooArgSet(mass),qcdDataSet, 1.0);
  RooHistPdf qcdPdf("qcdPdf","",RooArgSet(mass),qcdDataHist);
  mass.setBins( nBins_ );

  //////////////////////////////////////
  //      bkg
  //////////////////////////////////////

  RooDataSet bkgDataSet("bkgDataSet","dataset for bkg-pass template", RooArgSet(mass), Import( *fullTreeBkgCut ) );
  bkgDataSet.append( *(qcdPdf.generate( mass, (int)(fractionQCD*bkgDataSet.numEntries()) )) );

  RooDataHist bkgDataHist("bkgDataHist","",RooArgSet(mass),bkgDataSet, 1.0);
  
  // pdf0: RooDataHist
  RooHistPdf  bkgPdf("bkgPdf","",RooArgSet(mass),bkgDataHist);
  
  // pdf1: gaussian plus pedestal
  //RooRealVar alpha("alpha","",0,-10,10);
  //RooRealVar c0("c0","",10,0,1000);
  RooPolynomial bkgPdf1_tail_C("bkgPdf1_tail_C","",mass,RooArgList());
  RooRealVar mean_C("mean_C","",55,40,65);
  RooRealVar sigma_C("sigma_C","",10,5,50);
  RooGaussian bkgPdf1_core_C("bkgPdf1_core_C","",mass,mean_C,sigma_C);
  RooRealVar f_C("f_C","",0.5,0,1);
  RooAddPdf bkgPdf1_C("bkgPdf1_C","",RooArgList( bkgPdf1_tail_C, bkgPdf1_core_C),f_C);
  RooFitResult* fitRes1_C = bkgPdf1_C.fitTo(bkgDataSet, Minos(1), Save(1), NumCPU(4) );
  RooArgSet fitParam1_C(fitRes1_C->floatParsFinal());
  RooRealVar* mean_CFit      = (RooRealVar*)(&fitParam1_C["mean_C"]);

  RooPolynomial bkgPdf1_tail("bkgPdf1_tail","",mass,RooArgList());
  RooRealVar mean("mean","",mean_CFit->getVal());
  RooRealVar sigma("sigma","",10,5,50);
  RooGaussian bkgPdf1_core("bkgPdf1_core","",mass,mean,sigma);
  RooRealVar f("f","",0.5,0,1);
  RooAddPdf bkgPdf1("bkgPdf1","",RooArgList( bkgPdf1_tail, bkgPdf1_core),f);

  // pdf2: exponential
  RooRealVar c_C("c_C","",0,-10,10);
  RooExponential bkgPdf2_C("bkgPdf2_C","",mass,c_C);
  bkgPdf2_C.fitTo(bkgDataSet, Minos(1), Save(1), NumCPU(4) );
  RooRealVar c("c","",0,-10,10);
  RooExponential bkgPdf2("bkgPdf2","",mass,c);

  // pdf3: RooKeysPdf
  RooKeysPdf bkgPdf3("bkgPdf3","",mass,bkgDataSet);

  // pdf4: CB
  RooRealVar m1_C("m1_C","m1",60,50,70);
  RooRealVar sigma1_C("sigma1_C","sigma",13,5,30);
  RooRealVar alfa1_C("alfa1_C","", -1.,-10.,0.);
  RooRealVar n_C("n_C","n", 1,1e-06,10);
  RooCBShape bkgPdf4_C("bkgPdf4_C","",mass,m1_C,sigma1_C,alfa1_C,n_C);
  RooFitResult* fitRes4_C = bkgPdf4_C.fitTo(bkgDataSet, Minos(1), Save(1), NumCPU(4) );
  RooArgSet fitParam4_C(fitRes4_C->floatParsFinal());
  RooRealVar* m1_CFit      = (RooRealVar*)(&fitParam4_C["m1_C"]);
  RooRealVar* n_CFit      = (RooRealVar*)(&fitParam4_C["n_C"]);
  RooRealVar* alfa1_CFit  = (RooRealVar*)(&fitParam4_C["alfa1_C"]);
  
  RooRealVar m1("m1","m1",m1_CFit->getVal()*(1+deltaM_),45,65);
  RooRealVar sigma1("sigma1","sigma",13,5,30);
  RooRealVar alfa1("alfa1","",alfa1_CFit->getVal(),-10.,0.);
  RooRealVar n("n","n", n_CFit->getVal(), 0,10);
  RooCBShape bkgPdf4("bkgPdf4","",mass,m1,sigma1,alfa1,n);

  //////////////////////////////////////
  //      sgn
  //////////////////////////////////////

  RooDataSet sgnDataSet("sgnDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeSgnCut ) );
  RooDataHist sgnDataHist("sgnDataHist","",RooArgSet(mass),sgnDataSet, 1.0);
  RooHistPdf  sgnTemplatePdf("sgnTemplatePdf","",RooArgSet(mass),sgnDataHist);

  RooRealVar meanRes("meanRes","",0,-10,10);
  RooRealVar sigmaRes("sigmaRes","",0.5,0,10);

  RooGaussian resolMod("resolModP","",mass,meanRes,sigmaRes);
  RooFFTConvPdf sgnPdf("sgnPdf","",mass,sgnTemplatePdf,resolMod);

  /////////////////////////////////////////
  //      mix
  /////////////////////////////////////////

  RooDataSet mixDataSet("sgnDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeMixCut ) );
  mass.setBins( nBins_ );
  RooDataHist mixDataHist("sgnDataHist","",RooArgSet(mass),mixDataSet, 1.0);

  /////////////////////////////////////////
  //     sum
  /////////////////////////////////////////

  RooRealVar Nsgn("Nsgn","",1000,0,100000);
  RooRealVar Nbkg("Nbkg","",1000,0,100000);

  RooAddPdf sum("sum","",RooArgList(sgnPdf,bkgPdf4),RooArgList(Nsgn,Nbkg));

  RooPlot* frameBkg = mass.frame(  Bins(nBins_) );
  RooPlot* frameSgn = mass.frame(  Bins(nBins_) );
  RooPlot* frameMix = mass.frame(  Bins(nBins_) );

  RooFitResult* fitRes = 0;
  if(doBinned_){
    mass.setBins( 10000,"fft");   
    //mass.setRange(xLow_,xHigh_);
    fitRes = sum.fitTo( mixDataHist,/* Range(xLow_,xHigh_),*/ Minos(1), Save(1), NumCPU(4) );
  }
  else{
    mass.setBins( 10000,"fft");   
    //mass.setRange(xLow_,xHigh_);
    fitRes = sum.fitTo( mixDataSet,/* Range(xLow_,xHigh_),*/ Minos(1), Save(1), NumCPU(4) );
  }
  //RooArgSet fitParam(fitRes->floatParsFinal());

  c2->Divide(2,2);

  c2->cd(1);
  bkgDataSet.plotOn(frameBkg);
  bkgPdf4_C.plotOn(frameBkg,LineColor(kRed));
  //bkgPdf4_C.plotOn(frameBkg);
  frameBkg->Draw();
  c2->cd(2);
  mixDataSet.plotOn(frameMix);
  sum.plotOn(frameMix,LineColor(kBlue));
  sum.plotOn(frameMix,LineColor(kRed),  Components("sgnPdf"));
  sum.plotOn(frameMix,LineColor(kGreen),Components("bkgPdf4"));
  frameMix->Draw();
  c2->cd(3);
  qcdDataSet.plotOn(frameSgn);
  qcdPdf.plotOn(frameSgn);
  frameSgn->Draw();
  c2->cd(4);
  frameSgn->Draw();

  c2->Update();
  c2->Draw();

  
}



void fitStudyTemplatesFromMC(const string tnp_      = "etoTauSCMargNoCracks80",
			     const string category_ = "tauAntiEMVA",
			     const string condition_ = ">=",
			     double cutValue_       = 0.5,
			     const string bin_      = "abseta<1.5",
			     const string additionalCut_ = "abseta>-1",
			     double nBins_          = 18,
			     double xLow_           = 40,
			     double xHigh_          = 120,
			     float deltaM_          = 0.0,
			     bool doBinned_         = true
			     ){
  
  TCanvas *c2 = new TCanvas("fitCanvasTemplate","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  // signal
  TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA-PILEUP-NOHLT.root");
  TTree *fullTreeSgn = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());
  
  // bkg
  TFile fbkg("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_bkg_tauAntiEMVA_PILEUP.root");
  TTree *fullTreeBkg = (TTree*)fbkg.Get((tnp_+"/fitter_tree").c_str());
  
  // mix
  //TFile fmix("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA_PILEUP.root");
  TFile fmix("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Data.root");
  TTree *fullTreeMix = (TTree*)fmix.Get((tnp_+"/fitter_tree").c_str());

  // QCD 33 pb
  TFile fqcd("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_QCD.root");
  TTree *fullTreeQcd = (TTree*)fqcd.Get((tnp_+"/fitter_tree").c_str());
  // bkg with QCD 33 pb
  TFile fbkgQcd("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA_33pb_withQCD.root");
  TTree *fullTreeBkgQcd = (TTree*)fbkgQcd.Get((tnp_+"/fitter_tree").c_str());

  // Wjets
  TFile fWen("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_WToENu-PILEUP.root");
  TTree *fullTreeWen = (TTree*)fWen.Get((tnp_+"/fitter_tree").c_str());

  // Z-> tau tau
  TFile fZtt("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToTauTau-PYTHIA-PILEUP.root");
  TTree *fullTreeZtt = (TTree*)fZtt.Get((tnp_+"/fitter_tree").c_str());

  // TTb
  TFile fTTb("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_TT-PILEUP.root");
  TTree *fullTreeTTb = (TTree*)fTTb.Get((tnp_+"/fitter_tree").c_str());
 
  TFile *McP = new TFile("dummy1.root","RECREATE");
  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeBkgCut = fullTreeBkg->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeMixCut = fullTreeMix->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  
  TTree* fullTreeQcdCut = fullTreeQcd->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeBkgQcdCut = fullTreeBkgQcd->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );

  TTree* fullTreeWenCut = fullTreeWen->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZttCut = fullTreeZtt->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeTTbCut = fullTreeTTb->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeSgnFakeCut = fullTreeSgn->CopyTree( Form("(!mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );

  float Lumi_ = 33.;

  // compute fraction of QCD
  TH1F* hBkgQcd = new TH1F("hBkgQcd","",1,40,140);
  TH1F* hQcd    = new TH1F("hQcd","",1,40,140);
  fullTreeBkgQcdCut->Draw("mass>>hBkgQcd","weight");
  fullTreeQcdCut->Draw("mass>>hQcd","weight");
  float fractionQCD = (float)hQcd->Integral()/(float)hBkgQcd->Integral();
  float expQCD = (float)hQcd->Integral()*(Lumi_/33.);

  fWen.cd("allEventsFilter");
  TH1F* totalEventsWen = (TH1F*)gDirectory->Get("totalEvents");
  float readEventsWen = totalEventsWen->GetBinContent(1);
  float expWen = fullTreeWenCut->GetEntries()*Lumi_/(readEventsWen/(7899.*1.32));

  fZtt.cd("allEventsFilter");
  TH1F* totalEventsZtt = (TH1F*)gDirectory->Get("totalEvents");
  float readEventsZtt = totalEventsZtt->GetBinContent(1);
  float expZtt = fullTreeZttCut->GetEntries()*Lumi_/(readEventsZtt/(1300.*1.33));

  fTTb.cd("allEventsFilter");
  TH1F* totalEventsTTb = (TH1F*)gDirectory->Get("totalEvents");
  float readEventsTTb = totalEventsTTb->GetBinContent(1);
  float expTTb = fullTreeTTbCut->GetEntries()*Lumi_/(readEventsTTb/(94.*1.75));

  fsgn.cd("allEventsFilter");
  TH1F* totalEventsSgnFake = (TH1F*)gDirectory->Get("totalEvents");
  float readEventsSgnFake = totalEventsSgnFake->GetBinContent(1);
  float expSgnFake = fullTreeSgnFakeCut->GetEntries()*Lumi_/(readEventsSgnFake/(1300.*1.33));

  McP->cd();
  // mass variable

  RooRealVar mass("mass","mass",xLow_,xHigh_);
  mass.setBins( 10000, "fft" );
  mass.setBins( nBins_ );

  //////////////////////////////////////
  //      qcd
  //////////////////////////////////////
  mass.setBins( 6 );
  RooDataSet qcdDataSet("qcdDataSet","dataset for qcd", RooArgSet(mass), Import( *fullTreeQcdCut ) );
  RooDataHist qcdDataHist("bkgDataHist","",RooArgSet(mass),qcdDataSet, 1.0);
  RooHistPdf qcdPdf("qcdPdf","",RooArgSet(mass),qcdDataHist);

  //////////////////////////////////////
  //      bkg + qcd
  //////////////////////////////////////

  mass.setBins( nBins_ );
  RooDataSet bkgDataSet("bkgDataSet","dataset for bkg-pass template", RooArgSet(mass), Import( *fullTreeBkgCut ) );
  bkgDataSet.append( *(qcdPdf.generate( mass, (int)(fractionQCD*bkgDataSet.numEntries()) )) );
  RooDataHist bkgDataHist("bkgDataHist","",RooArgSet(mass),bkgDataSet, 1.0);

  //////////////////////////////////////
  //      sgn
  //////////////////////////////////////

  mass.setBins( 50 );
  RooDataSet sgnDataSet("sgnDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeSgnCut ) );
  RooDataHist sgnDataHist("sgnDataHist","",RooArgSet(mass),sgnDataSet, 1.0);
  RooHistPdf  sgnTemplatePdf("sgnTemplatePdf","",RooArgSet(mass),sgnDataHist);
  //RooHistPdf  sgnPdf("sgnPdf","",RooArgSet(mass),sgnDataHist);

  // smearing
  RooRealVar meanRes("meanRes","",0,-10,10);
  RooRealVar sigmaRes("sigmaRes","",0.5,0,10);
  RooGaussian resolMod("resolModP","",mass,meanRes,sigmaRes);
  RooFFTConvPdf sgnPdf("sgnPdf","",mass,sgnTemplatePdf,resolMod);

  /////////////////////////////////////////
  //      mix
  /////////////////////////////////////////

  mass.setBins( nBins_ );
  RooDataSet mixDataSet("sgnDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeMixCut ) );
  RooDataHist mixDataHist("sgnDataHist","",RooArgSet(mass),mixDataSet, 1.0);
  RooHistPdf  mixPdf("mixPdf","",RooArgSet(mass),mixDataHist);

  /////////////////////////////////////////
  //      Wen
  /////////////////////////////////////////

  mass.setBins( 8 );
  RooDataSet wenDataSet("wenDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeWenCut ) );
  RooDataHist wenDataHist("wenDataHist","",RooArgSet(mass),wenDataSet, 1.0);
  RooHistPdf  wenPdf("wenPdf","",RooArgSet(mass),wenDataHist,4);

  /////////////////////////////////////////
  //      ZTT
  /////////////////////////////////////////

  mass.setBins( 20 );
  RooDataSet zttDataSet("zttDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeZttCut ) );
  RooDataHist zttDataHist("zttDataHist","",RooArgSet(mass),zttDataSet, 1.0);
  RooHistPdf  zttPdf("zttPdf","",RooArgSet(mass),zttDataHist);

  /////////////////////////////////////////
  //      TTb
  /////////////////////////////////////////

  mass.setBins( 6 );
  RooDataSet ttbDataSet("ttbDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeTTbCut ) );
  RooDataHist ttbDataHist("ttbDataHist","",RooArgSet(mass),ttbDataSet, 1.0);
  RooHistPdf  ttbPdf("ttbPdf","",RooArgSet(mass),ttbDataHist);

  /////////////////////////////////////////
  //      sgn fake
  /////////////////////////////////////////

  mass.setBins( 8 );
  RooDataSet sgnFakeDataSet("sgnFakeDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeSgnFakeCut ) );
  RooDataHist sgnFakeDataHist("sgnFakeDataHist","",RooArgSet(mass),sgnFakeDataSet, 1.0);
  RooHistPdf  sgnFakePdf("sgnFakePdf","",RooArgSet(mass),sgnFakeDataHist,4);

  // bkg pdf
  RooRealVar Nqcd("Nqcd","",100,0,10000);
  RooRealVar Nztt("Nztt","",100,0,10000);
  RooRealVar Nwen("Nwen","",100,0,10000);
  RooRealVar Nttb("Nttb","",100,0,10000);
  RooRealVar NsgnFake("sgnFake","",100,0,10000);

  RooConstVar expQCD_cv("expQCD_cv","",expQCD);
  RooConstVar expZtt_cv("expZtt_cv","",expZtt);
  RooConstVar expWen_cv("expWen_cv","",expWen);
  RooConstVar expTTb_cv("expTTb_cv","",expTTb);
  RooConstVar expSgnFake_cv("expSgnFake_cv","",expSgnFake);
  RooConstVar expQCD_err_cv("expQCD_err_cv","",expQCD*(1+4.0));
  RooConstVar expZtt_err_cv("expZtt_err_cv","",expZtt*(1+.2));
  RooConstVar expWen_err_cv("expWen_err_cv","",expWen*(1+4.0));
  RooConstVar expTTb_err_cv("expTTb_err_cv","",expTTb*(1+.5));
  RooConstVar expSgnFake_err_cv("expSgnFake_err_cv","",expSgnFake*(1+.5));

  RooGaussian NqcdConstraint("NqcdConstraint","",Nqcd,expQCD_cv,expQCD_err_cv) ;
  RooGaussian NzttConstraint("NzttConstraint","",Nztt,expZtt_cv,expZtt_err_cv) ;
  RooGaussian NwenConstraint("NwenConstraint","",Nwen,expWen_cv,expWen_err_cv) ;
  RooGaussian NttbConstraint("NttbConstraint","",Nttb,expTTb_cv,expTTb_err_cv) ;
  RooGaussian NsgnFakeConstraint("NsgnFakeConstraint","",NsgnFake,expSgnFake_cv,expSgnFake_err_cv) ;

  RooAddPdf bkgPdf("bkgPdf","",RooArgList(qcdPdf,zttPdf,wenPdf,ttbPdf,sgnFakePdf),RooArgList(Nqcd,Nztt,Nwen,Nttb,NsgnFake));

  /////////////////////////////////////////
  //     sum
  /////////////////////////////////////////

  RooRealVar Nsgn("Nsgn","",1000,0,100000);
  RooRealVar Nbkg("Nbkg","",1000,0,100000);

  RooAddPdf sum("sum","",RooArgList(sgnPdf,bkgPdf),RooArgList(Nsgn,Nbkg));

  RooPlot* frameBkg     = mass.frame(  Bins(nBins_) );
  RooPlot* frameSgn     = mass.frame(  Bins(nBins_) );
  RooPlot* frameMix     = mass.frame(  Bins(nBins_) );
  RooPlot* frameQcd     = mass.frame(  Bins(nBins_) );
  RooPlot* frameWen     = mass.frame(  Bins(nBins_) );
  RooPlot* frameZtt     = mass.frame(  Bins(nBins_) );
  RooPlot* frameTTb     = mass.frame(  Bins(nBins_) );
  RooPlot* frameSgnFake = mass.frame(  Bins(nBins_) );

  mass.setBins( nBins_ );

  RooFitResult* fitRes = 0;
  if(doBinned_){
    fitRes = sum.fitTo( mixDataHist,ExternalConstraints( RooArgSet(NqcdConstraint,NzttConstraint,NwenConstraint,NttbConstraint,NsgnFakeConstraint) ), /* Range(xLow_,xHigh_),*/ Minos(1), Save(1), NumCPU(4) );
  }
  else{
    fitRes = sum.fitTo( mixDataSet,ExternalConstraints( RooArgSet(NqcdConstraint,NzttConstraint,NwenConstraint,NttbConstraint,NsgnFakeConstraint) ), /* Range(xLow_,xHigh_),*/ Minos(1), Save(1), NumCPU(4)  );
  }
  //RooArgSet fitParam(fitRes->floatParsFinal());

  c2->Divide(2,4);

  c2->cd(1);
  bkgDataSet.plotOn(frameBkg);
  bkgPdf.plotOn(frameBkg,LineColor(kRed));
  frameBkg->Draw();
  c2->cd(2);
  sgnDataSet.plotOn(frameSgn);
  sgnPdf.plotOn(frameSgn);
  frameSgn->Draw();
  c2->cd(3);
  mixDataSet.plotOn(frameMix);
  sum.plotOn(frameMix,LineColor(kBlue));
  sum.plotOn(frameMix,LineColor(kRed),  Components("sgnPdf"));
  sum.plotOn(frameMix,LineColor(kGreen),Components("bkgPdf"));
  frameMix->Draw();
  c2->cd(4);
  qcdDataSet.plotOn(frameQcd);
  qcdPdf.plotOn(frameQcd);
  frameQcd->Draw();
  c2->cd(5);
  wenDataSet.plotOn(frameWen);
  wenPdf.plotOn(frameWen);
  frameWen->Draw();
  c2->cd(6);
  zttDataSet.plotOn(frameZtt);
  zttPdf.plotOn(frameZtt);
  frameZtt->Draw();
  c2->cd(7);
  ttbDataSet.plotOn(frameTTb);
  ttbPdf.plotOn(frameTTb);
  frameTTb->Draw();
  c2->cd(8);
  sgnFakeDataSet.plotOn(frameSgnFake);
  sgnFakePdf.plotOn(frameSgnFake);
  frameSgnFake->Draw();

  c2->Update();
  c2->Draw();

  
  std::cout << "qcd " << expQCD << ", Ztattau " <<expZtt << ", Wenu " << expWen << ", TTbar " << expTTb << ", sgn fake " << expSgnFake << endl; 

}



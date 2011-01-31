#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
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
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooAbsCategory.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "RooTruthModel.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"

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
			 const float xHigh_     = 113,
			 const float nBins_     = 24,
			 bool doBinned_         = true,
			 float deltaAlpha_      = 0.0,
			 float deltaN_          = 0.0,
			 float scale_           = 0.0
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


  // signal
  TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA-PILEUP-NOHLT.root");
  TTree *fullTreeSgn  = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());

  // QCD ~33 pb
  TFile fQcd("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_QCD.root");
  TTree *fullTreeQcd = (TTree*)fQcd.Get((tnp_+"/fitter_tree").c_str());

  // W->e nu
  TFile fWen("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_WToENu-PILEUP.root");
  TTree *fullTreeWen = (TTree*)fWen.Get((tnp_+"/fitter_tree").c_str());

  // Z-> tau tau
  TFile fZtt("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToTauTau-PYTHIA-PILEUP.root");
  TTree *fullTreeZtt = (TTree*)fZtt.Get((tnp_+"/fitter_tree").c_str());

  // TTb
  TFile fTTb("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_TT-PILEUP.root");
  TTree *fullTreeTTb = (TTree*)fTTb.Get((tnp_+"/fitter_tree").c_str());

  // LS/anti-iso data 
  TFile fLS( "/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Run2010B-LSantiISO.root");
  TTree *fullTreeLS = (TTree*)fLS.Get("etoTauMargLooseNoCracks90/fitter_tree");

  // soup
  TFile fmix("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA_PILEUP.root");
  TTree *fullTreeMix = (TTree*)fmix.Get((tnp_+"/fitter_tree").c_str());

  // data
  TFile fdat("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Run2010AB-HLT.root");
  //TFile fdat("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Run2010B_39X.root");
  TTree *fullTreeData = (TTree*)fdat.Get((tnp_+"/fitter_tree").c_str());
  TTree *fullTreeDataForTemplate = (TTree*)fdat.Get("etoTauMargTightNoCracks60/fitter_tree");
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

  // variables to constrct the pdfs
  RooRealVar mcTrue("mcTrue","",0,1);
  RooRealVar matchedID("matchedID","",0,1);
  RooRealVar tauAntiEMVA("tauAntiEMVA","",0,1);
  RooRealVar electronPreIDOutput("electronPreIDOutput","",-1000,1.0);
  RooRealVar signalPFChargedHadrCands("signalPFChargedHadrCands","",0,10);
  RooRealVar leadPFChargedHadrCandTrackPt("leadPFChargedHadrCandTrackPt","",0,200);
  RooRealVar pt("pt","",0,200);
  RooRealVar abseta("abseta","",0,10);

  RooRealVar mass("mass","m_{tp} (GeV/c^{2})",xLow_,xHigh_);
  mass.setBins( 10000, "fft" );
  mass.setBins( nBins_ );
  ////////////////////////////////////////////////


  ////////////////////////////////////////////////
  //    pdfs for bkg/data in passing region:
  ////////////////////////////////////////////////

  // file to copy the trees
  TFile *templFile = new TFile("dummyTempl.root","RECREATE");
  
  TTree* fullTreeSgnCutP = fullTreeSgn->CopyTree( Form("(mcTrue && %s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) );
  TTree* fullTreeSgnCutF = fullTreeSgn->CopyTree( Form("(mcTrue && %s<%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) );
  TTree* fullTreeQcdCutP = fullTreeQcd->CopyTree( Form("(%s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) );
  TTree* fullTreeWenCutP = fullTreeWen->CopyTree( Form("(%s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) );
  TTree* fullTreeZttCutP = fullTreeZtt->CopyTree( Form("(%s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) );
  TTree* fullTreeTTbCutP = fullTreeTTb->CopyTree( Form("(%s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) );
  TTree* fullTreeSgnFakeCutP = fullTreeSgn->CopyTree( Form("(!mcTrue && %s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) );
  TTree* fullTreeLSCutP = fullTreeLS->CopyTree( Form("(%s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) );

  TTree* fullTreeMixCutP = fullTreeMix->CopyTree( Form("(%s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) );
  TTree* fullTreeMixCutF = fullTreeMix->CopyTree( Form("(%s<%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) );
  TTree* fullTreeDataCutP = fullTreeData->CopyTree( Form("(%s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) );
  TTree* fullTreeDataCutTemplP = fullTreeDataForTemplate->CopyTree( Form("(%s>=%f && %s && (leadPFChargedHadrCandTrackPt>25 && leadPFCandPt>15 && signalPFChargedHadrCands<1.5))",category_.c_str(),cutValue_,bin_.c_str()) );
  TTree* fullTreeDataCutF = fullTreeData->CopyTree( Form("(%s<%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) );


  // expected yields

  fQcd.cd();
  float expQcd = fullTreeQcdCutP->GetEntries()/33.;

  fWen.cd("allEventsFilter");
  TH1F* totalEventsWen = (TH1F*)gDirectory->Get("totalEvents");
  float readEventsWen = totalEventsWen->GetBinContent(1);
  float expWen = fullTreeWenCutP->GetEntries()/(readEventsWen/(7899.*1.32));

  fZtt.cd("allEventsFilter");
  TH1F* totalEventsZtt = (TH1F*)gDirectory->Get("totalEvents");
  float readEventsZtt = totalEventsZtt->GetBinContent(1);
  float expZtt = fullTreeZttCutP->GetEntries()/(readEventsZtt/(1300.*1.33));

  fTTb.cd("allEventsFilter");
  TH1F* totalEventsTTb = (TH1F*)gDirectory->Get("totalEvents");
  float readEventsTTb = totalEventsTTb->GetBinContent(1);
  float expTTb = fullTreeTTbCutP->GetEntries()/(readEventsTTb/(94.*1.75));

  fsgn.cd("allEventsFilter");
  TH1F* totalEventsSgnFake = (TH1F*)gDirectory->Get("totalEvents");
  float readEventsSgnFake = totalEventsSgnFake->GetBinContent(1);
  float expSgnFake = fullTreeSgnFakeCutP->GetEntries()/(readEventsSgnFake/(1300.*1.33));
  float readEventsSgn = totalEventsSgnFake->GetBinContent(1);
  float expSgn = fullTreeSgnCutP->GetEntries()/(readEventsSgn/(1300.*1.33));

  templFile->cd();
  //////////////////////////////////////


  //////////////////////////////////////
  //    extraction of the templates:
  //////////////////////////////////////
  


  //////////////////////////////////////
  //      qcd
  //////////////////////////////////////

  mass.setBins( 20 );
  RooDataSet qcdDataSetP("qcdDataSetP","dataset for qcd", RooArgSet(mass), Import( *fullTreeLSCutP ) );
  RooDataHist qcdDataHistP("bkgDataHistP","",RooArgSet(mass),qcdDataSetP, 1.0);
  //RooHistPdf qcdPdfP("qcdPdfP","",RooArgSet(mass),qcdDataHistP);

  // Landau 
  RooRealVar meanLQcd("meanLQcd","",59,40,70);
  RooRealVar sigmaLQcd("sigmaLQcd","",11,5,30);
  RooLandau LQcd("LQcd","",mass,meanLQcd,sigmaLQcd);

  // fit to qcd dataset
  RooFitResult* ResQcdFit = LQcd.fitTo(qcdDataHistP, Minos(1), Save(1), NumCPU(4) );
  RooArgSet FitParamQcd(ResQcdFit->floatParsFinal());
  RooRealVar* meanLQcdFit   = (RooRealVar*)(&FitParamQcd["meanLQcd"]);
  RooRealVar* sigmaLQcdFit  = (RooRealVar*)(&FitParamQcd["sigmaLQcd"]);
  RooConstVar meanLQcd_C("meanLQcd_C","",meanLQcdFit->getVal() );
  RooConstVar sigmaLQcd_C("sigmaLQcd_C","",sigmaLQcdFit->getVal());

  // fitted Landau
  RooLandau qcdPdfP("qcdPdfP","",mass, meanLQcd_C,sigmaLQcd_C);

  /////////////////////////////////////////
  //      Wen
  /////////////////////////////////////////
  
  mass.setBins( 10 );
  RooDataSet wenDataSetP("wenDataSetP","dataset for Wenu", RooArgSet(mass), Import( *fullTreeWenCutP ) );
  RooDataHist wenDataHistP("wenDataHistP","",RooArgSet(mass),wenDataSetP, 1.0);
  //RooHistPdf  wenPdfP("wenPdfP","",RooArgSet(mass),wenDataHist,4);

  // Crystall Ball
  RooRealVar m1Wen("m1Wen","m1",61,50,70);
  RooRealVar sigmaWen("sigmaWen","sigma",12,5,25);
  RooRealVar alfaWen("alfaWen","alfa", -0.5,-5,5);
  RooRealVar nWen("nWen","n", 1,1e-06,10);
  RooCBShape cbWen("cbWen","",mass,m1Wen,sigmaWen,alfaWen,nWen);

  // fit to Wenu dataset
  RooFitResult* ResWenFit = cbWen.fitTo(wenDataSetP, Minos(1), Save(1), NumCPU(4) );
  RooArgSet FitParamWen(ResWenFit->floatParsFinal());
  RooRealVar* m1WenFit    = (RooRealVar*)(&FitParamWen["m1Wen"]);
  RooRealVar* sigmaWenFit = (RooRealVar*)(&FitParamWen["sigmaWen"]);
  RooRealVar* alfaWenFit  = (RooRealVar*)(&FitParamWen["alfaWen"]);
  RooRealVar* nWenFit     = (RooRealVar*)(&FitParamWen["nWen"]);
 
  RooConstVar m1Wen_C("m1Wen_C","",m1WenFit->getVal()*(1+scale_));
  RooConstVar sigmaWen_C("sigmaWen_C","",sigmaWenFit->getVal()*(1+scale_));
  RooConstVar alfaWen_C("alfaWen_C","",alfaWenFit->getVal());
  RooConstVar nWen_C("nWen_C","",nWenFit->getVal());

  // fitted Crystall Ball
  RooCBShape wenPdfP("wenPdfP","",mass,m1Wen_C,sigmaWen_C,alfaWen_C,nWen_C);

  /////////////////////////////////////////
  //      ZTT
  /////////////////////////////////////////

  mass.setBins( 15 );
  RooDataSet zttDataSetP("zttDataSetP","dataset for Ztautau", RooArgSet(mass), Import( *fullTreeZttCutP ) );
  RooDataHist zttDataHistP("zttDataHistP","",RooArgSet(mass),zttDataSetP, 1.0);
  RooHistPdf  zttPdfP("zttPdfP","",RooArgSet(mass),zttDataHistP,4);

  /////////////////////////////////////////
  //      TTb
  /////////////////////////////////////////

  mass.setBins( 1 );
  RooDataSet ttbDataSetP("ttbDataSetP","dataset for TT~", RooArgSet(mass), Import( *fullTreeTTbCutP ) );
  RooDataHist ttbDataHistP("ttbDataHistP","",RooArgSet(mass),ttbDataSetP, 1.0);
  RooHistPdf  ttbPdfP("ttbPdfP","",RooArgSet(mass),ttbDataHistP);

  /////////////////////////////////////////
  //      sgn fake
  /////////////////////////////////////////

  mass.setBins( 8 );
  RooDataSet sgnFakeDataSetP("sgnFakeDataSetP","dataset for fakes in Z->ee", RooArgSet(mass), Import( *fullTreeSgnFakeCutP ) );
  RooDataHist sgnFakeDataHistP("sgnFakeDataHistP","",RooArgSet(mass),sgnFakeDataSetP, 1.0);
  //RooHistPdf  sgnFakePdfP("sgnFakePdfP","",RooArgSet(mass),sgnFakeDataHistP,4);

  // Crystall Ball
  RooRealVar m1sgnFake("m1sgnFake","m1",61,50,70);
  RooRealVar sigmasgnFake("sigmasgnFake","sigma",12,5,25);
  RooRealVar alfasgnFake("alfasgnFake","alfa", -0.5,-5,5);
  RooRealVar nsgnFake("nsgnFake","n", 1,1e-06,10);
  RooCBShape cbsgnFake("cbsgnFake","",mass,m1sgnFake,sigmasgnFake,alfasgnFake,nsgnFake);

  RooFitResult* RessgnFakeFit = cbsgnFake.fitTo(sgnFakeDataSetP, Minos(1), Save(1), NumCPU(4) );
  RooArgSet FitParamsgnFake(RessgnFakeFit->floatParsFinal());
  RooRealVar* m1sgnFakeFit    = (RooRealVar*)(&FitParamsgnFake["m1sgnFake"]);
  RooRealVar* sigmasgnFakeFit = (RooRealVar*)(&FitParamsgnFake["sigmasgnFake"]);
  RooRealVar* alfasgnFakeFit  = (RooRealVar*)(&FitParamsgnFake["alfasgnFake"]);
  RooRealVar* nsgnFakeFit     = (RooRealVar*)(&FitParamsgnFake["nsgnFake"]);
 
  RooConstVar m1sgnFake_C("m1sgnFake_C","",m1sgnFakeFit->getVal()*(1+scale_));
  RooConstVar sigmasgnFake_C("sigmasgnFake_C","",sigmasgnFakeFit->getVal()*(1+scale_));
  RooConstVar alfasgnFake_C("alfasgnFake_C","",alfasgnFakeFit->getVal());
  RooConstVar nsgnFake_C("nsgnFake_C","",nsgnFakeFit->getVal());

  // fitted Crystall Ball
  RooCBShape sgnFakePdfP("sgnFakePdfP","",mass,m1sgnFake_C,sigmasgnFake_C,alfasgnFake_C,nsgnFake_C);


  //////////////////////////////////////
  //      sgn
  //////////////////////////////////////

  mass.setBins( 50 );

  // failing:

  RooDataSet sgnDataSetF("sgnDataSetF","dataset for signal fail", RooArgSet(mass), Import( *fullTreeSgnCutF ) );
  RooDataHist sgnDataHistF("sgnDataHistF","",RooArgSet(mass),sgnDataSetF, 1.0);
  RooHistPdf  sgnPdfF_raw("sgnPdfF_raw","",RooArgSet(mass),sgnDataHistF);

  RooRealVar sgnMeanResF("sgnMeanResF","",0,-10,10);
  RooRealVar sgnSigmaResF("sgnSigmaResF","",0.5,0,10);
  RooGaussian resolModF("sgnResolModF","",mass,sgnMeanResF,sgnSigmaResF);
  RooFFTConvPdf sgnPdfF("sgnPdfF","",mass,sgnPdfF_raw,resolModF);

  // passing:

  RooDataSet sgnDataSetP("sgnDataSetP","dataset for signal", RooArgSet(mass), Import( *fullTreeSgnCutP ) );
  RooDataHist sgnDataHistP("sgnDataHistP","",RooArgSet(mass),sgnDataSetP, 1.0);
  RooHistPdf  sgnTemplatePdfP("sgnTemplatePdfP","",RooArgSet(mass),sgnDataHistP);
  //RooHistPdf  sgnPdfP("sgnPdfP","",RooArgSet(mass),sgnDataHistP);

  mass.setBins( 25 );
  RooDataSet sgnDataSetTemplateFromDataP("sgnDataSetTemplateFromDataP","dataset for signal", RooArgSet(mass), Import( *fullTreeDataCutTemplP ) );
  RooDataHist sgnDataHistTemplateFromDataP("sgnDataHistTemplateFromDataP","",RooArgSet(mass),sgnDataSetTemplateFromDataP, 1.0);
  //RooHistPdf  sgnPdfP("sgnPdfP","",RooArgSet(mass),sgnDataHistTemplateFromDataP);

  mass.setBins( 50 );

  // Breit-Wigner
  RooConstVar meanSgn("meanSgn","mean",91.19);
  RooConstVar widthSgn("widthSgn","width",2.49);
  RooBreitWigner bwSgn("bwSgn","bw",mass,meanSgn,widthSgn);

  // Crystall Ball
  RooRealVar m1Sgn("m1Sgn","m1",0,-20,20);
  RooRealVar sigmaSgn("sigmaSgn","sigma",0.5,0,20);
  RooRealVar alfaSgn("alfaSgn","alfa", 0.5,0,20);
  RooRealVar nSgn("nSgn","n", 1,1e-06,50);
  RooCBShape cbSgn("cbSgn","",mass,m1Sgn,sigmaSgn,alfaSgn,nSgn);

  // BW (X) CB
  RooFFTConvPdf bvcbSgn("bvcbSgn","",mass,bwSgn, cbSgn);
  RooFitResult* ResSgnFit = bvcbSgn.fitTo(sgnDataSetP, Minos(1), Save(1), NumCPU(4) );
  RooArgSet FitParamSgn(ResSgnFit->floatParsFinal());
  RooRealVar* m1SgnFit    = (RooRealVar*)(&FitParamSgn["m1Sgn"]);
  RooRealVar* sigmaSgnFit = (RooRealVar*)(&FitParamSgn["sigmaSgn"]);
  RooRealVar* alfaSgnFit  = (RooRealVar*)(&FitParamSgn["alfaSgn"]);
  RooRealVar* nSgnFit     = (RooRealVar*)(&FitParamSgn["nSgn"]);

  RooRealVar m1Sgn_C("m1Sgn_C","m1",m1SgnFit->getVal(),-10,10);
  RooRealVar sigmaSgn_C("sigmaSgn_C","sigma",sigmaSgnFit->getVal(),0,20);
  // choose to let it float or not
  RooConstVar alfaSgn_C("alfaSgn_C","alfa",alfaSgnFit->getVal()*(1+deltaAlpha_)/*,0,20*/);
  RooConstVar nSgn_C("nSgn_C","n",nSgnFit->getVal()*(1+deltaN_)/*,0,50*/);
  RooCBShape cbSgn_C("cbSgn_C","",mass,m1Sgn_C,sigmaSgn_C,alfaSgn_C,nSgn_C);

  // fitted BW (X) CB 
  RooFFTConvPdf sgnPdfP("sgnPdfP","",mass,bwSgn, cbSgn_C);

  //////////////////////////////////////////
  //  soup -- data
  //////////////////////////////////////////

  mass.setBins( nBins_ );
  RooDataSet mixDataSetP("mixDataSetP","dataset for mix pass", RooArgSet(mass), Import( *fullTreeMixCutP ) );
  mixDataSetP.append( *(qcdPdfP.generate( mass, (int)expQcd*500 )) );
  RooDataHist mixDataHistP("mixDataHistP","",RooArgSet(mass),mixDataSetP, 1.0);
  RooDataSet mixDataSetF("mixDataSetF","dataset for mix fail", RooArgSet(mass), Import( *fullTreeMixCutF ) );
  RooDataHist mixDataHistF("mixDataHistF","",RooArgSet(mass),mixDataSetF, 1.0);

  mass.setBins( nBins_ );
  RooDataSet DataDataSetP("DataDataSetP","dataset for Data pass", RooArgSet(mass), Import( *fullTreeDataCutP ) );
  RooDataHist DataDataHistP("DataDataHistP","",RooArgSet(mass),DataDataSetP, 1.0);
  RooDataSet DataDataSetF("DataDataSetF","dataset for Data fail", RooArgSet(mass), Import( *fullTreeDataCutF ) );
  RooDataHist DataDataHistF("DataDataHistF","",RooArgSet(mass),DataDataSetF, 1.0);

  // intermadiate yields
  RooRealVar Nqcd("Nqcd","",        100,0,10000);
  RooRealVar Nztt("Nztt","",        100,0,10000);
  RooRealVar Nwen("Nwen","",        100,0,10000);
  RooRealVar Nttb("Nttb","",        100,0,10000);
  RooRealVar NsgnFake("NsgnFake","",100,0,10000);

  // all numbers centered around the expected
  RooConstVar Nztt_m("Nztt_m","",1.0);
  RooConstVar Nqcd_m("Nqcd_m","",1.0);
  RooConstVar Nwen_m("Nwen_m","",1.0);
  RooConstVar Nttb_m("Nttb_m","",1.0);
  RooConstVar NsgnFake_m("NsgnFake_m","",1.0);

  // variances are the same for mix and data
  RooConstVar Nqcd_s("Nqcd_s","",         4.0 );    // 4.0 // 4.0
  RooConstVar Nztt_s("Nztt_s","",         0.4 );    // 0.4 // 0.4
  RooConstVar Nwen_s("Nwen_s","",         2.0 );    // 1.0 // 2.0
  RooConstVar Nttb_s("Nttb_s","",         2.0 );    // 2.0 // 1.0 
  RooConstVar NsgnFake_s("NsgnFake_s","", 2.0 );    // 1.0 // 1.0

  // mix 
  RooConstVar McExpQCD_cv("McExpQCD_cv","",expQcd*500.);
  RooConstVar McExpZtt_cv("McExpZtt_cv","",expZtt*500.);
  RooConstVar McExpWen_cv("McExpWen_cv","",expWen*500.);
  RooConstVar McExpTTb_cv("McExpTTb_cv","",expTTb*500.);
  RooConstVar McExpSgnFake_cv("McExpSgnFake_cv","",expSgnFake*500.);

  RooFormulaVar McNzttN("McNzttN","","Nztt/McExpZtt_cv",RooArgSet(Nztt,McExpZtt_cv));
  RooFormulaVar McNqcdN("McNqcdN","Nqcd/McExpQCD_cv",   RooArgSet(Nqcd,McExpQCD_cv));
  RooFormulaVar McNwenN("McNwen","Nwen/McExpWen_cv",    RooArgSet(Nwen,McExpWen_cv));
  RooFormulaVar McNttbN("McNttbN","Nttb/McExpTTb_cv",   RooArgSet(Nttb,McExpTTb_cv));
  RooFormulaVar McNsgnFakeN("McNsgnFakeN","NsgnFake/McExpSgnFake_cv",RooArgSet(NsgnFake,McExpSgnFake_cv));

  RooGaussian McNqcdConstraint("McNqcdConstraint","",McNqcdN,Nqcd_m,Nqcd_s) ;
  RooGaussian McNzttConstraint("McNzttConstraint","",McNzttN,Nztt_m,Nztt_s) ;
  RooGaussian McNwenConstraint("McNwenConstraint","",McNwenN,Nwen_m,Nwen_s) ;
  RooGaussian McNttbConstraint("McNttbConstraint","",McNttbN,Nttb_m,Nttb_s) ;
  RooGaussian McNsgnFakeConstraint("McNsgnFakeConstraint","",McNsgnFakeN,NsgnFake_m,NsgnFake_s) ;


  // data
  RooConstVar DataExpQCD_cv("DataExpQCD_cv","",expQcd*33.);
  RooConstVar DataExpZtt_cv("DataExpZtt_cv","",expZtt*33.);
  RooConstVar DataExpWen_cv("DataExpWen_cv","",expWen*33.);
  RooConstVar DataExpTTb_cv("DataExpTTb_cv","",expTTb*33.);
  RooConstVar DataExpSgnFake_cv("DataExpSgnFake_cv","",expSgnFake*33.);

  RooFormulaVar DataNzttN("DataNzttN","","Nztt/DataExpZtt_cv",RooArgList(Nztt,DataExpZtt_cv));
  RooFormulaVar DataNqcdN("DataNqcdN","Nqcd/DataExpQCD_cv",RooArgSet(Nqcd,DataExpQCD_cv));
  RooFormulaVar DataNwenN("DataNwen","Nwen/DataExpWen_cv",RooArgSet(Nwen,DataExpWen_cv));
  RooFormulaVar DataNttbN("DataNttbN","Nttb/DataExpTTb_cv",RooArgSet(Nttb,DataExpTTb_cv));
  RooFormulaVar DataNsgnFakeN("DataNsgnFakeN","NsgnFake/DataExpSgnFake_cv",RooArgSet(NsgnFake,DataExpSgnFake_cv));

  RooGaussian DataNqcdConstraint("DataNqcdConstraint","",DataNqcdN,Nqcd_m,Nqcd_s) ;
  RooGaussian DataNzttConstraint("DataNzttConstraint","",DataNzttN,Nztt_m,Nztt_s) ;
  RooGaussian DataNwenConstraint("DataNwenConstraint","",DataNwenN,Nwen_m,Nwen_s) ;
  RooGaussian DataNttbConstraint("DataNttbConstraint","",DataNttbN,Nttb_m,Nttb_s) ;
  RooGaussian DataNsgnFakeConstraint("DataNsgnFakeConstraint","",DataNsgnFakeN,NsgnFake_m,NsgnFake_s) ;


  RooRealVar McCF("McCF","",0);
  RooExponential McBackgroundPdfF("McBackgroundPdfF","",mass,McCF);

  ////////////////////////////////////////////////

  RooRealVar DataCF("DataCF","",0);
  RooExponential DataBackgroundPdfF("DataBackgroundPdfF","",mass,DataCF);

  ////////////////////////////////////////////////

 

  RooPlot* TemplateFrameP = mass.frame(Bins(nBins_),Title("Template passing"));
  sgnDataSetP.plotOn(TemplateFrameP);
  bvcbSgn.plotOn(TemplateFrameP);
  
  RooPlot* TemplateFrameF = mass.frame(Bins(nBins_),Title("Template failing"));
  sgnDataSetF.plotOn(TemplateFrameF);
  sgnPdfF.plotOn(TemplateFrameF);

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


  RooRealVar McNumBkgF("McNumBkgF","",0);
  RooRealVar McNumSgn("McNumSgn","",0,1000000);
  RooRealVar McEfficiency("McEfficiency","",0.04,0,1);
  RooFormulaVar McNumSgnP("McNumSgnP","McEfficiency*McNumSgn",RooArgSet(McEfficiency,McNumSgn));
  RooFormulaVar McNumSgnF("McNumSgnF","(1-McEfficiency)*McNumSgn",RooArgSet(McEfficiency,McNumSgn));
 

  RooAddPdf McModelP("McModelP","",RooArgList(sgnPdfP,qcdPdfP,zttPdfP,wenPdfP,ttbPdfP,sgnFakePdfP),RooArgList(McNumSgnP,Nqcd,Nztt,Nwen,Nttb,NsgnFake));
  RooAddPdf McModelF("McModelF","",RooArgList(sgnPdfF,McBackgroundPdfF),RooArgList(McNumSgnF,McNumBkgF));

  RooDataHist McCombData("McCombData","combined data",mass,Index(category),Import("pass", *(mixDataSetP.createHistogram("histoSoupP",mass))) ,Import("fail", *(mixDataSetF.createHistogram("histoSoupF",mass))) ) ;

  RooSimultaneous McSimPdf("McSimPdf","simultaneous pdf",category) ;
  McSimPdf.addPdf(McModelP,"pass") ;
  McSimPdf.addPdf(McModelF,"fail") ;

  RooFitResult* ResMcCombinedFit = McSimPdf.fitTo(McCombData , Extended(1), Minos(1), Save(1), NumCPU(4), ExternalConstraints( RooArgSet(McNqcdConstraint,McNzttConstraint,McNwenConstraint,McNttbConstraint,McNsgnFakeConstraint) ) );
  outFile_->cd(Form("bin%.2f",binCenter_));
  ResMcCombinedFit->Write("McFitResults_Combined");

  
  RooArgSet McFitParam(ResMcCombinedFit->floatParsFinal());
  RooRealVar* McEffFit     = (RooRealVar*)(&McFitParam["McEfficiency"]);
  RooRealVar* McNumSigFit  = (RooRealVar*)(&McFitParam["McNumSgn"]);

  RooPlot* McFrameP = mass.frame(Bins(nBins_),Title("MC: passing sample"));
  McCombData.plotOn(McFrameP,Cut("category==category::pass"));
  McSimPdf.plotOn(McFrameP,Slice(category,"pass"), ProjWData(category,McCombData), LineColor(kBlue));
  McSimPdf.plotOn(McFrameP,Slice(category,"pass"), ProjWData(category,McCombData), Components("sgnPdfP"), LineColor(kRed), LineStyle(kSolid));
  McSimPdf.plotOn(McFrameP,Slice(category,"pass"), ProjWData(category,McCombData), Components("qcdPdfP"), LineColor(kBlack), LineStyle(4));
  McSimPdf.plotOn(McFrameP,Slice(category,"pass"), ProjWData(category,McCombData), Components("zttPdfP"), LineColor(kYellow), LineStyle(5));
  McSimPdf.plotOn(McFrameP,Slice(category,"pass"), ProjWData(category,McCombData), Components("wenPdfP"), LineColor(kGreen), LineStyle(6));
  McSimPdf.plotOn(McFrameP,Slice(category,"pass"), ProjWData(category,McCombData), Components("ttbPdfP"), LineColor(kMagenta), LineStyle(7));
  McSimPdf.plotOn(McFrameP,Slice(category,"pass"), ProjWData(category,McCombData), Components("sgnFakePdfP"), LineColor(kRed), LineStyle(8));

  RooPlot* McFrameF = mass.frame(Bins(nBins_),Title("MC: failing sample"));
  McCombData.plotOn(McFrameF,Cut("category==category::fail"));
  McSimPdf.plotOn(McFrameF,Slice(category,"fail"), ProjWData(category,McCombData), LineColor(kBlue));
  McSimPdf.plotOn(McFrameF,Slice(category,"fail"), ProjWData(category,McCombData), Components("sgnPdfF"), LineColor(kRed), LineStyle(kSolid));

  ////////////////////////////////////////////////
  //                 Data
  ////////////////////////////////////////////////

  RooRealVar DataNumBkgF("DataNumBkgF","",0.01); // was 0.0
  RooRealVar DataNumSgn("DataNumSgn","",0,1000000);
  RooRealVar DataEfficiency("DataEfficiency","",0.04,0,1);

  RooFormulaVar DataNumSgnP("DataNumSgnP","DataEfficiency*DataNumSgn",RooArgSet(DataEfficiency,DataNumSgn));
  RooFormulaVar DataNumSgnF("DataNumSgnF","(1-DataEfficiency)*DataNumSgn",RooArgSet(DataEfficiency,DataNumSgn));
 
  RooAddPdf DataModelP("DataModelP","",RooArgList(sgnPdfP,qcdPdfP,zttPdfP,wenPdfP,ttbPdfP,sgnFakePdfP),RooArgList(DataNumSgnP,Nqcd,Nztt,Nwen,Nttb,NsgnFake));
  RooAddPdf DataModelF("DataModelF","",RooArgList(sgnPdfF,DataBackgroundPdfF),RooArgList(DataNumSgnF,DataNumBkgF));

  mass.setBins( nBins_ );
  // binned combined dataset
  RooDataHist DataCombData("DataCombData","combined data",mass,Index(category),Import("pass", *(DataDataSetP.createHistogram("histoDataP",mass)) ) ,Import("fail", *(DataDataSetF.createHistogram("histoDataF",mass))) ) ;
  // unbinned combined dataset
  RooDataSet DataCombDataUnBinned("DataCombDataUnBinned","combined data",mass,Index(category),Import("pass", DataDataSetP ) ,Import("fail",DataDataSetF) ) ;

  RooSimultaneous DataSimPdf("DataSimPdf","simultaneous pdf",category) ;
  DataSimPdf.addPdf(DataModelP,"pass") ;
  DataSimPdf.addPdf(DataModelF,"fail") ;

  //mass.setBins( 10000, "fft" );
  RooFitResult* ResDataCombinedFit =  0;
  if(doBinned_)  ResDataCombinedFit = DataSimPdf.fitTo(DataCombData , Extended(1), Minos(1), Save(1), NumCPU(4), ExternalConstraints( RooArgSet(DataNqcdConstraint,DataNzttConstraint,DataNwenConstraint,DataNttbConstraint,DataNsgnFakeConstraint) )  );
  else ResDataCombinedFit = DataSimPdf.fitTo(DataCombDataUnBinned , Extended(1), Minos(1), Save(1), NumCPU(4),  ExternalConstraints( RooArgSet(DataNqcdConstraint,DataNzttConstraint,DataNwenConstraint,DataNttbConstraint,DataNsgnFakeConstraint) ) );
  outFile_->cd(Form("bin%.2f",binCenter_));
  ResDataCombinedFit->Write("DataFitResults_Combined");

  
  RooArgSet DataFitParam(ResDataCombinedFit->floatParsFinal());
  RooRealVar* DataEffFit      = (RooRealVar*)(&DataFitParam["DataEfficiency"]);
  RooRealVar* DataNumSigFit   = (RooRealVar*)(&DataFitParam["DataNumSgn"]);

  RooPlot* DataFrameP = mass.frame(Bins(nBins_),Title("Data: passing sample"));
  DataCombData.plotOn(DataFrameP,Cut("category==category::pass"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), LineColor(kBlue));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("sgnPdfP"), LineColor(kRed), LineStyle(kSolid));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("qcdPdfP"), LineColor(kBlack), LineStyle(4));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("zttPdfP"), LineColor(kYellow), LineStyle(5));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("wenPdfP"), LineColor(kGreen), LineStyle(6));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("ttbPdfP"), LineColor(kMagenta), LineStyle(7));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("sgnFakePdfP"), LineColor(kRed), LineStyle(8));


  RooPlot* DataFrameF = mass.frame(Bins(nBins_),Title("Data: failing sample"));
  DataCombData.plotOn(DataFrameF,Cut("category==category::fail"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), LineColor(kBlue));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("sgnPdfF"), LineColor(kRed), LineStyle(kSolid));


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
	      const float xLow_      = 40, 
	      const float xHigh_     = 120,
	      const float nBins_     = 20,
	      bool doBinned_         = false,
	      double bin1_           = 0.0,
	      double bin2_           = 1.5,
	      double bin3_           = 2.5,
	      double yHigh_          = 0.06
	      ){
  
    // output file
  TFile *outFile = new TFile( Form("EtoTauPlots_%s_%s_%s.root",tnp_.c_str(),category_.c_str(),var_.c_str()),"RECREATE");

  cout << "******************** bin1" << endl;
  vector<Double_t*> bin1Results = simFit(outFile,tnp_, category_,cutValue_,
					 string(Form("%s && %s>%f && %s<%f",otherCuts_.c_str(),var_.c_str(),bin1_,var_.c_str(),bin2_)),
					 (bin2_+bin1_)/2,(bin2_-bin1_)/2 
					 ,xLow_, xHigh_, nBins_,doBinned_);

  cout << "******************** bin2" << endl;
  vector<Double_t*> bin2Results = simFit(outFile,tnp_, category_, cutValue_,
					 string(Form("(%s && %s > %f && %s<%f)",otherCuts_.c_str(),var_.c_str(),bin2_,var_.c_str(),bin3_)),
					 (bin3_+bin2_)/2,(bin3_-bin2_)/2 ,
					 xLow_, xHigh_, nBins_,doBinned_);
  
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

  TLegend* leg = new TLegend(0.18,0.6,0.45,0.89,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);

  double binsEdges[3] = {bin1_,bin2_,bin3_};
  TH1F* h1 = new TH1F("h1","",2, binsEdges);
  h1->SetAxisRange(0,yHigh_,"Y");
  h1->SetXTitle( "|#eta|" );
  string YTitle="fake-rate";
  h1->SetYTitle( YTitle.c_str() );
  
  TGraphAsymmErrors* graph_truthMC = new TGraphAsymmErrors(2,truthMC_x,truthMC_y, truthMC_xL,truthMC_xH,truthMC_yL,truthMC_yH);
  graph_truthMC->SetMarkerColor(kBlue);
  graph_truthMC->SetMarkerStyle(kFullCircle);
  graph_truthMC->SetMarkerSize(1.5);
  ///////////////////////////////////////////////////////
  TGraphAsymmErrors* graph_tnpMC = new TGraphAsymmErrors(2,tnpMC_x,tnpMC_y,tnpMC_xL,tnpMC_xH,tnpMC_yL,tnpMC_yH);
  graph_tnpMC->SetMarkerColor(kRed);
  graph_tnpMC->SetMarkerStyle(kFullTriangleDown);
  graph_tnpMC->SetMarkerSize(1.2);
  ///////////////////////////////////////////////////////
  TGraphAsymmErrors* graph_tnpDATA = new TGraphAsymmErrors(2,tnpDATA_x, tnpDATA_y,tnpDATA_xL,tnpDATA_xH,tnpDATA_yL,tnpDATA_yH);
  graph_tnpDATA->SetMarkerColor(kBlack);
  graph_tnpDATA->SetMarkerStyle(kFullTriangleUp);
  graph_tnpDATA->SetMarkerSize(1.5);

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
  if(category_.find("tauAntiEMVA")!=string::npos) discr="passing #xi<-0.1";
  else if(category_.find("electronPreIDOutput")!=string::npos) discr="passing #xi<0.6";
  else if( category_.find("matchedID")!=string::npos) discr="failing WP95";
  string tau = "" ;
  if(tnp_.find("SC")==string::npos) tau="HPS #tau-candidate";
  else tau="Shrinking-Cone #tau-candidate";
  leg->SetHeader(Form("#splitline{CMS Preliminary L=33 pb^{-1}}{%s %s}",tau.c_str(),discr.c_str()));
  leg->AddEntry(h1mcTruth,"MC-truth");
  leg->AddEntry(h1tnpMC,"t&p: simulation");
  leg->AddEntry(h1tnpDATA,"t&p: 7 TeV Data");
  leg->Draw();

  outFile->cd();
  c1->Write();
  outFile->Close();

}



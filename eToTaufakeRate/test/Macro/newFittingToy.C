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
#include "RooChiSquarePdf.h"

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

void fitStudyTemplatesFromMC(const string tnp_      = "etoTauMargLooseNoCracks70",
			     const string category_ = "tauAntiEMVA",
			     const string condition_= ">=",
			     double cutValue_       = 0.5,
			     const string bin_      = "abseta<1.5",
			     const string additionalCut_ = "abseta>-1",
			     int nToys_             = 1, 
			     double nBins_          = 20,
			     double xLow_           = 40,
			     double xHigh_          = 120,
			     float deltaAlpha_      = 0.0,
			     float deltaN_          = 0.0,
			     float scale_           = 0.0, 
			     bool doBinned_         = true
			     ){
  
  TCanvas *c2 = new TCanvas("canvas","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  // signal
  TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA-PILEUP-NOHLT.root");
  //TFile fsgn("/data_CMS/cms/lbianchini/35pb/Htt/testNewWriteFromPAT_DYToEE-PYTHIA-PILEUP-Htt.root");
  //TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");
  TTree *fullTreeSgn = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());
  
  // bkg
  TFile fbkg("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_bkg_tauAntiEMVA_PILEUP.root");
  //TFile fbkg("/data_CMS/cms/lbianchini/35pb/Htt/testNewWriteFromPAT_soup_Htt.root");
  TTree *fullTreeBkg = (TTree*)fbkg.Get((tnp_+"/fitter_tree").c_str());
  
  // mix
  //TFile fmix("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA_PILEUP.root");
  TFile fmix("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Data.root");
  //TFile fmix("/data_CMS/cms/lbianchini/35pb/Htt/testNewWriteFromPAT_Run2010B-HLT-Htt.root");
  TTree *fullTreeMix = (TTree*)fmix.Get((tnp_+"/fitter_tree").c_str());
  TTree *fullTreeMixForTemplate = (TTree*)fmix.Get("etoTauMargTightNoCracks60/fitter_tree");

  // QCD 33 pb
  TFile fqcd("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_QCD.root");
  TTree *fullTreeQcd = (TTree*)fqcd.Get((tnp_+"/fitter_tree").c_str());
  // bkg with QCD 33 pb
  TFile fbkgQcd("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA_33pb_withQCD.root");
  TTree *fullTreeBkgQcd = (TTree*)fbkgQcd.Get((tnp_+"/fitter_tree").c_str());

  // Wjets
  TFile fWen("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_WToENu-PILEUP.root");
  //TFile fWen("/data_CMS/cms/lbianchini/35pb/Htt/testNewWriteFromPAT_WToENu-PILEUP-Htt.root");
  TTree *fullTreeWen = (TTree*)fWen.Get((tnp_+"/fitter_tree").c_str());

  // Z-> tau tau
  TFile fZtt("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToTauTau-PYTHIA-PILEUP.root");
  //TFile fZtt("/data_CMS/cms/lbianchini/35pb/Htt/testNewWriteFromPAT_DYToTauTau-PYTHIA-PILEUP-Htt.root");
  TTree *fullTreeZtt = (TTree*)fZtt.Get((tnp_+"/fitter_tree").c_str());

  // TTb
  TFile fTTb("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_TT-PILEUP.root");
  //TFile fTTb("/data_CMS/cms/lbianchini/35pb/Htt/testNewWriteFromPAT_TT-PILEUP-Htt.root");
  TTree *fullTreeTTb = (TTree*)fTTb.Get((tnp_+"/fitter_tree").c_str());
 
  // LS DATA
  TFile fLS( "/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Run2010B-LSantiISO.root" ,"READ");
  TTree *fullTreeLS = (TTree*)fLS.Get("etoTauMargLooseNoCracks90/fitter_tree");

  TFile *McP = new TFile("dummy1.root","RECREATE");
  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeBkgCut = fullTreeBkg->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeMixCut = fullTreeMix->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeMixCutTempl = fullTreeMixForTemplate->CopyTree( Form("(%s%s%f && %s && (leadPFChargedHadrCandTrackPt>25 && leadPFCandPt>15 && signalPFChargedHadrCands<1.5))",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str()) );  

  TTree* fullTreeQcdCut = fullTreeQcd->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  //TTree* fullTreeQcdCut = fullTreeQcd->CopyTree( Form("(tauAntiEMVA%s%f && %s && %s)",condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeBkgQcdCut = fullTreeBkgQcd->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  //TTree* fullTreeBkgQcdCut = fullTreeBkgQcd->CopyTree( Form("(tauAntiEMVA%s%f && %s && %s)",condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );

  TTree* fullTreeWenCut = fullTreeWen->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZttCut = fullTreeZtt->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeTTbCut = fullTreeTTb->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeSgnFakeCut = fullTreeSgn->CopyTree( Form("(!mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeLSCut = fullTreeLS->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  //TTree* fullTreeLSCut = fullTreeLS->CopyTree( Form("(tauAntiEMVA%s%f && %s && %s)",condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );

  float Lumi_ = 30;

  // compute fraction of QCD
  TH1F* hBkgQcd = new TH1F("hBkgQcd","",1,20,140);
  TH1F* hQcd    = new TH1F("hQcd","",1,20,140);
  fullTreeBkgQcdCut->Draw("mass>>hBkgQcd","weight");
  fullTreeQcdCut->Draw("mass>>hQcd","weight");

  float fractionQCD = (float)hQcd->Integral()/(float)hBkgQcd->Integral();
  float expQCD = (float) (hQcd->Integral()*(Lumi_/33.)*(1+0.2));

  cout<< "fraction " <<  fractionQCD << " of  " << hBkgQcd->Integral() << endl;

  fWen.cd("allEventsFilter");
  TH1F* totalEventsWen = (TH1F*)gDirectory->Get("totalEvents");
  float readEventsWen = totalEventsWen->GetBinContent(1);
  float expWen = fullTreeWenCut->GetEntries()*Lumi_/(readEventsWen/(7899.*1.32));

  fZtt.cd("allEventsFilter");
  TH1F* totalEventsZtt = (TH1F*)gDirectory->Get("totalEvents");
  float readEventsZtt = totalEventsZtt->GetBinContent(1);
  float expZtt = fullTreeZttCut->GetEntries()*Lumi_/(readEventsZtt/(1300.*1.28));

  fTTb.cd("allEventsFilter");
  TH1F* totalEventsTTb = (TH1F*)gDirectory->Get("totalEvents");
  float readEventsTTb = totalEventsTTb->GetBinContent(1);
  float expTTb = fullTreeTTbCut->GetEntries()*Lumi_/(readEventsTTb/(94.*1.75));


  fsgn.cd("allEventsFilter");
  TH1F* totalEventsSgnFake = (TH1F*)gDirectory->Get("totalEvents");
  float readEventsSgnFake = totalEventsSgnFake->GetBinContent(1);
  float expSgnFake = fullTreeSgnFakeCut->GetEntries()*Lumi_/(readEventsSgnFake/(1300.*1.28));

  float readEventsSgn = totalEventsSgnFake->GetBinContent(1);
  float expSgn = fullTreeSgnCut->GetEntries()*Lumi_/(readEventsSgn/(1300.*1.28));

  McP->cd();
  // mass variable

  RooRealVar mass("mass","m_{tp} (GeV/c^{2})",xLow_,xHigh_);
  mass.setBins( 10000, "fft" );
  mass.setBins( nBins_ );

  //////////////////////////////////////
  //      qcd
  //////////////////////////////////////
  mass.setBins( 20 );

  RooDataSet qcdDataSet("qcdDataSet","dataset for qcd", RooArgSet(mass), Import( *fullTreeLSCut ) );
  RooDataHist qcdDataHist("bkgDataHist","",RooArgSet(mass),qcdDataSet, 1.0);
  //RooHistPdf qcdPdf("qcdPdf","",RooArgSet(mass),qcdDataHist);
  RooRealVar meanLQcd("meanLQcd","",59,40,70);
  RooRealVar sigmaLQcd("sigmaLQcd","",11,5,30);
  RooLandau LQcd("LQcd","",mass,meanLQcd,sigmaLQcd);
  RooFitResult* ResQcdFit = LQcd.fitTo(qcdDataHist, Minos(1), Save(1), NumCPU(4) );
  RooArgSet FitParamQcd(ResQcdFit->floatParsFinal());
  RooRealVar* meanLQcdFit   = (RooRealVar*)(&FitParamQcd["meanLQcd"]);
  RooRealVar* sigmaLQcdFit  = (RooRealVar*)(&FitParamQcd["sigmaLQcd"]);
  RooConstVar meanLQcd_C("meanLQcd_C","",meanLQcdFit->getVal() );
  RooConstVar sigmaLQcd_C("sigmaLQcd_C","",sigmaLQcdFit->getVal());

  RooLandau qcdPdf("qcdPdf","",mass, meanLQcd_C,sigmaLQcd_C);

  //////////////////////////////////////
  //      sgn
  //////////////////////////////////////

  mass.setBins( 50 );
  RooDataSet sgnDataSet("sgnDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeSgnCut ) );
  RooDataHist sgnDataHist("sgnDataHist","",RooArgSet(mass),sgnDataSet, 1.0);
  RooHistPdf  sgnTemplatePdf("sgnTemplatePdf","",RooArgSet(mass),sgnDataHist);

  mass.setBins( 25 );
  RooDataSet sgnTemplateDataSet("sgnTemplateDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeMixCutTempl ) );
  RooDataHist sgnTemplateDataHist("sgnTemplateDataHist","",RooArgSet(mass),sgnTemplateDataSet, 1.0);
  //RooHistPdf  sgnPdf("sgnPdf","",RooArgSet(mass),sgnTemplateDataHist);

  mass.setBins( 50 );

  // define the functional form for the passing signal in the template
  // breit-wigner
  RooConstVar meanSgn("meanSgn","mean",91.19);
  RooConstVar widthSgn("widthSgn","width",2.49);
  RooBreitWigner bwSgn("bwSgn","bw",mass,meanSgn,widthSgn);

  // constrained crystall ball
  RooRealVar m1Sgn("m1Sgn","m1",0,-20,20);
  RooRealVar sigmaSgn("sigmaSgn","sigma",0.5,0,20);
  RooRealVar alfaSgn("alfaSgn","alfa", 0.5,0,20);
  RooRealVar nSgn("nSgn","n", 1,1e-06,50);
  RooCBShape cbSgn("cbSgn","",mass,m1Sgn,sigmaSgn,alfaSgn,nSgn);

  RooConstVar sigmaRSgn("sigmaRSgn","sigmaR",0.0);
  RooRealVar sigmaLSgn("sigmaLSgn","sigmaL",0.5,0.,20.);
  RooBifurGauss bifurcSgn("bifurcSgn","bifurc",mass,m1Sgn,sigmaLSgn,sigmaRSgn); 
  RooRealVar fSgn("fSgn","",0.5,0,1);
  RooAddPdf cbPlusBifurcSgn("cbPlusBifurcSgn","",RooArgList(cbSgn,bifurcSgn),fSgn);
  
  // convolute
  RooFFTConvPdf bvcbSgn("bvcbSgn","",mass,bwSgn,/*cbPlusBifurcSgn*/ cbSgn);
  RooFitResult* ResSgnFit = bvcbSgn.fitTo(sgnDataSet, Minos(1), Save(1), NumCPU(4) );
  RooArgSet FitParamSgn(ResSgnFit->floatParsFinal());
  RooRealVar* m1SgnFit    = (RooRealVar*)(&FitParamSgn["m1Sgn"]);
  RooRealVar* sigmaSgnFit = (RooRealVar*)(&FitParamSgn["sigmaSgn"]);
  RooRealVar* alfaSgnFit  = (RooRealVar*)(&FitParamSgn["alfaSgn"]);
  RooRealVar* nSgnFit     = (RooRealVar*)(&FitParamSgn["nSgn"]);

  RooRealVar m1Sgn_C("m1Sgn_C","m1",m1SgnFit->getVal(),-10,10);
  RooRealVar sigmaSgn_C("sigmaSgn_C","sigma",sigmaSgnFit->getVal(),0,20);
  RooRealVar alfaSgn_C("alfaSgn_C","alfa",alfaSgnFit->getVal(),0,20); //0
  RooRealVar nSgn_C("nSgn_C","n",nSgnFit->getVal(),0,50);             //0

  RooLognormal alfaSgn_CPdf("alfaSgn_CPdf","",alfaSgn_C,RooConst(alfaSgnFit->getVal()),RooConst(1.5));//1.5
  RooLognormal nSgn_CPdf("nSgn_CPdf","",nSgn_C,RooConst(nSgnFit->getVal()),RooConst(1.5));//1.5

  RooCBShape cbSgnMc_C("cbSgnMc_C","",mass,/*RooConst(m1SgnFit->getVal())*/m1Sgn_C,/*RooConst(sigmaSgnFit->getVal())*/sigmaSgn_C,RooConst(alfaSgnFit->getVal()),RooConst(nSgnFit->getVal()));
  RooCBShape cbSgn_C("cbSgn_C","",mass,m1Sgn_C,sigmaSgn_C,alfaSgn_C,nSgn_C);
  RooFFTConvPdf sgnPdf("sgnPdf","",mass,bwSgn, cbSgn_C);
  RooFFTConvPdf sgnMcPdf("sgnMcPdf","",mass,bwSgn, cbSgnMc_C);
  
  // smearing
  RooRealVar meanRes("meanRes","",0,-10,10);
  RooRealVar sigmaRes("sigmaRes","",0.5,0,10);
  RooGaussian resolMod("resolModP","",mass,meanRes,sigmaRes);
  //RooFFTConvPdf sgnPdf("sgnPdf","",mass,sgnTemplatePdf,resolMod);


  /////////////////////////////////////////
  //      Wen
  /////////////////////////////////////////

  mass.setBins( 10 );
  RooDataSet wenDataSet("wenDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeWenCut ) );
  RooDataHist wenDataHist("wenDataHist","",RooArgSet(mass),wenDataSet, 1.0);
  //RooHistPdf  wenPdf("wenPdf","",RooArgSet(mass),wenDataHist,4);

  // CB for Wenu
  RooRealVar m1Wen("m1Wen","m1",61,50,70);
  RooRealVar sigmaWen("sigmaWen","sigma",12,5,25);
  RooRealVar alfaWen("alfaWen","alfa", -0.5,-5,5);
  RooRealVar nWen("nWen","n", 1,1e-06,10);
  RooCBShape cbWen("cbWen","",mass,m1Wen,sigmaWen,alfaWen,nWen);

  RooFitResult* ResWenFit = cbWen.fitTo(wenDataSet, Minos(1), Save(1), NumCPU(4) );
  RooArgSet FitParamWen(ResWenFit->floatParsFinal());
  RooRealVar* m1WenFit    = (RooRealVar*)(&FitParamWen["m1Wen"]);
  RooRealVar* sigmaWenFit = (RooRealVar*)(&FitParamWen["sigmaWen"]);
  RooRealVar* alfaWenFit  = (RooRealVar*)(&FitParamWen["alfaWen"]);
  RooRealVar* nWenFit     = (RooRealVar*)(&FitParamWen["nWen"]);
 
  RooConstVar m1Wen_C("m1Wen_C","",m1WenFit->getVal()*(1+scale_));
  RooConstVar sigmaWen_C("sigmaWen_C","",sigmaWenFit->getVal()*(1+scale_));
  RooConstVar alfaWen_C("alfaWen_C","",alfaWenFit->getVal());
  RooConstVar nWen_C("nWen_C","",nWenFit->getVal());

  RooCBShape wenPdf("wenPdf","",mass,m1Wen_C,sigmaWen_C,alfaWen_C,nWen_C);


  /////////////////////////////////////////
  //      ZTT
  /////////////////////////////////////////

  mass.setBins( 15 );
  RooDataSet zttDataSet("zttDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeZttCut ) );

  RooDataSet zttDataSetBias("zttDataSetBias","", RooArgSet(mass) );
  for(int i = 0; i < zttDataSet.numEntries() ; i++){
    const RooArgSet* massSet_i = zttDataSet.get(i);
    RooRealVar* mass_i = (RooRealVar*)massSet_i->find("mass");
    mass_i->setVal(mass_i->getVal()*(1.+scale_));
    zttDataSetBias.add(RooArgSet(*mass_i));
  }
  RooDataHist zttDataHist("zttDataHist","",RooArgSet(mass),zttDataSetBias, 1.0);
  RooHistPdf  zttPdf("zttPdf","",RooArgSet(mass),zttDataHist,4);

  // CB for Ztautau
  RooRealVar m1Ztt("m1Ztt","m1",61,50,70);
  RooRealVar sigmaZtt("sigmaZtt","sigma",12,5,25);
  RooRealVar alfaZtt("alfaZtt","alfa", -0.5,-5,5);
  RooRealVar nZtt("nZtt","n", 1,1e-06,10);
  RooCBShape cbZtt("cbZtt","",mass,m1Ztt,sigmaZtt,alfaZtt,nZtt);

  RooFitResult* ResZttFit = cbZtt.fitTo(zttDataSet, Minos(1), Save(1), NumCPU(4) );
  RooArgSet FitParamZtt(ResZttFit->floatParsFinal());
  RooRealVar* m1ZttFit    = (RooRealVar*)(&FitParamZtt["m1Ztt"]);
  RooRealVar* sigmaZttFit = (RooRealVar*)(&FitParamZtt["sigmaZtt"]);
  RooRealVar* alfaZttFit  = (RooRealVar*)(&FitParamZtt["alfaZtt"]);
  RooRealVar* nZttFit     = (RooRealVar*)(&FitParamZtt["nZtt"]);
 
  RooConstVar m1Ztt_C("m1Ztt_C","",m1ZttFit->getVal()*(1+scale_));
  RooConstVar sigmaZtt_C("sigmaZtt_C","",sigmaZttFit->getVal()*(1+scale_));
  RooConstVar alfaZtt_C("alfaZtt_C","",alfaZttFit->getVal());
  RooConstVar nZtt_C("nZtt_C","",nZttFit->getVal());

  //RooCBShape zttPdf("zttPdf","",mass,m1Ztt_C,sigmaZtt_C,alfaZtt_C,nZtt_C);
  

  /////////////////////////////////////////
  //      TTb
  /////////////////////////////////////////

  mass.setBins( 1 );
  RooDataSet ttbDataSet("ttbDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeTTbCut ) );
  RooDataHist ttbDataHist("ttbDataHist","",RooArgSet(mass),ttbDataSet, 1.0);
  RooHistPdf  ttbPdf("ttbPdf","",RooArgSet(mass),ttbDataHist);

  /////////////////////////////////////////
  //      sgn fake
  /////////////////////////////////////////

  mass.setBins( 8 );
  RooDataSet sgnFakeDataSet("sgnFakeDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeSgnFakeCut ) );
  RooDataHist sgnFakeDataHist("sgnFakeDataHist","",RooArgSet(mass),sgnFakeDataSet, 1.0);
  //RooHistPdf  sgnFakePdf("sgnFakePdf","",RooArgSet(mass),sgnFakeDataHist,4);


  // CB for sgnFake
  RooRealVar m1sgnFake("m1sgnFake","m1",61,50,70);
  RooRealVar sigmasgnFake("sigmasgnFake","sigma",12,5,25);
  RooRealVar alfasgnFake("alfasgnFake","alfa", -0.5,-5,5);
  RooRealVar nsgnFake("nsgnFake","n", 1,1e-06,10);
  RooCBShape cbsgnFake("cbsgnFake","",mass,m1sgnFake,sigmasgnFake,alfasgnFake,nsgnFake);

  RooFitResult* RessgnFakeFit = cbsgnFake.fitTo(sgnFakeDataSet, Minos(1), Save(1), NumCPU(4) );
  RooArgSet FitParamsgnFake(RessgnFakeFit->floatParsFinal());
  RooRealVar* m1sgnFakeFit    = (RooRealVar*)(&FitParamsgnFake["m1sgnFake"]);
  RooRealVar* sigmasgnFakeFit = (RooRealVar*)(&FitParamsgnFake["sigmasgnFake"]);
  RooRealVar* alfasgnFakeFit  = (RooRealVar*)(&FitParamsgnFake["alfasgnFake"]);
  RooRealVar* nsgnFakeFit     = (RooRealVar*)(&FitParamsgnFake["nsgnFake"]);
 
  RooConstVar m1sgnFake_C("m1sgnFake_C","",m1sgnFakeFit->getVal()*(1+scale_));
  RooConstVar sigmasgnFake_C("sigmasgnFake_C","",sigmasgnFakeFit->getVal()*(1+scale_));
  RooConstVar alfasgnFake_C("alfasgnFake_C","",alfasgnFakeFit->getVal());
  RooConstVar nsgnFake_C("nsgnFake_C","",nsgnFakeFit->getVal());

  RooCBShape sgnFakePdf("sgnFakePdf","",mass,m1sgnFake_C,sigmasgnFake_C,alfasgnFake_C,nsgnFake_C);

  ////////////////////////////////////////////////////////////////


  RooRealVar Nsgn("Nsgn","",1000,0,1000000);

  RooRealVar Nqcd("Nqcd","",100,0,10000);
  RooRealVar Nztt("Nztt","",100,0,10000);
  RooRealVar Nwen("Nwen","",100,0,10000);
  RooRealVar Nttb("Nttb","",100,0,10000);
  RooRealVar NsgnFake("NsgnFake","",100,0,10000);

  RooConstVar expQCD_cv("expQCD_cv","",expQCD);
  RooConstVar expZtt_cv("expZtt_cv","",expZtt);
  RooConstVar expWen_cv("expWen_cv","",expWen);
  RooConstVar expTTb_cv("expTTb_cv","",expTTb);
  RooConstVar expSgnFake_cv("expSgnFake_cv","",expSgnFake);

  RooLognormal NqcdConstraint("NqcdConstraint","",Nqcd,expQCD_cv,RooConst(1.5)) ; //1.5
  RooLognormal NzttConstraint("NzttConstraint","",Nztt,expZtt_cv,RooConst(1.2)) ; //1.2
  RooLognormal NwenConstraint("NwenConstraint","",Nwen,expWen_cv,RooConst(2.0)) ; //2.0
  RooLognormal NttbConstraint("NttbConstraint","",Nttb,expTTb_cv,RooConst(2.0)) ; //2.0
  RooLognormal NsgnFakeConstraint("NsgnFakeConstraint","",NsgnFake,expSgnFake_cv,RooConst(2.0)) ; //2.0

  RooAddPdf sum("sum","",RooArgList(sgnPdf,qcdPdf,zttPdf,wenPdf,ttbPdf,sgnFakePdf),RooArgList(Nsgn,Nqcd,Nztt,Nwen,Nttb,NsgnFake));
  RooProdPdf sumTimesConstr("sumTimesConstr","",RooArgList(sum,NqcdConstraint,NzttConstraint,NwenConstraint,NttbConstraint,NsgnFakeConstraint,alfaSgn_CPdf,nSgn_CPdf));
  
  TH1F* h1 = new TH1F("h1","pull Nsgn",160,-4,4);
  TH1F* h3 = new TH1F("h3","pull alfaSgn",80,-4,4);
  TH1F* h4 = new TH1F("h4","pull nSgn",80,-4,4);
  TH1F* h2 = new TH1F("h2","-2ln(L/L^{sat})",100,0,100);

  RooDataSet mixDataSet("mixDataSet", "", RooArgSet(mass) );

  RooFitResult* fitRes = 0;
  for(unsigned int iToy = 1; iToy <= nToys_; iToy++){

    mixDataSet.reset();
    mixDataSet.append( *(bvcbSgn.generate( mass,(int)expSgn , Extended() ))        );
    int nSgnGen = mixDataSet.numEntries();
    mixDataSet.append( *(wenPdf.generate( mass, (int)expWen , Extended()  ))       );
    mixDataSet.append( *(zttPdf.generate( mass ,(int)expZtt , Extended()  ))       );
    mixDataSet.append( *(sgnFakePdf.generate( mass, (int)expSgnFake , Extended() )));
    mixDataSet.append( *(ttbPdf.generate( mass, (int)expTTb , Extended() ))        );
    mixDataSet.append( *(qcdPdf.generate( mass, (int)expQCD , Extended() ))        );

    mass.setBins(nBins_);
    RooDataHist mixDataHist("mixDataHist","",RooArgSet(mass),mixDataSet, 1.0);
    int nEntries = mixDataHist.sumEntries();
    RooHistPdf  mixSaturatedPdf("mixSaturatedPdf","",RooArgSet(mass),mixDataHist);

    // restore intial values
    Nqcd.setVal(expQCD);
    Nwen.setVal(expWen);
    Nztt.setVal(expZtt);
    NsgnFake.setVal(expSgnFake);
    Nttb.setVal(expTTb);
    Nsgn.setVal(expSgn);
    m1Sgn_C.setVal(m1SgnFit->getVal());
    sigmaSgn_C.setVal(sigmaSgnFit->getVal());
    alfaSgn_C.setVal(alfaSgnFit->getVal());
    nSgn_C.setVal(nSgnFit->getVal());

    if(doBinned_){
      fitRes = sum.fitTo( mixDataHist,ExternalConstraints( RooArgSet(NqcdConstraint,NzttConstraint,NwenConstraint,NttbConstraint,NsgnFakeConstraint,alfaSgn_CPdf,nSgn_CPdf) ), Extended(),  Minos(1), Save(1), NumCPU(4) );
      //fitRes = sumTimesConstr.fitTo( mixDataHist, Extended(),  Minos(1), Save(1), NumCPU(4) );
      RooNLLVar nLogLike("nLogLike","",sumTimesConstr,mixDataHist,kTRUE);
      RooNLLVar nLogLikeSaturated("nLogLikeSaturated","",mixSaturatedPdf,mixDataHist);

      h2->Fill(2*(nLogLike.getVal() - (nLogLikeSaturated.getVal()+ nEntries - nEntries*log(nEntries) ) ));

      //std::cout << "- Log likelihhod from RooNLLVar after fit " << nLogLike.getVal() 
      //	<< "- Log likelihhod from RooNLLVarSaturated  " << (nLogLikeSaturated.getVal()+nEntries - nEntries*log(nEntries)) 
      //	<< "----- min - Log likelihood " << fitRes->minNll() << std::endl;
    }
    else{
      fitRes = sum.fitTo( mixDataSet,ExternalConstraints( RooArgSet(NqcdConstraint,NzttConstraint,NwenConstraint,NttbConstraint,NsgnFakeConstraint,alfaSgn_CPdf,nSgn_CPdf) ),  Minos(1), Save(1), Extended(), NumCPU(4)  );
    }

    RooArgSet fitParam(fitRes->floatParsFinal());
    RooRealVar* numSgnFit_i  = (RooRealVar*)(&fitParam["Nsgn"]);
    RooRealVar* alfaSgnFit_i = (RooRealVar*)(&fitParam["alfaSgn_C"]);
    RooRealVar* nSgnFit_i    = (RooRealVar*)(&fitParam["nSgn_C"]);
    h1->Fill( (numSgnFit_i->getVal()  - expSgn)/expSgn/*numSgnFit_i->getError()*/  );
    h3->Fill( (alfaSgnFit_i->getVal() - alfaSgnFit->getVal())/alfaSgnFit_i->getError()  );
    h4->Fill( (nSgnFit_i->getVal()    - nSgnFit->getVal())/nSgnFit_i->getError()  );
      
  }

  RooRealVar chi2("chi2","",0,100);
  RooRealVar ndof("ndof","",0,100);
  RooDataHist chi2Hist("chi2Hist","",chi2,h2,1.0);
  RooChiSquarePdf chi2Pdf("chi2Pdf","",chi2,ndof);
  RooFitResult* chi2fit = chi2Pdf.fitTo(chi2Hist, Save(1));
  RooArgSet fitParamChi2(chi2fit->floatParsFinal());
  RooRealVar* ndofFit  = (RooRealVar*)(&fitParamChi2["ndof"]);
  cout << "fitted ndof " << ndofFit->getVal() << endl;


  c2->Divide(2,2);
  c2->cd(1);
  h1->Draw();
  c2->cd(2);
  h2->Draw();
  c2->cd(3);
  h3->Draw();
  c2->cd(4);
  h4->Draw();


}




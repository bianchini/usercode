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

void fitStudyTemplatesFromMC(const string tnp_      = "etoTauSCMargNoCracks80",
			     const string category_ = "tauAntiEMVA",
			     const string condition_= ">=",
			     double cutValue_       = 0.5,
			     const string bin_      = "abseta<1.5",
			     const string additionalCut_ = "abseta>-1",
			     double nBins_          = 18,
			     double xLow_           = 40,
			     double xHigh_          = 120,
			     float deltaAlpha_      = 0.0,
			     float deltaN_          = 0.0,
			     float scale_           = 0.0, 
			     bool doBinned_         = true,
			     bool isData_           = true,
			     bool fitInFail_        = false
			     ){
  
  TCanvas *c2 = new TCanvas("fitCanvasTemplate2","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  TCanvas *c1 = new TCanvas("fitCanvasTemplate1","canvas",10,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  // signal
  TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA-PILEUP-NOHLT.root");
  //TFile fsgn("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");
  TTree *fullTreeSgn = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());
  
  // bkg
  TFile fbkg("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_bkg_tauAntiEMVA_PILEUP.root");
  TTree *fullTreeBkg = (TTree*)fbkg.Get((tnp_+"/fitter_tree").c_str());
  
  // mix
  //TFile fmix("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA_PILEUP.root");
  TFile fmix("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Data.root");
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
  TTree *fullTreeWen = (TTree*)fWen.Get((tnp_+"/fitter_tree").c_str());

  // Z-> tau tau
  TFile fZtt("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToTauTau-PYTHIA-PILEUP.root");
  TTree *fullTreeZtt = (TTree*)fZtt.Get((tnp_+"/fitter_tree").c_str());

  // TTb
  TFile fTTb("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_TT-PILEUP.root");
  TTree *fullTreeTTb = (TTree*)fTTb.Get((tnp_+"/fitter_tree").c_str());
 
  // LS DATA
  TFile fLS( "/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_Run2010B-LSantiISO.root" ,"READ");
  TTree *fullTreeLS = (TTree*)fLS.Get("etoTauMargLooseNoCracks90/fitter_tree");

  TFile *McP = new TFile("dummy1.root","RECREATE");
  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeBkgCut = fullTreeBkg->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeMixCut = fullTreeMix->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeMixCutTempl = fullTreeMixForTemplate->CopyTree( Form("(%s%s%f && %s && (leadPFChargedHadrCandTrackPt>25 && leadPFCandPt>15 && signalPFChargedHadrCands<1.5))",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str()) );  

  //fullTreeMixCutTempl->Draw("mass");
  //return;

  TTree* fullTreeQcdCut = fullTreeQcd->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeBkgQcdCut = fullTreeBkgQcd->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );

  TTree* fullTreeWenCut = fullTreeWen->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZttCut = fullTreeZtt->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeTTbCut = fullTreeTTb->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeSgnFakeCut = fullTreeSgn->CopyTree( Form("(!mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeLSCut = fullTreeLS->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );

  float Lumi_ = 33.;

  // compute fraction of QCD
  TH1F* hBkgQcd = new TH1F("hBkgQcd","",1,20,140);
  TH1F* hQcd    = new TH1F("hQcd","",1,20,140);
  fullTreeBkgQcdCut->Draw("mass>>hBkgQcd","weight");
  fullTreeQcdCut->Draw("mass>>hQcd","weight");

  float fractionQCD = (float)hQcd->Integral()/(float)hBkgQcd->Integral();
  float expQCD = (float)hQcd->Integral()*(Lumi_/33.);

  cout<< "fraction " <<  fractionQCD << " of  " << hBkgQcd->Integral() << endl;

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

  float readEventsSgn = totalEventsSgnFake->GetBinContent(1);
  float expSgn = fullTreeSgnCut->GetEntries()*Lumi_/(readEventsSgn/(1300.*1.33));

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
  //RooRealVar* meanFit    = (RooRealVar*)(&FitParamSgn["meanSgn"]);
  //RooRealVar* widthFit   = (RooRealVar*)(&FitParamSgn["widthSgn"]);
  RooRealVar* m1SgnFit    = (RooRealVar*)(&FitParamSgn["m1Sgn"]);
  RooRealVar* sigmaSgnFit = (RooRealVar*)(&FitParamSgn["sigmaSgn"]);
  RooRealVar* alfaSgnFit  = (RooRealVar*)(&FitParamSgn["alfaSgn"]);
  RooRealVar* nSgnFit     = (RooRealVar*)(&FitParamSgn["nSgn"]);
  //RooRealVar* fSgnFit     = (RooRealVar*)(&FitParamSgn["fSgn"]);
  //RooRealVar* sigmaLSgnFit= (RooRealVar*)(&FitParamSgn["sigmaLSgn"]);

  RooRealVar m1Sgn_C("m1Sgn_C","m1",/*m1SgnFit->getVal()*/0,-10,10);
  RooRealVar sigmaSgn_C("sigmaSgn_C","sigma",/*sigmaSgnFit->getVal()*/0.5,0,20);
  RooRealVar alfaSgn_C("alfaSgn_C","alfa",alfaSgnFit->getVal()*(1+deltaAlpha_),0,20);
  RooRealVar nSgn_C("nSgn_C","n",nSgnFit->getVal()*(1+deltaN_),0,50);
  //RooConstVar fSgn_C("fSgn_C","n",fSgnFit->getVal());
  //RooConstVar sigmaLSgn_C("sigmaLSgn_C","n",sigmaLSgnFit->getVal());

  RooLognormal alfaSgn_CPdf("alfaSgn_CPdf","",alfaSgn_C,RooConst(alfaSgnFit->getVal()),RooConst(1.5));
  RooLognormal nSgn_CPdf("nSgn_CPdf","",nSgn_C,RooConst(nSgnFit->getVal()),RooConst(1.5));


  RooCBShape cbSgn_C("cbSgn_C","",mass,m1Sgn_C,sigmaSgn_C,alfaSgn_C,nSgn_C);
  //RooBifurGauss bifurcSgn_C("bifurcSgn_C","bifurc",mass,m1Sgn_C,sigmaLSgn_C,sigmaRSgn); 
  //RooAddPdf cbPlusBifurcSgn_C("cbPlusBifurcSgn_C","",RooArgList(cbSgn_C,bifurcSgn_C),fSgn_C);

  RooFFTConvPdf sgnPdf("sgnPdf","",mass,bwSgn,  /*cbPlusBifurcSgn_C*/ cbSgn_C);
  

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

  /*
  // gauss (X) landau  for Wen
  RooRealVar meanWen("meanWen","",  50,30,70);
  RooRealVar sigmaWen("sigmaWen","",20,5,50);
  RooGaussian gausWen("gausWen","",mass,meanWen,sigmaWen);
  RooRealVar meanLWen("meanLWen","",5,0,20);
  RooRealVar sigmaLWen("sigmaLWen","",10,0,50);
  RooLandau LWen("LWen","",mass,meanLWen,sigmaLWen);

  RooFFTConvPdf gausLWen("gausLWen","",mass,gausWen,LWen);

  RooFitResult* ResWenFit = gausLWen.fitTo(wenDataHist, Minos(1), Save(1), NumCPU(4) );
  RooArgSet FitParamWen(ResWenFit->floatParsFinal());
  RooRealVar* meanWenFit    = (RooRealVar*)(&FitParamWen["meanWen"]);
  RooRealVar* sigmaWenFit   = (RooRealVar*)(&FitParamWen["sigmaWen"]);
  RooRealVar* meanLWenFit   = (RooRealVar*)(&FitParamWen["meanLWen"]);
  RooRealVar* sigmaLWenFit  = (RooRealVar*)(&FitParamWen["sigmaLWen"]);
 
  RooConstVar meanWen_C("meanWen_C","",meanWenFit->getVal());
  RooConstVar sigmaWen_C("sigmaWen_C","",sigmaWenFit->getVal());
  RooConstVar meanLWen_C("meanLWen_C","",meanLWenFit->getVal());
  RooConstVar sigmaLWen_C("sigmaLWen_C","",sigmaLWenFit->getVal());

  RooGaussian gausWen_C("gausWen_C","",mass,meanWen_C,sigmaWen_C);
  RooGaussian LWen_C("LWen_C","",mass,meanLWen_C,sigmaLWen_C);

  RooFFTConvPdf wenPdf("wenPdf","",mass,gausWen_C,LWen_C);
  */

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

  /////////////////////////////////////////
  //      mix
  /////////////////////////////////////////
  
  mass.setBins( nBins_ );
  RooDataSet mixDataSet("sgnDataSet","dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeMixCut ) );
  if(!isData_){
    // this is not necessary, but good x-check
    mixDataSet.reset();
    mixDataSet.append( *((RooDataSet*)sgnDataSet.reduce(EventRange(1,(int)expSgn)))  );
    mixDataSet.append( *((RooDataSet*)wenDataSet.reduce(EventRange(1,(int)expWen)) ));
    mixDataSet.append( *((RooDataSet*)zttDataSet.reduce(EventRange(1 ,(int)expZtt)) ));
    mixDataSet.append( *((RooDataSet*)sgnFakeDataSet.reduce(EventRange(1 , (int)expSgnFake)))  );
    mixDataSet.append( *((RooDataSet*)ttbDataSet.reduce(EventRange(1 , (int)expTTb)))  );
    
    mixDataSet.append( *(qcdPdf.generate( mass, (int)(fractionQCD*mixDataSet.numEntries()) )) );
  }
  RooDataHist mixDataHist("sgnDataHist","",RooArgSet(mass),mixDataSet, 1.0);
  RooHistPdf  mixPdf("mixPdf","",RooArgSet(mass),mixDataHist);


  //////////////////////////////////////
  //      bkg + qcd
  //////////////////////////////////////

  mass.setBins( nBins_ );
  RooDataSet bkgDataSet("bkgDataSet","dataset for bkg-pass template", RooArgSet(mass), Import( *fullTreeBkgCut ) );
  bkgDataSet.append( *(qcdPdf.generate( mass, (int)(fractionQCD*mixDataSet.numEntries()) )) );
  int qcdYield = (int)(fractionQCD*mixDataSet.numEntries());
  RooDataHist bkgDataHist("bkgDataHist","",RooArgSet(mass),bkgDataSet, 1.0);

  // bkg pdf
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

  RooFormulaVar NzttN("NzttN","","Nztt/expZtt_cv",RooArgList(Nztt,expZtt_cv));
  RooConstVar NzttN_mean("NzttN_mean","",1.0);
  RooFormulaVar NqcdN("NqcdN","Nqcd/expQCD_cv",RooArgSet(Nqcd,expQCD_cv));
  RooConstVar NqcdN_mean("NqcdN_mean","",1.0);
  RooFormulaVar NwenN("Nwen","Nwen/expWen_cv",RooArgSet(Nwen,expWen_cv));
  RooConstVar NwenN_mean("NwenN_mean","",1.0);
  RooFormulaVar NttbN("NttbN","Nttb/expTTb_cv",RooArgSet(Nttb,expTTb_cv));
  RooConstVar NttbN_mean("NttbN_mean","",1.0);
  RooFormulaVar NsgnFakeN("NsgnFakeN","NsgnFake/expSgnFake_cv",RooArgSet(NsgnFake,expSgnFake_cv));
  RooConstVar NsgnFakeN_mean("NsgnFakeN_mean","",1.0);

  RooConstVar expQCD_err_cv("expQCD_err_cv","",         4.0 );    // 4.0 // 4.0
  RooConstVar expZtt_err_cv("expZtt_err_cv","",         0.4 );    // 0.4 // 0.4
  RooConstVar expWen_err_cv("expWen_err_cv","",         1.0 );    // 1.0 // 2.0
  RooConstVar expTTb_err_cv("expTTb_err_cv","",         2.0 );    // 2.0 // 1.0 
  RooConstVar expSgnFake_err_cv("expSgnFake_err_cv","", 1.0 );    // 1.0 // 1.0

  /*
  RooGaussian NqcdConstraint("NqcdConstraint","",NqcdN,NqcdN_mean,expQCD_err_cv) ;
  RooGaussian NzttConstraint("NzttConstraint","",NzttN, NzttN_mean , expZtt_err_cv) ;
  RooGaussian NwenConstraint("NwenConstraint","",NwenN,NwenN_mean,expWen_err_cv) ;
  RooGaussian NttbConstraint("NttbConstraint","",NttbN,NttbN_mean,expTTb_err_cv) ;
  RooGaussian NsgnFakeConstraint("NsgnFakeConstraint","",NsgnFakeN,NsgnFakeN_mean,expSgnFake_err_cv) ;
  */

  RooLognormal NqcdConstraint("NqcdConstraint","",Nqcd,expQCD_cv,RooConst(4)) ;
  RooLognormal NzttConstraint("NzttConstraint","",Nztt,expZtt_cv,RooConst(1.4)) ;
  RooLognormal NwenConstraint("NwenConstraint","",Nwen,expWen_cv,RooConst(2)) ;
  RooLognormal NttbConstraint("NttbConstraint","",Nttb,expTTb_cv,RooConst(2)) ;
  RooLognormal NsgnFakeConstraint("NsgnFakeConstraint","",NsgnFake,expSgnFake_cv,RooConst(2)) ;





  RooAddPdf bkgPdf("bkgPdf","",RooArgList(qcdPdf,zttPdf,wenPdf,ttbPdf,sgnFakePdf),RooArgList(Nqcd,Nztt,Nwen,Nttb,NsgnFake));

  /////////////////////////////////////////
  //     sum
  /////////////////////////////////////////

  RooRealVar Nsgn("Nsgn","",1000,0,1000000);
  RooRealVar Nbkg("Nbkg","",1000,0,1000000);
  RooRealVar purity("purity","",0.5,0,1);
  
  //RooAddPdf sum("sum","",RooArgList(sgnPdf,bkgPdf),RooArgList(Nsgn,Nbkg));
  RooAddPdf sum("sum","",RooArgList(sgnPdf,qcdPdf,zttPdf,wenPdf,ttbPdf,sgnFakePdf),RooArgList(Nsgn,Nqcd,Nztt,Nwen,Nttb,NsgnFake));
  //RooAddPdf sum("sum","",RooArgList(sgnPdf,bkgPdf),purity);
  
  RooAddPdf sumFail("sumFail","",RooArgList(sgnPdf,cbZtt),RooArgList(Nsgn,Nbkg));

  RooPlot* frameBkg     = mass.frame(  Bins(nBins_), Title("Sum of backgrounds (MC) + fit to data") );
  RooPlot* frameSgn     = mass.frame(  Bins(nBins_), Title("Signal template (MC)") );
  RooPlot* frameMix     = mass.frame(  Bins(nBins_), Title("Data set") );
  RooPlot* frameQcd     = mass.frame(  Bins(nBins_), Title("QCD template (data)") );
  RooPlot* frameWen     = mass.frame(  Bins(nBins_), Title("Wenu template (MC)") );
  RooPlot* frameZtt     = mass.frame(  Bins(nBins_), Title("Z#tau#tau template (MC)") );
  RooPlot* frameTTb     = mass.frame(  Bins(nBins_), Title("t#bar{t} template (MC)") );
  RooPlot* frameSgnFake = mass.frame(  Bins(nBins_), Title("Z->ee,jet->#tau  template (MC)") );

  mass.setBins( nBins_ );
  
  
  RooFitResult* fitRes = 0;
  if(doBinned_){
    fitRes = sum.fitTo( mixDataHist,ExternalConstraints( RooArgSet(NqcdConstraint,NzttConstraint,NwenConstraint,NttbConstraint,NsgnFakeConstraint,alfaSgn_CPdf,nSgn_CPdf) ),  Minos(1), Save(1), NumCPU(4) );
    if(fitInFail_) sumFail.fitTo( mixDataHist,  Minos(1), Save(1), NumCPU(4) );
  }
  else{
    fitRes = sum.fitTo( mixDataSet,ExternalConstraints( RooArgSet(NqcdConstraint,NzttConstraint,NwenConstraint,NttbConstraint,NsgnFakeConstraint) ),  Minos(1), Save(1), NumCPU(4)  );
    if(fitInFail_) sumFail.fitTo( mixDataSet,  Minos(1), Save(1), NumCPU(4)  );
  }
  //RooArgSet fitParam(fitRes->floatParsFinal());
  
  c1->cd();
  mixDataSet.plotOn(frameMix);
  sum.plotOn(frameMix,LineColor(kBlue));
  sum.plotOn(frameMix,LineColor(kRed),  Components("sgnPdf"));
  sum.plotOn(frameMix,LineColor(kGreen),Components("zttPdf"),LineStyle(kDashed));
  sum.plotOn(frameMix,LineColor(kBlack),Components("qcdPdf"),LineStyle(kDashed));
  sum.plotOn(frameMix,LineColor(kMagenta),Components("ttbPdf"),LineStyle(kDashed));
  sum.plotOn(frameMix,LineColor(kYellow),Components("wenPdf"),LineStyle(kDashed));
  sum.plotOn(frameMix,LineColor(kRed),Components("sgnFakePdf"),LineStyle(kDotted));
  frameMix->Draw();
  c1->Draw();

  /*
  c2->Divide(2,4);

  c2->cd(1);
  bkgDataSet.plotOn(frameBkg);
  bkgPdf.plotOn(frameBkg,LineColor(kBlue),LineStyle(kDashed));
  frameBkg->Draw();
  c2->cd(2);
  sgnDataSet.plotOn(frameSgn);
  sgnPdf.plotOn(frameSgn,LineColor(kRed),LineStyle(kSolid));
  frameSgn->Draw();
  c2->cd(3);
  mixDataSet.plotOn(frameMix);
  sum.plotOn(frameMix,LineColor(kBlue));
  sum.plotOn(frameMix,LineColor(kRed),  Components("sgnPdf"));
  sum.plotOn(frameMix,LineColor(kGreen),Components("zttPdf"),LineStyle(kDashed));
  sum.plotOn(frameMix,LineColor(kBlack),Components("qcdPdf"),LineStyle(kDashed));
  sum.plotOn(frameMix,LineColor(kMagenta),Components("ttbPdf"),LineStyle(kDashed));
  sum.plotOn(frameMix,LineColor(kYellow),Components("wenPdf"),LineStyle(kDashed));
  sum.plotOn(frameMix,LineColor(kRed),Components("sgnFakePdf"),LineStyle(kDotted));
  frameMix->Draw();
  c2->cd(4);
  qcdDataHist.plotOn(frameQcd);
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
  */
  
  c2->Divide(2,3);
  
  c2->cd(1);
  //bkgDataSet.plotOn(frameBkg);
  //bkgPdf.plotOn(frameBkg,LineColor(kBlue),LineStyle(kDashed));
  //frameBkg->Draw();
  sgnDataSet.plotOn(frameSgn);
  bvcbSgn.plotOn(frameSgn,LineColor(kRed),LineStyle(kSolid));
  frameSgn->Draw();
  c2->cd(2);
  sgnFakeDataSet.plotOn(frameSgnFake);
  sgnFakePdf.plotOn(frameSgnFake,LineColor(kRed),LineStyle(kDashed));
  frameSgnFake->Draw();
  c2->cd(3);
  qcdDataHist.plotOn(frameQcd);
  qcdPdf.plotOn(frameQcd,LineColor(kBlack),LineStyle(kDashed));
  frameQcd->Draw();
  c2->cd(4);
  wenDataSet.plotOn(frameWen);
  wenPdf.plotOn(frameWen,LineColor(kYellow),LineStyle(kDashed));
  frameWen->Draw();
  c2->cd(5);
  zttDataSet.plotOn(frameZtt);
  zttPdf.plotOn(frameZtt,LineColor(kGreen),LineStyle(kDashed));
  frameZtt->Draw();
  c2->cd(6);
  ttbDataSet.plotOn(frameTTb);
  ttbPdf.plotOn(frameTTb,LineColor(kMagenta),LineStyle(kDashed));
  frameTTb->Draw();
  
  c2->Update();
  c2->Draw();

  
  std::cout << "sgn " << expSgn  << ", qcd " << expQCD << ", Ztattau " <<expZtt << ", Wenu " << expWen << ", TTbar " << expTTb << ", sgn fake " << expSgnFake << endl; 

  std::cout << "QCD yield " << qcdYield << std::endl;

}




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
#include "TRandom3.h"

#include <vector>

using namespace std;
using namespace RooFit;


void toy( 
	 unsigned int nToys_    = 100,
	 double mean_H          = -1.0,        
	 double sigma_H         = 2.0,        
	 const string tnp_      = "etoTauMargLooseNoCracks70",
	 const string category_ = "tauAntiEMVA",
	 double cutValue_       = 0.5,
	 const string bin_      = "abseta<1.5",
	 double nBins_          = 25,
	 double xLow_           = 65,
	 double xHigh_          = 120
	 )
{

  TCanvas *c2 = new TCanvas("fitCanvasTemplate","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  TH1F* histoTempl   = new TH1F("histoTempl",  "Resolution model is null",20,-4,4);
  TH1F* histoResol   = new TH1F("histoResol",  "Resolution model left floating",20,-4,4);
  TH1F* histoResol_H = new TH1F("histoResol_H",Form("Resolution mis-modelled by (%.1f,%.1f)",mean_H,sigma_H),20,-4,4);

  TFile fsup("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_soup_tauAntiEMVA.root");
  //TFile fsup("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");
  TTree *fullTreeSoup = (TTree*)fsup.Get((tnp_+"/fitter_tree").c_str());

  TFile fsgnHighStat("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_DYToEE-PYTHIA.root");
  TTree *fullTreeSgnHighStat = (TTree*)fsgnHighStat.Get((tnp_+"/fitter_tree").c_str());

  int all = (int)fullTreeSgnHighStat->GetEntries(Form("(mcTrue && %s && mass>%f && mass<%f)",bin_.c_str(),xLow_,xHigh_));
  int passing = (int)fullTreeSgnHighStat->GetEntries(Form("(mcTrue && %s>=%f && %s && mass>%f && mass<%f)",category_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_));
  float efficiency      = (float)passing/(float)all;
  float efficiencyError = TMath::Sqrt(efficiency*(1-efficiency)/(float)all);

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
  mass.setBins(10000,"fft");
  mass.setBins( nBins_ );

  RooDataSet templateP("templateP","dataset for signal-pass template", RooArgSet(mass,mcTrue,tauAntiEMVA,signalPFChargedHadrCands,abseta,matchedID,leadPFChargedHadrCandTrackPt,pt), Import( *fullTreeSgnHighStat ), Cut(Form("(mcTrue && %s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) ) );
  RooDataHist templateHistP("templateHistP","",RooArgSet(mass), templateP, 1.0);
  RooHistPdf TemplateSignalPdfP("TemplateSignalPdfP","",RooArgSet(mass),templateHistP);

  RooRealVar McMeanResP("McMeanResP","",0,-10,10);
  RooRealVar McSigmaResP("McSigmaResP","",0.5,0,10);
  RooGaussian McResolModP("McResolModP","",mass,McMeanResP,McSigmaResP);
  RooFFTConvPdf McSignalPdfP("McSignalPdfP","",mass,TemplateSignalPdfP,McResolModP);

  RooRealVar McMeanResP_H("McMeanResP_H","",  mean_H);
  RooRealVar McSigmaResP_H("McSigmaResP_H","",sigma_H);
  RooGaussian McResolModP_H("McResolModP_H","",mass,McMeanResP_H,McSigmaResP_H);
  RooFFTConvPdf McSignalPdfP_H("McSignalPdfP_H","",mass,TemplateSignalPdfP,McResolModP_H);

  RooRealVar McCP("McCP","",0,-10,10);
  RooExponential McBackgroundPdfP("McBackgroundPdfP","",mass,McCP);
  
  RooRealVar McNumBkgP("McNumBkgP","",0,1000000);
  RooRealVar McNumSgnP("McNumSgnP","",0,1000000);
  
  RooAddPdf McModelPResol("McModelPResol","",RooArgList(McBackgroundPdfP,McSignalPdfP),RooArgList(McNumBkgP,McNumSgnP));
  RooAddPdf McModelPResol_H("McModelPResol_H","",RooArgList(McBackgroundPdfP,McSignalPdfP_H),RooArgList(McNumBkgP,McNumSgnP));
  RooAddPdf McModelPTempl("McModelPTempl","",RooArgList(McBackgroundPdfP,TemplateSignalPdfP),RooArgList(McNumBkgP,McNumSgnP));


  TFile *McP = new TFile("dummy1.root","RECREATE");
  TTree* treeSoupP = fullTreeSoup->CopyTree( Form("(%s>=%f && %s)",category_.c_str(),cutValue_,bin_.c_str()) );
  int SOUPPass = (int)treeSoupP->GetEntries(Form("mass>%f && mass<%f",xLow_,xHigh_));
  int SOUPFail = (int)fullTreeSoup->GetEntries(Form("(%s<%f && %s && mass>%f && mass<%f)",category_.c_str(),cutValue_,bin_.c_str(),xLow_,xHigh_));

  
  RooDataSet McDataP("McDataP","dataset pass for the soup pass", RooArgSet(mass), Import( *treeSoupP ) );
  RooDataHist McDataHistP("McDataHistP","",RooArgSet(mass), McDataP, 1.0);
  RooHistPdf McDataPdfP("McDataPdfP","",RooArgSet(mass),McDataHistP);

  int eventsPass   = (int)SOUPPass*(33./500.);
  int eventsFail   = (int)SOUPFail*(33./500.);

  TRandom3* random = new TRandom3();

  /*
  for(unsigned int toy = 0; toy< 100; toy++){
    
    RooDataSet* toySet  = McDataPdfP.generate(mass, eventsPass, Extended() );
    RooDataHist toyHist("toyHist","",RooArgSet(mass), *toySet, 1.0);
    
    // resolution is null
    RooFitResult* ResMcCombinedTemplFit = McModelPTempl.fitTo(toyHist , Extended(1), Minos(1), Save(1), NumCPU(4) );
    RooArgSet McFitParamTempl(ResMcCombinedTemplFit->floatParsFinal());
    RooRealVar* McNumSigFitTempl  = (RooRealVar*)(&McFitParamTempl["McNumSgnP"]);
    RooRealVar* McNumBkgPFitTempl = (RooRealVar*)(&McFitParamTempl["McNumBkgP"]);

    histoTempl->Fill(McNumSigFitTempl->getVal()/(toySet->sumEntries()));
 }

 histoTempl->Sumw2();
 TF1* f1 = new TF1("f1","gaus",-1,1);
 histoTempl->Fit(f1);
 float meanPurity =  f1->GetParameter(1);
  */

  RooPlot* McFrameP = mass.frame(Bins(nBins_),Title("MC: passing sample"));
  for(unsigned int toy = 0; toy< nToys_; toy++){

    //failing sample:
    int failing = random->Poisson(eventsFail);

    float i_efficiency = random->Gaus(efficiency,efficiencyError);

    RooDataSet* toySet = McDataPdfP.generate(mass, eventsPass, Extended() );
    RooDataHist toyHist("toyHist","",RooArgSet(mass), *toySet, 1.0);

    RooFitResult* ResMcCombinedTemplFit = McModelPTempl.fitTo(toyHist , Extended(1), Minos(1), Save(1), NumCPU(4) );
    RooArgSet McFitParamTempl(ResMcCombinedTemplFit->floatParsFinal());
    RooRealVar* McNumSigFitTempl  = (RooRealVar*)(&McFitParamTempl["McNumSgnP"]);
    RooRealVar* McNumBkgPFitTempl = (RooRealVar*)(&McFitParamTempl["McNumBkgP"]);
    float fittedEff = McNumSigFitTempl->getVal()/(failing+McNumSigFitTempl->getVal());
    histoTempl->Fill( (fittedEff-i_efficiency)/McNumSigFitTempl->getError()*(failing+McNumSigFitTempl->getVal()) );
   
    // resolution left floating
    RooFitResult* ResMcCombinedResolFit = McModelPResol.fitTo(toyHist , Extended(1), Minos(1), Save(1), NumCPU(4) );
    RooArgSet McFitParamResol(ResMcCombinedResolFit->floatParsFinal());
    RooRealVar* McNumSigFitResol  = (RooRealVar*)(&McFitParamResol["McNumSgnP"]);
    RooRealVar* McNumBkgPFitResol = (RooRealVar*)(&McFitParamResol["McNumBkgP"]);
    fittedEff = McNumSigFitResol->getVal()/(failing+McNumSigFitResol->getVal());
    histoResol->Fill( (fittedEff-i_efficiency)/McNumSigFitResol->getError()*(failing+McNumSigFitResol->getVal()) );

    // resolution fixed to some wrong value H
    RooFitResult* ResMcCombinedResolFit_H = McModelPResol_H.fitTo(toyHist , Extended(1), Minos(1), Save(1), NumCPU(4) );
    RooArgSet McFitParamResol_H(ResMcCombinedResolFit_H->floatParsFinal());
    RooRealVar* McNumSigFitResol_H  = (RooRealVar*)(&McFitParamResol_H["McNumSgnP"]);
    RooRealVar* McNumBkgPFitResol_H = (RooRealVar*)(&McFitParamResol_H["McNumBkgP"]);
    fittedEff = McNumSigFitResol_H->getVal()/(failing+McNumSigFitResol_H->getVal());
    histoResol_H->Fill( (fittedEff-i_efficiency)/McNumSigFitResol_H->getError()*(failing+McNumSigFitResol_H->getVal()) );

    
    /*   
    toyHist.plotOn(McFrameP,Bins(nBins_),Title("MC: passing sample"));
    McModelPResol.plotOn(McFrameP, LineColor(kBlue));
    McModelPResol.plotOn(McFrameP, LineColor(kGreen),Components("McBackgroundPdfP"));
    McModelPResol.plotOn(McFrameP, LineColor(kRed),Components("McSignalPdfP"));
    McModelPTempl.plotOn(McFrameP, LineColor(kRed),LineStyle(kDashed),Components("TemplateSignalPdfP"));
    cout << "########################################" << endl;
    break;
    */

  }

  //McFrameP->Draw();

  histoTempl->Sumw2();
  histoResol->Sumw2();
  histoResol_H->Sumw2();

  TF1* gauss = new TF1("gauss","gaus",-2,2);

  c2->Divide(2,2);
  c2->cd(1);
  histoResol->Fit(gauss);
  histoResol->Draw("P");
  c2->cd(2);
  histoResol_H->Fit(gauss);
  histoResol_H->Draw("P");
  c2->cd(3);
  histoTempl->Fit(gauss);
  histoTempl->Draw("P");
  //c2->cd(4);
  //histoTempl->Fit(f1);
  //histoTempl->Draw("HIST");


  c2->cd();
  c2->Draw();

  delete random;

}

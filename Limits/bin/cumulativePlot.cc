#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TKey.h"
#include "TObjString.h"
#include "TObject.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TList.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCut.h"
#include "TArrayF.h"
#include "TStyle.h"
#include "TCollection.h"
#include "TGraphAsymmErrors.h"
#include "TGraphPainter.h"
#include "TMultiGraph.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>

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
#include "RooLognormal.h"


using namespace RooFit;
using namespace std;



void mergeDifferentHistos(TH1F* hTarget=0, TH1F* hInput=0, float weight = 1.0, bool isError = false, bool verbose = false){

  
  for(unsigned int i = 1; i <= (unsigned int)(hTarget->GetNbinsX()); i++){
 
    float centralValueInput = hTarget->GetBinCenter(i);
    float widthInput        = hTarget->GetXaxis()->GetBinUpEdge(i)-hTarget->GetXaxis()->GetBinLowEdge(i);
    float minusHalfInput    = (centralValueInput-widthInput/2.) ;//+ 0.001;
    float plusHalfInput     = (centralValueInput+widthInput/2.) ;//- 0.001;

    if(verbose){
      cout << "Target Bin " << i << endl;
      cout << " ==> centralValueInput = " << centralValueInput << endl;
      cout << " ==> widthInput = " << widthInput << endl; 
      cout << " ==> minusHalfInput = " << minusHalfInput << "; plusHalfInput = " << plusHalfInput << endl; 
    }

    int binLowInput = 0; float widthBinLowInput =0. ; float binLowInputHedgeHigh =0. ; float binLowInputFracHigh=0.;
    if(minusHalfInput<=hInput->GetXaxis()->GetXmax()){
      binLowInput          = hInput->FindBin(minusHalfInput);
      widthBinLowInput     = hInput->GetXaxis()->GetBinUpEdge(binLowInput)-hInput->GetXaxis()->GetBinLowEdge(binLowInput);
      binLowInputHedgeHigh = hInput->GetXaxis()->GetBinUpEdge(binLowInput);
      binLowInputFracHigh  = (binLowInputHedgeHigh-minusHalfInput)/widthBinLowInput;
      if(binLowInputFracHigh<0 || binLowInputFracHigh>1){
	cout << "Low bin fraction is not ok... " << binLowInputFracHigh << endl; 
      }
      else
	if(verbose){
	  cout << "binLowInput = " << binLowInput << ", widthBinLowInput = " << widthBinLowInput << ", binLowInputHedgeHigh = " << binLowInputHedgeHigh <<  endl;
	}
    }
    else
      binLowInput = hInput->GetNbinsX()+1;

    int binHighInput  = 0; float widthBinHighInput=0.; float binHighInputHedgeLow=0.; float binHighInputFracLow=0.;
    if(plusHalfInput<=hInput->GetXaxis()->GetXmax()){
      binHighInput         = hInput->FindBin(plusHalfInput);
      widthBinHighInput    = hInput->GetXaxis()->GetBinUpEdge(binHighInput)-hInput->GetXaxis()->GetBinLowEdge(binHighInput);
      binHighInputHedgeLow = hInput->GetXaxis()->GetBinLowEdge(binHighInput);
      binHighInputFracLow  = (plusHalfInput-binHighInputHedgeLow)/widthBinHighInput;
      if(binHighInputFracLow<0 || binHighInputFracLow>1){
	cout << "Low High fraction is not ok... " << binHighInputFracLow << endl; 
      }
      else
	if(verbose){
	  cout << "binHighInput = " << binHighInput << ", widthBinHighInput = " << widthBinHighInput << ", binHighInputHedgeLow = " << binHighInputHedgeLow <<  endl;
	}
    }
    else
      binHighInput = hInput->GetNbinsX();
      
    float entryi = 0.; float errori2 = 0.;

    for(int j = binLowInput; j <= binHighInput; j++){
      float widthj = hInput->GetBinWidth(j);
      float inputj = hInput->GetBinContent(j);

      float errorj = hInput->GetBinError(j)*widthj;
      float entryj = inputj*widthj;
      if(j==binLowInput){
	entryi  += binLowInputFracHigh*entryj;
	errori2 += (binLowInputFracHigh*errorj)*(binLowInputFracHigh*errorj);
	if(verbose) cout << entryi << endl;
      }
      else if(j==binHighInput){
	entryi  += binHighInputFracLow*entryj;
	errori2 += (binHighInputFracLow*errorj)*(binHighInputFracLow*errorj);
	if(verbose) cout << "+" << entryi << endl;
      }
      else{
	entryi  += entryj;
	errori2 += (errorj*errorj);
	if(verbose) cout << "+" << entryj << endl;
      }
    }

    float oldContenti = hTarget->GetBinContent(i);
    float oldErrori   = hTarget->GetBinError(i);
    oldContenti *= widthInput;
    oldErrori   *= widthInput;
    float newContenti = oldContenti + entryi*weight;
    float newErrori   = TMath::Sqrt(oldErrori*oldErrori + errori2*weight*weight);
    hTarget->SetBinContent(i, newContenti/widthInput);
    hTarget->SetBinError(i,   newErrori/widthInput);

  }

  if(verbose)  cout << "#### Initial ==> " << hInput->Integral("width") << endl;
  if(verbose)  cout << "#### Final   ==> " << hTarget->Integral("width") << endl;

}



void produce(int mH                   = 120,
	     float signalRescale      = 0.20, 
	     vector<string>& channels = *(new  vector<string>()),
	     vector<string>& com      = *(new  vector<string>()),
	     Int_t nBins_ = 80, 
	     Float_t xMin_=  0, Float_t xMax_=400,
	     TString title            = ""
	     ){

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(0);

  TPad* pad1 = new TPad("pad1","",0.05,0.22,0.96,0.97);
  TPad* pad2 = new TPad("pad2","",0.05,0.02,0.96,0.20);
  pad1->SetFillColor(0);
  pad2->SetFillColor(0);
  pad1->Draw();
  pad2->Draw();

  TCanvas *c1W = new TCanvas("c1W","",5,30,650,600);
  c1W->SetGrid(0,0);
  c1W->SetFillStyle(4000);
  c1W->SetFillColor(10);
  c1W->SetTicky();
  c1W->SetObjectStat(0);
  c1W->SetLogy(0);

  TPad* pad1W = new TPad("pad1W","",0.05,0.22,0.96,0.97);
  TPad* pad2W = new TPad("pad2W","",0.05,0.02,0.96,0.20);
  pad1W->SetFillColor(0);
  pad2W->SetFillColor(0);
  pad1W->Draw();
  pad2W->Draw();

  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.1,"y");

  TLegend* leg = new TLegend(0.50,0.65,0.75,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);
  //leg->SetHeader(title);

  for(unsigned int en = 0; en<com.size(); en++){
    cout << com[en] << endl;
  }

  typedef std::map<string, std::pair<vector<string>,vector<string> > >::iterator IT;
  std::map<string, std::pair<vector<string>,vector<string> > > channelToCatMap; 

  vector<string> categoriesLT; vector<string> categoriesIndexLT;
  vector<string> categoriesTT; vector<string> categoriesIndexTT;
  vector<string> categoriesMM; vector<string> categoriesIndexMM;

  //categoriesLT.push_back("0jet_low");   categoriesIndexLT.push_back("0");
  //categoriesLT.push_back("0jet_high");  categoriesIndexLT.push_back("1");
  //categoriesLT.push_back("boost_low");  categoriesIndexLT.push_back("2");
  //categoriesLT.push_back("boost_high"); categoriesIndexLT.push_back("3");
  categoriesLT.push_back("vbf");        categoriesIndexLT.push_back("5");

  categoriesTT.push_back("boost");      categoriesIndexTT.push_back("0");
  categoriesTT.push_back("vbf");        categoriesIndexTT.push_back("1");

  categoriesMM.push_back("dummy");      categoriesIndexTT.push_back("0");

  vector<string> channelsDir;
  for(unsigned int ch = 0; ch<channels.size(); ch++){
    cout << channels[ch] << endl;
    if(channels[ch].find("mt")!=string::npos){
      channelToCatMap.insert( make_pair("mt", make_pair(categoriesLT, categoriesIndexLT) )  );
      channelsDir.push_back("muTau");
    }
    if(channels[ch].find("et")!=string::npos){
      channelToCatMap.insert( make_pair("et", make_pair(categoriesLT, categoriesIndexLT) )  );
      channelsDir.push_back("eleTau");
    }
    if(channels[ch].find("tt")!=string::npos){
      channelToCatMap.insert( make_pair("tt", make_pair(categoriesTT, categoriesIndexTT) )  );
      channelsDir.push_back("tauTau");
    }
    if(channels[ch].find("em")!=string::npos){
      channelToCatMap.insert( make_pair("em", make_pair(categoriesLT, categoriesIndexLT) )  );
      channelsDir.push_back("emu");
    }
    if(channels[ch].find("mm")!=string::npos){
      channelToCatMap.insert( make_pair("mm", make_pair(categoriesMM, categoriesIndexMM) )  );
      channelsDir.push_back("muMu");
    }
  }


  // input txt file with bins
  ifstream is;

  char* c = new char[10];
  is.open("../test/bins_coarse.txt"); 
  if(nBins_<0 &&  !is.good()){
    cout << "Bins file not found" << endl;
    return;
  }

  int nBinsFromFile = 0;
  while (is.good())     
    {
      is.getline(c,999,',');     
      if (is.good()){
	nBinsFromFile++;
      }
    }

  int nBins =  nBins_>0 ? nBins_ : nBinsFromFile-1 ;
  TArrayF bins(nBins+1);
  cout << "Making histograms with " << nBins << " bins:" << endl;

  is.close();
  is.open("../test/bins_coarse.txt");
  
  nBinsFromFile = 0;

  if(nBins_>0){
    for( ; nBinsFromFile <= nBins ; nBinsFromFile++){
      bins[nBinsFromFile] =  xMin_ + nBinsFromFile*(xMax_-xMin_)/nBins_;
    }
  }
  else{
    while (is.good())  
      {
	is.getline(c,999,',');     
	if (is.good() && nBinsFromFile<=nBins) {
	  bins[nBinsFromFile] = atof(c);
	  cout << bins[nBinsFromFile] << ", " ;
	}
	nBinsFromFile++;
      }
    cout << endl;
  }

  THStack* aStack  = new THStack("aStack","");
  THStack* aStackW = new THStack("aStack","");

  TH1F* hData  = new TH1F( "hData" ,"Data; mass (GeV); events/GeV"            , nBins , bins.GetArray());
  TH1F* hDataW = new TH1F( "hDataW","Data weighted; mass (GeV); weighted events/GeV"   , nBins , bins.GetArray());
  TH1F* hBkg   = new TH1F( "hBkg"   ,"MC"              , nBins , bins.GetArray());
  TH1F* hBkgW  = new TH1F( "hBkgW"  ,"MC weighted"     , nBins , bins.GetArray());
  TH1F* hSgn   = new TH1F( "hSgn"  ,"Signal"          , nBins , bins.GetArray());
  TH1F* hSgnW  = new TH1F( "hSgnW" ,"Signal weightes" , nBins , bins.GetArray());
  TH1F* hErr   = new TH1F( "hErr"  ,"Signal"          , nBins , bins.GetArray());
  TH1F* hErrW  = new TH1F( "hErrW" ,"Signal weightes" , nBins , bins.GetArray());

  hBkg->SetFillColor(kBlue);
  hBkg->SetLineColor(kBlue);
  hBkg->SetLineWidth(1);
  hBkg->SetFillStyle(3003);
  hBkgW->SetFillColor(kBlue);
  hBkgW->SetFillColor(kBlue);
  hBkg->SetLineWidth(1);
  hBkgW->SetFillStyle(3003);

  hData->Sumw2();
  hData->SetLabelSize(0.04,"X");
  hData->SetTitleSize(0.06,"X");
  hData->SetLabelSize(0.04,"Y");
  hData->SetTitleSize(0.06,"Y");
  hData->SetMarkerStyle(20);
  hData->SetMarkerSize(1.0);
  hData->SetMarkerColor(kBlack);
  hData->SetLineColor(kBlack);
  hDataW->Sumw2();
  hData->SetLabelSize(0.04,"X");
  hData->SetTitleSize(0.06,"X");
  hData->SetLabelSize(0.04,"Y");
  hData->SetTitleSize(0.06,"Y");
  hDataW->SetMarkerStyle(20);
  hDataW->SetMarkerSize(1.0);
  hDataW->SetMarkerColor(kBlack);
  hDataW->SetLineColor(kBlack);

  hSgn->SetFillColor(kRed);
  hSgn->SetLineColor(kRed);
  hSgn->SetFillStyle(1001);
  hSgn->SetLineWidth(1);
  hSgn->SetLineStyle(kSolid);
  hSgnW->SetFillColor(kRed);
  hSgnW->SetFillStyle(1001);
  hSgnW->SetLineColor(kRed);
  hSgnW->SetLineWidth(1);
  hSgnW->SetLineStyle(kSolid);

  hErr->SetMarkerSize(0);
  hErr->SetFillColor(1);
  hErr->SetFillStyle(3013);
  hErrW->SetMarkerSize(0);
  hErrW->SetFillColor(1);
  hErrW->SetFillStyle(3013);

  leg->AddEntry(hData,"Observed","P");
  leg->AddEntry(hBkg,"Background","F");
  leg->AddEntry(hSgn,Form("SM Higgs (m_{H}=%d GeV)",mH),"F");

  TList *listData = new TList();
  TList *listMC   = new TList();
  TList *listSgn  = new TList();
  TList *listDataW = new TList();
  TList *listMCW   = new TList();
  TList *listSgnW  = new TList();

  float totalSgn = 0.;
  
  for(unsigned int en = 0; en<com.size(); en++){
    string iEn = com[en];

    for(unsigned int ch = 0; ch<channels.size(); ch++){

      string iCh   = channels[ch];
      string iChIn = channelsDir[ch];

      vector<string> categories      = channelToCatMap[iCh].first;
      vector<string> categoriesIndex = channelToCatMap[iCh].second;
      
      for(unsigned int cat = 0; cat< categories.size(); cat++){

	string iCat   = categories[cat];
	string iCatIn = categoriesIndex[cat];

	cout << "Now doing " << string(Form("%s_%s_rescaled_%s_.root", iChIn.c_str(),iCat.c_str(),iEn.c_str() )) << endl;
	
	TFile* fInput = new TFile(Form("../test/%s_%s_rescaled_%s_.root", iChIn.c_str(),iCat.c_str(),iEn.c_str() ),"READ");
	if(!fInput || fInput->IsZombie() ) continue;
	
	TH1F* ihSgn = 0;  TH1F* ihBkg  = 0; TH1F* ihData  = 0; TH1F* ihErr  = 0;
	TH1F* ihSgnW = 0; TH1F* ihBkgW = 0; TH1F* ihDataW = 0; TH1F* ihErrW = 0;

	TIter nextkey(gDirectory->GetListOfKeys());
	TH1F *key;
	while ( (key = (TH1F*)nextkey()) ) {
	  string name(key->GetName());

	  TH1F* h = (TH1F*)gDirectory->FindObjectAny( key->GetName());

	  if(  name.find("ggH")!=string::npos ){
	    //cout << "This is signal ... " << name << endl;
	    if( !ihSgn ){
	      ihSgn  = (TH1F*)h->Clone("ihSgn");
	      ihSgnW = (TH1F*)h->Clone("ihSgnW");
	      ihSgn->Add(h);
	      ihSgnW->Add(h);
	    }
	    else{
	      ihSgn->Add(h);
	      ihSgnW->Add(h);
	    }
	  }
	  else if(  name.find("data_obs")==string::npos &&  name.find("Ztt")!=string::npos) {
	    //cout << "This is background ... " << name << endl;
	    if( !ihBkg ){
	      ihBkg  = (TH1F*)h->Clone("ihBkg");
	      ihBkgW = (TH1F*)h->Clone("ihBkgW");
	      ihBkg->Add(h);
	      ihBkgW->Add(h);
	    }
	    else{
	      ihBkg->Add(h);
	      ihBkgW->Add(h);
	    }
	  }
	  else if( name.find("data_obs")!=string::npos){
	    //cout << "This is data ... " << name << endl;
	    if( !ihData ){
	      ihData  = (TH1F*)h->Clone("ihData");
	      ihDataW = (TH1F*)h->Clone("ihDataW");
	      ihData->Add(h);
	      ihDataW->Add(h);
	    }
	    else{
	      ihData->Add(h);
	      ihDataW->Add(h);
	    }
	  }
	  else if( name.find("errorBand")!=string::npos ){
	    //cout << "This is error band ... " << name << endl;
	    if( !ihErr ){
	      ihErr  = (TH1F*)h->Clone("ihErr");
	      ihErrW = (TH1F*)h->Clone("ihErrW");
	      ihErr->Add(h);
	      ihErrW->Add(h);
	    }
	    else{
	      ihErr->Add(h);
	      ihErrW->Add(h);
	    }
	  }
	  else{
	    //cout << "This histo will not be considered!" << endl;
	  }
	}

	if( !ihSgn || !ihBkg || !ihData){
	  cout << "Could not find some histos... continue" << endl;
	  continue;
	}

	if(ihSgn->Integral() > ihBkg->Integral()){
	  ihSgn->Add(ihBkg, -1.0);
	  ihSgnW->Add(ihBkg, -1.0);
	}
	ihSgn->Scale(signalRescale);
	ihSgnW->Scale(signalRescale);

	totalSgn += ihSgn->Integral("width");

	//cout << "Sgn  = " << ihSgn->Integral() << endl;
	//cout << "Bkg  = " << ihBkg->Integral() << endl;
	//cout << "Data = " << ihData->Integral() << endl;

	double quantiles[] = {0.05,0.16, 0.50, 0.84, 0.95};
	double positions[] = {0.0,0.0, 0.0, 0.0, 0.0};
	ihSgn->GetQuantiles(5, positions, quantiles);
	int binLow  = ihSgn->FindBin(positions[1]);
	int binHigh = ihSgn->FindBin(positions[3]);

	float totalSig = ihSgn->Integral(binLow,binHigh);
	float totalBkg = ihBkg->Integral(binLow,binHigh);

	float SoB     = totalSig/(totalBkg+totalSig);
	float SoSqrtB = totalSig/TMath::Sqrt(totalBkg);
	float SoBpS   = totalSig/(totalBkg+totalSig);

	cout << "Weight S/B = "       << SoB << endl;
	//cout << "Weight S/sqrt(B) = " << SoSqrtB << endl;
	//cout << "Weight S/(S+B) = "   << SoBpS << endl;

	//ihDataW->Scale(SoB);
	//ihSgnW->Scale(SoB);
	//ihBkgW->Scale(SoB);
	//ihErrW->Scale(SoB);

	mergeDifferentHistos(hSgn,  ihSgn,   1.0);
	mergeDifferentHistos(hSgnW, ihSgnW , SoB);
	mergeDifferentHistos(hBkg,  ihBkg,   1.0);
	mergeDifferentHistos(hBkgW, ihBkgW , SoB);
	mergeDifferentHistos(hData, ihData,  1.0);
	mergeDifferentHistos(hDataW,ihDataW, SoB);
	mergeDifferentHistos(hErr, ihErr,    1.0, true);
	mergeDifferentHistos(hErrW,ihErrW,   SoB, true);

	listData->Add(ihData);
	listDataW->Add(ihDataW);
	listMC->Add(ihBkg);
	listMCW->Add(ihBkgW);
	listSgn->Add(ihSgn);
	listSgnW->Add(ihSgnW);

	cout << "SgnW ==> " << hSgnW->Integral() << endl;

      }
      
    }
     
  }
    

//   hSgn->Merge(listSgn);
//   hSgnW->Merge(listSgnW);
//   hData->Merge(listData);
//   hDataW->Merge(listDataW);
//   hBkg->Merge(listMC);
//   hBkgW->Merge(listMCW);
//   hSgnW->Scale(totalSgn/hSgnW->Integral("width"));
//   hDataW->Scale(totalSgn/hSgnW->Integral("width"));
//   hBkgW->Scale(totalSgn/hSgnW->Integral("width"));

  TH1F* hErrSub = (TH1F*)hErr->Clone("hErrSub");
  hErrSub->Add(hBkg,-1.0);

  TH1F* hErrSubW = (TH1F*)hErrW->Clone("hErrSubW");
  hErrSubW->Add(hBkgW,-1.0);

  aStack->Add(hBkg);
  aStack->Add(hSgn);

  aStackW->Add(hBkgW);
  aStackW->Add(hSgnW);

  cout << "Unweighted = " << endl;
  cout << " Sgn TOT = " << hSgn->Integral() << endl;
  cout << " Bkg TOT = " << hBkg->Integral() << endl;
  cout << " Data TOT = " << hData->Integral() << endl;
  cout << "Weighted = " << endl;
  cout << " SgnW TOT = " << hSgnW->Integral() << endl;
  cout << " BkgW TOT = " << hBkgW->Integral() << endl;
  cout << " DataW TOT = " << hDataW->Integral() << endl;

  c1->cd();
  pad1->cd();
  hData->Draw("P");
  aStack->Draw("HISTSAME");
  hData->Draw("PSAME");

  leg->Draw();

  pad2->cd();

  TH1F* hRatio = (TH1F*)hData->Clone("hRatio");
  hRatio->SetTitle("");
  hRatio->SetYTitle("");
  hRatio->SetXTitle("");
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerSize(1.0);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetLineColor(kBlack);
  hRatio->SetLabelSize(0.12,"Y");
  hRatio->SetLabelSize(0.12,"X");
  hRatio->Add(hBkg,-1.0);
  for(unsigned int i = 1; i <= (unsigned int)(hRatio->GetNbinsX()); i++)
    hRatio->SetBinError(i, hData->GetBinError(i));
    

  TF1* line = new TF1("line","0",hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());
  line->SetLineStyle(3);
  line->SetLineWidth(1.5);
  line->SetLineColor(kBlack);

  hRatio->Draw("PE");
  hSgn->Draw("HISTSAME");
  hRatio->Draw("PESAME");
  hErrSub->Draw("e2SAME");
  line->Draw("SAME");

  c1->SaveAs("unweighted.png");


  c1W->cd();

  pad1W->cd();
  hDataW->Draw("P");
  aStackW->Draw("HISTSAME");
  hDataW->Draw("PSAME");
  hErrW->Draw("e2SAME");
  leg->Draw();

  pad2W->cd();

  TH1F* hRatioW = (TH1F*)hDataW->Clone("hRatioW");
  hRatioW->SetTitle("");
  hRatioW->SetYTitle("");
  hRatioW->SetXTitle("");
  hRatioW->SetMarkerStyle(20);
  hRatioW->SetMarkerSize(1.0);
  hRatioW->SetMarkerColor(kBlack);
  hRatioW->SetLineColor(kBlack);
  hRatioW->SetLabelSize(0.10,"Y");
  hRatioW->SetLabelSize(0.10,"X");
  hRatioW->Add(hBkgW,-1.0);
  for(unsigned int i = 1; i <= (unsigned int)(hRatioW->GetNbinsX()); i++)
    hRatioW->SetBinError(i, hDataW->GetBinError(i));


  hRatioW->Draw("PE");
  hSgnW->Draw("HISTSAME");
  hRatioW->Draw("PESAME");
  hErrSubW->Draw("e2SAME");
  line->Draw("SAME");

  c1W->SaveAs("weighted.png");

  delete hData; delete hDataW; delete hBkg; delete hBkgW; delete hSgn; delete hSgnW;

  return;
}


  

void produceAll(){

  vector<string> sevenTeV;
  sevenTeV.push_back("7TeV");

  vector<string> eightTeV;
  eightTeV.push_back("8TeV");
  
  vector<string> allTeV;
  allTeV.push_back("7TeV"); allTeV.push_back("8TeV");

  vector<string> mt;
  mt.push_back("mt");

  vector<string> et;
  et.push_back("et");

  vector<string> tt;
  tt.push_back("tt");

  vector<string> em;
  em.push_back("em");

  vector<string> lt;
  lt.push_back("mt"); lt.push_back("et");

  vector<string> lx;
  lx.push_back("em"); 
  lx.push_back("mt"); /*lx.push_back("et");*/  

  vector<string> xt;
  xt.push_back("mt"); xt.push_back("et");  xt.push_back("tt");

  vector<string> xx;
  xx.push_back("mt"); xx.push_back("et");  xx.push_back("tt");  xx.push_back("em"); 

  produce(125, 0.20,lx, allTeV, -1, 0, 400, "");

  return;

}



int main(int argc, const char* argv[])
{

  std::cout << "produce()" << std::endl;
  gROOT->SetBatch(true);
 

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  produceAll();

}



// 	for(unsigned int i = 1; i <= (unsigned int)(hData->GetNbinsX()); i++){
// 	  float centralValue = hData->GetBinCenter(i);
// 	  //cout <<  " ==> " << centralValue << endl;

// 	  float totData = ihData->GetBinContent(ihData->FindBin(centralValue));
// 	  float totSgn  = ihSgn->GetBinContent(ihSgn->FindBin(centralValue));
// 	  float totBkg  = ihBkg->GetBinContent(ihBkg->FindBin(centralValue));

// 	  float binWidth= ihData->GetBinWidth(ihData->FindBin(centralValue));
// 	  //cout << totData << endl;cout << totSgn << endl;cout << totBkg << endl;

// 	  float oldDataContent = hData->GetBinContent(i);
// 	  float oldDataError   = hData->GetBinError(i);
// 	  hData->SetBinContent(i,  totData+oldDataContent);
// 	  hData->SetBinError(i, TMath::Sqrt( oldDataError*oldDataError + totData*totData) ); // to be checked!
// 	  hDataW->SetBinContent(i,  SoB*totData+oldDataContent);
// 	  hDataW->SetBinError(i, TMath::Sqrt( oldDataError*oldDataError + SoB*SoB*totData*totData) ); // to be checked!

// 	  float oldSgnContent = hSgn->GetBinContent(i);
// 	  float oldSgnError   = hSgn->GetBinError(i);
// 	  hSgn->SetBinContent(i,  totSgn+oldSgnContent);
// 	  hSgn->SetBinError(i, TMath::Sqrt( oldSgnError*oldSgnError + totSgn*totSgn) ); // to be checked!
// 	  hSgnW->SetBinContent(i,  SoB*totSgn+oldSgnContent);
// 	  hSgnW->SetBinError(i, TMath::Sqrt( oldSgnError*oldSgnError + SoB*SoB*totSgn*totSgn) ); // to be checked!
	  
// 	  float oldBkgContent = hBkg->GetBinContent(i);
// 	  float oldBkgError   = hBkg->GetBinError(i);
// 	  hBkg->SetBinContent(i,  totBkg+oldBkgContent);
// 	  hBkg->SetBinError(i, TMath::Sqrt( oldBkgError*oldBkgError + totBkg*totBkg) ); // to be checked!
// 	  hBkgW->SetBinContent(i,  SoB*totBkg+oldBkgContent);
// 	  hBkgW->SetBinError(i, TMath::Sqrt( oldBkgError*oldBkgError + SoB*SoB*totBkg*totBkg) ); // to be checked!
// 	}

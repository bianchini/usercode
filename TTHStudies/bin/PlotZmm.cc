#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>

#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCut.h"
#include "TArrayF.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "Bianchi/TTHStudies/interface/Samples.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"

#define VERBOSE 1

using namespace std;


void ClearAllHisto(map<string, TH1F*> aMap){

  for(  std::map<string,TH1F*>::iterator it = aMap.begin(); it!= aMap.end() ; it++){
    delete (it->second);
  }

  return;

}


void CanvasAndLegend(TCanvas* c1,  TLegend* leg, int logy=0){

  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(logy);

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
  gStyle->SetTitleOffset(1.3,"y");

  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  
}


void DrawHistogramMC(TTree* tree = 0, 
		     TString variable = "", 
		     float& normalization      = *(new float()), 
		     float& normalizationError = *(new float()), 
		     float scaleFactor = 0., 
		     TH1F* h = 0, 
		     TCut cut = TCut(""),
		     int verbose = 0 ){
 
  if(tree!=0 && h!=0){
    h->Reset();
    tree->Draw(variable+">>"+TString(h->GetName()),"(PUweight*weightTrig2012)"*cut); //lheWeight
    h->Scale(scaleFactor);
    normalization      = h->Integral();
    normalizationError = TMath::Sqrt(h->GetEntries())*(normalization/h->GetEntries());
    if(verbose==0) h->Reset();
  }
  else{
    cout << "Function DrawHistogramMC has raised an error" << endl;
    return;
  }
}

void DrawHistogramData(TTree* tree = 0, 
		       TString variable = "", 
		       float& normalization      = *(new float()), 
		       float& normalizationError = *(new float()), 
		       float scaleFactor = 0., 
		       TH1F* h = 0, 
		       TCut cut = TCut(""),
		       int verbose = 0 ){
 
  if(tree!=0 && h!=0){
    h->Reset();
    tree->Draw(variable+">>"+TString(h->GetName()),"((EVENT.json == 1 || EVENT.run < 196532) && H.HiggsFlag==1 && (triggerFlags[14]>0 || triggerFlags[21]>0 || triggerFlags[22]>0 || triggerFlags[23]>0))"*cut);
    h->Scale(scaleFactor);
    normalization      = h->Integral();
    normalizationError = TMath::Sqrt(h->GetEntries())*(normalization/h->GetEntries());
    if(verbose==0) h->Reset();
  }
  else{
    cout << "Function DrawHistogramData has raised an error" << endl;
    return;
  }
}


int main(int argc, const char* argv[])
{

  std::cout << "PlotZmm" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");
  AutoLibraryLoader::enable();

  PythonProcessDesc builder(argv[1]);
 
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;
  const edm::VParameterSet& plots   = in.getParameter<edm::VParameterSet>("plots") ;

  std::string pathToFile( in.getParameter<std::string>("pathToFile" ) );
  std::string ordering(   in.getParameter<std::string>("ordering" ) );
  double lumi(       in.getParameter<double>("lumi" ) );


  Samples* mySamples = new Samples(pathToFile, ordering, samples, lumi, VERBOSE);
  map<string, TH1F*> mapHist;
  vector<string> mySampleFiles;

  if(mySamples->IsOk()){

    cout << "Ok!" << endl;
    mySampleFiles = mySamples->Files();

    for( unsigned int i = 0 ; i < mySampleFiles.size(); i++){
      string sampleName       = mySampleFiles[i];

      if(VERBOSE){
	cout << mySampleFiles[i] << " ==> " << mySamples->GetXSec(sampleName) 
	     << " pb, Num. events  "        << mySamples->GetTree(sampleName, "tree")->GetEntries() 
	     << " ==> weight = "            << mySamples->GetWeight(sampleName) << endl;
      }
    }
  }
  else{
    cout << "Problems... leaving" << endl;
    return 0;
  }

 

  BOOST_FOREACH(edm::ParameterSet pset, plots){

    float totalMC   = 0.;
    float totalData = 0.;
    float normalization      = 0;
    float normalizationError = 0;

    bool skip        = pset.getParameter<bool>("skip");
    float xLow       = pset.getParameter<double>("xLow");
    float xHigh      = pset.getParameter<double>("xHigh");
    int nBins        = pset.getParameter<int>("nBins");
    string var       = pset.getParameter<string>("variable");
    string xTitle    = pset.getParameter<string>("xTitle");
    string yTitle    = pset.getParameter<string>("yTitle");
    string histoName = pset.getParameter<string>("histoName");
    int logy         = pset.getParameter<int>("logy");

    if(skip) continue;

    for(unsigned int i = 0 ; i < mySampleFiles.size(); i++){

      string currentName       = mySampleFiles[i];

      TH1F* h1 = new TH1F( ("h1_"+currentName).c_str(), (" ;"+xTitle+";"+yTitle).c_str(), nBins, xLow, xHigh );
      h1->SetLineColor( mySamples->GetColor(currentName) );
      h1->SetFillColor( mySamples->GetColor(currentName) );

      if( currentName.find("Data")!=string::npos) 
	h1->SetMarkerStyle(20);h1->SetMarkerSize(1.2);h1->SetMarkerColor(kBlack);h1->SetLineColor(kBlack);

      mapHist[currentName] = h1;


      TH1F*  currentHisto      = mapHist[currentName];
      TTree* currentTree       = mySamples->GetTree( currentName, "tree");
      
      if( currentHisto==0 || currentTree==0){ cout << "No histo/file..." << endl; continue; }
      
      TString variable(var.c_str());
     
      float scaleFactor        = mySamples->GetWeight(currentName);

      TCut myCut("Vtype==0 && H.HiggsFlag==1 && vLepton_pt[0]>24 && vLepton_pt[1]>20 && vLepton_charge[0]*vLepton_charge[1]<0 && vLepton_pfCombRelIso[0]<0.10 && vLepton_pfCombRelIso[1]<0.30 && V.mass>60 && hJet_pt[0]>30 && hJet_pt[1]>30");
      
      if( currentName.find("Data")==string::npos ){
	DrawHistogramMC(    currentTree, variable, normalization, normalizationError, scaleFactor, currentHisto, myCut, 1);
	totalMC += normalization;
	if(VERBOSE) cout << "Histo for " << currentName << " has integral " << currentHisto->Integral() << endl;
      }
      else{
	DrawHistogramData(  currentTree, variable, normalization, normalizationError, scaleFactor, currentHisto, myCut, 1);
	totalData += normalization;
	if(VERBOSE) cout << "Histo for " << currentName << " has integral " << currentHisto->Integral() << endl;
      }
      
    }
    
    
    if(VERBOSE) cout << "Toatl MC = "   << totalMC << endl;
    if(VERBOSE) cout << "Toatl Data = " << totalData << endl;
    
    
    TCanvas *c1  = new TCanvas("c1","",5,30,650,600);
    TLegend* leg = new TLegend(0.65,0.45,0.90,0.90,NULL,"brNDC");
    
    CanvasAndLegend(c1, leg, logy);
    THStack* aStack = new THStack("aStack","");

    bool flagSingleTop = false;    
    for(  std::map<string,TH1F*>::iterator it = mapHist.begin(); it!= mapHist.end() ; it++){
      if( (it->first).find("Data")==string::npos  ){
	if( (it->first).find("ZH")==string::npos ||
	    ((it->first).find("ZH")!=string::npos && (it->first).find("125")!=string::npos )  ){
	  aStack->Add((it->second));

	  bool isSingleTop = ((it->first).find("Tt")!=string::npos || (it->first).find("TtW")!=string::npos  || (it->first).find("Ts")!=string::npos ||
			      (it->first).find("Tbar")!=string::npos );

	  if( isSingleTop && !flagSingleTop){
	    leg->AddEntry(it->second, "Single top", "F");
	    flagSingleTop = true;
	  }
	  else if(!isSingleTop){
	    leg->AddEntry(it->second, (it->first).c_str(), "F");
	  }

	}
      }
    }
    
    c1->cd();
    if(mapHist["DataZmm"]!=0){
      leg->AddEntry(mapHist["DataZmm"], "Data", "P");
      mapHist["DataZmm"]->Draw("P");
    }
    aStack->Draw("HISTSAME");
    mapHist["DataZmm"]->Draw("PSAME");
    leg->Draw();
    c1->SaveAs(("Plots/"+histoName+".png").c_str());
    c1->SaveAs(("Plots/"+histoName+".root").c_str());
   


    /// clear
    ClearAllHisto(mapHist); delete c1; delete leg;

  }


  return 0;

}

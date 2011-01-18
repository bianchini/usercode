#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "TList.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TObject.h"
#include "THStack.h"
#include "TLegend.h"
#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <sstream> 
#include <strstream> 
#include <algorithm>

using namespace std;

void setUpHisto(TH1F* h = 0, float Lumi = 1.0, string selection = ""){

  if(h==0) return;
  h->SetAxisRange(30.,130,"X");
  h->SetXTitle("M_{t&p} [GeV/c^{2}]");
  h->SetYTitle( Form("Events / (%.1f GeV/c^{2}) / %.1f pb^{-1}", h->GetBinWidth(1), Lumi ) );
  h->SetMinimum(0.0);
  h->SetTitle(("Tag&probe pair mass for selection "+selection).c_str());

}



void makeTree(string sample = "DY", string selection = "", float Lumi = 1){

  vector<pair<string,float> > files;
  vector< TH1F* > histos;

  vector<string> directories;
  
  if(sample.find("DY")!=string::npos){
    files.push_back( make_pair("DY-madgraph-10to50.root",      (310*0.19*Lumi)  ) );
    files.push_back( make_pair("DY-madgraph-50.root",          (2289*0.56*Lumi)  ) );
  }
  
  if(sample.find("Wjets")!=string::npos){
    files.push_back( make_pair("W1Jets_ptW-0to100.root",       (3.693e+03*0.44*Lumi)  ) );
    files.push_back( make_pair("W1Jets_ptW-100to300.root",     (7.197e+01*0.58*Lumi)  ) );
    files.push_back( make_pair("W1Jets_ptW-300to800.root",     (5.658e-01*0.73*Lumi)  ) );
    files.push_back( make_pair("W2Jets_ptW-0to100.root",       (9.434e+02*0.46*Lumi)  ) );
    files.push_back( make_pair("W2Jets_ptW-100to300.root",     (6.718e+01*0.59*Lumi)  ) );
    files.push_back( make_pair("W2Jets_ptW-300to800.root",     (8.553e-01*0.75*Lumi)  ) );
    files.push_back( make_pair("W3Jets_ptW-0to100.root",       (2.087e+02*0.48*Lumi)  ) );
    files.push_back( make_pair("W3Jets_ptW-100to300.root",     (3.243e+01*0.60*Lumi)  ) );
    files.push_back( make_pair("W3Jets_ptW-300to800.root",     (6.229e-01*0.76*Lumi)  ) );
    files.push_back( make_pair("W4Jets_ptW-0to100.root",       (4.446e+01*0.50*Lumi)  ) );
    files.push_back( make_pair("W4Jets_ptW-100to300.root",     (1.138e+01*0.61*Lumi)  ) );
    files.push_back( make_pair("W4Jets_ptW-300to800.root",     (2.950e-01*0.77*Lumi)  ) );
    files.push_back( make_pair("W5Jets_ptW-0to100.root",       (1.111e+01*0.53*Lumi)  ) );
    files.push_back( make_pair("W5Jets_ptW-100to300.root",     (3.789e+00*0.65*Lumi)  ) );
    files.push_back( make_pair("W5Jets_ptW-300to800.root",     (1.565e-01*0.79*Lumi)  ) );
  }
  if(sample.find("TT")!=string::npos){
    files.push_back( make_pair("TT.root",                      (94*Lumi)  ) );
  }
  if(sample.find("QCD")!=string::npos){
    files.push_back( make_pair("QCD_Pt-0to5.root",             (4.84e+10*6.07E-07*Lumi)) );
    files.push_back( make_pair("QCD_Pt-5to15.root",            (3.68e+10*1.82E-06*Lumi)) );
    files.push_back( make_pair("QCD_Pt-15to30.root",           (8.16e+08*0.000113*Lumi)) );
    files.push_back( make_pair("QCD_Pt-30to50.root",           (5.31e+07*0.00388126*Lumi)) );
    files.push_back( make_pair("QCD_Pt-50to80.root",           (6.36e+06*0.02321727*Lumi)) );
    files.push_back( make_pair("QCD_Pt-80to120.root",          (7.84e+05*0.06100585*Lumi)) );
    files.push_back( make_pair("QCD_Pt-120to170.root",         (1.15e+05*0.11389810*Lumi)) );
    files.push_back( make_pair("QCD_Pt-170to300.root",         (2.43e+04*0.13732413*Lumi)) );
    files.push_back( make_pair("QCD_Pt-300to470.root",         (1.17e+03*0.317390358*Lumi)) );
    files.push_back( make_pair("QCD_Pt_470to600.root",         (7.02e+01*0.431763719*Lumi)) );
    files.push_back( make_pair("QCD_Pt_600to800.root",         (1.56e+01*0.508048972*Lumi)) );
    files.push_back( make_pair("QCD_Pt_800to1000.root",        (1.84e+00*0.385363968*Lumi)) );
    files.push_back( make_pair("QCD_Pt_1000to1400.root",       (3.32e-01*0.661857989*Lumi)) );
    //files.push_back( make_pair("QCD_Pt_1400to1800.root",       (1.09e-02*0.784767648*Lumi)) );
  }

  
  directories.push_back("etoTauMargLooseNoCracks80");
  directories.push_back("etoTauMargLooseNoCracks70");
  directories.push_back("etoTauMargLooseNoCracks60");
  directories.push_back("etoTauMargMediumNoCracks80");
  directories.push_back("etoTauMargMediumNoCracks70");
  directories.push_back("etoTauMargMediumNoCracks60");
  directories.push_back("etoTauMargTightNoCracks80");
  directories.push_back("etoTauMargTightNoCracks70");
  directories.push_back("etoTauMargTightNoCracks60");

  directories.push_back("etoTauSCMargNoCracks80");
  directories.push_back("etoTauSCMargNoCracks70");
  directories.push_back("etoTauSCMargNoCracks60");
    

  TFile *outFile = new TFile(("testNewWriteFromPAT_soup"+selection+".root").c_str(),"RECREATE");
  
  int numFiles = 0;

  for(unsigned int j = 0; j<directories.size(); j++){

    numFiles = 0;

    histos.clear();

    for(unsigned int i = 0; i<files.size();i++){
      histos.push_back( new TH1F(Form("h_%d",i),Form("file %s",(files[i].first).c_str()),20,40,120) );
    }

    THStack* aStack = new THStack("aStack","");
    TLegend* leg = new TLegend(0.1331269,0.5926573,0.3622291,0.8671329,NULL,"brNDC");
    leg->SetFillStyle(4000);
    leg->SetBorderSize(0);
    leg->SetFillColor(10);
    leg->SetTextSize(0.03);

    cout << "Directory " <<  directories[j] << endl;

    TChain* outChain = new TChain((directories[j]+"/fitter_tree").c_str());

    int counter = 0;
    for(unsigned int i = 0; i<files.size();i++){

      TFile *file = new TFile(("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_"+files[i].first).c_str(),"READ");
      if(file->IsZombie()){
	cout << "No file " << files[i].first << " is found" << endl;
	continue;
      }

      file->cd("allEventsFilter");
      TH1F* totalEvents = (TH1F*)gDirectory->Get("totalEvents");
      float readEvents = totalEvents->GetBinContent(1);

      file->cd(directories[j].c_str());
      TTree *oldTree = (TTree*) gDirectory->Get("fitter_tree");

      float scaleToNNLO = 1.0;
      if((files[i].first).find("DY")!=string::npos)    scaleToNNLO = 1.33;
      if((files[i].first).find("W")!=string::npos)     scaleToNNLO = 1.23;
      if((files[i].first).find("TT")!=string::npos)    scaleToNNLO = 1.75;
     
      int entries = std::min( (int)oldTree->GetEntries(),  
			      (int)(((files[i].second*scaleToNNLO)/readEvents)*oldTree->GetEntries()) );

      TFile *newFile = new TFile( ("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_"+files[i].first+"new").c_str(), "RECREATE");
      TDirectory* dir = (TDirectory*) newFile->mkdir(directories[j].c_str());
      TTree *tree = oldTree->CloneTree( entries );

      // weight for the tree: 1 if you require less than the overall tree entries, otherwise set it to L*sigma/processed 
      //float weight = std::max((double)files[i].second,1.0);
      //tree->SetWeight( weight );
      double weight;
      tree->Branch("weight", &weight,"weight/D");
      weight = std::max((double)(files[i].second/readEvents),1.0);
      tree->Fill();

      TH1F* h = new TH1F("h","",20,40,120);
      tree->Draw("mass>>h",("weight*("+selection+")").c_str());
      histos[i]->Add(h,1);
      histos[i]->SetFillColor(i+1);
      histos[i]->SetLineColor(i+1);
      //setUpHisto(histos[i],Lumi, directories[j]);
      aStack->Add( histos[i] );
      leg->AddEntry( histos[i], (files[i].first).c_str(),"F");

      dir->cd();
      tree->Write();
      newFile->Write();

      counter+=entries;

      cout << files[i].first 
	   << " ==> total entries " << oldTree->GetEntries() 
	   << " --- Required " <<  (int)(files[i].second/readEvents*oldTree->GetEntries())
	   << " --- Provided " << entries
	   <<  endl;
      std::cout << "The tree will have weight: " << weight << std::endl; 
      //outTree->CopyEntries(tree,entries);

      if(entries>0){
	int isOK =  outChain->AddFile( ("/data_CMS/cms/lbianchini/35pb/testNewWriteFromPAT_"+files[i].first+"new").c_str());
	if( isOK == 0) cout << "No linked file" << endl;
	numFiles += isOK;
      }
      cout << " outChain has " << outChain->GetEntries() 
	   << " and attached files " << numFiles
	   << endl;

      file->Close();
      newFile->Close();
      //delete h;
      delete file;
      delete newFile;
    } // file

    //outFile->cd();
    TDirectory* dir = (TDirectory*) outFile->mkdir(directories[j].c_str());
    
    TTree* outTree = (TTree*) outChain->CopyTree("");
  
    //outTree->Draw("mass");

    cout<< "Total events copied in output directory " << outTree->GetEntries() << " from requested " << counter << endl;

    dir->cd();
    outTree->Write("", TObject::kOverwrite);

    TCanvas *c1 = new TCanvas("c1","stacked mass",10,30,650,600);
    c1->SetGrid(0,0);
    c1->SetFillStyle(4000);
    c1->SetFillColor(10);
    c1->SetTicky();
    c1->SetObjectStat(0);

    aStack->Draw("HIST");

    setUpHisto( (TH1F*)aStack->GetHistogram(), Lumi, directories[j] );
    leg->SetHeader(directories[j].c_str());
    leg->Draw(); 
    c1->Write();

    for(int i = 0; i<histos.size(); i++){
      delete histos[i];
    }
    delete aStack;
    delete leg;

  }//directories

  


  outFile->Write();
  outFile->Close();
  delete outFile;

}




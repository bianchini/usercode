#include <iostream> 
#include <map>
#include <string>

#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TCut.h"
#include "TGraph.h"


void plot(int steps_ = 30){


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.55,0.65,0.80,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);

  TFile *fSgn = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleVBFH115-Mu-powheg-PUS1_Open_MuTauStream.root","READ");
  TFile *fBkg = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleWJets-Mu-madgraph-PUS1_Open_MuTauStream.root","READ");

  if((fSgn->IsZombie()) || (fBkg->IsZombie()) ){
    cout << "No file" << endl;
    exit(1);
  }

  float vx_Mt[steps_];
  float vy_Mt[steps_];
  float vx_pZeta[steps_];
  float vy_pZeta[steps_];

  TH2F* hMaster = new TH2F("hMaster"," ; qqH efficiency ; W+jets efficiency", 100,0,1,100,0,1);

  TTree* treeSgn = (TTree*)fSgn->Get("outTreePtOrd");
  TTree* treeBkg = (TTree*)fBkg->Get("outTreePtOrd");

  TFile* dummy = new TFile("dummy.root","RECREATE");
  TTree* treeSgnCut = (TTree*)treeSgn->CopyTree("tightestHPSWP>0 && pt1>25 && pt2>20 && eta1*eta2<0");  
  TTree* treeBkgCut = (TTree*)treeBkg->CopyTree("tightestHPSWP>0 && pt1>25 && pt2>20 && eta1*eta2<0");  
  //TTree* treeSgnCut = (TTree*)treeSgn->CopyTree("tightestHPSWP>0");  
  //TTree* treeBkgCut = (TTree*)treeBkg->CopyTree("tightestHPSWP>0");  

  cout << treeSgnCut->GetEntries() << endl;
  cout << treeBkgCut->GetEntries() << endl;

  std::vector<string> antiW;
  antiW.push_back("MtLeg1");
  antiW.push_back("pZetaCutVar");

  TH1F* h1 = new TH1F("h1","",1,-10,10);

  for(unsigned int m = 0; m < antiW.size() ; m++){
    
    treeSgnCut->Draw( "eta1>>h1" ,"puWeight");
    float sgnDen = (float)h1->Integral();
    cout << sgnDen << endl;
    h1->Reset();
    treeBkgCut->Draw( "eta1>>h1" ,"puWeight");
    float bkgDen = (float)h1->Integral();
    cout << bkgDen << endl;
    h1->Reset();
    
    for(int i = 0; i < steps_; i++){
	
      cout << antiW[m] << " --> " << i << endl;
      
      float step_i = antiW[m].find("Mt")!=string::npos ? 
	i*(60./steps_) : i*(60/steps_);

      TCut cut1( Form("puWeight*(%s<=%f)",  antiW[m].c_str(), step_i ) );
      TCut cut2( Form("puWeight*(%s>=(20-%f))", antiW[m].c_str(), step_i ) );

      TCut cut = antiW[m].find("Mt")!=string::npos ? cut1 : cut2;

      treeSgnCut->Draw( "eta1>>h1" ,cut );
      float sgnNum = (float)h1->Integral();
      cout << sgnNum << endl;
      h1->Reset();
      treeBkgCut->Draw( "eta1>>h1" ,cut );
      float bkgNum = (float)h1->Integral();
      cout << bkgNum << endl;
      h1->Reset();
 
      if(m==0){
	vx_Mt[i] = sgnNum/sgnDen;
	vy_Mt[i] = bkgNum/bkgDen;
      }
      if(m==1){
	vx_pZeta[i] = sgnNum/sgnDen;
	vy_pZeta[i] = bkgNum/bkgDen;
      }
    }
  }

  leg->SetHeader("ROC for anti-W cuts after VBF-preselection");

  hMaster->Draw();
  
  TGraph* graph_Mt    = new TGraph(steps_,vx_Mt,vy_Mt);
  TGraph* graph_pZeta = new TGraph(steps_,vx_pZeta,vy_pZeta);

  std::map<std::string, TGraph*> roc;
  roc["MtLeg1"]      = graph_Mt;
  roc["pZetaCutVar"] = graph_pZeta;

  for(unsigned int k = 0; k < antiW.size() ; k++){
    if(k==0){
      roc[antiW[k]]->SetMarkerColor(kBlue);
      roc[antiW[k]]->SetMarkerSize(1.4);
      roc[antiW[k]]->SetMarkerStyle(kOpenCircle);
      roc[antiW[k]]->GetYaxis()->SetRangeUser(0.80,1.1);
      leg->AddEntry(roc[antiW[k]],"Mt","P");
      roc[antiW[k]]->Draw("P");
    }
    if(k==1){
      roc[antiW[k]]->SetMarkerColor(kRed);
      roc[antiW[k]]->SetMarkerSize(1.4);
      roc[antiW[k]]->SetMarkerStyle(kOpenSquare);
      leg->AddEntry(roc[antiW[k]],"pZeta","P");
      roc[antiW[k]]->Draw("PSAME");
    }
  }

  leg->Draw();
  
}

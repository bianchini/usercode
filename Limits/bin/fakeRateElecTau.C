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
#include "TFitResultPtr.h"


void fakeRate(){

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);


  float Lumi   = (-47.4 + 215.6 + 955.3 + 389.9 + 706.7 + 2714);
  float lumiCorrFactor = (1-0.056);
  Lumi *= lumiCorrFactor;


  TFile *fData              
    = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamFall11_04May2012_PreApproval/nTupleRun2011-ElecTau-All_run_Open_ElecTauStream.root", "READ");  

  TFile *fBackgroundWJets   
    = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamFall11_04May2012_PreApproval/nTupleWJets-ElecTau-madgraph-PUS6_run_Open_ElecTauStream.root","READ"); 
  TFile *fBackgroundDY
    = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamFall11_04May2012_PreApproval/nTupleDYJets-ElecTau-50-madgraph-PUS6_run_Open_ElecTauStream.root","READ"); 
  TFile *fBackgroundTTbar  
    = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamFall11_04May2012_PreApproval/nTupleTTJets-ElecTau-madgraph-PUS6_run_Open_ElecTauStream.root","READ"); 
  TFile *fBackgroundOthers  
    = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/ElecTauStreamFall11_04May2012_PreApproval/nTupleOthers-ElecTau-PUS6_run_Open_ElecTauStream.root","READ"); 
  
  TString tree         = "outTreePtOrd";

  TTree *data          = (TTree*)fData->Get(tree);

  TTree *backgroundTTbar     = (TTree*)fBackgroundTTbar->Get(tree);
  TTree *backgroundWJets     = (TTree*)fBackgroundWJets->Get(tree);
  TTree *backgroundOthers    = (TTree*)fBackgroundOthers->Get(tree);
 
  TFile* dummy1 = new TFile("dummy1.root","RECREATE");
  cout << "Now copying g/Z -> tau+ tau- " << endl;
  TTree *backgroundDYTauTau  = ((TTree*)fBackgroundDY->Get(tree))->CopyTree("abs(genDecay)==(23*15)");                 // g/Z -> tau+ tau-
  cout << "Now copying g/Z -> e+e- e->tau" << endl;
  TTree *backgroundDYEtoTau = ((TTree*)fBackgroundDY->Get(tree))->CopyTree("abs(genDecay)!=(23*15) &&  leptFakeTau"); // g/Z -> mu+mu- mu->tau
  cout << "Now copying g/Z -> e+e- jet->tau" << endl;
  TTree *backgroundDYJtoTau  = ((TTree*)fBackgroundDY->Get(tree))->CopyTree("abs(genDecay)!=(23*15) && !leptFakeTau"); // g/Z -> mu+mu- jet->tau


  bool useMt      = true;
  string antiWcut = useMt ? "MtLeg1MVA" : "-(pZetaMVA-1.5*pZetaVisMVA)" ;
  float antiWsgn  = useMt ? 40. :  20. ; 
  float antiWsdb  = useMt ? 60. :  40. ; 

  TCut lpt("ptL1>20 && TMath::Abs(etaL1)<2.1");
  TCut tpt("ptL2>20 && ptL2<999 && TMath::Abs(etaL2)<2.3");

  TCut lID("((TMath::Abs(etaL1)<0.925 && mvaPOGNonTrig>0.85) || (TMath::Abs(etaL1)<1.479 && TMath::Abs(etaL1)>0.80 && mvaPOGNonTrig>0.975) || (TMath::Abs(etaL1)>1.479 && mvaPOGNonTrig>0.985))");
  tpt = tpt&&lID;

  TCut tiso("tightestHPSMVAWP>=0 && tightestAntiECutWP>1"); 
  TCut liso("combRelIsoLeg1DBeta<0.10");
  TCut laiso("combRelIsoLeg1DBetav2>0.20 && combRelIsoLeg1DBetav2<0.50");
  TCut lliso("combRelIsoLeg1DBetav2<0.20");
  TCut lveto("elecFlag==0");
  TCut SS("diTauCharge!=0");
  TCut OS("diTauCharge==0");
  TCut hltevent("pairIndex<1 && HLTx==1 && (run>=163269 || run==1)");
  TCut hltmatch("HLTmatch==1");

  TCut pZ( Form("((%s)<%f)",antiWcut.c_str(),antiWsgn));
  TCut apZ(Form("((%s)>%f)",antiWcut.c_str(),antiWsdb));

 
  vector<int> bins;
  bins.push_back(20);
  bins.push_back(22);
  bins.push_back(24);
  bins.push_back(26);
  bins.push_back(28);
  bins.push_back(30);
  bins.push_back(32);
  bins.push_back(34);
  bins.push_back(36);
  bins.push_back(40);
  bins.push_back(45);
  bins.push_back(50);
  bins.push_back(60); 
  bins.push_back(80); 
  bins.push_back(100); 

  int nBins =  bins.size() ;
  TArrayF binsT(nBins);

  for(unsigned int k = 0 ; k < nBins ; k++ )
    binsT[k] = bins[k];

  TH1F* hFakeRate  = new TH1F("hFakeRate", "", nBins-1, binsT.GetArray() );

  TH1F* hFakeRateErrUp    = new TH1F("hFakeRateErrUp",   "", nBins-1, binsT.GetArray() );
  TH1F* hFakeRateErrDown  = new TH1F("hFakeRateErrDown", "", nBins-1, binsT.GetArray() );

  TH1F* hPuritySdb = new TH1F("hPuritySdb","", nBins-1, binsT.GetArray() );
  TH1F* hPurityQCD = new TH1F("hPurityQCD","", nBins-1, binsT.GetArray() );
  TH1F* hPurityAIs = new TH1F("hPurityAIs","", nBins-1, binsT.GetArray() );

  //TH1F* hFakeRate  = new TH1F("hFakeRate", "", nBins-1, 20,100 );
  //TH1F* hPuritySdb = new TH1F("hPuritySdb","", nBins-1, 20,100 );
  //TH1F* hPurityQCD = new TH1F("hPurityQCD","", nBins-1, 20,100 );

  TH1F* hExtrap        = new TH1F("hExtrap","",1, -10,10);
  TH1F* hSVfit         = new TH1F("hSVfit",     "",40, 0,400);
  TH1F* hSVfitFake     = new TH1F("hSVfitFake", "",40, 0,400);
  TH1F* hSVfitAIso     = new TH1F("hSVfitAIso", "",40, 0,400);
  TH1F* hSVfitHelpAdd  = new TH1F("hSVfitHelpAdd", "",40, 0,400);
  TH1F* hSVfitHelp     = new TH1F("hSVfitHelp", "",40, 0,400);

  TCut sbinSSaIsoInclusiveAllPt = lpt && tpt && tiso && laiso && lveto && SS && pZ  && hltevent && hltmatch;

  for(int i = 0; i <  bins.size()-1 ; i++){

    float SSIsoToSSAIsoRatioQCD = 1.0;

    float min = bins[i];
    float max = bins[i+1]  ;

    TCut lpt_i = lpt && TCut(Form("ptL1>%f && ptL1<=%f", min, max));

    TCut sbinInclusive;
    sbinInclusive            = lpt_i && tpt && tiso && liso && lveto && OS && pZ  && hltevent && hltmatch;
    TCut sbinPZetaRelSSInclusive;
    sbinPZetaRelSSInclusive  = lpt_i && tpt && tiso && liso && lveto && SS        && hltevent && hltmatch;
    TCut sbinSSInclusive;
    sbinSSInclusive          = lpt_i && tpt && tiso && liso && lveto && SS && pZ  && hltevent && hltmatch;
    TCut sbinSSaIsoInclusive;
    sbinSSaIsoInclusive      = lpt_i && tpt && tiso && laiso&& lveto && SS && pZ  && hltevent && hltmatch;


    cout << "******** Extrapolation factors for QCD normalization: " << " bin " << min << "," << max << " " << "********" << endl;

    hExtrap->Reset();
    backgroundWJets->Draw("etaL1>>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelSSInclusive&&pZ));
    float ExtrapSSWinSignalRegionMC   = hExtrap->Integral();
    hExtrap->Reset();
    backgroundWJets->Draw("etaL1>>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelSSInclusive&&apZ));
    float ExtrapSSWinSidebandRegionMC = hExtrap->Integral();
    float ExtrapscaleFactorSS         = ExtrapSSWinSignalRegionMC>0 ? ExtrapSSWinSidebandRegionMC/ExtrapSSWinSignalRegionMC : 1.0;
    cout << " Extrapolation factor W SS (inclusive) " << ExtrapscaleFactorSS << endl;
    
    hExtrap->Reset();
    backgroundTTbar->Draw("etaL1>>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelSSInclusive&&apZ));
    float ExtrapttbarExtrSS    = hExtrap->Integral()*Lumi/1000;
    hExtrap->Reset();
    backgroundOthers->Draw("etaL1>>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelSSInclusive&&apZ));
    float ExtrapothersExtrSS   = hExtrap->Integral()*Lumi/1000;
    hExtrap->Reset();
    backgroundDYJtoTau->Draw("etaL1>>hExtrap","(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*(sbinPZetaRelSSInclusive&&apZ));
    float ExtrapdyjtotauExtrSS = hExtrap->Integral()*Lumi/1000;


    hExtrap->Reset();
    data->Draw("etaL1>>hExtrap", sbinPZetaRelSSInclusive&&apZ);
    float ExtrapSSWinSignalRegionDATA = hExtrap->Integral();
    cout << "Extrapolation for QCD (inclusive): total data events in sideband " << ExtrapSSWinSignalRegionDATA << endl;

    float totalBkgSdb = 
      ExtrapttbarExtrSS +
      ExtrapothersExtrSS +
      ExtrapdyjtotauExtrSS;

    hPuritySdb->SetBinContent(i+1, totalBkgSdb/ExtrapSSWinSignalRegionDATA);


    ExtrapSSWinSignalRegionDATA -= ExtrapttbarExtrSS;
    ExtrapSSWinSignalRegionDATA -= ExtrapothersExtrSS;
    ExtrapSSWinSignalRegionDATA -= ExtrapdyjtotauExtrSS;
    ExtrapSSWinSignalRegionDATA /= ExtrapscaleFactorSS;
    cout << "Extrapolation for QCD (inclusive): W+jets in SS signal region (inclusive) is estimated to be " << ExtrapSSWinSignalRegionDATA << endl;
    

    hSVfitHelpAdd->Reset(); 

    float totalBkg =  0.;
    hExtrap->Reset();
    data->Draw("etaL1>>hExtrap", sbinSSInclusive);
    data->Draw("diTauNSVfitMass>>hSVfitHelpAdd", sbinSSInclusive);
    hSVfit->Add(hSVfitHelpAdd,+1);

    float SSeventsExtrapTot = hExtrap->Integral();
    float SSeventsExtrap    = hExtrap->Integral();

    cout << "Extrapolation for SS events in data (inclusive) " << hExtrap->GetEntries() << endl;
    cout << "Subtracting W+jets (SS)..." << endl;
    SSeventsExtrap  -= ExtrapSSWinSignalRegionDATA;
    totalBkg        += ExtrapSSWinSignalRegionDATA;

    hSVfitHelp->Reset();

    backgroundWJets->Draw("diTauNSVfitMass>>hSVfitHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*sbinSSInclusive);
    if(hSVfitHelp->Integral()>0) hSVfitHelp->Scale(ExtrapSSWinSignalRegionDATA/hSVfitHelp->Integral());
    hSVfit->Add(hSVfitHelp,-1);
    hSVfitHelp->Reset();

    hExtrap->Reset();
    backgroundTTbar->Draw("etaL1>>hExtrap", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*sbinSSInclusive);
    backgroundTTbar->Draw("diTauNSVfitMass>>hSVfitHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*sbinSSInclusive);
    SSeventsExtrap  -= hExtrap->Integral()*Lumi/1000;
    totalBkg        += hExtrap->Integral()*Lumi/1000;
    hSVfitHelp->Scale(Lumi/1000);
    hSVfit->Add(hSVfitHelp,-1);
    hSVfitHelp->Reset();

    hExtrap->Reset();
    backgroundDYEtoTau->Draw("etaL1>>hExtrap", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*sbinSSInclusive);
    backgroundDYEtoTau->Draw("diTauNSVfitMass>>hSVfitHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*sbinSSInclusive);
    SSeventsExtrap  -= hExtrap->Integral()*Lumi/1000;
    totalBkg        += hExtrap->Integral()*Lumi/1000;
    hSVfitHelp->Scale(Lumi/1000);
    hSVfit->Add(hSVfitHelp,-1);
    hSVfitHelp->Reset();

    hExtrap->Reset();
    backgroundDYJtoTau->Draw("etaL1>>hExtrap", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*sbinSSInclusive);
    backgroundDYJtoTau->Draw("diTauNSVfitMass>>hSVfitHelp", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*sbinSSInclusive);
    SSeventsExtrap  -= hExtrap->Integral()*Lumi/1000;
    totalBkg        += hExtrap->Integral()*Lumi/1000;
    hExtrap->Reset();
    hSVfitHelp->Scale(Lumi/1000);
    hSVfit->Add(hSVfitHelp,-1);
    hSVfitHelp->Reset();

    data->Draw("etaL1>>hExtrap", sbinSSaIsoInclusive);
    float SSeventsExtrapAiso = hExtrap->GetEntries();

    backgroundWJets->Draw("etaL1>>hExtrap", "(sampleWeight*puWeight*HLTweightTau*HLTweightElec*SFTau*SFElec)"*sbinSSaIsoInclusive);
    hExtrap->Scale(Lumi/1000);
    hPurityAIs->SetBinContent(i+1, hExtrap->Integral()/SSeventsExtrapAiso );
    hExtrap->Reset();
    SSeventsExtrapAiso -= hExtrap->Integral();

    SSIsoToSSAIsoRatioQCD = SSeventsExtrap>0 ? 1./(SSeventsExtrap/SSeventsExtrapAiso) : 1.;
    cout << "The extrapolation factor Iso<0.1 / 0.3<Iso<0.5 is " << SSIsoToSSAIsoRatioQCD << endl;
    

    hPurityQCD->SetBinContent(i+1,totalBkg/(SSeventsExtrap+totalBkg));


    float error_iDen = 0.;
    if(SSeventsExtrapAiso>0) error_iDen += 1./TMath::Sqrt(SSeventsExtrapAiso);

    float error_iNum2 = 0.;
    error_iNum2 += SSeventsExtrap;
    error_iNum2 += (totalBkg*0.06)*(totalBkg*0.06);
    float error_iNum = SSeventsExtrap>0 ? TMath::Sqrt(error_iNum2)/SSeventsExtrap : 0.0;

    hFakeRate->SetBinContent(i+1, SSIsoToSSAIsoRatioQCD);
    hFakeRate->SetBinError(i+1,   (error_iDen+error_iNum)*SSIsoToSSAIsoRatioQCD );

    cout <<  hFakeRate->GetBinContent(i+1) << endl;
    cout <<  hFakeRate->GetBinError(i+1) << endl;
    cout << "************** END extrapolation *******************" << endl;

  }

  //return;

  TF1* fit = new TF1("fit","[0]*TMath::Exp([1]*x)+[2]",20, bins[bins.size()-1]);
  fit->SetLineColor(kRed);
  fit->SetParLimits(0,0,20);
  fit->SetParLimits(1,-5,0);
  fit->SetParLimits(2,0,1);

  //TF1* fit = new TF1("fit","[0]+[1]*x+[2]*x*x + [3]*x*x*x",20, bins[bins.size()-2]);
  //fit->SetLineColor(kRed);
  //fit->SetParLimits(0,0,10);
  //fit->SetParLimits(1,-5,5);
  //fit->SetParLimits(2,0,1);
  //fit->SetParLimits(3,-1,1);


  //TF1 *fit = new TF1("fit","[0]+[1]*TMath::Erf([2]*x+[3])",20, bins[bins.size()-1]);
  //fit->SetParLimits(0,-200,200);
  //fit->SetParLimits(1,0,500);
  //fit->SetParLimits(2,0,0.5);
  //fit->SetParLimits(3,-100,1);


  hFakeRate->Fit("fit", "", "", 20, 100);

  float par0 = fit->GetParameter(0);
  float par1 = fit->GetParameter(1);
  float par2 = fit->GetParameter(2);

  float parE0 = fit->GetParError(0);
  float parE1 = fit->GetParError(1);
  float parE2 = fit->GetParError(2);



  for(int i = 0; i <  bins.size()-1 ; i++){

    float min = bins[i];
    float max = bins[i+1]  ;

    double integE    = fit->IntegralError(min,max );
    double integ     = fit->Integral(min,max);
    double binWidth  = max-min;
    cout << "bin " << min << "," << max << " ==> " << integ/binWidth << endl;
    hFakeRateErrUp->SetBinContent( i+1,   (integ + 0.5*integE)/binWidth );
    hFakeRateErrDown->SetBinContent( i+1, (integ - 0.5*integE)/binWidth );
  }
  hFakeRateErrDown->SetLineColor(kRed);
  hFakeRateErrDown->SetLineStyle(kDashed);
  hFakeRateErrUp->SetLineColor(kRed);
  hFakeRateErrUp->SetLineStyle(kDotted);

  cout << "error @ 30 = " << TMath::Sqrt(parE0*parE0 + parE1*parE1 + parE2*parE2) << endl;

  TF1* fitE1 = new TF1("fitE1",Form("(%f)*TMath::Exp((%f)*x)+(%f)", par0+parE0, par1, par2),20, bins[bins.size()-1]);
  TF1* fitE2 = new TF1("fitE2",Form("(%f)*TMath::Exp((%f)*x)+(%f)", par0, par1+parE1, par2),20, bins[bins.size()-1]);
  TF1* fitE3 = new TF1("fitE3",Form("(%f)*TMath::Exp((%f)*x)+(%f)", par0, par1, par2+parE2),20, bins[bins.size()-1]);
  TF1* fitE4 = new TF1("fitE4",Form("(%f)*TMath::Exp((%f)*x)+(%f)", par0-parE0, par1, par2),20, bins[bins.size()-1]);
  TF1* fitE5 = new TF1("fitE5",Form("(%f)*TMath::Exp((%f)*x)+(%f)", par0, par1-parE1, par2),20, bins[bins.size()-1]);
  TF1* fitE6 = new TF1("fitE6",Form("(%f)*TMath::Exp((%f)*x)+(%f)", par0, par1, par2-parE2),20, bins[bins.size()-1]);

  fitE1->SetLineColor(kBlue);
  fitE2->SetLineColor(kGreen);
  fitE3->SetLineColor(kMagenta);
  fitE4->SetLineColor(kYellow);
  fitE5->SetLineColor(kBlack);
  fitE6->SetLineColor(kYellow-6);

  string scaleFact = "( ";
  for(int i = 0; i < bins.size()-1; i++){
    
    float min = bins[i];
    float max = bins[i+1];
    cout << min << "," << max << endl;
    float bin = hFakeRate->FindBin((max+min)/2.);
    //float bin = i+1;
    cout << bin << endl;
    float weightBin_i =  fit->Eval( (max+min)/2.);
    //float weightBin_i =  hFakeRate->GetBinContent( bin );
    
    scaleFact += string( Form("(ptL1>=%f && ptL1<%f)*%f", min , max, 1./weightBin_i ) );
    if(i < bins.size() - 2 ) scaleFact += " + ";
  }
  scaleFact += " )";
  cout << scaleFact << endl;
  
  data->Draw("diTauNSVfitMass>>hSVfitAIso", sbinSSaIsoInclusiveAllPt);
  data->Draw("diTauNSVfitMass>>hSVfitFake", (TCut(scaleFact.c_str()))*sbinSSaIsoInclusiveAllPt);
  hSVfitFake->SetMarkerStyle(kOpenCircle);
  hSVfitFake->SetMarkerSize(1.2);
  hSVfitFake->SetMarkerColor(kRed);

  hSVfitAIso->SetLineColor(kMagenta);

  hSVfit->Print("all");
  hSVfitFake->Print("all");


  hSVfit->DrawNormalized();
  hSVfitFake->DrawNormalized("PSAME");
  hSVfitAIso->DrawNormalized("SAME");
  //return;

  TCanvas *c2 = new TCanvas("c2","",5,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  hFakeRate->SetLineColor(kRed);
  hFakeRate->Draw("E");

  hFakeRateErrDown->Draw("SAME");
  hFakeRateErrUp->Draw("SAME");

//   fitE1->Draw("SAME");
//   fitE2->Draw("SAME");
//   fitE3->Draw("SAME");
//   fitE4->Draw("SAME");
//   fitE5->Draw("SAME");
//   fitE6->Draw("SAME");

  hPurityQCD->SetLineColor(kBlack);
  hPurityQCD->Draw("SAME");
  hPuritySdb->SetLineColor(kBlue);
  hPuritySdb->Draw("SAME");
  hPurityAIs->SetLineColor(kGreen);
  hPurityAIs->Draw("SAME");


}

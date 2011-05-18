


{


  float Lumi = 24.86+159.15;

  TCut hlt("( ((HLTmu==1 && run<=163261) || (HLTx==1 && run>163261)) )");
  
  TCut signal_region_data("sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==0 && MtLeg1<40 && diTauCharge==0 )");
  TCut signal_region("sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==0 && MtLeg1<40 && diTauCharge==0)");
  TCut signal_SS_region("sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==0 && MtLeg1<40 && diTauCharge!=0)");
  TCut Wenrich_OS_region("sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==0 && MtLeg1>60 && diTauCharge==0)");
  TCut Wenrich_SS_region("sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==0 && MtLeg1>60 && diTauCharge!=0)");
  TCut Wrelaxed_OS_region("sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==0 && diTauCharge==0)");
  TCut Wrelaxed_SS_region("sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==0 && diTauCharge!=0)");
  TCut QCDenrich_OS_region("sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta>0.15 && muFlag==0 && MtLeg1<40 && diTauCharge==0)");
  TCut QCDenrich_SS_region("sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta>0.15 && muFlag==0 && MtLeg1<40 && diTauCharge!=0)");
  TCut ZMuMuFakeJet_OS_region("sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==1 && MtLeg1<40 && diTauCharge==0)");
  TCut ZMuMuFakeJet_SS_region("sampleWeight*puWeight*(tightestHPSWP>0 && combRelIsoLeg1DBeta<0.1 && muFlag==1 && MtLeg1<40 && diTauCharge!=0) ");

  TFile *fData = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2//Inclusive/nTupleRun2011-Mu_Open_MuTauStream.root","READ");
  TTree* treeData = (TTree*)fData->Get("outTreePtOrd");

  TFile *fDiBoson = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleDiBoson-Mu_Open_MuTauStream.root","READ");
  TTree* treeDiBoson = (TTree*)fDiBoson->Get("outTreePtOrd");

  TFile *fTTbar = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleTTJets-Mu-madgraph-PUS1_Open_MuTauStream.root","READ");
  TTree* treeTTbar = (TTree*)fTTbar->Get("outTreePtOrd");

  TFile *fSingleTop = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleSingleTop-Mu_Open_MuTauStream.root","READ");
  TTree* treeSingleTop = (TTree*)fSingleTop->Get("outTreePtOrd");

  TFile *fWJets = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleWJets-Mu-madgraph-PUS1_Open_MuTauStream.root","READ");
  TTree* treeWJets = (TTree*)fWJets->Get("outTreePtOrd");

  TFile *fZTauTau = new TFile("/data_CMS/cms/lbianchini/VbfJetsStudy/OpenNtuples/MuTauStream2011_iter2/Inclusive/nTupleDYJets-Mu-50-madgraph-PUS1_Open_MuTauStream.root","READ");
  TTree* treeZTauTauAll = (TTree*)fZTauTau->Get("outTreePtOrd");
  TFile* dummy1 = new TFile("dummy1.root","RECREATE");
  TTree *treeZTauTau = (TTree*)treeZTauTauAll->CopyTree("isTauLegMatched>0.5");
  TTree *treeZFakes  = (TTree*)treeZTauTauAll->CopyTree("isTauLegMatched<0.5");


  TH1F* h1 = new TH1F("h1","",1,-10,10);

  //////////////// DiBoson + TTbar + SingleTop
 
  treeDiBoson->Draw("etaL1>>h1",signal_region);
  float N_DiBoson_signal_region = (float)h1->Integral()*(Lumi/1000.);
  h1->Reset();

  treeDiBoson->Draw("etaL1>>h1",Wenrich_OS_region);
  float N_DiBoson_Wenrich_OS_region = (float)h1->Integral()*(Lumi/1000.);
  h1->Reset();

  treeTTbar->Draw("etaL1>>h1",signal_region);
  float N_TTbar_signal_region = (float)h1->Integral()*(Lumi/1000.);
  h1->Reset();

  treeTTbar->Draw("etaL1>>h1",Wenrich_OS_region);
  float N_TTbar_Wenrich_OS_region = (float)h1->Integral()*(Lumi/1000.);
  h1->Reset();

  treeSingleTop->Draw("etaL1>>h1",signal_region);
  float N_SingleTop_signal_region = (float)h1->Integral()*(Lumi/1000.);
  h1->Reset();

  treeSingleTop->Draw("etaL1>>h1",Wenrich_OS_region);
  float N_SingleTop_Wenrich_OS_region = (float)h1->Integral()*(Lumi/1000.);
  h1->Reset();


  cout << "N_DiBoson_signal_region = " <<  N_DiBoson_signal_region << endl;
  cout << "N_DiBoson_Wenrich_OS_region = " << N_DiBoson_Wenrich_OS_region << endl;
  cout << "N_TTbar_signal_region = " << N_TTbar_signal_region << endl;
  cout << "N_TTbar_Wenrich_OS_region = " << N_TTbar_Wenrich_OS_region << endl;
  cout << "N_SingleTop_signal_region = " << N_SingleTop_signal_region << endl;
  cout << "N_SingleTop_Wenrich_OS_region = " << N_SingleTop_Wenrich_OS_region << endl;

  cout << "//////////////////////////////////////" << endl;

  //////////////// Z->mumu, jet->tau, mu->tau
  
  treeData->Draw("etaL1>>h1",ZMuMuFakeJet_SS_region&&hlt);
  float N_ZMuMuFakeJet_signal_region = (float)h1->Integral();
  h1->Reset();
  // correct for A and eff (from AN-10-430)
  float eA = 0.57;
  float eI = 0.976;
  N_ZMuMuFakeJet_signal_region *= (  (1-eA+eA*(1-eI))/(eA*eI) );
  float N_ZMuMuFakeJet_signal_SS_region =  N_ZMuMuFakeJet_signal_region;

  // apply factor 2 ratio between Z->mumu, jet->tau and Z->mumu, mu->tau
  float N_ZMuMuFakeMu_signal_region = N_ZMuMuFakeJet_signal_region * 2;

  treeZFakes->Draw("etaL1>>h1",signal_region);
  float Exp_ZFake_signal_region = (float)h1->Integral();
  Exp_ZFake_signal_region *= (Lumi/1000.*0.92);
  h1->Reset();

  cout << "N_ZMuMuFakeJet_signal_region = " << N_ZMuMuFakeJet_signal_region << endl;
  cout << "N_ZMuMuFakeJet_signal_SS_region = " << N_ZMuMuFakeJet_signal_SS_region << endl;
  cout << "N_ZMuMuFakeMu_signal_region = "  << N_ZMuMuFakeMu_signal_region << endl;
  cout << "Exp_ZFake_signal_region = "  << Exp_ZFake_signal_region << endl;

  cout << "//////////////////////////////////////" << endl;


  /////////////// W for OS
  treeData->Draw("etaL1>>h1", Wenrich_OS_region&&hlt);
  float N_WMuNu_Wenrich_OS_region = (float)h1->Integral();
  h1->Reset();

  TH1F* h1_Mt = new TH1F("h1_Mt","",200,0,200);
  treeWJets->Draw("MtLeg1>>h1_Mt", Wrelaxed_OS_region);
  float WratioOS = (h1_Mt->Integral(61,120))/(h1_Mt->Integral(0,40));
  h1_Mt->Reset();

  N_WMuNu_Wenrich_OS_region -= N_DiBoson_Wenrich_OS_region;
  N_WMuNu_Wenrich_OS_region -= N_TTbar_Wenrich_OS_region;
  N_WMuNu_Wenrich_OS_region -= N_SingleTop_Wenrich_OS_region;

  float N_WMuNu_signal_region = N_WMuNu_Wenrich_OS_region*(1./WratioOS);

  float WTauWMuRatio = 0.267;
  float N_WTauNu_signal_region = N_WMuNu_signal_region*WTauWMuRatio;


  cout << "WratioOS = " << WratioOS << endl;
  cout << "N_WMuNu_Wenrich_OS_region = " << N_WMuNu_Wenrich_OS_region << endl;
  cout << "N_WMuNu_signal_region = "  <<N_WMuNu_signal_region  << endl;
  cout << "N_WTauNu_signal_region = " << N_WTauNu_signal_region << endl;
  cout << "N_WTauNu_signal_region = " << N_WTauNu_signal_region << endl;

  cout << "//////////////////////////////////////" << endl;

  /////////////// W for SS
  treeData->Draw("etaL1>>h1", Wenrich_SS_region&&hlt);
  float N_WMuNu_Wenrich_SS_region = (float)h1->Integral();
  h1->Reset();

  treeWJets->Draw("MtLeg1>>h1_Mt", Wrelaxed_SS_region);
  float WratioSS = (h1_Mt->Integral(61,120))/(h1_Mt->Integral(0,40));
  h1_Mt->Reset();

  float N_WMuNu_signal_SS_region  = N_WMuNu_Wenrich_SS_region*(1./WratioSS);
  float N_WTauNu_signal_SS_region = N_WMuNu_signal_SS_region*WTauWMuRatio;

  cout << "WratioSS = " << WratioSS << endl;
  cout << "N_WMuNu_signal_SS_region = " << N_WMuNu_signal_SS_region << endl;
  cout << "N_WTauNu_signal_SS_region = " <<  N_WTauNu_signal_SS_region << endl;

  cout << "//////////////////////////////////////" << endl;


  //QCD
  treeData->Draw("etaL1>>h1",signal_SS_region&&hlt);
  float N_QCD_signal_SS_region = (float)h1->Integral();
  h1->Reset();

  N_QCD_signal_SS_region -= N_WMuNu_signal_SS_region; 
  N_QCD_signal_SS_region -= N_WTauNu_signal_SS_region;
  N_QCD_signal_SS_region -= N_ZMuMuFakeJet_signal_SS_region;

  treeData->Draw("etaL1>>h1",QCDenrich_OS_region&&hlt);
  float N_QCD_QCDenrich_OS_region = (float)h1->Integral();
  h1->Reset();
  treeData->Draw("etaL1>>h1",QCDenrich_SS_region&&hlt);
  float N_QCD_QCDenrich_SS_region = (float)h1->Integral();
  h1->Reset();
  
  float SSOSratio = N_QCD_QCDenrich_OS_region/N_QCD_QCDenrich_SS_region;
  float N_QCD_signal_region = N_QCD_signal_SS_region * SSOSratio;

  cout << "SSOSratio = " << SSOSratio  << endl;
  cout << "N_QCD_signal_SS_region = " << N_QCD_signal_SS_region << endl;
  cout << "N_QCD_signal_region = " << N_QCD_signal_region << endl;

  cout << "//////////////////////////////////////" << endl;

 
  treeData->Draw("etaL1>>h1",signal_region&&hlt);
  float N_Z_signal_region = (float)h1->Integral();
  h1->Reset();

  N_Z_signal_region -= N_QCD_signal_region;
  N_Z_signal_region -= N_DiBoson_signal_region;
  N_Z_signal_region -= N_TTbar_signal_region;
  N_Z_signal_region -= N_SingleTop_signal_region;
  N_Z_signal_region -= N_ZMuMuFakeJet_signal_region;
  N_Z_signal_region -= N_ZMuMuFakeMu_signal_region;
  N_Z_signal_region -= N_WMuNu_signal_region;
  N_Z_signal_region -= N_WTauNu_signal_region;

  treeZTauTau->Draw("etaL1>>h1",signal_region);
  float Exp_Z_signal_region = (float)h1->Integral();
  Exp_Z_signal_region *= (Lumi/1000.*0.92);
  h1->Reset();

  cout << "Estimated ==> " << N_Z_signal_region << " --- Expected ---"   << Exp_Z_signal_region << endl;

}

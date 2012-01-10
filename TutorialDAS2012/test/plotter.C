
{

  TFile fData("treeSkimmedMuTau_Data.root");
  TFile fDYJets("treeSkimmedMuTau_DYJets.root");
  TFile fWJets("treeSkimmedMuTau_WJets.root");
  TFile fTTJets("treeSkimmedMuTau_TTJets.root");

  TTree* tData   = (TTree*)fData.Get("outTreePtOrd");
  TTree* tDYJets = (TTree*)fDYJets.Get("outTreePtOrd");
  TTree* tWJets  = (TTree*)fWJets.Get("outTreePtOrd");
  TTree* tTTJets = (TTree*)fTTJets.Get("outTreePtOrd");

  int nBins = 10;
  TArrayF bins(nBins+1);
  bins[0] = 0;    bins[1] = 30;    bins[2] = 40;    bins[3] = 50;    bins[4] = 60;    bins[5] = 70;
  bins[6] = 80;   bins[7] = 100;   bins[8] = 120;   bins[9] = 150;   bins[10] = 200;


  TH1F* hData    = new TH1F("hData","",   nBins, bins.GetArray());
  TH1F* hQCD     = new TH1F("hQCD","",    nBins, bins.GetArray());
  TH1F* hDYJets  = new TH1F("hDYJets","", nBins, bins.GetArray());
  TH1F* hWJets   = new TH1F("hWJets","",  nBins, bins.GetArray());
  TH1F* hWJetsSS = new TH1F("hWJetsSS","",nBins, bins.GetArray());
  TH1F* hTTJets  = new TH1F("hTTJets","", nBins, bins.GetArray());

  hData->SetMarkerStyle(kFullCircle);
  hQCD->SetFillColor(kMagenta-9);
  hDYJets->SetFillColor(kYellow-9);
  hWJets->SetFillColor(kRed-3);
  hTTJets->SetFillColor(kBlue-2);

  TCut sbin(       "ptL1>17 && ptL2>20 && tightestHPSDBWP>0 && diTauCharge==0 && MtLeg1<40 && muFlag==0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10");
  TCut sbinAntiW(  "ptL1>17 && ptL2>20 && tightestHPSDBWP>0 && diTauCharge==0 && MtLeg1>60 && muFlag==0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10");
  TCut sbinSS(     "ptL1>17 && ptL2>20 && tightestHPSDBWP>0 && diTauCharge!=0 && MtLeg1<40 && muFlag==0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10");
  TCut sbinSSAntiW("ptL1>17 && ptL2>20 && tightestHPSDBWP>0 && diTauCharge!=0 && MtLeg1>60 && muFlag==0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10");

  float lumiFact    = 53./100;
  float OStoSSRatio = 1.07;
  // estimation of W+jets
  TH1F* h1 = new TH1F("h1","",1,-10,10);
  tWJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*SFMu*SFTau*sampleWeight"*sbinAntiW);
  float WsbinAntiW  = h1->Integral()*lumiFact;
  h1->Reset();
  tWJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*SFMu*SFTau*sampleWeight"*sbin);
  float Wsbin       = h1->Integral()*lumiFact;
  h1->Reset();
  tTTJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*SFMu*SFTau*sampleWeight"*sbinAntiW);
  float TTsbinAntiW = h1->Integral()*lumiFact;
  h1->Reset();
  tData->Draw("etaL1>>h1",   sbinAntiW);
  float DatasbinAntiW = h1->Integral();
  h1->Reset();

  cout << "Wsbin (MC) = " << Wsbin << endl;
  float Wsbin = (DatasbinAntiW - TTsbinAntiW)*(Wsbin/WsbinAntiW);
  cout << "Wsbin = (" << DatasbinAntiW << " - " << TTsbinAntiW << " )*" << Wsbin/WsbinAntiW << " = " << Wsbin << endl;

  // estimation of QCD
  TH1F* h1 = new TH1F("h1","",1,-10,10);
  tWJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*SFMu*SFTau*sampleWeight"*sbinSSAntiW);
  float WsbinSSAntiW  = h1->Integral()*lumiFact;
  h1->Reset();
  tWJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*SFMu*SFTau*sampleWeight"*sbinSS);
  float WsbinSS       = h1->Integral()*lumiFact;
  h1->Reset();
  tTTJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*SFMu*SFTau*sampleWeight"*sbinSSAntiW);
  float TTsbinSSAntiW = h1->Integral()*lumiFact;
  h1->Reset();
  tData->Draw("etaL1>>h1",   sbinSSAntiW);
  float DatasbinAntiW = h1->Integral();
  h1->Reset();
  tData->Draw("etaL1>>h1",   sbinSS);
  float DatasbinSS = h1->Integral();
  h1->Reset();

  float WsbinSS   = (DatasbinAntiW - TTsbinAntiW)*(Wsbin/WsbinAntiW);
  float QCDsbinSS = DatasbinSS - WsbinSS;
  float QCDsbin   = QCDsbinSS*OStoSSRatio;

  tData->Draw("diTauVisMass>>hData", sbin);
  tData->Draw("diTauVisMass>>hQCD",  sbinSS);
  tDYJets->Draw("diTauVisMass>>hDYJets", "puWeight*SFMu*SFTau*sampleWeight"*sbin);
  tWJets->Draw("diTauVisMass>>hWJets",   "puWeight*HLTweightMu*SFMu*SFTau*sampleWeight"*sbin);
  tWJets->Draw("diTauVisMass>>hWJetsSS",   "puWeight*HLTweightMu*SFMu*SFTau*sampleWeight"*sbinSS);
  tTTJets->Draw("diTauVisMass>>hTTJets", "puWeight*HLTweightMu*SFMu*SFTau*sampleWeight"*sbin);

  hDYJets->Scale( lumiFact );
  hWJets->Scale(     Wsbin/hWJets->Integral());
  hWJetsSS->Scale( WsbinSS/hWJetsSS->Integral());
  hTTJets->Scale( lumiFact );
  hQCD->Add(   hWJetsSS, -1);
  hQCD->Scale( QCDsbin/hQCD->Integral());

  THStack* aStack = new THStack("aStack","");
  aStack->Add(hTTJets);
  aStack->Add(hQCD);
  aStack->Add(hWJets);
  aStack->Add(hDYJets);

  hData->Sumw2();
  hData->Draw("P");
  aStack->Draw("HISTSAME");
  hData->Draw("PSAME");
}

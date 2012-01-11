
{

  TFile fData("treeSkimmedMuTau_Data.root");
  TFile fDYJets("treeSkimmedMuTau_DYJets.root");
  TFile fWJets("treeSkimmedMuTau_WJets.root");
  TFile fTTJets("treeSkimmedMuTau_TTJets.root");
  TFile fVBFH130("treeSkimmedMuTau_VBFH130.root");
  TFile fGGFH130("treeSkimmedMuTau_GGFH130.root");

  TTree* tData   = (TTree*)fData.Get("outTreePtOrd");
  TTree* tDYJets = (TTree*)fDYJets.Get("outTreePtOrd");
  TTree* tWJets  = (TTree*)fWJets.Get("outTreePtOrd");
  TTree* tTTJets = (TTree*)fTTJets.Get("outTreePtOrd");
  TTree* tVBFH130= (TTree*)fVBFH130.Get("outTreePtOrd");
  TTree* tGGFH130= (TTree*)fGGFH130.Get("outTreePtOrd");

  int nBins = 11;
  TArrayF bins(nBins+1);
  bins[0] = 0;    bins[1] = 30;    bins[2] = 40;    bins[3] = 50;    bins[4]  = 60;    bins[5] = 70;
  bins[6] = 80;   bins[7] = 90;    bins[8] = 100;   bins[9] = 120;   bins[10] = 150;  bins[11] = 200;
  //for(int k = 0 ; k <= nBins ; k++) bins[k] = 25*k; 
  TString variable("diTauVisMass");
  TString labels(" ; mass (GeV) ; Events");

  TH1F* hData    = new TH1F("hData"   ,labels, nBins, bins.GetArray());
  TH1F* hQCD     = new TH1F("hQCD"    ,labels, nBins, bins.GetArray());
  TH1F* hDYJets  = new TH1F("hDYJets" ,labels, nBins, bins.GetArray());
  TH1F* hWJets   = new TH1F("hWJets"  ,labels, nBins, bins.GetArray());
  TH1F* hWJetsSS = new TH1F("hWJetsSS",labels, nBins, bins.GetArray());
  TH1F* hTTJets  = new TH1F("hTTJets" ,labels, nBins, bins.GetArray());
  TH1F* hVBFH130 = new TH1F("hVBFH130",labels, nBins, bins.GetArray());
  TH1F* hGGFH130 = new TH1F("hGGFH130",labels, nBins, bins.GetArray());

  hData->SetMarkerStyle(kFullCircle);
  hQCD->SetFillColor(kMagenta-9);
  hDYJets->SetFillColor(kYellow-9);
  hWJets->SetFillColor(kRed-3);
  hTTJets->SetFillColor(kBlue-2);
  hVBFH130->SetLineWidth(3);
  hVBFH130->SetLineStyle(kDashed);
  hVBFH130->SetLineColor(kBlue);
  hGGFH130->SetLineWidth(3);
  hGGFH130->SetLineStyle(kDashed);
  hGGFH130->SetLineColor(kBlue);

  TCut sbin(         "ptL1>17 && ptL2>20 && tightestHPSDBWP>0 && diTauCharge==0 && MtLeg1<40 && muFlag==0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10 && muFlag==0");
  TCut sbinAntiW(    "ptL1>17 && ptL2>20 && tightestHPSDBWP>0 && diTauCharge==0 && MtLeg1>60 && muFlag==0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10 && muFlag==0");
  TCut sbinSS(       "ptL1>17 && ptL2>20 && tightestHPSDBWP>0 && diTauCharge!=0 && MtLeg1<40 && muFlag==0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10 && muFlag==0");
  TCut sbinSSRelIso( "ptL1>17 && ptL2>20 && tightestHPSDBWP>0 && diTauCharge!=0 && MtLeg1<40 && muFlag==0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.30 && muFlag==0");
  TCut sbinSSAntiW(  "ptL1>17 && ptL2>20 && tightestHPSDBWP>0 && diTauCharge!=0 && MtLeg1>60 && muFlag==0 && HLTx && HLTmatch && combRelIsoLeg1DBeta<0.10 && muFlag==0");

  float lumiFact    = 53./100;
  float OStoSSRatio = 1.07;
  // estimation of W+jets
  TH1F* h1 = new TH1F("h1","",1,-10,10);
  tWJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbinAntiW);
  float WsbinAntiW  = h1->Integral()*lumiFact;
  h1->Reset();
  tWJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbin);
  float Wsbin       = h1->Integral()*lumiFact;
  h1->Reset();
  tTTJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbinAntiW);
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
  tWJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbinSSAntiW);
  float WsbinSSAntiW  = h1->Integral()*lumiFact;
  h1->Reset();
  tWJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbinSS);
  float WsbinSS       = h1->Integral()*lumiFact;
  h1->Reset();
  tTTJets->Draw("etaL1>>h1",   "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbinSSAntiW);
  float TTsbinSSAntiW = h1->Integral()*lumiFact;
  h1->Reset();
  tData->Draw("etaL1>>h1",   sbinSSAntiW);
  float DatasbinSSAntiW = h1->Integral();
  h1->Reset();
  tData->Draw("etaL1>>h1",   sbinSS);
  float DatasbinSS = h1->Integral();
  h1->Reset();
  
  cout << "WsbinSS (MC) = " << WsbinSS << endl; 
  float WsbinSS   = (DatasbinSSAntiW - TTsbinSSAntiW)*(WsbinSS/WsbinSSAntiW);
  cout << "WsbinSS = (" << DatasbinSSAntiW << " - " << TTsbinSSAntiW << " )*" << WsbinSS/WsbinSSAntiW << " = " << WsbinSS << endl;
  float QCDsbinSS = DatasbinSS - WsbinSS;
  float QCDsbin   = QCDsbinSS*OStoSSRatio;

  tData->Draw(  variable+">>hData",   sbin);
  tData->Draw(  variable+">>hQCD",    sbinSSRelIso);
  tDYJets->Draw(variable+">>hDYJets", "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbin);
  tWJets->Draw( variable+">>hWJets",  "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbin);
  tWJets->Draw( variable+">>hWJetsSS","puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbinSS);
  tTTJets->Draw(variable+">>hTTJets", "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbin);
  tVBFH130->Draw(variable+">>hVBFH130", "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbin);
  tGGFH130->Draw(variable+">>hGGFH130", "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight*HqTWeight"*sbin);
  

  hDYJets->Scale( lumiFact );
  hWJets->Scale(     Wsbin/hWJets->Integral());
  hWJetsSS->Scale( WsbinSS/hWJetsSS->Integral());
  hTTJets->Scale( lumiFact );
  hQCD->Add(   hWJetsSS, -1);
  hQCD->Scale( QCDsbin/hQCD->Integral());

  hVBFH130->Scale( lumiFact*100 );
  hGGFH130->Scale( lumiFact*100 );
  hVBFH130->Add(hGGFH130,1.0);

  THStack* aStack = new THStack("aStack","");
  aStack->Add(hTTJets);
  aStack->Add(hQCD);
  aStack->Add(hWJets);
  aStack->Add(hDYJets);
  aStack->Add(hVBFH130);

  hData->Sumw2();
  hData->Draw("P");
  aStack->Draw("HISTSAME");
  hData->Draw("PSAME");
  //hVBFH130->Draw("HISTSAME");

  TLegend* leg = new TLegend(0.52,0.50,0.75,0.87,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  leg->SetHeader("CMS Preliminary 2011 #sqrt{s}=7 TeV, L=53 pb^{-1}");
  leg->AddEntry(hData,"Observed","P");
  leg->AddEntry(hDYJets,"Z#rightarrow#tau#tau","F");
  leg->AddEntry(hWJets,"W+jets","F");
  leg->AddEntry(hQCD,"QCD","F");
  leg->AddEntry(hTTJets,"t#bar{t}","F");
  leg->AddEntry(hVBFH130,"100 X H#rightarrow#tau#tau, m_{H}=130 GeV","F");
  leg->Draw();

}


{

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetStatColor(0);
  gStyle->SetOptStat(0000000);
  gStyle->SetOptFit(0111);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetPalette(1);

  TFile fData("treeSkimmedMuTau_Data.root");
  TFile fDYJets("treeSkimmedMuTau_DYJets.root");
  TFile fWJets("treeSkimmedMuTau_WJets.root");
  TFile fTTJets("treeSkimmedMuTau_TTJets.root");
  TFile fVBFH130("treeSkimmedMuTau_VBFH130.root");
  TFile fGGFH130("treeSkimmedMuTau_GGFH130.root");
  TFile fQCD("treeSkimmedMuTau_QCD.root");

  TTree* tData   = (TTree*)fData.Get("outTreePtOrd");
  TTree* tDYJets = (TTree*)fDYJets.Get("outTreePtOrd");
  TTree* tWJets  = (TTree*)fWJets.Get("outTreePtOrd");
  TTree* tTTJets = (TTree*)fTTJets.Get("outTreePtOrd");
  TTree* tVBFH130= (TTree*)fVBFH130.Get("outTreePtOrd");
  TTree* tGGFH130= (TTree*)fGGFH130.Get("outTreePtOrd");
  TTree* tQCD    = (TTree*)fQCD.Get("outTreePtOrd");

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

  // full signal selection
  TCut sbin(         "ptL1>15 && ptL2>20 && tightestHPSDBWP>0 && diTauCharge==0&& HLTx && HLTmatch");
  // luminosity of data is 53 pb
  float lumiFact    = 53./100;

  // Draw with cuts and weights !!!
  tData->Draw(  variable+">>hData",   sbin);
  tQCD->Draw(  variable+">>hQCD",     "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbin);
  tDYJets->Draw(variable+">>hDYJets", "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbin);
  tWJets->Draw( variable+">>hWJets",  "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbin);
  tTTJets->Draw(variable+">>hTTJets", "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbin);
  tVBFH130->Draw(variable+">>hVBFH130", "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight"*sbin);
  tGGFH130->Draw(variable+">>hGGFH130", "puWeight*HLTweightMu*HLTweightTau*SFMu*SFTau*sampleWeight*HqTWeight"*sbin);
  
  // Scale histograms
  hDYJets->Scale( lumiFact );
  hWJets->Scale( lumiFact );
  hTTJets->Scale( lumiFact );
  hQCD->Scale( lumiFact );
  hVBFH130->Scale( lumiFact*100 );
  hGGFH130->Scale( lumiFact*100 );
  hVBFH130->Add(hGGFH130,1.0);

  // Add all together
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



  // Legend
  TLegend* leg = new TLegend(0.52,0.50,0.75,0.87,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  leg->SetHeader("#splitline{CMS Preliminary 2011}{#sqrt{s}=7 TeV, L=53 pb^{-1}}");
  leg->AddEntry(hData,"Observed","P");
  leg->AddEntry(hDYJets,"Z#rightarrow#tau#tau","F");
  leg->AddEntry(hWJets,"W+jets","F");
  leg->AddEntry(hQCD,"QCD","F");
  leg->AddEntry(hTTJets,"t#bar{t}","F");
  leg->AddEntry(hVBFH130,"100 X H#rightarrow#tau#tau, m_{H}=130 GeV","F");
  leg->Draw();

}

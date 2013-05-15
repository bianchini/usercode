#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>

#include "TSystem.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"
#include "TVectorD.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"


#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooExponential.h"
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
#include "RooGaussian.h"
#include "RooAbsCategory.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "RooTruthModel.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooBernstein.h"
#include "RooNDKeysPdf.h"
#include "RooChiSquarePdf.h"

#include "Bianchi/TTHStudies/interface/MEIntegratorNew.h"
#include "Bianchi/TTHStudies/interface/Samples.h"


#define GENJETDR 0.3
#define VERBOSE  false

using namespace std;
using namespace RooFit;

typedef struct 
{
  float et; 
  float sumet;
  float sig;
  float phi;
} metInfo;


typedef struct
{
  float bmass;
  float bpt;
  float beta;
  float bphi;
  float bstatus;
  float wdau1mass;
  float wdau1pt;
  float wdau1eta;
  float wdau1phi;
  float wdau1id;
  float wdau2mass;
  float wdau2pt;
  float wdau2eta;
  float wdau2phi;
  float wdau2id;
} genTopInfo;

typedef struct 
{
  float mass; 
  float pt;
  float eta;
  float phi;
  float status;
  float charge;
  float momid;
} genParticleInfo;


int main(int argc, const char* argv[])
{

  std::cout << "TestMENew" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");
  //gSystem->Load("BianchiTTHStudiesPlugins");

  AutoLibraryLoader::enable();

  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;
  std::string outFileName  ( in.getParameter<std::string>  ("outFileName" ) );
  std::string pathToFile   ( in.getParameter<std::string>("pathToFile" ) );
  std::string ordering     ( in.getParameter<std::string>("ordering" ) );
  std::string pathToTF     ( in.getParameter<std::string> ("pathToTF"   ) );
  int vegasPoints          ( in.getParameter<int> ("vegasPoints"   ) );
  bool verbose             ( in.getParameter<bool>("verbose") );
  double lumi              ( in.getParameter<double>("lumi") );

  float met       ( in.getParameter<double> ("met") );
  float pertW1    ( in.getParameter<double> ("pertW1") );
  float pertW2    ( in.getParameter<double> ("pertW2") );
  float pertBHad  ( in.getParameter<double> ("pertBHad") );
  float pertBLep  ( in.getParameter<double> ("pertBLep") );
  float enlargeE1 ( in.getParameter<double>  ("enlargeE1") );
  float enlargeEh1( in.getParameter<double> ("enlargeEh1") );
  float enlargePt ( in.getParameter<double>  ("enlargePt") );
  vector<int> evLimits( in.getParameter<vector<int> >  ("evLimits") );
  int evLow = evLimits[0];
  int evHigh = evLimits[1];

  //TFile* fout = 0;
  //fout = new TFile(outFileName.c_str(),"RECREATE");
  //fout->Close();
  //cout << "Deleting file..." << endl;
  //delete fout
  gSystem->Exec(("rm "+outFileName).c_str());


  int printP4 = 0;
  int par = 4;
  MEIntegratorNew* meIntegrator = new MEIntegratorNew( pathToTF , par, int(verbose));
  //meIntegrator->setFile(fout);

  ROOT::Math::GSLMCIntegrator ig2("vegas", 1.e-12, 1.e-5, vegasPoints);
  ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par); 
  ig2.SetFunction(toIntegrate);
  
  const int nMassPoints = 21;
  double mH[nMassPoints] = {140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240};
  TH1F*  hMass     = new TH1F("hMass","",21, 137.5, 242.5);
  TH1F*  hBestMass = new TH1F("hBestMass","",100, 50, 250);

  ////////////////////////////////////////////////////////


  bool openAllFiles  = false;
  Samples* mySamples = new Samples(openAllFiles, pathToFile, ordering, samples, lumi, verbose);
  vector<string> mySampleFiles;

  if(mySamples->IsOk()){

    cout << "Ok!" << endl;
    mySampleFiles = mySamples->Files();

    for( unsigned int i = 0 ; i < mySampleFiles.size(); i++){
      string sampleName       = mySampleFiles[i];

      if(verbose){
	cout << mySampleFiles[i] << " ==> " << mySamples->GetXSec(sampleName) 
	     << " pb,"
	     << " ==> weight = "            << mySamples->GetWeight(sampleName) << endl;
      }
    }
  }
  else{
    cout << "Problems... leaving" << endl;
    return 0;
  }


  ////////////////////////////////////////////////////////


  for(unsigned int i = 0 ; i < mySampleFiles.size(); i++){
    
    string currentName       = mySampleFiles[i];

    mySamples->OpenFile( currentName );
    cout << "Opening file " << currentName << endl;
    TTree* currentTree       = mySamples->GetTree( currentName, "tree");
    cout << "Done!!" << endl;
    
    genParticleInfo genB, genBbar;
    genTopInfo genTop, genTbar;
    metInfo METtype1p2corr;
 
    currentTree->SetBranchAddress("genB",   &genB);
    currentTree->SetBranchAddress("genBbar",&genBbar);
    currentTree->SetBranchAddress("genTop", &genTop);
    currentTree->SetBranchAddress("genTbar",&genTbar);
    currentTree->SetBranchAddress("METtype1p2corr",&METtype1p2corr);
 
    int counter = 0;
    Long64_t nentries = currentTree->GetEntries();
    for (Long64_t i = 0; i < nentries ; i++){
      
      if(i%5000==0) cout << i << endl;
      currentTree->GetEntry(i);
      
      //continue;
      
      TLorentzVector topBLV   ;
      TLorentzVector topW1LV  ; 
      TLorentzVector topW2LV  ;
      TLorentzVector atopBLV  ; 
      TLorentzVector atopW1LV ; 
      TLorentzVector atopW2LV ;
      TLorentzVector genBLV   ;
      TLorentzVector genBbarLV;
      if(genTop.bmass>0){
	topBLV.SetPtEtaPhiM(  genTop.bpt,genTop.beta,genTop.bphi,genTop.bmass );
	topW1LV.SetPtEtaPhiM( genTop.wdau1pt,genTop.wdau1eta, genTop.wdau1phi,genTop.wdau1mass);
	topW2LV.SetPtEtaPhiM( genTop.wdau2pt,genTop.wdau2eta, genTop.wdau2phi,genTop.wdau2mass);
      }
      if(genTbar.bmass>0){
	atopBLV.SetPtEtaPhiM(  genTbar.bpt,genTbar.beta,genTbar.bphi,genTbar.bmass );
	atopW1LV.SetPtEtaPhiM( genTbar.wdau1pt,genTbar.wdau1eta, genTbar.wdau1phi,genTbar.wdau1mass);
	atopW2LV.SetPtEtaPhiM( genTbar.wdau2pt,genTbar.wdau2eta,genTbar.wdau2phi,genTbar.wdau2mass);
      }
      if(genB.mass>0 && genB.momid==25){
	genBLV.SetPtEtaPhiM(genB.pt,genB.eta ,genB.phi, genB.mass );
      }
      if(genBbar.mass>0 && genBbar.momid==25){
	genBbarLV.SetPtEtaPhiM(genBbar.pt,genBbar.eta ,genBbar.phi, genBbar.mass );
      }
  
      TLorentzVector TOPLEP;
      TLorentzVector TOPHAD;
      TLorentzVector TOPHADW1;
      TLorentzVector TOPHADW2;
      TLorentzVector TOPHADB;
      TLorentzVector TOPLEPW1;
      TLorentzVector TOPLEPW2;
      TLorentzVector TOPLEPB;
      TLorentzVector HIGGS;

      float neutCosTheta;

      bool properEvent = (genBLV.Pt()>0 && genBbarLV.Pt()>0 && topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0);

      HIGGS.SetPxPyPzE( (genBLV+genBbarLV).Px(), (genBLV+genBbarLV).Py(),(genBLV+genBbarLV).Pz(),(genBLV+genBbarLV).E());
      if( abs(genTop.wdau1id)>6 && abs(genTbar.wdau1id)<6){
	TOPLEP.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	TOPHAD.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	TOPHADW1.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	TOPHADW2.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	TOPHADB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	if( abs(genTop.wdau1id)==11 || abs(genTop.wdau1id)==13 || abs(genTop.wdau1id)==15 ){
	  TOPLEPW1.SetPxPyPzE  ( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	  TOPLEPW2.SetPtEtaPhiM( topW2LV.Pt(), 0.0, topW2LV.Phi(), 0.0);
	  neutCosTheta = TMath::Cos( topW2LV.Vect().Theta() );
	}
	else{
	  TOPLEPW2.SetPtEtaPhiM( topW1LV.Pt(), 0.0, topW1LV.Phi(), 0.0);
	  TOPLEPW1.SetPxPyPzE  ( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	  neutCosTheta = TMath::Cos( topW1LV.Vect().Theta() );
	}
	TOPLEPB.SetPxPyPzE(  topBLV.Px(),  topBLV.Py(),   topBLV.Pz(), topBLV.E());
      }
      else if(abs(genTop.wdau1id)<6 && abs(genTbar.wdau1id)>6){
	TOPHAD.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	TOPLEP.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	TOPHADW1.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	TOPHADW2.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	TOPHADB.SetPxPyPzE( topBLV.Px(),  topBLV.Py(),   topBLV.Pz(),  topBLV.E());
	if( abs(genTbar.wdau1id)==11 || abs(genTbar.wdau1id)==13 || abs(genTbar.wdau1id)==15 ){
	  TOPLEPW1.SetPxPyPzE  ( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	  TOPLEPW2.SetPtEtaPhiM( atopW2LV.Pt(), 0.0, atopW2LV.Phi(), 0.0);
	  neutCosTheta = TMath::Cos( atopW2LV.Vect().Theta() );
	}
	else{
	  TOPLEPW2.SetPtEtaPhiM( atopW1LV.Pt(), 0.0, atopW1LV.Phi(), 0.0);
	  TOPLEPW1.SetPxPyPzE  ( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	  neutCosTheta = TMath::Cos( atopW1LV.Vect().Theta() );
	}
	TOPLEPB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
      }      
      else{
	properEvent=false;
      }

      properEvent = ( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
		      TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
		      TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
		      TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
		      TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
		      genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
		      genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
		      );

      if(!properEvent) continue;
      counter++;

      if( !(counter>=evLow && counter<=evHigh) ) continue;
      cout << "Processing event # " << counter << endl;

      //TLorentzVector jet1,jet2,jet3,jet4,jet5,jet6,jet7,jet8;
      //jet1.SetPtEtaPhiM(19.8454,            0.822914,-2.3529,0.10566);     // lep
      //jet2.SetPtEtaPhiM(91.934,             0.0,-0.0161266,   0.);         // MET    0.543608 <=> cosTheta = 0.495714;
      //jet3.SetPtEtaPhiM(175.197 * pertBLep, 0.0936224,-0.969521,4.8);      // b from top lep
      //jet4.SetPtEtaPhiM(122.914 * pertW1,   1.17562,2.17839,0.33);         // j1 from W
      //jet5.SetPtEtaPhiM(34.9939 * pertW2,   1.06994,0.879248,0.33);        // j2 from W
      //jet6.SetPtEtaPhiM(95.5899 * pertBHad, 0.829954,3.10669,4.8);         // b from top hadr
      //jet7.SetPtEtaPhiM(32.805,             0.678347,-0.0951163,4.8);      // b1 from H
      //jet8.SetPtEtaPhiM(75.4309,            2.08021,2.75503,4.8);          // b2 from H

      vector<TLorentzVector> jets;
      jets.push_back( TOPLEPW1 );  
      jets.push_back( TOPLEPW2 );  
      jets.push_back( TOPLEPB );
      jets.push_back( TOPHADW1 );  
      jets.push_back( TOPHADW2 );  
      jets.push_back( TOPHADB );
      jets.push_back( genBLV );
      jets.push_back( genBbarLV );

      if(printP4){
	cout << "*******START******" << endl;      
	cout << "SumEt = " << METtype1p2corr.sumet  << endl;
	cout << "lep:  jet1.SetPtEtaPhiM(" << TOPLEPW1.Pt() << "," <<  TOPLEPW1.Eta() << "," << TOPLEPW1.Phi() << "," << TOPLEPW1.M() << ")" << endl;
	cout << "met:  jet2.SetPtEtaPhiM(" << TOPLEPW2.Pt() << "," <<  TOPLEPW2.Eta() << "," << TOPLEPW2.Phi() << "," << TOPLEPW2.M() << ")" << endl;
	cout << "blep: jet3.SetPtEtaPhiM(" << TOPLEPB.Pt() << "," <<  TOPLEPB.Eta() << "," << TOPLEPB.Phi() << "," << TOPLEPB.M() << ")" << endl;
	cout << "w1:   jet4.SetPtEtaPhiM(" << TOPHADW1.Pt() << "," <<  TOPHADW1.Eta() << "," << TOPHADW1.Phi() << "," << TOPHADW1.M() << ")" << endl;
	cout << "w2:   jet5.SetPtEtaPhiM(" << TOPHADW2.Pt() << "," <<  TOPHADW2.Eta() << "," << TOPHADW2.Phi() << "," << TOPHADW2.M() << ")" << endl;
	cout << "bhad: jet6.SetPtEtaPhiM(" << TOPHADB.Pt() << "," <<  TOPHADB.Eta() << "," << TOPHADB.Phi() << "," << TOPHADB.M() << ")" << endl;
	cout << "h1:   jet7.SetPtEtaPhiM(" << genBLV.Pt() << "," <<  genBLV.Eta() << "," << genBLV.Phi() << "," << genBLV.M() << ")" << endl;
	cout << "h2:   jet8.SetPtEtaPhiM(" << genBbarLV.Pt() << "," <<  genBbarLV.Eta() << "," << genBbarLV.Phi() << "," << genBbarLV.M() << ")" << endl;	
	cout << "Top Lep mass = " << (TOPLEPW1+TOPLEPW2+TOPLEPB).M() << " <=> neutrino eta=0!!!" << endl;
	cout << "Top Had mass = " << (TOPHADW1+TOPHADW2+TOPHADB).M() << endl;
	cout << "Higgs mass = "   << (genBLV+genBbarLV).M() << endl;
	cout << "Neut cosTheta = " << neutCosTheta << endl;
      }


      /////////////////////////////////////////////////////////////////////////////

      if(printP4) cout << ">>>" << endl;
      meIntegrator->SetPar(par);
      if(printP4) cout << "Setup jets..." << endl;
      meIntegrator->setJets(&jets);
      if(printP4) cout << "Setup H mass..." << endl;
      meIntegrator->setMass( met );
      meIntegrator->setSumEt( METtype1p2corr.sumet );
      if(printP4) cout << "Setup versors..." << endl;
      meIntegrator->initVersors(1);
      if(printP4) cout << "Setup TFs..." << endl;
      meIntegrator->initTF();
      if(printP4) cout << ">>>" << endl;
      
      /////////////////////////////////////////////////////////////////////////////
    
      double E1low   = (meIntegrator->getCachedTF("tfWjet1"))->GetXaxis()->GetXmin() * (1-enlargeE1);
      double E1high  = (meIntegrator->getCachedTF("tfWjet1"))->GetXaxis()->GetXmax() * (1+enlargeE1);
      
      double Ptlow   = (meIntegrator->getCachedTF("tfMetPt"))->GetXaxis()->GetXmin() * (1-enlargePt);
      double Pthigh  = (meIntegrator->getCachedTF("tfMetPt"))->GetXaxis()->GetXmax() * (1+enlargePt);
      
      Ptlow  = -1; //0.781794 - 0.05;
      Pthigh = +1; //0.781794 + 0.05;
      
      double Philow   = -(meIntegrator->getCachedTF("tfMetPhi"))->GetYaxis()->GetXmax();
      double Phihigh  =  (meIntegrator->getCachedTF("tfMetPhi"))->GetYaxis()->GetXmax();
            
      //Philow  = -0.04;
      //Phihigh = +0.04;

      double Eh1low   = (meIntegrator->getCachedTF("tfHiggs1"))->GetXaxis()->GetXmin() * (1-enlargeEh1);
      double Eh1high  = (meIntegrator->getCachedTF("tfHiggs1"))->GetXaxis()->GetXmax() * (1+enlargeEh1);
      
      /////////////////////////////////////////////////////////////////////////////
      
      if(printP4){
	cout << "E1 :       ["  << E1low  << "," << E1high  << "]" << endl;
	cout << "cosTheta : ["  << Ptlow  << "," << Pthigh  << "]" << endl;
	cout << "Phi :      ["  << Philow << "," << Phihigh << "]" << endl;
	cout << "Eh1 :      ["  << Eh1low << "," << Eh1high << "]" << endl;
      }

      /////////////////////////////////////////////////////////////////////////////
      
      double xL[4] = {  E1low,  Ptlow,   Philow,   Eh1low};
      double xU[4] = {  E1high, Pthigh , Phihigh,  Eh1high};
      
      hMass->Reset();
      for(int m = 0; m < nMassPoints ; m++){
	//meIntegrator->setMass( mH[m] );
	meIntegrator->setTopMass( mH[m], 80.19);
	double p = ig2.Integral(xL, xU);
	 if(printP4) cout << "Mass " << mH[m] << " => prob  = " << p << endl;
	hMass->Fill(mH[m], p);
      }

      hBestMass->Fill( hMass->GetBinCenter(hMass->GetMaximumBin())  );

       if(printP4) cout << "Opening file to save histos" << endl;
      TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");
      TDirectory *dir = fout_tmp->mkdir(Form("Event_%d", counter));
      dir->cd();
      (meIntegrator->getCachedPdf("pdfGammaWHad"))->Write("pdfGammaWHad",TObject::kOverwrite);
      (meIntegrator->getCachedPdf("pdfBetaWHad")) ->Write("pdfBetaWHad",TObject::kOverwrite);
      (meIntegrator->getCachedPdf("pdfGammaWLep"))->Write("pdfGammaWLep",TObject::kOverwrite);
      (meIntegrator->getCachedPdf("pdfBetaWLep")) ->Write("pdfBetaWLep",TObject::kOverwrite);
      (meIntegrator->getCachedPdf("pdfGammaTTH")) ->Write("pdfGammaTTH",TObject::kOverwrite);
      (meIntegrator->getCachedPdf("pdf3D"))       ->Write("pdf3D",TObject::kOverwrite);
      (meIntegrator->getCachedTF("tfWjet1"))      ->Write("tfWjet1",TObject::kOverwrite);
      (meIntegrator->getCachedTF("tfWjet2"))      ->Write("tfWjet2",TObject::kOverwrite);  
      (meIntegrator->getCachedTF("tfbHad"))       ->Write("tfbHad",TObject::kOverwrite);
      (meIntegrator->getCachedTF("tfbLep"))       ->Write("tfbLep",TObject::kOverwrite);
      (meIntegrator->getCachedTF("tfHiggs1"))     ->Write("tfHiggs1",TObject::kOverwrite);
      (meIntegrator->getCachedTF("tfHiggs2"))     ->Write("tfHiggs2",TObject::kOverwrite);
      (meIntegrator->getCachedTF("tfMetPt"))      ->Write("tfMetPt",TObject::kOverwrite);
      (meIntegrator->getCachedTF("tfMetPhi"))     ->Write("tfMetPhi",TObject::kOverwrite);
      hMass->Write("hMass",TObject::kOverwrite);
       if(printP4) cout << "Closing file to save histos" << endl;
      fout_tmp->Close();


      meIntegrator->deleteTF();
      if(printP4) cout << "*******END********" << endl;      

    }
    
  }
  

  TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");
  hBestMass->Write("hBestMass",TObject::kOverwrite);
  fout_tmp->Close();

  //cout << "Closing file..." << endl;
  //fout->Close();
  //cout << "Deleting file..." << endl;
  //delete fout;
  
  cout << "Delete meIntegrator..." << endl;
  delete meIntegrator;
  cout << "Finished!!!" << endl;
  return 0;
  

}

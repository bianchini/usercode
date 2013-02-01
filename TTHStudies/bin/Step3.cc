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
#include "TList.h"
#include "THStack.h"
#include "TCut.h"
#include "TArrayF.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TKey.h"
#include "TMultiGraph.h"
#include "Bianchi/TTHStudies/interface/Samples.h"
#include "Bianchi/TTHStudies/interface/Test.h"
#include "Math/GenVector/LorentzVector.h"

#include <algorithm>

#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"

#define VERBOSE 1

using namespace std;

typedef  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LV;

struct sorterByPt {
  bool operator() (float i,float j) const { return (i>j);}
};


int main(int argc, const char* argv[])
{

  std::cout << "Step3" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();

  PythonProcessDesc builder(argv[1]);
 
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;

  std::string pathToFile( in.getParameter<std::string>("pathToFile" ) );
  std::string ordering(   in.getParameter<std::string>("ordering" ) );
  double lumi(            in.getParameter<double>("lumi" ) );

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


  for(unsigned int i = 0 ; i < mySampleFiles.size(); i++){
    
    string currentName       = mySampleFiles[i];
    
    //fwlite::TFileService fs = fwlite::TFileService( ("TThNTuples_"+currentName+".root").c_str());
    //TTree* outTree = fs.make<TTree>("outTree","TTH Tree");
    //TIter nextkey((mySamples->GetFile( currentName ))->GetListOfKeys());
    //TH1F *key;
    //while ( (key = (TH1F*)nextkey()) ) {
    //string name(key->GetName());
    //if( name.find("tree")!=string::npos) continue;
    //TH1F* h = (TH1F*)(mySamples->GetFile( currentName ))->FindObjectAny( key->GetName());
    //h->Write(key->GetName());
    //}
    //outTree = (mySamples->GetTree( currentName, "tree"))->CopyTree("","");


    cout << "Start copying..." << endl;
    TTree* outTree = mySamples->GetTree( currentName, "tree");
    cout << "Done!!" << endl;

    //////////////////////////////////NEW VARIABLES///////////////////////////////
    //////////////////////////////////////////////////////////////////////////////

    int hJetRank_[100];
    int aJetRank_[100];
    float pt1_,  pt2_,  pt3_,  pt4_,  pt5_,  pt6_;
    float eta1_, eta2_, eta3_, eta4_, eta5_, eta6_;
    float phi1_, phi2_, phi3_, phi4_, phi5_, phi6_;
    float csv1_, csv2_, csv3_, csv4_, csv5_, csv6_;
    int numJets30_; int numJets30bTag_;
    int numJets20_; int numJets20bTag_;
    float isotropy_, circularity_, sphericity_, aplanarity_, Cparam_, Dparam_;
    float aveCsv_;
    float aveDeltaRbTag_;
    float aveMbTag_, aveMunTag_;
    float closestJJbTagMass_;
    float varCsv_;
    float maxCsv_;  float minCsv_;
    float sumAllPt_;
    float massAll_;
    float massLJ_;
    float M3_;

    TBranch *hJetRankBR = outTree->Branch("hJetRank",hJetRank_,"jhJetRank[nhJets]/I");
    TBranch *aJetRankBR = outTree->Branch("aJetRank",aJetRank_,"ahJetRank[naJets]/I");

    TBranch *pt1BR = outTree->Branch("pt1",&pt1_,"pt1/F");
    TBranch *pt2BR = outTree->Branch("pt2",&pt2_,"pt2/F");
    TBranch *pt3BR = outTree->Branch("pt3",&pt3_,"pt3/F");
    TBranch *pt4BR = outTree->Branch("pt4",&pt4_,"pt4/F");
    TBranch *pt5BR = outTree->Branch("pt5",&pt5_,"pt5/F");
    TBranch *pt6BR = outTree->Branch("pt6",&pt6_,"pt6/F");

    TBranch *eta1BR = outTree->Branch("eta1",&eta1_,"eta1/F");
    TBranch *eta2BR = outTree->Branch("eta2",&eta2_,"eta2/F");
    TBranch *eta3BR = outTree->Branch("eta3",&eta3_,"eta3/F");
    TBranch *eta4BR = outTree->Branch("eta4",&eta4_,"eta4/F");
    TBranch *eta5BR = outTree->Branch("eta5",&eta5_,"eta5/F");
    TBranch *eta6BR = outTree->Branch("eta6",&eta6_,"eta6/F");
 
    TBranch *phi1BR = outTree->Branch("phi1",&phi1_,"phi1/F");
    TBranch *phi2BR = outTree->Branch("phi2",&phi2_,"phi2/F");
    TBranch *phi3BR = outTree->Branch("phi3",&phi3_,"phi3/F");
    TBranch *phi4BR = outTree->Branch("phi4",&phi4_,"phi4/F");
    TBranch *phi5BR = outTree->Branch("phi5",&phi5_,"phi5/F");
    TBranch *phi6BR = outTree->Branch("phi6",&phi6_,"phi6/F");

    TBranch *csv1BR = outTree->Branch("csv1",&csv1_,"csv1/F");
    TBranch *csv2BR = outTree->Branch("csv2",&csv2_,"csv2/F");
    TBranch *csv3BR = outTree->Branch("csv3",&csv3_,"csv3/F");
    TBranch *csv4BR = outTree->Branch("csv4",&csv4_,"csv4/F");
    TBranch *csv5BR = outTree->Branch("csv5",&csv5_,"csv5/F");
    TBranch *csv6BR = outTree->Branch("csv6",&csv6_,"csv6/F");

    TBranch *numJets30BR     = outTree->Branch("numJets30",&numJets30_,"numJets30/I");
    TBranch *numJets20BR     = outTree->Branch("numJets20",&numJets20_,"numJets20/I");
    TBranch *numJets30bTagBR = outTree->Branch("numJets30bTag",&numJets30bTag_,"numJets30bTag/I");
    TBranch *numJets20bTagBR = outTree->Branch("numJets20bTag",&numJets20bTag_,"numJets20bTag/I");

    TBranch *isotropyBR      = outTree->Branch("isotropy",&isotropy_,"isotropy/F");
    TBranch *circularityBR   = outTree->Branch("circularity",&circularity_,"circularity/F");
    TBranch *sphericityBR    = outTree->Branch("sphericity",&sphericity_,"sphericity/F");
    TBranch *aplanarityBR    = outTree->Branch("aplanarity",&aplanarity_,"aplanarity/F");
    TBranch *CparamBR        = outTree->Branch("Cparam",&Cparam_,"Cparam/F");
    TBranch *DparamBR        = outTree->Branch("Dparam",&Dparam_,"Dparam/F");

    TBranch *aveCsvBR            = outTree->Branch("aveCsv",&aveCsv_,"aveCsv/F");
    TBranch *aveDeltaRbTagBR     = outTree->Branch("aveDeltaRbTag",&aveDeltaRbTag_,"aveDeltaRbTag/F");
    TBranch *aveMbTagBR          = outTree->Branch("aveMbTag",&aveMbTag_,"aveMbTag/F");
    TBranch *aveMunTagBR         = outTree->Branch("aveMunTag",&aveMunTag_,"aveMunTag/F");
    TBranch *closestJJbTagMassBR = outTree->Branch("closestJJbTagMass",&closestJJbTagMass_,"closestJJbTagMass/F");

    TBranch *varCsvBR            = outTree->Branch("varCsv",&varCsv_,"varCsv/F");
    TBranch *maxCsvBR            = outTree->Branch("maxCsv",&maxCsv_,"maxCsv/F");
    TBranch *minCsvBR            = outTree->Branch("minCsv",&minCsv_,"minCsv/F");

    TBranch *sumAllPtBR          = outTree->Branch("sumAllPt",&sumAllPt_,"sumAllPt/F");

    TBranch *massAllBR           = outTree->Branch("massAll",&massAll_,"massAll/F");
    TBranch *massLJBR            = outTree->Branch("massLJ",&massLJ_,"massLJ/F");
    TBranch *M3BR                = outTree->Branch("M3",&M3_,"M3/F");
    //TBranch *BR = outTree->Branch("",&_,"/F");

    //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////





    //////////////////////////////////INPUT VARIABLES/////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
   
    int nhJets;
    int naJets;
    int nvlep;
    float met[999];
    float hJetspt[999]; float hJetseta[999]; float hJetsphi[999]; float hJetse[999]; float hJetscsv[999]; 
    float aJetspt[999]; float aJetseta[999]; float aJetsphi[999]; float aJetse[999]; float aJetscsv[999]; 
    float vLeptonpt[999]; float vLeptoneta[999]; float vLeptonphi[999];

    outTree->SetBranchAddress("nhJets",   &nhJets);
    outTree->SetBranchAddress("naJets",   &naJets);

    outTree->SetBranchAddress("hJet_pt",   hJetspt);
    outTree->SetBranchAddress("hJet_eta",  hJetseta);
    outTree->SetBranchAddress("hJet_phi",  hJetsphi);
    outTree->SetBranchAddress("hJet_e",    hJetse);
    outTree->SetBranchAddress("hJet_csv",  hJetscsv);

    outTree->SetBranchAddress("aJet_pt",   aJetspt);
    outTree->SetBranchAddress("aJet_eta",  aJetseta);
    outTree->SetBranchAddress("aJet_phi",  aJetsphi);
    outTree->SetBranchAddress("aJet_e",    aJetse);
    outTree->SetBranchAddress("aJet_csv",  aJetscsv);

    outTree->SetBranchAddress("vLepton_pt", vLeptonpt);
    outTree->SetBranchAddress("vLepton_phi",vLeptonphi);
    outTree->SetBranchAddress("vLepton_eta",vLeptoneta);
    outTree->SetBranchAddress("nvlep",      &nvlep);

    outTree->SetBranchAddress("METtype1p2corr",  met);
    //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////



    Long64_t nentries = outTree->GetEntries();
 
    for (Long64_t i = 0; i < nentries; i++){

      if(i%5000==0) cout << i << endl;

      pt1_  = -99; pt2_  = -99; pt3_  = -99; pt4_  = -99; pt5_  = -99; pt6_  = -99;
      eta1_ = -99; eta2_ = -99; eta3_ = -99; eta4_ = -99; eta5_ = -99; eta6_ = -99;
      phi1_ = -99; phi2_ = -99; phi3_ = -99; phi4_ = -99; phi5_ = -99; phi6_ = -99;
      csv1_ = -99; csv2_ = -99; csv3_ = -99; csv4_ = -99; csv5_ = -99; csv6_ = -99;
      sumAllPt_ = 0.; massLJ_ = 0.; varCsv_ = 0.;

      outTree->GetEntry(i);

      std::vector<math::RhoEtaPhiVector> allJetsForShapeVar;

      std::vector<LV> allLeptons;
      std::vector<LV> allBtagJetsForShapeVar;
      std::vector<LV> allUntagJetsForShapeVar;
      std::vector<float> csvBtag;

      std::map<float, int, sorterByPt> hMapPt;
      std::map<float, int, sorterByPt> aMapPt;
      std::map<float, int, sorterByPt> allMapPt30;

      numJets20_     = 0;  numJets30_     = 0;
      numJets30bTag_ = 0;  numJets20bTag_ = 0;
      float sumBtag  = 0.;


      for(int i = 0; i < nhJets; i++){
	hMapPt[ hJetspt[i]]   = i;

	if( hJetspt[i] > 20.){
	  numJets20_++;
	  if( hJetscsv[i] > 0.898) numJets20bTag_++;
	}
	if( hJetspt[i]>30.){
	  numJets30_++;
	  float jetMass2 = hJetse[i]*hJetse[i] -  TMath::Power(hJetspt[i]*TMath::CosH(hJetseta[i]) ,2);
	  if( hJetscsv[i] > 0.898){
	    numJets30bTag_++;
	    sumBtag += hJetscsv[i];
	    allBtagJetsForShapeVar.push_back( LV(hJetspt[i], hJetseta[i], hJetsphi[i], TMath::Sqrt(jetMass2)) );
	    csvBtag.push_back( hJetscsv[i] );
	  }
	  else{
	    allUntagJetsForShapeVar.push_back( LV(hJetspt[i], hJetseta[i], hJetsphi[i], TMath::Sqrt(jetMass2)) );
	  }
	  allMapPt30[ hJetspt[i]] = i;

	  allJetsForShapeVar.push_back( math::RhoEtaPhiVector(hJetspt[i], hJetseta[i], hJetsphi[i]) );

	}
      }


      for(int i = 0; i < naJets; i++){
	aMapPt[ aJetspt[i]] = i;

	if( aJetspt[i] > 20.){
	  numJets20_++;
	  if( aJetscsv[i] > 0.898) numJets20bTag_++;
	}
	if( aJetspt[i]>30.){
	  numJets30_++;
	  float jetMass2 = aJetse[i]*aJetse[i] -  TMath::Power(aJetspt[i]*TMath::CosH(aJetseta[i]),2);
	  if( aJetscsv[i] > 0.898){
	    numJets30bTag_++;
	    sumBtag += aJetscsv[i];
	    allBtagJetsForShapeVar.push_back( LV(aJetspt[i], aJetseta[i], aJetsphi[i],  TMath::Sqrt(jetMass2)) );
	    csvBtag.push_back( aJetscsv[i] );
	  }
	  else{
	    allUntagJetsForShapeVar.push_back( LV(aJetspt[i], aJetseta[i], aJetsphi[i],  TMath::Sqrt(jetMass2)) );
	  }
	  allMapPt30[ aJetspt[i]] = -i-1; // distinguish jet collection

	  allJetsForShapeVar.push_back( math::RhoEtaPhiVector(aJetspt[i], aJetseta[i], aJetsphi[i]) );

	}
      }



      aveCsv_ = numJets30bTag_>0 ? sumBtag/numJets30bTag_ : -99;

      float sumDeltaRbTag = 0.;
      float sumMassBtag   = 0.;
      float sumMassUntag  = 0.;
      float minDeltaRBtag = 999.;

      if(allBtagJetsForShapeVar.size()>1){
	for(unsigned k = 0; k< allBtagJetsForShapeVar.size()-1 ; k++){
	  for(unsigned l = k+1; l<allBtagJetsForShapeVar.size() ; l++){
	    LV first  = allBtagJetsForShapeVar[k];
	    LV second = allBtagJetsForShapeVar[l];
	    float deltaR = TMath::Sqrt( TMath::Power(first.Eta()-second.Eta(),2) + TMath::Power(first.Phi()-second.Phi(),2)  );
	    sumDeltaRbTag += deltaR;
	    if(deltaR < minDeltaRBtag){
	      minDeltaRBtag = deltaR;
	      closestJJbTagMass_ = (first+second).M();
	    }
	    sumMassBtag   += (first+second).M();
	  }
	}
      }

      aveDeltaRbTag_ = numJets30bTag_>1 ? sumDeltaRbTag/(numJets30bTag_/2*(numJets30bTag_-1)) : -99;
      aveMbTag_      = numJets30bTag_>1 ? sumMassBtag/(numJets30bTag_/2*(numJets30bTag_-1))   : -99;


      float maxCsv = -999.; float minCsv = 999.;
      for(unsigned k = 0; k< csvBtag.size() ; k++){
	varCsv_ +=  TMath::Power(csvBtag[k]-aveCsv_,2);
	if( csvBtag[k]>maxCsv ) maxCsv = csvBtag[k];
	if( csvBtag[k]<minCsv ) minCsv = csvBtag[k];
      }
      maxCsv_ = maxCsv;
      minCsv_ = minCsv;
      varCsv_ = numJets30bTag_>1 ? varCsv_/(numJets30bTag_-1) : -99;
  
      if(allUntagJetsForShapeVar.size()>1){
	for(unsigned k = 0; k< allUntagJetsForShapeVar.size()-1 ; k++){
	  for(unsigned l = k+1; l<allUntagJetsForShapeVar.size() ; l++){
	    LV first  = allUntagJetsForShapeVar[k];
	    LV second = allUntagJetsForShapeVar[l];
	    sumMassUntag   += (first+second).M();
	  }
	}
      }
      aveMunTag_     = (numJets30_-numJets30bTag_)>1 ? sumMassUntag/((numJets30_-numJets30bTag_)/2*((numJets30_-numJets30bTag_)-1)) : -99;
      

      LV allP4(0,0,0,0);
      LV MEtP4(met[0], 0., met[3] , 0.);      
    
      for(int k = 0; k<nvlep ; k++){
	sumAllPt_ += vLeptonpt[k];
	allLeptons.push_back( LV(vLeptonpt[k], vLeptoneta[k], vLeptonphi[k], 0.) );
      }


      if(allBtagJetsForShapeVar.size()>1){
	float minDeltaR = 999.;
	for(unsigned k = 0; k< allBtagJetsForShapeVar.size() ; k++){
	  for(unsigned l = 0; l<allLeptons.size() ; l++){
	    LV first  = allBtagJetsForShapeVar[k];
	    LV second = allLeptons[l];
	    float deltaR = TMath::Sqrt( TMath::Power(first.Eta()-second.Eta(),2) + TMath::Power(first.Phi()-second.Phi(),2)  );
	    if(deltaR < minDeltaR){
	      massLJ_ = (first+second).M();
	      minDeltaR = deltaR;
	    }
	  }
	}
      }

      if(allBtagJetsForShapeVar.size()>0 && allUntagJetsForShapeVar.size()>1){
	float highestPt = -999;
	for(unsigned k = 0; k< allBtagJetsForShapeVar.size() ; k++){
	  LV zero = allBtagJetsForShapeVar[k];
	  for(unsigned l = 0; l<allUntagJetsForShapeVar.size()-1 ; l++){
	    for(unsigned m = l+1; m<allUntagJetsForShapeVar.size() ; m++){
	      LV first  = allUntagJetsForShapeVar[l];
	      LV second = allUntagJetsForShapeVar[m];
	      if( (zero+first+second).Pt()> highestPt){
		highestPt = (zero+first+second).Pt();
		M3_ = (zero+first+second).M();
	      }
	    }
	  }

	}
      }


      for(unsigned k = 0; k< allLeptons.size() ; k++){
	allP4 += allLeptons[k];
      }
      for(unsigned k = 0; k< allUntagJetsForShapeVar.size() ; k++){
	sumAllPt_ += allUntagJetsForShapeVar[k].Pt();
	allP4     += allUntagJetsForShapeVar[k];
      }
      for(unsigned k = 0; k< allBtagJetsForShapeVar.size() ; k++){
	sumAllPt_ += allBtagJetsForShapeVar[k].Pt();
	allP4     += allBtagJetsForShapeVar[k];
      }

      sumAllPt_ += MEtP4.Pt();
      allP4     += MEtP4;

      massAll_ = allP4.M();


      int h = 0; int a = 0;
      for(std::map<float, int>::iterator it = hMapPt.begin(); it!=hMapPt.end(); it++, h++) hJetRank_[h] = it->second ;
      for(std::map<float, int>::iterator it = aMapPt.begin(); it!=aMapPt.end(); it++, a++) aJetRank_[a] = it->second ;

      int all = 0;
      for(std::map<float, int>::iterator it = allMapPt30.begin(); it!=allMapPt30.end(); it++, all++){


	int index = it->second;

	if(all==0){
	  pt1_  = (index >= 0 ) ?  hJetspt[index]  : aJetspt[ -index-1];
	  eta1_ = (index >= 0 ) ?  hJetseta[index] : aJetseta[-index-1];
	  phi1_ = (index >= 0 ) ?  hJetsphi[index] : aJetsphi[-index-1];
	  csv1_ = (index >= 0 ) ?  hJetscsv[index] : aJetscsv[-index-1];
	}
	else if(all==1){
	  pt2_  = (index >= 0 ) ?  hJetspt[index]  : aJetspt[ -index-1];
	  eta2_ = (index >= 0 ) ?  hJetseta[index] : aJetseta[-index-1];
	  phi2_ = (index >= 0 ) ?  hJetsphi[index] : aJetsphi[-index-1];
	  csv2_ = (index >= 0 ) ?  hJetscsv[index] : aJetscsv[-index-1];
	}
	else if(all==2){
	  pt3_  = (index >= 0 ) ?  hJetspt[index]  : aJetspt[ -index-1];
	  eta3_ = (index >= 0 ) ?  hJetseta[index] : aJetseta[-index-1];
	  phi3_ = (index >= 0 ) ?  hJetsphi[index] : aJetsphi[-index-1];
	  csv3_ = (index >= 0 ) ?  hJetscsv[index] : aJetscsv[-index-1];
	}
	else if(all==3){
	  pt4_  = (index >= 0 ) ?  hJetspt[index]  : aJetspt[ -index-1];
	  eta4_ = (index >= 0 ) ?  hJetseta[index] : aJetseta[-index-1];
	  phi4_ = (index >= 0 ) ?  hJetsphi[index] : aJetsphi[-index-1];
	  csv4_ = (index >= 0 ) ?  hJetscsv[index] : aJetscsv[-index-1];
	}
	else if(all==4){
	  pt5_  = (index >= 0 ) ?  hJetspt[index]  : aJetspt[ -index-1];
	  eta5_ = (index >= 0 ) ?  hJetseta[index] : aJetseta[-index-1];
	  phi5_ = (index >= 0 ) ?  hJetsphi[index] : aJetsphi[-index-1];
	  csv5_ = (index >= 0 ) ?  hJetscsv[index] : aJetscsv[-index-1];
	}
	else if(all==5){
	  pt6_  = (index >= 0 ) ?  hJetspt[index]  : aJetspt[ -index-1];
	  eta6_ = (index >= 0 ) ?  hJetseta[index] : aJetseta[-index-1];
	  phi6_ = (index >= 0 ) ?  hJetsphi[index] : aJetsphi[-index-1];
	  csv6_ = (index >= 0 ) ?  hJetscsv[index] : aJetscsv[-index-1];
	}
	else{}


      }

      EventShapeVariables* eventShapes = new EventShapeVariables(allJetsForShapeVar);
      isotropy_   = eventShapes->isotropy();
      circularity_= eventShapes->circularity();
      sphericity_ = eventShapes->sphericity();
      aplanarity_ = eventShapes->aplanarity();
      Cparam_     = eventShapes->C();
      Dparam_     = eventShapes->D();


      //////////////////////////////////FILL////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////

      hJetRankBR->Fill();
      aJetRankBR->Fill();
      numJets30BR->Fill();
      numJets20BR->Fill();
      numJets30bTagBR->Fill();
      numJets20bTagBR->Fill();
      pt1BR->Fill();
      pt2BR->Fill();
      pt3BR->Fill();
      pt4BR->Fill();
      pt5BR->Fill();
      pt6BR->Fill();
      eta1BR->Fill();
      eta2BR->Fill();
      eta3BR->Fill();
      eta4BR->Fill();
      eta5BR->Fill();
      eta6BR->Fill();
      phi1BR->Fill();
      phi2BR->Fill();
      phi3BR->Fill();
      phi4BR->Fill();
      phi5BR->Fill();
      phi6BR->Fill();
      csv1BR->Fill();
      csv2BR->Fill();
      csv3BR->Fill();
      csv4BR->Fill();
      csv5BR->Fill();
      csv6BR->Fill();
      isotropyBR->Fill();
      circularityBR->Fill();
      sphericityBR->Fill();
      aplanarityBR->Fill();
      CparamBR->Fill();
      DparamBR->Fill();
      aveCsvBR->Fill();
      aveDeltaRbTagBR->Fill();
      aveMbTagBR->Fill();
      aveMunTagBR->Fill();
      closestJJbTagMassBR->Fill();
      varCsvBR->Fill();
      maxCsvBR->Fill();
      minCsvBR->Fill();
      sumAllPtBR->Fill();
      massAllBR->Fill();
      massLJBR->Fill();
      M3BR->Fill();


      delete eventShapes;

    }


    outTree->Write("",TObject::kOverwrite );

   }



  return 0;


}

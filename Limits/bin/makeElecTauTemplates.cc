#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TSystem.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>


#define RESCALETO1PB true



using namespace std;



float higgsXsection(int mH = 120, string process = "ggH"){

  float xsection = 1.0;

  if(process.find("ggH")!=string::npos){
    if(mH == 110) xsection = 8.02e-02*19.84;
    if(mH == 115) xsection = 7.65e-02*18.13;
    if(mH == 120) xsection = 7.10e-02*16.63;
    if(mH == 125) xsection = 6.37e-02*15.31;
    if(mH == 130) xsection = 5.48e-02*14.12;
    if(mH == 135) xsection = 4.52e-02*13.08;
    if(mH == 140) xsection = 3.54e-02*12.13;
    if(mH == 145) xsection = 2.61e-02*11.27;
  }
  if(process.find("qqH")!=string::npos){
    if(mH == 110) xsection = 8.02e-02*1.398;
    if(mH == 115) xsection = 7.65e-02*1.332;
    if(mH == 120) xsection = 7.10e-02*1.269;
    if(mH == 125) xsection = 6.37e-02*1.211;
    if(mH == 130) xsection = 5.48e-02*1.154;
    if(mH == 135) xsection = 4.52e-02*1.100;
    if(mH == 140) xsection = 3.54e-02*1.052;
    if(mH == 145) xsection = 2.61e-02*1.004;
  }
  if(process.find("VH")!=string::npos){
    if(mH == 110) xsection = 8.02e-02*(0.8754+0.4721+0.1257  );
    if(mH == 115) xsection = 7.65e-02*(0.7546+0.4107+0.1106  );
    if(mH == 120) xsection = 7.10e-02*(0.6561+0.3598+0.09756 );
    if(mH == 125) xsection = 6.37e-02*(0.5729+0.3158+0.08634 );
    if(mH == 130) xsection = 5.48e-02*(0.4942+0.2778+0.07658 );
    if(mH == 135) xsection = 4.52e-02*(0.4390+0.2453+0.06810 );
    if(mH == 140) xsection = 3.54e-02*(0.3857+0.2172+0.06072 );
    if(mH == 145) xsection = 2.61e-02*(0.3406+0.1930+0.05435 );
  }

  return xsection;
}



void produce(   
	     int mH_=120,
	     string variable_  = "diTauVisMass",
	     string analysis_  = "",
	     string bin_       = "inclusive",
	     TString outputDir = "Dec2011"
	     ){


  cout << "Now doing mass mH=" << mH_ << ", for variable " << variable_ << " analysis " << analysis_ << " and bin " << bin_ << endl;
  TFile* fin = new TFile(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/histograms/%s/eTau_mH%d_%s_%s_%s.root", outputDir.Data(), 120, bin_.c_str() , analysis_.c_str(), variable_.c_str()), "READ");


  ///////////////////////////////////////////////
  TFile* fin_jUp     = new TFile(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/histograms/%s/eTau_mH%d_%s_JetUp_%s.root",   outputDir.Data(), 120, bin_.c_str() , variable_.c_str()), "READ");
  TFile* fin_jDown   = new TFile(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/histograms/%s/eTau_mH%d_%s_JetDown_%s.root", outputDir.Data(), 120, bin_.c_str() , variable_.c_str()), "READ");
  TFile* fin_tUp     = new TFile(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/histograms/%s/eTau_mH%d_%s_TauUp_%s.root",   outputDir.Data(), 120, bin_.c_str() , variable_.c_str()), "READ");
  TFile* fin_tDown   = new TFile(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/histograms/%s/eTau_mH%d_%s_TauDown_%s.root", outputDir.Data(), 120, bin_.c_str() , variable_.c_str()), "READ");
  TFile* fin_nominal = new TFile(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/histograms/%s/eTau_mH%d_%s__%s.root",        outputDir.Data(), 120, bin_.c_str() , variable_.c_str()), "READ");
  ///////////////////////////////////////////////

  float rescaleggH = RESCALETO1PB ? higgsXsection(mH_,"ggH") : 1.0;
  float rescaleqqH = RESCALETO1PB ? higgsXsection(mH_,"qqH") : 1.0;
  float rescaleVH  = RESCALETO1PB ? higgsXsection(mH_,"VH")  : 1.0;


  string binNameSpace = "";
  if(bin_.find("inclusive")!=string::npos)      
    binNameSpace =  "inclusive";

  else if(bin_.find("novbf")!=string::npos) 
    binNameSpace =  "SM0";

  if(bin_.find("novbfLow")  !=string::npos) 
    binNameSpace =  "0jet_low";

  if(bin_.find("novbfHigh") !=string::npos) 
    binNameSpace =  "0jet_high";


  else if(bin_.find("boost")!=string::npos) 
    binNameSpace =  "SM1";

  if(bin_.find("boostLow")!=string::npos) 
    binNameSpace =  "boost_low";

  if(bin_.find("boostHigh")!=string::npos) 
    binNameSpace =  "boost_high";


  else if(bin_.find("vbf")!=string::npos && bin_.find("novbf")==string::npos) 
    binNameSpace =  "vbf";

  else if(bin_.find("vh")!=string::npos) 
    binNameSpace =  "2jet";


  else if(bin_.find("bTag")!=string::npos) 
    binNameSpace =  "SM4";

  if(bin_.find("bTagLow")!=string::npos) 
    binNameSpace =  "btag_low";

  if(bin_.find("bTagHigh")!=string::npos) 
    binNameSpace =  "btag_high";


  else if(bin_.find("twoJets")!=string::npos) 
    binNameSpace =  "SMpre2";

  else if(bin_.find("oneJet")!=string::npos) 
    binNameSpace =  "SMpre2a";

  //TFile* fTemplOut = new TFile(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/datacards/%s/eTauSM_%s.root",outputDir.Data(), variable_.c_str()),"UPDATE");
  TFile* fTemplOut = new TFile(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/datacards/%s/eTauSM.root",outputDir.Data()),"UPDATE");
  
  string suffix = "";
  if(analysis_.find("TauUp")!=string::npos)
    suffix = "_CMS_scale_tUp";
  else if(analysis_.find("TauDown")!=string::npos)
    suffix = "_CMS_scale_tDown";
  else if(analysis_.find("MuUp")!=string::npos)
    suffix = "_CMS_scale_mUp";
  else if(analysis_.find("MuDown")!=string::npos)
    suffix = "_CMS_scale_mDown";
  else if(analysis_.find("JetUp")!=string::npos)
    suffix = "_CMS_scale_jUp";
  else if(analysis_.find("JetDown")!=string::npos)
    suffix = "_CMS_scale_jDown";
  
  cout << "Adding histos with suffix " << suffix << endl;
  TString dirName( Form("eTau_%s",binNameSpace.c_str()) );
  
  if(! (fTemplOut->cd( dirName.Data()  ) ) ){
    
    cout << "Editing the directory for bin and variable " << binNameSpace << ", " << variable_ << endl;

    TDirectory* dir = fTemplOut->mkdir( dirName.Data() );
    dir->cd();


    ((TH1F*)fin->Get("hData"))->Write("data_obs");

    TH1F* hSgn2 = (TH1F*)fin->Get(Form("hggH%d",mH_));
    hSgn2->Scale(1./rescaleggH);
    hSgn2->Write(Form("ggH%d%s" ,mH_,suffix.c_str()));

    TH1F* hSgn1 = (TH1F*)fin->Get(Form("hqqH%d",mH_));
    hSgn1->Scale(1./rescaleqqH);
    hSgn1->Write(Form("qqH%d%s" ,mH_,suffix.c_str()));

    TH1F* hSgn3 = (TH1F*)fin->Get(Form("hVH%d",mH_));
    hSgn3->Scale(1./rescaleVH);
    hSgn3->Write(Form("VH%d%s" ,mH_,suffix.c_str()));
				  

    if(bin_.find("novbf")!=string::npos){

      ((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s",suffix.c_str()));

      TH1F* hQCD         = (TH1F*)fin->Get("hQCD");
      TH1F* hAntiIsoKeys = (TH1F*)fin->Get("hAntiIsoKeys") ;
      //hAntiIsoKeys->Scale(hQCD->Integral()/hAntiIsoKeys->Integral());

      if( bin_.find("Low")!=string::npos ){
	//hQCD->Write(Form("QCD%s"    ,suffix.c_str()));
	hAntiIsoKeys->Scale(hQCD->Integral()/hAntiIsoKeys->Integral());
	hAntiIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));
      }
      else{
	//hAntiIsoKeys->Scale(hQCD->Integral()/hAntiIsoKeys->Integral());
	hAntiIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));
      }
   
      //((TH1F*)fin->Get("hW"))->Write(Form("W%s"           ,suffix.c_str()));
      //((TH1F*)fin->Get("hWKeys"))->Write(Form("W%s"         ,suffix.c_str()));
      TH1F* hW     = (TH1F*)fin->Get("hW");
      TH1F* hWKeys = (TH1F*)fin->Get("hWKeys");
      hWKeys->Scale(hW->Integral()/hWKeys->Integral());
      hWKeys->Write(Form("W%s"           ,suffix.c_str()));


      ((TH1F*)fin->Get("hZmj"))->Write(Form("ZJ%s"        ,suffix.c_str()));
      //((TH1F*)fin->Get("hZmm"))->Write(Form("ZL%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hZmmKeys"))->Write(Form("ZL%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hZfakes"))->Write(Form("ZLL%s"    ,suffix.c_str()));
      ((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hVV"))->Write(Form("VV%s"         ,suffix.c_str()));

    }
    else if(bin_.find("boost")!=string::npos){

    
      ((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s",suffix.c_str()));

      TH1F* hQCD         = (TH1F*)fin->Get("hQCD");
      TH1F* hAntiIsoKeys = (TH1F*)fin->Get("hAntiIsoKeys");
      //hAntiIsoKeys->Scale(hQCD->Integral()/hAntiIsoKeys->Integral());
      hAntiIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));
      //hQCD->Write(Form("QCD%s"    ,suffix.c_str()));

      //TH1F* hAntiIsoKeys = (TH1F*)fin->Get("hAntiIsoKeys");
      //hAntiIsoKeys->Scale(hQCD->Integral()/hAntiIsoKeys->Integral());
      //hAntiIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));


      ((TH1F*)fin->Get("hWKeys"))->Write(Form("W%s"           ,suffix.c_str()));
      //((TH1F*)fin->Get("hW"))->Write(Form("W%s"           ,suffix.c_str()));

      //TH1F* hZmjKeys    = (TH1F*)fin->Get("hZmjKeys");
      //TH1F* hZmmKeys    = (TH1F*)fin->Get("hZmmKeys");
      TH1F* hZmjKeys    = (TH1F*)fin->Get("hZmj");
      TH1F* hZmmKeys    = (TH1F*)fin->Get("hZmm");
      TH1F* hZfakesKeys = (TH1F*)fin->Get("hZfakes");
      hZfakesKeys->Reset();
      hZfakesKeys->Add(hZmjKeys,1.0);
      hZfakesKeys->Add(hZmmKeys,1.0);
      hZfakesKeys->Write(Form("ZLL%s"        ,suffix.c_str()));

      ((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hVV"))->Write(Form("VV%s"         ,suffix.c_str()));
    }

    else if(bin_.find("bTag")!=string::npos){

      ((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s",suffix.c_str()));

      TH1F* hQCD         = (TH1F*)fin->Get("hQCD");
      //hQCD->Write(Form("QCD%s"    ,suffix.c_str()));
      TH1F* hAntiIsoKeys = (TH1F*)fin->Get("hAntiIsoKeys");
      //TH1F* hAntiIso      = (TH1F*)fin->Get("hAntiIso");
      //TH1F* hLooseIsoKeys = (TH1F*)fin->Get("hLooseIsoKeys");
      //hLooseIsoKeys->Scale(hAntiIso->Integral()/hLooseIsoKeys->Integral());
      //hLooseIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));

      if( bin_.find("Low")!=string::npos )
	hQCD->Write(Form("QCD%s"    ,suffix.c_str()));
	//hAntiIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));
      else
	hAntiIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));

      ((TH1F*)fin->Get("hWKeys"))->Write(Form("W%s"           ,suffix.c_str()));
      //((TH1F*)fin->Get("hW"))->Write(Form("W%s"           ,suffix.c_str()));

      //TH1F* hZmjKeys    = (TH1F*)fin->Get("hZmjKeys");
      //TH1F* hZmmKeys    = (TH1F*)fin->Get("hZmmKeys");
      TH1F* hZmjKeys    = (TH1F*)fin->Get("hZmj");
      TH1F* hZmmKeys    = (TH1F*)fin->Get("hZmm");
      TH1F* hZfakesKeys = (TH1F*)fin->Get("hZfakes");
      hZfakesKeys->Reset();
      hZfakesKeys->Add(hZmjKeys,1.0);
      hZfakesKeys->Add(hZmmKeys,1.0);
      hZfakesKeys->Write(Form("ZLL%s"        ,suffix.c_str()));

      ((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hVV"))->Write(Form("VV%s"         ,suffix.c_str()));
    }
    else if(bin_.find("vbf")!=string::npos && bin_.find("novbf")==string::npos){

      ((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s",suffix.c_str()));

      //((TH1F*)fin->Get("hAntiIsoKeys"))->Write(Form("QCD%s"    ,suffix.c_str()));
      TH1F* hAntiIsoKeys  = (TH1F*)fin->Get("hAntiIsoKeys");
      TH1F* hAntiIso      = (TH1F*)fin->Get("hAntiIso");
      TH1F* hLooseIsoKeys = (TH1F*)fin->Get("hLooseIsoKeys");
      hLooseIsoKeys->Scale(hAntiIso->Integral()/hLooseIsoKeys->Integral());
      hLooseIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));
      //((TH1F*)fin->Get("hW3JetsKeys"))->Write(Form("W%s"           ,suffix.c_str()));

      //((TH1F*)fin->Get("hWKeys"))->Write(Form("W%s"           ,suffix.c_str()));
      ((TH1F*)fin->Get("hW3Jets"))->Write(Form("W%s"           ,suffix.c_str()));

      TH1F* hZmjKeys    = (TH1F*)fin->Get("hZmjKeys");
      TH1F* hZmmKeys    = (TH1F*)fin->Get("hZmmKeys");
      TH1F* hZmj        = (TH1F*)fin->Get("hZmj");
      TH1F* hZmm        = (TH1F*)fin->Get("hZmm");
      TH1F* hZfakesKeys = (TH1F*)hZmjKeys->Clone("hZfakesKeys");
      hZfakesKeys->Reset();
      if(hZmjKeys->Integral()>0) hZfakesKeys->Add(hZmjKeys,1.0);
      else hZfakesKeys->Add(hZmj,1.0);
      if(hZmmKeys->Integral()>0) hZfakesKeys->Add(hZmmKeys,1.0);  
      else hZfakesKeys->Add(hZmm,1.0);

      hZfakesKeys->Write(Form("ZLL%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hVVKeys"))->Write(Form("VV%s"         ,suffix.c_str()));
    }
    else if(bin_.find("vh")!=string::npos){

      ((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s",suffix.c_str()));

      //((TH1F*)fin->Get("hAntiIsoKeys"))->Write(Form("QCD%s"    ,suffix.c_str()));
      TH1F* hAntiIsoKeys  = (TH1F*)fin->Get("hAntiIsoKeys");
      TH1F* hAntiIso      = (TH1F*)fin->Get("hAntiIso");
      TH1F* hLooseIsoKeys = (TH1F*)fin->Get("hLooseIsoKeys");
      hLooseIsoKeys->Scale(hAntiIso->Integral()/hLooseIsoKeys->Integral());
      hLooseIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));
      //((TH1F*)fin->Get("hW3JetsKeys"))->Write(Form("W%s"           ,suffix.c_str()));

      ((TH1F*)fin->Get("hWKeys"))->Write(Form("W%s"           ,suffix.c_str()));

      TH1F* hZmjKeys    = (TH1F*)fin->Get("hZmjKeys");
      TH1F* hZmmKeys    = (TH1F*)fin->Get("hZmmKeys");
      TH1F* hZmj        = (TH1F*)fin->Get("hZmj");
      TH1F* hZmm        = (TH1F*)fin->Get("hZmm");
      TH1F* hZfakesKeys = (TH1F*)hZmjKeys->Clone("hZfakesKeys");
      hZfakesKeys->Reset();
      if(hZmjKeys->Integral()>0) hZfakesKeys->Add(hZmjKeys,1.0);
      else hZfakesKeys->Add(hZmj,1.0);
      if(hZmmKeys->Integral()>0) hZfakesKeys->Add(hZmmKeys,1.0);  
      else hZfakesKeys->Add(hZmm,1.0);

      hZfakesKeys->Write(Form("ZLL%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hVVKeys"))->Write(Form("VV%s"         ,suffix.c_str()));
    }
    else{


      ((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s",suffix.c_str()));
      ((TH1F*)fin->Get("hQcd"))->Write(Form("QCD%s"    ,suffix.c_str()));
      ((TH1F*)fin->Get("hW"))->Write(Form("W%s"           ,suffix.c_str()));
      ((TH1F*)fin->Get("hZmj"))->Write(Form("ZJ%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hZmm"))->Write(Form("ZL%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hZfakes"))->Write(Form("ZLL%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hVV"))->Write(Form("VV%s"         ,suffix.c_str()));
    }


           
  }else{
    cout << "Directory is there, filling new histos..." << endl;

    TDirectory* dir = (TDirectory*)fTemplOut->Get( dirName.Data() );

    fTemplOut->cd( dirName.Data() );
    

    if(dir->FindObjectAny(Form("qqH%d%s"         ,mH_,suffix.c_str()))==0 ){
      TH1F* hSgn2 = (TH1F*)fin->Get(Form("hqqH%d",mH_));
      hSgn2->Scale(1./rescaleqqH);
      hSgn2->Write(Form("qqH%d%s" ,mH_,suffix.c_str()));
    }
    if(dir->FindObjectAny(Form("ggH%d%s"          , mH_,suffix.c_str()))==0 ){
      TH1F* hSgn1 = (TH1F*)fin->Get(Form("hggH%d",mH_));
      hSgn1->Scale(1./rescaleggH);
      hSgn1->Write(Form("ggH%d%s" ,mH_,suffix.c_str()));
    }
    if(dir->FindObjectAny(Form("VH%d%s"          , mH_,suffix.c_str()))==0 ){
      TH1F* hSgn3 = (TH1F*)fin->Get(Form("hVH%d",mH_));
      hSgn3->Scale(1./rescaleVH);
      hSgn3->Write(Form("VH%d%s" ,mH_,suffix.c_str()));
    }


    if(bin_.find("novbf")!=string::npos){

      if(dir->FindObjectAny(Form("ZTT%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s"       ,suffix.c_str()));
      if(dir->FindObjectAny(Form("QCD%s"       ,suffix.c_str()))==0 ){
	TH1F* hQCD         = (TH1F*)fin->Get("hQCD");
	TH1F* hAntiIsoKeys = (TH1F*)fin->Get("hAntiIsoKeys") ;
	hAntiIsoKeys->Scale(hQCD->Integral()/hAntiIsoKeys->Integral());
	if( bin_.find("Low")!=string::npos ){
	  //hQCD->Write(Form("QCD%s"    ,suffix.c_str()));
	  hAntiIsoKeys->Scale(hQCD->Integral()/hAntiIsoKeys->Integral());
	  hAntiIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));
	}
	else{
	  //hAntiIsoKeys->Scale(hQCD->Integral()/hAntiIsoKeys->Integral());
	  hAntiIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));
	}

      }
      if(dir->FindObjectAny(Form("W%s"       ,suffix.c_str()))==0 ){

	//((TH1F*)fin->Get("hWKeys"))->Write(Form("W%s"         ,suffix.c_str()));
	TH1F* hW     = (TH1F*)fin->Get("hW");
	TH1F* hWKeys = (TH1F*)fin->Get("hWKeys");
	hWKeys->Scale(hW->Integral()/hWKeys->Integral());
	hWKeys->Write(Form("W%s"           ,suffix.c_str()));

      }
      if(dir->FindObjectAny(Form("ZJ%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hZmj"))->Write(Form("ZJ%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZL%s"       ,suffix.c_str()))==0 ) 
	//((TH1F*)fin->Get("hZmm"))->Write(Form("ZL%s"        ,suffix.c_str()));
	((TH1F*)fin->Get("hZmmKeys"))->Write(Form("ZL%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZLL%s"       ,suffix.c_str()))==0 ) 
	((TH1F*)fin->Get("hZfakes"))->Write(Form("ZLL%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("TT%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("VV%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hVV"))->Write(Form("VV%s"         ,suffix.c_str()));
     
      

    }
    else if(bin_.find("boost")!=string::npos){

      if(dir->FindObjectAny(Form("ZTT%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s"       ,suffix.c_str()));
      if(dir->FindObjectAny(Form("QCD%s"       ,suffix.c_str()))==0 ){
	TH1F* hQCD         = (TH1F*)fin->Get("hQCD");
	TH1F* hAntiIsoKeys = (TH1F*)fin->Get("hAntiIsoKeys");
	//hAntiIsoKeys->Scale(hQCD->Integral()/hAntiIsoKeys->Integral());
	hAntiIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));
	//hQCD->Write(Form("QCD%s"    ,suffix.c_str()));
      }
      if(dir->FindObjectAny(Form("W%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hWKeys"))->Write(Form("W%s"           ,suffix.c_str()));
	//((TH1F*)fin->Get("hW"))->Write(Form("W%s"           ,suffix.c_str()));

      if(dir->FindObjectAny(Form("ZJ%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hZmj"))->Write(Form("ZJ%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZL%s"       ,suffix.c_str()))==0 ) 
	((TH1F*)fin->Get("hZmm"))->Write(Form("ZL%s"        ,suffix.c_str()));

      if(dir->FindObjectAny(Form("ZLL%s"       ,suffix.c_str()))==0 ){
	//TH1F* hZmjKeys    = (TH1F*)fin->Get("hZmjKeys");
	//TH1F* hZmmKeys    = (TH1F*)fin->Get("hZmmKeys");
	TH1F* hZmjKeys    = (TH1F*)fin->Get("hZmj");
	TH1F* hZmmKeys    = (TH1F*)fin->Get("hZmm");
	TH1F* hZfakesKeys = (TH1F*)hZmjKeys->Clone("hZfakesKeys");
	hZfakesKeys->Reset();
	hZfakesKeys->Add(hZmjKeys,1.0);
	hZfakesKeys->Add(hZmmKeys,1.0);  
	hZfakesKeys->Write(Form("ZLL%s"        ,suffix.c_str()));
      }
      if(dir->FindObjectAny(Form("TT%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("VV%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hVV"))->Write(Form("VV%s"         ,suffix.c_str()));
    
    }
    else if(bin_.find("bTag")!=string::npos){

      if(dir->FindObjectAny(Form("ZTT%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s"       ,suffix.c_str()));
      if(dir->FindObjectAny(Form("QCD%s"       ,suffix.c_str()))==0 ){
	TH1F* hQCD         = (TH1F*)fin->Get("hQCD");
	TH1F* hAntiIsoKeys  = (TH1F*)fin->Get("hAntiIsoKeys");

	if(  bin_.find("Low")!=string::npos )
	  hQCD->Write(Form("QCD%s"    ,suffix.c_str()));
	  //hAntiIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));
	else 
	  hAntiIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));
      }
      if(dir->FindObjectAny(Form("W%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hWKeys"))->Write(Form("W%s"           ,suffix.c_str()));
	//((TH1F*)fin->Get("hW"))->Write(Form("W%s"           ,suffix.c_str()));

      if(dir->FindObjectAny(Form("ZJ%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hZmj"))->Write(Form("ZJ%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZL%s"       ,suffix.c_str()))==0 ) 
	((TH1F*)fin->Get("hZmm"))->Write(Form("ZL%s"        ,suffix.c_str()));

      if(dir->FindObjectAny(Form("ZLL%s"       ,suffix.c_str()))==0 ){
	//TH1F* hZmjKeys    = (TH1F*)fin->Get("hZmjKeys");
	//TH1F* hZmmKeys    = (TH1F*)fin->Get("hZmmKeys");
	TH1F* hZmjKeys    = (TH1F*)fin->Get("hZmj");
	TH1F* hZmmKeys    = (TH1F*)fin->Get("hZmm");
	TH1F* hZfakesKeys = (TH1F*)hZmjKeys->Clone("hZfakes");
	hZfakesKeys->Reset();
	hZfakesKeys->Add(hZmjKeys,1.0);
	hZfakesKeys->Add(hZmmKeys,1.0);  
	hZfakesKeys->Write(Form("ZLL%s"        ,suffix.c_str()));
      }
      if(dir->FindObjectAny(Form("TT%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("VV%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hVV"))->Write(Form("VV%s"         ,suffix.c_str()));

   

    }
    else if(bin_.find("vbf")!=string::npos && bin_.find("novbf")==string::npos){

      if(dir->FindObjectAny(Form("ZTT%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s"       ,suffix.c_str()));

      if(dir->FindObjectAny(Form("QCD%s"       ,suffix.c_str()))==0 ){

	TH1F* hAntiIsoKeys  = (TH1F*)fin->Get("hAntiIsoKeys");
	TH1F* hAntiIso      = (TH1F*)fin->Get("hAntiIso");
	TH1F* hLooseIsoKeys = (TH1F*)fin->Get("hLooseIsoKeys");
	hLooseIsoKeys->Scale(hAntiIso->Integral()/hLooseIsoKeys->Integral());
	hLooseIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));
    
	hLooseIsoKeys->Write(Form("QCD%s"       ,suffix.c_str()));

      }

      if(dir->FindObjectAny(Form("W%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hW3Jets"))->Write(Form("W%s"           ,suffix.c_str()));
	//((TH1F*)fin->Get("hWKeys"))->Write(Form("W%s"           ,suffix.c_str()));

      if(dir->FindObjectAny(Form("ZJ%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hZmj"))->Write(Form("ZJ%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZL%s"       ,suffix.c_str()))==0 ) 
	((TH1F*)fin->Get("hZmm"))->Write(Form("ZL%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZLL%s"       ,suffix.c_str()))==0 ){
	TH1F* hZmjKeys    = (TH1F*)fin->Get("hZmjKeys");
	TH1F* hZmmKeys    = (TH1F*)fin->Get("hZmmKeys");
	TH1F* hZmj        = (TH1F*)fin->Get("hZmj");
	TH1F* hZmm        = (TH1F*)fin->Get("hZmm");
	TH1F* hZfakesKeys = (TH1F*)hZmjKeys->Clone("hZfakesKeys");
	hZfakesKeys->Reset();
	if(hZmjKeys->Integral()>0) hZfakesKeys->Add(hZmjKeys,1.0);
	else hZfakesKeys->Add(hZmj,1.0);
	if(hZmmKeys->Integral()>0) hZfakesKeys->Add(hZmmKeys,1.0);  
	else hZfakesKeys->Add(hZmm,1.0);
	hZfakesKeys->Write(Form("ZLL%s"        ,suffix.c_str()));
      }  
      if(dir->FindObjectAny(Form("TT%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("VV%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hVVKeys"))->Write(Form("VV%s"         ,suffix.c_str()));

     

    }
    else if(bin_.find("vh")!=string::npos){

      if(dir->FindObjectAny(Form("ZTT%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s"       ,suffix.c_str()));

      if(dir->FindObjectAny(Form("QCD%s"       ,suffix.c_str()))==0 ){

	TH1F* hAntiIsoKeys  = (TH1F*)fin->Get("hAntiIsoKeys");
	TH1F* hAntiIso      = (TH1F*)fin->Get("hAntiIso");
	TH1F* hLooseIsoKeys = (TH1F*)fin->Get("hLooseIsoKeys");
	hLooseIsoKeys->Scale(hAntiIso->Integral()/hLooseIsoKeys->Integral());
	hLooseIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));
    
	hLooseIsoKeys->Write(Form("QCD%s"       ,suffix.c_str()));

      }

      if(dir->FindObjectAny(Form("W%s"       ,suffix.c_str()))==0 )
	//((TH1F*)fin->Get("hW3JetsKeys"))->Write(Form("W%s"           ,suffix.c_str()));
	((TH1F*)fin->Get("hWKeys"))->Write(Form("W%s"           ,suffix.c_str()));

      if(dir->FindObjectAny(Form("ZJ%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hZmj"))->Write(Form("ZJ%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZL%s"       ,suffix.c_str()))==0 ) 
	((TH1F*)fin->Get("hZmm"))->Write(Form("ZL%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZLL%s"       ,suffix.c_str()))==0 ){
	TH1F* hZmjKeys    = (TH1F*)fin->Get("hZmjKeys");
	TH1F* hZmmKeys    = (TH1F*)fin->Get("hZmmKeys");
	TH1F* hZmj        = (TH1F*)fin->Get("hZmj");
	TH1F* hZmm        = (TH1F*)fin->Get("hZmm");
	TH1F* hZfakesKeys = (TH1F*)hZmjKeys->Clone("hZfakesKeys");
	hZfakesKeys->Reset();
	if(hZmjKeys->Integral()>0) hZfakesKeys->Add(hZmjKeys,1.0);
	else hZfakesKeys->Add(hZmj,1.0);
	if(hZmmKeys->Integral()>0) hZfakesKeys->Add(hZmmKeys,1.0);  
	else hZfakesKeys->Add(hZmm,1.0);
	hZfakesKeys->Write(Form("ZLL%s"        ,suffix.c_str()));
      }  
      if(dir->FindObjectAny(Form("TT%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("VV%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hVVKeys"))->Write(Form("VV%s"         ,suffix.c_str()));



    }
    else{
      if(dir->FindObjectAny(Form("ZTT%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s"       ,suffix.c_str()));
      if(dir->FindObjectAny(Form("QCD%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hAntiIso"))->Write(Form("QCD%s"       ,suffix.c_str()));
      if(dir->FindObjectAny(Form("W%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hW"))->Write(Form("W%s"           ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZJ%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hZmj"))->Write(Form("ZJ%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZL%s"       ,suffix.c_str()))==0 ) 
	((TH1F*)fin->Get("hZmm"))->Write(Form("ZL%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZLL%s"       ,suffix.c_str()))==0 ) 
	((TH1F*)fin->Get("hZfakes"))->Write(Form("ZLL%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("TT%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("VV%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hVV"))->Write(Form("VV%s"         ,suffix.c_str()));

    }


  }


  
  fTemplOut->Close();
  
  // edit the datacards only for the nominal analysis
  if(analysis_.find("Up")!=string::npos || analysis_.find("Down")!=string::npos) 
    return;


  ifstream in;

  char* c = new char[1000];
  in.open(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/templates/eTau_%s_template_v3.txt",     binNameSpace.c_str()));
  //ofstream out(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/datacards/%s/eTau_%s_mH%d_%s.txt", outputDir.Data(), binNameSpace.c_str(), mH_, variable_.c_str()));
  ofstream out(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/datacards/%s/eTau_%s_mH%d.txt", outputDir.Data(), binNameSpace.c_str(), mH_ ));
  out.precision(8);

  while (in.good())
    {
      in.getline(c,1000,'\n');
      if (in.good()){

	string line(c);
	if(line.find("observation")!=string::npos){
	  line.replace( line.find("XXX") , 3 , string(Form("%.0f", ((TH1F*)fin->Get("hData"))->Integral()) )  );
	  out << line << endl;
	}
	else if(line.find("shapes")!=string::npos){
          //line.replace( line.find("XXX") , 3 , string(Form("eTauSM_%s",variable_.c_str()))  );
	  line.replace( line.find("XXX") , 3 , string("eTauSM")  );
	  out << line << endl;
	}
	else if(line.find("process")!=string::npos && line.find("qqH")!=string::npos){

	  //if(!RESCALETO1PB){
	  //line.replace( line.find("XXX") , 3 , string(Form("%d",mH_))  );
	  //line.replace( line.find("YYY") , 3 , string(Form("%d",mH_))  );
	  //line.replace( line.find("KKK") , 3 , string(Form("%d",mH_))  );
	  //}
	  //else{
	  line.replace( line.find("XXX") , 3 , string(Form("   "))  );
	  line.replace( line.find("YYY") , 3 , string(Form("   "))  );
	  line.replace( line.find("KKK") , 3 , string(Form("   "))  );
	  //}
          out << line << endl;
        }
	else if(line.find("rate")!=string::npos){
	  
	  float QCDyield = 0;
	  if(bin_.find("boost")!=string::npos){
	    //QCDyield = ((TH1F*)fin->Get("hQCD"))->Integral(); 
	    QCDyield = ((TH1F*)fin->Get("hAntiIso"))->Integral(); 
	  }
	  else if(bin_.find("novbf")!=string::npos ){
	    if(   bin_.find("Low")!=string::npos )
	      QCDyield = ((TH1F*)fin->Get("hQCD"))->Integral(); 	   
	    else
	      QCDyield = ((TH1F*)fin->Get("hAntiIso"))->Integral();
	  }
	  else if(bin_.find("vbf")!=string::npos && bin_.find("novbf")==string::npos || bin_.find("vh")!=string::npos ){
	    QCDyield = ((TH1F*)fin->Get("hAntiIso"))->Integral(); 	   
	  }
	  else if( bin_.find("bTag")!=string::npos  ){
	    if(   bin_.find("Low")!=string::npos )
	      QCDyield = ((TH1F*)fin->Get("hQCD"))->Integral();
	      //QCDyield = ((TH1F*)fin->Get("hAntiIso"))->Integral(); 	   
	    else
	      QCDyield = ((TH1F*)fin->Get("hAntiIso"))->Integral(); 	   
	  }
	  else{
	    QCDyield = ((TH1F*)fin->Get("hQCD"))->Integral(); 
	  }

	  /////////////////////////////////////////////
	  /////////////////////////////////////////////
	  
	  TH1F* hSgn2 = (TH1F*)fin->Get(Form("hggH%d",mH_));
	  //hSgn2->Scale(1./rescaleggH);
	  
	  TH1F* hSgn1 = (TH1F*)fin->Get(Form("hqqH%d",mH_));
	  //hSgn1->Scale(1./rescaleqqH);
	  
	  TH1F* hSgn3 = (TH1F*)fin->Get(Form("hVH%d",mH_));
	  //hSgn3->Scale(1./rescaleVH);


	  string rate = "rate                                           ";
	  string space = "              ";
	  out << rate ;
	  out <<          hSgn3->Integral()
	      << space << hSgn2->Integral()
	      << space << hSgn1->Integral()
	      << space << ((TH1F*)fin->Get("hDataEmb"))->Integral()
	      << space << QCDyield
	      << space << ((TH1F*)fin->Get("hW"))->Integral();
	  if(binNameSpace.find("0jet_")!=string::npos){
	    out << space << ((TH1F*)fin->Get("hZmj"))->Integral()
		<< space << ((TH1F*)fin->Get("hZmm"))->Integral();
	      }
	  else
	    out << space << ((TH1F*)fin->Get("hZfakes"))->Integral();
	  out << space << ((TH1F*)fin->Get("hTTb"))->Integral()
	      << space << ((TH1F*)fin->Get("hVV"))->Integral()
	      << endl;
	}
	else if(line.find("CMS_htt_ztt_extrap")!=string::npos){

	  float extrapFactor = ((TH1F*)fin->Get("hParameters"))->GetBinContent(8);
	  float extrapError  = ((TH1F*)fin->Get("hParameters"))->GetBinContent(9);

	  float relUnc = 1.+extrapError/extrapFactor;

	  line.replace( line.find("XXX") , 3 , string(Form("%.5f",relUnc))  );
	  line.replace( line.find("YYY") , 3 , string(Form("%.5f",relUnc))  );
	  if( line.find("KKK") != string::npos )
	    line.replace( line.find("KKK") , 3 , string(Form("%.5f",relUnc))  );
	  out << line << endl;

	}
	else if(line.find("CMS_scale_j")!=string::npos){
	  float VBFrel = ((TH1F*)fin->Get(Form("hVH%d",mH_)))->Integral()>0 ? 
	    TMath::Max(TMath::Abs((((TH1F*)fin_jUp->Get(Form("hVH%d",mH_)))->Integral()/((TH1F*)fin->Get(Form("hVH%d",mH_)))->Integral())),
		       TMath::Abs((((TH1F*)fin_jDown->Get(Form("hVH%d",mH_)))->Integral()/((TH1F*)fin->Get(Form("hVH%d",mH_)))->Integral()))
		       ): 2.0;

	  float VBFrelUp   = ((TH1F*)fin->Get(Form("hVH%d",mH_)))->Integral()>0 ?
	    (((TH1F*)fin_jUp->Get(Form("hVH%d",mH_)))->Integral()/((TH1F*)fin->Get(Form("hVH%d",mH_)))->Integral())  : 2.0;

	  float VBFrelDown = ((TH1F*)fin->Get(Form("hVH%d",mH_)))->Integral()>0 ?
	    (((TH1F*)fin_jDown->Get(Form("hVH%d",mH_)))->Integral()/((TH1F*)fin->Get(Form("hVH%d",mH_)))->Integral()) : 2.0;

	  float SMrel  = ((TH1F*)fin->Get(Form("hggH%d",mH_)))->Integral()>0 ? 
	    TMath::Max(TMath::Abs((((TH1F*)fin_jUp->Get(Form("hggH%d",mH_)))->Integral()/((TH1F*)fin->Get(Form("hggH%d",mH_)))->Integral())),
		       TMath::Abs((((TH1F*)fin_jDown->Get(Form("hggH%d",mH_)))->Integral()/((TH1F*)fin->Get(Form("hggH%d",mH_)))->Integral()))
		       ) : 2.0;

	  float SMrelUp   = ((TH1F*)fin->Get(Form("hggH%d",mH_)))->Integral()>0 ?
	    (((TH1F*)fin_jUp->Get(Form("hggH%d",mH_)))->Integral()/((TH1F*)fin->Get(Form("hggH%d",mH_)))->Integral()) : 2.0;

	  float SMrelDown = ((TH1F*)fin->Get(Form("hggH%d",mH_)))->Integral()>0 ?
	    (((TH1F*)fin_jDown->Get(Form("hggH%d",mH_)))->Integral()/((TH1F*)fin->Get(Form("hggH%d",mH_)))->Integral()) : 2.0;

	  float VHrel  = TMath::Max(TMath::Abs((((TH1F*)fin_jUp->Get(Form("hqqH%d",mH_)))->Integral()/((TH1F*)fin->Get(Form("hqqH%d",mH_)))->Integral())),
				    TMath::Abs((((TH1F*)fin_jDown->Get(Form("hqqH%d",mH_)))->Integral()/((TH1F*)fin->Get(Form("hqqH%d",mH_)))->Integral()))
				    );

	  float VHrelUp   = ((TH1F*)fin->Get(Form("hqqH%d",mH_)))->Integral()>0 ?
	    (((TH1F*)fin_jUp->Get(Form("hqqH%d",mH_)))->Integral()/((TH1F*)fin->Get(Form("hqqH%d",mH_)))->Integral()) : 2.0;
	  float VHrelDown = ((TH1F*)fin->Get(Form("hqqH%d",mH_)))->Integral()>0 ?
	    (((TH1F*)fin_jDown->Get(Form("hqqH%d",mH_)))->Integral()/((TH1F*)fin->Get(Form("hqqH%d",mH_)))->Integral()): 2.0;


	  float ZTTrel = ((TH1F*)fin->Get("hZtt"))->Integral()>0 ? 
	    TMath::Max(TMath::Abs((((TH1F*)fin_jUp->Get("hZtt"))->Integral()/((TH1F*)fin->Get("hZtt"))->Integral())),
		       TMath::Abs((((TH1F*)fin_jDown->Get("hZtt"))->Integral()/((TH1F*)fin->Get("hZtt"))->Integral()))
		       ) : 2.0;

	  float ZTTrelUp   = ((TH1F*)fin->Get("hZtt"))->Integral()>0 ? 
	    (((TH1F*)fin_jUp->Get("hZtt"))->Integral()/((TH1F*)fin->Get("hZtt"))->Integral()) : 2.0;
	  float ZTTrelDown = ((TH1F*)fin->Get("hZtt"))->Integral()>0 ?
	    (((TH1F*)fin_jDown->Get("hZtt"))->Integral()/((TH1F*)fin->Get("hZtt"))->Integral()) : 2.0;


	  float TTrel  = ((TH1F*)fin->Get("hTTb"))->Integral()>0 ?
	    TMath::Max(TMath::Abs((((TH1F*)fin_jUp->Get("hTTb"))->Integral()/((TH1F*)fin->Get("hTTb"))->Integral())),
		       TMath::Abs((((TH1F*)fin_jDown->Get("hTTb"))->Integral()/((TH1F*)fin->Get("hTTb"))->Integral()))
		       ) : 2.0;

	  float TTrelUp   = ((TH1F*)fin->Get("hTTb"))->Integral()>0 ?
	    (((TH1F*)fin_jUp->Get("hTTb"))->Integral()/((TH1F*)fin->Get("hTTb"))->Integral()) : 2.0;
	  float TTrelDown = ((TH1F*)fin->Get("hTTb"))->Integral()>0 ?
	    (((TH1F*)fin_jDown->Get("hTTb"))->Integral()/((TH1F*)fin->Get("hTTb"))->Integral()) : 2.0;


	  float VVrel  = ((TH1F*)fin->Get("hVV"))->Integral()>0 ?
	    TMath::Max(TMath::Abs((((TH1F*)fin_jUp->Get("hVV"))->Integral()/((TH1F*)fin->Get("hVV"))->Integral())),
		       TMath::Abs((((TH1F*)fin_jDown->Get("hVV"))->Integral()/((TH1F*)fin->Get("hVV"))->Integral()))
		       ) : 2.0;

	  float VVrelUp   = ((TH1F*)fin->Get("hVV"))->Integral()>0 ?
	    (((TH1F*)fin_jUp->Get("hVV"))->Integral()/((TH1F*)fin->Get("hVV"))->Integral()) : 2.0;
	  float VVrelDown = ((TH1F*)fin->Get("hVV"))->Integral()>0 ?
	    (((TH1F*)fin_jDown->Get("hVV"))->Integral()/((TH1F*)fin->Get("hVV"))->Integral()) : 2.0;


	  if(TMath::Abs(VBFrelUp-1)<0.001)   VBFrelUp  *=(1+0.001);
	  if(TMath::Abs(VBFrelDown-1)<0.001) VBFrelDown*=(1+0.001);
	  if(TMath::Abs(SMrelUp-1)<0.001)    SMrelUp   *=(1+0.001);
	  if(TMath::Abs(SMrelDown-1)<0.001)  SMrelDown *=(1+0.001);
	  if(TMath::Abs(VHrelUp-1)<0.001)    VHrelUp   *=(1+0.001);
	  if(TMath::Abs(VHrelDown-1)<0.001)  VHrelDown *=(1+0.001);
	  if(TMath::Abs(ZTTrelUp-1)<0.001)   ZTTrelUp  *=(1+0.001);
	  if(TMath::Abs(ZTTrelDown-1)<0.001) ZTTrelDown*=(1+0.001);
	  if(TMath::Abs(TTrelUp-1)<0.001)    TTrelUp   *=(1+0.001);
	  if(TMath::Abs(TTrelDown-1)<0.001)  TTrelDown *=(1+0.001);
	  if(TMath::Abs(VVrelUp-1)<0.001)    VVrelUp   *=(1+0.001);
	  if(TMath::Abs(VVrelDown-1)<0.001)  VVrelDown *=(1+0.001);

	  string space      = "                   ";
	  string longspace  = "                      ";
	  string shortspace = "     ";
	  out << "CMS_scale_j";
	  out << longspace << "lnN" << shortspace;
	  if(binNameSpace.find("0jet_")!=string::npos) 
	    out <<          string(Form("%.2f",1+TMath::Max(TMath::Abs(-1+VHrelDown), TMath::Abs(-1+VHrelUp) ) * TMath::Abs(-1+VHrelUp)/(-1+VHrelUp)   )); // VH
	  else
	    out << space << string(Form("%.2f",1+TMath::Max(TMath::Abs(-1+SMrelDown), TMath::Abs(-1+SMrelUp))  * TMath::Abs(-1+SMrelUp)/(-1+SMrelUp)   )); // VH = GGF
	  out << space << string(Form("%.2f",1+TMath::Max(  TMath::Abs(-1+SMrelDown), TMath::Abs(-1+SMrelUp))  * TMath::Abs(-1+SMrelUp)/(-1+SMrelUp)   ));   // GGF
	  out << space << string(Form("%.2f",1+TMath::Max(  TMath::Abs(-1+VBFrelDown),TMath::Abs(-1+VBFrelUp))* TMath::Abs(-1+VBFrelUp)/(-1+VBFrelUp) ));   // VBF
	  out << space << "-" ; //string(Form("%.2f",1+(-1+ZTTrelDown)*1.0));  // ZTT
	  out << space << "-" ;                                                // QCD
	  out << space << "-" ; //string(Form("%.2f",1+(-1+ZTTrelDown)*1.0));  // W
	  out << space << "-" ; //string(Form("%.2f",1+(-1+ZTTrelDown)*1.0));  // Z+fakes

	  if(binNameSpace.find("0jet_")!=string::npos) 
	    out << space << "-" ; //string(Form("%.2f",1+(-1+ZTTrelDown)*1.0));

	  out << space << string(Form("%.2f",1+TMath::Max(  TMath::Abs(-1+TTrelDown),TMath::Abs(-1+TTrelUp))  * TMath::Abs(-1+TTrelUp)/(-1+TTrelUp) ));  // TT

	  //if(binNameSpace.find("0jet_")!=string::npos) 
	  out << space << string(Form("%.2f",1+TMath::Max(TMath::Abs(-1+VVrelDown),TMath::Abs(-1+VVrelUp) ) * TMath::Abs(-1+VVrelUp)/(-1+VVrelUp)));   // VV
	  //else
	  //out << space << string(Form("%.2f",1+TMath::Max(TMath::Abs(-1+TTrelDown),TMath::Abs(-1+TTrelUp))  * TMath::Abs(-1+TTrelUp)/(-1+TTrelUp) ));  // TT

	  out << space << "JEC uncertainty";
	  out << endl;
	}
	else if(line.find( Form("CMS_htt_eTau_%s_QCDNorm",binNameSpace.c_str()) )!=string::npos){
	  line.replace( line.find("XXX") , 3 , string(Form("%.0f",((TH1F*)fin->Get("hParameters"))->GetBinContent(5)))  );
	  line.replace( line.find("YYY") , 3 , string(Form("%.4f",((TH1F*)fin->Get("hParameters"))->GetBinContent(6)))  );
	  out << line << endl;
	}
	else if(line.find( Form("CMS_htt_eTau_%s_WNorm",binNameSpace.c_str()) )!=string::npos){
	  line.replace( line.find("XXX") , 3 , string(Form("%.0f",((TH1F*)fin->Get("hParameters"))->GetBinContent(2)))  );
	  line.replace( line.find("YYY") , 3 , string(Form("%.4f",((TH1F*)fin->Get("hParameters"))->GetBinContent(1)))  );
	  out << line << endl;
	}
	else{
	  out << c << endl;
	}
      
      }
 
   }


  return;

}



void produceAll(  TString outputDir = "May2012/Reload_PreApproval" ){

  vector<string> variables;
  vector<int> mH;

  //variables.push_back("diTauVisMass");
  variables.push_back("diTauNSVfitMass");
  //variables.push_back("diTauSVFitMassCal0");
  //variables.push_back("diTauSVFitMassCal1");
  //variables.push_back("diTauSVFitMassCal2");


  //mH.push_back(105);
  //mH.push_back(110);
  //mH.push_back(115);
  //mH.push_back(120);
  //mH.push_back(125);
  //mH.push_back(130);
  //mH.push_back(135);
  mH.push_back(140);
  mH.push_back(145);

  for(unsigned int i = 0 ; i < variables.size(); i++){
    for(unsigned j = 0; j < mH.size(); j++){

      produce(mH[j],variables[i], ""        , "novbfLow", outputDir);
      produce(mH[j],variables[i], "TauUp"   , "novbfLow", outputDir);
      produce(mH[j],variables[i], "TauDown" , "novbfLow", outputDir);
      produce(mH[j],variables[i], "JetUp"   , "novbfLow", outputDir);
      produce(mH[j],variables[i], "JetDown" , "novbfLow", outputDir);
      
      produce(mH[j],variables[i], ""        , "novbfHigh", outputDir);
      produce(mH[j],variables[i], "TauUp"   , "novbfHigh", outputDir);
      produce(mH[j],variables[i], "TauDown" , "novbfHigh", outputDir);
      produce(mH[j],variables[i], "JetUp"   , "novbfHigh", outputDir);
      produce(mH[j],variables[i], "JetDown" , "novbfHigh", outputDir);

      produce(mH[j],variables[i], ""        , "boostLow", outputDir);
      produce(mH[j],variables[i], "TauUp"   , "boostLow", outputDir);
      produce(mH[j],variables[i], "TauDown" , "boostLow", outputDir);
      produce(mH[j],variables[i], "JetUp"   , "boostLow", outputDir);
      produce(mH[j],variables[i], "JetDown" , "boostLow", outputDir);

      produce(mH[j],variables[i], ""        , "boostHigh", outputDir);
      produce(mH[j],variables[i], "TauUp"   , "boostHigh", outputDir);
      produce(mH[j],variables[i], "TauDown" , "boostHigh", outputDir);
      produce(mH[j],variables[i], "JetUp"   , "boostHigh", outputDir);
      produce(mH[j],variables[i], "JetDown" , "boostHigh", outputDir);

      produce(mH[j],variables[i], ""        , "bTagLow", outputDir);
      produce(mH[j],variables[i], "TauUp"   , "bTagLow", outputDir);
      produce(mH[j],variables[i], "TauDown" , "bTagLow", outputDir);
      produce(mH[j],variables[i], "JetUp"   , "bTagLow", outputDir);
      produce(mH[j],variables[i], "JetDown" , "bTagLow", outputDir);

      produce(mH[j],variables[i], ""        , "bTagHigh", outputDir);
      produce(mH[j],variables[i], "TauUp"   , "bTagHigh", outputDir);
      produce(mH[j],variables[i], "TauDown" , "bTagHigh", outputDir);
      produce(mH[j],variables[i], "JetUp"   , "bTagHigh", outputDir);
      produce(mH[j],variables[i], "JetDown" , "bTagHigh", outputDir);

      produce(mH[j],variables[i], ""        , "vbf", outputDir);
      produce(mH[j],variables[i], "TauUp"   , "vbf", outputDir);
      produce(mH[j],variables[i], "TauDown" , "vbf", outputDir);
      produce(mH[j],variables[i], "JetUp"   , "vbf", outputDir);
      produce(mH[j],variables[i], "JetDown" , "vbf", outputDir);

      produce(mH[j],variables[i], ""        , "vh", outputDir);
      produce(mH[j],variables[i], "TauUp"   , "vh", outputDir);
      produce(mH[j],variables[i], "TauDown" , "vh", outputDir);
      produce(mH[j],variables[i], "JetUp"   , "vh", outputDir);
      produce(mH[j],variables[i], "JetDown" , "vh", outputDir);

     
  
      //produce(mH[j],variables[i], ""        , "novbf", outputDir);
      //produce(mH[j],variables[i], "TauUp"   , "novbf", outputDir);
      //produce(mH[j],variables[i], "TauDown" , "novbf", outputDir);
      //produce(mH[j],variables[i], "JetUp"   , "novbf", outputDir);
      //produce(mH[j],variables[i], "JetDown" , "novbf", outputDir);
  
    }

  }

}




int main(int argc, const char* argv[])
{

  std::cout << "produce()" << std::endl;
  gROOT->SetBatch(true);
 

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  produceAll();

}

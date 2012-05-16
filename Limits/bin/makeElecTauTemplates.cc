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

using namespace std;

void produce(   
	     int mH_=120,
	     string variable_  = "diTauVisMass",
	     string analysis_  = "",
	     string bin_       = "inclusive",
	     TString outputDir = "Nov2011"
	     ){


  cout << "Now doing mass mH=" << mH_ << ", for variable " << variable_ << " analysis " << analysis_ << " and bin " << bin_ << endl;
  TFile* fin = new TFile(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/histograms/%s/eTau_mH%d_%s_%s_%s.root", outputDir.Data(), mH_, bin_.c_str() , analysis_.c_str(), variable_.c_str()), "READ");

  TFile* fin_jUp   = new TFile(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/histograms/%s/eTau_mH%d_%s_JetUp_%s.root",   outputDir.Data(), mH_, bin_.c_str() , variable_.c_str()), "READ");
  TFile* fin_jDown = new TFile(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/histograms/%s/eTau_mH%d_%s_JetDown_%s.root", outputDir.Data(), mH_, bin_.c_str() , variable_.c_str()), "READ");
  TFile* fin_tUp   = new TFile(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/histograms/%s/eTau_mH%d_%s_TauUp_%s.root",   outputDir.Data(), mH_, bin_.c_str() , variable_.c_str()), "READ");
  TFile* fin_tDown = new TFile(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/histograms/%s/eTau_mH%d_%s_TauDown_%s.root", outputDir.Data(), mH_, bin_.c_str() , variable_.c_str()), "READ");
  TFile* fin_nominal = new TFile(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/histograms/%s/eTau_mH%d_%s__%s.root", outputDir.Data(), mH_, bin_.c_str() , variable_.c_str()), "READ");

  string binNameSpace = "";
  if(bin_.find("inclusive")!=string::npos)      
    binNameSpace =  "inclusive";
  else if(bin_.find("novbf")!=string::npos) 
    binNameSpace =  "SM0";
  else if(bin_.find("boost")!=string::npos) 
    binNameSpace =  "SM1";
  else if(bin_.find("vbf")!=string::npos && bin_.find("novbf")==string::npos) 
    binNameSpace =  "SM2";
  else if(bin_.find("twoJets")!=string::npos) 
    binNameSpace =  "SMpre2";
  else if(bin_.find("oneJet")!=string::npos) 
    binNameSpace =  "SMpre2a";

  TFile* fTemplOut = new TFile(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/datacards/%s/eTauSM_%s.root",outputDir.Data(), variable_.c_str()),"UPDATE");
  
  string suffix = "";
  if(analysis_.find("TauUp")!=string::npos)
    suffix = "_CMS_scale_tUp";
  else if(analysis_.find("TauDown")!=string::npos)
    suffix = "_CMS_scale_tDown";
  else if(analysis_.find("ElecUp")!=string::npos)
    suffix = "_CMS_scale_eUp";
  else if(analysis_.find("ElecDown")!=string::npos)
    suffix = "_CMS_scale_eDown";
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

    if(bin_.find("novbf")!=string::npos){
      ((TH1F*)fin->Get("hData"))->Write("data_obs");
      ((TH1F*)fin->Get("hSgn1"))->Write(Form("VBF%d%s" ,mH_,suffix.c_str()));
      ((TH1F*)fin->Get("hSgn2"))->Write(Form("SM%d%s"  ,mH_,suffix.c_str()));
      ((TH1F*)fin->Get("hSgn3"))->Write(Form("VH%d%s"  ,mH_,suffix.c_str()));
      ((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s",suffix.c_str()));
      TH1F* hQCD         = (TH1F*)fin->Get("hQCD");
      TH1F* hAntiIsoKeys = (TH1F*)fin->Get("hAntiIsoKeys");
      hAntiIsoKeys->Scale(hQCD->Integral()/hAntiIsoKeys->Integral());
      hAntiIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));
      //hQCD->Write(Form("QCD%s"    ,suffix.c_str()));
      ((TH1F*)fin->Get("hW"))->Write(Form("W%s"           ,suffix.c_str()));
      ((TH1F*)fin->Get("hZmj"))->Write(Form("ZJ%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hZmm"))->Write(Form("ZL%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hZfakes"))->Write(Form("ZLL%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hVV"))->Write(Form("VV%s"         ,suffix.c_str()));

    }
    else if(bin_.find("boost")!=string::npos){
      ((TH1F*)fin->Get("hData"))->Write("data_obs");
      ((TH1F*)fin->Get("hSgn1"))->Write(Form("VBF%d%s" ,mH_,suffix.c_str()));
      ((TH1F*)fin->Get("hSgn2"))->Write(Form("SM%d%s"  ,mH_,suffix.c_str()));
      ((TH1F*)fin->Get("hSgn3"))->Write(Form("VH%d%s"  ,mH_,suffix.c_str()));
      ((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s",suffix.c_str()));
      ((TH1F*)fin->Get("hAntiIsoKeys"))->Write(Form("QCD%s"    ,suffix.c_str()));
      ((TH1F*)fin->Get("hWKeys"))->Write(Form("W%s"           ,suffix.c_str()));
      TH1F* hZmjKeys    = (TH1F*)fin->Get("hZmjKeys");
      TH1F* hZmmKeys    = (TH1F*)fin->Get("hZmmKeys");
      TH1F* hZfakesKeys = (TH1F*)fin->Get("hZfakes");
      hZfakesKeys->Reset();
      hZfakesKeys->Add(hZmjKeys,1.0);
      hZfakesKeys->Add(hZmmKeys,1.0);
      hZfakesKeys->Write(Form("ZLL%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hVV"))->Write(Form("VV%s"         ,suffix.c_str()));
    }
    else if(bin_.find("vbf")!=string::npos && bin_.find("novbf")==string::npos){
      ((TH1F*)fin->Get("hData"))->Write("data_obs");
      ((TH1F*)fin->Get("hSgn1"))->Write(Form("VBF%d%s" ,mH_,suffix.c_str()));
      ((TH1F*)fin->Get("hSgn2"))->Write(Form("SM%d%s"  ,mH_,suffix.c_str()));
      ((TH1F*)fin->Get("hSgn3"))->Write(Form("VH%d%s"  ,mH_,suffix.c_str()));
      ((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s",suffix.c_str()));
      ((TH1F*)fin->Get("hAntiIsoKeys"))->Write(Form("QCD%s"    ,suffix.c_str()));
      //((TH1F*)fin->Get("hW3JetsKeys"))->Write(Form("W%s"           ,suffix.c_str()));
      ((TH1F*)fin->Get("hW"))->Write(Form("W%s"           ,suffix.c_str()));
      TH1F* hZmjKeys    = (TH1F*)fin->Get("hZmjKeys");
      TH1F* hZmmKeys    = (TH1F*)fin->Get("hZmmKeys");
      TH1F* hZfakesKeys = (TH1F*)hZmjKeys->Clone("hZfakesKeys");
      hZfakesKeys->Reset();
      hZfakesKeys->Add(hZmjKeys,1.0);
      hZfakesKeys->Add(hZmmKeys,1.0);  
      hZfakesKeys->Write(Form("ZLL%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
      ((TH1F*)fin->Get("hVVKeys"))->Write(Form("VV%s"         ,suffix.c_str()));
    }
    else{
      ((TH1F*)fin->Get("hData"))->Write("data_obs");
      ((TH1F*)fin->Get("hSgn1"))->Write(Form("VBF%d%s" ,mH_,suffix.c_str()));
      ((TH1F*)fin->Get("hSgn2"))->Write(Form("SM%d%s"  ,mH_,suffix.c_str()));
      ((TH1F*)fin->Get("hSgn3"))->Write(Form("VH%d%s"  ,mH_,suffix.c_str()));
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

    if(bin_.find("novbf")!=string::npos){

      if(dir->FindObjectAny(Form("ZTT%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s"       ,suffix.c_str()));
      if(dir->FindObjectAny(Form("QCD%s"       ,suffix.c_str()))==0 ){
	TH1F* hQCD         = (TH1F*)fin->Get("hQCD");
	TH1F* hAntiIsoKeys = (TH1F*)fin->Get("hAntiIsoKeys");
	hAntiIsoKeys->Scale(hQCD->Integral()/hAntiIsoKeys->Integral());
	hAntiIsoKeys->Write(Form("QCD%s"    ,suffix.c_str()));
	//hQCD->Write(Form("QCD%s"    ,suffix.c_str()));     
      }
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
      if(dir->FindObjectAny(Form("VBF%d%s"         ,mH_,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hSgn1"))->Write(Form("VBF%d%s",mH_  ,suffix.c_str()));
      if(dir->FindObjectAny(Form("SM%d%s"          , mH_,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hSgn2"))->Write(Form("SM%d%s" ,mH_,suffix.c_str()));
      if(dir->FindObjectAny(Form("VH%d%s"          , mH_,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hSgn3"))->Write(Form("VH%d%s" ,mH_,suffix.c_str()));
    }
    else if(bin_.find("boost")!=string::npos){
      if(dir->FindObjectAny(Form("ZTT%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s"       ,suffix.c_str()));
      if(dir->FindObjectAny(Form("QCD%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hAntiIsoKeys"))->Write(Form("QCD%s"       ,suffix.c_str()));
      if(dir->FindObjectAny(Form("W%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hWKeys"))->Write(Form("W%s"           ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZJ%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hZmj"))->Write(Form("ZJ%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZL%s"       ,suffix.c_str()))==0 ) 
	((TH1F*)fin->Get("hZmm"))->Write(Form("ZL%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZLL%s"       ,suffix.c_str()))==0 ){
	TH1F* hZmjKeys    = (TH1F*)fin->Get("hZmjKeys");
	TH1F* hZmmKeys    = (TH1F*)fin->Get("hZmmKeys");
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
      if(dir->FindObjectAny(Form("VBF%d%s"         ,mH_,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hSgn1"))->Write(Form("VBF%d%s",mH_  ,suffix.c_str()));
      if(dir->FindObjectAny(Form("SM%d%s"          , mH_,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hSgn2"))->Write(Form("SM%d%s" ,mH_,suffix.c_str()));
      if(dir->FindObjectAny(Form("VH%d%s"          , mH_,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hSgn3"))->Write(Form("VH%d%s" ,mH_,suffix.c_str()));
    }
    else if(bin_.find("vbf")!=string::npos && bin_.find("novbf")==string::npos){
       if(dir->FindObjectAny(Form("ZTT%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hDataEmb"))->Write(Form("ZTT%s"       ,suffix.c_str()));
      if(dir->FindObjectAny(Form("QCD%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hAntiIsoKeys"))->Write(Form("QCD%s"       ,suffix.c_str()));
      if(dir->FindObjectAny(Form("W%s"       ,suffix.c_str()))==0 )
	//((TH1F*)fin->Get("hW3JetsKeys"))->Write(Form("W%s"           ,suffix.c_str()));
	((TH1F*)fin->Get("hW"))->Write(Form("W%s"           ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZJ%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hZmj"))->Write(Form("ZJ%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZL%s"       ,suffix.c_str()))==0 ) 
	((TH1F*)fin->Get("hZmm"))->Write(Form("ZL%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("ZLL%s"       ,suffix.c_str()))==0 ){
	TH1F* hZmjKeys    = (TH1F*)fin->Get("hZmjKeys");
	TH1F* hZmmKeys    = (TH1F*)fin->Get("hZmmKeys");
	TH1F* hZfakesKeys = (TH1F*)hZmjKeys->Clone("hZfakesKeys");
	hZfakesKeys->Reset();
	hZfakesKeys->Add(hZmjKeys,1.0);
	hZfakesKeys->Add(hZmmKeys,1.0);  
	hZfakesKeys->Write(Form("ZLL%s"        ,suffix.c_str()));
      }  
      if(dir->FindObjectAny(Form("TT%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
      if(dir->FindObjectAny(Form("VV%s"       ,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hVVKeys"))->Write(Form("VV%s"         ,suffix.c_str()));
      if(dir->FindObjectAny(Form("VBF%d%s"         ,mH_,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hSgn1"))->Write(Form("VBF%d%s",mH_  ,suffix.c_str()));
      if(dir->FindObjectAny(Form("SM%d%s"          , mH_,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hSgn2"))->Write(Form("SM%d%s" ,mH_,suffix.c_str()));
      if(dir->FindObjectAny(Form("VH%d%s"          , mH_,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hSgn3"))->Write(Form("VH%d%s" ,mH_,suffix.c_str()));
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
      if(dir->FindObjectAny(Form("VBF%d%s"         ,mH_,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hSgn1"))->Write(Form("VBF%d%s",mH_  ,suffix.c_str()));
      if(dir->FindObjectAny(Form("SM%d%s"          , mH_,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hSgn2"))->Write(Form("SM%d%s" ,mH_,suffix.c_str()));
      if(dir->FindObjectAny(Form("VH%d%s"          , mH_,suffix.c_str()))==0 )
	((TH1F*)fin->Get("hSgn3"))->Write(Form("VH%d%s" ,mH_,suffix.c_str()));
    }


  }


  fTemplOut->Close();
  
  // edit the datacards only for the nominal analysis
  if(analysis_.find("Up")!=string::npos || analysis_.find("Down")!=string::npos) 
    return;

  ifstream in;

  char* c = new char[1000];
  in.open(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/templates/eTau_%s_template_v2.txt",     binNameSpace.c_str()));
  ofstream out(Form("/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/htautau/datacards/%s/eTau_%s_mH%d_%s.txt", outputDir.Data(), binNameSpace.c_str(), mH_, variable_.c_str()));
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
          line.replace( line.find("XXX") , 3 , string(Form("eTauSM_%s",variable_.c_str()))  );
	  out << line << endl;
	}
	else if(line.find("process")!=string::npos && line.find("VBF")!=string::npos){
          line.replace( line.find("XXX") , 3 , string(Form("%d",mH_))  );
          line.replace( line.find("YYY") , 3 , string(Form("%d",mH_))  );
	  line.replace( line.find("KKK") , 3 , string(Form("%d",mH_))  );
          out << line << endl;
        }
	else if(line.find("rate")!=string::npos){

	  float QCDyield = 0;
	  if(bin_.find("novbf")!=string::npos){
	    QCDyield = ((TH1F*)fin->Get("hQCD"))->Integral(); 
	  }
	  else if(bin_.find("vbf")!=string::npos && bin_.find("novbf")==string::npos){
	    QCDyield = ((TH1F*)fin->Get("hAntiIsoKeys"))->Integral(); 
	  }
	  else if(bin_.find("boost")!=string::npos){
	    QCDyield = ((TH1F*)fin->Get("hAntiIsoKeys"))->Integral(); 
	  }
	  else{
	    QCDyield = ((TH1F*)fin->Get("hQCD"))->Integral(); 
	  }
	  string rate = "rate                                           ";
	  string space = "              ";
	  out << rate ;
	  out <<          ((TH1F*)fin->Get("hSgn3"))->Integral()
	      << space << ((TH1F*)fin->Get("hSgn2"))->Integral()
	      << space << ((TH1F*)fin->Get("hSgn1"))->Integral()
	      << space << ((TH1F*)fin->Get("hDataEmb"))->Integral()
	      << space << QCDyield
	      << space << ((TH1F*)fin->Get("hW"))->Integral();
	  if(binNameSpace.find("SM0")!=string::npos){
	    out << space << ((TH1F*)fin->Get("hZmj"))->Integral()
		<< space << ((TH1F*)fin->Get("hZmm"))->Integral();
	      }
	  else
	    out << space << ((TH1F*)fin->Get("hZmj"))->Integral()+((TH1F*)fin->Get("hZmm"))->Integral();
	  out << space << ((TH1F*)fin->Get("hTTb"))->Integral()
	      << space << ((TH1F*)fin->Get("hVV"))->Integral()
	      << endl;
	}
	else if(line.find("CMS_scale_j")!=string::npos){
	  float VBFrel = TMath::Max(TMath::Abs((((TH1F*)fin_jUp->Get("hSgn1"))->Integral()/((TH1F*)fin->Get("hSgn1"))->Integral())),
				    TMath::Abs((((TH1F*)fin_jDown->Get("hSgn1"))->Integral()/((TH1F*)fin->Get("hSgn1"))->Integral()))
				    );
	  float VBFrelUp   = (((TH1F*)fin_jUp->Get("hSgn1"))->Integral()/((TH1F*)fin->Get("hSgn1"))->Integral());
	  float VBFrelDown = (((TH1F*)fin_jDown->Get("hSgn1"))->Integral()/((TH1F*)fin->Get("hSgn1"))->Integral());

	  float SMrel  = TMath::Max(TMath::Abs((((TH1F*)fin_jUp->Get("hSgn2"))->Integral()/((TH1F*)fin->Get("hSgn2"))->Integral())),
				    TMath::Abs((((TH1F*)fin_jDown->Get("hSgn2"))->Integral()/((TH1F*)fin->Get("hSgn2"))->Integral()))
				    );

	  float SMrelUp   = (((TH1F*)fin_jUp->Get("hSgn2"))->Integral()/((TH1F*)fin->Get("hSgn2"))->Integral());
	  float SMrelDown = (((TH1F*)fin_jDown->Get("hSgn2"))->Integral()/((TH1F*)fin->Get("hSgn2"))->Integral());

	  float VHrel  = TMath::Max(TMath::Abs((((TH1F*)fin_jUp->Get("hSgn3"))->Integral()/((TH1F*)fin->Get("hSgn3"))->Integral())),
				    TMath::Abs((((TH1F*)fin_jDown->Get("hSgn3"))->Integral()/((TH1F*)fin->Get("hSgn3"))->Integral()))
				    );

	  float VHrelUp   = (((TH1F*)fin_jUp->Get("hSgn3"))->Integral()/((TH1F*)fin->Get("hSgn3"))->Integral());
	  float VHrelDown = (((TH1F*)fin_jDown->Get("hSgn3"))->Integral()/((TH1F*)fin->Get("hSgn3"))->Integral());


	  float ZTTrel = TMath::Max(TMath::Abs((((TH1F*)fin_jUp->Get("hZtt"))->Integral()/((TH1F*)fin->Get("hZtt"))->Integral())),
				    TMath::Abs((((TH1F*)fin_jDown->Get("hZtt"))->Integral()/((TH1F*)fin->Get("hZtt"))->Integral()))
				    );

	  float ZTTrelUp   = (((TH1F*)fin_jUp->Get("hZtt"))->Integral()/((TH1F*)fin->Get("hZtt"))->Integral());
	  float ZTTrelDown = (((TH1F*)fin_jDown->Get("hZtt"))->Integral()/((TH1F*)fin->Get("hZtt"))->Integral());


	  float TTrel  = TMath::Max(TMath::Abs((((TH1F*)fin_jUp->Get("hTTb"))->Integral()/((TH1F*)fin->Get("hTTb"))->Integral())),
				    TMath::Abs((((TH1F*)fin_jDown->Get("hTTb"))->Integral()/((TH1F*)fin->Get("hTTb"))->Integral()))
				    );

	  float TTrelUp   = (((TH1F*)fin_jUp->Get("hTTb"))->Integral()/((TH1F*)fin->Get("hTTb"))->Integral());
	  float TTrelDown = (((TH1F*)fin_jDown->Get("hTTb"))->Integral()/((TH1F*)fin->Get("hTTb"))->Integral());


	  float VVrel  = TMath::Max(TMath::Abs((((TH1F*)fin_jUp->Get("hVV"))->Integral()/((TH1F*)fin->Get("hVV"))->Integral())),
				    TMath::Abs((((TH1F*)fin_jDown->Get("hVV"))->Integral()/((TH1F*)fin->Get("hVV"))->Integral()))
				    );

	  float VVrelUp   = (((TH1F*)fin_jUp->Get("hVV"))->Integral()/((TH1F*)fin->Get("hVV"))->Integral());
	  float VVrelDown = (((TH1F*)fin_jDown->Get("hVV"))->Integral()/((TH1F*)fin->Get("hVV"))->Integral());


	  string space      = "                   ";
	  string longspace  = "                      ";
	  string shortspace = "     ";
	  out << "CMS_scale_j";
	  out << longspace << "lnN" << shortspace;
	  if(binNameSpace.find("SM0")!=string::npos) 
	    out <<          string(Form("%.2f",1+TMath::Max(TMath::Abs(-1+VHrelDown), TMath::Abs(-1+VHrelUp) ) * TMath::Abs(-1+VHrelUp)/(-1+VHrelUp)   )); // VH
	  else
	    out << space << string(Form("%.2f",1+TMath::Max(TMath::Abs(-1+SMrelDown), TMath::Abs(-1+SMrelUp))  * TMath::Abs(-1+SMrelUp)/(-1+SMrelUp)   )); // VH = GGF
	  out << space << string(Form("%.2f",1+TMath::Max(  TMath::Abs(-1+SMrelDown), TMath::Abs(-1+SMrelUp))  * TMath::Abs(-1+SMrelUp)/(-1+SMrelUp)   ));   // GGF
	  out << space << string(Form("%.2f",1+TMath::Max(  TMath::Abs(-1+VBFrelDown),TMath::Abs(-1+VBFrelUp))* TMath::Abs(-1+VBFrelUp)/(-1+VBFrelUp) ));   // VBF
	  out << space << "-" ; //string(Form("%.2f",1+(-1+ZTTrelDown)*1.0));  // ZTT
	  out << space << "-" ;                                                // QCD
	  out << space << "-" ; //string(Form("%.2f",1+(-1+ZTTrelDown)*1.0));  // W
	  out << space << "-" ; //string(Form("%.2f",1+(-1+ZTTrelDown)*1.0));  // Z+fakes
	  if(binNameSpace.find("SM0")!=string::npos) 
	    out << space << "-" ; //string(Form("%.2f",1+(-1+ZTTrelDown)*1.0));
	  out << space << string(Form("%.2f",1+TMath::Max(  TMath::Abs(-1+TTrelDown),TMath::Abs(-1+TTrelUp))  * TMath::Abs(-1+TTrelUp)/(-1+TTrelUp) ));  // TT
	  if(binNameSpace.find("SM0")!=string::npos) 
	    out << space << string(Form("%.2f",1+TMath::Max(TMath::Abs(-1+VVrelDown),TMath::Abs(-1+VVrelUp) ) * TMath::Abs(-1+VVrelUp)/(-1+VVrelUp)));   // VV
	  else
	    out << space << string(Form("%.2f",1+TMath::Max(TMath::Abs(-1+TTrelDown),TMath::Abs(-1+TTrelUp))  * TMath::Abs(-1+TTrelUp)/(-1+TTrelUp) ));  // TT
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



void produceAll(  TString outputDir = "May2012/Reload_NewCategories" ){

  vector<string> variables;
  vector<int> mH;

  //variables.push_back("diTauVisMass");
  variables.push_back("diTauNSVfitMass");

  //mH.push_back(105);
  //mH.push_back(110);
  mH.push_back(115);
  //mH.push_back(120);
  mH.push_back(125);
  //mH.push_back(130);
  mH.push_back(135);
  mH.push_back(140);
  //mH.push_back(145);

  for(unsigned int i = 0 ; i < variables.size(); i++){
    for(unsigned j = 0; j < mH.size(); j++){

      //produce(mH[j],variables[i], ""        , "vbf", outputDir);
      //produce(mH[j],variables[i], "TauUp"   , "vbf", outputDir);
      //produce(mH[j],variables[i], "TauDown" , "vbf", outputDir);
      //produce(mH[j],variables[i], "JetUp"   , "vbf", outputDir);
      //produce(mH[j],variables[i], "JetDown" , "vbf", outputDir);

      //produce(mH[j],variables[i], ""        , "boost", outputDir);
      //produce(mH[j],variables[i], "TauUp"   , "boost", outputDir);
      //produce(mH[j],variables[i], "TauDown" , "boost", outputDir);
      //produce(mH[j],variables[i], "JetUp"   , "boost", outputDir);
      //produce(mH[j],variables[i], "JetDown" , "boost", outputDir);
    
      produce(mH[j],variables[i], ""        , "novbf", outputDir);
      produce(mH[j],variables[i], "TauUp"   , "novbf", outputDir);
      produce(mH[j],variables[i], "TauDown" , "novbf", outputDir);
      produce(mH[j],variables[i], "JetUp"   , "novbf", outputDir);
      produce(mH[j],variables[i], "JetDown" , "novbf", outputDir);

      //produce(mH[j],variables[i], ""        , "twoJets", outputDir);
      //produce(mH[j],variables[i], "TauUp"   , "twoJets", outputDir);
      //produce(mH[j],variables[i], "TauDown" , "twoJets", outputDir);
      //produce(mH[j],variables[i], "JetUp"   , "twoJets", outputDir);
      //produce(mH[j],variables[i], "JetDown" , "twoJets", outputDir);
      
      //produce(mH[j],variables[i], ""        , "oneJet", outputDir);
      //produce(mH[j],variables[i], "TauUp"   , "oneJet", outputDir);
      //produce(mH[j],variables[i], "TauDown" , "oneJet", outputDir);
      //produce(mH[j],variables[i], "JetUp"   , "oneJet", outputDir);
      //produce(mH[j],variables[i], "JetDown" , "oneJet", outputDir);
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

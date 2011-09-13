
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TDirectory.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>


void produce(   
	     int mH_=120,
	     string variable_ = "diTauVisMass",
	     string analysis_ = "",
	     string bin_= "inclusive"
	     ){

  cout << "Now doing mass mH=" << mH_ << ", for variable " << variable_ << " analysis " << analysis_ << " and bin " << bin_ << endl;


  string binNameSpace = "";
  if(bin_.find("inclusive")!=string::npos)      
    binNameSpace =  "inclusive";
  else if(bin_.find("novbf")!=string::npos) 
    binNameSpace =  "SM0";
  else if(bin_.find("vbf")!=string::npos && bin_.find("novbf")==string::npos) 
    binNameSpace =  "SM2";

  ifstream in;

  char* c = new char[1000];
  in.open(Form("templates/muTau_%s_template.txt",  binNameSpace.c_str()));
  ofstream out(Form("muTau_%s_mH%d.txt", binNameSpace.c_str(), mH_));
  out.precision(8);

  TFile* fin = new TFile(Form("histograms/muTau_mH%d_%s_%s_%s.root", mH_, bin_.c_str() , analysis_.c_str(), variable_.c_str()), "READ");

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
          line.replace( line.find("muTauSM") , 7 , string(Form("muTauSM_%s",variable_.c_str()))  );
	  out << line << endl;
	}
	else if(line.find("process")!=string::npos && line.find("VBF")!=string::npos){
          line.replace( line.find("XXX") , 3 , string(Form("%d",mH_))  );
          line.replace( line.find("YYY") , 3 , string(Form("%d",mH_))  );
          out << line << endl;
        }
	else if(line.find("rate")!=string::npos){
	  string rate = "rate                      ";
	  string space = "                ";
	  out << rate ;
	  out << space << ((TH1F*)fin->Get("hSgn2"))->Integral() 
	      << space << ((TH1F*)fin->Get("hSgn1"))->Integral()
	      << space << ((TH1F*)fin->Get("hZtt"))->Integral()
	      << space << ((TH1F*)fin->Get("hQCD"))->Integral()
	      << space << ((TH1F*)fin->Get("hW"))->Integral()
	      << space << ((TH1F*)fin->Get("hZmj"))->Integral()
	      << space << ((TH1F*)fin->Get("hZmm"))->Integral()
	      << space << ((TH1F*)fin->Get("hTTb"))->Integral()
	      << space << ((TH1F*)fin->Get("hVV"))->Integral()
	      << endl;
	}
	else if(line.find(Form("CMS_htt_muTau_%s_QCDNorm", binNameSpace.c_str()))!=string::npos){
	  line.replace( line.find("XXX") , 3 , string(Form("%.0f",((TH1F*)fin->Get("hParameters"))->GetBinContent(5)))  );
	  line.replace( line.find("YYY") , 3 , string(Form("%.2f",((TH1F*)fin->Get("hParameters"))->GetBinContent(6)))  );
	  out << line << endl;
	}
	else if(line.find(Form("CMS_htt_muTau_%s_WNorm",binNameSpace.c_str()))!=string::npos && binNameSpace!="SM2"){
          line.replace( line.find("XXX") , 3 , string(Form("%.0f",((TH1F*)fin->Get("hParameters"))->GetBinContent(2)))  );
          line.replace( line.find("YYY") , 3 , string(Form("%.2f",((TH1F*)fin->Get("hParameters"))->GetBinContent(1)))  );
          out << line << endl;
        }
	else{
	  out << c << endl;
	}
      
      }
 
   }

  TFile* fTemplOut = new TFile(Form("muTauSM_%s.root",variable_.c_str()),"UPDATE");

  string suffix = "";
  if(analysis_.find("TauUp")!=string::npos)
    suffix = "CMS_scale_tUp";
  else if(analysis_.find("TauDown")!=string::npos)
    suffix = "CMS_scale_tDown";
  else if(analysis_.find("MuUp")!=string::npos)
    suffix = "CMS_scale_mUp";
  else if(analysis_.find("MuDown")!=string::npos)
    suffix = "CMS_scale_mDown";
  else if(analysis_.find("JetUp")!=string::npos)
    suffix = "CMS_scale_jUp";
  else if(analysis_.find("JetDown")!=string::npos)
    suffix = "CMS_scale_jDown";
  
  cout << "Adding histos with suffix " << suffix << endl;
  
  if(! (fTemplOut->cd(Form("muTau%s_%s",binNameSpace.c_str(),variable_.c_str())))  ){
    
    cout << "Editing the directory for bin and variable " << binNameSpace << ", " << variable_ << endl;

    TDirectory* dir = fTemplOut->mkdir( Form("muTau%s_%s",binNameSpace.c_str(),variable_.c_str()) );
    dir->cd();

    ((TH1F*)fin->Get("hSgn1"))->Write(Form("VBF%d%s",mH_,suffix.c_str()));
    ((TH1F*)fin->Get("hSgn2"))->Write(Form("SM%d%s" ,mH_,suffix.c_str()));
    ((TH1F*)fin->Get("hZtt"))->Write(Form("ZTT%s"       ,suffix.c_str()));
    ((TH1F*)fin->Get("hQCD"))->Write(Form("QCD%s"       ,suffix.c_str()));
    ((TH1F*)fin->Get("hW"))->Write(Form("W%s"           ,suffix.c_str()));
    ((TH1F*)fin->Get("hZmj"))->Write(Form("ZJ%s"        ,suffix.c_str()));
    ((TH1F*)fin->Get("hZmm"))->Write(Form("ZL%s"        ,suffix.c_str()));
    ((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
    ((TH1F*)fin->Get("hVV"))->Write(Form("VV%s"         ,suffix.c_str()));
           
  }else{
    cout << "Directory is there, filling only higgs histos..." << endl;
    //write the higgs only
    fTemplOut->cd(Form("muTau%s_%s",binNameSpace.c_str(),variable_.c_str()));

    if(fTemplOut->FindObjectAny(Form("ZTT%s"       ,suffix.c_str()))==0 )
      ((TH1F*)fin->Get("hZtt"))->Write(Form("ZTT%s"       ,suffix.c_str()));
    if(fTemplOut->FindObjectAny(Form("QCD%s"       ,suffix.c_str()))==0 )
      ((TH1F*)fin->Get("hQCD"))->Write(Form("QCD%s"       ,suffix.c_str()));
    if(fTemplOut->FindObjectAny(Form("W%s"       ,suffix.c_str()))==0 )
      ((TH1F*)fin->Get("hW"))->Write(Form("W%s"           ,suffix.c_str()));
    if(fTemplOut->FindObjectAny(Form("ZJ%s"       ,suffix.c_str()))==0 )
      ((TH1F*)fin->Get("hZmj"))->Write(Form("ZJ%s"        ,suffix.c_str()));
    if(fTemplOut->FindObjectAny(Form("ZL%s"       ,suffix.c_str()))==0 ) 
      ((TH1F*)fin->Get("hZmm"))->Write(Form("ZL%s"        ,suffix.c_str()));
    if(fTemplOut->FindObjectAny(Form("TT%s"       ,suffix.c_str()))==0 )
      ((TH1F*)fin->Get("hTTb"))->Write(Form("TT%s"        ,suffix.c_str()));
    if(fTemplOut->FindObjectAny(Form("VV%s"       ,suffix.c_str()))==0 )
      ((TH1F*)fin->Get("hVV"))->Write(Form("VV%s"         ,suffix.c_str()));
    if(fTemplOut->FindObjectAny(Form("VBF%d%s"         ,mH_,suffix.c_str()))==0 )
      ((TH1F*)fin->Get("hSgn1"))->Write(Form("VBF%d%s",mH_  ,suffix.c_str()));
    if(fTemplOut->FindObjectAny(Form("SM%d%s"          , mH_,suffix.c_str()))==0 )
      ((TH1F*)fin->Get("hSgn2"))->Write(Form("SM%d%s" ,mH_,suffix.c_str()));
  }

  //fTemplOut->Write();
  fTemplOut->Close();

  return;

}



void produceAll(){

  produce(120,"diTauVisMass", ""        , "vbf");
  produce(120,"diTauVisMass", "TauUp"   , "vbf");
  produce(120,"diTauVisMass", "TauDown" , "vbf");

  produce(120,"diTauVisMass", ""        , "novbf");
  produce(120,"diTauVisMass", "TauUp"   , "novbf");
  produce(120,"diTauVisMass", "TauDown" , "novbf");

  produce(130,"diTauVisMass", ""        , "vbf");
  produce(130,"diTauVisMass", "TauUp"   , "vbf");
  produce(130,"diTauVisMass", "TauDown" , "vbf");

  produce(130,"diTauVisMass", ""        , "novbf");
  produce(130,"diTauVisMass", "TauUp"   , "novbf");
  produce(130,"diTauVisMass", "TauDown" , "novbf");


}

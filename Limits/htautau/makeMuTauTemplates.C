
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>



void fun(){

  int mH_=120;

  ifstream in;

  char* c = new char[1000];
  in.open("muTau_SM0_template.txt");

  ofstream out(Form("muTau_SM0_mH%d_v2.txt", mH_));
  out.precision(8);

  TFile* fin = new TFile("muTau_mH120_inclusive__diTauVisMass.root");


  while (in.good())
    {
      in.getline(c,1000,'\n');
      if (in.good()){

	string line(c);
	if(line.find("rate")!=string::npos){
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
	else if(line.find("CMS_htt_muTau_QCDNorm")!=string::npos){
	  line.replace( line.find("XXX") , 3 , string(Form("%.0f",((TH1F*)fin->Get("hParameters"))->GetBinContent(5)))  );
	  line.replace( line.find("YYY") , 3 , string(Form("%.2f",((TH1F*)fin->Get("hParameters"))->GetBinContent(6)))  );
	  out << line << endl;
	}
	else if(line.find("CMS_htt_muTau_SM0_WNorm")!=string::npos){
          line.replace( line.find("XXX") , 3 , string(Form("%.0f",((TH1F*)fin->Get("hParameters"))->GetBinContent(2)))  );
          line.replace( line.find("YYY") , 3 , string(Form("%.2f",((TH1F*)fin->Get("hParameters"))->GetBinContent(1)))  );
          out << line << endl;
        }
	else{
	  out << c << endl;
	}
      
      }
 
   }

}





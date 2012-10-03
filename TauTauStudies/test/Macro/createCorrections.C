#include "TFile.h"
#include "TF1.h"
#include "TSystem.h"

//#include "ratioEfficiencyElec.C"
#include "ratioEfficiencyTest.C"


Double_t myFuncRatioElec20BL(Double_t* x, Double_t *par) {  //modified
  gSystem->Load("ratioEfficiencyTest_C.so");  
  ratioEfficiencyTest* fitEffElec20BL = new ratioEfficiencyTest(20.97643939, 1.15196354, 2.27544602, 1.01743868, 2.04391816);
  ratioEfficiencyTest* fitEffElec20BLMC = new ratioEfficiencyTest(20.58604584, -1.89456806, 3.69311772, 1.05480046, 1.28655181);  
  
  Float_t xx = x[0];  
  return fitEffElec20BL->turnOn(xx)/fitEffElec20BLMC->turnOn(xx);  

}

Double_t myFuncTurnOnElec20BL(Double_t* x, Double_t *par) {  //modified

  gSystem->Load("ratioEfficiencyTest_C.so");   
  ratioEfficiencyTest* fitEffElec20BL = new ratioEfficiencyTest(20.97643939, 1.15196354, 2.27544602, 1.01743868, 2.04391816);
  
  Float_t xx = x[0];   
  return fitEffElec20BL->turnOn(xx);
}

Double_t myFuncRatioElec20EC(Double_t* x, Double_t *par) { //modified

  gSystem->Load("ratioEfficiencyTest_C.so");   
  ratioEfficiencyTest* fitEffElec20EC = new ratioEfficiencyTest(20.59874300, 1.25425435, 1.61098921, 1.00146962, 60.35067579);
  ratioEfficiencyTest* fitEffElec20ECMC = new ratioEfficiencyTest(20.15425918, 0.75449122, 1.06027513, 1.01106686, 7.01956561);
   
  Float_t xx = x[0];   
  return fitEffElec20EC->turnOn(xx)/fitEffElec20ECMC->turnOn(xx);   
 
}

Double_t myFuncTurnOnElec20EC(Double_t* x, Double_t *par) { //modified
  
  gSystem->Load("ratioEfficiencyTest_C.so");    
  ratioEfficiencyTest* fitEffElec20EC = new ratioEfficiencyTest(20.59874300, 1.25425435, 1.61098921, 1.00146962, 60.35067579);

  Float_t xx = x[0];    
  return fitEffElec20EC->turnOn(xx);
}

Double_t myFuncRatioEle22BL(Double_t* x, Double_t *par) {  //modified

  gSystem->Load("ratioEfficiencyTest_C.so");   
  ratioEfficiencyTest* fitEffElec22BL = new ratioEfficiencyTest(22.90752344, 1.32376429, 2.17813319, 1.03674051, 2.15454768);
  ratioEfficiencyTest* fitEffElec20BLMC = new ratioEfficiencyTest(20.58604584, -1.89456806, 3.69311772, 1.05480046, 1.28655181);   
  
  Float_t xx = x[0];   
  return fitEffElec22BL->turnOn(xx)/fitEffElec20BLMC->turnOn(xx);   
 
}

Double_t myFuncTurnOnEle22BL(Double_t* x, Double_t *par) { //modified

  gSystem->Load("ratioEfficiencyTest_C.so");    
  ratioEfficiencyTest* fitEffElec22BL = new ratioEfficiencyTest(22.90752344, 1.32376429, 2.17813319, 1.03674051, 2.15454768);

  Float_t xx = x[0];    
  return fitEffElec22BL->turnOn(xx);
}


Double_t myFuncRatioEle22EC(Double_t* x, Double_t *par) { //modified

  gSystem->Load("ratioEfficiencyTest_C.so");    
  ratioEfficiencyTest* fitEffElec22EC = new ratioEfficiencyTest(22.14553261, 1.19913124, 1.75642067, 1.00826962, 9.04331617);
  ratioEfficiencyTest* fitEffElec20ECMC = new ratioEfficiencyTest(20.15425918, 0.75449122, 1.06027513, 1.01106686, 7.01956561);

   
  Float_t xx = x[0];    
  return fitEffElec22EC->turnOn(xx)/fitEffElec20ECMC->turnOn(xx);    
  
}


Double_t myFuncTurnOnEle22EC(Double_t* x, Double_t *par) { //modified

  gSystem->Load("ratioEfficiencyTest_C.so");     
  ratioEfficiencyTest* fitEffElec22EC = new ratioEfficiencyTest(22.14553261, 1.19913124, 1.75642067, 1.00826962, 9.04331617);

  Float_t xx = x[0];     
  return fitEffElec22EC->turnOn(xx);
}

Double_t myFuncRatioEleAllBL(Double_t* x, Double_t *par) { //modified

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle20BL = new ratioEfficiencyTest(20.97643939, 1.15196354, 2.27544602, 1.01743868, 2.04391816);
  ratioEfficiencyTest* fitEffEle22BL = new ratioEfficiencyTest(22.90752344, 1.32376429, 2.17813319, 1.03674051, 2.15454768);

  ratioEfficiencyTest* fitEffElec20MCBL = new ratioEfficiencyTest(20.58604584, -1.89456806, 3.69311772, 1.05480046, 1.28655181);

  float weightElec20 = 0.14; //change
  float weightElec22 = 0.86; //change


  float total = weightElec20+weightElec22;
  weightElec20/=total;
  weightElec22/=total;


  Float_t xx = x[0];

  return (fitEffEle20BL->turnOn(xx) * weightElec20 +
	  fitEffEle22BL->turnOn(xx) * weightElec22 
	  )/fitEffElec20MCBL->turnOn(xx);
}

Double_t myFuncTurnOnEleAllBL(Double_t* x, Double_t *par) { //modified

  gSystem->Load("ratioEfficiencyTest_C.so"); 
  ratioEfficiencyTest* fitEffEle20BL = new ratioEfficiencyTest(20.97643939, 1.15196354, 2.27544602, 1.01743868, 2.04391816); 
  ratioEfficiencyTest* fitEffEle22BL = new ratioEfficiencyTest(22.90752344, 1.32376429, 2.17813319, 1.03674051, 2.15454768); 
 
  float weightElec20 = 0.14; //change 
  float weightElec22 = 0.86; //change 
 
 
  float total = weightElec20+weightElec22; 
  weightElec20/=total; 
  weightElec22/=total; 
 
  Float_t xx = x[0];

  return (fitEffEle20BL->turnOn(xx) * weightElec20 +
	  fitEffEle22BL->turnOn(xx) * weightElec22 
	  );
 
}



Double_t myFuncRatioEleAllEC(Double_t* x, Double_t *par) { //modified

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle20EC = new ratioEfficiencyTest(20.59874300, 1.25425435, 1.61098921, 1.00146962, 60.35067579);
  ratioEfficiencyTest* fitEffEle22EC = new ratioEfficiencyTest(22.14553261, 1.19913124, 1.75642067, 1.00826962, 9.04331617);

  ratioEfficiencyTest* fitEffElec20MCEC = new ratioEfficiencyTest(20.15425918, 0.75449122, 1.06027513, 1.01106686, 7.01956561);

  float weightElec20 = 0.14;  //change
  float weightElec22 = 0.86;  //change

  float total = weightElec20+weightElec22;
  weightElec20/=total;
  weightElec22/=total;

  Float_t xx = x[0];

  return (fitEffEle20EC->turnOn(xx) * weightElec20 +
	  fitEffEle22EC->turnOn(xx) * weightElec22 
	  )/fitEffElec20MCEC->turnOn(xx);
}

Double_t myFuncTurnOnEleAllEC(Double_t* x, Double_t *par) { //modified

  gSystem->Load("ratioEfficiencyTest_C.so"); 
  ratioEfficiencyTest* fitEffEle20EC = new ratioEfficiencyTest(20.59874300, 1.25425435, 1.61098921, 1.00146962, 60.35067579); 
  ratioEfficiencyTest* fitEffEle22EC = new ratioEfficiencyTest(22.14553261, 1.19913124, 1.75642067, 1.00826962, 9.04331617); 
 
  float weightElec20 = 0.14;  //change 
  float weightElec22 = 0.86;  //change 
 
  float total = weightElec20+weightElec22; 
  weightElec20/=total; 
  weightElec22/=total;
  
  Float_t xx = x[0];

  return (fitEffEle20EC->turnOn(xx) * weightElec20 +
	  fitEffEle22EC->turnOn(xx) * weightElec22
	  );
 
}


///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
Double_t myFuncRatioElecIDBL(Double_t* x, Double_t *par) {  //modified
  double ratio = 1.; 
 
  Float_t xx = x[0]; 
   
  if( xx < 20 )  
    ratio = 0.922; 
  else if( xx >= 20 && xx< 30)  
    ratio = 0.922  ; 
  else if( xx >= 30)
    ratio = 0.964 ; 
  
  return ratio; 
  
}
Double_t myFuncTurnOnElecIDBL(Double_t* x, Double_t *par) {  //modified
  double ratio = 1.;  
  
  Float_t xx = x[0];  
    
  if( xx < 20 )   
    ratio = 0.733;  
  else if( xx >= 20 && xx< 30)   
    ratio = 0.733 ;  
  else if( xx >= 30) 
    ratio = 0.876 ;  
   
  return ratio;  
  
}
Double_t myFuncRatioElecIDEC(Double_t* x, Double_t *par) { //modified
  double ratio = 1.;   
   
  Float_t xx = x[0];   
     
  if( xx < 20 )    
    ratio = 0.898 ;   
  else if( xx >= 20 && xx< 30)    
    ratio = 0.898 ;   
  else if( xx >= 30)  
    ratio = 0.958 ;   
    
  return ratio;   

}
Double_t myFuncTurnOnElecIDEC(Double_t* x, Double_t *par) { //modified
  double ratio = 1.;    
    
  Float_t xx = x[0];    
      
  if( xx < 20 )     
    ratio = 0.389;    
  else if( xx >= 20 && xx< 30)     
    ratio = 0.389 ;    
  else if( xx >= 30)   
    ratio = 0.592 ;    
     
  return ratio;    

}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

Double_t myFuncRatioElecIsoBL(Double_t* x, Double_t *par) { //modified
  
  double ratio = 1.;

  Float_t xx = x[0];
  
  if( xx < 20 ) 
    ratio = 0.974 ;
  else if( xx >= 20 && xx< 30) 
    ratio = 0.974 ;
  else if( xx >= 30)
    ratio = 0.997 ;
  
  return ratio;
}

Double_t myFuncTurnOnElecIsoBL(Double_t* x, Double_t *par) { //modified

  double ratio = 1.;

  Float_t xx = x[0];
  
 if( xx < 20 ) 
    ratio = 0.715;
  else if( xx >= 20 && xx< 30) 
    ratio = 0.715 ;
  else if( xx >= 30)
    ratio = 0.893 ;

  return ratio;

}

Double_t myFuncRatioElecIsoEC(Double_t* x, Double_t *par) { //modified

  double ratio = 1.;

  Float_t xx = x[0];

  if( xx < 20 )  
    ratio = 1.008 ; 
  else if( xx >= 20 && xx< 30)  
    ratio = 1.008  ; 
  else if( xx >= 30) 
    ratio = 0.983  ; 

  return ratio;
}

Double_t myFuncTurnOnElecIsoEC(Double_t* x, Double_t *par) {  //modified
 
  double ratio = 1.;

  Float_t xx = x[0];
 
  if( xx < 20 )   
    ratio = 0.745 ;  
  else if( xx >= 20 && xx< 30)   
    ratio = 0.745 ;  
  else if( xx >= 30)  
    ratio = 0.896 ;  


  return ratio;

}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

Double_t myFuncRatioMuIDBL(Double_t* x, Double_t *par) { //modified
  double ratio = 1.; 
 
  Float_t xx = x[0]; 
   
  if( xx < 15 )  
    ratio = 0.989 ; 
  else if( xx >= 15 && xx< 20)  
    ratio = 0.989 ; 
  else if( xx >= 20 && xx< 30)  
    ratio = 0.991 ; 
  else if( xx >= 30)
    ratio = 0.989 ; 
 
  return ratio; 

}
Double_t myFuncTurnOnMuIDBL(Double_t* x, Double_t *par) { //modified
  double ratio = 1.;  
  
  Float_t xx = x[0];  
    
  if( xx < 15 )   
    ratio = 0.968  ;  
  else if( xx >= 15 && xx< 20)   
    ratio = 0.968  ;  
  else if( xx >= 20 && xx< 30)   
    ratio = 0.961  ;  
  else if( xx >= 30) 
    ratio = 0.962  ;  
  
  return ratio;  

}
Double_t myFuncRatioMuIDEC(Double_t* x, Double_t *par) {  //modified
  double ratio = 1.;   
   
  Float_t xx = x[0];   
     
  if( xx < 15 )    
    ratio = 0.977  ;   
  else if( xx >= 15 && xx< 20)    
    ratio = 0.977  ;   
  else if( xx >= 20 && xx< 30)    
    ratio = 0.974  ;   
  else if( xx >= 30)  
    ratio = 0.989  ;   
   
  return ratio;   

}
Double_t myFuncTurnOnMuIDEC(Double_t* x, Double_t *par) {  //modified
  double ratio = 1.;    
    
  Float_t xx = x[0];    
      
  if( xx < 15 )     
    ratio = 0.959  ;    
  else if( xx >= 15 && xx< 20)     
    ratio = 0.959  ;    
  else if( xx >= 20 && xx< 30)     
    ratio = 0.951  ;    
  else if( xx >= 30)   
    ratio = 0.955  ;    
    
  return ratio;    
 
}

///////////////////////////////////////////////////////


Double_t myFuncRatioMuIsoBL(Double_t* x, Double_t *par) { //modified
  
  double ratio = 1.;     
     
  Float_t xx = x[0];     
       
  if( xx < 15 )      
    ratio = 0.945  ;     
  else if( xx >= 15 && xx< 20)      
    ratio = 0.945  ;     
  else if( xx >= 20 && xx< 30)      
    ratio = 1.005  ;     
  else if( xx >= 30)    
    ratio = 0.993  ;     
     
  return ratio;     

}

Double_t myFuncTurnOnMuIsoBL(Double_t* x, Double_t *par) {  //modified

     
  double ratio = 1.;      
      
  Float_t xx = x[0];      
        
  if( xx < 15 )       
    ratio = 0.688  ;      
  else if( xx >= 15 && xx< 20)       
    ratio = 0.688  ;      
  else if( xx >= 20 && xx< 30)       
    ratio = 0.775  ;      
  else if( xx >= 30)     
    ratio = 0.915 ;      
      
  return ratio;      

}

Double_t myFuncRatioMuIsoEC(Double_t* x, Double_t *par) {  //modified

  double ratio = 1.;       
       
  Float_t xx = x[0];       
         
  if( xx < 15 )        
    ratio = 1.047   ;       
  else if( xx >= 15 && xx< 20)        
    ratio = 1.047   ;       
  else if( xx >= 20 && xx< 30)        
    ratio = 0.992   ;       
  else if( xx >= 30)      
    ratio = 1.005  ;       
       
  return ratio;       

}

Double_t myFuncTurnOnMuIsoEC(Double_t* x, Double_t *par) { //modified
  double ratio = 1.;        
        
  Float_t xx = x[0];        
          
  if( xx < 15 )         
    ratio = 0.750    ;        
  else if( xx >= 15 && xx< 20)         
    ratio = 0.750    ;        
  else if( xx >= 20 && xx< 30)         
    ratio = 0.834    ;        
  else if( xx >= 30)       
    ratio = 0.935   ;        
        
  return ratio;        
 
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

Double_t myFuncRatioMu18BL(Double_t* x, Double_t *par) { //modified
  gSystem->Load("ratioEfficiencyTest_C.so"); 
  ratioEfficiencyTest* fitEffMuBL = new ratioEfficiencyTest(15.99983195, -0.39072829, 0.28256338, 1.72861719, 0.95769408); 
  ratioEfficiencyTest* fitEffMuBLMC = new ratioEfficiencyTest(16.99389526, -0.04080190, 0.00794730, 1.60377906, 0.99626161); 
 
  Float_t xx = x[0]; 
  return fitEffMuBL->turnOn(xx)/fitEffMuBLMC->turnOn(xx); 

}

Double_t myFuncTurnOnMu18BL(Double_t* x, Double_t *par) { //modified
  gSystem->Load("ratioEfficiencyTest_C.so");  
  ratioEfficiencyTest* fitEffMuBL = new ratioEfficiencyTest(15.99983195, -0.39072829, 0.28256338, 1.72861719, 0.95769408);
  Float_t xx = x[0];  
  return fitEffMuBL->turnOn(xx);
}

Double_t myFuncRatioMu18EC(Double_t* x, Double_t *par) { //modified
  gSystem->Load("ratioEfficiencyTest_C.so");  
  ratioEfficiencyTest* fitEffMuEC = new ratioEfficiencyTest(18.49754887, -0.16941614, 0.26076717, 1.05494469, 1.53819978);
  ratioEfficiencyTest* fitEffMuECMC = new ratioEfficiencyTest(16.99065795, -0.11993730, 0.01384991, 2.38867304, 0.86552275);
  
  Float_t xx = x[0];  
  return fitEffMuEC->turnOn(xx)/fitEffMuECMC->turnOn(xx);
  
}

Double_t myFuncTurnOnMu18EC(Double_t* x, Double_t *par) { //modified
  gSystem->Load("ratioEfficiencyTest_C.so");   
  ratioEfficiencyTest* fitEffMuEC = new ratioEfficiencyTest(18.49754887, -0.16941614, 0.26076717, 1.05494469, 1.53819978); 

  Float_t xx = x[0];   
  return fitEffMuEC->turnOn(xx);
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

Double_t myFuncRatioMu17BL(Double_t* x, Double_t *par) { //modified
  gSystem->Load("ratioEfficiencyTest_C.so");  
  ratioEfficiencyTest* fitEffMu17BL = new ratioEfficiencyTest(17.21270264, 0.54997112, 1.02874912, 1.29646487, 0.96724273);
  ratioEfficiencyTest* fitEffMu18BLMC = new ratioEfficiencyTest(16.99389526, -0.04080190, 0.00794730, 1.60377906, 0.99626161);  
  
  Float_t xx = x[0];  
  return fitEffMu17BL->turnOn(xx)/fitEffMu18BLMC->turnOn(xx);  
}

Double_t myFuncTurnOnMu17BL(Double_t* x, Double_t *par) {  //modified

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffMu17BL = new ratioEfficiencyTest(17.21270264, 0.54997112, 1.02874912, 1.29646487, 0.96724273);

  Float_t xx = x[0];
  return fitEffMu17BL->turnOn(xx);
}

Double_t myFuncRatioMu17EC(Double_t* x, Double_t *par) { //modified
  gSystem->Load("ratioEfficiencyTest_C.so");   
  ratioEfficiencyTest* fitEffMu17EC = new ratioEfficiencyTest(15.98037640, 0.12062946, 0.02183977, 2.84751010, 0.83985656);
  ratioEfficiencyTest* fitEffMu18ECMC = new ratioEfficiencyTest(16.99065795, -0.11993730, 0.01384991, 2.38867304, 0.86552275);
   
  Float_t xx = x[0];   
  return fitEffMu17EC->turnOn(xx)/fitEffMu18ECMC->turnOn(xx);   
   
}

Double_t myFuncTurnOnMu17EC(Double_t* x, Double_t *par) { //modified 

  gSystem->Load("ratioEfficiencyTest_C.so");    
  ratioEfficiencyTest* fitEffMu17EC = new ratioEfficiencyTest(15.98037640, 0.12062946, 0.02183977, 2.84751010, 0.83985656); 
  Float_t xx = x[0];    
  return fitEffMu17EC->turnOn(xx);
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

Double_t myFuncRatioMuAllBL(Double_t* x, Double_t *par) { //modified

  TF1* ratioMu18BL = new TF1("ratioMu18BL", myFuncRatioMu18BL,0,400,0);
  TF1* ratioMu17BL = new TF1("ratioMu17BL", myFuncRatioMu17BL,0,400,0);
  
  float weightIsoMu18 = 0.14;  //change
  float weightIsoMu17 = 0.86;  //change

  float total = weightIsoMu18+weightIsoMu17;
  weightIsoMu18 /= total;
  weightIsoMu17 /= total;

  Float_t xx = x[0];

  return (ratioMu18BL->Eval(xx) * weightIsoMu18 +
	  ratioMu17BL->Eval(xx) * weightIsoMu17 );

}

Double_t myFuncTurnOnMuAllBL(Double_t* x, Double_t *par) { //modified

  TF1* turnOnMu18BL = new TF1("turnOnMu18BL", myFuncTurnOnMu18BL,0,400,0);
  TF1* turnOnMu17BL = new TF1("turnOnMu17BL", myFuncTurnOnMu17BL,0,400,0);

  float weightIsoMu18 = 0.14;  //change
  float weightIsoMu17 = 0.86; //change

  float total = weightIsoMu18+weightIsoMu17;
  weightIsoMu18 /=total;
  weightIsoMu17 /=total;

  Float_t xx = x[0];
  
  return (turnOnMu18BL->Eval(xx) * weightIsoMu18 +
	  turnOnMu17BL->Eval(xx) * weightIsoMu17 );
}


Double_t myFuncRatioMuAllEC(Double_t* x, Double_t *par) { //modified

  TF1* ratioMu18EC = new TF1("ratioMu18EC", myFuncRatioMu18EC,0,400,0);
  TF1* ratioMu17EC = new TF1("ratioMu17EC", myFuncRatioMu17EC,0,400,0);
  
  float weightIsoMu18 = 0.14;  //change
  float weightIsoMu17 = 0.86;  //change

  float total = weightIsoMu18+weightIsoMu17;
  weightIsoMu18 /= total;
  weightIsoMu17 /= total;

  Float_t xx = x[0];

  return (ratioMu18EC->Eval(xx) * weightIsoMu18 +
	  ratioMu17EC->Eval(xx) * weightIsoMu17 );

}

Double_t myFuncTurnOnMuAllEC(Double_t* x, Double_t *par) {  //modified

  TF1* turnOnMu18EC = new TF1("turnOnMu18EC", myFuncTurnOnMu18EC,0,400,0);
  TF1* turnOnMu17EC = new TF1("turnOnMu17EC", myFuncTurnOnMu17EC,0,400,0);

  float weightIsoMu18 = 0.14;  //change
  float weightIsoMu17 = 0.86;  //change

  float total = weightIsoMu18+weightIsoMu17;
  weightIsoMu18 /=total;
  weightIsoMu17 /=total;

  Float_t xx = x[0];

  return (turnOnMu18EC->Eval(xx) * weightIsoMu18 +
	  turnOnMu17EC->Eval(xx) * weightIsoMu17 );
}






/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

// HLT PFTauLoose20 / MC Summer12 / mu+tau

Double_t myFuncTurnOnTauLoose20MuTauMCBL(Double_t* x, Double_t *par) {  //modified
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MC = new ratioEfficiencyTest(18.86257072, 0.25680380, 0.16916101, 2.42931257, 0.89590264);
  
  Float_t xx = x[0];

  return ratioEffTauLoose20MC->turnOn(xx);

}

Double_t myFuncTurnOnTauLoose20MuTauMCEC(Double_t* x, Double_t *par) {  //modified
   
  gSystem->Load("ratioEfficiencyTest_C.so"); 
  ratioEfficiencyTest* ratioEffTauLoose20MC = new ratioEfficiencyTest(18.74764561, 1.82036845, 701.46994969, 101.57913480, 0.82547043);
  
  Float_t xx = x[0]; 
  
  return ratioEffTauLoose20MC->turnOn(xx); 
 
} 

/////////////////////////////////////////////////


/////////////////////////////////////////////////
/////////////////////////////////////////////////

// HLT PFTauLoose20 / Run2012A (Run<=193686)    / mu+tau

Double_t myFuncRatioTauLoose20MuTauRunABL(Double_t* x, Double_t *par) { //modified
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunABL = new ratioEfficiencyTest(18.52262128, 1.85879597, 3.48843815, 1.15491294, 1.02489024);
  ratioEfficiencyTest* ratioEffTauLoose20MuTauMCBL = new ratioEfficiencyTest(18.86257072, 0.25680380, 0.16916101, 2.42931257, 0.89590264);

  Float_t xx = x[0];

  return ratioEffTauLoose20MuTauRunABL->turnOn(xx)/ratioEffTauLoose20MuTauMCBL->turnOn(xx);

}

Double_t myFuncTurnOnTauLoose20MuTauRunABL(Double_t* x, Double_t *par) { //modified
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunABL = new ratioEfficiencyTest(18.52262128, 1.85879597, 3.48843815, 1.15491294, 1.02489024);
 
  Float_t xx = x[0];

  return ratioEffTauLoose20MuTauRunABL->turnOn(xx);

}

Double_t myFuncRatioTauLoose20MuTauRunAEC(Double_t* x, Double_t *par) { //modified 
   
  gSystem->Load("ratioEfficiencyTest_C.so"); 
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunAEC = new ratioEfficiencyTest(18.90119559, 0.14025596, 0.14482632, 1.56126508, 0.81188198);
  ratioEfficiencyTest* ratioEffTauLoose20MuTauMCEC = new ratioEfficiencyTest(18.74764561, 1.82036845, 701.46994969, 101.57913480, 0.82547043);
 
  Float_t xx = x[0]; 
 
  return ratioEffTauLoose20MuTauRunAEC->turnOn(xx)/ratioEffTauLoose20MuTauMCEC->turnOn(xx); 
 
} 

Double_t myFuncTurnOnTauLoose20MuTauRunAEC(Double_t* x, Double_t *par) { //modified

  gSystem->Load("ratioEfficiencyTest_C.so");  
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunAEC = new ratioEfficiencyTest(18.90119559, 0.14025596, 0.14482632, 1.56126508, 0.81188198);

  Float_t xx = x[0];  
  
  return ratioEffTauLoose20MuTauRunAEC->turnOn(xx);
}

/////////////////////////////////////////////////

// HLT PFTauLoose20 / Run2012B    / mu+tau

Double_t myFuncRatioTauLoose20MuTauRunBBL(Double_t* x, Double_t *par) { //modified
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunBBL = new ratioEfficiencyTest(17.92648563, 1.96846742, 4.46406075, 1.02023992, 1.52260575);
  ratioEfficiencyTest* ratioEffTauLoose20MuTauMCBL = new ratioEfficiencyTest(18.86257072, 0.25680380, 0.16916101, 2.42931257, 0.89590264);


  Float_t xx = x[0];

  return ratioEffTauLoose20MuTauRunBBL->turnOn(xx)/ ratioEffTauLoose20MuTauMCBL->turnOn(xx);

}

Double_t myFuncTurnOnTauLoose20MuTauRunBBL(Double_t* x, Double_t *par) { //modified
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunBBL = new ratioEfficiencyTest(17.92648563, 1.96846742, 4.46406075, 1.02023992, 1.52260575);

  Float_t xx = x[0];

  return ratioEffTauLoose20MuTauRunBBL->turnOn(xx);

}
 
Double_t myFuncRatioTauLoose20MuTauRunBEC(Double_t* x, Double_t *par) { //modified  
    
  gSystem->Load("ratioEfficiencyTest_C.so");  
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunBEC = new ratioEfficiencyTest(18.59856420, 2.49132550, 10.99643595, 1.50651123, 0.87952970);
  ratioEfficiencyTest* ratioEffTauLoose20MuTauMCEC = new ratioEfficiencyTest(18.74764561, 1.82036845, 701.46994969, 101.57913480, 0.82547043); 
  
  Float_t xx = x[0];  
  
  return ratioEffTauLoose20MuTauRunBEC->turnOn(xx)/ratioEffTauLoose20MuTauMCEC->turnOn(xx);  
  
}  

Double_t myFuncTurnOnTauLoose20MuTauRunBEC(Double_t* x, Double_t *par) { //modified 

  gSystem->Load("ratioEfficiencyTest_C.so");   
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunBEC = new ratioEfficiencyTest(18.59856420, 2.49132550, 10.99643595, 1.50651123, 0.87952970);

  Float_t xx = x[0];   
   
  return ratioEffTauLoose20MuTauRunBEC->turnOn(xx);
}

/////////////////////////////////////////////////


/////////////////////////////////////////////////


Double_t myFuncRatioTauMuTauAllBL(Double_t* x, Double_t *par) { //modified

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunABL = new ratioEfficiencyTest(18.52262128, 1.85879597, 3.48843815, 1.15491294, 1.02489024); 
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunBBL = new ratioEfficiencyTest(17.92648563, 1.96846742, 4.46406075, 1.02023992, 1.52260575);

  ratioEfficiencyTest* ratioEffTauLoose20MuTauMCBL = new ratioEfficiencyTest(18.86257072, 0.25680380, 0.16916101, 2.42931257, 0.89590264);

  float weightLoose20RunA =    0.14; 
  float weightLoose20RunB =    0.86; 
   
  float total = weightLoose20RunA+weightLoose20RunB; 
  weightLoose20RunA /= total; 
  weightLoose20RunB /= total; 
 
  Float_t xx = x[0];

  return (ratioEffTauLoose20MuTauRunABL->turnOn(xx) * weightLoose20RunA + 
	  ratioEffTauLoose20MuTauRunBBL->turnOn(xx) * weightLoose20RunB
	  )/ratioEffTauLoose20MuTauMCBL->turnOn(xx) ;

}

Double_t myFuncTurnOnTauMuTauAllBL(Double_t* x, Double_t *par) { //modified

  gSystem->Load("ratioEfficiencyTest_C.so"); 
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunABL = new ratioEfficiencyTest(18.52262128, 1.85879597, 3.48843815, 1.15491294, 1.02489024);  
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunBBL = new ratioEfficiencyTest(17.92648563, 1.96846742, 4.46406075, 1.02023992, 1.52260575);
  
  float weightLoose20RunA =    0.14;
  float weightLoose20RunB =    0.86;
  
  float total = weightLoose20RunA+weightLoose20RunB;
  weightLoose20RunA /= total;
  weightLoose20RunB /= total;
  
  Float_t xx = x[0];
  
  return (ratioEffTauLoose20MuTauRunABL->turnOn(xx) * weightLoose20RunA + 
	  ratioEffTauLoose20MuTauRunBBL->turnOn(xx) * weightLoose20RunB
	  );
  
 

}



////////////////////////////////////////////////////////

Double_t myFuncRatioTauMuTauAllEC(Double_t* x, Double_t *par) { //modified

  gSystem->Load("ratioEfficiencyTest_C.so");  
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunAEC = new ratioEfficiencyTest(18.90119559, 0.14025596, 0.14482632, 1.56126508, 0.81188198); 
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunBEC = new ratioEfficiencyTest(18.59856420, 2.49132550, 10.99643595, 1.50651123, 0.87952970);
  
  ratioEfficiencyTest* ratioEffTauLoose20MuTauMCEC = new ratioEfficiencyTest(18.74764561, 1.82036845, 701.46994969, 101.57913480, 0.82547043);

  
  float weightLoose20RunA =    0.14;  
  float weightLoose20RunB =    0.86;  
    
  float total = weightLoose20RunA+weightLoose20RunB;  
  weightLoose20RunA /= total;  
  weightLoose20RunB /= total;  
  
  Float_t xx = x[0]; 
 
  return (ratioEffTauLoose20MuTauRunAEC->turnOn(xx) * weightLoose20RunA +  
          ratioEffTauLoose20MuTauRunBEC->turnOn(xx) * weightLoose20RunB 
          )/ratioEffTauLoose20MuTauMCEC->turnOn(xx) ;  
}

Double_t myFuncTurnOnTauMuTauAllEC(Double_t* x, Double_t *par) { //modified 
  
  gSystem->Load("ratioEfficiencyTest_C.so");  
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunAEC = new ratioEfficiencyTest(18.90119559, 0.14025596, 0.14482632, 1.56126508, 0.81188198);  
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunBEC = new ratioEfficiencyTest(18.59856420, 2.49132550, 10.99643595, 1.50651123, 0.87952970);
   
  float weightLoose20RunA =    0.14; 
  float weightLoose20RunB =    0.86; 
  
  float total = weightLoose20RunA+weightLoose20RunB; 
  weightLoose20RunA /= total; 
  weightLoose20RunB /= total; 
  
  Float_t xx = x[0]; 
  
  return (ratioEffTauLoose20MuTauRunAEC->turnOn(xx) * weightLoose20RunA +  
          ratioEffTauLoose20MuTauRunBEC->turnOn(xx) * weightLoose20RunB 
          ); 
 
} 

/////////////////////////////////////////////////////////

// HLT PFTauLoose20 / MC Summer12 / elec+tau

Double_t myFuncTurnOnTauLoose20ElecTauMC(Double_t* x, Double_t *par) { //modified
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MC = new ratioEfficiencyTest(18.77448606, 0.45765507, 0.26077509, 13.43372485, 0.88037836);
  
  Float_t xx = x[0];

  return ratioEffTauLoose20MC->turnOn(xx);

}

/////////////////////////////////////////////////

// HLT PFTauLoose20 / Run201BA    / e+tau

Double_t myFuncRatioTauLoose20ElecTauRunA(Double_t* x, Double_t *par) { //modified
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunA = new ratioEfficiencyTest(18.84658959, 0.25958704, 0.17300958, 2.43491208, 0.85872017);
  ratioEfficiencyTest* ratioEffTauLoose20MC = new ratioEfficiencyTest(18.77448606, 0.45765507, 0.26077509, 13.43372485, 0.88037836);

  
  Float_t xx = x[0];

  return ratioEffTauLoose20ElecTauRunA->turnOn(xx)/ ratioEffTauLoose20MC->turnOn(xx);

}

Double_t myFuncTurnOnTauLoose20ElecTauRunA(Double_t* x, Double_t *par) { //modified
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunA = new ratioEfficiencyTest(18.84658959, 0.25958704, 0.17300958, 2.43491208, 0.85872017);

  Float_t xx = x[0];

  return ratioEffTauLoose20ElecTauRunA->turnOn(xx);
  
}

/////////////////////////////////////////////////

// HLT PFTauLoose20 / Run2012B / e+tau

Double_t myFuncRatioTauLoose20ElecTauRunB(Double_t* x, Double_t *par) { //modified
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunB = new ratioEfficiencyTest(18.48663118, 1.63417147, 20.25695815, 138.55422224, 0.89456038);
  ratioEfficiencyTest* ratioEffTauLoose20MC = new ratioEfficiencyTest(18.77448606, 0.45765507, 0.26077509, 13.43372485, 0.88037836);
  
  Float_t xx = x[0];

  return ratioEffTauLoose20ElecTauRunB->turnOn(xx)/ ratioEffTauLoose20MC->turnOn(xx);

}

Double_t myFuncTurnOnTauLoose20ElecTauRunB(Double_t* x, Double_t *par) { //modified
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunB = new ratioEfficiencyTest(18.48663118, 1.63417147, 20.25695815, 138.55422224, 0.89456038);
  
  Float_t xx = x[0];
  
  return ratioEffTauLoose20ElecTauRunB->turnOn(xx);
  
}

/////////////////////////////////////////////////

// HLT PFTau All / e+tau

Double_t myFuncRatioTauElecTauAll(Double_t* x, Double_t *par) {  //modified
  
  gSystem->Load("ratioEfficiencyTest_C.so");
 
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunA = new ratioEfficiencyTest(18.84658959, 0.25958704, 0.17300958, 2.43491208, 0.85872017); 
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunB = new ratioEfficiencyTest(18.48663118, 1.63417147, 20.25695815, 138.55422224, 0.89456038);
  ratioEfficiencyTest* ratioEffTauLoose20MC = new ratioEfficiencyTest(18.77448606, 0.45765507, 0.26077509, 13.43372485, 0.88037836); 
  

  float weightLoose20RunA =    0.14;  
  float weightLoose20RunB =    0.86;  
   
  float total = weightLoose20RunA+weightLoose20RunB;  
  weightLoose20RunA /= total;  
  weightLoose20RunB /= total;
  
  Float_t xx = x[0];

  float ratio = ratioEffTauLoose20MC->turnOn(xx) > 0 ? 
    (ratioEffTauLoose20ElecTauRunA->turnOn(xx)  * weightLoose20RunA +
     ratioEffTauLoose20ElecTauRunB->turnOn(xx)  * weightLoose20RunB
     ) / ratioEffTauLoose20MC->turnOn(xx) : 0.;
  
  return ratio;

}


Double_t myFuncTurnOnTauElecTauAll(Double_t* x, Double_t *par) { //modified
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunA = new ratioEfficiencyTest(18.84658959, 0.25958704, 0.17300958, 2.43491208, 0.85872017);  
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunB = new ratioEfficiencyTest(18.48663118, 1.63417147, 20.25695815, 138.55422224, 0.89456038); 
  
  float weightLoose20RunA =    0.14;   
  float weightLoose20RunB =    0.86;   
    
  float total = weightLoose20RunA+weightLoose20RunB;   
  weightLoose20RunA /= total;   
  weightLoose20RunB /= total; 
   
  Float_t xx = x[0]; 
  
  return (ratioEffTauLoose20ElecTauRunA->turnOn(xx)  * weightLoose20RunA +
	  ratioEffTauLoose20ElecTauRunB->turnOn(xx)  * weightLoose20RunB
	  );

}


/////////////////////////////////////////////////

void makeFile(){

  TFile* fout = new TFile("llrCorrections.root","RECREATE");

  TF1 *ratioElec20BL        = new TF1("ratioElec20BL",           myFuncRatioElec20BL,       24,800,0);
  TF1 *turnOnElec20BL       = new TF1("turnOnElec20BL",          myFuncTurnOnElec20BL,      24,800,0);
  TF1 *ratioElec20EC        = new TF1("ratioElec20EC",           myFuncRatioElec20EC,       24,800,0);
  TF1 *turnOnElec20EC       = new TF1("turnOnElec20EC",          myFuncTurnOnElec20EC,      24,800,0);

  TF1 *ratioElec22BL        = new TF1("ratioElec22BL",           myFuncRatioEle22BL,        24,800,0);
  TF1 *turnOnElec22BL       = new TF1("turnOnElec22BL",          myFuncTurnOnEle22BL,       24,800,0);
  TF1 *ratioElec22EC        = new TF1("ratioElec22EC",           myFuncRatioEle22EC,        24,800,0);
  TF1 *turnOnElec22EC       = new TF1("turnOnElec22EC",          myFuncTurnOnEle22EC,       24,800,0);

  TF1 *ratioElecAllBL       = new TF1("ratioElecAllBL",          myFuncRatioEleAllBL,       24,800,0);
  TF1 *turnOnElecAllBL      = new TF1("turnOnElecAllBL",         myFuncTurnOnEleAllBL,      24,800,0);
  TF1 *ratioElecAllEC       = new TF1("ratioElecAllEC",          myFuncRatioEleAllEC,       24,800,0);
  TF1 *turnOnElecAllEC      = new TF1("turnOnElecAllEC",         myFuncTurnOnEleAllEC,      24,800,0);

  TF1 *ratioElecIDBL        = new TF1("ratioElecIDBL",           myFuncRatioElecIDBL ,      14,800,0);
  TF1 *turnOnElecIDBL       = new TF1("turnOnElecIDBL",          myFuncTurnOnElecIDBL ,     14,800,0);
  TF1 *ratioElecIDEC        = new TF1("ratioElecIDEC",           myFuncRatioElecIDEC ,      14,800,0);
  TF1 *turnOnElecIDEC       = new TF1("turnOnElecIDEC",          myFuncTurnOnElecIDEC ,     14,800,0);

  TF1 *ratioElecIsoBL       = new TF1("ratioElecIsoBL",          myFuncRatioElecIsoBL ,     14,800,0);
  TF1 *turnOnElecIsoBL      = new TF1("turnOnElecIsoBL",         myFuncTurnOnElecIsoBL ,    14,800,0);
  TF1 *ratioElecIsoEC       = new TF1("ratioElecIsoEC",          myFuncRatioElecIsoEC ,     14,800,0);
  TF1 *turnOnElecIsoEC      = new TF1("turnOnElecIsoEC",         myFuncTurnOnElecIsoEC ,    14,800,0);

  TF1 *ratioMu18BL      = new TF1("ratioMu18BL",         myFuncRatioMu18BL,     20,800,0);
  TF1 *turnOnMu18BL     = new TF1("turnOnMu18BL",        myFuncTurnOnMu18BL,    20,800,0);

  TF1 *ratioMu18EC      = new TF1("ratioMu18EC",         myFuncRatioMu18EC,     20,800,0);
  TF1 *turnOnMu18EC     = new TF1("turnOnMu18EC",        myFuncTurnOnMu18EC,    20,800,0);

  TF1 *ratioMu17BL      = new TF1("ratioMu17BL",         myFuncRatioMu17BL,     20,800,0);
  TF1 *turnOnMu17BL     = new TF1("turnOnMu17BL",        myFuncTurnOnMu17BL,    20,800,0);

  TF1 *ratioMu17EC      = new TF1("ratioMu17EC",         myFuncRatioMu17EC,     20,800,0);
  TF1 *turnOnMu17EC     = new TF1("turnOnMu17EC",        myFuncTurnOnMu17EC,    20,800,0);
 
  TF1 *ratioMuAllBL         = new TF1("ratioMuAllBL",            myFuncRatioMuAllBL,        20,800,0);
  TF1 *turnOnMuAllBL        = new TF1("turnOnMuAllBL",           myFuncTurnOnMuAllBL,       20,800,0);

  TF1 *ratioMuAllEC         = new TF1("ratioMuAllEC",            myFuncRatioMuAllEC,        20,800,0);
  TF1 *turnOnMuAllEC        = new TF1("turnOnMuAllEC",           myFuncTurnOnMuAllEC,       20,800,0);
 
  TF1 *ratioMuIDBL          = new TF1("ratioMuIDBL",             myFuncRatioMuIDBL ,        14,800,0);
  TF1 *turnOnMuIDBL         = new TF1("turnOnMuIDBL",            myFuncTurnOnMuIDBL ,       14,800,0);

  TF1 *ratioMuIDEC          = new TF1("ratioMuIDEC",             myFuncRatioMuIDEC ,        14,800,0);
  TF1 *turnOnMuIDEC         = new TF1("turnOnMuIDEC",            myFuncTurnOnMuIDEC ,       14,800,0);

  TF1 *ratioMuIsoBL         = new TF1("ratioMuIsoBL",            myFuncRatioMuIsoBL ,       14,800,0);
  TF1 *turnOnMuIsoBL        = new TF1("turnOnMuIsoBL",           myFuncTurnOnMuIsoBL ,      14,800,0);

  TF1 *ratioMuIsoEC         = new TF1("ratioMuIsoEC",            myFuncRatioMuIsoEC ,       14,800,0);
  TF1 *turnOnMuIsoEC        = new TF1("turnOnMuIsoEC",           myFuncTurnOnMuIsoEC ,      14,800,0);

  TF1 *turnOnTauLoose20MuTauMCBL      = new TF1("turnOnTauLoose20MuTauMCBL",   myFuncTurnOnTauLoose20MuTauMCBL,      20,800,0);
  TF1 *turnOnTauLoose20MuTauMCEC      = new TF1("turnOnTauLoose20MuTauMCEC",   myFuncTurnOnTauLoose20MuTauMCEC,	     20,800,0); 

  TF1 *ratioTauLoose20MuTauRunABL    = new TF1("ratioTauLoose20MuTauRunABL",   myFuncRatioTauLoose20MuTauRunABL,     20,800,0);
  TF1 *turnOnTauLoose20MuTauRunABL   = new TF1("turnOnTauLoose20MuTauRunABL",  myFuncTurnOnTauLoose20MuTauRunABL,    20,800,0);
  TF1 *ratioTauLoose20MuTauRunAEC    = new TF1("ratioTauLoose20MuTauRunAEC",   myFuncRatioTauLoose20MuTauRunAEC,     20,800,0); 
  TF1 *turnOnTauLoose20MuTauRunAEC   = new TF1("turnOnTauLoose20MuTauRunAEC",  myFuncTurnOnTauLoose20MuTauRunAEC,    20,800,0);
  
  TF1 *turnOnTauLoose20ElecTauMC      = new TF1("turnOnTauLoose20ElecTauMC",   myFuncTurnOnTauLoose20ElecTauMC,      20,800,0);
  
  TF1 *ratioTauLoose20ElecTauRunA    = new TF1("ratioTauLoose20ElecTauRunA",   myFuncRatioTauLoose20ElecTauRunA ,  20,800,0);
  TF1 *turnOnTauLoose20ElecTauRunA   = new TF1("turnOnTauLoose20ElecTauRunA",  myFuncTurnOnTauLoose20ElecTauRunA,  20,800,0);
  
  TF1 *ratioTauLoose20MuTauRunBBL    = new TF1("ratioTauLoose20MuTauRunBBL",   myFuncRatioTauLoose20MuTauRunBBL,     20,800,0);
  TF1 *turnOnTauLoose20MuTauRunBBL   = new TF1("turnOnTauLoose20MuTauRunBBL",  myFuncTurnOnTauLoose20MuTauRunBBL,    20,800,0);
  TF1 *ratioTauLoose20MuTauRunBEC    = new TF1("ratioTauLoose20MuTauRunBEC",   myFuncRatioTauLoose20MuTauRunBEC,     20,800,0); 
  TF1 *turnOnTauLoose20MuTauRunBEC   = new TF1("turnOnTauLoose20MuTauRunBEC",  myFuncTurnOnTauLoose20MuTauRunBEC,    20,800,0);

  TF1 *ratioTauLoose20ElecTauRunB    = new TF1("ratioTauLoose20ElecTauRunB",   myFuncRatioTauLoose20ElecTauRunB ,  20,800,0);
  TF1 *turnOnTauLoose20ElecTauRunB   = new TF1("turnOnTauLoose20ElecTauRunB",  myFuncTurnOnTauLoose20ElecTauRunB,  20,800,0);

 
  TF1 *ratioTauElecTauAll    = new TF1("ratioTauElecTauAll", myFuncRatioTauElecTauAll ,20,800,0);
  TF1 *turnOnTauElecTauAll   = new TF1("turnOnTauElecTauAll",myFuncTurnOnTauElecTauAll,20,800,0);
 
  TF1 *ratioTauMuTauAllBL      = new TF1("ratioTauMuTauAllBL", myFuncRatioTauMuTauAllBL ,20,800,0);
  TF1 *turnOnTauMuTauAllBL     = new TF1("turnOnTauMuTauAllBL",myFuncTurnOnTauMuTauAllBL,20,800,0);
  TF1 *ratioTauMuTauAllEC      = new TF1("ratioTauMuTauAllEC", myFuncRatioTauMuTauAllEC ,20,800,0); 
  TF1 *turnOnTauMuTauAllEC     = new TF1("turnOnTauMuTauAllEC",myFuncTurnOnTauMuTauAllEC,20,800,0);
  

  fout->cd();
  
  

  ratioElecIsoBL->SetNpx(25600);
  turnOnElecIsoBL->SetNpx(25600);
  ratioElecIsoEC->SetNpx(25600);
  turnOnElecIsoEC->SetNpx(25600);
  ratioMuIsoBL->SetNpx(25600);
  ratioMuIsoEC->SetNpx(25600);
  turnOnMuIsoBL->SetNpx(25600);
  turnOnMuIsoEC->SetNpx(25600);
  
  ratioElec20BL->SetNpx(25600);
  turnOnElec20BL->SetNpx(25600);
  ratioElec20EC->SetNpx(25600);
  turnOnElec20EC->SetNpx(25600);
  
  ratioElec22BL->SetNpx(25600); 
  turnOnElec22BL->SetNpx(25600); 
  ratioElec22EC->SetNpx(25600); 
  turnOnElec22EC->SetNpx(25600);

  ratioElecAllBL->SetNpx(25600);
  turnOnElecAllBL->SetNpx(25600);
  ratioElecAllEC->SetNpx(25600);
  turnOnElecAllEC->SetNpx(25600);

  
  ratioMu17BL->SetNpx(25600);
  turnOnMu17BL->SetNpx(25600);
  ratioMu17EC->SetNpx(25600);
  turnOnMu17EC->SetNpx(25600);

  ratioMu18BL->SetNpx(25600); 
  turnOnMu18BL->SetNpx(25600); 
  ratioMu18EC->SetNpx(25600); 
  turnOnMu18EC->SetNpx(25600);

  ratioMuAllBL->SetNpx(25600);
  turnOnMuAllBL->SetNpx(25600);
  ratioMuAllEC->SetNpx(25600);
  turnOnMuAllEC->SetNpx(25600);


  ratioTauElecTauAll->SetNpx(25600);
  turnOnTauElecTauAll->SetNpx(25600);

  ratioTauMuTauAllBL->SetNpx(25600);
  turnOnTauMuTauAllBL->SetNpx(25600);
  ratioTauMuTauAllEC->SetNpx(25600); 
  turnOnTauMuTauAllEC->SetNpx(25600);

  turnOnTauLoose20MuTauMCBL->SetNpx(25600);
  turnOnTauLoose20MuTauMCEC->SetNpx(25600);
  ratioTauLoose20MuTauRunABL->SetNpx(25600);
  turnOnTauLoose20MuTauRunABL->SetNpx(25600);
  ratioTauLoose20MuTauRunAEC->SetNpx(25600);
  turnOnTauLoose20MuTauRunAEC->SetNpx(25600);
  turnOnTauLoose20ElecTauMC->SetNpx(25600);
  ratioTauLoose20ElecTauRunA->SetNpx(25600);
  turnOnTauLoose20ElecTauRunA->SetNpx(25600);
  ratioTauLoose20MuTauRunBBL->SetNpx(25600);
  turnOnTauLoose20MuTauRunBBL->SetNpx(25600);
  ratioTauLoose20MuTauRunBEC->SetNpx(25600);
  turnOnTauLoose20MuTauRunBEC->SetNpx(25600);
  ratioTauLoose20ElecTauRunB->SetNpx(25600);
  turnOnTauLoose20ElecTauRunB->SetNpx(25600);


  ratioElecIDBL->SetNpx(25600);
  turnOnElecIDBL->SetNpx(25600);
  ratioElecIDEC->SetNpx(25600);
  turnOnElecIDEC->SetNpx(25600);
  ratioMuIDBL->SetNpx(25600);
  turnOnMuIDBL->SetNpx(25600);
  ratioMuIDEC->SetNpx(25600);
  turnOnMuIDEC->SetNpx(25600);

  //==============================

  ratioElecIsoBL->Write();
  turnOnElecIsoBL->Write();
  ratioElecIsoEC->Write();
  turnOnElecIsoEC->Write();
  ratioMuIsoBL->Write();
  ratioMuIsoEC->Write();
  turnOnMuIsoBL->Write();
  turnOnMuIsoEC->Write();
  
  ratioElec20BL->Write();
  turnOnElec20BL->Write();
  ratioElec20EC->Write();
  turnOnElec20EC->Write();
  
  ratioElec22BL->Write(); 
  turnOnElec22BL->Write(); 
  ratioElec22EC->Write(); 
  turnOnElec22EC->Write();

  ratioElecAllBL->Write();
  turnOnElecAllBL->Write();
  ratioElecAllEC->Write();
  turnOnElecAllEC->Write();
  
  ratioMu17BL->Write();
  turnOnMu17BL->Write();
  ratioMu17EC->Write();
  turnOnMu17EC->Write();

  ratioMu18BL->Write(); 
  turnOnMu18BL->Write(); 
  ratioMu18EC->Write(); 
  turnOnMu18EC->Write();

  ratioMuAllBL->Write();
  turnOnMuAllBL->Write();
  ratioMuAllEC->Write();
  turnOnMuAllEC->Write();


  ratioTauElecTauAll->Write();
  turnOnTauElecTauAll->Write();

  ratioTauMuTauAllBL->Write();
  turnOnTauMuTauAllBL->Write();
  ratioTauMuTauAllEC->Write(); 
  turnOnTauMuTauAllEC->Write();

  turnOnTauLoose20MuTauMCBL->Write();
  turnOnTauLoose20MuTauMCEC->Write();
  ratioTauLoose20MuTauRunABL->Write();
  turnOnTauLoose20MuTauRunABL->Write();
  ratioTauLoose20MuTauRunAEC->Write();
  turnOnTauLoose20MuTauRunAEC->Write();
  turnOnTauLoose20ElecTauMC->Write();
  ratioTauLoose20ElecTauRunA->Write();
  turnOnTauLoose20ElecTauRunA->Write();
  ratioTauLoose20MuTauRunBBL->Write();
  turnOnTauLoose20MuTauRunBBL->Write();
  ratioTauLoose20MuTauRunBEC->Write();
  turnOnTauLoose20MuTauRunBEC->Write();
  ratioTauLoose20ElecTauRunB->Write();
  turnOnTauLoose20ElecTauRunB->Write();


  ratioElecIDBL->Write();
  turnOnElecIDBL->Write();
  ratioElecIDEC->Write();
  turnOnElecIDEC->Write();
  ratioMuIDBL->Write();
  turnOnMuIDBL->Write();
  ratioMuIDEC->Write();
  turnOnMuIDEC->Write();


  fout->Write();
  fout->Close();


}

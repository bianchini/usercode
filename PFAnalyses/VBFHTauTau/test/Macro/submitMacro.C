

{

  //gSystem->Load("templateStudy_C.so");

  //bkgTemplateStudy(1000,440,0.32,1,"etoTauSCMargNoCracks80",true,"tauAntiEMVA", ">=",0.5,"abseta<1.5","abseta>-1",18,65,120);
  //bkgTemplateStudy(500,440,0.32,2,"etoTauSCMargNoCracks80",true,"tauAntiEMVA", ">=",0.5,"abseta<1.5","abseta>-1",18,65,120);
  //bkgTemplateStudy(500,440,0.32,3,"etoTauSCMargNoCracks80",true,"tauAntiEMVA", ">=",0.5,"abseta<1.5","abseta>-1",18,65,120);
  //bkgTemplateStudy(1000,220,0.40,1,"etoTauSCMargNoCracks80",true,"tauAntiEMVA", ">=",0.5,"abseta>1.5","abseta>-1",18,65,120);
  //bkgTemplateStudy(500,220,0.40,2,"etoTauSCMargNoCracks80",true,"tauAntiEMVA", ">=",0.5,"abseta>1.5","abseta>-1",18,65,120);
  //bkgTemplateStudy(500,220,0.40,3,"etoTauSCMargNoCracks80",true,"tauAntiEMVA", ">=",0.5,"abseta>1.5","abseta>-1",18,65,120);

  //bkgTemplateStudy(1000,220,0.75,1,"etoTauMargLooseNoCracks70",true,"tauAntiEMVA", ">=",0.5,"abseta<1.5","abseta>-1",18,65,120);
  //bkgTemplateStudy(500,220,0.75,2,"etoTauMargLooseNoCracks70",true,"tauAntiEMVA", ">=",0.5,"abseta<1.5","abseta>-1",18,65,120);
  //bkgTemplateStudy(500,220,0.75,3,"etoTauMargLooseNoCracks70",true,"tauAntiEMVA", ">=",0.5,"abseta<1.5","abseta>-1",18,65,120);
  //bkgTemplateStudy(1000,120,0.70,1,"etoTauMargLooseNoCracks70",true,"tauAntiEMVA", ">=",0.5,"abseta>1.5","abseta>-1",18,65,120);
  //bkgTemplateStudy(500,120,0.70,2,"etoTauMargLooseNoCracks70",true,"tauAntiEMVA", ">=",0.5,"abseta>1.5","abseta>-1",18,65,120);
  //bkgTemplateStudy(500,120,0.70,3,"etoTauMargLooseNoCracks70",true,"tauAntiEMVA", ">=",0.5,"abseta>1.5","abseta>-1",18,65,120);


  gSystem->Load("fitHisto_C.so");

  makePlot("etoTauSCMargNoCracks80","tauAntiEMVA",0.5,"abseta","abseta>-1",65,120,18,true,0.0,1.5,2.5,0.06);
  makePlot("etoTauSCMargNoCracks70","matchedID",0.975,"abseta","abseta>-1",65,120,18,false,0.0,1.5,2.5,0.15);
  makePlot("etoTauSCMargNoCracks60","-electronPreIDOutput",-0.6,"abseta","abseta>-1",65,120,18,false,0.0,1.5,2.5,0.30);

  makePlot("etoTauMargLooseNoCracks70","matchedID",0.975,"abseta","abseta>-1",65,120,18,false,0.0,1.5,2.5,0.15);
  makePlot("etoTauMargLooseNoCracks70","-electronPreIDOutput",-0.6,"abseta","abseta>-1",65,120,18,false,0.0,1.5,2.5,0.30);
  makePlot("etoTauMargLooseNoCracks70","tauAntiEMVA",0.5,"abseta","abseta>-1",65,120,18,true,0.0,1.5,2.5,0.06);
  

  //makePlot("etoTauMargLooseNoCracks70","tauAntiEMVA",0.5,"abseta","leadPFChargedHadrCandTrackPt>5 && pt>15",65,120,18,true,0.0,1.5,2.5);
 
  //makePlot("etoTauSCMargNoCracks70",   "matchedID",0.975,"abseta", 65,120,25,0.0,1.5,2.5);
  //makePlot("etoTauSCMargNoCracks70",   "tauAntiEMVA",0.5,"abseta","pt>15",65,120,18,false,0.0,1.5,2.5);

  //makePlot("etoTauMargMediumNoCracks70","tauAntiEMVA",0.5,"abseta", 65,120,18,0.0,1.5,2.5);
  //makePlot("etoTauMargTightNoCracks70", "tauAntiEMVA",0.5,"abseta", 65,120,25,0.0,1.5,2.5);

 

}

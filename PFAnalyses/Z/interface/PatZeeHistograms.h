#ifndef PatZeeHistograms_h
#define PatZeeHistograms_h

// Base class for histogram managing.
//
// Original Author:  Artur Kalinowski
//         Created:  Wed Jul 22 12:56:54 CEST 2009
// $Id: PatZeeHistograms.h,v 1.1 2010/04/16 08:39:55 bianchi Exp $
//
//
#include "PFAnalyses/CommonTools/interface/AnalysisHistograms.h"

class PatZeeHistograms: public AnalysisHistograms {
 public:

  PatZeeHistograms(std::string fileName="TestHistos.root", int opt=0);

  PatZeeHistograms( TFileDirectory *myDir);

  PatZeeHistograms( TFileDirectory *myDir, const std::string & fileName = "");

  virtual ~PatZeeHistograms();

  void fillHistograms(float et);

 private:

  virtual void defineHistograms();


};

#endif

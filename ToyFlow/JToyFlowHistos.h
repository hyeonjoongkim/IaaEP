//===========================================================
// JToyFlowHistos.h
//===========================================================

#ifndef JTOYFLOWHISTOS_H
#define JTOYFLOWHISTOS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TFile.h>
#include "JConst.h"


class JToyFlowHistos {

public:
  JToyFlowHistos(); 
  virtual ~JToyFlowHistos(){;}	  //destructor
  JToyFlowHistos(const JToyFlowHistos& obj);
  JToyFlowHistos& operator=(const JToyFlowHistos& obj); 
  // ALICE methods =====================================================
  void CreateHistos();

private:

public:
  //===================================================
  // Toy Flow histograms
  //===================================================
  TH1D *pah[R_COUNT][NC], *pbh[R_COUNT][NC], *pch[R_COUNT][NC];
  TH1D *evph[R_COUNT][NC];
  TH2D *contami2d[R_COUNT][NC];
  TH2D *highcontami2d[R_COUNT][NC];
  TH2D *evpcorr2d[R_COUNT][NC];
  TH2D *evpcorrvsdet2d[NC];
  TH1D *evpdifference[R_COUNT][NC];
  TH1D *hPhiEbE[NC];

};

#endif



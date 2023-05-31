#ifndef NURATEGENERATOR_H
#define NURATEGENERATOR_H

#include "TH2.h"
#include "TGraph.h"

//TH2D*  GetTrueEvents(int flvf, double Energy[], double CosT[] ,int nbinsx, int nbinsy);
TH2D*  ObservedEvents(int flvf, double Energy[], double CosT[] ,int Ebins, int Tbins);


//Generate energy dependent Oscillation Probabilities  (fixed Baseline) 
//TGraph*  Matter_Osc(int flvi, int flvf, double Emin , double Emax,double rho , double L);

#endif
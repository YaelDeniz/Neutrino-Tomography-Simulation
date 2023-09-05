#ifndef OBSEVENTGENERATOR_H
#define OBSEVENTGENERATOR_H

#include "TH2.h"
#include <string.h>
#include <stdlib.h>


TH2D*  GetObservedEvents( std::string model, int flvf, double E_GeV[], double Eta[],double dAz ,int Ebins, int Tbins, double Det_par[] );


//Generate energy dependent Oscillation Probabilities  (fixed Baseline) 
//TGraph*  Matter_Osc(int flvi, int flvf, double Emin , double Emax,double rho , double L);

#endif
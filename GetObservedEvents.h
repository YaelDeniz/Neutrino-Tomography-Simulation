#ifndef OBSEVENTGENERATOR_H
#define OBSEVENTGENERATOR_H

#include "TH2.h"
#include "TObjArray.h"
#include <string.h>
#include <stdlib.h>

TH2D*  AsimovObservedEvents(std::string modelname, int flvf, double Region[], int Bins[],double Det_par[], double NnT);

//TH2D*  GetObservedEvents(std::string modelname, int flvf, double Region[], int Bins[],double Det_par[], double NnT);

//TObjArray*  SmearEvents( TObjArray* TrueEvents, int flvf, double Region[], int Bins[],double Det_par[], double NnT, int K );


//Generate energy dependent Oscillation Probabilities  (fixed Baseline) 
//TGraph*  Matter_Osc(int flvi, int flvf, double Emin , double Emax,double rho , double L);

#endif
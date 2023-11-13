#ifndef TRUEEVENTGENERATOR_H
#define TRUEEVENTGENERATOR_H

#include "TH2.h"
#include "TObjArray.h"
#include <string.h>
#include <stdlib.h>

TH2D*  AsimovTrueEvents(std::string modelname, int flvf, double Region[], int Bins[], double NnT);

//TObjArray*  GetTrueEventsMC(std::string modelname, int flvf , double Region[], int Bins[], double NnT, int K);
//TH2D*  GetTrueEvents(std::string modelname, int flvf , double Region[], int Bins[], double NnT);
/*
List of Variables 

Region[]={Emin, Emax, Etamin, Etamax, dAz}

Bins[]={Angular bins, Energy bins}

NnT = Gigaton-year exposure

K = Number of Pseudo Experiments

*/


#endif
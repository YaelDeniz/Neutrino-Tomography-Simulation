#ifndef TRUEEVENTGENERATOR_H
#define TRUEEVENTGENERATOR_H

#include "TH2.h"
#include <string.h>
#include <stdlib.h>


TH2D*  GetTrueEvents(std::string model, int flvf, double E_GeV[], double Eta[], double dAz ,int Ebins, int Tbins, double NnT);
//TGraph*  Matter_Osc(int flvi, int flvf, double Emin , double Emax,double rho , double L);

#endif
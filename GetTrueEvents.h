#ifndef TRUEEVENTGENERATOR_H
#define TRUEEVENTGENERATOR_H

#include "TH2.h"


TH2D*  GetTrueEvents(int flvf, double E_GeV[], double Eta[], double dAz ,int Ebins, int Tbins);
//TGraph*  Matter_Osc(int flvi, int flvf, double Emin , double Emax,double rho , double L);

#endif
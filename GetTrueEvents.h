#ifndef TRUEEVENTGENERATOR_H
#define TRUEEVENTGENERATOR_H

#include "TH2.h"


TH2D*  GetTrueEvents(int flvf, double Energy[], double CosT[] ,int nbinsx, int nbinsy);
//TGraph*  Matter_Osc(int flvi, int flvf, double Emin , double Emax,double rho , double L);

#endif
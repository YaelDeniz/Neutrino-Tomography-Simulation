#ifndef TRUEEVENTGENERATOR_H
#define TRUEEVENTGENERATOR_H

//C++
#include <string.h>
#include <stdlib.h>
#include <vector>

//CERN ROOT
#include "TH2.h"
#include "TH3.h"
#include "TObjArray.h"


using namespace std;

//TH2D*  AsimovTrueEvents(std::string modelname,bool LLVP,  std::vector<int> layers, int flvf, double Region[], int Bins[], double NnT);

class AsimovSimulation
{
    public:

    //Earth Settings
    std::string PremModel="prem_44layers";
    bool MantleAnomaly = false;
    std::string AnomalyShape ="pancake";
    //std::vector<int> AnomalousLayers;

    //Modify specific Layers

    int PremTableNumber = 44 ;
    double DensityContrast = 0 ;
    double ChemicalContrats = 0 ;

    //Neutrino Settings
    int flvf;
    double EnuMin = 1.00; //GeV
    double EnuMax = 10.00; //GeV

    double ZenMin; // Min 90;
    double ZenMax; // Max 180;

    double AziMin=0; // Min 0;
    double AziMax=360; //Max 360;

    //Simulation Settings
    int nbinsZen; //Bins in Zenith
    int nbinsAzi; //Bins in Azimuth
    int nbinsE; //Bins in Energy

    double NnT; //Detector Exposure

    void ModifyLayer (int n , double rhopct , double zoapct)
    {

        PremTableNumber = n; // Index of layer
        DensityContrast = rhopct; // Density percentage difference
        ChemicalContrats= zoapct; // zoa(compositional) perctengae differnce

    }

    void SetIntervals(double zlow, double zup, double alow, double aup, double elow, double eup)
    {
        ZenMin = zlow;
        ZenMax = zup;

        AziMin = alow;
        AziMax = aup;
        
        EnuMin = elow;
        EnuMax = eup;

    }

    void SetBinning(int zbins, int abins, int ebins)
    {
        nbinsZen = zbins;
        nbinsAzi = abins;
        nbinsE   = ebins;
    }

    void SetExposure(double exposure) 
    {

        NnT = exposure;

    }

    TH3D * GetTrueEvents3D();

    TH2D * GetTrueEvents2D(); //Azimuth independet

    TH2D * TestTrueEvents2D(std::string model_std , std::string model_alt);




};




TH2D*  OscProbEarth(std::string modelname, int flvf, double Region[], int Bins[], double NnT);

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
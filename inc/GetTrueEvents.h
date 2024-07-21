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
#include "TGraph.h"



using namespace std;
# define Rocean 6368.0 // Earth radius in (km)
# define Rearth 6371.0 // Earth radius in (km)

//TH2D*  AsimovTrueEvents(std::string modelname,bool LLVP,  std::vector<int> layers, int flvf, double Region[], int Bins[], double NnT);

class AsimovSimulation
{
    public:
    

    double NnT; //Detector Exposure
  
    
    double Rdet[3] = {0.0, 0.0, -1.0*Rearth};
    std::string FluxFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/";
    std::string FluxFile = "SP_AziAveraged_solmin/spl-nu-20-01-000.d";
    std::string HondaTable = FluxFolder + FluxFile; //Class Assumes Sout Pole flux as default

    //Earth Settings

    //The Earth model
    std::string PremFolder = "/home/dehy0499/OscProb/PremTables/";
    std::string PremName = "prem_44layers.txt";
    std::string PremTable = PremFolder+PremName;
 
    //The LLVP
    bool MantleAnomaly = false;
    double PileHeight = 1000;
    double aperture = 45.0;
    std::string AnomalyShape ="cake";
    double PileDensityContrast = 2.0;
    double PileChemContrast = 0.0;
    

    //Modify specific Layers
    int PremTableNumber = 44 ;
    double DensityContrast = 0 ;
    double ChemicalContrast = 0 ;

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

    void SetDetectorPosition( double x[3])
    {

        Rdet[0] = x[0];
        Rdet[1] = x[1];
        Rdet[2] = x[2];

    }


    void ModifyLayer (int n , double rhopct , double zoapct)
    {

        PremTableNumber = n; // Index of layer
        DensityContrast = rhopct; // Density percentage difference
        ChemicalContrast = zoapct; // zoa(compositional) perctengae differnce

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

    TH2D * GetOscProb2D(int flvi, int flvf, bool nunubar  ); //Azimuth independet
    
    TGraph * GetOscProb( int flvi, int flvf, bool nunubar, double cth ); //To be Deleted
    
    TH2D * SensitivityTrueEvents2D( int n, double pct ); //To be Del

    //TH2D * TestTrueEvents2D(std::string model_std , std::string model_alt);

    std::vector< std::vector<double> > GetPremMatrix( std::string path2table  );




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
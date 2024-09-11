#ifndef OBSEVENTGENERATOR_H
#define OBSEVENTGENERATOR_H

//C++
#include <string.h>
#include <stdlib.h>
#include <vector>

//CERN ROOT
#include "TH2.h"
#include "TH3.h"
#include "TObjArray.h"
#include "TGraph.h"

# define Rocean 6368.0 // Earth radius in (km)
# define Rearth 6371.0 // Earth radius in (km)

class AsimovObsSimulation
{
    public:

    //Detector
    double xyzTelescope[3] = {0.0, 0.0, -1.0*Rearth};

    std::string FluxFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/";
    std::string FluxFile = "SP_AziAveraged_solmin/spl-nu-20-01-000.d";
    std::string TheHondaTable = FluxFolder + FluxFile; //Class Assumes Sout Pole flux as default

    //The Earth model
    std::string PremFolder = "/home/dehy0499/OscProb/PremTables/";
    std::string PremName = "prem_44layers.txt";
    std::string ThePremTable = PremFolder+PremName;

    //Modify specific Layers
    int PremLayer = 44 ;
    double DensityContrast = 0;
    double ChemicalContrats = 0;

    //LLVP

    bool PileInModel = false; 
    std::string ShapeOfPile = "cake";
    double ThePileHight = 1000.0;
    double ThePileAperture=45.0;
    double ThePileDensityContrast = 10.0;
    double ThePileChemicalContrast = 10.0;

    //Neutrino Settings
    int flvf;
    double EnuMin = 1.00; //GeV
    double EnuMax = 10.00; //GeV

    double ZenMin; // Min 90;
    double ZenMax; // Max 180;

    double AziMin=0; // Min 0;
    double AziMax=360; //Max 360;

    //Simulation Settings
    int truebinsZen; //Bins in Zenith
    int truebinsAzi; //Bins in Azimuth
    int truebinsE; //Bins in Energy

    int recobinsZen; //Bins in Zenith
    int recobinsAzi; //Bins in Azimuth
    int recobinsE; //Bins in Energy


    double NnT; //Detector Exposure //0.05*e +0.1*sqrt(e)
    
    //Detector Resolution-------------------------------------------------------

    double Ae = 0.05;
    double Be = 0.1;

    double Ath = 0.04;  //0.04 + 0.18/(sqrt(e))
    double Bth = 0.18;

    void SetEnergyResolution(double A, double B);

	double PDFE(double Ereco , double Etrue );//Gaussian

    void SetAngularResolution(double A, double B);

    double PDFth(double thReco, double Etrue, double thTrue); //Gaussian

    double PDFcth(double cthreco, double Etrue, double cthtrue); //Von Mises Fisher






    void ModifyLayer (int LayerNumber , double rhopct , double zoapct)
    {

        PremLayer = LayerNumber; // Index of layer
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

    void SetTrueBinning(int zbins, int abins, int ebins)
    {
        truebinsZen = zbins;
        truebinsAzi = abins;
        truebinsE   = ebins;
    }

    void SetRecoBinning(int zbins, int abins, int ebins) 
    {
        recobinsZen = zbins;
        recobinsAzi = abins+1;
        recobinsE   = ebins;
    }


    void SetExposure(double exposure) 
    {

        NnT = exposure;

    }

    void SetDetectorXYZ(double xyz[3])
    {
        xyzTelescope[0]=xyz[0];
        xyzTelescope[1]=xyz[1];
        xyzTelescope[2]=xyz[2];

    }

    //TH3D * GetObsEvents3D();

    TH2D * GetObsEvents2Dcth(); //Azimuth independet

    std::vector < TH2D* >  GetObsEvents3Dcth();


    //TH2D * TestObsEvents2D(std::string model_std , std::string model_alt);




};

//TH2D*  AsimovObservedEvents(std::string modelname, int flvf, double Region[], int Bins[],double Det_par[], double NnT);

//TH2D*  DetectorRes( double Region[], int Bins[],double Det_par[] );

//TH2D*  GetObservedEvents(std::string modelname, int flvf, double Region[], int Bins[],double Det_par[], double NnT);

//TObjArray*  SmearEvents( TObjArray* TrueEvents, int flvf, double Region[], int Bins[],double Det_par[], double NnT, int K );


//Generate energy dependent Oscillation Probabilities  (fixed Baseline) 
//TGraph*  Matter_Osc(int flvi, int flvf, double Emin , double Emax,double rho , double L);

#endif
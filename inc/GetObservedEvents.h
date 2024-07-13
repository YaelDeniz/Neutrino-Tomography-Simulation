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

class AsimovObsSimulation
{
    public:

    //Earth Settings
    std::string PremModel="prem_44layers";
    bool Pile = false;
    std::string PileShape ="pancake";

    //Modify specific Layers
    int PremLayer = 44 ;
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






    void ModifyLayer (int LayerNumber , double PctInDensity , double PctInZoA)
    {

        PremTableNumber = LayerNumber; // Index of layer
        DensityContrast = PctInDensity; // Density percentage difference
        ChemicalContrats= PctInZoA; // zoa(compositional) perctengae differnce

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
        recobinsAzi = abins;
        recobinsE   = ebins;
    }


    void SetExposure(double exposure) 
    {

        NnT = exposure;

    }

    TH3D * GetObsEvents3D();

    TH2D * GetObsEvents2D(); //Azimuth independet

    TH2D * TestObsEvents2D(std::string model_std , std::string model_alt);




};

TH2D*  AsimovObservedEvents(std::string modelname, int flvf, double Region[], int Bins[],double Det_par[], double NnT);

TH2D*  DetectorRes( double Region[], int Bins[],double Det_par[] );

//TH2D*  GetObservedEvents(std::string modelname, int flvf, double Region[], int Bins[],double Det_par[], double NnT);

//TObjArray*  SmearEvents( TObjArray* TrueEvents, int flvf, double Region[], int Bins[],double Det_par[], double NnT, int K );


//Generate energy dependent Oscillation Probabilities  (fixed Baseline) 
//TGraph*  Matter_Osc(int flvi, int flvf, double Emin , double Emax,double rho , double L);

#endif
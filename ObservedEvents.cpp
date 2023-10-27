#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

//Cern ROOT
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "Math/Interpolator.h"
//My Libraries
#include "MathTools.h"
#include "PhyTools.h"
#include "GetObservedEvents.h"


# define myPi 3.14159265358979323846  /* pi */
// Some Constants i need
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds
# define R_earth 6368.0 //km



using namespace std;


int main()
{
    // Plot and Histogram2D settings

    std::cout << " Neutrino Oscillation tomography. " << std::endl;
    std::cout << " Observed Events: track-like only. "<< std::endl; 

    int Ebins=100; // # of Bins of Energy
    int Etabins=100; // # of Bins of cosEta

    int Ebins_o=20; // # of Bins of Energy
    int Etabins_o=20; // # of Bins of cosEta

    int Bins[]={Etabins, Ebins, Etabins_o, Ebins_o};


    //Energy interval (in GeV):
    double Emin = 1.0; 
    double Emax = 10.0;

    //------------------

    // Alternative Earth Model:

    double h_llsvp = 700; //km
    double R_cmb = 3500.0;
    double R_llsvp = R_cmb + h_llsvp; //km

    double R_min = 3000.0; // Distance from the center of the Earth , 3500 Km CMB
    //double R_llsvp = R_cmb + h_llsvp; //Km
    double R_max = 4500; // Distance from the center of the Earth

    double Etamin = TMath::ASin( (R_min)/R_earth )*(180.0/TMath::Pi()) ;
    //double Etamax_LLSVP = TMath::ASin( (R_cmb + h_llsvp)/R_earth )*(180.0/TMath::Pi()) ;
    double Etamax = TMath::ASin( R_max/R_earth )*(180.0/TMath::Pi()) ;


    //------------------

    double Phim = 0.0;
    double PhiM = 80.0 ; 

    double dAz = PhiM-Phim;

    double Region[] = {Emin,Emax,Etamin,Etamax,dAz};

    
    //Detector Properties:

    //Resolution

    double a_E = 0.2;

    double a_Eta= 0.25;

    double Det_par[] = {a_E,a_Eta};
    
    //Exposure
    
    double T = 10; // Exposure time in years 

    double DetMass = 10; // Detetor mass in Mtons

    double NnT      = (DetMass*MTon)*(T*years2sec)/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)

    std::cout <<"Simulation for T= "<<T<<" years and M= "<<DetMass <<" Mton"<<std::endl; 

    std::cout <<"a_E= "<< a_E<<" a_Eta= "<< a_Eta <<std::endl; 
    


    int flvf = 1;

    std::string  modelname;
            
    //prem_default   = "/home/dehy0499/OscProb/PremTables/prem_default.txt"; //Specify PREM table from OscProb

    modelname = "prem_default";

    //TH2D* EventOsc = GetObservedEvents(prem_default, flvf, E, Eta , dAz , Ebins, Tbins, Det_par);




    //-------------------------------------------------------------------------------------------------------------------

    //Density contrast
    double drho_dp = 2.0; // 2 % density contrats

    std::string prem_llsvp, prem_default;
    
    std::ifstream PREM_DEFAULT;

    //Density profile PREM altered--------------------------------------------------------------------------

    PREM_DEFAULT.open("/home/dehy0499/OscProb/PremTables/prem_default.txt");
               
    std::ofstream PREM_LLSVP("/home/dehy0499/OscProb/PremTables/prem_llsvp.txt"); //Create a file with Alternative Earth density profile 
    
    // Variables for storing table rows
    float radius, density, zoa, layer;

    
    // Loop over table rows
    while(PREM_DEFAULT >> radius >> density >> zoa >> layer)
    {

        if( radius >=  R_cmb  && radius <= R_llsvp )
            {

                PREM_LLSVP << radius << " " << density*(1.0 + drho_dp/100) << " " << zoa << " " << layer << "\n";
                            

            }
                        
        else 
            {

                PREM_LLSVP << radius << " " << density << " " << zoa << " " << layer << "\n";

            }

    }

    PREM_DEFAULT.close();
    PREM_LLSVP.close();


    //-------------------------------------------------------------------------------------------------------------------

    prem_default   = "prem_default"; //Specify PREM table from OscProb

    TH2D* EventOsc_null = GetObservedEvents(prem_default, flvf,  Region, Bins , Det_par,  NnT);

    prem_llsvp   = "prem_llsvp"; //Specify PREM table from OscProb

    TH2D* EventOsc_alt = GetObservedEvents(prem_llsvp, flvf,  Region, Bins , Det_par,  NnT);


    /*
    std::cout << "Observed data is stored inside './SimulationResults' " << std::endl;

    TCanvas *c = new TCanvas();
    
    EventOsc->Draw("SURF1 Z");
    
    gStyle->SetPalette(55);

    gPad->SetTheta(30.0); // default is 30
    gPad->SetPhi(330.0); // default is 30 --Azimuthal
    gPad->Update();
    gPad->SetRightMargin(0.18);

    c->Print("SimulationResults/ObsEventsResults/test_poimeans.png");

    */

    return 0;

}
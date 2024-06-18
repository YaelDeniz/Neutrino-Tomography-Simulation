
/* Direction of neutrinos a re given in degrees*/
#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

//Cern ROOT
#include "TCanvas.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "Math/Interpolator.h"
//My Libraries
#include "MathTools.h"
#include "PhyTools.h"
#include "GetTrueEvents.h"


# define myPi 3.14159265358979323846  /* pi */
// Some Constants i need
# define R_earth 6368.0 //km
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds

using namespace std;


int main()
{
    // Plot and Histogram2D settings

    std::cout << " Neutrino Oscillation tomography. " << std::endl;

    // Number of bins
    
    int zbins=100 ; // # Bins in zenith/cos(zenith)
    int abins=100 ; // # Bins in azimuth
    int ebins=100 ; // bins in energy
    

    //Interval of integration in Zenith, Azimuth and E for Neutrino Events

    //Energy interval (in GeV)--------------------------------------------
    double Emin=1.0 ; 
    double Emax=10.0 ;

    //Zenith Angle Interval-----------------------------------------------

    //double Etamin = 33.0;
    //double Etamax = 38.0;

    // Alternative Earth Model:

    double h_llsvp = 800; //km
    double R_cmb = 3480.0;
    double R_llsvp = R_cmb + h_llsvp; //km

    double R_min = R_cmb-500.0;
    //double R_min = 3500.0; // Distance from the center of the Earth , 3500 Km CMB
    //double R_llsvp = R_cmb + h_llsvp; //Km
    double R_max = R_llsvp+500.0; // Distance from the center of the Earth

    
    double zenmax = 180-TMath::ASin( (R_min)/R_earth )*(180.0/TMath::Pi()) ; // max 180
    //double Etamax_LLSVP = TMath::ASin( (R_cmb + h_llsvp)/R_earth )*(180.0/TMath::Pi()) ;
    double zenmin = 180-TMath::ASin( R_max/R_earth )*(180.0/TMath::Pi()) ; // min 90

    double cz_min = cos(zenmax*TMath::Pi()/180.0);
    
    double cz_max = cos(zenmin*TMath::Pi()/180.0);

    std::cout << zenmin << " " << zenmax <<  std::endl;

    

    //double Etamin = 10 ;
    //double Etamax_LLSVP = TMath::ASin( (R_cmb + h_llsvp)/R_earth )*(180.0/TMath::Pi()) ;
    //double Etamax = 30 ;


    //MOdify prem to fit H_LLSVP

    //--------------------------------------------------------------------------------------------

    //double Etamin =10 ;
    //double Etamax =30 ;

    //Azimuthal Interal-------------------------------------------------------------------------
    double Phim = 0.0;
    double PhiM = 80.0 ;
    //double PhiM = 360.0 ; 

    //Detector size-----------------------------------------------------------------------------

    double DetMass = 10.0*MTon; //Mass in megaton units
    double Nn      = DetMass/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)
    double T       = 10.0*years2sec; //Detector Exposure time in sec: One Year

    double NnT = Nn*T; // Exposure Mton*years

    //Density contrast
    double drho_dp = 2.0; // 2 % density contrats

    
    

    //Simulation Information
      std::cout << "Energy Range in GeV: [" << Emin << " - " << Emax << "]" << "Angular Range(zenith): [" << zenmin << " - " <<zenmax << "]" <<std::endl;
      std::cout << "LLSVP information-  density contrats(%): " << drho_dp << " Height(km): "<< h_llsvp << std::endl;
      std::cout << "Simulation set up- Angular bins: " << Tbins << " Energy bins: "  << Ebins<<std::endl;
      std::cout << "PREM tables located in /OscProb/PremTables"<< std::endl;

    
    
    std::string prem_llsvp, prem_default;
    
    std::ifstream PREM_DEFAULT;

    //Density profile PREM altered--------------------------------------------------------------------------

    PREM_DEFAULT.open("/home/dehy0499/OscProb/PremTables/prem_44layers.txt");
               
   // std::ofstream PREM_LLSVP("/home/dehy0499/OscProb/PremTables/prem_llsvp.txt"); //Create a file with Alternative Earth density profile 
    
    // Variables for storing table rows
    float radius, density, zoa, layer;

    
    // Loop over table rows

    int nlayer = 1;
    std::vector <int> layers; 

    while(PREM_DEFAULT >> radius >> density >> zoa >> layer)
    {
            std::cout << " #layer: " << nlayer << std::endl;

        if( radius >=  R_cmb  && radius <= R_llsvp )
            {

                std::cout <<" LLVP segment:" << radius << " " << density << " " << zoa << " " << layer << "\n";
                std::cout <<" Remaning height " << R_llsvp - radius << std::endl;
                std::cout <<" Segment height : " << radius - 3480 << std::endl;
                layers.push_back(nlayer);   
                            

            }
                        
        else 
            {

//PREM_LLSVP << radius << " " << density << " " << zoa << " " << layer << "\n";

            }

        nlayer += 1;

    }

    PREM_DEFAULT.close();
   


   //Event generation-----------------------------------------------------------------------------------------
   int nuflv = 1; // neutrino  final state options nue (0), numu (1) or nutau (2)
    
   // Standard Earth model
    
    /*
    std::string prem_default   = "prem_default"; //Specify PREM table from OscProb
    TH2D* nullhist = AsimovTrueEvents(prem_default,false, layers ,flvf , Region,  Bins,  NnT);
    */

   AsimovSimulation StandardEarth;

   StandartEarth.PremModel = "prem_default";
   StandardEarth.MantleAnomaly = false;
   StandardEarth.SetIntervals(zenmin,zenmax,Phim,PhiM,EnuMin,EnuMax);
   StandardEarth.SetBinning(zbins,abins,ebins);
   StandardEarth.SetExposure(NnT);
   StandardEarth.flvf(nuflv);

   TH3D * TrueStd = StandardEarth.GetTrueEvents3D();

    // Alternative Earth Model

   AsimovSimulation AlternativeEarth;

   AlternativeEarth.PremModel = "prem_default";
   AlternativeEarth.MantleAnomaly = true;
   AlternativeEarth.AnomalousLayers = layers;
   AlternativeEarth.SetIntervals(zenmin,zenmax,Phim,PhiM,EnuMin,EnuMax);
   AlternativeEarth.SetBinning(zbins,abins,ebins);
   AlternativeEarth.SetExposure(NnT);
   AlternativeEarth.flvf(nuflv);

   TH3D * TrueAlt = AlternativeEarth.GetTrueEvents3D();



    //Difference in events

   TH3D* TrueDiff = new TH3D("TrueDiff","Percentage difference in  neutrino events",zbins,cz_min,cz_max,abins,Phim,PhiM,ebins,EnuMin,EnuMax); 

   
   // Data visualization
   
   std::ofstream EventDiff("SimulationResults/TrueEventsResults/3DSimulation.csv"); 
   double zen, azi, e, nexp, nobs, dn;

    for (int j = 1; i <= abins; j++)
    {
        phi = TrueDiff->GetYaxis()->GetBinCenter(j);
        
        for(int i=1; i<= zbins  ; i++) 
        {    
            cth = TrueDiff->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small

                        
                        for (int k=1; k <= ebins ; k++)
                        { 
                            e = TrueDiff->GetZaxis()->GetBinCenter(k); //< This will defined a constant L por different values of ct provided Dct is Small
                            
                            nexp = TrueStd->GetBinContent(i,j,k); // expected
                            nobs = TrueAlt->GetBinContent(i,j,k); // obserbed

                            dn = 100.0*(nobs-nexp)/nexp;

                            diffhist->SetBinContent(i,j, dn); //Create histogram for  kth Pseudo-Experiment

                            EventDiff<<  zen << ", " << e << ", "<< nexp << ", "<< nobs<< ", " << dn << "\n";

                            std::cout<<  zen << ", " << e << ", "<< nexp << ", "<< nobs<< ", " << dn << std::endl;


                            
                        } // loop energy

        } // Loop eta

    }

    EventDiff.close();


    TCanvas *c = new TCanvas();
    diffhist->Draw("COLZ");
    c->Print("SimulationResults/Histograms/3DModel.png");
    

    

    return 0;

}

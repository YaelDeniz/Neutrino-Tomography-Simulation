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
#include "Math/ProbFuncMathCore.h"
//My Libraries
#include "MathTools.h"
#include "PhyTools.h"
#include "GetObservedEvents.h"




# define myPi 3.14159265358979323846  /* pi */

# define R_earth 6368.0 //km

// Some Constants i need
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds

using namespace std;


int main()
{
    
    std::cout << " Neutrino Oscillation tomography. " << std::endl;
    std::cout << " Log-Likelihood Ratio Test (LLRT) / Chi^2 test for Observed events. "<< std::endl; 

    
     // Number of bins
    int Ebins=10; // # of Bins of Energy
    int Etabins=10; // # of Bins of cosEta

    int Ebins_o=10; // # of Bins of Energy
    int Etabins_o=10; // # of Bins of cosEta

    int Bins[]={Etabins, Ebins, Etabins_o, Ebins_o};
    //-----------------------------------------------------------------------------------------------------------------------------------------
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

    //---------------------------------------------------------------------------------------------------------------------------------------
    //LLSVP hypothesis


    double rho_min = -2; //More thermal hypothesis
    double rho_max = 2;  //More thermochemical hypothesis

   

    int m = 10; // # Number of points

    double rho_arr [m+1]; 
   

    double drho_arr = (rho_max-rho_min)/(m);



    for (int i = 0; i <= m; ++i)
    {

        rho_arr [i] = rho_min + i*drho_arr;
        
    }




    //--------------------------------------------------------------------------------------------


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
    
   


    // neutrino state options nue (0), numu (1) or nutau (2)
    
    int flvf = 1;


    //STATS--------------------------------------------------------------------------------------------

    // Get the model table paths ------------
    std::string Path2Original;
    std::ifstream StdPrem;



    std::string location = "SimulationResults/chi2results/" ;
    
    std::string title = "ObsDiff_Para_"+std::to_string(flvf)+"_"+std::to_string(Bins[0])+"_"+std::to_string(Bins[1])+"_"+std::to_string(Bins[0])+"_"+std::to_string(Bins[1])+"_"+std::to_string(Region[0])+"_"+std::to_string(Region[1])+".csv";
    std::string filename = location+title;




    ofstream ChiSqrtVals(filename, std::ofstream::trunc); //Opens a file and rewrite content, if files does not exist it Creates new file


    //Loops---------------------------------------------------------------------------------------------

        for (int i = 0; i <= m; ++i)
        {
                
                StdPrem.open("/home/dehy0499/OscProb/PremTables/prem_default.txt");
               
                std::ofstream AltPrem("/home/dehy0499/OscProb/PremTables/prem_alt.txt"); //Create a file with Alternative Earth density profile 
    

                double drho_dp = rho_arr[i];
                //double drho_dp = 2.0;


                std::cout << "density contrats:" << drho_dp << std::endl;


                // Variables for storing table rows
                float radius, density, zoa, layer;

                // Keep track of previous radius
                double rprev = 0;


                // Loop over table rows
                while(StdPrem >> radius >> density >> zoa >> layer){


                    if( radius > R_cmb  && radius < R_llsvp )
                    {


                        AltPrem << radius << " " << density*(1.0 + drho_dp/100) << " " << zoa << " " << layer << "\n";
                        

                    }
                    
                    else 
                    {

                        AltPrem << radius << " " << density << " " << zoa << " " << layer << "\n";

                    }

                }

                AltPrem.close();
                StdPrem.close();

                //---------------------------------------

                
                std::string prem_default, prem_alt;

                // Null hypothesis: Standard Earth
                
                prem_default   = "prem_default"; //Specify PREM table from OscProb

                // Alternative hypothesis

                prem_alt   = "prem_alt"; //Specify PREM table from OscProb

                TH2D* EventsExp = GetObservedEvents(prem_default, flvf,  Region, Bins , Det_par,  NnT); //Null hypothesis

                TH2D* EventsObs = GetObservedEvents(prem_alt, flvf,  Region, Bins , Det_par,  NnT); //Observed data from Alternative
                   
                //Chi2 Analysis: Log-likelihood Analysis.

                double dChi2 = 0;

                for(int m=1; m <= Etabins  ; m++) 
                {    

                    
                    for (int n = 1; n <= Ebins ; ++n)
                    { 
                        double nexp = EventsExp->GetBinContent(m,n); // expected
                        double nobs = EventsObs->GetBinContent(m,n); // obserbed

                        dChi2 +=  2*( nexp - nobs + nobs*TMath::Log(nobs/nexp) ); // LLRT

                        //dChi2 +=  (nexp - nobs)*(nexp - nobs)/(nexp); // Pearson's Chi^2

                        //dChi2 +=  (nexp - nobs)*(nexp - nobs)/(nobs); // Neyman's Chi^2


                        
                    } // loop energy

                } // Loop eta

                 //P-Value calculation

                int r = 1; //Degrees of freedorm

                //double pval_c = ROOT::Math::chisquared_cdf_c  (dChi2, r) ;
                //double pval = ROOT::Math::chisquared_cdf   (dChi2, r);  

                //std::cout << "exposure time" << T <<std::endl;

                std::cout << "Chi^2 Result = " << dChi2 << std::endl;       

                //std::cout << "Pvalue = " << pval_c << std::endl;

                std::cout << " " << std::endl;

                ChiSqrtVals << T <<" " << Det_par[0] << " " << Det_par[1]<< " "<< drho_dp << " "<< dChi2 << " "<< Ebins << " "<< Etabins <<"\n";

                delete EventsExp;
                delete EventsObs;

        }//End Loop



    ChiSqrtVals.close();

    return 0;

}
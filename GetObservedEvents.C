#include "GetObservedEvents.h"

//C++
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

//Cern Root
#include <math.h>
#include "TF1.h"
#include "TMath.h"
#include "TH2.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/Interpolator.h"
#include "Math/SpecFuncMathCore.h"

//#include "TGraphTrueD.h"

//OSCPROB
#ifndef __CINT__
#include "PremModel.h"
#include "PMNS_Fast.h"
bool isCINT = false;
#else
bool isCINT = true;
#endif

// My Math Tools
#include "MathTools.h"
#include "PhyTools.h"

using namespace std;


// Some Constants i need
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds

//Resolution functions

double w_E(double E , std::vector<double> Eo ); 

double w_th(double th, double E,std::vector<double> tho);



// Make oscillogram for given final flavour and MH
TH2D*  ObservedEvents(int flvf, double Energy[], double CosT[] ,int Ebins, int Tbins)
{

    ofstream ObservedData("SimulationResults/ObservedEvents.csv", std::ofstream::trunc); //Opens a file and rewrite content, if files does not exist it Creates new file
    
    //DETECTOR PROPERTIES
    double DetMass = 10*MTon; //Mass in megaton units
    double Nn      = DetMass/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)
    double T       = 10*years2sec; //Detector Exposure time in sec: One Year

    //Simulation Variables
    double N_ij    = 0;   //Mean number of events at bin ij. True distribution.
    double No_mn = 0;     //Mean number of events at bin mn. Observed distribution.
    double Ntot;

    //INTERVAL LIMITS FOR INTEGRATION:

    //Limits of Phi/Azimuthal angle
    double Phim    = 0.0;
    double PhiM    = 360 ; //DOlivo Azimuthal? for LLSVPs set 80
    
    //Integration limits of E
    double Emin      = Energy[0];
    double Emax      = Energy[1];
    
    //Range in Theta for Event Oscillogram
    double thmin   = 10.0;
    double thmax   = 30.0;
    
    double ctmin   = cos(thmax*TMath::Pi()/180); 
    double ctmax   = cos(thmin*TMath::Pi()/180);


    //True Bins:
    int ibins = 50; // Number of  angular bins of True event distribution
    int jbins = 50; // Number of  energy bins of True event distribution
    double dth = (thmax - thmin)/(2.0*ibins); //< Bin width/2
    double dE = (Emax - Emin)/(2.0*jbins); //< Bin width/2
    TH2D* hTrue = new TH2D("hTrue","True Events",jbins,thmin,thmax,ibins,Emin,Emax); // xbins correspond to energy values and ybins to zenith angle cost

    //Observed Bins:
    int mbins = Tbins; //Number of angular bins of Observed event distribution
    int nbins = Ebins; //Number of Energy bins of Observed event distribution
    double dtho = (thmax - thmin)/(2.0*mbins); //< Bin width/2
    double dEo = (Emax - Emin)/(2.0*nbins); //< Bin width/2
    TH2D* hObs = new TH2D("hObs","Observed Events",mbins,thmin,thmax,nbins,Emin,Emax); // xbins correspond to energy values and ybins to zenith angle cost

    
    
    //NEUTRINO OSCILLATION PROB

    // Neutrino flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino
    // Particle type
    int nu         = 1; //neutrino
    int nubar      = -1; //antineutrino

    // Model Set up: 
    OscProb::PMNS_Fast PMNS_Ho; // Create PMNS objects
    
    /*
    // Set specific Oscillation Parameters
    
    double dm21 = 7.42e-5;
    double dm31 = 2.533e-3;
    double th12 = 33.7712*TMath::Pi()/180;
    double th13 = 8.588*TMath::Pi()/180;
    double th23 = 48.504*TMath::Pi()/180;
    double dcp  = 214.2*TMath::Pi()/180;

    // Set PMNS parameters
    
    PMNS_Ho.SetDm(2, dm21);
    PMNS_Ho.SetDm(3, dm31);
    PMNS_Ho.SetAngle(1,2, th12);
    PMNS_Ho.SetAngle(1,3, th13);
    PMNS_Ho.SetAngle(2,3, th23);
    PMNS_Ho.SetDelta(1,3, dcp);
    */

    PMNS_Ho.SetStdPars(); // Set PDG 3-flavor parameters

    // Create PREM Model
    string test_profile;

    test_profile   = "/home/dehy0499/OscProb/PremTables/prem_test.txt"; //Specify PREM table from OscProb
    
    OscProb::PremModel prem(test_profile);

    // Atmospheric Neutrino Flux data:
    
    TFile *HF = new TFile("./NuFlux/Honda2014_spl-solmin-allavg.root","read"); //South Pole (IC telescope) Averaged all directions
    
    //TFile *HF = new TFile("./NuFlux/TestFlux.root","read"); // Specific flux V. Agrawal, et al. (1996)
    TDirectory *Zen;
    Zen = (TDirectory*) HF->Get("CosZ_all"); // Avg Flux for -0.9 <~ CosT < 0
    TTree *flux= (TTree*) Zen->Get("FluxData"); //Opens data file      
    
    double Enu,NuMu,NuMubar,NuE,NuEbar;
    flux->SetBranchAddress("Enu",&  Enu  );
    flux->SetBranchAddress("NuMu",&  NuMu );
    flux->SetBranchAddress("NuMubar",& NuMubar );
    flux->SetBranchAddress("NuE",& NuE );
    flux->SetBranchAddress("NuEbar",& NuEbar );

    //Interpolation of Neutrino flux data
    int n = flux->GetEntries();
    std::vector <double> EPsi,PsiNuMu,PsiNuMubar,PsiNuE,PsiNuEbar;
    for (int i = 0; i < n-1 ; ++i)
    {
        flux->GetEntry(i);
        EPsi.push_back(Enu);
        PsiNuMu.push_back(NuMu);
        PsiNuMubar.push_back(NuMubar);
        PsiNuE.push_back(NuE);
        PsiNuEbar.push_back(NuEbar);;
    }

    //Interpolate Muon-neutrino flux
    ROOT::Math::Interpolator dPsiMudE(EPsi,PsiNuMu, ROOT::Math::Interpolation::kLINEAR   );       // Muon Neutrino Flux Interpolation
    ROOT::Math::Interpolator dPsiMubardE(EPsi,PsiNuMubar, ROOT::Math::Interpolation::kLINEAR   ); // Muon Antineutrino Flux Interpolation

    //Interplate Electron-neutrinos flux
    ROOT::Math::Interpolator dPsiEdE(EPsi,PsiNuE, ROOT::Math::Interpolation::kLINEAR   );         // Electron Neutrino Flux Interpolation
    ROOT::Math::Interpolator dPsiEbardE(EPsi,PsiNuEbar, ROOT::Math::Interpolation::kLINEAR   );   // Electron Antineutrino Flux Interpolation
    
    // Calculation of Observed events:

   // std::vector<double> Th_o, Th, E_o, Et, dN;

    for (int m = 1; m <= mbins; ++m) //Observed Angular bins
    {

        double tho = hObs->GetXaxis()->GetBinCenter(m); // Select angular bin from the "Observed" histogram
        std::vector<double> Th_o;
        Th_o.push_back(tho - dtho);
        Th_o.push_back(tho);
        Th_o.push_back(tho + dtho);

        for (int n = 1; n <=nbins ; ++n) //Observed Energy bins
        {
     
            double Eo = hObs->GetYaxis()->GetBinCenter(n); // Select Energy bin from the "Observed" histogram
            std::vector<double> E_o;
            E_o.push_back(Eo - dEo);
            E_o.push_back(Eo);
            E_o.push_back(Eo + dEo);

            Ntot = 0;   // Reset total number of events.

            for(int i=1; i<= ibins ; ++i) //True Angular bins
            {    

                // Get cos(theta_z) from bin center, It fix a value of L
                double th = hTrue->GetXaxis()->GetBinCenter(i); // Select angular bin of the "True" histogram
                std::vector<double> Th;
                Th.push_back(th - dth);
                Th.push_back(th);
                Th.push_back(th + dth);

                double cosT =cos( (180-th)*TMath::Pi()/180 );

                if(cosT < -1 || cosT > 1) break; // Skip if cosT is unphysical 

                double L_Ho = prem.GetTotalL(cosT); // Set total path length L  

                prem.FillPath(cosT); // Fill paths from PREM model

                // Set paths in OscProb  
                PMNS_Ho.SetPath(prem.GetNuPath()); //STANDARD EARTH (DESCRIBED BY PREM MODEL)
            

                for (int j = 1; j <= jbins; ++j) //True Energy bins
                { 

                    double E = hTrue->GetYaxis()->GetBinCenter(j); // Select Energy bin from the "True" histogram
                    std::vector<double> Et;
                    Et.push_back(E - dE);
                    Et.push_back(E);
                    Et.push_back(E + dE);

                    // Numerical Integration of Energy dependet part

                    double hE = dE; //Step size for Integration of E

                    //std::vector<double> Th, Et, dN; // R= event rate d^2( N )/ (dtheta dE)
                    std::vector<double>  dN; // R= event rate d^2( N )/ (dtheta dE)

                    
                    for (int l = 0; l <= 2 ; ++l)
                    {   
                        
                        //Et.push_back(E + (l-1)*dE);
                        //Th.push_back(th + (l-1)*dth);
                        
                        //Resolution factor
                        //double res = w_E(Et[l],E_o)*( (-2.0*dth/3.0)*( w_th(Th[0],Et[l],Th_o)*sin(Th[0]) + 4*w_th(Th[1],Et[l],Th_o)*sin(Th[1]) +  w_th(Th[2],Et[l],Th_o)*sin(Th[2]) ) );
                        double res = w_E(Et[l],E_o)*( (2.0*dth/3.0)*( w_th(Th[0],Et[l],Th_o) + 4*w_th(Th[1],Et[l],Th_o) +  w_th(Th[2],Et[l],Th_o) ) );

                    
                        //Neutrino Contribution 
                        PMNS_Ho.SetIsNuBar(false); 
                        double r_ho = XSec(Et[l],nu)*( PMNS_Ho.Prob(numu, flvf, Et[l])*dPsiMudE.Eval(Et[l]) + PMNS_Ho.Prob(nue,flvf,Et[l])*dPsiEdE.Eval(Et[l]) );

                        //Antineutrino contribution
                        PMNS_Ho.SetIsNuBar(true); 
                        double r_hobar = XSec(Et[l],nubar)*( PMNS_Ho.Prob(numu,flvf, Et[l])*dPsiMubardE.Eval(Et[l]) + PMNS_Ho.Prob(nue,flvf,Et[l])*dPsiEbardE.Eval(Et[l]) ); 

                        //store data
                        dN.push_back( res*(r_ho + r_hobar) );

                    }
                    

                    //Integration in negery variables using NC quadrature
                    double NdTdM  =  simpson(Et,dN); //Integration of Event Rate Interaction for neutrinos passing Homogeneus mantle over E

                    //Number of events at bin(i,j) 
                    double N_ij=Nn*T*NdTdM;
                                    
                    Ntot += N_ij; 
                
                    //Et.clear();
                
                } // loop e
                
                //Th.clear();

            } // Loop th

                    
            No_mn = Ntot;

            ObservedData << tho << ", " << Eo << ", "<< No_mn << "\n"; //Write im file.
            
            hTrue->SetBinContent(m,n,No_mn); // Set events in Observed histogram. 

          //  E_o.clear();

        } // loop in eo

        //Th_o.clear();

    } // loop in tho

    //E_o.clear();
    //Et.clear();
    //Th_o.clear();
    //Th.clear();


    ObservedData.close();

    hObs->SetTitle("Observed events for #nu_{#mu} and #bar{#nu}_#mu ");
    
    hObs->SetStats(0);
            
    return hObs;
}


double w_E(double E , std::vector<double> Eo ) 
{

double aE = 0.1;

double sdE = aE*E;

double wE = (0.5)*( ROOT::Math::erf (  (Eo[2] - E)/( sqrt(2)*sdE )  ) - ROOT::Math::erf (  (Eo[0] - E)/( sqrt(2)*sdE )  ) );    

return wE;

}

double w_th(double th, double E, std::vector<double> tho)
{

double ath = 0.1;

double sdth = ath/(sqrt(E));

double wth = (0.5)*( ROOT::Math::erf (  (tho[2] - th)/( sqrt(2)*sdth )  ) - ROOT::Math::erf (  (tho[0] - th)/( sqrt(2)*sdth )  ) );   

return wth;
    
}
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

double w_E(double aE, double E , std::vector<double> Eo, double E_GeV[] );

double w_th(double aTh, double th, double E,std::vector<double> tho, double Eta[]);

// Make oscillogram for given final flavour and MH
TH2D*  GetObservedEvents(std::string model, int flvf, double E_GeV[], double Eta[],double dAz ,int Ebins, int Tbins, double Det_par[] )
{   

    //ofstream ObservedData("SimulationResults/NmuObs_extend.csv", std::ofstream::trunc); //Opens a file and rewrite content, if files does not exist it Creates new file
    
    //DETECTOR PROPERTIES
    //double DetMass = 10.0*MTon; //Mass in megaton units
    //double Nn      = DetMass/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)
    //double T       = 10.0*years2sec; //Detector Exposure time in sec: One Year

    double aE =  Det_par[0] ;
    double aTh = Det_par[1] ;
    double NnT = Det_par[2] ; //Exposure in MegaTon*Years

    //std::cout << "Detector expousre :" << DetMass*T << "MTon-Year" << std::endl;

    //Simulation Variables
    double N_ij    = 0.0;   //Mean number of events at bin ij. True distribution.
    double No_mn = 0.0;     //Mean number of events at bin mn. Observed distribution.
    double Ntot;
    double Ntrue;

    //INTERVAL LIMITS FOR INTEGRATION:
    
    //Integration limits of E
    double Emin      = E_GeV[0];
    double Emax      = E_GeV[1];
    
    double thmin   = Eta[0];
    double thmax   = Eta[1];


    double ctmin   = cos((180.0-thmax)*TMath::Pi()/180); 
    double ctmax   = cos((180.0-thmin)*TMath::Pi()/180);


    //True Bins:
    int ibins = Tbins; // Number of  angular bins of True event distribution
    int jbins = Ebins; // Number of  energy bins of True event distribution
    double dth = (thmax - thmin)/(2.0*ibins); //< Bin width/2
    double dE = (Emax - Emin)/(2.0*jbins); //< Bin width/2
    TH2D* hTrue = new TH2D("hTrue","True Events",jbins,thmin,thmax,ibins,Emin,Emax); // xbins correspond to energy values and ybins to zenith angle cost

    //Observed Bins:
    int mbins = Tbins; //Number of angular bins of Observed event distribution
    int nbins = Ebins; //Number of Energy bins of Observed event distribution
    double dtho = (thmax - thmin)/(2.0*mbins); //< Bin width/2
    double dEo = (Emax - Emin)/(2.0*nbins); //< Bin width/2
    TH2D* hObs = new TH2D("hObs","Observed Events",mbins,thmin,thmax,nbins,Emin,Emax); // xbins correspond to energy values and ybins to zenith angle cost

    //std::cout << "deltaE " << 2*dE <<" deltaTh"<< 2*dth <<std::endl; 
    
    //NEUTRINO OSCILLATION PROB

    // Neutrino flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino
    // Particle type
    int nu         = 1; //neutrino
    int nubar      = -1; //antineutrino

    // Model Set up: 
    OscProb::PMNS_Fast PMNS_Ho; // Create PMNS objects
    
    
    // Set specific Oscillation Parameters
    
    /*
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

    OscProb::PremModel prem(model);

    // Atmospheric Neutrino Flux data:
    
    TFile *HF = new TFile("./NuFlux/Honda2014_spl-solmin-allavg.root","read"); //South Pole (IC telescope) Averaged all directions
    
    //TFile *HF = new TFile("./NuFlux/TestFlux.root","read"); // Specific flux V. Agrawal, et al. (1996)
    //TDirectory *Zen;
    TDirectory *Zen = (TDirectory*) HF->Get("CosZ_all"); // Avg Flux for -0.9 <~ CosT < 0
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
    ROOT::Math::Interpolator dPsiMudE(EPsi,PsiNuMu, ROOT::Math::Interpolation::kCSPLINE);       // Muon Neutrino Flux Interpolation
    ROOT::Math::Interpolator dPsiMubardE(EPsi,PsiNuMubar, ROOT::Math::Interpolation::kCSPLINE); // Muon Antineutrino Flux Interpolation

    //Interplate Electron-neutrinos flux
    ROOT::Math::Interpolator dPsiEdE(EPsi,PsiNuE, ROOT::Math::Interpolation::kCSPLINE);         // Electron Neutrino Flux Interpolation
    ROOT::Math::Interpolator dPsiEbardE(EPsi,PsiNuEbar, ROOT::Math::Interpolation::kCSPLINE);   // Electron Antineutrino Flux Interpolation


   // Observed Events 
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
  

//----------------------------------------------------------------------------------------------
            
            for(int i=1; i<= ibins ; i++) //Loop in cosT
            {    

                // Get cos(theta_z) from bin center, It fix a value of L

                double th = hTrue->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
                
                double cosT =cos( (180.0-th)*TMath::Pi()/180.0 );

                if(cosT < -1 || cosT > 1) break; // Skip if cosT is unphysical 
                
                double L_Ho = prem.GetTotalL(cosT); // Set total path length L  
                
                prem.FillPath(cosT); // Fill paths from PREM model
                
                // Set paths in OscProb  
                PMNS_Ho.SetPath(prem.GetNuPath()); //STANDARD EARTH (DESCRIBED BY PREM MODEL)

                std::vector<double> Th;
                Th.push_back(th - dth);
                Th.push_back(th);
                Th.push_back(th + dth);


                for (int j = 1; j <= jbins; ++j)
                { 
                    //NUMERICAL INTEGRATION

                    double ene = hTrue->GetYaxis()->GetBinCenter(j); //< This will defined a constant L por different values of ct provided Dct is Small

                    //INTEGRATION IN ENERGY
                    int nstep = 2; //Steps for Integration of E

                    //double deltaE =2*dE;

                    //double hE = deltaE/nstep; //Step size for Integration of E

                    std::vector<double>  Et, R_Ho; // R= event rate d^2( N )/ (dtheta dE)

                    for (int l = 0; l <= nstep ; ++l)
                    {   
                        Et.push_back( ene + (l-1)*dE );
                        

                        double res = w_E( aE, Et[l],E_o,E_GeV)*(dth/3.0)*(  w_th( aTh, Th[0],Et[l],Th_o,Eta)*sin( (180.0-Th[0])*TMath::Pi()/180.0 )+ 4.0*w_th( aTh, Th[1],Et[l],Th_o,Eta)*sin( (180.0-Th[1])*TMath::Pi()/180.0 ) +  w_th( aTh, Th[2], Et[l], Th_o, Eta)*sin( (180.0-Th[2])*TMath::Pi()/180.0 )  );
                        
                        //Calculate the Interaction Rate at E for neutrinos passing trough lower mantle of Standar Earth

                        //Neutrino contribution
                        PMNS_Ho.SetIsNuBar(false); 
                        double r_ho = XSec(Et[l],nu)*( PMNS_Ho.Prob(numu, flvf, Et[l])*dPsiMudE.Eval(Et[l]) + PMNS_Ho.Prob(nue,flvf,Et[l])*dPsiEdE.Eval(Et[l]) );
                        
                        //Antineutrino contribution
                        PMNS_Ho.SetIsNuBar(true); 
                        double r_hobar = XSec(Et[l],nubar)*( PMNS_Ho.Prob(numu,flvf, Et[l])*dPsiMubardE.Eval(Et[l]) + PMNS_Ho.Prob(nue,flvf,Et[l])*dPsiEbardE.Eval(Et[l]) ); 
                        
                        //store data
                        R_Ho.push_back(res*(r_ho + r_hobar));

                    }
                    
            
                    //Integration in negery variables using NC quadrature
                    double dN_Ho_dOm = simpson(Et, R_Ho); //Integration of Event Rate Interaction for neutrinos passing Homogeneus mantle over E

                    //NUMBER OF EVENTS FOR [Emin Emax]&[cTmin cTmax] or bin(i,j) 
                    double N_ij=(NnT*(dAz)*(TMath::Pi()/180.0) )*dN_Ho_dOm;
                    Ntot += N_ij;

                    //ObservedData << th << ", " << ene << ", "<< N_ij << "\n";
                    
                } // loop e

            } // Loop ct
            
            No_mn = Ntot;



            //ObservedData << tho << ", " << Eo << ", "<< No_mn << "\n"; //Write im file.
            
            hObs->SetBinContent(m,n,No_mn); // Set events in Observed histogram. 
            

        } // loop in eo

    } // loop in tho

   //ObservedData.close();
       
   hObs->SetTitle("Observed events for #nu_{#mu} and #bar{#nu}_#mu ");
    
   hObs->SetStats(0);

   //Delete 

   delete HF;
   delete hTrue;

   return hObs;
}


double w_E(double aE, double E , std::vector<double> Eo, double E_GeV[] ) 
{

//double aE = 0.2;

double sdE = aE*E;

double AE = (0.5)*( ROOT::Math::erf (  (E - E_GeV[0] )/( sqrt(2)*sdE )  ) - ROOT::Math::erf (  (E - E_GeV[1] )/( sqrt(2)*sdE )  ) );  

double wE = (0.5)*( ROOT::Math::erf (  (E - Eo[0] )/( sqrt(2)*sdE )  ) - ROOT::Math::erf (  (E - Eo[2])/( sqrt(2)*sdE )  ) );    

return wE/AE;

}

double w_th(double aTh, double th, double E, std::vector<double> tho, double Eta[])
{

//double ath = 0.25;

double sdth = aTh/(sqrt(E));

double ATh = (180.0/TMath::Pi())*(0.5)*( ROOT::Math::erf (  (TMath::Pi()/180.0)*(th - Eta[0])/( sqrt(2)*sdth )  ) - ROOT::Math::erf (  (TMath::Pi()/180.0)*(th - Eta[1])/( sqrt(2)*sdth )  ) );

double wth = (180.0/TMath::Pi())*(0.5)*( ROOT::Math::erf (  (TMath::Pi()/180.0)*(th - tho[0])/( sqrt(2)*sdth )  ) - ROOT::Math::erf (  (TMath::Pi()/180.0)*(th - tho[2])/( sqrt(2)*sdth )  ) );   

return wth/ATh;
    
}

/*
 Changes:
    1.-HondaFlux was chaged to Flux Avaraged for al Azimuthal directions
    2.-Events account for neutrinos and Antineutrinos

*/
#include "GetTrueEvents.h"

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
#include <Math/Interpolator.h>
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

// Make oscillogram for given final flavour and MH
TH2D*  GetTrueEvents(std::string model, int flvf, double E_GeV[], double Eta[], double dAz ,int Ebins, int Tbins, double NnT)
{  



    ofstream TrueEvents("SimulationResults/test2.csv", std::ofstream::trunc); //Opens a file and rewrite content, if files does not exist it Creates new file
    
    double u       = 0;
    double N_ij    = 0;   //Mean number of events at bin ij.
    double Ntot    = 0;   // Total number of events.

   


    //std::cout << "Detector expousre :" << DetMass*T << "MTon-Year" << std::endl; 


    //INTERVAL LIMITS FOR INTEGRATION:

    //Limits of Phi/Azimuthal angle
    /*
    double Phim    = 0.0;
    double PhiM    = 80.0;
    double PhiM    = 360 ; //DOlivo Azimuthal?
    */

    //Integration limits of E
    double Emin      = E_GeV[0];
    double Emax      = E_GeV[1];

    //Range in Theta and E for Event Oscillogram
    //double thmin   = 10.0;
    //double thmax   = 30.0;

    double thmin   = Eta[0];
    double thmax   = Eta[1];

    double ctmin   = cos((180-thmax)*TMath::Pi()/180); 
    double ctmax   = cos((180-thmin)*TMath::Pi()/180);

    int ibins = Tbins; // Number of  angular bins of True event distribution
    int jbins = Ebins; // Number of  energy bins of True event distribution



    // Create 2D histogram for Event Oscillogram
    double dth = (thmax - thmin)/(2.0*ibins); //< Bin width/2
    double dE = (Emax - Emin)/(2.0*jbins); //< Bin width/2
    TH2D* hTrue = new TH2D("hTrue","True Events",ibins,thmin,thmax,jbins,Emin,Emax); // xbins correspond to energy values and ybins to zenith angle cost
    
    //NEUTRINO OSCILLATION PROB

    // Neutrino flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino
    // Particle type
    int nu         = 1; //neutrino
    int nubar      = -1; //antineutrino

    // Model Set up: 
    OscProb::PMNS_Fast PMNS_Ho; // Create PMNS objects
    
    
    // Set parameters to PDG
    
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

    // Create PREM Model
    //string file_test;
    //file_test   = "/home/dehy0499/OscProb/PremTables/prem_lbl.txt"; //Specify PREM table from OscProb
    //OscProb::PremModel prem(file_test);
    
    

    OscProb::PremModel prem(model);
    

    //OscProb::PremModel prem; //Default PREM table

    // Atmospheric Neutrino Flux data:
    TFile *HF = new TFile("./NuFlux/Honda2014_spl-solmin-allavg.root","read"); //South Pole (IC telescope)
    //TFile *HF = new TFile("./NuFlux/TestFlux.root","read"); //Flux Test
   
    // Set Avarege flux for a range of cosT
    TDirectory *Zen;
    Zen = (TDirectory*) HF->Get("CosZ_all"); // Avg Flux for -0.9 <~ CosT < -0.8
    TTree *flux= (TTree*) Zen->Get("FluxData"); //Opens data file      
    
    double Enu,NuMu,NuMubar,NuE,NuEbar;
    flux->SetBranchAddress("Enu",&  Enu  );
    flux->SetBranchAddress("NuMu",&  NuMu );
    flux->SetBranchAddress("NuMubar",& NuMubar );
    flux->SetBranchAddress("NuE",& NuE );
    flux->SetBranchAddress("NuEbar",& NuEbar );

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
    
    //Interpolate Neutrino flux data

    //Interpolate Muon-neutrino flux
    ROOT::Math::Interpolator dPsiMudE(EPsi,PsiNuMu, ROOT::Math::Interpolation::kCSPLINE    );       // Muon Neutrino Flux Interpolation
    ROOT::Math::Interpolator dPsiMubardE(EPsi,PsiNuMubar, ROOT::Math::Interpolation::kCSPLINE    ); // Muon Antineutrino Flux Interpolation

    //Interplate Electron-neutrinos flux
    ROOT::Math::Interpolator dPsiEdE(EPsi,PsiNuE, ROOT::Math::Interpolation::kCSPLINE    );         // Electron Neutrino Flux Interpolation
    ROOT::Math::Interpolator dPsiEbardE(EPsi,PsiNuEbar, ROOT::Math::Interpolation::kCSPLINE    );   // Electron Antineutrino Flux Interpolation
  

    for(int i=1; i<= ibins ; i++) //Loop in cosT
    {    

        // Get cos(theta_z) from bin center, It fix a value of L

        //double cosT = hTrue->GetYaxis()->GetBinCenter(ct); //< This will defined a constant L por different values of ct provided Dct is Small
        double th = hTrue->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
        double cosT =cos( (180-th)*TMath::Pi()/180 );

        if(cosT < -1 || cosT > 1) break; // Skip if cosT is unphysical 
        double L_Ho = prem.GetTotalL(cosT); // Set total path length L  
        prem.FillPath(cosT); // Fill paths from PREM model
        // Set paths in OscProb  
        PMNS_Ho.SetPath(prem.GetNuPath()); //STANDARD EARTH (DESCRIBED BY PREM MODEL)
    

        for (int j = 1; j <=jbins; ++j)
        { 
            //NUMERICAL INTEGRATION

            double e = hTrue->GetYaxis()->GetBinCenter(j); //< This will defined a constant L por different values of ct provided Dct is Small

            //INTEGRATION IN ENERGY
            int nstep = 2; //Steps for Integration of E

            double deltaE =2*dE;

            double hE = deltaE/nstep; //Step size for Integration of E

            std::vector<double>  E, R_Ho; // R= event rate d^2( N )/ (dtheta dE)
            for (int l = 0; l <= nstep ; ++l)
            {   
                
                E.push_back( e + (l-1)*dE);
               
                //Calculate the Interaction Rate at E for neutrinos passing trough lower mantle of Standar Earth
               
                //Neutrino contribution
                PMNS_Ho.SetIsNuBar(false); 
                double r_ho = XSec(E[l],nu)*( PMNS_Ho.Prob(numu, flvf, E[l])*dPsiMudE.Eval(E[l]) + PMNS_Ho.Prob(nue,flvf,E[l])*dPsiEdE.Eval(E[l]) );
                //Antineutrino contribution
                PMNS_Ho.SetIsNuBar(true); 
                double r_hobar = XSec(E[l],nubar)*( PMNS_Ho.Prob(numu,flvf, E[l])*dPsiMubardE.Eval(E[l]) + PMNS_Ho.Prob(nue,flvf,E[l])*dPsiEbardE.Eval(E[l]) ); 
                //store data
                R_Ho.push_back(r_ho + r_hobar);

            }
            

            //Integration in negery variables using NC quadrature
            double dN_Ho_dOm =  simpson(E, R_Ho); //Integration of Event Rate Interaction for neutrinos passing Homogeneus mantle over E

            // INTEGRATION IN THETA AND PHI (Solid Angle):
            double DOm = (180/TMath::Pi())*( cos( (180-(th+dth))*TMath::Pi()/180 ) - cos( (180-(th-dth))*TMath::Pi()/180 ) )*(dAz)*(TMath::Pi()/180) ; //Solid Angle =  DeltaCosT*DeltaPhi

            //NUMBER OF EVENTS FOR [Emin Emax]&[cTmin cTmax] or bin(i,j) 
            double N_ij=NnT*dN_Ho_dOm*DOm;
            
            TrueEvents << th << ", " << e << ", "<< N_ij << "\n";

            hTrue->SetBinContent(i,j,N_ij);      //Expected Events
        
            
        } // loop energy

    } // Loop eta

    TrueEvents.close();

    delete HF;
            
             
    return hTrue;
}

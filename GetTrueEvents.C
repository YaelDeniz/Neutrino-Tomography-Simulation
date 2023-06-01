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


//Delete later

// Define the PREM tables path
//#include "../prem_default.hpp"

// Macro to load OscProb library
//#include "LoadOscProb.C"

// Some functions to make nice plots
//#include "SetNiceStyle.C"

// My Math Tools
#include "MathTools.h"
#include "PhyTools.h"

using namespace std;


// Some Constants i need
//# define myPi 3.14159265358979323846  /* pi */
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds

// Make oscillogram for given final flavour and MH
TH2D*  GetTrueEvents(int flvf, double Energy[], double CosT[] ,int nbinsx, int nbinsy)
{  

    ofstream TrueEvents("SimulationResults/TrueEventsNCotes.csv", std::ofstream::trunc); //Opens a file and rewrite content, if files does not exist it Creates new file
    
    //Simulation Variables
    srand ( time(NULL) ); // Random Uniform Draw
    double u       = 0;
    double N_ij    = 0;   //Mean number of events at bin ij.
    double Nij_sim = 0;   //Random Sample from Poisson dist.
    double Ntot    = 0;   // Total number of events.
    double fij_exp = 0;   // Probability Density function.

    //INTERVAL LIMITS FOR INTEGRATION:

    //Limits of Phi/Azimuthal angle
    double Phim    = 0.0;
    //double PhiM    = 80.0;
    double PhiM    = 360 ; //DOlivo Azimuthal?
    //Integration limits of E
    double Emin      = Energy[0];
    double Emax      = Energy[1];
    double dE = (Emax - Emin)/(2.0*nbinsy); //< Bin width/2
    //Range in Theta and E for Event Oscillogram
    //double ctmin = CosT[0]; 
    //double ctmax = CosT[1];
    double thmin   = 10.0;
    double thmax   = 30.0;
    double dth = (thmax - thmin)/(2.0*nbinsx); //< Bin width/2
    double ctmin   = cos(thmax*TMath::Pi()/180); 
    double ctmax   = cos(thmin*TMath::Pi()/180);
    //double cosT = (ctmax + ctmin)/(2.0); // Mind Point
    //double dcT = (ctmax - ctmin)/(2.0*nbinsx); //< Bin width/2
    std::cout << "Region- E=( " << Emin <<"-"<< Emax <<" ), th=( "<< thmin <<"-"<< thmax <<" )" <<std::endl; 
    
    //DETECTOR PROPERTIES
    double DetMass = 10*MTon; //Mass in megaton units
    double Nn      = DetMass/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)
    double T       = 10*years2sec; //Detector Exposure time in sec: One Year
    
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
    string file_test;
    file_test   = "/home/dehy0499/OscProb/PremTables/prem_test.txt"; //Specify PREM table from OscProb
    OscProb::PremModel prem(file_test);

    //OscProb::PremModel prem; //Default PREM table
    
     //Create Energy vectors with  Logarithmic spaced points from Emin GeV to Emax GeV
   // std::vector<double> xbins = logspace(log10(Em),log10(EM),nbinsx + 1);
    //std::cout << " --- " <<  log10(Em) << "   " << log10(EM) << std::endl;
    //double Energies[nbinsx + 1];

    //for (int i = 0; i < nbinsx + 1; ++i)
    //{ 
        //Energies[i] = xbins[i]; 

    //    Energies[i] = Em + i*(EM-Em)/(nbinsx + 1); 

        // std::cout << "Check" << Energies[i] << std::endl; 
    //}

    std::cout<<" Azimuth angle averaged flux : Integration using Simpsons rule 1/3"<<std::endl;

    // Create 2D histogram for Event Oscillogram
    //TH2D* hTrue = new TH2D("hTrue","True Events",nbinsx, Energies ,nbinsy,ctmin,ctmax); // xbins correspond to energy values and ybins to zenith angle cost
    TH2D* hTrue = new TH2D("hTrue","True Events",nbinsx,thmin,thmax,nbinsy,Emin,Emax); // xbins correspond to energy values and ybins to zenith angle cost
    
    TFile *HF = new TFile("./NuFlux/Honda2014_spl-solmin-allavg.root","read"); //South Pole (IC telescope)
    
    //TFile *HF = new TFile("TestFlux.root","read"); //Flux Test
   
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
    ROOT::Math::Interpolator dPsiMudE(EPsi,PsiNuMu, ROOT::Math::Interpolation::kCSPLINE   );       // Muon Neutrino Flux Interpolation
    ROOT::Math::Interpolator dPsiMubardE(EPsi,PsiNuMubar, ROOT::Math::Interpolation::kCSPLINE   ); // Muon Antineutrino Flux Interpolation
    //Interplate Electron-neutrinos flux
    ROOT::Math::Interpolator dPsiEdE(EPsi,PsiNuE, ROOT::Math::Interpolation::kCSPLINE   );         // Electron Neutrino Flux Interpolation
    ROOT::Math::Interpolator dPsiEbardE(EPsi,PsiNuEbar, ROOT::Math::Interpolation::kCSPLINE   );   // Electron Antineutrino Flux Interpolation
    //END OF INTERPOLATION
    
    //Pseucdo Experiment  
    for(int ct=1; ct<= nbinsx ; ct++) //Loop in cosT
    {    

        // Get cos(theta_z) from bin center, It fix a value of L

        //double cosT = hTrue->GetYaxis()->GetBinCenter(ct); //< This will defined a constant L por different values of ct provided Dct is Small
        double t = hTrue->GetXaxis()->GetBinCenter(ct); //< This will defined a constant L por different values of ct provided Dct is Small
        double cosT =cos( (180-t)*TMath::Pi()/180 );

        if(cosT < -1 || cosT > 1) break; // Skip if cosT is unphysical 
        double L_Ho = prem.GetTotalL(cosT); // Set total path length L  
        prem.FillPath(cosT); // Fill paths from PREM model
        // Set paths in OscProb  
        PMNS_Ho.SetPath(prem.GetNuPath()); //STANDARD EARTH (DESCRIBED BY PREM MODEL)
    
        //NEUTRINO FLUX: Avarged over all Azimuth angles
        // Set Avarege flux for a range of cosT
        /*
        TDirectory *Zen;
        Zen = (TDirectory*) HF->Get("cosZ_-9_-8"); // Avg Flux for -0.9 <~ CosT < -0.8
        if(cosT >= -0.8 && cosT < -0.7) { Zen = (TDirectory*) HF->Get("cosZ_-8_-7"); } // Avg Flux for -0.8 <~ CosT < -0.7
        */

        //NEUTRINO FLUX: Avarged over all directions
        //TFile *HF = new TFile("Honda2014_spl-solmin-allavg.root","read"); //South Pole (IC telescope)

        for (int e = 1; e <=nbinsy; ++e)
        { 
            //NUMERICAL INTEGRATION

            double ene = hTrue->GetYaxis()->GetBinCenter(e); //< This will defined a constant L por different values of ct provided Dct is Small

            //INTEGRATION IN ENERGY
            int nstep = 100; //Steps for Integration of E

            double deltaE =2*dE;

            double hE = deltaE/nstep; //Step size for Integration of E

            std::vector<double>  Ei, R_Ho; // R= event rate d^2( N )/ (dtheta dE)
            for (int i = 0; i <= nstep ; ++i)
            {   
                
                Ei.push_back( (ene-dE) + i*hE);
                //Calculate the Interaction Rate at E for neutrinos passing trough lower mantle of Standar Earth
                //Neutrino contribution
                PMNS_Ho.SetIsNuBar(false); 
                double r_ho = XSec(Ei[i],nu)*( PMNS_Ho.Prob(numu, flvf, Ei[i])*dPsiMudE.Eval(Ei[i]) + PMNS_Ho.Prob(nue,flvf,Ei[i])*dPsiEdE.Eval(Ei[i]) );
                //Antineutrino contribution
                PMNS_Ho.SetIsNuBar(true); 
                double r_hobar = XSec(Ei[i],nubar)*( PMNS_Ho.Prob(numu,flvf, Ei[i])*dPsiMubardE.Eval(Ei[i]) + PMNS_Ho.Prob(nue,flvf,Ei[i])*dPsiEbardE.Eval(Ei[i]) ); 
                //store data
                R_Ho.push_back(r_ho + r_hobar);

            }
            
            /*

            PMNS_Ho.SetIsNuBar(false);
            double r_nobar = XSec(ene,nu)*( PMNS_Ho.Prob(numu, flvf, ene)*dPsiMudE.Eval(ene) + PMNS_Ho.Prob(nue,flvf,ene)*dPsiEdE.Eval(ene) );
            //Antineutrino contribution
            PMNS_Ho.SetIsNuBar(true);
            double r_bar = XSec(ene,nubar)*( PMNS_Ho.Prob(numu,flvf,ene)*dPsiMubardE.Eval(ene) + PMNS_Ho.Prob(nue,flvf,ene)*dPsiEbardE.Eval(ene) ); 
            double f = r_nobar + r_bar;
            
            */

            //Integration in negery variables using NC quadrature
            double dN_Ho_dOm =  ncquad(Ei, R_Ho); //Integration of Event Rate Interaction for neutrinos passing Homogeneus mantle over E

            // INTEGRATION IN THETA AND PHI (Solid Angle):
            double cTmin =cos( (180-(t-dth))*TMath::Pi()/180 );
            double cTmax =cos( (180-(t+dth))*TMath::Pi()/180 );
            double DOm = (cTmax-cTmin)*(PhiM - Phim)*(TMath::Pi()/180) ; //Solid Angle =  DeltaCosT*DeltaPhi

            //NUMBER OF EVENTS FOR [Emin Emax]&[cTmin cTmax] or bin(i,j) 
            double N_ij=Nn*T*dN_Ho_dOm*DOm;
            Ntot += N_ij;
            fij_exp = N_ij;
            //fij_exp = N_ij;
            //double N_ij=f;


            //Random Sample from Pois(N_ij)
            //u = static_cast<double>(rand()) / RAND_MAX; //Uniform Distributed?
            //Nij_sim = Sample_Pois(N_ij, u);

            //hTrue->SetBinContent(e,ct, Nij_sim ); //Sample Events
            //std::cout <<"Bin id (" << ct << "," << e << "): E:("<< Energies[e-1] << "-" << Energies[e] << ") N= "  << N_ij << " sim: "<< Nij_sim<<std::endl;

            //if(N_ij < 0) { N_ij = 0; }

            TrueEvents << t << ", " << ene << ", "<< fij_exp << "\n";

            hTrue->SetBinContent(ct,e,fij_exp);      //Expected Events
            //std::cout <<"Bin id (" << ct << "," << e << "): E[ " << ene-dE <<"-"<< ene+dE <<" ] th[" << t-dth <<"-"<< t+dth << "] " << " DeltaC: "<< cTmax-cTmin << "- DeltaE: " << deltaE << "| N= "  << N_ij << std::endl;
            //std::cout <<"Bin id (" << ct << "," << e << "): E:("<< Energies[e-1] << "-" << Energies[e] << ") N= "  << N_ij << std::endl;
            
        } // loop e

    } // Loop ct

    TrueEvents.close();



    //Total Events

    //INTEGRATION IN ENERGY
            int nstep = 1000; //Steps for Integration of E
            double deltaE =Emax - Emin;
            double hE = deltaE/nstep; //Step size for Integration of E
            std::vector<double>  Etot, R_tot; // R= event rate d^2( N )/ (dtheta dE)
            for (int i = 0; i <= nstep ; ++i)
            {   
                Etot.push_back(Emin + i*hE);
                //Calculate the Interaction Rate at E for neutrinos passing trough lower mantle of Standar Earth
                //Neutrino contribution
                PMNS_Ho.SetIsNuBar(false); 
                double r_tot = XSec(Etot[i],nu)*( PMNS_Ho.Prob(numu, flvf, Etot[i])*dPsiMudE.Eval(Etot[i]) + PMNS_Ho.Prob(nue,flvf,Etot[i])*dPsiEdE.Eval(Etot[i]) );
                //Antineutrino contribution
                PMNS_Ho.SetIsNuBar(true); 
                double r_totbar = XSec(Etot[i],nubar)*( PMNS_Ho.Prob(numu,flvf, Etot[i])*dPsiMubardE.Eval(Etot[i]) + PMNS_Ho.Prob(nue,flvf,Etot[i])*dPsiEbardE.Eval(Etot[i]) ); 
                //store data
                R_tot.push_back(r_tot + r_totbar);
            }


            //Integration in negery variables using NC quadrature
            double dN_tot_dOm =  ncquad(Etot, R_tot); //Integration of Event Rate Interaction for neutrinos passing Homogeneus mantle over E

            // INTEGRATION IN THETA AND PHI (Solid Angle):
            double ctmim =cos(170*TMath::Pi()/180 );
            double ctmaxm =cos(150*TMath::Pi()/180 );
            double Dotot = (ctmaxm-ctmim)*(PhiM - Phim)*(TMath::Pi()/180) ; //Solid Angle =  DeltaCosT*DeltaPhi





    std::cout << "N total= " << Ntot<< " -- " << Nn*T*dN_tot_dOm*Dotot<< std::endl; 
    //hTrue->Scale(1. / hTrue->Integral(), "width");
    hTrue->SetTitle("Expected events for #nu_{#mu} and #bar{#nu}_#mu ");
    //hTrue->GetYaxis()->SetTitle("E (GeV)");
    //hTrue->GetXaxis()->SetTitle("eta");
    //hTrue->GetZaxis()->SetTitle("Events per bin");
    hTrue->SetStats(0);
            
             
    return hTrue;
}

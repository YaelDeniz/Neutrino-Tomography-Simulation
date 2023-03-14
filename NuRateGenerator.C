/*
 Changes:
    1.-HondaFlux was chaged to Flux Avaraged for al Azimuthal directions
    2.-Events account for neutrinos and Antineutrinos

*/
#include "NuRateGenerator.h"

#include <iostream>
#include "TF1.h"
#include "TMath.h"
//#include "TGraph.h"
//#include "TMultiGraph.h"
//#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include <Math/Interpolator.h>
//#include "TGraph2D.h"

//#ifndef __CINT__
#include "PremModel.h"
#include "PMNS_Fast.h"
//bool isCINT = false;
//#else
//bool isCINT = true;
//#endif

// Define the PREM tables path
//#include "../prem_default.hpp"

// Macro to load OscProb library
//#include "LoadOscProb.C"

// Some functions to make nice plots
//#include "SetNiceStyle.C"

// My Math Tools
#include "MathTools.h"
#include "PhyTools.h"



// Some Constants i need
# define myPi 3.14159265358979323846  /* pi */
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds

 
// Make oscillogram for given final flavour and MH
double Nflvf_ij(int flvf, double CosT[] , double Energy[])
{
    
    // Load OscProb library
  //  if(isCINT) LoadOscProb();   

    //INTERVAL LIMITS FOR INTEGRATION:

    //Integration limits of Phi/Azimuthal angle
    double Phim= 0.0;
    double PhiM =80.0;
    //Integration limits of E
    double Em = Energy[0];
    double EM = Energy[1];
    //Range in Theta and E for Event Oscillogram
    double ctmin = CosT[0]; 
    double ctmax = CosT[1];
    double cTmid = (ctmax + ctmin)/(2.0); // Mind Point

    

    // CONSTANTS
    //Detector properties
    double DetMass = 1*MTon;
    double Nn = DetMass/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)
    double T = 1*years2sec; //Detector Exposure time in sec: One Year
    // Neutrino flavour
    int nue = 0;  // electron neutrino  
    int numu = 1; // muon neutrino
    // Particle type
    int nu = 1; //neutrino
    int nubar= -1; //antineutrino
    // END CONSTANTS




    //OSCILLATION PROBABILITY CALCULATOR

    //Model Set up: 
    OscProb::PMNS_Fast PMNS_Ho; // Create PMNS objects
    PMNS_Ho.SetStdPars(); // Set PDG 3-flavor parameters
    
    // Create default PREM Model
    //std::string file_prem;
    //file_prem = prem; //Use prem_default.txt as Earth Model
    OscProb::PremModel prem;
    //OscProb::PremModel prem(file_prem); //Set eviroment for neutrino propagation: Neutrinos traveling through Earth
    double L_Ho = prem.GetTotalL(cTmid); // Set total path length L  
    // Skip if cosT is unphysical  
    //if(cTmid < -1 || cTmid > 1) break;
    prem.FillPath(cTmid); // Fill paths from PREM model
    // Set paths in OscProb  
    PMNS_Ho.SetPath(prem.GetNuPath()); //STANDARD EARTH (DESCRIBED BY PREM MODEL)
        
    



    //NEUTRINO FLUX: Avarged over all Azimuth angles

    TFile *HF = new TFile("Honda2014_spl-solmin-avgaz.root","read"); //South Pole (IC telescope)
    // Set Avarege flux for a range of cosT
    TDirectory *Zen;
    Zen = (TDirectory*) HF->Get("cosZ_-9_-8"); // Avg Flux for -0.9 <~ CosT < -0.8
    if(cTmid >= -0.8 && cTmid < -0.7) { Zen = (TDirectory*) HF->Get("cosZ_-8_-7"); /* Avg Flux for -0.8 <~ CosT < -0.7*/ } 
    
    //INTERPOLATE FLUX DATA POINTS

    TTree *flux= (TTree*) Zen->Get("FluxData"); //Opens data file
            
    // File data
    double Enu,NuMu,NuMubar,NuE,NuEbar;
                        
    flux->SetBranchAddress("Enu",&  Enu  );
    flux->SetBranchAddress("NuMu",&  NuMu );
    flux->SetBranchAddress("NuMubar",& NuMubar );
    flux->SetBranchAddress("NuE",& NuE );
    flux->SetBranchAddress("NuEbar",& NuEbar );

    int n = flux->GetEntries();
        
    std::vector <double> EPsi,PsiNuMu,PsiNuMubar,PsiNuE,PsiNuEbar;
        
    for (int i = 0; i < n ; ++i)
    {
            flux->GetEntry(i);
            EPsi.push_back(Enu);
            PsiNuMu.push_back(NuMu);
            PsiNuMubar.push_back(NuMubar);
            PsiNuE.push_back(NuE);
            PsiNuEbar.push_back(NuEbar);
    }

    //Interpolate Muon-neutrino flux
    ROOT::Math::Interpolator dPsiMudE(EPsi,PsiNuMu, ROOT::Math::Interpolation::kCSPLINE );       // Muon Neutrino Flux Interpolation
    ROOT::Math::Interpolator dPsiMubardE(EPsi,PsiNuMubar, ROOT::Math::Interpolation::kCSPLINE ); // Muon Antineutrino Flux Interpolation

    //Interplate Electron-neutrinos flux
    ROOT::Math::Interpolator dPsiEdE(EPsi,PsiNuE, ROOT::Math::Interpolation::kCSPLINE );         // Electron Neutrino Flux Interpolation
    ROOT::Math::Interpolator dPsiEbardE(EPsi,PsiNuEbar, ROOT::Math::Interpolation::kCSPLINE );   // Electron Antineutrino Flux Interpolation
       
    //END OF INTERPOLATION

    
    //NUMERICAL INTEGRATION

    //INTEGRATION IN ENERGY
    int nstep = 200; //Steps for Integration of E
    double hE = (EM - Em)/nstep; //Step size for Integration of E
    std::vector<double>  Ei, R_Ho; // R= event rate d^2( N )/ (dtheta dE)
    for (int i = 0; i <= nstep ; ++i)
    {   
        
        Ei.push_back(Em + i*hE);
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
    //Integration in negery variables using NC quadrature
    double dN_Ho_dOm =  ncquad(Ei, R_Ho); //Integration of Event Rate Interaction for neutrinos passing Homogeneus mantle over E
    // INTEGRATION IN THETA AND PHI (Solid Angle):
    double DOm = (ctmax - ctmin) * (PhiM - Phim)*(myPi/180.0) ; //Solid Angle =  DeltaCosT*DeltaPhi


    //NUMBER OF EVENTS FOR [Emin Emax]&[cTmin cTmax] or bin(i,j) 

    double Nij_Ho= Nn*T*(dN_Ho_dOm)*DOm;
                       
    return Nij_Ho;
}

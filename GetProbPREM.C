/*
 Changes:
    1.-HondaFlux was chaged to Flux Avaraged for al Azimuthal directions
    2.-Events account for neutrinos and Antineutrinos

*/
#include "GetProbPREM.h"

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

//OSCPROB
#ifndef __CINT__
#include "PremModel.h"
#include "PMNS_Fast.h"
bool isCINT = false;
#else
bool isCINT = true;
#endif


using namespace std;

// Make oscillogram for given final flavour and MH
TH2D*  GetProbPREM(int flvi, int flvf, double E_GeV[], double Eta[],int Ebins, int Tbins)
{  



    ofstream Pab("SimulationResults/Delta_Prob_LLSVP.csv", std::ofstream::trunc); //Opens a file and rewrite content, if files does not exist it Creates new file
    

    double P_ij    = 0;   //Mean number of events at bin ij.
    double P_ijbar    = 0;   //Mean number of events at bin ij.

    //Integration limits of E
    double Emin      = E_GeV[0];
    double Emax      = E_GeV[1];

    double thmin   = Eta[0];
    double thmax   = Eta[1];

    double ctmin   = cos((180-thmax)*TMath::Pi()/180); 
    double ctmax   = cos((180-thmin)*TMath::Pi()/180);

    int ibins = Tbins; // Number of  angular bins of True event distribution
    int jbins = Ebins; // Number of  energy bins of True event distribution



    // Create 2D histogram for Event Oscillogram
    double dth = (thmax - thmin)/(2.0*ibins); //< Bin width/2
    double dE = (Emax - Emin)/(2.0*jbins); //< Bin width/2
    TH2D* dPab = new TH2D("dPab","Oscillation Probabilities",ibins,thmin,thmax,jbins,Emin,Emax); // xbins correspond to energy values and ybins to zenith angle cost
    
    //NEUTRINO OSCILLATION PROB

    // Model Set up: 

    OscProb::PMNS_Fast PMNS_Ho;
    OscProb::PMNS_Fast PMNS_H; // Create PMNS objects
    
    
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
    PMNS_H.SetStdPars(); // Set PDG 3-flavor parameters

    //Create PREM Model
    //string file_test;
    //file_test   = "/home/dehy0499/OscProb/PremTables/prem_lbl.txt"; //Specify PREM table from OscProb
    //OscProb::PremModel prem(file_test);
    
    OscProb::PremModel PREM; //Default PREM table
    
    std::string prem_llsvp;
    prem_llsvp   = "/home/dehy0499/OscProb/PremTables/Prem_LLSVPlow.txt"; //Specify PREM table from OscProb
    OscProb::PremModel PREM_llsvp(prem_llsvp);
    

    for(int i=1; i<= ibins ; i++) //Loop in cosT
    {    

        // Get cos(theta_z) from bin center, It fix a value of L

        //double cosT = dPab->GetYaxis()->GetBinCenter(ct); //< This will defined a constant L por different values of ct provided Dct is Small
        double th = dPab->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
        double cosT =cos( (180-th)*TMath::Pi()/180 );

        if(cosT < -1 || cosT > 1) break; // Skip if cosT is unphysical 
        
        double L_Ho = PREM.GetTotalL(cosT); // Set total path length L
        PREM.FillPath(cosT); // Fill paths from PREM model

        double L_H  = PREM_llsvp.GetTotalL(cosT);  
        PREM_llsvp.FillPath(cosT);

        // Set paths in OscProb  
        PMNS_Ho.SetPath(PREM.GetNuPath()); //STANDARD EARTH (DESCRIBED BY PREM MODEL)
        PMNS_H.SetPath(PREM_llsvp.GetNuPath()); //STANDARD EARTH (DESCRIBED BY PREM MODEL)
    

        for (int j = 1; j <=jbins; ++j)
        { 
            //NUMERICAL INTEGRATION

            double e = dPab->GetYaxis()->GetBinCenter(j); //< This will defined a constant L por different values of ct provided Dct is Small

            PMNS_Ho.SetIsNuBar(false); 
            PMNS_H.SetIsNuBar(false); 

            P_ij = PMNS_H.Prob(flvi,flvf,e) - PMNS_Ho.Prob(flvi,flvf,e);

            PMNS_Ho.SetIsNuBar(true); 
            PMNS_H.SetIsNuBar(true); 

            P_ijbar = PMNS_H.Prob(flvi,flvf,e) - PMNS_Ho.Prob(flvi,flvf,e);
            
            Pab << th << ", " << e << ", "<< P_ij << "," << P_ijbar <<"\n";

            dPab->SetBinContent(i,j,P_ij);      //Expected Events
        
            
        } // loop energy

    } // Loop eta

    Pab.close();
            
             
    return dPab;
}

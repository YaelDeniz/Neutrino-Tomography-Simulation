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

//Random number generator
#include <random>


//Cern Root
#include <math.h>
#include "TF1.h"
#include "TMath.h"
#include "TH2.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
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
TObjArray*  GetTrueEventsMC(std::string modelname, int flvf, double Region[], int Bins[], double NnT, int K)
{  

    double total = 0;

    //Data Storage -----------------------------------------------------------------------------------------------------
    std::string model = "/home/dehy0499/OscProb/PremTables/"+modelname+".txt";
    
    std::cout << modelname <<std::endl;

    std::string location = "SimulationResults/PoisMeans/" ;
    
    std::string Earthmodel= modelname;
    
    std::string title = "_PoiMeansT_"+std::to_string(flvf)+"_"+std::to_string(Bins[0])+"_"+std::to_string(Bins[1])+"_"+std::to_string(K)+".csv";
    
    std::string filename = location+Earthmodel+title;
    
    ofstream TrueEvents(filename, std::ofstream::trunc); //Opens a file and rewrite content, if files does not exist it Creates new file

    TObjArray *PseudoExp = new TObjArray(K+2);

    // Monte Carlo Simulation - Random Sampling-------------------------------------------------------------------------

    std::random_device rd; //produce non-deterministic random numbers.
    
    std::mt19937 gen(rd()); //high quality, but not cryptographically secure, unsigned integer random numbers of type UIntType on the interval
    
    double N_ij = 0   ;   //Poisson mean for bin ij.
    double n_ij  = 0   ;   //Random Sample for bin ij
    double n_ijk  =0   ;   //Random Sample for bin ij; Pseudo experiment k
    double Ntot  = 0  ;   //Total number of events.
    int psudoexp   = K;   //Montecarlo simulation.


    //Binnig scheme and Oscillogram-------------------------------------------------------------------------------------

    //Energy Intervals
    double Emin      = Region[0];//Lower limit for Energy
    double Emax      = Region[1];//Upper limit for Energy

    //Zenith Angle Interval
    double etamin   = Region[2];//Lower limit for angle
    double etamax   = Region[3];//Upperlimit for angle
    double cetamin   = cos((180-etamax)*TMath::Pi()/180); 
    double cetamax   = cos((180-etamin)*TMath::Pi()/180);

    //Phi/Azimuthal Intervals
    double dAz = Region[4];

    //Bins
    int ibins = Bins[0]; // Number of  angular bins of True event distribution
    int jbins = Bins[1]; // Number of  energy bins of True event distribution
    double deta = (etamax - etamin)/(2.0*ibins); //< Bin width/2
    double dE = (Emax - Emin)/(2.0*jbins); //< Bin width/2

    /* Create 2D histogram for Event Oscillogram,
     xbins correspond to energy values and ybins to zenith angle cosEta*/
    TH2D* hTrue = new TH2D("hTrue","True Events",ibins,etamin,etamax,jbins,Emin,Emax); 
    TH2D* hTrue_k= new TH2D("hTrue_k","True Events; #eta ; #E",ibins,etamin,etamax,jbins,Emin,Emax); //Store data for each pseudo-experiment.
    TH2D* hTrue_means= new TH2D("hTrue_means","True Events( Means ); #eta ; #E",ibins,etamin,etamax,jbins,Emin,Emax); //Store data for each pseudo-experiment.
    
    //Neutrino event generation-----------------------------------------------------------------------------------------

    // Neutrino flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino
    // Particle type
    int nu         = 1; //neutrino
    int nubar      = -1; //antineutrino

    //Neutrino Oscillation Probabilities calculation--------------------------------------------------------------------
    OscProb::PMNS_Fast PMNS_H; // Create PMNS objects
    

    /*
    // Set parameters to PDG
    double dm21 = 7.42e-5;
    double dm31 = 2.533e-3;
    double th12 = 33.7712*TMath::Pi()/180;
    double th13 = 8.588*TMath::Pi()/180;
    double th23 = 48.504*TMath::Pi()/180;
    double dcp  = 214.2*TMath::Pi()/180;
    // Set PMNS parameters
    PMNS_H.SetDm(2, dm21);
    PMNS_H.SetDm(3, dm31);
    PMNS_H.SetAngle(1,2, th12);
    PMNS_H.SetAngle(1,3, th13);
    PMNS_H.SetAngle(2,3, th23);
    PMNS_H.SetDelta(1,3, dcp);
    */
    
    PMNS_H.SetStdPars(); // Set PDG 3-flavor parameters

    // Create PREM Model
    //string file_test;
    //file_test   = "/home/dehy0499/OscProb/PremTables/prem_lbl.txt"; //Specify PREM table from OscProb
    //OscProb::PremModel prem(file_test);
    
    
    //Earth model for neutrino propagation 
    OscProb::PremModel prem(model);
    //OscProb::PremModel prem; //Default PREM table

    // Atmospheric Neutrino Flux data:
    TFile *HF = new TFile("./NuFlux/Honda2014_spl-solmin-allavg.root","read"); //South Pole (IC telescope)
    //TFile *HF = new TFile("./NuFlux/TestFlux.root","read"); //Flux Test
   
    // Set Avarege flux for a range of cosEta
    TDirectory *Zen;
    Zen = (TDirectory*) HF->Get("CosZ_all"); // Avg Flux for -0.9 <~ CosEta < -0.8
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
        PsiNuEbar.push_back(NuEbar);
    }
    
    //Interpolate Neutrino flux data:
    //Muon-neutrino flux
    ROOT::Math::Interpolator dPsiMudE(EPsi,PsiNuMu, ROOT::Math::Interpolation::kCSPLINE);// Muon Neutrino Flux Interpolation
    ROOT::Math::Interpolator dPsiMubardE(EPsi,PsiNuMubar, ROOT::Math::Interpolation::kCSPLINE);// Muon Antineutrino Flux Interpolation
    //Electron-neutrinos flux
    ROOT::Math::Interpolator dPsiEdE(EPsi,PsiNuE, ROOT::Math::Interpolation::kCSPLINE);// Electron Neutrino Flux Interpolation
    ROOT::Math::Interpolator dPsiEbardE(EPsi,PsiNuEbar, ROOT::Math::Interpolation::kCSPLINE);// Electron Antineutrino Flux Interpolation

    //Monte Carlo Simulation--------------------------------------------------------------------------------------------  
    for (int k = 1; k <= psudoexp; ++k) //kth sudo Experiments for the MC simulation
    {
        //std::cout << "Pseudo-exp #" << k<<std::endl;

        total = 0;

        for(int i=1; i<= ibins ; i++) //Loop in Angular bins
        {    

            // Get cos(eta) from bin center, This is used to calculate the baseline.

            double eta = hTrue->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
            

            double cosEta =cos( (180-eta)*TMath::Pi()/180 );

            if(cosEta < -1 || cosEta > 1) break; // Skip if cosEta is unphysical 
            double L_Ho = prem.GetTotalL(cosEta); // Set total path length L  
            prem.FillPath(cosEta); // Fill paths from PREM model
            PMNS_H.SetPath(prem.GetNuPath()); // Set paths in OscProb  
        
            for (int j = 1; j <=jbins; ++j)
            { 
                double e = hTrue->GetYaxis()->GetBinCenter(j); //< This will defined a constant L por different values of ct provided Dct is Small
                
                
                //NUMERICAL INTEGRATION
                
                int nstep = 2; //Steps for Integration of E
                double deltaE =2*dE;
                double hE = deltaE/nstep; //Step size for Integration of E
                std::vector<double>  E, R_nu; // R= event rate d^2( N )/ (dtheta dE)
                
                for (int l = 0; l <= nstep ; ++l)
                {   
                    
                    E.push_back( e + (l-1)*dE);
                   
                    //Calculate the Interaction Rate at E for neutrinos passing trough lower mantle of Standar Earth
                   
                    //Neutrino contribution
                    PMNS_H.SetIsNuBar(false); 
                    double r_nu = XSec(E[l],nu)*( PMNS_H.Prob(numu, flvf, E[l])*dPsiMudE.Eval(E[l]) + PMNS_H.Prob(nue,flvf,E[l])*dPsiEdE.Eval(E[l]) );
                    //Antineutrino contribution
                    PMNS_H.SetIsNuBar(true); 
                    double r_nubar = XSec(E[l],nubar)*( PMNS_H.Prob(numu,flvf, E[l])*dPsiMubardE.Eval(E[l]) + PMNS_H.Prob(nue,flvf,E[l])*dPsiEbardE.Eval(E[l]) ); 
                    //store data
                    R_nu.push_back(r_nu + r_nubar);

                }
                

                //Integration in negery variables using NC quadrature
                double dN_nu_dOm =  simpson(E, R_nu); //Integration of Event Rate Interaction for neutrinos passing Homogeneus mantle over E

                // INTEGRATION IN THETA AND PHI (Solid Angle):
                double DOm = (180/TMath::Pi())*( cos( (180-(eta+deta))*TMath::Pi()/180 ) - cos( (180-(eta-deta))*TMath::Pi()/180 ) )*(dAz)*(TMath::Pi()/180) ; //Solid Angle =  DeltaCosEta*DeltaPhi

                //NUMBER OF EVENTS FOR [Emin Emax]&[cetamin cetamax] or bin(i,j) 

                N_ij=NnT*dN_nu_dOm*DOm; //Mean of the Monte carlo simulation 

                
                
                //Appoximation for large bin size
                
                /*
                double r_nu = XSec(e,nu)*( PMNS_H.Prob(numu, flvf, e)*dPsiMudE.Eval(e) + PMNS_H.Prob(nue,flvf,e)*dPsiEdE.Eval(e) );
                double r_nubar = XSec(e,nubar)*( PMNS_H.Prob(numu,flvf, e)*dPsiMubardE.Eval(e) + PMNS_H.Prob(nue,flvf,e)*dPsiEbardE.Eval(e) ); 
                N_ij=NnT*(r_nu*r_nubar)*(2*dE)*(2*deta*TMath::Pi()/180)*dAz*TMath::Pi()/180; //Mean of the Monte carlo simulation 
                */

                if ( k == 1)
                {
                    TrueEvents << eta << ", " << e << ", "<< N_ij << "\n";
                    //std::cout << eta << ", " << e << ", "<< N_ij << "\n" << std::endl;
                   hTrue_means->SetBinContent(i,j, N_ij); //Create histogram for  kth Pseudo-Experiment
                }
                
                //Inverse Tranfrom Sampling.

                //double sam_ij=Sample_Pois(N_ij,u);

                std::poisson_distribution<> d(N_ij);
                
                n_ijk = d(gen); //Random sample from a Poisson distribution with mean N_ij

                total += n_ijk;

                hTrue_k->SetBinContent(i,j, n_ijk); //Create histogram for  kth Pseudo-Experiment

                

                int gbin = hTrue->GetBin(i,j); //Generalized mean



                n_ij = hTrue->GetBinContent(gbin)+ n_ijk; //Total Events histogram

    

               // std::cout << n_ijk << " " << n_ij << " " << " " << N_ij << std::endl;
                

                hTrue->SetBinContent(i,j, n_ij);      //Expected Events



                //std::cout << n_ij << " " << nij<< " " << N_ij << std::endl;            
                
            } // loop energy

        } // Loop eta

        std::cout << " total events kth PseudoExp: " << total << std::endl;

        PseudoExp->Add(hTrue_k);

    } // # Of PseudoExperiments.

    PseudoExp -> Add(hTrue);
    PseudoExp -> Add(hTrue_means);

    TrueEvents.close();

    delete HF;
            
             
    return PseudoExp;
}

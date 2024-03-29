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

# define R_earth 6368.0 //km


// Some Constants i need
//# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
//# define MTon  1E9  //Metric MegaTon
//# define years2sec 3.154E7 // Years in Seconds

// Make oscillogram for given final flavour and MH
TH2D*  AsimovTrueEvents(std::string modelname, int flvf, double Region[], int Bins[], double NnT)
{  
    std::cout << "Generating Asimov data set: Reimann Integration Method" << std::endl;

    double total = 0;

    //Data Storage -----------------------------------------------------------------------------------------------------
    std::string model = "/home/dehy0499/OscProb/PremTables/"+modelname+".txt";
    std::string model_default = "/home/dehy0499/OscProb/PremTables/prem_default.txt";
    
    std::cout << modelname <<std::endl;

    std::string location = "SimulationResults/AsimovData/" ;
    
    std::string Earthmodel= modelname;
    
    std::string title = "_Asimov_"+std::to_string(flvf)+"_"+std::to_string(Bins[0])+"_"+std::to_string(Bins[1])+"_"+std::to_string(Region[0])+"_"+std::to_string(Region[1])+".csv";
    
    std::string filename = location+Earthmodel+title;
    
    ofstream TrueEvents(filename, std::ofstream::trunc); //Opens a file and rewrite content, if files does not exist it Creates new file

    
    double N_ij = 0   ;   //Poisson mean for bin ij.
  

    //Binnig scheme and Oscillogram-------------------------------------------------------------------------------------

    //Energy Intervals
    double Emin      = Region[0];//Lower limit for Energy
    double Emax      = Region[1];//Upper limit for Energy

    double etamin   = ( 180-Region[3] )*TMath::Pi()/180;
    double etamax   = ( 180-Region[2] )*TMath::Pi()/180;;

    //Zenith Angle Interval
    //double etamax   = ( 180-Region[2] )*TMath::Pi()/180;//Lower limit for angle
    //double etamin   = ( 180-Region[3])*TMath::Pi()/180;//Upperlimit for angle
    //double cetamin   = cos((180-etamax)*TMath::Pi()/180); 
    //double cetamax   = cos((180-etamin)*TMath::Pi()/180);

    //Phi/Azimuthal Intervals
    double dAz = Region[4]*TMath::Pi()/180;

    //Bins
    int ibins = Bins[0]; // Number of  angular bins of True event distribution
    int jbins = Bins[1]; // Number of  energy bins of True event distribution
    double deta = (etamax - etamin)/(ibins); //< Bin width/2
    double dE = (Emax - Emin)/(jbins); //< Bin width/2

    /* Create 2D histogram for Event Oscillogram,
     xbins correspond to energy values and ybins to zenith angle cosEta*/
    TH2D* hTrue = new TH2D("hTrue","True Events",ibins,etamin,etamax,jbins,Emin,Emax); 
    //TH2D* hTrue_k= new TH2D("hTrue_k","True Events; #eta ; #E",ibins,etamin,etamax,jbins,Emin,Emax); //Store data for each pseudo-experiment.
    //TH2D* hTrue_means= new TH2D("hTrue_means","True Events( Means ); #eta ; #E",ibins,etamin,etamax,jbins,Emin,Emax); //Store data for each pseudo-experiment.
    
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
    OscProb::PremModel prem_model(model);
    OscProb::PremModel prem(model_default); //Default PREM table


    double R_lim = 3500.0; // Distance from the center of the Earth , 3500 Km CMB
       
    double EtaLim = TMath::Pi() - TMath::ASin( (R_lim)/R_earth ) ;
    //double Etamax_LLSVP = TMath::ASin( (R_cmb + h_llsvp)/R_earth )*(180.0/TMath::Pi()) ;




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



        for(int i=1; i<= ibins ; i++) //Loop in Angular bins
        {    

            // Get cos(eta) from bin center, This is used to calculate the baseline.

            double eta = hTrue->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
            
            double cosEta =cos(eta);

            if(cosEta < -1 || cosEta > 1) break; // Skip if cosEta is unphysical 

             

            //double L_Ho = prem.GetTotalL(cosEta); // Set total path length L

            //Select Earth model dependeing on the neutrino path


            if (eta <= EtaLim ) 
             {

                std::cout << eta << " " << EtaLim << " --- " << std::endl;
                double cosEta =cos(eta);
                 
                prem_model.FillPath(cosEta); // Fill paths from PREM model

                PMNS_H.SetPath(prem_model.GetNuPath()); // Set paths in OscProb  

              }
              else

              {

                std::cout << eta << " " << EtaLim << " *** " << std::endl;

                prem.FillPath(cosEta); // Fill paths from PREM model
            
                PMNS_H.SetPath(prem.GetNuPath()); // Set paths in OscProb  
              }

              std::cout << " GET PREM LAYERS- - - - - - - - - - - - - - - - - - - " << std::endl; 

              vector<OscProb::NuPath> p;

              p = prem.GetNuPath();

              int nsteps = p.size();
              
              double totL = prem.GetTotalL(cosEta);

              double sumL = prem.GetTotalL(cosEta);

              for (int i = 0; i < nsteps; ++i)
              {
                  double NuL = p[i].length;  
                  
                  double rho = p[i].density;
                  
                  double idx = p[i].layer;

                  std::cout << "*****" << NuL << " " << rho << " " << idx << " " << std::endl;
              }
              
              std::cout << " - - - - - - - - - - - - - - - - - - - " << std::endl; 
            
            for (int j = 1; j <=jbins; ++j)
            { 
                double e = hTrue->GetYaxis()->GetBinCenter(j); //< This will defined a constant L por different values of ct provided Dct is Small
                
                //Neutrino
                PMNS_H.SetIsNuBar(false); 
                double Ri_nu = XSec(e,nu)*( PMNS_H.Prob(numu, flvf, e)*dPsiMudE.Eval(e) + PMNS_H.Prob(nue,flvf,e)*dPsiEdE.Eval(e) );
                
                //Antineutrino contribution
                PMNS_H.SetIsNuBar(true); 
                double Ri_nubar = XSec(e,nubar)*( PMNS_H.Prob(numu,flvf, e)*dPsiMubardE.Eval(e) + PMNS_H.Prob(nue,flvf,e)*dPsiEbardE.Eval(e) ); 

                double N_ij = NnT*(Ri_nu + Ri_nubar)*dE*deta*dAz;


                TrueEvents << eta << ", " << e << ", "<< N_ij  << "\n";

                hTrue->SetBinContent(i,j, N_ij); //Create histogram for  kth Pseudo-Experimens

            } // loop energy

        } // Loop eta


    TrueEvents.close();

    delete HF;
            
             
    return hTrue;
}

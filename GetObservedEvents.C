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


double w_E(double e_o , double e, double dE_o, double a_E, double E[] ); 

double w_eta(double eta_o, double e, double eta, double deta_o,double a_eta ,double Eta[]);

// Make oscillogram for given final flavour and MH


TH2D*  GetObservedEvents(std::string modelname, int flvf, double Region[], int Bins[],double Det_par[], double NnT)
{   
    //Data Storage -----------------------------------------------------------------------------------------------------
    std::string model = "/home/dehy0499/OscProb/PremTables/"+modelname+".txt";

    std::string location = "SimulationResults/ObsEventsResults/" ;
    std::string Earthmodel= modelname;
    std::string title = "_PO_"+std::to_string(flvf)+"_"+std::to_string(Bins[2])+"_"+std::to_string(Bins[3])+"_"+std::to_string(Region[0])+"_"+std::to_string(Region[1])+".csv";
    std::string filename = location+Earthmodel+title;
    
    ofstream ObservedData(filename, std::ofstream::trunc); //Opens a file and rewrite content, if files does not exist it Creates new file
    

    //Detector Response-------------------------------------------------------------------------------------------------
    
    double a_E =  Det_par[0] ;
    
    double a_eta = Det_par[1] ;


    //Observed Variables:-----------------------------------------------------------------------------------------------
    
    //Observed bins
    int mbins = Bins[2]; //Number of angular bins of Observed event distribution
    
    int nbins = Bins[3]; //Number of Energy bins of Observed event distribution
    
    //Energy[GeV]
    double EOmin      = Region[0];
    
    double EOmax      = Region[1];
    
    double dE_o = (EOmax - EOmin)/(2.0*nbins); //< Bin width/2
    
    double E_itv[] = {EOmin, EOmax};
    
    //Observed Angle
    double etaOmin   = Region[2];
    
    double etaOmax   = Region[3];
    
    double deta_o = (etaOmax - etaOmin)/(2.0*mbins); //< Bin width/2

    double Eta_itv[] = {etaOmin, etaOmax};
    
    //Azimuthal 
    double dAz = Region[4];

    //Observed distribution:
    TH2D* hObs = new TH2D("hObs","Observed Events",mbins,etaOmin,etaOmax,nbins,EOmin,EOmax); // xbins correspond to energy values and ybins to zenith angle cost
    

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = Bins[0]; // Number of  angular bins of True event distribution

    int jbins = Bins[1]; // Number of  energy bins of True event distribution
    
    //True Energy[GeV]
    double Emin=Region[0];

    double Emax=Region[1];
    
    double dE = (Emax - Emin)/(2.0*jbins); //< Bin width/2

    //True Angle
    double etamin= Region[2];

    double etamax= Region[3];

    double deta = (etamax - etamin)/(2.0*ibins); //< Bin width/2

    
    //True distribution        
    TH2D* hTrue = new TH2D("hTrue","True Events",jbins,etamin,etamax,ibins,Emin,Emax); // xbins correspond to energy values and ybins to zenith angle cost





    std::cout << "True variable region ["<< Emin << "-"<<Emax<< "]*["<<etamin<< "-"<<etamax <<"]" << std::endl;

    std::cout << "Observed variable region ["<< EOmin << "-"<<EOmax<< "]*["<<etaOmin<< "-"<<etaOmax <<"]" << std::endl;
    
    //NEUTRINO OSCILLATION PROB-----------------------------------------------------------------------------------------

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

    // Atmospheric Neutrino Flux data-----------------------------------------------------------------------------------
    
    TFile *HF = new TFile("./NuFlux/Honda2014_spl-solmin-allavg.root","read"); //South Pole (IC telescope) Averaged all directions
    
    //TFile *HF = new TFile("./NuFlux/TestFlux.root","read"); // Specific flux V. Agrawal, et al. (1996)
    
    
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



   //Calculation of observed events------------------------------------------------------------------------------------- 
    for (int m = 1; m <= mbins; ++m) //Observed Angular bins
    {

        double eta_o = hObs->GetXaxis()->GetBinCenter(m); // Select angular bin from the "Observed" histogram
       
        for (int n = 1; n <=nbins ; ++n) //Observed Energy bins
        {
     
            double e_o = hObs->GetYaxis()->GetBinCenter(n); // Select Energy bin from the "Observed" histogram
       
            double Ntot = 0;   // Reset total number of events.
            
            for(int i=1; i<= ibins ; i++) //Loop in cosT
            {    

                // Get cos(theta_z) from bin center, It fix a value of L

                double eta = hTrue->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
                
                double cosEta =cos( (180.0-eta)*TMath::Pi()/180.0 );

                if(cosEta < -1 || cosEta > 1) break; // Skip if cosT is unphysical 
                
                double L_Ho = prem.GetTotalL(cosEta); // Set total path length L  
                
                prem.FillPath(cosEta); // Fill paths from PREM model
                
                // Set paths in OscProb  
                PMNS_Ho.SetPath(prem.GetNuPath()); //STANDARD EARTH (DESCRIBED BY PREM MODEL)


                for (int j = 1; j <= jbins; ++j)
                { 
                    //NUMERICAL INTEGRATION

                    double e = hTrue->GetYaxis()->GetBinCenter(j); //< This will defined a constant L por different values of ct provided Dct is Small
                    
                    int nstep = 2; //Steps for Integration of E

                    std::vector<double>  Et, R_Ho; // R= event rate d^2( N )/ (dtheta dE)

                    for (int l = 0; l <= nstep ; ++l)
                    {   
                        Et.push_back( e + (l-1)*dE );

                        //Integration in angular variable using simpson's rule.

                        double etai_min=eta - deta; //Lower bound
                        
                        double etai= eta; //Center
                        
                        double etai_max=eta + deta; //Upper bound
                        
                        double res_eta = (deta/3.0)*( w_eta(eta_o,e,etai_min,deta_o,a_eta,Eta_itv)*sin( (180.0-etai_min)*TMath::Pi()/180.0 )+ 4.0*w_eta(eta_o,e,eta,deta_o,a_eta,Eta_itv)*sin( (180.0-eta)*TMath::Pi()/180.0 ) +  w_eta(eta_o,e,etai_max,deta_o,a_eta,Eta_itv)*sin( (180.0-etai_max)*TMath::Pi()/180.0 ) );

                        double res_e = w_E(e_o, e, dE_o, a_E, E_itv);

                        //double res = w_E( aE, Et[l],E_o,E_GeV)*(dth/3.0)*(  w_th( aTh, Th[0],Et[l],Eta_o,Eta)*sin( (180.0-Th[0])*TMath::Pi()/180.0 )+ 4.0*w_th( aTh, Th[1],Et[l],Eta_o,Eta)*sin( (180.0-Th[1])*TMath::Pi()/180.0 ) +  w_th( aTh, Th[2], Et[l], Eta_o, Eta)*sin( (180.0-Th[2])*TMath::Pi()/180.0 )  );
                        
                        //Calculate the Interaction Rate at bin center E for neutrinos passing trough lower mantle of Standar Earth

                        //Neutrino contribution
                        PMNS_Ho.SetIsNuBar(false); 
                        double r_ho = XSec(Et[l],nu)*( PMNS_Ho.Prob(numu, flvf, Et[l])*dPsiMudE.Eval(Et[l]) + PMNS_Ho.Prob(nue,flvf,Et[l])*dPsiEdE.Eval(Et[l]) );
                        
                        //Antineutrino contribution
                        PMNS_Ho.SetIsNuBar(true); 
                        double r_hobar = XSec(Et[l],nubar)*( PMNS_Ho.Prob(numu,flvf, Et[l])*dPsiMubardE.Eval(Et[l]) + PMNS_Ho.Prob(nue,flvf,Et[l])*dPsiEbardE.Eval(Et[l]) ); 
                        

                        R_Ho.push_back( res_eta*res_e*(r_ho + r_hobar) );


                    }
                    
            
                    //Integration in negery variables using NC quadrature
                    double dNHo_dAzdNnT = simpson(Et, R_Ho); //Integration of Event Rate Interaction for neutrinos passing Homogeneus mantle over E

                    double N_ij=(NnT*(dAz)*(TMath::Pi()/180.0) )*dNHo_dAzdNnT;

                    Ntot += N_ij;
                    
                } // loop e

            } // Loop ct
            
            double No_mn = Ntot;

            ObservedData << eta_o << ", " << e_o << ", "<< No_mn << "\n"; //Write im file.
            
            hObs->SetBinContent(m,n,No_mn); // Set events in Observed histogram. 

         
            

        } // loop in eo

    } // loop in eta_o


   ObservedData.close();
       
   hObs->SetTitle("Observed events for #nu_{#mu} and #bar{#nu}_#mu ");
    
   hObs->SetStats(0);

   //Delete 
   delete hTrue;
   delete HF;
   

   return hObs;
}


//double w_E(double a_E, double E , std::vector<double> eo, double E_GeV[] ) 

double w_E(double e_o , double e, double dE_o, double a_E, double E[] ) 
{

//True

double e_min = E[0];

double e_max = E[1]; 

//Observed

double e_omin = e_o - dE_o;

double e_omax = e_o + dE_o; 

//Resolution function 

double sdE = a_E*e;

double AE = (0.5)*( ROOT::Math::erf (  (e - e_min )/( sqrt(2)*sdE )  ) - ROOT::Math::erf (  (e - e_max )/( sqrt(2)*sdE )  ) );  

double wE = (0.5)*( ROOT::Math::erf (  (e - e_omin )/( sqrt(2)*sdE )  ) - ROOT::Math::erf (  (e - e_omax)/( sqrt(2)*sdE )  ) );    

return wE;

}

//double w_eta(double a_eta, double eta, double E, std::vector<double> eta_o, double Eta[])

double w_eta(double eta_o, double e, double eta, double deta_o,double a_eta ,double Eta[])
{


//True Variable

double eta_min = Eta[0];

double eta_max = Eta[1];

//Observed variable

double eta_omin = eta_o - deta_o;

double eta_omax = eta_o + deta_o;

//Resolution function

double s_eta = a_eta/(sqrt(e));

double A_eta = (180.0/TMath::Pi())*(0.5)*( ROOT::Math::erf (  (TMath::Pi()/180.0)*(eta - eta_min)/( sqrt(2)*s_eta )  ) - ROOT::Math::erf (  (TMath::Pi()/180.0)*(eta - eta_max)/( sqrt(2)*s_eta )  ) );

double w_eta = (TMath::Pi()/180)*(0.5)*( ROOT::Math::erf (  (TMath::Pi()/180.0)*(eta - eta_omin)/( sqrt(2)*s_eta )  ) - ROOT::Math::erf (  (TMath::Pi()/180.0)*(eta - eta_omax)/( sqrt(2)*s_eta )  ) );   

return w_eta;
    
}
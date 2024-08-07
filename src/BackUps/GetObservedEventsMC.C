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
#include "Math/PdfFuncMathCore.h"

//#include "TGraphTrueD.h"

//OSCPROB
//#ifndef __CINT__
//#include "PremModel.h"
//#include "PMNS_Fast.h"
//bool isCINT = false;
//#else
//bool isCINT = true;
//#endif

// My Math Tools
#include "MathTools.h"
#include "PhyTools.h"

using namespace std;


// Some Constants i need
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds

//Resolution functions

double PDF_energy(double e_o, double e, double E_itv[], double a_E);

double PDF_angle(double eta_o,double e, double eta, double Eta_itv[], double a_eta);

// Make oscillogram for given final flavour and MH

//TH2D*  GetTrueEvents(std::string modelname, int flvf, double Region[], int Bins[], double NnT, int K)

TObjArray*  SmearEvents( TObjArray* TrueEvents, int flvf, double Region[], int Bins[],double Det_par[], double NnT, int K )
{   

    double total = 0;

    //Binning Scheme

     //Energy Intervals
    double Emin      = Region[0];//Lower limit for Energy

    double Emax      = Region[1];//Upper limit for Energy

    double E_itv [] = {Emin, Emax}; //For normalizing porpuse 

    //Zenith Angle Interval
    double etamin   = Region[2];//Lower limit for angle

    double etamax   = Region[3];//Upperlimit for angle

    double Eta_itv [] = {etamin, etamax}; //For normalizing porpuse 



    double cetamin   = cos((180-etamax)*TMath::Pi()/180); 

    double cetamax   = cos((180-etamin)*TMath::Pi()/180);

    // True Bins:
    int ibins = Bins[0]; // Number of  angular bins of True event distribution
    
    int jbins = Bins[1]; // Number of  energy bins of True event distribution
    
    double deta = (etamax - etamin)/(ibins); //< Bin width/2

    double dE = (Emax - Emin)/(jbins); //< Bin width/2

    
    // Observed Bins:
    
    int mbins = Bins[2]; //Number of angular bins of Observed event distribution
    
    int nbins = Bins[3]; //Number of Energy bins of Observed event distribution
    
    double deta_o = (etamax - etamin)/(mbins); //< Bin width/2
    
    double dE_o = (Emax - Emin)/(nbins); //< Bin width/2

    
    //Detector resolution.
    
    double a_eta = Det_par[0]; 
    
    double a_E = Det_par[1];

    



    //Observed distribution 
    
    TH2D* hObs_k = new TH2D("hObs","Observed Events: kth experiment",mbins,etamin,etamax,nbins,Emin,Emax); // xbins correspond to energy values and ybins to zenith angle cost
    
    TH2D* hObs = new TH2D("hObs","Observed Events",mbins,etamin,etamax,nbins,Emin,Emax); // xbins correspond to energy values and ybins to zenith angle cost

    TH2D* hObs_means = new TH2D("hObs","Observed Events: Convoluted poisson means",mbins,etamin,etamax,nbins,Emin,Emax); // xbins correspond to energy values and ybins to zenith angle cost

    TObjArray* PseudoExp_O = new TObjArray(K+2);

    TH2D* hTrue;

    double e, e_o, eta, eta_o, Nmn_k, Nmn, Ntrue;

    int tbin , obin; //Generalized bins


    for (int k = 0; k < K; ++k)
    {
        hTrue=(TH2D*)TrueEvents->At(k); //Select the true event histogram generated in the kth pseudo experiment.
    
        total = 0;

       // Observed Events 
        for (int m = 1; m <= mbins; ++m) //Observed Angular bins
        {

             eta_o = hObs->GetXaxis()->GetBinCenter(m); // Select the center of a observed angular bin 
       

            for (int n = 1; n <=nbins ; ++n) //Observed Energy bins
            {
         
                 e_o = hObs->GetYaxis()->GetBinCenter(n); // Select the center of a energy angular bin

                
                for(int i=1; i<= ibins ; i++) //Loop in cosT
                {    

                     eta = hTrue->GetXaxis()->GetBinCenter(i); // Select the center of original/true angular bin 
                    

                    for (int j = 1; j <= jbins; ++j)
                    { 
                       
                        e = hTrue->GetYaxis()->GetBinCenter(j); // Select the center of original/true energy bin

                        tbin = hTrue->GetBin(i,j); //Generalized bin number for bin (i,j) in the True Event distribution.

                        Ntrue = hTrue->GetBinContent(tbin); //Total Events histogram

                        Nmn_k  += PDF_angle(eta_o, e, eta, Eta_itv, a_eta)*PDF_energy(e_o, e, E_itv, a_E)*(Ntrue)*deta*dE*deta_o*dE_o;

                        //Convolution of events
                        
                    } // loop e

                } // Loop ct

                hObs_k->SetBinContent(m,n,Nmn_k);

                obin = hObs->GetBin(m,n); //Generalized bin number for bin (i,j) in the True Event distribution.

                Nmn = hObs->GetBinContent(obin) + Nmn_k; 

                hObs->SetBinContent(m,n, Nmn);   

                total += Nmn;         

            } // loop in e_o

        } // loop in eta_o

        PseudoExp_O->Add(hObs_k);// Smearing of pseudo-experiment.

        std::cout << " total events kth PseudoExp (obs): " << total << std::endl;

    }

    PseudoExp_O->Add(hObs); //Monte Carlo simulation result.

    return PseudoExp_O;
}



double PDF_energy(double e_o, double e, double E_itv[], double a_E)
{
 
//double e_min = e - dE/2;

//double e_max = e + dE/2;

double e_min = E_itv[0];

double e_max = E_itv[1];

double s_e = a_E*e;

double A_E = (0.5)*( ROOT::Math::erf (  (e - e_min )/( sqrt(2)*s_e )  ) - ROOT::Math::erf (  (e - e_max )/( sqrt(2)*s_e )  ) );  

//double f_E = (0.5)*( ROOT::Math::erf (  (E - eo[0] )/( sqrt(2)*sdE )  ) - ROOT::Math::erf (  (E - eo[2])/( sqrt(2)*sdE )  ) );  

double f_E = ROOT::Math::gaussian_pdf (e_o, s_e, e);    

return f_E/A_E;
//return f_E;

}

double PDF_angle(double eta_o,double e,double eta,double Eta_itv[],double a_eta)
{

//double eta_min = eta - deta/2;

//double eta_max = eta + deta/2;

double eta_min = Eta_itv[0];

double eta_max = Eta_itv[1];

double s_eta = a_eta/(sqrt(e));

double A_eta = (180.0/TMath::Pi())*(0.5)*( ROOT::Math::erf (  (TMath::Pi()/180.0)*(eta - eta_min)/( sqrt(2)*s_eta )  ) - ROOT::Math::erf (  (TMath::Pi()/180.0)*(eta - eta_max)/( sqrt(2)*s_eta )  ) ); //Normalization Constant

//double f_eta = (180.0/TMath::Pi())*(0.5)*( ROOT::Math::erf (  (TMath::Pi()/180.0)*(th - eta_o[0])/( sqrt(2)*sdth )  ) - ROOT::Math::erf (  (TMath::Pi()/180.0)*(th - eta_o[2])/( sqrt(2)*sdth )  ) );   

double f_eta = ROOT::Math::gaussian_pdf ( (TMath::Pi()/180.0)*eta_o, s_eta, (TMath::Pi()/180.0)*eta );  

return f_eta/A_eta;
//return f_eta;
    
}

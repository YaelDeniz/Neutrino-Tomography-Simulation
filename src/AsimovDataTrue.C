/*
Notes:

All angular inputs are in degrees

All angular variables are in radias

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
#include "TH3.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include <Math/Interpolator.h>
//#include "TGraphEventsD.h"

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
#include "Earth3DModel.h"

using namespace std;

# define R_earth 6368.0 //km


// Some Constants i need
//# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
//# define MTon  1E9  //Metric MegaTon
//# define years2sec 3.154E7 // Years in Seconds

// Make oscillogram for given final flavour and MH

/*
class AsimovSimulation()
{
    Public:

    //Earth Settings
    std::string PremModel;
    bool MantleAnomaly;
    std::vector<int> AnomalousLayers;

    //Neutrino Settings
    int flvf;
    double EnuMin = 1.00; //GeV
    double EnuMax = 10.00; //GeV

    double ZenMin; // Min 90;
    double ZenMax; // Max 180;

    double AziMin=0; // Min 0;
    double AziMax=360; //Max 360;

    //Simulation Settings
    int nbinsZen; //Bins in Zenith
    int nbinsAzi; //Bins in Azimuth
    int nbinsE; //Bins in Energy

    double NnT; //Detector Exposure

    TH3D * GetTrueEvents3D();



};

*/

TH3D* AsimovSimulation::GetTrueEvents3D()
{  

//Saving information-----------------------------------------------------------------------------------------------------------
    
    std::cout << "Simulation of True events assuming Asimov data set" << std::endl;

    std::string PremFile = PremModel+".txt";

    std::cout << PremModel <<std::endl;

    std::string location = "SimulationResults/AsimovData/" ;
       
    std::string nudetails = "nu"+std::to_string(flvf)+"E"+std::to_string(EnuMin)+std::to_string(EnuMax);
   
    std::string simdetails = "asmvtrue"+std::to_string(nbinsZen)+std::to_string(nbinsAzi)+std::to_string(nbinsE);
   
    std::string filename = location+PremModel+nudetails+simdetails+".csv";
    
    ofstream TrueEvents(filename, std::ofstream::trunc); //Opens a file and rewrite content, if files does not exist it Creates new file

//Binnig scheme and Oscillogram-------------------------------------------------------------------------------------

    //Energy Intervals
    double Emin      = EnuMin;//Lower limit for Energy
    double Emax      = EnuMax;//Upper limit for Energy

    //Angular Intervals

    //Zenith (given  in degrees -> transformed to radians)
    double thmin = ZenMin*TMath::Pi()/180; //[min pi/2]
    double thmax = ZenMax*TMath::Pi()/180;  //[max pi]

    double cthmin = cos(thmax); //min cth = -1 
    double cthmax = cos(thmin); // max cth = 0

    //Azimuth (given  in degrees -> transformed to radians)
    double phimin = AziMin*TMath::Pi()/180; //[min 0]
    double phimax = AziMax*TMath::Pi()/180;  //[max 2pi]

    //Bins
    int ibins = nbinsZen; //Bins in Zenith
    int jbins = nbinsAzi; //Bins in Azimuth
    int kbins = nbinsE;   //Bins in Energy

    double N_ijk = 0   ;   //Poisson mean for bin ijk.

//Histrogram--------------------------------------------------------------------------------------------------------------------

    // TH3D * TrueHist("TrueHist","True Event Histrogram", ibins,thmin,thmax,jbins,phimin,phimax,kbins,Emin,Emax) //binning in th

    TH3D * TrueHist = new TH3D("TrueHist","True Event Histrogram", ibins,cthmin,cthmax,jbins,phimin,phimax,kbins,Emin,Emax); //binning in cth 
    
//Neutrino event generation-----------------------------------------------------------------------------------------

    // Neutrino final flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino

    // Particle type
    int nu         = 1; //neutrino
    int nubar      = -1; //antineutrino

//Neutrino Oscillation Probabilities calculation--------------------------------------------------------------------
    OscProb::PMNS_Fast PMNS_H; // Create PMNS objects
    

    /*
    // Get parameters to PDG
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

//Honda flux distribution-----------------------------------------------------------------------------------------------------

    NuFlux SPflux;

    std::vector< std::vector<double> > FluxData = SPflux.SetFluxData();

    //Matrix for Histogram & Histogram Draw

    TH2D* muflux =  SPflux.GetFluxHist(1,FluxData); //MuFlux
    TH2D* mubflux =  SPflux.GetFluxHist(2,FluxData); //MuBarFlux
    TH2D* eflux =  SPflux.GetFluxHist(3,FluxData); //EFlux
    TH2D* ebflux =  SPflux.GetFluxHist(4,FluxData); //EBarFlux

//Set earth model -------------------------------------------------------------------------------------------------------------
    
     Earth3DModel MyEarthModel;

     MyEarthModel.SetModel(PremFile);

     MyEarthModel.ActiveHeterogeneity( MantleAnomaly );

     MyEarthModel.WhichLayersLLVPs = AnomalousLayers;

//Event Calculation
    
    double l,d,z,ly;

    for (int j = 1; j <= jbins; j++) //Loop In Azimuth
    {
        double phi = TrueHist->GetYaxis()->GetBinCenter(j); //< This will defined a constant L por different values of ct provided Dct is Small
        
        double dphi = (TrueHist -> GetYaxis()->GetBinWidth(j))*(TMath::Pi()/180); //In Radians
        
      
        for(int i=1; i<= ibins ; i++) //Loop in Zenith
        {    

            // Get cos(eta) from bin center, This is used to calculate the baseline.

            double cth = TrueHist->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
            double dcth = (TrueHist -> GetXaxis()->GetBinWidth(i))*(TMath::Pi()/180);
            double th =acos(cth);

            if(cth < -1 || cth > 1) break; // Skip if cosEta is unphysical 
            
            //std::vector< std::vector<double> >  PathMatrix = NuPATHS3D (PREM_MODELTEST, eta , 0.0, LLVP);

            MyEarthModel.SetDirection(th,phi); 

            std::vector<std::vector<double>> EarthPath = MyEarthModel.Create3DPath();

            l = EarthPath[0][0];
            d = EarthPath[0][1];
            z = EarthPath[0][2];
            //ly =  PathMatrix[0][3]; 
            
            PMNS_H.SetPath(l,d,z);
            
           
            for (int i = 1; i < EarthPath.size(); i++) 
            { 
        
                l = EarthPath[i][0];
                d = EarthPath[i][1];
                z = EarthPath[i][2];
               // ly =  PathMatrix[i][3]; 
                
                PMNS_H.AddPath(l,d,z);
            
            } 
            
            for (int k = 1; k <=kbins; ++k)
            { 
                double e = TrueHist->GetZaxis()->GetBinCenter(k); //< This will defined a constant L por different values of ct provided Dct is Small
                double dE = TrueHist->GetZaxis()->GetBinWidth(k);
                double logEi = log10(e);
                

                //Neutrino Flux interpolation

                 double logdPsiMu = muflux->Interpolate(logEi,cth);
                 double logdPsiMub = mubflux->Interpolate(logEi,cth);
                 double logdPsiE = eflux->Interpolate(logEi,cth);
                 double logdPsiEb = ebflux->Interpolate(logEi,cth);

         

                 double dPsiMudEdct = pow(10,logdPsiMu);     //Muon neutrino flux
                 double dPsiMubardEdct = pow(10,logdPsiMub); //Muon anti-neutrino flux
                 double dPsiEdEdct = pow(10,logdPsiE);        //Electron neutrino flux
                 double dPsiEbardEdct = pow(10,logdPsiEb);    //Electron anti-neutrino flux



                //Neutrino Contribution;

                PMNS_H.SetIsNuBar(false); 

                double Ri_e = XSec(e,nu)*( PMNS_H.Prob(nue,flvf,e)*dPsiEdEdct); //Electron neutrino contribution

                double  Ri_mu = XSec(e,nu)*(PMNS_H.Prob(numu, flvf, e)*dPsiMudEdct); //Muon neutrino contribution  
                
                double Ri_nu = Ri_e + Ri_mu;

                //double Ri_nu = XSec(e,nu)*( PMNS_H.Prob(numu, flvf, e)*dPsiMudE.Eval(e) + PMNS_H.Prob(nue,flvf,e)*dPsiEdE.Eval(e) );
                
                //Antineutrino contribution
                PMNS_H.SetIsNuBar(true); 

                double Ri_eb=XSec(e,nubar)*(PMNS_H.Prob(nue,flvf,e)*dPsiEbardEdct ); //Electron anti-neutrino contribution

                double Ri_mub=XSec(e,nubar)*( PMNS_H.Prob(numu,flvf, e)*dPsiMubardEdct ); //Muon anti-neutrino contribution

                double Ri_nubar = Ri_eb + Ri_mub;

                //double Ri_nubar = XSec(e,nubar)*( PMNS_H.Prob(numu,flvf, e)*dPsiMubardE.Eval(e) + PMNS_H.Prob(nue,flvf,e)*dPsiEbardE.Eval(e) ); 

                //Events at bin
                double N_ijk = NnT*(Ri_nu + Ri_nubar)*dE*dcth*dphi;
                //double N_ijk = NnT*(Ri_nu + Ri_nubar)*dE*dth*dphi;


                TrueEvents << cth << "," << phi <<" , " << e << ", "<< N_ijk  << "\n";

                TrueHist->SetBinContent(i,j,k, N_ijk); //Create histogram for  kth Pseudo-Experimens

            } // loop energy

        } // Loop zenith

    } //Loop in Azimuth

    TrueEvents.close();
           
    return TrueHist;
}

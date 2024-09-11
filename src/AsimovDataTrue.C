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
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include <Math/Interpolator.h>
#include "TGraph.h"
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


// Some Constants i need
//# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
//# define MTon  1E9  //Metric MegaTon
//# define years2sec 3.154E7 // Years in Seconds

# define mN   1.67492749804E-27  // Nucleons mass (~ Neutron mass) in kg

// Make oscillogram for given final flavour and MH

std::vector< TH2D* > AsimovSimulation::GetTrueEvents3D()
{  
    std::ofstream IndexAzi("SimulationResults/TrueIndexTable.csv"); 

    // Neutrino final flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino

    // Particle type
    int nu         = 1; //neutrino
    int nubar      = -1; //antineutrino

    //Binning Scheme 
    double cthmin = cos(thmax); // Min value for cth (cth = -1) 
    double cthmax = cos(thmin); // Max value for cth (cth = 1)
    //ibins - Total number of zenith("theta") bins | cos(theta)
    //jbins - Total number of azimuth bins
    //kbins - Total number of Energy bins 
    double N_ijk = 0   ;        // Poisson mean for bin (i,j,k).

    // Histogram definition (3D: cos(theta), azimuth, energy)
    TH3D * EventHist3D = new TH3D("TrueHist","True Event Histrogram", 
                                  ibins,cthmin,cthmax,
                                  jbins,phimin,phimax,
                                  kbins,Emin,Emax); 

    //Neutrino Oscillation Probabilities calculation
    OscProb::PMNS_Fast PMNS_H; // PMNS object
    PMNS_H.SetStdPars();       // Set PDG 3-flavor parameters

    //Load Honda flux data
    NuFlux HondaFlux;
    std::vector< std::vector<double> > FluxData = HondaFlux.SetFluxData(HondaTable);

    //Generate flux histograms for different neutrino flavours
    TH2D* muflux =  HondaFlux.GetFluxHist(1,FluxData); // Muon neutrino
    TH2D* mubflux =  HondaFlux.GetFluxHist(2,FluxData);// Muon antineutrino
    TH2D* eflux =  HondaFlux.GetFluxHist(3,FluxData);  // Electron neutrino
    TH2D* ebflux =  HondaFlux.GetFluxHist(4,FluxData); // Electron antineutrino

    // Set up the Earth's 3D model    
    Earth3DModel Earth3D;
    Earth3D.SetModel(PremTable);
    Earth3D.SetDetector(pos);
    Earth3D.PileThickness = PileHeight;  
    Earth3D.aWidth = aperture;
    Earth3D.SetPile( MantleAnomaly, AnomalyShape, PileDensityContrast, PileChemContrast);
    Earth3D.SetLayerProp(PremTableNumber, DensityContrast, ChemicalContrast);
    
    // Initialize vector for 2D histograms and other variables
    std::vector< TH2D* > HistVec;
    TH2D * EventHist2D[jbins];
    double l,d,z;

    //Loop over azimuth angles
    for (int j = 1; j <= jbins; j++) 
    {
        // Get azimuth center and width for the current bin
        double phi = EventHist3D->GetYaxis()->GetBinCenter(j); 
        double dphi = EventHist3D-> GetYaxis()->GetBinWidth(j); 

        // Create a new 2D histogram for each azimuth bin
        //char histchar[64];
        //strcat(histchar, "truehist%d" , j);
        //const char *histname = histchar;

        TString histname = Form("OscHist%d",j);
        EventHist2D[j] = new TH2D(histname,Form("Oscillogram%d",j), ibins,cthmin,cthmax,kbins,Emin,Emax); 

        // Loop over zenith bins (cos(theta))
        for(int i=1; i<= ibins ; i++) 
        {    
            // Get cos(theta) center and width for the current bin.
            // The distance traveled through the earth "Baseline" is L=-2 * R_earth * cos(theta)

            double cth = EventHist2D[j]->GetXaxis()->GetBinCenter(i);
            double dcth = EventHist2D[j] -> GetXaxis()->GetBinWidth(i);

            if(cth < -1 || cth > 1) break; // Skip unphysical cos(theta)

            // Set Earth model direction and calculate Earth path
            
            Earth3D.SetDirection(cth,phi); 
            std::vector<std::vector<double>> EarthPath = Earth3D.Create3DPath();

            // Set PMNS neutrino oscillation path
            l = EarthPath[0][0];
            d = EarthPath[0][1];
            z = EarthPath[0][2];
            
            PMNS_H.SetPath(l,d,z); // Initialize path
            
            for (int np = 1; np < EarthPath.size(); np++) 
            { 
                l = EarthPath[np][0];
                d = EarthPath[np][1];
                z = EarthPath[np][2];    
                PMNS_H.AddPath(l,d,z);
            } 
            
            // Loop over energy bins
            for (int k = 1; k <=kbins; ++k)
            { 
                double e = EventHist2D[j]->GetYaxis()->GetBinCenter(k); 
                double dE = EventHist2D[j]->GetYaxis()->GetBinWidth(k);
                double logEi = log10(e);
                

                // Interpolate neutrino flux values
                double logdPsiMu = muflux->Interpolate(logEi,cth);
                double logdPsiMub = mubflux->Interpolate(logEi,cth);
                double logdPsiE = eflux->Interpolate(logEi,cth);
                double logdPsiEb = ebflux->Interpolate(logEi,cth);

                // Convert logarithmic flux values to linear
                double dPsiMudEdct = pow(10,logdPsiMu);     // Muon neutrino flux
                double dPsiMubardEdct = pow(10,logdPsiMub); // Muon anti-neutrino flux
                double dPsiEdEdct = pow(10,logdPsiE);       // Electron neutrino flux
                double dPsiEbardEdct = pow(10,logdPsiEb);   // Electron anti-neutrino flux

                PMNS_H.SetIsNuBar(false);// neutrino
                double Ri_e = XSec(e,nu)*( PMNS_H.Prob(nue,flvf,e)*dPsiEdEdct); 
                double  Ri_mu = XSec(e,nu)*(PMNS_H.Prob(numu, flvf, e)*dPsiMudEdct); 
                double Ri_nu = Ri_e + Ri_mu;
                
                PMNS_H.SetIsNuBar(true); //Antineutrino 
                double Ri_eb=XSec(e,nubar)*(PMNS_H.Prob(nue,flvf,e)*dPsiEbardEdct ); 
                double Ri_mub=XSec(e,nubar)*( PMNS_H.Prob(numu,flvf, e)*dPsiMubardEdct );
                double Ri_nubar = Ri_eb + Ri_mub;

                // Calculate Poisson mean for the current bins (N_ijk)
                double N_ijk = NT*(Ri_nu + Ri_nubar)*dE*dcth*dphi; 

                // Set the bin content in the 2D histogram
                EventHist2D[j]->SetBinContent(i,k, N_ijk); 

            } // End energy loop 

        } // End zenith loop

        // Add the current 2D histogram to the vector

        HistVec.push_back(EventHist2D[j]);

    } //End  Azimuth Loop 


    // Close index file and return the vector of histograms      
    IndexAzi.close();
    return HistVec;
}


TH2D* AsimovSimulation::GetTrueEvents2D( ) //To be Deleted
{   


    std::string filename = "ControlData_"+label+"_.csv";
    
    std::ofstream Control("SimulationResults/Control/ControlData/"+filename); 

    std::ofstream ControlFlux("SimulationResults/Control/ControlFlux/ControlFlux_"+label+"_.csv"); 

    std::ofstream ControlProb("SimulationResults/Control/ControlProb/ControlProb_"+label+"_.csv");

    //Binnig scheme and Oscillogram-------------------------------------------------------------------------------------

    //Energy Intervals
//    double Emin      = EnuMin;//Lower limit for Energy
//    double Emax      = EnuMax;//Upper limit for Energy

    //Angular Intervals

    //Zenith (given  in degrees -> transformed to radians)
    //double thmin = ZenMin*TMath::Pi()/180; //[min pi/2]
    //double thmax = ZenMax*TMath::Pi()/180;  //[max pi]

    double cthmin = cos(thmax); //min cth = -1 
    double cthmax = cos(thmin); // max cth = 0

    //Azimuth (given  in degrees -> transformed to radians)
    //double phimin = AziMin*TMath::Pi()/180; //[min 0]
    //double phimax = AziMax*TMath::Pi()/180;  //[max 2pi]
    
    double dphi = phimax - phimin;

    //Bins
    //int ibins = nbinsZen; //Bins in Zenith
    //int jbins = nbinsAzi; //Bins in Azimuth
    //int kbins = nbinsE;   //Bins in Energy



//Histrogram--------------------------------------------------------------------------------------------------------------------

     TH2D * TrueHist = new TH2D("TrueHist","True Event Histrogram", ibins,cthmin,cthmax,kbins,Emin,Emax); //binning in cth 
    
//Neutrino event generation-----------------------------------------------------------------------------------------

    // Neutrino final flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino

    // Particle type
    int nu         = 1; //neutrino
    int nubar      = -1; //antineutrino

//Neutrino Oscillation Probabilities calculation--------------------------------------------------------------------
    OscProb::PMNS_Fast PMNS_H; // Create PMNS objects

    PMNS_H.SetStdPars(); // Set PDG 3-flavor parameters

//Honda flux distribution-----------------------------------------------------------------------------------------------------

    NuFlux HondaFlux;
    std::vector< std::vector<double> > FluxData = HondaFlux.SetFluxData(HondaTable);

    //Matrix for Histogram & Histogram Draw

    TH2D* muflux =  HondaFlux.GetFluxHist(1,FluxData); //MuFlux
    TH2D* mubflux =  HondaFlux.GetFluxHist(2,FluxData); //MuBarFlux
    TH2D* eflux =  HondaFlux.GetFluxHist(3,FluxData); //EFlux
    TH2D* ebflux =  HondaFlux.GetFluxHist(4,FluxData); //EBarFlux

//Set earth model -------------------------------------------------------------------------------------------------------------
    
     Earth3DModel Earth3D;

     Earth3D.SetModel(PremTable);

     Earth3D.SetDetector(pos);

     Earth3D.PileThickness = PileHeight;  

     Earth3D.aWidth = aperture;

     Earth3D.SetPile( MantleAnomaly, AnomalyShape, PileDensityContrast, PileChemContrast);

     Earth3D.SetLayerProp(PremTableNumber, DensityContrast, ChemicalContrast);
    

     //OscProb::PremModel prem(PremTable);
    
     //Event Calculation
    
        double l,d,z,ly;

        double phi = 0.0; //< This will defined a constant L for different values of ct provided Dct is Small
    
        for(int i=1; i<= ibins ; i++) //Loop in Zenith
        {    

            // Get cos(eta) from bin center, This is used to calculate the baseline.

            double cth  = TrueHist->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
            double dcth = TrueHist -> GetXaxis()->GetBinWidth(i);

            if(cth < -1 || cth > 1) break; // Skip if cosEta is unphysical 
            
             
            Earth3D.SetDirection(cth,phi); 

            std::vector<std::vector<double>> EarthPath = Earth3D.Create3DPath();

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
              //ly =  PathMatrix[i][3]; 
                
              PMNS_H.AddPath(l,d,z);
            
            } 
             
          
            

            //prem.FillPath(cth); // Fill paths from PREM model

            //PMNS_H.SetPath(prem.GetNuPath()); // Set paths in OscProb  
            
            for (int k = 1; k <=kbins; ++k)
            { 
                double e = TrueHist->GetYaxis()->GetBinCenter(k); //< This will defined a constant L por different values of ct provided Dct is Small
                double dE = TrueHist->GetYaxis()->GetBinWidth(k);
                double logEi = log10(e);
                
                //Neutrino Flux interpolation

                 double logdPsiMu = muflux->Interpolate(logEi,cth);
                 double logdPsiMub = mubflux->Interpolate(logEi,cth);
                 double logdPsiE = eflux->Interpolate(logEi,cth);
                 double logdPsiEb = ebflux->Interpolate(logEi,cth);

                //Evaluate differential Flux
                 double dPsiMudEdct = pow(10,logdPsiMu);     //Muon neutrino flux
                 double dPsiMubardEdct = pow(10,logdPsiMub); //Muon anti-neutrino flux
                 double dPsiEdEdct = pow(10,logdPsiE);        //Electron neutrino flux
                 double dPsiEbardEdct = pow(10,logdPsiEb);    //Electron anti-neutrino flux


                //Neutrino Contribution;

                PMNS_H.SetIsNuBar(false); 

                double Pef = PMNS_H.Prob(nue,flvf,e);

                double Pmuf = PMNS_H.Prob(numu,flvf, e);
                 
                double R_ef = (XSec(e,nu)/mN)*( Pef*dPsiEdEdct); //Electron neutrino to final-neutrino contribution

                double  R_muf = (XSec(e,nu)/mN)*(Pmuf*dPsiMudEdct); //Muon neutrino to final-neutrino contribution

                double R_f = R_ef + R_muf; // Events in the final-neutrino channel

                
                //Antineutrino contribution
                PMNS_H.SetIsNuBar(true); 

                double Pefb = PMNS_H.Prob(nue,flvf,e);

                double Pmufb = PMNS_H.Prob(numu,flvf, e);

                double R_efb=(XSec(e,nubar)/mN)*(  Pefb*dPsiEbardEdct ); //Electron anti-neutrino to final-anti-neutrino contribution

                double R_mufb=(XSec(e,nubar)/mN)*( Pmufb *dPsiMubardEdct ); //Muon anti-neutrino to final-anti-neutrino contribution

                double R_fbar = R_efb + R_mufb; //// Events in the final-antineutrino channel

                 
                double R = R_f+R_fbar;
                
                double N_f = NT*(R)*dE*(dcth)*dphi;  //Total events in the final-neutrino channel

                //Dispaly events 

                std::cout << i << " " <<  k << " " << ibins << " " << kbins <<  " | " << R << " " << N_f << " " << e << " " << cth << " " << dcth << " " << dphi << " " << dE << std::endl;

                Control << cth << " " << e << " " << dcth << " " << dE << " " <<  R_ef << " " << R_muf << " " << R_efb << " " << R_mufb << " " << R << " " << N_f << std::endl; 

                ControlProb << cth << " " << e << " " << Pef << " " << Pefb << " " << Pmuf << " " << Pmuf << " " << XSec(e,nu) << " " << XSec(e,nubar) << std::endl;

                ControlFlux << cth << " " << e << " " << dPsiEdEdct << " " << dPsiEbardEdct << " " << dPsiMudEdct << " " << dPsiMubardEdct << std::endl;

                TrueHist->SetBinContent(i,k, N_f ); //Create histogram for  kth Pseudo-Experimens

            } // loop energy

            PMNS_H.ClearPath();

        } // Loop zenith

        //Simulation Summary

    std::cout << "Detector location" << pos[2] << std::endl;
    std::cout << "PREM model: " << PremTable  <<std::endl;
    std::cout << "Flux Table: " << HondaTable  <<std::endl;
    std::cout << "Pile| h:" <<  PileHeight << " aWidth:" << aperture << " Shape: " << AnomalyShape << " rhopct: " << PileDensityContrast << std::endl;
    std::cout <<" Simulation Summary| flvf:" << flvf << " IsNuBar:" <<  (nunubar <= 0) << " E: " << Emin <<  " " <<  Emax << " cth:" << cthmin << " " << cthmax << " th:"<< thmin*180.0/TMath::Pi() << " " << thmax*180.0/TMath::Pi() <<  std::endl;


    Control.close();
    ControlProb.close();
    ControlFlux.close();

    return TrueHist;
}


TGraph * AsimovSimulation::GetOscProb( int flvi, int flvf, bool nunubar, double cth ) //To be Deleted
{  
    
    //Energy Intervals
    //double Emin      = EnuMin;//Lower limit for Energy
    //double Emax      = EnuMax;//Upper limit for Energy
    int n = 1000;   //Bins in Energy
    
//Neutrino Oscillation Probabilities calculation--------------------------------------------------------------------
    OscProb::PMNS_Fast PMNS_H; // Create PMNS objects

    PMNS_H.SetStdPars(); // Set PDG 3-flavor parameters


///----------------

    NuFlux HondaFlux;
    std::vector< std::vector<double> > FluxData = HondaFlux.SetFluxData(HondaTable);

//Set earth model -------------------------------------------------------------------------------------------------------------
    
        

     std::cout << "3D" << std::endl;
    /*
     Earth3DModel Earth3D;
     Earth3D.SetModel(PremTable);
     Earth3D.SetDetector(pos);
     Earth3D.SetPile( MantleAnomaly, AnomalyShape, PileDensityContrast, PileChemContrast);
     Earth3D.SetLayerProp(PremTableNumber, DensityContrast, ChemicalContrast);
    
    double l,d,z,ly;
    double phi = 0.0; //< This will defined a constant L for different values of ct provided Dct is Small    
    Earth3D.SetDirection(cth,phi); 
    std::vector<std::vector<double>> EarthPath = Earth3D.Create3DPath();

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
        PMNS_H.AddPath(l,d,z);
        
    } 
    
    */
    
    std::cout << "OscProb" << std::endl;
    OscProb::PremModel prem(PremTable);
    prem.FillPath(cth); // Fill paths from PREM model
    prem.FillPath(cth); // Fill paths from PREM model
    PMNS_H.SetPath(prem.GetNuPath()); // Set paths in OscProb  
    
    PMNS_H.SetIsNuBar(nunubar);
    
    //Oscillation Probabilities

    double dE = (Emax-Emin)/1000.0;

    //TH1D * EHist = new TH1D("EHist", "Energy Value", 1000 , Emin, Emax);

    std::cout <<  Emin << " " << Emax <<  std::endl; 

    TGraph * g = new TGraph(1001);

    for (int i = 0; i < 1001; ++i)
        {
            //double e=EHist->GetXaxis()->GetBinCenter(i);
            double E = Emin + i*dE;
            double pab = PMNS_H.Prob(flvi,flvf,E);
            g->SetPoint(i,E,pab);

          //  std::cout << E << " " << pab << std::endl;
            
        }    

   
    return g;
}



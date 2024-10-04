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
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds


// Make oscillogram for given final flavour and MH

//Simulation with binning in cos(theta)
std::vector< TH2D* > AsimovSimulation::GetTrueEvents3D()
{  
    std::string NuTomoPath= "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation";
    std::string IntFolder = "/SimulationResults/PreliminaryResults/IntEvents/";
    std::string AziIndex = NuTomoPath+IntFolder+"TrueIndexTable.txt";
    std::ofstream IndexAzi(AziIndex); 

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

        IndexAzi << j << " " << phi*180/TMath::Pi() << " " << dphi*180/TMath::Pi() << std::endl; 

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

                double  mN = 1.67492749804E-27; // Approximate mass of nucleon
                // Calculate Poisson mean for the current bins (N_ijk)
                // Note : MT/mN  refers to target nucleons per year
                double N_ijk = (MT/mN)*(Ri_nu + Ri_nubar)*dE*dcth*dphi; 

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

//Simulation with binning in theta
std::vector< TH2D* > AsimovSimulation::GetTrueEvents3Dth()
{  
    std::string NuTomoPath= "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation";
    std::string IntFolder = "/SimulationResults/PreliminaryResults/IntEvents/";
    std::string AziIndex = NuTomoPath+IntFolder+"TrueIndexTable.txt";
    std::ofstream IndexAzi(AziIndex); 

    // Neutrino final flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino

    // Particle type
    int nu         = 1; //neutrino
    int nubar      = -1; //antineutrino

    //Binning Scheme 

    //double cthmin = cos(thmax); // Min value for cth (cth = -1) 
    //double cthmax = cos(thmin); // Max value for cth (cth = 1)
    //ibins - Total number of zenith("theta") bins | cos(theta)
    //jbins - Total number of azimuth bins
    //kbins - Total number of Energy bins 
    double N_ijk = 0   ;        // Poisson mean for bin (i,j,k).

    // Histogram definition (3D: cos(theta), azimuth, energy)
    TH3D * EventHist3D = new TH3D("TrueHist","True Event Histrogram", 
                                  ibins,thmin,thmax,
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

        IndexAzi << j << " " << phi*180/TMath::Pi() << " " << dphi*180/TMath::Pi() << std::endl; 

        // Create a new 2D histogram for each azimuth bin
        //char histchar[64];
        //strcat(histchar, "truehist%d" , j);
        //const char *histname = histchar;

        TString histname = Form("OscHist%d",j);
        EventHist2D[j] = new TH2D(histname,Form("Oscillogram%d",j), ibins,thmin,thmax,kbins,Emin,Emax); 

        // Loop over zenith bins (cos(theta))
        for(int i=1; i<= ibins ; i++) 
        {    
            // Get cos(theta) center and width for the current bin.
            // The distance traveled through the earth "Baseline" is L=-2 * R_earth * cos(theta)

            double th = EventHist2D[j]->GetXaxis()->GetBinCenter(i);
            double cth = cos(th);
            double sth = sin(th);
            double dth = EventHist2D[j] -> GetXaxis()->GetBinWidth(i);

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
                double dPsiMudEdth    = sth*pow(10,logdPsiMu);     // Muon neutrino flux
                double dPsiMubardEdth = sth*pow(10,logdPsiMub); // Muon anti-neutrino flux
                double dPsiEdEdth     = sth*pow(10,logdPsiE);       // Electron neutrino flux
                double dPsiEbardEdth  = sth*pow(10,logdPsiEb);   // Electron anti-neutrino flux

                PMNS_H.SetIsNuBar(false);// neutrino
                double Ri_e = XSec(e,nu)*( PMNS_H.Prob(nue,flvf,e)*dPsiEdEdth); 
                double  Ri_mu = XSec(e,nu)*(PMNS_H.Prob(numu, flvf, e)*dPsiMudEdth); 
                double Ri_nu = Ri_e + Ri_mu;
                
                PMNS_H.SetIsNuBar(true); //Antineutrino 
                double Ri_eb=XSec(e,nubar)*(PMNS_H.Prob(nue,flvf,e)*dPsiEbardEdth ); 
                double Ri_mub=XSec(e,nubar)*( PMNS_H.Prob(numu,flvf, e)*dPsiMubardEdth );
                double Ri_nubar = Ri_eb + Ri_mub;

                double  mN = 1.67492749804E-27; // Approximate mass of nucleon
                // Calculate Poisson mean for the current bins (N_ijk)
                // Note : MT/mN  refers to target nucleons per year
                double N_ijk = (MT/mN)*(Ri_nu + Ri_nubar)*dE*dth*dphi; 

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


std::vector< TH2D* > AsimovSimulation::GetTrueEvents2D( ) //To be Deleted
{   
    //Define cos(theta) range and azimuthal angle range
    double cthmin = cos(thmax); // min cos(theta) = -1 
    double cthmax = cos(thmin); // max cos(theta) = 0
    double dphi = phimax - phimin; //Azimutham angle range

    //Bin labels
    // ibins - Bins in Zenith
    // jbins - Bins in Azimuth
    // kbins - Bins in Energy

    //Initialize histograms for neutrino events
    TH2D * EventHist = new TH2D("EventHist","True Event Histrogram", 
                                 ibins,cthmin,cthmax,kbins,Emin,Emax); 
    TH2D * Evtseflv = new TH2D("Evtseflv","#nu_{e} contribution to True Event", 
                                 ibins,cthmin,cthmax,kbins,Emin,Emax); 
    TH2D * Evtseflvb = new TH2D("Evtseflvb","#bar{#nu}_{e} contribution to True Event", 
                                 ibins,cthmin,cthmax,kbins,Emin,Emax); 
    TH2D * Evtsmuflv = new TH2D("Evtsmuflv","#nu_{#mu} contribution to True Event", 
                                 ibins,cthmin,cthmax,kbins,Emin,Emax);
    TH2D * Evtsmuflvb = new TH2D("Evtsmuflvb","#bar{#nu}_{#mu} contribution to True Event", 
                                 ibins,cthmin,cthmax,kbins,Emin,Emax); 
                            
     std::vector < TH2D * >  Histflv; //Store Oscillograms for (muflv, muflvbar, eflv, eflvbar, flv+ flvbar)
     
     TH2D * Hist[5];
    
    //Neutrino properties

    // Neutrino initial flavour
    int nue        = 0;  // Electron neutrino  
    int numu       = 1;  // Muon neutrino
    // Particle type
    int nu         = 1;  // Neutrino
    int nubar      = -1; // Antineutrino

    // PMNS matrix setup
    OscProb::PMNS_Fast PMNS_H; 
    PMNS_H.SetStdPars(); // Set PDG 3-flavor parameters

    // Load Honda flux data and create histograms for different neutrino types
    NuFlux HondaFlux;
    std::vector< std::vector<double> > FluxData = HondaFlux.SetFluxData(HondaTable);
    //Generate flux histograms for different neutrino flavours
    TH2D* muflux =  HondaFlux.GetFluxHist(1,FluxData); // Muon neutrino
    TH2D* mubflux =  HondaFlux.GetFluxHist(2,FluxData);// Muon antineutrino
    TH2D* eflux =  HondaFlux.GetFluxHist(3,FluxData);  // Electron neutrino
    TH2D* ebflux =  HondaFlux.GetFluxHist(4,FluxData); // Electron antineutrino

    OscProb::PremModel prem(PremTable);


    // Set Earth model and Detector    
    Earth3DModel Earth3D;
    Earth3D.SetModel(PremTable); // Earth model
    Earth3D.SetDetector(pos);    // Detector position
    // Modify Specific layer properties
    Earth3D.SetLayerProp(PremTableNumber, 
                         DensityContrast, ChemicalContrast); 
    // Set heterogeneities in lower mantle (e.g., LLVPs)
    Earth3D.PileThickness = PileHeight;  // Height of LLVP
    Earth3D.aWidth = aperture;           // Aperture of LLVP
    Earth3D.SetPile( MantleAnomaly, 
                   AnomalyShape, PileDensityContrast, PileChemContrast);
        
        // Loop over zenith bins (cos(theta))
        for(int i=1; i<= ibins ; i++) //Loop in Zenith
        {    
            // Get the center and width of a bin in cos(theta)
            double cth  = EventHist->GetXaxis()->GetBinCenter(i);
            double dcth = EventHist -> GetXaxis()->GetBinWidth(i);

            if(cth < -1 || cth > 1) break; // Skip unphysical values 
            
            // Set the path through the Earth for direction (cos(theta),phi)
            prem.FillPath(cth); // Fill paths from PREM model
            prem.FillPath(cth); // Fill paths from PREM model
            PMNS_H.SetPath(prem.GetNuPath()); // Set paths in OscProb  
            //Earth3D.SetDirection(cth,0);  
            //std::vector<std::vector<double>> EarthPath = Earth3D.Create3DPath();

            // Distance through "ith" layer : l - EarthPath[i][0];
            // Densty of "ith" layer : rho - EarthPath[i][2];
            // z/a of "ith" layer : zoa - EarthPath[i][3];
            
            // Set initial Earth path in the PMNS model
            //PMNS_H.SetPath(EarthPath[0][0],EarthPath[0][1],EarthPath[0][2]);
            
            // Add subsequent segments of the path
            //for (int i = 1; i < EarthPath.size(); i++) 
            //{ 
              //PMNS_H.AddPath(EarthPath[i][0],EarthPath[i][1],EarthPath[i][2]);
            //} 

            // Loop through energy bins 
            for (int k = 1; k <=kbins; ++k)
            { 
                double e = EventHist->GetYaxis()->GetBinCenter(k); 
                double dE = EventHist->GetYaxis()->GetBinWidth(k);
                double logEi = log10(e);
                
                //Interpolate neutrino fluxes
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
                PMNS_H.SetIsNuBar(false); // Set neutrino probabilities 
                // iflv -> fflv oscillation probabilities
                double Pef = PMNS_H.Prob(nue,flvf,e);   // nue -> nufinal
                double Pmuf = PMNS_H.Prob(numu,flvf, e);// numu -> nufinal

                // Neutrino Interacting rate
                double R_ef = XSec(e,nu)*( Pef*dPsiEdEdct);
                double  R_muf = XSec(e,nu)*(Pmuf*dPsiMudEdct); 
                double R_f = R_ef + R_muf; // Events in the final-neutrino channel

                //Antineutrino contribution
                PMNS_H.SetIsNuBar(true); // Set antineutrino probabilities 
                // iflv -> fflv oscillation probabilities
                double Pefb = PMNS_H.Prob(nue,flvf,e);   // nuebar -> nufinalbar
                double Pmufb = PMNS_H.Prob(numu,flvf, e);// numubar -> nufinalbar
                
                // Antineutrino Interacting rate
                double R_efb=XSec(e,nubar)*(  Pefb*dPsiEbardEdct ); 
                double R_mufb=XSec(e,nubar)*( Pmufb *dPsiMubardEdct ); 
                double R_fbar = R_efb + R_mufb; //// Events in the final-antineutrino channel

                // Total Event rate
                double R = R_f+R_fbar;

                double  mN = 1.67492749804E-27; // Approximate mass of nucleon
                //double  Meff = 14.6*pow(log10(e),1.8 )*MTon;
                double    Meff = 1*MTon;
                //double  N_A = 6.02214E23;
                double  T = 1*years2sec;
                //double  N_f = (2*TMath::Pi())*(Meff*T/mN)*R*dE*dcth;  // Number of events for bin (cth,e)

                  
                //Note : MT/mN  refers to target nucleons per year
                //double  N_f = (2*TMath::Pi())*(MT/mN)*R*dE*dcth;  // Number of events for bin (cth,e)
                double  N_f = (e*e)*(2*TMath::Pi())*(Meff*T/mN)*R*dcth;  // Number of events for bin (cth,e) 

                Evtseflv->SetBinContent(i,k,  ((Meff*T/mN)*dE*dcth*2*TMath::Pi())*R_ef );
                Evtseflvb->SetBinContent(i,k, ((Meff*T/mN)*dE*dcth*2*TMath::Pi())*R_efb );
                Evtsmuflv->SetBinContent(i,k, ((Meff*T/mN)*dE*dcth*2*TMath::Pi())*R_muf );
                Evtsmuflvb->SetBinContent(i,k,((Meff*T/mN)*dE*dcth*2*TMath::Pi())*R_mufb );
                EventHist->SetBinContent(i,k, N_f ); //Create histogram for  kth Pseudo-Experimens

            } // End of loop in energy

            // Clear PMNS path for the next iteration
            PMNS_H.ClearPath();

        }   //   End of loop in zenith

    //Store all histrogram in a vector
    
    Histflv.push_back(Evtsmuflv); // Muon neutrino contribution to events
    Histflv.push_back(Evtsmuflvb);// Muon antineutrino contribution to events
    Histflv.push_back(Evtseflv);  // Electron neutrino contribution to events
    Histflv.push_back(Evtseflvb); // Electron antineutrino contribution to events
    Histflv.push_back(EventHist); // flvf neutrino events

    // Simulation summary output
    std::cout << "Detector location" << pos[2] << std::endl;
    std::cout << "PREM model: " << PremTable  <<std::endl;
    std::cout << "Flux Table: " << HondaTable  <<std::endl;
    std::cout << "Pile| h:" <<  PileHeight << " aWidth:" << aperture 
              << " Shape: " << AnomalyShape << " rhopct: " << PileDensityContrast << std::endl;
    std::cout <<" Simulation Summary| flvf:" << flvf << " IsNuBar:" <<  (nunubar <= 0) << " E: " 
              << Emin <<  " " <<  Emax << " cth:" << cthmin << " " << cthmax << " th:"<< thmin*180.0/TMath::Pi() 
              << " " << thmax*180.0/TMath::Pi() <<  std::endl;

    return Histflv;
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



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

#define MyPremTables "/home/dehy0499/OscProb/PremTables/"

// Make oscillogram for given final flavour and MH

TH3D* AsimovSimulation::GetTrueEvents3D()
{  

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
    
     Earth3DModel MyEarthModel;

     MyEarthModel.SetModel(PremTable);

     MyEarthModel.SetDetector(Rdet);

     MyEarthModel.SetPile( MantleAnomaly, AnomalyShape, PileDensityContrast, PileChemContrast);

     MyEarthModel.SetLayerProp(PremTableNumber, DensityContrast, ChemicalContrast);

//Event Calculation
    
    double l,d,z,ly;

    for (int j = 1; j <= jbins; j++) //Loop In Azimuth
    {
        double phi = TrueHist->GetYaxis()->GetBinCenter(j); //< This will defined a constant L por different values of ct provided Dct is Small
        
        double dphi = TrueHist -> GetYaxis()->GetBinWidth(j); //In Radians
        
      
        for(int i=1; i<= ibins ; i++) //Loop in Zenith
        {    

            // Get cos(eta) from bin center, This is used to calculate the baseline.

            double cth = TrueHist->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
            double dcth = (TrueHist -> GetXaxis()->GetBinWidth(i));
            //double th =acos(cth);

            if(cth < -1 || cth > 1) break; // Skip if cosEta is unphysical 
            
            //std::vector< std::vector<double> >  PathMatrix = NuPATHS3D (PREM_MODELTEST, eta , 0.0, LLVP);

            MyEarthModel.SetDirection(cth,phi); 

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
                

                //Neutrino Flux bilinear interpolation


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
                //double N_ijk = (Ri_nu + Ri_nubar); //Rates
                double N_ijk = NnT*(Ri_nu + Ri_nubar)*dE*dcth*dphi; //Binning in cth
                //double N_ijk = NnT*(Ri_nu + Ri_nubar)*dE*dth*dphi; //Binning in th

                TrueHist->SetBinContent(i,j,k, N_ijk); //Create histogram for  kth Pseudo-Experimens

            } // loop energy

        } // Loop zenith

    } //Loop in Azimuth

    std::cout << "  3D true simulation summary: " << std::endl;
    std::cout << "  Energy range [" << EnuMin << " - " << EnuMax << "] Binning " << kbins  << std::endl;
    std::cout << "  Zenith range [" << cthmin << "(" << ZenMax <<") - "<< cthmax << "(" << ZenMin <<")]" << " Binning " << ibins  << std::endl;
    std::cout << "  Azimuth range [" << AziMin << " - " << AziMax << "] Binning " << jbins  << std::endl;

          
    return TrueHist;
}


TH2D* AsimovSimulation::GetTrueEvents2D( ) //To be Deleted
{  
    std::cout<< " ** * * * * * * * * * * * Internal Confirmation "<< HondaTable << std::endl;
    
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
    double dphi = phimax - phimin;

    //Bins
    int ibins = nbinsZen; //Bins in Zenith
    int jbins = nbinsAzi; //Bins in Azimuth
    int kbins = nbinsE;   //Bins in Energy

    double N_ik = 0   ;   //Poisson mean for bin ik.

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
    
     Earth3DModel MyEarthModel;

     MyEarthModel.SetModel(PremTable);

     MyEarthModel.SetDetector(Rdet);

     MyEarthModel.PileThickness = PileHeight;  

     MyEarthModel.aWidth = aperture;

     MyEarthModel.SetPile( MantleAnomaly, AnomalyShape, PileDensityContrast, PileChemContrast);

     MyEarthModel.SetLayerProp(PremTableNumber, DensityContrast, ChemicalContrast);
    
     
     //std::string model = "/home/dehy0499/OscProb/PremTables/"+PremFile;

     //OscProb::PremModel prem(model);
    
     //Event Calculation
    
    double l,d,z,ly;

        double phi = 0.0; //< This will defined a constant L for different values of ct provided Dct is Small
    
        for(int i=1; i<= ibins ; i++) //Loop in Zenith
        {    

            // Get cos(eta) from bin center, This is used to calculate the baseline.

            double cth  = TrueHist->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
            double dcth = TrueHist -> GetXaxis()->GetBinWidth(i);

            if(cth < -1 || cth > 1) break; // Skip if cosEta is unphysical 
            
            MyEarthModel.SetDirection(cth,phi); 

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
                //double N_ik = NnT*(Ri_nu + Ri_nubar)*dE*dcth*dphi; //Binning in cos(th)
                //double N_ik = NnT*(Ri_nu + Ri_nubar)*dE*(-sin(th)*dth)*dphi;  //Binning in th
                
                N_ik = NnT*(Ri_nu + Ri_nubar)*dE*(dcth)*dphi;  //Binning in th
            
                TrueHist->SetBinContent(i,k, N_ik); //Create histogram for  kth Pseudo-Experimens

            } // loop energy

        } // Loop zenith
   
    return TrueHist;
}

TH2D* AsimovSimulation::GetOscProb2D( int flvi, int flvf, bool nunubar ) //To be Deleted
{  
    
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

//Histrogram--------------------------------------------------------------------------------------------------------------------

     TH2D * ProbHist = new TH2D("ProbHist","Probabilities of Neutrinos Oscillation inside the Earth", ibins,cthmin,cthmax,kbins,Emin,Emax); //binning in cth 
    
//Neutrino Oscillation Probabilities calculation--------------------------------------------------------------------
    OscProb::PMNS_Fast PMNS_H; // Create PMNS objects

    PMNS_H.SetStdPars(); // Set PDG 3-flavor parameters

//Set earth model -------------------------------------------------------------------------------------------------------------
    
    Earth3DModel MyEarthModel;

    MyEarthModel.SetModel(PremTable);

    MyEarthModel.SetDetector(Rdet);

    MyEarthModel.SetPile( MantleAnomaly, AnomalyShape, PileDensityContrast, PileChemContrast);

    MyEarthModel.SetLayerProp(PremTableNumber, DensityContrast, ChemicalContrast);
    
    OscProb::PremModel prem(PremTable);
    
//Event Calculation
    
    double l,d,z,ly;

    double phi = 0.0; //< This will defined a constant L for different values of ct provided Dct is Small
    
    for(int i=1; i<= ibins ; i++) //Loop in Zenith
    {    

        // Get cos(eta) from bin center, This is used to calculate the baseline.

        double cth  = ProbHist->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
        double dcth = ProbHist -> GetXaxis()->GetBinWidth(i);

        if(cth < -1 || cth > 1) break; // Skip if cosEta is unphysical 

            
        MyEarthModel.SetDirection(cth,phi); 

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
            


        prem.FillPath(cth); // Fill paths from PREM model

        std::vector<OscProb::NuPath> paths = prem.GetNuPath();

   //     PMNS_H.SetPath(prem.GetNuPath()); // Set paths in OscProb  
            
            for (int k = 1; k <=kbins; ++k)
            { 
                double e = ProbHist->GetYaxis()->GetBinCenter(k); //< This will defined a constant L por different values of ct provided Dct is Small

                PMNS_H.SetIsNuBar(nunubar); 

                double  Pab = PMNS_H.Prob(flvi, flvf, e); //Muon neutrino contribution  
                
                ProbHist->SetBinContent(i,k, Pab); //Create histogram for  kth Pseudo-Experimens

            } // loop energy

            std::vector<OscProb::NuPath> CheckPath = PMNS_H.GetPath();

            std::cout << "Check if Paths are being filled correctly" << CheckPath.size() - paths.size() <<  std::endl;

            PMNS_H.ClearPath();

        } // Loop zenith
   
    return ProbHist;
}

TGraph * AsimovSimulation::GetOscProb( int flvi, int flvf, bool nunubar, double cth ) //To be Deleted
{  
    
    //Energy Intervals
    double Emin      = EnuMin;//Lower limit for Energy
    double Emax      = EnuMax;//Upper limit for Energy
    int n = 1000;   //Bins in Energy
    
//Neutrino Oscillation Probabilities calculation--------------------------------------------------------------------
    OscProb::PMNS_Fast PMNS_H; // Create PMNS objects

    PMNS_H.SetStdPars(); // Set PDG 3-flavor parameters


///----------------

    NuFlux HondaFlux;
    std::vector< std::vector<double> > FluxData = HondaFlux.SetFluxData(HondaTable);

//Set earth model -------------------------------------------------------------------------------------------------------------
    
        

     std::cout << "3D" << std::endl;
     Earth3DModel MyEarthModel;
     MyEarthModel.SetModel(PremTable);
     MyEarthModel.SetDetector(Rdet);
     MyEarthModel.SetPile( MantleAnomaly, AnomalyShape, PileDensityContrast, PileChemContrast);
     MyEarthModel.SetLayerProp(PremTableNumber, DensityContrast, ChemicalContrast);
    
    double l,d,z,ly;
    double phi = 0.0; //< This will defined a constant L for different values of ct provided Dct is Small    
    MyEarthModel.SetDirection(cth,phi); 
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
        PMNS_H.AddPath(l,d,z);
        
    } 
    
    
    /*
    std::cout << "OscProb" << std::endl;
    OscProb::PremModel prem(PremTable);
    prem.FillPath(cth); // Fill paths from PREM model
    prem.FillPath(cth); // Fill paths from PREM model
    PMNS_H.SetPath(prem.GetNuPath()); // Set paths in OscProb  
    */
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


std::vector< std::vector<double> > AsimovSimulation::GetPremMatrix( std::string path2table )
{
   

   std::vector< std::vector<double> > PremMatrix; // Matrix data  form of Prem tables

   std::vector<double> PremRow(4);               // PremRow = [radius density z/a, layer ID]

   double radius, density, zoa, layer; // Variables for storing Prem table values

   std::ifstream PREM_DATA;

   PREM_DATA.open(path2table);

      // Loop over table rows
   while(PREM_DATA >> radius >> density >> zoa >> layer)
   {

      PremRow[0] = radius;
      PremRow[1] = density;
      PremRow[2] = zoa;
      PremRow[3] = layer;

      PremMatrix.push_back(PremRow);
    }

   return PremMatrix;
}


TH2D* AsimovSimulation::SensitivityTrueEvents2D( int n, double pct ) //To be Deleted
{  
    int layer = n - 1;
   

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
    double dphi = phimax - phimin;

    //Bins
    int ibins = nbinsZen; //Bins in Zenith
    int jbins = nbinsAzi; //Bins in Azimuth
    int kbins = nbinsE;   //Bins in Energy

    double N_ik = 0   ;   //Poisson mean for bin ik.

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
    
 //    Earth3DModel MyEarthModel;

 //    MyEarthModel.SetModel(PremTable);

 //    MyEarthModel.SetDetector(Rdet);

 //    MyEarthModel.SetPile( MantleAnomaly, AnomalyShape, PileDensityContrast, PileChemContrast);

 //    MyEarthModel.SetLayerProp(PremTableNumber, DensityContrast, ChemicalContrast);
    
     //std::vector< std::vector<double> > PremData = GetPremMatrix( PremFile );

     OscProb::PremModel prem(PremTable);
    
//Event Calculation
    
    double l,d,z,ly;

        double phi = 0.0; //< This will defined a constant L for different values of ct provided Dct is Small
    
        for(int i=1; i<= ibins ; i++) //Loop in Zenith
        {    

            // Get cos(eta) from bin center, This is used to calculate the baseline.

            double cth  = TrueHist->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
            double dcth = TrueHist -> GetXaxis()->GetBinWidth(i);

            if(cth < -1 || cth > 1) break; // Skip if cosEta is unphysical 
          
          /*  
            MyEarthModel.SetDirection(cth,phi); 

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
            
        */
            prem.FillPath(cth); // Fill paths from PREM model

            PMNS_H.SetPath(prem.GetNuPath()); // Set paths in OscProb  
            
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
                //double N_ik = NnT*(Ri_nu + Ri_nubar)*dE*dcth*dphi; //Binning in cos(th)
                //double N_ik = NnT*(Ri_nu + Ri_nubar)*dE*(-sin(th)*dth)*dphi;  //Binning in th
                
                N_ik = NnT*(Ri_nu + Ri_nubar)*dE*(dcth)*dphi;  //Binning in th
            
                TrueHist->SetBinContent(i,k, N_ik); //Create histogram for  kth Pseudo-Experimens

            } // loop energy

        } // Loop zenith
   
    return TrueHist;
}


/*
TH2D* AsimovSimulation::TestTrueEvents2D(std::string model_std , std::string model_alt)
{  

//Saving information-----------------------------------------------------------------------------------------------------------
    
    std::cout << "TEST OF NEUTRINO EVENT CALCULATION" << std::endl;

    std::string PremFile = PremModel+".txt";

    std::cout << PremModel <<std::endl;

    std::string location = "SimulationResults/AsimovData/" ;
       
    std::string nudetails = "nu"+std::to_string(flvf)+"E"+std::to_string(EnuMin)+std::to_string(EnuMax);
   
    std::string simdetails = "asmvtrue2D"+std::to_string(nbinsZen)+std::to_string(nbinsE);
   
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
    double dphi = phimax - phimin;


     //SECTION TO BE DELETED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double R_cmb = 3480.0;
    double thtest = TMath::Pi() - TMath::ASin( (R_cmb)/6368.0 ) ; // max 180


    std::cout << "Region: " << thmin*180.0/TMath::Pi() << " " <<  thmax*180.0/TMath::Pi() << " " << thtest*180.0/TMath::Pi() << " " << Emin << " " << Emax << std::endl;


    //Bins
    int ibins = nbinsZen; //Bins in Zenith
    int jbins = nbinsAzi; //Bins in Azimuth
    int kbins = nbinsE;   //Bins in Energy

    double N_ik = 0   ;   //Poisson mean for bin ijk.

//Histrogram--------------------------------------------------------------------------------------------------------------------

    // TH3D * TrueHist("TrueHist","True Event Histrogram", ibins,thmin,thmax,jbins,phimin,phimax,kbins,Emin,Emax) //binning in th

     //TH2D * TrueHist = new TH2D("TrueHist","True Event Histrogram", ibins,cthmin,cthmax,kbins,Emin,Emax); //binning in cth 
     TH2D * TrueHist = new TH2D("TrueHist","True Event Histrogram", ibins,thmin,thmax,kbins,Emin,Emax); //binning in cth 
    
//Neutrino event generation-----------------------------------------------------------------------------------------

    // Neutrino final flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino

    // Particle type
    int nu         = 1; //neutrino
    int nubar      = -1; //antineutrino

//Neutrino Oscillation Probabilities calculation--------------------------------------------------------------------
    OscProb::PMNS_Fast PMNS_H; // Create PMNS objects
    

    
    // Get parameters to PDG
    //double dm21 = 7.42e-5;
    //double dm31 = 2.533e-3;
    //double th12 = 33.7712*TMath::Pi()/180;
    //double th13 = 8.588*TMath::Pi()/180;
    //double th23 = 48.504*TMath::Pi()/180;
    //double dcp  = 214.2*TMath::Pi()/180;

    // Set PMNS parameters
    //PMNS_H.SetDm(2, dm21);
    //PMNS_H.SetDm(3, dm31);
    //PMNS_H.SetAngle(1,2, th12);
    //PMNS_H.SetAngle(1,3, th13);
    //PMNS_H.SetAngle(2,3, th23);
    //PMNS_H.SetDelta(1,3, dcp);
    
    
    PMNS_H.SetStdPars(); // Set PDG 3-flavor parameters

//Honda flux distribution-----------------------------------------------------------------------------------------------------

    NuFlux HondaFlux;
    HondaFlux.MediterraneanSeaFlux(); // Gran Sasso
    //HondaFlux.SouthPoleFlux(); 
    std::vector< std::vector<double> > FluxData = HondaFlux.SetFluxData(HondaFlux.FluxFileName);

    //Matrix for Histogram & Histogram Draw

    TH2D* muflux =  HondaFlux.GetFluxHist(1,FluxData); //MuFlux
    TH2D* mubflux =  HondaFlux.GetFluxHist(2,FluxData); //MuBarFlux
    TH2D* eflux =  HondaFlux.GetFluxHist(3,FluxData); //EFlux
    TH2D* ebflux =  HondaFlux.GetFluxHist(4,FluxData); //EBarFlux



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









//Set earth model -------------------------------------------------------------------------------------------------------------
    
     Earth3DModel MyEarthModel;

     MyEarthModel.SetModel(PremFile);

     MyEarthModel.SetPile( MantleAnomaly, AnomalyShape, PileDensityContrast, PileChemContrast);

     //MyEarthModel.WhichLayersLLVPs = AnomalousLayers;

     

     double l,d,z,ly;
//Event Calculation
//Previous approach------------------------------------------------------------------------------------------------------------
//Earth model for neutrino propagation 
    
    OscProb::PremModel prem_std(model_std);
    
    OscProb::PremModel prem_alt(model_alt); //Default PREM table

     

        double phi = 0.0; //< This will defined a constant L por different values of ct provided Dct is Small
    
        for(int i=1; i<= ibins ; i++) //Loop in Zenith
        {    

            // Get cos(eta) from bin center, This is used to calculate the baseline.

            //double cth = TrueHist->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
            //double dcth = TrueHist -> GetXaxis()->GetBinWidth(i);
            

            double th  = TrueHist->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
            double cth = cos(th);
            double dth = (TrueHist -> GetXaxis()->GetBinWidth(i))*(TMath::Pi()/180);

            if(cth < -1 || cth > 1) break; // Skip if cosEta is unphysical 
            
            //std::vector< std::vector<double> >  PathMatrix = NuPATHS3D (PREM_MODELTEST, eta , 0.0, LLVP);

            MyEarthModel.SetDirection(cth,phi); 

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
            
            
             if (th > thtest ) 
             {

                std::cout << th << " " << thtest << " --- " << std::endl;

                prem_std.FillPath(cth); // Fill paths from PREM model

                PMNS_H.SetPath(prem_std.GetNuPath()); // Set paths in OscProb  

              }
              else

              {

                std::cout << th << " " << thtest << " *** " << std::endl;

                prem_alt.FillPath(cth); // Fill paths from PREM model

                PMNS_H.SetPath(prem_alt.GetNuPath()); // Set paths in OscProb  
              }

           
            
            for (int k = 1; k <=kbins; ++k)
            { 
                double e = TrueHist->GetYaxis()->GetBinCenter(k); //< This will defined a constant L por different values of ct provided Dct is Small
                double dE = TrueHist->GetYaxis()->GetBinWidth(k);
                double logEi = log10(e);
                

                //Neutrino Flux interpolation

                 std::cout << " ////////////////////////////////////////////////////////////////////Points to interpolate " << Emin << " " << logEi << " " << cth << std::endl; 
                 double logdPsiMu = muflux->Interpolate(logEi,cth);
                 double logdPsiMub = mubflux->Interpolate(logEi,cth);
                 double logdPsiE = eflux->Interpolate(logEi,cth);
                 double logdPsiEb = ebflux->Interpolate(logEi,cth);

                
                 std::cout << "Interpolation of the 2d honda flux "<< std::endl;
                 double dPsiMudEdct = pow(10,logdPsiMu);     //Muon neutrino flux
                 double dPsiMubardEdct = pow(10,logdPsiMub); //Muon anti-neutrino flux
                 double dPsiEdEdct = pow(10,logdPsiE);        //Electron neutrino flux
                 double dPsiEbardEdct = pow(10,logdPsiEb);    //Electron anti-neutrino flux
                 
                 
                 
                 std::cout << "Interpolation of the 1d honda flux "<< std::endl;
                 double dPsiMudEdct = dPsiMudE.Eval(e);     //Muon neutrino flux
                 double dPsiMubardEdct = dPsiMubardE.Eval(e); //Muon anti-neutrino flux
                 double dPsiEdEdct = dPsiEdE.Eval(e);        //Electron neutrino flux
                 double dPsiEbardEdct = dPsiEbardE.Eval(e);    //Electron anti-neutrino flux
                 

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
                //double N_ik = NnT*(Ri_nu + Ri_nubar)*dE*dcth*dphi; //Binning in cos(th)
                //double N_ik = NnT*(Ri_nu + Ri_nubar)*dE*(-sin(th)*dth)*dphi;  //Binning in th
                double N_ik = NnT*(Ri_nu + Ri_nubar)*dE*(dth)*dphi;  //Binning in th

                //double N_ijk = NnT*(Ri_nu + Ri_nubar)*dE*dth*dphi;


                TrueEvents << th << "," << phi <<" , " << e << ", "<< N_ik  << "\n";

                std::cout << "Simulation: "<< AziMin << " " << AziMax << " " << th*180.0/TMath::Pi() << "," << e << ", " << dphi << " , "<< N_ik  << "\n";

                TrueHist->SetBinContent(i,k, N_ik); //Create histogram for  kth Pseudo-Experimens

            } // loop energy

        } // Loop zenith

    

    TrueEvents.close();
           
    return TrueHist;
}
*/

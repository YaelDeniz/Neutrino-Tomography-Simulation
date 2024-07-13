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
#include "Math/RootFinder.h"
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
#include "Earth3DModel.h"


using namespace std;


// Some Constants i need
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds

# define R_earth 6368.0 //km

//RESOLUTION FUNCITONS----------------------------------------------------------

//Energy Resolution Parameters
void AsimovObsSimulation::SetEnergyResolution(double A, double B){
    Ae = A;
    Be = B;}

//Energy Resolution function
double AsimovObsSimulation::PDFE(double Ereco , double Etrue ) {
    double SigmaE = Ae*Etrue +Be*sqrt(Etrue); //Detector Energy Resolution
    double pdfe = ROOT::Math::gaussian_pdf( Ereco, SigmaE, Etrue) ;
    return pdfe; } //Gaussian Distribution

//Angular (Zenith only) Resolution Parameters
void AsimovObsSimulation::SetAngularResolution(double A, double B){
    Ath = A;
    Bth = B;}

//Angular(Zenith only) Resolution function
double AsimovObsSimulation::PDFth(double thReco, double Etrue, double thTrue){
    double Sigmath = Ath + Bth/(sqrt(Etrue)); //Detector Angular Resolution
    double pdfth = ROOT::Math::gaussian_pdf( thReco, Sigmath, thTrue) ;
    return pdfth; } //Gaussian Distribution

double AsimovObsSimulation::PDFcth(double cthreco, double Etrue, double cthtrue){

    /// @param sigma - std. dev. in radians. Angular resolution
    /// kappa = 1/sigma^2 - VMF concentrarion parameter. 

    double sigmath = Ath + Bth/(sqrt(Etrue)); // Angular resolution in radians

    double cos_z = cthreco;
    double cos_z0 = cthtrue;

    double sin_z = sqrt(1 - pow(cos_z, 2));
    double sin_z0 = sqrt(1 - pow(cos_z0, 2));

    double kappa = 1/(sigmath*sigmath);

    double exp_arg = kappa * cos_z * cos_z0;

    double norm = kappa;
    
    if (kappa < 50)
        norm /= 2 * sinh(kappa);
    else
        exp_arg -= kappa;

    double arg = kappa * sin_z * sin_z0;
    double bess = 1;
    if (arg < 50)
        bess = TMath::BesselI0(arg);
    else 
    {
        exp_arg += arg;
        bess *=
            (1 + 1. / (8 * arg) + 4.5 / pow(8 * arg, 2) + 37.5 / pow(8 * arg, 3));
        bess /= sqrt(2 * M_PI * arg);
    }

    double expo = exp(exp_arg);

    return norm * expo * bess * sin_z;} //Von Mises Fisher

TH2D*  GetObsEvents2D(){  

    std::cout << "Asimov Approach Simulation for Observed Events" << std::endl;
    std::string PremFile = PremModel+".txt";    

    //Observed Variables:-----------------------------------------------------------------------------------------------
    
    //Observed bins
    int mbins = recobinsZen; //Number of angular bins of Observed event distribution
    int nbins = recobinsE; //Number of Energy bins of Observed event distribution
    
    //Energy[GeV]
    double ErecoMin      = EnuMin;
    double ErecoMax      = EnuMax;
    
    //Observed Angle
    double threcoMin   = ZenMin*TMath::Pi()/180.0;
    double threcoMax   = ZenMax*TMath::Pi()/180.0;

    double CthrecoMin = cos(threcoMax); //min cth = -1 
    double CthrecoMax = cos(threcoMin); // max cth = 0
    

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = truebinsZen; // Number of  angular bins of True event distribution
    int jbins = truebinsAzi;
    int kbins = truebinsE; // Number of  energy bins of True event distribution
    
    //True Energy[GeV]
    TF1 *fmin = new TF1("fmin","-1*x + [0]-4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmin->SetParameter(0,ErecoMin);
    fmin->SetParameter(1,Ae);
    fmin->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_min;
    finder_min.Solve(*fmin,0,100);

    TF1 *fmax = new TF1("fmax","-1*x + [0] + 4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmax->SetParameter(0,ErecoMax);
    fmax->SetParameter(1,Ae);
    fmax->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_max;
    finder_max.Solve(*fmax,0,100);

    double EtrueMin = finder_min.Root();
    double EtrueMax = finder_max.Root();

    //True Angle-
    double thtrueMin_a= (threcoMin)-4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMin_b= 0;

    double thtrueMax_a= (threcoMax)+4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMax_b= TMath::Pi();

    // True Angular Range
    double thtrueMin = max(thtrueMin_b,thtrueMin_a);
    double thtrueMax = min(thtrueMax_b,thtrueMax_a);

    double cthtrueMin = cos(thtrueMax);
    double cthtrueMax = cos(thtrueMin);

    //Azimuthal 
    double phiTrueMin = AziMin*TMath::Pi()/180.0;
    double phiTrueMax = AziMin*TMath::Pi()/180.0;

    //HISTOGRAMS--------------------------------------------------------------------------------------------------------
     //True Histogram        
    AsimovSimulation TrueSimulation;

    TrueSimulation.PremModel = PremModel;
    TrueSimulation.MantleAnomaly = Pile;
    TrueSimulation.SetIntervals(thtrueMin,thtrueMax,phiTrueMin,phiTrueMax,EtrueMin,EtrueMax);
    TrueSimulation.SetBinning(ibins,jbins,kbins);
    TrueSimulation.SetExposure(NnT);
    TrueSimulation.flvf=nuflv;
    TH2D * TrueHist = TrueSimulation.GetTrueEvents2D();

    //Observed distribution:
    TH2D* RecoHist = new TH2D("RecoHist","Observed Events Histogram",mbins,CthrecoMin,CthrecoMax,nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

    std::cout << "Observed variable region ["<< EOmin << "-"<<EOmax<< "]*["<<etaOmin<< "-"<<etaOmax <<"]" << std::endl;

    //std::cout << "Extended variable region ["<< Emin << "-"<<Emax<< "]*["<<etamin_a<< "-"<<etamax_a <<"]" << std::endl;

    std::cout << "True variable region ["<< Emin << "-"<<Emax<< "]*["<<etamin<< "-"<<etamax <<"]" << std::endl;

   //Calculation of observed events------------------------------------------------------------------------------------- 
    for (int m = 1; m <= mbins; ++m) //Observed Angular bins
    {

        double cthreco = RecoHist->GetXaxis()->GetBinCenter(m); // Angular bin from the "Observed" histogram
        double threco = acos(cthreco)
        double dcthreco = RecoHist->GetXaxis()->GetBinWidth(m); // Angular bin width of the "Observed" histogram
       
        for (int n = 1; n <=nbins ; ++n) //Observed Energy bins
        {
     
            double Ereco = RecoHist->GetYaxis()->GetBinCenter(n); //  Energy bin from the "Observed" histogram
            double dEreco = RecoHist->GetYaxis()->GetBinWidth(n); //  Energy bin width of the "Observed" histogram
       
            double Nsmear = 0;   // Reset total number of events.
            
            for(int i=1; i<= ibins ; i++) //Loop in cosT
            {    
                // Get cos(theta_z) from bin center, It fix a value of L
                double cthtrue = TrueHist->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
                double thtrue = acos(cthtrue);
                double dcthtrue = TrueHist->GetXaxis()->GetBinWidth(i); 

                for (int k = 1; k <= kbins; ++k)
                { 
                    //NUMERICAL INTEGRATION

                    double Etrue = TrueHist->GetYaxis()->GetBinCenter(k); //< This will defined a constant L por different values of ct provided Dct is Small
                    double dEtrue = TrueHist->GetYaxis()->GetBinWidth(k); 

                    double NikTrue = TrueHist->GetBinContente(i,k);

                    double dEdcthReco = dEreco*dcthreco;
                    double Nsmear += PDFcth(cthreco,Etrue,cthtrue)*PDFE(Ereco,Etrue)*(NikTrue)*(dEdcthReco);
    
                } // loop e

            } // Loop ct
            
            double Nmn = Nsmear;

            RecoHist->SetBinContent(m,n,Nmn); // Set events in Observed histogram. 

         
        } // loop in Ereco

    } // loop in CthReco

   delete TrueHist;
   
   return RecoHist;}


TH2D*  GetObsEvents2D(){  

    std::cout << "Asimov Approach Simulation for Observed Events" << std::endl;
    std::string PremFile = PremModel+".txt";    

    //Observed Variables:-----------------------------------------------------------------------------------------------
    
    //Observed bins
    int mbins = recobinsZen; //Number of angular bins of Observed event distribution
    int nbins = recobinsE; //Number of Energy bins of Observed event distribution
    
    //Energy[GeV]
    double ErecoMin      = EnuMin;
    double ErecoMax      = EnuMax;
    
    //Observed Angle
    double threcoMin   = ZenMin*TMath::Pi()/180.0;
    double threcoMax   = ZenMax*TMath::Pi()/180.0;

    double CthrecoMin = cos(threcoMax); //min cth = -1 
    double CthrecoMax = cos(threcoMin); // max cth = 0
    
    //Azimuthal 
    double dAz = (AziMax-AziMin)*TMath::Pi()/180.0;

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = truebinsZen; // Number of  angular bins of True event distribution
    int jbins = truebinsE; // Number of  energy bins of True event distribution
    
    //True Energy[GeV]
    TF1 *fmin = new TF1("fmin","-1*x + [0]-4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmin->SetParameter(0,ErecoMin);
    fmin->SetParameter(1,Ae);
    fmin->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_min;
    finder_min.Solve(*fmin,0,100);

    TF1 *fmax = new TF1("fmax","-1*x + [0] + 4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmax->SetParameter(0,ErecoMax);
    fmax->SetParameter(1,Ae);
    fmax->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_max;
    finder_max.Solve(*fmax,0,100);

    double EtrueMin = finder_min.Root();
    double EtrueMax = finder_max.Root();

    //True Angle-
    double thtrueMin_a= (threcoMin)-4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMin_b= 0;

    double thtrueMax_a= (threcoMax)+4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMax_b= TMath::Pi();

    // True Angular Range
    double thtrueMin = max(thtrueMin_b,thtrueMin_a);
    double thtrueMax = min(thtrueMax_b,thtrueMax_a);

    double cthtrueMin = cos(thtrueMax);
    double cthtrueMax = cos(thtrueMin);

    //HISTOGRAMS--------------------------------------------------------------------------------------------------------
    //True Histogram        
    TH2D* TrueHist = new TH2D("TrueHist","True Event Histogram",ibins,cthtrueMin,cthtrueMax,jbins,EtrueMin,EtrueMax); // xbins correspond to energy values and ybins to zenith angle cost
    
    //Observed distribution:
    TH2D* RecoHist = new TH2D("RecoHist","Observed Events Histogram",mbins,CthrecoMin,CthrecoMax,nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

    std::cout << "Observed variable region ["<< EOmin << "-"<<EOmax<< "]*["<<etaOmin<< "-"<<etaOmax <<"]" << std::endl;

    //std::cout << "Extended variable region ["<< Emin << "-"<<Emax<< "]*["<<etamin_a<< "-"<<etamax_a <<"]" << std::endl;

    std::cout << "True variable region ["<< Emin << "-"<<Emax<< "]*["<<etamin<< "-"<<etamax <<"]" << std::endl;

    
    //NEUTRINO OSCILLATION PROB-----------------------------------------------------------------------------------------

    // Neutrino flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino
    // Particle type
    int nu         = 1; //neutrino
    int nubar      = -1; //antineutrino

    // Model Set up: 
    OscProb::PMNS_Fast PMNS_H; // Create PMNS objects
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
     MyEarthModel.SetPile( Pile, PileShape );
     MyEarthModel.SetLayerProp(PremLayer, DensityContrast, ChemicalContrats);


   //Calculation of observed events------------------------------------------------------------------------------------- 
    for (int m = 1; m <= mbins; ++m) //Observed Angular bins
    {

        double cthreco = RecoHist->GetXaxis()->GetBinCenter(m); // Angular bin from the "Observed" histogram
        double dcthreco = RecoHist->GetXaxis()->GetBinWidth(m); // Angular bin width of the "Observed" histogram
       
        for (int n = 1; n <=nbins ; ++n) //Observed Energy bins
        {
     
            double Ereco = RecoHist->GetYaxis()->GetBinCenter(n); //  Energy bin from the "Observed" histogram
            double dEreco = RecoHist->GetYaxis()->GetBinWidth(n); //  Energy bin width of the "Observed" histogram
       
            double Ntot = 0;   // Reset total number of events.
            
            for(int i=1; i<= ibins ; i++) //Loop in cosT
            {    
                // Get cos(theta_z) from bin center, It fix a value of L
                double cthtrue = TrueHist->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
                double dcthtrue = TrueHist->GetXaxis()->GetBinWidth(i); 

                if(cosEta < -1 || cosEta > 1) break; // Skip if cosT is unphysical 
                
                MyEarthModel.SetDirection(cth,phi); 
                std::vector<std::vector<double>> EarthPath = MyEarthModel.Create3DPath();

                l = EarthPath[0][0]; //Length
                d = EarthPath[0][1]; //Density
                z = EarthPath[0][2]; //ZoA 
            
                PMNS_H.SetPath(l,d,z);
                for (int seg = 1; seg < EarthPath.size(); seg++) 
                { 
                    l = EarthPath[seg][0];
                    d = EarthPath[seg][1];
                    z = EarthPath[seg][2];
                    PMNS_H.AddPath(l,d,z);
                 }

                for (int j = 1; j <= jbins; ++j)
                { 
                    //NUMERICAL INTEGRATION

                    double Etrue = TrueHist->GetYaxis()->GetBinCenter(j); //< This will defined a constant L por different values of ct provided Dct is Small
                    double dEtrue = TrueHist->GetYaxis()->GetBinWidth(j); 

                    double logEi = log10(Etrue);
                
                    //Neutrino Flux interpolation
                    double logdPsiMu = muflux->Interpolate(logEi,cthtrue);
                    double logdPsiMub = mubflux->Interpolate(logEi,cthtrue);
                    double logdPsiE = eflux->Interpolate(logEi,cthtrue);
                    double logdPsiEb = ebflux->Interpolate(logEi,cthtrue);

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
                    
                    //Antineutrino contribution
                    PMNS_H.SetIsNuBar(true); 
                    double Ri_eb=XSec(e,nubar)*(PMNS_H.Prob(nue,flvf,e)*dPsiEbardEdct ); //Electron anti-neutrino contribution
                    double Ri_mub=XSec(e,nubar)*( PMNS_H.Prob(numu,flvf, e)*dPsiMubardEdct ); //Muon anti-neutrino contribution
                    double Ri_nubar = Ri_eb + Ri_mub;

                    //Convolution

                    double dEdcthTrue =dEtrue*dcthtrue;
                    double dEdcthReco =dEreco*dcthreco;
                    double N_ij = NnT*PDFth(thReco,Etrue,thTrue)*PDFE(Ereco,Etrue)*(Ri_nu+Ri_nubar)*(dEdcthTrue)*(dEdcthReco)*dAz;

                    //std::cout << Ri_nu << " " << Ri_nubar << " " << N_ij << std::endl;

                    Ntot += N_ij;
                    
                } // loop e

            } // Loop ct
            
            double N_mn = Ntot;

            RecoHist->SetBinContent(m,n,No_mn); // Set events in Observed histogram. 

         
        } // loop in Ereco

    } // loop in CthReco

   delete TrueHist;

   return RecoHist;}

TH2D*  OriginalApproach(std::string modelname, int flvf, double Region[], int Bins[],double Det_par[], double NnT)
{   
    std::cout << "Generating Asimov data set: Reimann Integration + Extended region" << std::endl;
    //Data Storage -----------------------------------------------------------------------------------------------------
    std::string model = "/home/dehy0499/OscProb/PremTables/"+modelname+".txt";
    std::string model_default = "/home/dehy0499/OscProb/PremTables/prem_default.txt";

    std::string location = "SimulationResults/AsimovData/" ;
    std::string Earthmodel= modelname;
    std::string title = "_AsiObsNextGen_"+std::to_string(flvf)+"_"+std::to_string(Bins[0])+"_"+std::to_string(Bins[2])+"_"+std::to_string(Region[0])+"_"+std::to_string(Region[1])+".csv";
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
    
    double dE_o = (EOmax - EOmin)/(nbins); //< Bin width
    
    double E_itv[] = {EOmin, EOmax};
    
    //Observed Angle
    double etaOmax   = (180.0-Region[2])*TMath::Pi()/180.0;
    
    double etaOmin   = (180.0-Region[3])*TMath::Pi()/180.0;
    
    double deta_o = (etaOmax - etaOmin)/(mbins); //< Bin width

    double Eta_itv[] = {etaOmin, etaOmax};
    
    //Azimuthal 
    double dAz = Region[4]*TMath::Pi()/180.0;

    //Observed distribution:
    TH2D* hObs = new TH2D("hObs","Observed Events",mbins,etaOmin,etaOmax,nbins,EOmin,EOmax); // xbins correspond to energy values and ybins to zenith angle cost
    

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = Bins[0]; // Number of  angular bins of True event distribution

    int jbins = Bins[1]; // Number of  energy bins of True event distribution
/*    
    //True Energy[GeV]
    
    double Emin=Region[0];

    double Emax=Region[1];
    
    double dE = (Emax - Emin)/(jbins); //< Bin width

    //True Angle
    double etamin= etaOmin;

    double etamax= etaOmax;

    double deta = (etamax - etamin)/(ibins); //< Bin width
    */

    
      //True Energy[GeV]
    //double Emin=(1/(1+4*a_E))*(EOmin);

    //double Emax=(1/(1-4*a_E))*(EOmax);
    TF1 *fmin = new TF1("fmin","-1*x + [0]-4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmin->SetParameter(0,EOmin);
    fmin->SetParameter(1,0.05);
    fmin->SetParameter(2,0.1);

    ROOT::Math::RootFinder finder_min;
    finder_min.Solve(*fmin,0,100);

    double Emin = finder_min.Root();

    TF1 *fmax = new TF1("fmax","-1*x + [0] + 4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmax->SetParameter(0,EOmax);
    fmax->SetParameter(1,0.05);
    fmax->SetParameter(2,0.1);

    ROOT::Math::RootFinder finder_max;
    finder_max.Solve(*fmax,0,100);

    double Emax = finder_max.Root();

    
    //NextGen Energy

   // double Emin=0.6;

   // double Emax=14.5;
    
    double dE = (Emax - Emin)/(jbins); //< Bin width

    //True Angle
    /*
    double etamin_a= (etaOmin)-4*(a_eta/sqrt(Emin));
    double etamin_b= 0;

    double etamax_a= (etaOmax)+4*(a_eta/sqrt(Emin));
    double etamax_b= TMath::Pi();
    */


    double etamin_a= (etaOmin)-4*(0.04 + 0.18/sqrt(Emin));
    double etamin_b= 0;

    double etamax_a= (etaOmax)+4*(0.04 + 0.18/sqrt(Emin));
    double etamax_b= TMath::Pi();


    double etamin = max(etamin_b,etamin_a);

    double etamax = min(etamax_b,etamax_a);

    double deta = (etamax - etamin)/(ibins); //< Bin width
    
    
    //True distribution        
    TH2D* hTrue = new TH2D("hTrue","True Events",jbins,etamin,etamax,ibins,Emin,Emax); // xbins correspond to energy values and ybins to zenith angle cost


    std::cout << "Observed variable region ["<< EOmin << "-"<<EOmax<< "]*["<<etaOmin<< "-"<<etaOmax <<"]" << std::endl;

    //std::cout << "Extended variable region ["<< Emin << "-"<<Emax<< "]*["<<etamin_a<< "-"<<etamax_a <<"]" << std::endl;

    std::cout << "True variable region ["<< Emin << "-"<<Emax<< "]*["<<etamin<< "-"<<etamax <<"]" << std::endl;

    
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

    //Earth model for neutrino propagation 
    OscProb::PremModel prem_model(model);
    OscProb::PremModel prem(model_default); //Default PREM table

    //Limit region to turn off model
    double R_lim = 3500.0; // Distance from the center of the Earth , 3500 Km CMB
    double EtaLim = TMath::Pi() - TMath::ASin( (R_lim)/R_earth ) ;


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
                
                double cosEta =cos(eta);

                if(cosEta < -1 || cosEta > 1) break; // Skip if cosT is unphysical 
                
                

                
                if (eta <= EtaLim ) 
                {

                    //std::cout << eta << " " << EtaLim << " --- " << std::endl;
                    double cosEta =cos(eta);
                     
                    prem_model.FillPath(cosEta); // Fill paths from PREM model

                    PMNS_Ho.SetPath(prem_model.GetNuPath()); // Set paths in OscProb  

                }
                else

                {

                    //std::cout << eta << " " << EtaLim << " *** " << std::endl;

                    prem.FillPath(cosEta); // Fill paths from PREM model
                
                    PMNS_Ho.SetPath(prem.GetNuPath()); // Set paths in OscProb  
                }
/*
                std::cout << " Things is off " << std::endl;

                prem_model.FillPath(cosEta); // Fill paths from PREM model
                
                // Set paths in OscProb  
                PMNS_Ho.SetPath(prem_model.GetNuPath()); //STANDARD EARTH (DESCRIBED BY PREM MODEL)
*/

                for (int j = 1; j <= jbins; ++j)
                { 
                    //NUMERICAL INTEGRATION

                    double e = hTrue->GetYaxis()->GetBinCenter(j); //< This will defined a constant L por different values of ct provided Dct is Small

                    //Neutrino
                    PMNS_Ho.SetIsNuBar(false); 
                    double Ri_nu = XSec(e,nu)*( PMNS_Ho.Prob(numu, flvf, e)*dPsiMudE.Eval(e) + PMNS_Ho.Prob(nue,flvf,e)*dPsiEdE.Eval(e) );
                
                    //Antineutrino contribution
                    PMNS_Ho.SetIsNuBar(true); 
                    double Ri_nubar = XSec(e,nubar)*( PMNS_Ho.Prob(numu,flvf, e)*dPsiMubardE.Eval(e) + PMNS_Ho.Prob(nue,flvf,e)*dPsiEbardE.Eval(e) ); 


                    double N_ij = NnT*PDF_angle(eta_o,e,eta,a_eta)*PDF_energy(e_o,e,a_E)*(Ri_nu+Ri_nubar)*(dE*deta)*(dE_o*deta_o)*dAz;

                    //std::cout << Ri_nu << " " << Ri_nubar << " " << N_ij << std::endl;
                    
                   

                    /*
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
                    */





                    //double N_ij=(NnT*(dAz)*(TMath::Pi()/180.0) )*dNHo_dAzdNnT;

                    Ntot += N_ij;
                    
                } // loop e

            } // Loop ct
            
            double No_mn = Ntot;

            //std::cout << eta_o << " " << e_o << " ****************************** " << No_mn << std::endl;
                    

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
/*
double PDF_energy(double e_o , double e, double a_e ) 
{

//double sigma_ene = a_e*e;

double sigma_ene = 0.05*e +0.1*sqrt(e);

double pdf_ene = ROOT::Math::gaussian_pdf( e_o, sigma_ene, e) ;
 
return pdf_ene; 



}

//double w_eta(double a_eta, double eta, double E, std::vector<double> eta_o, double Eta[])

double PDF_angle(double eta_o, double e, double eta, double a_eta)
{

//double sigma_angle = a_eta/(sqrt(e));

double sigma_angle = 0.04 + 0.18/(sqrt(e)); //next gen - Maderer

double pdf_angle = ROOT::Math::gaussian_pdf     ( eta_o, sigma_angle, eta) ;
 
return pdf_angle; 


}
*/
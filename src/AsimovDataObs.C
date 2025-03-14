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

//#include "PremModel.h"
//#include "PMNS_Fast.h"


// My Math Tools
#include "MathTools.h"
#include "PhyTools.h"
#include "GetTrueEvents.h"
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

//Angular Resolution Parameters (rad)
void AsimovObsSimulation::SetAngularResolution(double A, double B){
    Ath = A;
    Bth = B;}

//Angular(Zenith only) Resolution function
double AsimovObsSimulation::Gaussth(double threco, double Etrue, double thtrue){
    double Sigmath = Ath + Bth/(sqrt(Etrue)); //Detector Angular Resolution
    double pdfth = ROOT::Math::gaussian_pdf( threco, Sigmath, thtrue) ;
    return pdfth; } //Gaussian Distribution

double AsimovObsSimulation::VMFth(double threco, double Etrue, double thtrue){

    /// @param sigma - std. dev. in radians. Angular resolution
    /// kappa = 1/sigma^2 - VMF concentrarion parameter. 

    double sigmath = Ath + Bth/(sqrt(Etrue)); // Angular resolution in radians

    //double cos_z = cthreco;
    //double cos_z0 = cthtrue;

    double cos_z  = cos(threco);
    double cos_z0 = cos(thtrue);

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

// Observed events in cth binning
std::vector < TH2D* >  AsimovObsSimulation::GetObsEvents3Dcth(){  

    std::cout << "Asimov Approach Simulation for Observed Events" << std::endl;


    //Observed Variables:-----------------------------------------------------------------------------------------------
    
    //Observed bins
    int mbins = recobinsZen; //Number of angular bins of Observed event distribution
    int nbins = recobinsE; //Number of Energy bins of Observed event distribution 
    
    //Energy[GeV]
    double ErecoMin    = EnuMin;
    double ErecoMax    = EnuMax;
    
    //Observed Angle
    double threcoMin   = ZenMin;
    double threcoMax   = ZenMax;

    double CthrecoMin  = cos(threcoMax*TMath::Pi()/180.0); //min cth = -1 
    double CthrecoMax  = cos(threcoMin*TMath::Pi()/180.0); // max cth = 0
    

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = truebinsZen; // Number of  angular bins of True event distribution
    int jbins = truebinsAzi;
    int kbins = truebinsE; // Number of  energy bins of True event distribution
    
    //True Energy[GeV]
    TF1 *fmin = new TF1("fmin","-1*x + [0] -4*( [1]*x+[2]*sqrt(x) )", 0, 100);
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
    double thtrueMin_a= (threcoMin*TMath::Pi()/180.0)-4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMin_b= 0;

    double thtrueMax_a= (threcoMax*TMath::Pi()/180.0)+4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMax_b= TMath::Pi();

    // True Angular Range
    double thtrueMin = max(thtrueMin_b,thtrueMin_a)*180.0/TMath::Pi(); //From Radians to degrees
    double thtrueMax = min(thtrueMax_b,thtrueMax_a)*180.0/TMath::Pi(); // From radias to degrees

    //Azimuthal 
    double phiTrueMin = AziMin;
    double phiTrueMax = AziMax;


    //std::cout << "TrueVariable | E=[ " << EtrueMin << " , " << EtrueMax << " ] Zenith= [" << thtrueMin << " , " << thtrueMax << "]" << std::endl;
    //std::cout << "RecoVariable | E=[ " << ErecoMin << " , " << ErecoMax << " ] Zenith= [" << threcoMin << " , " << threcoMax << "]" << std::endl;


    //HISTOGRAMS--------------------------------------------------------------------------------------------------------
     //True Histogram        
    AsimovSimulation TrueSimulation;

    TrueSimulation.PremTable = ThePremTable; //Prem Model
    TrueSimulation.HondaTable = TheHondaTable; // Honda Table
    TrueSimulation.SetDetectorPosition(xyzTelescope); //Detector location



    //If LLVP set geometry and properties
    TrueSimulation.MantleAnomaly = PileInModel;
    TrueSimulation.AnomalyShape= ShapeOfPile;
    TrueSimulation.PileHeight = ThePileHight; // LLVP mas height 
    TrueSimulation.aperture=ThePileAperture; //Aperture
    TrueSimulation.PileDensityContrast = ThePileDensityContrast;
    TrueSimulation.PileChemContrast = ThePileChemicalContrast;




    TrueSimulation.ModifyLayer(PremLayer,DensityContrast,ChemicalContrats);
    TrueSimulation.SetIntervals(thtrueMin,thtrueMax,phiTrueMin,phiTrueMax,EtrueMin,EtrueMax);
    TrueSimulation.SetBinning(ibins,jbins,kbins);
    TrueSimulation.SetExposure(MT);
    TrueSimulation.flvf=flvf;


    std::vector<TH2D*> TrueHist= TrueSimulation.GetTrueEvents3D(); //cos(theta) binning



    std::vector<TH2D*> HistVec;
    TH2D* RecoHist[jbins]; 

    //Observed distribution:
    //TH2D* RecoHist = new TH2D("RecoHist","Observed Events Histogram",mbins,CthrecoMin,CthrecoMax,nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

    //Loop Over Histogram in Azimuth
    for (int nphi = 0 ; nphi < TrueHist.size(); ++nphi)
    {

        RecoHist[nphi] = new TH2D(Form("RecoHist%d",nphi),Form("ObsOscillogram%d",nphi),
                                  mbins,CthrecoMin,CthrecoMax,
                                  nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

        //Loop over Reconstructed cth
        for (int m = 1; m <= mbins; ++m) //Observed Angular bins
        {

            double cthreco = RecoHist[nphi]->GetXaxis()->GetBinCenter(m); // Angular bin from the "Observed" histogram
            double threco  =  acos(cthreco);
            double sthreco =  sin(threco); 
            double dcthreco= RecoHist[nphi]->GetXaxis()->GetBinWidth(m); // Angular bin width of the "Observed" histogram

            //Loop over Reconstructed energy
            for (int n = 1; n <=nbins ; ++n) //Observed Energy bins
            {
         
                double Ereco = RecoHist[nphi]->GetYaxis()->GetBinCenter(n); //  Energy bin from the "Observed" histogram
                double dEreco = RecoHist[nphi]->GetYaxis()->GetBinWidth(n); //  Energy bin width of the "Observed" histogram
           
                double Nsmear = 0;   // Reset total number of events.
                
                for(int i=1; i<= ibins ; i++) //Loop in cosT
                {    
                    // Get cos(theta_z) from bin center, It fix a value of L
                    double cthtrue = TrueHist[nphi]->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
                    double thtrue = acos(cthtrue);
                    double dcthtrue = TrueHist[nphi]->GetXaxis()->GetBinWidth(i); 

                    for (int k = 1; k <= kbins; ++k)
                    { 
                        //NUMERICAL INTEGRATION


                        double Etrue = TrueHist[nphi]->GetYaxis()->GetBinCenter(k); //< This will defined a constant L por different values of ct provided Dct is Small
                        double dEtrue = TrueHist[nphi]->GetYaxis()->GetBinWidth(k); 

                        double NikTrue = TrueHist[nphi]->GetBinContent(i,k);

                        double dEdcthReco = dEreco*dcthreco;


                         ZenithBinSize = TrueHist[nphi]->GetXaxis()->GetXmax();
                         AzimuthBinSize =0;
                         EnergyBinSize = TrueHist[nphi]->GetYaxis()->GetXmax();
                       
                        Nsmear += (VMFth(threco,Etrue,thtrue)/sthreco)*PDFE(Ereco,Etrue)*(NikTrue)*(dEdcthReco);
        
                    } // loop e

                } // Loop ct
                
                double Nmn = Nsmear;

                RecoHist[nphi]->SetBinContent(m,n,Nmn); // Set events in Observed histogram. 
             
            } // End reconstructed E loop 

        } // End reconstructed cth loop

        HistVec.push_back(RecoHist[nphi]); 

    } // End Histogram loop
   
   return HistVec;}


// Observed events in th binning
std::vector < TH2D* >  AsimovObsSimulation::GetObsEvents3Dth(){  

    std::cout << "Asimov Approach Simulation for Observed Events" << std::endl;


    //Observed Variables:-----------------------------------------------------------------------------------------------
    
    //Observed bins
    int mbins = recobinsZen; //Number of angular bins of Observed event distribution
    int nbins = recobinsE; //Number of Energy bins of Observed event distribution 
    
    //Energy[GeV]
    double ErecoMin      = EnuMin;
    double ErecoMax      = EnuMax;
    
    //Observed Angle (rad)
    double threcoMin   = ZenMin*TMath::Pi()/180.0;
    double threcoMax   = ZenMax*TMath::Pi()/180.0;

    //double CthrecoMin = cos(threcoMax*TMath::Pi()/180.0); //min cth = -1 
    //double CthrecoMax = cos(threcoMin*TMath::Pi()/180.0); // max cth = 0
    

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = truebinsZen; // Number of  angular bins of True event distribution
    int jbins = truebinsAzi;
    int kbins = truebinsE; // Number of  energy bins of True event distribution
    
    //True Energy[GeV]
    TF1 *fmin = new TF1("fmin","-1*x + [0] -4*( [1]*x+[2]*sqrt(x) )", 0, 100);
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
    double thtrueMin_a= threcoMin - 4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMin_b= 0;

    double thtrueMax_a= threcoMax + 4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMax_b= TMath::Pi();

    // True Angular Range
    double thtrueMin = max(thtrueMin_b,thtrueMin_a)*180.0/TMath::Pi(); //From Radians to degrees
    double thtrueMax = min(thtrueMax_b,thtrueMax_a)*180.0/TMath::Pi(); // From radias to degrees

    //Azimuthal 
    double phiTrueMin = AziMin;
    double phiTrueMax = AziMax;
    

    //HISTOGRAMS--------------------------------------------------------------------------------------------------------
     //True Histogram        
    AsimovSimulation TrueSimulation;

    TrueSimulation.PremTable = ThePremTable; //Prem Model
    TrueSimulation.HondaTable = TheHondaTable; // Honda Table
    TrueSimulation.SetDetectorPosition(xyzTelescope); //Detector location

    //If LLVP set geometry and properties
    TrueSimulation.MantleAnomaly = PileInModel;
    TrueSimulation.AnomalyShape= ShapeOfPile;
    TrueSimulation.PileHeight = ThePileHight; // LLVP mas height 
    TrueSimulation.aperture=ThePileAperture; //Aperture
    TrueSimulation.PileDensityContrast = ThePileDensityContrast;
    TrueSimulation.PileChemContrast = ThePileChemicalContrast;


    TrueSimulation.ModifyLayer(PremLayer,DensityContrast,ChemicalContrats);
    TrueSimulation.SetIntervals(thtrueMin,thtrueMax,phiTrueMin,phiTrueMax,EtrueMin,EtrueMax);
    TrueSimulation.SetBinning(ibins,jbins,kbins);
    TrueSimulation.SetExposure(MT);
    TrueSimulation.flvf=flvf; 

    //Binning in zenith angle
    std::vector<TH2D*> TrueHist= TrueSimulation.GetTrueEvents3Dth();


    std::vector<TH2D*> HistVec;
    TH2D* RecoHist[jbins]; 

    //Observed distribution:
    //TH2D* RecoHist = new TH2D("RecoHist","Observed Events Histogram",mbins,CthrecoMin,CthrecoMax,nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

    //Loop Over Histogram in Azimuth
    for (int nphi = 0 ; nphi < TrueHist.size(); ++nphi)
    {

        RecoHist[nphi] = new TH2D(Form("RecoHist%d",nphi),Form("ObsOscillogram%d",nphi),
                                  mbins,threcoMin,threcoMax,
                                  nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

        //Loop over Reconstructed cth
        for (int m = 1; m <= mbins; ++m) //Observed Angular bins
        {

            double threco = RecoHist[nphi]->GetXaxis()->GetBinCenter(m); // Angular bin from the "Observed" histogram
            double dthreco= RecoHist[nphi]->GetXaxis()->GetBinWidth(m); // Angular bin width of the "Observed" histogram

            double cthreco  =  cos(threco);

            //Loop over Reconstructed energy
            for (int n = 1; n <=nbins ; ++n) //Observed Energy bins
            {
         
                double Ereco = RecoHist[nphi]->GetYaxis()->GetBinCenter(n); //  Energy bin from the "Observed" histogram
                double dEreco = RecoHist[nphi]->GetYaxis()->GetBinWidth(n); //  Energy bin width of the "Observed" histogram
           
                double Nsmear = 0;   // Reset total number of events.
                
                for(int i=1; i<= ibins ; i++) //Loop in cosT
                {    
                    // Get cos(theta_z) from bin center, It fix a value of L
                    double thtrue = TrueHist[nphi]->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
                    double dthtrue = TrueHist[nphi]->GetXaxis()->GetBinWidth(i); 

                    double cthtrue = cos(thtrue);

                    for (int k = 1; k <= kbins; ++k)
                    { 
                        //NUMERICAL INTEGRATION

                        double Etrue = TrueHist[nphi]->GetYaxis()->GetBinCenter(k); //< This will defined a constant L por different values of ct provided Dct is Small
                        double dEtrue = TrueHist[nphi]->GetYaxis()->GetBinWidth(k); 

                        double NikTrue = TrueHist[nphi]->GetBinContent(i,k);

                        double dEdthReco = dEreco*dthreco;
                       
                        Nsmear += VMFth(threco,Etrue,thtrue)*PDFE(Ereco,Etrue)*(NikTrue)*(dEdthReco);

                    } // loop e

                } // Loop th
                
                double Nmn = Nsmear;

                RecoHist[nphi]->SetBinContent(m,n,Nmn); // Set events in Observed histogram. 
             
            } // End reconstructed E loop 

        } // End reconstructed th loop

        HistVec.push_back(RecoHist[nphi]); 

    } // End Histogram loop

    std::cout << "TrueVariable | E=[ " << EtrueMin << " , " << EtrueMax << " ] Zenith= [" << thtrueMin << " , " << thtrueMax << "]" << std::endl;
    std::cout << "RecoVariable | E=[ " << ErecoMin << " , " << ErecoMax << " ] Zenith= [" << threcoMin << " , " << threcoMax << "]" << std::endl;
   
   return HistVec;}

// Observed events in th binning
std::vector < TH2D* >  AsimovObsSimulation::GetObsEvents3Dthlog10E(){  

    std::cout << "Asimov Approach Simulation for Observed Events" << std::endl;


    //Observed Variables:-----------------------------------------------------------------------------------------------
    
    //Observed bins
    int mbins = recobinsZen; //Number of angular bins of Observed event distribution
    int nbins = recobinsE; //Number of Energy bins of Observed event distribution 
    
    //Energy[GeV]
    double ErecoMin      = EnuMin;
    double ErecoMax      = EnuMax;
    
    //Observed Angle (rad)
    double threcoMin   = ZenMin*TMath::Pi()/180.0;
    double threcoMax   = ZenMax*TMath::Pi()/180.0;

    //double CthrecoMin = cos(threcoMax*TMath::Pi()/180.0); //min cth = -1 
    //double CthrecoMax = cos(threcoMin*TMath::Pi()/180.0); // max cth = 0
    

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = truebinsZen; // Number of  angular bins of True event distribution
    int jbins = truebinsAzi;
    int kbins = truebinsE; // Number of  energy bins of True event distribution
    
    //True Energy[GeV]
    TF1 *fmin = new TF1("fmin","-1*x + [0] -4*( [1]*x+[2]*sqrt(x) )", 0, 100);
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
    double thtrueMin_a= threcoMin - 4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMin_b= 0;

    double thtrueMax_a= threcoMax + 4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMax_b= TMath::Pi();

    // True Angular Range
    double thtrueMin = max(thtrueMin_b,thtrueMin_a)*180.0/TMath::Pi(); //From Radians to degrees
    double thtrueMax = min(thtrueMax_b,thtrueMax_a)*180.0/TMath::Pi(); // From radias to degrees

    //Azimuthal 
    double phiTrueMin = AziMin;
    double phiTrueMax = AziMax;
    

    //HISTOGRAMS--------------------------------------------------------------------------------------------------------
     //True Histogram        
    AsimovSimulation TrueSimulation;

    TrueSimulation.PremTable = ThePremTable; //Prem Model
    TrueSimulation.HondaTable = TheHondaTable; // Honda Table
    TrueSimulation.SetDetectorPosition(xyzTelescope); //Detector location

    //If LLVP set geometry and properties
    TrueSimulation.MantleAnomaly = PileInModel;
    TrueSimulation.AnomalyShape= ShapeOfPile;
    TrueSimulation.PileHeight = ThePileHight; // LLVP mas height 
    TrueSimulation.aperture=ThePileAperture; //Aperture
    TrueSimulation.PileDensityContrast = ThePileDensityContrast;
    TrueSimulation.PileChemContrast = ThePileChemicalContrast;


    TrueSimulation.ModifyLayer(PremLayer,DensityContrast,ChemicalContrats);
    TrueSimulation.SetIntervals(thtrueMin,thtrueMax,phiTrueMin,phiTrueMax,EtrueMin,EtrueMax);
    TrueSimulation.SetBinning(ibins,jbins,kbins);
    TrueSimulation.SetExposure(MT);
    TrueSimulation.flvf=flvf; 

    //Binning in zenith angle
    std::vector<TH2D*> TrueHist= TrueSimulation.GetTrueEvents3Dthlog10E();


    std::vector<TH2D*> HistVec;
    TH2D* RecoHist[jbins]; 

    //Observed distribution:
    //TH2D* RecoHist = new TH2D("RecoHist","Observed Events Histogram",mbins,CthrecoMin,CthrecoMax,nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

    //Loop Over Histogram in Azimuth
    for (int nphi = 0 ; nphi < TrueHist.size(); ++nphi)
    {

        RecoHist[nphi] = new TH2D(Form("RecoHist%d",nphi),Form("ObsOscillogram%d",nphi),
                                  mbins,threcoMin,threcoMax,
                                  nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

        //Loop over Reconstructed cth
        for (int m = 1; m <= mbins; ++m) //Observed Angular bins
        {

            double threco = RecoHist[nphi]->GetXaxis()->GetBinCenter(m); // Angular bin from the "Observed" histogram
            double dthreco= RecoHist[nphi]->GetXaxis()->GetBinWidth(m); // Angular bin width of the "Observed" histogram

            double cthreco  =  cos(threco);

            //Loop over Reconstructed energy
            for (int n = 1; n <=nbins ; ++n) //Observed Energy bins
            {
         
                double log10Ereco = RecoHist[nphi]->GetYaxis()->GetBinCenter(n); //  Energy bin from the "Observed" histogram
                double Ereco = pow(10, log10Ereco);
                double dlog10Ereco = RecoHist[nphi]->GetYaxis()->GetBinWidth(n); //  Energy bin width of the "Observed" histogram
           
                double Nsmear = 0;   // Reset total number of events.
                
                for(int i=1; i<= ibins ; i++) //Loop in cosT
                {    
                    // Get cos(theta_z) from bin center, It fix a value of L
                    double thtrue = TrueHist[nphi]->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
                    double dthtrue = TrueHist[nphi]->GetXaxis()->GetBinWidth(i); 

                    double cthtrue = cos(thtrue);

                    for (int k = 1; k <= kbins; ++k)
                    { 
                        //NUMERICAL INTEGRATION

                        double log10Etrue = TrueHist[nphi]->GetYaxis()->GetBinCenter(k); //< This will defined a constant L por different values of ct provided Dct is Small
                        double Etrue = pow(10, log10Etrue);
                        double dlog10Etrue = TrueHist[nphi]->GetYaxis()->GetBinWidth(k); 

                        double NikTrue = TrueHist[nphi]->GetBinContent(i,k);

                        double dlog10EdthReco = dlog10Ereco*dthreco;
                       
                        Nsmear += VMFth(threco,Etrue,thtrue)*PDFE(Ereco,Etrue)*(NikTrue)*(log(10)*Ereco)*(dlog10EdthReco);

                    } // loop e

                } // Loop th
                
                double Nmn = Nsmear;

                RecoHist[nphi]->SetBinContent(m,n,Nmn); // Set events in Observed histogram. 
             
            } // End reconstructed E loop 

        } // End reconstructed th loop

        HistVec.push_back(RecoHist[nphi]); 

    } // End Histogram loop

    std::cout << "TrueVariable | E=[ " << EtrueMin << " , " << EtrueMax << " ] Zenith= [" << thtrueMin << " , " << thtrueMax << "]" << std::endl;
    std::cout << "RecoVariable | E=[ " << ErecoMin << " , " << ErecoMax << " ] Zenith= [" << threcoMin << " , " << threcoMax << "]" << std::endl;
   
   return HistVec;}

// Observed events in th binning - Variable binning
std::vector < TH2D* >  AsimovObsSimulation::GetObsEvents3Dthvb(){  

   
    //Observed Variables:-----------------------------------------------------------------------------------------------
    
    //Observed bins
    int mbins = recobinsZen; //Number of angular bins of Observed event distribution
    //int nbins = recobinsE; //Number of Energy bins of Observed event distribution 
    
    //Energy[GeV]
    double ErecoMin      = EnuMin;
    double ErecoMax      = EnuMax;
    
    //Observed Angle (rad)
    double threcoMin   = ZenMin*TMath::Pi()/180.0;
    double threcoMax   = ZenMax*TMath::Pi()/180.0;

    //double CthrecoMin = cos(threcoMax*TMath::Pi()/180.0); //min cth = -1 
    //double CthrecoMax = cos(threcoMin*TMath::Pi()/180.0); // max cth = 0

    //Create bins of variable size

    std::vector<double> ErecoEdges; 
    ErecoEdges.push_back(ErecoMin);
    double Erecoi = ErecoMin;

    while (ErecoEdges.back() < ErecoMax) 
    {
        //Find center
        TF1 *frecoC = new TF1("fmin","-1*x + [0] + 0.5*( [1]*x+[2]*sqrt(x) )", 0, 100);
        frecoC->SetParameter(0,Erecoi);
        frecoC->SetParameter(1,Ae);
        frecoC->SetParameter(2,Be);
        ROOT::Math::RootFinder finder_center;
        finder_center.Solve(*frecoC ,0,100);
        double ErecoC = finder_center.Root();
        
        //Find next edge
        double Erecoii = ErecoC + 0.5*(Ae*ErecoC + Be*sqrt(ErecoC));
        
        ErecoEdges.push_back(Erecoii);
        
        Erecoi = Erecoii;
        
    }

    double ErecobinEdges [ErecoEdges.size()];
        
    for (int i=0 ; i < ErecoEdges.size(); i++ )
    {
          ErecobinEdges[i] = ErecoEdges[i];
            
    }
    
    int nbins = sizeof(ErecobinEdges)/sizeof(ErecobinEdges[0]) - 1; // Number of bins

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = truebinsZen; // Number of  angular bins of True event distribution
    int jbins = truebinsAzi;
    int kbins = truebinsE; // Number of  energy bins of True event distribution
    
    //True Energy[GeV]
    TF1 *fmin = new TF1("fmin","-1*x + [0] -4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmin->SetParameter(0,ErecoMin);
    fmin->SetParameter(1,Ae);
    fmin->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_min;
    finder_min.Solve(*fmin,0,100);

    TF1 *fmax = new TF1("fmax","-1*x + [0] + 4*( [1]*x+[2]*sqrt(x) )", 0, 100);
    fmax->SetParameter(0,ErecoEdges.back());
    fmax->SetParameter(1,Ae);
    fmax->SetParameter(2,Be);
    ROOT::Math::RootFinder finder_max;
    finder_max.Solve(*fmax,0,100);

    double EtrueMin = finder_min.Root();
    double EtrueMax = finder_max.Root();

    //True Angle-
    double thtrueMin_a= threcoMin - 4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMin_b= 0;

    double thtrueMax_a= threcoMax + 4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMax_b= TMath::Pi();

    // True Angular Range
    double thtrueMin = max(thtrueMin_b,thtrueMin_a)*180.0/TMath::Pi(); //From Radians to degrees
    double thtrueMax = min(thtrueMax_b,thtrueMax_a)*180.0/TMath::Pi(); // From radias to degrees

    //Azimuthal 
    double phiTrueMin = AziMin;
    double phiTrueMax = AziMax;
    

    //HISTOGRAMS--------------------------------------------------------------------------------------------------------
     //True Histogram        
    AsimovSimulation TrueSimulation;

    TrueSimulation.PremTable = ThePremTable; //Prem Model
    TrueSimulation.HondaTable = TheHondaTable; // Honda Table
    TrueSimulation.SetDetectorPosition(xyzTelescope); //Detector location

    //If LLVP set geometry and properties
    TrueSimulation.MantleAnomaly = PileInModel;
    TrueSimulation.AnomalyShape= ShapeOfPile;
    TrueSimulation.PileHeight = ThePileHight; // LLVP mas height 
    TrueSimulation.aperture=ThePileAperture; //Aperture
    TrueSimulation.PileDensityContrast = ThePileDensityContrast;
    TrueSimulation.PileChemContrast = ThePileChemicalContrast;


    TrueSimulation.ModifyLayer(PremLayer,DensityContrast,ChemicalContrats);
    TrueSimulation.SetIntervals(thtrueMin,thtrueMax,phiTrueMin,phiTrueMax,EtrueMin,EtrueMax);
    TrueSimulation.SetBinning(ibins,jbins,kbins);
    TrueSimulation.SetExposure(MT);
    TrueSimulation.flvf=flvf; 

    //Binning in zenith angle
    std::vector<TH2D*> TrueHist= TrueSimulation.GetTrueEvents3Dth();


    std::vector<TH2D*> HistVec;
    TH2D* RecoHist[jbins]; 

    //Observed distribution:
    //TH2D* RecoHist = new TH2D("RecoHist","Observed Events Histogram",mbins,CthrecoMin,CthrecoMax,nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

    //Loop Over Histogram in Azimuth
    for (int nphi = 0 ; nphi < TrueHist.size(); ++nphi)
    {

        RecoHist[nphi] = new TH2D(Form("RecoHist%d",nphi),Form("ObsOscillogram%d",nphi),
                                  mbins,threcoMin,threcoMax,
                                  nbins,ErecobinEdges); // xbins correspond to energy values and ybins to zenith angle cost

        //Loop over Reconstructed cth
        for (int m = 1; m <= mbins; ++m) //Observed Angular bins
        {

            double threco = RecoHist[nphi]->GetXaxis()->GetBinCenter(m); // Angular bin from the "Observed" histogram
            double dthreco= RecoHist[nphi]->GetXaxis()->GetBinWidth(m); // Angular bin width of the "Observed" histogram

            double cthreco  =  cos(threco);

            //Loop over Reconstructed energy
            for (int n = 1; n <=nbins ; ++n) //Observed Energy bins
            {
         
                double Ereco = RecoHist[nphi]->GetYaxis()->GetBinCenter(n); //  Energy bin from the "Observed" histogram
                double dEreco = RecoHist[nphi]->GetYaxis()->GetBinWidth(n); //  Energy bin width of the "Observed" histogram
           
                double Nsmear = 0;   // Reset total number of events.
                
                for(int i=1; i<= ibins ; i++) //Loop in cosT
                {    
                    // Get cos(theta_z) from bin center, It fix a value of L
                    double thtrue = TrueHist[nphi]->GetXaxis()->GetBinCenter(i); //< This will defined a constant L por different values of ct provided Dct is Small
                    double dthtrue = TrueHist[nphi]->GetXaxis()->GetBinWidth(i); 

                    double cthtrue = cos(thtrue);

                    for (int k = 1; k <= kbins; ++k)
                    { 
                        //NUMERICAL INTEGRATION

                        double Etrue = TrueHist[nphi]->GetYaxis()->GetBinCenter(k); //< This will defined a constant L por different values of ct provided Dct is Small
                        double dEtrue = TrueHist[nphi]->GetYaxis()->GetBinWidth(k); 

                        double NikTrue = TrueHist[nphi]->GetBinContent(i,k);

                        double dEdthReco = dEreco*dthreco;
                       
                        Nsmear += VMFth(threco,Etrue,thtrue)*PDFE(Ereco,Etrue)*(NikTrue)*(dEdthReco);

                    } // loop e

                } // Loop th
                
                double Nmn = Nsmear;

                RecoHist[nphi]->SetBinContent(m,n,Nmn); // Set events in Observed histogram. 
             
            } // End reconstructed E loop 

        } // End reconstructed th loop

        HistVec.push_back(RecoHist[nphi]); 

    } // End Histogram loop

    std::cout << "TrueVariable | E=[ " << EtrueMin << " , " << EtrueMax << " ] Zenith= [" << thtrueMin << " , " << thtrueMax << "]" << std::endl;
    std::cout << "RecoVariable | E=[ " << ErecoMin << " , " << ErecoMax << " ] Zenith= [" << threcoMin << " , " << threcoMax << "]" << std::endl;
   
   return HistVec;}

/*
TH2D*  AsimovObsSimulation::GetObsEvents2Dcth(){  

    std::cout << "Asimov Approach Simulation for Observed Events" << std::endl;


    //Observed Variables:-----------------------------------------------------------------------------------------------
    
    //Observed bins
    int mbins = recobinsZen; //Number of angular bins of Observed event distribution
    int nbins = recobinsE; //Number of Energy bins of Observed event distribution 
    
    //Energy[GeV]
    double ErecoMin      = EnuMin;
    double ErecoMax      = EnuMax;
    
    //Observed Angle
    double threcoMin   = ZenMin;
    double threcoMax   = ZenMax;

    double CthrecoMin = cos(threcoMax*TMath::Pi()/180.0); //min cth = -1 
    double CthrecoMax = cos(threcoMin*TMath::Pi()/180.0); // max cth = 0
    

    //True Variables:---------------------------------------------------------------------------------------------------

    //True bins
    int ibins = truebinsZen; // Number of  angular bins of True event distribution
    int jbins = truebinsAzi;
    int kbins = truebinsE; // Number of  energy bins of True event distribution
    
    //True Energy[GeV]
    TF1 *fmin = new TF1("fmin","-1*x + [0] -4*( [1]*x+[2]*sqrt(x) )", 0, 100);
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
    double thtrueMin_a= (threcoMin*TMath::Pi()/180.0)-4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMin_b= 0;

    double thtrueMax_a= (threcoMax*TMath::Pi()/180.0)+4*(Ath + Bth/sqrt(EtrueMin));
    double thtrueMax_b= TMath::Pi();

    // True Angular Range
    double thtrueMin = max(thtrueMin_b,thtrueMin_a)*180.0/TMath::Pi(); //From Radians to degrees
    double thtrueMax = min(thtrueMax_b,thtrueMax_a)*180.0/TMath::Pi(); // From radias to degrees

    //Azimuthal 
    double phiTrueMin = AziMin;
    double phiTrueMax = AziMax;


    std::cout << "TrueVariable | E=[ " << EtrueMin << " , " << EtrueMax << " ] Zenith= [" << thtrueMin << " , " << thtrueMax << "]" << std::endl;
    std::cout << "RecoVariable | E=[ " << ErecoMin << " , " << ErecoMax << " ] Zenith= [" << threcoMin << " , " << threcoMax << "]" << std::endl;


    //HISTOGRAMS--------------------------------------------------------------------------------------------------------
     //True Histogram        
    AsimovSimulation TrueSimulation;

    TrueSimulation.PremTable = ThePremTable; //Prem Model
    TrueSimulation.HondaTable = TheHondaTable; // Honda Table
    TrueSimulation.SetDetectorPosition(xyzTelescope); //Detector location

    //If LLVP set geometry and properties
    TrueSimulation.MantleAnomaly = PileInModel;
    TrueSimulation.AnomalyShape= ShapeOfPile;
    TrueSimulation.PileHeight = ThePileHight; // LLVP mas height 
    TrueSimulation.aperture=ThePileAperture; //Aperture
    TrueSimulation.PileDensityContrast = ThePileDensityContrast;
    TrueSimulation.PileChemContrast = ThePileChemicalContrast;


    TrueSimulation.ModifyLayer(PremLayer,DensityContrast,ChemicalContrats);
    TrueSimulation.SetIntervals(thtrueMin,thtrueMax,phiTrueMin,phiTrueMax,EtrueMin,EtrueMax);
    TrueSimulation.SetBinning(ibins,jbins,kbins);
    TrueSimulation.SetExposure(NnT);
    TrueSimulation.flvf=flvf;
    TH2D * TrueHist = TrueSimulation.GetTrueEvents2D();

    //Observed distribution:
    TH2D* RecoHist = new TH2D("RecoHist","Observed Events Histogram",mbins,CthrecoMin,CthrecoMax,nbins,ErecoMin,ErecoMax); // xbins correspond to energy values and ybins to zenith angle cost

   //Calculation of observed events------------------------------------------------------------------------------------- 
    for (int m = 1; m <= mbins; ++m) //Observed Angular bins
    {

        double cthreco = RecoHist->GetXaxis()->GetBinCenter(m); // Angular bin from the "Observed" histogram
        double threco = acos(cthreco);
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

                    double NikTrue = TrueHist->GetBinContent(i,k);

                  //  std::cout << "True bin conent: " <<  NikTrue << std::endl;

                    double dEdcthReco = dEreco*dcthreco;
                   
                    //Nsmear += PDFcth(cthreco,Etrue,cthtrue)*PDFE(Ereco,Etrue)*(NikTrue)*(dEdcthReco);

                    Nsmear += PDFE(Ereco,Etrue)*(NikTrue);

    
                } // loop e

            } // Loop ct
            
            double Nmn = Nsmear;

            RecoHist->SetBinContent(m,n,Nmn); // Set events in Observed histogram. 

         
        } // loop in Ereco

    } // loop in CthReco

   std::cout << "TrueVariable | E=[ " << EtrueMin << " , " << EtrueMax << " ] Zenith= [" << thtrueMin << " , " << thtrueMax << "] " << phiTrueMin << " " <<  phiTrueMax << std::endl;
   std::cout << "RecoVariable | E=[ " << ErecoMin << " , " << ErecoMax << " ] Zenith= [" << threcoMin << " , " << threcoMax << "]" << std::endl;
   
   delete TrueHist;
   
   return RecoHist;}
*/

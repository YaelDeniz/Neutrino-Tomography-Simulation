
/* Direction of neutrinos a re given in degrees*/
#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

//Cern ROOT
#include "TApplication.h"
#include "TRootCanvas.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "Math/Interpolator.h"
//My Libraries
#include "MathTools.h"
#include "PhyTools.h"
#include "GetTrueEvents.h"


# define myPi 3.14159265358979323846  /* pi */
// Some Constants i need
# define R_earth 6368.0 //km
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds

using namespace std;


void GetDiff3D( TH3D * histstd , TH3D * histalt, TH3D * diff);

int main(int argc, char **argv)
{
    // Plot and Histogram2D settings

    std::cout << " Neutrino Oscillation tomography. " << std::endl;

    // Number of bins
    
    int czbins=20 ; // # Bins in zenith/cos(zenith)
    int abins=20 ; // # Bins in azimuth
    int ebins=20 ; // bins in energy
    

    //Interval of integration in Zenith, Azimuth and E for Neutrino Events

    //Energy interval (in GeV)--------------------------------------------
    double EnuMin=1.0 ; 
    double EnuMax=10.0 ;

    //Zenith Angle Interval-----------------------------------------------

    //double Etamin = 33.0;
    //double Etamax = 38.0;

    // Alternative Earth Model:

    double h_llsvp = 800; //km
    double R_cmb = 3480.0;
    double R_llsvp = R_cmb + h_llsvp; //km

    double R_min = R_cmb-500.0;
    //double R_min = 3500.0; // Distance from the center of the Earth , 3500 Km CMB
    //double R_llsvp = R_cmb + h_llsvp; //Km
    double R_max = R_llsvp+500.0; // Distance from the center of the Earth

    
    double zenmax = 180-TMath::ASin( (R_min)/R_earth )*(180.0/TMath::Pi()) ; // max 180
    //double Etamax_LLSVP = TMath::ASin( (R_cmb + h_llsvp)/R_earth )*(180.0/TMath::Pi()) ;
    double zenmin = 180-TMath::ASin( R_max/R_earth )*(180.0/TMath::Pi()) ; // min 90

    double czmin = cos(zenmax*TMath::Pi()/180.0);
    
    double czmax = cos(zenmin*TMath::Pi()/180.0);

    std::cout << zenmin << " " << zenmax <<  std::endl;

    //Azimuthal Interal-------------------------------------------------------------------------
    double phimin = 0.0;
    double phimax = 360.0 ;
   

    //Detector size-----------------------------------------------------------------------------

    double DetMass = 10.0*MTon; //Mass in megaton units
    double Nn      = DetMass/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)
    double T       = 10.0*years2sec; //Detector Exposure time in sec: One Year

    double NnT = Nn*T; // Exposure Mton*years

    //Density contrast
    double drho_dp = 2.0; // 2 % density contrats

    std::vector<vector<double>> chi2data;
    double chi2;

    //Event generation-----------------------------------------------------------------------------------------
    int nuflv = 1; // neutrino  final state options nue (0), numu (1) or nutau (2)
    
    // Standard Earth model
    
    /*
    std::string prem_default   = "prem_default"; //Specify PREM table from OscProb
    TH2D* nullhist = AsimovTrueEvents(prem_default,false, layers ,flvf , Region,  Bins,  NnT);
    */

    AsimovSimulation StandardEarth;

    std::string PremName = "prem_44layers";

    StandardEarth.PremModel = PremName;
    StandardEarth.SetIntervals(zenmin,zenmax,phimin,phimax,EnuMin,EnuMax);
    StandardEarth.SetBinning(czbins,abins,ebins);
    StandardEarth.SetExposure(NnT);
    StandardEarth.flvf=nuflv;

   //TH2D * TrueStd = StandardEarth.TestTrueEvents2D("/home/dehy0499/OscProb/PremTables/prem_default.txt","/home/dehy0499/OscProb/PremTables/prem_default.txt");

     TH3D * TrueStd = StandardEarth.GetTrueEvents3D();

    // Alternative Earth Model

   AsimovSimulation AlternativeEarth;

   AlternativeEarth.PremModel = PremName;
   AlternativeEarth.SetLayerProp(23,5.0,0.0);
   AlternativeEarth.SetIntervals(zenmin,zenmax,phimin,phimax,EnuMin,EnuMax);
   AlternativeEarth.SetBinning(czbins,abins,ebins);
   AlternativeEarth.SetExposure(NnT);
   AlternativeEarth.flvf=nuflv;

   TH3D * TrueAlt = AlternativeEarth.GetTrueEvents3D();


   chi2data.push_back(layer,chi2)



   
  TH3D * TrueDiff3D = new TH3D("TrueHist","True Event Histrogram", czbins,czmin,czmax,abins,phimin,phimax,ebins,EnuMin,EnuMax); //binning in cth 

  GetDiff3D( TrueStd , TrueAlt, TrueDiff3D );

   TApplication app("app", &argc, argv);

   TCanvas *c = new TCanvas();
   gStyle->SetPalette(kBird);
   TrueDiff3D->Draw("ISO");
   TrueDiff3D->SetStats(0);
   gPad->Update();

   c->Modified(); c->Update();

  
    TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
    rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    app.Run();

    return 0;

}

void GetDiff3D( TH3D * histstd , TH3D * histalt, TH3D * diff)
{


   std::ofstream TrueDiffFile("SimulationResults/TrueEventsResults/dN_true_3D.csv"); 

   double cth, azi, e, Nexp, Nobs, dN;


    for (int j = 1; j <= diff->GetYaxis()->GetNbins(); ++j)
    {
        azi=diff->GetYaxis()->GetBinCenter(j);
        
        for (int i = 1; i <= diff->GetXaxis()->GetNbins(); ++i)
        {
            cth=diff->GetXaxis()->GetBinCenter(i);
         
            for (int k = 1; k <= diff->GetZaxis()->GetNbins(); ++k)
            {

                e=diff->GetZaxis()->GetBinCenter(k);
                Nexp = histstd->GetBinContent(i,j,k);
                Nobs = histalt->GetBinContent(i,j,k);
                dN = 100*(Nobs-Nexp)/Nexp;

                diff->SetBinContent(i,j,k,dN);

                TrueDiffFile << cth << " , " << azi << " , " << e << " , " << dN << std::endl;

            }
        }
        
    }

    TrueDiffFile.close();

}

double GetChi2( TH3D * histstd, TH3D * histalt)
{

    double = Nexp, Nobs;

    double chi2 = 0;

    for (int j = 1; j <= diff->GetYaxis()->GetNbins(); ++j) //loop in azimuth
    {
        
        for (int i = 1; i <= diff->GetXaxis()->GetNbins(); ++i) //loop in zenith
        {
         
            for (int k = 1; k <= diff->GetZaxis()->GetNbins(); ++k) //loop in energy
            {

                Nexp = histstd->GetBinContent(i,j,k);
                Nobs = histalt->GetBinContent(i,j,k);
                chi2 +=  2*( Nexp - Nobs + Nobs*TMath::Log(Nobs/Nexp) ); // LLRT

            }
        }
        
    }

    return chi2;

}